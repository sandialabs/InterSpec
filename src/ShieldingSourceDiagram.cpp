/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include "InterSpec_config.h"

#include <sstream>
#include <iomanip>

#include <Wt/WString>
#include <Wt/WLogger>
#include <Wt/WComboBox>
#include <Wt/WGridLayout>
#include <Wt/WLayoutItem>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>
#include <Wt/WWidget>

#include <nlohmann/json.hpp>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceDiagram.h"

using namespace Wt;
using namespace std;

Shielding2DView::Shielding2DView( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                                  GammaInteractionCalc::GeometryType geometry,
                                  double detectorDistance,
                                  double detectorDiameter,
                                  Wt::WContainerWidget *parent )
  : Wt::WContainerWidget( parent ),
    m_shieldings( shieldings ),
    m_sources( sources ),
    m_geometry( geometry ),
    m_detectorDistance( detectorDistance ),
    m_detectorDiameter( detectorDiameter )
{
  setStyleClass("Shielding2DView");
  
  // Load JS and CSS
  Wt::WApplication *app = Wt::WApplication::instance();
  app->require("InterSpec_resources/d3.v3.min.js", "d3.v3.js");
  app->require("InterSpec_resources/Shielding2DView.js");
  app->useStyleSheet("InterSpec_resources/Shielding2DView.css");
  
  defineJavaScript();
}

void Shielding2DView::defineJavaScript()
{
  std::string jsonData = createJsonData();
  
  // Create options object with default padding
  string options = "{ padding: { top: 40, right: 40, bottom: 40, left: 40 } }";
  string js = "var c=" + jsRef() + ";if(c){c.chart=new Shielding2DView('" + id() + "', " + jsonData + ", " + options + ");}";
  
  doJavaScript( js );
  
  // Set up ResizeObserver
  setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + id() + "') ){"
          "const c=" + jsRef() + ";"
          "if(c && c.chart)"
            "c.chart.handleResize();"
        "}"
      "}"
    "});"
  );
  
  callJavaScriptMember( "resizeObserver.observe", jsRef() );
}

// JSON creation code (moved from ShieldingJsonUtils)
std::string createShieldingDiagramJson(
    const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
    const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
    GammaInteractionCalc::GeometryType geometry,
    double detectorDistance,
    double detectorDiameter
)
{
  nlohmann::json j;
    
  // Geometry
  std::string geomStr;
  switch( geometry )
  {
    case GammaInteractionCalc::GeometryType::Spherical:
      geomStr = "Spherical";
      break;
      
    case GammaInteractionCalc::GeometryType::CylinderEndOn:
      geomStr = "CylinderEndOn";
      break;
      
    case GammaInteractionCalc::GeometryType::CylinderSideOn:
      geomStr = "CylinderSideOn";
      break;
      
    case GammaInteractionCalc::GeometryType::Rectangular:
      geomStr = "Rectangular";
      break;
      
    case GammaInteractionCalc::GeometryType::NumGeometryType:
      assert( 0 );
      geomStr = "Unknown";
      break;
  }//switch( geometry )
  
  j["geometry"] = geomStr;
  
  // Distance in mm
  j["distance"] = detectorDistance / PhysicalUnits::mm;
  
  // Detector Diameter in mm
  j["detectorDiameter"] = detectorDiameter / PhysicalUnits::mm;
  
  j["fitSources"] = nlohmann::json::array();
  for( const ShieldingSourceFitCalc::IsoFitStruct &src : sources )
  {
    nlohmann::json src_json;
    std::string symbol = src.nuclide ? src.nuclide->symbol : "Unknown";
    src_json["nuclide"] = symbol;
    
    // Type
    std::string typeStr;
    switch( src.sourceType ) {
      case ShieldingSourceFitCalc::ModelSourceType::Point: typeStr = "Point"; break;
      case ShieldingSourceFitCalc::ModelSourceType::Intrinsic: typeStr = "Intrinsic"; break;
      case ShieldingSourceFitCalc::ModelSourceType::Trace: typeStr = "Trace"; break;
    }
    src_json["type"] = typeStr;
    
    // Activity
    src_json["activity"] = src.activity;
    src_json["activityPretty"] = PhysicalUnits::printToBestActivityUnits( src.activity );
    
    if( src.activityUncertainty > 0 ) {
      src_json["activityUncert"] = src.activityUncertainty;
    }
    
    // Age
    src_json["age"] = src.age;
    
    j["fitSources"].push_back( src_json );
  }
  
  j["shieldings"] = nlohmann::json::array();
  
  // Logic to track inner dimensions for volume calculation
  double inner_rad = 0.0;
  double inner_cyl_rad = 0.0;
  double inner_cyl_half_len = 0.0;
  double inner_half_width = 0.0;
  double inner_half_height = 0.0;
  double inner_half_depth = 0.0;
  
  for( const ShieldingSourceFitCalc::ShieldingInfo &s : shieldings )
  {
    nlohmann::json s_json;
    
    std::string matName = "Generic";
    double density = 0.0;
    
    if( s.m_material ) {
      matName = s.m_material->name;
      density = s.m_material->density / (PhysicalUnits::g / PhysicalUnits::cm3);
    } else if( s.m_isGenericMaterial ) {
      matName = "Generic(Z=" + std::to_string(int(s.m_dimensions[0])) + ")";
      density = 0.0;
      s_json["atomicNumber"] = s.m_dimensions[0];
      s_json["arealDensity"] = s.m_dimensions[1] / (PhysicalUnits::g / (PhysicalUnits::cm2));
    }
    
    s_json["material"] = matName;
    s_json["density"] = density;
    
    s_json["dimensions"] = {
      s.m_dimensions[0] / PhysicalUnits::mm,
      s.m_dimensions[1] / PhysicalUnits::mm,
      s.m_dimensions[2] / PhysicalUnits::mm
    };
    
    std::vector<std::string> sourcesList;
    
    for( const auto &ts : s.m_traceSources )
    {
      string ss_src = "Trace: ";
      if( ts.m_nuclide )
        ss_src += ts.m_nuclide->symbol;
      
      ss_src += " (" + PhysicalUnits::printToBestActivityUnits( ts.m_activity );
      
      switch( ts.m_type )
      {
        case GammaInteractionCalc::TraceActivityType::TotalActivity:
        case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
          break;
          
        case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
          ss_src += "/cm3";
          break;
          
        case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
          ss_src += "/g";
          break;
          
        case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
          ss_src += "/m2";
          break;
      }//switch( ts.m_type )
      
      ss_src += ")";
      
      sourcesList.push_back( std::move(ss_src) );
    }//for( const auto &ts : s.m_traceSources )
    
    for( const auto &elemPair : s.m_nuclideFractions_ )
    {
      for( const auto &tup : elemPair.second )
      {
        if( std::get<0>(tup) )
          sourcesList.push_back( std::string("Self: ") + std::get<0>(tup)->symbol );
      }
    }
    
    s_json["sources"] = sourcesList;
    s_json["hasSource"] = !sourcesList.empty();
    
    double volume = 0.0;
    double mass = 0.0;
    
    if( !s.m_isGenericMaterial )
    {
      const double pi = PhysicalUnits::pi;
      
      switch( geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
        {
          const double thickness = s.m_dimensions[0];
          const double outer_rad = inner_rad + thickness;
          const double innerVol = (4.0/3.0) * pi * inner_rad * inner_rad * inner_rad;
          const double outerVol = (4.0/3.0) * pi * outer_rad * outer_rad * outer_rad;
          volume = outerVol - innerVol;
          inner_rad = outer_rad;
          break;
        }
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
        {
          const double radThickness = s.m_dimensions[0];
          const double lenThickness = s.m_dimensions[1];
          const double rad = inner_cyl_rad + radThickness;
          const double half_len = inner_cyl_half_len + lenThickness;
          const double innerVol = pi * inner_cyl_rad * inner_cyl_rad * (2.0 * inner_cyl_half_len);
          const double outerVol = pi * rad * rad * (2.0 * half_len);
          volume = outerVol - innerVol;
          inner_cyl_rad = rad;
          inner_cyl_half_len = half_len;
          break;
        }
          
        case GammaInteractionCalc::GeometryType::Rectangular:
        {
          const double wThick = s.m_dimensions[0];
          const double hThick = s.m_dimensions[1];
          const double dThick = s.m_dimensions[2];
          const double half_width = inner_half_width + wThick;
          const double half_height = inner_half_height + hThick;
          const double half_depth = inner_half_depth + dThick;
          const double innerVol = 8.0 * inner_half_width * inner_half_height * inner_half_depth;
          const double outerVol = 8.0 * half_width * half_height * half_depth;
          volume = outerVol - innerVol;
          inner_half_width = half_width;
          inner_half_height = half_height;
          inner_half_depth = half_depth;
          break;
        }
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( geometry )
      
      if( s.m_material )
        mass = volume * s.m_material->density;
    }//if( !s.m_isGenericMaterial )
    
    s_json["mass"] = mass / PhysicalUnits::gram;
    s_json["volume"] = volume / PhysicalUnits::cm3;
    
    j["shieldings"].push_back(s_json);
  }//for( const ShieldingSourceFitCalc::ShieldingInfo &s : shieldings )
  
  return j.dump();
}//std::string createShieldingDiagramJson(...)


std::string Shielding2DView::createJsonData() const
{
  return createShieldingDiagramJson( m_shieldings, m_sources, m_geometry, m_detectorDistance, m_detectorDiameter );
}

void Shielding2DView::updateData( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                                  GammaInteractionCalc::GeometryType geometry,
                                  double detectorDistance,
                                  double detectorDiameter )
{
  m_shieldings = shieldings;
  m_sources = sources;
  m_geometry = geometry;
  m_detectorDistance = detectorDistance;
  m_detectorDiameter = detectorDiameter;
  
  // Update JavaScript chart with new data
  if( isRendered() )
  {
    string jsonData = createJsonData();
    string js = "var c=" + jsRef() + ";if(c && c.chart && typeof c.chart.setData === 'function'){c.chart.setData(" + jsonData + ");}";
    doJavaScript( js );
  }
}//void Shielding2DView::updateData(...)

// Shielding3DView implementation (moved from Shielding3DView.cpp)
Shielding3DView::Shielding3DView( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                                  GammaInteractionCalc::GeometryType geometry,
                                  double detectorDistance,
                                  double detectorDiameter,
                                  Wt::WContainerWidget *parent )
  : Wt::WContainerWidget( parent ),
    m_shieldings( shieldings ),
    m_sources( sources ),
    m_geometry( geometry ),
    m_detectorDistance( detectorDistance ),
    m_detectorDiameter( detectorDiameter )
{
  setStyleClass("Shielding3DView");
  
  Wt::WApplication *app = Wt::WApplication::instance();
  app->require("InterSpec_resources/Shielding3DView.js?v=2");
  app->useStyleSheet("InterSpec_resources/Shielding3DView.css");
  
  defineJavaScript();
}

void Shielding3DView::defineJavaScript()
{
  string jsonData = createJsonData();
  string js = "var c=" + jsRef() + ";if(c){c.chart=new Shielding3DView('" + id() + "', " + jsonData + ");}";
  doJavaScript( js );
}

std::string Shielding3DView::createJsonData() const
{
  return createShieldingDiagramJson( m_shieldings, m_sources, m_geometry, m_detectorDistance, m_detectorDiameter );
}

void Shielding3DView::updateData( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                                  GammaInteractionCalc::GeometryType geometry,
                                  double detectorDistance,
                                  double detectorDiameter )
{
  m_shieldings = shieldings;
  m_sources = sources;
  m_geometry = geometry;
  m_detectorDistance = detectorDistance;
  m_detectorDiameter = detectorDiameter;
  
  // Update JavaScript 3D view with new data
  if( isRendered() )
  {
    string jsonData = createJsonData();
    string js = "var c=" + jsRef() + ";if(c && c.chart && typeof c.chart.setData === 'function'){c.chart.setData(" + jsonData + ");}";
    doJavaScript( js );
  }
}//void Shielding3DView::updateData(...)

// ShieldingDiagramDialog implementation
ShieldingDiagramDialog::ShieldingDiagramDialog(
  const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
  GammaInteractionCalc::GeometryType geometry,
  double detectorDistance,
  double detectorDiameter
)
  : SimpleDialog( Wt::WString::tr("ssd-shield-source-diagram") ),
    m_2DView( nullptr ),
    m_3DView( nullptr ),
    m_select( nullptr ),
    m_layout( nullptr ),
    m_shieldings( shieldings ),
    m_sources( sources ),
    m_geometry( geometry ),
    m_detectorDistance( detectorDistance ),
    m_detectorDiameter( detectorDiameter )
{
  InterSpec *viewer = InterSpec::instance();
  if( viewer )
    viewer->useMessageResourceBundle( "ShieldingSourceDisplay" );
  
  addStyleClass( "ShieldingDiagramDialog" );
  resize( WLength(95,WLength::Percentage), WLength(95,WLength::Percentage) );
  
  // Create initial 2D view
  m_2DView = new Shielding2DView( m_shieldings, m_sources, m_geometry, m_detectorDistance, m_detectorDiameter );
  
  m_layout = new WGridLayout();
  m_layout->addWidget( m_2DView, 0, 0 );
  m_layout->setColumnStretch( 0, 1 );
  m_layout->setRowStretch( 0, 1 );
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  contents()->setLayout( m_layout );
  contents()->setOverflow( Wt::WContainerWidget::Overflow::OverflowHidden );
  
  WContainerWidget *type_row = new WContainerWidget();
  m_layout->addWidget( type_row, 1, 0 );
  m_select = new WComboBox( type_row );
  m_select->setFloatSide( Wt::Side::Right );
  m_select->addItem( "2D View" );
  m_select->addItem( "3D View" );
  m_select->setCurrentIndex( 0 );
  m_select->activated().connect( this, &ShieldingDiagramDialog::handleViewTypeToggle );
  
  static const std::string chart_css_rule_size_name = "ShieldingDiagramDialog-size";
  Wt::WCssStyleSheet &style = wApp->styleSheet();
  if( !style.isDefined(chart_css_rule_size_name) )
  {
    style.addRule( ".simple-dialog.ShieldingDiagramDialog", "max-width: 95vw; max-height: 95vh; width: 95vw; height: 95vh;", chart_css_rule_size_name );
    style.addRule( ".simple-dialog.ShieldingDiagramDialog .body", "max-width: 95vw;", chart_css_rule_size_name + "-body");
  }
  
  addButton( WString::tr("Close") );
}


void ShieldingDiagramDialog::handleViewTypeToggle()
{
  const bool show3D = (m_select->currentIndex() == 1);
  switchView( show3D );
}//void handleViewTypeToggle()

void ShieldingDiagramDialog::updateData( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                         const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                                         GammaInteractionCalc::GeometryType geometry,
                                         double detectorDistance,
                                         double detectorDiameter )
{
  // Update stored data
  m_shieldings = shieldings;
  m_sources = sources;
  m_geometry = geometry;
  m_detectorDistance = detectorDistance;
  m_detectorDiameter = detectorDiameter;
  
  // Update the currently visible view
  const bool show3D = (m_select && m_select->currentIndex() == 1);
  
  if( show3D && m_3DView )
  {
    m_3DView->updateData( shieldings, sources, geometry, detectorDistance, detectorDiameter );
  }
  else if( !show3D && m_2DView )
  {
    m_2DView->updateData( shieldings, sources, geometry, detectorDistance, detectorDiameter );
  }
  
  // Also update the other view if it exists (so switching views will show updated data)
  if( m_2DView )
    m_2DView->updateData( shieldings, sources, geometry, detectorDistance, detectorDiameter );
  
  if( m_3DView )
    m_3DView->updateData( shieldings, sources, geometry, detectorDistance, detectorDiameter );
}//void ShieldingDiagramDialog::updateData(...)


void ShieldingDiagramDialog::switchView( bool show3D )
{
  if( show3D && !m_3DView )
  {
    if( m_2DView )
      delete m_2DView;
    m_2DView = nullptr;
    
    m_3DView = new Shielding3DView( m_shieldings, m_sources, m_geometry, m_detectorDistance, m_detectorDiameter );
    m_layout->addWidget( m_3DView, 0, 0 );
  }//if( show3D && !m_3DView )
  
  if( !show3D && !m_2DView )
  {
    if( m_3DView )
      delete m_3DView;
    m_3DView = nullptr;
    
    m_2DView = new Shielding2DView( m_shieldings, m_sources, m_geometry, m_detectorDistance, m_detectorDiameter );
    m_layout->addWidget( m_2DView, 0, 0 );
  }//if( !show3D && !m_2DView )
}

// Static factory method
ShieldingDiagramDialog *ShieldingDiagramDialog::createShieldingDiagram(
  const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
  GammaInteractionCalc::GeometryType geometry,
  double detectorDistance,
  double detectorDiameter
)
{
  return new ShieldingDiagramDialog( shieldings, sources, geometry, detectorDistance, detectorDiameter );
}

