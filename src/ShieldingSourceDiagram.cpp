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
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include <nlohmann/json.hpp>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceDiagram.h"

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
  std::stringstream options;
  options << "{ padding: { top: 40, right: 40, bottom: 40, left: 40 } }";
  
  std::stringstream ss;
  ss << "new Shielding2DView('" << id() << "', " << jsonData << ", " << options.str() << ");";
  
  setJavaScriptMember( "chart", ss.str() );
  
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
    switch(geometry) {
        case GammaInteractionCalc::GeometryType::Spherical: geomStr = "Spherical"; break;
        case GammaInteractionCalc::GeometryType::CylinderEndOn: geomStr = "CylinderEndOn"; break;
        case GammaInteractionCalc::GeometryType::CylinderSideOn: geomStr = "CylinderSideOn"; break;
        case GammaInteractionCalc::GeometryType::Rectangular: geomStr = "Rectangular"; break;
        default: geomStr = "Unknown"; break;
    }
    j["geometry"] = geomStr;
    
    // Distance in mm
    j["distance"] = detectorDistance / PhysicalUnits::mm;
    
    // Detector Diameter in mm
    j["detectorDiameter"] = detectorDiameter / PhysicalUnits::mm;
    
    j["fitSources"] = nlohmann::json::array();
    for( const auto &src : sources ) {
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

    for( const auto &s : shieldings ) {
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
        
        for( const auto &ts : s.m_traceSources ) {
            std::stringstream ss_src;
            ss_src << "Trace: ";
            if( ts.m_nuclide ) ss_src << ts.m_nuclide->symbol;
            
            ss_src << " (" << PhysicalUnits::printToBestActivityUnits( ts.m_activity );
            
            switch( ts.m_type ) {
                case GammaInteractionCalc::TraceActivityType::TotalActivity:
                    break;
                case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
                    ss_src << "/cm3";
                    break;
                case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
                    ss_src << "/g";
                    break;
                case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
                    ss_src << "/m2";
                    break;
                default:
                    break;
            }
            ss_src << ")";
            sourcesList.push_back( ss_src.str() );
        }
        
        for( const auto &elemPair : s.m_nuclideFractions_ ) {
             for( const auto &tup : elemPair.second ) {
                 const auto *nuc = std::get<0>(tup);
                 if( nuc ) {
                     sourcesList.push_back( std::string("Self: ") + nuc->symbol );
                 }
             }
        }
        s_json["sources"] = sourcesList;
        s_json["hasSource"] = !sourcesList.empty();
        
        double volume = 0.0;
        double mass = 0.0;
        
        if( !s.m_isGenericMaterial ) {
            const double pi = PhysicalUnits::pi;
            
            switch( geometry ) {
                case GammaInteractionCalc::GeometryType::Spherical: {
                    const double thickness = s.m_dimensions[0];
                    const double outer_rad = inner_rad + thickness;
                    const double innerVol = (4.0/3.0) * pi * inner_rad * inner_rad * inner_rad;
                    const double outerVol = (4.0/3.0) * pi * outer_rad * outer_rad * outer_rad;
                    volume = outerVol - innerVol;
                    inner_rad = outer_rad;
                    break;
                }
                case GammaInteractionCalc::GeometryType::CylinderEndOn:
                case GammaInteractionCalc::GeometryType::CylinderSideOn: {
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
                case GammaInteractionCalc::GeometryType::Rectangular: {
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
                default:
                    break;
            }
            
            if( s.m_material ) {
                mass = volume * s.m_material->density;
            }
        }
        
        s_json["mass"] = mass / PhysicalUnits::gram;
        s_json["volume"] = volume / PhysicalUnits::cm3;
        
        j["shieldings"].push_back(s_json);
    }
    
    return j.dump();
}

std::string Shielding2DView::createJsonData() const
{
  return createShieldingDiagramJson( m_shieldings, m_sources, m_geometry, m_detectorDistance, m_detectorDiameter );
}

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
  std::string jsonData = createJsonData();
  
  std::stringstream ss;
  ss << "new Shielding3DView('" << id() << "', " << jsonData << ");";
  
  doJavaScript( ss.str() );
}

std::string Shielding3DView::createJsonData() const
{
  return createShieldingDiagramJson( m_shieldings, m_sources, m_geometry, m_detectorDistance, m_detectorDiameter );
}

// Dialog creation function (moved from ShieldingSourceDisplay.cpp)
SimpleDialog *createShieldingDiagram(
  const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
  GammaInteractionCalc::GeometryType geometry,
  double detectorDistance,
  double detectorDiameter
)
{
  using namespace Wt;
  
  InterSpec *viewer = InterSpec::instance();
  if( viewer )
    viewer->useMessageResourceBundle( "ShieldingSourceDisplay" ); //
  
  SimpleDialog *dialog = new SimpleDialog( WString::tr("ssd-shield-source-diagram") );
  dialog->addStyleClass( "ShieldingDiagramDialog" );
  dialog->resize( WLength(95,WLength::Percentage), WLength(95,WLength::Percentage) );
  
  Shielding2DView *view = new Shielding2DView( shieldings, sources, geometry, detectorDistance, detectorDiameter );
  
  WGridLayout *layout = new WGridLayout();
  layout->addWidget( view, 0, 0 );
  layout->setColumnStretch( 0, 1 );
  layout->setRowStretch( 0, 1 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  dialog->contents()->setLayout( layout );
  dialog->contents()->setOverflow( Wt::WContainerWidget::Overflow::OverflowHidden );

  
  WContainerWidget *type_row = new WContainerWidget();
  layout->addWidget( type_row, 1, 0 );
  WComboBox *select = new WComboBox( type_row );
  select->setFloatSide( Wt::Side::Right );
  select->addItem( "2D View" );
  select->addItem( "3D View" );
  select->setCurrentIndex( 0 );
  select->activated().connect( std::bind([select, layout, shieldings, sources, geometry, detectorDistance, detectorDiameter](){
    
    bool found_old = false;
    for( int i = 0; !found_old && (i < layout->count()); ++i )
    {
      WLayoutItem *item = layout->itemAt(i);
      WWidget *w = item ? item->widget() : nullptr;
      if( Shielding2DView *view2D = dynamic_cast<Shielding2DView *>(w) )
      {
        found_old = true;
        delete view2D;
      }else if( Shielding3DView *view3D = dynamic_cast<Shielding3DView *>(w) )
      {
        found_old = true;
        delete view3D;
      }
    }
    
    assert( found_old );
    if( !found_old )
      return;
    
    if( select->currentIndex() == 0 )
    {
      Shielding2DView *view = new Shielding2DView( shieldings, sources, geometry, detectorDistance, detectorDiameter );
      layout->addWidget( view, 0, 0 );
    }else
    {
      Shielding3DView *view = new Shielding3DView( shieldings, sources, geometry, detectorDistance, detectorDiameter );
      layout->addWidget( view, 0, 0 );
    }
  }) );
  
  static const std::string chart_css_rule_size_name = "ShieldingDiagramDialog-size";
  Wt::WCssStyleSheet &style = wApp->styleSheet();
  if( !style.isDefined(chart_css_rule_size_name) )
  {
    style.addRule( ".simple-dialog.ShieldingDiagramDialog", "max-width: 95vw; max-height: 95vh; width: 95vw; height: 95vh;", chart_css_rule_size_name );
    style.addRule( ".simple-dialog.ShieldingDiagramDialog .body", "max-width: 95vw;", chart_css_rule_size_name + "-body");
  }
  
  dialog->addButton( WString::tr("Close") );
  
  return dialog;
}

