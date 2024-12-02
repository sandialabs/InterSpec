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


#include <map>
#include <memory>
#include <vector>


#include <Wt/Utils>
#include <Wt/WColor>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>


#include <Wt/Json/Array>
#include <Wt/Json/Value>
#include <Wt/Json/Object>
#include <Wt/Json/Serializer>


#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/Chi2Graphic.h"
#include "InterSpec/GammaInteractionCalc.h"


#if( SHOW_Chi2Graphic_PEAK_PREVIEW )
#include <Wt/WSvgImage>

#include "SpecUtils/SpecFile.h"

#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#endif //SHOW_Chi2Graphic_PEAK_PREVIEW

using namespace Wt;
using namespace std;


Chi2Graphic::Chi2Graphic( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_jsgraph( jsRef() + ".chart" ),
  m_xAxisTitle( WString::tr("Energy (keV)") ),
  m_yAxisChiTitle( WString::tr("x2g-yaxis-title-chi") ),    //From ShieldingSourceDisplay.xml
  m_yAxisScaleTitle( WString::tr("x2g-yaxis-title-mult") ), //From ShieldingSourceDisplay.xml
  m_displayType( Chi2Graphic::ChiDispType::Pull ),
  m_topMargin( 5 ),
  m_rightMargin( 5 ),
  m_bottomMargin( 0 ),
  m_leftMargin( 5 ),
  m_titlePadding( 1 ),
  m_pendingJs{},
  m_renderFlags( 0 ),
  m_dataPointMousedOver( this, "dataPointMousedOver", false ),
  m_ndof( 0 ),
  m_used_points{}
{
  addStyleClass( "Chi2Graphic" );
  setOverflow(Overflow::OverflowHidden);
  
  wApp->useStyleSheet( "InterSpec_resources/Chi2Graphic.css" );
  
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/Chi2Graphic.js" );
  
  m_renderFlags |= Chi2RenderActions::SetData;
  
#if( SHOW_Chi2Graphic_PEAK_PREVIEW )
  m_dataPointMousedOver.connect( boost::bind( &Chi2Graphic::dataPointMousedOver, this,
                                             boost::placeholders::_1, boost::placeholders::_2 ) );
#endif
}//Chi2Graphic constructor


void Chi2Graphic::defineJavaScript()
{
  string options = "{ "
  "margins: {"
  " top: " + std::to_string(m_topMargin) + ","
  " right: " + std::to_string(m_rightMargin) + ","
  " bottom: " + std::to_string(m_bottomMargin) + ","
  " left: " + std::to_string(m_leftMargin) +
  " }, "
  " titleToAxisPadding: " + std::to_string(m_titlePadding);
  
  
  string jsDispType;
  switch( m_displayType )
  {
    case ChiDispType::Pull:  jsDispType = "0"; break;
    case ChiDispType::Scale: jsDispType += "1"; break;
  }

  if( !m_xAxisTitle.empty() )
    options += ", xAxisTitle: " + m_xAxisTitle.jsStringLiteral('\"');
  if( !m_yAxisChiTitle.empty() )
    options += ", yAxisChiTitle: " + m_yAxisChiTitle.jsStringLiteral('\"');
  if( !m_yAxisScaleTitle.empty() )
    options += ", yAxisScaleTitle: " + m_yAxisScaleTitle.jsStringLiteral('\"');
  options += ", displayType: " + jsDispType;
  options += " }";
  
  setJavaScriptMember( "chart", "new Chi2Graphic(" + jsRef() + ", " + options + ");");
  
  setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + id() + "') ){"
          // When we "Clear Session", jsRef() will give a null result temporarily, so we'll protect against that
          "const c=" + jsRef() + ";"
          "if(c && c.chart)"
            "c.chart.handleResize();"
        "}"
      "}"
    "});"
  );
  
  callJavaScriptMember( "resizeObserver.observe", jsRef() );
  
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void Chi2Graphic::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
  
  if( m_renderFlags.testFlag(Chi2RenderActions::SetData) )
    setDataToClient();
  
  m_renderFlags = 0;
}//void render( flags )


void Chi2Graphic::setData( const int ndof, const std::vector<GammaInteractionCalc::PeakResultPlotInfo> &points )
{
  m_ndof = ndof;
  m_used_points = points;
  m_renderFlags |= Chi2RenderActions::SetData;
  scheduleRender();
}//void setData(ndof, info)


string Chi2Graphic::to_json( const int ndof,
                      const std::vector<GammaInteractionCalc::PeakResultPlotInfo> &used_points )
{
  Wt::Json::Object json;
  
  json["NDOF"] = ndof;
  Wt::Json::Array &points = json["Points"] = Wt::Json::Array{};
  
  for( const GammaInteractionCalc::PeakResultPlotInfo &point : used_points )
  {
    Wt::Json::Object jpoint;
    jpoint["energy"] = point.energy;
    jpoint["numSigmaOff"] = point.numSigmaOff;
    jpoint["observedOverExpected"] = point.observedOverExpected;
    jpoint["observedOverExpectedUncert"] = point.observedOverExpectedUncert;
    const auto peak = point.m_user_info ? point.m_user_info->m_peak : nullptr;
    WColor color;
    if( peak )
      color = peak->lineColor();
    if( color.isDefault() )
    {
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        color = viewer->getColorTheme()->defaultPeakLine;
      else
        color = WColor("blue");
    }//if( color.isDefault() )
    color.setRgb( color.red(), color.green(), color.blue(), 255 );
    
    jpoint["color"] = WString::fromUTF8(color.cssText());
    string nuclide = "n/a";
    if( peak && peak->parentNuclide() )
      nuclide = peak->parentNuclide()->symbol;
    jpoint["nuclide"] = WString::fromUTF8(nuclide);
    
    const shared_ptr<const GammaInteractionCalc::PeakResultPlotInfo::ExtraInfo> &user_info = point.m_user_info;
    if( user_info )
    {
      jpoint["numObserved"] = user_info->m_observed;
      jpoint["numObservedUncert"] = user_info->m_observed_uncert;
      jpoint["numExpected"] = user_info->m_expected;
      if( user_info->m_back_peak_area > 0 )
      {
        jpoint["numBackground"] = user_info->m_back_peak_area;
        jpoint["numBackgroundUncert"] = user_info->m_back_peak_area_uncert;
      }
      
      if( user_info->m_det_eff > 0 )
        jpoint["detEff"] = user_info->m_det_eff;
      
      if( user_info->m_det_intrinsic_eff > 0 )
        jpoint["detIntrinsicEff"] = user_info->m_det_intrinsic_eff;
      
      if( user_info->m_det_solid_angle > 0 )
        jpoint["detSolidAngle"] = user_info->m_det_solid_angle;
      
      if( !user_info->m_calc_details.m_sources.empty() )
      {
        Wt::Json::Array &srcs = jpoint["sourcesInfo"] = Wt::Json::Array{};
        
        for( const GammaInteractionCalc::PeakDetailSrc &src : user_info->m_calc_details.m_sources )
        {
          Wt::Json::Object jsonsrc;
          jsonsrc["name"] = WString::fromUTF8( src.nuclide ? src.nuclide->symbol : std::string{} );
          jsonsrc["energy"] = src.energy;
          jsonsrc["br"] = src.br;
          jsonsrc["cps"] = src.modelContribToPeak / user_info->m_live_time;
          jsonsrc["counts"] = src.modelContribToPeak;
          
          //src.age;
          //src.calcActivity
          //... there are other quantities we could put here, but its not clear we would use
          if( src.decayCorrection > 0.0 )
            jsonsrc["decayCorrection"] = src.decayCorrection;
          
          srcs.push_back( jsonsrc );
        }//for( loop over point.m_calc_details.m_sources )
      }//if( !user_info->m_sources_info.empty() )
      
      if( !user_info->m_calc_details.m_volumetric_srcs.empty() )
      {
        Wt::Json::Array &srcs = jpoint["volumetricSources"] = Wt::Json::Array{};
        
        for( const GammaInteractionCalc::PeakDetail::VolumeSrc &src : user_info->m_calc_details.m_volumetric_srcs )
        {
          Wt::Json::Object jsonsrc;
          jsonsrc["name"] = WString::fromUTF8( src.sourceName );
          jsonsrc["type"] = WString::fromUTF8( src.trace ? "Trace" : "Self-Attenuating" );
          jsonsrc["volumetricActivity"] = src.srcVolumetricActivity;
          jsonsrc["detectorIntrinsicEff"] = src.detIntrinsicEff;
          if( src.trace )
          {
            jsonsrc["inSituExponential"] = src.inSituExponential;
            jsonsrc["relaxationLength_mm"] = src.inSituRelaxationLength; //TODO: should print to proper length, or divide by `PhysicalUnits::mm`
          }
          srcs.push_back( jsonsrc );
        }//for( loop over user_info->m_volumetric_srcs )
      }//if( !user_info->m_volumetric_srcs.empty() )
    }//if( point.m_user_info )
    
    points.push_back( std::move(jpoint) );
  }//for( loop over used_points )
  
  return Wt::Json::serialize(json);
}//to_json(...)


void Chi2Graphic::setDataToClient()
{
  const string json = to_json(m_ndof, m_used_points);
  const string js = m_jsgraph + ".setData(" + json + ");";
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setData( std::shared_ptr<Measurement> data_hist )


#if( SHOW_Chi2Graphic_PEAK_PREVIEW )
void Chi2Graphic::dataPointMousedOver( const std::string &tooltip_id, const double energy )
{
  InterSpec *viewer = InterSpec::instance();
  PeakModel *peakModel = viewer ? viewer->peakModel() : nullptr;
  assert( peakModel );
  
  shared_ptr<const SpecUtils::Measurement> measurement = viewer ? viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground) : nullptr;
  shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = peakModel ? peakModel->peaks() : nullptr;
  
  if( !measurement || !peaks )
    return;
  
  double lower_energy = -1.0, upper_energy = -1.0, delta = 99999.9;
  
  for( const shared_ptr<const PeakDef> &p : *peaks )
  {
    try
    {
      const double d = fabs(p->gammaParticleEnergy() - energy);
      if( (d < 1.0) && (d < delta) )
      {
        delta = d;
        const double roi_width = p->upperX() - p->lowerX();
        lower_energy = p->lowerX() - roi_width;
        upper_energy = p->upperX() + roi_width;
      }
    }catch( std::exception & )
    {
      // We get here if a peak doesnt have a nuclide/x-ray/reaction associated with it.
    }
  }//for( const shared_ptr<const PeakDef> &p : *peaks )
  
  if( delta > 2 )
    return;
  
  const bool compact = true;
  shared_ptr<const ColorTheme> theme = viewer->getColorTheme();
  const vector<std::shared_ptr<const ReferenceLineInfo>> reflines;
  const int width_px = 300, height_px = 160;
  
  shared_ptr<Wt::WSvgImage> svg = PeakSearchGuiUtils::renderChartToSvg( measurement, peaks,
                                                 reflines, lower_energy, upper_energy,
                                                width_px, height_px, theme, compact );
  if( svg )
  {
    stringstream strm;
    svg->write( strm );
    
    const string svgstr = Wt::WWebWidget::jsStringLiteral(strm.str(),'\'');
    string js = "$('#" + tooltip_id + " > .tt-chart').html(" + svgstr + ");";
    doJavaScript(js);
  }
}//void dataPointMousedOver( const std::string &tooltip_id, const double energy )
#endif //SHOW_Chi2Graphic_PEAK_PREVIEW

void Chi2Graphic::setXAxisTitle( const Wt::WString &title )
{
  m_xAxisTitle = title;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setXAxisTitle(" + m_xAxisTitle.jsStringLiteral('\'') + ",true);" );
}//


void Chi2Graphic::setYAxisTitle( const Wt::WString &chi_title, const Wt::WString &scale_title )
{
  m_yAxisChiTitle = chi_title;
  m_yAxisScaleTitle = scale_title;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setYAxisTitles(" + m_yAxisChiTitle.jsStringLiteral('\'') + ","
                 + m_yAxisScaleTitle.jsStringLiteral('\'') + ",true);" );
}


void Chi2Graphic::setContentMargins( int top, int right, int bottom, int left )
{
  m_topMargin = top;
  m_rightMargin = right;
  m_bottomMargin = bottom;
  m_leftMargin = left;
  
  if( isRendered() )
  {
    const string margins = "{"
    " top: " + std::to_string(m_topMargin) + ","
    " right: " + std::to_string(m_rightMargin) + ","
    " bottom: " + std::to_string(m_bottomMargin) + ","
    " left: " + std::to_string(m_leftMargin) +
    " }";
    doJavaScript( m_jsgraph + ".setMargins(" + margins + ");" );
  }//if( isRendered() )
}//void setContentMargins( int top, int right, int bottom, int left );


Chi2Graphic::~Chi2Graphic()
{
}//~Chi2Graphic()


void Chi2Graphic::setChartType( const Chi2Graphic::ChiDispType type )
{
  m_displayType = type;
  
  string jsDispType = m_jsgraph + ".ChartType.";
  switch( m_displayType )
  {
    case ChiDispType::Pull:  jsDispType += "Pull"; break;
    case ChiDispType::Scale: jsDispType += "Scale"; break;
  }
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setDisplayType(" + jsDispType + ");" );
}//void setChartType( type )
