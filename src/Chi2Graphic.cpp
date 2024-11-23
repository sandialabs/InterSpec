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
#include "InterSpec/Chi2Graphic.h"


using namespace Wt;
using namespace std;


Chi2Graphic::Chi2Graphic( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_jsgraph( jsRef() + ".chart" ),
  m_xAxisTitle{},
  m_yAxisTitle{},
  m_displayType( Chi2Graphic::ChiDispType::Pull ),
  m_topMargin( 5 ),
  m_rightMargin( 5 ),
  m_bottomMargin( 0 ),
  m_leftMargin( 5 ),
  m_titlePadding( -3 ),
  m_pendingJs{},
  m_renderFlags( 0 ),
  m_ndof( 0 ),
  m_used_points{}
{
  addStyleClass( "Chi2Graphic" );
  setOverflow(Overflow::OverflowHidden);
  
  wApp->useStyleSheet( "InterSpec_resources/Chi2Graphic.css" );
  
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/Chi2Graphic.js" );
  
  m_renderFlags |= Chi2RenderActions::SetData;
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
  if( !m_yAxisTitle.empty() )
    options += ", yAxisTitle: " + m_yAxisTitle.jsStringLiteral('\"');
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


void Chi2Graphic::setData( const int ndof, const std::vector<Chi2Graphic::PeakFitInfo> &points )
{
  m_ndof = ndof;
  m_used_points = points;
  m_renderFlags |= Chi2RenderActions::SetData;
  scheduleRender();
}//void setData(ndof, info)


void Chi2Graphic::setDataToClient()
{
  Wt::Json::Object json;
  
  json["NDOF"] = m_ndof;
  Wt::Json::Array &points = json["Points"] = Wt::Json::Array{};
  
  for( const Chi2Graphic::PeakFitInfo &point : m_used_points )
  {
    Wt::Json::Object jpoint;
    jpoint["energy"] = point.energy;
    jpoint["numSigmaOff"] = point.numSigmaOff;
    jpoint["observedOverExpected"] = point.observedOverExpected;
    jpoint["observedOverExpectedUncert"] = point.observedOverExpectedUncert;
    jpoint["color"] = WString::fromUTF8( point.peakColor.isDefault() ? string("blue") : point.peakColor.cssText() );
    jpoint["nuclide"] = point.nuclidename;
    
    points.push_back( jpoint );
  }//for( loop over m_used_points )
  
  const string js = m_jsgraph + ".setData(" + Wt::Json::serialize(json) + ");";
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setData( std::shared_ptr<Measurement> data_hist )


void Chi2Graphic::setXAxisTitle( const Wt::WString &title )
{
  m_xAxisTitle = title;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setXAxisTitle(" + m_xAxisTitle.jsStringLiteral('\'') + ",true);" );
}//


void Chi2Graphic::setYAxisTitle( const Wt::WString &title )
{
  m_yAxisTitle = title;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setYAxisTitle(" + m_yAxisTitle.jsStringLiteral('\'') + ",true);" );
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
