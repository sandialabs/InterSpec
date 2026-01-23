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


#include <memory>
#include <sstream>
#include <string>

#include <Wt/Utils>
#include <Wt/WColor>
#include <Wt/WLength>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WCssStyleSheet>
#include <Wt/WContainerWidget>

#include "InterSpec/DrfChart.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;


DrfChart::DrfChart( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_detector( nullptr ),
  m_jsgraph( jsRef() + ".chart" )
{
  addStyleClass( "DrfChart" );
  setOverflow( Overflow::OverflowHidden );
  
  // Require JavaScript resources
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/DetectorPeakResponseJS.js" );
  wApp->require( "InterSpec_resources/DrfChart.js" );
  
  wApp->useStyleSheet( "InterSpec_resources/DrfChart.css" );
}//DrfChart constructor


void DrfChart::updateChart( std::shared_ptr<const DetectorPeakResponse> det )
{
  m_detector = det;
  
  // Generate and send detector data using DetectorPeakResponse JSON generation
  string detectorData = (!det || !det->isValid()) ? string("null") : det->toJSON();
  const string detectorJs = m_jsgraph + ".setDetectorData(" + detectorData + ");";
  if( isRendered() )
    doJavaScript( detectorJs );
  else
    m_pendingJs.push_back( detectorJs );
} //DrfChart::updateChart()


void DrfChart::defineJavaScript()
{
  string options = "{ "
    "margins: {"
    " top: 5,"
    " right: 50,"
    " bottom: 40,"
    " left: 60"
    " } }";
  
  setJavaScriptMember( "chart", "new DrfChart(" + jsRef() + ", " + options + ");");
  
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
  
  // Execute any pending JS calls
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void DrfChart::defineJavaScript()


void DrfChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
}//void DrfChart::render(...)


void DrfChart::setXAxisRange( double minEnergy, double maxEnergy )
{
  const string js = m_jsgraph + ".setXRange(" + std::to_string(minEnergy) + ", " + std::to_string(maxEnergy) + ");";
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void DrfChart::setXAxisRange(...)



DrfChart::~DrfChart()
{
  Wt::WApplication *app = Wt::WApplication::instance();
  if( app )
  {
    WCssStyleSheet &style = app->styleSheet();
    for( const auto &rule : m_cssRules )
      style.removeRule( rule.second );
    m_cssRules.clear();
  }//if( app )
}//~DrfChart()

