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

#include <Wt/Utils.h>
#include <Wt/WColor.h>
#include <Wt/WLength.h>
#include <Wt/WJavaScript.h>
#include <Wt/WApplication.h>
#include <Wt/WStringStream.h>
#include <Wt/WCssStyleSheet.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/DrfChart.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;


DrfChart::DrfChart()
: WContainerWidget(),
  m_detector( nullptr ),
  m_showAngles( false ),
  m_sourceDistance( 25.0 * PhysicalUnits::cm ),
  m_intrinsic( false ),
  m_showFwhm( true ),
  m_jsgraph( jsRef() + ".chart" ),
  m_jsDefined( false ),
  m_xRangeSet( false ),
  m_xRangeMin( 0.0 ),
  m_xRangeMax( 0.0 )
{
  addStyleClass( "DrfChart" );
  setOverflow( Overflow::Hidden );
  
  // Require JavaScript resources
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/DetectorPeakResponseJS.js" );
  wApp->require( "InterSpec_resources/DrfChart.js" );
  
  wApp->useStyleSheet( "InterSpec_resources/DrfChart.css" );
}//DrfChart constructor


void DrfChart::doChartJs( const std::string &method_call )
{
  // See the note on this function in DrfChart.h: the element and the chart object may both be
  //  absent, so the call has to check before it runs.
  doJavaScript( "{const c=" + jsRef() + ";if(c&&c.chart){c.chart." + method_call + ";}}" );
}//DrfChart::doChartJs(...)


void DrfChart::updateChart( std::shared_ptr<const DetectorPeakResponse> det )
{
  m_detector = det;

  // Generate and send detector data using DetectorPeakResponse JSON generation
  const string detectorData = (!det || !det->isValid()) ? string("null") : det->toJSON();
  doChartJs( "setDetectorData(" + detectorData + ")" );

  pushAngleSeries();
} //DrfChart::updateChart()


void DrfChart::pushAngleSeries()
{
  const string series = (m_showAngles && m_detector && m_detector->isValid())
                          ? m_detector->responseAngleSeriesJSON( m_sourceDistance )
                          : string("null");
  doChartJs( "setResponseSeries(" + series + ")" );
}//DrfChart::pushAngleSeries()


void DrfChart::setShowResponseAngles( const bool show )
{
  m_showAngles = show;
  pushAngleSeries();
}//DrfChart::setShowResponseAngles(...)


void DrfChart::setSourceDistance( const double distance )
{
  m_sourceDistance = (distance > 0.0) ? distance : m_sourceDistance;
  pushAngleSeries();
}//DrfChart::setSourceDistance(...)


void DrfChart::setIntrinsicEfficiency( const bool intrinsic )
{
  m_intrinsic = intrinsic;
  doChartJs( "setEfficiencyMode('" + string(intrinsic ? "intrinsic" : "absolute") + "')" );
}//DrfChart::setIntrinsicEfficiency(...)


void DrfChart::setShowFwhm( const bool show )
{
  m_showFwhm = show;
  doChartJs( "setShowFwhm(" + string(show ? "true" : "false") + ")" );
}//DrfChart::setShowFwhm(...)


void DrfChart::defineJavaScript()
{
  m_jsDefined = true;

  string options = "{ "
    "margins: {"
    " top: 5,"
    " right: 50,"
    " bottom: 40,"
    " left: 3"
    " } }";
  
  // Guarded: `Wt4_x_y.$(id)` is a getElementById, which is null while this widget is only a stub
  //  on a tab that has not been shown.  Wt re-applies JavaScript members (and issues a full render)
  //  when the real element replaces the stub, so the object gets built then instead.
  setJavaScriptMember( "chart",
      "(function(){const e=" + jsRef() + ";return e ? new DrfChart(e, " + options + ") : null;})()" );
  
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

  // Re-send everything the chart is currently showing, rather than replaying a one-shot queue: this
  //  runs again whenever the client-side object is (re)built - when a stub becomes a real element,
  //  say - and the object that is built is empty.
  const string detectorData = (!m_detector || !m_detector->isValid())
                                ? string("null") : m_detector->toJSON();
  doChartJs( "setDetectorData(" + detectorData + ")" );
  doChartJs( "setShowFwhm(" + string(m_showFwhm ? "true" : "false") + ")" );
  doChartJs( "setEfficiencyMode('" + string(m_intrinsic ? "intrinsic" : "absolute") + "')" );
  pushAngleSeries();
  if( m_xRangeSet )
    doChartJs( "setXRange(" + std::to_string(m_xRangeMin) + ", " + std::to_string(m_xRangeMax) + ")" );
}//void DrfChart::defineJavaScript()


void DrfChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = flags.test( Wt::RenderFlag::Full );
  
  WContainerWidget::render( flags );
  
  // `!m_jsDefined` as well as the Full flag: a first render does not always carry it (see
  //  m_jsDefined), and without the client-side object every queued call hits `undefined`.
  if( renderFull || !m_jsDefined )
    defineJavaScript();
}//void DrfChart::render(...)


void DrfChart::setXAxisRange( double minEnergy, double maxEnergy )
{
  m_xRangeSet = true;
  m_xRangeMin = minEnergy;
  m_xRangeMax = maxEnergy;

  doChartJs( "setXRange(" + std::to_string(minEnergy) + ", " + std::to_string(maxEnergy) + ")" );
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

