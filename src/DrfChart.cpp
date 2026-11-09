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


#include <atomic>
#include <memory>
#include <string>
#include <sstream>
#include <iostream>

#include <Wt/Utils.h>
#include <Wt/WColor.h>
#include <Wt/WLength.h>
#include <Wt/WJavaScript.h>
#include <Wt/WApplication.h>
#include <Wt/WWebWidget.h>
#include <Wt/WStringStream.h>
#include <Wt/WCssStyleSheet.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/AppUtils.h"
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
  m_reference( nullptr ),
  m_showAngles( false ),
  m_sourceDistance( 25.0 * PhysicalUnits::cm ),
  m_intrinsic( false ),
  m_showFwhm( true ),
  m_showEfficiency( true ),
  m_jsgraph( jsRef() + ".chart" ),
  m_jsDefined( false ),
  m_xRangeSet( false ),
  m_xRangeMin( 0.0 ),
  m_xRangeMax( 0.0 ),
  m_xRangeSinceUpdate( false ),
  m_xRangeAfterNextData( false ),
  m_detectorSerial( 1 ),
  m_anglesSerial( 1 ),
  m_referenceSerial( 1 ),
  m_detectorJson( "null" ),
  m_anglesJson( "null" ),
  m_referenceJson( "null" ),
  m_detectorJsonSerial( 1 ),
  m_anglesJsonSerial( 1 ),
  m_referenceJsonSerial( 1 ),
  m_detectorRequestedSerial( 1 ),
  m_anglesRequestedSerial( 1 ),
  m_referenceRequestedSerial( 1 ),
  m_computeGeneration( 0 ),
  m_pending( false ),
  m_alive( std::make_shared<std::atomic<bool>>( true ) ),
  m_dataPointsJs(),
  m_showEffPoints( true ),
  m_showFwhmPoints( true ),
  m_keepZoom( false ),
  m_dataRangeLow( 0.0 ),
  m_dataRangeHigh( 0.0 )
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
  // Always a new serial, even for the same pointer: a caller may have changed the DRF it points to.
  m_detector = det;
  ++m_detectorSerial;
  ++m_anglesSerial;
  m_xRangeSinceUpdate = false;
  m_xRangeAfterNextData = false;   //a range restored for the previous detector does not carry over
  scheduleRender();
} //DrfChart::updateChart()


void DrfChart::setShowResponseAngles( const bool show )
{
  m_showAngles = show;
  ++m_anglesSerial;
  scheduleRender();
}//DrfChart::setShowResponseAngles(...)


void DrfChart::setSourceDistance( const double distance )
{
  m_sourceDistance = (distance > 0.0) ? distance : m_sourceDistance;
  ++m_anglesSerial;
  scheduleRender();
}//DrfChart::setSourceDistance(...)


void DrfChart::refreshData()
{
  // Out of date, and not already being computed (a render while a computation runs - the
  //  completion's own push causes one - must not start the same work again).
  const bool need_detector = (m_detectorJsonSerial != m_detectorSerial)
                             && (m_detectorRequestedSerial != m_detectorSerial);
  const bool need_angles = (m_anglesJsonSerial != m_anglesSerial)
                           && (m_anglesRequestedSerial != m_anglesSerial);
  const bool need_reference = (m_referenceJsonSerial != m_referenceSerial)
                              && (m_referenceRequestedSerial != m_referenceSerial);
  if( !need_detector && !need_angles && !need_reference )
    return;

  if( need_detector )
    m_detectorRequestedSerial = m_detectorSerial;
  if( need_angles )
    m_anglesRequestedSerial = m_anglesSerial;
  if( need_reference )
    m_referenceRequestedSerial = m_referenceSerial;

  const shared_ptr<const DetectorPeakResponse> det = (m_detector && m_detector->isValid())
                                                       ? m_detector : nullptr;
  const shared_ptr<const DetectorPeakResponse> ref = (m_reference && m_reference->isValid())
                                                       ? m_reference : nullptr;
  const bool show_angles = m_showAngles;
  const double distance = m_sourceDistance;
  const uint64_t detector_serial = m_detectorSerial, angles_serial = m_anglesSerial;
  const uint64_t reference_serial = m_referenceSerial;
  const uint64_t generation = ++m_computeGeneration;

  struct Result
  {
    string detector_json = "null", angles_json = "null", reference_json = "null";
  };

  const auto compute = [show_angles, distance, need_detector, need_angles, need_reference](
                                const shared_ptr<const DetectorPeakResponse> &drf,
                                const shared_ptr<const DetectorPeakResponse> &ref_drf, Result &result ){
    try
    {
      if( need_detector && drf )
        result.detector_json = drf->toJSON();
      if( need_angles && drf && show_angles )
        result.angles_json = drf->responseAngleSeriesJSON( distance );
      if( need_reference && ref_drf )
        result.reference_json = ref_drf->toJSON();
    }catch( std::exception &e )
    {
      cerr << "DrfChart: failed to compute chart data: " << e.what() << endl;
    }
  };//compute lambda

  // Only a Monte-Carlo / transfer response makes this expensive (it is evaluated a few hundred
  //  times); for anything else, a worker thread would only add a round trip before the curve shows.
  const bool slow = (det && det->ceeloResponse() && (need_detector || need_angles))
                    || (ref && ref->ceeloResponse() && need_reference);
  if( !slow )
  {
    Result result;
    compute( det, ref, result );
    if( need_detector )
    {
      m_detectorJson = std::move( result.detector_json );
      m_detectorJsonSerial = detector_serial;
    }
    if( need_angles )
    {
      m_anglesJson = std::move( result.angles_json );
      m_anglesJsonSerial = angles_serial;
    }
    if( need_reference )
    {
      m_referenceJson = std::move( result.reference_json );
      m_referenceJsonSerial = reference_serial;
    }
    pushData( need_detector, need_angles, need_reference );
    setPending( false );
    return;
  }//if( cheap to compute )

  setPending( true );

  // Only snapshots of the DRFs (const objects, whose queries are thread safe) and plain values
  //  cross to the worker; the chart itself is touched only in the completion, on the session
  //  thread.  Snapshots, because nothing stops a caller changing the DRF it passed to #updateChart
  //  (it then calls #updateChart again); copies are cheap - they share the curves and the response.
  const shared_ptr<const DetectorPeakResponse> snapshot
                                  = det ? make_shared<const DetectorPeakResponse>( *det ) : nullptr;
  const shared_ptr<const DetectorPeakResponse> ref_snapshot
                                  = ref ? make_shared<const DetectorPeakResponse>( *ref ) : nullptr;
  const shared_ptr<Result> result = make_shared<Result>();
  const shared_ptr<std::atomic<bool>> alive = m_alive;
  DrfChart * const self = this;

  AppUtils::run_on_ioservice( [compute, snapshot, ref_snapshot, result](){
      compute( snapshot, ref_snapshot, *result );
    },
    [self, alive, result, need_detector, need_angles, need_reference,
     detector_serial, angles_serial, reference_serial, generation](){
      if( !alive->load() )
        return;

      // A result is still good if what it was computed for has not changed since, even when a
      //  later computation (for something else that did change) is on its way.
      const bool use_detector = need_detector && (detector_serial == self->m_detectorSerial);
      const bool use_angles = need_angles && (angles_serial == self->m_anglesSerial);
      const bool use_reference = need_reference && (reference_serial == self->m_referenceSerial);
      if( use_detector )
      {
        self->m_detectorJson = std::move( result->detector_json );
        self->m_detectorJsonSerial = detector_serial;
      }
      if( use_angles )
      {
        self->m_anglesJson = std::move( result->angles_json );
        self->m_anglesJsonSerial = angles_serial;
      }
      if( use_reference )
      {
        self->m_referenceJson = std::move( result->reference_json );
        self->m_referenceJsonSerial = reference_serial;
      }
      self->pushData( use_detector, use_angles, use_reference );

      if( generation == self->m_computeGeneration )
        self->setPending( false );
    } );
}//void refreshData()


void DrfChart::pushData( const bool detector, const bool angles, const bool reference )
{
  if( detector )
  {
    doChartJs( "setDetectorData(" + m_detectorJson + ")" );

    // The client just reset the x-axis to the detector's range; restore one set to override it.
    if( m_xRangeSet && (m_xRangeSinceUpdate || m_xRangeAfterNextData) )
      doChartJs( "setXRange(" + std::to_string(m_xRangeMin) + ", " + std::to_string(m_xRangeMax) + ")" );
    m_xRangeAfterNextData = false;
  }//if( detector )

  if( angles )
    doChartJs( "setResponseSeries(" + m_anglesJson + ")" );

  if( reference )
    doChartJs( "setReferenceData(" + m_referenceJson + ")" );
}//void pushData(...)


void DrfChart::setPending( const bool pending )
{
  if( pending == m_pending )
    return;

  m_pending = pending;
  doChartJs( "setPending(" + string(pending ? "true" : "false") + ")" );
}//void setPending( const bool pending )


void DrfChart::setIntrinsicEfficiency( const bool intrinsic )
{
  m_intrinsic = intrinsic;
  doChartJs( "setEfficiencyMode('" + string(intrinsic ? "intrinsic" : "absolute") + "')" );
}//DrfChart::setIntrinsicEfficiency(...)


void DrfChart::setDataPoints( const std::vector<DrfChartPoint> &points )
{
  WStringStream js;
  js << "[";
  for( size_t i = 0; i < points.size(); ++i )
  {
    const DrfChartPoint &p = points[i];
    js << (i ? "," : "") << "{e:" << p.energy
       << ",eff:" << p.efficiency << ",effSig:" << p.efficiencyUncert
       << ",fwhm:" << p.fwhm << ",fwhmSig:" << p.fwhmUncert
       << ",label:" << WWebWidget::jsStringLiteral(p.label)
       << ",color:" << WWebWidget::jsStringLiteral(p.color) << "}";
  }
  js << "]";

  m_dataPointsJs = js.str();
  doChartJs( "setDataPoints(" + m_dataPointsJs + ")" );
}//setDataPoints(...)


void DrfChart::setShowEfficiencyPoints( const bool show )
{
  m_showEffPoints = show;
  doChartJs( "setShowEffPoints(" + string(show ? "true" : "false") + ")" );
}//setShowEfficiencyPoints(...)


void DrfChart::setShowFwhmPoints( const bool show )
{
  m_showFwhmPoints = show;
  doChartJs( "setShowFwhmPoints(" + string(show ? "true" : "false") + ")" );
}//setShowFwhmPoints(...)


void DrfChart::setDataRange( const double lowKeV, const double highKeV )
{
  m_dataRangeLow = lowKeV;
  m_dataRangeHigh = highKeV;
  doChartJs( "setDataRange(" + std::to_string(lowKeV) + "," + std::to_string(highKeV) + ")" );
}//setDataRange(...)


void DrfChart::setKeepZoomOnUpdate( const bool keep )
{
  m_keepZoom = keep;
  doChartJs( "setKeepZoom(" + string(keep ? "true" : "false") + ")" );
}//setKeepZoomOnUpdate(...)


void DrfChart::setShowFwhm( const bool show )
{
  m_showFwhm = show;
  doChartJs( "setShowFwhm(" + string(show ? "true" : "false") + ")" );
}//DrfChart::setShowFwhm(...)


void DrfChart::setShowEfficiency( const bool show )
{
  m_showEfficiency = show;
  doChartJs( "setShowEfficiency(" + string(show ? "true" : "false") + ")" );
}//DrfChart::setShowEfficiency(...)


void DrfChart::setReferenceDetector( std::shared_ptr<const DetectorPeakResponse> det )
{
  if( det == m_reference )
    return;
  m_reference = det;
  ++m_referenceSerial;
  scheduleRender();
}//DrfChart::setReferenceDetector(...)


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
  //  say - and the object that is built is empty.  The detector data and angle series are the ones
  //  last computed; if they are out of date, #render computes new ones right after this.
  doChartJs( "setDetectorData(" + m_detectorJson + ")" );
  doChartJs( "setShowFwhm(" + string(m_showFwhm ? "true" : "false") + ")" );
  doChartJs( "setShowEfficiency(" + string(m_showEfficiency ? "true" : "false") + ")" );
  doChartJs( "setReferenceData(" + m_referenceJson + ")" );
  doChartJs( "setEfficiencyMode('" + string(m_intrinsic ? "intrinsic" : "absolute") + "')" );
  doChartJs( "setShowEffPoints(" + string(m_showEffPoints ? "true" : "false") + ")" );
  doChartJs( "setShowFwhmPoints(" + string(m_showFwhmPoints ? "true" : "false") + ")" );
  doChartJs( "setKeepZoom(" + string(m_keepZoom ? "true" : "false") + ")" );
  if( m_dataRangeHigh > m_dataRangeLow )
    doChartJs( "setDataRange(" + std::to_string(m_dataRangeLow) + ","
                               + std::to_string(m_dataRangeHigh) + ")" );
  if( !m_dataPointsJs.empty() )
    doChartJs( "setDataPoints(" + m_dataPointsJs + ")" );
  doChartJs( "setResponseSeries(" + m_anglesJson + ")" );
  if( m_xRangeSet )
    doChartJs( "setXRange(" + std::to_string(m_xRangeMin) + ", " + std::to_string(m_xRangeMax) + ")" );
  m_xRangeAfterNextData = (m_xRangeSet && (m_detectorJsonSerial != m_detectorSerial));
  if( m_pending )
    doChartJs( "setPending(true)" );
}//void DrfChart::defineJavaScript()


void DrfChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = flags.test( Wt::RenderFlag::Full );
  
  WContainerWidget::render( flags );
  
  // `!m_jsDefined` as well as the Full flag: a first render does not always carry it (see
  //  m_jsDefined), and without the client-side object every queued call hits `undefined`.
  if( renderFull || !m_jsDefined )
    defineJavaScript();

  // Here rather than in #updateChart and friends, so a detector pushed several times in one event
  //  (a common pattern) is computed once.
  refreshData();
}//void DrfChart::render(...)


void DrfChart::setXAxisRange( double minEnergy, double maxEnergy )
{
  m_xRangeSet = true;
  m_xRangeSinceUpdate = true;
  m_xRangeMin = minEnergy;
  m_xRangeMax = maxEnergy;

  doChartJs( "setXRange(" + std::to_string(minEnergy) + ", " + std::to_string(maxEnergy) + ")" );
}//void DrfChart::setXAxisRange(...)



DrfChart::~DrfChart()
{
  m_alive->store( false );

  Wt::WApplication *app = Wt::WApplication::instance();
  if( app )
  {
    WCssStyleSheet &style = app->styleSheet();
    for( const auto &rule : m_cssRules )
      style.removeRule( rule.second );
    m_cssRules.clear();
  }//if( app )
}//~DrfChart()

