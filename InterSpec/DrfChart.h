#ifndef DrfChart_h
#define DrfChart_h
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
#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <cstdint>

#include <Wt/WContainerWidget.h>

struct ColorTheme;
class DetectorPeakResponse;

namespace Wt
{
  class WColor;
  class WCssTextRule;
}

/** A measured data point drawn on a DrfChart: an efficiency marker with error bar on the left
 axis, and/or a FWHM marker on the right axis (a value <= 0 draws nothing on that axis). */
struct DrfChartPoint
{
  double energy = 0.0;               //keV
  double efficiency = 0.0, efficiencyUncert = 0.0;
  double fwhm = 0.0, fwhmUncert = 0.0;
  std::string label;                 //tooltip text, e.g. "Cs137 661.7 keV: 12345 counts"
  std::string color;                 //CSS color, or empty for the default
};//struct DrfChartPoint


/** A D3 chart of a DRF's efficiency (with its uncertainty band) and FWHM, optionally with per-angle
 response curves.

 The chart data is not computed when it is asked for, but once per render (see #render): callers
 commonly push the same detector two or three times in one event.  For a DRF with a Monte-Carlo /
 transfer response, computing it means evaluating the response a few hundred times, so that is done
 on a worker thread, and the chart dims until the result arrives; a DRF without one is cheap and is
 computed in place.
 */
class DrfChart : public Wt::WContainerWidget
{
protected:
  void defineJavaScript();
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );

  std::shared_ptr<const DetectorPeakResponse> m_detector;

  /** Whether to draw the per-angle response curves (only meaningful when the
   detector has a Monte-Carlo / transfer parameterized response). */
  bool m_showAngles;

  /** Source distance (PhysicalUnits) the absolute per-angle curves are drawn
   at; also where the near-field validity flag is evaluated. */
  double m_sourceDistance;

  /** True = show intrinsic (per-gamma-striking-face) efficiency; false =
   absolute (per-emitted-gamma) at #m_sourceDistance. */
  bool m_intrinsic;

  /** Whether the DRFs FWHM curve (and its right-hand axis) is drawn. */
  bool m_showFwhm;

  /** The javascript variable name used to refer to the DrfChart object.
   Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;
  
  /** Whether #defineJavaScript has run at least once. */
  bool m_jsDefined;

  /** The x-axis range #setXAxisRange was last given, so #defineJavaScript can restore it. */
  bool m_xRangeSet;
  double m_xRangeMin, m_xRangeMax;

  /** Whether #setXAxisRange has been called since the last #updateChart.  The client resets the
   x-axis to the detector's energy range whenever detector data arrives, so a range set after the
   detector it is meant to override is sent again once that detector's data has been pushed.
   */
  bool m_xRangeSinceUpdate;

  /** Set when #defineJavaScript restores the x-axis range, but the detector data it sent is not yet
   up to date: the range has to survive that data's arrival too, as it did when the data was
   computed in place.  Cleared by the next detector-data push.
   */
  bool m_xRangeAfterNextData;

  /** Bumped by every change to what the chart data depends on: #m_detectorSerial by #updateChart;
   #m_anglesSerial by that and by the angle settings.  #render computes data for the current serials
   (once, however many changes there were), and a result is kept only if it was computed for the
   serials still current when it is ready.
   */
  uint64_t m_detectorSerial, m_anglesSerial;

  /** The client-side detector data (`DetectorPeakResponse::toJSON()`) and per-angle series last
   computed, and the serials they were computed for - so a full re-render re-sends them rather than
   computing them again.
   */
  std::string m_detectorJson, m_anglesJson;
  uint64_t m_detectorJsonSerial, m_anglesJsonSerial;

  /** The serials the most recent computation was started for, so a render while it is still running
   does not start the same work again. */
  uint64_t m_detectorRequestedSerial, m_anglesRequestedSerial;

  /** Bumped by every computation #render starts; the one with the latest value clears the chart's
   "computing" state when it lands. */
  uint64_t m_computeGeneration;

  /** Whether the chart is waiting on a background computation. */
  bool m_pending;

  /** Cleared by the destructor.  A background computation's completion runs on the session thread,
   as does the destructor, and checks this before touching the chart - rather than looking the
   chart up by id, which misses one sitting on a tab that has not been shown.
   */
  std::shared_ptr<std::atomic<bool>> m_alive;

  /** The measured-point overlay and its options, held for the same reason as #m_xRangeSet: the
   client-side object is rebuilt empty whenever a stub becomes a real element, so everything the
   chart is showing has to be re-sendable.  #m_dataPointsJs is the serialized array literal
   #setDataPoints last built. */
  std::string m_dataPointsJs;
  bool m_showEffPoints, m_showFwhmPoints;
  bool m_keepZoom;
  double m_dataRangeLow, m_dataRangeHigh;
  
  std::map<std::string,Wt::WCssTextRule *> m_cssRules;

  
public:
  DrfChart();
  virtual ~DrfChart();
  
  void updateChart( std::shared_ptr<const DetectorPeakResponse> det );

  /** Set the x-axis range from C++ */
  void setXAxisRange( double minEnergy, double maxEnergy );

  /** Enables/disables the per-angle response curves (0, 22.5, 45, 62.5, 90
   degrees) plus the angle-flat far-field reference.  No effect for detectors
   without a Monte-Carlo / transfer parameterized response. */
  void setShowResponseAngles( const bool show );

  /** Source distance (in PhysicalUnits) the absolute per-angle curves are
   evaluated at; the angle series is re-sampled for it at the next render. */
  void setSourceDistance( const double distance );

  /** Selects intrinsic (per-gamma-striking-face) vs absolute
   (per-emitted-gamma) efficiency for the per-angle curves. */
  void setIntrinsicEfficiency( const bool intrinsic );

  /** Enables/disables drawing the DRFs FWHM curve and its right-hand axis;
   defaults to enabled.  Turned off where the chart is previewing a geometry /
   efficiency response, which the FWHM has nothing to say about. */
  void setShowFwhm( const bool show );

  /** Sets (or, with an empty vector, clears) measured data points to draw over the curves - the
   efficiency points the Create DRF tool fits, with error bars, and their FWHM markers. */
  void setDataPoints( const std::vector<DrfChartPoint> &points );

  /** Show/hide the efficiency markers, and the FWHM markers, of #setDataPoints. */
  void setShowEfficiencyPoints( const bool show );
  void setShowFwhmPoints( const bool show );

  /** Shades the energy regions outside [lowKeV, highKeV] - where the curve is extrapolated beyond
   the data.  A range with high <= low clears the shading. */
  void setDataRange( const double lowKeV, const double highKeV );

  /** When true, #updateChart keeps the current x-axis zoom rather than resetting it to the new
   detectors energy range - for a chart whose detector is re-pushed after every re-fit. */
  void setKeepZoomOnUpdate( const bool keep );

protected:
  /** Computes whatever chart data is out of date (see #m_detectorSerial) and pushes it to the client
   - in place, or on a worker thread for a DRF with a Monte-Carlo / transfer response.
   */
  void refreshData();

  /** Pushes #m_detectorJson and/or #m_anglesJson to the client. */
  void pushData( const bool detector, const bool angles );

  /** Shows (or clears) the client's "computing" state. */
  void setPending( const bool pending );

  /** Calls one method on the client-side chart object, e.g. `setShowFwhm(true)`.

   The emitted statement checks that the element and its chart object exist before calling, because
   neither is guaranteed: a widget sitting on a tab that has never been shown is sent to the client
   as a *stub*, so its element is simply not in the DOM.  Emitting a bare call in that state threw a
   TypeError which aborted the rest of Wt's update block - taking the surrounding dialog's event
   wiring with it.  #defineJavaScript re-sends the whole state once the element does exist, so
   nothing is lost by a call that no-ops.
   */
  void doChartJs( const std::string &method_call );
};//class DrfChart

#endif //DrfChart_h

