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

#include <memory>

#include <map>
#include <vector>

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
  
  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
   */
  std::vector<std::string> m_pendingJs;
  
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
   evaluated at; re-samples and pushes the angle series. */
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
  /** Re-samples the per-angle series (at #m_sourceDistance) and pushes it to
   the client, or clears it when angle curves are disabled / unavailable. */
  void pushAngleSeries();
};//class DrfChart

#endif //DrfChart_h

