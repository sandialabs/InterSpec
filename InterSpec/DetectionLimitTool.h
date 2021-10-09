#ifndef DetectionLimitTool_h
#define DetectionLimitTool_h
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

#include <tuple>
#include <vector>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

/** TODO:
 - Every update seems to trigger a layout resize that causes the chart to grow.
 - Could allow entering a scale factor for spectrum; this is like you have a 5 minute background, but are interested in a 60s dwell
 - The likelihood based estimate does not seem to be reliable yet
 - Add in checkbox to allow accounting for attenuation in the air (or just always do this)
 - Add in option to scale time period results are for.  Like if the spectrum is a five-minute spectrum, but you want a 1 second detection limit.
 - Add in allowing to calculate the maximum detection distance; not sure if this should be individual peaks, or all peaks; maybe a toggle for the whole screen
 - Add a drop-box to allow selecting confidence for limit (e.g., 80%, 90%, 95%, 99%, etc)
 - Allow asserting you know there isnt a peak in the spectrum, and affix the peak continuums to the observed spectrum
 - Allow adjusting the simple MDA side-widths
 - Check numerical accuracies of calculations
 */


class PeakDef;
class PeakModel;
class InterSpec;
class MaterialDB;
class MdaPeakRow;
class ShieldingSelect;
class DetectorDisplay;
class NativeFloatSpinBox;
class D3SpectrumDisplayDiv;
class DetectionLimitTool;

namespace SandiaDecay
{
  struct Nuclide;
};//namespace SandiaDecay

namespace Wt
{
  class WText;
  class WLineEdit;
  class WSuggestionPopup;
  class WStandardItemModel;
  namespace Chart
  {
    class WCartesianChart;
  }
}//namespace Wt

class DetectionLimitWindow : public AuxWindow
{
public:
  DetectionLimitWindow( InterSpec *viewer,
                             MaterialDB *materialDB,
                             Wt::WSuggestionPopup *materialSuggest );
  virtual ~DetectionLimitWindow();
  
protected:
  DetectionLimitTool *m_tool;
};//class DetectionLimitWindow


class DetectionLimitTool : public Wt::WContainerWidget
{
public:
  DetectionLimitTool( InterSpec *viewer,
                          MaterialDB *materialDB,
                          Wt::WSuggestionPopup *materialSuggest,
                          Wt::WContainerWidget *parent = 0 );
  
  virtual ~DetectionLimitTool();
  
  void do_development();
  
  void setRefLinesAndGetLineInfo();
  
  void scheduleCalcUpdate();
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  void initChi2Chart();
  
  void roiDraggedCallback( double new_roi_lower_energy,
                   double new_roi_upper_energy,
                   double new_roi_lower_px,
                   double new_roi_upper_px,
                   double original_roi_lower_energy,
                   bool is_final_range );
  
  void doCalc();
  void updateShownPeaks();
  void computeForAcivity( const double activity,
                          std::vector<PeakDef> &peaks,
                          double &chi2, int &numDOF );
  
  void handleUserAgeChange();
  void handleUserNuclideChange();
  void handleNuclideChange( const bool update_to_default_age );
  
  /** Looks at the current nuclide/DRF/shielding, and finds the minimum relative intensity for there to be 20 or less gamma line entries.
   
   Does not hide GUI widgets, or schedule calculation update or anything.
   */
  void calcAndSetDefaultMinRelativeIntensity();
  
  void handleUserMinRelativeIntensityChange();
  void updateShownGammaLinesFromMinIntensity();
  
  
  void handleInputChange();
  
  
  /** Returns the current {energy, observed_br}.
   Where observed_br accounts for B.R. after aging and attenuation through shielding, and detection efficiency at the given distance; e.g.
   how many counts per second you would expect to detect per Bq of activity.
   
   May throw exception if error calculating attenuation, or there is no DRF, or no current valid nuclide.
   */
  std::vector<std::tuple<double,double>> gammaLines() const;
  
  
  InterSpec *m_interspec;
  
  /** When inputs change will mark the widget that it needs update, and delay computatations until
   Wt calls render(WFlags).
   */
  bool m_needsUpdate;
  
  D3SpectrumDisplayDiv *m_chart;
  PeakModel *m_peakModel;
  
  Wt::WLineEdit *m_nuclideEdit;
  Wt::WLineEdit *m_ageEdit;
  
  const SandiaDecay::Nuclide *m_currentNuclide;
  double m_currentAge;
  Wt::WSuggestionPopup *m_nuclideSuggest;
  DetectorDisplay *m_detectorDisplay;
  
  Wt::WLineEdit *m_distanceEdit;
  
  MaterialDB *m_materialDB;                 //not owned by this object
  Wt::WSuggestionPopup *m_materialSuggest;  //not owned by this object
  ShieldingSelect *m_shieldingSelect;
  
  NativeFloatSpinBox *m_minRelIntensity;
  
  
  Wt::WLineEdit *m_displayActivity;
  
  /** Holds m_chi2Chart, m_bestChi2Act, and m_upperLimit. */
  Wt::WContainerWidget *m_results;
  
  /** Holds the D3.js based MdaChi2Chart */
  Wt::WContainerWidget *m_chi2Chart;
  Wt::WText *m_bestChi2Act;
  Wt::WText *m_upperLimit;
  
  Wt::WText *m_errorMsg;
  
  std::shared_ptr<SpecMeas> m_our_meas;
  Wt::WContainerWidget *m_peaks;
  
  
  // We want to preserve user changes to ROIs across nuclide changes, so we'll cache all
  //  previous values every time we update.
  struct PreviousRoiValue
  {
    float energy;
    bool use_for_likelihood;
    float roi_start;
    float roi_end;
    
    PreviousRoiValue();
    PreviousRoiValue( const MdaPeakRow * const row );
  };//struct PreviousRoiValue
  
  std::map<float,PreviousRoiValue> m_previousRoiValues;
};//class DetectionLimitTool





#endif //DetectionLimitTool_h
