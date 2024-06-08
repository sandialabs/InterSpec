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
 - [x] Every update seems to trigger a layout resize that causes the chart to grow.
 - [ ] The MdaPeakRow (i.e., each gamma) should probably be like a WPanel that has just the "Use for Likilihood" checkbox in the title, as well as the currie MDA summary text, but then can be expanded to show the other settings
 - [ ] Could allow entering a scale factor for spectrum; this is like you have a 5 minute background, but are interested in a 60s dwell
 - [ ] Have a little pop-up window or something that gives full information for Curie calculation
 - [ ] Add in a Curie style MDA calculator where you enter limits, BR, and stuff and get out counts and such
 - [ ] The likelihood based estimate does not seem to be reliable yet
 - [x] Add in checkbox to allow accounting for attenuation in the air (or just always do this)
 - [x] Add in allowing to calculate the maximum detection distance; not sure if this should be individual peaks, or all peaks; maybe a toggle for the whole screen
 - [x] Add a drop-box to allow selecting confidence for limit (e.g., 80%, 90%, 95%, 99%, etc)
 - [ ] Allow asserting you know there isnt a peak in the spectrum, and affix the peak continuums to the observed spectrum
 - [x] Allow adjusting the simple MDA side-widths
 - [ ] Allow fitting for energy calibration, or accounting for it anyway
 - [ ] If spectrum already has a peak defined for a gamma, re-use its definition of roi lower and upper energies
 - [ ] When deconvolution method is used along with continuum determined from channels above/below the ROI, need to probably take into account the stat uncertainty of the channels...
 - [ ] Check numerical accuracies of calculations
 - [ ] Report mass as well as activity
 - [ ] Make sure when presence is reported on Currie style limit (e.g., both lower and upper limits given), the coverage is actually correct; e.g. 95% of times in given interval, and not 97.5% or 90% or something.
 - [ ] Allow volumetric sources; either trace or self-attenuating.
 - [ ] Add in pass-by (e.g., at a fixed speed/distance) calculation
  - [ ] Have C-api to call from node/python
 
 - [x] Switch to using  a CSS grid layout for the options section
 - [x] put in BR (after DRF and shield) to limit total number of peaks; have it be zero by default if less than 10 peaks, or 0.1 or something otherwise
 - probably also filter lines that have essentually a zero efficiency per bq
 - [ ] Make a "by eye" equivalent CL
 
 - [ ] Make test cases that will quickly iterate through, to test things
 - [ ] Give the user a choice about using continuum fixed at null hypothesis
 - [ ] For each row show plot of current peak in that row
 - [ ] Allow combining ROI with neighboring peaks
    - Actually should just include all gammas in an energy range, and have them all share a continuum - so a user could widen the ROI to include multiple peaks, and all would be fine
 - [ ] In addition to error messages, have an area for warning messages (like if scaling spectrum more than 1.0, etc)
 - [ ] Allow user to pick Currie limit ranges, and improve clarity of this stuff, like maybe have each peak be a WPanel or something
 - [ ] Have the energy rows fold down to show more information, similar to Steves tool, for each energy
 - [ ] If user clicks on a result row, have chart zoom to that general region
 
 - [ ] Allow minor gamma-lines overlapping with primary gamma lines to contribute to peak area
 - [ ] Allow users to double click on the spectrum to add a peak to the limit, or similarly for erasing a peak
 - [ ] Remember ROI properties for all user changed ROIs, for full use, not just if nuclide changes
 
 - [x] Make the Chi2 plot a D3 based plot
 - [x] Put the Chi2 chart to the right of the spectrum, when it should exist
 - [x] Make it so when user change ROI on chart, it updates the text input
 - [x] Add in allowing to age nuclide (didnt I generalize inputting a nuclide somewhere?  Hopefully just re-use that)
 - [x] Default fill in reference lines/shielding/age as user has in Reference PhotoPeak tool
 */


class PeakDef;
class PeakModel;
class InterSpec;
class MaterialDB;
class MdaPeakRow;
class SwitchCheckbox;
class ShieldingSelect;
class DetectorDisplay;
class NativeFloatSpinBox;
class DetectionLimitTool;
class D3SpectrumDisplayDiv;
class DetectorPeakResponse;

namespace DetectionLimitCalc
{
  struct CurieMdaInput;
  struct CurieMdaResult;
  struct DeconComputeInput;
  enum class DeconContinuumNorm : int;
}

namespace SpecUtils
{
  class Measurement;
};//namespace SpecUtils

namespace SandiaDecay
{
  struct Nuclide;
};//namespace SandiaDecay

namespace Wt
{
  class WText;
  class WLabel;
  class WCheckBox;
  class WLineEdit;
  class WComboBox;
  class WPushButton;
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
  
  DetectionLimitTool *tool();
protected:
  DetectionLimitTool *m_tool;
};//class DetectionLimitWindow


class DetectionLimitTool : public Wt::WContainerWidget
{
public:
  /** The limits type the user can select to determine. */
  enum class LimitType
  {
    Activity,
    Distance
  };
  
public:
  DetectionLimitTool( InterSpec *viewer,
                          MaterialDB *materialDB,
                          Wt::WSuggestionPopup *materialSuggest,
                          Wt::WContainerWidget *parent = 0 );
  
  virtual ~DetectionLimitTool();
  
  void do_development();
  
  void setRefLinesAndGetLineInfo();
  
  void scheduleCalcUpdate();
  
  
  /** Handles receiving a "deep-link" url starting with "interspec://detection-limit/...".
   
   TODO: NOT CURRENTLY IMPLEMENTED
   
   Example URIs:
   - "interspec://detection-limit?V=1&..."
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
   "interspec://detection-limit?V=1...", then this string would be "V=1...".
   Assumes the string passed in has alaready been url-decoded.
   If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string does not include the
   "interspec://" protocol, or "detection-limit" path; so will look something like "V=1&...",
   and it will not be url-encoded.
   
   TODO: NOT CURRENTLY IMPLEMENTED - returns empty string
   */
  std::string encodeStateToUrl() const;
  
  /** Updates the decorative region, visible region, and text for a Currie limit result.
   
   @param chart Chart to update; must not be nullptr.
   @param PeakModel Peak model to add resulting peak to; may be nullptr if you dont want a peak on chart.
   @param input The Currie MDA input; must be filled out.
   @param result The Currie limit results.  May be nullptr if it hasnt/couldnt be computed (in which case no peak/text will be drawn)
   @param drf Only use for peak FWHM; if nullptr or no FWHM info, will use 0.25 ROI for peak sigma width.
   @param limitType The limit type - currently only activity limit is supported
   @param gammas_per_bq The expected detected gammas per Bq.  If zero or negative, will put text in terms of counts
   */
  static void update_spectrum_for_currie_result( D3SpectrumDisplayDiv *chart,
                                         PeakModel *pmodel,
                                         const DetectionLimitCalc::CurieMdaInput &input,
                                         const DetectionLimitCalc::CurieMdaResult * const result,
                                         std::shared_ptr<const DetectorPeakResponse> drf,
                                         DetectionLimitTool::LimitType limitType,
                                         const double gammas_per_bq );
  
  static void createCurrieRoiMoreInfoWindow( const SandiaDecay::Nuclide *const nuclide,
                                  const DetectionLimitCalc::CurieMdaResult &result,
                                  std::shared_ptr<const DetectorPeakResponse> drf,
                                  DetectionLimitTool::LimitType limitType,
                                            const double distance,
                                            const bool do_air_attenuation,
                                            const double branch_ratio,
                                            const double counts_per_bq_into_4pi  );
  
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
  
  /** Returns the input for computing Chi2 
   
   Returns nullptr if cant prepare input (e.g., no peaks, do displayed histogram, etc).
   
   You can simple change `DetectionLimitCalc::DeconComputeInput::activity` and
   `DetectionLimitCalc::DeconComputeInput::distance` of the returned object, and the
   rest of the fields will remain valid (e.g., for scanning over activity or distance, you only need to
   prepare the input once).
   */
  std::shared_ptr<DetectionLimitCalc::DeconComputeInput> getComputeForActivityInput(
                                                                    const double activity,
                                                                    const double distance,
                                                                    std::vector<PeakDef> &peaks );
  
  void computeForActivity( const double activity,
                          const double distance,
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
  
  void handleUserChangedConfidenceLevel();
  
  void handleUserChangedUseAirAttenuate();
  
  void handleUserChangedToComputeActOrDist();
  
  void handleDrfSelected( std::shared_ptr<DetectorPeakResponse> new_drf );
  
  void handleFitFwhmRequested();
  
  void handleInputChange();
  
  /** Gets DRF from GUI widget, and if that isnt valid, gets it from the spectrum. */
  std::shared_ptr<const DetectorPeakResponse> detector();
  
  float currentConfidenceLevel();
  
  /** Returns either the current user entered distance (if determining activity limit), or the current display distance (if determining
   distance limit).
   
   Throws exception if invalid or zero distance.
   */
  double currentDisplayDistance() const;
  
  /** Returns either the current user entered activity (if determining activity limit), or the current display activity (if determining
   distance limit).
   
   Throws exception if invalid or zero distance.
   */
  double currentDisplayActivity() const;
  
  /** Returns the type of limit, either activity or distance, the user currently has selected t determine. */
  LimitType limitType() const;
  
  /** A struct to hold a little information about gamma lines for a nuclide. */
  struct GammaLineInfo
  {
    /** Energy of the gamma line this information is for. */
    double energy;
    
    /** The branching ratio of the gamma, per decay of the parent isotope. */
    double branching_ratio;
    
    /** Accounts for B.R. after aging and attenuation through shielding, but not attenuation in the air. */
    double gammas_into_4pi;
    
    /** Adds in attenuation in the air, at the currently specified distance (if distance input is invalid, will use 1m) */
    double gammas_4pi_after_air_attenuation;
    
    /** The fraction of gammas that make it through the air without interacting; e.g. near 1.0 for short distances. */
    double air_transmission;
    
    /** The fraction of gammas that make it through the specified shielding without interacting. Will 1.0 if no shielding. */
    double shield_transmission;
  };//struct GammaLineInfo
  
  /** Returns the current gamma lines information.

   May throw exception if error calculating attenuation, or no current valid nuclide.
   */
  std::vector<GammaLineInfo> gammaLines() const;
  
  
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
  
  /** The `DetectorDisplay` only holds a weak pointer to the selected DRF; and when the user selects
   the DRF in this widget, we dont assign it to the foreground, so we will keep the selected detector alive
   using this shared pointer.
   */
  std::shared_ptr<DetectorPeakResponse> m_drf_cache;
  
  /** The switch to toggle if the user wants to compute activity limit, or distance limit.
   m_distOrActivity->isChecked() indicates compute distance.
   */
  SwitchCheckbox *m_distOrActivity;
  Wt::WLabel *m_activityLabel;
  Wt::WLabel *m_distanceLabel;
  
  Wt::WLineEdit *m_distanceForActivityLimit;
  
  Wt::WLineEdit *m_activityForDistanceLimit;
  
  MaterialDB *m_materialDB;                 //not owned by this object
  Wt::WSuggestionPopup *m_materialSuggest;  //not owned by this object
  ShieldingSelect *m_shieldingSelect;
  
  NativeFloatSpinBox *m_minRelIntensity;
  
  Wt::WCheckBox *m_attenuateForAir;
  
  Wt::WLabel *m_displayActivityLabel;
  Wt::WLabel *m_displayDistanceLabel;
  Wt::WLineEdit *m_displayActivity; //!< Used for peak plotting when activity limit is being determined
  Wt::WLineEdit *m_displayDistance; //!< Used for peak plotting when distance limit is being determined
  
  enum ConfidenceLevel { OneSigma, TwoSigma, ThreeSigma, FourSigma, FiveSigma, NumConfidenceLevel };
  Wt::WComboBox *m_confidenceLevel;
  
  /** Holds m_chi2Chart, m_bestChi2Act, and m_upperLimit. */
  Wt::WContainerWidget *m_results;
  
  /** Holds the D3.js based MdaChi2Chart */
  Wt::WContainerWidget *m_chi2Chart;
  Wt::WText *m_bestChi2Act;
  Wt::WText *m_upperLimit;
  
  Wt::WText *m_errorMsg;
  
  Wt::WPushButton *m_fitFwhmBtn;
  
  std::shared_ptr<SpecMeas> m_our_meas;
  Wt::WContainerWidget *m_peaks;
  
public:
  
  // We want to preserve user changes to ROIs across nuclide changes, so we'll cache all
  //  previous values every time we update.
  //
  // This struct could/should maybe be composed into parts that can be changed by the MdaPeakRow
  //  widget, vs outside of it, vs data, but for the development (and the moment), we'll just
  //  cram it into one struct until it becomes clear the best way to organize things, or even if
  //  it matters.
  struct MdaPeakRowInput
  {
    bool use_for_likelihood;
    LimitType limit_type;
    bool do_air_attenuation;
    
    float energy;
    double branch_ratio; //!< Decay branching ratio of nuclide
    double counts_per_bq_into_4pi; //!< This includes spectrum live time, shielding, and gamma BR.  This will need to be turned into some vector to cover multiple gamma lines
    double counts_per_bq_into_4pi_with_air;
    double trans_through_air;
    double trans_through_shielding;
    double distance;
    double activity;
    float roi_start;
    float roi_end;
    size_t num_side_channels;
    float confidence_level;
    
    /** How we should normalize the continuum.
     E.g., allowed to float as normal, fixed using the edges of the continuum, or fixed using the whole ROI, assuming no signal
     (the gaussian component will always just sit on top of the continuum, no matter the activity, for the later two options).
     */
    DetectionLimitCalc::DeconContinuumNorm decon_cont_norm_method;
    
    /** The continuum type to use for the deconvolution limit type.
     Currently only linear or quadratic, but maybe we can expand this.
     */
    PeakContinuum::OffsetType decon_continuum_type;
    
    std::shared_ptr<const DetectorPeakResponse> drf;
    std::shared_ptr<const SpecUtils::Measurement> measurement;
  };//struct MdaPeakRowInput
  
protected:
  std::map<float,MdaPeakRowInput> m_previousRoiValues;
};//class DetectionLimitTool





#endif //DetectionLimitTool_h
