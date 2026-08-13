#ifndef DetectionLimitSimple_h
#define DetectionLimitSimple_h
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

#include <string>
#include <vector>
#include <memory>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

class DetectorDisplay;
class NuclideSourceEnter;
class DetectionLimitSimple;
class DetectorPeakResponse;
class NuclideSourceEnterController;

namespace Wt
{
  class WMenu;
  class WText;
  class WCheckBox;
  class WComboBox;
  class WLineEdit;
  class WTabWidget;
  class WPushButton;
  class WGridLayout;
  class WButtonGroup;
  class WRadioButton;
  class WStackedWidget;
  class WSuggestionPopup;
}//namespace Wt

namespace SandiaDecay
{
  class Nuclide;
  struct Transition;
}

namespace SpecUtils
{
  class Measurement;
  enum class SpectrumType : int;
}


namespace DetectionLimitCalc
{
  struct CurrieMdaInput;
  struct CurrieMdaResult;
  struct DeconComputeInput;
  struct DeconActivityOrDistanceLimitResult;
}//namespace DetectionLimitCalc

class DetectionLimitSimpleWindow : public AuxWindow
{
public:
  DetectionLimitSimpleWindow( Wt::WSuggestionPopup *materialSuggestion,
                  InterSpec* viewer );
  
  virtual ~DetectionLimitSimpleWindow();

  DetectionLimitSimple *tool();
  
protected:
  DetectionLimitSimple *m_tool;
};//class DetectionLimitSimple


/**
 */
class DetectionLimitSimple : public Wt::WContainerWidget
{
public:
  
  DetectionLimitSimple( Wt::WSuggestionPopup *materialSuggestion,
                  InterSpec *specViewer,
                  Wt::WContainerWidget *parent = 0 );
  
  
  virtual ~DetectionLimitSimple();
  
  /**
   
   @param age Age of nuclide - if negative, use default age for the nuclide
   */
  void setNuclide( const SandiaDecay::Nuclide *nuc, const double age, const double energy );
  
  float photopeakEnergy() const;
  const SandiaDecay::Nuclide *nuclide() const;
  
  /** Handles receiving a "deep-link" url starting with "interspec://simple-mda/...".
   
   Example URIs:
   - "interspec://simple-mda/convoluted?ver=1&nuc=u235&energy=185&dist=100cm&..."
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://simple-mda/currie?nuc=u238...", then this string would be "curie?nuc=u238...",
          showing the Curie-style limit.
          This string is is in standard URL format of "key1=value1&key2=value2&..." with ordering not mattering.
          Capitalization is not important.
          Assumes the string passed in has already been url-decoded.
          If not a valid path or query_str, throws exception.
   */
  void handleAppUrl( std::string uri );
  
  /** Encodes current tool state to app-url format.  Returned string does not include the
   "interspec://" protocol, or "simple-mda" authority; so will look something like "decon?nuc=Cs137&energy=661&dist=100 cm&...",
   The path part of the URI specifies tab the tool is on.
   and it will not be url-encoded (will likely contain spaces).
   */
  std::string encodeStateToUrl() const;
protected:
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  void init();
  
  void roiDraggedCallback( double new_roi_lower_energy,
                   double new_roi_upper_energy,
                   double new_roi_px,
                   double original_roi_lower_energy,
                   const std::string &spec_type,
                   bool is_final_range );
  
  void handleMethodChanged();
  
  void handleNuclideChanged();
  
  void handleGammaChanged();
  
  void handleDistanceChanged();
  
  void handleConfidenceLevelChanged();

  void handleAdvancedToggled();
  void handleAlphaChanged();
  void handleBetaChanged();
  void handleSystematicUncertChanged();

  /** The combined relative systematic uncertainty to hand to
   `DetectionLimitCalc::CurrieMdaInput::additional_uncertainty`.

   Zero when "Advanced" is unchecked or both fields are empty, so that with the section off the tool
   computes exactly what it computed before the section existed.

   @param [out] note Set to a user-facing note when the combination had to be qualified - a distance
          uncertainty dropped because the detector response is fixed-geometry.  Left untouched
          otherwise.  Cannot be `m_warningTxt` directly: this runs from `updateResult()`, and
          `updateSpectrumDecorationsAndResultText()` clears that text afterwards.
   @returns The relative (not percent) 1-sigma systematic uncertainty, in [0,1).

   Throws `std::runtime_error` carrying a localized message when the entered strings cannot be used -
   an unparseable field, a distance uncertainty with no distance to divide by, or a combination at or
   above 100% (which `currie_mda_calc` rejects).
   */
  float currentSystematicUncertainty( Wt::WString &note ) const;

  void handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf );
  
  void handleFitFwhmRequested();
  
  void handleSelectDetectorRequested();
  
  void handleSpectrumChanged( const SpecUtils::SpectrumType type );
  
  void handleUserChangedRoi();
  
  void handleUserChangedNumSideChannel();
  
  void setFwhmFromEstimate();
  void handleUserChangedFwhm();
  
  void handleDeconPriorChange();
  void handleNoSignalPresentChanged();
  void handleDeconContinuumTypeChange();

  /** Toggles enabled state of the planned-measurement-time input, repopulates the field with
   the current foreground real time when the checkbox is off OR when the field is
   empty (so the displayed value is always meaningful), and schedules a recompute. */
  void handlePlanTimeChanged();

  /** The planned measurement time, in seconds, or zero when none was requested.

   Zero when the control is hidden or unchecked, or the field is blank - all of which mean "not
   asked for".  Throws only when the control is active and the field holds a non-empty but
   unparseable or non-positive value. */
  double currentPlanTimeSeconds() const;

  /** Which measurement the deconvolution limit describes.

   The "background spectrum" checkbox means the same assertion for both methods, but drives
   different machinery: zero Currie side channels, or this measurement model.  Always
   `CurrentSpectrum` under the Currie method, which has no notion of a future measurement. */
  DetectionLimitCalc::DeconMeasurementModel currentMeasurementModel() const;

  /** The spectra and exposure bookkeeping the current method / model combination needs.

   Thin wrapper over `DetectionLimitCalc::plan_measurement`, which is where the subtlety lives:
   under a background reference the deconvolution must see the UNSCALED spectrum with the planned
   time carried as an exposure, or the projection is applied twice. */
  DetectionLimitCalc::PlannedMeasurement currentEffectiveForeground() const;
  
  void updateSpectrumDecorationsAndResultText();
  
  /** Returns the confidence level (ex., 0.95, 0.9973, etc) that is selected by the GUI. */
  double currentConfidenceLevel() const;
  
  SimpleDialog *createDeconvolutionLimitMoreInfo();
  
  void createMoreInfoWindow();
  void handleMoreInfoWindowClose( SimpleDialog *dialog );
  void programmaticallyCloseMoreInfoWindow();
  
  void updateResult();
  
  /** Currently only update or not */
  enum RenderActions
  {
    UpdateLimit = 0x01,
    AddUndoRedoStep = 0x02,
    UpdateDisplayedSpectrum = 0x04,
    UpdateSpectrumDecorations = 0x08
  };//enum RenderActions
  
  Wt::WFlags<RenderActions> m_renderFlags;
  
  
  InterSpec *m_viewer;
  Wt::WSuggestionPopup *m_materialSuggest;
  // ShieldingSelect *m_enterShieldingSelect;
  
  Wt::WStackedWidget *m_chartErrMsgStack;
  Wt::WText *m_errMsg;
  Wt::WPushButton *m_fitFwhmBtn;
  
  D3SpectrumDisplayDiv *m_spectrum;
  PeakModel *m_peakModel;

  Wt::WText *m_resultTxt;
  
  Wt::WPushButton *m_moreInfoButton;
  
  Wt::WLineEdit *m_nuclideEdit;
  Wt::WLineEdit *m_nuclideAgeEdit;
  NuclideSourceEnterController *m_nucEnterController;
  
  Wt::WComboBox *m_photoPeakEnergy;
  
  /** We'll keep a copy of the energies in `m_photoPeakEnergy`, for ease of access.
   
   We could instead use our own `WAbstractItemModel` to power `m_photoPeakEnergy`,
   but just doing the dumb thing for the moment.
   */
  std::vector<std::pair<double,double>> m_photoPeakEnergiesAndBr;
  
  
  Wt::WLineEdit *m_distance;
  
  enum ConfidenceLevel 
  {
    NinetyFivePercent, //0.95
    NinetyNinePercent, //0.99
    OneSigma,          //0.682689492137086
    TwoSigma,          //0.954499736103642
    ThreeSigma,        //0.997300203936740
    FourSigma,         //0.999936657516334
    FiveSigma,         //0.999999426696856
    NumConfidenceLevel
  };//enum ConfidenceLevel
  
  Wt::WComboBox *m_confidenceLevel;
  
  DetectorDisplay *m_detectorDisplay;
  
  /** Enum to track ID of `m_methodGroup` */
  enum class MethodIds : int
  {
    Currie = 0,
    Deconvolution = 1
  };
  
  Wt::WButtonGroup *m_methodGroup;
  Wt::WText *m_methodDescription;

  /** Reveals the advanced statistical inputs below the rest of the tool.  Unchecked by default;
   while unchecked none of those values reach the calculation, so the answers are exactly what they
   were before this section existed. */
  Wt::WCheckBox *m_advancedCb;

  /** Container for the advanced inputs.  A sibling of the main "GeneralInput" grid rather than an
   eleventh row of it: `GridLayoutHelpers.css` stops at `GridTenthRow`, the phone-only row-shift table
   in `DetectionLimitSimple.css` renumbers `GridThirdRow`..`GridTenthRow` inside `.GeneralInput`, and
   the phone input-width rules are `.GeneralInput`-scoped as well - so an eleventh row would mean
   growing all three.  A sibling also gets its own, saner phone layout. */
  Wt::WContainerWidget *m_advancedDiv;

  /** Probability of deciding a signal is present when there is none.  Sets k_alpha, and so the
   decision threshold L_c.  A probability, not a percent, in (0, 0.5); 0.05 is usual.
   \sa DetectionLimitCalc::CurrieMdaInput::alpha */
  NativeFloatSpinBox *m_alpha;

  /** Probability of failing to detect a signal whose true size is the detection limit.  Sets k_beta,
   and so L_d.  A probability in (0, 0.5); 0.05 is usual.
   \sa DetectionLimitCalc::CurrieMdaInput::beta */
  NativeFloatSpinBox *m_beta;

  /** False until the user edits the corresponding field.  While false the field tracks the
   confidence-level combo as 1 - CL, and - importantly - the *sentinel* rather than the field's
   contents is handed to the calculation, so the answers stay bit-for-bit what they were before
   alpha/beta were separable.  (The field shows a rounding of 1 - CL; feeding that back would move
   the last digits for the sigma-style confidence levels.)  The flags are encoded in the state URI by
   the presence or absence of the ALPHA / BETA tokens. */
  bool m_alphaUserSet;
  bool m_betaUserSet;

  /** 1-sigma uncertainty of the source-to-detector distance, as a length string in the same grammar
   as the distance field above it ("1 cm").  Empty, or a zero length, means none. */
  Wt::WLineEdit *m_distanceUncert;

  /** Most recent valid distance-uncertainty text, so an invalid entry can be reverted - same pattern,
   and same reason, as `m_prevDistance`. */
  Wt::WString m_prevDistanceUncert;

  /** 1-sigma *relative* uncertainty of the counts expected per unit activity, as a PERCENT ("5").
   Empty or zero means none.

   Covers everything that scales the expected counts and is not counting statistics: the detector
   efficiency curve and the gamma branching ratio, which enter identically and so cannot be usefully
   separated here.

   TODO: `DetectorPeakResponse` carries no uncertainty on its efficiency today.  When it does,
         pre-fill this field from the DRF (combined in quadrature with the branching ratio's own
         uncertainty, if `SandiaDecay` ever carries one), leaving the user able to override. */
  Wt::WLineEdit *m_effUncert;

  /** Note under the advanced inputs saying they apply to the Currie method; shown only while the
   Deconvolution method is selected. */
  Wt::WText *m_advancedNote;

  /** A note about the systematic-uncertainty combination that has to survive from `updateResult()`
   to `updateSpectrumDecorationsAndResultText()`, which runs afterwards and clears `m_warningTxt`. */
  Wt::WString m_systematicNote;


  float m_numFwhmWide;
  NativeFloatSpinBox *m_lowerRoi;
  NativeFloatSpinBox *m_upperRoi;
  Wt::WLabel *m_numSideChannelLabel;
  Wt::WSpinBox *m_numSideChannel;
  
  NativeFloatSpinBox *m_fwhm;
  Wt::WText *m_fwhmSuggestTxt;
  Wt::WPushButton *m_addFwhmBtn;
  Wt::WPushButton *m_selectDetectorBtn;
  
  WContainerWidget *m_isBackgroundDiv;
  Wt::WCheckBox *m_isBackgroundSpectrum;
  Wt::WLabel *m_continuumPriorLabel;
  Wt::WComboBox *m_continuumPrior;
  Wt::WLabel *m_continuumTypeLabel;
  Wt::WComboBox *m_continuumType;

  /** The planned measurement time (T_s) - the dwell being asked about.

   Under the Currie method, and under the deconvolution method describing the current spectrum, the
   spectrum's counts are projected to this real time.  Under a background reference the spectrum is
   left alone and this becomes the exposure of the predicted future measurement.  Hidden for
   deconvolution + current spectrum, where projecting the spectrum in hand and then bounding signal
   in it is circular.  When unchecked, `m_planTimeEdit` shows the current foreground's real time for
   reference and is disabled.  \sa DetectionLimitCalc::plan_measurement */
  Wt::WCheckBox        *m_planTimeCb;
  Wt::WLineEdit        *m_planTimeEdit;
  Wt::WContainerWidget *m_planTimeDiv;

  /** The help icon beside the background checkbox; its tooltip is swapped by method, because the
   checkbox means the same assertion but drives different machinery for each. */
  Wt::WImage           *m_isBackgroundHelpImg;

  /** Notes that qualify a limit that WAS produced.  A sibling of `m_resultTxt` rather than a third
   page of `m_chartErrMsgStack`, which is either-or: an error page OR the chart, with nowhere to
   put a warning beside a successful result.  These were dropped entirely before Increment C. */
  Wt::WText            *m_warningTxt;
  
  SimpleDialog *m_moreInfoWindow;
  
  const SandiaDecay::Nuclide *m_currentNuclide;
  double m_currentAge;
  double m_currentEnergy;
  
  /** If false, then only gamma specified will be used in the limit.
   Otherwise, all gammas within the ROI will be added to gamma specified.
   */
  bool m_allGammasInRoi;
  
  /** The most  recent valid distance - used to reset the distance field, if user enters an invalid distance. */
  Wt::WString m_prevDistance;
  
  // I'm still not sure how to handle undo-redo.
  //  Maybe at first we'll use the URI string, but maybe it would be easier/better
  //  to just capture input values, and use those.
  
  /** For tracking undo/redo, we will keep track of the widgets state, as a URI.
   \sa encodeStateToUrl
   \sa handleAppUrl
   */
  std::string m_stateUri;
  
  std::shared_ptr<const DetectionLimitCalc::CurrieMdaInput> m_currentCurrieInput;
  std::shared_ptr<const DetectionLimitCalc::CurrieMdaResult> m_currentCurrieResults;
  
  std::shared_ptr<const DetectionLimitCalc::DeconComputeInput> m_currentDeconInput;
  std::shared_ptr<const DetectionLimitCalc::DeconActivityOrDistanceLimitResult> m_currentDeconResults;

  /** The predicted spread of the limit, when a measurement time other than the spectrum's is being
   asked about; invalid otherwise.

   In counts, matching the limit the active method reports.  A projection is a prediction about a
   measurement nobody has taken, and quoting only its middle hides how far the answer can move - by
   about `sqrt(1+k)` more than a plain scaling implies, `k` being the projection factor.
   \sa DetectionLimitCalc::currie_projected_limit, DetectionLimitCalc::decon_projected_limit
   */
  DetectionLimitCalc::ProjectedLimit m_currentProjectedLimit;
};//class DoseCalcWidget

#endif //DetectionLimitSimple_h
