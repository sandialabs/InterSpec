#ifndef EnergyCalGainMatch_h
#define EnergyCalGainMatch_h
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

#include <set>
#include <map>
#include <deque>
#include <tuple>
#include <memory>
#include <string>
#include <vector>

#include <Wt/WDialog.h>
#include <Wt/WContainerWidget.h>

// Forward Declarations
class SpecMeas;
class InterSpec;
class AuxWindow;
class EnergyCalTool;
struct PeakFitDetPrefs;
class D3SpectrumDisplayDiv;
class NativeFloatSpinBox;

namespace Wt
{
  class WText;
  class WComboBox;
  class WCheckBox;
  class WPushButton;
  class WButtonGroup;
  class WRadioButton;
}

namespace SpecUtils
{
  class Measurement;
  struct EnergyCalibration;
  enum class SpectrumType : int;
}


/** GUI-free gain-matching computation, so unit tests
 (target/testing/test_EnergyCalGainMatch.cpp) can exercise each stage without a WApplication.

 The overall approach to make a "test" spectrum consistent with a "reference" spectrum is to
 find the linear energy transform, E' = gain*E + offset, of the test spectrums energy
 calibration that makes the two spectra as similar as possible:
   1) A coarse gain estimate from cross-correlating the two spectra on a log-energy grid
      (a pure gain difference is a constant shift in log-energy), which is robust to even
      wildly wrong starting gains.
   2) A refinement minimizing the chi2 between the spectra (with the overall normalization
      of the test spectrum floating freely), using 1-D Brent minimizations.
   3) Optionally, a photopeak is fit in both spectra and the test spectrums gain is scaled
      so its fitted mean coincides with the references fitted mean.
 */
namespace GainMatchCalc
{
  enum class WarningType : int
  {
    LowStatistics,      //Too few counts in the energy range; spectrum excluded from matching
    WeakCorrelation,    //Best cross-correlation was weak; result may be unreliable
    CorrelationAtEdge,  //Best correlation at edge of gain search window; true gain may be beyond it
    PeakFitFailed,      //Peak refinement requested, but peak couldnt be fit; chi2 result kept
    AmbiguousCorrelation, //A competing correlation maximum was nearly as good; result may be wrong
    MatchFailed,        //Matching this spectrum failed entirely; its calibration left unchanged
    NumWarningType
  };//enum class WarningType


  struct MatchOptions
  {
    /** Lower energy (keV) of range used to compare spectra. */
    float lower_energy = 200.0f;

    /** Upper energy (keV) of comparison range; a value <= lower_energy means no upper limit. */
    float upper_energy = 0.0f;

    /** If true, fit the offset term in addition to gain. */
    bool fit_offset = false;

    /** If greater than zero, after the chi2 match, fit a photopeak at this energy (keV) in the
     reference and each test spectrum, and adjust each test gain so its fitted mean matches the
     references fitted mean.
     */
    double ref_peak_energy = 0.0;
  };//struct MatchOptions


  struct SpectrumInput
  {
    /** Label for messages (detector name, or file description). */
    std::string name;

    /** The spectrum for this input; e.g., a single detectors data summed over the samples of
     interest, or a whole files displayed sum.  Must have a valid energy calibration.
     */
    std::shared_ptr<const SpecUtils::Measurement> spectrum;

    /** Peak-fit preferences used for the (optional) peak-refinement stage; may be nullptr, in
     which case defaults appropriate to the spectrums resolution are used.
     */
    std::shared_ptr<const PeakFitDetPrefs> fit_prefs;
  };//struct SpectrumInput


  struct SpectrumResult
  {
    /** False if this spectrum was excluded (invalid, low statistics), or was not matched. */
    bool used = false;

    /** Multiplicative term of the energy transform E' = gain*E + offset. */
    double gain = 1.0;

    /** Additive term, in keV, of the energy transform. */
    double offset = 0.0;

    /** Peak normalized cross-correlation value from stage 1, in [-1,1]. */
    double correlation = 0.0;

    /** Reduced chi2 of (reference vs scaled test) after stage 2; qualitative only since
     neighboring channels are correlated by the rebinning.
     */
    double chi2_dof = 0.0;

    /** Fitted peak mean (keV, in the stage-2 adjusted energy scale) if peak refinement ran for
     this spectrum; <= 0.0 otherwise.
     */
    double fit_peak_mean = 0.0;

    /** The updated calibration for this spectrum; nullptr for the reference spectrum. */
    std::shared_ptr<const SpecUtils::EnergyCalibration> updated_cal;

    std::vector<WarningType> warnings;
  };//struct SpectrumResult


  struct MatchResults
  {
    /** Index into the input vector of the reference spectrum used. */
    size_t reference_index = 0;

    /** The references fitted peak mean (keV) when peak refinement ran; <= 0.0 otherwise. */
    double ref_peak_mean = 0.0;

    /** Result for each input spectrum; parallel to the input vector. */
    std::vector<SpectrumResult> results;
  };//struct MatchResults


  /** Returns a new calibration transformed so E' = gain*E + offset.

   Polynomial (incl. UnspecifiedUsingDefaultPolynomial) and FullRangeFraction: energy is linear
   in the coefficients (including the FRF low-energy 1/(1+60x) term), so all coefficients are
   multiplied by gain, then offset added to the zeroth - this is exact.
   LowerChannelEdge: each channel edge is transformed.
   Deviation pairs are copied unchanged: they are small nonlinear corrections applied after the
   polynomial, so the error of not scaling them is second order for the few-percent gain changes
   this tool produces (documented limitation for very large gain changes).

   Throws std::exception if cal is null/invalid, gain is not positive, or the result is invalid.
   */
  std::shared_ptr<SpecUtils::EnergyCalibration>
  transform_calibration( const std::shared_ptr<const SpecUtils::EnergyCalibration> &cal,
                         const double gain, const double offset );

  /** Sum of gamma counts in channels whose lower edge is in [lower,upper);
   upper <= lower means no upper limit.  Returns 0.0 for null/invalid spectra.
   */
  double counts_in_range( const std::shared_ptr<const SpecUtils::Measurement> &spec,
                          const float lower_energy, const float upper_energy );

  /** Stage 1: coarse pure-gain estimate, by maximizing the normalized (Pearson) correlation of
   sqrt(counts/keV) between the two spectra resampled onto uniform-in-log(energy) grids, as a
   function of the log-energy shift (i.e., of gain).  The count normalizations of the spectra do
   not matter.  Searches gains in [1/8,8] in ~0.2% steps, with parabolic interpolation of the
   correlation peak.

   @param correlation [out] the best correlation value found
   @param warnings [out] appended with WeakCorrelation and/or CorrelationAtEdge as applicable
   @returns the gain estimate

   Throws std::exception if either spectrum is null/invalid or the usable range is too narrow.
   */
  double coarse_gain_xcorr( const std::shared_ptr<const SpecUtils::Measurement> &test,
                            const std::shared_ptr<const SpecUtils::Measurement> &reference,
                            const float lower_energy, const float upper_energy,
                            double &correlation,
                            std::vector<WarningType> &warnings );

  /** Stage 2 objective: chi2 of the reference spectrum vs the test spectrum transformed by
   (gain,offset) and rebinned onto the references channel edges (restricted to the energy range),
   with the overall normalization of the test spectrum solved for analytically.
   Reference channels in the energy range that the transformed test spectrum does not cover
   count as zeros (i.e., are penalized, not dropped), so shrinking the overlap cant lower chi2.
   Exposed primarily for unit tests.
   */
  double chi2_for_gain_offset( const std::shared_ptr<const SpecUtils::Measurement> &test,
                               const std::shared_ptr<const SpecUtils::Measurement> &reference,
                               const float lower_energy, const float upper_energy,
                               const double gain, const double offset,
                               size_t &num_channels_used );

  /** Stage 2: refines gain (and offset, when fit_offset) starting from the passed-in values.

   The comparison window is fixed from the starting transform (see #chi2_for_gain_offset), then
   gain (and offset, alternating) are minimized by a coarse-to-fine direct scan followed by a
   Brent polish - a plain bracketed minimization can entirely miss the very narrow (~0.1% in
   gain) chi2 valleys that sharp HPGe photopeaks produce.

   Updates gain/offset in place; returns the reduced chi2.
   */
  double refine_gain_offset( const std::shared_ptr<const SpecUtils::Measurement> &test,
                             const std::shared_ptr<const SpecUtils::Measurement> &reference,
                             const float lower_energy, const float upper_energy,
                             const bool fit_offset,
                             double &gain, double &offset );

  /** Stage 3 helper: fits a photopeak near expected_energy (headless), returning the fitted mean.

   @param fit_prefs may be nullptr, in which case defaults for the spectrums resolution are used.

   Throws std::exception if no acceptable peak is found (fitted mean must be within
   max(2*expected_fwhm, 3% of expected_energy) of expected_energy).
   */
  double fit_peak_mean( const std::shared_ptr<const SpecUtils::Measurement> &spec,
                        const double expected_energy,
                        const std::shared_ptr<const PeakFitDetPrefs> &fit_prefs );

  /** Top-level driver: stages 1-3 for each usable non-reference input.

   @param inputs the spectra to match
   @param reference_index index of the reference spectrum, or negative to auto-select the
          spectrum with the most counts in the energy range
   @param options see #MatchOptions

   Throws std::runtime_error if fewer than two usable spectra.
   */
  MatchResults match( const std::vector<SpectrumInput> &inputs,
                      const int reference_index,
                      const MatchOptions &options );


  /** A ready-to-apply gain-match result, produced by #analyzeForAutoMatch when a
   multi-detector spectrum is loaded, and consumed by #applyDetectorGains.
   */
  struct DetectorMatchSuggestion
  {
    std::weak_ptr<SpecMeas> file;
    SpecUtils::SpectrumType type;
    std::set<int> displayed_samples;
    /** The largest implied energy discrepancy between detectors (keV) - the measure used
     to decide the match is worth suggesting.
     */
    double max_shift_kev = 0.0;
    /** Per detector: (detector name, gain, offset). */
    std::vector<std::tuple<std::string,double,double>> per_detector;
  };//struct DetectorMatchSuggestion


  /** Output of #matchDetectorInputs: the per-detector corrections and whether applying them
   would meaningfully improve consistency.
   */
  struct DetectorMatchResult
  {
    bool beneficial = false;
    double max_shift_kev = 0.0;
    std::vector<std::tuple<std::string,double,double>> per_detector;  //(name, gain, offset)
  };//struct DetectorMatchResult


  /** Builds one summed input spectrum per gamma detector with data in \p displayed_samples,
   applying the cost gate (returns empty if <2 usable detectors, or the file is too large to
   analyze quickly).  Reads the SpecMeas, so must be called on the session thread; the summed
   Measurements it returns are fresh and immutable, so the result can then be handed to a
   worker thread.  Never throws.
   */
  std::vector<SpectrumInput>
  buildDetectorInputs( const std::shared_ptr<SpecMeas> &meas,
                       const std::set<int> &displayed_samples );

  /** Pure gain-match analysis on pre-built inputs (from #buildDetectorInputs): matches the
   detectors and decides, via a peak-width-relative measure, whether matching would visibly
   sharpen the summed data (largest implied inter-detector shift exceeds both \p min_shift_kev
   and ~1/3 of the reference peak FWHM).  GUI- and SpecMeas-free, so safe to run on a worker
   thread.  Never throws.
   */
  DetectorMatchResult
  matchDetectorInputs( const std::vector<SpectrumInput> &inputs,
                       const double min_shift_kev = 3.0 );

  /** Synchronous convenience wrapper: #buildDetectorInputs + #matchDetectorInputs, returning a
   #DetectorMatchSuggestion (or nullptr if not beneficial).  Used by tests and simple callers;
   the on-load path in InterSpec splits the two halves across threads instead.  Never throws.
   */
  std::shared_ptr<DetectorMatchSuggestion>
  analyzeForAutoMatch( const SpecUtils::SpectrumType type,
                       const std::shared_ptr<SpecMeas> &meas,
                       const std::set<int> &displayed_samples,
                       const double min_shift_kev = 3.0 );


  /** Applies a per-detector gain/offset map to a files measurements: transforms each
   in-scope measurements calibration, registers one undo/redo step, refreshes the display,
   and relaunches the automated-search ("hint") peaks for the displayed spectrum.

   Shared by the Gain Match dialogs "within file" apply path and the auto-match "Accept"
   toast button.

   @param apply_samples the sample numbers to modify (all samples in the file, or just the
          displayed ones)
   @param per_detector (detector name, gain, offset) for each detector to change
   @returns the number of detectors whose calibration was changed
   */
  size_t applyDetectorGains( InterSpec *interspec,
                             const SpecUtils::SpectrumType type,
                             const std::shared_ptr<SpecMeas> &meas,
                             const std::set<int> &apply_samples,
                             const std::vector<std::tuple<std::string,double,double>> &per_detector );


  /** A peak found in most/all detectors, used for the multi-peak refinement table.  The
   detectors are first aligned by the stage-1/2 (gain,offset) match, then the same peak is
   located in each; #detector_channels gives the channel it sits at in each input detector
   (calibration-invariant), so a polynomial can be fit per detector.
   */
  struct SharedPeak
  {
    /** The peaks energy (keV) in the reference detector - its location in the aligned data. */
    double energy = 0.0;

    /** True if this peak was matched to a nuclide-assigned peak, so #target_energy is the
     known gamma energy rather than the reference detectors fitted position.
     */
    bool identified = false;

    /** The energy (keV) each detectors channel for this peak should map to: the known gamma
     energy when #identified, otherwise the reference detectors position (#energy).
     */
    double target_energy = 0.0;

    /** A short label for the table, e.g. "356.0 keV" or "Ba133 356.0 keV". */
    std::string label;

    /** Number of input detectors this peak was located in. */
    int num_detectors = 0;

    /** Per input detector, the (fractional) channel this peak sits at, or < 0 if not found in
     that detector.  Parallel to the inputs vector passed to #findSharedPeaks.
     */
    std::vector<double> detector_channels;
  };//struct SharedPeak


  /** Finds peaks common to the detectors, for the multi-peak refinement.  Aligns each detector
   by its stage-1/2 (gain,offset) result, runs a fast candidate search per detector, clusters
   the candidates across detectors by energy, and keeps clusters present in a healthy fraction
   of detectors.  Clusters whose energy matches one of \p id_peaks (a nuclide-assigned peak) are
   marked identified, with #target_energy set to the assigned gamma energy.

   Pure and GUI-free (safe on a worker thread), but can take up to a few seconds for many
   detectors.  Never throws.
   */
  std::vector<SharedPeak>
  findSharedPeaks( const std::vector<SpectrumInput> &inputs,
                   const MatchResults &stage2,
                   const std::deque<std::shared_ptr<const PeakDef>> &id_peaks,
                   const float lower_energy, const float upper_energy );

  /** Multi-peak refinement: for each detector, fits a polynomial energy calibration of the given
   order to the selected shared peaks (channel -> target energy), returning the new calibration
   per input detector (nullptr where unchanged or too few peaks).

   @param order 1 => offset+gain, 2 => +quadratic, 3 => +cubic.  A detector needs at least
          order+1 selected peaks located in it to be fit.
   @param use_peak parallel to \p peaks; only peaks with a true entry are used.

   The reference detector is only re-fit when a selected peak is identified (so it can move
   toward the known energy); otherwise it stays fixed (nullptr).  Never throws.
   */
  std::vector<std::shared_ptr<const SpecUtils::EnergyCalibration>>
  refineWithSharedPeaks( const std::vector<SpectrumInput> &inputs,
                         const MatchResults &stage2,
                         const std::vector<SharedPeak> &peaks,
                         const std::vector<bool> &use_peak,
                         const int order );

  /** Applies a per-detector calibration directly (used by the multi-peak refinement, where the
   correction is a full polynomial rather than a linear gain/offset): sets each in-scope
   measurements calibration, registers one undo/redo step, refreshes, and relaunches the
   automated-search peaks.  Measurements whose channel count differs from the fitted calibration
   are skipped.  Returns the number of detectors changed.
   */
  size_t applyDetectorCals( InterSpec *interspec,
                            const SpecUtils::SpectrumType type,
                            const std::shared_ptr<SpecMeas> &meas,
                            const std::set<int> &apply_samples,
                            const std::vector<std::pair<std::string,
                                   std::shared_ptr<const SpecUtils::EnergyCalibration>>> &per_detector );
}//namespace GainMatchCalc


/** Content widget of the "Gain Match..." dialog opened from the Energy Calibration tabs
 "More Actions" column; created by EnergyCalAddActionsWindow, with buttons in the AuxWindow
 footer (analogous to EnergyCalMultiFile).

 Lets the user match the energy calibrations of the individual detectors within a spectrum file
 (or of the foreground/background/secondary files to each other), previewing each spectrum as a
 separate line on a chart.  No SpecMeas is modified until the user accepts the dialog.
 */
class EnergyCalGainMatch : public Wt::WContainerWidget
{
public:
  EnergyCalGainMatch( EnergyCalTool *cal, AuxWindow *parent );
  virtual ~EnergyCalGainMatch();

  void handleFinish( Wt::DialogCode result );

protected:
  enum class MatchMode : int
  {
    WithinFile,   //match individual detectors of one file
    BetweenFiles  //match displayed foreground/background/secondary (summed) to each other
  };//enum class MatchMode

  enum GainMatchRenderActions
  {
    UpdateRows    = 0x01,  //rebuild the row list (mode or file changed)
    UpdateCalc    = 0x02,  //re-run GainMatchCalc::match (+ find shared peaks if refining)
    UpdateRefine  = 0x04,  //re-apply the multi-peak refinement using the cached shared peaks
    UpdatePreview = 0x08   //re-send preview spectra to the chart
  };//enum GainMatchRenderActions

  /** One selectable spectrum: a detector of the chosen file, or a displayed file. */
  struct Row
  {
    Wt::WCheckBox *use = nullptr;
    Wt::WRadioButton *reference = nullptr;
    Wt::WText *name = nullptr;
    Wt::WText *gain = nullptr;
    Wt::WText *offset = nullptr;
    Wt::WText *status = nullptr;
    std::string det_name;                //WithinFile mode
    SpecUtils::SpectrumType spec_type;   //file this row belongs to
    std::shared_ptr<const SpecUtils::Measurement> spectrum;  //summed input (a copy)
    std::shared_ptr<const SpecUtils::EnergyCalibration> orig_cal;
    GainMatchCalc::SpectrumResult result;
    bool refined = false;                //true if result.updated_cal came from multi-peak refinement
  };//struct Row

  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;

  MatchMode currentMode() const;
  SpecUtils::SpectrumType selectedFileType() const;
  void handleModeChange();
  void handleOptionsChange();
  void rebuildRows();
  bool refineWithPeaksEnabled() const;  //true when the "refine with peaks" box is on and applicable
  int selectedFitOrder() const;         //1 = offset+gain, 2 = +quadratic, 3 = +cubic
  void rebuildPeaksTable();             //fills m_peaksTable from m_sharedPeaks
  void applyRefineToRows();             //refineWithSharedPeaks -> each rows updated_cal
  void launchSharedPeakSearch();        //runs findSharedPeaks on a worker thread
  //Session-thread continuation of launchSharedPeakSearch(); public so the worker post-back
  //  (which re-looks-up the widget by id) can call it.  `generation` guards against stale results.
public:
  void onSharedPeaksFound( int generation, const std::vector<GainMatchCalc::SharedPeak> &peaks );
protected:
  float upperEnergyLimit();             // <= lower means no limit
  void doCalcUpdate();
  void doRefineUpdate();
  void updatePreview();
  void updateResultColumns();
  void applyChanges();

  InterSpec *m_interspec;
  EnergyCalTool *m_calibrator;
  AuxWindow *m_parent;

  std::shared_ptr<Wt::WButtonGroup> m_modeGroup;
  Wt::WRadioButton *m_withinFileBtn;
  Wt::WRadioButton *m_betweenFilesBtn;
  Wt::WComboBox *m_fileSelect;
  std::vector<SpecUtils::SpectrumType> m_fileSelectTypes;  //parallel to m_fileSelect items
  NativeFloatSpinBox *m_lowerEnergy;
  NativeFloatSpinBox *m_upperEnergy;   //empty text means no upper limit
  Wt::WCheckBox *m_fitOffset;
  Wt::WCheckBox *m_refineWithPeaks;    //enable the multi-peak polynomial refinement (WithinFile)
  Wt::WComboBox *m_fitOrderCombo;      //polynomial order for the refinement (offset+gain / +quad / +cubic)
  Wt::WContainerWidget *m_peaksRow;    //container holding the refine controls + peaks table
  Wt::WContainerWidget *m_peaksTable;  //one row per shared peak (Use checkbox + label + found-in)
  Wt::WText *m_peaksStatus;            //e.g. "8 shared peaks" or "not enough shared peaks"
  Wt::WCheckBox *m_allSamples;
  Wt::WCheckBox *m_showOriginal;
  Wt::WContainerWidget *m_rowTable;
  std::shared_ptr<Wt::WButtonGroup> m_refGroup;
  std::vector<Row> m_rows;

  // Cached inputs/results of the last stage-1/2 match, so the multi-peak refinement can be
  //  re-applied (on peak selection / order changes) without re-running the match or peak search.
  std::vector<GainMatchCalc::SpectrumInput> m_lastInputs;
  std::vector<size_t> m_lastInputRows;   //row index for each entry in m_lastInputs
  GainMatchCalc::MatchResults m_lastStage2;
  std::vector<GainMatchCalc::SharedPeak> m_sharedPeaks;
  std::vector<Wt::WCheckBox *> m_peakUseCbs;  //parallel to m_sharedPeaks
  bool m_sharedPeaksValid;               //false when the match changed and peaks must be re-found
  bool m_findingSharedPeaks;             //true while the worker-thread search is in flight
  int m_refineGeneration;                //bumped per search; the post-back discards stale results
  D3SpectrumDisplayDiv *m_chart;
  size_t m_rangeHighlightId;  //decorative highlight region id; 0 means none
  bool m_chartXRangeSet;      //whether initial x-domain has been sent to client
  Wt::WText *m_status;
  Wt::WPushButton *m_cancel;
  Wt::WPushButton *m_use;

  Wt::WFlags<GainMatchRenderActions> m_renderFlags;

  /** Files snapshotted at row build; verified still displayed in applyChanges(). */
  std::shared_ptr<SpecMeas> m_withinFileMeas;
  std::map<SpecUtils::SpectrumType,std::shared_ptr<SpecMeas>> m_betweenFileMeas;
};//class EnergyCalGainMatch

#endif //EnergyCalGainMatch_h
