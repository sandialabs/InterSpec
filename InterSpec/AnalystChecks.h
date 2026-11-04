#ifndef AnalystChecks_h
#define AnalystChecks_h
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

#include <deque>
#include <string>
#include <vector>
#include <memory>
#include <variant>
#include <optional>
#include <functional>

#include "SpecUtils/SpecFile.h"

#include <Wt/WFlags.h>

#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/ReactionGamma.h" //for ReactionGamma::Reaction

// Forward declarations
class PeakDef;
class InterSpec;

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}

namespace AnalystChecks
{
  /** Options for peak detection analysis. */
  struct DetectedPeaksOptions {
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
    bool nonBackgroundPeaksOnly = false;
    std::optional<double> lowerEnergy;
    std::optional<double> upperEnergy;
  };
  
  /** Results of peak detection analysis. */
  struct DetectedPeakStatus {
    std::vector<std::shared_ptr<const PeakDef>> peaks;

    /** Subset of `peaks` that are user-defined analysis peaks (not auto-detected peaks). */
    std::vector<std::shared_ptr<const PeakDef>> analysis_peaks;
  };
  
  /** Perform automated peak detection on the specified spectrum.
   
   This function combines user-defined peaks with automatically detected peaks,
   avoiding duplicates by checking peak overlap.
   
   @param options Specifies which spectrum to analyze and optional user session
   @param interspec Pointer to the InterSpec session containing the spectrum
   @return DetectedPeakStatus containing the found peaks and session info
   @throws std::runtime_error if spectrum is not available or other errors occur
   */
  InterSpec_API DetectedPeakStatus detected_peaks(const DetectedPeaksOptions& options, InterSpec* interspec);

  /** The amplitude uncertainty the elevated-above-background tests use: the detection-equivalent one,
   the amplitude over its detection significance (see PeakFitLM::peak_detection_significance, chi2
   weights), which their thresholds were tuned on, rather than the larger, marginal, uncertainty a fit
   reports; the reported one where that cannot be had.  `roi_peaks` should hold the other peaks of
   `peak`'s ROI; `data` is the spectrum it was fit to.
   */
  InterSpec_API double elevation_test_uncert( const PeakDef &peak,
                                              const std::vector<std::shared_ptr<const PeakDef>> &roi_peaks,
                                              const std::shared_ptr<const SpecUtils::Measurement> &data );
  
  
  struct FitPeakOptions {
    bool doNotAddToAnalysisPeaks = false;
    double energy = 0.0;
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
    std::optional<std::string> source;

    /** A ROI the caller fully specifies: exactly these peaks are fit together in one ROI. */
    struct ExplicitRoi
    {
      /** Starting means (keV) of the peaks to fit; at most #sm_max_peaks. */
      std::vector<double> energies;

      /** Source for each peak (empty string for none); either empty, or the same size as `energies`. */
      std::vector<std::string> sources;

      /** ROI range (keV); a bound not given is chosen from the data around the outermost peaks. */
      std::optional<double> lowerEnergy;
      std::optional<double> upperEnergy;

      /** A `PeakContinuum::offset_type_str(...)` name; if not given, chosen from the number of peaks. */
      std::optional<std::string> continuumType;

      static constexpr size_t sm_max_peaks = 10;
    };//struct ExplicitRoi

    /** If set, `energy` and `source` are ignored, and the requested peaks are fit in one ROI,
     instead of searching for the best peak near `energy`.
     */
    std::optional<ExplicitRoi> explicitRoi;
  };

  /** What the data show near an energy where the (single-energy) automated-threshold peak search
   accepted no peak; see `FitPeakStatus::diagnostic`.
   */
  struct UnfitPeakDiagnostic
  {
    /** Why no peak was accepted, e.g., "significance 2.1 sigma is below the 5 sigma the automated
     search requires"; intended to be read by an LLM or a user. */
    std::string reason;

    /** The best candidate peak (and its ROI) found using the more lenient thresholds of a user
     double-click; null if not even that found a peak.  Never added to the session. */
    std::shared_ptr<const PeakDef> candidate;
    std::vector<std::shared_ptr<const PeakDef>> candidateRoiPeaks;

    /** Expected FWHM (keV) at the energy, and a description of where it came from. */
    double expectedFwhm = 0.0;
    std::string expectedFwhmSource;

    /** A Gaussian of the expected FWHM at exactly the requested energy, fit with a quadratic
     continuum over +-2 FWHM (see `FitPeaksForNuclides::detail::fixed_shape_peak_fit`). */
    bool fixedFitValid = false;
    double fixedFitArea = 0.0;        // signed net counts
    double fixedFitAreaUncert = 0.0;
    double fixedFitChi2Dof = -1.0;
    double fixedFitLowerEnergy = 0.0;
    double fixedFitUpperEnergy = 0.0;
  };//struct UnfitPeakDiagnostic


  struct FitPeakStatus {
    /** The fit peak; for a `FitPeakOptions::explicitRoi` fit, the peak nearest the first requested
     energy.  When `significant` is false, no peak was accepted, and this is the rejected candidate
     (`UnfitPeakDiagnostic::candidate`; may be null) - so test `significant`, not this, for success. */
    std::shared_ptr<const PeakDef> fitPeak;
    std::vector<std::shared_ptr<const PeakDef>> peaksInRoi;

    /** False when the single-energy search found no peak to add, and the caller asked to only
     probe the energy (`FitPeakOptions::doNotAddToAnalysisPeaks`); `diagnostic` then says why,
     and what the data show at the energy. */
    bool significant = true;
    std::optional<UnfitPeakDiagnostic> diagnostic;

    /** If this was a `FitPeakOptions::explicitRoi` fit; the members below are only used for these. */
    bool explicitRoi = false;

    /** Previously existing peaks replaced by the fit (a requested peak was at about their energy). */
    std::vector<std::shared_ptr<const PeakDef>> removedPeaks;

    /** Notes for the caller, e.g., weak peaks, or peaks that moved a lot. */
    std::vector<std::string> warnings;
  };

  /** Fits a peak, and unless `options.doNotAddToAnalysisPeaks`, adds it to the analysis peaks.

   Normally searches for the best peak near `options.energy`, with the fitter choosing the ROI and the
   peaks in it (like a user double-click, but with the automated-search significance requirements).

   If `options.explicitRoi` is set, instead creates one ROI with the given (or a default) range and
   continuum, holding one Gaussian peak per requested energy, all fit together, even if weak.
   Existing ROIs are combined into the new ROI (their peaks fit along with the requested ones) as a
   user double-click would combine them - or, if ROI bounds are given, if at least half of an existing
   ROI, or all its peaks, are inside the bounds.  An existing peak at about a requested energy is
   replaced by the requested peak.  Overlapping, but not combined, ROIs must not have peaks inside the
   new ROI; that, or overlapping a data-defined peak, is an error.  Records an undo/redo step.

   Throws std::exception, with a description intended to be actionable by an LLM, on error.
   */
  InterSpec_API FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec );

  /** Async version of fit_user_peak that runs the expensive peak search (or explicit ROI fit) in a
   background thread.
   The callback is called on the GUI thread with either a FitPeakStatus (success) or a string (error message).
   */
  using FitPeakCallback = std::function<void( std::variant<FitPeakStatus, std::string> )>;
  InterSpec_API void fit_user_peak_async( const FitPeakOptions &options, InterSpec *interspec,
                                          FitPeakCallback callback );

  struct ExplicitRoiSetup;
  struct ExplicitRoiFitResult;

  /** The stages of a `FitPeakOptions::explicitRoi` fit, which `fit_user_peak_async` runs on the GUI
   thread, a background thread, and the GUI thread, respectively.  Exposed so tests can interleave
   the stages of several requests, as happens with parallel tool calls.
   */
  struct ExplicitRoiPreparedFit
  {
    std::shared_ptr<ExplicitRoiSetup> setup;
    std::shared_ptr<ExplicitRoiFitResult> fit;
  };

  /** Validates the request, and works out what to do with existing peaks (GUI thread); throws on error. */
  InterSpec_API ExplicitRoiPreparedFit prepare_explicit_roi_fit( const FitPeakOptions &options, InterSpec *interspec );

  /** Fits the peaks; touches no session state, so may be run on any thread.  Throws on error. */
  InterSpec_API void run_explicit_roi_fit( ExplicitRoiPreparedFit &prepared );

  /** Updates the analysis peaks (GUI thread); throws if the peaks it depends on have changed. */
  InterSpec_API FitPeakStatus commit_explicit_roi_fit( const ExplicitRoiPreparedFit &prepared, InterSpec *interspec );

  struct GetUserPeakOptions {
    SpecUtils::SpectrumType specType;
    std::optional<double> lowerEnergy;
    std::optional<double> upperEnergy;
  };
  
  struct GetUserPeakStatus {
    std::vector<std::shared_ptr<const PeakDef>> peaks;
  };
  
  InterSpec_API GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec );

  /** Get a list of unique sources assigned to peaks in the specified spectrum.

   Returns a vector of source strings (nuclide symbols, element symbols, or reaction names)
   that are assigned to user peaks in the specified spectrum. Each source appears only once
   in the result, even if multiple peaks are assigned to that source.

   @param specType Which spectrum to get sources from (Foreground, Background, or SecondForeground)
   @param interspec Pointer to the InterSpec session
   @return Vector of unique source strings assigned to peaks
   @throws std::runtime_error if InterSpec is null or spectrum not loaded
   */
  InterSpec_API std::vector<std::string> get_identified_sources( const SpecUtils::SpectrumType specType, InterSpec *interspec );

  struct FitPeaksForNuclideOptions {
    std::vector<std::string> sources;
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
    bool doNotAddPeaksToUserSession = false;
    Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> fitSrcPeaksOptions;
  };
  
  struct FitPeaksForNuclideStatus
  {
    std::vector<std::shared_ptr<const PeakDef>> fitPeaks;
    std::vector<std::shared_ptr<const PeakDef>> peaksToRemove;
    bool peaksWereRemovedFromSession = false;
    bool noObservablePeaksFit = false;
    std::vector<std::string> warnings;
  };
  
  InterSpec_API FitPeaksForNuclideStatus fit_peaks_for_nuclides( const FitPeaksForNuclideOptions &options, InterSpec *interspec );

  /** Async version of fit_peaks_for_nuclides that runs the expensive fitting in a background thread.
   The callback is called on the GUI thread with either a FitPeaksForNuclideStatus (success) or a string (error message).
   */
  using FitPeaksCallback = std::function<void( std::variant<FitPeaksForNuclideStatus, std::string> )>;
  InterSpec_API void fit_peaks_for_nuclides_async( const FitPeaksForNuclideOptions &options, InterSpec *interspec,
                                                    FitPeaksCallback callback );

  /** Calculates an approximate importance that a peak in a spectrum will have for a given nuclide.
   * 
   * Importance is defined by `yield(i)*sqrt(energy(i))/sum(yield*sqrt(energy))`.
   * A very rough approximation, but can be useful for matching things.
   * This is what PhotoPeak.lis uses for its fourth column.
   * 
   * @param gamma_energies_and_yields A vector of tuples containing the energy and yield of each gamma.
   * @return A tuple containing the energy, yield, and estimated importance.
   */
  std::vector<std::tuple<float,float,float>> cacl_estimated_gamma_importance( const std::vector<std::tuple<float,float>> &gamma_energies_and_yields );

  /** Returns up to three of the most prominent photon energies of a nuclide, x-ray element, or
   reaction; empty if it emits no photons (e.g., a pure beta emitter such as H3 or P32).
   */
  InterSpec_API std::vector<float> get_characteristic_gammas( const std::string &nuclide );

  /** Fewest 511 keV annihilation photons per decay for a nuclide to count as a positron emitter -
   for its characteristic lines (#get_characteristic_gammas), and for 511 keV lookups. */
  constexpr double sm_min_annihilation_photons_per_decay = 0.1;

  /** A positron emitter whose annihilation photons are a significant 511 keV source. */
  struct AnnihilationSource
  {
    const SandiaDecay::Nuclide *nuclide = nullptr;
    double photons511PerDecay = 0.0;  // 2 x positrons per decay of the parent, progeny included
  };//struct AnnihilationSource

  /** Positron emitters a 511 keV annihilation peak may come from: the nuclide-search category
   nuclides, the characteristic-gamma nuclides, and common beta+ emitters (Na22, F18, Ga68, Ge68,
   Cu64, ...), that give at least `min_photons_per_decay` 511 keV photons per decay, at their default
   age; sorted by that yield, largest first.
   */
  InterSpec_API std::vector<AnnihilationSource> positron_emitters_for_annihilation( InterSpec *interspec,
                                      const double min_photons_per_decay = sm_min_annihilation_photons_per_decay );
  
  /** Returns the Nuclides, x-rays, and reaction gammas from CharacteristicGammas.txt, in the energy range. */
  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_nuclides_with_characteristics_in_energy_range( double lower_energy, double upper_energy, InterSpec *interspec );
  
  /** Gets the approximate FWHM of the foreground spectrum for the provided energy, and returns
   `get_nuclides_with_characteristics_in_energy_range(energy-fwhm, energy+fwhm)`
   
   Throws exception on error (no foreground, invalid InterSpec pointer, etc).
   */
  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_characteristics_near_energy( const double energy, InterSpec *interspec );

  /** Get the expected Full Width at Half Maximum (FWHM) for a peak at the specified energy.
   
   Uses the detector response function if available, otherwise fits detected peaks to estimate
   resolution, or falls back to detector type-based estimation.
   
   @param energy The energy (in keV) to calculate the expected FWHM for
   @param interspec Pointer to the InterSpec session
   @param source If non-null, set to a short description of where the FWHM came from.
   @return The expected FWHM in keV
   @throws std::runtime_error if no InterSpec session, no foreground loaded, energy not positive,
           or FWHM cannot be determined
   */
  InterSpec_API float get_expected_fwhm( const double energy, InterSpec *interspec, std::string *source = nullptr );

  /** Returns a callable giving the expected FWHM (keV) at an energy (keV).

   The resolution source is resolved once, up front, so the returned function is cheap to evaluate
   repeatedly, and does not touch InterSpec session state when called.  In order of preference:
   - the loaded DRF's resolution, unless the spectrum's own peaks show it is off by more than ~33%
     (a generic DRF on a better or worse detector);
   - a resolution curve fit to the detected peaks, anchored to the shape of the DRF's (or, without
     one, the detector class's) curve, so a few peaks still give a sensible curve;
   - the nearest detected peak (within 20% of the requested energy);
   - a generic HPGe or low-resolution estimate.
   Gets the detected peaks (which may run the automated peak search, once per spectrum).

   @param interspec Pointer to the InterSpec session
   @param source If non-null, set to a short description of where the FWHM comes from.
   @return Function mapping energy (keV) to expected FWHM (keV)
   @throws std::runtime_error if no InterSpec session or no foreground loaded
   */
  InterSpec_API std::function<float(double)> make_expected_fwhm_fcn( InterSpec *interspec, std::string *source = nullptr );

  /** Worker-safe implementation of `make_expected_fwhm_fcn(InterSpec *)`.

   @param drf The loaded detector response; used if it has resolution info consistent with `peaks`.  May be null.
   @param peaks Peaks of `spectrum` to check the DRF against, or fit a resolution to.  May be null or empty.
   @param spectrum The spectrum the peaks are from.
   @param det_type The detector class; selects the resolution shape and the final generic fallback.
   @param source If non-null, set to a short description of where the FWHM comes from.
   */
  InterSpec_API std::function<float(double)> make_expected_fwhm_fcn( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                     const std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &peaks,
                                     const std::shared_ptr<const SpecUtils::Measurement> &spectrum,
                                     const PeakFitUtils::CoarseResolutionType det_type,
                                     std::string *source = nullptr );

  struct SpectrumCountsInEnergyRange
  {
    double lower_energy;  // actual lower energy, snapped to channel boundary
    double upper_energy;  // actual upper energy, snapped to channel boundary

    double foreground_counts;  // total counts in range

    /** Counts divided by live time.  If live time of spectrum is not valid, this will be set to NaN. */
    double foreground_cps;

    // Per-channel data (empty when max_channels <= 1)
    std::vector<float> channel_lower_energies;    // lower edge of each (possibly combined) bin
    std::vector<float> foreground_channel_counts;  // counts in each bin
    size_t channel_combine_factor = 1;             // the 2^N combining factor used

    struct BackgroundSecondaryInfo
    {
      double counts;

      /** Counts divided by live time.  If live time of spectrum is not valid, this will be set to NaN. */
      double cps;

      /** The number of statistical sigma the foreground is elevated,
       * relative to this spectrum (the background or secondary), when live-time normalized.
       * Positive numbers indicate foreground is elevated, negative numbers indicate
       * deficit in foreground.
       */
      double num_sigma_rel_foreground;

      /** Scale factor: foreground_live_time / this_live_time.
       * Set to 1.0 if either live time is zero or invalid.
       */
      double live_time_scale_factor;

      /** Per-channel counts, rebinned to foreground energy calibration and scaled by live-time ratio.
       * Empty when max_channels <= 1.
       */
      std::vector<float> channel_counts;
    };//struct BackgroundSecondaryInfo

    std::optional<BackgroundSecondaryInfo> background_info;
    std::optional<BackgroundSecondaryInfo> secondary_info;
  };//struct SpectrumCountsInEnergyRange

  SpectrumCountsInEnergyRange get_counts_in_energy_range( double lower_energy, double upper_energy,
                                                          int max_channels, InterSpec *interspec );

  /** Structural actions for peak editing that are mutually exclusive with property changes. */
  enum class EditPeakStructuralAction
  {
    DeletePeak,
    SplitFromRoi,
    MergeWithLeft,
    MergeWithRight
  };//enum class EditPeakStructuralAction

  const char *to_string( EditPeakStructuralAction action );
  EditPeakStructuralAction structural_action_from_string( const std::string &str );


  /** Options for editing a peak at a specified energy.

   Supports setting multiple peak properties in a single call.  When a structural action is
   specified (Delete, Split, Merge), all property fields are ignored.

   Shape-affecting property changes (mean, fwhm, amplitude, skew type/params, continuum
   type/coefficients, ROI bounds) trigger an automatic refit unless suppressed via #refit.
   Fit-for flags and metadata changes do not trigger auto-refit.

   Important: setting a parameter value (e.g., skewPar0) triggers a refit that may adjust
   that value.  To lock a value, also set its corresponding fitFor flag to false.
   */
  struct EditAnalysisPeakOptions
  {
    /** Energy (keV) used to find the nearest analysis peak (must be within 1 FWHM). */
    double energy;
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;

    /** When set, performs a structural action and all property fields below are ignored. */
    std::optional<EditPeakStructuralAction> structuralAction;

    // --- Gaussian peak properties ---
    std::optional<double> meanEnergy;
    std::optional<double> fwhm;
    std::optional<double> amplitude;

    std::optional<double> meanUncertainty;
    std::optional<double> fwhmUncertainty;
    std::optional<double> amplitudeUncertainty;

    std::optional<bool> fitForMean;
    std::optional<bool> fitForFwhm;       // maps to PeakDef::CoefficientType::Sigma
    std::optional<bool> fitForAmplitude;

    // --- Skew ---
    std::optional<std::string> skewType;
    std::optional<double> skewPar0;
    std::optional<double> skewPar1;
    std::optional<double> skewPar2;
    std::optional<double> skewPar3;
    std::optional<double> skewPar4;
    std::optional<double> skewPar5;
    std::optional<bool> fitForSkewPar0;
    std::optional<bool> fitForSkewPar1;
    std::optional<bool> fitForSkewPar2;
    std::optional<bool> fitForSkewPar3;
    std::optional<bool> fitForSkewPar4;
    std::optional<bool> fitForSkewPar5;

    // --- Continuum ---
    std::optional<std::string> continuumType;
    std::optional<double> continuumCoef0;
    std::optional<double> continuumCoef1;
    std::optional<double> continuumCoef2;
    std::optional<double> continuumCoef3;
    std::optional<bool> fitForContinuumCoef0;
    std::optional<bool> fitForContinuumCoef1;
    std::optional<bool> fitForContinuumCoef2;
    std::optional<bool> fitForContinuumCoef3;

    // --- ROI bounds ---
    std::optional<double> roiLowerEnergy;
    std::optional<double> roiUpperEnergy;

    // --- Source / metadata ---
    std::optional<std::string> source;
    std::optional<std::string> color;
    std::optional<std::string> userLabel;

    // --- Usage flags ---
    std::optional<bool> useForEnergyCalibration;
    std::optional<bool> useForShieldingSourceFit;
    std::optional<bool> useForManualRelEff;

    /** Refit control.
     nullopt: auto (refit if shape-affecting values changed).
     true: force refit even if only metadata changed.
     false: suppress refit even if shape values changed.
     */
    std::optional<bool> refit;
  };//struct EditAnalysisPeakOptions

  /** Results of peak edit operation. */
  struct EditAnalysisPeakStatus
  {
    bool success;
    std::string message;
    bool refitPerformed = false;
    std::optional<std::shared_ptr<const PeakDef>> modifiedPeak;
    std::vector<std::shared_ptr<const PeakDef>> peaksInRoi;
  };//struct EditAnalysisPeakStatus

  /** Edit or delete a peak at the specified energy.

   Supports modifying multiple peak properties in a single call, with optional automatic
   refitting when shape-affecting properties change.

   @param options Specifies the energy and property changes (or structural action) to apply
   @param interspec Pointer to the InterSpec session containing the peak
   @return EditAnalysisPeakStatus containing success/failure, refit info, and modified peak
   */
  InterSpec_API EditAnalysisPeakStatus edit_analysis_peak( const EditAnalysisPeakOptions &options, InterSpec *interspec );

  /** Enum specifying the type of escape peak. */
  enum class EscapePeakType
  {
    SingleEscape,  // Parent energy minus 511 keV (one positron escapes)
    DoubleEscape   // Parent energy minus 1022 keV (both positrons escape)
  };//enum class EscapePeakType

  /** Converts EscapePeakType enum to string representation.

   @param type The escape peak type to convert
   @return String representation of the type
   */
  const char* to_string( EscapePeakType type );

  /** Information about the parent peak if the queried energy is an escape peak. */
  struct ParentPeakInfo
  {
    double parentPeakEnergy;
    EscapePeakType escapeType;
    std::optional<std::string> sourceLabel;  // e.g., "Th232 S.E. 2614.52 keV" or "Co56 D.E. 3451.15 keV"
    std::optional<std::string> userLabel;    // e.g., "2614.53 S.E." if no source assigned
  };//struct ParentPeakInfo

  /** Options for checking if a peak is an escape peak. */
  struct EscapePeakCheckOptions
  {
    double energy = 0.0;
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
  };//struct EscapePeakCheckOptions

  /** A peak found at an energy the escape check looked at. */
  struct EscapeCheckPeak
  {
    std::shared_ptr<const PeakDef> peak;
    std::optional<std::string> source;  // the peak's assigned nuclide, x-ray element, or reaction
  };//struct EscapeCheckPeak

  /** A possible detector-material K x-ray escape relation between the queried peak and a
   higher-energy line: a photon photo-absorbed near the crystal surface, whose fluorescence
   K x-ray escapes, deposits (parent energy - x-ray energy).
   */
  struct XrayEscapeInfo
  {
    std::string element;          // detector element the x-ray comes from, e.g., "Ge", "I", "Cd"
    std::string line;             // "Ka" or "Kb"
    double xrayEnergy = 0.0;      // keV
    double kEdgeEnergy = 0.0;     // keV; a parent must be above this to produce the x-ray
    double parentEnergy = 0.0;    // queried energy + xrayEnergy: where the parent line would be

    /** A peak found at `parentEnergy`, if any. */
    std::optional<EscapeCheckPeak> parentPeak;

    /** Queried peak area / parent peak area; only set if both peaks exist. */
    std::optional<double> areaRatio;

    /** A generous upper bound on (escape peak area / parent peak area) for this line, at the
     parent energy; zero if the parent is below the K edge. */
    double maxEscapeFraction = 0.0;

    /** True if a parent peak exists, and `areaRatio` is no larger than `maxEscapeFraction`. */
    bool plausible = false;
  };//struct XrayEscapeInfo

  /** Results of escape peak analysis. */
  struct EscapePeakCheckStatus
  {
    // If the specified peak IS an escape peak of another peak in the spectrum
    std::optional<ParentPeakInfo> parentPeak;

    // Where this peak's own escape peaks would be, if this peak were a parent; only set when the
    //  energy is above the 1022 keV pair-production threshold, and the escape energy is inside
    //  the spectrum.
    std::optional<double> potentialSingleEscapePeakEnergy;     // input - 511 keV
    std::optional<double> potentialDoubleEscapePeakEnergy;     // input - 1022 keV

    // Peaks found at those energies, if any.
    std::optional<EscapeCheckPeak> singleEscapePeakFound;
    std::optional<EscapeCheckPeak> doubleEscapePeakFound;

    // Where the parent would be if this peak were its single or double escape peak (always
    //  above 1022 keV for a double escape; for a single escape, below 1022 keV when the input is
    //  below 511 keV - in which case no such parent can pair-produce).
    double potentialParentPeakSingleEscape = 0.0;     // input + 511 keV
    double potentialParentPeakDoubleEscape = 0.0;     // input + 1022 keV

    // Peaks found at those energies, whether or not they can be the parent.
    std::optional<EscapeCheckPeak> peakAtPotentialSingleEscapeParent;
    std::optional<EscapeCheckPeak> peakAtPotentialDoubleEscapeParent;

    /** True when the input energy is below 1022 keV, so it can have no escape peaks of its own. */
    bool belowPairProductionThreshold = false;

    /** Detector material the x-ray escape checks assume, e.g., "Ge (HPGe)", "NaI", or "unknown
     (all common materials checked)". */
    std::string detectorMaterialAssumed;

    /** Detector K x-ray escape relations of the queried peak to higher-energy lines. */
    std::vector<XrayEscapeInfo> xrayEscapes;

    /** Plain-language notes, e.g., that a would-be parent is below the pair-production threshold. */
    std::vector<std::string> notes;

    // Search windows (±keV) for finding peaks at those energies
    // Calculated as max(1.0, min(0.5*fwhm, 15.0)) - fairly arbitrary and untested
    double singleEscapeSearchWindow = 0.0;            // window for finding S.E. peak of input energy
    double doubleEscapeSearchWindow = 0.0;            // window for finding D.E. peak of input energy
    double singleEscapeParentSearchWindow = 0.0;      // window for finding parent if input is S.E.
    double doubleEscapeParentSearchWindow = 0.0;      // window for finding parent if input is D.E.
  };//struct EscapePeakCheckStatus

  /** Check if a peak at a specified energy is a single or double escape peak.

   This function analyzes whether the peak at the given energy is an escape peak of another
   peak in the spectrum. It checks for parent peaks at energy+511 keV (single escape) and
   energy+1022 keV (double escape), verifying that parent peaks are above the pair production
   threshold (~1255 keV for HPGe). The function also handles the edge case where the parent
   peak is above the detector's energy range by checking nuclide photon energies.

   Note: The function checks for double escape parents (1022 keV above) BEFORE single escape
   parents (511 keV above) to avoid misidentifying a double escape peak as a single escape.

   @param options Specifies the energy to check and which spectrum to use
   @param interspec Pointer to the InterSpec session containing the spectrum
   @return EscapePeakCheckStatus containing parent peak info (if found) and calculated potential energies
   @throws std::runtime_error if InterSpec is null, no spectrum loaded, or other errors occur
   */
  InterSpec_API EscapePeakCheckStatus escape_peak_check( const EscapePeakCheckOptions &options, InterSpec *interspec );

  /** Enum specifying the type of sum peak. */
  enum class SumPeakType
  {
    NotASumPeak,          // Not identified as a sum peak
    RandomSum,            // Sum from random coincidence (pile-up)
    CascadeSum,           // Sum from cascade decay photons emitted nearly simultaneously
    CascadeAndRandomSum,  // Both: a cascade pair, with random summing a sizable part of the rate
    Unknown               // Could be a sum peak but type cannot be determined
  };//enum class SumPeakType

  /** Converts SumPeakType enum to string representation.

   @param type The sum peak type to convert
   @return String representation of the type
   */
  const char* to_string( SumPeakType type );

  /** Information about a sum peak if the queried energy is identified as one. */
  struct SumPeakInfo
  {
    SumPeakType sumType;
    std::shared_ptr<const PeakDef> firstPeak;   // First contributing peak
    std::shared_ptr<const PeakDef> secondPeak;  // Second contributing peak
    std::optional<std::string> userLabel;       // e.g., "Sum 661.66 + 661.66 keV" or "Co60 cascade-sum 1173.23+1332.49 keV"

    /** For cascade sums: probability per decay that both photons are emitted (the two photon
     yields times their coincidence probability), for the nuclide in `cascadeNuclide`. */
    std::optional<double> coincidenceFraction;
    std::optional<std::string> cascadeNuclide;
    std::optional<std::string> coincidenceType;  // "gamma-gamma", "x-ray-gamma", "x-ray-x-ray", or "annihilation-gamma"

    /** True if the pair is a cascade, but random summing alone explains the observed rate. */
    bool cascadeAlsoPossible = false;

    /** Count rate of the sum peak (net, per live-time second), if a peak is there. */
    std::optional<double> observedSumCps;

    /** Random-summing rate expected for this pair: tau*R^2 for a line summing with itself,
     2*tau*R1*R2 for two different lines; see `SumPeakCheckStatus::pileupTau`.  When tau is only
     an upper bound, so is this. */
    std::optional<double> expectedRandomSumCps;

    /** Why this classification was made, for an LLM or user to read. */
    std::string reason;
  };//struct SumPeakInfo

  /** Options for checking if a peak is a sum peak. */
  struct SumPeakCheckOptions
  {
    double energy = 0.0;
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
    std::optional<double> distance;     // Optional: source-to-detector distance (in PhysicalUnits, e.g., PhysicalUnits::cm). If > 15 cm, cascade-sums are not considered
  };//struct SumPeakCheckOptions

  /** Results of sum peak analysis. */
  struct SumPeakCheckStatus
  {
    // If the specified peak IS a sum peak
    std::optional<SumPeakInfo> sumPeakInfo;

    // Energy windows used for searching for contributing peaks
    double searchWindow = 0.0;  // calculated as max(1.0, min(0.5*fwhm, 15.0))

    /** Pile-up (random summing) time constant, in seconds, so a line of rate R sums with itself at
     tau*R^2, and two lines at 2*tau*R1*R2.  Measured from self-sum peaks of strong lines when there
     are any (other than the queried one), otherwise only an upper bound from the dead time. */
    std::optional<double> pileupTau;
    bool pileupTauIsUpperBound = false;
    std::string pileupTauSource;  // e.g., "558.30 keV = 2 x 279.15 keV self-sum peak", or "dead time (upper bound)"

    /** Other sources with a characteristic line at about the queried energy, which could explain
     the peak instead of summing. */
    std::vector<std::string> otherSourcesAtThisEnergy;
  };//struct SumPeakCheckStatus

  /** Check if a peak at a specified energy is a sum peak (random-sum or cascade-sum).

   This function analyzes whether the peak at the given energy is the result of summing
   two lower-energy photons. There are two types of sum peaks:

   1. Random-sum peaks: Occur when unrelated photons arrive at the detector nearly
      simultaneously, typically when dead time is high (>20% for HPGe, >10% for low-res).
      The function validates by checking if the highest-amplitude peaks produce expected
      sum peaks.

   2. Cascade-sum peaks: Occur when a nuclide emits two photons in quick succession and
      the detector is close enough (<~15 cm) to detect both. Both contributing peaks must
      have the same nuclide source assigned, and the nuclide's decay scheme must list
      those gammas as cascade coincidences. If a distance parameter is specified and is
      greater than 15 cm, cascade-sum detection is disabled.

   When multiple candidate peak pairs are found, the function selects the best match by:
   - First prioritizing the closest energy match to the target
   - If multiple pairs have similar energies (within ~0.5 keV), selecting the pair with
     the largest combined area (for cascade sums: area1 × area2 × coincidenceFraction)

   @param options Specifies the energy to check, which spectrum to use, and optionally
                  the source-to-detector distance (in PhysicalUnits)
   @param interspec Pointer to the InterSpec session containing the spectrum
   @return SumPeakCheckStatus containing sum peak info (if identified) and search parameters
   @throws std::runtime_error if InterSpec is null, no spectrum loaded, or other errors occur
   */
  InterSpec_API SumPeakCheckStatus sum_peak_check( const SumPeakCheckOptions &options, InterSpec *interspec );


  /** Confidence that a candidate peak is a Compton-scatter feature of another peak. */
  enum class ComptonScatterConfidence
  {
    No,
    Unlikely,
    Likely,
    VeryLikely
  };//enum class ComptonScatterConfidence

  const char* to_string( ComptonScatterConfidence c );

  /** Type of Compton-scatter feature identified.

   ComptonPeak corresponds to a "backscatter peak" in standard gamma-spectroscopy
   literature (Knoll Ch. 10.2; Gilmore §2.2.3): a photon scatters off material
   in front of the detector at angle theta and the reduced-energy photon
   photopeaks in the detector at E' = E_p / ( 1 + (E_p/m_e c^2)(1 - cos theta) ).
   Distinct peaks form near 180 deg due to kinematic focusing.

   ComptonEdge is the maximum energy *deposited* in the detector when an
   incoming photon Compton-scatters near 180 deg inside the crystal and the
   scattered photon escapes: E_edge = E_p * 2 alpha / ( 1 + 2 alpha ), with
   alpha = E_p / m_e c^2.  The edge is a step/shoulder rather than a true peak,
   so a Gaussian fit on the edge is necessarily an artifact.

   ComptonPeakOrEdge is used near E_p ~ 256 keV (m_e c^2 / 2), where the
   backscatter peak and the Compton-edge fall close together (E_bs + E_edge = E_p).
   */
  enum class ComptonScatterType
  {
    None,
    ComptonPeak,
    ComptonEdge,
    ComptonPeakOrEdge
  };//enum class ComptonScatterType

  const char* to_string( ComptonScatterType t );

  /** Information about the matched parent peak when a candidate peak is
   identified as a Compton-scatter feature.  All numeric fields are zero and
   userLabel is empty when no scatter match was found.
   */
  struct ComptonScatterParentInfo
  {
    std::shared_ptr<const PeakDef> parentPeak;
    double parentEnergy = 0.0;   // keV
    double parentFwhm = 0.0;     // keV
    double scatterAngle = 0.0;   // degrees; 0 for pure ComptonEdge matches
    double areaFraction = 0.0;   // candidateArea / parentArea
    std::optional<std::string> userLabel;  // e.g. "Compton scatter from Co60 1332.50 keV peak"
  };//struct ComptonScatterParentInfo

  /** Options for compton_peak_check. */
  struct ComptonPeakCheckOptions
  {
    double energy = 0.0;                                                       // keV
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
  };//struct ComptonPeakCheckOptions

  /** A higher-energy peak whose single-Compton-scatter features cover an energy - found from the
   kinematics alone, so no peak needs to have been fit at that energy.

   A photon of energy Ep that scatters by 180 deg in material around the detector, then is absorbed,
   deposits E_bs = Ep / (1 + 2 Ep/511) - the backscatter peak, a broad bump from scatters of ~120-180
   deg.  One that scatters by 180 deg in the detector and escapes deposits E_edge = Ep - E_bs - the
   Compton edge, a step down in the continuum.
   */
  struct ComptonKinematicParent
  {
    std::shared_ptr<const PeakDef> parentPeak;
    double parentEnergy = 0.0;        // keV
    double parentArea = 0.0;          // counts
    double parentCps = -1.0;          // counts per live-time second; -1 if no live time
    std::optional<std::string> parentSource;
    double backscatterEnergy = 0.0;   // E_bs (keV)
    double comptonEdgeEnergy = 0.0;   // E_edge (keV)
    double window = 0.0;              // keV the energy may be outside the backscatter range, or from E_edge

    /** For a backscatter match: the angle (degrees) a photon from the parent must scatter by to have
     the queried energy (180 if the energy is just below E_bs).  -1 for an edge-only match. */
    double scatterAngle = -1.0;

    /** ComptonPeak (backscatter, angle >= ~120 deg), ComptonEdge (within the window of E_edge), or
     ComptonPeakOrEdge (both); ComptonPeak with a smaller angle is less likely to be a visible feature. */
    ComptonScatterType feature = ComptonScatterType::None;

    /** Parent area / area of the largest peak above the queried energy. */
    double prominence = 0.0;

    /** Significance of the excess of the counts in [E - w, E + w] over a straight line
     interpolated between the regions just outside the parent's [E_bs - w, E_edge + w];
     0 if not computable. */
    double excessSigma = 0.0;

    /** How strongly the region alone (`prominence` and `excessSigma`) suggests a scatter feature
     from this parent.  A fitted candidate peak is judged on more criteria (its width and area), so
     the overall verdict can be lower than this.
     */
    ComptonScatterConfidence confidence = ComptonScatterConfidence::No;
  };//struct ComptonKinematicParent

  /** Finds the higher-energy peaks in `peaks` with a scatter feature at `energy`: either the energy
   is in the parent's backscatter range ([E_bs - window, E_bs(120 deg) + window], where scatters of
   ~120-180 degrees bunch up into the backscatter peak), or it is within `window` of the parent's
   Compton edge.  (Scatters at smaller angles spread over the Compton continuum, and make no
   feature.)  Most confident first.

   Only Gaussian peaks above ~50 keV, with significance of at least 5, are considered parents.
   Thread-safe: touches no session state.

   @param window keV beyond E_bs and E_edge that `energy` may be; typically the FWHM at `energy`.
   */
  InterSpec_API std::vector<ComptonKinematicParent> compton_kinematic_parents( const double energy,
                                    const double window,
                                    const std::vector<std::shared_ptr<const PeakDef>> &peaks,
                                    const std::shared_ptr<const SpecUtils::Measurement> &spectrum );

  /** Results of Compton-scatter analysis for a candidate peak. */
  struct ComptonPeakCheckStatus
  {
    ComptonScatterConfidence confidence = ComptonScatterConfidence::No;
    ComptonScatterType type = ComptonScatterType::None;

    /** True if no fitted peak was found within `searchWindow` of the energy; the confidence and
     type then come from `kinematicParents` alone. */
    bool noCandidatePeak = false;

    /** The peak nearest the queried energy, and its distance (keV), if there are any peaks. */
    std::optional<double> nearestPeakEnergy;
    std::optional<double> nearestPeakDistance;

    /** Why the result is what it is, e.g., why a candidate peak was rejected. */
    std::string reason;

    /** Parents whose scatter features cover the energy, from kinematics alone; see
     #compton_kinematic_parents.  Only filled in by the `InterSpec *` overload. */
    std::vector<ComptonKinematicParent> kinematicParents;

    // Populated only when confidence != No.
    std::optional<ComptonScatterParentInfo> parentInfo;

    // Candidate-peak diagnostics; zero if no peak was found near `energy`.
    double candidateMean = 0.0;
    double candidateFwhm = 0.0;

    // FWHM the detector would give at candidateMean.  -1 if the cheap-path
    // (DRF resolution info or cached automatedSearchPeaks) was unavailable
    // and a sqrt(E)-vs-parent FWHM fallback was used instead.
    double expectedDetectorFwhm = 0.0;

    // Parent peak area / continuum integral over its ROI, when a parent matched.
    double parentSignalToContinuum = 0.0;

    // Tolerance used to find candidate peak near input energy (keV)
    double searchWindow = 0.0;
  };//struct ComptonPeakCheckStatus

  /** Check if `candidate` is a Compton-scatter feature (backscatter peak or
   Compton edge) of a higher-energy peak in `user_peaks`.

   Worker-thread-safe overload: takes already-fetched parameters and does not
   touch any InterSpec session state, so it may be called from a worker thread
   (e.g., from inside `IsotopeId::populateCandidateNuclides`).

   For each `user_peaks` entry above `candidate->mean()`, inverts the
   Compton-scattering kinematic for cos(theta) and computes the Compton-edge
   energy; picks the parent with the smallest energy-match residual.  The FWHM
   criterion uses `detector->peakResolutionFWHM(...)` when
   `detector->hasResolutionInfo()` is true, otherwise falls back to a sqrt(E)
   scaling of the matched parent's measured FWHM (HPGe ~ sqrt(a + bE),
   NaI ~ sqrt(E)).

   Confidence is the conservative minimum across five buckets: FWHM ratio,
   area fraction (candidate/parent), parent peak signal-to-continuum, scatter
   angle, and energy-match tightness.  Hard-No guards include parent area
   <= candidate area, area_fraction > 0.5, and fwhm_ratio < 1.15.

   @param candidate   The peak being tested.  Must be a Gaussian-fit user peak.
   @param user_peaks  All user peaks in the spectrum (parents are searched here).
   @param spectrum    The displayed histogram, used to integrate the parent's continuum.
   @param detector    Optional; if non-null and carrying resolution info, used for FWHM check.
   @param expected_fwhm Optional; if given, used for the FWHM check instead of `detector` (e.g.,
                      from #make_expected_fwhm_fcn, which checks the DRF against the spectrum).
   @return ComptonPeakCheckStatus with confidence, scatter type, and parent info.
           confidence == No / type == None when no kinematic match is found or
           any hard-No guard fails; `reason` then says which.
   */
  InterSpec_API ComptonPeakCheckStatus compton_peak_check(
      const std::shared_ptr<const PeakDef> &candidate,
      const std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > &user_peaks,
      const std::shared_ptr<const SpecUtils::Measurement> &spectrum,
      const std::shared_ptr<const DetectorPeakResponse> &detector,
      const std::function<float(double)> &expected_fwhm = nullptr );

  /** Main-thread wrapper around the worker-safe overload above.

   Fetches the user peaks, spectrum, and detector for `options.specType` from
   `interspec`, locates the user peak nearest `options.energy`, and delegates
   to the worker-safe overload.

   @param options    Candidate energy and which spectrum to use.
   @param interspec  Pointer to the InterSpec session containing the spectrum.
   @throws std::runtime_error if InterSpec is null or no spectrum is loaded.
   */
  InterSpec_API ComptonPeakCheckStatus compton_peak_check( const ComptonPeakCheckOptions &options, InterSpec *interspec );


  /** Graded verdict of the pure-beta-emitter (bremsstrahlung) continuum check.

   The check background-subtracts the foreground and analyzes the net continuum shape;
   a smooth quasi-exponential hump with no gamma-peak evidence and a gradual (non-cliff)
   high-energy fade is characteristic of bremsstrahlung from a pure beta emitter
   (e.g. P-32, Sr-90/Y-90, Tl-204).
   */
  enum class BetaContinuumVerdict
  {
    BackgroundNotLoaded,   // check requires a displayed Background spectrum to subtract
    NoContinuumElevation,  // net continuum not significantly above background
    BremLike,              // all gates passed: smooth brem-shaped continuum, no gamma peaks
    Ambiguous,             // elevation present, but a soft gate failed (smoothness marginal,
                           //  or too little tail to analyze)
    NotBremLike,           // a hard gate failed (gamma peaks present, cliff-like termination,
                           //  shape too hard/soft, or termination beyond plausible beta endpoint)
    AnnihilationDominated  // strong 511 keV peak with continuum ending just above it -
                           //  suggests a positron source (e.g. F-18), not a pure beta-minus emitter
  };//enum class BetaContinuumVerdict

  const char *to_string( BetaContinuumVerdict verdict );

  /** A peak the beta-continuum check considered as gamma-source evidence (or an x-ray report). */
  struct BetaContinuumPeakInfo
  {
    double energy = 0.0;             // keV
    double amplitude = 0.0;          // peak area, counts
    double numSigmaElevated = 0.0;   // significance of elevation over the (live-time scaled)
                                     //  matching background peak; over zero if no match
    bool isAnnihilation = false;     // within tolerance of 511 keV
  };//struct BetaContinuumPeakInfo

  /** A pure-beta nuclide whose (equilibrium-chain) endpoint is consistent with the observed
   continuum termination energy.  A consistency candidate - NOT an identification.
   */
  struct BetaContinuumCandidateNuclide
  {
    std::string nuclide;             // e.g. "Sr90"
    double betaEndpoint = 0.0;       // keV; max endpoint over the parent and any shorter-lived
                                     //  (in-equilibrium) descendants
    std::string note;                // e.g. "endpoint is from in-equilibrium daughter Y90"
  };//struct BetaContinuumCandidateNuclide

  /** Options for beta_continuum_check. */
  struct BetaContinuumCheckOptions
  {
    /** Spectrum to test as the candidate beta/bremsstrahlung source; Foreground or
     SecondForeground only (the displayed Background is always the subtraction baseline). */
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
  };//struct BetaContinuumCheckOptions

  /** Worker-thread-safe inputs for beta_continuum_check - no InterSpec session state. */
  struct BetaContinuumCheckInput
  {
    std::shared_ptr<const SpecUtils::Measurement> foreground;    // required
    std::shared_ptr<const SpecUtils::Measurement> background;    // null => BackgroundNotLoaded
    std::vector<std::shared_ptr<const PeakDef>> foregroundPeaks; // user + auto-search peaks
    std::vector<std::shared_ptr<const PeakDef>> backgroundPeaks; // user + auto-search peaks
    std::function<float(double)> expectedFwhm;                   // FWHM(energy in keV) in keV; required
  };//struct BetaContinuumCheckInput

  /** Results of the beta/bremsstrahlung continuum analysis. */
  struct BetaContinuumCheckStatus
  {
    BetaContinuumVerdict verdict = BetaContinuumVerdict::BackgroundNotLoaded;

    /** Total net-continuum significance (sigma) over the analysis range, peak regions excluded. */
    double continuumElevationSigma = 0.0;

    // Shape metrics; populated only once continuum elevation was established:
    std::optional<double> modeEnergy;           // keV; maximum of the smoothed net continuum (turn-on hump)
    std::optional<double> terminationEnergy;    // keV; last coarse bin > 3 sigma.  A statistics-limited
                                                //  LOWER BOUND on the beta endpoint energy.
    std::optional<double> terminationRatio;     // smoothed net at termination / smoothed max; small
                                                //  (<~0.05) for a gradual brem fade, large for a cliff
    std::optional<double> characteristicEnergy; // keV; -1/(dln(net)/dE) at fit midrange (continuum softness)
    std::optional<double> fitChi2Dof;           // robust log-fit chi2/dof (3% systematic floor included)
    std::optional<double> fitOutlierFraction;   // fraction of fit bins rejected as outliers
    std::optional<double> fitMaxSegmentResidual;// max |mean residual| over 5 fit segments, in sigma

    /** Peaks >~120 keV (non-511) significantly elevated over background - gamma-source evidence
     that vetoes a BremLike verdict. */
    std::vector<BetaContinuumPeakInfo> elevatedGammaPeaks;

    /** Peaks at <~120 keV - candidate fluorescence/decay x-rays (e.g. Hg K x-rays from Tl-204's
     EC branch).  Reported for context; these do NOT veto a BremLike verdict. */
    std::vector<BetaContinuumPeakInfo> xrayRegionPeaks;

    /** True when a strong 511 keV peak dominates and the continuum terminates just above it. */
    bool annihilationDominant = false;

    // Set when a weak net tail (beyond the count-quantile used for shape fitting - i.e., holding
    // <~1% of the net counts) extends past terminationEnergy.  Commonly random-summing (pile-up)
    // at elevated count rates; not used for shape grading or the endpoint bound.
    std::optional<double> highEnergyTailUpperEnergy;  // keV
    std::optional<double> highEnergyTailSigma;        // significance of the tail region net counts

    /** Pure-beta nuclides consistent with terminationEnergy; only filled for
     BremLike/Ambiguous/AnnihilationDominated verdicts. */
    std::vector<BetaContinuumCandidateNuclide> candidateNuclides;

    std::vector<std::string> warnings;
  };//struct BetaContinuumCheckStatus

  /** Check whether the net (background-subtracted) continuum is consistent with bremsstrahlung
   from a pure beta emitter.

   Worker-thread-safe overload: pure function of the input struct; touches no InterSpec session
   state.  See BetaContinuumVerdict and BetaContinuumCheckStatus for the semantics of the result.

   @param input Spectra, peak lists, and FWHM function; see BetaContinuumCheckInput.
   @return BetaContinuumCheckStatus with graded verdict and shape metrics.
   @throws std::runtime_error if the foreground or FWHM function is missing/unusable.
   */
  InterSpec_API BetaContinuumCheckStatus beta_continuum_check( const BetaContinuumCheckInput &input );

  /** Main-thread wrapper around the worker-safe overload above.

   Gathers the displayed foreground (or second foreground) and background spectra, their
   user + automated-search peaks (running background-peak recovery first), and the expected-FWHM
   function from the session, then delegates.

   @param options Which spectrum to test (Foreground or SecondForeground).
   @param interspec Pointer to the InterSpec session.
   @throws std::runtime_error if InterSpec is null or the candidate spectrum is not loaded.
   */
  InterSpec_API BetaContinuumCheckStatus beta_continuum_check( const BetaContinuumCheckOptions &options,
                                                               InterSpec *interspec );

} // namespace AnalystChecks

#endif // AnalystChecks_h 
