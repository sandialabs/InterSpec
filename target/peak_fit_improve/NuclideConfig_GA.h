#ifndef NuclideConfig_GA_h
#define NuclideConfig_GA_h
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

#include <memory>
#include <string>
#include <vector>
#include <functional>

#define OPENGA_EXTERN_LOCAL_VARS 1

#include "openGA.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/FitPeaksForNuclides.h"

#include "PeakFitImprove.h"

namespace SpecUtils
{
  class Measurement;
}

class DetectorPeakResponse;


namespace NuclideConfig_GA
{

using PeakFitResult = FitPeaksForNuclides::PeakFitResult;
using PeakFitForNuclideConfig = FitPeaksForNuclides::PeakFitForNuclideConfig;


/** How background is handled during the fit. */
enum class BackgroundMode : int
{
  /** Pass long_background to fit_peaks_for_nuclides for background subtraction. */
  BackgroundSubtracted,

  /** Pass nullptr as background - no subtraction. */
  NoBackground,

  /** Pass nullptr as background, but fit NORM peaks simultaneously using a second rel-eff curve. */
  NoBackgroundFitNorm
};//enum class BackgroundMode

/** Compile-time selection of which background mode to use for the GA optimization.
 Change this and rebuild to switch modes.
 */
constexpr BackgroundMode sm_background_mode = BackgroundMode::BackgroundSubtracted;

/** Weight applied to the background-fit-false-positive penalty before it is
 added to a per-spectrum's combined score.  0.25 makes background false
 positives ~1/4 as important as foreground accuracy. */
constexpr double sm_background_fit_penalty_weight = 0.25;

/** How the manual-rel-eff chi2/dof cap is treated by the GA. */
enum class RelEffChi2CapMode : int { Disabled, Fixed, GAOptimized };

/** Whether to run the bg-fit-false-positive trial in addition to the
 foreground fit.  Set from CLI in main(). */
extern bool sm_do_background_fit_trial;

/** Whether tuning/evaluation/report fits pass DisableAutoInterfererFit.  Set from the
 `--disable-auto-interferer-fit` CLI switch; false preserves the production default. */
extern bool sm_disable_auto_interferer_fit;

/** Mode for `PeakFitForNuclideConfig::initial_manual_rel_eff_max_chi2_dof`.
 Disabled    - cap is set to a sentinel large value (gate never fires)
 Fixed       - cap = sm_rel_eff_chi2_cap_fixed_value
 GAOptimized - cap is a gene the GA optimizes (range 5..1000 log-uniform) */
extern RelEffChi2CapMode sm_rel_eff_chi2_cap_mode;

/** Cap value used when sm_rel_eff_chi2_cap_mode == Fixed. */
extern double sm_rel_eff_chi2_cap_fixed_value;

/** Detector-resolution class this GA run is optimizing for.  `genes_to_settings` starts from
 `PeakFitForNuclideConfig::default_config(sm_base_det_type)` so that the non-GA-optimized fields
 (e.g. skew_type) match the production default for this detector type.  Set from the CLI-selected
 det-type in main() before `do_nuclide_config_ga(...)`.  Defaults to Low (non-HPGe). */
extern PeakFitUtils::CoarseResolutionType sm_base_det_type;

/** Checkpoint/resume for NuclideConfigGA (all opt-in via the CLI).

 `sm_checkpoint_name`: when non-empty (--checkpoint <name>), the full population's genes are written
 every generation to `<sm_output_file_prefix><name>_nuclide_config_ga_checkpoint.tsv` (atomic
 temp+rename).  The det-type prefix plus the user-supplied name keep concurrent runs from
 clobbering one another.

 `sm_resume_path`: when non-empty (--resume <file>), the GA's initial population is seeded from this
 checkpoint (genes only; costs are recomputed; the RNG/stall-counter are not restored).

 `sm_checkpoint_options_summary`: one-line summary of the options that affect the objective and gene
 meaning (det-type, dataset filters, chi2-cap mode/value, background-fit-trial).  Written into the
 checkpoint header; on resume it is compared against the current run so an incompatible continuation
 is refused.  Thread counts are intentionally excluded as they do not affect results. */
extern std::string sm_checkpoint_name;
extern std::string sm_resume_path;
extern std::string sm_checkpoint_options_summary;

/** Cap on the raw (pre-weight) per-spectrum background-fit penalty so a
 single spectrum's false-positive cluster (or fit failure) cannot dominate
 total_score.  30.0 is roughly the typical magnitude of a per-spectrum
 foreground score contribution; combined with the 0.25 weight, this bounds
 the per-spectrum penalty contribution at 7.5. */
constexpr double sm_background_fit_penalty_per_spectrum_cap = 30.0;

/** Penalty assigned to a single spectrum's foreground cost when the fit fails (non-Success status,
 or throws).

 A "failure" is dominated by GRACEFUL EMPTIES: on the NaI eval-holdout, 40 of 49 failures are
 FailedToSetupProblem "No ROIs are defined" - the exact same zero-peaks-recovered OUTCOME as a
 Success that returns 0 observable peaks (which scores at most the total-miss ceiling,
 sm_miss_penalty_weight).  The old value (100) put failures on a wildly different scale: a handful of
 status flips dominated the summed objective and set the A/B noise floor to ~100 per flip - ~40x any
 real algorithm effect (e.g. R1.1 moved the NaI total by ~2.4) - so measuring fitter changes against
 it was impossible.  It is now scored just above the total-miss ceiling (a small premium so an
 outright error is still ranked slightly worse than a graceful empty).  Failure RATE is not lost:
 the per-spectrum TSV emits the status enum, so a failure-prone config stays visible without a cliff
 that swamps the accuracy signal (the two-axis view of GrossScoring_Review.md R2).
 TODO: normalize per-spectrum costs by expected-peak count so multi-line spectra do not dominate the
 summed objective; the failure penalty would then move to that normalized scale.
 [gross-scoring review R2, 2026-07-18] */
constexpr double sm_fit_failure_penalty = 7.5;  // ~= 1.5 x total-miss ceiling (sm_miss_penalty_weight)


/** Weight for the per-spectrum "missed expected area" penalty term, added to the foreground cost:
 `sm_miss_penalty_weight * (fraction of definitely-wanted expected peak-area not detected)`.

 The fraction is in [0,1], so this is the maximum per-spectrum penalty (a total miss).  It exists
 because the rest of the objective is asymmetric: finding a peak is rewarded and its area error is
 penalized, but a *missed* peak only forgoes the reward (its area error is never counted), so the GA
 can lower its cost by dropping hard-to-fit peaks and total-miss spectra score ~0.  Area-weighting
 makes missing a dominant peak hurt far more than a marginal one.  MUST be a fixed constant, not a GA
 gene (the GA would drive its own objective weight to zero).  See missed_def_wanted_area_fraction(). */
constexpr double sm_miss_penalty_weight = 5.0;


/** Holds precomputed data for each spectrum, so the expensive search_for_peaks call is done once. */
struct PrecomputedNuclideData
{
  /** Pointer to the original DataSrcInfo (not owned). */
  const DataSrcInfo *src_info = nullptr;

  /** The source nuclides/elements for this entry. */
  std::vector<RelActCalcAuto::SrcVariant> sources;

  /** Precomputed auto-search peaks (the expensive call). */
  std::vector<std::shared_ptr<const PeakDef>> auto_search_peaks;

  /** The foreground spectrum. */
  std::shared_ptr<const SpecUtils::Measurement> foreground;

  /** Background (may be nullptr depending on BackgroundMode). */
  std::shared_ptr<const SpecUtils::Measurement> background;

  /** Precomputed auto-search peaks on the long_background spectrum, for the
   GA background-false-positive penalty path.  Empty when all of `sources`
   are NORM-like (per `is_norm_like_for_ga`) or when long_background is null. */
  std::vector<std::shared_ptr<const PeakDef>> background_auto_search_peaks;

  /** DRF (may be nullptr). */
  std::shared_ptr<const DetectorPeakResponse> drf;

  /** The coarse resolution type for this detector. */
  PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::Unknown;

  /** PeakFitDetPrefs with the correct m_det_type set. */
  std::shared_ptr<PeakFitDetPrefs> peak_fit_prefs;
};//struct PrecomputedNuclideData


/** Per-spectrum diagnostic detail for the background-fit penalty path.
 Populated by `compute_background_fit_penalty` when a non-null pointer is
 supplied; used by the HTML reporter to show which sources fit peaks in
 their backgrounds and at what significance. */
struct BackgroundFitDetail
{
  /** All peaks attributed to a source in the background fit.  Suppressed
   ones are recorded so the HTML can render them for sanity but they do
   not contribute to `penalty`. */
  std::vector<PeakDef> source_attributed_peaks;

  /** Per-peak normalized significance (post livetime normalization,
   pre-7.5 cap), parallel to `source_attributed_peaks`. */
  std::vector<double> normalized_significances;

  /** Per-peak suppression reason, parallel to `source_attributed_peaks`.
   Empty string means the peak contributed to `penalty`.  Non-empty
   values describe why the peak was excluded:
     "norm-like"       - the peak's source is in the NORM-like list
     "below 30 keV"    - low-energy filter (xrays / scatter)
     "near 511 keV"    - within 2 keV of annihilation energy
     "near NORM line"  - within ~1 FWHM of a strong NORM gamma */
  std::vector<std::string> suppression_reasons;

  /** Final per-spectrum penalty, post-cap.  Zero if all sources are
   NORM-like or background_auto_search_peaks is empty. */
  double penalty = 0.0;

  /** Set if the background fit failed; penalty is set to the per-spectrum
   cap in that case. */
  std::string error_message;
};//struct BackgroundFitDetail


/** The exact fit-quality components used by NuclideConfigEval and the GA objective.
 The objective is `area_cost - find_reward - candidate_reward
 + sm_miss_penalty_weight*miss_fraction`; the background term is recorded separately. */
struct FitAccuracyBreakdown
{
  double cost = 0.0;
  double area_cost = 0.0;
  double find_reward = 0.0;
  double candidate_reward = 0.0;
  double miss_fraction = 0.0;
  size_t missed_definitely_wanted = 0;
  size_t extra_peaks = 0;
};


/** One authoritative NuclideConfigEval result.  The HTML reporter consumes these cached records
 instead of fitting the spectrum a second time, so status, errors, peaks, and scores cannot drift
 from the command-line evaluation. */
struct SpectrumEvaluation
{
  std::string spectrum_id;
  std::string anchor_id;
  bool has_fit_result = false;
  bool exception = false;
  bool legitimate_empty = false;
  bool mechanical_failure = false;
  std::string status;
  std::string error_message;
  PeakFitResult fit_result;
  FitAccuracyBreakdown accuracy;
  BackgroundFitDetail background_detail;
  double background_penalty = 0.0;
};


struct ConfigEvaluation
{
  std::vector<SpectrumEvaluation> spectra;
  double total_fg = 0.0;
  double total_bg_raw = 0.0;
  double total_cost = 0.0;
  FitAccuracyBreakdown accuracy_totals;
  size_t successes = 0;
  size_t legitimate_empties = 0;
  size_t mechanical_failures = 0;
};


namespace ReportDetail
{
  /** Stable across thread counts and filesystem enumeration order for the same corpus paths. */
  std::string canonical_spectrum_key( const DataSrcInfo &info );
  std::string canonical_spectrum_key( const PrecomputedNuclideData &pd );
  std::string stable_spectrum_id( const PrecomputedNuclideData &pd );
  std::string html_escape( const std::string &text );
  size_t roi_channel_count( const SpecUtils::Measurement &measurement,
                            double lower_energy, double upper_energy );
  double roi_width_in_fwhm( double lower_energy, double upper_energy,
                            double representative_fwhm );
  void validate_evaluation_coverage( size_t selected_spectra,
                                     const ConfigEvaluation &evaluation );

  /** The shared objective helper used by both CLI scoring and the HTML reporter. */
  std::vector<ExpectedPhotopeakInfo> filter_truth_for_spectroscopic_extent(
    const std::vector<ExpectedPhotopeakInfo> &truth,
    const PrecomputedNuclideData &pd );
  FitAccuracyBreakdown score_observable_peaks(
    const PrecomputedNuclideData &pd,
    const std::vector<PeakDef> &observable_peaks );

  /** Runs inexpensive pure reporter checks.  Throws on failure. */
  void run_self_tests();
}//namespace ReportDetail


/** Export the same cached evaluation records and report-only ROI diagnostics consumed by the
 manual-review HTML.  Throws rather than overwriting an existing artifact. */
void write_config_evaluation_tsv(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const ConfigEvaluation &evaluation,
  const std::string &filename );


/** Computes the background-fit penalty for a single PrecomputedNuclideData.

 Runs `fit_peaks_for_nuclides` against the long_background spectrum (as the
 foreground argument, with no further background subtraction) using the
 given `config`, and sums livetime-normalized per-peak significances of
 source-attributed peaks for any source that is not NORM-like.  The result
 is capped at `sm_background_fit_penalty_per_spectrum_cap`.

 Returns 0.0 if all sources are NORM-like or `background_auto_search_peaks`
 is empty.  On fit failure or exception, returns the per-spectrum cap.

 If `detail_out` is non-null, populates it with the per-peak diagnostic
 detail for HTML rendering.
 */
double compute_background_fit_penalty(
    const PrecomputedNuclideData &pd,
    const PeakFitForNuclideConfig &config,
    BackgroundMode bg_mode,
    BackgroundFitDetail *detail_out );


/** Resolve a source name (from DataSrcInfo) to a vector of SrcVariant (nuclides, elements, reactions).
 Handles edge cases like Tl201, Pu238, U233, etc.
 Returns empty if the source should be skipped (e.g., Cf252, Am241Li).
 */
std::vector<RelActCalcAuto::SrcVariant> resolve_sources( const std::string &src_name );


/** Precompute the expensive data (search_for_peaks, source resolution, detector type) for all spectra.
 Skips sources that cannot be resolved (Cf252, Am241Li, etc.).
 */
std::vector<PrecomputedNuclideData> precompute_nuclide_data(
  const std::vector<DataSrcInfo> &srcs_info,
  const BackgroundMode bg_mode );


struct NuclideConfigSolution;  // defined below; needed for the `genes` parameter

/** Write the HTML and N42 output files for a given config, using the precomputed data.
 This is the shared logic extracted from eval_peaks_for_nuclide, and can be called
 from both the GA reporting and the PeaksForNuclide action.

 `genes` is the chromosome that produced `config`; it is rendered as a table in the report so the
 reader can see exactly which GA parameters generated the plots.
 */
void write_results_html_and_n42(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const PeakFitForNuclideConfig &config,
  const NuclideConfigSolution &genes,
  const BackgroundMode bg_mode,
  const std::string &html_filename,
  const std::string &n42_output_dir,
  const ConfigEvaluation *evaluation = nullptr,
  const std::string &run_metadata = std::string() );


/** The GA chromosome: maps ~50 genes to PeakFitForNuclideConfig fields. */
struct NuclideConfigSolution
{
  // --- FWHM functional form (0=kGadrasResolutionFcn, 1=kSqrtPolynomial, 2=kSqrtEnergyPlusInverse, 3=kConstantPlusSqrtEnergy) ---
  int fwhm_functional_form;

  // --- Manual RelEff parameters ---
  double rel_eff_manual_base_rel_eff_uncert;
  double initial_nuc_match_cluster_num_sigma;
  double manual_eff_cluster_num_sigma;

  // AICc complexity-penalty scales: manual rel-eff form/order and per-ROI continuum order are
  // selected per spectrum/per ROI at runtime; only these kappa scales are optimized.  Replaces
  // the former per-peak-count form/order genes and the width-threshold quadratic-continuum genes.
  double manual_releff_aicc_penalty;
  double cont_order_aicc_penalty;

  // Manual clustering keep-gate: min z = S_est/sqrt(S_est + B_est) vs the sideband continuum.
  // (Replaces the former absolute-count min_data_area/min_est_peak_area gates + S/sqrt(gross) sig.)
  double manual_keep_significance_z;
  double manual_rel_eff_sol_min_fwhm_roi;
  double manual_rel_eff_sol_max_fwhm;
  double manual_roi_core_num_fwhm;  // adaptive-extent core half-width, manual stage

  // --- Auto RelEff parameters ---
  int fwhm_form;  // RelActCalcAuto::FwhmForm (0..12, excluding NotApplicable)
  double rel_eff_auto_base_rel_eff_uncert;
  double auto_rel_eff_cluster_num_sigma;
  double auto_keep_significance_z;  // keep-gate z for the auto stage (see manual_keep_significance_z)
  double auto_roi_core_num_fwhm;    // adaptive-extent core half-width, auto stage
  double roi_extend_z;              // sideband block-consistency z (shared across stages)
  double roi_max_num_fwhm;          // extension cap beyond outermost gamma (shared)
  double auto_rel_eff_sol_max_fwhm;
  double auto_rel_eff_sol_min_fwhm_roi;
  int auto_roi_partition_overwide;  // 0/1: admit measured-data partition of over-wide components
  int auto_roi_final_fitted_partition;  // 0/1: admit a separately validated final fitted-ROI partition
  double auto_roi_partition_min_gap_fwhm;
  int auto_roi_partition_allow_clean_gap_override;
  double auto_roi_partition_residual_valley_max_excess_z;
  int auto_roi_final_partition_max_proposals;
  int auto_roi_final_partition_max_atoms;
  double auto_roi_final_partition_min_width_fwhm;
  int auto_roi_partition_max_children;
  double auto_roi_partition_force_gap_fwhm;

  // --- RelActAuto model parameters ---
  int rel_eff_eqn_type;   // RelActCalc::RelEffEqnForm (0..4, excluding FramPhysicalModel=4 is fine to include)
  int rel_eff_eqn_order;  // [0, 6]

  // --- Desperation physical model ---
  double desperation_phys_model_atomic_number;
  double desperation_phys_model_areal_density_g_per_cm2;

  // --- Physical model booleans ---
  int nucs_of_el_same_age;   // 0 or 1
  int phys_model_use_hoerl;  // 0 or 1
  int fit_energy_cal;         // 0 or 1

  // --- ROI merge/split: clean-gap anchoring test (shared across stages) ---
  double merge_tail_z;           // predicted-tail-vs-continuum-noise z per block
  double merge_clean_gap_fwhm;   // required clean-window width, in FWHM

  // --- ROI significance and observable peak thresholds ---
  // Single likelihood-ratio equivalent-z (replaces min_chi2_reduction / min_peak_sig /
  // min_quad_cont_chi2_dof, which were redundant parameterizations of one test).
  double roi_significance_z;
  double observable_peak_initial_significance_threshold;
  double observable_peak_final_significance_threshold;

  // --- Step continuum parameters ---
  double step_cont_min_peak_significance;
  double step_trial_chi2_margin;  // chi2 margin the step trial fit must win by (replaced step_cont_left_right_nsigma)

  // Note: skew_type is NOT optimized by GA - it is manually chosen per detector type
  // (e.g., NoSkew for most detectors, DoubleSidedCrystalBall for CZT)

  // --- Manual rel-eff chi2/dof cap (only varied when sm_rel_eff_chi2_cap_mode == GAOptimized) ---
  // Otherwise this field is overridden by resolve_chi2_dof_cap() at genes_to_settings time.
  double initial_manual_rel_eff_max_chi2_dof;

  std::string to_string( const std::string &separator ) const;

  /** Parse a line produced by `to_string()` back into a solution.  Name-keyed, so it is
   order-independent and tolerant of unknown keys; returns false if any known gene is missing or
   unparseable.  NOTE this means a checkpoint from before a gene RENAME does NOT resume (its old
   key is ignored as unknown, its new key is missing -> the individual is discarded, and if every
   line fails the resume exits) - fail-loud is intentional, since a checkpoint tuned against
   different gene semantics should not silently seed a new run.  E.g. checkpoints predating the
   step_cont_left_right_nsigma -> step_trial_chi2_margin replacement (2026-07) must be re-tuned
   from scratch.  Used by --resume. */
  static bool from_string( const std::string &line, const std::string &separator,
                           NuclideConfigSolution &out );
};//struct NuclideConfigSolution


struct NuclideConfigCost
{
  /** Combined optimization objective: `objective_fg + sm_background_fit_penalty_weight * objective_bg`.
   This is what openGA actually minimizes. */
  double objective1;

  /** Foreground-only contribution to the cost (summed across all spectra).
   Equals `objective1` when sm_do_background_fit_trial is false. */
  double objective_fg;

  /** Raw, un-weighted background-fit penalty summed across all spectra.
   Zero when sm_do_background_fit_trial is false. */
  double objective_bg;
};//struct NuclideConfigCost


typedef EA::Genetic<NuclideConfigSolution,NuclideConfigCost> GA_Type;
typedef EA::GenerationType<NuclideConfigSolution,NuclideConfigCost> Generation_Type;


/** Convert GA genes to a PeakFitForNuclideConfig. */
PeakFitForNuclideConfig genes_to_settings( const NuclideConfigSolution &solution );

/** Inverse of `genes_to_settings`: extract the GA-optimized fields of a config into a gene set.
 Round-trip exactness caveat: `genes_to_settings` resolves `initial_manual_rel_eff_max_chi2_dof`
 through `resolve_chi2_dof_cap()`, so the cap gene only round-trips when
 `sm_rel_eff_chi2_cap_mode == GAOptimized` (or the Fixed value happens to match).
 Used by NuclideConfigEval's `--config-genes default` to score/dump the production defaults. */
NuclideConfigSolution settings_to_genes( const PeakFitForNuclideConfig &config );


void init_genes( NuclideConfigSolution &p, const std::function<double(void)> &rnd01 );

bool eval_solution( const NuclideConfigSolution &p, NuclideConfigCost &c );

NuclideConfigSolution mutate( const NuclideConfigSolution &X_base,
                              const std::function<double(void)> &rnd01,
                              double shrink_scale );

NuclideConfigSolution crossover( const NuclideConfigSolution &X1,
                                 const NuclideConfigSolution &X2,
                                 const std::function<double(void)> &rnd01 );

double calculate_SO_total_fitness( const GA_Type::thisChromosomeType &X );

void SO_report_generation( int generation_number,
                           const EA::GenerationType<NuclideConfigSolution,NuclideConfigCost> &last_generation,
                           const NuclideConfigSolution &best_genes );


/** Run the GA optimization. Returns the best PeakFitForNuclideConfig found.

 @param precomputed The precomputed data for all spectra.
 @param ga_eval_fcn Function that takes a PeakFitForNuclideConfig and two
        optional out-pointers (foreground-only and raw-bg-penalty totals).
        Returns the combined cost (lower is better).  The out-pointers are
        used to surface the breakdown in reporting.  If `sm_do_background_fit_trial`
        is false, the bg component is always 0 and combined == foreground.
 */
PeakFitForNuclideConfig do_nuclide_config_ga(
  const std::vector<PrecomputedNuclideData> &precomputed,
  std::function<double( const PeakFitForNuclideConfig &, double *, double * )> ga_eval_fcn );

}//namespace NuclideConfig_GA

#endif //NuclideConfig_GA_h
