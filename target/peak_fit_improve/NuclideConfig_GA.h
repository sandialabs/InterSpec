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

/** Cap on the raw (pre-weight) per-spectrum background-fit penalty so a
 single spectrum's false-positive cluster (or fit failure) cannot dominate
 total_score.  30.0 is roughly the typical magnitude of a per-spectrum
 foreground score contribution; combined with the 0.25 weight, this bounds
 the per-spectrum penalty contribution at 7.5. */
constexpr double sm_background_fit_penalty_per_spectrum_cap = 30.0;

/** Penalty assigned to a single spectrum's foreground cost when the fit fails
 (non-Success status, or throws).  The GA MINIMIZES the objective, and a successful
 per-spectrum foreground cost is `total_weight - (find_weight + candidate_score)`,
 which ranges roughly from ~ -30 (all source peaks found with good areas) up to ~ +40
 (poor areas and/or many spurious peaks).  This penalty must be strictly worse (larger)
 than any plausible successful cost so the GA is driven away from configs that cannot
 produce a fit at all.
 TODO: normalize per-spectrum costs by expected-peak count so multi-line spectra do not
 dominate the summed objective; the failure penalty would then move to that normalized scale. */
constexpr double sm_fit_failure_penalty = 100.0;


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


/** Write the HTML and N42 output files for a given config, using the precomputed data.
 This is the shared logic extracted from eval_peaks_for_nuclide, and can be called
 from both the GA reporting and the PeaksForNuclide action.
 */
void write_results_html_and_n42(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const PeakFitForNuclideConfig &config,
  const BackgroundMode bg_mode,
  const std::string &html_filename,
  const std::string &n42_output_dir );


/** The GA chromosome: maps ~50 genes to PeakFitForNuclideConfig fields. */
struct NuclideConfigSolution
{
  // --- FWHM functional form (0=kGadrasResolutionFcn, 1=kSqrtPolynomial, 2=kSqrtEnergyPlusInverse, 3=kConstantPlusSqrtEnergy) ---
  int fwhm_functional_form;

  // --- Manual RelEff parameters ---
  double rel_eff_manual_base_rel_eff_uncert;
  double initial_nuc_match_cluster_num_sigma;
  double manual_eff_cluster_num_sigma;

  // Manual RelEff equation forms/orders per peak count
  int initial_manual_relEff_1peak_eqn_order;
  int initial_manual_relEff_1peak_form;
  int initial_manual_relEff_2peak_eqn_order;
  int initial_manual_relEff_2peak_form;
  int initial_manual_relEff_3peak_eqn_order;
  int initial_manual_relEff_3peak_form;
  int initial_manual_relEff_4peak_physical_use_hoerl;  // 0 or 1
  int initial_manual_relEff_4peak_eqn_order;
  int initial_manual_relEff_4peak_form;
  int initial_manual_relEff_many_peak_physical_use_hoerl;  // 0 or 1
  int initial_manual_relEff_many_peak_eqn_order;
  int initial_manual_relEff_manypeak_form;

  // Manual clustering thresholds
  double manual_rel_eff_sol_min_data_area_keep;
  double manual_rel_eff_sol_min_est_peak_area_keep;
  double manual_rel_eff_sol_min_est_significance_keep;
  double manual_rel_eff_sol_min_fwhm_roi;
  double manual_rel_eff_sol_min_fwhm_quad_cont;
  double manual_rel_eff_sol_max_fwhm;
  double manual_rel_eff_roi_width_num_fwhm_lower;
  double manual_rel_eff_roi_width_num_fwhm_upper;

  // --- Auto RelEff parameters ---
  int fwhm_form;  // RelActCalcAuto::FwhmForm (0..12, excluding NotApplicable)
  double rel_eff_auto_base_rel_eff_uncert;
  double auto_rel_eff_cluster_num_sigma;
  double auto_rel_eff_sol_min_data_area_keep;
  double auto_rel_eff_sol_min_est_peak_area_keep;
  double auto_rel_eff_sol_min_est_significance_keep;
  double auto_rel_eff_roi_width_num_fwhm_lower;
  double auto_rel_eff_roi_width_num_fwhm_upper;
  double auto_rel_eff_sol_max_fwhm;
  double auto_rel_eff_sol_min_fwhm_roi;
  double auto_rel_eff_sol_min_fwhm_quad_cont;

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

  // --- Tail contribution merge prevention ---
  double manual_rel_eff_min_tail_contribution;
  double manual_rel_eff_tail_width_scale_fwhm;
  double auto_rel_eff_min_tail_contribution;
  double auto_rel_eff_tail_width_scale_fwhm;

  // --- ROI significance and observable peak thresholds ---
  double roi_significance_min_chi2_reduction;
  double roi_significance_min_peak_sig;
  double roi_significance_min_quad_cont_chi2_dof;
  double observable_peak_initial_significance_threshold;
  double observable_peak_final_significance_threshold;

  // --- Step continuum parameters ---
  double step_cont_min_peak_area;
  double step_cont_min_peak_significance;
  double step_cont_left_right_nsigma;

  // Note: skew_type is NOT optimized by GA - it is manually chosen per detector type
  // (e.g., NoSkew for most detectors, DoubleSidedCrystalBall for CZT)

  // --- Manual rel-eff chi2/dof cap (only varied when sm_rel_eff_chi2_cap_mode == GAOptimized) ---
  // Otherwise this field is overridden by resolve_chi2_dof_cap() at genes_to_settings time.
  double initial_manual_rel_eff_max_chi2_dof;

  std::string to_string( const std::string &separator ) const;
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
