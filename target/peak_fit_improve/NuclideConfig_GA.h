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

  /** DRF (may be nullptr). */
  std::shared_ptr<const DetectorPeakResponse> drf;

  /** The coarse resolution type for this detector. */
  PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::Unknown;

  /** PeakFitDetPrefs with the correct m_det_type set. */
  std::shared_ptr<PeakFitDetPrefs> peak_fit_prefs;
};//struct PrecomputedNuclideData


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

  std::string to_string( const std::string &separator ) const;
};//struct NuclideConfigSolution


struct NuclideConfigCost
{
  double objective1;
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
 @param ga_eval_fcn Function that takes a PeakFitForNuclideConfig and returns a cost (lower is better).
 */
PeakFitForNuclideConfig do_nuclide_config_ga(
  const std::vector<PrecomputedNuclideData> &precomputed,
  std::function<double( const PeakFitForNuclideConfig & )> ga_eval_fcn );

}//namespace NuclideConfig_GA

#endif //NuclideConfig_GA_h
