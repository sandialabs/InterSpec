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

#ifndef FitPeaksForNuclides_h
#define FitPeaksForNuclides_h

#include <string>
#include <vector>
#include <memory>

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DetectorPeakResponse.h"

namespace SpecUtils
{
  class Measurement;
}

namespace FitPeaksForNuclides
{

// Settings for the gamma clustering algorithm - different values may be used
// for the initial RelActManual stage vs subsequent RelActAuto refinement stages
struct GammaClusteringSettings
{
  double cluster_num_sigma;         // How many sigma to use for clustering gamma lines
  double min_data_area_keep;        // Minimum data area to keep a cluster
  double min_est_peak_area_keep;    // Minimum estimated peak area to keep
  double min_est_significance_keep; // Minimum significance (est_counts / sqrt(data_area))
  double roi_width_num_fwhm_lower;  // FWHM below lowest gamma line for ROI extent
  double roi_width_num_fwhm_upper;  // FWHM above highest gamma line for ROI extent
  double max_fwhm_width;            // Maximum ROI width in FWHM before breaking up
  double min_fwhm_roi;              // Minimum ROI width in FWHM to keep
  double min_fwhm_quad_cont;        // Width threshold to use quadratic continuum

  // Parameters for synthetic spectrum-based ROI breaking
  // Region around minimum/maximum to compute significance (in FWHM units)
  double break_check_fwhm_fraction = 0.5;

  // Threshold for considering a peak "significant" between breakpoints (sigma)
  // Must have a peak exceeding this between any two breakpoints (or ROI edge and breakpoint)
  double break_peak_significance_threshold = 2.0;

  // When multiple breakpoint candidates have similar significance, use depth_score as tiebreaker
  // This defines "similar" - candidates within this sigma of each other are considered tied
  double break_significance_tie_threshold = 0.5;

  // Step continuum decision thresholds
  // Minimum peak area (counts) to consider checking for step continuum
  double step_cont_min_peak_area = 1000.0;
  // Minimum peak significance (peak_area / sqrt(data_area)) to consider step continuum
  double step_cont_min_peak_significance = 30.0;
  // How many sigma higher the left side must be than the right side to use step continuum
  double step_cont_left_right_nsigma = 3.0;
};


// Result from RelActAuto peak fitting
struct PeakFitResult
{
  RelActCalcAuto::RelActAutoSolution::Status status;
  std::string error_message;
  std::vector<std::string> warnings;  // warnings that don't prevent success

  // Peaks after combining overlapping peaks within ROIs.
  // Peaks that are close together (within 1.5 sigma) or where a smaller peak
  // is not statistically distinguishable from a larger peak's tail are merged
  // into a single peak with combined properties.
  std::vector<PeakDef> fit_peaks;

  // Original uncombined peaks from the fit - preserves all individual peak information.
  // This is the raw output from RelActAuto before peak combination.
  std::vector<PeakDef> uncombined_fit_peaks;

  // Peaks that are statistically observable in the spectrum.
  // Computed from fit_peaks by:
  // 1. Removing peaks with initial significance < threshold (using raw data area)
  // 2. Refitting each ROI with PeakFitLM::refitPeaksThatShareROI_LM (SmallRefinementOnly)
  // 3. Iteratively removing peaks with final significance < threshold and refitting
  std::vector<PeakDef> observable_peaks;

  RelActCalcAuto::RelActAutoSolution solution;  // only valid if status == Success
};//struct PeakFitResult


// Configuration parameters for peak fitting for nuclides
// These values will be optimized via genetic algorithm for different detector types
struct PeakFitForNuclideConfig
{
  // FWHM functional form
  DetectorPeakResponse::ResolutionFnctForm fwhm_functional_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;

  // RelActManual parameters for initial relative efficiency estimation
  double rel_eff_manual_base_rel_eff_uncert = 0.0;
  double initial_nuc_match_cluster_num_sigma = 1.5;
  double manual_eff_cluster_num_sigma = 1.5;

  // RelActManual equation form and order based on number of matched peaks
  size_t initial_manual_relEff_1peak_eqn_order = 0;
  RelActCalc::RelEffEqnForm initial_manual_relEff_1peak_form = RelActCalc::RelEffEqnForm::LnXLnY;

  size_t initial_manual_relEff_2peak_eqn_order = 0;
  RelActCalc::RelEffEqnForm initial_manual_relEff_2peak_form = RelActCalc::RelEffEqnForm::LnXLnY;

  size_t initial_manual_relEff_3peak_eqn_order = 1;
  RelActCalc::RelEffEqnForm initial_manual_relEff_3peak_form = RelActCalc::RelEffEqnForm::LnXLnY;

  bool initial_manual_relEff_4peak_physical_use_hoerl = false;
  size_t initial_manual_relEff_4peak_eqn_order = 2;
  RelActCalc::RelEffEqnForm initial_manual_relEff_4peak_form = RelActCalc::RelEffEqnForm::LnXLnY;

  bool initial_manual_relEff_many_peak_physical_use_hoerl = false;
  size_t initial_manual_relEff_many_peak_eqn_order = 3;
  RelActCalc::RelEffEqnForm initial_manual_relEff_manypeak_form = RelActCalc::RelEffEqnForm::LnXLnY;

  // ROI clustering thresholds for manual RelEff stage
  double manual_rel_eff_sol_min_data_area_keep = 10.0;
  double manual_rel_eff_sol_min_est_peak_area_keep = 5.0;
  double manual_rel_eff_sol_min_est_significance_keep = 2.0;
  double manual_rel_eff_sol_min_fwhm_roi = 1.0;
  double manual_rel_eff_sol_min_fwhm_quad_cont = 8.0;
  double manual_rel_eff_sol_max_fwhm = 15.0;
  double manual_rel_eff_roi_width_num_fwhm_lower = 3.0;
  double manual_rel_eff_roi_width_num_fwhm_upper = 3.0;

  // RelActAuto parameters
  RelActCalcAuto::FwhmForm fwhm_form = RelActCalcAuto::FwhmForm::Berstein_3;
  double num_sigma_contribution = 1.5;  // Peak area contribution sigma range
  double rel_eff_auto_base_rel_eff_uncert = 0.1;  // BR uncertainty for RelActAuto

  // ROI clustering thresholds for auto RelEff refinement stage
  double auto_rel_eff_cluster_num_sigma = 2.0;  // Slightly wider clustering with better rel-eff
  double auto_rel_eff_sol_min_data_area_keep = 10.0;
  double auto_rel_eff_sol_min_est_peak_area_keep = 5.0;
  double auto_rel_eff_sol_min_est_significance_keep = 3.0;
  double auto_rel_eff_roi_width_num_fwhm_lower = 3.5;  // Slightly more generous for refined fit
  double auto_rel_eff_roi_width_num_fwhm_upper = 3.5;
  double auto_rel_eff_sol_max_fwhm = 12.0;  // Tighter constraint as solution improves
  double auto_rel_eff_sol_min_fwhm_roi = 1.0;
  double auto_rel_eff_sol_min_fwhm_quad_cont = 8.0;


  // RelActAuto relative efficiency model parameters
  RelActCalc::RelEffEqnForm rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnXLnY;
  size_t rel_eff_eqn_order = 2;

  // Physical model shielding (empty vectors mean no shielding)
  std::vector<std::shared_ptr<RelActCalc::PhysicalModelShieldInput>> phys_model_self_atten;
  std::vector<std::shared_ptr<RelActCalc::PhysicalModelShieldInput>> phys_model_external_atten;

  // Desperation Physical Model shielding configuration
  // Atomic number for desperation external shielding (0 = don't use, 1-98 = valid elements)
  double desperation_phys_model_atomic_number = 26.0;  // Default: iron

  // Areal density starting value for desperation external shielding (in g/cm2)
  double desperation_phys_model_areal_density_g_per_cm2 = 1.0;  // Default: 1.0 g/cm2

  // Physical model options (only used when rel_eff_eqn_type == FramPhysicalModel)
  bool nucs_of_el_same_age = true;
  bool phys_model_use_hoerl = true;

  // Fields for RelActAuto options configuration
  RelActCalcAuto::EnergyCalFitType energy_cal_type = RelActCalcAuto::EnergyCalFitType::NonLinearFit;
  
  // ROI significance threshold for iterative refinement
  // Minimum total chi2 reduction required for peaks in a ROI to be considered significant
  // The chi2 with peaks must be at least this much lower than chi2 with continuum-only
  double roi_significance_min_chi2_reduction = 10.0;

  // Minimum peak significance (peak_area / sqrt(continuum)) for a peak to be considered significant
  // If any peak in a ROI has significance above this threshold, the ROI is considered significant
  // (alternative to chi2 reduction test - ROI passes if EITHER test passes)
  double roi_significance_min_peak_sig = 3.5;

  // Threshold for initial peak significance filter before refitting for observable_peaks.
  // Significance = (peak_amplitude * 0.7607) / sqrt(data_area_in_pm_1_fwhm)
  // 0.7607 is fraction of Gaussian area within +/-1 FWHM
  double observable_peak_initial_significance_threshold = 2.25;

  // Threshold for final peak significance after refitting for observable_peaks.
  // Significance = peak_area / peak_area_uncertainty
  double observable_peak_final_significance_threshold = 2.0;

  // Step continuum decision parameters
  // Minimum peak area (counts) to consider checking for step continuum
  double step_cont_min_peak_area = 1000.0;
  // Minimum peak significance (peak_area / sqrt(data_area)) to consider step continuum
  double step_cont_min_peak_significance = 30.0;
  // How many sigma higher the left side must be than the right side to use step continuum
  double step_cont_left_right_nsigma = 3.0;



  // Get GammaClusteringSettings for manual RelEff stage
  GammaClusteringSettings get_manual_clustering_settings() const
  {
    GammaClusteringSettings settings;
    settings.cluster_num_sigma = manual_eff_cluster_num_sigma;
    settings.min_data_area_keep = manual_rel_eff_sol_min_data_area_keep;
    settings.min_est_peak_area_keep = manual_rel_eff_sol_min_est_peak_area_keep;
    settings.min_est_significance_keep = manual_rel_eff_sol_min_est_significance_keep;
    settings.roi_width_num_fwhm_lower = manual_rel_eff_roi_width_num_fwhm_lower;
    settings.roi_width_num_fwhm_upper = manual_rel_eff_roi_width_num_fwhm_upper;
    settings.max_fwhm_width = manual_rel_eff_sol_max_fwhm;
    settings.min_fwhm_roi = manual_rel_eff_sol_min_fwhm_roi;
    settings.min_fwhm_quad_cont = manual_rel_eff_sol_min_fwhm_quad_cont;
    settings.step_cont_min_peak_area = step_cont_min_peak_area;
    settings.step_cont_min_peak_significance = step_cont_min_peak_significance;
    settings.step_cont_left_right_nsigma = step_cont_left_right_nsigma;
    return settings;
  }

  // Get GammaClusteringSettings for auto refinement stage
  GammaClusteringSettings get_auto_clustering_settings() const
  {
    GammaClusteringSettings settings;
    settings.cluster_num_sigma = auto_rel_eff_cluster_num_sigma;
    settings.min_data_area_keep = auto_rel_eff_sol_min_data_area_keep;
    settings.min_est_peak_area_keep = auto_rel_eff_sol_min_est_peak_area_keep;
    settings.min_est_significance_keep = auto_rel_eff_sol_min_est_significance_keep;
    settings.roi_width_num_fwhm_lower = auto_rel_eff_roi_width_num_fwhm_lower;
    settings.roi_width_num_fwhm_upper = auto_rel_eff_roi_width_num_fwhm_upper;
    settings.max_fwhm_width = auto_rel_eff_sol_max_fwhm;
    settings.min_fwhm_roi = auto_rel_eff_sol_min_fwhm_roi;
    settings.min_fwhm_quad_cont = auto_rel_eff_sol_min_fwhm_quad_cont;
    settings.step_cont_min_peak_area = step_cont_min_peak_area;
    settings.step_cont_min_peak_significance = step_cont_min_peak_significance;
    settings.step_cont_left_right_nsigma = step_cont_left_right_nsigma;
    return settings;
  }
};//struct PeakFitForNuclideConfig


/** Comprehensive peak fitting function for nuclides.

 This function performs the complete peak fitting workflow:
 1. Determines FWHM functional form from auto-search peaks or DRF
 2. Matches auto-search peaks to source nuclides
 3. Performs RelActManual fit to get initial relative efficiency estimate
 4. Clusters gamma lines into ROIs based on manual rel-eff
 5. Calls fit_peaks_for_nuclide_relactauto with three different configurations
 6. Returns the best result based on chi2

 @param auto_search_peaks Initial peaks found by automatic search
 @param foreground Foreground spectrum to fit
 @param sources Vector of source nuclides to fit
 @param long_background Long background spectrum (can be nullptr)
 @param drf_input Detector response function (can be nullptr, will use generic if needed)
 @param config Configuration for peak fitting parameters
 @param isHPGe Whether this is an HPGe detector
 @return PeakFitResult with status, error message, fit peaks, and solution
 */
  
PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::shared_ptr<const SpecUtils::Measurement> &long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const PeakFitForNuclideConfig &config,
  const bool isHPGe );
  
PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const SpecUtils::Measurement> &long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const PeakFitForNuclideConfig &config,
  const bool isHPGe );
  

// Helper function for estimating initial ROIs when no peaks are available
std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_without_peaks(
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const bool isHPGe,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const double min_valid_energy,
  const double max_valid_energy,
  const PeakFitForNuclideConfig &config );

}//namespace FitPeaksForNuclides

#endif //FitPeaksForNuclides_h
