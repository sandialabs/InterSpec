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

#include <chrono>
#include <random>
#include <string>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp> // For boost::asio::post

#include <Wt/WColor>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "FinalFit_GA.h"
#include "PeakFitImprove.h"
#include "InitialFit_GA.h"
#include "CandidatePeak_GA.h"

using namespace std;
namespace FitPeaksForNuclideDev
{

// Struct to store cluster info including gamma line energies and amplitudes
struct ClusteredGammaInfo {
  double lower;
  double upper;
  vector<double> gamma_energies;    // energies of gamma lines in this cluster
  vector<double> gamma_amplitudes;  // expected peak areas/amplitudes
};


// Settings for the gamma clustering algorithm - different values may be used
// for the initial RelActManual stage vs subsequent RelActAuto refinement stages
struct GammaClusteringSettings {
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


namespace
{
  /** Get photon energies and intensities for a source at a given age.

   For nuclides: Uses NuclideMixture with the specified age
   For elements: Returns xrays (age parameter ignored)
   For reactions: Returns gammas (age parameter ignored)

   \param src The source variant
   \param activity Activity in becquerels (for nuclides) or scaling factor (for elements/reactions)
   \param age Age in PhysicalUnits seconds (for nuclides only, ignored for others)
   \return Vector of EnergyRatePair with photon energies and rates
   \throws runtime_error if source is null or invalid
   */
  std::vector<SandiaDecay::EnergyRatePair> get_source_photons(
    const RelActCalcAuto::SrcVariant &src,
    const double activity,
    const double age )
  {
    std::vector<SandiaDecay::EnergyRatePair> result;

    const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( src );
    if( nuc )
    {
      // Use NuclideMixture for nuclides
      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( nuc, activity, age );
      result = mix.photons( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy );
      return result;
    }

    const SandiaDecay::Element *el = RelActCalcAuto::element( src );
    if( el )
    {
      // Convert xrays to EnergyRatePair format
      // Xrays have relative intensities, scale by activity for consistency
      for( const SandiaDecay::EnergyIntensityPair &xray : el->xrays )
      {
        // intensity is relative to strongest line (1.0 = strongest)
        // Scale by activity to get rate
        result.emplace_back( activity * xray.intensity, xray.energy );
      }

      return result;
    }

    const ReactionGamma::Reaction *rxn = RelActCalcAuto::reaction( src );
    if( rxn )
    {
      // Convert reaction gammas to EnergyRatePair format
      for( const ReactionGamma::Reaction::EnergyYield &gamma : rxn->gammas )
      {
        // abundance is yield per reaction
        // Scale by activity to get rate
        result.emplace_back( activity * gamma.abundance, gamma.energy );
      }

      return result;
    }

    throw runtime_error( "get_source_photons: invalid or null source" );
  }//get_source_photons(...)


  /** Get appropriate age for a source.

   For nuclides: Returns PeakDef::defaultDecayTime() if age < 0, else returns age
   For elements/reactions: Returns 0.0 (age not applicable)

   \param src The source variant
   \param age_input User-specified age (negative means use default for nuclides)
   \return Age in PhysicalUnits seconds
   */
  double get_source_age( const RelActCalcAuto::SrcVariant &src, const double age_input )
  {
    const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( src );
    if( nuc )
    {
      if( age_input >= 0.0 )
        return age_input;

      return PeakDef::defaultDecayTime( nuc );
    }

    // Elements and reactions don't have age
    return 0.0;
  }//get_source_age(...)


  /** Determines if external shielding should be used for desperation Physical Model attempts.

   Shielding is disabled if:
   - Atomic number is invalid (0, < 1, or > 98)
   - There are fewer than 2 ROIs
   - Energy range is too small (< 60 keV if any ROI extends below 120 keV, otherwise < 100 keV)

   @param atomic_number The configured atomic number for desperation shielding
   @param rois The vector of ROI ranges to check
   @return true if shielding should be used, false otherwise
   */
  bool should_use_desperation_shielding( const double atomic_number,
                                         const std::vector<RelActCalcAuto::RoiRange> &rois )
  {
    // Check 1: Valid atomic number
    if( atomic_number < 1.0 || atomic_number > 98.0 )
      return false;

    // Check 2: At least 2 ROIs required
    if( rois.size() < 2 )
      return false;

    // Check 3: Sufficient energy range
    // Find min lower_energy and max upper_energy across all ROIs
    double min_lower = std::numeric_limits<double>::max();
    double max_upper = std::numeric_limits<double>::lowest();
    bool has_low_energy_roi = false;

    for( const RelActCalcAuto::RoiRange &roi : rois )
    {
      min_lower = std::min( min_lower, roi.lower_energy );
      max_upper = std::max( max_upper, roi.upper_energy );

      if( roi.lower_energy < 120.0 )
        has_low_energy_roi = true;
    }

    const double energy_range = max_upper - min_lower;
    const double required_range = has_low_energy_roi ? 60.0 : 100.0;

    if( energy_range < required_range )
      return false;

    return true;
  }//should_use_desperation_shielding(...)


  /** Creates a PhysicalModelShieldInput for desperation attempts with the specified atomic number.

   @param atomic_number The atomic number to use (must be in range 1-98)
   @param starting_areal_density The starting areal density in g/cm2
   @return Shared pointer to configured PhysicalModelShieldInput
   @throws std::invalid_argument if atomic_number is out of valid range
   */
  std::shared_ptr<RelActCalc::PhysicalModelShieldInput>
  create_desperation_shielding( const double atomic_number, const double starting_areal_density )
  {
    if( atomic_number < 1.0 || atomic_number > 98.0 )
      throw std::invalid_argument( "create_desperation_shielding: atomic_number must be in range [1, 98]" );

    std::shared_ptr<RelActCalc::PhysicalModelShieldInput> shield
      = std::make_shared<RelActCalc::PhysicalModelShieldInput>();

    // Configure as generic element (not material)
    shield->atomic_number = atomic_number;
    shield->material = nullptr;

    // Set starting areal density (convert from g/cm2 to PhysicalUnits)
    shield->areal_density = starting_areal_density * PhysicalUnits::g_per_cm2;

    // Do NOT fit atomic number - only fit areal density
    shield->fit_atomic_number = false;

    // Configure areal density fitting with reasonable bounds
    shield->fit_areal_density = true;
    shield->lower_fit_areal_density = 0.1 * PhysicalUnits::g_per_cm2;
    shield->upper_fit_areal_density = 10.0 * PhysicalUnits::g_per_cm2;

    return shield;
  }//create_desperation_shielding(...)
}//namespace


// Forward declaration
vector<RelActCalcAuto::RoiRange> cluster_gammas_to_rois(
  const function<double(double)> &rel_eff_fcn,
  const vector<pair<RelActCalcAuto::SrcVariant, double>> &sources_and_activities,
  const shared_ptr<const SpecUtils::Measurement> &foreground,
  const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
  const vector<float> &fwhm_coefficients,
  const double fwhm_lower_energy,
  const double fwhm_upper_energy,
  const double lowest_energy,
  const double highest_energy,
  const GammaClusteringSettings &settings );


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


/** Combined scoring struct that evaluates both final fit quality and peak finding completeness.

 This struct provides a comprehensive assessment by combining:
 - FinalFitScore: Area, width, and position accuracy metrics
 - PeakFindAndFitWeights: Peak detection completeness and area accuracy
 - final_weight: Combined score (simple sum of find_weight + total_weight)
 */
struct CombinedPeakFitScore
{
  FinalFit_GA::FinalFitScore final_fit_score;
  InitialFit_GA::PeakFindAndFitWeights initial_fit_weights;
  CandidatePeak_GA::CandidatePeakScore candidate_peak_score;

  /** Combined weight score (find_weight + total_weight).
   * find_weight: higher is better (more correct peaks found)
   * total_weight: lower is better (more accurate areas)
   */
  double final_weight = 0.0;

  std::string print( const std::string &varname ) const
  {
    std::string answer;

    answer += "=== " + varname + " - Final Fit Metrics ===\n";
    answer += final_fit_score.print( varname + ".final_fit_score" );

    answer += "\n=== " + varname + " - Peak Finding Metrics ===\n";
    answer += varname + ".initial_fit_weights.find_weight =                     " + std::to_string(initial_fit_weights.find_weight) + "\n";
    answer += varname + ".initial_fit_weights.def_wanted_area_weight =          " + std::to_string(initial_fit_weights.def_wanted_area_weight) + "\n";
    answer += varname + ".initial_fit_weights.maybe_wanted_area_weight =        " + std::to_string(initial_fit_weights.maybe_wanted_area_weight) + "\n";
    answer += varname + ".initial_fit_weights.not_wanted_area_weight =          " + std::to_string(initial_fit_weights.not_wanted_area_weight) + "\n";
    answer += varname + ".initial_fit_weights.def_wanted_area_median_weight =   " + std::to_string(initial_fit_weights.def_wanted_area_median_weight) + "\n";
    answer += varname + ".initial_fit_weights.maybe_wanted_area_median_weight = " + std::to_string(initial_fit_weights.maybe_wanted_area_median_weight) + "\n";
    answer += varname + ".initial_fit_weights.not_wanted_area_median_weight =   " + std::to_string(initial_fit_weights.not_wanted_area_median_weight) + "\n";

    answer += "\n=== " + varname + " - Candidate Peak Metrics ===\n";
    answer += varname + ".candidate_peak_score.score =                           " + std::to_string(candidate_peak_score.score) + "\n";
    answer += varname + ".candidate_peak_score.num_peaks_found =                 " + std::to_string(candidate_peak_score.num_peaks_found) + "\n";
    answer += varname + ".candidate_peak_score.num_def_wanted_not_found =         " + std::to_string(candidate_peak_score.num_def_wanted_not_found) + "\n";
    answer += varname + ".candidate_peak_score.num_def_wanted_peaks_found =      " + std::to_string(candidate_peak_score.num_def_wanted_peaks_found) + "\n";
    answer += varname + ".candidate_peak_score.num_possibly_accepted_peaks_not_found = " + std::to_string(candidate_peak_score.num_possibly_accepted_peaks_not_found) + "\n";
    answer += varname + ".candidate_peak_score.num_extra_peaks =                " + std::to_string(candidate_peak_score.num_extra_peaks) + "\n";

    answer += "\n=== " + varname + " - Combined Metrics ===\n";
    answer += varname + ".final_weight = " + std::to_string(final_weight) + "\n";

    return answer;
  }//print(...)
};//struct CombinedPeakFitScore


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
  bool fit_energy_cal = true;
  std::vector<RelActCalcAuto::NucInputInfo> base_nuclides;  // Set from outside

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


/** Determines if two peaks should be combined based on proximity and relative contribution.

Peaks are combined if:
1. They are within `always_combine_nsigma` sigma of each other (using larger peak's sigma), OR
2. The smaller peak is not statistically distinguishable from the larger peak's tail
   (smaller peak area < 4 * sqrt(larger peak contribution at smaller peak's location))

@param larger_peak The peak with larger area (must have >= area of smaller_peak)
@param smaller_peak The peak with smaller area
@param always_combine_nsigma Peaks within this many sigma are always combined (default 1.5)
@return true if peaks should be combined, false otherwise
*/
bool should_combine_peaks( const PeakDef &larger_peak,
                           const PeakDef &smaller_peak,
                           const double always_combine_nsigma = 1.5 )
{
  assert( larger_peak.peakArea() >= smaller_peak.peakArea() );

  const double dist = std::fabs( larger_peak.mean() - smaller_peak.mean() );
  const double sigma_large = larger_peak.sigma();

  // Criterion 1: Always combine if within threshold sigma
  if( dist < always_combine_nsigma * sigma_large )
    return true;

  // Criterion 2: Check if smaller peak is overwhelmed by larger peak's contribution
  // Calculate larger peak's contribution at smaller peak's location
  const double sigma_small = smaller_peak.sigma();
  const double mean_small = smaller_peak.mean();

  // Integrate larger peak over +/- 0.5 sigma of smaller peak's mean
  const double x0 = mean_small - 0.5 * sigma_small;
  const double x1 = mean_small + 0.5 * sigma_small;
  const double contribution = larger_peak.gauss_integral( x0, x1 );

  // If smaller peak area < 4 * sqrt(contribution), combine
  // This means the smaller peak is not statistically distinguishable from the tail of the larger peak
  const double smaller_area = smaller_peak.peakArea();
  if( contribution > 0.0 && smaller_area < 4.0 * std::sqrt( contribution ) )
    return true;

  return false;
}//should_combine_peaks


/** Combines multiple peaks into a single peak using area-weighted averaging.

Creates a new peak with:
- Mean: area-weighted average of input peak means
- Sigma: area-weighted average of input peak sigmas
- Area: sum of input peak areas
- Uncertainty: quadrature sum of input peak uncertainties
- Source assignment: from the peak with largest area (dominant peak)
- Continuum: shared from input peaks (must all share same continuum)

@param peaks_to_combine Vector of peaks to combine (must share same continuum, must not be empty)
@return Combined peak
@throws std::invalid_argument if peaks vector is empty or peaks dont share continuum
*/
PeakDef combine_peaks( const std::vector<const PeakDef *> &peaks_to_combine )
{
  assert( !peaks_to_combine.empty() );
  if( peaks_to_combine.empty() )
    throw std::invalid_argument( "combine_peaks: empty peaks vector" );

  // Find dominant peak (largest area) and verify all share same continuum
  const PeakDef *dominant = peaks_to_combine[0];
  const std::shared_ptr<const PeakContinuum> cont = peaks_to_combine[0]->continuum();

  double total_area = 0.0;
  double sum_area_mean = 0.0;
  double sum_area_sigma = 0.0;
  double sum_uncert_sq = 0.0;

  for( const PeakDef *peak : peaks_to_combine )
  {
    assert( peak->continuum() == cont );
    if( peak->continuum() != cont )
      throw std::invalid_argument( "combine_peaks: peaks must share same continuum" );

    const double area = peak->peakArea();
    const double uncert = peak->peakAreaUncert();

    total_area += area;
    sum_area_mean += area * peak->mean();
    sum_area_sigma += area * peak->sigma();
    sum_uncert_sq += uncert * uncert;

    if( area > dominant->peakArea() )
      dominant = peak;
  }//for( const PeakDef *peak : peaks_to_combine )

  // Start with a copy of the dominant peak - this copies the continuum, source assignment,
  // skew type, line color, and other settings
  PeakDef combined = *dominant;

  // Update the gaussian parameters to the combined values
  combined.setMean( sum_area_mean / total_area );
  combined.setSigma( sum_area_sigma / total_area );
  combined.setPeakArea( total_area );
  combined.setPeakAreaUncert( std::sqrt( sum_uncert_sq ) );

  return combined;
}//combine_peaks


/** Combines overlapping peaks within each ROI based on proximity and contribution criteria.

For each ROI (group of peaks sharing a continuum):
1. Identifies peaks that should be combined using proximity and contribution checks
2. Uses transitive closure (Union-Find) to group connected peaks
3. Creates combined peaks with area-weighted properties

@param uncombined_peaks The original fit peaks, potentially with multiple peaks per ROI
@return Vector of peaks with overlapping peaks combined within each ROI
*/
std::vector<PeakDef> combine_overlapping_peaks_in_rois(
    const std::vector<PeakDef> &uncombined_peaks )
{
  if( uncombined_peaks.empty() )
    return {};

  // Phase 1: Group peaks by continuum (ROI)
  std::map<std::shared_ptr<const PeakContinuum>, std::vector<size_t>> continuum_to_peak_indices;
  for( size_t i = 0; i < uncombined_peaks.size(); ++i )
  {
    const std::shared_ptr<const PeakContinuum> cont = uncombined_peaks[i].continuum();
    continuum_to_peak_indices[cont].push_back( i );
  }

  std::vector<PeakDef> combined_peaks;
  combined_peaks.reserve( uncombined_peaks.size() );

  if( PeakFitImprove::debug_printout )
  {
    cout << "combine_overlapping_peaks_in_rois: " << continuum_to_peak_indices.size()
         << " unique ROIs from " << uncombined_peaks.size() << " peaks" << endl;
  }

  // Phase 2 & 3: Process each ROI independently
  for( const auto &roi_entry : continuum_to_peak_indices )
  {
    const std::vector<size_t> &peak_indices = roi_entry.second;

    if( peak_indices.size() == 1 )
    {
      // Single peak in ROI - no combination needed
      combined_peaks.push_back( uncombined_peaks[peak_indices[0]] );
      continue;
    }

    if( PeakFitImprove::debug_printout )
    {
      cout << "  ROI with " << peak_indices.size() << " peaks at energies:";
      for( size_t idx : peak_indices )
        cout << " " << uncombined_peaks[idx].mean();
      cout << endl;
    }

    // Sort indices by peak area (descending) - we'll process largest peaks first
    std::vector<size_t> sorted_indices = peak_indices;
    std::sort( sorted_indices.begin(), sorted_indices.end(),
      [&uncombined_peaks]( size_t a, size_t b ) {
        return uncombined_peaks[a].peakArea() > uncombined_peaks[b].peakArea();
      });

    // Greedy clustering: start with largest peak, cluster nearby peaks into it,
    // then move to next largest unclustered peak
    std::set<size_t> clustered;  // Tracks which indices have been assigned to a cluster
    std::vector<std::vector<const PeakDef *>> clusters;

    for( size_t i = 0; i < sorted_indices.size(); ++i )
    {
      const size_t idx_i = sorted_indices[i];

      // Skip if this peak has already been clustered with a larger peak
      if( clustered.count( idx_i ) )
        continue;

      // Start a new cluster with this peak (the largest remaining unclustered peak)
      const PeakDef &anchor_peak = uncombined_peaks[idx_i];
      std::vector<const PeakDef *> cluster;
      cluster.push_back( &anchor_peak );
      clustered.insert( idx_i );

      // Find all smaller peaks that should be combined with this anchor peak
      for( size_t j = i + 1; j < sorted_indices.size(); ++j )
      {
        const size_t idx_j = sorted_indices[j];

        // Skip if already in another cluster
        if( clustered.count( idx_j ) )
          continue;

        const PeakDef &candidate_peak = uncombined_peaks[idx_j];

        // Check if this smaller peak should be combined with the anchor peak
        if( should_combine_peaks( anchor_peak, candidate_peak, 1.5 ) )
        {
          cluster.push_back( &candidate_peak );
          clustered.insert( idx_j );
        }
      }//for( size_t j = i + 1; j < sorted_indices.size(); ++j )

      clusters.push_back( std::move( cluster ) );
    }//for( size_t i = 0; i < sorted_indices.size(); ++i )

    if( PeakFitImprove::debug_printout )
    {
      cout << "    -> " << clusters.size() << " clusters after greedy clustering" << endl;
    }

    // Create combined peaks for each cluster
    std::vector<PeakDef> roi_peaks;  // Peaks for this ROI
    for( const std::vector<const PeakDef *> &cluster : clusters )
    {
      if( cluster.size() == 1 )
      {
        roi_peaks.push_back( *cluster[0] );
      }
      else
      {
        if( PeakFitImprove::debug_printout )
        {
          cout << "    Combining " << cluster.size() << " peaks at energies:";
          for( const PeakDef *p : cluster )
            cout << " " << p->mean();
          cout << endl;
        }
        roi_peaks.push_back( combine_peaks( cluster ) );
      }
    }//for( const std::vector<const PeakDef *> &cluster : clusters )

    // Make all peaks in this ROI share a new continuum
    if( !roi_peaks.empty() )
    {
      auto roi_continuum = make_shared<PeakContinuum>( *roi_peaks.front().continuum() );
      for( PeakDef &peak : roi_peaks )
        peak.setContinuum( roi_continuum );
    }

    // Add this ROI's peaks to the combined output
    combined_peaks.insert( combined_peaks.end(), roi_peaks.begin(), roi_peaks.end() );
  }//for( const auto &roi_entry : continuum_to_peak_indices )

  // Sort combined peaks by mean energy
  std::sort( combined_peaks.begin(), combined_peaks.end(),
    []( const PeakDef &lhs, const PeakDef &rhs ) {
      return lhs.mean() < rhs.mean();
    });

  if( PeakFitImprove::debug_printout )
  {
    cout << "combine_overlapping_peaks_in_rois: " << uncombined_peaks.size()
         << " input peaks -> " << combined_peaks.size() << " output peaks" << endl;
  }

  return combined_peaks;
}//combine_overlapping_peaks_in_rois


/** Computes observable peaks from fit_peaks by filtering for significance and refitting.

Algorithm:
1. Initial filter: Remove peaks where (amplitude * 0.7607) / sqrt(data_area_pm_1_fwhm) < initial_threshold
   - 0.7607 = fraction of Gaussian within +/-1 FWHM
   - data_area = total spectrum counts in [mean - fwhm, mean + fwhm]
2. Adjust ROI bounds if edge peaks were removed, create new continuum for affected ROIs
3. Group remaining peaks by shared PeakContinuum (ROI)
4. For each ROI (in parallel):
   a. Refit using PeakFitLM::refitPeaksThatShareROI_LM with SmallRefinementOnly option
   b. Check each peak: significance = area / area_uncertainty
   c. Remove peaks with significance < final_threshold
   d. If edge peaks removed, adjust ROI bounds and create new continuum
   e. If any peaks removed, repeat from (a) until all remaining peaks are significant
5. Return all peaks that passed filtering

@param fit_peaks Input peaks from RelActAuto fit
@param foreground The spectrum measurement
@param config Configuration containing thresholds and ROI width settings
@return Vector of observable peaks
*/
std::vector<PeakDef> compute_observable_peaks(
  const std::vector<PeakDef> &fit_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const PeakFitForNuclideConfig &config )
{
  
#if( PERFORM_DEVELOPER_CHECKS )
  // Debug: Check if input fit_peaks have overlapping ROIs
  if( PeakFitImprove::debug_printout )
  {
    std::map<std::shared_ptr<const PeakContinuum>, std::vector<double>> continuum_to_peak_means_input;
    for( const PeakDef &peak : fit_peaks )
      continuum_to_peak_means_input[peak.continuum()].push_back( peak.mean() );

    struct RoiDebugInfoInput
    {
      double lower_energy;
      double upper_energy;
      std::vector<double> peak_means;
    };

    std::vector<RoiDebugInfoInput> rois_debug_input;
    for( const auto &entry : continuum_to_peak_means_input )
    {
      RoiDebugInfoInput info;
      info.lower_energy = entry.first->lowerEnergy();
      info.upper_energy = entry.first->upperEnergy();
      info.peak_means = entry.second;
      rois_debug_input.push_back( info );
    }

    std::sort( begin(rois_debug_input), end(rois_debug_input),
      []( const RoiDebugInfoInput &a, const RoiDebugInfoInput &b ) { return a.lower_energy < b.lower_energy; } );

    bool found_overlap_input = false;
    for( size_t i = 1; i < rois_debug_input.size(); ++i )
    {
      const RoiDebugInfoInput &prev_roi = rois_debug_input[i - 1];
      const RoiDebugInfoInput &curr_roi = rois_debug_input[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap_input )
          cerr << "compute_observable_peaks: INPUT fit_peaks ALREADY HAVE OVERLAPPING ROIs:" << endl;
        found_overlap_input = true;

        cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : prev_roi.peak_means )
          cerr << mean << " ";
        cerr << "keV" << endl;

        cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : curr_roi.peak_means )
          cerr << mean << " ";
        cerr << "keV" << endl;
        cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << endl;
      }
    }

    if( !found_overlap_input && PeakFitImprove::debug_printout )
      cout << "compute_observable_peaks: Input fit_peaks have no overlapping ROIs" << endl;
  }//if( PeakFitImprove::debug_printout )
#endif

  // Fraction of Gaussian area within +/-1 FWHM
  // erf(sqrt(ln(2))) = 0.7607
  const double fwhm_fraction = 0.7607;
  const double initial_significance_threshold = config.observable_peak_initial_significance_threshold;
  const double final_significance_threshold = config.observable_peak_final_significance_threshold;
  const double roi_width_num_fwhm_lower = config.auto_rel_eff_roi_width_num_fwhm_lower;
  const double roi_width_num_fwhm_upper = config.auto_rel_eff_roi_width_num_fwhm_upper;

  // Lambda to adjust ROI bounds when edge peaks are removed.
  // Returns true if bounds were adjusted, false otherwise.
  // Takes peaks sorted by mean energy.
  const auto reduce_roi_bounds_if_needed = [roi_width_num_fwhm_lower, roi_width_num_fwhm_upper](
    std::vector<PeakDef> &peaks,
    const double orig_left_mean,
    const double orig_right_mean ) -> bool
  {
    if( peaks.empty() )
      return false;

    // Sort peaks by mean for edge detection
    std::sort( begin(peaks), end(peaks), &PeakDef::lessThanByMean );

    const double new_left_mean = peaks.front().mean();
    const double new_right_mean = peaks.back().mean();
    const bool left_edge_changed = (fabs(new_left_mean - orig_left_mean) > 0.1);
    const bool right_edge_changed = (fabs(new_right_mean - orig_right_mean) > 0.1);

    if( !left_edge_changed && !right_edge_changed )
      return false;

    // Get current continuum from first peak
    std::shared_ptr<const PeakContinuum> old_continuum = peaks.front().continuum();

    // Calculate new ROI bounds based on new edge peaks
    const double new_left_fwhm = peaks.front().fwhm();
    const double new_right_fwhm = peaks.back().fwhm();
    double new_lower_energy = new_left_mean - roi_width_num_fwhm_lower * new_left_fwhm;
    double new_upper_energy = new_right_mean + roi_width_num_fwhm_upper * new_right_fwhm;

    // Constrain new bounds to not expand beyond original ROI bounds
    // This function is only called when edge peaks are removed, so ROI should only shrink, never expand
    const double orig_lower = old_continuum->lowerEnergy();
    const double orig_upper = old_continuum->upperEnergy();
    new_lower_energy = std::max( new_lower_energy, orig_lower );
    new_upper_energy = std::min( new_upper_energy, orig_upper );

    // Create new continuum as copy with updated energy range
    std::shared_ptr<PeakContinuum> new_continuum = std::make_shared<PeakContinuum>( *old_continuum );
    new_continuum->setRange( new_lower_energy, new_upper_energy );

    // Set new continuum to all peaks
    for( PeakDef &p : peaks )
      p.setContinuum( new_continuum );

    if( PeakFitImprove::debug_printout )
    {
      cout << "  Observable filter: adjusted ROI bounds from ["
           << old_continuum->lowerEnergy() << ", " << old_continuum->upperEnergy()
           << "] to [" << new_lower_energy << ", " << new_upper_energy << "] keV" << endl;
    }

    return true;
  };//reduce_roi_bounds_if_needed lambda

  // Step 1: Group input peaks by ROI (shared continuum)
  std::map<std::shared_ptr<const PeakContinuum>, std::vector<PeakDef>> input_rois;
  for( const PeakDef &peak : fit_peaks )
    input_rois[peak.continuum()].push_back( peak );

  // Step 2: Initial significance filter and ROI adjustment per ROI
  std::vector<PeakDef> filtered_peaks;
  for( auto &roi_entry : input_rois )
  {
    std::vector<PeakDef> &roi_peaks = roi_entry.second;

    // Sort by mean for edge tracking
    std::sort( begin(roi_peaks), end(roi_peaks), &PeakDef::lessThanByMean );

    const double orig_left_mean = roi_peaks.front().mean();
    const double orig_right_mean = roi_peaks.back().mean();

    // Filter peaks by initial significance
    std::vector<PeakDef> kept_peaks;
    for( const PeakDef &peak : roi_peaks )
    {
      const double mean = peak.mean();
      const double fwhm = peak.fwhm();
      const double lower_energy = mean - fwhm;
      const double upper_energy = mean + fwhm;

      // Get total data counts in +/-1 FWHM range
      const double data_area = foreground->gamma_integral( lower_energy, upper_energy );

      // Calculate initial significance
      const double peak_contrib = peak.amplitude() * fwhm_fraction;
      const double significance = peak_contrib / std::sqrt( std::max(data_area, 1.0) );

      if( significance >= initial_significance_threshold )
      {
        kept_peaks.push_back( peak );
      }
      else if( PeakFitImprove::debug_printout )
      {
        cout << "  Observable filter (initial sig=" << significance << " < " << initial_significance_threshold
             << "): peak at " << mean << " keV" << endl;
      }
    }//for( const PeakDef &peak : roi_peaks )

    // Adjust ROI bounds if edge peaks were removed
    if( !kept_peaks.empty() )
    {
      reduce_roi_bounds_if_needed( kept_peaks, orig_left_mean, orig_right_mean );
      filtered_peaks.insert( end(filtered_peaks), begin(kept_peaks), end(kept_peaks) );
    }
  }//for( auto &roi_entry : input_rois )

  if( filtered_peaks.empty() )
    return filtered_peaks;

  // Step 3: Group peaks by shared continuum (ROI) - may have new continuums from adjustment
  std::map<std::shared_ptr<const PeakContinuum>, std::vector<PeakDef>> rois;
  for( const PeakDef &peak : filtered_peaks )
    rois[peak.continuum()].push_back( peak );

  // Debug: Check if rois already have overlaps before parallel processing
#if( PERFORM_DEVELOPER_CHECKS )
  if( PeakFitImprove::debug_printout )
  {
    struct RoiDebugInfo2
    {
      double lower_energy;
      double upper_energy;
      size_t num_peaks;
    };

    std::vector<RoiDebugInfo2> rois_debug2;
    for( const auto &entry : rois )
    {
      RoiDebugInfo2 info;
      info.lower_energy = entry.first->lowerEnergy();
      info.upper_energy = entry.first->upperEnergy();
      info.num_peaks = entry.second.size();
      rois_debug2.push_back( info );
    }

    std::sort( begin(rois_debug2), end(rois_debug2),
      []( const RoiDebugInfo2 &a, const RoiDebugInfo2 &b ) { return a.lower_energy < b.lower_energy; } );

    bool found_overlap2 = false;
    for( size_t i = 1; i < rois_debug2.size(); ++i )
    {
      const RoiDebugInfo2 &prev_roi = rois_debug2[i - 1];
      const RoiDebugInfo2 &curr_roi = rois_debug2[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap2 )
          cerr << "compute_observable_peaks: OVERLAPS BEFORE PARALLEL PROCESSING:" << endl;
        found_overlap2 = true;

        cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy
             << "] keV, " << prev_roi.num_peaks << " peaks" << endl;
        cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy
             << "] keV, " << curr_roi.num_peaks << " peaks" << endl;
        cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << endl;
      }
    }

    if( !found_overlap2 && PeakFitImprove::debug_printout )
      cout << "compute_observable_peaks: No overlaps before parallel processing" << endl;
  }
#endif

  // Step 4: Iteratively refit each ROI and remove insignificant peaks - in parallel
  const size_t num_rois = rois.size();
  std::vector<std::vector<PeakDef>> roi_results( num_rois );

  // Copy ROI data to vector for parallel processing
  std::vector<std::pair<std::shared_ptr<const PeakContinuum>, std::vector<PeakDef>>> roi_vec( rois.begin(), rois.end() );

  SpecUtilsAsync::ThreadPool pool;

  for( size_t roi_index = 0; roi_index < num_rois; ++roi_index )
  {
    pool.post( [roi_index, &roi_vec, &roi_results, &foreground,
                final_significance_threshold, &reduce_roi_bounds_if_needed]()
    {
      std::vector<PeakDef> roi_peaks = roi_vec[roi_index].second;

      // Sort peaks by mean energy for edge detection
      std::sort( begin(roi_peaks), end(roi_peaks), &PeakDef::lessThanByMean );

      const size_t max_iterations = 3;
      size_t iteration = 0;
      bool changed = true;
      while( changed && !roi_peaks.empty() && (iteration < max_iterations) )
      {
        changed = false;

        // Track original edge peaks before filtering
        const double orig_left_mean = roi_peaks.front().mean();
        const double orig_right_mean = roi_peaks.back().mean();

        // Convert to shared_ptr for refitPeaksThatShareROI_LM
        std::vector<std::shared_ptr<const PeakDef>> input_peaks;
        for( const PeakDef &p : roi_peaks )
          input_peaks.push_back( std::make_shared<PeakDef>(p) );

        // Refit with SmallRefinementOnly option
        std::vector<std::shared_ptr<const PeakDef>> refit_result =
          PeakFitLM::refitPeaksThatShareROI_LM( foreground, nullptr, input_peaks,
            PeakFitLM::PeakFitLMOptions::SmallRefinementOnly );

        // If refit failed, keep original peaks
        if( refit_result.empty() )
          break;

        // Check significance and remove insignificant peaks
        std::vector<PeakDef> kept_peaks;
        for( const std::shared_ptr<const PeakDef> &peak : refit_result )
        {
          const double mean = peak->mean();
          const double area = peak->peakArea();
          const double area_uncert = peak->peakAreaUncert();
          const double final_sig = (area_uncert > 0.0) ? (area / area_uncert) : 0.0;

          if( final_sig >= final_significance_threshold )
          {
            kept_peaks.push_back( *peak );
          }
          else
          {
            changed = true;
            if( PeakFitImprove::debug_printout )
            {
              cout << "  Observable filter post-refit (final sig=" << final_sig << " < " << final_significance_threshold
                   << "): peak at " << mean << " keV" << endl;
            }
          }
        }//for( const std::shared_ptr<const PeakDef> &peak : refit_result )

        // Adjust ROI bounds if edge peaks were removed
        if( !kept_peaks.empty() )
        {
          const bool bounds_adjusted = reduce_roi_bounds_if_needed( kept_peaks, orig_left_mean, orig_right_mean );
          if( bounds_adjusted )
            changed = true;  // Need to refit with new continuum bounds
        }

        roi_peaks = std::move( kept_peaks );
        iteration += 1;
      }//while( changed && !roi_peaks.empty() && (iteration < max_iterations) )

      roi_results[roi_index] = std::move( roi_peaks );
    });//pool.post lambda
  }//for( size_t roi_index = 0; roi_index < num_rois; ++roi_index )

  pool.join();

  // Collect all results
  std::vector<PeakDef> observable_peaks;
  for( const std::vector<PeakDef> &roi_result : roi_results )
    observable_peaks.insert( end(observable_peaks), begin(roi_result), end(roi_result) );

  // Sort by energy
  std::sort( begin(observable_peaks), end(observable_peaks), &PeakDef::lessThanByMean );

#if( PERFORM_DEVELOPER_CHECKS )
  // Debug: Check for overlapping ROI bounds
  if( PeakFitImprove::debug_printout )
  {
    // Group peaks by continuum to identify unique ROIs
    std::map<std::shared_ptr<const PeakContinuum>, std::vector<double>> continuum_to_peak_means;
    for( const PeakDef &peak : observable_peaks )
      continuum_to_peak_means[peak.continuum()].push_back( peak.mean() );

    // Extract ROIs sorted by lower energy
    struct RoiDebugInfo
    {
      double lower_energy;
      double upper_energy;
      std::vector<double> peak_means;
    };

    std::vector<RoiDebugInfo> rois_debug;
    for( const auto &entry : continuum_to_peak_means )
    {
      RoiDebugInfo info;
      info.lower_energy = entry.first->lowerEnergy();
      info.upper_energy = entry.first->upperEnergy();
      info.peak_means = entry.second;
      rois_debug.push_back( info );
    }

    std::sort( begin(rois_debug), end(rois_debug),
      []( const RoiDebugInfo &a, const RoiDebugInfo &b ) { return a.lower_energy < b.lower_energy; } );

    // Check for overlaps
    bool found_overlap = false;
    for( size_t i = 1; i < rois_debug.size(); ++i )
    {
      const RoiDebugInfo &prev_roi = rois_debug[i - 1];
      const RoiDebugInfo &curr_roi = rois_debug[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap )
          cerr << "compute_observable_peaks: FOUND OVERLAPPING ROIs:" << endl;
        found_overlap = true;

        cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : prev_roi.peak_means )
          cerr << mean << " ";
        cerr << "keV" << endl;

        cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : curr_roi.peak_means )
          cerr << mean << " ";
        cerr << "keV" << endl;
        cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << endl;
      }
    }

    if( !found_overlap && PeakFitImprove::debug_printout )
      cout << "compute_observable_peaks: No overlapping ROIs detected" << endl;
  }//if( PeakFitImprove::debug_printout )
#endif

  return observable_peaks;
}//compute_observable_peaks

  
/** Finds the valid energy range of a foreground spectrum based on spectroscopic extent
   
   @param foreground The spectrum to analyze
   
   @return A pair of (min_valid_energy, max_valid_energy). Returns (0.0, 0.0) if foreground is
   null or has no data.
*/
pair<double,double> find_valid_energy_range( const shared_ptr<const SpecUtils::Measurement> &meas )
{
  // This implementation is an updated implementation of `find_spectroscopic_extent(...)`,
  //  and that function should probably be updated
  if( !meas 
    || !meas->energy_calibration() 
    || !meas->energy_calibration()->valid()
    || (meas->num_gamma_channels() < 7) )
  {
    return {0.0, 0.0};
  }
  
  size_t lower_channel = 0, upper_channel = 0;

  const vector<float> &channel_counts = *meas->gamma_counts();
  const size_t nbin = channel_counts.size();
  
  //First detect where spectrum begins
  const int side_bins = 3;
  const int order = 2;
  const int derivative = 2;
  vector<float> smoothed_2nd, smoothed_2nd_variance;
  SavitzyGolayCoeffs sgcoeffs( side_bins, side_bins, order, derivative );
  sgcoeffs.smooth_with_variance( channel_counts, smoothed_2nd, smoothed_2nd_variance );

  // Find where spectrum begins by looking for significant negative curvature (detector turn-on)
  // followed by the curvature leveling off, using variance to set thresholds
  size_t channel = 0;

  // First, find where we have significant negative curvature (downturn after detector turn-on)
  // Use variance to set a threshold - we're looking for a signal that's statistically significant
  while( channel < nbin )
  {
    const float sigma = (smoothed_2nd_variance[channel] > 0.0f) ? std::sqrt(smoothed_2nd_variance[channel]) : 1.0f;
    // Look for curvature more negative than -2 sigma (statistically significant downturn)
    if( smoothed_2nd[channel] < -2.0f * sigma )
      break;
    ++channel;
  }

  // Now find where the curvature levels off (spectrum becomes relatively flat/linear)
  // Again use variance to determine when we're close enough to zero
  while( channel < nbin )
  {
    const float sigma = (smoothed_2nd_variance[channel] > 0.0f) ? std::sqrt(smoothed_2nd_variance[channel]) : 1.0f;
    // Look for curvature within 1 sigma of zero (spectrum has leveled off)
    if( std::abs(smoothed_2nd[channel]) < sigma )
      break;
    ++channel;
  }
  
  lower_channel = std::min(channel-0,smoothed_2nd.size()-1);
  
  // Find the upper energy limit by looking for the last channel with meaningful signal
  // Start from the end and work backwards
  size_t upperlastchannel = nbin - 1;

  // First, skip trailing channels with very low counts (< 3 counts)
  while( (upperlastchannel > 0) && (channel_counts[upperlastchannel] < 3.0f) )
    --upperlastchannel;

  // Now increase upperlastchannel until we find at least two consecutive zero-count channels
  bool found_two_zeros = false;
  while( (upperlastchannel < nbin - 2) && !found_two_zeros )
  {
    ++upperlastchannel;

    // Check if we have two consecutive zeros starting at upperlastchannel
    if( (channel_counts[upperlastchannel] == 0.0f)
        && (channel_counts[upperlastchannel + 1] == 0.0f) )
    {
      found_two_zeros = true;
      ++upperlastchannel; // Set to one more than the second zero
    }
  }

  // Ensure upperlastchannel is valid
  upperlastchannel = std::min( upperlastchannel, nbin - 1 );


  //  upperlastchannel = (std::min)( upperlastchannel + 5, nbin - 1 );
  
  // Search backwards from the end to find the last peak-like structure
  // Look for significant negative curvature, then find where it becomes positive on the right side
  size_t lastchannel = upperlastchannel;

  bool found_last_signal = false;
  for( size_t i = upperlastchannel; i > lower_channel; --i )
  {
    // Look for statistically significant negative curvature (peak shape)
    if( (i >= smoothed_2nd.size()) || (smoothed_2nd_variance[i] <= 0.0f) )
      continue;

    const float sigma_2nd = std::sqrt( smoothed_2nd_variance[i] );
    if( smoothed_2nd[i] < (-2.0f * sigma_2nd) )
    {
      // Found significant negative curvature - this is a potential peak
      // Now find where the second derivative becomes positive on the right side
      size_t right_extent = i;
      for( size_t j = i + 1; j < smoothed_2nd.size(); ++j )
      {
        if( smoothed_2nd[j] > 0.0f )
        {
          right_extent = j;
          break;
        }
      }

      // Calculate peak width and add ~25% buffer
      const size_t peak_width = (right_extent > i) ? (right_extent - i) : 1;
      const size_t buffer = std::max( size_t(4), peak_width / 4 );
      lastchannel = std::min( right_extent + buffer, nbin - 1 );
      found_last_signal = true;
      break;
    }
  }

  // If we didn't find peak structure in the backwards search, use the simple upperlastchannel
  if( !found_last_signal )
    lastchannel = upperlastchannel;

  // If the sophisticated detection failed, fall back to a simple algorithm
  if( (lower_channel > (nbin/3)) || (lastchannel <= lower_channel) )
  {
    // Find first channel where this channel and the next three all have non-zero counts
    size_t first_nonzero = 0;
    while( first_nonzero + 3 < nbin )
    {
      if( (channel_counts[first_nonzero] > 0.0f)
         && (channel_counts[first_nonzero + 1] > 0.0f)
         && (channel_counts[first_nonzero + 2] > 0.0f)
         && (channel_counts[first_nonzero + 3] > 0.0f) )
      {
        break;
      }
      ++first_nonzero;
    }

    // Find last channel where this channel and the previous channel both have non-zero counts
    size_t last_nonzero = nbin - 1;
    while( last_nonzero > 0 )
    {
      if( (channel_counts[last_nonzero] > 0.0f) && (channel_counts[last_nonzero - 1] > 0.0f) )
        break;
      --last_nonzero;
    }

    // Make sure we found valid channels
    if( first_nonzero < nbin && last_nonzero > first_nonzero )
    {
      lower_channel = first_nonzero;
      upper_channel = last_nonzero;
    }
    else
    {
      // Ultimate fallback - use entire spectrum
      lower_channel = 0;
      upper_channel = nbin - 1;
    }

    return { meas->gamma_channel_lower(lower_channel), meas->gamma_channel_upper(upper_channel) };
  }//

  // Add some buffer beyond the last detected signal to ensure we don't cut off peak tails
  const size_t buffer_channels = std::max( size_t(10), size_t(0.005 * nbin) );
  upper_channel = std::min( lastchannel + buffer_channels, nbin - 1 );

  // Make sure we don't cut off too early - at minimum use upperlastchannel
  upper_channel = std::max( upper_channel, upperlastchannel );
  upper_channel = std::min( upper_channel, nbin - 1 );

  return { meas->gamma_channel_lower(lower_channel), meas->gamma_channel_upper(upper_channel) };
}//find_valid_energy_range
  

// Helper to check if two ROI sets are similar enough to stop iterating
bool rois_are_similar( const vector<RelActCalcAuto::RoiRange> &a,
                       const vector<RelActCalcAuto::RoiRange> &b,
                       const double energy_tolerance = 1.0 )
{
  if( a.size() != b.size() )
    return false;

  for( size_t i = 0; i < a.size(); ++i )
  {
    if( fabs(a[i].lower_energy - b[i].lower_energy) > energy_tolerance )
      return false;
    if( fabs(a[i].upper_energy - b[i].upper_energy) > energy_tolerance )
      return false;
  }
  return true;
}//rois_are_similar


/** Result of computing chi2 significance for a single ROI.
 */
struct RoiSignificanceResult
{
  double chi2_with_peaks = 0.0;
  double chi2_continuum_only = 0.0;
  double chi2_reduction = 0.0;  // chi2_continuum_only - chi2_with_peaks
  double max_peak_significance = 0.0;  // Maximum peak_area / sqrt(continuum) for any peak
  size_t num_channels = 0;
  bool passes_chi2_test = false;       // chi2_reduction >= threshold
  bool passes_peak_sig_test = false;   // max_peak_significance >= threshold
  bool has_significant_peaks = false;  // passes_chi2_test OR passes_peak_sig_test
};//struct RoiSignificanceResult


/** Computes chi2 significance for a single ROI by comparing chi2 with peaks vs continuum-only.

 This function:
 1. Finds peaks within the ROI
 2. Computes chi2 with those peaks using chi2_for_region()
 3. Fits a continuum-only model using fit_amp_and_offset()
 4. Computes total chi2 reduction
 5. Computes peak significance (peak_area / sqrt(continuum)) for each peak
 6. Determines if peaks are significant based on either threshold

 @param roi The ROI range to evaluate
 @param all_peaks All peaks in the solution
 @param data The spectrum measurement
 @param min_chi2_reduction Minimum chi2 reduction to consider peaks significant
 @param min_peak_significance Minimum peak significance to consider a peak significant
 @return RoiSignificanceResult with chi2 values and significance determination
 */
RoiSignificanceResult compute_roi_chi2_significance(
  const RelActCalcAuto::RoiRange &roi,
  const std::vector<PeakDef> &all_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const double min_chi2_reduction,
  const double min_peak_significance )
{
  RoiSignificanceResult result;

  // Find peaks in this ROI
  std::vector<std::shared_ptr<const PeakDef>> peaks_in_roi;
  for( const PeakDef &peak : all_peaks )
  {
    if( peak.mean() >= roi.lower_energy && peak.mean() <= roi.upper_energy )
      peaks_in_roi.push_back( std::make_shared<PeakDef>( peak ) );
  }

  if( peaks_in_roi.empty() )
    return result;

  // Get channel range
  const size_t start_channel = data->find_gamma_channel( static_cast<float>( roi.lower_energy ) );
  const size_t end_channel = data->find_gamma_channel( static_cast<float>( roi.upper_energy ) );
  result.num_channels = (end_channel > start_channel) ? (end_channel - start_channel) : 0;

  if( result.num_channels < 3 )
    return result;


  // Get continuum parameters from first peak
  const std::shared_ptr<const PeakContinuum> continuum = peaks_in_roi[0]->continuum();
  if( !continuum )
    return result;

  const int num_poly_terms = static_cast<int>( continuum->parameters().size() );
  const bool step_continuum = (continuum->type() == PeakContinuum::OffsetType::FlatStep
                            || continuum->type() == PeakContinuum::OffsetType::LinearStep
                            || continuum->type() == PeakContinuum::OffsetType::BiLinearStep);

  // Prepare data arrays for fit_amp_and_offset
  const std::vector<float> &channel_energies = *data->channel_energies();
  std::vector<float> channel_counts( result.num_channels );
  for( size_t i = 0; i < result.num_channels; ++i )
    channel_counts[i] = data->gamma_channel_content( start_channel + i );

  std::vector<double> empty_means, empty_sigmas;
  std::vector<double> continuum_coeffs, continuum_uncerts;
  std::vector<double> dummy_amps, dummy_amp_uncerts;

  // To compare apples to apples, call `fit_amp_and_offset(...)` for `peaks_in_roi`,
  // but with the amplitudes of the peaks fixed, so the continuum is re-fit, but not the peak amplitudes.
  // This way we can use the chi2 computed by fit_amp_and_offset, to consistently compare against
  // the chi2 computed for the continuum fit for no peaks.
  std::vector<PeakDef> fixed_peaks;
  for( const std::shared_ptr<const PeakDef> &peak : peaks_in_roi )
    fixed_peaks.push_back( *peak );

  result.chi2_with_peaks = fit_amp_and_offset(
    &channel_energies[start_channel],
    channel_counts.data(),
    result.num_channels,
    num_poly_terms,
    step_continuum,
    continuum->referenceEnergy(),
    empty_means,
    empty_sigmas,
    fixed_peaks,  // Pass peaks as fixed - their amplitudes won't be fit
    PeakDef::SkewType::NoSkew,
    nullptr,
    dummy_amps,
    continuum_coeffs,
    dummy_amp_uncerts,
    continuum_uncerts
  );

  // Fit continuum only (no peaks)
  std::vector<PeakDef> empty_fixed_peaks;

  result.chi2_continuum_only = fit_amp_and_offset(
    &channel_energies[start_channel],
    channel_counts.data(),
    result.num_channels,
    num_poly_terms,
    step_continuum,
    continuum->referenceEnergy(),
    empty_means,
    empty_sigmas,
    empty_fixed_peaks,
    PeakDef::SkewType::NoSkew,
    nullptr,
    dummy_amps,
    continuum_coeffs,
    dummy_amp_uncerts,
    continuum_uncerts
  );

  // Compute chi2 reduction
  result.chi2_reduction = result.chi2_continuum_only - result.chi2_with_peaks;
  result.passes_chi2_test = (result.chi2_reduction >= min_chi2_reduction);

  // Compute peak significance for each peak: peak_area / sqrt(continuum)
  // For each peak, integrate the continuum between mean ± 1 FWHM
  // Peak amplitude in this range is ~97.93% of total (for a Gaussian, erf(2.355/sqrt(2)) ≈ 0.9793)
  const double fwhm_coverage_fraction = 0.9793;

  for( const std::shared_ptr<const PeakDef> &peak : peaks_in_roi )
  {
    const double peak_mean = peak->mean();
    const double peak_fwhm = peak->fwhm();
    const double peak_lower = peak_mean - peak_fwhm;
    const double peak_upper = peak_mean + peak_fwhm;

    // Integrate continuum between peak_lower and peak_upper
    const std::shared_ptr<const PeakContinuum> peak_cont = peak->continuum();
    if( !peak_cont )
      continue;

    // Continuum integral, with a minimum of 1.0 to avoid division issues
    const double continuum_integral = std::max( 1.0, peak_cont->offset_integral( peak_lower, peak_upper, data ) );

    // Peak amplitude in the ±1 FWHM range
    const double peak_amp_in_range = peak->amplitude() * fwhm_coverage_fraction;

    // Significance = peak_area / sqrt(continuum)
    const double peak_sig = peak_amp_in_range / std::sqrt( continuum_integral );
    
    if( peak_sig > result.max_peak_significance )
      result.max_peak_significance = peak_sig;
  }//for( loop over peaks in ROI )

  result.passes_peak_sig_test = (result.max_peak_significance >= min_peak_significance);
  result.has_significant_peaks = (result.passes_chi2_test || result.passes_peak_sig_test);

  return result;
}//compute_roi_chi2_significance


/** Computes filtered chi2/dof by only including ROIs with significant peaks.

 This function loops through all ROIs in the solution, evaluates each for significance,
 and returns the chi2/dof considering only the significant ROIs.

 @param solution The RelActAuto solution containing ROIs and peaks
 @param data The spectrum measurement
 @param min_chi2_reduction Minimum chi2 reduction for ROI to be significant
 @param min_peak_significance Minimum peak significance for ROI to be significant
 @param insignificant_roi_indices Output vector of ROI indices that were insignificant
 @return Chi2/DOF considering only significant ROIs
 */
double compute_filtered_chi2_dof(
  const RelActCalcAuto::RelActAutoSolution &solution,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const double min_chi2_reduction,
  const double min_peak_significance,
  std::vector<size_t> &insignificant_roi_indices )
{
  insignificant_roi_indices.clear();

  double total_chi2 = 0.0;
  size_t total_channels = 0;

  for( size_t roi_idx = 0; roi_idx < solution.m_final_roi_ranges.size(); ++roi_idx )
  {
    const RelActCalcAuto::RoiRange &roi = solution.m_final_roi_ranges[roi_idx];

    const RoiSignificanceResult sig_result = compute_roi_chi2_significance(
      roi, solution.m_peaks_without_back_sub, data, min_chi2_reduction, min_peak_significance );

    if( sig_result.has_significant_peaks )
    {
      total_chi2 += sig_result.chi2_with_peaks;
      total_channels += sig_result.num_channels;

      if( PeakFitImprove::debug_printout )
      {
        cout << "  Significant ROI " << roi_idx << " [" << roi.lower_energy << ", " << roi.upper_energy << "] keV: "
             << "chi2_test=" << (sig_result.passes_chi2_test ? "PASS" : "fail")
             << " (reduction=" << sig_result.chi2_reduction << ", thresh=" << min_chi2_reduction << "), "
             << "peak_sig_test=" << (sig_result.passes_peak_sig_test ? "PASS" : "fail")
             << " (max_sig=" << sig_result.max_peak_significance << ", thresh=" << min_peak_significance << ")" << endl;

        // Print peaks in this ROI
        for( const PeakDef &peak : solution.m_peaks_without_back_sub )
        {
          const double peak_roi_lower = peak.continuum()->lowerEnergy();
          const double peak_roi_upper = peak.continuum()->upperEnergy();
          if( (fabs( peak_roi_lower - roi.lower_energy ) < 1.0)
             && (fabs( peak_roi_upper - roi.upper_energy ) < 1.0) )
          {
            cout << "    Peak: mean=" << peak.mean() << " keV, area=" << peak.peakArea() << endl;
          }
        }//for( loop over peaks )
      }//if( debug_printout )
    }
    else
    {
      insignificant_roi_indices.push_back( roi_idx );

      if( PeakFitImprove::debug_printout )
      {
        cout << "  Insignificant ROI " << roi_idx << " [" << roi.lower_energy << ", " << roi.upper_energy << "] keV: "
             << "chi2_test=" << (sig_result.passes_chi2_test ? "PASS" : "fail")
             << " (reduction=" << sig_result.chi2_reduction << ", thresh=" << min_chi2_reduction << "), "
             << "peak_sig_test=" << (sig_result.passes_peak_sig_test ? "PASS" : "fail")
             << " (max_sig=" << sig_result.max_peak_significance << ", thresh=" << min_peak_significance << ")" << endl;

        // Print peaks in this ROI
        for( const PeakDef &peak : solution.m_peaks_without_back_sub )
        {
          const double peak_roi_lower = peak.continuum()->lowerEnergy();
          const double peak_roi_upper = peak.continuum()->upperEnergy();
          if( (fabs( peak_roi_lower - roi.lower_energy ) < 1.0)
             && (fabs( peak_roi_upper - roi.upper_energy ) < 1.0) )
          {
            cout << "    Peak: mean=" << peak.mean() << " keV, area=" << peak.peakArea() << endl;
          }
        }//for( loop over peaks )
      }//if( debug_printout )
    }
  }//for( loop over ROIs )

  if( total_channels == 0 )
    return std::numeric_limits<double>::max();

  return total_chi2 / static_cast<double>( total_channels );
}//compute_filtered_chi2_dof


/** Fits peaks for nuclides using RelActAuto.

 This function encapsulates the peak fitting workflow using RelActAuto.
 It takes configuration parameters and performs the fit, including iterative refinement.

 @param auto_search_peaks Initial peaks found by automatic search
 @param foreground Foreground spectrum to fit
 @param sources Vector of source nuclides to fit
 @param options RelActAuto options configuration (created from config)
 @param long_background Long background spectrum (can be nullptr)
 @param drf Detector response function
 @param config Configuration for peak fitting parameters
 @return PeakFitResult with status, error message, fit peaks, and solution
 */
PeakFitResult fit_peaks_for_nuclide_relactauto(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::vector<RelActCalcAuto::RoiRange> &initial_rois,
  const std::shared_ptr<const SpecUtils::Measurement> &long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitForNuclideConfig &config,
  const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
  const std::vector<float> &fwhm_coefficients,
  const double fwhm_lower_energy,
  const double fwhm_upper_energy )
{
  PeakFitResult result;

  // Validate sources
  if( sources.empty() )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "No sources provided";
    return result;
  }

  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    if( RelActCalcAuto::is_null( src ) )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Null source in sources vector";
      return result;
    }
  }

  // Build base_nuclides from sources
  std::vector<RelActCalcAuto::NucInputInfo> base_nuclides;
  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    RelActCalcAuto::NucInputInfo nuc_info;
    nuc_info.source = src;
    nuc_info.age = get_source_age( src, -1.0 );
    nuc_info.fit_age = false;
    base_nuclides.push_back( nuc_info );
  }

  // Define skew type to use
  const PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;

  // Create relative efficiency curve from config parameters
  RelActCalcAuto::RelEffCurveInput rel_eff_curve;
  rel_eff_curve.rel_eff_eqn_type = config.rel_eff_eqn_type;
  rel_eff_curve.rel_eff_eqn_order = config.rel_eff_eqn_order;

  // FramPhysicalModel requires rel_eff_eqn_order to be 0
  assert( (config.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel)
         || (rel_eff_curve.rel_eff_eqn_order == 0) );
  if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    rel_eff_curve.rel_eff_eqn_order = 0;

  rel_eff_curve.nucs_of_el_same_age = config.nucs_of_el_same_age;
  rel_eff_curve.phys_model_use_hoerl = config.phys_model_use_hoerl;
  rel_eff_curve.nuclides = base_nuclides;

  // Copy shielding inputs from config (converting to const shared_ptr)
  if( !config.phys_model_self_atten.empty() )
    rel_eff_curve.phys_model_self_atten = config.phys_model_self_atten.front();

  rel_eff_curve.phys_model_external_atten.clear();
  for( const auto &shield : config.phys_model_external_atten )
    rel_eff_curve.phys_model_external_atten.push_back( shield );

  // Generate name based on equation type and shielding (just for informational purposes)
  if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    rel_eff_curve.name = "Physical Model";
    if( config.phys_model_external_atten.empty() )
    {
      rel_eff_curve.name += " (no shielding) Peak Fit";
    }else
    {
      const int z = static_cast<int>( config.phys_model_external_atten.front()->atomic_number + 0.5 );
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      const SandiaDecay::Element *el = db->element(z);
      assert( el );
      rel_eff_curve.name += " " + (el ? el->symbol : ("z=" + std::to_string(z)) ) + string(" Peak Fit");
    }
  }else
  {
    rel_eff_curve.name = RelActCalc::to_str(config.rel_eff_eqn_type) + std::string(" Peak Fit");
  }

  // Create RelActAuto options from config
  RelActCalcAuto::Options options;
  options.rel_eff_curves.push_back( rel_eff_curve );
  options.rois = initial_rois;
  options.fit_energy_cal = config.fit_energy_cal;
  options.fwhm_form = config.fwhm_form;
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
  options.skew_type = skew_type;
  options.additional_br_uncert = config.rel_eff_auto_base_rel_eff_uncert;

  // Find valid energy range based on contiguous channels with data
  const auto [min_valid_energy, max_valid_energy] = find_valid_energy_range( foreground );

  try
  {
    // Call RelActAuto::solve with provided options
    RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
      options, foreground, long_background, drf, auto_search_peaks
    );

    // As of 20260103, energy calibration adjustments may cause failure to fit the correct solution sometimes,
    //  so if our current solution failed, or is really bad, we'll try without fitting energy cal
    if( options.fit_energy_cal
      && ((solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success)
        || ((solution.m_chi2 / solution.m_dof) > 10.0)) ) //10.0 arbitrary - and un-explored
    {
      RelActCalcAuto::Options no_ecal_opts = options;
      no_ecal_opts.fit_energy_cal = false;

      RelActCalcAuto::RelActAutoSolution trial_solution = RelActCalcAuto::solve(
        no_ecal_opts, foreground, long_background, drf, auto_search_peaks
      );
      
      // If the solution is still really bad - we'll try a Physical Model solution
      // Optionally with external shielding if configured
      if( (trial_solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success)
          || ((solution.m_chi2 / solution.m_dof) > 10.0) )
      {
        RelActCalcAuto::Options desperation_opts = options;
        RelActCalcAuto::RelEffCurveInput &curve = desperation_opts.rel_eff_curves.front();
        curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
        curve.rel_eff_eqn_order = 0;
        curve.phys_model_self_atten = nullptr;
        curve.phys_model_external_atten.clear();

        // Apply external shielding if conditions are met
        if( should_use_desperation_shielding( config.desperation_phys_model_atomic_number, options.rois ) )
        {
          try
          {
            std::shared_ptr<RelActCalc::PhysicalModelShieldInput> shield
              = create_desperation_shielding( config.desperation_phys_model_atomic_number,
                                              config.desperation_phys_model_areal_density_g_per_cm2 );
            curve.phys_model_external_atten.push_back( shield );

            if( PeakFitImprove::debug_printout )
            {
              std::cerr << "First desperation attempt: using external shielding with AN="
                        << config.desperation_phys_model_atomic_number
                        << ", starting AD=" << config.desperation_phys_model_areal_density_g_per_cm2
                        << " g/cm2" << std::endl;
            }
          }
          catch( const std::exception &e )
          {
            if( PeakFitImprove::debug_printout )
              std::cerr << "Failed to create desperation shielding: " << e.what() << std::endl;
            // Continue without shielding
          }
        }
        else
        {
          if( PeakFitImprove::debug_printout )
            std::cerr << "First desperation attempt: not using external shielding" << std::endl;
        }

        curve.phys_model_use_hoerl = (options.rois.size() > 2);

        trial_solution = RelActCalcAuto::solve(
          desperation_opts, foreground, long_background, drf, auto_search_peaks
        );
      }//If( still a bad solution )
      
      if( (trial_solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success)
        && ( (solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success)
            || ((solution.m_chi2 / solution.m_dof) > (trial_solution.m_chi2 / trial_solution.m_dof)) ) )
      {
        if( PeakFitImprove::debug_printout )
          cerr << "Abandoning fitting e-cal for nuclide" << endl;
        solution = trial_solution;
      }
    }

    // Check if initial solve failed
    if( solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      result.status = solution.m_status;
      result.error_message = solution.m_error_message;
      return result;
    }

    cout << "Initial RelActAuto solution (" << options.rel_eff_curves.front().name << "):" << endl;
    solution.print_summary( std::cout );
    cout << "Chi2/DOF = " << solution.m_chi2 << "/" << solution.m_dof << " = " << (solution.m_chi2 / solution.m_dof) << endl;

    // Iteratively refine ROIs using RelActAuto solutions
    // The idea is that each iteration provides a better relative efficiency estimate,
    // which allows us to better identify significant gamma lines and create better ROIs.
    {
      const size_t max_iterations = 3;
      size_t num_extra_allowed = 0; //If we switch to our "desperation" model type retry - we will increment this to 1.
      for( size_t iter = 0; iter < (max_iterations + num_extra_allowed); ++iter )
      {
        // Use the first rel-eff curve (index 0)
        const size_t rel_eff_index = 0;

        // Create rel_eff lambda from current RelActAuto solution
        const auto auto_rel_eff = [&solution, rel_eff_index]( double energy ) -> double {
          return solution.relative_efficiency( energy, rel_eff_index );
        };

        // Collect sources and activities from the current solution
        // m_rel_activities is a 2D vector: [rel_eff_curve_index][source_index]
        vector<pair<RelActCalcAuto::SrcVariant, double>> sources_and_acts;
        if( rel_eff_index < solution.m_rel_activities.size() )
        {
          if( PeakFitImprove::debug_printout )
            cout << "Collecting " << solution.m_rel_activities[rel_eff_index].size() << " sources from RelActAuto solution:" << endl;

          for( const RelActCalcAuto::NuclideRelAct &nuc_act : solution.m_rel_activities[rel_eff_index] )
          {
            if( RelActCalcAuto::is_null( nuc_act.source ) )
              continue;

            const double live_time_seconds = foreground->live_time();
            // RelActAuto's rel_activity is per second, need to multiply by live time for clustering
            const double activity_for_clustering = nuc_act.rel_activity * live_time_seconds;

            if( PeakFitImprove::debug_printout )
              cout << "  " << nuc_act.name() << ": rel_activity=" << nuc_act.rel_activity
                   << ", live_time=" << live_time_seconds
                   << "s, activity_for_clustering=" << activity_for_clustering << endl;

            sources_and_acts.emplace_back( nuc_act.source, activity_for_clustering );
          }
        }

        // Use the FWHM coefficients passed to this function (computed from auto-search peaks
        // or DRF in fit_peaks_for_nuclides), rather than relying on solution.m_drf which may
        // not have valid FWHM info or may have incorrect values

        // Get auto clustering settings from config
        const GammaClusteringSettings auto_settings = config.get_auto_clustering_settings();

        // Cluster gammas using current solution's relative efficiency
        vector<RelActCalcAuto::RoiRange> refined_rois = cluster_gammas_to_rois(
            auto_rel_eff, sources_and_acts, foreground,
            fwhm_form, fwhm_coefficients,
            fwhm_lower_energy, fwhm_upper_energy,
            min_valid_energy, max_valid_energy,
            auto_settings );
        
        if( refined_rois.empty() )
        {
          // If we lost all ROIs and are not already using a PhysicalModel, try re-fitting with one
          // This is similar to the "desperation" approach used above
          const bool using_physical_model = (solution.m_options.rel_eff_curves.front().rel_eff_eqn_type
                                             == RelActCalc::RelEffEqnForm::FramPhysicalModel);

          if( !using_physical_model )
          {
            if( PeakFitImprove::debug_printout )
              cerr << "Lost all ROIs, trying PhysicalModel as desperation attempt..." << endl;

            RelActCalcAuto::Options desperation_opts = options;
            RelActCalcAuto::RelEffCurveInput &curve = desperation_opts.rel_eff_curves.front();
            curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
            curve.rel_eff_eqn_order = 0;
            curve.phys_model_self_atten = nullptr;
            curve.phys_model_external_atten.clear();

            // Apply external shielding if conditions are met
            // Note: We use desperation_opts.rois (not options.rois) since they may differ at this point
            if( should_use_desperation_shielding( config.desperation_phys_model_atomic_number, desperation_opts.rois ) )
            {
              try
              {
                std::shared_ptr<RelActCalc::PhysicalModelShieldInput> shield
                  = create_desperation_shielding( config.desperation_phys_model_atomic_number,
                                                  config.desperation_phys_model_areal_density_g_per_cm2 );
                curve.phys_model_external_atten.push_back( shield );

                if( PeakFitImprove::debug_printout )
                {
                  std::cerr << "Second desperation attempt: using external shielding with AN="
                            << config.desperation_phys_model_atomic_number
                            << ", starting AD=" << config.desperation_phys_model_areal_density_g_per_cm2
                            << " g/cm2" << std::endl;
                }
              }
              catch( const std::exception &e )
              {
                if( PeakFitImprove::debug_printout )
                  std::cerr << "Failed to create desperation shielding: " << e.what() << std::endl;
                // Continue without shielding
              }
            }
            else
            {
              if( PeakFitImprove::debug_printout )
                std::cerr << "Second desperation attempt: not using external shielding" << std::endl;
            }

            curve.phys_model_use_hoerl = (desperation_opts.rois.size() > 2);

            RelActCalcAuto::RelActAutoSolution desperation_solution = RelActCalcAuto::solve(
              desperation_opts, foreground, long_background, drf, auto_search_peaks
            );

            if( (desperation_solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success)
                && ((desperation_solution.m_chi2 / desperation_solution.m_dof)
                    < (solution.m_chi2 / solution.m_dof)) )
            {
              if( PeakFitImprove::debug_printout )
                cerr << "PhysicalModel desperation solution succeeded and improved chi2/dof" << endl;
              
              if( !num_extra_allowed ) //Allow an extra iteration since we changed
                num_extra_allowed += 1;
              
              solution = desperation_solution;
              // Continue with another iteration using the new solution
              continue;
            }
            else
            {
              if( PeakFitImprove::debug_printout )
                cerr << "PhysicalModel desperation solution did not improve result" << endl;
            }
          }//if( !using_physical_model )

          result.warnings.push_back( "Lost all ROIs while iterationing to refine solution - stopped early." );
          cerr << "Have lost all ROIs!  Halting iterations to refine solution." << endl;
          break;
        }//if( refined_rois.empty() )
        
        // Debug output: print refined ROIs with expected counts
        if( PeakFitImprove::debug_printout )
        {
          cout << "Iteration " << iter << " refined ROIs:" << endl;
          for( size_t roi_idx = 0; roi_idx < refined_rois.size(); ++roi_idx )
          {
            const RelActCalcAuto::RoiRange &roi = refined_rois[roi_idx];
            const double roi_data_counts = foreground->gamma_integral(
                static_cast<float>(roi.lower_energy), static_cast<float>(roi.upper_energy) );
            cout << "  ROI " << roi_idx << ": [" << roi.lower_energy << " - " << roi.upper_energy
                 << "] keV, width=" << (roi.upper_energy - roi.lower_energy)
                 << " keV, data_counts=" << roi_data_counts << endl;
          }
        }

        // Check if ROIs changed significantly - if not, stop iterating
        if( rois_are_similar( refined_rois, solution.m_options.rois ) )
        {
          cout << "Iteration " << iter << ": ROIs are similar, stopping refinement" << endl;
          break;
        }

        cout << "Iteration " << iter << ": trying " << refined_rois.size() << " refined ROIs" << endl;

        // Re-run RelActAuto with refined ROIs
        RelActCalcAuto::Options refined_options = solution.m_options;
        refined_options.rois = refined_rois;

        RelActCalcAuto::RelActAutoSolution refined_solution
            = RelActCalcAuto::solve( refined_options, foreground, long_background, drf, auto_search_peaks );

        if( refined_solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        {
          cout << "Iteration " << iter << " failed: " << refined_solution.m_error_message << endl;
          break;
        }

        // Developer check: Validate refined solution's final ROIs don't overlap
#if( PERFORM_DEVELOPER_CHECKS )
        for( size_t i = 1; i < refined_solution.m_final_roi_ranges.size(); ++i )
        {
          const RelActCalcAuto::RoiRange &prev_roi = refined_solution.m_final_roi_ranges[i - 1];
          const RelActCalcAuto::RoiRange &curr_roi = refined_solution.m_final_roi_ranges[i];
          if( curr_roi.lower_energy < prev_roi.upper_energy )
          {
            cerr << "ERROR: RelActAuto returned overlapping ROIs[" << (i-1) << "] and [" << i << "]: "
                 << "[" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] vs "
                 << "[" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "]" << endl;
            assert( curr_roi.lower_energy >= prev_roi.upper_energy );
          }
        }
#endif

        // Compute filtered chi2/dof that only includes ROIs with significant peaks.
        // This avoids the problem where adding a ROI in a flat region (with no real peaks)
        // would artificially reduce chi2/dof.
        std::vector<size_t> old_insignificant_rois, new_insignificant_rois;
        const double old_chi2_dof = compute_filtered_chi2_dof(
          solution, foreground, config.roi_significance_min_chi2_reduction,
          config.roi_significance_min_peak_sig, old_insignificant_rois );
        const double new_chi2_dof = compute_filtered_chi2_dof(
          refined_solution, foreground, config.roi_significance_min_chi2_reduction,
          config.roi_significance_min_peak_sig, new_insignificant_rois );

        // Check if chi2/dof improved
        if( new_chi2_dof >= old_chi2_dof )
        {
          cout << "Iteration " << iter << " did not improve filtered chi2/dof ("
               << old_chi2_dof << " -> " << new_chi2_dof << "), stopping" << endl;
          if( !new_insignificant_rois.empty() )
            cout << "  (" << new_insignificant_rois.size() << " ROIs had insignificant peaks)" << endl;
          break;
        }

        solution = std::move( refined_solution );
        cout << "Iteration " << iter << " improved: chi2/dof=" << new_chi2_dof
             << " (was " << old_chi2_dof << ")" << endl;
      }//for( size_t iter = 0; iter < max_iterations; ++iter )

      cout << "Final solution after refinement:" << endl;
      solution.print_summary( std::cout );
      cout << endl;

      // Print ROIs and sum fit peak areas for each ROI
      cout << "Solution ROIs and fit peak areas:" << endl;
      for( size_t roi_index = 0; roi_index < solution.m_final_roi_ranges.size(); ++roi_index )
      {
        const RelActCalcAuto::RoiRange &roi = solution.m_final_roi_ranges[roi_index];
        double sum_peak_area = 0.0;
        size_t num_peaks_in_roi = 0;

        for( const PeakDef &peak : solution.m_fit_peaks )
        {
          const double peak_roi_lower = peak.continuum()->lowerEnergy();
          const double peak_roi_upper = peak.continuum()->upperEnergy();
          // Match if peak's ROI bounds are within 1 keV of the solution ROI bounds
          if( (fabs( peak_roi_lower - roi.lower_energy ) < 1.0)
             && (fabs( peak_roi_upper - roi.upper_energy ) < 1.0) )
          {
            sum_peak_area += peak.peakArea();
            ++num_peaks_in_roi;
          }
        }//for( loop over fit peaks )

        cout << "  ROI " << roi_index << ": [" << roi.lower_energy << ", " << roi.upper_energy << "] keV"
             << ", " << num_peaks_in_roi << " peaks, sum area = " << sum_peak_area << endl;
      }//for( loop over ROIs )
      cout << endl;
    }//iterative refinement

    // Identify ROIs without significant peaks for filtering
    std::vector<size_t> final_insignificant_rois;
    compute_filtered_chi2_dof( solution, foreground,
      config.roi_significance_min_chi2_reduction, config.roi_significance_min_peak_sig,
      final_insignificant_rois );

    // Build set of insignificant ROI ranges for filtering
    std::vector<std::pair<double,double>> insignificant_roi_ranges;
    for( const size_t roi_idx : final_insignificant_rois )
    {
      const RelActCalcAuto::RoiRange &roi = solution.m_final_roi_ranges[roi_idx];
      insignificant_roi_ranges.emplace_back( roi.lower_energy, roi.upper_energy );
    }

    // Populate result, filtering out peaks from insignificant ROIs
    result.status = solution.m_status;
    result.error_message = solution.m_error_message;

    // Filter peaks - only include those NOT in insignificant ROIs
    result.fit_peaks.clear();
    if( PeakFitImprove::debug_printout && !insignificant_roi_ranges.empty() )
      cout << "Peak filtering by ROI significance:" << endl;

    for( const PeakDef &peak : solution.m_peaks_without_back_sub )
    {
      const double mean = peak.mean();
      const double peak_roi_lower = peak.continuum()->lowerEnergy();
      const double peak_roi_upper = peak.continuum()->upperEnergy();
      bool in_insignificant_roi = false;

      for( const std::pair<double,double> &roi_range : insignificant_roi_ranges )
      {
        // Match if peak's ROI bounds are within 1 keV of the insignificant ROI bounds
        if( (fabs( peak_roi_lower - roi_range.first ) < 1.0)
           && (fabs( peak_roi_upper - roi_range.second ) < 1.0) )
        {
          in_insignificant_roi = true;
          break;
        }
      }

      if( in_insignificant_roi )
      {
        if( PeakFitImprove::debug_printout )
        {
          cout << "  Filtered (insignificant ROI [" << peak.continuum()->lowerEnergy() << ", " << peak.continuum()->upperEnergy()
               << "] keV): peak at " << mean << " keV, area = " << peak.peakArea() << endl;
        }
      }else
      {
        bool mean_in_roi = false;
        for( size_t roi_index = 0; !mean_in_roi && (roi_index < solution.m_final_roi_ranges_in_spectrum_cal.size()); ++roi_index )
        {
          const auto pos = std::find( begin(final_insignificant_rois), end(final_insignificant_rois), roi_index );
          if( pos != end(final_insignificant_rois) )
            continue;
          const RelActCalcAuto::RoiRange &roi = solution.m_final_roi_ranges_in_spectrum_cal[roi_index];
          mean_in_roi = ((mean >= roi.lower_energy) && (mean <= roi.upper_energy));
        }

        if( PeakFitImprove::debug_printout && !insignificant_roi_ranges.empty() )
        {
          cout << "  Kept (significant ROI [" << peak.continuum()->lowerEnergy() << ", " << peak.continuum()->upperEnergy()
               << "] keV): peak at " << mean << " keV, area = " << peak.peakArea() << (mean_in_roi ? " (was in ROI)" : " (skipping peak, not in a ROI)") << endl;
        }
        
        if( mean_in_roi )
          result.fit_peaks.push_back( peak );
      }
    }

    if( !insignificant_roi_ranges.empty() )
    {
      const size_t num_filtered = solution.m_peaks_without_back_sub.size() - result.fit_peaks.size();
      cout << "Filtered out " << num_filtered << " peaks from "
           << insignificant_roi_ranges.size() << " ROIs without significant chi2 improvement" << endl;
    }

    result.solution = std::move( solution );

    // Combine overlapping peaks within ROIs
    // First, preserve the uncombined peaks, then create combined version
    result.uncombined_fit_peaks = result.fit_peaks;
    result.fit_peaks = combine_overlapping_peaks_in_rois( result.uncombined_fit_peaks );

    if( PeakFitImprove::debug_printout && (result.fit_peaks.size() != result.uncombined_fit_peaks.size()) )
    {
      cout << "Combined " << result.uncombined_fit_peaks.size() << " peaks into "
           << result.fit_peaks.size() << " peaks" << endl;
    }

    // Compute observable peaks - peaks that users can visually see in the spectrum
    result.observable_peaks = compute_observable_peaks( result.fit_peaks, foreground, config );

    if( PeakFitImprove::debug_printout && (result.observable_peaks.size() != result.fit_peaks.size()) )
    {
      cout << "Observable peaks: " << result.observable_peaks.size() << " of "
           << result.fit_peaks.size() << " fit_peaks" << endl;
    }

    // Developer check: Look for duplicate peaks (same mean, different ROI)
#if( PERFORM_DEVELOPER_CHECKS )
    for( size_t i = 0; i < result.observable_peaks.size(); ++i )
    {
      for( size_t j = i + 1; j < result.observable_peaks.size(); ++j )
      {
        const PeakDef &peak_i = result.observable_peaks[i];
        const PeakDef &peak_j = result.observable_peaks[j];
        const double mean_diff = std::fabs( peak_i.mean() - peak_j.mean() );
        if( mean_diff < 0.5 )  // Same energy within 0.5 keV
        {
          const bool same_continuum = (peak_i.continuum() == peak_j.continuum());
          if( !same_continuum )
          {
            cerr << "WARNING: Duplicate peaks at " << peak_i.mean() << " keV with different ROIs: "
                 << "[" << peak_i.continuum()->lowerEnergy() << ", " << peak_i.continuum()->upperEnergy() << "] vs "
                 << "[" << peak_j.continuum()->lowerEnergy() << ", " << peak_j.continuum()->upperEnergy() << "]" << endl;
          }
        }
      }
    }
#endif

  }catch( const std::exception &e )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    result.error_message = e.what();
  }

  return result;
}//fit_peaks_for_nuclide_relactauto


/** Estimate initial ROIs when no auto-search peaks matched source nuclides.

 Creates ROIs directly from the highest BR*efficiency gamma lines for each source.
 For each source, selects up to 4 gammas with highest BR*efficiency scores.
 Each gamma defines an ROI as gamma_energy ± 2.5*FWHM.

 ROIs are merged if they overlap AND the combined width doesn't exceed
 config.auto_rel_eff_sol_max_fwhm FWHMs. If ROIs overlap but can't merge due to
 width constraint, they are split to eliminate overlap (critical requirement).

 All ROIs use Linear continuum and Fixed range limits.

 @param sources Vector of source nuclides/elements/reactions to consider
 @param drf Detector response function (may be nullptr, will use generic)
 @param isHPGe True for HPGe detector, false for NaI
 @param fwhmFnctnlForm FWHM functional form
 @param fwhm_coefficients FWHM coefficients
 @param lower_fwhm_energy Lower energy bound for FWHM validity
 @param upper_fwhm_energy Upper energy bound for FWHM validity
 @param min_valid_energy Minimum valid energy for spectrum
 @param max_valid_energy Maximum valid energy for spectrum
 @param config Configuration with merging thresholds
 @return Vector of ROI ranges, empty if no valid ROIs created
 */
std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_without_peaks(
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const bool isHPGe,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const double min_valid_energy,
  const double max_valid_energy,
  const PeakFitForNuclideConfig &config )
{
  // Step 1: Get or create valid DRF (use generic if nullptr)
  std::shared_ptr<const DetectorPeakResponse> drf_to_use = drf;
  if( !drf_to_use || !drf_to_use->isValid() )
  {
    drf_to_use = isHPGe
        ? DetectorPeakResponse::getGenericHPGeDetector()
        : DetectorPeakResponse::getGenericNaIDetector();
  }

  if( !drf_to_use || !drf_to_use->isValid() )
    return {};

  // Step 2: Collect top gammas per source
  // Data structure to hold gamma info
  struct GammaInfo
  {
    double energy;
    double br_times_eff;
    RelActCalcAuto::SrcVariant source;  // For debugging
  };

  std::vector<GammaInfo> selected_gammas;

  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    if( RelActCalcAuto::is_null( src ) )
      continue;

    // Get source age and photons
    const double age = get_source_age( src, -1.0 );
    const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src, 1.0, age );

    // Compute BR*eff scores for valid gammas
    std::vector<GammaInfo> candidates;
    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      if( photon.energy < min_valid_energy || photon.energy > max_valid_energy )
        continue;

      const double br = photon.numPerSecond;  // BR since we used unit activity
      const double eff = drf_to_use->intrinsicEfficiency( static_cast<float>(photon.energy) );
      const double score = br * eff;

      if( score > 0.0 )
        candidates.push_back( {photon.energy, score, src} );
    }

    // Sort by score (descending) and take top 4
    std::sort( candidates.begin(), candidates.end(),
      [](const GammaInfo &a, const GammaInfo &b) { return a.br_times_eff > b.br_times_eff; } );

    const size_t num_to_take = std::min( candidates.size(), size_t(4) );
    for( size_t i = 0; i < num_to_take; ++i )
      selected_gammas.push_back( candidates[i] );

    // Debug output
    if( PeakFitImprove::debug_printout )
    {
      cerr << "estimate_initial_rois_without_peaks: Source "
           << RelActCalcAuto::to_name( src ) << " - selected " << num_to_take << " gammas" << endl;
      for( size_t i = 0; i < num_to_take; ++i )
      {
        cerr << "  " << candidates[i].energy << " keV, BR*eff=" << candidates[i].br_times_eff << endl;
      }
    }
  }

  if( selected_gammas.empty() )
  {
    if( PeakFitImprove::debug_printout )
      cerr << "estimate_initial_rois_without_peaks: No valid gammas found" << endl;
    return {};
  }

  // Step 3: Create initial ROIs
  struct InitialRoi
  {
    RelActCalcAuto::RoiRange roi;
    double center_energy;
    double fwhm;
  };

  std::vector<InitialRoi> initial_rois;

  for( const GammaInfo &gamma : selected_gammas )
  {
    // Compute FWHM
    const double fwhm = DetectorPeakResponse::peakResolutionFWHM(
      static_cast<float>(gamma.energy), fwhmFnctnlForm, fwhm_coefficients );

    // Validate FWHM
    if( !std::isfinite( fwhm ) || fwhm <= 0.0 )
    {
      if( PeakFitImprove::debug_printout )
        cerr << "Warning: Invalid FWHM at " << gamma.energy << " keV, skipping" << endl;
      continue;
    }

    // Create ROI: energy ± 2.5 FWHM, clamped to valid range
    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = std::max( min_valid_energy, gamma.energy - 2.5 * fwhm );
    roi.upper_energy = std::min( max_valid_energy, gamma.energy + 2.5 * fwhm );
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;

    initial_rois.push_back( {roi, gamma.energy, fwhm} );
  }

  if( initial_rois.empty() )
  {
    if( PeakFitImprove::debug_printout )
      cerr << "estimate_initial_rois_without_peaks: No valid ROIs created" << endl;
    return {};
  }

  // Sort by lower_energy for merging
  std::sort( initial_rois.begin(), initial_rois.end(),
    [](const InitialRoi &a, const InitialRoi &b)
    {
      return a.roi.lower_energy < b.roi.lower_energy;
    } );

  // Step 4: Merge/split overlapping ROIs
  // Critical: ensure NO overlapping ROIs in final result
  std::vector<RelActCalcAuto::RoiRange> merged_rois;
  std::vector<std::vector<double>> merged_centers;  // Track center energies for each merged ROI
  std::vector<double> merged_fwhms;  // Track FWHM for validation

  for( const InitialRoi &current : initial_rois )
  {
    if( merged_rois.empty() )
    {
      merged_rois.push_back( current.roi );
      merged_centers.push_back( {current.center_energy} );
      merged_fwhms.push_back( current.fwhm );
      continue;
    }

    RelActCalcAuto::RoiRange &last = merged_rois.back();
    std::vector<double> &last_centers = merged_centers.back();
    const double last_fwhm = merged_fwhms.back();

    // Check if ROIs overlap
    const bool overlaps = (current.roi.lower_energy < last.upper_energy);

    if( !overlaps )
    {
      // No overlap - add new ROI
      merged_rois.push_back( current.roi );
      merged_centers.push_back( {current.center_energy} );
      merged_fwhms.push_back( current.fwhm );
      continue;
    }

    // ROIs overlap - check width constraint
    const double combined_upper = std::max( last.upper_energy, current.roi.upper_energy );
    const double combined_width = combined_upper - last.lower_energy;
    const double mid_energy = 0.5 * (last.lower_energy + combined_upper);
    const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
      static_cast<float>(mid_energy), fwhmFnctnlForm, fwhm_coefficients );

    const bool width_ok = (combined_width <= config.auto_rel_eff_sol_max_fwhm * mid_fwhm);

    if( width_ok )
    {
      // MERGE: Extend last ROI to encompass both
      last.upper_energy = combined_upper;
      last_centers.push_back( current.center_energy );
      // Update FWHM to use the average
      merged_fwhms.back() = 0.5 * (last_fwhm + current.fwhm);

      if( PeakFitImprove::debug_printout )
      {
        cerr << "Merged overlapping ROI [" << current.roi.lower_energy << ", "
             << current.roi.upper_energy << "] into [" << last.lower_energy
             << ", " << last.upper_energy << "]" << endl;
      }
    }
    else
    {
      // SPLIT: Width constraint prevents merge - give each ROI half the overlap
      const double overlap = last.upper_energy - current.roi.lower_energy;
      const double half_overlap = overlap / 2.0;

      // Adjust boundaries to split overlap
      const double original_last_upper = last.upper_energy;
      last.upper_energy = last.upper_energy - half_overlap;
      RelActCalcAuto::RoiRange adjusted_current = current.roi;
      adjusted_current.lower_energy = current.roi.lower_energy + half_overlap;

      // Validate adjusted last ROI still contains all its center energies and is wide enough
      const double last_width = last.upper_energy - last.lower_energy;
      bool last_valid = (last_width >= last_fwhm);
      for( const double center : last_centers )
      {
        if( center < last.lower_energy || center > last.upper_energy )
        {
          last_valid = false;
          break;
        }
      }

      // Validate adjusted current ROI is still useful
      const double adjusted_width = adjusted_current.upper_energy - adjusted_current.lower_energy;
      const bool current_contains_center = (current.center_energy >= adjusted_current.lower_energy)
                                        && (current.center_energy <= adjusted_current.upper_energy);
      const bool current_wide_enough = (adjusted_width >= current.fwhm);
      const bool current_valid = current_contains_center && current_wide_enough;

      if( !last_valid )
      {
        // Last ROI is no longer valid after split - remove it and add current as-is
        merged_rois.pop_back();
        merged_centers.pop_back();
        merged_fwhms.pop_back();

        merged_rois.push_back( current.roi );
        merged_centers.push_back( {current.center_energy} );
        merged_fwhms.push_back( current.fwhm );

        if( PeakFitImprove::debug_printout )
        {
          cerr << "Removed last ROI (invalid after split), kept current ROI at "
               << current.center_energy << " keV" << endl;
        }
      }
      else if( current_valid )
      {
        // Both ROIs valid - add the split current ROI
        merged_rois.push_back( adjusted_current );
        merged_centers.push_back( {current.center_energy} );
        merged_fwhms.push_back( current.fwhm );

        if( PeakFitImprove::debug_printout )
        {
          cerr << "Split overlapping ROIs (width constraint): overlap=" << overlap
               << " keV, last=[" << last.lower_energy << ", " << last.upper_energy
               << "], current=[" << adjusted_current.lower_energy << ", "
               << adjusted_current.upper_energy << "]" << endl;
        }
      }
      else
      {
        // Current ROI invalid after split, but last is valid - just keep last as adjusted
        if( PeakFitImprove::debug_printout )
        {
          cerr << "Skipping adjusted current ROI at " << current.center_energy << " keV: ";
          if( !current_contains_center )
            cerr << "doesn't contain source energy";
          else
            cerr << "too narrow (" << adjusted_width << " keV < " << current.fwhm << " keV FWHM)";
          cerr << endl;
        }
      }
    }
  }

  // Validate no overlaps (developer check)
#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < merged_rois.size(); ++i )
  {
    assert( merged_rois[i].lower_energy >= merged_rois[i-1].upper_energy );
  }
#endif

  if( PeakFitImprove::debug_printout )
  {
    cerr << "estimate_initial_rois_without_peaks: Created " << merged_rois.size()
         << " final ROIs" << endl;
  }

  return merged_rois;
}//estimate_initial_rois_without_peaks


/** Fallback function to estimate initial ROIs when RelActManual fails.

 Uses detector efficiency and spectrum integration to estimate activities for each source nuclide,
 then creates ROIs using the same clustering function as the successful RelActManual case.

 For each nuclide:
 1. Find the gamma with largest expected yield (branching_ratio * detector_efficiency)
 2. If an auto-fit peak matches this gamma (within 0.5*FWHM), use its area to estimate activity
 3. Otherwise, integrate spectrum counts around the gamma energy and use 1/4 as estimated peak area
 4. Use these activities with cluster_gammas_to_rois() to create ROIs
 */
std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_fallback(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const bool isHPGe,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const GammaClusteringSettings &settings )
{
  // Step 1: Get or create valid DRF
  std::shared_ptr<const DetectorPeakResponse> drf_to_use = drf;
  if( !drf_to_use || !drf_to_use->isValid() )
  {
    drf_to_use = isHPGe
        ? DetectorPeakResponse::getGenericHPGeDetector()
        : DetectorPeakResponse::getGenericNaIDetector();
  }

  if( !drf_to_use || !drf_to_use->isValid() )
    return {};

  const double live_time = foreground->live_time();
  if( live_time <= 0.0 )
    return {};

  // Find valid energy range based on contiguous channels with data
  const auto [min_valid_energy, max_valid_energy] = find_valid_energy_range( foreground );

  // Step 2: Estimate activity for each source
  std::vector<std::pair<RelActCalcAuto::SrcVariant, double>> sources_and_acts;

  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    if( RelActCalcAuto::is_null( src ) )
      continue;

    // Get gamma lines at default age
    const double age = get_source_age( src, -1.0 );
    const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src, 1.0, age );

    // Find gamma with largest expected yield (br * efficiency)
    double best_yield = 0.0;
    double best_energy = 0.0;
    double best_br = 0.0;
    double best_eff = 0.0;

    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      if( photon.energy < min_valid_energy || photon.energy > max_valid_energy )
        continue;

      const double br = photon.numPerSecond;  // BR since we used unit activity
      const double eff = drf_to_use->intrinsicEfficiency( static_cast<float>(photon.energy) );
      const double yield = br * eff;

      if( yield > best_yield )
      {
        best_yield = yield;
        best_energy = photon.energy;
        best_br = br;
        best_eff = eff;
      }
    }

    if( best_yield <= 0.0 || best_br <= 0.0 || best_eff <= 0.0 )
      continue;

    // Determine energy tolerance for peak matching (0.5 * FWHM)
    const double fwhm_at_energy = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(best_energy), fwhmFnctnlForm, fwhm_coefficients );
    const double energy_tolerance = 0.5 * fwhm_at_energy;

    // Try to find matching auto-fit peak
    std::shared_ptr<const PeakDef> matched_peak = nullptr;
    double min_distance = std::numeric_limits<double>::max();

    for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
    {
      if( !peak || !peak->gausPeak() )
        continue;

      const double distance = std::abs( peak->mean() - best_energy );
      if( distance < energy_tolerance && distance < min_distance )
      {
        min_distance = distance;
        matched_peak = peak;
      }
    }

    // Estimate activity
    double estimated_activity = 0.0;

    if( matched_peak )
    {
      // Use peak area: activity = peak_area / (br * eff) - we wont divide by live_time, to be consistent with `cluster_gammas_to_rois(...)` convention
      const double peak_area = matched_peak->peakArea();
      estimated_activity = peak_area / (best_br * best_eff);

      cout << "Fallback: " << RelActCalcAuto::to_name( src ) << " matched peak at " << matched_peak->mean()
           << " keV (gamma at " << best_energy << " keV), area=" << peak_area
           << ", estimated activity=" << estimated_activity << endl;
    }
    else
    {
      // Integrate spectrum ±0.5 FWHM, use 1/4 as estimated peak area
      const double integration_half_width = 0.5 * fwhm_at_energy;
      const float lower_e = static_cast<float>(best_energy - integration_half_width);
      const float upper_e = static_cast<float>(best_energy + integration_half_width);

      const double total_counts = foreground->gamma_integral( lower_e, upper_e );
      const double estimated_peak_area = total_counts / 4.0;

      estimated_activity = estimated_peak_area / (best_br * best_eff * live_time);

      cout << "Fallback: " << RelActCalcAuto::to_name( src ) << " no matching peak for gamma at " << best_energy
           << " keV, integrated counts=" << total_counts << ", est. peak area=" << estimated_peak_area
           << ", estimated activity=" << estimated_activity << endl;
    }

    if( estimated_activity > 0.0 )
      sources_and_acts.emplace_back( src, estimated_activity );
  }

  if( sources_and_acts.empty() )
  {
    cerr << "Fallback: Could not estimate activity for any source" << endl;
    return {};
  }

  // Step 3: Create relative efficiency function from DRF
  const auto fallback_rel_eff = [&drf_to_use]( double energy ) -> double {
    return drf_to_use->intrinsicEfficiency( static_cast<float>(energy) );
  };

  // Step 4: Call cluster_gammas_to_rois with estimated activities
  return cluster_gammas_to_rois(
    fallback_rel_eff,
    sources_and_acts,
    foreground,
    fwhmFnctnlForm,
    fwhm_coefficients,
    lower_fwhm_energy,
    upper_fwhm_energy,
    min_valid_energy,
    max_valid_energy,
    settings
  );
}//estimate_initial_rois_fallback


/** Estimates initial ROIs using RelActManual with multiple fallback strategies.

 This function attempts to estimate initial ROIs by:
 1. Converting auto_search_peaks to RelActManual format
 2. Matching peaks to source nuclides
 3. If no peaks matched -> falls back to estimate_initial_rois_without_peaks()
 4. Fitting a relative efficiency curve using RelActManual with matched peaks
 5. Extracting relative activities for each source from the solution
 6. Clustering gammas into ROIs using the relative efficiency and activities
 7. If RelActManual fails -> falls back to estimate_initial_rois_fallback()

 @param auto_search_peaks Peaks found by automatic search
 @param foreground Foreground measurement
 @param sources Source variants to fit
 @param drf Detector response function (may be nullptr)
 @param isHPGe Whether detector is HPGe (vs NaI)
 @param fwhmFnctnlForm FWHM functional form
 @param fwhm_coefficients FWHM coefficients
 @param lower_fwhm_energy Lower energy bound for FWHM validity
 @param upper_fwhm_energy Upper energy bound for FWHM validity
 @param min_valid_energy Minimum valid energy for clustering
 @param max_valid_energy Maximum valid energy for clustering
 @param manual_settings Clustering settings for manual mode
 @param config Peak fitting configuration
 @param fallback_warning Output parameter for warning message if fallback is used

 @return Vector of initial ROI ranges (never empty - will use fallback if needed)
 */
std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_using_relactmanual(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const bool isHPGe,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const double min_valid_energy,
  const double max_valid_energy,
  const GammaClusteringSettings &manual_settings,
  const PeakFitForNuclideConfig &config,
  std::string &fallback_warning )
{
  using namespace std;

  vector<RelActCalcAuto::RoiRange> initial_rois;

  // Step 1: Convert auto_search_peaks to RelActManual format
  vector<RelActCalcManual::GenericPeakInfo> rel_act_manual_peaks;

  for( const shared_ptr<const PeakDef> &peak : auto_search_peaks )
  {
    assert( peak && peak->gausPeak() );
    if( !peak || !peak->gausPeak() )
      continue;

    RelActCalcManual::GenericPeakInfo peak_info;
    peak_info.m_energy = peak_info.m_mean = peak->mean();
    peak_info.m_fwhm = peak->fwhm();
    peak_info.m_counts = peak->amplitude();
    peak_info.m_counts_uncert = peak->amplitudeUncert();
    peak_info.m_base_rel_eff_uncert = config.rel_eff_manual_base_rel_eff_uncert;

    rel_act_manual_peaks.push_back( peak_info );
  }

  // Step 2: Match peaks to source nuclides
  vector<RelActCalcManual::PeakCsvInput::NucAndAge> rel_act_manual_srcs;
  for( const RelActCalcAuto::SrcVariant &src : sources )
    rel_act_manual_srcs.emplace_back( RelActCalcAuto::to_name( src ), -1.0, false );

  const std::vector<std::pair<float,float>> energy_ranges{};
  const std::vector<float> excluded_peak_energies{};
  const float real_time = foreground->real_time();

  const RelActCalcManual::PeakCsvInput::NucMatchResults peak_match_results
   = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( rel_act_manual_peaks,
                                                          RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay,
                                                          energy_ranges,
                                                          rel_act_manual_srcs,
                                                          config.initial_nuc_match_cluster_num_sigma, excluded_peak_energies,
                                                          real_time );

  vector<RelActCalcManual::GenericPeakInfo> peaks_matched = peak_match_results.peaks_matched;
  std::sort( begin(peaks_matched), end(peaks_matched),
    []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ){
      return lhs.m_energy < rhs.m_energy;
  });

  if( PeakFitImprove::debug_printout )
  {
    const vector<string> &used_isotopes = peak_match_results.used_isotopes;
    const vector<string> &unused_isotopes = peak_match_results.unused_isotopes;
    if( unused_isotopes.empty() )
      cout << "Matched up all source nuclides to initial peak fit" << endl;
    cout << "Failed to match up nuclides: {";
    for( const string &nuc : unused_isotopes )
      cout << nuc << ", ";
    cout << "} to initial auto-fit peaks" << endl;
    if( !used_isotopes.empty() )
    {
      cout << "Matched up nuclides: {";
      for( const string &nuc : used_isotopes )
        cout << nuc << ", ";
      cout << "} to initial auto-fit peaks to a total of " << peaks_matched.size() << " peaks." << endl;

      if( !peaks_matched.empty() )
      {
        cout << "Matched peak energies (keV): {";
        for( size_t i = 0; i < peaks_matched.size(); ++i )
        {
          cout << peaks_matched[i].m_energy;
          if( i < peaks_matched.size() - 1 )
            cout << ", ";
        }
        cout << "}" << endl;
      }
    }
  }

  // If no matched peaks, fall back to estimate_initial_rois_without_peaks
  if( peaks_matched.empty() )
  {
    return estimate_initial_rois_without_peaks(
      sources, drf, isHPGe,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      min_valid_energy, max_valid_energy, config );
  }

  // Step 3: Configure RelActManual input based on number of matched peaks
  RelActCalcManual::RelEffInput manual_input;
  manual_input.peaks = peaks_matched;

  manual_input.eqn_order = 0;
  manual_input.use_ceres_to_fit_eqn = false;
  manual_input.phys_model_use_hoerl = false;

  if( peaks_matched.size() == 1 )
  {
    manual_input.eqn_form = config.initial_manual_relEff_1peak_form;

    // With only one peak, we can only have a zeroth order (constant) relative efficiency
    manual_input.eqn_order = 0;
    if( config.initial_manual_relEff_1peak_eqn_order != 0 )
    {
      fallback_warning = "Only 1 peak matched; forcing RelEff equation order to 0 (was "
                       + std::to_string(config.initial_manual_relEff_1peak_eqn_order) + ")";
    }

    if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      manual_input.use_ceres_to_fit_eqn = true;
      // With only one peak and physical model, we cannot fit any shielding parameters
      // The shielding vectors are left empty (set to empty below in the FramPhysicalModel block)
      if( fallback_warning.empty() )
        fallback_warning = "Only 1 peak matched with FramPhysicalModel; shielding parameters will not be fit";
      else
        fallback_warning += "; shielding parameters will not be fit";
    }
  }
  else if( peaks_matched.size() == 2 )
  {
    manual_input.eqn_order = config.initial_manual_relEff_2peak_eqn_order;
    manual_input.eqn_form = config.initial_manual_relEff_2peak_form;
    if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      manual_input.eqn_order = 0;
      manual_input.use_ceres_to_fit_eqn = true;
    }
  }
  else if( peaks_matched.size() == 3 )
  {
    manual_input.eqn_order = config.initial_manual_relEff_3peak_eqn_order;
    manual_input.eqn_form = config.initial_manual_relEff_3peak_form;
    if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      manual_input.eqn_order = 0;
      manual_input.use_ceres_to_fit_eqn = true;
    }
  }
  else if( peaks_matched.size() == 4 )
  {
    manual_input.eqn_order = config.initial_manual_relEff_4peak_eqn_order;
    manual_input.eqn_form = config.initial_manual_relEff_4peak_form;
    if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      manual_input.eqn_order = 0;
      manual_input.use_ceres_to_fit_eqn = true;
      manual_input.phys_model_use_hoerl = config.initial_manual_relEff_4peak_physical_use_hoerl;
    }
  }
  else
  {
    assert( peaks_matched.size() > 4 );
    manual_input.eqn_order = config.initial_manual_relEff_many_peak_eqn_order;
    manual_input.eqn_form = config.initial_manual_relEff_manypeak_form;
    if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      manual_input.eqn_order = 0;
      manual_input.use_ceres_to_fit_eqn = true;
      manual_input.phys_model_use_hoerl = config.initial_manual_relEff_many_peak_physical_use_hoerl;
    }
  }

  if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    manual_input.phys_model_detector = drf;
    if( !manual_input.phys_model_detector )
    {
      if( isHPGe )
        manual_input.phys_model_detector = DetectorPeakResponse::getGenericHPGeDetector();
      else
        manual_input.phys_model_detector = DetectorPeakResponse::getGenericNaIDetector();
    }

    manual_input.phys_model_self_atten = shared_ptr<const RelActCalc::PhysicalModelShieldInput>{};
    manual_input.phys_model_external_attens = vector<shared_ptr<const RelActCalc::PhysicalModelShieldInput>>{};
  }

  // Step 2: Solve for relative efficiency with retry logic
  try
  {
    RelActCalcManual::RelEffSolution manual_solution
        = RelActCalcManual::solve_relative_efficiency( manual_input );

    double chi2_dof = manual_solution.m_chi2 / (std::max)(manual_solution.m_dof, 1);

    if( (manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success)
       || (chi2_dof > 20.0) )
    {
      if( PeakFitImprove::debug_printout )
      {
        cout << "Initial manual solution failed: status=";
        switch( manual_solution.m_status )
        {
          case RelActCalcManual::ManualSolutionStatus::NotInitialized:
            cout << "NotInitialized";
            break;
          case RelActCalcManual::ManualSolutionStatus::ErrorInitializing:
            cout << "ErrorInitializing";
            break;
          case RelActCalcManual::ManualSolutionStatus::ErrorFindingSolution:
            cout << "ErrorFindingSolution";
            break;
          case RelActCalcManual::ManualSolutionStatus::ErrorGettingSolution:
            cout << "ErrorGettingSolution";
            break;
          case RelActCalcManual::ManualSolutionStatus::Success:
            cout << "Success";
            break;
        }
        cout << ", form=" << RelActCalc::to_str( manual_solution.m_input.eqn_form )
             << ", order=" << manual_solution.m_input.eqn_order
             << ", chi2=" << manual_solution.m_chi2
             << ", dof=" << manual_solution.m_dof
             << ", chi2/dof=" << chi2_dof;
        if( !manual_solution.m_error_message.empty() )
          cout << ", error: " << manual_solution.m_error_message;
        cout << endl;
      }//if( PeakFitImprove::debug_printout )
      RelActCalcManual::RelEffInput retry_input = manual_input;
      retry_input.eqn_form = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      retry_input.phys_model_use_hoerl = (manual_input.peaks.size() > 3);
      retry_input.use_ceres_to_fit_eqn = true;
      retry_input.eqn_order = 0;
      if( isHPGe )
        retry_input.phys_model_detector = DetectorPeakResponse::getGenericHPGeDetector();
      else
        retry_input.phys_model_detector = DetectorPeakResponse::getGenericNaIDetector();

      if( auto_search_peaks.size() > 3 )
      {

      }

      RelActCalcManual::RelEffSolution retry_solution
        = RelActCalcManual::solve_relative_efficiency( retry_input );

      double retry_chi2_dof = retry_solution.m_chi2 / (std::max)(retry_solution.m_dof, 1);

      if( PeakFitImprove::debug_printout )
      {
        cout << "Retry manual solution fit comparison:" << endl;
        cout << "  Original: form=" << RelActCalc::to_str( manual_solution.m_input.eqn_form )
        << ", order=" << manual_solution.m_input.eqn_order
        << ", chi2/dof=" << manual_solution.m_chi2 << "/" << manual_solution.m_dof
        << "=" << chi2_dof
        << ", success=" << (manual_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success)
        << endl;
        if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
          cout << "  Original error: " << manual_solution.m_error_message << endl;
        cout << "  Retry:    form=" << RelActCalc::to_str( retry_solution.m_input.eqn_form )
        << ", order=" << retry_solution.m_input.eqn_order
        << ", chi2/dof=" << retry_solution.m_chi2 << "/" << retry_solution.m_dof
        << "=" << retry_chi2_dof
        << ", success=" << (retry_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success)
        << endl;
        if( retry_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
          cout << "  Retry error: " << retry_solution.m_error_message << endl;

        if( retry_solution.m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        {
          cout << "  Retry shielding:";
          if( retry_solution.m_input.phys_model_self_atten )
          {
            const auto &self = retry_solution.m_input.phys_model_self_atten;
            if( self->material )
              cout << " self-atten=" << self->material->name;
            else
              cout << " self-atten AN=" << self->atomic_number << " AD=" << self->areal_density;
          }
          else
          {
            cout << " no self-atten";
          }

          if( !retry_solution.m_input.phys_model_external_attens.empty() )
          {
            cout << ", external-atten=" << retry_solution.m_input.phys_model_external_attens.size();
            for( size_t i = 0; i < retry_solution.m_input.phys_model_external_attens.size(); ++i )
            {
              const auto &ext = retry_solution.m_input.phys_model_external_attens[i];
              if( ext->material )
                cout << "[" << i << "]=" << ext->material->name;
              else
                cout << "[" << i << "] AN=" << ext->atomic_number << " AD=" << ext->areal_density;
            }
          }
          else
          {
            cout << ", no external-atten";
          }
          cout << ", use_hoerl=" << (retry_solution.m_input.phys_model_use_hoerl ? "true" : "false") << endl;
        }
        cout << endl;
      }


      if( (retry_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success)
         && (retry_chi2_dof < chi2_dof) )
      {
        if( PeakFitImprove::debug_printout )
          cout << "Will use retry solution!" << endl;
        chi2_dof = retry_chi2_dof;
        manual_solution = std::move(retry_solution);
      }else
      {
        if( PeakFitImprove::debug_printout )
          cout << "Will use original solution!" << endl;
      }
    }


    if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
      throw std::runtime_error( "Failed to fit initial RelActCalcManual::RelEffSolution: " + manual_solution.m_error_message );

    if( PeakFitImprove::debug_printout )
    {
      cout << "Successfully fitted initial RelActCalcManual::RelEffSolution: chi2/dof="
      << manual_solution.m_chi2 << "/" << manual_solution.m_dof << "="
      << manual_solution.m_chi2 / manual_solution.m_dof
      << " using " << peaks_matched.size() << " peaks"
      << endl;
      cout << endl;

      // Print expected vs observed peak areas and efficiency for each peak
      cout << "Peak areas and efficiency:" << endl;
      cout << "Energy(keV)  Observed    Expected    Efficiency  Obs/Exp" << endl;
      cout << "--------------------------------------------------------" << endl;
      for( const RelActCalcManual::GenericPeakInfo &peak : manual_solution.m_input.peaks )
      {
        const double efficiency = manual_solution.relative_efficiency( peak.m_energy );

        // Calculate expected source counts from relative activities and yields
        double expected_src_counts = 0.0;
        for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
        {
          const double rel_activity = manual_solution.relative_activity( line.m_isotope );
          expected_src_counts += rel_activity * line.m_yield;
        }

        const double expected_counts = expected_src_counts * efficiency;
        const double obs_over_exp = (expected_counts > 0.0) ? (peak.m_counts / expected_counts) : 0.0;

        cout << fixed;
        cout << setw(10) << setprecision(2) << peak.m_energy
             << setw(12) << setprecision(1) << peak.m_counts
             << setw(12) << setprecision(1) << expected_counts
             << setw(12) << setprecision(4) << efficiency
             << setw(10) << setprecision(3) << obs_over_exp
             << endl;
      }
      cout << endl;
    }//if( PeakFitImprove::debug_printout )

    // Step 3: Create rel_eff lambda from manual solution - handles extrapolation clamping
    const auto manual_rel_eff = [&manual_solution, &peaks_matched]( double energy ) -> double {
      // Extrapolation is terrible for rel-eff, so clamp to the lowest/highest peak energy
      if( energy < peaks_matched.front().m_energy )
        return manual_solution.relative_efficiency( peaks_matched.front().m_energy );
      else if( energy > peaks_matched.back().m_energy )
        return manual_solution.relative_efficiency( peaks_matched.back().m_energy );
      else
        return manual_solution.relative_efficiency( energy );
    };

    // Step 4: Collect sources and activities from the manual solution
    vector<pair<RelActCalcAuto::SrcVariant, double>> sources_and_acts;
    for( const RelActCalcManual::IsotopeRelativeActivity &rel_act : manual_solution.m_rel_activities )
    {
      const RelActCalcAuto::SrcVariant src = RelActCalcAuto::source_from_string( rel_act.m_isotope );
      if( RelActCalcAuto::is_null( src ) )
        throw std::logic_error( "Failed to get source from RelAct isotope '" + rel_act.m_isotope + "'" );
      sources_and_acts.emplace_back( src, manual_solution.relative_activity( rel_act.m_isotope ) );
    }

    // Step 5: Use the reusable clustering function to create ROIs
    initial_rois = cluster_gammas_to_rois( manual_rel_eff, sources_and_acts, foreground,
                                          fwhmFnctnlForm, fwhm_coefficients,
                                          lower_fwhm_energy, upper_fwhm_energy,
                                          min_valid_energy, max_valid_energy,
                                          manual_settings );

    if( PeakFitImprove::debug_printout )
    {
      cout << "Initial ROIs from RelActManual: ";
      for( const auto &roi : initial_rois )
        cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
      cout << endl;
    }
  }catch( std::exception &e )
  {
    cerr << "Error trying to fit initial manual rel-eff solution: " << e.what() << endl;
    cerr << "Using fallback activity estimation..." << endl;

    initial_rois = estimate_initial_rois_fallback(
      auto_search_peaks, foreground, sources, drf, isHPGe,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      manual_settings );

    fallback_warning = "RelActManual fitting failed (" + std::string(e.what())
                     + "); used simplified activity estimation fallback";

    cout << "Fallback ROIs: ";
    for( const RelActCalcAuto::RoiRange &roi : initial_rois )
      cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
    cout << endl;
  }

  return initial_rois;
}//estimate_initial_rois_using_relactmanual


/** Test function for estimate_initial_rois_without_peaks.

 Tests the function with Pu238, Pu239, and Pu241 sources and validates:
 - ROIs are created for each source
 - No ROIs overlap
 - All ROIs are at least 1 FWHM wide
 - All ROIs contain their source gamma energies

 @return true if all tests pass, false otherwise
 */
bool test_estimate_initial_rois_without_peaks()
{
  cout << "\n=== Testing estimate_initial_rois_without_peaks ===" << endl;

  bool all_tests_passed = true;

  try
  {
    // Setup test sources: Pu238, Pu239, Pu241
    std::vector<RelActCalcAuto::SrcVariant> test_sources;

    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
    {
      cerr << "TEST FAILED: Could not get decay database" << endl;
      return false;
    }

    const SandiaDecay::Nuclide * const pu238 = db->nuclide( "Pu238" );
    const SandiaDecay::Nuclide * const pu239 = db->nuclide( "Pu239" );
    const SandiaDecay::Nuclide * const pu241 = db->nuclide( "Pu241" );

    if( !pu238 || !pu239 || !pu241 )
    {
      cerr << "TEST FAILED: Could not find Pu238, Pu239, or Pu241 in database" << endl;
      return false;
    }

    test_sources.push_back( pu238 );
    test_sources.push_back( pu239 );
    test_sources.push_back( pu241 );

    // Setup test parameters
    const bool isHPGe = true;
    const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm
        = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;

    // Typical HPGe FWHM coefficients: FWHM = sqrt(a + b*E + c*E^2) where E is in keV
    const std::vector<float> fwhm_coefficients = {0.5f, 0.0f, 0.0f};  // ~0.7 keV FWHM constant
    const double lower_fwhm_energy = 0.0;
    const double upper_fwhm_energy = 3000.0;
    const double min_valid_energy = 50.0;
    const double max_valid_energy = 2000.0;

    // Create test config
    PeakFitForNuclideConfig config;
    config.auto_rel_eff_sol_max_fwhm = 12.0;

    // Get generic HPGe DRF
    std::shared_ptr<const DetectorPeakResponse> drf = DetectorPeakResponse::getGenericHPGeDetector();
    if( !drf || !drf->isValid() )
    {
      cerr << "TEST FAILED: Could not get generic HPGe detector" << endl;
      return false;
    }

    // Call the function under test
    const std::vector<RelActCalcAuto::RoiRange> rois = estimate_initial_rois_without_peaks(
      test_sources, drf, isHPGe,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      min_valid_energy, max_valid_energy, config );

    cout << "Created " << rois.size() << " ROIs" << endl;

    // Test 1: Check that we have ROIs
    if( rois.empty() )
    {
      cerr << "TEST FAILED: No ROIs created" << endl;
      all_tests_passed = false;
    }
    else
    {
      cout << "PASS: ROIs were created" << endl;
    }

    // Test 2: Check no overlaps
    for( size_t i = 1; i < rois.size(); ++i )
    {
      if( rois[i].lower_energy < rois[i-1].upper_energy )
      {
        cerr << "TEST FAILED: ROIs " << (i-1) << " and " << i << " overlap: "
             << "[" << rois[i-1].lower_energy << ", " << rois[i-1].upper_energy << "] "
             << "and [" << rois[i].lower_energy << ", " << rois[i].upper_energy << "]" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: No overlapping ROIs" << endl;

    // Test 3: Check that each ROI is at least 1 FWHM wide
    for( size_t i = 0; i < rois.size(); ++i )
    {
      const double roi_width = rois[i].upper_energy - rois[i].lower_energy;
      const double mid_energy = 0.5 * (rois[i].lower_energy + rois[i].upper_energy);
      const double fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(mid_energy), fwhmFnctnlForm, fwhm_coefficients );

      if( roi_width < fwhm )
      {
        cerr << "TEST FAILED: ROI " << i << " is narrower than 1 FWHM: "
             << roi_width << " keV < " << fwhm << " keV FWHM at " << mid_energy << " keV" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: All ROIs are at least 1 FWHM wide" << endl;

    // Test 4: Check that we have ROIs for each source
    // Collect expected gamma energies for each source (top 4 by BR*eff)
    std::map<const SandiaDecay::Nuclide*, std::vector<double>> expected_gammas;

    for( const RelActCalcAuto::SrcVariant &src : test_sources )
    {
      const double age = get_source_age( src, -1.0 );
      const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src, 1.0, age );

      struct GammaScore
      {
        double energy;
        double score;
      };
      std::vector<GammaScore> candidates;

      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        if( photon.energy < min_valid_energy || photon.energy > max_valid_energy )
          continue;

        const double br = photon.numPerSecond;
        const double eff = drf->intrinsicEfficiency( static_cast<float>(photon.energy) );
        const double score = br * eff;

        if( score > 0.0 )
          candidates.push_back( {photon.energy, score} );
      }

      // Sort and take top 4
      std::sort( candidates.begin(), candidates.end(),
        [](const GammaScore &a, const GammaScore &b) { return a.score > b.score; } );

      const size_t num_to_take = std::min( candidates.size(), size_t(4) );
      const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( src );

      for( size_t i = 0; i < num_to_take; ++i )
        expected_gammas[nuc].push_back( candidates[i].energy );
    }

    // Check that ROIs cover expected gammas
    for( const auto &[nuc, gammas] : expected_gammas )
    {
      cout << "Checking " << nuc->symbol << " gammas:" << endl;

      for( const double gamma_energy : gammas )
      {
        bool found = false;
        for( const RelActCalcAuto::RoiRange &roi : rois )
        {
          if( gamma_energy >= roi.lower_energy && gamma_energy <= roi.upper_energy )
          {
            found = true;
            cout << "  " << gamma_energy << " keV: FOUND in ROI [" << roi.lower_energy
                 << ", " << roi.upper_energy << "]" << endl;
            break;
          }
        }

        if( !found )
        {
          cerr << "TEST FAILED: " << nuc->symbol << " gamma at " << gamma_energy
               << " keV not found in any ROI" << endl;
          all_tests_passed = false;
        }
      }
    }

    // Test 5: Verify all ROIs have Linear continuum type
    for( size_t i = 0; i < rois.size(); ++i )
    {
      if( rois[i].continuum_type != PeakContinuum::OffsetType::Linear )
      {
        cerr << "TEST FAILED: ROI " << i << " does not have Linear continuum type" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: All ROIs have Linear continuum type" << endl;

    // Test 6: Verify all ROIs have Fixed range limits
    for( size_t i = 0; i < rois.size(); ++i )
    {
      if( rois[i].range_limits_type != RelActCalcAuto::RoiRange::RangeLimitsType::Fixed )
      {
        cerr << "TEST FAILED: ROI " << i << " does not have Fixed range limits" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: All ROIs have Fixed range limits" << endl;

  }
  catch( const std::exception &e )
  {
    cerr << "TEST FAILED: Exception thrown: " << e.what() << endl;
    all_tests_passed = false;
  }

  if( all_tests_passed )
    cout << "=== ALL TESTS PASSED ===" << endl;
  else
    cout << "=== SOME TESTS FAILED ===" << endl;

  cout << endl;

  return all_tests_passed;
}//test_estimate_initial_rois_without_peaks


PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const SpecUtils::Measurement> &long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const PeakFitForNuclideConfig &config,
  const bool isHPGe )
{
  PeakFitResult result;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "Failed to open SandiaDecayDataBase";
    return result;
  }

  // Validate sources
  if( sources.empty() )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "No sources provided";
    return result;
  }

  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    if( RelActCalcAuto::is_null( src ) )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Null source in sources vector";
      return result;
    }
  }

  // Use input DRF or create a copy we can modify
  std::shared_ptr<const DetectorPeakResponse> drf = drf_input;

  std::string fallback_warning;  // Set if we use the fallback activity estimation

  try
  {
    // Step 1: Determine FWHM functional form from auto-search peaks or DRF
    DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm = config.fwhm_functional_form;
    double lower_fwhm_energy = -1.0, upper_fwhm_energy = -1.0;
    vector<float> fwhm_coefficients, fwhm_uncerts;

    if( !drf || !drf->isValid() || !drf->hasResolutionInfo() || (auto_search_peaks.size() > 6) )
    {
      const int num_auto_peaks = static_cast<int>(auto_search_peaks.size());
      int sqrtEqnOrder = (std::min)( 6, num_auto_peaks / (1 + (num_auto_peaks > 3)) );
      if( auto_search_peaks.size() < 3 )
        sqrtEqnOrder = static_cast<int>( auto_search_peaks.size() );

      auto auto_search_peaks_dq
        = make_shared<const deque<shared_ptr<const PeakDef>>>( begin(auto_search_peaks), end(auto_search_peaks) );

      MakeDrfFit::performResolutionFit( auto_search_peaks_dq, fwhmFnctnlForm, sqrtEqnOrder, fwhm_coefficients, fwhm_uncerts );
      auto_search_peaks_dq = MakeDrfFit::removeOutlyingWidthPeaks( auto_search_peaks_dq, fwhmFnctnlForm, fwhm_coefficients );
      MakeDrfFit::performResolutionFit( auto_search_peaks_dq, fwhmFnctnlForm, sqrtEqnOrder, fwhm_coefficients, fwhm_uncerts );

      // Set energy range based on peaks used for FWHM fit
      if( !auto_search_peaks_dq->empty() )
      {
        lower_fwhm_energy = auto_search_peaks_dq->front()->mean();
        upper_fwhm_energy = auto_search_peaks_dq->back()->mean();
      }
    }
    else
    {
      fwhmFnctnlForm = drf->resolutionFcnType();
      fwhm_coefficients = drf->resolutionFcnCoefficients();

      // Get energy range from detector response function
      lower_fwhm_energy = drf->lowerEnergy();
      upper_fwhm_energy = drf->upperEnergy();
    }

    // Validate that we have valid FWHM coefficients
    if( fwhm_coefficients.empty() )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Failed to determine FWHM coefficients - unable to proceed with peak fitting";
      return result;
    }

    // Check that coefficients are finite
    for( size_t i = 0; i < fwhm_coefficients.size(); ++i )
    {
      if( !std::isfinite( fwhm_coefficients[i] ) )
      {
        result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
        result.error_message = "FWHM coefficient[" + std::to_string(i) + "] is not finite (value="
                               + std::to_string(fwhm_coefficients[i]) + ")";
        return result;
      }
    }

    if( PeakFitImprove::debug_printout )
    {
      cout << "FWHM function form: " << static_cast<int>(fwhmFnctnlForm) << ", coefficients: [";
      for( size_t i = 0; i < fwhm_coefficients.size(); ++i )
        cout << (i > 0 ? ", " : "") << fwhm_coefficients[i];
      cout << "]" << endl;
    }

    // Find valid energy range based on contiguous channels with data
    const auto [min_valid_energy, max_valid_energy] = find_valid_energy_range( foreground );

    // Determine energy range for gamma lines
    double highest_energy_gamma = 0.0, lowest_energy_gamma = std::numeric_limits<double>::max();

    vector<RelActCalcAuto::NucInputInfo> base_nuclides;
    for( const RelActCalcAuto::SrcVariant &src : sources )
    {
      RelActCalcAuto::NucInputInfo nuc_info;
      nuc_info.source = src;
      nuc_info.age = get_source_age( src, -1.0 );
      nuc_info.fit_age = false;
      base_nuclides.push_back(nuc_info);

      const vector<SandiaDecay::EnergyRatePair> photons
          = get_source_photons( src, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits, nuc_info.age );
      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        highest_energy_gamma = (std::max)( highest_energy_gamma, photon.energy );
        lowest_energy_gamma = (std::min)( lowest_energy_gamma, photon.energy );
      }
    }

    lowest_energy_gamma = (std::max)( lowest_energy_gamma - (isHPGe ? 5 : 25), (double)foreground->gamma_energy_min() );
    highest_energy_gamma = (std::min)( highest_energy_gamma + (isHPGe ? 5 : 25), (double)foreground->gamma_energy_max() );

    // Step 2, 3 & 4: Estimate initial ROIs using RelActManual with multiple fallbacks
    // This function internally:
    // - Converts auto_search_peaks to RelActManual format and matches to sources
    // - Falls back to estimate_initial_rois_without_peaks() if no peaks match
    // - Fits relative efficiency curve and clusters gammas into ROIs
    // - Falls back to estimate_initial_rois_fallback() if RelActManual fails
    const GammaClusteringSettings manual_settings = config.get_manual_clustering_settings();

    const vector<RelActCalcAuto::RoiRange> initial_rois = estimate_initial_rois_using_relactmanual(
      auto_search_peaks, foreground, sources, drf, isHPGe,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      min_valid_energy, max_valid_energy, manual_settings,
      config, fallback_warning );

    // Call RelActAuto with initial_rois
    result = fit_peaks_for_nuclide_relactauto(
      auto_search_peaks, foreground, sources,
      initial_rois, long_background, drf, config,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy );

    /*
    // If RelActAuto failed completely, try a desperation Physical Model approach
    // This is a last-ditch effort to get some result
    if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      if( PeakFitImprove::debug_printout )
        std::cerr << "fit_peaks_for_nuclides: RelActAuto failed, trying desperation Physical Model..." << std::endl;

      // Create a new config with Physical Model settings
      PeakFitForNuclideConfig desperation_config = config;
      desperation_config.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      desperation_config.rel_eff_eqn_order = 0;
      desperation_config.phys_model_self_atten.clear();
      desperation_config.phys_model_external_atten.clear();

      // Apply external shielding if conditions are met
      if( should_use_desperation_shielding( config.desperation_phys_model_atomic_number, initial_rois ) )
      {
        try
        {
          std::shared_ptr<RelActCalc::PhysicalModelShieldInput> shield
            = create_desperation_shielding( config.desperation_phys_model_atomic_number,
                                            config.desperation_phys_model_areal_density_g_per_cm2 );
          desperation_config.phys_model_external_atten.push_back( shield );

          if( PeakFitImprove::debug_printout )
          {
            std::cerr << "fit_peaks_for_nuclides desperation: using external shielding with AN="
                      << config.desperation_phys_model_atomic_number
                      << ", starting AD=" << config.desperation_phys_model_areal_density_g_per_cm2
                      << " g/cm2" << std::endl;
          }
        }
        catch( const std::exception &e )
        {
          if( PeakFitImprove::debug_printout )
            std::cerr << "Failed to create desperation shielding: " << e.what() << std::endl;
          // Continue without shielding
        }
      }
      else
      {
        if( PeakFitImprove::debug_printout )
          std::cerr << "fit_peaks_for_nuclides desperation: not using external shielding" << std::endl;
      }

      desperation_config.phys_model_use_hoerl = (initial_rois.size() > 2);

      PeakFitResult desperation_result = fit_peaks_for_nuclide_relactauto(
        auto_search_peaks, foreground, sources,
        initial_rois, long_background, drf, desperation_config,
        fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy );

      // Use desperation result if it succeeded
      if( desperation_result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
      {
        if( PeakFitImprove::debug_printout )
          std::cerr << "fit_peaks_for_nuclides: desperation Physical Model succeeded!" << std::endl;

        desperation_result.warnings.push_back( "Original fit failed; used desperation Physical Model approach" );
        result = std::move( desperation_result );
      }
      else
      {
        if( PeakFitImprove::debug_printout )
          std::cerr << "fit_peaks_for_nuclides: desperation Physical Model also failed" << std::endl;
      }
    }
     */
  }catch( const std::exception &e )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    result.error_message = e.what();
  }

  // Propagate fallback warning if we used it
  if( !fallback_warning.empty() )
    result.warnings.push_back( fallback_warning );

  return result;
}//fit_peaks_for_nuclides


// Struct to hold local minimum info for ROI breakpoint selection
struct LocalMinimum {
  size_t channel;                  // Absolute channel number of minimum
  double synthetic_value;
  double depth_score;              // For tiebreaking (higher = better)
  double statistical_significance; // Primary criterion (lower = better breakpoint)
};


/** Build a synthetic spectrum from expected Gaussian peak shapes.

 Uses PeakDists::photopeak_function_integral() for consistency with other peak calculations.
 The synthetic spectrum uses the same binning as the actual data.

 @param gamma_energies Energies of gamma lines
 @param gamma_amplitudes Expected peak areas/amplitudes
 @param fwhm_at_energy Function returning FWHM at a given energy
 @param channel_energies Pointer to channel energy edges from data->channel_energies()->data()
 @param start_channel The starting channel for the synthetic spectrum
 @param num_channels Number of channels to compute
 @return Vector of synthetic counts for each channel
 */
vector<double> build_synthetic_spectrum(
  const vector<double> &gamma_energies,
  const vector<double> &gamma_amplitudes,
  const function<double(double)> &fwhm_at_energy,
  const float *channel_energies,
  size_t start_channel,
  size_t num_channels )
{
  assert( gamma_energies.size() == gamma_amplitudes.size() );

  // Zero-initialize the synthetic spectrum (same binning as data)
  vector<double> synthetic( num_channels, 0.0 );

  // Add each gamma line's Gaussian contribution using photopeak_function_integral
  for( size_t g = 0; g < gamma_energies.size(); ++g )
  {
    const double mean = gamma_energies[g];
    const double amplitude = gamma_amplitudes[g];
    const double sigma = fwhm_at_energy( mean ) / 2.35482;  // FWHM to sigma

    // photopeak_function_integral adds the peak contribution to synthetic[]
    // Pass pointer to channel energies starting at start_channel
    PeakDists::photopeak_function_integral<double>(
      mean,
      sigma,
      amplitude,
      PeakDef::SkewType::NoSkew,  // No skew for synthetic spectrum
      nullptr,                     // No skew parameters
      num_channels,
      &channel_energies[start_channel],
      synthetic.data()
    );
  }

  return synthetic;
}//build_synthetic_spectrum


/** Computes how significant the peak in the synthetic spectrum is compared to the noise level in the data.

 This metric checks whether the expected peak counts (from the synthetic spectrum)
 would be detectable above the statistical noise in the actual data.
 Significance = sum(synthetic) / sqrt(sum(data))
 This is analogous to signal / noise.

 @param synthetic Pre-computed synthetic spectrum, starting at start_channel
 @param start_channel The channel that synthetic[0] corresponds to
 @param check_start_ch First channel to check (absolute)
 @param check_end_ch One past last channel to check (absolute)
 @param data The spectrum measurement to compare against
 @return Significance metric: sum(synthetic) / sqrt(sum(data))
 */
double compute_significance_in_region(
  const vector<double> &synthetic,
  size_t start_channel,
  size_t check_start_ch,
  size_t check_end_ch,
  const shared_ptr<const SpecUtils::Measurement> &data )
{
  if( check_end_ch <= check_start_ch )
    return 0.0;  // Too few channels, return neutral

  // Ensure check range is within synthetic range
  if( check_start_ch < start_channel )
    check_start_ch = start_channel;
  if( check_end_ch > start_channel + synthetic.size() )
    check_end_ch = start_channel + synthetic.size();

  if( check_end_ch <= check_start_ch )
    return 0.0;

  // Sum the expected peak counts from synthetic and the data counts
  double sum_synthetic = 0.0;
  double sum_data = 0.0;

  for( size_t ch = check_start_ch; ch < check_end_ch; ++ch )
  {
    const size_t syn_idx = ch - start_channel;
    sum_synthetic += synthetic[syn_idx];
    sum_data += data->gamma_channel_content( ch );
  }

  // Significance = expected_signal / noise
  // where noise = sqrt(data_counts) (Poisson statistics)
  if( sum_data <= 0.0 )
    return 0.0;

  return sum_synthetic / std::sqrt( sum_data );
}//compute_significance_in_region


/** Checks if there's a statistically significant local maximum in the synthetic spectrum.

 Uses the pre-computed synthetic spectrum and compares to actual data.

 @param lower_channel Search range start (absolute channel)
 @param upper_channel Search range end (absolute channel)
 @param synthetic Pre-computed synthetic spectrum
 @param start_channel Channel that synthetic[0] corresponds to
 @param data Spectrum measurement to compare against
 @param channel_energies Pointer to channel energy edges
 @param fwhm_at_energy Function returning FWHM at a given energy
 @param check_fwhm_fraction Fraction of FWHM to use for check region
 @param significance_threshold Minimum significance to consider a peak significant
 @return true if a significant peak exists between the channels
 */
bool has_significant_peak_between(
  size_t lower_channel,
  size_t upper_channel,
  const vector<double> &synthetic,
  size_t start_channel,
  const shared_ptr<const SpecUtils::Measurement> &data,
  const float *channel_energies,
  const function<double(double)> &fwhm_at_energy,
  double check_fwhm_fraction,
  double significance_threshold )
{
  const size_t num_channels = synthetic.size();

  // Convert to indices within synthetic
  if( lower_channel < start_channel )
    lower_channel = start_channel;
  if( upper_channel > start_channel + num_channels )
    upper_channel = start_channel + num_channels;

  // Find local maxima in synthetic between lower_channel and upper_channel
  size_t num_local_max = 0;
  double max_significance = 0.0;

  for( size_t ch = lower_channel; ch < upper_channel; ++ch )
  {
    const size_t i = ch - start_channel;
    if( i == 0 || i >= num_channels - 1 )
      continue;

    // Check if this is a local maximum in the synthetic spectrum
    if( synthetic[i] > synthetic[i-1] && synthetic[i] > synthetic[i+1] )
    {
      num_local_max++;

      // Get energy at this channel for FWHM calculation
      const double ch_energy = 0.5 * (channel_energies[ch] + channel_energies[ch + 1]);
      const double fwhm = fwhm_at_energy( ch_energy );

      // Convert FWHM fraction to channel range for significance check
      const double half_width = check_fwhm_fraction * fwhm;
      const size_t check_start = data->find_gamma_channel( static_cast<float>(ch_energy - half_width) );
      const size_t check_end = data->find_gamma_channel( static_cast<float>(ch_energy + half_width) );

      const double significance = compute_significance_in_region(
        synthetic, start_channel, check_start, check_end, data );

      if( significance > max_significance )
        max_significance = significance;

      if( significance >= significance_threshold )
        return true;  // Found a significant peak
    }
  }

  // Debug output - only print occasionally to avoid overwhelming output
  static size_t call_count = 0;
  if( PeakFitImprove::debug_printout && (call_count++ % 10 == 0) )
  {
    const double lower_energy = channel_energies[lower_channel];
    const double upper_energy = channel_energies[upper_channel];
    cerr << "  has_significant_peak_between [" << lower_energy << ", " << upper_energy
         << "] keV: found " << num_local_max << " local maxima, max_sig=" << max_significance
         << ", threshold=" << significance_threshold;

    // Also show sample synthetic and data values
    double max_synth = 0.0;
    double data_at_max_synth = 0.0;
    for( size_t ch = lower_channel; ch < upper_channel; ++ch )
    {
      const size_t i = ch - start_channel;
      if( i >= num_channels )
        break;
      if( synthetic[i] > max_synth )
      {
        max_synth = synthetic[i];
        data_at_max_synth = data->gamma_channel_content( ch );
      }
    }
    cerr << ", max_synth=" << max_synth << ", data_at_max_synth=" << data_at_max_synth << endl;
  }

  return false;  // No significant peak found
}//has_significant_peak_between


/** Finds local minima in the synthetic spectrum and computes their significance.

 Computes both depth_score (tiebreaker) and statistical significance (primary criterion).
 Uses the pre-computed synthetic spectrum for all calculations.

 @param synthetic Pre-computed synthetic spectrum
 @param start_channel Channel that synthetic[0] corresponds to
 @param data Spectrum measurement to compare against
 @param channel_energies Pointer to channel energy edges
 @param fwhm_at_energy Function returning FWHM at a given energy
 @param check_fwhm_fraction Fraction of FWHM to use for check region
 @return Vector of LocalMinimum structs for each local minimum found
 */
vector<LocalMinimum> find_synthetic_minima(
  const vector<double> &synthetic,
  size_t start_channel,
  const shared_ptr<const SpecUtils::Measurement> &data,
  const float *channel_energies,
  const function<double(double)> &fwhm_at_energy,
  double check_fwhm_fraction )
{
  vector<LocalMinimum> minima;
  const size_t num_channels = synthetic.size();

  for( size_t i = 1; i < num_channels - 1; ++i )
  {
    // Check if this is a local minimum
    if( synthetic[i] < synthetic[i-1] && synthetic[i] < synthetic[i+1] )
    {
      const double min_value = synthetic[i];
      const size_t abs_channel = start_channel + i;

      // Find left maximum (scan left until we find a local max or boundary)
      double left_max = min_value;
      for( size_t j = i; j > 0; --j )
      {
        if( synthetic[j] > left_max )
          left_max = synthetic[j];
        if( j > 0 && synthetic[j] > synthetic[j-1] && synthetic[j] > synthetic[j+1] )
          break;  // Found left local maximum
      }

      // Find right maximum
      double right_max = min_value;
      for( size_t j = i; j < num_channels - 1; ++j )
      {
        if( synthetic[j] > right_max )
          right_max = synthetic[j];
        if( synthetic[j] > synthetic[j-1] && synthetic[j] > synthetic[j+1] )
          break;  // Found right local maximum
      }

      // Compute relative depth score (for tiebreaking)
      const double smaller_max = min( left_max, right_max );
      const double larger_max = max( left_max, right_max );
      const double depth = smaller_max - min_value;
      const double depth_score = (larger_max > 0) ? depth / larger_max : 0.0;

      // Get energy at this channel for FWHM-based check region
      const double ch_energy = 0.5 * (channel_energies[abs_channel] + channel_energies[abs_channel + 1]);
      const double fwhm = fwhm_at_energy( ch_energy );
      const double half_width = check_fwhm_fraction * fwhm;
      const size_t check_start = data->find_gamma_channel( static_cast<float>(ch_energy - half_width) );
      const size_t check_end = data->find_gamma_channel( static_cast<float>(ch_energy + half_width) );

      const double significance = compute_significance_in_region(
        synthetic, start_channel, check_start, check_end, data );

      LocalMinimum lm;
      lm.channel = abs_channel;
      lm.synthetic_value = min_value;
      lm.depth_score = depth_score;
      lm.statistical_significance = significance;
      minima.push_back( lm );
    }
  }

  return minima;
}//find_synthetic_minima


/** Determines if a step continuum should be used for a ROI based on left/right continuum comparison.

 The decision is made by:
 1. Finding the gamma line with the largest amplitude (largest BR * eff * activity) in the ROI
 2. Integrating data 1.5 FWHM away from this reference gamma on each side (0.25 FWHM width)
 3. Comparing left vs right using statistical significance

 @param cluster The clustered gamma info containing gamma energies and amplitudes
 @param foreground The spectrum measurement data
 @param fwhm_form The FWHM functional form
 @param fwhm_coefficients FWHM coefficients
 @param fwhm_lower_energy Lower energy bound for FWHM clamping (-1 if no limit)
 @param fwhm_upper_energy Upper energy bound for FWHM clamping (-1 if no limit)
 @param roi_lower Lower energy of the ROI
 @param roi_upper Upper energy of the ROI
 @param step_cont_left_right_nsigma Threshold for left/right difference (in sigma)
 @return true if step continuum should be used, false otherwise
 */
bool should_use_step_continuum(
  const ClusteredGammaInfo &cluster,
  const shared_ptr<const SpecUtils::Measurement> &foreground,
  const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
  const vector<float> &fwhm_coefficients,
  const double fwhm_lower_energy,
  const double fwhm_upper_energy,
  const double roi_lower,
  const double roi_upper,
  const double step_cont_left_right_nsigma )
{
  // Fixed integration parameters
  static constexpr double INTEGRATION_OFFSET_FWHM = 1.5;
  static constexpr double INTEGRATION_WIDTH_FWHM = 0.25;

  if( cluster.gamma_amplitudes.empty() || !foreground )
    return false;

  // Find the gamma with the largest amplitude (which corresponds to largest BR * eff * activity)
  const auto max_it = max_element( begin(cluster.gamma_amplitudes), end(cluster.gamma_amplitudes) );
  const size_t max_idx = static_cast<size_t>( distance( begin(cluster.gamma_amplitudes), max_it ) );
  const double ref_gamma_energy = cluster.gamma_energies[max_idx];

  // Check for valid FWHM range
  const bool have_fwhm_range = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0);
  double fwhm_eval_energy = ref_gamma_energy;
  if( have_fwhm_range )
    fwhm_eval_energy = std::max( fwhm_lower_energy, std::min( fwhm_upper_energy, ref_gamma_energy ) );

  const float fwhm = DetectorPeakResponse::peakResolutionFWHM(
      static_cast<float>(fwhm_eval_energy), fwhm_form, fwhm_coefficients );

  if( !std::isfinite(fwhm) || (fwhm <= 0.0f) )
    return false;

  // Calculate integration regions
  // Left side: 1.5 FWHM below reference gamma, with 0.25 FWHM width
  const double left_center = ref_gamma_energy - INTEGRATION_OFFSET_FWHM * fwhm;
  const double left_lower = left_center - 0.5 * INTEGRATION_WIDTH_FWHM * fwhm;
  const double left_upper = left_center + 0.5 * INTEGRATION_WIDTH_FWHM * fwhm;

  // Right side: 1.5 FWHM above reference gamma, with 0.25 FWHM width
  const double right_center = ref_gamma_energy + INTEGRATION_OFFSET_FWHM * fwhm;
  const double right_lower = right_center - 0.5 * INTEGRATION_WIDTH_FWHM * fwhm;
  const double right_upper = right_center + 0.5 * INTEGRATION_WIDTH_FWHM * fwhm;

  // Ensure integration regions are within the ROI
  if( (left_lower < roi_lower) || (right_upper > roi_upper) )
    return false;  // Integration regions extend beyond ROI, cannot determine

  // Sum counts in each region
  const double left_sum = foreground->gamma_integral( static_cast<float>(left_lower),
                                                       static_cast<float>(left_upper) );
  const double right_sum = foreground->gamma_integral( static_cast<float>(right_lower),
                                                        static_cast<float>(right_upper) );

  // Correct Poisson uncertainty for difference: sigma = sqrt(left + right)
  const double combined_uncert = std::sqrt( left_sum + right_sum );

  if( combined_uncert <= 0.0 )
    return false;

  const double nsigma = (left_sum - right_sum) / combined_uncert;

  if( PeakFitImprove::debug_printout )
  {
    cerr << "should_use_step_continuum: ref_gamma=" << ref_gamma_energy << " keV, "
         << "left_sum=" << left_sum << ", right_sum=" << right_sum
         << ", nsigma=" << nsigma
         << " (threshold=" << step_cont_left_right_nsigma << ")" << endl;
  }

  // If left side is significantly higher than right side, suggest step continuum
  return (nsigma >= step_cont_left_right_nsigma);
}//should_use_step_continuum



/** Clusters gamma lines and creates ROIs based on a relative efficiency function.

 This function takes gamma lines from the specified nuclides, estimates their relative
 amplitudes using the provided relative efficiency function, and clusters them into
 ROIs (Regions of Interest) for fitting.

 The clustering algorithm:
 1. Calculates expected counts for each gamma line using rel_eff_fcn
 2. Orders gamma lines by expected counts (highest first)
 3. For each gamma line, creates a cluster covering +/- cluster_num_sigma * sigma
 4. Merges overlapping clusters
 5. Breaks up clusters that are too wide (> max_fwhm_width FWHMs)
 6. Filters clusters that are too small or have insignificant expected signal
 7. Creates ROIs with appropriate bounds and continuum types

 @param rel_eff_fcn Lambda returning relative efficiency at a given energy
 @param nuclides_and_activities Vector of (nuclide*, activity) pairs
 @param foreground The spectrum to use for data area calculations
 @param fwhm_form FWHM functional form
 @param fwhm_coefficients FWHM function coefficients
 @param fwhm_lower_energy Lower energy bound for FWHM validity (use -1 if no limit)
 @param fwhm_upper_energy Upper energy bound for FWHM validity (use -1 if no limit)
 @param lowest_energy Lower energy bound for gamma lines
 @param highest_energy Upper energy bound for gamma lines
 @param settings Clustering and ROI parameters

 @return Vector of RoiRange objects
 */
vector<RelActCalcAuto::RoiRange> cluster_gammas_to_rois(
    const function<double(double)> &rel_eff_fcn,
    const vector<pair<RelActCalcAuto::SrcVariant, double>> &sources_and_activities,
    const shared_ptr<const SpecUtils::Measurement> &foreground,
    const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
    const vector<float> &fwhm_coefficients,
    const double fwhm_lower_energy,
    const double fwhm_upper_energy,
    const double lowest_energy,
    const double highest_energy,
    const GammaClusteringSettings &settings )
{
  vector<RelActCalcAuto::RoiRange> result_rois;

  // Collect all gamma lines with their expected counts
  vector<pair<double,double>> gammas_by_counts;  // (energy, expected_counts)

  for( const auto &src_act : sources_and_activities )
  {
    const RelActCalcAuto::SrcVariant src = src_act.first;
    const double activity = src_act.second;

    if( RelActCalcAuto::is_null( src ) || (activity <= 0.0) )
      continue;

    const double age = get_source_age( src, -1.0 );

    if( PeakFitImprove::debug_printout )
    {
      cerr << "cluster_gammas_to_rois: Source " << RelActCalcAuto::to_name( src ) << ", activity=" << activity
           << ", age=" << (age / PhysicalUnits::second) << " seconds ("
           << (age / PhysicalUnits::year) << " years)" << endl;
    }

    const vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src, activity, age );

    if( PeakFitImprove::debug_printout )
    {
      cerr << "  " << photons.size() << " photons from " << RelActCalcAuto::to_name( src ) << ", energy range ["
           << lowest_energy << ", " << highest_energy << "] keV" << endl;
      // Show photons near 807 keV
      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        //if( photon.energy >= 800.0 && photon.energy <= 820.0 )
        //{
        //  const double rel_eff = rel_eff_fcn( photon.energy );
        //  cerr << "    *** Photon near 807 keV: " << photon.energy << " keV, rate=" << photon.numPerSecond
        //       << " /s, rel_eff=" << rel_eff << ", est_counts=" << (photon.numPerSecond * rel_eff) << endl;
        //}
      }
    }

    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      if( (photon.energy < lowest_energy)
         || (photon.energy > highest_energy)
         || (photon.numPerSecond <= std::numeric_limits<double>::epsilon()) )
      {
        continue;
      }

      const double rel_eff = rel_eff_fcn( photon.energy );
      if( rel_eff <= 0.0 )
        continue;

      gammas_by_counts.emplace_back( photon.energy, photon.numPerSecond * rel_eff );
    }//for( const SandiaDecay::EnergyRatePair &photon : photons )
  }//for( const auto &src_act : sources_and_activities )

  if( gammas_by_counts.empty() )
    return result_rois;

  // Create a copy sorted by energy for efficient lookup
  vector<pair<double,double>> gammas_by_energy = gammas_by_counts;
  const auto lessThanByEnergy = []( const pair<double,double> &lhs, const pair<double,double> &rhs ) {
    return lhs.first < rhs.first;
  };
  std::sort( begin(gammas_by_energy), end(gammas_by_energy), lessThanByEnergy );

  // Sort gammas by expected counts (highest first)
  std::sort( begin(gammas_by_counts), end(gammas_by_counts),
    []( const pair<double,double> &lhs, const pair<double,double> &rhs ) {
      return lhs.second > rhs.second;
  } );

  if( PeakFitImprove::debug_printout )
  {
    cerr << "cluster_gammas_to_rois: Input gammas (" << gammas_by_counts.size() << " total):" << endl;
    //for( const auto &gc : gammas_by_energy )
    //{
      // Show all gammas, but highlight those in 800-820 keV range
    //  if( gc.first >= 800.0 && gc.first <= 820.0 )
    //    cerr << "  *** " << gc.first << " keV, est_counts=" << gc.second << " ***" << endl;
    //}
    // Also show top 10 by counts
    cerr << "  Top 20 by expected counts:" << endl;
    for( size_t i = 0; i < std::min( gammas_by_counts.size(), size_t(20) ); ++i )
      cerr << "    " << gammas_by_counts[i].first << " keV, est_counts=" << gammas_by_counts[i].second << endl;
  }

  // Cluster gamma lines
  vector<ClusteredGammaInfo> clustered_gammas;

  // Check if we have a valid FWHM energy range for clamping
  const bool have_fwhm_range = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0);

  for( const pair<double,double> &energy_counts : gammas_by_counts )
  {
    auto ene_pos = std::lower_bound( begin(gammas_by_energy), end(gammas_by_energy),
                                    energy_counts, lessThanByEnergy );
    if( ene_pos == end(gammas_by_energy) )
      continue;

    if( ene_pos->first != energy_counts.first )
    {
      // Already removed from gammas_by_energy (absorbed into another cluster)
      if( PeakFitImprove::debug_printout && (energy_counts.first >= 800.0) && (energy_counts.first <= 820.0) )
      {
        cerr << "cluster_gammas_to_rois: Gamma at " << energy_counts.first
             << " keV already absorbed into another cluster" << endl;
      }
      continue;
    }

    const double energy = energy_counts.first;
    const double counts = energy_counts.second;

    // Clamp energy to valid FWHM range to avoid extrapolation issues
    const double fwhm_eval_energy = have_fwhm_range
        ? std::clamp( energy, fwhm_lower_energy, fwhm_upper_energy )
        : energy;

    const float fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(fwhm_eval_energy), fwhm_form, fwhm_coefficients );

    // Check for invalid FWHM (NaN, zero, or negative)
    if( !std::isfinite(fwhm) || (fwhm <= 0.0f) )
    {
      if( PeakFitImprove::debug_printout )
        cerr << "Warning: Invalid FWHM=" << fwhm << " at energy=" << energy << " keV, skipping gamma line" << endl;
      continue;
    }

    const double sigma = fwhm / PhysicalUnits::fwhm_nsigma;

    const double lower = std::max( lowest_energy, energy - settings.cluster_num_sigma * sigma );
    const double upper = std::min( highest_energy, energy + settings.cluster_num_sigma * sigma );

    // Find all gammas in this range
    const auto start_remove = std::lower_bound( begin(gammas_by_energy), end(gammas_by_energy),
                                               make_pair(lower, 0.0), lessThanByEnergy );
    const auto end_remove = std::upper_bound( begin(gammas_by_energy), end(gammas_by_energy),
                                             make_pair(upper, 0.0), lessThanByEnergy );

    const double counts_in_region = std::accumulate( start_remove, end_remove, 0.0,
        []( const double &sum, const pair<double,double> &el ) {
          return sum + el.second;
    } );

    // Capture the gamma lines before erasing them - as separate energy and amplitude arrays
    vector<double> gamma_energies_in_cluster;
    vector<double> gamma_amplitudes_in_cluster;
    for( auto it = start_remove; it != end_remove; ++it )
    {
      gamma_energies_in_cluster.push_back( it->first );
      gamma_amplitudes_in_cluster.push_back( it->second );
    }

    const double data_area = foreground->gamma_integral( static_cast<float>(lower), static_cast<float>(upper) );

    gammas_by_energy.erase( start_remove, end_remove );

    const double signif = (data_area > 0.0) ? (counts_in_region / sqrt(data_area)) : (1.0 / sqrt(counts_in_region));

    // Additional safety check - ensure lower and upper are finite and valid
    if( !std::isfinite(lower) || !std::isfinite(upper) || (lower >= upper) )
    {
      if( PeakFitImprove::debug_printout )
        cerr << "Warning: Invalid cluster bounds [" << lower << ", " << upper << "], skipping" << endl;
      continue;
    }

    const bool passes_data_area = (data_area > settings.min_data_area_keep);
    const bool passes_counts = (counts_in_region > settings.min_est_peak_area_keep);
    const bool passes_signif = (signif > settings.min_est_significance_keep);

    if( passes_data_area && passes_counts && passes_signif )
    {
      ClusteredGammaInfo cluster_info;
      cluster_info.lower = lower;
      cluster_info.upper = upper;
      cluster_info.gamma_energies = std::move( gamma_energies_in_cluster );
      cluster_info.gamma_amplitudes = std::move( gamma_amplitudes_in_cluster );
      clustered_gammas.push_back( std::move( cluster_info ) );
    }
    
    if( PeakFitImprove::debug_printout )
    {
      const string status_str = (passes_data_area && passes_counts && passes_signif) ? "Accepted" : "Rejected";
      // Print why this cluster was rejected
      cerr << "cluster_gammas_to_rois: " << status_str << " [" << std::fixed << std::setprecision(1) << lower << ", " << upper << "] keV (e="
           << energy << " keV): "
           << "data=" << data_area << (passes_data_area ? " > " : " < ") << settings.min_data_area_keep << "; "
           << "est_counts=" << counts_in_region << (passes_counts ? " > " : " < ") << settings.min_est_peak_area_keep << "; "
         << "sig=" << signif << (passes_signif ? " > " : " < ") << settings.min_est_significance_keep << "; "
         << endl;
    }
  }//for( const pair<double,double> &energy_counts : gammas_by_counts )

  // Sort by lower energy
  std::sort( begin(clustered_gammas), end(clustered_gammas),
    []( const ClusteredGammaInfo &lhs, const ClusteredGammaInfo &rhs ) {
      return lhs.lower < rhs.lower;
  } );

  if( PeakFitImprove::debug_printout )
  {
    cerr << "cluster_gammas_to_rois: Initial " << clustered_gammas.size() << " clustered ROIs:" << endl;
    for( size_t i = 0; i < clustered_gammas.size(); ++i )
    {
      const ClusteredGammaInfo &c = clustered_gammas[i];
      cerr << "  [" << i << "] range=[" << c.lower << ", " << c.upper << "] keV, "
           << c.gamma_energies.size() << " gammas: ";
      for( size_t j = 0; j < std::min( c.gamma_energies.size(), size_t(5) ); ++j )
        cerr << c.gamma_energies[j] << (j + 1 < c.gamma_energies.size() ? ", " : "");
      if( c.gamma_energies.size() > 5 )
        cerr << "... (" << c.gamma_energies.size() - 5 << " more)";
      cerr << endl;
    }
  }

  // At this point we have significant peaks, roughly as would be observed (the dominant peak, and then any smaller
  // peaks within cluster_num_sigma of the dominant peak).
  // After merging overlapping clusters and breaking up overly-wide clusters, we calculate weighted mean and effective
  // FWHM (accounting for both individual peak widths and the spread of gamma energies), then set the final ROI bounds.

  
  // Update cluster bounds based on weighted mean and effective FWHM
  // The effective FWHM accounts for both individual peak widths and the spread of gamma line energies
  // Uses the law of total variance for a mixture of Gaussians:
  //   σ²_total = (weighted avg of individual variances) + (weighted variance of the means)
  for( ClusteredGammaInfo &cluster : clustered_gammas )
  {
    if( cluster.gamma_energies.empty() )
      continue;
    
    // Calculate weighted mean energy (weighted by expected amplitude)
    double sum_weighted_energy = 0.0;
    double sum_weights = 0.0;
    for( size_t i = 0; i < cluster.gamma_energies.size(); ++i )
    {
      sum_weighted_energy += cluster.gamma_energies[i] * cluster.gamma_amplitudes[i];
      sum_weights += cluster.gamma_amplitudes[i];
    }
    
    if( sum_weights <= 0.0 )
      continue;
    
    const double weighted_mean = sum_weighted_energy / sum_weights;
    
    // Calculate effective variance using law of total variance for mixture of Gaussians
    // σ²_total = E[Var(X|I)] + Var(E[X|I])
    //          = (weighted avg of individual σ²) + (weighted variance of means)
    double sum_weighted_var = 0.0;        // For E[Var(X|I)] - weighted average of individual variances
    double sum_weighted_sq_dev = 0.0;     // For Var(E[X|I]) - weighted variance of the means
    
    for( size_t i = 0; i < cluster.gamma_energies.size(); ++i )
    {
      const double energy = cluster.gamma_energies[i];
      const double amplitude = cluster.gamma_amplitudes[i];
      
      // Calculate FWHM at this gamma energy (clamped to valid range)
      const double energy_clamped = have_fwhm_range
      ? std::clamp( energy, fwhm_lower_energy, fwhm_upper_energy )
      : energy;
      const double fwhm_i = DetectorPeakResponse::peakResolutionFWHM(
                                                                     static_cast<float>(energy_clamped), fwhm_form, fwhm_coefficients );
      
      if( !std::isfinite(fwhm_i) || (fwhm_i <= 0.0) )
        continue;
      
      const double sigma_i = fwhm_i / PhysicalUnits::fwhm_nsigma;
      const double var_i = sigma_i * sigma_i;
      
      // Weighted average of individual variances
      sum_weighted_var += amplitude * var_i;
      
      // Weighted squared deviation from weighted mean
      const double dev = energy - weighted_mean;
      sum_weighted_sq_dev += amplitude * dev * dev;
    }
    
    // Total variance = weighted avg variance + weighted variance of means
    const double total_var = (sum_weighted_var + sum_weighted_sq_dev) / sum_weights;
    const double effective_sigma = std::sqrt( total_var );
    const double effective_fwhm = effective_sigma * PhysicalUnits::fwhm_nsigma;
    
    if( !std::isfinite(effective_fwhm) || (effective_fwhm <= 0.0) )
      continue;
    
    // Set new bounds based on weighted mean and effective FWHM
    const double old_lower = cluster.lower;
    const double old_upper = cluster.upper;
    cluster.lower = std::max( lowest_energy, weighted_mean - settings.roi_width_num_fwhm_lower * effective_fwhm );
    cluster.upper = std::min( highest_energy, weighted_mean + settings.roi_width_num_fwhm_upper * effective_fwhm );
    
    if( PeakFitImprove::debug_printout )
    {
      cerr << "cluster_gammas_to_rois: Setting effective FWHM bounds for cluster with "
      << cluster.gamma_energies.size() << " gammas:" << endl
      << "  weighted_mean=" << weighted_mean << " keV, effective_fwhm=" << effective_fwhm << " keV" << endl
      << "  old range=[" << old_lower << ", " << old_upper << "] keV"
      << " -> new range=[" << cluster.lower << ", " << cluster.upper << "] keV" << endl;
    }
  }//for( ClusteredGammaInfo &cluster : final_clusters )
  
  
  // Merge overlapping clusters
  vector<ClusteredGammaInfo> merged_clusters;
  for( const ClusteredGammaInfo &cluster : clustered_gammas )
  {
    if( merged_clusters.empty() || (cluster.lower > merged_clusters.back().upper) )
    {
      merged_clusters.push_back( cluster );
    }
    else
    {
      if( PeakFitImprove::debug_printout )
      {
        cerr << "cluster_gammas_to_rois: Merging cluster [" << cluster.lower << ", " << cluster.upper
             << "] into [" << merged_clusters.back().lower << ", " << merged_clusters.back().upper << "]"
             << " -> new upper=" << std::max( merged_clusters.back().upper, cluster.upper ) << endl;
      }

      merged_clusters.back().upper = std::min( highest_energy, std::max( merged_clusters.back().upper, cluster.upper ) );
      merged_clusters.back().gamma_energies.insert(
        end(merged_clusters.back().gamma_energies),
        begin(cluster.gamma_energies),
        end(cluster.gamma_energies)
      );
      merged_clusters.back().gamma_amplitudes.insert(
        end(merged_clusters.back().gamma_amplitudes),
        begin(cluster.gamma_amplitudes),
        end(cluster.gamma_amplitudes)
      );
    }
  }//for( const ClusteredGammaInfo &cluster : clustered_gammas )

  // Validate merged clusters - ensure merge didn't introduce NaN values
  merged_clusters.erase(
    std::remove_if( begin(merged_clusters), end(merged_clusters),
      []( const ClusteredGammaInfo &cluster ) {
        const bool invalid = !std::isfinite(cluster.lower) || !std::isfinite(cluster.upper) || (cluster.lower >= cluster.upper);
        if( invalid && PeakFitImprove::debug_printout )
          cerr << "Warning: Removing invalid merged cluster with bounds [" << cluster.lower << ", " << cluster.upper << "]" << endl;
        return invalid;
      } ),
    end(merged_clusters)
  );

  if( PeakFitImprove::debug_printout )
  {
    cerr << "cluster_gammas_to_rois: After merging, " << merged_clusters.size() << " clusters:" << endl;
    for( size_t i = 0; i < merged_clusters.size(); ++i )
    {
      const ClusteredGammaInfo &c = merged_clusters[i];
      cerr << "  [" << i << "] range=[" << c.lower << ", " << c.upper << "] keV ("
           << (c.upper - c.lower) << " keV wide), " << c.gamma_energies.size() << " gammas" << endl;
    }
  }

  // Break up ROIs that are too wide
  vector<ClusteredGammaInfo> final_clusters;
  for( ClusteredGammaInfo &cluster : merged_clusters )
  {
    const double mid_energy = 0.5 * (cluster.lower + cluster.upper);

    // Clamp mid_energy to valid FWHM range to avoid extrapolation
    const double mid_fwhm_eval_energy = have_fwhm_range
        ? std::clamp( mid_energy, fwhm_lower_energy, fwhm_upper_energy )
        : mid_energy;

    const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(mid_fwhm_eval_energy), fwhm_form, fwhm_coefficients );
    const double max_width = settings.max_fwhm_width * mid_fwhm;
    const double current_width = cluster.upper - cluster.lower;

    if( PeakFitImprove::debug_printout )
    {
      cerr << "cluster_gammas_to_rois: Cluster [" << cluster.lower << ", " << cluster.upper
           << "] keV: width=" << current_width << " keV, mid_fwhm=" << mid_fwhm
           << " keV, max_fwhm_width=" << settings.max_fwhm_width
           << ", max_width=" << max_width << " keV, "
           << (current_width <= max_width ? "NOT breaking" : "NEEDS breaking") << endl;
    }

    if( current_width <= max_width )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    // Need to break up this cluster using synthetic spectrum-based breakpoint selection
    // Sort gamma energies and amplitudes together by energy
    {
      // Create index array for sorting
      vector<size_t> indices( cluster.gamma_energies.size() );
      std::iota( begin(indices), end(indices), 0 );
      std::sort( begin(indices), end(indices),
        [&cluster]( size_t a, size_t b ) {
          return cluster.gamma_energies[a] < cluster.gamma_energies[b];
        } );

      // Reorder both arrays based on sorted indices
      vector<double> sorted_energies( cluster.gamma_energies.size() );
      vector<double> sorted_amplitudes( cluster.gamma_amplitudes.size() );
      for( size_t i = 0; i < indices.size(); ++i )
      {
        sorted_energies[i] = cluster.gamma_energies[indices[i]];
        sorted_amplitudes[i] = cluster.gamma_amplitudes[indices[i]];
      }
      cluster.gamma_energies = std::move( sorted_energies );
      cluster.gamma_amplitudes = std::move( sorted_amplitudes );
    }

    if( cluster.gamma_energies.size() <= 1 )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    // Get channel range for this cluster
    const size_t start_channel = foreground->find_gamma_channel( static_cast<float>(cluster.lower) );
    const size_t end_channel = foreground->find_gamma_channel( static_cast<float>(cluster.upper) );
    const size_t num_channels = end_channel - start_channel;

    if( num_channels < 3 )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    const float *channel_energies = foreground->channel_energies()->data();

    // Create fwhm_at_energy lambda for the helper functions
    const auto fwhm_at_energy = [&]( double energy ) -> double {
      const double clamped_energy = have_fwhm_range
        ? std::clamp( energy, fwhm_lower_energy, fwhm_upper_energy )
        : energy;
      return DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(clamped_energy), fwhm_form, fwhm_coefficients );
    };

    // Build synthetic spectrum from expected Gaussians (same binning as data)
    vector<double> synthetic = build_synthetic_spectrum(
      cluster.gamma_energies,
      cluster.gamma_amplitudes,
      fwhm_at_energy,
      channel_energies,
      start_channel,
      num_channels );

    // Find local minima in synthetic spectrum (with significance computed)
    const double start_energy = channel_energies[start_channel];
    const double end_energy = channel_energies[end_channel];
    const bool debug_this_cluster = PeakFitImprove::debug_printout &&
      ((start_energy >= 60.0 && end_energy <= 140.0) ||
       (start_energy >= 140.0 && start_energy <= 160.0));

    if( debug_this_cluster )
    {
      cerr << "DEBUG: Breaking cluster [" << cluster.lower << ", " << cluster.upper
           << "] keV: start_channel=" << start_channel
           << ", end_channel=" << end_channel
           << ", num_channels=" << num_channels << endl;
    }

    vector<LocalMinimum> minima = find_synthetic_minima(
      synthetic,
      start_channel,
      foreground,
      channel_energies,
      fwhm_at_energy,
      settings.break_check_fwhm_fraction );

    if( debug_this_cluster )
    {
      cerr << "DEBUG: Found " << minima.size() << " minima in cluster ["
           << cluster.lower << ", " << cluster.upper << "] keV:" << endl;
      for( size_t i = 0; i < minima.size(); ++i )
      {
        const LocalMinimum &m = minima[i];
        cerr << "  [" << i << "] channel=" << m.channel
             << ", energy=" << channel_energies[m.channel] << " keV"
             << ", synthetic_value=" << m.synthetic_value
             << ", depth_score=" << m.depth_score
             << ", stat_sig=" << m.statistical_significance << endl;
      }
    }

    // If no minima found, don't break the ROI
    if( minima.empty() )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    //Currently debugging the 150.1 keV ROI - maybe things are fixed so that it wont be split badly, but maybe not - extending the clusters lower/upper according to effective mean and FWHM not tested, and breaking up large ROIs also totally untested.
    
    // Sort by statistical_significance (ascending - least significant first)
    // Use depth_score as tiebreaker (descending - deepest first)
    // Note: We cannot use a threshold-based tie because it violates strict weak ordering
    //       (transitivity of equivalence fails when A~B and B~C but A!~C)
    std::sort( begin(minima), end(minima),
      []( const LocalMinimum &a, const LocalMinimum &b ) {
        if( a.statistical_significance != b.statistical_significance )
          return a.statistical_significance < b.statistical_significance;
        return a.depth_score > b.depth_score;
      });

    // Select breakpoints, validating significant peaks exist between them
    vector<size_t> selected_breakpoint_channels;
    size_t current_lower_ch = start_channel;

    for( const LocalMinimum &candidate : minima )
    {
      // Check if significant peak exists between current_lower_ch and this candidate
      const bool has_left_peak = has_significant_peak_between(
        current_lower_ch,
        candidate.channel,
        synthetic,
        start_channel,
        foreground,
        channel_energies,
        fwhm_at_energy,
        settings.break_check_fwhm_fraction,
        settings.break_peak_significance_threshold );

      // Check if significant peak exists between this candidate and ROI upper
      const bool has_right_peak = has_significant_peak_between(
        candidate.channel,
        end_channel,
        synthetic,
        start_channel,
        foreground,
        channel_energies,
        fwhm_at_energy,
        settings.break_check_fwhm_fraction,
        settings.break_peak_significance_threshold );

      if( debug_this_cluster )
      {
        cerr << "DEBUG: Evaluating breakpoint candidate at channel " << candidate.channel
             << " (" << channel_energies[candidate.channel] << " keV): "
             << "has_left_peak=" << (has_left_peak ? "YES" : "NO")
             << ", has_right_peak=" << (has_right_peak ? "YES" : "NO")
             << ", accepted=" << ((has_left_peak && has_right_peak) ? "YES" : "NO") << endl;
      }

      if( has_left_peak && has_right_peak )
      {
        selected_breakpoint_channels.push_back( candidate.channel );
        current_lower_ch = candidate.channel;

        // Check if we've satisfied the max_fwhm_width constraint for all sub-ROIs
        bool all_rois_ok = true;
        size_t prev_ch = start_channel;
        for( size_t bp_ch : selected_breakpoint_channels )
        {
          const double sub_lower = channel_energies[prev_ch];
          const double sub_upper = channel_energies[bp_ch];
          const double sub_width = sub_upper - sub_lower;
          const double sub_center = 0.5 * (sub_lower + sub_upper);
          const double sub_fwhm = fwhm_at_energy( sub_center );
          if( sub_width > settings.max_fwhm_width * sub_fwhm )
            all_rois_ok = false;
          prev_ch = bp_ch;
        }
        // Check last sub-ROI
        {
          const double sub_lower = channel_energies[prev_ch];
          const double sub_upper = channel_energies[end_channel];
          const double sub_width = sub_upper - sub_lower;
          const double sub_center = 0.5 * (sub_lower + sub_upper);
          const double sub_fwhm = fwhm_at_energy( sub_center );
          if( sub_width > settings.max_fwhm_width * sub_fwhm )
            all_rois_ok = false;
        }

        if( all_rois_ok )
          break;
      }
    }//for( const LocalMinimum &candidate : minima )

    // Create sub-clusters based on selected breakpoint channels
    if( selected_breakpoint_channels.empty() )
    {
      final_clusters.push_back( std::move( cluster ) );
    }
    else
    {
      std::sort( begin(selected_breakpoint_channels), end(selected_breakpoint_channels) );

      // Convert breakpoint channels to energies
      vector<double> breakpoint_energies;
      for( size_t bp_ch : selected_breakpoint_channels )
        breakpoint_energies.push_back( channel_energies[bp_ch] );

      if( PeakFitImprove::debug_printout )
      {
        cerr << "cluster_gammas_to_rois: Breaking up cluster [" << cluster.lower << ", " << cluster.upper
             << "] at " << breakpoint_energies.size() << " breakpoints: ";
        for( size_t i = 0; i < breakpoint_energies.size(); ++i )
          cerr << breakpoint_energies[i] << (i + 1 < breakpoint_energies.size() ? ", " : "");
        cerr << " keV" << endl;
      }

      double seg_start = cluster.lower;
      size_t gamma_start_idx = 0;

      for( const double bp_energy : breakpoint_energies )
      {
        ClusteredGammaInfo sub_cluster;
        sub_cluster.lower = std::max( lowest_energy, seg_start );
        sub_cluster.upper = std::min( highest_energy, bp_energy );

        // Add gamma lines that fall within this segment
        for( size_t j = gamma_start_idx; j < cluster.gamma_energies.size(); ++j )
        {
          if( cluster.gamma_energies[j] >= bp_energy )
          {
            gamma_start_idx = j;
            break;
          }
          sub_cluster.gamma_energies.push_back( cluster.gamma_energies[j] );
          sub_cluster.gamma_amplitudes.push_back( cluster.gamma_amplitudes[j] );
          if( j == cluster.gamma_energies.size() - 1 )
            gamma_start_idx = j + 1;
        }

        if( !sub_cluster.gamma_energies.empty() )
          final_clusters.push_back( std::move( sub_cluster ) );

        seg_start = bp_energy;
      }

      // Add the last segment
      ClusteredGammaInfo sub_cluster;
      sub_cluster.lower = std::max( lowest_energy, seg_start );
      sub_cluster.upper = std::min( highest_energy, cluster.upper );
      for( size_t j = gamma_start_idx; j < cluster.gamma_energies.size(); ++j )
      {
        sub_cluster.gamma_energies.push_back( cluster.gamma_energies[j] );
        sub_cluster.gamma_amplitudes.push_back( cluster.gamma_amplitudes[j] );
      }

      if( !sub_cluster.gamma_energies.empty() )
        final_clusters.push_back( std::move( sub_cluster ) );
    }
  }//for( ClusteredGammaInfo &cluster : merged_clusters )


  // Developer check: Validate final_clusters don't overlap
#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < final_clusters.size(); ++i )
  {
    const ClusteredGammaInfo &prev = final_clusters[i - 1];
    const ClusteredGammaInfo &curr = final_clusters[i];
    if( curr.lower < prev.upper )
    {
      cerr << "ERROR: final_clusters[" << (i-1) << "] and [" << i << "] overlap: "
           << "[" << prev.lower << ", " << prev.upper << "] vs "
           << "[" << curr.lower << ", " << curr.upper << "]" << endl;
      assert( curr.lower >= prev.upper );
    }
  }
#endif

  // Create ROIs from final clusters
  double previous_roi_upper = 0.0;

  for( const ClusteredGammaInfo &cluster : final_clusters )
  {
    RelActCalcAuto::RoiRange roi;

    // Use the pre-calculated bounds from cluster (based on weighted mean and effective FWHM)
    // Constrain lower bound to not overlap with previous ROI, and clamp to valid energy range
    roi.lower_energy = std::max( lowest_energy, std::max( cluster.lower, previous_roi_upper ) );
    roi.upper_energy = std::min( highest_energy, cluster.upper );

    const double mid_energy = 0.5 * (roi.upper_energy + roi.lower_energy);
    const double mid_energy_clamped = have_fwhm_range
        ? std::clamp( mid_energy, fwhm_lower_energy, fwhm_upper_energy )
        : mid_energy;
    const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(mid_energy_clamped), fwhm_form, fwhm_coefficients );
    const double num_fwhm_wide = (roi.upper_energy - roi.lower_energy) / mid_fwhm;

    if( num_fwhm_wide < settings.min_fwhm_roi )
    {
      if( PeakFitImprove::debug_printout )
      {
        cerr << "cluster_gammas_to_rois: Rejected ROI [" << roi.lower_energy << ", " << roi.upper_energy
             << "] keV: too narrow (" << num_fwhm_wide << " FWHM < " << settings.min_fwhm_roi << " min)" << endl;
      }
      // Don't update previous_roi_upper since we're not adding this ROI
      continue;
    }

    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    if( num_fwhm_wide > settings.min_fwhm_quad_cont )
      roi.continuum_type = PeakContinuum::OffsetType::Quadratic;

    // Check if step continuum should be used based on peak area, significance, and left/right comparison
    if( !cluster.gamma_amplitudes.empty() )
    {
      const double max_amplitude = *max_element( begin(cluster.gamma_amplitudes),
                                                  end(cluster.gamma_amplitudes) );
      const double data_area = foreground->gamma_integral( static_cast<float>(roi.lower_energy),
                                                            static_cast<float>(roi.upper_energy) );
      const double est_significance = (data_area > 0.0)
          ? (max_amplitude / sqrt(data_area))
          : 0.0;

      // Must pass both area and significance thresholds before checking left/right comparison
      if( (max_amplitude >= settings.step_cont_min_peak_area)
          && (est_significance >= settings.step_cont_min_peak_significance) )
      {
        if( should_use_step_continuum( cluster, foreground, fwhm_form, fwhm_coefficients,
                                        fwhm_lower_energy, fwhm_upper_energy,
                                        roi.lower_energy, roi.upper_energy,
                                        settings.step_cont_left_right_nsigma ) )
        {
          // Use FlatStep for narrower ROIs, LinearStep for wider ones
          if( num_fwhm_wide > settings.min_fwhm_quad_cont )
            roi.continuum_type = PeakContinuum::OffsetType::LinearStep;
          else
            roi.continuum_type = PeakContinuum::OffsetType::FlatStep;
        }
      }
    }//if( !cluster.gamma_amplitudes.empty() )

    result_rois.push_back( roi );
    previous_roi_upper = roi.upper_energy;
  }//for( const ClusteredGammaInfo &cluster : final_clusters )

  if( PeakFitImprove::debug_printout )
  {
    cerr << "cluster_gammas_to_rois: Final " << result_rois.size() << " ROIs:" << endl;
    for( size_t i = 0; i < result_rois.size(); ++i )
    {
      const RelActCalcAuto::RoiRange &roi = result_rois[i];
      const double width = roi.upper_energy - roi.lower_energy;
      const double mid = 0.5 * (roi.lower_energy + roi.upper_energy);
      const double mid_clamped = have_fwhm_range
          ? std::clamp( mid, fwhm_lower_energy, fwhm_upper_energy ) : mid;
      const double fwhm = DetectorPeakResponse::peakResolutionFWHM(
          static_cast<float>(mid_clamped), fwhm_form, fwhm_coefficients );
      const char *cont_str = "Unknown";
      switch( roi.continuum_type )
      {
        case PeakContinuum::OffsetType::Linear:     cont_str = "Linear";     break;
        case PeakContinuum::OffsetType::Quadratic:  cont_str = "Quadratic";  break;
        case PeakContinuum::OffsetType::FlatStep:   cont_str = "FlatStep";   break;
        case PeakContinuum::OffsetType::LinearStep: cont_str = "LinearStep"; break;
        default: break;
      }
      cerr << "  [" << i << "] range=[" << roi.lower_energy << ", " << roi.upper_energy << "] keV ("
           << width << " keV, " << (width / fwhm) << " FWHM), cont=" << cont_str << endl;
    }
  }

  // Developer check: Validate result ROIs don't overlap
#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < result_rois.size(); ++i )
  {
    const RelActCalcAuto::RoiRange &prev_roi = result_rois[i - 1];
    const RelActCalcAuto::RoiRange &curr_roi = result_rois[i];
    if( curr_roi.lower_energy < prev_roi.upper_energy )
    {
      cerr << "ERROR: result_rois[" << (i-1) << "] and [" << i << "] overlap: "
           << "[" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] vs "
           << "[" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "]" << endl;
      assert( curr_roi.lower_energy >= prev_roi.upper_energy );
    }
  }
#endif

  return result_rois;
}//cluster_gammas_to_rois


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
 @param source_nuclides Vector of source nuclides to fit
 @param long_background Long background spectrum (can be nullptr)
 @param drf Detector response function (can be nullptr, will use generic if needed)
 @param config Configuration for peak fitting parameters
 @param isHPGe Whether this is an HPGe detector
 @return PeakFitResult with status, error message, fit peaks, and solution
 */
PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const SpecUtils::Measurement> &long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitForNuclideConfig &config,
  const bool isHPGe );


// This is development playground for `AnalystChecks::fit_peaks_for_nuclides(...)`
void eval_peaks_for_nuclide( const std::vector<DataSrcInfo> &srcs_info )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    throw runtime_error( "Failed to open SandiaDecayDataBase" );

  // Configuration for peak fitting - these values will eventually be optimized via GA
  PeakFitForNuclideConfig config;
  // TODO: lets better assign FWHM form that should be used based on detector type
  // TODO: add an enum for how to handle existing ROIs here {ignore, replace_source, no_overlap}

  // Get clustering settings from config
  const GammaClusteringSettings manual_settings = config.get_manual_clustering_settings();
  const GammaClusteringSettings auto_settings = config.get_auto_clustering_settings();


  double total_score = 0.0;

  const bool write_n42 = true;

  // Initialize HTML output file (only if we're writing N42 files)
  std::unique_ptr<ofstream> html_output;
  size_t chart_counter = 0;  // For unique div IDs

  if( write_n42 )
  {
    html_output = std::make_unique<ofstream>( "relact_auto_peakfit_dev.html" );
    try
    {
      D3SpectrumExport::write_html_page_header( *html_output, "RelActAuto Peak Fits", "InterSpec_resources" );
    }catch( std::exception &e )
    {
      cerr << "\n\nError: " << e.what() << endl;
      cerr << "You probably need to symbolically link InterSpec_resources to your CWD." << endl;
      exit(1);
    }

    *html_output << "<body>" << endl;
    *html_output << "<style>"
      << ".TopLinesTable{ margin-top: 25px; margin-left: auto; margin-right: auto; border-collapse: collapse; border: 1px solid black; }" << endl
      << "table, th, td{ border: 1px solid black; }" << endl
      << "fieldset{width: 90vw; margin-left: auto; margin-right: auto; margin-top: 20px;}" << endl
      << "</style>" << endl;
  }

  for( const DataSrcInfo &src_info : srcs_info )
  {
    const InjectSourceInfo &src = src_info.src_info;
    assert( src.src_spectra.size() >= 12 );
    if( src.src_spectra.size() < 12 )
      throw runtime_error( "Unexpected number of measurements." );

    string src_name = src.src_name;
    const auto underscore_pos = src_name.find( '_' );
    if( underscore_pos != string::npos )
    {
      const string last_part = src_name.substr(underscore_pos);
      assert( (last_part == "_Sh")  || (last_part == "_Sh-Point") || (last_part == "_Sh_BPE") || (last_part == "_Sh_PE")
             || (last_part == "_Sh_100g") || (last_part == "_ShHeavy") || (last_part == "_Unsh")
             || (last_part == "_Unsh-Point") || (last_part == "_Phantom")
             || (last_part == "_Unsh_0200") || (last_part == "_Unsh_0300") || (last_part == "_Unsh_0400")
             || (last_part == "_Unsh_0500") || (last_part == "_Unsh_1000") || (last_part == "_Unsh_2000")
             || (last_part == "_Unsh_4000") || (last_part == "_Unsh_5000") || (last_part == "_Unsh_6000")
             || (last_part == "_Unsh_8000") || (last_part == "_Unsh_9000") || (last_part == "_Unsh_9330")
      );
      src_name = src_name.substr(0,underscore_pos);
    }//if( underscore_pos != string::npos )

    vector<RelActCalcAuto::SrcVariant> sources;
    if( src_name == "Tl201woTl202" )
    {
      sources.push_back( db->nuclide("Tl201"));
      sources.push_back( db->nuclide("Tl202")); //I think it does actually have Tl202 in it
    }else if( src_name == "Tl201wTl202" )
    {
      sources.push_back( db->nuclide("Tl201"));
      sources.push_back( db->nuclide("Tl202"));
    }else if( src_name == "Tl201" )
    {
      sources.push_back( db->nuclide("Tl201"));
      sources.push_back( db->nuclide("Tl202"));
    }else if( src_name == "I125" )
    {
      sources.push_back( db->nuclide("I125"));
      sources.push_back( db->nuclide("I126"));
    }else if( src_name == "U233" )
    {
      sources.push_back( db->element( "U" ) );
      sources.push_back( db->nuclide("U232"));
      sources.push_back( db->nuclide("U233"));
    }else if( src_name == "Pu238" )
    {
      sources.push_back( db->element( "Pu" ) );
      sources.push_back( db->nuclide("Pu238"));
      sources.push_back( db->nuclide("Pu239"));
      sources.push_back( db->nuclide("Pu241"));
    }else if( src_name == "Pu239" )
    {
      sources.push_back( db->element( "Pu" ) );
      sources.push_back( db->nuclide("Pu239"));
      sources.push_back( db->nuclide("Pu241"));
    }else if( src_name == "Uore" )
    {
      sources.push_back( db->element( "U" ) );
      sources.push_back( db->nuclide("U235") );
      sources.push_back( db->nuclide("U238") );
      sources.push_back( db->nuclide("Ra226") );
    }else if( src_name == "U235" )
    {
      sources.push_back( db->element( "U" ) );
      sources.push_back( db->nuclide("U235") );
      sources.push_back( db->nuclide("U238") );
    }else if( src_name == "Np237" )
    {
      sources.push_back( db->element( "Np" ) );
      sources.push_back( db->nuclide("Np237") );
    }else if( src_name == "Xe133" )
    {
      sources.push_back( db->nuclide("Xe133") );
      sources.push_back( db->nuclide("Xe133m") );
    }else if( src_name == "Cf252" )
    {
      sources.push_back( db->nuclide("Cf252") );
      
      const ReactionGamma * const rctn_db = ReactionGammaServer::database();
      vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
      if( rctn_db )
        rctn_db->gammas( "H(n,g)", possible_rctns );  //May throw exception if reaction name is invalid
      assert( !possible_rctns.empty() );
      if( !possible_rctns.empty() )
        sources.push_back( possible_rctns[0].reaction );
      
      continue;
    }else if( src_name == "Am241Li" )
    {
      continue;
    }else
    {
      const SandiaDecay::Nuclide * const nuc = db->nuclide( src_name );
      if( !nuc )
      {
        cout << "Unable to get nuclide for '" << src.src_name << "' - so will skip" << endl;
        continue;
      }
      sources.push_back( nuc );
      assert( nuc->symbol == src_name );
    }

    assert( !sources.empty() && !RelActCalcAuto::is_null( sources.front() ) );
    if( sources.empty() || RelActCalcAuto::is_null( sources.front() ) )
      throw runtime_error( "Failed to get a source" );

    const shared_ptr<SpecUtils::SpecFile> &spec_file = src.spec_file;
    const shared_ptr<const SpecUtils::Measurement> &foreground = src.src_spectra.front(); // TODO: we cold loop over all 16 of these histograms
    const shared_ptr<const SpecUtils::Measurement> &long_background = src.long_background;

    const bool isHPGe = true;
    const bool singleThreaded = false;
    shared_ptr<const DetectorPeakResponse> drf = nullptr; //Probably

    try
    {
      shared_ptr<const deque<shared_ptr<const PeakDef>>> dummy_origpeaks;

      const vector<shared_ptr<const PeakDef> > auto_search_peaks
          = ExperimentalAutomatedPeakSearch::search_for_peaks( foreground, drf, dummy_origpeaks, singleThreaded, isHPGe );

      const PeakFitResult curve_results = fit_peaks_for_nuclides(
        auto_search_peaks, foreground, sources, long_background, drf, config, isHPGe );

      if( curve_results.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        throw std::runtime_error( "fit_peaks_for_nuclides failed for all RelEff curve types: " + curve_results.error_message );

      const RelActCalcAuto::RelActAutoSolution &solution = curve_results.solution;
      const std::vector<PeakDef> &fit_peaks = curve_results.observable_peaks;

#if( PERFORM_DEVELOPER_CHECKS )
      {// Begin check for overlapping ROIs in fit_peaks
        // Collect unique continuums (ROIs) and their bounds
        std::map<std::shared_ptr<const PeakContinuum>, std::vector<double>> continuum_to_peak_means;
        for( const PeakDef &peak : fit_peaks )
          continuum_to_peak_means[peak.continuum()].push_back( peak.mean() );

        struct RoiCheckInfo
        {
          double lower_energy;
          double upper_energy;
          std::vector<double> peak_means;
        };

        std::vector<RoiCheckInfo> rois_to_check;
        for( const auto &entry : continuum_to_peak_means )
        {
          RoiCheckInfo info;
          info.lower_energy = entry.first->lowerEnergy();
          info.upper_energy = entry.first->upperEnergy();
          info.peak_means = entry.second;
          rois_to_check.push_back( info );
        }

        // Sort ROIs by lower energy
        std::sort( begin(rois_to_check), end(rois_to_check),
          []( const RoiCheckInfo &a, const RoiCheckInfo &b ) { return a.lower_energy < b.lower_energy; } );

        // Check for overlaps
        bool found_overlap = false;
        for( size_t i = 1; i < rois_to_check.size(); ++i )
        {
          const RoiCheckInfo &prev_roi = rois_to_check[i - 1];
          const RoiCheckInfo &curr_roi = rois_to_check[i];

          if( curr_roi.lower_energy < prev_roi.upper_energy )
          {
            if( !found_overlap )
              cerr << "\n*** OVERLAPPING ROIs DETECTED in fit_peaks (observable_peaks) ***" << endl;
            found_overlap = true;

            cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
            for( double mean : prev_roi.peak_means )
              cerr << mean << " ";
            cerr << "keV" << endl;

            cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
            for( double mean : curr_roi.peak_means )
              cerr << mean << " ";
            cerr << "keV" << endl;
            cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << endl;
          }
        }

        if( found_overlap )
          cerr << "*** END OVERLAPPING ROIs ***\n" << endl;
      }// End check for overlapping ROIs in fit_peaks
#endif

      cout << "Best RelEff curve type solution (" << solution.m_options.rel_eff_curves.front().name << ") selected with chi2/dof="
           << solution.m_chi2 << "/" << solution.m_dof << endl;

      // Score the fit results using only signal photopeaks
      CombinedPeakFitScore combined_score;

      // Calculate FinalFitScore using refactored helper
      combined_score.final_fit_score = FinalFit_GA::calculate_final_fit_score(
        fit_peaks,
        src_info.expected_signal_photopeaks,
        config.num_sigma_contribution
      );

      // Calculate PeakFindAndFitWeights using refactored helper
      combined_score.initial_fit_weights = InitialFit_GA::calculate_peak_find_weights(
        fit_peaks,
        src_info.expected_signal_photopeaks,
        config.num_sigma_contribution
      );

      // Calculate CandidatePeakScore using refactored helper
      combined_score.candidate_peak_score = CandidatePeak_GA::calculate_candidate_peak_score_for_source(
        fit_peaks,
        src_info.expected_signal_photopeaks
      );

      // Correct score to exclude escape peaks from penalty counts
      CandidatePeak_GA::correct_score_for_escape_peaks(
        combined_score.candidate_peak_score,
        src_info.expected_signal_photopeaks
      );

      // Calculate combined final weight (simple sum)
      combined_score.final_weight = combined_score.initial_fit_weights.find_weight
                                    + combined_score.final_fit_score.total_weight 
                                    + combined_score.candidate_peak_score.score;

      // Add to cumulative total score
      total_score += combined_score.final_weight;

      if( PeakFitImprove::debug_printout )
      {
        cout << "Combined score for " << src.src_name << ":" << endl;
        cout << combined_score.print( "combined_score" ) << endl;
      }

      // Write N42 file with foreground, background, and fit peaks (only if fit succeeded)
      if( curve_results.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
      {
        if( write_n42 )
        {
        string outdir = "output_n42_relactauto";
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        outdir = SpecUtils::append_path( outdir, src_info.detector_name );
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        outdir = SpecUtils::append_path( outdir, src_info.location_name );
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        outdir = SpecUtils::append_path( outdir, src_info.live_time_name );
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        // Use the source name for the output file
        const string out_n42 = SpecUtils::append_path( outdir, src.src_name ) + "_relactauto_fit.n42";

        SpecMeas output;
        output.set_manufacturer( src.spec_file->manufacturer() );
        output.set_detector_type( src.spec_file->detector_type() );
        output.set_instrument_model( src.spec_file->instrument_model() );
        output.set_instrument_type( src.spec_file->instrument_type() );
        output.set_filename( src.spec_file->filename() );

        output.remove_measurements( output.measurements() );

        // Generate a unique UUID based on current time + random component
        {
          const auto now_time = std::chrono::system_clock::now();
          const auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
              now_time.time_since_epoch() ).count();
          std::random_device rd;
          std::mt19937 gen( rd() );
          std::uniform_int_distribution<uint32_t> dis( 0, 0xFFFFFFFF );
          std::ostringstream uuid_oss;
          uuid_oss << std::hex << std::setfill('0')
                   << std::setw(12) << now_ms << "-"
                   << std::setw(8) << dis(gen) << "-"
                   << std::setw(8) << dis(gen);
          output.set_uuid( uuid_oss.str() );
        }

        output.add_remark( combined_score.print( "combined_score" ) );
        output.add_remark( "RelActAuto solution chi2/dof=" + std::to_string( solution.m_chi2 ) + "/" + std::to_string( solution.m_dof ) );

        output.set_instrument_model( "RelActAuto Fit Test" );

        // Add foreground spectrum
        shared_ptr<SpecUtils::Measurement> out_foreground = make_shared<SpecUtils::Measurement>( *foreground );
        out_foreground->set_sample_number( 1 );
        out_foreground->set_title( src_name + " - Foreground" );
        const auto now = std::chrono::time_point_cast<std::chrono::microseconds>( std::chrono::system_clock::now() );
        out_foreground->set_start_time( now );
        output.add_measurement( out_foreground, false );

        // Add background spectrum if present
        if( long_background )
        {
          shared_ptr<SpecUtils::Measurement> out_background = make_shared<SpecUtils::Measurement>( *long_background );
          out_background->set_sample_number( 2 );
          out_background->set_title( src_name + " - Background" );
          out_background->set_start_time( now );
          output.add_measurement( out_background, false );
        }

        // Add fit peaks associated with foreground
        deque<shared_ptr<const PeakDef>> peaks;
        for( const PeakDef &p : fit_peaks )
          peaks.push_back( make_shared<PeakDef>( p ) );

        output.setPeaks( peaks, {1} );

        // Add RelActAuto GUI state so the fit can be loaded in InterSpec
        RelActCalcAuto::RelActAutoGuiState gui_state;
        gui_state.options = solution.m_options;
        gui_state.background_subtract = (long_background != nullptr);
        gui_state.show_ref_lines = true;

        // Set display energy range from ROIs or foreground spectrum
        if( !solution.m_final_roi_ranges.empty() )
        {
          gui_state.lower_display_energy = solution.m_final_roi_ranges.front().lower_energy;
          gui_state.upper_display_energy = solution.m_final_roi_ranges.back().upper_energy;
        }
        else
        {
          gui_state.lower_display_energy = foreground->gamma_energy_min();
          gui_state.upper_display_energy = foreground->gamma_energy_max();
        }

        gui_state.note = "RelActAuto fit from PeakFitImprove";
        gui_state.description = "Source: " + src.src_name;

        output.setRelActAutoGuiState( &gui_state );

        output.save2012N42File( out_n42, [=](){
          cerr << "Failed to write '" << out_n42 << "'" << endl;
        });

        if( PeakFitImprove::debug_printout )
          cout << "Wrote N42 file: " << out_n42 << endl;

        // Generate reference lines for the nuclides that were fit
        std::map<std::string, std::string> reference_lines_json;

        // Define colors to rotate through for different nuclides
        const std::vector<Wt::WColor> colors = {
          Wt::WColor(0, 0, 139),      // darkBlue
          Wt::WColor(0, 139, 139),    // darkCyan
          Wt::WColor(0, 100, 0),      // darkGreen
          Wt::WColor(139, 0, 139),    // darkMagenta
          Wt::WColor(255, 140, 0),    // darkOrange
          Wt::WColor(139, 0, 0),      // darkRed
        };

        size_t color_index = 0;

        for( size_t rel_eff_index = 0; rel_eff_index < solution.m_rel_activities.size(); ++rel_eff_index )
        {
          const std::vector<RelActCalcAuto::NuclideRelAct> &nuclides = solution.m_rel_activities[rel_eff_index];

          for( const RelActCalcAuto::NuclideRelAct &nuc_info : nuclides )
          {
            const std::string nuc_name = nuc_info.name();

            // Create RefLineInput for this nuclide
            RefLineInput input;
            input.m_input_txt = nuc_name;

            // Set age if it was fit
            if( nuc_info.age_was_fit && nuc_info.age > 0 )
            {
              input.m_age = std::to_string(nuc_info.age) + " s";
            }

            // Rotate through colors
            input.m_color = colors[color_index % colors.size()];
            color_index++;

            input.m_showGammas = true;
            input.m_showXrays = true;
            input.m_showAlphas = false;
            input.m_showBetas = false;
            input.m_promptLinesOnly = false;
            input.m_lower_br_cutt_off = 0.0; // Only show lines with BR > 1%

            // Generate reference line info
            std::shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );

            if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string ref_json;
              ref_info->toJson( ref_json );
              reference_lines_json[nuc_name] = ref_json;
            }
          }
        }

        {// Begin add reference lines for valid energy range
          const auto [min_valid_energy, max_valid_energy] = find_valid_energy_range( foreground );

          if( (min_valid_energy > 0.0) && (max_valid_energy > min_valid_energy) )
          {
            // Create RefLineInput for minimum energy
            RefLineInput min_input;
            min_input.m_input_txt = std::to_string(min_valid_energy) + " keV";
            min_input.m_color = Wt::WColor(255, 165, 0); // Orange
            min_input.m_showGammas = true;
            min_input.m_showXrays = false;
            min_input.m_showAlphas = false;
            min_input.m_showBetas = false;
            min_input.m_promptLinesOnly = false;
            min_input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> min_ref_info = ReferenceLineInfo::generateRefLineInfo( min_input );
            if( min_ref_info && min_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string min_ref_json;
              min_ref_info->toJson( min_ref_json );
              reference_lines_json["Valid Range Min"] = min_ref_json;
            }

            // Create RefLineInput for maximum energy
            RefLineInput max_input;
            max_input.m_input_txt = std::to_string(max_valid_energy) + " keV";
            max_input.m_color = Wt::WColor(255, 165, 0); // Orange
            max_input.m_showGammas = true;
            max_input.m_showXrays = false;
            max_input.m_showAlphas = false;
            max_input.m_showBetas = false;
            max_input.m_promptLinesOnly = false;
            max_input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> max_ref_info = ReferenceLineInfo::generateRefLineInfo( max_input );
            if( max_ref_info && max_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string max_ref_json;
              max_ref_info->toJson( max_ref_json );
              reference_lines_json["Valid Range Max"] = max_ref_json;
            }
          }
        }// End add reference lines for valid energy range

        {// Begin add red reference lines for each energy in `combined_score.candidate_peak_score.def_expected_but_not_detected`
          const std::vector<ExpectedPhotopeakInfo> &not_detected = combined_score.candidate_peak_score.def_expected_but_not_detected;
          for( size_t i = 0; i < not_detected.size(); ++i )
          {
            const double energy = not_detected[i].effective_energy;

            RefLineInput red_line_input;
            red_line_input.m_input_txt = std::to_string(energy) + " keV";
            red_line_input.m_color = Wt::WColor(255, 0, 0); // Red
            red_line_input.m_showGammas = true;
            red_line_input.m_showXrays = false;
            red_line_input.m_showAlphas = false;
            red_line_input.m_showBetas = false;
            red_line_input.m_promptLinesOnly = false;
            red_line_input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> red_ref_info = ReferenceLineInfo::generateRefLineInfo( red_line_input );
            if( red_ref_info && red_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string red_ref_json;
              red_ref_info->toJson( red_ref_json );
              reference_lines_json["Missing Peak " + std::to_string(i+1)] = red_ref_json;
            }
          }
        }// End add red reference lines for each energy in `combined_score.candidate_peak_score.def_expected_but_not_detected`

        // Setup chart options
        std::string title = src_name + " - RelActAuto Fit";
        std::string dataTitle = "";
        bool useLogYAxis = true;
        bool showVerticalGridLines = false;
        bool showHorizontalGridLines = false;
        bool legendEnabled = true;
        bool compactXAxis = true;
        bool showPeakUserLabels = false;
        bool showPeakEnergyLabels = false;
        bool showPeakNuclideLabels = false;
        bool showPeakNuclideEnergyLabels = false;
        bool showEscapePeakMarker = false;
        bool showComptonPeakMarker = false;
        bool showComptonEdgeMarker = false;
        bool showSumPeakMarker = false;
        bool backgroundSubtract = false;

        // Set energy range
        float xMin = 0.0f;
        float xMax = 3000.0f;
        if( !solution.m_final_roi_ranges.empty() )
        {
          xMin = static_cast<float>( solution.m_final_roi_ranges.front().lower_energy );
          xMax = static_cast<float>( solution.m_final_roi_ranges.back().upper_energy );
        }
        else if( foreground )
        {
          xMin = static_cast<float>( foreground->gamma_energy_min() );
          xMax = static_cast<float>( foreground->gamma_energy_max() );
        }

        D3SpectrumExport::D3SpectrumChartOptions options(
          title, "Energy (keV)", "Counts/Channel",
          dataTitle, useLogYAxis,
          showVerticalGridLines, showHorizontalGridLines,
          legendEnabled, compactXAxis,
          showPeakUserLabels, showPeakEnergyLabels, showPeakNuclideLabels,
          showPeakNuclideEnergyLabels, showEscapePeakMarker, showComptonPeakMarker,
          showComptonEdgeMarker, showSumPeakMarker, backgroundSubtract,
          xMin, xMax, reference_lines_json
        );

        // Setup spectrum data with peaks
        D3SpectrumExport::D3SpectrumOptions foreground_opts;
        foreground_opts.line_color = "black";
        foreground_opts.title = src_name;
        foreground_opts.display_scale_factor = 1.0;
        foreground_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;

        // Convert fit_peaks to shared_ptr vector for peak_json
        vector<shared_ptr<const PeakDef>> fit_peaks_ptrs;
        for( const PeakDef &p : fit_peaks )
          fit_peaks_ptrs.push_back( make_shared<PeakDef>( p ) );

        foreground_opts.peaks_json = PeakDef::peak_json( fit_peaks_ptrs, foreground, Wt::WColor(), false );

        // Create unique div ID
        const string div_id = "chart_" + std::to_string(chart_counter);
        chart_counter++;

        // Write HTML fieldset and div
        *html_output << "<fieldset style=\"\">" << endl
          << "<legend>" << src.src_name << " (chi2/dof=" << solution.m_chi2 << "/" << solution.m_dof << ")</legend>" << endl;
        *html_output << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>" << endl;
        *html_output << "<script>" << endl;

        // Initialize chart
        D3SpectrumExport::write_js_for_chart( *html_output, div_id, options.m_dataTitle, options.m_xAxisTitle, options.m_yAxisTitle );

        // Set spectrum data
        std::vector< std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
        measurements.emplace_back( foreground.get(), foreground_opts );

        // Add background if present
        if( long_background )
        {
          D3SpectrumExport::D3SpectrumOptions background_opts;
          background_opts.line_color = "steelblue";
          background_opts.title = "Background";
          background_opts.display_scale_factor = foreground->live_time() / long_background->live_time();
          background_opts.spectrum_type = SpecUtils::SpectrumType::Background;
          measurements.emplace_back( long_background.get(), background_opts );
        }

        D3SpectrumExport::write_and_set_data_for_chart( *html_output, div_id, measurements );

        // Write resize handler
        *html_output << R"delim(
const resizeChart)delim" << chart_counter-1 << R"delim( = function(){
  let height = window.innerHeight;
  let width = document.documentElement.clientWidth;
  let el = spec_chart_)delim" << div_id << R"delim(.chart;
  el.style.width = 0.8*width + "px";
  el.style.height = Math.min(500,Math.max(250, Math.min(0.4*width,height-175))) + "px";
  el.style.marginLeft = 0.05*width + "px";
  el.style.marginRight = 0.05*width + "px";
  )delim"
          << "  spec_chart_" << div_id << R"delim(.handleResize();
};

window.addEventListener('resize', resizeChart)delim" << chart_counter-1 << R"delim();
)delim" << endl;

        // Set chart options (this creates reference_lines_<div_id> variable)
        D3SpectrumExport::write_set_options_for_chart( *html_output, div_id, options );

        // Apply reference lines
        *html_output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );" << endl;

        // Trigger initial resize
        *html_output << "resizeChart" << chart_counter-1 << "();" << endl;
        *html_output << "</script>" << endl;

        // Add fit score information with three sections
        *html_output << "<div style=\"margin: 10px auto; max-width: 800px;\">" << endl;
        *html_output << "<h4>Fit Quality Metrics</h4>" << endl;

        // === Combined Score (prominent display) ===
        *html_output << "<div style=\"background: #f0f8ff; padding: 10px; margin-bottom: 15px; border: 2px solid #4682b4; border-radius: 5px;\">" << endl;
        *html_output << "<h5 style=\"margin: 0 0 10px 0;\">Combined Score</h5>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\"><strong>Final Weight (Combined)</strong></td>"
                     << "<td style=\"text-align: right; padding: 5px;\"><strong>"
                     << combined_score.final_weight << "</strong></td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\">Component 0: Peak Find Score</td>"
                     << "<td style=\"text-align: right; padding: 5px;\">"
                     << combined_score.candidate_peak_score.score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\">Component 1: Find Weight</td>"
                     << "<td style=\"text-align: right; padding: 5px;\">"
                     << combined_score.initial_fit_weights.find_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\">Component 2: Total Weight (Avrg nsigma off)</td>"
                     << "<td style=\"text-align: right; padding: 5px;\">"
                     << combined_score.final_fit_score.total_weight << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</div>" << endl;

        // === Final Fit Metrics (collapsible) ===
        *html_output << "<details>" << endl;
        *html_output << "<summary style=\"cursor: pointer; font-weight: bold; padding: 5px; background: #f5f5f5; border: 1px solid #ddd;\">"
                     << "Final Fit Quality Metrics</summary>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em; margin-top: 5px;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Area Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.area_score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Width Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.width_score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Position Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.position_score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Ignored Unexpected Peaks</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.ignored_unexpected_peaks << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Unexpected Peaks Sum Significance</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.unexpected_peaks_sum_significance << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Number of Peaks Used</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.num_peaks_used << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Total Weight (Area)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.total_weight << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</details>" << endl;

        // === Peak Finding Metrics (collapsible) ===
        *html_output << "<details style=\"margin-top: 10px;\">" << endl;
        *html_output << "<summary style=\"cursor: pointer; font-weight: bold; padding: 5px; background: #f5f5f5; border: 1px solid #ddd;\">"
                     << "Peak Finding & Area Accuracy Metrics</summary>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em; margin-top: 5px;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Find Weight (Score)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.find_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Def Wanted Area Weight (Mean)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.def_wanted_area_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Def Wanted Area Weight (Median)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.def_wanted_area_median_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Maybe Wanted Area Weight (Mean)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.maybe_wanted_area_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Maybe Wanted Area Weight (Median)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.maybe_wanted_area_median_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Not Wanted Area Weight (Mean)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.not_wanted_area_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Not Wanted Area Weight (Median)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.not_wanted_area_median_weight << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</details>" << endl;

        // === Candidate Peak Metrics (collapsible) ===
        *html_output << "<details style=\"margin-top: 10px;\">" << endl;
        *html_output << "<summary style=\"cursor: pointer; font-weight: bold; padding: 5px; background: #f5f5f5; border: 1px solid #ddd;\">"
                     << "Candidate Peak Detection Metrics</summary>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em; margin-top: 5px;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Number of Peaks Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_peaks_found << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Definitely Wanted Peaks Not Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_def_wanted_not_found << "</td></tr>" << endl;
        if( !combined_score.candidate_peak_score.def_expected_but_not_detected.empty() )
        {
          *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Definitely Wanted Peaks Not Found (Energies)</td>"
                       << "<td style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">";
          for( size_t i = 0; i < combined_score.candidate_peak_score.def_expected_but_not_detected.size(); ++i )
          {
            if( i > 0 )
              *html_output << ", ";
            *html_output << std::fixed << std::setprecision(2) 
                         << combined_score.candidate_peak_score.def_expected_but_not_detected[i].effective_energy << " keV";
          }
          *html_output << "</td></tr>" << endl;
        }
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Definitely Wanted Peaks Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_def_wanted_peaks_found << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Possibly Accepted Peaks Not Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_possibly_accepted_peaks_not_found << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Extra Peaks (Not Expected)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_extra_peaks << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</details>" << endl;
        
        *html_output << "<p>Definitely Wanted Peaks Not Found: "
                     << combined_score.candidate_peak_score.num_def_wanted_not_found << "</p>" << endl;

        *html_output << "</div>" << endl;

        *html_output << "</fieldset>" << endl;
        }//if( write_n42 )
      }//if( curve_results.status == Success )

    }catch( const std::exception &e )
    {
      //throw std::runtime_error("Error in fit_peaks_for_nuclides: " + string(e.what()));
      cerr << e.what() << endl;
    }
  }//for( const DataSrcInfo &src_info : srcs_info )

  // Print total score summary
  cout << endl << "========================================" << endl;
  cout << "Total score across all " << srcs_info.size() << " sources: " << total_score << endl;
  cout << "========================================" << endl << endl;

  // Close HTML file
  if( html_output )
  {
    *html_output << "</body>" << endl;
    *html_output << "</html>" << endl;
    html_output->close();

    if( PeakFitImprove::debug_printout )
      cout << "Wrote HTML file: relact_auto_peakfit_dev.html with " << chart_counter << " charts" << endl;
  }
}//void eval_peaks_for_nuclide( const DataSrcInfo &src_info )


}//namespace FitPeaksForNuclideDev
