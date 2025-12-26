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
#include <string>
#include <fstream>
#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp> // For boost::asio::post

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
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
  vector<pair<double,double>> gamma_lines; // (energy, amplitude)
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
};


// Forward declaration
vector<RelActCalcAuto::RoiRange> cluster_gammas_to_rois(
  const function<double(double)> &rel_eff_fcn,
  const vector<pair<const SandiaDecay::Nuclide *, double>> &nuclides_and_activities,
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
  std::vector<PeakDef> fit_peaks;  // from m_peaks_without_back_sub
  RelActCalcAuto::RelActAutoSolution solution;  // only valid if status == Success
};//struct PeakFitResult


// Configuration parameters for peak fitting for nuclides
// These values will be optimized via genetic algorithm for different detector types
struct PeakFitForNuclideConfig {
  // FWHM functional form
  DetectorPeakResponse::ResolutionFnctForm fwhm_functional_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;

  // RelActManual parameters for initial relative efficiency estimation
  double rel_eff_manual_base_rel_eff_uncert = 0.0;
  double initial_nuc_match_cluster_num_sigma = 1.5;
  double manual_eff_cluster_num_sigma = 1.5;

  // RelActManual equation form and order based on number of matched peaks
  size_t initial_manual_relEff_1peak_eqn_order = 0;
  RelActCalc::RelEffEqnForm initial_manual_relEff_1peak_form = RelActCalc::RelEffEqnForm::LnX;

  size_t initial_manual_relEff_2peak_eqn_order = 0;
  RelActCalc::RelEffEqnForm initial_manual_relEff_2peak_form = RelActCalc::RelEffEqnForm::LnX;

  size_t initial_manual_relEff_3peak_eqn_order = 1;
  RelActCalc::RelEffEqnForm initial_manual_relEff_3peak_form = RelActCalc::RelEffEqnForm::LnX;

  bool initial_manual_relEff_4peak_physical_use_hoerl = false;
  size_t initial_manual_relEff_4peak_eqn_order = 2;
  RelActCalc::RelEffEqnForm initial_manual_relEff_4peak_form = RelActCalc::RelEffEqnForm::LnX;

  bool initial_manual_relEff_many_peak_physical_use_hoerl = false;
  size_t initial_manual_relEff_many_peak_eqn_order = 3;
  RelActCalc::RelEffEqnForm initial_manual_relEff_manypeak_form = RelActCalc::RelEffEqnForm::LnX;

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
  double rel_eff_auto_base_rel_eff_uncert = 0.01;  // BR uncertainty for RelActAuto

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
    return settings;
  }

  // RelActAuto relative efficiency model parameters
  RelActCalc::RelEffEqnForm rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  size_t rel_eff_eqn_order = 3;

  // Physical model shielding (empty vectors mean no shielding)
  std::vector<std::shared_ptr<RelActCalc::PhysicalModelShieldInput>> phys_model_self_atten;
  std::vector<std::shared_ptr<RelActCalc::PhysicalModelShieldInput>> phys_model_external_atten;

  // Physical model options (only used when rel_eff_eqn_type == FramPhysicalModel)
  bool nucs_of_el_same_age = true;
  bool phys_model_use_hoerl = true;

  // Fields for RelActAuto options configuration
  bool fit_energy_cal = true;
  std::vector<RelActCalcAuto::NucInputInfo> base_nuclides;  // Set from outside
};//struct PeakFitForNuclideConfig


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
  const std::vector<const SandiaDecay::Nuclide *> &sources,
  const RelActCalcAuto::Options &options,
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
    result.error_message = "No source nuclides provided";
    return result;
  }

  for( const SandiaDecay::Nuclide *nuc : sources )
  {
    if( !nuc )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Null nuclide in sources vector";
      return result;
    }
  }

  try
  {
    // Call RelActAuto::solve with provided options
    RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
      options, foreground, long_background, drf, auto_search_peaks
    );

    // Check if initial solve failed
    if( solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      result.status = solution.m_status;
      result.error_message = solution.m_error_message;
      return result;
    }

    cout << "Initial RelActAuto solution (" << options.rel_eff_curves.front().name << "):" << endl;
    solution.print_summary( std::cout );
    cout << endl;

    // Iteratively refine ROIs using RelActAuto solutions
    // The idea is that each iteration provides a better relative efficiency estimate,
    // which allows us to better identify significant gamma lines and create better ROIs.
    {
      const size_t max_iterations = 3;
      for( size_t iter = 0; iter < max_iterations; ++iter )
      {
        // Use the first rel-eff curve (index 0)
        const size_t rel_eff_index = 0;

        // Create rel_eff lambda from current RelActAuto solution
        const auto auto_rel_eff = [&solution, rel_eff_index]( double energy ) -> double {
          return solution.relative_efficiency( energy, rel_eff_index );
        };

        // Collect nuclides and activities from the current solution
        // m_rel_activities is a 2D vector: [rel_eff_curve_index][nuclide_index]
        vector<pair<const SandiaDecay::Nuclide *, double>> nucs_and_acts;
        if( rel_eff_index < solution.m_rel_activities.size() )
        {
          if( PeakFitImprove::debug_printout )
            cout << "Collecting " << solution.m_rel_activities[rel_eff_index].size() << " nuclides from RelActAuto solution:" << endl;

          for( const RelActCalcAuto::NuclideRelAct &nuc_act : solution.m_rel_activities[rel_eff_index] )
          {
            // source is a SrcVariant which can be a nuclide or element
            const SandiaDecay::Nuclide * const * nuc_ptr_local = std::get_if<const SandiaDecay::Nuclide *>( &nuc_act.source );
            if( nuc_ptr_local && *nuc_ptr_local )
            {
              const double live_time_seconds = foreground->live_time();
              // RelActAuto's rel_activity is per second, need to multiply by live time for clustering
              const double activity_for_clustering = nuc_act.rel_activity * live_time_seconds;

              if( PeakFitImprove::debug_printout )
                cout << "  " << nuc_act.name() << ": rel_activity=" << nuc_act.rel_activity
                     << ", live_time=" << live_time_seconds
                     << "s, activity_for_clustering=" << activity_for_clustering << endl;

              nucs_and_acts.emplace_back( *nuc_ptr_local, activity_for_clustering );
            }
          }
        }

        // Use the FWHM coefficients passed to this function (computed from auto-search peaks
        // or DRF in fit_peaks_for_nuclides), rather than relying on solution.m_drf which may
        // not have valid FWHM info or may have incorrect values
        
        

        // Get energy range from foreground
        const double lowest_energy_gamma = foreground->gamma_energy_min();
        const double highest_energy_gamma = foreground->gamma_energy_max();

        // Get auto clustering settings from config
        const GammaClusteringSettings auto_settings = config.get_auto_clustering_settings();

        // Cluster gammas using current solution's relative efficiency
        const vector<RelActCalcAuto::RoiRange> refined_rois = cluster_gammas_to_rois(
            auto_rel_eff, nucs_and_acts, foreground,
            fwhm_form, fwhm_coefficients,
            fwhm_lower_energy, fwhm_upper_energy,
            lowest_energy_gamma, highest_energy_gamma,
            auto_settings );

        // Debug output: print refined ROIs with expected counts
        if( PeakFitImprove::debug_printout )
        {
          cout << "Iteration " << iter << " refined ROIs:" << endl;
          for( size_t roi_idx = 0; roi_idx < refined_rois.size(); ++roi_idx )
          {
            const RelActCalcAuto::RoiRange &roi = refined_rois[roi_idx];
            const double roi_expected_counts = foreground->gamma_integral(
                static_cast<float>(roi.lower_energy), static_cast<float>(roi.upper_energy) );
            cout << "  ROI " << roi_idx << ": [" << roi.lower_energy << " - " << roi.upper_energy
                 << "] keV, width=" << (roi.upper_energy - roi.lower_energy)
                 << " keV, expected_counts=" << roi_expected_counts << endl;
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

        const double old_chi2_dof = solution.m_chi2 / std::max( 1.0, static_cast<double>(solution.m_dof) );
        const double new_chi2_dof = refined_solution.m_chi2 / std::max( 1.0, static_cast<double>(refined_solution.m_dof) );

        // Check if chi2/dof improved
        if( new_chi2_dof >= old_chi2_dof )
        {
          cout << "Iteration " << iter << " did not improve chi2/dof ("
               << old_chi2_dof << " -> " << new_chi2_dof << "), stopping" << endl;
          break;
        }

        solution = std::move( refined_solution );
        cout << "Iteration " << iter << " improved: chi2/dof=" << new_chi2_dof
             << " (was " << old_chi2_dof << ")" << endl;
      }//for( size_t iter = 0; iter < max_iterations; ++iter )

      cout << "Final solution after refinement:" << endl;
      solution.print_summary( std::cout );
      cout << endl;
    }//iterative refinement

    // Populate result with successful solution
    result.status = solution.m_status;
    result.error_message = solution.m_error_message;
    result.fit_peaks = solution.m_peaks_without_back_sub;
    result.solution = std::move( solution );

  } catch( const std::exception &e )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    result.error_message = e.what();
  }

  return result;
}//fit_peaks_for_nuclide_relactauto


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
  const std::vector<const SandiaDecay::Nuclide *> &source_nuclides,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const bool isHPGe,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const double lowest_energy_gamma,
  const double highest_energy_gamma,
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

  // Step 2: Estimate activity for each nuclide
  std::vector<std::pair<const SandiaDecay::Nuclide *, double>> nucs_and_acts;

  for( const SandiaDecay::Nuclide * const nuc : source_nuclides )
  {
    if( !nuc )
      continue;

    // Get gamma lines at default age
    const double age = PeakDef::defaultDecayTime( nuc );
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc, 1.0, age );
    const std::vector<SandiaDecay::EnergyRatePair> photons
        = mix.photons( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy );

    // Find gamma with largest expected yield (br * efficiency)
    double best_yield = 0.0;
    double best_energy = 0.0;
    double best_br = 0.0;
    double best_eff = 0.0;

    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      if( photon.energy < lowest_energy_gamma || photon.energy > highest_energy_gamma )
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
      // Use peak area: activity = peak_area / (br * eff * live_time)
      const double peak_area = matched_peak->peakArea();
      estimated_activity = peak_area / (best_br * best_eff * live_time);

      cout << "Fallback: " << nuc->symbol << " matched peak at " << matched_peak->mean()
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

      cout << "Fallback: " << nuc->symbol << " no matching peak for gamma at " << best_energy
           << " keV, integrated counts=" << total_counts << ", est. peak area=" << estimated_peak_area
           << ", estimated activity=" << estimated_activity << endl;
    }

    if( estimated_activity > 0.0 )
      nucs_and_acts.emplace_back( nuc, estimated_activity );
  }

  if( nucs_and_acts.empty() )
  {
    cerr << "Fallback: Could not estimate activity for any nuclide" << endl;
    return {};
  }

  // Step 3: Create relative efficiency function from DRF
  const auto fallback_rel_eff = [&drf_to_use]( double energy ) -> double {
    return drf_to_use->intrinsicEfficiency( static_cast<float>(energy) );
  };

  // Step 4: Call cluster_gammas_to_rois with estimated activities
  return cluster_gammas_to_rois(
    fallback_rel_eff,
    nucs_and_acts,
    foreground,
    fwhmFnctnlForm,
    fwhm_coefficients,
    lower_fwhm_energy,
    upper_fwhm_energy,
    lowest_energy_gamma,
    highest_energy_gamma,
    settings
  );
}//estimate_initial_rois_fallback


PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<const SandiaDecay::Nuclide *> &source_nuclides,
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

  // Validate source nuclides
  if( source_nuclides.empty() )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "No source nuclides provided";
    return result;
  }

  for( const SandiaDecay::Nuclide *nuc : source_nuclides )
  {
    if( !nuc )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Null nuclide in source_nuclides vector";
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

    // Determine energy range for gamma lines
    double highest_energy_gamma = 0.0, lowest_energy_gamma = std::numeric_limits<double>::max();

    vector<RelActCalcAuto::NucInputInfo> base_nuclides;
    for( const SandiaDecay::Nuclide *nuc : source_nuclides )
    {
      RelActCalcAuto::NucInputInfo nuc_info;
      nuc_info.source = nuc;
      nuc_info.age = PeakDef::defaultDecayTime(nuc);
      nuc_info.fit_age = false;
      base_nuclides.push_back(nuc_info);

      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( nuc, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits, nuc_info.age );
      const vector<SandiaDecay::EnergyRatePair> photons = mix.photons( 0.0 );
      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        highest_energy_gamma = (std::max)( highest_energy_gamma, photon.energy );
        lowest_energy_gamma = (std::min)( lowest_energy_gamma, photon.energy );
      }
    }

    lowest_energy_gamma = (std::max)( lowest_energy_gamma - (isHPGe ? 5 : 25), (double)foreground->gamma_energy_min() );
    highest_energy_gamma = (std::min)( highest_energy_gamma + (isHPGe ? 5 : 25), (double)foreground->gamma_energy_max() );

    // Step 2: Match auto-search peaks to source nuclides
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

    vector<RelActCalcManual::PeakCsvInput::NucAndAge> rel_act_manual_srcs;
    for( const SandiaDecay::Nuclide *nuc : source_nuclides )
      rel_act_manual_srcs.emplace_back( nuc->symbol, -1.0, false );

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
      }
    }

    // Get clustering settings from config
    const GammaClusteringSettings manual_settings = config.get_manual_clustering_settings();

    vector<RelActCalcAuto::RoiRange> initial_rois;

    // Step 3 & 4: Perform RelActManual fit and create initial ROIs
    if( !peaks_matched.empty() )
    {
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

      try
      {
        const RelActCalcManual::RelEffSolution manual_solution
            = RelActCalcManual::solve_relative_efficiency( manual_input );

        if( manual_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success )
        {
          cout << "Successfully fitted initial RelActCalcManual::RelEffSolution: chi2/dof="
          << manual_solution.m_chi2 << "/" << manual_solution.m_dof << "="
          << manual_solution.m_chi2 / manual_solution.m_dof
          << " using " << peak_match_results.peaks_matched.size() << " peaks"
          << endl;
          cout << endl;

          // Create rel_eff lambda from manual solution - handles extrapolation clamping
          const auto manual_rel_eff = [&manual_solution, &peaks_matched]( double energy ) -> double {
            // Extrapolation is terrible for rel-eff, so clamp to the lowest/highest peak energy
            if( energy < peaks_matched.front().m_energy )
              return manual_solution.relative_efficiency( peaks_matched.front().m_energy );
            else if( energy > peaks_matched.back().m_energy )
              return manual_solution.relative_efficiency( peaks_matched.back().m_energy );
            else
              return manual_solution.relative_efficiency( energy );
          };

          // Collect nuclides and activities from the manual solution
          vector<pair<const SandiaDecay::Nuclide *, double>> nucs_and_acts;
          for( const RelActCalcManual::IsotopeRelativeActivity &rel_act : manual_solution.m_rel_activities )
          {
            const SandiaDecay::Nuclide * const nuc = db->nuclide( rel_act.m_isotope );
            assert( nuc );
            if( !nuc )
              throw std::logic_error( "Failed to get nuclide from RelAct nuc '" + rel_act.m_isotope + "'" );
            nucs_and_acts.emplace_back( nuc, manual_solution.relative_activity( rel_act.m_isotope ) );
          }

          // Use the reusable clustering function to create ROIs
          initial_rois = cluster_gammas_to_rois( manual_rel_eff, nucs_and_acts, foreground,
                                                 fwhmFnctnlForm, fwhm_coefficients,
                                                 lower_fwhm_energy, upper_fwhm_energy,
                                                 lowest_energy_gamma, highest_energy_gamma,
                                                 manual_settings );

          cout << "Initial ROIs from RelActManual: ";
          for( const auto &roi : initial_rois )
            cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
          cout << endl;
        }
        else
        {
          cout << "Failed to fit initial RelActCalcManual::RelEffSolution: " << manual_solution.m_error_message << endl;
          cout << endl;
        }
      }
      catch( std::exception &e )
      {
        cerr << "Error trying to fit initial manual rel-eff solution: " << e.what() << endl;
        cerr << "Using fallback activity estimation..." << endl;

        initial_rois = estimate_initial_rois_fallback(
          auto_search_peaks, foreground, source_nuclides, drf, isHPGe,
          fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
          lowest_energy_gamma, highest_energy_gamma, manual_settings );

        fallback_warning = "RelActManual fitting failed (" + std::string(e.what())
                         + "); used simplified activity estimation fallback";

        cout << "Fallback ROIs: ";
        for( const RelActCalcAuto::RoiRange &roi : initial_rois )
          cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
        cout << endl;
      }
    }

    // Step 5: Prepare configuration for RelActAuto
    PeakFitForNuclideConfig config_with_nucs = config;
    config_with_nucs.base_nuclides = base_nuclides;
    // initial_rois is now passed as parameter to create_*_options() methods

    // Define skew type to try
    const PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;

    // Step 6: Create RelActAuto options from config
    RelActCalcAuto::Options options;

    // Create relative efficiency curve from config parameters
    RelActCalcAuto::RelEffCurveInput rel_eff_curve;
    rel_eff_curve.rel_eff_eqn_type = config.rel_eff_eqn_type;
    rel_eff_curve.rel_eff_eqn_order = config.rel_eff_eqn_order;

    // FramPhysicalModel requires rel_eff_eqn_order to be 0
    if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      rel_eff_curve.rel_eff_eqn_order = 0;

    rel_eff_curve.nucs_of_el_same_age = config.nucs_of_el_same_age;
    rel_eff_curve.phys_model_use_hoerl = config.phys_model_use_hoerl;
    rel_eff_curve.nuclides = config_with_nucs.base_nuclides;

    // Copy shielding inputs from config (converting to const shared_ptr)
    if( !config.phys_model_self_atten.empty() )
      rel_eff_curve.phys_model_self_atten = config.phys_model_self_atten.front();

    rel_eff_curve.phys_model_external_atten.clear();
    for( const auto &shield : config.phys_model_external_atten )
      rel_eff_curve.phys_model_external_atten.push_back( shield );

    // Generate name based on equation type and shielding
    if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::LnX )
    {
      rel_eff_curve.name = "LnX Peak Fit";
    }
    else if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      rel_eff_curve.name = "Physical Model";
      if( !config.phys_model_external_atten.empty() )
      {
        const int z = static_cast<int>( config.phys_model_external_atten.front()->atomic_number + 0.5 );
        if( z == 13 )
          rel_eff_curve.name += " Al";
        else if( z == 82 )
          rel_eff_curve.name += " Pb";
        else
          rel_eff_curve.name += " Z=" + std::to_string(z);
      }
      rel_eff_curve.name += " Peak Fit";
    }

    options.rel_eff_curves.push_back( rel_eff_curve );
    options.rois = initial_rois;
    options.fit_energy_cal = config.fit_energy_cal;
    options.fwhm_form = config.fwhm_form;
    options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
    options.skew_type = skew_type;
    options.additional_br_uncert = config.rel_eff_auto_base_rel_eff_uncert;

    // Call RelActAuto with the configured options
    result = fit_peaks_for_nuclide_relactauto(
      auto_search_peaks, foreground, source_nuclides,
      options, long_background, drf, config_with_nucs,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy );

  }
  catch( const std::exception &e )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    result.error_message = e.what();
  }

  // Propagate fallback warning if we used it
  if( !fallback_warning.empty() )
    result.warnings.push_back( fallback_warning );

  return result;
}//fit_peaks_for_nuclides


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
    const vector<pair<const SandiaDecay::Nuclide *, double>> &nuclides_and_activities,
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

  for( const auto &nuc_act : nuclides_and_activities )
  {
    const SandiaDecay::Nuclide * const nuc = nuc_act.first;
    const double activity = nuc_act.second;

    if( !nuc || (activity <= 0.0) )
      continue;

    const double age = PeakDef::defaultDecayTime( nuc );

    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc, activity, age );
    const vector<SandiaDecay::EnergyRatePair> photons
        = mix.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );

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
  }//for( const auto &nuc_act : nuclides_and_activities )

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
      continue; // Already removed from gammas_by_energy

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

    const double lower = energy - settings.cluster_num_sigma * sigma;
    const double upper = energy + settings.cluster_num_sigma * sigma;

    // Find all gammas in this range
    const auto start_remove = std::lower_bound( begin(gammas_by_energy), end(gammas_by_energy),
                                               make_pair(lower, 0.0), lessThanByEnergy );
    const auto end_remove = std::upper_bound( begin(gammas_by_energy), end(gammas_by_energy),
                                             make_pair(upper, 0.0), lessThanByEnergy );

    const double counts_in_region = std::accumulate( start_remove, end_remove, 0.0,
        []( const double &sum, const pair<double,double> &el ) {
          return sum + el.second;
    } );

    // Capture the gamma lines before erasing them
    vector<pair<double,double>> gamma_lines_in_cluster( start_remove, end_remove );

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

    if( (data_area > settings.min_data_area_keep)
       && (counts_in_region > settings.min_est_peak_area_keep)
       && (signif > settings.min_est_significance_keep) )
    {
      ClusteredGammaInfo cluster_info;
      cluster_info.lower = lower;
      cluster_info.upper = upper;
      cluster_info.gamma_lines = std::move( gamma_lines_in_cluster );
      clustered_gammas.push_back( std::move( cluster_info ) );
    }
  }//for( const pair<double,double> &energy_counts : gammas_by_counts )

  // Sort by lower energy
  std::sort( begin(clustered_gammas), end(clustered_gammas),
    []( const ClusteredGammaInfo &lhs, const ClusteredGammaInfo &rhs ) {
      return lhs.lower < rhs.lower;
  } );

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
      merged_clusters.back().upper = std::max( merged_clusters.back().upper, cluster.upper );
      merged_clusters.back().gamma_lines.insert(
        end(merged_clusters.back().gamma_lines),
        begin(cluster.gamma_lines),
        end(cluster.gamma_lines)
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

    if( current_width <= max_width )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    // Need to break up this cluster - sort gamma lines by energy first
    std::sort( begin(cluster.gamma_lines), end(cluster.gamma_lines),
      []( const pair<double,double> &a, const pair<double,double> &b ) {
        return a.first < b.first;
    } );

    if( cluster.gamma_lines.size() <= 1 )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    // Calculate potential breakpoints between consecutive gamma lines
    // See algorithm description in the comments near line 607-627
    struct BreakpointInfo {
      double gap_center;
      double score;
      size_t left_index;
    };
    vector<BreakpointInfo> potential_breaks;

    const double min_distance = 0.1; // Minimum distance to avoid division by zero (0.1 keV)

    for( size_t i = 0; i < cluster.gamma_lines.size() - 1; ++i )
    {
      const double left_energy = cluster.gamma_lines[i].first;
      const double right_energy = cluster.gamma_lines[i+1].first;

      const double gap_center = 0.5 * (left_energy + right_energy);
      const double gap_width = right_energy - left_energy;

      // Calculate weighted sum for left side (inverse distance weighting, quadrature sum)
      double left_weight_sq_sum = 0.0;
      for( size_t j = 0; j <= i; ++j )
      {
        const double gamma_energy = cluster.gamma_lines[j].first;
        const double gamma_amp = cluster.gamma_lines[j].second;
        const double distance = std::max( min_distance, std::abs(gap_center - gamma_energy) );
        const double weighted_amp = gamma_amp / distance;
        left_weight_sq_sum += weighted_amp * weighted_amp;
      }
      const double left_weight = std::sqrt( left_weight_sq_sum );

      // Calculate weighted sum for right side
      double right_weight_sq_sum = 0.0;
      for( size_t j = i + 1; j < cluster.gamma_lines.size(); ++j )
      {
        const double gamma_energy = cluster.gamma_lines[j].first;
        const double gamma_amp = cluster.gamma_lines[j].second;
        const double distance = std::max( min_distance, std::abs(gap_center - gamma_energy) );
        const double weighted_amp = gamma_amp / distance;
        right_weight_sq_sum += weighted_amp * weighted_amp;
      }
      const double right_weight = std::sqrt( right_weight_sq_sum );

      const double score = gap_width / (left_weight + right_weight);
      potential_breaks.push_back( { gap_center, score, i } );
    }//for( size_t i = 0; i < cluster.gamma_lines.size() - 1; ++i )

    // Sort by score descending (best breakpoints first)
    std::sort( begin(potential_breaks), end(potential_breaks),
      []( const BreakpointInfo &a, const BreakpointInfo &b ) {
        return a.score > b.score;
    } );

    // Greedily select breakpoints to ensure no sub-ROI exceeds max width
    vector<size_t> selected_break_indices;

    for( const BreakpointInfo &bp : potential_breaks )
    {
      selected_break_indices.push_back( bp.left_index );
      std::sort( begin(selected_break_indices), end(selected_break_indices) );

      // Check if all sub-ROIs are now within max width
      bool all_ok = true;
      double seg_start = cluster.lower;

      for( size_t break_idx : selected_break_indices )
      {
        const double seg_end_energy = 0.5 * (cluster.gamma_lines[break_idx].first
                                            + cluster.gamma_lines[break_idx + 1].first);
        const double seg_mid = 0.5 * (seg_start + seg_end_energy);
        const double seg_mid_clamped = have_fwhm_range
            ? std::clamp( seg_mid, fwhm_lower_energy, fwhm_upper_energy )
            : seg_mid;
        const double seg_fwhm = DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(seg_mid_clamped), fwhm_form, fwhm_coefficients );
        const double seg_max = settings.max_fwhm_width * seg_fwhm;

        if( (seg_end_energy - seg_start) > seg_max )
        {
          all_ok = false;
          break;
        }
        seg_start = seg_end_energy;
      }

      // Check the last segment
      if( all_ok )
      {
        const double seg_mid = 0.5 * (seg_start + cluster.upper);
        const double seg_mid_clamped = have_fwhm_range
            ? std::clamp( seg_mid, fwhm_lower_energy, fwhm_upper_energy )
            : seg_mid;
        const double seg_fwhm = DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(seg_mid_clamped), fwhm_form, fwhm_coefficients );
        const double seg_max = settings.max_fwhm_width * seg_fwhm;

        if( (cluster.upper - seg_start) > seg_max )
          all_ok = false;
      }

      if( all_ok )
        break;
    }//for( const BreakpointInfo &bp : potential_breaks )

    // Create sub-clusters based on selected breakpoints
    if( selected_break_indices.empty() )
    {
      final_clusters.push_back( std::move( cluster ) );
    }
    else
    {
      std::sort( begin(selected_break_indices), end(selected_break_indices) );

      double seg_start = cluster.lower;
      size_t gamma_start_idx = 0;

      for( size_t break_idx : selected_break_indices )
      {
        const double seg_end = 0.5 * (cluster.gamma_lines[break_idx].first
                                     + cluster.gamma_lines[break_idx + 1].first);

        ClusteredGammaInfo sub_cluster;
        sub_cluster.lower = seg_start;
        sub_cluster.upper = seg_end;
        for( size_t j = gamma_start_idx; j <= break_idx; ++j )
          sub_cluster.gamma_lines.push_back( cluster.gamma_lines[j] );

        final_clusters.push_back( std::move( sub_cluster ) );

        seg_start = seg_end;
        gamma_start_idx = break_idx + 1;
      }

      // Add the last segment
      ClusteredGammaInfo sub_cluster;
      sub_cluster.lower = seg_start;
      sub_cluster.upper = cluster.upper;
      for( size_t j = gamma_start_idx; j < cluster.gamma_lines.size(); ++j )
        sub_cluster.gamma_lines.push_back( cluster.gamma_lines[j] );

      final_clusters.push_back( std::move( sub_cluster ) );
    }
  }//for( ClusteredGammaInfo &cluster : merged_clusters )

  // Create ROIs from final clusters
  double previous_roi_upper = 0.0;

  for( const ClusteredGammaInfo &cluster : final_clusters )
  {
    RelActCalcAuto::RoiRange roi;

    if( cluster.gamma_lines.empty() )
    {
      roi.lower_energy = cluster.lower;
      roi.upper_energy = cluster.upper;
    }
    else
    {
      // Find min/max gamma energies in this cluster
      double min_gamma_energy = cluster.gamma_lines.front().first;
      double max_gamma_energy = cluster.gamma_lines.front().first;
      for( const pair<double,double> &gamma : cluster.gamma_lines )
      {
        min_gamma_energy = std::min( min_gamma_energy, gamma.first );
        max_gamma_energy = std::max( max_gamma_energy, gamma.first );
      }

      // Calculate FWHM at the min and max gamma energies (clamp to valid range)
      const double min_gamma_clamped = have_fwhm_range
          ? std::clamp( min_gamma_energy, fwhm_lower_energy, fwhm_upper_energy )
          : min_gamma_energy;
      const double max_gamma_clamped = have_fwhm_range
          ? std::clamp( max_gamma_energy, fwhm_lower_energy, fwhm_upper_energy )
          : max_gamma_energy;

      const double fwhm_at_lower = DetectorPeakResponse::peakResolutionFWHM(
          static_cast<float>(min_gamma_clamped), fwhm_form, fwhm_coefficients );
      const double fwhm_at_upper = DetectorPeakResponse::peakResolutionFWHM(
          static_cast<float>(max_gamma_clamped), fwhm_form, fwhm_coefficients );

      // Calculate desired ROI bounds based on FWHM distance from outer gamma lines
      const double desired_lower = min_gamma_energy - settings.roi_width_num_fwhm_lower * fwhm_at_lower;
      const double desired_upper = max_gamma_energy + settings.roi_width_num_fwhm_upper * fwhm_at_upper;

      // Constrain lower bound to not overlap with previous ROI
      roi.lower_energy = std::max( desired_lower, previous_roi_upper );
      roi.upper_energy = desired_upper;
    }

    const double mid_energy = 0.5 * (roi.upper_energy + roi.lower_energy);
    const double mid_energy_clamped = have_fwhm_range
        ? std::clamp( mid_energy, fwhm_lower_energy, fwhm_upper_energy )
        : mid_energy;
    const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(mid_energy_clamped), fwhm_form, fwhm_coefficients );
    const double num_fwhm_wide = (roi.upper_energy - roi.lower_energy) / mid_fwhm;

    if( num_fwhm_wide < settings.min_fwhm_roi )
      continue;

    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    if( num_fwhm_wide > settings.min_fwhm_quad_cont )
      roi.continuum_type = PeakContinuum::OffsetType::Quadratic;

    result_rois.push_back( roi );
    previous_roi_upper = roi.upper_energy;
  }//for( const ClusteredGammaInfo &cluster : final_clusters )

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
  const std::vector<const SandiaDecay::Nuclide *> &source_nuclides,
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

    vector<const SandiaDecay::Nuclide *> source_nuclides;
    if( src_name == "Tl201woTl202" )
    {
      source_nuclides.push_back( db->nuclide("Tl201"));
    }else if( src_name == "Tl201wTl202" )
    {
      source_nuclides.push_back( db->nuclide("Tl201"));
      source_nuclides.push_back( db->nuclide("Tl202"));
    }else if( src_name == "Uore" )
    {
      source_nuclides.push_back( db->nuclide("U235") );
      source_nuclides.push_back( db->nuclide("U238") );
      source_nuclides.push_back( db->nuclide("Ra226") );
    }else if( src_name == "Am241Li" )
    {
      continue;
    }else
    {
      const SandiaDecay::Nuclide * const nuc = db->nuclide( src_name );
      if( !nuc )
      {
        cout << "Unable to get nuclide for '" << src_name << "' - so will skip" << endl;
        continue;
      }
      source_nuclides.push_back( nuc );
      assert( nuc->symbol == src_name );
    }

    assert( !source_nuclides.empty() && source_nuclides.front() );
    if( source_nuclides.empty() || !source_nuclides.front() )
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

      // Try all three RelEff curve types
      std::vector<PeakFitResult> curve_results;

      // Try LnX
      PeakFitForNuclideConfig config_lnx = config;
      config_lnx.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
      config_lnx.rel_eff_eqn_order = 3;
      config_lnx.phys_model_self_atten.clear();
      config_lnx.phys_model_external_atten.clear();
      curve_results.push_back( fit_peaks_for_nuclides(
        auto_search_peaks, foreground, source_nuclides, long_background, drf, config_lnx, isHPGe ) );

      // Try Aluminum shielding
      PeakFitForNuclideConfig config_aluminum = config;
      config_aluminum.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      config_aluminum.phys_model_use_hoerl = true;
      config_aluminum.phys_model_self_atten.clear();
      config_aluminum.phys_model_external_atten.clear();

      std::shared_ptr<RelActCalc::PhysicalModelShieldInput> al_shield
        = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
      al_shield->atomic_number = 13.0; // Aluminum
      al_shield->fit_atomic_number = false;
      al_shield->areal_density = 0.1 * PhysicalUnits::g_per_cm2; // Starting value in g/cm2
      al_shield->fit_areal_density = true;
      al_shield->lower_fit_areal_density = 0.0;
      al_shield->upper_fit_areal_density = 10.0 * PhysicalUnits::g_per_cm2;
      config_aluminum.phys_model_external_atten.push_back( al_shield );

      curve_results.push_back( fit_peaks_for_nuclides(
        auto_search_peaks, foreground, source_nuclides, long_background, drf, config_aluminum, isHPGe ) );

      // Try Lead shielding
      PeakFitForNuclideConfig config_lead = config;
      config_lead.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      config_lead.phys_model_use_hoerl = true;
      config_lead.phys_model_self_atten.clear();
      config_lead.phys_model_external_atten.clear();

      std::shared_ptr<RelActCalc::PhysicalModelShieldInput> pb_shield
        = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
      pb_shield->atomic_number = 82.0; // Lead
      pb_shield->fit_atomic_number = false;
      pb_shield->areal_density = 0.1 * PhysicalUnits::g_per_cm2; // Starting value in g/cm2
      pb_shield->fit_areal_density = true;
      pb_shield->lower_fit_areal_density = 0.0;
      pb_shield->upper_fit_areal_density = 10.0 * PhysicalUnits::g_per_cm2;
      config_lead.phys_model_external_atten.push_back( pb_shield );

      curve_results.push_back( fit_peaks_for_nuclides(
        auto_search_peaks, foreground, source_nuclides, long_background, drf, config_lead, isHPGe ) );

      // Select best result based on chi2
      PeakFitResult best_result;
      double best_chi2 = std::numeric_limits<double>::max();
      bool found_valid = false;
      std::string last_error;
      const char *curve_names[] = { "LnX", "Aluminum", "Lead" };

      for( size_t i = 0; i < curve_results.size(); ++i )
      {
        const PeakFitResult &curve_result = curve_results[i];

        cout << "RelEff curve type " << curve_names[i] << " status=" << static_cast<int>(curve_result.status);
        if( curve_result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
          cout << " chi2/dof=" << curve_result.solution.m_chi2 << "/" << curve_result.solution.m_dof << endl;
        else
          cout << " error: " << curve_result.error_message << endl;

        if( curve_result.status == RelActCalcAuto::RelActAutoSolution::Status::Success
            && curve_result.solution.m_chi2 < best_chi2 )
        {
          best_chi2 = curve_result.solution.m_chi2;
          best_result = curve_result;
          found_valid = true;
        }
        else if( curve_result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        {
          last_error = curve_result.error_message;
        }
      }

      if( !found_valid )
      {
        throw std::runtime_error( "fit_peaks_for_nuclides failed for all RelEff curve types: " + last_error );
      }

      const RelActCalcAuto::RelActAutoSolution &solution = best_result.solution;
      const std::vector<PeakDef> &fit_peaks = best_result.fit_peaks;

      cout << "Best RelEff curve type solution (" << solution.m_options.rel_eff_curves.front().name << ") selected with chi2/dof="
           << solution.m_chi2 << "/" << solution.m_dof << endl;

      // Old manually-inlined code was here - now using fit_peaks_for_nuclides()
      /*
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
          = make_shared<const deque<shared_ptr<const PeakDef>>>( begin(auto_search_peaks) , end(auto_search_peaks) );

        MakeDrfFit::performResolutionFit( auto_search_peaks_dq, fwhmFnctnlForm, sqrtEqnOrder, fwhm_coefficients, fwhm_uncerts );
        auto_search_peaks_dq = MakeDrfFit::removeOutlyingWidthPeaks( auto_search_peaks_dq, fwhmFnctnlForm, fwhm_coefficients );
        MakeDrfFit::performResolutionFit( auto_search_peaks_dq, fwhmFnctnlForm, sqrtEqnOrder, fwhm_coefficients, fwhm_uncerts );

        // Set energy range based on peaks used for FWHM fit
        if( !auto_search_peaks_dq->empty() )
        {
          lower_fwhm_energy = auto_search_peaks_dq->front()->mean();
          upper_fwhm_energy = auto_search_peaks_dq->back()->mean();
        }
      }else
      {
        fwhmFnctnlForm = drf->resolutionFcnType();
        fwhm_coefficients = drf->resolutionFcnCoefficients();

        // Get energy range from detector response function
        lower_fwhm_energy = drf->lowerEnergy();
        upper_fwhm_energy = drf->upperEnergy();
      }//if( we will fit FWHM functional form from auto-fit peaks ) / else (use detector eff FWHM )

      // Validate that we have valid FWHM coefficients
      if( fwhm_coefficients.empty() )
      {
        throw std::runtime_error( "Failed to determine FWHM coefficients - unable to proceed with peak fitting" );
      }

      // Check that coefficients are finite
      for( size_t i = 0; i < fwhm_coefficients.size(); ++i )
      {
        if( !std::isfinite( fwhm_coefficients[i] ) )
        {
          throw std::runtime_error( "FWHM coefficient[" + std::to_string(i) + "] is not finite (value="
                                   + std::to_string(fwhm_coefficients[i]) + ")" );
        }
      }

      if( PeakFitImprove::debug_printout )
      {
        cout << "FWHM function form: " << static_cast<int>(fwhmFnctnlForm) << ", coefficients: [";
        for( size_t i = 0; i < fwhm_coefficients.size(); ++i )
          cout << (i > 0 ? ", " : "") << fwhm_coefficients[i];
        cout << "]" << endl;
      }

      double highest_energy_gamma = 0.0, lowest_energy_gamma = std::numeric_limits<double>::max();
      
      vector<RelActCalcAuto::NucInputInfo> base_nuclides;
      for( const SandiaDecay::Nuclide *nuc : source_nuclides )
      {
        RelActCalcAuto::NucInputInfo nuc_info;
        nuc_info.source = nuc;
        nuc_info.age = PeakDef::defaultDecayTime(nuc);
        nuc_info.fit_age = false; // Don't fit age by default
        base_nuclides.push_back(nuc_info);
        
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nuc, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits, nuc_info.age );
        const vector<SandiaDecay::EnergyRatePair> photons = mix.photons( 0.0 );
        for( const SandiaDecay::EnergyRatePair &photon : photons )
        {
          highest_energy_gamma = (std::max)( highest_energy_gamma, photon.energy );
          lowest_energy_gamma = (std::min)( lowest_energy_gamma, photon.energy );
        }
      }//for( const SandiaDecay::Nuclide *nuc : source_nuclides )

      lowest_energy_gamma = (std::max)( lowest_energy_gamma - (isHPGe ? 5 : 25), (double)foreground->gamma_energy_min() );
      highest_energy_gamma = (std::min)( highest_energy_gamma + (isHPGe ? 5 : 25), (double)foreground->gamma_energy_max() );
    
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
        // peak_info.m_source_gammas = std::vector<GenericLineInfo>{}; // `fill_in_nuclide_info(...)` will fill this in
        
        rel_act_manual_peaks.push_back( peak_info );
      }//for( const shared_ptr<const PeakDef> &peak : auto_search_peaks )
      
      vector<RelActCalcManual::PeakCsvInput::NucAndAge> rel_act_manual_srcs;
      for( const SandiaDecay::Nuclide *nuc : source_nuclides )
        rel_act_manual_srcs.emplace_back( nuc->symbol, -1.0, false );
      
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
      //const vector<RelActCalcManual::GenericPeakInfo> &peaks_not_matched = peak_match_results.peaks_not_matched;
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
        }
      }//if( PeakFitImprove::debug_printout )
      
      
      vector<RelActCalcAuto::RoiRange> initial_rois;
      
      // If we matched any peaks, we will do a manual rel-eff fit to roughly get the scale of what peaks we should look for
      if( !peaks_matched.empty() )
      {
        RelActCalcManual::RelEffInput manual_input;
        manual_input.peaks = peaks_matched;
        
        manual_input.eqn_order = 0;
        manual_input.use_ceres_to_fit_eqn = false;
        manual_input.phys_model_use_hoerl = false;
        if( peaks_matched.size() == 1 )
        {
          manual_input.eqn_order = config.initial_manual_relEff_1peak_eqn_order;
          manual_input.eqn_form = config.initial_manual_relEff_1peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
          }
        }else if( peaks_matched.size() == 2 )
        {
          manual_input.eqn_order = config.initial_manual_relEff_2peak_eqn_order;
          manual_input.eqn_form = config.initial_manual_relEff_2peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
          }
        }else if( peaks_matched.size() == 3 )
        {
          manual_input.eqn_order = config.initial_manual_relEff_3peak_eqn_order;
          manual_input.eqn_form = config.initial_manual_relEff_3peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
          }
        }else if( peaks_matched.size() == 4 )
        {
          manual_input.eqn_order = config.initial_manual_relEff_4peak_eqn_order;
          manual_input.eqn_form = config.initial_manual_relEff_4peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
            manual_input.phys_model_use_hoerl = config.initial_manual_relEff_4peak_physical_use_hoerl;
          }
        }else
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
            
            // Other generic detectors we could use:
            //DetectorPeakResponse::getGenericLaBrDetector();
            //DetectorPeakResponse::getGenericCZTGeneralDetector();
            //DetectorPeakResponse::getGenericCZTGoodDetector();
          }//if( !manual_input.phys_model_detector )
          
          manual_input.phys_model_self_atten = shared_ptr<const RelActCalc::PhysicalModelShieldInput>{};
          manual_input.phys_model_external_attens = vector<shared_ptr<const RelActCalc::PhysicalModelShieldInput>>{};
        }//if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        
        try
        {
          const RelActCalcManual::RelEffSolution manual_solution
              = RelActCalcManual::solve_relative_efficiency( manual_input );

          if( manual_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success )
          {
            cout << "Successfully fitted initial RelActCalcManual::RelEffSolution: chi2/dof="
            << manual_solution.m_chi2 << "/" << manual_solution.m_dof << "="
            << manual_solution.m_chi2 / manual_solution.m_dof
            << " using " << peak_match_results.peaks_matched.size() << " peaks"
            << endl;
            cout << endl;

            // Create rel_eff lambda from manual solution - handles extrapolation clamping
            const auto manual_rel_eff = [&manual_solution, &peaks_matched]( double energy ) -> double {
              // Extrapolation is terrible for rel-eff, so clamp to the lowest/highest peak energy
              if( energy < peaks_matched.front().m_energy )
                return manual_solution.relative_efficiency( peaks_matched.front().m_energy );
              else if( energy > peaks_matched.back().m_energy )
                return manual_solution.relative_efficiency( peaks_matched.back().m_energy );
              else
                return manual_solution.relative_efficiency( energy );
            };

            // Collect nuclides and activities from the manual solution
            vector<pair<const SandiaDecay::Nuclide *, double>> nucs_and_acts;
            for( const RelActCalcManual::IsotopeRelativeActivity &rel_act : manual_solution.m_rel_activities )
            {
              const SandiaDecay::Nuclide * const nuc = db->nuclide( rel_act.m_isotope );
              assert( nuc );
              if( !nuc )
                throw std::logic_error( "Failed to get nuclide from RelAct nuc '" + rel_act.m_isotope + "'" );
              nucs_and_acts.emplace_back( nuc, manual_solution.relative_activity( rel_act.m_isotope ) );
            }

            // Use the reusable clustering function to create ROIs
            initial_rois = cluster_gammas_to_rois( manual_rel_eff, nucs_and_acts, foreground,
                                                   fwhmFnctnlForm, fwhm_coefficients,
                                                   lower_fwhm_energy, upper_fwhm_energy,
                                                   lowest_energy_gamma, highest_energy_gamma,
                                                   manual_settings );

            cout << "Initial ROIs from RelActManual: ";
            for( const auto &roi : initial_rois )
              cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
            cout << endl;
          }else
          {
            cout << "Failed to fit initial RelActCalcManual::RelEffSolution: " << manual_solution.m_error_message << endl;
            cout << endl;
          }
        }catch( std::exception &e )
        {
          cerr << "Error trying to fit initial manual rel-eff solution: " << e.what() << endl;
          cerr << endl;
        }
      }//if( !peaks_matched.empty() )

      // Populate config fields needed for options creation
      //  Convert source_nuclides to NucInputInfo for base_nuclides
      config.base_nuclides.clear();
      for( const SandiaDecay::Nuclide *nuc : source_nuclides )
      {
        RelActCalcAuto::NucInputInfo nuc_info;
        nuc_info.source = nuc;
        nuc_info.age = PeakDef::defaultDecayTime( nuc );
        nuc_info.fit_age = false;
        config.base_nuclides.push_back( nuc_info );
      }
      // initial_rois is now passed as parameter to create_*_options() methods

      // Create source variant for the first nuclide
      // (for now we only support single-nuclide sources)
      RelActCalcAuto::SrcVariant source_variant = source_nuclides.front();

      // Define skew types to try
      std::vector<PeakDef::SkewType> skew_types = { PeakDef::SkewType::NoSkew };

      // Try three configurations and collect results
      std::vector<PeakFitResult> results;
      for( PeakDef::SkewType skew_type : skew_types )
      {
        results.push_back( fit_peaks_for_nuclide_relactauto(
          auto_search_peaks, foreground, source_variant,
          config.create_lnx_options( skew_type, initial_rois ), long_background, drf, config ) );

        results.push_back( fit_peaks_for_nuclide_relactauto(
          auto_search_peaks, foreground, source_variant,
          config.create_aluminum_options( skew_type, initial_rois ), long_background, drf, config ) );

        results.push_back( fit_peaks_for_nuclide_relactauto(
          auto_search_peaks, foreground, source_variant,
          config.create_lead_options( skew_type, initial_rois ), long_background, drf, config ) );
      }

      // Find best solution by chi2
      PeakFitResult best_result;
      double best_chi2 = std::numeric_limits<double>::max();
      bool found_valid = false;
      std::string last_error;

      for( size_t i = 0; i < results.size(); ++i )
      {
        const PeakFitResult &result = results[i];

        cout << "RelActAuto trial " << i << " status=" << static_cast<int>(result.status);
        if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
          cout << " chi2/dof=" << result.solution.m_chi2 << "/" << result.solution.m_dof << endl;
        else
          cout << " error: " << result.error_message << endl;

        if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success
            && result.solution.m_chi2 < best_chi2 )
        {
          best_chi2 = result.solution.m_chi2;
          best_result = result;
          found_valid = true;
        } else if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        {
          last_error = result.error_message;
        }
      }

      if( !found_valid )
        throw std::runtime_error( "RelActAuto failed for all configs: " + last_error );

      // Use best result (OLD CODE - now handled above)
      const RelActCalcAuto::RelActAutoSolution &solution = best_result.solution;
      const std::vector<PeakDef> &fit_peaks = best_result.fit_peaks;

      cout << "Best solution (" << solution.m_options.rel_eff_curves.front().name << ") selected" << endl;
      */
      // End of old code

      // Score the fit results using only signal photopeaks
      FinalFit_GA::FinalFitScore fit_score;

      for( const PeakDef &found_peak : fit_peaks )
      {
        const double found_energy = found_peak.mean();
        const double found_fwhm = found_peak.fwhm();
        const double peak_lower_contrib = found_energy - config.num_sigma_contribution * found_peak.sigma();
        const double peak_upper_contrib = found_energy + config.num_sigma_contribution * found_peak.sigma();

        // Find the nearest expected peak from signal-only photopeaks
        const ExpectedPhotopeakInfo *nearest_signal_peak = nullptr;

        for( const ExpectedPhotopeakInfo &expected_peak : src_info.expected_signal_photopeaks )
        {
          const double expected_energy = expected_peak.effective_energy;

          if( (expected_energy >= peak_lower_contrib) && (expected_energy <= peak_upper_contrib) )
          {
            if( !nearest_signal_peak || (fabs(expected_energy - found_energy) < fabs(nearest_signal_peak->effective_energy - found_energy)) )
              nearest_signal_peak = &expected_peak;
          }
        }//for( const ExpectedPhotopeakInfo &expected_peak : src_info.expected_signal_photopeaks )

        if( nearest_signal_peak )
        {
          fit_score.num_peaks_used += 1;

          const double expected_energy = nearest_signal_peak->effective_energy;
          const double expected_fwhm = nearest_signal_peak->effective_fwhm;
          const double expected_sigma = expected_fwhm / 2.35482;
          const double expected_area = nearest_signal_peak->peak_area;

          const double area_diff = fabs( found_peak.amplitude() - expected_area );
          const double area_score = area_diff / sqrt( (expected_area < 1.0) ? 1.0 : expected_area );
          const double width_score = fabs( expected_fwhm - found_fwhm ) / expected_sigma;
          const double position_score = fabs( found_energy - expected_energy ) / expected_sigma;

          fit_score.area_score += std::min( area_score, 20.0 );
          fit_score.width_score += std::min( width_score, 1.0 );
          fit_score.position_score += std::min( position_score, 1.5 );
        }else
        {
          // Found an extra peak we didn't expect
          if( found_peak.amplitude() < 1.0 ) //Not a real peak, so ignore it
            continue;

          const double sqrt_area = sqrt( found_peak.amplitude() ); //The uncertainty that comes from the RelEff fit may be way small
          const double area_uncert = found_peak.amplitudeUncert() > 0.0 ? (std::max)(found_peak.amplitudeUncert(), sqrt_area) : sqrt_area;

          // If this is a significant peak that we didn't expect, penalize - the 8.0 is arbitrary - and we dont really ever expect to encounter having to assign this penalty
          if( found_peak.amplitude() > 8.0 * area_uncert )
          {
            if( PeakFitImprove::debug_printout )
              cout << "Found unexpected peak {mean: " << found_peak.mean()
                   << ", fwhm: " << found_peak.fwhm()
                   << ", area: " << found_peak.amplitude()
                   << ", area_uncert: " << found_peak.amplitudeUncert() << "}" << endl;
          }

          fit_score.ignored_unexpected_peaks += 1;
          fit_score.unexpected_peaks_sum_significance += std::min( 7.5, found_peak.amplitude() / area_uncert ); //Cap at 7.5 sigma (arbitrary)
        }//if( nearest_signal_peak ) / else
      }//for( const PeakDef &found_peak : fit_peaks )

      // Calculate total weight for scoring
      if( fit_score.num_peaks_used <= 1 )
        fit_score.total_weight = fit_score.area_score;
      else
        fit_score.total_weight = fit_score.area_score / fit_score.num_peaks_used;

      // Add to cumulative total score
      total_score += fit_score.total_weight;

      if( PeakFitImprove::debug_printout )
      {
        cout << "Fit score for " << src_name << ":" << endl;
        cout << fit_score.print( "fit_score" ) << endl;
      }

      // Write N42 file with foreground, background, and fit peaks (only if fit succeeded)
      if( best_result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
      {
        const bool write_n42 = true;
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
        const string out_n42 = SpecUtils::append_path( outdir, src_name ) + "_relactauto_fit.n42";

        SpecMeas output;

        output.add_remark( fit_score.print( "fit_score" ) );
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

        output.save2012N42File( out_n42, [=](){
          cerr << "Failed to write '" << out_n42 << "'" << endl;
        });

        if( PeakFitImprove::debug_printout )
          cout << "Wrote N42 file: " << out_n42 << endl;
        }//if( write_n42 )
      }//if( best_result.status == Success )

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
}//void eval_peaks_for_nuclide( const DataSrcInfo &src_info )


}//namespace FitPeaksForNuclideDev
