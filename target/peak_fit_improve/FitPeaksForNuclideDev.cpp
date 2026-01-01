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

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDists.h"
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

  // Peaks after combining overlapping peaks within ROIs.
  // Peaks that are close together (within 1.5 sigma) or where a smaller peak
  // is not statistically distinguishable from a larger peak's tail are merged
  // into a single peak with combined properties.
  std::vector<PeakDef> fit_peaks;

  // Original uncombined peaks from the fit - preserves all individual peak information.
  // This is the raw output from RelActAuto before peak combination.
  std::vector<PeakDef> uncombined_fit_peaks;

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

  }catch( const std::exception &e )
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
      // Use peak area: activity = peak_area / (br * eff) - we wont divide by live_time, to be consistent with `cluster_gammas_to_rois(...)` convention
      const double peak_area = matched_peak->peakArea();
      estimated_activity = peak_area / (best_br * best_eff);

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

        if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
          throw std::runtime_error( "Failed to fit initial RelActCalcManual::RelEffSolution: " + manual_solution.m_error_message );
        
        
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
      }catch( std::exception &e )
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
    assert( (config.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel)
           || (rel_eff_curve.rel_eff_eqn_order == 0) );
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
    if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
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
    }else
    {
      rel_eff_curve.name = RelActCalc::to_str(config.rel_eff_eqn_type) + string(" Peak Fit");
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

    if( PeakFitImprove::debug_printout )
    {
      cerr << "cluster_gammas_to_rois: Nuclide " << nuc->symbol << ", activity=" << activity
           << ", age=" << (age / PhysicalUnits::second) << " seconds ("
           << (age / PhysicalUnits::year) << " years)" << endl;
    }

    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc, activity, age );
    const vector<SandiaDecay::EnergyRatePair> photons
        = mix.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );

    if( PeakFitImprove::debug_printout )
    {
      cerr << "  " << photons.size() << " photons from " << nuc->symbol << ", energy range ["
           << lowest_energy << ", " << highest_energy << "] keV" << endl;
      // Show photons near 807 keV
      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        if( photon.energy >= 800.0 && photon.energy <= 820.0 )
        {
          const double rel_eff = rel_eff_fcn( photon.energy );
          cerr << "    *** Photon near 807 keV: " << photon.energy << " keV, rate=" << photon.numPerSecond
               << " /s, rel_eff=" << rel_eff << ", est_counts=" << (photon.numPerSecond * rel_eff) << endl;
        }
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

  if( PeakFitImprove::debug_printout )
  {
    cerr << "cluster_gammas_to_rois: Input gammas (" << gammas_by_counts.size() << " total):" << endl;
    for( const auto &gc : gammas_by_energy )
    {
      // Show all gammas, but highlight those in 800-820 keV range
      if( gc.first >= 800.0 && gc.first <= 820.0 )
        cerr << "  *** " << gc.first << " keV, est_counts=" << gc.second << " ***" << endl;
    }
    // Also show top 10 by counts
    cerr << "  Top 10 by expected counts:" << endl;
    for( size_t i = 0; i < std::min( gammas_by_counts.size(), size_t(10) ); ++i )
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
    else if( PeakFitImprove::debug_printout )
    {
      // Print why this cluster was rejected
      cerr << "cluster_gammas_to_rois: Rejected cluster [" << lower << ", " << upper << "] keV (energy="
           << energy << " keV): ";
      if( !passes_data_area )
        cerr << "data_area=" << data_area << " < " << settings.min_data_area_keep << "; ";
      if( !passes_counts )
        cerr << "est_counts=" << counts_in_region << " < " << settings.min_est_peak_area_keep << "; ";
      if( !passes_signif )
        cerr << "signif=" << signif << " < " << settings.min_est_significance_keep << "; ";
      cerr << endl;
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
    cluster.lower = weighted_mean - settings.roi_width_num_fwhm_lower * effective_fwhm;
    cluster.upper = weighted_mean + settings.roi_width_num_fwhm_upper * effective_fwhm;
    
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

      merged_clusters.back().upper = std::max( merged_clusters.back().upper, cluster.upper );
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
        sub_cluster.lower = seg_start;
        sub_cluster.upper = bp_energy;

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
      sub_cluster.lower = seg_start;
      sub_cluster.upper = cluster.upper;
      for( size_t j = gamma_start_idx; j < cluster.gamma_energies.size(); ++j )
      {
        sub_cluster.gamma_energies.push_back( cluster.gamma_energies[j] );
        sub_cluster.gamma_amplitudes.push_back( cluster.gamma_amplitudes[j] );
      }

      if( !sub_cluster.gamma_energies.empty() )
        final_clusters.push_back( std::move( sub_cluster ) );
    }
  }//for( ClusteredGammaInfo &cluster : merged_clusters )

  

  // Create ROIs from final clusters
  double previous_roi_upper = 0.0;

  for( const ClusteredGammaInfo &cluster : final_clusters )
  {
    RelActCalcAuto::RoiRange roi;

    // Use the pre-calculated bounds from cluster (based on weighted mean and effective FWHM)
    // Constrain lower bound to not overlap with previous ROI
    roi.lower_energy = std::max( cluster.lower, previous_roi_upper );
    roi.upper_energy = cluster.upper;

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
        cout << "Unable to get nuclide for '" << src.src_name << "' - so will skip" << endl;
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
      config_lnx.rel_eff_eqn_type = config.rel_eff_eqn_type;
      config_lnx.rel_eff_eqn_order = config.rel_eff_eqn_order;
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

      //curve_results.push_back( fit_peaks_for_nuclides(
      //  auto_search_peaks, foreground, source_nuclides, long_background, drf, config_aluminum, isHPGe ) );

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

      //curve_results.push_back( fit_peaks_for_nuclides(
      //  auto_search_peaks, foreground, source_nuclides, long_background, drf, config_lead, isHPGe ) );

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
        cout << "Fit score for " << src.src_name << ":" << endl;
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
        const string out_n42 = SpecUtils::append_path( outdir, src.src_name ) + "_relactauto_fit.n42";

        SpecMeas output;

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
