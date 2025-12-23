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

    const float fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(energy), fwhm_form, fwhm_coefficients );
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

    const double signif = (data_area > 0.0) ? (counts_in_region / sqrt(data_area)) : 0.0;

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

  // Break up ROIs that are too wide
  vector<ClusteredGammaInfo> final_clusters;
  for( ClusteredGammaInfo &cluster : merged_clusters )
  {
    const double mid_energy = 0.5 * (cluster.lower + cluster.upper);
    const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(mid_energy), fwhm_form, fwhm_coefficients );
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
        const double seg_fwhm = DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(seg_mid), fwhm_form, fwhm_coefficients );
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
        const double seg_fwhm = DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(seg_mid), fwhm_form, fwhm_coefficients );
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

      // Calculate FWHM at the min and max gamma energies
      const double fwhm_at_lower = DetectorPeakResponse::peakResolutionFWHM(
          static_cast<float>(min_gamma_energy), fwhm_form, fwhm_coefficients );
      const double fwhm_at_upper = DetectorPeakResponse::peakResolutionFWHM(
          static_cast<float>(max_gamma_energy), fwhm_form, fwhm_coefficients );

      // Calculate desired ROI bounds based on FWHM distance from outer gamma lines
      const double desired_lower = min_gamma_energy - settings.roi_width_num_fwhm_lower * fwhm_at_lower;
      const double desired_upper = max_gamma_energy + settings.roi_width_num_fwhm_upper * fwhm_at_upper;

      // Constrain lower bound to not overlap with previous ROI
      roi.lower_energy = std::max( desired_lower, previous_roi_upper );
      roi.upper_energy = desired_upper;
    }

    const double mid_energy = 0.5 * (roi.upper_energy + roi.lower_energy);
    const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(mid_energy), fwhm_form, fwhm_coefficients );
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


// This is development playground for `AnalystChecks::fit_peaks_for_nuclides(...)`
void eval_peaks_for_nuclide( const std::vector<DataSrcInfo> &srcs_info )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    throw runtime_error( "Failed to open SandiaDecayDataBase" );

  
  // Parameters we will try to optimize
  DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
  const double rel_eff_manual_base_rel_eff_uncert = 0.0;
  const double initial_nuc_match_cluster_num_sigma = 1.5;
  const double manual_eff_cluster_num_sigma = 1.5;
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
  const double manual_rel_eff_sol_min_data_area_keep = 10;
  const double manual_rel_eff_sol_min_est_peak_area_keep = 5.0;
  const double manual_rel_eff_sol_min_est_significance_keep = 0.5;
  const double manual_rel_eff_sol_min_fwhm_roi = 1.0;
  const double manual_rel_eff_sol_min_fwhm_quad_cont = 8;
  const double manual_rel_eff_sol_max_fwhm = 15.0;
  const double manual_rel_eff_roi_width_num_fwhm_lower = 3.0;
  const double manual_rel_eff_roi_width_num_fwhm_upper = 3.0;
  // TODO: lets better assign FWHM form that should be used
  RelActCalcAuto::FwhmForm fwhm_form = RelActCalcAuto::FwhmForm::Berstein_3;
  
  // TODO: add an enum for how to handle existing ROIs here {ignore, replace_source, no_overlap}
  
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
      }else
      {
        fwhmFnctnlForm = drf->resolutionFcnType();
        fwhm_coefficients = drf->resolutionFcnCoefficients();
      }//if( we will fit FWHM functional form from auto-fit peaks ) / else (use detector eff FWHM )
      
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



      // TODO: Implement the ROI breaking up logic here that `RelActCalcAuto::solve` and `RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres` will eventually be udpated to
      //throw runtime_error( "TODO: Implement the ROI breaking up logic here that `RelActCalcAuto::solve` and `RelActCalcAutoImp::RelActAutoCostFcn::solve_ceres` will eventually be udpated to" );
      // I think first round should be based on peaks in the spectrum...
      // Next round use RelEff fit
      
    
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
        peak_info.m_base_rel_eff_uncert = rel_eff_manual_base_rel_eff_uncert;
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
                                                              initial_nuc_match_cluster_num_sigma, excluded_peak_energies,
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
      
      // If we matched any peaks, we will do a maual rel-eff fit to roughly get the scale of what peaks we should look for
      if( !peaks_matched.empty() )
      {
        RelActCalcManual::RelEffInput manual_input;
        manual_input.peaks = peaks_matched;
        
        manual_input.eqn_order = 0;
        manual_input.use_ceres_to_fit_eqn = false;
        manual_input.phys_model_use_hoerl = false;
        if( peaks_matched.size() == 1 )
        {
          manual_input.eqn_order = initial_manual_relEff_1peak_eqn_order;
          manual_input.eqn_form = initial_manual_relEff_1peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
          }
        }else if( peaks_matched.size() == 2 )
        {
          manual_input.eqn_order = initial_manual_relEff_2peak_eqn_order;
          manual_input.eqn_form = initial_manual_relEff_2peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
          }
        }else if( peaks_matched.size() == 3 )
        {
          manual_input.eqn_order = initial_manual_relEff_3peak_eqn_order;
          manual_input.eqn_form = initial_manual_relEff_3peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
          }
        }else if( peaks_matched.size() == 4 )
        {
          manual_input.eqn_order = initial_manual_relEff_4peak_eqn_order;
          manual_input.eqn_form = initial_manual_relEff_4peak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
            manual_input.phys_model_use_hoerl = initial_manual_relEff_4peak_physical_use_hoerl = false;
          }
        }else
        {
          assert( peaks_matched.size() > 4 );
          manual_input.eqn_order = initial_manual_relEff_many_peak_eqn_order;
          manual_input.eqn_form = initial_manual_relEff_manypeak_form;
          if( manual_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
          {
            manual_input.eqn_order = 0;
            manual_input.use_ceres_to_fit_eqn = true;
            manual_input.phys_model_use_hoerl = initial_manual_relEff_many_peak_physical_use_hoerl;
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

            // Create settings for the manual RelEff clustering stage
            GammaClusteringSettings manual_settings;
            manual_settings.cluster_num_sigma = manual_eff_cluster_num_sigma;
            manual_settings.min_data_area_keep = manual_rel_eff_sol_min_data_area_keep;
            manual_settings.min_est_peak_area_keep = manual_rel_eff_sol_min_est_peak_area_keep;
            manual_settings.min_est_significance_keep = manual_rel_eff_sol_min_est_significance_keep;
            manual_settings.roi_width_num_fwhm_lower = manual_rel_eff_roi_width_num_fwhm_lower;
            manual_settings.roi_width_num_fwhm_upper = manual_rel_eff_roi_width_num_fwhm_upper;
            manual_settings.max_fwhm_width = manual_rel_eff_sol_max_fwhm;
            manual_settings.min_fwhm_roi = manual_rel_eff_sol_min_fwhm_roi;
            manual_settings.min_fwhm_quad_cont = manual_rel_eff_sol_min_fwhm_quad_cont;

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

      // Create all trial options combinations to run in parallel
      std::vector<RelActCalcAuto::Options> trial_options;

      // Define skew types to try
      std::vector<PeakDef::SkewType> skew_types = {
        PeakDef::SkewType::NoSkew,
        //PeakDef::SkewType::GaussExp
      };

      for (PeakDef::SkewType skew_type : skew_types) {
        // LnX equation form
        {
          RelActCalcAuto::Options solve_options;

          // Create relative efficiency curve with LnX
          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "LnX Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
          rel_eff_curve.rel_eff_eqn_order = 3;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;

          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois = initial_rois;

          // Set other options
          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = fwhm_form;
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;

          trial_options.push_back(std::move(solve_options));
        }

        // FramPhysicalModel with Aluminum shielding
        {
          RelActCalcAuto::Options solve_options;

          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "Physical Model Al Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;
          rel_eff_curve.phys_model_use_hoerl = true;

          // Create aluminum external shielding
          auto al_shield = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
          al_shield->atomic_number = 13.0; // Aluminum
          al_shield->fit_atomic_number = false;
          al_shield->areal_density = 0.1; // Starting value in g/cm2
          al_shield->fit_areal_density = true;
          al_shield->lower_fit_areal_density = 0.001;
          al_shield->upper_fit_areal_density = 10.0;

          rel_eff_curve.phys_model_external_atten.push_back(al_shield);

          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois = initial_rois;

          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = fwhm_form;
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;

          trial_options.push_back(std::move(solve_options));
        }

        // FramPhysicalModel with Lead shielding
        {
          RelActCalcAuto::Options solve_options;

          RelActCalcAuto::RelEffCurveInput rel_eff_curve;
          rel_eff_curve.name = "Physical Model Pb Peak Fit";
          rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
          rel_eff_curve.nucs_of_el_same_age = true;
          rel_eff_curve.nuclides = base_nuclides;
          rel_eff_curve.phys_model_use_hoerl = true;

          // Create lead external shielding
          auto pb_shield = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
          pb_shield->atomic_number = 82.0; // Lead
          pb_shield->fit_atomic_number = false;
          pb_shield->areal_density = 0.1; // Starting value in g/cm2
          pb_shield->fit_areal_density = true;
          pb_shield->lower_fit_areal_density = 0.0;
          pb_shield->upper_fit_areal_density = 10.0;

          rel_eff_curve.phys_model_external_atten.push_back(pb_shield);

          solve_options.rel_eff_curves.push_back(rel_eff_curve);
          solve_options.rois = initial_rois;

          solve_options.fit_energy_cal = true;
          solve_options.fwhm_form = fwhm_form; //RelActCalcAuto::FwhmForm::Gadras;
          solve_options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
          solve_options.skew_type = skew_type;
          solve_options.additional_br_uncert = 0.0;

          trial_options.push_back(std::move(solve_options));
        }
      }

      // Execute all trials in parallel using ThreadPool
      std::vector<RelActCalcAuto::RelActAutoSolution> solutions(trial_options.size());
      std::mutex solutions_mutex;

      SpecUtilsAsync::ThreadPool pool;
      for (size_t i = 0; i < trial_options.size(); ++i) {
        pool.post([&, i]() {
          try {
            RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
              trial_options[i], foreground, long_background, drf, auto_search_peaks
            );

            std::lock_guard<std::mutex> lock(solutions_mutex);
            solutions[i] = std::move(solution);
          } catch (const std::exception &e) {
            std::lock_guard<std::mutex> lock(solutions_mutex);
            solutions[i].m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
            solutions[i].m_error_message = e.what();
          }
        });
      }

      pool.join();

      // Find the best solution from all parallel trials
      RelActCalcAuto::RelActAutoSolution best_solution;
      double best_chi2 = std::numeric_limits<double>::max();
      bool found_valid_solution = false;
      string last_error_message;

      for (size_t i = 0; i < solutions.size(); ++i) {
        const RelActCalcAuto::RelActAutoSolution &solution = solutions[i];

        cout << "Initial RelActAuto solution " << i << " chi2/dof=" << solution.m_chi2 << "/" << solution.m_dof
        << "=" << solution.m_chi2 / solution.m_dof << endl;
        
        if (solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success &&
            solution.m_chi2 < best_chi2) {
          best_chi2 = solution.m_chi2;
          best_solution = solution;
          found_valid_solution = true;
        } else if (solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success) {
          last_error_message = solution.m_error_message;
          cerr << "RelActCalcAuto::solve(trial " << i << "): failed: " << solution.m_error_message << endl;
        }
      }

      if (!found_valid_solution) {
        throw std::runtime_error("RelActCalcAuto::solve failed for all attempted configurations ('" + last_error_message + "')");
      }

      // Use the best solution
      RelActCalcAuto::RelActAutoSolution solution = std::move(best_solution);

      cout << "Best initial solution:" << endl;
      solution.print_summary( std::cout );
      cout << endl;

      // Iteratively refine ROIs using RelActAuto solutions
      // The idea is that each iteration provides a better relative efficiency estimate,
      // which allows us to better identify significant gamma lines and create better ROIs.
      {
        // Settings for the auto refinement stage - may want different values than manual stage
        GammaClusteringSettings auto_settings;
        auto_settings.cluster_num_sigma = 2.0;  // Slightly wider clustering with better rel-eff
        auto_settings.min_data_area_keep = 10;
        auto_settings.min_est_peak_area_keep = 5.0;
        auto_settings.min_est_significance_keep = 0.5;
        auto_settings.roi_width_num_fwhm_lower = 3.5;  // Slightly more generous for refined fit
        auto_settings.roi_width_num_fwhm_upper = 3.5;
        auto_settings.max_fwhm_width = 12.0;  // Tighter constraint as solution improves
        auto_settings.min_fwhm_roi = 1.0;
        auto_settings.min_fwhm_quad_cont = 8;

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
              const SandiaDecay::Nuclide * const * nuc_ptr = std::get_if<const SandiaDecay::Nuclide *>( &nuc_act.source );
              if( nuc_ptr && *nuc_ptr )
              {
                const double live_time_seconds = foreground->live_time();
                // RelActAuto's rel_activity is per second, need to multiply by live time for clustering
                const double activity_for_clustering = nuc_act.rel_activity * live_time_seconds;

                if( PeakFitImprove::debug_printout )
                  cout << "  " << nuc_act.name() << ": rel_activity=" << nuc_act.rel_activity
                       << ", live_time=" << live_time_seconds
                       << "s, activity_for_clustering=" << activity_for_clustering << endl;

                nucs_and_acts.emplace_back( *nuc_ptr, activity_for_clustering );
              }
            }
          }

          // Cluster gammas using current solution's relative efficiency
          const vector<RelActCalcAuto::RoiRange> refined_rois = cluster_gammas_to_rois(
              auto_rel_eff, nucs_and_acts, foreground,
              fwhmFnctnlForm, fwhm_coefficients,
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


      // Score the fit results using only signal photopeaks
      FinalFit_GA::FinalFitScore fit_score;

      // Use m_peaks_without_back_sub for scoring
      const vector<PeakDef> &fit_peaks = solution.m_peaks_without_back_sub;

      // The number of sigma away from expected, that we will use to calculate the peak area contribution
      const double num_sigma_contribution = 1.5;

      for( const PeakDef &found_peak : fit_peaks )
      {
        const double found_energy = found_peak.mean();
        const double found_fwhm = found_peak.fwhm();
        const double peak_lower_contrib = found_energy - num_sigma_contribution * found_peak.sigma();
        const double peak_upper_contrib = found_energy + num_sigma_contribution * found_peak.sigma();

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

          const double area_uncert = found_peak.amplitudeUncert() > 0.0 ? found_peak.amplitudeUncert() : 1.0 / sqrt( found_peak.amplitude() );

          // If this is a significant peak that we didn't expect, penalize
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


      // Write N42 file with foreground, background, and fit peaks
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
