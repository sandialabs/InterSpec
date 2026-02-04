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

#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"

#include "InterSpec/FitPeaksForNuclides.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"
#include "SpecUtils/SpecUtilsAsync.h"

// Note: PeakFitImprove::debug_printout is used in the original code but is in target/peak_fit_improve/
// For now, we'll conditionally compile debug output based on PERFORM_DEVELOPER_CHECKS
// This can be refined later if needed
#ifndef PERFORM_DEVELOPER_CHECKS
#define PERFORM_DEVELOPER_CHECKS 0
#endif

namespace
{
  // Local debug flag - can be set externally if needed
  // Note: In the original code, this was PeakFitImprove::debug_printout from target/peak_fit_improve/
  // For the moved code, we use a local flag that can be controlled
  bool local_debug_printout = false;
  
  // Helper to check if debug output should be printed
  // Replaces PeakFitImprove::debug_printout from the original code
  inline bool should_debug_print()
  {
#if PERFORM_DEVELOPER_CHECKS
    return local_debug_printout;
#else
    return false;
#endif
  }
  
  // Macro to replace PeakFitImprove::debug_printout throughout the code
  // This allows easy replacement of all instances
#define PEAK_FIT_DEBUG_PRINTOUT should_debug_print()
  
  // Helper functions that need access to should_use_desperation_shielding and create_desperation_shielding
  // These are used by fit_peaks_for_nuclide_relactauto but were in the anonymous namespace in the original
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
}

using namespace std;

namespace FitPeaksForNuclides
{

// Anonymous namespace for helper functions
namespace
{
  // Forward declarations for structs used in helper functions
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
  };

  struct LocalMinimum
  {
    size_t channel;                  // Absolute channel number of minimum
    double synthetic_value;
    double depth_score;              // For tiebreaking (higher = better)
    double statistical_significance; // Primary criterion (lower = better breakpoint)
  };

  struct ClusteredGammaInfo
  {
    double lower;
    double upper;
    std::vector<double> gamma_energies;    // energies of gamma lines in this cluster
    std::vector<double> gamma_amplitudes;  // expected peak areas/amplitudes
  };
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

    throw std::runtime_error( "get_source_photons: invalid or null source" );
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
}//namespace


std::pair<double,double> find_valid_energy_range( const std::shared_ptr<const SpecUtils::Measurement> &meas )
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

  const std::vector<float> &channel_counts = *meas->gamma_counts();
  const size_t nbin = channel_counts.size();
  
  //First detect where spectrum begins
  const int side_bins = 3;
  const int order = 2;
  const int derivative = 2;
  std::vector<float> smoothed_2nd, smoothed_2nd_variance;
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


bool rois_are_similar( const std::vector<RelActCalcAuto::RoiRange> &a,
                       const std::vector<RelActCalcAuto::RoiRange> &b,
                       const double energy_tolerance = 1.0 )
{
  if( a.size() != b.size() )
    return false;

  for( size_t i = 0; i < a.size(); ++i )
  {
    if( std::fabs(a[i].lower_energy - b[i].lower_energy) > energy_tolerance )
      return false;
    if( std::fabs(a[i].upper_energy - b[i].upper_energy) > energy_tolerance )
      return false;
  }
  return true;
}//rois_are_similar


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
    }
    else
    {
      insignificant_roi_indices.push_back( roi_idx );
    }
  }//for( loop over ROIs )

  if( total_channels == 0 )
    return std::numeric_limits<double>::max();

  return total_chi2 / static_cast<double>( total_channels );
}//compute_filtered_chi2_dof


bool should_combine_peaks( const PeakDef &larger_peak,
                           const PeakDef &smaller_peak,
                           const double always_combine_nsigma )
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
        roi_peaks.push_back( combine_peaks( cluster ) );
      }
    }//for( const std::vector<const PeakDef *> &cluster : clusters )

    // Make all peaks in this ROI share a new continuum
    if( !roi_peaks.empty() )
    {
      std::shared_ptr<PeakContinuum> roi_continuum = std::make_shared<PeakContinuum>( *roi_peaks.front().continuum() );
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

  return combined_peaks;
}//combine_overlapping_peaks_in_rois


std::vector<PeakDef> compute_observable_peaks(
  const std::vector<PeakDef> &fit_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const PeakFitForNuclideConfig &config )
{
  
#if( PERFORM_DEVELOPER_CHECKS )
  // Debug: Check if input fit_peaks have overlapping ROIs
  if( should_debug_print() )
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

    std::sort( std::begin(rois_debug_input), std::end(rois_debug_input),
      []( const RoiDebugInfoInput &a, const RoiDebugInfoInput &b ) { return a.lower_energy < b.lower_energy; } );

    bool found_overlap_input = false;
    for( size_t i = 1; i < rois_debug_input.size(); ++i )
    {
      const RoiDebugInfoInput &prev_roi = rois_debug_input[i - 1];
      const RoiDebugInfoInput &curr_roi = rois_debug_input[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap_input )
          std::cerr << "compute_observable_peaks: INPUT fit_peaks ALREADY HAVE OVERLAPPING ROIs:" << std::endl;
        found_overlap_input = true;

        std::cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : prev_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;

        std::cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : curr_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;
        std::cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << std::endl;
      }
    }

    if( !found_overlap_input && should_debug_print() )
      std::cout << "compute_observable_peaks: Input fit_peaks have no overlapping ROIs" << std::endl;
  }//if( should_debug_print() )
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
    std::sort( std::begin(peaks), std::end(peaks), &PeakDef::lessThanByMean );

    const double new_left_mean = peaks.front().mean();
    const double new_right_mean = peaks.back().mean();
    const bool left_edge_changed = (std::fabs(new_left_mean - orig_left_mean) > 0.1);
    const bool right_edge_changed = (std::fabs(new_right_mean - orig_right_mean) > 0.1);

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

    if( should_debug_print() )
    {
      std::cout << "  Observable filter: adjusted ROI bounds from ["
           << old_continuum->lowerEnergy() << ", " << old_continuum->upperEnergy()
           << "] to [" << new_lower_energy << ", " << new_upper_energy << "] keV" << std::endl;
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
    std::sort( std::begin(roi_peaks), std::end(roi_peaks), &PeakDef::lessThanByMean );

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
      else if( should_debug_print() )
      {
        std::cout << "  Observable filter (initial sig=" << significance << " < " << initial_significance_threshold
             << "): peak at " << mean << " keV" << std::endl;
      }
    }//for( const PeakDef &peak : roi_peaks )

    // Adjust ROI bounds if edge peaks were removed
    if( !kept_peaks.empty() )
    {
      reduce_roi_bounds_if_needed( kept_peaks, orig_left_mean, orig_right_mean );
      filtered_peaks.insert( std::end(filtered_peaks), std::begin(kept_peaks), std::end(kept_peaks) );
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
  if( should_debug_print() )
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

    std::sort( std::begin(rois_debug2), std::end(rois_debug2),
      []( const RoiDebugInfo2 &a, const RoiDebugInfo2 &b ) { return a.lower_energy < b.lower_energy; } );

    bool found_overlap2 = false;
    for( size_t i = 1; i < rois_debug2.size(); ++i )
    {
      const RoiDebugInfo2 &prev_roi = rois_debug2[i - 1];
      const RoiDebugInfo2 &curr_roi = rois_debug2[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap2 )
          std::cerr << "compute_observable_peaks: OVERLAPS BEFORE PARALLEL PROCESSING:" << std::endl;
        found_overlap2 = true;

        std::cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy
             << "] keV, " << prev_roi.num_peaks << " peaks" << std::endl;
        std::cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy
             << "] keV, " << curr_roi.num_peaks << " peaks" << std::endl;
        std::cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << std::endl;
      }
    }

    if( !found_overlap2 && should_debug_print() )
      std::cout << "compute_observable_peaks: No overlaps before parallel processing" << std::endl;
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
      std::sort( std::begin(roi_peaks), std::end(roi_peaks), &PeakDef::lessThanByMean );

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
            if( should_debug_print() )
            {
              std::cout << "  Observable filter post-refit (final sig=" << final_sig << " < " << final_significance_threshold
                   << "): peak at " << mean << " keV" << std::endl;
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
    observable_peaks.insert( std::end(observable_peaks), std::begin(roi_result), std::end(roi_result) );

  // Sort by energy
  std::sort( std::begin(observable_peaks), std::end(observable_peaks), &PeakDef::lessThanByMean );

#if( PERFORM_DEVELOPER_CHECKS )
  // Debug: Check for overlapping ROI bounds
  if( should_debug_print() )
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

    std::sort( std::begin(rois_debug), std::end(rois_debug),
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
          std::cerr << "compute_observable_peaks: FOUND OVERLAPPING ROIs:" << std::endl;
        found_overlap = true;

        std::cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : prev_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;

        std::cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : curr_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;
        std::cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << std::endl;
      }
    }

    if( !found_overlap && should_debug_print() )
      std::cout << "compute_observable_peaks: No overlapping ROIs detected" << std::endl;
  }//if( should_debug_print() )
#endif

  return observable_peaks;
}//compute_observable_peaks


std::vector<double> build_synthetic_spectrum(
  const std::vector<double> &gamma_energies,
  const std::vector<double> &gamma_amplitudes,
  const std::function<double(double)> &fwhm_at_energy,
  const float *channel_energies,
  size_t start_channel,
  size_t num_channels )
{
  assert( gamma_energies.size() == gamma_amplitudes.size() );

  // Zero-initialize the synthetic spectrum (same binning as data)
  std::vector<double> synthetic( num_channels, 0.0 );

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


double compute_significance_in_region(
  const std::vector<double> &synthetic,
  size_t start_channel,
  size_t check_start_ch,
  size_t check_end_ch,
  const std::shared_ptr<const SpecUtils::Measurement> &data )
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


bool has_significant_peak_between(
  size_t lower_channel,
  size_t upper_channel,
  const std::vector<double> &synthetic,
  size_t start_channel,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const float *channel_energies,
  const std::function<double(double)> &fwhm_at_energy,
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
  for( size_t ch = lower_channel; ch < upper_channel; ++ch )
  {
    const size_t i = ch - start_channel;
    if( i == 0 || i >= num_channels - 1 )
      continue;

    // Check if this is a local maximum in the synthetic spectrum
    if( synthetic[i] > synthetic[i-1] && synthetic[i] > synthetic[i+1] )
    {
      // Get energy at this channel for FWHM calculation
      const double ch_energy = 0.5 * (channel_energies[ch] + channel_energies[ch + 1]);
      const double fwhm = fwhm_at_energy( ch_energy );

      // Convert FWHM fraction to channel range for significance check
      const double half_width = check_fwhm_fraction * fwhm;
      const size_t check_start = data->find_gamma_channel( static_cast<float>(ch_energy - half_width) );
      const size_t check_end = data->find_gamma_channel( static_cast<float>(ch_energy + half_width) );

      const double significance = compute_significance_in_region(
        synthetic, start_channel, check_start, check_end, data );

      if( significance >= significance_threshold )
        return true;  // Found a significant peak
    }
  }

  return false;  // No significant peak found
}//has_significant_peak_between


std::vector<LocalMinimum> find_synthetic_minima(
  const std::vector<double> &synthetic,
  size_t start_channel,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const float *channel_energies,
  const std::function<double(double)> &fwhm_at_energy,
  double check_fwhm_fraction )
{
  std::vector<LocalMinimum> minima;
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
      const double smaller_max = std::min( left_max, right_max );
      const double larger_max = std::max( left_max, right_max );
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


bool should_use_step_continuum(
  const ClusteredGammaInfo &cluster,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
  const std::vector<float> &fwhm_coefficients,
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
  const auto max_it = std::max_element( std::begin(cluster.gamma_amplitudes), std::end(cluster.gamma_amplitudes) );
  const size_t max_idx = static_cast<size_t>( std::distance( std::begin(cluster.gamma_amplitudes), max_it ) );
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

  // If left side is significantly higher than right side, suggest step continuum
  return (nsigma >= step_cont_left_right_nsigma);
}//should_use_step_continuum


// Due to the large size of the remaining helper functions, they will be added
// by copying from FitPeaksForNuclideDev.cpp. For now, adding stubs that will
// be replaced with full implementations.

std::vector<RelActCalcAuto::RoiRange> cluster_gammas_to_rois(
    const std::vector<std::function<double(double)>> &rel_eff_fcns,
    const std::vector<std::vector<std::tuple<RelActCalcAuto::SrcVariant, double /*age*/, double/*act*/>>> &sources_age_activity_sets,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
    const std::vector<float> &fwhm_coefficients,
    const double fwhm_lower_energy,
    const double fwhm_upper_energy,
    const double lowest_energy,
    const double highest_energy,
    const GammaClusteringSettings &settings )
{
  assert( rel_eff_fcns.size() == sources_age_activity_sets.size() );
  if( rel_eff_fcns.size() != sources_age_activity_sets.size() )
    throw runtime_error( "cluster_gammas_to_rois: there is a different number of relative efficiency functions and sets of sources" );

  vector<RelActCalcAuto::RoiRange> result_rois;

  // Collect all gamma lines with their expected counts
  vector<std::pair<double,double>> gammas_by_counts;  // (energy, expected_counts)

  for( size_t rel_eff_index = 0; rel_eff_index < rel_eff_fcns.size(); ++rel_eff_index )
  {
    const function<double(double)> &rel_eff_fcn = rel_eff_fcns[rel_eff_index];
    const vector<tuple<RelActCalcAuto::SrcVariant,double,double>> &src_age_and_activities = sources_age_activity_sets[rel_eff_index];

    for( const tuple<RelActCalcAuto::SrcVariant,double,double> &src_age_act : src_age_and_activities )
    {
      const RelActCalcAuto::SrcVariant &src = get<0>(src_age_act);
      const double age = get<1>(src_age_act);
      const double activity = get<2>(src_age_act);

      if( RelActCalcAuto::is_null(src) || (activity <= 0.0) )
        continue;

      if( should_debug_print() )
      {
        std::cerr << "cluster_gammas_to_rois: Source " << RelActCalcAuto::to_name( src ) << ", activity=" << activity
        << ", age=" << (age / PhysicalUnits::second) << " seconds ("
        << (age / PhysicalUnits::year) << " years)" << std::endl;
      }

      const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src, activity, age );

      if( should_debug_print() )
      {
        std::cerr << "  " << photons.size() << " photons from " << RelActCalcAuto::to_name( src ) << ", energy range ["
        << lowest_energy << ", " << highest_energy << "] keV" << std::endl;
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
  }//for( size_t rel_eff_index = 0; rel_eff_index < rel_eff_fcns.size(); ++rel_eff_index )

  if( gammas_by_counts.empty() )
    return result_rois;

  // Create a copy sorted by energy for efficient lookup
  std::vector<std::pair<double,double>> gammas_by_energy = gammas_by_counts;
  const auto lessThanByEnergy = []( const std::pair<double,double> &lhs, const std::pair<double,double> &rhs ) {
    return lhs.first < rhs.first;
  };
  std::sort( std::begin(gammas_by_energy), std::end(gammas_by_energy), lessThanByEnergy );

  // Sort gammas by expected counts (highest first)
  std::sort( std::begin(gammas_by_counts), std::end(gammas_by_counts),
    []( const std::pair<double,double> &lhs, const std::pair<double,double> &rhs ) {
      return lhs.second > rhs.second;
  } );

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: Input gammas (" << gammas_by_counts.size() << " total):" << std::endl;
    std::cerr << "  Top 20 by expected counts:" << std::endl;
    for( size_t i = 0; i < std::min( gammas_by_counts.size(), size_t(20) ); ++i )
      std::cerr << "    " << gammas_by_counts[i].first << " keV, est_counts=" << gammas_by_counts[i].second << std::endl;
  }

  // Cluster gamma lines
  std::vector<ClusteredGammaInfo> clustered_gammas;

  // Check if we have a valid FWHM energy range for clamping
  const bool have_fwhm_range = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0);

  for( const std::pair<double,double> &energy_counts : gammas_by_counts )
  {
    auto ene_pos = std::lower_bound( std::begin(gammas_by_energy), std::end(gammas_by_energy),
                                    energy_counts, lessThanByEnergy );
    if( ene_pos == std::end(gammas_by_energy) )
      continue;

    if( ene_pos->first != energy_counts.first )
    {
      // Already removed from gammas_by_energy (absorbed into another cluster)
      if( should_debug_print() && (energy_counts.first >= 800.0) && (energy_counts.first <= 820.0) )
      {
        std::cerr << "cluster_gammas_to_rois: Gamma at " << energy_counts.first
             << " keV already absorbed into another cluster" << std::endl;
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
      if( should_debug_print() )
        std::cerr << "Warning: Invalid FWHM=" << fwhm << " at energy=" << energy << " keV, skipping gamma line" << std::endl;
      continue;
    }

    const double sigma = fwhm / PhysicalUnits::fwhm_nsigma;

    const double lower = std::max( lowest_energy, energy - settings.cluster_num_sigma * sigma );
    const double upper = std::min( highest_energy, energy + settings.cluster_num_sigma * sigma );

    // Find all gammas in this range
    const auto start_remove = std::lower_bound( std::begin(gammas_by_energy), std::end(gammas_by_energy),
                                               std::make_pair(lower, 0.0), lessThanByEnergy );
    const auto end_remove = std::upper_bound( std::begin(gammas_by_energy), std::end(gammas_by_energy),
                                             std::make_pair(upper, 0.0), lessThanByEnergy );

    const double counts_in_region = std::accumulate( start_remove, end_remove, 0.0,
        []( const double &sum, const std::pair<double,double> &el ) {
          return sum + el.second;
    } );

    // Capture the gamma lines before erasing them - as separate energy and amplitude arrays
    std::vector<double> gamma_energies_in_cluster;
    std::vector<double> gamma_amplitudes_in_cluster;
    for( auto it = start_remove; it != end_remove; ++it )
    {
      gamma_energies_in_cluster.push_back( it->first );
      gamma_amplitudes_in_cluster.push_back( it->second );
    }

    const double data_area = foreground->gamma_integral( static_cast<float>(lower), static_cast<float>(upper) );

    gammas_by_energy.erase( start_remove, end_remove );

    const double signif = (data_area > 0.0) ? (counts_in_region / std::sqrt(data_area)) : (1.0 / std::sqrt(counts_in_region));

    // Additional safety check - ensure lower and upper are finite and valid
    if( !std::isfinite(lower) || !std::isfinite(upper) || (lower >= upper) )
    {
      if( should_debug_print() )
        std::cerr << "Warning: Invalid cluster bounds [" << lower << ", " << upper << "], skipping" << std::endl;
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
    
    if( should_debug_print() )
    {
      const std::string status_str = (passes_data_area && passes_counts && passes_signif) ? "Accepted" : "Rejected";
      std::cerr << "cluster_gammas_to_rois: " << status_str << " [" << std::fixed << std::setprecision(1) << lower << ", " << upper << "] keV (e="
           << energy << " keV): "
           << "data=" << data_area << (passes_data_area ? " > " : " < ") << settings.min_data_area_keep << "; "
           << "est_counts=" << counts_in_region << (passes_counts ? " > " : " < ") << settings.min_est_peak_area_keep << "; "
         << "sig=" << signif << (passes_signif ? " > " : " < ") << settings.min_est_significance_keep << "; "
         << std::endl;
    }
  }//for( const std::pair<double,double> &energy_counts : gammas_by_counts )

  // Sort by lower energy
  std::sort( std::begin(clustered_gammas), std::end(clustered_gammas),
    []( const ClusteredGammaInfo &lhs, const ClusteredGammaInfo &rhs ) {
      return lhs.lower < rhs.lower;
  } );

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: Initial " << clustered_gammas.size() << " clustered ROIs:" << std::endl;
    for( size_t i = 0; i < clustered_gammas.size(); ++i )
    {
      const ClusteredGammaInfo &c = clustered_gammas[i];
      std::cerr << "  [" << i << "] range=[" << c.lower << ", " << c.upper << "] keV, "
           << c.gamma_energies.size() << " gammas: ";
      for( size_t j = 0; j < std::min( c.gamma_energies.size(), size_t(5) ); ++j )
        std::cerr << c.gamma_energies[j] << (j + 1 < c.gamma_energies.size() ? ", " : "");
      if( c.gamma_energies.size() > 5 )
        std::cerr << "... (" << c.gamma_energies.size() - 5 << " more)";
      std::cerr << std::endl;
    }
  }

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
    
    if( should_debug_print() )
    {
      std::cerr << "cluster_gammas_to_rois: Setting effective FWHM bounds for cluster with "
      << cluster.gamma_energies.size() << " gammas:" << std::endl
      << "  weighted_mean=" << weighted_mean << " keV, effective_fwhm=" << effective_fwhm << " keV" << std::endl
      << "  old range=[" << old_lower << ", " << old_upper << "] keV"
      << " -> new range=[" << cluster.lower << ", " << cluster.upper << "] keV" << std::endl;
    }
  }//for( ClusteredGammaInfo &cluster : clustered_gammas )
  
  // Merge overlapping clusters
  std::vector<ClusteredGammaInfo> merged_clusters;
  for( const ClusteredGammaInfo &cluster : clustered_gammas )
  {
    if( merged_clusters.empty() || (cluster.lower > merged_clusters.back().upper) )
    {
      merged_clusters.push_back( cluster );
    }
    else
    {
      if( should_debug_print() )
      {
        std::cerr << "cluster_gammas_to_rois: Merging cluster [" << cluster.lower << ", " << cluster.upper
             << "] into [" << merged_clusters.back().lower << ", " << merged_clusters.back().upper << "]"
             << " -> new upper=" << std::max( merged_clusters.back().upper, cluster.upper ) << std::endl;
      }

      merged_clusters.back().upper = std::min( highest_energy, std::max( merged_clusters.back().upper, cluster.upper ) );
      merged_clusters.back().gamma_energies.insert(
        std::end(merged_clusters.back().gamma_energies),
        std::begin(cluster.gamma_energies),
        std::end(cluster.gamma_energies)
      );
      merged_clusters.back().gamma_amplitudes.insert(
        std::end(merged_clusters.back().gamma_amplitudes),
        std::begin(cluster.gamma_amplitudes),
        std::end(cluster.gamma_amplitudes)
      );
    }
  }//for( const ClusteredGammaInfo &cluster : clustered_gammas )

  // Validate merged clusters - ensure merge didn't introduce NaN values
  merged_clusters.erase(
    std::remove_if( std::begin(merged_clusters), std::end(merged_clusters),
      []( const ClusteredGammaInfo &cluster ) {
        const bool invalid = !std::isfinite(cluster.lower) || !std::isfinite(cluster.upper) || (cluster.lower >= cluster.upper);
        if( invalid && should_debug_print() )
          std::cerr << "Warning: Removing invalid merged cluster with bounds [" << cluster.lower << ", " << cluster.upper << "]" << std::endl;
        return invalid;
      } ),
    std::end(merged_clusters)
  );

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: After merging, " << merged_clusters.size() << " clusters:" << std::endl;
    for( size_t i = 0; i < merged_clusters.size(); ++i )
    {
      const ClusteredGammaInfo &c = merged_clusters[i];
      std::cerr << "  [" << i << "] range=[" << c.lower << ", " << c.upper << "] keV ("
           << (c.upper - c.lower) << " keV wide), " << c.gamma_energies.size() << " gammas" << std::endl;
    }
  }

  // Break up ROIs that are too wide
  std::vector<ClusteredGammaInfo> final_clusters;
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

    if( should_debug_print() )
    {
      std::cerr << "cluster_gammas_to_rois: Cluster [" << cluster.lower << ", " << cluster.upper
           << "] keV: width=" << current_width << " keV, mid_fwhm=" << mid_fwhm
           << " keV, max_fwhm_width=" << settings.max_fwhm_width
           << ", max_width=" << max_width << " keV, "
           << (current_width <= max_width ? "NOT breaking" : "NEEDS breaking") << std::endl;
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
      std::vector<size_t> indices( cluster.gamma_energies.size() );
      std::iota( std::begin(indices), std::end(indices), 0 );
      std::sort( std::begin(indices), std::end(indices),
        [&cluster]( size_t a, size_t b ) {
          return cluster.gamma_energies[a] < cluster.gamma_energies[b];
        } );

      // Reorder both arrays based on sorted indices
      std::vector<double> sorted_energies( cluster.gamma_energies.size() );
      std::vector<double> sorted_amplitudes( cluster.gamma_amplitudes.size() );
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
    std::vector<double> synthetic = build_synthetic_spectrum(
      cluster.gamma_energies,
      cluster.gamma_amplitudes,
      fwhm_at_energy,
      channel_energies,
      start_channel,
      num_channels );

    // Find local minima in synthetic spectrum (with significance computed)
    const double start_energy = channel_energies[start_channel];
    const double end_energy = channel_energies[end_channel];
    const bool debug_this_cluster = should_debug_print() &&
      ((start_energy >= 60.0 && end_energy <= 140.0) ||
       (start_energy >= 140.0 && start_energy <= 160.0));

    if( debug_this_cluster )
    {
      std::cerr << "DEBUG: Breaking cluster [" << cluster.lower << ", " << cluster.upper
           << "] keV: start_channel=" << start_channel
           << ", end_channel=" << end_channel
           << ", num_channels=" << num_channels << std::endl;
    }

    std::vector<LocalMinimum> minima = find_synthetic_minima(
      synthetic,
      start_channel,
      foreground,
      channel_energies,
      fwhm_at_energy,
      settings.break_check_fwhm_fraction );

    if( debug_this_cluster )
    {
      std::cerr << "DEBUG: Found " << minima.size() << " minima in cluster ["
           << cluster.lower << ", " << cluster.upper << "] keV:" << std::endl;
      for( size_t i = 0; i < minima.size(); ++i )
      {
        const LocalMinimum &m = minima[i];
        std::cerr << "  [" << i << "] channel=" << m.channel
             << ", energy=" << channel_energies[m.channel] << " keV"
             << ", synthetic_value=" << m.synthetic_value
             << ", depth_score=" << m.depth_score
             << ", stat_sig=" << m.statistical_significance << std::endl;
      }
    }

    // If no minima found, don't break the ROI
    if( minima.empty() )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }
    
    // Sort by statistical_significance (ascending - least significant first)
    // Use depth_score as tiebreaker (descending - deepest first)
    std::sort( std::begin(minima), std::end(minima),
      []( const LocalMinimum &a, const LocalMinimum &b ) {
        if( a.statistical_significance != b.statistical_significance )
          return a.statistical_significance < b.statistical_significance;
        return a.depth_score > b.depth_score;
      });

    // Select breakpoints, validating significant peaks exist between them
    std::vector<size_t> selected_breakpoint_channels;
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
        std::cerr << "DEBUG: Evaluating breakpoint candidate at channel " << candidate.channel
             << " (" << channel_energies[candidate.channel] << " keV): "
             << "has_left_peak=" << (has_left_peak ? "YES" : "NO")
             << ", has_right_peak=" << (has_right_peak ? "YES" : "NO")
             << ", accepted=" << ((has_left_peak && has_right_peak) ? "YES" : "NO") << std::endl;
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
      std::sort( std::begin(selected_breakpoint_channels), std::end(selected_breakpoint_channels) );

      // Convert breakpoint channels to energies
      std::vector<double> breakpoint_energies;
      for( size_t bp_ch : selected_breakpoint_channels )
        breakpoint_energies.push_back( channel_energies[bp_ch] );

      if( should_debug_print() )
      {
        std::cerr << "cluster_gammas_to_rois: Breaking up cluster [" << cluster.lower << ", " << cluster.upper
             << "] at " << breakpoint_energies.size() << " breakpoints: ";
        for( size_t i = 0; i < breakpoint_energies.size(); ++i )
          std::cerr << breakpoint_energies[i] << (i + 1 < breakpoint_energies.size() ? ", " : "");
        std::cerr << " keV" << std::endl;
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
      std::cerr << "ERROR: final_clusters[" << (i-1) << "] and [" << i << "] overlap: "
           << "[" << prev.lower << ", " << prev.upper << "] vs "
           << "[" << curr.lower << ", " << curr.upper << "]" << std::endl;
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
      if( should_debug_print() )
      {
        std::cerr << "cluster_gammas_to_rois: Rejected ROI [" << roi.lower_energy << ", " << roi.upper_energy
             << "] keV: too narrow (" << num_fwhm_wide << " FWHM < " << settings.min_fwhm_roi << " min)" << std::endl;
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
      const double max_amplitude = *std::max_element( std::begin(cluster.gamma_amplitudes),
                                                  std::end(cluster.gamma_amplitudes) );
      const double data_area = foreground->gamma_integral( static_cast<float>(roi.lower_energy),
                                                            static_cast<float>(roi.upper_energy) );
      const double est_significance = (data_area > 0.0)
          ? (max_amplitude / std::sqrt(data_area))
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

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: Final " << result_rois.size() << " ROIs:" << std::endl;
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
      std::cerr << "  [" << i << "] range=[" << roi.lower_energy << ", " << roi.upper_energy << "] keV ("
           << width << " keV, " << (width / fwhm) << " FWHM), cont=" << cont_str << std::endl;
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
      std::cerr << "ERROR: result_rois[" << (i-1) << "] and [" << i << "] overlap: "
           << "[" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] vs "
           << "[" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "]" << std::endl;
      assert( curr_roi.lower_energy >= prev_roi.upper_energy );
    }
  }
#endif

  return result_rois;
}//cluster_gammas_to_rois


struct InitialRoi
{
  RelActCalcAuto::RoiRange roi;
  double center_energy;
  double fwhm;
};

std::vector<RelActCalcAuto::RoiRange> merge_rois( std::vector<InitialRoi> initial_rois, const PeakFitForNuclideConfig &config )
{
  // Sort by lower_energy for merging
  std::sort( initial_rois.begin(), initial_rois.end(), [](const InitialRoi &a, const InitialRoi &b){
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
    const double last_mid = 0.5 * (last.lower_energy + last.upper_energy);
    const double current_mid = 0.5 * (current.roi.lower_energy + current.roi.upper_energy);
    const double mid_dist = (fabs(current_mid - last_mid) < 0.1) ? 0.5 : ((mid_energy - last_mid) / (current_mid - last_mid));
    const double mid_fwhm = last_fwhm + mid_dist*(current.fwhm - last_fwhm);
    assert( ((mid_fwhm >= current.fwhm) && (mid_fwhm <= last_fwhm)) || ((mid_fwhm <= current.fwhm) && (mid_fwhm >= last_fwhm)) );
    //const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
    //  static_cast<float>(mid_energy), fwhmFnctnlForm, fwhm_coefficients );

    const bool width_ok = (combined_width <= config.auto_rel_eff_sol_max_fwhm * mid_fwhm);

    if( width_ok )
    {
      // MERGE: Extend last ROI to encompass both
      last.upper_energy = combined_upper;
      last_centers.push_back( current.center_energy );
      // Update FWHM to use the average
      merged_fwhms.back() = 0.5 * (last_fwhm + current.fwhm);

      if( should_debug_print() )
      {
        std::cerr << "Merged overlapping ROI [" << current.roi.lower_energy << ", "
             << current.roi.upper_energy << "] into [" << last.lower_energy
             << ", " << last.upper_energy << "]" << std::endl;
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

        if( should_debug_print() )
        {
          std::cerr << "Removed last ROI (invalid after split), kept current ROI at "
               << current.center_energy << " keV" << std::endl;
        }
      }
      else if( current_valid )
      {
        // Both ROIs valid - add the split current ROI
        merged_rois.push_back( adjusted_current );
        merged_centers.push_back( {current.center_energy} );
        merged_fwhms.push_back( current.fwhm );

        if( should_debug_print() )
        {
          std::cerr << "Split overlapping ROIs (width constraint): overlap=" << overlap
               << " keV, last=[" << last.lower_energy << ", " << last.upper_energy
               << "], current=[" << adjusted_current.lower_energy << ", "
               << adjusted_current.upper_energy << "]" << std::endl;
        }
      }
      else
      {
        // Current ROI invalid after split, but last is valid - just keep last as adjusted
        if( should_debug_print() )
        {
          std::cerr << "Skipping adjusted current ROI at " << current.center_energy << " keV: ";
          if( !current_contains_center )
            std::cerr << "doesn't contain source energy";
          else
            std::cerr << "too narrow (" << adjusted_width << " keV < " << current.fwhm << " keV FWHM)";
          std::cerr << std::endl;
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

  if( should_debug_print() )
  {
    std::cerr << "estimate_initial_rois_without_peaks: Created " << merged_rois.size()
         << " final ROIs" << std::endl;
  }

  return merged_rois;
}//std::vector<InitialRoi> merge_rois(...)


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

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
      continue;

    // Get source age and photons
    const double age = src.age;
    const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src.source, 1.0, age );

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
        candidates.push_back( {photon.energy, score, src.source} );
    }

    // Sort by score (descending) and take top 4
    std::sort( candidates.begin(), candidates.end(),
      [](const GammaInfo &a, const GammaInfo &b) { return a.br_times_eff > b.br_times_eff; } );

    const size_t num_to_take = std::min( candidates.size(), size_t(4) );
    for( size_t i = 0; i < num_to_take; ++i )
      selected_gammas.push_back( candidates[i] );

    // Debug output
    if( should_debug_print() )
    {
      std::cerr << "estimate_initial_rois_without_peaks: Source "
           << RelActCalcAuto::to_name( src.source ) << " - selected " << num_to_take << " gammas" << std::endl;
      for( size_t i = 0; i < num_to_take; ++i )
      {
        std::cerr << "  " << candidates[i].energy << " keV, BR*eff=" << candidates[i].br_times_eff << std::endl;
      }
    }
  }

  if( selected_gammas.empty() )
  {
    if( should_debug_print() )
      std::cerr << "estimate_initial_rois_without_peaks: No valid gammas found" << std::endl;
    return {};
  }

  // Step 3: Create initial ROIs
  std::vector<InitialRoi> initial_rois;

  for( const GammaInfo &gamma : selected_gammas )
  {
    // Compute FWHM
    const double fwhm = DetectorPeakResponse::peakResolutionFWHM(
      static_cast<float>(gamma.energy), fwhmFnctnlForm, fwhm_coefficients );

    // Validate FWHM
    if( !std::isfinite( fwhm ) || fwhm <= 0.0 )
    {
      if( should_debug_print() )
        std::cerr << "Warning: Invalid FWHM at " << gamma.energy << " keV, skipping" << std::endl;
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
    if( should_debug_print() )
      std::cerr << "estimate_initial_rois_without_peaks: No valid ROIs created" << std::endl;
    return {};
  }

  return merge_rois( initial_rois, config );
}//estimate_initial_rois_without_peaks

std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_fallback(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
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
  std::vector<tuple<RelActCalcAuto::SrcVariant, double, double>> source_age_and_acts;

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
      continue;

    // Get gamma lines at default age
    const double age = src.age;
    const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src.source, 1.0, age );

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

      std::cout << "Fallback: " << RelActCalcAuto::to_name( src.source ) << " matched peak at " << matched_peak->mean()
           << " keV (gamma at " << best_energy << " keV), area=" << peak_area
           << ", estimated activity=" << estimated_activity << std::endl;
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

      std::cout << "Fallback: " << RelActCalcAuto::to_name( src.source ) << " no matching peak for gamma at " << best_energy
           << " keV, integrated counts=" << total_counts << ", est. peak area=" << estimated_peak_area
           << ", estimated activity=" << estimated_activity << std::endl;
    }

    if( estimated_activity > 0.0 )
      source_age_and_acts.emplace_back( src.source, age, estimated_activity );
  }

  if( source_age_and_acts.empty() )
  {
    std::cerr << "Fallback: Could not estimate activity for any source" << std::endl;
    return {};
  }

  // Step 3: Create relative efficiency function from DRF
  const auto fallback_rel_eff = [&drf_to_use]( double energy ) -> double {
    return drf_to_use->intrinsicEfficiency( static_cast<float>(energy) );
  };

  // Step 4: Call cluster_gammas_to_rois with estimated activities
  return cluster_gammas_to_rois(
    {fallback_rel_eff},
    {source_age_and_acts},
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

std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_using_relactmanual(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
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
  std::vector<RelActCalcAuto::RoiRange> initial_rois;

  // Step 1: Convert auto_search_peaks to RelActManual format
  std::vector<RelActCalcManual::GenericPeakInfo> rel_act_manual_peaks;

  for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
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
  std::vector<RelActCalcManual::PeakCsvInput::NucAndAge> rel_act_manual_srcs;
  for( const RelActCalcAuto::NucInputInfo &src : sources )
    rel_act_manual_srcs.emplace_back( RelActCalcAuto::to_name( src.source ), -1.0, false );

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

  std::vector<RelActCalcManual::GenericPeakInfo> peaks_matched = peak_match_results.peaks_matched;
  std::sort( std::begin(peaks_matched), std::end(peaks_matched),
    []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ){
      return lhs.m_energy < rhs.m_energy;
  });

  if( should_debug_print() )
  {
    const std::vector<std::string> &used_isotopes = peak_match_results.used_isotopes;
    const std::vector<std::string> &unused_isotopes = peak_match_results.unused_isotopes;
    if( unused_isotopes.empty() )
      std::cout << "Matched up all source nuclides to initial peak fit" << std::endl;
    std::cout << "Failed to match up nuclides: {";
    for( const std::string &nuc : unused_isotopes )
      std::cout << nuc << ", ";
    std::cout << "} to initial auto-fit peaks" << std::endl;
    if( !used_isotopes.empty() )
    {
      std::cout << "Matched up nuclides: {";
      for( const std::string &nuc : used_isotopes )
        std::cout << nuc << ", ";
      std::cout << "} to initial auto-fit peaks to a total of " << peaks_matched.size() << " peaks." << std::endl;

      if( !peaks_matched.empty() )
      {
        std::cout << "Matched peak energies (keV): {";
        for( size_t i = 0; i < peaks_matched.size(); ++i )
        {
          std::cout << peaks_matched[i].m_energy;
          if( i < peaks_matched.size() - 1 )
            std::cout << ", ";
        }
        std::cout << "}" << std::endl;
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

    manual_input.phys_model_self_atten = std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>{};
    manual_input.phys_model_external_attens = std::vector<std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>>{};
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
      if( should_debug_print() )
      {
        std::cout << "Initial manual solution failed: status=";
        switch( manual_solution.m_status )
        {
          case RelActCalcManual::ManualSolutionStatus::NotInitialized:
            std::cout << "NotInitialized";
            break;
          case RelActCalcManual::ManualSolutionStatus::ErrorInitializing:
            std::cout << "ErrorInitializing";
            break;
          case RelActCalcManual::ManualSolutionStatus::ErrorFindingSolution:
            std::cout << "ErrorFindingSolution";
            break;
          case RelActCalcManual::ManualSolutionStatus::ErrorGettingSolution:
            std::cout << "ErrorGettingSolution";
            break;
          case RelActCalcManual::ManualSolutionStatus::Success:
            std::cout << "Success";
            break;
        }
        std::cout << ", form=" << RelActCalc::to_str( manual_solution.m_input.eqn_form )
             << ", order=" << manual_solution.m_input.eqn_order
             << ", chi2=" << manual_solution.m_chi2
             << ", dof=" << manual_solution.m_dof
             << ", chi2/dof=" << chi2_dof;
        if( !manual_solution.m_error_message.empty() )
          std::cout << ", error: " << manual_solution.m_error_message;
        std::cout << std::endl;
      }//if( should_debug_print() )
      RelActCalcManual::RelEffInput retry_input = manual_input;
      retry_input.eqn_form = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      retry_input.phys_model_use_hoerl = (manual_input.peaks.size() > 3);
      retry_input.use_ceres_to_fit_eqn = true;
      retry_input.eqn_order = 0;
      if( isHPGe )
        retry_input.phys_model_detector = DetectorPeakResponse::getGenericHPGeDetector();
      else
        retry_input.phys_model_detector = DetectorPeakResponse::getGenericNaIDetector();

      RelActCalcManual::RelEffSolution retry_solution
        = RelActCalcManual::solve_relative_efficiency( retry_input );

      double retry_chi2_dof = retry_solution.m_chi2 / (std::max)(retry_solution.m_dof, 1);

      if( should_debug_print() )
      {
        std::cout << "Retry manual solution fit comparison:" << std::endl;
        std::cout << "  Original: form=" << RelActCalc::to_str( manual_solution.m_input.eqn_form )
        << ", order=" << manual_solution.m_input.eqn_order
        << ", chi2/dof=" << manual_solution.m_chi2 << "/" << manual_solution.m_dof
        << "=" << chi2_dof
        << ", success=" << (manual_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success)
        << std::endl;
        if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
          std::cout << "  Original error: " << manual_solution.m_error_message << std::endl;
        std::cout << "  Retry:    form=" << RelActCalc::to_str( retry_solution.m_input.eqn_form )
        << ", order=" << retry_solution.m_input.eqn_order
        << ", chi2/dof=" << retry_solution.m_chi2 << "/" << retry_solution.m_dof
        << "=" << retry_chi2_dof
        << ", success=" << (retry_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success)
        << std::endl;
        if( retry_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
          std::cout << "  Retry error: " << retry_solution.m_error_message << std::endl;
        std::cout << std::endl;
      }

      if( (retry_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success)
         && (retry_chi2_dof < chi2_dof) )
      {
        if( should_debug_print() )
          std::cout << "Will use retry solution!" << std::endl;
        chi2_dof = retry_chi2_dof;
        manual_solution = std::move(retry_solution);
      }else
      {
        if( should_debug_print() )
          std::cout << "Will use original solution!" << std::endl;
      }
    }

    if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
      throw std::runtime_error( "Failed to fit initial RelActCalcManual::RelEffSolution: " + manual_solution.m_error_message );

    if( should_debug_print() )
    {
      std::cout << "Successfully fitted initial RelActCalcManual::RelEffSolution: chi2/dof="
      << manual_solution.m_chi2 << "/" << manual_solution.m_dof << "="
      << manual_solution.m_chi2 / manual_solution.m_dof
      << " using " << peaks_matched.size() << " peaks"
      << std::endl;
      std::cout << std::endl;
    }

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
    vector<tuple<RelActCalcAuto::SrcVariant, double, double>> source_age_and_acts;
    for( const RelActCalcManual::IsotopeRelativeActivity &rel_act : manual_solution.m_rel_activities )
    {
      const RelActCalcAuto::SrcVariant src = RelActCalcAuto::source_from_string( rel_act.m_isotope );
      if( RelActCalcAuto::is_null( src ) )
        throw std::logic_error( "Failed to get source from RelAct isotope '" + rel_act.m_isotope + "'" );

      double age = 0.0;
      if( RelActCalcAuto::nuclide(src) )
      {
        bool found_src = false;
        for( const RelActCalcAuto::NucInputInfo &input_src : sources )
        {
          if( RelActCalcAuto::to_name(input_src.source) == rel_act.m_isotope )
          {
            age = input_src.age;
            found_src = true;
            break;
          }
        }//for( const RelActCalcAuto::NucInputInfo &input_src : sources )

        assert( found_src );
      }//if( RelActCalcAuto::nuclide(src) )

      const double act = manual_solution.relative_activity( rel_act.m_isotope );
      source_age_and_acts.emplace_back( src, age, act );
    }//for( loop over manual_solution.m_rel_activities )


    // Step 5: Use the reusable clustering function to create ROIs
    initial_rois = cluster_gammas_to_rois( {manual_rel_eff}, {source_age_and_acts}, foreground,
                                          fwhmFnctnlForm, fwhm_coefficients,
                                          lower_fwhm_energy, upper_fwhm_energy,
                                          min_valid_energy, max_valid_energy,
                                          manual_settings );

    if( should_debug_print() )
    {
      std::cout << "Initial ROIs from RelActManual: ";
      for( const auto &roi : initial_rois )
        std::cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
      std::cout << std::endl;
    }
  }catch( std::exception &e )
  {
    std::cerr << "Error trying to fit initial manual rel-eff solution: " << e.what() << std::endl;
    std::cerr << "Using fallback activity estimation..." << std::endl;

    initial_rois = estimate_initial_rois_fallback(
      auto_search_peaks, foreground, sources, drf, isHPGe,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      manual_settings );

    fallback_warning = "RelActManual fitting failed (" + std::string(e.what())
                     + "); used simplified activity estimation fallback";

    std::cout << "Fallback ROIs: ";
    for( const RelActCalcAuto::RoiRange &roi : initial_rois )
      std::cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
    std::cout << std::endl;
  }

  return initial_rois;
}//estimate_initial_rois_using_relactmanual

PeakFitResult fit_peaks_for_nuclide_relactauto(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &orig_foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::vector<RelActCalcAuto::RoiRange> &input_rois,
  const std::shared_ptr<const SpecUtils::Measurement> &orig_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const Wt::WFlags<FitSrcPeaksOptions> user_options,
  const PeakFitForNuclideConfig &config,
  const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
  const std::vector<float> &fwhm_coefficients,
  const bool isHPGe,
  const double fwhm_lower_energy,
  const double fwhm_upper_energy )
{
  PeakFitResult result;

  const bool fit_norm_peaks = user_options.testFlag(FitSrcPeaksOptions::FitNormBkgrndPeaks);
  const bool apply_energy_cal_between = config.fit_energy_cal
                                        && !user_options.testFlag(FitSrcPeaksOptions::DoNotVaryEnergyCal)
                                        && !user_options.testFlag(FitSrcPeaksOptions::DoNotRefineEnergyCal);
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );

  // Validate sources
  if( sources.empty() && !fit_norm_peaks )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "No sources provided";
    return result;
  }

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Null source in sources vector";
      return result;
    }
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
  rel_eff_curve.nuclides = sources;

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
      rel_eff_curve.name += " " + (el ? el->symbol : ("z=" + std::to_string(z)) ) + std::string(" Peak Fit");
    }
  }else
  {
    rel_eff_curve.name = RelActCalc::to_str(config.rel_eff_eqn_type) + std::string(" Peak Fit");
  }

  // Create RelActAuto options from config
  RelActCalcAuto::Options options;
  options.rel_eff_curves.push_back( rel_eff_curve );
  options.rois = input_rois;
  options.fit_energy_cal = config.fit_energy_cal;
  options.fwhm_form = config.fwhm_form;
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
  options.skew_type = skew_type;
  options.additional_br_uncert = config.rel_eff_auto_base_rel_eff_uncert;

  // Find valid energy range based on contiguous channels with data
  const auto [min_valid_energy, max_valid_energy] = find_valid_energy_range( orig_foreground );


  if( fit_norm_peaks )
  {
    try
    {
      // See `getBackgroundRefLines()` in ReferenceLineInfo.cpp for
      const SandiaDecay::Nuclide * const u238  = db->nuclide( "U238" );
      const SandiaDecay::Nuclide * const ra226 = db->nuclide( "Ra226" );
      const SandiaDecay::Nuclide * const u235  = db->nuclide( "U235" );
      const SandiaDecay::Nuclide * const th232 = db->nuclide( "Th232" );
      const SandiaDecay::Nuclide * const k40   = db->nuclide( "K40" );

      assert( u238 && ra226 && u235 && th232 && k40 );

      // Denominators are for observation distance of 200cm, 5cm detector radius, and a 1 m radius sphere
      const vector<tuple<const SandiaDecay::Nuclide *,double, double>> nuc_activity{
        { u238,   0.0004653/410.2892,  5.0*u238->promptEquilibriumHalfLife() },   //make the 1001 keV have amp 0.0004653, before norm; the denom is integral at 1001 keV times BR of 1001.
        { ra226,  0.02515/17990.5430,  5.0*ra226->promptEquilibriumHalfLife() },  //make 609 keV have amp 0.02515, before norm
        { u235,   0.001482/14603.0156, 5.0*u235->promptEquilibriumHalfLife() },  //make 185 keV have amp 0.001482, before norm
        { th232,  0.02038/27897.2617,  5.0*th232->secularEquilibriumHalfLife() }, //make 2614 keV have amp 0.02038, before norm
        { k40,    0.1066/6523.8994,    0.0 }                                        //make 1460 keV have amp 0.1066, before norm
      };//nuc_activity

      // TODO: we coud find the 609, and if not the 1460, and if not 2614 keV peaks to use to normalize the starting
      //  relative activity, but we'll skip for the moment.
      //for( const shared_ptr<const PeakDef> &peak : auto_search_peaks )
      //{
      //}

      vector<RelActCalcAuto::NucInputInfo> norm_sources;

      for( const tuple<const SandiaDecay::Nuclide *,double, double> &info : nuc_activity )
      {
        RelActCalcAuto::NucInputInfo norm_src;
        norm_src.source  = get<0>(info);
        norm_src.age     = get<2>(info);
        norm_src.fit_age = false;
        //norm_src.starting_rel_act = ;
        norm_sources.push_back( norm_src );
      }

      PeakFitForNuclideConfig norm_config = config;
      norm_config.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      norm_config.rel_eff_eqn_order = 0;
      norm_config.phys_model_use_hoerl = false;

      norm_config.phys_model_self_atten.clear();
      auto self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>();
      //const char *soil_chem_formula = "H0.022019C0.009009O0.593577Al0.066067Si0.272289K0.01001Fe0.027029 d=1.6";
      //self_atten->material = make_share<Material>( MaterialDB::materialFromChemicalFormula( soil_chem_formula, db ) );
      self_atten->atomic_number = 10.4;
      self_atten->areal_density = (1.6 * PhysicalUnits::g / PhysicalUnits::cm3) * (100 * PhysicalUnits::cm);  //Transmission frac @2614 keV, though 100 cm soil is 0.2%
      self_atten->fit_atomic_number = false;
      self_atten->fit_areal_density = false;
      norm_config.phys_model_self_atten.push_back( self_atten );

      norm_config.phys_model_external_atten.clear();

      const bool many_peaks = true; //TODO: actually try to match the epaks up to estimate this!
      shared_ptr<RelActCalc::PhysicalModelShieldInput> ext_atten;
      if( many_peaks )
      {
        norm_config.phys_model_use_hoerl = true;
        ext_atten = make_shared<RelActCalc::PhysicalModelShieldInput>();

        // Have the external attenuator be concrete; we'll just use the AN/AD approx for the actual material, for the moment
        ext_atten->atomic_number = 11.3;
        ext_atten->areal_density = (2.3 * PhysicalUnits::g / PhysicalUnits::cm3) * (1.0 * PhysicalUnits::mm);
        ext_atten->fit_atomic_number = false;
        ext_atten->fit_areal_density = true;

        norm_config.phys_model_external_atten.push_back( ext_atten );
      }//if( many_peaks )

      auto norm_drf = drf;
      if( !norm_drf )
        norm_drf = isHPGe ? DetectorPeakResponse::getGenericHPGeDetector() : DetectorPeakResponse::getGenericNaIDetector();
      // We also have getGenericLaBrDetector(), getGenericCZTGeneralDetector(), getGenericCZTGoodDetector();

      const GammaClusteringSettings manual_settings = config.get_manual_clustering_settings();

      string norm_fallback_warning;
      const vector<RelActCalcAuto::RoiRange> norm_rois = estimate_initial_rois_using_relactmanual( auto_search_peaks,
        orig_foreground, norm_sources, norm_drf, isHPGe,
        fwhm_form, fwhm_coefficients,
        fwhm_lower_energy, fwhm_upper_energy,
        min_valid_energy, max_valid_energy,
        manual_settings, norm_config,
        norm_fallback_warning
      );

      vector<RelActCalcAuto::RoiRange> initial_src_norm_rois = input_rois;
      initial_src_norm_rois.insert( end(initial_src_norm_rois), begin(norm_rois), end(norm_rois) );
      vector<InitialRoi> initial_src_norm_info_rois;
      for( const RelActCalcAuto::RoiRange &roi : initial_src_norm_rois )
      {
        InitialRoi roi_info;
        roi_info.roi = roi;
        roi_info.center_energy = 0.5*(roi.upper_energy + roi.lower_energy);
        roi_info.fwhm = DetectorPeakResponse::peakResolutionFWHM(
                                                                 static_cast<float>(roi_info.center_energy), fwhm_form, fwhm_coefficients );

        float min_sigma_width, max_sigma_width;
        expected_peak_width_limits( roi_info.fwhm, isHPGe, orig_foreground, min_sigma_width, max_sigma_width );

        if( roi_info.fwhm < (min_sigma_width*PhysicalUnits::fwhm_nsigma) )
          roi_info.fwhm = min_sigma_width*PhysicalUnits::fwhm_nsigma;
        if( roi_info.fwhm > (max_sigma_width*PhysicalUnits::fwhm_nsigma) )
          roi_info.fwhm = max_sigma_width*PhysicalUnits::fwhm_nsigma;

        initial_src_norm_info_rois.push_back( roi_info );
      }//for( const RelActCalcAuto::RoiRange &roi : initial_src_norm_rois )

      options.rois = merge_rois( initial_src_norm_info_rois, config );

      RelActCalcAuto::RelEffCurveInput norm_rel_eff_curve;
      norm_rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      norm_rel_eff_curve.rel_eff_eqn_order = 0;
      norm_rel_eff_curve.nucs_of_el_same_age = false;
      norm_rel_eff_curve.phys_model_use_hoerl = false;
      norm_rel_eff_curve.phys_model_self_atten = self_atten;

      if( many_peaks )
      {
        norm_config.phys_model_use_hoerl = true;
        if( ext_atten )
          norm_config.phys_model_external_atten.push_back( ext_atten );
      }

      norm_rel_eff_curve.nuclides = norm_sources;

      options.rel_eff_curves.push_back( norm_rel_eff_curve );
    }catch( const std::exception &e )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
      result.error_message = "Error performing initial estimation of background peaks: " + string(e.what());
    }
  }//if( fit_norm_peaks )

  try
  {
    // Call RelActAuto::solve with provided options
    RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
      options, orig_foreground, orig_background, drf, auto_search_peaks
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
        no_ecal_opts, orig_foreground, orig_background, drf, auto_search_peaks
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

            if( should_debug_print() )
            {
              std::cerr << "First desperation attempt: using external shielding with AN="
                        << config.desperation_phys_model_atomic_number
                        << ", starting AD=" << config.desperation_phys_model_areal_density_g_per_cm2
                        << " g/cm2" << std::endl;
            }
          }
          catch( const std::exception &e )
          {
            if( should_debug_print() )
              std::cerr << "Failed to create desperation shielding: " << e.what() << std::endl;
            // Continue without shielding
          }
        }
        else
        {
          if( should_debug_print() )
            std::cerr << "First desperation attempt: not using external shielding" << std::endl;
        }

        curve.phys_model_use_hoerl = (options.rois.size() > 2);

        trial_solution = RelActCalcAuto::solve(
          desperation_opts, orig_foreground, orig_background, drf, auto_search_peaks
        );
      }//If( still a bad solution )
      
      if( (trial_solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success)
        && ( (solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success)
            || ((solution.m_chi2 / solution.m_dof) > (trial_solution.m_chi2 / trial_solution.m_dof)) ) )
      {
        if( should_debug_print() )
          std::cerr << "Abandoning fitting e-cal for nuclide" << std::endl;
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

    std::cout << "Initial RelActAuto solution (" << options.rel_eff_curves.front().name << "):" << std::endl;
    solution.print_summary( std::cout );
    std::cout << "Chi2/DOF = " << solution.m_chi2 << "/" << solution.m_dof << " = " << (solution.m_chi2 / solution.m_dof) << std::endl;

    // We may adjust energy calibration - sow we'll use `background` and `foreground` for this.
    shared_ptr<const SpecUtils::Measurement> foreground = orig_foreground;
    shared_ptr<const SpecUtils::Measurement> background = orig_background;

    // Iteratively refine ROIs using RelActAuto solutions
    // The idea is that each iteration provides a better relative efficiency estimate,
    // which allows us to better identify significant gamma lines and create better ROIs.
    {
      const size_t max_iterations = 3;
      size_t num_extra_allowed = 0; //If we switch to our "desperation" model type retry - we will increment this to 1.
      for( size_t iter = 0; iter < (max_iterations + num_extra_allowed); ++iter )
      {
        if( apply_energy_cal_between && config.fit_energy_cal )
        {
          const shared_ptr<SpecUtils::EnergyCalibration> fitted_cal = solution.get_adjusted_energy_cal();
          const shared_ptr<const SpecUtils::EnergyCalibration> orig_fg_cal = foreground->energy_calibration();

          shared_ptr<SpecUtils::Measurement> new_foreground = make_shared<SpecUtils::Measurement>( *foreground );
          new_foreground->set_energy_calibration( fitted_cal );
          foreground = new_foreground;

          if( background )
          {
            // Propagate the foreground cal change (orig_fg_cal -> fitted_cal) to
            // the background.  propogate_energy_cal_change requires orig and new
            // cals to be Polynomial or FRF; if LowerChannelEdge just use fitted_cal directly.
            shared_ptr<const SpecUtils::EnergyCalibration> new_bkg_cal;
            if( orig_fg_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge )
              new_bkg_cal = fitted_cal;
            else
              new_bkg_cal = EnergyCal::propogate_energy_cal_change( orig_fg_cal, fitted_cal, background->energy_calibration() );

            shared_ptr<SpecUtils::Measurement> new_background = make_shared<SpecUtils::Measurement>( *background );
            new_background->set_energy_calibration( new_bkg_cal );
            background = new_background;
          }
        }//if( apply_energy_cal_between && config.fit_energy_cal )


        vector<function<double(double)>> auto_rel_effs;
        vector<vector<tuple<RelActCalcAuto::SrcVariant,double,double>>> source_age_and_acts;
        for( size_t rel_eff_index = 0; rel_eff_index < solution.m_rel_activities.size(); ++rel_eff_index )
        {
          source_age_and_acts.push_back( {} );

          // Create rel_eff lambda from current RelActAuto solution
          auto_rel_effs.push_back( [&solution, rel_eff_index]( double energy ) -> double {
            return solution.relative_efficiency( energy, rel_eff_index );
          } );

          // Collect sources and activities from the current solution
          // m_rel_activities is a 2D vector: [rel_eff_curve_index][source_index]

          if( should_debug_print() )
            std::cout << "Collecting " << solution.m_rel_activities[rel_eff_index].size() << " sources from RelActAuto solution:" << std::endl;

          for( const RelActCalcAuto::NuclideRelAct &nuc_act : solution.m_rel_activities[rel_eff_index] )
          {
            if( RelActCalcAuto::is_null( nuc_act.source ) )
              continue;

            const double live_time_seconds = foreground->live_time();
            // RelActAuto's rel_activity is per second, need to multiply by live time for clustering
            const double activity_for_clustering = nuc_act.rel_activity * live_time_seconds;

            if( should_debug_print() )
              std::cout << "  " << nuc_act.name() << ": rel_activity=" << nuc_act.rel_activity
                   << ", live_time=" << live_time_seconds
                   << "s, activity_for_clustering=" << activity_for_clustering << std::endl;

            source_age_and_acts.back().emplace_back( nuc_act.source, nuc_act.age, activity_for_clustering );
          }//for( loop over solution.m_rel_activities[rel_eff_index] )
        }//for( size_t rel_eff_index = 0; rel_eff_index < ; ++rel_eff_index )

        // Use the FWHM coefficients passed to this function (computed from auto-search peaks
        // or DRF in fit_peaks_for_nuclides), rather than relying on solution.m_drf which may
        // not have valid FWHM info or may have incorrect values

        // Get auto clustering settings from config
        const GammaClusteringSettings auto_settings = config.get_auto_clustering_settings();

        // Cluster gammas using current solution's relative efficiency
        std::vector<RelActCalcAuto::RoiRange> refined_rois = cluster_gammas_to_rois(
            auto_rel_effs, source_age_and_acts, foreground,
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
            if( should_debug_print() )
              std::cerr << "Lost all ROIs, trying PhysicalModel as desperation attempt..." << std::endl;

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

                if( should_debug_print() )
                {
                  std::cerr << "Second desperation attempt: using external shielding with AN="
                            << config.desperation_phys_model_atomic_number
                            << ", starting AD=" << config.desperation_phys_model_areal_density_g_per_cm2
                            << " g/cm2" << std::endl;
                }
              }
              catch( const std::exception &e )
              {
                if( should_debug_print() )
                  std::cerr << "Failed to create desperation shielding: " << e.what() << std::endl;
                // Continue without shielding
              }
            }
            else
            {
              if( should_debug_print() )
                std::cerr << "Second desperation attempt: not using external shielding" << std::endl;
            }

            curve.phys_model_use_hoerl = (desperation_opts.rois.size() > 2);

            RelActCalcAuto::RelActAutoSolution desperation_solution = RelActCalcAuto::solve(
              desperation_opts, foreground, background, drf, auto_search_peaks
            );

            if( (desperation_solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success)
                && ((desperation_solution.m_chi2 / desperation_solution.m_dof)
                    < (solution.m_chi2 / solution.m_dof)) )
            {
              if( should_debug_print() )
                std::cerr << "PhysicalModel desperation solution succeeded and improved chi2/dof" << std::endl;
              
              if( !num_extra_allowed ) //Allow an extra iteration since we changed
                num_extra_allowed += 1;
              
              solution = desperation_solution;
              // Continue with another iteration using the new solution
              continue;
            }
            else
            {
              if( should_debug_print() )
                std::cerr << "PhysicalModel desperation solution did not improve result" << std::endl;
            }
          }//if( !using_physical_model )

          result.warnings.push_back( "Lost all ROIs while iterationing to refine solution - stopped early." );
          std::cerr << "Have lost all ROIs!  Halting iterations to refine solution." << std::endl;
          break;
        }//if( refined_rois.empty() )
        
        // Debug output: print refined ROIs with expected counts
        if( should_debug_print() )
        {
          std::cout << "Iteration " << iter << " refined ROIs:" << std::endl;
          for( size_t roi_idx = 0; roi_idx < refined_rois.size(); ++roi_idx )
          {
            const RelActCalcAuto::RoiRange &roi = refined_rois[roi_idx];
            const double roi_data_counts = foreground->gamma_integral(
                static_cast<float>(roi.lower_energy), static_cast<float>(roi.upper_energy) );
            std::cout << "  ROI " << roi_idx << ": [" << roi.lower_energy << " - " << roi.upper_energy
                 << "] keV, width=" << (roi.upper_energy - roi.lower_energy)
                 << " keV, data_counts=" << roi_data_counts << std::endl;
          }
        }

        // Check if ROIs changed significantly - if not, stop iterating
        if( rois_are_similar( refined_rois, solution.m_options.rois ) )
        {
          std::cout << "Iteration " << iter << ": ROIs are similar, stopping refinement" << std::endl;
          break;
        }

        std::cout << "Iteration " << iter << ": trying " << refined_rois.size() << " refined ROIs" << std::endl;

        // Re-run RelActAuto with refined ROIs
        RelActCalcAuto::Options refined_options = solution.m_options;
        refined_options.rois = refined_rois;

        RelActCalcAuto::RelActAutoSolution refined_solution
            = RelActCalcAuto::solve( refined_options, foreground, background, drf, auto_search_peaks );

        if( refined_solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        {
          std::cout << "Iteration " << iter << " failed: " << refined_solution.m_error_message << std::endl;
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
            std::cerr << "ERROR: RelActAuto returned overlapping ROIs[" << (i-1) << "] and [" << i << "]: "
                 << "[" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] vs "
                 << "[" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "]" << std::endl;
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
          std::cout << "Iteration " << iter << " did not improve filtered chi2/dof ("
               << old_chi2_dof << " -> " << new_chi2_dof << "), stopping" << std::endl;
          if( !new_insignificant_rois.empty() )
            std::cout << "  (" << new_insignificant_rois.size() << " ROIs had insignificant peaks)" << std::endl;
          break;
        }

        solution = std::move( refined_solution );
        std::cout << "Iteration " << iter << " improved: chi2/dof=" << new_chi2_dof
             << " (was " << old_chi2_dof << ")" << std::endl;
      }//for( size_t iter = 0; iter < max_iterations; ++iter )

      std::cout << "Final solution after refinement:" << std::endl;
      solution.print_summary( std::cout );
      std::cout << std::endl;

      // Print ROIs and sum fit peak areas for each ROI
      std::cout << "Solution ROIs and fit peak areas:" << std::endl;
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
          if( (std::fabs( peak_roi_lower - roi.lower_energy ) < 1.0)
             && (std::fabs( peak_roi_upper - roi.upper_energy ) < 1.0) )
          {
            sum_peak_area += peak.peakArea();
            ++num_peaks_in_roi;
          }
        }//for( loop over fit peaks )

        std::cout << "  ROI " << roi_index << ": [" << roi.lower_energy << ", " << roi.upper_energy << "] keV"
             << ", " << num_peaks_in_roi << " peaks, sum area = " << sum_peak_area << std::endl;
      }//for( loop over ROIs )
      std::cout << std::endl;
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
    if( should_debug_print() && !insignificant_roi_ranges.empty() )
      std::cout << "Peak filtering by ROI significance:" << std::endl;

    for( const PeakDef &peak : solution.m_peaks_without_back_sub )
    {
      const double mean = peak.mean();
      const double peak_roi_lower = peak.continuum()->lowerEnergy();
      const double peak_roi_upper = peak.continuum()->upperEnergy();
      bool in_insignificant_roi = false;

      for( const std::pair<double,double> &roi_range : insignificant_roi_ranges )
      {
        // Match if peak's ROI bounds are within 1 keV of the insignificant ROI bounds
        if( (std::fabs( peak_roi_lower - roi_range.first ) < 1.0)
           && (std::fabs( peak_roi_upper - roi_range.second ) < 1.0) )
        {
          in_insignificant_roi = true;
          break;
        }
      }

      if( in_insignificant_roi )
      {
        if( should_debug_print() )
        {
          std::cout << "  Filtered (insignificant ROI [" << peak.continuum()->lowerEnergy() << ", " << peak.continuum()->upperEnergy()
               << "] keV): peak at " << mean << " keV, area = " << peak.peakArea() << std::endl;
        }
      }else
      {
        bool mean_in_roi = false;
        for( size_t roi_index = 0; !mean_in_roi && (roi_index < solution.m_final_roi_ranges_in_spectrum_cal.size()); ++roi_index )
        {
          const auto pos = std::find( std::begin(final_insignificant_rois), std::end(final_insignificant_rois), roi_index );
          if( pos != std::end(final_insignificant_rois) )
            continue;
          const RelActCalcAuto::RoiRange &roi = solution.m_final_roi_ranges_in_spectrum_cal[roi_index];
          mean_in_roi = ((mean >= roi.lower_energy) && (mean <= roi.upper_energy));
        }

        if( should_debug_print() && !insignificant_roi_ranges.empty() )
        {
          std::cout << "  Kept (significant ROI [" << peak.continuum()->lowerEnergy() << ", " << peak.continuum()->upperEnergy()
               << "] keV): peak at " << mean << " keV, area = " << peak.peakArea() << (mean_in_roi ? " (was in ROI)" : " (skipping peak, not in a ROI)") << std::endl;
        }
        
        if( mean_in_roi )
          result.fit_peaks.push_back( peak );
      }
    }

    if( !insignificant_roi_ranges.empty() )
    {
      const size_t num_filtered = solution.m_peaks_without_back_sub.size() - result.fit_peaks.size();
      std::cout << "Filtered out " << num_filtered << " peaks from "
           << insignificant_roi_ranges.size() << " ROIs without significant chi2 improvement" << std::endl;
    }
    
    

    result.solution = std::move( solution );

    // Combine overlapping peaks within ROIs
    // First, preserve the uncombined peaks, then create combined version
    result.uncombined_fit_peaks = result.fit_peaks;
    result.fit_peaks = combine_overlapping_peaks_in_rois( result.uncombined_fit_peaks );

    if( should_debug_print() && (result.fit_peaks.size() != result.uncombined_fit_peaks.size()) )
    {
      std::cout << "Combined " << result.uncombined_fit_peaks.size() << " peaks into "
           << result.fit_peaks.size() << " peaks" << std::endl;
    }

    // Compute observable peaks - peaks that users can visually see in the spectrum
    result.observable_peaks = compute_observable_peaks( result.fit_peaks, foreground, config );


    if( apply_energy_cal_between && config.fit_energy_cal )
    {
      // Peaks are currently in foreground's (fitted) energy cal; translate them
      // back to the original foreground's energy cal for display.
      const shared_ptr<const SpecUtils::EnergyCalibration> fitted_cal = foreground->energy_calibration();
      const shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = orig_foreground->energy_calibration();

      if( fitted_cal && orig_cal && (*fitted_cal != *orig_cal) )
      {
        auto translate_peaks = [&fitted_cal, &orig_cal]( vector<PeakDef> &peaks )
        {
          deque<shared_ptr<const PeakDef>> tmp_peaks;
          for( const PeakDef &p : peaks )
            tmp_peaks.push_back( make_shared<const PeakDef>( p ) );

          const deque<shared_ptr<const PeakDef>> translated
              = EnergyCal::translatePeaksForCalibrationChange( tmp_peaks, fitted_cal, orig_cal );

          peaks.clear();
          for( const shared_ptr<const PeakDef> &p : translated )
            peaks.push_back( *p );
        };

        translate_peaks( result.fit_peaks );
        translate_peaks( result.uncombined_fit_peaks );
        translate_peaks( result.observable_peaks );
      }
    }//if( apply_energy_cal_between && config.fit_energy_cal )

    if( should_debug_print() && (result.observable_peaks.size() != result.fit_peaks.size()) )
    {
      std::cout << "Observable peaks: " << result.observable_peaks.size() << " of "
           << result.fit_peaks.size() << " fit_peaks" << std::endl;
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
            std::cerr << "WARNING: Duplicate peaks at " << peak_i.mean() << " keV with different ROIs: "
                 << "[" << peak_i.continuum()->lowerEnergy() << ", " << peak_i.continuum()->upperEnergy() << "] vs "
                 << "[" << peak_j.continuum()->lowerEnergy() << ", " << peak_j.continuum()->upperEnergy() << "]" << std::endl;
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

  
PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const SpecUtils::Measurement> &background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const Wt::WFlags<FitSrcPeaksOptions> options,
  const PeakFitForNuclideConfig &config,
  const bool isHPGe )
{
  
  std::vector<RelActCalcAuto::NucInputInfo> base_nuclides;
  base_nuclides.reserve( sources.size() );
  
  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    RelActCalcAuto::NucInputInfo nuc_info;
    nuc_info.age = get_source_age( src, -1.0 );
    nuc_info.source = src;
    nuc_info.fit_age = false;  // not currently exposed in UI
    base_nuclides.push_back( nuc_info );
  }
  
  return fit_peaks_for_nuclides( auto_search_peaks, foreground, base_nuclides,
                                background, drf_input, options, config, isHPGe );
}
  
  
PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::shared_ptr<const SpecUtils::Measurement> &long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const Wt::WFlags<FitSrcPeaksOptions> options,
  const PeakFitForNuclideConfig &config,
  const bool isHPGe )
{
  /* These functions needs the following options addresses/implemented from FitSrcPeaksOptions:
   -[ ] DoNotReplaceExistingPeaksForSource = 0x01,
   -[ ] ExistingPeaksAsFreePeak = 0x02,
   -[ ] DoNotVaryEnergyCal = 0x04,
   -[x] DoNotRefineEnergyCal = 0x08,
   -[x] FitNormBkgrndPeaks = 0x10,
   -[ ] FitNormBkgrndPeaksDontUse = 0x20
   
   Also need to take care of case where one of the sources is a NORM nuclide, and we're being asked to fit norm
   */
  
  
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

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
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
    std::vector<float> fwhm_coefficients, fwhm_uncerts;

    if( !drf || !drf->isValid() || !drf->hasResolutionInfo() || (auto_search_peaks.size() > 6) )
    {
      // If we have peaks, estimate FWHM from peak widths. Otherwise fall back to a generic detector.
      if( !auto_search_peaks.empty() )
      {
        const int num_auto_peaks = static_cast<int>(auto_search_peaks.size());
        int sqrtEqnOrder = (std::min)( 6, num_auto_peaks / (1 + (num_auto_peaks > 3)) );
        if( auto_search_peaks.size() < 3 )
          sqrtEqnOrder = static_cast<int>( auto_search_peaks.size() );

        std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> auto_search_peaks_dq
          = std::make_shared<const std::deque<std::shared_ptr<const PeakDef>>>( begin(auto_search_peaks), end(auto_search_peaks) );

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
        // With no peaks and no DRF resolution info, use generic detector resolution coefficients.
        if( isHPGe )
          drf = DetectorPeakResponse::getGenericHPGeDetector();
        else
          drf = DetectorPeakResponse::getGenericNaIDetector();

        if( drf && drf->isValid() && drf->hasResolutionInfo() )
        {
          fwhmFnctnlForm = drf->resolutionFcnType();
          fwhm_coefficients = drf->resolutionFcnCoefficients();
          lower_fwhm_energy = drf->lowerEnergy();
          upper_fwhm_energy = drf->upperEnergy();
          fallback_warning = "No peaks were available to estimate resolution; using generic detector FWHM parameters.";
        }
      }//if( !auto_search_peaks.empty() ) / else
    }else
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

    // Find valid energy range based on contiguous channels with data
    const auto [min_valid_energy, max_valid_energy] = find_valid_energy_range( foreground );

    /*
    // Determine energy range for gamma lines
    double highest_energy_gamma = 0.0, lowest_energy_gamma = std::numeric_limits<double>::max();

    
    for( const RelActCalcAuto::NucInputInfo &src : sources )
    {
      const std::vector<SandiaDecay::EnergyRatePair> photons
          = get_source_photons( src.source, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits, src.age );
      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        highest_energy_gamma = (std::max)( highest_energy_gamma, photon.energy );
        lowest_energy_gamma = (std::min)( lowest_energy_gamma, photon.energy );
      }
    }

    lowest_energy_gamma = (std::max)( lowest_energy_gamma - (isHPGe ? 5 : 25), (double)foreground->gamma_energy_min() );
    highest_energy_gamma = (std::min)( highest_energy_gamma + (isHPGe ? 5 : 25), (double)foreground->gamma_energy_max() );
    */

    // Step 2, 3 & 4: Estimate initial ROIs using RelActManual with multiple fallbacks
    // This function internally:
    // - Converts auto_search_peaks to RelActManual format and matches to sources
    // - Falls back to estimate_initial_rois_without_peaks() if no peaks match
    // - Fits relative efficiency curve and clusters gammas into ROIs
    // - Falls back to estimate_initial_rois_fallback() if RelActManual fails
    const GammaClusteringSettings manual_settings = config.get_manual_clustering_settings();

    const std::vector<RelActCalcAuto::RoiRange> initial_rois = estimate_initial_rois_using_relactmanual(
      auto_search_peaks, foreground, sources, drf, isHPGe,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      min_valid_energy, max_valid_energy, manual_settings,
      config, fallback_warning );

    // Call RelActAuto with initial_rois
    result = fit_peaks_for_nuclide_relactauto(
      auto_search_peaks, foreground, sources,
      initial_rois, long_background, drf, options, config,
      fwhmFnctnlForm, fwhm_coefficients, isHPGe,
      lower_fwhm_energy, upper_fwhm_energy
    );

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

}//namespace FitPeaksForNuclides
