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

#include <map>
#include <cmath>
#include <mutex>
#include <chrono>
#include <string>
#include <fstream>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitSpecImp.h"
#include "InterSpec/FitPeaksForNuclides.h"

#include "PeakFitImprove.h"
#include "CandidatePeak_GA.h"

using namespace std;


std::mutex EA::mtx_rand;


namespace
{
std::function<double(const FindCandidateSettings &)> ns_ga_eval_fcn;

const std::vector<DataSrcInfo> *ns_input_srcs = nullptr;

// Cache for find_valid_energy_range results, keyed by Measurement pointer.
// The same spectra are evaluated for every GA individual, so caching avoids
// re-running the expensive range computation ~500 * N_sources times per generation.
std::mutex ns_valid_range_cache_mutex;
std::map<const SpecUtils::Measurement *, std::pair<size_t,size_t>> ns_valid_range_cache;

std::ofstream sm_output_file;

bool sm_set_best_genes = false;
CandidatePeak_GA::CandidatePeakSolution sm_best_genes;
double sm_best_total_cost = 1.0E99;

bool sm_has_been_called = false;
}


namespace CandidatePeak_GA
{

/** Estimates peak sigma by fitting a parabola in log-space (i.e., a Gaussian) to the raw
 spectrum data near the candidate peak center, after subtracting a linear continuum.

 Returns the estimated sigma in keV, or -1.0 if the fit fails.

 Algorithm:
 1. Estimate a linear continuum from channels just outside the ROI
 2. Find the peak channel (max net counts near the candidate mean)
 3. Determine the fit window: channels where net counts >= 15% of peak height
    This focuses on the peak top where the Gaussian shape is best defined
 4. Fit ln(net_counts) = a + b*x + c*x^2 (weighted by net counts)
 5. Extract sigma = sqrt(-1/(2c))
 */
/** Define to 1 to use two-pass log-parabola Gaussian fit for candidate peak sigma estimation,
 *  or 0 to use the original zero-crossing method.
 *  The Gaussian fit method reduces FWHM estimation error by ~49% vs zero-crossings.
 */
#define USE_TWO_PASS_SIMPLE_GAUSS_FIT 1


/** Fits a log-parabola (ln(y) = a + b*x + c*x^2) to continuum-subtracted spectrum data
 *  near a candidate peak to estimate sigma.
 *
 *  The curvature coefficient c (which determines sigma) is robust to continuum errors because
 *  it depends on the relative shape of the peak, not its absolute height.  The intercept a
 *  (which determines amplitude/area via exp(a)) is sensitive to even small continuum biases,
 *  so this method gives much better FWHM estimates than the SG zero-crossing method, but
 *  worse area estimates than the second-derivative amplitude method.
 *
 *  @param spectrum The raw spectrum channel counts
 *  @param energy_cal The energy calibration
 *  @param nchannel Total number of channels
 *  @param peak_ch Channel of the peak maximum
 *  @param peak_energy Energy of the peak maximum (keV)
 *  @param cont_ref_ch Reference channel for continuum interpolation (the lower bound)
 *  @param lower_cont Continuum level (counts/channel) at the lower reference
 *  @param cont_slope Continuum slope (counts/channel per channel)
 *  @param area_out If non-null, receives the estimated Gaussian peak area (counts)
 *  @return Estimated sigma in keV, or -1.0 on failure
 */
static double fit_log_parabola( const std::vector<float> &spectrum,
                                const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                                const size_t nchannel,
                                const size_t peak_ch,
                                const double peak_energy,
                                const size_t cont_ref_ch,
                                const double lower_cont,
                                const double cont_slope,
                                double *area_out = nullptr )
{
  // Find net peak height with current continuum
  const double cont_at_peak = lower_cont + cont_slope * static_cast<double>( peak_ch - cont_ref_ch );
  const double peak_net = spectrum[peak_ch] - cont_at_peak;

  if( peak_net <= 1.0 )
    return -1.0;

  // Determine fit window: channels with net counts >= 5% of peak net counts.
  // This focuses on the peak top where the Gaussian shape is well-defined.
  // We use 5% (not higher) because HPGe peaks can be as narrow as 3-4 channels.
  const double threshold = 0.05 * peak_net;

  size_t fit_lower = peak_ch, fit_upper = peak_ch;

  while( fit_lower > 0 )
  {
    const double cont_at = lower_cont + cont_slope * static_cast<double>( fit_lower - 1 - cont_ref_ch );
    const double net = spectrum[fit_lower - 1] - cont_at;
    if( net < threshold )
      break;
    --fit_lower;
  }

  while( (fit_upper + 1) < nchannel )
  {
    const double cont_at = lower_cont + cont_slope * static_cast<double>( fit_upper + 1 - cont_ref_ch );
    const double net = spectrum[fit_upper + 1] - cont_at;
    if( net < threshold )
      break;
    ++fit_upper;
  }

  // Need at least 3 channels for a 3-parameter fit
  if( (fit_upper - fit_lower) < 2 )
    return -1.0;

  // Weighted least-squares fit of ln(y) = a + b*x + c*x^2
  // where x is energy relative to peak energy, y is continuum-subtracted counts
  // Weight by net counts (variance of ln(y) ~ 1/y for Poisson data)
  double S = 0, Sx = 0, Sx2 = 0, Sx3 = 0, Sx4 = 0;
  double Sy = 0, Sxy = 0, Sx2y = 0;
  size_t n_valid = 0;

  for( size_t ch = fit_lower; ch <= fit_upper; ++ch )
  {
    const double cont_at = lower_cont + cont_slope * static_cast<double>( ch - cont_ref_ch );
    const double net_counts = spectrum[ch] - cont_at;

    if( net_counts <= 0.5 )
      continue;

    const double energy_ch = energy_cal->energy_for_channel( static_cast<double>( ch ) + 0.5 );
    const double x = energy_ch - peak_energy;
    const double ln_y = std::log( net_counts );
    const double w = net_counts;

    const double x2 = x * x;
    S += w;
    Sx += w * x;
    Sx2 += w * x2;
    Sx3 += w * x2 * x;
    Sx4 += w * x2 * x2;
    Sy += w * ln_y;
    Sxy += w * x * ln_y;
    Sx2y += w * x2 * ln_y;
    ++n_valid;
  }

  if( n_valid < 3 || S <= 0.0 )
    return -1.0;

  // Solve 3x3 normal equations via Cramer's rule for a, b, c
  const double d00 = S,   d01 = Sx,  d02 = Sx2;
  const double d10 = Sx,  d11 = Sx2, d12 = Sx3;
  const double d20 = Sx2, d21 = Sx3, d22 = Sx4;

  const double det = d00 * (d11 * d22 - d12 * d21)
                   - d01 * (d10 * d22 - d12 * d20)
                   + d02 * (d10 * d21 - d11 * d20);

  if( std::fabs( det ) < 1.0e-30 )
    return -1.0;

  const double det_c = d00 * (d11 * Sx2y - Sxy * d21)
                     - d01 * (d10 * Sx2y - Sxy * d20)
                     + Sy  * (d10 * d21  - d11 * d20);

  const double c = det_c / det;

  // For a Gaussian: ln(y) = ln(A) - x^2 / (2*sigma^2), so c = -1/(2*sigma^2)
  if( c >= 0.0 )
    return -1.0;

  const double sigma_keV = std::sqrt( -1.0 / (2.0 * c) );

  if( sigma_keV < 0.1 || sigma_keV > 100.0 )
    return -1.0;

  // Optionally compute peak area: a = ln(peak_height), area = exp(a) * sigma * sqrt(2*pi) / keV_per_ch
  // The fit is in counts-per-channel vs energy, so exp(a) is counts/channel at peak center.
  // Total counts = integral of Gaussian / channel_width = height * sigma * sqrt(2*pi) / keV_per_ch
  if( area_out )
  {
    const double det_a = Sy  * (d11 * d22 - d12 * d21)
                       - d01 * (Sxy * d22 - Sx2y * d20)
                       + d02 * (Sxy * d21 - Sx2y * d11);
    const double a = det_a / det;
    const double peak_height = std::exp( a );  // counts at peak center (per channel)

    // Estimate keV per channel at the peak
    const double keV_per_ch = energy_cal->energy_for_channel( static_cast<double>( peak_ch ) + 1.0 )
                            - energy_cal->energy_for_channel( static_cast<double>( peak_ch ) );

    if( keV_per_ch > 0.0 )
      *area_out = peak_height * sigma_keV * std::sqrt( 2.0 * boost::math::constants::pi<double>() ) / keV_per_ch;
    else
      *area_out = -1.0;
  }

  return sigma_keV;
}//fit_log_parabola(...)


double estimate_sigma_from_raw_data( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                     const double peak_mean_keV,
                                     const double roi_lower_energy,
                                     const double roi_upper_energy,
                                     double *area_out = nullptr )
{
  const std::shared_ptr<const SpecUtils::EnergyCalibration> energy_cal = data->energy_calibration();
  if( !energy_cal || !energy_cal->valid() || !data->gamma_counts() )
    return -1.0;

  const std::vector<float> &spectrum = *data->gamma_counts();
  const size_t nchannel = spectrum.size();
  if( nchannel < 20 )
    return -1.0;

  // Find channels for ROI boundaries
  const size_t roi_lower_ch = static_cast<size_t>(
      std::max( 0.0, std::floor( energy_cal->channel_for_energy( roi_lower_energy ) ) ) );
  const size_t roi_upper_ch = static_cast<size_t>(
      std::min( static_cast<double>( nchannel - 1 ),
                std::ceil( energy_cal->channel_for_energy( roi_upper_energy ) ) ) );

  if( roi_upper_ch <= (roi_lower_ch + 2) )
    return -1.0;

  // --- Pass 1: Estimate continuum from channels just outside the ROI ---
  const size_t cont_channels = std::min( size_t(3), (roi_upper_ch - roi_lower_ch) / 4 + 1 );

  double lower_cont = 0.0, upper_cont = 0.0;
  size_t lower_cont_n = 0, upper_cont_n = 0;
  for( size_t i = 0; i < cont_channels; ++i )
  {
    if( roi_lower_ch >= (cont_channels - i) )
    {
      lower_cont += spectrum[roi_lower_ch - cont_channels + i];
      ++lower_cont_n;
    }

    const size_t ui = roi_upper_ch + 1 + i;
    if( ui < nchannel )
    {
      upper_cont += spectrum[ui];
      ++upper_cont_n;
    }
  }

  if( lower_cont_n == 0 || upper_cont_n == 0 )
    return -1.0;

  lower_cont /= lower_cont_n;
  upper_cont /= upper_cont_n;

  const double roi_span = static_cast<double>( roi_upper_ch - roi_lower_ch );
  const double cont_slope = (upper_cont - lower_cont) / roi_span;

  // Find the actual peak channel (max net counts) within a reasonable range of the candidate mean
  const double mean_channel = energy_cal->channel_for_energy( peak_mean_keV );
  const size_t search_half = std::max( size_t(3), (roi_upper_ch - roi_lower_ch) / 3 );
  const size_t search_lower = (static_cast<size_t>( mean_channel ) > search_half)
                            ? (static_cast<size_t>( mean_channel ) - search_half) : 0;
  const size_t search_upper = std::min( static_cast<size_t>( mean_channel ) + search_half, nchannel - 1 );

  size_t peak_ch = static_cast<size_t>( std::max( 0.0, std::round( mean_channel ) ) );
  double peak_net = -1.0e9;

  for( size_t ch = search_lower; ch <= search_upper; ++ch )
  {
    const double cont_at_ch = lower_cont + cont_slope * static_cast<double>( ch - roi_lower_ch );
    const double net = spectrum[ch] - cont_at_ch;
    if( net > peak_net )
    {
      peak_net = net;
      peak_ch = ch;
    }
  }

  if( peak_net <= 1.0 )
    return -1.0;

  const double peak_energy = energy_cal->energy_for_channel( static_cast<double>( peak_ch ) + 0.5 );

  // First-pass fit using ROI-edge continuum
  double pass1_area = -1.0;
  const double pass1_sigma = fit_log_parabola( spectrum, energy_cal, nchannel,
                                                peak_ch, peak_energy,
                                                roi_lower_ch, lower_cont, cont_slope,
                                                area_out ? &pass1_area : nullptr );

  if( pass1_sigma <= 0.0 )
    return -1.0;

  // --- Pass 2: Symmetric fit window to mitigate neighbor contamination ---
  // Use pass1 sigma to define a symmetric window of ~2.5*sigma on each side.
  // Only accept pass2 if it gives a smaller sigma (indicating pass1 was inflated
  // by a neighbor) and the reduction is not too drastic (> 60% of pass1, otherwise
  // it's likely noise-driven for weak peaks).
  const double keV_per_ch = (roi_upper_energy - roi_lower_energy) / roi_span;
  const double sigma_ch = pass1_sigma / keV_per_ch;

  const size_t sym_half = static_cast<size_t>( std::max( 1.0, std::round( 2.5 * sigma_ch ) ) );
  const size_t sym_lower = (peak_ch > sym_half) ? (peak_ch - sym_half) : 0;
  const size_t sym_upper = std::min( peak_ch + sym_half, nchannel - 1 );

  if( (sym_upper - sym_lower) >= 2 )
  {
    // Estimate local continuum from channels just outside the symmetric window
    const size_t sym_cont_ch = std::min( size_t(3), sym_half );

    double sym_lower_cont = 0.0, sym_upper_cont = 0.0;
    size_t sym_lower_n = 0, sym_upper_n = 0;

    for( size_t i = 0; i < sym_cont_ch; ++i )
    {
      if( sym_lower >= (sym_cont_ch - i) )
      {
        sym_lower_cont += spectrum[sym_lower - sym_cont_ch + i];
        ++sym_lower_n;
      }

      const size_t ui = sym_upper + 1 + i;
      if( ui < nchannel )
      {
        sym_upper_cont += spectrum[ui];
        ++sym_upper_n;
      }
    }

    if( sym_lower_n > 0 && sym_upper_n > 0 )
    {
      sym_lower_cont /= sym_lower_n;
      sym_upper_cont /= sym_upper_n;

      const double sym_span = static_cast<double>( sym_upper - sym_lower );
      const double sym_slope = (sym_upper_cont - sym_lower_cont) / sym_span;

      double pass2_area = -1.0;
      const double pass2_sigma = fit_log_parabola( spectrum, energy_cal, nchannel,
                                                    peak_ch, peak_energy,
                                                    sym_lower, sym_lower_cont, sym_slope,
                                                    area_out ? &pass2_area : nullptr );

      // Accept pass2 if it gives a smaller sigma (indicating pass1 was inflated
      // by neighbor contamination extending the fit window asymmetrically)
      if( (pass2_sigma > 0.0) && (pass2_sigma < pass1_sigma) )
      {
        if( area_out )
          *area_out = pass2_area;
        return pass2_sigma;
      }
    }
  }

  if( area_out )
    *area_out = pass1_area;

  return pass1_sigma;
}//estimate_sigma_from_raw_data(...)


std::vector<PeakDef> find_candidate_peaks( const std::shared_ptr<const SpecUtils::Measurement> &data,
                          size_t start_channel,
                          size_t end_channel,
                          const FindCandidateSettings &settings )
{
  // Use find_valid_energy_range to determine the valid energy bounds, replacing the
  //  ad-hoc turn-on detection in PeakFitSpec::find_candidate_peaks when start_channel == 0.
  //  Results are cached since the same spectra are evaluated for every GA individual.
  if( (start_channel == 0) && data && (data->num_gamma_channels() > 16) )
  {
    const SpecUtils::Measurement *key = data.get();
    std::pair<size_t,size_t> cached_channels;
    bool found_in_cache = false;

    {
      std::lock_guard<std::mutex> lock( ns_valid_range_cache_mutex );
      const auto it = ns_valid_range_cache.find( key );
      if( it != end(ns_valid_range_cache) )
      {
        cached_channels = it->second;
        found_in_cache = true;
      }
    }

    if( !found_in_cache )
    {
      try
      {
        const auto [min_energy, max_energy] = FitPeaksForNuclides::find_valid_energy_range( data );

        if( min_energy < max_energy )
        {
          cached_channels.first = data->find_gamma_channel( static_cast<float>( min_energy ) );
          cached_channels.second = data->find_gamma_channel( static_cast<float>( max_energy ) );
        }else
        {
          cached_channels = { 0, 0 };
        }
      }catch( std::exception & )
      {
        // find_valid_energy_range or find_gamma_channel may throw for unusual energy
        //  calibrations (e.g., full range fraction that fails to converge); fall back
        //  to the ad-hoc turn-on detection in PeakFitSpec::find_candidate_peaks.
        cached_channels = { 0, 0 };
      }

      std::lock_guard<std::mutex> lock( ns_valid_range_cache_mutex );
      ns_valid_range_cache[key] = cached_channels;
    }//if( !found_in_cache )

    if( cached_channels.first > 0 || cached_channels.second > 0 )
    {
      start_channel = cached_channels.first;

      if( (end_channel == 0) || (end_channel >= (data->num_gamma_channels() - 1)) )
        end_channel = cached_channels.second;
    }
  }//if( start_channel == 0 )

  return PeakFitSpec::find_candidate_peaks( data, start_channel, end_channel, settings );
}//find_candidate_peaks(...)


CandidatePeakScore calculate_candidate_peak_score_for_source( const vector<PeakDef> &detected_peaks,
                                                              const vector<ExpectedPhotopeakInfo> &expected_photopeaks )
{
  vector<tuple<float,float,float>> detected_peaks_tuples; //{mean, sigma, amplitude}
  for( const PeakDef &p : detected_peaks )
    detected_peaks_tuples.emplace_back( p.mean(), p.sigma(), p.amplitude() );

  double score = 0.0;
  size_t num_detected_expected = 0, num_detected_not_expected = 0;
  size_t num_possibly_accepted_but_not_detected = 0, num_def_wanted_but_not_detected = 0;
  size_t num_def_wanted_detected = 0;

  // First, go through expected peaks, and match up to candidates
  vector<ExpectedPhotopeakInfo> possibly_expected_but_not_detected, expected_and_was_detected, def_expected_but_not_detected;
  possibly_expected_but_not_detected.reserve( expected_photopeaks.size() );
  def_expected_but_not_detected.reserve( expected_photopeaks.size() );
  expected_and_was_detected.reserve( expected_photopeaks.size() );

  vector<bool> detected_matched_to_expected( detected_peaks_tuples.size(), false );

  for( size_t expected_index = 0; expected_index < expected_photopeaks.size(); ++expected_index )
  {
    const ExpectedPhotopeakInfo &expected = expected_photopeaks[expected_index];

    vector<pair<tuple<float,float,float>,size_t>> detected_matching_expected; //{mean, sigma, amplitude}
    for( size_t det_index = 0; det_index < detected_peaks_tuples.size(); ++det_index )
    {
      const tuple<float,float,float> &det_peak = detected_peaks_tuples[det_index];
      const float mean = get<0>(det_peak);
      const float det_sigma = get<1>(det_peak);

      // Match if detected peak mean is within expected ROI, or if the detected
      //  peak's ±1 sigma range overlaps with the expected peak's ROI.
      //  This handles smoothing shifting detected peak means.
      const bool mean_in_roi = (mean > expected.roi_lower) && (mean < expected.roi_upper);
      const bool det_roi_overlaps = ((mean + det_sigma) > expected.roi_lower)
                                 && ((mean - det_sigma) < expected.roi_upper);
      if( mean_in_roi || det_roi_overlaps )
      {
        detected_matching_expected.push_back( make_pair(det_peak,det_index) );
        detected_matched_to_expected[det_index] = true;
      }
    }//for( size_t det_index = 0; det_index < detected_peaks_tuples.size(); ++det_index )

    if( detected_matching_expected.empty() )
    {
      num_possibly_accepted_but_not_detected += 1;

      if( (expected.nsigma_over_background > JudgmentFactors::def_want_nsigma)
         && (expected.peak_area > JudgmentFactors::min_def_wanted_counts) )
      {
        num_def_wanted_but_not_detected += 1;
        def_expected_but_not_detected.push_back( expected );
        if( PeakFitImprove::debug_printout )
        {
          cerr << "Def. wanted m=" << expected.effective_energy
               << " keV, nsigma=" << expected.nsigma_over_background
               << ", area=" << expected.peak_area
               << ", FWHM=" << expected.effective_fwhm
               << ", ROI=[" << expected.roi_lower << ", " << expected.roi_upper << "]"
               << ", but didnt find." << endl;
        }
      }else
      {
        possibly_expected_but_not_detected.push_back( expected );
      }
    }else
    {
      // we got a match
      num_detected_expected += 1;
      expected_and_was_detected.push_back( expected );

      if( (expected.nsigma_over_background > JudgmentFactors::def_want_nsigma)
         && (expected.peak_area > JudgmentFactors::min_def_wanted_counts) )
      {
        num_def_wanted_detected += 1;
      }

      if( expected.nsigma_over_background > JudgmentFactors::def_want_nsigma )
      {
        score += 1.0;
      }else if( expected.nsigma_over_background > JudgmentFactors::lower_want_nsigma )
      {
        const double amount_short = JudgmentFactors::def_want_nsigma - expected.nsigma_over_background;
        const double fraction_short = amount_short / (JudgmentFactors::def_want_nsigma - JudgmentFactors::lower_want_nsigma);
        assert( (fraction_short >= 0.0) && (fraction_short <= 1.0) );
        score += (1 - fraction_short);
      }else
      {
        score += 0.0; //No punishment, but no reward.
      }
    }//if( detected_matching_expected.empty() ) / else
  }//for( loop over expected_photopeaks )

  assert( expected_and_was_detected.size() == num_detected_expected );

  vector<tuple<float,float,float>> detected_expected, detected_not_expected;
  detected_expected.reserve( expected_photopeaks.size() );
  detected_not_expected.reserve( detected_peaks_tuples.size() );

  for( size_t i = 0; i < detected_matched_to_expected.size(); ++i )
  {
    if( detected_matched_to_expected[i] )
      detected_expected.push_back( detected_peaks_tuples[i] );
    else
      detected_not_expected.push_back( detected_peaks_tuples[i] );
  }//for( size_t i = 0; i < detected_matched_to_expected.size(); ++i )

  num_detected_not_expected = detected_not_expected.size();

  // The 511 Annih. line isnt in our truth line, unless the source makes them, so lets
  //  remove this one from the score
  for( const auto &p : detected_not_expected )
  {
    if( (get<0>(p) > 508) && (get<0>(p) < 514) )
    {
      num_detected_not_expected -= 1;
      break;
    }
  }

  if( PeakFitImprove::debug_printout && !detected_not_expected.empty() )
  {
    cerr << "Extra (unexpected) peaks detected:" << endl;
    for( const auto &p : detected_not_expected )
    {
      cerr << "  Extra peak: m=" << get<0>(p) << " keV, sigma=" << get<1>(p)
           << ", amp=" << get<2>(p) << endl;
    }
  }

  score -= JudgmentFactors::found_extra_punishment * num_detected_not_expected;

  CandidatePeakScore result;
  result.score = score;
  result.num_peaks_found = num_detected_expected;
  result.num_def_wanted_not_found = num_def_wanted_but_not_detected;
  result.num_def_wanted_peaks_found = num_def_wanted_detected;
  result.num_possibly_accepted_peaks_not_found = num_possibly_accepted_but_not_detected;
  result.num_extra_peaks = num_detected_not_expected;
  result.detected_expected = std::move(detected_expected);
  result.detected_not_expected = std::move(detected_not_expected);
  result.possibly_expected_but_not_detected = std::move(possibly_expected_but_not_detected);
  result.expected_and_was_detected = std::move(expected_and_was_detected);
  result.def_expected_but_not_detected = std::move(def_expected_but_not_detected);
  return result;
}//calculate_candidate_peak_score_for_source


void correct_score_for_escape_peaks( CandidatePeakScore &score,
                                     const vector<ExpectedPhotopeakInfo> &expected_photopeaks )
{
  // Single and double escape energies (electron rest mass energy)
  const double single_escape_offset = 510.9989; // keV
  const double double_escape_offset = 2.0 * single_escape_offset; // 1021.9978 keV

  // Tolerance for matching escape peak to parent
  const double energy_tolerance = 1.0; // keV

  // Minimum parent-to-escape area ratio to consider it a likely escape peak (under a few MeV, anyway)
  const double min_parent_ratio = 2.0;

  // Iterate through def_expected_but_not_detected to find escape peaks
  vector<ExpectedPhotopeakInfo> corrected_def_expected;
  corrected_def_expected.reserve( score.def_expected_but_not_detected.size() );

  size_t num_escape_peaks_removed = 0;

  for( const ExpectedPhotopeakInfo &not_found : score.def_expected_but_not_detected )
  {
    const double escape_energy = not_found.effective_energy;
    bool is_escape_peak = false;

    // Check for single escape (parent at +511 keV)
    const double se_parent_energy = escape_energy + single_escape_offset;
    for( const ExpectedPhotopeakInfo &potential_parent : expected_photopeaks )
    {
      if( fabs(potential_parent.effective_energy - se_parent_energy) < energy_tolerance )
      {
        // Check if parent peak is significantly larger
        if( potential_parent.peak_area > (min_parent_ratio * not_found.peak_area) )
        {
          is_escape_peak = true;
          if( PeakFitImprove::debug_printout )
          {
            cerr << "Identified S.E. peak at " << escape_energy << " keV "
                 << "(parent at " << potential_parent.effective_energy << " keV)" << endl;
          }
          break;
        }
      }
    }

    // Check for double escape (parent at +1022 keV) if not already identified as S.E.
    if( !is_escape_peak )
    {
      const double de_parent_energy = escape_energy + double_escape_offset;
      for( const ExpectedPhotopeakInfo &potential_parent : expected_photopeaks )
      {
        if( fabs(potential_parent.effective_energy - de_parent_energy) < energy_tolerance )
        {
          // Check if parent peak is significantly larger
          if( potential_parent.peak_area > (min_parent_ratio * not_found.peak_area) )
          {
            is_escape_peak = true;
            if( PeakFitImprove::debug_printout )
            {
              cerr << "Identified D.E. peak at " << escape_energy << " keV "
                   << "(parent at " << potential_parent.effective_energy << " keV)" << endl;
            }
            break;
          }
        }
      }
    }

    // Only keep non-escape peaks in the corrected list
    if( !is_escape_peak )
    {
      corrected_def_expected.push_back( not_found );
    }
    else
    {
      num_escape_peaks_removed++;
    }
  }

  // Update the score structure
  score.def_expected_but_not_detected = std::move(corrected_def_expected);
  score.num_def_wanted_not_found -= num_escape_peaks_removed;

  if( PeakFitImprove::debug_printout && (num_escape_peaks_removed > 0) )
  {
    cerr << "Removed " << num_escape_peaks_removed
         << " escape peak(s) from def_wanted_not_found count" << endl;
  }
}//correct_score_for_escape_peaks


CandidatePeakScore eval_candidate_settings( const FindCandidateSettings settings, const vector<DataSrcInfo> &input_srcs, const bool write_n42 )
{
  double sum_score = 0.0;
  size_t num_possibly_accepted_peaks_not_found = 0, num_def_wanted_not_found = 0;
  size_t num_extra_peaks = 0, num_peaks_found = 0;
  size_t num_def_wanted_peaks_found = 0;

  for( const DataSrcInfo &info : input_srcs )
  {
    shared_ptr<const SpecUtils::Measurement> src_spectrum = info.src_info.src_spectra.front();
    const vector<PeakDef> peaks = CandidatePeak_GA::find_candidate_peaks( src_spectrum, 0, 0, settings );

    // Use the refactored helper function to calculate score for this source
    const CandidatePeakScore source_score = calculate_candidate_peak_score_for_source( peaks, info.expected_photopeaks );

    // Accumulate scores and counts across all sources
    sum_score += source_score.score;
    num_peaks_found += source_score.num_peaks_found;
    num_possibly_accepted_peaks_not_found += source_score.num_possibly_accepted_peaks_not_found;
    num_def_wanted_not_found += source_score.num_def_wanted_not_found;
    num_def_wanted_peaks_found += source_score.num_def_wanted_peaks_found;
    num_extra_peaks += source_score.num_extra_peaks;

    if( PeakFitImprove::debug_printout )
    {
      cerr << "\n--- Source: " << info.src_info.src_name
           << " (" << info.detector_name << "/" << info.location_name << "/" << info.live_time_name << ")"
           << " ---\n"
           << "  Score: " << source_score.score
           << ", Found: " << source_score.num_peaks_found
           << ", Def wanted found: " << source_score.num_def_wanted_peaks_found
           << ", Def missed: " << source_score.num_def_wanted_not_found
           << ", Extra: " << source_score.num_extra_peaks
           << ", Total expected: " << info.expected_photopeaks.size()
           << ", Total detected: " << peaks.size()
           << endl;

      // Compare FWHM and area accuracy: zero-crossing vs Gaussian fit
      size_t zc_matched = 0, gf_matched = 0;
      double zc_sum_ratio = 0, zc_sum_abs_err = 0, zc_sum_sq_err = 0;
      double gf_sum_ratio = 0, gf_sum_abs_err = 0, gf_sum_sq_err = 0;
      size_t gf_fail_count = 0;

      // Area comparison accumulators
      double sd_area_sum_ratio = 0, gf_area_sum_ratio = 0;
      size_t sd_area_matched = 0, gf_area_matched = 0;

      for( const PeakDef &p : peaks )
      {
        // Find the best-matching expected photopeak
        const ExpectedPhotopeakInfo *best_match = nullptr;
        double best_dist = 1.0e9;
        for( const ExpectedPhotopeakInfo &epi : info.expected_photopeaks )
        {
          const double dist = std::fabs( p.mean() - epi.effective_energy );
          const double match_tol = std::max( 2.0 * p.sigma(), epi.effective_fwhm );
          if( dist < match_tol && dist < best_dist )
          {
            best_dist = dist;
            best_match = &epi;
          }
        }

        if( !best_match || best_match->effective_fwhm < 0.01 )
          continue;

        const double true_fwhm = best_match->effective_fwhm;
        const double true_area = best_match->peak_area;

        // Current zero-crossing FWHM
        const double zc_fwhm = p.fwhm();
        const double zc_ratio = zc_fwhm / true_fwhm;
        zc_sum_ratio += zc_ratio;
        zc_sum_abs_err += std::fabs( zc_fwhm - true_fwhm );
        zc_sum_sq_err += (zc_fwhm - true_fwhm) * (zc_fwhm - true_fwhm);
        ++zc_matched;

        // Second-derivative area (what PeakDef stores as "amplitude"/area)
        const double sd_area = p.amplitude();
        if( true_area > 1.0 && sd_area > 0.0 )
        {
          sd_area_sum_ratio += sd_area / true_area;
          ++sd_area_matched;
        }

        // Gaussian fit FWHM and area
        const double roi_lower = p.continuum()->lowerEnergy();
        const double roi_upper = p.continuum()->upperEnergy();
        double gf_area = -1.0;
        const double gf_sigma = estimate_sigma_from_raw_data( src_spectrum, p.mean(), roi_lower, roi_upper, &gf_area );

        if( gf_sigma > 0.0 )
        {
          const double gf_fwhm = gf_sigma * 2.35482;
          const double gf_ratio = gf_fwhm / true_fwhm;
          gf_sum_ratio += gf_ratio;
          gf_sum_abs_err += std::fabs( gf_fwhm - true_fwhm );
          gf_sum_sq_err += (gf_fwhm - true_fwhm) * (gf_fwhm - true_fwhm);
          ++gf_matched;

          if( true_area > 1.0 && gf_area > 0.0 )
          {
            gf_area_sum_ratio += gf_area / true_area;
            ++gf_area_matched;
          }

          // Check if this peak has a close neighbor in the expected photopeaks
          bool has_neighbor = false;
          for( const ExpectedPhotopeakInfo &epi : info.expected_photopeaks )
          {
            if( &epi == best_match )
              continue;
            if( std::fabs( epi.effective_energy - best_match->effective_energy ) < 3.0 * true_fwhm )
            {
              has_neighbor = true;
              break;
            }
          }

          // Per-peak output for bias investigation
          const double nsigma = best_match->nsigma_over_background;
          const double roi_width = roi_upper - roi_lower;
          const double roi_to_fwhm = roi_width / true_fwhm;
          cerr << "  PEAK_DIAG: E=" << p.mean()
               << " true_fwhm=" << true_fwhm
               << " gf_fwhm=" << gf_fwhm
               << " gf_ratio=" << gf_ratio
               << " sd_area=" << sd_area
               << " gf_area=" << gf_area
               << " true_area=" << true_area
               << " sd_area_ratio=" << (true_area > 1.0 ? sd_area / true_area : -1.0)
               << " gf_area_ratio=" << (true_area > 1.0 ? gf_area / true_area : -1.0)
               << " nsigma=" << nsigma
               << " neighbor=" << has_neighbor
               << endl;
        }
        else
        {
          ++gf_fail_count;
        }
      }

      if( zc_matched > 0 )
      {
        cerr << "  ZeroCross FWHM: " << zc_matched << " matched, "
             << "mean_ratio=" << (zc_sum_ratio / zc_matched)
             << ", mean_abs_err=" << (zc_sum_abs_err / zc_matched) << " keV"
             << ", rms_err=" << std::sqrt( zc_sum_sq_err / zc_matched ) << " keV" << endl;
      }
      if( gf_matched > 0 )
      {
        cerr << "  GaussFit  FWHM: " << gf_matched << " matched, "
             << "mean_ratio=" << (gf_sum_ratio / gf_matched)
             << ", mean_abs_err=" << (gf_sum_abs_err / gf_matched) << " keV"
             << ", rms_err=" << std::sqrt( gf_sum_sq_err / gf_matched ) << " keV";
        if( gf_fail_count > 0 )
          cerr << " (" << gf_fail_count << " fits failed)";
        cerr << endl;
      }

      // Area comparison summary
      if( sd_area_matched > 0 )
        cerr << "  SecDeriv AREA: " << sd_area_matched << " matched, mean_ratio=" << (sd_area_sum_ratio / sd_area_matched) << endl;
      if( gf_area_matched > 0 )
        cerr << "  GaussFit AREA: " << gf_area_matched << " matched, mean_ratio=" << (gf_area_sum_ratio / gf_area_matched) << endl;
    }

    // For N42 output, use the vectors from source_score (already computed by helper function)
    const vector<tuple<float,float,float>> &detected_expected = source_score.detected_expected;
    const vector<tuple<float,float,float>> &detected_not_expected = source_score.detected_not_expected;
    const vector<ExpectedPhotopeakInfo> &possibly_expected_but_not_detected = source_score.possibly_expected_but_not_detected;
    const vector<ExpectedPhotopeakInfo> &expected_and_was_detected = source_score.expected_and_was_detected;
    const vector<ExpectedPhotopeakInfo> &def_expected_but_not_detected = source_score.def_expected_but_not_detected;

    // For N42 output, we also need all detected peaks as tuples for visualization
    vector<tuple<float,float,float>> detected_peaks; //{mean, sigma, amplitude}
    for( const PeakDef &p : peaks )
      detected_peaks.emplace_back( p.mean(), p.sigma(), p.amplitude() );

    const vector<tuple<float,float,float>> orig_peak_candidates = detected_peaks;

    if( write_n42 )
    {
      //const string src_dir = SpecUtils::parent_path( info.src_info.file_base_path );

      string outdir = "output_n42";
      if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      outdir = SpecUtils::append_path( outdir, info.detector_name );
      if( !SpecUtils::is_directory(outdir) && SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      outdir = SpecUtils::append_path( outdir, info.location_name );
      if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      outdir = SpecUtils::append_path( outdir, info.live_time_name );
      if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      const string out_n42 = SpecUtils::append_path( outdir, info.src_info.src_name ) + "_peak_candidates.n42";

      SpecMeas output;

      output.add_remark( "Failed to find " + std::to_string(possibly_expected_but_not_detected.size() + def_expected_but_not_detected.size()) + " peaks that are possibly accepted." );
      output.add_remark( "Found " + std::to_string(detected_expected.size()) + " peaks that were expected." );
      output.add_remark( "Found " + std::to_string(detected_not_expected.size()) + " peaks that were NOT expected." );
      output.add_remark( settings.print("settings") );
      output.set_instrument_model( info.detector_name );
      if( SpecUtils::icontains(info.detector_name, "Detective" ) )
        output.set_manufacturer( "ORTEC" );
      else if( SpecUtils::icontains(info.detector_name, "Falcon 5000" ) )
        output.set_manufacturer( "Canberra" );
      else if( SpecUtils::icontains(info.detector_name, "Fulcrum" ) )
        output.set_manufacturer( "PHDS" );

      output.set_measurement_location_name( info.location_name );

      shared_ptr<SpecUtils::Measurement> out_with_cand = make_shared<SpecUtils::Measurement>( *src_spectrum );
      out_with_cand->set_sample_number( 1 );

      const string title = info.detector_name + "/" + info.src_info.src_name;
      out_with_cand->set_title( title + " + candidate peaks" );
      auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
      out_with_cand->set_start_time( now );
      output.add_measurement( out_with_cand, false );

      deque<shared_ptr<const PeakDef>> peaks;

      vector<tuple<float,float,float>> all_candidates = detected_expected;
      all_candidates.insert( end(all_candidates), begin(detected_not_expected), end(detected_not_expected) );
      for( size_t peak_index = 0; peak_index < all_candidates.size(); ++peak_index )
      {
        const tuple<float,float,float> &p = all_candidates[peak_index]; //{mean, sigma, amplitude}
        const float mean = get<0>(p);
        const float sigma = get<1>(p);
        const float amp = get<2>(p);

        const bool is_expected = (peak_index < detected_expected.size());

        auto peak = make_shared<PeakDef>( mean, sigma, amp );
        peak->setFitFor( PeakDef::CoefficientType::Mean, false );
        peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
        peak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
        peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
        peak->continuum()->setRange( mean - 3*sigma, mean + 3*sigma );
        peak->continuum()->calc_linear_continuum_eqn( src_spectrum, mean, mean - 4*sigma, mean + 4*sigma, 5, 5 );
        peak->setLineColor( is_expected ? Wt::GlobalColor::darkGreen : Wt::GlobalColor::red );

        if( is_expected )
        {
          //Put truth-level info to user-label field
          string info_str;
          for( const ExpectedPhotopeakInfo &exp_info : info.expected_photopeaks )
          {
            if( mean > exp_info.roi_lower && mean < exp_info.roi_upper )
            {
              if( !info_str.empty() )
                info_str += ".\n";

              info_str += "E: A=" + SpecUtils::printCompact(exp_info.peak_area, 3)
              + ", W=" + SpecUtils::printCompact(exp_info.gamma_lines.front().fwhm, 3)
              + ", #s=" + SpecUtils::printCompact(exp_info.nsigma_over_background, 3);
            }
          }

          peak->setUserLabel( info_str );
        }else
        {
          peak->setUserLabel( "Not Expected" );
        }

        peaks.push_back( peak );
      }//for( const tuple<float,float,float> &p : all_candidates )


      vector<ExpectedPhotopeakInfo> missing_peaks = def_expected_but_not_detected;
      missing_peaks.insert( end(missing_peaks), begin(possibly_expected_but_not_detected), end(possibly_expected_but_not_detected) );

      for( size_t missing_index = 0; missing_index < missing_peaks.size(); ++missing_index )
      {
        const ExpectedPhotopeakInfo &missing = missing_peaks[missing_index];
        const bool def_wanted = (missing_index < def_expected_but_not_detected.size());

        const double mean = missing.effective_energy;
        const double sigma = missing.gamma_lines.front().fwhm/2.35482f;
        auto peak = make_shared<PeakDef>( mean, sigma, missing.peak_area );
        peak->setFitFor( PeakDef::CoefficientType::Mean, false );
        peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
        peak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
        peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
        peak->continuum()->setRange( mean - 3*sigma, mean + 3*sigma );
        peak->continuum()->calc_linear_continuum_eqn( src_spectrum, mean, mean - 4*sigma, mean + 4*sigma, 5, 5 );
        const Wt::WColor orange(252, 94, 3), yellow(168, 150, 50);
        peak->setLineColor( def_wanted ? orange : yellow );
        string label = def_wanted ? "Not det. - wanted" : "Not det. - poss.";
        label += ", N_Sigma=" + SpecUtils::printCompact(missing.nsigma_over_background, 3);
        peak->setUserLabel( label );

        peaks.push_back( peak );
      }//for( const ExpectedPhotopeakInfo &missing : not_detected )





      output.setPeaks( peaks, {1} );

      {// begin add second derivative to N42
        auto second_deriv = make_shared<vector<float>>();
        const int side_bins = settings.num_smooth_side_channels;
        const int poly_order = settings.smooth_polynomial_order;
        smoothSpectrum( src_spectrum, side_bins, poly_order, 2, *second_deriv );

        shared_ptr<SpecUtils::Measurement> second_deriv_meas = make_shared<SpecUtils::Measurement>();
        second_deriv_meas->set_gamma_counts( second_deriv, 1.0, 1.0 );
        second_deriv_meas->set_energy_calibration( src_spectrum->energy_calibration() );
        second_deriv_meas->set_sample_number( 2 );
        second_deriv_meas->set_title( title + " + smoothed second derivative." );
        output.add_measurement( second_deriv_meas, true );

        deque<shared_ptr<const PeakDef>> candidates;
        for( const tuple<float,float,float> &p : orig_peak_candidates )
        {
          auto peak = make_shared<PeakDef>( get<0>(p), get<1>(p), get<2>(p) );
          peak->continuum()->setType( PeakContinuum::OffsetType::Constant );
          peak->continuum()->setPolynomialCoef(0, 0.0);
          candidates.push_back( peak );
        }

        output.setPeaks( candidates, {2} );
      }// end add second derivative to N42


      {
        // I think this will be good for estimating sigma and area
        const int side_bins = 3;
        const int poly_order = 3;
        vector<float> rougher_second_deriv;
        smoothSpectrum( src_spectrum, side_bins, poly_order, 2, rougher_second_deriv );
        shared_ptr<SpecUtils::Measurement> rough_second_deriv_meas = make_shared<SpecUtils::Measurement>();
        auto rough_second_deriv_counts = make_shared<vector<float>>( rougher_second_deriv );
        rough_second_deriv_meas->set_gamma_counts( rough_second_deriv_counts, 1.0, 1.0 );
        rough_second_deriv_meas->set_energy_calibration( src_spectrum->energy_calibration() );
        rough_second_deriv_meas->set_sample_number( 3 );
        rough_second_deriv_meas->set_title( title + " + rougher smoothed second derivative." );
        output.add_measurement( rough_second_deriv_meas, true );
      }


      {
        const int side_bins = settings.num_smooth_side_channels;
        const int poly_order = 1;
        vector<float> smoothed_data;
        smoothSpectrum( src_spectrum, side_bins, poly_order, 0, smoothed_data );
        shared_ptr<SpecUtils::Measurement> meas = make_shared<SpecUtils::Measurement>();
        auto rough_second_deriv_counts = make_shared<vector<float>>( smoothed_data );
        meas->set_gamma_counts( rough_second_deriv_counts, out_with_cand->live_time(), out_with_cand->real_time() );
        meas->set_energy_calibration( src_spectrum->energy_calibration() );
        meas->set_sample_number( 4 );
        meas->set_title( title + " + smoothed data." );
        output.add_measurement( meas, true );
      }

      //ofstream outstrm( out_n42.c_str(), ios::out | ios::binary );
      //if( !outstrm )
      //  cerr << "Failed to open '" << out_n42 << "'" << endl;
      //output.write( outstrm, output.sample_numbers(), output.detector_names(), SpecUtils::SaveSpectrumAsType::N42_2012 );

      output.save2012N42File( out_n42, [=](){
        cerr << "Failed to write '" << out_n42 << "'" << endl;
      });

      //cout << "Wrote '" << out_n42 << endl;
    }//if( write_n42 )
  }//for( const DataSrcInfo &info : input_srcs )

  //cout << "Avrg score: " << sum_score/input_srcs.size() << endl;

  CandidatePeakScore result;
  result.score = sum_score / input_srcs.size();
  result.num_peaks_found = num_peaks_found;
  result.num_def_wanted_not_found = num_def_wanted_not_found;
  result.num_def_wanted_peaks_found = num_def_wanted_peaks_found;
  result.num_possibly_accepted_peaks_not_found = num_possibly_accepted_peaks_not_found;
  result.num_extra_peaks = num_extra_peaks;
  return result;
};//eval_candidate_settings lambda


string CandidatePeakSolution::to_string( const string &separator ) const
{
  return
  string("num_smooth_side_channels: ") + std::to_string(num_smooth_side_channels)
  + separator + "smooth_polynomial_order: " + std::to_string(smooth_polynomial_order)
  + separator + "threshold_FOM: " + std::to_string(threshold_FOM)
  + separator + "more_scrutiny_FOM_threshold_delta: " + std::to_string(more_scrutiny_FOM_threshold_delta)
  + separator + "pos_sum_threshold_sf: " + std::to_string(pos_sum_threshold_sf)
  + separator + "num_chan_fluctuate: " + std::to_string(num_chan_fluctuate)
  + separator + "more_scrutiny_coarser_FOM_delta: " + std::to_string(more_scrutiny_coarser_FOM_delta)
  + separator + "more_scrutiny_min_dev_from_line: " + std::to_string(more_scrutiny_min_dev_from_line)
  + separator + "amp_to_apply_line_test_below: " + std::to_string(amp_to_apply_line_test_below)
  + separator + "smooth_ref_fraction: " + std::to_string(smooth_ref_fraction)
  + separator + "smooth_scale_power: " + std::to_string(smooth_scale_power)
  + separator + "min_second_deriv_significance: " + std::to_string(min_second_deriv_significance)
  + separator + "compton_next_ratio_max: " + std::to_string(compton_next_ratio_max)
  + separator + "compton_prev_ratio_min: " + std::to_string(compton_prev_ratio_min)
  + separator + "compton_total_ratio_min: " + std::to_string(compton_total_ratio_min)
  + separator + "low_energy_test_max_keV: " + std::to_string(low_energy_test_max_keV)
  + separator + "low_energy_drop_fraction: " + std::to_string(low_energy_drop_fraction)
  + separator + "pcgap_feature_nsigma: " + std::to_string(pcgap_feature_nsigma)
  + separator + "pcgap_max_extent_nsigma: " + std::to_string(pcgap_max_extent_nsigma)
  + separator + "pcgap_roi_blend_weight: " + std::to_string(pcgap_roi_blend_weight)
  + separator + "pcgap_fom_blend_threshold: " + std::to_string(pcgap_fom_blend_threshold)
  ;
}


void init_genes(CandidatePeakSolution& p,const std::function<double(void)> &rnd01)
{
  // rnd01() gives a random number in 0~1
  p.num_smooth_side_channels          = 3+10*rnd01();
  p.smooth_polynomial_order           = 2;
  p.threshold_FOM                     = 0.75+1.25*rnd01();
  p.more_scrutiny_FOM_threshold_delta = -0.1+2.1*rnd01();
  p.pos_sum_threshold_sf              = -0.15+0.3*rnd01();
  p.num_chan_fluctuate                = 1;
  p.more_scrutiny_coarser_FOM_delta   = -0.1+4.1*rnd01();
  p.more_scrutiny_min_dev_from_line   = 0 + 7*rnd01();
  p.amp_to_apply_line_test_below      = 0.0+100*rnd01();

  // Energy-adaptive smoothing
  p.smooth_ref_fraction               = 0.05+0.20*rnd01();  // [0.05, 0.25]
  p.smooth_scale_power                = 0.6*rnd01();         // [0, 0.6]; 0.5 = sqrt (NaI physics)
  p.min_second_deriv_significance     = 8.0*rnd01();         // [0, 8]; 0 = disabled

  // Compton backscatter test
  p.compton_next_ratio_max            = 2.0+18.0*rnd01();   // [2, 20]; old algo: 4.0
  p.compton_prev_ratio_min            = 0.5*rnd01();         // [0, 0.5]; old algo: 0.2
  p.compton_total_ratio_min           = 0.8*rnd01();         // [0, 0.8]; old algo: 0.3

  // Low-energy drop-off test
  p.low_energy_test_max_keV           = 200*rnd01();         // [0, 200]; old algo: 130
  p.low_energy_drop_fraction          = 0.1+0.7*rnd01();    // [0.1, 0.8]; old algo: 0.45

  // PCGAP ROI
  p.pcgap_feature_nsigma              = 1.5+3.5*rnd01();    // [1.5, 5.0]; old algo: 3.0
  p.pcgap_max_extent_nsigma           = 2.5+2.5*rnd01();    // [2.5, 5.0]
  p.pcgap_roi_blend_weight            = rnd01();             // [0, 1]; 0 = disabled
  p.pcgap_fom_blend_threshold         = 1.0+9.0*rnd01();    // [1, 10]
}

FindCandidateSettings genes_to_settings( const CandidatePeakSolution &p )
{
  FindCandidateSettings settings;

  settings.num_smooth_side_channels = p.num_smooth_side_channels;
  settings.smooth_polynomial_order = p.smooth_polynomial_order;
  settings.threshold_FOM = p.threshold_FOM;
  settings.more_scrutiny_FOM_threshold = p.threshold_FOM + p.more_scrutiny_FOM_threshold_delta;
  settings.pos_sum_threshold_sf = p.pos_sum_threshold_sf;
  settings.num_chan_fluctuate = p.num_chan_fluctuate;
  settings.more_scrutiny_coarser_FOM = p.threshold_FOM + p.more_scrutiny_coarser_FOM_delta;
  settings.more_scrutiny_min_dev_from_line = p.more_scrutiny_min_dev_from_line;
  settings.amp_to_apply_line_test_below = p.amp_to_apply_line_test_below;

  // Energy-adaptive smoothing
  settings.smooth_ref_fraction = static_cast<float>( p.smooth_ref_fraction );
  settings.smooth_scale_power = static_cast<float>( p.smooth_scale_power );
  settings.min_second_deriv_significance = static_cast<float>( p.min_second_deriv_significance );

  // Compton backscatter test
  settings.compton_next_ratio_max = static_cast<float>( p.compton_next_ratio_max );
  settings.compton_prev_ratio_min = static_cast<float>( p.compton_prev_ratio_min );
  settings.compton_total_ratio_min = static_cast<float>( p.compton_total_ratio_min );

  // Low-energy drop-off test
  settings.low_energy_test_max_keV = static_cast<float>( p.low_energy_test_max_keV );
  settings.low_energy_drop_fraction = static_cast<float>( p.low_energy_drop_fraction );

  // PCGAP ROI
  settings.pcgap_feature_nsigma = static_cast<float>( p.pcgap_feature_nsigma );
  settings.pcgap_max_extent_nsigma = static_cast<float>( p.pcgap_max_extent_nsigma );
  settings.pcgap_roi_blend_weight = static_cast<float>( p.pcgap_roi_blend_weight );
  settings.pcgap_fom_blend_threshold = static_cast<float>( p.pcgap_fom_blend_threshold );

  return settings;
}

bool eval_solution( const CandidatePeakSolution &p, CandidatePeakCost &c )
{
  const FindCandidateSettings settings = genes_to_settings( p );

  assert(ns_ga_eval_fcn);

  try
  {
    c.objective1 = -1.0 * ns_ga_eval_fcn( settings );
  }catch( std::exception & )
  {
    return false; //reject solution
  }

  return true; // solution is accepted
}


CandidatePeakSolution mutate(
    const CandidatePeakSolution& X_base,
    const std::function<double(void)> &rnd01,
    double shrink_scale)
{
  CandidatePeakSolution X_new;
  const double mu = 0.2*shrink_scale; // mutation radius (adjustable)

  const double mutate_threshold = PeakFitImprove::sm_ga_mutate_threshold;

  bool in_range;

  size_t num_tries = 0;

  do{
    num_tries += 1;
    if( num_tries > 1000 )
    {
      cerr << "Has taken over " << num_tries << " tries to find some genes.\nX_new={\n"
      << X_new.to_string( "\n\t" )
      << "}\nX_base={\n"
      << X_base.to_string("\n\t" )
      << "}\nWill return new randomly inited gene.\n"
      << endl;

      std::random_device rng;
      std::uniform_real_distribution<double> unif_dist( 0.0, 1.0 );

      init_genes( X_new, [&](){return unif_dist(rng);} );
      return X_new;
    }

    in_range = true;
    X_new = X_base;

    if( rnd01() > mutate_threshold )
      X_new.num_smooth_side_channels += shrink_scale*(rnd01()-rnd01()); //not multiplying by `mu`, because we can get stuck in a single int
    if( rnd01() > mutate_threshold )
      X_new.num_smooth_side_channels = std::max( X_new.num_smooth_side_channels, 3 );
    if( rnd01() > mutate_threshold )
    {
      X_new.num_smooth_side_channels = std::min( X_new.num_smooth_side_channels, 13 );
      in_range=in_range&&(X_new.num_smooth_side_channels>=3 && X_new.num_smooth_side_channels<=13);
    }

    //X_new.smooth_polynomial_order+=mu*(rnd01()-rnd01()); //This is an int we could get stuck in...
    //in_range=in_range&&(X_new.smooth_polynomial_order>=2 && X_new.smooth_polynomial_order<3);

    if( rnd01() > mutate_threshold )
    {
      X_new.threshold_FOM += mu * (rnd01() - rnd01());
      //in_range=in_range&&(X_new.threshold_FOM>=0.75 && X_new.threshold_FOM<3.5);
      in_range=in_range&&(X_new.threshold_FOM>=0.75 && X_new.threshold_FOM<=2);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.more_scrutiny_FOM_threshold_delta += mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_FOM_threshold_delta>=-0.5 && X_new.more_scrutiny_FOM_threshold_delta<3.5);
      in_range=in_range&&(X_new.more_scrutiny_FOM_threshold_delta>=-0.1 && X_new.more_scrutiny_FOM_threshold_delta<=2);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.pos_sum_threshold_sf += 0.05 * mu * (rnd01() - rnd01()); //mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.pos_sum_threshold_sf>=-0.1 && X_new.pos_sum_threshold_sf<0.1);
      in_range=in_range&&(X_new.pos_sum_threshold_sf>=-0.15 && X_new.pos_sum_threshold_sf<=0.15);
    }

    //X_new.num_chan_fluctuate+=mu*(rnd01()-rnd01());  //THis is an int we could get stuck in...
    //in_range=in_range&&(X_new.num_chan_fluctuate>=1 && X_new.num_chan_fluctuate<4);

    if( rnd01() > mutate_threshold )
    {
      X_new.more_scrutiny_coarser_FOM_delta += mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_coarser_FOM_delta>=-0.1 && X_new.more_scrutiny_coarser_FOM_delta<5);
      in_range=in_range&&(X_new.more_scrutiny_coarser_FOM_delta>=-0.1 && X_new.more_scrutiny_coarser_FOM_delta<4);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.more_scrutiny_min_dev_from_line += 2 * mu*(rnd01()-rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_min_dev_from_line>=0.0 && X_new.more_scrutiny_min_dev_from_line<10.0);
      in_range=in_range&&(X_new.more_scrutiny_min_dev_from_line>=0.0 && X_new.more_scrutiny_min_dev_from_line<7.0);
    }

    if( rnd01() > mutate_threshold )
      X_new.amp_to_apply_line_test_below += 5 * mu*(rnd01()-rnd01());
    if( rnd01() > mutate_threshold )
      X_new.amp_to_apply_line_test_below = std::max( X_new.amp_to_apply_line_test_below, 0 );
    if( rnd01() > mutate_threshold )
    {
      X_new.amp_to_apply_line_test_below = std::min( X_new.amp_to_apply_line_test_below, 100 );
      in_range=in_range&&(X_new.amp_to_apply_line_test_below>=0.0 && X_new.amp_to_apply_line_test_below<=100);
      //in_range=in_range&&(X_new.amp_to_apply_line_test_below>=0.0 && X_new.amp_to_apply_line_test_below<=75);
    }

    // Energy-adaptive smoothing
    if( rnd01() > mutate_threshold )
    {
      X_new.smooth_ref_fraction += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.smooth_ref_fraction >= 0.05 && X_new.smooth_ref_fraction <= 0.25);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.smooth_scale_power += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.smooth_scale_power >= 0.0 && X_new.smooth_scale_power <= 0.6);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.min_second_deriv_significance += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.min_second_deriv_significance >= 0.0 && X_new.min_second_deriv_significance <= 8.0);
    }

    // Compton backscatter test
    if( rnd01() > mutate_threshold )
    {
      X_new.compton_next_ratio_max += 2.0 * mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.compton_next_ratio_max >= 2.0 && X_new.compton_next_ratio_max <= 20.0);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.compton_prev_ratio_min += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.compton_prev_ratio_min >= 0.0 && X_new.compton_prev_ratio_min <= 0.5);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.compton_total_ratio_min += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.compton_total_ratio_min >= 0.0 && X_new.compton_total_ratio_min <= 0.8);
    }

    // Low-energy drop-off test
    if( rnd01() > mutate_threshold )
    {
      X_new.low_energy_test_max_keV += 10.0 * mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.low_energy_test_max_keV >= 0.0 && X_new.low_energy_test_max_keV <= 200.0);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.low_energy_drop_fraction += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.low_energy_drop_fraction >= 0.1 && X_new.low_energy_drop_fraction <= 0.8);
    }

    // PCGAP ROI
    if( rnd01() > mutate_threshold )
    {
      X_new.pcgap_feature_nsigma += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.pcgap_feature_nsigma >= 1.5 && X_new.pcgap_feature_nsigma <= 5.0);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.pcgap_max_extent_nsigma += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.pcgap_max_extent_nsigma >= 2.5 && X_new.pcgap_max_extent_nsigma <= 5.0);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.pcgap_roi_blend_weight += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.pcgap_roi_blend_weight >= 0.0 && X_new.pcgap_roi_blend_weight <= 1.0);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.pcgap_fom_blend_threshold += mu * (rnd01() - rnd01());
      in_range = in_range && (X_new.pcgap_fom_blend_threshold >= 1.0 && X_new.pcgap_fom_blend_threshold <= 10.0);
    }

  } while(!in_range);
  return X_new;
}


CandidatePeakSolution crossover(
    const CandidatePeakSolution& X1,
    const CandidatePeakSolution& X2,
    const std::function<double(void)> &rnd01)
{
  const double crossover_threshold = PeakFitImprove::sm_ga_crossover_threshold;


  CandidatePeakSolution X_new = X1;


  double r;
  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.num_smooth_side_channels=r*X1.num_smooth_side_channels+(1.0-r)*X2.num_smooth_side_channels;
    //X_new.num_smooth_side_channels = X1.num_smooth_side_channels;
  }else if( rnd01() < 0.5 )
  {
    X_new.num_smooth_side_channels = X2.num_smooth_side_channels;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    //X_new.smooth_polynomial_order=r*X1.smooth_polynomial_order+(1.0-r)*X2.smooth_polynomial_order;
    X_new.smooth_polynomial_order = X1.smooth_polynomial_order;
  }else if( rnd01() < 0.5 )
  {
    X_new.smooth_polynomial_order = X2.smooth_polynomial_order;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.threshold_FOM=r*X1.threshold_FOM+(1.0-r)*X2.threshold_FOM;
  }else if( rnd01() < 0.5 )
  {
    X_new.threshold_FOM = X2.threshold_FOM;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.more_scrutiny_FOM_threshold_delta=r*X1.more_scrutiny_FOM_threshold_delta+(1.0-r)*X2.more_scrutiny_FOM_threshold_delta;
  }else if( rnd01() < 0.5 )
  {
    X_new.more_scrutiny_FOM_threshold_delta = X2.more_scrutiny_FOM_threshold_delta;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.pos_sum_threshold_sf=r*X1.pos_sum_threshold_sf+(1.0-r)*X2.pos_sum_threshold_sf;
  }else if( rnd01() < 0.5 )
  {
    X_new.pos_sum_threshold_sf = X2.pos_sum_threshold_sf;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    //X_new.num_chan_fluctuate=r*X1.num_chan_fluctuate+(1.0-r)*X2.num_chan_fluctuate;
    X_new.num_chan_fluctuate = X1.num_chan_fluctuate;
  }else if( rnd01() < 0.5 )
  {
    X_new.num_chan_fluctuate = X2.num_chan_fluctuate;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.more_scrutiny_coarser_FOM_delta=r*X1.more_scrutiny_coarser_FOM_delta+(1.0-r)*X2.more_scrutiny_coarser_FOM_delta;
  }else if( rnd01() < 0.5 )
  {
    X_new.more_scrutiny_coarser_FOM_delta = X2.more_scrutiny_coarser_FOM_delta;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.more_scrutiny_min_dev_from_line=r*X1.more_scrutiny_min_dev_from_line+(1.0-r)*X2.more_scrutiny_min_dev_from_line;
  }else if( rnd01() < 0.5 )
  {
    X_new.more_scrutiny_min_dev_from_line = X2.more_scrutiny_min_dev_from_line;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.amp_to_apply_line_test_below=r*X1.amp_to_apply_line_test_below+(1.0-r)*X2.amp_to_apply_line_test_below;
  }else if( rnd01() < 0.5 )
  {
    X_new.amp_to_apply_line_test_below = X2.amp_to_apply_line_test_below;
  }

  // Energy-adaptive smoothing
  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.smooth_ref_fraction = r * X1.smooth_ref_fraction + (1.0 - r) * X2.smooth_ref_fraction;
  }else if( rnd01() < 0.5 )
  {
    X_new.smooth_ref_fraction = X2.smooth_ref_fraction;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.smooth_scale_power = r * X1.smooth_scale_power + (1.0 - r) * X2.smooth_scale_power;
  }else if( rnd01() < 0.5 )
  {
    X_new.smooth_scale_power = X2.smooth_scale_power;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.min_second_deriv_significance = r * X1.min_second_deriv_significance + (1.0 - r) * X2.min_second_deriv_significance;
  }else if( rnd01() < 0.5 )
  {
    X_new.min_second_deriv_significance = X2.min_second_deriv_significance;
  }

  // Compton backscatter test
  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.compton_next_ratio_max = r * X1.compton_next_ratio_max + (1.0 - r) * X2.compton_next_ratio_max;
  }else if( rnd01() < 0.5 )
  {
    X_new.compton_next_ratio_max = X2.compton_next_ratio_max;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.compton_prev_ratio_min = r * X1.compton_prev_ratio_min + (1.0 - r) * X2.compton_prev_ratio_min;
  }else if( rnd01() < 0.5 )
  {
    X_new.compton_prev_ratio_min = X2.compton_prev_ratio_min;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.compton_total_ratio_min = r * X1.compton_total_ratio_min + (1.0 - r) * X2.compton_total_ratio_min;
  }else if( rnd01() < 0.5 )
  {
    X_new.compton_total_ratio_min = X2.compton_total_ratio_min;
  }

  // Low-energy drop-off test
  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.low_energy_test_max_keV = r * X1.low_energy_test_max_keV + (1.0 - r) * X2.low_energy_test_max_keV;
  }else if( rnd01() < 0.5 )
  {
    X_new.low_energy_test_max_keV = X2.low_energy_test_max_keV;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.low_energy_drop_fraction = r * X1.low_energy_drop_fraction + (1.0 - r) * X2.low_energy_drop_fraction;
  }else if( rnd01() < 0.5 )
  {
    X_new.low_energy_drop_fraction = X2.low_energy_drop_fraction;
  }

  // PCGAP ROI
  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.pcgap_feature_nsigma = r * X1.pcgap_feature_nsigma + (1.0 - r) * X2.pcgap_feature_nsigma;
  }else if( rnd01() < 0.5 )
  {
    X_new.pcgap_feature_nsigma = X2.pcgap_feature_nsigma;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.pcgap_max_extent_nsigma = r * X1.pcgap_max_extent_nsigma + (1.0 - r) * X2.pcgap_max_extent_nsigma;
  }else if( rnd01() < 0.5 )
  {
    X_new.pcgap_max_extent_nsigma = X2.pcgap_max_extent_nsigma;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.pcgap_roi_blend_weight = r * X1.pcgap_roi_blend_weight + (1.0 - r) * X2.pcgap_roi_blend_weight;
  }else if( rnd01() < 0.5 )
  {
    X_new.pcgap_roi_blend_weight = X2.pcgap_roi_blend_weight;
  }

  if( rnd01() < crossover_threshold )
  {
    r = rnd01();
    X_new.pcgap_fom_blend_threshold = r * X1.pcgap_fom_blend_threshold + (1.0 - r) * X2.pcgap_fom_blend_threshold;
  }else if( rnd01() < 0.5 )
  {
    X_new.pcgap_fom_blend_threshold = X2.pcgap_fom_blend_threshold;
  }

  return X_new;
}


double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
{
  // finalize the cost
  double final_cost=0.0;
  final_cost+=X.middle_costs.objective1;
  return final_cost;
}


void SO_report_generation(
    int generation_number,
    const EA::GenerationType<CandidatePeakSolution,CandidatePeakCost> &last_generation,
    const CandidatePeakSolution& best_genes)
{
  bool best_yet = false;
  if( !sm_set_best_genes || (last_generation.best_total_cost < sm_best_total_cost) )
  {
    best_yet = true;
    sm_set_best_genes = true;
    sm_best_genes = best_genes;
    sm_best_total_cost = last_generation.best_total_cost;
  }

  cout << "Generation [" << generation_number << "], "
       << "Best=" << last_generation.best_total_cost << ", "
       << "Average=" << last_generation.average_cost << ", "
       << "Exe_time=" << last_generation.exe_time << "\n";

  // Print score details for the best individual
  if( ns_input_srcs )
  {
    const FindCandidateSettings settings = genes_to_settings( best_genes );
    const CandidatePeakScore score = eval_candidate_settings( settings, *ns_input_srcs, false );
    cout << "  Score=" << score.score
         << ", Found=" << score.num_peaks_found
         << ", DefWantedFound=" << score.num_def_wanted_peaks_found
         << ", DefWantedMissed=" << score.num_def_wanted_not_found
         << ", MaybeWantedMissed=" << score.num_possibly_accepted_peaks_not_found
         << ", Extra=" << score.num_extra_peaks << "\n";
  }

  cout << "Best genes: {\n\t" << best_genes.to_string( "\n\t" ) << "\n}\n" << endl;

  sm_output_file
  << generation_number << "\t"
  << last_generation.average_cost << "\t"
  << last_generation.best_total_cost << "\t"
  << "{" << best_genes.to_string(", ") << "}\n\n";
}


FindCandidateSettings do_ga_eval( std::function<double(const FindCandidateSettings &)> ga_eval_fcn,
                                  const std::vector<DataSrcInfo> *input_srcs_for_report )
{
  assert( !sm_has_been_called );
  if( sm_has_been_called )
  {
    cerr << "You should only call CandidatePeak_GA::do_ga_eval(...) once per program execution!!!" << endl;
    exit(1);
  }

  sm_has_been_called = true;

  assert( !!ga_eval_fcn );
  if( !ga_eval_fcn )
    throw runtime_error( "Invalid eval function passed in." );

  ns_ga_eval_fcn = ga_eval_fcn;
  ns_input_srcs = input_srcs_for_report;


  sm_output_file.open("results.txt");
  sm_output_file<<"step"<<"\t"<<"cost_avg"<<"\t"<<"cost_best"<<"\t"<<"solution_best"<<"\n";

  EA::Chronometer timer;
  timer.tic();

  GA_Type ga_obj;
  ga_obj.problem_mode=EA::GA_MODE::SOGA; //Single objective genetic algorithm
  ga_obj.multi_threading=true;
  ga_obj.idle_delay_us=1; // switch between threads quickly
  ga_obj.dynamic_threading = true; //If false,  thread responsibilities are divided at the beginning
  ga_obj.verbose=false;
  ga_obj.population = static_cast<unsigned int>(PeakFitImprove::sm_ga_population);
  ga_obj.generation_max = static_cast<int>(PeakFitImprove::sm_ga_generation_max);
  ga_obj.calculate_SO_total_fitness=calculate_SO_total_fitness;
  ga_obj.init_genes=init_genes;
  ga_obj.eval_solution=eval_solution;
  ga_obj.mutate=mutate;
  ga_obj.crossover=crossover;
  ga_obj.SO_report_generation=SO_report_generation;
  ga_obj.crossover_fraction=PeakFitImprove::sm_ga_crossover_fraction;
  ga_obj.mutation_rate=PeakFitImprove::sm_ga_mutation_rate;
  ga_obj.best_stall_max = static_cast<int>(PeakFitImprove::sm_ga_best_stall_max);
  ga_obj.elite_count = static_cast<int>(PeakFitImprove::sm_ga_elite_count);
  ga_obj.N_threads = static_cast<int>(PeakFitImprove::sm_num_optimization_threads);
  EA::StopReason stop_reason = ga_obj.solve();

  cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason) << endl;
  cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;

  sm_output_file.close();

  return genes_to_settings( sm_best_genes );
}
}//namespace CandidatePeak_GA
