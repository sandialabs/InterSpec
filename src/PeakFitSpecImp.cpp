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
#include <limits>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>

#include <boost/math/constants/constants.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitSpecImp.h"

using namespace std;


// ============================================================================
// FindCandidateSettings print/to_json
// ============================================================================

namespace PeakFitSpec
{

string FindCandidateSettings::print( const string &var_name ) const
{
  string result;
  result += var_name + ".num_smooth_side_channels = "   + to_string( num_smooth_side_channels )       + ";\n";
  result += var_name + ".smooth_polynomial_order = "    + to_string( smooth_polynomial_order )        + ";\n";
  result += var_name + ".threshold_FOM = "              + to_string( threshold_FOM )                  + ";\n";
  result += var_name + ".more_scrutiny_FOM_threshold = "+ to_string( more_scrutiny_FOM_threshold )    + ";\n";
  result += var_name + ".pos_sum_threshold_sf = "       + to_string( pos_sum_threshold_sf )           + ";\n";
  result += var_name + ".num_chan_fluctuate = "          + to_string( num_chan_fluctuate )             + ";\n";
  result += var_name + ".more_scrutiny_coarser_FOM = "  + to_string( more_scrutiny_coarser_FOM )      + ";\n";
  result += var_name + ".more_scrutiny_min_dev_from_line = " + to_string( more_scrutiny_min_dev_from_line ) + ";\n";
  result += var_name + ".amp_to_apply_line_test_below = "    + to_string( amp_to_apply_line_test_below )    + ";\n";
  result += var_name + ".smooth_ref_fraction = "             + to_string( smooth_ref_fraction )              + ";\n";
  result += var_name + ".smooth_scale_power = "             + to_string( smooth_scale_power )               + ";\n";
  result += var_name + ".min_second_deriv_significance = "  + to_string( min_second_deriv_significance )    + ";\n";
  result += var_name + ".compton_next_ratio_max = "         + to_string( compton_next_ratio_max )           + ";\n";
  result += var_name + ".compton_prev_ratio_min = "         + to_string( compton_prev_ratio_min )           + ";\n";
  result += var_name + ".compton_total_ratio_min = "        + to_string( compton_total_ratio_min )          + ";\n";
  result += var_name + ".low_energy_test_max_keV = "        + to_string( low_energy_test_max_keV )          + ";\n";
  result += var_name + ".low_energy_drop_fraction = "       + to_string( low_energy_drop_fraction )         + ";\n";
  result += var_name + ".pcgap_feature_nsigma = "           + to_string( pcgap_feature_nsigma )             + ";\n";
  result += var_name + ".pcgap_max_extent_nsigma = "        + to_string( pcgap_max_extent_nsigma )          + ";\n";
  result += var_name + ".pcgap_roi_blend_weight = "         + to_string( pcgap_roi_blend_weight )           + ";\n";
  result += var_name + ".pcgap_fom_blend_threshold = "      + to_string( pcgap_fom_blend_threshold )        + ";\n";
  return result;
}//FindCandidateSettings::print


string FindCandidateSettings::to_json() const
{
  string result = "{\n";
  result += "  \"num_smooth_side_channels\": "   + to_string( num_smooth_side_channels )       + ",\n";
  result += "  \"smooth_polynomial_order\": "    + to_string( smooth_polynomial_order )        + ",\n";
  result += "  \"threshold_FOM\": "              + to_string( threshold_FOM )                  + ",\n";
  result += "  \"more_scrutiny_FOM_threshold\": "+ to_string( more_scrutiny_FOM_threshold )    + ",\n";
  result += "  \"pos_sum_threshold_sf\": "       + to_string( pos_sum_threshold_sf )           + ",\n";
  result += "  \"num_chan_fluctuate\": "          + to_string( num_chan_fluctuate )             + ",\n";
  result += "  \"more_scrutiny_coarser_FOM\": "  + to_string( more_scrutiny_coarser_FOM )      + ",\n";
  result += "  \"more_scrutiny_min_dev_from_line\": " + to_string( more_scrutiny_min_dev_from_line ) + ",\n";
  result += "  \"amp_to_apply_line_test_below\": "    + to_string( amp_to_apply_line_test_below )    + ",\n";
  result += "  \"smooth_ref_fraction\": "             + to_string( smooth_ref_fraction )              + ",\n";
  result += "  \"smooth_scale_power\": "             + to_string( smooth_scale_power )               + ",\n";
  result += "  \"min_second_deriv_significance\": "  + to_string( min_second_deriv_significance )    + ",\n";
  result += "  \"compton_next_ratio_max\": "          + to_string( compton_next_ratio_max )           + ",\n";
  result += "  \"compton_prev_ratio_min\": "          + to_string( compton_prev_ratio_min )           + ",\n";
  result += "  \"compton_total_ratio_min\": "         + to_string( compton_total_ratio_min )          + ",\n";
  result += "  \"low_energy_test_max_keV\": "         + to_string( low_energy_test_max_keV )          + ",\n";
  result += "  \"low_energy_drop_fraction\": "        + to_string( low_energy_drop_fraction )         + ",\n";
  result += "  \"pcgap_feature_nsigma\": "            + to_string( pcgap_feature_nsigma )             + ",\n";
  result += "  \"pcgap_max_extent_nsigma\": "         + to_string( pcgap_max_extent_nsigma )          + ",\n";
  result += "  \"pcgap_roi_blend_weight\": "          + to_string( pcgap_roi_blend_weight )           + ",\n";
  result += "  \"pcgap_fom_blend_threshold\": "       + to_string( pcgap_fom_blend_threshold )        + "\n";
  result += "}";
  return result;
}//FindCandidateSettings::to_json


// ============================================================================
// LowHighResClassifySettings print/to_json
// ============================================================================

string LowHighResClassifySettings::print( const string &var_name ) const
{
  string result;
  result += var_name + ".hpge_fwhm_bias_mult = " + to_string( hpge_fwhm_bias_mult ) + ";\n";
  result += var_name + ".nai_fwhm_bias_mult = " + to_string( nai_fwhm_bias_mult ) + ";\n";
  result += var_name + ".hpge_distance_weight = " + to_string( hpge_distance_weight ) + ";\n";
  result += var_name + ".nai_distance_weight = " + to_string( nai_distance_weight ) + ";\n";
  result += var_name + ".amp_clamp_denom = " + to_string( amp_clamp_denom ) + ";\n";
  result += var_name + ".weight_clamp_max = " + to_string( weight_clamp_max ) + ";\n";
  result += var_name + ".min_peak_energy = " + to_string( min_peak_energy ) + ";\n";
  result += var_name + ".max_peak_energy = " + to_string( max_peak_energy ) + ";\n";
  result += var_name + ".unknown_threshold = " + to_string( unknown_threshold ) + ";\n";
  result += var_name + ".min_peaks_for_classify = " + to_string( min_peaks_for_classify ) + ";\n";
  result += var_name + ".min_peak_significance = " + to_string( min_peak_significance ) + ";\n";
  result += var_name + ".sig_fraction_of_max = " + to_string( sig_fraction_of_max ) + ";\n";
  result += var_name + ".narrow_penalty_mult = " + to_string( narrow_penalty_mult ) + ";\n";
  result += var_name + ".narrow_sig_threshold = " + to_string( narrow_sig_threshold ) + ";\n";
  result += var_name + ".best_peak_czt_penalty = " + to_string( best_peak_czt_penalty ) + ";\n";
  result += var_name + ".best_peak_conf_threshold = " + to_string( best_peak_conf_threshold ) + ";\n";
  return result;
}//LowHighResClassifySettings::print


string LowHighResClassifySettings::to_json() const
{
  string result = "{\n";
  result += "  \"hpge_fwhm_bias_mult\": " + to_string( hpge_fwhm_bias_mult ) + ",\n";
  result += "  \"nai_fwhm_bias_mult\": " + to_string( nai_fwhm_bias_mult ) + ",\n";
  result += "  \"hpge_distance_weight\": " + to_string( hpge_distance_weight ) + ",\n";
  result += "  \"nai_distance_weight\": " + to_string( nai_distance_weight ) + ",\n";
  result += "  \"amp_clamp_denom\": " + to_string( amp_clamp_denom ) + ",\n";
  result += "  \"weight_clamp_max\": " + to_string( weight_clamp_max ) + ",\n";
  result += "  \"min_peak_energy\": " + to_string( min_peak_energy ) + ",\n";
  result += "  \"max_peak_energy\": " + to_string( max_peak_energy ) + ",\n";
  result += "  \"unknown_threshold\": " + to_string( unknown_threshold ) + ",\n";
  result += "  \"min_peaks_for_classify\": " + to_string( min_peaks_for_classify ) + ",\n";
  result += "  \"min_peak_significance\": " + to_string( min_peak_significance ) + ",\n";
  result += "  \"sig_fraction_of_max\": " + to_string( sig_fraction_of_max ) + ",\n";
  result += "  \"narrow_penalty_mult\": " + to_string( narrow_penalty_mult ) + ",\n";
  result += "  \"narrow_sig_threshold\": " + to_string( narrow_sig_threshold ) + ",\n";
  result += "  \"best_peak_czt_penalty\": " + to_string( best_peak_czt_penalty ) + ",\n";
  result += "  \"best_peak_conf_threshold\": " + to_string( best_peak_conf_threshold ) + "\n";
  result += "}";
  return result;
}//LowHighResClassifySettings::to_json

}//namespace PeakFitSpec


// ============================================================================
// File-static helpers for candidate peak finding
// ============================================================================

namespace
{

// Two-pass log-parabola Gaussian fit for candidate peak sigma estimation.
// Reduces FWHM estimation error by ~49% vs zero-crossing method.
#define USE_TWO_PASS_SIMPLE_GAUSS_FIT 1


/** Fits a log-parabola (ln(y) = a + b*x + c*x^2) to continuum-subtracted spectrum data
 *  near a candidate peak to estimate sigma.
 *
 *  The curvature coefficient c (which determines sigma) is robust to continuum errors because
 *  it depends on the relative shape of the peak, not its absolute height.
 */
static double fit_log_parabola( const vector<float> &spectrum,
                                const shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                                const size_t nchannel,
                                const size_t peak_ch,
                                const double peak_energy,
                                const size_t cont_ref_ch,
                                const double lower_cont,
                                const double cont_slope,
                                double *area_out = nullptr )
{
  const double cont_at_peak = lower_cont + cont_slope * static_cast<double>( peak_ch - cont_ref_ch );
  const double peak_net = spectrum[peak_ch] - cont_at_peak;

  if( peak_net <= 1.0 )
    return -1.0;

  // Determine fit window: channels with net counts >= 5% of peak net counts
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
    const double ln_y = log( net_counts );
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

  // Solve 3x3 normal equations via Cramer's rule
  const double d00 = S,   d01 = Sx,  d02 = Sx2;
  const double d10 = Sx,  d11 = Sx2, d12 = Sx3;
  const double d20 = Sx2, d21 = Sx3, d22 = Sx4;

  const double det = d00 * (d11 * d22 - d12 * d21)
                   - d01 * (d10 * d22 - d12 * d20)
                   + d02 * (d10 * d21 - d11 * d20);

  if( fabs( det ) < 1.0e-30 )
    return -1.0;

  const double det_c = d00 * (d11 * Sx2y - Sxy * d21)
                     - d01 * (d10 * Sx2y - Sxy * d20)
                     + Sy  * (d10 * d21  - d11 * d20);

  const double c = det_c / det;

  // For a Gaussian: ln(y) = ln(A) - x^2 / (2*sigma^2), so c = -1/(2*sigma^2)
  if( c >= 0.0 )
    return -1.0;

  const double sigma_keV = sqrt( -1.0 / (2.0 * c) );

  if( sigma_keV < 0.1 || sigma_keV > 100.0 )
    return -1.0;

  if( area_out )
  {
    const double det_a = Sy  * (d11 * d22 - d12 * d21)
                       - d01 * (Sxy * d22 - Sx2y * d20)
                       + d02 * (Sxy * d21 - Sx2y * d11);
    const double a = det_a / det;
    const double peak_height = exp( a );

    const double keV_per_ch = energy_cal->energy_for_channel( static_cast<double>( peak_ch ) + 1.0 )
                            - energy_cal->energy_for_channel( static_cast<double>( peak_ch ) );

    if( keV_per_ch > 0.0 )
      *area_out = peak_height * sigma_keV * sqrt( 2.0 * boost::math::constants::pi<double>() ) / keV_per_ch;
    else
      *area_out = -1.0;
  }

  return sigma_keV;
}//fit_log_parabola(...)


static double estimate_sigma_from_raw_data( const shared_ptr<const SpecUtils::Measurement> &data,
                                            const double peak_mean_keV,
                                            const double roi_lower_energy,
                                            const double roi_upper_energy,
                                            double *area_out = nullptr )
{
  const shared_ptr<const SpecUtils::EnergyCalibration> energy_cal = data->energy_calibration();
  if( !energy_cal || !energy_cal->valid() || !data->gamma_counts() )
    return -1.0;

  const vector<float> &spectrum = *data->gamma_counts();
  const size_t nchannel = spectrum.size();
  if( nchannel < 20 )
    return -1.0;

  const size_t roi_lower_ch = static_cast<size_t>(
      max( 0.0, floor( energy_cal->channel_for_energy( roi_lower_energy ) ) ) );
  const size_t roi_upper_ch = static_cast<size_t>(
      min( static_cast<double>( nchannel - 1 ),
                ceil( energy_cal->channel_for_energy( roi_upper_energy ) ) ) );

  if( roi_upper_ch <= (roi_lower_ch + 2) )
    return -1.0;

  // Pass 1: Estimate continuum from channels just outside the ROI
  const size_t cont_channels = min( size_t(3), (roi_upper_ch - roi_lower_ch) / 4 + 1 );

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
  const size_t search_half = max( size_t(3), (roi_upper_ch - roi_lower_ch) / 3 );
  const size_t search_lower = (static_cast<size_t>( mean_channel ) > search_half)
                            ? (static_cast<size_t>( mean_channel ) - search_half) : 0;
  const size_t search_upper = min( static_cast<size_t>( mean_channel ) + search_half, nchannel - 1 );

  size_t peak_ch = static_cast<size_t>( max( 0.0, round( mean_channel ) ) );
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

  // Pass 2: Symmetric fit window to mitigate neighbor contamination
  const double keV_per_ch = (roi_upper_energy - roi_lower_energy) / roi_span;
  const double sigma_ch = pass1_sigma / keV_per_ch;

  const size_t sym_half = static_cast<size_t>( max( 1.0, round( 2.5 * sigma_ch ) ) );
  const size_t sym_lower = (peak_ch > sym_half) ? (peak_ch - sym_half) : 0;
  const size_t sym_upper = min( peak_ch + sym_half, nchannel - 1 );

  if( (sym_upper - sym_lower) >= 2 )
  {
    const size_t sym_cont_ch = min( size_t(3), sym_half );

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

      // Accept pass2 if it gives a smaller sigma (pass1 was inflated by neighbor)
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

}//anonymous namespace


// ============================================================================
// find_candidate_peaks
// ============================================================================

namespace PeakFitSpec
{

vector<PeakDef> find_candidate_peaks( const shared_ptr<const SpecUtils::Measurement> &data,
                                      size_t start_channel,
                                      size_t end_channel,
                                      const FindCandidateSettings &settings )
{
  vector<PeakDef> result_peaks;

  const size_t nchannel = data->num_gamma_channels();
  const size_t nFluxuate = settings.num_chan_fluctuate;

  if( nchannel < (2*nFluxuate + 2) || (nchannel < 4) )
    return {};

  if( start_channel >= nchannel )
    start_channel = 0;
  if( end_channel <= start_channel || end_channel >= (nchannel-1) )
    end_channel = nchannel - nFluxuate - 1;


  const double threshold_FOM = settings.threshold_FOM;
  const float pos_sum_threshold_sf = settings.pos_sum_threshold_sf;


  const int order = settings.smooth_polynomial_order;
  const size_t side_bins = settings.num_smooth_side_channels;

  const vector<float> &spectrum = *data->gamma_counts();

  vector<float> second_deriv, second_deriv_var;
  vector<float> smoothed_data, rougher_first_deriv, rougher_second_deriv, rougher_second_deriv_var;

  // Energy-adaptive SG 2nd derivative: window size scales with channel position to match
  //  peak resolution vs energy. Pre-computes smoothings for each integer window size needed,
  //  then linearly interpolates per-channel.
  if( (settings.smooth_scale_power > 0.001f)
     && (settings.smooth_ref_fraction > 0.001f)
     && (nchannel >= 64)
     && (side_bins >= 3) )
  {
    const float ref_channel = settings.smooth_ref_fraction * nchannel;
    const int base_side = static_cast<int>( side_bins );
    const int max_allowed_side = std::min( static_cast<int>( nchannel / 4 ), 30 );

    // Compute the max side bins needed at the highest channel
    int max_side_needed = base_side;
    {
      const float scale = std::pow( static_cast<float>( nchannel - 1 ) / ref_channel,
                                   settings.smooth_scale_power );
      max_side_needed = std::min( static_cast<int>( std::ceil( base_side * scale ) ),
                                 max_allowed_side );
    }

    // Pre-compute SG smoothings (with variance) for each integer window size
    map<int, pair<vector<float>, vector<float>>> sg_results;
    for( int s = base_side; s <= max_side_needed; ++s )
    {
      SavitzyGolayCoeffs sg( s, s, order, 2 );
      sg.smooth_with_variance( spectrum, sg_results[s].first, sg_results[s].second );
    }

    // Build per-channel second_deriv and second_deriv_var by interpolating between
    //  the two nearest integer window sizes.
    second_deriv.resize( nchannel, 0.0f );
    second_deriv_var.resize( nchannel, 0.0f );

    for( size_t ch = 0; ch < nchannel; ++ch )
    {
      float ideal_side;
      if( ch <= static_cast<size_t>( ref_channel ) )
        ideal_side = static_cast<float>( base_side );
      else
        ideal_side = base_side * std::pow( static_cast<float>( ch ) / ref_channel,
                                          settings.smooth_scale_power );

      ideal_side = std::max( ideal_side, static_cast<float>( base_side ) );
      ideal_side = std::min( ideal_side, static_cast<float>( max_side_needed ) );

      const int lo = static_cast<int>( std::floor( ideal_side ) );
      const int hi = std::min( lo + 1, max_side_needed );
      const float frac = ideal_side - lo;

      const vector<float> &lo_smooth = sg_results[lo].first;
      const vector<float> &lo_var = sg_results[lo].second;
      const vector<float> &hi_smooth = sg_results[hi].first;
      const vector<float> &hi_var = sg_results[hi].second;

      second_deriv[ch] = (1.0f - frac) * lo_smooth[ch] + frac * hi_smooth[ch];
      second_deriv_var[ch] = (1.0f - frac) * lo_var[ch] + frac * hi_var[ch];
    }//for( each channel )
  }else
  {
    // Fixed window (backward compatible, e.g. for HPGe or when scale_power == 0)
    SavitzyGolayCoeffs sg( static_cast<int>( side_bins ), static_cast<int>( side_bins ), order, 2 );
    sg.smooth_with_variance( spectrum, second_deriv, second_deriv_var );
  }

  // Smooth the data lightly for stat significance summing and PCGAP walk
  smoothSpectrum( spectrum, static_cast<int>( side_bins ), 1, 0, smoothed_data );

  smoothSpectrum( spectrum, 3, 3, 1, rougher_first_deriv );

  // Rougher 2nd derivative (3-side-channel) used for sigma refinement and ROI extent checks
  {
    SavitzyGolayCoeffs sgcoeffs( 3, 3, 3, 2 );
    sgcoeffs.smooth_with_variance( spectrum, rougher_second_deriv, rougher_second_deriv_var );
  }


  const double amp_fake_factor = 1.0;

  const vector<float> &energies = *data->gamma_channel_energies();

  size_t minbin = 0, firstzero = 0, secondzero = 0;
  float secondsum = 0.0f, minval = 9999999999.9f;

  // Usually there is some turn-on for the spectrum, so lets skip past that
  if( start_channel == 0 )
  {
    while( (start_channel+nFluxuate) < end_channel )
    {
      bool all_neg = true;
      for( size_t ind = start_channel; ind < (start_channel+nFluxuate); ++ind )
        all_neg = (second_deriv[ind] < -1.0E-9);

      if( all_neg )
        break;
      ++start_channel;
    }

    while( (start_channel+nFluxuate) < end_channel )
    {
      bool all_pos = true;
      for( size_t ind = start_channel; ind < (start_channel+nFluxuate); ++ind )
        all_pos = (second_deriv[ind] > 1.0E-9);

      if( all_pos )
        break;
      ++start_channel;
    }

    // Add a buffer past the turn-on equal to the smoothing radius, so that the
    //  first accepted peak isnt contaminated by the edge of the spectrum.
    start_channel = std::min( start_channel + side_bins, end_channel - nFluxuate - 1 );
  }//if( start_channel == 0 )


  for( size_t channel = start_channel; channel <= end_channel; ++channel )
  {
    double rougher_amplitude = -1.0;
    const double channel_energy = data->gamma_channel_center(channel);

    const float secondDeriv = second_deriv[channel];

    bool secondSumPositive = true;
    float positive2ndDerivSum = 0.0f, positive2ndDataSum = 0.0f;
    for( size_t i = 0; (i < nFluxuate) && ((channel+i) <= end_channel); ++i )
    {
      const bool above = (second_deriv[channel+i] > 0.0f);
      if( above )
      {
        positive2ndDerivSum += second_deriv[channel+i];
        positive2ndDataSum += smoothed_data[channel+i];
      }
      secondSumPositive &= above;
    }

    secondSumPositive &= (positive2ndDerivSum > pos_sum_threshold_sf*secondsum);

    if( secondSumPositive && (minval < 99999999999.9)
       && (secondsum!=0.0) && (firstzero>0)
       && ((channel-firstzero)>2) )
    {
      if( pos_sum_threshold_sf < 0.0 )
      {
        size_t trial_second_zero = channel;
        while( (trial_second_zero > minbin)
              && (trial_second_zero > 1)
              && (second_deriv[trial_second_zero-1]  > 0.0) )
        {
          trial_second_zero -= 1;
        }

        secondzero = ((trial_second_zero > minbin) && second_deriv[trial_second_zero-1] < 0.0)
                        ? trial_second_zero
                        : channel;
      }else
      {
        secondzero = channel;
      }

      double mean = data->gamma_channel_center(minbin);

      const double sigma_smoothed = 0.5*(data->gamma_channel_center(secondzero)
                                         - data->gamma_channel_center(firstzero));
      double sigma = sigma_smoothed;


      assert( second_deriv[firstzero+1] <= 0.0 );

      // The more-coarse smoothed value at `minbin` may not be negative
      size_t rough_index = (rougher_second_deriv[minbin] < rougher_second_deriv[minbin-1]) ? minbin : minbin-1;
      rough_index = (rougher_second_deriv[rough_index] < rougher_second_deriv[minbin+1]) ? rough_index : minbin+1;

      double rougher_FOM = 0.0;

      size_t roi_begin_index = 0, roi_end_index = 0;

      bool rougher_confident_multi_roi = false;
      const double rougher_scnd_drv = rougher_second_deriv[rough_index];


      if( rougher_scnd_drv < 0.0 )
      {
        size_t lower_pos_index = rough_index, upper_pos_index = rough_index;
        while( lower_pos_index > nFluxuate )
        {
          bool all_pos = true;
          for( size_t ind = lower_pos_index; all_pos && (ind > (lower_pos_index - nFluxuate)); --ind )
            all_pos = (rougher_second_deriv[ind] >= -numeric_limits<float>::epsilon());
          if( all_pos )
            break;

          lower_pos_index -= 1;
        }


        while( upper_pos_index < (end_channel - nFluxuate) )
        {
          bool all_pos = true;
          for( size_t ind = upper_pos_index; all_pos && (ind < (upper_pos_index + nFluxuate)); ++ind )
            all_pos = (rougher_second_deriv[ind] >= -numeric_limits<float>::epsilon());
          if( all_pos )
            break;

          upper_pos_index += 1;
        }

        assert( (rougher_second_deriv[lower_pos_index] >= 0)   || (lower_pos_index <= (nFluxuate+1)) );
        assert( (rougher_second_deriv[lower_pos_index+1] <= 0) || (lower_pos_index <= (nFluxuate+1)) );
        assert( (rougher_second_deriv[upper_pos_index] >= 0)   || ((upper_pos_index + nFluxuate + 1) >= end_channel) );
        assert( (rougher_second_deriv[upper_pos_index-1] <= 0) || (upper_pos_index <= (nFluxuate+1)) );

        const float lower_crossing = lower_pos_index
          + (-rougher_second_deriv[lower_pos_index+1] / (rougher_second_deriv[lower_pos_index] - rougher_second_deriv[lower_pos_index+1]));
        const float upper_crossing = upper_pos_index
          + (rougher_second_deriv[upper_pos_index-1] / (rougher_second_deriv[upper_pos_index] - rougher_second_deriv[upper_pos_index-1]));

        const double lower_sigma_energy = data->energy_calibration()->energy_for_channel(lower_crossing);
        const double upper_sigma_energy = data->energy_calibration()->energy_for_channel(upper_crossing);

        const double rougher_sigma = 0.5*(upper_sigma_energy - lower_sigma_energy);


        double rougher_secondsum = 0.0;
        size_t min_rough_index = minbin;
        for( size_t index = lower_pos_index + 1; index < upper_pos_index; ++index )
        {
          if( rougher_second_deriv[index] < 0.0 )
            rougher_secondsum += rougher_second_deriv[index];
          if( rougher_second_deriv[index] < rougher_second_deriv[min_rough_index] )
            min_rough_index = index;
        }

        const double rougher_deriv_sigma = 0.5*( upper_crossing - lower_crossing );
        const double rougher_part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                         boost::math::constants::e<double>() ) )
                              / ( rougher_deriv_sigma * rougher_deriv_sigma );
        rougher_amplitude = -amp_fake_factor * rougher_secondsum / rougher_part;


        double rough_data_area = 0.0;
        for( size_t i = lower_pos_index+1; i <= upper_pos_index-1; ++i )
          rough_data_area += smoothed_data[i];
        const double rough_est_std_dev = sqrt( max(rough_data_area,1.0) );


        rougher_FOM = 0.68*rougher_amplitude / rough_est_std_dev;


        roi_begin_index = lower_pos_index;
        roi_end_index = upper_pos_index;

        assert( rougher_second_deriv[lower_pos_index] >= 0 || (lower_pos_index < (2*nFluxuate + 1) ) );
        assert( (rougher_second_deriv[upper_pos_index] >= 0) || (upper_pos_index > (end_channel - 2*nFluxuate - 1) ) );

        for( ; (roi_begin_index > 0) && (rougher_second_deriv[roi_begin_index] > 0.0); --roi_begin_index )
        {
          if( (roi_begin_index > 1) && (rougher_first_deriv[roi_begin_index] < 0) && (rougher_first_deriv[roi_begin_index-1] < 0) )
            break;
        }

        for( ; ((roi_end_index+1) < end_channel) && (rougher_second_deriv[roi_end_index] > 0.0); ++roi_end_index )
        {
          if( ((roi_end_index + 1) < end_channel) && (rougher_first_deriv[roi_end_index] > 0.0) && (rougher_first_deriv[roi_end_index+1] > 0.0) )
            break;
        }


        if( (rougher_sigma < sigma)
           && !IsNan(rougher_sigma)
           && !IsInf(rougher_sigma)
           && ((upper_crossing - lower_crossing) > 1.5f) )
        {
          sigma = rougher_sigma;
          if( rougher_second_deriv[min_rough_index] < 0.0 )
          {
            mean = data->gamma_channel_center(min_rough_index);
          }
        }
      }else //if( rougher_second_deriv[minbin] < 0 )
      {
        roi_begin_index = firstzero;
        roi_end_index = secondzero;

        for( ; (roi_begin_index > 0) && (second_deriv[roi_begin_index] < 0.0); --roi_begin_index )
        {
        }

        for( ; ((roi_end_index+1) < end_channel) && (second_deriv[roi_end_index] < 0.0); ++roi_end_index )
        {
        }

        assert( second_deriv[roi_begin_index] >= 0.0 );
        assert( second_deriv[roi_end_index] >= 0.0 );

        for( ; (roi_begin_index > 0) && (second_deriv[roi_begin_index] > 0.0); --roi_begin_index )
        {
        }

        for( ; ((roi_end_index+1) < end_channel) && (second_deriv[roi_end_index] > 0.0); ++roi_end_index )
        {
        }


        const auto second_deriv_sig = [minbin, &rougher_second_deriv, &rougher_second_deriv_var]( const int i ) -> double {
          if( (minbin < i) || ((minbin + i) >= rougher_second_deriv_var.size()) || (rougher_second_deriv_var[minbin + i] <= 0.0) )
            return 0.0;
          return rougher_second_deriv[minbin+i]/sqrt(rougher_second_deriv_var[minbin+i]);
        };

        int num_sig_pos = 0;
        const double second_deriv_sig_thresh = 1.0;
        num_sig_pos += (second_deriv_sig(-1) > second_deriv_sig_thresh);
        num_sig_pos += (second_deriv_sig(0) > second_deriv_sig_thresh);
        num_sig_pos += (second_deriv_sig(1) > second_deriv_sig_thresh);

        const bool roi_center_is_dip = (num_sig_pos >= 2);

        if( roi_center_is_dip )
        {
          bool peak_lower = false, peak_upper = false;
          const double required_convex_limit = -1.5;

          for( size_t i = roi_begin_index; !peak_lower && (i < minbin); ++i )
          {
            const double val = rougher_second_deriv[i];
            const double var = (rougher_second_deriv_var[i] > 0.0) ? sqrt(rougher_second_deriv_var[i]) : 1.0;
            peak_lower = ((val / var) < required_convex_limit);
          }

          for( size_t i = minbin + 1; !peak_upper && (i <= roi_end_index); ++i )
          {
            const double val = rougher_second_deriv[i];
            const double var = (rougher_second_deriv_var[i] > 0.0) ? sqrt(rougher_second_deriv_var[i]) : 1.0;
            peak_upper = ((val / var) < required_convex_limit);
          }

          rougher_confident_multi_roi = (peak_lower && peak_upper);
        }

        if( roi_center_is_dip )
        {
          size_t min_2nd_deriv_index = firstzero;
          for( size_t i = firstzero + 1; i < secondzero; ++i )
          {
            if( rougher_second_deriv[i] < rougher_second_deriv[min_2nd_deriv_index] )
              min_2nd_deriv_index = i;
          }

          size_t rough_firstzero = firstzero, rough_secondzero = secondzero;
          for( size_t i = min_2nd_deriv_index + 1; i < secondzero; ++i )
          {
            if( signbit(rougher_second_deriv[i-1]) != signbit(rougher_second_deriv[i]) )
            {
              rough_secondzero = i;
              break;
            }
          }

          for( size_t i = min_2nd_deriv_index - 1; (i > 0) && (i > firstzero); ++i )
          {
            if( signbit(rougher_second_deriv[i]) != signbit(rougher_second_deriv[i+1]) )
            {
              rough_firstzero = i;
              break;
            }
          }

          const double rough_mean = data->gamma_channel_center(min_2nd_deriv_index);
          const double rough_sigma = 0.5*(data->gamma_channel_center(rough_secondzero)
                                             - data->gamma_channel_center(rough_firstzero));
          if( rough_sigma > 0.01 )
          {
            mean = rough_mean;
            sigma = rough_sigma;

            double rougher_secondsum = 0.0;
            for( size_t index = rough_firstzero; index <= rough_secondzero; ++index )
            {
              if( rougher_second_deriv[index] < 0.0 )
                rougher_secondsum += rougher_second_deriv[index];
            }
            const double rougher_part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                             boost::math::constants::e<double>() ) )
                                  / ( rough_sigma * rough_sigma );
            rougher_amplitude = -amp_fake_factor * rougher_secondsum / rougher_part;
          }
        }//if( roi_center_is_dip )

      }//if( rougher_second_deriv can be used ) / else

      assert( roi_begin_index != roi_end_index );

      // PCGAP-style ROI refinement: walk outward from the peak tracking background,
      //  and blend with derivative-based ROI boundaries.
      if( (settings.pcgap_roi_blend_weight > 0.01f)
         && (sigma > 0.01) )
      {
        const size_t max_low_ch = data->find_gamma_channel(
          static_cast<float>( std::max( 0.0, mean - settings.pcgap_max_extent_nsigma * sigma ) ) );
        const size_t max_high_ch = std::min(
          data->find_gamma_channel(
            static_cast<float>( mean + settings.pcgap_max_extent_nsigma * sigma ) ),
          end_channel );

        // Walk lower boundary from ~1.5 sigma below mean
        size_t pcgap_lower = max_low_ch;
        {
          const size_t walk_start = data->find_gamma_channel(
            static_cast<float>( mean - 1.5 * sigma ) );
          float bg_sum = 0.0f;
          size_t bg_count = 0;
          for( size_t j = 0; (j < 3) && ((walk_start + j) < nchannel); ++j )
          {
            bg_sum += smoothed_data[walk_start + j];
            bg_count++;
          }

          for( size_t ch = walk_start; (ch > max_low_ch) && (ch > 0); --ch )
          {
            const float bg = std::max( bg_sum / std::max( bg_count, size_t(1) ), 1.0f );
            const float threshold = bg + settings.pcgap_feature_nsigma * std::sqrt( bg );
            if( smoothed_data[ch] > threshold )
            {
              pcgap_lower = std::min( ch + 3, walk_start );
              break;
            }
            bg_sum += smoothed_data[ch];
            bg_count++;
          }
        }

        // Walk upper boundary from ~1.5 sigma above mean
        size_t pcgap_upper = max_high_ch;
        {
          const size_t walk_start = std::min(
            data->find_gamma_channel( static_cast<float>( mean + 1.5 * sigma ) ),
            end_channel );
          float bg_sum = 0.0f;
          size_t bg_count = 0;
          for( size_t j = 0; (j < 3) && (walk_start >= j); ++j )
          {
            bg_sum += smoothed_data[walk_start - j];
            bg_count++;
          }

          for( size_t ch = walk_start; (ch < max_high_ch) && (ch < nchannel); ++ch )
          {
            const float bg = std::max( bg_sum / std::max( bg_count, size_t(1) ), 1.0f );
            const float threshold = bg + settings.pcgap_feature_nsigma * std::sqrt( bg );
            if( smoothed_data[ch] > threshold )
            {
              pcgap_upper = (ch > 3) ? std::max( ch - 3, walk_start ) : walk_start;
              break;
            }
            bg_sum += smoothed_data[ch];
            bg_count++;
          }
        }

        // Blend derivative-based ROI with PCGAP ROI
        const float w = settings.pcgap_roi_blend_weight;
        roi_begin_index = static_cast<size_t>( (1.0f - w) * roi_begin_index + w * pcgap_lower + 0.5f );
        roi_end_index = static_cast<size_t>( (1.0f - w) * roi_end_index + w * pcgap_upper + 0.5f );

        roi_begin_index = std::max( roi_begin_index, size_t(0) );
        roi_end_index = std::min( roi_end_index, nchannel - 1 );
        if( roi_end_index <= roi_begin_index )
          roi_end_index = roi_begin_index + 1;
      }//if( PCGAP blending enabled )

      const float roi_start_energy = data->gamma_channel_lower( roi_begin_index );
      const float roi_end_energy = data->gamma_channel_upper( roi_end_index );

      const double deriv_sigma = 0.5*( secondzero - firstzero );
      const double part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                       boost::math::constants::e<double>() ) )
                            / ( deriv_sigma * deriv_sigma );
      const double amplitude = -amp_fake_factor * secondsum / part;


      double data_area = 0.0;
      for( size_t i = firstzero; i <= secondzero; ++i )
        data_area += smoothed_data[i];


      double est_std_dev = sqrt( max(data_area,1.0) );

      const double figure_of_merit = 0.68*amplitude / est_std_dev;


      bool passed_higher_scrutiny = true;
      if( figure_of_merit < settings.more_scrutiny_FOM_threshold
         && (settings.more_scrutiny_FOM_threshold >= threshold_FOM) )
      {
        passed_higher_scrutiny = ( (rougher_scnd_drv < 0.0) && (rougher_FOM >= settings.more_scrutiny_coarser_FOM) );
      }//if( figure_of_merit < settings.more_scrutiny_FOM_threshold )


      // For peaks in low-stat area, check straight-line chi2
      if( passed_higher_scrutiny && ((figure_of_merit < settings.more_scrutiny_FOM_threshold)
                                     || (amplitude < settings.amp_to_apply_line_test_below)) )
      {
        const size_t cont_est_channels = min( size_t(3), 1 + roi_end_index - roi_begin_index );
        const size_t lower_edge_start = (cont_est_channels < roi_begin_index)
                                      ? (roi_begin_index - cont_est_channels) : 0;
        const size_t upper_edge_end = ((roi_end_index + cont_est_channels) < end_channel)
                                      ? (roi_end_index + cont_est_channels)
                                      : ((end_channel - cont_est_channels) > roi_end_index) ? (end_channel - cont_est_channels)
                                                                                            : min( roi_end_index + side_bins, end_channel - 1 );

        size_t lower_nchannel = 0, upper_nchannel = 0;
        double lower_cnts_per_chnl = 0.0, upper_cnts_per_chnl = 0.0;
        for( size_t i = lower_edge_start; i < roi_begin_index; ++i, ++lower_nchannel  )
          lower_cnts_per_chnl += spectrum[i];
        for( size_t i = roi_end_index + 1; i <= upper_edge_end; ++i, ++upper_nchannel )
          upper_cnts_per_chnl += spectrum[i];

        assert( lower_nchannel == (roi_begin_index - lower_edge_start) );
        assert( (upper_nchannel == (upper_edge_end - roi_end_index)) || (upper_edge_end < roi_end_index) );

        lower_cnts_per_chnl /= (lower_nchannel ? lower_nchannel : size_t(1));
        upper_cnts_per_chnl /= (upper_nchannel ? upper_nchannel : size_t(1));
        const double diff_per_chnl = (upper_cnts_per_chnl - lower_cnts_per_chnl) / (1 + roi_end_index - roi_begin_index);

        double chi2_dof = 0.0, max_chi2 = 0.0;
        for( size_t i = roi_begin_index; i <= roi_end_index; ++i )
        {
          const double pred_chnl_nts = lower_cnts_per_chnl + ((i-roi_begin_index) + 0.5) * diff_per_chnl;
          const double dc = spectrum[i] - pred_chnl_nts;

          const double uncert2 = max( (pred_chnl_nts > 0.0) ? pred_chnl_nts : static_cast<double>(spectrum[i]), 1.0 );
          const double val = dc*dc / uncert2;

          if( spectrum[i] > pred_chnl_nts )
            max_chi2 = max( max_chi2, val );

          chi2_dof += val;
        }
        chi2_dof /= (roi_end_index - roi_begin_index);

        passed_higher_scrutiny = (max_chi2 > settings.more_scrutiny_min_dev_from_line);
      }//if( check Chi2 of region )


      if( (figure_of_merit < 5) && (rougher_scnd_drv >= 0.0) && !rougher_confident_multi_roi )
      {
        passed_higher_scrutiny = false;
      }


      // Compton backscatter asymmetry test: check that positive 2nd-derivative sums
      //  flanking the peak are reasonably symmetric (not a Compton edge).
      if( passed_higher_scrutiny
         && (settings.compton_next_ratio_max < 50.0f) )
      {
        const size_t nflux = std::max( nFluxuate, ((secondzero - firstzero) / side_bins) + size_t(1) );

        // Sum positive 2nd-derivative forward from secondzero
        float nextpositivesum = 0.0f;
        for( size_t fwd = secondzero; fwd <= end_channel; ++fwd )
        {
          bool all_neg = true;
          for( size_t j = 0; (j < nflux) && ((fwd + j) <= end_channel); ++j )
            all_neg &= (second_deriv[fwd + j] < 0.0f);
          if( all_neg )
            break;
          if( second_deriv[fwd] > 0.0f )
            nextpositivesum += second_deriv[fwd];
        }

        // Sum positive 2nd-derivative backward from firstzero
        float prevpositivesum = 0.0f;
        for( size_t bwd = firstzero; bwd > 0; --bwd )
        {
          bool all_neg = true;
          for( size_t j = 0; (j < nflux) && (bwd >= j + 1); ++j )
            all_neg &= (second_deriv[bwd - j] < 0.0f);
          if( all_neg )
            break;
          if( second_deriv[bwd] > 0.0f )
            prevpositivesum += second_deriv[bwd];
        }

        if( secondsum < 0.0f )
        {
          const float nextratio = -nextpositivesum / secondsum;
          const float prevratio = -prevpositivesum / secondsum;
          const bool compton_ok = (nextratio < settings.compton_next_ratio_max
                                   || prevratio > settings.compton_prev_ratio_min)
                                  && ((nextratio + prevratio) > settings.compton_total_ratio_min);
          if( !compton_ok )
            passed_higher_scrutiny = false;
        }
      }//if( Compton backscatter test enabled )


      // Low-energy drop-off validation: check that peaks below a threshold energy
      //  have a Gaussian-like shape, not a detector efficiency edge.
      if( passed_higher_scrutiny
         && (settings.low_energy_test_max_keV > 0.0f)
         && (mean < settings.low_energy_test_max_keV)
         && (sigma > 0.01) && (minbin > 0) )
      {
        const size_t p15sig_ch = data->find_gamma_channel(
          static_cast<float>( mean + 1.5 * sigma ) );
        if( (p15sig_ch < nchannel) && (p15sig_ch > minbin) )
        {
          // Expected Gaussian drop from center to +1.5 sigma
          //  gaussian(x) = exp(-0.5 * ((x-mean)/sigma)^2), so at 1.5 sigma: exp(-0.5*2.25) = exp(-1.125) ≈ 0.325
          //  The ratio of Gaussian at +1.5σ to center is ~0.325, so expected drop is ~0.675 of center value.
          const double gauss_ratio_at_15sig = std::exp( -0.5 * 1.5 * 1.5 );
          const double expected_drop_fraction = 1.0 - gauss_ratio_at_15sig;

          const float peak_counts = smoothed_data[minbin];
          const float side_counts = smoothed_data[p15sig_ch];

          if( peak_counts > 1.0f )
          {
            const double actual_drop_fraction = (peak_counts - side_counts) / peak_counts;
            if( actual_drop_fraction < (settings.low_energy_drop_fraction * expected_drop_fraction) )
              passed_higher_scrutiny = false;
          }
        }
      }//if( low-energy drop-off test enabled )


      // Second-derivative significance cut: reject if the 2nd derivative minimum
      //  is not statistically significant relative to Poisson noise propagated
      //  through the SG smoothing coefficients.
      if( (settings.min_second_deriv_significance > 0.0f)
         && (minbin < second_deriv_var.size())
         && (second_deriv_var[minbin] > 0.0f) )
      {
        const float sig = std::abs( second_deriv[minbin] )
                        / std::sqrt( second_deriv_var[minbin] );
        if( sig < settings.min_second_deriv_significance )
        {
          passed_higher_scrutiny = false;
          // Also skip if even high-FOM: a statistically insignificant dip is noise
          if( figure_of_merit <= threshold_FOM )
          {
            secondsum = 0.0;
            minval = 9999999999.9f;
            minbin = secondzero = firstzero = 0;
            continue;
          }
        }
      }//if( 2nd derivative significance check )


      if( (figure_of_merit > threshold_FOM) || passed_higher_scrutiny )
      {
#if( USE_TWO_PASS_SIMPLE_GAUSS_FIT )
        {
          const double gf_sigma = estimate_sigma_from_raw_data( data, mean, roi_start_energy, roi_end_energy );
          if( gf_sigma > 0.0 )
            sigma = gf_sigma;
        }
#endif

        PeakDef peak( mean, sigma, amplitude );
        peak.continuum()->setRange( roi_start_energy, roi_end_energy );

        vector<double> cont_equation{0.0, 0.0};
        PeakContinuum::eqn_from_offsets( roi_begin_index, roi_end_index, mean, data,
                                        settings.num_smooth_side_channels,
                                        settings.num_smooth_side_channels,
                                        cont_equation[1], cont_equation[0] );
        peak.continuum()->setType( PeakContinuum::OffsetType::Linear );
        peak.continuum()->setParameters( mean, cont_equation, {} );

        result_peaks.push_back( std::move(peak) );
      }//if( passed )

      secondsum = 0.0;
      minval = 9999999999.9f;
      minbin = secondzero = firstzero = 0;
    }else
    {
      bool belowzero = true, goingnegative = true, abovezero = true;
      for( size_t i = 0; i < nFluxuate; ++i )
      {
        if( (channel+i+1) < end_channel )
          goingnegative &= (second_deriv[channel+i+1] < 0.0f);
        if( channel >= i )
        {
          belowzero &= (second_deriv[channel-i] <= 0.0f);
          abovezero &= (second_deriv[channel-i] > 0.0f);
        }
      }//for( size_t i = 0; i < nFluxuate; ++i )

      if( channel && !firstzero && goingnegative )
      {
        firstzero = channel;
        minbin = channel;
        minval = secondDeriv;

        for( size_t i = 1; i < nFluxuate; ++i )
          if( channel >= i )
            secondsum += second_deriv[channel-i];
      }else if( secondSumPositive )
      {
        secondsum = 0.0;
        minval = 9999999999.9f;
        minbin = secondzero = firstzero = 0;
      }

      if( firstzero > 0 )
      {
        secondsum += secondDeriv;

        if( secondDeriv < minval )
        {
          minbin = channel;
          minval = secondDeriv;
        }
      }//if( firstzero > 0 )
    }//if( we are out of region of interest) / else( in region of )
  }//for( loop over bins )

  return result_peaks;
}//find_candidate_peaks(...)


// ============================================================================
// classify_lowres_highres_from_peak_candidates
// ============================================================================

SpecClassType classify_lowres_highres_from_peak_candidates(
  const vector<PeakDef> &candidate_peaks,
  const LowHighResClassifySettings &settings,
  double *confidence_out )
{

  struct PeakVote
  {
    double w;
    double hpge_frac;
    double nai_frac;
  };

  // First pass: find max significance among energy-filtered peaks
  double max_significance = 0.0;
  for( const PeakDef &p : candidate_peaks )
  {
    if( !p.gausPeak() )
      continue;
    const double energy = p.mean();
    if( (energy < settings.min_peak_energy) || (energy > settings.max_peak_energy) )
      continue;
    const double area = max( p.peakArea(), 0.0 );
    const double significance = sqrt( area );
    max_significance = max( max_significance, significance );
  }

  const double effective_min_sig = max( settings.min_peak_significance,
                                        settings.sig_fraction_of_max * max_significance );

  vector<PeakVote> peak_votes;
  peak_votes.reserve( candidate_peaks.size() );

  double hpge_w = 0.0, nai_w = 0.0;
  int num_passing = 0;

  // Track the most significant peak for the best-peak CZT check
  double best_peak_significance = 0.0;
  double best_peak_energy = 0.0;
  double best_peak_fwhm = 0.0;

  for( const PeakDef &p : candidate_peaks )
  {
    if( !p.gausPeak() )
      continue;

    const double energy = p.mean();
    if( (energy < settings.min_peak_energy) || (energy > settings.max_peak_energy) )
      continue;

    const double area = max( p.peakArea(), 0.0 );
    const double significance = sqrt( area );
    if( significance < effective_min_sig )
      continue;

    num_passing += 1;

    // Track most significant peak for best-peak CZT check
    if( significance > best_peak_significance )
    {
      best_peak_significance = significance;
      best_peak_energy = energy;
      best_peak_fwhm = p.fwhm();
    }

    // Significance-based weighting
    const double w = min( significance / settings.amp_clamp_denom, settings.weight_clamp_max );

    // Compute biased reference FWHMs
    const double ref_hpge = settings.hpge_fwhm_bias_mult * PeakFitUtils::hpge_fwhm_fcn( static_cast<float>( energy ) );
    const double ref_nai  = settings.nai_fwhm_bias_mult  * PeakFitUtils::nai_fwhm_fcn( static_cast<float>( energy ) );

    const double measured_fwhm = p.fwhm();

    // Log-space distance to HPGe reference, with asymmetric penalty for narrow peaks
    const double log_ratio_hpge = log( measured_fwhm / ref_hpge );
    double dist_hpge = fabs( log_ratio_hpge ) * settings.hpge_distance_weight;

    if( log_ratio_hpge < 0.0 )
    {
      const double sig_mod = max( 0.0, 1.0 - significance / settings.narrow_sig_threshold );
      dist_hpge *= (1.0 + settings.narrow_penalty_mult * sig_mod);
    }

    const double dist_nai = fabs( log( measured_fwhm / ref_nai ) ) * settings.nai_distance_weight;

    // Soft voting: inverse distance squared
    const double eps = 1.0e-6;
    const double inv_hpge = 1.0 / (dist_hpge * dist_hpge + eps);
    const double inv_nai  = 1.0 / (dist_nai  * dist_nai  + eps);
    const double inv_sum = inv_hpge + inv_nai;

    const double hpge_frac = inv_hpge / inv_sum;
    const double nai_frac  = inv_nai  / inv_sum;

    hpge_w += w * hpge_frac;
    nai_w  += w * nai_frac;

    peak_votes.push_back( { w, hpge_frac, nai_frac } );
  }//for( candidate_peaks )


  if( num_passing < settings.min_peaks_for_classify )
  {
    if( confidence_out )
      *confidence_out = 0.0;
    return SpecClassType::Unknown;
  }

  const double total_w = hpge_w + nai_w;
  if( total_w <= 0.0 )
  {
    if( confidence_out )
      *confidence_out = 0.0;
    return SpecClassType::Unknown;
  }

  // First-pass winner (0=HPGe, 1=NaI)
  const int winner = (hpge_w >= nai_w) ? 0 : 1;

  // Two-pass outlier rejection: only when >4 peaks, exclude at most 25%
  if( num_passing > 4 )
  {
    const size_t max_exclude = static_cast<size_t>( num_passing ) / 4;

    vector<pair<double, size_t>> outliers;

    for( size_t i = 0; i < peak_votes.size(); ++i )
    {
      const PeakVote &pv = peak_votes[i];
      const double winner_frac = (winner == 0) ? pv.hpge_frac : pv.nai_frac;
      const double peak_best_frac = max( pv.hpge_frac, pv.nai_frac );

      if( winner_frac < peak_best_frac )
        outliers.push_back( { winner_frac, i } );
    }

    if( !outliers.empty() )
    {
      sort( outliers.begin(), outliers.end() );

      const size_t num_exclude = min( outliers.size(), max_exclude );

      hpge_w = 0.0;
      nai_w = 0.0;

      vector<bool> excluded( peak_votes.size(), false );
      for( size_t i = 0; i < num_exclude; ++i )
        excluded[outliers[i].second] = true;

      for( size_t i = 0; i < peak_votes.size(); ++i )
      {
        if( excluded[i] )
          continue;

        const PeakVote &pv = peak_votes[i];
        hpge_w += pv.w * pv.hpge_frac;
        nai_w  += pv.w * pv.nai_frac;
      }

    }
  }//if( num_passing > 4 ) -- outlier rejection

  // Best-peak CZT check
  if( (best_peak_significance > 0.0) && (settings.best_peak_czt_penalty > 0.0) && (hpge_w > nai_w) )
  {
    const double hpge_frac = hpge_w / (hpge_w + nai_w);
    if( hpge_frac < settings.best_peak_conf_threshold )
    {
      const double ref_hpge = settings.hpge_fwhm_bias_mult
        * PeakFitUtils::hpge_fwhm_fcn( static_cast<float>( best_peak_energy ) );
      const double ref_czt = PeakFitUtils::czt_fwhm_fcn( static_cast<float>( best_peak_energy ) );

      const double dist_hpge = fabs( log( best_peak_fwhm / ref_hpge ) );
      const double dist_czt  = fabs( log( best_peak_fwhm / ref_czt ) );

      if( dist_czt < dist_hpge )
      {
        const double transfer = settings.best_peak_czt_penalty * (dist_hpge - dist_czt) * best_peak_significance;
        const double actual_transfer = min( transfer, hpge_w * 0.5 );
        hpge_w -= actual_transfer;
        nai_w  += actual_transfer;
      }
    }
  }//best-peak CZT check

  const double final_total = hpge_w + nai_w;
  if( final_total <= 0.0 )
  {
    if( confidence_out )
      *confidence_out = 0.0;
    return SpecClassType::Unknown;
  }

  const double final_max = max( hpge_w, nai_w );
  const double vote_fraction = final_max / final_total;

  if( confidence_out )
    *confidence_out = vote_fraction;

  if( vote_fraction < settings.unknown_threshold )
    return SpecClassType::Unknown;

  if( final_max == hpge_w )
    return SpecClassType::High;

  return SpecClassType::LowOrMedRes;
}//classify_lowres_highres_from_peak_candidates


// ============================================================================
// initial_lowres_highres_classify
// ============================================================================

SpecClassType initial_lowres_highres_classify(
  const shared_ptr<const SpecUtils::Measurement> &measurement,
  double *confidence_out )
{
  if( !measurement || (measurement->num_gamma_channels() < 16) )
  {
    if( confidence_out )
      *confidence_out = 0.0;
    return SpecClassType::Unknown;
  }

  // GA-optimized settings for HPGe vs LowOrMedRes classification.
  // Optimized 2026-03-01 over 119,812 spectra from diverse detector types.
  // Performance: 97.0% correct, 0.5% incorrect, 2.5% unknown
  //   (correct=HPGe classified as High, or non-HPGe classified as LowOrMedRes;
  //    incorrect=HPGe classified as LowOrMedRes, or non-HPGe classified as High)
  // See target/peak_fit_improve/ for GA optimization code (ClassifyDetType_GA).

  FindCandidateSettings cand_settings;
  cand_settings.num_smooth_side_channels = 9;
  cand_settings.smooth_polynomial_order = 2;
  cand_settings.threshold_FOM = 0.758621;
  cand_settings.more_scrutiny_FOM_threshold = 1.598265;
  cand_settings.pos_sum_threshold_sf = 0.119178f;
  cand_settings.num_chan_fluctuate = 1;
  cand_settings.more_scrutiny_coarser_FOM = 3.001943f;
  cand_settings.more_scrutiny_min_dev_from_line = 6.816465;
  cand_settings.amp_to_apply_line_test_below = 6.000000;

  LowHighResClassifySettings cls_settings;
  cls_settings.hpge_fwhm_bias_mult = 0.493;
  cls_settings.nai_fwhm_bias_mult = 1.121;
  cls_settings.hpge_distance_weight = 4.644;
  cls_settings.nai_distance_weight = 2.266;
  cls_settings.amp_clamp_denom = 7.419;
  cls_settings.weight_clamp_max = 21.120;
  cls_settings.min_peak_energy = 44.267;
  cls_settings.max_peak_energy = 2973.993;
  cls_settings.unknown_threshold = 0.559;
  cls_settings.min_peaks_for_classify = 1;
  cls_settings.min_peak_significance = 2.251;
  cls_settings.sig_fraction_of_max = 0.009;
  cls_settings.narrow_penalty_mult = 1.208;
  cls_settings.narrow_sig_threshold = 65.150;
  cls_settings.best_peak_czt_penalty = 4.545;
  cls_settings.best_peak_conf_threshold = 0.642;

  const size_t nchannel = measurement->num_gamma_channels();
  const vector<PeakDef> candidates = find_candidate_peaks( measurement, 0, nchannel, cand_settings );

  return classify_lowres_highres_from_peak_candidates( candidates, cls_settings, confidence_out );
}//initial_lowres_highres_classify

}//namespace PeakFitSpec
