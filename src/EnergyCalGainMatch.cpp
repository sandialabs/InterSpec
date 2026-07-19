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
#include <set>
#include <cmath>
#include <deque>
#include <tuple>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <utility>
#include <algorithm>
#include <functional>
#include <stdexcept>

#include <boost/math/tools/minima.hpp>

#include <Wt/WText.h>
#include <Wt/WColor.h>
#include <Wt/WLabel.h>
#include <Wt/WCheckBox.h>
#include <Wt/WComboBox.h>
#include <Wt/WPushButton.h>
#include <Wt/WServer.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>
#include <Wt/WButtonGroup.h>
#include <Wt/WRadioButton.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/EnergyCalGainMatch.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

using namespace Wt;
using namespace std;

namespace
{
  /** Log-energy grid step; a 0.2% gain step per bin. */
  const double sm_log_step = std::log( 1.002 );

  /** The cross-correlation searches gains in [1/sm_max_gain_search, sm_max_gain_search]. */
  const double sm_max_gain_search = 8.0;

  /** Log-energy grids never extend below this energy (keV). */
  const float sm_min_grid_energy = 10.0f;

  /** Correlations below this get a WarningType::WeakCorrelation.
   Genuine matches between same-system detectors typically correlate above ~0.9 (after the
   high-pass), while spurious locks (e.g., true gain outside the search window) come in
   well below ~0.5.
   */
  const double sm_weak_correlation = 0.5;

  /** Spectra with fewer counts than this in the energy range are excluded from matching. */
  const double sm_min_counts_in_range = 100.0;

  /** Minimum number of overlapping log-grid bins for a gain hypothesis to be considered. */
  const size_t sm_min_overlap_bins = 20;

  /** Chi2 value returned for invalid/unusable (gain,offset) hypotheses. */
  const double sm_invalid_chi2 = 1.0E30;


  /** Rebins the spectrum onto the given (ascending) energy edges, returning
   sqrt( counts-per-keV ) for each of the (edges.size()-1) resulting bins.
   Bins outside the spectrums energy range come back as zero.
   */
  vector<float> sqrt_density_on_edges( const SpecUtils::Measurement &spec,
                                       const vector<float> &edges )
  {
    const shared_ptr<const SpecUtils::EnergyCalibration> cal = spec.energy_calibration();
    const shared_ptr<const vector<float>> &counts = spec.gamma_counts();
    assert( cal && cal->valid() && counts );

    vector<float> rebinned;
    SpecUtils::rebin_by_lower_edge( *cal->channel_energies(), *counts, edges, rebinned );

    vector<float> density( edges.size() - 1, 0.0f );
    const size_t nbin = std::min( density.size(), rebinned.size() );
    for( size_t i = 0; i < nbin; ++i )
    {
      const float width = edges[i+1] - edges[i];
      const float d = (width > 0.0f) ? (std::max(rebinned[i], 0.0f) / width) : 0.0f;
      density[i] = std::sqrt( d );
    }

    return density;
  }//sqrt_density_on_edges(...)


  /** Subtracts a centered moving average (half-width hw bins) from vals, leaving only
   structure narrower than the window.  Without this, smoothly falling spectra (e.g., plastic
   scintillator, or continuum-dominated data) are nearly self-similar under log-energy shifts,
   and their overall trend correlates strongly at any lag, swamping the actual features.
   */
  void highpass_baseline_subtract( vector<float> &vals, const int hw )
  {
    const int n = static_cast<int>( vals.size() );
    if( (n < 3) || (hw < 1) )
      return;

    vector<double> cumulative( n + 1, 0.0 );
    for( int i = 0; i < n; ++i )
      cumulative[i+1] = cumulative[i] + vals[i];

    vector<float> answer( n );
    for( int i = 0; i < n; ++i )
    {
      const int lo = std::max( 0, i - hw );
      const int hi = std::min( n - 1, i + hw );
      const double avg = (cumulative[hi+1] - cumulative[lo]) / (hi - lo + 1);
      answer[i] = vals[i] - static_cast<float>( avg );
    }

    vals.swap( answer );
  }//highpass_baseline_subtract(...)


  /** Evaluates fcn at (2*half_steps + 1) points, center +- i*step, returning the x with the
   smallest value.  Used to localize narrow chi2 valleys before polishing with Brent.
   */
  template <class Fcn>
  double scan_minimum( const Fcn &fcn, const double center, const double step,
                       const int half_steps )
  {
    double best_x = center;
    double best_val = fcn( center );

    for( int i = -half_steps; i <= half_steps; ++i )
    {
      if( !i )
        continue;

      const double x = center + i*step;
      const double val = fcn( x );
      if( val < best_val )
      {
        best_val = val;
        best_x = x;
      }
    }

    return best_x;
  }//scan_minimum(...)
}//namespace


namespace GainMatchCalc
{

std::shared_ptr<SpecUtils::EnergyCalibration>
transform_calibration( const std::shared_ptr<const SpecUtils::EnergyCalibration> &cal,
                       const double gain, const double offset )
{
  if( !cal || !cal->valid() )
    throw runtime_error( "transform_calibration: invalid input calibration" );

  if( (gain <= 0.0) || std::isnan(gain) || std::isinf(gain)
      || std::isnan(offset) || std::isinf(offset) )
    throw runtime_error( "transform_calibration: invalid gain or offset" );

  const size_t nchannel = cal->num_channels();
  auto answer = make_shared<SpecUtils::EnergyCalibration>();

  switch( cal->type() )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    {
      // Energy is linear in the coefficients for both polynomial and FRF (including the FRF
      //  low-energy 1/(1+60x) term), so E'=gain*E+offset is exactly a coefficient transform.
      vector<float> coefs = cal->coefficients();
      assert( !coefs.empty() );
      for( float &c : coefs )
        c = static_cast<float>( gain * c );
      coefs[0] += static_cast<float>( offset );

      if( cal->type() == SpecUtils::EnergyCalType::FullRangeFraction )
        answer->set_full_range_fraction( nchannel, coefs, cal->deviation_pairs() );
      else if( cal->type() == SpecUtils::EnergyCalType::Polynomial )
        answer->set_polynomial( nchannel, coefs, cal->deviation_pairs() );
      else
        answer->set_default_polynomial( nchannel, coefs, cal->deviation_pairs() );
      break;
    }//case polynomial or FRF

    case SpecUtils::EnergyCalType::LowerChannelEdge:
    {
      const shared_ptr<const vector<float>> &energies = cal->channel_energies();
      assert( energies && !energies->empty() );

      vector<float> new_energies( energies->size() );
      for( size_t i = 0; i < energies->size(); ++i )
        new_energies[i] = static_cast<float>( gain * (*energies)[i] + offset );

      answer->set_lower_channel_energy( nchannel, std::move(new_energies) );
      break;
    }//case LowerChannelEdge

    case SpecUtils::EnergyCalType::InvalidEquationType:
      throw runtime_error( "transform_calibration: invalid equation type" );
  }//switch( cal->type() )

  return answer;
}//transform_calibration(...)


double counts_in_range( const std::shared_ptr<const SpecUtils::Measurement> &spec,
                        const float lower_energy, const float upper_energy )
{
  if( !spec )
    return 0.0;

  const shared_ptr<const SpecUtils::EnergyCalibration> cal = spec->energy_calibration();
  const shared_ptr<const vector<float>> &counts = spec->gamma_counts();
  if( !cal || !cal->valid() || !counts || counts->empty() )
    return 0.0;

  const vector<float> &edges = *cal->channel_energies();
  const bool no_upper = (upper_energy <= lower_energy);

  double sum = 0.0;
  const size_t nchannel = std::min( counts->size(), edges.size() - 1 );
  for( size_t i = 0; i < nchannel; ++i )
  {
    const float e = edges[i];
    if( (e >= lower_energy) && (no_upper || (e < upper_energy)) )
      sum += (*counts)[i];
  }

  return sum;
}//counts_in_range(...)


double coarse_gain_xcorr( const std::shared_ptr<const SpecUtils::Measurement> &test,
                          const std::shared_ptr<const SpecUtils::Measurement> &reference,
                          const float lower_energy, const float upper_energy,
                          double &correlation,
                          std::vector<WarningType> &warnings )
{
  correlation = 0.0;

  if( !test || !reference )
    throw runtime_error( "coarse_gain_xcorr: null spectrum" );

  const shared_ptr<const SpecUtils::EnergyCalibration> ref_cal = reference->energy_calibration();
  const shared_ptr<const SpecUtils::EnergyCalibration> test_cal = test->energy_calibration();
  if( !ref_cal || !ref_cal->valid() || !reference->gamma_counts()
      || !test_cal || !test_cal->valid() || !test->gamma_counts() )
    throw runtime_error( "coarse_gain_xcorr: invalid spectrum" );

  const bool no_upper = (upper_energy <= lower_energy);
  const float win_lo = std::max( { lower_energy, sm_min_grid_energy, ref_cal->lower_energy() } );
  const float win_hi = no_upper ? ref_cal->upper_energy()
                                : std::min( upper_energy, ref_cal->upper_energy() );

  if( (win_lo <= 0.0f) || (win_hi <= (1.15f * win_lo)) )
    throw runtime_error( "coarse_gain_xcorr: usable energy range too narrow" );

  const int nwin = static_cast<int>( std::ceil( std::log(win_hi / win_lo) / sm_log_step ) );
  const int nsearch = static_cast<int>( std::ceil( std::log(sm_max_gain_search) / sm_log_step ) );

  // High-pass window half-width: structure wider than ~15% in energy is treated as smooth
  //  continuum trend and removed before correlating.
  const int highpass_hw = static_cast<int>( std::ceil( std::log(1.15) / sm_log_step ) );

  if( (nwin - 2*highpass_hw) < static_cast<int>(sm_min_overlap_bins) )
    throw runtime_error( "coarse_gain_xcorr: usable energy range too narrow" );

  // Reference bins i = 0..nwin-1, with lower edges win_lo*exp(i*step)
  vector<float> ref_edges( nwin + 1 );
  for( int i = 0; i <= nwin; ++i )
    ref_edges[i] = win_lo * static_cast<float>( std::exp( i * sm_log_step ) );
  vector<float> ref_vals = sqrt_density_on_edges( *reference, ref_edges );
  highpass_baseline_subtract( ref_vals, highpass_hw );

  // Test lattice: index j has lower edge at log-energy log(win_lo) + (j - nsearch)*step,
  //  j = 0..(nwin + 2*nsearch - 1), clamped to the range the test spectrum covers.
  const int nlattice = nwin + 2*nsearch;
  const double test_lo = test_cal->lower_energy();
  const double test_hi = test_cal->upper_energy();

  auto lattice_energy = [win_lo,nsearch]( const int j ) -> double {
    return win_lo * std::exp( (j - nsearch) * sm_log_step );
  };

  int j_first = 0, j_last = nlattice - 1;  //inclusive lattice bin indices
  while( (j_first < nlattice) && (lattice_energy(j_first) < test_lo) )
    ++j_first;
  while( (j_last >= 0) && (lattice_energy(j_last + 1) > test_hi) )
    --j_last;

  if( (j_last - j_first) < static_cast<int>(sm_min_overlap_bins) )
    throw runtime_error( "coarse_gain_xcorr: no overlap between spectra energy ranges" );

  vector<float> test_edges( j_last - j_first + 2 );
  for( int j = j_first; j <= (j_last + 1); ++j )
    test_edges[j - j_first] = static_cast<float>( lattice_energy(j) );
  vector<float> test_vals = sqrt_density_on_edges( *test, test_edges );
  highpass_baseline_subtract( test_vals, highpass_hw );

  // The high-pass leaves systematic artifacts within a half-window of each array end (the
  //  moving average is one-sided there), and those artifacts correlate between the two
  //  arrays at whatever lag aligns the array ends - so trim them from the comparison.
  const int i_first = highpass_hw;
  const int i_last = nwin - 1 - highpass_hw;
  const int jv_first = j_first + highpass_hw;
  const int jv_last = j_last - highpass_hw;

  // For gain hypothesis g = exp(lag*step), reference bin i pairs with test lattice
  //  bin j = i + nsearch - lag (a feature at reference energy E sits at E/g in the test).
  // Requiring a healthy fraction of the window to overlap keeps a few continuum bins from
  //  producing spuriously high correlations.
  const size_t min_overlap = std::max( sm_min_overlap_bins, static_cast<size_t>(nwin)/3 );

  int best_lag = 0;
  double best_corr = -2.0;
  vector<double> corrs( 2*nsearch + 1, -2.0 );

  for( int lag = -nsearch; lag <= nsearch; ++lag )
  {
    double sum_r = 0.0, sum_t = 0.0, sum_rr = 0.0, sum_tt = 0.0, sum_rt = 0.0;
    size_t num = 0;

    for( int i = i_first; i <= i_last; ++i )
    {
      const int j = i + nsearch - lag;
      if( (j < jv_first) || (j > jv_last) )
        continue;

      const double r = ref_vals[i];
      const double t = test_vals[j - j_first];
      sum_r += r;
      sum_t += t;
      sum_rr += r*r;
      sum_tt += t*t;
      sum_rt += r*t;
      ++num;
    }//for( loop over reference bins )

    if( num < min_overlap )
      continue;

    const double n = static_cast<double>( num );
    const double cov = sum_rt - sum_r*sum_t/n;
    const double var_r = sum_rr - sum_r*sum_r/n;
    const double var_t = sum_tt - sum_t*sum_t/n;
    if( (var_r <= 0.0) || (var_t <= 0.0) )
      continue;

    const double corr = cov / std::sqrt( var_r * var_t );
    corrs[lag + nsearch] = corr;

    if( corr > best_corr )
    {
      best_corr = corr;
      best_lag = lag;
    }
  }//for( loop over lags )

  if( best_corr < -1.5 )
    throw runtime_error( "coarse_gain_xcorr: no gain hypothesis had enough overlapping data" );

  correlation = best_corr;
  if( best_corr < sm_weak_correlation )
    warnings.push_back( WarningType::WeakCorrelation );
  if( (best_lag == -nsearch) || (best_lag == nsearch) )
    warnings.push_back( WarningType::CorrelationAtEdge );

  // If a competing maximum well away from the best lag is nearly as good, the match is
  //  ambiguous (e.g., the true gain is outside the search window, or the spectra are
  //  self-similar), so flag it.
  const int exclusion = static_cast<int>( std::ceil( std::log(1.05) / sm_log_step ) );
  double second_best = -2.0;
  for( size_t idx = 0; idx < corrs.size(); ++idx )
  {
    if( (std::abs(static_cast<int>(idx) - (best_lag + nsearch)) > exclusion)
        && (corrs[idx] > second_best) )
      second_best = corrs[idx];
  }
  if( (second_best > -1.5) && ((best_corr - second_best) < 0.05) )
    warnings.push_back( WarningType::AmbiguousCorrelation );

  // Parabolic interpolation of the correlation peak for sub-bin gain resolution
  double lag_refined = best_lag;
  const int c = best_lag + nsearch;
  if( (c > 0) && (c < 2*nsearch) && (corrs[c-1] > -1.5) && (corrs[c+1] > -1.5) )
  {
    const double denom = corrs[c-1] - 2.0*corrs[c] + corrs[c+1];
    if( denom < 0.0 )
    {
      double delta = 0.5 * (corrs[c-1] - corrs[c+1]) / denom;
      delta = std::max( -0.5, std::min(0.5, delta) );
      lag_refined += delta;
    }
  }//if( can parabolic interpolate )

  return std::exp( lag_refined * sm_log_step );
}//coarse_gain_xcorr(...)


double chi2_for_gain_offset( const std::shared_ptr<const SpecUtils::Measurement> &test,
                             const std::shared_ptr<const SpecUtils::Measurement> &reference,
                             const float lower_energy, const float upper_energy,
                             const double gain, const double offset,
                             size_t &num_channels_used )
{
  num_channels_used = 0;

  if( !test || !reference )
    throw runtime_error( "chi2_for_gain_offset: null spectrum" );

  const shared_ptr<const SpecUtils::EnergyCalibration> ref_cal = reference->energy_calibration();
  const shared_ptr<const vector<float>> &ref_counts = reference->gamma_counts();
  const shared_ptr<const vector<float>> &test_counts = test->gamma_counts();
  if( !ref_cal || !ref_cal->valid() || !ref_counts
      || !test->energy_calibration() || !test->energy_calibration()->valid() || !test_counts )
    throw runtime_error( "chi2_for_gain_offset: invalid spectrum" );

  shared_ptr<SpecUtils::EnergyCalibration> trans_cal;
  vector<float> scaled;  //test counts rebinned onto the reference grid
  try
  {
    trans_cal = transform_calibration( test->energy_calibration(), gain, offset );
    SpecUtils::rebin_by_lower_edge( *trans_cal->channel_energies(), *test_counts,
                                    *ref_cal->channel_energies(), scaled );
  }catch( std::exception & )
  {
    return sm_invalid_chi2;
  }

  const vector<float> &ref_edges = *ref_cal->channel_energies();
  const bool no_upper = (upper_energy <= lower_energy);

  const size_t nchannel = std::min( ref_counts->size(), scaled.size() );

  // First pass: analytically solve the free normalization of the test spectrum.
  // Channels in the window that the transformed test doesnt cover have scaled[i]==0 (from the
  //  rebin), so they dont bias the normalization, but do penalize chi2 below - meaning the
  //  optimizer cant lower chi2 by shrinking the overlap.
  double sum_wrs = 0.0, sum_wss = 0.0;
  size_t nused = 0;
  for( size_t i = 0; i < nchannel; ++i )
  {
    const float e_lo = ref_edges[i];
    if( (e_lo < lower_energy) || (!no_upper && (e_lo >= upper_energy)) )
      continue;

    const double r = (*ref_counts)[i];
    const double s = scaled[i];
    const double w = 1.0 / std::max( r, 1.0 );
    sum_wrs += w * r * s;
    sum_wss += w * s * s;
    ++nused;
  }//for( first pass )

  if( (nused < 10) || (sum_wss <= 0.0) )
    return sm_invalid_chi2;

  const double norm = sum_wrs / sum_wss;

  double chi2 = 0.0;
  for( size_t i = 0; i < nchannel; ++i )
  {
    const float e_lo = ref_edges[i];
    if( (e_lo < lower_energy) || (!no_upper && (e_lo >= upper_energy)) )
      continue;

    const double r = (*ref_counts)[i];
    const double diff = r - norm * scaled[i];
    chi2 += diff * diff / std::max( r, 1.0 );
  }//for( second pass )

  num_channels_used = nused;

  return chi2;
}//chi2_for_gain_offset(...)


double refine_gain_offset( const std::shared_ptr<const SpecUtils::Measurement> &test,
                           const std::shared_ptr<const SpecUtils::Measurement> &reference,
                           const float lower_energy, const float upper_energy,
                           const bool fit_offset,
                           double &gain, double &offset )
{
  using boost::math::tools::brent_find_minima;

  if( (gain <= 0.0) || std::isnan(gain) || std::isinf(gain) )
    throw runtime_error( "refine_gain_offset: invalid starting gain" );

  if( !test || !reference || !reference->energy_calibration()
      || !reference->energy_calibration()->valid() )
    throw runtime_error( "refine_gain_offset: invalid spectrum" );

  // Fix the comparison window from the starting transform, so the set of compared channels
  //  cant change during optimization; channels a trial transform doesnt cover then count as
  //  zeros against the reference (see chi2_for_gain_offset), instead of dropping out.
  const shared_ptr<const SpecUtils::EnergyCalibration> ref_cal = reference->energy_calibration();
  const shared_ptr<const SpecUtils::EnergyCalibration> start_cal
                    = transform_calibration( test->energy_calibration(), gain, offset );

  const bool no_upper = (upper_energy <= lower_energy);
  const float eff_lower = std::max( lower_energy, start_cal->lower_energy() );
  float eff_upper = no_upper ? ref_cal->upper_energy()
                             : std::min( upper_energy, ref_cal->upper_energy() );
  eff_upper = std::min( eff_upper, start_cal->upper_energy() );

  if( eff_upper <= eff_lower )
    throw runtime_error( "refine_gain_offset: no overlap between spectra in the energy range" );

  const int bits = 24;

  auto chi2_fcn = [&]( const double g, const double b ) -> double {
    size_t nchan;
    return chi2_for_gain_offset( test, reference, eff_lower, eff_upper, g, b, nchan );
  };

  // Minimize in x = log(gain).  Sharp HPGe photopeaks make chi2 valleys as narrow as ~0.1%
  //  in gain, which a plain bracketed Brent can miss entirely, so localize the valley with a
  //  coarse-to-fine direct scan and only then polish within one fine step.
  double x = std::log( gain );

  auto fx = [&]( const double xv ) -> double { return chi2_fcn( std::exp(xv), offset ); };

  x = scan_minimum( fx, x, 1.25E-3, 24 );  // +-3.0% in 0.125% steps
  x = scan_minimum( fx, x, 1.25E-4, 10 );  // +-0.125% in 0.0125% steps
  {
    std::uintmax_t max_iter = 100;
    const pair<double,double> found = brent_find_minima( fx, x - 1.25E-4, x + 1.25E-4,
                                                         bits, max_iter );
    x = found.first;
  }

  if( fit_offset )
  {
    for( int round = 0; round < 5; ++round )
    {
      const double prev_x = x, prev_b = offset;

      {
        auto fb = [&]( const double bv ) -> double { return chi2_fcn( std::exp(x), bv ); };
        double b = scan_minimum( fb, offset, 2.0, 25 );  // +-50 keV in 2 keV steps
        b = scan_minimum( fb, b, 0.2, 10 );              // +-2 keV in 0.2 keV steps
        std::uintmax_t max_iter = 100;
        const pair<double,double> found = brent_find_minima( fb, b - 0.2, b + 0.2,
                                                             bits, max_iter );
        offset = found.first;
      }

      {
        x = scan_minimum( fx, x, 2.5E-4, 20 );  // re-polish gain within +-0.5%
        std::uintmax_t max_iter = 100;
        const pair<double,double> found = brent_find_minima( fx, x - 2.5E-4, x + 2.5E-4,
                                                             bits, max_iter );
        x = found.first;
      }

      if( (std::fabs(x - prev_x) < 1.0E-6) && (std::fabs(offset - prev_b) < 0.01) )
        break;
    }//for( alternate between offset and gain )
  }//if( fit_offset )

  gain = std::exp( x );

  size_t nchan = 0;
  const double chi2 = chi2_for_gain_offset( test, reference, eff_lower, eff_upper,
                                            gain, offset, nchan );
  const size_t npar = fit_offset ? 3 : 2;  //gain, normalization, and possibly offset

  return chi2 / static_cast<double>( std::max<size_t>( 1, (nchan > npar) ? (nchan - npar) : 1 ) );
}//refine_gain_offset(...)


double fit_peak_mean( const std::shared_ptr<const SpecUtils::Measurement> &spec,
                      const double expected_energy,
                      const std::shared_ptr<const PeakFitDetPrefs> &fit_prefs )
{
  if( !spec || !spec->energy_calibration() || !spec->energy_calibration()->valid()
      || !spec->gamma_counts() )
    throw runtime_error( "fit_peak_mean: invalid spectrum" );

  if( expected_energy <= 0.0 )
    throw runtime_error( "fit_peak_mean: invalid expected energy" );

  const bool highres = PeakFitUtils::is_high_res( spec );

  shared_ptr<const PeakFitDetPrefs> prefs = fit_prefs;
  if( !prefs )
  {
    auto defprefs = make_shared<PeakFitDetPrefs>();
    defprefs->m_det_type = highres ? PeakFitUtils::CoarseResolutionType::High
                                   : PeakFitUtils::CoarseResolutionType::Low;
    prefs = defprefs;
  }//if( !prefs )

  const vector<shared_ptr<const PeakDef>> no_peaks;
  const pair<PeakShrdVec,PeakShrdVec> found
             = searchForPeakFromUser( expected_energy, -1.0, spec, no_peaks, nullptr, nullptr, prefs );

  if( found.first.empty() )
    throw runtime_error( "fit_peak_mean: no peak found near "
                         + to_string(expected_energy) + " keV" );

  double mean = found.first[0]->mean();
  for( const shared_ptr<const PeakDef> &p : found.first )
  {
    if( std::fabs(p->mean() - expected_energy) < std::fabs(mean - expected_energy) )
      mean = p->mean();
  }

  const float energy = static_cast<float>( expected_energy );
  const float fwhm_est = highres ? PeakFitUtils::hpge_fwhm_fcn( energy )
                                 : PeakFitUtils::nai_fwhm_fcn( energy );
  const double tolerance = std::max( 2.0*fwhm_est, 0.03*expected_energy );

  if( std::fabs(mean - expected_energy) > tolerance )
    throw runtime_error( "fit_peak_mean: fitted peak at " + to_string(mean)
                         + " keV is too far from expected " + to_string(expected_energy) + " keV" );

  return mean;
}//fit_peak_mean(...)


MatchResults match( const std::vector<SpectrumInput> &inputs,
                    const int reference_index,
                    const MatchOptions &options )
{
  MatchResults answer;
  answer.results.resize( inputs.size() );

  const float lower = options.lower_energy;
  const float upper = options.upper_energy;

  // Screen out invalid or low-statistics spectra
  vector<double> counts( inputs.size(), 0.0 );
  size_t num_usable = 0;
  for( size_t i = 0; i < inputs.size(); ++i )
  {
    SpectrumResult &res = answer.results[i];
    const shared_ptr<const SpecUtils::Measurement> &spec = inputs[i].spectrum;
    const bool valid = spec && spec->energy_calibration() && spec->energy_calibration()->valid()
                       && spec->gamma_counts() && (spec->num_gamma_channels() >= 16);

    counts[i] = valid ? counts_in_range( spec, lower, upper ) : 0.0;
    res.used = (valid && (counts[i] >= sm_min_counts_in_range));
    if( res.used )
      ++num_usable;
    else
      res.warnings.push_back( WarningType::LowStatistics );
  }//for( screen inputs )

  if( num_usable < 2 )
    throw runtime_error( "Gain matching requires at least two usable spectra in the energy range" );

  // Choose the reference spectrum
  size_t ref_index = 0;
  if( reference_index >= 0 )
  {
    if( (static_cast<size_t>(reference_index) >= inputs.size())
        || !answer.results[reference_index].used )
      throw runtime_error( "The requested reference spectrum is not usable" );
    ref_index = static_cast<size_t>( reference_index );
  }else
  {
    double max_counts = -1.0;
    for( size_t i = 0; i < inputs.size(); ++i )
    {
      if( answer.results[i].used && (counts[i] > max_counts) )
      {
        max_counts = counts[i];
        ref_index = i;
      }
    }
  }//if( user specified reference ) / else

  answer.reference_index = ref_index;
  const shared_ptr<const SpecUtils::Measurement> &ref_spec = inputs[ref_index].spectrum;

  // The reference keeps its calibration; updated_cal left null
  answer.results[ref_index].gain = 1.0;
  answer.results[ref_index].offset = 0.0;

  // Fit the refinement peak in the reference first; a failure there disables stage 3 entirely
  bool do_peak = (options.ref_peak_energy > 0.0);
  if( do_peak )
  {
    try
    {
      answer.ref_peak_mean = fit_peak_mean( ref_spec, options.ref_peak_energy,
                                            inputs[ref_index].fit_prefs );
    }catch( std::exception & )
    {
      do_peak = false;
      answer.ref_peak_mean = 0.0;
      answer.results[ref_index].warnings.push_back( WarningType::PeakFitFailed );
    }
  }//if( do_peak )

  for( size_t i = 0; i < inputs.size(); ++i )
  {
    if( (i == ref_index) || !answer.results[i].used )
      continue;

    SpectrumResult &res = answer.results[i];
    const shared_ptr<const SpecUtils::Measurement> &spec = inputs[i].spectrum;

    try
    {
      // The log-energy cross-correlation assumes a pure gain difference; an unknown offset
      //  shifts low-energy features by a different log-amount than high-energy ones and can
      //  derail it.  So when offset fitting is enabled, try a small grid of offset
      //  hypotheses (pre-shifting the test spectrum by each) and keep the best-correlating.
      double gain = 1.0, offset = 0.0;

      if( options.fit_offset )
      {
        double best_corr = -2.0;
        for( int i = -5; i <= 5; ++i )
        {
          const double cand_offset = 12.0 * i;

          shared_ptr<const SpecUtils::Measurement> shifted = spec;
          if( i != 0 )
          {
            auto shifted_meas = make_shared<SpecUtils::Measurement>( *spec );
            shifted_meas->set_energy_calibration(
                    transform_calibration( spec->energy_calibration(), 1.0, cand_offset ) );
            shifted = shifted_meas;
          }

          try
          {
            double corr = 0.0;
            vector<WarningType> cand_warnings;
            const double cand_gain = coarse_gain_xcorr( shifted, ref_spec, lower, upper,
                                                        corr, cand_warnings );
            if( corr > best_corr )
            {
              best_corr = corr;
              res.correlation = corr;
              res.warnings = cand_warnings;
              gain = cand_gain;
              offset = cand_gain * cand_offset;  //composed: g*(E + b_cand) = g*E + g*b_cand
            }
          }catch( std::exception & )
          {
          }
        }//for( offset hypotheses )

        if( best_corr < -1.5 )
          throw runtime_error( "No offset hypothesis produced a usable correlation" );
      }else
      {
        gain = coarse_gain_xcorr( spec, ref_spec, lower, upper, res.correlation, res.warnings );
      }//if( options.fit_offset ) / else

      res.chi2_dof = refine_gain_offset( spec, ref_spec, lower, upper, options.fit_offset,
                                         gain, offset );

      if( do_peak )
      {
        // Fit the peak in a copy adjusted by the chi2 result, then scale the gain so the
        //  fitted mean lands exactly on the references fitted mean.
        try
        {
          auto trial = make_shared<SpecUtils::Measurement>( *spec );
          trial->set_energy_calibration( transform_calibration( spec->energy_calibration(),
                                                                gain, offset ) );

          const double mu = fit_peak_mean( trial, answer.ref_peak_mean, inputs[i].fit_prefs );
          res.fit_peak_mean = mu;

          if( ((mu - offset) > 1.0) && ((answer.ref_peak_mean - offset) > 1.0) )
            gain *= (answer.ref_peak_mean - offset) / (mu - offset);
          else
            res.warnings.push_back( WarningType::PeakFitFailed );
        }catch( std::exception & )
        {
          res.warnings.push_back( WarningType::PeakFitFailed );
        }
      }//if( do_peak )

      res.gain = gain;
      res.offset = offset;
      res.updated_cal = transform_calibration( spec->energy_calibration(), gain, offset );
    }catch( std::exception & )
    {
      res.used = false;
      res.warnings.push_back( WarningType::MatchFailed );
    }
  }//for( loop over non-reference inputs )

  return answer;
}//match(...)

}//namespace GainMatchCalc


namespace
{
  using SpecUtils::Measurement;
  using SpecUtils::SpectrumType;
  using SpecUtils::EnergyCalibration;

  const char *warning_txt_key( const GainMatchCalc::WarningType type )
  {
    switch( type )
    {
      case GainMatchCalc::WarningType::LowStatistics:        return "ecgm-warn-low-stats";
      case GainMatchCalc::WarningType::WeakCorrelation:      return "ecgm-warn-weak-corr";
      case GainMatchCalc::WarningType::CorrelationAtEdge:    return "ecgm-warn-corr-edge";
      case GainMatchCalc::WarningType::PeakFitFailed:        return "ecgm-warn-peak-fit-failed";
      case GainMatchCalc::WarningType::AmbiguousCorrelation: return "ecgm-warn-ambiguous";
      case GainMatchCalc::WarningType::MatchFailed:          return "ecgm-warn-match-failed";
      case GainMatchCalc::WarningType::NumWarningType:       break;
    }
    return "";
  }//warning_txt_key(...)


  /** Line color for the i'th spectrum, from the color themes reference-line colors. */
  WColor row_color( InterSpec *viewer, const size_t index )
  {
    static const char * const s_default_colors[] = {
      "#0000ff", "#f44336", "#4caf50", "#9c27b0", "#ff9800", "#00bcd4", "#795548", "#607d8b"
    };

    vector<WColor> colors;
    const shared_ptr<const ColorTheme> theme = viewer ? viewer->getColorTheme() : nullptr;
    if( theme )
      colors = theme->referenceLineColor;
    if( colors.empty() )
    {
      for( const char *color : s_default_colors )
        colors.emplace_back( color );
    }

    return colors[index % colors.size()];
  }//row_color(...)


  /** Names of detectors with gamma data in the displayed samples of the given file
   (mirrors EnergyCalTool::gammaDetectorsForDisplayedSamples).
   */
  vector<string> gamma_detectors( InterSpec *viewer, const SpectrumType type )
  {
    vector<string> answer;

    const shared_ptr<const SpecMeas> meas = viewer->measurment( type );
    if( !meas )
      return answer;

    const set<int> &samples = viewer->displayedSamples( type );

    for( const string &name : meas->gamma_detector_names() )
    {
      for( const int sample : samples )
      {
        const shared_ptr<const Measurement> m = meas->measurement( sample, name );
        if( m && (m->num_gamma_channels() > 4)
            && m->energy_calibration() && m->energy_calibration()->valid() )
        {
          answer.push_back( name );
          break;
        }
      }
    }//for( loop over detector names )

    return answer;
  }//gamma_detectors(...)


  // For undo/redo: per-Measurement old/new calibrations, and per-sample-set old/new peaks
  //  (same shapes as the anonymous-namespace types in EnergyCalTool.cpp).
  typedef vector<tuple<weak_ptr<const Measurement>,
                       shared_ptr<const EnergyCalibration>,
                       shared_ptr<const EnergyCalibration>>> meas_old_new_cal_t;

  typedef vector<tuple<set<int>,
                       deque<shared_ptr<const PeakDef>>,
                       deque<shared_ptr<const PeakDef>>>> meas_old_new_peaks_t;

  struct GainMatchFileChange
  {
    weak_ptr<SpecMeas> file;
    SpectrumType type;
    meas_old_new_cal_t cals;
    meas_old_new_peaks_t peaks;
    meas_old_new_peaks_t hint_peaks;

    /** If true, after (re)applying the cal change, the files stale automated-search
     ("hint") peaks are invalidated and a fresh search is relaunched for the displayed
     spectrum.  Used for within-file matching, where the individual detectors realign onto
     the (unchanged) reference energy scale and so the summed data changes shape, but the
     fit peaks are left in place (see applyChanges).
     */
    bool relaunch_hints = false;
  };//struct GainMatchFileChange


  /** Applies (is_undo==false), or reverts (is_undo==true), the gain-match calibration and
   peak changes; used both for the initial commit and by the undo/redo step.
   Follows the same file-still-displayed guards as the energy-cal undo/redo machinery
   (see do_undo_or_redo in EnergyCalTool.cpp).
   */
  void apply_gain_match_direction( const bool is_undo, const vector<GainMatchFileChange> &changes )
  {
    InterSpec * const viewer = InterSpec::instance();
    assert( viewer );
    if( !viewer )
      return;

    const shared_ptr<SpecMeas> foreground = viewer->measurment( SpectrumType::Foreground );
    const set<int> &foreSamples = viewer->displayedSamples( SpectrumType::Foreground );
    const shared_ptr<SpecMeas> background = viewer->measurment( SpectrumType::Background );
    const set<int> &backSamples = viewer->displayedSamples( SpectrumType::Background );
    const shared_ptr<SpecMeas> secondary = viewer->measurment( SpectrumType::SecondForeground );
    const set<int> &secoSamples = viewer->displayedSamples( SpectrumType::SecondForeground );

    const shared_ptr<PeakModel> peakModel = viewer->peakModel();

    for( const GainMatchFileChange &change : changes )
    {
      const shared_ptr<SpecMeas> specfile = change.file.lock();
      const shared_ptr<SpecMeas> specfile_now = viewer->measurment( change.type );

      if( !specfile || (specfile != specfile_now) )
      {
        viewer->logMessage( WString::tr("ecgm-file-changed"), 2 );
        continue;
      }

      for( const auto &m_o_n : change.cals )
      {
        const shared_ptr<const Measurement> m = get<0>(m_o_n).lock();
        const shared_ptr<const EnergyCalibration> &to_cal = is_undo ? get<1>(m_o_n) : get<2>(m_o_n);
        assert( to_cal );
        if( !m || !to_cal )
          continue;

        try
        {
          specfile->set_energy_calibration( to_cal, m );
        }catch( std::exception & )
        {
          // Measurement no longer in file (e.g., file was modified) - nothing sensible to do
        }
      }//for( loop over calibration changes )

      for( const auto &m_o_n : change.peaks )
      {
        const set<int> &samples = get<0>(m_o_n);
        const deque<shared_ptr<const PeakDef>> &to_peaks = is_undo ? get<1>(m_o_n) : get<2>(m_o_n);

        specfile->setPeaks( to_peaks, samples );
        if( peakModel && (specfile == foreground) && (samples == foreSamples) )
          peakModel->setPeakFromSpecMeas( foreground, foreSamples, SpectrumType::Foreground );
        else if( peakModel && (specfile == background) && (samples == backSamples) )
          peakModel->setPeakFromSpecMeas( background, backSamples, SpectrumType::Background );
        else if( peakModel && (specfile == secondary) && (samples == secoSamples) )
          peakModel->setPeakFromSpecMeas( secondary, secoSamples, SpectrumType::SecondForeground );
      }//for( loop over peak changes )

      for( const auto &m_o_n : change.hint_peaks )
      {
        const set<int> &samples = get<0>(m_o_n);
        const deque<shared_ptr<const PeakDef>> &to_peaks = is_undo ? get<1>(m_o_n) : get<2>(m_o_n);
        auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( to_peaks );
        specfile->setAutomatedSearchPeaks( samples, peaks );
      }//for( loop over hint peak changes )
    }//for( loop over changed files )

    viewer->refreshDisplayedCharts();
    EnergyCalTool * const tool = viewer->energyCalTool();
    if( tool )
      tool->refreshGuiFromFiles();

    // Relaunch the automated-search ("hint") peaks for any displayed spectrum whose data
    //  changed shape - done after refreshDisplayedCharts() so the search runs on the new
    //  summed data.  Applies on the initial commit as well as on undo/redo.
    for( const GainMatchFileChange &change : changes )
    {
      if( !change.relaunch_hints )
        continue;

      const shared_ptr<SpecMeas> specfile = change.file.lock();
      if( !specfile || (specfile != viewer->measurment(change.type)) )
        continue;

      const set<int> disp_samples = viewer->displayedSamples( change.type );
      const shared_ptr<const SpecUtils::Measurement> disp_hist
                                          = viewer->displayedHistogram( change.type );
      if( disp_hist )
      {
        specfile->setAutomatedSearchPeaks( disp_samples,
                                     make_shared<deque<shared_ptr<const PeakDef>>>() );
        viewer->searchForHintPeaks( specfile, disp_samples, disp_hist,
                                    specfile->peakFitDetPrefs() );
      }
    }//for( relaunch hint peaks where requested )
  }//apply_gain_match_direction(...)
}//namespace


namespace GainMatchCalc
{

size_t applyDetectorGains( InterSpec *interspec,
                           const SpecUtils::SpectrumType type,
                           const std::shared_ptr<SpecMeas> &meas,
                           const std::set<int> &apply_samples,
                           const std::vector<std::tuple<std::string,double,double>> &per_detector )
{
  if( !interspec || !meas || per_detector.empty() )
    return 0;

  GainMatchFileChange change;
  change.file = meas;
  change.type = type;
  change.relaunch_hints = true;

  size_t num_updated = 0;
  for( const auto &det_gain_offset : per_detector )
  {
    const string &det = get<0>( det_gain_offset );
    const double gain = get<1>( det_gain_offset );
    const double offset = get<2>( det_gain_offset );

    bool changed_this_det = false;
    for( const int sample : apply_samples )
    {
      const shared_ptr<const Measurement> m = meas->measurement( sample, det );
      if( !m || (m->num_gamma_channels() < 5) )
        continue;

      const shared_ptr<const EnergyCalibration> old_cal = m->energy_calibration();
      if( !old_cal || !old_cal->valid() )
        continue;

      try
      {
        const shared_ptr<const EnergyCalibration> new_cal
                                      = transform_calibration( old_cal, gain, offset );
        change.cals.emplace_back( m, old_cal, new_cal );
        changed_this_det = true;
      }catch( std::exception & )
      {
      }
    }//for( loop over samples )

    if( changed_this_det )
      num_updated += 1;
  }//for( loop over detectors )

  if( change.cals.empty() )
    return 0;

  auto changes = make_shared<vector<GainMatchFileChange>>();
  changes->push_back( std::move(change) );

  apply_gain_match_direction( false, *changes );

  UndoRedoManager * const undoManager = interspec->undoRedoManager();
  if( undoManager )
  {
    const shared_ptr<const vector<GainMatchFileChange>> applied = changes;
    auto undo = [applied](){ apply_gain_match_direction( true, *applied ); };
    auto redo = [applied](){ apply_gain_match_direction( false, *applied ); };
    undoManager->addUndoRedoStep( undo, redo, "Gain-match energy calibration" );
  }//if( undoManager )

  return num_updated;
}//applyDetectorGains(...)


std::vector<SpectrumInput>
buildDetectorInputs( const std::shared_ptr<SpecMeas> &meas,
                     const std::set<int> &displayed_samples )
{
  vector<SpectrumInput> inputs;

  try
  {
    if( !meas || displayed_samples.empty() )
      return inputs;

    // Cost gate: keep the analysis quick.  Skip large passthrough/search files where summing
    //  many samples would be slow, or systems with an unusual number of detectors.
    if( displayed_samples.size() > 64 )
      return inputs;

    // Gather gamma detectors that have data in the displayed samples
    vector<pair<string,shared_ptr<const EnergyCalibration>>> dets;
    for( const string &name : meas->gamma_detector_names() )
    {
      for( const int sample : displayed_samples )
      {
        const shared_ptr<const Measurement> m = meas->measurement( sample, name );
        if( m && (m->num_gamma_channels() > 4)
            && m->energy_calibration() && m->energy_calibration()->valid() )
        {
          dets.emplace_back( name, m->energy_calibration() );
          break;
        }
      }
    }//for( loop over detector names )

    if( (dets.size() < 2) || (dets.size() > 32)
        || ((dets.size() * displayed_samples.size()) > 512) )
      return inputs;

    for( const auto &name_cal : dets )
    {
      SpectrumInput input;
      input.name = name_cal.first;
      input.fit_prefs = meas->peakFitDetPrefs();
      try
      {
        // sum_measurements makes fresh, immutable Measurement objects, so the returned inputs
        //  are safe to hand to a worker thread (no further SpecMeas access needed).
        input.spectrum = meas->sum_measurements( displayed_samples, {name_cal.first},
                                                 name_cal.second );
      }catch( std::exception & )
      {
      }
      if( input.spectrum && (input.spectrum->num_gamma_channels() <= 65536) )
        inputs.push_back( input );
    }//for( build one input per detector )

    if( inputs.size() < 2 )
      inputs.clear();
  }catch( std::exception & )
  {
    inputs.clear();
  }

  return inputs;
}//buildDetectorInputs(...)


DetectorMatchResult
matchDetectorInputs( const std::vector<SpectrumInput> &inputs, const double min_shift_kev )
{
  DetectorMatchResult answer;

  try
  {
    if( inputs.size() < 2 )
      return answer;

    MatchOptions options;
    options.lower_energy = 200.0f;
    options.upper_energy = 0.0f;   //no upper limit
    options.fit_offset = false;    //gain only, for speed and robustness on load

    const MatchResults results = match( inputs, -1, options );
    if( results.reference_index >= inputs.size() )
      return answer;

    // Require every detector to have matched confidently - if any is unreliable, dont
    //  silently suggest a change; the user can still open the tool by hand.
    double min_gain = 1.0, max_gain = 1.0;
    for( size_t i = 0; i < results.results.size(); ++i )
    {
      const SpectrumResult &res = results.results[i];
      if( !res.used )
        return answer;

      for( const WarningType w : res.warnings )
      {
        if( (w == WarningType::WeakCorrelation) || (w == WarningType::CorrelationAtEdge)
            || (w == WarningType::AmbiguousCorrelation) || (w == WarningType::MatchFailed) )
          return answer;
      }

      if( i != results.reference_index )
      {
        if( !res.updated_cal )
          return answer;
        min_gain = std::min( min_gain, res.gain );
        max_gain = std::max( max_gain, res.gain );
      }
    }//for( check every detector matched confidently )

    // Consistency measure: compare the largest implied energy discrepancy between detectors to
    //  the intrinsic width of a real peak, so we only suggest matching when it would visibly
    //  sharpen the summed data.  A detector-type-agnostic estimate of the resolution is taken
    //  by fitting the tallest peak in the reference spectrum (channel-count heuristics like
    //  is_high_res misclassify high-channel-count NaI portals, so we measure the data instead).
    const shared_ptr<const Measurement> ref_spec = inputs[results.reference_index].spectrum;
    const shared_ptr<const EnergyCalibration> ref_cal = ref_spec->energy_calibration();
    const shared_ptr<const vector<float>> ref_counts = ref_spec->gamma_counts();
    const vector<float> &ref_edges = *ref_cal->channel_energies();

    // Energy of the tallest channel above the lower energy - the seed for the resolution fit
    double tallest_energy = 0.0, tallest_counts = -1.0;
    for( size_t i = 0; (i < ref_counts->size()) && ((i+1) < ref_edges.size()); ++i )
    {
      if( (ref_edges[i] >= options.lower_energy) && ((*ref_counts)[i] > tallest_counts) )
      {
        tallest_counts = (*ref_counts)[i];
        tallest_energy = 0.5*( ref_edges[i] + ref_edges[i+1] );
      }
    }

    const double gain_spread = max_gain - min_gain;

    // Try to measure the peak width; the implied inter-detector spread at that peak is compared
    //  to it.  A spread of ~1/3 of the FWHM corresponds to roughly a 5% broadening of the
    //  summed peak (added in quadrature) - the users "shrinks FWHM by >5%" criterion.
    double eval_energy = tallest_energy;
    double resolution_floor = 0.0;
    if( tallest_energy > 0.0 )
    {
      try
      {
        const shared_ptr<const PeakFitDetPrefs> prefs = inputs[results.reference_index].fit_prefs;
        const vector<shared_ptr<const PeakDef>> no_peaks;
        const pair<PeakShrdVec,PeakShrdVec> found
              = searchForPeakFromUser( tallest_energy, -1.0, ref_spec, no_peaks, nullptr, nullptr, prefs );
        if( !found.first.empty() )
        {
          shared_ptr<const PeakDef> peak = found.first.front();
          for( const shared_ptr<const PeakDef> &p : found.first )
          {
            if( std::fabs(p->mean() - tallest_energy) < std::fabs(peak->mean() - tallest_energy) )
              peak = p;
          }
          eval_energy = peak->mean();
          resolution_floor = 0.33 * peak->fwhm();
        }
      }catch( std::exception & )
      {
      }
    }//if( tallest_energy > 0.0 )

    const double max_shift = gain_spread * eval_energy;

    // Flag only when the implied spread exceeds both an absolute floor (mostly relevant for
    //  high-resolution detectors) and the resolution-based floor (quiet on low-resolution
    //  detectors whose broad peaks tolerate a few keV of misalignment).
    if( (max_shift <= min_shift_kev) || (max_shift <= resolution_floor) )
      return answer;

    answer.max_shift_kev = max_shift;
    for( size_t i = 0; i < inputs.size(); ++i )
    {
      if( i != results.reference_index )
        answer.per_detector.emplace_back( inputs[i].name, results.results[i].gain,
                                          results.results[i].offset );
    }
    answer.beneficial = !answer.per_detector.empty();
  }catch( std::exception & )
  {
    answer = DetectorMatchResult();
  }

  return answer;
}//matchDetectorInputs(...)


std::shared_ptr<DetectorMatchSuggestion>
analyzeForAutoMatch( const SpecUtils::SpectrumType type,
                     const std::shared_ptr<SpecMeas> &meas,
                     const std::set<int> &displayed_samples,
                     const double min_shift_kev )
{
  const vector<SpectrumInput> inputs = buildDetectorInputs( meas, displayed_samples );
  const DetectorMatchResult result = matchDetectorInputs( inputs, min_shift_kev );
  if( !result.beneficial )
    return nullptr;

  auto suggestion = make_shared<DetectorMatchSuggestion>();
  suggestion->file = meas;
  suggestion->type = type;
  suggestion->displayed_samples = displayed_samples;
  suggestion->max_shift_kev = result.max_shift_kev;
  suggestion->per_detector = result.per_detector;

  return suggestion;
}//analyzeForAutoMatch(...)


std::vector<SharedPeak>
findSharedPeaks( const std::vector<SpectrumInput> &inputs,
                 const MatchResults &stage2,
                 const std::deque<std::shared_ptr<const PeakDef>> &id_peaks,
                 const float lower_energy, const float upper_energy )
{
  vector<SharedPeak> answer;

  try
  {
    const size_t ndet = inputs.size();
    if( (ndet < 2) || (stage2.results.size() != ndet) )
      return answer;

    const bool no_upper = (upper_energy <= lower_energy);
    const size_t ref = stage2.reference_index;

    // A fast candidate-peak search per detector, on its stage-2-aligned spectrum (so candidates
    //  are in a common energy scale).  Keep only the strongest few candidates per detector - the
    //  raw finder returns many noise bumps, which would otherwise cluster into spurious "peaks".
    struct Candidate{ double energy; double channel; double area; };
    vector<vector<Candidate>> det_cands( ndet );
    const size_t max_per_det = 12;

    for( size_t d = 0; d < ndet; ++d )
    {
      const SpectrumResult &res = stage2.results[d];
      const shared_ptr<const Measurement> spec = inputs[d].spectrum;
      if( !res.used || !spec || !spec->energy_calibration() || !spec->energy_calibration()->valid() )
        continue;

      shared_ptr<const EnergyCalibration> adj_cal;
      if( d == ref )
        adj_cal = spec->energy_calibration();
      else
        adj_cal = transform_calibration( spec->energy_calibration(), res.gain, res.offset );

      auto adjusted = make_shared<Measurement>( *spec );
      adjusted->set_energy_calibration( adj_cal );

      vector<tuple<float,float,float>> cands;  //(mean, sigma, area)
      secondDerivativePeakCanidates( adjusted, inputs[d].fit_prefs, 0, 0, cands );

      vector<Candidate> keep;
      for( const tuple<float,float,float> &c : cands )
      {
        const double mean = get<0>( c ), area = get<2>( c );
        if( (mean < lower_energy) || (!no_upper && (mean > upper_energy)) || (area <= 0.0) )
          continue;
        // channel is calibration-invariant (adj_cal is a linear transform of the original)
        keep.push_back( Candidate{ mean, adj_cal->channel_for_energy(mean), area } );
      }

      // Keep the strongest few by area, then order by energy for the matching below.
      std::sort( begin(keep), end(keep),
                 [](const Candidate &a, const Candidate &b){ return a.area > b.area; } );
      if( keep.size() > max_per_det )
        keep.resize( max_per_det );
      std::sort( begin(keep), end(keep),
                 [](const Candidate &a, const Candidate &b){ return a.energy < b.energy; } );
      det_cands[d] = keep;
    }//for( per-detector candidate search )

    const size_t min_dets = std::max<size_t>( 2, (ndet + 1)/2 );  //present in >= half the detectors

    // Seed clusters from the reference detector's peaks (so every kept peak has a reference
    //  target), then find the nearest candidate in each other detector.
    double last_energy = -1.0e9;
    for( const Candidate &seed : det_cands[ref] )
    {
      const double window = std::max( 0.01*seed.energy, 2.0 );  //~1% of energy, at least 2 keV
      if( (seed.energy - last_energy) < window )
        continue;  //too close to the previous kept peak - avoid duplicates/overlap

      SharedPeak peak;
      peak.detector_channels.assign( ndet, -1.0 );
      peak.detector_channels[ref] = seed.channel;
      int found = 1;

      for( size_t d = 0; d < ndet; ++d )
      {
        if( (d == ref) || !stage2.results[d].used )
          continue;

        double best_diff = window, best_channel = -1.0;
        for( const Candidate &c : det_cands[d] )
        {
          const double diff = std::fabs( c.energy - seed.energy );
          if( diff < best_diff )
          {
            best_diff = diff;
            best_channel = c.channel;
          }
        }
        if( best_channel >= 0.0 )
        {
          peak.detector_channels[d] = best_channel;
          found += 1;
        }
      }//for( find this peak in each other detector )

      if( static_cast<size_t>(found) < min_dets )
        continue;

      last_energy = seed.energy;
      peak.num_detectors = found;
      peak.energy = seed.energy;
      peak.target_energy = seed.energy;

      // Mark identified if a nuclide-assigned peak sits at this energy
      string nuclide;
      for( const shared_ptr<const PeakDef> &p : id_peaks )
      {
        if( !p || !p->gausPeak() || !p->hasSourceGammaAssigned() )
          continue;
        if( std::fabs(p->mean() - peak.energy) <= window )
        {
          const double gamma_energy = p->gammaParticleEnergy();
          if( gamma_energy > 0.0 )
          {
            peak.identified = true;
            peak.target_energy = gamma_energy;
            if( p->parentNuclide() )
              nuclide = p->parentNuclide()->symbol;
          }
          break;
        }
      }//for( match against identified peaks )

      char buf[64] = { '\0' };
      snprintf( buf, sizeof(buf), "%.1f keV", peak.target_energy );
      peak.label = nuclide.empty() ? string(buf) : (nuclide + " " + buf);

      answer.push_back( peak );
    }//for( reference-seeded clustering )
  }catch( std::exception & )
  {
    answer.clear();
  }

  return answer;
}//findSharedPeaks(...)


std::vector<std::shared_ptr<const SpecUtils::EnergyCalibration>>
refineWithSharedPeaks( const std::vector<SpectrumInput> &inputs,
                       const MatchResults &stage2,
                       const std::vector<SharedPeak> &peaks,
                       const std::vector<bool> &use_peak,
                       const int order )
{
  const size_t ndet = inputs.size();
  vector<shared_ptr<const EnergyCalibration>> answer( ndet, nullptr );

  try
  {
    if( stage2.results.size() != ndet )
      return answer;

    const int ncoef = std::max( 2, std::min(4, order + 1) );  //order 1 -> 2 coefs (offset+gain)

    bool any_identified = false;
    for( size_t j = 0; (j < peaks.size()) && (j < use_peak.size()); ++j )
      any_identified = (any_identified || (use_peak[j] && peaks[j].identified));

    for( size_t d = 0; d < ndet; ++d )
    {
      if( !stage2.results[d].used )
        continue;

      // In relative mode (no identified peaks) the reference defines the target energies, so it
      //  is left unchanged; with identified peaks it is re-fit toward the known energies too.
      const bool is_reference = (d == stage2.reference_index);
      if( is_reference && !any_identified )
        continue;

      const shared_ptr<const Measurement> spec = inputs[d].spectrum;
      const shared_ptr<const EnergyCalibration> orig = spec ? spec->energy_calibration() : nullptr;
      if( !orig || !orig->valid() )
        continue;

      const size_t nchannel = orig->num_channels();

      vector<EnergyCal::RecalPeakInfo> recalinfos;
      for( size_t j = 0; (j < peaks.size()) && (j < use_peak.size()); ++j )
      {
        if( !use_peak[j] || (d >= peaks[j].detector_channels.size()) )
          continue;
        const double channel = peaks[j].detector_channels[d];
        if( channel < 0.0 )
          continue;

        EnergyCal::RecalPeakInfo info;
        info.peakMeanBinNumber = channel;
        info.peakMean = orig->energy_for_channel( channel );
        info.peakMeanUncert = 1.0;   //equal weighting
        info.photopeakEnergy = peaks[j].target_energy;
        recalinfos.push_back( info );
      }//for( gather this detectors selected peaks )

      if( static_cast<int>(recalinfos.size()) < ncoef )
        continue;  //not enough peaks located in this detector for the requested order

      const vector<bool> fitfor( ncoef, true );
      vector<float> coefs( ncoef, 0.0f ), coefs_uncert;

      try
      {
        EnergyCal::fit_energy_cal_poly( recalinfos, fitfor, nchannel, orig->deviation_pairs(),
                                        coefs, coefs_uncert );
        auto newcal = make_shared<EnergyCalibration>();
        newcal->set_polynomial( nchannel, coefs, orig->deviation_pairs() );
        answer[d] = newcal;
      }catch( std::exception & )
      {
      }
    }//for( loop over detectors )
  }catch( std::exception & )
  {
    answer.assign( ndet, nullptr );
  }

  return answer;
}//refineWithSharedPeaks(...)


size_t applyDetectorCals( InterSpec *interspec,
                          const SpecUtils::SpectrumType type,
                          const std::shared_ptr<SpecMeas> &meas,
                          const std::set<int> &apply_samples,
                          const std::vector<std::pair<std::string,
                                 std::shared_ptr<const SpecUtils::EnergyCalibration>>> &per_detector )
{
  if( !interspec || !meas || per_detector.empty() )
    return 0;

  GainMatchFileChange change;
  change.file = meas;
  change.type = type;
  change.relaunch_hints = true;

  size_t num_updated = 0;
  for( const auto &name_cal : per_detector )
  {
    const string &det = name_cal.first;
    const shared_ptr<const EnergyCalibration> newcal = name_cal.second;
    if( !newcal || !newcal->valid() )
      continue;

    bool changed_this_det = false;
    for( const int sample : apply_samples )
    {
      const shared_ptr<const Measurement> m = meas->measurement( sample, det );
      if( !m || (m->num_gamma_channels() < 5) )
        continue;

      const shared_ptr<const EnergyCalibration> old_cal = m->energy_calibration();
      if( !old_cal || !old_cal->valid() )
        continue;

      // The fitted calibration is a polynomial on this detectors channels; apply it directly
      //  when the channel count matches (the common case - a detectors measurements share the
      //  same channelization).  Skip any measurement whose channel count differs.
      if( m->num_gamma_channels() != newcal->num_channels() )
        continue;

      change.cals.emplace_back( m, old_cal, newcal );
      changed_this_det = true;
    }//for( loop over samples )

    if( changed_this_det )
      num_updated += 1;
  }//for( loop over detectors )

  if( change.cals.empty() )
    return 0;

  auto changes = make_shared<vector<GainMatchFileChange>>();
  changes->push_back( std::move(change) );

  apply_gain_match_direction( false, *changes );

  UndoRedoManager * const undoManager = interspec->undoRedoManager();
  if( undoManager )
  {
    const shared_ptr<const vector<GainMatchFileChange>> applied = changes;
    auto undo = [applied](){ apply_gain_match_direction( true, *applied ); };
    auto redo = [applied](){ apply_gain_match_direction( false, *applied ); };
    undoManager->addUndoRedoStep( undo, redo, "Gain-match energy calibration" );
  }//if( undoManager )

  return num_updated;
}//applyDetectorCals(...)

}//namespace GainMatchCalc


EnergyCalGainMatch::EnergyCalGainMatch( EnergyCalTool *cal, AuxWindow *parent )
  : WContainerWidget(),
    m_interspec( InterSpec::instance() ),
    m_calibrator( cal ),
    m_parent( parent ),
    m_modeGroup( nullptr ),
    m_withinFileBtn( nullptr ),
    m_betweenFilesBtn( nullptr ),
    m_fileSelect( nullptr ),
    m_fileSelectTypes{},
    m_lowerEnergy( nullptr ),
    m_upperEnergy( nullptr ),
    m_fitOffset( nullptr ),
    m_refineWithPeaks( nullptr ),
    m_fitOrderCombo( nullptr ),
    m_peaksRow( nullptr ),
    m_peaksTable( nullptr ),
    m_peaksStatus( nullptr ),
    m_allSamples( nullptr ),
    m_showOriginal( nullptr ),
    m_rowTable( nullptr ),
    m_refGroup( nullptr ),
    m_rows{},
    m_lastInputs{},
    m_lastInputRows{},
    m_lastStage2{},
    m_sharedPeaks{},
    m_peakUseCbs{},
    m_sharedPeaksValid( false ),
    m_findingSharedPeaks( false ),
    m_refineGeneration( 0 ),
    m_chart( nullptr ),
    m_rangeHighlightId( 0 ),
    m_chartXRangeSet( false ),
    m_status( nullptr ),
    m_cancel( nullptr ),
    m_use( nullptr ),
    m_renderFlags(),
    m_withinFileMeas( nullptr ),
    m_betweenFileMeas{}
{
  assert( m_interspec && cal && parent );

  wApp->useStyleSheet( "InterSpec_resources/EnergyCalGainMatch.css" );
  m_interspec->useMessageResourceBundle( "EnergyCalGainMatch" );

  addStyleClass( "EnergyCalGainMatch" );

  m_chart = addNew<D3SpectrumDisplayDiv>();
  m_chart->addStyleClass( "GainMatchChart" );
  m_chart->setCompactAxis( true );
  m_chart->setYAxisLog( true );
  m_chart->enableLegend();
  m_chart->showYAxisScalers( false );

  WContainerWidget * const controls = addNew<WContainerWidget>();
  controls->addStyleClass( "GainMatchControls" );

  // Mode row: "Match detectors within [file]"  /  "Match between displayed spectra"
  WContainerWidget * const modeRow = controls->addNew<WContainerWidget>();
  modeRow->addStyleClass( "GainMatchRow" );

  m_modeGroup = std::make_shared<WButtonGroup>();
  m_withinFileBtn = modeRow->addNew<WRadioButton>( WString::tr("ecgm-mode-within") );
  m_modeGroup->addButton( m_withinFileBtn );
  m_fileSelect = modeRow->addNew<WComboBox>();
  m_betweenFilesBtn = modeRow->addNew<WRadioButton>( WString::tr("ecgm-mode-between") );
  m_modeGroup->addButton( m_betweenFilesBtn );

  const SpectrumType types[3] = { SpectrumType::Foreground, SpectrumType::Background,
                                  SpectrumType::SecondForeground };
  const char * const type_keys[3] = { "Foreground", "Background", "Secondary" };

  for( int i = 0; i < 3; ++i )
  {
    if( m_interspec->measurment(types[i]) && (gamma_detectors(m_interspec, types[i]).size() > 1) )
    {
      m_fileSelect->addItem( WString::tr(type_keys[i]) );
      m_fileSelectTypes.push_back( types[i] );
    }
  }//for( loop over spectrum types )

  set<shared_ptr<SpecMeas>> distinct_files;
  for( int i = 0; i < 3; ++i )
  {
    const shared_ptr<SpecMeas> meas = m_interspec->measurment( types[i] );
    if( meas )
      distinct_files.insert( meas );
  }

  const bool can_within = !m_fileSelectTypes.empty();
  const bool can_between = (distinct_files.size() > 1);
  m_withinFileBtn->setEnabled( can_within );
  m_betweenFilesBtn->setEnabled( can_between );
  m_fileSelect->setEnabled( can_within );
  m_modeGroup->setCheckedButton( can_within ? m_withinFileBtn : m_betweenFilesBtn );

  m_modeGroup->checkedChanged().connect( this, [this]( WRadioButton * ){ handleModeChange(); } );
  m_fileSelect->activated().connect( this, [this]( int ){ handleModeChange(); } );

  // Options row: energy range, offset, refinement peak
  WContainerWidget * const optionsRow = controls->addNew<WContainerWidget>();
  optionsRow->addStyleClass( "GainMatchRow" );

  optionsRow->addNew<WLabel>( WString::tr("ecgm-energy-range") );
  m_lowerEnergy = optionsRow->addNew<NativeFloatSpinBox>();
  m_lowerEnergy->setSpinnerHidden( true );
  m_lowerEnergy->setValue( 200.0f );
  optionsRow->addNew<WLabel>( WString::tr("ecgm-to") );
  m_upperEnergy = optionsRow->addNew<NativeFloatSpinBox>();
  m_upperEnergy->setSpinnerHidden( true );
  m_upperEnergy->setText( "" );
  m_upperEnergy->setPlaceholderText( WString::tr("ecgm-no-limit") );
  optionsRow->addNew<WLabel>( WString::tr("ecgm-kev") );

  m_fitOffset = optionsRow->addNew<WCheckBox>( WString::tr("ecgm-fit-offset") );

  // Multi-peak refinement row (WithinFile mode only): fit higher-order coefficients per detector
  //  from peaks common to all detectors.
  m_peaksRow = controls->addNew<WContainerWidget>();
  m_peaksRow->addStyleClass( "GainMatchRow GainMatchPeaksRow" );
  m_refineWithPeaks = m_peaksRow->addNew<WCheckBox>( WString::tr("ecgm-refine-peaks") );
  m_peaksRow->addNew<WLabel>( WString::tr("ecgm-fit-order") );
  m_fitOrderCombo = m_peaksRow->addNew<WComboBox>();
  m_fitOrderCombo->addItem( WString::tr("ecgm-order-linear") );     //offset + gain
  m_fitOrderCombo->addItem( WString::tr("ecgm-order-quadratic") );  //+ quadratic
  m_fitOrderCombo->addItem( WString::tr("ecgm-order-cubic") );      //+ cubic
  m_fitOrderCombo->setCurrentIndex( 1 );                           //default quadratic
  m_fitOrderCombo->disable();
  m_peaksStatus = m_peaksRow->addNew<WText>();
  m_peaksStatus->addStyleClass( "GainMatchPeaksStatus" );

  m_peaksTable = controls->addNew<WContainerWidget>();
  m_peaksTable->addStyleClass( "GainMatchPeaksTable" );
  m_peaksTable->setHidden( true );

  // Second options row: sample scope and preview toggle
  WContainerWidget * const scopeRow = controls->addNew<WContainerWidget>();
  scopeRow->addStyleClass( "GainMatchRow" );
  m_allSamples = scopeRow->addNew<WCheckBox>( WString::tr("ecgm-apply-all-samples") );
  m_allSamples->setChecked( true );
  m_showOriginal = scopeRow->addNew<WCheckBox>( WString::tr("ecgm-show-original") );

  // Table of spectra
  m_rowTable = controls->addNew<WContainerWidget>();
  m_rowTable->addStyleClass( "GainMatchTable HideOffset" );

  m_status = controls->addNew<WText>();
  m_status->addStyleClass( "GainMatchStatus" );

  // Signal wiring
  m_lowerEnergy->valueChanged().connect( this, [this]( float ){ handleOptionsChange(); } );
  m_upperEnergy->valueChanged().connect( this, [this]( float ){ handleOptionsChange(); } );
  m_fitOffset->checked().connect( this, [this](){
    m_rowTable->removeStyleClass( "HideOffset" );
    handleOptionsChange();
  } );
  m_fitOffset->unChecked().connect( this, [this](){
    m_rowTable->addStyleClass( "HideOffset" );
    handleOptionsChange();
  } );
  // Toggling refinement needs the shared-peak search (an UpdateCalc); changing the order only
  //  re-fits from the cached peaks (an UpdateRefine).
  m_refineWithPeaks->checked().connect( this, [this](){
    m_fitOrderCombo->enable();
    m_sharedPeaksValid = false;
    handleOptionsChange();
  } );
  m_refineWithPeaks->unChecked().connect( this, [this](){
    m_fitOrderCombo->disable();
    m_peaksTable->setHidden( true );
    handleOptionsChange();
  } );
  m_fitOrderCombo->activated().connect( this, [this]( int ){
    m_renderFlags |= UpdateRefine;
    scheduleRender();
  } );
  m_showOriginal->checked().connect( this, [this](){
    m_renderFlags |= UpdatePreview;
    scheduleRender();
  } );
  m_showOriginal->unChecked().connect( this, [this](){
    m_renderFlags |= UpdatePreview;
    scheduleRender();
  } );

  // Footer buttons
  WContainerWidget * const footer = parent->footer();
  AuxWindow::addHelpInFooter( footer, "energy-calibration" );
  m_cancel = footer->addNew<WPushButton>( WString::tr("Cancel") );
  m_cancel->clicked().connect( this, [this](){ handleFinish( Wt::DialogCode::Rejected ); } );
  m_use = footer->addNew<WPushButton>( WString::tr("Use") );
  m_use->clicked().connect( this, [this](){ handleFinish( Wt::DialogCode::Accepted ); } );
  m_use->disable();

  m_renderFlags |= UpdateRows;
  m_renderFlags |= UpdateCalc;
  scheduleRender();

  parent->resizeScaledWindow( 0.85, 0.85 );
  parent->centerWindow();
}//EnergyCalGainMatch constructor


EnergyCalGainMatch::~EnergyCalGainMatch()
{
}


void EnergyCalGainMatch::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool update_rows = m_renderFlags.test( UpdateRows );
  const bool update_calc = update_rows || m_renderFlags.test( UpdateCalc );
  const bool update_refine = m_renderFlags.test( UpdateRefine );
  const bool update_preview = update_calc || update_refine || m_renderFlags.test( UpdatePreview );

  m_renderFlags = Wt::WFlags<GainMatchRenderActions>();

  if( update_rows )
    rebuildRows();

  if( update_calc )
    doCalcUpdate();       //runs the match (+ shared-peak search & refinement when enabled)
  else if( update_refine )
    doRefineUpdate();     //re-applies the refinement using the cached shared peaks

  if( update_preview )
    updatePreview();

  WContainerWidget::render( flags );
}//render( flags )


EnergyCalGainMatch::MatchMode EnergyCalGainMatch::currentMode() const
{
  const bool between = m_modeGroup && (m_modeGroup->checkedButton() == m_betweenFilesBtn);
  return between ? MatchMode::BetweenFiles : MatchMode::WithinFile;
}


SpecUtils::SpectrumType EnergyCalGainMatch::selectedFileType() const
{
  const int index = m_fileSelect ? m_fileSelect->currentIndex() : -1;
  if( (index < 0) || (index >= static_cast<int>(m_fileSelectTypes.size())) )
    return SpectrumType::Foreground;
  return m_fileSelectTypes[index];
}


void EnergyCalGainMatch::handleModeChange()
{
  const bool within = (currentMode() == MatchMode::WithinFile);
  if( m_fileSelect )
    m_fileSelect->setEnabled( within && m_fileSelect->count() );

  // Multi-peak refinement is a within-file (per-detector) concept; hide it between files.
  if( m_peaksRow )
    m_peaksRow->setHidden( !within );
  if( !within && m_peaksTable )
    m_peaksTable->setHidden( true );

  m_sharedPeaksValid = false;
  m_renderFlags |= UpdateRows;
  m_renderFlags |= UpdateCalc;
  scheduleRender();
}


void EnergyCalGainMatch::handleOptionsChange()
{
  m_sharedPeaksValid = false;  //the stage-1/2 result may change, invalidating the shared peaks
  m_renderFlags |= UpdateCalc;
  scheduleRender();
}


bool EnergyCalGainMatch::refineWithPeaksEnabled() const
{
  return m_refineWithPeaks && m_refineWithPeaks->isChecked()
         && (currentMode() == MatchMode::WithinFile);
}


int EnergyCalGainMatch::selectedFitOrder() const
{
  const int index = m_fitOrderCombo ? m_fitOrderCombo->currentIndex() : 1;
  return std::max( 1, std::min(3, index + 1) );  //index 0 -> order 1 (offset+gain)
}


void EnergyCalGainMatch::rebuildRows()
{
  m_rowTable->clear();
  m_rows.clear();
  m_refGroup = std::make_shared<WButtonGroup>();
  m_withinFileMeas = nullptr;
  m_betweenFileMeas.clear();

  // Header row
  const char * const header_keys[7] = { "ecgm-col-use", "ecgm-col-ref", "", "ecgm-col-spectrum",
                                        "ecgm-col-gain", "ecgm-col-offset", "ecgm-col-status" };
  for( int col = 0; col < 7; ++col )
  {
    WText * const txt = m_rowTable->addNew<WText>( *header_keys[col]
                                                   ? WString::tr(header_keys[col]) : WString() );
    txt->addStyleClass( "GainMatchHeaderCell" );
    if( col == 5 )
      txt->addStyleClass( "GainMatchOffsetCell" );
  }

  vector<WString> row_names;

  if( currentMode() == MatchMode::WithinFile )
  {
    const SpectrumType type = selectedFileType();
    const shared_ptr<SpecMeas> meas = m_interspec->measurment( type );
    m_withinFileMeas = meas;

    if( meas )
    {
      const set<int> &samples = m_interspec->displayedSamples( type );

      for( const string &det : gamma_detectors( m_interspec, type ) )
      {
        shared_ptr<const EnergyCalibration> cal;
        for( const int sample : samples )
        {
          const shared_ptr<const Measurement> m = meas->measurement( sample, det );
          if( m && (m->num_gamma_channels() > 4)
              && m->energy_calibration() && m->energy_calibration()->valid() )
          {
            cal = m->energy_calibration();
            break;
          }
        }//for( find a calibration for this detector )

        if( !cal )
          continue;

        shared_ptr<const Measurement> spectrum;
        try
        {
          spectrum = meas->sum_measurements( samples, {det}, cal );
        }catch( std::exception & )
        {
        }

        if( !spectrum )
          continue;

        Row row;
        row.det_name = det;
        row.spec_type = type;
        row.spectrum = spectrum;
        row.orig_cal = cal;
        m_rows.push_back( row );
        row_names.emplace_back( WString::fromUTF8(det) );
      }//for( loop over gamma detectors )
    }//if( meas )
  }else
  {
    const SpectrumType types[3] = { SpectrumType::Foreground, SpectrumType::Background,
                                    SpectrumType::SecondForeground };
    const char * const type_keys[3] = { "Foreground", "Background", "Secondary" };

    set<shared_ptr<SpecMeas>> seen;
    for( int i = 0; i < 3; ++i )
    {
      const shared_ptr<SpecMeas> meas = m_interspec->measurment( types[i] );
      if( !meas || seen.count(meas) )
        continue;
      seen.insert( meas );

      const shared_ptr<const Measurement> spectrum = m_interspec->displayedHistogram( types[i] );
      if( !spectrum || !spectrum->energy_calibration() || !spectrum->energy_calibration()->valid() )
        continue;

      m_betweenFileMeas[types[i]] = meas;

      Row row;
      row.spec_type = types[i];
      row.spectrum = spectrum;
      row.orig_cal = spectrum->energy_calibration();
      m_rows.push_back( row );
      row_names.emplace_back( WString::tr(type_keys[i]) );
    }//for( loop over spectrum types )
  }//if( within-file mode ) / else

  // Widgets for each row
  for( size_t i = 0; i < m_rows.size(); ++i )
  {
    Row &row = m_rows[i];

    row.use = m_rowTable->addNew<WCheckBox>();
    row.use->setChecked( true );
    row.use->checked().connect( this, &EnergyCalGainMatch::handleOptionsChange );
    row.use->unChecked().connect( this, &EnergyCalGainMatch::handleOptionsChange );

    row.reference = m_rowTable->addNew<WRadioButton>();
    m_refGroup->addButton( row.reference, static_cast<int>(i) );

    WText * const swatch = m_rowTable->addNew<WText>();
    swatch->addStyleClass( "GainMatchSwatch" );
    swatch->setAttributeValue( "style",
                    "background-color: " + row_color(m_interspec, i).cssText() + ";" );

    row.name = m_rowTable->addNew<WText>( row_names[i] );

    row.gain = m_rowTable->addNew<WText>( "--" );
    row.gain->addStyleClass( "GainMatchNumCell" );

    row.offset = m_rowTable->addNew<WText>( "--" );
    row.offset->addStyleClass( "GainMatchNumCell GainMatchOffsetCell" );

    row.status = m_rowTable->addNew<WText>();
    row.status->addStyleClass( "GainMatchStatusCell" );
  }//for( loop over rows )

  m_refGroup->checkedChanged().connect( this, [this]( WRadioButton * ){ handleOptionsChange(); } );

  // Default reference: the spectrum with the most counts in the current energy range
  //  (the foreground in between-files mode, since rows are in fore/back/secondary order
  //  and the foreground will normally have the most counts - but max-counts either way)
  const float lower = m_lowerEnergy->value();
  const float upper = upperEnergyLimit();
  size_t best_row = 0;
  double best_counts = -1.0;
  for( size_t i = 0; i < m_rows.size(); ++i )
  {
    const double counts = GainMatchCalc::counts_in_range( m_rows[i].spectrum, lower, upper );
    if( counts > best_counts )
    {
      best_counts = counts;
      best_row = i;
    }
  }
  if( best_row < m_rows.size() )
    m_refGroup->setCheckedButton( m_rows[best_row].reference );
}//rebuildRows()


void EnergyCalGainMatch::rebuildPeaksTable()
{
  m_peaksTable->clear();
  m_peakUseCbs.clear();

  if( m_sharedPeaks.empty() )
  {
    m_peaksTable->setHidden( true );
    return;
  }

  m_peaksTable->setHidden( false );

  const size_t ndet = m_lastInputs.size();
  for( const GainMatchCalc::SharedPeak &peak : m_sharedPeaks )
  {
    WContainerWidget * const rowdiv = m_peaksTable->addNew<WContainerWidget>();
    rowdiv->addStyleClass( "GainMatchPeakRow" );

    WCheckBox * const cb = rowdiv->addNew<WCheckBox>( WString::fromUTF8(peak.label) );
    cb->setChecked( true );
    cb->changed().connect( this, [this](){
      m_renderFlags |= UpdateRefine;
      scheduleRender();
    } );
    m_peakUseCbs.push_back( cb );

    WText * const found = rowdiv->addNew<WText>(
              WString::tr("ecgm-peak-found-in").arg(peak.num_detectors).arg(static_cast<int>(ndet)) );
    found->addStyleClass( "GainMatchPeakFound" );
  }//for( loop over shared peaks )
}//rebuildPeaksTable()


void EnergyCalGainMatch::applyRefineToRows()
{
  for( Row &row : m_rows )
    row.refined = false;

  if( !refineWithPeaksEnabled() || m_sharedPeaks.empty() || m_lastInputs.empty() )
    return;

  vector<bool> use( m_sharedPeaks.size(), true );
  for( size_t j = 0; (j < m_peakUseCbs.size()) && (j < use.size()); ++j )
    use[j] = m_peakUseCbs[j] && m_peakUseCbs[j]->isChecked();

  const vector<shared_ptr<const SpecUtils::EnergyCalibration>> cals
      = GainMatchCalc::refineWithSharedPeaks( m_lastInputs, m_lastStage2, m_sharedPeaks,
                                              use, selectedFitOrder() );

  for( size_t k = 0; (k < cals.size()) && (k < m_lastInputRows.size()); ++k )
  {
    if( !cals[k] )
      continue;
    const size_t rowidx = m_lastInputRows[k];
    if( rowidx < m_rows.size() )
    {
      m_rows[rowidx].result.updated_cal = cals[k];
      m_rows[rowidx].refined = true;
    }
  }
}//applyRefineToRows()


void EnergyCalGainMatch::doRefineUpdate()
{
  applyRefineToRows();

  bool have_result = false;
  for( const Row &row : m_rows )
    have_result = (have_result || (row.result.used && row.result.updated_cal));
  m_use->setEnabled( have_result );

  updateResultColumns();
}//doRefineUpdate()


float EnergyCalGainMatch::upperEnergyLimit()
{
  if( !m_upperEnergy || m_upperEnergy->text().empty() )
    return 0.0f;
  return m_upperEnergy->value();
}


void EnergyCalGainMatch::doCalcUpdate()
{
  m_status->setText( WString() );

  const WRadioButton * const checked_ref = m_refGroup ? m_refGroup->checkedButton() : nullptr;

  vector<GainMatchCalc::SpectrumInput> inputs;
  vector<size_t> input_rows;
  int ref_index = -1;

  for( size_t i = 0; i < m_rows.size(); ++i )
  {
    Row &row = m_rows[i];
    row.result = GainMatchCalc::SpectrumResult();

    // The reference spectrum always participates; other rows only when their Use is checked
    const bool is_ref = row.reference && (checked_ref == row.reference);
    if( !row.spectrum || (!is_ref && (!row.use || !row.use->isChecked())) )
      continue;

    if( is_ref )
      ref_index = static_cast<int>( inputs.size() );

    shared_ptr<SpecMeas> meas;
    if( currentMode() == MatchMode::WithinFile )
    {
      meas = m_withinFileMeas;
    }else
    {
      const auto pos = m_betweenFileMeas.find( row.spec_type );
      meas = (pos == end(m_betweenFileMeas)) ? nullptr : pos->second;
    }

    GainMatchCalc::SpectrumInput input;
    input.name = row.name ? row.name->text().toUTF8() : row.det_name;
    input.spectrum = row.spectrum;
    input.fit_prefs = meas ? meas->peakFitDetPrefs() : nullptr;

    input_rows.push_back( i );
    inputs.push_back( input );
  }//for( loop over rows )

  m_lastInputs = inputs;
  m_lastInputRows = input_rows;
  m_lastStage2 = GainMatchCalc::MatchResults();

  bool have_result = false;

  try
  {
    GainMatchCalc::MatchOptions options;
    options.lower_energy = m_lowerEnergy->value();
    options.upper_energy = upperEnergyLimit();
    // When refining with peaks, keep the coarse stage as gain-only; the polynomial fit handles
    //  the offset (and higher orders).  The single-peak refinement is not used by the dialog.
    options.fit_offset = (m_fitOffset->isChecked() && !refineWithPeaksEnabled());
    options.ref_peak_energy = 0.0;

    const GainMatchCalc::MatchResults results = GainMatchCalc::match( inputs, ref_index, options );
    m_lastStage2 = results;

    assert( results.results.size() == input_rows.size() );
    for( size_t k = 0; k < results.results.size() && k < input_rows.size(); ++k )
    {
      m_rows[ input_rows[k] ].result = results.results[k];
      m_rows[ input_rows[k] ].refined = false;
      if( results.results[k].updated_cal )
        have_result = true;
    }

    // If the reference was auto-selected (no radio checked somehow), reflect it in the GUI
    if( (ref_index < 0) && (results.reference_index < input_rows.size()) )
    {
      const Row &ref_row = m_rows[ input_rows[results.reference_index] ];
      if( ref_row.reference )
        m_refGroup->setCheckedButton( ref_row.reference );
    }

    // Multi-peak refinement: find peaks common to the detectors (on a worker thread, cached
    //  until the match changes), then fit a per-detector polynomial from the selected ones.
    if( refineWithPeaksEnabled() )
    {
      if( !m_sharedPeaksValid )
        launchSharedPeakSearch();

      applyRefineToRows();  //uses whatever peaks we already have (none while a search is running)
      have_result = false;
      for( const Row &row : m_rows )
        have_result = (have_result || (row.result.used && row.result.updated_cal));

      // Don't let the user apply the coarse (linear) result while the refinement is still
      //  computing - they explicitly asked for the peak refinement.
      if( m_findingSharedPeaks )
        have_result = false;
    }//if( refineWithPeaksEnabled() )
  }catch( std::exception &e )
  {
    m_status->setText( WString::tr("ecgm-err-calc").arg( string(e.what()) ) );
  }

  m_use->setEnabled( have_result );
  updateResultColumns();
}//doCalcUpdate()


void EnergyCalGainMatch::launchSharedPeakSearch()
{
  // Peaks already fit (and possibly nuclide-identified) in the displayed spectrum, so a shared
  //  peak that matches one can be snapped to the true gamma energy.  Gathered here on the
  //  session thread; only immutable data crosses to the worker.
  const SpectrumType type = selectedFileType();
  const shared_ptr<SpecMeas> fmeas = m_interspec->measurment( type );
  deque<shared_ptr<const PeakDef>> id_peaks;
  if( fmeas )
  {
    const shared_ptr<deque<shared_ptr<const PeakDef>>> p
                              = fmeas->peaks( m_interspec->displayedSamples(type) );
    if( p )
      id_peaks = *p;
  }

  m_sharedPeaksValid = true;   //don't relaunch for this same match while the search is in flight
  m_findingSharedPeaks = true;
  m_sharedPeaks.clear();
  rebuildPeaksTable();
  m_peaksStatus->setText( WString::tr("ecgm-finding-peaks") );

  const int generation = ++m_refineGeneration;
  const string widgetId = id();
  const string sessionId = wApp ? wApp->sessionId() : string();
  const vector<GainMatchCalc::SpectrumInput> inputs = m_lastInputs;
  const GainMatchCalc::MatchResults stage2 = m_lastStage2;
  const float lower = m_lowerEnergy->value();
  const float upper = upperEnergyLimit();

  Wt::WServer * const server = Wt::WServer::instance();
  if( !server || sessionId.empty() )
  {
    // No server (e.g. unusual test context) - fall back to a synchronous search.
    onSharedPeaksFound( generation,
                        GainMatchCalc::findSharedPeaks( inputs, stage2, id_peaks, lower, upper ) );
    return;
  }

  server->ioService().boost::asio::io_service::post( [=](){
    // --- worker thread: the CPU-heavy peak search on the (immutable) summed spectra ---
    const vector<GainMatchCalc::SharedPeak> peaks
          = GainMatchCalc::findSharedPeaks( inputs, stage2, id_peaks, lower, upper );

    // --- back on the session thread: only the inert widget id crosses the boundary; the widget
    //     tree lookup happens here, under the applications UpdateLock (see CLAUDE.md) ---
    server->post( sessionId, [widgetId, generation, peaks](){
      WApplication * const app = WApplication::instance();
      if( !app || !app->domRoot() )
        return;
      WWidget * const w = app->domRoot()->findById( widgetId );
      EnergyCalGainMatch * const self = dynamic_cast<EnergyCalGainMatch *>( w );
      if( self )
      {
        self->onSharedPeaksFound( generation, peaks );
        app->triggerUpdate();
      }
    } );
  } );
}//launchSharedPeakSearch()


void EnergyCalGainMatch::onSharedPeaksFound( int generation,
                                             const std::vector<GainMatchCalc::SharedPeak> &peaks )
{
  if( generation != m_refineGeneration )
    return;  //a newer search (option change) superseded this result

  m_findingSharedPeaks = false;

  // The user may have turned refinement off (or switched modes) while the search ran.
  if( !refineWithPeaksEnabled() )
  {
    m_peaksTable->setHidden( true );
    m_peaksStatus->setText( WString() );
    return;
  }

  m_sharedPeaks = peaks;
  rebuildPeaksTable();
  m_peaksStatus->setText( m_sharedPeaks.empty() ? WString::tr("ecgm-no-shared-peaks")
                 : WString::tr("ecgm-num-shared-peaks").arg(static_cast<int>(m_sharedPeaks.size())) );

  applyRefineToRows();

  bool have_result = false;
  for( const Row &row : m_rows )
    have_result = (have_result || (row.result.used && row.result.updated_cal));
  m_use->setEnabled( have_result );

  updateResultColumns();
  updatePreview();
}//onSharedPeaksFound()


void EnergyCalGainMatch::updateResultColumns()
{
  const WRadioButton * const checked_ref = m_refGroup ? m_refGroup->checkedButton() : nullptr;

  for( Row &row : m_rows )
  {
    if( !row.gain || !row.offset || !row.status )
      continue;

    string warning_txt;
    for( const GainMatchCalc::WarningType warning : row.result.warnings )
    {
      const char * const key = warning_txt_key( warning );
      if( *key )
        warning_txt += (warning_txt.empty() ? "" : ", ") + WString::tr(key).toUTF8();
    }

    if( row.reference && (checked_ref == row.reference) )
    {
      row.gain->setText( "1" );
      row.offset->setText( "0" );
      row.status->setText( warning_txt.empty() ? WString::tr("ecgm-ref-label")
                                               : WString::fromUTF8(warning_txt) );
    }else if( row.result.updated_cal )
    {
      double gain = row.result.gain, offset = row.result.offset;

      // For a peak-refined (possibly higher-order) row, show the effective linear gain/offset
      //  of the new calibration relative to the original, evaluated across the data.
      if( row.refined && row.orig_cal && row.orig_cal->valid() )
      {
        try
        {
          const double nchan = static_cast<double>( row.orig_cal->num_channels() );
          const double clo = 0.05*nchan, chi = 0.95*nchan;
          const double eo_lo = row.orig_cal->energy_for_channel( clo );
          const double eo_hi = row.orig_cal->energy_for_channel( chi );
          const double en_lo = row.result.updated_cal->energy_for_channel( clo );
          const double en_hi = row.result.updated_cal->energy_for_channel( chi );
          if( std::fabs(eo_hi - eo_lo) > 1.0E-6 )
          {
            gain = (en_hi - en_lo) / (eo_hi - eo_lo);
            offset = en_lo - gain*eo_lo;
          }
        }catch( std::exception & )
        {
        }
      }//if( refined )

      char buffer[64];
      snprintf( buffer, sizeof(buffer), "%.4f", gain );
      row.gain->setText( buffer );
      snprintf( buffer, sizeof(buffer), "%.1f", offset );
      row.offset->setText( buffer );

      WString status = row.refined ? WString::tr("ecgm-status-refined") : WString::tr("ecgm-status-ok");
      if( !warning_txt.empty() )
        status = WString::fromUTF8( warning_txt );
      row.status->setText( status );
    }else
    {
      row.gain->setText( "--" );
      row.offset->setText( "--" );
      row.status->setText( WString::fromUTF8(warning_txt) );
    }
  }//for( loop over rows )
}//updateResultColumns()


void EnergyCalGainMatch::updatePreview()
{
  const bool show_original = m_showOriginal->isChecked();
  const float lower = m_lowerEnergy->value();
  const float upper = upperEnergyLimit();

  const WRadioButton * const checked_ref = m_refGroup ? m_refGroup->checkedButton() : nullptr;

  double ref_counts = 0.0;
  float max_energy = 0.0f;
  for( const Row &row : m_rows )
  {
    if( !row.spectrum )
      continue;
    if( row.reference && (checked_ref == row.reference) )
      ref_counts = GainMatchCalc::counts_in_range( row.spectrum, lower, upper );
    if( row.spectrum->energy_calibration() && row.spectrum->energy_calibration()->valid() )
      max_energy = std::max( max_energy, row.spectrum->energy_calibration()->upper_energy() );
  }

  vector<pair<shared_ptr<const Measurement>,D3SpectrumExport::D3SpectrumOptions>> display;

  for( size_t i = 0; i < m_rows.size(); ++i )
  {
    const Row &row = m_rows[i];
    if( !row.spectrum )
      continue;

    shared_ptr<const Measurement> spectrum = row.spectrum;
    if( !show_original && row.result.updated_cal )
    {
      try
      {
        auto adjusted = make_shared<Measurement>( *row.spectrum );
        adjusted->set_energy_calibration( row.result.updated_cal );
        spectrum = adjusted;
      }catch( std::exception & )
      {
      }
    }//if( show the matched calibration )

    D3SpectrumExport::D3SpectrumOptions options;
    options.line_color = row_color( m_interspec, i ).cssText();
    options.spectrum_type = SpectrumType::Foreground;

    string title = row.name ? row.name->text().toUTF8() : row.det_name;
    if( row.reference && (checked_ref == row.reference) )
      title += " (" + WString::tr("ecgm-ref-label").toUTF8() + ")";
    options.title = title;

    const double row_counts = GainMatchCalc::counts_in_range( row.spectrum, lower, upper );
    options.display_scale_factor = ((ref_counts > 0.0) && (row_counts > 0.0))
                                     ? (ref_counts / row_counts) : 1.0;

    display.emplace_back( spectrum, options );
  }//for( loop over rows )

  m_chart->setMultipleSpectra( std::move(display), !m_chartXRangeSet );
  m_chartXRangeSet = true;

  // Show the energy range being matched over
  m_chart->removeAllDecorativeHighlightRegions();
  m_rangeHighlightId = 0;
  const float region_upper = (upper > lower) ? upper : max_energy;
  if( region_upper > lower )
    m_rangeHighlightId = m_chart->addDecorativeHighlightRegion( lower, region_upper,
                                WColor(125, 175, 225, 65),
                                D3SpectrumDisplayDiv::HighlightRegionFill::BelowData,
                                WString::tr("ecgm-match-range") );
}//updatePreview()


void EnergyCalGainMatch::applyChanges()
{
  try
  {
    if( currentMode() == MatchMode::WithinFile )
    {
      const SpectrumType type = selectedFileType();
      const shared_ptr<SpecMeas> meas = m_interspec->measurment( type );
      if( !meas || (meas != m_withinFileMeas) )
        throw runtime_error( WString::tr("ecgm-file-changed").toUTF8() );

      const set<int> samples = m_allSamples->isChecked() ? meas->sample_numbers()
                                             : m_interspec->displayedSamples( type );

      // Peaks are intentionally NOT translated: the references energy scale is unchanged, and
      //  the (previously smeared) summed features consolidate onto it, so peaks stay at
      //  (approximately) the correct energies.  The stale automated-search peaks are
      //  invalidated and re-searched by the apply helper; any user-fit peaks are worth a
      //  re-fit, which we hint at below.
      size_t num_updated = 0;
      if( refineWithPeaksEnabled() )
      {
        // Multi-peak refinement produces a full per-detector calibration (higher order), so apply
        //  the calibration directly rather than a linear gain/offset transform.
        vector<pair<string,shared_ptr<const EnergyCalibration>>> per_detector;
        for( const Row &row : m_rows )
        {
          if( row.result.updated_cal )
            per_detector.emplace_back( row.det_name, row.result.updated_cal );
        }
        num_updated = GainMatchCalc::applyDetectorCals( m_interspec, type, meas, samples, per_detector );
      }else
      {
        vector<tuple<string,double,double>> per_detector;
        for( const Row &row : m_rows )
        {
          if( row.result.updated_cal )
            per_detector.emplace_back( row.det_name, row.result.gain, row.result.offset );
        }
        num_updated = GainMatchCalc::applyDetectorGains( m_interspec, type, meas, samples, per_detector );
      }//if( refineWithPeaksEnabled() ) / else

      if( num_updated && !meas->sampleNumsWithPeaks().empty() )
        m_interspec->logMessage( WString::tr("ecgm-note-refit-peaks"), 1 );
      if( num_updated )
        m_interspec->logMessage( WString::tr("ecgm-applied-toast")
                                  .arg( static_cast<int>(num_updated) ), 1 );
      return;
    }//if( within-file mode )

    // Between-files mode: each files features move with its recalibration, so its peaks are
    //  translated (not re-searched) - handled below.
    size_t num_updated = 0;
    auto changes = make_shared<vector<GainMatchFileChange>>();

    {
      for( const Row &row : m_rows )
      {
        if( !row.result.updated_cal )
          continue;

        const shared_ptr<SpecMeas> meas = m_interspec->measurment( row.spec_type );
        const auto pos = m_betweenFileMeas.find( row.spec_type );
        if( !meas || (pos == end(m_betweenFileMeas)) || (meas != pos->second) )
          throw runtime_error( WString::tr("ecgm-file-changed").toUTF8() );

        GainMatchFileChange change;
        change.file = meas;
        change.type = row.spec_type;

        const set<int> scope = m_allSamples->isChecked() ? meas->sample_numbers()
                                       : m_interspec->displayedSamples( row.spec_type );

        for( const shared_ptr<const Measurement> &m : meas->measurements() )
        {
          if( !m || !scope.count(m->sample_number()) || (m->num_gamma_channels() < 5) )
            continue;

          const shared_ptr<const EnergyCalibration> old_cal = m->energy_calibration();
          if( !old_cal || !old_cal->valid() )
            continue;

          const auto new_cal = GainMatchCalc::transform_calibration( old_cal, row.result.gain,
                                                                     row.result.offset );
          change.cals.emplace_back( m, old_cal, new_cal );
        }//for( loop over measurements )

        // This whole files features move with the recalibration, so its peaks must move too
        //  (same semantics as EnergyCalTool::applyCalChange)
        const vector<string> disp_dets = m_interspec->detectorsToDisplay( row.spec_type );
        for( const set<int> &peak_samples : meas->sampleNumsWithPeaks() )
        {
          bool contained = true;
          for( const int sample : peak_samples )
            contained = (contained && scope.count(sample));
          if( !contained )
            continue;

          const shared_ptr<deque<shared_ptr<const PeakDef>>> old_peaks = meas->peaks( peak_samples );
          const shared_ptr<const deque<shared_ptr<const PeakDef>>> old_hints
                                          = meas->automatedSearchPeaks( peak_samples );
          const bool have_peaks = (old_peaks && !old_peaks->empty());
          const bool have_hints = (old_hints && !old_hints->empty());
          if( !have_peaks && !have_hints )
            continue;

          try
          {
            const shared_ptr<const EnergyCalibration> old_sum
                       = meas->suggested_sum_energy_calibration( peak_samples, disp_dets );
            if( !old_sum || !old_sum->valid() )
              throw runtime_error( "no summed energy calibration" );

            const auto new_sum = GainMatchCalc::transform_calibration( old_sum, row.result.gain,
                                                                       row.result.offset );

            if( have_peaks )
              change.peaks.emplace_back( peak_samples, *old_peaks,
                  EnergyCal::translatePeaksForCalibrationChange(*old_peaks, old_sum, new_sum) );

            if( have_hints )
              change.hint_peaks.emplace_back( peak_samples, *old_hints,
                  EnergyCal::translatePeaksForCalibrationChange(*old_hints, old_sum, new_sum) );
          }catch( std::exception & )
          {
            m_interspec->logMessage( WString::tr("ecgm-note-refit-peaks"), 2 );
          }
        }//for( loop over sample sets with peaks )

        if( !change.cals.empty() )
        {
          num_updated += 1;
          changes->push_back( std::move(change) );
        }
      }//for( loop over rows )
    }//( between-files mode block )

    if( changes->empty() )
      return;

    apply_gain_match_direction( false, *changes );

    UndoRedoManager * const undoManager = m_interspec->undoRedoManager();
    if( undoManager )
    {
      const shared_ptr<const vector<GainMatchFileChange>> applied = changes;
      auto undo = [applied](){ apply_gain_match_direction( true, *applied ); };
      auto redo = [applied](){ apply_gain_match_direction( false, *applied ); };
      undoManager->addUndoRedoStep( undo, redo, "Gain-match energy calibration" );
    }//if( undoManager )

    m_interspec->logMessage( WString::tr("ecgm-applied-toast")
                              .arg( static_cast<int>(num_updated) ), 1 );
  }catch( std::exception &e )
  {
    m_interspec->logMessage( WString::tr("ecgm-err-apply").arg( string(e.what()) ), 2 );
  }
}//applyChanges()


void EnergyCalGainMatch::handleFinish( Wt::DialogCode result )
{
  switch( result )
  {
    case Wt::DialogCode::Rejected:
    {
      UndoRedoManager * const undoManager = m_interspec->undoRedoManager();
      if( m_parent && undoManager )
      {
        auto undo = [](){
          InterSpec * const viewer = InterSpec::instance();
          EnergyCalTool * const tool = viewer ? viewer->energyCalTool() : nullptr;
          if( tool )
            tool->moreActionBtnClicked( MoreActionsIndex::GainMatch );
        };

        auto redo = [](){
          InterSpec * const viewer = InterSpec::instance();
          EnergyCalTool * const tool = viewer ? viewer->energyCalTool() : nullptr;
          if( tool )
            tool->cancelMoreActionWindow();
        };

        undoManager->addUndoRedoStep( undo, redo, "Cancel energy cal gain match" );
      }//if( m_parent && undoManager )
      break;
    }//case Rejected

    case Wt::DialogCode::Accepted:
      applyChanges();
      break;
  }//switch( result )

  if( m_parent )
    m_parent->hide();
}//handleFinish(...)
