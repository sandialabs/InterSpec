#ifndef PeakFitLMObjective_imp_hpp
#define PeakFitLMObjective_imp_hpp
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

/** The linear-parameter (amplitude + continuum) solve, and the Poisson deviance residual, for
 PeakFitLM's Poisson-likelihood fits (see `SPARSE_DATA_LIKELIHOOD_USE_CERES` in src/PeakFitLM.cpp).
 Included only by src/PeakFitLM.cpp (and its unit test), after "InterSpec/PeakFit_imp.hpp".

 `PeakFit::fit_amp_and_offset_imp(...)` is the solve for the default modified-Neyman chi2, and is
 deliberately left byte-identical (it has many callers, some with other Jet sizes).  The IRLS solve
 here builds the same design matrix, but:
   - takes arbitrary per-channel variances, and
   - does not clip `data - fixed contributions` at zero, which would bias a Poisson fit.

 TODO: fold `build_roi_linear_basis(...)` back into `fit_amp_and_offset_imp(...)`, behind a
       bit-compare of the default fit on a test corpus.
 */

#include <cmath>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <type_traits>

#include "Eigen/Dense"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakDists_imp.hpp"


namespace PeakFitLMObjective
{

/** The ROI's linear model: `counts = A * [poly coefs..., peak amps...] + fixed_counts + fixed_step_counts`.

 Mirrors the design matrix of `PeakFit::fit_amp_and_offset_imp(...)`, but un-weighted.
 */
template<typename ScalarType>
struct RoiLinearBasis
{
  size_t nbin = 0;
  size_t num_poly = 0;
  size_t npeaks = 0;
  bool cdf_step = false;

  /** Counts per unit coefficient, `nbin x (num_poly + npeaks)`.  For the CDF step types each peak's
   column also carries its step term `(SUM_k s_k*E'^k) * CDFbar_i * dx`.
   */
  Eigen::MatrixX<ScalarType> A;

  /** Fixed-amplitude peaks' counts (Gaussian part); `nbin` entries. */
  std::vector<ScalarType> fixed_counts;

  /** For the CDF step types, the fixed-amplitude peaks' step contribution (part of the continuum);
   empty otherwise. */
  std::vector<ScalarType> fixed_step_counts;

  /** For the CDF step types, `[npeaks][nbin]` step part of each peak column; empty otherwise. */
  std::vector<std::vector<ScalarType>> peak_step_counts;
};//struct RoiLinearBasis


/** Builds the linear model of a ROI; see `PeakFit::fit_amp_and_offset_imp(...)` for the meaning of
 every argument.  `basis_data` is only used by the data-step continua.
 */
template<typename PeakType, typename ScalarType>
RoiLinearBasis<ScalarType> build_roi_linear_basis( const float *x,
                                                   const float *basis_data,
                                                   const size_t nbin,
                                                   const PeakContinuum::OffsetType cont_type,
                                                   const std::type_identity_t<ScalarType> * const step_coeffs,
                                                   const ScalarType ref_energy,
                                                   const std::vector<ScalarType> &means,
                                                   const std::vector<ScalarType> &sigmas,
                                                   const std::vector<PeakType> &fixedAmpPeaks,
                                                   const PeakDef::SkewType skew_type,
                                                   const ScalarType *skew_parameters )
{
  using namespace std;

  RoiLinearBasis<ScalarType> basis;

  const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_linear_fit_pars( cont_type ) );
  const bool step_continuum = PeakContinuum::is_step_continuum( cont_type );
  const bool cdf_step = PeakContinuum::is_peak_cdf_step_continuum( cont_type );
  const size_t num_step_terms = PeakContinuum::num_cdf_step_pars( cont_type );

  if( sigmas.size() != means.size() )
    throw runtime_error( "build_roi_linear_basis: invalid input" );

  if( !skew_parameters && (skew_type != PeakDef::SkewType::NoSkew) )
    throw std::logic_error( "build_roi_linear_basis: skew pars not provided" );

  const bool have_step_coeffs = (num_step_terms > 0) && step_coeffs;
  const size_t npeaks = sigmas.size();
  const Eigen::Index num_poly_terms = static_cast<Eigen::Index>( num_polynomial_terms );

  basis.nbin = nbin;
  basis.num_poly = static_cast<size_t>( num_polynomial_terms );
  basis.npeaks = npeaks;
  basis.cdf_step = cdf_step;
  basis.A.resize( static_cast<Eigen::Index>(nbin), num_poly_terms + static_cast<Eigen::Index>(npeaks) );

  // The data-step basis, exactly as in `fit_amp_and_offset_imp(...)` and `PeakContinuum::offset_integral(...)`.
  ScalarType roi_data_sum( 0.0 ), step_cumulative_data( 0.0 );
#if( PEAK_CONTINUUM_DATA_STEP_SUBTRACT )
  double min_data_val = static_cast<double>( basis_data[0] );
  for( size_t row = 0; row < nbin; ++row )
  {
    roi_data_sum += (std::max)( static_cast<double>( basis_data[row] ), 0.0 );
    min_data_val = (std::min)( min_data_val, static_cast<double>( basis_data[row] ) );
  }
  if( step_continuum && !cdf_step )
  {
    min_data_val = (std::max)( min_data_val, 0.0 );
    roi_data_sum -= ScalarType( min_data_val * static_cast<double>( nbin ) );
  }else
  {
    min_data_val = 0.0;
  }
#else
  const double min_data_val = 0.0;
  for( size_t row = 0; row < nbin; ++row )
    roi_data_sum += (std::max)( static_cast<double>(basis_data[row]), 0.0 );
#endif

  const auto step_coeff_for_channel = [&]( const size_t channel ) -> ScalarType {
    if( !have_step_coeffs )
      return ScalarType( 0.0 );

    const ScalarType center_rel = ScalarType(0.5) * ((ScalarType(x[channel]) - ref_energy)
                                                      + (ScalarType(x[channel+1]) - ref_energy));
    ScalarType step_coeff = step_coeffs[0];
    ScalarType energy_pow = center_rel;
    for( size_t k = 1; k < num_step_terms; ++k )
    {
      step_coeff += step_coeffs[k] * energy_pow;
      energy_pow *= center_rel;
    }

    return step_coeff;
  };//step_coeff_for_channel

  const size_t nfixedpeak = fixedAmpPeaks.size();
  basis.fixed_counts.assign( nbin, ScalarType(0.0) );
  for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
    fixedAmpPeaks[peak_index].gauss_integral( x, basis.fixed_counts.data(), nbin );

  if( nfixedpeak && cdf_step )
  {
    basis.fixed_step_counts.assign( nbin, ScalarType(0.0) );
    vector<ScalarType> pdf_per_channel( nbin ), cdf_per_channel( nbin );

    for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
    {
      const ScalarType &amp_val = fixedAmpPeaks[peak_index].amplitude();
      const PeakDef::SkewType skew = fixedAmpPeaks[peak_index].skewType();

      const ScalarType *fixed_skew_pars = nullptr;
      if( skew != PeakDef::SkewType::NoSkew )
      {
        if constexpr ( std::is_same_v<PeakType, PeakDef> )
          fixed_skew_pars = fixedAmpPeaks[peak_index].coefficients() + PeakDef::CoefficientType::SkewPar0;
        else
          fixed_skew_pars = fixedAmpPeaks[peak_index].skew_parameters();
      }

      std::fill( begin(pdf_per_channel), end(pdf_per_channel), ScalarType(0.0) );
      PeakDists::photopeak_function_integral( fixedAmpPeaks[peak_index].mean(),
                                              fixedAmpPeaks[peak_index].sigma(), ScalarType(1.0),
                                              skew, fixed_skew_pars, nbin, x, pdf_per_channel.data() );
      PeakDists::unit_pdf_to_cdf( pdf_per_channel.data(), cdf_per_channel.data(), nbin );

      for( size_t channel = 0; channel < nbin; ++channel )
      {
        const ScalarType dx = ScalarType( x[channel+1] - x[channel] );
        basis.fixed_step_counts[channel] += step_coeff_for_channel(channel)
                                            * amp_val * cdf_per_channel[channel] * dx;
      }
    }//for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
  }//if( nfixedpeak && cdf_step )

  for( size_t row = 0; row < nbin; ++row )
  {
    const ScalarType dataval = ScalarType( static_cast<double>(basis_data[row]) );

    const double x0 = x[row];
    const double x1 = x[row+1];

    const ScalarType x0_rel = x0 - ref_energy;
    const ScalarType x1_rel = x1 - ref_energy;

    if( step_continuum && !cdf_step )
      step_cumulative_data += (dataval - ScalarType( min_data_val ));

    for( Eigen::Index col = 0; col < num_poly_terms; ++col )
    {
      const ScalarType exp = ScalarType(col + 1.0);

      if( !cdf_step && step_continuum
         && ((num_polynomial_terms == 2) || (num_polynomial_terms == 3))
         && (col == (num_polynomial_terms - 1)) )
      {
        const ScalarType frac_data = (roi_data_sum > 0.0)
            ? (step_cumulative_data - ScalarType( 0.5*(basis_data[row] - min_data_val) )) / roi_data_sum
            : ScalarType( 0.5 );
        basis.A(row,col) = ScalarType( frac_data * (x1 - x0) );
      }else if( !cdf_step && step_continuum && (num_polynomial_terms == 4) )
      {
        const ScalarType frac_data = (roi_data_sum > 0.0)
            ? (step_cumulative_data - ScalarType( 0.5*(basis_data[row] - min_data_val) )) / roi_data_sum
            : ScalarType( 0.5 );

        ScalarType contrib( 0.0 );
        switch( col )
        {
          case 0: contrib = (1.0 - frac_data) * (x1_rel - x0_rel);                     break;
          case 1: contrib = 0.5 * (1.0 - frac_data) * (x1_rel*x1_rel - x0_rel*x0_rel); break;
          case 2: contrib = frac_data * (x1_rel - x0_rel);                             break;
          case 3: contrib = 0.5 * frac_data * (x1_rel*x1_rel - x0_rel*x0_rel);         break;
          default: assert( 0 ); break;
        }//switch( col )

        basis.A(row,col) = contrib;
      }else
      {
        basis.A(row,col) = (1.0/exp) * (pow(x1_rel,exp) - pow(x0_rel,exp));
      }
    }//for( Eigen::Index col = 0; col < num_poly_terms; ++col )
  }//for( size_t row = 0; row < nbin; ++row )

  vector<ScalarType> peak_counts( nbin ), cdf_at_centers( cdf_step ? nbin : size_t(0) );
  if( cdf_step )
    basis.peak_step_counts.assign( npeaks, vector<ScalarType>(nbin,ScalarType(0.0)) );

  for( size_t i = 0; i < npeaks; ++i )
  {
    std::fill( begin(peak_counts), end(peak_counts), ScalarType(0.0) );
    PeakDists::photopeak_function_integral( means[i], sigmas[i], ScalarType(1.0),
                                            skew_type, skew_parameters, nbin, x, peak_counts.data() );

    if( cdf_step )
    {
      PeakDists::unit_pdf_to_cdf( peak_counts.data(), cdf_at_centers.data(), nbin );

      for( size_t channel = 0; channel < nbin; ++channel )
      {
        const ScalarType dx = ScalarType( x[channel+1] - x[channel] );
        basis.peak_step_counts[i][channel] = step_coeff_for_channel(channel) * cdf_at_centers[channel] * dx;
        peak_counts[channel] += basis.peak_step_counts[i][channel];
      }
    }//if( cdf_step )

    for( size_t channel = 0; channel < nbin; ++channel )
      basis.A(static_cast<Eigen::Index>(channel), num_poly_terms + static_cast<Eigen::Index>(i)) = peak_counts[channel];
  }//for( size_t i = 0; i < npeaks; ++i )

  return basis;
}//build_roi_linear_basis(...)


/** Adds the model counts for the given coefficients to `peak_counts`; the continuum (polynomial plus
 any CDF step) is clamped at zero per channel, exactly as `fit_amp_and_offset_imp(...)` and
 `PeakContinuum::offset_integral(...)` do.  `extra_fixed` (may be null) is added as is.
 */
template<typename ScalarType>
void add_model_counts( const RoiLinearBasis<ScalarType> &basis,
                       const Eigen::VectorX<ScalarType> &coeffs,
                       const double * const extra_fixed,
                       ScalarType * const peak_counts )
{
  const Eigen::Index num_poly = static_cast<Eigen::Index>( basis.num_poly );

  for( size_t bin = 0; bin < basis.nbin; ++bin )
  {
    const Eigen::Index row = static_cast<Eigen::Index>( bin );

    ScalarType continuum_bin( 0.0 );
    for( Eigen::Index col = 0; col < num_poly; ++col )
      continuum_bin += coeffs(col) * basis.A(row,col);

    if( basis.cdf_step )
    {
      for( size_t i = 0; i < basis.npeaks; ++i )
        continuum_bin += coeffs(num_poly + static_cast<Eigen::Index>(i)) * basis.peak_step_counts[i][bin];
      if( !basis.fixed_step_counts.empty() )
        continuum_bin += basis.fixed_step_counts[bin];
    }

    if( continuum_bin < 0.0 )
      continuum_bin = ScalarType(0.0);

    ScalarType y_pred = continuum_bin;
    for( size_t i = 0; i < basis.npeaks; ++i )
    {
      const Eigen::Index col = num_poly + static_cast<Eigen::Index>(i);
      y_pred += coeffs(col) * (basis.cdf_step ? (basis.A(row,col) - basis.peak_step_counts[i][bin])
                                              : basis.A(row,col));
    }

    y_pred += basis.fixed_counts[bin];
    if( extra_fixed )
      y_pred += ScalarType( extra_fixed[bin] );

    peak_counts[bin] += y_pred;
  }//for( size_t bin = 0; bin < basis.nbin; ++bin )
}//add_model_counts(...)


/** Weighted linear least squares of `data ~ A*coeffs + fixed + extra_fixed`, with weights
 `1/variances`, no clipping of the data; coefficient uncertainties from `(A^T W A)^-1`.

 `data` (which also defines any data-step continuum), `variances` and (optional) `extra_fixed` have
 `nbin` entries; all variances must be positive.  The model counts are added to `peak_counts` (see
 `add_model_counts(...)`).  Arguments otherwise as `PeakFit::fit_amp_and_offset_imp(...)`.
 Throws on failure.
 */
template<typename PeakType, typename ScalarType>
void fit_amp_and_offset_weighted( const float *x,
                                  const float *data,
                                  const double *variances,
                                  const double *extra_fixed,
                                  const size_t nbin,
                                  const PeakContinuum::OffsetType cont_type,
                                  const std::type_identity_t<ScalarType> * const step_coeffs,
                                  const ScalarType ref_energy,
                                  const std::vector<ScalarType> &means,
                                  const std::vector<ScalarType> &sigmas,
                                  const std::vector<PeakType> &fixedAmpPeaks,
                                  const PeakDef::SkewType skew_type,
                                  const ScalarType *skew_parameters,
                                  std::vector<ScalarType> &amplitudes,
                                  std::vector<ScalarType> &continuum_coeffs,
                                  std::vector<ScalarType> &amplitudes_uncerts,
                                  std::vector<ScalarType> &continuum_coeffs_uncerts,
                                  ScalarType * const peak_counts )
{
  const RoiLinearBasis<ScalarType> basis
         = build_roi_linear_basis<PeakType,ScalarType>( x, data, nbin, cont_type, step_coeffs, ref_energy,
                                                        means, sigmas, fixedAmpPeaks, skew_type, skew_parameters );

  const Eigen::Index nrow = static_cast<Eigen::Index>( nbin );
  const Eigen::Index ncol = basis.A.cols();
  Eigen::MatrixX<ScalarType> Aw( nrow, ncol );
  Eigen::VectorX<ScalarType> yw( nrow );

  for( Eigen::Index row = 0; row < nrow; ++row )
  {
    const size_t bin = static_cast<size_t>( row );
    assert( variances[bin] > 0.0 );
    const double inv_uncert = 1.0 / std::sqrt( variances[bin] );

    ScalarType y = ScalarType( static_cast<double>( data[bin] ) ) - basis.fixed_counts[bin];
    if( !basis.fixed_step_counts.empty() )
      y -= basis.fixed_step_counts[bin];
    if( extra_fixed )
      y -= ScalarType( extra_fixed[bin] );

    yw(row) = y * inv_uncert;
    for( Eigen::Index col = 0; col < ncol; ++col )
      Aw(row,col) = basis.A(row,col) * inv_uncert;
  }//for( rows )

#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<ScalarType>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd( Aw );
#else
  const Eigen::BDCSVD<Eigen::MatrixX<ScalarType>> svd( Aw, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  const Eigen::VectorX<ScalarType> coeffs = svd.solve( yw );
  check_jet_array_for_NaN( coeffs.data(), coeffs.size() );

  const Eigen::MatrixX<ScalarType> covariance = (Aw.transpose() * Aw).inverse();

  const size_t num_poly = basis.num_poly;
  continuum_coeffs.resize( num_poly );
  continuum_coeffs_uncerts.resize( num_poly );
  for( size_t i = 0; i < num_poly; ++i )
  {
    const Eigen::Index k = static_cast<Eigen::Index>( i );
    continuum_coeffs[i] = coeffs(k);
    continuum_coeffs_uncerts[i] = sqrt( covariance(k,k) );
  }

  amplitudes.resize( basis.npeaks );
  amplitudes_uncerts.resize( basis.npeaks );
  for( size_t i = 0; i < basis.npeaks; ++i )
  {
    const Eigen::Index k = static_cast<Eigen::Index>( num_poly + i );
    amplitudes[i] = coeffs(k);
    amplitudes_uncerts[i] = sqrt( covariance(k,k) );
  }

  if( peak_counts )
    add_model_counts( basis, coeffs, extra_fixed, peak_counts );
}//fit_amp_and_offset_weighted(...)


/** Poisson deviance of one channel; `expected` must be positive. */
inline double poisson_deviance_term( const double observed, const double expected )
{
  if( observed <= 0.0 )
    return 2.0*expected;
  return 2.0*( expected - observed + observed*std::log(observed/expected) );
}


/** Signed square root of the Poisson deviance of one channel: `r*r` is the deviance, and `r` has
 the sign of `observed - expected`, so a least-squares minimizer of `r` minimizes the deviance.

 Written as `r = (n - m)/sqrt(n) * psi(m/n)`, `psi(t) = sqrt(2*(t - 1 - ln t))/|t - 1|`, using the
 series of `psi` near `t == 1`, so the derivative stays finite where the data and model agree (a
 plain `sqrt(deviance)` has an infinite derivative there).  `psi == 1` would be the Neyman residual.
 Expected counts below `floor` (> 0) are replaced by the smooth, positive continuation
 `floor^2/(2*floor - m)` (equal value and slope at `floor`, and never underflowing to zero), so the
 residual stays finite, and C1, for a model that touches or crosses zero.
 */
template<typename T>
T poisson_deviance_residual( const double observed, const T &expected, const double floor )
{
  using std::log;
  using std::sqrt;

  const T m = (expected >= floor) ? expected : T( (floor*floor) / (2.0*floor - expected) );

  if( observed <= 0.0 )
    return -sqrt( 2.0*m );

  const double n = observed;
  const T u = m/n - 1.0;
  if( (u < 1.0E-3) && (u > -1.0E-3) )
  {
    const T psi = 1.0 + u*(-1.0/3.0 + u*(7.0/36.0 - u*0.13518518518518519));
    return -sqrt( n ) * u * psi;
  }

  const T d = 2.0*n*( u - log( m/n ) );   // = deviance; positive away from u == 0
  return (u < 0.0) ? sqrt( d ) : -sqrt( d );
}//poisson_deviance_residual(...)


}//namespace PeakFitLMObjective

#endif //PeakFitLMObjective_imp_hpp
