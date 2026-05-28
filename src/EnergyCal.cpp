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
#include <deque>
#include <limits>
#include <vector>
#include <numeric>
#include <cassert>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "Eigen/Dense"

#include "ceres/ceres.h"

#include "InterSpec/PeakDef.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/EnergyCal.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/RelActCalcAuto_EnergyCal_imp.hpp"

using namespace std;


namespace
{
double poly_coef_fcn( size_t order, double channel, size_t nchannel )
{
  return pow( channel, static_cast<double>(order) );
}

double frf_coef_fcn( size_t order, double channel, size_t nchannel )
{
  const double x = channel / nchannel;
  if( order == 4 )
    return 1.0 / (1.0 + 60.0*x);
  return pow( x, static_cast<double>(order) );
}


double fit_energy_cal_imp( const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                                      const std::vector<bool> &fitfor,
                                      const size_t nchannels,
                                      const std::vector<std::pair<float,float>> &dev_pairs,
                                      std::vector<float> &coefs,
                                      std::vector<float> &coefs_uncert,
                                      std::function<double(size_t,double,size_t)> coeffcn )
{
  assert( coeffcn );
  
  const size_t npeaks = peakinfos.size();
  const size_t nparsfit = static_cast<size_t>( std::count(begin(fitfor),end(fitfor),true) );
  
  if( npeaks < 1 )
    throw runtime_error( "Must have at least one peak" );
  
  if( nparsfit < 1 )
    throw runtime_error( "Must fit for at least one coefficient" );
  
  if( nparsfit > npeaks )
    throw runtime_error( "Must have at least as many peaks as coefficients fitting for" );
  
  if( (nparsfit != fitfor.size()) && (coefs.size() != fitfor.size()) )
    throw runtime_error( "You must supply input coefficient when any of the coefficients are fixed" );
  
  //Energy = P0 + P1*x + P2*x^2 + P3*x^3, where x is bin number
  //  However, some of the coeffeicents may not be being fit for.
  vector<float> mean_bin( npeaks ), true_energies( npeaks ), energy_uncerts( npeaks );
  for( size_t i = 0; i < npeaks; ++i )
  {
    mean_bin[i] = peakinfos[i].peakMeanBinNumber;
    true_energies[i] = peakinfos[i].photopeakEnergy;
    energy_uncerts[i] = true_energies[i] * peakinfos[i].peakMeanUncert / std::max(peakinfos[i].peakMean,1.0);
  }
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation using Eigen SVD for numerical stability
  //Energy_i = P0 + P1*pow(i,1) + P2*pow(i,2) + P3*pow(i,3)  (for polynomial)
  
  Eigen::MatrixX<double> A( npeaks, nparsfit );
  Eigen::VectorX<double> b( npeaks );
  
  for( size_t row = 0; row < npeaks; ++row )
  {
    double data_y = true_energies[row];
    const double data_y_uncert = fabs( energy_uncerts[row] );
    
    data_y -= SpecUtils::correction_due_to_dev_pairs( true_energies[row], dev_pairs );
    
    for( size_t col = 0, coef_index = 0; coef_index < fitfor.size(); ++coef_index )
    {
      if( fitfor[coef_index] )
      {
        assert( col < nparsfit );
        A(row,col) = coeffcn( coef_index, mean_bin[row], nchannels ) / data_y_uncert; //std::pow( mean_bin[row], double(coef_index)) / data_y_uncert;
        ++col;
      }else
      {
        data_y -= coefs[coef_index] * coeffcn( coef_index, mean_bin[row], nchannels);
      }
    }//
    
    b(row) = data_y / data_y_uncert;
  }//for( int col = 0; col < order; ++col )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  
  // Compute uncertainties using normal equations for compatibility
  // Note: For maximum numerical stability, we could extract uncertainties from the SVD,
  // but for consistency with existing code behavior, we use the normal equations approach
  const Eigen::MatrixX<double> A_transpose = A.transpose();
  const Eigen::MatrixX<double> alpha = A_transpose * A;
  const Eigen::MatrixX<double> C = alpha.inverse();
  
  coefs.resize( fitfor.size(), 0.0 );
  coefs_uncert.resize( fitfor.size(), 0.0 );
  
  for( size_t col = 0, coef_index = 0; coef_index < fitfor.size(); ++coef_index )
  {
    if( fitfor[coef_index] )
    {
      assert( col < nparsfit );
      coefs[coef_index] = static_cast<float>( a(col) );
      coefs_uncert[coef_index] = static_cast<float>( std::sqrt( C(col,col) ) );
      ++col;
    }else
    {
      coefs_uncert[coef_index] = 0.0;
    }
  }//for( int coef = 0; coef < order; ++coef )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < npeaks; ++bin )
  {
    double y_pred = 0.0;
    for( size_t i = 0; i < fitfor.size(); ++i )
      y_pred += coefs[i] * coeffcn( i, mean_bin[bin], nchannels );
    y_pred += SpecUtils::deviation_pair_correction( y_pred, dev_pairs );
    chi2 += std::pow( (y_pred - true_energies[bin]) / energy_uncerts[bin], 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_energy_cal_imp


double fit_from_channel_energies_imp( const size_t ncoeffs, const vector<float> &channel_energies,
                                      std::function<double(size_t,double,size_t)> coeffcn,
                                      vector<float> &coefs )
{
  if( ncoeffs < 2 )
    throw runtime_error( "fit_from_channel_energies_imp: You must request at least two coefficients" );
    
  //For polynomial this isnt a programming or math limitation, just a sanity limitation; FRF should
  //  have max of 5 coefficients
  if( ncoeffs >= 6 )
    throw runtime_error( "fit_from_channel_energies_imp: You must request less than 6 coefficients" );
  
  const size_t nenergies = channel_energies.size();
  if( nenergies <= 6 )
    throw runtime_error( "fit_from_channel_energies_imp: Input energies must have at least 6 entries" );
  
  const size_t nchannel = nenergies - 1;
  
  for( size_t i = 1; i < nenergies; ++i )
    if( channel_energies[i-1] >= channel_energies[i] )
      throw runtime_error( "fit_from_channel_energies_imp: Input energies must be stricktly increasing" );
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation using Eigen SVD for numerical stability
  
  Eigen::MatrixX<double> A( nenergies, ncoeffs );
  Eigen::VectorX<double> b( nenergies );
  
  for( size_t row = 0; row < nenergies; ++row )
  {
    //Energy_i = P0 + P1*pow(i,1) + P2*pow(i,2) + P3*pow(i,3)  //for polynomial
    const double uncert = 1.0;
    for( size_t col = 0; col < ncoeffs; ++col )
      A(row,col) = coeffcn(col, row, nchannel) / uncert;
    b(row) = channel_energies[row] / uncert;
  }//for( int col = 0; col < order; ++col )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  
  coefs.resize( ncoeffs, 0.0 );
  for( size_t col = 0; col < ncoeffs; ++col )
    coefs[col] = static_cast<float>( a(col) );
  
  double avrgdiffs = 0.0;
  for( size_t i = 0; i < nenergies; ++i )
  {
    double energy = 0.0;
    for( size_t col = 0; col < ncoeffs; ++col )
      energy += coefs[col] * coeffcn(col, i, nchannel);
    avrgdiffs += fabs( energy - channel_energies[i] );
  }
  
  return avrgdiffs / nenergies;
}//fit_from_channel_energies_imp



/** Stride for Ceres autodiff (max expected number of fit parameters).
 If exceeded, Ceres handles it by chunking - just slightly less efficient.
 Energy cal fits are small (poly/FRF order ≤ 5, LowerChannelEdge fits 2 params).
 */
constexpr int sm_energy_cal_ceres_stride = 8;


/** Linear interpolation of the baseline lower-channel-edge vector at fractional channel `channel`,
 transformed by the relative affine model
   new[i] = baseline[i] + offset + (gain - 1) * (baseline[i] - baseline[0])
          = baseline[0] + offset + gain * (baseline[i] - baseline[0])
 where `offset` is the *delta* applied at channel 0 (so 0 means no shift) and `gain` is a
 multiplier on the energy span (so 1 means identity). This matches the values shown to the
 user in the EnergyCalTool UI.

 `baseline` has `nchannel + 1` entries (last entry is upper edge of last channel).
 Extrapolation outside `[0, nchannel]` uses the first/last channel width.
 Auto-diff gradient flows through `offset` and `gain`; the baseline values are constants.
 */
template<typename T>
T predicted_lce_energy( const T &offset,
                        const T &gain,
                        const double channel,
                        const std::vector<float> &baseline )
{
  assert( baseline.size() >= 2 );
  const size_t nchannel = baseline.size() - 1;
  const double base0 = static_cast<double>( baseline[0] );

  double base_at_channel;
  if( channel < 0.0 )
  {
    const double e0 = static_cast<double>( baseline[0] );
    const double width = static_cast<double>( baseline[1] ) - e0;
    base_at_channel = e0 + width * channel;
  }
  else if( channel >= static_cast<double>(nchannel) )
  {
    const double e_last = static_cast<double>( baseline[nchannel] );
    const double width  = e_last - static_cast<double>( baseline[nchannel - 1] );
    base_at_channel = e_last + width * (channel - static_cast<double>(nchannel));
  }
  else
  {
    const size_t low = static_cast<size_t>( std::floor(channel) );
    assert( (low + 1) < baseline.size() );
    const double e_lo = static_cast<double>( baseline[low] );
    const double e_hi = static_cast<double>( baseline[low + 1] );
    const double frac = channel - static_cast<double>(low);
    base_at_channel = e_lo + (e_hi - e_lo) * frac;
  }

  return T(base0) + offset + gain * T( base_at_channel - base0 );
}//predicted_lce_energy(...)


/** Ceres autodiff cost functor for fitting energy calibration coefficients.

 Handles Polynomial, FullRangeFraction (with deviation pairs via the templated cubic spline in
 RelActCalcAuto_EnergyCal_imp.hpp), and LowerChannelEdge (with offset/gain affine transform of a
 baseline lower-channel-edge vector — `params[0] = offset`, `params[1] = gain`).

 Each residual i is `(predicted_energy_i - photopeak_energy_i) / uncert_i`, so Ceres' final cost
 `0.5 * sum(r_i^2)` equals half the Minuit-style chi-squared. The caller returns `2*final_cost`
 to preserve the chi2 return contract.
 */
struct EnergyCalFitCost
{
  SpecUtils::EnergyCalType m_eqnType;
  size_t m_nchannel;
  size_t m_npars;
  std::vector<EnergyCal::RecalPeakInfo> m_peakInfo;
  std::vector<std::pair<float,float>> m_dev_pairs;
  std::vector<float> m_baseline_lce_edges; // only populated when m_eqnType == LowerChannelEdge

  template<typename T>
  bool operator()( T const * const * parameters, T *residuals ) const
  {
    const T * const pars = parameters[0];

    if( m_eqnType == SpecUtils::EnergyCalType::LowerChannelEdge )
    {
      for( size_t i = 0; i < m_peakInfo.size(); ++i )
      {
        const T predE = predicted_lce_energy<T>( pars[0], pars[1],
                                                 m_peakInfo[i].peakMeanBinNumber,
                                                 m_baseline_lce_edges );
        const double w = (m_peakInfo[i].peakMeanUncert <= 0.0) ? 1.0 : m_peakInfo[i].peakMeanUncert;
        residuals[i] = (predE - T(m_peakInfo[i].photopeakEnergy)) / T(w);
      }
      return true;
    }

    // Build templated dev-pair spline. Anchors are constant doubles, so the offsets
    // are wrapped as T(0)-derivative jets — gradient flows only through `pars`.
    std::vector<std::pair<double,T>> dev_pairs_T;
    dev_pairs_T.reserve( m_dev_pairs.size() );
    for( const auto &dp : m_dev_pairs )
      dev_pairs_T.emplace_back( static_cast<double>(dp.first), T(static_cast<double>(dp.second)) );
    const std::vector<RelActCalcAutoImp::CubicSplineNodeT<T>> spline
      = RelActCalcAutoImp::create_cubic_spline( dev_pairs_T );

    std::vector<T> coefs( pars, pars + m_npars );

    if( m_eqnType == SpecUtils::EnergyCalType::FullRangeFraction )
    {
      for( size_t i = 0; i < m_peakInfo.size(); ++i )
      {
        const T predE = RelActCalcAutoImp::fullrangefraction_energy<T>(
          T(m_peakInfo[i].peakMeanBinNumber), coefs, m_nchannel, spline );
        const double w = (m_peakInfo[i].peakMeanUncert <= 0.0) ? 1.0 : m_peakInfo[i].peakMeanUncert;
        residuals[i] = (predE - T(m_peakInfo[i].photopeakEnergy)) / T(w);
      }
    }
    else
    {
      // Polynomial or UnspecifiedUsingDefaultPolynomial
      for( size_t i = 0; i < m_peakInfo.size(); ++i )
      {
        const T predE = RelActCalcAutoImp::polynomial_energy<T>(
          T(m_peakInfo[i].peakMeanBinNumber), coefs, spline );
        const double w = (m_peakInfo[i].peakMeanUncert <= 0.0) ? 1.0 : m_peakInfo[i].peakMeanUncert;
        residuals[i] = (predE - T(m_peakInfo[i].photopeakEnergy)) / T(w);
      }
    }

    return true;
  }//operator()
};//struct EnergyCalFitCost

}//namespace





double EnergyCal::fit_energy_cal_frf( const std::vector<EnergyCal::RecalPeakInfo> &peaks,
                                      const std::vector<bool> &fitfor,
                                      const size_t nchannels,
                                      const std::vector<std::pair<float,float>> &dev_pairs,
                                      std::vector<float> &coefs,
                                      std::vector<float> &uncert )
{
  return fit_energy_cal_imp( peaks, fitfor, nchannels, dev_pairs, coefs, uncert, &frf_coef_fcn );
}


double EnergyCal::fit_energy_cal_poly( const std::vector<EnergyCal::RecalPeakInfo> &peaks,
                            const vector<bool> &fitfor,
                            const size_t nchannels,
                            const std::vector<std::pair<float,float>> &dev_pairs,
                            vector<float> &coefs,
                            vector<float> &uncert )
{
  return fit_energy_cal_imp( peaks, fitfor, nchannels, dev_pairs, coefs, uncert, &poly_coef_fcn );
}//double fit_energy_cal_poly(...)


double EnergyCal::fit_poly_from_channel_energies( const size_t ncoeffs,
                                             const std::vector<float> &channel_energies,
                                             std::vector<float> &coefs )
{
  return fit_from_channel_energies_imp( ncoeffs, channel_energies, &poly_coef_fcn, coefs );
}//fit_poly_from_channel_energies(...)



double EnergyCal::fit_full_range_fraction_from_channel_energies( const size_t ncoeffs,
                                                  const std::vector<float> &channel_energies,
                                                  std::vector<float> &coefs )
{
  if( ncoeffs >= 5 )
    throw runtime_error( "fit_full_range_fraction_from_channel_energies:"
                         " You must request less than 5 coefficients" );
  
  return fit_from_channel_energies_imp( ncoeffs, channel_energies, &frf_coef_fcn, coefs );
}//fit_full_range_fraction_from_channel_energies(...)


double EnergyCal::fit_energy_cal_iterative( const std::vector<EnergyCal::RecalPeakInfo> &peakInfo,
                                           const size_t nbin,
                                           const SpecUtils::EnergyCalType eqnType,
                                           const std::vector<bool> fitfor,
                                           std::vector<float> &startingCoefs,
                                           const std::vector<std::pair<float,float>> &devpair,
                                           std::vector<float> &coefs,
                                           std::vector<float> &coefs_uncert,
                                           std::string &warning_msg,
                                           const std::vector<float> &lower_channel_edges )
{
  coefs.clear();
  coefs_uncert.clear();
  warning_msg.clear();

  if( peakInfo.size() < 1 )
    throw runtime_error( "fit_energy_cal_iterative: no peaks specified" );

  const size_t npars = startingCoefs.size();

  if( npars < 2 )
    throw runtime_error( "fit_energy_cal_iterative: must be at least two coefficients" );

  if( fitfor.size() != npars )
    throw runtime_error( "fit_energy_cal_iterative: fitfor.size() != startingCoefs.size()" );

  switch( eqnType )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      if( !lower_channel_edges.empty() )
        throw runtime_error( "fit_energy_cal_iterative: lower_channel_edges must be empty for"
                             " Polynomial/FullRangeFraction calibrations" );
      break;

    case SpecUtils::EnergyCalType::LowerChannelEdge:
      if( npars != 2 )
        throw runtime_error( "fit_energy_cal_iterative: LowerChannelEdge requires exactly 2"
                             " coefficients (offset, gain)" );
      if( lower_channel_edges.size() != (nbin + 1) )
        throw runtime_error( "fit_energy_cal_iterative: lower_channel_edges size must equal"
                             " nchannels + 1 for LowerChannelEdge fits" );
      break;

    case SpecUtils::EnergyCalType::InvalidEquationType:
      throw runtime_error( "fit_energy_cal_iterative: invalid calibration type" );
      break;
  }//switch( eqnType )

  if( nbin < 7 )
    throw runtime_error( "fit_energy_cal_iterative: you must have at least 7 channels" );

  const size_t nfit = std::accumulate( begin(fitfor), end(fitfor), size_t(0) );

  if( nfit < 1 )
    throw runtime_error( "fit_energy_cal_iterative: must fit for at least one coefficient" );

  if( nfit > peakInfo.size() )
    throw runtime_error( "fit_energy_cal_iterative: must have at least as many peaks as"
                         " coefficients being fit for" );

  EnergyCalFitCost *cost_functor = new EnergyCalFitCost();
  cost_functor->m_eqnType = eqnType;
  cost_functor->m_nchannel = nbin;
  cost_functor->m_npars = npars;
  cost_functor->m_peakInfo = peakInfo;
  cost_functor->m_dev_pairs = devpair;
  if( eqnType == SpecUtils::EnergyCalType::LowerChannelEdge )
    cost_functor->m_baseline_lce_edges = lower_channel_edges;

  ceres::DynamicAutoDiffCostFunction<EnergyCalFitCost,sm_energy_cal_ceres_stride> *cost_function
    = new ceres::DynamicAutoDiffCostFunction<EnergyCalFitCost,sm_energy_cal_ceres_stride>(
        cost_functor, ceres::TAKE_OWNERSHIP );
  cost_function->AddParameterBlock( static_cast<int>(npars) );
  cost_function->SetNumResiduals( static_cast<int>(peakInfo.size()) );

  std::vector<double> parameters( npars, 0.0 );
  for( size_t i = 0; i < npars; ++i )
    parameters[i] = static_cast<double>( startingCoefs[i] );

  ceres::Problem problem;
  problem.AddResidualBlock( cost_function, nullptr, parameters.data() );

  // Hold fixed parameters constant via SubsetManifold.
  std::vector<int> constant_indices;
  for( size_t i = 0; i < npars; ++i )
  {
    if( !fitfor[i] )
      constant_indices.push_back( static_cast<int>(i) );
  }
  if( !constant_indices.empty() && (constant_indices.size() < npars) )
  {
    ceres::Manifold *subset = new ceres::SubsetManifold( static_cast<int>(npars), constant_indices );
    problem.SetManifold( parameters.data(), subset );
  }

  // Bound the linear/gain term (index 1) to be strictly positive when it is being fit.
  // For polynomial/FRF this preserves the Minuit2 lower-limit at 0.0; for LowerChannelEdge the
  // 'gain' multiplier must also stay positive to avoid an inverted energy axis.
  if( fitfor[1] )
    problem.SetParameterLowerBound( parameters.data(), 1, 0.0 );

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_type = ceres::TRUST_REGION;
  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
  options.use_nonmonotonic_steps = true;
  options.max_num_iterations = 200;
  options.function_tolerance = 1e-10;
  options.parameter_tolerance = 1e-12;
  options.logging_type = ceres::SILENT;
  options.minimizer_progress_to_stdout = false;

  ceres::Solver::Summary summary;
  ceres::Solve( options, &problem, &summary );

  if( summary.termination_type == ceres::FAILURE
      || summary.termination_type == ceres::USER_FAILURE )
  {
    throw runtime_error( "Fit for calibration parameters failed: " + summary.message );
  }

  for( size_t i = 0; i < npars; ++i )
  {
    if( std::isinf(parameters[i]) || std::isnan(parameters[i]) )
      throw runtime_error( "Invalid calibration parameter from fit :(" );
  }

  if( summary.termination_type == ceres::NO_CONVERGENCE )
  {
    warning_msg = "Warning: calibration coefficient fit did not fully converge; please check"
                  " the result, and revert this calibration if necessary.";
    bool fithigher = false;
    for( size_t i = 2; i < npars; ++i )
      fithigher |= fitfor[i];
    if( fithigher )
      warning_msg += " You might try not fitting for quadratic or higher terms.";
  }

  // Extract parameter uncertainties from the covariance matrix.
  std::vector<double> uncerts( npars, 0.0 );
  if( nfit < npars || nfit > 0 )
  {
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::DENSE_SVD;
    cov_options.min_reciprocal_condition_number = 1e-14;
    cov_options.null_space_rank = -1;

    ceres::Covariance covariance( cov_options );
    std::vector<std::pair<const double*,const double*>> cov_blocks;
    cov_blocks.emplace_back( parameters.data(), parameters.data() );

    if( covariance.Compute( cov_blocks, &problem ) )
    {
      std::vector<double> row_major( npars * npars, 0.0 );
      const std::vector<const double*> blocks( 1, parameters.data() );
      if( covariance.GetCovarianceMatrix( blocks, row_major.data() ) )
      {
        for( size_t i = 0; i < npars; ++i )
        {
          const double v = row_major[i * npars + i];
          uncerts[i] = (v > 0.0) ? std::sqrt(v) : 0.0;
        }
      }
    }
    else
    {
      if( !warning_msg.empty() ) warning_msg += " ";
      warning_msg += "Could not compute parameter uncertainties (covariance was rank-deficient).";
    }
  }

  // Zero uncertainties on fixed parameters (they were not fit).
  for( size_t i = 0; i < npars; ++i )
  {
    if( !fitfor[i] )
      uncerts[i] = 0.0;
  }

  coefs.resize( npars, 0.0f );
  coefs_uncert.resize( npars, 0.0f );
  for( size_t i = 0; i < npars; ++i )
  {
    coefs[i] = static_cast<float>( parameters[i] );
    coefs_uncert[i] = static_cast<float>( uncerts[i] );
  }

  // Ceres minimizes 0.5 * sum(r_i^2); chi^2 = sum(r_i^2) = 2 * final_cost.
  return 2.0 * summary.final_cost;
}//fit_energy_cal_iterative(...)



deque<shared_ptr<const PeakDef>>
EnergyCal::translatePeaksForCalibrationChange( const std::deque<std::shared_ptr<const PeakDef>> &inputPeaks,
                              const std::shared_ptr<const SpecUtils::EnergyCalibration> &old_cal,
                              const std::shared_ptr<const SpecUtils::EnergyCalibration> &new_cal )
{
//#if( PERFORM_DEVELOPER_CHECKS )
//  const double preGausArea = peak.gauss_integral(peak.lowerX(), peak.upperX() );
//  const double preContArea = peak.continuum()->offset_integral(peak.lowerX(), peak.upperX());
//#endif
  if( !old_cal || !new_cal )
    throw runtime_error( "translatePeaksForCalibrationChange: null calibration passed in" );
  
  if( old_cal->type() == SpecUtils::EnergyCalType::InvalidEquationType )
    throw runtime_error( "translatePeaksForCalibrationChange: old calibration invalid" );
  
  if( new_cal->type() == SpecUtils::EnergyCalType::InvalidEquationType )
    throw runtime_error( "translatePeaksForCalibrationChange: new calibration invalid" );
  
  if( old_cal->num_channels() < 5 )
    throw runtime_error( "translatePeaksForCalibrationChange: old calibration has less than 5 channels" );
  
  if( new_cal->num_channels() < 5 )
    throw runtime_error( "translatePeaksForCalibrationChange: new calibration has less than 5 channels" );
  
  deque<shared_ptr<const PeakDef>> answer;
  
  // First go through and map peaks that share a continuum
  map<shared_ptr<const PeakContinuum>,deque<shared_ptr<const PeakDef>>> conts_to_peaks;
  for( const auto &p : inputPeaks )
    conts_to_peaks[p->continuum()].push_back( p );
  
  // Incase any peaks use an external continuum, we will only create a single copy of it and use
  //  for all peaks.
  std::shared_ptr<SpecUtils::Measurement> data_continuum;
  
  for( auto roi : conts_to_peaks )
  {
    const shared_ptr<const PeakContinuum> &oldcont = roi.first;
    const deque<shared_ptr<const PeakDef>> &peaks = roi.second;
    
    // Lets first translate the continuum
    assert( oldcont );
    auto newcont = make_shared<PeakContinuum>();
    newcont->setType( oldcont->type() );
    
    float strech = 1.0f;
    if( oldcont->energyRangeDefined() )
    {
      const float oldLowEnergy = static_cast<float>( oldcont->lowerEnergy() );
      const float oldlowbin = old_cal->channel_for_energy( oldLowEnergy );
      const float new_lowenergy = new_cal->energy_for_channel( oldlowbin );
      
      const float oldHighEnergy = static_cast<float>( oldcont->upperEnergy() );
      const float oldhighbin = old_cal->channel_for_energy( oldHighEnergy );
      const float new_highenergy = new_cal->energy_for_channel( oldhighbin );
        
      strech = (new_highenergy - new_lowenergy) / (oldHighEnergy - oldLowEnergy);
      newcont->setRange( new_lowenergy, new_highenergy );
    }//if( peak.continuum().energyRangeDefined() )
      
    if( oldcont->isPolynomial() )
    {
      const double oldref = oldcont->referenceEnergy();
      const float oldrefbin = old_cal->channel_for_energy( oldref );
      const float newref = new_cal->energy_for_channel( oldrefbin );
        
      if( !oldcont->energyRangeDefined() )
        strech = static_cast<float>( newref / oldref );
        
      if( IsNan(strech) || IsInf(strech) )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__,
                              "Found an invalid stretch value when calculated "
                              "by the reference energy" );
#endif
        //strech = newMean / oldMean;
      }//if( IsNan(strech) || IsInf(strech) )
 
      vector<double> vars = oldcont->parameters();
      vector<double> uncerts = oldcont->uncertainties();

      for( size_t i = 0; i < vars.size(); ++i )
      {
        vars[i] = vars[i] / std::pow( strech, static_cast<float>(i+1.0f) );
        uncerts[i] = uncerts[i] / std::pow( strech, static_cast<float>(i+1.0f) );
      }//for( size_t i = 0; i < contvars.size(); ++i )
        
      newcont->setParameters( newref, vars, uncerts );
    }else if( oldcont->externalContinuum() )
    {
      if( !data_continuum )
      {
        data_continuum = std::make_shared<SpecUtils::Measurement>( *oldcont->externalContinuum() );
        data_continuum->set_energy_calibration( new_cal );
      }
      
      newcont->setExternalContinuum( data_continuum );
    }//if( peak.continuum().isPolynomial() ) / else if(
      
    
    // Now that we have he continuum modified, lets modify all the peaks for this ROI.
    for( const shared_ptr<const PeakDef> &oldpeak : peaks )
    {
      assert( oldpeak );
      shared_ptr<PeakDef> newpeak = make_shared<PeakDef>( *oldpeak );
      newpeak->setContinuum( newcont );
      
      if( !oldpeak->gausPeak() )
      {
        const float oldMean = static_cast<float>( oldpeak->mean() );
        const float meanbin = old_cal->channel_for_energy( oldMean );
        const float newMean = new_cal->energy_for_channel(  meanbin );
        newpeak->set_coefficient( newMean, PeakDef::Mean );
        
        // The new continuums range should have already been set, but if not could do it as:
        //const float oldlow = static_cast<float>( oldpeak->lowerX() );
        //const float oldhigh = static_cast<float>( oldpeak->upperX() );
        //const float lowbin = old_cal->channel_for_energy( oldlow );
        //const float highbin = old_cal->channel_for_energy( oldhigh );
        //const float newLower = new_cal->energy_for_channel( lowbin );
        //const float newUpper = new_cal->energy_for_channel( highbin );
        //newcont->setRange( newLower, newUpper );
        
        answer.push_back( newpeak );
        continue;
      }//if( !peak.gausPeak() )
      
      
      const float oldMean = static_cast<float>( oldpeak->mean() );
      const float oldSigma = static_cast<float>( oldpeak->sigma() );
      const float meanbin = old_cal->channel_for_energy( oldMean );
      const float newMean = new_cal->energy_for_channel( meanbin );
      
      const float oldneg2sigmabin = old_cal->channel_for_energy( oldMean - 2.0*oldSigma );
      const float oldpos2sigmabin = old_cal->channel_for_energy( oldMean + 2.0*oldSigma );
      const float newneg2sigma = new_cal->energy_for_channel( oldneg2sigmabin );
      const float newpos2sigma = new_cal->energy_for_channel( oldpos2sigmabin );
      
      float strech = 0.25f*(newpos2sigma - newneg2sigma) / oldSigma;
      
      if( IsNan(strech) || IsInf(strech) )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        const char * msg = "Found an invalid stretch value when claculated from the"
        " mean";
        log_developer_error( __func__, msg );
#endif
        //throw runtime_error( "translatePeaksForCalibrationChange: found an invalid stretch" );
        strech = 1.0f; //JIC - dont expect to happen unless `oldSigma` is 0, or the input peak is invalid.
      }//if( IsNan(strech) || IsInf(strech) )
      
#if( PERFORM_DEVELOPER_CHECKS )
      const float newbin = new_cal->channel_for_energy( newMean );
      if( fabs(newbin - meanbin) > 0.025 )  //0.025 arbitrary
      {
        stringstream msg;
        msg.precision( 9 );
        msg << "When recalibrating from coefs={";
        for( size_t i = 0; i < old_cal->coefficients().size(); ++i )
          msg << (i ? ", " : " ") << old_cal->coefficients()[i];
        msg << " }";
        if( old_cal->deviation_pairs().size() )
        {
          msg << ", devpairs={";
          for( size_t i = 0; i < old_cal->deviation_pairs().size(); ++i )
            msg << (i ? ", {" : " {") << old_cal->deviation_pairs()[i].first
            << "," << old_cal->deviation_pairs()[i].second << "}";
          msg << "}";
        }
        
        msg << " to {";
        for( size_t i = 0; i < new_cal->coefficients().size(); ++i )
          msg << (i ? ", " : " ") << new_cal->coefficients()[i];
        msg << " }";
        if( new_cal->deviation_pairs().size() )
        {
          msg << ", devpairs={";
          for( size_t i = 0; i < new_cal->deviation_pairs().size(); ++i )
            msg << (i ? ", {" : " {") << new_cal->deviation_pairs()[i].first << ","
            << new_cal->deviation_pairs()[i].second << "}";
          msg << "}";
        }
        msg << ", peak at mean=" << oldMean << ", moved to mean=" << newMean
        << ", but some error caused it to move from bin " << meanbin << " to "
        << newbin << ", which shouldnt have happend (should have stayed same "
        << "bin number).";
        
        log_developer_error( __func__, msg.str().c_str() );
      }//if( fabs(newbin - oldbin) > 0.025 )
#endif
      
      newpeak->setMean( newMean );
      newpeak->setSigma( strech * oldSigma );
      newpeak->setMeanUncert( strech * oldpeak->meanUncert() );
      newpeak->setSigmaUncert( strech * oldpeak->sigmaUncert() );
      
      
      answer.push_back( newpeak );
    }//for( loop over old peaks )
  }//for( auto roi : conts_to_peaks )
  
  return answer;
}//translatePeakForCalibrationChange(...)


shared_ptr<const SpecUtils::EnergyCalibration>
EnergyCal::propogate_energy_cal_change( const shared_ptr<const SpecUtils::EnergyCalibration> &orig_cal,
                             const shared_ptr<const SpecUtils::EnergyCalibration> &new_cal,
                             const shared_ptr<const SpecUtils::EnergyCalibration> &other_cal )
{
  using namespace SpecUtils;
  
  if( !orig_cal || !new_cal || !other_cal
     || !orig_cal->valid() || !new_cal->valid() || !other_cal->valid() )
    throw runtime_error( "EnergyCal::propogate_energy_cal_change invalid input" );
  
  if( orig_cal == new_cal )
    return other_cal;
  
  auto answer = make_shared<EnergyCalibration>();
  
  const vector<float> &prev_disp_coefs = orig_cal->coefficients();
  const vector<float> &new_disp_coefs = new_cal->coefficients();
  const vector<float> &other_coeffs = other_cal->coefficients();
  
  const size_t orig_num_channel = orig_cal->num_channels();
  const size_t new_num_channel = orig_cal->num_channels();
  const size_t other_num_channel = other_cal->num_channels();
  
  // \TODO: we currently arent using deviation pairs when converting between channel number and
  //        energy in this function.  This isnt actually correct; the aprehentions I have are:
  //        1) I havent tested that with deviation is numerically stable enough.
  //        2) What if its only deviation pairs that have changed, do we then want to correct for
  //           this using the coefficients?
  //        3) I'm not a hundred percent certain we actually want to correct for deviation pairs;
  //           needs more thought;
  const vector<pair<float,float>> prev_disp_devs; // = orig_cal->deviation_pairs();
  const vector<pair<float,float>> new_disp_devs;  // = new_cal->deviation_pairs();
  const vector<pair<float,float>> other_devs;     // = other_cal->deviation_pairs();
  
  // Deal with the easy case of other_cal being lower channel energies.
  if( other_cal->type() == EnergyCalType::LowerChannelEdge )
  {
    const size_t nchannel = other_cal->num_channels();
    const shared_ptr<const vector<float>> &old_lower_ptr = other_cal->channel_energies();
    assert( old_lower_ptr && !old_lower_ptr->empty() );
    
    const vector<float> &old_lower = *old_lower_ptr;
    assert( old_lower.size() >= (nchannel + 1) ); //actually should always be equal
    if( nchannel >= old_lower.size() )
      throw runtime_error( "EnergyCal::propogate_energy_cal_change: really unexpected programing error" );
    
    vector<float> new_lower( nchannel + 1 );
    for( size_t i = 0; i <= nchannel; ++i )
    {
      const double equiv_channel = orig_cal->channel_for_energy( old_lower[i] );
      new_lower[i] = new_cal->energy_for_channel( equiv_channel );
    }
    
    answer->set_lower_channel_energy( nchannel, std::move(new_lower) );
    return answer;
  }//if( other_cal->type() == EnergyCalType::LowerChannelEdge )

  
  // LowerChannelEdge foreground propagating to a poly / FRF measurement: the displayed LCE
  // change is always a 2-parameter affine of the loaded edges (offset + gain — that's what the
  // UI exposes), so we don't need the full channel-matching dance below. Just shift the other
  // cal's offset by the same energy delta and multiply its linear term by the same gain ratio.
  // Higher-order coefficients are left alone (they represent intrinsic non-linearity of the
  // other detector that shouldn't change when the displayed detector's offset/gain is nudged).
  if( orig_cal->type() == EnergyCalType::LowerChannelEdge )
  {
    const auto orig_edges = orig_cal->channel_energies();
    const auto new_edges  = (new_cal->type() == EnergyCalType::LowerChannelEdge)
                              ? new_cal->channel_energies()
                              : shared_ptr<const vector<float>>{};
    if( !orig_edges || orig_edges->size() < 2 || !new_edges || new_edges->size() < 2
        || (orig_edges->size() != new_edges->size()) )
      throw runtime_error( "EnergyCal::propogate_energy_cal_change: LowerChannelEdge edges missing or size mismatch" );

    const double orig_span = static_cast<double>(orig_edges->back())
                             - static_cast<double>(orig_edges->front());
    const double new_span = static_cast<double>(new_edges->back())
                            - static_cast<double>(new_edges->front());
    const double offset_delta = static_cast<double>(new_edges->front())
                                - static_cast<double>(orig_edges->front());
    const double gain_factor = (std::fabs(orig_span) > 1e-6) ? (new_span / orig_span) : 1.0;

    vector<float> new_other_coefs = other_cal->coefficients();
    if( new_other_coefs.empty() )
      new_other_coefs.push_back( 0.0f );
    new_other_coefs[0] = static_cast<float>( static_cast<double>(new_other_coefs[0]) + offset_delta );
    if( new_other_coefs.size() >= 2 )
      new_other_coefs[1] = static_cast<float>( static_cast<double>(new_other_coefs[1]) * gain_factor );

    const auto &dev_pairs = other_cal->deviation_pairs();
    switch( other_cal->type() )
    {
      case EnergyCalType::Polynomial:
      case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        answer->set_polynomial( other_num_channel, new_other_coefs, dev_pairs );
        break;

      case EnergyCalType::FullRangeFraction:
        answer->set_full_range_fraction( other_num_channel, new_other_coefs, dev_pairs );
        break;

      case EnergyCalType::LowerChannelEdge:
      case EnergyCalType::InvalidEquationType:
        assert( 0 );
        break;
    }//switch( other_cal->type() )

    return answer;
  }//if( orig_cal is LowerChannelEdge )

  // At this point both `orig_cal` and `new_cal` are poly / FRF, and `other_cal` is poly / FRF
  // too (LowerChannelEdge other_cal returned above).
  assert( (orig_cal->type() == EnergyCalType::FullRangeFraction)
          || (orig_cal->type() == EnergyCalType::Polynomial)
          || (orig_cal->type() == EnergyCalType::UnspecifiedUsingDefaultPolynomial) );
  assert( (new_cal->type() == EnergyCalType::FullRangeFraction)
          || (new_cal->type() == EnergyCalType::Polynomial)
          || (new_cal->type() == EnergyCalType::UnspecifiedUsingDefaultPolynomial) );
  assert( (other_cal->type() == EnergyCalType::FullRangeFraction)
          || (other_cal->type() == EnergyCalType::Polynomial)
          || (other_cal->type() == EnergyCalType::UnspecifiedUsingDefaultPolynomial) );

  const double accuracy = 0.00001;
  const size_t order = std::max( other_coeffs.size(),
                                 std::max(prev_disp_coefs.size(), new_disp_coefs.size()) );
  
  vector<pair<double,double>> channels_energies;  //this gives <channel number,energy it should be>
  for( size_t i = 0; i < order; ++i )
  {
    const size_t display_channel = ((order - i - 1) * orig_num_channel) / (order - 1);
     
    double old_disp_energy = std::numeric_limits<double>::quiet_NaN(),
           new_disp_energy = std::numeric_limits<double>::quiet_NaN(),
           other_channel = std::numeric_limits<double>::quiet_NaN();
    
    switch( orig_cal->type() )
    {
      case EnergyCalType::FullRangeFraction:
      {
        old_disp_energy = fullrangefraction_energy( display_channel, prev_disp_coefs, orig_num_channel, prev_disp_devs );
        new_disp_energy = fullrangefraction_energy( display_channel, new_disp_coefs, new_num_channel, new_disp_devs );
        break;
      }//case orig_cal was FRF

      case EnergyCalType::Polynomial:
      case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      {
        old_disp_energy = polynomial_energy( display_channel, prev_disp_coefs, prev_disp_devs );
        new_disp_energy = polynomial_energy( display_channel, new_disp_coefs, new_disp_devs );
        break;
      }//case: orig_cal was Poly

      case EnergyCalType::LowerChannelEdge:
      case EnergyCalType::InvalidEquationType:
        // LowerChannelEdge handled by the early-return block above; InvalidEquationType
        // rejected by the up-front validity check.
        assert( 0 );
        break;
    }//switch( orig_cal->type() )
    
    assert( !IsNan(old_disp_energy) && !IsNan(new_disp_energy) );
    
    double check_energy = std::numeric_limits<double>::quiet_NaN();  //Just for development test
    switch( other_cal->type() )
    {
      case EnergyCalType::Polynomial:
      case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        other_channel = find_polynomial_channel( old_disp_energy, other_coeffs, other_num_channel,
                                                 other_devs, accuracy );
        check_energy = polynomial_energy( other_channel, other_coeffs, other_devs );
        break;
        
      case EnergyCalType::FullRangeFraction:
        other_channel = find_fullrangefraction_channel( old_disp_energy, other_coeffs,
                                                       other_num_channel, other_devs, accuracy );
        check_energy = fullrangefraction_energy( other_channel, other_coeffs, other_num_channel,
                                                 other_devs );
        break;
        
      case EnergyCalType::LowerChannelEdge:
      case EnergyCalType::InvalidEquationType:
        assert( 0 );
        break;
    }//switch( other_cal->type() )
    
    assert( !IsNan(other_channel) && !IsNan(check_energy) );
    
    channels_energies.push_back( {other_channel, new_disp_energy} );
  
#if( PERFORM_DEVELOPER_CHECKS )
    const double diff = fabs(check_energy - old_disp_energy);
    const double max_energy = std::max( fabs(check_energy), fabs(old_disp_energy) );
    
    if( (diff > (max_energy * 1.0E-5)) && (diff > 0.00001) )
    {
      char buffer[256];
      snprintf( buffer, sizeof(buffer), "Found case going from energy-->channel-->energy"
               " gave seconf energy too different than initial energy by %f with check_energy=%f"
               " and old_disp_energy=%f", diff, check_energy, old_disp_energy );
      log_developer_error( __func__, buffer );
    }//if( diff is larger than we wanted )

    assert( diff <= max_energy*1.0E-6 || (diff < 0.00001) );
#endif
  }//for( size_t i = 0; i < order; ++i )
  
  const auto &dev_pairs = other_cal->deviation_pairs();
  switch( other_cal->type() )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    {
      const vector<float> new_other_coefs = EnergyCal::fit_for_poly_coefs( channels_energies, static_cast<int>(order) );
      answer->set_polynomial( other_num_channel, new_other_coefs, dev_pairs );
      break;
    }
      
    case SpecUtils::EnergyCalType::FullRangeFraction:
    {
      const vector<float> new_other_coefs = EnergyCal::fit_for_fullrangefraction_coefs( channels_energies,
                                                                        other_num_channel, static_cast<int>(order) );
      answer->set_full_range_fraction( other_num_channel, new_other_coefs, dev_pairs );
      break;
    }
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
      break;
  }//switch( other_cal->type() )
  
  assert( answer->valid() );
  
  return answer;
}//propogate_energy_cal_change(...)


std::vector<float> EnergyCal::fit_for_poly_coefs( const std::vector<std::pair<double,double>> &channels_energies,
                            const int poly_terms )
{
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation using Eigen SVD for numerical stability
  
  const size_t npoints = channels_energies.size();
  Eigen::MatrixX<double> A( npoints, poly_terms );
  Eigen::VectorX<double> b( npoints );
  
  for( size_t row = 0; row < npoints; ++row )
  {
    b(row) = channels_energies[row].second;
    for( int col = 0; col < poly_terms; ++col )
      A(row,col) = std::pow( channels_energies[row].first, double(col) );
  }//for( int col = 0; col < poly_terms; ++col )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  
  vector<float> poly_coeffs( poly_terms );
  for( int coef = 0; coef < poly_terms; ++coef )
    poly_coeffs[coef] = static_cast<float>( a(coef) );
  
  return poly_coeffs;
}//void fit_for_poly_coefs(...)


std::vector<float> EnergyCal::fit_for_fullrangefraction_coefs( const std::vector<std::pair<double,double>> &channels_energies,
                            const size_t nchannels, const int nterms )
{
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation using Eigen SVD for numerical stability
  
  const int polyterms = std::min( nterms, 5 );
  const size_t npoints = channels_energies.size();
  
  Eigen::MatrixX<double> A( npoints, polyterms );
  Eigen::VectorX<double> b( npoints );
  
  for( size_t row = 0; row < npoints; ++row )
  {
    const double x = channels_energies[row].first / nchannels;
    
    b(row) = channels_energies[row].second;
    for( int col = 0; col < std::min(4,polyterms); ++col )
      A(row,col) = std::pow( x, static_cast<double>(col) );
    if( polyterms > 4 )
      A(row,4) = 1.0 / (1.0 + 60.0*x);
  }//for( int col = 0; col < poly_terms; ++col )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  
  vector<float> frf_coeffs( polyterms );
  for( int coef = 0; coef < polyterms; ++coef )
    frf_coeffs[coef] = static_cast<float>( a(coef) );
  
  return frf_coeffs;
}//void fit_for_fullrangefraction_coefs(...)


