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
#include <algorithm>
#include <stdexcept>

#include "Eigen/Dense"


//Ceres includes; note: "ceres/ceres.h" must come before
//  "InterSpec/RelActCalcAuto_EnergyCal_imp.hpp", so the real ceres::Jet is used
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


/** Solves the (already row-weighted) linear least squares problem `A*x = b` via SVD of the
 column-equilibrated matrix, optionally also computing the coefficient variances.

 Raw polynomial basis columns (e.g., channel^k, with channel up to ~64k) differ by many orders of
 magnitude, which loses accuracy in the SVD - and did so catastrophically in the coefficient
 uncertainties, which were previously computed as inverse(A^T*A) (squaring the condition number).
 Scaling column j of A by s_j is equivalent to fitting x_j/s_j, so the solution is un-scaled
 before returning; in exact arithmetic the answer is unchanged.  The covariance is computed from
 the SVD factors as the pseudo-inverse V*diag(1/sigma_i^2)*V^T, with singular values below
 max(sigma)*1e-10 contributing zero variance, rather than near-singular directions blowing the
 variances up.
 */
Eigen::VectorXd solve_lls_equilibrated( Eigen::MatrixX<double> A,
                                        const Eigen::VectorX<double> &b,
                                        Eigen::VectorXd * const coef_variances )
{
  const Eigen::Index ncols = A.cols();

  Eigen::VectorXd col_scale( ncols );
  for( Eigen::Index col = 0; col < ncols; ++col )
  {
    const double mag = A.col(col).cwiseAbs().maxCoeff();
    col_scale(col) = (mag > 0.0) ? (1.0 / mag) : 1.0;
    A.col(col) *= col_scale(col);
  }//for( loop over columns of A )

#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif

  const Eigen::VectorXd scaled_answer = svd.solve(b);

  if( coef_variances )
  {
    const Eigen::VectorXd &sigmas = svd.singularValues();
    const double sigma_min_allowed = (sigmas.size() ? sigmas(0) : 0.0) * 1.0E-10;

    Eigen::VectorXd inv_sigma2 = Eigen::VectorXd::Zero( sigmas.size() );
    for( Eigen::Index i = 0; i < sigmas.size(); ++i )
      inv_sigma2(i) = (sigmas(i) > sigma_min_allowed) ? 1.0/(sigmas(i)*sigmas(i)) : 0.0;

    const Eigen::MatrixX<double> &V = svd.matrixV();

    coef_variances->resize( ncols );
    for( Eigen::Index j = 0; j < ncols; ++j )
    {
      double scaled_var = 0.0;
      for( Eigen::Index i = 0; i < sigmas.size(); ++i )
        scaled_var += V(j,i) * V(j,i) * inv_sigma2(i);

      (*coef_variances)(j) = col_scale(j) * col_scale(j) * scaled_var;
    }//for( loop over coefficients )
  }//if( coef_variances )

  return col_scale.cwiseProduct( scaled_answer );
}//solve_lls_equilibrated(...)


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
  
  Eigen::VectorXd coef_variances;
  const Eigen::VectorXd a = solve_lls_equilibrated( A, b, &coef_variances );

  coefs.resize( fitfor.size(), 0.0 );
  coefs_uncert.resize( fitfor.size(), 0.0 );

  for( size_t col = 0, coef_index = 0; coef_index < fitfor.size(); ++coef_index )
  {
    if( fitfor[coef_index] )
    {
      assert( col < nparsfit );
      coefs[coef_index] = static_cast<float>( a(col) );
      coefs_uncert[coef_index] = static_cast<float>( std::sqrt( coef_variances(col) ) );
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

  const Eigen::VectorXd a = solve_lls_equilibrated( A, b, nullptr );

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



/** Cost functor for #EnergyCal::fit_energy_cal_ceres.

 One parameter block, laid out as [coefficients..., deviation-pair offsets...]; parameters not
 being fit are held constant with a ceres::SubsetManifold.  Polynomial coefficients are scaled by
 nchannel^order internally, so all parameters (and their gradients) are of roughly keV magnitude.

 The residuals are one per peak, `(E_predicted(channel) - photopeakEnergy)/sigma` - the same
 model and weighting as the linear fits - plus two soft hinge residuals that punish the
 calibration becoming non-monotonic (replacing the hard chi2 cliff of the old Minuit based
 implementation, which gave the minimizer no gradient to recover along).
 */
struct EnergyCalCeresCostFcn
{
  const EnergyCal::EnergyCalCeresFitSetup &m_setup;

  /** If no deviation pair offsets are being fit, the (fixed) deviation pairs are subtracted from
   the target peak energies up-front - exactly like the linear fits do - instead of evaluating
   the spline inside the fit. */
  bool m_fit_any_dev_offsets;

  /** Multiplicative scale from parameter-space to coefficient-space. */
  std::vector<double> m_par_scales;

  std::vector<double> m_channels, m_sigmas, m_targets;

  /** Per-peak {original energy, fraction of original energy range}; only for LowerChannelEdge,
   whose predicted energy is `orig_energy + offset + gain*frac` (see
   #EnergyCal::adjust_lower_channel_energy_cal). */
  /** For a LowerChannelEdge fit, each peaks ORIGINAL energy (E_orig); the model is
   E_new = offset + gain*E_orig (offset keV, gain dimensionless, nominal 1.0). */
  std::vector<double> m_lower_chan_orig;

  /** The original energy range of a LowerChannelEdge calibration; its channel energies stay
   monotonic exactly when `gain > -range`. */
  double m_lower_chan_range;

  /** How hard to punish each keV of non-monotonicity at the spectrum ends. */
  static constexpr double sm_mono_hinge_weight = 100.0;

  EnergyCalCeresCostFcn( const std::vector<EnergyCal::RecalPeakInfo> &peaks,
                         const EnergyCal::EnergyCalCeresFitSetup &setup )
  : m_setup( setup ),
    m_fit_any_dev_offsets( false ),
    m_lower_chan_range( 0.0 )
  {
    for( const bool fit : setup.fit_dev_pair_offsets )
      m_fit_any_dev_offsets = (m_fit_any_dev_offsets || fit);

    const size_t ncoefs = setup.fitfor.size();
    m_par_scales.resize( ncoefs + setup.dev_pairs.size(), 1.0 );

    switch( setup.cal_type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        for( size_t i = 0; i < ncoefs; ++i )
          m_par_scales[i] = 1.0 / std::pow( static_cast<double>(setup.num_channels),
                                            static_cast<double>(i) );
        break;

      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        break;  //coefficients are already keV magnitude

      case SpecUtils::EnergyCalType::InvalidEquationType:
        assert( 0 );
        break;
    }//switch( setup.cal_type )

    if( setup.cal_type == SpecUtils::EnergyCalType::LowerChannelEdge )
    {
      assert( setup.lower_channel_energies.size() >= 2 );
      m_lower_chan_range = setup.lower_channel_energies.back() - setup.lower_channel_energies.front();
    }

    const std::vector<RelActCalcAutoImp::CubicSplineNodeT<double>> no_spline;

    for( const EnergyCal::RecalPeakInfo &peak : peaks )
    {
      m_channels.push_back( peak.peakMeanBinNumber );

      const double sigma = peak.photopeakEnergy * peak.peakMeanUncert
                           / (std::max)( peak.peakMean, 1.0 );
      m_sigmas.push_back( (std::max)( fabs(sigma), 1.0e-6 ) );

      double target = peak.photopeakEnergy;
      if( !m_fit_any_dev_offsets && !setup.dev_pairs.empty() )
        target -= SpecUtils::correction_due_to_dev_pairs( peak.photopeakEnergy, setup.dev_pairs );
      m_targets.push_back( target );

      if( setup.cal_type == SpecUtils::EnergyCalType::LowerChannelEdge )
      {
        // E_new = offset + gain*E_orig; store the peaks original energy
        const double orig_energy = RelActCalcAutoImp::lowerchannel_energy(
                              peak.peakMeanBinNumber, setup.lower_channel_energies, no_spline );
        m_lower_chan_orig.push_back( orig_energy );
      }
    }//for( const EnergyCal::RecalPeakInfo &peak : peaks )
  }//EnergyCalCeresCostFcn constructor

  size_t num_coefs() const { return m_setup.fitfor.size(); }
  size_t num_parameters() const { return m_setup.fitfor.size() + m_setup.dev_pairs.size(); }
  size_t num_residuals() const { return m_channels.size() + 2; }

  template<typename T>
  bool operator()( T const *const *parameters, T *residuals ) const
  {
    try
    {
      const T * const pars = parameters[0];
      const size_t ncoefs = num_coefs();
      const size_t npeaks = m_channels.size();

      std::vector<T> coefs( ncoefs );
      for( size_t i = 0; i < ncoefs; ++i )
        coefs[i] = pars[i] * m_par_scales[i];

      // The deviation pair spline; only built (and only needed) when offsets are being fit
      std::vector<RelActCalcAutoImp::CubicSplineNodeT<T>> spline;
      if( m_fit_any_dev_offsets && (m_setup.dev_pairs.size() > 1) )
      {
        std::vector<std::pair<double,T>> pairs;
        pairs.reserve( m_setup.dev_pairs.size() );
        for( size_t i = 0; i < m_setup.dev_pairs.size(); ++i )
          pairs.emplace_back( static_cast<double>(m_setup.dev_pairs[i].first), pars[ncoefs + i] );

        spline = RelActCalcAutoImp::create_cubic_spline( pairs );
      }//if( fitting deviation pair offsets )

      // Predicted (observed) energy at a channel; for peaks the per-peak lower-channel values
      //  are precomputed, but the monotonicity checks below need arbitrary channels, so
      //  lower-channel is special-cased where this lambda is used.
      const auto base_energy = [&]( const double channel ) -> T {
        switch( m_setup.cal_type )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            return RelActCalcAutoImp::polynomial_energy( T(channel), coefs, spline );

          case SpecUtils::EnergyCalType::FullRangeFraction:
            return RelActCalcAutoImp::fullrangefraction_energy( T(channel), coefs,
                                                            m_setup.num_channels, spline );

          case SpecUtils::EnergyCalType::LowerChannelEdge:
          case SpecUtils::EnergyCalType::InvalidEquationType:
            assert( 0 );  //lower-channel handled separately; invalid never gets here
            break;
        }//switch( m_setup.cal_type )

        return T(0.0);
      };//base_energy lambda

      const bool lower_channel = (m_setup.cal_type == SpecUtils::EnergyCalType::LowerChannelEdge);

      for( size_t i = 0; i < npeaks; ++i )
      {
        T pred;
        if( lower_channel )
          pred = coefs[0] + coefs[1]*m_lower_chan_orig[i];  //offset + gain*E_orig
        else
          pred = base_energy( m_channels[i] );

        residuals[i] = (pred - m_targets[i]) / m_sigmas[i];
      }//for( loop over peaks )

      // Soft punishments for the calibration becoming non-monotonic at the spectrum ends (for a
      //  lower-channel calibration the monotonicity condition is uniform: gain > -range)
      T lower_diff, upper_diff;
      if( lower_channel )
      {
        //Monotonic iff the new energy range (gain*orig_range) is positive, i.e., gain > 0
        lower_diff = upper_diff = coefs[1] * T(m_lower_chan_range);
      }else
      {
        const double nchan = static_cast<double>( m_setup.num_channels );
        lower_diff = base_energy(1.0) - base_energy(0.0);
        upper_diff = base_energy(nchan) - base_energy(nchan - 1.0);
      }

      residuals[npeaks] = (lower_diff < T(0.0)) ? T(-sm_mono_hinge_weight)*lower_diff : T(0.0);
      residuals[npeaks + 1] = (upper_diff < T(0.0)) ? T(-sm_mono_hinge_weight)*upper_diff : T(0.0);
    }catch( std::exception & )
    {
      return false;  //Ceres treats this as an invalid step, and tries a different one
    }

    return true;
  }//operator()
};//struct EnergyCalCeresCostFcn

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


std::shared_ptr<const SpecUtils::EnergyCalibration>
EnergyCal::adjust_lower_channel_energy_cal( const std::shared_ptr<const SpecUtils::EnergyCalibration> &orig,
                                            const double offset, const double gain )
{
  using SpecUtils::EnergyCalibration;

  if( !orig || !orig->valid() || (orig->type() != SpecUtils::EnergyCalType::LowerChannelEdge) )
    throw runtime_error( "adjust_lower_channel_energy_cal: input must be a valid lower channel"
                         " energy calibration" );

  //A non-positive gain would make the (increasing) original energies non-increasing
  if( gain <= 0.0 )
    throw runtime_error( "adjust_lower_channel_energy_cal: gain must be positive" );

  const shared_ptr<const vector<float>> &orig_energies = orig->channel_energies();
  assert( orig_energies && (orig_energies->size() > 1) );
  if( !orig_energies || (orig_energies->size() < 2) )
    throw runtime_error( "adjust_lower_channel_energy_cal: invalid input channel energies" );

  const size_t nchannel = orig->num_channels();

  //  E_new[i] = offset + gain*E_orig[i]   (offset in keV, gain dimensionless, nominal 1.0)
  vector<float> new_energies( *orig_energies );
  for( size_t i = 0; i < new_energies.size(); ++i )
    new_energies[i] = static_cast<float>( offset + gain*new_energies[i] );

  auto answer = make_shared<EnergyCalibration>();
  //set_lower_channel_energy(...) throws if the result isnt monotonically increasing
  answer->set_lower_channel_energy( nchannel, std::move(new_energies) );

  return answer;
}//adjust_lower_channel_energy_cal(...)


double EnergyCal::fit_energy_cal_lower_channel( const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                              const std::shared_ptr<const SpecUtils::EnergyCalibration> &orig_cal,
                                                const std::vector<bool> &fitfor,
                                                std::vector<float> &coefs,
                                                std::vector<float> &coefs_uncert )
{
  if( !orig_cal || !orig_cal->valid()
      || (orig_cal->type() != SpecUtils::EnergyCalType::LowerChannelEdge) )
    throw runtime_error( "fit_energy_cal_lower_channel: input must be a valid lower channel"
                         " energy calibration" );

  if( fitfor.size() != 2 )
    throw runtime_error( "fit_energy_cal_lower_channel: fitfor must have exactly two ({offset,"
                         " gain}) entries" );

  const shared_ptr<const vector<float>> &orig_energies = orig_cal->channel_energies();
  assert( orig_energies && (orig_energies->size() > 1) );
  if( !orig_energies || (orig_energies->size() < 2) )
    throw runtime_error( "fit_energy_cal_lower_channel: invalid input channel energies" );

  const size_t nchannels = orig_cal->num_channels();
  const double range = orig_energies->back() - orig_energies->front();
  if( range <= 0.0 )
    throw runtime_error( "fit_energy_cal_lower_channel: invalid original energy range" );

  // We model the true peak energies as:
  //   E_true = offset + gain*E_orig(channel)   (offset keV, gain dimensionless, nominal 1.0)
  //  which is a plain 2-parameter linear fit ({offset, gain}), so we reuse the linear-least-
  //  squares fitter directly: coef 0 = offset (basis 1), coef 1 = gain (basis E_orig(channel)).
  const auto basis = [orig_cal]( size_t order, double channel, size_t nchan ) -> double {
    switch( order )
    {
      case 0: return 1.0;
      case 1: return orig_cal->energy_for_channel( channel );
    }
    assert( 0 );
    return 0.0;
  };

  const vector<pair<float,float>> dev_pairs;  //lower channel energy cals dont have deviation pairs
  const double chi2 = fit_energy_cal_imp( peakinfos, fitfor, nchannels, dev_pairs,
                                          coefs, coefs_uncert, basis );

  assert( (coefs.size() == 2) && (coefs_uncert.size() == 2) );

  return chi2;
}//fit_energy_cal_lower_channel(...)


EnergyCal::EnergyCalCeresFitResult EnergyCal::fit_energy_cal_ceres(
                                          const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                                          const EnergyCal::EnergyCalCeresFitSetup &setup_in )
{
  using SpecUtils::EnergyCalType;

  EnergyCalCeresFitSetup setup = setup_in;  //we may prepend an implicit {0,0} deviation pair

  const size_t ncoefs = setup.fitfor.size();

  if( peakinfos.empty() )
    throw runtime_error( "fit_energy_cal_ceres: no peaks specified" );

  if( setup.starting_coefs.size() != ncoefs )
    throw runtime_error( "fit_energy_cal_ceres: fitfor and starting_coefs must be the same size" );

  if( !setup.fit_dev_pair_offsets.empty()
      && (setup.fit_dev_pair_offsets.size() != setup.dev_pairs.size()) )
    throw runtime_error( "fit_energy_cal_ceres: fit_dev_pair_offsets must be empty, or the same"
                         " size as dev_pairs" );

  if( setup.num_channels < 2 )
    throw runtime_error( "fit_energy_cal_ceres: must have at least two channels" );

  switch( setup.cal_type )
  {
    case EnergyCalType::Polynomial:
    case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      if( ncoefs < 2 )
        throw runtime_error( "fit_energy_cal_ceres: must have at least two coefficients" );
      break;

    case EnergyCalType::FullRangeFraction:
      if( (ncoefs < 2) || (ncoefs > 5) )
        throw runtime_error( "fit_energy_cal_ceres: full range fraction must have two to five"
                             " coefficients" );
      break;

    case EnergyCalType::LowerChannelEdge:
      if( ncoefs != 2 )
        throw runtime_error( "fit_energy_cal_ceres: lower channel energy must have exactly the"
                             " two {offset, gain} coefficients" );
      if( !setup.dev_pairs.empty() )
        throw runtime_error( "fit_energy_cal_ceres: lower channel energy calibrations cant have"
                             " deviation pairs" );
      if( setup.lower_channel_energies.size() < (setup.num_channels + 1) )
        throw runtime_error( "fit_energy_cal_ceres: not enough lower channel energies" );
      break;

    case EnergyCalType::InvalidEquationType:
      throw runtime_error( "fit_energy_cal_ceres: invalid calibration type" );
  }//switch( setup.cal_type )

  bool fit_any_dev = false;
  size_t num_fit_coefs = 0, num_fit_dev = 0;
  for( const bool fit : setup.fitfor )
    num_fit_coefs += fit;
  for( const bool fit : setup.fit_dev_pair_offsets )
  {
    num_fit_dev += fit;
    fit_any_dev = (fit_any_dev || fit);
  }

  if( (num_fit_coefs + num_fit_dev) < 1 )
    throw runtime_error( "fit_energy_cal_ceres: must fit at least one parameter" );

  if( (num_fit_coefs + num_fit_dev) > peakinfos.size() )
    throw runtime_error( "fit_energy_cal_ceres: must have at least as many peaks as parameters"
                         " being fit" );

  // If fitting deviation pair offsets with only a single input pair, prepend the implicit
  //  (fixed) {0,0} pair, matching the convention SpecUtils uses when applying deviation pairs
  if( fit_any_dev && (setup.dev_pairs.size() == 1) && (setup.dev_pairs[0].first >= 0.1f) )
  {
    setup.dev_pairs.insert( begin(setup.dev_pairs), {0.0f, 0.0f} );
    setup.fit_dev_pair_offsets.insert( begin(setup.fit_dev_pair_offsets), false );
  }

  EnergyCalCeresFitResult result;

  const auto append_warning = [&result]( const std::string &msg ){
    if( !result.warning_msg.empty() )
      result.warning_msg += "\n";
    result.warning_msg += msg;
  };

  // If, even after the implicit {0,0} handling, there arent enough pairs to form a spline (i.e.,
  //  a lone pair below 0.1 keV, which SpecUtils treats as no deviation at all), a fit offset
  //  would enter no residual - Ceres would silently leave it untouched, with a zero Jacobian
  //  column poisoning the covariance - so dont pretend to fit it.
  if( fit_any_dev && (setup.dev_pairs.size() < 2) )
  {
    append_warning( "The deviation pair offset was not fit: a lone pair below 0.1 keV has no"
                    " effect on the energy calibration." );
    std::fill( begin(setup.fit_dev_pair_offsets), end(setup.fit_dev_pair_offsets), false );
    fit_any_dev = false;
    num_fit_dev = 0;

    if( num_fit_coefs < 1 )
      throw runtime_error( "fit_energy_cal_ceres: must fit at least one parameter" );
  }//if( fitting dev pair offsets, but they cant form a spline )

  // Bound (per fitted deviation pair) applied below; 0 means unbounded/not-fit.
  std::vector<double> dev_offset_limit( setup.dev_pairs.size(), 0.0 );

  if( fit_any_dev )
  {
    // A deviation pair offset is only meaningfully constrained if a calibration peak falls in the
    //  energy region where that pair is the closest deviation pair - otherwise fitting it lets
    //  the optimizer overfit through spline interpolation, dragging the coefficients (the offset
    //  especially) to non-physical values.  Hold such pairs fixed, and warn.  For the pairs we do
    //  fit, bound the offset to a fraction of the spacing to the nearest neighbor (capped), so
    //  even a weakly-conditioned pair cant produce a non-physical spike.
    const size_t ndev = setup.dev_pairs.size();
    const double inf = std::numeric_limits<double>::infinity();

    bool warned_unconstrained = false;
    for( size_t i = 0; i < ndev; ++i )
    {
      if( !setup.fit_dev_pair_offsets[i] )
        continue;

      // The energy region where deviation pair i is the closest pair (endpoints reach to infinity)
      const double lo = (i == 0)      ? -inf : 0.5*(setup.dev_pairs[i-1].first + setup.dev_pairs[i].first);
      const double hi = (i+1 == ndev) ?  inf : 0.5*(setup.dev_pairs[i].first + setup.dev_pairs[i+1].first);

      bool has_peak = false;
      for( const EnergyCal::RecalPeakInfo &peak : peakinfos )
        has_peak = (has_peak || ((peak.photopeakEnergy >= lo) && (peak.photopeakEnergy <= hi)));

      if( !has_peak )
      {
        setup.fit_dev_pair_offsets[i] = false;  //hold fixed - nothing constrains it
        if( !warned_unconstrained )
          append_warning( "A deviation pair offset was not fit: no calibration peak falls in the"
                          " energy region it governs, so its value cant be determined." );
        warned_unconstrained = true;
        continue;
      }//if( no peak governs this deviation pair )

      double nearest_spacing = inf;
      if( i > 0 )
        nearest_spacing = std::min( nearest_spacing,
                              static_cast<double>(setup.dev_pairs[i].first - setup.dev_pairs[i-1].first) );
      if( (i+1) < ndev )
        nearest_spacing = std::min( nearest_spacing,
                              static_cast<double>(setup.dev_pairs[i+1].first - setup.dev_pairs[i].first) );

      const double spacing_limit = std::isinf(nearest_spacing) ? 50.0 : (0.15 * nearest_spacing);
      dev_offset_limit[i] = std::max( 0.5, std::min( spacing_limit, 50.0 ) );
    }//for( loop over deviation pairs )

    // Re-tally after possibly holding some fixed
    num_fit_dev = 0;
    fit_any_dev = false;
    for( const bool fit : setup.fit_dev_pair_offsets )
    {
      num_fit_dev += fit;
      fit_any_dev = (fit_any_dev || fit);
    }

    // Fitting the offset coefficient together with ALL the (remaining) deviation pair offsets
    //  leaves a constant-energy-shift degeneracy; the SVD covariance absorbs it, but warn.
    if( fit_any_dev && !setup.fitfor.empty() && setup.fitfor[0] )
    {
      bool all_dev_fit = true;
      for( const bool fit : setup.fit_dev_pair_offsets )
        all_dev_fit = (all_dev_fit && fit);

      if( all_dev_fit )
        append_warning( "Fitting the offset coefficient together with all deviation pair offsets"
                        " is degenerate (a constant energy shift can go into either), so the"
                        " result may be poorly determined." );
    }
  }//if( fit_any_dev )

  // Seed the polynomial/FRF coefficients with a quick linear fit (deviation pairs held at their
  //  current values), so Ceres starts near the solution.  The non-linear deviation-pair fit is
  //  otherwise sensitive to a poor starting calibration and can diverge from it - the linear fit
  //  cant fit the deviation-pair offsets, but it gets the offset/gain right, which is the hard
  //  part for the non-linear solver to recover from a bad start.
  std::vector<float> seed_coefs = setup.starting_coefs;
  if( (num_fit_coefs > 0) && (setup.cal_type != EnergyCalType::LowerChannelEdge) )
  {
    try
    {
      std::vector<float> lls_coefs = setup.starting_coefs, lls_uncert;
      if( setup.cal_type == EnergyCalType::FullRangeFraction )
        fit_energy_cal_frf( peakinfos, setup.fitfor, setup.num_channels, setup.dev_pairs,
                            lls_coefs, lls_uncert );
      else
        fit_energy_cal_poly( peakinfos, setup.fitfor, setup.num_channels, setup.dev_pairs,
                             lls_coefs, lls_uncert );

      bool valid = (lls_coefs.size() == seed_coefs.size());
      for( const float c : lls_coefs )
        valid = (valid && !std::isnan(c) && !std::isinf(c));
      if( valid )
        seed_coefs = lls_coefs;
    }catch( std::exception & )
    {
      //keep setup.starting_coefs as the seed
    }
  }//if( can seed coefficients with a linear fit )

  // Setup and solve the Ceres problem
  EnergyCalCeresCostFcn cost_functor( peakinfos, setup );

  const size_t num_pars = cost_functor.num_parameters();
  assert( num_pars == (ncoefs + setup.dev_pairs.size()) );

  vector<double> pars( num_pars, 0.0 );
  for( size_t i = 0; i < ncoefs; ++i )
    pars[i] = seed_coefs[i] / cost_functor.m_par_scales[i];
  for( size_t i = 0; i < setup.dev_pairs.size(); ++i )
    pars[ncoefs + i] = setup.dev_pairs[i].second;

  vector<int> constant_pars;
  for( size_t i = 0; i < ncoefs; ++i )
  {
    if( !setup.fitfor[i] )
      constant_pars.push_back( static_cast<int>(i) );
  }
  for( size_t i = 0; i < setup.dev_pairs.size(); ++i )
  {
    if( setup.fit_dev_pair_offsets.empty() || !setup.fit_dev_pair_offsets[i] )
      constant_pars.push_back( static_cast<int>(ncoefs + i) );
  }

  auto cost_function = new ceres::DynamicAutoDiffCostFunction<EnergyCalCeresCostFcn,4>(
                              &cost_functor, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP );
  cost_function->AddParameterBlock( static_cast<int>(num_pars) );
  cost_function->SetNumResiduals( static_cast<int>(cost_functor.num_residuals()) );

  ceres::Problem problem;
  problem.AddResidualBlock( cost_function, nullptr, pars.data() );  //problem owns cost_function

  if( !constant_pars.empty() )
  {
    ceres::SubsetManifold * const manifold
                    = new ceres::SubsetManifold( static_cast<int>(num_pars), constant_pars );
    problem.SetManifold( pars.data(), manifold );  //problem owns manifold
  }

  // Bound the fitted deviation-pair offsets, so even a weakly-conditioned pair cant produce a
  //  non-physical spike (which the coefficients would then compensate for, wrecking the fit)
  for( size_t i = 0; i < setup.dev_pairs.size(); ++i )
  {
    if( setup.fit_dev_pair_offsets.empty() || !setup.fit_dev_pair_offsets[i] )
      continue;

    assert( dev_offset_limit[i] > 0.0 );
    const int par_index = static_cast<int>( ncoefs + i );
    const double start = pars[par_index];
    problem.SetParameterLowerBound( pars.data(), par_index, start - dev_offset_limit[i] );
    problem.SetParameterUpperBound( pars.data(), par_index, start + dev_offset_limit[i] );
  }//for( bound the fitted deviation pair offsets )

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_type = ceres::TRUST_REGION;
  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
  options.use_nonmonotonic_steps = true;
  options.max_num_iterations = 250;
  options.function_tolerance = 1.0e-9;
  options.parameter_tolerance = 1.0e-11;
  options.num_threads = 1;
  options.logging_type = ceres::SILENT;
  options.minimizer_progress_to_stdout = false;

  ceres::Solver::Summary summary;
  ceres::Solve( options, &problem, &summary );

  switch( summary.termination_type )
  {
    case ceres::CONVERGENCE:
    case ceres::USER_SUCCESS:
      break;

    case ceres::NO_CONVERGENCE:
      append_warning( "The fit did not fully converge, so the calibration may not be optimal;"
                      " please check the result." );
      break;

    case ceres::FAILURE:
    case ceres::USER_FAILURE:
      throw runtime_error( "fit_energy_cal_ceres: minimization failed: " + summary.message );
  }//switch( summary.termination_type )

  // Extract the results
  result.coefs.resize( ncoefs );
  result.coef_uncerts.assign( ncoefs, 0.0f );
  for( size_t i = 0; i < ncoefs; ++i )
    result.coefs[i] = static_cast<float>( pars[i] * cost_functor.m_par_scales[i] );

  result.dev_pairs = setup.dev_pairs;
  result.dev_pair_offset_uncerts.assign( setup.dev_pairs.size(), 0.0f );
  for( size_t i = 0; i < setup.dev_pairs.size(); ++i )
    result.dev_pairs[i].second = static_cast<float>( pars[ncoefs + i] );

  // Parameter uncertainties, from the SVD pseudo-inverse based covariance
  {
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::DENSE_SVD;
    cov_options.null_space_rank = -1;
    cov_options.min_reciprocal_condition_number = 1.0e-14;

    ceres::Covariance covariance( cov_options );
    vector<pair<const double *, const double *>> cov_blocks;
    cov_blocks.emplace_back( pars.data(), pars.data() );

    vector<double> cov( num_pars * num_pars, 0.0 );
    if( covariance.Compute( cov_blocks, &problem )
        && covariance.GetCovarianceBlock( pars.data(), pars.data(), cov.data() ) )
    {
      for( size_t i = 0; i < ncoefs; ++i )
      {
        const double var = cov[i*num_pars + i];
        if( setup.fitfor[i] && (var > 0.0) )
          result.coef_uncerts[i] = static_cast<float>( cost_functor.m_par_scales[i] * std::sqrt(var) );
      }

      for( size_t i = 0; i < setup.dev_pairs.size(); ++i )
      {
        const size_t index = ncoefs + i;
        const double var = cov[index*num_pars + index];
        if( !setup.fit_dev_pair_offsets.empty() && setup.fit_dev_pair_offsets[i] && (var > 0.0) )
          result.dev_pair_offset_uncerts[i] = static_cast<float>( std::sqrt(var) );
      }
    }else
    {
      append_warning( "Parameter uncertainties could not be computed." );
    }
  }//end code-block to compute uncertainties

  // Compute chi2 using the SpecUtils forward functions, consistent with the linear fits
  result.chi2 = 0.0;
  for( size_t i = 0; i < peakinfos.size(); ++i )
  {
    const EnergyCal::RecalPeakInfo &peak = peakinfos[i];

    double pred = 0.0;
    switch( setup.cal_type )
    {
      case EnergyCalType::Polynomial:
      case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        pred = SpecUtils::polynomial_energy( peak.peakMeanBinNumber, result.coefs, result.dev_pairs );
        break;

      case EnergyCalType::FullRangeFraction:
        pred = SpecUtils::fullrangefraction_energy( peak.peakMeanBinNumber, result.coefs,
                                                    setup.num_channels, result.dev_pairs );
        break;

      case EnergyCalType::LowerChannelEdge:
        pred = result.coefs[0] + result.coefs[1]*cost_functor.m_lower_chan_orig[i];
        break;

      case EnergyCalType::InvalidEquationType:
        assert( 0 );
        break;
    }//switch( setup.cal_type )

    result.chi2 += std::pow( (pred - peak.photopeakEnergy) / cost_functor.m_sigmas[i], 2.0 );
  }//for( loop over peaks )

  // Hard validity check: make sure the result actually forms a valid calibration
  try
  {
    auto check_cal = make_shared<SpecUtils::EnergyCalibration>();
    switch( setup.cal_type )
    {
      case EnergyCalType::Polynomial:
      case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        check_cal->set_polynomial( setup.num_channels, result.coefs, result.dev_pairs );
        break;

      case EnergyCalType::FullRangeFraction:
        check_cal->set_full_range_fraction( setup.num_channels, result.coefs, result.dev_pairs );
        break;

      case EnergyCalType::LowerChannelEdge:
      {
        auto orig_cal = make_shared<SpecUtils::EnergyCalibration>();
        orig_cal->set_lower_channel_energy( setup.num_channels, setup.lower_channel_energies );
        adjust_lower_channel_energy_cal( orig_cal, result.coefs[0], result.coefs[1] );
        break;
      }

      case EnergyCalType::InvalidEquationType:
        assert( 0 );
        break;
    }//switch( setup.cal_type )
  }catch( std::exception &e )
  {
    throw runtime_error( "fit_energy_cal_ceres: fit produced an invalid calibration: "
                         + std::string(e.what()) );
  }//try / catch

  return result;
}//fit_energy_cal_ceres(...)



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

  
  assert( (orig_cal->type() != EnergyCalType::InvalidEquationType) );
  assert( (new_cal->type() != EnergyCalType::InvalidEquationType) );
  assert( (other_cal->type() == EnergyCalType::FullRangeFraction)
          || (other_cal->type() == EnergyCalType::Polynomial)
          || (other_cal->type() == EnergyCalType::UnspecifiedUsingDefaultPolynomial) );

  const double accuracy = 0.00001;

  // For a lower-channel-energy calibration coefficients() returns all the channel energies, so
  //  dont let it drive the number of sampled points (its offset/gain style adjustments are
  //  linear, so two points would do; other_cal's order still gets honored).
  const size_t prev_order = (orig_cal->type() == EnergyCalType::LowerChannelEdge)
                            ? size_t(2) : prev_disp_coefs.size();
  const size_t new_order = (new_cal->type() == EnergyCalType::LowerChannelEdge)
                            ? size_t(2) : new_disp_coefs.size();
  const size_t order = std::max( other_coeffs.size(), std::max(prev_order, new_order) );
  
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
      {
        // Lower channel energy calibrations dont have deviation pairs, so using the member
        //  functions here is consistent with the dev-pair-ignoring convention above
        old_disp_energy = orig_cal->energy_for_channel( static_cast<double>(display_channel) );
        new_disp_energy = new_cal->energy_for_channel( static_cast<double>(display_channel) );
        break;
      }//case: orig_cal was LowerChannelEdge

      case EnergyCalType::InvalidEquationType:
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

  const Eigen::VectorXd a = solve_lls_equilibrated( A, b, nullptr );

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

  const Eigen::VectorXd a = solve_lls_equilibrated( A, b, nullptr );

  vector<float> frf_coeffs( polyterms );
  for( int coef = 0; coef < polyterms; ++coef )
    frf_coeffs[coef] = static_cast<float>( a(coef) );
  
  return frf_coeffs;
}//void fit_for_fullrangefraction_coefs(...)


