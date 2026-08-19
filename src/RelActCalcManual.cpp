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

#include <set>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <chrono>
#include <vector>
#include <memory>
#include <optional>
#include <iostream>

#include <boost/math/distributions/chi_squared.hpp>

#include <Wt/WApplication>

#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#undef ERROR
#endif

#include "Eigen/Dense"
#include "ceres/ceres.h"
#include "ceres/loss_function.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "InterSpec/RelActCalc_imp.hpp"
#include "InterSpec/RelActCalc_CeresJetTraits.hpp"

using namespace std;


namespace
{
struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct


/** Console-chatter gate for the solver: set the INTERSPEC_RELACT_DEBUG environment variable to a
 non-"0" value to re-enable Ceres iteration progress and the solver diagnostic printouts.  An
 env-var (rather than a compile flag) so it can be flipped on a production GUI server without a
 rebuild - and so debug builds stay quiet for tests.
 */
bool debug_printout()
{
  static const bool s_debug = [](){
    const char * const env_val = std::getenv( "INTERSPEC_RELACT_DEBUG" );
    return env_val && env_val[0] && strcmp(env_val, "0");
  }();

  return s_debug;
}//bool debug_printout()


/** Whether this Physical Model shielding contributes any Ceres parameters: it has a material,
 or its atomic number is being fit, or it has a valid fixed atomic number in [1,98].

 These are the canonical thresholds - previously several call sites re-implemented this
 predicate with slightly different cutoffs (0.99 / 0.999 / 1.0), which only agreed because
 `PhysicalModelShieldInput::check_valid()` forbids AN in (0,1).  Any drift between call sites
 is a parameter-index-shift bug, so use these helpers everywhere.
 */
bool shield_is_present( const RelActCalc::PhysicalModelShieldInput * const opt )
{
  if( !opt )
    return false;

  return opt->material || opt->fit_atomic_number
         || ((opt->atomic_number >= 1.0) && (opt->atomic_number <= 98.0));
}//bool shield_is_present(...)

bool shield_is_present( const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt )
{
  return shield_is_present( opt.get() );
}

/** Whether the shields atomic number is a fitted Ceres parameter (no material, and AN fit). */
bool shield_fits_an( const RelActCalc::PhysicalModelShieldInput * const opt )
{
  return opt && !opt->material && opt->fit_atomic_number;
}

bool shield_fits_an( const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt )
{
  return shield_fits_an( opt.get() );
}

/** We can either use a residual to force the normalization of the R.E. curve to 1.0 at the lowest
 * energy.  Or we can manually force the average measured relative efficiency to 1.0.
 * This does not apply to `RelEffEqnForm::FramPhysicalModel`.
 */
#define USE_RESIDUAL_TO_BREAK_DEGENERACY 0

  
/** Computes the linear-least-squares parameter covariance `C = (A^T A)^{-1}` from the (thin) SVD of
 `A = U*S*V^T`, using the numerically-stable pseudo-inverse `C = V * S^{-2} * V^T`.  Singular values
 below `tol*sigma_max` are dropped (treated as 1/s^2 -> 0), so this avoids forming `A^T A` (which
 squares the condition number) and yields a finite, minimum-norm covariance for rank-deficient `A`
 instead of a blow-up or a slightly-negative diagonal.

 The covariance is only ever needed in double precision: the auto-differentiation (`ceres::Jet`)
 residual path passes a null covariance pointer (it is unused there, and a Jet matrix pseudo-inverse
 would be expensive and is not needed).  The `static_assert` makes that a compile-time requirement for
 any future caller. */
template<typename T>
Eigen::MatrixX<T> lls_covariance_from_svd( const Eigen::MatrixX<T> &V,
                                           const Eigen::VectorX<T> &singular_values,
                                           const Eigen::Index num_rows )
{
  static_assert( std::is_same_v<T, double>,
    "lls_covariance_from_svd: covariance is computed in double precision only.  The auto-diff"
    " (ceres::Jet) path passes a null covariance pointer (it is unused there); the SVD pseudo-inverse"
    " is neither needed nor implemented for Jet types." );

  const Eigen::Index p = V.rows();          // parameter dimension (covariance is p x p)
  const Eigen::Index num_sv = singular_values.size();

  // Eigen returns singular values in descending order, so element 0 is the largest.
  const double sigma_max = (num_sv > 0) ? singular_values(0) : 0.0;
  const double tol = std::numeric_limits<double>::epsilon()
                       * static_cast<double>( std::max<Eigen::Index>( num_rows, p ) )
                       * sigma_max;

  Eigen::MatrixX<T> C = Eigen::MatrixX<T>::Zero( p, p );
  for( Eigen::Index k = 0; k < num_sv; ++k )
  {
    const double s = singular_values(k);
    if( s > tol )
      C.noalias() += (1.0 / (s*s)) * (V.col(k) * V.col(k).transpose());
  }//for( loop over singular values )

  return C;
}//lls_covariance_from_svd(...)


/** Smooth, gradient-preserving lower bound: ~`value` for `value >> floor_val`, ~`floor_val` for
 `value << floor_val`, with slope always in (0,1).  Used to keep a quantity strictly positive before
 a log()/division without zeroing a `ceres::Jet` derivative (a hard clamp `value = floor_val` would
 zero it, re-creating the bound-trap the physical-model AD code avoids). */
template<typename T>
T smooth_lower_bound( const T &value, const double floor_val )
{
  using namespace std;
  using namespace ceres;
  const T d = value - floor_val;
  return floor_val + 0.5*( d + sqrt( d*d + floor_val*floor_val ) );
}//smooth_lower_bound(...)


template<typename T>
void fit_rel_eff_eqn_lls_imp( const RelActCalc::RelEffEqnForm fcn_form,
                              const size_t order,
                              const vector<double> &energies,
                              const vector<T> &data_values,
                              const vector<T> &data_uncertainties_orig,
                              vector<T> &fit_pars,
                              vector<vector<T>> *covariance )
{
  using namespace std;
  using namespace ceres;
  
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  const vector<T> &data_uncertainties = data_uncertainties_orig;
  
  assert( !data_values.empty() );
  assert( energies.size() == data_values.size() );
  assert( energies.size() == data_uncertainties.size() );
  
  if( data_values.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no data points provided." );
  
  const int poly_terms = static_cast<int>(order) + 1;
  const int num_peaks = static_cast<int>( data_values.size() );
  
  
  Eigen::MatrixX<T> A( num_peaks, poly_terms );
  Eigen::VectorX<T> b( num_peaks );
  
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const double energy = energies[row];
    const T measured_rel_eff = data_values[row];
    T uncertainty( data_uncertainties[row] );
    
    // Letting negative relative efficiencies through doesnt feel right, but I guess we'll do
    //  it to not mess up the L-M fitting of relative activities....
    //if( measured_rel_eff <= 0.0 )
    //  throw runtime_error( "fit_rel_eff_eqn_lls: Measured relative efficiency for energy "
    //                      + to_string(energy) + " is invalid ("
    //                      + to_string(measured_rel_eff) + ")" );
    
    //  But we'll put our foot down for negative or zero uncertainties.
    if( uncertainty <= 0.0 )
      throw runtime_error( "fit_rel_eff_eqn_lls: Uncertainty for energy " + to_string(energy)
                          + " is invalid." );
    
    switch( fcn_form )
    {
      case RelActCalc::RelEffEqnForm::LnX:
      {
        //y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
        b(row) = measured_rel_eff / uncertainty;
        break;
      }
        
      case RelActCalc::RelEffEqnForm::LnY:
      case RelActCalc::RelEffEqnForm::LnXLnY:
      case RelActCalc::RelEffEqnForm::FramEmpirical:
      {
        //LnY:           y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
        //LnXLnY:        y = exp (a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
        //FramEmpirical: y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
        //
        // We'll take the log of each side of the equation, and then solve for the parameters
        //
        // Note: when `f = ln(x)`, then `uncert(f) = uncert(x)/x`
        
        uncertainty = uncertainty / measured_rel_eff;
        
        // Note that we get the same answer (for a few problems I checked) if we use the following
        //  approximation to estimate uncertainty.
        //if( measured_rel_eff < 2.0*uncertainty )
        //  uncertainty = (2.0*uncertainty/measured_rel_eff) * fabs( log(0.75*measured_rel_eff) - log(1.25*measured_rel_eff) );
        //else
        //  uncertainty = 2.0*fabs( std::log(measured_rel_eff - 0.25*uncertainty)
        //                             - std::log(measured_rel_eff + 0.25*uncertainty) );
        
        b(row) = log(measured_rel_eff) / uncertainty;
        
        break;
      }
        
      case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      {
        assert( 0 );
        throw runtime_error( "fit_rel_eff_eqn_lls: FramPhysicalModel not supported." );
      }
    }//switch( fcn_form )
    
    
    for( int col = 0; col < poly_terms; ++col )
    {
      switch( fcn_form )
      {
        case RelActCalc::RelEffEqnForm::LnX:
        case RelActCalc::RelEffEqnForm::LnXLnY:
          //y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
          // and
          //ln(y) = a + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ...
          A(row,col) = pow(log(energy), static_cast<double>(col)) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::LnY:
          //ln(y) = a + b*x + c/x + d/x^2 + e/x^3 + ...
          if( col == 0 )
            A(row,col) = 1.0 / uncertainty;
          else if( col == 1 )
            A(row,col) = energy / uncertainty;
          else
            A(row,col) = pow(energy, 1.0 - col) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::FramEmpirical:
          //ln(y) = a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3
          if( col == 0 )
            A(row,col) = 1.0 / uncertainty;
          else if( col == 1 )
            A(row,col) = (1.0 / (energy*energy)) / uncertainty;
          else
            A(row,col) = pow(log(energy), col - 1.0) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::FramPhysicalModel:
          assert( 0 );
          throw runtime_error( "fit_rel_eff_eqn_lls: FramPhysicalModel not supported." );
      }//switch( fcn_form )
    }//for( int col = 0; col < poly_terms; ++col )
  }//for( int col = 0; col < poly_terms; ++col )
  
  // TODO: determine if HouseholderQr or BDC SVD is better/more-stable/faster/whatever
  //const Eigen::VectorX<T> solution = A.colPivHouseholderQr().solve(b);

  // deprecated way to compute the BDCSVD matrix
  //const Eigen::BDCSVD<Eigen::MatrixX<T>> bdc = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  // What I think is the updated way.
  Eigen::BDCSVD<Eigen::MatrixX<T>,Eigen::ComputeThinU | Eigen::ComputeThinV> bdc;
  bdc.compute(A);
  const Eigen::VectorX<T> solution = bdc.solve(b);
  
  assert( solution.size() == (order + 1) );
  
  fit_pars.resize( solution.size() );
  for( size_t i = 0; i <= order; ++i )
    fit_pars[i] = solution(i);
  
  // Only compute covariance if it is wanted
  if( covariance )
  {
    if constexpr ( std::is_same_v<T, double> )
    {
      // Covariance C = (A^T A)^{-1}, computed as the SVD pseudo-inverse V*S^{-2}*V^T (see helper);
      //  this avoids forming A^T A (which squares the condition number) and stays finite for
      //  rank-deficient A.
      const Eigen::MatrixX<T> C = lls_covariance_from_svd<T>( bdc.matrixV(), bdc.singularValues(),
                                                            static_cast<Eigen::Index>(num_peaks) );

      assert( C.rows() == solution.size() );
      assert( C.cols() == solution.size() );

      covariance->resize( solution.size() );

      for( size_t i = 0; i <= order; ++i )
      {
        vector<T> &row = (*covariance)[i];
        row.resize( solution.size() );
        for( size_t j = 0; j <= order; ++j )
          row[j] = C(i,j);
      }//for( loop over coefficients index )
    }else
    {
      // The auto-diff (ceres::Jet) residual path passes covariance==nullptr (it is unused there),
      //  so this is unreachable; we never want to differentiate through a matrix pseudo-inverse.
      assert( 0 );
      throw std::logic_error( "fit_rel_eff_eqn_lls_imp: covariance output is not supported for"
                              " ceres::Jet types - pass a null covariance pointer." );
    }//if constexpr( T is double ) / else
  }//if( covariance )
}//fit_rel_eff_eqn_lls_imp(...)


template<typename T>
void fit_rel_eff_eqn_lls_imp( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<std::string> &isotopes,
                           const std::vector<T> &rel_acts,
                           const std::vector<RelActCalcManual::GenericPeakInfo> &peak_infos,
                           std::vector<T> &fit_pars,
                           std::vector<std::vector<T>> *covariance )
{
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  assert( !isotopes.empty() );
  if( isotopes.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no isotopes specified." );
  
  const int poly_terms = static_cast<int>(order) + 1;
  const int num_peaks = static_cast<int>( peak_infos.size() );
  
  vector<double> energies( num_peaks, 0.0 );
  vector<T> meas_rel_eff( num_peaks, T(0.0) ), meas_rel_eff_uncert( num_peaks, T(0.0) );
  
  bool unweighted_fit = false;
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const RelActCalcManual::GenericPeakInfo &peak = peak_infos[row];
    
    // A basic sanity check that the uncertainty in counts isnt garbage.
    if( peak.m_counts_uncert < 0.0 )
      throw runtime_error( "fit_rel_eff_eqn_lls: peak counts uncertainty can not be <0" );
    
    // Check there is a non-zero peak counts uncertainty; if its zero, we'll (arbitrarily) restrict
    //  to doing an un-weighted fit.  We could accept any non-zero peak.m_base_rel_eff_uncert
    //  and compute things just fine, but this would be highly suspect that the user has messed
    //  up filling out peak information, so we'll throw an exception.
    if( (peak.m_counts_uncert == 0.0) && (peak.m_base_rel_eff_uncert != -1.0) )
      throw runtime_error( "fit_rel_eff_eqn_lls: you must either provide a non-zero peak counts"
                          " uncertainty, or perform a unweighted fit" );
    
    const double energy = peak.m_energy;
    const double counts = peak.m_counts;
    const double counts_uncert = peak.m_counts_uncert;
    
    T raw_rel_counts( 0.0 );
    
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
    {
      const auto iso_pos = std::lower_bound( std::begin(isotopes), std::end(isotopes), line.m_isotope );

      // `lower_bound` requires `isotopes` sorted (the public overload sorts caller input), and
      //  without the equality check a missing isotope would silently take the next entry's
      //  activity.  Reachable from the public API with bad input, so throw rather than assert.
      if( (iso_pos == std::end(isotopes)) || ((*iso_pos) != line.m_isotope) )
        throw std::logic_error( "fit_rel_eff_eqn_lls: peak source '" + line.m_isotope
                                + "' is not in the isotope list" );

      const size_t iso_index = static_cast<size_t>( iso_pos - std::begin(isotopes) );
      const T rel_act_value = rel_acts[iso_index];
      
      raw_rel_counts += line.m_yield * rel_act_value;
    }//for( const GenericLineInfo &line : peak.m_source_gammas )

    // The optimizer is free to probe zero relative activity (activities are bounded below by 0),
    //  which drives `raw_rel_counts` to zero for any peak fed only by that source.  This used to
    //  throw, which the Ceres functor turned into a failed evaluation and hence a trust-region
    //  shrink - an expensive way to say "that point is bad", and it happens on ordinary LnX fits.
    //  Instead floor it smoothly: the point stays evaluable with a finite, correctly-signed
    //  gradient, so the optimizer walks away from it on its own.
    //
    //  The floor is scale-free (tied to the peak's own counts) and matches the clamp the Ceres
    //  twin already uses in `eval_internal_nl_rel_eff`.  Note the weighted least-squares row is
    //  well behaved in the limit: `measured_rel_eff` and `measured_rel_eff_uncert` both scale as
    //  1/raw_rel_counts, so their ratio - the pull this peak exerts - stays at counts/counts_uncert
    //  while the design-matrix row shrinks toward zero, i.e. the peak simply stops constraining
    //  the curve rather than blowing up the normal equations.
    const double raw_rel_counts_floor = 1.0E-6 * (std::max)( 1.0, counts );
    raw_rel_counts = smooth_lower_bound( raw_rel_counts, raw_rel_counts_floor );

    T measured_rel_eff = counts / raw_rel_counts;
    T measured_rel_eff_uncert = counts_uncert / raw_rel_counts;

    // Keep rel-eff in a usable range (workaround since Eigen's LM doesn't constrain parameter ranges).
    //  The log-based forms take log(measured_rel_eff) and divide by it downstream, so a non-positive
    //  or non-finite value poisons the fit; use a smooth, gradient-preserving floor for those (a hard
    //  clamp to a constant would zero the ceres::Jet derivative w.r.t. the activities - the same
    //  bound-trap the physical-model AD code deliberately avoids).  LnX uses rel-eff linearly, so a
    //  hard clamp to zero remains valid there.
    const double rel_eff_floor = static_cast<double>(numeric_limits<float>::epsilon());
    const bool log_form = (fcn_form == RelActCalc::RelEffEqnForm::LnY)
                          || (fcn_form == RelActCalc::RelEffEqnForm::LnXLnY)
                          || (fcn_form == RelActCalc::RelEffEqnForm::FramEmpirical);

    if( log_form )
    {
      if( isnan(measured_rel_eff) || isinf(measured_rel_eff) )
      {
        // raw_rel_counts ~ 0: no usable measurement (no meaningful derivative either) - floor and down-weight.
        measured_rel_eff = T( rel_eff_floor );
        measured_rel_eff_uncert = T( 1.0 );
      }else
      {
        measured_rel_eff = smooth_lower_bound( measured_rel_eff, rel_eff_floor );
      }
    }else if( (measured_rel_eff <= rel_eff_floor) || isinf(measured_rel_eff) || isnan(measured_rel_eff) )
    {
      measured_rel_eff = T( 0.0 );
      if( peak.m_base_rel_eff_uncert > rel_eff_floor )
        measured_rel_eff_uncert = T( 0.0 );
      else
        measured_rel_eff_uncert = T( 1.0 );
    }
    
    if( peak.m_base_rel_eff_uncert == -1.0 )
    {
      if( row && !unweighted_fit )
        throw runtime_error( "fit_rel_eff_eqn_lls: for unweighted fit, all peaks must specify m_base_rel_eff_uncert == -1" );
      unweighted_fit = true;
      measured_rel_eff_uncert = T( 1.0 );
    }else
    {
      if( unweighted_fit )
        throw runtime_error( "fit_rel_eff_eqn_lls: for unweighted fit, all peaks must specify m_base_rel_eff_uncert == -1" );
      
      if( (peak.m_base_rel_eff_uncert < 0.0) || (peak.m_base_rel_eff_uncert > 1.0) )
        throw runtime_error( "fit_rel_eff_eqn_lls: m_base_rel_eff_uncert must be in range [0,1]" );
      
      if( peak.m_base_rel_eff_uncert > 0.0 )
      {
        // We should to be consistent with #ManualGenericRelActFunctor::eval in how we compute the
        // uncertainty
        
        const T add_uncert( counts * peak.m_base_rel_eff_uncert );
        measured_rel_eff_uncert = sqrt( pow(counts_uncert,2.0) + pow(add_uncert, 2.0) );
        measured_rel_eff_uncert /= raw_rel_counts;
        assert( !isinf(counts) && !isnan(counts) );
        assert( !isinf(measured_rel_eff_uncert) && !isnan(measured_rel_eff_uncert) );
        assert( !isinf(raw_rel_counts) && !isnan(raw_rel_counts) );
        assert( !isinf(raw_rel_counts) && !isnan(raw_rel_counts) );
        assert( !isinf(measured_rel_eff_uncert) && !isnan(measured_rel_eff_uncert) );
      }//if( peak.m_base_rel_eff_uncert > 0.0 )
      
      // else keep as counts_uncert / raw_rel_counts
    }//if( do unweighted fit ) / else
    
    
    energies[row] = energy;
    meas_rel_eff[row] = measured_rel_eff;
    meas_rel_eff_uncert[row] = measured_rel_eff_uncert;
    //cout << "energy(" << energy << "): meas_rel_eff[" << row << "] = " << meas_rel_eff[row] << ", meas_rel_eff_uncert[" << row << "] = " << meas_rel_eff_uncert[row] << endl;
    assert( !isinf(measured_rel_eff_uncert) && !isnan(measured_rel_eff_uncert) );
  }//for( int col = 0; col < poly_terms; ++col )
  
  
#if( !USE_RESIDUAL_TO_BREAK_DEGENERACY )
  // Pin the measured rel-eff points to an average of 1.0; this is what breaks the overall
  //  activities-vs-curve scale degeneracy in the LLS fit mode (the counts-space residuals then
  //  penalize any deviation of the average from 1).  The uncertainties are divided by the same
  //  factor so absolute (non-log) uncertainty use stays consistent with the normalized values.
  //  See also ManualGenericRelActFunctor::average_measured_rel_eff(), which must compute the
  //  identical normalization.
  const T sum_re = std::accumulate( begin(meas_rel_eff), end(meas_rel_eff), T(0.0) ); //Previous to 20250110, a value of 1,0 was used to initialize accumulate - not sure want that was, should probably check on this again
  const T average_re = sum_re / static_cast<double>( meas_rel_eff.size() );
  for( T &re : meas_rel_eff )
    re /= average_re;
  for( T &re_uncert : meas_rel_eff_uncert )
    re_uncert /= average_re;
#endif
  
  
  fit_rel_eff_eqn_lls_imp<T>( fcn_form, order,
                          energies, meas_rel_eff, meas_rel_eff_uncert,
                          fit_pars, covariance );
}//fit_rel_eff_eqn_lls_imp(...)

/*
template<typename T>
void fit_rel_eff_eqn_lls_imp( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<SandiaDecayNucRelAct<T>> &nuclides,
                           const double base_rel_eff_uncert,
                           const std::vector<std::shared_ptr<const PeakDef>> &peak_infos,
                           vector<T> &fit_pars,
                           std::vector<std::vector<T>> *covariance )
{
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  assert( !nuclides.empty() );
  if( nuclides.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no nuclides specified." );
  
  assert( !peak_infos.empty() );
  if( peak_infos.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no peaks specified." );
  
  // We will map from the peaks mean, to the total number of gammas that contribute to that peak
  map<double,double> energy_gammas_map;
  
  // Get peak energies and widths (normally width this is just 'sigma', but for non-Gaussian peaks
  //  its 0.25 of the ROI)
  set<double> energies_seen; //a very poor check that there arent duplicate peaks
  vector< pair<double,double> > energy_widths, energy_obs_counts, energy_obs_counts_uncert;
  for( const auto &p : peak_infos )
  {
    // To use non-Gaussian peaks we would need to pass in the shared_ptr<const Measurement> data...
    //  maybe later if it ever matters
    if( !p->gausPeak() )
      throw runtime_error( "fit_rel_eff_eqn_lls: non-Gaussian peaks not supported yet" );
    
    // Note that in GammaInteractionCalc::ShieldingSourceChi2Fcn::observedPeakEnergyWidths
    //  we use the assigned nuclides gamma energy, as the energy - here we are using the peak mean.
    //  TODO: - revisit ether to use peak mean or its nuclide gamma as the energy - after implementing the rest of the manual RelAct calc stuff.
    const double energy = p->mean();
    const double sigma = p->gausPeak() ? p->sigma() : 0.25*p->roiWidth();
    const double amp = p->amplitude();
    //const double amp = p->gausPeak() ? p->amplitude() : p->areaFromData(data);
    const double ampUncert = p->amplitudeUncert();
    
    if( energies_seen.count(energy) )
      throw runtime_error( "fit_rel_eff_eqn_lls: multiple peaks with same energy - not allowed." );
    energies_seen.insert( energy );
    
    energy_widths.push_back( {energy, sigma} );
    energy_obs_counts.push_back( {energy, amp} );
    energy_obs_counts_uncert.push_back( {energy, ampUncert} );
  }//for( const PeakDef &peak : peaks )
  
  // JIC the peaks werent sorted, sort by just energies (although we did check no duplicate energies
  //  but we'll play it safe)
  auto sortByFirstOnly = []( const pair<double, double> &lhs, const pair<double, double> &rhs ){
    return lhs.first < rhs.first;
  };
  
  std::stable_sort( begin(energy_widths), end(energy_widths), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts), end(energy_obs_counts), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts_uncert), end(energy_obs_counts_uncert), sortByFirstOnly );
  
  
  // Now we will go through and get the amplitude of gammas we expect to contribute to a single peak
  //  (there may be multiple gammas from the same nuclide, as well as multiple nuclides that
  //   contribute to a single observable peaks).
  //  We will select a 'cluster' sigma of 1.5; this is what Activity/Shielding fit uses, but I dont
  //  think this value was derived by anything more than "that seems about right", and I havent
  //  run into an obvious case where this is not correct.
  const double photopeakClusterSigma = 1.5;
  set<const SandiaDecay::Nuclide *> nuclides_seen;
  for( const auto &n : nuclides )
  {
    if( nuclides_seen.count(n.nuclide) )
      throw runtime_error( "fit_rel_eff_eqn_lls: input nuclides must be unique" );
    
    nuclides_seen.insert( n.nuclide );
    
    SandiaDecay::NuclideMixture mixture;
    mixture.addNuclideByActivity(n.nuclide, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits);
    
    const double energyToCluster = -1;
    // TODO: we could account for decays during the measurement, but would need realTime here
    const bool accountForDecayDuringMeas = false;
    const double realTime = -1;
    GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map,
                                                                          energy_widths, mixture, n.rel_activity, n.age,
                                                                          photopeakClusterSigma, energyToCluster,
                                                                          accountForDecayDuringMeas, realTime, nullptr, nullptr );
  }//for( const auto &n : nuclides )
  
  // Convert energy_gammas_map to a vector for convenience
  vector<pair<double,double>> energy_gammas;
  for( const auto &ec : energy_gammas_map )
    energy_gammas.push_back( ec );
  
  assert( energy_gammas.size() == peak_infos.size() );
  assert( energy_gammas.size() == energy_widths.size() );
  assert( energy_gammas.size() == energy_obs_counts.size() );
  assert( energy_gammas.size() == energy_obs_counts_uncert.size() );
  
  
  // Now put all this info onto a form so we can call into fit_rel_eff_eqn_lls(...); there is a
  //  commented out implementation of not having to do this, yet another, transformation of
  //  information.
  
  double max_pred_counts = 0.0;
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
    max_pred_counts = std::max(max_pred_counts, energy_gammas[peak_index].second);
  
  // Instead of keeping counts from each nuclide for each peak separate, we summed all nuclides
  //  together for each peak - so here we'll only use a single "Effective" isotope.
  const vector<string> isotopes{ "Effective" };
  const vector<double> rel_acts( 1, max_pred_counts );
  vector<GenericPeakInfo> generic_peak_infos;
  
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
  {
    GenericPeakInfo peak;
    peak.m_mean = peak_infos[peak_index]->mean();
    peak.m_energy = energy_gammas[peak_index].first;
    peak.m_fwhm = 2.35482*energy_widths[peak_index].second;
    peak.m_counts = energy_obs_counts[peak_index].second;
    peak.m_counts_uncert = (energy_obs_counts_uncert[peak_index].second > 0.0)
    ? energy_obs_counts_uncert[peak_index].second
    : sqrt(peak.m_counts);
    peak.m_base_rel_eff_uncert = base_rel_eff_uncert;
    
    const double yield = energy_gammas[peak_index].second / max_pred_counts;
    peak.m_source_gammas.emplace_back( yield, "Generic" );
    
    generic_peak_infos.push_back( peak );
  }//for( size_t row = 0; row < num_peaks; ++row )
  
  return fit_rel_eff_eqn_lls_imp<T>( fcn_form, order, isotopes, rel_acts, generic_peak_infos,
                             fit_pars, covariance );
}//fit_rel_eff_eqn_lls_imp(...)
*/
  

/** Setups the parameters for a shielding.
  
 Note that the number of parameters used is variable; when AN is being fit, a parameter will be used, otherwise not.
 A parameter for AD is always used.
 This was a result of the development path (specifically trouble getting auto-differentiation to work, and holding and
 manifolds to behave), and I think could be removed so number of parameter is fixed, but maybe we'll wait until
 we make another change (like being able to fix ratios of nuclides, or something) to change this
*/
void setup_physical_model_shield_par_manual( vector<int> &constant_parameters,
                                            double * const pars,
                                            size_t &index,
                                            const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt )
{
  // Note: deliberately narrower than `shield_is_present()` - an *invalid* fixed AN (e.g. 0.5)
  //  falls through to check_valid() below for a clean error, instead of being silently skipped.
  //  For inputs check_valid() accepts, this condition is exactly `!shield_is_present(opt)`.
  if( !opt || (!opt->material && (opt->atomic_number == 0.0) && !opt->fit_atomic_number) )
    return;

  opt->check_valid();
  
  if( opt->material )
  {
    if( opt->fit_atomic_number )
      throw runtime_error( "You can not fit AN when defining a material" );
  }else if( opt->fit_atomic_number )
  {
    double lower_an = opt->lower_fit_atomic_number;
    double upper_an = opt->upper_fit_atomic_number;
    if( (lower_an == upper_an) && (lower_an == 0.0) )
    {
      lower_an = 1.0;
      upper_an = 98.0;
    }
    
    double an = opt->atomic_number;
    if( an < 1.0 )
    {
      assert( an == 0.0 );
      an = 0.5*(opt->lower_fit_atomic_number + opt->upper_fit_atomic_number);
    }
    
    if( (an < 1) || (an > 98) || (an < lower_an) || (an > upper_an) )
      throw runtime_error( "Self-atten AN is invalid" );
    
    pars[index] = an / RelActCalc::ns_an_ceres_mult;
    index += 1; //Add parameter for AN
    
    if( (lower_an < 1) || (upper_an > 98) || (upper_an <= lower_an) )
      throw runtime_error( "Self-atten AN limits is invalid" );
  }else
  {
    double an = opt->atomic_number;
    if( (an < 1) || (an > 98) )
      throw runtime_error( "Self-atten fixed AN is invalid" );
    
    // We wont actually add/use a parameter for fixed AN
  }//if( opt.material ) / else if( opt->fit_atomic_number )
  
  const double max_ad = RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2;
  double ad = opt->areal_density / PhysicalUnits::g_per_cm2;
  double lower_ad = opt->lower_fit_areal_density / PhysicalUnits::g_per_cm2;
  double upper_ad = opt->upper_fit_areal_density / PhysicalUnits::g_per_cm2;
  
  if( (lower_ad == upper_ad) && (lower_ad == 0.0) )
  {
    lower_ad = 0.0;
    upper_ad = max_ad;
  }
  
  if( (ad == 0.0) && opt->fit_areal_density )
  {
    //ad = 0.5*(lower_ad + upper_ad); //Something like 250 would be way too much
    ad = 2.5; // We want something away from zero, because Ceres doesnt like zero values much - 2.5 is arbitrary
  }
  
  if( (ad < 0.0) || (ad > max_ad) )
    throw runtime_error( "Self-atten AD is invalid" );
  
  pars[index] = ad;
  
  
  if( opt->fit_areal_density )
  {
    // Check for limits
    if( (lower_ad < 0.0) || (upper_ad > max_ad) || (lower_ad >= upper_ad) )
      throw runtime_error( "Self-atten AD limits is invalid" );
  }else
  {
    constant_parameters.push_back( static_cast<int>(index) );
  }
  
  index += 1; //Add parameter for AD, always
}//void setup_physical_model_shield_par_manual( ceres::Problem... )


  
/** Functor for minimizing the relative activities; relative efficiency is fit for each set of
 activities via the #fit_rel_eff_eqn_lls function.
 */
struct ManualGenericRelActFunctor  /* : ROOT::Minuit2::FCNBase() */
{
  /** The form of relative efficiency equation to use. */
  //const RelActCalc::RelEffEqnForm m_eqn_form;
  
  /** The order of relative efficiency equation to use (not equation will have one more than this
   value coefficients)
   */
  //const int m_eqn_order;
  
  /** All isotope relative activities will be fit for using LM, with the relative eff curve forced
   to have a value near 1.0 for the first peak.
   */
  std::vector<string> m_isotopes;
  
  /** We will first normalize relative efficiency for each isotope to a flat line at y = 1.0.
   We will then use L-M to fit the multiples of these values that yield the best answer; this is
   to keep the values being fit for roughly around 1.0.
   
   i.e., This vector contains the relative activities for the relative efficiency line of y = 1.0
  (independent of energy), except for constrained nuclides; the entries for these will be -1.0.
   */
  vector<double> m_rel_act_norms;
  
  /** The input peak information.
   
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
   There will be one more residual than the number of peaks (with the very last residual being the
   difference of the relative efficiency, at the lowest energy, from 1.0).
#endif
   */
  //std::vector<RelActCalcManual::GenericPeakInfo> m_peak_infos;
  
  /** The input for the manual relative efficiency calculation. */
  RelActCalcManual::RelEffInput m_input;

  /** Warnings from setting up the problem - does not include problems evaluating things. */
  std::vector<std::string> m_setup_warnings;

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
  /** Per-element sigma-block for mass-fraction constraints - the Manual mirror of
   RelActCalcAutoImp::RelActAutoCostFcn::MassFracBlock, sharing the exact decode in
   RelActCalc_imp.hpp (RelActCalc::decode_mass_frac_block; this replaced the former soft-cap).

   The CARRIER (first range-constrained isotope, roster order) parameter holds `t in [0,1]`
   (Manual parameters use offset 0.5, i.e. bounds [0.5, 1.5]) mapping to the range-constrained
   TOTAL fraction `sigma`; remaining range slots hold the `g_k` distribution values.  When every
   isotope of the element is constrained, `sigma` is constant and the carrier slot instead holds
   the elements total relative mass (`scale_multiple * x`, `x >= 0`).
   */
  struct ManualMassFracBlock
  {
    RelActCalc::MassFracBlockSpec spec;

    /** Index into m_isotopes (== parameter index) of the carrier; max() when no carrier is
     needed (a fixed-only element that still has unconstrained isotopes). */
    size_t carrier_index = std::numeric_limits<size_t>::max();

    /** Parameter indices of the g_k distribution values (range isotopes 1..). */
    std::vector<size_t> dist_indices;

    /** Range-constrained isotopes, in block order (carrier first; parallel to spec.lower/upper). */
    std::vector<std::string> range_isos;

    /** Fixed (lower == upper) constrained isotopes, and their pinned fractions. */
    std::vector<std::string> fixed_isos;
    std::vector<double> fixed_fractions;

    /** The elements isotope roster: isotope -> specific activity (from the constraints). */
    std::map<std::string,double> specific_activities;

    /** For `spec.all_constrained` only: element total rel mass per unit of x[carrier_index]. */
    double scale_multiple = 1.0;
  };//struct ManualMassFracBlock

  std::vector<ManualMassFracBlock> m_mass_frac_blocks;
#endif //USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT

  /** just for debug purposes, we'll keep track of how many times the eval function gets called. */
  mutable std::atomic<size_t> m_ncalls;

  /** Profile-likelihood penalty channel state; see `number_residuals()` and #arm_profile.
   Index into `m_input.profile_targets` of the quantity currently being driven, or -1 for
   "disarmed" (the extra residual is then identically zero, so the nominal fit is untouched).

   These are mutable so they can be changed on the const functor Ceres holds - but they are
   written ONLY between `ceres::Solve` calls, never during one: Ceres may evaluate the functor
   from several threads within a single solve.
   */
  mutable int m_profile_active = -1;
  mutable double m_profile_target = 0.0;
  mutable double m_profile_weight = 0.0;
  mutable double m_profile_denom_floor = 0.0;

  /** Point the penalty channel at `target_index`'s quantity, pulling it toward `target_value`
   with weight `weight` (residual `weight*(g(x) - target_value)`).  `denom_floor` smoothly floors
   the quantity's denominator so a vanishing element/denominator activity cannot make it NaN.
   Call with a negative index to disarm.  Not thread safe - see the note above.
   */
  void arm_profile( const int target_index, const double target_value,
                    const double weight, const double denom_floor ) const
  {
    assert( target_index < static_cast<int>(m_input.profile_targets.size()) );
    m_profile_active = target_index;
    m_profile_target = target_value;
    m_profile_weight = weight;
    m_profile_denom_floor = denom_floor;
  }//void arm_profile(...)

  /** Whether reported relative activities carry the LLS-gauge multiple `average_measured_rel_eff`
   - i.e. the Ceres-fit empirical forms, where only the product (curve x activities) is determined
   by the data.

   Seeded from the equation form in the constructor, but `solve_relative_efficiency` OVERWRITES it
   with the decision it actually used: that decision has a fallback (a non-finite average, or a
   peak left with no fitted source activity) which the form alone cannot predict.  Re-deriving it
   here instead would let the profiled quantity apply a gauge multiple the reported activities do
   not, so the two would silently disagree.
   */
  mutable bool m_gauge_normalizes_activities = false;

  bool gauge_normalizes_activities() const
  {
    return m_gauge_normalizes_activities;
  }//bool gauge_normalizes_activities() const

  void set_gauge_normalizes_activities( const bool normalizes ) const
  {
    m_gauge_normalizes_activities = normalizes;
  }


  /** The value of a profile target's quantity at parameters `x`.

   Every form is invariant under the overall activity<->curve scale gauge: the ratio forms because
   the scale cancels, and `RelativeActivity` because it carries the same `average_measured_rel_eff`
   multiple the reported activities do.  All of them go through `relative_activity()`, so
   act-ratio chains and mass-fraction block decodes are automatically respected.
   */
  template<typename T>
  T profile_quantity( const RelActCalcManual::ProfileTarget &target, const std::vector<T> &x,
                      const double denom_floor ) const
  {
    if( target.m_type == RelActCalcManual::ProfileTarget::Type::RelativeActivity )
    {
      const T act = relative_activity( target.m_nuclide, x );
      return gauge_normalizes_activities() ? (average_measured_rel_eff( x ) * act) : act;
    }//if( RelativeActivity )

    if( target.m_type == RelActCalcManual::ProfileTarget::Type::ActivityRatio )
    {
      const T numer = relative_activity( target.m_nuclide, x );
      const T denom = relative_activity( target.m_denom_nuclide, x );
      return numer / smooth_lower_bound( denom, denom_floor );
    }//if( ActivityRatio )

    // Mass fraction of the element roster: (A_t/s_t) / sum_roster(A_i/s_i)
    T total_mass( 0.0 ), target_mass( 0.0 );
    for( const std::map<std::string,double>::value_type &iso_sa : target.m_specific_activities )
    {
      const T mass = relative_activity( iso_sa.first, x ) / iso_sa.second;
      total_mass += mass;
      if( iso_sa.first == target.m_nuclide )
        target_mass = mass;
    }

    return target_mass / smooth_lower_bound( total_mass, denom_floor );
  }//T profile_quantity(...)


  /** Constructor for this functior.
   
   Will throw exception on error.
   */
  ManualGenericRelActFunctor( const RelActCalcManual::RelEffInput &input )
  : 
    m_isotopes{},
    m_rel_act_norms{},
  //m_eqn_form( input.eqn_form ),
  //m_eqn_order( static_cast<int>(input.eqn_order) ),
  //m_peak_infos( input.peaks ),
    m_input( input ),
    m_setup_warnings{},
    m_ncalls( 0 )
  {
    if( m_input.peaks.size() < 1 )
      throw runtime_error( "You must use at least one peak." );

    m_gauge_normalizes_activities = (m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel)
                                     && m_input.use_ceres_to_fit_eqn;

    size_t num_rel_eff_pars_fit = 0;
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      if( m_input.eqn_order != 0 )
        throw runtime_error( "ManualGenericRelActFunctor: equation order must be 0 for FramPhysicalModel." );
      if( !m_input.phys_model_detector )
        throw runtime_error( "ManualGenericRelActFunctor: detector must be specified for FramPhysicalModel." );

      if( m_input.phys_model_use_hoerl )
        num_rel_eff_pars_fit += 2;

      try
      {
        if( m_input.phys_model_self_atten )
        {
          m_input.phys_model_self_atten->check_valid();

          if( shield_is_present(m_input.phys_model_self_atten) )
          {
            if( m_input.phys_model_self_atten->fit_areal_density )
              num_rel_eff_pars_fit += 1;
            if( shield_fits_an(m_input.phys_model_self_atten) )
              num_rel_eff_pars_fit += 1;
          }
        }//if( m_input.phys_model_self_atten )

        for( const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt : m_input.phys_model_external_attens )
        {
          if( !opt )
            throw runtime_error( "ManualGenericRelActFunctor: external attenuation may not be nullptr for FramPhysicalModel." );
          opt->check_valid();

          if( shield_is_present(opt) )
          {
            if( opt->fit_areal_density )
              num_rel_eff_pars_fit += 1;
            if( shield_fits_an(opt) )
              num_rel_eff_pars_fit += 1;
          }
        }//for( loop over m_input.phys_model_external_attens )
      }catch( const std::exception &e )
      {
        throw runtime_error( "ManualGenericRelActFunctor: attenuation input is invalid: " + std::string(e.what()) );
      }
    }else
    {
      // Apply some sanity checks to the eqn_order.  Realistically eqn_order should probably be
      //  between 3 and 6, but we'll allow an arbitrary amount of slop here.
      if( m_input.eqn_order >= 10 )
        throw runtime_error( "ManualGenericRelActFunctor: equation order must be 9 or less." );

#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
      std::sort( begin(m_input.peaks), end(m_input.peaks),
                []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ) -> bool {
        return lhs.m_energy < rhs.m_energy;
        } );
#endif

      num_rel_eff_pars_fit += m_input.eqn_order;  //We actually fit one more paramater than `eqn_order`, but since we force to an average of 1, we effectively have one less.
      if( static_cast<int>(m_input.peaks.size()) < m_input.eqn_order )
        throw runtime_error( "You must use at least as many peaks as the equation order." );

    }//if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else

    const size_t num_peaks = m_input.peaks.size();
    
    // We will check that we arent being passed in a ridiculous number of peak.
    const size_t max_allowed_peaks = 2000; // Arbitrarily chosen value
    if( num_peaks > max_allowed_peaks )
      throw std::runtime_error( "ManualGenericRelActFunctor: too many peaks ("
                                + std::to_string(num_peaks) + ") for the manual rel. eff. solver (max "
                                + std::to_string(max_allowed_peaks) + ")." );
    
    bool unweighted_fit = false;
    
    // We'll go through and make a unique list of isotopes being fit for in m_isotopes, sorting them
    //  alphabetically (alphabetical ordering is arbitrary, other than for the order of isotope
    //  parameters to not be dependent on input ordering; i.e., for debugging purposes)
    //
    //  TODO: we should maybe order the nuclides by a niave guess as to the highest activity nuclide
    //        maybe just take the largest peak, or something...
    for( size_t peak_index = 0; peak_index < m_input.peaks.size(); ++peak_index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = m_input.peaks[peak_index];
      
      if( peak.m_source_gammas.empty() )
        throw std::runtime_error( "ManualGenericRelActFunctor: Peak at " + std::to_string(peak.m_mean)
                                 + " keV has no source gammas defined." );
      
      if( peak.m_counts <= 0.0 )
        throw std::runtime_error( "ManualGenericRelActFunctor: peak counts must be >0.0." );
      
      if( (peak.m_counts_uncert < 0.0)
         || ((peak.m_counts_uncert == 0.0) && (peak.m_base_rel_eff_uncert != -1.0)) )
        throw std::runtime_error( "ManualGenericRelActFunctor: peak count uncertainty must be >0.0." );
      
      if( peak.m_base_rel_eff_uncert == -1.0 )
      {
        if( peak_index && !unweighted_fit )
          throw runtime_error( "ManualGenericRelActFunctor: RelActCalcManual::GenericPeakInfo::m_base_rel_eff_uncert must consistently either have a value [0,1], or -1" );
        
        unweighted_fit = true;
      }else if( unweighted_fit || (peak.m_base_rel_eff_uncert < 0.0) || (peak.m_base_rel_eff_uncert > 1.0) )
      {
        throw runtime_error( "ManualGenericRelActFunctor: RelActCalcManual::GenericPeakInfo::m_base_rel_eff_uncert must consistently either have a value [0,1], or -1" );
      }//
      
      set<string> nuclides_in_peak;
      
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        const std::string &iso = line.m_isotope;
        
        // We probably wont usually have more than one line from a nuclide for a peak -- but we
        //  could so we'll just print out a warning, rather than throwing a hard error or something.
        if( nuclides_in_peak.count(iso) )
        {
          m_setup_warnings.push_back( "peak contained multiple lines from same nuclide." );
          std::cerr << "Warning: ManualGenericRelActFunctor: peak contained multiple lines from same nuclide.\n";
        }
        nuclides_in_peak.insert( iso );
        
        // lower_bound returns iterator to the that is not less than (i.e. greater or equal to)
        const auto pos = std::lower_bound( std::begin(m_isotopes), std::end(m_isotopes), iso );
        
        if( (pos == std::end(m_isotopes)) || ((*pos) != iso) )
          m_isotopes.insert( pos, iso );
        
        //  TODO: look at the various decay chains to determine better minimal yield value.
        //    Note that numeric_limits<float> is 1.19209e-07 on my system, which the BR for
        //    Pu that are important can be significantly smaller than this - so for the moment we'll
        //    essentially use zero.
        const double min_allowable_yield = std::numeric_limits<float>::min();
        if( line.m_yield <= min_allowable_yield )
          throw std::runtime_error( "ManualGenericRelActFunctor: yields must be greater than zero." );
      }//for( loop over nuclides that contribute to this peak )
    }//for( loop over fit peaks )
    
    const size_t num_isotopes = m_isotopes.size();
    
    // We should be guaranteed to have at least one isotope by here, but we'll throw in a check to
    //  avoid static analysis warnings, or jic, or something
    assert( num_isotopes );
    if( num_isotopes < 1 )
      throw std::runtime_error( "ManualGenericRelActFunctor: no isotopes specified." );
  
    // We'll first estimate the relative activities using a flat line
    vector<double> dummy_rel_act_norm_uncerts;
    const RelActCalc::RelEffEqnForm est_eqn_form = m_input.eqn_form;
    
    // If we have no act ratio constraints, we can just fit the relative activities directly,
    //  otherwise we need to re-define the isotopes and BRs to eliminate constrained isotopes.
    if( m_input.act_ratio_constraints.empty() 
#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
        && m_input.mass_fraction_constraints.empty() 
#endif // USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT
    )
    {
      if( est_eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        RelActCalcManual::fit_act_to_phys_rel_eff( m_input, m_isotopes, m_input.peaks,
                          m_rel_act_norms, dummy_rel_act_norm_uncerts );
      }else
      {
        const vector<double> flat_rel_eff_coefs{ (est_eqn_form==RelActCalc::RelEffEqnForm::LnX ? 1.0 : 0.0), 0.0 };
        RelActCalcManual::fit_act_to_rel_eff( est_eqn_form, flat_rel_eff_coefs,
                         m_isotopes, m_input.peaks,
                         m_rel_act_norms, dummy_rel_act_norm_uncerts );
      }
    }else
    {
      // If we have act ratio constraints, we will re-define the isotopes of the
      //  peaks, so that only non-constrained isotopes are fit for; 
      vector<string> mod_isotopes = m_isotopes;
      vector<RelActCalcManual::GenericPeakInfo> mod_peaks = m_input.peaks;
      
      for( RelActCalcManual::GenericPeakInfo &peak : mod_peaks )
      {
        for( RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
        {
          size_t num_iter = 0; //I think the logic is solid - but add in a protection anyway.
          double act_scale_factor = 1.0;
          string parent_iso, line_iso = line.m_isotope;

          while( parent_iso != line_iso )
          {
            parent_iso = line_iso;

            for( const RelActCalcManual::ManualActRatioConstraint &c : m_input.act_ratio_constraints )
            {
              if( c.m_constrained_nuclide == line_iso )
              {
                act_scale_factor *= c.m_constrained_to_controlled_activity_ratio;
                line_iso = c.m_controlling_nuclide;
                break;
              }
            }//for( const ManualActRatioConstraint &c : m_input.act_ratio_constraints )

            ++num_iter;
            assert( num_iter < 100 );
            if( num_iter > 100 )
              throw runtime_error( "ManualGenericRelActFunctor: act ratio constraint loop." );
          }//while( parent_iso != line_iso )
          
          assert( parent_iso == line_iso );

          line.m_yield *= act_scale_factor;
          line.m_isotope = parent_iso;
        }//for( GenericLineInfo &line : peak.m_source_gammas )

        // Later on there is a check that each peak only has each source a single time,
        //  we will lump them together here.
        std::map<string,RelActCalcManual::GenericLineInfo> line_map;
        for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
        {
          auto pos = line_map.find(line.m_isotope);
          if( pos != line_map.end() )
            pos->second.m_yield += line.m_yield;
          else
            line_map[line.m_isotope] = line;
        }
        peak.m_source_gammas.clear();
        for( const auto &[iso, line] : line_map )
          peak.m_source_gammas.push_back( line );
      }//for( GenericPeakInfo &peak : mod_peaks )

      for( const RelActCalcManual::ManualActRatioConstraint &constraint : m_input.act_ratio_constraints )
      {
        const auto pos = std::find( begin(mod_isotopes), end(mod_isotopes), constraint.m_constrained_nuclide );
        if( pos != std::end(mod_isotopes) )
          mod_isotopes.erase( pos );
      }

    vector<double> mod_activities;

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
      // Mass-fraction constraints solve through the per-element sigma-block (see
      //  #ManualMassFracBlock and RelActCalc::decode_mass_frac_block) - the fit itself handles
      //  fixed and range constraints exactly, including all-constrained elements.  This PRE-fit
      //  estimate stays simpler: range-constrained isotopes are fit freely (their windows are
      //  ignored for the start), and fixed ones are excluded below.
      if( !m_input.mass_fraction_constraints.empty() )
      {
        for( const RelActCalcManual::MassFractionConstraint &constraint : m_input.mass_fraction_constraints )
        {
          // If we are not fixing the mass-fraction to a certain value, then we will freely fit things, totally
          //  ignoring the constraints - again, not ideal, but this is just a first guess.
          if( constraint.m_mass_fraction_lower != constraint.m_mass_fraction_upper )
            continue;

          const auto mod_isotope_pos = std::find( begin(mod_isotopes), end(mod_isotopes), constraint.m_nuclide );
          assert( mod_isotope_pos != end(mod_isotopes) );
          if( mod_isotope_pos != std::end(mod_isotopes) )
            mod_isotopes.erase( mod_isotope_pos );

          for( RelActCalcManual::GenericPeakInfo &peak : mod_peaks )
          {
            // Remove all lines from the peak that are the mass-fraction-constrained nuclide.
            peak.m_source_gammas.erase( std::remove_if( begin(peak.m_source_gammas), end(peak.m_source_gammas),
                            [&constraint]( const RelActCalcManual::GenericLineInfo &line ) -> bool {
              return line.m_isotope == constraint.m_nuclide;
            } ), end(peak.m_source_gammas) );
          }//for( RelActCalcManual::GenericPeakInfo &peak : mod_peaks )
        }//for( RelActCalcManual::MassFractionConstraint &constraint : m_input.mass_fraction_constraints )

        mod_peaks.erase( std::remove_if( begin(mod_peaks), end(mod_peaks),
                            []( const RelActCalcManual::GenericPeakInfo &peak ) -> bool {
              return peak.m_source_gammas.empty();
            } ), end(mod_peaks) );
      }//if( !m_input.mass_fraction_constraints.empty() )
#endif  //USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT

      // Every isotope may have been excluded from the free pre-fit (e.g., an all-FIXED
      //  mass-fraction-constrained element, where only the element total is fit) - then there is
      //  nothing to estimate here, and the LLS below would reduce an empty matrix; the norms
      //  simply fall back to 1.0 in the loop below.
      if( !mod_isotopes.empty() )
      {
        if( est_eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
        {
          RelActCalcManual::fit_act_to_phys_rel_eff( m_input, mod_isotopes, mod_peaks,
                              mod_activities, dummy_rel_act_norm_uncerts );
        }else
        {
          const vector<double> flat_rel_eff_coefs{ (est_eqn_form==RelActCalc::RelEffEqnForm::LnX ? 1.0 : 0.0), 0.0 };
          RelActCalcManual::fit_act_to_rel_eff( est_eqn_form, flat_rel_eff_coefs,
                             mod_isotopes, mod_peaks,
                             mod_activities, dummy_rel_act_norm_uncerts );
        }
      }//if( !mod_isotopes.empty() )

      assert( mod_activities.size() == mod_isotopes.size() );
      m_rel_act_norms.resize( m_isotopes.size() );
      for( size_t i = 0; i < m_isotopes.size(); ++i )
      {
        const string &iso = m_isotopes[i];
        const auto pos = std::find( begin(mod_isotopes), end(mod_isotopes), iso );
        if( pos == end(mod_isotopes) )  
        {
          // This nuclide was removed from the free fit because it is constrained.  Act-ratio-controlled
          //  nuclides keep the -1.0 sentinel; a *fixed* (lower==upper) mass-fraction constraint is also removed
          //  here, but takes a norm of 1.0 - the same value range mass-fraction constraints get, and the value
          //  release already settles on via the normalization loop below - so the norm invariant in
          //  mass_fraction() stays consistent.  (Range constraints are not removed from mod_isotopes, so any
          //  mass-fraction-constrained nuclide reaching this branch is necessarily a fixed one.)
          bool is_mass_constrained = false;
#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
          for( size_t j = 0; !is_mass_constrained && (j < m_input.mass_fraction_constraints.size()); ++j )
            is_mass_constrained = (m_input.mass_fraction_constraints[j].m_nuclide == iso);
#endif
          m_rel_act_norms[i] = is_mass_constrained ? 1.0 : -1.0;
#ifndef NDEBUG
          bool has_constrained = is_mass_constrained;
          for( size_t j = 0; !has_constrained && (j < m_input.act_ratio_constraints.size()); ++j )
            has_constrained = (m_input.act_ratio_constraints[j].m_constrained_nuclide == iso);
          assert( has_constrained );
#endif
        }else
        {
          m_rel_act_norms[i] = mod_activities[pos - begin(mod_isotopes)];
#ifndef NDEBUG
          bool is_controlled = false;
          for( size_t j = 0; !is_controlled && (j < m_input.act_ratio_constraints.size()); ++j )
            is_controlled = (m_input.act_ratio_constraints[j].m_constrained_nuclide == iso);
          assert( !is_controlled );

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
          bool is_mass_constrained = false;
          for( size_t j = 0; !is_mass_constrained && (j < m_input.mass_fraction_constraints.size()); ++j )
            is_mass_constrained = (m_input.mass_fraction_constraints[j].m_nuclide == iso);
          if( is_mass_constrained )
            m_rel_act_norms[i] = 1.0;
#endif
#endif
        }//if( pos == end(mod_isotopes) ) / else
      }//for( size_t i = 0; i < m_isotopes.size(); ++i )
    }//if( !m_input.act_ratio_constraints.empty() )

    assert( m_isotopes.size() == m_rel_act_norms.size() );

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
    // Group the mass-fraction constraints into per-element sigma-blocks (see #ManualMassFracBlock).
    //  The roster (m_specific_activities) identifies an element; check_nuclide_constraints() already
    //  validated the rosters are consistent and the windows feasible.
    for( const RelActCalcManual::MassFractionConstraint &constraint : m_input.mass_fraction_constraints )
    {
      const auto already_pos = std::find_if( begin(m_mass_frac_blocks), end(m_mass_frac_blocks),
        [&constraint]( const ManualMassFracBlock &b ) -> bool {
          return b.specific_activities.count( constraint.m_nuclide );
        } );
      if( already_pos != end(m_mass_frac_blocks) )
        continue; //this constraints element already has its block

      ManualMassFracBlock block;
      block.specific_activities = constraint.m_specific_activities;

      size_t num_constrained = 0;
      std::vector<std::pair<double,double>> windows;
      for( const std::map<std::string,double>::value_type &iso_sa : block.specific_activities )
      {
        const string &iso = iso_sa.first;
        const auto c_pos = std::find_if( begin(m_input.mass_fraction_constraints), end(m_input.mass_fraction_constraints),
          [&iso]( const RelActCalcManual::MassFractionConstraint &c ) -> bool { return c.m_nuclide == iso; } );
        if( c_pos == end(m_input.mass_fraction_constraints) )
          continue;

        num_constrained += 1;
        if( c_pos->m_mass_fraction_lower == c_pos->m_mass_fraction_upper )
        {
          block.fixed_isos.push_back( iso );
          block.fixed_fractions.push_back( c_pos->m_mass_fraction_lower );
        }else
        {
          block.range_isos.push_back( iso );
          windows.push_back( {c_pos->m_mass_fraction_lower, c_pos->m_mass_fraction_upper} );
        }
      }//for( const auto &iso_sa : block.specific_activities )

      double fixed_sum = 0.0;
      for( const double f : block.fixed_fractions )
        fixed_sum += f;

      const bool all_constrained = (num_constrained >= block.specific_activities.size());

      block.spec = RelActCalc::make_mass_frac_block_spec( windows, fixed_sum, all_constrained );

      if( block.range_isos.empty() )
      {
        assert( !block.fixed_isos.empty() );
        if( all_constrained )
          block.carrier_index = iso_index( block.fixed_isos[0] );
      }else
      {
        block.carrier_index = iso_index( block.range_isos[0] );
        for( size_t k = 1; k < block.range_isos.size(); ++k )
          block.dist_indices.push_back( iso_index( block.range_isos[k] ) );
      }

      // TODO: for all-constrained elements a better `scale_multiple` (element rel mass) start
      //  could come from the pre-fit activities; the solver refines it from 1.0 fine so far.
      m_mass_frac_blocks.push_back( std::move(block) );
    }//for( const RelActCalcManual::MassFractionConstraint &constraint : ... )
#endif //USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT

    size_t num_fit_activities = 0;
    for( size_t i = 0; i < m_rel_act_norms.size(); ++i )
    {
      bool is_controlled = false;
      for( size_t j = 0; !is_controlled && (j < m_input.act_ratio_constraints.size()); ++j )
        is_controlled = (m_input.act_ratio_constraints[j].m_constrained_nuclide == m_isotopes[i]);
      assert( !is_controlled || (m_rel_act_norms[i] == -1.0) );
      if( is_controlled )
        m_rel_act_norms[i] = -1.0;// JIC

      bool is_mass_constrained = false; //declared outside the #if guard - also used below it
#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
      for( size_t j = 0; !is_mass_constrained && (j < m_input.mass_fraction_constraints.size()); ++j )
        is_mass_constrained = (m_input.mass_fraction_constraints[j].m_nuclide == m_isotopes[i]);
      assert( !is_mass_constrained || (m_rel_act_norms[i] == 1.0) );
      if( is_mass_constrained )
        m_rel_act_norms[i] = 1.0; // JIC
#endif

      // Only replace unusable (non-positive / non-finite) initial estimates: a legitimately
      //  small pre-fit activity (CPS-scaled data, trace isotopes) is a good starting scale, and
      //  flooring it at 1.0 used to leave that nuclides Ceres parameter far from O(1).
      if( !is_controlled && ((m_rel_act_norms[i] <= 0.0) || IsInf(m_rel_act_norms[i]) || IsNan(m_rel_act_norms[i])) )
      {
        m_setup_warnings.push_back( "The initial activity estimate for " + m_isotopes[i]
                                    + " was " + std::to_string(m_rel_act_norms[i])
                                    + ", so will use 1.0 instead.");
        m_rel_act_norms[i] = 1.0;
      }//if( initial activity estimate is unusable )

      if( !is_mass_constrained && !is_controlled )
        num_fit_activities += 1;
    }//for( size_t i = 0; i < m_rel_act_norms.size(); ++i )

    //assert( mod_activities.size() == num_fit_activities );

    if( (num_fit_activities + num_rel_eff_pars_fit) > m_input.peaks.size() )
      throw std::runtime_error( "ManualGenericRelActFunctor: you must have at least as many peaks as"
                               " parameters you are fitting for." );
  }//ManualGenericRelActFunctor constructor
  
  size_t number_residuals() const
  {
    // One extra channel, present for the life of the functor, when profile-likelihood intervals
    //  were asked for; it is inert (identically zero) until armed between solves.  Ceres needs a
    //  constant residual count, hence "allocate up front, arm later".
    const size_t num_profile = m_input.profile_targets.empty() ? 0 : 1;

#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      return m_input.peaks.size() + num_profile;
    return m_input.peaks.size() + 1 + num_profile;
#else
    return m_input.peaks.size() + num_profile;
#endif
  }//size_t number_residuals() const
  
  size_t iso_index( const std::string &iso ) const
  {
    const auto iso_pos = std::lower_bound( std::begin(m_isotopes), std::end(m_isotopes), iso );
    assert( (iso_pos != std::end(m_isotopes)) && ((*iso_pos) == iso) );

    if( (iso_pos == std::end(m_isotopes)) || ((*iso_pos) != iso) )
      throw std::logic_error( "ManualGenericRelActFunctor: missing nuclide '" + iso + "'" );

    return static_cast<size_t>( iso_pos - std::begin(m_isotopes) );
  }
  
  template<typename T>
  T relative_activity( const std::string &iso, const vector<T> &x ) const
  {
    const size_t index = iso_index( iso );
    assert( index < x.size() );

    int constraint_index = -1;
    for( size_t i = 0; ((constraint_index < 0) && (i < m_input.act_ratio_constraints.size())); ++i )  
      if( m_input.act_ratio_constraints[i].m_constrained_nuclide == iso )
        constraint_index = static_cast<int>(i);

    if( constraint_index >= 0 )
    {
#ifndef NDEBUG
      const size_t index = iso_index( iso );
      assert( (abs(x[index] - -1.0) < 1.0E-6) || (abs(x[index] - 0.0) < 1.0E-6) );
#endif

      const RelActCalcManual::ManualActRatioConstraint &constraint = m_input.act_ratio_constraints[constraint_index];
      return relative_activity( constraint.m_controlling_nuclide, x ) * constraint.m_constrained_to_controlled_activity_ratio;
    }//if( constraint_index >= 0 )

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
    for( const ManualMassFracBlock &block : m_mass_frac_blocks )
    {
      // Is `iso` one of this blocks constrained isotopes (fixed or range)?
      int range_pos = -1;
      bool is_fixed = false;
      T this_frac( 0.0 );

      for( size_t i = 0; !is_fixed && (i < block.fixed_isos.size()); ++i )
      {
        if( block.fixed_isos[i] == iso )
        {
          is_fixed = true;
          this_frac = T( block.fixed_fractions[i] ); //fixed constraints decode to exactly their pinned fraction
        }
      }//for( loop over fixed-constrained isotopes )

      for( size_t i = 0; !is_fixed && (range_pos < 0) && (i < block.range_isos.size()); ++i )
      {
        if( block.range_isos[i] == iso )
          range_pos = static_cast<int>( i );
      }

      if( !is_fixed && (range_pos < 0) )
        continue; //`iso` is not part of this elements block

      const auto this_sa_pos = block.specific_activities.find( iso );
      assert( this_sa_pos != end(block.specific_activities) );
      if( this_sa_pos == end(block.specific_activities) )
        throw std::logic_error( "ManualGenericRelActFunctor: missing nuclide in constraint.m_specific_activities???" );

      const size_t num_range = block.range_isos.size();

      // The range-constrained total mass fraction, from the carrier parameter; a hard box bound
      //  keeps `fixed_sum + sigma <= 1 - delta` (this replaced the former soft-cap decode, whose
      //  0.95-of-budget knee compressed the top of every window - see RelActCalc::MassFracBlockSpec).
      T sigma( block.spec.sig_lo );
      if( block.spec.all_constrained )
      {
        sigma = T( block.spec.sig_hi ); //== 1 - fixed_sum; the carrier slot holds the element scale instead
      }else if( num_range > 0 )
      {
        const T t = x[block.carrier_index] - 0.5; //box-bounded to [0.5, 1.5]
        assert( (t > -0.02) && (t < 1.02) ); //leave some room for numerical differentiation
        sigma = block.spec.sig_lo + t*(block.spec.sig_hi - block.spec.sig_lo);
      }

      if( range_pos >= 0 )
      {
        // Decode all the blocks range fractions through the shared exact sigma-block decode.
        vector<T> gs( (num_range > 1) ? (num_range - 1) : size_t(0), T(0.0) );
        vector<T> fractions( num_range, T(0.0) );
        for( size_t k = 1; k < num_range; ++k )
        {
          gs[k-1] = x[ block.dist_indices[k-1] ] - 0.5; //box-bounded to [0.5, 1.5]
          assert( (gs[k-1] > -0.02) && (gs[k-1] < 1.02) );
        }

        RelActCalc::decode_mass_frac_block( block.spec, sigma, gs.data(), fractions.data() );
        this_frac = fractions[range_pos];
      }//if( range_pos >= 0 )

      if( block.spec.all_constrained )
      {
        // Every isotope of the element is constrained: no unconstrained isotopes carry the
        //  elements absolute scale, so the carrier slot holds the total relative mass (x >= 0).
        const T el_rel_mass = block.scale_multiple * x[block.carrier_index];
        return el_rel_mass * this_frac * this_sa_pos->second;
      }//if( block.spec.all_constrained )

      // Sum the relative masses of the elements unconstrained isotopes (rel act / specific
      //  activity); together they make up the `1 - sum_constrained` remainder of the element
      //  mass - strictly positive, by the hard bound on the carrier.
      const T sum_constrained = T(block.spec.fixed_sum) + sigma;
      T sum_unconstrained_rel_mass_of_el( 0.0 );
      for( const std::map<std::string, double>::value_type &iso_sa : block.specific_activities )
      {
        if( std::count( begin(block.fixed_isos), end(block.fixed_isos), iso_sa.first )
            || std::count( begin(block.range_isos), end(block.range_isos), iso_sa.first ) )
          continue; //constrained isotopes are folded into sum_constrained

        sum_unconstrained_rel_mass_of_el += (relative_activity( iso_sa.first, x ) / iso_sa.second);
      }//for( const auto &iso_sa : block.specific_activities )

      const T unconstrained_rel_mass_frac_of_el = 1.0 - sum_constrained;
      const T total_rel_mass = sum_unconstrained_rel_mass_of_el / unconstrained_rel_mass_frac_of_el;
      const T this_rel_mass = total_rel_mass * this_frac;

      return this_rel_mass * this_sa_pos->second;
    }//for( const ManualMassFracBlock &block : m_mass_frac_blocks )
#endif


    return m_rel_act_norms[index] * x[index];
  }
  
  
  template<typename T>
  struct PhysModelRelEqnDef
  {
    shared_ptr<const DetectorPeakResponse> det;
    
    std::optional<RelActCalc::PhysModelShield<T>> self_atten;
    std::vector<RelActCalc::PhysModelShield<T>> external_attens;
    
    std::optional<T> hoerl_b, hoerl_c;
  };//struct PhysModelRelEqnDef
  
  
  template<typename T>
  static PhysModelRelEqnDef<T> make_phys_eqn_input( const RelActCalcManual::RelEffInput &input,
                                         const std::vector<T> &eqn_coefficients )
  {
    // See `setup_physical_model_shield_par_manual(...)` for details of equation coefficients, but
    //  the short of it is there will only be a parameter for AN if material is nullptr, and
    //  atomic number is being fit for.
    // AD will always have a parameter dedicated to it.
    // if a RelActCalc::PhysicalModelShieldInput is nullptr, or material is null and
    //  opt->atomic_number == 0.0, and fit AN is false, then there will be no parameters
    //  dedicated to the PhysicalModelShieldInput.
    //  (this design is a bit vestigial, and we *could* change things to always have a
    //   consistent parameter definition, but maybe its worth waiting until we make another
    //   change to do this)
    
    assert( input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel );
    if( input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      throw std::logic_error( "make_phys_eqn_input called for non-physical model." );
    
    PhysModelRelEqnDef<T> answer;
    
    answer.det = input.phys_model_detector;
    
    size_t shield_index = 0;
    if( shield_is_present( input.phys_model_self_atten ) )
    {
      RelActCalc::PhysModelShield<T> self_atten;

      self_atten.material = input.phys_model_self_atten->material;
      if( !self_atten.material )
      {
        if( input.phys_model_self_atten->fit_atomic_number )
        {
          assert( (shield_index + 1) <= eqn_coefficients.size() );
          if( (shield_index + 1) > eqn_coefficients.size() )
            throw logic_error( "make_phys_eqn_input: not enough input coefficients (1)" );
          
          const T an = eqn_coefficients[shield_index] * RelActCalc::ns_an_ceres_mult;
          shield_index += 1;
          assert( an >= 1.0 && an <= 98.0 );
          self_atten.atomic_number = an;
        }else
        {
          self_atten.atomic_number = T( input.phys_model_self_atten->atomic_number );
        }
      }//if( !self_atten.material )
      
      assert( (shield_index + 1) <= eqn_coefficients.size() );
      if( (shield_index + 1) > eqn_coefficients.size() )
        throw logic_error( "make_phys_eqn_input: not enough input coefficients (2)" );
      
      assert( eqn_coefficients[shield_index] >= 0.0 && eqn_coefficients[shield_index] <= 500.0 );
      self_atten.areal_density = eqn_coefficients[shield_index] * PhysicalUnits::g_per_cm2;
      shield_index += 1;
      
      answer.self_atten = std::move(self_atten);
    }//m_options.phys_model_self_atten
      
      
    for( size_t i = 0; i < input.phys_model_external_attens.size(); ++i )
    {
      const auto &a = input.phys_model_external_attens[i];
      if( !shield_is_present(a) )
        continue;

      RelActCalc::PhysModelShield<T> atten;
      atten.material = a->material;
      if( !atten.material )
      {
        if( a->fit_atomic_number )
        {
          assert( (shield_index + 1) <= eqn_coefficients.size() );
          if( (shield_index + 1) > eqn_coefficients.size() )
            throw logic_error( "make_phys_eqn_input: not enough input coefficients (3)" );
          
          const T an = eqn_coefficients[shield_index]  * RelActCalc::ns_an_ceres_mult;
          assert( (an >= 1.0) && (an <= 98.0) );
          atten.atomic_number = an;
          shield_index += 1;
        }else
        {
          atten.atomic_number = T( a->atomic_number );
        }
      }//if( !atten.material )
      
      assert( (shield_index + 1) <= eqn_coefficients.size() );
      if( (shield_index + 1) > eqn_coefficients.size() )
        throw logic_error( "make_phys_eqn_input: not enough input coefficients (4)" );
      
      assert( eqn_coefficients[shield_index] >= 0.0 && eqn_coefficients[shield_index] <= 500.0 );
      atten.areal_density = eqn_coefficients[shield_index] * PhysicalUnits::g_per_cm2;
      shield_index += 1;
      
      answer.external_attens.push_back( std::move(atten) );
    }//for( loop over input.phys_model_external_attens )
      
    if( input.phys_model_use_hoerl )
    {
      assert( (shield_index + 2) <= eqn_coefficients.size() );
      if( (shield_index + 2) > eqn_coefficients.size() )
        throw logic_error( "make_phys_eqn_input: not enough input coefficients (5)" );
      
      answer.hoerl_b = (eqn_coefficients[shield_index] - RelActCalc::ns_decay_hoerl_b_offset) * RelActCalc::ns_decay_hoerl_b_multiple;
      shield_index += 1;
      answer.hoerl_c = (eqn_coefficients[shield_index] - RelActCalc::ns_decay_hoerl_c_offset) * RelActCalc::ns_decay_hoerl_c_multiple;
      shield_index += 1;
    }//if( input.phys_model_use_hoerl )
    
    assert( shield_index == eqn_coefficients.size() );
    if( shield_index != eqn_coefficients.size() )
      throw logic_error( "make_phys_eqn_input: number of equation coefficients mismatch" );
    
    return answer;
  }//PhysModelRelEqnDef make_phys_eqn_input( const std::vector<double> &eqn_coefficients ) const
  
  
  template<typename T>
  std::function<T(double)> make_rel_eff_fcn( const std::vector<T> &eqn_coefficients ) const
  {
    const RelActCalc::RelEffEqnForm eqn_form = m_input.eqn_form;
    if( eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      // Note: will take copy of `eqn_coefficients` to avoid potential life-time issues, but we
      //       could maybe instead take a reference, for all our current use cases.
      return [eqn_form,eqn_coefficients]( double energy ){
        return eval_eqn_imp( energy, eqn_form, eqn_coefficients.data(), eqn_coefficients.size() );
      };
    }
      
    PhysModelRelEqnDef input = make_phys_eqn_input( m_input, eqn_coefficients );
    
    
    return [input]( double energy ){
      return RelActCalc::eval_physical_model_eqn_imp<T>( energy, input.self_atten, input.external_attens,
                                            input.det.get(), input.hoerl_b, input.hoerl_c );
    };
  }//std::function<T(double)> rel_eff_fcn( const std::vector<T> &x ) const


  /** The plain average, over all peaks, of the measured relative efficiency implied by the
   activities `x`: (1/P) * sum_p( C_p / sum_i( A_i(x)*y_pi ) ).

   This is the normalization `fit_rel_eff_eqn_lls_imp(...)` divides out of the measured rel-eff
   points in the LLS fit mode; it is used to place Ceres-fit empirical-form solutions in that
   same gauge (average measured rel. eff. == 1), so the two fit modes report directly comparable
   activities.  Note m(x)*A_i(x) is exactly invariant under the scale orbit (activities times k,
   curve divided by k).

   The denominator carries the same smooth floor `fit_rel_eff_eqn_lls_imp` applies: this is not
   only evaluated once at a converged solution - the profile-likelihood penalty channel calls it
   for `Type::RelativeActivity` targets, so it sees whatever trial point the optimizer probes,
   including ones where a source's activity has been driven to its lower bound of zero.  Without
   the floor that is a division by zero, i.e. a non-finite residual.
   */
  template<typename T>
  T average_measured_rel_eff( const std::vector<T> &x, bool *any_peak_floored = nullptr ) const
  {
    assert( !m_input.peaks.empty() );

    if( any_peak_floored )
      *any_peak_floored = false;

    T sum( 0.0 );
    for( const RelActCalcManual::GenericPeakInfo &peak : m_input.peaks )
    {
      T rel_src_counts( 0.0 );
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
        rel_src_counts += relative_activity( line.m_isotope, x ) * line.m_yield;

      const double floor_val = 1.0E-6 * (std::max)( 1.0, peak.m_counts );

      // Report, but do not hide, a floored peak.  The floor keeps this evaluable, but a peak whose
      //  only source has been driven to zero contributes ~1e6 rather than the infinity it used to,
      //  which is finite enough to sail through a caller's isinf() check while being just as
      //  meaningless - see where the reporting gauge is computed.
      if( any_peak_floored && (rel_src_counts < 10.0*floor_val) )
        *any_peak_floored = true;

      sum += peak.m_counts / smooth_lower_bound( rel_src_counts, floor_val );
    }//for( loop over peaks )

    return sum / static_cast<double>( m_input.peaks.size() );
  }//T average_measured_rel_eff( const std::vector<T> &x, bool *any_peak_floored ) const


  template<typename T>
  void eval_internal_lls_rel_eff( const std::vector<T> &x, T *residuals ) const
  {
    using namespace std;
    using namespace ceres;
    
    m_ncalls += 1;
    
    assert( residuals );
    assert( m_input.eqn_order >= 0 );
    assert( !m_input.use_ceres_to_fit_eqn );
    assert( x.size() == m_isotopes.size() );
    assert( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel );
    
    const size_t num_eqn_pars = static_cast<size_t>( m_input.eqn_order + 1 );
    const size_t num_isotopes = m_isotopes.size();
    assert( num_isotopes >= 1 );
    
    const T * const pars = x.data();
    vector<T> rel_activities( num_isotopes );
    for( size_t i = 0; i < num_isotopes; ++i )
    {
      rel_activities[i] = this->relative_activity(m_isotopes[i], x);
      if( isnan(rel_activities[i]) || isinf(rel_activities[i]) )
        cerr << "Got inf/nan rel act for iso " << i << endl;
      assert( !isnan(rel_activities[i]) && !isinf(rel_activities[i]) );
    }

    vector<T> eqn_coefficients;
    // Note: we deliberately pass nullptr for the covariance: it is unused here, and requesting it
    //  would force a (ceres::Jet) matrix pseudo-inverse on every residual evaluation.  The final
    //  covariance is computed once, in double precision, where the solution is extracted.
    fit_rel_eff_eqn_lls_imp<T>( m_input.eqn_form, m_input.eqn_order,
                                m_isotopes,
                                rel_activities,
                                m_input.peaks,
                                eqn_coefficients, nullptr );

    assert( eqn_coefficients.size() == (m_input.eqn_order + 1) );
    
    
    std::function<T(double)> rel_eff_curve = make_rel_eff_fcn<T>( eqn_coefficients );
    
    
    for( size_t index = 0; index < m_input.peaks.size(); ++index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = m_input.peaks[index];
      
      const T curve_val = rel_eff_curve( peak.m_energy );
      
      if( isnan(curve_val) || isinf(curve_val) )
        throw std::runtime_error( "Perpective Rel. Eff. curve value invalid at " + std::to_string(peak.m_energy) + " keV" );
      
      T rel_src_counts( 0.0 );
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        const T rel_activity = relative_activity( line.m_isotope, x );
        rel_src_counts += rel_activity * line.m_yield;
      }//for( const GenericLineInfo &line : peak.m_source_gammas )
      
      if( isnan(rel_src_counts) || isinf(rel_src_counts) )
        throw std::runtime_error( "Perpective source counts value invalid at " + std::to_string(peak.m_energy) + " keV" );
      
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
      if( (index == 0) && (m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel) )
      {
        // Pin the first relative efficiency to 1.0; this is to prevent degeneracy of LM and LLS
        //  both essentially fitting for the overall normalization
        const double weight = 10.0; //Fairly arbitrary - not sure what to actually use - need something so the uncerts on activities dont get huge (which they will with no weight)
        residuals[m_peak_infos.size()] = weight*((peak.m_counts / rel_src_counts) - 1.0);
      }
#endif
      
      
      // TODO: should we fold in the uncertainty from the relative efficiency equation into the below as well - I dont think so since we are treating them as nuisance parameters?
      
      // Note: we want to compute:
      //  `(rel_efficiency - curve_val) / rel_eff_uncertainty`,
      //  But this has some divide by zero issues when rel act is zero, but the below gives the same
      //  residual, but avoiding these issues.
      if( peak.m_base_rel_eff_uncert == -1.0 )
      {
        // We are doing an unweighted fit
        
        // Avoid dividing by zero, so make sure rel_src_counts isnt really close to zero.
        // TODO: - need to evaluate using a different formulation where \c rel_src_counts being zero isnt a problem (I think its fine, but need to check before making the change - and make sure it wont effect how we use the covariances).
        if( ((rel_src_counts < 1.0E-8) && (rel_src_counts < (1.0E-6*peak.m_counts)))
           || (rel_src_counts < numeric_limits<T>::epsilon()) ) //keep consistent with the eval_internal_nl_rel_eff twin
        {
          rel_src_counts = T(1.0E-6) * peak.m_counts;
        }
        
        residuals[index] = (peak.m_counts / rel_src_counts) - curve_val;
      }else if( peak.m_base_rel_eff_uncert == 0.0 )
      {
        // We are not using a m_base_rel_eff_uncert value
        const T pred_counts = curve_val * rel_src_counts;
        residuals[index] = (peak.m_counts - pred_counts) / peak.m_counts_uncert;
      }else
      {
        // We are using a m_base_rel_eff_uncert value
        assert( peak.m_base_rel_eff_uncert <= 1.0 );
        
        const T pred_counts = curve_val * rel_src_counts;
        // Note: for `add_uncert` below, we are using peak.m_counts, but it *could* also be
        //       reasonable to use `rel_src_counts` (which I was doing pre 20220720) or
        //       even `pred_counts`.  This is maybe worth revisiting.  Note that if you change
        //       how things are calculated here, you should also be consistent with
        //       #fit_rel_eff_eqn_lls.
        const double add_uncert = peak.m_counts * peak.m_base_rel_eff_uncert;
        const double uncert = sqrt( pow(peak.m_counts_uncert,2.0) + pow(add_uncert,2.0) );
        
        residuals[index] = (peak.m_counts - pred_counts) / uncert;
        //cout << "residuals[index] = (peak.m_counts - pred_counts) / uncert --> " << residuals[index] << " = (" << peak.m_counts << " - " << pred_counts << ") / " << uncert << endl;
        //cout << "curve_val=" << curve_val << ", rel_src_counts=" << rel_src_counts << endl;
      }
    }//for( loop over energies to evaluate at )
  }//eval_internal_lls_rel_eff(...)


  template<typename T>
  void eval_internal_nl_rel_eff( const std::vector<T> &x, T *residuals ) const
  {
    using namespace std;
    using namespace ceres;
    
    m_ncalls += 1;
    
    assert( residuals );

    const size_t num_isos = num_isotopes();
  
    std::function<T (double)> rel_eff_fcn;
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      size_t par_num = num_isos;
      
      assert( m_input.phys_model_detector );
      const DetectorPeakResponse * const det = m_input.phys_model_detector.get();
      
      std::optional<RelActCalc::PhysModelShield<T>> self_atten;

      assert( !m_input.phys_model_self_atten || shield_is_present(m_input.phys_model_self_atten) );

      if( shield_is_present( m_input.phys_model_self_atten ) )
      {
        RelActCalc::PhysModelShield<T> att;
        att.material = m_input.phys_model_self_atten->material;
        
        if( att.material )
        {
          att.atomic_number = T(0.0);
        }else if( m_input.phys_model_self_atten->fit_atomic_number )
        {
          T an = x[par_num] * RelActCalc::ns_an_ceres_mult;
          // Do NOT clamp `an` with fmin/fmax: that zeroes the ceres::Jet derivative at the [1,98]
          //  bound, re-creating the bound-trap the areal-density AD path deliberately avoids (see
          //  RelActCalc_imp.hpp).  The AN parameter is already Ceres-bounded (to [lower_an,upper_an]/
          //  ns_an_ceres_mult), and get_atten_coef_for_an interpolates smoothly - keeping the gradient
          //  non-zero - for the tiny margin the optimizer can reach past the bound.
          att.atomic_number = an;
          par_num += 1;
        }else
        {
          att.atomic_number = T( fmin( fmax(m_input.phys_model_self_atten->atomic_number, 1.0), 98.0) );
        }
        
        att.areal_density = x[par_num] * PhysicalUnits::g_per_cm2;
        par_num += 1;
        
        self_atten = std::move(att);
      }//if( use internal attenuation shielding )
      
      vector<RelActCalc::PhysModelShield<T>> external_attens;
      for( const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &ext_atten : m_input.phys_model_external_attens )
      {
        assert( ext_atten );
        if( !ext_atten )
          continue;
        
        assert( shield_is_present(ext_atten) );

        if( !shield_is_present(ext_atten) )
        {
          assert( 0 );
          continue;
        }

        RelActCalc::PhysModelShield<T> att;
        att.material = ext_atten->material;

        if( att.material )
        {
          att.atomic_number = T(0.0);
        }else if( ext_atten->fit_atomic_number )
        {
          T an = x[par_num] * RelActCalc::ns_an_ceres_mult;
          // Do NOT clamp `an` with fmin/fmax: that zeroes the ceres::Jet derivative at the [1,98]
          //  bound, re-creating the bound-trap the areal-density AD path deliberately avoids (see
          //  RelActCalc_imp.hpp).  The AN parameter is already Ceres-bounded (to [lower_an,upper_an]/
          //  ns_an_ceres_mult), and get_atten_coef_for_an interpolates smoothly - keeping the gradient
          //  non-zero - for the tiny margin the optimizer can reach past the bound.
          att.atomic_number = an;
          par_num += 1;
        }else
        {
          assert( (ext_atten->atomic_number >= 1.0) && (ext_atten->atomic_number <= 98.0) );
          att.atomic_number = T( ext_atten->atomic_number );
        }
        
        att.areal_density = x[par_num] * PhysicalUnits::g_per_cm2;
        par_num += 1;
        
        external_attens.push_back( std::move(att) );
      }//for( loop over external attenuators )
      
      
      std::optional<T> b, c;
      if( m_input.phys_model_use_hoerl )
      {
        b = (x[par_num] - RelActCalc::ns_decay_hoerl_b_offset) * RelActCalc::ns_decay_hoerl_b_multiple;
        par_num += 1;
        c = (x[par_num] - RelActCalc::ns_decay_hoerl_c_offset) * RelActCalc::ns_decay_hoerl_c_multiple;
        par_num += 1;
      }
      
      rel_eff_fcn = [self_atten,external_attens, det, b, c]( double energy ) -> T {
        return RelActCalc::eval_physical_model_eqn_imp<T>( energy, self_atten, external_attens, det, b, c );
      };
      
      assert( par_num == x.size() );
    }else
    {
      assert( x.size() == (num_isos + m_input.eqn_order + 1) );
      rel_eff_fcn = [&x, this, num_isos]( double energy ){
        return eval_eqn_imp( energy, m_input.eqn_form, &(x[num_isos]), m_input.eqn_order + 1 );
      };
    }
    
    for( size_t index = 0; index < m_input.peaks.size(); ++index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = m_input.peaks[index];
      
      T curve_val = rel_eff_fcn( peak.m_energy );
      
      T rel_src_counts{}; //scalar part of Jet will be default constructed to 0.0
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        const T rel_activity = relative_activity( line.m_isotope, x );
        rel_src_counts += rel_activity * line.m_yield;
      }//for( const GenericLineInfo &line : peak.m_source_gammas )
      
      
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
      if( (index == 0) && (m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel) )
      {
        // Pin the first relative efficiency to 1.0; this is to prevent degeneracy of LM and LLS
        //  both essentially fitting for the overall normalization
        const double weight = 10.0; //Fairly arbitrary - not sure what to actually use - need something so the uncerts on activities dont get huge (which they will with no weight)
        residuals[m_peak_infos.size()] = weight*((peak.m_counts / rel_src_counts) - 1.0);
      }
#endif
      
      
      // TODO: should we fold in the uncertainty from the relative efficiency equation into the below as well - I dont think so since we are treating them as nuisance parameters?
      
      // Note: we want to compute:
      //  `(rel_efficiency - curve_val) / rel_eff_uncertainty`,
      //  But this has some divide by zero issues when rel act is zero, but the below gives the same
      //  residual, but avoiding these issues.
      if( peak.m_base_rel_eff_uncert == -1.0 )
      {
        // We are doing an unweighted fit
        
        // Avoid dividing by zero, so make sure rel_src_counts isnt really close to zero.
        // TODO: - need to evaluate using a different formulation where \c rel_src_counts being zero isnt a problem (I think its fine, but need to check before making the change - and make sure it wont effect how we use the covariances).
        if( ((rel_src_counts < 1.0E-8) && (rel_src_counts < (1.0E-6*peak.m_counts)))
           || (rel_src_counts < numeric_limits<T>::epsilon()) )
        {
          rel_src_counts = T(1.0E-6) * peak.m_counts;
        }
        
        residuals[index] = (peak.m_counts / rel_src_counts) - curve_val;
      }else if( peak.m_base_rel_eff_uncert == 0.0 )
      {
        // We are not using a m_base_rel_eff_uncert value
        const T pred_counts = curve_val * rel_src_counts;
        residuals[index] = (peak.m_counts - pred_counts) / peak.m_counts_uncert;
      }else
      {
        // We are using a m_base_rel_eff_uncert value
        assert( peak.m_base_rel_eff_uncert <= 1.0 );
        
        const T pred_counts = curve_val * rel_src_counts;
        // Note: for `add_uncert` below, we are using peak.m_counts, but it *could* also be
        //       reasonable to use `rel_src_counts` (which I was doing pre 20220720) or
        //       even `pred_counts`.  This is maybe worth revisiting.  Note that if you change
        //       how things are calculated here, you should also be consistent with
        //       #fit_rel_eff_eqn_lls.
        const double add_uncert = peak.m_counts * peak.m_base_rel_eff_uncert;
        const double uncert = sqrt( pow(peak.m_counts_uncert,2.0) + pow(add_uncert,2.0) );

        residuals[index] = (peak.m_counts - pred_counts) / uncert;
      }
    }//for( loop over energies to evaluate at )
  }//eval_internal_nl_rel_eff(...)

/*
  template<typename T, int N>
  void eval( std::vector<ceres::Jet<T,N>> x, ceres::Jet<T,N> *residuals ) const
  {
    
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      const size_t num_isos = num_isotopes();
      // TODO: make sure scalar portion of Jet is valid range - need to use proper logic for parameter location...
      //We'll give the minimizer a little wiggle room to go negative for AD, but not to much.
      //  However, we'll make sure the AD is non-negative.
      //assert( x[num_isos + 1] > -1.0E-6 );
      //x[num_isos + 1].a = std::max( 0.0, x[num_isos + 1].a );
      //for( size_t i = 0; i < m_input.phys_model_external_attens.size(); ++i )
      //{
      //  assert( x[num_isos + 2 + 2*i + 1] > -1.0E-6 );
      //  x[num_isos + 2 + 2*i + 1].a = std::max( 0.0, x[num_isos + 2 + 2*i + 1].a );
      //}
    }//if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    
    if( !m_input.use_ceres_to_fit_eqn )
      throw std::logic_error( "ManualGenericRelActFunctor: relative efficiency equation should be fit using Ceres when using auto-differentiation." );

    eval_internal_nl_rel_eff( x, residuals );
  }
  */

  template<typename T>
  void eval( std::vector<T> x, T *residuals ) const
  {
    if( m_input.use_ceres_to_fit_eqn )
      eval_internal_nl_rel_eff<T>( x, residuals );
    else
      eval_internal_lls_rel_eff<T>( x, residuals );

    // The profile-likelihood penalty channel (last residual); zero unless armed, so the nominal
    //  fit is bit-for-bit what it would be without profiling.
    if( !m_input.profile_targets.empty() )
    {
      const size_t profile_index = number_residuals() - 1;

      if( m_profile_active < 0 )
      {
        residuals[profile_index] = T(0.0);
      }else
      {
        assert( static_cast<size_t>(m_profile_active) < m_input.profile_targets.size() );
        const RelActCalcManual::ProfileTarget &target = m_input.profile_targets[m_profile_active];
        residuals[profile_index] = m_profile_weight
                    * (profile_quantity( target, x, m_profile_denom_floor ) - m_profile_target);
      }
    }//if( the penalty channel exists )


#ifndef NDEBUG
    for( size_t index = 0; index < number_residuals(); ++index )
    {
      assert( !isnan(residuals[index]) && !isinf(residuals[index]) );
      if( isnan(residuals[index]) || isinf(residuals[index]) )
        cerr << "Residual " << index << " has inf or nan value." << endl;

      if constexpr ( !std::is_same_v<T, double> )
      {
        for( size_t jac_index = 0; jac_index < residuals[index].v.rows(); ++jac_index )
        {
          assert( !IsNan(residuals[index].v[jac_index]) && !IsInf(residuals[index].v[jac_index]) );
          if( isnan(residuals[index].v[jac_index]) || isinf(residuals[index].v[jac_index]) )
            cerr << "Residual " << index << " has Jacobian " << jac_index << " value if inf or nan" << endl;
        }
        //cout << "Residual " << index << " has value " << residuals[index].a << endl;
      }else
      {
        //cout << "Residual " << index << " has value " << residuals[index] << endl;
      }//if constexpr ( !std::is_same_v<T, double> ) / else
    }//for( size_t index = 0; index < number_residuals(); ++index )
#endif
  }
  

  /*
  virtual double operator()( const std::vector<double> &x ) const
  {
    vector<double> residuals( number_residuals(), 0.0 );
    try
    {
      eval( x, residuals.data() );
    }catch( std::exception &e )
    {
      cerr << "ManualGenericRelActFunctor::operator() caught: " << e.what() << endl;
      return std::numeric_limits<double>::max();
    }
    
    double chi2 = 0.0;
    for( size_t i = 0; i < m_input.peaks.size(); ++i )
      chi2 += residuals[i]*residuals[i];
    
    return chi2;
  }//operator() - for minuit
  
  
  // For Minuit2
  virtual double Up() const
  {
    return 1.0;
  }
   */

  size_t num_isotopes() const
  {
    return m_isotopes.size();
  }

  size_t num_parameters() const
  {
    size_t num_pars = m_isotopes.size();
    if( !m_input.use_ceres_to_fit_eqn )
      return num_pars;
    
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      if( shield_is_present( m_input.phys_model_self_atten ) )
      {
        if( shield_fits_an( m_input.phys_model_self_atten ) )
          num_pars += 1; //AN is only a parameter if being fit
        num_pars += 1; //AD is always a paramater
      }

      for( const auto &atten : m_input.phys_model_external_attens )
      {
        assert( atten );
        if( shield_is_present( atten ) )
        {
          if( shield_fits_an( atten ) )
            num_pars += 1;
          num_pars += 1;
        }
      }
      
      if( m_input.phys_model_use_hoerl )
        num_pars += 2;
    }else
    {
      num_pars += m_input.eqn_order + 1;
    }
    
    return num_pars;
  }
  
  // Walk to controlling nuclide similar to solution.walk_to_controlling_nuclide()
  bool walk_to_controlling_nuclide( size_t &iso_index, double &multiple ) const
  {
    assert( multiple == 1.0 );
    multiple = 1.0;
    assert( iso_index < m_isotopes.size() );
    
    if( m_input.act_ratio_constraints.empty() )
      return false;
    
    size_t controller_index = std::numeric_limits<size_t>::max();
    bool found_controller = false;
    size_t sentinel = 0; // Safety counter to prevent infinite loops
    
    while( controller_index != iso_index )
    {
      sentinel += 1;
      assert( sentinel < 100 );
      if( sentinel > 1000 )
        throw std::logic_error( "ManualGenericRelActFunctor::walk_to_controlling_nuclide: possible infinite loop - logic error" );
      
      controller_index = iso_index;
      for( const RelActCalcManual::ManualActRatioConstraint &constraint : m_input.act_ratio_constraints )
      {
        if( constraint.m_constrained_nuclide == m_isotopes[iso_index] )
        {
          // Find the index of the controlling nuclide
          iso_index = this->iso_index( constraint.m_controlling_nuclide );
          multiple *= constraint.m_constrained_to_controlled_activity_ratio;
          found_controller = true;
          break;
        }
      }
    }//while( controller_index != iso_index )
    
    return found_controller;
  }//bool walk_to_controlling_nuclide( size_t &iso_index, double &multiple ) const
  
  // The return value indicates whether the computation of the
  // residuals and/or jacobians was successful or not.
  template<typename T>
  bool operator()(T const* const* parameters, T* residuals) const 
  {
    try
    {
      vector<T> pars( parameters[0], parameters[0] + num_parameters() );
     
      eval( pars, residuals );
    }catch( std::exception &e )
    {
      cerr << "ManualGenericRelActFunctor::operator() caught: " << e.what() << endl;
      return false;
    }
    
    return true;
  };//bool operator() - for Ceres

};//class ManualGenericRelActFunctor
}//namespace




namespace RelActCalcManual
{
GenericLineInfo::GenericLineInfo()
: m_yield( std::numeric_limits<double>::quiet_NaN() ),
  m_isotope( "InvalidIsotope" )
{
}

GenericLineInfo::GenericLineInfo( const double yield, const std::string &isotope )
 : m_yield( yield ),
   m_isotope( std::move(isotope) )
{
}


GenericPeakInfo::GenericPeakInfo()
 : m_energy( 0.0 ),
   m_mean( 0.0 ),
   m_fwhm( 0.0 ),
   m_counts( 0.0 ),
   m_counts_uncert( 0.0 ),
   m_base_rel_eff_uncert( 0.0 ),
   m_source_gammas{}
{
}

/*
void ManualActRatioConstraint::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  if( m_constrained_to_controlled_activity_ratio <= 0.0 )
    throw logic_error( "ManualActRatioConstraint::toXml: Constrained to controlled activity ratio is less than or equal to 0.0." );

  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "ManualActRatioConstraint::toXml: invalid parent." );

  rapidxml::xml_document<char> *doc = parent->document();
  rapidxml::xml_node<char> *base_node = doc->allocate_node( node_element, "ManualActRatioConstraint", nullptr, 24, 0 );
  parent->append_node( base_node );
  append_version_attrib( base_node, ManualActRatioConstraint::sm_xmlSerializationVersion );
  
  append_string_node( base_node, "ControllingNuclide", m_controlling_nuclide );
  append_string_node( base_node, "ConstrainedNuclide", m_constrained_nuclide );
  append_float_node( base_node, "ActivityRatio", m_constrained_to_controlled_activity_ratio );
}//ManualActRatioConstraint::toXml(...)

void ManualActRatioConstraint::fromXml( const ::rapidxml::xml_node<char> *constraint_node )
{
  if( !constraint_node )
    throw runtime_error( "ManualActRatioConstraint::fromXml: invalid input" );
    
  if( !rapidxml::internal::compare( constraint_node->name(), constraint_node->name_size(), "ManualActRatioConstraint", 24, false ) )
    throw std::logic_error( "invalid input node name" );
    
  // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
  static_assert( ManualActRatioConstraint::sm_xmlSerializationVersion == 0,
                  "ManualActRatioConstraint::fromXml: needs to be updated for new serialization version." );
    
  XmlUtils::check_xml_version( constraint_node, ManualActRatioConstraint::sm_xmlSerializationVersion );

  const rapidxml::xml_node<char> *controlling_node = XmlUtils::get_required_node( constraint_node, "ControllingNuclide" );
  const rapidxml::xml_node<char> *constrained_node = XmlUtils::get_required_node( constraint_node, "ConstrainedNuclide" );
  
  m_controlling_nuclide = SpecUtils::xml_value_str( controlling_node );
  m_constrained_nuclide = SpecUtils::xml_value_str( constrained_node );
  m_constrained_to_controlled_activity_ratio = XmlUtils::get_float_node_value( constraint_node, "ActivityRatio" );
  if( m_constrained_to_controlled_activity_ratio <= 0.0 )
    throw runtime_error( "ManualActRatioConstraint::fromXml: Activity ratio is less than or equal to 0.0." );
}//ManualActRatioConstraint::fromXml(...)

#if( PERFORM_DEVELOPER_CHECKS )
void ManualActRatioConstraint::equalEnough( const ManualActRatioConstraint &lhs, const ManualActRatioConstraint &rhs )
{
  if( fabs(lhs.m_constrained_to_controlled_activity_ratio - rhs.m_constrained_to_controlled_activity_ratio) > 1e-6 )
    throw logic_error( "ManualActRatioConstraint: Constrained to controlled activity ratio is not equal." );

  if( lhs.m_controlling_nuclide != rhs.m_controlling_nuclide )
    throw logic_error( "ManualActRatioConstraint: Controlling nuclide is not equal." );

  if( lhs.m_constrained_nuclide != rhs.m_constrained_nuclide )
    throw logic_error( "ManualActRatioConstraint: Constrained nuclide is not equal." );
}
#endif
*/


void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form, const size_t order,
                             const vector<double> &energies,
                             const vector<double> &data_values,
                             const vector<double> &data_uncertainties_orig,
                             vector<double> &fit_pars,
                             vector<vector<double>> *covariance )
{
  fit_rel_eff_eqn_lls_imp( fcn_form, order, energies, data_values,
                            data_uncertainties_orig, fit_pars, covariance );
}
  
  
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<std::string> &isotopes,
                           const std::vector<double> &rel_acts,
                           const std::vector<GenericPeakInfo> &peak_infos,
                           std::vector<double> &fit_pars,
                           std::vector<std::vector<double>> *covariance )
{
  // The implementation looks isotopes up with std::lower_bound, so the (index-paired) isotope
  //  and activity lists must be sorted together; dont rely on the caller for that.
  assert( isotopes.size() == rel_acts.size() );
  if( isotopes.size() != rel_acts.size() )
    throw std::logic_error( "fit_rel_eff_eqn_lls: isotopes and rel_acts sizes dont match" );

  if( std::is_sorted( std::begin(isotopes), std::end(isotopes) ) )
  {
    fit_rel_eff_eqn_lls_imp( fcn_form, order, isotopes, rel_acts, peak_infos, fit_pars, covariance );
  }else
  {
    std::vector<size_t> order_index( isotopes.size() );
    for( size_t i = 0; i < order_index.size(); ++i )
      order_index[i] = i;
    std::sort( std::begin(order_index), std::end(order_index),
               [&isotopes]( const size_t lhs, const size_t rhs ){ return isotopes[lhs] < isotopes[rhs]; } );

    std::vector<std::string> sorted_isotopes( isotopes.size() );
    std::vector<double> sorted_rel_acts( rel_acts.size() );
    for( size_t i = 0; i < order_index.size(); ++i )
    {
      sorted_isotopes[i] = isotopes[order_index[i]];
      sorted_rel_acts[i] = rel_acts[order_index[i]];
    }

    fit_rel_eff_eqn_lls_imp( fcn_form, order, sorted_isotopes, sorted_rel_acts, peak_infos, fit_pars, covariance );
  }//if( already sorted ) / else
}
  
  
vector<GenericPeakInfo> add_nuclides_to_peaks( const std::vector<GenericPeakInfo> &peaks,
                                              const std::vector<SandiaDecayNuc> &nuclides,
                                              const double real_time,
                                              const double cluster_sigma )
{
  vector<GenericPeakInfo> answer = peaks;
  
  vector< pair<double,double> > energy_widths, energy_obs_counts, energy_obs_counts_uncert;
  for( const auto &p : peaks )
  {
    // The logic to assign nuclides below will fail if two input peaks have identical energies.
    const auto prev_pos = std::find_if(begin(energy_widths),end(energy_widths), [&p]( const pair<double,double> &pe ){
      return p.m_energy == pe.first;
    });
    assert( prev_pos == end(energy_widths) );
    if( prev_pos != end(energy_widths) )
      throw runtime_error( "add_nuclides_to_peaks: two peaks have exactly the same energies - not allowed." );
    
    energy_widths.push_back( {p.m_energy, p.m_fwhm / 2.35482} );
  }

  // `cluster_peak_activities(...)` looks peaks up with std::lower_bound, so the energy/width
  //  list must be sorted by energy - the caller's peak ordering is not otherwise required to be.
  //  (The per-peak yield write-back below is keyed by energy, so `answer` keeps caller order.)
  std::sort( begin(energy_widths), end(energy_widths),
             []( const pair<double,double> &lhs, const pair<double,double> &rhs ){
               return lhs.first < rhs.first;
             } );


  set<const void *> nuclides_seen;
  for( const auto &n : nuclides )
  {
    string name;
    const void *src_ptr = nullptr;
    if( n.nuclide )
    {
      src_ptr = static_cast<const void *>(n.nuclide);
      name = n.nuclide->symbol;
    }else if( n.element )
    {
      src_ptr = static_cast<const void *>(n.element);
      name = n.element->symbol;
    }else if( n.reaction )
    {
      src_ptr = static_cast<const void *>(n.reaction);
      name = n.reaction->name();
    }else
    {
      throw runtime_error( "add_nuclides_to_peaks: null input" );
    }
    
    if( nuclides_seen.count( src_ptr ) )
      throw runtime_error( "add_nuclides_to_peaks: input nuclides must be unique" );
    nuclides_seen.insert( src_ptr );
    
    // We will map from the peaks mean, to the total number of gammas that contribute to that peak,
    //  for this nuclide
    map<double,double> energy_gammas_map;
    
    if( n.nuclide )
    {
      SandiaDecay::NuclideMixture mixture;
      mixture.addNuclideByActivity(n.nuclide, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits);
      
      if( n.correct_for_decay_during_meas && (real_time <= 0) )
        throw runtime_error( "add_nuclides_to_peaks: measurement time must be specified if"
                            " correcting activities for nuclide decays during measurement.");
      
      const double activity = 1.0;
      const bool decay_correct = n.correct_for_decay_during_meas;
      
      GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map,
                                                                            energy_widths, mixture, activity, n.age, cluster_sigma, -1,
                                                                            decay_correct, real_time, nullptr, nullptr );
    }else if( n.element )
    {
      for( const pair<double,double> &p : energy_widths )
      {
        double src_yield = 0.0;
        const double lower_x = p.first - cluster_sigma*p.second;
        const double upper_x = p.first + cluster_sigma*p.second;
        for( const SandiaDecay::EnergyIntensityPair &xray : n.element->xrays )
        {
          if( (xray.energy >= lower_x) && (xray.energy <= upper_x) )
            src_yield += xray.intensity;
        }
        
        energy_gammas_map[p.first] = src_yield;
      }//for( const auto &p : peaks )
    }else if( n.reaction )
    {
      for( const pair<double,double> &p : energy_widths )
      {
        double src_yield = 0.0;
        const double lower_x = p.first - cluster_sigma*p.second;
        const double upper_x = p.first + cluster_sigma*p.second;
        for( const ReactionGamma::Reaction::EnergyYield &rxctn : n.reaction->gammas )
        {
          if( (rxctn.energy >= lower_x) && (rxctn.energy <= upper_x) )
            src_yield += rxctn.abundance;
        }
        
        energy_gammas_map[p.first] = src_yield;
      }//for( const auto &p : peaks )
    }
    
    // Write the per-peak yields back keyed by ENERGY: `energy_gammas_map` iterates in ascending
    //  energy order, while `answer` preserves the caller's peak ordering, so a positional zip
    //  would silently mis-assign yields whenever the input peaks are not energy-sorted.
    //  (Duplicate energies were rejected above, so the lookup is unambiguous; keys were inserted
    //  as the exact same peak-energy doubles, so exact find() is safe.)
    assert( energy_gammas_map.size() == answer.size() );

    for( GenericPeakInfo &peak : answer )
    {
      const map<double,double>::const_iterator pos = energy_gammas_map.find( peak.m_energy );
      assert( pos != end(energy_gammas_map) );
      const double yield = (pos == end(energy_gammas_map)) ? 0.0 : pos->second;
      if( yield > numeric_limits<float>::min() )
        peak.m_source_gammas.emplace_back( yield, name );
    }//for( GenericPeakInfo &peak : answer )
  }//for( const auto &n : nuclides )
  
  // an alternate way to do this, but they dont match exactly
  //vector<RelActCalcManual::PeakCsvInput::NucAndAge> isotopes;
        //for( const auto &n : nuc_sources )
        //{
        //  if( n.nuclide )
        //    isotopes.emplace_back( n.nuclide->symbol, n.age, correct_for_decay );
        //}
        //RelActCalcManual::PeakCsvInput::NucMatchResults matched_res
        //  = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( peaks_in_range,
        //                                    RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay,
        //                                    {}, isotopes, cluster_num_sigma, {}, real_time );
        // const auto peaks_with_nucs = matched_res.peaks_matched;
      
        /*
        // Code to help debug difference between matching stuff...
        for( auto &p : peaks_with_nucs )
        {
          std::sort( begin(p.m_source_gammas), end(p.m_source_gammas), []( auto &lhs, auto &rhs ){
            return lhs.m_isotope < rhs.m_isotope;
          } );
        }
      
        for( auto &p : matched_res.peaks_matched )
        {
          std::sort( begin(p.m_source_gammas), end(p.m_source_gammas), []( auto &lhs, auto &rhs ){
            return lhs.m_isotope < rhs.m_isotope;
          } );
        }
      
        assert( matched_res.peaks_matched.size() == peaks_with_nucs.size() );
      
      
        for( size_t i = 0; i < std::max(matched_res.peaks_matched.size(), peaks_with_nucs.size()); ++i )
        {
          const auto newp = matched_res.peaks_matched[i];
          const auto oldp = peaks_with_nucs[i];
          assert( newp.m_energy == oldp.m_energy );
          assert( newp.m_counts == oldp.m_counts );
          assert( newp.m_counts_uncert == oldp.m_counts_uncert );
          assert( newp.m_fwhm == oldp.m_fwhm );
          assert( newp.m_base_rel_eff_uncert == oldp.m_base_rel_eff_uncert );
          assert( newp.m_source_gammas.size() == oldp.m_source_gammas.size() );
          for( size_t j = 0; j < newp.m_source_gammas.size(); ++j )
          {
            assert( newp.m_source_gammas[j].m_isotope == oldp.m_source_gammas[j].m_isotope );
          
            double diff = fabs( newp.m_source_gammas[j].m_yield - oldp.m_source_gammas[j].m_yield );
            assert( diff <= 0.00001*newp.m_source_gammas[j].m_yield );
            assert( diff <= 0.00001*oldp.m_source_gammas[j].m_yield );
            if( newp.m_source_gammas[j].m_yield != oldp.m_source_gammas[j].m_yield )
            {
              double brnew = newp.m_source_gammas[j].m_yield;
              double brold = oldp.m_source_gammas[j].m_yield;
              cout << "Mismatcht BR: " << brnew << " vs " << brold << " for " << newp.m_energy << " keV" << endl;
              cout << endl;
            }
            //assert( newp.m_source_gammas[j].m_yield == oldp.m_source_gammas[j].m_yield );
          }
        
          if( i < matched_res.peaks_matched.size() )
          {
            const auto p = matched_res.peaks_matched[i];
            cout << "new " << i << ": e=" << p.m_energy << ", fwhm=" << p.m_fwhm << endl;
            for( const auto g : p.m_source_gammas )
              cout << "\tsource: " << g.m_isotope << ": " << g.m_yield << endl;
          }
        
          if( i < peaks_with_nucs.size() )
          {
            const auto p = peaks_with_nucs[i];
            cout << "old " << i << ": e=" << p.m_energy << ", fwhm=" << p.m_fwhm << endl;
            for( const auto g : p.m_source_gammas )
              cout << "\tsource: " << g.m_isotope << ": " << g.m_yield << endl;
          }
        }
        cout << endl << endl;
        //peaks_with_nucs = matched_res.peaks_matched;
      */

  return answer;
}//add_nuclides_to_peaks(...)


/*
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                         const size_t order,
                         const std::vector<SandiaDecayNucRelAct<double>> &nuclides,
                         const double base_rel_eff_uncert,
                         const std::vector<std::shared_ptr<const PeakDef>> &peak_infos,
                         vector<double> &fit_pars,
                         std::vector<std::vector<double>> *covariance )
{
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  assert( !nuclides.empty() );
  if( nuclides.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no nuclides specified." );
  
  assert( !peak_infos.empty() );
  if( peak_infos.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no peaks specified." );
  
  // We will map from the peaks mean, to the total number of gammas that contribute to that peak
  map<double,double> energy_gammas_map;
  
  // Get peak energies and widths (normally width this is just 'sigma', but for non-Gaussian peaks
  //  its 0.25 of the ROI)
  set<double> energies_seen; //a very poor check that there arent duplicate peaks
  vector< pair<double,double> > energy_widths, energy_obs_counts, energy_obs_counts_uncert;
  for( const auto &p : peak_infos )
  {
    // To use non-Gaussian peaks we would need to pass in the shared_ptr<const Measurement> data...
    //  maybe later if it ever matters
    if( !p->gausPeak() )
      throw runtime_error( "fit_rel_eff_eqn_lls: non-Gaussian peaks not supported yet" );
    
    // Note that in GammaInteractionCalc::ShieldingSourceChi2Fcn::observedPeakEnergyWidths
    //  we use the assigned nuclides gamma energy, as the energy - here we are using the peak mean.
    //  TODO: - revisit ether to use peak mean or its nuclide gamma as the energy - after implementing the rest of the manual RelAct calc stuff.
    const double energy = p->mean();
    const double sigma = p->gausPeak() ? p->sigma() : 0.25*p->roiWidth();
    const double amp = p->amplitude();
    //const double amp = p->gausPeak() ? p->amplitude() : p->areaFromData(data);
    const double ampUncert = p->amplitudeUncert();

    if( energies_seen.count(energy) )
      throw runtime_error( "fit_rel_eff_eqn_lls: multiple peaks with same energy - not allowed." );
    energies_seen.insert( energy );
    
    energy_widths.push_back( {energy, sigma} );
    energy_obs_counts.push_back( {energy, amp} );
    energy_obs_counts_uncert.push_back( {energy, ampUncert} );
  }//for( const PeakDef &peak : peaks )
  
  // JIC the peaks werent sorted, sort by just energies (although we did check no duplicate energies
  //  but we'll play it safe)
  auto sortByFirstOnly = []( const pair<double, double> &lhs, const pair<double, double> &rhs ){
    return lhs.first < rhs.first;
  };
  
  std::stable_sort( begin(energy_widths), end(energy_widths), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts), end(energy_obs_counts), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts_uncert), end(energy_obs_counts_uncert), sortByFirstOnly );
  
  
  // Now we will go through and get the amplitude of gammas we expect to contribute to a single peak
  //  (there may be multiple gammas from the same nuclide, as well as multiple nuclides that
  //   contribute to a single observable peaks).
  //  We will select a 'cluster' sigma of 1.5; this is what Activity/Shielding fit uses, but I dont
  //  think this value was derived by anything more than "that seems about right", and I havent
  //  run into an obvious case where this is not correct.
  const double photopeakClusterSigma = 1.5;
  set<const SandiaDecay::Nuclide *> nuclides_seen;
  for( const auto &n : nuclides )
  {
    if( nuclides_seen.count(n.nuclide) )
      throw runtime_error( "fit_rel_eff_eqn_lls: input nuclides must be unique" );
  
    nuclides_seen.insert( n.nuclide );
    
    SandiaDecay::NuclideMixture mixture;
    mixture.addNuclideByActivity(n.nuclide, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits);
    
    const double energyToCluster = -1;
    // TODO: we could account for decays during the measurement, but would need realTime here
    const bool accountForDecayDuringMeas = false;
    const double realTime = -1;
    GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map,
                                            energy_widths, mixture, n.rel_activity, n.age,
                                            photopeakClusterSigma, energyToCluster,
                                            accountForDecayDuringMeas, realTime, nullptr, nullptr );
  }//for( const auto &n : nuclides )
  
  // Convert energy_gammas_map to a vector for convenience
  vector<pair<double,double>> energy_gammas;
  for( const auto &ec : energy_gammas_map )
    energy_gammas.push_back( ec );

  assert( energy_gammas.size() == peak_infos.size() );
  assert( energy_gammas.size() == energy_widths.size() );
  assert( energy_gammas.size() == energy_obs_counts.size() );
  assert( energy_gammas.size() == energy_obs_counts_uncert.size() );
  
  
  // Now put all this info onto a form so we can call into fit_rel_eff_eqn_lls(...); there is a
  //  commented out implementation of not having to do this, yet another, transformation of
  //  information.
  
  double max_pred_counts = 0.0;
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
    max_pred_counts = std::max(max_pred_counts, energy_gammas[peak_index].second);
  
  // Instead of keeping counts from each nuclide for each peak separate, we summed all nuclides
  //  together for each peak - so here we'll only use a single "Effective" isotope.
  const vector<string> isotopes{ "Effective" };
  const vector<double> rel_acts( 1, max_pred_counts );
  vector<GenericPeakInfo> generic_peak_infos;
  
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
  {
    GenericPeakInfo peak;
    peak.m_mean = peak_infos[peak_index]->mean();
    peak.m_energy = energy_gammas[peak_index].first;
    peak.m_fwhm = 2.35482*energy_widths[peak_index].second;
    peak.m_counts = energy_obs_counts[peak_index].second;
    peak.m_counts_uncert = (energy_obs_counts_uncert[peak_index].second > 0.0)
                             ? energy_obs_counts_uncert[peak_index].second
                             : sqrt(peak.m_counts);
    peak.m_base_rel_eff_uncert = base_rel_eff_uncert;
    
    const double yield = energy_gammas[peak_index].second / max_pred_counts;
    peak.m_source_gammas.emplace_back( yield, "Generic" );
    
    generic_peak_infos.push_back( peak );
  }//for( size_t row = 0; row < num_peaks; ++row )
  
  return fit_rel_eff_eqn_lls( fcn_form, order, isotopes, rel_acts, generic_peak_infos,
                              fit_pars, covariance );
}//fit_rel_eff_eqn_lls(...)
*/

void fit_act_to_rel_eff( const RelActCalc::RelEffEqnForm eqn_form,
                        const std::vector<double> &eqn_pars,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_uncerts )
{
  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    throw runtime_error( "fit_act_to_rel_eff: FramPhysicalModel not supported." );
  
  if( eqn_pars.empty() || eqn_pars.size() > 10 )
    throw runtime_error( "fit_act_to_rel_eff: invalid equation passed in." );

  auto eff_eqn = [eqn_form, eqn_pars]( double energy ){
    return eval_eqn( energy, eqn_form, eqn_pars );
  };
  
  fit_act_to_rel_eff( eff_eqn, isotopes, peak_infos, fit_rel_acts, fit_uncerts );
}//fit_act_to_rel_eff(...)


void fit_act_to_phys_rel_eff( const RelEffInput &input,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_uncerts )
{
  if( !input.phys_model_detector || !input.phys_model_detector->isValid() )
    throw runtime_error( "fit_act_to_phys_rel_eff: invalid detector passed in." );
  
  
  vector<int> dummy_const_pars;
  vector<double> rel_eff_pars( 2 + 2*input.phys_model_external_attens.size() + 2 ); //Likely too big - we'll resize later
  size_t rel_eff_index = 0;
  setup_physical_model_shield_par_manual( dummy_const_pars, rel_eff_pars.data(),
                                         rel_eff_index, input.phys_model_self_atten );
  
  for( const auto &opt : input.phys_model_external_attens )
  {
    assert( rel_eff_index < rel_eff_pars.size() );
    setup_physical_model_shield_par_manual( dummy_const_pars, rel_eff_pars.data(),
                                           rel_eff_index, opt );
  }
  
  if( input.phys_model_use_hoerl )
  {
    assert( rel_eff_index < rel_eff_pars.size() );
    rel_eff_pars[rel_eff_index] = (0.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset; //(energy/1000)^b
    rel_eff_index += 1;
    assert( rel_eff_index < rel_eff_pars.size() );
    rel_eff_pars[rel_eff_index] = (1.0/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset; //c^(1000/energy)
    rel_eff_index += 1;
  }//if( input.phys_model_use_hoerl )
  
  assert( rel_eff_index <= rel_eff_pars.size() );
  
  rel_eff_pars.resize( rel_eff_index );

  const ManualGenericRelActFunctor::PhysModelRelEqnDef eqn_input
        = ManualGenericRelActFunctor::make_phys_eqn_input( input, rel_eff_pars );
  
  const function<double(double)> eff_eqn
         = RelActCalc::physical_model_eff_function( eqn_input.self_atten, eqn_input.external_attens,
                                        eqn_input.det, eqn_input.hoerl_b, eqn_input.hoerl_c );

  fit_act_to_rel_eff( eff_eqn, isotopes, peak_infos, fit_rel_acts, fit_uncerts );
}//fit_act_to_rel_eff(...)


void fit_act_to_rel_eff( const std::function<double(double)> &eff_fcn,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_uncerts )
{
  // We want to fit the relative activities of the isotopes, given a relative efficiency curve
  
  // We'll start off with some basic sanity checks of the input.
  if( !eff_fcn )
    throw runtime_error( "fit_act_to_rel_eff: invalid rel eff fcn passed in." );
  
  if( peak_infos.size() < isotopes.size() )
    throw runtime_error( "fit_act_to_rel_eff: less peaks than isotopes." );
  
  std::vector<bool> used_isotopes( isotopes.size(), false );
  for( const GenericPeakInfo &info : peak_infos )
  {
    // We could blissfully ignore peaks with no source gammas (and assign them an activity of zero),
    //  but we'll be strict to maybe prevent some input mistakes, or something.
    if( info.m_source_gammas.empty() )
      throw runtime_error( "fit_act_to_rel_eff: peak at "
                          + to_string(info.m_mean) + " keV has no source gammas." );
    
    for( const GenericLineInfo &line : info.m_source_gammas )
    {
      if( line.m_yield <= 0.0 || IsInf(line.m_yield) || IsNan(line.m_yield) )
        throw runtime_error( "fit_act_to_rel_eff: invalid yield." );
      
      const auto pos = find( begin(isotopes), end(isotopes), line.m_isotope );
      if( pos == end(isotopes) )
        throw runtime_error( "fit_act_to_rel_eff: peak source isotope '" + line.m_isotope
                            + "' is not in list of wanted isotopes." );
      
      used_isotopes[pos - begin(isotopes)] = true;
    }//for( const GenericLineInfo &line : info.m_source_gammas )
    
    // Check peaks dont have the same nuclide multiple times
    for( size_t i = 1; i < info.m_source_gammas.size(); ++i )
    {
      for( size_t j = 0; j < i; ++j )
      {
        if( info.m_source_gammas[i].m_isotope == info.m_source_gammas[j].m_isotope )
          throw runtime_error( "fit_act_to_rel_eff: peak uses same isotope twice." );
      }//for( size_t j = 0; j < i; ++j )
    }//for( size_t i = 1; i < info.m_source_gammas.size(); ++i )
  }//for( const RelActCalc::GenericPeakInfo &info : peak_infos )
  
  assert( used_isotopes.size() == isotopes.size() );
  for( size_t i = 0; i < used_isotopes.size(); ++i )
  {
    if( !used_isotopes[i] )
      throw runtime_error( "fit_act_to_rel_eff: no peak with isotope '" + isotopes[i] + "' passed in." );
  }//for( size_t i = 0; i < used_isotopes.size(); ++i )
  
  // Checks are done, get to work
  
  //  We want to solve Ax = b, where
  //    Elements of A are the branching ratios for each isotope (e.g., where column correspond to
  //      isotopes), divided by the uncertainty of the peak - which corresponds to the row)
  //    x is the activities we are solving for
  //    b is the counts in each peak, divided by the relative efficiency for that energy (and all divided by the uncertainty of the peak).
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  
  const int num_isotopes = static_cast<int>( isotopes.size() );
  const int num_peaks = static_cast<int>( peak_infos.size() );
  
  Eigen::MatrixX<double> A = Eigen::MatrixX<double>::Zero( num_peaks, num_isotopes );
  Eigen::VectorX<double> b( num_peaks );
  
  
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const GenericPeakInfo &peak = peak_infos[row];
    const double energy = peak.m_energy;
    const double counts = peak.m_counts;
    const double counts_uncert = peak.m_counts_uncert;
    
    const double rel_eff = eff_fcn( energy );
    
    if( IsNan(rel_eff) || IsInf(rel_eff) )
      throw runtime_error( "fit_act_to_rel_eff: invalid rel eff at " + to_string(energy) + " keV." );

    if( counts_uncert <= 0.0 )
      throw runtime_error( "fit_act_to_rel_eff: peak counts uncertainty must be > 0 at "
                          + to_string(energy) + " keV." );

    // The per-peak weighted least-squares row is (sum_iso yield*rel_eff*act - counts)/counts_uncert.
    //  Writing b and A directly - instead of via rel_act = counts/rel_eff and
    //  rel_act_uncert = rel_act*counts_uncert/counts - avoids dividing by rel_eff and by counts; the
    //  algebra is identical (b = counts/counts_uncert, A = yield*rel_eff/counts_uncert), so the
    //  solution and covariance are unchanged, but a zero rel_eff or zero counts no longer produce NaN.
    b(row) = counts / counts_uncert;

    for( const GenericLineInfo &info : peak.m_source_gammas )
    {
      const auto pos = std::find( begin(isotopes), end(isotopes), info.m_isotope );
      assert( pos != end(isotopes) ); //we already checked this.

      const int column = static_cast<int>( pos - begin(isotopes) );
      assert( (column >= 0) && (column < static_cast<int>(isotopes.size())) ); //sometimes re-assurance is good

      assert( A(row,column) == 0.0 );
      A(row,column) = (info.m_yield * rel_eff) / counts_uncert;
    }//for( const GenericLineInfo &info : peak.m_source_gammas )
  }//for( size_t row = 0; row < num_peaks; ++row )
  
  // TODO: determine if HouseholderQr or BDC SVD is better/more-stable/faster/whatever
  //const Eigen::VectorXd solution = A.colPivHouseholderQr().solve(b);
  
  //const Eigen::BDCSVD<Eigen::MatrixX<double>> bdc = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV); //depreciated.
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
    const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> bdc(A); //slow, but best decomposition
#else
    const Eigen::BDCSVD<Eigen::MatrixX<double>> bdc(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif

  const Eigen::VectorXd solution = bdc.solve(b);
  
  assert( solution.size() == num_isotopes );
  
  // Covariance C = (A^T A)^{-1} via the SVD pseudo-inverse V*S^{-2}*V^T (see helper): avoids forming
  //  A^T A (which squares the condition number) and stays finite for rank-deficient A.
  const Eigen::MatrixX<double> C = lls_covariance_from_svd<double>( bdc.matrixV(), bdc.singularValues(),
                                                            static_cast<Eigen::Index>(num_peaks) );

  fit_rel_acts.resize( solution.size() );
  fit_uncerts.resize( solution.size() );

  for( size_t i = 0; i < num_isotopes; ++i )
  {
    fit_rel_acts[i] = solution(i);
    fit_uncerts[i] = std::sqrt( std::max( 0.0, C(i,i) ) );
  }
}//fit_act_to_rel_eff(...)


void RelEffInput::check_nuclide_constraints() const
{
  // Make sure nuclides in constraints are non-null, not the same nuclide, and are in the nuclides 
  //  list that we are fitting for.
  for( const ManualActRatioConstraint &nuc_constraint : act_ratio_constraints )
  {
    if( nuc_constraint.m_constrained_to_controlled_activity_ratio <= 0.0 )
      throw logic_error( "RelEffInput: Constrained to controlled activity ratio is less than or equal to 0.0." );

    if( nuc_constraint.m_constrained_nuclide.empty() )
      throw logic_error( "RelEffInput: Constrained nuclide is empty." );

    if( nuc_constraint.m_controlling_nuclide.empty() )
      throw logic_error( "RelEffInput: Controlling nuclide is empty." );

    if( nuc_constraint.m_constrained_nuclide == nuc_constraint.m_controlling_nuclide )
      throw logic_error( "RelEffInput: Constrained and controlling nuclides are the same." );

    // Check that the constrained nuclide is a nuclide in this RelEffCurve
    bool is_constrained_nuclide_in_curve = false, is_controlling_nuclide_in_curve = false;
    for( const GenericPeakInfo &peak : peaks )
    {
      for( const GenericLineInfo &line : peak.m_source_gammas )
      {
        if( nuc_constraint.m_constrained_nuclide == line.m_isotope )
          is_constrained_nuclide_in_curve = true;

        if( nuc_constraint.m_controlling_nuclide == line.m_isotope )
          is_controlling_nuclide_in_curve = true;
      }//for( const GenericLineInfo &line : peak.m_source_gammas )

      if( is_constrained_nuclide_in_curve && is_controlling_nuclide_in_curve )
        break;
    }//for( const GenericPeakInfo &peak : peaks )
        
    if( !is_constrained_nuclide_in_curve )
      throw logic_error( "RelEffInput: Constrained nuclide is not in any peak." );

    if( !is_controlling_nuclide_in_curve )
      throw logic_error( "RelEffInput: Controlling nuclide is not in any peak." );


#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
    //Check `nuc_constraint.m_constrained_nuclide` is not also in mass_fraction_constraints
    for( const RelActCalcManual::MassFractionConstraint &mass_frac_constraint : mass_fraction_constraints )
    {
      if( mass_frac_constraint.m_nuclide == nuc_constraint.m_constrained_nuclide )
        throw logic_error( "RelEffInput: Constrained nuclide is also in mass fraction constraints." );
    }
#endif // USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT
  }//for( const RelEffCurveInput::ActRatioConstraint &nuc_constraint : act_ratio_constraints )

  // Check that the constrained nuclide is only controlled by one nuclide
  for( const ManualActRatioConstraint &nuc_constraint : act_ratio_constraints )
  {
    size_t count = 0;
    for( const ManualActRatioConstraint &other_constraint : act_ratio_constraints )
    {
      if( other_constraint.m_constrained_nuclide == nuc_constraint.m_constrained_nuclide )
        ++count;
    }

    if( count > 1 )
      throw logic_error( "RelEffInput: Constrained nuclide is controlled by more than one nuclide." );
  }//for( const auto &nuc_constraint : act_ratio_constraints )

  // Make sure no duplicate constraints
  for( size_t i = 1; i < act_ratio_constraints.size(); ++i )
  {
    const ManualActRatioConstraint &outer_constraint = act_ratio_constraints[i];
    for( size_t j = 0; j < i; ++j )
    {
      const ManualActRatioConstraint &inner_constraint = act_ratio_constraints[j];
      if( (outer_constraint.m_constrained_nuclide == inner_constraint.m_constrained_nuclide)
        && (outer_constraint.m_controlling_nuclide == inner_constraint.m_controlling_nuclide) )
        {
          throw logic_error( "RelEffInput: Duplicate nuclide constraints." );
        }
    }
  }//for( size_t i = 0; i < act_ratio_constraints.size(); ++i )

  // Now we need to walk the chain of constraints to make sure we dont have a cycle
  // e.g. {constrained: U235, controlling: U238} -> {constrained: U238, controlling: U234} -> {constrained: U234, controlling: U235}
  for( size_t outer_index = 0; outer_index < act_ratio_constraints.size(); ++outer_index )
  { 
    const ManualActRatioConstraint &outer_constraint = act_ratio_constraints[outer_index];
    
    set<size_t> visited_constraints;
    visited_constraints.insert( outer_index );

    bool found_controller = true;
    const string *current_controller = &(outer_constraint.m_controlling_nuclide);  // e.g. U238

    while( found_controller )
    {
      found_controller = false;
      
      for( size_t inner_index = 0; inner_index < act_ratio_constraints.size(); ++inner_index )
      {
        const ManualActRatioConstraint &inner_constraint = act_ratio_constraints[inner_index];
        if( (*current_controller) == inner_constraint.m_constrained_nuclide )
        {
          if( visited_constraints.count( inner_index ) )
            throw logic_error( "Cycle in nuclide constraints." );

          found_controller = true;
          current_controller = &(inner_constraint.m_controlling_nuclide); // e.g. U234
          visited_constraints.insert( inner_index );
          break;
        }
      }//for( size_t inner_index = 0; inner_index < rel_eff_curve.act_ratio_constraints.size(); ++inner_index )
    }//while( found_constroller )
  }//for( size_t outer_index = 0; outer_index < act_ratio_constraints.size(); ++outer_index )

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
  // Guard against unbounded recursion between act-ratio chains and mass-fraction blocks:
  //  relative_activity() of an act-ratio constrained nuclide recurses to its (terminal)
  //  controller; if that controller is mass-fraction constrained, its sigma-block decode sums
  //  the element's non-mass-fraction-constrained isotopes by calling relative_activity() on
  //  each - and if such an isotope's own act-ratio chain terminates back on a nuclide of the
  //  same (or a mutually-referencing) block, the evaluation never terminates (a stack-overflow
  //  crash, not an exception).  Model each distinct element roster as a node, add an edge
  //  A -> B when decoding block A requires a chain that terminates on a nuclide mass-fraction
  //  constrained in block B, and reject any cycle.
  if( !act_ratio_constraints.empty() && !mass_fraction_constraints.empty() )
  {
    // Distinct element rosters ("blocks"); constraints sharing a nuclide are validated (below)
    //  to carry identical rosters, so keying on the sorted isotope set is well-defined.
    vector<set<string>> rosters;
    const auto roster_index_of_constraint = [&rosters]( const MassFractionConstraint &mfc ) -> size_t {
      set<string> roster;
      for( const auto &iso_act : mfc.m_specific_activities )
        roster.insert( iso_act.first );
      for( size_t i = 0; i < rosters.size(); ++i )
        if( rosters[i] == roster )
          return i;
      rosters.push_back( roster );
      return rosters.size() - 1;
    };

    const auto mass_frac_constraint_for = [this]( const string &iso ) -> const MassFractionConstraint * {
      for( const MassFractionConstraint &mfc : mass_fraction_constraints )
        if( mfc.m_nuclide == iso )
          return &mfc;
      return nullptr;
    };

    const auto terminal_of_chain = [this]( const string &iso ) -> const string * {
      // Walk act-ratio links until the current nuclide is not itself constrained; cycles among
      //  the act-ratio constraints were already rejected above, so this terminates.
      const string *current = &iso;
      bool moved = true;
      while( moved )
      {
        moved = false;
        for( const ManualActRatioConstraint &link : act_ratio_constraints )
        {
          if( link.m_constrained_nuclide == (*current) )
          {
            current = &(link.m_controlling_nuclide);
            moved = true;
            break;
          }
        }
      }//while( moved )
      return current;
    };

    vector<size_t> constraint_roster( mass_fraction_constraints.size() );
    for( size_t i = 0; i < mass_fraction_constraints.size(); ++i )
      constraint_roster[i] = roster_index_of_constraint( mass_fraction_constraints[i] );

    // Edges between blocks, plus the member/terminal names that created them (for the message).
    vector<set<size_t>> edges( rosters.size() );
    vector<map<size_t,pair<string,string>>> edge_reason( rosters.size() );
    for( size_t a = 0; a < rosters.size(); ++a )
    {
      for( const string &member : rosters[a] )
      {
        if( mass_frac_constraint_for(member) )
          continue; //the block decode does not recurse into mass-fraction constrained members

        const string * const terminal = terminal_of_chain( member );
        if( (*terminal) == member )
          continue; //not act-ratio constrained

        const MassFractionConstraint * const term_mfc = mass_frac_constraint_for( *terminal );
        if( !term_mfc )
          continue; //chain ends on a plain nuclide - evaluation terminates

        // Look the roster up WITHOUT inserting: growing `rosters` here would invalidate the
        //  range-for above (and leave `edges`, sized from `rosters`, too short).  Every
        //  constraint was already seeded into `rosters` by the loop above, so a miss is
        //  impossible; skip defensively rather than risk it.
        size_t b = rosters.size();
        {
          set<string> terminal_roster;
          for( const auto &iso_specact : term_mfc->m_specific_activities )
            terminal_roster.insert( iso_specact.first );
          for( size_t r = 0; r < rosters.size(); ++r )
            if( rosters[r] == terminal_roster )
              b = r;
        }
        assert( b < rosters.size() );
        if( b >= rosters.size() )
          continue;
        edges[a].insert( b );
        edge_reason[a].emplace( b, make_pair(member, *terminal) );
      }//for( const string &member : rosters[a] )
    }//for( size_t a = 0; a < rosters.size(); ++a )

    // Reject any cycle (including a self-loop) via a simple bounded walk from each node.
    for( size_t start = 0; start < rosters.size(); ++start )
    {
      set<size_t> visited{ start };
      vector<size_t> to_visit( begin(edges[start]), end(edges[start]) );
      while( !to_visit.empty() )
      {
        const size_t node = to_visit.back();
        to_visit.pop_back();
        if( node == start )
        {
          // The cycle may close through an intermediate block, so report the edge that leaves
          //  `start` (there is at least one, or we would not be walking).
          assert( !edge_reason[start].empty() );
          const pair<string,string> &reason = edge_reason[start].count(node)
                                                ? edge_reason[start][node]
                                                : begin(edge_reason[start])->second;
          throw logic_error( "RelEffInput: evaluating mass-fraction constrained element containing '"
              + reason.first + "' requires the activity of '" + reason.second
              + "' (through an activity-ratio constraint chain), which is mass-fraction"
              " constrained in a way that recurses back to the same element - this configuration"
              " is not supported." );
        }
        if( visited.count(node) )
          continue;
        visited.insert( node );
        to_visit.insert( end(to_visit), begin(edges[node]), end(edges[node]) );
      }//while( !to_visit.empty() )
    }//for( size_t start = 0; start < rosters.size(); ++start )
  }//if( have both act-ratio and mass-fraction constraints )
#endif // USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT


#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
  // Check that mass fraction constraints are valid
  for( const RelActCalcManual::MassFractionConstraint &constraint : mass_fraction_constraints )
  {
    if( constraint.m_nuclide.empty() )
      throw logic_error( "RelEffInput: Mass fraction constraint nuclide is empty." );

    if( constraint.m_mass_fraction_lower < 0.0 )
      throw logic_error( "RelEffInput: Mass fraction lower constraint is less than 0." );

    if( constraint.m_mass_fraction_lower > 1.0 )
      throw logic_error( "RelEffInput: Mass fraction lower constraint is greater than 1." );

    if( constraint.m_mass_fraction_upper < 0.0 )
      throw logic_error( "RelEffInput: Mass fraction upper constraint is less than 0." );

    if( constraint.m_mass_fraction_upper > 1.0 )
      throw logic_error( "RelEffInput: Mass fraction upper constraint is greater than 1." );

    if( (constraint.m_mass_fraction_lower == constraint.m_mass_fraction_upper)
       && ((constraint.m_mass_fraction_lower == 0.0) || (constraint.m_mass_fraction_lower == 1.0)))
      throw logic_error( "RelEffInput: Mass fraction is fixed and equal to 0 or 1 - not allowed." );

    if( constraint.m_mass_fraction_lower > constraint.m_mass_fraction_upper )
      throw logic_error( "RelEffInput: Mass fraction constraint upper value is less than lower value." );

    // Check that the nuclide is in the nuclides list
    // Check that the constrained nuclide is a nuclide in this RelEffCurve
    const auto check_iso_in_curve = [&]( const string &iso ) -> bool
    {
      for( const GenericPeakInfo &peak : peaks )
      {
        for( const GenericLineInfo &line : peak.m_source_gammas )
        {
          if( iso == line.m_isotope )
            return true;
        }//for( const GenericLineInfo &line : peak.m_source_gammas )
      }//for( const GenericPeakInfo &peak : peaks )

      return false;
    };//const auto check_iso_in_curve = [&]( const string &iso ) -> bool

    if( !check_iso_in_curve( constraint.m_nuclide ) )
      throw logic_error( "RelEffInput: Mass fraction constraint nuclide is not in any peak." );

    if( constraint.m_specific_activities.size() < 2 )
      throw logic_error( "RelEffInput: Mass fraction constraint has less than 2 specific activities." );

    // Check each nuclide with a specific activity is in the curve, and has a positive specific activity,
    //  and that this nuclide has a specific activity.
    bool has_this_nuc = false;
    for( const auto &specific_activity : constraint.m_specific_activities )
    {
      if( (specific_activity.second <= 0.0)
          || isnan(specific_activity.second) || isinf(specific_activity.second) )
        throw logic_error( "RelEffInput: Mass fraction constraint specific activity must be a"
                           " positive, finite number." );

      if( !check_iso_in_curve( specific_activity.first ) )
        throw logic_error( "RelEffInput: Mass fraction constraint specific activity nuclide is not in any peak." );

      has_this_nuc |= (specific_activity.first == constraint.m_nuclide);
    }//for( const auto &specific_activity : constraint.m_specific_activities )

    if( !has_this_nuc )
      throw logic_error( "RelEffInput: Mass fraction constraint nuclide is not in the specific activity list." );

    // Check that any other constraint for nuclides in `constraint.m_specific_activities`, has the same `m_specific_activities`.
    for( const RelActCalcManual::MassFractionConstraint &other_constraint : mass_fraction_constraints )
    {
      if( other_constraint.m_nuclide == constraint.m_nuclide )
        continue;
      
      const auto other_pos = other_constraint.m_specific_activities.find( constraint.m_nuclide );
      if( other_pos != end(other_constraint.m_specific_activities) )
      {
        //`other_constraint.m_specific_activities` should be the same as `constraint.m_specific_activities`
        if( other_constraint.m_specific_activities.size() != constraint.m_specific_activities.size() )
          throw logic_error( "RelEffInput: Mass fraction constraint nuclide is in another mass fraction constraint with a different number of specific activities." );

        for( const auto &specific_activity : constraint.m_specific_activities )
        {
          const auto other_const_nuc_pos = other_constraint.m_specific_activities.find( specific_activity.first );
          if( other_const_nuc_pos == end(other_constraint.m_specific_activities) )
            throw logic_error( "RelEffInput: Mass fraction constraint nuclide is in another mass fraction constraint with a different list of specific activities." );

          if( fabs(other_const_nuc_pos->second - specific_activity.second) > 1.0e-6*std::max(other_const_nuc_pos->second, specific_activity.second) )
            throw logic_error( "RelEffInput: Mass fraction constraint nuclide is in another mass fraction constraint with a different specific activity." );
        }//for( const auto &specific_activity : constraint.m_specific_activities )
      }//if( other_pos != end(other_constraint.m_specific_activities) )
    }//for( const RelActCalcManual::MassFractionConstraint &other_constraint : mass_fraction_constraints )


    // Feasibility of the elements constrained windows: with unconstrained isotopes remaining, the
    //  constrained lower bounds must leave them a positive mass remainder; an all-constrained
    //  element (supported via the sigma-block element-scale parameter) instead needs its windows
    //  able to sum to exactly 1.  (Mirrors RelActCalcAuto::RelEffCurveInput::check_nuclide_constraints.)
    size_t num_not_mass_frac_constrained = 0;
    double sum_lower_mass_frac_constrained = 0.0, sum_upper_mass_frac_constrained = 0.0;
    for( const auto &specific_activity : constraint.m_specific_activities )
    {
      const auto mass_frac_pos = std::find_if( begin(mass_fraction_constraints), end(mass_fraction_constraints),
        [&]( const RelActCalcManual::MassFractionConstraint &mfc )
        {
          return mfc.m_nuclide == specific_activity.first;
        } );

      if( mass_frac_pos == end(mass_fraction_constraints) )
      {
        ++num_not_mass_frac_constrained;
      }else
      {
        sum_lower_mass_frac_constrained += mass_frac_pos->m_mass_fraction_lower;
        sum_upper_mass_frac_constrained += mass_frac_pos->m_mass_fraction_upper;
      }
    }//for( const auto &specific_activity : constraint.m_specific_activities )

    if( num_not_mass_frac_constrained == 0 )
    {
      if( (sum_lower_mass_frac_constrained > (1.0 + 1.0E-6))
          || (sum_upper_mass_frac_constrained < (1.0 - 1.0E-6)) )
        throw logic_error( "RelEffInput: All of an elements isotopes are mass-fraction constrained,"
                           " but the windows can not sum to exactly 1 (sum of lower limits="
                           + std::to_string(sum_lower_mass_frac_constrained) + ", sum of upper limits="
                           + std::to_string(sum_upper_mass_frac_constrained) + ")." );
    }else
    {
      if( sum_lower_mass_frac_constrained >= (1.0 - 1.0E-6) )
        throw logic_error( "RelEffInput: The sum of lower mass fraction constraints leaves no mass"
                           " for the elements unconstrained isotopes." );
    }//if( all-constrained element ) / else
  }//for( const RelActCalcManual::MassFractionConstraint &constraint : mass_fraction_constraints )

#endif // USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT


}//RelEffInput::check_nuclide_constraints()

double RelEffSolution::relative_activity( const std::string &nuclide ) const
{
  for( const IsotopeRelativeActivity &i : m_rel_activities )
  {
    if( i.m_isotope == nuclide )
      return i.m_rel_activity;
  }
  
  throw runtime_error( "RelEffSolution::relative_activity: no nuclide '" + nuclide + "'" );
  return 0.0;
}//double relative_activity( const std::string &nuclide ) const


double RelEffSolution::relative_activity_uncertainty( const std::string &nuclide ) const
{
  for( const IsotopeRelativeActivity &i : m_rel_activities )
  {
    if( i.m_isotope == nuclide )
    {
      if( i.m_rel_activity_uncert < 0.0 )
        throw runtime_error( "RelEffSolution::relative_activity_uncertainty: uncertainties not able to be computed." );
      return i.m_rel_activity_uncert;
    }
  }

  throw runtime_error( "RelEffSolution::relative_activity_uncertainty: no nuclide '" + nuclide + "'" );
  return 0.0;
}//double relative_activity_uncertainty( const std::string &nuclide ) const


double RelEffSolution::relative_efficiency( const double energy ) const
{
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return eval_eqn( energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );

  const ManualGenericRelActFunctor::PhysModelRelEqnDef eqn_input
        = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );

  return RelActCalc::eval_physical_model_eqn( energy, eqn_input.self_atten, 
                              eqn_input.external_attens, eqn_input.det.get(), 
                              eqn_input.hoerl_b, eqn_input.hoerl_c );
}//double relative_efficiency( const double energy ) const


size_t RelEffSolution::nuclide_index( const std::string &nuc ) const
{
  for( size_t i = 0; i < m_rel_activities.size(); ++i )
  {
    if( m_rel_activities[i].m_isotope == nuc )
      return i;
  }
  
  throw runtime_error( "RelEffSolution::nuclide_index: invalid nuclide '" + nuc + "'" );
  return m_rel_activities.size(); //prevent warnings.
}//nuclide_index(...)


double RelEffSolution::activity_ratio( const size_t iso1, const size_t iso2 ) const
{
  assert( iso1 < m_rel_activities.size() );
  assert( iso2 < m_rel_activities.size() );
  return m_rel_activities[iso1].m_rel_activity / m_rel_activities[iso2].m_rel_activity;
}

double RelEffSolution::activity_ratio( const std::string &nuc1, const std::string &nuc2 ) const
{
  return activity_ratio( nuclide_index(nuc1), nuclide_index(nuc2) );
}

bool RelEffSolution::walk_to_controlling_nuclide( size_t &iso_index, double &multiple ) const
{
  assert( multiple == 1.0 );
  multiple = 1.0;
  assert( iso_index < m_rel_activities.size() );

#ifndef NDEBUG
  const size_t original_iso_index = iso_index;
#endif

  if( m_input.act_ratio_constraints.empty() )
    return false;

  size_t controller_index = std::numeric_limits<size_t>::max();

  bool found_controller = false;
  size_t sentinel = 0; //Dont need, but just to check the logic for development

  while( controller_index != iso_index )
  {
    sentinel += 1;
    assert( sentinel < 100 );
    if( sentinel > 1000 )
      throw logic_error( "RelEffSolution::activity_ratio_uncert: possible infinite loop - logic error" );

    controller_index = iso_index;
    for( const ManualActRatioConstraint &constraint : m_input.act_ratio_constraints )
    {
      if( constraint.m_constrained_nuclide == m_rel_activities[iso_index].m_isotope )
      {
        iso_index = nuclide_index(constraint.m_controlling_nuclide );
        multiple *= constraint.m_constrained_to_controlled_activity_ratio;
        found_controller = true;
        break;
      }
    }
  }//while( controller_index != iso_index )

  
#ifndef NDEBUG
  // A controlled nuclide's uncertainty is its controller's, times the chain multiple.  This now
  //  comes out of two separate J*C*J^T accumulations rather than being exact by construction, so
  //  compare relatively; and it says nothing when no covariance was computed (both are -1).
  if( (m_rel_activities[iso_index].m_rel_activity_uncert > 0.0)
     && (m_rel_activities[original_iso_index].m_rel_activity_uncert > 0.0) )
  {
    const double expected = multiple * m_rel_activities[iso_index].m_rel_activity_uncert;
    const double actual = m_rel_activities[original_iso_index].m_rel_activity_uncert;
    assert( fabs(expected - actual) <= 1.0E-6*(std::max)(fabs(expected), fabs(actual)) );
  }
#endif

  return found_controller;
}//bool walk_to_controlling_nuclide( size_t &iso_index, double &multiple ) const;


double RelEffSolution::activity_ratio_uncert( const size_t iso1_index, const size_t iso2_index ) const
{
  if( (iso1_index == iso2_index)
     || (iso1_index >= m_rel_activities.size())
     || (iso2_index >= m_rel_activities.size()) )
    throw runtime_error( "RelEffSolution::activity_ratio_uncert: invalid isotope index" );

  if( m_rel_act_covariance.empty() )
    throw runtime_error( "RelEffSolution::activity_ratio_uncert: covariance not available" );

  assert( m_rel_act_covariance.size() == m_rel_activities.size() );

  // Nuclides tied (possibly through a chain) to the same ultimate controller have an exactly
  //  fixed ratio, with zero uncertainty.
  if( !m_input.act_ratio_constraints.empty() )
  {
    double iso1_mult = 1.0, iso2_mult = 1.0;
    size_t iso1 = iso1_index, iso2 = iso2_index;
    walk_to_controlling_nuclide( iso1, iso1_mult );
    walk_to_controlling_nuclide( iso2, iso2_mult );

    if( iso1 == iso2 )
      return 0.0;
  }//if( !m_input.act_ratio_constraints.empty() )

  // `m_rel_act_covariance` is the full-Jacobian covariance of the final relative activities, so
  //  activity norms, act-ratio constraint chains, and mass-fraction constraints are all already
  //  accounted for; first-order propagation of the ratio is all that is left.
  const double act_1 = m_rel_activities[iso1_index].m_rel_activity;
  const double act_2 = m_rel_activities[iso2_index].m_rel_activity;

  // Relative activities are bounded below by zero, so a nuclide with no real signal can land
  //  exactly at 0 - the relative-uncertainty form below would then produce inf/NaN rather than
  //  failing, and NaN would flow all the way to the displayed "+-" (callers guard by catching).
  if( (act_1 <= 0.0) || (act_2 <= 0.0) )
    throw runtime_error( "RelEffSolution::activity_ratio_uncert: a relative activity is zero,"
                         " so their ratio's uncertainty is not defined" );

  const double ratio = act_1 / act_2;

  const double cov_1_1 = m_rel_act_covariance[iso1_index][iso1_index];
  const double cov_1_2 = m_rel_act_covariance[iso1_index][iso2_index];
  const double cov_2_2 = m_rel_act_covariance[iso2_index][iso2_index];

  // The correlation term can drive the radicand negative for strongly correlated isotopes;
  // clamp to zero rather than producing NaN.
  return fabs(ratio) * sqrt( std::max( 0.0,
                    (cov_1_1 / (act_1 * act_1))
                    + (cov_2_2 / (act_2 * act_2))
                    - (2.0 * cov_1_2 / (act_1 * act_2)) ) );
}


double RelEffSolution::activity_ratio_uncert( const std::string &iso1, const std::string &iso2 ) const
{
  return activity_ratio_uncert( nuclide_index(iso1), nuclide_index(iso2) );
}


double RelEffSolution::mass_fraction( const std::string &nuclide ) const
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *wanted_nuc = db->nuclide( nuclide );
  if( !wanted_nuc )
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "'): invalid nuclide" );
  
  double sum_rel_mass = 0.0, nuc_rel_mas = -1.0;
  for( size_t index = 0; index < m_rel_activities.size(); ++index )
  {
    const IsotopeRelativeActivity &act = m_rel_activities[index];
    const SandiaDecay::Nuclide * const nuc = db->nuclide( act.m_isotope );
    if( !nuc )
      continue;
    
    const double rel_mass = act.m_rel_activity / nuc->activityPerGram();
    sum_rel_mass += rel_mass;
    
    if( nuc == wanted_nuc )
    {
      assert( nuc_rel_mas == -1.0 );
      nuc_rel_mas = rel_mass;
    }
  }//for( size_t index = 0; index < m_rel_activities.size(); ++index )
  
  if( nuc_rel_mas < 0.0 )
    throw runtime_error( "mass_fraction: invalid nuclide: " + nuclide );
  
  return nuc_rel_mas / sum_rel_mass;
}//double mass_fraction( const std::string &nuclide ) const

  
double RelEffSolution::mass_fraction( const std::string &nuclide, const double num_sigma ) const
{
  // The `double RelActAutoSolution::mass_enrichment_fraction(...)` function is similar to
  //  this one, so if any issues are found in this function, please also check that function.

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *wanted_nuc = db->nuclide( nuclide );
  
  if( !wanted_nuc ) //Will be nullptr for reactions and x-rays
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "', num_sigma): invalid nuclide" );
  
  //assert( !m_rel_act_covariance.empty() ); // Failing to compute the activity covarances can happen sometimes
  if( m_rel_act_covariance.empty() )
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "', num_sigma): no valid covariance." );

  assert( m_rel_act_covariance.size() == m_rel_activities.size() );

  const size_t nuc_index = std::find_if( begin(m_rel_activities), end(m_rel_activities),
                          [&nuclide]( const IsotopeRelativeActivity &val ) {
    return val.m_isotope == nuclide;
  }) - begin(m_rel_activities);

  // Note: asking for a nuclide that isnt in the solution is a legitimate query (callers ask for
  //  e.g. "U235" without knowing what was fit, inside a try/catch), so this must throw - dont
  //  assert on it ahead of the check, which would abort debug builds.
  if( nuc_index >= m_rel_activities.size() )
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "', "
                        + std::to_string(num_sigma) + "): nuclide not in solution set" );

  assert( nuc_index < m_rel_act_covariance.size() );

  // Vary the target nuclides activity by num_sigma of its uncertainty, moving every other
  //  activity along its regression slope cov_kn/var_n (the conditional mean of the joint
  //  Gaussian) - a linearized profile.  `m_rel_act_covariance` is the full-Jacobian covariance
  //  of the final relative activities, so all constraint co-movements (act-ratio chains,
  //  mass-fraction block decodes) are already encoded in the slopes: e.g., a fixed
  //  mass-fraction nuclide automatically stays at its pinned within-element fraction, since its
  //  activity co-varies proportionally with the rest of its element.
  const double var_n = m_rel_act_covariance[nuc_index][nuc_index];
  if( var_n <= 0.0 ) // an exactly-pinned (or degenerate) activity: no variation to propagate
    return mass_fraction( nuclide );

  const double sqrt_var_n = sqrt( var_n );

  double sum_rel_mass = 0.0, nuc_rel_mas = -1.0;
  for( size_t index = 0; index < m_rel_activities.size(); ++index )
  {
    const IsotopeRelativeActivity &act = m_rel_activities[index];
    const SandiaDecay::Nuclide * const nuc = db->nuclide( act.m_isotope );
    if( !nuc ) //For example when an x-ray or reaction
      continue;

    const double shifted_act = act.m_rel_activity
                               + (num_sigma * m_rel_act_covariance[index][nuc_index] / sqrt_var_n);
    const double rel_mass = shifted_act / nuc->activityPerGram();

    sum_rel_mass += (std::max)( rel_mass, 0.0 );

    if( index == nuc_index )
    {
      assert( nuc == wanted_nuc );

      nuc_rel_mas = rel_mass;
    }//if( index == nuc_index )
  }//for( size_t index = 0; index < m_rel_activities.size(); ++index )

  if( nuc_rel_mas < 0.0 ) // This can happen when we go down a couple sigma
    return 0.0;

  return nuc_rel_mas / sum_rel_mass;
}//double mass_fraction( const std::string &iso, const double num_sigma ) const
    
  
std::string RelEffSolution::parameter_name( const size_t par_num ) const
{
  if( par_num >= m_fit_parameters.size() )
    throw runtime_error( "RelEffSolution::parameter_name: invalid parameter." );
  
  if( par_num < m_rel_activities.size() )
    return "Act(" + m_rel_activities[par_num].m_isotope + ")";
  
  if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    size_t working_par_num = m_rel_activities.size();
    if( working_par_num )
      working_par_num -= 1;
    
    if( shield_is_present( m_input.phys_model_self_atten ) )
    {
      if( shield_fits_an( m_input.phys_model_self_atten ) )
      {
        working_par_num += 1;
        if( working_par_num == par_num )
          return "SAtt(AN)";
      }

      working_par_num += 1;
      if( working_par_num == par_num )
        return "SAtt(AD)";
    }//if( use internal attenuation shielding )

    for( size_t i = 0; i < m_input.phys_model_external_attens.size(); ++i )
    {
      const auto &ext_atten = m_input.phys_model_external_attens[i];
      if( !shield_is_present( ext_atten ) )
        continue;

      if( shield_fits_an( ext_atten ) )
      {
        working_par_num += 1;
        if( working_par_num == par_num )
          return "EAtt" + std::to_string(i) + "(AN)";
      }
      
      working_par_num += 1;
      if( working_par_num == par_num )
        return "EAtt" + std::to_string(i) + "(AD)";
    }//for( loop over external attenuators )
    
    if( m_input.phys_model_use_hoerl )
    {
      working_par_num += 1;
      if( working_par_num == par_num )
        return "Hoerl(b)";
      working_par_num += 1;
      if( working_par_num == par_num )
        return "Hoerl(c)";
    }
    
    assert( 0 );
    throw std::logic_error( "Logic for determining Physical Model coefficient name is whack." );
  }//if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  
  
  return "RE_" + std::to_string( par_num - m_rel_activities.size() );
}//string parameter_name( const size_t par_num ) const
  

std::ostream &RelEffSolution::print_summary( std::ostream &strm ) const
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  switch( m_status )
  {
    case RelActCalcManual::ManualSolutionStatus::NotInitialized:
      strm << "Status: NotInitialized\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::ErrorInitializing:
      strm << "Status: ErrorInitializing\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::ErrorFindingSolution:
      strm << "Status: ErrorFindingSolution\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::ErrorGettingSolution:
      strm << "Status: ErrorGettingSolution\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::Success:
      strm << "Status: Success\n";
      break;
  }//switch( solution.m_status )
  
  if( !m_error_message.empty() )
  {
    strm << "--------------------------------------------------------------------------------\n";
    strm << "Error: " << m_error_message << "\n";
    strm << "--------------------------------------------------------------------------------\n";
  }
  
  if( !m_warnings.empty() )
  {
    strm << "--------------------------------------------------------------------------------\n";
    for( const string &warning : m_warnings )
      strm << "Warning: " << warning << "\n";
    strm << "--------------------------------------------------------------------------------\n";
  }//if( !m_warnings.empty() )
  
  strm << "Eqn coefficients: ";
  for( size_t i = 0; i < m_rel_eff_eqn_coefficients.size(); ++i )
  strm << (!i ? "" : ", ") << m_rel_eff_eqn_coefficients[i];
  strm << endl;
  strm << "Eqn coefficient covariance:\n";
  for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
  {
    for( size_t j = 0; j < m_rel_eff_eqn_covariance[i].size(); ++j )
    strm << std::setw(14) << m_rel_eff_eqn_covariance[i][j];
    strm << endl;
  }
  strm << endl;
  
  strm << "Relative activities:" << endl;
  for( const RelActCalcManual::IsotopeRelativeActivity &i : m_rel_activities )
    strm << "\t" << i.m_isotope << ": " << i.m_rel_activity << " +- " << i.m_rel_activity_uncert << endl;

  
  strm << "Relative masses:" << endl;
  for( const RelActCalcManual::IsotopeRelativeActivity &i : m_rel_activities )
  {
    const SandiaDecay::Nuclide *nuclide = db->nuclide( i.m_isotope );
    
    if( nuclide )
    {
      const double rel_mass = mass_fraction( i.m_isotope );
      strm << "\t" << i.m_isotope << ": " << 100.0*rel_mass;

      try
      {
        const double rel_frac_plus = mass_fraction( i.m_isotope, 1.0 );
        const double rel_frac_minus = mass_fraction( i.m_isotope, -1.0 );
        const double uncert = 0.5*(rel_frac_plus - rel_frac_minus);
        strm << " +- " << uncert;
      }catch( std::exception &e )
      {

      }

      strm << endl;
    }
  }
  
  
  strm << "Chi2: " << m_chi2 << endl;
  strm << "Num Fcnt Evals: " << m_num_function_eval_total << "\n\n" << endl;
  
  
  for( const auto &peak: m_input.peaks )
  {
    double function_val = std::numeric_limits<double>::quiet_NaN();
    if( !m_rel_eff_eqn_coefficients.empty() )
    {
      if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        function_val = RelActCalc::eval_eqn( peak.m_energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );
      }else
      {
        const ManualGenericRelActFunctor::PhysModelRelEqnDef input
              = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
        
        function_val = RelActCalc::eval_physical_model_eqn( peak.m_energy, input.self_atten,
                            input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
      }//if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else
    }//if( !m_rel_eff_eqn_coefficients.empty() )

    strm << "For energy " << peak.m_energy << " (";
    for( size_t i = 0; i < peak.m_source_gammas.size(); ++i )
    strm << (i ? ", " : "") << peak.m_source_gammas[i].m_isotope;
    strm << ") function value is " << function_val << ", and:" << endl;
    
    //peak.m_counts
    strm << "\t";
    double total_contrib_counts = 0.0;
    for( size_t i = 0; i < peak.m_source_gammas.size(); ++i )
    {
      const RelActCalcManual::GenericLineInfo &line = peak.m_source_gammas[i];
      
      const string &iso = line.m_isotope;
      const double yield = line.m_yield;
      double rel_act = std::numeric_limits<double>::quiet_NaN();
      if( !m_rel_activities.empty() )
        rel_act = relative_activity(iso);
      
      const double contrib_counts = rel_act * yield * function_val;
      total_contrib_counts += contrib_counts;
      strm << (i ? ", " : "") << iso << ": fit " << contrib_counts << " counts";
    }//for( const RelActCalc::GammaLineInfo &line : peak.m_source_gammas )
    
    strm << " and observed " << peak.m_counts << "+-" << peak.m_counts_uncert << " (off by "
    << ((total_contrib_counts - peak.m_counts) / peak.m_counts_uncert)
    << " sigma)" << endl;
    
    
  }//for( const auto &peak: peak_infos )
  
  return strm;
}//std::ostream &print_summary( std::ostream &strm ) const


void RelEffSolution::get_mass_fraction_table( std::ostream &results_html ) const
{
  const int nsigfig = 4;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    return;
  
  results_html << "<table class=\"nuctable resulttable\">\n";
  results_html << "  <caption>Relative activities and mass fractions</caption>\n";
  results_html << "  <thead><tr>"
                 " <th scope=\"col\">Nuclide</th>"
                 " <th scope=\"col\">Rel. Act.</th>"
                 " <th scope=\"col\">Mass Frac.</th>"
                 //" <th scope=\"col\" title=\"This is percentage uncertainty - i.e., the percent error on the reported percent.\">Uncert.</th>"
                 " <th scope=\"col\" >2&sigma; interval</th>"
                 " </tr></thead>\n";
  results_html << "  <tbody>\n";
  
  // We will put the normalized relative activity (i.e., activity divided by largest activity) in
  //  a tool-tip.
  //  However - we should probably just do this for the actual display - I'm just nervous of messing something up somewhere
  //  If we had more room, we could add a "Norm Act." column.
  double largest_rel_act = 0.0;
  for( const auto &act : m_rel_activities )
    largest_rel_act = std::max( largest_rel_act, act.m_rel_activity );
  
  for( size_t index = 0; index < m_rel_activities.size(); ++index )
  {
    const IsotopeRelativeActivity &act = m_rel_activities[index];

    results_html << "  <tr><td>" << act.m_isotope << "</td>"
    << "<td title=\"The normalized relative activity (i.e., divided by largest rel. act.) is "
    << SpecUtils::printCompact( act.m_rel_activity / largest_rel_act, nsigfig + 1)
    << "\">" << SpecUtils::printCompact( act.m_rel_activity, nsigfig ) << "</td>";
    
    try
    {
      const double frac_mass = mass_fraction(act.m_isotope);
      results_html << "<td>" << SpecUtils::printCompact(100.0*frac_mass, nsigfig)      << "%</td>";
    }catch( std::exception & )
    {
      results_html << "<td>N.A.</td>";
    }
    
    // Convention: report the asymmetric 1-sigma interval, with the symmetric sigma defined as
    //  the interval half-width, 0.5*(|delta_plus| + |delta_minus|) - the same convention the
    //  chart-title enrichment uses.  (Check the "Profile uncert." checkbox of a nuclide for a
    //  profile-likelihood interval, which is also valid near the physical bounds.)
    string error_tt;
    try
    {
      const double frac_mass = mass_fraction(act.m_isotope);

      const double frac_mass_plus1 = mass_fraction(act.m_isotope, 1.0);
      const double frac_mass_minus1 = mass_fraction(act.m_isotope, -1.0);

      const double delta_plus1 = frac_mass_plus1 - frac_mass;
      const double delta_minus1 = frac_mass_minus1 - frac_mass;
      const double half_width = 0.5*( fabs(delta_plus1) + fabs(delta_minus1) );
      const double uncert_percent = 100.0 * half_width / frac_mass;

      error_tt = "The 1-sigma mass fraction interval is ["
      + SpecUtils::printCompact(100.0*frac_mass_minus1, nsigfig+1)
      + "% to "
      + SpecUtils::printCompact(100.0*frac_mass_plus1, nsigfig+1)
      + "%]. \nThe symmetric 1-sigma percentage uncertainty (interval half-width, as percent of"
      " the reported percent) is "
      + SpecUtils::printCompact(uncert_percent, nsigfig+1)
      + "%";
    }catch( std::exception & )
    {
      // If we dont have covariance matrix, we will end up here
      error_tt = "No covariance available.";
    }
    
    results_html << "<td title=\"" << error_tt << "\">" ;
    //<< SpecUtils::printCompact(uncert_percent, nsigfig-1)
    try
    {
      const double frac_mass_plus2 = mass_fraction(act.m_isotope, 2.0);
      const double frac_mass_minus2 = mass_fraction(act.m_isotope, -2.0);

      results_html << SpecUtils::printCompact(100.0*frac_mass_minus2, nsigfig-1)
      << "%, "
      << SpecUtils::printCompact( 100.0*frac_mass_plus2, nsigfig-1)
      << "%";
    }catch( std::exception & )
    {
      // We will get here for reactions and x-rays.
    }
    
    results_html << "</td>"
    << "</tr>\n";
  }
  results_html << "  </tbody>\n"
  << "</table>\n\n";

  // Profile-likelihood intervals, when they were asked for (see RelEffInput::profile_targets);
  //  these lines serve both the in-app results tab and the HTML report.
  for( const ProfileResult &profile : m_profile_results )
  {
    if( profile.intervals.empty() )
      continue;

    const bool is_mass_frac = (profile.target.m_type == ProfileTarget::Type::MassFraction);
    const double disp_mult = is_mass_frac ? 100.0 : 1.0;
    const char * const units = is_mass_frac ? "%" : "";

    results_html << "<div class=\"profileuncert\">Profile likelihood " << profile.target.m_nuclide;
    switch( profile.target.m_type )
    {
      case ProfileTarget::Type::MassFraction:
        results_html << " mass fraction (of its element): ";
        break;
      case ProfileTarget::Type::ActivityRatio:
        results_html << "/" << profile.target.m_denom_nuclide << " activity ratio: ";
        break;
      case ProfileTarget::Type::RelativeActivity:
        results_html << " relative activity: ";
        break;
    }//switch( profile.target.m_type )
    results_html << SpecUtils::printCompact(disp_mult*profile.nominal_value, 4) << units;

    for( const ProfileInterval &interval : profile.intervals )
    {
      // Mark the individual end that is not a chi2 crossing, rather than tagging the whole
      //  interval: that end is only as far as the scan got, so the true end-point lies beyond it -
      //  hence the inequality - and which of the two it is changes how the number should be read.
      results_html << "; " << SpecUtils::printCompact(100.0*interval.confidence_level, 4)
                   << "% CL in ["
                   << (interval.lower_at_bound ? "&lt;" : "")
                   << SpecUtils::printCompact(disp_mult*interval.lower_value, 4) << units
                   << (interval.lower_at_bound ? " (limit)" : "")
                   << ", "
                   << (interval.upper_at_bound ? "&gt;" : "")
                   << SpecUtils::printCompact(disp_mult*interval.upper_value, 4) << units
                   << (interval.upper_at_bound ? " (limit)" : "")
                   << "]";
    }//for( const ProfileInterval &interval : profile.intervals )

    results_html << "</div>\n";

    for( const std::string &warning : profile.warnings )
      results_html << "<div class=\"profileuncert profileuncertwarn\">" << warning << "</div>\n";
  }//for( const ProfileResult &profile : m_profile_results )
}//void get_mass_fraction_table( std::ostream &strm ) const


void RelEffSolution::get_phys_model_shield_text( std::ostream &strm ) const
{
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return;

  // The shield AD/AN uncertainties (computed from the parameter covariance, including the
  //  chi2/dof inflation) were previously computed but shown nowhere.
  const auto print_shield = [&strm]( const PhysModelShieldFit &shield, const std::string &label ){
    strm << "<div class=\"shieldfit\">" << label << ": ";
    if( shield.m_material )
      strm << shield.m_material->name;
    else
    {
      strm << "AN=" << SpecUtils::printCompact( shield.m_atomic_number, 4 );
      if( shield.m_atomic_number_uncert > 0.0 )
        strm << " &plusmn; " << SpecUtils::printCompact( shield.m_atomic_number_uncert, 3 );
    }

    strm << ", AD=" << SpecUtils::printCompact( shield.m_areal_density/PhysicalUnits::g_per_cm2, 4 );
    if( shield.m_areal_density_uncert > 0.0 )
      strm << " &plusmn; " << SpecUtils::printCompact( shield.m_areal_density_uncert/PhysicalUnits::g_per_cm2, 3 );
    strm << " g/cm<sup>2</sup></div>\n";
  };//print_shield

  if( m_phys_model_self_atten_shield )
    print_shield( *m_phys_model_self_atten_shield, "Self-atten" );

  for( size_t i = 0; i < m_phys_model_external_atten_shields.size(); ++i )
  {
    if( m_phys_model_external_atten_shields[i] )
      print_shield( *m_phys_model_external_atten_shields[i], "Ext-atten " + std::to_string(i+1) );
  }
}//void get_phys_model_shield_text( std::ostream &strm ) const


void RelEffSolution::get_mass_ratio_table( std::ostream &results_html ) const
{
  const size_t nsigfig = 4;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    return;
  
  // Make the table of mass and activity ratios
  results_html << "<table class=\"massratiotable resulttable\">\n";
  results_html << "  <caption>Mass and Activity Ratios.</caption>\n";
  results_html << "  <thead><tr>"
  "<th scope=\"col\">Nuclides</th>"
  "<th scope=\"col\">Mass Ratio</th>"
  "<th scope=\"col\">Activity Ratio</th>"
  "</tr></thead>\n";
  results_html << "  <tbody>\n";
  
  for( size_t i = 1; i < m_rel_activities.size(); ++i )
  {
    for( size_t j = 0; j < i; ++j )
    {
      const string &nuc_i_str = m_rel_activities[i].m_isotope;
      const string &nuc_j_str = m_rel_activities[j].m_isotope;
      const SandiaDecay::Nuclide * const nuc_i = db->nuclide( nuc_i_str );
      const SandiaDecay::Nuclide * const nuc_j = db->nuclide( nuc_j_str );
      
      if( !nuc_i || !nuc_j )
      {
        results_html << "<tr><td>" << nuc_i_str << "/" << nuc_j_str << "</td><td>--</td><td>--</td></tr>\n";
        results_html << "<tr><td>" << nuc_j_str << "/" << nuc_i_str << "</td><td>--</td><td>--</td></tr>\n";
        continue;
      }
      
      const double act_i = relative_activity( nuc_i_str );
      const double act_j = relative_activity( nuc_j_str );
      
      const double mass_i = act_i / nuc_i->activityPerGram();
      const double mass_j = act_j / nuc_j->activityPerGram();
      
      const double i_to_j_specific_act = nuc_i->activityPerGram() / nuc_j->activityPerGram();
      const double j_to_i_specific_act = 1.0 / i_to_j_specific_act;
      
      const double i_to_j_act_ratio = activity_ratio( nuc_i_str, nuc_j_str );
      const double j_to_i_act_ratio = activity_ratio( nuc_j_str, nuc_i_str );
      
      const double i_to_j_mass_ratio = i_to_j_act_ratio * j_to_i_specific_act;
      const double j_to_i_mass_ratio = j_to_i_act_ratio * i_to_j_specific_act;
      
      // activity_ratio_uncert(...) throws when the covariance was not successfully computed
      //  (or an index is invalid); fall back to the "n/a" cells in that case.
      bool have_uncerts = false;
      double i_to_j_act_ratio_uncert = -1.0, j_to_i_act_ratio_uncert = -1.0;
      try
      {
        i_to_j_act_ratio_uncert = activity_ratio_uncert( nuc_i_str, nuc_j_str );
        j_to_i_act_ratio_uncert = activity_ratio_uncert( nuc_j_str, nuc_i_str );
        have_uncerts = true;
      }catch( std::exception & )
      {
      }

      if( have_uncerts )
      {
        const double i_to_j_mass_ratio_uncert = i_to_j_act_ratio_uncert * j_to_i_specific_act;
        const double j_to_i_mass_ratio_uncert = j_to_i_act_ratio_uncert * i_to_j_specific_act;

        results_html << "<tr><td>" << nuc_i->symbol << "/" << nuc_j->symbol
        << "</td><td>" << PhysicalUnits::printValueWithUncertainty( i_to_j_mass_ratio, i_to_j_mass_ratio_uncert, nsigfig )
        << "</td><td>" << PhysicalUnits::printValueWithUncertainty( i_to_j_act_ratio, i_to_j_act_ratio_uncert, nsigfig )
        << "</td></tr>\n";
        
        results_html << "<tr><td>" << nuc_j->symbol << "/" << nuc_i->symbol
        << "</td><td>" << PhysicalUnits::printValueWithUncertainty( j_to_i_mass_ratio, j_to_i_mass_ratio_uncert, nsigfig )
        << "</td><td>"<< PhysicalUnits::printValueWithUncertainty( j_to_i_act_ratio, j_to_i_act_ratio_uncert, nsigfig )
        << "</td></tr>\n";
      }else
      {
        results_html << "<tr><td>" << nuc_i->symbol << "/" << nuc_j->symbol
        << "</td><td>" << SpecUtils::printCompact(i_to_j_mass_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td><td>" << SpecUtils::printCompact(i_to_j_act_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td></tr>\n";
        
        results_html << "<tr><td>" <<nuc_j->symbol << "/" << nuc_i->symbol
        << "</td><td>" << SpecUtils::printCompact(j_to_i_mass_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td><td>" << SpecUtils::printCompact(j_to_i_act_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td></tr>\n";
      }//if( we have covariance ) / else
    }//for( size_t j = 0; j < i; ++j )
  }//for( size_t i = 0; i < used_isotopes.size(); ++i )
  results_html << "  </tbody>\n"
  << "</table>\n\n";
}//void get_mass_ratio_table( std::ostream &strm ) const


double RelEffSolution::rel_eff_eqn_value( const double energy ) const
{
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return RelActCalc::eval_eqn( energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );
  
  const ManualGenericRelActFunctor::PhysModelRelEqnDef input 
                = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );

  return RelActCalc::eval_physical_model_eqn_imp<double>( energy, input.self_atten, input.external_attens,
                                            input.det.get(), input.hoerl_b, input.hoerl_c );
}


double RelEffSolution::rel_eff_eqn_uncert( const double energy ) const
{
  if( m_rel_eff_eqn_covariance.empty() )
    throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: Rel. Eff. Eqn. covariances not available." );

  assert( m_rel_eff_eqn_covariance.size() == m_rel_eff_eqn_coefficients.size() );
  if( m_rel_eff_eqn_covariance.size() != m_rel_eff_eqn_coefficients.size() )
    throw std::logic_error( "RelEffSolution::rel_eff_eqn_uncert: covariance matrix does not match expected." );

  // I think we would be safe skipping this following check, at least on non-debug builds, but whatever
  for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
  {
    assert( m_rel_eff_eqn_covariance[i].size() == m_rel_eff_eqn_covariance.size() );
    if( m_rel_eff_eqn_covariance[i].size() != m_rel_eff_eqn_covariance.size() )  //JIC for release builds
      throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: covariance not a square matrix." );
  }//for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )

  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    // We can delegate to RelActCalc::eval_eqn_uncertainty() no matter whether the covariance
    //  came from the log-space LLS fit (use_ceres_to_fit_eqn == false), or from Ceres fitting
    //  the coefficients directly: LnX is linear in its coefficients, and for the exponential
    //  forms (LnY, LnXLnY, FramEmpirical) ln(y) is linear in the *same* coefficients in both
    //  parameterizations, so sigma_y = y*sqrt(t^T*C*t) is the identical quadratic form either
    //  way (the non-log Jacobian factors, y*t_i, just pull y^2 outside the double sum).
    const double val = RelActCalc::eval_eqn_uncertainty( energy, m_input.eqn_form,
                                            m_rel_eff_eqn_coefficients, m_rel_eff_eqn_covariance );
    if( isnan(val) )
      throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: NaN value for uncertainty." );

    return val;
  }//if( an empirical eqn form )

  // FramPhysicalModel: compute the gradient of the rel. eff. curve with respect to the
  //  Ceres-space parameters (i.e., AN/ns_an_ceres_mult, AD as the g/cm2 number, and the
  //  offset Hoerl b/c values), using automatic differentiation, then contract it with
  //  `m_rel_eff_eqn_covariance` - which is the raw Ceres covariance sub-block in those same
  //  units.  `make_phys_eqn_input` applies all unit transforms internally, exactly like the
  //  cost functor does during the fit, so the chain rule to parameter space is automatic.
  //  Parameters held constant during the fit have all-zero covariance rows/columns, so their
  //  gradient components contribute nothing.
  assert( m_input.use_ceres_to_fit_eqn );

  const size_t num_pars = m_rel_eff_eqn_coefficients.size();
  vector<double> gradient( num_pars, 0.0 );

  for( size_t chunk_start = 0; chunk_start < num_pars; chunk_start += 8 )
  {
    typedef ceres::Jet<double,8> JetType;

    vector<JetType> jet_coefs( num_pars );
    for( size_t i = 0; i < num_pars; ++i )
      jet_coefs[i] = JetType( m_rel_eff_eqn_coefficients[i] );

    const size_t num_seed = std::min( num_pars - chunk_start, size_t(8) );
    for( size_t k = 0; k < num_seed; ++k )
      jet_coefs[chunk_start + k].v[k] = 1.0;

    const ManualGenericRelActFunctor::PhysModelRelEqnDef<JetType> eqn_input
                        = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, jet_coefs );

    const JetType eval_val = RelActCalc::eval_physical_model_eqn_imp<JetType>( energy,
                        eqn_input.self_atten, eqn_input.external_attens, eqn_input.det.get(),
                        eqn_input.hoerl_b, eqn_input.hoerl_c );

    for( size_t k = 0; k < num_seed; ++k )
      gradient[chunk_start + k] = eval_val.v[k];
  }//for( loop over stride-8 chunks of parameters )

  double uncert_sq = 0.0;
  for( size_t i = 0; i < num_pars; ++i )
  {
    for( size_t j = 0; j < num_pars; ++j )
      uncert_sq += gradient[i] * m_rel_eff_eqn_covariance[i][j] * gradient[j];
  }//for( size_t i = 0; i < num_pars; ++i )

  if( uncert_sq < 0.0 )
    throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: negative squared uncertainty." );

  return sqrt( uncert_sq );
}//double rel_eff_eqn_uncert( const double energy ) const


string RelEffSolution::rel_eff_eqn_js_uncert_fcn() const
{
  vector<double> energies;
  double current_energy = 20.0;
  const double upper_energy = m_input.peaks.empty() ? 3000.0 
                                                    : std::max(3000.0, m_input.peaks.back().m_energy);

  for( const GenericPeakInfo &peak : m_input.peaks )
  {
    double min_dx = 1.0;
    if( current_energy < 130 )
      min_dx = 1.0;
    else if( current_energy < 300 )
      min_dx = 5.0;
    else
      min_dx = 15.0;

    // We'll try to get in at least ~10 points between each peak
    if( peak.m_energy > current_energy )
      min_dx = std::min( min_dx, 0.1*(peak.m_energy - current_energy) );
    min_dx = std::max( min_dx, 1.0 ); //but less than a keV between points is just to small.

    for( ; current_energy < peak.m_energy; current_energy += min_dx )
      energies.push_back( current_energy );
    
    if( !energies.empty() && (energies.back() < peak.m_energy) )
      current_energy = peak.m_energy;
  }//for( const GenericPeakInfo &peak : m_input.peaks )
  
  for( ; current_energy < upper_energy; current_energy += 15 )
    energies.push_back( current_energy );
  
  size_t num_points = 0;
  string fcn = "function(x){\n"
  "  const points = [";
  bool is_first_point = true;
  for( double x : energies )
  {
    try
    {
      double y = rel_eff_eqn_uncert( x );
      //assert( (y >= 0.0) && !IsNan(y) && !IsInf(y) ); //Can happen when we are out of bounds of the physical model
      if( isnan(y) || isinf(y) )
        continue;
    
      fcn += is_first_point ? "" : ",";
      fcn += "[" + SpecUtils::printCompact(x, 4) + "," + SpecUtils::printCompact(y, 4) + "]";

      num_points += 1;
      is_first_point = false;
    }catch( std::exception &e )
    {
      // This can happennty when we are out of bounds of the physical model
    }
  }//for( double x : energies )

if( num_points < 2 )
 return "null";

  fcn += "];\n"
  "  if( x <= points[0][0] )\n"
  "    return points[0][1];\n"
  "  if( x >= points[points.length - 1][0] )\n"
  "    return points[points.length - 1][1];\n"
  "  for (let i = 0; i < points.length - 1; i++) {\n"
  "    const [x1, y1] = points[i];\n"
  "    const [x2, y2] = points[i + 1];\n"
  "    if( x >= x1 && x <= x2) {\n"
  "      const t = (x - x1) / (x2 - x1);\n"
  "      return y1 + t * (y2 - y1);\n"
  "    }\n"
  "  }\n"
  "console.assert(0,'Shouldnt get here in interpolating');\n"
  "return points[points.length - 1][1];"
  "}";
  
  return fcn;
}//string RelEffSolution::rel_eff_eqn_js_uncert_fcn() const


string RelEffSolution::rel_eff_eqn_txt( const bool html_format ) const
{
  if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    const ManualGenericRelActFunctor::PhysModelRelEqnDef input
          = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
    return RelActCalc::physical_model_rel_eff_eqn_text( input.self_atten,
                input.external_attens, input.det, input.hoerl_b, input.hoerl_c, html_format );
  }
  
  return RelActCalc::rel_eff_eqn_text( m_input.eqn_form, m_rel_eff_eqn_coefficients );
}//std::string RelEffSolution::rel_eff_eqn_txt( const bool html_format ) const
  
  
string RelEffSolution::rel_eff_eqn_js_function() const
{
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return RelActCalc::rel_eff_eqn_js_function( m_input.eqn_form, m_rel_eff_eqn_coefficients );
  

  const ManualGenericRelActFunctor::PhysModelRelEqnDef input
        = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
    
  return RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                          input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
}//string RelEffSolution::rel_eff_eqn_js_function() const
  

void RelEffSolution::print_html_report( ostream &output_html_file,
                                       string spectrum_title,
                                       shared_ptr<const SpecUtils::Measurement> spectrum,
                                       vector<shared_ptr<const PeakDef>> display_peaks,
                                       shared_ptr<const SpecUtils::Measurement> background,
                                       double background_normalization,
                                       vector<shared_ptr<const PeakDef>> background_peaks ) const
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  const size_t nsigfig = 4;
  
  char buffer[512] = { '\0' };
  
  
  stringstream results_html;
  results_html << "<div>&chi;<sup>2</sup>=" << SpecUtils::printCompact(m_chi2, nsigfig)
  << " and there were " << m_dof << " DOF";
  if( m_dof > 0 ) //m_dof can legitimately be zero
    results_html << "; &chi;<sup>2</sup>/dof=" << SpecUtils::printCompact(m_chi2/m_dof, nsigfig);
  results_html << " </div>\n";
  
  results_html << "<div class=\"releffeqn\">Rel. Eff. Eqn: y = ";
  results_html << rel_eff_eqn_txt( true );

  results_html << "</div>\n";

  get_phys_model_shield_text( results_html );

  get_mass_fraction_table( results_html );
  get_mass_ratio_table( results_html );
  
  const bool has_decay_corr = !m_input.peaks_before_decay_correction.empty();

  bool any_peak_has_multiple_srcs = false;
  for( const GenericPeakInfo &info : m_input.peaks )
    any_peak_has_multiple_srcs |= (info.m_source_gammas.size() > 1);

  // Make table giving info on each of the _used_ peaks
  results_html << "<table class=\"peaktable resulttable\">\n";
  results_html << "  <caption>Peaks used for analysis.</caption>\n";
  results_html << "  <thead><tr>"
  "<th scope=\"col\">Energy (keV)</th>"
  "<th scope=\"col\">Nuclide</th>"
  "<th scope=\"col\">Yield</th>"
  "<th scope=\"col\">Net Area</th>"
  "<th scope=\"col\">Net Area Uncert</th>"
  "<th scope=\"col\">Counts/Yield</th>"
  "<th scope=\"col\">Counts/Yield Unc</th>"
  "<th scope=\"col\">Add. Unc.</th>"
  "<th scope=\"col\">Meas. Rel Eff</th>"
  "<th scope=\"col\">Meas. Rel Eff Unct</th>"
  << (any_peak_has_multiple_srcs ? "<th scope=\"col\">Peak Frac</th>" : "")
  << (has_decay_corr ? "<th scope=\"col\">Decay Corr.</th>" : "")
  << "</tr></thead>\n"
  "  <tbody>\n";
  
  
  for( const GenericPeakInfo &info : m_input.peaks )
  {
    snprintf(buffer, sizeof(buffer), "%.2f", info.m_mean );
    results_html << "  <tr><td>" << buffer << "</td>";
    for( size_t i = 0; i < info.m_source_gammas.size(); ++i )
    {
      const GenericLineInfo &line = info.m_source_gammas[i];
      
      if( i )
        results_html << "<tr><td></td>";
      
      const double rel_act = relative_activity(line.m_isotope);
      const double counts_over_yield = info.m_counts / line.m_yield;
      const double counts_uncert_percent = 100.0*(info.m_counts_uncert / info.m_counts);
      const double meas_rel_eff = info.m_counts / (line.m_yield * rel_act);
      const double meas_rel_eff_uncert = 100* info.m_counts_uncert / info.m_counts;
      
      results_html << "<td>" << line.m_isotope
      << "</td><td>" << SpecUtils::printCompact( line.m_yield, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( info.m_counts, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( info.m_counts_uncert, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( counts_over_yield, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( counts_uncert_percent, nsigfig ) << "%"
      << "</td><td>" << SpecUtils::printCompact( info.m_base_rel_eff_uncert, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( meas_rel_eff, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( meas_rel_eff_uncert, nsigfig ) << "%";

      if( any_peak_has_multiple_srcs )
      {
        results_html << "</td><td>";
        if( info.m_source_gammas.size() >= 2 )
        {
          double sum_counts = 0.0;
          for( const GenericLineInfo &sum_line : info.m_source_gammas )
            sum_counts += relative_activity(sum_line.m_isotope) * sum_line.m_yield;
          const double contrib_percent = 100.0*(rel_act*line.m_yield) / sum_counts;
          results_html << SpecUtils::printCompact( contrib_percent, nsigfig ) << "%";
        }
      }//if( any_peak_has_multiple_srcs )

      if( has_decay_corr )
      {
        results_html << "</td><td>";
        
        for( const GenericPeakInfo &un_corr_peak : m_input.peaks_before_decay_correction )
        {
          if( fabs(info.m_energy - un_corr_peak.m_energy) > 0.001 )
            continue;
          
          for( const GenericLineInfo &un_corr_line : un_corr_peak.m_source_gammas )
          {
            if( un_corr_line.m_isotope != line.m_isotope )
              continue;
            
            const double ratio = line.m_yield / un_corr_line.m_yield;
            if( !IsInf(ratio) && !IsNan(ratio) )
              results_html << " " << SpecUtils::printCompact( ratio, nsigfig );
            else
              results_html << "--";
          }//for( loop over un_corr_peak.m_source_gammas )
        }//for( loop over m_input.peaks_before_decay_correction )
      }//if( has_decay_corr )
      
      results_html << "</td></tr>\n";
    }//for( size_t i = 0; i < info.m_source_gammas.size(); ++i )
  }//for( const RelEff::PeakInfo &info : input_peaks )
  
  results_html << "  </tbody>\n"
  << "</table>\n\n";
  
  
  auto html_sanitize = []( string &val ){
    // We'll do some really stupid/simple sanitization
    SpecUtils::ireplace_all( val, "&", "&amp;" );
    SpecUtils::ireplace_all( val, "<", "&lt;" );
    SpecUtils::ireplace_all( val, ">", "&gt;" );
    SpecUtils::ireplace_all( val, "'", "&#39;" );
    SpecUtils::ireplace_all( val, "\"", "&quot;" );
  };
  
  results_html << "</div>\n";
  
  
  if( !m_warnings.empty() )
  {
    results_html << "<div class=\"warnings\">\n"
    << "<h3>Warnings</h3>\n";
    for( string warning : m_warnings )
    {
      html_sanitize( warning );
      
      results_html << "<div class=\"warningline\">" << warning << "</div>\n";
    }//for( string warning : warnings )
    
    results_html << "</div>\n";
  }//if( !warnings.empty() )
  
  
  time_t rawtime;
  struct tm *timeinfo;
  time( &rawtime );
  timeinfo = localtime (&rawtime);
  results_html << "<div class=\"anatime\">Analysis performed " << asctime(timeinfo)
  << " with " << spectrum_title << " compiled " __TIMESTAMP__
  << "</div>\n";
  
  
  //Write out the data JSON and CSS data
  stringstream rel_eff_plot_values;
  
  // Previous to 20250520, we used to use CSS to color data markers.
  //  Then we switched to just doing it in the JSON/JS directly, to make it
  //  easier to keep things in sync acros doing HTML reports and within interactive
  //  InterSpec - leaving the code in, but commented out for for the moment incase we
  //  want to go back
  //stringstream add_rel_eff_plot_css;
  
  rel_eff_plot_values << "[";
  for( size_t index = 0; index < m_input.peaks.size(); ++index )
  {
    const GenericPeakInfo &peak = m_input.peaks[index];
    
    string isotopes_json;
    double src_counts = 0.0;
    for( const GenericLineInfo &line : peak.m_source_gammas )
    {
      const double rel_act = relative_activity( line.m_isotope );
      src_counts += rel_act * line.m_yield;
      
      //const double meas_rel_eff = info.m_counts / (info.m_source_gammas[i].m_yield * rel_act);
      
      snprintf( buffer, sizeof(buffer), "%s{\"nuc\": \"%s\", \"br\": %1.6G, \"rel_act\": %1.6G}",
               (isotopes_json.empty() ? "" : ", "), line.m_isotope.c_str(), line.m_yield, rel_act );
      
      isotopes_json += buffer;
    }//for( const RelEff::GammaLineInfo &line : peak.m_source_gammas )
    
    const double eff = peak.m_counts / src_counts;
    const double eff_uncert = peak.m_counts_uncert / src_counts;
    
    snprintf( buffer, sizeof(buffer),
             "%s"
             "{\"energy\": %.2f, \"counts\": %1.7g,"
             " \"counts_uncert\": %1.7g, \"eff\": %1.6g,"
             " \"eff_uncert\": %1.6g",
             (index ? ", " : ""), 
             peak.m_mean, peak.m_counts,
             peak.m_counts_uncert, eff,
             eff_uncert );
    
    rel_eff_plot_values << buffer;
    rel_eff_plot_values << ", \"nuc_info\": [" << isotopes_json << "]}";
  }//for( size_t index = 0; index < input_peaks.size(); ++index )
  
  rel_eff_plot_values << "]";
  
  
  set<const SandiaDecay::Nuclide *> nuclides_with_colors;
  for( const shared_ptr<const PeakDef> &p : display_peaks )
  {
    assert( p );
    const SandiaDecay::Nuclide * const nuc = p ? p->parentNuclide() : nullptr;
    if( nuc && !nuclides_with_colors.count(nuc) && !p->lineColor().isDefault() )
    {
      //add_rel_eff_plot_css << "        .RelEffPlot circle." << nuc->symbol
      //                     << "{ fill: " << p->lineColor().cssText() << "; }\n";
      nuclides_with_colors.insert( nuc );
    }
  }//for( const shared_ptr<const PeakDef> &p : display_peaks )
  
  size_t unseen_nuc_index = 0;
  const vector<string> default_nuc_colors{ "#003f5c", "#ffa600", "#7a5195", "#ff764a", "#ef5675", "#374c80" };
  for( const IsotopeRelativeActivity &act : m_rel_activities )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( act.m_isotope );
    if( !nuclides_with_colors.count(nuc) )
    {
      string nucstr = act.m_isotope;
      if( nucstr.empty() )
        nucstr = "default";
      
      for( size_t i = 0; i < nucstr.size(); ++i )
        nucstr[i] = std::isalpha( static_cast<unsigned char>(nucstr[i]) ) ? nucstr[i] : '_';
      
      //add_rel_eff_plot_css << "        .RelEffPlot circle." << nucstr << "{ fill: "
      //<< default_nuc_colors[unseen_nuc_index % default_nuc_colors.size()]
      //<< "; }\n";
      
      unseen_nuc_index += 1;
    }
  }//for( const IsotopeRelativeActivity &act : m_rel_activities )
  
  
  auto load_file_contents = []( string filename ) -> string {
    Wt::WApplication *app = Wt::WApplication::instance();
    string filepath = app ? app->docRoot() : string("");
    filepath = SpecUtils::append_path( filepath, "InterSpec_resources" );
    filepath = SpecUtils::append_path( filepath, filename );
    
    vector<char> file_data;
    try
    {
      SpecUtils::load_file_data( filepath.c_str(), file_data );
    }catch( std::exception & )
    {
      throw std::runtime_error( "Failed to read " + filename );
    }
    
    return string( begin( file_data ), end( file_data ) );
  };//load_file_contents(...)
  
  
  string html = load_file_contents( "static_text/manual_rel_act_report.tmplt.html" );
  const string d3_js = load_file_contents( "d3.v3.min.js" );
  
  SpecUtils::ireplace_all( html, "\\;", ";" );
  
  SpecUtils::ireplace_all( html, "${D3_SCRIPT}", d3_js.c_str() );
  
  SpecUtils::ireplace_all( html, "${TITLE}", spectrum_title.c_str() );
  
  SpecUtils::ireplace_all( html, "${REL_EFF_DATA_VALS}", rel_eff_plot_values.str().c_str() );
  
  string rel_eff_fcn;
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    rel_eff_fcn = RelActCalc::rel_eff_eqn_js_function( m_input.eqn_form, m_rel_eff_eqn_coefficients );
  }else
  {
    const ManualGenericRelActFunctor::PhysModelRelEqnDef input
           = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
    
    rel_eff_fcn = RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                                                                     input.external_attens,
                                                                     input.det.get(),
                                                                     input.hoerl_b,
                                                                     input.hoerl_c );
  }

  string unc_fcn = rel_eff_eqn_js_uncert_fcn();

  SpecUtils::ireplace_all( html, "${FIT_REL_EFF_EQUATION}", rel_eff_fcn.c_str() );
  SpecUtils::ireplace_all( html, "${FIT_REL_EFF_EQUATION_UNCERTAINTY}", unc_fcn.c_str() );
  SpecUtils::ireplace_all( html, "${RESULTS_TXT}", results_html.str().c_str() );
  
  
  const string rel_eff_plot_js = load_file_contents( "RelEffPlot.js" );
  const string rel_eff_plot_css = load_file_contents( "RelEffPlot.css" );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_JS}", rel_eff_plot_js.c_str() );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_CSS}", rel_eff_plot_css.c_str() );
  //SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_ADDITIONAL_CSS}", add_rel_eff_plot_css.str().c_str() );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_ADDITIONAL_CSS}", "" );
  
  
  if( spectrum )
  {
    stringstream set_js_str;
    
    D3SpectrumExport::write_js_for_chart( set_js_str, "specchart", "", "Energy (keV)", "Counts"  );
    
    set_js_str <<
    "  let spec_observer = new ResizeObserver(entries => {\n"
    "    for (let entry of entries) {\n"
    "      if (entry.target && (entry.target.id === \"specchart\")) {\n"
    "        spec_chart_specchart.handleResize(false);\n"
    "      }\n"
    "    }\n"
    "  });\n"
    "  spec_observer.observe( document.getElementById(\"specchart\") );\n"
    ;
    
    D3SpectrumExport::D3SpectrumChartOptions chart_options;
    chart_options.m_useLogYAxis = true;
    chart_options.m_legendEnabled = false;
    chart_options.m_compactXAxis = true;
    chart_options.m_allowDragRoiExtent = false;
    write_set_options_for_chart( set_js_str, "specchart", chart_options );
    set_js_str << "  spec_chart_specchart.setShowLegend(false);\n";
    set_js_str << "  spec_chart_specchart.setXAxisRange("
               << spectrum->gamma_energy_min()
               << ", " << spectrum->gamma_energy_max() << ", false);\n";
    
    D3SpectrumExport::D3SpectrumOptions spec_options;
    spec_options.spectrum_type = SpecUtils::SpectrumType::Foreground;
    
    vector<shared_ptr<const PeakDef>> peaks;
    for( const auto &p : display_peaks )
      peaks.push_back( make_shared<PeakDef>(*p) );
    spec_options.peaks_json = PeakDef::peak_json( peaks, spectrum, Wt::WColor(0,51,255), 255 );
    
    
    vector<pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > meas_to_plot;
    
    const SpecUtils::Measurement * const meas_ptr = spectrum.get();
    meas_to_plot.emplace_back(meas_ptr, spec_options);
    
    if( background )
    {
      D3SpectrumExport::D3SpectrumOptions spec_options;
      spec_options.spectrum_type = SpecUtils::SpectrumType::Background;
      spec_options.line_color = "steelblue";
      spec_options.display_scale_factor = 1.0;
      if( !background_peaks.empty() )
        spec_options.peaks_json = PeakDef::peak_json( background_peaks, background, Wt::WColor(0,51,255), 255 );


      if( (background_normalization <= 0.0f)
         || IsNan(background_normalization)
         || IsInf(background_normalization) )
      {
        float back_lt = background->live_time();
        if( (back_lt <= 0.0f) || IsNan(back_lt) || IsNan(back_lt) )
          back_lt = background->real_time();
        
        float fore_lt = spectrum->live_time();
        if( (fore_lt <= 0.0f) || IsNan(fore_lt) || IsNan(fore_lt) )
          fore_lt = spectrum->real_time();
        
        background_normalization = fore_lt / back_lt;
      }//if( background normalization wasnt provided )
      
      if( (background_normalization > 0.0f)
         && !IsNan(background_normalization)
         && !IsInf(background_normalization) )
      {
        spec_options.display_scale_factor = background_normalization;
      }
        
      spec_options.title = "Background";
      meas_to_plot.emplace_back(background.get(), spec_options);
    }//if( background )
    
    
    D3SpectrumExport::write_and_set_data_for_chart( set_js_str, "specchart", meas_to_plot );
    
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_DIV}",
                            "<div id=\"specchart\" style=\"height: 30vw; flex: 1 2; overflow: hidden;\" class=\"SpecChart\"></div>" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_INIT_JS}", set_js_str.str().c_str() );
    
    
    const string spectrum_chart_d3_js = load_file_contents( "SpectrumChartD3.js" );
    const string spectrum_chart_d3_css = load_file_contents( "SpectrumChartD3.css" );
    
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_JS}", spectrum_chart_d3_js.c_str() );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_CSS}", spectrum_chart_d3_css.c_str() );
    
    SpecUtils::ireplace_all( html, "${CHART_SPACER_LEFT}", "" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_RIGHT}", "" );
  }else
  {
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_DIV}", "" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_JS}", "" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_CSS}", "" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_INIT_JS}", "" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_LEFT}", "<div style=\"width: 10%\"> </div>" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_RIGHT}", "<div style=\"width: 15%\"> </div>" );
  }//if( spectrum ) / else
  
  
  output_html_file << html;
}//void print_html_report( std::ostream &strm ) const


RelEffSolution solve_relative_efficiency( const RelEffInput &input_orig )
{
  // When fitting the AN using the Physical Model, we can easily get caught in a local-minimum, and also
  // we dont currently have great control over the step sizes for AN/AD, so here is 
  // a work-around, to get a decent starting point for AN.
  // Right now we are scanning on AN, but scanning on AD as well would be better.
  // I assume that we could avoid this with a proper implementation of DynamicCostFunction.
#define SCAN_AN_FOR_BEST_FIT 1

  // A mutable copy: the AN scan replaces the self-attenuator with its best-fit one, and the
  //  profile-target resolution below fills in element rosters.
  RelEffInput input = input_orig;

  input.check_nuclide_constraints();

  // Resolve the profile targets (see RelEffInput::profile_targets) against the nuclides actually
  //  in the problem, before the cost functor - which allocates the penalty residual channel for
  //  them - is built.  A target we cannot make sense of is dropped with a warning rather than
  //  failing the whole fit; the covariance-based uncertainties remain as the fallback.
  vector<string> profile_setup_warnings;
  if( !input.profile_targets.empty() )
  {
    set<string> problem_isotopes;
    for( const GenericPeakInfo &peak : input.peaks )
      for( const GenericLineInfo &line : peak.m_source_gammas )
        problem_isotopes.insert( line.m_isotope );

    // `database()` THROWS on failure rather than returning null, and this runs outside the
    //  try/catch that turns problems into a solution status - so catch it here and let the
    //  per-target `db` guards below drop what they cannot resolve.
    const SandiaDecay::SandiaDecayDataBase *db = nullptr;
    try
    {
      db = DecayDataBaseServer::database();
    }catch( std::exception & )
    {
      db = nullptr;
    }

    vector<ProfileTarget> resolved;
    for( ProfileTarget target : input.profile_targets )
    {
      const string what = "profile of " + target.m_nuclide;

      if( !problem_isotopes.count(target.m_nuclide) )
      {
        profile_setup_warnings.push_back( "Skipped the " + what + ": it is not in the fit." );
        continue;
      }

      if( target.m_type == ProfileTarget::Type::RelativeActivity )
      {
        resolved.push_back( std::move(target) );
        continue;
      }//if( RelativeActivity )

      if( target.m_type == ProfileTarget::Type::ActivityRatio )
      {
        if( !problem_isotopes.count(target.m_denom_nuclide) )
        {
          profile_setup_warnings.push_back( "Skipped the " + what + "/" + target.m_denom_nuclide
                                            + " activity ratio: the denominator is not in the fit." );
          continue;
        }

        if( target.m_denom_nuclide == target.m_nuclide )
        {
          profile_setup_warnings.push_back( "Skipped the " + what + " activity ratio: the"
                                            " numerator and denominator are the same nuclide." );
          continue;
        }

        resolved.push_back( std::move(target) );
        continue;
      }//if( ActivityRatio )

      // A mass fraction pinned by a fixed constraint cannot move, so there is nothing to profile.
#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
      bool is_pinned = false;
      for( const MassFractionConstraint &mfc : input.mass_fraction_constraints )
        is_pinned |= ((mfc.m_nuclide == target.m_nuclide)
                      && (mfc.m_mass_fraction_lower == mfc.m_mass_fraction_upper));
      if( is_pinned )
      {
        profile_setup_warnings.push_back( "Skipped the " + what + " mass fraction: it is fixed by"
                                          " a mass-fraction constraint." );
        continue;
      }

      // A constraint on the target already names the element roster the fraction is relative to.
      if( target.m_specific_activities.empty() )
      {
        for( const MassFractionConstraint &mfc : input.mass_fraction_constraints )
          if( mfc.m_nuclide == target.m_nuclide )
            target.m_specific_activities = mfc.m_specific_activities;
      }
#endif //USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT

      const SandiaDecay::Nuclide * const target_nuc = db ? db->nuclide( target.m_nuclide ) : nullptr;
      if( target.m_specific_activities.empty() && target_nuc )
      {
        for( const string &iso : problem_isotopes )
        {
          const SandiaDecay::Nuclide * const nuc = db->nuclide( iso );
          if( nuc && (nuc->atomicNumber == target_nuc->atomicNumber) )
            target.m_specific_activities[iso] = nuc->activityPerGram();
        }
      }//if( no roster from the caller or a constraint )

      // The roster must be a subset of the fitted nuclides: an isotope that is not fit has no
      //  activity to contribute, and asking the functor for one would throw.
      for( map<string,double>::iterator iter = begin(target.m_specific_activities);
           iter != end(target.m_specific_activities); /**/ )
      {
        if( problem_isotopes.count(iter->first) && (iter->second > 0.0) )
          ++iter;
        else
          iter = target.m_specific_activities.erase( iter );
      }

      if( target.m_specific_activities.size() < 2 )
      {
        profile_setup_warnings.push_back( "Skipped the " + what + " mass fraction: fewer than two"
                    " isotopes of its element are in the fit, so the fraction is 1 by definition." );
        continue;
      }

      if( !target.m_specific_activities.count(target.m_nuclide) )
      {
        profile_setup_warnings.push_back( "Skipped the " + what + " mass fraction: it is not in"
                                          " the roster the fraction would be relative to." );
        continue;
      }

      resolved.push_back( std::move(target) );
    }//for( ProfileTarget target : input.profile_targets )

    input.profile_targets.swap( resolved );
  }//if( !input.profile_targets.empty() )

  // Estimate the additional per-peak fractional uncertainty from the data, if asked to.  The
  //  estimate needs a fitted model to measure deviations against, so: solve, estimate, apply,
  //  and repeat (the weights change the fit a little, which changes the estimate a little).
  //  Converges in a couple of rounds; sub-solves skip the covariance work.
  // Widen the peak uncertainties by one common factor estimated from the data, if asked to.
  //  The estimate needs a fitted model to measure deviations against, so solve once first.  A
  //  single pass suffices: scaling every peak by the same factor cannot move the fit, so
  //  re-estimating afterwards would only ever return 1.
  double auto_stat_multiple = -1.0;
  if( input.auto_estimate_add_uncert )
  {
    RelEffInput trial_input = input;
    trial_input.auto_estimate_add_uncert = false;
    trial_input.point_estimate_only = true;
    trial_input.profile_targets.clear(); //only the final solve profiles

    const RelEffSolution trial_solution = solve_relative_efficiency( trial_input );
    if( trial_solution.m_status == ManualSolutionStatus::Success )
      auto_stat_multiple = estimate_stat_uncert_multiple( trial_solution );

    if( auto_stat_multiple >= 1.0 )
    {
      for( GenericPeakInfo &peak : input.peaks )
        peak.m_counts_uncert *= auto_stat_multiple;
    }
  }//if( input.auto_estimate_add_uncert )

  const std::vector<GenericPeakInfo> &peak_infos = input.peaks;
  const RelActCalc::RelEffEqnForm eqn_form = input.eqn_form;
  const size_t eqn_order = input.eqn_order;

  const auto start_time = std::chrono::high_resolution_clock::now();
  
  RelEffSolution solution;
  
  DoWorkOnDestruct setFinalTime( [&solution,start_time](){
    const auto end_time = std::chrono::high_resolution_clock::now();
    solution.m_num_microseconds_eval = std::chrono::duration<double, std::micro>(end_time - start_time).count();
  });
  
  solution.m_input = input;
  solution.m_status = ManualSolutionStatus::NotInitialized;

  solution.m_warnings.insert( begin(solution.m_warnings),
                      begin(input.prep_warnings), end(input.prep_warnings) );

  solution.m_warnings.insert( end(solution.m_warnings),
                      begin(profile_setup_warnings), end(profile_setup_warnings) );

  if( input.auto_estimate_add_uncert )
  {
    solution.m_auto_stat_uncert_multiple = auto_stat_multiple;
    if( auto_stat_multiple < 0.0 )
      solution.m_warnings.push_back( "Could not estimate how much to widen the peak uncertainties"
                                     " from the data; using them as they are." );
    else if( auto_stat_multiple >= 24.9 ) //i.e. it ran into estimate_stat_uncert_multiple's cap
      solution.m_warnings.push_back( "The peaks disagree with any achievable relative efficiency"
                " curve so badly that the estimated uncertainty widening hit its limit - treat this"
                " fit, and its uncertainties, with suspicion." );
  }//if( input.auto_estimate_add_uncert )

#if( SCAN_AN_FOR_BEST_FIT )
  const bool scan_an_for_best_fit = (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel
                                      && input.phys_model_self_atten
                                      && input.phys_model_self_atten->fit_atomic_number
                                      && !input.skip_an_scan );
  if( scan_an_for_best_fit )
  {
    SpecUtilsAsync::ThreadPool threadpool;
    std::mutex best_chi2_mutex;
    double best_chi2 = std::numeric_limits<double>::max();
    double best_an = 0.0, best_ad = 0.0;
    
    const auto do_work = [&input, &best_chi2, &best_an, &best_ad, &best_chi2_mutex](
                                        const double an, const double start_ad_gcm2 ){
      RelEffInput new_input = input;
      auto new_self_atten = std::make_shared<RelActCalc::PhysicalModelShieldInput>(*input.phys_model_self_atten);
      new_self_atten->fit_atomic_number = false;
      new_self_atten->atomic_number = an;
      if( (start_ad_gcm2 > 0.0) && new_self_atten->fit_areal_density )
      {
        // Only use the multi-start AD if it is inside the callers allowed fit range (Ceres
        //  requires the initial value within bounds).
        double lower_ad = new_self_atten->lower_fit_areal_density / PhysicalUnits::g_per_cm2;
        double upper_ad = new_self_atten->upper_fit_areal_density / PhysicalUnits::g_per_cm2;
        if( (lower_ad == 0.0) && (upper_ad == 0.0) )
          upper_ad = RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2;
        if( (start_ad_gcm2 >= lower_ad) && (start_ad_gcm2 <= upper_ad) )
          new_self_atten->areal_density = start_ad_gcm2 * PhysicalUnits::g_per_cm2;
      }
      new_input.phys_model_self_atten = new_self_atten;
      new_input.point_estimate_only = true; //scan only reads m_chi2 and the fitted AD - skip covariance work
      new_input.profile_targets.clear(); //only the final solve profiles
      // The peaks already carry any auto-estimated additional uncertainty (it is applied before
      //  this scan runs); re-estimating it per scan point would both cost several solves each and
      //  make the scan compare chi2 values computed with different weights.
      new_input.auto_estimate_add_uncert = false;
      RelEffSolution solution = solve_relative_efficiency( new_input );

      std::lock_guard<std::mutex> lock(best_chi2_mutex);

      if( (solution.m_status==ManualSolutionStatus::Success) && (solution.m_chi2 < best_chi2) )
      {
        best_chi2 = solution.m_chi2;
        best_an = an;
        assert( solution.m_phys_model_self_atten_shield );
        best_ad = solution.m_phys_model_self_atten_shield->m_areal_density / PhysicalUnits::g_per_cm2;
      }
    };//do_work

    double an_step = 5.0;
    double min_an = input.phys_model_self_atten->lower_fit_atomic_number;
    double max_an = input.phys_model_self_atten->upper_fit_atomic_number;

    // A (0,0) fit range means "use the default range" - the same convention the bounds-setting
    //  lambda applies; without this the scan would degenerate to the single point AN = 0.
    if( (min_an == 0.0) && (max_an == 0.0) )
    {
      min_an = 1.0;
      max_an = 98.0;
    }

    // AD is the other known local-minimum axis, so the coarse pass can additionally be
    //  multi-started over areal-density decades (add e.g. 0.5, 5.0, 50.0 below).  Measured on
    //  the spec184 problem with a fit-AN self attenuator (20260816): the decade multi-start
    //  tripled the scan cost (3.7 -> 11.1 CPU-seconds) and reached the identical chi2 and
    //  enrichment, so only the caller's own starting AD is used until a problem shows a benefit.
    //  (A start of 0 keeps the caller's input AD.)
    const double coarse_start_ads[] = { 0.0 }; //g/cm2
    for( double an = min_an; an <= max_an; an += an_step )
    {
      for( const double start_ad : coarse_start_ads )
        threadpool.post( [&do_work, an, start_ad](){ do_work(an, start_ad); } );
    }
    threadpool.join();

    if( best_chi2 != std::numeric_limits<double>::max() )
    {
      if( debug_printout() )
        cout << "Initial best AN = " << best_an << " AD = " << best_ad << endl;
      const double refine_start_ad = best_ad;
      best_chi2 = std::numeric_limits<double>::max();
      min_an = std::max( min_an, best_an - an_step );
      max_an = std::min( max_an, best_an + an_step );
      an_step = 1.0;
      for( double an = min_an; an <= max_an; an += an_step )
        threadpool.post( [&do_work, an, refine_start_ad](){ do_work(an, refine_start_ad); } );
      threadpool.join();
    }
    
    if( best_chi2 != std::numeric_limits<double>::max() )
    {
      if( debug_printout() )
        cout << "Final best AN = " << best_an << " AD = " << best_ad << endl;
      auto new_self_atten = std::make_shared<RelActCalc::PhysicalModelShieldInput>(*input.phys_model_self_atten);
      //new_self_atten->fit_atomic_number = input.use_ceres_to_fit_eqn ? false : true;
      new_self_atten->atomic_number = best_an;
      new_self_atten->areal_density = best_ad * PhysicalUnits::g_per_cm2;
      input.phys_model_self_atten = new_self_atten;
      
      // // How many parameters, and how they are defined will change based on if we are fitting
      // //  atomic number - so we need to update the input stored in the solution to reflect what
      // //  we actually used.
      solution.m_input = input;
    }//if( best_chi2 != std::numeric_limits<double>::max() )
  }//if( scan_an_for_best_fit )
#endif //SCAN_AN_FOR_BEST_FIT

  try
  {
    if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {  
      if( !input.phys_model_detector || !input.phys_model_detector->isValid() )
        throw runtime_error( "You must specify a detector for the Physical Model." );
      
      // Count the number of parameters being fit to ensure we have enough peaks
      // Count unique nuclides in peaks
      std::set<std::string> unique_nuclides;
      for( const GenericPeakInfo &peak : peak_infos )
      {
        for( const GenericLineInfo &line : peak.m_source_gammas )
          unique_nuclides.insert( line.m_isotope );
      }

      size_t num_fit_parameters = unique_nuclides.size();

      // Count self-attenuation parameters
      if( input.phys_model_self_atten )
      {
        if( input.phys_model_self_atten->fit_areal_density )
          num_fit_parameters += 1; // atomic number
        if( !input.phys_model_self_atten->material && input.phys_model_self_atten->fit_atomic_number )
          num_fit_parameters += 1; 
      }

      // Count external attenuation parameters
      for( const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &ext_atten : input.phys_model_external_attens )
      {
        if( ext_atten )
        {
          if( ext_atten->fit_areal_density )
            num_fit_parameters += 1; // atomic number
          if( !ext_atten->material && ext_atten->fit_atomic_number )
            num_fit_parameters += 1;
        }
      }

      // Count Hoerl equation parameters (b and c)
      if( input.phys_model_use_hoerl )
        num_fit_parameters += 2;

      // Ensure we have at least as many peaks as parameters being fit
      if( peak_infos.size() < num_fit_parameters )
      {
        string msg = "You must specify at least " + std::to_string(num_fit_parameters) + " peak"
          + string(num_fit_parameters > 1 ? "s" : "")
          + " for the Physical Model (currently have " + std::to_string(peak_infos.size())
          + "), since you are fitting " + std::to_string(num_fit_parameters) + " parameter"
          + string(num_fit_parameters > 1 ? "s" : "")
          + ".";
        throw runtime_error( msg );
      }

      if( !input.use_ceres_to_fit_eqn )
        throw logic_error( "You must specify to use Ceres to fit Rel. Eff. equation for the Physical Model." );
      
      //if( input.use_ceres_to_fit_eqn && input.phys_model_self_atten && input.phys_model_self_atten->fit_atomic_number )
      //  throw logic_error( "You can not use Ceres to fit atomic number of Physical relative efficiency equation." );
      
      if( input.phys_model_self_atten && !input.phys_model_self_atten->material
         && ((input.phys_model_self_atten->atomic_number < 0.999) || (input.phys_model_self_atten->atomic_number > 98.001)) )
        throw logic_error( "A self-attenuator is specified, but with no material or atomic number." );
      
      // TODO: add checks for external attenuators
    }else
    {
      if( input.phys_model_self_atten || !input.phys_model_external_attens.empty() )
        throw runtime_error( "Attenuations can only be specified for FramPhysicalModel." );
    }//if( RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorInitializing;
    solution.m_error_message = e.what();
    return solution;
  }//try / catch to check input


  ManualGenericRelActFunctor *cost_functor = nullptr;
  try
  {
    cost_functor = new ManualGenericRelActFunctor( input );
    solution.m_status = ManualSolutionStatus::ErrorFindingSolution;
    
    solution.m_warnings.insert( end(solution.m_warnings),
                               begin(cost_functor->m_setup_warnings),
                               end(cost_functor->m_setup_warnings) );
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorInitializing;
    solution.m_error_message = e.what();
    
    return solution;
  }//try / catch to setup the problem
  
  const size_t num_peaks = cost_functor->m_input.peaks.size();
  const size_t num_nuclides = cost_functor->m_isotopes.size();
  const size_t num_parameters = cost_functor->num_parameters();

  // The chi2 of the DATA channels alone, at the given parameters.  Everything that reports or
  //  differences a chi2 goes through here, so an armed profile-likelihood penalty channel (see
  //  RelEffInput::profile_targets) can never leak into a reported number: 2*Ceres' final cost
  //  would include it.  With the channel disarmed this is exactly 2*final_cost.
  const auto peaks_only_chi2 = [cost_functor, num_peaks]( const vector<double> &pars ) -> double {
    vector<double> residuals( cost_functor->number_residuals(), 0.0 );
    cost_functor->eval( pars, residuals.data() );

    double chi2 = 0.0;
    for( size_t i = 0; i < num_peaks; ++i )
      chi2 += residuals[i]*residuals[i];
    return chi2;
  };//peaks_only_chi2 lambda

  solution.m_activity_norms = cost_functor->m_rel_act_norms;

  // Relative activities multiples start out as 1.0 because ManualGenericRelActFunctor constructor
  //   estimates the activities for a flat rel eff = 1.0; see
  //   #ManualGenericRelActFunctor::m_rel_act_norms.
  vector<double> parameters( num_parameters, 1.0 );
  double *pars = &parameters[0];


  ceres::CostFunction *cost_function = nullptr;

  // From a few example cases inspected by hand, it looks like auto and numerical differentiation
  //  get the same answers, but auto diff is faster, and requires a lot fewer evaluations, so
  //  we'll just always use auto differentiation.
  const bool use_auto_diff = true;
  const int auto_diff_stride = 8;
  if( use_auto_diff )
  {
    ceres::DynamicAutoDiffCostFunction<ManualGenericRelActFunctor,auto_diff_stride> *dyn_auto_diff_cost_function
          = new ceres::DynamicAutoDiffCostFunction<ManualGenericRelActFunctor,auto_diff_stride>( cost_functor,
                                                                          ceres::TAKE_OWNERSHIP );
    // The number of residuals is the number of peaks, unless USE_RESIDUAL_TO_BREAK_DEGENERACY then
    //  we add one more residual to clamp the relative efficiency curve to 1.0 at the lowest energy.
 
    dyn_auto_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
    dyn_auto_diff_cost_function->AddParameterBlock( static_cast<int>(num_parameters) );

    cost_function = dyn_auto_diff_cost_function;
  }else
  {
    ceres::NumericDiffOptions num_diff_options;

#if( SCAN_AN_FOR_BEST_FIT )
    if( scan_an_for_best_fit )
    {
      // TODO: need to more closely evaluate the step sizes
      num_diff_options.ridders_relative_initial_step_size = 0.1;
      num_diff_options.relative_step_size = 1.0E-3;
    }
#endif //SCAN_AN_FOR_BEST_FIT

    ceres::DynamicNumericDiffCostFunction<ManualGenericRelActFunctor> *dyn_num_diff_cost_function
          = new ceres::DynamicNumericDiffCostFunction<ManualGenericRelActFunctor>( cost_functor,
                                                        ceres::TAKE_OWNERSHIP, num_diff_options );
    dyn_num_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );        
    dyn_num_diff_cost_function->AddParameterBlock( static_cast<int>(num_parameters) );     

    cost_function = dyn_num_diff_cost_function;
  }//if( use_auto_diff ) / else
    
  
  ceres::Problem problem;
  
  // Robust/outlier-tolerant fitting has been investigated repeatedly (most recently 20260816,
  //  in the detail below) and is deliberately NOT adopted.  Earlier passes reached the same
  //  conclusion from the other direction: the stock Huber and Cauchy losses change little on any
  //  single problem, which the mechanism in (1) explains.  Do not re-open without new evidence.
  //
  //  1) A ceres::LossFunction cannot do it here.  Ceres applies the loss to the squared norm of a
  //     whole RESIDUAL BLOCK, and this problem adds one block holding every peak, so rho() acts
  //     on the total chi2 - a monotone transform of the cost, which leaves the minimum exactly
  //     where it was.  Measured: HuberLoss(3.5) moved the enrichment by 1e-4 percentage points
  //     (i.e., nothing), while perturbing the reported uncertainties.
  //  2) Genuine per-peak robustification (Huber influence applied to each residual inside the
  //     cost functor, so peaks past N sigma contribute linearly) was implemented and measured at
  //     N = 4, 3, 2, and 1.5 sigma.  It moved the enrichment by at most 0.08 percentage points -
  //     and slightly AWAY from the reference value - while shrinking the reported uncertainties,
  //     which is the wrong direction: those outlying peaks reflect real peak-fit/model error that
  //     belongs in the error budget.
  //
  //  Spectra used: `manual_rel_eff/spec184_235U_12.9543.n42` (IAEA IDB reference spectrum,
  //  12.9543 wt% U235, HPGe, 13.2M counts / 693 s live time, 24 fit peaks spanning 72-1001 keV),
  //  fit both with the Physical Model and with an empirical LnX order-4 curve.  In both cases the
  //  chi2 excess comes from a handful of individually mis-fit peaks; that is handled where it
  //  belongs - in the uncertainty scale (see the outlier-insensitive `m_cov_scale` below) and by
  //  telling the user which peaks to re-fit - rather than by down-weighting data in the fit.
  ceres::LossFunction *lossfcn = nullptr;
  //ceres::LossFunction *lossfcn = new ceres::CauchyLoss(1.0);

  problem.AddResidualBlock( cost_function, lossfcn, pars );
  problem.AddParameterBlock( pars, static_cast<int>(num_parameters) );

  vector<int> constant_parameters;

  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    assert( input.use_ceres_to_fit_eqn );
    
    size_t par_num = num_nuclides;
    assert( par_num <= num_parameters );
    
    const auto &self_atten_opt = input.phys_model_self_atten;
    setup_physical_model_shield_par_manual( constant_parameters, pars, par_num, self_atten_opt );
    
    assert( par_num <= num_parameters );
    
    for( size_t ext_ind = 0; ext_ind < input.phys_model_external_attens.size(); ++ext_ind )
    {
      assert( par_num <= num_parameters );
      const auto &opt = input.phys_model_external_attens[ext_ind];
      setup_physical_model_shield_par_manual( constant_parameters, pars, par_num, opt );
    }//for( size_t ext_ind = 0; ext_ind < options.phys_model_external_atten.size(); ++ext_ind )

    assert( par_num <= num_parameters );
    
    if( input.phys_model_use_hoerl )
    {
      // set the b and c parameters for the relative efficiency equation
      pars[par_num] = (0.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset;  //(energy/1000)^b - start b at 0, so term is 1.0
      par_num += 1;
      pars[par_num] = (1.0/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset;  //c^(1000/energy) - start c at 1, so term is 1
      par_num += 1;
    }//if( input.phys_model_use_hoerl )
    
    assert( par_num == num_parameters );
    if( par_num != num_parameters )
      throw logic_error( "Num paramaters doesnt match expected" );
  }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  
  for( const ManualActRatioConstraint &constraint : input.act_ratio_constraints )
  {
    assert( num_nuclides == cost_functor->m_isotopes.size() );

    for( int i = 0; i < static_cast<int>(num_nuclides); ++i )
    {
      if( constraint.m_constrained_nuclide == cost_functor->m_isotopes[i] )
      {
        assert( std::find( begin(constant_parameters), end(constant_parameters), i ) == end(constant_parameters) );
        constant_parameters.push_back( i );
        pars[i] = -1.0; //so we can assert on this later to make sure things are reasonable
        break;
      }
    }//for( int i = 0; i < static_cast<int>(num_nuclides); ++i )
  }//for( const ManualActRatioConstraint &constraint : input.act_ratio_constraints )

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
  for( const MassFractionConstraint &constraint : input.mass_fraction_constraints )
  {
    assert( num_nuclides == cost_functor->m_isotopes.size() );

    if( constraint.m_mass_fraction_lower == constraint.m_mass_fraction_upper )
    {
      for( int i = 0; i < static_cast<int>(num_nuclides); ++i )
      {
        if( constraint.m_nuclide == cost_functor->m_isotopes[i] )
        {
          // Fixed constraints decode to their pinned fraction, so the slot is constant - EXCEPT
          //  when it is an all-constrained elements carrier, which holds the element scale (a
          //  live parameter; see #ManualMassFracBlock) - configured in the bounds loop below.
          bool is_all_const_carrier = false;
          for( const auto &block : cost_functor->m_mass_frac_blocks )
            is_all_const_carrier |= (block.spec.all_constrained && (block.carrier_index == static_cast<size_t>(i)));

          if( !is_all_const_carrier )
          {
            assert( std::find( begin(constant_parameters), end(constant_parameters), i ) == end(constant_parameters) );
            constant_parameters.push_back( i );
            pars[i] = -1.0; //so we can assert on this later to make sure things are reasonable
          }
          break;
        }
      }//for( int i = 0; i < static_cast<int>(num_nuclides); ++i )
    }else
    {

    }//if( fixed mass fraction ) / else
  }//for( const ManualMassFractionConstraint &constraint : input.mass_fraction_constraints )
#endif // USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT
  
  // Note: the SubsetManifold holding `constant_parameters` fixed is created below, after the
  //  Ceres-fit empirical-form coefficient seeding, which adds a gauge-pinned coefficient.

  // Set a lower bound on relative activities to be 0, unless it is constrained
  for( size_t i = 0; i < num_nuclides; ++i )
  {
    bool is_fixed = false, is_mass_frac_constrained = false;
    for( const auto &constraint : input.act_ratio_constraints )
    {
      if( constraint.m_constrained_nuclide == cost_functor->m_isotopes[i] )
        is_fixed = true;
    }

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
    // An all-constrained elements carrier holds the element total-rel-mass scale: start 1.0,
    //  bounded below by 0, no upper bound (regardless of whether its own constraint is fixed).
    for( const auto &block : cost_functor->m_mass_frac_blocks )
    {
      if( block.spec.all_constrained && (block.carrier_index == i) )
      {
        pars[i] = 1.0;
        is_mass_frac_constrained = true;
        problem.SetParameterLowerBound( pars, static_cast<int>(i), 0.0 );
        break;
      }
    }//for( const auto &block : cost_functor->m_mass_frac_blocks )

    for( const auto &constraint : input.mass_fraction_constraints )
    {
      if( is_mass_frac_constrained )
        break;

      if( constraint.m_nuclide == cost_functor->m_isotopes[i] )
      {
        if( constraint.m_mass_fraction_lower == constraint.m_mass_fraction_upper )
        {
          is_fixed = true;
        }else
        {
          // A sigma-block `t` (carrier) or `g_k` (distribution) value - box [0.5, 1.5], start at
          //  the midpoint (t = g = 0.5).
          pars[i] = 1.0;
          is_mass_frac_constrained = true;
          problem.SetParameterLowerBound( pars, static_cast<int>(i), 0.5 );
          problem.SetParameterUpperBound( pars, static_cast<int>(i), 1.5 );
        }

        break;
      }//if( constraint.m_nuclide == cost_functor->m_isotopes[i] )
    }//for( const auto &constraint : input.mass_fraction_constraints )
#endif // USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT
    if( is_fixed )
      pars[i] = -1.0; //so we can assert on this later to make sure things are reasonable
    else if( !is_mass_frac_constrained )
      problem.SetParameterLowerBound( pars, static_cast<int>(i), 0.0 );
  }//for( size_t i = 0; i < num_nuclides; ++i )
  
  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    const auto set_bounds = [&problem, num_nuclides, pars]( const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt, size_t &index ){
      // Note: the previous inline predicate here skipped a fit-AN shield whose atomic_number
      //  was 0 (valid per check_valid), mis-aligning every subsequent parameter bound.
      if( !shield_is_present(opt) )
        return;

      if( opt->fit_atomic_number )
      {
        double lower_an = opt->lower_fit_atomic_number;
        double upper_an = opt->upper_fit_atomic_number;
        if( (lower_an == upper_an) && (lower_an == 0.0) )
        {
          lower_an = 1.0;
          upper_an = 98.0;
        }

        problem.SetParameterLowerBound( pars, static_cast<int>(index), lower_an / RelActCalc::ns_an_ceres_mult );
        problem.SetParameterUpperBound( pars, static_cast<int>(index), upper_an / RelActCalc::ns_an_ceres_mult );
        index += 1; //Add parameter for AN, only if fitting it
      }
      
      if( opt->fit_areal_density )
      {
        double lower_ad = opt->lower_fit_areal_density / PhysicalUnits::g_per_cm2;
        double upper_ad = opt->upper_fit_areal_density / PhysicalUnits::g_per_cm2;
        if( (lower_ad == upper_ad) && (lower_ad == 0.0) )
        {
          lower_ad = 0.0;
          upper_ad = RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2;
        }
        
        problem.SetParameterLowerBound( pars, static_cast<int>(index), lower_ad );
        problem.SetParameterUpperBound( pars, static_cast<int>(index), upper_ad );
      }
      
      index += 1; //Add parameter for AD, always
    };//set_bounds lambda

    size_t index = num_nuclides;
    set_bounds( input.phys_model_self_atten, index );
    for( size_t i = 0; i < input.phys_model_external_attens.size(); ++i )
      set_bounds( input.phys_model_external_attens[i], index );
    
    if( input.phys_model_use_hoerl )
    {
      const double b_lower = (0.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset;
      const double b_upper = (2.0/RelActCalc::ns_decay_hoerl_b_multiple) + RelActCalc::ns_decay_hoerl_b_offset;
      const double c_lower = (1.0E-6/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset;  //e.x, pow(-0.1889,1000/124.8) is NaN
      const double c_upper = (3.0/RelActCalc::ns_decay_hoerl_c_multiple) + RelActCalc::ns_decay_hoerl_c_offset;
      
      assert( num_parameters > 2 );
      problem.SetParameterLowerBound( pars, static_cast<int>(index), b_lower );
      problem.SetParameterUpperBound( pars, static_cast<int>(index), b_upper );
      index += 1;
      assert( index < num_parameters );
      problem.SetParameterLowerBound( pars, static_cast<int>(index), c_lower );
      problem.SetParameterUpperBound( pars, static_cast<int>(index), c_upper );
      index += 1;
    }
    
    assert( index == num_parameters );
  }else if( input.use_ceres_to_fit_eqn )
  {
    try
    {
      vector<double> rel_activities( num_nuclides ), dummy_parameters( num_parameters, 1.0 );
      for( size_t i = 0; i < num_nuclides; ++i )
        rel_activities[i] = cost_functor->relative_activity( cost_functor->m_isotopes[i], parameters );

      vector<double> fit_pars;
      fit_rel_eff_eqn_lls( eqn_form, eqn_order, cost_functor->m_isotopes, rel_activities, peak_infos, fit_pars, nullptr );
      assert( fit_pars.size() == (eqn_order + 1) );
      assert( parameters.size() == (cost_functor->m_isotopes.size() + eqn_order + 1) );
      for( size_t i = 0; i < (eqn_order + 1); ++i )
        parameters[num_nuclides + i] = fit_pars[i];

      // The counts-space residuals are exactly invariant under scaling all activities by k while
      //  dividing the curve by k (LnX: all coefficients divided by k; exponential forms:
      //  c0 -> c0 - ln(k)), leaving the Jacobian rank-deficient by one - and the unweighted
      //  residual variant is actively driven toward k -> infinity.  Pin the gauge by holding one
      //  coefficient at its LLS-seeded value; after the solve the reported activities and
      //  coefficients are re-expressed in the LLS-mode convention (average measured rel. eff.
      //  == 1; see `average_measured_rel_eff`), so both fit methods agree in values and
      //  uncertainties.
      size_t held_coef_index = 0; // exponential forms: only c0 is on the scale orbit
      if( eqn_form == RelActCalc::RelEffEqnForm::LnX )
      {
        // For LnX the orbit scales ALL coefficients, so hold the largest-magnitude one.  (If
        //  every seeded coefficient were zero the data are degenerate and the solve fails anyway.)
        for( size_t i = 0; i < (eqn_order + 1); ++i )
          if( fabs(parameters[num_nuclides + i]) > fabs(parameters[num_nuclides + held_coef_index]) )
            held_coef_index = i;
      }
      constant_parameters.push_back( static_cast<int>(num_nuclides + held_coef_index) );
    }catch( std::exception &e )
    {
      solution.m_status = ManualSolutionStatus::ErrorInitializing;
      solution.m_error_message = "Failed to fit initial relative efficiency equation: " + string(e.what());
      return solution;
    }
  }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else if( input.use_ceres_to_fit_eqn )

  if( !constant_parameters.empty() )
  {
    ceres::Manifold *subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_parameters), constant_parameters );
    problem.SetManifold( pars, subset_manifold ); //Looks like it takes ownership of subset_manifold
  }

  // Okay - we've set our problem up
  ceres::Solver::Options ceres_options;
  ceres_options.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH
  ceres_options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG
  ceres_options.linear_solver_type = ceres::DENSE_QR;
  ceres_options.logging_type = debug_printout() ? ceres::PER_MINIMIZER_ITERATION : ceres::SILENT;
  ceres_options.minimizer_progress_to_stdout = debug_printout();
  ceres_options.max_num_iterations = 1000;
  ceres_options.max_solver_time_in_seconds = 60.0;
  //ceres_options.min_trust_region_radius = 1e-10;
  
  ceres_options.use_nonmonotonic_steps = true;
  ceres_options.max_consecutive_nonmonotonic_steps = 5;
  // There are a lot more options that could be useful here!

  //parameter_tolerance = 1e-8;
  //double gradient_tolerance = 1e-10;
  //double function_tolerance = 1e-6;
  //int max_num_consecutive_invalid_steps = 5;
  
  // Setting ceres_options.num_threads >1 doesnt seem to do much (any?) good
  ceres_options.num_threads = std::thread::hardware_concurrency();
  assert( ceres_options.num_threads );
  if( !ceres_options.num_threads )
    ceres_options.num_threads = 4;
  
  //cout << "Starting parameter values: {";
  //for( size_t i = 0; i < num_parameters; ++i )
  //  cout << (i ? ", " : "") << parameters[i];
  //cout << "}\n";
  
  ceres::Solver::Summary summary;
  
  try
  {
    ceres::Solve(ceres_options, &problem, &summary);
    
    solution.m_num_function_eval_solution = static_cast<int>( cost_functor->m_ncalls );
    
    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        // good deal, all is well with our little fit
        break;
        
      case ceres::NO_CONVERGENCE:
      case ceres::FAILURE:
      case ceres::USER_FAILURE:
        throw runtime_error( "The L-M solving failed." );
        break;
    }//switch( summary.termination_type )
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorFindingSolution;
    solution.m_error_message += e.what();
    
    cerr << "RelActCalcManual::solve_relative_efficiency: Failed in solving solution: "
         << e.what() << endl;
    
    return solution;
  }//try / catch to fit the solution
  
  
  if( debug_printout() )
  {
    std::cout << summary.BriefReport() << "\n";
    //std::cout << summary.FullReport() << "\n";
    const auto nmicro = std::chrono::duration<double, std::micro>(std::chrono::high_resolution_clock::now() - start_time).count();
    cout << "Took " << solution.m_num_function_eval_solution << " calls and " << setprecision(6) << nmicro << " us to solve." << endl;
    cout << "Final parameter values: {";
    for( size_t i = 0; i < num_parameters; ++i )
      cout << (i ? ", " : "") << parameters[i];
    cout << "}\n";
    cout << "Chi2=" << summary.final_cost << " (from initial value " << summary.initial_cost << ")\n\n";
  }//if( debug_printout() )


  solution.m_fit_parameters = parameters;

  // Chi2 under the weights actually minimized (i.e., including any
  //  GenericPeakInfo::m_base_rel_eff_uncert).  Note the gauge pin for Ceres-fit empirical forms
  //  is a constant parameter, not an extra residual, so the data channels are the whole cost
  //  here; `peaks_only_chi2` (== 2*summary.final_cost on this path) is used so this and the
  //  profile scan's chi2 values are produced by the identical expression.
  solution.m_chi2_fit_weights = peaks_only_chi2( parameters );
#if( !USE_RESIDUAL_TO_BREAK_DEGENERACY )
  //  Only equal to Ceres' cost when every residual is a peak; the degeneracy-breaking channel
  //  would be part of `final_cost` but not of `peaks_only_chi2`.
  assert( fabs(solution.m_chi2_fit_weights - 2.0*summary.final_cost)
          <= 1.0E-6*(std::max)(1.0, fabs(solution.m_chi2_fit_weights)) );
#endif

  // Count the effective number of fitted parameters once, from the same bookkeeping the solver
  //  used: `constant_parameters` holds act-ratio controlled and fixed mass-fraction slots,
  //  non-fit shield parameters, and (for Ceres-fit empirical forms) the gauge-pinned
  //  coefficient.  In the LLS fit mode the (eqn_order + 1) curve coefficients are profiled by
  //  the inner LLS fit, with one combination fixed by the average-measured-rel-eff
  //  normalization convention, so they add a further (eqn_order + 1) - 1 effective parameters.
  //  The residual count deliberately does NOT enter here: `number_residuals()` may include
  //  channels that are not data (the profile-likelihood penalty channel), and the DOF of the fit
  //  is set by the peaks.
  {
    size_t num_effective_pars = num_parameters - constant_parameters.size();
    if( (eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel) && !input.use_ceres_to_fit_eqn )
      num_effective_pars += eqn_order;  // == (eqn_order + 1) - 1
    assert( cost_functor->m_input.peaks.size() == num_peaks );
    solution.m_dof = static_cast<int>( num_peaks ) - static_cast<int>( num_effective_pars );
  }

  // For the Ceres-fit empirical forms, results are reported in the LLS-mode gauge (average
  //  measured rel. eff. == 1): activities are multiplied by m(x), and the curve coefficients
  //  compensated, below.  m(x)*A_i(x) is exactly invariant along the scale orbit, so the
  //  derived rel-act covariance is well-defined and directly comparable between fit methods.
  bool gauge_normalize = (eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel)
                         && input.use_ceres_to_fit_eqn;
  double gauge_mult = 1.0;
  if( gauge_normalize )
  {
    bool any_peak_floored = false;
    gauge_mult = cost_functor->average_measured_rel_eff( parameters, &any_peak_floored );

    // `any_peak_floored` matters as much as the finiteness test: a nuclide fitted to zero activity
    //  that owns a peak makes that peak's term ~1e6 instead of infinite, so the average comes back
    //  finite and would silently scale every reported activity by ~1e5.
    if( any_peak_floored || (gauge_mult <= 0.0) || isnan(gauge_mult) || isinf(gauge_mult) )
    {
      //shouldnt happen; leave results in the solver gauge if it somehow does
      gauge_mult = 1.0;
      gauge_normalize = false;
      if( any_peak_floored )
        solution.m_warnings.push_back( "At least one peak has no fitted source activity left to"
              " explain it, so the relative activities are reported in the solver's own"
              " normalization rather than the usual one; compare ratios, not absolute values." );
    }
  }

  // Whatever was decided above is what "reported activity" means from here on - including for a
  //  ProfileTarget::Type::RelativeActivity scan.
  cost_functor->set_gauge_normalizes_activities( gauge_normalize );

  try
  {
    ceres::Covariance::Options cov_options;
    //cov_options.algorithm_type = ceres::CovarianceAlgorithmType::SPARSE_QR; //faster, but not capable of computing the covariance if the Jacobian is rank deficient.
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::DENSE_SVD;
    cov_options.null_space_rank = -1;
    // Drop near-null directions (sigma_i/sigma_max < 1e-6) from the covariance instead of
    //  letting their ~1/sigma^2 variances blow up the reported uncertainties; these are the
    //  structurally near-degenerate directions (e.g., self-AD vs ext-AD vs Hoerl tilt in the
    //  physical model).  Same remedy as the RelActCalcAuto solver; the manual problem's
    //  parameters are already near-unit-scaled (activity multiples, AN/50, AD in g/cm2,
    //  offset Hoerl), so a raw-Jacobian cutoff is appropriate here without the column
    //  equilibration the Auto solver needs.
    cov_options.min_reciprocal_condition_number = 1.0E-12;
    cov_options.num_threads = ceres_options.num_threads;
    
    vector<double> uncertainties( num_nuclides, 0.0 ), uncerts_squared( num_nuclides, 0.0 );
    
    ceres::Covariance covariance(cov_options);
    vector<const double*> parameter_blocks;
    vector<pair<const double*, const double*> > covariance_blocks;

    parameter_blocks.push_back( pars );
    covariance_blocks.push_back( {pars, pars} );
    
    solution.m_nonlin_covariance.clear();
    if( input.point_estimate_only )
    {
      // Skip all covariance work: derived uncertainties carry their "not available" sentinels
      //  (empty matrices, -1 uncertainties), but point estimates, chi2, DOF, and the shield fit
      //  values below are still produced.  Used by the AN scan and profile-likelihood sub-solves.
    }else if( !covariance.Compute(covariance_blocks, &problem) )
    {
      cerr << "Failed to compute final covariances!" << endl;
      solution.m_warnings.push_back( "Failed to compute final covariances." );
    }else
    {
      // row-major order: the elements of the first row are consecutively given, followed by second
      //                  row contents, etc
      vector<double> row_major_covariance( num_parameters * num_parameters );
      
      const bool success = covariance.GetCovarianceMatrix( parameter_blocks, row_major_covariance.data() );
      assert( success );
      if( !success )
        throw runtime_error( "Failed to get covariance matrix - maybe didnt add all covariance blocks?" );
      
      solution.m_nonlin_covariance.resize( num_parameters, vector<double>(num_parameters,0.0) );
      for( size_t row = 0; row < num_parameters; ++row )
      {
        for( size_t col = 0; col < num_parameters; ++col )
          solution.m_nonlin_covariance[row][col] = row_major_covariance[row*num_parameters + col];
      }//for( size_t row = 0; row < num_nuclides; ++row )

      // Inflate the covariance when the data scatter about the model exceeds the assumed
      //  uncertainties (the classic "variance of unit weight" / Birge rescale): the raw
      //  covariance then systematically understates the parameter uncertainties.  Never shrink
      //  (a scale below 1 is clamped away), and skip it entirely for unweighted fits.
      //
      //  We deliberately do NOT use the usual chi2/dof for the scale.  In this problem the
      //  excess is usually NOT a uniform underestimate of the peak uncertainties - it is a
      //  handful of individually mis-fit peaks (wrong continuum type, insufficient skew, an
      //  unmodeled interference).  Those land many sigma off when the peak is high-statistics,
      //  yet barely move the enrichment; using the mean of the squared pulls would let them
      //  inflate EVERY reported uncertainty.  Measured on the spec184 problem 20260816: with a
      //  LnX order-4 fit the median |pull| is 0.93 (i.e. the typical peak fits fine) while the
      //  top 3 of 24 peaks carry 55% of the chi2 - chi2/dof = 6.1 would inflate sigmas by 2.5x,
      //  where the robust estimator below gives 1.7x.  When the excess really is broad, the two
      //  agree: the same spectrum's physical-model fit gives 2.12x (mean) vs 2.06x (robust).
      //
      //  The robust estimator is the median of the squared pulls, scaled by the median of a
      //  chi2(1) variate (0.4549) so it equals 1 for well-behaved data, and by n/dof to undo
      //  the variance the fit itself removes from the residuals (average leverage; makes the
      //  estimator unbiased in the no-outlier Gaussian case).
      //  \sa RelActCalcAuto::RelActAutoSolution::m_cov_scale, which uses the chi2/dof form.
      {
        bool any_unweighted = false;
        for( const GenericPeakInfo &peak : cost_functor->m_input.peaks )
          any_unweighted |= (peak.m_base_rel_eff_uncert < -1.0E-9);

        solution.m_cov_scale = 1.0;

        if( (solution.m_dof > 0) && !any_unweighted && input.widen_uncerts_for_scatter )
        {
          // Peak channels only - a profile-likelihood penalty channel, or the
          //  USE_RESIDUAL_TO_BREAK_DEGENERACY normalization channel, is not data and would drag
          //  the median of the squared pulls.
          const size_t num_resids = num_peaks;
          vector<double> residuals( cost_functor->number_residuals(), 0.0 );
          cost_functor->eval( parameters, residuals.data() );

          vector<double> sq_pulls( num_resids );
          for( size_t i = 0; i < num_resids; ++i )
            sq_pulls[i] = residuals[i]*residuals[i];

          vector<double> sorted_sq = sq_pulls;
          std::sort( begin(sorted_sq), end(sorted_sq) );
          const double median_sq = (num_resids % 2)
                  ? sorted_sq[num_resids/2]
                  : 0.5*(sorted_sq[num_resids/2 - 1] + sorted_sq[num_resids/2]);

          const double chi2_1_median = 0.4549364; //median of a chi-squared(1) distribution
          const double leverage_corr = static_cast<double>(num_resids) / solution.m_dof;
          solution.m_cov_scale = (std::max)( 1.0, median_sq * leverage_corr / chi2_1_median );

          // Point at the peaks driving the misfit - the actionable cure is re-fitting those
          //  peaks, not accepting a wider error bar.
          vector<pair<double,double>> abs_pull_and_energy; //|pull|, energy
          for( size_t i = 0; (i < num_resids) && (i < cost_functor->m_input.peaks.size()); ++i )
            abs_pull_and_energy.emplace_back( fabs(residuals[i]),
                                              cost_functor->m_input.peaks[i].m_energy );
          std::sort( begin(abs_pull_and_energy), end(abs_pull_and_energy),
                     []( const pair<double,double> &l, const pair<double,double> &r ){
                       return l.first > r.first;
                     } );

          if( !abs_pull_and_energy.empty() && (abs_pull_and_energy[0].first > 3.0) )
          {
            string msg = "Some peaks are much further from the fit than their uncertainties allow:";
            for( size_t i = 0; (i < 3) && (i < abs_pull_and_energy.size()); ++i )
            {
              if( abs_pull_and_energy[i].first <= 3.0 )
                break;
              msg += " " + SpecUtils::printCompact(abs_pull_and_energy[i].second, 5) + " keV ("
                     + SpecUtils::printCompact(abs_pull_and_energy[i].first, 2) + "&sigma;)";
            }
            msg += ".  This usually indicates a peak-fit issue (continuum type, peak skew, or an"
                   " unmodeled interference) rather than a source-composition effect - checking"
                   " those peak fits is worthwhile.  Reported uncertainties use an"
                   " outlier-insensitive estimate of the scatter, so a few such peaks do not"
                   " inflate every uncertainty.";
            solution.m_warnings.push_back( msg );
          }//if( there are strongly deviating peaks )
        }//if( we can estimate the scatter )

        if( solution.m_cov_scale > 1.0 )
        {
          for( vector<double> &row : solution.m_nonlin_covariance )
            for( double &val : row )
              val *= solution.m_cov_scale;
        }

        if( solution.m_cov_scale > 2.25 ) // i.e., sigmas inflated by more than 1.5x
          solution.m_warnings.push_back( "The scatter of the data about the fitted model exceeds"
                " the statistical (plus any additional) uncertainties, in a way that is spread"
                " across the peaks rather than confined to a few: reported uncertainties have"
                " been inflated by " + SpecUtils::printCompact( std::sqrt(solution.m_cov_scale), 3 )
                + "x (chi2/dof = "
                + SpecUtils::printCompact( solution.m_chi2_fit_weights / solution.m_dof, 3 )
                + " under the fit weights)." );
      }

      // Compute the full Jacobian of the reported relative activities with respect to the fit
      //  parameters, d(RelAct_i)/d(par_j), by automatic differentiation, and from it the
      //  relative-activity covariance C_relact = J * C * J^T.  Every derived uncertainty
      //  (per-isotope sigmas, activity ratios, mass-fraction variations) is computed from this
      //  one matrix, so activity norms, act-ratio constraint chains, and mass-fraction block
      //  decodes (a constrained nuclides activity co-varying with the rest of its element) are
      //  all consistently accounted for.  `m_nonlin_covariance` itself stays pristine in
      //  parameter space.  Parameters held constant (SubsetManifold) have all-zero covariance
      //  rows, so seeding them costs nothing.
      if constexpr ( use_auto_diff )
      {
        solution.m_rel_act_jacobian.assign( num_nuclides, vector<double>(num_parameters, 0.0) );

        for( size_t chunk_start = 0; chunk_start < num_parameters; chunk_start += auto_diff_stride )
        {
          vector<ceres::Jet<double,auto_diff_stride>> input_jets( begin(parameters), end(parameters) );
          const size_t num_seed = std::min( num_parameters - chunk_start, static_cast<size_t>(auto_diff_stride) );
          for( size_t k = 0; k < num_seed; ++k )
            input_jets[chunk_start + k].v[k] = 1.0;

          // In the LLS-reporting gauge the reported activity is m(x)*A_i(x) (see
          //  `gauge_normalize` above), so differentiate that product.
          const ceres::Jet<double,auto_diff_stride> chunk_gauge = gauge_normalize
                     ? cost_functor->average_measured_rel_eff( input_jets )
                     : ceres::Jet<double,auto_diff_stride>( 1.0 );

          for( size_t i = 0; i < num_nuclides; ++i )
          {
            const ceres::Jet<double,auto_diff_stride> rel_act
                     = chunk_gauge * cost_functor->relative_activity( cost_functor->m_isotopes[i], input_jets );
            for( size_t k = 0; k < num_seed; ++k )
              solution.m_rel_act_jacobian[i][chunk_start + k] = rel_act.v[k];
          }
        }//for( loop over stride-sized chunks of parameters )

        // C_relact = J * C * J^T  (problem sizes are tiny, so plain loops are fine)
        vector<vector<double>> jac_cov( num_nuclides, vector<double>(num_parameters, 0.0) );
        for( size_t i = 0; i < num_nuclides; ++i )
          for( size_t l = 0; l < num_parameters; ++l )
            for( size_t k = 0; k < num_parameters; ++k )
              jac_cov[i][k] += solution.m_rel_act_jacobian[i][l] * solution.m_nonlin_covariance[l][k];

        solution.m_rel_act_covariance.assign( num_nuclides, vector<double>(num_nuclides, 0.0) );
        for( size_t i = 0; i < num_nuclides; ++i )
          for( size_t j = 0; j < num_nuclides; ++j )
            for( size_t k = 0; k < num_parameters; ++k )
              solution.m_rel_act_covariance[i][j] += jac_cov[i][k] * solution.m_rel_act_jacobian[j][k];
      }else
      {
        // dont expect to ever not use auto-diff, so wont worry about this case
        static_assert( use_auto_diff, "Numeric diff not implemented for computing the manual solver rel. act. covariance" );
      }
    }//if( we failed to get covariance ) / else
    
    // Compute the Jacobian (the rank-deficiency check below auto-skips when this is skipped)
    if( !input.point_estimate_only )
    try
    {
      const size_t num_residuals = cost_functor->number_residuals();
      vector<double> residuals( num_residuals );
      vector<double> jacobian( num_parameters * num_residuals );

      const double * const parameters_ptr = parameters.data();
      double * const residuals_ptr = &(residuals[0]);
      double * jacobians_ptr = &(jacobian[0]);  //We only have a single paramater block; the pointer passed into Ceres expects first index to index into parameter block

      const bool success = cost_function->Evaluate( &parameters_ptr, residuals_ptr, &jacobians_ptr );
      if( !success )
        throw std::runtime_error( "Failed to evaluate the cost function." );

      // Keep only the peak rows, so the k-index stays the peak index no matter which extra
      //  channels the functor carries (profile-likelihood penalty, degeneracy pin).
      solution.m_nonlin_jacobian.resize( num_peaks, vector<double>(num_parameters, 0.0) );
      for( size_t k = 0; k < num_peaks; ++k )
      {
        for( size_t i = 0; i < num_parameters; ++i )
          solution.m_nonlin_jacobian[k][i] = jacobian[k*num_parameters + i];
      }
    }catch(const std::exception& e)
    {
      cerr << "Failed to compute final Jacobian! - " << e.what() << endl;
      solution.m_warnings.push_back( "Failed to compute Jacobian: " + string(e.what()) );
    }//try / catch to compute Jacobian

    // Detect effectively-unconstrained (near-degenerate) parameter directions, and warn.  With
    //  DENSE_SVD + null_space_rank = -1 + the min_reciprocal_condition_number above, such
    //  directions are dropped from the covariance (truncated pseudo-inverse) rather than
    //  blowing it up - but the per-parameter uncertainties then omit those directions.  The
    //  parameters are near-unit-scaled, so a raw-Jacobian singular value cutoff is meaningful.
    if( !solution.m_nonlin_jacobian.empty() && !solution.m_nonlin_covariance.empty() )
    {
      vector<size_t> free_pars;
      for( size_t i = 0; i < num_parameters; ++i )
      {
        if( !std::count( begin(constant_parameters), end(constant_parameters), static_cast<int>(i) ) )
          free_pars.push_back( i );
      }

      const size_t num_resids = solution.m_nonlin_jacobian.size();
      if( !free_pars.empty() && (num_resids >= free_pars.size()) )
      {
        Eigen::MatrixXd jacobian_mat( num_resids, free_pars.size() );
        for( size_t r = 0; r < num_resids; ++r )
          for( size_t c = 0; c < free_pars.size(); ++c )
            jacobian_mat(r,c) = solution.m_nonlin_jacobian[r][free_pars[c]];

        const Eigen::JacobiSVD<Eigen::MatrixXd> svd( jacobian_mat );
        const Eigen::VectorXd &sing_vals = svd.singularValues(); //sorted decreasing

        size_t num_deficient = 0;
        for( Eigen::Index i = 0; i < sing_vals.size(); ++i )
          num_deficient += ((sing_vals.size() > 0) && (sing_vals(i) < 1.0E-6*sing_vals(0)));

        if( num_deficient > 0 )
          solution.m_warnings.push_back( "The fit is rank-deficient: "
                + std::to_string(num_deficient) + " of " + std::to_string(free_pars.size())
                + " fitted parameters are effectively unconstrained by the data"
                " (near-degenerate directions).  Per-parameter uncertainties omit these"
                " directions, so uncertainties of quantities that depend on them may be"
                " understated." );
      }//if( sensible matrix dimensions )
    }//if( have jacobian and covariance )


    // Compute the relative activities
    vector<double> rel_activities( num_nuclides );
    
    for( size_t i = 0; i < num_nuclides; ++i )
      rel_activities[i] = cost_functor->relative_activity( cost_functor->m_isotopes[i], parameters );
    
    if( (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel) || input.use_ceres_to_fit_eqn )
    {
      solution.m_rel_eff_eqn_coefficients = {pars + num_nuclides, pars + num_parameters };

      if( gauge_normalize )
      {
        // Re-express the curve in the LLS-mode gauge: the reported activities are multiplied by
        //  `gauge_mult` (below), so the curve is divided by it.
        if( eqn_form == RelActCalc::RelEffEqnForm::LnX )
        {
          for( double &coef : solution.m_rel_eff_eqn_coefficients )
            coef /= gauge_mult;
        }else
        {
          solution.m_rel_eff_eqn_coefficients[0] -= std::log( gauge_mult );
        }
      }//if( gauge_normalize )

      if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        if( !solution.m_nonlin_covariance.empty() )
        {
          // The covariance matrix for the relative efficiency equation is the lower-right
          //  submatrix of the full covariance matrix.  For the Physical Model this marginal
          //  covariance is the right thing: the DRF fixes the curve's absolute scale, so the
          //  shield/Hoerl parameters are identifiable on their own.
          const size_t num_rel_eff_params = num_parameters - num_nuclides;
          assert( solution.m_nonlin_covariance.size() == num_parameters );

          solution.m_rel_eff_eqn_covariance.resize( num_rel_eff_params, vector<double>(num_rel_eff_params, 0.0) );
          for( size_t i = num_nuclides; i < num_parameters; ++i )
          {
            assert( solution.m_nonlin_covariance[i].size() == num_parameters );
            for( size_t j = num_nuclides; j < num_parameters; ++j )
              solution.m_rel_eff_eqn_covariance[i-num_nuclides][j-num_nuclides] = solution.m_nonlin_covariance[i][j];
          }
        }//if( we have the covariance matrix )
      }else if( !input.point_estimate_only )
      {
        // Empirical forms: use the coefficient covariance CONDITIONAL on the fitted activities -
        //  i.e., the same one the LLS fit mode produces - rather than the marginal sub-block of
        //  the joint covariance.
        //
        //  For these forms only the PRODUCT (rel. eff. curve x activities) is determined by the
        //  data; the split between them is a gauge choice (fixed here by normalizing the average
        //  measured rel. eff. to 1).  The marginal coefficient covariance therefore contains the
        //  curve<->activity trade, and using it for the plotted band gives absurd widths -
        //  measured on spec184 20260816, at 186 keV the band went from 0.047 (conditional) to
        //  4.5 (LnX), 25 (LnXLnY) and 1.8 (FRAM Empirical), on a curve whose value is ~1.1.
        //  It would also double-count: the plotted data points are themselves computed with the
        //  fitted activities, so they move WITH the curve under that trade.  The band the chart
        //  wants is "how well does this curve describe these points", which is the conditional
        //  one.  (A refinement would be to take the conditional covariance from the coefficient
        //  block of the Ceres Jacobian instead of re-fitting; that is exact for the counts-space
        //  objective, but changes band values, so it is left for a deliberate change.)
        vector<double> band_coefficients;
        fit_rel_eff_eqn_lls( eqn_form, eqn_order, cost_functor->m_isotopes, rel_activities, peak_infos,
                            band_coefficients, &(solution.m_rel_eff_eqn_covariance) );

        if( solution.m_cov_scale > 1.0 )
        {
          for( vector<double> &row : solution.m_rel_eff_eqn_covariance )
            for( double &val : row )
              val *= solution.m_cov_scale;
        }
      }//if( FramPhysicalModel ) / else
    }else
    {
      fit_rel_eff_eqn_lls( eqn_form, eqn_order, cost_functor->m_isotopes, rel_activities, peak_infos,
                          solution.m_rel_eff_eqn_coefficients, &(solution.m_rel_eff_eqn_covariance) );
      assert( solution.m_rel_eff_eqn_coefficients.size() == (eqn_order + 1) );

      // Propagate the chi2/dof covariance inflation to the LLS-fit curve covariance too, so the
      //  plotted uncertainty band is consistent with the (inflated) activity uncertainties.
      if( solution.m_cov_scale > 1.0 )
      {
        for( vector<double> &row : solution.m_rel_eff_eqn_covariance )
          for( double &val : row )
            val *= solution.m_cov_scale;
      }
    }
    
    for( size_t i = 0; i < cost_functor->m_isotopes.size(); ++i )
    {
      const string &iso = cost_functor->m_isotopes[i];
      
      IsotopeRelativeActivity rel_act;
      rel_act.m_isotope = iso;

      // `gauge_mult` is 1.0 except for Ceres-fit empirical forms (LLS-gauge reporting)
      rel_act.m_rel_activity = gauge_mult * cost_functor->relative_activity( iso, parameters );
      
      if( solution.m_rel_act_covariance.empty() )
      {
        rel_act.m_rel_activity_uncert = -1.0;
      }else
      {
        assert( i < solution.m_rel_act_covariance.size() );
        assert( i < solution.m_rel_act_covariance[i].size() );

        // The full-Jacobian relative-activity covariance already accounts for activity norms,
        //  act-ratio constraint chains, and mass-fraction block decodes (a constrained nuclide
        //  co-varying with the rest of its element), so this one expression covers every
        //  constraint case - including a fixed mass-fraction nuclide in a mixed element, whose
        //  relative activity still varies with the element total.
        rel_act.m_rel_activity_uncert = std::sqrt( std::max( 0.0, solution.m_rel_act_covariance[i][i] ) );
      }//if( no covariance ) / else
      
      solution.m_rel_activities.push_back( std::move(rel_act) );
    }//for( loop over relative activities )
    
    if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      const auto get_res = [num_nuclides, &parameters, &solution]( const RelActCalc::PhysicalModelShieldInput &orig_opt, 
                                                              size_t &shield_index ) -> unique_ptr<RelEffSolution::PhysModelShieldFit> {
        if( !shield_is_present(&orig_opt) )
          return nullptr;
        
        assert( solution.m_nonlin_covariance.empty()
                || (solution.m_nonlin_covariance.size() > shield_index) );

        auto shield_result = std::make_unique<RelEffSolution::PhysModelShieldFit>();
        shield_result->m_material = orig_opt.material;
        if( !shield_result->m_material )
        {
          if( orig_opt.fit_atomic_number )
          {
            shield_result->m_atomic_number = parameters[shield_index] * RelActCalc::ns_an_ceres_mult;
            
            if( !solution.m_nonlin_covariance.empty()  )
            {
              assert( shield_index < solution.m_nonlin_covariance.size() );
              assert( orig_opt.fit_atomic_number
                     || (solution.m_nonlin_covariance[shield_index][shield_index] == 0.0) );
              
              if( orig_opt.fit_atomic_number )
                shield_result->m_atomic_number_uncert = sqrt( std::max(0.0, solution.m_nonlin_covariance[shield_index][shield_index]) ) * RelActCalc::ns_an_ceres_mult;
            }//if( !solution.m_nonlin_covariance.empty()  )
          
            shield_index += 1;
          }else
          {
            shield_result->m_atomic_number = orig_opt.atomic_number;
          }
        }//if( !shield_result->m_material )
        
        const double ad_g_cm2 = parameters[shield_index];
        shield_result->m_areal_density = ad_g_cm2 * PhysicalUnits::g_per_cm2;

        if( !solution.m_nonlin_covariance.empty()  )
        {
          //assert( orig_opt.fit_areal_density
          //        || (solution.m_nonlin_covariance[shield_index][shield_index] == 0.0) );
          
          if( orig_opt.fit_areal_density )
            shield_result->m_areal_density_uncert = sqrt( std::max(0.0, solution.m_nonlin_covariance[shield_index][shield_index]) ) * PhysicalUnits::g_per_cm2;
        }//if( !solution.m_nonlin_covariance.empty()  )
        
        shield_index += 1; //increment for areal density
        
        return std::move(shield_result);
      };//get_res lamda

      size_t shield_index = num_nuclides;
      if( input.phys_model_self_atten )
        solution.m_phys_model_self_atten_shield = get_res( *input.phys_model_self_atten, shield_index );

      for( const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt : input.phys_model_external_attens )
      {
        assert( opt );
        if( opt )
          solution.m_phys_model_external_atten_shields.push_back( get_res( *opt, shield_index ) );
      }//for( size_t i = 0; i < input.phys_model_external_attens.size(); ++i )
      
      assert( (shield_index + (solution.m_input.phys_model_use_hoerl ? 2 : 0)) == num_parameters );
    }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    
    // We'll manually compute the Chi2 here, so we only take into account statistical uncertainties,
    //  (i.e., we ignore #GenericPeakInfo::m_base_rel_eff_uncert)
    solution.m_chi2 = 0.0;
    
    // Note: `m_dof` was computed from the solver parameter bookkeeping right after the solve
    //  (num peaks minus effective fitted parameters); it can legitimately be 0.
    if( solution.m_dof < 0 )
      throw runtime_error( "There are only " + std::to_string(cost_functor->m_input.peaks.size())
                          + " peaks, but the fit has "
                          + std::to_string( static_cast<int>(cost_functor->m_input.peaks.size()) - solution.m_dof )
                          + " effective parameters." );


    const vector<double> &rel_eff_coefs = solution.m_rel_eff_eqn_coefficients;
    ManualGenericRelActFunctor::PhysModelRelEqnDef<double> phys_mode_rel_eqn_input;
    if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      phys_mode_rel_eqn_input = ManualGenericRelActFunctor::make_phys_eqn_input( input, rel_eff_coefs );
    
    bool used_unweighted = false, used_add_uncert = false;
    for( const GenericPeakInfo &peak : cost_functor->m_input.peaks )
    {
      double curve_val;
      if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        curve_val = RelActCalc::eval_physical_model_eqn( peak.m_energy,
                                                  phys_mode_rel_eqn_input.self_atten,
                                                  phys_mode_rel_eqn_input.external_attens,
                                                  phys_mode_rel_eqn_input.det.get(),
                                                  phys_mode_rel_eqn_input.hoerl_b,
                                                  phys_mode_rel_eqn_input.hoerl_c );
      }else
      {
        curve_val = RelActCalc::eval_eqn( peak.m_energy, eqn_form, rel_eff_coefs );
      }
      
      used_add_uncert |= (peak.m_base_rel_eff_uncert > 1.0E-9);
      used_unweighted |= (peak.m_base_rel_eff_uncert < -1.0E-9);
      
      double expected_src_counts = 0.0;
      for( const GenericLineInfo &line : peak.m_source_gammas )
      {
        // `gauge_mult` (1.0 except Ceres-fit empirical forms) keeps the activities consistent
        //  with the gauge-rescaled curve coefficients; their product is gauge-invariant.
        const double rel_activity = gauge_mult * cost_functor->relative_activity( line.m_isotope, parameters );
        expected_src_counts += rel_activity * line.m_yield;
      }//for( const RelActCalc::GammaLineInfo &line : peak.m_source_gammas )
      
      const double expected_counts = expected_src_counts * curve_val;
      if( peak.m_counts_uncert > 0.0 ) //unweighted fits can carry zero statistical uncert
        solution.m_chi2 += std::pow( (expected_counts - peak.m_counts) / peak.m_counts_uncert, 2.0 );
    }//for( loop over energies to evaluate at )
    
    // `m_chi2` is documented as being in terms of the peaks' STATISTICAL uncertainties.  When the
    //  automatic estimate widened those uncertainties by a common multiple, undo it here so the
    //  displayed chi2/dof still reports how far the data sit from the model in counting-statistics
    //  terms - that number is the diagnostic that the model does not describe the data, and
    //  ballooning the uncertainties should not hide it.
    if( auto_stat_multiple > 1.0 )
      solution.m_chi2 *= (auto_stat_multiple * auto_stat_multiple);

    if( used_add_uncert )
      solution.m_warnings.push_back( "Additional uncertainties were applied to peaks"
                                     " - the result uncertainties include these, so may not be"
                                     " reliable to interpret. The &chi;<sup>2</sup>/DOF does"
                                     " not include the add. uncerts." );
    
    if( used_unweighted )
      solution.m_warnings.push_back( "Fit to rel. eff. was unweighted, so uncertainties may not have much meaning." );
    
    
    solution.m_status = ManualSolutionStatus::Success;
    solution.m_num_function_eval_total = static_cast<int>( cost_functor->m_ncalls );
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorGettingSolution;
    solution.m_error_message = e.what();
    cerr << "RelActCalcManual::solve_relative_efficiency: Failed to get solution after solving: "
         << e.what() << endl;

    return solution;
  }//try / catch to get the solution


  // Profile-likelihood (asymmetric) intervals - see RelEffInput::profile_targets.
  //
  //  Each scan point re-solves this same `problem`, from the parameters the nominal solve landed
  //  on, with one extra residual `w*(g(x) - target)` armed in the cost functor.  At the optimum
  //  of `C(x) + w^2*(g(x) - target)^2` we have `grad C = -lambda*grad g`, which is exactly the
  //  stationarity condition of "minimise C subject to g == g(x*)" - so the point is an EXACT
  //  point of the profile curve at the ACHIEVED g, whatever `w` was.  We therefore record the
  //  achieved value and discard the requested one, which is what lets `w` stay small (no
  //  penalty-method ill-conditioning) and removes any need to iterate onto exact targets.
  //
  //  Fixing the parameter instead would NOT work: that gives `grad C` parallel to a coordinate
  //  axis, which is stationary for the wrong constraint, and yields intervals that are too narrow.
  if( (solution.m_status == ManualSolutionStatus::Success) && !input.profile_targets.empty() )
  {
    const vector<double> nominal_parameters = parameters;

    // Restoring `parameters` is load-bearing: `pars` points into it, so without this everything
    //  downstream would see the last trial point rather than the solution.
    DoWorkOnDestruct restore_nominal( [&parameters, &nominal_parameters, cost_functor](){
      std::copy( begin(nominal_parameters), end(nominal_parameters), begin(parameters) );
      cost_functor->arm_profile( -1, 0.0, 0.0, 0.0 );
    } );

    // Same iteration budget as the nominal solve: a scan point is only usable if it CONVERGED (a
    //  point that ran out of iterations is not stationary, so it sits above the profile curve),
    //  and the Ceres-fit empirical forms genuinely need a few hundred iterations on some steps -
    //  cutting the budget just turns those steps into "at bound".  Warm starts mean the usual
    //  step takes a handful.
    ceres::Solver::Options profile_options = ceres_options;

    // The covariance-based uncertainties carry the max(1, scatter) inflation (see m_cov_scale);
    //  scale the chi2 thresholds by the same factor so the two agree in the Gaussian limit -
    //  i.e. the profile also reads over-scatter as under-stated per-peak errors, rather than as
    //  a sharper likelihood.
    const double chi2_scale = (std::max)( 1.0, solution.m_cov_scale );

    double max_quantile = 0.0;
    vector<double> quantiles;
    for( const double cl : input.profile_confidence_levels )
    {
      // Central two-sided interval: the chi2(1 dof) CDF at delta_chi2 already folds both Gaussian
      //  tails, so quantile(chi_squared(1), CL) is the threshold (a one-sided limit would instead
      //  use 2*CL - 1; see DetectionLimitCalc::decon_limit_delta for that convention).
      const double q = ((cl > 0.0) && (cl < 1.0))
              ? boost::math::quantile( boost::math::chi_squared_distribution<double>(1.0), cl )
              : -1.0;
      quantiles.push_back( q );
      max_quantile = (std::max)( max_quantile, q );
    }

    if( max_quantile <= 0.0 )
    {
      solution.m_warnings.push_back( "No valid confidence level was given for the"
                                     " profile-likelihood scan." );
    }else
    {
      for( size_t target_index = 0; target_index < input.profile_targets.size(); ++target_index )
      {
        const ProfileTarget &target = input.profile_targets[target_index];
        const bool is_mass_frac = (target.m_type == ProfileTarget::Type::MassFraction);
        const bool is_rel_act = (target.m_type == ProfileTarget::Type::RelativeActivity);

        ProfileResult result;
        result.target = target;
        result.nominal_chi2 = solution.m_chi2_fit_weights;

        try
        {
          std::copy( begin(nominal_parameters), end(nominal_parameters), begin(parameters) );

          // The floor keeps a vanishing denominator (whole element, or the ratio denominator,
          //  going to zero) from making the quantity - and its gradient - non-finite.  A plain
          //  relative activity has no denominator; its "denom" below is just the scale the
          //  first-order sigma is expressed against, i.e. the activity itself.
          double denom_nominal = 0.0;
          if( is_mass_frac )
          {
            for( const map<string,double>::value_type &iso_sa : target.m_specific_activities )
              denom_nominal += cost_functor->relative_activity( iso_sa.first, parameters ) / iso_sa.second;
          }else
          {
            denom_nominal = cost_functor->relative_activity(
                          is_rel_act ? target.m_nuclide : target.m_denom_nuclide, parameters );
          }

          if( !(denom_nominal > 0.0) )
            throw runtime_error( "the quantity is zero at the solution, so it cannot be profiled"
                                 " (the scan is multiplicative about the fitted value)" );

          const double denom_floor = 1.0E-8 * denom_nominal;
          result.nominal_value = cost_functor->profile_quantity( target, parameters, denom_floor );

          if( !(result.nominal_value > 0.0) || isinf(result.nominal_value) )
            throw runtime_error( "the quantity is not positive at the solution" );

          // With a reporting gauge multiple and only ONE nuclide in the fit, the reported activity
          //  is `(1/P)*sum_p C_p/y_p` exactly - every `A_t` cancels out of `m(x)*A_t` - so it is a
          //  constant of the data with no dependence on the fit at all.  (Its covariance-based
          //  sigma is ~0 for the same reason: the derivative row is identically zero.)  There is
          //  nothing to profile, and saying so beats reporting a zero-width interval.
          if( is_rel_act && cost_functor->gauge_normalizes_activities()
              && (cost_functor->m_isotopes.size() == 1) )
            throw runtime_error( "with an empirical equation fit by Ceres, the reported activity of"
                    " the only nuclide in the problem is fixed by the reporting normalization and"
                    " carries no information from the fit; use the LnX form or the physical model"
                    " for an absolute activity uncertainty" );

          // A RelativeActivity scan must walk exactly the number the solution reports - if the
          //  functor's gauge decision ever drifted from the reporting one, this is where it shows.
          assert( !is_rel_act
                  || (fabs(result.nominal_value - gauge_mult*cost_functor->relative_activity(target.m_nuclide, parameters))
                      <= 1.0E-9*result.nominal_value) );

          // First-order (covariance) sigma of the quantity, used only to size the first step - the
          //  scan re-estimates it from the first solved point, so a crude value is fine.  Note this
          //  carries the m_cov_scale inflation, which is divided back out to get the RAW-chi2 sigma
          //  the thresholds are expressed in.
          double sigma_raw = 0.0;
          if( !solution.m_rel_act_covariance.empty() )
          {
            // `m_rel_act_covariance` is the covariance of the REPORTED activities, which for
            //  Ceres-fit empirical forms are `gauge_mult` times the functor's raw ones.  The two
            //  ratio forms are unaffected by that factor, so only the denominator their
            //  derivatives divide by needs it; the RelativeActivity form does not use this at all
            //  (its derivative is just 1).
            const double denom_reported = gauge_mult * denom_nominal;

            // d(g)/d(A_i): mass fraction g = m_t/M with m_i = A_i/s_i gives (delta_it - g)/(s_i*M);
            //  activity ratio g = A_n/A_d gives 1/A_d for the numerator and -g/A_d for the
            //  denominator; a plain relative activity is just the identity.
            vector<double> dg( num_nuclides, 0.0 );
            for( size_t i = 0; i < num_nuclides; ++i )
            {
              const string &iso = cost_functor->m_isotopes[i];
              if( is_mass_frac )
              {
                const map<string,double>::const_iterator pos = target.m_specific_activities.find( iso );
                if( pos != end(target.m_specific_activities) )
                  dg[i] = ((iso == target.m_nuclide) ? (1.0 - result.nominal_value) : -result.nominal_value)
                          / (pos->second * denom_reported);
              }else if( iso == target.m_nuclide )
              {
                dg[i] = is_rel_act ? 1.0 : (1.0 / denom_reported);
              }else if( !is_rel_act && (iso == target.m_denom_nuclide) )
              {
                dg[i] = -result.nominal_value / denom_reported;
              }
            }//for( size_t i = 0; i < num_nuclides; ++i )

            double var = 0.0;
            for( size_t i = 0; i < num_nuclides; ++i )
              for( size_t j = 0; j < num_nuclides; ++j )
                var += dg[i] * dg[j] * solution.m_rel_act_covariance[i][j];

            if( (var > 0.0) && !isinf(var) )
              sigma_raw = std::sqrt( var / chi2_scale );
          }//if( we have the relative-activity covariance )

          if( !(sigma_raw > 0.0) || isinf(sigma_raw) )
            sigma_raw = is_mass_frac ? (std::max)( 0.02, 0.25*result.nominal_value )
                                     : 0.25*result.nominal_value;

          const double max_delta_chi2 = chi2_scale * max_quantile;
          const double s_max = std::sqrt( max_delta_chi2 ); //target displacement, in sigma_raw units

          // Weight: with `w = K/sigma`, a request `t` away from nominal lands at `K^2/(1+K^2)` of
          //  the way there (from `grad C = -grad penalty` with a locally quadratic C), so K = 4
          //  tracks requests to ~94% while adding only ~K^2 to the curvature in that one direction.
          const double weight_k = 4.0;

          vector<pair<double,double>> scan; //(achieved value, peaks-only chi2)
          scan.emplace_back( result.nominal_value, result.nominal_chi2 );

          // Why each side stopped short of the widest threshold, if it did; used to explain the
          //  `*_at_bound` ends below.  Indexed [0] = lower, [1] = upper.
          string stop_reasons[2];

          for( const int direction : { -1, 1 } )
          {
            //each side starts from the nominal basin
            std::copy( begin(nominal_parameters), end(nominal_parameters), begin(parameters) );
            double sigma_est = sigma_raw, prev_value = result.nominal_value, prev_s = 0.0;
            string &stop_reason = stop_reasons[(direction < 0) ? 0 : 1];
            stop_reason = "the allowed number of re-fits ran out";

            // A step the solver cannot use - it did not converge or went non-finite - is RETRIED
            //  by bisecting toward the last accepted request rather than ending the side.
            //  Abandoning instead throws away the whole remaining range: measured on a
            //  Br82/K42/Na24 spectrum, one 0.85-sigma overshoot ended a walk at 2.64 when the
            //  data constrain the quantity near 0.5.
            // Requests are anchored on the previous REQUEST, not on the previous achieved value,
            //  and a rejected request is retried by bisecting toward the last accepted one.  The
            //  achieved value lags its request by an amount comparable to the whole remaining
            //  range near a domain edge, so a ladder anchored on it can step BACKWARDS after a
            //  retry - which then reads as the profile going flat when nothing of the sort
            //  happened.  Bisection is monotone by construction.
            double good_request = result.nominal_value;
            double bad_request = 0.0;
            bool have_bad = false;

            // Retries are counted separately from productive steps: a retry buys no ground, so
            //  charging it to the step budget would let a couple of overshoots end the walk short
            //  of the range the data actually constrain.  The second counter is only a
            //  runaway-loop backstop.
            size_t retries = 0;
            const size_t max_retries = 2*input.profile_max_solves_per_side;

            for( size_t step = 0; (step < input.profile_max_solves_per_side) && (retries <= max_retries); )
            {
              // Step outward from where the LAST solve landed, not from the nominal.  `sigma_est`
              //  is revised as the walk proceeds, and a ladder built on the nominal plus a
              //  multiple of the revised sigma can ask for a point INSIDE the previous one (any
              //  shrink of more than 2x does it) - which then looks exactly like the profile
              //  going flat, when all that moved was the ladder.
              //
              //  Deliberately NOT clamped to the quantity's physical range: the target is only a
              //  knob - the achieved value is what gets recorded - and clamping caps how hard the
              //  penalty can pull, so the walk stalls short of the real crossing and misreports it
              //  as a bound.  The bound is enforced where it belongs, by the parameter bounds.
              const double requested = have_bad
                          ? 0.5*(good_request + bad_request)
                          : (good_request + direction * 0.4 * s_max * sigma_est);

              // The bracket has closed without finding a usable point.
              if( have_bad && (fabs(bad_request - good_request) < 1.0E-3*fabs(sigma_est)) )
                break;

              cost_functor->arm_profile( static_cast<int>(target_index), requested,
                                         weight_k / sigma_est, denom_floor );

              ceres::Solver::Summary profile_summary;
              ceres::Solve( profile_options, &problem, &profile_summary );

              // Only a CONVERGED penalized solve is a point of the profile curve: the whole
              //  argument above is a statement about optima, and a solve that merely ran out of
              //  iterations or time sits ABOVE the curve at its achieved g - which would pull the
              //  interpolated crossing inward and report a too-narrow interval as if it were a real
              //  one.  (The nominal solve treats NO_CONVERGENCE as an error for the same reason.)
              if( (profile_summary.termination_type != ceres::CONVERGENCE)
                  && (profile_summary.termination_type != ceres::USER_SUCCESS) )
              {
                if( debug_printout() )
                  cout << "Profile " << target.m_nuclide << " step " << step << " dir " << direction
                       << ": " << profile_summary.BriefReport() << endl;

                bad_request = requested;
                have_bad = true;
                retries += 1;
                stop_reason = "a re-fit did not converge";
                continue;
                break;
              }

              const double achieved = cost_functor->profile_quantity( target, parameters, denom_floor );
              const double chi2 = peaks_only_chi2( parameters );

              if( isnan(achieved) || isinf(achieved) || isnan(chi2) || isinf(chi2) )
              {
                if( debug_printout() )
                  cout << "Profile " << target.m_nuclide << " step " << step << " dir " << direction
                       << ": non-finite at requested " << requested << endl;

                bad_request = requested;
                have_bad = true;
                retries += 1;
                stop_reason = "a re-fit gave a non-finite result";
                continue;
              }

              // For the Ceres-fit empirical forms the reported activity carries the
              //  average-measured-rel-eff multiple, and that average divides by each peak's source
              //  counts.  Exactly, this makes a reported activity unable to fall below
              //  `(1/P)*sum over peaks owned by t of (C_p/y_p)` - the `A_t` cancels.  Only the
              //  smooth floor inside `average_measured_rel_eff` lets it through that wall, and
              //  then the end-point is decided by the floor constant rather than by the data.
              //  Stop instead: the floor exists to keep the fit evaluable, not to produce answers.
              if( is_rel_act && cost_functor->gauge_normalizes_activities() )
              {
                bool any_peak_floored = false;
                cost_functor->average_measured_rel_eff( parameters, &any_peak_floored );
                if( any_peak_floored )
                {
                  if( debug_printout() )
                    cout << "Profile " << target.m_nuclide << " step " << step << " dir " << direction
                         << ": at the structural floor, requested " << requested << endl;

                  // A request below the structural floor has no minimum to find: the optimizer
                  //  just drives the activity to zero, where the smooth floor - not the data -
                  //  decides the answer.  Stop; the zero-activity anchor solve below is the
                  //  correct end-point for this side, and it is exact.
                  stop_reason = "the reported activity reached its structural floor, where this"
                                " nuclide's peaks are entirely unexplained";
                  break;
                }
              }//if( the reported activity carries the gauge multiple )

              scan.emplace_back( achieved, chi2 );

              if( debug_printout() )
                cout << "Profile " << target.m_nuclide << " step " << step << " dir " << direction
                     << ": requested " << requested << ", achieved " << achieved
                     << ", chi2 " << chi2 << " (nominal " << result.nominal_chi2 << ")" << endl;

              const double delta_chi2 = chi2 - result.nominal_chi2;

              // Calibrate the step scale from THIS step before diagnosing anything with it.  For a
              //  locally quadratic C the penalized optimum lands at `K^2/(K^2 + r^2)` of the way to
              //  the request, with `r = sigma_est/sigma_true` - so the achieved fraction measures
              //  `r` directly, and one badly seeded step can be corrected instead of being read as
              //  a verdict.  Without this, both diagnoses below fire on step 0 for a seed off by
              //  more than ~100x (which the covariance-free fallback seed can be) and report a
              //  confident "unbounded" or "at a bound" that is purely an artifact of the seed.
              //  Measured against a real step: predicted 16/17 = 0.9412, achieved 0.9409.
              const double asked = requested - prev_value;
              if( fabs(asked) > 0.0 )
              {
                const double frac = (achieved - prev_value) / asked;
                if( (frac > 1.0E-6) && (frac < 1.0) )
                {
                  const double r = weight_k * std::sqrt( 1.0/frac - 1.0 );
                  if( (r > 0.0) && !isinf(r) )
                    sigma_est = (std::max)( 0.02*sigma_est, (std::min)( 50.0*sigma_est, sigma_est/r ) );
                }
              }

              // The quantity stopped responding: it is against a bound (activity at zero, a mass
              //  fraction at 0/1, a constraint-window edge), so no stronger pull will move it.
              //  This is NOT the flat case below - the fit is still resisting, something else is
              //  simply holding the quantity - so it is diagnosed separately and first.
              if( fabs(achieved - prev_value) < 1.0E-3*sigma_est )
              {
                stop_reason = "it ran into a bound before the chi2 rose that far";
                break;
              }

              // Has the profile gone flat?  Each step aims to buy a fixed amount of sqrt(chi2)
              //  - `0.4*s_max` - so a step that buys almost none of it means the likelihood has
              //  stopped resisting and the quantity is simply running away.  A "crossing"
              //  interpolated from points out there is decided by numerical noise rather than by
              //  the data: measured on a Br82/K42/Na24 spectrum the last two points differed by
              //  0.05 in chi2 across a 2.8x change in activity, so an asymptote at 4.02 reported
              //  a finite 2-sigma bound where one at 3.99 would have reported none.  (Healthy
              //  profiles are nowhere near this: on that same spectrum the K42 and Na24 walks
              //  bought 57-94% of the intended amount every step, against 4% for flat Br82.)
              //  Skipped when the fit got BETTER than nominal - `s` is pinned at 0 there, which
              //  would read as flat while what actually happened is the re-anchoring case below.
              const double s_now = std::sqrt( (std::max)( 0.0, delta_chi2 ) );
              const bool went_flat = (delta_chi2 > 0.0)
                                     && ((s_now - prev_s) < 0.10*0.4*s_max);
              prev_s = s_now;

              // The flat verdict is taken BEFORE "we reached the threshold", not after.  The whole
              //  point is the knife-edge: an asymptote settling at 4.02 would otherwise be reported
              //  as a finite bound while one settling at 3.99 is reported as unbounded, which is a
              //  1% difference in the data deciding the entire character of the answer.  A step
              //  that crosses the threshold while buying almost none of the sqrt(chi2) it aimed for
              //  has crossed on noise, so it is a limit either way.  A healthy crossing is nowhere
              //  near this: measured, the well-behaved walks bought 38-94% of the aim on the step
              //  that crossed, against 4% for a genuinely flat direction.
              if( went_flat )
              {
                stop_reason = "the fit stopped resisting further change, so it is unbounded on"
                              " this side";
                break;
              }

              if( delta_chi2 >= max_delta_chi2 )
              {
                stop_reason.clear();  //this side is properly bracketed
                break;
              }

              // Also fold in the chi2 actually bought, which measures the sigma of the profile
              //  itself rather than the local response to the penalty; the two agree for a
              //  quadratic profile and this one keeps the walk moving on a flattening one.
              if( delta_chi2 > 0.05 )
              {
                const double measured = fabs(achieved - result.nominal_value) / std::sqrt(delta_chi2);
                if( (measured > 0.0) && !isinf(measured) )
                  sigma_est = (std::max)( 0.2*sigma_est, (std::min)( 5.0*sigma_est, measured ) );
              }

              prev_value = achieved;
              good_request = requested;
              have_bad = false;
              step += 1;
            }//for( loop over steps out from the nominal )
          }//for( const int direction : { -1, 1 } )

          cost_functor->arm_profile( -1, 0.0, 0.0, 0.0 );
          std::copy( begin(nominal_parameters), end(nominal_parameters), begin(parameters) );

          // The lower end of a relative-activity profile IS the zero-activity hypothesis - "none of
          //  this nuclide" - and the walk cannot reach it.  With the reporting gauge the quantity
          //  has a hard floor there (`m(x)*A_t` tends to `(1/P)*sum over peaks owned by t of
          //  C_p/y_p` as `A_t` goes to zero), which the penalty can only approach; past it the
          //  smooth floor inside `average_measured_rel_eff` decides the answer instead of the data.
          //  So solve it directly: cap the activity just above zero and re-fit everything else.
          //  What comes back is an ordinary point of the profile curve, so the crossing logic
          //  needs no special case - and it is the point that decides whether "absent" is inside
          //  the interval, i.e. whether the nuclide is detected at that confidence level.
          //  (chi2 at the floor is essentially the nuclide's detection significance squared, so
          //  for a well-detected nuclide the anchor sits far outside every threshold and changes
          //  nothing.)
          //  Only where the reporting gauge applies.  Without it the profiled quantity is the raw
          //  activity, which has no hard floor - the walk reaches the soft one the LLS
          //  normalization imposes on its own, and forcing the activity to zero instead sails past
          //  it into the region where that normalization blows the chi2 up (measured: 47401
          //  against the 4.04 this nuclide's significance predicts).
          if( is_rel_act && cost_functor->gauge_normalizes_activities() )
          {
            const size_t act_index = cost_functor->iso_index( target.m_nuclide );
            const bool is_free_activity =
                  (act_index < cost_functor->m_rel_act_norms.size())
                  && (cost_functor->m_rel_act_norms[act_index] > 0.0)
                  && !std::count( begin(constant_parameters), end(constant_parameters),
                                  static_cast<int>(act_index) )
                  && (nominal_parameters[act_index] > 0.0);

            if( is_free_activity )
            {
              const double saved_upper = problem.GetParameterUpperBound( pars, static_cast<int>(act_index) );
              const double cap = 1.0E-4 * nominal_parameters[act_index];

              parameters[act_index] = cap;  //Ceres requires a feasible starting point
              problem.SetParameterUpperBound( pars, static_cast<int>(act_index), cap );

              ceres::Solver::Summary anchor_summary;
              ceres::Solve( profile_options, &problem, &anchor_summary );

              problem.SetParameterUpperBound( pars, static_cast<int>(act_index), saved_upper );

              const bool anchor_ok = (anchor_summary.termination_type == ceres::CONVERGENCE)
                                     || (anchor_summary.termination_type == ceres::USER_SUCCESS);
              const double anchor_value = anchor_ok
                        ? cost_functor->profile_quantity( target, parameters, denom_floor ) : 0.0;
              const double anchor_chi2 = anchor_ok ? peaks_only_chi2( parameters ) : 0.0;

              if( anchor_ok && !isnan(anchor_value) && !isinf(anchor_value)
                  && !isnan(anchor_chi2) && !isinf(anchor_chi2)
                  && (anchor_value < result.nominal_value) )
              {
                scan.emplace_back( anchor_value, anchor_chi2 );
                stop_reasons[0] = "the lower end is the zero-activity hypothesis: the data do not"
                        " exclude this nuclide being absent at this level (with an empirical"
                        " equation the reporting normalization expresses 'absent' as that value,"
                        " not as 0)";

                if( debug_printout() )
                  cout << "Profile " << target.m_nuclide << " zero-activity anchor: value "
                       << anchor_value << ", chi2 " << anchor_chi2
                       << " (nominal " << result.nominal_chi2 << ")" << endl;
              }

              std::copy( begin(nominal_parameters), end(nominal_parameters), begin(parameters) );
            }//if( is_free_activity )
          }//if( is_rel_act )

          // Nothing moved the quantity at all, on either side: it is not a free function of the
          //  fitted parameters (an act-ratio chain that ties the whole roster together, an element
          //  whose other isotopes are all pinned, a zero-width constraint window, ...).  Reporting a
          //  zero-width "(limit)" interval at the nominal would read as a confident answer, so say
          //  what happened instead and leave the covariance-based uncertainty as the answer.
          //  Note the stalled point IS recorded before the walk gives up, so this cannot be a test
          //  on how many points there are - it has to be on how far any of them got.
          double largest_move = 0.0;
          for( const pair<double,double> &pt : scan )
            largest_move = (std::max)( largest_move, fabs(pt.first - result.nominal_value) );

          if( largest_move <= 1.0E-6*result.nominal_value )
            throw runtime_error( "the quantity is fixed by the problem's constraints, so there is"
                                 " nothing to profile" );

          const pair<double,double> nominal_point( result.nominal_value, result.nominal_chi2 );

          std::sort( begin(scan), end(scan),
                     []( const pair<double,double> &l, const pair<double,double> &r ){
                       return l.first < r.first;
                     } );
          result.scan_points = scan;

          // If a scan point fit the data better than the nominal solve did, the nominal was a local
          //  minimum; reference the interval to the better point.  The reference is always the scan
          //  minimum (so `nominal_chi2` really is what the delta-chi2 is measured from, as its
          //  documentation says); only the warning is gated on the improvement being meaningful.
          double ref_chi2 = result.nominal_chi2;
          for( const pair<double,double> &pt : scan )
            ref_chi2 = (std::min)( ref_chi2, pt.second );

          if( ref_chi2 < (result.nominal_chi2 - 0.5) )
            result.warnings.push_back( "The profile scan found a better fit (chi2 "
                + SpecUtils::printCompact(ref_chi2, 5) + " vs nominal "
                + SpecUtils::printCompact(result.nominal_chi2, 5) + ") - the nominal solution may be"
                " a local minimum; the profile interval is referenced to the scan minimum." );
          result.nominal_chi2 = ref_chi2;

          // Interval end-points: the crossings of sqrt(chi2 - ref_chi2) = sqrt(delta_chi2), found by
          //  linear interpolation in (value, sqrt(delta chi2)) - which is exact for a quadratic
          //  likelihood, so a handful of points suffices.
          for( size_t cl_index = 0; cl_index < input.profile_confidence_levels.size(); ++cl_index )
          {
            if( quantiles[cl_index] <= 0.0 )
              continue;

            ProfileInterval interval;
            interval.confidence_level = input.profile_confidence_levels[cl_index];
            interval.delta_chi2 = chi2_scale * quantiles[cl_index];

            const double s_target = std::sqrt( interval.delta_chi2 );

            bool degenerate = false;

            for( const int direction : { -1, 1 } )
            {
              // The scan points on this side, ordered from the nominal outward.  The nominal is
              //  prepended rather than selected by value, so a stalled step that landed exactly on
              //  it cannot end up ahead of it (`scan` is sorted on the value alone).
              vector<pair<double,double>> side{ nominal_point };
              for( const pair<double,double> &pt : scan )
              {
                if( (direction < 0) ? (pt.first < result.nominal_value)
                                    : (pt.first > result.nominal_value) )
                  side.push_back( pt );
              }
              if( direction < 0 )
                std::reverse( begin(side) + 1, end(side) ); //keep the nominal first

              double crossing = result.nominal_value;
              bool crossed = false;

              // The nominal point already past the threshold means the scan re-anchored to a better
              //  minimum elsewhere, and there is no interval around the nominal to report.
              if( std::sqrt( (std::max)( 0.0, side[0].second - ref_chi2 ) ) >= s_target )
              {
                degenerate = true;
              }else for( size_t k = 0; !crossed && ((k + 1) < side.size()); ++k )
              {
                const double s_i = std::sqrt( (std::max)( 0.0, side[k].second - ref_chi2 ) );
                const double s_j = std::sqrt( (std::max)( 0.0, side[k+1].second - ref_chi2 ) );

                if( (s_i <= s_target) && (s_j >= s_target) )
                {
                  crossing = (s_j > s_i)
                       ? (side[k].first + (side[k+1].first - side[k].first)*(s_target - s_i)/(s_j - s_i))
                       : side[k+1].first;
                  crossed = true;
                }
              }//for( walk outward looking for a bracketing pair )

              if( !crossed && !degenerate )
                crossing = side.back().first; //furthest the quantity actually got

              if( direction < 0 )
              {
                interval.lower_value = crossing;
                interval.lower_at_bound = !crossed;
              }else
              {
                interval.upper_value = crossing;
                interval.upper_at_bound = !crossed;
              }
            }//for( const int direction : { -1, 1 } )

            if( degenerate )
              result.warnings.push_back( "The value the fit landed on is itself excluded at the "
                      + SpecUtils::printCompact(100.0*interval.confidence_level, 4) + "% level by the"
                      " profile scan (a better fit was found elsewhere), so no interval could be"
                      " formed around it." );

            // Explain each end that came back a limit, at the confidence level it actually applies
            //  to.  Emitting this from the walk instead would over-claim: a side can stop short of
            //  the 2-sigma threshold having crossed the 1-sigma one perfectly well, and a blanket
            //  "this side is unbounded" would then contradict the 1-sigma interval beside it.
            if( !degenerate )
            {
              const bool ends_at_bound[2] = { interval.lower_at_bound, interval.upper_at_bound };
              for( size_t side = 0; side < 2; ++side )
              {
                if( ends_at_bound[side] && !stop_reasons[side].empty() )
                  result.warnings.push_back( string("The ") + (side ? "upper" : "lower") + " "
                          + SpecUtils::printCompact(100.0*interval.confidence_level, 4)
                          + "% end for " + target.m_nuclide + " is a limit, not a measured bound: "
                          + stop_reasons[side] + "." );
              }
            }//if( !degenerate )

            result.intervals.push_back( interval );
          }//for( loop over confidence levels )
        }catch( std::exception &e )
        {
          cost_functor->arm_profile( -1, 0.0, 0.0, 0.0 );
          std::copy( begin(nominal_parameters), end(nominal_parameters), begin(parameters) );
          result.intervals.clear();
          solution.m_warnings.push_back( "Profile-likelihood scan for " + target.m_nuclide
                                         + " failed: " + string(e.what()) );
          continue;
        }//try / catch around one profile target

        solution.m_profile_results.push_back( std::move(result) );
      }//for( loop over profile targets )
    }//else (we have at least one valid confidence level)
  }//if( profile-likelihood intervals were asked for )


  /*
  if( (solution.m_status == ManualSolutionStatus::Success)
     && !solution.m_nonlin_covariance.empty() && !solution.m_fit_parameters.empty() )
  {
    cout << "Covariance matrix:" << endl;
    cout << "  " << setw(10) << " " << " ";
    for( size_t row = 0; row < num_parameters; ++row )
      cout << setw(10) << solution.parameter_name(row) << " ";
    cout << endl;
    
    for( size_t row = 0; row < num_parameters; ++row )
    {
      cout << "  " << setw(10) << solution.parameter_name(row) << " ";
      for( size_t col = 0; col < num_parameters; ++col )
        cout << setw(10) << setprecision(2) << solution.m_nonlin_covariance[row][col] << " ";
      cout << endl;
    }//for( size_t row = 0; row < num_nuclides; ++row )
  }//if( we have the covariance matrix )
  */
  
  return solution;
}//solve_relative_efficiency(...)


double estimate_stat_uncert_multiple( const RelEffSolution &solution, const double max_multiple )
{
  if( solution.m_status != ManualSolutionStatus::Success )
    return -1.0;

  const vector<GenericPeakInfo> &peaks = solution.m_input.peaks;
  if( peaks.size() < 4 )  //a median over fewer than this is not meaningful
    return -1.0;

  vector<double> sq_pulls;
  for( const GenericPeakInfo &peak : peaks )
  {
    if( peak.m_base_rel_eff_uncert < -1.0E-9 ) //unweighted fit - no statistical scale to work from
      return -1.0;

    if( (peak.m_counts <= 0.0) || (peak.m_counts_uncert <= 0.0) )
      return -1.0;

    double src_counts = 0.0;
    for( const GenericLineInfo &line : peak.m_source_gammas )
      src_counts += solution.relative_activity( line.m_isotope ) * line.m_yield;

    const double predicted = src_counts * solution.rel_eff_eqn_value( peak.m_energy );

    // Deliberately against the peak's STATISTICAL uncertainty: `k` is the multiple of that
    //  quantity we are solving for, whatever else may have been folded into the fit weights.
    const double pull = (predicted - peak.m_counts) / peak.m_counts_uncert;
    sq_pulls.push_back( pull * pull );
  }//for( const GenericPeakInfo &peak : peaks )

  const size_t num_peaks = sq_pulls.size();
  std::sort( begin(sq_pulls), end(sq_pulls) );
  const double median_sq = (num_peaks % 2) ? sq_pulls[num_peaks/2]
                                    : 0.5*(sq_pulls[num_peaks/2 - 1] + sq_pulls[num_peaks/2]);

  const double chi2_1_median = 0.4549364; //median of a chi-squared(1) distribution
  const double leverage_corr = static_cast<double>(num_peaks) / (std::max)( solution.m_dof, 1 );
  const double variance_multiple = median_sq * leverage_corr / chi2_1_median;

  return (std::min)( max_multiple, (std::max)( 1.0, std::sqrt(variance_multiple) ) );
}//double estimate_stat_uncert_multiple( const RelEffSolution &solution, const double max_multiple )






namespace PeakCsvInput
{
NuclideInfo::NuclideInfo( const char *p, const char *nuc, bool opt, float kev, float br )
: parent(p), source_nuclide(nuc), energy(kev), yield(br), optional(opt)
{
}

const char *to_str( const NucDataSrc src )
{
  switch( src )
  {
    case NucDataSrc::Icrp107_U:        return "Icrp107_U";
    case NucDataSrc::Lanl_U:           return "Lanl_U";
    case NucDataSrc::IcrpLanlGadras_U: return "IcrpLanlGadras_U";
    case NucDataSrc::SandiaDecay:      return "SandiaDecay";
    case NucDataSrc::Undefined:        return "Undefined";
  }
  assert( 0 );
  return "";
}//to_str( NucDataSrc )



NucMatchResults fill_in_nuclide_info( const vector<RelActCalcManual::GenericPeakInfo> peaks,
                                     const NucDataSrc nuc_data_src,
                                     const vector<pair<float,float>> energy_ranges,
                                     std::vector<NucAndAge> isotopes,
                                     const float energy_tolerance_sigma,
                                     const vector<float> excluded_peak_energies,
                                     const float measurement_duration )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  //https://journals.sagepub.com/doi/pdf/10.1177/ANIB_38_3
  // and more specifically the supplementary material at
  // https://journals.sagepub.com/doi/suppl/10.1177/ANIB_38_3
  const vector<NuclideInfo> icrp107{{//18 entries
    {"U235", "U235",         false, 143.76f,      0.11f},
    {"U235", "U235",         false, 163.33f,      0.0508f},
    {"U235", "U235",         false, 185.715f,     0.572f},
    {"U235", "U235",         true,  202.11f,      0.0108f},
    {"U235", "U235",         false, 205.311f,     0.0501f},
    
    {"U232", "Pb212",        false, 238.632f,     0.433f},
    {"U232", "Tl208",        false, 583.191f,     0.304f},
    {"U232", "Bi212",        false, 727.33f,      0.0658f},
    {"U232", "Tl208",        false, 860.564f,     0.0447f},
    
    {"U238", "Pa234m",       false, 258.26f,      0.000726833f},
    {"U238", "Pa234",        true,  569.3173913f, 0.00018952f},
    {"U238", "Pa234/Pa234m", true,  742.81f,      0.000831678f},
    {"U238", "Pa234m",       false, 766.36f,      0.002935286f},
    {"U238", "Pa234/Pa234m", true,  880.5742495f, 0.000204387f},
    {"U238", "Pa234/Pa234m", true,  883.24f,      0.000188208f},
    {"U238", "Pa234/Pa234m", true,  945.9684295f, 0.000313081f},
    {"U238", "Pa234m",       false, 1001.03f,     0.008356588f},
    
    {"U234", "U234",         false, 120.9f,       0.000397362f},
  }};//icrp107
  
  
  const vector<NuclideInfo> lanl{{ //17 entries
    {"U235", "U235",          false, 143.76f,  0.11f},
    {"U235", "U235",          false, 163.36f,  0.0505f},
    {"U235", "U235",          false, 185.715f, 0.57f},
    {"U235", "U235",          true,  202.11f,  0.01098f},
    {"U235", "U235",          false, 205.311f, 0.0503f},
    
    {"U232", "Pb212",         false, 238.625f, 0.48f},
    {"U232", "Tl208",         false, 583.187f, 0.306f},
    {"U232", "Bi212",         false, 727.3f,   0.0676f},
    {"U232", "Tl208",         false, 860.56f,  0.046f},
    
    {"U238", "Pa234m",        false, 258.26f,  0.000754f},
    {"U238", "Pa234/Pa234m",  true,  742.83f,  0.000907f},
    {"U238", "Pa234m",        false, 766.4f,   0.003074f},
    {"U238", "Pa234/Pa234m",  true,  880.47f,  0.000213f},
    {"U238", "Pa234/Pa234m",  true,  883.24f,  0.000213f},
    {"U238", "Pa234/Pa234m",  true,  945.95f,  0.000347f},
    {"U238", "Pa234m",        false, 1001.03f, 0.008371f},
    
    {"U234", "U234",          false, 120.905f, 0.00035f}
  }};
  
  
  const vector<NuclideInfo> lanl_icrp_gad{{ //21 entries
    {"U235", "U235",          false, 143.76f,      0.11f},
    {"U235", "U235",          false, 163.36f,      0.0505f},
    {"U235", "U235",          false, 185.715f,     0.57f},
    {"U235", "U235",          true,  202.11f,      0.01098f},
    {"U235", "U235",          false, 205.311f,     0.0503f},
    {"U235", "U235",          true,  221.38f,      0.0012f},
    {"U235", "U235",          true,  246.84f,      0.00053f},
    {"U235", "U235",          true,  345.9f,       0.0003f},
    
    {"U232", "Pb212",         false, 238.625f,     0.48f},
    {"U232", "Tl208",         false, 583.187f,     0.306f},
    {"U232", "Bi212",         false, 727.3f,       0.0676f},
    {"U232", "Tl208",         false, 860.56f,      0.046f},
    
    {"U238", "Pa234m",        false, 258.26f,      0.000754f},
    {"U238", "Pa234",         true,  569.3173913f, 0.00018952f},
    {"U238", "Pa234/Pa-234m", true,  742.83f,      0.000907f},
    {"U238", "Pa234m",        false, 766.4f,       0.003074f},
    {"U238", "Pa234/Pa-234m", true,  880.47f,      0.000213f},
    {"U238", "Pa234/Pa-234m", true,  883.24f,      0.000213f},
    {"U238", "Pa234/Pa-234m", true,  945.95f,      0.000347f},
    {"U238", "Pa234m",        false, 1001.03f,     0.008371f},
    
    {"U234", "U234",          false, 120.905f,     0.00035f}
  }};
  
  // Check isotopes are valid, and normalize their name, and if there is an age, get it.
  if( isotopes.empty() )
  {
    if( (nuc_data_src == NucDataSrc::SandiaDecay) || (nuc_data_src == NucDataSrc::Undefined) )
      throw runtime_error( "fill_in_nuclide_info: No nuclides, or specialized nuclear data sources specified." );
       
    vector<NuclideInfo> specialize_src;
    switch( nuc_data_src )
    {
      case NucDataSrc::Icrp107_U:
        specialize_src = icrp107;
        break;
        
      case NucDataSrc::Lanl_U:
        specialize_src = lanl;
        break;
        
      case NucDataSrc::IcrpLanlGadras_U:
        specialize_src = lanl_icrp_gad;
        break;
        
      case NucDataSrc::SandiaDecay:
      case NucDataSrc::Undefined:
        assert( 0 );
        break;
    }//switch( nuc_data_src )
    
    
    for( const auto &info : specialize_src )
    {
      bool has_iso = false;
      for( const auto &iso : isotopes )
        has_iso |= (iso.nuclide == info.parent);
      
      if( !has_iso )
      {
        const SandiaDecay::Nuclide *nuc = db->nuclide(info.parent);
        assert( nuc );
        if( !nuc )
          throw runtime_error( "Some how '" + info.parent + "' isnt a valid nuclide." );
        
        isotopes.emplace_back( nuc->symbol, -1.0, false );
      }//
    }//for( const auto &info : specialize_src )
  }//if( isotopes.empty() )
  
  
  // We will put raw source info into 'initial_nucs_info', and then filter this down to the
  //  nuclides and energy ranges we will actually use
  //  If we are correcting activities for nuclides decay during the measurement, then we will
  //  place the uncorrected info into `non_decay_corr_initial_nucs_info` (that is,
  //  `non_decay_corr_initial_nucs_info` will not have any entries if we arent decay correcting)
  vector<NuclideInfo> initial_nucs_info, non_decay_corr_initial_nucs_info;
  for( size_t i = 0; i < isotopes.size(); ++i )
  {
    string &iso = isotopes[i].nuclide;
    double &age = isotopes[i].age;
    
    const ReactionGamma::Reaction *rctn = nullptr;
    const SandiaDecay::Nuclide *nuc = db->nuclide( iso );
    const SandiaDecay::Element *element = !nuc ? db->element( iso ) : nullptr;
    if( !nuc && !element )
    {
      const ReactionGamma *reactiondb = ReactionGammaServer::database();
      
      try
      {
        vector<ReactionGamma::ReactionPhotopeak> reactions;
        reactiondb->gammas( iso, reactions );
        if( reactions.empty() )
          throw runtime_error( "unknown reaction" );
        
        // TODO: currently just using first possible reaction; need to implement retrieving reactions by name from ReactionGamma.
        rctn = reactions[0].reaction;
      }catch( std::exception &e )
      {
        cerr << "Invalid reaction ("<< iso << "): " << e.what() << endl;
      }
    }//if( !nuc )
    
    if( !nuc && !rctn && !element )
      throw runtime_error( "Invalid nuclide '" + iso + "' specified." );
    
    iso = nuc ? nuc->symbol : (element ? element->symbol : rctn->name());
    
    // Make sure we try to use a U-data source, only for U, and the nuclides of it we have,
    //  otherwise we'll default back to using SandiaDecay
    NucDataSrc src = nuc_data_src;
    if( (nuc && (nuc->atomicNumber != 92)) || rctn || element )
    {
      src = NucDataSrc::SandiaDecay;
    }else if( nuc->atomicNumber == 92 )
    {
      auto has_U_nuc = [db]( const SandiaDecay::Nuclide *nuc, const vector<NuclideInfo> &input ) -> bool {
        for( const NuclideInfo &inputnuc : input )
        {
          const SandiaDecay::Nuclide *inputnuc_ptr = db->nuclide( inputnuc.parent );
          assert( inputnuc_ptr );
          if( inputnuc_ptr == nuc )
            return true;
        }//for( const NuclideInfo &inputnuc : input )
        return false;
      };//filter_for_nuc(...)
      
      // TODO: we should probably give a warning that we are using SandiaDecay data if the
      //       uranium dataset doesn't have this nuclide.
      switch( src )
      {
        case NucDataSrc::Icrp107_U:
          if( !has_U_nuc(nuc, icrp107) )
            src = NucDataSrc::SandiaDecay;
        break;
          
        case NucDataSrc::Lanl_U:
          if( !has_U_nuc(nuc, lanl) )
            src = NucDataSrc::SandiaDecay;
        break;
          
        case NucDataSrc::IcrpLanlGadras_U:
          if( !has_U_nuc(nuc, lanl_icrp_gad) )
            src = NucDataSrc::SandiaDecay;
        break;
          
        case NucDataSrc::SandiaDecay:
        case NucDataSrc::Undefined:
          break;
      }//switch( src )
    }//if( non-U nuclide ) / else
    
    
    switch( src )
    {
      case NucDataSrc::Icrp107_U:
      case NucDataSrc::Lanl_U:
      case NucDataSrc::IcrpLanlGadras_U:
      case NucDataSrc::Undefined:
        age = -1.0;
        break;
        
      case NucDataSrc::SandiaDecay:
        if( nuc && (age < 0.0) )
          age = PeakDef::defaultDecayTime(nuc, nullptr);
        break;
    }//switch( nucdatasrc )
    
    
    auto filter_for_nuc = [db]( const SandiaDecay::Nuclide *nuc, const vector<NuclideInfo> &input ) -> vector<NuclideInfo> {
      vector<NuclideInfo> results;
      
      for( const NuclideInfo &inputnuc : input )
      {
        const SandiaDecay::Nuclide *inputnuc_ptr = db->nuclide( inputnuc.parent );
        assert( inputnuc_ptr );
        if( inputnuc_ptr == nuc )
          results.push_back( inputnuc );
      }//for( const NuclideInfo &inputnuc : input )
      
      return results;
    };//filter_for_nuc(...)
    
    
    switch( src )
    {
      case NucDataSrc::Icrp107_U:
      {
        const auto these_infos = filter_for_nuc( nuc, icrp107 );
        initial_nucs_info.insert( end(initial_nucs_info), begin(these_infos), end(these_infos) );
        break;
      }
        
      case NucDataSrc::Lanl_U:
      {
        const auto these_infos = filter_for_nuc( nuc, lanl );
        initial_nucs_info.insert( end(initial_nucs_info), begin(these_infos), end(these_infos) );
        break;
      }
        
      case NucDataSrc::IcrpLanlGadras_U:
      {
        const auto these_infos = filter_for_nuc( nuc, lanl_icrp_gad );
        initial_nucs_info.insert( end(initial_nucs_info), begin(these_infos), end(these_infos) );
        break;
      }
        
      case NucDataSrc::SandiaDecay:
      {
        if( isotopes.empty() )
          throw runtime_error( "You must specify the isotopes to use when using SandiaDecay as your source data." );
        
        if( rctn )
        {
          for( const auto &g : rctn->gammas )
          {
            if( g.abundance > std::numeric_limits<float>::min() )
              initial_nucs_info.emplace_back( rctn->name().c_str(), "", true, g.energy, g.abundance );
          }
        }else if( nuc )
        {
          assert( nuc );
          const double ref_act = 1.0*SandiaDecay::MBq;
          const double decrease_factor = std::exp( -age * nuc->decayConstant() );
          const double initial_activity = ref_act / decrease_factor;
          SandiaDecay::NuclideMixture mix;
          mix.addNuclideByActivity( nuc, initial_activity );
          
          assert( fabs(ref_act - mix.activity(age, nuc)) < 0.001*ref_act );
          
          vector<SandiaDecay::EnergyRatePair> un_decay_corrected_photons;
          vector<SandiaDecay::EnergyRatePair> photons = mix.photons( age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
          
          if( isotopes[i].decay_during_measurement )
          {
            vector<SandiaDecay::EnergyRatePair> corr_photons
                  = GammaInteractionCalc::decay_during_meas_corrected_gammas( mix, age, measurement_duration );
            
            // See what the decay correction is here, and note to then include in the results.
            assert( corr_photons.size() == photons.size() );
            assert( std::is_sorted(begin(corr_photons), end(corr_photons),
                   []( const SandiaDecay::EnergyRatePair &lhs, const SandiaDecay::EnergyRatePair &rhs ){
              return lhs.energy < rhs.energy;
            }) );
            assert( std::is_sorted(begin(photons), end(photons),
                   []( const SandiaDecay::EnergyRatePair &lhs, const SandiaDecay::EnergyRatePair &rhs ){
              return lhs.energy < rhs.energy;
            }) );
            //for( size_t i = 0; i < std::min(photons.size(),corr_photons.size()); ++i ){
            //  assert( fabs(photons[i].energy - corr_photons[i].energy) < 0.001 );
            //}
            
            un_decay_corrected_photons = std::move(photons);
            photons = std::move(corr_photons);
          }//if( isotopes[i].decay_during_measurement ) / else
          
          for( size_t info_index = 0; info_index < photons.size(); ++info_index )
          {
            const SandiaDecay::EnergyRatePair &rate_info = photons[info_index];
              
            const char * const parent = nuc->symbol.c_str();
            const char * const resp_nuc = ""; //TODO: bother to get the nuclide thats actually giving off this gamma; see decay_gammas(...) in RelActCalcAuto.cpp
            const bool optional_use = true;
            const float energy = rate_info.energy;
            const float br = static_cast<float>( rate_info.numPerSecond / ref_act );
            
            if( br > std::numeric_limits<float>::min() )
            {
              initial_nucs_info.emplace_back( parent, resp_nuc, optional_use, energy, br );
              
              if( info_index < un_decay_corrected_photons.size() )
              {
                assert( un_decay_corrected_photons.size() == photons.size() );
                const SandiaDecay::EnergyRatePair &noncorr_rate_info = un_decay_corrected_photons[info_index];
                assert( fabs(noncorr_rate_info.energy - energy) < 0.001 );
                
                const float non_corr_br = static_cast<float>( noncorr_rate_info.numPerSecond / ref_act );
                non_decay_corr_initial_nucs_info.emplace_back( parent, resp_nuc, optional_use, energy, non_corr_br );
              }//if( decay correct nuclides )
            }//if( br > std::numeric_limits<float>::min() )
          }//for( const SandiaDecay::EnergyRatePair info : photons )
        }else if( element )
        {
          for( const SandiaDecay::EnergyIntensityPair &g : element->xrays )
          {
            if( g.intensity > std::numeric_limits<float>::min() )
              initial_nucs_info.emplace_back( element->symbol.c_str(), "", true, g.energy, g.intensity );
          }
        }else
        {
          assert( 0 );
        }//if( rctn ) / else( nuc )
        
        
        break;
      }//case NucDataSrc::SandiaDecay:
        
      case NucDataSrc::Undefined:
        assert(0);
        throw std::logic_error( "shouldnt be here" );
        break;
    }//switch( nucdatasrc )
  }//for( size_t i = 0; i < isotopes.size(); ++i )
  
  
  std::sort( begin(initial_nucs_info), end(initial_nucs_info), []( const NuclideInfo &lhs, const NuclideInfo &rhs ){
    return lhs.energy < rhs.energy;
  });

  std::sort( begin(non_decay_corr_initial_nucs_info), end(non_decay_corr_initial_nucs_info),
            []( const NuclideInfo &lhs, const NuclideInfo &rhs ){
    return lhs.energy < rhs.energy;
  });
  
  
  // Now we'll put the entries from 'initial_nucs_info' we *might* actually use, into
  //  'candidate_nucs_info'.  Note that we're doing this in stages to provide warnings to the
  //  user about peaks not used, or matched, or whatever.
  vector<NuclideInfo> candidate_nucs_info;
  for( const auto &info : initial_nucs_info )
  {
    // Check if isotope is wanted
    bool nuc_wanted = false;
    for( const auto &n : isotopes )
      nuc_wanted |= (n.nuclide == info.parent);
    
    if( !nuc_wanted )
      continue;
    
    // Check if in energy range  (energy_ranges)
    bool in_energy_range = energy_ranges.empty();
    for( const pair<float,float> &r : energy_ranges )
      in_energy_range |= ((info.energy >= r.first) && (info.energy <= r.second));
    
    if( !in_energy_range )
      continue;
    
    candidate_nucs_info.push_back( info );
  }//for( const auto &info : initial_nucs_info )
  
  // Now go through and actually match up the info
  vector<RelActCalcManual::GenericPeakInfo> matched_peaks = peaks;
  vector<bool> used_isotope( isotopes.size(), false );
  vector<bool> used_peak( matched_peaks.size(), false );
  vector<bool> peak_was_excluded( matched_peaks.size(), false );
  vector<bool> used_candidate_nucs_info( candidate_nucs_info.size(), false );
  
  // Lets keep track of peaks with gammas that have received decay-during-measurement
  //  corrections; following variable will only have un-corrected values for peaks that
  //  received corrections and only entries in #GenericPeakInfo::m_source_gammas that
  //  have been corrected.
  vector<RelActCalcManual::GenericPeakInfo> un_corrected_peaks;
  
  for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  {
    RelActCalcManual::GenericPeakInfo &peak = matched_peaks[peak_index];
    
    const double peak_sigma = (peak.m_fwhm > 0.0) ? (peak.m_fwhm / 2.35482) : 1.0;
    
    // Check if this peak is specifically excluded via the 'exclude-peak' command line argument
    bool exclude = false;
    for( size_t i = 0; !exclude && (i < excluded_peak_energies.size()); ++i )
      exclude = (fabs(excluded_peak_energies[i] - peak.m_energy) <= (peak_sigma * energy_tolerance_sigma) );
    peak_was_excluded[peak_index] = exclude;
    if( exclude )
      continue;
    
    RelActCalcManual::GenericPeakInfo un_decay_corr_peak = peak;
    assert( un_decay_corr_peak.m_source_gammas.empty() );
    un_decay_corr_peak.m_source_gammas.clear();
    
    // Try to match this peak to a source gamma line
    for( size_t nuc_index = 0; nuc_index < candidate_nucs_info.size(); ++nuc_index )
    {
      NuclideInfo &nuc = candidate_nucs_info[nuc_index];
      
      const double peak_sigma = (peak.m_fwhm > 0.0) ? (peak.m_fwhm / 2.35482) : 1.0;
      
      if( fabs(nuc.energy - peak.m_energy) <= (peak_sigma * energy_tolerance_sigma) )
      {
        used_peak[peak_index] = true;
        used_candidate_nucs_info[nuc_index] = true;
        
        size_t iso_pos = isotopes.size();
        for( size_t i = 0; i < isotopes.size(); ++i )
        {
          if( isotopes[i].nuclide == nuc.parent )
          {
            iso_pos = i;
            break;
          }
        }

        assert( iso_pos != isotopes.size() );
        if( iso_pos == isotopes.size() )//Shouldnt ever happen
          throw std::logic_error( "Failed to find source nuclide in isotopes to use, after filtering" );
        
        used_isotope[iso_pos] = true;
        
        bool nuc_already_used = false;
        for( RelActCalcManual::GenericLineInfo &gamma : peak.m_source_gammas )
        {
          if( gamma.m_isotope == nuc.parent )
          {
            nuc_already_used = true;
            gamma.m_yield += nuc.yield;
            break;
          }
        }//for( loop over existing source gammas for `peak` )
        
        if( !nuc_already_used )
          peak.m_source_gammas.emplace_back( nuc.yield, nuc.parent );
        
        
        //look through `non_decay_corr_initial_nucs_info` and try to match to `isotopes[i]`
        //  TODO: `non_decay_corr_initial_nucs_info` is sorted by energy, so we could be much
        //         smarter here
        for( const NuclideInfo &un_corr_nuc : non_decay_corr_initial_nucs_info )
        {
          if( (nuc.parent == un_corr_nuc.parent)
             && (nuc.source_nuclide == un_corr_nuc.source_nuclide)
             && (fabs(nuc.energy - un_corr_nuc.energy) < 0.001) )
          {
            bool uncorr_nuc_already_used = false;
            for( RelActCalcManual::GenericLineInfo &gamma : un_decay_corr_peak.m_source_gammas )
            {
              if( gamma.m_isotope == nuc.parent )
              {
                uncorr_nuc_already_used = true;
                gamma.m_yield += un_corr_nuc.yield;
                break;
              }
            }//for( loop over existing source gammas for `peak` )
            
            if( !uncorr_nuc_already_used )
              un_decay_corr_peak.m_source_gammas.emplace_back( un_corr_nuc.yield, un_corr_nuc.parent );
          }//if( we found un-corrected line corresponding to `nuc` )
          
          if( un_corr_nuc.energy > nuc.energy )
            break;
        }//for( loop over non_decay_corr_initial_nucs_info )
      }//if( peak matched source data within tolerance )
    }//for( size_t nuc_index = 0; nuc_index < candidate_nucs_info.size(); ++candidate_nucs_info )
    
    if( !un_decay_corr_peak.m_source_gammas.empty() )
      un_corrected_peaks.push_back( un_decay_corr_peak );
  }//for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  
  
  NucMatchResults results;
  results.data_source = nuc_data_src;
  
  results.match_sigma_tolerance = energy_tolerance_sigma;
  results.energy_ranges = energy_ranges;
  results.not_decay_corrected_peaks = un_corrected_peaks;
  
  for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  {
    if( peak_was_excluded[peak_index] )
    {
      results.peaks_excluded.push_back( matched_peaks[peak_index] );
    }else
    {
      auto &place = used_peak[peak_index] ? results.peaks_matched : results.peaks_not_matched;
      place.push_back( matched_peaks[peak_index] );
    }
  }//for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  
  
  assert( used_isotope.size() == isotopes.size() );
  for( size_t index = 0; index < used_isotope.size(); ++index )
  {
    if( used_isotope[index] )
    {
      results.used_isotopes.push_back( isotopes[index].nuclide );
      results.used_isotope_ages.push_back( isotopes[index].age );
    }else
    {
      results.unused_isotopes.push_back( isotopes[index].nuclide );
    }
  }//for( size_t index = 0; index < used_isotope.size(); ++index )
  
  
  assert( used_candidate_nucs_info.size() == candidate_nucs_info.size() );
  for( size_t index = 0; index < candidate_nucs_info.size(); ++index )
  {
    auto &place = used_candidate_nucs_info[index] ? results.source_gammas_used : results.source_gammas_not_used;
    place.push_back( candidate_nucs_info[index] );
  }//for( size_t index = 0; index < candidate_nucs_info.size(); ++index )
  
  
  return results;
}//fill_in_nuclide_info(...)



vector<RelActCalcManual::GenericPeakInfo> peak_csv_to_peaks( istream &csv )
{
  //Gadras: "Energy(keV),sigma,Rate(cps),sigma,FWHM(keV),sigma,Leakage(1/s),Centroid,FileName,RecordIdx,Title,DateTime"
  //PeakEasy: "Centroid,  Net_Area,   Net_Area,      Peak, FWHM,   FWHM,Reduced, ROI_Total,ROI"
  //          "keV,    Counts,Uncertainty,       CPS,  keV,Percent,Chi_Sqr,    Counts,ID#,  File, LiveTime, Date, Time"
  
  enum class PeakCsvFormat{ PeakEasy, Gadras, Unknown };
  
  // If the CSV file is from InterSpec, and there are any nuclides provided, we will use them.
  bool contained_nuc_ids = false;
  PeakCsvFormat csv_format = PeakCsvFormat::Unknown;
  
  
  auto split_and_trim = []( vector<string> &fields, const string &line ){
    SpecUtils::split_no_delim_compress( fields, line, "," );
    for( string &s : fields )
      SpecUtils::trim( s );
  };
  
  string line;
  while( std::getline( csv, line ) )
  {
    SpecUtils::trim( line );
    if( line.empty() || line[0] == '#' )
      continue;
    
    // Line should either be "Energy(keV),sigma,Rate(cps)...", or "Centroid,  Net_Area,   Net_Area"
    vector<string> fields;
    split_and_trim( fields, line );
    if( fields.empty() )
      continue;
    
    if( fields.size() < 9 )
      throw runtime_error( "Invalid Peak CSV header line: '" + line + "'" );
    
    if( (fields[0] == "Energy(keV)") && (fields[1] == "sigma")
       && (fields[2] == "Rate(cps)") && (fields[3] == "sigma") )
    {
      csv_format = PeakCsvFormat::Gadras;
      break;
    }else if( (fields[0] == "Centroid") && (fields[1] == "Net_Area")
             && (fields[2] == "Net_Area") && (fields[3] == "Peak") )
    {
      if( !std::getline( csv, line ) )
        throw runtime_error( "Failed to get second line of PeakEasy CSV" );
      
      split_and_trim( fields, line );
      
      if( (fields.size() < 9) || (fields[0] != "keV") || (fields[1] != "Counts")
         || (fields[2] != "Uncertainty") || (fields[3] != "CPS") )
      {
        throw runtime_error( "Second line of PeakEasy CSV file is not correct: '" + line + "'" );
      }
      
      csv_format = PeakCsvFormat::PeakEasy;
      break;
    }else
    {
      throw runtime_error( "Invalid peak CSV line: '" + line + "'" );
    }
  }//while( std::getline( csv, line ) )
  
  
  size_t energy_index, amplitude_index, amplitude_sigma_index, nuclide_index, nuclide_energy_index, fwhm_kev_index;
  
  switch( csv_format )
  {
    case PeakCsvFormat::PeakEasy:
      energy_index = 0;
      amplitude_index= 1;
      amplitude_sigma_index = 2;
      nuclide_index = 13;
      nuclide_energy_index = 14;
      fwhm_kev_index = 4;
      //const size_t peak_cps_index = 3, fwhm_percent_index = 5;
      //const size_t roi_total_counts_index = 6, roi_id_index = 7, filename_index = 8;
      //const size_t live_time_index = 9, date_index = 10, time_index = 11;
      break;
      
    case PeakCsvFormat::Gadras:
      energy_index = 0;
      amplitude_index = 2;
      amplitude_sigma_index = 3;
      nuclide_index = nuclide_energy_index = 0;
      
      //const size_t sigma_index = 1; //uncert in energy
      fwhm_kev_index = 4;//, fwhm_sigma_index = 5, leakage_per_second_index = 6;
      //const size_t centroid_index = 7, filename_index = 8, record_idx_index = 9;
      //const size_t title_index = 10, date_time_index = 11;
      break;
      
    case PeakCsvFormat::Unknown:
      throw runtime_error( "Not a peak CSV file." );
  }//switch( csv_format )


  // Validate the column indices we'll use are reachable in the minimum
  //  field count we'll later require for each row (at least 9 fields).
  //  Indices that are 0 (unused / not-applicable) are allowed.
  const size_t required_fields = std::max( static_cast<size_t>(9),
        std::max({ energy_index, amplitude_index, amplitude_sigma_index, fwhm_kev_index })
        + 1 );

  vector<RelActCalcManual::GenericPeakInfo> answer;
  
  while( std::getline(csv, line) )
  {
    SpecUtils::trim(line);
    if( line.empty() || line[0]=='#' || (!isdigit(line[0]) && line[0]!='+' && line[0]!='-') )
      continue;
    
    vector<string> fields;
    split_and_trim( fields, line );
    
    const size_t nfields = fields.size();
    if( nfields == 0 )
      continue;

    if( nfields < required_fields )
      throw runtime_error( "Encountered line in GADRAS CSV file with only "
                          + std::to_string(nfields) + " fields.\n\tLine: \"" + line + "\"" );
    
    try
    {
      auto parse_field = []( const std::string &field, const char *what ) -> double {
        double val = 0.0;
        if( !SpecUtils::parse_double( field.c_str(), field.size(), val ) )
          throw runtime_error( std::string(what) + ": expected a number, but got '" + field + "'" );
        return val;
      };

      RelActCalcManual::GenericPeakInfo info;
      info.m_mean = info.m_energy = parse_field( fields[energy_index], "energy" );
      info.m_fwhm = parse_field( fields[fwhm_kev_index], "FWHM" );
      info.m_counts = parse_field( fields[amplitude_index], "amplitude" );
      info.m_counts_uncert = parse_field( fields[amplitude_sigma_index], "amplitude uncertainty" );
      
      //cout << "Found peak at " << info.m_energy << " keV"
      //     << ", with Amp=" << info.m_counts << " +- " << info.m_counts_uncert
      //     << endl;
      
      if( nuclide_index
         && nuclide_energy_index
         && (fields.size() > nuclide_energy_index)
         && !fields[nuclide_index].empty()
         && !fields[nuclide_energy_index].empty() )
      {
        contained_nuc_ids = true;
        
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        
        const SandiaDecay::Nuclide *nuc = db->nuclide( fields[nuclide_index] );
        if( !nuc )
          throw runtime_error( "Invalid Nuclide ID: " + fields[nuclide_index] );
        
        double nuc_energy = 0.0;
        if( !SpecUtils::parse_double( fields[nuclide_energy_index].c_str(),
                                      fields[nuclide_energy_index].size(), nuc_energy ) )
          throw runtime_error( "nuclide energy: expected a number, but got '"
                              + fields[nuclide_energy_index] + "'" );
        info.m_energy = nuc_energy;
        
        const double age = PeakDef::defaultDecayTime( nuc, nullptr );
        
        const double fwhm = info.m_fwhm;
        
        //size_t transition_index = 0;
        //const SandiaDecay::Transition *transition = nullptr;
        //PeakDef::SourceGammaType sourceGammaType;
        //PeakDef::findNearestPhotopeak( nuc, nuc_energy, -1.0, false, -1.0,
        //                               transition, transition_index, sourceGammaType );
        
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nuc, 1.0, age);
        const vector<SandiaDecay::EnergyRatePair> gammas = mix.gammas(0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, false);
        double yield = 0.0;
        for( const auto &erp : gammas )
        {
          // We'll sum all gammas within 1 FWHM of the specified gamma energy; not perfect, but good enough for now.
          if( fabs(erp.energy - nuc_energy) < fwhm )
            yield += erp.numPerSecond;
        }
        
        info.m_source_gammas.push_back( {yield, nuc->symbol} );
      }//if( we have nuclide ID from file )
      
      answer.push_back( info );
    }catch( std::exception &e )
    {
      throw runtime_error( "Invalid value on line '" + line + "', " + string(e.what()) );
    }//try / catch to parse a line into a peak
  }//while( SpecUtils::safe_get_line(csv, line, 2048) )
  
  if( answer.empty() )
    throw runtime_error( "No peak rows found in file." );
  
  return answer;
}//vector<GenericPeakInfo> peak_csv_to_peaks( istream &csv )

}//namespace PeakCsvInput

}//namespace RelActCalcManual
