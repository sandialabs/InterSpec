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

#include <cmath>
#include <deque>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "Eigen/Dense"

#include "ceres/ceres.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;


namespace
{
  /** The value part of a double or a ceres::Jet<>. */
  template<typename T>
  double scalar_of( const T &v )
  {
    if constexpr( std::is_arithmetic_v<T> )
      return static_cast<double>( v );
    else
      return static_cast<double>( v.a );
  }


  /** Solver options shared by the FWHM and efficiency fits. */
  ceres::Solver::Options make_solver_options()
  {
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.use_nonmonotonic_steps = true;
    options.max_num_iterations = 500;
    options.function_tolerance = 1.0e-10;
    options.gradient_tolerance = 1.0e-12;
    options.parameter_tolerance = 1.0e-12;
    options.num_threads = 1;
    options.logging_type = ceres::SILENT;
    options.minimizer_progress_to_stdout = false;
    return options;
  }//make_solver_options()


  /** Computes the (ambient-space) covariance of the single parameter block of `problem`, using
   the SVD based estimator.  Returns an empty vector if it cannot be computed.

   `ceres::Covariance` evaluates the Jacobian and runs Eigen's SVD on it, which is not safe on
   a non-finite Jacobian, so the residuals and Jacobian are checked first.
   */
  vector<double> problem_covariance( ceres::Problem &problem, double *pars, const size_t num_pars )
  {
    vector<double> residuals;
    ceres::CRSMatrix jacobian;
    double cost = 0.0;
    if( !problem.Evaluate( ceres::Problem::EvaluateOptions(), &cost, &residuals, nullptr, &jacobian ) )
      return {};

    for( const double r : residuals )
      if( std::isnan(r) || std::isinf(r) )
        return {};
    for( const double j : jacobian.values )
      if( std::isnan(j) || std::isinf(j) )
        return {};

    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::DENSE_SVD;
    cov_options.null_space_rank = -1;
    cov_options.min_reciprocal_condition_number = 1.0e-14;

    ceres::Covariance covariance( cov_options );
    vector<pair<const double *, const double *>> cov_blocks;
    cov_blocks.emplace_back( pars, pars );

    vector<double> cov( num_pars * num_pars, 0.0 );
    if( !covariance.Compute( cov_blocks, &problem )
        || !covariance.GetCovarianceBlock( pars, pars, cov.data() ) )
      return {};

    for( const double c : cov )
      if( std::isnan(c) || std::isinf(c) )
        return {};

    return cov;
  }//problem_covariance(...)


  /** The 1-sigma width uncertainty every FWHM fit uses - see MakeDrfFit::peak_width_chi2. */
  double sigma_uncert_for_fit( const PeakDef &peak )
  {
    const double sigma = peak.sigma();
    return (peak.sigmaUncert() > 0.0) ? std::max( peak.sigmaUncert(), 0.01*sigma ) : 0.05*sigma;
  }


  /** Residuals (sigma_pred - sigma_meas)/sigma_uncert for each Gaussian peak. */
  struct FwhmCostFunctor
  {
    vector<double> m_energies, m_sigmas, m_sigma_uncerts;
    DetectorPeakResponse::ResolutionFnctForm m_form;
    size_t m_num_pars;

    FwhmCostFunctor( const deque<shared_ptr<const PeakDef>> &peaks,
                     const DetectorPeakResponse::ResolutionFnctForm form,
                     const size_t num_pars )
      : m_form( form ), m_num_pars( num_pars )
    {
      for( const shared_ptr<const PeakDef> &peak : peaks )
      {
        if( !peak || !peak->gausPeak() )
          continue;
        m_energies.push_back( peak->mean() );
        m_sigmas.push_back( peak->sigma() );
        m_sigma_uncerts.push_back( sigma_uncert_for_fit(*peak) );
      }
    }//FwhmCostFunctor(...)

    size_t num_residuals() const { return m_energies.size(); }

    template<typename T>
    bool operator()( T const *const *parameters, T *residuals ) const
    {
      const T * const pars = parameters[0];
      try
      {
        for( size_t i = 0; i < m_energies.size(); ++i )
        {
          const T fwhm = DetectorPeakResponse::peakResolutionFWHM( T(m_energies[i]), m_form,
                                                                   pars, m_num_pars );
          const double val = scalar_of( fwhm );
          if( std::isnan(val) || std::isinf(val) )
            return false;  //Ceres treats this as an invalid trial and shrinks the step
          residuals[i] = (fwhm / 2.35482 - m_sigmas[i]) / m_sigma_uncerts[i];
        }
      }catch( std::exception & )
      {
        return false;
      }
      return true;
    }//operator()

    /** chi2 at `pars`; +infinity if the form cant be evaluated there. */
    double chi2( const vector<double> &pars ) const
    {
      vector<double> residuals( num_residuals(), 0.0 );
      const double * const blocks[1] = { pars.data() };
      if( !(*this)( blocks, residuals.data() ) )
        return std::numeric_limits<double>::infinity();
      double chi2 = 0.0;
      for( const double r : residuals )
        chi2 += r*r;
      return chi2;
    }//chi2(...)
  };//struct FwhmCostFunctor


  /** Parameter bounds, NaN meaning unbounded. */
  struct ParBounds
  {
    double lower = std::numeric_limits<double>::quiet_NaN();
    double upper = std::numeric_limits<double>::quiet_NaN();
  };


  /** Builds the Ceres problem for the FWHM fit around `pars` (which must outlive the problem). */
  void setup_fwhm_problem( ceres::Problem &problem, const FwhmCostFunctor &functor,
                           vector<double> &pars, const vector<int> &constant_pars,
                           const vector<ParBounds> &bounds )
  {
    const size_t num_pars = pars.size();
    // Ceres wants a non-const functor pointer, though it only ever calls the const operator()
    auto cost_function = new ceres::DynamicAutoDiffCostFunction<FwhmCostFunctor,4>(
                              const_cast<FwhmCostFunctor *>(&functor),
                              ceres::Ownership::DO_NOT_TAKE_OWNERSHIP );
    cost_function->AddParameterBlock( static_cast<int>(num_pars) );
    cost_function->SetNumResiduals( static_cast<int>(functor.num_residuals()) );
    problem.AddResidualBlock( cost_function, nullptr, pars.data() );  //problem owns cost_function

    if( !constant_pars.empty() )
    {
      if( constant_pars.size() == num_pars )
      {
        problem.SetParameterBlockConstant( pars.data() );
      }else
      {
        ceres::SubsetManifold * const manifold
                        = new ceres::SubsetManifold( static_cast<int>(num_pars), constant_pars );
        problem.SetManifold( pars.data(), manifold );  //problem owns manifold
      }
    }//if( !constant_pars.empty() )

    for( size_t i = 0; i < num_pars; ++i )
    {
      const bool is_const = (std::find( begin(constant_pars), end(constant_pars),
                                        static_cast<int>(i) ) != end(constant_pars));
      if( is_const || (i >= bounds.size()) )
        continue;
      if( !std::isnan(bounds[i].lower) )
      {
        pars[i] = std::max( pars[i], bounds[i].lower );
        problem.SetParameterLowerBound( pars.data(), static_cast<int>(i), bounds[i].lower );
      }
      if( !std::isnan(bounds[i].upper) )
      {
        pars[i] = std::min( pars[i], bounds[i].upper );
        problem.SetParameterUpperBound( pars.data(), static_cast<int>(i), bounds[i].upper );
      }
    }//for( size_t i = 0; i < num_pars; ++i )
  }//setup_fwhm_problem(...)


  struct FwhmSolution
  {
    vector<double> pars;
    double chi2 = std::numeric_limits<double>::infinity();
    bool converged = false;
  };

  /** One LM solve from `start`; never throws (a failed solve keeps chi2 = infinity). */
  FwhmSolution solve_fwhm( const FwhmCostFunctor &functor, vector<double> start,
                           const vector<int> &constant_pars, const vector<ParBounds> &bounds )
  {
    FwhmSolution answer;
    ceres::Problem problem;
    setup_fwhm_problem( problem, functor, start, constant_pars, bounds );

    ceres::Solver::Summary summary;
    ceres::Solve( make_solver_options(), &problem, &summary );

    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        answer.converged = true;
        break;
      case ceres::NO_CONVERGENCE:
        break;
      case ceres::FAILURE:
      case ceres::USER_FAILURE:
        return answer;
    }//switch( summary.termination_type )

    answer.pars = start;
    answer.chi2 = functor.chi2( start );
    return answer;
  }//solve_fwhm(...)


  /** Residuals L^{-1}*(model/measured - 1), with L the Cholesky factor of the fractional data
   covariance, so correlated (per-source) errors are weighted correctly. */
  struct EffCostFunctor
  {
    vector<double> m_energies, m_meas;
    Eigen::MatrixXd m_whiten;   //L^{-1}; lower triangular
    vector<double> m_par_scales;

    size_t num_residuals() const { return m_meas.size(); }
    size_t num_parameters() const { return m_par_scales.size(); }

    template<typename T>
    bool operator()( T const *const *parameters, T *residuals ) const
    {
      const size_t n = m_meas.size(), ncoef = m_par_scales.size();
      const T * const pars = parameters[0];

      vector<T> coefs( ncoef );
      for( size_t k = 0; k < ncoef; ++k )
        coefs[k] = pars[k] * m_par_scales[k];

      vector<T> raw( n );
      for( size_t i = 0; i < n; ++i )
      {
        const T model = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( T(m_energies[i]),
                                                                            coefs.data(), ncoef );
        const double val = scalar_of( model );
        if( std::isnan(val) || std::isinf(val) )
          return false;
        raw[i] = model / m_meas[i] - 1.0;
      }

      for( size_t i = 0; i < n; ++i )
      {
        T sum( 0.0 );
        for( size_t j = 0; j <= i; ++j )
          sum += m_whiten(i,j) * raw[j];
        residuals[i] = sum;
      }

      return true;
    }//operator()

    /** chi2 at coefficient values `coefs` (un-scaled); +infinity if not evaluable. */
    double chi2_of_coefs( const vector<double> &coefs ) const
    {
      vector<double> pars( coefs.size() );
      for( size_t k = 0; k < coefs.size(); ++k )
        pars[k] = coefs[k] / m_par_scales[k];
      vector<double> residuals( num_residuals(), 0.0 );
      const double * const blocks[1] = { pars.data() };
      if( !(*this)( blocks, residuals.data() ) )
        return std::numeric_limits<double>::infinity();
      double chi2 = 0.0;
      for( const double r : residuals )
        chi2 += r*r;
      return chi2;
    }//chi2_of_coefs(...)
  };//struct EffCostFunctor


  /** L^{-1} for C = L*L^T.  If C is not numerically positive definite a small diagonal jitter is
   added (up to a few times); as a last resort the correlations are dropped. */
  Eigen::MatrixXd whitening_matrix( const vector<double> &cov_row_major, const size_t n,
                                    std::string &warnings )
  {
    Eigen::MatrixXd C( n, n );
    double trace = 0.0;
    for( size_t i = 0; i < n; ++i )
    {
      for( size_t j = 0; j < n; ++j )
        C(i,j) = cov_row_major[i*n + j];
      trace += C(i,i);
    }

    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity( n, n );
    double jitter = 1.0E-12 * trace / static_cast<double>( n );
    for( int attempt = 0; attempt < 4; ++attempt )
    {
      Eigen::LLT<Eigen::MatrixXd> llt( C );
      if( llt.info() == Eigen::Success )
      {
        const Eigen::MatrixXd L = llt.matrixL();
        return L.triangularView<Eigen::Lower>().solve( identity );
      }
      C += jitter * identity;
      jitter *= 100.0;
    }//for( int attempt = 0; attempt < 4; ++attempt )

    warnings += "Data covariance was not positive definite; correlations between points ignored. ";
    Eigen::MatrixXd answer = Eigen::MatrixXd::Zero( n, n );
    for( size_t i = 0; i < n; ++i )
      answer(i,i) = 1.0 / std::sqrt( std::max( cov_row_major[i*n + i], 1.0E-12 ) );
    return answer;
  }//whitening_matrix(...)
}//namespace


namespace MakeDrfFit
{
  
double peak_width_chi2( double predicted_sigma, const PeakDef &peak )
{
  if( !peak.gausPeak() )
    return 0.0;  //shouldnt ever get here
  
  const double measured_sigma = peak.sigma();
  //double measured_sigma_uncert = peak.sigmaUncert();
  //if( measured_sigma_uncert <= 0.01 )
  //  measured_sigma_uncert = 1.0;
  const double measured_sigma_uncert = ((peak.sigmaUncert() > 0.0) ? std::max( peak.sigmaUncert(), 0.01*measured_sigma) : 0.05*measured_sigma);
  
  const double chi = (measured_sigma - predicted_sigma)/measured_sigma_uncert;
  return chi*chi;
}//peak_width_chi2(...)



FwhmFitResult performResolutionFitEx( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                                      const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                                      const int sqrtEqnOrder,
                                      const std::vector<float> &starting_coefs )
{
  if( !peaks || peaks->empty() )
    throw runtime_error( "MakeDrfFit::performResolutionFit(...): no input peaks" );
  
  const vector<shared_ptr<const PeakDef>> peakv( begin(*peaks), end(*peaks) );
  const PeakFitUtils::CoarseResolutionType coarse_type
                            = PeakFitUtils::coarse_resolution_from_peaks( peakv );
  const bool highres = (coarse_type == PeakFitUtils::CoarseResolutionType::High);
  const size_t npeaks = peaks->size();

  // Default starting values and bounds; the bounds span both the high- and low-resolution
  //  detector ranges, because we could be wrong about the detector being high-resolution.
  double a_initial = 0.0, b_initial = 0.0, c_initial = 0.0;
  double lowerA = 0.0, upperA = 0.0, lowerB = 0.0, upperB = 0.0, lowerC = 0.0, upperC = 0.0;
  
  // Closed-form linear least squares seed, where the form allows it
  bool fit_using_lls = false;
  vector<float> lls_coefs, lls_uncerts;
  
  switch( fnctnlForm )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
    {
      lowerA = std::min( 0.75*1.0, 1.5*-7.0);
      upperA = std::max( 2.0*1.77, 1.5*7.44000);
      lowerB = std::min( 0.75*0.20028, 0.75*2.13000 );
      upperB = std::max( 2.0*0.27759, 1.5*8.50000 );
      lowerC = std::min( 0.75*0.31, 0.75*0.20000 );
      upperC = std::max( 1.5*0.57223, 1.25*0.70000 );
      
      if( starting_coefs.size() == 3 )
      { //caller has provided some default values, lets use them
        a_initial = static_cast<double>( starting_coefs[0] );
        b_initial = static_cast<double>( starting_coefs[1] );
        c_initial = static_cast<double>( starting_coefs[2] );
      }else if( highres )
      {
        a_initial = 0.5*(1.0 + 1.77);
        b_initial = 0.5*(0.20028 + 0.27759);
        c_initial = 0.5*(0.31 + 0.57223);
      }else
      {
        a_initial = 0.0;
        b_initial = 0.5*(2.13000 + 8.50000);
        c_initial = 0.5*(0.20000 + 0.70000);
      }
      
      break;
    }//case kGadrasResolutionFcn:
      
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
    {
      a_initial = highres ? 2.6 : 100;
      lowerA = std::min( -2.5, -100.0 );
      upperA = std::max( 7.5, 400.0 );
      
      b_initial = highres ? 1.0 : (3600.0 / 0.661);
      lowerB = std::min( -5.0, 0.0);
      upperB = std::max( 5.0, (180*180 / 0.661) );
      
      c_initial = 0.0;
      lowerC = std::min( -5.0, -10000.0);
      upperC = std::max( 5.0, 10000.0);
      
      try
      {
        fit_sqrt_poly_fwhm_lls( *peaks, 3, true, lls_coefs, lls_uncerts );
        fit_using_lls = (lls_coefs.size() == 3);
      }catch( std::exception & )
      {
      }
      
      break;
    }//case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      
    case DetectorPeakResponse::kConstantPlusSqrtEnergy:
    {
      //Based on pretty much nothing
      a_initial = highres ? 1 : 0;
      lowerA = std::min( -10.0, -25.0 );
      upperA = std::max( 10.0, 50.0 );
      
      b_initial = highres ? 0.035 : 2.0;
      lowerB = std::min( 0.0, 0.0 );
      upperB = std::max( 5.0, 50.0 );
      
      try
      {
        fit_constant_plus_sqrt_fwhm_lls( *peaks, lls_coefs, lls_uncerts );
        fit_using_lls = (lls_coefs.size() == 2);
      }catch( std::exception & )
      {
      }
      
      break;
    }//case DetectorPeakResponse::kConstantPlusSqrtEnergy:
      
    case DetectorPeakResponse::kSqrtPolynomial:
    {
      if( sqrtEqnOrder < 1 )
        throw runtime_error( "performResolutionFit: sqrt eqn order should be at least 1" );
       
      //Based on pretty much nothing
      a_initial = highres ? 2.6 : 100.0;
      lowerA = std::min( -2.5, -100.0 );
      upperA = std::max( 7.5, 400.0 );
      
      b_initial = highres ? 1.0 : (3600.0 / 0.661);
      lowerB = std::min( -5.0, 0.0 );
      upperB = std::max( 5.0, (180*180 / 0.661) );
      
      c_initial = 0;
      lowerC = std::min( -5.0, (-3600.0/(0.661*0.661)) );
      upperC = std::max( 5.0, (3600.0/(0.661*0.661)) );
      
      try
      {
        fit_sqrt_poly_fwhm_lls( *peaks, sqrtEqnOrder, false, lls_coefs, lls_uncerts );
        fit_using_lls = (static_cast<int>(lls_coefs.size()) == sqrtEqnOrder);
      }catch( std::exception & )
      {
      }
      
      break;
    }//case kSqrtPolynomial:
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      throw runtime_error( "MakeDrfFit::performResolutionFit(...):"
                          " invalid ResolutionFnctForm" );
  }//switch( fnctnlForm )
  
  // Starting parameters, which are held fixed, and bounds - depending on how many peaks we have
  vector<double> start;
  vector<int> constant_pars;
  vector<ParBounds> bounds;
  auto bounded = []( const double lower, const double upper ){
    ParBounds b;
    b.lower = lower;
    b.upper = upper;
    return b;
  };
  
  bool search_gadras_a = false;
  
  switch( fnctnlForm )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
    {
      start = { a_initial, b_initial, c_initial };
      bounds = { bounded(lowerA,upperA), bounded(lowerB,upperB), bounded(lowerC,upperC) };
      if( npeaks == 1 )
        constant_pars = { 0, 2 };
      else if( npeaks == 2 )
        constant_pars = { 0 };
      else
        search_gadras_a = true;
      break;
    }//case kGadrasResolutionFcn:
      
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
    {
      if( fit_using_lls )
      {
        start.assign( begin(lls_coefs), end(lls_coefs) );
      }else
      {
        bounds = { bounded(lowerA,upperA), bounded(lowerB,upperB), bounded(lowerC,upperC) };
        if( npeaks == 1 )
        {
          start = { 0.0, b_initial, 0.0 };
          constant_pars = { 0, 2 };
        }else if( npeaks == 2 )
        {
          start = { a_initial, b_initial, 0.0 };
          constant_pars = { 2 };
        }else
        {
          start = { a_initial, b_initial, c_initial };
        }
      }//if( fit_using_lls ) / else
      break;
    }//case kSqrtEnergyPlusInverse:
      
    case DetectorPeakResponse::kConstantPlusSqrtEnergy:
    {
      if( fit_using_lls )
      {
        start.assign( begin(lls_coefs), end(lls_coefs) );
      }else
      {
        bounds = { bounded(lowerA,upperA), bounded(lowerB,upperB) };
        if( npeaks == 1 )
        {
          start = { 0.0, b_initial };
          constant_pars = { 0 };
        }else
        {
          start = { a_initial, b_initial };
        }
      }//if( fit_using_lls ) / else
      break;
    }//case kConstantPlusSqrtEnergy:
      
    case DetectorPeakResponse::kSqrtPolynomial:
    {
      if( fit_using_lls )
      {
        start.assign( begin(lls_coefs), end(lls_coefs) );
      }else
      {
        // Note: may be fewer coefficients than `sqrtEqnOrder` (callers rely on this)
        if( npeaks <= 3 )
        {
          start = { a_initial, b_initial };
          bounds = { bounded(lowerA,upperA), bounded(lowerB,upperB) };
          if( npeaks < 3 )
            constant_pars = { 1 };
        }else
        {
          start = { a_initial, b_initial, c_initial };
          bounds = { bounded(lowerA,upperA), bounded(lowerB,upperB), bounded(lowerC,upperC) };
        }
      }//if( fit_using_lls ) / else
      break;
    }//case kSqrtPolynomial:
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      assert( 0 );
      break;
  }//switch( fnctnlForm )
  
  const size_t num_pars = start.size();
  const FwhmCostFunctor functor( *peaks, fnctnlForm, num_pars );
  if( !functor.num_residuals() )
    throw runtime_error( "MakeDrfFit::performResolutionFit(...): no Gaussian peaks" );
  
  FwhmSolution best;
  
  if( search_gadras_a )
  {
    // The sign of A selects a different functional form below 661 keV, so a gradient method
    //  cant cross zero on its own: solve each branch from several starts and keep the best.
    const double a_neg_max = -1.0E-6;
    for( const bool negative : { true, false } )
    {
      vector<ParBounds> region_bounds = bounds;
      region_bounds[0] = negative ? bounded(lowerA, a_neg_max) : bounded(0.0, upperA);
      
      vector<double> a_starts = negative ? vector<double>{ -8.0, -4.0, -1.5, -0.3 }
                                         : vector<double>{ 0.0, 0.5, 1.4, 4.0 };
      if( (negative && (a_initial < 0.0)) || (!negative && (a_initial >= 0.0)) )
        a_starts.push_back( a_initial );
      
      for( const double a0 : a_starts )
      {
        const FwhmSolution trial = solve_fwhm( functor, { a0, b_initial, c_initial },
                                               constant_pars, region_bounds );
        if( trial.chi2 < best.chi2 )
          best = trial;
      }
    }//for( const bool negative : { true, false } )
  }else
  {
    best = solve_fwhm( functor, start, constant_pars, bounds );
  }//if( search_gadras_a ) / else
  
  if( fit_using_lls )
  {
    // Keep the closed-form answer if the nonlinear refinement didnt improve on it
    const vector<double> lls_pars( begin(lls_coefs), end(lls_coefs) );
    const double lls_chi2 = functor.chi2( lls_pars );
    if( lls_chi2 <= best.chi2 )
    {
      best.pars = lls_pars;
      best.chi2 = lls_chi2;
      best.converged = true;
    }
  }//if( fit_using_lls )
  
  if( best.pars.empty() || std::isinf(best.chi2) || std::isnan(best.chi2) )
    throw runtime_error( "FWHM response function fit failed to find a valid solution" );
  
  FwhmFitResult result;
  result.chi2 = best.chi2;
  result.dof = static_cast<int>( functor.num_residuals() ) - static_cast<int>( num_pars - constant_pars.size() );
  if( !best.converged )
    result.warnings += "The FWHM fit did not fully converge. ";
  
  result.coefs.resize( num_pars );
  for( size_t i = 0; i < num_pars; ++i )
    result.coefs[i] = static_cast<float>( best.pars[i] );
  
  // Covariance at the solution (no further minimization; the problem is only built to evaluate it)
  {
    vector<double> pars = best.pars;
    ceres::Problem problem;
    // Bounds only matter for a solve; leave them out so a solution sitting on a bound isnt moved.
    setup_fwhm_problem( problem, functor, pars, constant_pars, {} );
    const vector<double> cov = problem_covariance( problem, pars.data(), num_pars );
    
    result.uncerts.assign( num_pars, 0.0f );
    if( cov.size() == num_pars*num_pars )
    {
      result.covRowMajor.assign( begin(cov), end(cov) );
      for( size_t i = 0; i < num_pars; ++i )
        result.uncerts[i] = static_cast<float>( std::sqrt( std::max( 0.0, cov[i*num_pars + i] ) ) );
    }else
    {
      result.warnings += "FWHM coefficient uncertainties could not be computed. ";
    }
  }
  
  return result;
}//performResolutionFitEx(...)


double performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                           const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                           const int sqrtEqnOrder,
                           std::vector<float> &answer,
                           std::vector<float> &uncerts )
{
  const FwhmFitResult result = performResolutionFitEx( peaks, fnctnlForm, sqrtEqnOrder, answer );
  answer = result.coefs;
  uncerts = result.uncerts;
  return result.chi2;
}//performResolutionFit(...)

  
  //removeOutlyingWidthPeaks(...): removes peaks whos width does not agree
  //  well with the functional form passed in.  Returns surviving peaks.
std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>>
removeOutlyingWidthPeaks( const std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &peaks,
                          const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                          const std::vector<float> &coefficients )
{
    const double npeaks = static_cast<double>( peaks->size() );
    if( npeaks < 5 )
      return peaks;
    const size_t ndel_max = static_cast<size_t>( floor(0.2*npeaks) );
    
    vector<double> weights;
    double mean_weight = 0.0;
    for( const PeakModel::PeakShrdPtr peak : *peaks )
    {
      double predicted_sigma = DetectorPeakResponse::peakResolutionSigma( peak->mean(), fnctnlForm, coefficients );
      if( IsNan(predicted_sigma) || IsInf(predicted_sigma) )
        predicted_sigma = 0.0;
      
      const double chi2 = MakeDrfFit::peak_width_chi2( predicted_sigma, *peak );
      
      mean_weight += chi2/npeaks;
      weights.push_back( chi2 );
    }//for( const EnergySigma &es : m_energies_and_sigmas )
    
    vector<size_t> indices( npeaks );
    for( size_t i = 0; i < npeaks; ++i )
      indices[i] = i;
    
    std::sort( begin(indices), end(indices), [&weights](const size_t &a, const size_t &b){
      return weights[b] < weights[a];
    } );
    
    size_t lastInd = 0;
    while( lastInd <= ndel_max && weights[indices[lastInd]] > 2.5*mean_weight )
      ++lastInd;
    indices.erase( indices.begin() + lastInd, indices.end() );
    
    std::shared_ptr< deque< std::shared_ptr<const PeakDef> > > reduced_peaks( new deque< PeakModel::PeakShrdPtr >() );
    for( size_t i = 0; i < npeaks; ++i )
      if( find( indices.begin(), indices.end(), i) == indices.end() )
        reduced_peaks->push_back( peaks->operator[](i) );
    
    return reduced_peaks;
}//removeOutlyingWidthPeaks(...)
  
  
  
double fit_sqrt_poly_fwhm_lls( const std::deque< std::shared_ptr<const PeakDef> > &peaks,
                               const int num_fit_coefficients,
                               const bool include_inv_term,
                               std::vector<float> &coeffs,
                               std::vector<float> &coeff_uncerts )
{
  size_t nbin = 0;
  for( size_t i = 0; i < peaks.size(); ++i )
    nbin += (peaks[i] && peaks[i]->gausPeak() ? 1 : 0);
    
  if( num_fit_coefficients < 1 )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls: num_fit_coefficients must be >= 1" );
  
  if( !nbin )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls: must have at least 1 input peak" );
  
  if( nbin < static_cast<size_t>(num_fit_coefficients) )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls: must have at least as many peaks as coefficients to fit to" );
  
  //log(eff(x)) = A0 + A1*logx + A2*logx^2 + A3*logx^3, where x is energy in MeV
  vector<float> x( nbin, 0.0f ), widths( nbin, 0.0f ), widths_uncert( nbin, 0.0f );
  
  for( size_t i = 0, peak_num = 0; i < peaks.size(); ++i )
  {
    const shared_ptr<const PeakDef> &peak = peaks[i];
    if( peak && peak->gausPeak() )
    {
      x[peak_num] = peak->mean() / (include_inv_term ? 1.0f : 1000.0f);
      widths[peak_num] = peak->fwhm();
      widths_uncert[peak_num] = 2.35482*((peak->sigmaUncert() > 0.0) 
                                          ? std::max( peak->sigmaUncert(), 0.01*widths[peak_num])
                                          : 0.05*widths[peak_num]);
      peak_num += 1;
      
      assert( peak_num <= nbin );
    }//if( valid peak we can use )
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation using Eigen SVD for numerical stability
  
  Eigen::MatrixX<double> A( nbin, num_fit_coefficients );
  Eigen::VectorX<double> b( nbin );
  
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const double data_y = widths[row] * widths[row];
    const double data_y_uncert = widths_uncert[row] * widths_uncert[row];
    
    b(row) = data_y / data_y_uncert;
    
    if( include_inv_term )
    {
      // Nominally we will only ever call this function to fit for the FRAM style FWHM  when num_fit_coefficients == 3,
      //  but we'll go ahead and code in possibility to fit for different num_fit_coefficients.  For num_fit_coefficients==1, we'll
      //  fit the energy dependent term, for num_fit_coefficients==2, the constant + energy dependent, and num_fit_coefficients==3
      //  the inverse term as well.  Higher num_fit_coefficients terms we'll do as power series in energy.
      //  Before we return, we'll swap the zeroth and first coefficients so things are as expected.
      assert( num_fit_coefficients == 3 );
      for( int col = 0; col < num_fit_coefficients; ++col )
      {
        if( col == 0 )
          A(row,col) = x[row] / data_y_uncert;
        else if( col == 1 )
          A(row,col) = 1.0 / data_y_uncert;
        else if( col == 2 )
          A(row,col) = (1.0/x[row]) / data_y_uncert;
        else
          A(row,col) = std::pow( x[row], double(col-1)) / data_y_uncert;
      }
    }else
    {
      for( int col = 0; col < num_fit_coefficients; ++col )
        A(row,col) = std::pow( x[row], static_cast<double>(col)) / data_y_uncert;
    }
  }//for( int col = 0; col < num_fit_coefficients; ++col )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  
  const Eigen::MatrixX<double> A_transpose = A.transpose();
  const Eigen::MatrixX<double> alpha = A_transpose * A;
  const Eigen::MatrixX<double> C = alpha.inverse();
  
  coeffs.resize( num_fit_coefficients );
  coeff_uncerts.resize( num_fit_coefficients );
  for( int coef = 0; coef < num_fit_coefficients; ++coef )
  {
    coeffs[coef] = static_cast<float>( a(coef) );
    coeff_uncerts[coef] = static_cast<float>( std::sqrt( C(coef,coef) ) );
  }//for( int coef = 0; coef < num_fit_coefficients; ++coef )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    double y_pred = 0.0;
    if( include_inv_term )
    {
      for( int i = 0; i < num_fit_coefficients; ++i )
      {
        if( i == 0 )
          y_pred += a(i) * x[bin];
        else if( i == 1 )
          y_pred += a(i);
        else if( i == 2 )
          y_pred += a(i) / x[bin];
        else
          y_pred += a(i) * std::pow( x[bin], static_cast<double>(i-1) );
      }
    }else
    {
      for( int i = 0; i < num_fit_coefficients; ++i )
        y_pred += a(i) * std::pow( x[bin], static_cast<double>(i) );
    }//if( include_inv_term ) / else
    
    y_pred = sqrt( y_pred );
    chi2 += std::pow( (y_pred - widths[bin]) / widths_uncert[bin], 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  if( include_inv_term && (coeffs.size() > 1) )
  {
    // Swap the zeroth and first coefficients so things are as expected.
    //  (i.e., constant term is zeroth, and energy dependent term is at index==1)
    std::swap( coeffs[0], coeffs[1] );
    std::swap( coeff_uncerts[0], coeff_uncerts[1] );
  }
  
  return chi2;
}//double fit_sqrt_poly_fwhm_lls(...)
  
  
double fit_constant_plus_sqrt_fwhm_lls( const std::deque< std::shared_ptr<const PeakDef> > &peaks,
                                  std::vector<float> &coeffs,
                                  std::vector<float> &coeff_uncerts )
{
  const size_t nbin = peaks.size();
  
  if( nbin < 2 )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls: must have at least 2 input peak" );
  
  //log(eff(x)) = A0 + A1*logx + A2*logx^2 + A3*logx^3, where x is energy in MeV
  vector<float> x, widths, widths_uncert;
  x.resize( peaks.size() );
  widths.resize( peaks.size() );
  widths_uncert.resize( peaks.size() );
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    if( peaks[i]->gausPeak() )
    {
      x[i] = peaks[i]->mean();
      widths[i] = peaks[i]->fwhm();
      widths_uncert[i] = 2.35482*((peaks[i]->sigmaUncert() > 0.0) ? std::max( peaks[i]->sigmaUncert(), 0.01*widths[i]) : 0.05*widths[i]);
    }
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation using Eigen SVD for numerical stability
  
  const size_t num_fit_coefficients = 2;
  Eigen::MatrixX<double> A( nbin, num_fit_coefficients );
  Eigen::VectorX<double> b( nbin );
  
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const double data_y = widths[row];
    const double data_y_uncert = widths_uncert[row];
    
    b(row) = data_y / data_y_uncert;
    
    A(row,0) = 1.0 / data_y_uncert;
    A(row,1) = std::sqrt( x[row] ) / data_y_uncert;
  }//for( int col = 0; col < num_fit_coefficients; ++col )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  
  const Eigen::MatrixX<double> A_transpose = A.transpose();
  const Eigen::MatrixX<double> alpha = A_transpose * A;
  const Eigen::MatrixX<double> C = alpha.inverse();
  
  coeffs.resize( num_fit_coefficients );
  coeff_uncerts.resize( num_fit_coefficients );
  for( int coef = 0; coef < num_fit_coefficients; ++coef )
  {
    coeffs[coef] = static_cast<float>( a(coef) );
    coeff_uncerts[coef] = static_cast<float>( std::sqrt( C(coef,coef) ) );
  }//for( int coef = 0; coef < num_fit_coefficients; ++coef )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    const double y_pred = a(0) + a(1) * sqrt( x[bin] );
    chi2 += std::pow( (y_pred - widths[bin]) / widths_uncert[bin], 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_constant_plus_sqrt_fwhm_lls(...)
  
  
std::vector<double> effDataCovariance( const std::vector<EffFitPoint> &data )
{
  const size_t n = data.size();
  vector<double> cov( n*n, 0.0 );
  for( size_t i = 0; i < n; ++i )
  {
    const EffFitPoint &a = data[i];
    cov[i*n + i] += static_cast<double>(a.fracStatUncert) * a.fracStatUncert;
    if( a.sourceKey.empty() )
    {
      cov[i*n + i] += static_cast<double>(a.fracCertUncert) * a.fracCertUncert
                      + static_cast<double>(a.fracDistUncert) * a.fracDistUncert;
      continue;
    }
    
    for( size_t j = 0; j < n; ++j )
    {
      const EffFitPoint &b = data[j];
      if( a.sourceKey == b.sourceKey )
        cov[i*n + j] += static_cast<double>(a.fracCertUncert) * b.fracCertUncert
                        + static_cast<double>(a.fracDistUncert) * b.fracDistUncert;
    }
  }//for( size_t i = 0; i < n; ++i )
  
  return cov;
}//effDataCovariance(...)


double fit_intrinsic_eff_least_linear_squares( const std::vector<EffFitPoint> &data,
                                               const int order,
                                               std::vector<float> &coeffs,
                                               std::vector<float> &coeff_uncerts,
                                               std::vector<double> *cov_row_major )
{
  const size_t nbin = data.size();
  if( (order < 1) || (nbin < static_cast<size_t>(order)) )
    throw runtime_error( "fit_intrinsic_eff_least_linear_squares: invalid order" );
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //log(eff(x)) = A0 + A1*logx + A2*logx^2 + A3*logx^3 + ...
  //  The weight of each point is its total fractional uncertainty (= sigma of log(eff)).
  Eigen::MatrixX<double> A( nbin, order );
  Eigen::VectorX<double> b( nbin );
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const EffFitPoint &p = data[row];
    if( (p.energy <= 0.0f) || (p.efficiency <= 0.0f) )
      throw runtime_error( "fit_intrinsic_eff_least_linear_squares: non-positive energy or efficiency" );
    
    const double frac_uncert = std::max( 1.0E-6, std::sqrt( static_cast<double>(p.fracStatUncert)*p.fracStatUncert
                                                          + static_cast<double>(p.fracCertUncert)*p.fracCertUncert
                                                          + static_cast<double>(p.fracDistUncert)*p.fracDistUncert ) );
    const double logx = log( static_cast<double>(p.energy) );
    
    b(row) = log( static_cast<double>(p.efficiency) ) / frac_uncert;
    for( int col = 0; col < order; ++col )
      A(row,col) = std::pow( logx, double(col) ) / frac_uncert;
  }//for( size_t row = 0; row < nbin; ++row )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  const Eigen::MatrixX<double> C = (A.transpose() * A).inverse();
  
  coeffs.resize( order );
  coeff_uncerts.resize( order );
  for( int coef = 0; coef < order; ++coef )
  {
    coeffs[coef] = static_cast<float>( a(coef) );
    coeff_uncerts[coef] = static_cast<float>( std::sqrt( std::max( 0.0, C(coef,coef) ) ) );
  }
  
  if( cov_row_major )
  {
    cov_row_major->assign( order*order, 0.0 );
    for( int i = 0; i < order; ++i )
      for( int j = 0; j < order; ++j )
        (*cov_row_major)[i*order + j] = C(i,j);
  }//if( cov_row_major )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    const EffFitPoint &p = data[bin];
    const double frac_uncert = std::max( 1.0E-6, std::sqrt( static_cast<double>(p.fracStatUncert)*p.fracStatUncert
                                                          + static_cast<double>(p.fracCertUncert)*p.fracCertUncert
                                                          + static_cast<double>(p.fracDistUncert)*p.fracDistUncert ) );
    const double y_pred = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( static_cast<double>(p.energy),
                                                                              a.data(), a.size() );
    chi2 += std::pow( (y_pred - p.efficiency) / (frac_uncert * p.efficiency), 2.0 );
  }//for( size_t bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_intrinsic_eff_least_linear_squares(...)
  
  
EffFitResult performEfficiencyFit( const std::vector<EffFitPoint> &data, const int fcnOrder )
{
  if( data.empty() )
    throw runtime_error( "MakeDrfFit::performEfficiencyFit(...): no input peaks" );
  
  if( fcnOrder < 1 )
    throw runtime_error( "MakeDrfFit::performEfficiencyFit(...): requested fit order " + std::to_string(fcnOrder) + " (must be at least 1)" );
  
  if( fcnOrder > static_cast<int>(data.size()) )
    throw runtime_error( "MakeDrfFit::performEfficiencyFit(...): requested fit order " + std::to_string(fcnOrder)
                         + ", with only " + std::to_string(data.size()) + " data points." );
  
  for( const EffFitPoint &p : data )
  {
    if( (p.energy <= 0.0f) || (p.efficiency <= 0.0f) || std::isnan(p.efficiency) || std::isinf(p.efficiency)
        || (p.fracStatUncert < 0.0f) || (p.fracCertUncert < 0.0f) || (p.fracDistUncert < 0.0f) )
      throw runtime_error( "MakeDrfFit::performEfficiencyFit(...): invalid data point" );
  }
  
  const size_t n = data.size(), ncoef = static_cast<size_t>( fcnOrder );
  EffFitResult result;
  result.dof = static_cast<int>( n ) - fcnOrder;
  
  // Seed from the closed-form log-space fit (falls back to a flat curve if that fails)
  vector<float> seed_coefs, seed_uncerts;
  vector<double> seed_cov;
  bool have_seed = false;
  try
  {
    fit_intrinsic_eff_least_linear_squares( data, fcnOrder, seed_coefs, seed_uncerts, &seed_cov );
    have_seed = true;
    for( const float c : seed_coefs )
      have_seed = (have_seed && !std::isnan(c) && !std::isinf(c));
  }catch( std::exception & )
  {
  }
  
  if( !have_seed )
  {
    double mean_log_eff = 0.0;
    for( const EffFitPoint &p : data )
      mean_log_eff += log( static_cast<double>(p.efficiency) ) / static_cast<double>( n );
    seed_coefs.assign( ncoef, 0.0f );
    seed_coefs[0] = static_cast<float>( mean_log_eff );
    seed_cov.clear();
  }//if( !have_seed )
  
  // Problem setup
  EffCostFunctor functor;
  for( const EffFitPoint &p : data )
  {
    functor.m_energies.push_back( p.energy );
    functor.m_meas.push_back( p.efficiency );
  }
  functor.m_whiten = whitening_matrix( effDataCovariance(data), n, result.warnings );
  functor.m_par_scales.resize( ncoef );
  for( size_t k = 0; k < ncoef; ++k )
    functor.m_par_scales[k] = (fabs(seed_coefs[k]) > 1.0E-6) ? fabs(static_cast<double>(seed_coefs[k])) : 1.0;
  
  vector<double> pars( ncoef );
  for( size_t k = 0; k < ncoef; ++k )
    pars[k] = seed_coefs[k] / functor.m_par_scales[k];
  
  auto cost_function = new ceres::DynamicAutoDiffCostFunction<EffCostFunctor,4>(
                              &functor, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP );
  cost_function->AddParameterBlock( static_cast<int>(ncoef) );
  cost_function->SetNumResiduals( static_cast<int>(n) );
  
  ceres::Problem problem;
  problem.AddResidualBlock( cost_function, nullptr, pars.data() );  //problem owns cost_function
  
  ceres::Solver::Summary summary;
  ceres::Solve( make_solver_options(), &problem, &summary );
  
  switch( summary.termination_type )
  {
    case ceres::CONVERGENCE:
    case ceres::USER_SUCCESS:
      break;
      
    case ceres::NO_CONVERGENCE:
      result.warnings += "The efficiency fit did not fully converge. ";
      break;
      
    case ceres::FAILURE:
    case ceres::USER_FAILURE:
      if( !have_seed )
        throw runtime_error( "Efficiency function fit failed: " + summary.message );
      result.warnings += "The efficiency fit failed (" + summary.message + "); using the linear seed. ";
      for( size_t k = 0; k < ncoef; ++k )
        pars[k] = seed_coefs[k] / functor.m_par_scales[k];
      break;
  }//switch( summary.termination_type )
  
  vector<double> coefs( ncoef );
  for( size_t k = 0; k < ncoef; ++k )
    coefs[k] = pars[k] * functor.m_par_scales[k];
  
  double chi2 = functor.chi2_of_coefs( coefs );
  
  // Keep the closed-form answer if the nonlinear refinement didnt improve on it
  vector<double> cov;
  if( have_seed )
  {
    const vector<double> seed_dbl( begin(seed_coefs), end(seed_coefs) );
    const double seed_chi2 = functor.chi2_of_coefs( seed_dbl );
    if( seed_chi2 < chi2 )
    {
      coefs = seed_dbl;
      chi2 = seed_chi2;
      for( size_t k = 0; k < ncoef; ++k )
        pars[k] = coefs[k] / functor.m_par_scales[k];
    }
  }//if( have_seed )
  
  if( std::isnan(chi2) || std::isinf(chi2) )
    throw runtime_error( "Efficiency function fit did not find a valid solution" );
  
  // Coefficient covariance at the solution.  The residuals are whitened, so the SVD based
  //  estimate is the covariance of the (scaled) parameters directly.
  {
    const vector<double> scaled_cov = problem_covariance( problem, pars.data(), ncoef );
    if( scaled_cov.size() == ncoef*ncoef )
    {
      cov.resize( ncoef*ncoef );
      for( size_t i = 0; i < ncoef; ++i )
        for( size_t j = 0; j < ncoef; ++j )
          cov[i*ncoef + j] = functor.m_par_scales[i] * functor.m_par_scales[j] * scaled_cov[i*ncoef + j];
    }else if( seed_cov.size() == ncoef*ncoef )
    {
      cov = seed_cov;
      result.warnings += "Coefficient covariance taken from the linear seed. ";
    }else
    {
      result.warnings += "Coefficient uncertainties could not be computed. ";
    }
  }
  
  // Birge (variance of unit weight) inflation: when the points scatter more than their stated
  //  uncertainties allow, the reported curve uncertainty is scaled up to match the scatter.
  //  chi2/dof rather than a robust estimator, since there are only ever a handful of points; see
  //  the "variance of unit weight" discussion in RelActCalcManual.cpp for the trade-offs.
  result.chi2 = chi2;
  result.birgeScale = (result.dof > 0) ? std::max( 1.0, chi2 / result.dof ) : 1.0;
  
  result.coefs.resize( ncoef );
  for( size_t k = 0; k < ncoef; ++k )
    result.coefs[k] = static_cast<float>( coefs[k] );
  
  result.uncerts.assign( ncoef, 0.0f );
  if( cov.size() == ncoef*ncoef )
  {
    result.covRowMajor.resize( ncoef*ncoef );
    for( size_t i = 0; i < ncoef*ncoef; ++i )
      result.covRowMajor[i] = static_cast<float>( result.birgeScale * cov[i] );
    for( size_t k = 0; k < ncoef; ++k )
      result.uncerts[k] = static_cast<float>( std::sqrt( std::max( 0.0, result.birgeScale * cov[k*ncoef + k] ) ) );

    // Decide here, where the one warning channel is, whether the covariance can actually be stored -
    //  DetectorEfficiencyUncert refuses one that describes no possible set of errors, and downstream
    //  (MakeDrfCalc::assembleDrf) that refusal would drop it with nothing said to the user.
    const vector<double> cov_dbl( begin(result.covRowMajor), end(result.covRowMajor) );
    string why;
    if( !DetectorEfficiencyUncert::covarianceIsUsable( cov_dbl, &why ) )
    {
      result.covRowMajor.clear();
      result.warnings += "The efficiency fit's coefficient covariance could not be used (" + why
                         + "), so the equation will carry only the uncertainty the measured points"
                           " imply. ";
    }
  }//if( cov.size() == ncoef*ncoef )
  
  return result;
}//performEfficiencyFit(...)


double performEfficiencyFit( const std::vector<DetEffDataPoint> data,
                            const int fcnOrder,
                            std::vector<float> &result,
                            std::vector<float> &uncerts )
{
  vector<EffFitPoint> points;
  for( const DetEffDataPoint &d : data )
  {
    EffFitPoint p;
    p.energy = d.energy;
    p.efficiency = d.efficiency;
    p.fracStatUncert = ((d.efficiency_uncert > 0.0f) && (d.efficiency > 0.0f))
                        ? (d.efficiency_uncert / d.efficiency) : 0.05f;
    points.push_back( std::move(p) );
  }
  
  const EffFitResult fit = performEfficiencyFit( points, fcnOrder );
  result = fit.coefs;
  uncerts = fit.uncerts;
  
  return (fit.dof > 0) ? (fit.chi2 / fit.dof) : fit.chi2;
}//performEfficiencyFit(...)
  
}//namespace MakeDrfFit
