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
#include <vector>
#include <cassert>
#include <algorithm>

#include "ceres/ceres.h"

#include "InterSpec/BersteinPolynomial.hpp"

namespace BersteinPolynomial
{

namespace
{
  // Number of sample points used for constrained fitting
  constexpr size_t sm_num_fit_points = 50;

  // Stride for Ceres autodiff (maximum number of parameters we expect)
  //. But if we go over this it will just be a ineffiecny, and Ceres should take acare of it
  constexpr int sm_ceres_stride = 8;


  /** Cost functor for Ceres autodiff to fit Bernstein coefficients to target y-values.

   Minimizes the sum of squared differences between the Bernstein polynomial evaluated
   at normalized x-points and the target y-values.
   */
  struct ConstrainedBernsteinFitCost
  {
    // Normalized x-values in [0,1] where we evaluate the polynomial
    std::vector<double> m_x_normalized;

    // Target y-values (clamped to bounds) that we want to fit
    std::vector<double> m_y_target;

    // Number of Bernstein coefficients (degree + 1)
    size_t m_num_coeffs;


    ConstrainedBernsteinFitCost( std::vector<double> x_normalized,
                                 std::vector<double> y_target,
                                 size_t num_coeffs )
      : m_x_normalized( std::move( x_normalized ) ),
        m_y_target( std::move( y_target ) ),
        m_num_coeffs( num_coeffs )
    {
      assert( m_x_normalized.size() == m_y_target.size() );
      assert( m_num_coeffs > 0 );
    }


    template<typename T>
    bool operator()( T const* const* parameters, T* residuals ) const
    {
      const T* const coeffs = parameters[0];

      for( size_t i = 0; i < m_x_normalized.size(); ++i )
      {
        const T x = T( m_x_normalized[i] );
        const T y_eval = BersteinPolynomial::evaluate( x, coeffs, m_num_coeffs );
        residuals[i] = y_eval - T( m_y_target[i] );
      }

      return true;
    }
  };//struct ConstrainedBernsteinFitCost

}//namespace


std::vector<double> constrained_power_series_to_bernstein( const double* power_coeffs,
                                                           const size_t num_coefficients,
                                                           const double x_min,
                                                           const double x_max,
                                                           const double lower_bound,
                                                           const double upper_bound )
{
  if( num_coefficients == 0 )
    throw std::invalid_argument( "BersteinPolynomial::constrained_power_series_to_bernstein:"
                                 " coefficients cannot be empty" );

  if( x_max <= x_min )
    throw std::invalid_argument( "BersteinPolynomial::constrained_power_series_to_bernstein:"
                                 " x_max must be greater than x_min" );

  if( upper_bound <= lower_bound )
    throw std::invalid_argument( "BersteinPolynomial::constrained_power_series_to_bernstein:"
                                 " upper_bound must be greater than lower_bound" );

  // First, try the exact conversion
  std::vector<double> bernstein_coeffs = power_series_to_bernstein( power_coeffs, num_coefficients,
                                                                    x_min, x_max );

  // Check if all coefficients are already within bounds
  bool all_in_bounds = true;
  for( const double coeff : bernstein_coeffs )
  {
    if( (coeff < lower_bound) || (coeff > upper_bound) )
    {
      all_in_bounds = false;
      break;
    }
  }

  // Fast path: if all coefficients are within bounds, return them directly
  if( all_in_bounds )
    return bernstein_coeffs;

  // Slow path: use Ceres to find constrained Bernstein coefficients

  // Sample the original polynomial at uniform points
  const double x_range = x_max - x_min;
  std::vector<double> x_normalized( sm_num_fit_points );
  std::vector<double> y_target( sm_num_fit_points );

  for( size_t i = 0; i < sm_num_fit_points; ++i )
  {
    // Normalized x in [0, 1]
    const double t = static_cast<double>(i) / static_cast<double>(sm_num_fit_points - 1);
    x_normalized[i] = t;

    // Evaluate original polynomial at corresponding x in [x_min, x_max]
    const double x = x_min + t * x_range;
    double y = 0.0;
    double x_pow = 1.0;
    for( size_t j = 0; j < num_coefficients; ++j )
    {
      y += power_coeffs[j] * x_pow;
      x_pow *= x;
    }

    // Clamp y to bounds - this is what we want to fit
    y_target[i] = std::clamp( y, lower_bound, upper_bound );
  }

  // Set up initial values: clamp the exact conversion coefficients to bounds
  std::vector<double> fit_coeffs( num_coefficients );
  for( size_t i = 0; i < num_coefficients; ++i )
  {
    fit_coeffs[i] = std::clamp( bernstein_coeffs[i], lower_bound, upper_bound );
  }

  // Set up Ceres problem
  ceres::Problem problem;

  ConstrainedBernsteinFitCost* cost_functor
    = new ConstrainedBernsteinFitCost( std::move( x_normalized ),
                                       std::move( y_target ),
                                       num_coefficients );

  ceres::DynamicAutoDiffCostFunction<ConstrainedBernsteinFitCost, sm_ceres_stride>* cost_function
    = new ceres::DynamicAutoDiffCostFunction<ConstrainedBernsteinFitCost, sm_ceres_stride>(
        cost_functor, ceres::TAKE_OWNERSHIP );

  cost_function->AddParameterBlock( static_cast<int>(num_coefficients) );
  cost_function->SetNumResiduals( static_cast<int>(sm_num_fit_points) );

  problem.AddResidualBlock( cost_function, nullptr, fit_coeffs.data() );

  // Set bounds on all coefficients
  for( size_t i = 0; i < num_coefficients; ++i )
  {
    problem.SetParameterLowerBound( fit_coeffs.data(), static_cast<int>(i), lower_bound );
    problem.SetParameterUpperBound( fit_coeffs.data(), static_cast<int>(i), upper_bound );
  }

  // Configure solver
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;
  options.max_num_iterations = 100;
  options.function_tolerance = 1e-10;
  options.parameter_tolerance = 1e-10;

  // Solve
  ceres::Solver::Summary summary;
  ceres::Solve( options, &problem, &summary );

  // Check for convergence
  if( summary.termination_type == ceres::FAILURE )
  {
    throw std::runtime_error( "BersteinPolynomial::constrained_power_series_to_bernstein:"
                              " Ceres optimization failed: " + summary.message );
  }

  return fit_coeffs;
}//constrained_power_series_to_bernstein(...)


std::vector<double> constrained_power_series_to_bernstein( const std::vector<double>& power_coeffs,
                                                           const double x_min,
                                                           const double x_max,
                                                           const double lower_bound,
                                                           const double upper_bound )
{
  return constrained_power_series_to_bernstein( power_coeffs.data(), power_coeffs.size(),
                                                x_min, x_max, lower_bound, upper_bound );
}

}//namespace BersteinPolynomial
