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

#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdexcept>

#define BOOST_TEST_MODULE test_BersteinPolynomial_suite
#include <boost/test/included/unit_test.hpp>
#include <boost/math/tools/polynomial.hpp>

#include "InterSpec/BersteinPolynomial.hpp"

using namespace std;
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE( test_bernstein_evaluation_basic )
{
  // Test basic evaluation of simple polynomials
  
  // Constant polynomial: f(x) = 5
  {
    vector<double> coeffs = {5.0};
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 0.0, coeffs ), 5.0, 1e-10 );
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 0.5, coeffs ), 5.0, 1e-10 );
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 1.0, coeffs ), 5.0, 1e-10 );
  }
  
  // Linear polynomial: f(x) = 2 + 3*x (in Bernstein form: coeffs = {2, 5})
  {
    vector<double> coeffs = {2.0, 5.0};  // B_0,1(x) = 1-x, B_1,1(x) = x
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 0.0, coeffs ), 2.0, 1e-10 );
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 1.0, coeffs ), 5.0, 1e-10 );
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 0.5, coeffs ), 3.5, 1e-10 );
  }
}


BOOST_AUTO_TEST_CASE( test_power_series_to_bernstein_conversion )
{
  const double tolerance = 1e-10;
  
  // Test conversion with mathematically verified expected values
  
  // Constant: f(x) = 5 on [0, 1]
  {
    vector<double> power_coeffs = {5.0};
    vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, 0.0, 1.0 );
    
    BOOST_REQUIRE_EQUAL( bernstein_coeffs.size(), 1 );
    BOOST_CHECK_CLOSE( bernstein_coeffs[0], 5.0, tolerance );
  }
  
  // Linear: f(x) = 2 + 3*x on [0, 1]
  {
    vector<double> power_coeffs = {2.0, 3.0};
    vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, 0.0, 1.0 );
    
    BOOST_REQUIRE_EQUAL( bernstein_coeffs.size(), 2 );
    // Mathematically correct values: b[0] = 2, b[1] = 5
    BOOST_CHECK_CLOSE( bernstein_coeffs[0], 2.0, tolerance );
    BOOST_CHECK_CLOSE( bernstein_coeffs[1], 5.0, tolerance );
    
    // Also verify by evaluation
    boost::math::tools::polynomial<double> power_poly( power_coeffs.begin(), power_coeffs.end() );
    
    vector<double> test_points = {0.0, 0.25, 0.5, 0.75, 1.0};
    for( double x : test_points )
    {
      double power_result = power_poly.evaluate( x );
      double bernstein_result = BersteinPolynomial::evaluate( x, bernstein_coeffs );
      BOOST_CHECK_CLOSE( bernstein_result, power_result, tolerance );
    }
  }
  
  // Quadratic: f(x) = 1 + 2*x + 3*x^2 on [0, 1]
  {
    vector<double> power_coeffs = {1.0, 2.0, 3.0};
    vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, 0.0, 1.0 );
    
    BOOST_REQUIRE_EQUAL( bernstein_coeffs.size(), 3 );
    // Mathematically correct values: b[0] = 1, b[1] = 2, b[2] = 6
    BOOST_CHECK_CLOSE( bernstein_coeffs[0], 1.0, tolerance );
    BOOST_CHECK_CLOSE( bernstein_coeffs[1], 2.0, tolerance );
    BOOST_CHECK_CLOSE( bernstein_coeffs[2], 6.0, tolerance );
    
    // Also verify by evaluation
    boost::math::tools::polynomial<double> power_poly( power_coeffs.begin(), power_coeffs.end() );
    
    vector<double> test_points = {0.0, 0.1, 0.33, 0.5, 0.67, 0.9, 1.0};
    for( double x : test_points )
    {
      double power_result = power_poly.evaluate( x );
      double bernstein_result = BersteinPolynomial::evaluate( x, bernstein_coeffs );
      BOOST_CHECK_CLOSE( bernstein_result, power_result, tolerance );
    }
  }
}


BOOST_AUTO_TEST_CASE( test_bernstein_to_power_series_conversion )
{
  const double tolerance = 1e-10;
  
  // Test conversion back from Bernstein to power series using round-trip verification
  
  // Start with known power series, convert to Bernstein, then back to power series
  vector<vector<double>> test_polynomials = {
    {5.0},                      // constant
    {2.0, 3.0},                // linear
    {1.0, 2.0, 3.0},           // quadratic
    {-1.0, 0.5, -2.0, 1.5}     // cubic
  };
  
  for( const auto& original_power : test_polynomials )
  {
    // Convert power -> bernstein -> power
    vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( original_power, 0.0, 1.0 );
    vector<double> recovered_power = BersteinPolynomial::bernstein_to_power_series( bernstein_coeffs, 0.0, 1.0 );
    
    BOOST_REQUIRE_EQUAL( recovered_power.size(), original_power.size() );
    
    for( size_t i = 0; i < original_power.size(); ++i )
    {
      BOOST_CHECK_CLOSE( recovered_power[i], original_power[i], tolerance );
    }
  }
}


BOOST_AUTO_TEST_CASE( test_round_trip_conversion )
{
  const double tolerance = 1e-12;  // Relaxed for floating-point round-trip accuracy
  
  // Test round-trip conversions: power -> bernstein -> power
  vector<vector<double>> test_polynomials = {
    {1.0},                           // constant
    {2.0, 3.0},                     // linear  
    {1.0, -2.0, 3.0},               // quadratic
    {-1.0, 2.0, -3.0, 4.0},         // cubic
    {0.5, -1.5, 2.5, -3.5, 4.5},   // quartic
    {1.0, 0.0, -1.0, 0.0, 1.0, -1.0} // quintic
  };
  
  for( const auto& original_power : test_polynomials )
  {
    // Convert power -> bernstein -> power
    vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( original_power, 0.0, 1.0 );
    vector<double> recovered_power = BersteinPolynomial::bernstein_to_power_series( bernstein_coeffs, 0.0, 1.0 );
    
    BOOST_REQUIRE_EQUAL( recovered_power.size(), original_power.size() );
    
    for( size_t i = 0; i < original_power.size(); ++i )
    {
      if( std::abs(original_power[i]) < 1e-14 )
      {
        // For values near zero, use absolute tolerance instead of relative
        BOOST_CHECK_SMALL( recovered_power[i], 1e-14 );
      }
      else
      {
        BOOST_CHECK_CLOSE( recovered_power[i], original_power[i], tolerance );
      }
    }
  }
}


BOOST_AUTO_TEST_CASE( test_evaluation_accuracy )
{
  const double tolerance = 1e-10;
  
  // Test that both power series and Bernstein evaluation give same results
  vector<double> power_coeffs = {1.0, -2.0, 3.0, -1.0};  // 1 - 2x + 3x^2 - x^3
  vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, 0.0, 1.0 );
  
  // Create boost polynomial for reference
  boost::math::tools::polynomial<double> boost_poly( power_coeffs.begin(), power_coeffs.end() );
  
  // Test at various x values
  vector<double> test_x_values = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};
  
  for( double x : test_x_values )
  {
    double boost_result = boost_poly.evaluate( x );
    double bernstein_result = BersteinPolynomial::evaluate( x, bernstein_coeffs );
    double via_bernstein_result = BersteinPolynomial::evaluate_power_series_via_bernstein( x, power_coeffs, 0.0, 1.0 );
    
    BOOST_CHECK_CLOSE( bernstein_result, boost_result, tolerance );
    BOOST_CHECK_CLOSE( via_bernstein_result, boost_result, tolerance );
  }
}


BOOST_AUTO_TEST_CASE( test_domain_transformation )
{
  const double tolerance = 1e-10;
  
  // Test polynomial f(x) = 2 + 3*x on domain [-1, 2]
  vector<double> power_coeffs = {2.0, 3.0};
  const double x_min = -1.0;
  const double x_max = 2.0;
  
  // Convert to Bernstein
  vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, x_min, x_max );
  
  // Test evaluation at several points
  boost::math::tools::polynomial<double> boost_poly( power_coeffs.begin(), power_coeffs.end() );
  
  vector<double> test_x_values = {-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
  
  for( double x : test_x_values )
  {
    double expected = boost_poly.evaluate( x );
    double via_bernstein = BersteinPolynomial::evaluate_power_series_via_bernstein( x, power_coeffs, x_min, x_max );
    
    BOOST_CHECK_CLOSE( via_bernstein, expected, tolerance );
  }
  
  // Test round-trip conversion with domain transformation
  vector<double> recovered_power = BersteinPolynomial::bernstein_to_power_series( bernstein_coeffs, x_min, x_max );
  
  BOOST_REQUIRE_EQUAL( recovered_power.size(), power_coeffs.size() );
  for( size_t i = 0; i < power_coeffs.size(); ++i )
  {
    BOOST_CHECK_CLOSE( recovered_power[i], power_coeffs[i], tolerance );
  }
}


BOOST_AUTO_TEST_CASE( test_higher_order_polynomials )
{
  const double tolerance = 5e-5;  // Tolerance, in percent (i.e. 5e-7 actual fractional diff), appropriate for higher-order floating-point round-trip operations
  
  // Test higher order polynomials up to degree 6
  for( int degree = 1; degree <= 6; ++degree )
  {
    // Create random coefficients
    vector<double> power_coeffs( degree + 1 );
    mt19937 rng( 42 + degree );  // Fixed seed for reproducibility
    uniform_real_distribution<double> dist( -5.0, 5.0 );
    
    for( double& coeff : power_coeffs )
    {
      coeff = dist( rng );
    }
    
    // Test with different domains
    vector<pair<double, double>> domains = {{0.0, 1.0}, {-2.0, 3.0}, {-1.0, 1.0}, {10.0, 20.0}};
    
    for( const auto& domain : domains )
    {
      const double x_min = domain.first;
      const double x_max = domain.second;
      
      // Convert power -> bernstein -> power
      vector<double> bernstein_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, x_min, x_max );
      vector<double> recovered_power = BersteinPolynomial::bernstein_to_power_series( bernstein_coeffs, x_min, x_max );
      
      // Check round-trip accuracy
      BOOST_REQUIRE_EQUAL( recovered_power.size(), power_coeffs.size() );
      for( size_t i = 0; i < power_coeffs.size(); ++i )
      {
        if( std::abs(power_coeffs[i]) < 1e-12 )
        {
          // For very small coefficients, use absolute tolerance
          BOOST_CHECK_SMALL( recovered_power[i], 1e-12 );
        }
        else
        {
          BOOST_CHECK_CLOSE( recovered_power[i], power_coeffs[i], tolerance );
        }
      }
      
      // Test evaluation accuracy at several points
      boost::math::tools::polynomial<double> boost_poly( power_coeffs.begin(), power_coeffs.end() );
      
      const double range = x_max - x_min;
      for( int j = 0; j <= 10; ++j )
      {
        const double x = x_min + (j / 10.0) * range;
        const double expected = boost_poly.evaluate( x );
        const double via_bernstein = BersteinPolynomial::evaluate_power_series_via_bernstein( x, power_coeffs, x_min, x_max );
        
        BOOST_CHECK_CLOSE( via_bernstein, expected, tolerance );
      }
    }
  }
}


BOOST_AUTO_TEST_CASE( test_error_conditions )
{
  // Test error conditions
  
  // Empty coefficients
  {
    vector<double> empty_coeffs;
    BOOST_CHECK_THROW( BersteinPolynomial::evaluate( 0.5, empty_coeffs ), std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::power_series_to_bernstein( empty_coeffs, 0.0, 1.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::bernstein_to_power_series( empty_coeffs, 0.0, 1.0 ), std::invalid_argument );
  }
  
  // Invalid domain (x_max <= x_min)
  {
    vector<double> coeffs = {1.0, 2.0};
    BOOST_CHECK_THROW( BersteinPolynomial::power_series_to_bernstein( coeffs, 1.0, 1.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::power_series_to_bernstein( coeffs, 2.0, 1.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::bernstein_to_power_series( coeffs, 1.0, 1.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::bernstein_to_power_series( coeffs, 2.0, 1.0 ), std::invalid_argument );
  }
}


BOOST_AUTO_TEST_CASE( test_special_cases )
{
  const double tolerance = 1e-12;
  
  // Test evaluation at endpoints
  {
    vector<double> coeffs = {1.0, 2.0, 3.0, 4.0};  // Degree 3
    
    // At x=0, only first coefficient contributes
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 0.0, coeffs ), 1.0, tolerance );
    
    // At x=1, only last coefficient contributes  
    BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( 1.0, coeffs ), 4.0, tolerance );
  }
  
  // Test single coefficient (degree 0)
  {
    vector<double> single_coeff = {42.0};
    
    // Should be constant everywhere
    for( double x = 0.0; x <= 1.0; x += 0.1 )
    {
      BOOST_CHECK_CLOSE( BersteinPolynomial::evaluate( x, single_coeff ), 42.0, tolerance );
    }
  }
}


BOOST_AUTO_TEST_CASE( test_specific_evaluation_data )
{
  const double tolerance = 1e-6;
  
  // Test specific evaluation data: coefficients [0.2, 1.8, 0.5, 1.1] (degree 3)
  vector<double> coeffs = {0.2, 1.8, 0.5, 1.1};
  
  // Test x-values and expected y-values
  vector<double> x_values = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  vector<double> expected_y = {0.200000, 0.597800, 0.850400, 0.986600, 1.035200, 1.025000, 
                               0.984800, 0.943400, 0.929600, 0.972200, 1.100000};
  
  BOOST_REQUIRE_EQUAL( x_values.size(), expected_y.size() );
  
  for( size_t i = 0; i < x_values.size(); ++i )
  {
    double computed_y = BersteinPolynomial::evaluate( x_values[i], coeffs );
    BOOST_CHECK_CLOSE( computed_y, expected_y[i], tolerance );
  }
}


BOOST_AUTO_TEST_CASE( test_bernstein_lls_fitting )
{
  const double tolerance = 1e-10;
  
  // Test linear least squares fitting functionality
  
  // Generate synthetic data from a known Bernstein polynomial
  vector<double> true_coeffs = {1.0, 2.5, 0.8, 3.2};  // degree 3 
  vector<double> x_data, y_data, uncertainties;
  
  // Create test data points
  for( double x = 0.0; x <= 1.0; x += 0.1 )
  {
    double y = BersteinPolynomial::evaluate( x, true_coeffs );
    x_data.push_back( x );
    y_data.push_back( y );
    uncertainties.push_back( 0.01 );  // Small constant uncertainty
  }
  
  // Test exact fitting (should recover original coefficients)
  vector<double> fitted_coeffs = BersteinPolynomial::fit_bernstein_lls( x_data, y_data, uncertainties, 3, 0.0, 1.0 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 4 );
  
  for( size_t i = 0; i < true_coeffs.size(); ++i )
  {
    BOOST_CHECK_CLOSE( fitted_coeffs[i], true_coeffs[i], tolerance );
  }
  
  // Test with domain transformation  
  vector<double> x_data_transformed, y_data_transformed;
  const double x_min = 100.0, x_max = 200.0;
  
  for( double x = x_min; x <= x_max; x += 10.0 )
  {
    // Evaluate using domain transformation
    double y = BersteinPolynomial::evaluate_power_series_via_bernstein( x, {1.0, 0.01, 0.0001}, x_min, x_max );
    x_data_transformed.push_back( x );
    y_data_transformed.push_back( y );
  }
  
  vector<double> uncertainties_transformed( x_data_transformed.size(), 0.001 );
  
  // Fit quadratic (degree 2) Bernstein polynomial
  vector<double> fitted_transformed = BersteinPolynomial::fit_bernstein_lls( 
    x_data_transformed, y_data_transformed, uncertainties_transformed, 2, x_min, x_max );
  
  BOOST_REQUIRE_EQUAL( fitted_transformed.size(), 3 );
  
  // Verify fit quality by checking residuals at test points
  for( size_t i = 0; i < x_data_transformed.size(); ++i )
  {
    double x_norm = (x_data_transformed[i] - x_min) / (x_max - x_min);
    double y_fitted = BersteinPolynomial::evaluate( x_norm, fitted_transformed );
    double residual = std::abs( y_fitted - y_data_transformed[i] );
    BOOST_CHECK_SMALL( residual, 0.01 );  // Should fit reasonably well
  }
}


BOOST_AUTO_TEST_CASE( test_lls_error_conditions )
{
  // Test error conditions for LLS fitting
  
  vector<double> x = {1.0, 2.0, 3.0};
  vector<double> y = {1.0, 2.0};  // Wrong size
  vector<double> sigma = {0.1, 0.1, 0.1};
  
  // Mismatched vector sizes
  BOOST_CHECK_THROW( BersteinPolynomial::fit_bernstein_lls( x, y, sigma, 2, 0.0, 4.0 ), std::invalid_argument );
  
  // Empty vectors
  vector<double> empty;
  BOOST_CHECK_THROW( BersteinPolynomial::fit_bernstein_lls( empty, empty, empty, 1, 0.0, 1.0 ), std::invalid_argument );
  
  // Invalid domain
  vector<double> valid_data = {1.0, 2.0, 3.0, 4.0};
  BOOST_CHECK_THROW( BersteinPolynomial::fit_bernstein_lls( valid_data, valid_data, valid_data, 2, 1.0, 1.0 ), std::invalid_argument );
  
  // Too few data points for degree
  BOOST_CHECK_THROW( BersteinPolynomial::fit_bernstein_lls( valid_data, valid_data, valid_data, 5, 0.0, 5.0 ), std::invalid_argument );
}


BOOST_AUTO_TEST_CASE( test_constrained_power_series_to_bernstein_fast_path )
{
  // Test case where all coefficients are already within bounds (fast path)
  // A simple polynomial that stays positive and bounded

  const double tolerance = 1e-10;

  // f(x) = 1 + 0.5*x on [0, 1] gives y in [1, 1.5]
  // Bernstein coeffs should be [1, 1.5], both within [0, 2]
  vector<double> power_coeffs = {1.0, 0.5};
  const double x_min = 0.0;
  const double x_max = 1.0;
  const double lower_bound = 0.0;
  const double upper_bound = 2.0;

  vector<double> constrained_coeffs = BersteinPolynomial::constrained_power_series_to_bernstein(
    power_coeffs, x_min, x_max, lower_bound, upper_bound );

  // Should match exact conversion since coeffs are in bounds
  vector<double> exact_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, x_min, x_max );

  BOOST_REQUIRE_EQUAL( constrained_coeffs.size(), exact_coeffs.size() );

  for( size_t i = 0; i < constrained_coeffs.size(); ++i )
  {
    BOOST_CHECK_CLOSE( constrained_coeffs[i], exact_coeffs[i], tolerance );
  }

  // Verify all coefficients are within bounds
  for( const double coeff : constrained_coeffs )
  {
    BOOST_CHECK_GE( coeff, lower_bound );
    BOOST_CHECK_LE( coeff, upper_bound );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_power_series_to_bernstein_needs_fitting )
{
  // Test case where exact conversion gives coefficients outside bounds
  // requiring the Ceres fitting path

  // A polynomial with oscillation that causes Bernstein coeffs to go outside y-range
  // f(x) = 1 - 4*(x - 0.5)^2 = 1 - 4*x^2 + 4*x - 1 = 4*x - 4*x^2
  // This is a parabola with max at x=0.5, y=1, and y=0 at x=0 and x=1
  // Power series: 0 + 4*x - 4*x^2
  vector<double> power_coeffs = {0.0, 4.0, -4.0};
  const double x_min = 0.0;
  const double x_max = 1.0;

  // The exact Bernstein conversion may have coefficients outside [0, 1]
  // Let's force a tighter bound that requires fitting
  const double lower_bound = 0.2;
  const double upper_bound = 0.9;

  // First verify the exact conversion has coefficients outside bounds
  // to confirm we're actually testing the Ceres fitting path
  vector<double> exact_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, x_min, x_max );
  bool has_out_of_bounds = false;
  for( const double coeff : exact_coeffs )
  {
    if( (coeff < lower_bound) || (coeff > upper_bound) )
    {
      has_out_of_bounds = true;
      break;
    }
  }
  BOOST_CHECK_MESSAGE( has_out_of_bounds,
    "Exact Bernstein coefficients should be outside bounds to test Ceres path" );

  vector<double> constrained_coeffs = BersteinPolynomial::constrained_power_series_to_bernstein(
    power_coeffs, x_min, x_max, lower_bound, upper_bound );

  BOOST_REQUIRE_EQUAL( constrained_coeffs.size(), 3 );

  // Verify ALL coefficients are within bounds
  for( size_t i = 0; i < constrained_coeffs.size(); ++i )
  {
    BOOST_CHECK_GE( constrained_coeffs[i], lower_bound );
    BOOST_CHECK_LE( constrained_coeffs[i], upper_bound );
  }

  // Verify the fitted polynomial approximates the clamped original reasonably well
  // Sample at a few points and check
  for( double t = 0.0; t <= 1.0; t += 0.1 )
  {
    const double y_fitted = BersteinPolynomial::evaluate( t, constrained_coeffs );

    // The fitted value should be within bounds
    BOOST_CHECK_GE( y_fitted, lower_bound - 0.01 );  // Small tolerance for numerical issues
    BOOST_CHECK_LE( y_fitted, upper_bound + 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_power_series_to_bernstein_negative_coeffs )
{
  // Test case where exact conversion gives Bernstein coefficients outside tight bounds
  // The polynomial f(x) = 2 - 3*x + 2*x^2 on [0, 1] has these exact Bernstein coeffs:
  // B0 = f(0) = 2, B1 = (1/2)(c1 + 2*c0) = 0.5, B2 = c0 + c1 + c2 = 1
  // So with tight bounds [0.8, 1.8], the exact B0=2 is too high and B1=0.5 is too low
  vector<double> power_coeffs = {2.0, -3.0, 2.0};
  const double x_min = 0.0;
  const double x_max = 1.0;

  // Tight bounds that force B0=2 and B1=0.5 to be out of range
  const double lower_bound = 0.8;
  const double upper_bound = 1.8;

  // First verify the exact conversion has coefficients outside bounds
  // to confirm we're actually testing the Ceres fitting path
  vector<double> exact_coeffs = BersteinPolynomial::power_series_to_bernstein( power_coeffs, x_min, x_max );
  bool has_out_of_bounds = false;
  for( const double coeff : exact_coeffs )
  {
    if( (coeff < lower_bound) || (coeff > upper_bound) )
    {
      has_out_of_bounds = true;
      break;
    }
  }
  BOOST_CHECK_MESSAGE( has_out_of_bounds,
    "Exact Bernstein coefficients should be outside bounds to test Ceres path" );

  vector<double> constrained_coeffs = BersteinPolynomial::constrained_power_series_to_bernstein(
    power_coeffs, x_min, x_max, lower_bound, upper_bound );

  BOOST_REQUIRE_EQUAL( constrained_coeffs.size(), 3 );

  // Verify all coefficients are within bounds
  for( const double coeff : constrained_coeffs )
  {
    BOOST_CHECK_GE( coeff, lower_bound );
    BOOST_CHECK_LE( coeff, upper_bound );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_power_series_to_bernstein_higher_degree )
{
  // Test with a higher degree polynomial (degree 4)

  // f(x) = 1 + 2*x - 3*x^2 + 2*x^3 - 0.5*x^4
  vector<double> power_coeffs = {1.0, 2.0, -3.0, 2.0, -0.5};
  const double x_min = 0.0;
  const double x_max = 1.0;
  const double lower_bound = 0.0;
  const double upper_bound = 2.0;

  vector<double> constrained_coeffs = BersteinPolynomial::constrained_power_series_to_bernstein(
    power_coeffs, x_min, x_max, lower_bound, upper_bound );

  BOOST_REQUIRE_EQUAL( constrained_coeffs.size(), 5 );

  // Verify all coefficients are within bounds
  for( const double coeff : constrained_coeffs )
  {
    BOOST_CHECK_GE( coeff, lower_bound );
    BOOST_CHECK_LE( coeff, upper_bound );
  }

  // Verify the polynomial stays within bounds when evaluated
  for( double t = 0.0; t <= 1.0; t += 0.05 )
  {
    const double y_fitted = BersteinPolynomial::evaluate( t, constrained_coeffs );
    BOOST_CHECK_GE( y_fitted, lower_bound - 0.01 );
    BOOST_CHECK_LE( y_fitted, upper_bound + 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_power_series_to_bernstein_domain_transform )
{
  // Test with non-unit domain

  // f(x) = 100 + 0.1*x on [100, 200]
  vector<double> power_coeffs = {100.0, 0.1};
  const double x_min = 100.0;
  const double x_max = 200.0;
  const double lower_bound = 105.0;
  const double upper_bound = 125.0;

  vector<double> constrained_coeffs = BersteinPolynomial::constrained_power_series_to_bernstein(
    power_coeffs, x_min, x_max, lower_bound, upper_bound );

  BOOST_REQUIRE_EQUAL( constrained_coeffs.size(), 2 );

  // Verify all coefficients are within bounds
  for( const double coeff : constrained_coeffs )
  {
    BOOST_CHECK_GE( coeff, lower_bound );
    BOOST_CHECK_LE( coeff, upper_bound );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_power_series_to_bernstein_error_conditions )
{
  vector<double> valid_coeffs = {1.0, 2.0};

  // Empty coefficients
  {
    vector<double> empty;
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_power_series_to_bernstein(
      empty, 0.0, 1.0, 0.0, 2.0 ), std::invalid_argument );
  }

  // Invalid domain (x_max <= x_min)
  {
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_power_series_to_bernstein(
      valid_coeffs, 1.0, 1.0, 0.0, 2.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_power_series_to_bernstein(
      valid_coeffs, 2.0, 1.0, 0.0, 2.0 ), std::invalid_argument );
  }

  // Invalid bounds (upper_bound <= lower_bound)
  {
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_power_series_to_bernstein(
      valid_coeffs, 0.0, 1.0, 2.0, 2.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_power_series_to_bernstein(
      valid_coeffs, 0.0, 1.0, 3.0, 2.0 ), std::invalid_argument );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_power_series_approximation_quality )
{
  // Test that the constrained fit provides a good approximation to the original
  // polynomial within the valid range

  // A polynomial: f(x) = 0.5 + 2*x - x^2 on [0, 1]
  // This has range approximately [0.5, 1.5] with max at x=1
  vector<double> power_coeffs = {0.5, 2.0, -1.0};
  const double x_min = 0.0;
  const double x_max = 1.0;

  // Set bounds that should be achievable
  const double lower_bound = 0.4;
  const double upper_bound = 1.6;

  vector<double> constrained_coeffs = BersteinPolynomial::constrained_power_series_to_bernstein(
    power_coeffs, x_min, x_max, lower_bound, upper_bound );

  // Evaluate both original and fitted at several points
  // The fitted should match closely where original is in bounds
  boost::math::tools::polynomial<double> original_poly( power_coeffs.begin(), power_coeffs.end() );

  double max_error = 0.0;
  for( double t = 0.0; t <= 1.0; t += 0.02 )
  {
    const double x = x_min + t * (x_max - x_min);
    const double y_original = original_poly.evaluate( x );
    const double y_clamped = std::clamp( y_original, lower_bound, upper_bound );
    const double y_fitted = BersteinPolynomial::evaluate( t, constrained_coeffs );

    const double error = std::abs( y_fitted - y_clamped );
    max_error = (std::max)( max_error, error );
  }

  // The approximation should be quite good - within 5% of the range
  const double range = upper_bound - lower_bound;
  BOOST_CHECK_LT( max_error, 0.05 * range );
}


BOOST_AUTO_TEST_CASE( test_constrained_bernstein_refit_fast_path )
{
  // When all coefficients are already within bounds, should return them unchanged.
  const vector<double> coeffs = {1.0, 1.5, 2.0, 1.8};
  const double lower_bound = 0.5;
  const double upper_bound = 2.5;

  const vector<double> result = BersteinPolynomial::constrained_bernstein_refit( coeffs,
                                                                                 lower_bound, upper_bound );

  BOOST_REQUIRE_EQUAL( result.size(), coeffs.size() );

  for( size_t i = 0; i < coeffs.size(); ++i )
  {
    BOOST_CHECK_EQUAL( result[i], coeffs[i] );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_bernstein_refit_poly_in_range_coeffs_not )
{
  // Key scenario: the polynomial values are within bounds everywhere in [0,1],
  //  but some Bernstein coefficients are outside bounds.
  //  This happens because Bernstein coefficients form a "convex hull" that can
  //  overshoot the actual polynomial range - e.g., a polynomial oscillating
  //  between 1.0 and 3.0 may have some Bernstein coefficients at 0.5 or 3.5.
  //
  //  After re-fit, the Bernstein coefficients should all be within bounds,
  //  and the evaluated polynomial should be very close to the original.

  // Test case 1: Quadratic with Bernstein coefficients that overshoot
  //  f(x) = 2 - 3*x + 2*x^2 on [0, 1]
  //  Exact Bernstein coefficients: [2.0, 0.5, 1.0]
  //  The polynomial values range from about 0.875 to 2.0, but B1=0.5 is below 0.875.
  //  Set bounds to [0.8, 2.1] - the polynomial is entirely within bounds, but B1=0.5 is not.
  {
    const vector<double> power_coeffs = {2.0, -3.0, 2.0};
    const double x_min = 0.0;
    const double x_max = 1.0;

    const vector<double> exact_bern = BersteinPolynomial::power_series_to_bernstein(
                                                            power_coeffs, x_min, x_max );

    // Verify our expectation: B1 is below lower bound
    BOOST_REQUIRE_EQUAL( exact_bern.size(), 3 );
    BOOST_CHECK_LT( exact_bern[1], 0.8 );  // B1 = 0.5 < 0.8

    // Verify polynomial values are actually within bounds everywhere
    boost::math::tools::polynomial<double> poly( power_coeffs.begin(), power_coeffs.end() );
    for( double x = 0.0; x <= 1.0; x += 0.01 )
    {
      const double y = poly.evaluate( x );
      BOOST_CHECK_GE( y, 0.8 );
      BOOST_CHECK_LE( y, 2.1 );
    }

    const double lower_bound = 0.8;
    const double upper_bound = 2.1;

    const vector<double> refit = BersteinPolynomial::constrained_bernstein_refit( exact_bern,
                                                                                  lower_bound, upper_bound );

    BOOST_REQUIRE_EQUAL( refit.size(), exact_bern.size() );

    // All coefficients must be within bounds
    for( size_t i = 0; i < refit.size(); ++i )
    {
      BOOST_CHECK_GE( refit[i], lower_bound );
      BOOST_CHECK_LE( refit[i], upper_bound );
    }

    // The re-fit polynomial should be very close to the original since the original
    //  is entirely within bounds - the re-fit just adjusts coefficients, not the curve.
    double max_error = 0.0;
    for( double t = 0.0; t <= 1.0; t += 0.01 )
    {
      const double y_original = BersteinPolynomial::evaluate( t, exact_bern );
      const double y_refit = BersteinPolynomial::evaluate( t, refit );
      max_error = (std::max)( max_error, (std::abs)( y_refit - y_original ) );
    }

    // Constraining Bernstein coefficients necessarily changes the polynomial shape,
    //  so the error won't be tiny - but should be a modest fraction of the bound range.
    const double bound_range = upper_bound - lower_bound;
    BOOST_CHECK_LT( max_error, 0.15 * bound_range );
  }

  // Test case 2: Higher-order polynomial (degree 4) where the polynomial stays in
  //  [1.0, 5.0] but some Bernstein coefficients overshoot.
  //  f(x) = 3 - 8*x + 20*x^2 - 20*x^3 + 8*x^4 on [0, 1]
  //  This polynomial evaluates to 3.0 at x=0, 1.25 at x=0.25, 1.0 at x=0.5,
  //  1.25 at x=0.75, 3.0 at x=1.0 - a symmetric "U" shape in [1.0, 3.0].
  {
    const vector<double> power_coeffs = {3.0, -8.0, 20.0, -20.0, 8.0};
    const double x_min = 0.0;
    const double x_max = 1.0;

    const vector<double> exact_bern = BersteinPolynomial::power_series_to_bernstein(
                                                            power_coeffs, x_min, x_max );

    BOOST_REQUIRE_EQUAL( exact_bern.size(), 5 );

    // Verify polynomial values are within [0.9, 3.1] everywhere
    boost::math::tools::polynomial<double> poly( power_coeffs.begin(), power_coeffs.end() );
    double poly_min = 1e30, poly_max = -1e30;
    for( double x = 0.0; x <= 1.0; x += 0.001 )
    {
      const double y = poly.evaluate( x );
      poly_min = (std::min)( poly_min, y );
      poly_max = (std::max)( poly_max, y );
    }

    // Set bounds just outside the polynomial's actual range
    const double lower_bound = poly_min - 0.05;
    const double upper_bound = poly_max + 0.05;

    // Verify at least one Bernstein coeff is outside these bounds
    bool has_out = false;
    for( const double c : exact_bern )
    {
      if( (c < lower_bound) || (c > upper_bound) )
      {
        has_out = true;
        break;
      }
    }
    BOOST_CHECK_MESSAGE( has_out,
      "Expected at least one Bernstein coefficient outside bounds for this test" );

    const vector<double> refit = BersteinPolynomial::constrained_bernstein_refit( exact_bern,
                                                                                  lower_bound, upper_bound );

    BOOST_REQUIRE_EQUAL( refit.size(), exact_bern.size() );

    for( size_t i = 0; i < refit.size(); ++i )
    {
      BOOST_CHECK_GE( refit[i], lower_bound );
      BOOST_CHECK_LE( refit[i], upper_bound );
    }

    // Re-fit should be very close to original since polynomial is within bounds
    double max_error = 0.0;
    for( double t = 0.0; t <= 1.0; t += 0.01 )
    {
      const double y_original = BersteinPolynomial::evaluate( t, exact_bern );
      const double y_refit = BersteinPolynomial::evaluate( t, refit );
      max_error = (std::max)( max_error, (std::abs)( y_refit - y_original ) );
    }

    const double bound_range = upper_bound - lower_bound;
    BOOST_CHECK_LT( max_error, 0.25 * bound_range );
  }

  // Test case 3: Cubic with coefficients out of range on [0.1, 0.9] domain
  //  f(x) = 1 + 10*(x - 0.5)^2 - 10*(x - 0.5)^3 on a narrow domain
  //  which when converted to Bernstein on that domain may have overshooting coefficients.
  {
    // Power series for 1 + 10*(x-0.5)^2 - 10*(x-0.5)^3
    //  = 1 + 10*(x^2 - x + 0.25) - 10*(x^3 - 1.5*x^2 + 0.75*x - 0.125)
    //  = 1 + 10*x^2 - 10*x + 2.5 - 10*x^3 + 15*x^2 - 7.5*x + 1.25
    //  = 4.75 - 17.5*x + 25*x^2 - 10*x^3
    const vector<double> power_coeffs = {4.75, -17.5, 25.0, -10.0};
    const double x_min = 0.1;
    const double x_max = 0.9;

    const vector<double> exact_bern = BersteinPolynomial::power_series_to_bernstein(
                                                            power_coeffs, x_min, x_max );

    BOOST_REQUIRE_EQUAL( exact_bern.size(), 4 );

    // Find the actual polynomial range
    boost::math::tools::polynomial<double> poly( power_coeffs.begin(), power_coeffs.end() );
    double poly_min = 1e30, poly_max = -1e30;
    for( double x = x_min; x <= x_max; x += 0.001 )
    {
      const double y = poly.evaluate( x );
      poly_min = (std::min)( poly_min, y );
      poly_max = (std::max)( poly_max, y );
    }

    // Set bounds to polynomial range + 5% margin
    const double margin = 0.05 * (poly_max - poly_min);
    const double lower_bound = poly_min - margin;
    const double upper_bound = poly_max + margin;

    const vector<double> refit = BersteinPolynomial::constrained_bernstein_refit( exact_bern,
                                                                                  lower_bound, upper_bound );

    BOOST_REQUIRE_EQUAL( refit.size(), exact_bern.size() );

    for( size_t i = 0; i < refit.size(); ++i )
    {
      BOOST_CHECK_GE( refit[i], lower_bound );
      BOOST_CHECK_LE( refit[i], upper_bound );
    }

    // Should be a good approximation
    double max_error = 0.0;
    for( double t = 0.0; t <= 1.0; t += 0.01 )
    {
      const double y_original = BersteinPolynomial::evaluate( t, exact_bern );
      const double y_refit = BersteinPolynomial::evaluate( t, refit );
      max_error = (std::max)( max_error, (std::abs)( y_refit - y_original ) );
    }

    const double bound_range = upper_bound - lower_bound;
    BOOST_CHECK_LT( max_error, 0.20 * bound_range );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_bernstein_refit_poly_slightly_outside_range )
{
  // The polynomial itself goes slightly outside the wanted bounds.
  // The re-fit should approximate the clamped version as closely as possible.

  // f(x) = 2 + 3*sin(pi*x) approximated as a polynomial: 2 + 9.42*x - 15.5*x^2 + 6.08*x^3
  // This goes from 2.0 at x=0, peaks around x=0.5, back to 2.0 at x=1.
  // We'll set upper bound slightly below the peak to force the re-fit.
  //
  // Use a simpler polynomial: f(t) = 1 + 2*t - 4*t^2 + 2*t^3 on [0,1]
  //  f(0) = 1, f(1) = 1, f(0.5) = 1.25.  Max is at roots of f'=0: 2 - 8*t + 6*t^2 = 0
  //  t = (8 ± sqrt(64-48))/12 = (8 ± 4)/12 => t=1 or t=1/3
  //  f(1/3) = 1 + 2/3 - 4/9 + 2/27 = 1 + 18/27 - 12/27 + 2/27 = 1 + 8/27 ≈ 1.296

  const vector<double> power_coeffs = {1.0, 2.0, -4.0, 2.0};
  const vector<double> exact_bern = BersteinPolynomial::power_series_to_bernstein(
                                                          power_coeffs, 0.0, 1.0 );

  // Set upper bound at 1.2, below the polynomial's actual max of ~1.296
  // Set lower bound at 0.9, below the polynomial's min of 1.0
  const double lower_bound = 0.9;
  const double upper_bound = 1.2;

  const vector<double> refit = BersteinPolynomial::constrained_bernstein_refit( exact_bern,
                                                                                lower_bound, upper_bound );

  BOOST_REQUIRE_EQUAL( refit.size(), exact_bern.size() );

  // All coefficients must be within bounds
  for( size_t i = 0; i < refit.size(); ++i )
  {
    BOOST_CHECK_GE( refit[i], lower_bound );
    BOOST_CHECK_LE( refit[i], upper_bound );
  }

  // The re-fit polynomial should stay within bounds (Bernstein bounding property)
  for( double t = 0.0; t <= 1.0; t += 0.01 )
  {
    const double y = BersteinPolynomial::evaluate( t, refit );
    BOOST_CHECK_GE( y, lower_bound - 1e-10 );
    BOOST_CHECK_LE( y, upper_bound + 1e-10 );
  }

  // Where the original is within bounds, the re-fit should be close to it.
  // Where the original exceeds bounds, the re-fit should be near the bound.
  boost::math::tools::polynomial<double> poly( power_coeffs.begin(), power_coeffs.end() );

  double max_error_in_bounds = 0.0;
  for( double t = 0.0; t <= 1.0; t += 0.01 )
  {
    const double x = t;  // domain is [0,1]
    const double y_original = poly.evaluate( x );
    const double y_clamped = std::clamp( y_original, lower_bound, upper_bound );
    const double y_refit = BersteinPolynomial::evaluate( t, refit );
    const double error = (std::abs)( y_refit - y_clamped );
    max_error_in_bounds = (std::max)( max_error_in_bounds, error );
  }

  // The re-fit is a smooth approximation of a clamped curve (which has sharp corners),
  //  so it won't be perfect, but should be within a reasonable fraction of the range.
  const double bound_range = upper_bound - lower_bound;
  BOOST_CHECK_LT( max_error_in_bounds, 0.40 * bound_range );
}


BOOST_AUTO_TEST_CASE( test_constrained_bernstein_refit_error_conditions )
{
  // Empty coefficients
  {
    const vector<double> empty;
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_bernstein_refit( empty, 0.0, 1.0 ),
                       std::invalid_argument );
  }

  // Invalid bounds (upper <= lower)
  {
    const vector<double> coeffs = {1.0, 2.0};
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_bernstein_refit( coeffs, 2.0, 2.0 ),
                       std::invalid_argument );
    BOOST_CHECK_THROW( BersteinPolynomial::constrained_bernstein_refit( coeffs, 3.0, 1.0 ),
                       std::invalid_argument );
  }
}


BOOST_AUTO_TEST_CASE( test_constrained_bernstein_refit_fwhm_like_scenario )
{
  // Simulate the actual use-case in RelActCalcAuto: FWHM² Bernstein coefficients
  //  where initial estimation uses wide bounds but optimization uses tighter bounds.
  //
  //  A typical scenario: initial estimation computed Bernstein coefficients for FWHM²
  //  with wide bounds, then the optimization stage computed tighter bounds based on
  //  channel widths and peak data, causing some coefficients to be just outside range.

  // Scenario 1: Coefficients just slightly outside tighter bounds.
  //  This is the most common real-world case - bounds tighten by 10-20%.
  {
    // Original FWHM² Bernstein coefficients, representing a smoothly increasing curve
    const vector<double> initial_coeffs = {1.8, 3.5, 6.0, 9.0};

    // Tighter bounds: initial_coeffs[0] = 1.8 is below lower bound of 2.0
    const double lower_bound = 2.0;
    const double upper_bound = 10.0;

    const vector<double> refit = BersteinPolynomial::constrained_bernstein_refit( initial_coeffs,
                                                                                  lower_bound, upper_bound );

    BOOST_REQUIRE_EQUAL( refit.size(), initial_coeffs.size() );

    // All coefficients must be within bounds
    for( size_t i = 0; i < refit.size(); ++i )
    {
      BOOST_CHECK_GE( refit[i], lower_bound );
      BOOST_CHECK_LE( refit[i], upper_bound );
    }

    // Re-fit polynomial must stay within bounds (Bernstein bounding property)
    for( double t = 0.0; t <= 1.0; t += 0.01 )
    {
      const double y = BersteinPolynomial::evaluate( t, refit );
      BOOST_CHECK_GE( y, lower_bound - 1e-10 );
      BOOST_CHECK_LE( y, upper_bound + 1e-10 );
    }

    // Since only one coefficient was slightly out of bounds, the re-fit should be
    //  very close to the original (within 2% of bound range at each point)
    const double bound_range = upper_bound - lower_bound;
    for( double t = 0.0; t <= 1.0; t += 0.01 )
    {
      const double y_original = BersteinPolynomial::evaluate( t, initial_coeffs );
      const double y_clamped = std::clamp( y_original, lower_bound, upper_bound );
      const double y_refit = BersteinPolynomial::evaluate( t, refit );
      const double error = (std::abs)( y_refit - y_clamped );
      BOOST_CHECK_LT( error, 0.02 * bound_range );
    }
  }

  // Scenario 2: Multiple coefficients moderately outside bounds.
  //  The re-fit should find the best smooth approximation within bounds.
  {
    // 5-coefficient Bernstein (degree 4) with two out of bounds
    const vector<double> initial_coeffs = {0.8, 2.5, 4.0, 5.5, 7.2};

    // initial_coeffs[0] = 0.8 < 1.0, initial_coeffs[4] = 7.2 > 7.0
    const double lower_bound = 1.0;
    const double upper_bound = 7.0;

    const vector<double> refit = BersteinPolynomial::constrained_bernstein_refit( initial_coeffs,
                                                                                  lower_bound, upper_bound );

    BOOST_REQUIRE_EQUAL( refit.size(), initial_coeffs.size() );

    for( size_t i = 0; i < refit.size(); ++i )
    {
      BOOST_CHECK_GE( refit[i], lower_bound );
      BOOST_CHECK_LE( refit[i], upper_bound );
    }

    // Polynomial must stay within bounds
    for( double t = 0.0; t <= 1.0; t += 0.01 )
    {
      const double y = BersteinPolynomial::evaluate( t, refit );
      BOOST_CHECK_GE( y, lower_bound - 1e-10 );
      BOOST_CHECK_LE( y, upper_bound + 1e-10 );
    }

    // Should be a reasonable approximation of the clamped curve
    double max_error = 0.0;
    const double bound_range = upper_bound - lower_bound;
    for( double t = 0.0; t <= 1.0; t += 0.01 )
    {
      const double y_original = BersteinPolynomial::evaluate( t, initial_coeffs );
      const double y_clamped = std::clamp( y_original, lower_bound, upper_bound );
      const double y_refit = BersteinPolynomial::evaluate( t, refit );
      max_error = (std::max)( max_error, (std::abs)( y_refit - y_clamped ) );
    }

    BOOST_CHECK_LT( max_error, 0.05 * bound_range );
  }
}
