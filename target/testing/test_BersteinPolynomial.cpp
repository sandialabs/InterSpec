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
  const double tolerance = 1e-6;  // Tolerance appropriate for higher-order floating-point round-trip operations
  
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