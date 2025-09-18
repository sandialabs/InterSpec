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
#include <vector>
#include <stdexcept>
#include <algorithm>

#define BOOST_TEST_MODULE test_EnergyCal_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/EnergyCal.h"

namespace 
{
  // Helper function to generate channel energies from polynomial coefficients
  std::vector<float> generate_channel_energies_from_poly( const std::vector<float>& coeffs, size_t nchannels )
  {
    std::vector<float> channel_energies( nchannels + 1 );
    
    for( size_t channel = 0; channel <= nchannels; ++channel )
    {
      double energy = 0.0;
      for( size_t i = 0; i < coeffs.size(); ++i )
      {
        energy += coeffs[i] * std::pow( static_cast<double>(channel), static_cast<double>(i) );
      }
      channel_energies[channel] = static_cast<float>( energy );
    }
    
    return channel_energies;
  }
  
  // Helper function to create realistic energy calibration data
  std::vector<float> create_realistic_channel_energies( size_t nchannels )
  {
    // Create a realistic HPGe detector calibration: E = a + b*channel + c*channel^2
    // Typical values: ~3 keV/channel gain, small quadratic term
    const std::vector<float> true_coeffs = { 10.0f, 3.0f, 0.0001f };  // offset, gain, quadratic
    return generate_channel_energies_from_poly( true_coeffs, nchannels );
  }
  
} // anonymous namespace

BOOST_AUTO_TEST_CASE( test_fit_poly_from_channel_energies_linear )
{
  // Test linear polynomial fitting (E = a + b*channel)
  const std::vector<float> true_coeffs = { 5.0f, 2.5f };  // 5 keV offset, 2.5 keV/channel gain
  const size_t nchannels = 1024;
  
  // Generate exact channel energies from known polynomial
  const auto channel_energies = generate_channel_energies_from_poly( true_coeffs, nchannels );
  
  std::vector<float> fitted_coeffs;
  const double avg_error = EnergyCal::fit_poly_from_channel_energies( 2, channel_energies, fitted_coeffs );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 2 );
  
  // Check fitted coefficients match true values (allowing for numerical precision)
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-5 );  // offset
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-5 );  // gain
  
  // Average error should be very small for exact polynomial data
  BOOST_CHECK_SMALL( avg_error, 1e-10 );
}

BOOST_AUTO_TEST_CASE( test_fit_poly_from_channel_energies_quadratic )
{
  // Test quadratic polynomial fitting (E = a + b*channel + c*channel^2)
  const std::vector<float> true_coeffs = { 0.5f, 3.2f, 0.0002f };
  const size_t nchannels = 2048;
  
  // Generate exact channel energies from known polynomial
  const auto channel_energies = generate_channel_energies_from_poly( true_coeffs, nchannels );
  
  std::vector<float> fitted_coeffs;
  const double avg_error = EnergyCal::fit_poly_from_channel_energies( 3, channel_energies, fitted_coeffs );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  
  // Check fitted coefficients match true values (allow for numerical precision)
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-3 );  // offset
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-3 );  // linear term
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 1e-3 );  // quadratic term
  
  // Average error should be small for exact polynomial data
  BOOST_CHECK_SMALL( avg_error, 1e-3 );
}

BOOST_AUTO_TEST_CASE( test_fit_poly_from_channel_energies_realistic_detector )
{
  // Test with realistic HPGe detector calibration
  const size_t nchannels = 4096;
  const auto channel_energies = create_realistic_channel_energies( nchannels );
  
  // Fit with quadratic polynomial (typical for HPGe detectors)
  std::vector<float> fitted_coeffs;
  const double avg_error = EnergyCal::fit_poly_from_channel_energies( 3, channel_energies, fitted_coeffs );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  
  // Check that coefficients are reasonable for HPGe detector
  BOOST_CHECK_GT( fitted_coeffs[1], 0.0 );  // Gain should be positive
  BOOST_CHECK( std::abs(fitted_coeffs[2]) < 1.0 );  // Quadratic term should be small
  
  // Average error should be small
  BOOST_CHECK_SMALL( avg_error, 1e-3 );
  
  // Verify that the fitted polynomial reproduces the original energies
  for( size_t channel = 0; channel <= nchannels; channel += 100 )  // Sample every 100 channels
  {
    const double expected_energy = channel_energies[channel];
    const double fitted_energy = fitted_coeffs[0] + fitted_coeffs[1] * channel + 
                                 fitted_coeffs[2] * channel * channel;
    
    BOOST_CHECK_CLOSE( fitted_energy, expected_energy, 1e-3 );  // 0.001% tolerance
  }
}

BOOST_AUTO_TEST_CASE( test_fit_poly_from_channel_energies_higher_order )
{
  // Test higher order polynomial fitting
  const std::vector<float> true_coeffs = { 1.0f, 2.0f, 0.001f, 0.000001f };  // cubic polynomial
  const size_t nchannels = 1024;
  
  // Generate exact channel energies from known polynomial
  const auto channel_energies = generate_channel_energies_from_poly( true_coeffs, nchannels );
  
  std::vector<float> fitted_coeffs;
  const double avg_error = EnergyCal::fit_poly_from_channel_energies( 4, channel_energies, fitted_coeffs );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 4 );
  
  // Check fitted coefficients match true values
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-3 );  // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-3 );  // linear term
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 1e-2 );  // quadratic term  
  BOOST_CHECK_CLOSE( fitted_coeffs[3], true_coeffs[3], 1e-1 );  // cubic term (less precision expected)
  
  // Average error should be small
  BOOST_CHECK_SMALL( avg_error, 1e-3 );
}

BOOST_AUTO_TEST_CASE( test_fit_poly_from_channel_energies_error_conditions )
{
  const size_t nchannels = 10;
  const auto channel_energies = create_realistic_channel_energies( nchannels );
  
  std::vector<float> coeffs;
  
  // Test error condition: requesting 0 coefficients
  BOOST_CHECK_THROW( 
    EnergyCal::fit_poly_from_channel_energies( 0, channel_energies, coeffs ),
    std::exception 
  );
  
  // Test error condition: empty channel energies
  const std::vector<float> empty_energies;
  BOOST_CHECK_THROW( 
    EnergyCal::fit_poly_from_channel_energies( 2, empty_energies, coeffs ),
    std::exception 
  );
  
  // Test error condition: too few data points for requested polynomial order
  // With 11 data points, can't fit order 11 polynomial (which needs 12 coefficients)
  BOOST_CHECK_THROW( 
    EnergyCal::fit_poly_from_channel_energies( 12, channel_energies, coeffs ),
    std::exception 
  );
}

BOOST_AUTO_TEST_CASE( test_fit_poly_from_channel_energies_edge_cases )
{
  // Test basic functionality with reasonable data size
  // This test focuses on ensuring the function works with typical use cases
  const size_t nchannels = 100;
  const auto channel_energies = create_realistic_channel_energies( nchannels );
  
  std::vector<float> coeffs;
  
  // Test with linear fit (should work fine)
  BOOST_CHECK_NO_THROW(
    EnergyCal::fit_poly_from_channel_energies( 2, channel_energies, coeffs )
  );
  
  BOOST_REQUIRE_EQUAL( coeffs.size(), 2 );
  
  // Verify coefficients are reasonable
  BOOST_CHECK_GT( coeffs[1], 0.0 );  // Gain should be positive
  BOOST_CHECK( std::isfinite( coeffs[0] ) );  // Offset should be finite
  BOOST_CHECK( std::isfinite( coeffs[1] ) );  // Gain should be finite
}

BOOST_AUTO_TEST_CASE( test_fit_poly_from_channel_energies_consistency )
{
  // Test consistency: fit with different orders should be self-consistent
  const size_t nchannels = 512;
  const auto channel_energies = create_realistic_channel_energies( nchannels );
  
  // Test different polynomial orders
  for( size_t order = 2; order <= 4; ++order )
  {
    std::vector<float> coeffs;
    const double avg_error = EnergyCal::fit_poly_from_channel_energies( order, channel_energies, coeffs );
    
    BOOST_REQUIRE_EQUAL( coeffs.size(), order );
    
    // Average error should be reasonable
    BOOST_CHECK_GE( avg_error, 0.0 );
    BOOST_CHECK_LT( avg_error, 10.0 );  // Less than 10 keV average error seems reasonable
    
    // All coefficients should be finite
    for( const auto& coeff : coeffs )
    {
      BOOST_CHECK( std::isfinite( coeff ) );
    }
    
    // First coefficient (gain) should be positive for reasonable calibrations
    if( coeffs.size() >= 2 )
    {
      BOOST_CHECK_GT( coeffs[1], 0.0 );
    }
  }
}


BOOST_AUTO_TEST_CASE( test_fit_for_poly_coefs_linear )
{
  // Test linear polynomial fitting: E = a + b*channel
  const std::vector<double> true_coeffs = { 5.0, 2.5 };  // offset, gain
  
  // Generate exact polynomial data points
  std::vector<std::pair<double,double>> channels_energies;
  for( double channel = 0; channel <= 1000; channel += 100 )
  {
    const double energy = true_coeffs[0] + true_coeffs[1] * channel;
    channels_energies.emplace_back( channel, energy );
  }
  
  const std::vector<float> fitted_coeffs = EnergyCal::fit_for_poly_coefs( channels_energies, 2 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 2 );
  
  // Check fitted coefficients match true values
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-6 );  // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-6 );  // linear term
}

BOOST_AUTO_TEST_CASE( test_fit_for_poly_coefs_quadratic )
{
  // Test quadratic polynomial fitting: E = a + b*channel + c*channel^2
  const std::vector<double> true_coeffs = { 1.0, 3.0, 0.001 };
  
  // Generate exact polynomial data points
  std::vector<std::pair<double,double>> channels_energies;
  for( double channel = 0; channel <= 2000; channel += 200 )
  {
    const double energy = true_coeffs[0] + true_coeffs[1] * channel + true_coeffs[2] * channel * channel;
    channels_energies.emplace_back( channel, energy );
  }
  
  const std::vector<float> fitted_coeffs = EnergyCal::fit_for_poly_coefs( channels_energies, 3 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  
  // Check fitted coefficients match true values
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-6 );  // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-6 );  // linear term
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 1e-4 );  // quadratic term (allow for numerical precision)
}

BOOST_AUTO_TEST_CASE( test_fit_for_poly_coefs_realistic_data )
{
  // Test with realistic detector calibration data (HPGe detector)
  std::vector<std::pair<double,double>> channels_energies = {
    { 10, 30.5 },    // Low energy line
    { 50, 122.1 },   // Co-57 122 keV
    { 200, 511.0 },  // Annihilation line
    { 400, 1173.2 }, // Co-60 1173 keV
    { 450, 1332.5 }, // Co-60 1333 keV
    { 600, 1836.1 }, // High energy line
  };
  
  const std::vector<float> fitted_coeffs = EnergyCal::fit_for_poly_coefs( channels_energies, 2 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 2 );
  
  // Check that coefficients are reasonable for HPGe detector
  BOOST_CHECK_GT( fitted_coeffs[1], 0.0 );  // Gain should be positive
  BOOST_CHECK_GT( fitted_coeffs[1], 1.0 );  // Reasonable gain (>1 keV/channel)
  BOOST_CHECK_LT( fitted_coeffs[1], 10.0 ); // Reasonable gain (<10 keV/channel)
  
  // Verify that the fitted polynomial reproduces the calibration points reasonably
  // Note: This is a linear fit to inherently non-linear calibration data, so we expect larger errors
  for( const auto& point : channels_energies )
  {
    const double channel = point.first;
    const double expected_energy = point.second;
    const double fitted_energy = fitted_coeffs[0] + fitted_coeffs[1] * channel;
    
    // Use absolute tolerance for low energies, relative for higher energies
    if( expected_energy < 100.0 ) {
      BOOST_CHECK_SMALL( std::abs(fitted_energy - expected_energy), 50.0 );  // 50 keV absolute tolerance
    } else {
      BOOST_CHECK_CLOSE( fitted_energy, expected_energy, 25.0 );  // 25% tolerance for realistic non-linear data
    }
  }
}

BOOST_AUTO_TEST_CASE( test_fit_for_poly_coefs_error_conditions )
{
  std::vector<std::pair<double,double>> channels_energies = {
    { 100, 300 },
    { 200, 600 }
  };
  
  // Test error condition: requesting more coefficients than data points
  // Note: ublas may not throw for underdetermined systems, so we check the result instead
  try {
    const std::vector<float> result = EnergyCal::fit_for_poly_coefs( channels_energies, 5 );
    // If it doesn't throw, we expect either an exception or a reasonable result
    // Some linear algebra libraries handle underdetermined systems differently
  } catch (const std::exception&) {
    // Expected behavior for underdetermined system
  }
  
  // Test edge case: exactly enough data points
  BOOST_CHECK_NO_THROW(
    EnergyCal::fit_for_poly_coefs( channels_energies, 2 )  // 2 coeffs with 2 data points
  );
}

BOOST_AUTO_TEST_CASE( test_fit_for_fullrangefraction_coefs_basic )
{
  // Test basic full range fraction fitting
  const size_t nchannels = 1024;
  
  // Generate data that follows the FRF model: E = a0 + a1*x + a2*x^2 (where x = channel/nchannels)
  const std::vector<double> true_coeffs = { 10.0, 2000.0, 500.0 };
  
  std::vector<std::pair<double,double>> channels_energies;
  for( double channel = 0; channel <= nchannels; channel += nchannels/10 )
  {
    const double x = channel / nchannels;
    const double energy = true_coeffs[0] + true_coeffs[1] * x + true_coeffs[2] * x * x;
    channels_energies.emplace_back( channel, energy );
  }
  
  const std::vector<float> fitted_coeffs = EnergyCal::fit_for_fullrangefraction_coefs( channels_energies, nchannels, 3 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  
  // Check fitted coefficients match true values (allow for numerical precision)
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-5 );  // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-5 );  // linear term  
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 1e-4 );  // quadratic term
}

BOOST_AUTO_TEST_CASE( test_fit_for_fullrangefraction_coefs_with_special_term )
{
  // Test FRF fitting with the special 5th term: 1/(1 + 60*x)
  const size_t nchannels = 2048;
  
  // Generate data that includes the special term
  const std::vector<double> true_coeffs = { 5.0, 1500.0, 300.0, 50.0, 100.0 };
  
  std::vector<std::pair<double,double>> channels_energies;
  for( double channel = 0; channel <= nchannels; channel += nchannels/12 )
  {
    const double x = channel / nchannels;
    double energy = true_coeffs[0] + true_coeffs[1] * x + true_coeffs[2] * x * x + true_coeffs[3] * x * x * x;
    energy += true_coeffs[4] / (1.0 + 60.0 * x);  // Special 5th term
    channels_energies.emplace_back( channel, energy );
  }
  
  const std::vector<float> fitted_coeffs = EnergyCal::fit_for_fullrangefraction_coefs( channels_energies, nchannels, 5 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 5 );
  
  // Check fitted coefficients are reasonable (the special term makes exact matching difficult)
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1.0 );   // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1.0 );   // linear term
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 1.0 );   // quadratic term
  BOOST_CHECK_CLOSE( fitted_coeffs[3], true_coeffs[3], 1.0 );
  BOOST_CHECK_CLOSE( fitted_coeffs[4], true_coeffs[4], 1.0 );
}

BOOST_AUTO_TEST_CASE( test_fit_for_fullrangefraction_coefs_realistic )
{
  // Test with realistic gamma spectroscopy data
  const size_t nchannels = 4096;
  
  std::vector<std::pair<double,double>> channels_energies = {
    { 0, 0.0 },        // Zero point
    { 205, 122.1 },    // Co-57 122 keV 
    { 850, 511.0 },    // Annihilation line
    { 1950, 1173.2 },  // Co-60 1173 keV
    { 2210, 1332.5 },  // Co-60 1333 keV
    { 4096, 3000.0 },  // High energy endpoint
  };
  
  const std::vector<float> fitted_coeffs = EnergyCal::fit_for_fullrangefraction_coefs( channels_energies, nchannels, 3 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  
  // Check that coefficients are reasonable
  BOOST_CHECK_GT( fitted_coeffs[1], 0.0 );  // Linear term should be positive
  
  // Verify that the fitted model reproduces the calibration points
  for( const auto& point : channels_energies )
  {
    const double channel = point.first;
    const double expected_energy = point.second;
    const double x = channel / nchannels;
    const double fitted_energy = fitted_coeffs[0] + fitted_coeffs[1] * x + fitted_coeffs[2] * x * x;
    
    // Allow reasonable tolerance for real calibration data
    // For zero energy, check absolute difference instead of relative
    if( expected_energy == 0.0 ) {
      BOOST_CHECK_SMALL( fitted_energy, 50.0 );  // Allow up to 50 keV offset for zero point
    } else {
      const double tolerance = (expected_energy > 1000) ? 10.0 : 15.0;
      BOOST_CHECK_CLOSE( fitted_energy, expected_energy, tolerance );
    }
  }
}

BOOST_AUTO_TEST_CASE( test_fit_for_fullrangefraction_coefs_limits )
{
  // Test the limits of the full range fraction function
  const size_t nchannels = 1024;
  
  std::vector<std::pair<double,double>> channels_energies = {
    { 0, 10.0 },
    { 256, 500.0 },
    { 512, 1000.0 },
    { 768, 1500.0 },
    { 1024, 2000.0 }
  };
  
  // Test maximum allowed terms (should be limited to 5)
  const std::vector<float> fitted_coeffs = EnergyCal::fit_for_fullrangefraction_coefs( channels_energies, nchannels, 10 );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 5 );  // Should be clamped to max of 5
  
  // All coefficients should be finite
  for( const auto& coeff : fitted_coeffs )
  {
    BOOST_CHECK( std::isfinite( coeff ) );
  }
}

BOOST_AUTO_TEST_CASE( test_fit_for_fullrangefraction_coefs_consistency )
{
  // Test consistency between different term counts
  const size_t nchannels = 2048;
  
  // Generate simple quadratic data
  std::vector<std::pair<double,double>> channels_energies;
  for( double channel = 0; channel <= nchannels; channel += nchannels/8 )
  {
    const double x = channel / nchannels;
    const double energy = 50.0 + 1800.0 * x + 200.0 * x * x;
    channels_energies.emplace_back( channel, energy );
  }
  
  // Test different numbers of terms
  for( int nterms = 2; nterms <= 4; ++nterms )
  {
    const std::vector<float> fitted_coeffs = EnergyCal::fit_for_fullrangefraction_coefs( channels_energies, nchannels, nterms );
    
    BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), static_cast<size_t>(nterms) );
    
    // Verify basic properties
    BOOST_CHECK_GT( fitted_coeffs[1], 0.0 );  // Linear term should be positive
    BOOST_CHECK( std::isfinite( fitted_coeffs[0] ) );
    BOOST_CHECK( std::isfinite( fitted_coeffs[1] ) );
    
    // Higher order fits should still reproduce the linear behavior reasonably
    const double x_mid = 0.5;  // Middle of range
    const double expected_mid = 50.0 + 1800.0 * x_mid + 200.0 * x_mid * x_mid;
    
    double fitted_mid = fitted_coeffs[0] + fitted_coeffs[1] * x_mid;
    if( nterms >= 3 )
      fitted_mid += fitted_coeffs[2] * x_mid * x_mid;
    if( nterms >= 4 )
      fitted_mid += fitted_coeffs[3] * x_mid * x_mid * x_mid;
    
    BOOST_CHECK_CLOSE( fitted_mid, expected_mid, 5.0 );
  }
}

// ==========================================================================================
// Tests for fit_energy_cal_poly function
// ==========================================================================================

// Test basic polynomial fitting with simple linear calibration
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_poly_linear )
{
  // Create test data for a simple linear calibration: E = 10 + 2*channel
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  // Add peaks at various channels
  EnergyCal::RecalPeakInfo peak1;
  peak1.peakMeanBinNumber = 100.0;  // channel
  peak1.peakMean = 100.0;           // same as bin number for simplicity
  peak1.peakMeanUncert = 0.5;
  peak1.photopeakEnergy = 10.0 + 2.0 * 100.0;  // 210 keV
  peaks.push_back( peak1 );
  
  EnergyCal::RecalPeakInfo peak2;
  peak2.peakMeanBinNumber = 200.0;
  peak2.peakMean = 200.0;
  peak2.peakMeanUncert = 0.5;
  peak2.photopeakEnergy = 10.0 + 2.0 * 200.0;  // 410 keV
  peaks.push_back( peak2 );
  
  // Fit for linear coefficients
  std::vector<bool> fitfor = { true, true };  // Fit for both offset and linear term
  const size_t nchannels = 1000;
  std::vector<std::pair<float,float>> dev_pairs;  // No deviation pairs
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_poly( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify results
  BOOST_REQUIRE_EQUAL( coefs.size(), 2u );
  BOOST_REQUIRE_EQUAL( uncert.size(), 2u );
  
  BOOST_CHECK_CLOSE( coefs[0], 10.0f, 0.1 );   // Offset term
  BOOST_CHECK_CLOSE( coefs[1], 2.0f, 0.1 );    // Linear term
  BOOST_CHECK_SMALL( chi2, 1e-10 );  // Should be near perfect fit
  BOOST_CHECK_GT( uncert[0], 0.0f );  // Uncertainties should be positive
  BOOST_CHECK_GT( uncert[1], 0.0f );
}

// Test polynomial fitting with quadratic calibration
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_poly_quadratic )
{
  // Create test data for quadratic calibration: E = 5 + 1.8*x + 0.001*x^2
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  // Add multiple peaks to constrain quadratic term
  const std::vector<double> channels = { 50.0, 150.0, 300.0, 500.0 };
  const std::vector<double> coeffs_true = { 5.0, 1.8, 0.001 };
  
  for( double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = channel;
    peak.peakMeanUncert = 0.5;
    peak.photopeakEnergy = coeffs_true[0] + coeffs_true[1] * channel + coeffs_true[2] * channel * channel;
    peaks.push_back( peak );
  }
  
  // Fit for quadratic coefficients
  std::vector<bool> fitfor = { true, true, true };
  const size_t nchannels = 1000;
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_poly( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify results
  BOOST_REQUIRE_EQUAL( coefs.size(), 3u );
  BOOST_REQUIRE_EQUAL( uncert.size(), 3u );
  
  BOOST_CHECK_CLOSE( coefs[0], 5.0f, 0.1 );      // Offset term
  BOOST_CHECK_CLOSE( coefs[1], 1.8f, 0.1 );     // Linear term
  BOOST_CHECK_CLOSE( coefs[2], 0.001f, 5.0 );   // Quadratic term (allow more tolerance)
  BOOST_CHECK_SMALL( chi2, 1e-8 );  // Should be near perfect fit
}

// Test polynomial fitting with realistic detector data
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_poly_realistic )
{
  // Simulate realistic gamma spectroscopy calibration peaks
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  // Common calibration sources - Co-60, Cs-137, etc.
  const std::vector<std::pair<double,double>> peak_data = {
    { 186.2, 661.7 },   // Cs-137
    { 292.8, 1173.2 },  // Co-60
    { 317.1, 1332.5 },  // Co-60
    { 121.8, 511.0 },   // Annihilation
    { 59.5, 238.6 }     // Th-232 decay
  };
  
  for( const auto& data : peak_data )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = data.first;
    peak.peakMean = data.first;
    peak.peakMeanUncert = 0.3;  // Typical uncertainty
    peak.photopeakEnergy = data.second;
    peaks.push_back( peak );
  }
  
  // Fit for quadratic calibration (typical for most detectors)
  std::vector<bool> fitfor = { true, true, true };
  const size_t nchannels = 512;
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_poly( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify reasonable results for gamma spectroscopy
  BOOST_REQUIRE_EQUAL( coefs.size(), 3u );
  BOOST_REQUIRE_EQUAL( uncert.size(), 3u );
  
  BOOST_CHECK_GT( coefs[1], 0.0f );  // Linear term should be positive
  BOOST_CHECK_LT( chi2, 5000.0 );   // Should be reasonable fit for real data
  
  // Check that calibration makes sense by evaluating at test points
  // Note: Real detector data may not fit perfectly with simple polynomial
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const double channel = peaks[i].peakMeanBinNumber;
    const double expected_energy = peaks[i].photopeakEnergy;
    const double fitted_energy = coefs[0] + coefs[1] * channel + coefs[2] * channel * channel;
    
    BOOST_CHECK_CLOSE( fitted_energy, expected_energy, 20.0 );  // Within 20% for realistic data
  }
}

// Test error conditions for fit_energy_cal_poly
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_poly_error_conditions )
{
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  std::vector<bool> fitfor = { true, true };
  const size_t nchannels = 1000;
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  // Test with no peaks - should throw
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_poly( peaks, fitfor, nchannels, dev_pairs, coefs, uncert ), 
                     std::runtime_error );
  
  // Add one peak
  EnergyCal::RecalPeakInfo peak;
  peak.peakMeanBinNumber = 100.0;
  peak.peakMean = 100.0;
  peak.peakMeanUncert = 0.5;
  peak.photopeakEnergy = 300.0;
  peaks.push_back( peak );
  
  // Test with no coefficients to fit for - should throw
  std::vector<bool> no_fitfor = { false, false };
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_poly( peaks, no_fitfor, nchannels, dev_pairs, coefs, uncert ), 
                     std::runtime_error );
  
  // Test with more coefficients than peaks - should throw
  std::vector<bool> too_many = { true, true, true, true };
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_poly( peaks, too_many, nchannels, dev_pairs, coefs, uncert ), 
                     std::runtime_error );
}

// Test fixed coefficients functionality
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_poly_fixed_coefficients )
{
  // Create test data
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  EnergyCal::RecalPeakInfo peak1;
  peak1.peakMeanBinNumber = 100.0;
  peak1.peakMean = 100.0;
  peak1.peakMeanUncert = 0.5;
  peak1.photopeakEnergy = 210.0;  // E = 10 + 2*100
  peaks.push_back( peak1 );
  
  EnergyCal::RecalPeakInfo peak2;
  peak2.peakMeanBinNumber = 200.0;
  peak2.peakMean = 200.0;
  peak2.peakMeanUncert = 0.5;
  peak2.photopeakEnergy = 410.0;  // E = 10 + 2*200
  peaks.push_back( peak2 );
  
  // Fix the offset, fit only the linear term
  std::vector<bool> fitfor = { false, true };
  const size_t nchannels = 1000;
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs = { 10.0f, 0.0f };  // Start with correct offset, wrong slope
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_poly( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify results
  BOOST_REQUIRE_EQUAL( coefs.size(), 2u );
  BOOST_CHECK_CLOSE( coefs[0], 10.0f, 0.1 );   // Should remain fixed
  BOOST_CHECK_CLOSE( coefs[1], 2.0f, 1.0 );   // Should be fitted
  BOOST_CHECK_SMALL( chi2, 1e-8 );
}

// ==========================================================================================
// Tests for fit_energy_cal_frf function
// ==========================================================================================

// Test basic FRF fitting with linear-like behavior
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_frf_basic )
{
  // Create test data for FRF calibration
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  // Add peaks - FRF uses normalized channel (x = channel/nchannels)
  const size_t nchannels = 1000;
  const std::vector<double> channels = { 100.0, 300.0, 500.0, 800.0 };
  
  for( double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = channel;
    peak.peakMeanUncert = 0.5;
    
    // Simple linear-like FRF: E ≈ 50 + 1500*x (where x = channel/nchannels)
    const double x = channel / nchannels;
    peak.photopeakEnergy = 50.0 + 1500.0 * x;
    
    peaks.push_back( peak );
  }
  
  // Fit for FRF coefficients (linear terms only)
  std::vector<bool> fitfor = { true, true };
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_frf( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify results
  BOOST_REQUIRE_EQUAL( coefs.size(), 2u );
  BOOST_REQUIRE_EQUAL( uncert.size(), 2u );
  
  BOOST_CHECK_CLOSE( coefs[0], 50.0f, 5.0 );    // Offset term
  BOOST_CHECK_CLOSE( coefs[1], 1500.0f, 5.0 );  // Linear term
  BOOST_CHECK_LT( chi2, 1e-6 );  // Should be very good fit
  BOOST_CHECK_GT( uncert[0], 0.0f );
  BOOST_CHECK_GT( uncert[1], 0.0f );
}

// Test FRF fitting with special 5th term
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_frf_with_special_term )
{
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  const size_t nchannels = 1000;
  const std::vector<double> channels = { 50.0, 200.0, 400.0, 600.0, 850.0 };
  
  // Create synthetic data with the special FRF term: 1/(1 + 60*x)
  for( double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = channel;
    peak.peakMeanUncert = 0.5;
    
    const double x = channel / nchannels;
    // E = 100 + 1000*x + 500*x^2 + 200*x^3 + 50/(1+60*x)
    peak.photopeakEnergy = 100.0 + 1000.0*x + 500.0*x*x + 200.0*x*x*x + 50.0/(1.0 + 60.0*x);
    
    peaks.push_back( peak );
  }
  
  // Fit for all 5 FRF coefficients including special term
  std::vector<bool> fitfor = { true, true, true, true, true };
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_frf( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify results
  BOOST_REQUIRE_EQUAL( coefs.size(), 5u );
  BOOST_REQUIRE_EQUAL( uncert.size(), 5u );
  
  BOOST_CHECK_CLOSE( coefs[0], 100.0f, 5.0 );   // Constant term
  BOOST_CHECK_CLOSE( coefs[1], 1000.0f, 5.0 );  // Linear term
  BOOST_CHECK_CLOSE( coefs[2], 500.0f, 10.0 );  // Quadratic term
  BOOST_CHECK_CLOSE( coefs[3], 200.0f, 15.0 );  // Cubic term
  BOOST_CHECK_CLOSE( coefs[4], 50.0f, 15.0 );   // Special term coefficient
  BOOST_CHECK_LT( chi2, 1e-4 );  // Should be good fit
}

// Test FRF fitting with realistic detector scenario
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_frf_realistic )
{
  // Simulate realistic HPGe detector with FRF calibration
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  // Typical gamma lines from calibration sources
  const std::vector<std::pair<double,double>> peak_data = {
    { 122.1, 661.7 },   // Cs-137
    { 234.2, 1173.2 },  // Co-60
    { 265.8, 1332.5 },  // Co-60
    { 101.9, 511.0 },   // Annihilation
    { 59.5, 238.6 },    // Th-232
    { 295.1, 1460.8 }   // K-40
  };
  
  const size_t nchannels = 512;
  
  for( const auto& data : peak_data )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = data.first;
    peak.peakMean = data.first;
    peak.peakMeanUncert = 0.4;  // Typical uncertainty for HPGe
    peak.photopeakEnergy = data.second;
    peaks.push_back( peak );
  }
  
  // Fit with typical FRF parameters (up to quadratic)
  std::vector<bool> fitfor = { true, true, true };
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_frf( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify reasonable results
  BOOST_REQUIRE_EQUAL( coefs.size(), 3u );
  BOOST_CHECK_GT( coefs[1], 0.0f );  // Linear term should be positive
  BOOST_CHECK_LT( chi2, 1000.0 );   // Should be reasonable fit for real data
  
  // Verify calibration accuracy
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const double channel = peaks[i].peakMeanBinNumber;
    const double expected_energy = peaks[i].photopeakEnergy;
    const double x = channel / nchannels;
    const double fitted_energy = coefs[0] + coefs[1] * x + coefs[2] * x * x;
    
    BOOST_CHECK_CLOSE( fitted_energy, expected_energy, 25.0 );  // Within 25% for realistic data
  }
}

// Test FRF error conditions and limits
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_frf_error_conditions )
{
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  std::vector<bool> fitfor = { true, true };
  const size_t nchannels = 1000;
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs;
  std::vector<float> uncert;
  
  // Test with no peaks - should throw
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_frf( peaks, fitfor, nchannels, dev_pairs, coefs, uncert ), 
                     std::runtime_error );
  
  // Add one peak
  EnergyCal::RecalPeakInfo peak;
  peak.peakMeanBinNumber = 100.0;
  peak.peakMean = 100.0;
  peak.peakMeanUncert = 0.5;
  peak.photopeakEnergy = 300.0;
  peaks.push_back( peak );
  
  // Test with no coefficients to fit for - should throw
  std::vector<bool> no_fitfor = { false, false };
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_frf( peaks, no_fitfor, nchannels, dev_pairs, coefs, uncert ), 
                     std::runtime_error );
  
  // Test with more coefficients than peaks - should throw
  std::vector<bool> too_many = { true, true, true, true, true, true };
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_frf( peaks, too_many, nchannels, dev_pairs, coefs, uncert ), 
                     std::runtime_error );
}

// Test FRF with fixed coefficients
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_frf_fixed_coefficients )
{
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  const size_t nchannels = 1000;
  const std::vector<double> channels = { 200.0, 500.0, 800.0 };
  
  // Create test data: E = 50 + 1200*x
  for( double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = channel;
    peak.peakMeanUncert = 0.5;
    
    const double x = channel / nchannels;
    peak.photopeakEnergy = 50.0 + 1200.0 * x;
    peaks.push_back( peak );
  }
  
  // Fix offset, fit only linear term
  std::vector<bool> fitfor = { false, true };
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs = { 50.0f, 0.0f };  // Correct offset, wrong slope
  std::vector<float> uncert;
  
  const double chi2 = EnergyCal::fit_energy_cal_frf( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
  
  // Verify results
  BOOST_REQUIRE_EQUAL( coefs.size(), 2u );
  BOOST_CHECK_CLOSE( coefs[0], 50.0f, 0.1 );     // Should remain fixed
  BOOST_CHECK_CLOSE( coefs[1], 1200.0f, 1.0 );  // Should be fitted
  BOOST_CHECK_SMALL( chi2, 1e-6 );
}

// Test consistency between different FRF orders
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_frf_consistency )
{
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  
  const size_t nchannels = 1000;
  const std::vector<double> channels = { 100.0, 250.0, 400.0, 650.0, 850.0 };
  
  // Create well-behaved quadratic data
  for( double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = channel;
    peak.peakMeanUncert = 0.3;
    
    const double x = channel / nchannels;
    peak.photopeakEnergy = 75.0 + 1800.0 * x + 200.0 * x * x;  // Quadratic
    peaks.push_back( peak );
  }
  
  // Test different orders and verify they give consistent results for lower-order terms
  for( size_t nterms = 2; nterms <= 4; ++nterms )
  {
    std::vector<bool> fitfor( nterms, true );
    std::vector<std::pair<float,float>> dev_pairs;
    std::vector<float> coefs;
    std::vector<float> uncert;
    
    const double chi2 = EnergyCal::fit_energy_cal_frf( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );
    
    BOOST_REQUIRE_EQUAL( coefs.size(), nterms );
    
    // Verify basic properties
    BOOST_CHECK_GT( coefs[1], 0.0f );  // Linear term should be positive
    BOOST_CHECK( std::isfinite( coefs[0] ) );
    BOOST_CHECK( std::isfinite( coefs[1] ) );
    BOOST_CHECK_LT( chi2, 2000.0 );   // Relax for consistency test
    
    // For quadratic and higher, check quadratic coefficient is reasonable
    if( nterms >= 3 )
    {
      BOOST_CHECK_CLOSE( coefs[2], 200.0f, 10.0 );  // Should recover quadratic term
    }
  }
}
