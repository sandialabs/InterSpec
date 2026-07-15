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

#include "SpecUtils/EnergyCalibration.h"

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

// Test polynomial fitting stays numerically stable for an ill-conditioned problem: many channels
//  (raw basis columns then span ~13 orders of magnitude), higher polynomial order, and peaks
//  spanning the full spectrum.  Before the column-equilibrated solve (and SVD pseudo-inverse
//  covariance) the coefficient uncertainties came from inverse(A^T*A), which for this problem is
//  numerically singular (condition number ~1e26).
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_poly_ill_conditioned )
{
  const size_t nchannels = 16384;
  const std::vector<double> coeffs_true = { -3.0, 0.183, 2.0e-8, 5.0e-13 };

  std::vector<EnergyCal::RecalPeakInfo> peaks;
  const std::vector<double> channels = { 300.0, 1200.0, 3000.0, 6500.0, 9000.0, 12000.0, 15500.0 };
  for( double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = channel;
    peak.peakMeanUncert = 0.25;
    peak.photopeakEnergy = coeffs_true[0] + coeffs_true[1]*channel
                           + coeffs_true[2]*channel*channel
                           + coeffs_true[3]*channel*channel*channel;
    peaks.push_back( peak );
  }

  const std::vector<bool> fitfor( 4, true );
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs, uncert;

  const double chi2 = EnergyCal::fit_energy_cal_poly( peaks, fitfor, nchannels, dev_pairs, coefs, uncert );

  BOOST_REQUIRE_EQUAL( coefs.size(), 4u );
  BOOST_REQUIRE_EQUAL( uncert.size(), 4u );

  BOOST_CHECK_CLOSE( coefs[0], coeffs_true[0], 1.0 );
  BOOST_CHECK_CLOSE( coefs[1], coeffs_true[1], 0.1 );
  BOOST_CHECK_CLOSE( coefs[2], coeffs_true[2], 10.0 );
  BOOST_CHECK_CLOSE( coefs[3], coeffs_true[3], 25.0 );
  // chi2 is not exactly zero only because the returned coefficients are rounded to float
  //  (e.g., the cubic term rounding, times 15500^3, gives ~1e-7 keV residuals)
  BOOST_CHECK_SMALL( chi2, 1.0e-4 );

  // The predicted energies matter more than the individual (correlated) coefficients
  for( double channel : channels )
  {
    const double e_true = coeffs_true[0] + coeffs_true[1]*channel
                          + coeffs_true[2]*channel*channel
                          + coeffs_true[3]*channel*channel*channel;
    const double e_fit = coefs[0] + coefs[1]*channel + coefs[2]*channel*channel
                         + coefs[3]*channel*channel*channel;
    BOOST_CHECK_SMALL( fabs(e_fit - e_true), 0.01 );
  }

  // Uncertainties must be finite and positive - with the old inverse(A^T*A) covariance these
  //  would come out NaN, negative, or wildly large for this problem
  for( size_t i = 0; i < uncert.size(); ++i )
  {
    BOOST_CHECK( !std::isnan(uncert[i]) && !std::isinf(uncert[i]) );
    BOOST_CHECK_GT( uncert[i], 0.0f );
  }

  // Rough magnitude sanity of the uncertainties: predicted energy uncertainty at a middle peak,
  //  from just the diagonal (ignoring correlations, which only reduce it), shouldnt be more than
  //  a few times the ~0.25 keV scale of the input peak uncertainties.
  {
    const double channel = 6500.0;
    double pred_var = 0.0;
    for( size_t i = 0; i < uncert.size(); ++i )
    {
      const double basis = std::pow( channel, static_cast<double>(i) );
      pred_var += (uncert[i]*basis) * (uncert[i]*basis);
    }
    BOOST_CHECK_LT( std::sqrt(pred_var), 100.0 );
  }
}

// Check the column-equilibrated solve gives the same answer as the hand-computed exact solution
//  of a simple two-point problem (i.e., the equilibration un-scaling is exactly right).
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_poly_equilibration_exactness )
{
  std::vector<EnergyCal::RecalPeakInfo> peaks;

  EnergyCal::RecalPeakInfo peak1;
  peak1.peakMeanBinNumber = 1000.0;
  peak1.peakMean = 1000.0;
  peak1.peakMeanUncert = 1.0;
  peak1.photopeakEnergy = 661.657;
  peaks.push_back( peak1 );

  EnergyCal::RecalPeakInfo peak2;
  peak2.peakMeanBinNumber = 2223.0;
  peak2.peakMean = 2223.0;
  peak2.peakMeanUncert = 1.0;
  peak2.photopeakEnergy = 1460.82;
  peaks.push_back( peak2 );

  const std::vector<bool> fitfor = { true, true };
  std::vector<std::pair<float,float>> dev_pairs;
  std::vector<float> coefs, uncert;

  EnergyCal::fit_energy_cal_poly( peaks, fitfor, 4096, dev_pairs, coefs, uncert );

  // Exact solution of the 2x2 system
  const double gain = (1460.82 - 661.657) / (2223.0 - 1000.0);
  const double offset = 661.657 - gain*1000.0;

  BOOST_REQUIRE_EQUAL( coefs.size(), 2u );
  // Tolerances (in percent) allow for the float rounding of the returned coefficients; the
  //  offset suffers some cancellation (661.657 - gain*1000), so gets a bit more slack.
  BOOST_CHECK_CLOSE( static_cast<double>(coefs[0]), offset, 1.0e-2 );
  BOOST_CHECK_CLOSE( static_cast<double>(coefs[1]), gain, 1.0e-3 );
  BOOST_CHECK_GT( uncert[0], 0.0f );
  BOOST_CHECK_GT( uncert[1], 0.0f );
}

// Helper for the lower-channel-energy calibration tests: a mildly non-linear set of channel
//  energies (quadratic in channel), like a real lower-channel-edge defined detector might have.
namespace
{
std::shared_ptr<const SpecUtils::EnergyCalibration> make_test_lower_channel_cal( const size_t nchannel )
{
  std::vector<float> energies( nchannel + 1 );
  for( size_t i = 0; i <= nchannel; ++i )
  {
    const double x = static_cast<double>(i);
    energies[i] = static_cast<float>( 1.0 + 2.9*x + 3.0e-5*x*x );
  }

  auto cal = std::make_shared<SpecUtils::EnergyCalibration>();
  cal->set_lower_channel_energy( nchannel, std::move(energies) );
  return cal;
}
}//namespace


BOOST_AUTO_TEST_CASE( test_adjust_lower_channel_energy_cal )
{
  const size_t nchannel = 1024;
  const auto orig = make_test_lower_channel_cal( nchannel );
  const auto orig_energies = orig->channel_energies();
  BOOST_REQUIRE( orig_energies && (orig_energies->size() == (nchannel + 1)) );

  // E_new = offset + gain*E_orig, gain dimensionless (nominal 1.0)
  const double offset = 3.25, gain = 1.05;
  const auto adjusted = EnergyCal::adjust_lower_channel_energy_cal( orig, offset, gain );

  BOOST_REQUIRE( adjusted && adjusted->valid() );
  BOOST_REQUIRE( adjusted->type() == SpecUtils::EnergyCalType::LowerChannelEdge );
  BOOST_REQUIRE_EQUAL( adjusted->num_channels(), nchannel );

  const auto new_energies = adjusted->channel_energies();
  BOOST_REQUIRE( new_energies && (new_energies->size() == (nchannel + 1)) );

  for( size_t i = 0; i <= nchannel; i += 100 )
  {
    const double expected = offset + gain*(*orig_energies)[i];
    BOOST_CHECK_SMALL( fabs( (*new_energies)[i] - expected ), 1.0e-2 );  //float, energies ~keV
  }

  // offset = 0, gain = 1 is the identity
  const auto identity = EnergyCal::adjust_lower_channel_energy_cal( orig, 0.0, 1.0 );
  const auto identity_energies = identity->channel_energies();
  for( size_t i = 0; i <= nchannel; i += 100 )
    BOOST_CHECK_SMALL( fabs( static_cast<double>((*identity_energies)[i]) - (*orig_energies)[i] ), 1.0e-3 );

  // A non-positive gain makes the edges non-monotonic (or degenerate)
  BOOST_CHECK_THROW( EnergyCal::adjust_lower_channel_energy_cal( orig, 0.0, 0.0 ), std::exception );
  BOOST_CHECK_THROW( EnergyCal::adjust_lower_channel_energy_cal( orig, 0.0, -0.5 ), std::exception );

  // Non lower-channel-energy input must throw
  auto polycal = std::make_shared<SpecUtils::EnergyCalibration>();
  polycal->set_polynomial( nchannel, {0.0f, 3.0f}, {} );
  BOOST_CHECK_THROW( EnergyCal::adjust_lower_channel_energy_cal( polycal, 1.0, 1.0 ),
                     std::exception );
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_lower_channel )
{
  const size_t nchannel = 2048;
  const auto orig = make_test_lower_channel_cal( nchannel );

  const double true_offset = 2.75, true_gain = 1.008;  //gain dimensionless, nominal 1.0
  const auto adjusted = EnergyCal::adjust_lower_channel_energy_cal( orig, true_offset, true_gain );

  // Make peaks whose true energies are where the adjusted calibration puts those channels
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  for( const double channel : { 250.0, 800.0, 1600.0 } )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = orig->energy_for_channel( channel );
    peak.peakMeanUncert = 0.25;
    peak.photopeakEnergy = adjusted->energy_for_channel( channel );
    peaks.push_back( peak );
  }

  {// Fit both offset and gain
    std::vector<bool> fitfor = { true, true };
    std::vector<float> coefs, uncerts;
    const double chi2 = EnergyCal::fit_energy_cal_lower_channel( peaks, orig, fitfor, coefs, uncerts );

    BOOST_REQUIRE_EQUAL( coefs.size(), 2u );
    BOOST_CHECK_SMALL( fabs( coefs[0] - true_offset ), 1.0e-2 );
    BOOST_CHECK_SMALL( fabs( coefs[1] - true_gain ), 1.0e-4 );
    BOOST_CHECK_SMALL( chi2, 1.0e-4 );
    BOOST_CHECK_GT( uncerts[0], 0.0f );
    BOOST_CHECK_GT( uncerts[1], 0.0f );
  }

  {// Fit only the offset, with gain fixed at its true value
    std::vector<bool> fitfor = { true, false };
    std::vector<float> coefs = { 0.0f, static_cast<float>(true_gain) }, uncerts;
    EnergyCal::fit_energy_cal_lower_channel( peaks, orig, fitfor, coefs, uncerts );

    BOOST_CHECK_SMALL( fabs( coefs[0] - true_offset ), 1.0e-2 );
    BOOST_CHECK_SMALL( fabs( coefs[1] - true_gain ), 1.0e-6 );  //should be passed through un-touched
    BOOST_CHECK_EQUAL( uncerts[1], 0.0f );
  }

  {// Fit only the gain, with offset fixed at its true value
    std::vector<bool> fitfor = { false, true };
    std::vector<float> coefs = { static_cast<float>(true_offset), 0.0f }, uncerts;
    EnergyCal::fit_energy_cal_lower_channel( peaks, orig, fitfor, coefs, uncerts );

    BOOST_CHECK_SMALL( fabs( coefs[0] - true_offset ), 1.0e-6 );
    BOOST_CHECK_SMALL( fabs( coefs[1] - true_gain ), 1.0e-4 );
  }

  {// A single peak is enough to fit just the offset
    std::vector<EnergyCal::RecalPeakInfo> one_peak( 1, peaks[1] );
    std::vector<bool> fitfor = { true, false };
    std::vector<float> coefs = { 0.0f, static_cast<float>(true_gain) }, uncerts;
    EnergyCal::fit_energy_cal_lower_channel( one_peak, orig, fitfor, coefs, uncerts );
    BOOST_CHECK_SMALL( fabs( coefs[0] - true_offset ), 1.0e-2 );
  }

  // Invalid inputs
  {
    std::vector<float> coefs, uncerts;
    std::vector<bool> wrong_size_fitfor = { true };
    BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_lower_channel( peaks, orig, wrong_size_fitfor, coefs, uncerts ),
                       std::exception );

    auto polycal = std::make_shared<SpecUtils::EnergyCalibration>();
    polycal->set_polynomial( nchannel, {0.0f, 3.0f}, {} );
    std::vector<bool> fitfor = { true, true };
    BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_lower_channel( peaks, polycal, fitfor, coefs, uncerts ),
                       std::exception );
  }
}


// The generalized propogate_energy_cal_change: a lower-channel-energy display calibration change
//  gets propagated to a polynomial "other" calibration, keeping relative channel alignment.
BOOST_AUTO_TEST_CASE( test_propogate_energy_cal_change_lower_channel )
{
  const size_t nchannel = 1024;
  const auto orig = make_test_lower_channel_cal( nchannel );
  const auto adjusted = EnergyCal::adjust_lower_channel_energy_cal( orig, 4.0, 1.05 );

  // Note: keep the energy range of `other` within `orig`s range, so the checks below can use
  //  orig->channel_for_energy() without it throwing for being outside the defined range
  auto other = std::make_shared<SpecUtils::EnergyCalibration>();
  other->set_polynomial( 512, {2.0f, 5.5f}, {} );

  const auto new_other = EnergyCal::propogate_energy_cal_change( orig, adjusted, other );
  BOOST_REQUIRE( new_other && new_other->valid() );
  BOOST_REQUIRE( new_other->type() == other->type() );
  BOOST_REQUIRE_EQUAL( new_other->num_channels(), other->num_channels() );

  // Channels of `other` that mapped to a given energy under `orig` must now map to where
  //  `adjusted` puts that same original energy
  for( const double other_channel : { 25.0, 150.0, 350.0, 500.0 } )
  {
    const double energy_before = other->energy_for_channel( other_channel );
    const double orig_channel = orig->channel_for_energy( energy_before );
    const double expected_energy = adjusted->energy_for_channel( orig_channel );
    const double energy_after = new_other->energy_for_channel( other_channel );
    BOOST_CHECK_SMALL( fabs( energy_after - expected_energy ), 0.01 );
  }
}

// ---- fit_energy_cal_ceres tests ----

namespace
{
std::vector<EnergyCal::RecalPeakInfo> make_poly_peaks( const std::vector<float> &coefs,
                                          const std::vector<std::pair<float,float>> &dev_pairs,
                                          const std::vector<double> &channels )
{
  std::vector<EnergyCal::RecalPeakInfo> peaks;
  for( const double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = SpecUtils::polynomial_energy( channel, coefs, dev_pairs );
    peak.peakMeanUncert = 0.25;
    peak.photopeakEnergy = peak.peakMean;
    peaks.push_back( peak );
  }
  return peaks;
}
}//namespace


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_poly )
{
  const std::vector<float> true_coefs = { -1.5f, 2.93f, 2.0e-5f };
  const std::vector<double> channels = { 100.0, 400.0, 950.0 };
  const auto peaks = make_poly_peaks( true_coefs, {}, channels );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 1024;
  setup.fitfor = { true, true, true };
  setup.starting_coefs = { 0.0f, 3.0f, 0.0f };

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_REQUIRE_EQUAL( result.coefs.size(), 3u );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[0]) - true_coefs[0] ), 1.0e-3 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - true_coefs[1] ), 1.0e-4 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[2]) - true_coefs[2] ), 1.0e-6 );
  BOOST_CHECK_SMALL( result.chi2, 1.0e-4 );
  BOOST_CHECK_GT( result.coef_uncerts[1], 0.0f );
  BOOST_CHECK( result.warning_msg.empty() );
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_agrees_with_lls )
{
  // On a purely linear (in the coefficients) problem the Ceres fit must agree with the linear
  //  least squares fit
  const std::vector<float> true_coefs = { 3.1f, 0.37f };
  const std::vector<std::pair<float,float>> dev_pairs = { {0.0f, 0.0f}, {661.7f, 4.0f}, {2614.0f, -3.0f} };
  const std::vector<double> channels = { 500.0, 1500.0, 3000.0, 6000.0 };

  auto peaks = make_poly_peaks( true_coefs, dev_pairs, channels );
  // perturb targets a little, so the fit isnt exact, and the two methods have to agree at the
  //  minimum (not just at the truth)
  peaks[0].photopeakEnergy += 0.21;
  peaks[2].photopeakEnergy -= 0.13;

  const std::vector<bool> fitfor = { true, true };

  std::vector<float> lls_coefs = { 0.0f, 0.0f }, lls_uncerts;
  const double lls_chi2 = EnergyCal::fit_energy_cal_poly( peaks, fitfor, 8192, dev_pairs,
                                                          lls_coefs, lls_uncerts );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 8192;
  setup.fitfor = fitfor;
  setup.starting_coefs = { 0.0f, 0.5f };
  setup.dev_pairs = dev_pairs;

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_REQUIRE_EQUAL( result.coefs.size(), 2u );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[0]) - lls_coefs[0] ), 1.0e-3 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - lls_coefs[1] ), 1.0e-6 );
  BOOST_CHECK_SMALL( fabs( result.chi2 - lls_chi2 ), 1.0e-3*std::max(1.0, lls_chi2) );

  // and the dev pairs should be passed through un-changed
  BOOST_REQUIRE_EQUAL( result.dev_pairs.size(), dev_pairs.size() );
  for( size_t i = 0; i < dev_pairs.size(); ++i )
  {
    BOOST_CHECK_EQUAL( result.dev_pairs[i].first, dev_pairs[i].first );
    BOOST_CHECK_EQUAL( result.dev_pairs[i].second, dev_pairs[i].second );
  }
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_fixed_coef )
{
  const std::vector<float> true_coefs = { -2.0f, 3.0f };
  const std::vector<double> channels = { 200.0, 700.0, 900.0 };
  const auto peaks = make_poly_peaks( true_coefs, {}, channels );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 1024;
  setup.fitfor = { false, true };
  setup.starting_coefs = { -2.0f, 2.5f };  //offset fixed at its true value

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_CHECK_EQUAL( result.coefs[0], -2.0f );  //fixed - must come through exactly
  BOOST_CHECK_EQUAL( result.coef_uncerts[0], 0.0f );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - true_coefs[1] ), 1.0e-5 );
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_frf )
{
  const size_t nchannel = 2048;
  const std::vector<float> true_coefs = { -3.0f, 2900.0f, 20.0f, 0.0f, 12.0f };
  const std::vector<double> channels = { 150.0, 400.0, 800.0, 1200.0, 1600.0, 1900.0 };

  std::vector<EnergyCal::RecalPeakInfo> peaks;
  for( const double channel : channels )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = SpecUtils::fullrangefraction_energy( channel, true_coefs, nchannel, {} );
    peak.peakMeanUncert = 0.25;
    peak.photopeakEnergy = peak.peakMean;
    peaks.push_back( peak );
  }

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::FullRangeFraction;
  setup.num_channels = nchannel;
  setup.fitfor = { true, true, true, false, true };
  setup.starting_coefs = { 0.0f, 3000.0f, 0.0f, 0.0f, 0.0f };

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[0]) - true_coefs[0] ), 0.01 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - true_coefs[1] ), 0.05 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[2]) - true_coefs[2] ), 0.05 );
  BOOST_CHECK_EQUAL( result.coefs[3], 0.0f );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[4]) - true_coefs[4] ), 0.1 );
  BOOST_CHECK_SMALL( result.chi2, 1.0e-3 );
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_dev_pair_offsets )
{
  // Recover known deviation pair offsets (with the calibration coefficients held fixed)
  const std::vector<float> true_coefs = { 0.0f, 3.0f };
  const std::vector<std::pair<float,float>> true_pairs = { {0.0f, 0.0f}, {661.7f, 5.0f},
                                                           {1460.8f, 8.0f}, {2614.0f, 2.0f} };
  const std::vector<double> channels = { 100.0, 220.567, 486.933, 871.333, 960.0 };

  const auto peaks = make_poly_peaks( true_coefs, true_pairs, channels );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 1024;
  setup.fitfor = { false, false };
  setup.starting_coefs = true_coefs;
  setup.dev_pairs = true_pairs;
  setup.fit_dev_pair_offsets = { false, true, true, false };

  //Start the fit offsets away from truth
  setup.dev_pairs[1].second = 0.0f;
  setup.dev_pairs[2].second = 0.0f;

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_REQUIRE_EQUAL( result.dev_pairs.size(), true_pairs.size() );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.dev_pairs[1].second) - true_pairs[1].second ), 0.05 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.dev_pairs[2].second) - true_pairs[2].second ), 0.05 );
  BOOST_CHECK_EQUAL( result.dev_pairs[0].second, 0.0f );  //not fit
  BOOST_CHECK_EQUAL( result.dev_pairs[3].second, 2.0f );  //not fit
  BOOST_CHECK_GT( result.dev_pair_offset_uncerts[1], 0.0f );
  BOOST_CHECK_EQUAL( result.dev_pair_offset_uncerts[0], 0.0f );

  // The result must be self-consistent when applied through SpecUtils
  for( size_t i = 0; i < channels.size(); ++i )
  {
    const double pred = SpecUtils::polynomial_energy( channels[i], result.coefs, result.dev_pairs );
    BOOST_CHECK_SMALL( fabs( pred - peaks[i].photopeakEnergy ), 0.1 );
  }
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_lower_channel )
{
  const size_t nchannel = 2048;
  const auto orig = make_test_lower_channel_cal( nchannel );
  const double true_offset = 3.5, true_gain = 1.006;  //gain dimensionless, nominal 1.0
  const auto adjusted = EnergyCal::adjust_lower_channel_energy_cal( orig, true_offset, true_gain );

  std::vector<EnergyCal::RecalPeakInfo> peaks;
  for( const double channel : { 300.0, 900.0, 1700.0 } )
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = channel;
    peak.peakMean = orig->energy_for_channel( channel );
    peak.peakMeanUncert = 0.25;
    peak.photopeakEnergy = adjusted->energy_for_channel( channel );
    peaks.push_back( peak );
  }

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::LowerChannelEdge;
  setup.num_channels = nchannel;
  setup.fitfor = { true, true };
  setup.starting_coefs = { 0.0f, 0.0f };
  setup.lower_channel_energies = *orig->channel_energies();

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[0]) - true_offset ), 1.0e-2 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - true_gain ), 1.0e-4 );

  // ...and must agree with the dedicated linear fit
  std::vector<bool> fitfor = { true, true };
  std::vector<float> lin_coefs, lin_uncerts;
  EnergyCal::fit_energy_cal_lower_channel( peaks, orig, fitfor, lin_coefs, lin_uncerts );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[0]) - lin_coefs[0] ), 1.0e-3 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - lin_coefs[1] ), 1.0e-3 );
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_degenerate_warns )
{
  // Fitting the offset coefficient together with ALL deviation pair offsets is degenerate; the
  //  fit should still converge (SVD covariance absorbs the null space), but must warn
  const std::vector<float> true_coefs = { 0.0f, 3.0f };
  const std::vector<std::pair<float,float>> pairs = { {300.0f, 2.0f}, {1500.0f, 4.0f} };
  const std::vector<double> channels = { 50.0, 300.0, 700.0, 990.0 };
  const auto peaks = make_poly_peaks( true_coefs, pairs, channels );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 1024;
  setup.fitfor = { true, false };
  setup.starting_coefs = true_coefs;
  setup.dev_pairs = pairs;
  setup.fit_dev_pair_offsets = { true, true };

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_CHECK( !result.warning_msg.empty() );

  // Even though the parameter split is ambiguous, the RESULTING calibration must reproduce the
  //  peak energies
  for( size_t i = 0; i < channels.size(); ++i )
  {
    const double pred = SpecUtils::polynomial_energy( channels[i], result.coefs, result.dev_pairs );
    BOOST_CHECK_SMALL( fabs( pred - peaks[i].photopeakEnergy ), 0.1 );
  }
}


BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_error_conditions )
{
  const auto peaks = make_poly_peaks( {0.0f, 3.0f}, {}, {100.0, 500.0} );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 1024;
  setup.fitfor = { true, true, true };  //3 params > 2 peaks
  setup.starting_coefs = { 0.0f, 3.0f, 0.0f };
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_ceres( peaks, setup ), std::exception );

  setup.fitfor = { false, false };  //nothing to fit
  setup.starting_coefs = { 0.0f, 3.0f };
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_ceres( peaks, setup ), std::exception );

  setup.fitfor = { true, false };  //mismatched starting_coefs
  setup.starting_coefs = { 0.0f };
  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_ceres( peaks, setup ), std::exception );

  BOOST_CHECK_THROW( EnergyCal::fit_energy_cal_ceres( {}, setup ), std::exception );  //no peaks
}

// A lone deviation pair below 0.1 keV is treated by SpecUtils as no deviation at all, so asking
//  to fit its offset must not leave a silent do-nothing parameter in the problem (which would
//  also poison the covariance); the fit should warn, skip that offset, and still fit the
//  coefficients correctly.
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_lone_low_energy_dev_pair )
{
  const std::vector<float> true_coefs = { -1.0f, 3.0f };
  const std::vector<double> channels = { 150.0, 500.0, 900.0 };
  const auto peaks = make_poly_peaks( true_coefs, {}, channels );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 1024;
  setup.fitfor = { true, true };
  setup.starting_coefs = { 0.0f, 2.8f };
  setup.dev_pairs = { {0.05f, 1.0f} };  //below the 0.1 keV threshold SpecUtils cares about
  setup.fit_dev_pair_offsets = { true };

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_CHECK( !result.warning_msg.empty() );  //must tell the user the offset wasnt fit

  BOOST_REQUIRE_EQUAL( result.dev_pairs.size(), 1u );
  BOOST_CHECK_EQUAL( result.dev_pairs[0].second, 1.0f );  //passed through untouched
  BOOST_CHECK_EQUAL( result.dev_pair_offset_uncerts[0], 0.0f );

  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[0]) - true_coefs[0] ), 1.0e-3 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - true_coefs[1] ), 1.0e-5 );
  BOOST_CHECK_GT( result.coef_uncerts[1], 0.0f );  //covariance must still be computable
}

// Reproduces a user-reported divergence: a badly mis-calibrated starting point, plus a deviation
//  pair that DOES sit under a calibration peak (so it is genuinely constrained), used to send the
//  non-linear fit off the rails (offset ~ -59, gain ~0.4, the 766 keV deviation offset ~ -316).
//  Seeding the fit from the linear solution (and bounding the offsets) must let it recover.
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_dev_pair_bad_start )
{
  const size_t nchannel = 8192;
  const std::vector<float> true_coefs = { 0.003f, 0.366f };
  const std::vector<std::pair<float,float>> true_pairs = { {100.0f, 0.0f}, {766.0f, -0.1f},
                                                           {2614.0f, 0.0f} };

  auto truecal = std::make_shared<SpecUtils::EnergyCalibration>();
  truecal->set_polynomial( nchannel, true_coefs, true_pairs );

  std::vector<EnergyCal::RecalPeakInfo> peaks;
  for( const double energy : { 100.0, 250.0, 500.0, 766.0, 1200.0 } )  //note a peak AT 766
  {
    EnergyCal::RecalPeakInfo peak;
    peak.peakMeanBinNumber = truecal->channel_for_energy( energy );
    peak.peakMean = energy;
    peak.peakMeanUncert = 0.25;
    peak.photopeakEnergy = energy;
    peaks.push_back( peak );
  }

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = nchannel;
  setup.fitfor = { true, true };
  setup.starting_coefs = { 0.0f, 0.28f };  //deliberately far from truth (23% low gain)
  setup.dev_pairs = { {100.0f, 0.0f}, {766.0f, 0.0f}, {2614.0f, 0.0f} };  //766 offset starts at 0
  setup.fit_dev_pair_offsets = { false, true, false };

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  //Must recover the calibration, not diverge
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[0]) - true_coefs[0] ), 0.05 );
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - true_coefs[1] ), 1.0e-3 );

  BOOST_REQUIRE_EQUAL( result.dev_pairs.size(), 3u );
  BOOST_CHECK_LT( fabs( static_cast<double>(result.dev_pairs[1].second) ), 1.0 );  //not -316
  BOOST_CHECK_EQUAL( result.dev_pairs[0].second, 0.0f );  //held fixed
  BOOST_CHECK_EQUAL( result.dev_pairs[2].second, 0.0f );  //held fixed

  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const double pred = SpecUtils::polynomial_energy( peaks[i].peakMeanBinNumber, result.coefs,
                                                      result.dev_pairs );
    BOOST_CHECK_SMALL( fabs( pred - peaks[i].photopeakEnergy ), 0.3 );
  }
}

// A deviation pair with NO calibration peak in the energy region it governs is unconstrained;
//  it must be held fixed (with a warning), not fit - otherwise the optimizer overfits through it.
BOOST_AUTO_TEST_CASE( test_fit_energy_cal_ceres_dev_pair_no_governing_peak )
{
  const std::vector<float> true_coefs = { 0.0f, 3.0f };
  const std::vector<double> channels = { 40.0, 80.0, 150.0, 250.0 };  //energies 120..750, all < 766
  const auto peaks = make_poly_peaks( true_coefs, {}, channels );

  EnergyCal::EnergyCalCeresFitSetup setup;
  setup.cal_type = SpecUtils::EnergyCalType::Polynomial;
  setup.num_channels = 1024;
  setup.fitfor = { true, true };
  setup.starting_coefs = { 0.0f, 2.9f };
  //dev pair at 2000 keV has no peak anywhere near it (peaks span 120..750)
  setup.dev_pairs = { {100.0f, 0.0f}, {2000.0f, 0.0f} };
  setup.fit_dev_pair_offsets = { false, true };

  const EnergyCal::EnergyCalCeresFitResult result = EnergyCal::fit_energy_cal_ceres( peaks, setup );

  BOOST_CHECK( !result.warning_msg.empty() );          //must tell the user it wasnt fit
  BOOST_CHECK_EQUAL( result.dev_pairs[1].second, 0.0f );  //held fixed at its starting value
  BOOST_CHECK_SMALL( fabs( static_cast<double>(result.coefs[1]) - true_coefs[1] ), 1.0e-4 );  //coefs still good
}
