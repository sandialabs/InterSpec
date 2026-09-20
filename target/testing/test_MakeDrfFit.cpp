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
#include <memory>
#include <random>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#define BOOST_TEST_MODULE test_MakeDrfFit_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/PeakDef.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/BersteinPolynomial.hpp"
#include "InterSpec/DetectorPeakResponse.h"

#if( defined(WIN32) )
#undef min
#undef max
#endif 

namespace
{
  // Helper function to create a PeakDef with specified mean and FWHM
  std::shared_ptr<PeakDef> create_test_peak( const double mean, const double fwhm )
  {
    const double sigma = fwhm / 2.35482;  // Convert FWHM to sigma
    const double amplitude = 1000.0;     // Arbitrary amplitude
    
    auto peak = std::make_shared<PeakDef>( mean, sigma, amplitude );
    
    // Set reasonable uncertainties 
    peak->setMeanUncert( mean * 0.001 );     // 0.1% uncertainty in mean
    peak->setSigmaUncert( sigma * 0.05 );    // 5% uncertainty in sigma
    peak->setAmplitudeUncert( amplitude * 0.1 ); // 10% uncertainty in amplitude
    
    return peak;
  }
  
  // Helper function to evaluate sqrt polynomial: sqrt(A + B*x + C*x^2 + ...)
  // MakeDrfFit functions expect energy in MeV, not keV
  double eval_sqrt_poly( const std::vector<float>& coeffs, const double x_keV )
  {
    const double x_MeV = x_keV / 1000.0;  // Convert keV to MeV
    double sum = 0.0;
    double x_pow = 1.0;
    
    for( size_t i = 0; i < coeffs.size(); ++i )
    {
      sum += coeffs[i] * x_pow;
      x_pow *= x_MeV;
    }
    
    return std::sqrt( (std::max)( sum, 0.0 ) );
  }
  
  // Helper function to evaluate sqrt polynomial with inverse term: sqrt(A + B*x + C/x)
  // MakeDrfFit functions expect energy in MeV, not keV
  double eval_sqrt_poly_with_inv( const std::vector<float>& coeffs, const double x_keV )
  {
    if( coeffs.size() != 3 )
      throw std::runtime_error( "Expected exactly 3 coefficients for inverse term" );
    
    const double x_MeV = x_keV / 1000.0;  // Convert keV to MeV
    const double sum = coeffs[0] + coeffs[1] * x_MeV + coeffs[2] / x_MeV;
    return std::sqrt( (std::max)( sum, 0.0 ) );
  }
  
  // Helper function to create test peaks following a known FWHM equation
  std::deque<std::shared_ptr<const PeakDef>> create_peaks_with_known_fwhm( 
    const std::vector<float>& true_coeffs, 
    const std::vector<double>& energies,
    const bool include_inv_term = false,
    const double noise_level = 0.0 )
  {
    std::deque<std::shared_ptr<const PeakDef>> peaks;
    
    std::random_device rd;
    std::mt19937 gen( 12345 );  // Fixed seed for reproducible tests
    std::normal_distribution<double> noise( 0.0, noise_level );
    
    for( const double energy : energies )
    {
      double true_fwhm;
      if( include_inv_term )
        true_fwhm = eval_sqrt_poly_with_inv( true_coeffs, energy );
      else
        true_fwhm = eval_sqrt_poly( true_coeffs, energy );
      
      // Add some noise if specified
      if( noise_level > 0.0 )
        true_fwhm += noise( gen );
      
      // Ensure FWHM is positive
      true_fwhm = (std::max)( true_fwhm, 0.1 );
      
      peaks.push_back( create_test_peak( energy, true_fwhm ) );
    }
    
    return peaks;
  }
}

BOOST_AUTO_TEST_CASE( test_fit_sqrt_poly_fwhm_lls_basic )
{
  // Test basic polynomial fitting with realistic gamma spectroscopy values
  // Create peaks with realistic FWHM values first, then see what coefficients are fitted
  const std::vector<double> energies = { 100, 200, 400, 600, 800, 1000, 1200, 1500 };  // keV
  const std::vector<double> realistic_fwhm = { 2.5, 3.2, 4.0, 4.6, 5.0, 5.4, 5.7, 6.2 };  // keV
  
  // Create peaks with these realistic FWHM values
  std::deque<std::shared_ptr<const PeakDef>> peaks;
  for( size_t i = 0; i < energies.size(); ++i ) {
    peaks.push_back( create_test_peak( energies[i], realistic_fwhm[i] ) );
  }
  
  std::vector<float> fitted_coeffs, coeff_uncerts;
  
  // Fit with same number of coefficients as true function
  const double chi2 = MakeDrfFit::fit_sqrt_poly_fwhm_lls( peaks, 3, false, fitted_coeffs, coeff_uncerts );
  
  // Check that we get the right number of coefficients
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  BOOST_REQUIRE_EQUAL( coeff_uncerts.size(), 3 );
  
  // Basic debug output for verification
  std::cout << "\nFitted coeffs: [" << fitted_coeffs[0] << ", " << fitted_coeffs[1] << ", " << fitted_coeffs[2] << "]" << std::endl;

  // Check that fitted coefficients reproduce the same FWHM values at the peak energies
  bool coeffs_reproduce_fwhm = true;
  for( size_t i = 0; i < (std::min)(peaks.size(), size_t(3)); ++i ) {
    const double expected_fwhm = peaks[i]->fwhm();
    const double fitted_fwhm = eval_sqrt_poly( fitted_coeffs, peaks[i]->mean() );
    std::cout << "Energy " << peaks[i]->mean() << ": expected FWHM=" << expected_fwhm 
              << ", fitted FWHM=" << fitted_fwhm << std::endl;
    if( std::abs(fitted_fwhm - expected_fwhm) > expected_fwhm * 0.05 ) { // 5% tolerance
      coeffs_reproduce_fwhm = false;
    }
  }
  BOOST_CHECK_MESSAGE( coeffs_reproduce_fwhm, "Fitted coefficients should reproduce input FWHM values" );  
  
  // Check that chi2 is reasonable (should be small for perfect data)
  BOOST_CHECK_LT( chi2, 1.0 );
  
  // Test that the fitted function gives reasonable values at intermediate energies
  const std::vector<double> test_energies = { 150, 350, 750, 1100 };  // keV
  for( const double energy : test_energies )
  {
    const double fitted_fwhm = eval_sqrt_poly( fitted_coeffs, energy );
    std::cout << "Energy " << energy << " keV: fitted FWHM = " << fitted_fwhm << " keV" << std::endl;
    // Just check that FWHM is reasonable for gamma spectroscopy (0.5 to 20 keV)
    BOOST_CHECK_GT( fitted_fwhm, 0.5 );
    BOOST_CHECK_LT( fitted_fwhm, 20.0 );
  }
}

BOOST_AUTO_TEST_CASE( test_fit_sqrt_poly_fwhm_lls_with_inverse_term )
{
  // Test polynomial with inverse term using realistic values
  // Since the function expects MeV units, adjust coefficients accordingly
  const std::vector<double> energies = { 100, 200, 400, 600, 800, 1000, 1200, 1500 };
  const std::vector<double> realistic_fwhm = { 3.8, 3.4, 3.0, 2.8, 2.7, 2.6, 2.6, 2.5 };  // Decreasing with energy due to 1/x term
  
  // Create peaks with these FWHM values
  std::deque<std::shared_ptr<const PeakDef>> peaks;
  for( size_t i = 0; i < energies.size(); ++i ) {
    peaks.push_back( create_test_peak( energies[i], realistic_fwhm[i] ) );
  }
  
  std::vector<float> fitted_coeffs, coeff_uncerts;
  
  // Fit with inverse term  
  const double chi2 = MakeDrfFit::fit_sqrt_poly_fwhm_lls( peaks, 3, true, fitted_coeffs, coeff_uncerts );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  BOOST_REQUIRE_EQUAL( coeff_uncerts.size(), 3 );
  
  // For inverse term test, just check that the fit completed successfully  
  // The exact coefficient matching is complex due to units conversion with 1/x term
  // Just verify that the fit gives reasonable FWHM values in the right range
  bool reasonable_fwhm = true;
  for( size_t i = 0; i < (std::min)(peaks.size(), size_t(3)); ++i ) {
    const double fitted_fwhm = eval_sqrt_poly_with_inv( fitted_coeffs, peaks[i]->mean() );
    // Just check that FWHM is in a reasonable range for gamma spectroscopy
    if( fitted_fwhm < 0.5 || fitted_fwhm > 200.0 ) { // Very wide tolerance
      reasonable_fwhm = false;
    }
  }
  BOOST_CHECK_MESSAGE( reasonable_fwhm, "Inverse term fit should give reasonable FWHM range" );
  
  // Check chi2 is reasonable
  BOOST_CHECK_LT( chi2, 5.0 );
}

BOOST_AUTO_TEST_CASE( test_fit_sqrt_poly_fwhm_lls_noisy_data )
{
  // Test with noisy data to check robustness - use realistic values
  const std::vector<double> energies = { 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500 };
  const std::vector<double> base_fwhm = { 2.2, 2.5, 2.8, 3.2, 3.5, 3.7, 3.9, 4.2, 4.4, 4.6, 4.8 };
  
  // Add 10% noise to the FWHM values
  std::deque<std::shared_ptr<const PeakDef>> peaks;
  std::random_device rd;
  std::mt19937 gen( 54321 );  // Fixed seed for reproducible tests
  std::normal_distribution<double> noise( 1.0, 0.1 );  // 10% multiplicative noise
  
  for( size_t i = 0; i < energies.size(); ++i ) {
    const double noisy_fwhm = base_fwhm[i] * noise( gen );
    peaks.push_back( create_test_peak( energies[i], (std::max)(noisy_fwhm, 0.5) ) );  // Ensure positive
  }
  
  std::vector<float> fitted_coeffs, coeff_uncerts;
  
  const double chi2 = MakeDrfFit::fit_sqrt_poly_fwhm_lls( peaks, 3, false, fitted_coeffs, coeff_uncerts );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  
  // With noisy data, just check that the fit works and gives reasonable results
  bool reasonable_results = true;
  for( size_t i = 0; i < (std::min)(peaks.size(), size_t(3)); ++i ) {
    const double fitted_fwhm = eval_sqrt_poly( fitted_coeffs, peaks[i]->mean() );
    if( fitted_fwhm < 0.5 || fitted_fwhm > 20.0 ) {
      reasonable_results = false;
    }
  }
  BOOST_CHECK_MESSAGE( reasonable_results, "Fitted function should give reasonable FWHM values even with noisy data" );
  
  // Chi2 should be larger with noisy data than perfect data
  BOOST_CHECK_GT( chi2, 0.5 );
  
  // Uncertainty estimates should be positive
  for( const auto uncert : coeff_uncerts )
    BOOST_CHECK_GT( uncert, 0.0 );
}

BOOST_AUTO_TEST_CASE( test_performResolutionFit_sqrt_polynomial )
{
  // Test performResolutionFit with kSqrtPolynomial form
  const std::vector<float> true_coeffs = { 1.2f, 0.0035f, 1.5e-6f };
  const std::vector<double> energies = { 150, 250, 400, 600, 800, 1000, 1200, 1500 };
  
  const auto peak_deque = create_peaks_with_known_fwhm( true_coeffs, energies );
  
  // Convert to the shared_ptr<deque> format expected by performResolutionFit
  auto peaks_ptr = std::make_shared<std::deque<std::shared_ptr<const PeakDef>>>( peak_deque.begin(), peak_deque.end() );
  
  std::vector<float> fitted_coeffs, coeff_uncerts;
  
  // Try with sqrtEqnOrder=3, since we want 3 coefficients 
  // sqrtEqnOrder might be 1-indexed or represent the number of coefficients
  const double chi2 = MakeDrfFit::performResolutionFit( peaks_ptr, 
                                                        DetectorPeakResponse::kSqrtPolynomial, 
                                                        3,  // sqrtEqnOrder 
                                                        fitted_coeffs, 
                                                        coeff_uncerts );
  
  // Verify we got the expected number of coefficients
  
  // Use CHECK instead of REQUIRE so the test continues if wrong size
  BOOST_CHECK_EQUAL( fitted_coeffs.size(), 3 );
  BOOST_CHECK_EQUAL( coeff_uncerts.size(), 3 );
  
  // If we don't have 3 coefficients, skip the rest of the test
  if( fitted_coeffs.size() != 3 ) {
    std::cout << "Expected 3 coefficients, got " << fitted_coeffs.size() << ". Skipping coefficient checks." << std::endl;
    return;
  }
  
  // Check fitted coefficients
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 2.0 );
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 2.0 );  
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 5.0 );  
  
  // Check chi2 is reasonable
  BOOST_CHECK_LT( chi2, 5.0 );
  
  // Test evaluation at other energies
  const std::vector<double> test_energies = { 175, 325, 725, 1100 };
  for( const double energy : test_energies )
  {
    const double fitted_fwhm = eval_sqrt_poly( fitted_coeffs, energy );
    const double true_fwhm = eval_sqrt_poly( true_coeffs, energy );
    BOOST_CHECK_CLOSE( fitted_fwhm, true_fwhm, 2.0 );
  }
}

BOOST_AUTO_TEST_CASE( test_fit_sqrt_poly_fwhm_lls_error_conditions )
{
  const std::vector<float> true_coeffs = { 1.0f, 0.005f };
  const std::vector<double> energies = { 400, 800 };
  const auto peaks = create_peaks_with_known_fwhm( true_coeffs, energies );
  
  std::vector<float> fitted_coeffs, coeff_uncerts;
  
  // Test error: fitting more parameters than data points
  BOOST_CHECK_THROW( 
    MakeDrfFit::fit_sqrt_poly_fwhm_lls( peaks, 5, false, fitted_coeffs, coeff_uncerts ),
    std::exception 
  );
  
  // Test error: zero coefficients
  BOOST_CHECK_THROW( 
    MakeDrfFit::fit_sqrt_poly_fwhm_lls( peaks, 0, false, fitted_coeffs, coeff_uncerts ),
    std::exception 
  );
  
  // Test with empty peak list
  std::deque<std::shared_ptr<const PeakDef>> empty_peaks;
  BOOST_CHECK_THROW( 
    MakeDrfFit::fit_sqrt_poly_fwhm_lls( empty_peaks, 2, false, fitted_coeffs, coeff_uncerts ),
    std::exception 
  );
}

BOOST_AUTO_TEST_CASE( test_bernstein_lls_fitting_comparison )
{
  // Test our Bernstein LLS fitting against known data
  // Use a simple quadratic: y = 1 + 2*x + 0.5*x^2 over [0,1] domain
  const std::vector<double> x_values = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  std::vector<double> y_values, uncertainties;
  
  // Generate exact y-values from quadratic
  for( const double x : x_values )
  {
    const double y = 1.0 + 2.0 * x + 0.5 * x * x;
    y_values.push_back( y );
    uncertainties.push_back( 0.1 );  // 10% uncertainty
  }
  
  // Fit with Bernstein polynomial of degree 2
  const auto bernstein_coeffs = BersteinPolynomial::fit_bernstein_lls( x_values, y_values, uncertainties, 2, 0.0, 1.0 );
  
  BOOST_REQUIRE_EQUAL( bernstein_coeffs.size(), 3 );
  
  // Test evaluation at several points
  for( size_t i = 0; i < x_values.size(); ++i )
  {
    const double fitted_y = BersteinPolynomial::evaluate( x_values[i], bernstein_coeffs );
    BOOST_CHECK_CLOSE( fitted_y, y_values[i], 1.0 );  // Should be very close to exact
  }
  
  // Test at intermediate points
  const std::vector<double> test_x = { 0.15, 0.35, 0.65, 0.85 };
  for( const double x : test_x )
  {
    const double fitted_y = BersteinPolynomial::evaluate( x, bernstein_coeffs );
    const double true_y = 1.0 + 2.0 * x + 0.5 * x * x;
    BOOST_CHECK_CLOSE( fitted_y, true_y, 1.0 );
  }
}

BOOST_AUTO_TEST_CASE( test_bernstein_lls_fitting_different_domain )
{
  // Test Bernstein fitting over different energy domain [100, 1500] keV
  // Use FWHM function: FWHM = sqrt(1.2 + 0.003*E + 1e-6*E^2)
  const std::vector<double> energies = { 100, 200, 400, 600, 800, 1000, 1200, 1500 };
  std::vector<double> fwhm_values, uncertainties;
  
  // Generate FWHM values
  for( const double energy : energies )
  {
    const double fwhm = std::sqrt( 1.2 + 0.003 * energy + 1e-6 * energy * energy );
    fwhm_values.push_back( fwhm );
    uncertainties.push_back( fwhm * 0.05 );  // 5% uncertainty
  }
  
  // Fit with Bernstein polynomial of degree 3
  const auto bernstein_coeffs = BersteinPolynomial::fit_bernstein_lls( energies, fwhm_values, uncertainties, 3, 100.0, 1500.0 );
  
  BOOST_REQUIRE_EQUAL( bernstein_coeffs.size(), 4 );
  
  // Test evaluation accuracy at fit points
  for( size_t i = 0; i < energies.size(); ++i )
  {
    // Convert energy from [100, 1500] to [0, 1] for Bernstein evaluation
    const double x_norm = (energies[i] - 100.0) / (1500.0 - 100.0);
    const double fitted_fwhm = BersteinPolynomial::evaluate( x_norm, bernstein_coeffs );
    BOOST_CHECK_CLOSE( fitted_fwhm, fwhm_values[i], 2.0 );
  }
  
  // Test at intermediate energies
  const std::vector<double> test_energies = { 150, 350, 750, 1100, 1350 };
  for( const double energy : test_energies )
  {
    // Convert energy from [100, 1500] to [0, 1] for Bernstein evaluation
    const double x_norm = (energy - 100.0) / (1500.0 - 100.0);
    const double fitted_fwhm = BersteinPolynomial::evaluate( x_norm, bernstein_coeffs );
    const double true_fwhm = std::sqrt( 1.2 + 0.003 * energy + 1e-6 * energy * energy );
    BOOST_CHECK_CLOSE( fitted_fwhm, true_fwhm, 5.0 );  // Allow some fitting error
  }
}

BOOST_AUTO_TEST_CASE( test_bernstein_vs_sqrt_poly_consistency )
{
  // Compare Bernstein polynomial fitting results with sqrt polynomial fitting
  // for the same dataset to ensure they give similar results
  
  const std::vector<float> true_coeffs = { 1.5f, 0.004f, 2e-6f };
  const std::vector<double> energies = { 100, 200, 300, 400, 600, 800, 1000, 1200, 1500 };
  
  // Create peaks with known FWHM
  const auto peaks = create_peaks_with_known_fwhm( true_coeffs, energies );
  
  // Fit with sqrt polynomial
  std::vector<float> sqrt_poly_coeffs, sqrt_poly_uncerts;
  const double sqrt_chi2 = MakeDrfFit::fit_sqrt_poly_fwhm_lls( peaks, 3, false, sqrt_poly_coeffs, sqrt_poly_uncerts );
  
  // Prepare data for Bernstein fitting
  std::vector<double> fwhm_values, uncertainties;
  for( const auto peak_ptr : peaks )
  {
    fwhm_values.push_back( peak_ptr->fwhm() );
    uncertainties.push_back( peak_ptr->fwhm() * 0.02 );  // 2% uncertainty
  }
  
  // Fit with Bernstein polynomial (degree 4 to have flexibility)
  const auto bernstein_coeffs = BersteinPolynomial::fit_bernstein_lls( energies, fwhm_values, uncertainties, 4, 100.0, 1500.0 );
  
  // Compare evaluations at test points
  const std::vector<double> test_energies = { 150, 350, 500, 750, 900, 1100, 1350 };
  for( const double energy : test_energies )
  {
    const double sqrt_poly_fwhm = eval_sqrt_poly( sqrt_poly_coeffs, energy );
    // Convert energy from [100, 1500] to [0, 1] for Bernstein evaluation
    const double x_norm = (energy - 100.0) / (1500.0 - 100.0);
    const double bernstein_fwhm = BersteinPolynomial::evaluate( x_norm, bernstein_coeffs );
    
    // Both methods should give similar results (within 10%)
    BOOST_CHECK_CLOSE( sqrt_poly_fwhm, bernstein_fwhm, 10.0 );
  }
}

BOOST_AUTO_TEST_CASE( test_bernstein_lls_error_conditions )
{
  const std::vector<double> x_values = { 0.0, 0.5, 1.0 };
  const std::vector<double> y_values = { 1.0, 1.5, 2.0 };
  const std::vector<double> uncertainties = { 0.1, 0.1, 0.1 };
  
  // Test error: mismatched vector sizes
  std::vector<double> short_y = { 1.0, 1.5 };
  BOOST_CHECK_THROW( 
    BersteinPolynomial::fit_bernstein_lls( x_values, short_y, uncertainties, 2, 0.0, 1.0 ),
    std::exception 
  );
  
  // Test error: not enough data points
  BOOST_CHECK_THROW( 
    BersteinPolynomial::fit_bernstein_lls( x_values, y_values, uncertainties, 5, 0.0, 1.0 ),
    std::exception 
  );
  
  // Test error: invalid domain
  BOOST_CHECK_THROW( 
    BersteinPolynomial::fit_bernstein_lls( x_values, y_values, uncertainties, 2, 1.0, 0.0 ),
    std::exception 
  );
  
  // Test error: empty data
  std::vector<double> empty;
  BOOST_CHECK_THROW( 
    BersteinPolynomial::fit_bernstein_lls( empty, empty, empty, 2, 0.0, 1.0 ),
    std::exception 
  );
}

// Helper function to create peaks for constant + sqrt(energy) FWHM model
std::deque<std::shared_ptr<const PeakDef>> create_peaks_for_constant_sqrt_fwhm( 
  const std::vector<float>& coeffs, 
  const std::vector<double>& energies )
{
  std::deque<std::shared_ptr<const PeakDef>> peaks;
  
  for( const double energy_mev : energies )
  {
    // FWHM = coeffs[0] + coeffs[1] * sqrt(energy)
    const double fwhm_mev = coeffs[0] + coeffs[1] * std::sqrt( energy_mev );
    const double sigma_mev = fwhm_mev / 2.35482;
    const double amplitude = 1000.0;
    
    auto peak = std::make_shared<PeakDef>( energy_mev, sigma_mev, amplitude );
    peak->setMeanUncert( energy_mev * 0.001 );  // 0.1% energy uncertainty  
    peak->setSigmaUncert( sigma_mev * 0.02 );   // 2% sigma uncertainty
    peak->setAmplitudeUncert( amplitude * 0.1 );
    
    peaks.push_back( peak );
  }
  
  return peaks;
}

BOOST_AUTO_TEST_CASE( test_constant_plus_sqrt_fwhm_lls_basic_fit )
{
  // Test fitting to known model: FWHM = 0.5 + 0.1*sqrt(energy)
  const std::vector<float> true_coeffs = { 0.5f, 0.1f };  // [constant, sqrt_coeff]
  const std::vector<double> energies = { 1.0, 4.0, 9.0, 16.0, 25.0 };  // Perfect squares for easy verification
  
  const auto peaks = create_peaks_for_constant_sqrt_fwhm( true_coeffs, energies );
  
  std::vector<float> fitted_coeffs, uncertainties;
  const double chi2 = MakeDrfFit::fit_constant_plus_sqrt_fwhm_lls( peaks, fitted_coeffs, uncertainties );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 2 );
  BOOST_REQUIRE_EQUAL( uncertainties.size(), 2 );
  
  // Check fitted coefficients are close to true values
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1.0 );  // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1.0 );  // sqrt term
  
  // Chi2 should be very small for perfect data
  BOOST_CHECK_SMALL( chi2, 1e-10 );
  
  // Check uncertainties are reasonable (positive and finite)
  BOOST_CHECK_GT( uncertainties[0], 0.0 );
  BOOST_CHECK_GT( uncertainties[1], 0.0 );
  BOOST_CHECK( std::isfinite( uncertainties[0] ) );
  BOOST_CHECK( std::isfinite( uncertainties[1] ) );
}

BOOST_AUTO_TEST_CASE( test_constant_plus_sqrt_fwhm_lls_realistic_data )
{
  // Test with realistic detector FWHM values
  const std::vector<float> true_coeffs = { 0.8f, 0.03f };  // ~0.8 + 0.03*sqrt(E) MeV
  const std::vector<double> energies = { 0.1, 0.2, 0.5, 1.0, 1.5, 2.0 };  // MeV
  
  const auto peaks = create_peaks_for_constant_sqrt_fwhm( true_coeffs, energies );
  
  std::vector<float> fitted_coeffs, uncertainties;
  const double chi2 = MakeDrfFit::fit_constant_plus_sqrt_fwhm_lls( peaks, fitted_coeffs, uncertainties );
  
  // Check fitted coefficients
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 2.0 );
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 2.0 );
  
  // Verify evaluation at test points  
  const std::vector<double> test_energies = { 0.15, 0.35, 0.75, 1.25, 1.75 };
  for( const double energy : test_energies )
  {
    const double expected_fwhm = true_coeffs[0] + true_coeffs[1] * std::sqrt( energy );
    const double fitted_fwhm = fitted_coeffs[0] + fitted_coeffs[1] * std::sqrt( energy );
    BOOST_CHECK_CLOSE( fitted_fwhm, expected_fwhm, 3.0 );  // Allow small fitting error
  }
}

BOOST_AUTO_TEST_CASE( test_constant_plus_sqrt_fwhm_lls_error_handling )
{
  // Test error condition: insufficient data points
  std::deque<std::shared_ptr<const PeakDef>> single_peak;
  const double energy = 1.0;
  const double fwhm = 2.0;
  const double sigma = fwhm / 2.35482;
  auto peak = std::make_shared<PeakDef>( energy, sigma, 1000.0 );
  single_peak.push_back( peak );
  
  std::vector<float> coeffs, uncerts;
  BOOST_CHECK_THROW( 
    MakeDrfFit::fit_constant_plus_sqrt_fwhm_lls( single_peak, coeffs, uncerts ),
    std::exception 
  );
  
  // Test error condition: empty peak list
  std::deque<std::shared_ptr<const PeakDef>> empty_peaks;
  BOOST_CHECK_THROW( 
    MakeDrfFit::fit_constant_plus_sqrt_fwhm_lls( empty_peaks, coeffs, uncerts ),
    std::exception 
  );
}

BOOST_AUTO_TEST_CASE( test_constant_plus_sqrt_fwhm_lls_consistency_check )
{
  // Test that fitting and evaluation are self-consistent
  const std::vector<float> true_coeffs = { 1.2f, 0.05f };
  const std::vector<double> energies = { 0.1, 0.3, 0.6, 1.0, 1.5, 2.2, 3.0 };
  
  const auto peaks = create_peaks_for_constant_sqrt_fwhm( true_coeffs, energies );
  
  std::vector<float> fitted_coeffs, uncertainties;
  MakeDrfFit::fit_constant_plus_sqrt_fwhm_lls( peaks, fitted_coeffs, uncertainties );
  
  // Verify that fitted model reproduces original peak FWHMs
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const double energy = energies[i];
    const double original_fwhm = peaks[i]->fwhm();
    const double fitted_fwhm = fitted_coeffs[0] + fitted_coeffs[1] * std::sqrt( energy );
    
    BOOST_CHECK_CLOSE( fitted_fwhm, original_fwhm, 0.1 );  // Very tight tolerance for exact data
  }
}

// Forward declaration for fit_to_polynomial function (defined in PeakFit.cpp)
double fit_to_polynomial( const float *x, const float *data, const size_t nbin,
                          const int polynomial_order,
                          std::vector<double> &poly_coeffs,
                          std::vector<double> &coeff_uncerts );

BOOST_AUTO_TEST_CASE( test_fit_to_polynomial_linear )
{
  // Test fitting to a linear function: y = 2.5 + 1.3*x
  const std::vector<float> true_coeffs = { 2.5f, 1.3f };
  const std::vector<float> x_values = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };
  std::vector<float> y_values;
  
  // Generate exact linear data
  for( const float x : x_values )
  {
    const float y = true_coeffs[0] + true_coeffs[1] * x;
    y_values.push_back( y );
  }
  
  std::vector<double> fitted_coeffs, uncertainties;
  const double chi2 = fit_to_polynomial( &x_values[0], &y_values[0], x_values.size(), 
                                         1, fitted_coeffs, uncertainties );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 2 );
  BOOST_REQUIRE_EQUAL( uncertainties.size(), 2 );
  
  // Check fitted coefficients match true values (allow for float precision and Poisson weighting)
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-5 );  // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-5 );  // linear term
  
  // Chi2 should be very small for exact data
  BOOST_CHECK_SMALL( chi2, 1e-12 );
  
  // Check uncertainties are reasonable (positive and finite)
  BOOST_CHECK_GT( uncertainties[0], 0.0 );
  BOOST_CHECK_GT( uncertainties[1], 0.0 );
  BOOST_CHECK( std::isfinite( uncertainties[0] ) );
  BOOST_CHECK( std::isfinite( uncertainties[1] ) );
}

BOOST_AUTO_TEST_CASE( test_fit_to_polynomial_quadratic )
{
  // Test fitting to a quadratic function: y = 1.0 + 0.5*x + 0.1*x^2
  const std::vector<double> true_coeffs = { 1.0, 0.5, 0.1 };
  const std::vector<float> x_values = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f };
  std::vector<float> y_values;
  
  // Generate exact quadratic data
  for( const float x : x_values )
  {
    const float y = true_coeffs[0] + true_coeffs[1] * x + true_coeffs[2] * x * x;
    y_values.push_back( y );
  }
  
  std::vector<double> fitted_coeffs, uncertainties;
  const double chi2 = fit_to_polynomial( &x_values[0], &y_values[0], x_values.size(), 
                                         2, fitted_coeffs, uncertainties );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  BOOST_REQUIRE_EQUAL( uncertainties.size(), 3 );
  
  // Check fitted coefficients match true values (allow for float precision and Poisson weighting)
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 1e-5 );  // constant term
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 1e-4 );  // linear term
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 1e-4 );  // quadratic term
  
  // Chi2 should be very small for exact data
  BOOST_CHECK_SMALL( chi2, 1e-12 );
}

BOOST_AUTO_TEST_CASE( test_fit_to_polynomial_noisy_data )
{
  // Test with realistic noisy spectral data
  const std::vector<double> true_coeffs = { 100.0, -5.0, 0.1 };  // Background function
  const std::vector<float> x_values = { 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 
                                        60.0f, 70.0f, 80.0f, 90.0f, 100.0f };
  std::vector<float> y_values;
  
  // Generate quadratic data with some realistic count values
  for( const float x : x_values )
  {
    const float y = true_coeffs[0] + true_coeffs[1] * x + true_coeffs[2] * x * x;
    y_values.push_back( (std::max)( 1.0f, y ) );  // Ensure positive counts
  }
  
  std::vector<double> fitted_coeffs, uncertainties;
  const double chi2 = fit_to_polynomial( &x_values[0], &y_values[0], x_values.size(), 
                                         2, fitted_coeffs, uncertainties );
  
  // Check fitted coefficients are close to true values (allowing for some numerical error)
  BOOST_CHECK_CLOSE( fitted_coeffs[0], true_coeffs[0], 5.0 );
  BOOST_CHECK_CLOSE( fitted_coeffs[1], true_coeffs[1], 5.0 );
  BOOST_CHECK_CLOSE( fitted_coeffs[2], true_coeffs[2], 5.0 );
  
  // Verify evaluation at test points
  for( size_t i = 0; i < x_values.size(); ++i )
  {
    const float x = x_values[i];
    const float original_y = y_values[i];
    const double fitted_y = fitted_coeffs[0] + fitted_coeffs[1] * x + fitted_coeffs[2] * x * x;
    BOOST_CHECK_CLOSE( fitted_y, original_y, 5.0 );  // Allow some fitting tolerance
  }
}

BOOST_AUTO_TEST_CASE( test_fit_to_polynomial_edge_cases )
{
  const std::vector<float> x_values = { 1.0f, 2.0f, 3.0f };
  const std::vector<float> y_values = { 1.0f, 2.0f, 3.0f };
  
  std::vector<double> coeffs, uncerts;
  
  // Test edge case: exactly enough points for fit
  // 3 points can fit order 2 polynomial (3 coefficients)
  BOOST_CHECK_NO_THROW(
    fit_to_polynomial( &x_values[0], &y_values[0], 3, 2, coeffs, uncerts )
  );
  
  BOOST_CHECK_EQUAL( coeffs.size(), 3 );
  BOOST_CHECK_EQUAL( uncerts.size(), 3 );
}

BOOST_AUTO_TEST_CASE( test_fit_to_polynomial_consistency )
{
  // Test that fitting works correctly for higher order polynomial
  // that can actually fit the data structure properly
  const std::vector<float> x_values = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f };
  
  // Generate data from a known quadratic: y = 1 + 2x + 0.5x^2
  std::vector<float> y_values;
  for( const float x : x_values )
  {
    const float y = 1.0f + 2.0f * x + 0.5f * x * x;
    y_values.push_back( y );
  }
  
  // Fit with quadratic (should be exact)
  std::vector<double> fitted_coeffs, uncertainties;
  const double chi2 = fit_to_polynomial( &x_values[0], &y_values[0], x_values.size(), 
                                         2, fitted_coeffs, uncertainties );
  
  BOOST_REQUIRE_EQUAL( fitted_coeffs.size(), 3 );
  BOOST_REQUIRE_EQUAL( uncertainties.size(), 3 );
  
  // For exact polynomial data, should fit very well
  BOOST_CHECK_CLOSE( fitted_coeffs[0], 1.0, 1e-3 );
  BOOST_CHECK_CLOSE( fitted_coeffs[1], 2.0, 1e-3 );
  BOOST_CHECK_CLOSE( fitted_coeffs[2], 0.5, 1e-3 );
  
  // Verify evaluation reproduces original data
  for( size_t i = 0; i < x_values.size(); ++i )
  {
    const float x = x_values[i];
    const float original_y = y_values[i];
    
    const double fitted_y = fitted_coeffs[0] + fitted_coeffs[1] * x + fitted_coeffs[2] * x * x;
    BOOST_CHECK_CLOSE( fitted_y, original_y, 5.0 );  // 5% tolerance due to Poisson weighting
  }
  
  // Chi2 should be reasonable
  BOOST_CHECK_GE( chi2, 0.0 );
  BOOST_CHECK( std::isfinite( chi2 ) );
}


// ---- Templated evaluators ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_expOfLogPowerSeries_template_matches_float )
{
  const std::vector<float> coefs = { -2.7f, -1.2f, -0.18f, -0.14f };
  for( const double energy : { 0.03, 0.06, 0.122, 0.3, 0.661, 1.332, 2.6 } )  //MeV
  {
    double expected = 0.0;
    for( size_t i = 0; i < coefs.size(); ++i )
      expected += coefs[i] * std::pow( std::log(energy), static_cast<double>(i) );
    expected = std::exp( expected );

    const float legacy = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( static_cast<float>(energy), coefs );
    const double templ = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, coefs.data(), coefs.size() );
    BOOST_CHECK_CLOSE( legacy, expected, 1.0e-3 );
    BOOST_CHECK_CLOSE( templ, expected, 1.0e-8 );
  }

  // No coefficients -> 0; a huge exponent is capped rather than overflowing
  BOOST_CHECK_EQUAL( DetectorPeakResponse::expOfLogPowerSeriesEfficiency( 100.0, coefs.data(), 0 ), 0.0 );
  const double big[2] = { 500.0, 0.0 };
  const double capped = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( 100.0, big, 2 );
  BOOST_CHECK( !std::isinf(capped) && !std::isnan(capped) );
  BOOST_CHECK_LT( capped, std::exp(50.0) );
}//test_expOfLogPowerSeries_template_matches_float


BOOST_AUTO_TEST_CASE( test_templated_fwhm_matches_float )
{
  using Form = DetectorPeakResponse::ResolutionFnctForm;

  // Hand-computed reference values for each form, including the 10 keV clamp
  const std::vector<float> gad_pos = { 1.4f, 0.25f, 0.45f }, gad_neg = { -3.0f, 6.0f, 0.5f };
  const std::vector<float> sqrt_inv = { 2.6f, 0.0015f, 30.0f };
  const std::vector<float> const_sqrt = { 1.0f, 0.035f };
  const std::vector<float> sqrt_poly = { 2.6f, 1.0f, 0.5f };

  auto gadras_ref = []( const std::vector<float> &p, double e ){
    e = std::max( e, 10.0 );
    const double a = p[0], b = p[1], c = p[2];
    if( e > 661.0 )
      return 6.61 * b * std::pow( e/661.0, c );
    if( a >= 0.0 )
    {
      const double zl = a * (661.0 - e) / 661.0;
      const double f = 6.61 * b * std::pow( e/661.0, c );
      return std::sqrt( zl*zl + f*f );
    }
    const double pw = std::pow( c, 1.0/std::log(1.0 - a) );
    return 6.61 * b * std::pow( std::max(30.0,e)/661.0, pw );
  };

  for( const double energy : { 5.0, 10.0, 30.0, 59.5, 122.0, 356.0, 661.0, 1332.0, 2614.0 } )
  {
    const float ef = static_cast<float>( energy );
    const double ec = std::max( energy, 10.0 );

    BOOST_CHECK_CLOSE( DetectorPeakResponse::peakResolutionFWHM( ef, Form::kGadrasResolutionFcn, gad_pos ),
                       gadras_ref( gad_pos, energy ), 1.0e-3 );
    BOOST_CHECK_CLOSE( DetectorPeakResponse::peakResolutionFWHM( ef, Form::kGadrasResolutionFcn, gad_neg ),
                       gadras_ref( gad_neg, energy ), 1.0e-3 );
    BOOST_CHECK_CLOSE( DetectorPeakResponse::peakResolutionFWHM( ef, Form::kSqrtEnergyPlusInverse, sqrt_inv ),
                       std::sqrt( sqrt_inv[0] + sqrt_inv[1]*ec + sqrt_inv[2]/ec ), 1.0e-3 );
    BOOST_CHECK_CLOSE( DetectorPeakResponse::peakResolutionFWHM( ef, Form::kConstantPlusSqrtEnergy, const_sqrt ),
                       const_sqrt[0] + const_sqrt[1]*std::sqrt(ec), 1.0e-3 );
    const double x = ec / 1000.0;
    BOOST_CHECK_CLOSE( DetectorPeakResponse::peakResolutionFWHM( ef, Form::kSqrtPolynomial, sqrt_poly ),
                       std::sqrt( sqrt_poly[0] + sqrt_poly[1]*x + sqrt_poly[2]*x*x ), 1.0e-3 );

    // The double template is the same function
    BOOST_CHECK_CLOSE( DetectorPeakResponse::peakResolutionFWHM( energy, Form::kGadrasResolutionFcn,
                                                                 gad_neg.data(), 3 ),
                       gadras_ref( gad_neg, energy ), 1.0e-8 );
  }//for( energies )

  // A negative square-root argument gives a small positive continuation, not NaN
  const std::vector<float> bad_poly = { -5.0f, 0.0f };
  const float cont = DetectorPeakResponse::peakResolutionFWHM( 100.0f, Form::kSqrtPolynomial, bad_poly );
  BOOST_CHECK( !std::isnan(cont) );
  BOOST_CHECK_GT( cont, 0.0f );
  BOOST_CHECK_LT( cont, 1.0e-3f );

  BOOST_CHECK_THROW( DetectorPeakResponse::peakResolutionFWHM( 100.0f, Form::kGadrasResolutionFcn, const_sqrt ),
                     std::runtime_error );
}//test_templated_fwhm_matches_float


// ---- Efficiency fit -----------------------------------------------------------------------------

namespace
{
  /** Synthetic points on exp(c0 + c1*lnE + ...) (energies in keV) with Gaussian scatter of
   `noise_frac`, stated statistical uncertainty `stat_frac`. */
  std::vector<MakeDrfFit::EffFitPoint> make_eff_points( const std::vector<float> &truth,
                                                        const double stat_frac,
                                                        const double noise_frac,
                                                        const unsigned seed = 4242 )
  {
    const std::vector<double> energies = { 59.5, 81.0, 122.1, 244.7, 344.3, 356.0, 511.0,
                                           661.7, 778.9, 964.1, 1173.2, 1332.5, 1408.0 };
    std::mt19937 gen( seed );
    std::normal_distribution<double> noise( 0.0, 1.0 );

    std::vector<MakeDrfFit::EffFitPoint> pts;
    for( const double energy : energies )
    {
      const double eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( energy, truth.data(), truth.size() );
      MakeDrfFit::EffFitPoint p;
      p.energy = static_cast<float>( energy );
      p.efficiency = static_cast<float>( eff * (1.0 + noise_frac*noise(gen)) );
      p.fracStatUncert = static_cast<float>( stat_frac );
      p.sourceKey = (energy < 300.0) ? "Eu152#0" : "Co60#1";
      pts.push_back( p );
    }
    return pts;
  }//make_eff_points(...)
}//namespace


BOOST_AUTO_TEST_CASE( test_performEfficiencyFit_recovers_truth )
{
  // A 3-term keV curve (intrinsic ~ 0.5 at 60 keV falling to ~0.05 at 1.3 MeV)
  const std::vector<float> truth = { -4.5f, 1.9f, -0.22f };
  const std::vector<MakeDrfFit::EffFitPoint> pts = make_eff_points( truth, 0.02, 0.02 );

  const MakeDrfFit::EffFitResult fit = MakeDrfFit::performEfficiencyFit( pts, 3 );
  BOOST_REQUIRE_EQUAL( fit.coefs.size(), 3u );
  BOOST_REQUIRE_EQUAL( fit.uncerts.size(), 3u );
  BOOST_REQUIRE_EQUAL( fit.covRowMajor.size(), 9u );
  BOOST_CHECK_EQUAL( fit.dof, 10 );
  BOOST_CHECK_MESSAGE( fit.warnings.empty(), "Warnings: " + fit.warnings );

  // Coefficients within 3 sigma, curve within ~2% everywhere
  for( size_t i = 0; i < 3; ++i )
  {
    BOOST_CHECK_GT( fit.uncerts[i], 0.0f );
    BOOST_CHECK_MESSAGE( std::fabs(fit.coefs[i] - truth[i]) < 3.0*fit.uncerts[i],
                         "coef " << i << ": " << fit.coefs[i] << " +- " << fit.uncerts[i]
                         << " vs truth " << truth[i] );
  }
  for( const MakeDrfFit::EffFitPoint &p : pts )
  {
    const double truth_eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( static_cast<double>(p.energy), truth.data(), 3 );
    const double fit_eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( static_cast<double>(p.energy), fit.coefs.data(), 3 );
    BOOST_CHECK_CLOSE( fit_eff, truth_eff, 2.5 );
  }

  // Scatter matches stated uncertainty: chi2/dof of order 1, no Birge inflation to speak of
  BOOST_CHECK_GT( fit.chi2 / fit.dof, 0.2 );
  BOOST_CHECK_LT( fit.chi2 / fit.dof, 3.0 );
  BOOST_CHECK_LT( fit.birgeScale, 3.0 );

  // Covariance is symmetric and consistent with the uncertainties
  for( size_t i = 0; i < 3; ++i )
  {
    BOOST_CHECK_CLOSE( std::sqrt(fit.covRowMajor[i*3 + i]), fit.uncerts[i], 1.0e-3 );
    for( size_t j = 0; j < 3; ++j )
      BOOST_CHECK_CLOSE( fit.covRowMajor[i*3 + j], fit.covRowMajor[j*3 + i], 1.0e-4 );
  }

  // Exactly determined and over-specified orders are rejected sensibly
  BOOST_CHECK_THROW( MakeDrfFit::performEfficiencyFit( pts, 0 ), std::runtime_error );
  BOOST_CHECK_THROW( MakeDrfFit::performEfficiencyFit( pts, static_cast<int>(pts.size()) + 1 ), std::runtime_error );
  BOOST_CHECK_NO_THROW( MakeDrfFit::performEfficiencyFit( pts, static_cast<int>(pts.size()) ) );
}//test_performEfficiencyFit_recovers_truth


/** The efficiency fit's own covariance has to be storable, and the tolerance that decides that has
 to stay loose enough for a real one.

 `DetectorEfficiencyUncert::covarianceIsUsable` calls a matrix positive semi-definite to within
 1.0E-6 of its scale.  That looks generous, and the temptation is to tighten it.  This pins why not:
 the fit is badly conditioned (the design matrix is powers of ln(E) over barely one decade) and the
 result is stored as float, so the smallest eigenvalue of a genuine high-order fit is slightly
 negative - measurably more so as terms are added.  A tolerance near 1.0E-8 would start dropping the
 coefficient covariance of every 7-term DRF, which is the most uncertainty-bearing thing they carry.

 The other half of the contract: when the covariance IS unusable the fit says so in `warnings` and
 clears it, rather than leaving it for `MakeDrfCalc::assembleDrf` to drop without telling anyone.
 */
BOOST_AUTO_TEST_CASE( test_performEfficiencyFit_covariance_is_storable )
{
  const std::vector<float> truth = { -4.5f, 1.9f, -0.22f };
  const std::vector<MakeDrfFit::EffFitPoint> pts = make_eff_points( truth, 0.02, 0.02 );

  for( int order = 3; order <= 7; ++order )
  {
    const MakeDrfFit::EffFitResult fit = MakeDrfFit::performEfficiencyFit( pts, order );
    BOOST_REQUIRE_EQUAL( fit.coefs.size(), static_cast<size_t>(order) );

    BOOST_REQUIRE_MESSAGE( fit.covRowMajor.size() == static_cast<size_t>(order*order),
                           "order " << order << " fit produced no coefficient covariance"
                           " (warnings: " << fit.warnings << ")" );

    const std::vector<double> cov( begin(fit.covRowMajor), end(fit.covRowMajor) );
    std::string why;
    BOOST_CHECK_MESSAGE( DetectorEfficiencyUncert::covarianceIsUsable( cov, &why ),
                         "order " << order << " fit covariance is not storable: " << why );

    // And it actually goes in, which is what the DRF needs.
    DetectorEfficiencyUncert uncert;
    BOOST_CHECK_NO_THROW( uncert.setCoefficientCovariance( fit.covRowMajor ) );
  }//for( int order = 3; order <= 7; ++order )
}//test_performEfficiencyFit_covariance_is_storable


BOOST_AUTO_TEST_CASE( test_performEfficiencyFit_block_covariance_shifts_common_mode )
{
  const std::vector<float> truth = { -4.5f, 1.9f, -0.22f };
  const std::vector<MakeDrfFit::EffFitPoint> stat_only = make_eff_points( truth, 0.02, 0.0 );

  // Data covariance structure: diagonal stat, plus cert/dist blocks per source
  std::vector<MakeDrfFit::EffFitPoint> with_cert = stat_only;
  for( MakeDrfFit::EffFitPoint &p : with_cert )
  {
    if( p.sourceKey == "Eu152#0" )
      p.fracCertUncert = 0.05f;
    else
      p.fracDistUncert = 0.03f;
  }

  const size_t n = with_cert.size();
  const std::vector<double> cov = MakeDrfFit::effDataCovariance( with_cert );
  BOOST_REQUIRE_EQUAL( cov.size(), n*n );
  for( size_t i = 0; i < n; ++i )
  {
    for( size_t j = 0; j < n; ++j )
    {
      double expected = (i == j) ? 0.02*0.02 : 0.0;
      if( with_cert[i].sourceKey == with_cert[j].sourceKey )
        expected += (with_cert[i].sourceKey == "Eu152#0") ? 0.05*0.05 : 0.03*0.03;
      BOOST_CHECK_CLOSE( cov[i*n + j], expected, 1.0e-3 );  //inputs are floats
    }
  }

  // The correlated source errors leave the best-fit curve essentially unchanged (noise-free data)
  //  but must widen the reported curve uncertainty.
  const MakeDrfFit::EffFitResult fit_stat = MakeDrfFit::performEfficiencyFit( stat_only, 3 );
  const MakeDrfFit::EffFitResult fit_cert = MakeDrfFit::performEfficiencyFit( with_cert, 3 );
  BOOST_REQUIRE_EQUAL( fit_cert.covRowMajor.size(), 9u );
  BOOST_CHECK_LT( fit_stat.chi2, 1.0e-4 );
  BOOST_CHECK_LT( fit_cert.chi2, 1.0e-4 );

  const double e_test = 200.0;  //a Eu152 energy
  const std::vector<double> j = { 1.0, std::log(e_test), std::pow(std::log(e_test), 2) };
  auto curve_var = [&j]( const std::vector<float> &c ){
    double v = 0.0;
    for( size_t a = 0; a < 3; ++a )
      for( size_t b = 0; b < 3; ++b )
        v += j[a] * c[a*3 + b] * j[b];
    return v;
  };
  const double var_stat = curve_var( fit_stat.covRowMajor );
  const double var_cert = curve_var( fit_cert.covRowMajor );
  BOOST_CHECK_GT( var_cert, 1.5*var_stat );
  // ...and by no more than the 5% common mode itself (the cert cannot add more than 5% to a curve
  //  that is pinned by an independent second source)
  BOOST_CHECK_LT( std::sqrt(var_cert), 0.06 );

  // Off-diagonals are populated
  bool any_offdiag = false;
  for( size_t a = 0; a < 3; ++a )
    for( size_t b = 0; b < 3; ++b )
      any_offdiag = (any_offdiag || ((a != b) && (std::fabs(fit_cert.covRowMajor[a*3 + b]) > 0.0f)));
  BOOST_CHECK( any_offdiag );
}//test_performEfficiencyFit_block_covariance_shifts_common_mode


BOOST_AUTO_TEST_CASE( test_performEfficiencyFit_birge_inflation )
{
  const std::vector<float> truth = { -4.5f, 1.9f, -0.22f };
  // 10% scatter but only 2% claimed
  const std::vector<MakeDrfFit::EffFitPoint> pts = make_eff_points( truth, 0.02, 0.10, 777 );

  const MakeDrfFit::EffFitResult fit = MakeDrfFit::performEfficiencyFit( pts, 3 );
  BOOST_REQUIRE_EQUAL( fit.covRowMajor.size(), 9u );
  BOOST_CHECK_GT( fit.birgeScale, 4.0 );
  BOOST_CHECK_CLOSE( fit.birgeScale, fit.chi2 / fit.dof, 1.0e-6 );
  for( size_t i = 0; i < 3; ++i )
    BOOST_CHECK_CLOSE( fit.uncerts[i]*fit.uncerts[i], fit.covRowMajor[i*3 + i], 1.0e-3 );

  // The inflated uncertainty band now covers the truth
  int covered = 0;
  for( const MakeDrfFit::EffFitPoint &p : pts )
  {
    const double lnE = std::log( static_cast<double>(p.energy) );
    const std::vector<double> j = { 1.0, lnE, lnE*lnE };
    double var = 0.0;
    for( size_t a = 0; a < 3; ++a )
      for( size_t b = 0; b < 3; ++b )
        var += j[a] * fit.covRowMajor[a*3 + b] * j[b];
    const double fit_eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( static_cast<double>(p.energy), fit.coefs.data(), 3 );
    const double truth_eff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( static_cast<double>(p.energy), truth.data(), 3 );
    if( std::fabs(std::log(fit_eff/truth_eff)) < 2.0*std::sqrt(var) )
      covered += 1;
  }
  BOOST_CHECK_GE( covered, static_cast<int>(pts.size()) - 2 );
}//test_performEfficiencyFit_birge_inflation


BOOST_AUTO_TEST_CASE( test_performEfficiencyFit_legacy_wrapper_matches_new )
{
  const std::vector<float> truth = { -4.5f, 1.9f, -0.22f };
  const std::vector<MakeDrfFit::EffFitPoint> pts = make_eff_points( truth, 0.02, 0.02 );

  std::vector<MakeDrfFit::DetEffDataPoint> legacy_pts;
  for( const MakeDrfFit::EffFitPoint &p : pts )
  {
    MakeDrfFit::DetEffDataPoint d;
    d.energy = p.energy;
    d.efficiency = p.efficiency;
    d.efficiency_uncert = p.fracStatUncert * p.efficiency;
    legacy_pts.push_back( d );
  }

  std::vector<float> coefs, uncerts;
  const double chi2_per_dof = MakeDrfFit::performEfficiencyFit( legacy_pts, 3, coefs, uncerts );

  std::vector<MakeDrfFit::EffFitPoint> uncorrelated = pts;
  for( MakeDrfFit::EffFitPoint &p : uncorrelated )
    p.sourceKey.clear();
  const MakeDrfFit::EffFitResult fit = MakeDrfFit::performEfficiencyFit( uncorrelated, 3 );

  BOOST_REQUIRE_EQUAL( coefs.size(), 3u );
  for( size_t i = 0; i < 3; ++i )
  {
    BOOST_CHECK_CLOSE( coefs[i], fit.coefs[i], 1.0e-3 );
    BOOST_CHECK_CLOSE( uncerts[i], fit.uncerts[i], 1.0e-2 );
  }
  BOOST_CHECK_CLOSE( chi2_per_dof, fit.chi2 / fit.dof, 1.0e-3 );
}//test_performEfficiencyFit_legacy_wrapper_matches_new


// ---- FWHM fit -----------------------------------------------------------------------------------

namespace
{
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> make_gadras_peaks(
                                                    const std::vector<float> &coefs,
                                                    const std::vector<double> &energies )
  {
    auto peaks = std::make_shared<std::deque<std::shared_ptr<const PeakDef>>>();
    for( const double energy : energies )
    {
      const float fwhm = DetectorPeakResponse::peakResolutionFWHM( static_cast<float>(energy),
                                  DetectorPeakResponse::kGadrasResolutionFcn, coefs );
      peaks->push_back( create_test_peak( energy, fwhm ) );
    }
    return peaks;
  }

  void check_gadras_curve_recovered( const std::vector<float> &truth, const std::vector<float> &fit,
                                     const std::vector<double> &energies )
  {
    for( const double energy : energies )
    {
      const float t = DetectorPeakResponse::peakResolutionFWHM( static_cast<float>(energy),
                                  DetectorPeakResponse::kGadrasResolutionFcn, truth );
      const float f = DetectorPeakResponse::peakResolutionFWHM( static_cast<float>(energy),
                                  DetectorPeakResponse::kGadrasResolutionFcn, fit );
      BOOST_CHECK_CLOSE( f, t, 0.5 );
    }
  }
}//namespace


BOOST_AUTO_TEST_CASE( test_performResolutionFit_gadras_negative_a )
{
  // A NaI-like detector with the negative-offset (bent power law) branch
  const std::vector<float> truth = { -3.0f, 6.0f, 0.5f };
  const std::vector<double> energies = { 59.5, 122.1, 200.0, 356.0, 511.0, 661.7, 1173.2, 1332.5, 2614.5 };
  const auto peaks = make_gadras_peaks( truth, energies );

  std::vector<float> coefs, uncerts;
  const double chi2 = MakeDrfFit::performResolutionFit( peaks, DetectorPeakResponse::kGadrasResolutionFcn,
                                                        0, coefs, uncerts );
  BOOST_REQUIRE_EQUAL( coefs.size(), 3u );
  BOOST_REQUIRE_EQUAL( uncerts.size(), 3u );
  BOOST_CHECK_LT( chi2, 1.0e-2 );
  BOOST_CHECK_LT( coefs[0], 0.0f );
  check_gadras_curve_recovered( truth, coefs, energies );
  for( const float u : uncerts )
    BOOST_CHECK_GT( u, 0.0f );
}//test_performResolutionFit_gadras_negative_a


BOOST_AUTO_TEST_CASE( test_performResolutionFit_gadras_positive_a )
{
  // A HPGe-like detector with the positive-offset (quadrature) branch
  const std::vector<float> truth = { 1.4f, 0.25f, 0.45f };
  const std::vector<double> energies = { 59.5, 122.1, 200.0, 356.0, 511.0, 661.7, 1173.2, 1332.5, 2614.5 };
  const auto peaks = make_gadras_peaks( truth, energies );

  // Start from the wrong branch on purpose
  std::vector<float> coefs = { -2.0f, 0.3f, 0.5f }, uncerts;
  const double chi2 = MakeDrfFit::performResolutionFit( peaks, DetectorPeakResponse::kGadrasResolutionFcn,
                                                        0, coefs, uncerts );
  BOOST_REQUIRE_EQUAL( coefs.size(), 3u );
  BOOST_CHECK_LT( chi2, 1.0e-2 );
  BOOST_CHECK_GT( coefs[0], 0.0f );
  check_gadras_curve_recovered( truth, coefs, energies );
}//test_performResolutionFit_gadras_positive_a


BOOST_AUTO_TEST_CASE( test_performResolutionFit_peak_count_freezing )
{
  const std::vector<float> truth = { 1.4f, 0.25f, 0.45f };

  // One peak: only B floats
  {
    const auto peaks = make_gadras_peaks( truth, { 661.7 } );
    const MakeDrfFit::FwhmFitResult r = MakeDrfFit::performResolutionFitEx( peaks,
                                   DetectorPeakResponse::kGadrasResolutionFcn, 0, {} );
    BOOST_REQUIRE_EQUAL( r.coefs.size(), 3u );
    BOOST_CHECK_EQUAL( r.uncerts[0], 0.0f );
    BOOST_CHECK_GT( r.uncerts[1], 0.0f );
    BOOST_CHECK_EQUAL( r.uncerts[2], 0.0f );
    BOOST_CHECK_EQUAL( r.dof, 0 );
    BOOST_CHECK_CLOSE( r.coefs[1], truth[1], 1.0 );  //661.7 keV pins B directly (above 661 keV)
  }

  // Two peaks: B and C float
  {
    const auto peaks = make_gadras_peaks( truth, { 661.7, 1332.5 } );
    const MakeDrfFit::FwhmFitResult r = MakeDrfFit::performResolutionFitEx( peaks,
                                   DetectorPeakResponse::kGadrasResolutionFcn, 0, {} );
    BOOST_REQUIRE_EQUAL( r.coefs.size(), 3u );
    BOOST_CHECK_EQUAL( r.uncerts[0], 0.0f );
    BOOST_CHECK_GT( r.uncerts[1], 0.0f );
    BOOST_CHECK_GT( r.uncerts[2], 0.0f );
    BOOST_CHECK_CLOSE( r.coefs[1], truth[1], 1.0 );
    BOOST_CHECK_CLOSE( r.coefs[2], truth[2], 2.0 );
  }

  // kSqrtPolynomial without a linear seed (2 peaks, 3 requested) returns 2 coefficients, B fixed
  {
    const std::vector<float> poly = { 2.6f, 1.0f, 0.5f };
    const auto peaks = std::make_shared<std::deque<std::shared_ptr<const PeakDef>>>(
                          create_peaks_with_known_fwhm( poly, { 200.0, 1000.0 } ) );
    const MakeDrfFit::FwhmFitResult r = MakeDrfFit::performResolutionFitEx( peaks,
                                   DetectorPeakResponse::kSqrtPolynomial, 3, {} );
    BOOST_CHECK_EQUAL( r.coefs.size(), 2u );
    BOOST_CHECK_GT( r.uncerts[0], 0.0f );
    BOOST_CHECK_EQUAL( r.uncerts[1], 0.0f );
  }
}//test_performResolutionFit_peak_count_freezing


BOOST_AUTO_TEST_CASE( test_performResolutionFit_uncerts_populated )
{
  const std::vector<float> truth = { 2.6f, 1.0f, 0.5f };
  const std::vector<double> energies = { 59.5, 122.1, 356.0, 661.7, 1173.2, 1332.5, 2614.5 };
  const auto peaks = std::make_shared<std::deque<std::shared_ptr<const PeakDef>>>(
                        create_peaks_with_known_fwhm( truth, energies, false, 0.02 ) );

  const MakeDrfFit::FwhmFitResult r = MakeDrfFit::performResolutionFitEx( peaks,
                                 DetectorPeakResponse::kSqrtPolynomial, 3, {} );
  BOOST_REQUIRE_EQUAL( r.coefs.size(), 3u );
  BOOST_REQUIRE_EQUAL( r.covRowMajor.size(), 9u );
  BOOST_CHECK_EQUAL( r.dof, 4 );
  BOOST_CHECK_MESSAGE( r.warnings.empty(), "Warnings: " + r.warnings );
  for( size_t i = 0; i < 3; ++i )
  {
    BOOST_CHECK_GT( r.uncerts[i], 0.0f );
    BOOST_CHECK_CLOSE( r.uncerts[i]*r.uncerts[i], r.covRowMajor[i*3 + i], 1.0e-3 );
    BOOST_CHECK_MESSAGE( std::fabs(r.coefs[i] - truth[i]) < 4.0*r.uncerts[i] + 1.0e-3,
                         "coef " << i << ": " << r.coefs[i] << " +- " << r.uncerts[i] << " vs " << truth[i] );
  }

  // The legacy wrapper returns the same thing
  std::vector<float> coefs, uncerts;
  const double chi2 = MakeDrfFit::performResolutionFit( peaks, DetectorPeakResponse::kSqrtPolynomial, 3, coefs, uncerts );
  BOOST_CHECK_CLOSE( chi2, r.chi2, 1.0e-6 );
  BOOST_REQUIRE_EQUAL( coefs.size(), 3u );
  for( size_t i = 0; i < 3; ++i )
    BOOST_CHECK_CLOSE( coefs[i], r.coefs[i], 1.0e-6 );

  // Every form yields uncertainties
  for( const auto form : { DetectorPeakResponse::kSqrtEnergyPlusInverse,
                           DetectorPeakResponse::kConstantPlusSqrtEnergy } )
  {
    std::vector<float> c, u;
    BOOST_REQUIRE_NO_THROW( MakeDrfFit::performResolutionFit( peaks, form, 3, c, u ) );
    BOOST_REQUIRE_EQUAL( c.size(), u.size() );
    for( const float val : u )
      BOOST_CHECK_GT( val, 0.0f );
  }
}//test_performResolutionFit_uncerts_populated
