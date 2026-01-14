#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

//#if __has_include("ceres/jet.h") && __has_include("eigen3/Eigen/Core")
#  define HAS_CERES_JET 1
#  include "eigen3/Eigen/Core"
#  include "ceres/jet.h"
//#else
//#  define HAS_CERES_JET 0
//#endif

#include "../pseudo_voigt_exp_tail.hpp"
#include "scipy_voigt_test_data.hpp"

namespace {

using namespace pseudo_voigt;

BOOST_AUTO_TEST_SUITE(PseudoVoigtExpTailTestSuite)

BOOST_AUTO_TEST_CASE(TestNormalization) {
    const double mean = 94.67;
    const double sigma = 0.5;
    const double gamma = 0.0876;
    const double tail_ratio = 0.10;
    const double tail_slope = 1.0;

    // Integrate from -infinity to +infinity using a very wide range to cover long tails
    double x_min = mean - 200.0 * sigma;
    double x_max = mean + 200.0 * sigma;

    double cdf_min = voigt_exp_indefinite(x_min, mean, sigma, gamma, tail_ratio, tail_slope);
    double cdf_max = voigt_exp_indefinite(x_max, mean, sigma, gamma, tail_ratio, tail_slope);
    double integral = cdf_max - cdf_min;

    // Note: Tolerance is 0.01 because pseudo-Voigt approximation may have slight normalization differences
    BOOST_CHECK_CLOSE(integral, 1.0, 1.0);
}

BOOST_AUTO_TEST_CASE(TestPureVoigt) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;  // Pure Gaussian
    const double tail_ratio = 0.0;
    const double tail_slope = 1.0;

    // Test at mean
    double pdf_at_mean = voigt_pdf(mean, mean, sigma, gamma);
    double expected_gaussian = 1.0 / (sigma * std::sqrt(2.0 * M_PI));

    BOOST_CHECK_CLOSE(pdf_at_mean, expected_gaussian, 1e-6);
}

BOOST_AUTO_TEST_CASE(TestPureGaussExp) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;
    const double tail_ratio = 1.0;
    const double tail_slope = 1.0;

    // Test normalization
    double x_min = mean - 20.0 * sigma;
    double x_max = mean + 20.0 * sigma;
    double cdf_min = voigt_exp_indefinite(x_min, mean, sigma, gamma, tail_ratio, tail_slope);
    double cdf_max = voigt_exp_indefinite(x_max, mean, sigma, gamma, tail_ratio, tail_slope);
    double integral = cdf_max - cdf_min;

    BOOST_CHECK(std::abs(integral - 1.0) < 1e-4); // tolerance for exponential tails
}

BOOST_AUTO_TEST_CASE(TestSymmetry) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.0;  // Pure Voigt, should be symmetric
    const double tail_slope = 1.0;

    double offset = 2.0;
    double pdf_left = voigt_pdf(mean - offset, mean, sigma, gamma);
    double pdf_right = voigt_pdf(mean + offset, mean, sigma, gamma);

    // Pseudo-Voigt should be symmetric (within approximation error)
    BOOST_CHECK_CLOSE(pdf_left, pdf_right, 1e-3);
}

BOOST_AUTO_TEST_CASE(TestCoverageLimits) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.1;
    const double tail_slope = 1.0;
    const double p = 0.05;  // 95% coverage

    auto limits = voigt_exp_coverage_limits(mean, sigma, gamma, tail_ratio, tail_slope, p);

    double cdf_lower = voigt_exp_indefinite(limits.first, mean, sigma, gamma, tail_ratio, tail_slope);
    double cdf_upper = voigt_exp_indefinite(limits.second, mean, sigma, gamma, tail_ratio, tail_slope);
    double coverage = cdf_upper - cdf_lower;

    BOOST_CHECK_CLOSE(coverage, (1.0 - p), 1e-2);
}

BOOST_AUTO_TEST_CASE(TestPureGaussian) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;  // No Lorentzian
    const double tail_ratio = 0.0;  // No tail
    const double tail_slope = 1.0;

    // Test PDF at mean
    double pdf_at_mean = voigt_pdf(mean, mean, sigma, gamma);
    double expected_gaussian = 1.0 / (sigma * std::sqrt(2.0 * M_PI));

    BOOST_CHECK_CLOSE(pdf_at_mean, expected_gaussian, 1e-6);

    // Test CDF
    double cdf_at_mean = voigt_cdf(mean, mean, sigma, gamma);
    double expected_cdf = 0.5;

    BOOST_CHECK_CLOSE(cdf_at_mean, expected_cdf, 1e-6);

    // Test symmetry
    double offset = 2.0;
    double pdf_left = voigt_pdf(mean - offset, mean, sigma, gamma);
    double pdf_right = voigt_pdf(mean + offset, mean, sigma, gamma);

    BOOST_CHECK_CLOSE(pdf_left, pdf_right, 1e-6);
}

BOOST_AUTO_TEST_CASE(TestGammaSpectroscopyParameters) {
    // Test with realistic gamma spectroscopy parameters
    const double mean = 94.67;      // Uranium 94.67 keV x-ray
    const double sigma = 0.5;       // 500 eV Gaussian width (typical HPGe)
    const double gamma = 0.0876;    // 87.6 eV Lorentzian width (realistic)
    const double tail_ratio = 0.05; // 5% tailing
    const double tail_slope = 1.0;  // Tail slope

    // Test PDF values at peak
    double pdf_at_peak = voigt_pdf(mean, mean, sigma, gamma);
    BOOST_CHECK(pdf_at_peak > 0.0);

    // Test CDF increases monotonically
    double cdf1 = voigt_exp_indefinite(mean - 1.0, mean, sigma, gamma, tail_ratio, tail_slope);
    double cdf2 = voigt_exp_indefinite(mean, mean, sigma, gamma, tail_ratio, tail_slope);
    double cdf3 = voigt_exp_indefinite(mean + 1.0, mean, sigma, gamma, tail_ratio, tail_slope);

    BOOST_CHECK(cdf1 < cdf2);
    BOOST_CHECK(cdf2 < cdf3);
    BOOST_CHECK(cdf1 >= 0.0);
    BOOST_CHECK(cdf3 <= 1.0);
}

BOOST_AUTO_TEST_CASE(TestEdgeCases) {
    const double mean = 100.0;
    const double sigma = 1.0;

    // Test with very small gamma (should behave like Gaussian)
    double gamma_small = 1e-10;
    double pdf_gauss = voigt_pdf(mean, mean, sigma, 0.0);
    double pdf_small_gamma = voigt_pdf(mean, mean, sigma, gamma_small);
    BOOST_CHECK_CLOSE(pdf_gauss, pdf_small_gamma, 1e-6);

    // Test with very small sigma (should behave like Lorentzian)
    double sigma_small = 1e-10;
    double gamma_test = 0.1;
    double x_test = mean + 1.0;

    // For Lorentzian, PDF should be gamma/(pi * (x-mean)^2 + gamma^2)
    double expected_lorentz = gamma_test / (M_PI * (gamma_test * gamma_test));
    double pdf_small_sigma = voigt_pdf(mean, mean, sigma_small, gamma_test);
    BOOST_CHECK_CLOSE(pdf_small_sigma, expected_lorentz, 1e-3);

    // Test with zero tail ratio (should be pure Voigt)
    double tail_ratio_zero = 0.0;
    double pdf_voigt = voigt_pdf(mean, mean, sigma, gamma_test);
    double cdf_voigt_with_tail = voigt_exp_indefinite(mean, mean, sigma, gamma_test, tail_ratio_zero, 1.0);
    double cdf_pure_voigt = voigt_cdf(mean, mean, sigma, gamma_test);
    BOOST_CHECK_CLOSE(cdf_voigt_with_tail, cdf_pure_voigt, 1e-6);
}

BOOST_AUTO_TEST_CASE(TestArrayIntegration) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double amplitude = 1234.5;
    const double gamma = 0.1;
    const double tail_ratio = 0.1;
    const double tail_slope = 2.0;

    const size_t num_channels = 100;
    const double x_min = mean - 5*sigma;
    const double x_max = mean + 5*sigma;

    std::vector<float> energies(num_channels + 1);
    for (size_t i = 0; i <= num_channels; ++i) {
        energies[i] = static_cast<float>(x_min + i * (x_max - x_min) / num_channels);
    }

    std::vector<double> counts(num_channels, 0.0);
    voigt_exp_integral(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                      energies.data(), counts.data(), num_channels);

    // Check that total counts equals amplitude (within tolerance)
    // Pseudo-Voigt approximation may have slight normalization differences
    double total_counts = std::accumulate(counts.begin(), counts.end(), 0.0);
    double error_fraction = std::abs(total_counts - amplitude) / amplitude;

    BOOST_CHECK_MESSAGE(error_fraction < 0.05,  // 5% tolerance for pseudo-Voigt
        "Array integration should give correct amplitude. "
        "Expected: " << amplitude << ", Actual: " << total_counts <<
        ", Error: " << (error_fraction * 100.0) << "%");
}

BOOST_AUTO_TEST_CASE(TestExtremeParameters) {
    // Test case for extreme parameters (gamma >> sigma)
    const double mean = 14.0;
    const double sigma = 0.0005;
    const double gamma = 0.006;  // gamma/sigma = 12 (Lorentzian-dominated)
    const double tail_ratio = 0.15;
    const double tail_slope = 1.0;

    // Test 1: voigt_cdf should return values in [0, 1] for all x
    double x_near_mean = 13.9970874786377;  // Slightly below mean
    double cdf_near_mean = voigt_cdf(x_near_mean, mean, sigma, gamma);

    BOOST_CHECK_MESSAGE(cdf_near_mean >= 0.0 && cdf_near_mean <= 1.0,
        "voigt_cdf returns " << cdf_near_mean << " for gamma >> sigma case. "
        "Expected [0,1].");

    // Test 2: voigt_exp_indefinite should also return values in [0, 1]
    double cdf_combined = voigt_exp_indefinite(x_near_mean, mean, sigma, gamma, tail_ratio, tail_slope);
    BOOST_CHECK_MESSAGE(cdf_combined >= 0.0 && cdf_combined <= 1.0,
        "voigt_exp_indefinite returns " << cdf_combined <<
        " for extreme parameters. Expected [0,1].");

    // Test 3: Integration over wide range should give unit area
    std::pair<double,double> limits = voigt_exp_coverage_limits(mean, sigma, gamma, tail_ratio, tail_slope, 1.0e-9);
    double x0 = limits.first;
    double x1 = limits.second;

    double integral = voigt_exp_integral(mean, sigma, gamma, tail_ratio, tail_slope, x0, x1);

    BOOST_CHECK_MESSAGE(std::abs(integral - 1.0) < 0.01,
        "voigt_exp_integral over coverage limits returns " << integral <<
        " instead of ~1.0 for extreme parameters (gamma=" << gamma << ", sigma=" << sigma <<
        ", tail_ratio=" << tail_ratio << ", tail_slope=" << tail_slope << ").");
}

BOOST_AUTO_TEST_CASE(TestInterSpecFailures) {
    // Same extreme parameters as TestExtremeParameters
    const double mean = 14.0;
    const double sigma = 0.0005;
    const double gamma = 0.006;  // gamma/sigma = 12 (Lorentzian-dominated)
    const double tail_ratio = 0.15;
    const double tail_slope = 1.0;
    const double amplitude = 1.23456;

    // Test 1: Array integration (voigt_exp_integral with energy bins) should give correct amplitude
    {
        std::pair<double,double> limits = voigt_exp_coverage_limits(mean, sigma, gamma, tail_ratio, tail_slope, 1.0e-9);
        double x0 = limits.first;
        double x1 = limits.second;

        size_t num_channels = 2048;
        std::vector<double> counts(num_channels, 0.0);
        std::vector<float> energies(num_channels + 1, 0.0f);
        for (size_t i = 0; i < energies.size(); ++i)
            energies[i] = static_cast<float>(x0 + (i * (x1 - x0) / num_channels));

        voigt_exp_integral(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                          energies.data(), counts.data(), num_channels);

        double answer_sum = std::accumulate(counts.begin(), counts.end(), 0.0);
        double error_fraction = std::abs(answer_sum - amplitude) / amplitude;

        BOOST_CHECK_MESSAGE(error_fraction < 1.0e-2,
            "Array integration should give correct amplitude. "
            "Expected: " << amplitude << ", Actual: " << answer_sum <<
            ", Error: " << (error_fraction * 100.0) << "%.");
    }

    // Test 2: Unit area test (integral over coverage limits should equal 1.0)
    {
        std::pair<double,double> limits = voigt_exp_coverage_limits(mean, sigma, gamma, tail_ratio, tail_slope, 1.0e-9);
        double x0 = limits.first;
        double x1 = limits.second;

        double total_area = voigt_exp_integral(mean, sigma, gamma, tail_ratio, tail_slope, x0, x1);
        double error_fraction = std::abs(total_area - 1.0);

        BOOST_CHECK_MESSAGE(error_fraction < 1.0e-2,
            "Integral over coverage limits should be 1.0. "
            "Expected: 1.0, Actual: " << total_area <<
            ", Error: " << (error_fraction * 100.0) << "%.");
    }

    // Test 3: voigt_exp_coverage_limits should return correct fraction outside limits
    {
        double prob = 1.0e-6;
        std::pair<double,double> limits = voigt_exp_coverage_limits(mean, sigma, gamma,
                                                                    tail_ratio, tail_slope, prob);
        double fraction_inside = voigt_exp_integral(mean, sigma, gamma, tail_ratio,
                                                   tail_slope, limits.first, limits.second);
        double fraction_outside = 1.0 - fraction_inside;
        double error_ratio = fraction_outside / prob;

        BOOST_CHECK_MESSAGE(std::abs(fraction_outside - prob) / prob < 0.1,
            "voigt_exp_coverage_limits should return correct fraction outside. "
            "Expected fraction outside: " << prob << ", Actual: " << fraction_outside <<
            ", Ratio: " << error_ratio << "× (should be ~1.0).");
    }
}

// ============================================================================
// Tests against scipy-generated Voigt PDF values
// Note: Pseudo-Voigt is an approximation, so tolerances are more relaxed
// ============================================================================

BOOST_AUTO_TEST_CASE(TestAgainstScipy_GaussianDominated) {
    using namespace ScipyVoigtTestData;
    using namespace pseudo_voigt;
    
    // Get 95% coverage limits (p=0.05 means 95% coverage)
    auto limits = voigt_exp_coverage_limits(Gaussian_dominated_mean, 
                                          Gaussian_dominated_sigma, 
                                          Gaussian_dominated_gamma,
                                          0.0, 1.0, 0.05);
    double x_min = limits.first;
    double x_max = limits.second;
    
    size_t points_tested = 0;
    size_t failures = 0;
    for (size_t i = 0; i < Gaussian_dominated_npoints; ++i) {
        double x = Gaussian_dominated_x_values[i];
        
        // Only test points within 95% coverage range
        if (x < x_min || x > x_max) continue;
        
        points_tested++;
        double scipy_pdf = Gaussian_dominated_pdf_values[i];
        double computed_pdf = voigt_pdf(x, Gaussian_dominated_mean, 
                                        Gaussian_dominated_sigma, 
                                        Gaussian_dominated_gamma);
        
        // Require 5% relative error or small absolute error
        // For very small PDF values, use absolute error threshold
        double abs_error = std::abs(computed_pdf - scipy_pdf);
        double rel_error = (scipy_pdf > 1e-10) ? abs_error / scipy_pdf : 0.0;
        // Use slightly more lenient absolute error for very small values
        double abs_threshold = (scipy_pdf < 1e-3) ? 1e-6 : 1e-5;
        if (rel_error > 0.05 && abs_error > abs_threshold) {
            failures++;
        }
    }
    
    // Require 100% of points within coverage range to pass
    BOOST_CHECK_MESSAGE(failures == 0,
        "PDF mismatch: " << failures << " out of " << points_tested 
        << " points failed (within 95% coverage range [" << x_min << ", " << x_max << "])");
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_LorentzianDominated) {
    using namespace ScipyVoigtTestData;
    using namespace pseudo_voigt;
    
    // Get 95% coverage limits (p=0.05 means 95% coverage)
    auto limits = voigt_exp_coverage_limits(Lorentzian_dominated_mean, 
                                          Lorentzian_dominated_sigma, 
                                          Lorentzian_dominated_gamma,
                                          0.0, 1.0, 0.05);
    double x_min = limits.first;
    double x_max = limits.second;
    
    size_t points_tested = 0;
    size_t failures = 0;
    for (size_t i = 0; i < Lorentzian_dominated_npoints; ++i) {
        double x = Lorentzian_dominated_x_values[i];
        
        // Only test points within 95% coverage range
        if (x < x_min || x > x_max) continue;
        
        points_tested++;
        double scipy_pdf = Lorentzian_dominated_pdf_values[i];
        double computed_pdf = voigt_pdf(x, Lorentzian_dominated_mean, 
                                        Lorentzian_dominated_sigma, 
                                        Lorentzian_dominated_gamma);
        
        // Require 5% relative error or small absolute error
        // For very small PDF values, use absolute error threshold
        double abs_error = std::abs(computed_pdf - scipy_pdf);
        double rel_error = (scipy_pdf > 1e-10) ? abs_error / scipy_pdf : 0.0;
        // Use slightly more lenient absolute error for very small values
        double abs_threshold = (scipy_pdf < 1e-3) ? 1e-6 : 1e-5;
        if (rel_error > 0.05 && abs_error > abs_threshold) {
            failures++;
        }
    }
    
    // Require 100% of points within coverage range to pass
    BOOST_CHECK_MESSAGE(failures == 0,
        "PDF mismatch: " << failures << " out of " << points_tested 
        << " points failed (within 95% coverage range [" << x_min << ", " << x_max << "])");
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_Balanced) {
    using namespace ScipyVoigtTestData;
    using namespace pseudo_voigt;
    
    // Get 95% coverage limits (p=0.05 means 95% coverage)
    auto limits = voigt_exp_coverage_limits(Balanced_mean, 
                                          Balanced_sigma, 
                                          Balanced_gamma,
                                          0.0, 1.0, 0.05);
    double x_min = limits.first;
    double x_max = limits.second;
    
    size_t points_tested = 0;
    size_t failures = 0;
    for (size_t i = 0; i < Balanced_npoints; ++i) {
        double x = Balanced_x_values[i];
        
        // Only test points within 95% coverage range
        if (x < x_min || x > x_max) continue;
        
        points_tested++;
        double scipy_pdf = Balanced_pdf_values[i];
        double computed_pdf = voigt_pdf(x, Balanced_mean, 
                                        Balanced_sigma, 
                                        Balanced_gamma);
        
        // Require 6% relative error for balanced case (pseudo-Voigt approximation)
        // For very small PDF values, use absolute error threshold
        double abs_error = std::abs(computed_pdf - scipy_pdf);
        double rel_error = (scipy_pdf > 1e-10) ? abs_error / scipy_pdf : 0.0;
        // Use slightly more lenient absolute error for very small values
        double abs_threshold = (scipy_pdf < 1e-3) ? 1e-6 : 1e-5;
        if (rel_error > 0.06 && abs_error > abs_threshold) {
            failures++;
        }
    }
    
    // Require 100% of points within coverage range to pass
    BOOST_CHECK_MESSAGE(failures == 0,
        "PDF mismatch: " << failures << " out of " << points_tested 
        << " points failed (within 95% coverage range [" << x_min << ", " << x_max << "])");
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_GammaSpectroscopy) {
    using namespace ScipyVoigtTestData;
    using namespace pseudo_voigt;
    
    // Get 95% coverage limits (p=0.05 means 95% coverage)
    auto limits = voigt_exp_coverage_limits(Gamma_spectroscopy_mean, 
                                          Gamma_spectroscopy_sigma, 
                                          Gamma_spectroscopy_gamma,
                                          0.0, 1.0, 0.05);
    double x_min = limits.first;
    double x_max = limits.second;
    
    size_t points_tested = 0;
    size_t failures = 0;
    for (size_t i = 0; i < Gamma_spectroscopy_npoints; ++i) {
        double x = Gamma_spectroscopy_x_values[i];
        
        // Only test points within 95% coverage range
        if (x < x_min || x > x_max) continue;
        
        points_tested++;
        double scipy_pdf = Gamma_spectroscopy_pdf_values[i];
        double computed_pdf = voigt_pdf(x, Gamma_spectroscopy_mean, 
                                        Gamma_spectroscopy_sigma, 
                                        Gamma_spectroscopy_gamma);
        
        // Require 5% relative error or small absolute error
        // For very small PDF values, use absolute error threshold
        double abs_error = std::abs(computed_pdf - scipy_pdf);
        double rel_error = (scipy_pdf > 1e-10) ? abs_error / scipy_pdf : 0.0;
        // Use slightly more lenient absolute error for very small values
        double abs_threshold = (scipy_pdf < 1e-3) ? 1e-6 : 1e-5;
        if (rel_error > 0.05 && abs_error > abs_threshold) {
            failures++;
        }
    }
    
    // Require 100% of points within coverage range to pass
    BOOST_CHECK_MESSAGE(failures == 0,
        "PDF mismatch: " << failures << " out of " << points_tested 
        << " points failed (within 95% coverage range [" << x_min << ", " << x_max << "])");
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_SmallSigmaLargeGamma) {
    using namespace ScipyVoigtTestData;
    using namespace pseudo_voigt;
    
    // Get 95% coverage limits (p=0.05 means 95% coverage)
    auto limits = voigt_exp_coverage_limits(Small_sigma_large_gamma_mean, 
                                          Small_sigma_large_gamma_sigma, 
                                          Small_sigma_large_gamma_gamma,
                                          0.0, 1.0, 0.05);
    double x_min = limits.first;
    double x_max = limits.second;
    
    size_t points_tested = 0;
    size_t failures = 0;
    for (size_t i = 0; i < Small_sigma_large_gamma_npoints; ++i) {
        double x = Small_sigma_large_gamma_x_values[i];
        
        // Only test points within 95% coverage range
        if (x < x_min || x > x_max) continue;
        
        points_tested++;
        double scipy_pdf = Small_sigma_large_gamma_pdf_values[i];
        double computed_pdf = voigt_pdf(x, Small_sigma_large_gamma_mean, 
                                        Small_sigma_large_gamma_sigma, 
                                        Small_sigma_large_gamma_gamma);
        
        // Require 5% relative error or small absolute error
        // For very small PDF values, use absolute error threshold
        double abs_error = std::abs(computed_pdf - scipy_pdf);
        double rel_error = (scipy_pdf > 1e-10) ? abs_error / scipy_pdf : 0.0;
        // Use slightly more lenient absolute error for very small values
        double abs_threshold = (scipy_pdf < 1e-3) ? 1e-6 : 1e-5;
        if (rel_error > 0.05 && abs_error > abs_threshold) {
            failures++;
        }
    }
    
    // Require 100% of points within coverage range to pass
    BOOST_CHECK_MESSAGE(failures == 0,
        "PDF mismatch: " << failures << " out of " << points_tested 
        << " points failed (within 95% coverage range [" << x_min << ", " << x_max << "])");
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_LargeSigmaSmallGamma) {
    using namespace ScipyVoigtTestData;
    using namespace pseudo_voigt;
    
    // Get 95% coverage limits (p=0.05 means 95% coverage)
    auto limits = voigt_exp_coverage_limits(Large_sigma_small_gamma_mean, 
                                          Large_sigma_small_gamma_sigma, 
                                          Large_sigma_small_gamma_gamma,
                                          0.0, 1.0, 0.05);
    double x_min = limits.first;
    double x_max = limits.second;
    
    size_t points_tested = 0;
    size_t failures = 0;
    for (size_t i = 0; i < Large_sigma_small_gamma_npoints; ++i) {
        double x = Large_sigma_small_gamma_x_values[i];
        
        // Only test points within 95% coverage range
        if (x < x_min || x > x_max) continue;
        
        points_tested++;
        double scipy_pdf = Large_sigma_small_gamma_pdf_values[i];
        double computed_pdf = voigt_pdf(x, Large_sigma_small_gamma_mean, 
                                        Large_sigma_small_gamma_sigma, 
                                        Large_sigma_small_gamma_gamma);
        
        // Require 5% relative error or small absolute error
        // For very small PDF values, use absolute error threshold
        double abs_error = std::abs(computed_pdf - scipy_pdf);
        double rel_error = (scipy_pdf > 1e-10) ? abs_error / scipy_pdf : 0.0;
        // Use slightly more lenient absolute error for very small values
        double abs_threshold = (scipy_pdf < 1e-3) ? 1e-6 : 1e-5;
        if (rel_error > 0.05 && abs_error > abs_threshold) {
            failures++;
        }
    }
    
    // Require 100% of points within coverage range to pass
    BOOST_CHECK_MESSAGE(failures == 0,
        "PDF mismatch: " << failures << " out of " << points_tested 
        << " points failed (within 95% coverage range [" << x_min << ", " << x_max << "])");
}

// ============================================================================
// Ceres Jet compatibility and automatic differentiation tests
// ============================================================================

#if HAS_CERES_JET
BOOST_AUTO_TEST_CASE(TestJetCompatibility) {
    std::cout << "\n[PseudoVoigt] Running TestJetCompatibility with Ceres Jet types...\n";
    using ceres::Jet;
    using T = Jet<double, 1>;

    // Realistic gamma spectroscopy parameters (HPGe around 95 keV)
    const double mean = 94.67;
    const double sigma = 0.5;
    const double gamma = 0.0876;
    const double tail_ratio = 0.05;
    const double tail_slope = 1.5;
    const double x = 95.2; // near peak shoulder

    const double h = 1e-5;
    const double abs_tol = 1e-7;
    const double rel_tol = 5e-5; // 50 ppm relative tolerance

    auto check_derivative = [&](auto fn, double param_base) {
        // Central finite difference reference
        double fwd = fn(param_base + h);
        double bwd = fn(param_base - h);
        double fd = (fwd - bwd) / (2 * h);

        // Jet-based derivative (v[0])
        T p;
        p.a = param_base;
        p.v[0] = 1.0;
        T val = fn(p);
        double jet_deriv = val.v[0];

        double mag = std::max(1.0, std::abs(fd));
        double tol = std::max(abs_tol, rel_tol * mag);
        BOOST_CHECK_SMALL(jet_deriv - fd, tol);
    };

    // Derivative w.r.t. mean for PDF
    check_derivative(
        [&](auto m) {
            return voigt_pdf(static_cast<decltype(m)>(x),
                             m,
                             static_cast<decltype(m)>(sigma),
                             static_cast<decltype(m)>(gamma));
        },
        mean);

    // Derivative w.r.t. sigma for CDF (indefinite integral with tail)
    check_derivative(
        [&](auto s) {
            return voigt_exp_indefinite(static_cast<decltype(s)>(x),
                                        static_cast<decltype(s)>(mean),
                                        s,
                                        static_cast<decltype(s)>(gamma),
                                        static_cast<decltype(s)>(tail_ratio),
                                        static_cast<decltype(s)>(tail_slope));
        },
        sigma);

    // Derivative w.r.t. gamma for PDF
    check_derivative(
        [&](auto g) {
            return voigt_pdf(static_cast<decltype(g)>(x),
                             static_cast<decltype(g)>(mean),
                             static_cast<decltype(g)>(sigma),
                             g);
        },
        gamma);

    // Derivative w.r.t. tail_ratio for the skewed CDF
    check_derivative(
        [&](auto t_ratio) {
            return voigt_exp_indefinite(static_cast<decltype(t_ratio)>(x),
                                        static_cast<decltype(t_ratio)>(mean),
                                        static_cast<decltype(t_ratio)>(sigma),
                                        static_cast<decltype(t_ratio)>(gamma),
                                        t_ratio,
                                        static_cast<decltype(t_ratio)>(tail_slope));
        },
        tail_ratio);

    // Derivative w.r.t. tail_slope for the skewed CDF
    check_derivative(
        [&](auto t_slope) {
            return voigt_exp_indefinite(static_cast<decltype(t_slope)>(x),
                                        static_cast<decltype(t_slope)>(mean),
                                        static_cast<decltype(t_slope)>(sigma),
                                        static_cast<decltype(t_slope)>(gamma),
                                        static_cast<decltype(t_slope)>(tail_ratio),
                                        t_slope);
        },
        tail_slope);
}

// Comprehensive Ceres Jet derivative tests for all parameters
BOOST_AUTO_TEST_CASE(TestJetDerivativesAllParameters) {
    std::cout << "\n[PseudoVoigt] Running TestJetDerivativesAllParameters with Ceres Jet types...\n";
    using ceres::Jet;
    using T = Jet<double, 1>;

    const double h = 1e-6;
    const double abs_tol = 1e-7;
    const double rel_tol = 1e-4; // 0.01% relative tolerance

    auto check_derivative = [&](const std::string& name, auto fn_double, auto fn_jet, double param_base) {
        // Central finite difference reference
        double fwd = fn_double(param_base + h);
        double bwd = fn_double(param_base - h);
        double fd = (fwd - bwd) / (2 * h);

        // Jet-based derivative
        T p;
        p.a = param_base;
        p.v[0] = 1.0;
        T val = fn_jet(p);
        double jet_deriv = val.v[0];

        double mag = std::max(1.0, std::abs(fd));
        double tol = std::max(abs_tol, rel_tol * mag);
        double error = std::abs(jet_deriv - fd);

        BOOST_CHECK_MESSAGE(error < tol,
            name << " derivative: Jet=" << jet_deriv << ", FD=" << fd <<
            ", error=" << error << ", tol=" << tol);

        return error < tol;
    };

    // Test suite 1: Normal parameter ranges
    {
        const double mean = 100.0;
        const double sigma = 1.0;
        const double gamma = 0.1;
        const double tail_ratio = 0.1;
        const double tail_slope = 2.0;
        const double x = 101.0;

        // Test voigt_pdf derivatives
        check_derivative("voigt_pdf d/d(mean)",
            [&](double m) { return voigt_pdf(x, m, sigma, gamma); },
            [&](T m) { return voigt_pdf(T(x), m, T(sigma), T(gamma)); },
            mean);

        check_derivative("voigt_pdf d/d(sigma)",
            [&](double s) { return voigt_pdf(x, mean, s, gamma); },
            [&](T s) { return voigt_pdf(T(x), T(mean), s, T(gamma)); },
            sigma);

        check_derivative("voigt_pdf d/d(gamma)",
            [&](double g) { return voigt_pdf(x, mean, sigma, g); },
            [&](T g) { return voigt_pdf(T(x), T(mean), T(sigma), g); },
            gamma);

        // Test voigt_cdf derivatives
        check_derivative("voigt_cdf d/d(mean)",
            [&](double m) { return voigt_cdf(x, m, sigma, gamma); },
            [&](T m) { return voigt_cdf(T(x), m, T(sigma), T(gamma)); },
            mean);

        check_derivative("voigt_cdf d/d(sigma)",
            [&](double s) { return voigt_cdf(x, mean, s, gamma); },
            [&](T s) { return voigt_cdf(T(x), T(mean), s, T(gamma)); },
            sigma);

        check_derivative("voigt_cdf d/d(gamma)",
            [&](double g) { return voigt_cdf(x, mean, sigma, g); },
            [&](T g) { return voigt_cdf(T(x), T(mean), T(sigma), g); },
            gamma);

        // Test voigt_exp_indefinite derivatives (all 5 parameters)
        check_derivative("voigt_exp_indefinite d/d(mean)",
            [&](double m) { return voigt_exp_indefinite(x, m, sigma, gamma, tail_ratio, tail_slope); },
            [&](T m) { return voigt_exp_indefinite(T(x), m, T(sigma), T(gamma), T(tail_ratio), T(tail_slope)); },
            mean);

        check_derivative("voigt_exp_indefinite d/d(sigma)",
            [&](double s) { return voigt_exp_indefinite(x, mean, s, gamma, tail_ratio, tail_slope); },
            [&](T s) { return voigt_exp_indefinite(T(x), T(mean), s, T(gamma), T(tail_ratio), T(tail_slope)); },
            sigma);

        check_derivative("voigt_exp_indefinite d/d(gamma)",
            [&](double g) { return voigt_exp_indefinite(x, mean, sigma, g, tail_ratio, tail_slope); },
            [&](T g) { return voigt_exp_indefinite(T(x), T(mean), T(sigma), g, T(tail_ratio), T(tail_slope)); },
            gamma);

        check_derivative("voigt_exp_indefinite d/d(tail_ratio)",
            [&](double tr) { return voigt_exp_indefinite(x, mean, sigma, gamma, tr, tail_slope); },
            [&](T tr) { return voigt_exp_indefinite(T(x), T(mean), T(sigma), T(gamma), tr, T(tail_slope)); },
            tail_ratio);

        check_derivative("voigt_exp_indefinite d/d(tail_slope)",
            [&](double ts) { return voigt_exp_indefinite(x, mean, sigma, gamma, tail_ratio, ts); },
            [&](T ts) { return voigt_exp_indefinite(T(x), T(mean), T(sigma), T(gamma), T(tail_ratio), ts); },
            tail_slope);
    }

    // Test suite 2: Extreme parameters (gamma >> sigma, as used in InterSpec)
    {
        const double mean = 14.0;
        const double sigma = 0.0005;
        const double gamma = 0.006;  // gamma/sigma = 12 (Lorentzian-dominated)
        const double tail_ratio = 0.15;
        const double tail_slope = 1.0;
        const double x = 14.001;

        // Critical test: derivatives should still be accurate when using pseudo-Voigt approximation
        check_derivative("voigt_cdf d/d(mean) [gamma>>sigma]",
            [&](double m) { return voigt_cdf(x, m, sigma, gamma); },
            [&](T m) { return voigt_cdf(T(x), m, T(sigma), T(gamma)); },
            mean);

        check_derivative("voigt_cdf d/d(gamma) [gamma>>sigma]",
            [&](double g) { return voigt_cdf(x, mean, sigma, g); },
            [&](T g) { return voigt_cdf(T(x), T(mean), T(sigma), g); },
            gamma);

        check_derivative("voigt_exp_indefinite d/d(sigma) [gamma>>sigma]",
            [&](double s) { return voigt_exp_indefinite(x, mean, s, gamma, tail_ratio, tail_slope); },
            [&](T s) { return voigt_exp_indefinite(T(x), T(mean), s, T(gamma), T(tail_ratio), T(tail_slope)); },
            sigma);

        check_derivative("voigt_exp_indefinite d/d(gamma) [gamma>>sigma]",
            [&](double g) { return voigt_exp_indefinite(x, mean, sigma, g, tail_ratio, tail_slope); },
            [&](T g) { return voigt_exp_indefinite(T(x), T(mean), T(sigma), g, T(tail_ratio), T(tail_slope)); },
            gamma);
    }

    // Test suite 3: Edge cases
    {
        // Very small gamma (nearly Gaussian)
        const double mean = 100.0;
        const double sigma = 1.0;
        const double gamma_tiny = 1e-8;
        const double x = 100.5;

        check_derivative("voigt_pdf d/d(sigma) [gamma~0]",
            [&](double s) { return voigt_pdf(x, mean, s, gamma_tiny); },
            [&](T s) { return voigt_pdf(T(x), T(mean), s, T(gamma_tiny)); },
            sigma);

        check_derivative("voigt_cdf d/d(mean) [gamma~0]",
            [&](double m) { return voigt_cdf(x, m, sigma, gamma_tiny); },
            [&](T m) { return voigt_cdf(T(x), m, T(sigma), T(gamma_tiny)); },
            mean);
    }

    // Test suite 4: Array integration with Jet types
    {
        const double mean = 100.0;
        const double sigma = 1.0;
        const double amplitude = 1000.0;
        const double gamma = 0.1;
        const double tail_ratio = 0.1;
        const double tail_slope = 2.0;

        const size_t num_channels = 10;
        std::vector<float> energies(num_channels + 1);
        for (size_t i = 0; i <= num_channels; ++i) {
            energies[i] = 95.0f + i * 1.0f;
        }

        // Test that array integration works with Jet types (basic functionality)
        std::vector<T> channels_jet(num_channels, T(0.0));
        voigt_exp_integral(T(mean), T(sigma), T(amplitude), T(gamma), 
                          T(tail_ratio), T(tail_slope),
                          energies.data(), channels_jet.data(), num_channels);

        // Verify all channels have valid values
        for (size_t i = 0; i < num_channels; ++i) {
            BOOST_CHECK(std::isfinite(channels_jet[i].a));  // Check scalar part
            BOOST_CHECK(channels_jet[i].a >= 0.0);
            // Verify derivative part is also finite (may be zero if energies are fixed)
            BOOST_CHECK(std::isfinite(channels_jet[i].v[0]));
        }

        // Test derivative w.r.t. mean using array integration
        // Note: Since energies are fixed floats, derivatives w.r.t. mean should propagate
        // through the CDF evaluations
        T mean_jet;
        mean_jet.a = mean;
        mean_jet.v[0] = 1.0;
        
        std::vector<T> channels_jet_mean(num_channels, T(0.0));
        voigt_exp_integral(mean_jet, T(sigma), T(amplitude), T(gamma), 
                          T(tail_ratio), T(tail_slope),
                          energies.data(), channels_jet_mean.data(), num_channels);

        const double h = 1e-5;
        std::vector<double> channels_fwd(num_channels, 0.0);
        std::vector<double> channels_bwd(num_channels, 0.0);
        
        voigt_exp_integral(mean + h, sigma, amplitude, gamma, tail_ratio, tail_slope,
                          energies.data(), channels_fwd.data(), num_channels);
        voigt_exp_integral(mean - h, sigma, amplitude, gamma, tail_ratio, tail_slope,
                          energies.data(), channels_bwd.data(), num_channels);

        // Compare finite difference with Jet derivative for a few channels
        // Use more lenient tolerance for pseudo-Voigt approximation
        for (size_t i = 2; i < num_channels - 2; ++i) {
            double fd_deriv = (channels_fwd[i] - channels_bwd[i]) / (2 * h);
            double jet_deriv = channels_jet_mean[i].v[0];
            double error = std::abs(jet_deriv - fd_deriv);
            // More lenient tolerance for pseudo-Voigt (1% relative or absolute)
            double tol = std::max(1e-3, 0.01 * std::abs(fd_deriv));
            
            BOOST_CHECK_MESSAGE(error < tol,
                "Array integration derivative w.r.t. mean at channel " << i <<
                ": Jet=" << jet_deriv << ", FD=" << fd_deriv <<
                ", error=" << error << ", tol=" << tol);
        }
    }
    std::cout << "[PseudoVoigt] TestJetDerivativesAllParameters completed successfully.\n";
}
#else
BOOST_AUTO_TEST_CASE(TestJetCompatibility) {
    std::cout << "\n[PseudoVoigt] SKIPPING TestJetCompatibility (ceres/jet.h not available).\n";
}

BOOST_AUTO_TEST_CASE(TestJetDerivativesAllParameters) {
    std::cout << "\n[PseudoVoigt] SKIPPING TestJetDerivativesAllParameters (ceres/jet.h not available).\n";
}
#endif // HAS_CERES_JET

BOOST_AUTO_TEST_SUITE_END()

} // namespace
