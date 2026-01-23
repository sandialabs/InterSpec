#include <cmath>

#define BOOST_TEST_MODULE VoigtDistributionTests
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#if __has_include("ceres/jet.h") && __has_include("eigen3/Eigen/Core")
#  define HAS_CERES_JET 1
#  include "eigen3/Eigen/Core"
#  include "ceres/jet.h"
#else
#  define HAS_CERES_JET 0
#endif


#include "../voigt_exp_tail.hpp"
#include "../Faddeeva.hpp"
#include <numeric>
#include <vector>

using namespace voigt_exp_tail;

// Forward declaration of MIT Faddeeva-based MIT implementation from test_voigt_against_mit_imp.cpp
extern void voigt_exp_integral_mit_gauss_legendre(double peak_mean, double sigma_gauss,
                                                  double peak_amplitude, double gamma_lor,
                                                  double tail_ratio, double tail_slope,
                                                  const float * const energies, double *channels,
                                                  size_t nchannel);

// Include scipy-generated Voigt PDF test data
#include "scipy_voigt_test_data.hpp"

namespace {

BOOST_AUTO_TEST_SUITE(VoigtDistributionTestSuite)


BOOST_AUTO_TEST_CASE(TestPureVoigt) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;  // Pure Gaussian
    const double tail_ratio = 0.0;
    const double tail_slope = 1.0;

    // Test at mean
    double pdf_at_mean = voigt_pdf<double>(mean, mean, sigma, gamma);
    double expected_gaussian = 1.0 / (sigma * std::sqrt(2.0 * M_PI));

    BOOST_CHECK_CLOSE(pdf_at_mean, expected_gaussian, 1e-6);
}


BOOST_AUTO_TEST_CASE(TestSymmetry) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.0;  // Pure Voigt, should be symmetric
    const double tail_slope = 1.0;

    double offset = 2.0;
    double pdf_left = voigt_pdf<double>(mean - offset, mean, sigma, gamma);
    double pdf_right = voigt_pdf<double>(mean + offset, mean, sigma, gamma);

    BOOST_CHECK_CLOSE(pdf_left, pdf_right, 1e-6);
}


BOOST_AUTO_TEST_CASE(TestPureGaussian) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;  // No Lorentzian
    const double tail_ratio = 0.0;  // No tail
    const double tail_slope = 1.0;

    // Test PDF at mean
    double pdf_at_mean = voigt_pdf<double>(mean, mean, sigma, gamma);
    double expected_gaussian = 1.0 / (sigma * std::sqrt(2.0 * M_PI));

    BOOST_CHECK_CLOSE(pdf_at_mean, expected_gaussian, 1e-6);

    // Test symmetry
    double offset = 2.0;
    double pdf_left = voigt_pdf<double>(mean - offset, mean, sigma, gamma);
    double pdf_right = voigt_pdf<double>(mean + offset, mean, sigma, gamma);

    BOOST_CHECK_CLOSE(pdf_left, pdf_right, 1e-6);
}

BOOST_AUTO_TEST_CASE(TestAgainstExactFaddeeva1) {
    // Test 1: Full VoigtExpTail (main test case)
    const double mean = 94.67;
    const double sigma = 0.5;
    const double gamma = 0.0876;
    const double tail_ratio = 0.1;
    const double tail_slope = 1.0;
    const double amplitude = 1.0;
    const double energy_min = 80.0;
    const double energy_max = 110.0;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((energy_max - energy_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(energy_min + i * bin_width);
    }

    // Generate spectrum using templated implementation
    std::vector<double> computed_counts(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                       energies.data(), computed_counts.data(), nchannel);

    // Generate exact reference using MIT Faddeeva
    std::vector<double> exact_counts(nchannel, 0.0);
    voigt_exp_integral_mit_gauss_legendre(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                                         energies.data(), exact_counts.data(), nchannel);

    // Check total counts match
    double total_computed = 0.0;
    double total_exact = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_computed += computed_counts[i];
        total_exact += exact_counts[i];
    }

    BOOST_CHECK_CLOSE(total_computed, total_exact, 0.01); // 0.01% tolerance
}

BOOST_AUTO_TEST_CASE(TestAgainstExactFaddeeva2) {
    // Test 2: Pure Voigt (no tail)
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.0;
    const double tail_slope = 1.0;
    const double amplitude = 1.0;
    const double energy_min = 90.0;
    const double energy_max = 110.0;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((energy_max - energy_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(energy_min + i * bin_width);
    }

    std::vector<double> computed_counts(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                       energies.data(), computed_counts.data(), nchannel);

    std::vector<double> exact_counts(nchannel, 0.0);
    voigt_exp_integral_mit_gauss_legendre(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                                         energies.data(), exact_counts.data(), nchannel);

    double total_computed = 0.0;
    double total_exact = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_computed += computed_counts[i];
        total_exact += exact_counts[i];
    }

    BOOST_CHECK_CLOSE(total_computed, total_exact, 0.01); // 0.01% tolerance
}

BOOST_AUTO_TEST_CASE(TestAgainstExactFaddeeva3) {
    // Test 3: Gaussian with tail (no Lorentzian)
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;
    const double tail_ratio = 0.15;
    const double tail_slope = 1.0;
    const double amplitude = 1.0;
    const double energy_min = 90.0;
    const double energy_max = 110.0;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((energy_max - energy_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(energy_min + i * bin_width);
    }

    std::vector<double> computed_counts(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                       energies.data(), computed_counts.data(), nchannel);

    std::vector<double> exact_counts(nchannel, 0.0);
    voigt_exp_integral_mit_gauss_legendre(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                                         energies.data(), exact_counts.data(), nchannel);

    double total_computed = 0.0;
    double total_exact = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_computed += computed_counts[i];
        total_exact += exact_counts[i];
    }

    BOOST_CHECK_CLOSE(total_computed, total_exact, 0.01); // 0.01% tolerance
}

BOOST_AUTO_TEST_CASE(TestAgainstExactFaddeeva4) {
    // Test 4: Large Lorentzian, no tail
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.5;
    const double tail_ratio = 0.0;
    const double tail_slope = 1.0;
    const double amplitude = 1.0;
    const double energy_min = 90.0;
    const double energy_max = 110.0;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((energy_max - energy_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(energy_min + i * bin_width);
    }

    std::vector<double> computed_counts(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                       energies.data(), computed_counts.data(), nchannel);

    std::vector<double> exact_counts(nchannel, 0.0);
    voigt_exp_integral_mit_gauss_legendre(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                                         energies.data(), exact_counts.data(), nchannel);

    double total_computed = 0.0;
    double total_exact = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_computed += computed_counts[i];
        total_exact += exact_counts[i];
    }

    BOOST_CHECK_CLOSE(total_computed, total_exact, 0.01); // 0.01% tolerance
}

BOOST_AUTO_TEST_CASE(TestAgainstExactFaddeeva5) {
    // Test 5: Pure Gaussian
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;
    const double tail_ratio = 0.0;
    const double tail_slope = 1.0;
    const double amplitude = 1.0;
    const double energy_min = 90.0;
    const double energy_max = 110.0;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((energy_max - energy_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(energy_min + i * bin_width);
    }

    std::vector<double> computed_counts(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                       energies.data(), computed_counts.data(), nchannel);

    std::vector<double> exact_counts(nchannel, 0.0);
    voigt_exp_integral_mit_gauss_legendre(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                                         energies.data(), exact_counts.data(), nchannel);

    double total_computed = 0.0;
    double total_exact = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_computed += computed_counts[i];
        total_exact += exact_counts[i];
    }

    BOOST_CHECK_CLOSE(total_computed, total_exact, 0.01); // 0.01% tolerance
}

BOOST_AUTO_TEST_CASE(TestGammaSpectroscopyParameters) {
    // Test with realistic gamma spectroscopy parameters
    const double mean = 94.67;      // Uranium 94.67 keV x-ray
    const double sigma = 0.5;       // 500 eV Gaussian width (typical HPGe)
    const double gamma = 0.0876;    // 87.6 eV Lorentzian width (realistic)
    const double tail_ratio = 0.05; // 5% tailing
    const double tail_slope = 1.0;  // Tail slope

    // Test PDF values at peak
    double pdf_at_peak = voigt_pdf<double>(mean, mean, sigma, gamma);
    BOOST_CHECK(pdf_at_peak > 0.0);

    // Test PDF values are non-negative and decrease away from peak
    double pdf_left = voigt_exp_pdf<double>(mean - 1.0, mean, sigma, gamma, tail_ratio, tail_slope);
    double pdf_center = voigt_exp_pdf<double>(mean, mean, sigma, gamma, tail_ratio, tail_slope);
    double pdf_right = voigt_exp_pdf<double>(mean + 1.0, mean, sigma, gamma, tail_ratio, tail_slope);

    BOOST_CHECK(pdf_left >= 0.0);
    BOOST_CHECK(pdf_center >= 0.0);
    BOOST_CHECK(pdf_right >= 0.0);
    BOOST_CHECK(pdf_center > pdf_left);  // PDF should be higher at peak than 1 keV away
    BOOST_CHECK(pdf_center > pdf_right); // PDF should be higher at peak than 1 keV away
}


BOOST_AUTO_TEST_CASE(TestJetCompatibility) {
#if !HAS_CERES_JET
    std::cout << "Skipping Jet compatibility test (ceres/jet.h not available).\n";
    return;
#else
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
    const double abs_tol = 1e-8;
    const double rel_tol = 1e-5; // 10 ppm relative tolerance (0.001%)

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

    // Test voigt_exp_pdf derivatives
    // Derivative w.r.t. mean
    check_derivative(
        [&](auto m) {
            return voigt_exp_tail::voigt_exp_pdf(static_cast<decltype(m)>(x),
                                                 m,
                                                 static_cast<decltype(m)>(sigma),
                                                 static_cast<decltype(m)>(gamma),
                                                 static_cast<decltype(m)>(tail_ratio),
                                                 static_cast<decltype(m)>(tail_slope));
        },
        mean);

    // Derivative w.r.t. sigma
    check_derivative(
        [&](auto s) {
            return voigt_exp_tail::voigt_exp_pdf(static_cast<decltype(s)>(x),
                                                 static_cast<decltype(s)>(mean),
                                                 s,
                                                 static_cast<decltype(s)>(gamma),
                                                 static_cast<decltype(s)>(tail_ratio),
                                                 static_cast<decltype(s)>(tail_slope));
        },
        sigma);

    // Derivative w.r.t. gamma
    check_derivative(
        [&](auto g) {
            return voigt_exp_tail::voigt_exp_pdf(static_cast<decltype(g)>(x),
                                                 static_cast<decltype(g)>(mean),
                                                 static_cast<decltype(g)>(sigma),
                                                 g,
                                                 static_cast<decltype(g)>(tail_ratio),
                                                 static_cast<decltype(g)>(tail_slope));
        },
        gamma);

    // Derivative w.r.t. tail_ratio
    check_derivative(
        [&](auto t_ratio) {
            return voigt_exp_tail::voigt_exp_pdf(static_cast<decltype(t_ratio)>(x),
                                                 static_cast<decltype(t_ratio)>(mean),
                                                 static_cast<decltype(t_ratio)>(sigma),
                                                 static_cast<decltype(t_ratio)>(gamma),
                                                 t_ratio,
                                                 static_cast<decltype(t_ratio)>(tail_slope));
        },
        tail_ratio);

    // Derivative w.r.t. tail_slope
    check_derivative(
        [&](auto t_slope) {
            return voigt_exp_tail::voigt_exp_pdf(static_cast<decltype(t_slope)>(x),
                                                 static_cast<decltype(t_slope)>(mean),
                                                 static_cast<decltype(t_slope)>(sigma),
                                                 static_cast<decltype(t_slope)>(gamma),
                                                 static_cast<decltype(t_slope)>(tail_ratio),
                                                 t_slope);
        },
        tail_slope);
#endif // HAS_CERES_JET
}


BOOST_AUTO_TEST_CASE(TestJetCompatibilityIntegral) {
#if !HAS_CERES_JET
    std::cout << "Skipping Jet compatibility test for voigt_exp_integral (ceres/jet.h not available).\n";
    return;
#else
    using ceres::Jet;
    using T = Jet<double, 1>;

    // Realistic gamma spectroscopy parameters (HPGe around 95 keV)
    const double mean = 94.67;
    const double sigma = 0.5;
    const double gamma = 0.0876;
    const double tail_ratio = 0.05;
    const double tail_slope = 1.5;
    const double amplitude = 1.0;

    // Create energy bins around the peak
    const double energy_min = 90.0;
    const double energy_max = 100.0;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((energy_max - energy_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(energy_min + i * bin_width);
    }

    const double h = 1e-5;
    const double abs_tol = 1e-7;
    const double rel_tol = 1e-5; // 10 ppm relative tolerance (0.001%)

    // Helper to compute total integral (sum of all channels) for double parameters
    auto compute_total_integral_double = [&](double p_mean, double p_sigma, double p_gamma, 
                                             double p_tail_ratio, double p_tail_slope) -> double {
        std::vector<double> channels(nchannel, 0.0);
        voigt_exp_tail::voigt_exp_integral(p_mean, p_sigma, amplitude, p_gamma, p_tail_ratio, p_tail_slope,
                                           energies.data(), channels.data(), nchannel);
        double total = 0.0;
        for (size_t i = 0; i < nchannel; ++i) {
            total += channels[i];
        }
        return total;
    };

    // Helper lambdas to compute total integral for Jet parameters (one per parameter)
    auto compute_total_integral_jet_mean = [&](const T& param_val) -> T {
        std::vector<T> channels(nchannel, T(0));
        voigt_exp_tail::voigt_exp_integral(param_val, T(sigma), T(amplitude), T(gamma), T(tail_ratio), T(tail_slope),
                                           energies.data(), channels.data(), nchannel);
        T total = T(0);
        for (size_t i = 0; i < nchannel; ++i) total += channels[i];
        return total;
    };
    
    auto compute_total_integral_jet_sigma = [&](const T& param_val) -> T {
        std::vector<T> channels(nchannel, T(0));
        voigt_exp_tail::voigt_exp_integral(T(mean), param_val, T(amplitude), T(gamma), T(tail_ratio), T(tail_slope),
                                           energies.data(), channels.data(), nchannel);
        T total = T(0);
        for (size_t i = 0; i < nchannel; ++i) total += channels[i];
        return total;
    };
    
    auto compute_total_integral_jet_gamma = [&](const T& param_val) -> T {
        std::vector<T> channels(nchannel, T(0));
        voigt_exp_tail::voigt_exp_integral(T(mean), T(sigma), T(amplitude), param_val, T(tail_ratio), T(tail_slope),
                                           energies.data(), channels.data(), nchannel);
        T total = T(0);
        for (size_t i = 0; i < nchannel; ++i) total += channels[i];
        return total;
    };
    
    auto compute_total_integral_jet_tail_ratio = [&](const T& param_val) -> T {
        std::vector<T> channels(nchannel, T(0));
        voigt_exp_tail::voigt_exp_integral(T(mean), T(sigma), T(amplitude), T(gamma), param_val, T(tail_slope),
                                           energies.data(), channels.data(), nchannel);
        T total = T(0);
        for (size_t i = 0; i < nchannel; ++i) total += channels[i];
        return total;
    };
    
    auto compute_total_integral_jet_tail_slope = [&](const T& param_val) -> T {
        std::vector<T> channels(nchannel, T(0));
        voigt_exp_tail::voigt_exp_integral(T(mean), T(sigma), T(amplitude), T(gamma), T(tail_ratio), param_val,
                                           energies.data(), channels.data(), nchannel);
        T total = T(0);
        for (size_t i = 0; i < nchannel; ++i) total += channels[i];
        return total;
    };

    auto check_derivative = [&](auto compute_jet, double param_base, const char* param_name,
                                double p_mean_fwd, double p_sigma_fwd, double p_gamma_fwd, 
                                double p_tail_ratio_fwd, double p_tail_slope_fwd,
                                double p_mean_bwd, double p_sigma_bwd, double p_gamma_bwd, 
                                double p_tail_ratio_bwd, double p_tail_slope_bwd) {
        // Central finite difference reference
        double fwd = compute_total_integral_double(p_mean_fwd, p_sigma_fwd, p_gamma_fwd, p_tail_ratio_fwd, p_tail_slope_fwd);
        double bwd = compute_total_integral_double(p_mean_bwd, p_sigma_bwd, p_gamma_bwd, p_tail_ratio_bwd, p_tail_slope_bwd);
        double fd = (fwd - bwd) / (2 * h);

        // Jet-based derivative
        T p;
        p.a = param_base;
        p.v[0] = 1.0;
        T total_jet = compute_jet(p);
        double jet_deriv = total_jet.v[0];

        double mag = std::max(1.0, std::abs(fd));
        double tol = std::max(abs_tol, rel_tol * mag);
        
        BOOST_TEST_MESSAGE("Testing derivative w.r.t. " << param_name 
                          << ": FD=" << fd << ", Jet=" << jet_deriv 
                          << ", error=" << std::abs(jet_deriv - fd) << ", tol=" << tol);
        BOOST_CHECK_SMALL(jet_deriv - fd, tol);
    };

    // Test derivatives for each parameter
    check_derivative(compute_total_integral_jet_mean, mean, "mean",
                     mean + h, sigma, gamma, tail_ratio, tail_slope,
                     mean - h, sigma, gamma, tail_ratio, tail_slope);
    check_derivative(compute_total_integral_jet_sigma, sigma, "sigma",
                     mean, sigma + h, gamma, tail_ratio, tail_slope,
                     mean, sigma - h, gamma, tail_ratio, tail_slope);
    check_derivative(compute_total_integral_jet_gamma, gamma, "gamma",
                     mean, sigma, gamma + h, tail_ratio, tail_slope,
                     mean, sigma, gamma - h, tail_ratio, tail_slope);
    check_derivative(compute_total_integral_jet_tail_ratio, tail_ratio, "tail_ratio",
                     mean, sigma, gamma, tail_ratio + h, tail_slope,
                     mean, sigma, gamma, tail_ratio - h, tail_slope);
    check_derivative(compute_total_integral_jet_tail_slope, tail_slope, "tail_slope",
                     mean, sigma, gamma, tail_ratio, tail_slope + h,
                     mean, sigma, gamma, tail_ratio, tail_slope - h);
#endif // HAS_CERES_JET
}


BOOST_AUTO_TEST_CASE(TestEdgeCases) {
    const double mean = 100.0;
    const double sigma = 1.0;

    // Test with very small gamma (should behave like Gaussian)
    double gamma_small = 1e-10;
    double pdf_gauss = voigt_pdf<double>(mean, mean, sigma, 0.0);
    double pdf_small_gamma = voigt_pdf<double>(mean, mean, sigma, gamma_small);
    BOOST_CHECK_CLOSE(pdf_gauss, pdf_small_gamma, 1e-6);

    // Test with very small sigma (should behave like Lorentzian)
    double sigma_small = 1e-10;
    double gamma_test = 0.1;
    double x_test = mean + 1.0;

    // For Lorentzian, PDF should be gamma/(pi * (x-mean)^2 + gamma^2)
    double expected_lorentz = gamma_test / (3.14159265358979323846264338327950288 * (gamma_test * gamma_test));
    double pdf_small_sigma = voigt_pdf<double>(mean, mean, sigma_small, gamma_test);
    BOOST_CHECK_CLOSE(pdf_small_sigma, expected_lorentz, 1e-3);

    // Test with zero tail ratio (should be pure Voigt)
    double tail_ratio_zero = 0.0;
    double pdf_voigt = voigt_pdf<double>(mean, mean, sigma, gamma_test);
    double pdf_voigt_with_tail = voigt_exp_pdf<double>(mean, mean, sigma, gamma_test, tail_ratio_zero, 1.0);
    BOOST_CHECK_CLOSE(pdf_voigt_with_tail, pdf_voigt, 1e-6);
}

// Test case from InterSpec that tests PDF functions with extreme parameters
// This test verifies that PDF functions work correctly when gamma >> sigma
BOOST_AUTO_TEST_CASE(TestExtremeParametersFromInterSpec) {
    // Test case for extreme parameters from InterSpec test_PeakDists.cpp:
    // - sigma = 0.0005 keV (0.5 eV) - very narrow Gaussian
    // - gamma = 0.006 keV (6 eV) --> gamma/sigma = 12 (Lorentzian-dominated)
    // - tail_slope = 1.0 keV (1000 eV) --> tail_slope/sigma = 2000 (extreme tail)
    //
    // When gamma >> sigma, the complex argument z = ((x - mean) + i*gamma) / (sigma * sqrt(2))
    // has a large imaginary component (z_imag ≈ 8.49), which can cause numerical instability
    // in the Faddeeva function computation. PDF functions should handle this correctly.
    const double mean = 14.0;
    const double sigma = 0.0005;
    const double gamma = 0.006;  // gamma >> sigma (Lorentzian-dominated regime)
    const double tail_ratio = 0.15;
    const double tail_slope = 1.0;  // tail_slope >> sigma (extreme tail)

    // Test 1: voigt_pdf should return non-negative, finite values
    double x_near_mean = 13.9970874786377;  // Slightly below mean
    double pdf_near_mean = voigt_pdf<double>(x_near_mean, mean, sigma, gamma);

    BOOST_CHECK_MESSAGE(pdf_near_mean >= 0.0 && std::isfinite(pdf_near_mean),
        "voigt_pdf returns " << pdf_near_mean << " for gamma >> sigma case. "
        "Expected non-negative, finite value.");

    // Test 2: voigt_exp_pdf should also return non-negative, finite values
    double pdf_combined = voigt_exp_pdf<double>(x_near_mean, mean, sigma, gamma, tail_ratio, tail_slope);
    BOOST_CHECK_MESSAGE(pdf_combined >= 0.0 && std::isfinite(pdf_combined),
        "voigt_exp_pdf returns " << pdf_combined <<
        " for extreme parameters. Expected non-negative, finite value.");

    // Test 3: PDF should be highest at the peak
    double pdf_at_peak = voigt_exp_pdf<double>(mean, mean, sigma, gamma, tail_ratio, tail_slope);
    double pdf_left = voigt_exp_pdf<double>(mean - 1.0, mean, sigma, gamma, tail_ratio, tail_slope);
    double pdf_right = voigt_exp_pdf<double>(mean + 1.0, mean, sigma, gamma, tail_ratio, tail_slope);

    BOOST_CHECK_MESSAGE(pdf_at_peak > pdf_left,
        "PDF should be higher at peak than 1 keV to the left. "
        "pdf(peak) = " << pdf_at_peak << ", pdf(peak-1) = " << pdf_left);
    BOOST_CHECK_MESSAGE(pdf_at_peak > pdf_right,
        "PDF should be higher at peak than 1 keV to the right. "
        "pdf(peak) = " << pdf_at_peak << ", pdf(peak+1) = " << pdf_right);

    // Test 4: Verify that voigt_pdf handles Lorentzian-dominated case correctly
    // The complex argument z = ((x - mean) + i*gamma) / (sigma * sqrt(2)) has:
    // - z_real ≈ -4.12 (x slightly below mean)
    // - z_imag ≈ 8.49 (gamma/sigma * sqrt(2) = 12 * sqrt(2) ≈ 17, but divided by sqrt(2) gives ~8.49)
    double sqrt2 = 1.41421356237309504880;
    double z_real = (x_near_mean - mean) / (sigma * sqrt2);
    double z_imag = gamma / (sigma * sqrt2);
    
    // Verify z_imag is large enough to potentially cause numerical issues
    // For this test case: z_imag = gamma/(sigma*sqrt(2)) = 0.006/(0.0005*1.414) ≈ 8.49
    BOOST_CHECK_MESSAGE(z_imag > 5.0,
        "z_imag = " << z_imag << " should be large enough to potentially cause numerical issues. "
        "For gamma/sigma = " << (gamma/sigma) << ", z_imag = " << z_imag);
    
    // Verify that PDF is close to Lorentzian PDF for Lorentzian-dominated case
    // Since gamma/sigma = 12 > 5, the PDF should be close to Lorentzian
    double lorentzian_pdf_val = lorentzian_pdf<double>(x_near_mean, mean, gamma);
    // For Lorentzian-dominated case, voigt_pdf should be close to lorentzian_pdf
    // Allow some tolerance since it's not exactly Lorentzian
    BOOST_CHECK_MESSAGE(std::abs(pdf_near_mean - lorentzian_pdf_val) / lorentzian_pdf_val < 0.1,
        "voigt_pdf should be close to Lorentzian PDF when gamma/sigma > 5. "
        "voigt_pdf = " << pdf_near_mean << ", lorentzian_pdf = " << lorentzian_pdf_val);
}

// This test verifies PDF-based integration for extreme parameters (gamma >> sigma with tail_ratio > 0).
// Tests use PDF functions since CDF-based functions (voigt_exp_indefinite, voigt_exp_coverage_limits)
// have been removed. The array version of voigt_exp_integral uses PDF integration via Gauss-Legendre quadrature.
BOOST_AUTO_TEST_CASE(TestInterSpecFailures) {
    // Same extreme parameters as TestExtremeParametersFromInterSpec
    const double mean = 14.0;
    const double sigma = 0.0005;
    const double gamma = 0.006;  // gamma/sigma = 12 (Lorentzian-dominated)
    const double tail_ratio = 0.15;
    const double tail_slope = 1.0;
    const double amplitude = 1.23456;

    // Test 1: Array integration (voigt_exp_integral with energy bins) should give correct amplitude
    // This test verifies that PDF-based integration correctly captures the full distribution.
    // Uses a wide range to ensure full coverage of the distribution.
    {
        // Use a wide range that should cover the full distribution
        // For Lorentzian-dominated: need ~20*gamma on each side, plus tail extent
        double range = std::max(20.0 * gamma, tail_slope * 10.0);
        double x0 = mean - range;
        double x1 = mean + range;

        size_t num_channels = 2048;
        std::vector<double> counts(num_channels, 0.0);
        std::vector<float> energies(num_channels + 1, 0.0f);
        for (size_t i = 0; i < energies.size(); ++i)
            energies[i] = static_cast<float>(x0 + (i * (x1 - x0) / num_channels));

        voigt_exp_integral(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                          &energies[0], &counts[0], num_channels);

        double answer_sum = std::accumulate(counts.begin(), counts.end(), 0.0);
        double error_fraction = std::abs(answer_sum - amplitude) / amplitude;

        BOOST_CHECK_MESSAGE(error_fraction < 1.0e-3,
            "Array integration should give correct amplitude. "
            "Expected: " << amplitude << ", Actual: " << answer_sum <<
            ", Error: " << (error_fraction * 100.0) << "%. "
            "PDF-based integration should capture full distribution.");
    }

    // Test 2: PDF values should be non-negative and finite
    // This test verifies that PDF functions work correctly for extreme parameters
    {
        double pdf_at_peak = voigt_exp_pdf<double>(mean, mean, sigma, gamma, tail_ratio, tail_slope);
        double pdf_left = voigt_exp_pdf<double>(mean - 1.0, mean, sigma, gamma, tail_ratio, tail_slope);
        double pdf_right = voigt_exp_pdf<double>(mean + 1.0, mean, sigma, gamma, tail_ratio, tail_slope);

        BOOST_CHECK(pdf_at_peak > 0.0);
        BOOST_CHECK(pdf_left >= 0.0);
        BOOST_CHECK(pdf_right >= 0.0);
        BOOST_CHECK(std::isfinite(pdf_at_peak));
        BOOST_CHECK(std::isfinite(pdf_left));
        BOOST_CHECK(std::isfinite(pdf_right));
    }

    // Test 3: PDF should decrease away from peak
    // This test verifies basic PDF properties for extreme parameters
    {
        double pdf_at_peak = voigt_exp_pdf<double>(mean, mean, sigma, gamma, tail_ratio, tail_slope);
        double pdf_far_left = voigt_exp_pdf<double>(mean - 5.0, mean, sigma, gamma, tail_ratio, tail_slope);
        double pdf_far_right = voigt_exp_pdf<double>(mean + 5.0, mean, sigma, gamma, tail_ratio, tail_slope);

        BOOST_CHECK_MESSAGE(pdf_at_peak > pdf_far_left,
            "PDF should be higher at peak than far left. "
            "pdf(peak) = " << pdf_at_peak << ", pdf(peak-5) = " << pdf_far_left);
        BOOST_CHECK_MESSAGE(pdf_at_peak > pdf_far_right,
            "PDF should be higher at peak than far right. "
            "pdf(peak) = " << pdf_at_peak << ", pdf(peak+5) = " << pdf_far_right);
    }
}

// Comprehensive Ceres Jet derivative tests for all parameters using PDF functions
BOOST_AUTO_TEST_CASE(TestJetDerivativesAllParameters) {
#if !HAS_CERES_JET
    std::cout << "Skipping comprehensive Jet derivative tests (ceres/jet.h not available).\n";
    return;
#else
    using ceres::Jet;
    using T = Jet<double, 1>;

    const double h = 1e-6;
    const double abs_tol = 1e-8;
    const double rel_tol = 1e-5; // 10 ppm relative tolerance (0.001%)

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

        // Test voigt_exp_pdf derivatives (all 5 parameters)
        check_derivative("voigt_exp_pdf d/d(mean)",
            [&](double m) { return voigt_exp_pdf(x, m, sigma, gamma, tail_ratio, tail_slope); },
            [&](T m) { return voigt_exp_pdf(T(x), m, T(sigma), T(gamma), T(tail_ratio), T(tail_slope)); },
            mean);

        check_derivative("voigt_exp_pdf d/d(sigma)",
            [&](double s) { return voigt_exp_pdf(x, mean, s, gamma, tail_ratio, tail_slope); },
            [&](T s) { return voigt_exp_pdf(T(x), T(mean), s, T(gamma), T(tail_ratio), T(tail_slope)); },
            sigma);

        check_derivative("voigt_exp_pdf d/d(gamma)",
            [&](double g) { return voigt_exp_pdf(x, mean, sigma, g, tail_ratio, tail_slope); },
            [&](T g) { return voigt_exp_pdf(T(x), T(mean), T(sigma), g, T(tail_ratio), T(tail_slope)); },
            gamma);

        check_derivative("voigt_exp_pdf d/d(tail_ratio)",
            [&](double tr) { return voigt_exp_pdf(x, mean, sigma, gamma, tr, tail_slope); },
            [&](T tr) { return voigt_exp_pdf(T(x), T(mean), T(sigma), T(gamma), tr, T(tail_slope)); },
            tail_ratio);

        check_derivative("voigt_exp_pdf d/d(tail_slope)",
            [&](double ts) { return voigt_exp_pdf(x, mean, sigma, gamma, tail_ratio, ts); },
            [&](T ts) { return voigt_exp_pdf(T(x), T(mean), T(sigma), T(gamma), T(tail_ratio), ts); },
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

        // Critical test: derivatives should still be accurate when using Lorentzian approximation
        check_derivative("voigt_exp_pdf d/d(mean) [gamma>>sigma]",
            [&](double m) { return voigt_exp_pdf(x, m, sigma, gamma, tail_ratio, tail_slope); },
            [&](T m) { return voigt_exp_pdf(T(x), m, T(sigma), T(gamma), T(tail_ratio), T(tail_slope)); },
            mean);

        check_derivative("voigt_exp_pdf d/d(gamma) [gamma>>sigma]",
            [&](double g) { return voigt_exp_pdf(x, mean, sigma, g, tail_ratio, tail_slope); },
            [&](T g) { return voigt_exp_pdf(T(x), T(mean), T(sigma), g, T(tail_ratio), T(tail_slope)); },
            gamma);

        check_derivative("voigt_exp_pdf d/d(sigma) [gamma>>sigma]",
            [&](double s) { return voigt_exp_pdf(x, mean, s, gamma, tail_ratio, tail_slope); },
            [&](T s) { return voigt_exp_pdf(T(x), T(mean), s, T(gamma), T(tail_ratio), T(tail_slope)); },
            sigma);
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

        // Very small sigma (nearly Lorentzian)
        const double sigma_tiny = 1e-8;
        const double gamma_large = 1.0;

        check_derivative("voigt_pdf d/d(gamma) [sigma~0]",
            [&](double g) { return voigt_pdf(x, mean, sigma_tiny, g); },
            [&](T g) { return voigt_pdf(T(x), T(mean), T(sigma_tiny), g); },
            gamma_large);

        // Zero tail_ratio (pure Voigt, no exponential tail)
        const double tail_ratio_zero = 0.0;

        check_derivative("voigt_exp_pdf d/d(mean) [tail_ratio=0]",
            [&](double m) { return voigt_exp_pdf(x, m, sigma, gamma_large, tail_ratio_zero, 1.0); },
            [&](T m) { return voigt_exp_pdf(T(x), m, T(sigma), T(gamma_large), T(tail_ratio_zero), T(1.0)); },
            mean);

        // Large tail_ratio (tail-dominated)
        const double tail_ratio_large = 0.9;

        check_derivative("voigt_exp_pdf d/d(tail_ratio) [tail_ratio=0.9]",
            [&](double tr) { return voigt_exp_pdf(x, mean, sigma, gamma_large, tr, 1.0); },
            [&](T tr) { return voigt_exp_pdf(T(x), T(mean), T(sigma), T(gamma_large), tr, T(1.0)); },
            tail_ratio_large);
    }
#endif
}

// Test array-filling function with Jets (critical for RelActCalcAuto fitting)
BOOST_AUTO_TEST_CASE(TestJetArrayIntegral) {
#if !HAS_CERES_JET
    std::cout << "Skipping Jet array integral test (ceres/jet.h not available).\n";
    return;
#else
    using ceres::Jet;
    using T = Jet<double, 1>;

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

    // Test derivative w.r.t. amplitude
    {
        const double h = 1e-6;

        // Finite difference
        std::vector<double> counts_fwd(num_channels);
        std::vector<double> counts_bwd(num_channels);
        voigt_exp_integral(mean, sigma, amplitude + h, gamma, tail_ratio, tail_slope,
                          energies.data(), counts_fwd.data(), num_channels);
        voigt_exp_integral(mean, sigma, amplitude - h, gamma, tail_ratio, tail_slope,
                          energies.data(), counts_bwd.data(), num_channels);

        // Jet version
        std::vector<T> counts_jet(num_channels);
        T amp_jet;
        amp_jet.a = amplitude;
        amp_jet.v[0] = 1.0;
        voigt_exp_integral(T(mean), T(sigma), amp_jet, T(gamma), T(tail_ratio), T(tail_slope),
                          energies.data(), counts_jet.data(), num_channels);

        // Check derivatives for each channel
        bool all_passed = true;
        for (size_t i = 0; i < num_channels; ++i) {
            double fd = (counts_fwd[i] - counts_bwd[i]) / (2 * h);
            double jet_deriv = counts_jet[i].v[0];
            double error = std::abs(jet_deriv - fd);
            double tol = std::max(1e-8, 1e-5 * std::abs(fd)); // 0.001% relative tolerance

            if (error > tol) {
                all_passed = false;
                BOOST_TEST_MESSAGE("Channel " << i << ": Jet=" << jet_deriv <<
                                   ", FD=" << fd << ", error=" << error);
            }
        }
        BOOST_CHECK_MESSAGE(all_passed,
            "Array integral derivatives w.r.t. amplitude should match finite difference");
    }

    // Test derivative w.r.t. mean (peak position)
    {
        const double h = 1e-6;

        std::vector<double> counts_fwd(num_channels);
        std::vector<double> counts_bwd(num_channels);
        voigt_exp_integral(mean + h, sigma, amplitude, gamma, tail_ratio, tail_slope,
                          energies.data(), counts_fwd.data(), num_channels);
        voigt_exp_integral(mean - h, sigma, amplitude, gamma, tail_ratio, tail_slope,
                          energies.data(), counts_bwd.data(), num_channels);

        std::vector<T> counts_jet(num_channels);
        T mean_jet;
        mean_jet.a = mean;
        mean_jet.v[0] = 1.0;
        voigt_exp_integral(mean_jet, T(sigma), T(amplitude), T(gamma), T(tail_ratio), T(tail_slope),
                          energies.data(), counts_jet.data(), num_channels);

        bool all_passed = true;
        int failed_count = 0;
        int skipped_count = 0;
        for (size_t i = 0; i < num_channels; ++i) {
            // Skip channels very close to peak mean where derivative w.r.t. mean changes sign sharply
            // Finite differences are numerically unstable at these points
            double channel_center = x_min + (i + 0.5) * (x_max - x_min) / num_channels;
            if (std::abs(channel_center - mean) < 0.5 * sigma) {
                skipped_count++;
                continue; // Skip this channel
            }

            double fd = (counts_fwd[i] - counts_bwd[i]) / (2 * h);
            double jet_deriv = counts_jet[i].v[0];
            double error = std::abs(jet_deriv - fd);
            double tol = std::max(1e-8, 1e-5 * std::abs(fd)); // 0.001% relative tolerance

            if (error > tol) {
                all_passed = false;
                if (failed_count < 5) {  // Only print first 5 failures
                    std::cout << "Channel " << i << " (center=" << channel_center << "): Jet=" << jet_deriv <<
                              ", FD=" << fd << ", error=" << error << ", tol=" << tol << "\n";
                }
                failed_count++;
            }
        }
        if (failed_count > 0) {
            std::cout << "Total failures: " << failed_count << " out of " << (num_channels - skipped_count) << " tested channels ";
            std::cout << "(" << skipped_count << " channels near peak skipped)\n";
        }
        BOOST_CHECK_MESSAGE(all_passed,
            "Array integral derivatives w.r.t. mean should match finite difference");
    }

    // Test derivative w.r.t. sigma in extreme case (gamma >> sigma)
    {
        const double mean_extreme = 14.0;
        const double sigma_extreme = 0.0005;
        const double gamma_extreme = 0.006;
        const double tail_ratio_extreme = 0.15;
        const double tail_slope_extreme = 1.0;
        const double amp_extreme = 100.0;

        const double x_min_extreme = mean_extreme - 0.1;
        const double x_max_extreme = mean_extreme + 0.1;

        std::vector<float> energies_extreme(num_channels + 1);
        for (size_t i = 0; i <= num_channels; ++i) {
            energies_extreme[i] = static_cast<float>(x_min_extreme +
                                  i * (x_max_extreme - x_min_extreme) / num_channels);
        }

        const double h = 1e-8;

        std::vector<double> counts_fwd(num_channels);
        std::vector<double> counts_bwd(num_channels);
        voigt_exp_integral(mean_extreme, sigma_extreme + h, amp_extreme, gamma_extreme,
                          tail_ratio_extreme, tail_slope_extreme,
                          energies_extreme.data(), counts_fwd.data(), num_channels);
        voigt_exp_integral(mean_extreme, sigma_extreme - h, amp_extreme, gamma_extreme,
                          tail_ratio_extreme, tail_slope_extreme,
                          energies_extreme.data(), counts_bwd.data(), num_channels);

        std::vector<T> counts_jet(num_channels);
        T sigma_jet;
        sigma_jet.a = sigma_extreme;
        sigma_jet.v[0] = 1.0;
        voigt_exp_integral(T(mean_extreme), sigma_jet, T(amp_extreme), T(gamma_extreme),
                          T(tail_ratio_extreme), T(tail_slope_extreme),
                          energies_extreme.data(), counts_jet.data(), num_channels);

        bool all_passed = true;
        for (size_t i = 0; i < num_channels; ++i) {
            double fd = (counts_fwd[i] - counts_bwd[i]) / (2 * h);
            double jet_deriv = counts_jet[i].v[0];
            double error = std::abs(jet_deriv - fd);
            double tol = std::max(1e-6, 1e-4 * std::abs(fd)); // 0.01% relative tolerance for extreme case

            if (error > tol) {
                all_passed = false;
            }
        }
        BOOST_CHECK_MESSAGE(all_passed,
            "Array integral derivatives w.r.t. sigma [gamma>>sigma] should match finite difference");
    }
#endif
}

// Multi-parameter Jets test (for Hessian computation) using PDF functions
BOOST_AUTO_TEST_CASE(TestMultiParameterJets) {
#if !HAS_CERES_JET
    std::cout << "Skipping multi-parameter Jet test (ceres/jet.h not available).\n";
    return;
#else
    using ceres::Jet;

    // Use 5-parameter Jets to compute all partial derivatives simultaneously
    using T = Jet<double, 5>;

    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.1;
    const double tail_slope = 2.0;
    const double x = 101.0;

    // Set up Jets with one-hot vectors
    T mean_jet, sigma_jet, gamma_jet, tail_ratio_jet, tail_slope_jet, x_jet;
    mean_jet.a = mean; mean_jet.v[0] = 1.0;
    sigma_jet.a = sigma; sigma_jet.v[1] = 1.0;
    gamma_jet.a = gamma; gamma_jet.v[2] = 1.0;
    tail_ratio_jet.a = tail_ratio; tail_ratio_jet.v[3] = 1.0;
    tail_slope_jet.a = tail_slope; tail_slope_jet.v[4] = 1.0;
    x_jet.a = x; // x is not a fit parameter, so all v[i] = 0

    // Compute voigt_exp_pdf with all parameter derivatives
    T result = voigt_exp_pdf(x_jet, mean_jet, sigma_jet, gamma_jet,
                              tail_ratio_jet, tail_slope_jet);

    // Check that all derivatives are non-zero and finite
    BOOST_CHECK(std::isfinite(result.a));
    BOOST_CHECK(result.v[0] != 0.0); // d/d(mean)
    BOOST_CHECK(std::isfinite(result.v[0]));
    BOOST_CHECK(result.v[1] != 0.0); // d/d(sigma)
    BOOST_CHECK(std::isfinite(result.v[1]));
    BOOST_CHECK(result.v[2] != 0.0); // d/d(gamma)
    BOOST_CHECK(std::isfinite(result.v[2]));
    BOOST_CHECK(std::isfinite(result.v[3])); // d/d(tail_ratio) - could be zero at x far from tail
    BOOST_CHECK(std::isfinite(result.v[4])); // d/d(tail_slope) - could be zero at x far from tail

    // Verify gradient consistency with finite differences
    const double h = 1e-6;

    double fd_mean = (voigt_exp_pdf(x, mean + h, sigma, gamma, tail_ratio, tail_slope) -
                      voigt_exp_pdf(x, mean - h, sigma, gamma, tail_ratio, tail_slope)) / (2*h);
    BOOST_CHECK_CLOSE(result.v[0], fd_mean, 0.01); // 0.01% tolerance

    double fd_sigma = (voigt_exp_pdf(x, mean, sigma + h, gamma, tail_ratio, tail_slope) -
                       voigt_exp_pdf(x, mean, sigma - h, gamma, tail_ratio, tail_slope)) / (2*h);
    BOOST_CHECK_CLOSE(result.v[1], fd_sigma, 0.01); // 0.01% tolerance

    double fd_gamma = (voigt_exp_pdf(x, mean, sigma, gamma + h, tail_ratio, tail_slope) -
                       voigt_exp_pdf(x, mean, sigma, gamma - h, tail_ratio, tail_slope)) / (2*h);
    BOOST_CHECK_CLOSE(result.v[2], fd_gamma, 0.01); // 0.01% tolerance
#endif
}

// ============================================================================
// Tests against scipy-generated Voigt PDF values
// ============================================================================

BOOST_AUTO_TEST_CASE(TestAgainstScipy_GaussianDominated) {
    using namespace ScipyVoigtTestData;
    
    for (size_t i = 0; i < Gaussian_dominated_npoints; ++i) {
        double x = Gaussian_dominated_x_values[i];
        double scipy_pdf = Gaussian_dominated_pdf_values[i];
        double computed_pdf = voigt_pdf<double>(x, Gaussian_dominated_mean, 
                                                Gaussian_dominated_sigma, 
                                                Gaussian_dominated_gamma);
        
        // For Gaussian-dominated case, expect very close agreement
        double rel_error = std::abs(computed_pdf - scipy_pdf) / scipy_pdf;
        BOOST_CHECK_MESSAGE(rel_error < 1e-3 || std::abs(computed_pdf - scipy_pdf) < 1e-6,
            "PDF mismatch at x=" << x << ": scipy=" << scipy_pdf 
            << ", computed=" << computed_pdf << ", rel_error=" << (rel_error * 100) << "%");
    }
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_LorentzianDominated) {
    using namespace ScipyVoigtTestData;
    
    for (size_t i = 0; i < Lorentzian_dominated_npoints; ++i) {
        double x = Lorentzian_dominated_x_values[i];
        double scipy_pdf = Lorentzian_dominated_pdf_values[i];
        double computed_pdf = voigt_pdf<double>(x, Lorentzian_dominated_mean, 
                                                Lorentzian_dominated_sigma, 
                                                Lorentzian_dominated_gamma);
        
        // For Lorentzian-dominated case, allow slightly more tolerance
        double rel_error = std::abs(computed_pdf - scipy_pdf) / scipy_pdf;
        BOOST_CHECK_MESSAGE(rel_error < 5e-3 || std::abs(computed_pdf - scipy_pdf) < 1e-5,
            "PDF mismatch at x=" << x << ": scipy=" << scipy_pdf 
            << ", computed=" << computed_pdf << ", rel_error=" << (rel_error * 100) << "%");
    }
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_Balanced) {
    using namespace ScipyVoigtTestData;
    
    for (size_t i = 0; i < Balanced_npoints; ++i) {
        double x = Balanced_x_values[i];
        double scipy_pdf = Balanced_pdf_values[i];
        double computed_pdf = voigt_pdf<double>(x, Balanced_mean, 
                                                Balanced_sigma, 
                                                Balanced_gamma);
        
        double rel_error = std::abs(computed_pdf - scipy_pdf) / scipy_pdf;
        BOOST_CHECK_MESSAGE(rel_error < 1e-3 || std::abs(computed_pdf - scipy_pdf) < 1e-6,
            "PDF mismatch at x=" << x << ": scipy=" << scipy_pdf 
            << ", computed=" << computed_pdf << ", rel_error=" << (rel_error * 100) << "%");
    }
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_GammaSpectroscopy) {
    using namespace ScipyVoigtTestData;
    
    for (size_t i = 0; i < Gamma_spectroscopy_npoints; ++i) {
        double x = Gamma_spectroscopy_x_values[i];
        double scipy_pdf = Gamma_spectroscopy_pdf_values[i];
        double computed_pdf = voigt_pdf<double>(x, Gamma_spectroscopy_mean, 
                                                Gamma_spectroscopy_sigma, 
                                                Gamma_spectroscopy_gamma);
        
        double rel_error = std::abs(computed_pdf - scipy_pdf) / scipy_pdf;
        BOOST_CHECK_MESSAGE(rel_error < 1e-3 || std::abs(computed_pdf - scipy_pdf) < 1e-6,
            "PDF mismatch at x=" << x << ": scipy=" << scipy_pdf 
            << ", computed=" << computed_pdf << ", rel_error=" << (rel_error * 100) << "%");
    }
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_SmallSigmaLargeGamma) {
    using namespace ScipyVoigtTestData;
    
    for (size_t i = 0; i < Small_sigma_large_gamma_npoints; ++i) {
        double x = Small_sigma_large_gamma_x_values[i];
        double scipy_pdf = Small_sigma_large_gamma_pdf_values[i];
        double computed_pdf = voigt_pdf<double>(x, Small_sigma_large_gamma_mean, 
                                                Small_sigma_large_gamma_sigma, 
                                                Small_sigma_large_gamma_gamma);
        
        // For extreme cases, allow more tolerance
        double rel_error = std::abs(computed_pdf - scipy_pdf) / scipy_pdf;
        BOOST_CHECK_MESSAGE(rel_error < 5e-3 || std::abs(computed_pdf - scipy_pdf) < 1e-5,
            "PDF mismatch at x=" << x << ": scipy=" << scipy_pdf 
            << ", computed=" << computed_pdf << ", rel_error=" << (rel_error * 100) << "%");
    }
}

BOOST_AUTO_TEST_CASE(TestAgainstScipy_LargeSigmaSmallGamma) {
    using namespace ScipyVoigtTestData;
    
    for (size_t i = 0; i < Large_sigma_small_gamma_npoints; ++i) {
        double x = Large_sigma_small_gamma_x_values[i];
        double scipy_pdf = Large_sigma_small_gamma_pdf_values[i];
        double computed_pdf = voigt_pdf<double>(x, Large_sigma_small_gamma_mean, 
                                                Large_sigma_small_gamma_sigma, 
                                                Large_sigma_small_gamma_gamma);
        
        double rel_error = std::abs(computed_pdf - scipy_pdf) / scipy_pdf;
        BOOST_CHECK_MESSAGE(rel_error < 1e-3 || std::abs(computed_pdf - scipy_pdf) < 1e-6,
            "PDF mismatch at x=" << x << ": scipy=" << scipy_pdf 
            << ", computed=" << computed_pdf << ", rel_error=" << (rel_error * 100) << "%");
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace
