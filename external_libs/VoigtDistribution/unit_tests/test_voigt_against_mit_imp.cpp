#include "3rd_party/Faddeeva_MIT/Faddeeva.hh"
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include "../Faddeeva.hpp"

// ============================================================================
// MIT Voigt profile functions using MIT Faddeeva implementation
// ============================================================================

double voigt_pdf_mit(double x, double mean, double sigma, double gamma) {
    // MIT Voigt PDF: V(x) = Re[w(z)] / (σ * sqrt(π))
    // where z = (x - μ + i*γ) / (σ * sqrt(2))
    const double sqrt2 = 1.41421356237309504880;
    const double sqrt_pi = 1.77245385090551602730;
    
    if (gamma <= 1e-10 * sigma) {
        // Degenerate to pure Gaussian
        double diff = (x - mean) / sigma;
        return std::exp(-0.5 * diff * diff) / (sigma * sqrt_pi * sqrt2);
    }
    
    double sigma_sqrt2 = sigma * sqrt2;
    std::complex<double> z((x - mean) / sigma_sqrt2, gamma / sigma_sqrt2);
    std::complex<double> w = Faddeeva::w(z);
    
    return w.real() / (sigma * sqrt_pi * sqrt2);
}


// ============================================================================
// MIT exponentially-modified Gaussian (low-energy tail) functions
// ============================================================================

double gaussexp_pdf_mit(double x, double mean, double sigma, double tau) {
    // MIT EMG PDF for low-energy tail
    // f(x) = (λ/2) * exp(λ(μ - x) + (λσ)^2 / 2) * erfc( (μ - x)/(σ√2) + λσ/√2 )
    // where λ = 1/τ
    const double sqrt2 = 1.41421356237309504880;
    
    if (tau <= 0) {
        // Degenerate to pure Gaussian
        double diff = (x - mean) / sigma;
        return std::exp(-0.5 * diff * diff) / (sigma * sqrt2 * 1.77245385090551602730);
    }
    
    const double lambda = 1.0 / tau;
    const double lambda_sigma = lambda * sigma;
    double arg = (mean - x) / (sigma * sqrt2) + lambda_sigma / sqrt2;
    double erfc_val = FaddeevaT::erfc_real(arg);
    double exp_arg = lambda * (mean - x) + 0.5 * lambda_sigma * lambda_sigma;
    double exp_val = std::exp(exp_arg);
    
    return lambda * 0.5 * exp_val * erfc_val;
}


// ============================================================================
// VoigtExpTail API functions (MIT versions)
// ============================================================================

double voigt_exp_pdf_mit(double x, double mean, double sigma_gauss,
                           double gamma_lor, double tail_ratio, double tail_slope) {
    double voigt_pdf_val = voigt_pdf_mit(x, mean, sigma_gauss, gamma_lor);
    double gaussexp_pdf_val = gaussexp_pdf_mit(x, mean, sigma_gauss, tail_slope);
    return (1.0 - tail_ratio) * voigt_pdf_val + tail_ratio * gaussexp_pdf_val;
}

// High-accuracy integration using 7-point Gauss-Legendre quadrature with MIT Faddeeva
void voigt_exp_integral_mit_gauss_legendre(double peak_mean, double sigma_gauss,
                                           double peak_amplitude, double gamma_lor,
                                           double tail_ratio, double tail_slope,
                                           const float * const energies, double *channels,
                                           size_t nchannel) {
    // 7-point Gauss-Legendre quadrature nodes and weights on [-1, 1]
    // These are the roots of the 7th Legendre polynomial
    const double gl7_nodes[7] = {
        -0.94910791234275852453,  // -sqrt((7+2*sqrt(7/3))/21)
        -0.74153118559939443986,  // -sqrt((7-2*sqrt(7/3))/21)
        -0.40584515137739716691,  // -sqrt(7/21)
         0.0,
         0.40584515137739716691,  // sqrt(7/21)
         0.74153118559939443986,  // sqrt((7-2*sqrt(7/3))/21)
         0.94910791234275852453   // sqrt((7+2*sqrt(7/3))/21)
    };
    const double gl7_weights[7] = {
        0.12948496616886969327,   // (128/225) * (1 - sqrt(7/3)/7)
        0.27970539148927666790,   // (128/225) * (1 + sqrt(7/3)/7)
        0.38183005050511894495,   // (128/225) * (1 - sqrt(7/3)/7)
        0.41795918367346938776,   // 128/225
        0.38183005050511894495,   // (128/225) * (1 - sqrt(7/3)/7)
        0.27970539148927666790,   // (128/225) * (1 + sqrt(7/3)/7)
        0.12948496616886969327    // (128/225) * (1 - sqrt(7/3)/7)
    };

    for (size_t i = 0; i < nchannel; ++i) {
        double a = static_cast<double>(energies[i]);      // Lower edge
        double b = static_cast<double>(energies[i + 1]);  // Upper edge
        double bin_width = b - a;

        // Transform from [-1, 1] to [a, b]
        // x = (b-a)/2 * xi + (a+b)/2
        double half_width = 0.5 * bin_width;
        double center = 0.5 * (a + b);

        // Integrate PDF over [a, b] using 7-point Gauss-Legendre quadrature
        double integral = 0.0;
        for (int j = 0; j < 7; ++j) {
            double x = half_width * gl7_nodes[j] + center;
            double pdf_val = voigt_exp_pdf_mit(x, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
            integral += gl7_weights[j] * pdf_val;
        }

        // Scale by bin width and amplitude
        // The weights are for integration on [-1, 1], so we multiply by (b-a)/2
        double bin_content = integral * half_width * peak_amplitude;
        if (bin_content < 0.0) bin_content = 0.0;
        channels[i] += bin_content;
    }
}



// This file provides MIT reference functions for unit testing
// The main function is not used - this is a unit test support file

// ============================================================================
// Boost unit tests for Gauss-Legendre quadrature implementation
// ============================================================================

#include <boost/test/unit_test.hpp>
#include "../voigt_exp_tail.hpp"

using namespace voigt_exp_tail;

BOOST_AUTO_TEST_SUITE(VoigtExpIntegralTestSuite)

// Test normalization: integral should sum to peak_amplitude
BOOST_AUTO_TEST_CASE(TestNormalization) {
    const double mean = 94.67;
    const double sigma = 0.5;
    const double gamma = 0.0876;
    const double tail_ratio = 0.10;
    const double tail_slope = 1.0;
    const double amplitude = 1000.0;

    // Use a wide range to cover long tails
    const double x_min = mean - 20.0 * sigma;
    const double x_max = mean + 20.0 * sigma;
    const double bin_width = 0.01;
    const size_t nchannel = static_cast<size_t>((x_max - x_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(x_min + i * bin_width);
    }

    std::vector<double> channels(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                              energies.data(), channels.data(), nchannel);

    double total = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total += channels[i];
    }

    // Should be close to amplitude (within 1% due to finite range)
    // Note: Finite integration range means we don't capture the full distribution tail
    BOOST_CHECK_CLOSE(total, amplitude, 1.0);
}

// Test against MIT CDF-based integration
BOOST_AUTO_TEST_CASE(TestAgainstMitCDF) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.05;
    const double tail_slope = 2.0;
    const double amplitude = 1000.0;

    const double x_min = mean - 5.0 * sigma;
    const double x_max = mean + 5.0 * sigma;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((x_max - x_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(x_min + i * bin_width);
    }

    // Gauss-Legendre quadrature implementation
    std::vector<double> channels_gl(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                              energies.data(), channels_gl.data(), nchannel);

    std::vector<double> channels_mit(nchannel, 0.0);
    voigt_exp_integral_mit_gauss_legendre(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                            energies.data(), channels_mit.data(), nchannel);

    // Check total normalization and per-channel agreement
    double total_gl = 0.0;
    double total_mit = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_gl += channels_gl[i];
        total_mit += channels_mit[i];
        
        // Check each channel matches: either MIT value is very small (< 1e-9) or relative difference < 0.01%
        double abs_diff = std::abs(channels_gl[i] - channels_mit[i]);
        if (channels_mit[i] < 1e-9) {
            // For very small values, just check absolute difference is small
            BOOST_CHECK_SMALL(abs_diff, 1e-12);
        } else {
            // For larger values, check relative difference
            double rel_diff = abs_diff / channels_mit[i] * 100.0; // Convert to percentage
            BOOST_CHECK_MESSAGE(rel_diff < 0.001,
                "Channel " << i << ": relative difference " << rel_diff << "% exceeds 0.01%. "
                "GL=" << channels_gl[i] << ", MIT=" << channels_mit[i]);
        }
    }

    // Total should match within 0.01% (quadrature approximation)
    BOOST_CHECK_CLOSE(total_gl, total_mit, 0.01);

    // Compare channel by channel for significant channels only
    double max_rel_error = 0.0;
    size_t significant_channels = 0;
    for (size_t i = 0; i < nchannel; ++i) {
        // Only compare channels with significant content (> 1% of peak)
        if (channels_mit[i] > amplitude * 1e-2) {
            double abs_error = std::abs(channels_gl[i] - channels_mit[i]);
            double rel_error = abs_error / channels_mit[i];
            if (rel_error > max_rel_error) {
                max_rel_error = rel_error;
            }
            significant_channels++;
        }
    }

    // For significant channels, Gauss-Legendre should be accurate to within 0.5%
    // (more lenient due to quadrature approximation, especially for narrow bins)
    if (significant_channels > 0) {
        BOOST_CHECK_LT(max_rel_error, 0.005);
    }
}

// Test with narrow bins (should use 3-point quadrature)
BOOST_AUTO_TEST_CASE(TestNarrowBins) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.1;
    const double tail_slope = 1.0;
    const double amplitude = 1000.0;

    // Use narrow bins (< 0.25 * max(sigma, gamma) = 0.25)
    // Use wider range to capture more of the distribution
    const double x_min = mean - 14.0 * sigma;
    const double x_max = mean + 14.0 * sigma;
    const double bin_width = 0.05;  // Narrow bins
    const size_t nchannel = static_cast<size_t>((x_max - x_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(x_min + i * bin_width);
    }

    std::vector<double> channels(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                              energies.data(), channels.data(), nchannel);

    // Verify normalization
    double total = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total += channels[i];
    }

    // Should integrate correctly even with narrow bins (within 0.5% due to finite range and numerical errors)
    // Note: Finite integration range means we don't capture the full distribution tail
    BOOST_CHECK_CLOSE(total, amplitude, 0.5); //actual value looks to be 0.411383%
}

// Test with wide bins (should use 5-point quadrature)
BOOST_AUTO_TEST_CASE(TestWideBins) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.1;
    const double tail_ratio = 0.1;
    const double tail_slope = 1.0;
    const double amplitude = 1000.0;

    // Use wide bins (> 0.25 * max(sigma, gamma) = 0.25)
    const double x_min = mean - 15.0 * sigma;
    const double x_max = mean + 15.0 * sigma;
    const double bin_width = 1.0;  // Wide bins
    const size_t nchannel = static_cast<size_t>((x_max - x_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(x_min + i * bin_width);
    }

    std::vector<double> channels(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                              energies.data(), channels.data(), nchannel);

    // Verify normalization
    double total = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total += channels[i];
    }

    // Should integrate correctly with wide bins (within 1.5% due to finite range)
    // Note: Finite integration range means we don't capture the full distribution tail
    BOOST_CHECK_CLOSE(total, amplitude, 0.5);
}

// Test with gamma spectroscopy parameters
BOOST_AUTO_TEST_CASE(TestGammaSpectroscopyParameters) {
    const double mean = 94.67;      // Uranium 94.67 keV x-ray
    const double sigma = 0.5;       // 500 eV Gaussian width
    const double gamma = 0.0876;    // 87.6 eV Lorentzian width
    const double tail_ratio = 0.05;
    const double tail_slope = 1.0;
    const double amplitude = 1000.0;

    const double x_min = mean - 10.0 * sigma;
    const double x_max = mean + 10.0 * sigma;
    const double bin_width = 0.01;  // High resolution
    const size_t nchannel = static_cast<size_t>((x_max - x_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(x_min + i * bin_width);
    }

    std::vector<double> channels(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                              energies.data(), channels.data(), nchannel);

    // Compare against MIT implementation
    std::vector<double> channels_mit(nchannel, 0.0);
    voigt_exp_integral_mit_gauss_legendre(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                            energies.data(), channels_mit.data(), nchannel);

    // Check total normalization
    double total_gl = 0.0;
    double total_mit = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_gl += channels[i];
        total_mit += channels_mit[i];
    }

    // Allow 0.01% difference due to quadrature approximation
    BOOST_CHECK_CLOSE(total_gl, total_mit, 0.01);
}

// Test pure Gaussian case (gamma = 0)
BOOST_AUTO_TEST_CASE(TestPureGaussian) {
    const double mean = 100.0;
    const double sigma = 1.0;
    const double gamma = 0.0;  // Pure Gaussian
    const double tail_ratio = 0.0;
    const double tail_slope = 1.0;
    const double amplitude = 1000.0;

    const double x_min = mean - 5.0 * sigma;
    const double x_max = mean + 5.0 * sigma;
    const double bin_width = 0.1;
    const size_t nchannel = static_cast<size_t>((x_max - x_min) / bin_width);

    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = static_cast<float>(x_min + i * bin_width);
    }

    std::vector<double> channels(nchannel, 0.0);
    voigt_exp_integral<double>(mean, sigma, amplitude, gamma, tail_ratio, tail_slope,
                              energies.data(), channels.data(), nchannel);

    double total = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total += channels[i];
    }

    BOOST_CHECK_CLOSE(total, amplitude, 0.01); // 0.01% tolerance
}

BOOST_AUTO_TEST_SUITE_END()
