#include "3rd_party/Faddeeva_MIT/Faddeeva.hh"
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include "../Faddeeva.hpp"

// ============================================================================
// Exact Voigt profile functions using MIT Faddeeva implementation
// ============================================================================

double voigt_pdf_exact(double x, double mean, double sigma, double gamma) {
    // Exact Voigt PDF: V(x) = Re[w(z)] / (σ * sqrt(π))
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

double voigt_cdf_exact(double x, double mean, double sigma, double gamma) {
    // Exact Voigt CDF: F_V(x) = 0.5 + 0.5 * Re[erf(z)]
    // where z = (x - μ + i*γ) / (σ * sqrt(2))
    const double sqrt2 = 1.41421356237309504880;
    
    if (gamma <= 1e-10 * sigma) {
        // Degenerate to pure Gaussian CDF
        double diff = (x - mean) / (sigma * sqrt2);
        return 0.5 * (1.0 + FaddeevaT::erf_real(diff));
    }
    
    double sigma_sqrt2 = sigma * sqrt2;
    std::complex<double> z((x - mean) / sigma_sqrt2, gamma / sigma_sqrt2);
    std::complex<double> erf_z = Faddeeva::erf(z);
    
    return 0.5 + 0.5 * erf_z.real();
}

// ============================================================================
// Exact exponentially-modified Gaussian (low-energy tail) functions
// ============================================================================

double gaussexp_pdf_exact(double x, double mean, double sigma, double tau) {
    // Exact EMG PDF for low-energy tail
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

double gaussexp_cdf_exact(double x, double mean, double sigma, double tau) {
    // Exact EMG CDF for low-energy tail
    // F(x) = 1 - Φ(a) + exp(-λ(μ - x) + (λσ)^2 / 2) * Φ(a - λσ)
    // where a = (μ - x)/σ and λ = 1/τ
    const double sqrt2 = 1.41421356237309504880;
    
    if (tau <= 0) {
        // Degenerate to pure Gaussian CDF
        double diff = (x - mean) / (sigma * sqrt2);
        return 0.5 * (1.0 + FaddeevaT::erf_real(diff));
    }
    
    const double lambda = 1.0 / tau;
    const double lambda_sigma = lambda * sigma;
    
    double a = (mean - x) / sigma;
    double phi_a = 0.5 * (1.0 + FaddeevaT::erf_real(a / sqrt2));
    double phi_shift = 0.5 * (1.0 + FaddeevaT::erf_real((a - lambda_sigma) / sqrt2));
    double exp_term = std::exp(-lambda * (mean - x) + 0.5 * lambda_sigma * lambda_sigma);
    
    double result = 1.0 - phi_a + exp_term * phi_shift;
    
    // Guard against tiny numerical drift
    if (result < 0.0) result = 0.0;
    if (result > 1.0) result = 1.0;
    
    return result;
}

// ============================================================================
// VoigtExpTail API functions (exact versions)
// ============================================================================

double voigt_exp_indefinite_exact(double x, double mean, double sigma_gauss,
                                   double gamma_lor, double tail_ratio, double tail_slope) {
    double voigt_cdf_val = voigt_cdf_exact(x, mean, sigma_gauss, gamma_lor);
    double gaussexp_cdf_val = gaussexp_cdf_exact(x, mean, sigma_gauss, tail_slope);
    return (1.0 - tail_ratio) * voigt_cdf_val + tail_ratio * gaussexp_cdf_val;
}

void voigt_exp_integral_exact(double peak_mean, double sigma_gauss,
                              double peak_amplitude, double gamma_lor,
                              double tail_ratio, double tail_slope,
                              const float * const energies, double *channels,
                              size_t nchannel) {
    // Evaluate CDF at all bin edges
    std::vector<double> cdf_values(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        double x = static_cast<double>(energies[i]);
        cdf_values[i] = voigt_exp_indefinite_exact(x, peak_mean, sigma_gauss, 
                                                    gamma_lor, tail_ratio, tail_slope);
    }
    
    // Compute bin contents as differences
    for (size_t i = 0; i < nchannel; ++i) {
        double bin_content = (cdf_values[i + 1] - cdf_values[i]) * peak_amplitude;
        if (bin_content < 0.0) bin_content = 0.0;  // Ensure non-negative
        channels[i] += bin_content;
    }
}

// ============================================================================
// Test distribution generation
// ============================================================================

struct TestDistribution {
    std::string name;
    double peak_mean;
    double sigma_gauss;
    double gamma_lor;
    double tail_ratio;
    double tail_slope;
    double peak_amplitude;
    double energy_min;
    double energy_max;
    double bin_width;
};

void generate_test_distribution(const TestDistribution& test, const std::string& output_name) {
    std::cout << "\n=== Generating " << test.name << " ===" << std::endl;
    std::cout << "  Peak mean: " << test.peak_mean << " keV" << std::endl;
    std::cout << "  Gaussian width: " << test.sigma_gauss << " keV" << std::endl;
    std::cout << "  Lorentzian width: " << test.gamma_lor << " keV" << std::endl;
    std::cout << "  Tail ratio: " << test.tail_ratio << std::endl;
    std::cout << "  Tail slope: " << test.tail_slope << " keV" << std::endl;
    std::cout << "  Energy range: " << test.energy_min << " - " << test.energy_max << " keV" << std::endl;
    std::cout << "  Bin width: " << test.bin_width << " keV" << std::endl;
    
    size_t nchannel = static_cast<size_t>((test.energy_max - test.energy_min) / test.bin_width);
    std::cout << "  Number of channels: " << nchannel << std::endl;
    
    // Create energy array (nchannel+1 entries for bin edges)
    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = test.energy_min + i * test.bin_width;
    }
    
    // Create channel array
    std::vector<double> channels(nchannel, 0.0);
    
    // Fill channels using exact implementation
    voigt_exp_integral_exact(test.peak_mean, test.sigma_gauss, test.peak_amplitude,
                             test.gamma_lor, test.tail_ratio, test.tail_slope,
                             energies.data(), channels.data(), nchannel);
    
    // Verify normalization
    double total_counts = 0.0;
    for (size_t i = 0; i < nchannel; ++i) {
        total_counts += channels[i];
    }
    
    std::cout << "  Total counts: " << total_counts << " (expected: " << test.peak_amplitude << ")" << std::endl;
    
    // Write to CSV
    std::ofstream csv_file(output_name + ".csv");
    csv_file << std::fixed << std::setprecision(6);
    csv_file << "energy (keV),counts\n";
    for (size_t i = 0; i < nchannel; ++i) {
        csv_file << energies[i] << "," << channels[i] * 100.0 << "\n";
    }
    csv_file.close();
    std::cout << "  Written to " << output_name << ".csv" << std::endl;
    
    // Generate C++ code for hard-coding into main.cpp
    std::ofstream cpp_file(output_name + "_test_data.cpp");
    cpp_file << "// Exact test data for: " << test.name << "\n";
    cpp_file << "// Generated using MIT Faddeeva implementation\n";
    cpp_file << "// Parameters: mean=" << test.peak_mean << ", sigma=" << test.sigma_gauss 
             << ", gamma=" << test.gamma_lor << ", tail_ratio=" << test.tail_ratio 
             << ", tail_slope=" << test.tail_slope << "\n";
    cpp_file << "// Energy range: " << test.energy_min << " - " << test.energy_max 
             << " keV, bin_width=" << test.bin_width << " keV\n\n";
    cpp_file << "const double " << output_name << "_energy_min = " << test.energy_min << ";\n";
    cpp_file << "const double " << output_name << "_energy_max = " << test.energy_max << ";\n";
    cpp_file << "const double " << output_name << "_bin_width = " << test.bin_width << ";\n";
    cpp_file << "const size_t " << output_name << "_nchannel = " << nchannel << ";\n\n";
    cpp_file << "const double " << output_name << "_exact_counts[" << nchannel << "] = {\n";
    cpp_file << std::scientific << std::setprecision(15);
    for (size_t i = 0; i < nchannel; ++i) {
        cpp_file << "  " << channels[i];
        if (i < nchannel - 1) cpp_file << ",";
        if ((i + 1) % 5 == 0) cpp_file << "\n";
    }
    cpp_file << "\n};\n";
    cpp_file.close();
    std::cout << "  Generated C++ test data in " << output_name << "_test_data.cpp" << std::endl;
}

// This file provides exact reference functions for unit testing
// The main function is not used - this is a unit test support file
