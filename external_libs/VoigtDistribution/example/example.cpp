#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "../voigt_exp_tail.hpp"

/**
 * Example program demonstrating VoigtDistribution library usage
 *
 * This program generates various gamma spectroscopy spectra using the
 * VoigtExpTail distribution and saves them to CSV files for analysis.
 */

void generate_spectrum_to_csv(const std::string& filename,
                               const double peak_mean, const double sigma_gauss,
                               const double gamma_lor, const double tail_ratio,
                               const double tail_slope, const double peak_amplitude,
                               const double energy_min, const double energy_max,
                               const double bin_width) {
    const size_t nchannel = static_cast<size_t>((energy_max - energy_min) / bin_width);

    // Create energy array (nchannel+1 entries for bin edges)
    std::vector<float> energies(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        energies[i] = energy_min + i * bin_width;
    }

    // Create channel array
    std::vector<double> channels(nchannel, 0.0);

    // Fill channels using the VoigtExpTail distribution
    // Use CDF-based integration (accurate but numerically problematic)
    voigt_exp_tail::voigt_exp_integral<double>(peak_mean, sigma_gauss, peak_amplitude, gamma_lor,
                       tail_ratio, tail_slope, energies.data(), channels.data(), nchannel);


    // Write spectrum to CSV file
    std::ofstream csv_file(filename);
    csv_file << std::fixed << std::setprecision(6);
    csv_file << "energy (keV),counts\n";
    for (size_t i = 0; i < nchannel; ++i) {
        csv_file << energies[i] << "," << channels[i] << "\n";
    }
    csv_file.close();

    // Report statistics
    double total_counts = 0.0;
    double max_counts = 0.0;
    size_t peak_channel = 0;

    for (size_t i = 0; i < nchannel; ++i) {
        total_counts += channels[i];
        if (channels[i] > max_counts) {
            max_counts = channels[i];
            peak_channel = i;
        }
    }

    double peak_energy = energies[peak_channel] + 0.5 * bin_width;

    std::cout << "  Generated " << filename << ":" << std::endl;
    std::cout << "    Total counts: " << total_counts << std::endl;
    std::cout << "    Peak energy: " << peak_energy << " keV" << std::endl;
    std::cout << "    Peak counts: " << max_counts << std::endl;
    std::cout << "    Channels: " << nchannel << std::endl;
}

int main() {
    std::cout << "VoigtDistribution Example - Gamma Spectroscopy Spectra" << std::endl;
    std::cout << "======================================================" << std::endl;
    std::cout << std::endl;

    // Example 1: Uranium 94.67 keV x-ray peak (realistic gamma spectroscopy parameters)
    std::cout << "Example 1: Uranium 94.67 keV x-ray peak" << std::endl;
    std::cout << "Parameters: HPGe detector (σ=0.5 keV), natural linewidth (γ=87.6 eV), 5% tailing" << std::endl;
    generate_spectrum_to_csv("uranium_xray.csv",
                           94.67,    // peak mean (keV)
                           0.5,      // Gaussian sigma (keV) - typical HPGe resolution
                           0.0876,   // Lorentzian gamma (keV) - natural linewidth
                           0.05,     // tail ratio (5% tailing)
                           1.0,      // tail slope
                           1000.0,   // peak amplitude (counts)
                           85.0,     // energy min (keV)
                           105.0,    // energy max (keV)
                           0.01);    // bin width (keV) - high resolution

    std::cout << std::endl;

    // Example 2: Pure Gaussian (no Lorentzian, no tail)
    std::cout << "Example 2: Pure Gaussian peak" << std::endl;
    std::cout << "Parameters: Only detector resolution, no natural broadening or tailing" << std::endl;
    generate_spectrum_to_csv("pure_gaussian.csv",
                           100.0,    // peak mean (keV)
                           1.0,      // Gaussian sigma (keV)
                           0.0,      // no Lorentzian
                           0.0,      // no tail
                           1.0,      // tail slope (unused)
                           1000.0,   // peak amplitude
                           90.0,     // energy min
                           110.0,    // energy max
                           0.1);     // bin width

    std::cout << std::endl;

    // Example 3: Gaussian with exponential tail (incomplete charge collection)
    std::cout << "Example 3: Gaussian with exponential low-energy tail" << std::endl;
    std::cout << "Parameters: Detector resolution + incomplete charge collection tailing" << std::endl;
    generate_spectrum_to_csv("gaussian_with_tail.csv",
                           100.0,    // peak mean
                           1.0,      // Gaussian sigma
                           0.0,      // no Lorentzian
                           0.15,     // 15% tail ratio
                           1.0,      // tail slope
                           1000.0,   // peak amplitude
                           90.0,     // energy min
                           110.0,    // energy max
                           0.1);     // bin width

    std::cout << std::endl;

    // Example 4: Voigt profile with large Lorentzian component
    std::cout << "Example 4: Voigt profile (significant natural broadening)" << std::endl;
    std::cout << "Parameters: Both detector resolution and large natural linewidth" << std::endl;
    generate_spectrum_to_csv("voigt_large_lorentzian.csv",
                           100.0,    // peak mean
                           1.0,      // Gaussian sigma
                           0.5,      // large Lorentzian gamma
                           0.0,      // no tail
                           1.0,      // tail slope (unused)
                           1000.0,   // peak amplitude
                           90.0,     // energy min
                           110.0,    // energy max
                           0.1);     // bin width

    std::cout << std::endl;
    std::cout << "======================================================" << std::endl;
    std::cout << "Generated CSV files:" << std::endl;
    std::cout << "  - uranium_xray.csv        (realistic gamma spectroscopy example)" << std::endl;
    std::cout << "  - pure_gaussian.csv       (detector resolution only)" << std::endl;
    std::cout << "  - gaussian_with_tail.csv  (with incomplete charge collection)" << std::endl;
    std::cout << "  - voigt_large_lorentzian.csv (with natural linewidth)" << std::endl;
    std::cout << std::endl;
    std::cout << "These CSV files can be plotted to visualize the different peak shapes." << std::endl;

    return 0;
}
