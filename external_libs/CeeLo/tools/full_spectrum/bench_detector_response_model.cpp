/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
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

#include "full_spectrum/FullSpectrumModel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string read_text(const std::string& path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open " + path);
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc < 2 || argc > 4) {
            throw std::invalid_argument(
                "usage: bench_detector_response_model PAYLOAD "
                "[iterations=100] [bins=4096]");
        }
        const int iterations = argc >= 3 ? std::stoi(argv[2]) : 100;
        const int bins = argc >= 4 ? std::stoi(argv[3]) : 4096;
        if (iterations < 1 || bins < 1)
            throw std::invalid_argument("iterations and bins must be positive");
        const auto model =
            ceelo::full_spectrum::DetectorPulseResponseModel::deserialize(
                read_text(argv[1]));
        std::vector<double> edges(static_cast<std::size_t>(bins + 1));
        constexpr double maximum_keV = 12000.0;
        for (int index = 0; index <= bins; ++index) {
            edges[static_cast<std::size_t>(index)] =
                maximum_keV * index / bins;
        }
        // Warm all code and coefficient pages before starting the CPU clock.
        (void)model.evaluate(1332.492, 0.35, edges);
        volatile double checksum = 0.0;
        const std::clock_t start = std::clock();
        for (int index = 0; index < iterations; ++index) {
            const double fraction =
                (index + 0.5) / static_cast<double>(iterations);
            const double energy =
                10.0 * std::exp(std::log(1200.0) * fraction);
            const double cosine =
                -1.0 + 2.0 * std::fmod(index * 0.61803398875, 1.0);
            const auto measure = model.evaluate(energy, cosine, edges);
            checksum = checksum + measure.integrated_yield();
        }
        const double cpu_seconds =
            static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;
        std::cout
            << "iterations,bins,cpu_seconds,cpu_ms_per_response,checksum\n"
            << iterations << ',' << bins << ',' << cpu_seconds << ','
            << 1000.0 * cpu_seconds / iterations << ','
            << checksum << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "bench_detector_response_model: "
                  << error.what() << '\n';
        return 2;
    }
}
