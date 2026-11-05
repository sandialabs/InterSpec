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

// Microbenchmark for ChargedParticleSourceFluxModel::evaluate.
// Times the double path only: the fitted-scalar count is fixed, so double
// timing is the comparison metric between candidate representations.

#include "full_spectrum/FullSpectrumModel.h"

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
                "usage: bench_charged_particle_source PAYLOAD "
                "[iterations=200] [bins=1024]");
        }
        const int iterations = argc >= 3 ? std::stoi(argv[2]) : 200;
        const int bins = argc >= 4 ? std::stoi(argv[3]) : 1024;
        if (iterations < 1 || bins < 1)
            throw std::invalid_argument(
                "iterations and bins must be positive");
        const auto model =
            ceelo::full_spectrum::ChargedParticleSourceFluxModel::
                deserialize(read_text(argv[1]));
        std::vector<double> edges(static_cast<std::size_t>(bins + 1));
        constexpr double maximum_keV = 3000.0;
        for (int index = 0; index <= bins; ++index) {
            edges[static_cast<std::size_t>(index)] =
                maximum_keV * index / bins;
        }
        for (bool positron : {false, true}) {
            (void)model.evaluate(2290.0, 1.0, 35.0, positron, edges);
            volatile double checksum = 0.0;
            const std::clock_t start = std::clock();
            for (int index = 0; index < iterations; ++index) {
                const double fraction =
                    (index + 0.5) / static_cast<double>(iterations);
                const double kinetic_energy =
                    50.0 * std::exp(std::log(200.0) * fraction);
                const double areal_density =
                    0.05 + 5.0 * std::fmod(index * 0.61803398875, 1.0);
                const double effective_Z =
                    4.0 + 84.0 * std::fmod(index * 0.3819660113, 1.0);
                const auto measure = model.evaluate(
                    kinetic_energy, areal_density, effective_Z,
                    positron, edges);
                checksum = checksum + measure.integrated_yield();
            }
            const std::clock_t stop = std::clock();
            const double microseconds_per_call =
                1.0e6 * static_cast<double>(stop - start)
                / CLOCKS_PER_SEC / iterations;
            std::cout << (positron ? "positron" : "electron")
                      << " calls=" << iterations
                      << " bins=" << bins
                      << " us_per_call=" << microseconds_per_call
                      << " checksum=" << checksum << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "bench_charged_particle_source: "
                  << error.what() << '\n';
        return 2;
    }
}
