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

// Emit Geant4's natural-element photonuclear mass attenuation coefficient.
// This is a fast, deterministic feature-bank generator; no run manager or
// transport histories are required.

#include "cross_sections/CrossSectionData.h"

#include "G4Gamma.hh"
#include "G4GammaNuclearXS.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct Config {
    double minimum_energy_keV = 1000.0;
    double maximum_energy_keV = 12000.0;
    int energy_nodes = 96;
    int minimum_Z = 1;
    int maximum_Z = 92;
    std::string output;
};

Config parse(int argc, char** argv) {
    Config config;
    auto value = [&](int& index) -> std::string {
        if (++index >= argc)
            throw std::invalid_argument("missing option value");
        return argv[index];
    };
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--min-energy-kev")
            config.minimum_energy_keV = std::stod(value(index));
        else if (argument == "--max-energy-kev")
            config.maximum_energy_keV = std::stod(value(index));
        else if (argument == "--energy-nodes")
            config.energy_nodes = std::stoi(value(index));
        else if (argument == "--min-z")
            config.minimum_Z = std::stoi(value(index));
        else if (argument == "--max-z")
            config.maximum_Z = std::stoi(value(index));
        else if (argument == "--output")
            config.output = value(index);
        else
            throw std::invalid_argument(
                "unknown option: " + argument);
    }
    if (!(config.minimum_energy_keV > 0.0)
        || config.maximum_energy_keV < config.minimum_energy_keV
        || config.energy_nodes < 2
        || config.minimum_Z < 1
        || config.maximum_Z > 92
        || config.maximum_Z < config.minimum_Z) {
        throw std::invalid_argument("invalid grid");
    }
    return config;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse(argc, argv);
        G4GammaNuclearXS cross_section;
        cross_section.BuildPhysicsTable(*G4Gamma::GammaDefinition());
        constexpr double avogadro = 6.02214076e23;
        std::ofstream output_file;
        std::ostream* output = &std::cout;
        if (!config.output.empty()) {
            output_file.open(config.output);
            if (!output_file)
                throw std::runtime_error(
                    "cannot open " + config.output);
            output = &output_file;
        }
        *output << std::setprecision(12)
                << "energy_keV,Z,sigma_barn,mu_over_rho_cm2_g\n";
        for (int energy_index = 0;
             energy_index < config.energy_nodes; ++energy_index) {
            const double fraction =
                static_cast<double>(energy_index)
                / static_cast<double>(config.energy_nodes - 1);
            const double energy_keV =
                config.minimum_energy_keV
                * std::pow(
                    config.maximum_energy_keV
                        / config.minimum_energy_keV,
                    fraction);
            for (int Z = config.minimum_Z; Z <= config.maximum_Z; ++Z) {
                const double sigma_barn =
                    cross_section.ElementCrossSection(
                        energy_keV * keV, Z)
                    / barn;
                const double atomic_weight =
                    ceelo::CrossSectionData::instance().atomic_weight(Z);
                const double mass_attenuation =
                    sigma_barn * 1e-24 * avogadro / atomic_weight;
                *output << energy_keV << ',' << Z << ','
                        << sigma_barn << ',' << mass_attenuation
                        << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "g4_photonuclear_xs: "
                  << error.what() << '\n';
        return 2;
    }
}
