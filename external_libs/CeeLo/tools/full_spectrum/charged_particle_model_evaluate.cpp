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

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ceelo::full_spectrum::SpectralComponent;

const char* component_name(SpectralComponent component, bool positron)
{
    switch (component) {
    case SpectralComponent::SourceFluorescenceKAlpha:
        return positron ? "positron_fluorescence_k_alpha" : "electron_fluorescence_k_alpha";
    case SpectralComponent::SourceFluorescenceKBeta:
        return positron ? "positron_fluorescence_k_beta" : "electron_fluorescence_k_beta";
    case SpectralComponent::SourceFluorescenceL:
        return positron ? "positron_fluorescence_l" : "electron_fluorescence_l";
    case SpectralComponent::PrimaryElectronBrems:
        return "electron_brems";
    case SpectralComponent::PrimaryPositronBrems:
        return "positron_brems";
    case SpectralComponent::PrimaryPositronAnnihilation:
        return "positron_511_atom";
    case SpectralComponent::PrimaryPositronAnnihilationContinuum:
        return "positron_annihilation_continuum";
    case SpectralComponent::PrimaryPositronAnnihilationInFlight:
        return "positron_annihilation_in_flight_direct";
    case SpectralComponent::
            PrimaryPositronAnnihilationInFlightContinuum:
        return "positron_annihilation_in_flight_continuum";
    default:
        return "unsupported";
    }
}

struct Config {
    std::string payload;
    std::string output = "charged_particle_model.csv";
    double energy_keV = 2290.0;
    double effective_Z = 92.0;
    double areal_density_g_cm2 = 1.0;
    int bins = 1024;
    double maximum_energy_keV = 3000.0;
    bool positron = false;
    bool fluorescence_bases = false;
};

Config parse(int argc, char** argv)
{
    Config config;
    auto value = [&](int& index) {
        if (++index >= argc)
            throw std::invalid_argument("missing option value");
        return std::string(argv[index]);
    };
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--payload")
            config.payload = value(index);
        else if (option == "--output")
            config.output = value(index);
        else if (option == "--energy-kev")
            config.energy_keV = std::stod(value(index));
        else if (option == "--effective-z")
            config.effective_Z = std::stod(value(index));
        else if (option == "--areal-density")
            config.areal_density_g_cm2 = std::stod(value(index));
        else if (option == "--bins")
            config.bins = std::stoi(value(index));
        else if (option == "--max-energy-kev")
            config.maximum_energy_keV = std::stod(value(index));
        else if (option == "--positron")
            config.positron = true;
        else if (option == "--fluorescence-bases")
            config.fluorescence_bases = true;
        else
            throw std::invalid_argument("unknown option: " + option);
    }
    if (config.fluorescence_bases)
        return config;
    if (config.payload.empty() || !(config.energy_keV > 0.0)
        || !(config.effective_Z >= 1.0 && config.effective_Z <= 92.0)
        || config.areal_density_g_cm2 < 0.0 || config.bins < 1
        || !(config.maximum_energy_keV > 0.0)) {
        throw std::invalid_argument("invalid model-evaluation configuration");
    }
    return config;
}

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
        const Config config = parse(argc, argv);
        if (config.fluorescence_bases) {
            // The fitter's fluorescence bases, from the evaluator's own
            // function: stdin rows "KE Z areal-density f_H positron(0/1)";
            // per row and family, a (impact ionization and absorbed brems) and
            // b (absorbed annihilation photons per unit containment).
            const ceelo::full_spectrum::SmoothElementCrossSection cross_sections(0.42);
            std::cout << std::setprecision(17)
                      << "energy_keV,Z,areal_density_g_cm2,hydrogen_mass_fraction,"
                         "positron,family,a,b\n";
            static constexpr const char* kFamilies[] = {"k_alpha", "k_beta", "l"};
            double energy = 0.0, z = 0.0, thickness = 0.0, fraction = 0.0;
            int positron = 0;
            while (std::cin >> energy >> z >> thickness >> fraction >> positron) {
                const auto bases = ceelo::full_spectrum::detail::charged_fluorescence_bases(
                    energy, thickness, z, fraction, positron != 0, cross_sections);
                for (int family = 0; family < 3; ++family)
                    std::cout << energy << ',' << z << ',' << thickness << ',' << fraction
                              << ',' << positron << ',' << kFamilies[family] << ','
                              << bases.a[family] << ',' << bases.b[family] << '\n';
            }
            return 0;
        }
        const auto model =
            ceelo::full_spectrum::ChargedParticleSourceFluxModel::
                deserialize(read_text(config.payload));
        std::vector<double> edges(
            static_cast<std::size_t>(config.bins + 1));
        for (int index = 0; index <= config.bins; ++index) {
            edges[static_cast<std::size_t>(index)] =
                config.maximum_energy_keV * index / config.bins;
        }
        const auto measure = model.evaluate(
            config.energy_keV, config.areal_density_g_cm2,
            config.effective_Z, config.positron, edges);

        std::ofstream output(config.output);
        if (!output)
            throw std::runtime_error("cannot open " + config.output);
        output << std::setprecision(17)
            << "campaign,source,primary_energy_keV,Z,"
               "areal_density_g_cm2,particle,component,bin_lo_keV,"
               "bin_hi_keV,sum_w,sum_w2,normalization_events,seed,physics\n";
        const char* source = config.positron ? "e+" : "e-";
        for (const auto& atom : measure.atoms) {
            output << "charged_particle_functional," << source << ','
                   << config.energy_keV << ',' << config.effective_Z << ','
                   << config.areal_density_g_cm2 << ",photon,"
                   << component_name(atom.component, config.positron) << ','
                   << atom.energy_keV << ',' << atom.energy_keV << ','
                   << atom.yield << ",0,1,0,"
                   << "positive-log-energy-bspline-rbf\n";
        }
        for (const auto& bin : measure.continuous) {
            output << "charged_particle_functional," << source << ','
                   << config.energy_keV << ',' << config.effective_Z << ','
                   << config.areal_density_g_cm2 << ",photon,"
                   << component_name(bin.component, config.positron) << ','
                   << bin.lower_keV << ',' << bin.upper_keV << ','
                   << bin.integral << ",0,1,0,"
                   << "positive-log-energy-bspline-rbf\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_charged_particle_model_evaluate: "
                  << error.what() << '\n';
        return 2;
    }
}
