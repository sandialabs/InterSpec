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

#include "CampaignCommon.h"
#include "ElectronCampaignOutput.h"

#include "cross_sections/CrossSectionData.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {

using PhotonOrigin =
    ceelo::SourceGeometry::SourceSecondaryPhoton::Origin;

/// The family label of an unscattered shield x-ray of element Z (exact line
/// energy; lines below 10 keV are outside the model), or nullptr.
const char* fluorescence_family(int Z, double energy_keV)
{
    if (energy_keV < 10.0) return nullptr;
    const auto& xs = ceelo::CrossSectionData::instance();
    if (const ceelo::FluorescenceData* k = xs.fluorescence(Z))
        for (int i = 0; i < k->num_lines; ++i)
            if (std::abs(k->line_energy_keV[i] - energy_keV) < 1e-6)
                return ceelo::k_line_family(k->line_filling_subshell[i])
                        == ceelo::KLineFamily::Alpha
                    ? "fluorescence_k_alpha" : "fluorescence_k_beta";
    if (const ceelo::LFluorescenceData* l = xs.l_fluorescence(Z))
        for (const auto& sub : l->sub)
            for (int i = 0; i < sub.num_lines; ++i)
                if (std::abs(sub.line_energy_keV[i] - energy_keV) < 1e-6)
                    return "fluorescence_l";
    return nullptr;
}

const char* escaped_photon_component(
    const ceelo::SourceGeometry::SourceSecondaryPhoton& photon,
    bool primary_is_positron, int Z)
{
    // Shield K and L x-rays (photoabsorption or impact ionization) that leave
    // unscattered are the fluorescence atoms; scattered ones stay with brems.
    if ((photon.origin == PhotonOrigin::Fluorescence
         || photon.origin == PhotonOrigin::ImpactIonization)
        && photon.compton_interactions == 0) {
        if (const char* family = fluorescence_family(Z, photon.energy_keV))
            return family;
    }
    if (primary_is_positron
        && photon.origin == PhotonOrigin::PrimaryPositronAnnihilation) {
        return photon.compton_interactions == 0
            && std::abs(photon.energy_keV - 510.99895) < 1e-6
            ? "positron_511_atom"
            : "positron_annihilation_continuum";
    }
    if (primary_is_positron
        && photon.origin
            == PhotonOrigin::PrimaryPositronAnnihilationInFlight) {
        return photon.compton_interactions == 0
            ? "positron_annihilation_in_flight_direct"
            : "positron_annihilation_in_flight_continuum";
    }
    if (photon.origin == PhotonOrigin::SecondaryElectronBrems)
        return "secondary_charged_brems";
    return primary_is_positron
        ? "positron_brems" : "electron_brems";
}

struct Config {
    double energy_keV = 2290.0;
    int Z = 92;
    double areal_density_g_cm2 = 1.0;
    std::uint64_t events = 100000;
    std::uint64_t seed = 0x454c454354524f4eULL;
    int bins = 512;
    double maximum_energy_keV = 3000.0;
    double minimum_photon_keV = 10.0;
    bool positron = false;
    std::string output = "ceelo_source_electron.csv";
};

Config parse(int argc, char** argv)
{
    Config config;
    auto value = [&](int& index) {
        if (++index >= argc)
            throw std::invalid_argument("missing option value");
        return std::string(argv[index]);
    };
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--energy-kev")
            config.energy_keV = std::stod(value(i));
        else if (option == "--z")
            config.Z = std::stoi(value(i));
        else if (option == "--areal-density")
            config.areal_density_g_cm2 = std::stod(value(i));
        else if (option == "--events")
            config.events = std::stoull(value(i));
        else if (option == "--seed")
            config.seed = std::stoull(value(i));
        else if (option == "--bins")
            config.bins = std::stoi(value(i));
        else if (option == "--max-energy-kev")
            config.maximum_energy_keV = std::stod(value(i));
        else if (option == "--min-photon-kev")
            config.minimum_photon_keV = std::stod(value(i));
        else if (option == "--positron")
            config.positron = true;
        else if (option == "--output")
            config.output = value(i);
        else
            throw std::invalid_argument("unknown option: " + option);
    }
    if (!(config.energy_keV > 0.0 && config.energy_keV <= 12000.0)
        || config.Z < 1 || config.Z > 92
        || config.areal_density_g_cm2 < 0.0
        || config.events == 0 || config.bins < 8
        || config.maximum_energy_keV <= 0.0
        || config.minimum_photon_keV < 0.0)
        throw std::invalid_argument("invalid campaign configuration");
    return config;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const Config config = parse(argc, argv);
        ceelo::Material material(
            "equivalent-element", 1.0,
            {{static_cast<std::uint8_t>(config.Z), 1.0}});
        ceelo::SourceGeometry source;
        source.configure_point(Eigen::Vector3d::Zero());
        if (config.areal_density_g_cm2 > 0.0)
            source.add_shield(&material, config.areal_density_g_cm2);
        source.set_source_brems(true);
        source.set_source_positron_transport(true);
        // Shield K and L x-rays, from the charged particle's bremsstrahlung
        // and annihilation photons absorbed and from impact ionization: the
        // fluorescence atoms.
        source.set_source_fluorescence(true);
        source.set_source_impact_ionization(true);
        source.set_rayleigh_enabled(true);
        source.set_record_truth_labels(true);

        ceelo::full_spectrum::electron_campaign::Spectra spectra(
            config.bins, config.maximum_energy_keV);
        std::mt19937_64 rng(config.seed);
        for (std::uint64_t event = 0; event < config.events; ++event) {
            const auto result = source.transport_source_electron(
                Eigen::Vector3d::Zero(),
                ceelo::full_spectrum::campaign::isotropic_direction(rng),
                config.energy_keV, rng, config.positron);
            spectra.tally(
                "diagnostic", "primary_path_g_cm2",
                0.0, result.traversed_gcm2);
            if (result.escaped) {
                spectra.tally(
                    config.positron ? "positron" : "electron",
                    "escaped_primary", result.escaped_electron.energy_keV, 1.0);
            }
            for (const auto& photon : result.created_brems_photons) {
                if (photon.energy_keV >= config.minimum_photon_keV) {
                    spectra.tally(
                        "created_photon",
                        photon.origin == PhotonOrigin::SecondaryElectronBrems
                            ? "secondary_charged_created"
                            : "primary_electron_created",
                        photon.energy_keV, 1.0);
                }
            }
            for (const auto& photon : result.photons) {
                if (photon.energy_keV >= config.minimum_photon_keV) {
                    spectra.tally(
                        "photon",
                        escaped_photon_component(
                            photon, config.positron, config.Z),
                        photon.energy_keV, 1.0);
                }
            }
        }
        spectra.write(
            config.output, "mono_electron",
            config.positron ? "e+" : "e-", config.energy_keV, config.Z,
            config.areal_density_g_cm2, config.events, config.seed,
            "CeeLo-CSDA-SeltzerBerger-MollerBhabha-KLfluorescence-impactionization");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_source_electron_campaign: "
                  << error.what() << '\n';
        return 2;
    }
}
