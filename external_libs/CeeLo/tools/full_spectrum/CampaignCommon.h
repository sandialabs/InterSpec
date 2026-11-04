#pragma once
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

#include <Eigen/Core>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace ceelo::full_spectrum::campaign {

struct Options {
    double energy_keV = 661.657;
    int Z = 82;
    double tau = 1.0;
    double areal_density_g_cm2 = -1.0;
    // Shield hydrogen mass fraction; the remainder is the elemental Z. Zero
    // is a single-element shield.
    double hydrogen_mass_fraction = 0.0;
    double cos_theta = 1.0;
    // An arrival class (full_spectrum::arrival_band: iso, front, back, side,
    // mu-band:LO,HI) in place of the fixed beam angle. Empty: --cos-theta.
    std::string incidence;
    std::uint64_t seed = 0x46554c4c53504543ULL;
    std::uint64_t min_events = 20000;
    std::uint64_t max_events = 250000;
    int bins = 256;
    // A built-in detector (nai3x3, hpge_coax), or, with geometry_file, the
    // label written in the CSV for the detector that file describes.
    std::string detector = "nai3x3";
    // A GeometryDescriptor in its "CEELO_GEOMETRY 1" text form.
    std::string geometry_file;
    // Write the detector's geometry text to stdout instead of running.
    bool print_geometry = false;
    bool interaction_orders = false;
    // Opt-in finer pair-history labels, emitted as extra CSV rows.
    bool refined_pair_labels = false;
    // Detector-side secondary physics, opt-in. The bake enables both, as does
    // any comparison against Geant4 (which always has both).
    bool detector_brems = false;
    bool detector_fluorescence = false;
    // Spectrum-fidelity detector physics (EfficiencyCalculator::
    // enable_spectrum_fidelity): Compton K-vacancy x-rays and the Detailed
    // electron walk. The bake enables it.
    bool detector_spectrum_fidelity = false;
    bool fine_interaction_orders = false;
    bool source_rayleigh = true;
    bool source_compton_doppler = true;
    std::string multiple_surrogate_file;
    std::string secondary_surrogate_file;
    std::string photonuclear_xs_file;
    std::string photonuclear_surrogate_file;
    // Rayleigh-only ("elastic") line constants. Negative means "leave the
    // payload default alone", so an unpassed flag never silently overwrites
    // what a loaded payload already carries.
    double elastic_escape_scale = -1.0;
    double elastic_compton_damping = -1.0;
    std::string payload_file;
    std::string write_payload_file;
};

inline Options parse_options(int argc, char** argv, bool detector_mode) {
    Options options;
    bool cos_theta_given = false;
    auto value = [&](int& i) -> std::string {
        if (++i >= argc)
            throw std::invalid_argument("Missing value after command-line option");
        return argv[i];
    };
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--refined-pair-labels")
            options.refined_pair_labels = true;
        else if (arg == "--detector-brems")
            options.detector_brems = true;
        else if (arg == "--detector-fluorescence")
            options.detector_fluorescence = true;
        else if (arg == "--detector-spectrum-fidelity")
            options.detector_spectrum_fidelity = true;
        else if (arg == "--energy-kev")
            options.energy_keV = std::stod(value(i));
        else if (arg == "--z")
            options.Z = std::stoi(value(i));
        else if (arg == "--tau")
            options.tau = std::stod(value(i));
        else if (arg == "--areal-density")
            options.areal_density_g_cm2 = std::stod(value(i));
        else if (arg == "--hydrogen-mass-fraction")
            options.hydrogen_mass_fraction = std::stod(value(i));
        else if (arg == "--cos-theta") {
            options.cos_theta = std::stod(value(i));
            cos_theta_given = true;
        }
        else if (arg == "--incidence")
            options.incidence = value(i);
        else if (arg == "--seed")
            options.seed = std::stoull(value(i));
        else if (arg == "--min-events")
            options.min_events = std::stoull(value(i));
        else if (arg == "--max-events")
            options.max_events = std::stoull(value(i));
        else if (arg == "--events")
            options.min_events = options.max_events = std::stoull(value(i));
        else if (arg == "--bins")
            options.bins = std::stoi(value(i));
        else if (arg == "--detector")
            options.detector = value(i);
        else if (arg == "--geometry")
            options.geometry_file = value(i);
        else if (arg == "--print-geometry")
            options.print_geometry = true;
        else if (arg == "--interaction-orders")
            options.interaction_orders = true;
        else if (arg == "--fine-interaction-orders") {
            options.interaction_orders = true;
            options.fine_interaction_orders = true;
        }
        else if (arg == "--source-rayleigh")
            options.source_rayleigh = true;
        else if (arg == "--no-source-rayleigh")
            options.source_rayleigh = false;
        else if (arg == "--source-compton-doppler")
            options.source_compton_doppler = true;
        else if (arg == "--no-source-compton-doppler")
            options.source_compton_doppler = false;
        else if (arg == "--multiple-surrogate")
            options.multiple_surrogate_file = value(i);
        else if (arg == "--secondary-surrogate")
            options.secondary_surrogate_file = value(i);
        else if (arg == "--photonuclear-xs")
            options.photonuclear_xs_file = value(i);
        else if (arg == "--photonuclear-surrogate")
            options.photonuclear_surrogate_file = value(i);
        else if (arg == "--elastic-escape-scale")
            options.elastic_escape_scale = std::stod(value(i));
        else if (arg == "--elastic-compton-damping")
            options.elastic_compton_damping = std::stod(value(i));
        else if (arg == "--payload")
            options.payload_file = value(i);
        else if (arg == "--write-payload")
            options.write_payload_file = value(i);
        else if (arg == "--help") {
            std::cout
                << "Cell driver options:\n"
                << "  --energy-kev E --seed N --events N --bins N\n";
            if (detector_mode)
                std::cout
                    << "  --cos-theta C --detector nai3x3|hpge_coax\n"
                    << "  --incidence iso|front|back|side|mu-band:LO,HI  "
                       "photons arriving\n"
                    << "      uniformly over a band of direction cosines "
                       "(area-weighted;\n"
                    << "      rows are cell_kind detector_arrival), instead "
                       "of --cos-theta\n"
                    << "  --geometry FILE --detector LABEL  any detector, from "
                       "a CEELO_GEOMETRY text file\n"
                    << "  --print-geometry  write the detector's geometry text "
                       "and exit\n";
            else
                std::cout
                    << "  --z Z (--tau T | --areal-density G_CM2)\n"
                    << "  --hydrogen-mass-fraction F  two-component shield: "
                       "F hydrogen + (1-F) element Z\n"
                    << "  --interaction-orders  split multiple-Compton truth\n"
                    << "  --fine-interaction-orders  additionally split "
                       "3, 4, and 5+ Compton truth\n"
                    << "  --source-rayleigh  enable coherent source scatter "
                       "(the default)\n"
                    << "  --no-source-rayleigh  disable coherent source "
                       "scatter for an ablation\n"
                    << "  --source-compton-doppler  enable bound-electron "
                       "Doppler broadening (the source-campaign default)\n"
                    << "  --no-source-compton-doppler  free-electron energy "
                       "kinematics ablation\n"
                    << "  --multiple-surrogate FILE  residual payload\n"
                    << "  --secondary-surrogate FILE  pair/brems residual "
                       "payload\n"
                    << "  --photonuclear-xs FILE  Geant4 photonuclear "
                       "mass-attenuation table\n"
                    << "  --elastic-escape-scale S       Rayleigh-only "
                       "line scale (default: payload's)\n"
                    << "  --elastic-compton-damping D    Compton damping of "
                       "that line (default: payload's)\n"
                    << "  --photonuclear-surrogate FILE  fitted (E,Z,tau) "
                       "photonuclear yield/shape nodes\n"
                    << "  --payload FILE  load a serialized "
                       "CEELO_SOURCE_ESCAPE payload instead of the "
                       "surrogate files (evaluator only)\n"
                    << "  --write-payload FILE  serialize the assembled "
                       "source-escape payload (evaluator only)\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }
    if (!(options.energy_keV >= 10.0 && options.energy_keV <= 12000.0))
        throw std::invalid_argument("energy must be in [10,12000] keV");
    if (options.bins < 8 || options.bins > 16384)
        throw std::invalid_argument("bins must be in [8,16384]");
    if (options.min_events == 0 || options.max_events < options.min_events)
        throw std::invalid_argument("invalid event limits");
    if (options.hydrogen_mass_fraction < 0.0
        || options.hydrogen_mass_fraction >= 1.0) {
        throw std::invalid_argument(
            "hydrogen mass fraction must be in [0,1)");
    }
    if (detector_mode) {
        if (options.cos_theta < -1.0 || options.cos_theta > 1.0)
            throw std::invalid_argument("cos-theta must be in [-1,1]");
        if (cos_theta_given && !options.incidence.empty())
            throw std::invalid_argument(
                "--incidence and --cos-theta are exclusive");
    } else {
        if (options.Z < 1 || options.Z > 92 || options.tau < 0.0
            || (options.areal_density_g_cm2 < 0.0
                && options.areal_density_g_cm2 != -1.0))
            throw std::invalid_argument(
                "Z must be 1..92 and shield amount nonnegative");
    }
    return options;
}

struct Tally {
    double sum_w = 0.0;
    double sum_w2 = 0.0;
    void add(double weight) {
        sum_w += weight;
        sum_w2 += weight * weight;
    }
    double relative_uncertainty() const {
        return sum_w > 0.0 ? std::sqrt(sum_w2) / sum_w : 0.0;
    }
};

inline bool precision_reached(
    std::uint64_t events,
    const std::vector<Tally>& component_areas,
    const std::vector<Tally>& coarse_continuum)
{
    if (events == 0)
        return false;
    for (const Tally& area : component_areas) {
        if (area.sum_w >= 100.0 * events / 100000.0
            && area.relative_uncertainty() > 0.003)
            return false;
    }
    for (const Tally& region : coarse_continuum) {
        if (region.sum_w >= 25.0
            && region.relative_uncertainty() > 0.02)
            return false;
    }
    return true;
}

inline Eigen::Vector3d isotropic_direction(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double z = 2.0 * uniform(rng) - 1.0;
    const double phi = 2.0 * 3.14159265358979323846 * uniform(rng);
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    return {r * std::cos(phi), r * std::sin(phi), z};
}

} // namespace ceelo::full_spectrum::campaign
