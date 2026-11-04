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

#include "cross_sections/CrossSectionData.h"
#include "full_spectrum/ArrivalResponse.h"
#include "geometry/Geometry.h"
#include "geometry/GeometryDescriptor.h"
#include "materials/Material.h"
#include "transport/PhotonTransport.h"
#include "transport/TransportUtils.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

using ceelo::InteractionType;
using ceelo::TransportResult;
using ceelo::full_spectrum::campaign::Options;
using ceelo::full_spectrum::campaign::Tally;

enum Component : std::size_t {
    ZeroDeposit,
    FullEnergy,
    PairFullEnergy,
    SingleEscape,
    DoubleEscape,
    Annihilation511,
    KAlphaEscape,
    KBetaEscape,
    KAlphaAbsorbed,
    KBetaAbsorbed,
    PhotoelectricEscape,
    SingleCompton,
    MultipleCompton,
    PairNoCompton,
    PairSingleCompton,
    PairMultipleCompton,
    NumComponents
};

const char* component_name(Component component) {
    static constexpr const char* names[] = {
        "zero_deposit", "full_energy", "pair_full_energy",
        "single_escape", "double_escape", "annihilation_511",
        "k_alpha_escape", "k_beta_escape",
        "k_alpha_absorbed", "k_beta_absorbed",
        "photoelectric_escape", "single_compton", "multiple_compton",
        "pair_no_compton", "pair_single_compton",
        "pair_multiple_compton"};
    return names[static_cast<std::size_t>(component)];
}

// Probability-weighted mean energy of an element's K-alpha and K-beta lines
// (families from the lines' filling subshells). Returns false when the
// element has no tabulated K fluorescence or lacks either family.
bool k_line_family_means(int Z, double* alpha_mean_keV,
                         double* beta_mean_keV) {
    const ceelo::FluorescenceData* fluorescence =
        ceelo::CrossSectionData::instance().fluorescence(Z);
    if (!fluorescence || fluorescence->num_lines == 0)
        return false;
    const ceelo::KLineFamilyMeans means =
        ceelo::k_line_family_means(*fluorescence);
    if (!(means.alpha_probability > 0.0) || !(means.beta_probability > 0.0))
        return false;
    *alpha_mean_keV = means.alpha_keV;
    *beta_mean_keV = means.beta_keV;
    return true;
}

Component classify(const TransportResult& result, double energy_keV,
                   int crystal_Z) {
    const double deposit = result.energy_deposited_scoring;
    if (deposit <= 0.0)
        return ZeroDeposit;
    int compton = 0;
    int pair_compton = 0;
    bool pair = false;
    bool photoelectric = false;
    bool after_pair = false;
    for (const auto& interaction : result.interactions) {
        if (interaction.type == InteractionType::PairProduction) {
            pair = true;
            after_pair = true;
        } else if (interaction.type == InteractionType::Compton) {
            ++compton;
            if (after_pair)
                ++pair_compton;
        }
        if (interaction.type == InteractionType::Photoelectric)
            photoelectric = true;
    }
    if (pair && std::abs(deposit - energy_keV) < 0.05)
        return PairFullEnergy;
    if (!pair && std::abs(deposit - energy_keV) < 1e-6)
        return FullEnergy;
    if (pair && energy_keV > 1021.9979
        && std::abs(deposit - (energy_keV - 510.99895)) < 0.05)
        return SingleEscape;
    if (pair && energy_keV > 1021.9979
        && std::abs(deposit - (energy_keV - 1021.9979)) < 0.05)
        return DoubleEscape;
    // This is a physical pair-history atom rather than a generic energy
    // window. It primarily represents pair production in passive detector
    // material followed by full absorption of one at-rest annihilation
    // photon in the active crystal. Keeping it out of the pair continuum
    // prevents a smooth basis from having to synthesize a delta function.
    if (pair && std::abs(deposit - 510.99895) < 0.05)
        return Annihilation511;
    if (pair) {
        if (pair_compton == 0)
            return PairNoCompton;
        return pair_compton == 1
            ? PairSingleCompton : PairMultipleCompton;
    }
    // A crystal K X-ray that escapes deposits Ein - E_line — a discrete
    // line, not part of the smooth photoelectric-escape or Compton continua
    // (same sharp-features-are-atoms rule as the pair escape peaks and the
    // 511 atom above). That holds however the vacancy was made (photoelectric
    // absorption, or a Compton collision with spectrum fidelity on) and
    // however many Compton scatters came first, since everything but the
    // X-ray is absorbed. The deposit window plus the escaped flag reject
    // dead-layer/attenuator fluorescence that scattered into the crystal.
    for (const auto& record : result.fluorescence_records) {
        // Only the scoring crystal's own K x-ray makes a K-escape line
        // at the energies this campaign writes; a dead-layer or endcap
        // fluorescence record would be labelled with the wrong element's
        // line energy, and (being gated on the crystal's edge) could be
        // dropped from the output entirely.
        if (record.element_Z != crystal_Z)
            continue;
        if (!record.was_k_shell || !record.fluorescence_emitted)
            continue;
        double alpha_mean = 0.0, beta_mean = 0.0;
        if (!k_line_family_means(record.element_Z, &alpha_mean, &beta_mean))
            break;
        const bool is_alpha = ceelo::k_line_family(record.filling_subshell)
            == ceelo::KLineFamily::Alpha;
        if (record.fluorescence_escaped) {
            if (std::abs(deposit
                         - (energy_keV
                            - record.fluorescence_energy_keV))
                < 0.05)
                return is_alpha ? KAlphaEscape : KBetaEscape;
            continue;
        }
        // The complementary population: the vacancy is made in a dead
        // layer, its electrons and Auger cascade are lost there, and
        // only the K x-ray reaches the active volume, where
        // it is absorbed. The deposit is then the LINE energy itself,
        // independent of Ein -- a second pair of discrete lines that
        // carries most of the photoelectric-escape continuum below
        // ~100 keV (67% of it at 42.5 keV).
        if (std::abs(deposit - record.fluorescence_energy_keV) < 0.05)
            return is_alpha ? KAlphaAbsorbed : KBetaAbsorbed;
    }
    if (photoelectric && compton == 0)
        return PhotoelectricEscape;
    return compton <= 1 ? SingleCompton : MultipleCompton;
}

// Optional finer partition of the pair-history continua, written as extra
// rows alongside the standard labels, which are unaffected.
// The post-pair Compton count mixes physically distinct spectral sources: an
// annihilation photon that escapes cleanly, one that escapes after scattering,
// and one that is absorbed all produce different deposited-energy structure
// but can share a Compton count. These labels separate them by fate instead.
enum RefinedPairComponent : std::size_t {
    PairEscape0,
    PairEscape1Unscattered,
    PairEscape1Scattered,
    PairEscape2Unscattered,
    PairEscape2Scattered,
    PairInFlight,
    NumRefinedPairComponents
};

const char* refined_pair_name(RefinedPairComponent component) {
    static constexpr const char* names[] = {
        "pair_escape0", "pair_escape1_unscattered",
        "pair_escape1_scattered", "pair_escape2_unscattered",
        "pair_escape2_scattered", "pair_in_flight"};
    return names[static_cast<std::size_t>(component)];
}

RefinedPairComponent classify_refined_pair(const TransportResult& result) {
    const auto& history = result.pair_history;
    // In-flight annihilation produces a continuum photon pair rather than two
    // 511 keV photons, so its escape structure is not comparable.
    if (history.annihilated_in_flight)
        return PairInFlight;
    // A positron that left the geometry takes both photons with it: neither
    // deposits, as when both escape unscattered.
    if (history.positron_escaped)
        return PairEscape2Unscattered;
    const bool scattered = history.any_scattered_escape();
    switch (history.escape_count()) {
    case 0:
        return PairEscape0;
    case 1:
        return scattered ? PairEscape1Scattered : PairEscape1Unscattered;
    default:
        return scattered ? PairEscape2Scattered : PairEscape2Unscattered;
    }
}

// The two built-in detectors, as descriptors: built through
// GeometryDescriptor::build_geometry they trace exactly as the hand-built
// geometries they replace (same calls, same order, same materials).
ceelo::GeometryDescriptor builtin_detector(const std::string& name) {
    ceelo::GeometryDescriptor gd;
    gd.crystal_material_index = 0;
    if (name == "hpge_coax") {
        gd.set_dimensions(ceelo::CylinderDims{4.0, 8.0});
        gd.materials = {ceelo::MaterialSpec::from(ceelo::make_HPGe()),
                        ceelo::MaterialSpec::from(ceelo::make_Aluminum())};
        gd.bore = ceelo::BoreHoleConfig{0.5, 5.5};
        gd.dead_layer = ceelo::DeadLayerConfig{0.05, 0.08, 0.08};
        gd.layers = {{1, 0.10, 0.10, -0.10, 8.10}};
    } else if (name == "nai3x3") {
        gd.set_dimensions(ceelo::CylinderDims{3.81, 7.62});
        gd.materials = {ceelo::MaterialSpec::from(ceelo::make_NaI()),
                        ceelo::MaterialSpec::from(ceelo::make_Aluminum())};
        gd.layers = {{1, 0.05, 0.05, -0.05, 7.67}};
    } else {
        throw std::invalid_argument(
            "detector must be nai3x3 or hpge_coax, or come from --geometry FILE");
    }
    return gd;
}

ceelo::GeometryDescriptor read_geometry_file(const std::string& path) {
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot read geometry file " + path);
    std::ostringstream text;
    text << in.rdbuf();
    return ceelo::GeometryDescriptor::from_text(text.str());
}

// The element whose K lines the escape and absorbed atoms are labelled with:
// the crystal's heaviest constituent (Ge in HPGe, I in NaI).
int k_line_element(const ceelo::GeometryDescriptor& gd) {
    int z = 0;
    for (const ceelo::MaterialComponent& c :
         gd.materials.at(static_cast<std::size_t>(gd.crystal_material_index))
             .composition)
        z = std::max(z, static_cast<int>(c.Z));
    return z;
}

Eigen::Vector3d perpendicular(const Eigen::Vector3d& direction) {
    const Eigen::Vector3d trial =
        std::abs(direction.z()) < 0.8
        ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitX();
    return direction.cross(trial).normalized();
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options =
            ceelo::full_spectrum::campaign::parse_options(
                argc, argv, true);

        const ceelo::GeometryDescriptor descriptor =
            options.geometry_file.empty()
                ? builtin_detector(options.detector)
                : read_geometry_file(options.geometry_file);
        if (options.detector.empty()
            || options.detector.find_first_of(", \t\n") != std::string::npos)
            throw std::invalid_argument(
                "--detector must be a label without commas or whitespace");
        if (options.print_geometry) {
            std::cout << descriptor.to_text();
            return 0;
        }
        std::vector<std::unique_ptr<ceelo::Material>> materials;
        const ceelo::Geometry geometry = descriptor.build_geometry(materials);
        const int crystal_Z = k_line_element(descriptor);

        // Detector-side bremsstrahlung and fluorescence are opt-in (the bake
        // enables both). They matter above ~2 MeV: with brems off a
        // contained electron deposits all of its kinetic energy and therefore
        // lands exactly on a full-energy or escape line, which inflates those
        // lines at the expense of the continuum. Enable them to compare
        // against a Geant4 run, which always has both.

        ceelo::TransportConfig config;
        config.record_interactions = true;
        config.enable_fluorescence = options.detector_fluorescence;
        config.enable_pair_production = true;
        config.enable_rayleigh = true;
        config.disable_brems = !options.detector_brems;
        if (options.detector_spectrum_fidelity) {
            config.enable_compton_fluorescence = true;
            config.electron_transport = ceelo::ElectronTransport::Detailed;
            config.enable_impact_ionization = true;
        }

        const double sine =
            std::sqrt(std::max(0.0, 1.0 - options.cos_theta
                                       * options.cos_theta));
        const Eigen::Vector3d direction(
            sine, 0.0, options.cos_theta);
        const Eigen::Vector3d e1 = perpendicular(direction);
        const Eigen::Vector3d e2 = direction.cross(e1).normalized();
        const auto z_extent = geometry.outer_z_extent();
        const Eigen::Vector3d center(
            0.0, 0.0, 0.5 * (z_extent.first + z_extent.second));
        const double half_z = 0.5 * (z_extent.second - z_extent.first);
        const double bound = 1.01 * std::sqrt(
            geometry.outer_bounding_radius()
                * geometry.outer_bounding_radius()
            + half_z * half_z);

        std::array<Tally, NumComponents> areas;
        std::array<std::vector<Tally>, NumComponents> spectra;
        for (auto& spectrum : spectra)
            spectrum.resize(static_cast<std::size_t>(options.bins));
        std::array<Tally, NumRefinedPairComponents> refined_areas;
        std::array<std::vector<Tally>, NumRefinedPairComponents>
            refined_spectra;
        if (options.refined_pair_labels) {
            for (auto& spectrum : refined_spectra)
                spectrum.resize(static_cast<std::size_t>(options.bins));
        }
        std::vector<Tally> coarse(32);
        std::mt19937_64 rng(options.seed);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        std::uint64_t events = 0;
        constexpr std::uint64_t batch = 5000;

        // An arrival class replaces the fixed beam: every photon draws its own
        // direction from the band (full_spectrum::ArrivalBeam), so the cell is
        // the area-weighted average of the fixed-angle response over it.
        const bool arrival = !options.incidence.empty();
        const ceelo::full_spectrum::ArrivalBand band = arrival
            ? ceelo::full_spectrum::arrival_band(options.incidence)
            : ceelo::full_spectrum::ArrivalBand{};
        ceelo::full_spectrum::ArrivalBeam arrival_beam(geometry, band);

        while (events < options.max_events) {
            const std::uint64_t stop =
                std::min(options.max_events, events + batch);
            while (events < stop) {
                Eigen::Vector3d origin;
                Eigen::Vector3d photon_direction = direction;
                if (arrival) {
                    arrival_beam.sample(rng, origin, photon_direction);
                } else {
                    do {
                        const double r = bound * std::sqrt(uniform(rng));
                        const double phi =
                            2.0 * 3.14159265358979323846 * uniform(rng);
                        origin = center - direction * (2.0 * bound)
                            + e1 * (r * std::cos(phi))
                            + e2 * (r * std::sin(phi));
                    } while (
                        !geometry.ray_hits_outer_boundary(origin, direction));
                }

                const TransportResult result = ceelo::transport_photon(
                    origin, photon_direction, options.energy_keV,
                    geometry, config, rng);
                const Component component =
                    classify(result, options.energy_keV, crystal_Z);
                const double weight = result.weight;
                areas[component].add(weight);
                const double deposit = result.energy_deposited_scoring;
                const bool atom = component == ZeroDeposit
                    || component == FullEnergy
                    || component == PairFullEnergy
                    || component == SingleEscape
                    || component == DoubleEscape
                    || component == Annihilation511
                    || component == KAlphaEscape
                    || component == KBetaEscape
                    || component == KAlphaAbsorbed
                    || component == KBetaAbsorbed;
                if (!atom && deposit >= 0.0
                    && deposit < options.energy_keV) {
                    const int bin = std::min(
                        options.bins - 1,
                        static_cast<int>(
                            deposit / options.energy_keV * options.bins));
                    spectra[component][static_cast<std::size_t>(bin)]
                        .add(weight);
                    if (options.refined_pair_labels
                        && result.pair_history.occurred) {
                        const RefinedPairComponent refined =
                            classify_refined_pair(result);
                        refined_areas[refined].add(weight);
                        refined_spectra[refined]
                            [static_cast<std::size_t>(bin)].add(weight);
                    }
                    const int coarse_bin = std::min(
                        31, static_cast<int>(
                            deposit / options.energy_keV * 32));
                    coarse[static_cast<std::size_t>(coarse_bin)].add(weight);
                }
                ++events;
            }
            std::vector<Tally> area_vector(areas.begin(), areas.end());
            if (events >= options.min_events
                && ceelo::full_spectrum::campaign::precision_reached(
                    events, area_vector, coarse))
                break;
        }

        std::cout << std::setprecision(12);
        // Arrival cells say so in cell_kind (so a fixed-angle reader rejects
        // them), leave cos_theta empty, and append the band and the boundary's
        // mean projected area over it (bounding disk x acceptance, +- sigma).
        std::cout
            << "cell_kind,energy_keV,Z,tau,areal_density_g_cm2,cos_theta,"
               "detector,component,bin_lo_keV,bin_hi_keV,sum_w,sum_w2,"
               "events,seed,analytic_direct"
            << (arrival ? ",incidence,mu_lo,mu_hi,mean_projected_area_cm2,"
                          "mean_projected_area_sigma_cm2"
                        : "")
            << '\n';
        std::ostringstream suffix_text;
        suffix_text << std::setprecision(12);
        if (arrival) {
            suffix_text << ',' << options.incidence << ',' << band.mu_lo << ','
                << band.mu_hi << ',' << arrival_beam.mean_projected_area_cm2()
                << ',' << arrival_beam.mean_projected_area_sigma_cm2();
        }
        const std::string row_suffix = suffix_text.str();
        auto write_row = [&](const char* name, double lo, double hi,
                             const Tally& tally) {
            std::cout << (arrival ? "detector_arrival," : "detector,")
                << options.energy_keV << ",,,,";
            if (!arrival)
                std::cout << options.cos_theta;
            std::cout << ',' << options.detector << ',' << name
                << ',' << lo << ',' << hi << ','
                << tally.sum_w << ',' << tally.sum_w2 << ','
                << events << ',' << options.seed << ',' << row_suffix << '\n';
        };
        // Crystal K-escape atom energies: the probability-weighted family
        // means of the crystal element's K lines. Rows are emitted only
        // above the K edge.
        double k_alpha_mean = 0.0, k_beta_mean = 0.0;
        const bool k_escape_possible =
            k_line_family_means(crystal_Z, &k_alpha_mean, &k_beta_mean)
            && options.energy_keV
                > ceelo::CrossSectionData::instance()
                      .fluorescence(crystal_Z)->k_edge_keV;
        for (std::size_t c = 0; c < NumComponents; ++c) {
            const auto component = static_cast<Component>(c);
            const bool line = component == ZeroDeposit
                || component == FullEnergy
                || component == PairFullEnergy
                || component == SingleEscape
                || component == DoubleEscape
                || component == Annihilation511
                || component == KAlphaEscape
                || component == KBetaEscape
                || component == KAlphaAbsorbed
                || component == KBetaAbsorbed;
            if (line) {
                if (component == SingleEscape
                    && options.energy_keV <= 1021.9979)
                    continue;
                if (component == DoubleEscape
                    && options.energy_keV <= 1021.9979)
                    continue;
                if ((component == KAlphaEscape
                     || component == KBetaEscape
                     || component == KAlphaAbsorbed
                     || component == KBetaAbsorbed)
                    && !k_escape_possible)
                    continue;
                double energy = 0.0;
                if (component == FullEnergy
                    || component == PairFullEnergy)
                    energy = options.energy_keV;
                else if (component == SingleEscape)
                    energy = options.energy_keV - 510.99895;
                else if (component == DoubleEscape)
                    energy = options.energy_keV - 1021.9979;
                else if (component == Annihilation511)
                    energy = 510.99895;
                else if (component == KAlphaEscape)
                    energy = options.energy_keV - k_alpha_mean;
                else if (component == KBetaEscape)
                    energy = options.energy_keV - k_beta_mean;
                else if (component == KAlphaAbsorbed)
                    energy = k_alpha_mean;
                else if (component == KBetaAbsorbed)
                    energy = k_beta_mean;
                write_row(component_name(component), energy, energy,
                          areas[c]);
                continue;
            }
            for (int b = 0; b < options.bins; ++b) {
                write_row(component_name(component),
                          options.energy_keV * b / options.bins,
                          options.energy_keV * (b + 1) / options.bins,
                          spectra[c][static_cast<std::size_t>(b)]);
            }
        }

        // Refined pair rows are additive: the twelve standard labels above are
        // unchanged, so a consumer that does not know these names can ignore
        // them.
        if (options.refined_pair_labels) {
            for (std::size_t c = 0; c < NumRefinedPairComponents; ++c) {
                const auto refined =
                    static_cast<RefinedPairComponent>(c);
                for (int b = 0; b < options.bins; ++b) {
                    write_row(refined_pair_name(refined),
                              options.energy_keV * b / options.bins,
                              options.energy_keV * (b + 1) / options.bins,
                              refined_spectra[c][static_cast<std::size_t>(b)]);
                }
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_detector_response_campaign: "
                  << error.what() << '\n';
        return 2;
    }
}
