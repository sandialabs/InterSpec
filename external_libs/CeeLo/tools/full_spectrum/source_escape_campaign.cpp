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
#include "full_spectrum/FullSpectrumModel.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <utility>
#include <vector>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

using ceelo::SourceGeometry;
using ceelo::full_spectrum::campaign::Options;
using ceelo::full_spectrum::campaign::Tally;

enum Component : std::size_t {
    Direct,
    Elastic,
    SingleCompton,
    DoubleCompton,
    HigherCompton,
    SingleComptonRayleigh,
    DoubleComptonRayleigh,
    HigherComptonRayleigh,
    TripleCompton,
    QuadrupleCompton,
    FivePlusCompton,
    TripleComptonRayleigh,
    QuadrupleComptonRayleigh,
    FivePlusComptonRayleigh,
    Pair511Atom,
    PairContinuum,
    PairInFlight,
    PhotoelectronBrems,
    ComptonElectronBrems,
    PairElectronBrems,
    FluorescenceKAlpha,
    FluorescenceKBeta,
    FluorescenceL,
    FluorescenceContinuum,
    Absorbed,
    NumComponents
};

const char* component_name(Component component, bool interaction_orders) {
    switch (component) {
    case Direct:
        return "uncollided";
    case Elastic:
        return "elastic_line";
    case SingleCompton:
        return "single_compton";
    case DoubleCompton:
        return interaction_orders
            ? "double_compton" : "multiple_compton";
    case HigherCompton:
        return interaction_orders
            ? "higher_compton" : "multiple_compton";
    case SingleComptonRayleigh:
        return interaction_orders
            ? "single_compton_rayleigh" : "multiple_compton";
    case DoubleComptonRayleigh:
        return interaction_orders
            ? "double_compton_rayleigh" : "multiple_compton";
    case HigherComptonRayleigh:
        return interaction_orders
            ? "higher_compton_rayleigh" : "multiple_compton";
    case TripleCompton:
        return "triple_compton";
    case QuadrupleCompton:
        return "quadruple_compton";
    case FivePlusCompton:
        return "five_plus_compton";
    case TripleComptonRayleigh:
        return "triple_compton_rayleigh";
    case QuadrupleComptonRayleigh:
        return "quadruple_compton_rayleigh";
    case FivePlusComptonRayleigh:
        return "five_plus_compton_rayleigh";
    case Pair511Atom:
        return "pair_511";
    case PairContinuum:
        return "pair_annihilation_continuum";
    case PairInFlight:
        return "pair_annihilation_in_flight";
    case PhotoelectronBrems:
        return "photoelectron_brems";
    case ComptonElectronBrems:
        return "compton_electron_brems";
    case PairElectronBrems:
        return "pair_electron_brems";
    case FluorescenceKAlpha:
        return "fluorescence_k_alpha";
    case FluorescenceKBeta:
        return "fluorescence_k_beta";
    case FluorescenceL:
        return "fluorescence_l";
    case FluorescenceContinuum:
        return "fluorescence_continuum";
    case Absorbed:
        return "absorbed";
    case NumComponents:
        break;
    }
    return "unknown";
}

Component classify_primary(
    const SourceGeometry::SourceFullTransportResult& result,
    bool fine_interaction_orders)
{
    if (!result.survived)
        return Absorbed;
    if (!result.interacted)
        return Direct;
    if (result.truth.compton_interactions == 0
        && result.truth.rayleigh_interactions > 0)
        return Elastic;
    const auto compton = result.truth.compton_interactions;
    if (fine_interaction_orders && compton >= 3) {
        if (compton == 3)
            return TripleCompton;
        if (compton == 4)
            return QuadrupleCompton;
        return FivePlusCompton;
    }
    // Interaction order in this study means energy-changing Compton order.
    // Rayleigh remains active in transport, but does not promote a history
    // from single- to multiple-Compton or change its Compton-order label.
    if (compton == 1)
        return SingleCompton;
    if (compton == 2)
        return DoubleCompton;
    return HigherCompton;
}

/// The shield element's fluorescence lines (exact float energies, as the
/// transport emits them) and their families: K-alpha / K-beta from each K
/// line's filling subshell, L for the retained L lines.
struct LineFamilies {
    std::vector<std::pair<double, Component>> lines;
    std::array<double, 3> mean_keV{};   // probability-weighted, informational
    Component family_of(double energy_keV) const {
        for (const auto& [line, family] : lines)
            if (std::abs(line - energy_keV) < 1e-6) return family;
        return FluorescenceContinuum;
    }
};

LineFamilies shield_line_families(int Z)
{
    LineFamilies out;
    std::array<double, 3> sum{}, weight{};
    const auto& xs = ceelo::CrossSectionData::instance();
    if (const ceelo::FluorescenceData* k = xs.fluorescence(Z)) {
        for (int i = 0; i < k->num_lines; ++i) {
            const int f = static_cast<int>(ceelo::k_line_family(k->line_filling_subshell[i]));
            out.lines.emplace_back(k->line_energy_keV[i], f == 0 ? FluorescenceKAlpha
                                                                 : FluorescenceKBeta);
            sum[f] += k->line_energy_keV[i] * k->line_probability[i];
            weight[f] += k->line_probability[i];
        }
    }
    if (const ceelo::LFluorescenceData* l = xs.l_fluorescence(Z)) {
        for (const auto& sub : l->sub)
            for (int i = 0; i < sub.num_lines; ++i) {
                out.lines.emplace_back(sub.line_energy_keV[i], FluorescenceL);
                sum[2] += sub.line_energy_keV[i] * sub.fluorescence_yield * sub.line_probability[i];
                weight[2] += sub.fluorescence_yield * sub.line_probability[i];
            }
    }
    for (int f = 0; f < 3; ++f) out.mean_keV[f] = weight[f] > 0.0 ? sum[f] / weight[f] : 0.0;
    return out;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options =
            ceelo::full_spectrum::campaign::parse_options(
                argc, argv, false);
        // Two-component shield: hydrogen mass fraction + elemental remainder.
        // f_H = 0 constructs a single-element material.
        const double f_H = options.hydrogen_mass_fraction;
        std::vector<ceelo::MaterialComponent> composition;
        if (f_H > 0.0)
            composition.push_back({std::uint8_t{1}, f_H});
        composition.push_back(
            {static_cast<std::uint8_t>(options.Z), 1.0 - f_H});
        ceelo::Material element("equivalent-element", 1.0, composition);
        const auto incident_xs =
            element.macroscopic_xs(options.energy_keV * 1e-3);
        const double mu_mass = (
            incident_xs.mu_pe + incident_xs.mu_cs + incident_xs.mu_pp
            + (options.source_rayleigh ? incident_xs.mu_rs : 0.0))
            / element.density();
        const double areal_density =
            options.areal_density_g_cm2 >= 0.0
            ? options.areal_density_g_cm2
            : (options.tau == 0.0 ? 0.0 : options.tau / mu_mass);
        const double optical_depth = areal_density * mu_mass;

        SourceGeometry source;
        source.configure_point(Eigen::Vector3d::Zero());
        if (areal_density > 0.0)
            source.add_shield(&element, areal_density);
        // Gamma-produced source-electron bremsstrahlung is part of the source
        // escape kernel. "Gamma-only" excludes beta-decay primaries, not the
        // secondary electrons created by photon interactions.
        source.set_source_brems(true);
        source.set_source_electron_transport(false);
        source.set_source_positron_transport(true);
        // Shield K and L x-rays, from photoabsorption, Compton vacancies (with
        // Doppler sampling) and electron-impact ionization: the fluorescence
        // atoms. The rows record both switches; the fitter refuses a bank
        // without them.
        source.set_source_fluorescence(true);
        source.set_source_impact_ionization(true);
        const LineFamilies line_families = shield_line_families(options.Z);
        source.set_record_truth_labels(true);
        source.set_rayleigh_enabled(options.source_rayleigh);
        source.set_compton_doppler_enabled(
            options.source_compton_doppler);

        std::array<Tally, NumComponents> areas;
        std::array<std::vector<Tally>, NumComponents> spectra;
        for (auto& spectrum : spectra)
            spectrum.resize(static_cast<std::size_t>(options.bins));
        std::vector<Tally> coarse(32);
        std::mt19937_64 rng(options.seed);
        std::uint64_t events = 0;
        constexpr std::uint64_t batch = 5000;

        auto tally_photon = [&](double energy, double weight,
                                Component component, bool atom) {
            areas[component].add(weight);
            if (!atom && energy >= 0.0 && energy < options.energy_keV) {
                const int bin = std::min(
                    options.bins - 1,
                    static_cast<int>(energy / options.energy_keV * options.bins));
                spectra[component][static_cast<std::size_t>(bin)].add(weight);
                const int coarse_bin = std::min(
                    31, static_cast<int>(energy / options.energy_keV * 32));
                coarse[static_cast<std::size_t>(coarse_bin)].add(weight);
            }
        };

        while (events < options.max_events) {
            const std::uint64_t stop =
                std::min(options.max_events, events + batch);
            while (events < stop) {
                const auto result = source.transport_source_photon(
                    Eigen::Vector3d::Zero(),
                    ceelo::full_spectrum::campaign::isotropic_direction(rng),
                    options.energy_keV, rng);
                const Component primary_component = classify_primary(
                    result, options.fine_interaction_orders);
                if (primary_component == Absorbed) {
                    areas[Absorbed].add(result.bias_weight);
                } else {
                    const bool atom = primary_component == Direct
                        || primary_component == Elastic;
                    tally_photon(
                        result.energy_keV, result.bias_weight,
                        primary_component, atom);
                }
                for (const auto& secondary : result.secondaries) {
                    Component component = ComptonElectronBrems;
                    bool atom = false;
                    using Origin =
                        SourceGeometry::SourceSecondaryPhoton::Origin;
                    switch (secondary.origin) {
                    case Origin::PairAnnihilation:
                        atom = std::abs(
                            secondary.energy_keV - 510.99895) < 1e-3;
                        component = atom ? Pair511Atom : PairContinuum;
                        break;
                    case Origin::PairAnnihilationInFlight:
                        component = PairInFlight;
                        break;
                    case Origin::PhotoelectronBrems:
                        component = PhotoelectronBrems;
                        break;
                    case Origin::ComptonElectronBrems:
                        component = ComptonElectronBrems;
                        break;
                    case Origin::SecondaryElectronBrems:
                        // Delta-electron brems is a collision descendant of
                        // the tracked source electron.  Keep it in the broad
                        // recoil-electron continuum until the payload gains a
                        // dedicated component.
                        component = ComptonElectronBrems;
                        break;
                    case Origin::PairElectronBrems:
                        component = PairElectronBrems;
                        break;
                    case Origin::PrimaryElectronBrems:
                    case Origin::RadioactiveBetaBrems:
                    case Origin::ConversionElectronBrems:
                    case Origin::PrimaryPositronAnnihilation:
                    case Origin::PrimaryPositronAnnihilationInFlight:
                    case Origin::Fluorescence:
                    case Origin::ImpactIonization:
                        // Lines below 10 keV are outside the model. Unscattered
                        // (a Rayleigh-only x-ray keeps its line energy) -> its
                        // family's atom; Compton-scattered -> the fluorescence
                        // continuum.
                        if (secondary.energy_keV
                            < ceelo::full_spectrum::kFluorescenceMinimumKeV)
                            continue;
                        component = secondary.compton_interactions == 0
                            ? line_families.family_of(secondary.energy_keV)
                            : FluorescenceContinuum;
                        atom = component != FluorescenceContinuum;
                        break;
                    case Origin::Unknown:
                        continue;
                    }
                    tally_photon(
                        secondary.energy_keV, result.bias_weight,
                        component, atom);
                }
                ++events;
            }
            std::vector<Tally> area_vector(areas.begin(), areas.end() - 1);
            if (events >= options.min_events
                && ceelo::full_spectrum::campaign::precision_reached(
                    events, area_vector, coarse))
                break;
        }

        const double analytic_direct = std::exp(-optical_depth);
        std::cout << std::setprecision(12);
        std::cout
            << "cell_kind,energy_keV,Z,tau,areal_density_g_cm2,cos_theta,"
               "detector,component,bin_lo_keV,bin_hi_keV,sum_w,sum_w2,"
               "events,seed,analytic_direct,rayleigh_enabled,"
               "compton_doppler_enabled,positron_transport_enabled,"
               "fluorescence_enabled,hydrogen_mass_fraction,"
               "impact_ionization_enabled\n";
        for (std::size_t c = 0; c < NumComponents; ++c) {
            const auto component = static_cast<Component>(c);
            const bool line = component == Direct || component == Elastic
                || component == Pair511Atom || component == Absorbed
                || component == FluorescenceKAlpha || component == FluorescenceKBeta
                || component == FluorescenceL;
            if (line) {
                double energy = component == Pair511Atom ? 510.99895
                    : (component == Absorbed ? -1.0 : options.energy_keV);
                if (component == FluorescenceKAlpha) energy = line_families.mean_keV[0];
                if (component == FluorescenceKBeta) energy = line_families.mean_keV[1];
                if (component == FluorescenceL) energy = line_families.mean_keV[2];
                const Tally& tally = areas[c];
                std::cout << "source," << options.energy_keV << ','
                    << options.Z << ',' << optical_depth << ','
                    << areal_density << ",,,"
                    << component_name(component, options.interaction_orders)
                    << ',' << energy << ',' << energy << ','
                    << tally.sum_w << ',' << tally.sum_w2 << ','
                    << events << ',' << options.seed << ','
                    << analytic_direct << ','
                    << options.source_rayleigh << ','
                    << options.source_compton_doppler << ','
                    << source.source_positron_transport_enabled() << ','
                    << source.source_fluorescence_enabled() << ','
                    << options.hydrogen_mass_fraction << ','
                    << source.source_impact_ionization_enabled() << '\n';
                continue;
            }
            for (int b = 0; b < options.bins; ++b) {
                const Tally& tally = spectra[c][static_cast<std::size_t>(b)];
                const double lo = options.energy_keV * b / options.bins;
                const double hi = options.energy_keV * (b + 1) / options.bins;
                std::cout << "source," << options.energy_keV << ','
                    << options.Z << ',' << optical_depth << ','
                    << areal_density << ",,,"
                    << component_name(component, options.interaction_orders)
                    << ',' << lo << ',' << hi << ','
                    << tally.sum_w << ',' << tally.sum_w2 << ','
                    << events << ',' << options.seed << ','
                    << analytic_direct << ','
                    << options.source_rayleigh << ','
                    << options.source_compton_doppler << ','
                    << source.source_positron_transport_enabled() << ','
                    << source.source_fluorescence_enabled() << ','
                    << options.hydrogen_mass_fraction << ','
                    << source.source_impact_ionization_enabled() << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_source_escape_campaign: "
                  << error.what() << '\n';
        return 2;
    }
}
