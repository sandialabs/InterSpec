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

#include "full_spectrum/FullSpectrumModel.h"
#include "materials/Material.h"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

using ceelo::full_spectrum::SpectralComponent;

const char* component_name(SpectralComponent component) {
    switch (component) {
    case SpectralComponent::Uncollided:
        return "uncollided";
    case SpectralComponent::ElasticLine:
        return "elastic_line";
    case SpectralComponent::SingleCompton:
        return "single_compton";
    case SpectralComponent::MultipleCompton:
        return "multiple_compton";
    case SpectralComponent::DoubleCompton:
        return "double_compton";
    case SpectralComponent::HigherCompton:
        return "higher_compton";
    case SpectralComponent::ComptonRayleigh:
        return "compton_rayleigh";
    case SpectralComponent::SingleComptonRayleigh:
        return "single_compton_rayleigh";
    case SpectralComponent::DoubleComptonRayleigh:
        return "double_compton_rayleigh";
    case SpectralComponent::HigherComptonRayleigh:
        return "higher_compton_rayleigh";
    case SpectralComponent::MultipleComptonRayleigh:
        return "multiple_compton_rayleigh";
    case SpectralComponent::PairAnnihilation:
        return "pair_511";
    case SpectralComponent::PairAnnihilationContinuum:
        return "pair_annihilation_continuum";
    case SpectralComponent::PairAnnihilationInFlight:
        return "pair_annihilation_in_flight";
    case SpectralComponent::SourceElectronBrems:
        return "source_electron_brems";
    case SpectralComponent::PhotoelectronBrems:
        return "photoelectron_brems";
    case SpectralComponent::ComptonElectronBrems:
        return "compton_electron_brems";
    case SpectralComponent::PairElectronBrems:
        return "pair_electron_brems";
    case SpectralComponent::PhotonuclearSecondary:
        return "photonuclear_secondary";
    case SpectralComponent::SourceFluorescenceKAlpha:
        return "fluorescence_k_alpha";
    case SpectralComponent::SourceFluorescenceKBeta:
        return "fluorescence_k_beta";
    case SpectralComponent::SourceFluorescenceL:
        return "fluorescence_l";
    default:
        return "unsupported";
    }
}

void read_component(
    std::istream& input,
    const std::string& expected_label,
    ceelo::full_spectrum::MultipleScatterRbfComponent& component)
{
    std::string label;
    int degree = 0;
    std::size_t count = 0;
    if (!(input >> label >> degree) || label != expected_label
        || degree + 1
            != static_cast<int>(
                ceelo::full_spectrum::kMultipleScatterBasisSize)) {
        throw std::runtime_error(
            "invalid multiple-scatter component header: "
            + expected_label);
    }
    // The caller sets feature_dimension before calling: 7 for the secondary
    // surrogate's process coordinate, the default 3 (raw) otherwise.
    const std::size_t dimension = std::min<std::size_t>(
        component.feature_dimension,
        ceelo::full_spectrum::kSourceRbfMaxFeatures);
    for (std::size_t index = 0; index < dimension; ++index)
        input >> component.yield_bandwidth[index];
    for (std::size_t index = 0; index < dimension; ++index)
        input >> component.shape_bandwidth[index];
    input >> count;
    if (!input)
        throw std::runtime_error(
            "truncated multiple-scatter component header: "
            + expected_label);
    component.nodes.clear();
    component.nodes.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        std::string node_key;
        std::string node_label;
        ceelo::full_spectrum::MultipleScatterRbfNode node;
        if (!(input >> node_key >> node_label)
            || node_key != "node" || node_label != expected_label) {
            throw std::runtime_error(
                "invalid multiple-scatter node: " + expected_label);
        }
        for (std::size_t axis = 0; axis < dimension; ++axis)
            input >> node.feature[axis];
        input >> node.log_conditional_yield;
        for (double& value : node.shape_coefficients)
            input >> value;
        if (!input)
            throw std::runtime_error(
                "truncated multiple-scatter node: " + expected_label);
        component.nodes.push_back(std::move(node));
    }
}

void read_multiple_surrogate(
    const std::string& path,
    ceelo::full_spectrum::SourceModelPayload& payload)
{
    std::ifstream input(path);
    std::string magic;
    int version = 0;
    if (!input || !(input >> magic >> version)
        || magic != "CEELO_MULTIPLE_SURROGATE" || version != 1) {
        throw std::runtime_error(
            "invalid multiple-scatter surrogate: " + path);
    }
    std::string key;
    if (!(input >> key
          >> payload.multiple_scatter_surrogate
                .inverse_energy_warp_scale)
        || key != "inverse_energy_warp_scale") {
        throw std::runtime_error(
            "missing multiple-scatter inverse-energy warp scale: "
            + path);
    }
    if (!(input >> key
          >> payload.multiple_scatter_surrogate
                .absorption_envelope_strength
          >> payload.multiple_scatter_surrogate
                .absorption_envelope_power
          >> payload.multiple_scatter_surrogate
                .absorption_envelope_max_ratio)
        || key != "absorption_envelope") {
        throw std::runtime_error(
            "missing multiple-scatter absorption envelope: " + path);
    }
    {
        std::string basis;
        if (!(input >> key >> basis) || key != "basis"
            || (basis != "cubic_bspline"
                && basis != "bernstein")) {
            throw std::runtime_error(
                "invalid multiple-scatter basis: " + path);
        }
        payload.multiple_scatter_surrogate.use_cubic_bspline =
            basis == "cubic_bspline";
    }
    std::string mode;
    if (!(input >> key >> mode) || key != "higher_mode"
        || mode != "merged") {
        throw std::runtime_error(
            "unsupported multiple-scatter higher mode: " + path);
    }
    payload.multiple_scatter_surrogate.merge_higher_rayleigh =
        payload.enable_rayleigh;
    std::size_t yield_count = 0;
    auto& yield_model =
        payload.multiple_scatter_surrogate.higher_yield;
    if (!(input >> key >> yield_model.degree)
        || key != "yield_model"
        || yield_model.degree + 1
            > ceelo::full_spectrum::kMultipleScatterYieldBasisSize) {
        throw std::runtime_error(
            "invalid multiple-scatter yield model: " + path);
    }
    for (double& bandwidth : yield_model.bandwidth)
        input >> bandwidth;
    input >> yield_model.minimum_tau >> yield_model.maximum_tau
          >> yield_count;
    if (!input)
        throw std::runtime_error(
            "truncated multiple-scatter yield model: " + path);
    yield_model.nodes.reserve(yield_count);
    for (std::size_t index = 0; index < yield_count; ++index) {
        ceelo::full_spectrum::MultipleScatterYieldNode node;
        if (!(input >> key) || key != "yield_node")
            throw std::runtime_error(
                "invalid multiple-scatter yield node: " + path);
        for (double& value : node.feature)
            input >> value;
        for (std::size_t k = 0; k <= yield_model.degree; ++k)
            input >> node.coefficients[k];
        if (!input)
            throw std::runtime_error(
                "truncated multiple-scatter yield node: " + path);
        yield_model.nodes.push_back(std::move(node));
    }
    read_component(
        input, "higher_aggregate",
        payload.multiple_scatter_surrogate.higher_compton);
    payload.multiple_scatter_surrogate
        .higher_compton_rayleigh.nodes.clear();
    std::string end;
    if (!(input >> end) || end != "END")
        throw std::runtime_error(
            "missing END in multiple-scatter surrogate: " + path);
}

void read_secondary_surrogate(
    const std::string& path,
    ceelo::full_spectrum::SourceModelPayload& payload)
{
    std::ifstream input(path);
    std::string magic;
    int version = 0;
    if (!input || !(input >> magic >> version)
        || magic != "CEELO_SOURCE_SECONDARY_SURROGATE" || version != 1) {
        throw std::runtime_error(
            "invalid source-secondary surrogate: " + path);
    }
    auto& surrogate = payload.secondary_surrogate;
    std::string key;
    if (!(input >> key >> surrogate.inverse_energy_warp_scale)
        || key != "inverse_energy_warp_scale") {
        throw std::runtime_error(
            "missing secondary inverse-energy warp: " + path);
    }
    if (!(input >> key >> surrogate.absorption_envelope_strength
          >> surrogate.absorption_envelope_power
          >> surrogate.absorption_envelope_max_ratio)
        || key != "absorption_envelope") {
        throw std::runtime_error(
            "missing secondary absorption envelope: " + path);
    }
    std::string basis;
    if (!(input >> key >> basis) || key != "basis"
        || (basis != "cubic_bspline" && basis != "bernstein")) {
        throw std::runtime_error(
            "invalid secondary basis: " + path);
    }
    surrogate.use_cubic_bspline = basis == "cubic_bspline";
    // The 7-D "process" coordinate: (log E, Z/92, log tau) plus the
    // parent-process optical depth, the descendant's own absorption depth,
    // thickness in radiation lengths, and the hydrogen mass fraction. See
    // detail::source_process_features.
    const std::uint32_t feature_dimension = static_cast<std::uint32_t>(
        ceelo::full_spectrum::kSourceRbfMaxFeatures);
    double local_linear_ridge = 0.0;
    std::uint32_t local_linear_feature_index = 2;
    if (!(input >> key >> local_linear_ridge >> local_linear_feature_index)
        || key != "local_linear") {
        throw std::runtime_error(
            "missing secondary local_linear line: " + path);
    }
    if (!(input >> key >> surrogate.tail_decay_kappa)
        || key != "tail_decay") {
        throw std::runtime_error(
            "missing secondary tail_decay line: " + path);
    }
    for (auto* component : surrogate.components()) {
        component->feature_dimension = feature_dimension;
        component->local_linear_ridge = local_linear_ridge;
        component->local_linear_feature_index = local_linear_feature_index;
        component->tail_decay_kappa = surrogate.tail_decay_kappa;
    }
    // In SourceSecondarySurrogatePayload::components() order, which is the
    // fitter's ALL_COMPONENTS.
    static constexpr const char* kLabels[] = {
        "pair_511_atom", "pair_annihilation_continuum",
        "pair_annihilation_in_flight", "photoelectron_brems",
        "compton_electron_brems", "pair_electron_brems",
        "fluorescence_k_alpha_atom", "fluorescence_k_beta_atom",
        "fluorescence_l_atom"};
    static_assert(std::size(kLabels)
                  == ceelo::full_spectrum::SourceSecondarySurrogatePayload::kComponentCount);
    {
        const auto components = surrogate.components();
        for (std::size_t c = 0; c < components.size(); ++c)
            read_component(input, kLabels[c], *components[c]);
    }
    std::string end;
    if (!(input >> end) || end != "END")
        throw std::runtime_error(
            "missing END in source-secondary surrogate: " + path);
    // Derive each component's deepest populated optical depth from its own
    // nodes, matching what payload deserialization does.
    for (auto* component : surrogate.components())
        component->refresh_tail_decay_tau();
}

/// Fitted (E, Z, tau) photonuclear yield/shape nodes from fit_photonuclear.py.
void read_photonuclear_surrogate(
    const std::string& path,
    ceelo::full_spectrum::SourceModelPayload& payload)
{
    std::ifstream input(path);
    std::string magic;
    int version = 0;
    if (!input || !(input >> magic >> version)
        || magic != "CEELO_PHOTONUCLEAR_SURROGATE" || version != 1) {
        throw std::runtime_error(
            "invalid photonuclear surrogate: " + path);
    }
    auto& photonuclear = payload.photonuclear;
    std::string key;
    if (!(input >> key >> photonuclear.inverse_energy_warp_scale)
        || key != "inverse_energy_warp_scale") {
        throw std::runtime_error(
            "missing photonuclear inverse-energy warp scale: " + path);
    }
    std::string basis;
    if (!(input >> key >> basis) || key != "basis"
        || (basis != "cubic_bspline" && basis != "bernstein")) {
        throw std::runtime_error(
            "invalid photonuclear basis: " + path);
    }
    photonuclear.use_cubic_bspline = basis == "cubic_bspline";
    if (!(input >> key >> photonuclear.escape_kappa)
        || key != "escape_kappa") {
        throw std::runtime_error(
            "missing photonuclear escape kappa: " + path);
    }
    read_component(input, "photonuclear", photonuclear.secondary);
    std::string end;
    if (!(input >> end) || end != "END") {
        throw std::runtime_error(
            "missing END in photonuclear surrogate: " + path);
    }
}

void read_photonuclear_cross_sections(
    const std::string& path,
    ceelo::full_spectrum::SourceModelPayload& payload)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error(
            "cannot open photonuclear cross-section table: " + path);
    std::string line;
    if (!std::getline(input, line)
        || line != "energy_keV,Z,sigma_barn,mu_over_rho_cm2_g") {
        throw std::runtime_error(
            "invalid photonuclear cross-section header: " + path);
    }
    std::vector<ceelo::full_spectrum::PhotonuclearCrossSectionRow> rows;
    std::array<bool, 92> populated{};
    double current_energy = -1.0;
    auto finish_row = [&]() {
        if (current_energy < 0.0)
            return;
        for (std::size_t index = 0; index < populated.size(); ++index) {
            if (!populated[index]) {
                throw std::runtime_error(
                    "incomplete photonuclear Z grid at "
                    + std::to_string(current_energy) + " keV: " + path);
            }
        }
        rows.back().energy_keV = current_energy;
    };
    std::size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty())
            continue;
        for (char& character : line) {
            if (character == ',')
                character = ' ';
        }
        std::istringstream fields(line);
        double energy = 0.0;
        int Z = 0;
        double sigma_barn = 0.0;
        double mass_attenuation = 0.0;
        if (!(fields >> energy >> Z >> sigma_barn >> mass_attenuation)
            || Z < 1 || Z > 92 || !(energy > 0.0)
            || mass_attenuation < 0.0) {
            throw std::runtime_error(
                "invalid photonuclear row " + std::to_string(line_number)
                + ": " + path);
        }
        if (current_energy < 0.0 || energy != current_energy) {
            finish_row();
            if (current_energy >= 0.0 && energy <= current_energy) {
                throw std::runtime_error(
                    "photonuclear energies are not strictly ordered: "
                    + path);
            }
            rows.emplace_back();
            populated.fill(false);
            current_energy = energy;
        }
        const std::size_t index = static_cast<std::size_t>(Z - 1);
        if (populated[index]) {
            throw std::runtime_error(
                "duplicate photonuclear Z row: " + path);
        }
        populated[index] = true;
        rows.back().mass_attenuation_cm2_g[index] = mass_attenuation;
    }
    finish_row();
    if (rows.size() < 2)
        throw std::runtime_error(
            "photonuclear table needs at least two energy rows: " + path);
    payload.photonuclear.cross_sections = std::move(rows);
    payload.enable_photonuclear = true;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto options =
            ceelo::full_spectrum::campaign::parse_options(argc, argv, false);
        const double f_H = options.hydrogen_mass_fraction;
        std::vector<ceelo::MaterialComponent> composition;
        if (f_H > 0.0)
            composition.push_back({std::uint8_t{1}, f_H});
        composition.push_back(
            {static_cast<std::uint8_t>(options.Z), 1.0 - f_H});
        ceelo::Material element("equivalent-element", 1.0, composition);
        std::vector<double> edges(
            static_cast<std::size_t>(options.bins) + 1);
        for (int bin = 0; bin <= options.bins; ++bin) {
            edges[static_cast<std::size_t>(bin)] =
                options.energy_keV * bin / options.bins;
        }
        ceelo::full_spectrum::SourceModelPayload payload;
        if (!options.payload_file.empty()) {
            std::ifstream in(options.payload_file);
            if (!in)
                throw std::runtime_error(
                    "cannot open payload file: " + options.payload_file);
            std::ostringstream text;
            text << in.rdbuf();
            payload = ceelo::full_spectrum::SourceGeomFluxModel::deserialize(
                          text.str())
                          .payload();
        } else {
            payload.enable_rayleigh = options.source_rayleigh;
            payload.enable_compton_doppler =
                options.source_compton_doppler;
            payload.enable_source_electron_brems = true;
            if (!options.multiple_surrogate_file.empty())
                read_multiple_surrogate(
                    options.multiple_surrogate_file, payload);
            if (!options.secondary_surrogate_file.empty())
                read_secondary_surrogate(
                    options.secondary_surrogate_file, payload);
            if (!options.photonuclear_xs_file.empty())
                read_photonuclear_cross_sections(
                    options.photonuclear_xs_file, payload);
            if (!options.photonuclear_surrogate_file.empty())
                read_photonuclear_surrogate(
                    options.photonuclear_surrogate_file, payload);
        }
        // Applied after either load path, so they override a payload's
        // stored values as well as the struct defaults.
        if (options.elastic_escape_scale >= 0.0)
            payload.elastic_escape_scale = options.elastic_escape_scale;
        if (options.elastic_compton_damping >= 0.0) {
            payload.elastic_compton_damping =
                options.elastic_compton_damping;
        }
        const bool rayleigh_enabled = payload.enable_rayleigh;
        const bool doppler_enabled = payload.enable_compton_doppler;
        const auto incident_xs =
            element.macroscopic_xs(options.energy_keV * 1e-3);
        const double mu_mass = (
            incident_xs.mu_pe + incident_xs.mu_cs + incident_xs.mu_pp
            + (rayleigh_enabled ? incident_xs.mu_rs : 0.0))
            / element.density();
        const double areal_density =
            options.areal_density_g_cm2 >= 0.0
            ? options.areal_density_g_cm2
            : (options.tau == 0.0 ? 0.0 : options.tau / mu_mass);
        const double optical_depth = areal_density * mu_mass;
        const ceelo::full_spectrum::SourceGeomFluxModel model(payload);
        if (!options.write_payload_file.empty()) {
            std::ofstream out(options.write_payload_file);
            if (!out)
                throw std::runtime_error(
                    "cannot open payload output file: "
                    + options.write_payload_file);
            out << model.serialize();
            if (!out)
                throw std::runtime_error(
                    "failed writing payload file: "
                    + options.write_payload_file);
        }
        const auto measure = model.evaluate(
            options.energy_keV, areal_density,
            static_cast<double>(options.Z),
            options.hydrogen_mass_fraction, edges);

        std::cout << std::setprecision(12)
                  << "cell_kind,energy_keV,Z,tau,areal_density_g_cm2,"
                     "cos_theta,detector,component,bin_lo_keV,bin_hi_keV,"
                     "sum_w,sum_w2,events,seed,analytic_direct,"
                     "rayleigh_enabled,compton_doppler_enabled,"
                     "hydrogen_mass_fraction\n";
        for (const auto& atom : measure.atoms) {
            std::cout << "source," << options.energy_keV << ','
                      << options.Z << ',' << optical_depth << ','
                      << areal_density << ",,,"
                      << component_name(atom.component) << ','
                      << atom.energy_keV << ',' << atom.energy_keV << ','
                      << atom.yield << ",0,1,0,"
                      << std::exp(-optical_depth) << ','
                      << rayleigh_enabled << ','
                      << doppler_enabled << ','
                      << options.hydrogen_mass_fraction << '\n';
        }
        for (const auto& bin : measure.continuous) {
            std::cout << "source," << options.energy_keV << ','
                      << options.Z << ',' << optical_depth << ','
                      << areal_density << ",,,"
                      << component_name(bin.component) << ','
                      << bin.lower_keV << ',' << bin.upper_keV << ','
                      << bin.integral << ",0,1,0,"
                      << std::exp(-optical_depth) << ','
                      << rayleigh_enabled << ','
                      << doppler_enabled << ','
                      << options.hydrogen_mass_fraction << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_source_model_evaluate: "
                  << error.what() << '\n';
        return 2;
    }
}
