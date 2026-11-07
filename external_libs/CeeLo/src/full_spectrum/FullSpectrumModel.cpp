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

#include "cross_sections/CrossSectionData.h"
#include "materials/Material.h"
#include "physics/ElectronCsda.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ceelo::full_spectrum {
namespace {

PartialMassAttenuation compound_mass_cross_sections(
    double energy_keV,
    const std::vector<std::pair<int, double>>& mass_fractions)
{
    PartialMassAttenuation answer;
    for (const auto& component : mass_fractions) {
        const PartialMassAttenuation element =
            element_mass_cross_sections(component.first, energy_keV);
        answer.pe += component.second * element.pe;
        answer.cs += component.second * element.cs;
        answer.rs += component.second * element.rs;
        answer.pp += component.second * element.pp;
    }
    return answer;
}

template <std::size_t N>
void write_array(std::ostream& out, const char* name,
                 const std::array<double, N>& values) {
    out << name;
    for (double value : values)
        out << ' ' << std::setprecision(17) << value;
    out << '\n';
}

template <std::size_t N>
void read_array(std::istream& in, std::array<double, N>& values) {
    for (double& value : values) {
        if (!(in >> value))
            throw std::runtime_error("Truncated FullSpecApprox coefficient array");
    }
}

// Scattered-RBF components carry their own feature width. Only the first
// `feature_dimension` bandwidths and node coordinates are written.
// Components without that member (the detector response) keep their full
// fixed width; the exact-match overload wins for the ones that have it.
template <typename Component>
std::size_t rbf_feature_dimension(const Component& component) {
    return component.yield_bandwidth.size();
}

inline std::size_t rbf_feature_dimension(
    const ceelo::full_spectrum::MultipleScatterRbfComponent& component) {
    return std::min<std::size_t>(
        component.feature_dimension,
        ceelo::full_spectrum::kSourceRbfMaxFeatures);
}

/// Payload keys of the source secondary components, in
/// SourceSecondarySurrogatePayload::components() order; a node line is the key
/// plus "_node".
constexpr std::array<const char*, SourceSecondarySurrogatePayload::kComponentCount>
    kSecondaryComponentKeys{{
        "source_pair_atom_rbf", "source_pair_continuum_rbf",
        "source_pair_in_flight_rbf", "source_photoelectron_brems_rbf",
        "source_compton_electron_brems_rbf", "source_pair_electron_brems_rbf",
        "source_fluorescence_k_alpha_rbf", "source_fluorescence_k_beta_rbf",
        "source_fluorescence_l_rbf"}};

MultipleScatterRbfComponent* secondary_component_for_key(
    SourceSecondarySurrogatePayload& secondary, const std::string& key)
{
    const auto components = secondary.components();
    for (std::size_t c = 0; c < components.size(); ++c)
        if (key == kSecondaryComponentKeys[c]) return components[c];
    return nullptr;
}

MultipleScatterRbfComponent* secondary_component_for_node_key(
    SourceSecondarySurrogatePayload& secondary, const std::string& key)
{
    const auto components = secondary.components();
    for (std::size_t c = 0; c < components.size(); ++c)
        if (key == std::string(kSecondaryComponentKeys[c]) + "_node") return components[c];
    return nullptr;
}

template <typename Component>
void write_rbf_component(
    std::ostream& out,
    const char* name,
    const Component& component)
{
    const std::size_t dimension = rbf_feature_dimension(component);
    out << name;
    for (std::size_t index = 0; index < dimension; ++index)
        out << ' ' << std::setprecision(17)
            << component.yield_bandwidth[index];
    for (std::size_t index = 0; index < dimension; ++index)
        out << ' ' << std::setprecision(17)
            << component.shape_bandwidth[index];
    out << ' ' << component.nodes.size() << '\n';
    for (const auto& node : component.nodes) {
        out << name << "_node";
        for (std::size_t index = 0; index < dimension; ++index)
            out << ' ' << std::setprecision(17) << node.feature[index];
        out << ' ' << std::setprecision(17)
            << node.log_conditional_yield;
        for (double value : node.shape_coefficients)
            out << ' ' << std::setprecision(17) << value;
        out << '\n';
    }
}

template <typename Component>
void read_rbf_component_body(
    std::istream& in,
    const std::string& name,
    Component& component)
{
    std::string key;
    std::size_t count = 0;
    // The caller sets `feature_dimension` (from the section options line)
    // before getting here; it defaults to the 3-D raw coordinate.
    const std::size_t dimension = rbf_feature_dimension(component);
    const auto read_n = [&](auto& values) {
        for (std::size_t index = 0; index < dimension; ++index) {
            if (!(in >> values[index]))
                throw std::runtime_error(
                    "Truncated FullSpecApprox coefficient array");
        }
    };
    read_n(component.yield_bandwidth);
    read_n(component.shape_bandwidth);
    if (!(in >> count))
        throw std::runtime_error(
            "Truncated charged-particle component " + name);
    component.nodes.clear();
    component.nodes.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        typename Component::value_type node;
        if (!(in >> key) || key != name + "_node")
            throw std::runtime_error(
                "Missing charged-particle node for " + name);
        read_n(node.feature);
        if (!(in >> node.log_conditional_yield))
            throw std::runtime_error(
                "Truncated charged-particle node yield for " + name);
        read_array(in, node.shape_coefficients);
        component.nodes.push_back(std::move(node));
    }
}

template <typename Component>
void read_rbf_component(
    std::istream& in,
    const std::string& name,
    Component& component)
{
    std::string key;
    if (!(in >> key) || key != name)
        throw std::runtime_error(
            "Missing charged-particle component " + name);
    read_rbf_component_body(in, name, component);
}

// ---- key-dispatched (source-section) scattered-RBF readers ---------------
// The feature arrays are kSourceRbfMaxFeatures wide but only the component's
// own `feature_dimension` entries appear in the file (3 for the higher-Compton
// and photonuclear components). The secondary surrogate's width is set from
// its options line, which the writer emits before any of its components.

void read_rbf_header(
    std::istream& in,
    ceelo::full_spectrum::MultipleScatterRbfComponent& component)
{
    const std::size_t dimension = std::min<std::size_t>(
        component.feature_dimension,
        ceelo::full_spectrum::kSourceRbfMaxFeatures);
    for (std::size_t index = 0; index < dimension; ++index) {
        if (!(in >> component.yield_bandwidth[index]))
            throw std::runtime_error("Truncated RBF yield bandwidth");
    }
    for (std::size_t index = 0; index < dimension; ++index) {
        if (!(in >> component.shape_bandwidth[index]))
            throw std::runtime_error("Truncated RBF shape bandwidth");
    }
    std::size_t count = 0;
    if (!(in >> count))
        throw std::runtime_error("Truncated RBF node count");
    component.nodes.reserve(count);
}

void read_rbf_node(
    std::istream& in,
    ceelo::full_spectrum::MultipleScatterRbfComponent& component)
{
    const std::size_t dimension = std::min<std::size_t>(
        component.feature_dimension,
        ceelo::full_spectrum::kSourceRbfMaxFeatures);
    ceelo::full_spectrum::MultipleScatterRbfNode node;
    for (std::size_t index = 0; index < dimension; ++index) {
        if (!(in >> node.feature[index]))
            throw std::runtime_error("Truncated RBF node feature");
    }
    if (!(in >> node.log_conditional_yield))
        throw std::runtime_error("Truncated RBF node yield");
    read_array(in, node.shape_coefficients);
    component.nodes.push_back(std::move(node));
}

// Charged-particle payload structured-field IO. The format is documented on
// ChargedParticleSplineField in FullSpectrumModel.h.

// One validator for both codecs, so the embedded and standalone readers accept
// exactly the same anchors (azimuth_quadrature_points = 0 would silently zero
// the side term; a huge count would size a vector by it).
void validate_rate_surface_anchor(const DetectorRateSurfaceAnchor& anchor)
{
    if (!(anchor.radius_cm > 0.0) || !(anchor.length_cm > 0.0)
        || !(anchor.crystal_density_g_cm3 > 0.0)
        || anchor.crystal_composition.empty()
        || !(anchor.envelope_radius_cm >= anchor.radius_cm)
        || !(anchor.envelope_length_cm >= anchor.length_cm)
        || !std::isfinite(anchor.blend_keV) || !std::isfinite(anchor.blend_power)
        || anchor.azimuth_quadrature_points < 2
        || anchor.azimuth_quadrature_points > 256) {
        throw std::runtime_error(
            "Unsupported detector rate-surface anchor");
    }
}

// The embedded geometry block of both detector codecs: "<label> <n>" then the
// n lines of GeometryDescriptor::to_text().
void write_geometry_block(std::ostream& out, const char* label,
                          const GeometryDescriptor& geometry)
{
    const std::string text = geometry.to_text();
    out << label << ' ' << std::count(text.begin(), text.end(), '\n') << '\n'
        << text;
}

GeometryDescriptor read_geometry_block(std::istream& in)
{
    std::size_t count = 0;
    if (!(in >> count) || count < 3 || count > 4096)
        throw std::runtime_error("Invalid detector geometry block");
    std::string line;
    std::getline(in, line);   // the rest of the count line
    std::string text;
    for (std::size_t i = 0; i < count; ++i) {
        if (!std::getline(in, line))
            throw std::runtime_error("Truncated detector geometry block");
        text += line;
        text += '\n';
    }
    GeometryDescriptor geometry = GeometryDescriptor::from_text(text);
    if (!geometry.problems().empty())
        throw std::runtime_error(
            std::string("Invalid detector geometry: ")
            + to_string(geometry.problems().front()));
    return geometry;
}

// Family, radius, length and density follow the geometry; the pilot
// constants and the resolution window stay those of the family's preset.
void apply_payload_geometry(DetectorModelPayload& payload)
{
    if (!payload.geometry)
        return;
    const GeometryDescriptor& geometry = *payload.geometry;
    const MaterialSpec& crystal = geometry.materials.at(
        static_cast<std::size_t>(geometry.crystal_material_index));
    const CylinderDims dims = geometry.cylinder_dims();
    payload.radius_cm = dims.radius_cm;
    payload.length_cm = dims.full_length_cm;
    payload.density_g_cm3 = crystal.density_g_per_cm3;
    set_rate_surface_anchor_geometry(
        payload.detected_rate_surface.anchor, geometry);
}

void write_spline_field(
    std::ostream& out,
    const char* name,
    const ChargedParticleSplineField& field)
{
    const std::size_t axes = field.yield_grid.size();
    out << name << "_domain " << axes;
    for (double value : field.domain)
        out << ' ' << std::setprecision(17) << value;
    out << '\n';
    out << name << "_yield_grid " << axes;
    for (std::uint32_t size : field.yield_grid)
        out << ' ' << size;
    out << '\n';
    out << name << "_yield_coefficients";
    for (double value : field.yield_coefficients)
        out << ' ' << std::setprecision(17) << value;
    out << '\n';
    out << name << "_shape_rank " << field.shape_rank << '\n';
    for (const auto& archetype : field.shape_basis) {
        out << name << "_shape_basis";
        for (double value : archetype)
            out << ' ' << std::setprecision(17) << value;
        out << '\n';
    }
    if (field.shape_rank >= 2) {
        out << name << "_weight_grid " << axes;
        for (std::uint32_t size : field.weight_grid)
            out << ' ' << size;
        out << '\n';
        for (const auto& log_odds : field.weight_log_odds) {
            out << name << "_weight_log_odds";
            for (double value : log_odds)
                out << ' ' << std::setprecision(17) << value;
            out << '\n';
        }
    }
}

void read_spline_field(
    std::istream& in,
    const std::string& name,
    ChargedParticleSplineField& field)
{
    std::string key;
    const auto expect = [&](const std::string& suffix) {
        if (!(in >> key) || key != name + suffix)
            throw std::runtime_error(
                "Missing charged-particle spline record "
                + name + suffix);
    };
    const auto read_grid =
        [&](const std::string& suffix,
            std::vector<std::uint32_t>& grid,
            std::size_t axes) {
            expect(suffix);
            std::size_t declared = 0;
            if (!(in >> declared) || declared != axes)
                throw std::runtime_error(
                    "Unsupported charged-particle spline axis count for "
                    + name);
            grid.resize(axes);
            for (auto& size : grid) {
                if (!(in >> size) || size < 4 || size > 32)
                    throw std::runtime_error(
                        "Invalid charged-particle spline grid for " + name);
            }
        };
    // Three axes today; a future composition axis extends this count.
    constexpr std::size_t axes = 3;
    expect("_domain");
    std::size_t declared = 0;
    if (!(in >> declared) || declared != axes)
        throw std::runtime_error(
            "Unsupported charged-particle spline axis count for " + name);
    field.domain.resize(2 * axes);
    for (auto& value : field.domain) {
        if (!(in >> value) || !std::isfinite(value))
            throw std::runtime_error(
                "Invalid charged-particle spline domain for " + name);
    }
    for (std::size_t axis = 0; axis < axes; ++axis) {
        if (!(field.domain[2 * axis + 1] > field.domain[2 * axis]))
            throw std::runtime_error(
                "Invalid charged-particle spline domain for " + name);
    }
    read_grid("_yield_grid", field.yield_grid, axes);
    std::size_t yield_count = 1;
    for (std::uint32_t size : field.yield_grid)
        yield_count *= size;
    expect("_yield_coefficients");
    field.yield_coefficients.resize(yield_count);
    for (auto& value : field.yield_coefficients) {
        if (!(in >> value) || !std::isfinite(value))
            throw std::runtime_error(
                "Invalid charged-particle spline coefficients for " + name);
    }
    expect("_shape_rank");
    if (!(in >> field.shape_rank) || field.shape_rank > 16)
        throw std::runtime_error(
            "Invalid charged-particle spline shape rank for " + name);
    field.shape_basis.resize(field.shape_rank);
    for (auto& archetype : field.shape_basis) {
        expect("_shape_basis");
        for (auto& value : archetype) {
            if (!(in >> value) || !std::isfinite(value) || value < 0.0)
                throw std::runtime_error(
                    "Invalid charged-particle spline archetype for "
                    + name);
        }
    }
    field.weight_grid.clear();
    field.weight_log_odds.clear();
    if (field.shape_rank >= 2) {
        read_grid("_weight_grid", field.weight_grid, axes);
        std::size_t weight_count = 1;
        for (std::uint32_t size : field.weight_grid)
            weight_count *= size;
        field.weight_log_odds.resize(field.shape_rank - 1);
        for (auto& log_odds : field.weight_log_odds) {
            expect("_weight_log_odds");
            log_odds.resize(weight_count);
            for (auto& value : log_odds) {
                if (!(in >> value) || !std::isfinite(value))
                    throw std::runtime_error(
                        "Invalid charged-particle spline log-odds for "
                        + name);
            }
        }
    }
}

void write_detector_rbf_component(
    std::ostream& out,
    const char* name,
    const DetectorResponseRbfComponent& component)
{
    out << name;
    for (double bandwidth : component.yield_bandwidth)
        out << ' ' << std::setprecision(17) << bandwidth;
    for (double bandwidth : component.shape_bandwidth)
        out << ' ' << std::setprecision(17) << bandwidth;
    out << ' ' << component.basis_size
        << ' ' << component.nodes.size();
    if (component.uses_blocks()) {
        // Block partition of the deposited-energy support. The count comes
        // first so a reader can skip it without knowing the component's
        // landmark structure.
        std::size_t blocks = 0;
        for (std::uint32_t size : component.block_sizes)
            if (size != 0)
                ++blocks;
        out << ' ' << blocks;
        for (std::uint32_t size : component.block_sizes)
            if (size != 0)
                out << ' ' << size;
        out << '\n';
    } else {
        out << " 0\n";
    }
    for (const auto& node : component.nodes) {
        out << name << "_node";
        for (double value : node.feature)
            out << ' ' << std::setprecision(17) << value;
        out << ' ' << std::setprecision(17)
            << node.log_conditional_yield;
        for (std::size_t index = 0;
             index < component.basis_size; ++index) {
            out << ' ' << std::setprecision(17)
                << node.shape_coefficients[index];
        }
        out << '\n';
    }
}

void read_detector_rbf_component_body(
    std::istream& in,
    const std::string& name,
    DetectorResponseRbfComponent& component)
{
    std::string key;
    std::size_t count = 0;
    read_array(in, component.yield_bandwidth);
    read_array(in, component.shape_bandwidth);
    if (!(in >> component.basis_size >> count)
        || component.basis_size < 1
        || component.basis_size > kDetectorResponseBasisSize) {
        throw std::runtime_error(
            "Invalid compact detector component " + name);
    }
    component.block_sizes.fill(0);
    {
        std::size_t blocks = 0;
        if (!(in >> blocks)
            || blocks > detail::kDetectorMaxBlocks) {
            throw std::runtime_error(
                "Invalid detector block count for " + name);
        }
        std::uint32_t total = 0;
        for (std::size_t index = 0; index < blocks; ++index) {
            std::uint32_t size = 0;
            if (!(in >> size) || size < 4
                || size > detail::kDetectorMaxBlockSize) {
                throw std::runtime_error(
                    "Invalid detector block size for " + name);
            }
            component.block_sizes[index] = size;
            total += size;
        }
        if (blocks != 0 && total != component.basis_size) {
            throw std::runtime_error(
                "Detector block sizes must sum to basis size for " + name);
        }
    }
    component.nodes.clear();
    component.nodes.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        DetectorResponseRbfNode node;
        if (!(in >> key) || key != name + "_node")
            throw std::runtime_error(
                "Missing compact detector node for " + name);
        read_array(in, node.feature);
        if (!(in >> node.log_conditional_yield))
            throw std::runtime_error(
                "Truncated compact detector node yield for " + name);
        for (std::size_t coefficient = 0;
             coefficient < component.basis_size; ++coefficient) {
            if (!(in >> node.shape_coefficients[coefficient])) {
                throw std::runtime_error(
                    "Truncated compact detector shape for " + name);
            }
        }
        component.nodes.push_back(std::move(node));
    }
}

void read_detector_rbf_component(
    std::istream& in,
    const std::string& name,
    DetectorResponseRbfComponent& component)
{
    std::string key;
    if (!(in >> key) || key != name)
        throw std::runtime_error(
            "Missing compact detector component " + name);
    read_detector_rbf_component_body(in, name, component);
}

} // namespace

SmoothElementCrossSection::PhotonuclearInterpolation
SmoothElementCrossSection::photonuclear_interpolation(
    const PhotonuclearSourcePayload& payload, double energy_keV)
{
    PhotonuclearInterpolation result;
    if (payload.cross_sections.empty()
        || energy_keV < payload.cross_sections.front().energy_keV) {
        result.below_table = true;
        return result;
    }
    std::size_t upper = 0;
    while (upper < payload.cross_sections.size()
           && payload.cross_sections[upper].energy_keV < energy_keV) {
        ++upper;
    }
    if (upper == 0) {
        result.lower = result.upper = 0;
    } else if (upper >= payload.cross_sections.size()) {
        result.lower = result.upper = payload.cross_sections.size() - 1;
    } else {
        result.lower = upper - 1;
        result.upper = upper;
        const double denominator =
            payload.cross_sections[upper].energy_keV
            - payload.cross_sections[result.lower].energy_keV;
        result.fraction =
            (energy_keV - payload.cross_sections[result.lower].energy_keV)
            / denominator;
    }
    return result;
}

PartialMassAttenuation element_mass_cross_sections(
    int Z, double energy_keV)
{
    const auto& data = CrossSectionData::instance();
    const auto microscopic =
        data.all_cross_sections(Z, std::max(energy_keV, 0.001) * 1e-3);
    const double atoms_per_gram =
        kAvogadro / data.atomic_weight(Z) * kBarnToCm2;
    return {
        atoms_per_gram * microscopic.sigma_pe,
        atoms_per_gram * microscopic.sigma_cs,
        atoms_per_gram * microscopic.sigma_rs,
        atoms_per_gram * microscopic.sigma_pp};
}

void set_rate_surface_anchor_geometry(DetectorRateSurfaceAnchor& anchor,
                                      const GeometryDescriptor& geometry)
{
    if (geometry.shape != DetectorShape::Cylinder)
        throw std::invalid_argument(
            "the rate-surface anchor needs a cylindrical crystal");
    // build_geometry() runs problems() and gives the envelope exactly as the
    // detector campaign computes it.
    std::vector<std::unique_ptr<Material>> materials;
    const Geometry traced = geometry.build_geometry(materials);
    const CylinderDims dims = geometry.cylinder_dims();
    const MaterialSpec& crystal = geometry.materials.at(
        static_cast<std::size_t>(geometry.crystal_material_index));
    anchor.radius_cm = dims.radius_cm;
    anchor.length_cm = dims.full_length_cm;
    anchor.front_dead_cm = geometry.dead_layer ? geometry.dead_layer->front : 0.0;
    anchor.side_dead_cm = geometry.dead_layer ? geometry.dead_layer->side : 0.0;
    anchor.back_dead_cm = geometry.dead_layer ? geometry.dead_layer->back : 0.0;
    anchor.bore_radius_cm = geometry.bore ? geometry.bore->radius : 0.0;
    anchor.crystal_composition = crystal.composition;
    anchor.crystal_density_g_cm3 = crystal.density_g_per_cm3;
    anchor.front_layers.clear();
    anchor.side_layers.clear();
    anchor.back_layers.clear();
    const auto layer = [&](int index, double thickness) {
        const MaterialSpec& m =
            geometry.materials.at(static_cast<std::size_t>(index));
        return AnchorLayer{m.composition, m.density_g_per_cm3, thickness};
    };
    // As trace_ray() stacks the shells: each covers the front and side of the
    // shell inside it, and the back only where its z_end reaches past it.
    double inner_z_max = dims.full_length_cm;
    for (const LayerSpec& l : geometry.layers) {
        if (l.front_thickness_cm > 0.0)
            anchor.front_layers.push_back(
                layer(l.material_index, l.front_thickness_cm));
        if (l.side_thickness_cm > 0.0)
            anchor.side_layers.push_back(
                layer(l.material_index, l.side_thickness_cm));
        const double shell_z_max = std::max(inner_z_max, l.z_end_cm);
        if (shell_z_max > inner_z_max)
            anchor.back_layers.push_back(
                layer(l.material_index, shell_z_max - inner_z_max));
        inner_z_max = shell_z_max;
    }
    if (geometry.collimator && geometry.collimator->side_thickness_cm > 0.0)
        anchor.side_layers.push_back(layer(
            geometry.collimator->material_index,
            geometry.collimator->side_thickness_cm));
    const auto [z_min, z_max] = traced.outer_z_extent();
    anchor.envelope_radius_cm = traced.outer_bounding_radius();
    anchor.envelope_length_cm = z_max - z_min;
}

double passive_landmark_keV(const DetectorRateSurfaceAnchor& anchor)
{
    const auto depth = [&](double energy_keV) {
        return detail::layer_stack_optical_depth(anchor.front_layers, energy_keV)
            + detail::total_mass_attenuation(anchor.crystal_composition,
                                             energy_keV)
                * anchor.crystal_density_g_cm3 * anchor.front_dead_cm;
    };
    // Scan down from 1000 keV for the first (highest-energy) crossing, then
    // bisect it: an absorption edge in the stack can make the depth
    // non-monotone below that.
    constexpr int kSteps = 400;
    const double lo_keV = 10.0, hi_keV = 1000.0;
    const double ratio = std::pow(hi_keV / lo_keV, 1.0 / kSteps);
    double upper = hi_keV;
    if (depth(upper) >= kPassiveLandmarkOpticalDepth)
        return hi_keV;
    for (int step = 0; step < kSteps; ++step) {
        const double lower = upper / ratio;
        if (depth(lower) >= kPassiveLandmarkOpticalDepth) {
            double a = lower, b = upper;
            for (int i = 0; i < 60; ++i) {
                const double mid = std::sqrt(a * b);
                (depth(mid) >= kPassiveLandmarkOpticalDepth ? a : b) = mid;
            }
            return std::sqrt(a * b);
        }
        upper = lower;
    }
    return lo_keV;
}

void set_detector_geometry(DetectorModelPayload& payload,
                           GeometryDescriptor geometry)
{
    payload.geometry = std::move(geometry);
    apply_payload_geometry(payload);
    payload.passive_landmark_keV =
        passive_landmark_keV(payload.detected_rate_surface.anchor);
}

double element_csda_range_g_cm2(
    int Z, double kinetic_energy_keV, bool is_positron)
{
    return ElectronCsda::instance().range_gcm2(
        Z, kinetic_energy_keV, is_positron);
}

double element_inverse_radiation_length_cm2_g(int Z)
{
    const double radiation_length =
        ElectronCsda::radiation_length_gcm2_element(
            Z, CrossSectionData::instance().atomic_weight(Z));
    return radiation_length > 0.0 ? 1.0 / radiation_length : 0.0;
}

namespace {

ElementFluorescence build_element_fluorescence(int Z)
{
    ElementFluorescence out;
    const auto& xs = CrossSectionData::instance();
    const FluorescenceData* k = xs.fluorescence(Z);
    const LFluorescenceData* l = xs.l_fluorescence(Z);
    std::vector<FluorescenceLineYield> lines;

    // An L vacancy in subshell s ends radiating from subshell t with
    // probability P[s][t] (radiate with the retained yield first, then
    // Coster-Kronig f13 / f12 from L1 and f23 from L2), as sample_l_relaxation.
    std::array<std::array<double, 3>, 3> P{};
    if (l) {
        const double w1 = l->sub[0].fluorescence_yield, w2 = l->sub[1].fluorescence_yield,
                     w3 = l->sub[2].fluorescence_yield;
        const double f12 = l->coster_kronig[0], f13 = l->coster_kronig[1],
                     f23 = l->coster_kronig[2];
        P[0] = {w1, f12 * w2, f13 * w3 + f12 * f23 * w3};
        P[1] = {0.0, w2, f23 * w3};
        P[2] = {0.0, 0.0, w3};
    }
    // L lines per L-subshell vacancy, and the L vacancies a K vacancy leaves.
    std::array<double, 3> l_per_k{};
    if (k && k->num_lines > 0) {
        out.k_edge_keV = k->k_edge_keV;
        out.omega_K = k->fluorescence_yield;
        for (int i = 0; i < k->num_lines; ++i) {
            FluorescenceLineYield line;
            line.energy_keV = k->line_energy_keV[i];
            line.family = static_cast<int>(k_line_family(k->line_filling_subshell[i]));
            line.per_k = k->fluorescence_yield * k->line_probability[i];
            lines.push_back(line);
            const int filling = k->line_filling_subshell[i];
            if (filling >= kEadlL1 && filling <= kEadlL3)
                l_per_k[filling - kEadlL1] += line.per_k;
        }
    }
    if (l) {
        for (int t = 0; t < 3; ++t) {
            const LSubshellFluor& sub = l->sub[t];
            if (sub.num_lines > 0)
                for (int s = 0; s < 3; ++s) out.l_edge_keV[s] = l->sub[s].edge_keV;
            for (int i = 0; i < sub.num_lines; ++i) {
                FluorescenceLineYield line;
                line.energy_keV = sub.line_energy_keV[i];
                line.family = 2;
                for (int s = 0; s < 3; ++s)
                    line.per_l[s] = P[s][t] * sub.line_probability[i];
                for (int s = 0; s < 3; ++s)
                    line.per_k += l_per_k[s] * line.per_l[s];
                lines.push_back(line);
            }
        }
    }
    // Lines below 10 keV (light-element K lines) are outside the model's
    // validity and the photon tables; the campaign does not tally them either.
    lines.erase(std::remove_if(lines.begin(), lines.end(),
                               [](const FluorescenceLineYield& line) {
                                   return line.energy_keV < kFluorescenceMinimumKeV;
                               }),
                lines.end());
    // Compress each family: lines under 0.5% of the family's yield fold into
    // the nearest kept line, then lines within 0.05 keV merge.
    const auto yield = [](const FluorescenceLineYield& line) {
        return line.per_k + line.per_l[0] + line.per_l[1] + line.per_l[2];
    };
    const auto add_into = [](FluorescenceLineYield& to, const FluorescenceLineYield& from) {
        to.per_k += from.per_k;
        for (int s = 0; s < 3; ++s) to.per_l[s] += from.per_l[s];
    };
    for (int family = 0; family < kFluorescenceFamilies; ++family) {
        std::vector<FluorescenceLineYield> members;
        double total = 0.0;
        for (const auto& line : lines)
            if (line.family == family) {
                members.push_back(line);
                total += yield(line);
            }
        if (members.empty() || !(total > 0.0)) continue;
        std::vector<FluorescenceLineYield> kept, dropped;
        for (const auto& line : members)
            (yield(line) >= 0.005 * total ? kept : dropped).push_back(line);
        for (const auto& line : dropped) {
            auto nearest = std::min_element(kept.begin(), kept.end(),
                [&](const auto& a, const auto& b) {
                    return std::abs(a.energy_keV - line.energy_keV)
                        < std::abs(b.energy_keV - line.energy_keV);
                });
            add_into(*nearest, line);
        }
        std::sort(kept.begin(), kept.end(),
                  [](const auto& a, const auto& b) { return a.energy_keV < b.energy_keV; });
        std::vector<FluorescenceLineYield> merged;
        for (const auto& line : kept) {
            if (!merged.empty() && line.energy_keV - merged.back().energy_keV < 0.05) {
                FluorescenceLineYield& m = merged.back();
                const double wm = yield(m), wl = yield(line);
                m.energy_keV = (m.energy_keV * wm + line.energy_keV * wl) / (wm + wl);
                add_into(m, line);
            } else {
                merged.push_back(line);
            }
        }
        out.lines.insert(out.lines.end(), merged.begin(), merged.end());
    }
    return out;
}

}  // namespace

const ElementFluorescence& element_fluorescence(int Z)
{
    static const std::array<ElementFluorescence, 93> table = [] {
        std::array<ElementFluorescence, 93> t{};
        for (int z = 1; z <= 92; ++z) t[z] = build_element_fluorescence(z);
        return t;
    }();
    return table[(Z >= 1 && Z <= 92) ? Z : 0];
}

double element_impact_xray_yield(int Z, double kinetic_energy_keV)
{
    if (Z < 1 || Z > 92 || !(kinetic_energy_keV > 0.0)) return 0.0;
    static const std::array<std::optional<Material>, 93> materials = [] {
        std::array<std::optional<Material>, 93> out;
        for (int z = 1; z <= 92; ++z)
            out[z].emplace("element", 1.0,
                           std::vector<MaterialComponent>{{static_cast<std::uint8_t>(z), 1.0}});
        return out;
    }();
    return ElectronCsda::instance().impact_xray_yield(*materials[Z], kinetic_energy_keV);
}

std::array<double, 4> element_photoelectric_vacancy_shares(int Z, double energy_keV)
{
    std::array<double, 4> out{};
    const auto& xs = CrossSectionData::instance();
    const FluorescenceData* k = xs.fluorescence(Z);
    if (k && energy_keV > k->k_edge_keV) {
        const double total = xs.sigma_photoelectric(Z, energy_keV * 1e-3);
        if (total > 0.0)
            out[0] = std::min(1.0, xs.sigma_K_photoelectric(Z, energy_keV * 1e-3) / total);
    }
    if (const LFluorescenceData* l = xs.l_fluorescence(Z)) {
        double shares[3];
        if (l_photoelectric_shares(*l, energy_keV, shares))
            for (int s = 0; s < 3; ++s) out[1 + s] = shares[s];
    }
    return out;
}

std::array<std::pair<double, double>, 3> element_photoelectron_energies(
    int Z, double energy_keV)
{
    const auto share = element_photoelectric_vacancy_shares(Z, energy_keV);
    const auto& xs = CrossSectionData::instance();
    const FluorescenceData* k = xs.fluorescence(Z);
    const LFluorescenceData* l = xs.l_fluorescence(Z);
    std::array<std::pair<double, double>, 3> out{};
    out[0] = {k ? energy_keV - k->k_edge_keV : energy_keV, share[0]};
    double l_share = 0.0, l_binding = 0.0;
    if (l)
        for (int s = 0; s < 3; ++s) {
            l_share += share[1 + s];
            l_binding += share[1 + s] * l->sub[s].edge_keV;
        }
    out[1] = {l_share > 0.0 ? energy_keV - l_binding / l_share : energy_keV, l_share};
    out[2] = {energy_keV, std::max(0.0, 1.0 - share[0] - l_share)};
    return out;
}

namespace {

// Serialization order of the charged-particle payload components.
template <typename Payload, typename Visitor>
void visit_charged_components(Payload& payload, Visitor&& visit) {
    visit("electron_brems", payload.electron_brems_spline);
    visit("positron_brems", payload.positron_brems_spline);
    visit("positron_containment", payload.positron_containment_spline);
    visit("positron_511_atom", payload.positron_511_atom_spline);
    visit("positron_annihilation_continuum",
          payload.positron_annihilation_continuum_spline);
    visit("positron_annihilation_in_flight_direct",
          payload.positron_annihilation_in_flight_direct_spline);
    visit("positron_annihilation_in_flight_continuum",
          payload.positron_annihilation_in_flight_continuum_spline);
    visit("electron_fluorescence_k_alpha", payload.electron_fluorescence_k_alpha_spline);
    visit("electron_fluorescence_k_beta", payload.electron_fluorescence_k_beta_spline);
    visit("electron_fluorescence_l", payload.electron_fluorescence_l_spline);
    visit("positron_fluorescence_k_alpha", payload.positron_fluorescence_k_alpha_spline);
    visit("positron_fluorescence_k_beta", payload.positron_fluorescence_k_beta_spline);
    visit("positron_fluorescence_l", payload.positron_fluorescence_l_spline);
}

} // namespace

std::string ChargedParticleSourceFluxModel::serialize() const
{
    std::ostringstream out;
    out << "CEELO_CHARGED_PARTICLE_SOURCE "
        << payload_.model_version << '\n';
    out << "validity "
        << std::setprecision(17)
        << payload_.min_kinetic_energy_keV << ' '
        << payload_.max_kinetic_energy_keV << ' '
        << payload_.min_photon_energy_keV << ' '
        << payload_.min_effective_Z << ' '
        << payload_.max_effective_Z << ' '
        << payload_.min_areal_density_g_cm2 << ' '
        << payload_.max_areal_density_g_cm2 << '\n';
    out << "features "
        << std::setprecision(17)
        << payload_.z_blend_width << ' '
        << payload_.range_feature_scale << ' '
        << payload_.photon_absorption_path_fraction << '\n';
    visit_charged_components(
        payload_,
        [&](const char* name, const ChargedParticleSplineField& spline) {
            if (!spline.enabled()) {
                out << name << " none\n";
                return;
            }
            out << name << " spline\n";
            write_spline_field(out, name, spline);
        });
    out << "END\n";
    return out.str();
}

ChargedParticleSourceFluxModel
ChargedParticleSourceFluxModel::deserialize(const std::string& text)
{
    std::istringstream in(text);
    ChargedParticleSourcePayload payload;
    std::string key;
    if (!(in >> key >> payload.model_version)
        || key != "CEELO_CHARGED_PARTICLE_SOURCE"
        || payload.model_version != 1) {
        throw std::runtime_error(
            "Unsupported charged-particle source payload (expected version 1)");
    }
    if (!(in >> key)
        || key != "validity"
        || !(in >> payload.min_kinetic_energy_keV
                >> payload.max_kinetic_energy_keV
                >> payload.min_photon_energy_keV
                >> payload.min_effective_Z
                >> payload.max_effective_Z
                >> payload.min_areal_density_g_cm2
                >> payload.max_areal_density_g_cm2)) {
        throw std::runtime_error(
            "Invalid charged-particle source validity");
    }
    if (!(in >> key)
        || key != "features"
        || !(in >> payload.z_blend_width
                >> payload.range_feature_scale
                >> payload.photon_absorption_path_fraction)) {
        throw std::runtime_error(
            "Invalid charged-particle source features");
    }
    visit_charged_components(
        payload,
        [&](const char* name, ChargedParticleSplineField& spline) {
            std::string tag_name, tag;
            if (!(in >> tag_name >> tag) || tag_name != name
                || (tag != "spline" && tag != "none")) {
                throw std::runtime_error(
                    "Missing charged-particle component record for "
                    + std::string(name));
            }
            if (tag == "spline")
                read_spline_field(in, name, spline);
        });
    if (!(in >> key) || key != "END")
        throw std::runtime_error(
            "Missing END in charged-particle source payload");
    return ChargedParticleSourceFluxModel(std::move(payload));
}

DetectorModelPayload nai3x3_payload() {
    DetectorModelPayload payload;
    payload.zero_deposit.basis_size = 1;
    payload.full_energy.basis_size = 1;
    payload.pair_fraction.basis_size = 1;
    payload.pair_full_energy.basis_size = 1;
    payload.single_escape.basis_size = 1;
    payload.double_escape.basis_size = 1;
    payload.annihilation_511.basis_size = 1;
    return payload;
}

DetectorModelPayload coaxial_hpge_payload() {
    DetectorModelPayload payload = nai3x3_payload();
    payload.family = DetectorFamily::CoaxialHPGe;
    payload.radius_cm = 4.0;
    payload.length_cm = 8.0;
    payload.density_g_cm3 = 5.323;
    payload.back_incidence_scale = -0.03;
    payload.fep_base = 0.10;
    payload.fep_photoelectric_scale = 3.0;
    payload.fep_reinteraction_scale = 1.35;
    payload.single_escape_scale = 0.29;
    payload.double_escape_scale = 0.12;
    payload.resolution_half_window_channels = 48;
    payload.single_compton_shape =
        {0.18, 0.52, 0.92, 1.20, 1.13, 0.73, 0.31, 0.06};
    payload.multiple_compton_shape =
        {0.91, 1.18, 1.22, 1.02, 0.72, 0.43, 0.19, 0.05};
    return payload;
}

namespace {

// "<key> N" then N rows "E delta_full R_full delta_detected R_detected"
// (FiniteDistanceDepth).
void write_finite_distance(std::ostream& out, const char* key,
                           const FiniteDistanceDepth& depth)
{
    out << key << ' ' << depth.energies_keV.size() << '\n'
        << std::setprecision(17);
    for (std::size_t i = 0; i < depth.energies_keV.size(); ++i)
        out << depth.energies_keV[i] << ' ' << depth.full_energy_cm[i] << ' '
            << depth.full_energy_radius_cm[i] << ' ' << depth.detected_cm[i]
            << ' ' << depth.detected_radius_cm[i] << '\n';
}

void read_finite_distance(std::istream& in, FiniteDistanceDepth& depth)
{
    std::size_t count = 0;
    if (!(in >> count) || count > 4096)
        throw std::runtime_error("invalid finite_distance record");
    depth = FiniteDistanceDepth{};
    for (auto* column : {&depth.energies_keV, &depth.full_energy_cm,
                         &depth.full_energy_radius_cm, &depth.detected_cm,
                         &depth.detected_radius_cm})
        column->resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        if (!(in >> depth.energies_keV[i] >> depth.full_energy_cm[i]
                 >> depth.full_energy_radius_cm[i] >> depth.detected_cm[i]
                 >> depth.detected_radius_cm[i]))
            throw std::runtime_error("invalid finite_distance row");
        if (!(depth.energies_keV[i] > 0.0)
            || (i > 0 && !(depth.energies_keV[i] > depth.energies_keV[i - 1]))
            || !(depth.full_energy_cm[i] >= 0.0) || !(depth.detected_cm[i] >= 0.0)
            || !(depth.full_energy_radius_cm[i] >= 0.0)
            || !(depth.detected_radius_cm[i] >= 0.0))
            throw std::runtime_error(
                "finite_distance needs ascending energies and nonnegative "
                "depths and radii");
    }
}

} // namespace

std::string DetectorPulseResponseModel::serialize() const {
    std::ostringstream out;
    out << "CEELO_DETECTOR_RESPONSE " << kDetectorPayloadVersion << '\n';
    out << "family "
        << (payload_.family == DetectorFamily::CoaxialHPGe
                ? "hpge_coax" : "nai3x3")
        << '\n';
    if (payload_.geometry) {
        out << "detector " << payload_.detector_name << '\n';
        write_geometry_block(out, "geometry", *payload_.geometry);
    }
    out << "validity " << std::setprecision(17)
        << payload_.min_energy_keV << ' ' << payload_.max_energy_keV << ' '
        << payload_.min_cos_theta << ' ' << payload_.max_cos_theta << '\n';
    out << "feature_map log_energy_angle_passive_pair_threshold\n";
    out << "passive_landmark_keV " << std::setprecision(17)
        << payload_.passive_landmark_keV << '\n';
    out << "basis_size " << kDetectorResponseBasisSize << '\n';
    out << "rate_surface_enabled "
        << payload_.detected_rate_surface.enabled << '\n';
    if (payload_.detected_rate_surface.enabled) {
        const std::size_t energy_size =
            payload_.detected_rate_surface.energy_size;
        const std::size_t angle_size =
            payload_.detected_rate_surface.angle_size;
        out << "rate_surface "
            << energy_size << ' ' << angle_size << ' '
            << std::setprecision(17)
            << payload_.detected_rate_surface.energy_warp_alpha;
        out << ' '
            << payload_.detected_rate_surface.angle_warp_exponent << ' '
            << payload_.detected_rate_surface.angle_warp_epsilon << '\n';
        out << "rate_surface_coefficients";
        for (std::size_t index = 0; index < energy_size * angle_size; ++index) {
            out << ' ' << std::setprecision(17)
                << payload_.detected_rate_surface
                       .log_odds_coefficients[index];
        }
        out << '\n';
        const auto& anchor = payload_.detected_rate_surface.anchor;
        out << "rate_surface_anchor_enabled " << anchor.enabled << '\n';
        if (anchor.enabled) {
            // Its geometry is the payload's.
            out << "rate_surface_anchor " << std::setprecision(17)
                << anchor.blend_keV << ' ' << anchor.blend_power << ' '
                << anchor.azimuth_quadrature_points << '\n';
        }
    }
    const auto write_component =
        [&](const char* name,
            const DetectorResponseRbfComponent& component) {
            write_detector_rbf_component(out, name, component);
        };
    write_component("zero_deposit", payload_.zero_deposit);
    write_component("full_energy", payload_.full_energy);
    write_component("pair_fraction", payload_.pair_fraction);
    write_component(
        "pair_full_energy", payload_.pair_full_energy);

    write_component("single_escape", payload_.single_escape);
    write_component("double_escape", payload_.double_escape);
    write_component("annihilation_511", payload_.annihilation_511);
    out << "k_escape_lines " << std::setprecision(17)
        << payload_.k_escape_alpha_energy_keV << ' '
        << payload_.k_escape_beta_energy_keV << ' '
        << payload_.k_escape_edge_keV << '\n';
    out << "k_escape_support" << std::setprecision(17);
    for (double limit : payload_.k_escape_support_keV)
        out << ' ' << limit;
    out << '\n';
    write_component("k_alpha_escape", payload_.k_alpha_escape);
    write_component("k_beta_escape", payload_.k_beta_escape);
    write_component("k_alpha_absorbed", payload_.k_alpha_absorbed);
    write_component("k_beta_absorbed", payload_.k_beta_absorbed);
    write_component(
        "photoelectric_escape",
        payload_.detector_photoelectric_escape);
    write_component(
        "single_compton", payload_.detector_single_compton);
    write_component(
        "multiple_compton", payload_.detector_multiple_compton);
    write_component(
        "pair_no_compton", payload_.detector_pair_no_compton);
    write_component(
        "pair_single_compton",
        payload_.detector_pair_single_compton);
    write_component(
        "pair_multiple_compton",
        payload_.detector_pair_multiple_compton);
    write_finite_distance(out, "finite_distance", payload_.finite_distance);
    out << "END\n";
    return out.str();
}

DetectorPulseResponseModel DetectorPulseResponseModel::deserialize(
    const std::string& text)
{
    std::istringstream in(text);
    std::string key;
    std::uint32_t version = 0;
    if (!(in >> key >> version) || key != "CEELO_DETECTOR_RESPONSE"
        || version != kDetectorPayloadVersion) {
        throw std::runtime_error(
            "Unsupported detector-response payload header (expected version "
            + std::to_string(kDetectorPayloadVersion) + ")");
    }
    if (!(in >> key) || key != "family")
        throw std::runtime_error("Missing detector-response family");
    std::string family;
    if (!(in >> family))
        throw std::runtime_error("Truncated detector-response family");
    DetectorModelPayload payload =
        family == "hpge_coax" ? coaxial_hpge_payload() : nai3x3_payload();
    if (family != "hpge_coax" && family != "nai3x3")
        throw std::runtime_error(
            "Unsupported detector-response family " + family);
    if (!(in >> key))
        throw std::runtime_error("Truncated detector-response payload");
    if (key == "detector") {
        if (!(in >> payload.detector_name) || !(in >> key) || key != "geometry")
            throw std::runtime_error(
                "A detector-response label must precede its geometry");
        payload.geometry = read_geometry_block(in);
        if (!(in >> key))
            throw std::runtime_error("Truncated detector-response payload");
    }
    if (key != "validity"
        || !(in >> payload.min_energy_keV >> payload.max_energy_keV
                >> payload.min_cos_theta >> payload.max_cos_theta)) {
        throw std::runtime_error(
            "Missing detector-response validity range");
    }
    {
        std::string feature_map;
        if (!(in >> key >> feature_map)
            || key != "feature_map"
            || feature_map != "log_energy_angle_passive_pair_threshold") {
            throw std::runtime_error(
                "Unsupported detector-response feature map");
        }
        if (!(in >> key >> payload.passive_landmark_keV)
            || key != "passive_landmark_keV"
            || !(payload.passive_landmark_keV > 0.0)) {
            throw std::runtime_error(
                "Missing detector-response passive landmark");
        }
    }
    std::size_t basis_size = 0;
    if (!(in >> key >> basis_size) || key != "basis_size"
        || basis_size < 1 || basis_size > kDetectorResponseBasisSize) {
        throw std::runtime_error(
            "Unsupported detector-response basis size");
    }
    bool read_rate_surface = false;
    if (!(in >> key >> read_rate_surface)
        || key != "rate_surface_enabled") {
        throw std::runtime_error(
            "Missing detector rate-surface enable flag");
    }

    if (read_rate_surface) {
        std::size_t energy_basis_size = 0;
        std::size_t angle_basis_size = 0;
        if (!(in >> key >> energy_basis_size >> angle_basis_size
                >> payload.detected_rate_surface.energy_warp_alpha)
            || key != "rate_surface"
            || energy_basis_size < 2
            || energy_basis_size > kDetectorRateMaxEnergyBasisSize
            || angle_basis_size < 2
            || angle_basis_size > kDetectorRateMaxAngleBasisSize) {
            throw std::runtime_error(
                "Unsupported detector rate surface");
        }
        payload.detected_rate_surface.energy_size =
            static_cast<std::uint32_t>(energy_basis_size);
        payload.detected_rate_surface.angle_size =
            static_cast<std::uint32_t>(angle_basis_size);
        if (!(in >> payload.detected_rate_surface.angle_warp_exponent
                >> payload.detected_rate_surface.angle_warp_epsilon)) {
            throw std::runtime_error(
                "Missing detector rate-surface angle warp");
        }
        if (!(in >> key) || key != "rate_surface_coefficients")
            throw std::runtime_error(
                "Missing detector rate-surface coefficients");
        for (std::size_t index = 0;
             index < energy_basis_size * angle_basis_size; ++index) {
            if (!(in >> payload.detected_rate_surface
                          .log_odds_coefficients[index])) {
                throw std::runtime_error(
                    "Truncated detector rate-surface coefficients");
            }
        }
        payload.detected_rate_surface.enabled = true;
        auto& anchor = payload.detected_rate_surface.anchor;
        bool anchor_enabled = false;
        if (!(in >> key >> anchor_enabled)
            || key != "rate_surface_anchor_enabled") {
            throw std::runtime_error(
                "Missing detector rate-surface anchor flag");
        }
        if (anchor_enabled) {
            if (!(in >> key >> anchor.blend_keV >> anchor.blend_power
                     >> anchor.azimuth_quadrature_points)
                || key != "rate_surface_anchor") {
                throw std::runtime_error(
                    "Invalid detector rate-surface anchor");
            }
            if (!payload.geometry)
                throw std::runtime_error(
                    "A detector rate-surface anchor needs the detector geometry");
            anchor.enabled = true;
        }
    }
    const auto read_component =
        [&](const char* name,
            DetectorResponseRbfComponent& component) {
            read_detector_rbf_component(in, name, component);
        };
    read_component("zero_deposit", payload.zero_deposit);
    read_component("full_energy", payload.full_energy);
    read_component("pair_fraction", payload.pair_fraction);
    read_component(
        "pair_full_energy", payload.pair_full_energy);

    read_component("single_escape", payload.single_escape);
    read_component("double_escape", payload.double_escape);
    read_component(
        "annihilation_511", payload.annihilation_511);
    if (!(in >> key >> payload.k_escape_alpha_energy_keV
             >> payload.k_escape_beta_energy_keV
             >> payload.k_escape_edge_keV)
        || key != "k_escape_lines") {
        throw std::runtime_error(
            "Missing detector K-escape line energies");
    }
    // Upper end of each K component's fitted support (zero disables the
    // gate for that component); see DetectorModelPayload.
    if (!(in >> key) || key != "k_escape_support")
        throw std::runtime_error(
            "Missing detector K-escape support record");
    for (double& limit : payload.k_escape_support_keV) {
        if (!(in >> limit) || !std::isfinite(limit) || limit < 0.0)
            throw std::runtime_error(
                "Invalid detector K-escape support limit");
    }
    read_component("k_alpha_escape", payload.k_alpha_escape);
    read_component("k_beta_escape", payload.k_beta_escape);
    read_component("k_alpha_absorbed", payload.k_alpha_absorbed);
    read_component("k_beta_absorbed", payload.k_beta_absorbed);
    read_component(
        "photoelectric_escape",
        payload.detector_photoelectric_escape);

    read_component(
        "single_compton", payload.detector_single_compton);
    read_component(
        "multiple_compton", payload.detector_multiple_compton);
    read_component(
        "pair_no_compton",
        payload.detector_pair_no_compton);
    read_component(
        "pair_single_compton",
        payload.detector_pair_single_compton);
    read_component(
        "pair_multiple_compton",
        payload.detector_pair_multiple_compton);
    if (!(in >> key) || key != "finite_distance")
        throw std::runtime_error(
            "Missing detector-response finite_distance record");
    read_finite_distance(in, payload.finite_distance);
    if (!(in >> key) || key != "END")
        throw std::runtime_error("Missing detector-response END marker");
    apply_payload_geometry(payload);
    if (payload.detected_rate_surface.anchor.enabled)
        validate_rate_surface_anchor(payload.detected_rate_surface.anchor);
    payload.use_response_surrogate =
        payload.detected_rate_surface.enabled
        || !payload.full_energy.nodes.empty()
        || !payload.detector_single_compton.nodes.empty();
    return DetectorPulseResponseModel(std::move(payload));
}

PartialMassAttenuation DetectorPulseResponseModel::mass_cross_sections(
    double energy_keV) const
{
    if (payload_.geometry) {
        const GeometryDescriptor& geometry = *payload_.geometry;
        std::vector<std::pair<int, double>> composition;
        for (const MaterialComponent& c :
             geometry.materials.at(static_cast<std::size_t>(
                 geometry.crystal_material_index)).composition)
            composition.emplace_back(static_cast<int>(c.Z), c.mass_fraction);
        return compound_mass_cross_sections(energy_keV, composition);
    }
    if (payload_.family == DetectorFamily::CoaxialHPGe)
        return element_mass_cross_sections(32, energy_keV);

    const double a_na = CrossSectionData::instance().atomic_weight(11);
    const double a_i = CrossSectionData::instance().atomic_weight(53);
    const double denominator = a_na + a_i;
    return compound_mass_cross_sections(
        energy_keV, {{11, a_na / denominator}, {53, a_i / denominator}});
}

namespace {

// Shared source-section codec. The standalone CEELO_SOURCE_ESCAPE payload and
// the source block embedded in CEELO_FULLSPEC use identical keys, so a
// source-escape payload and a full-spectrum payload can never disagree on the
// source model's serialized form.
void write_source_section(std::ostream& out, const SourceModelPayload& s) {
    // The embedded FULLSPEC caller already runs at full precision; the
    // standalone container must not depend on ambient stream state.
    out << std::setprecision(17);
    out << "source_validity " << s.min_energy_keV << ' '
        << s.max_energy_keV << ' ' << s.min_effective_Z << ' '
        << s.max_effective_Z << ' ' << s.min_areal_density_g_cm2 << ' '
        << s.max_areal_density_g_cm2 << '\n';
    out << "source_scalars " << s.model_version << ' ' << s.z_blend_width
        << ' ' << s.single_escape_scale << ' ' << s.double_escape_scale
        << ' ' << s.multiple_escape_scale
        << ' ' << s.compton_rayleigh_scale
        << ' ' << s.double_rayleigh_scale
        << ' ' << s.multiple_rayleigh_scale
        << ' ' << s.elastic_escape_scale << ' ' << s.brems_yield_scale
        << ' ' << s.enable_pair_annihilation
        << ' ' << s.enable_source_electron_brems << '\n';
    out << "source_elastic_damping " << std::setprecision(17)
        << s.elastic_compton_damping << '\n';
    out << "source_rayleigh " << s.enable_rayleigh << '\n';
    out << "source_fluorescence " << s.enable_fluorescence << '\n';
    out << "source_compton_doppler "
        << s.enable_compton_doppler << '\n';
    write_array(out, "source_single", s.single_compton_shape);
    write_array(out, "source_multiple", s.multiple_compton_shape);
    write_array(out, "source_compton_rayleigh",
                s.compton_rayleigh_shape);
    write_array(out, "source_brems", s.brems_shape);
    out << "source_multiple_surrogate_version "
        << s.multiple_scatter_surrogate.model_version << '\n';
    out << "source_multiple_surrogate_options "
        << s.multiple_scatter_surrogate.inverse_energy_warp_scale << ' '
        << s.multiple_scatter_surrogate.absorption_envelope_strength << ' '
        << s.multiple_scatter_surrogate.absorption_envelope_power << ' '
        << s.multiple_scatter_surrogate.absorption_envelope_max_ratio << ' '
        << s.multiple_scatter_surrogate.use_cubic_bspline << ' '
        << s.multiple_scatter_surrogate.merge_higher_rayleigh << '\n';
    const auto& higher_yield =
        s.multiple_scatter_surrogate.higher_yield;
    out << "source_higher_yield " << higher_yield.degree;
    for (double bandwidth : higher_yield.bandwidth)
        out << ' ' << std::setprecision(17) << bandwidth;
    out << ' ' << higher_yield.minimum_tau
        << ' ' << higher_yield.maximum_tau
        << ' ' << higher_yield.nodes.size() << '\n';
    for (const auto& node : higher_yield.nodes) {
        out << "source_higher_yield_node";
        for (double value : node.feature)
            out << ' ' << std::setprecision(17) << value;
        for (double value : node.coefficients)
            out << ' ' << std::setprecision(17) << value;
        out << '\n';
    }
    write_rbf_component(
        out, "source_higher_rbf",
        s.multiple_scatter_surrogate.higher_compton);
    write_rbf_component(
        out, "source_higher_rayleigh_rbf",
        s.multiple_scatter_surrogate.higher_compton_rayleigh);
    out << "source_secondary_surrogate_version "
        << s.secondary_surrogate.model_version << '\n';
    out << "source_secondary_surrogate_options "
        << s.secondary_surrogate.inverse_energy_warp_scale << ' '
        << s.secondary_surrogate.absorption_envelope_strength << ' '
        << s.secondary_surrogate.absorption_envelope_power << ' '
        << s.secondary_surrogate.absorption_envelope_max_ratio << ' '
        << s.secondary_surrogate.use_cubic_bspline;
    // All six components share one feature width and one local-linear ridge,
    // so they travel here rather than on every component line.
    out << ' ' << s.secondary_surrogate.pair_511_atom.feature_dimension
        << ' ' << std::setprecision(17)
        << s.secondary_surrogate.pair_511_atom.local_linear_ridge
        << ' '
        << s.secondary_surrogate.pair_511_atom.local_linear_feature_index
        << ' ' << std::setprecision(17)
        << s.secondary_surrogate.tail_decay_kappa;
    out << '\n';
    {
        const auto components = s.secondary_surrogate.components();
        for (std::size_t c = 0; c < components.size(); ++c)
            write_rbf_component(out, kSecondaryComponentKeys[c], *components[c]);
    }
    out << "source_photonuclear_options "
        << s.photonuclear.model_version << ' '
        << s.enable_photonuclear << ' '
        << std::setprecision(17)
        << s.photonuclear.escaped_yield_scale << '\n';
    write_array(
        out, "source_photonuclear_shape",
        s.photonuclear.shape);
    if (!s.photonuclear.secondary.nodes.empty()) {
        // Written only when the photonuclear surrogate is fitted.
        out << "source_photonuclear_surrogate_options "
            << std::setprecision(17)
            << s.photonuclear.inverse_energy_warp_scale << ' '
            << s.photonuclear.use_cubic_bspline << ' '
            << s.photonuclear.escape_kappa << '\n';
        write_rbf_component(
            out, "source_photonuclear_rbf", s.photonuclear.secondary);
    }
    out << "source_photonuclear_table "
        << s.photonuclear.cross_sections.size() << '\n';
    for (const auto& row : s.photonuclear.cross_sections) {
        out << "source_photonuclear_row "
            << std::setprecision(17) << row.energy_keV;
        for (double value : row.mass_attenuation_cm2_g)
            out << ' ' << std::setprecision(17) << value;
        out << '\n';
    }
}

bool read_source_key(
    const std::string& key, std::istream& in, SourceModelPayload& source)
{
    if (key == "source_validity") {
        in >> source.min_energy_keV
           >> source.max_energy_keV
           >> source.min_effective_Z
           >> source.max_effective_Z
           >> source.min_areal_density_g_cm2
           >> source.max_areal_density_g_cm2;
    } else if (key == "source_elastic_damping") {
        in >> source.elastic_compton_damping;
    } else if (key == "source_scalars") {
        in >> source.model_version
           >> source.z_blend_width
           >> source.single_escape_scale
           >> source.double_escape_scale
           >> source.multiple_escape_scale
           >> source.compton_rayleigh_scale
           >> source.double_rayleigh_scale
           >> source.multiple_rayleigh_scale
           >> source.elastic_escape_scale
           >> source.brems_yield_scale
           >> source.enable_pair_annihilation
           >> source.enable_source_electron_brems;
        if (in && source.model_version != 1)
            throw std::runtime_error(
                "Unsupported source-escape model version (expected 1)");
    } else if (key == "source_rayleigh") {
        in >> source.enable_rayleigh;
    } else if (key == "source_fluorescence") {
        in >> source.enable_fluorescence;
    } else if (key == "source_compton_doppler") {
        in >> source.enable_compton_doppler;
    } else if (key == "source_single") {
        read_array(in, source.single_compton_shape);
    } else if (key == "source_multiple") {
        read_array(in, source.multiple_compton_shape);
    } else if (key == "source_compton_rayleigh") {
        read_array(in, source.compton_rayleigh_shape);
    } else if (key == "source_brems") {
        read_array(in, source.brems_shape);
    } else if (key == "source_multiple_surrogate_version") {
        in >> source.multiple_scatter_surrogate.model_version;
        if (in && source.multiple_scatter_surrogate.model_version != 1)
            throw std::runtime_error(
                "Unsupported multiple-scatter surrogate version (expected 1)");
    } else if (key == "source_multiple_surrogate_options") {
        auto& multiple = source.multiple_scatter_surrogate;
        in >> multiple.inverse_energy_warp_scale
           >> multiple.absorption_envelope_strength
           >> multiple.absorption_envelope_power
           >> multiple.absorption_envelope_max_ratio
           >> multiple.use_cubic_bspline
           >> multiple.merge_higher_rayleigh;
    } else if (key == "source_higher_yield") {
        std::size_t count = 0;
        auto& yield_model =
            source.multiple_scatter_surrogate.higher_yield;
        in >> yield_model.degree;
        if (yield_model.degree + 1
            > kMultipleScatterYieldBasisSize) {
            throw std::runtime_error(
                "Unsupported higher-yield polynomial degree");
        }
        read_array(in, yield_model.bandwidth);
        in >> yield_model.minimum_tau
           >> yield_model.maximum_tau >> count;
        yield_model.nodes.reserve(count);
    } else if (key == "source_higher_yield_node") {
        MultipleScatterYieldNode node;
        read_array(in, node.feature);
        read_array(in, node.coefficients);
        source.multiple_scatter_surrogate
            .higher_yield.nodes.push_back(std::move(node));
    } else if (key == "source_higher_rbf") {
        read_rbf_header(in, source.multiple_scatter_surrogate.higher_compton);
    } else if (key == "source_higher_rbf_node") {
        read_rbf_node(in, source.multiple_scatter_surrogate.higher_compton);
    } else if (key == "source_higher_rayleigh_rbf") {
        read_rbf_header(in, source.multiple_scatter_surrogate.higher_compton_rayleigh);
    } else if (key == "source_higher_rayleigh_rbf_node") {
        read_rbf_node(in, source.multiple_scatter_surrogate.higher_compton_rayleigh);
    } else if (key == "source_secondary_surrogate_version") {
        in >> source.secondary_surrogate.model_version;
        if (in && source.secondary_surrogate.model_version != 1)
            throw std::runtime_error(
                "Unsupported secondary surrogate version (expected 1)");
    } else if (key == "source_secondary_surrogate_options") {
        auto& secondary = source.secondary_surrogate;
        std::uint32_t feature_dimension = 3;
        double local_linear_ridge = 0.0;
        std::uint32_t local_linear_feature_index = 2;
        in >> secondary.inverse_energy_warp_scale
           >> secondary.absorption_envelope_strength
           >> secondary.absorption_envelope_power
           >> secondary.absorption_envelope_max_ratio
           >> secondary.use_cubic_bspline
           >> feature_dimension >> local_linear_ridge
           >> local_linear_feature_index
           >> secondary.tail_decay_kappa;
        if (in
            && (feature_dimension < 3
                || feature_dimension > kSourceRbfMaxFeatures)) {
            throw std::runtime_error(
                "Unsupported source secondary feature dimension");
        }
        // Applied before any component is read, so read_rbf_header /
        // read_rbf_node parse the right width.
        for (auto* component : secondary.components()) {
            component->feature_dimension = feature_dimension;
            component->local_linear_ridge = local_linear_ridge;
            component->local_linear_feature_index =
                local_linear_feature_index;
            component->tail_decay_kappa = secondary.tail_decay_kappa;
        }
    } else if (MultipleScatterRbfComponent* component =
                   secondary_component_for_key(source.secondary_surrogate, key)) {
        read_rbf_header(in, *component);
    } else if (MultipleScatterRbfComponent* component =
                   secondary_component_for_node_key(source.secondary_surrogate, key)) {
        read_rbf_node(in, *component);
    } else if (key == "source_photonuclear_options") {
        in >> source.photonuclear.model_version
           >> source.enable_photonuclear
           >> source.photonuclear.escaped_yield_scale;
        if (in && source.photonuclear.model_version != 1)
            throw std::runtime_error(
                "Unsupported photonuclear model version (expected 1)");
    } else if (key == "source_photonuclear_shape") {
        read_array(in, source.photonuclear.shape);
    } else if (key == "source_photonuclear_surrogate_options") {
        in >> source.photonuclear.inverse_energy_warp_scale
           >> source.photonuclear.use_cubic_bspline
           >> source.photonuclear.escape_kappa;
    } else if (key == "source_photonuclear_rbf") {
        read_rbf_header(in, source.photonuclear.secondary);
    } else if (key == "source_photonuclear_rbf_node") {
        read_rbf_node(in, source.photonuclear.secondary);
    } else if (key == "source_photonuclear_table") {
        std::size_t count = 0;
        in >> count;
        source.photonuclear.cross_sections.reserve(count);
    } else if (key == "source_photonuclear_row") {
        PhotonuclearCrossSectionRow row;
        in >> row.energy_keV;
        read_array(in, row.mass_attenuation_cm2_g);
        source.photonuclear.cross_sections.push_back(
            std::move(row));
    } else {
        return false;
    }
    return true;
}

} // namespace

std::string SourceGeomFluxModel::serialize() const {
    std::ostringstream out;
    out << "CEELO_SOURCE_ESCAPE 1\n";
    write_source_section(out, payload_);
    out << "END\n";
    return out.str();
}

// Derive each secondary component's deepest populated optical depth from its
// own node cloud, so the tail-decay envelope can never disagree with the nodes
// actually stored. A no-op when tail_decay_kappa is zero.
void finalize_source_payload(SourceModelPayload& payload) {
    for (auto* component : payload.secondary_surrogate.components())
        component->refresh_tail_decay_tau();
}

SourceGeomFluxModel SourceGeomFluxModel::deserialize(const std::string& text) {
    std::istringstream in(text);
    std::string key;
    std::uint32_t file_version = 0;
    if (!(in >> key >> file_version) || key != "CEELO_SOURCE_ESCAPE")
        throw std::runtime_error("Not a CeeLo source-escape payload");
    if (file_version != 1)
        throw std::runtime_error(
            "Unsupported source-escape payload version");
    SourceModelPayload payload;
    bool saw_end = false;
    while (in >> key) {
        if (key == "END") {
            saw_end = true;
            break;
        }
        if (!read_source_key(key, in, payload))
            throw std::runtime_error(
                "Unknown source-escape payload field: " + key);
        if (!in)
            throw std::runtime_error("Malformed source-escape payload");
    }
    if (!saw_end)
        throw std::runtime_error("Truncated source-escape payload");
    finalize_source_payload(payload);
    return SourceGeomFluxModel(std::move(payload));
}


std::string FastSpectrumModel::serialize() const {
    std::ostringstream out;
    out << "CEELO_FULLSPEC " << payload_.model_version << '\n';
    out << "model_name " << std::quoted(payload_.model_name) << '\n';
    out << "validity " << std::quoted(payload_.validity) << '\n';
    out << "training_hash " << std::quoted(payload_.training_hash) << '\n';
    out << "statistical_uncertainty "
        << std::quoted(payload_.statistical_uncertainty) << '\n';
    out << "physics_toggles " << std::quoted(payload_.physics_toggles) << '\n';
    out << "distance_cm " << std::setprecision(17)
        << payload_.distance_cm << '\n';
    out << "incidence_cos_theta " << payload_.incidence_cos_theta << '\n';
    out << "source_quadrature_bins " << payload_.source_quadrature_bins << '\n';

    write_source_section(out, payload_.source);

    const auto& d = payload_.detector;
    out << "detector_validity " << d.min_energy_keV << ' '
        << d.max_energy_keV << ' ' << d.min_cos_theta << ' '
        << d.max_cos_theta << '\n';
    out << "detector_scalars " << d.model_version << ' '
        << static_cast<int>(d.family) << ' ' << d.radius_cm << ' '
        << d.length_cm << ' ' << d.density_g_cm3 << ' '
        << d.back_incidence_scale << ' ' << d.fep_base << ' '
        << d.fep_photoelectric_scale << ' ' << d.fep_reinteraction_scale
        << ' ' << d.single_escape_scale << ' ' << d.double_escape_scale
        << ' ' << d.resolution_half_window_channels << '\n';
    if (d.geometry) {
        if (d.detector_name.empty()
            || d.detector_name.find_first_of(" \t\n") != std::string::npos)
            throw std::runtime_error(
                "a detector with a geometry needs a one-word label");
        out << "detector_name " << d.detector_name << '\n';
        write_geometry_block(out, "detector_geometry", *d.geometry);
    }
    out << "detector_passive_landmark_keV " << std::setprecision(17)
        << d.passive_landmark_keV << '\n';
    write_finite_distance(out, "detector_finite_distance", d.finite_distance);
    write_array(out, "detector_single", d.single_compton_shape);
    write_array(out, "detector_multiple", d.multiple_compton_shape);
    out << "detector_surrogate_enabled "
        << d.use_response_surrogate << '\n';
    if (d.use_response_surrogate) {
        out << "detector_rate_surface_enabled "
            << d.detected_rate_surface.enabled << '\n';
        if (d.detected_rate_surface.enabled) {
            out << "detector_rate_surface "
                << d.detected_rate_surface.energy_size << ' '
                << d.detected_rate_surface.angle_size << ' '
                << std::setprecision(17)
                << d.detected_rate_surface.energy_warp_alpha << ' '
                << d.detected_rate_surface.angle_warp_exponent << ' '
                << d.detected_rate_surface.angle_warp_epsilon << '\n';
            out << "detector_rate_surface_coefficients";
            for (std::size_t index = 0;
                 index < std::size_t(d.detected_rate_surface.energy_size)
                             * d.detected_rate_surface.angle_size;
                 ++index) {
                out << ' ' << std::setprecision(17)
                    << d.detected_rate_surface.log_odds_coefficients[index];
            }
            out << '\n';
            // Written only when the surface carries an anchor.
            const auto& anchor = d.detected_rate_surface.anchor;
            if (anchor.enabled) {
                // Its geometry is the detector's.
                out << "detector_rate_surface_anchor "
                    << std::setprecision(17)
                    << anchor.blend_keV << ' ' << anchor.blend_power << ' '
                    << anchor.azimuth_quadrature_points << '\n';
            }
        }
        out << "detector_response_basis "
            << kDetectorResponseBasisSize << '\n';
        out << "detector_component_basis_sizes "
            << d.zero_deposit.basis_size << ' '
            << d.full_energy.basis_size << ' '
            << d.single_escape.basis_size << ' '
            << d.double_escape.basis_size << ' '
            << d.annihilation_511.basis_size << ' '
            << d.detector_photoelectric_escape.basis_size << ' '
            << d.detector_single_compton.basis_size << ' '
            << d.detector_multiple_compton.basis_size << '\n';
        // Block partitions are written only for components that use them.
        {
            const std::array<
                std::pair<const char*, const DetectorResponseRbfComponent*>, 5>
                blocked{{
                    {"photoelectric_escape",
                     &d.detector_photoelectric_escape},
                    {"single_compton", &d.detector_single_compton},
                    {"multiple_compton", &d.detector_multiple_compton},
                    {"pair_single_compton",
                     &d.detector_pair_single_compton},
                    {"pair_multiple_compton",
                     &d.detector_pair_multiple_compton}}};
            for (const auto& entry : blocked) {
                if (!entry.second->uses_blocks())
                    continue;
                out << "detector_component_blocks " << entry.first;
                for (std::uint32_t size : entry.second->block_sizes)
                    if (size != 0)
                        out << ' ' << size;
                out << " 0\n";
            }
            if (d.detector_pair_no_compton.uses_blocks()) {
                out << "detector_component_blocks pair_no_compton";
                for (std::uint32_t size :
                     d.detector_pair_no_compton.block_sizes)
                    if (size != 0)
                        out << ' ' << size;
                out << " 0\n";
            }
        }
        if (!d.pair_fraction.nodes.empty()) {
            out << "detector_pair_component_basis_sizes "
                << d.pair_fraction.basis_size << ' '
                << d.pair_full_energy.basis_size << ' '
                << d.detector_pair_no_compton.basis_size << ' '
                << d.detector_pair_single_compton.basis_size << ' '
                << d.detector_pair_multiple_compton.basis_size << '\n';
            write_rbf_component(
                out, "detector_pair_fraction_rbf",
                d.pair_fraction);
            write_rbf_component(
                out, "detector_pair_fep_rbf",
                d.pair_full_energy);
        }
        write_rbf_component(
            out, "detector_zero_rbf", d.zero_deposit);
        write_rbf_component(
            out, "detector_fep_rbf", d.full_energy);
        write_rbf_component(
            out, "detector_single_escape_rbf", d.single_escape);
        write_rbf_component(
            out, "detector_double_escape_rbf", d.double_escape);
        write_rbf_component(
            out, "detector_annihilation_511_rbf",
            d.annihilation_511);
        // K-escape lines are written only when the detector has them.
        if (d.k_escape_alpha_energy_keV > 0.0) {
            out << "detector_k_escape_lines " << std::setprecision(17)
                << d.k_escape_alpha_energy_keV << ' '
                << d.k_escape_beta_energy_keV << ' '
                << d.k_escape_edge_keV;
            for (double limit : d.k_escape_support_keV)
                out << ' ' << limit;
            out << '\n';
            write_rbf_component(
                out, "detector_k_alpha_escape_rbf", d.k_alpha_escape);
            write_rbf_component(
                out, "detector_k_beta_escape_rbf", d.k_beta_escape);
            if (!d.k_alpha_absorbed.nodes.empty()
                || !d.k_beta_absorbed.nodes.empty()) {
                write_rbf_component(
                    out, "detector_k_alpha_absorbed_rbf",
                    d.k_alpha_absorbed);
                write_rbf_component(
                    out, "detector_k_beta_absorbed_rbf",
                    d.k_beta_absorbed);
            }
        }
        write_rbf_component(
            out, "detector_photoelectric_escape_rbf",
            d.detector_photoelectric_escape);
        write_rbf_component(
            out, "detector_single_compton_rbf",
            d.detector_single_compton);
        write_rbf_component(
            out, "detector_multiple_compton_rbf",
            d.detector_multiple_compton);
        if (!d.pair_fraction.nodes.empty()) {
            write_rbf_component(
                out, "detector_pair_no_compton_rbf",
                d.detector_pair_no_compton);
            write_rbf_component(
                out, "detector_pair_single_compton_rbf",
                d.detector_pair_single_compton);
            write_rbf_component(
                out, "detector_pair_multiple_compton_rbf",
                d.detector_pair_multiple_compton);
        }
    }
    out << "coefficient_standard_errors "
        << payload_.coefficient_standard_errors.size();
    for (double error : payload_.coefficient_standard_errors)
        out << ' ' << error;
    out << '\n';
    out << "END\n";
    return out.str();
}

FastSpectrumModel FastSpectrumModel::deserialize(const std::string& text) {
    std::istringstream in(text);
    std::string key;
    FullSpectrumPayload payload;
    std::uint32_t file_version = 0;
    if (!(in >> key >> file_version) || key != "CEELO_FULLSPEC")
        throw std::runtime_error("Not a CeeLo FullSpecApprox payload");
    if (file_version != 1)
        throw std::runtime_error("Unsupported FullSpecApprox payload version");
    payload.model_version = file_version;

    bool saw_end = false;
    while (in >> key) {
        if (key == "END") {
            saw_end = true;
            break;
        } else if (key == "model_name") {
            in >> std::quoted(payload.model_name);
        } else if (key == "validity") {
            in >> std::quoted(payload.validity);
        } else if (key == "training_hash") {
            in >> std::quoted(payload.training_hash);
        } else if (key == "statistical_uncertainty") {
            in >> std::quoted(payload.statistical_uncertainty);
        } else if (key == "physics_toggles") {
            in >> std::quoted(payload.physics_toggles);
        } else if (key == "distance_cm") {
            in >> payload.distance_cm;
        } else if (key == "incidence_cos_theta") {
            in >> payload.incidence_cos_theta;
        } else if (key == "source_quadrature_bins") {
            in >> payload.source_quadrature_bins;
        } else if (read_source_key(key, in, payload.source)) {
            // Handled by the shared source-section reader.
        } else if (key == "detector_validity") {
            in >> payload.detector.min_energy_keV
               >> payload.detector.max_energy_keV
               >> payload.detector.min_cos_theta
               >> payload.detector.max_cos_theta;
        } else if (key == "detector_scalars") {
            int family = 0;
            in >> payload.detector.model_version >> family
               >> payload.detector.radius_cm >> payload.detector.length_cm
               >> payload.detector.density_g_cm3
               >> payload.detector.back_incidence_scale
               >> payload.detector.fep_base
               >> payload.detector.fep_photoelectric_scale
               >> payload.detector.fep_reinteraction_scale
               >> payload.detector.single_escape_scale
               >> payload.detector.double_escape_scale
               >> payload.detector.resolution_half_window_channels;
            payload.detector.family =
                static_cast<DetectorFamily>(family);
        } else if (key == "detector_single") {
            read_array(in, payload.detector.single_compton_shape);
        } else if (key == "detector_multiple") {
            read_array(in, payload.detector.multiple_compton_shape);
        } else if (key == "detector_surrogate_enabled") {
            in >> payload.detector.use_response_surrogate;
        } else if (key == "detector_rate_surface_enabled") {
            in >> payload.detector.detected_rate_surface.enabled;
        } else if (key == "detector_rate_surface") {
            std::size_t energy_basis_size = 0;
            std::size_t angle_basis_size = 0;
            auto& surface = payload.detector.detected_rate_surface;
            in >> energy_basis_size >> angle_basis_size
                >> surface.energy_warp_alpha
                >> surface.angle_warp_exponent
                >> surface.angle_warp_epsilon;
            if (energy_basis_size < 2
                || energy_basis_size > kDetectorRateMaxEnergyBasisSize
                || angle_basis_size < 2
                || angle_basis_size > kDetectorRateMaxAngleBasisSize) {
                throw std::runtime_error(
                    "Unsupported embedded detector rate surface");
            }
            surface.energy_size =
                static_cast<std::uint32_t>(energy_basis_size);
            surface.angle_size =
                static_cast<std::uint32_t>(angle_basis_size);
        } else if (key == "detector_rate_surface_coefficients") {
            auto& surface = payload.detector.detected_rate_surface;
            for (std::size_t index = 0;
                 index < std::size_t(surface.energy_size)
                             * surface.angle_size;
                 ++index) {
                if (!(in >> surface.log_odds_coefficients[index])) {
                    throw std::runtime_error(
                        "Truncated embedded detector rate-surface "
                        "coefficients");
                }
            }
        } else if (key == "detector_name") {
            in >> payload.detector.detector_name;
        } else if (key == "detector_geometry") {
            payload.detector.geometry = read_geometry_block(in);
        } else if (key == "detector_passive_landmark_keV") {
            in >> payload.detector.passive_landmark_keV;
        } else if (key == "detector_finite_distance") {
            read_finite_distance(in, payload.detector.finite_distance);
        } else if (key == "detector_rate_surface_anchor") {
            auto& anchor = payload.detector.detected_rate_surface.anchor;
            if (!(in >> anchor.blend_keV >> anchor.blend_power
                     >> anchor.azimuth_quadrature_points)) {
                throw std::runtime_error(
                    "Truncated embedded detector rate-surface anchor");
            }
            // Validated once the geometry has filled it in, after the loop.
            anchor.enabled = true;
        } else if (key == "detector_response_basis") {
            std::size_t basis_size = 0;
            in >> basis_size;
            if (basis_size != kDetectorResponseBasisSize) {
                throw std::runtime_error(
                    "Unsupported embedded detector-response basis size");
            }
        } else if (key == "detector_component_basis_sizes") {
            auto& detector = payload.detector;
            in >> detector.zero_deposit.basis_size
               >> detector.full_energy.basis_size
               >> detector.single_escape.basis_size
               >> detector.double_escape.basis_size
               >> detector.annihilation_511.basis_size
               >> detector.detector_photoelectric_escape.basis_size
               >> detector.detector_single_compton.basis_size
               >> detector.detector_multiple_compton.basis_size;
            const std::array<std::uint32_t, 8> component_sizes{
                detector.zero_deposit.basis_size,
                detector.full_energy.basis_size,
                detector.single_escape.basis_size,
                detector.double_escape.basis_size,
                detector.annihilation_511.basis_size,
                detector.detector_photoelectric_escape.basis_size,
                detector.detector_single_compton.basis_size,
                detector.detector_multiple_compton.basis_size};
            for (const auto component_size : component_sizes) {
                if (component_size < 1
                    || component_size > kDetectorResponseBasisSize) {
                    throw std::runtime_error(
                        "Unsupported embedded detector component basis size");
                }
            }
        } else if (key == "detector_component_blocks") {
            auto& detector = payload.detector;
            std::string component_name;
            if (!(in >> component_name))
                throw std::runtime_error(
                    "Truncated embedded detector block partition");
            DetectorResponseRbfComponent* target = nullptr;
            if (component_name == "photoelectric_escape")
                target = &detector.detector_photoelectric_escape;
            else if (component_name == "single_compton")
                target = &detector.detector_single_compton;
            else if (component_name == "multiple_compton")
                target = &detector.detector_multiple_compton;
            else if (component_name == "pair_no_compton")
                target = &detector.detector_pair_no_compton;
            else if (component_name == "pair_single_compton")
                target = &detector.detector_pair_single_compton;
            else if (component_name == "pair_multiple_compton")
                target = &detector.detector_pair_multiple_compton;
            else
                throw std::runtime_error(
                    "Unknown embedded detector block component "
                    + component_name);
            target->block_sizes.fill(0);
            std::size_t index = 0;
            std::uint32_t size = 0;
            while ((in >> size) && size != 0) {
                if (index >= detail::kDetectorMaxBlocks || size < 4
                    || size > detail::kDetectorMaxBlockSize) {
                    throw std::runtime_error(
                        "Invalid embedded detector block partition for "
                        + component_name);
                }
                target->block_sizes[index++] = size;
            }
        } else if (key == "detector_pair_component_basis_sizes") {
            auto& detector = payload.detector;
            in >> detector.pair_fraction.basis_size
               >> detector.pair_full_energy.basis_size
               >> detector.detector_pair_no_compton.basis_size
               >> detector.detector_pair_single_compton.basis_size
               >> detector.detector_pair_multiple_compton.basis_size;
        } else if (key == "detector_pair_fraction_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.pair_fraction);
        } else if (key == "detector_pair_fep_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.pair_full_energy);
        } else if (key == "detector_zero_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.zero_deposit);
        } else if (key == "detector_fep_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.full_energy);
        } else if (key == "detector_single_escape_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.single_escape);
        } else if (key == "detector_double_escape_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.double_escape);
        } else if (key == "detector_annihilation_511_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.annihilation_511);
        } else if (key == "detector_k_escape_lines") {
            if (!(in >> payload.detector.k_escape_alpha_energy_keV
                     >> payload.detector.k_escape_beta_energy_keV
                     >> payload.detector.k_escape_edge_keV
                     >> payload.detector.k_escape_support_keV[0]
                     >> payload.detector.k_escape_support_keV[1]
                     >> payload.detector.k_escape_support_keV[2]
                     >> payload.detector.k_escape_support_keV[3])) {
                throw std::runtime_error(
                    "Truncated embedded detector K-escape line energies");
            }
        } else if (key == "detector_k_alpha_escape_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.k_alpha_escape);
        } else if (key == "detector_k_beta_escape_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.k_beta_escape);
        } else if (key == "detector_k_alpha_absorbed_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.k_alpha_absorbed);
        } else if (key == "detector_k_beta_absorbed_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.k_beta_absorbed);
        } else if (key == "detector_photoelectric_escape_rbf") {
            read_rbf_component_body(
                in, key,
                payload.detector.detector_photoelectric_escape);
        } else if (key == "detector_single_compton_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.detector_single_compton);
        } else if (key == "detector_multiple_compton_rbf") {
            read_rbf_component_body(
                in, key, payload.detector.detector_multiple_compton);
        } else if (key == "detector_pair_no_compton_rbf") {
            read_rbf_component_body(
                in, key,
                payload.detector.detector_pair_no_compton);
        } else if (key == "detector_pair_single_compton_rbf") {
            read_rbf_component_body(
                in, key,
                payload.detector.detector_pair_single_compton);
        } else if (key == "detector_pair_multiple_compton_rbf") {
            read_rbf_component_body(
                in, key,
                payload.detector.detector_pair_multiple_compton);
        } else if (key == "coefficient_standard_errors") {
            std::size_t count = 0;
            in >> count;
            payload.coefficient_standard_errors.resize(count);
            for (double& error : payload.coefficient_standard_errors)
                in >> error;
        } else {
            throw std::runtime_error(
                "Unknown FullSpecApprox payload field: " + key);
        }
        if (!in)
            throw std::runtime_error("Malformed FullSpecApprox payload");
    }
    if (!saw_end)
        throw std::runtime_error("Truncated FullSpecApprox payload");
    finalize_source_payload(payload.source);
    if (payload.detector.detected_rate_surface.anchor.enabled
        && !payload.detector.geometry)
        throw std::runtime_error(
            "An embedded detector rate-surface anchor needs the detector "
            "geometry");
    apply_payload_geometry(payload.detector);
    if (payload.detector.detected_rate_surface.anchor.enabled)
        validate_rate_surface_anchor(payload.detector.detected_rate_surface.anchor);
    return FastSpectrumModel(std::move(payload));
}

} // namespace ceelo::full_spectrum
