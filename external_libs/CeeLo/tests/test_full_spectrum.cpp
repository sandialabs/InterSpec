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

#define BOOST_TEST_MODULE FullSpectrumApproxTests
#include <boost/test/unit_test.hpp>

#include "cross_sections/CrossSectionData.h"
#include "full_spectrum/FullSpectrumModel.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

namespace fs = ceelo::full_spectrum;

namespace {

std::vector<double> uniform_edges(double lo, double hi, std::size_t bins) {
    std::vector<double> edges(bins + 1);
    for (std::size_t i = 0; i <= bins; ++i)
        edges[i] = lo + (hi - lo) * static_cast<double>(i) / bins;
    return edges;
}

double sum(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0);
}

template <typename Scalar>
Scalar exact_component_yield(
    const fs::SpectralMeasure<Scalar>& measure,
    fs::SpectralComponent component)
{
    Scalar answer = Scalar(0);
    for (const auto& atom : measure.atoms)
        if (atom.component == component)
            answer += atom.yield;
    for (const auto& bin : measure.continuous)
        if (bin.component == component)
            answer += bin.integral;
    return answer;
}

} // namespace

namespace {

fs::ChargedParticleSplineField make_test_spline_field(
    std::uint32_t shape_rank)
{
    fs::ChargedParticleSplineField field;
    field.domain = {0.1, 1.05, 0.0, 1.0, 0.02, 0.9};
    field.yield_grid = {6, 5, 8};
    field.yield_coefficients.resize(6 * 5 * 8);
    for (std::size_t index = 0;
         index < field.yield_coefficients.size(); ++index) {
        field.yield_coefficients[index] =
            -0.4 + 0.013 * static_cast<double>((7 * index + 3) % 17);
    }
    field.shape_rank = shape_rank;
    field.shape_basis.resize(shape_rank);
    for (std::uint32_t rank = 0; rank < shape_rank; ++rank) {
        double total = 0.0;
        for (std::size_t k = 0; k < fs::kMultipleScatterBasisSize; ++k) {
            field.shape_basis[rank][k] =
                1.0 + static_cast<double>((5 * k + 3 * rank) % 11);
            total += field.shape_basis[rank][k];
        }
        for (auto& value : field.shape_basis[rank])
            value /= total;
    }
    if (shape_rank >= 2) {
        field.weight_grid = {4, 4, 5};
        field.weight_log_odds.resize(shape_rank - 1);
        for (std::size_t rank = 0;
             rank + 1 < shape_rank; ++rank) {
            field.weight_log_odds[rank].resize(4 * 4 * 5);
            for (std::size_t index = 0;
                 index < field.weight_log_odds[rank].size(); ++index) {
                field.weight_log_odds[rank][index] =
                    (static_cast<double>((5 * index + rank) % 13) / 13.0)
                    - 0.5;
            }
        }
    }
    return field;
}

} // namespace

BOOST_AUTO_TEST_CASE(source_zero_density_and_direct_invariant) {
    fs::SourceGeomFluxModel model;
    const auto edges = uniform_edges(0.0, 1500.0, 256);
    const auto zero = model.evaluate(1332.5, 0.0, 26.0, edges);
    BOOST_TEST(
        zero.component_yield(fs::SpectralComponent::Uncollided) == 1.0,
        boost::test_tools::tolerance(1e-14));
    BOOST_TEST(
        zero.component_yield(fs::SpectralComponent::ElasticLine) == 0.0,
        boost::test_tools::tolerance(1e-14));
    BOOST_TEST(
        zero.component_yield(fs::SpectralComponent::SingleCompton) == 0.0,
        boost::test_tools::tolerance(1e-14));
    BOOST_TEST(
        zero.component_yield(fs::SpectralComponent::MultipleCompton) == 0.0,
        boost::test_tools::tolerance(1e-14));
    BOOST_TEST(
        zero.component_yield(fs::SpectralComponent::PairAnnihilation) == 0.0,
        boost::test_tools::tolerance(1e-14));

    const double areal_density = 2.5;
    const auto shielded =
        model.evaluate(661.657, areal_density, 82.0, edges);
    const auto xs = fs::element_mass_cross_sections(82, 661.657);
    const double analytic = std::exp(
        -areal_density * (xs.pe + xs.cs + xs.rs + xs.pp));
    const double approximated =
        shielded.component_yield(fs::SpectralComponent::Uncollided);
    BOOST_TEST(approximated == analytic, boost::test_tools::tolerance(0.005));
    BOOST_TEST(
        shielded.component_yield(fs::SpectralComponent::ElasticLine) > 0.0);
    BOOST_TEST(
        shielded.component_yield(
            fs::SpectralComponent::SingleComptonRayleigh) > 0.0);
    BOOST_TEST(
        shielded.component_yield(
            fs::SpectralComponent::DoubleComptonRayleigh) > 0.0);
    const double expected_single =
        exact_component_yield(
            shielded, fs::SpectralComponent::SingleCompton)
        + exact_component_yield(
            shielded, fs::SpectralComponent::ComptonRayleigh)
        + exact_component_yield(
            shielded, fs::SpectralComponent::SingleComptonRayleigh);
    BOOST_TEST(
        shielded.component_yield(fs::SpectralComponent::SingleCompton)
            == expected_single,
        boost::test_tools::tolerance(1e-13));
    const double expected_multiple =
        exact_component_yield(
            shielded, fs::SpectralComponent::MultipleCompton)
        + exact_component_yield(
            shielded, fs::SpectralComponent::DoubleCompton)
        + exact_component_yield(
            shielded, fs::SpectralComponent::HigherCompton)
        + exact_component_yield(
            shielded, fs::SpectralComponent::DoubleComptonRayleigh)
        + exact_component_yield(
            shielded, fs::SpectralComponent::HigherComptonRayleigh)
        + exact_component_yield(
            shielded, fs::SpectralComponent::MultipleComptonRayleigh);
    BOOST_TEST(
        shielded.component_yield(fs::SpectralComponent::MultipleCompton)
            == expected_multiple,
        boost::test_tools::tolerance(1e-13));

    fs::SourceModelPayload no_rayleigh_payload;
    no_rayleigh_payload.enable_rayleigh = false;
    const fs::SourceGeomFluxModel no_rayleigh_model(no_rayleigh_payload);
    const auto without_rayleigh =
        no_rayleigh_model.evaluate(661.657, areal_density, 82.0, edges);
    const double analytic_without_rayleigh = std::exp(
        -areal_density * (xs.pe + xs.cs + xs.pp));
    BOOST_TEST(
        without_rayleigh.component_yield(fs::SpectralComponent::Uncollided)
            == analytic_without_rayleigh,
        boost::test_tools::tolerance(0.005));
    BOOST_TEST(
        without_rayleigh.component_yield(
            fs::SpectralComponent::ElasticLine) == 0.0,
        boost::test_tools::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(source_outputs_are_nonnegative_and_supported) {
    fs::SourceModelPayload payload;
    payload.enable_source_electron_brems = true;
    fs::SourceGeomFluxModel model(payload);
    for (double energy : {10.0, 511.0, 1021.0, 1023.0, 12000.0}) {
        const auto edges = uniform_edges(0.0, energy, 80);
        for (double Z : {1.0, 13.5, 47.0, 82.0, 92.0}) {
            const auto measure = model.evaluate(energy, 3.0, Z, edges);
            for (const auto& atom : measure.atoms) {
                BOOST_TEST(std::isfinite(atom.yield));
                BOOST_TEST(atom.yield >= -1e-14);
                BOOST_TEST(atom.energy_keV >= 0.0);
                BOOST_TEST(atom.energy_keV <= energy);
            }
            for (const auto& bin : measure.continuous) {
                BOOST_TEST(std::isfinite(bin.integral));
                BOOST_TEST(bin.integral >= -1e-13);
                BOOST_TEST(bin.lower_keV >= 0.0);
                BOOST_TEST(bin.upper_keV <= energy);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(charged_particle_source_is_positive_reversible_and_zero_in_vacuum) {
    fs::ChargedParticleSourcePayload payload;
    payload.electron_brems_spline = make_test_spline_field(3);
    payload.positron_brems_spline = make_test_spline_field(3);
    payload.positron_containment_spline = make_test_spline_field(0);
    payload.positron_511_atom_spline = make_test_spline_field(0);
    payload.positron_annihilation_continuum_spline =
        make_test_spline_field(2);
    payload.positron_annihilation_in_flight_direct_spline =
        make_test_spline_field(3);
    payload.positron_annihilation_in_flight_continuum_spline =
        make_test_spline_field(1);

    const fs::ChargedParticleSourceFluxModel original(payload);
    const auto decoded =
        fs::ChargedParticleSourceFluxModel::deserialize(
            original.serialize());
    const auto bins = uniform_edges(0.0, 3000.0, 600);
    for (bool positron : {false, true}) {
        const auto vacuum =
            decoded.evaluate(2290.0, 0.0, 35.0, positron, bins);
        BOOST_TEST(
            vacuum.integrated_yield() == 0.0,
            boost::test_tools::tolerance(1e-14));

        const auto measure =
            decoded.evaluate(2290.0, 1.0, 35.0, positron, bins);
        BOOST_TEST(measure.integrated_yield() > 0.0);
        for (const auto& atom : measure.atoms) {
            BOOST_TEST(std::isfinite(atom.yield));
            BOOST_TEST(atom.yield >= 0.0);
            BOOST_TEST(atom.energy_keV >= 0.0);
            BOOST_TEST(atom.energy_keV <= 2800.99895);
        }
        for (const auto& bin : measure.continuous) {
            BOOST_TEST(std::isfinite(bin.integral));
            BOOST_TEST(bin.integral >= 0.0);
            BOOST_TEST(bin.lower_keV >= 0.0);
            BOOST_TEST(bin.lower_keV < 2800.99895);
            BOOST_TEST(bin.upper_keV <= bins.back());
        }
        if (positron) {
            BOOST_TEST(
                measure.component_yield(
                    fs::SpectralComponent::PrimaryPositronBrems) > 0.0);
            BOOST_TEST(
                measure.component_yield(
                    fs::SpectralComponent::
                        PrimaryPositronAnnihilation) > 0.0);
            BOOST_TEST(
                measure.component_yield(
                    fs::SpectralComponent::
                        PrimaryPositronAnnihilationContinuum) > 0.0);
            BOOST_TEST(
                measure.component_yield(
                    fs::SpectralComponent::
                        PrimaryPositronAnnihilationInFlight) > 0.0);
            BOOST_TEST(
                measure.component_yield(
                    fs::SpectralComponent::
                        PrimaryPositronAnnihilationInFlightContinuum)
                    > 0.0);
        } else {
            BOOST_TEST(
                measure.component_yield(
                    fs::SpectralComponent::PrimaryElectronBrems) > 0.0);
            BOOST_TEST(measure.atoms.empty());
        }
    }

    // Sub-threshold positrons still annihilate; only their bremsstrahlung is
    // below the model's photon cutoff. Electrons below that cutoff have no
    // supported photon source term.
    const auto low_energy_electron =
        decoded.evaluate(5.0, 1.0, 35.0, false, bins);
    const auto low_energy_positron =
        decoded.evaluate(5.0, 1.0, 35.0, true, bins);
    BOOST_TEST(
        low_energy_electron.integrated_yield() == 0.0,
        boost::test_tools::tolerance(1e-14));
    BOOST_TEST(
        low_energy_positron.component_yield(
            fs::SpectralComponent::PrimaryPositronAnnihilation) > 0.0);

    fs::FullSpectrumPayload full_payload;
    full_payload.source_quadrature_bins = 8;
    const fs::FastSpectrumModel composed(full_payload, decoded);
    const fs::ShieldParameters<double> shield{1.0, 35.0};
    const fs::ResolutionParameters<double> resolution{
        0.4, 0.001, 1.0e-6};
    const std::vector<fs::ChargedParticleLine> charged_lines{
        {5.0, 1.0, true}};
    const auto charged_spectrum = composed.predict(
        std::vector<fs::EmissionLine>{}, charged_lines, 1000.0,
        shield, resolution, bins);
    BOOST_TEST(sum(charged_spectrum) > 0.0);
}


BOOST_AUTO_TEST_CASE(charged_particle_payload_round_trip_is_exact) {
    // Some components fitted, the rest absent: both record kinds round-trip.
    fs::ChargedParticleSourcePayload payload;
    payload.electron_brems_spline = make_test_spline_field(3);
    payload.positron_containment_spline = make_test_spline_field(0);
    payload.positron_annihilation_continuum_spline =
        make_test_spline_field(2);

    const fs::ChargedParticleSourceFluxModel original(payload);
    const std::string text = original.serialize();
    BOOST_TEST(
        text.rfind("CEELO_CHARGED_PARTICLE_SOURCE 1\n", 0) == 0u);
    BOOST_TEST(text.find("electron_brems spline") != std::string::npos);
    BOOST_TEST(text.find("positron_brems none") != std::string::npos);
    BOOST_TEST(
        text.find("positron_containment spline") != std::string::npos);
    const auto decoded =
        fs::ChargedParticleSourceFluxModel::deserialize(text);
    // Serialization must round-trip exactly.
    BOOST_TEST(decoded.serialize() == text);
    const auto bins = uniform_edges(0.0, 3000.0, 600);
    for (bool positron : {false, true}) {
        const auto a = original.evaluate(2290.0, 1.0, 35.0, positron, bins);
        const auto b = decoded.evaluate(2290.0, 1.0, 35.0, positron, bins);
        BOOST_TEST(a.atoms.size() == b.atoms.size());
        BOOST_TEST(a.continuous.size() == b.continuous.size());
        for (std::size_t i = 0; i < a.atoms.size(); ++i)
            BOOST_TEST(
                a.atoms[i].yield == b.atoms[i].yield,
                boost::test_tools::tolerance(1e-14));
        for (std::size_t i = 0; i < a.continuous.size(); ++i)
            BOOST_TEST(
                a.continuous[i].integral == b.continuous[i].integral,
                boost::test_tools::tolerance(1e-14));
    }
}

BOOST_AUTO_TEST_CASE(
    charged_particle_spline_positivity_and_frozen_extrapolation) {
    // Constant-coefficient yield field: the tensor cubic B-spline basis is a
    // partition of unity, so conditional_yield == exp(c) everywhere and the
    // component yield is exactly base * exp(c). A non-constant field frozen
    // beyond its trained f2 domain must reproduce the same base ratio.
    const auto bins = uniform_edges(0.0, 3000.0, 600);
    const auto electron_yield = [&](
        const fs::ChargedParticleSourceFluxModel& model, double rho) {
        const auto measure = model.evaluate(2290.0, rho, 35.0, false, bins);
        return measure.component_yield(
            fs::SpectralComponent::PrimaryElectronBrems);
    };

    fs::ChargedParticleSourcePayload constant_payload;
    constant_payload.electron_brems_spline = make_test_spline_field(0);
    std::fill(
        constant_payload.electron_brems_spline.yield_coefficients.begin(),
        constant_payload.electron_brems_spline.yield_coefficients.end(),
        0.37);
    const fs::ChargedParticleSourceFluxModel constant_model(
        constant_payload);

    fs::ChargedParticleSourcePayload frozen_payload;
    frozen_payload.electron_brems_spline = make_test_spline_field(3);
    // Trained f2 domain ends well below the queried thicknesses, so both
    // queries evaluate the frozen (constant) boundary field.
    frozen_payload.electron_brems_spline.domain = {
        0.1, 1.05, 0.0, 1.0, 0.02, 0.15};
    const fs::ChargedParticleSourceFluxModel frozen_model(frozen_payload);

    // rho = 3.4 and 6.8 g/cm2 give range ratios of roughly 2 and 4 at
    // 2290 keV, Z = 35: f2 well beyond 0.15 with distinct mechanistic bases.
    const double constant_ratio =
        electron_yield(constant_model, 6.8)
        / electron_yield(constant_model, 3.4);
    const double frozen_ratio =
        electron_yield(frozen_model, 6.8)
        / electron_yield(frozen_model, 3.4);
    BOOST_TEST(
        frozen_ratio == constant_ratio,
        boost::test_tools::tolerance(1e-12));

    // Positivity and finiteness at extreme queries, including far outside
    // the validity box.
    for (double rho : {1.0e-4, 1.0, 60.0, 1000.0}) {
        for (double z : {1.0, 35.25, 92.0}) {
            for (bool positron : {false, true}) {
                const auto measure = frozen_model.evaluate(
                    2290.0, rho, z, positron, bins);
                for (const auto& atom : measure.atoms) {
                    BOOST_TEST(std::isfinite(atom.yield));
                    BOOST_TEST(atom.yield >= 0.0);
                }
                for (const auto& bin : measure.continuous) {
                    BOOST_TEST(std::isfinite(bin.integral));
                    BOOST_TEST(bin.integral >= 0.0);
                }
            }
        }
    }
    // Vacuum still yields exactly zero.
    const auto vacuum = frozen_model.evaluate(2290.0, 0.0, 35.0, false, bins);
    BOOST_TEST(
        vacuum.integrated_yield() == 0.0,
        boost::test_tools::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(photonuclear_correction_is_positive_and_conservative) {
    fs::SourceModelPayload payload;
    payload.enable_photonuclear = true;
    fs::PhotonuclearCrossSectionRow lower;
    lower.energy_keV = 1000.0;
    lower.mass_attenuation_cm2_g.fill(0.001);
    fs::PhotonuclearCrossSectionRow upper;
    upper.energy_keV = 2000.0;
    upper.mass_attenuation_cm2_g.fill(0.003);
    payload.photonuclear.cross_sections = {lower, upper};
    payload.photonuclear.escaped_yield_scale = 0.85;
    fs::SourceGeomFluxModel model(payload);
    const auto bins = uniform_edges(0.0, 1500.0, 256);
    const double areal_density = 20.0;
    const auto measure =
        model.evaluate(1500.0, areal_density, 26.0, bins);
    const fs::SmoothElementCrossSection smooth_cross_sections(
        payload.z_blend_width);
    const auto electromagnetic =
        smooth_cross_sections.evaluate(1500.0, 26.0);
    const double expected_direct = std::exp(
        -areal_density
        * (electromagnetic[0] + electromagnetic[1]
           + electromagnetic[2] + electromagnetic[3] + 0.002));
    BOOST_TEST(
        measure.component_yield(fs::SpectralComponent::Uncollided)
            == expected_direct,
        boost::test_tools::tolerance(1e-12));
    BOOST_TEST(
        measure.component_yield(
            fs::SpectralComponent::PhotonuclearSecondary)
            == 0.85 * (1.0 - std::exp(-areal_density * 0.002)),
        boost::test_tools::tolerance(1e-12));

    const auto vacuum = model.evaluate(1500.0, 0.0, 26.0, bins);
    BOOST_TEST(
        vacuum.component_yield(
            fs::SpectralComponent::PhotonuclearSecondary) == 0.0,
        boost::test_tools::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(single_compton_has_exact_kinematic_support) {
    fs::SourceModelPayload payload;
    payload.enable_compton_doppler = false;
    fs::SourceGeomFluxModel model(payload);
    const double energy_keV = 300.0;
    const double minimum_keV =
        energy_keV / (1.0 + 2.0 * energy_keV / 510.99895);
    const auto edges = uniform_edges(0.0, energy_keV, 600);
    const auto measure = model.evaluate(
        energy_keV, 7.874, 26.0, edges);
    double below_support = 0.0;
    double in_support = 0.0;
    for (const auto& bin : measure.continuous) {
        if (bin.component != fs::SpectralComponent::SingleCompton)
            continue;
        if (bin.upper_keV <= minimum_keV)
            below_support += bin.integral;
        else
            in_support += bin.integral;
    }
    BOOST_TEST(below_support == 0.0, boost::test_tools::tolerance(1e-14));
    BOOST_TEST(in_support > 0.0);
}

BOOST_AUTO_TEST_CASE(single_compton_doppler_populates_bound_electron_tail) {
    fs::SourceModelPayload payload;
    payload.enable_compton_doppler = true;
    fs::SourceGeomFluxModel model(payload);
    const double energy_keV = 100.0;
    const double free_minimum_keV =
        energy_keV / (1.0 + 2.0 * energy_keV / 510.99895);
    const auto edges = uniform_edges(0.0, energy_keV, 400);
    const auto measure = model.evaluate(
        energy_keV, 2.0, 82.0, edges);
    double below_free_support = 0.0;
    for (const auto& bin : measure.continuous) {
        if (bin.component == fs::SpectralComponent::SingleCompton
            && bin.upper_keV <= free_minimum_keV) {
            below_free_support += bin.integral;
        }
    }
    BOOST_TEST(below_free_support > 0.0);
}

BOOST_AUTO_TEST_CASE(double_compton_is_smooth_across_former_angle_nodes) {
    fs::SourceModelPayload payload;
    payload.enable_compton_doppler = false;
    fs::SourceGeomFluxModel model(payload);
    constexpr double energy_keV = 3000.0;
    constexpr double effective_Z = 35.0;
    const fs::SmoothElementCrossSection cross_sections(payload.z_blend_width);
    const auto mu = cross_sections.evaluate(energy_keV, effective_Z);
    const double areal_density = 1.0 / (mu[0] + mu[1] + mu[2] + mu[3]);

    const auto maximum_relative_adjacent_change =
        [&](double lower_keV, double upper_keV,
            fs::SpectralComponent component) {
            const auto edges =
                uniform_edges(lower_keV, upper_keV,
                              static_cast<std::size_t>(upper_keV - lower_keV));
            const auto measure =
                model.evaluate(energy_keV, areal_density, effective_Z, edges);
            std::vector<double> density(edges.size() - 1, 0.0);
            for (const auto& bin : measure.continuous) {
                if (bin.component != component)
                    continue;
                const auto index = static_cast<std::size_t>(
                    std::llround(bin.lower_keV - lower_keV));
                BOOST_REQUIRE(index < density.size());
                density[index] +=
                    bin.integral / (bin.upper_keV - bin.lower_keV);
            }
            double maximum = 0.0;
            for (std::size_t index = 1; index < density.size(); ++index) {
                const double scale =
                    std::max(density[index - 1], density[index]);
                BOOST_REQUIRE(scale > 0.0);
                maximum = std::max(
                    maximum,
                    std::abs(density[index] - density[index - 1]) / scale);
            }
            return maximum;
        };

    // A global eight-point first-angle rule produced 30-60% steps at
    // 792.49 and 1367.53 keV for this cell. Kinematically segmented panels
    // make the integral continuous across those former node energies.
    for (const auto component :
         {fs::SpectralComponent::DoubleCompton,
          fs::SpectralComponent::DoubleComptonRayleigh}) {
        BOOST_TEST(maximum_relative_adjacent_change(780.0, 805.0, component) <
                   0.02);
        BOOST_TEST(maximum_relative_adjacent_change(1355.0, 1380.0, component) <
                   0.02);
    }

    const auto coarse = model.evaluate(energy_keV, areal_density, effective_Z,
                                       uniform_edges(0.0, energy_keV, 256));
    const auto fine = model.evaluate(energy_keV, areal_density, effective_Z,
                                     uniform_edges(0.0, energy_keV, 1024));
    for (const auto component :
         {fs::SpectralComponent::DoubleCompton,
          fs::SpectralComponent::DoubleComptonRayleigh}) {
        const double coarse_yield = exact_component_yield(coarse, component);
        const double fine_yield = exact_component_yield(fine, component);
        BOOST_TEST(coarse_yield == fine_yield,
                   boost::test_tools::tolerance(0.002));
    }
}

BOOST_AUTO_TEST_CASE(detector_measure_normalizes_per_incident_photon) {
    for (const auto payload : {
             fs::nai3x3_payload(), fs::coaxial_hpge_payload()}) {
        fs::DetectorPulseResponseModel model(payload);
        for (double energy : {10.0, 511.0, 1023.0, 2000.0, 12000.0}) {
            const auto edges = uniform_edges(0.0, energy, 256);
            for (double cosine : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
                const auto response = model.evaluate(energy, cosine, edges);
                BOOST_TEST(
                    response.integrated_yield() == 1.0,
                    boost::test_tools::tolerance(2e-12));
                for (const auto& atom : response.atoms) {
                    BOOST_TEST(std::isfinite(atom.yield));
                    BOOST_TEST(atom.yield >= 0.0);
                }
                for (const auto& bin : response.continuous) {
                    BOOST_TEST(std::isfinite(bin.integral));
                    BOOST_TEST(bin.integral >= 0.0);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(detector_surrogate_round_trip_and_pair_threshold) {
    auto payload = fs::coaxial_hpge_payload();
    payload.use_response_surrogate = true;
    const auto populate =
        [](fs::DetectorResponseRbfComponent& component, double yield) {
            fs::DetectorResponseRbfNode node;
            node.feature = {0.5, 0.0, 1.0};
            node.log_conditional_yield = std::log(yield);
            node.shape_coefficients.fill(
                1.0 / fs::kDetectorResponseBasisSize);
            component.nodes.push_back(node);
        };
    // zero_deposit stores detected/non-detected odds; all other fields are
    // conditional fractions which the evaluator normalizes.
    populate(payload.zero_deposit, 4.0);
    populate(payload.full_energy, 0.25);
    populate(payload.pair_fraction, 0.30);
    populate(payload.pair_full_energy, 0.10);
    populate(payload.single_escape, 0.04);
    populate(payload.double_escape, 0.02);
    populate(payload.annihilation_511, 0.015);
    payload.k_escape_alpha_energy_keV = 9.8331279752;
    payload.k_escape_beta_energy_keV = 10.9426924714;
    payload.k_escape_edge_keV = 11.103;
    populate(payload.k_alpha_escape, 0.012);
    populate(payload.k_beta_escape, 0.003);
    populate(payload.k_alpha_absorbed, 0.020);
    populate(payload.k_beta_absorbed, 0.004);
    populate(payload.detector_photoelectric_escape, 0.035);
    populate(payload.detector_single_compton, 0.25);
    populate(payload.detector_multiple_compton, 0.38);
    populate(payload.detector_pair_no_compton, 0.02);
    populate(payload.detector_pair_single_compton, 0.025);
    populate(payload.detector_pair_multiple_compton, 0.015);
    payload.zero_deposit.basis_size = 1;
    payload.full_energy.basis_size = 1;
    payload.pair_fraction.basis_size = 1;
    payload.pair_full_energy.basis_size = 1;
    payload.single_escape.basis_size = 1;
    payload.double_escape.basis_size = 1;
    payload.annihilation_511.basis_size = 1;
    payload.k_alpha_escape.basis_size = 1;
    payload.k_beta_escape.basis_size = 1;
    payload.k_alpha_absorbed.basis_size = 1;
    payload.k_beta_absorbed.basis_size = 1;
    payload.detector_photoelectric_escape.basis_size = 16;
    payload.detector_single_compton.basis_size = 48;
    payload.detector_single_compton.block_sizes = {32, 16, 0, 0};
    payload.detector_multiple_compton.basis_size = 64;
    payload.detector_multiple_compton.block_sizes = {43, 21, 0, 0};
    payload.detector_pair_no_compton.basis_size = 32;
    payload.detector_pair_no_compton.block_sizes = {8, 8, 8, 8};
    payload.detector_pair_single_compton.basis_size = 48;
    payload.detector_pair_single_compton.block_sizes = {12, 12, 12, 12};
    payload.detector_pair_multiple_compton.basis_size = 64;
    payload.detector_pair_multiple_compton.block_sizes = {16, 16, 16, 16};
    payload.detected_rate_surface.enabled = true;
    for (std::size_t index = 0;
         index < payload.detected_rate_surface
                       .log_odds_coefficients.size(); ++index) {
        payload.detected_rate_surface.log_odds_coefficients[index] =
            -0.4 + 0.002 * index;
    }

    const fs::DetectorPulseResponseModel original(payload);
    const auto restored = fs::DetectorPulseResponseModel::deserialize(
        original.serialize());
    const auto edges = uniform_edges(0.0, 1332.492, 512);
    const auto response = restored.evaluate(1332.492, 0.35, edges);
    BOOST_TEST(
        response.integrated_yield() == 1.0,
        boost::test_tools::tolerance(2e-12));
    BOOST_TEST(
        exact_component_yield(
            response, fs::SpectralComponent::DetectorDoubleEscape) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorAnnihilation511) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorPhotoelectricEscape) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorKAlphaEscape) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorKBetaEscape) > 0.0);
    // The K-escape atoms sit at Ein minus the payload's line energies.
    {
        bool found_alpha = false;
        for (const auto& atom : response.atoms) {
            if (atom.component
                == fs::SpectralComponent::DetectorKAlphaEscape) {
                found_alpha = true;
                BOOST_TEST(
                    atom.energy_keV == 1332.492 - 9.8331279752,
                    boost::test_tools::tolerance(1e-12));
            }
        }
        BOOST_TEST(found_alpha);
    }
    // Absorbed dead-layer fluorescence deposits the line energy itself, so
    // its atoms sit at a fixed energy rather than tracking Ein.
    {
        bool found_absorbed = false;
        for (const auto& atom : response.atoms) {
            if (atom.component
                == fs::SpectralComponent::DetectorKAlphaAbsorbed) {
                found_absorbed = true;
                BOOST_TEST(
                    atom.energy_keV == 9.8331279752,
                    boost::test_tools::tolerance(1e-12));
            }
        }
        BOOST_TEST(found_absorbed);
    }
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorKAlphaAbsorbed) > 0.0);

    // Below the crystal K edge the atoms are gated off entirely.
    const auto below_k_edge = restored.evaluate(11.0, 0.35, edges);
    BOOST_TEST(
        exact_component_yield(
            below_k_edge,
            fs::SpectralComponent::DetectorKAlphaEscape) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            below_k_edge,
            fs::SpectralComponent::DetectorKBetaEscape) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            below_k_edge,
            fs::SpectralComponent::DetectorKAlphaAbsorbed) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorPairFullEnergy) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorPairNoCompton) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorPairSingleCompton) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            response,
            fs::SpectralComponent::DetectorPairMultipleCompton) > 0.0);

    const auto below_pair = restored.evaluate(1000.0, 0.35, edges);
    BOOST_TEST(
        exact_component_yield(
            below_pair, fs::SpectralComponent::DetectorSingleEscape) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            below_pair, fs::SpectralComponent::DetectorDoubleEscape) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            below_pair,
            fs::SpectralComponent::DetectorPairFullEnergy) == 0.0);
    const auto at_pair_threshold =
        restored.evaluate(1021.9979, 0.35, edges);
    BOOST_TEST(
        exact_component_yield(
            at_pair_threshold,
            fs::SpectralComponent::DetectorPairFullEnergy) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            at_pair_threshold,
            fs::SpectralComponent::DetectorPairNoCompton) == 0.0);

    fs::FullSpectrumPayload full_payload;
    full_payload.detector = payload;
    const auto embedded = fs::FastSpectrumModel::deserialize(
        fs::FastSpectrumModel(full_payload).serialize());
    BOOST_TEST(embedded.payload().detector.use_response_surrogate);
    BOOST_TEST(
        embedded.payload().detector.detected_rate_surface.enabled);
    BOOST_TEST(
        embedded.payload().detector.zero_deposit.nodes.size() == 1u);
    BOOST_TEST(
        embedded.payload().detector.zero_deposit.basis_size == 1u);
    BOOST_TEST(
        embedded.payload()
            .detector.detector_photoelectric_escape.basis_size == 16u);
    BOOST_TEST(
        embedded.payload().detector.detector_single_compton.basis_size
        == 48u);
    BOOST_TEST(
        embedded.payload().detector.detector_multiple_compton.basis_size
        == 64u);
    BOOST_TEST(
        embedded.payload().detector.pair_fraction.nodes.size() == 1u);
    BOOST_TEST(
        embedded.payload()
            .detector.detector_pair_no_compton.basis_size == 32u);
}

namespace {

// The detector campaign's built-in hpge_coax, as `--print-geometry` writes it.
const char* const kHpgeCoaxGeometry =
    "CEELO_GEOMETRY 1\n"
    "shape cylinder\n"
    "dimensions 4 8\n"
    "bullet_radius 0\n"
    "crystal_material 0\n"
    "bore 0.5 5.5 0\n"
    "dead_layer 0.05 0.08 0.08\n"
    "layer 1 0.1 0.1 -0.1 8.1\n"
    "symmetry axial\n"
    "reference_point crystal_face\n"
    "material 0 5.323 1 32 1 HPGe\n"
    "material 1 2.7 1 13 1 Al\n"
    "end\n";

// A populated hierarchical payload whose continuum components use the block
// partition, for the structural invariants below.
fs::DetectorModelPayload blocked_hpge_payload() {
    auto payload = fs::coaxial_hpge_payload();
    payload.detector_name = "hpge_coax";
    fs::set_detector_geometry(
        payload, ceelo::GeometryDescriptor::from_text(kHpgeCoaxGeometry));
    payload.use_response_surrogate = true;
    const auto populate =
        [](fs::DetectorResponseRbfComponent& component, double yield,
           std::uint32_t basis_size,
           std::array<std::uint32_t, 4> blocks = {}) {
            component.basis_size = basis_size;
            component.block_sizes = blocks;
            fs::DetectorResponseRbfNode node;
            node.feature = {0.5, 0.0, 1.0};
            node.log_conditional_yield = std::log(yield);
            // A varying, strictly positive shape exercises every basis
            // function rather than a flat vector that hides index errors.
            for (std::size_t index = 0; index < basis_size; ++index) {
                node.shape_coefficients[index] =
                    (1.0 + 0.5 * std::sin(0.7 * index)) / basis_size;
            }
            component.nodes.push_back(node);
        };
    populate(payload.zero_deposit, 4.0, 1);
    populate(payload.full_energy, 0.25, 1);
    populate(payload.pair_fraction, 0.30, 1);
    populate(payload.pair_full_energy, 0.10, 1);
    populate(payload.single_escape, 0.04, 1);
    populate(payload.double_escape, 0.02, 1);
    populate(payload.annihilation_511, 0.015, 1);
    populate(payload.detector_photoelectric_escape, 0.035, 16);
    populate(payload.detector_single_compton, 0.25, 48, {32, 16, 0, 0});
    populate(payload.detector_multiple_compton, 0.38, 96, {64, 32, 0, 0});
    populate(payload.detector_pair_no_compton, 0.02, 96, {24, 24, 24, 24});
    populate(
        payload.detector_pair_single_compton, 0.025, 64, {16, 16, 16, 16});
    populate(
        payload.detector_pair_multiple_compton, 0.015, 96, {24, 24, 24, 24});
    payload.detected_rate_surface.enabled = true;
    for (std::size_t index = 0;
         index < payload.detected_rate_surface
                       .log_odds_coefficients.size(); ++index) {
        payload.detected_rate_surface.log_odds_coefficients[index] =
            -0.4 + 0.002 * index;
    }
    return payload;
}

} // namespace

BOOST_AUTO_TEST_CASE(detector_rate_anchor_round_trip_and_identity) {
    auto payload = blocked_hpge_payload();
    // Disabled anchor must leave the detected rate bit-identical.
    const fs::DetectorPulseResponseModel plain(payload);
    payload.detected_rate_surface.anchor.enabled = true;
    const fs::DetectorPulseResponseModel anchored(payload);
    const auto text = anchored.serialize();
    BOOST_TEST(text.find("rate_surface_anchor_enabled 1")
               != std::string::npos);
    const auto restored =
        fs::DetectorPulseResponseModel::deserialize(text);
    BOOST_TEST(restored.serialize() == text);
    BOOST_TEST(restored.payload().detected_rate_surface.anchor.enabled);
    // The anchor's geometry is the payload's, rebuilt from the embedded text.
    BOOST_TEST(restored.payload().detector_name == "hpge_coax");
    BOOST_REQUIRE(restored.payload().geometry.has_value());
    BOOST_TEST(restored.payload().geometry->to_text()
               == payload.geometry->to_text());
    BOOST_TEST(restored.payload().passive_landmark_keV
               == payload.passive_landmark_keV);
    const auto& restored_anchor = restored.payload().detected_rate_surface.anchor;
    const auto& anchor = payload.detected_rate_surface.anchor;
    BOOST_TEST(restored_anchor.bore_radius_cm == 0.5);
    BOOST_TEST(restored_anchor.back_layers.size() == 1u);
    BOOST_TEST(restored_anchor.envelope_radius_cm == anchor.envelope_radius_cm);
    BOOST_TEST(restored_anchor.envelope_length_cm == anchor.envelope_length_cm);
    // The embedded (FULLSPEC) codec carries the same detector.
    fs::FullSpectrumPayload full_payload;
    full_payload.detector = payload;
    const auto embedded = fs::FastSpectrumModel::deserialize(
        fs::FastSpectrumModel(full_payload).serialize());
    BOOST_REQUIRE(embedded.payload().detector.geometry.has_value());
    BOOST_TEST(embedded.payload().detector.geometry->to_text()
               == payload.geometry->to_text());
    BOOST_TEST(embedded.payload().detector.detector_name == "hpge_coax");
    BOOST_TEST(embedded.payload().detector.passive_landmark_keV
               == payload.passive_landmark_keV);
    BOOST_TEST(
        embedded.payload().detector.detected_rate_surface.anchor.envelope_radius_cm
        == anchor.envelope_radius_cm);

    const auto edges = uniform_edges(0.0, 60.0, 128);
    const auto without = plain.evaluate(60.0, 0.85, edges);
    const auto with = restored.evaluate(60.0, 0.85, edges);
    // Both are normalized measures; the anchor moves the detected/zero split.
    BOOST_TEST(
        without.integrated_yield() == 1.0,
        boost::test_tools::tolerance(2e-12));
    BOOST_TEST(
        with.integrated_yield() == 1.0,
        boost::test_tools::tolerance(2e-12));
    const auto zero_yield = [](const auto& measure) {
        return exact_component_yield(
            measure, fs::SpectralComponent::DetectorZeroDeposit);
    };
    BOOST_TEST(zero_yield(with) != zero_yield(without));
    // The blend weight 1/(1+(E/150)^2) suppresses the anchor away from the
    // turn-on: by 6 MeV it moves the detected rate by only ~3e-4 relative
    // (and the fit absorbs even that, since it sees the same anchor).
    const auto high_edges = uniform_edges(0.0, 6000.0, 128);
    const auto high_plain = plain.evaluate(6000.0, 0.85, high_edges);
    const auto high_anchored = restored.evaluate(6000.0, 0.85, high_edges);
    BOOST_TEST(
        zero_yield(high_anchored) == zero_yield(high_plain),
        boost::test_tools::tolerance(1e-3));
}

BOOST_AUTO_TEST_CASE(detector_passive_landmark_follows_the_front_stack) {
    // hpge_coax's front entry is the reference: kPassiveLandmarkOpticalDepth
    // is its depth at 55 keV, the landmark the feature map was selected with.
    const auto payload = blocked_hpge_payload();
    BOOST_TEST(payload.passive_landmark_keV == 55.0,
               boost::test_tools::tolerance(1e-6));
    // A thicker window moves the turn-on, and the landmark with it.
    auto thicker = ceelo::GeometryDescriptor::from_text(kHpgeCoaxGeometry);
    thicker.layers.front().front_thickness_cm = 0.3;
    auto moved = fs::coaxial_hpge_payload();
    fs::set_detector_geometry(moved, thicker);
    BOOST_TEST(moved.passive_landmark_keV > 58.0);
    // The anchor normalizes per photon on the envelope, which the endcap
    // stack sets: a wider can is a larger target with the same crystal.
    BOOST_TEST(moved.detected_rate_surface.anchor.envelope_length_cm
               > payload.detected_rate_surface.anchor.envelope_length_cm);
    BOOST_TEST(moved.radius_cm == 4.0);
    BOOST_TEST(moved.density_g_cm3 == 5.323);
}

BOOST_AUTO_TEST_CASE(detector_rate_anchor_counts_the_bore_at_back_incidence) {
    // Backward photons that enter through the bore opening skip the back
    // stack; an anchor without the bore saw only the back face, and put the
    // turn-on of backward low-energy incidence far too low.
    auto payload = blocked_hpge_payload();
    auto& anchor = payload.detected_rate_surface.anchor;
    anchor.enabled = true;
    // The blended log odds are monotone in the entry probability.
    const auto log_odds = [&](double energy_keV, double cosine) {
        return fs::rate_surface_anchor_log_odds(anchor, energy_keV, cosine);
    };
    const double with_bore = log_odds(30.0, -1.0);
    anchor.bore_radius_cm = 0.0;
    const double without_bore = log_odds(30.0, -1.0);
    BOOST_TEST(with_bore > without_bore + 1.0);
    // Front incidence never sees the bore.
    anchor.bore_radius_cm = 0.5;
    const double front_with = log_odds(30.0, 1.0);
    anchor.bore_radius_cm = 0.0;
    BOOST_TEST(log_odds(30.0, 1.0) == front_with);
}

BOOST_AUTO_TEST_CASE(detector_rate_anchor_tolerates_cosine_overshoot) {
    // A caller whose cosine came from a dot product can hand us 1 + 1e-12.
    // Before the shrink, 1 - cos^2 went negative inside the anchor and NaN'd
    // every atom and every continuum bin of the measure.
    auto payload = blocked_hpge_payload();
    payload.detected_rate_surface.anchor.enabled = true;
    const fs::DetectorPulseResponseModel model(payload);
    const auto edges = uniform_edges(0.0, 40.0, 64);
    for (double cosine : {1.0, 1.0 + 1e-12, 1.0 + 1e-9, 1.0 + 1e-7,
                          -1.0, -1.0 - 1e-9}) {
        const auto measure = model.evaluate(40.0, cosine, edges);
        BOOST_TEST_CONTEXT("cos_theta " << cosine) {
            BOOST_TEST(std::isfinite(measure.integrated_yield()));
            BOOST_TEST(
                measure.integrated_yield() == 1.0,
                boost::test_tools::tolerance(2e-12));
            for (const auto& atom : measure.atoms)
                BOOST_TEST(std::isfinite(atom.yield));
            for (const auto& bin : measure.continuous)
                BOOST_TEST(std::isfinite(bin.integral));
        }
    }
}

BOOST_AUTO_TEST_CASE(detector_k_escape_support_limit_gates_atoms) {
    // The fitted K-escape yield fields are normalized kernel averages, so
    // above their node cloud they level off instead of decaying. The support
    // limit stops the model emitting a spurious escape satellite beside the
    // full-energy peak at energies the fit never saw.
    auto payload = blocked_hpge_payload();
    payload.k_escape_alpha_energy_keV = 9.8331279752;
    payload.k_escape_beta_energy_keV = 10.9426924714;
    payload.k_escape_edge_keV = 11.103;
    payload.k_escape_support_keV = {320.0, 175.0, 280.0, 150.0};
    const auto populate =
        [](fs::DetectorResponseRbfComponent& component, double yield) {
            fs::DetectorResponseRbfNode node;
            node.feature = {0.5, 0.0, 1.0};
            node.log_conditional_yield = std::log(yield);
            node.shape_coefficients.fill(
                1.0 / fs::kDetectorResponseBasisSize);
            component.basis_size = 1;
            component.nodes.push_back(node);
        };
    populate(payload.k_alpha_escape, 0.012);
    populate(payload.k_beta_escape, 0.003);
    populate(payload.k_alpha_absorbed, 0.020);
    populate(payload.k_beta_absorbed, 0.004);
    const fs::DetectorPulseResponseModel model(payload);
    const auto k_escape_yield = [](const auto& measure) {
        return exact_component_yield(
                   measure, fs::SpectralComponent::DetectorKAlphaEscape)
            + exact_component_yield(
                   measure, fs::SpectralComponent::DetectorKBetaEscape);
    };
    const auto inside =
        model.evaluate(120.0, 0.85, uniform_edges(0.0, 120.0, 128));
    BOOST_TEST(k_escape_yield(inside) > 0.0);
    for (double energy : {321.0, 1332.492, 11000.0}) {
        const auto outside = model.evaluate(
            energy, 0.85, uniform_edges(0.0, energy, 128));
        BOOST_TEST_CONTEXT("energy " << energy) {
            BOOST_TEST(k_escape_yield(outside) == 0.0);
            BOOST_TEST(
                outside.integrated_yield() == 1.0,
                boost::test_tools::tolerance(2e-12));
        }
    }
    // Each component is gated by its OWN support: k_beta_escape stops at
    // 175 keV while k_alpha_escape runs to 320, so a shared limit would let
    // the beta field extrapolate flat across most of a decade.
    const auto between = model.evaluate(
        250.0, 0.85, uniform_edges(0.0, 250.0, 128));
    BOOST_TEST(
        exact_component_yield(
            between, fs::SpectralComponent::DetectorKAlphaEscape) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            between, fs::SpectralComponent::DetectorKBetaEscape) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            between, fs::SpectralComponent::DetectorKAlphaAbsorbed) > 0.0);
    BOOST_TEST(
        exact_component_yield(
            between, fs::SpectralComponent::DetectorKBetaAbsorbed) == 0.0);

    // A payload written before the limit existed stores zero and keeps its
    // previous (ungated) behavior.
    payload.k_escape_support_keV = {};
    const fs::DetectorPulseResponseModel ungated(payload);
    BOOST_TEST(
        k_escape_yield(
            ungated.evaluate(1332.492, 0.85,
                             uniform_edges(0.0, 1332.492, 128))) > 0.0);
}

BOOST_AUTO_TEST_CASE(detector_block_basis_round_trip_and_normalization) {
    const auto payload = blocked_hpge_payload();
    const fs::DetectorPulseResponseModel original(payload);
    const std::string text = original.serialize();
    BOOST_TEST(text.rfind("CEELO_DETECTOR_RESPONSE 1\n", 0) == 0u);
    const auto restored =
        fs::DetectorPulseResponseModel::deserialize(text);
    BOOST_TEST(
        restored.payload().detector_single_compton.block_sizes[0] == 32u);
    BOOST_TEST(
        restored.payload().detector_pair_no_compton.block_sizes[3] == 24u);

    for (double energy : {42.5, 800.0, 1100.0, 2614.511, 5000.0, 12000.0}) {
        for (double cosine : {-1.0, -0.15, 0.5, 1.0}) {
            const auto edges = uniform_edges(0.0, energy, 256);
            const auto response = restored.evaluate(energy, cosine, edges);
            BOOST_TEST(
                response.integrated_yield() == 1.0,
                boost::test_tools::tolerance(2e-12));
            for (const auto& bin : response.continuous)
                BOOST_TEST(bin.integral >= 0.0);
            for (const auto& atom : response.atoms)
                BOOST_TEST(atom.yield >= 0.0);
        }
    }
    // The block boundaries move with incident energy; below threshold the
    // pair family must still vanish identically.
    const auto edges = uniform_edges(0.0, 1021.9979, 256);
    const auto at_threshold = restored.evaluate(1021.9979, 0.5, edges);
    BOOST_TEST(
        exact_component_yield(
            at_threshold,
            fs::SpectralComponent::DetectorPairNoCompton) == 0.0);
    BOOST_TEST(
        exact_component_yield(
            at_threshold,
            fs::SpectralComponent::DetectorPairFullEnergy) == 0.0);
}

BOOST_AUTO_TEST_CASE(detector_rate_surface_angle_warp_round_trip) {
    // The angle warp lives in both the C++ evaluator and the Python fitter and
    // has to agree exactly; a serialization slip that silently reverted the
    // warp to its identity setting is only visible off cos(theta) = 0, where
    // the warp does nothing, so check the shape of the detected rate across
    // angle rather than a single point.
    auto payload = blocked_hpge_payload();
    auto& surface = payload.detected_rate_surface;
    surface.energy_size = 24;
    surface.angle_size = 20;
    surface.angle_warp_exponent = 0.5;
    surface.angle_warp_epsilon = 0.03;
    for (std::size_t index = 0;
         index < std::size_t(surface.energy_size) * surface.angle_size;
         ++index) {
        surface.log_odds_coefficients[index] =
            -0.4 + 0.03 * std::sin(0.9 * index);
    }
    const fs::DetectorPulseResponseModel original(payload);
    const auto restored = fs::DetectorPulseResponseModel::deserialize(
        original.serialize());
    BOOST_TEST(restored.payload().detected_rate_surface.energy_size == 24u);
    BOOST_TEST(restored.payload().detected_rate_surface.angle_size == 20u);
    BOOST_TEST(
        restored.payload().detected_rate_surface.angle_warp_exponent == 0.5);

    const auto detected = [](const fs::DetectorPulseResponseModel& model,
                             double cosine) {
        const auto edges = uniform_edges(0.0, 2000.0, 8);
        const auto response = model.evaluate(2000.0, cosine, edges);
        return 1.0
            - exact_component_yield(
                response, fs::SpectralComponent::DetectorZeroDeposit);
    };
    for (double cosine : {-1.0, -0.97, -0.4, 0.0, 0.4, 0.97, 1.0}) {
        BOOST_TEST(
            detected(restored, cosine) == detected(original, cosine),
            boost::test_tools::tolerance(1e-12));
    }
    // With the warp active the detected rate must actually differ from the
    // unwarped surface away from cos(theta) = 0, otherwise a reverted warp
    // would pass unnoticed.
    auto identity = payload;
    identity.detected_rate_surface.angle_warp_exponent = 1.0;
    const fs::DetectorPulseResponseModel unwarped(identity);
    BOOST_TEST(
        std::abs(detected(restored, 0.97) - detected(unwarped, 0.97)) > 1e-6);
    BOOST_TEST(
        std::abs(detected(restored, 0.0) - detected(unwarped, 0.0)) < 1e-12);
}

BOOST_AUTO_TEST_CASE(detector_response_is_continuous_across_pair_landmarks) {
    // Two pair landmarks cross at Ein = 1533 keV (511 and Ein - 1022) and all
    // of them collapse at the 1022 keV threshold. The block partition is
    // defined on the sorted landmarks so the degenerate interval carries no
    // mass; the emitted spectrum must not jump as incident energy sweeps
    // through either point.
    const fs::DetectorPulseResponseModel model(blocked_hpge_payload());
    const auto shape_at = [&](double energy) {
        // Compare on a fixed fractional grid so the two spectra are directly
        // commensurable despite their different support.
        const auto edges = uniform_edges(0.0, energy, 256);
        const auto response = model.evaluate(energy, 0.35, edges);
        std::vector<double> profile(256, 0.0);
        for (const auto& bin : response.continuous) {
            const std::size_t index = static_cast<std::size_t>(
                256.0 * bin.lower_keV / energy + 0.5);
            if (index < profile.size())
                profile[index] += bin.integral;
        }
        double total = 0.0;
        for (double value : profile)
            total += value;
        if (total > 0.0)
            for (double& value : profile)
                value /= total;
        return profile;
    };
    for (double landmark : {1021.9979, 1532.99685}) {
        const auto below = shape_at(landmark * (1.0 - 1e-4));
        const auto above = shape_at(landmark * (1.0 + 1e-4));
        double variation = 0.0;
        for (std::size_t index = 0; index < below.size(); ++index)
            variation += 0.5 * std::abs(below[index] - above[index]);
        BOOST_TEST(variation < 0.02);
    }
}

BOOST_AUTO_TEST_CASE(detector_continuum_has_no_knot_scale_ripple) {
    // The emitted continuum should be smooth between physical landmarks. A
    // basis artifact - the failure mode that put a spurious feature near
    // 4650 keV into an earlier payload - shows up as curvature on the scale
    // of the knot spacing rather than of the spectrum.
    const fs::DetectorPulseResponseModel model(blocked_hpge_payload());
    const double energy = 5000.0;
    const std::size_t bins = 512;
    const auto edges = uniform_edges(0.0, energy, bins);
    const auto response = model.evaluate(energy, 0.97, edges);
    std::vector<double> profile(bins, 0.0);
    for (const auto& bin : response.continuous) {
        const std::size_t index = static_cast<std::size_t>(
            bins * bin.lower_keV / energy + 0.5);
        if (index < profile.size())
            profile[index] += bin.integral;
    }
    const double bin_width = energy / bins;
    const double landmarks[] = {
        510.99895, energy - 1021.9979, energy - 510.99895, energy,
        energy * (1.0 - 1.0 / (1.0 + 2.0 * energy / 510.99895))};
    double worst = 0.0;
    for (std::size_t index = 1; index + 1 < bins; ++index) {
        const double center = (index + 0.5) * bin_width;
        bool near_landmark = false;
        for (double landmark : landmarks)
            near_landmark |= std::abs(center - landmark) <= 2.0 * bin_width;
        if (near_landmark)
            continue;
        const double local =
            (profile[index - 1] + profile[index] + profile[index + 1]) / 3.0;
        if (!(local > 1e-12))
            continue;
        const double curvature = std::abs(
            profile[index - 1] - 2.0 * profile[index] + profile[index + 1]);
        worst = std::max(worst, curvature / local);
    }
    BOOST_TEST(worst < 0.35);
}

BOOST_AUTO_TEST_CASE(finite_distance_factor_limits) {
    using fs::detail::finite_distance_factor;
    // The far field is the payload's own normalization.
    BOOST_TEST(finite_distance_factor(1e6, 2.0, 4.0) == 1.0,
               boost::test_tools::tolerance(1e-5));
    // Far from the detector it is (d / (d + delta))^2.
    BOOST_TEST(finite_distance_factor(1000.0, 2.0, 4.0)
                   == std::pow(1000.0 / 1002.0, 2),
               boost::test_tools::tolerance(2e-5));
    // With no depth it is the disk's solid angle over its far-field value.
    const double d = 10.0, R = 4.0;
    const double omega = 0.5 * (1.0 - d / std::sqrt(d * d + R * R));
    BOOST_TEST(finite_distance_factor(d, 0.0, R) == omega / (R * R / (4.0 * d * d)),
               boost::test_tools::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(finite_distance_depths_round_trip) {
    auto detector = blocked_hpge_payload();
    detector.finite_distance = {{60.0, 662.0, 2614.5}, {0.24, 2.5, 3.2},
                                {3.7, 1.9, 0.7}, {0.235, 2.3, 2.9}, {3.6, 1.9, 0.0}};
    const auto restored = fs::DetectorPulseResponseModel::deserialize(
        fs::DetectorPulseResponseModel(detector).serialize()).payload();
    BOOST_TEST(restored.finite_distance.full_energy_radius_cm
                   == detector.finite_distance.full_energy_radius_cm,
               boost::test_tools::per_element());
    BOOST_TEST(restored.finite_distance.energies_keV
                   == detector.finite_distance.energies_keV,
               boost::test_tools::per_element());
    BOOST_TEST(restored.finite_distance.full_energy_cm
                   == detector.finite_distance.full_energy_cm,
               boost::test_tools::per_element());
    BOOST_TEST(restored.finite_distance.detected_cm
                   == detector.finite_distance.detected_cm,
               boost::test_tools::per_element());
    fs::FullSpectrumPayload composed;
    composed.detector = detector;
    const auto embedded = fs::FastSpectrumModel::deserialize(
        fs::FastSpectrumModel(composed).serialize()).payload();
    BOOST_TEST(embedded.detector.finite_distance.detected_cm
                   == detector.finite_distance.detected_cm,
               boost::test_tools::per_element());
    // None measured: an empty record, and predict's far-field limit.
    const auto none = fs::DetectorPulseResponseModel::deserialize(
        fs::DetectorPulseResponseModel(blocked_hpge_payload()).serialize()).payload();
    BOOST_TEST(!none.finite_distance.enabled());
}

BOOST_AUTO_TEST_CASE(finite_distance_depths_scale_predict) {
    const std::vector<fs::EmissionLine> lines{{661.657, 1.0}};
    const fs::ShieldParameters<double> bare{0.0, 26.0};
    const fs::ResolutionParameters<double> narrow{0.01, 0.0, 0.0};
    const auto channels = uniform_edges(0.0, 800.0, 8000);
    fs::FullSpectrumPayload payload;
    payload.distance_cm = 30.0;
    const auto far = fs::FastSpectrumModel(payload).predict(
        lines, 1.0, bare, narrow, channels);
    auto peak = [&](const std::vector<double>& counts) {
        double s = 0.0;
        for (std::size_t i = 0; i < counts.size(); ++i)
            if (std::abs(0.5 * (channels[i] + channels[i + 1]) - 661.657) < 0.3)
                s += counts[i];
        return s;
    };
    using fs::detail::finite_distance_factor;
    // One depth for everything: every count scales by its factor.
    payload.detector.finite_distance = {{100.0, 1000.0}, {2.0, 2.0}, {4.0, 4.0},
                                        {2.0, 2.0}, {4.0, 4.0}};
    const auto near = fs::FastSpectrumModel(payload).predict(
        lines, 1.0, bare, narrow, channels);
    BOOST_TEST(sum(near) == sum(far) * finite_distance_factor(30.0, 2.0, 4.0),
               boost::test_tools::tolerance(1e-9));
    // A deeper full-energy depth: the peak takes its factor, the detected
    // total keeps the detected depth's.
    payload.detector.finite_distance.full_energy_cm = {3.0, 3.0};
    const auto split = fs::FastSpectrumModel(payload).predict(
        lines, 1.0, bare, narrow, channels);
    BOOST_TEST(sum(split) == sum(far) * finite_distance_factor(30.0, 2.0, 4.0),
               boost::test_tools::tolerance(1e-9));
    BOOST_TEST(peak(split) == peak(far) * finite_distance_factor(30.0, 3.0, 4.0),
               boost::test_tools::tolerance(1e-3));
    // Measured on the axis: side-on, predict keeps the far-field acceptance.
    payload.incidence_cos_theta = 0.0;
    fs::FullSpectrumPayload side = payload;
    side.detector.finite_distance = fs::FiniteDistanceDepth{};
    BOOST_TEST(sum(fs::FastSpectrumModel(payload).predict(lines, 1.0, bare, narrow, channels))
                   == sum(fs::FastSpectrumModel(side).predict(lines, 1.0, bare, narrow, channels)),
               boost::test_tools::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(fixed_energy_detector_atoms_are_not_box_deposited) {
    // predict spreads an atom over a source-quadrature bin only when the atom
    // moves with the incident energy; the absorbed K lines are dead-layer K
    // x-rays at fixed energies, like the 511 keV line, and must stay points.
    using C = fs::SpectralComponent;
    for (C c : {C::DetectorFullEnergy, C::DetectorPairFullEnergy,
                C::DetectorSingleEscape, C::DetectorDoubleEscape,
                C::DetectorKAlphaEscape, C::DetectorKBetaEscape})
        BOOST_TEST(fs::detail::detector_atom_tracks_incident(c));
    for (C c : {C::DetectorZeroDeposit, C::DetectorAnnihilation511,
                C::DetectorKAlphaAbsorbed, C::DetectorKBetaAbsorbed})
        BOOST_TEST(!fs::detail::detector_atom_tracks_incident(c));
}

BOOST_AUTO_TEST_CASE(resolution_convolution_conserves_expected_counts) {
    fs::FullSpectrumPayload payload;
    payload.source_quadrature_bins = 8;
    fs::FastSpectrumModel model(payload);
    const std::vector<fs::EmissionLine> lines{{661.657, 0.851}};
    const fs::ShieldParameters<double> shield{0.5, 26.0};
    const auto channels = uniform_edges(-200.0, 1200.0, 700);
    const fs::ResolutionParameters<double> narrow{0.01, 0.0, 0.0};
    const fs::ResolutionParameters<double> realistic{4.0, 2.0, 2e-4};
    const auto a = model.predict(lines, 1e6, shield, narrow, channels);
    const auto b = model.predict(lines, 1e6, shield, realistic, channels);
    BOOST_TEST(sum(a) == sum(b), boost::test_tools::tolerance(2e-12));
}

BOOST_AUTO_TEST_CASE(source_quadrature_refinement_converges) {
    // Guards against depositing the escaping source continuum as nq delta
    // functions at quadrature bin centres: under heavy shielding the predicted
    // continuum shape then never converges in nq (TV(nq, 2nq) stuck at 15-23%
    // through nq = 64) even though total counts do. The bin yields must be
    // spread across the bin width.
    const std::vector<fs::EmissionLine> lines{{661.657, 0.851}};
    const fs::ShieldParameters<double> shield{7.874, 26.0};
    const auto channels = uniform_edges(0.0, 800.0, 800);
    const fs::ResolutionParameters<double> resolution{0.21, 0.0023, 0.0};
    auto spectrum_at = [&](int quadrature_bins) {
        fs::FullSpectrumPayload payload;
        payload.source_quadrature_bins = quadrature_bins;
        const fs::FastSpectrumModel model(payload);
        return model.predict(lines, 1e6, shield, resolution, channels);
    };
    const auto coarse = spectrum_at(8);
    const auto fine = spectrum_at(64);
    const double coarse_total = sum(coarse);
    const double fine_total = sum(fine);
    BOOST_TEST(coarse_total == fine_total,
               boost::test_tools::tolerance(5e-3));
    double variation = 0.0;
    for (std::size_t i = 0; i < coarse.size(); ++i)
        variation += 0.5 * std::abs(
            coarse[i] / coarse_total - fine[i] / fine_total);
    BOOST_TEST(variation < 0.03);
}

BOOST_AUTO_TEST_CASE(payload_round_trip_is_prediction_identical) {
    fs::FullSpectrumPayload payload;
    payload.model_name = "round-trip fixture";
    payload.training_hash = "sha256:test";
    payload.statistical_uncertainty = "component sumw2";
    payload.source.single_escape_scale = 0.812345678901234;
    auto& surrogate = payload.source.multiple_scatter_surrogate;
    surrogate.use_cubic_bspline = true;
    surrogate.higher_yield.degree = 3;
    surrogate.higher_yield.bandwidth = {0.08, 0.025};
    fs::MultipleScatterYieldNode yield_node;
    yield_node.feature = {0.67, 0.32};
    yield_node.coefficients = {-2.1, 0.3, -0.08, 0.01, 0.0, 0.0};
    surrogate.higher_yield.nodes.push_back(yield_node);
    surrogate.higher_compton.yield_bandwidth = {0.07, 0.03, 0.04};
    surrogate.higher_compton.shape_bandwidth = {0.10, 0.02, 0.03};
    fs::MultipleScatterRbfNode shape_node;
    shape_node.feature = {0.67, 0.32, 0.25};
    shape_node.log_conditional_yield = 0.0;
    shape_node.shape_coefficients.fill(
        1.0 / fs::kMultipleScatterBasisSize);
    surrogate.higher_compton.nodes.push_back(shape_node);
    auto& secondary = payload.source.secondary_surrogate;
    secondary.pair_511_atom.yield_bandwidth = {0.03, 0.05, 0.03};
    secondary.pair_511_atom.shape_bandwidth = {0.04, 0.07, 0.04};
    fs::MultipleScatterRbfNode pair_atom_node;
    pair_atom_node.feature = {0.77, 35.0 / 92.0, 0.22};
    pair_atom_node.log_conditional_yield = -0.08;
    pair_atom_node.shape_coefficients.fill(
        1.0 / fs::kMultipleScatterBasisSize);
    secondary.pair_511_atom.nodes.push_back(pair_atom_node);

    fs::MultipleScatterRbfNode pair_continuum_node = pair_atom_node;
    pair_continuum_node.log_conditional_yield = 0.13;
    for (std::size_t k = 0;
         k < fs::kMultipleScatterBasisSize; ++k) {
        pair_continuum_node.shape_coefficients[k] =
            0.1 + static_cast<double>(
                (3 * k + 1) % fs::kMultipleScatterBasisSize);
    }
    secondary.pair_annihilation_continuum.nodes.push_back(
        pair_continuum_node);
    fs::MultipleScatterRbfNode pair_in_flight_node =
        pair_continuum_node;
    pair_in_flight_node.log_conditional_yield = -1.3;
    secondary.pair_annihilation_in_flight.nodes.push_back(
        pair_in_flight_node);

    fs::MultipleScatterRbfNode compton_brems_node = pair_atom_node;
    compton_brems_node.log_conditional_yield = -3.2;
    secondary.compton_electron_brems.nodes.push_back(
        compton_brems_node);
    payload.source.enable_source_electron_brems = true;
    payload.source.enable_compton_doppler = false;
    payload.source.enable_photonuclear = true;
    fs::PhotonuclearCrossSectionRow photonuclear_lower;
    photonuclear_lower.energy_keV = 1000.0;
    fs::PhotonuclearCrossSectionRow photonuclear_upper;
    photonuclear_upper.energy_keV = 2000.0;
    for (std::size_t index = 0;
         index < photonuclear_lower.mass_attenuation_cm2_g.size();
         ++index) {
        photonuclear_lower.mass_attenuation_cm2_g[index] =
            1e-5 * static_cast<double>(index + 1);
        photonuclear_upper.mass_attenuation_cm2_g[index] =
            2e-5 * static_cast<double>(index + 1);
    }
    payload.source.photonuclear.cross_sections = {
        photonuclear_lower, photonuclear_upper};
    payload.source.photonuclear.escaped_yield_scale = 0.73;
    payload.detector = fs::coaxial_hpge_payload();
    fs::FastSpectrumModel original(payload);
    const fs::FastSpectrumModel decoded =
        fs::FastSpectrumModel::deserialize(original.serialize());
    const auto channels = uniform_edges(0.0, 1600.0, 256);
    const std::vector<fs::EmissionLine> lines{
        {1173.228, 0.999}, {1332.492, 0.9998}};
    const fs::ShieldParameters<double> shield{1.2, 29.4};
    const fs::ResolutionParameters<double> resolution{0.4, 0.001, 1e-6};
    const auto expected =
        original.predict(lines, 10000.0, shield, resolution, channels);
    const auto actual =
        decoded.predict(lines, 10000.0, shield, resolution, channels);
    BOOST_REQUIRE_EQUAL(expected.size(), actual.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        BOOST_TEST(actual[i] == expected[i],
                   boost::test_tools::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(hydrogen_mixture_cross_sections_are_exact_at_endpoints) {
    const fs::SmoothElementCrossSection cross_sections;
    for (const double energy : {30.0, 300.0, 1500.0, 6000.0}) {
        const auto blend = cross_sections.evaluate(energy, 26.0);
        const auto at_zero = cross_sections.evaluate(energy, 26.0, 0.0);
        for (std::size_t process = 0; process < blend.size(); ++process) {
            // f_H = 0 must reproduce the pure effective-Z blend exactly:
            // every existing payload and test rides on this identity.
            BOOST_TEST(at_zero[process] == blend[process]);
        }
        const auto hydrogen =
            fs::element_mass_cross_sections(1, energy);
        const std::array<double, 4> hydrogen_values{
            hydrogen.pe, hydrogen.cs, hydrogen.rs, hydrogen.pp};
        const auto at_one = cross_sections.evaluate(energy, 26.0, 1.0);
        for (std::size_t process = 0; process < blend.size(); ++process) {
            // The f_H = 1 endpoint carries the roundoff of the
            // blend + f*(H - blend) form (blend >> H cancellation); only the
            // f_H = 0 identity is exact by construction.
            BOOST_TEST(at_one[process] == hydrogen_values[process],
                       boost::test_tools::tolerance(1e-9));
        }
        // Water's hydrogen fraction sits strictly between the endpoints and
        // raises Compton per gram (hydrogen has ~2x the electrons per gram).
        const auto water_like =
            cross_sections.evaluate(energy, 8.0, 0.1119);
        const auto oxygen_only = cross_sections.evaluate(energy, 8.0, 0.0);
        BOOST_TEST(water_like[1] > oxygen_only[1]);
    }

    // Angular and Doppler densities: f_H = 0 exact, mixture stays a density.
    const double cosine = 0.3;
    for (const double energy : {100.0, 662.0}) {
        const auto blend_density =
            cross_sections.compton_angular_density(energy, 26.0, cosine);
        BOOST_TEST(
            cross_sections.compton_angular_density(
                energy, 26.0, 0.0, cosine) == blend_density);
        const auto mixed = cross_sections.compton_angular_density(
            energy, 26.0, 0.35, cosine);
        BOOST_TEST(mixed > 0.0);
        const double lower = 0.8 * energy;
        const double upper = 0.85 * energy;
        BOOST_TEST(
            cross_sections.compton_doppler_bin_probability(
                energy, 26.0, 0.0, cosine, lower, upper)
            == cross_sections.compton_doppler_bin_probability(
                   energy, 26.0, cosine, lower, upper));
        BOOST_TEST(
            cross_sections.compton_doppler_bin_probability(
                energy, 26.0, 0.35, cosine, lower, upper) >= 0.0);
    }
}

BOOST_AUTO_TEST_CASE(hydrogen_fraction_shifts_source_measure_physically) {
    // Oxygen-with-hydrogen (water-like) vs pure effective-Z=8 at the same
    // areal density: hydrogen adds electrons per gram, so more Compton
    // removal (smaller uncollided line) and more scattered yield.
    const fs::SourceGeomFluxModel model;
    const auto bins = uniform_edges(0.0, 662.0, 64);
    const auto plain = model.evaluate(661.657, 10.0, 8.0, bins);
    const auto watery = model.evaluate(661.657, 10.0, 8.0, 0.1119, bins);
    const auto forwarded = model.evaluate(661.657, 10.0, 8.0, 0.0, bins);

    auto atom_yield = [](const auto& measure, fs::SpectralComponent kind) {
        double total = 0.0;
        for (const auto& atom : measure.atoms)
            if (atom.component == kind)
                total += atom.yield;
        return total;
    };
    auto continuum_yield = [](const auto& measure) {
        double total = 0.0;
        for (const auto& bin : measure.continuous)
            total += bin.integral;
        return total;
    };

    // f_H = 0 through the explicit argument must match the 4-arg overload
    // exactly (the forwarding identity everything else relies on).
    BOOST_REQUIRE_EQUAL(plain.atoms.size(), forwarded.atoms.size());
    for (std::size_t i = 0; i < plain.atoms.size(); ++i)
        BOOST_TEST(plain.atoms[i].yield == forwarded.atoms[i].yield);
    BOOST_REQUIRE_EQUAL(
        plain.continuous.size(), forwarded.continuous.size());
    for (std::size_t i = 0; i < plain.continuous.size(); ++i)
        BOOST_TEST(plain.continuous[i].integral
                   == forwarded.continuous[i].integral);

    const double plain_direct =
        atom_yield(plain, fs::SpectralComponent::Uncollided);
    const double watery_direct =
        atom_yield(watery, fs::SpectralComponent::Uncollided);
    BOOST_TEST(watery_direct < plain_direct);
    BOOST_TEST(watery_direct > 0.0);
    BOOST_TEST(continuum_yield(watery) > continuum_yield(plain));
}

BOOST_AUTO_TEST_CASE(source_escape_payload_round_trip_is_exact) {
    fs::SourceModelPayload source;
    source.single_escape_scale = 0.812345678901234;
    source.enable_source_electron_brems = true;
    source.enable_compton_doppler = false;
    source.enable_photonuclear = true;
    auto& surrogate = source.multiple_scatter_surrogate;
    surrogate.use_cubic_bspline = true;
    surrogate.higher_yield.degree = 3;
    surrogate.higher_yield.bandwidth = {0.08, 0.025};
    fs::MultipleScatterYieldNode yield_node;
    yield_node.feature = {0.67, 0.32};
    yield_node.coefficients = {-2.1, 0.3, -0.08, 0.01, 0.0, 0.0};
    surrogate.higher_yield.nodes.push_back(yield_node);
    surrogate.higher_compton.yield_bandwidth = {0.07, 0.03, 0.04};
    surrogate.higher_compton.shape_bandwidth = {0.10, 0.02, 0.03};
    fs::MultipleScatterRbfNode shape_node;
    shape_node.feature = {0.67, 0.32, 0.25};
    shape_node.log_conditional_yield = 0.0;
    shape_node.shape_coefficients.fill(
        1.0 / fs::kMultipleScatterBasisSize);
    surrogate.higher_compton.nodes.push_back(shape_node);
    auto& secondary = source.secondary_surrogate;
    secondary.pair_511_atom.yield_bandwidth = {0.03, 0.05, 0.03};
    secondary.pair_511_atom.shape_bandwidth = {0.04, 0.07, 0.04};
    fs::MultipleScatterRbfNode pair_atom_node;
    pair_atom_node.feature = {0.77, 35.0 / 92.0, 0.22};
    pair_atom_node.log_conditional_yield = -0.08;
    pair_atom_node.shape_coefficients.fill(
        1.0 / fs::kMultipleScatterBasisSize);
    secondary.pair_511_atom.nodes.push_back(pair_atom_node);
    fs::MultipleScatterRbfNode brems_node = pair_atom_node;
    brems_node.log_conditional_yield = -3.2;
    secondary.compton_electron_brems.nodes.push_back(brems_node);
    fs::PhotonuclearCrossSectionRow photonuclear_lower;
    photonuclear_lower.energy_keV = 1000.0;
    fs::PhotonuclearCrossSectionRow photonuclear_upper;
    photonuclear_upper.energy_keV = 2000.0;
    for (std::size_t index = 0;
         index < photonuclear_lower.mass_attenuation_cm2_g.size();
         ++index) {
        photonuclear_lower.mass_attenuation_cm2_g[index] =
            1e-5 * static_cast<double>(index + 1);
        photonuclear_upper.mass_attenuation_cm2_g[index] =
            2e-5 * static_cast<double>(index + 1);
    }
    source.photonuclear.cross_sections = {
        photonuclear_lower, photonuclear_upper};
    source.photonuclear.escaped_yield_scale = 0.73;

    const fs::SourceGeomFluxModel original(source);
    const std::string text = original.serialize();
    BOOST_TEST(text.rfind("CEELO_SOURCE_ESCAPE 1", 0) == 0u);
    const fs::SourceGeomFluxModel decoded =
        fs::SourceGeomFluxModel::deserialize(text);
    // Serialization must round-trip exactly.
    BOOST_TEST(decoded.serialize() == text);

    const auto bins = uniform_edges(0.0, 2700.0, 128);
    const auto expected = original.evaluate(2614.5, 8.0, 26.0, bins);
    const auto actual = decoded.evaluate(2614.5, 8.0, 26.0, bins);
    BOOST_REQUIRE_EQUAL(expected.atoms.size(), actual.atoms.size());
    for (std::size_t i = 0; i < expected.atoms.size(); ++i) {
        BOOST_TEST(actual.atoms[i].yield == expected.atoms[i].yield,
                   boost::test_tools::tolerance(1e-14));
    }
    BOOST_REQUIRE_EQUAL(
        expected.continuous.size(), actual.continuous.size());
    for (std::size_t i = 0; i < expected.continuous.size(); ++i) {
        BOOST_TEST(
            actual.continuous[i].integral
                == expected.continuous[i].integral,
            boost::test_tools::tolerance(1e-14));
    }
}

namespace {

/// A one-node, log-yield-0, 7-D fluorescence surrogate: its family yield is
/// exactly the mechanistic base.
fs::MultipleScatterRbfComponent unit_fluorescence_surrogate() {
    fs::MultipleScatterRbfComponent component;
    component.feature_dimension = fs::kSourceRbfMaxFeatures;
    fs::MultipleScatterRbfNode node;
    node.feature.fill(0.5);
    node.log_conditional_yield = 0.0;
    node.shape_coefficients.fill(1.0 / fs::kMultipleScatterBasisSize);
    component.nodes.push_back(node);
    return component;
}

}  // namespace

BOOST_AUTO_TEST_CASE(lead_fluorescence_lines_conserve_their_yields) {
    const auto& pb = fs::element_fluorescence(82);
    std::array<int, 3> count{};
    double k_lines = 0.0, l3 = 0.0;
    for (const auto& line : pb.lines) {
        ++count[line.family];
        BOOST_TEST(line.energy_keV >= fs::kFluorescenceMinimumKeV);
        if (line.family == 0)
            BOOST_TEST((line.energy_keV > 72.0 && line.energy_keV < 75.5));
        if (line.family == 1)
            BOOST_TEST((line.energy_keV > 84.0 && line.energy_keV < 88.1));
        if (line.family == 2)
            BOOST_TEST((line.energy_keV > 9.0 && line.energy_keV < 16.0));
        k_lines += line.family < 2 ? line.per_k : 0.0;
        l3 += line.per_l[2];
    }
    BOOST_TEST(count[0] >= 2);
    BOOST_TEST(count[1] >= 1);
    BOOST_TEST(count[2] >= 3);
    // Compression folds the weak lines into kept ones: the K x-ray yield per
    // K vacancy is still omega_K, an L3 vacancy's still its retained yield.
    const auto& xs = ceelo::CrossSectionData::instance();
    BOOST_TEST(k_lines == xs.fluorescence(82)->fluorescence_yield,
               boost::test_tools::tolerance(1e-6));
    BOOST_TEST(l3 == xs.l_fluorescence(82)->sub[2].fluorescence_yield,
               boost::test_tools::tolerance(1e-6));
    BOOST_TEST(fs::element_fluorescence(26).lines.empty());  // Fe K lines < 10 keV
}

BOOST_AUTO_TEST_CASE(fluorescence_bases_follow_the_edges_and_become_atoms) {
    const fs::SmoothElementCrossSection xs(0.42);
    const fs::SourceSecondarySurrogatePayload options;
    for (double energy : {50.0, 200.0, 2614.5}) {
        const auto mu0 = xs.evaluate(energy, 82.0, 0.0);
        const double mu0_total = mu0[0] + mu0[1] + mu0[2] + mu0[3];
        const auto bases = fs::detail::source_component_bases(
            energy, 1.0, 82.0, 0.0, mu0, mu0_total, xs, options);
        const double k_alpha = bases.base[fs::detail::kFluorescenceKAlpha];
        const double k_beta = bases.base[fs::detail::kFluorescenceKBeta];
        const double l = bases.base[fs::detail::kFluorescenceL];
        if (energy < 88.0) {   // below the Pb K edge: L x-rays only
            BOOST_TEST(k_alpha == 0.0);
            BOOST_TEST(k_beta == 0.0);
        } else {
            BOOST_TEST(k_alpha > k_beta);
            BOOST_TEST(k_beta > 0.0);
        }
        BOOST_TEST(l > 0.0);
        for (int family = 0; family < fs::kFluorescenceFamilies; ++family) {
            double sum = 0.0;
            for (const auto& line : bases.fluorescence_lines[family])
                sum += line.second;
            BOOST_TEST(sum == bases.base[fs::detail::kFluorescenceKAlpha + family],
                       boost::test_tools::tolerance(1e-12));
        }
    }

    // A unit surrogate turns each family's base into its lines' atoms.
    fs::SourceModelPayload payload;
    // One feature width for every secondary component, as the payload format
    // carries it once.
    for (auto* component : payload.secondary_surrogate.components())
        component->feature_dimension = fs::kSourceRbfMaxFeatures;
    payload.secondary_surrogate.fluorescence_k_alpha_atom = unit_fluorescence_surrogate();
    payload.secondary_surrogate.fluorescence_k_beta_atom = unit_fluorescence_surrogate();
    payload.secondary_surrogate.fluorescence_l_atom = unit_fluorescence_surrogate();
    const fs::SourceGeomFluxModel model(payload);
    const auto bins = uniform_edges(0.0, 700.0, 64);
    const auto measure = model.evaluate(661.657, 1.0, 82.0, bins);
    const auto mu0 = xs.evaluate(661.657, 82.0, 0.0);
    const auto bases = fs::detail::source_component_bases(
        661.657, 1.0, 82.0, 0.0, mu0, mu0[0] + mu0[1] + mu0[2] + mu0[3], xs, options);
    BOOST_TEST(measure.component_yield(fs::SpectralComponent::SourceFluorescenceKAlpha)
                   == bases.base[fs::detail::kFluorescenceKAlpha],
               boost::test_tools::tolerance(1e-12));
    BOOST_TEST(measure.component_yield(fs::SpectralComponent::SourceFluorescence)
                   == bases.base[fs::detail::kFluorescenceKAlpha]
                      + bases.base[fs::detail::kFluorescenceKBeta]
                      + bases.base[fs::detail::kFluorescenceL],
               boost::test_tools::tolerance(1e-12));
    for (const auto& atom : measure.atoms)
        if (atom.component == fs::SpectralComponent::SourceFluorescenceKAlpha)
            BOOST_TEST((atom.energy_keV > 72.0 && atom.energy_keV < 75.5));

    // Off: no fluorescence atoms; and the payload round trip keeps them.
    payload.enable_fluorescence = false;
    BOOST_TEST(fs::SourceGeomFluxModel(payload)
                   .evaluate(661.657, 1.0, 82.0, bins)
                   .component_yield(fs::SpectralComponent::SourceFluorescence) == 0.0);
    payload.enable_fluorescence = true;
    const std::string text = model.serialize();
    BOOST_TEST(fs::SourceGeomFluxModel::deserialize(text).serialize() == text);
}

BOOST_AUTO_TEST_CASE(charged_fluorescence_escapes_from_where_photons_are_absorbed) {
    const fs::SmoothElementCrossSection xs(0.42);
    const double range = fs::detail::charged_blend_csda_range(
        5000.0, 82.0, 0.0, false, xs.blend_width_Z());
    const auto at = [&](double range_ratio, bool positron) {
        return fs::detail::charged_fluorescence_bases(
            5000.0, range_ratio * range, 82.0, 0.0, positron, xs);
    };
    for (const bool positron : {false, true}) {
        const auto bases = at(1.0, positron);
        std::array<double, fs::kFluorescenceFamilies> a{}, b{};
        for (const auto& line : bases.lines) {
            a[line.family] += line.a;
            b[line.family] += line.b;
        }
        for (int family = 0; family < fs::kFluorescenceFamilies; ++family) {
            BOOST_TEST(bases.a[family] > 0.0);
            BOOST_TEST(a[family] == bases.a[family], boost::test_tools::tolerance(1e-12));
            BOOST_TEST(b[family] == bases.b[family], boost::test_tools::tolerance(1e-12));
            BOOST_TEST((bases.b[family] > 0.0) == positron);
        }
        // Past the stop radius the photons born on the path still reach the
        // surface: from 2 to 4 ranges of lead the K x-rays fall by the
        // attenuation of the most penetrating bremsstrahlung, not by the
        // line's own, which would be ~exp(-2.4 x 10).
        const double two = at(2.0, positron).a[0], four = at(4.0, positron).a[0];
        BOOST_TEST(four < two);
        BOOST_TEST(four > 0.05 * two);
    }
    // Calcium's K lines are below the 10 keV floor: nothing to emit.
    const auto calcium = fs::detail::charged_fluorescence_bases(
        1000.0, 0.1, 20.0, 0.0, true, xs);
    BOOST_TEST(calcium.lines.empty());
}

BOOST_AUTO_TEST_CASE(source_truth_labels_do_not_change_transport_stream) {
    ceelo::Material lead = ceelo::make_Lead();
    ceelo::SourceGeometry plain;
    plain.configure_point(Eigen::Vector3d::Zero());
    plain.add_shield(&lead, 0.2);
    plain.set_source_brems(false);

    ceelo::SourceGeometry labelled = plain;
    labelled.set_record_truth_labels(true);
    std::mt19937_64 rng_plain(987654321);
    std::mt19937_64 rng_labelled(987654321);
    const Eigen::Vector3d direction(0.2, -0.3, std::sqrt(0.87));
    for (int event = 0; event < 1000; ++event) {
        const auto a = plain.transport_source_photon(
            Eigen::Vector3d::Zero(), direction, 1332.5, rng_plain);
        const auto b = labelled.transport_source_photon(
            Eigen::Vector3d::Zero(), direction, 1332.5, rng_labelled);
        BOOST_TEST(a.survived == b.survived);
        BOOST_TEST(a.energy_keV == b.energy_keV);
        BOOST_TEST(a.position.x() == b.position.x());
        BOOST_TEST(a.position.y() == b.position.y());
        BOOST_TEST(a.position.z() == b.position.z());
        BOOST_TEST(a.direction.x() == b.direction.x());
        BOOST_TEST(a.direction.y() == b.direction.y());
        BOOST_TEST(a.direction.z() == b.direction.z());
        BOOST_TEST(a.secondaries.size() == b.secondaries.size());
    }
    BOOST_TEST(rng_plain() == rng_labelled());
}

BOOST_AUTO_TEST_CASE(source_transport_can_remove_rayleigh_competition) {
    ceelo::Material iron = ceelo::make_Iron();
    ceelo::SourceGeometry source;
    BOOST_TEST(source.rayleigh_enabled());
    source.configure_point(Eigen::Vector3d::Zero());
    source.add_shield(&iron, 1.0);
    source.set_source_brems(false);
    source.set_record_truth_labels(true);
    source.set_rayleigh_enabled(false);

    constexpr double energy_keV = 300.0;
    const auto xs = iron.macroscopic_xs(energy_keV * 1e-3);
    const double expected = std::exp(
        -(xs.mu_pe + xs.mu_cs + xs.mu_pp));
    const double direct = source.no_interaction_probability(
        Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(), energy_keV);
    BOOST_TEST(direct == expected, boost::test_tools::tolerance(1e-13));

    std::mt19937_64 rng(20260728);
    for (int event = 0; event < 2000; ++event) {
        const auto result = source.transport_source_photon(
            Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(),
            energy_keV, rng);
        BOOST_TEST(result.truth.rayleigh_interactions == 0);
    }
}

BOOST_AUTO_TEST_CASE(source_positron_transport_labels_in_flight_annihilation) {
    ceelo::Material lead_equivalent(
        "lead-equivalent", 1.0,
        {{static_cast<std::uint8_t>(82), 1.0}});
    ceelo::SourceGeometry source;
    source.configure_point(Eigen::Vector3d::Zero());
    source.add_shield(&lead_equivalent, 20.134);
    source.set_source_brems(false);
    source.set_record_truth_labels(true);
    source.set_source_positron_transport(true);

    bool saw_at_rest = false;
    bool saw_in_flight = false;
    std::mt19937_64 rng(120008201);
    for (int event = 0; event < 5000; ++event) {
        const auto result = source.transport_source_photon(
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d::UnitZ(), 12000.0, rng);
        for (const auto& secondary : result.secondaries) {
            using Origin =
                ceelo::SourceGeometry::SourceSecondaryPhoton::Origin;
            saw_at_rest = saw_at_rest
                || secondary.origin == Origin::PairAnnihilation;
            saw_in_flight = saw_in_flight
                || secondary.origin == Origin::PairAnnihilationInFlight;
        }
    }
    BOOST_TEST(saw_at_rest);
    BOOST_TEST(saw_in_flight);
}
