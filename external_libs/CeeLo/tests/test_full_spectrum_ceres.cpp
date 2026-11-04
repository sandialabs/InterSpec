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

#define BOOST_TEST_MODULE FullSpectrumApproxCeresTests
#include <boost/test/unit_test.hpp>

#include "full_spectrum/FullSpectrumModel.h"

#include <ceres/jet.h>

#include <array>
#include <cmath>
#include <type_traits>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = ceelo::full_spectrum;

namespace {

std::vector<double> edges() {
    std::vector<double> result(129);
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = 1600.0 * static_cast<double>(i) / 128.0;
    return result;
}

template <typename Scalar>
Scalar objective(const std::array<Scalar, 6>& p) {
    fs::FullSpectrumPayload payload;
    payload.source_quadrature_bins = 6;
    payload.detector.resolution_half_window_channels = 64;
    auto& surrogate = payload.source.multiple_scatter_surrogate;
    surrogate.use_cubic_bspline = true;
    surrogate.higher_yield.degree = 3;
    surrogate.higher_yield.bandwidth = {0.20, 0.20};
    surrogate.higher_compton.yield_bandwidth = {0.20, 0.20, 0.20};
    surrogate.higher_compton.shape_bandwidth = {0.20, 0.20, 0.20};
    for (int index = 0; index < 3; ++index) {
        fs::MultipleScatterYieldNode yield_node;
        yield_node.feature = {
            0.35 + 0.22 * index, 0.12 + 0.25 * index};
        yield_node.coefficients = {
            -2.4 + 0.15 * index, 0.35 - 0.03 * index,
            -0.07 + 0.01 * index, 0.008, 0.0, 0.0};
        surrogate.higher_yield.nodes.push_back(yield_node);

        fs::MultipleScatterRbfNode shape_node;
        shape_node.feature = {
            0.35 + 0.22 * index, 0.12 + 0.25 * index,
            0.18 + 0.17 * index};
        shape_node.log_conditional_yield = 0.0;
        for (std::size_t k = 0;
             k < fs::kMultipleScatterBasisSize; ++k) {
            shape_node.shape_coefficients[k] =
                0.2 + static_cast<double>(
                    (k + 3 * index)
                    % fs::kMultipleScatterBasisSize);
        }
        surrogate.higher_compton.nodes.push_back(shape_node);
    }
    auto& secondary = payload.source.secondary_surrogate;
    secondary.pair_511_atom.yield_bandwidth = {0.12, 0.12, 0.12};
    secondary.pair_511_atom.shape_bandwidth = {0.12, 0.12, 0.12};
    secondary.pair_annihilation_continuum.yield_bandwidth =
        {0.12, 0.12, 0.12};
    secondary.pair_annihilation_continuum.shape_bandwidth =
        {0.12, 0.12, 0.12};
    secondary.pair_annihilation_in_flight.yield_bandwidth =
        {0.12, 0.12, 0.12};
    secondary.pair_annihilation_in_flight.shape_bandwidth =
        {0.12, 0.12, 0.12};
    secondary.compton_electron_brems.yield_bandwidth =
        {0.12, 0.12, 0.12};
    secondary.compton_electron_brems.shape_bandwidth =
        {0.12, 0.12, 0.12};
    fs::MultipleScatterRbfNode secondary_node;
    secondary_node.feature = {0.67, 0.30, 0.28};
    secondary_node.log_conditional_yield = -0.2;
    for (std::size_t k = 0;
         k < fs::kMultipleScatterBasisSize; ++k) {
        secondary_node.shape_coefficients[k] =
            0.1 + static_cast<double>(
                (5 * k + 2) % fs::kMultipleScatterBasisSize);
    }
    secondary.pair_511_atom.nodes.push_back(secondary_node);
    secondary.pair_annihilation_continuum.nodes.push_back(
        secondary_node);
    secondary.pair_annihilation_in_flight.nodes.push_back(
        secondary_node);
    secondary.compton_electron_brems.nodes.push_back(
        secondary_node);
    payload.source.enable_source_electron_brems = true;
    payload.source.enable_photonuclear = true;
    fs::PhotonuclearCrossSectionRow photonuclear_lower;
    photonuclear_lower.energy_keV = 300.0;
    fs::PhotonuclearCrossSectionRow photonuclear_upper;
    photonuclear_upper.energy_keV = 1600.0;
    for (std::size_t index = 0;
         index < photonuclear_lower.mass_attenuation_cm2_g.size();
         ++index) {
        const double z = static_cast<double>(index + 1);
        photonuclear_lower.mass_attenuation_cm2_g[index] = 1e-7 * z;
        photonuclear_upper.mass_attenuation_cm2_g[index] = 2e-6 * z;
    }
    payload.source.photonuclear.cross_sections = {
        photonuclear_lower, photonuclear_upper};
    fs::FastSpectrumModel model(payload);
    const std::vector<fs::EmissionLine> lines{
        {356.013, 0.6205}, {661.657, 0.851}, {1173.228, 0.999}};
    const fs::ShieldParameters<Scalar> shield{p[1], p[2]};
    const fs::ResolutionParameters<Scalar> resolution{p[3], p[4], p[5]};
    const auto spectrum =
        model.predict(lines, p[0], shield, resolution, edges());
    Scalar answer = Scalar(0);
    // A nonlinear channel weighting makes this observable sensitive to
    // resolution. A linear spectral moment is almost invariant under a
    // symmetric Gaussian and leaves finite differences dominated by roundoff.
    for (std::size_t i = 0; i < spectrum.size(); ++i)
        answer += spectrum[i]
            * (1.0 + 0.2 * std::sin(0.37 * static_cast<double>(i)));
    return answer;
}

// Charged-source objective. The evaluation point is chosen
// so the range-ratio coordinate lands EXACTLY on an interior knot of the
// 8-basis f2 grid for one thickness in the sweep (the C2 continuity of the
// clamped cubic spline is what makes the interval-select Jet-safe), and a
// second call probes the frozen beyond-domain region.
template <typename Scalar>
Scalar charged_source_spline_objective(
    const std::array<Scalar, 2>& parameters) {
    fs::ChargedParticleSourcePayload payload;
    const auto make_field = [](std::uint32_t rank) {
        fs::ChargedParticleSplineField field;
        field.domain = {0.1, 1.05, 0.0, 1.0, 0.02, 0.9};
        field.yield_grid = {6, 5, 8};
        field.yield_coefficients.resize(6 * 5 * 8);
        for (std::size_t index = 0;
             index < field.yield_coefficients.size(); ++index) {
            field.yield_coefficients[index] =
                -0.4 + 0.013 * static_cast<double>((7 * index + 3) % 17);
        }
        field.shape_rank = rank;
        field.shape_basis.resize(rank);
        for (std::uint32_t r = 0; r < rank; ++r) {
            double total = 0.0;
            for (std::size_t k = 0;
                 k < fs::kMultipleScatterBasisSize; ++k) {
                field.shape_basis[r][k] =
                    1.0 + static_cast<double>((5 * k + 3 * r) % 11);
                total += field.shape_basis[r][k];
            }
            for (auto& value : field.shape_basis[r])
                value /= total;
        }
        if (rank >= 2) {
            field.weight_grid = {4, 4, 5};
            field.weight_log_odds.resize(rank - 1);
            for (std::size_t r = 0; r + 1 < rank; ++r) {
                field.weight_log_odds[r].resize(4 * 4 * 5);
                for (std::size_t index = 0;
                     index < field.weight_log_odds[r].size(); ++index) {
                    field.weight_log_odds[r][index] =
                        (static_cast<double>((5 * index + r) % 13) / 13.0)
                        - 0.5;
                }
            }
        }
        return field;
    };
    payload.positron_brems_spline = make_field(3);
    payload.positron_containment_spline = make_field(0);
    payload.positron_511_atom_spline = make_field(0);
    payload.positron_annihilation_continuum_spline = make_field(2);
    payload.positron_annihilation_in_flight_direct_spline = make_field(3);
    payload.positron_annihilation_in_flight_continuum_spline =
        make_field(2);
    const fs::ChargedParticleSourceFluxModel model(payload);
    Scalar answer = Scalar(0);
    // Middle thickness scaled so f2 sits on the interior knot u = 2/5 of the
    // 8-basis grid; the last thickness is deep in the frozen region.
    for (double thickness_scale : {0.35, 1.0, 40.0}) {
        const auto measure = model.evaluate(
            2290.0, parameters[0] * Scalar(thickness_scale),
            parameters[1], true, edges());
        for (std::size_t index = 0;
             index < measure.continuous.size(); ++index) {
            answer += measure.continuous[index].integral
                * (1.0 + 0.2 * std::sin(0.13 * index));
        }
        for (const auto& atom : measure.atoms)
            answer += Scalar(1.17) * atom.yield;
    }
    return answer;
}

template <typename Scalar>
Scalar detector_response_objective(const Scalar& cosine) {
    auto payload = fs::coaxial_hpge_payload();
    payload.use_response_surrogate = true;
    const auto populate =
        [](fs::DetectorResponseRbfComponent& component,
           double yield, double phase, bool pair_feature = false) {
            component.yield_bandwidth = {0.08, 0.12, 0.18};
            component.shape_bandwidth = {0.06, 0.10, 0.16};
            for (int node_index = 0; node_index < 3; ++node_index) {
                const double node_cosine = -0.6 + 0.6 * node_index;
                fs::DetectorResponseRbfNode node;
                const double energy_coordinate = pair_feature
                    ? std::log1p((2290.0 - 1021.9979) / 50.0)
                        / std::log1p(
                            (12000.0 - 1021.9979) / 50.0)
                    : std::log(2290.0 / 10.0)
                        / std::log(1200.0);
                node.feature = {
                    energy_coordinate,
                    std::acos(0.99999999 * node_cosine)
                        / 3.14159265358979323846,
                    1.0 / (1.0 + std::pow(2290.0 / 55.0, 3.0))};
                node.log_conditional_yield =
                    std::log(yield * (1.0 + 0.08 * node_index));
                for (std::size_t index = 0;
                     index < node.shape_coefficients.size(); ++index) {
                    node.shape_coefficients[index] =
                        1.0 + std::fmod(
                            phase + 7.0 * index + 3.0 * node_index,
                            23.0);
                }
                component.nodes.push_back(node);
            }
        };
    populate(payload.zero_deposit, 2.5, 1.0);
    populate(payload.full_energy, 0.25, 2.0);
    populate(payload.detector_photoelectric_escape, 0.035, 4.7);
    populate(payload.detector_single_compton, 0.25, 5.0);
    populate(payload.detector_multiple_compton, 0.38, 6.0);
    populate(payload.pair_fraction, 0.05, 7.0, true);
    populate(payload.pair_full_energy, 0.12, 7.5, true);
    populate(payload.single_escape, 0.18, 8.0, true);
    populate(payload.double_escape, 0.06, 8.5, true);
    populate(payload.annihilation_511, 0.01, 9.0, true);
    populate(payload.detector_pair_no_compton, 0.08, 9.5, true);
    populate(payload.detector_pair_single_compton, 0.18, 10.0, true);
    populate(payload.detector_pair_multiple_compton, 0.37, 10.5, true);
    // Continuum components carry their deposited-energy block partition; the
    // angle warp below is what this case is really exercising through Jet.
    payload.detector_single_compton.block_sizes = {64, 32, 0, 0};
    payload.detector_multiple_compton.block_sizes = {64, 32, 0, 0};
    payload.detector_pair_no_compton.block_sizes = {24, 24, 24, 24};
    payload.detector_pair_single_compton.block_sizes = {24, 24, 24, 24};
    payload.detector_pair_multiple_compton.block_sizes = {24, 24, 24, 24};
    payload.detected_rate_surface.enabled = true;
    payload.detected_rate_surface.angle_warp_exponent = 0.5;
    payload.detected_rate_surface.angle_warp_epsilon = 0.03;
    for (std::size_t angle = 0;
         angle < fs::kDetectorRateAngleBasisSize; ++angle) {
        for (std::size_t energy = 0;
             energy < fs::kDetectorRateEnergyBasisSize; ++energy) {
            const std::size_t index =
                angle * fs::kDetectorRateEnergyBasisSize + energy;
            payload.detected_rate_surface
                .log_odds_coefficients[index] =
                -0.7 + 0.04 * angle + 0.01 * energy;
        }
    }
    const fs::DetectorPulseResponseModel model(payload);
    const auto measure = model.evaluate(2290.0, cosine, edges());
    Scalar answer = Scalar(0);
    for (std::size_t index = 0;
         index < measure.continuous.size(); ++index) {
        answer += measure.continuous[index].integral
            * (1.0 + 0.17 * std::sin(0.11 * index));
    }
    for (const auto& atom : measure.atoms)
        answer += atom.yield * (1.0 + 0.0003 * atom.energy_keV);
    return answer;
}

/// A fitted payload named by an environment variable, or empty. Fitted
/// payloads are generated locally (tools/full_spectrum/README.md) and are not
/// part of the source tree, so the cases that need one skip without it.
std::string payload_text_from_environment(const char* variable) {
    const char* path = std::getenv(variable);
    if (path == nullptr || *path == '\0')
        return {};
    std::ifstream input(path);
    if (!input) {
        BOOST_TEST_MESSAGE(
            std::string(variable) + " names an unreadable file: " + path);
        return {};
    }
    BOOST_TEST_MESSAGE(std::string("using ") + variable + " = " + path);
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

template <typename Scalar>
Scalar fitted_payload_objective(
    const fs::FastSpectrumModel& model, const std::array<Scalar, 7>& p)
{
    const std::vector<fs::EmissionLine> lines{
        {661.657, 0.851}, {1173.228, 0.999}, {2614.511, 0.9976}};
    const fs::ShieldParameters<Scalar> shield{p[1], p[2], p[6]};
    const fs::ResolutionParameters<Scalar> resolution{p[3], p[4], p[5]};
    const auto spectrum =
        model.predict(lines, p[0], shield, resolution, edges());
    Scalar answer = Scalar(0);
    for (std::size_t i = 0; i < spectrum.size(); ++i)
        answer += spectrum[i]
            * (1.0 + 0.2 * std::sin(0.37 * static_cast<double>(i)));
    return answer;
}

} // namespace

BOOST_AUTO_TEST_CASE(jet_jacobian_matches_central_differences) {
    constexpr int N = 6;
    using Jet = ceres::Jet<double, N>;
    const std::array<double, N> value{
        2.3, 1.7, 27.35, 3.0, 1.2, 2e-4};
    std::array<Jet, N> jet;
    for (int i = 0; i < N; ++i) {
        jet[i] = Jet(value[i]);
        jet[i].v[i] = 1.0;
    }
    const Jet autodiff = objective(jet);

    for (int parameter = 0; parameter < N; ++parameter) {
        auto plus = value;
        auto minus = value;
        auto plus2 = value;
        auto minus2 = value;
        static constexpr double steps[N] = {
            1e-4, 1e-3, 1e-3, 1e-3, 1e-3, 1e-6};
        const double h = steps[parameter];
        plus[parameter] += h;
        minus[parameter] -= h;
        plus2[parameter] += 2.0 * h;
        minus2[parameter] -= 2.0 * h;
        const double finite_difference =
            (-objective(plus2) + 8.0 * objective(plus)
             - 8.0 * objective(minus) + objective(minus2)) / (12.0 * h);
        const double scale =
            std::max(1e-10, std::abs(finite_difference));
        BOOST_TEST_CONTEXT("parameter " << parameter) {
            BOOST_TEST(
                std::abs(autodiff.v[parameter] - finite_difference) / scale
                    < 1e-5);
        }
    }
}

BOOST_AUTO_TEST_CASE(fluorescence_bases_jet_matches_central_differences) {
    // Areal density and effective Z, at Z = 82.1 (bismuth inside the 0.5-2%
    // weight taper) and 82.5 (lead and bismuth both fully weighted).
    using Jet = ceres::Jet<double, 2>;
    const fs::SmoothElementCrossSection xs(0.42);
    const fs::SourceSecondarySurrogatePayload options;
    const auto bases = [&](const auto& areal_density, const auto& z) {
        using S = std::decay_t<decltype(z)>;
        const auto mu0 = xs.evaluate(661.657, z, S(0));
        const S mu0_total = mu0[0] + mu0[1] + mu0[2] + mu0[3];
        const auto b = fs::detail::source_component_bases(
            661.657, areal_density, z, S(0), mu0, mu0_total, xs, options);
        return b.base[fs::detail::kFluorescenceKAlpha]
            + 3.0 * b.base[fs::detail::kFluorescenceKBeta]
            + 10.0 * b.base[fs::detail::kFluorescenceL];
    };
    for (double z : {82.1, 82.5}) {
        const double t = 1.3;
        const Jet autodiff = bases(Jet(t, 0), Jet(z, 1));
        for (int parameter = 0; parameter < 2; ++parameter) {
            const double h = 1e-4;
            const auto at = [&](double shift) {
                return parameter == 0 ? bases(t + shift, z) : bases(t, z + shift);
            };
            const double finite_difference =
                (-at(2 * h) + 8.0 * at(h) - 8.0 * at(-h) + at(-2 * h)) / (12.0 * h);
            BOOST_TEST_CONTEXT("Z " << z << " parameter " << parameter) {
                BOOST_TEST(std::abs(autodiff.v[parameter] - finite_difference)
                               / std::max(1e-12, std::abs(finite_difference))
                           < 1e-5);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(charged_fluorescence_bases_jet_matches_central_differences) {
    // Areal density and effective Z of a positron's shield, thin and thick,
    // at Z = 82.1 (bismuth inside the weight taper) and 82.5.
    using Jet = ceres::Jet<double, 2>;
    const fs::SmoothElementCrossSection xs(0.42);
    const auto bases = [&](const auto& areal_density, const auto& z) {
        const auto b = fs::detail::charged_fluorescence_bases(
            1000.0, areal_density, z, std::decay_t<decltype(z)>(0), true, xs);
        return b.a[0] + 3.0 * b.a[1] + 10.0 * b.a[2] + b.b[0] + 3.0 * b.b[1]
            + 10.0 * b.b[2];
    };
    for (double t : {0.1, 2.0}) {
        for (double z : {82.1, 82.5}) {
            const Jet autodiff = bases(Jet(t, 0), Jet(z, 1));
            for (int parameter = 0; parameter < 2; ++parameter) {
                const double h = 1e-4;
                const auto at = [&](double shift) {
                    return parameter == 0 ? bases(t + shift, z) : bases(t, z + shift);
                };
                const double finite_difference =
                    (-at(2 * h) + 8.0 * at(h) - 8.0 * at(-h) + at(-2 * h)) / (12.0 * h);
                BOOST_TEST_CONTEXT("t " << t << " Z " << z << " parameter " << parameter) {
                    BOOST_TEST(std::abs(autodiff.v[parameter] - finite_difference)
                                   / std::max(1e-12, std::abs(finite_difference))
                               < 1e-5);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(
    charged_source_spline_jet_matches_central_differences) {
    constexpr int N = 2;
    using Jet = ceres::Jet<double, N>;
    const std::array<double, N> value{1.0, 35.25};
    std::array<Jet, N> jet;
    for (int index = 0; index < N; ++index) {
        jet[index] = Jet(value[index]);
        jet[index].v[index] = 1.0;
    }
    const Jet autodiff = charged_source_spline_objective(jet);
    const std::array<double, N> steps{1e-4, 1e-4};
    for (int parameter = 0; parameter < N; ++parameter) {
        auto plus = value;
        auto minus = value;
        auto plus2 = value;
        auto minus2 = value;
        plus[parameter] += steps[parameter];
        minus[parameter] -= steps[parameter];
        plus2[parameter] += 2.0 * steps[parameter];
        minus2[parameter] -= 2.0 * steps[parameter];
        const double finite_difference =
            (-charged_source_spline_objective(plus2)
             + 8.0 * charged_source_spline_objective(plus)
             - 8.0 * charged_source_spline_objective(minus)
             + charged_source_spline_objective(minus2))
            / (12.0 * steps[parameter]);
        const double scale =
            std::max(1e-10, std::abs(finite_difference));
        BOOST_TEST_CONTEXT("parameter " << parameter) {
            BOOST_TEST(
                std::abs(
                    autodiff.v[parameter] - finite_difference) / scale
                    < 1e-5);
        }
    }
}

BOOST_AUTO_TEST_CASE(
    fitted_charged_payload_jet_matches_central_differences) {
    const std::string text =
        payload_text_from_environment("CEELO_CHARGED_SOURCE_PAYLOAD");
    if (text.empty()) {
        BOOST_TEST_MESSAGE(
            "no fitted charged-particle payload; case skipped "
            "(set CEELO_CHARGED_SOURCE_PAYLOAD to run it)");
        return;
    }
    const auto model =
        fs::ChargedParticleSourceFluxModel::deserialize(text);
    const auto objective = [&](const auto& parameters) {
        using Scalar =
            typename std::decay_t<decltype(parameters[0])>;
        Scalar answer = Scalar(0);
        for (bool positron : {false, true}) {
            const auto measure = model.evaluate(
                1650.0, parameters[0], parameters[1], positron, edges());
            for (std::size_t index = 0;
                 index < measure.continuous.size(); ++index) {
                answer += measure.continuous[index].integral
                    * (1.0 + 0.2 * std::sin(0.13 * index));
            }
            for (const auto& atom : measure.atoms)
                answer += Scalar(1.17) * atom.yield;
        }
        return answer;
    };
    constexpr int N = 2;
    using Jet = ceres::Jet<double, N>;
    const std::array<double, N> value{1.7, 47.5};
    std::array<Jet, N> jet;
    for (int index = 0; index < N; ++index) {
        jet[index] = Jet(value[index]);
        jet[index].v[index] = 1.0;
    }
    const Jet autodiff = objective(jet);
    const std::array<double, N> steps{1e-4, 1e-4};
    for (int parameter = 0; parameter < N; ++parameter) {
        auto plus = value;
        auto minus = value;
        auto plus2 = value;
        auto minus2 = value;
        plus[parameter] += steps[parameter];
        minus[parameter] -= steps[parameter];
        plus2[parameter] += 2.0 * steps[parameter];
        minus2[parameter] -= 2.0 * steps[parameter];
        const double finite_difference =
            (-objective(plus2) + 8.0 * objective(plus)
             - 8.0 * objective(minus) + objective(minus2))
            / (12.0 * steps[parameter]);
        const double scale =
            std::max(1e-10, std::abs(finite_difference));
        BOOST_TEST_CONTEXT("parameter " << parameter) {
            BOOST_TEST(
                std::abs(
                    autodiff.v[parameter] - finite_difference) / scale
                    < 1e-5);
        }
    }
}

BOOST_AUTO_TEST_CASE(fitted_source_payload_jet_matches_central_differences) {
    const std::string text =
        payload_text_from_environment("CEELO_SOURCE_ESCAPE_PAYLOAD");
    if (text.empty()) {
        BOOST_TEST_MESSAGE(
            "no fitted source-escape payload; case skipped "
            "(set CEELO_SOURCE_ESCAPE_PAYLOAD to run it)");
        return;
    }
    fs::FullSpectrumPayload payload;
    payload.source =
        fs::SourceGeomFluxModel::deserialize(text).payload();
    payload.detector = fs::coaxial_hpge_payload();
    payload.source_quadrature_bins = 6;
    payload.detector.resolution_half_window_channels = 64;
    const fs::FastSpectrumModel model(payload);

    constexpr int N = 7;
    using Jet = ceres::Jet<double, N>;
    const std::array<double, N> value{
        2.3, 8.1, 27.35, 3.0, 1.2, 2e-4, 0.08};
    std::array<Jet, N> jet;
    for (int i = 0; i < N; ++i) {
        jet[i] = Jet(value[i]);
        jet[i].v[i] = 1.0;
    }
    const Jet autodiff = fitted_payload_objective(model, jet);

    for (int parameter = 0; parameter < N; ++parameter) {
        auto plus = value;
        auto minus = value;
        auto plus2 = value;
        auto minus2 = value;
        static constexpr double steps[N] = {
            1e-4, 1e-3, 1e-3, 1e-3, 1e-3, 1e-6, 1e-3};
        const double h = steps[parameter];
        plus[parameter] += h;
        minus[parameter] -= h;
        plus2[parameter] += 2.0 * h;
        minus2[parameter] -= 2.0 * h;
        const double finite_difference =
            (-fitted_payload_objective(model, plus2)
             + 8.0 * fitted_payload_objective(model, plus)
             - 8.0 * fitted_payload_objective(model, minus)
             + fitted_payload_objective(model, minus2)) / (12.0 * h);
        const double scale =
            std::max(1e-10, std::abs(finite_difference));
        BOOST_TEST_CONTEXT("parameter " << parameter) {
            BOOST_TEST(
                std::abs(autodiff.v[parameter] - finite_difference) / scale
                    < 1e-5);
        }
    }
}

BOOST_AUTO_TEST_CASE(detector_angle_jet_matches_central_difference) {
    using Jet = ceres::Jet<double, 1>;
    constexpr double value = 0.27;
    Jet cosine(value);
    cosine.v[0] = 1.0;
    const Jet autodiff = detector_response_objective(cosine);
    constexpr double h = 1e-4;
    const double finite_difference =
        (-detector_response_objective(value + 2.0 * h)
         + 8.0 * detector_response_objective(value + h)
         - 8.0 * detector_response_objective(value - h)
         + detector_response_objective(value - 2.0 * h))
        / (12.0 * h);
    BOOST_TEST(
        std::abs(autodiff.v[0] - finite_difference)
            / std::max(1e-10, std::abs(finite_difference))
            < 1e-5);
}
