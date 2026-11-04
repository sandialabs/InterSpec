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

#define BOOST_TEST_MODULE EnvFitPayloadTest
#include <boost/test/included/unit_test.hpp>

#include "io/EnvFitArchive.h"
#include "io/EnvFitPayload.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

using namespace ceelo;

namespace {

std::vector<double> cubic_bernstein_knots()
{
    return {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
}

std::vector<double> geometry_knots(std::uint32_t level)
{
    const std::size_t intervals = std::size_t{1} << (level - 1U);
    std::vector<double> knots(4U, 0.0);
    for (std::size_t i = 1; i < intervals; ++i)
        knots.push_back(static_cast<double>(i) / static_cast<double>(intervals));
    knots.insert(knots.end(), 4U, 1.0);
    return knots;
}

std::vector<double> greville_abscissae(const std::vector<double>& knots)
{
    const std::size_t n_basis = knots.size() - 4U;
    std::vector<double> result(n_basis, 0.0);
    for (std::size_t i = 0; i < n_basis; ++i)
        result[i] = (knots[i + 1U] + knots[i + 2U] + knots[i + 3U]) / 3.0;
    return result;
}

std::uint64_t process_id()
{
#if defined(_WIN32)
    return static_cast<std::uint64_t>(_getpid());
#else
    return static_cast<std::uint64_t>(getpid());
#endif
}

EnvFitRuntimeModel model_at(const EnvFitRuntimeQuery& query,
                            const std::string& name,
                            EnvFitAxisKind kind,
                            double yield,
                            std::size_t energy_basis_index,
                            std::size_t mu_basis_index)
{
    std::array<double, 4> q{};
    std::string error;
    if (!env_fit_transform_query(query, q, &error))
        throw std::runtime_error(error);
    EnvFitRuntimeModel model;
    model.name = name;
    model.axis_kind = kind;
    model.rank = 1;
    model.level = 1;
    model.energy_knots = kind == EnvFitAxisKind::LineE0 ||
                         kind == EnvFitAxisKind::Line511
        ? std::vector<double>{} : cubic_bernstein_knots();
    model.mu_knots = cubic_bernstein_knots();
    for (std::size_t i = 0; i < 4; ++i) {
        model.q_min[i] = q[i] - 0.25;
        model.q_span[i] = 0.5;
    }
    model.q_min[0] = std::log(59.5);
    model.q_span[0] = std::log(2614.0) - model.q_min[0];
    // A clamped B-spline tensor is a partition of unity.  Equal coefficients
    // therefore produce an exactly constant log-yield field.
    model.log_yield_coefficients.assign(4U * 4U * 4U * 4U, std::log(yield));
    model.logit_coefficients.assign(model.log_yield_coefficients.size(), 0.0);
    const std::size_t n_energy = model.energy_knots.empty() ? 1U : 4U;
    model.archetype_coefficients.assign(n_energy * 4U, 0.0);
    model.archetype_coefficients[energy_basis_index * 4U + mu_basis_index] = 1.0;
    return model;
}

void declare_model_inventory(EnvFitPayload& payload)
{
    payload.model_inventory.clear();
    for (const auto& model : payload.continuum_models)
        payload.model_inventory.push_back({model.name, model.axis_kind, false});
    for (const auto& model : payload.line_models)
        payload.model_inventory.push_back({model.name, model.axis_kind, true});
}

EnvFitPayload make_payload(const EnvFitRuntimeQuery& query)
{
    EnvFitPayload payload;
    payload.config_sha256 = env_fit_sha256_hex("payload-test-configuration");
    payload.material_name = "Air";
    payload.geometry_tier = "air";
    payload.accepted_for_publication = true;
    payload.brems_activation_threshold_keV = 150.0;
    payload.physical_domain_min = {{59.5, 0.5, 0.0, 0.0}};
    payload.physical_domain_max = {{2614.0, 49.87, 4.16, 4.16}};
    payload.validation_metadata_json =
        "{\"final_nodes\":64,\"max_joint_tv\":0.008}";
    payload.continuum_models.push_back(model_at(
        query, "primary_order1", EnvFitAxisKind::Primary, 2.0, 0, 0));
    EnvFitRuntimeModel inactive = model_at(
        query, "annihilation_order1", EnvFitAxisKind::Annihilation, 10.0, 3, 3);
    inactive.has_energy_threshold = true;
    inactive.energy_threshold_keV = 2.0 * kEnvFitElectronRestEnergyKeV;
    payload.continuum_models.push_back(inactive);
    payload.line_models.push_back(model_at(
        query, "elastic_E0", EnvFitAxisKind::LineE0, 0.5, 0, 0));
    EnvFitRuntimeModel line511 = model_at(
        query, "annihilation_511", EnvFitAxisKind::Line511, 4.0, 0, 3);
    line511.has_energy_threshold = true;
    line511.energy_threshold_keV = 2.0 * kEnvFitElectronRestEnergyKeV;
    payload.line_models.push_back(line511);
    declare_model_inventory(payload);
    return payload;
}

double sum(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

EnvFitPayload canonical_cross_language_fixture()
{
    EnvFitPayload payload;
    payload.config_sha256.assign(64, '0');
    payload.material_name = "Air";
    payload.geometry_tier = "air";
    payload.accepted_for_publication = true;
    payload.brems_activation_threshold_keV = 150.0;
    payload.physical_domain_min = {{1.0, 1.0, 0.0, 0.0}};
    payload.physical_domain_max = {{2.0, 2.0, 1.0, 1.0}};
    payload.validation_metadata_json = "{}";
    EnvFitRuntimeModel model;
    model.name = "primary_order1";
    model.axis_kind = EnvFitAxisKind::Primary;
    model.e_cut_keV = 10.0;
    model.rank = 1;
    model.level = 1;
    model.energy_knots = cubic_bernstein_knots();
    model.mu_knots = cubic_bernstein_knots();
    model.q_min = {{0.0, 0.0, 0.0, 0.0}};
    model.q_span = {{1.0, 1.0, 1.0, 1.0}};
    model.log_yield_coefficients.assign(256, 0.0);
    model.logit_coefficients.assign(256, 0.0);
    model.archetype_coefficients.assign(16, 0.0);
    model.archetype_coefficients.front() = 1.0;
    payload.continuum_models.push_back(std::move(model));
    declare_model_inventory(payload);
    return payload;
}

EnvFitRuntimeModel parity_model(const EnvFitRuntimeQuery& query, bool line)
{
    std::array<double, 4> q{};
    std::string error;
    if (!env_fit_transform_query(query, q, &error))
        throw std::runtime_error(error);
    EnvFitRuntimeModel model;
    model.name = line ? "elastic_E0" : "brems_order0";
    model.axis_kind = line ? EnvFitAxisKind::LineE0 : EnvFitAxisKind::Brems;
    model.rank = 2;
    model.level = 1;
    if (!line) {
        model.has_energy_threshold = true;
        model.energy_threshold_keV = 150.0;
        model.energy_knots = {0.0, 0.0, 0.0, 0.0, 0.4,
                              1.0, 1.0, 1.0, 1.0};
    }
    model.mu_knots = {0.0, 0.0, 0.0, 0.0, 0.3, 0.7,
                      1.0, 1.0, 1.0, 1.0};
    const std::array<double, 4> offsets{{0.23, 0.37, 0.11, 0.19}};
    model.q_span = {{1.1, 0.9, 0.6, 0.8}};
    for (std::size_t i = 0; i < 4; ++i) model.q_min[i] = q[i] - offsets[i];

    for (int i0 = 0; i0 < 4; ++i0) {
        for (int i1 = 0; i1 < 4; ++i1) {
            for (int i2 = 0; i2 < 4; ++i2) {
                for (int i3 = 0; i3 < 4; ++i3) {
                    const double base = 0.4 + 0.01 * (i0 - 1.5) +
                        0.02 * (i1 - 1.5) - 0.015 * (i2 - 1.5) +
                        0.005 * (i3 - 1.5);
                    model.log_yield_coefficients.push_back(base - (line ? 1.1 : 0.0));
                    model.logit_coefficients.push_back(
                        -0.3 + 0.05 * i0 - 0.04 * i1 + 0.03 * i2 - 0.02 * i3);
                    model.logit_coefficients.push_back(
                        0.2 - 0.02 * i0 + 0.06 * i1 - 0.01 * i2 + 0.04 * i3);
                }
            }
        }
    }
    if (line) {
        for (int r = 0; r < 2; ++r) {
            double normalizer = 0.0;
            for (int j = 0; j < 6; ++j)
                normalizer += r == 0 ? j + 1.0 : 6.0 - j;
            for (int j = 0; j < 6; ++j)
                model.archetype_coefficients.push_back(
                    (r == 0 ? j + 1.0 : 6.0 - j) / normalizer);
        }
    } else {
        for (int r = 0; r < 2; ++r) {
            std::vector<double> raw;
            double normalizer = 0.0;
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 6; ++j) {
                    const double value = r == 0
                        ? (i + 1.0) * (j + 2.0)
                        : (6.0 - i) * (7.0 - j) + 0.3 * i * j;
                    raw.push_back(value);
                    normalizer += value;
                }
            }
            for (double value : raw)
                model.archetype_coefficients.push_back(value / normalizer);
        }
    }
    return model;
}

}  // namespace

BOOST_AUTO_TEST_CASE(cubic_archetypes_are_integrated_not_sampled_as_histograms)
{
    const EnvFitRuntimeQuery query{122.0, 4.0, 1.0, 1.2, 0.25};
    const EnvFitPayload payload = make_payload(query);
    const double backscatter = query.e0_keV /
        (1.0 + 2.0 * query.e0_keV / kEnvFitElectronRestEnergyKeV);
    EnvFitEvaluation result;
    std::string error;
    EnvFitPayloadEvaluator evaluator;
    BOOST_REQUIRE_MESSAGE(evaluator.set_payload(payload, &error), error);
    BOOST_REQUIRE_MESSAGE(evaluator.evaluate(
        query, {0.0, 10.0, backscatter, query.e0_keV, 200.0},
        {-1.0, 0.0, 1.0}, result, &error), error);

    BOOST_REQUIRE_EQUAL(result.components.size(), 2U);
    BOOST_TEST(result.direct_fluence_per_cm2 == 0.25);
    // First normalized Bernstein basis is 4(1-u)^3.  Its exact mass on
    // u=[0,1/2] is 15/16, and the same identity holds for t=(mu+1)/2.
    const auto& primary = result.components[0].s_hat_cell_masses;
    BOOST_TEST(primary[0] == 0.0, boost::test_tools::tolerance(1e-14));
    BOOST_TEST(primary[2] == 2.0 * (15.0 / 16.0) * (15.0 / 16.0),
               boost::test_tools::tolerance(2e-13));
    BOOST_TEST(primary[3] == 2.0 * (15.0 / 16.0) * (1.0 / 16.0),
               boost::test_tools::tolerance(2e-13));
    BOOST_TEST(sum(primary) == 2.0, boost::test_tools::tolerance(2e-13));
    BOOST_TEST(sum(result.components[1].s_hat_cell_masses) == 0.0);
    BOOST_TEST(sum(result.continuum_s_hat_cell_masses) == 2.0,
               boost::test_tools::tolerance(2e-13));
    BOOST_TEST(sum(result.incident_continuum_cell_masses_per_cm2) == 0.5,
               boost::test_tools::tolerance(2e-13));
    BOOST_TEST(sum(result.components[0].incident_cell_masses_per_cm2) == 0.5,
               boost::test_tools::tolerance(2e-13));

    BOOST_REQUIRE_EQUAL(result.lines.size(), 1U);
    BOOST_TEST(result.lines[0].energy_keV == query.e0_keV);
    BOOST_TEST(result.lines[0].s_hat_mu_bin_masses[0] == 0.5 * 15.0 / 16.0,
               boost::test_tools::tolerance(2e-13));
    BOOST_TEST(result.lines[0].s_hat_mu_bin_masses[1] == 0.5 / 16.0,
               boost::test_tools::tolerance(2e-13));
    BOOST_TEST(sum(result.lines[0].incident_mu_bin_masses_per_cm2) == 0.125,
               boost::test_tools::tolerance(2e-13));
}

BOOST_AUTO_TEST_CASE(nested_geometry_levels_match_python_affine_reproduction)
{
    // The Python hierarchy test evaluates the same point, intercept, slopes,
    // and Greville coefficients.  Exact affine reproduction simultaneously
    // checks the dyadic knots, feature ordering/count, and C++ evaluation for
    // every serialized geometry level.
    const EnvFitRuntimeQuery query{662.0, 8.0, 1.5, 0.8, 0.2};
    std::array<double, 4> q{};
    std::string error;
    BOOST_REQUIRE_MESSAGE(env_fit_transform_query(query, q, &error), error);
    const std::array<double, 4> normalized{{0.13, 0.37, 0.61, 0.83}};
    const std::array<double, 4> slopes{{0.2, -0.3, 0.4, 0.1}};
    const double expected_log_yield = 0.7 +
        normalized[0] * slopes[0] + normalized[1] * slopes[1] +
        normalized[2] * slopes[2] + normalized[3] * slopes[3];

    for (std::uint32_t level = 1; level <= 4; ++level) {
        EnvFitPayload payload;
        payload.config_sha256 = env_fit_sha256_hex("nested-geometry-level");
        payload.material_name = "Air";
        payload.geometry_tier = "air";
        payload.accepted_for_publication = true;
        payload.brems_activation_threshold_keV = 150.0;
        payload.physical_domain_min = {{
            std::exp(q[0] - normalized[0]), 0.5, 0.0, 0.0}};
        payload.physical_domain_max = {{
            std::exp(q[0] - normalized[0] + 1.0), 49.87, 4.16, 4.16}};
        payload.validation_metadata_json = "{}";

        EnvFitRuntimeModel model;
        model.name = "primary_order1";
        model.axis_kind = EnvFitAxisKind::Primary;
        model.rank = 1;
        model.level = level;
        model.energy_knots = cubic_bernstein_knots();
        model.mu_knots = cubic_bernstein_knots();
        for (std::size_t dim = 0; dim < 4; ++dim) {
            model.q_min[dim] = q[dim] - normalized[dim];
            model.q_span[dim] = 1.0;
        }
        const auto greville = greville_abscissae(geometry_knots(level));
        for (double x0 : greville) {
            for (double x1 : greville) {
                for (double x2 : greville) {
                    for (double x3 : greville) {
                        model.log_yield_coefficients.push_back(
                            0.7 + slopes[0] * x0 + slopes[1] * x1 +
                            slopes[2] * x2 + slopes[3] * x3);
                        model.logit_coefficients.push_back(0.0);
                    }
                }
            }
        }
        model.archetype_coefficients.assign(16U, 0.0);
        model.archetype_coefficients.front() = 1.0;
        payload.continuum_models.push_back(model);
        declare_model_inventory(payload);

        EnvFitEvaluation result;
        BOOST_REQUIRE_MESSAGE(evaluate_env_fit_payload(
            payload, query, {0.0, 10.0, query.e0_keV}, {-1.0, 1.0},
            result, &error), error);
        BOOST_TEST(sum(result.continuum_s_hat_cell_masses) ==
                       std::exp(expected_log_yield),
                   boost::test_tools::tolerance(2.0e-12));

        if (level >= 3U) {
            // Non-nested dimensions (6^4 at level three, 7^4 at level four)
            // are rejected, so a coefficient block laid out for a non-nested
            // basis cannot be silently reinterpreted under the nested
            // hierarchy.
            EnvFitPayload invalid = payload;
            const std::size_t nonnested_basis =
                static_cast<std::size_t>(level) + 3U;
            const std::size_t nonnested_features = nonnested_basis
                * nonnested_basis * nonnested_basis * nonnested_basis;
            invalid.continuum_models.front().log_yield_coefficients.assign(
                nonnested_features, 0.0);
            invalid.continuum_models.front().logit_coefficients.assign(
                nonnested_features, 0.0);
            BOOST_TEST(!validate_env_fit_payload(invalid, &error));
            BOOST_TEST(error.find("rank/level") != std::string::npos);
        }
    }
}

BOOST_AUTO_TEST_CASE(binary_round_trip_is_deterministic_and_checksummed)
{
    const EnvFitRuntimeQuery query{122.0, 4.0, 1.0, 1.2, 0.25};
    const EnvFitPayload expected = make_payload(query);
    const std::string path = "/tmp/ceelo-fit-payload-" +
                             std::to_string(process_id()) + ".bin";
    const std::string second = path + ".second";
    std::string error;
    BOOST_REQUIRE_MESSAGE(save_env_fit_payload(expected, path, &error), error);
    EnvFitPayload restored;
    BOOST_REQUIRE_MESSAGE(load_env_fit_payload(restored, path, &error), error);
    BOOST_REQUIRE_MESSAGE(save_env_fit_payload(restored, second, &error), error);

    std::ifstream first_input(path, std::ios::binary);
    std::ifstream second_input(second, std::ios::binary);
    const std::string first_bytes((std::istreambuf_iterator<char>(first_input)), {});
    const std::string second_bytes((std::istreambuf_iterator<char>(second_input)), {});
    BOOST_TEST(first_bytes == second_bytes);
    BOOST_TEST(restored.validation_metadata_json == expected.validation_metadata_json);
    BOOST_TEST(restored.accepted_for_publication);
    BOOST_TEST(restored.brems_activation_threshold_keV == 150.0);
    BOOST_TEST(restored.model_content_sha256 ==
               env_fit_payload_model_content_sha256(restored));
    const auto restored_primary = std::find_if(
        restored.continuum_models.begin(), restored.continuum_models.end(),
        [](const auto& model) { return model.name == "primary_order1"; });
    const auto expected_primary = std::find_if(
        expected.continuum_models.begin(), expected.continuum_models.end(),
        [](const auto& model) { return model.name == "primary_order1"; });
    BOOST_REQUIRE(restored_primary != restored.continuum_models.end());
    BOOST_REQUIRE(expected_primary != expected.continuum_models.end());
    BOOST_TEST(restored_primary->archetype_coefficients ==
                   expected_primary->archetype_coefficients,
               boost::test_tools::per_element());

    // Damage one payload byte without updating the stored SHA-256.
    BOOST_REQUIRE_GT(first_bytes.size(), 80U);
    std::string corrupt = first_bytes;
    corrupt.back() = static_cast<char>(corrupt.back() ^ 0x01);
    {
        std::ofstream output(second, std::ios::binary | std::ios::trunc);
        output.write(corrupt.data(), static_cast<std::streamsize>(corrupt.size()));
    }
    EnvFitPayload rejected;
    BOOST_TEST(!load_env_fit_payload(rejected, second, &error));
    BOOST_TEST(error.find("checksum") != std::string::npos);
    std::remove(path.c_str());
    std::remove(second.c_str());
}

BOOST_AUTO_TEST_CASE(candidate_loading_is_explicit_and_content_hash_excludes_validation)
{
    const EnvFitRuntimeQuery query{122.0, 4.0, 1.0, 1.2, 0.25};
    EnvFitPayload candidate = make_payload(query);
    candidate.accepted_for_publication = false;
    const std::string model_hash = env_fit_payload_model_content_sha256(candidate);
    candidate.validation_metadata_json = "{\"acceptance_passed\":false}";
    BOOST_TEST(env_fit_payload_model_content_sha256(candidate) == model_hash);

    const std::string path = "/tmp/ceelo-fit-candidate-" +
                             std::to_string(process_id()) + ".bin";
    std::string error;
    BOOST_REQUIRE_MESSAGE(save_env_fit_payload(candidate, path, &error), error);
    EnvFitPayload restored;
    BOOST_TEST(!load_env_fit_payload(restored, path, &error));
    BOOST_TEST(error.find("candidate") != std::string::npos);
    BOOST_REQUIRE_MESSAGE(load_env_fit_payload(restored, path, &error, true), error);
    BOOST_TEST(!restored.accepted_for_publication);
    BOOST_TEST(restored.model_content_sha256 == model_hash);

    restored.continuum_models.front().name = "primary_typo";
    restored.model_inventory.front().name = "primary_typo";
    BOOST_TEST(!validate_env_fit_payload(restored, &error));
    const bool taxonomy_error = error.find("canonical") != std::string::npos ||
                                error.find("unknown") != std::string::npos;
    BOOST_TEST(taxonomy_error);
    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(binary_layout_matches_the_stdlib_python_fixture)
{
    // This digest was independently generated with Python struct.pack using
    // the documented little-endian field order.  It catches layout drift even
    // on installations where the optional NumPy fitter is not available.
    const std::string path = "/tmp/ceelo-fit-payload-fixture-" +
                             std::to_string(process_id()) + ".bin";
    std::string error;
    BOOST_REQUIRE_MESSAGE(save_env_fit_payload(
        canonical_cross_language_fixture(), path, &error), error);
    std::ifstream input(path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(input)), {});
    BOOST_TEST(env_fit_sha256_hex(bytes) ==
        "7a947f72ddc89c54337afe8274f390b7443f1cc83a23b504bb32f4824b399808");
    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(domain_support_and_normalization_are_enforced)
{
    const EnvFitRuntimeQuery query{122.0, 4.0, 1.0, 1.2, 0.25};
    EnvFitPayload payload = make_payload(query);
    std::string error;
    BOOST_TEST(validate_env_fit_payload(payload, &error), error);

    EnvFitEvaluation ignored;
    EnvFitRuntimeQuery outside = query;
    outside.detector_height_m += 2.0;
    BOOST_TEST(!evaluate_env_fit_payload(payload, outside,
                                         {0.0, 200.0}, {-1.0, 1.0},
                                         ignored, &error));
    BOOST_TEST(error.find("outside q domain") != std::string::npos);

    EnvFitRuntimeQuery outside_campaign = query;
    outside_campaign.e0_keV = 40.0;
    BOOST_TEST(!evaluate_env_fit_payload(payload, outside_campaign,
                                         {0.0, 200.0}, {-1.0, 1.0},
                                         ignored, &error));
    BOOST_TEST(error.find("physical domain") != std::string::npos);

    EnvFitRuntimeQuery missing_direct = query;
    missing_direct.direct_fluence_per_cm2 = 0.0;
    BOOST_TEST(!evaluate_env_fit_payload(payload, missing_direct,
                                         {0.0, 200.0}, {-1.0, 1.0},
                                         ignored, &error));
    BOOST_TEST(error.find("analytic direct fluence") != std::string::npos);

    payload.continuum_models[0].archetype_coefficients[0] = -1.0;
    BOOST_TEST(!validate_env_fit_payload(payload, &error));
    BOOST_TEST(error.find("nonnegative") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(softmax_weights_and_yields_remain_positive)
{
    const EnvFitRuntimeQuery query{662.0, 8.0, 1.5, 0.8, 0.1};
    EnvFitPayload payload;
    payload.config_sha256 = env_fit_sha256_hex("softmax-test");
    payload.material_name = "soil";
    payload.geometry_tier = "infinite-field";
    payload.brems_activation_threshold_keV = 150.0;
    payload.physical_domain_min = {{59.5, 0.5, 0.0, 0.0}};
    payload.physical_domain_max = {{2614.0, 49.87, 4.16, 4.16}};
    payload.validation_metadata_json = "{}";
    EnvFitRuntimeModel model = model_at(
        query, "brems_order0", EnvFitAxisKind::Brems, 3.0, 0, 0);
    model.has_energy_threshold = true;
    model.energy_threshold_keV = 150.0;
    model.rank = 2;
    model.logit_coefficients.resize(model.log_yield_coefficients.size() * 2U);
    for (std::size_t f = 0; f < model.log_yield_coefficients.size(); ++f) {
        model.logit_coefficients[2U * f] = -1000.0;
        model.logit_coefficients[2U * f + 1U] = 1000.0;
    }
    model.archetype_coefficients.assign(2U * 4U * 4U, 0.0);
    model.archetype_coefficients[0] = 1.0;
    model.archetype_coefficients[16U + 15U] = 1.0;
    payload.continuum_models.push_back(model);
    declare_model_inventory(payload);

    EnvFitEvaluation result;
    std::string error;
    BOOST_REQUIRE_MESSAGE(evaluate_env_fit_payload(
        payload, query, {0.0, 10.0, 662.0, 1000.0}, {-1.0, 0.0, 1.0},
        result, &error), error);
    BOOST_TEST(sum(result.continuum_s_hat_cell_masses) == 3.0,
               boost::test_tools::tolerance(2e-13));
    BOOST_TEST(sum(result.incident_continuum_cell_masses_per_cm2) == 0.3,
               boost::test_tools::tolerance(2e-13));
    for (double value : result.continuum_s_hat_cell_masses) BOOST_TEST(value >= 0.0);
    model.energy_threshold_keV = 149.0;
    payload.continuum_models.front() = model;
    BOOST_TEST(!validate_env_fit_payload(payload, &error));
    BOOST_TEST(error.find("brems support") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(nontrivial_python_cpp_evaluator_golden_vector)
{
    // The same rank-2 fixture and expected vector are evaluated by the Python
    // RuntimePayload test.  It exercises all four q coordinates, nonconstant
    // logits/yields, interior spline knots, and non-knot requested edges.
    const EnvFitRuntimeQuery query{662.0, 8.0, 1.5, 0.8, 0.2};
    EnvFitPayload payload;
    payload.config_sha256 = env_fit_sha256_hex("python-cpp-parity-v1");
    payload.material_name = "Air";
    payload.geometry_tier = "air";
    payload.brems_activation_threshold_keV = 150.0;
    payload.physical_domain_min = {{550.0, 0.5, 0.0, 0.0}};
    payload.physical_domain_max = {{1500.0, 49.87, 4.16, 4.16}};
    payload.validation_metadata_json = "{\"parity_fixture\":1}";
    payload.continuum_models.push_back(parity_model(query, false));
    payload.line_models.push_back(parity_model(query, true));
    declare_model_inventory(payload);

    const std::vector<double> expected_continuum{
        0.0, 0.0, 0.0, 0.0,
        0.08332336753366282, 0.08478522302174381,
        0.046526312153207856, 0.07838387164298673,
        0.12131648531579858, 0.13280593285509656,
        0.0796843340401528, 0.15671859710585198,
        0.06857324135921093, 0.08270895357362763,
        0.05479029268231385, 0.12335331946729351,
        0.04509909119202966, 0.05932341261672018,
        0.0423212377411697, 0.10354784535377856,
        0.020488893946752323, 0.028309544066437688,
        0.020960043780111244, 0.05322364672284197,
        0.0, 0.0, 0.0, 0.0
    };
    const std::vector<double> expected_line{
        0.12409460886151082, 0.1348171647338879,
        0.08019496002549564, 0.15562079951936106
    };
    EnvFitEvaluation result;
    std::string error;
    BOOST_REQUIRE_MESSAGE(evaluate_env_fit_payload(
        payload, query, {0.0, 10.0, 17.0, 83.0, 251.0, 500.0, 662.0, 900.0},
        {-1.0, -0.73, -0.1, 0.37, 1.0}, result, &error), error);
    BOOST_REQUIRE_EQUAL(result.continuum_s_hat_cell_masses.size(), expected_continuum.size());
    for (std::size_t i = 0; i < expected_continuum.size(); ++i)
        BOOST_TEST(result.continuum_s_hat_cell_masses[i] == expected_continuum[i],
                   boost::test_tools::tolerance(1.0e-11));
    for (std::size_t i = 0; i < expected_continuum.size(); ++i)
        BOOST_TEST(result.incident_continuum_cell_masses_per_cm2[i] ==
                       0.2 * expected_continuum[i],
                   boost::test_tools::tolerance(1.0e-11));
    BOOST_REQUIRE_EQUAL(result.lines.size(), 1U);
    BOOST_REQUIRE_EQUAL(result.lines[0].s_hat_mu_bin_masses.size(), expected_line.size());
    for (std::size_t i = 0; i < expected_line.size(); ++i)
        BOOST_TEST(result.lines[0].s_hat_mu_bin_masses[i] == expected_line[i],
                   boost::test_tools::tolerance(1.0e-11));
    for (std::size_t i = 0; i < expected_line.size(); ++i)
        BOOST_TEST(result.lines[0].incident_mu_bin_masses_per_cm2[i] ==
                       0.2 * expected_line[i],
                   boost::test_tools::tolerance(1.0e-11));
}
