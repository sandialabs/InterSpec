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

#define BOOST_TEST_MODULE EnvFitBlockBuilderTest
#include <boost/test/unit_test.hpp>

#include "io/EnvFitBlockBuilder.h"
#include "materials/Material.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

using namespace ceelo;

namespace {

EnvFitManifest manifest()
{
    EnvFitManifest result;
    result.canonical_config_json = "{\"schema\":\"builder-test-v1\"}";
    result.config_sha256 = env_fit_sha256_hex(result.canonical_config_json);
    result.material_name = "Air";
    result.material_definition_json = "{}";
    result.geometry_tier = "air";
    result.rng_algorithm = "mt19937_64/sha256-tuple-stream-v1";
    result.estimator_flags = "next-event";
    result.units = "energy=keV,length=cm,fluence=cm^-2";
    result.observable = "bin-integrated scalar incident fluence";
    result.normalization = "raw";
    return result;
}

EnvFitBlock build(int threads, std::uint64_t block_id = 1,
                  int max_vertices = 50,
                  std::uint64_t histories = 32)
{
    const Material air = make_Air();
    const EnvironmentGeometry environment(air);
    EnvScatterConfig config;
    config.e_floor_keV = 10.0;
    config.axes = make_registered_landmark_axes(122.0, 12, 4, 10.0);
    config.scoring_radius_cm = 4.0;
    config.max_vertices = max_vertices;

    EnvFitBlockRequest request;
    request.manifest = manifest();
    request.block_id = block_id;
    request.role = "pilot";
    request.e0_keV = 122.0;
    request.distance_cm = 100.0;
    request.histories = histories;
    request.threads = threads;
    return run_env_fit_block(environment, config,
                             Eigen::Vector3d(0.0, 0.0, -100.0), request);
}

double sum(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

} // namespace

BOOST_AUTO_TEST_CASE(node_and_replica_identity_match_the_cross_language_spec)
{
    const std::string config_sha256(64, '0');
    const std::uint64_t node = env_fit_node_id(
        config_sha256, 122.0, 100.0, 0.0, 0.0, 0.0);
    BOOST_TEST(node == 0xc5dbb8b7b2b1ad9eULL);
    BOOST_TEST(env_fit_replica_seed(config_sha256, node, 1, 0) ==
               0x8fa2f48da2f95952ULL);
}

BOOST_AUTO_TEST_CASE(scheduling_does_not_change_fit_block_statistics)
{
    const EnvFitBlock serial = build(1);
    const EnvFitBlock parallel = build(4);

    BOOST_TEST(serial.node_id == parallel.node_id);
    BOOST_TEST(serial.continuum_total.s1 == parallel.continuum_total.s1,
               boost::test_tools::per_element());
    BOOST_TEST(serial.continuum_total.s2 == parallel.continuum_total.s2,
               boost::test_tools::per_element());
    BOOST_TEST(serial.continuum_components.s1 ==
                   parallel.continuum_components.s1,
               boost::test_tools::per_element());
    BOOST_TEST(serial.line_components.s1 == parallel.line_components.s1,
               boost::test_tools::per_element());
    BOOST_TEST(serial.stopping.s1 == parallel.stopping.s1,
               boost::test_tools::per_element());
    BOOST_TEST(serial.stopping.definition_json ==
               env_fit_stopping_definition_json(serial.axes));
    BOOST_TEST(serial.diagnostics.interaction_names ==
                   parallel.diagnostics.interaction_names,
               boost::test_tools::per_element());
    BOOST_TEST(serial.diagnostics.interaction_counts ==
                   parallel.diagnostics.interaction_counts,
               boost::test_tools::per_element());
    for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
        BOOST_TEST(serial.replicas[replica].rng_seed ==
                   parallel.replicas[replica].rng_seed);
        BOOST_TEST(serial.replicas[replica].history_count == 1U);
    }
}

BOOST_AUTO_TEST_CASE(formal_depth_statistics_are_thread_count_invariant)
{
    const EnvFitBlock serial = build(1, 17, 50, kEnvFitFormalBlockHistories);
    const EnvFitBlock parallel = build(8, 17, 50, kEnvFitFormalBlockHistories);

    BOOST_TEST(serial.continuum_total.s1 == parallel.continuum_total.s1,
               boost::test_tools::per_element());
    BOOST_TEST(serial.continuum_total.s2 == parallel.continuum_total.s2,
               boost::test_tools::per_element());
    BOOST_TEST(serial.continuum_components.s1 ==
                   parallel.continuum_components.s1,
               boost::test_tools::per_element());
    BOOST_TEST(serial.line_components.s1 == parallel.line_components.s1,
               boost::test_tools::per_element());
    for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica)
        BOOST_TEST(serial.replicas[replica].history_count == 8192U);
}

BOOST_AUTO_TEST_CASE(tiny_blocks_reuse_persistent_pool_without_completion_lifetime_race)
{
    for (std::uint64_t block_id = 100; block_id < 132; ++block_id) {
        const EnvFitBlock block = build(32, block_id);
        BOOST_TEST(block.node.histories == 32U);
    }
}

BOOST_AUTO_TEST_CASE(writer_rejects_diagnostic_history_count)
{
    const EnvFitBlock diagnostic = build(1);
    std::string error;
    BOOST_REQUIRE_MESSAGE(validate_env_fit_block(diagnostic, &error), error);
    BOOST_TEST(!save_env_fit_block_atomic(
        diagnostic, "/tmp/ceelo_env_fit_noncanonical.ceefit", &error));
    BOOST_TEST(error.find("exactly 262144 histories") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(transport_diagnostics_are_named_and_caps_are_not_silent)
{
    const EnvFitBlock normal = build(2);
    const std::vector<std::string> expected{
        "photoelectric", "compton", "rayleigh", "pair_production",
        "brems_emission"};
    BOOST_TEST(normal.diagnostics.interaction_names == expected,
               boost::test_tools::per_element());
    BOOST_TEST(normal.diagnostics.interaction_counts.size() == expected.size());
    BOOST_TEST(std::accumulate(normal.diagnostics.interaction_counts.begin(),
                               normal.diagnostics.interaction_counts.end(),
                               std::uint64_t{0}) > 0U);

    // max_vertices=0 deliberately truncates each primary before transport.
    // The worker block remains inspectable, but cannot masquerade as complete.
    const EnvFitBlock capped = build(1, 1, 0);
    BOOST_TEST(capped.diagnostics.capped);
    BOOST_TEST(capped.node.status == "capped");
    BOOST_TEST(capped.diagnostics.failure_reason.find("32 photon chains") !=
               std::string::npos);
    std::string error;
    BOOST_TEST(validate_env_fit_block(capped, &error), error);
}

BOOST_AUTO_TEST_CASE(native_components_and_exact_lines_close_to_common_total)
{
    const EnvFitBlock block = build(2);
    BOOST_REQUIRE_EQUAL(block.axes.energy_axes.size(), 3U);
    BOOST_REQUIRE_EQUAL(block.components.size(),
                        static_cast<std::size_t>(kFluxComponentCount));

    const double native_total = sum(block.continuum_components.s1) +
                                sum(block.line_components.s1);
    const double common_total = sum(block.continuum_total.s1) +
                                sum(block.line_total.s1);
    BOOST_TEST(native_total == common_total,
               boost::test_tools::tolerance(1e-12));

    std::string error;
    BOOST_TEST(validate_env_fit_block(block, &error), error);
}

BOOST_AUTO_TEST_CASE(block_ranges_are_disjoint_and_ids_are_stable)
{
    const Material air = make_Air();
    const EnvironmentGeometry environment(air);
    EnvScatterConfig config;
    config.e_floor_keV = 10.0;
    config.axes = make_registered_landmark_axes(122.0, 12, 4, 10.0);
    std::vector<EnvFitBlockRequest> requests(2);
    for (std::size_t i = 0; i < requests.size(); ++i) {
        requests[i].manifest = manifest();
        requests[i].block_id = i + 1;
        requests[i].role = "pilot";
        requests[i].e0_keV = 122.0;
        requests[i].distance_cm = 100.0;
        requests[i].histories = 32;
        requests[i].threads = 1;
    }
    const std::vector<EnvFitBlock> blocks = run_env_fit_blocks(
        environment, config, Eigen::Vector3d(0.0, 0.0, -100.0), requests);
    BOOST_REQUIRE_EQUAL(blocks.size(), 2U);
    const EnvFitBlock& first = blocks[0];
    const EnvFitBlock& second = blocks[1];
    BOOST_TEST(first.node_id == second.node_id);
    BOOST_TEST(first.block_id != second.block_id);
    for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
        const auto& a = first.replicas[replica];
        const auto& b = second.replicas[replica];
        BOOST_TEST(a.history_begin + a.history_count <= b.history_begin);
        BOOST_TEST(a.rng_seed != b.rng_seed);
    }

    EnvFitBlockReducer reducer;
    std::string error;
    BOOST_REQUIRE(reducer.add(first, &error));
    BOOST_TEST(reducer.add(second, &error));
    BOOST_TEST(reducer.aggregate().node.histories == 64U);
}
