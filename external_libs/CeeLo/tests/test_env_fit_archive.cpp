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

#define BOOST_TEST_MODULE EnvFitArchiveTest
#include <boost/test/unit_test.hpp>

#include "io/EnvFitArchive.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>

using namespace ceelo;

namespace {

EnvFitMomentArray constant_moments(std::size_t cells, std::uint64_t histories,
                                   double scale = 1.0)
{
    EnvFitMomentArray result;
    result.n.assign(cells, histories);
    result.s1.resize(cells);
    result.s2.resize(cells);
    for (std::size_t i = 0; i < cells; ++i) {
        const double contribution = scale * static_cast<double>(i + 1) * 1e-6;
        result.s1[i] = static_cast<double>(histories) * contribution;
        result.s2[i] = static_cast<double>(histories) * contribution * contribution;
    }
    return result;
}

std::vector<float> replica_sums(std::size_t cells,
                                std::uint64_t histories_per_replica,
                                double scale = 1.0)
{
    std::vector<float> result(kEnvFitReplicaCount * cells);
    for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
        for (std::size_t i = 0; i < cells; ++i) {
            const double contribution = scale * static_cast<double>(i + 1) * 1e-6;
            result[replica * cells + i] = static_cast<float>(
                static_cast<double>(histories_per_replica) * contribution);
        }
    }
    return result;
}

EnvFitBlock make_block(std::uint64_t block_id, std::uint64_t history_base,
                       std::uint64_t seed_base)
{
    EnvFitBlock block;
    block.block_id = block_id;

    block.manifest.canonical_config_json =
        "{\"energy_bins\":384,\"material\":\"air\",\"mu_bins\":64}";
    block.manifest.config_sha256 =
        env_fit_sha256_hex(block.manifest.canonical_config_json);
    block.manifest.material_name = "air";
    block.manifest.material_definition_json = "{\"density\":0.001205}";
    block.manifest.geometry_tier = "air";
    block.manifest.rng_algorithm = "mt19937_64";
    block.manifest.estimator_flags = "nee,rayleigh,brems";
    block.manifest.units = "keV,cm,photon/cm2/source-photon";
    block.manifest.observable = "incident scalar fluence";
    block.manifest.normalization = "raw per source photon";

    block.axes.energy_axes = {
        {"primary", 1, {0.0, 0.5, 1.0}, {10.0, 184.3, 662.0}},
        {"annihilation", 1, {0.0, 0.5, 1.0}, {10.0, 170.3, 511.0}}
    };
    block.axes.total_energy_axis = 0;
    block.axes.mu_edges = {-1.0, 0.0, 1.0};
    block.axes.lines = {
        {"elastic_e0", 662.0, 0, 0},
        {"annihilation_511", 511.0, 1, 0}
    };
    block.components = {
        {"primary_one", 0, 1, 0, {1, 1}},
        {"primary_two_plus", 0, 2, 0, {1, 1}},
        {"annihilation_one", 1, 1, 1, {1, 1}}
    };

    constexpr std::uint64_t histories_per_replica =
        kEnvFitFormalBlockHistories / kEnvFitReplicaCount;
    block.node.e0_keV = 662.0;
    block.node.distance_cm = 1000.0;
    block.node.h_src_cm = 100.0;
    block.node.h_det_cm = 100.0;
    block.node.transformed_coordinates = {
        std::log(662.0), std::log(10.0), std::atan(0.1), std::atan(0.1)
    };
    block.node.role = "train";
    block.node.status = "converged";
    block.node.direct_transmission = 0.8;
    block.node.direct_fluence = 6.4e-8;
    block.node.histories = histories_per_replica * kEnvFitReplicaCount;
    block.node.wall_seconds = 3.0;
    block.node.cpu_seconds = 10.0;
    block.node.target_total_rel_se = 0.0025;
    block.node.achieved_total_rel_se = 0.002;
    block.node_id = env_fit_node_id(block.manifest.config_sha256, block.node);

    for (std::size_t i = 0; i < kEnvFitReplicaCount; ++i) {
        block.replicas[i].replica_id = static_cast<std::uint32_t>(i);
        block.replicas[i].history_begin =
            history_base + i * histories_per_replica;
        block.replicas[i].history_count = histories_per_replica;
        block.replicas[i].rng_seed = seed_base + i;
    }

    const std::size_t component_cells =
        env_fit_continuum_component_cells(block);
    const std::size_t total_cells = env_fit_continuum_total_cells(block);
    const std::size_t line_component_cells =
        env_fit_line_component_cells(block);
    const std::size_t line_total_cells = env_fit_line_total_cells(block);
    block.continuum_components =
        constant_moments(component_cells, block.node.histories, 1.0);
    block.continuum_total =
        constant_moments(total_cells, block.node.histories, 2.0);
    block.line_components =
        constant_moments(line_component_cells, block.node.histories, 3.0);
    block.line_total.n.assign(line_total_cells, block.node.histories);
    block.line_total.s1.assign(line_total_cells, 0.0);
    block.line_total.s2.assign(line_total_cells, 0.0);
    for (std::size_t cell = 0; cell < line_total_cells; ++cell) {
        double contribution = 0.0;
        for (std::size_t component = 0; component < block.components.size(); ++component) {
            contribution += block.line_components.s1[
                component * line_total_cells + cell] / block.node.histories;
        }
        block.line_total.s1[cell] = block.node.histories * contribution;
        block.line_total.s2[cell] =
            block.node.histories * contribution * contribution;
    }
    block.replica_sums.continuum_components =
        replica_sums(component_cells, histories_per_replica, 1.0);
    block.replica_sums.continuum_total =
        replica_sums(total_cells, histories_per_replica, 2.0);
    block.replica_sums.line_components =
        replica_sums(line_component_cells, histories_per_replica, 3.0);
    block.replica_sums.line_total.assign(
        kEnvFitReplicaCount * line_total_cells, 0.0f);
    for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
        for (std::size_t cell = 0; cell < line_total_cells; ++cell) {
            float total = 0.0f;
            for (std::size_t component = 0;
                 component < block.components.size(); ++component) {
                total += block.replica_sums.line_components[
                    replica * line_component_cells +
                    component * line_total_cells + cell];
            }
            block.replica_sums.line_total[
                replica * line_total_cells + cell] = total;
        }
    }

    block.stopping.definition_json =
        env_fit_stopping_definition_json(block.axes);
    block.stopping.n = block.node.histories;
    for (std::size_t i = 0; i < kEnvFitStoppingDimension; ++i) {
        const double xi = static_cast<double>(i + 1) * 1e-5;
        block.stopping.s1[i] = block.node.histories * xi;
        for (std::size_t j = 0; j < kEnvFitStoppingDimension; ++j) {
            const double xj = static_cast<double>(j + 1) * 1e-5;
            block.stopping.outer_s2[i * kEnvFitStoppingDimension + j] =
                block.node.histories * xi * xj;
        }
    }

    EnvFitDiagnosticChannel air;
    air.name = "air_only";
    air.continuum = constant_moments(total_cells, block.node.histories, 0.5);
    air.lines = constant_moments(line_total_cells, block.node.histories, 0.25);
    block.diagnostics.path_channels.push_back(std::move(air));
    block.diagnostics.interaction_names = {"compton", "rayleigh", "pair"};
    block.diagnostics.interaction_counts = {100, 50, 2};
    block.diagnostics.scattered_total.n = block.node.histories;
    const double per_history_total = 78e-5; // Sum of stopping bands 1..12.
    block.diagnostics.scattered_total.s1 =
        block.node.histories * per_history_total;
    block.diagnostics.scattered_total.s2 =
        block.node.histories * per_history_total * per_history_total;
    block.diagnostics.dropped_weight.n = block.node.histories;
    block.diagnostics.dropped_weight.s1 = block.node.histories * 1e-12;
    block.diagnostics.dropped_weight.s2 = block.node.histories * 1e-24;
    block.diagnostics.maximum_history_contribution = 2e-4;
    return block;
}

} // namespace

BOOST_AUTO_TEST_CASE(sha256_and_stable_ids_are_standard_and_deterministic)
{
    BOOST_CHECK_EQUAL(
        env_fit_sha256_hex(""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    BOOST_CHECK_EQUAL(
        env_fit_sha256_hex("abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    BOOST_CHECK_EQUAL(env_fit_stable_id("abc"), 0xba7816bf8f01cfeaULL);
    BOOST_CHECK_EQUAL(env_fit_stable_id("abc"), env_fit_stable_id("abc"));
    BOOST_CHECK_NE(env_fit_stable_id("abc"), env_fit_stable_id("abd"));
    const EnvFitBlock block = make_block(1, 0, 1);
    BOOST_CHECK_EQUAL(block.node_id, 0xe00a5bb3240138d0ULL);
}

BOOST_AUTO_TEST_CASE(worker_block_round_trips_all_scientific_records)
{
    const EnvFitBlock source = make_block(17, 0, 1000);
    std::string error;
    BOOST_REQUIRE_MESSAGE(validate_env_fit_block(source, &error), error);

    const std::string path = "/tmp/ceelo_env_fit_roundtrip.ceefit";
    BOOST_REQUIRE_MESSAGE(save_env_fit_block_atomic(source, path, &error), error);

    EnvFitBlock loaded;
    BOOST_REQUIRE_MESSAGE(load_env_fit_block(loaded, path, &error), error);
    BOOST_CHECK_EQUAL(loaded.node_id, source.node_id);
    BOOST_CHECK_EQUAL(loaded.block_id, source.block_id);
    BOOST_CHECK_EQUAL(loaded.manifest.config_sha256,
                      source.manifest.config_sha256);
    BOOST_CHECK_EQUAL(loaded.manifest.canonical_config_json,
                      source.manifest.canonical_config_json);
    BOOST_CHECK_EQUAL(loaded.manifest.normalization,
                      source.manifest.normalization);
    BOOST_CHECK_EQUAL(loaded.axes.energy_axes[1].physical_edges_keV[1], 170.3);
    BOOST_CHECK_EQUAL(loaded.components[2].post_creation_order, 1U);
    BOOST_CHECK_EQUAL(loaded.replicas[31].rng_seed, 1031U);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        loaded.continuum_components.s1.begin(),
        loaded.continuum_components.s1.end(),
        source.continuum_components.s1.begin(),
        source.continuum_components.s1.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
        loaded.replica_sums.line_components.begin(),
        loaded.replica_sums.line_components.end(),
        source.replica_sums.line_components.begin(),
        source.replica_sums.line_components.end());
    BOOST_CHECK_EQUAL(loaded.stopping.outer_s2[15 * 16 + 14],
                      source.stopping.outer_s2[15 * 16 + 14]);
    BOOST_CHECK_EQUAL(loaded.stopping.definition_json,
                      source.stopping.definition_json);
    BOOST_CHECK_EQUAL(loaded.diagnostics.path_channels.front().name, "air_only");
    BOOST_CHECK_EQUAL(loaded.diagnostics.interaction_counts[1], 50U);
    BOOST_CHECK_EQUAL(loaded.diagnostics.scattered_total.s1,
                      source.diagnostics.scattered_total.s1);

    EnvFitBlock invalid = source;
    invalid.block_id = 0;
    BOOST_CHECK(!save_env_fit_block_atomic(invalid, path, &error));
    EnvFitBlock preserved;
    BOOST_REQUIRE(load_env_fit_block(preserved, path, &error));
    BOOST_CHECK_EQUAL(preserved.block_id, source.block_id);
    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(checksum_and_truncation_corruption_are_rejected)
{
    const EnvFitBlock source = make_block(18, 0, 2000);
    const std::string path = "/tmp/ceelo_env_fit_corrupt.ceefit";
    std::string error;
    BOOST_REQUIRE(save_env_fit_block_atomic(source, path, &error));

    {
        std::fstream file(path, std::ios::in | std::ios::out | std::ios::binary);
        BOOST_REQUIRE(file);
        file.seekg(100);
        char byte = 0;
        file.read(&byte, 1);
        byte ^= 0x40;
        file.seekp(100);
        file.write(&byte, 1);
    }
    EnvFitBlock sentinel;
    sentinel.node_id = 999;
    BOOST_CHECK(!load_env_fit_block(sentinel, path, &error));
    BOOST_CHECK_EQUAL(sentinel.node_id, 999U); // Failure does not modify output.
    BOOST_CHECK_NE(error.find("checksum"), std::string::npos);

    {
        std::ofstream file(path, std::ios::binary | std::ios::trunc);
        file.write("CEEFIT", 6);
    }
    BOOST_CHECK(!load_env_fit_block(sentinel, path, &error));
    BOOST_CHECK_NE(error.find("truncated"), std::string::npos);
    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(validation_rejects_inconsistent_scientific_metadata)
{
    std::string error;
    EnvFitBlock bad = make_block(19, 0, 3000);
    bad.manifest.config_sha256[0] = '0';
    BOOST_CHECK(!validate_env_fit_block(bad, &error));

    bad = make_block(19, 0, 3000);
    bad.replicas[1].rng_seed = bad.replicas[0].rng_seed;
    BOOST_CHECK(!validate_env_fit_block(bad, &error));

    bad = make_block(19, 0, 3000);
    bad.continuum_total.n[0] -= 1;
    BOOST_CHECK(!validate_env_fit_block(bad, &error));

    bad = make_block(19, 0, 3000);
    bad.components[0].support_mask.pop_back();
    BOOST_CHECK(!validate_env_fit_block(bad, &error));

    bad = make_block(19, 0, 3000);
    bad.replica_sums.line_total[0] += 1.0f;
    BOOST_CHECK(!validate_env_fit_block(bad, &error));

    bad = make_block(19, 0, 3000);
    bad.diagnostics.scattered_total.s1 *= 1.01;
    bad.diagnostics.scattered_total.s2 *= 1.01 * 1.01;
    BOOST_CHECK(!validate_env_fit_block(bad, &error));

    bad = make_block(19, 0, 3000);
    bad.manifest.observable.clear();
    BOOST_CHECK(!validate_env_fit_block(bad, &error));

    bad = make_block(19, 0, 3000);
    bad.stopping.definition_json = "{}";
    BOOST_CHECK(!validate_env_fit_block(bad, &error));
}

BOOST_AUTO_TEST_CASE(reducer_merges_disjoint_blocks_and_rejects_replays)
{
    const EnvFitBlock first = make_block(20, 0, 4000);
    const EnvFitBlock second =
        make_block(21, kEnvFitFormalBlockHistories, 5000);
    std::string error;
    EnvFitBlockReducer reducer;
    BOOST_REQUIRE_MESSAGE(reducer.add(first, &error), error);
    BOOST_REQUIRE_MESSAGE(reducer.add(second, &error), error);
    BOOST_CHECK_EQUAL(reducer.block_ids().size(), 2U);
    BOOST_CHECK_EQUAL(reducer.history_ranges().size(),
                      2U * kEnvFitReplicaCount);
    BOOST_CHECK_EQUAL(reducer.history_ranges().back().block_id, second.block_id);
    BOOST_CHECK_EQUAL(reducer.history_ranges().back().replica.replica_id, 31U);
    BOOST_CHECK_EQUAL(reducer.aggregate().block_id, 0U);
    BOOST_CHECK_EQUAL(reducer.aggregate().node.histories,
                      2U * first.node.histories);
    BOOST_CHECK_CLOSE(reducer.aggregate().continuum_total.s1[2],
                      2.0 * first.continuum_total.s1[2], 1e-10);
    BOOST_CHECK_CLOSE(reducer.aggregate().replica_sums.line_components[5],
                      2.0f * first.replica_sums.line_components[5], 1e-4);

    BOOST_CHECK(!reducer.add(second, &error));
    BOOST_CHECK_NE(error.find("duplicate block_id"), std::string::npos);

    EnvFitBlock overlap = make_block(22, 100, 6000);
    BOOST_CHECK(!reducer.add(overlap, &error));
    BOOST_CHECK_NE(error.find("overlaps"), std::string::npos);

    EnvFitBlock mismatch =
        make_block(23, 2 * kEnvFitFormalBlockHistories, 7000);
    mismatch.manifest.canonical_config_json = "{\"different\":true}";
    mismatch.manifest.config_sha256 =
        env_fit_sha256_hex(mismatch.manifest.canonical_config_json);
    mismatch.node_id =
        env_fit_node_id(mismatch.manifest.config_sha256, mismatch.node);
    BOOST_REQUIRE(validate_env_fit_block(mismatch, &error));
    BOOST_CHECK(!reducer.add(mismatch, &error));
    BOOST_CHECK_NE(error.find("manifest"), std::string::npos);

    BOOST_CHECK_EQUAL(reducer.block_ids().size(), 2U);
    BOOST_CHECK_EQUAL(reducer.aggregate().node.histories,
                      2U * first.node.histories);
}

BOOST_AUTO_TEST_CASE(mean_and_standard_error_use_raw_per_history_moments)
{
    const EnvFitEstimate estimate = env_fit_estimate(4, 10.0, 30.0);
    BOOST_CHECK_CLOSE(estimate.mean, 2.5, 1e-12);
    BOOST_CHECK_CLOSE(estimate.standard_error, std::sqrt(5.0 / 12.0), 1e-12);

    EnvFitMomentArray moments;
    moments.n = {4};
    moments.s1 = {10.0};
    moments.s2 = {30.0};
    BOOST_CHECK_CLOSE(env_fit_estimate(moments, 0).standard_error,
                      estimate.standard_error, 1e-12);
    BOOST_CHECK(std::isnan(env_fit_estimate(1, 2.0, 4.0).standard_error));
    BOOST_CHECK(std::isnan(env_fit_estimate(moments, 1).mean));
}
