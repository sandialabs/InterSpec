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

#include "io/EnvFitBlockBuilder.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstring>
#include <ctime>
#include <deque>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace ceelo {

namespace {

/// Process-local executor shared by successive blocks with the same requested
/// concurrency.  A fit campaign commonly appends many blocks, so constructing
/// and joining 32 replica threads for each 262k-history block is pure overhead.
/// Tasks own no pool state and every block waits for all of its tasks, making
/// cross-block scheduling irrelevant to RNG identity or reduction order.
class PersistentWorkerPool {
public:
    explicit PersistentWorkerPool(int workers)
    {
        threads_.reserve(static_cast<std::size_t>(workers));
        for (int i = 0; i < workers; ++i) {
            threads_.emplace_back([this] {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(mutex_);
                        ready_.wait(lock, [this] {
                            return stopping_ || !tasks_.empty();
                        });
                        if (stopping_ && tasks_.empty()) return;
                        task = std::move(tasks_.front());
                        tasks_.pop_front();
                    }
                    task();
                }
            });
        }
    }

    ~PersistentWorkerPool()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        ready_.notify_all();
        for (std::thread& thread : threads_) thread.join();
    }

    PersistentWorkerPool(const PersistentWorkerPool&) = delete;
    PersistentWorkerPool& operator=(const PersistentWorkerPool&) = delete;

    void submit(std::function<void()> task)
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            tasks_.push_back(std::move(task));
        }
        ready_.notify_one();
    }

private:
    std::mutex mutex_;
    std::condition_variable ready_;
    std::deque<std::function<void()>> tasks_;
    std::vector<std::thread> threads_;
    bool stopping_ = false;
};

PersistentWorkerPool& fit_worker_pool(int workers)
{
    static std::mutex pools_mutex;
    static std::map<int, std::unique_ptr<PersistentWorkerPool>> pools;
    std::lock_guard<std::mutex> lock(pools_mutex);
    std::unique_ptr<PersistentWorkerPool>& pool = pools[workers];
    if (!pool) pool.reset(new PersistentWorkerPool(workers));
    return *pool;
}

void append_u32_le(std::string& bytes, std::uint32_t value)
{
    for (unsigned i = 0; i < 4; ++i)
        bytes.push_back(static_cast<char>((value >> (8 * i)) & 0xffU));
}

void append_u64_le(std::string& bytes, std::uint64_t value)
{
    for (unsigned i = 0; i < 8; ++i)
        bytes.push_back(static_cast<char>((value >> (8 * i)) & 0xffU));
}

void append_f64_le(std::string& bytes, double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "double must be binary64");
    std::memcpy(&bits, &value, sizeof(bits));
    append_u64_le(bytes, bits);
}

std::uint64_t checked_history_begin(std::uint64_t block_id,
                                    std::uint64_t block_histories,
                                    std::uint32_t replica_id,
                                    std::uint64_t replica_histories)
{
    const std::uint64_t block_index = block_id - 1;
    if (block_index > std::numeric_limits<std::uint64_t>::max() /
                          block_histories)
        throw std::invalid_argument("fit block history range overflows");
    const std::uint64_t base = block_index * block_histories;
    const std::uint64_t offset =
        static_cast<std::uint64_t>(replica_id) * replica_histories;
    if (offset > std::numeric_limits<std::uint64_t>::max() - base)
        throw std::invalid_argument("fit replica history range overflows");
    return base + offset;
}

std::vector<double> physical_edges(const FluenceAxes& axes)
{
    if (!axes.edges.empty()) return axes.edges;
    std::vector<double> result(static_cast<std::size_t>(axes.n_energy) + 1);
    for (int i = 0; i < axes.n_energy; ++i)
        result[static_cast<std::size_t>(i)] = axes.energy_low(i);
    result.back() = axes.energy_high(axes.n_energy - 1);
    return result;
}

std::vector<double> registered_edges(int bins)
{
    std::vector<double> result(static_cast<std::size_t>(bins) + 1);
    for (int i = 0; i <= bins; ++i)
        result[static_cast<std::size_t>(i)] =
            static_cast<double>(i) / static_cast<double>(bins);
    return result;
}

std::vector<double> mu_edges(int bins)
{
    std::vector<double> result(static_cast<std::size_t>(bins) + 1);
    for (int i = 0; i <= bins; ++i)
        result[static_cast<std::size_t>(i)] =
            -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(bins);
    result.front() = -1.0;
    result.back() = 1.0;
    return result;
}

EnvFitMomentArray moments(const std::vector<double>& s1,
                          const std::vector<double>& s2,
                          std::uint64_t histories)
{
    if (s1.size() != s2.size())
        throw std::logic_error("fit tally S1/S2 dimension mismatch");
    EnvFitMomentArray result;
    result.n.assign(s1.size(), histories);
    result.s1 = s1;
    result.s2 = s2;
    return result;
}

const char* origin_name(FluxOrigin origin)
{
    switch (origin) {
    case FluxOrigin::Primary: return "primary";
    case FluxOrigin::Annihilation: return "annihilation";
    case FluxOrigin::Brems: return "brems";
    default: return "unknown";
    }
}

const char* order_name(PostCreationOrder order)
{
    switch (order) {
    case PostCreationOrder::Zero: return "order0";
    case PostCreationOrder::One: return "order1";
    case PostCreationOrder::TwoPlus: return "order2plus";
    default: return "unknown";
    }
}

bool component_has_support(FluxOrigin origin, PostCreationOrder order,
                           double e0_keV, double brems_floor_keV)
{
    if (origin == FluxOrigin::Primary)
        return order != PostCreationOrder::Zero;
    if (origin == FluxOrigin::Annihilation)
        return e0_keV >= 2.0 * 510.998950 &&
               order != PostCreationOrder::Zero;
    return e0_keV >= brems_floor_keV;
}

} // namespace

std::uint64_t env_fit_node_id(const std::string& config_sha256,
                              double e0_keV, double distance_cm,
                              double h_src_cm, double h_det_cm,
                              double transverse_extent_cm)
{
    std::string key("CEEFIT-NODE-v1\0", 15);
    key += config_sha256;
    append_f64_le(key, e0_keV);
    append_f64_le(key, distance_cm);
    append_f64_le(key, h_src_cm);
    append_f64_le(key, h_det_cm);
    append_f64_le(key, transverse_extent_cm);
    return env_fit_stable_id(key);
}

std::uint64_t env_fit_replica_seed(const std::string& config_sha256,
                                   std::uint64_t node_id,
                                   std::uint64_t block_id,
                                   std::uint32_t replica_id)
{
    std::string key("CEEFIT-RNG-v1\0", 14);
    key += config_sha256;
    append_u64_le(key, node_id);
    append_u64_le(key, block_id);
    append_u32_le(key, replica_id);
    return env_fit_stable_id(key);
}

EnvFitBlock run_env_fit_block(const EnvironmentGeometry& environment,
                              const EnvScatterConfig& config,
                              const Eigen::Vector3d& source,
                              const EnvFitBlockRequest& request)
{
    if (request.block_id == 0 || request.histories == 0 ||
        request.histories % kEnvFitReplicaCount != 0)
        throw std::invalid_argument(
            "fit block requires nonzero id and a history count divisible by 32");
    if (request.history_slot_size < request.histories ||
        request.history_slot_size % kEnvFitReplicaCount != 0)
        throw std::invalid_argument(
            "fit block history slot must cover histories and be divisible by 32");
    if (!(request.e0_keV > config.axes.e_min_keV) ||
        !(request.distance_cm > 0.0) || request.threads < 1)
        throw std::invalid_argument("invalid fit block node or thread count");
    if (std::abs(request.h_src_cm - request.h_det_cm) >= request.distance_cm)
        throw std::invalid_argument("invalid fit geometry: |h_src-h_det| >= d");
    if (request.manifest.config_sha256 !=
        env_fit_sha256_hex(request.manifest.canonical_config_json))
        throw std::invalid_argument("fit manifest configuration hash mismatch");
    std::string manifest_error;
    if (!validate_env_fit_manifest(request.manifest, &manifest_error))
        throw std::invalid_argument("invalid fit manifest: " + manifest_error);

    const auto component_axes = make_fit_component_axes(
        request.e0_keV, config.axes.n_energy, config.axes.n_mu,
        config.axes.e_min_keV);
    if (!config.axes.same_axis(component_axes[static_cast<int>(FluxOrigin::Primary)]))
        throw std::invalid_argument(
            "fit block common axis must be the registered primary axis");

    const std::uint64_t node_id = env_fit_node_id(
        request.manifest.config_sha256, request.e0_keV, request.distance_cm,
        request.h_src_cm, request.h_det_cm, request.transverse_extent_cm);
    const std::uint64_t histories_per_replica =
        request.histories / kEnvFitReplicaCount;

    std::vector<EnvScatterResult> replica_results;
    replica_results.reserve(kEnvFitReplicaCount);
    for (std::size_t i = 0; i < kEnvFitReplicaCount; ++i)
        replica_results.emplace_back(config.axes, request.e0_keV, component_axes);

    const NextEventEstimator estimator(environment, config, source);
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    const int worker_count = std::min<int>(
        request.threads, static_cast<int>(kEnvFitReplicaCount));
    const auto wall_start = std::chrono::steady_clock::now();
    const std::clock_t cpu_start = std::clock();
    std::size_t remaining = kEnvFitReplicaCount;
    std::mutex completion_mutex;
    std::condition_variable complete;
    PersistentWorkerPool& pool = fit_worker_pool(worker_count);
    std::size_t submitted = 0;
    std::exception_ptr submit_error;
    try {
        for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
            pool.submit([&, replica] {
                try {
                    const std::uint64_t seed = env_fit_replica_seed(
                        request.manifest.config_sha256, node_id,
                        request.block_id, static_cast<std::uint32_t>(replica));
                    estimator.run(request.e0_keV, histories_per_replica, seed,
                                  replica_results[replica],
                                  static_cast<int>(replica));
                } catch (...) {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!worker_error)
                        worker_error = std::current_exception();
                }
                {
                    std::lock_guard<std::mutex> lock(completion_mutex);
                    if (--remaining == 0) complete.notify_one();
                }
            });
            ++submitted;
        }
    } catch (...) {
        submit_error = std::current_exception();
        const std::size_t unscheduled = kEnvFitReplicaCount - submitted;
        {
            std::lock_guard<std::mutex> lock(completion_mutex);
            remaining -= unscheduled;
            if (remaining == 0) complete.notify_one();
        }
    }
    {
        std::unique_lock<std::mutex> lock(completion_mutex);
        complete.wait(lock, [&remaining] { return remaining == 0; });
    }
    {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (worker_error) std::rethrow_exception(worker_error);
    }
    if (submit_error) std::rethrow_exception(submit_error);

    const double wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    const double cpu_seconds = static_cast<double>(std::clock() - cpu_start) /
                               static_cast<double>(CLOCKS_PER_SEC);

    EnvScatterResult merged(config.axes, request.e0_keV, component_axes);
    for (const auto& result : replica_results) {
        merged.tally.merge(result.tally);
        merged.diagnostics.merge(result.diagnostics);
        merged.direct_transmission = result.direct_transmission;
        merged.direct_fluence = result.direct_fluence;
        merged.source_distance_cm = result.source_distance_cm;
        merged.histories += result.histories;
    }

    EnvFitBlock block;
    block.node_id = node_id;
    block.block_id = request.block_id;
    block.manifest = request.manifest;
    block.node.e0_keV = request.e0_keV;
    block.node.distance_cm = request.distance_cm;
    block.node.h_src_cm = request.h_src_cm;
    block.node.h_det_cm = request.h_det_cm;
    block.node.transverse_extent_cm = request.transverse_extent_cm;
    const double d_m = request.distance_cm * 0.01;
    const double h_src_m = request.h_src_cm * 0.01;
    const double h_det_m = request.h_det_cm * 0.01;
    const double dz_m = h_src_m - h_det_m;
    const double rho_m = std::sqrt(std::max(
        std::numeric_limits<double>::min(), d_m * d_m - dz_m * dz_m));
    block.node.transformed_coordinates = {{
        std::log(request.e0_keV), std::log(rho_m),
        std::atan2(h_src_m, rho_m), std::atan2(h_det_m, rho_m)}};
    block.node.role = request.role;
    block.node.status = "complete";
    block.node.direct_transmission = merged.direct_transmission;
    block.node.direct_fluence = merged.direct_fluence;
    block.node.histories = merged.tally.num_histories();
    block.node.wall_seconds = wall_seconds;
    block.node.cpu_seconds = cpu_seconds;
    block.node.target_total_rel_se = request.target_total_rel_se;

    const std::array<const char*, 3> axis_names{{
        "primary_registered", "annihilation_511_registered",
        "brems_log_relative"}};
    for (int origin = 0; origin < static_cast<int>(FluxOrigin::Count); ++origin) {
        EnvFitEnergyAxis axis;
        axis.name = axis_names[static_cast<std::size_t>(origin)];
        axis.transform_version = 1;
        axis.registered_edges = registered_edges(component_axes[origin].n_energy);
        axis.physical_edges_keV = physical_edges(component_axes[origin]);
        block.axes.energy_axes.push_back(std::move(axis));
    }
    block.axes.total_energy_axis = 0;
    block.axes.mu_edges = mu_edges(config.axes.n_mu);
    for (int line = 0; line < merged.tally.n_lines(); ++line) {
        const LineDef& definition = merged.tally.line(line);
        EnvFitLineDefinition archive_line;
        archive_line.name = definition.origin == LineOrigin::ElasticE0
                                ? "elastic_E0" : "annihilation_511";
        archive_line.energy_keV = definition.energy_keV;
        archive_line.origin = static_cast<std::uint32_t>(
            definition.origin == LineOrigin::ElasticE0
                ? FluxOrigin::Primary : FluxOrigin::Annihilation);
        // A line may contain several all-Rayleigh interaction orders.
        archive_line.post_creation_order =
            static_cast<std::uint32_t>(PostCreationOrder::Count);
        block.axes.lines.push_back(std::move(archive_line));
    }

    for (int component = 0; component < kFluxComponentCount; ++component) {
        const FluxOrigin origin = flux_component_origin(component);
        const PostCreationOrder order = flux_component_order(component);
        EnvFitComponent descriptor;
        descriptor.name = std::string(origin_name(origin)) + "_" + order_name(order);
        descriptor.origin = static_cast<std::uint32_t>(origin);
        descriptor.post_creation_order = static_cast<std::uint32_t>(order);
        descriptor.energy_axis = static_cast<std::uint32_t>(origin);
        const bool supported = component_has_support(
            origin, order, request.e0_keV, config.brems_electron_floor_keV);
        descriptor.support_mask.assign(
            static_cast<std::size_t>(component_axes[static_cast<int>(origin)].n_energy),
            supported ? 1U : 0U);
        block.components.push_back(std::move(descriptor));
    }

    const FluenceTally& tally = merged.tally;
    block.continuum_components = moments(
        tally.component_sum_w(), tally.component_sum_w2(), request.histories);
    block.continuum_total = moments(
        tally.continuum_sum_w(), tally.continuum_sum_w2(), request.histories);
    block.line_components = moments(
        tally.line_component_sum_w(), tally.line_component_sum_w2(),
        request.histories);
    block.line_total = moments(tally.line_sum_w(), tally.line_sum_w2(),
                               request.histories);

    const std::size_t component_cells = tally.component_sum_w().size();
    const std::size_t total_cells = tally.continuum_sum_w().size();
    const std::size_t line_component_cells =
        tally.line_component_sum_w().size();
    const std::size_t line_total_cells = tally.line_sum_w().size();
    block.replica_sums.continuum_components.resize(
        kEnvFitReplicaCount * component_cells);
    block.replica_sums.continuum_total.resize(kEnvFitReplicaCount * total_cells);
    block.replica_sums.line_components.resize(
        kEnvFitReplicaCount * line_component_cells);
    block.replica_sums.line_total.resize(
        kEnvFitReplicaCount * line_total_cells);
    for (std::size_t replica = 0; replica < kEnvFitReplicaCount; ++replica) {
        const FluenceTally& source_tally = replica_results[replica].tally;
        std::transform(source_tally.component_sum_w().begin(),
                       source_tally.component_sum_w().end(),
                       block.replica_sums.continuum_components.begin() +
                           static_cast<std::ptrdiff_t>(replica * component_cells),
                       [](double value) { return static_cast<float>(value); });
        std::transform(source_tally.continuum_sum_w().begin(),
                       source_tally.continuum_sum_w().end(),
                       block.replica_sums.continuum_total.begin() +
                           static_cast<std::ptrdiff_t>(replica * total_cells),
                       [](double value) { return static_cast<float>(value); });
        std::transform(source_tally.line_component_sum_w().begin(),
                       source_tally.line_component_sum_w().end(),
                       block.replica_sums.line_components.begin() +
                           static_cast<std::ptrdiff_t>(
                               replica * line_component_cells),
                       [](double value) { return static_cast<float>(value); });
        std::transform(source_tally.line_sum_w().begin(),
                       source_tally.line_sum_w().end(),
                       block.replica_sums.line_total.begin() +
                           static_cast<std::ptrdiff_t>(
                               replica * line_total_cells),
                       [](double value) { return static_cast<float>(value); });
        EnvFitReplicaRecord& record = block.replicas[replica];
        record.replica_id = static_cast<std::uint32_t>(replica);
        record.history_count = histories_per_replica;
        record.history_begin = checked_history_begin(
            request.block_id, request.history_slot_size,
            static_cast<std::uint32_t>(replica), histories_per_replica);
        record.rng_seed = env_fit_replica_seed(
            request.manifest.config_sha256, node_id, request.block_id,
            static_cast<std::uint32_t>(replica));
    }

    block.stopping.definition_json =
        env_fit_stopping_definition_json(block.axes);
    block.stopping.n = request.histories;
    block.stopping.s1 = tally.coarse_sum_w();
    block.stopping.outer_s2 = tally.coarse_outer_sum_w();

    const int n_classes = static_cast<int>(ChainClass::Count);
    const std::array<const char*, 4> class_names{{
        "air_only", "one_slab_bounce", "multi_same_slab",
        "multi_distinct_slab"}};
    for (int cls = 0; cls < n_classes; ++cls) {
        EnvFitDiagnosticChannel channel;
        channel.name = class_names[static_cast<std::size_t>(cls)];
        std::vector<double> continuum_s1(total_cells), continuum_s2(total_cells);
        for (std::size_t cell = 0; cell < total_cells; ++cell) {
            continuum_s1[cell] = tally.class_sum_w()[cell * n_classes + cls];
            continuum_s2[cell] = tally.class_sum_w2()[cell * n_classes + cls];
        }
        std::vector<double> line_s1(line_total_cells), line_s2(line_total_cells);
        for (std::size_t cell = 0; cell < line_total_cells; ++cell) {
            line_s1[cell] = tally.line_class_sum_w()[cell * n_classes + cls];
            line_s2[cell] = tally.line_class_sum_w2()[cell * n_classes + cls];
        }
        channel.continuum = moments(continuum_s1, continuum_s2, request.histories);
        channel.lines = moments(line_s1, line_s2, request.histories);
        block.diagnostics.path_channels.push_back(std::move(channel));
    }
    block.diagnostics.dropped_weight.n = request.histories;
    block.diagnostics.dropped_weight.s1 = tally.dropped_weight();
    block.diagnostics.dropped_weight.s2 = tally.dropped_weight_sum_w2();
    block.diagnostics.scattered_total = {
        request.histories, tally.total_sum_w(), tally.total_sum_w2()};
    block.diagnostics.maximum_history_contribution = tally.max_history_total();
    block.diagnostics.interaction_names = {
        "photoelectric", "compton", "rayleigh", "pair_production",
        "brems_emission"};
    block.diagnostics.interaction_counts.assign(
        merged.diagnostics.interaction_counts.begin(),
        merged.diagnostics.interaction_counts.end());
    block.diagnostics.transport_failures = merged.diagnostics.failures;
    block.diagnostics.capped =
        merged.diagnostics.max_vertex_terminations != 0;
    if (block.diagnostics.transport_failures != 0) {
        block.node.status = "failed";
        block.diagnostics.failure_reason =
            std::to_string(block.diagnostics.transport_failures) +
            " recoverable transport failures";
        if (block.diagnostics.capped) {
            block.diagnostics.failure_reason += "; max_vertices also terminated " +
                std::to_string(merged.diagnostics.max_vertex_terminations) +
                " photon chains";
        }
    } else if (block.diagnostics.capped) {
        block.node.status = "capped";
        block.diagnostics.failure_reason =
            "max_vertices terminated " +
            std::to_string(merged.diagnostics.max_vertex_terminations) +
            " photon chains";
    }

    const EnvFitEstimate total = env_fit_estimate(
        request.histories, tally.total_sum_w(), tally.total_sum_w2());
    block.node.achieved_total_rel_se =
        total.mean > 0.0 && std::isfinite(total.standard_error)
            ? total.standard_error / total.mean : 0.0;

    std::string validation_error;
    if (!validate_env_fit_block(block, &validation_error))
        throw std::runtime_error("constructed invalid fit block: " + validation_error);
    return block;
}

std::vector<EnvFitBlock> run_env_fit_blocks(
    const EnvironmentGeometry& environment,
    const EnvScatterConfig& config,
    const Eigen::Vector3d& source,
    const std::vector<EnvFitBlockRequest>& requests)
{
    std::vector<EnvFitBlock> blocks;
    blocks.reserve(requests.size());
    for (const EnvFitBlockRequest& request : requests)
        blocks.push_back(
            run_env_fit_block(environment, config, source, request));
    return blocks;
}

} // namespace ceelo
