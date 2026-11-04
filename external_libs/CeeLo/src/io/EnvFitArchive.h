#pragma once
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

/// @file EnvFitArchive.h
/// @brief Lossless, mergeable Monte-Carlo worker blocks for environment fits.
///
/// `EnvFitBlock` is the scientific record produced by one independently seeded
/// block of histories.  It is deliberately separate from EnvKernelTable, which
/// is a compact runtime product and is not rich enough to refit or continue a
/// simulation.  The binary codec has no third-party dependency.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ceelo {

constexpr std::uint32_t kEnvFitBlockFormatVersion = 1;
constexpr std::size_t kEnvFitReplicaCount = 32;
constexpr std::size_t kEnvFitStoppingDimension = 16;
constexpr std::uint64_t kEnvFitFormalBlockHistories = 262144;

/// Configuration and observable identity that must be identical for blocks
/// combined into one fit node.  `canonical_config_json` is supplied by the
/// campaign driver in canonical (stable key order/number formatting) form.  The
/// codec verifies that `config_sha256` is the SHA-256 of those exact UTF-8
/// bytes.
struct EnvFitManifest {
    std::string canonical_config_json;
    std::string config_sha256;

    std::string material_name;
    std::string material_definition_json;
    std::string geometry_tier;
    std::string rng_algorithm;
    std::string estimator_flags;
    std::string units;
    std::string observable;
    std::string normalization;
};

/// One registered/physical energy coordinate.  Different physical components
/// may use different coordinates (primary Compton landmark, fixed 511-keV
/// landmark, logarithmic brems coordinate) in the same block.
struct EnvFitEnergyAxis {
    std::string name;
    std::uint32_t transform_version = 1;
    std::vector<double> registered_edges;
    std::vector<double> physical_edges_keV;
};

/// Generic component descriptor.  Numeric tags intentionally do not hard-code
/// the current number of origins; the transport layer owns those enums.
struct EnvFitComponent {
    std::string name;
    std::uint32_t origin = 0;
    std::uint32_t post_creation_order = 0;
    std::uint32_t energy_axis = 0;
    /// One byte per energy bin.  Zero denotes structural absence.
    std::vector<std::uint8_t> support_mask;
};

struct EnvFitLineDefinition {
    std::string name;
    double energy_keV = 0.0;
    std::uint32_t origin = 0;
    std::uint32_t post_creation_order = 0;
};

struct EnvFitAxes {
    std::vector<EnvFitEnergyAxis> energy_axes;
    /// Axis on which the independently accumulated total continuum is stored.
    std::uint32_t total_energy_axis = 0;
    std::vector<double> mu_edges;
    std::vector<EnvFitLineDefinition> lines;
};

struct EnvFitNodeMetadata {
    double e0_keV = 0.0;
    double distance_cm = 0.0;
    double h_src_cm = 0.0;
    double h_det_cm = 0.0;
    double transverse_extent_cm = 0.0;
    /// (log(E0/keV), log(rho/m), atan(h_src/rho), atan(h_det/rho)).
    std::array<double, 4> transformed_coordinates{{0.0, 0.0, 0.0, 0.0}};

    std::string role;       ///< pilot, train, validation, final, or diagnostic.
    std::string status;     ///< running, converged, capped, or failed.
    double direct_transmission = 0.0;
    double direct_fluence = 0.0;
    std::uint64_t histories = 0;
    double wall_seconds = 0.0;
    double cpu_seconds = 0.0;
    double target_total_rel_se = 0.0;
    double achieved_total_rel_se = 0.0;
};

/// One independent logical replica/stream within a block.  Replica ids are
/// exactly 0..31 and index the leading dimension of EnvFitReplicaSums.
struct EnvFitReplicaRecord {
    std::uint32_t replica_id = 0;
    std::uint64_t history_begin = 0;  ///< Inclusive global history id.
    std::uint64_t history_count = 0;
    std::uint64_t rng_seed = 0;
};

struct EnvFitSourceRange {
    std::uint64_t block_id = 0;
    EnvFitReplicaRecord replica;
};

/// Raw sufficient statistics for per-primary-history contributions W_h.
/// Every vector has the same flattened tensor length.
struct EnvFitMomentArray {
    std::vector<std::uint64_t> n;
    std::vector<double> s1;
    std::vector<double> s2;
};

struct EnvFitScalarMoments {
    std::uint64_t n = 0;
    double s1 = 0.0;
    double s2 = 0.0;
};

/// Float32 raw S1 maps, indexed replica-major.  Counts and seeds live in the
/// replica records, so arbitrary rebinnings/folds can estimate covariance.
struct EnvFitReplicaSums {
    std::vector<float> continuum_components;
    std::vector<float> continuum_total;
    std::vector<float> line_components;
    std::vector<float> line_total;
};

/// Exact coarse observables used for stopping.  `outer_s2[i*16+j]` stores
/// sum_h x_h[i]*x_h[j], not a centered covariance.
struct EnvFitStoppingMoments {
    /// Canonical CEEFIT-STOPPING-v1 JSON defining all 16 vector entries.
    std::string definition_json;
    std::uint64_t n = 0;
    std::array<double, kEnvFitStoppingDimension> s1{};
    std::array<double,
               kEnvFitStoppingDimension * kEnvFitStoppingDimension> outer_s2{};
};

/// A diagnostic path-class map on the total energy axis plus exact line maps.
struct EnvFitDiagnosticChannel {
    std::string name;
    EnvFitMomentArray continuum;
    EnvFitMomentArray lines;
};

struct EnvFitDiagnostics {
    std::vector<EnvFitDiagnosticChannel> path_channels;
    std::vector<std::string> interaction_names;
    std::vector<std::uint64_t> interaction_counts;
    /// Exact integrated continuum-plus-lines per-history moments.
    EnvFitScalarMoments scattered_total;
    EnvFitScalarMoments dropped_weight;
    std::uint64_t transport_failures = 0;
    double maximum_history_contribution = 0.0;
    bool capped = false;
    std::string failure_reason;
};

/// One CEEFIT-BLOCK-v1 worker result.
///
/// Flattening:
///   continuum_components: concatenate component tensors in descriptor order;
///     each tensor is [energy_bin][mu_bin], mu fastest.
///   continuum_total: [energy_bin on total_energy_axis][mu_bin].
///   line_components: [component][line][mu_bin].
///   line_total: [line][mu_bin], accumulated independently so its per-history
///     S2 retains covariance between additive components.
/// Replica arrays add a leading [replica_id] dimension to those layouts.
struct EnvFitBlock {
    std::uint64_t node_id = 0;
    std::uint64_t block_id = 0;
    EnvFitManifest manifest;
    EnvFitAxes axes;
    std::vector<EnvFitComponent> components;
    EnvFitNodeMetadata node;
    std::array<EnvFitReplicaRecord, kEnvFitReplicaCount> replicas{};

    EnvFitMomentArray continuum_components;
    EnvFitMomentArray continuum_total;
    EnvFitMomentArray line_components;
    EnvFitMomentArray line_total;
    EnvFitReplicaSums replica_sums;
    EnvFitStoppingMoments stopping;
    EnvFitDiagnostics diagnostics;
};

struct EnvFitEstimate {
    double mean = 0.0;
    double standard_error = 0.0;
};

/// SHA-256 utilities used by campaign configuration and stable ids.
std::string env_fit_sha256_hex(const void* data, std::size_t size);
std::string env_fit_sha256_hex(const std::string& text);
std::uint64_t env_fit_stable_id(const std::string& canonical_key);
/// Cross-language node identity: SHA-256 domain separator, ASCII config hash,
/// and the five physical coordinates as little-endian binary64 values.
std::uint64_t env_fit_node_id(const std::string& config_sha256,
                              const EnvFitNodeMetadata& node);

/// Canonical self-describing definition for the 12 registered-index bands and
/// Legendre l=1..4 moments associated with an axis inventory.
std::string env_fit_stopping_definition_json(const EnvFitAxes& axes);

/// Validate the configuration hash and the required material/observable fields
/// before a worker spends time transporting a block.
bool validate_env_fit_manifest(const EnvFitManifest& manifest,
                               std::string* error = nullptr);

/// Expected flattened dimensions for a block's axes/component descriptors.
std::size_t env_fit_continuum_component_cells(const EnvFitBlock& block);
std::size_t env_fit_continuum_total_cells(const EnvFitBlock& block);
std::size_t env_fit_line_component_cells(const EnvFitBlock& block);
std::size_t env_fit_line_total_cells(const EnvFitBlock& block);

/// Validate all dimensions, the configuration hash, ids/ranges, axes, and
/// finite values.  On failure `error` receives a human-readable reason when
/// non-null.
bool validate_env_fit_block(const EnvFitBlock& block,
                            std::string* error = nullptr);

/// Atomically write an adjacent temporary file and rename it into place.
/// The file has a SHA-256 payload checksum; load rejects truncation, trailing
/// bytes, malformed dimensions, and checksum mismatches without modifying out.
bool save_env_fit_block_atomic(const EnvFitBlock& block,
                               const std::string& path,
                               std::string* error = nullptr);
bool load_env_fit_block(EnvFitBlock& out, const std::string& path,
                        std::string* error = nullptr);

/// Derive a mean and standard error from raw per-history sufficient statistics.
/// Returns NaN for the standard error when fewer than two histories are present.
EnvFitEstimate env_fit_estimate(std::uint64_t n, double s1, double s2);
EnvFitEstimate env_fit_estimate(const EnvFitMomentArray& moments,
                                std::size_t index);

/// Incremental in-memory reduction of disjoint blocks for one fit node.
/// The reducer retains every source range/id so continuation can reject replayed
/// work. `aggregate()` is a synthetic sum (block_id=0); save the original blocks,
/// not the aggregate, as CEEFIT-BLOCK-v1 worker records.
class EnvFitBlockReducer {
public:
    bool add(const EnvFitBlock& block, std::string* error = nullptr);
    bool empty() const { return block_ids_.empty(); }
    const EnvFitBlock& aggregate() const { return aggregate_; }
    const std::vector<std::uint64_t>& block_ids() const { return block_ids_; }
    const std::vector<EnvFitSourceRange>& history_ranges() const
    {
        return history_ranges_;
    }

private:
    EnvFitBlock aggregate_;
    std::vector<std::uint64_t> block_ids_;
    std::vector<EnvFitSourceRange> history_ranges_;
};

} // namespace ceelo
