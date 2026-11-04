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

/// @file EnvFitBlockBuilder.h
/// @brief Scheduling-independent production of one CEEFIT worker block.

#include "io/EnvFitArchive.h"
#include "transport/NextEventEstimator.h"

#include <cstdint>
#include <string>
#include <vector>

namespace ceelo {

struct EnvFitBlockRequest {
    EnvFitManifest manifest;
    std::uint64_t block_id = 1;
    std::string role = "train";

    double e0_keV = 0.0;
    double distance_cm = 0.0;
    double h_src_cm = 0.0;
    double h_det_cm = 0.0;
    double transverse_extent_cm = 0.0;

    std::uint64_t histories = 262144;
    /// Reserved global-history range per block. Partial smoke/pilot blocks may
    /// use fewer histories without changing later block identities.
    std::uint64_t history_slot_size = 262144;
    int threads = 1;
    double target_total_rel_se = 0.0025;
};

/// Stable identities used by both the worker and the Python campaign planner.
std::uint64_t env_fit_node_id(const std::string& config_sha256,
                              double e0_keV, double distance_cm,
                              double h_src_cm, double h_det_cm,
                              double transverse_extent_cm = 0.0);
std::uint64_t env_fit_replica_seed(const std::string& config_sha256,
                                   std::uint64_t node_id,
                                   std::uint64_t block_id,
                                   std::uint32_t replica_id);

/// Run exactly one 32-replica block. The component-native axes are derived from
/// `request.e0_keV`; `config.axes` is the common total registered axis.
/// Throws std::invalid_argument for a malformed request and std::runtime_error
/// if a worker fails.
EnvFitBlock run_env_fit_block(const EnvironmentGeometry& environment,
                              const EnvScatterConfig& config,
                              const Eigen::Vector3d& source,
                              const EnvFitBlockRequest& request);

/// Run consecutive or otherwise grouped blocks in one process.  Results retain
/// request order, while replica tasks use the same persistent worker pool as
/// run_env_fit_block().  Queue workers should prefer this entry point (or keep
/// their process alive) so thread construction is amortized across blocks.
std::vector<EnvFitBlock> run_env_fit_blocks(
    const EnvironmentGeometry& environment,
    const EnvScatterConfig& config,
    const Eigen::Vector3d& source,
    const std::vector<EnvFitBlockRequest>& requests);

} // namespace ceelo
