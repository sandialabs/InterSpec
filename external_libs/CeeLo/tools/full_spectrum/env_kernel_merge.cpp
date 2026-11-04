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

/// @file env_kernel_merge.cpp
/// @brief Reassemble a distributed environment-scatter bake from its shards.
///
/// A production bake is split by node index across machines
/// (`env_kernel_bake --shard I --nshards N`). Because every node's RNG seed is
/// derived from its GLOBAL grid index, a sharded bake produces bit-identical
/// nodes to a single-machine one -- the split is bookkeeping, not physics.
/// This tool puts the pieces back together and, critically, REFUSES the
/// failure modes that would otherwise ship silently in a write-once product:
///
///   * shards from a different grid, material, or scoring radius;
///   * duplicate nodes (overlapping shard assignment -- which implies some
///     other node was never baked at all);
///   * a merged table with holes, reported as the specific missing
///     coordinates (degenerate |h_src - h_det| >= d corners are excluded from
///     the check, since the geometry axes are coupled and those points have
///     no field geometry to bake).
///
///   env_kernel_merge --out kernel.bin shard0.bin shard1.bin ...
///   env_kernel_merge --check kernel.bin        (completeness report only)
///   env_kernel_merge --diff A.bin B.bin        (node-by-node value compare)
///
/// `--diff` is how the shard-equivalence claim is verified, and how a spot
/// check on returned worker output is validated: re-bake one node locally and
/// confirm it matches the shard bit for bit.

#include "io/EnvScatterKernel.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace ceelo;

int main(int argc, char** argv)
{
    std::string out_path, check_path;
    bool diff_mode = false;
    std::vector<std::string> inputs;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_path = argv[++i];
        else if (a == "--check" && i + 1 < argc) check_path = argv[++i];
        else if (a == "--diff") diff_mode = true;
        else if (a == "--help" || a == "-h") {
            std::cout << "Usage: env_kernel_merge --out OUT.bin SHARD.bin ...\n"
                         "       env_kernel_merge --check TABLE.bin\n"
                         "       env_kernel_merge --diff A.bin B.bin\n";
            return 0;
        }
        else inputs.push_back(a);
    }

    if (diff_mode) {
        if (inputs.size() != 2) {
            std::cerr << "--diff needs exactly two files\n";
            return 1;
        }
        EnvKernelTable a, b;
        if (!load_env_kernel(a, inputs[0]) || !load_env_kernel(b, inputs[1])) {
            std::cerr << "cannot read inputs\n";
            return 1;
        }
        int compared = 0, unmatched = 0;
        double worst = 0.0;
        std::string worst_at;
        for (const EnvKernelNode& na : a.nodes) {
            const EnvKernelNode* nb = nullptr;
            for (const EnvKernelNode& c : b.nodes) {
                if (c.e0_keV == na.e0_keV && c.distance_cm == na.distance_cm &&
                    c.h_src_cm == na.h_src_cm && c.h_det_cm == na.h_det_cm &&
                    c.transverse_extent_cm == na.transverse_extent_cm) {
                    nb = &c;
                    break;
                }
            }
            if (!nb) { ++unmatched; continue; }
            ++compared;
            if (na.s_hat.size() != nb->s_hat.size()) {
                std::cout << "  SIZE MISMATCH at E0=" << na.e0_keV
                          << " d=" << na.distance_cm << "\n";
                worst = 1e300;
                continue;
            }
            for (std::size_t k = 0; k < na.s_hat.size(); ++k) {
                const double x = na.s_hat[k], y = nb->s_hat[k];
                const double denom = std::max(std::abs(x), std::abs(y));
                if (denom <= 0.0) continue;
                const double rel = std::abs(x - y) / denom;
                if (rel > worst) {
                    worst = rel;
                    worst_at = "E0=" + std::to_string(na.e0_keV) + " d=" +
                               std::to_string(na.distance_cm) + " bin " +
                               std::to_string(k);
                }
            }
        }
        std::cout << "compared " << compared << " matching node(s), "
                  << unmatched << " unmatched\n";
        std::cout << "worst relative s_hat difference: " << worst;
        if (!worst_at.empty()) std::cout << "  (" << worst_at << ")";
        std::cout << "\n";
        if (worst == 0.0 && unmatched == 0)
            std::cout << "IDENTICAL\n";
        return (worst == 0.0 && unmatched == 0) ? 0 : 2;
    }

    if (!check_path.empty()) {
        EnvKernelTable t;
        if (!load_env_kernel(t, check_path)) {
            std::cerr << "cannot read " << check_path << "\n";
            return 1;
        }
        std::string missing;
        const bool ok = env_kernel_is_complete(t, &missing);
        std::cout << check_path << ": " << t.nodes.size() << " nodes, kind '"
                  << t.kind << "'\n";
        std::cout << (ok ? "COMPLETE\n" : "INCOMPLETE -- missing:\n" + missing);
        return ok ? 0 : 2;
    }

    if (out_path.empty() || inputs.empty()) {
        std::cerr << "need --out OUT.bin and at least one shard\n";
        return 1;
    }

    EnvKernelTable merged;
    for (const std::string& p : inputs) {
        EnvKernelTable s;
        if (!load_env_kernel(s, p)) {
            std::cerr << "cannot read " << p << "\n";
            return 1;
        }
        if (!merge_env_kernel_shard(merged, s)) {
            std::cerr << "REFUSED to merge " << p
                      << ": incompatible grid/material, or a duplicate node "
                         "(overlapping shard assignment)\n";
            return 1;
        }
        std::cout << "  + " << p << " (" << s.nodes.size() << " nodes)\n";
    }

    std::string missing;
    const bool complete = env_kernel_is_complete(merged, &missing);
    if (!complete)
        std::cout << "WARNING: merged table is INCOMPLETE -- missing:\n"
                  << missing;

    if (!save_env_kernel_binary(merged, out_path)) {
        std::cerr << "FAILED to write " << out_path << "\n";
        return 1;
    }
    std::cout << "Wrote " << merged.nodes.size() << " nodes to " << out_path
              << (complete ? "  [COMPLETE]\n" : "  [INCOMPLETE]\n");
    return complete ? 0 : 2;
}
