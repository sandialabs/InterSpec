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

/// @file env_fit_block_tool.cpp
/// @brief Verify and summarize one or more CEEFIT worker blocks.

#include "io/EnvFitArchive.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

EnvFitEstimate total_estimate(const EnvFitBlock& block)
{
    double s1 = 0.0;
    double s2 = 0.0;
    for (std::size_t i = 0; i < 12; ++i) {
        s1 += block.stopping.s1[i];
        for (std::size_t j = 0; j < 12; ++j)
            s2 += block.stopping.outer_s2[
                i * kEnvFitStoppingDimension + j];
    }
    return env_fit_estimate(block.stopping.n, s1, s2);
}

void print_summary(const EnvFitBlock& block, const std::string& label)
{
    const EnvFitEstimate total = total_estimate(block);
    std::cout << label << ": node=" << std::hex << std::setw(16)
              << std::setfill('0') << block.node_id << std::dec
              << std::setfill(' ') << " block=" << block.block_id
              << " role=" << block.node.role
              << " status=" << block.node.status
              << " histories=" << block.node.histories
              << " E0=" << block.node.e0_keV
              << "keV d=" << block.node.distance_cm << "cm"
              << " total=" << std::scientific << total.mean
              << " se=" << total.standard_error;
    if (total.mean > 0.0 && std::isfinite(total.standard_error))
        std::cout << " rel=" << total.standard_error / total.mean;
    std::cout << " config=" << block.manifest.config_sha256.substr(0, 16)
              << "...\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: env_fit_block_tool BLOCK.ceefit [BLOCK.ceefit ...]\n";
        return 2;
    }
    EnvFitBlockReducer reducer;
    for (int i = 1; i < argc; ++i) {
        EnvFitBlock block;
        std::string error;
        if (!load_env_fit_block(block, argv[i], &error)) {
            std::cerr << argv[i] << ": " << error << "\n";
            return 1;
        }
        print_summary(block, argv[i]);
        if (!reducer.add(block, &error)) {
            std::cerr << argv[i] << ": cannot merge: " << error << "\n";
            return 1;
        }
    }
    if (argc > 2) print_summary(reducer.aggregate(), "aggregate");
    return 0;
}
