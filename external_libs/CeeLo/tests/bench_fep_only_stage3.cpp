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

// bench_fep_only_stage3 -- unbiasedness + figure-of-merit bench for running
// the stage-3 (near-field) generation nodes in fep_only mode.
//
// The stage-3 tensor scan consumes ONLY full_energy_peak_efficiency (the
// eps_tot tally of those nodes is discarded), so FEP-only transport there
// discards nothing -- IF its eps_fep estimate is unbiased. This bench runs
// the matrix
//     {nai3x3, hpge_coax} x {122, 662, 1332, 2614 keV}
//   x {far d=10a ct=1, contact d=0.15a ct=1, graze d=0.5a ct=0.2}
//   x {full, fep_only}
// with identical per-point seeds and budgets, and reports per point:
//   eps ratio (fep_only/full), combined sigma, CPU-s ratio, and the FOM
//   gain, FOM = 1 / (rel_sigma^2 * cpu_s).
//
// Gate (spec C2): |eps_fep(fep_only)/eps_fep(full) - 1| < 2 * combined sigma
// at EVERY matrix point; exit status 1 if any point fails.
//
// Usage: bench_fep_only_stage3 [out_csv] [precision=0.003] [max_cpu_s=60]

#include "efficiency/EfficiencyCalculator.h"
#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "materials/Material.h"

#include "response_presets.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

// Deterministic per-point seed (never 0); mirrors ResponseGenerator's
// node_seed with a bench-private stage id.
uint64_t bench_seed(uint32_t point) {
    uint64_t s = 0xFEB0517A6E000000ULL ^ point;
    s ^= s >> 30;
    s *= 0xBF58476D1CE4E5B9ULL;
    s ^= s >> 27;
    return s | 1ULL;
}

Eigen::Vector3d src_pos(double d_cm, double cos_theta) {
    const double st = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    return Eigen::Vector3d(d_cm * st, 0.0, -d_cm * cos_theta);
}

}  // namespace

int main(int argc, char** argv) {
    const std::string out_csv = argc > 1 ? argv[1] : "";
    const double precision = argc > 2 ? std::atof(argv[2]) : 0.003;
    const double max_cpu_s = argc > 3 ? std::atof(argv[3]) : 60.0;

    const std::vector<std::string> presets = {"nai3x3", "hpge_coax"};
    const std::vector<double> energies = {122.0, 662.0, 1332.0, 2614.0};
    struct Geom { const char* name; double d_over_a, ct; };
    const std::vector<Geom> geoms = {
        {"far",     10.0, 1.0},   // far-field, on-axis
        {"contact", 0.15, 1.0},   // near-contact, on-axis
        {"graze",   0.5,  0.2},   // near-grazing (worst near-field corner)
    };

    FILE* csv = nullptr;
    if (!out_csv.empty()) {
        csv = std::fopen(out_csv.c_str(), "w");
        if (!csv) { std::printf("cannot write %s\n", out_csv.c_str()); return 1; }
        std::fprintf(csv, "preset,E_keV,geom,d_cm,cos_theta,"
                     "eps_full,unc_full,cpu_full,events_full,"
                     "eps_fep,unc_fep,cpu_fep,events_fep,"
                     "ratio,comb_rel_sigma,n_sigma,cpu_ratio,fom_gain,pass\n");
    }

    std::printf("bench_fep_only_stage3: precision %.4g, max_cpu %.0f CPU-s "
                "per point\n", precision, max_cpu_s);
    std::printf("%-10s %5s %-8s | %11s %11s | %7s %7s %7s %8s %5s\n",
                "preset", "E", "geom", "eps_full", "eps_fep",
                "dev%", "2sig%", "cpuX", "fomX", "gate");

    bool all_pass = true;
    uint32_t point = 0;
    for (const std::string& preset : presets) {
        const GeometryDescriptor gd = preset_descriptor(preset);
        const double a = gd.transverse_half_extent();
        for (const double E : energies) {
            for (const Geom& g : geoms) {
                const double d = g.d_over_a * a;
                const Eigen::Vector3d src = src_pos(d, g.ct);
                const uint64_t seed = bench_seed(point++);

                EfficiencyResult res[2];  // [0]=full, [1]=fep_only
                for (int mode = 0; mode < 2; ++mode) {
                    EfficiencyCalculator calc;
                    std::vector<std::unique_ptr<Material>> mats;
                    ResponseGenerator::configure_calculator(calc, gd, mats);
                    calc.set_point_source(src);
                    calc.enable_fep_only_mode(mode == 1);
                    SimulationConfig cfg;
                    cfg.energy_keV = E;
                    cfg.termination.target_fep_rel_precision = precision;
                    cfg.termination.max_events = 20000000;
                    cfg.termination.max_cpu_seconds = max_cpu_s;
                    cfg.termination.min_events = 20000;
                    cfg.batch_size = 2000;
                    cfg.seed = seed;  // same seed for the full/fep pair
                    res[mode] = calc.compute(cfg);
                }

                const double eps_f = res[0].full_energy_peak_efficiency;
                const double eps_o = res[1].full_energy_peak_efficiency;
                const double rel_f = eps_f > 0.0 ? res[0].fep_uncertainty / eps_f : 0.0;
                const double rel_o = eps_o > 0.0 ? res[1].fep_uncertainty / eps_o : 0.0;
                const double comb = std::sqrt(rel_f * rel_f + rel_o * rel_o);
                const double ratio = eps_f > 0.0 ? eps_o / eps_f : 0.0;
                const double n_sig = comb > 0.0 ? std::fabs(ratio - 1.0) / comb : 0.0;
                const double cpu_f = res[0].cpu_time_seconds;
                const double cpu_o = res[1].cpu_time_seconds;
                const double cpu_ratio = cpu_o > 0.0 ? cpu_f / cpu_o : 0.0;
                const double fom_f = (rel_f > 0.0 && cpu_f > 0.0)
                    ? 1.0 / (rel_f * rel_f * cpu_f) : 0.0;
                const double fom_o = (rel_o > 0.0 && cpu_o > 0.0)
                    ? 1.0 / (rel_o * rel_o * cpu_o) : 0.0;
                const double fom_gain = fom_f > 0.0 ? fom_o / fom_f : 0.0;
                const bool pass = std::fabs(ratio - 1.0) < 2.0 * comb;
                all_pass = all_pass && pass;

                std::printf("%-10s %5.0f %-8s | %11.4e %11.4e | %+6.2f%% %6.2f%% "
                            "%6.2fx %7.2fx %5s\n",
                            preset.c_str(), E, g.name, eps_f, eps_o,
                            100.0 * (ratio - 1.0), 200.0 * comb,
                            cpu_ratio, fom_gain, pass ? "pass" : "FAIL");
                std::fflush(stdout);
                if (csv)
                    std::fprintf(csv, "%s,%.1f,%s,%.4f,%.3f,"
                                 "%.8e,%.3e,%.3f,%llu,%.8e,%.3e,%.3f,%llu,"
                                 "%.6f,%.6f,%.3f,%.4f,%.4f,%d\n",
                                 preset.c_str(), E, g.name, d, g.ct,
                                 eps_f, res[0].fep_uncertainty, cpu_f,
                                 (unsigned long long)res[0].num_events_simulated,
                                 eps_o, res[1].fep_uncertainty, cpu_o,
                                 (unsigned long long)res[1].num_events_simulated,
                                 ratio, comb, n_sig, cpu_ratio, fom_gain,
                                 pass ? 1 : 0);
                if (csv) std::fflush(csv);
            }
        }
    }
    if (csv) std::fclose(csv);

    std::printf("bench_fep_only_stage3: %s\n",
                all_pass ? "ALL POINTS PASS (|ratio-1| < 2 sigma)"
                         : "GATE FAILED at one or more points");
    return all_pass ? 0 : 1;
}
