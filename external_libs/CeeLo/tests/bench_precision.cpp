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

/// @file bench_precision.cpp
/// @brief Benchmark: precision-based termination vs fixed event counts.
///
/// Compares the precision-based SimulationConfig (target_fep_rel_precision = 0.0005)
/// against fixed-count compute() at 20000 events for a 3"x3" NaI detector, point
/// source at z = -10 cm, across six energies.

#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace ceelo;

static void print_row(double E, const EfficiencyResult& result) {
    double rel_prec = (result.full_energy_peak_efficiency > 0.0)
        ? result.fep_uncertainty / result.full_energy_peak_efficiency
        : 0.0;

    const char* stop_str = "???";
    switch (result.stop_reason) {
        case StopReason::MaxEvents:       stop_str = "MaxEvents"; break;
        case StopReason::FepPrecision:    stop_str = "FepPrec"; break;
        case StopReason::TotalPrecision:  stop_str = "TotalPrec"; break;
        case StopReason::WallTime:        stop_str = "WallTime"; break;
        case StopReason::CpuTime:         stop_str = "CpuTime"; break;
    }

    std::printf("%-10.0f %12llu %14.6e %14.6e %12.6f %12.3f %10s\n",
                E,
                (unsigned long long)result.num_events_simulated,
                result.full_energy_peak_efficiency,
                result.fep_uncertainty,
                rel_prec,
                result.wall_time_seconds,
                stop_str);
    std::fflush(stdout);
}

static void print_header() {
    std::printf("%-10s %12s %14s %14s %12s %12s %10s\n",
                "E (keV)", "N_events", "FEP_eff", "FEP_unc",
                "rel_prec", "time (s)", "stop");
    std::fflush(stdout);
}

int main() {
    // Disable stdout buffering entirely
    std::setbuf(stdout, nullptr);

    // --- Detector setup: 3"x3" NaI cylinder ---
    // 3 inches = 7.62 cm; radius = 3.81 cm, length = 7.62 cm
    Material nai = make_NaI();

    const std::vector<double> energies = {100.0, 200.0, 351.0, 662.0, 1173.0, 2614.0};

    // =============================================
    // Part 1: Precision-based termination at various targets
    // =============================================
    for (double target : {0.05, 0.02, 0.01, 0.005}) {
        std::printf("\n=== Precision-based termination (target rel precision = %.1f%%) ===\n",
                    target * 100.0);
        print_header();

        for (double E : energies) {
            EfficiencyCalculator calc;
            calc.set_fep_window_keV(kTestFepWindowKeV);
            calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

            SimulationConfig config;
            config.energy_keV = E;
            config.termination.target_fep_rel_precision = target;
            config.termination.max_events = 50000000;
            config.termination.min_events = 1000;
            config.termination.max_wall_seconds = 60.0;
            config.num_threads = 1;
            config.batch_size = 5000;

            auto result = calc.compute(config);
            print_row(E, result);
        }
    }

    // =============================================
    // Part 2: Fixed event count (20,000 events)
    // =============================================
    std::printf("\n=== Fixed event count (N = 20000) ===\n");
    print_header();

    for (double E : energies) {
        EfficiencyCalculator calc;
        calc.set_fep_window_keV(kTestFepWindowKeV);
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

        SimulationConfig config;
        config.energy_keV = E;
        config.termination.max_events = 20000;
        config.termination.target_fep_rel_precision = 0.0; // disabled
        config.termination.min_events = 20000;
        config.num_threads = 1;
        config.batch_size = 5000;

        auto result = calc.compute(config);
        print_row(E, result);
    }

    return 0;
}
