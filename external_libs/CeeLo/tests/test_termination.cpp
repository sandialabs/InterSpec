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

#define BOOST_TEST_MODULE TerminationTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "materials/Material.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <vector>

using namespace ceelo;

// ============================================================
//  Feature 1: Flexible Termination + Progress Callback
// ============================================================

BOOST_AUTO_TEST_SUITE(TerminationCriteria)

BOOST_AUTO_TEST_CASE(max_events_terminates_at_count) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    SimulationConfig config;
    config.energy_keV = 662.0;
    config.termination.max_events = 50000;
    config.num_threads = 2;
    config.batch_size = 5000;

    auto result = calc.compute(config);

    // Should stop at approximately max_events (may slightly overshoot by up to batch_size * num_threads)
    BOOST_CHECK_GE(result.num_events_simulated, 50000u);
    BOOST_CHECK_LE(result.num_events_simulated, 50000u + 5000u * 2u);
    BOOST_CHECK(result.stop_reason == StopReason::MaxEvents);
}

BOOST_AUTO_TEST_CASE(fep_precision_converges) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    SimulationConfig config;
    config.energy_keV = 662.0;
    config.termination.target_fep_rel_precision = 0.02; // 2%
    config.termination.min_events = 10000;
    config.termination.max_events = 10000000; // safety cap
    config.num_threads = 2;
    config.batch_size = 10000;

    auto result = calc.compute(config);

    BOOST_CHECK(result.stop_reason == StopReason::FepPrecision);
    // Verify achieved precision is close to target
    double rel_prec = result.fep_uncertainty / result.full_energy_peak_efficiency;
    BOOST_CHECK_LE(rel_prec, 0.025); // within 25% of target (batch overshoot)
}

BOOST_AUTO_TEST_CASE(wall_time_terminates) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    SimulationConfig config;
    config.energy_keV = 662.0;
    config.termination.max_wall_seconds = 0.5; // 500 ms
    config.termination.max_events = 100000000; // huge cap
    config.num_threads = 0;  // auto-detect
    config.batch_size = 1000;

    auto result = calc.compute(config);

    BOOST_CHECK(result.stop_reason == StopReason::WallTime);
    BOOST_CHECK_GE(result.wall_time_seconds, 0.4); // at least close
    BOOST_CHECK_LE(result.wall_time_seconds, 2.0);  // not too long
}

BOOST_AUTO_TEST_CASE(cpu_time_terminates) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    SimulationConfig config;
    config.energy_keV = 662.0;
    config.termination.max_cpu_seconds = 0.5;
    config.termination.max_events = 100000000; // huge cap
    config.num_threads = 2;
    config.batch_size = 1000;

    auto result = calc.compute(config);

    BOOST_CHECK(result.stop_reason == StopReason::CpuTime);
    BOOST_CHECK_GE(result.cpu_time_seconds, 0.4);  // reached the cap
    BOOST_CHECK_LE(result.cpu_time_seconds, 4.0);  // no runaway overshoot
    // Summed thread CPU should not wildly exceed wall x threads.
    BOOST_CHECK_LE(result.cpu_time_seconds,
                   result.wall_time_seconds * 2.0 + 0.5);
}

BOOST_AUTO_TEST_CASE(cpu_time_reported_without_cap) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    SimulationConfig config;
    config.energy_keV = 662.0;
    config.termination.max_events = 50000;
    config.num_threads = 2;
    config.batch_size = 5000;

    auto result = calc.compute(config);

    BOOST_CHECK(result.stop_reason == StopReason::MaxEvents);
    BOOST_CHECK_GT(result.cpu_time_seconds, 0.0);
}

BOOST_AUTO_TEST_CASE(progress_callback_invoked) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    std::atomic<int> callback_count{0};
    uint64_t last_events = 0;
    Progress last;

    SimulationConfig config;
    config.energy_keV = 662.0;
    config.termination.max_events = 50000;
    config.num_threads = 0;  // auto-detect
    config.batch_size = 5000;
    config.progress_callback = [&](const Progress& p) {
        callback_count++;
        BOOST_CHECK_GE(p.num_events, last_events);
        last_events = p.num_events;
        BOOST_CHECK_GE(p.frac_complete, 0.0);
        BOOST_CHECK_LE(p.frac_complete, 1.0);
        BOOST_CHECK_GE(p.fep_uncertainty, 0.0);
        BOOST_CHECK_GE(p.total_uncertainty, 0.0);
        BOOST_CHECK_GE(p.total_efficiency, p.fep_efficiency - 1e-12);
        last = p;
    };

    auto result = calc.compute(config);

    BOOST_CHECK_GT(callback_count.load(), 0);
    // Completion fire: guaranteed, payload identical to the returned result.
    BOOST_CHECK_EQUAL(last.num_events, result.num_events_simulated);
    BOOST_CHECK_EQUAL(last.frac_complete, 1.0);
    BOOST_CHECK_EQUAL(last.fep_efficiency, result.full_energy_peak_efficiency);
    BOOST_CHECK_EQUAL(last.fep_uncertainty, result.fep_uncertainty);
    BOOST_CHECK_EQUAL(last.total_efficiency, result.total_efficiency);
    BOOST_CHECK_EQUAL(last.total_uncertainty, result.total_uncertainty);
}

BOOST_AUTO_TEST_CASE(old_compute_still_works) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    auto result = calc.compute(662.0, 10000, 1);

    BOOST_CHECK_EQUAL(result.num_events_simulated, 10000u);
    BOOST_CHECK_GT(result.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_GE(result.total_efficiency,
                   result.full_energy_peak_efficiency - 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()
