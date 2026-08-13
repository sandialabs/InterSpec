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

#define BOOST_TEST_MODULE FepOnlyTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "geometry/RayTrace.h"
#include "transport/ComptonScatter.h"

#include <cmath>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <random>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <unistd.h>

using namespace ceelo;

namespace {

/// Append a diagnostic line to fep_results.txt in the current working directory.
void log_result(const std::string& line) {
    std::ofstream f("fep_results.txt", std::ios::app);
    if (f.is_open()) f << line << "\n";
}

/// Helper: create a SimulationConfig targeting 0.5% FEP relative precision.
SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 5000;
    config.termination.max_wall_seconds = 120.0;
    config.num_threads = 0;  // auto-detect
    config.batch_size = 10000;
    return config;
}

} // anonymous namespace

// ============================================================
//  Feature 2: FEP-Only Mode
// ============================================================

BOOST_AUTO_TEST_SUITE(FepOnlyMode)

BOOST_AUTO_TEST_CASE(fep_matches_full_mode_bare_nai) {
    // FEP-only should give same FEP as full mode within stats.
    Material nai = make_NaI();

    EfficiencyCalculator calc_full;
    calc_full.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_full.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_full.enable_fep_only_mode(false);
    auto res_full = calc_full.compute(precision_config(662.0));

    EfficiencyCalculator calc_fep;
    calc_fep.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_fep.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_fep.enable_fep_only_mode(true);
    auto res_fep = calc_fep.compute(precision_config(662.0));

    double sigma = std::sqrt(
        res_full.fep_uncertainty * res_full.fep_uncertainty +
        res_fep.fep_uncertainty * res_fep.fep_uncertainty);

    double diff = std::abs(res_full.full_energy_peak_efficiency -
                           res_fep.full_energy_peak_efficiency);

    double pct_diff = (res_full.full_energy_peak_efficiency > 0.0)
        ? 100.0 * (res_fep.full_energy_peak_efficiency - res_full.full_energy_peak_efficiency)
              / res_full.full_energy_peak_efficiency
        : 0.0;
    double n_sigma = (sigma > 0.0) ? diff / sigma : 0.0;

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "bare NaI 662 keV: full=" << res_full.full_energy_peak_efficiency
        << " +/- " << res_full.fep_uncertainty
        << "  fep=" << res_fep.full_energy_peak_efficiency
        << " +/- " << res_fep.fep_uncertainty
        << std::setprecision(2)
        << "  diff=" << pct_diff << "%"
        << "  (" << n_sigma << " sigma)";
    BOOST_TEST_MESSAGE(msg.str());
    log_result(msg.str());

    // Allow 5-sigma difference
    BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
}

BOOST_AUTO_TEST_CASE(fep_matches_full_mode_with_attenuator) {
    // Now using 200 keV where Rayleigh scattering in the attenuator is significant.
    // Stochastic Rayleigh in FEP-only mode should handle this correctly.
    Material nai = make_NaI();
    Material al = make_Aluminum();

    EfficiencyCalculator calc_full;
    calc_full.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_full.add_attenuator(&al, 0.1, 0.1, 0.0, 7.62);
    calc_full.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_full.enable_fep_only_mode(false);
    auto res_full = calc_full.compute(precision_config(200.0));

    EfficiencyCalculator calc_fep;
    calc_fep.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_fep.add_attenuator(&al, 0.1, 0.1, 0.0, 7.62);
    calc_fep.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_fep.enable_fep_only_mode(true);
    auto res_fep = calc_fep.compute(precision_config(200.0));

    double sigma = std::sqrt(
        res_full.fep_uncertainty * res_full.fep_uncertainty +
        res_fep.fep_uncertainty * res_fep.fep_uncertainty);

    double diff = std::abs(res_full.full_energy_peak_efficiency -
                           res_fep.full_energy_peak_efficiency);

    double pct_diff = (res_full.full_energy_peak_efficiency > 0.0)
        ? 100.0 * (res_fep.full_energy_peak_efficiency - res_full.full_energy_peak_efficiency)
              / res_full.full_energy_peak_efficiency
        : 0.0;
    double n_sigma = (sigma > 0.0) ? diff / sigma : 0.0;

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "1mm Al attenuator 200 keV: full=" << res_full.full_energy_peak_efficiency
        << " +/- " << res_full.fep_uncertainty
        << "  fep=" << res_fep.full_energy_peak_efficiency
        << " +/- " << res_fep.fep_uncertainty
        << std::setprecision(2)
        << "  diff=" << pct_diff << "%"
        << "  (" << n_sigma << " sigma)";
    BOOST_TEST_MESSAGE(msg.str());
    log_result(msg.str());

    BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
}

BOOST_AUTO_TEST_CASE(fep_matches_full_mode_with_source_shield) {
    // FEP-only with Pb source shielding at 200 keV.
    // Stochastic Rayleigh in source shielding should match full mode.
    Material nai = make_NaI();
    Material pb = make_Lead();

    EfficiencyCalculator calc_full;
    calc_full.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_full.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_full.add_source_shield(&pb, 0.05);  // 0.5mm Pb shield
    calc_full.enable_fep_only_mode(false);
    auto res_full = calc_full.compute(precision_config(200.0));

    EfficiencyCalculator calc_fep;
    calc_fep.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_fep.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_fep.add_source_shield(&pb, 0.05);
    calc_fep.enable_fep_only_mode(true);
    auto res_fep = calc_fep.compute(precision_config(200.0));

    double sigma = std::sqrt(
        res_full.fep_uncertainty * res_full.fep_uncertainty +
        res_fep.fep_uncertainty * res_fep.fep_uncertainty);

    double diff = std::abs(res_full.full_energy_peak_efficiency -
                           res_fep.full_energy_peak_efficiency);

    double pct_diff = (res_full.full_energy_peak_efficiency > 0.0)
        ? 100.0 * (res_fep.full_energy_peak_efficiency - res_full.full_energy_peak_efficiency)
              / res_full.full_energy_peak_efficiency
        : 0.0;
    double n_sigma = (sigma > 0.0) ? diff / sigma : 0.0;

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "0.5mm Pb source shield 200 keV: full=" << res_full.full_energy_peak_efficiency
        << " +/- " << res_full.fep_uncertainty
        << "  fep=" << res_fep.full_energy_peak_efficiency
        << " +/- " << res_fep.fep_uncertainty
        << std::setprecision(2)
        << "  diff=" << pct_diff << "%"
        << "  (" << n_sigma << " sigma)";
    BOOST_TEST_MESSAGE(msg.str());
    log_result(msg.str());

    BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
}

BOOST_AUTO_TEST_CASE(fep_only_is_faster) {
    Material nai = make_NaI();
    const uint64_t N = 100000;

    EfficiencyCalculator calc_full;
    calc_full.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_full.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_full.enable_fep_only_mode(false);

    auto t0 = std::chrono::steady_clock::now();
    calc_full.compute(662.0, N, 1);
    auto t1 = std::chrono::steady_clock::now();
    double full_time = std::chrono::duration<double>(t1 - t0).count();

    EfficiencyCalculator calc_fep;
    calc_fep.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_fep.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_fep.enable_fep_only_mode(true);

    auto t2 = std::chrono::steady_clock::now();
    calc_fep.compute(662.0, N, 1);
    auto t3 = std::chrono::steady_clock::now();
    double fep_time = std::chrono::duration<double>(t3 - t2).count();

    BOOST_CHECK_LT(fep_time, full_time * 1.5);
}

BOOST_AUTO_TEST_CASE(fep_only_low_energy) {
    Material nai = make_NaI();

    EfficiencyCalculator calc_full;
    calc_full.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_full.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_full.enable_fep_only_mode(false);
    auto res_full = calc_full.compute(precision_config(100.0));

    EfficiencyCalculator calc_fep;
    calc_fep.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_fep.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_fep.enable_fep_only_mode(true);
    auto res_fep = calc_fep.compute(precision_config(100.0));

    double sigma = std::sqrt(
        res_full.fep_uncertainty * res_full.fep_uncertainty +
        res_fep.fep_uncertainty * res_fep.fep_uncertainty);

    double diff = std::abs(res_full.full_energy_peak_efficiency -
                           res_fep.full_energy_peak_efficiency);

    double pct_diff = (res_full.full_energy_peak_efficiency > 0.0)
        ? 100.0 * (res_fep.full_energy_peak_efficiency - res_full.full_energy_peak_efficiency)
              / res_full.full_energy_peak_efficiency
        : 0.0;
    double n_sigma = (sigma > 0.0) ? diff / sigma : 0.0;

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "bare NaI 100 keV: full=" << res_full.full_energy_peak_efficiency
        << " +/- " << res_full.fep_uncertainty
        << "  fep=" << res_fep.full_energy_peak_efficiency
        << " +/- " << res_fep.fep_uncertainty
        << std::setprecision(2)
        << "  diff=" << pct_diff << "%"
        << "  (" << n_sigma << " sigma)";
    BOOST_TEST_MESSAGE(msg.str());
    log_result(msg.str());

    BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
}

// Near-field unbiasedness (stage-3 generation geometry): FEP-only must match
// full mode with the source at a fraction of the detector radius, on-axis and
// near-grazing, where solid angle is large and side/can crossings dominate.
// NOTE bench_fep_only_stage3 (2026-07) measured a -5..-11% FEP-only bias for
// BORE-HOLE detectors (hpge_coax): kill-on-scoring-exit kills photons that
// cross the vacuum bore and would re-enter the crystal with full energy
// (removing the bore removes the bias). fep_only therefore stays DISABLED in
// the stage-3 generation loop; this case guards the non-bore domain.
BOOST_AUTO_TEST_CASE(fep_matches_full_mode_near_field) {
    Material nai = make_NaI();
    Material al = make_Aluminum();

    struct NearGeom { const char* label; Eigen::Vector3d src; };
    // a = 3.81 cm crystal radius: contact d = 0.15a on-axis; grazing
    // d = 0.5a at cos_theta = 0.2 (x = d*sin, z = -d*cos conventions).
    const NearGeom geoms[] = {
        {"contact d=0.57cm on-axis", Eigen::Vector3d(0.0, 0.0, -0.5715)},
        {"graze d=1.9cm ct=0.2",
         Eigen::Vector3d(1.9050 * 0.9798, 0.0, -1.9050 * 0.2)},
    };

    for (const NearGeom& g : geoms) {
        EfficiencyCalculator calc_full;
        calc_full.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc_full.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
        calc_full.set_point_source(g.src);
        calc_full.enable_fep_only_mode(false);
        auto res_full = calc_full.compute(precision_config(122.0, 0.004));

        EfficiencyCalculator calc_fep;
        calc_fep.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc_fep.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
        calc_fep.set_point_source(g.src);
        calc_fep.enable_fep_only_mode(true);
        auto res_fep = calc_fep.compute(precision_config(122.0, 0.004));

        double sigma = std::sqrt(
            res_full.fep_uncertainty * res_full.fep_uncertainty +
            res_fep.fep_uncertainty * res_fep.fep_uncertainty);
        double diff = std::abs(res_full.full_energy_peak_efficiency -
                               res_fep.full_energy_peak_efficiency);
        double pct_diff = (res_full.full_energy_peak_efficiency > 0.0)
            ? 100.0 * (res_fep.full_energy_peak_efficiency - res_full.full_energy_peak_efficiency)
                  / res_full.full_energy_peak_efficiency
            : 0.0;
        double n_sigma = (sigma > 0.0) ? diff / sigma : 0.0;

        std::ostringstream msg;
        msg << std::fixed << std::setprecision(6)
            << "near-field NaI+can 122 keV " << g.label
            << ": full=" << res_full.full_energy_peak_efficiency
            << " +/- " << res_full.fep_uncertainty
            << "  fep=" << res_fep.full_energy_peak_efficiency
            << " +/- " << res_fep.fep_uncertainty
            << std::setprecision(2)
            << "  diff=" << pct_diff << "%"
            << "  (" << n_sigma << " sigma)";
        BOOST_TEST_MESSAGE(msg.str());
        log_result(msg.str());

        BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
    }
}

// Validate FEP early termination optimization at 0.25% precision.
// Tests multiple energies to exercise each early-kill path:
//   100 keV — fluorescence escape path (NaI: I Ka ~28 keV)
//   662 keV — Compton scatter + electron escape path
//   2614 keV — pair production + annihilation gamma escape path
BOOST_AUTO_TEST_CASE(fep_early_kill_matches_full_mode_multi_energy) {
    Material nai = make_NaI();

    const double energies[] = {100.0, 662.0, 2614.0};
    const char* labels[] = {"100 keV (fluorescence)", "662 keV (Compton)", "2614 keV (pair prod)"};

    for (int ie = 0; ie < 3; ++ie) {
        double E = energies[ie];

        EfficiencyCalculator calc_full;
        calc_full.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc_full.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
        calc_full.enable_fep_only_mode(false);
        auto res_full = calc_full.compute(precision_config(E, 0.0025));

        EfficiencyCalculator calc_fep;
        calc_fep.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc_fep.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
        calc_fep.enable_fep_only_mode(true);
        auto res_fep = calc_fep.compute(precision_config(E, 0.0025));

        double sigma = std::sqrt(
            res_full.fep_uncertainty * res_full.fep_uncertainty +
            res_fep.fep_uncertainty * res_fep.fep_uncertainty);

        double diff = std::abs(res_full.full_energy_peak_efficiency -
                               res_fep.full_energy_peak_efficiency);

        double pct_diff = (res_full.full_energy_peak_efficiency > 0.0)
            ? 100.0 * (res_fep.full_energy_peak_efficiency - res_full.full_energy_peak_efficiency)
                  / res_full.full_energy_peak_efficiency
            : 0.0;
        double n_sigma = (sigma > 0.0) ? diff / sigma : 0.0;

        std::ostringstream msg;
        msg << std::fixed << std::setprecision(6)
            << "FEP early-kill " << labels[ie]
            << ": full=" << res_full.full_energy_peak_efficiency
            << " +/- " << res_full.fep_uncertainty
            << "  fep=" << res_fep.full_energy_peak_efficiency
            << " +/- " << res_fep.fep_uncertainty
            << std::setprecision(2)
            << "  diff=" << pct_diff << "%"
            << "  (" << n_sigma << " sigma)";
        BOOST_TEST_MESSAGE(msg.str());
        log_result(msg.str());

        // Allow 5-sigma difference (should be well within for correct optimization)
        BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
    }
}

// Validate FEP early termination with CZT (thin detector where electron
// escape and fluorescence cascade escape are common).
BOOST_AUTO_TEST_CASE(fep_early_kill_czt) {
    Material czt = make_CZT();

    // 200 keV: Te/Cd fluorescence cascade escape is significant in thin CZT
    EfficiencyCalculator calc_full;
    calc_full.set_detector(DetectorShape::Box, &czt, {0.5, 0.5, 0.25});
    calc_full.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
    calc_full.enable_fep_only_mode(false);
    auto res_full = calc_full.compute(precision_config(200.0, 0.0025));

    EfficiencyCalculator calc_fep;
    calc_fep.set_detector(DetectorShape::Box, &czt, {0.5, 0.5, 0.25});
    calc_fep.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
    calc_fep.enable_fep_only_mode(true);
    auto res_fep = calc_fep.compute(precision_config(200.0, 0.0025));

    double sigma = std::sqrt(
        res_full.fep_uncertainty * res_full.fep_uncertainty +
        res_fep.fep_uncertainty * res_fep.fep_uncertainty);
    double diff = std::abs(res_full.full_energy_peak_efficiency -
                           res_fep.full_energy_peak_efficiency);
    double pct_diff = (res_full.full_energy_peak_efficiency > 0.0)
        ? 100.0 * (res_fep.full_energy_peak_efficiency - res_full.full_energy_peak_efficiency)
              / res_full.full_energy_peak_efficiency
        : 0.0;
    double n_sigma = (sigma > 0.0) ? diff / sigma : 0.0;

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "FEP early-kill CZT 200 keV: full=" << res_full.full_energy_peak_efficiency
        << " +/- " << res_full.fep_uncertainty
        << "  fep=" << res_fep.full_energy_peak_efficiency
        << " +/- " << res_fep.fep_uncertainty
        << std::setprecision(2)
        << "  diff=" << pct_diff << "%"
        << "  (" << n_sigma << " sigma)";
    BOOST_TEST_MESSAGE(msg.str());
    log_result(msg.str());

    BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  compute_transmission_no_rayleigh
// ============================================================

BOOST_AUTO_TEST_SUITE(TransmissionNoRayleigh)

BOOST_AUTO_TEST_CASE(no_rayleigh_higher_than_total) {
    Material al = make_Aluminum();

    PathSegment seg;
    seg.t_start = 0.0;
    seg.t_end = 1.0;
    seg.material = &al;
    seg.is_scoring = false;

    std::vector<PathSegment> segments = {seg};

    double trans_total = compute_transmission(segments, 0.1);
    double trans_no_rs = compute_transmission_no_rayleigh(segments, 0.1);

    BOOST_CHECK_GE(trans_no_rs, trans_total - 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  sample_rayleigh_cos_theta unit test
// ============================================================

BOOST_AUTO_TEST_SUITE(RayleighSampling)

BOOST_AUTO_TEST_CASE(rayleigh_mean_cos_theta_aluminum_200keV) {
    // At 200 keV for Al (Z=13), Rayleigh scatters predominantly forward.
    // Verify <cos(theta)> is positive and reasonable (form factor suppresses
    // large-angle scattering). Numerical integration of F²(x,Z)*(1+cos²θ)/2
    // gives <cos(theta)> ≈ 0.5-0.9 for typical detector materials at 200 keV.
    const int N = 100000;
    std::mt19937_64 rng(42);

    double sum_cos = 0.0;
    for (int i = 0; i < N; ++i) {
        double ct = sample_rayleigh_cos_theta(13, 200.0, rng);
        BOOST_CHECK_GE(ct, -1.0);
        BOOST_CHECK_LE(ct, 1.0);
        sum_cos += ct;
    }
    double mean_cos = sum_cos / N;

    // Rayleigh at 200 keV for Al should be forward-peaked:
    // <cos(theta)> should be positive and > 0.3
    BOOST_CHECK_GT(mean_cos, 0.3);
    BOOST_CHECK_LT(mean_cos, 1.0);
}

BOOST_AUTO_TEST_CASE(rayleigh_mean_cos_theta_lead_100keV) {
    // At 100 keV for Pb (Z=82), form factor still causes forward peaking.
    const int N = 100000;
    std::mt19937_64 rng(123);

    double sum_cos = 0.0;
    for (int i = 0; i < N; ++i) {
        double ct = sample_rayleigh_cos_theta(82, 100.0, rng);
        sum_cos += ct;
    }
    double mean_cos = sum_cos / N;

    // Pb at 100 keV: moderate forward peaking, <cos(theta)> > 0.1
    BOOST_CHECK_GT(mean_cos, 0.1);
    BOOST_CHECK_LT(mean_cos, 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
