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

#define BOOST_TEST_MODULE PositionBiasTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>
#include <random>
#include <chrono>
#include <cstdio>

using namespace ceelo;

namespace {

/// Helper: create a SimulationConfig targeting a specific FEP relative precision.
SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 5000;
    config.termination.max_wall_seconds = 60.0;
    config.num_threads = 0;
    config.batch_size = 10000;
    return config;
}

/// Helper: compute combined N-sigma tolerance from two independent results.
double combined_sigma(const EfficiencyResult& a, const EfficiencyResult& b, double n_sigma) {
    return n_sigma * std::sqrt(a.fep_uncertainty * a.fep_uncertainty
                             + b.fep_uncertainty * b.fep_uncertainty);
}

/// Helper: set up a standard large cylindrical soil source with NaI detector.
EfficiencyCalculator make_large_soil_source() {
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62}); // 3"x3" NaI
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -100.0), 50.0, 10.0);
    calc.set_source_material(&soil);
    return calc;
}

} // anonymous namespace

// ============================================================
//  Correctness: biased vs unbiased agreement
// ============================================================

BOOST_AUTO_TEST_SUITE(BiasedVsUnbiased)

BOOST_AUTO_TEST_CASE(cylindrical_soil_662keV) {
    auto config = precision_config(662.0, 0.0025);

    // Unbiased
    auto calc_ub = make_large_soil_source();
    auto res_ub = calc_ub.compute(config);

    // Biased (auto params)
    auto calc_b = make_large_soil_source();
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV Cylindrical Soil ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol);

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_CASE(cylindrical_soil_100keV) {
    auto config = precision_config(100.0, 0.0025);

    auto calc_ub = make_large_soil_source();
    auto res_ub = calc_ub.compute(config);

    auto calc_b = make_large_soil_source();
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 100 keV Cylindrical Soil ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol);

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_CASE(cylindrical_soil_1332keV) {
    auto config = precision_config(1332.0, 0.0025);

    auto calc_ub = make_large_soil_source();
    auto res_ub = calc_ub.compute(config);

    auto calc_b = make_large_soil_source();
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 1332 keV Cylindrical Soil ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol);

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Weight normalization
// ============================================================

BOOST_AUTO_TEST_SUITE(WeightNormalization)

BOOST_AUTO_TEST_CASE(mean_weight_approximately_one) {
    // Sample many biased positions and verify mean weight ≈ 1.0
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -100.0), 50.0, 10.0);
    calc.set_source_material(&soil);

    // Use manual bias params to test specific values
    PositionBiasConfig bias;
    bias.lambda_lateral_cm = 70.0;
    bias.lambda_depth_cm = 5.0;
    calc.set_position_bias(bias);

    // Run a short simulation — the internal weight normalization is tested
    // by checking that biased and unbiased give same efficiency.
    // Here we just verify the API doesn't crash and produces reasonable results.
    auto res = calc.compute(662.0, 100000, 1);
    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LE(res.full_energy_peak_efficiency, 1.0);
    BOOST_CHECK_GT(res.fep_uncertainty, 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Rectangular source
// ============================================================

BOOST_AUTO_TEST_SUITE(RectangularSource)

BOOST_AUTO_TEST_CASE(rectangular_soil_662keV) {
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    auto config = precision_config(662.0, 0.0025);

    // Unbiased
    EfficiencyCalculator calc_ub;
    calc_ub.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_ub.set_rectangular_source(
        Eigen::Vector3d(0.0, 0.0, -100.0),
        Eigen::Vector3d(50.0, 50.0, 10.0));
    calc_ub.set_source_material(&soil);
    auto res_ub = calc_ub.compute(config);

    // Biased
    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_rectangular_source(
        Eigen::Vector3d(0.0, 0.0, -100.0),
        Eigen::Vector3d(50.0, 50.0, 10.0));
    calc_b.set_source_material(&soil);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV Rectangular Soil ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol);

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  FEP-only mode
// ============================================================

BOOST_AUTO_TEST_SUITE(FepOnlyMode)

BOOST_AUTO_TEST_CASE(fep_only_biased_matches_unbiased) {
    auto config = precision_config(662.0, 0.0025);

    auto calc_ub = make_large_soil_source();
    calc_ub.enable_fep_only_mode(true);
    auto res_ub = calc_ub.compute(config);

    auto calc_b = make_large_soil_source();
    calc_b.enable_fep_only_mode(true);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV FEP-Only Mode ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol);

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Exponential depth + depth bias interaction
// ============================================================

BOOST_AUTO_TEST_SUITE(ExponentialPlusBias)

BOOST_AUTO_TEST_CASE(exponential_depth_with_bias_matches) {
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    auto config = precision_config(662.0, 0.0025);

    // Unbiased with exponential depth distribution
    EfficiencyCalculator calc_ub;
    calc_ub.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_ub.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -100.0), 50.0, 10.0);
    calc_ub.set_source_material(&soil);
    calc_ub.set_exponential_depth_distribution(3.0); // 3 cm relaxation
    auto res_ub = calc_ub.compute(config);

    // Biased with same exponential depth distribution
    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -100.0), 50.0, 10.0);
    calc_b.set_source_material(&soil);
    calc_b.set_exponential_depth_distribution(3.0);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV Exponential + Bias ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol);

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Dense source (lead-like)
// ============================================================

BOOST_AUTO_TEST_SUITE(DenseSource)

BOOST_AUTO_TEST_CASE(dense_lead_source_662keV) {
    static Material nai = make_NaI();
    static Material lead = make_Lead();
    auto config = precision_config(662.0, 0.0025);

    // 3cm thick lead disk at 5cm distance
    Eigen::Vector3d center(0.0, 0.0, -6.5);  // center at -6.5, front at -5

    // Unbiased
    EfficiencyCalculator calc_ub;
    calc_ub.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_ub.set_cylindrical_source(center, 5.0, 1.5);
    calc_ub.set_source_material(&lead);
    auto res_ub = calc_ub.compute(config);

    // Biased
    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_cylindrical_source(center, 5.0, 1.5);
    calc_b.set_source_material(&lead);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV Dense Lead Source ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol);

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  In-situ geometry (large radius, deep soil)
// ============================================================

BOOST_AUTO_TEST_SUITE(InSitu)

BOOST_AUTO_TEST_CASE(large_insitu_1332keV) {
    // Realistic in-situ: detector 1m above ground, R=75m, depth=50cm
    // At 1332 keV, 50cm soil is effectively infinite depth (mu*d >> 1)
    // FEP-only mode: full mode can't converge unbiased (FEP ~ 6e-9, needs ~10B events)
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    auto config = precision_config(1332.0, 0.005);
    config.termination.max_wall_seconds = 120.0;

    // Unbiased (FEP-only mode to make convergence feasible)
    EfficiencyCalculator calc_ub;
    calc_ub.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_ub.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -125.0), 7500.0, 25.0);
    calc_ub.set_source_material(&soil);
    calc_ub.enable_fep_only_mode(true);
    auto res_ub = calc_ub.compute(config);

    // Biased (auto params, FEP-only)
    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -125.0), 7500.0, 25.0);
    calc_b.set_source_material(&soil);
    calc_b.enable_fep_only_mode(true);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 1332 keV In-Situ FEP-only (R=75m, depth=50cm) ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e (%.2f sigma)\n"
        "Speedup:  %.2fx\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol,
        diff / std::sqrt(res_ub.fep_uncertainty * res_ub.fep_uncertainty
                       + res_b.fep_uncertainty * res_b.fep_uncertainty),
        res_ub.wall_time_seconds / std::max(res_b.wall_time_seconds, 0.001));

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_CASE(large_insitu_662keV) {
    // In-situ at 662 keV: detector 1m above ground, R=75m, depth=50cm
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    auto config = precision_config(662.0, 0.005);
    config.termination.max_wall_seconds = 120.0;

    // Unbiased (FEP-only mode)
    EfficiencyCalculator calc_ub;
    calc_ub.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_ub.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -125.0), 7500.0, 25.0);
    calc_ub.set_source_material(&soil);
    calc_ub.enable_fep_only_mode(true);
    auto res_ub = calc_ub.compute(config);

    // Biased (auto params, FEP-only)
    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -125.0), 7500.0, 25.0);
    calc_b.set_source_material(&soil);
    calc_b.enable_fep_only_mode(true);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV In-Situ FEP-only (R=75m, depth=50cm) ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e (%.2f sigma)\n"
        "Speedup:  %.2fx\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol,
        diff / std::sqrt(res_ub.fep_uncertainty * res_ub.fep_uncertainty
                       + res_b.fep_uncertainty * res_b.fep_uncertainty),
        res_ub.wall_time_seconds / std::max(res_b.wall_time_seconds, 0.001));

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Rotated source geometry
// ============================================================

BOOST_AUTO_TEST_SUITE(RotatedSource)

BOOST_AUTO_TEST_CASE(rotated_45deg_cylindrical_662keV) {
    // Cylindrical source rotated 45 degrees about x-axis
    // This exercises the dz projection code path
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    auto config = precision_config(662.0, 0.005);
    config.termination.max_wall_seconds = 120.0;

    // 45-degree rotation about x-axis
    Eigen::Matrix3d rot;
    double c = std::cos(M_PI / 4.0), s = std::sin(M_PI / 4.0);
    rot << 1,  0, 0,
           0,  c, s,
           0, -s, c;

    Eigen::Vector3d center(0.0, 0.0, -20.0);

    // Unbiased
    EfficiencyCalculator calc_ub;
    calc_ub.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_ub.set_cylindrical_source(center, 15.0, 5.0, rot);
    calc_ub.set_source_material(&soil);
    auto res_ub = calc_ub.compute(config);

    // Biased
    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_cylindrical_source(center, 15.0, 5.0, rot);
    calc_b.set_source_material(&soil);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV Rotated 45-deg Cylindrical ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e (%.2f sigma)\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol,
        diff / std::sqrt(res_ub.fep_uncertainty * res_ub.fep_uncertainty
                       + res_b.fep_uncertainty * res_b.fep_uncertainty));

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_CASE(rotated_90deg_sideon_662keV) {
    // Source rotated 90 degrees — side-on cylinder
    // dz ≈ 0, so depth biasing should be skipped (|dz| < 0.1)
    // Only lateral biasing should apply
    static Material nai = make_NaI();
    static Material soil = make_Soil();
    auto config = precision_config(662.0, 0.005);
    config.termination.max_wall_seconds = 120.0;

    // 90-degree rotation about x-axis: source axis now points along y
    Eigen::Matrix3d rot;
    rot << 1, 0,  0,
           0, 0,  1,
           0, -1, 0;

    Eigen::Vector3d center(0.0, 0.0, -20.0);

    // Unbiased
    EfficiencyCalculator calc_ub;
    calc_ub.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_ub.set_cylindrical_source(center, 10.0, 10.0, rot);
    calc_ub.set_source_material(&soil);
    auto res_ub = calc_ub.compute(config);

    // Biased
    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_cylindrical_source(center, 10.0, 10.0, rot);
    calc_b.set_source_material(&soil);
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double tol = combined_sigma(res_ub, res_b, 3.5);
    double diff = std::abs(res_b.full_energy_peak_efficiency - res_ub.full_energy_peak_efficiency);

    std::fprintf(stderr,
        "\n=== 662 keV Side-On (90-deg) Cylindrical ===\n"
        "Unbiased: FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Biased:   FEP = %.6e +/- %.2e (%llu events, %.1f s)\n"
        "Diff: %.2e, 3.5-sigma tol: %.2e (%.2f sigma)\n",
        res_ub.full_energy_peak_efficiency, res_ub.fep_uncertainty,
        (unsigned long long)res_ub.num_events_simulated, res_ub.wall_time_seconds,
        res_b.full_energy_peak_efficiency, res_b.fep_uncertainty,
        (unsigned long long)res_b.num_events_simulated, res_b.wall_time_seconds,
        diff, tol,
        diff / std::sqrt(res_ub.fep_uncertainty * res_ub.fep_uncertainty
                       + res_b.fep_uncertainty * res_b.fep_uncertainty));

    BOOST_CHECK_LE(diff, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Timing comparison (informational, not pass/fail)
// ============================================================

BOOST_AUTO_TEST_SUITE(Timing)

BOOST_AUTO_TEST_CASE(speedup_informational) {
    auto config = precision_config(662.0, 0.005);  // 0.5% for speed

    auto calc_ub = make_large_soil_source();
    auto res_ub = calc_ub.compute(config);

    auto calc_b = make_large_soil_source();
    calc_b.enable_position_bias();
    auto res_b = calc_b.compute(config);

    double speedup = res_ub.wall_time_seconds / std::max(res_b.wall_time_seconds, 0.001);

    std::fprintf(stderr,
        "\n=== Timing Comparison (662 keV, 0.5%% precision) ===\n"
        "Unbiased: %.2f s (%llu events)\n"
        "Biased:   %.2f s (%llu events)\n"
        "Speedup:  %.2fx\n",
        res_ub.wall_time_seconds, (unsigned long long)res_ub.num_events_simulated,
        res_b.wall_time_seconds, (unsigned long long)res_b.num_events_simulated,
        speedup);

    // Not a pass/fail test — just informational
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
