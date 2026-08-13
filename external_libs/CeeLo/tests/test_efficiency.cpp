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

#define BOOST_TEST_MODULE EfficiencyTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>
#include <map>
#include <vector>

using namespace ceelo;

namespace {

/// Helper: create a SimulationConfig targeting 0.5% FEP relative precision.
SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 5000;
    config.termination.max_wall_seconds = 30.0;
    config.num_threads = 0;  // auto-detect
    config.batch_size = 10000;
    return config;
}

} // anonymous namespace

// ============================================================
//  Invariant Tests
//  These must hold for any valid simulation regardless of geometry
// ============================================================

BOOST_AUTO_TEST_SUITE(Invariants)

BOOST_AUTO_TEST_CASE(fep_leq_total_efficiency) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    for (double E : {100.0, 662.0, 1332.0, 2000.0}) {
        auto res = calc.compute(E, 5000, 1);
        BOOST_CHECK_GE(res.total_efficiency,
                       res.full_energy_peak_efficiency - 1e-9);
    }
}

BOOST_AUTO_TEST_CASE(efficiencies_non_negative) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    auto res = calc.compute(662.0, 5000, 1);
    BOOST_CHECK_GE(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_GE(res.total_efficiency, 0.0);
    BOOST_CHECK_GE(res.fep_uncertainty, 0.0);
    BOOST_CHECK_GE(res.total_uncertainty, 0.0);
}

BOOST_AUTO_TEST_CASE(efficiencies_leq_one) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    for (double E : {100.0, 662.0, 2000.0}) {
        auto res = calc.compute(E, 5000, 1);
        BOOST_CHECK_LE(res.full_energy_peak_efficiency, 1.0);
        BOOST_CHECK_LE(res.total_efficiency, 1.0);
    }
}

BOOST_AUTO_TEST_CASE(event_count_reported_correctly) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    const uint64_t N = 7777;
    auto res = calc.compute(662.0, N, 1);
    BOOST_CHECK_EQUAL(res.num_events_simulated, N);
}

BOOST_AUTO_TEST_CASE(fep_leq_total_hpge_with_bore) {
    Material hpge = make_HPGe();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &hpge, {5.0, 8.0});
    calc.set_bore_hole(0.5, 7.0);
    calc.set_dead_layer(0.03, 0.03);
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -15.0));

    auto res = calc.compute(662.0, 5000, 1);
    BOOST_CHECK_GE(res.total_efficiency,
                   res.full_energy_peak_efficiency - 1e-9);
    BOOST_CHECK_GE(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LE(res.total_efficiency, 1.0);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Physical Plausibility Tests
//  Verify qualitative physics trends with the MC simulator.
// ============================================================

BOOST_AUTO_TEST_SUITE(PhysicalPlausibility)

BOOST_AUTO_TEST_CASE(nonzero_efficiency_at_662keV) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    auto res = calc.compute(precision_config(662.0));

    BOOST_CHECK_GT(res.total_efficiency, 1e-5);
    BOOST_CHECK_LT(res.total_efficiency, 1.0);
    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
}

BOOST_AUTO_TEST_CASE(total_efficiency_decreases_with_energy) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    auto res_low  = calc.compute(precision_config(300.0));
    auto res_high = calc.compute(precision_config(3000.0));

    BOOST_CHECK_GT(res_low.total_efficiency, res_high.total_efficiency);
}

BOOST_AUTO_TEST_CASE(peak_to_total_ratio_higher_at_low_energy) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    auto res_low  = calc.compute(precision_config(100.0));
    auto res_high = calc.compute(precision_config(2000.0));

    double ratio_low  = res_low.full_energy_peak_efficiency  /
                        std::max(res_low.total_efficiency,  1e-10);
    double ratio_high = res_high.full_energy_peak_efficiency /
                        std::max(res_high.total_efficiency, 1e-10);

    BOOST_CHECK_GT(ratio_low, ratio_high);
}

BOOST_AUTO_TEST_CASE(dead_layer_reduces_efficiency) {
    Material nai = make_NaI();

    EfficiencyCalculator calc_nodl;
    calc_nodl.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_nodl.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    auto res_nodl = calc_nodl.compute(precision_config(100.0));

    EfficiencyCalculator calc_dl;
    calc_dl.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_dl.set_dead_layer(2.0, 0.0, 0.0);
    calc_dl.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    auto res_dl = calc_dl.compute(precision_config(100.0));

    BOOST_CHECK_LT(res_dl.total_efficiency, res_nodl.total_efficiency * 0.5);
}

BOOST_AUTO_TEST_CASE(batch_returns_correct_size) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    std::vector<double> energies = {100.0, 300.0, 662.0, 1000.0, 2000.0};
    auto results = calc.compute_batch(energies, 2000, 1);

    BOOST_REQUIRE_EQUAL(results.size(), energies.size());
    for (size_t i = 0; i < results.size(); ++i) {
        BOOST_CHECK_EQUAL(results[i].num_events_simulated, 2000u);
        BOOST_CHECK_GE(results[i].total_efficiency, 0.0);
        BOOST_CHECK_GE(results[i].total_efficiency,
                       results[i].full_energy_peak_efficiency - 1e-9);
    }
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Variance Reduction / Biasing Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(Biasing)

BOOST_AUTO_TEST_CASE(cone_sampling_consistent_with_isotropic) {
    Material nai = make_NaI();

    Eigen::Vector3d src(0.0, 0.0, -5.0);

    EfficiencyCalculator calc_cone;
    calc_cone.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_cone.set_point_source(src);
    calc_cone.enable_cone_sampling(true);
    auto res_cone = calc_cone.compute(precision_config(662.0));

    EfficiencyCalculator calc_iso;
    calc_iso.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_iso.set_point_source(src);
    calc_iso.enable_cone_sampling(false);
    auto res_iso = calc_iso.compute(precision_config(662.0));

    double sigma_combined = std::sqrt(
        res_cone.fep_uncertainty * res_cone.fep_uncertainty
      + res_iso.fep_uncertainty  * res_iso.fep_uncertainty);

    double diff = std::abs(res_cone.full_energy_peak_efficiency
                         - res_iso.full_energy_peak_efficiency);
    BOOST_CHECK_LT(diff, 5.0 * sigma_combined + 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Cascade Correction Tests (deterministic)
// ============================================================

BOOST_AUTO_TEST_SUITE(CascadeCorrection)

BOOST_AUTO_TEST_CASE(no_coincident_gammas_unchanged) {
    auto result = EfficiencyCalculator::cascade_correction(
        0.5, 0.8, {}, {});
    BOOST_CHECK_CLOSE(result.summing_out_factor, 1.0, 1e-9);
    BOOST_CHECK_CLOSE(result.corrected_fep,      0.5, 1e-9);
}

BOOST_AUTO_TEST_CASE(cascade_correction_single_gamma) {
    EfficiencyResult ej{};
    ej.total_efficiency = 0.1;

    std::map<double, EfficiencyResult> cache;
    cache[1173.2] = ej;

    auto result = EfficiencyCalculator::cascade_correction(
        0.5, 0.8, {{1173.2, 1.0}}, cache);

    BOOST_CHECK_CLOSE(result.summing_out_factor, 0.9,       1e-9);
    BOOST_CHECK_CLOSE(result.corrected_fep,      0.5 * 0.9, 1e-9);
}

BOOST_AUTO_TEST_CASE(cascade_correction_two_gammas) {
    EfficiencyResult e1{}, e2{};
    e1.total_efficiency = 0.2;
    e2.total_efficiency = 0.1;

    std::map<double, EfficiencyResult> cache;
    cache[1173.2] = e1;
    cache[1332.5] = e2;

    auto result = EfficiencyCalculator::cascade_correction(
        0.4, 0.8, {{1173.2, 0.5}, {1332.5, 1.0}}, cache);

    BOOST_CHECK_CLOSE(result.summing_out_factor, 0.81,       1e-9);
    BOOST_CHECK_CLOSE(result.corrected_fep,      0.4 * 0.81, 1e-9);
}

BOOST_AUTO_TEST_CASE(cascade_summing_out_leq_one) {
    EfficiencyResult ej{};
    ej.total_efficiency = 0.3;

    std::map<double, EfficiencyResult> cache;
    cache[511.0] = ej;

    auto result = EfficiencyCalculator::cascade_correction(
        0.5, 0.8, {{511.0, 0.5}}, cache);

    BOOST_CHECK_LE(result.summing_out_factor, 1.0);
    BOOST_CHECK_LE(result.corrected_fep, 0.5);
}

BOOST_AUTO_TEST_CASE(cascade_missing_energy_in_cache_skipped) {
    std::map<double, EfficiencyResult> empty_cache;

    auto result = EfficiencyCalculator::cascade_correction(
        0.5, 0.8, {{1173.2, 1.0}}, empty_cache);

    BOOST_CHECK_CLOSE(result.summing_out_factor, 1.0, 1e-9);
    BOOST_CHECK_CLOSE(result.corrected_fep,      0.5, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Pulse-Height Distribution Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(PulseHeight)

BOOST_AUTO_TEST_CASE(bin_values_non_negative) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    std::vector<float> edges;
    for (int i = 0; i <= 10; ++i)
        edges.push_back(static_cast<float>(i * 80.0f));

    auto res = calc.compute(662.0, 10000, 1, edges);
    BOOST_REQUIRE_EQUAL(res.pulse_height_distribution.size(), 10u);

    for (float v : res.pulse_height_distribution)
        BOOST_CHECK_GE(v, 0.0f);
}

BOOST_AUTO_TEST_CASE(bin_sum_consistent_with_total_eff) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    std::vector<float> edges;
    for (int i = 0; i <= 20; ++i)
        edges.push_back(static_cast<float>(i * 40.0f));

    auto res = calc.compute(662.0, 10000, 1, edges);
    BOOST_REQUIRE_EQUAL(res.pulse_height_distribution.size(), 20u);

    double sum = 0.0;
    for (float v : res.pulse_height_distribution)
        sum += static_cast<double>(v);

    BOOST_CHECK_CLOSE(sum, res.total_efficiency, 5.0);
}

BOOST_AUTO_TEST_CASE(fep_bin_has_counts) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    std::vector<float> edges = {0.0f, 600.0f, 700.0f, 800.0f};
    auto res = calc.compute(662.0, 20000, 1, edges);
    BOOST_REQUIRE_EQUAL(res.pulse_height_distribution.size(), 3u);

    BOOST_CHECK_GT(res.pulse_height_distribution[1], 0.0f);
}

// ------------------------------------------------------------
//  Per-bin statistical uncertainty (review fix D1)
//  pulse_height_uncertainty[i] = sqrt(Sum w_i^2) / sum_weights,
//  parallel to pulse_height_distribution[i] = (Sum w_i) / sum_weights.
//  sum_weights == num_events_simulated (the IS estimator denominator).
// ------------------------------------------------------------

BOOST_AUTO_TEST_CASE(uncertainty_array_well_formed) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -25.0));

    std::vector<float> edges;
    for (int i = 0; i <= 10; ++i)
        edges.push_back(static_cast<float>(i * 80.0f));

    auto res = calc.compute(662.0, 10000, 1, edges);

    // Parallel to the distribution: same length, every entry finite and >= 0.
    BOOST_REQUIRE_EQUAL(res.pulse_height_uncertainty.size(),
                        res.pulse_height_distribution.size());
    for (float s : res.pulse_height_uncertainty) {
        BOOST_CHECK(std::isfinite(s));
        BOOST_CHECK_GE(s, 0.0f);
    }
}

BOOST_AUTO_TEST_CASE(analog_unweighted_recovers_poisson) {
    // Analog (default BiasingConfig) + isotropic emission (cone off) gives
    // every event weight exactly 1.0, so Sum w^2 == Sum w == count and the
    // per-bin sigma must reduce to Poisson sqrt(count)/N. This is an algebraic
    // identity (not a statistical estimate), so the tolerance can be tight.
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));  // close -> hit rate
    calc.enable_cone_sampling(false);                        // isotropic, w == 1
    calc.set_biasing(BiasingConfig{});                       // force analog

    std::vector<float> edges;
    for (int i = 0; i <= 7; ++i)
        edges.push_back(static_cast<float>(i * 100.0f));     // 100 keV bins

    auto res = calc.compute(662.0, 200000, 0, edges);
    BOOST_REQUIRE_EQUAL(res.pulse_height_uncertainty.size(),
                        res.pulse_height_distribution.size());

    const double N = static_cast<double>(res.num_events_simulated);
    BOOST_REQUIRE_GT(N, 0.0);

    int populated = 0;
    for (size_t i = 0; i < res.pulse_height_distribution.size(); ++i) {
        const double dist = res.pulse_height_distribution[i];
        const double sigma = res.pulse_height_uncertainty[i];
        if (dist <= 0.0) continue;          // skip empty bins
        ++populated;
        const double count = dist * N;       // weighted bin count == sqrt argument
        // sigma * N must recover sqrt(count) for unweighted sampling.
        const double expected = std::sqrt(count) / N;
        BOOST_CHECK_CLOSE(sigma, expected, 1.0);  // 1% (float-rounding only)
    }
    BOOST_REQUIRE_GT(populated, 0);
}

BOOST_AUTO_TEST_CASE(weighted_cone_differs_from_poisson) {
    // Cone biasing (analog otherwise) gives every event the SAME constant
    // weight omega_frac < 1. Then R_i = sigma_i^2 * N / dist_i == Sum w^2 / Sum w
    // == omega_frac for every populated bin: clearly < 1 (so sigma differs from
    // the naive Poisson sqrt(count)/N) and constant across bins (the weight
    // dispersion that the bare counts cannot capture).
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
    calc.enable_cone_sampling(true);     // cone -> constant weight omega_frac
    calc.set_biasing(BiasingConfig{});   // analog otherwise (no forcing)

    std::vector<float> edges;
    for (int i = 0; i <= 7; ++i)
        edges.push_back(static_cast<float>(i * 100.0f));

    auto res = calc.compute(662.0, 200000, 0, edges);
    BOOST_REQUIRE_EQUAL(res.pulse_height_uncertainty.size(),
                        res.pulse_height_distribution.size());

    const double N = static_cast<double>(res.num_events_simulated);
    BOOST_REQUIRE_GT(N, 0.0);

    double ref_R = -1.0;
    int populated = 0;
    for (size_t i = 0; i < res.pulse_height_distribution.size(); ++i) {
        const double dist = res.pulse_height_distribution[i];
        const double sigma = res.pulse_height_uncertainty[i];
        if (dist <= 0.0) continue;
        ++populated;
        const double R = sigma * sigma * N / dist;  // == omega_frac
        // Weighted: R must be a constant strictly below the analog value of 1.
        BOOST_CHECK_LT(R, 0.99);
        BOOST_CHECK_GT(R, 0.0);
        if (ref_R < 0.0) ref_R = R;
        else BOOST_CHECK_CLOSE(R, ref_R, 1.0);  // same weight across all bins
    }
    BOOST_REQUIRE_GT(populated, 1);  // need >= 2 bins to check constancy
}

BOOST_AUTO_TEST_SUITE_END()
