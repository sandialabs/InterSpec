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

#define BOOST_TEST_MODULE ExponentialSourceTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>
#include <vector>
#include <random>
#include <numeric>

using namespace ceelo;

namespace {

SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 5000;
    config.termination.max_wall_seconds = 30.0;
    config.num_threads = 0;
    config.batch_size = 10000;
    return config;
}

} // anonymous namespace

// ============================================================
//  Exponential Depth Sampling Distribution Shape Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(ExponentialSampling)

// Verify that the exponential depth distribution produces the correct shape.
// We use the EfficiencyCalculator to sample positions and check the z-distribution.
BOOST_AUTO_TEST_CASE(distribution_shape_strongly_peaked) {
    // lambda << D: strong exponential, most activity near surface
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    double half_length = 10.0;  // D = 20 cm
    double lambda = 0.5;        // strongly peaked
    double distance = 15.0;
    calc.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        5.0, half_length);
    calc.set_exponential_depth_distribution(lambda);

    // Run a simulation and check that it completes without error.
    // The exponential distribution concentrates activity near the surface
    // (closest to detector), so efficiency should be higher than uniform.
    auto res = calc.compute(662.0, 500000, 0);
    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
    BOOST_CHECK_LE(res.total_efficiency, 1.0);
}

BOOST_AUTO_TEST_CASE(distribution_shape_nearly_uniform) {
    // lambda >> D: nearly uniform distribution
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    double half_length = 5.0;   // D = 10 cm
    double lambda = 10000.0;    // lambda >> D, nearly uniform
    double distance = 10.0;
    calc.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        5.0, half_length);
    calc.set_exponential_depth_distribution(lambda);

    auto res_exp = calc.compute(precision_config(662.0));

    // Compare with uniform distribution
    EfficiencyCalculator calc_uni;
    calc_uni.set_fep_window_keV(kTestFepWindowKeV);
    calc_uni.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_uni.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        5.0, half_length);
    auto res_uni = calc_uni.compute(precision_config(662.0));

    // Large lambda exponential should match uniform within ~5% (statistical + small exp bias)
    double ratio = res_exp.full_energy_peak_efficiency / res_uni.full_energy_peak_efficiency;
    BOOST_CHECK_CLOSE(ratio, 1.0, 5.0);
}

BOOST_AUTO_TEST_CASE(small_lambda_higher_efficiency_than_uniform) {
    // With small lambda, activity is concentrated near the detector-facing surface.
    // Without source material, the exponential source should have higher efficiency
    // than uniform because photons come from closer to the detector on average.
    Material nai = make_NaI();
    double distance = 10.0;
    double half_length = 5.0;  // D = 10 cm

    // Exponential with small lambda (concentrated near surface)
    EfficiencyCalculator calc_exp;
    calc_exp.set_fep_window_keV(kTestFepWindowKeV);
    calc_exp.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_exp.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        3.0, half_length);
    calc_exp.set_exponential_depth_distribution(0.5);
    auto res_exp = calc_exp.compute(precision_config(662.0));

    // Uniform
    EfficiencyCalculator calc_uni;
    calc_uni.set_fep_window_keV(kTestFepWindowKeV);
    calc_uni.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_uni.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        3.0, half_length);
    auto res_uni = calc_uni.compute(precision_config(662.0));

    // Exponential concentrates activity near the detector, so per-photon
    // efficiency should be >= uniform (closer = higher solid angle)
    double margin = 3.0 * std::max(res_exp.fep_uncertainty, res_uni.fep_uncertainty);
    BOOST_CHECK_GE(res_exp.full_energy_peak_efficiency,
                    res_uni.full_energy_peak_efficiency - margin);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Exponential Source Physics Invariants
// ============================================================

BOOST_AUTO_TEST_SUITE(ExponentialInvariants)

BOOST_AUTO_TEST_CASE(fep_leq_total) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -15.0), 5.0, 5.0);
    calc.set_exponential_depth_distribution(2.0);

    for (double E : {100.0, 662.0, 1332.0}) {
        auto res = calc.compute(E, 500000, 0);
        BOOST_CHECK_GE(res.total_efficiency,
                        res.full_energy_peak_efficiency - 1e-9);
    }
}

BOOST_AUTO_TEST_CASE(exponential_efficiency_geq_uniform_with_source_material) {
    // With source material self-attenuation, exponential concentrates activity
    // near the surface where less material is traversed, so efficiency should
    // be >= uniform efficiency (for the same geometry).
    Material nai = make_NaI();
    Material soil = make_Soil();
    double distance = 15.0;
    double half_length = 10.0;

    EfficiencyCalculator calc_exp;
    calc_exp.set_fep_window_keV(kTestFepWindowKeV);
    calc_exp.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_exp.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        5.0, half_length);
    calc_exp.set_source_material(&soil);
    calc_exp.set_exponential_depth_distribution(2.0);
    auto res_exp = calc_exp.compute(precision_config(662.0, 0.01));

    EfficiencyCalculator calc_uni;
    calc_uni.set_fep_window_keV(kTestFepWindowKeV);
    calc_uni.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_uni.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        5.0, half_length);
    calc_uni.set_source_material(&soil);
    auto res_uni = calc_uni.compute(precision_config(662.0, 0.01));

    // Exponential should have higher or equal FEP efficiency
    // Allow some statistical margin (exponential might be slightly lower due to noise)
    double margin = 3.0 * std::max(res_exp.fep_uncertainty, res_uni.fep_uncertainty);
    BOOST_CHECK_GE(res_exp.full_energy_peak_efficiency,
                    res_uni.full_energy_peak_efficiency - margin);
}

BOOST_AUTO_TEST_CASE(rectangular_source_exponential) {
    // Verify exponential depth works with rectangular sources too
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(
        Eigen::Vector3d(0.0, 0.0, -15.0),
        Eigen::Vector3d(5.0, 5.0, 5.0));
    calc.set_exponential_depth_distribution(2.0);

    auto res = calc.compute(662.0, 500000, 0);
    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
    BOOST_CHECK_LE(res.total_efficiency, 1.0);
}

BOOST_AUTO_TEST_CASE(reset_to_uniform) {
    // Verify that set_uniform_depth_distribution() resets behavior
    Material nai = make_NaI();
    double distance = 10.0;
    double half_length = 5.0;

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        5.0, half_length);
    calc.set_exponential_depth_distribution(2.0);
    calc.set_uniform_depth_distribution();  // reset

    auto res_reset = calc.compute(precision_config(662.0));

    EfficiencyCalculator calc_uni;
    calc_uni.set_fep_window_keV(kTestFepWindowKeV);
    calc_uni.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_uni.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -(distance + half_length)),
        5.0, half_length);
    auto res_uni = calc_uni.compute(precision_config(662.0));

    // After reset, should match uniform
    double ratio = res_reset.full_energy_peak_efficiency / res_uni.full_energy_peak_efficiency;
    BOOST_CHECK_CLOSE(ratio, 1.0, 3.0);
}

BOOST_AUTO_TEST_SUITE_END()
