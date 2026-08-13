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

#define BOOST_TEST_MODULE ExtendedSourceTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>
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
//  Cylindrical Source Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(CylindricalSource)

BOOST_AUTO_TEST_CASE(basic_invariants) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 2.0, 5.0);

    for (double E : {100.0, 662.0, 1332.0}) {
        auto res = calc.compute(E, 10000, 1);
        BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
        BOOST_CHECK_GE(res.full_energy_peak_efficiency, 0.0);
        BOOST_CHECK_LE(res.total_efficiency, 1.0);
        BOOST_CHECK_GT(res.total_efficiency, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(efficiency_bracketed_by_near_and_far_endpoints) {
    Material nai = make_NaI();

    EfficiencyCalculator calc_near;
    calc_near.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_near.set_point_source(Eigen::Vector3d(0.0, 0.0, -8.0));
    auto res_near = calc_near.compute(precision_config(662.0));

    EfficiencyCalculator calc_far;
    calc_far.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_far.set_point_source(Eigen::Vector3d(0.0, 0.0, -12.0));
    auto res_far = calc_far.compute(precision_config(662.0));

    EfficiencyCalculator calc_cyl;
    calc_cyl.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_cyl.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 0.1, 2.0);
    auto res_cyl = calc_cyl.compute(precision_config(662.0));

    double margin = 0.15;

    BOOST_CHECK_LE(res_cyl.total_efficiency,
                   res_near.total_efficiency * (1.0 + margin));
    BOOST_CHECK_GE(res_cyl.total_efficiency,
                   res_far.total_efficiency * (1.0 - margin));
}

BOOST_AUTO_TEST_CASE(point_like_cylinder_matches_point_source) {
    Material nai = make_NaI();
    Eigen::Vector3d center(0.0, 0.0, -10.0);

    EfficiencyCalculator calc_pt;
    calc_pt.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_pt.set_point_source(center);
    auto res_pt = calc_pt.compute(precision_config(662.0));

    EfficiencyCalculator calc_cyl;
    calc_cyl.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_cyl.set_cylindrical_source(center, 0.001, 0.001);
    auto res_cyl = calc_cyl.compute(precision_config(662.0));

    // With 0.5% precision on each, 5% tolerance is very conservative
    BOOST_CHECK_CLOSE(res_cyl.total_efficiency, res_pt.total_efficiency, 5.0);
}

BOOST_AUTO_TEST_CASE(tilted_source_basic_invariants) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    Eigen::Matrix3d rot;
    rot << 1, 0, 0,
           0, 0, -1,
           0, 1,  0;
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 1.0, 2.0, rot);

    auto res = calc.compute(662.0, 10000, 1);
    BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
    BOOST_CHECK_GE(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LE(res.total_efficiency, 1.0);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Rectangular Source Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(RectangularSource)

BOOST_AUTO_TEST_CASE(basic_invariants) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -10.0),
                                 Eigen::Vector3d(2.0, 2.0, 2.0));

    for (double E : {100.0, 662.0, 1332.0}) {
        auto res = calc.compute(E, 10000, 1);
        BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
        BOOST_CHECK_GE(res.full_energy_peak_efficiency, 0.0);
        BOOST_CHECK_LE(res.total_efficiency, 1.0);
        BOOST_CHECK_GT(res.total_efficiency, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(efficiency_bracketed_by_near_and_far_endpoints) {
    Material nai = make_NaI();

    EfficiencyCalculator calc_near;
    calc_near.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_near.set_point_source(Eigen::Vector3d(0.0, 0.0, -8.0));
    auto res_near = calc_near.compute(precision_config(662.0));

    EfficiencyCalculator calc_far;
    calc_far.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_far.set_point_source(Eigen::Vector3d(0.0, 0.0, -12.0));
    auto res_far = calc_far.compute(precision_config(662.0));

    EfficiencyCalculator calc_rect;
    calc_rect.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_rect.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -10.0),
                                      Eigen::Vector3d(0.05, 0.05, 2.0));
    auto res_rect = calc_rect.compute(precision_config(662.0));

    double margin = 0.15;

    BOOST_CHECK_LE(res_rect.total_efficiency,
                   res_near.total_efficiency * (1.0 + margin));
    BOOST_CHECK_GE(res_rect.total_efficiency,
                   res_far.total_efficiency * (1.0 - margin));
}

BOOST_AUTO_TEST_CASE(point_like_box_matches_point_source) {
    Material nai = make_NaI();
    Eigen::Vector3d center(0.0, 0.0, -10.0);

    EfficiencyCalculator calc_pt;
    calc_pt.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_pt.set_point_source(center);
    auto res_pt = calc_pt.compute(precision_config(662.0));

    EfficiencyCalculator calc_rect;
    calc_rect.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_rect.set_rectangular_source(center, Eigen::Vector3d(0.001, 0.001, 0.001));
    auto res_rect = calc_rect.compute(precision_config(662.0));

    BOOST_CHECK_CLOSE(res_rect.total_efficiency, res_pt.total_efficiency, 5.0);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Extended Source Consistency Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(ExtendedSourceConsistency)

BOOST_AUTO_TEST_CASE(cylindrical_source_fep_leq_total) {
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -8.0), 3.0, 3.0);

    for (double E : {60.0, 200.0, 662.0, 1332.0, 2000.0}) {
        auto res = calc.compute(E, 10000, 1);
        BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
    }
}

BOOST_AUTO_TEST_CASE(extended_source_efficiency_decreases_with_distance) {
    Material nai = make_NaI();

    EfficiencyCalculator calc_close;
    calc_close.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_close.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -5.0), 1.0, 1.0);
    auto res_close = calc_close.compute(precision_config(662.0));

    EfficiencyCalculator calc_far;
    calc_far.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_far.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -25.0), 1.0, 1.0);
    auto res_far = calc_far.compute(precision_config(662.0));

    BOOST_CHECK_GT(res_close.total_efficiency, res_far.total_efficiency);
}

BOOST_AUTO_TEST_CASE(rectangular_source_and_cylindrical_source_consistent) {
    Material nai = make_NaI();

    EfficiencyCalculator calc_cyl;
    calc_cyl.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_cyl.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 2.0, 2.0);
    auto res_cyl = calc_cyl.compute(precision_config(662.0));

    EfficiencyCalculator calc_rect;
    calc_rect.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_rect.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -10.0),
                                      Eigen::Vector3d(2.0, 2.0, 2.0));
    auto res_rect = calc_rect.compute(precision_config(662.0));

    BOOST_CHECK_GT(res_cyl.total_efficiency, 0.0);
    BOOST_CHECK_GT(res_rect.total_efficiency, 0.0);
    double ratio = res_cyl.total_efficiency / res_rect.total_efficiency;
    BOOST_CHECK_GT(ratio, 0.4);
    BOOST_CHECK_LT(ratio, 2.5);
}

BOOST_AUTO_TEST_CASE(hpge_with_extended_source) {
    Material hpge = make_HPGe();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &hpge, {4.0, 7.0});
    calc.set_bore_hole(0.5, 6.0);
    calc.set_dead_layer(0.1, 0.1);
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 3.0, 5.0);

    auto res = calc.compute(662.0, 10000, 1);
    BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
    BOOST_CHECK_GE(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LE(res.total_efficiency, 1.0);
    BOOST_CHECK_GT(res.total_efficiency, 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
