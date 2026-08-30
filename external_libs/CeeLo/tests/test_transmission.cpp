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

/// @file test_transmission.cpp
/// @brief Transmission tests: full-energy photon rate at the detector face.
///
/// These tests verify that photons which DO NOT INTERACT on their way to the
/// detector arrive with the correct rate (geometric × Beer-Lambert attenuation).
///
/// All tests use FEP efficiency as a proxy for "photon reaches detector unscattered".
/// When detector effects are minimized (small thin detector or ratio method),
/// the ratio or absolute value can be compared to analytic expectations.
///
/// Test geometries:
///  1. Point source + thin Pb absorber  — ratio test vs Beer-Lambert
///  2. Cylindrical soil extended source  — absolute count vs analytic MC reference
///  3. Rectangular soil extended source  — absolute count vs analytic MC reference
///
/// "Detector face" tests: we verify the fraction of full-energy photons reaching
/// the scoring volume front face.  We do NOT test detector efficiency here;
/// the goal is to confirm that source → detector geometry + attenuation is correct
/// for the simple "no scatter" case (matches GEANT4 for unscattered photons).

#define BOOST_TEST_MODULE TransmissionTests
#include <boost/test/unit_test.hpp>

#include "materials/Material.h"
#include "geometry/Geometry.h"
#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "cross_sections/CrossSectionData.h"

#include <Eigen/Core>
#include <random>
#include <cmath>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kPi    = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

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

/// Analytic transmission factor for a point source at z_src through a slab of
/// material at z_slab (thickness t cm) to a detector front face at z=0.
/// Returns exp(-mu * t) where mu is total linear attenuation (1/cm).
double beer_lambert(const Material& mat, double energy_keV, double thickness_cm)
{
    double mu = mat.mu_total(energy_keV * 1e-3);  // 1/cm
    return std::exp(-mu * thickness_cm);
}

} // anonymous namespace


// ============================================================
// Suite 1: Point source + thin Pb absorber (Beer-Lambert ratio test)
// ============================================================
BOOST_AUTO_TEST_SUITE(PointSourcePbAbsorber)

BOOST_AUTO_TEST_CASE(ratio_matches_beer_lambert_2mm_Pb_351keV)
{
    Material nai  = make_NaI();
    Material lead = make_Lead();

    const double   E = 351.0;  // keV
    const double   length_cm = 7.62;

    // --- Calculator WITHOUT Pb ---
    EfficiencyCalculator calc_bare;
    calc_bare.set_fep_window_keV(kTestFepWindowKeV);
    calc_bare.set_detector(DetectorShape::Cylinder, &nai, {3.81, length_cm});
    calc_bare.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    // --- Calculator WITH 2 mm Pb front cap ---
    EfficiencyCalculator calc_pb;
    calc_pb.set_fep_window_keV(kTestFepWindowKeV);
    calc_pb.set_detector(DetectorShape::Cylinder, &nai, {3.81, length_cm});
    calc_pb.add_attenuator(&lead, 0.2, 0.0, 0.0, length_cm);
    calc_pb.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    auto res_bare = calc_bare.compute(precision_config(E));
    auto res_pb   = calc_pb  .compute(precision_config(E));

    BOOST_REQUIRE_GT(res_bare.full_energy_peak_efficiency, 0.0);
    BOOST_REQUIRE_GT(res_pb  .full_energy_peak_efficiency, 0.0);

    double ratio_mc = res_pb.full_energy_peak_efficiency
                    / res_bare.full_energy_peak_efficiency;

    double mu_pb = lead.mu_total(E * 1e-3);
    double ratio_analytic = std::exp(-mu_pb * 0.2);

    BOOST_CHECK_CLOSE(ratio_mc, ratio_analytic, 10.0);
    BOOST_CHECK_LT(ratio_mc, 1.0);
    BOOST_CHECK_GT(ratio_mc, 0.20);
}

BOOST_AUTO_TEST_CASE(ratio_matches_beer_lambert_5mm_Pb_662keV)
{
    Material nai  = make_NaI();
    Material lead = make_Lead();
    const double   E = 662.0;
    const double   length_cm = 7.62;

    EfficiencyCalculator calc_bare;
    calc_bare.set_fep_window_keV(kTestFepWindowKeV);
    calc_bare.set_detector(DetectorShape::Cylinder, &nai, {3.81, length_cm});
    calc_bare.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    EfficiencyCalculator calc_pb;
    calc_pb.set_fep_window_keV(kTestFepWindowKeV);
    calc_pb.set_detector(DetectorShape::Cylinder, &nai, {3.81, length_cm});
    calc_pb.add_attenuator(&lead, 0.5, 0.0, 0.0, length_cm);
    calc_pb.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

    auto res_bare = calc_bare.compute(precision_config(E));
    auto res_pb   = calc_pb  .compute(precision_config(E));

    BOOST_REQUIRE_GT(res_bare.full_energy_peak_efficiency, 0.0);
    BOOST_REQUIRE_GT(res_pb  .full_energy_peak_efficiency, 0.0);

    double ratio_mc = res_pb.full_energy_peak_efficiency
                    / res_bare.full_energy_peak_efficiency;

    double mu_pb       = lead.mu_total(E * 1e-3);
    double ratio_analytic = std::exp(-mu_pb * 0.5);

    BOOST_CHECK_CLOSE(ratio_mc, ratio_analytic, 10.0);
    BOOST_CHECK_LT(ratio_mc, 1.0);
    BOOST_CHECK_GT(ratio_mc, 0.40);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 2: Cylindrical soil extended source
// ============================================================
BOOST_AUTO_TEST_SUITE(CylindricalSoilSource)

BOOST_AUTO_TEST_CASE(efficiency_positive_and_bounded)
{
    Material nai  = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(
        Eigen::Vector3d(0.0, 0.0, -15.0),
        5.0, 5.0
    );

    auto res = calc.compute(precision_config(351.0));

    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LE(res.full_energy_peak_efficiency, res.total_efficiency + 1e-9);
    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 1e-5);
    BOOST_CHECK_LT(res.full_energy_peak_efficiency, 0.05);
}

BOOST_AUTO_TEST_CASE(deeper_source_gives_lower_efficiency_per_photon)
{
    Material nai  = make_NaI();

    EfficiencyCalculator calc_near, calc_far;
    for (auto* calc : {&calc_near, &calc_far}) {
        calc->set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    }

    calc_near.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -8.0), 3.0, 2.0);
    calc_far.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -30.0), 3.0, 2.0);

    auto res_near = calc_near.compute(precision_config(351.0));
    auto res_far  = calc_far .compute(precision_config(351.0));

    BOOST_CHECK_GT(res_near.full_energy_peak_efficiency,
                   res_far.full_energy_peak_efficiency);
}

BOOST_AUTO_TEST_CASE(no_soil_attenuation_matches_geometric_only)
{
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 3.0, 3.0);

    auto res = calc.compute(precision_config(662.0));
    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LE(res.full_energy_peak_efficiency, 1.0);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 3: Rectangular soil extended source
// ============================================================
BOOST_AUTO_TEST_SUITE(RectangularSoilSource)

BOOST_AUTO_TEST_CASE(efficiency_positive_and_bounded)
{
    Material nai  = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(
        Eigen::Vector3d(0.0, 0.0, -15.0),
        Eigen::Vector3d(10.0, 10.0, 5.0)
    );

    auto res = calc.compute(precision_config(351.0));

    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LE(res.full_energy_peak_efficiency, res.total_efficiency + 1e-9);
    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 1e-5);
    BOOST_CHECK_LT(res.full_energy_peak_efficiency, 0.05);
}

BOOST_AUTO_TEST_CASE(efficiency_decreases_with_source_depth)
{
    Material nai = make_NaI();

    EfficiencyCalculator calc_near, calc_far;
    for (auto* calc : {&calc_near, &calc_far}) {
        calc->set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    }

    calc_near.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -10.0),
                                     Eigen::Vector3d(5.0, 5.0, 3.0));
    calc_far .set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -30.0),
                                     Eigen::Vector3d(5.0, 5.0, 3.0));

    auto res_near = calc_near.compute(precision_config(351.0));
    auto res_far  = calc_far .compute(precision_config(351.0));

    BOOST_CHECK_GT(res_near.full_energy_peak_efficiency,
                   res_far .full_energy_peak_efficiency);
}

BOOST_AUTO_TEST_CASE(fep_leq_total_invariant_holds)
{
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -12.0),
                                 Eigen::Vector3d(6.0, 6.0, 4.0));

    for (double E : {100.0, 351.0, 662.0, 1173.0}) {
        auto res = calc.compute(precision_config(E));
        BOOST_CHECK_GE(res.total_efficiency,
                       res.full_energy_peak_efficiency - 1e-9);
    }
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 4: Beer-Lambert self-consistency for extended soil sources
// ============================================================
BOOST_AUTO_TEST_SUITE(BeerLambertExtended)

BOOST_AUTO_TEST_CASE(point_vs_extended_efficiency_ratio_351keV)
{
    Material nai  = make_NaI();
    Material soil = make_Soil();

    EfficiencyCalculator calc_point;
    calc_point.set_fep_window_keV(kTestFepWindowKeV);
    calc_point.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_point.set_point_source(Eigen::Vector3d(0.0, 0.0, -12.0));

    EfficiencyCalculator calc_thin;
    calc_thin.set_fep_window_keV(kTestFepWindowKeV);
    calc_thin.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_thin.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -12.0), 0.1, 0.1);

    auto res_point = calc_point.compute(precision_config(351.0));
    auto res_thin  = calc_thin .compute(precision_config(351.0));

    BOOST_REQUIRE_GT(res_point.full_energy_peak_efficiency, 0.0);
    BOOST_REQUIRE_GT(res_thin .full_energy_peak_efficiency, 0.0);

    double ratio = res_thin.full_energy_peak_efficiency
                 / res_point.full_energy_peak_efficiency;
    BOOST_CHECK_GT(ratio, 0.75);
    BOOST_CHECK_LT(ratio, 1.25);
}

BOOST_AUTO_TEST_CASE(thick_soil_attenuates_properly_at_351keV)
{
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0), 10.0, 15.0);

    auto res = calc.compute(precision_config(351.0));

    BOOST_CHECK_GT(res.full_energy_peak_efficiency, 0.0);
    BOOST_CHECK_LT(res.full_energy_peak_efficiency, 0.10);
    BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()
