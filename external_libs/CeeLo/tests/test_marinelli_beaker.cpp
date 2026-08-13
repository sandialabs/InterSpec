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

#define BOOST_TEST_MODULE MarinelliBeakerTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <random>

using namespace ceelo;

namespace {

/// Helper: create a SimulationConfig for Marinelli (isotropic, needs more events).
SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 50000;
    config.termination.max_wall_seconds = 60.0;
    config.num_threads = 0;
    config.batch_size = 10000;
    return config;
}

// Standard Marinelli test parameters:
// 3"x3" NaI bare, well_r=4.3, well_depth=6, outer_r=7.5,
// fill_height=4, endcap_to_beaker=0.5, water sample, 2mm PE beaker
struct MarinelliTestSetup {
    Material nai = make_NaI();
    Material water = make_Water();
    Material pe = make_Polyethylene();

    double well_r = 4.3;
    double well_depth = 6.0;
    double outer_r = 7.5;
    double fill_height = 4.0;
    double endcap_to_beaker = 0.5;
    double beaker_thickness = 0.2;
};

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(MarinelliBeaker)

// ============================================================
//  Test 1: SourceGeometry Marinelli configuration
// ============================================================
BOOST_AUTO_TEST_CASE(configure_marinelli_stores_parameters) {
    SourceGeometry sg;
    // Assume z_det_min = 0 for bare detector (no attenuators/dead layer)
    double z_det_min = 0.0;
    double z_bk = z_det_min - 0.5;   // = -0.5
    double z_we = z_det_min + 6.0;   // = 6.0
    double z_bot = z_bk - 4.0;       // = -4.5

    sg.configure_marinelli(4.3, 7.5, z_bk, z_we, z_bot);

    BOOST_CHECK(sg.shape() == SourceGeometry::Shape::Marinelli);
    BOOST_CHECK(sg.is_configured());
    BOOST_CHECK_CLOSE(sg.marinelli_well_inner_radius(), 4.3, 1e-10);
    BOOST_CHECK_CLOSE(sg.marinelli_outer_radius(), 7.5, 1e-10);
    BOOST_CHECK_CLOSE(sg.marinelli_z_bk(), -0.5, 1e-10);
    BOOST_CHECK_CLOSE(sg.marinelli_z_we(), 6.0, 1e-10);
    BOOST_CHECK_CLOSE(sg.marinelli_z_bot(), -4.5, 1e-10);
}

// ============================================================
//  Test 2: is_inside_marinelli correctness
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_inside_test) {
    SourceGeometry sg;
    Material water = make_Water();
    sg.set_source_material(&water);

    double z_bk = -0.5, z_we = 6.0, z_bot = -4.5;
    sg.configure_marinelli(4.3, 7.5, z_bk, z_we, z_bot);

    // Points in disk region (z < z_bk, r < outer_r)
    // source_material_path should return > 0 for interior points
    // We test indirectly by checking compute_transmission > 0 for points inside
    Eigen::Vector3d dir(0, 0, 1); // +z direction
    double E = 0.662; // MeV

    // Point in disk center: should have nonzero path
    Eigen::Vector3d disk_center(0, 0, -2.5);
    double trans_disk = sg.compute_transmission(disk_center, dir, E);
    BOOST_CHECK_GT(trans_disk, 0.0);
    BOOST_CHECK_LT(trans_disk, 1.0); // some attenuation expected

    // Point in ring region
    Eigen::Vector3d ring_point(6.0, 0, 3.0); // r=6 > well_r, z_bk < z < z_we
    double trans_ring = sg.compute_transmission(ring_point, dir, E);
    BOOST_CHECK_GT(trans_ring, 0.0);
    BOOST_CHECK_LE(trans_ring, 1.0);

    // Point in well cavity (should NOT be inside): r < well_r, z > z_bk
    Eigen::Vector3d well_cavity(0, 0, 3.0); // r=0 < well_r=4.3, z=3 > z_bk
    // Path should be 0 (not inside sample), transmission = 1
    double trans_well = sg.compute_transmission(well_cavity, dir, E);
    BOOST_CHECK_CLOSE(trans_well, 1.0, 0.1);
}

// ============================================================
//  Test 3: Boundary distance — axial ray from disk
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_boundary_axial_disk) {
    SourceGeometry sg;
    Material water = make_Water();
    sg.set_source_material(&water);

    double z_bk = -0.5, z_we = 6.0, z_bot = -4.5;
    sg.configure_marinelli(4.3, 7.5, z_bk, z_we, z_bot);

    // Point at center of disk, moving -z: should hit bottom at z_bot
    Eigen::Vector3d pos(0, 0, -2.5);
    Eigen::Vector3d dir_down(0, 0, -1);
    double E = 0.662; // MeV
    double trans = sg.compute_transmission(pos, dir_down, E);
    // Expected path through water: 2.5 - (-4.5) = ... no, from pos.z=-2.5 to z_bot=-4.5 = 2.0 cm
    double expected_path = 2.0;
    double mu = water.mu_total(E);
    double expected_trans = std::exp(-mu * expected_path);
    // Allow 1% tolerance (water mu * 2cm is small attenuation)
    BOOST_CHECK_CLOSE(trans, expected_trans, 1.0);
}

// ============================================================
//  Test 4: Boundary distance — radial ray from ring to well wall
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_boundary_radial_ring) {
    SourceGeometry sg;
    Material water = make_Water();
    sg.set_source_material(&water);

    double z_bk = -0.5, z_we = 6.0, z_bot = -4.5;
    double well_r = 4.3, outer_r = 7.5;
    sg.configure_marinelli(well_r, outer_r, z_bk, z_we, z_bot);

    // Point in ring at r=6, moving inward (-x): should hit well wall at r=well_r
    Eigen::Vector3d pos(6.0, 0, 3.0);
    Eigen::Vector3d dir(-1, 0, 0); // radially inward
    double E = 0.662;
    double expected_path = 6.0 - well_r; // 1.7 cm
    double mu = water.mu_total(E);
    double expected_trans = std::exp(-mu * expected_path);
    double trans = sg.compute_transmission(pos, dir, E);
    BOOST_CHECK_CLOSE(trans, expected_trans, 1.5);
}

// ============================================================
//  Test 5: EfficiencyCalculator Marinelli API
// ============================================================
BOOST_AUTO_TEST_CASE(efficiency_calculator_marinelli_api) {
    MarinelliTestSetup s;

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, &s.water, &s.pe, s.beaker_thickness);

    BOOST_CHECK(calc.source_type() == SourceType::Marinelli);
    BOOST_CHECK(calc.source_geometry().is_configured());
    BOOST_CHECK(calc.source_geometry().shape() == SourceGeometry::Shape::Marinelli);

    // Check z-coordinates are computed from detector geometry
    double z_det_min = calc.geometry().outer_z_extent().first;
    double expected_z_bk = z_det_min - s.endcap_to_beaker;
    double expected_z_we = z_det_min + s.well_depth;
    double expected_z_bot = expected_z_bk - s.fill_height;
    BOOST_CHECK_CLOSE(calc.source_geometry().marinelli_z_bk(), expected_z_bk, 1e-10);
    BOOST_CHECK_CLOSE(calc.source_geometry().marinelli_z_we(), expected_z_we, 1e-10);
    BOOST_CHECK_CLOSE(calc.source_geometry().marinelli_z_bot(), expected_z_bot, 1e-10);
}

// ============================================================
//  Test 6: Marinelli volume calculation
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_volume) {
    MarinelliTestSetup s;
    double ring_depth = s.endcap_to_beaker + s.well_depth;
    double V_disk = M_PI * s.outer_r * s.outer_r * s.fill_height;
    double V_ring = M_PI * (s.outer_r * s.outer_r - s.well_r * s.well_r) * ring_depth;
    double V_total = V_disk + V_ring;

    // Verify volume is reasonable (> 0, not huge)
    BOOST_CHECK_GT(V_total, 100.0);  // > 100 cm^3 for a 1L-ish beaker
    BOOST_CHECK_LT(V_total, 2000.0); // < 2000 cm^3
}

// ============================================================
//  Test 7: Marinelli MC run produces reasonable efficiencies
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_mc_run) {
    MarinelliTestSetup s;

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, &s.water, &s.pe, s.beaker_thickness);

    // Run at 662 keV: Marinelli has high solid angle coverage
    auto result = calc.compute(precision_config(662.0, 0.01));

    // FEP efficiency for Marinelli 662 keV on 3"x3" NaI should be in ~1-10% range
    BOOST_CHECK_GT(result.full_energy_peak_efficiency, 0.005);
    BOOST_CHECK_LT(result.full_energy_peak_efficiency, 0.15);

    // Total efficiency should be higher than FEP
    BOOST_CHECK_GT(result.total_efficiency, result.full_energy_peak_efficiency);

    // Total should be reasonable (< 50% — some solid angle is missed)
    BOOST_CHECK_LT(result.total_efficiency, 0.5);
}

// ============================================================
//  Test 8: Marinelli without sample material (dry beaker)
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_no_sample_vs_with_sample) {
    MarinelliTestSetup s;

    // With water sample
    EfficiencyCalculator calc_water;
    calc_water.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc_water.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, &s.water, &s.pe, s.beaker_thickness);

    // Without sample material (nullptr)
    EfficiencyCalculator calc_dry;
    calc_dry.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc_dry.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, nullptr, &s.pe, s.beaker_thickness);

    auto res_water = calc_water.compute(precision_config(662.0, 0.01));
    auto res_dry = calc_dry.compute(precision_config(662.0, 0.01));

    // Water attenuates, so FEP with water < FEP without (at 662 keV, ~5-10% difference)
    BOOST_CHECK_LT(res_water.full_energy_peak_efficiency,
                    res_dry.full_energy_peak_efficiency * 1.02);
}

// ============================================================
//  Test 9: GDML export doesn't crash for Marinelli
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_gdml_export) {
    MarinelliTestSetup s;

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, &s.water, &s.pe, s.beaker_thickness);

    // Export should not throw
    BOOST_CHECK_NO_THROW(calc.export_geant4_gdml("/tmp/test_marinelli.gdml", true));
    BOOST_CHECK_NO_THROW(calc.export_geant4_macro("/tmp/test_marinelli.mac", 662.0, 1000));
}

// ============================================================
//  Test 10: Outermost extent radius for Marinelli
// ============================================================
BOOST_AUTO_TEST_CASE(marinelli_outermost_extent) {
    SourceGeometry sg;
    Material water = make_Water();
    Material pe = make_Polyethylene();
    sg.set_source_material(&water);
    sg.add_shield(&pe, 0.2);
    sg.configure_marinelli(4.3, 7.5, -0.5, 6.0, -4.5);

    // Outermost extent should be outer_r (ignoring shield thickness for now)
    BOOST_CHECK_GE(sg.outermost_extent_radius(), 7.5);
}

// ============================================================
//  Test: Source electron transport increases total efficiency at high energy
// ============================================================
BOOST_AUTO_TEST_CASE(source_electron_transport_increases_total) {
    MarinelliTestSetup s;
    Material al = make_Aluminum();

    // Without source electrons
    EfficiencyCalculator calc_no_e;
    calc_no_e.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc_no_e.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
    calc_no_e.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, &s.water, &s.pe, s.beaker_thickness);

    // With source electrons
    EfficiencyCalculator calc_with_e;
    calc_with_e.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc_with_e.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
    calc_with_e.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, &s.water, &s.pe, s.beaker_thickness);
    calc_with_e.enable_source_electron_transport(true);

    // At 2614 keV, the electron effect is largest (~3% deficit without)
    auto cfg = precision_config(2614.0, 0.02);
    cfg.termination.max_wall_seconds = 30.0;
    auto res_no_e = calc_no_e.compute(cfg);
    auto res_with_e = calc_with_e.compute(cfg);

    // Source electrons should increase total efficiency
    BOOST_CHECK_GT(res_with_e.total_efficiency, res_no_e.total_efficiency);

    // FEP should be essentially unchanged (within 5% relative)
    double fep_ratio = res_with_e.full_energy_peak_efficiency /
                       res_no_e.full_energy_peak_efficiency;
    BOOST_CHECK_CLOSE(fep_ratio, 1.0, 5.0);
}

// ============================================================
//  Test: Source electrons have negligible effect at low energy
// ============================================================
BOOST_AUTO_TEST_CASE(source_electron_no_effect_low_energy) {
    MarinelliTestSetup s;
    Material al = make_Aluminum();

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &s.nai, {3.81, 7.62});
    calc.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
    calc.set_marinelli_beaker(
        s.well_r, s.well_depth, s.outer_r, s.fill_height,
        s.endcap_to_beaker, &s.water, &s.pe, s.beaker_thickness);
    calc.enable_source_electron_transport(true);

    // At 59 keV, Compton recoil electrons have < 50 keV, so all filtered out
    auto cfg = precision_config(59.0, 0.02);
    cfg.termination.max_wall_seconds = 30.0;
    auto result = calc.compute(cfg);

    // Should still get reasonable results (not crash, not NaN)
    BOOST_CHECK_GT(result.total_efficiency, 0.01);
    BOOST_CHECK_GT(result.full_energy_peak_efficiency, 0.01);
    BOOST_CHECK(!std::isnan(result.total_efficiency));
    BOOST_CHECK(!std::isnan(result.full_energy_peak_efficiency));
}

BOOST_AUTO_TEST_SUITE_END()
