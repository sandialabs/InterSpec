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

#define BOOST_TEST_MODULE SourceShieldingTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <cmath>
#include <fstream>
#include <string>

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
//  Feature 3: Source Self-Attenuation & Shielding
// ============================================================

BOOST_AUTO_TEST_SUITE(SourceShielding)

BOOST_AUTO_TEST_CASE(point_source_pb_spherical_shell) {
    // Point source with Pb shielding: MC transmission should match exp(-mu*t).
    Material pb = make_Lead();
    Material nai = make_NaI();

    double thickness = 0.1; // 1 mm Pb
    double energy_keV = 662.0;
    double energy_MeV = energy_keV * 1e-3;
    double mu = pb.mu_total(energy_MeV);
    double expected_trans = std::exp(-mu * thickness);

    // Without shielding
    EfficiencyCalculator calc_bare;
    calc_bare.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_bare.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    auto res_bare = calc_bare.compute(precision_config(energy_keV));

    // With Pb shielding
    EfficiencyCalculator calc_pb;
    calc_pb.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_pb.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_pb.add_source_shield(&pb, thickness);
    auto res_pb = calc_pb.compute(precision_config(energy_keV));

    // The ratio should match the expected transmission factor
    double ratio = res_pb.full_energy_peak_efficiency /
                   res_bare.full_energy_peak_efficiency;

    // With 0.5% precision on each, allow ~1.5% tolerance on the ratio
    BOOST_CHECK_CLOSE(ratio, expected_trans, 3.0);
}

BOOST_AUTO_TEST_CASE(source_geometry_transmission_point) {
    Material pb = make_Lead();

    SourceGeometry sg;
    sg.configure_point(Eigen::Vector3d(0, 0, -10));
    sg.add_shield(&pb, 0.1); // 1 mm

    double energy_MeV = 0.662;
    double mu = pb.mu_total(energy_MeV);
    double expected = std::exp(-mu * 0.1);

    double trans = sg.point_source_transmission(energy_MeV);
    BOOST_CHECK_CLOSE(trans, expected, 0.01);
}

BOOST_AUTO_TEST_CASE(source_geometry_no_effects_by_default) {
    SourceGeometry sg;
    BOOST_CHECK(!sg.has_source_effects());
}

BOOST_AUTO_TEST_CASE(source_material_reduces_efficiency) {
    // Source material (water) should reduce efficiency due to self-attenuation.
    Material nai = make_NaI();
    Material water = make_Water();

    EfficiencyCalculator calc_bare;
    calc_bare.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_bare.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0),
                                     3.0, 3.0);
    auto res_bare = calc_bare.compute(precision_config(200.0));

    EfficiencyCalculator calc_water;
    calc_water.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_water.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0),
                                      3.0, 3.0);
    calc_water.set_source_material(&water);
    auto res_water = calc_water.compute(precision_config(200.0));

    // Water-filled source should have lower FEP
    BOOST_CHECK_LT(res_water.full_energy_peak_efficiency,
                   res_bare.full_energy_peak_efficiency);
}

BOOST_AUTO_TEST_CASE(multi_shield_layers) {
    Material pb = make_Lead();
    Material nai = make_NaI();
    double energy_keV = 200.0;

    EfficiencyCalculator calc_1;
    calc_1.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_1.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_1.add_source_shield(&pb, 0.05); // 0.5 mm
    auto res_1 = calc_1.compute(precision_config(energy_keV));

    EfficiencyCalculator calc_2;
    calc_2.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_2.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc_2.add_source_shield(&pb, 0.05);
    calc_2.add_source_shield(&pb, 0.05); // total 1 mm
    auto res_2 = calc_2.compute(precision_config(energy_keV));

    BOOST_CHECK_LT(res_2.full_energy_peak_efficiency,
                   res_1.full_energy_peak_efficiency);
}

BOOST_AUTO_TEST_CASE(gdml_export_point_source_shielding) {
    Material pb = make_Lead();
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&pb, 0.1);

    std::string fname = "test_gdml_src_shield_point.gdml";
    calc.export_geant4_gdml(fname);

    std::ifstream f(fname);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    BOOST_CHECK(content.find("<sphere name=\"SrcShieldSolid0\"") != std::string::npos);
    BOOST_CHECK(content.find("SrcShieldLV0") != std::string::npos);
    BOOST_CHECK(content.find("SrcShieldPV0") != std::string::npos);
    BOOST_CHECK(content.find("\"Pb\"") != std::string::npos);
    BOOST_CHECK(content.find("SrcMaterialSolid") == std::string::npos);

    std::remove(fname.c_str());
}

BOOST_AUTO_TEST_CASE(gdml_export_cylindrical_source_with_material) {
    Material pb = make_Lead();
    Material nai = make_NaI();
    Material water = make_Water();

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0), 3.0, 3.0);
    calc.set_source_material(&water);
    calc.add_source_shield(&pb, 0.2);

    std::string fname = "test_gdml_src_shield_cyl.gdml";
    calc.export_geant4_gdml(fname);

    std::ifstream f(fname);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    BOOST_CHECK(content.find("SrcMaterialSolid") != std::string::npos);
    BOOST_CHECK(content.find("SrcMaterialLV") != std::string::npos);
    BOOST_CHECK(content.find("SrcMaterialPV") != std::string::npos);
    BOOST_CHECK(content.find("SrcShieldSolid0") != std::string::npos);
    BOOST_CHECK(content.find("SrcShieldOuterSolid0") != std::string::npos);
    BOOST_CHECK(content.find("SrcShieldInnerSolid0") != std::string::npos);
    BOOST_CHECK(content.find("Water") != std::string::npos);
    BOOST_CHECK(content.find("<sphere") == std::string::npos);

    std::remove(fname.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
