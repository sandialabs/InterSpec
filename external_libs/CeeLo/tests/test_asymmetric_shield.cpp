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

#define BOOST_TEST_MODULE AsymmetricShieldTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <iomanip>
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
//  Per-dimension source shield thicknesses
// ============================================================

BOOST_AUTO_TEST_SUITE(AsymmetricShield)

BOOST_AUTO_TEST_CASE(analytic_cylinder_asymmetric) {
    // Cylindrical shield with different radial and end-cap thicknesses:
    // transmission from the center along the axis sees only t_end, and
    // radially only t_radial.
    Material pb = make_Lead();

    SourceGeometry sg;
    sg.configure_cylindrical(Eigen::Vector3d(0, 0, 0), 2.0, 3.0,
                             Eigen::Matrix3d::Identity());
    sg.add_shield(&pb, 0.3, 0.1);  // radial 3 mm, end caps 1 mm

    double energy_MeV = 0.662;
    double mu = pb.mu_total(energy_MeV);

    double trans_axial = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), energy_MeV);
    double trans_radial = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), energy_MeV);

    BOOST_CHECK_CLOSE(trans_axial, std::exp(-mu * 0.1), 1e-6);
    BOOST_CHECK_CLOSE(trans_radial, std::exp(-mu * 0.3), 1e-6);
}

BOOST_AUTO_TEST_CASE(analytic_cylinder_zero_end_cap) {
    // Side-wall-only shield: no attenuation along the axis.
    Material pb = make_Lead();

    SourceGeometry sg;
    sg.configure_cylindrical(Eigen::Vector3d(0, 0, 0), 2.0, 3.0,
                             Eigen::Matrix3d::Identity());
    sg.add_shield(&pb, 0.3, 0.0);

    double energy_MeV = 0.662;
    double mu = pb.mu_total(energy_MeV);

    double trans_axial = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, -1), energy_MeV);
    double trans_radial = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 0), energy_MeV);

    BOOST_CHECK_EQUAL(trans_axial, 1.0);
    BOOST_CHECK_CLOSE(trans_radial, std::exp(-mu * 0.3), 1e-6);
}

BOOST_AUTO_TEST_CASE(analytic_cylinder_zero_side_wall) {
    // End-caps-only shield: no attenuation radially.
    Material pb = make_Lead();

    SourceGeometry sg;
    sg.configure_cylindrical(Eigen::Vector3d(0, 0, 0), 2.0, 3.0,
                             Eigen::Matrix3d::Identity());
    sg.add_shield(&pb, 0.0, 0.2);

    double energy_MeV = 0.662;
    double mu = pb.mu_total(energy_MeV);

    double trans_axial = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), energy_MeV);
    double trans_radial = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), energy_MeV);

    BOOST_CHECK_CLOSE(trans_axial, std::exp(-mu * 0.2), 1e-6);
    BOOST_CHECK_EQUAL(trans_radial, 1.0);
}

BOOST_AUTO_TEST_CASE(fep_only_transmission_zero_dimensions) {
    // FEP-only stochastic-Rayleigh transmission: zero path through a
    // zero-thickness dimension gives weight exactly 1; a nonzero path
    // weights by exp(-mu_no_rayleigh * path) deterministically.
    Material pb = make_Lead();
    double energy_MeV = 0.662;
    MacroscopicXS xs = pb.macroscopic_xs(energy_MeV);
    double mu_no_rs = xs.mu_pe + xs.mu_cs + xs.mu_pp;

    std::mt19937_64 rng(12345);

    SourceGeometry sg_side;
    sg_side.configure_cylindrical(Eigen::Vector3d(0, 0, 0), 2.0, 3.0,
                                  Eigen::Matrix3d::Identity());
    sg_side.add_shield(&pb, 0.3, 0.0);
    auto res_axial = sg_side.compute_transmission_fep_only(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), energy_MeV, rng);
    BOOST_CHECK_EQUAL(res_axial.weight, 1.0);

    SourceGeometry sg_ends;
    sg_ends.configure_cylindrical(Eigen::Vector3d(0, 0, 0), 2.0, 3.0,
                                  Eigen::Matrix3d::Identity());
    sg_ends.add_shield(&pb, 0.0, 0.2);
    auto res_radial = sg_ends.compute_transmission_fep_only(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), energy_MeV, rng);
    BOOST_CHECK_EQUAL(res_radial.weight, 1.0);

    auto res_axial2 = sg_ends.compute_transmission_fep_only(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), energy_MeV, rng);
    BOOST_CHECK_CLOSE(res_axial2.weight, std::exp(-mu_no_rs * 0.2), 1e-6);
}

BOOST_AUTO_TEST_CASE(analytic_rectangular_per_axis) {
    // Box shield with per-axis thicknesses, including a zero dimension.
    Material pb = make_Lead();

    SourceGeometry sg;
    sg.configure_rectangular(Eigen::Vector3d(0, 0, 0),
                             Eigen::Vector3d(1.0, 2.0, 3.0),
                             Eigen::Matrix3d::Identity());
    sg.add_shield(&pb, 0.1, 0.0, 0.3);

    double energy_MeV = 0.662;
    double mu = pb.mu_total(energy_MeV);

    double trans_x = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), energy_MeV);
    double trans_y = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, -1, 0), energy_MeV);
    double trans_z = sg.compute_transmission(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), energy_MeV);

    BOOST_CHECK_CLOSE(trans_x, std::exp(-mu * 0.1), 1e-6);
    BOOST_CHECK_EQUAL(trans_y, 1.0);
    BOOST_CHECK_CLOSE(trans_z, std::exp(-mu * 0.3), 1e-6);
}

BOOST_AUTO_TEST_CASE(uniform_overload_equivalence) {
    // The per-dimension overloads with equal components must reproduce the
    // uniform overload exactly, for oblique directions too.
    Material pb = make_Lead();
    double t = 0.15;
    double energy_MeV = 0.2;

    SourceGeometry cyl_uniform, cyl_multi;
    for (SourceGeometry* sg : {&cyl_uniform, &cyl_multi}) {
        sg->configure_cylindrical(Eigen::Vector3d(0, 0, 0), 2.0, 3.0,
                                  Eigen::Matrix3d::Identity());
    }
    cyl_uniform.add_shield(&pb, t);
    cyl_multi.add_shield(&pb, t, t);

    SourceGeometry box_uniform, box_multi;
    for (SourceGeometry* sg : {&box_uniform, &box_multi}) {
        sg->configure_rectangular(Eigen::Vector3d(0, 0, 0),
                                  Eigen::Vector3d(1.0, 2.0, 3.0),
                                  Eigen::Matrix3d::Identity());
    }
    box_uniform.add_shield(&pb, t);
    box_multi.add_shield(&pb, t, t, t);

    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (int i = 0; i < 10; ++i) {
        Eigen::Vector3d dir(u(rng), u(rng), u(rng));
        if (dir.norm() < 1e-3) continue;
        dir.normalize();

        double tc_u = cyl_uniform.compute_transmission(
            Eigen::Vector3d(0.5, -0.3, 1.0), dir, energy_MeV);
        double tc_m = cyl_multi.compute_transmission(
            Eigen::Vector3d(0.5, -0.3, 1.0), dir, energy_MeV);
        BOOST_CHECK_EQUAL(tc_u, tc_m);

        double tb_u = box_uniform.compute_transmission(
            Eigen::Vector3d(0.2, 0.8, -1.5), dir, energy_MeV);
        double tb_m = box_multi.compute_transmission(
            Eigen::Vector3d(0.2, 0.8, -1.5), dir, energy_MeV);
        BOOST_CHECK_EQUAL(tb_u, tb_m);
    }
}

BOOST_AUTO_TEST_CASE(outermost_extent_radius_asymmetric) {
    Material pb = make_Lead();

    SourceGeometry cyl;
    cyl.configure_cylindrical(Eigen::Vector3d(0, 0, -15), 3.0, 4.0,
                              Eigen::Matrix3d::Identity());
    cyl.add_shield(&pb, 1.0, 2.0);
    BOOST_CHECK_CLOSE(cyl.outermost_extent_radius(),
                      std::sqrt(4.0 * 4.0 + 6.0 * 6.0), 1e-9);

    SourceGeometry box;
    box.configure_rectangular(Eigen::Vector3d(0, 0, -15),
                              Eigen::Vector3d(1.0, 2.0, 3.0),
                              Eigen::Matrix3d::Identity());
    box.add_shield(&pb, 0.5, 0.0, 1.0);
    BOOST_CHECK_CLOSE(box.outermost_extent_radius(),
                      Eigen::Vector3d(1.5, 2.0, 4.0).norm(), 1e-9);

    // Point-source extent semantics unchanged.
    SourceGeometry pt;
    pt.configure_point(Eigen::Vector3d(0, 0, -10));
    pt.add_shield(&pb, 0.1);
    pt.add_shield(&pb, 0.2);
    BOOST_CHECK_CLOSE(pt.outermost_extent_radius(), 0.3, 1e-9);
}

BOOST_AUTO_TEST_CASE(mc_end_cap_only_attenuation) {
    // MC integration: small on-axis cylindrical source facing a 3"x3" NaI.
    // An end-caps-only Pb shield should attenuate FEP by ~exp(-mu*t)
    // (slightly more due to oblique paths); a side-wall-only shield should
    // leave FEP essentially unchanged.
    Material pb = make_Lead();
    Material nai = make_NaI();

    double t = 0.1;  // 1 mm Pb
    double energy_keV = 662.0;
    double mu = pb.mu_total(energy_keV * 1e-3);
    double expected_trans = std::exp(-mu * t);

    auto make_calc = [&](EfficiencyCalculator& calc) {
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0), 0.5, 0.5);
    };

    EfficiencyCalculator calc_bare;
    calc_bare.set_fep_window_keV(kTestFepWindowKeV);
    make_calc(calc_bare);
    auto res_bare = calc_bare.compute(precision_config(energy_keV));

    EfficiencyCalculator calc_ends;
    calc_ends.set_fep_window_keV(kTestFepWindowKeV);
    make_calc(calc_ends);
    calc_ends.add_source_shield(&pb, 0.0, t);
    auto res_ends = calc_ends.compute(precision_config(energy_keV));

    EfficiencyCalculator calc_side;
    calc_side.set_fep_window_keV(kTestFepWindowKeV);
    make_calc(calc_side);
    calc_side.add_source_shield(&pb, t, 0.0);
    auto res_side = calc_side.compute(precision_config(energy_keV));

    double ratio_ends = res_ends.full_energy_peak_efficiency /
                        res_bare.full_energy_peak_efficiency;
    double ratio_side = res_side.full_energy_peak_efficiency /
                        res_bare.full_energy_peak_efficiency;

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "662 keV: bare=" << res_bare.full_energy_peak_efficiency
        << " +/- " << res_bare.fep_uncertainty
        << "  end-caps=" << res_ends.full_energy_peak_efficiency
        << " +/- " << res_ends.fep_uncertainty
        << "  side-only=" << res_side.full_energy_peak_efficiency
        << " +/- " << res_side.fep_uncertainty
        << std::setprecision(4)
        << "  ratio_ends=" << ratio_ends << " (expect ~" << expected_trans << ")"
        << "  ratio_side=" << ratio_side << " (expect ~1)";
    BOOST_TEST_MESSAGE(msg.str());

    BOOST_CHECK_CLOSE(ratio_ends, expected_trans, 3.0);
    BOOST_CHECK_CLOSE(ratio_side, 1.0, 3.0);
}

BOOST_AUTO_TEST_CASE(fep_only_matches_full_mode_asymmetric_shield) {
    // FEP-only vs full mode with an asymmetric cylindrical source shield.
    Material nai = make_NaI();
    Material pb = make_Lead();

    auto make_calc = [&](EfficiencyCalculator& calc) {
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 0.5, 0.5);
        calc.add_source_shield(&pb, 0.1, 0.025);  // 1 mm side, 0.25 mm ends
    };

    EfficiencyCalculator calc_full;
    calc_full.set_fep_window_keV(kTestFepWindowKeV);
    make_calc(calc_full);
    calc_full.enable_fep_only_mode(false);
    auto res_full = calc_full.compute(precision_config(200.0));

    EfficiencyCalculator calc_fep;
    calc_fep.set_fep_window_keV(kTestFepWindowKeV);
    make_calc(calc_fep);
    calc_fep.enable_fep_only_mode(true);
    auto res_fep = calc_fep.compute(precision_config(200.0));

    double sigma = std::sqrt(
        res_full.fep_uncertainty * res_full.fep_uncertainty +
        res_fep.fep_uncertainty * res_fep.fep_uncertainty);
    double diff = std::abs(res_full.full_energy_peak_efficiency -
                           res_fep.full_energy_peak_efficiency);

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "asymmetric Pb shield 200 keV: full=" << res_full.full_energy_peak_efficiency
        << " +/- " << res_full.fep_uncertainty
        << "  fep=" << res_fep.full_energy_peak_efficiency
        << " +/- " << res_fep.fep_uncertainty
        << "  (" << ((sigma > 0.0) ? diff / sigma : 0.0) << " sigma)";
    BOOST_TEST_MESSAGE(msg.str());

    BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-6);
}

BOOST_AUTO_TEST_CASE(gdml_export_asymmetric_cylinder) {
    Material pb = make_Lead();
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0), 3.0, 3.0);
    calc.add_source_shield(&pb, 0.3, 0.1);

    std::string fname = "test_gdml_asym_cyl.gdml";
    calc.export_geant4_gdml(fname);

    std::ifstream f(fname);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    // Outer solid grows by t_radial in r and t_end in half-z.
    BOOST_CHECK(content.find("SrcShieldOuterSolid0\" rmin=\"0\" rmax=\"3.300000\" z=\"6.200000\"")
                != std::string::npos);
    // Inner (subtracted) solid is inflated by kZeroDimEps = 1e-4 cm on BOTH
    // axes so its faces never coincide with the source-material surface it
    // wraps (the validated coincident-surface fix; thins the wall ~1 micron).
    BOOST_CHECK(content.find("SrcShieldInnerSolid0\" rmin=\"0\" rmax=\"3.000100\" z=\"6.000200\"")
                != std::string::npos);

    std::remove(fname.c_str());
}

BOOST_AUTO_TEST_CASE(gdml_export_zero_end_cap_epsilon) {
    Material pb = make_Lead();
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0), 3.0, 3.0);
    calc.add_source_shield(&pb, 0.3, 0.0);

    std::string fname = "test_gdml_zero_end_cyl.gdml";
    calc.export_geant4_gdml(fname);

    std::ifstream f(fname);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    // Outer half-z unchanged (t_end = 0); subtracted solid inflated by the
    // epsilon on both axes (radial and z) to avoid coincident boolean surfaces.
    BOOST_CHECK(content.find("SrcShieldOuterSolid0\" rmin=\"0\" rmax=\"3.300000\" z=\"6.000000\"")
                != std::string::npos);
    BOOST_CHECK(content.find("SrcShieldInnerSolid0\" rmin=\"0\" rmax=\"3.000100\" z=\"6.000200\"")
                != std::string::npos);

    std::remove(fname.c_str());
}

BOOST_AUTO_TEST_CASE(gdml_export_uniform_epsilon) {
    // Uniform shields are ALSO inflated by the epsilon: shield 0's inner face
    // coincides with the source-material surface regardless of symmetry, so the
    // subtracted solid grows by kZeroDimEps = 1e-4 cm on every axis.
    Material pb = make_Lead();
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -15.0), 3.0, 3.0);
    calc.add_source_shield(&pb, 0.2);

    std::string fname = "test_gdml_uniform_cyl.gdml";
    calc.export_geant4_gdml(fname);

    std::ifstream f(fname);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    BOOST_CHECK(content.find("SrcShieldOuterSolid0\" rmin=\"0\" rmax=\"3.200000\" z=\"6.400000\"")
                != std::string::npos);
    BOOST_CHECK(content.find("SrcShieldInnerSolid0\" rmin=\"0\" rmax=\"3.000100\" z=\"6.000200\"")
                != std::string::npos);
    // The un-inflated source dims must NOT appear as the subtracted solid.
    BOOST_CHECK(content.find("SrcShieldInnerSolid0\" rmin=\"0\" rmax=\"3.000000\" z=\"6.000000\"")
                == std::string::npos);

    std::remove(fname.c_str());
}

BOOST_AUTO_TEST_CASE(gdml_export_asymmetric_box) {
    Material pb = make_Lead();
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -15.0),
                                Eigen::Vector3d(1.0, 2.0, 3.0));
    calc.add_source_shield(&pb, 0.1, 0.0, 0.3);

    std::string fname = "test_gdml_asym_box.gdml";
    calc.export_geant4_gdml(fname);

    std::ifstream f(fname);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    BOOST_CHECK(content.find("SrcShieldOuterSolid0\" x=\"2.200000\" y=\"4.000000\" z=\"6.600000\"")
                != std::string::npos);
    // Subtracted solid inflated by epsilon on ALL axes (x, y, z), not just the
    // zero-thickness y face.
    BOOST_CHECK(content.find("SrcShieldInnerSolid0\" x=\"2.000200\" y=\"4.000200\" z=\"6.000200\"")
                != std::string::npos);

    std::remove(fname.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
