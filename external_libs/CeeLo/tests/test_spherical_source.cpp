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

#define BOOST_TEST_MODULE SphericalSourceTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "cascade/CascadeTypes.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>
#include <vector>

using namespace ceelo;

namespace {

/// SimulationConfig targeting `target` FEP relative precision (default 0.5%).
SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.target_total_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 5000;
    config.termination.max_wall_seconds = 30.0;
    config.num_threads = 0;  // auto-detect
    config.batch_size = 10000;
    return config;
}

/// z-score for the difference of two independent efficiencies.
double zdiff(double a, double sa, double b, double sb) {
    double s = std::sqrt(sa * sa + sb * sb);
    return (s > 0.0) ? std::abs(a - b) / s : 0.0;
}

/// SideOn cylinder rotation: cylinder local-z onto detector x (90° about y).
Eigen::Matrix3d side_on_rotation() {
    Eigen::Matrix3d R;
    R << 0, 0, 1,
         0, 1, 0,
        -1, 0, 0;
    return R;
}

} // anonymous namespace

// ============================================================
//  Spherical source — smoke + physics sanity
// ============================================================

BOOST_AUTO_TEST_SUITE(SphericalSource)

BOOST_AUTO_TEST_CASE(basic_invariants) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_spherical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 2.0);

    for (double E : {100.0, 662.0, 1332.0}) {
        auto res = calc.compute(E, 10000, 1);
        BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
        BOOST_CHECK_GE(res.full_energy_peak_efficiency, 0.0);
        BOOST_CHECK_LE(res.total_efficiency, 1.0);
        BOOST_CHECK_GT(res.total_efficiency, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(point_like_sphere_matches_point_source) {
    Material nai = make_NaI();
    Eigen::Vector3d center(0.0, 0.0, -10.0);

    EfficiencyCalculator calc_pt;
    calc_pt.set_fep_window_keV(kTestFepWindowKeV);
    calc_pt.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_pt.set_point_source(center);
    auto res_pt = calc_pt.compute(precision_config(662.0));

    EfficiencyCalculator calc_sph;
    calc_sph.set_fep_window_keV(kTestFepWindowKeV);
    calc_sph.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_sph.set_spherical_source(center, 0.001);
    auto res_sph = calc_sph.compute(precision_config(662.0));

    // Tiny sphere ≈ point at the same center (0.5% precision each → 5% is loose).
    BOOST_TEST_MESSAGE("point eff=" << res_pt.total_efficiency
        << " sphere eff=" << res_sph.total_efficiency);
    BOOST_CHECK_CLOSE(res_sph.total_efficiency, res_pt.total_efficiency, 5.0);
}

BOOST_AUTO_TEST_CASE(efficiency_bracketed_by_near_and_far_endpoints) {
    Material nai = make_NaI();

    EfficiencyCalculator calc_near;
    calc_near.set_fep_window_keV(kTestFepWindowKeV);
    calc_near.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_near.set_point_source(Eigen::Vector3d(0.0, 0.0, -8.0));
    auto res_near = calc_near.compute(precision_config(662.0));

    EfficiencyCalculator calc_far;
    calc_far.set_fep_window_keV(kTestFepWindowKeV);
    calc_far.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_far.set_point_source(Eigen::Vector3d(0.0, 0.0, -12.0));
    auto res_far = calc_far.compute(precision_config(662.0));

    EfficiencyCalculator calc_sph;
    calc_sph.set_fep_window_keV(kTestFepWindowKeV);
    calc_sph.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_sph.set_spherical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 1.0);
    auto res_sph = calc_sph.compute(precision_config(662.0));

    double margin = 0.15;
    BOOST_CHECK_LE(res_sph.total_efficiency, res_near.total_efficiency * (1.0 + margin));
    BOOST_CHECK_GE(res_sph.total_efficiency, res_far.total_efficiency * (1.0 - margin));
}

BOOST_AUTO_TEST_CASE(self_attenuating_sphere_below_bare) {
    // A self-attenuating (lead) sphere absorbs some of its own emission, so its
    // efficiency is below the same geometric sphere with no source material.
    Material nai = make_NaI();
    Material pb = make_Lead();
    Eigen::Vector3d center(0.0, 0.0, -8.0);

    EfficiencyCalculator calc_bare;
    calc_bare.set_fep_window_keV(kTestFepWindowKeV);
    calc_bare.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_bare.set_spherical_source(center, 3.0);
    auto res_bare = calc_bare.compute(precision_config(662.0));

    EfficiencyCalculator calc_self;
    calc_self.set_fep_window_keV(kTestFepWindowKeV);
    calc_self.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_self.set_spherical_source(center, 3.0);
    calc_self.set_source_material(&pb);
    auto res_self = calc_self.compute(precision_config(662.0));

    BOOST_TEST_MESSAGE("bare total=" << res_bare.total_efficiency << " +/- "
        << res_bare.total_uncertainty << "; self-atten total=" << res_self.total_efficiency
        << " +/- " << res_self.total_uncertainty);
    BOOST_CHECK_LT(res_self.total_efficiency, res_bare.total_efficiency);
    // The difference must be statistically significant.
    BOOST_CHECK_GT(zdiff(res_bare.total_efficiency, res_bare.total_uncertainty,
                         res_self.total_efficiency, res_self.total_uncertainty), 3.0);
}

BOOST_AUTO_TEST_CASE(spherical_shield_attenuates) {
    // Adding a spherical shield shell around the sphere lowers efficiency.
    Material nai = make_NaI();
    Material fe = make_Iron();
    Eigen::Vector3d center(0.0, 0.0, -8.0);

    EfficiencyCalculator calc0;
    calc0.set_fep_window_keV(kTestFepWindowKeV);
    calc0.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc0.set_spherical_source(center, 1.0);
    auto res0 = calc0.compute(precision_config(662.0));

    EfficiencyCalculator calc1;
    calc1.set_fep_window_keV(kTestFepWindowKeV);
    calc1.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc1.set_spherical_source(center, 1.0);
    calc1.add_source_shield(&fe, 1.0);   // 1 cm Fe shell
    auto res1 = calc1.compute(precision_config(662.0));

    BOOST_TEST_MESSAGE("no shield total=" << res0.total_efficiency
        << "; 1cm Fe total=" << res1.total_efficiency);
    BOOST_CHECK_LT(res1.total_efficiency, res0.total_efficiency);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Hollow shell / annular cylinder sanity
// ============================================================

BOOST_AUTO_TEST_SUITE(HollowAndAnnular)

BOOST_AUTO_TEST_CASE(hollow_shell_higher_than_solid_ball) {
    // Removing the deeply-buried (most self-absorbed) core activity raises the
    // per-photon efficiency: a void-center lead shell beats the solid lead ball
    // of the same outer radius and material.
    Material nai = make_NaI();
    Material pb = make_Lead();
    Eigen::Vector3d center(0.0, 0.0, -8.0);

    EfficiencyCalculator calc_solid;
    calc_solid.set_fep_window_keV(kTestFepWindowKeV);
    calc_solid.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_solid.set_spherical_source(center, 3.0);          // solid ball
    calc_solid.set_source_material(&pb);
    auto res_solid = calc_solid.compute(precision_config(662.0));

    EfficiencyCalculator calc_shell;
    calc_shell.set_fep_window_keV(kTestFepWindowKeV);
    calc_shell.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_shell.set_spherical_source(center, 3.0, Eigen::Matrix3d::Identity(), 2.0);  // shell [2,3]
    calc_shell.set_source_material(&pb);
    auto res_shell = calc_shell.compute(precision_config(662.0));

    BOOST_TEST_MESSAGE("solid ball total=" << res_solid.total_efficiency << " +/- "
        << res_solid.total_uncertainty << "; hollow shell total=" << res_shell.total_efficiency
        << " +/- " << res_shell.total_uncertainty);
    BOOST_CHECK_GT(res_shell.total_efficiency, res_solid.total_efficiency);
    BOOST_CHECK_GT(zdiff(res_shell.total_efficiency, res_shell.total_uncertainty,
                         res_solid.total_efficiency, res_solid.total_uncertainty), 3.0);
}

BOOST_AUTO_TEST_CASE(annular_cylinder_higher_than_solid_sideon) {
    // Side-on lead cylinder: the central bore is the deepest (most self-absorbed)
    // material in the viewing direction, so carving it out raises efficiency.
    Material nai = make_NaI();
    Material pb = make_Lead();
    Eigen::Vector3d center(0.0, 0.0, -8.0);
    Eigen::Matrix3d rot = side_on_rotation();

    EfficiencyCalculator calc_solid;
    calc_solid.set_fep_window_keV(kTestFepWindowKeV);
    calc_solid.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_solid.set_cylindrical_source(center, 3.0, 3.0, rot);   // solid
    calc_solid.set_source_material(&pb);
    auto res_solid = calc_solid.compute(precision_config(662.0));

    EfficiencyCalculator calc_tube;
    calc_tube.set_fep_window_keV(kTestFepWindowKeV);
    calc_tube.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_tube.set_cylindrical_source(center, 3.0, 3.0, rot, 2.0);  // bore r_in=2
    calc_tube.set_source_material(&pb);
    auto res_tube = calc_tube.compute(precision_config(662.0));

    BOOST_TEST_MESSAGE("solid cyl total=" << res_solid.total_efficiency << " +/- "
        << res_solid.total_uncertainty << "; annular total=" << res_tube.total_efficiency
        << " +/- " << res_tube.total_uncertainty);
    BOOST_CHECK_GT(res_tube.total_efficiency, res_solid.total_efficiency);
    BOOST_CHECK_GT(zdiff(res_tube.total_efficiency, res_tube.total_uncertainty,
                         res_solid.total_efficiency, res_solid.total_uncertainty), 3.0);
}

BOOST_AUTO_TEST_CASE(annular_cylinder_basic_invariants) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 2.0, 2.0,
                                Eigen::Matrix3d::Identity(), 1.0);  // bore r_in=1

    for (double E : {100.0, 662.0, 1332.0}) {
        auto res = calc.compute(E, 10000, 1);
        BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
        BOOST_CHECK_GT(res.total_efficiency, 0.0);
        BOOST_CHECK_LE(res.total_efficiency, 1.0);
    }
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  CylinderEndOn vs CylinderSideOn convention
// ============================================================

BOOST_AUTO_TEST_SUITE(CylinderOrientation)

BOOST_AUTO_TEST_CASE(endon_differs_from_sideon) {
    // A strongly anisotropic (long, thin) cylinder presents a very different
    // profile end-on vs side-on; the rotation must be threaded through, so the
    // two efficiencies must differ by many sigma.
    Material nai = make_NaI();
    Eigen::Vector3d center(0.0, 0.0, -10.0);

    EfficiencyCalculator calc_end;
    calc_end.set_fep_window_keV(kTestFepWindowKeV);
    calc_end.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_end.set_cylindrical_source(center, 0.5, 6.0);  // identity = EndOn
    auto res_end = calc_end.compute(precision_config(662.0));

    EfficiencyCalculator calc_side;
    calc_side.set_fep_window_keV(kTestFepWindowKeV);
    calc_side.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_side.set_cylindrical_source(center, 0.5, 6.0, side_on_rotation());
    auto res_side = calc_side.compute(precision_config(662.0));

    BOOST_TEST_MESSAGE("EndOn total=" << res_end.total_efficiency << " +/- "
        << res_end.total_uncertainty << "; SideOn total=" << res_side.total_efficiency
        << " +/- " << res_side.total_uncertainty);
    BOOST_CHECK_GT(zdiff(res_end.total_efficiency, res_end.total_uncertainty,
                         res_side.total_efficiency, res_side.total_uncertainty), 3.0);
    BOOST_CHECK_GT(res_end.total_efficiency, 0.0);
    BOOST_CHECK_GT(res_side.total_efficiency, 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Cascade summing on a spherical source (hand-built Co-60)
// ============================================================

BOOST_AUTO_TEST_SUITE(SphericalCascade)

BOOST_AUTO_TEST_CASE(co60_summing_out_on_sphere) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    // Close geometry maximises true-coincidence summing-out.
    calc.set_spherical_source(Eigen::Vector3d(0.0, 0.0, -2.0), 0.5);

    // Hand-built Co-60: two mutually coincident gammas (1173.2 + 1332.5 keV).
    DecayCascade co60;
    co60.branch_weight = 1.0;
    CascadeMember g1; g1.energy_keV = 1173.2; g1.intensity = 0.9985;
    CascadeMember g2; g2.energy_keV = 1332.5; g2.intensity = 0.9998;
    g1.coincident.push_back({1, 0.9998});  // partner index 1, P(g2|g1)
    g2.coincident.push_back({0, 0.9985});  // partner index 0, P(g1|g2)
    co60.members = {g1, g2};

    CascadeConfig cfg;
    cfg.cascades = {co60};
    // Explicit: this was written against the old 1.5 keV PeakWindow default, so
    //  spell it out rather than silently tracking whatever the default becomes.
    cfg.peaks = {PeakWindow{1173.2, kTestFepWindowKeV},
                 PeakWindow{1332.5, kTestFepWindowKeV}};
    cfg.num_events = 200000;
    cfg.method = CascadeMethod::Conditional;

    auto res = calc.compute_cascade(cfg);
    BOOST_REQUIRE_EQUAL(res.peaks.size(), 2u);
    for (const auto& p : res.peaks) {
        BOOST_TEST_MESSAGE("peak " << p.energy_keV << " eff_no_sum=" << p.eff_no_summing
            << " eff_sum=" << p.eff_with_summing << " factor=" << p.summing_factor
            << " +/- " << p.summing_factor_unc);
        BOOST_CHECK(p.found);
        BOOST_CHECK_GT(p.eff_no_summing, 0.0);
        BOOST_CHECK(std::isfinite(p.summing_factor));
        // Co-60: the partner gamma sums out the peak, so factor < 1.
        BOOST_CHECK_LT(p.summing_factor, 1.0);
        BOOST_CHECK_GT(p.summing_factor, 0.0);
    }
}

BOOST_AUTO_TEST_SUITE_END()
