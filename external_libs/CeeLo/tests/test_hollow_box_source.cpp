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

#define BOOST_TEST_MODULE HollowBoxSourceTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>
#include <random>
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

/// Source-material path length from trace_source_segments (0 if no segment).
double material_path(const SourceGeometry& sg, const Material* mat,
                     const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) {
    double path = 0.0;
    for (const auto& seg : sg.trace_source_segments(pos, dir, 662.0))
        if (seg.material == mat) path += seg.length;
    return path;
}

/// 90° rotation about y: source local z ← detector x (the SideOn convention).
Eigen::Matrix3d side_on_rotation() {
    Eigen::Matrix3d R;
    R << 0, 0, 1,
         0, 1, 0,
        -1, 0, 0;
    return R;
}

} // anonymous namespace

// ============================================================
//  Hollow box shell — analytic ray paths (MC-free)
// ============================================================

BOOST_AUTO_TEST_SUITE(HollowBoxPaths)

BOOST_AUTO_TEST_CASE(axis_aligned_paths_through_shell) {
    // Outer half-dims (3,3,3), inner void (2,2,2), centered at origin.
    Material pb = make_Lead();
    SourceGeometry sg;
    sg.configure_rectangular(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(3, 3, 3),
                             Eigen::Matrix3d::Identity(), Eigen::Vector3d(2, 2, 2));
    sg.set_source_material(&pb);

    // Emission in the +x wall, heading outward: only the remaining wall.
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {2.5, 0, 0}, {1, 0, 0}), 0.5, 1e-9);

    // Same point heading inward: near wall remnant (0.5), void crossing
    // (free), then the full far wall (1.0).
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {2.5, 0, 0}, {-1, 0, 0}), 1.5, 1e-9);

    // Query from the void center (e.g. a scattered photon): one wall thickness.
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {0, 0, 0}, {0, 0, 1}), 1.0, 1e-9);

    // Ray that never touches the inner box (|x| = 2.5 > 2 all along): the
    // full outer chord, exactly the solid-box path.
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {2.5, 2.5, 0}, {0, 0, -1}), 3.0, 1e-9);

    // compute_transmission must agree with the analytic path.
    double mu = pb.mu_total(0.662);
    double trans = sg.compute_transmission({2.5, 0, 0}, {-1, 0, 0}, 0.662);
    BOOST_CHECK_CLOSE(trans, std::exp(-mu * 1.5), 1e-6);
}

BOOST_AUTO_TEST_CASE(oblique_path_through_void) {
    Material pb = make_Lead();
    SourceGeometry sg;
    sg.configure_rectangular(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(3, 3, 3),
                             Eigen::Matrix3d::Identity(), Eigen::Vector3d(2, 2, 2));
    sg.set_source_material(&pb);

    // From the void center along (1,1,0)/sqrt(2): exits the inner box where
    // x = y = 2 (t = 2*sqrt(2)) and the outer box where x = y = 3
    // (t = 3*sqrt(2)); material path = sqrt(2).
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    Eigen::Vector3d dir(inv_sqrt2, inv_sqrt2, 0.0);
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {0, 0, 0}, dir), std::sqrt(2.0), 1e-9);

    // Oblique ray from within the wall crossing the void: from (2.5, 0, 0)
    // along (-1, 0.2, 0)/norm — enters the void at x = 2, exits at x = -2,
    // then crosses the far wall to x = -3. Path = (0.5 + 1.0) * norm/|dx|.
    Eigen::Vector3d d2(-1.0, 0.2, 0.0);
    double stretch = d2.norm() / 1.0;  // per unit x
    d2.normalize();
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {2.5, 0, 0}, d2), 1.5 * stretch, 1e-9);
}

BOOST_AUTO_TEST_CASE(rotated_shell_paths) {
    // Asymmetric hollow shell rotated side-on: the rotation must carry the
    // inner void with the outer box.
    // Per-axis wall thicknesses all differ (x: 0.3, y: 0.5, z: 0.6) so a
    // dropped/misapplied rotation would be caught.
    Material pb = make_Lead();
    Eigen::Matrix3d rot = side_on_rotation();
    SourceGeometry sg;
    sg.configure_rectangular(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 2, 3),
                             rot, Eigen::Vector3d(0.7, 1.5, 2.4));
    sg.set_source_material(&pb);

    // Detector-frame +z maps to source-local +x: wall [0.7, 1] → 0.3.
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {0, 0, 0}, {0, 0, 1}), 0.3, 1e-9);
    // Detector-frame +x maps to source-local -z: wall [2.4, 3] → 0.6.
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {0, 0, 0}, {1, 0, 0}), 0.6, 1e-9);
    // Detector-frame +y is unrotated (local y): wall [1.5, 2] → 0.5.
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {0, 0, 0}, {0, 1, 0}), 0.5, 1e-9);
}

BOOST_AUTO_TEST_CASE(solid_box_unchanged) {
    // inner = 0 must reduce exactly to the pre-hollow solid-box behavior,
    // whether defaulted or passed explicitly.
    Material pb = make_Lead();

    SourceGeometry sg_default;
    sg_default.configure_rectangular(Eigen::Vector3d(0, 0, 0),
                                     Eigen::Vector3d(3, 3, 3),
                                     Eigen::Matrix3d::Identity());
    sg_default.set_source_material(&pb);

    SourceGeometry sg_zero;
    sg_zero.configure_rectangular(Eigen::Vector3d(0, 0, 0),
                                  Eigen::Vector3d(3, 3, 3),
                                  Eigen::Matrix3d::Identity(),
                                  Eigen::Vector3d::Zero());
    sg_zero.set_source_material(&pb);

    const Eigen::Vector3d pos(1.0, -2.0, 0.5);
    const Eigen::Vector3d dir = Eigen::Vector3d(0.3, -0.2, 0.9).normalized();
    double p_default = material_path(sg_default, &pb, pos, dir);
    double p_zero = material_path(sg_zero, &pb, pos, dir);
    BOOST_CHECK_EQUAL(p_default, p_zero);

    // Hand value: from inside, path = distance to box exit along dir.
    double t_exit = (3.0 - pos.z()) / dir.z();  // z-face is the nearest exit
    t_exit = std::min(t_exit, (3.0 - pos.x()) / dir.x());
    t_exit = std::min(t_exit, (-3.0 - pos.y()) / dir.y());
    BOOST_CHECK_CLOSE(p_default, t_exit, 1e-9);
}

BOOST_AUTO_TEST_CASE(min_distance_to_boundary_inner_wall) {
    // With include_shields == false the inner-void wall is a material
    // boundary; the conservative per-axis candidate must kick in.
    SourceGeometry sg;
    sg.configure_rectangular(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(3, 3, 3),
                             Eigen::Matrix3d::Identity(), Eigen::Vector3d(2, 2, 2));

    // Deep in the +x wall, closer to the void than the outer surface.
    BOOST_CHECK_CLOSE(sg.min_distance_to_boundary({2.1, 0, 0}, false), 0.1, 1e-9);
    // Near the outer surface: outer face wins.
    BOOST_CHECK_CLOSE(sg.min_distance_to_boundary({2.9, 0, 0}, false), 0.1, 1e-9);
    // Mid-wall: equidistant.
    BOOST_CHECK_CLOSE(sg.min_distance_to_boundary({2.5, 0, 0}, false), 0.5, 1e-9);
    // Corner region (outside the inner box on all axes).
    BOOST_CHECK_CLOSE(sg.min_distance_to_boundary({2.5, 2.5, 2.5}, false), 0.5, 1e-9);
    // The outer envelope (include_shields == true) ignores the void.
    BOOST_CHECK_CLOSE(sg.min_distance_to_boundary({2.1, 0, 0}, true), 0.9, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Hollow box shell — position sampling (MC-free, fixed seed)
// ============================================================

BOOST_AUTO_TEST_SUITE(HollowBoxSampling)

BOOST_AUTO_TEST_CASE(samples_uniform_in_shell) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    const Eigen::Vector3d center(0.0, 0.0, -10.0);
    const Eigen::Vector3d outer(3.0, 3.0, 3.0);
    const Eigen::Vector3d inner(2.0, 2.0, 2.0);
    calc.set_rectangular_source(center, outer, Eigen::Matrix3d::Identity(), inner);

    std::mt19937_64 rng(987654321ULL);
    const int N = 200000;
    int n_front_slab = 0;   // local z in (2, 3]
    int n_x_wall = 0;       // |local x| > 2
    for (int i = 0; i < N; ++i) {
        const Eigen::Vector3d local = calc.sample_source_position_for_test(rng) - center;
        // Inside the outer box…
        BOOST_REQUIRE_LE(std::abs(local.x()), outer.x() + 1e-12);
        BOOST_REQUIRE_LE(std::abs(local.y()), outer.y() + 1e-12);
        BOOST_REQUIRE_LE(std::abs(local.z()), outer.z() + 1e-12);
        // …and never strictly inside the inner void.
        const bool in_void = std::abs(local.x()) < inner.x()
                          && std::abs(local.y()) < inner.y()
                          && std::abs(local.z()) < inner.z();
        BOOST_REQUIRE(!in_void);
        if (local.z() > inner.z()) ++n_front_slab;
        if (std::abs(local.x()) > inner.x()) ++n_x_wall;
    }

    // Uniform-in-shell volume fractions. V_shell = 6^3 - 4^3 = 152.
    // Front slab (z > 2): 6*6*1 = 36. |x| > 2 walls: 2*1*6*6 = 72.
    const double v_shell = 216.0 - 64.0;
    auto check_frac = [&](int n, double vol, const char* what) {
        const double p = vol / v_shell;
        const double frac = double(n) / N;
        const double sigma = std::sqrt(p * (1.0 - p) / N);
        BOOST_TEST_MESSAGE(what << ": frac=" << frac << " expected=" << p
                                << " sigma=" << sigma);
        BOOST_CHECK_LT(std::abs(frac - p), 5.0 * sigma);
    };
    check_frac(n_front_slab, 36.0, "front slab");
    check_frac(n_x_wall, 72.0, "x walls");
}

BOOST_AUTO_TEST_CASE(exponential_depth_respects_void) {
    // Exponential depth distribution + hollow shell: every sample still lands
    // in the shell (the rejection loop resamples the full proposal).
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    const Eigen::Vector3d center(0.0, 0.0, -10.0);
    const Eigen::Vector3d inner(1.5, 1.5, 1.5);
    calc.set_rectangular_source(center, Eigen::Vector3d(2.0, 2.0, 2.0),
                                Eigen::Matrix3d::Identity(), inner);
    calc.set_exponential_depth_distribution(1.0);

    std::mt19937_64 rng(24680);
    for (int i = 0; i < 20000; ++i) {
        const Eigen::Vector3d local = calc.sample_source_position_for_test(rng) - center;
        const bool in_void = std::abs(local.x()) < inner.x()
                          && std::abs(local.y()) < inner.y()
                          && std::abs(local.z()) < inner.z();
        BOOST_REQUIRE(!in_void);
    }
}

BOOST_AUTO_TEST_CASE(solid_sampling_identical_default_vs_zero) {
    // Solid boxes: defaulted inner and explicit zero inner must consume the
    // rng identically and produce bit-equal positions.
    Material nai = make_NaI();
    const Eigen::Vector3d center(0.0, 0.0, -10.0);

    EfficiencyCalculator calc_a;
    calc_a.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_a.set_rectangular_source(center, Eigen::Vector3d(3.0, 2.0, 1.0));

    EfficiencyCalculator calc_b;
    calc_b.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_b.set_rectangular_source(center, Eigen::Vector3d(3.0, 2.0, 1.0),
                                  Eigen::Matrix3d::Identity(),
                                  Eigen::Vector3d::Zero());

    std::mt19937_64 rng_a(1357), rng_b(1357);
    for (int i = 0; i < 1000; ++i) {
        const Eigen::Vector3d pa = calc_a.sample_source_position_for_test(rng_a);
        const Eigen::Vector3d pb = calc_b.sample_source_position_for_test(rng_b);
        BOOST_REQUIRE_EQUAL(pa.x(), pb.x());
        BOOST_REQUIRE_EQUAL(pa.y(), pb.y());
        BOOST_REQUIRE_EQUAL(pa.z(), pb.z());
    }
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Hollow box shell — MC sanity
// ============================================================

BOOST_AUTO_TEST_SUITE(HollowBoxMC)

BOOST_AUTO_TEST_CASE(basic_invariants) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -10.0),
                                Eigen::Vector3d(2.0, 2.0, 2.0),
                                Eigen::Matrix3d::Identity(),
                                Eigen::Vector3d(1.5, 1.5, 1.5));

    for (double E : {100.0, 662.0, 1332.0}) {
        auto res = calc.compute(E, 10000, 1);
        BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
        BOOST_CHECK_GT(res.total_efficiency, 0.0);
        BOOST_CHECK_LE(res.total_efficiency, 1.0);
    }
}

#ifdef CEELO_RUN_MC_TESTS
BOOST_AUTO_TEST_CASE(hollow_box_higher_than_solid_self_attenuating) {
    // Removing the deeply-buried (most self-absorbed) core activity raises
    // the per-photon efficiency — mirrors hollow_shell_higher_than_solid_ball.
    Material nai = make_NaI();
    Material pb = make_Lead();
    Eigen::Vector3d center(0.0, 0.0, -8.0);

    EfficiencyCalculator calc_solid;
    calc_solid.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_solid.set_rectangular_source(center, Eigen::Vector3d(3.0, 3.0, 3.0));
    calc_solid.set_source_material(&pb);
    auto res_solid = calc_solid.compute(precision_config(662.0));

    EfficiencyCalculator calc_shell;
    calc_shell.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_shell.set_rectangular_source(center, Eigen::Vector3d(3.0, 3.0, 3.0),
                                      Eigen::Matrix3d::Identity(),
                                      Eigen::Vector3d(2.0, 2.0, 2.0));
    calc_shell.set_source_material(&pb);
    auto res_shell = calc_shell.compute(precision_config(662.0));

    BOOST_TEST_MESSAGE("solid box total=" << res_solid.total_efficiency << " +/- "
        << res_solid.total_uncertainty << "; hollow shell total="
        << res_shell.total_efficiency << " +/- " << res_shell.total_uncertainty);
    BOOST_CHECK_GT(res_shell.total_efficiency, res_solid.total_efficiency);
    BOOST_CHECK_GT(zdiff(res_shell.total_efficiency, res_shell.total_uncertainty,
                         res_solid.total_efficiency, res_solid.total_uncertainty), 3.0);
}

BOOST_AUTO_TEST_CASE(biasing_consistency_hollow_box) {
    // Auto-biased (cone/two-stream) vs plain isotropic must agree — the
    // hollow-box arbiter for the worst-case-cone candidate handling (the
    // inner void adds no candidates, matching the cylinder-bore treatment).
    Material nai = make_NaI();
    Material wood = make_Cellulose();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -12.0),
                                Eigen::Vector3d(4.0, 4.0, 4.0),
                                Eigen::Matrix3d::Identity(),
                                Eigen::Vector3d(3.5, 3.5, 3.5));
    calc.set_source_material(&wood);

    auto res_biased = calc.compute(precision_config(662.0, 0.01));

    EfficiencyCalculator calc_iso;
    calc_iso.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_iso.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -12.0),
                                    Eigen::Vector3d(4.0, 4.0, 4.0),
                                    Eigen::Matrix3d::Identity(),
                                    Eigen::Vector3d(3.5, 3.5, 3.5));
    calc_iso.set_source_material(&wood);
    calc_iso.enable_cone_sampling(false);
    auto res_iso = calc_iso.compute(precision_config(662.0, 0.01));

    BOOST_TEST_MESSAGE("biased total=" << res_biased.total_efficiency << " +/- "
        << res_biased.total_uncertainty << "; isotropic total="
        << res_iso.total_efficiency << " +/- " << res_iso.total_uncertainty);
    BOOST_CHECK_LT(zdiff(res_biased.total_efficiency, res_biased.total_uncertainty,
                         res_iso.total_efficiency, res_iso.total_uncertainty), 4.0);
}
#endif // CEELO_RUN_MC_TESTS

BOOST_AUTO_TEST_SUITE_END()
