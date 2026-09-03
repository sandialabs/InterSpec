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

#define BOOST_TEST_MODULE PipeSourceTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>
#include <vector>

using namespace ceelo;

// Pipe modeling conventions verified here (see also stage-B spec):
//   - pipe CONTENTS = set_cylindrical_source(center, r_contents, half_len,
//     rot) + set_source_material(contents) + add_source_shield(steel, t_wall,
//     t_end); t_end = 0 models an open (uncapped) pipe.
//   - pipe WALL = set_cylindrical_source(..., inner_radius = r_bore) +
//     set_source_material(steel).
//   - KNOWN GAP (deferred): the annular-source bore is a NON-attenuating
//     void, so "wall source with attenuating contents" (e.g. steel pipe wall
//     activity with water inside) is not expressible yet.
//
// The (t_radial, t_end) shield is built in the SOURCE-LOCAL frame
// (SourceGeometry transforms position/direction by cyl_rotation_ before
// intersecting the shield tubes), so it composes with a rotated (side-on)
// source; the analytic checks below pin that down.

namespace {

/// SimulationConfig targeting `target` FEP relative precision.
SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.target_total_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 5000;
    config.termination.max_wall_seconds = 30.0;
    config.num_threads = 0;
    config.batch_size = 10000;
    return config;
}

/// z-score for the difference of two independent efficiencies.
double zdiff(double a, double sa, double b, double sb) {
    double s = std::sqrt(sa * sa + sb * sb);
    return (s > 0.0) ? std::abs(a - b) / s : 0.0;
}

/// Total path through `mat` from trace_source_segments.
double material_path(const SourceGeometry& sg, const Material* mat,
                     const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) {
    double path = 0.0;
    for (const auto& seg : sg.trace_source_segments(pos, dir, 662.0))
        if (seg.material == mat) path += seg.length;
    return path;
}

/// SideOn rotation: cylinder local z ← detector x (90° about y).
Eigen::Matrix3d side_on_rotation() {
    Eigen::Matrix3d R;
    R << 0, 0, 1,
         0, 1, 0,
        -1, 0, 0;
    return R;
}

} // anonymous namespace

// ============================================================
//  Pipe contents + wall shield: rotation composition (MC-free)
// ============================================================

BOOST_AUTO_TEST_SUITE(PipeShieldRotation)

BOOST_AUTO_TEST_CASE(side_on_open_pipe_paths) {
    // Side-on water-filled steel pipe: contents r = 2.1, wall t = 0.4
    // (r_outer = 2.5), half-length 15, open ends (t_end = 0), axis along
    // detector x, center 10 cm in front of the detector.
    Material water = make_Water();
    Material steel = make_StainlessSteel304();
    const Eigen::Vector3d center(0.0, 0.0, -10.0);

    SourceGeometry sg;
    sg.configure_cylindrical(center, 2.1, 15.0, side_on_rotation());
    sg.set_source_material(&water);
    sg.add_shield(&steel, 0.4, 0.0);  // side wall only — open pipe

    // From the pipe center toward the detector (+z detector = radial in the
    // pipe local frame): water radius, then one wall thickness of steel.
    BOOST_CHECK_CLOSE(material_path(sg, &water, center, {0, 0, 1}), 2.1, 1e-9);
    BOOST_CHECK_CLOSE(material_path(sg, &steel, center, {0, 0, 1}), 0.4, 1e-9);

    // Along the pipe axis (detector x = pipe local z): water half-length,
    // and NO steel — the open end has no cap.
    BOOST_CHECK_CLOSE(material_path(sg, &water, center, {1, 0, 0}), 15.0, 1e-9);
    BOOST_CHECK_EQUAL(material_path(sg, &steel, center, {1, 0, 0}), 0.0);

    // Transmission agrees with the analytic radial path.
    const double mu_w = water.mu_total(0.662);
    const double mu_s = steel.mu_total(0.662);
    BOOST_CHECK_CLOSE(sg.compute_transmission(center, {0, 0, 1}, 0.662),
                      std::exp(-mu_w * 2.1 - mu_s * 0.4), 1e-6);
}

BOOST_AUTO_TEST_CASE(side_on_capped_pipe_end_path) {
    // Capped pipe (t_end = t_wall): the axial ray now crosses one end cap.
    Material water = make_Water();
    Material steel = make_StainlessSteel304();
    const Eigen::Vector3d center(0.0, 0.0, -10.0);

    SourceGeometry sg;
    sg.configure_cylindrical(center, 2.1, 15.0, side_on_rotation());
    sg.set_source_material(&water);
    sg.add_shield(&steel, 0.4, 0.4);

    BOOST_CHECK_CLOSE(material_path(sg, &water, center, {1, 0, 0}), 15.0, 1e-9);
    BOOST_CHECK_CLOSE(material_path(sg, &steel, center, {1, 0, 0}), 0.4, 1e-9);
    // Radial path unchanged by the caps.
    BOOST_CHECK_CLOSE(material_path(sg, &steel, center, {0, 0, 1}), 0.4, 1e-9);
}

BOOST_AUTO_TEST_CASE(side_on_off_center_emission) {
    // Emission near one end of the side-on pipe (local (0, 0, +10) → detector
    // frame (-10, 0, -10)): the ray toward the detector is still radial in
    // the local frame, so it sees the same water radius + wall thickness.
    // This pins the shield to the source local frame at a non-center point.
    Material water = make_Water();
    Material steel = make_StainlessSteel304();
    const Eigen::Vector3d center(0.0, 0.0, -10.0);
    const Eigen::Matrix3d rot = side_on_rotation();

    SourceGeometry sg;
    sg.configure_cylindrical(center, 2.1, 15.0, rot);
    sg.set_source_material(&water);
    sg.add_shield(&steel, 0.4, 0.0);

    const Eigen::Vector3d pos = rot.transpose() * Eigen::Vector3d(0, 0, 10) + center;
    BOOST_CHECK_CLOSE(material_path(sg, &water, pos, {0, 0, 1}), 2.1, 1e-9);
    BOOST_CHECK_CLOSE(material_path(sg, &steel, pos, {0, 0, 1}), 0.4, 1e-9);
    // Along the axis toward the near open end: 5 cm of water, no steel.
    const Eigen::Vector3d axis_det = rot.transpose() * Eigen::Vector3d(0, 0, 1);
    BOOST_CHECK_CLOSE(material_path(sg, &water, pos, axis_det), 5.0, 1e-9);
    BOOST_CHECK_EQUAL(material_path(sg, &steel, pos, axis_det), 0.0);
}

BOOST_AUTO_TEST_CASE(end_on_matches_side_on_in_local_frame) {
    // The same pipe end-on (identity rotation): detector z is now the pipe
    // axis. The local-frame paths must equal the side-on case with the
    // detector-frame directions permuted accordingly.
    Material water = make_Water();
    Material steel = make_StainlessSteel304();
    const Eigen::Vector3d center(0.0, 0.0, -25.0);

    SourceGeometry sg;
    sg.configure_cylindrical(center, 2.1, 15.0, Eigen::Matrix3d::Identity());
    sg.set_source_material(&water);
    sg.add_shield(&steel, 0.4, 0.0);

    // Toward the detector = axial: 15 cm water, open end.
    BOOST_CHECK_CLOSE(material_path(sg, &water, center, {0, 0, 1}), 15.0, 1e-9);
    BOOST_CHECK_EQUAL(material_path(sg, &steel, center, {0, 0, 1}), 0.0);
    // Transverse = radial: water radius + wall.
    BOOST_CHECK_CLOSE(material_path(sg, &water, center, {1, 0, 0}), 2.1, 1e-9);
    BOOST_CHECK_CLOSE(material_path(sg, &steel, center, {1, 0, 0}), 0.4, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Pipe wall (annular source) paths (MC-free)
// ============================================================

BOOST_AUTO_TEST_SUITE(PipeWallSource)

BOOST_AUTO_TEST_CASE(side_on_wall_source_paths) {
    // Pipe WALL as the active region: annular cylinder [2.1, 2.5], steel
    // source material, side-on. The bore is a non-attenuating void.
    // KNOWN GAP (deferred): a wall source with ATTENUATING contents (e.g.
    // water inside) is not expressible — the bore is always empty.
    Material steel = make_StainlessSteel304();
    const Eigen::Vector3d center(0.0, 0.0, -10.0);
    const Eigen::Matrix3d rot = side_on_rotation();

    SourceGeometry sg;
    sg.configure_cylindrical(center, 2.5, 15.0, rot, 2.1);
    sg.set_source_material(&steel);

    // Mid-wall emission point: local (2.3, 0, 0) → detector (0, 0, -7.7).
    const Eigen::Vector3d pos = rot.transpose() * Eigen::Vector3d(2.3, 0, 0) + center;
    BOOST_CHECK_CLOSE(pos.z(), -7.7, 1e-9);

    // Outward toward the detector: remaining near-wall steel only.
    BOOST_CHECK_CLOSE(material_path(sg, &steel, pos, {0, 0, 1}), 0.2, 1e-9);
    // Inward: near-wall remnant (0.2), free bore crossing, far wall (0.4).
    BOOST_CHECK_CLOSE(material_path(sg, &steel, pos, {0, 0, -1}), 0.6, 1e-9);

    const double mu_s = steel.mu_total(0.662);
    BOOST_CHECK_CLOSE(sg.compute_transmission(pos, {0, 0, -1}, 0.662),
                      std::exp(-mu_s * 0.6), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Pipe MC sanity
// ============================================================

BOOST_AUTO_TEST_SUITE(PipeMC)

BOOST_AUTO_TEST_CASE(side_on_contents_basic_invariants) {
    Material nai = make_NaI();
    Material water = make_Water();
    Material steel = make_StainlessSteel304();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 2.1, 15.0,
                                side_on_rotation());
    calc.set_source_material(&water);
    calc.add_source_shield(&steel, 0.4, 0.0);

    for (double E : {186.0, 662.0, 1332.0}) {
        auto res = calc.compute(E, 10000, 1);
        BOOST_CHECK_GE(res.total_efficiency, res.full_energy_peak_efficiency - 1e-9);
        BOOST_CHECK_GT(res.total_efficiency, 0.0);
        BOOST_CHECK_LE(res.total_efficiency, 1.0);
    }
}

#ifdef CEELO_RUN_MC_TESTS
BOOST_AUTO_TEST_CASE(end_on_open_pipe_higher_than_capped) {
    // End-on pipe: photons toward the detector leave through the open end,
    // so capping it (t_end > 0) must lower the efficiency significantly.
    Material nai = make_NaI();
    Material water = make_Water();
    Material steel = make_StainlessSteel304();
    const Eigen::Vector3d center(0.0, 0.0, -20.0);

    EfficiencyCalculator calc_open;
    calc_open.set_fep_window_keV(kTestFepWindowKeV);
    calc_open.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_open.set_cylindrical_source(center, 2.1, 10.0);
    calc_open.set_source_material(&water);
    calc_open.add_source_shield(&steel, 0.4, 0.0);
    auto res_open = calc_open.compute(precision_config(186.0, 0.01));

    EfficiencyCalculator calc_capped;
    calc_capped.set_fep_window_keV(kTestFepWindowKeV);
    calc_capped.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_capped.set_cylindrical_source(center, 2.1, 10.0);
    calc_capped.set_source_material(&water);
    calc_capped.add_source_shield(&steel, 0.4, 0.4);
    auto res_capped = calc_capped.compute(precision_config(186.0, 0.01));

    BOOST_TEST_MESSAGE("open total=" << res_open.total_efficiency << " +/- "
        << res_open.total_uncertainty << "; capped total="
        << res_capped.total_efficiency << " +/- " << res_capped.total_uncertainty);
    BOOST_CHECK_GT(res_open.total_efficiency, res_capped.total_efficiency);
    BOOST_CHECK_GT(zdiff(res_open.total_efficiency, res_open.total_uncertainty,
                         res_capped.total_efficiency, res_capped.total_uncertainty), 3.0);
}
#endif // CEELO_RUN_MC_TESTS

BOOST_AUTO_TEST_SUITE_END()
