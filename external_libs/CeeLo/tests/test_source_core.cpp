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

#define BOOST_TEST_MODULE SourceCoreTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include "physics/ElectronCsda.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using namespace ceelo;

namespace {

/// Path through one material, summed over the traced segments.
double material_path(const SourceGeometry& sg, const Material* mat,
                     const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) {
    double path = 0.0;
    for (const auto& seg : sg.trace_source_segments(pos, dir, 662.0))
        if (seg.material == mat) path += seg.length;
    return path;
}

/// Geometric path through everything the trace returned, voids included.
double traced_extent(const SourceGeometry& sg,
                     const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) {
    double path = 0.0;
    for (const auto& seg : sg.trace_source_segments(pos, dir, 662.0))
        path += seg.length;
    return path;
}

/// A representative fan of rays, so a test is not passing on one lucky chord.
std::vector<Eigen::Vector3d> ray_fan() {
    std::vector<Eigen::Vector3d> dirs;
    for (double a = 0.0; a < 6.2; a += 0.37)
        for (double b = -1.0; b <= 1.0; b += 0.5)
            dirs.push_back(Eigen::Vector3d(std::cos(a), std::sin(a), b).normalized());
    return dirs;
}

} // anonymous namespace


// ============================================================
//  Analytic ray paths (MC-free)
// ============================================================

BOOST_AUTO_TEST_SUITE(SourceCorePaths)

BOOST_AUTO_TEST_CASE(core_of_same_material_is_a_solid_source) {
    // The invariant that pins the whole feature: a shell whose core is its OWN
    // material must be indistinguishable, ray by ray, from a solid source. Note
    // this is a statement about TRANSMISSION only - the samplers still differ,
    // since a shell source emits from the shell alone.
    Material pb = make_Lead();

    SourceGeometry solid;
    solid.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 0.0,
                              Eigen::Matrix3d::Identity());
    solid.set_source_material(&pb);

    SourceGeometry cored;
    cored.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                              Eigen::Matrix3d::Identity());
    cored.set_source_material(&pb);
    cored.add_core(&pb, 2.0);

    const std::vector<Eigen::Vector3d> starts = {
        {2.5, 0.0, 0.0}, {0.0, 2.2, 0.0}, {1.5, 1.5, 1.5}, {-2.9, 0.1, 0.0}
    };
    for (const Eigen::Vector3d& pos : starts) {
        for (const Eigen::Vector3d& dir : ray_fan()) {
            BOOST_CHECK_CLOSE(material_path(cored, &pb, pos, dir),
                              material_path(solid, &pb, pos, dir), 1e-8);
            BOOST_CHECK_CLOSE(cored.compute_transmission(pos, dir, 0.662),
                              solid.compute_transmission(pos, dir, 0.662), 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(no_core_leaves_the_hollow_path_untouched) {
    // Adding the feature must not move a source that does not use it: with no
    // add_core() call the original per-layer loops run, and a hollow shell still
    // sees a free (non-attenuating) centre.
    Material pb = make_Lead();
    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                           Eigen::Matrix3d::Identity());
    sg.set_source_material(&pb);

    // Straight through the middle: near wall 0.5 + far wall 1.0, void free.
    BOOST_CHECK_CLOSE(material_path(sg, &pb, {2.5, 0, 0}, {-1, 0, 0}), 1.5, 1e-9);
    BOOST_CHECK_CLOSE(sg.compute_transmission({2.5, 0, 0}, {-1, 0, 0}, 0.662),
                      std::exp(-pb.mu_total(0.662) * 1.5), 1e-6);
}

BOOST_AUTO_TEST_CASE(segments_are_ordered_and_geometrically_complete) {
    // The reason cores need their own tracer: transport locates an interaction
    // by accumulating segment lengths along the ray, so the list must be in
    // traversal order AND must account for every centimetre of the chord. The
    // pre-existing merged near+far shell path satisfies neither.
    Material pb = make_Lead();
    Material fe = make_Iron();
    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                           Eigen::Matrix3d::Identity());
    sg.set_source_material(&pb);
    sg.add_core(&fe, 2.0);

    const auto segs = sg.trace_source_segments({2.5, 0, 0}, {-1, 0, 0}, 662.0);
    BOOST_REQUIRE_EQUAL(segs.size(), 3u);
    BOOST_CHECK_EQUAL(segs[0].material, &pb);   //near wall
    BOOST_CHECK_EQUAL(segs[1].material, &fe);   //core
    BOOST_CHECK_EQUAL(segs[2].material, &pb);   //far wall
    BOOST_CHECK_CLOSE(segs[0].length, 0.5, 1e-9);
    BOOST_CHECK_CLOSE(segs[1].length, 4.0, 1e-9);
    BOOST_CHECK_CLOSE(segs[2].length, 1.0, 1e-9);

    // Nothing missing: the traced lengths sum to the whole chord to the surface.
    BOOST_CHECK_CLOSE(traced_extent(sg, {2.5, 0, 0}, {-1, 0, 0}), 5.5, 1e-9);
}

BOOST_AUTO_TEST_CASE(a_dense_core_attenuates_exactly_its_chord) {
    // A void core and an iron core differ by exactly exp(-mu_Fe * core chord).
    Material pb = make_Lead();
    Material fe = make_Iron();

    SourceGeometry voided;
    voided.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                               Eigen::Matrix3d::Identity());
    voided.set_source_material(&pb);

    SourceGeometry cored;
    cored.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                              Eigen::Matrix3d::Identity());
    cored.set_source_material(&pb);
    cored.add_core(&fe, 2.0);

    const Eigen::Vector3d pos(2.5, 0, 0), dir(-1, 0, 0);
    const double core_chord = 4.0;
    BOOST_CHECK_CLOSE(cored.compute_transmission(pos, dir, 0.662),
                      voided.compute_transmission(pos, dir, 0.662)
                        * std::exp(-fe.mu_total(0.662) * core_chord), 1e-6);

    // And a ray that misses the core entirely is unaffected.
    const Eigen::Vector3d tangent_pos(0.0, 2.5, 0.0), tangent_dir(1, 0, 0);
    BOOST_CHECK_CLOSE(cored.compute_transmission(tangent_pos, tangent_dir, 0.662),
                      voided.compute_transmission(tangent_pos, tangent_dir, 0.662),
                      1e-9);
}

BOOST_AUTO_TEST_CASE(cores_are_additive) {
    // Four 0.5 cm iron cores must equal one 2.0 cm iron core.
    Material pb = make_Lead();
    Material fe = make_Iron();

    SourceGeometry one;
    one.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                            Eigen::Matrix3d::Identity());
    one.set_source_material(&pb);
    one.add_core(&fe, 2.0);

    SourceGeometry many;
    many.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                             Eigen::Matrix3d::Identity());
    many.set_source_material(&pb);
    for (int i = 0; i < 4; ++i) many.add_core(&fe, 0.5);

    for (const Eigen::Vector3d& dir : ray_fan()) {
        BOOST_CHECK_CLOSE(one.compute_transmission({2.5, 0, 0}, dir, 0.662),
                          many.compute_transmission({2.5, 0, 0}, dir, 0.662), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(partly_filled_core_leaves_a_void_but_keeps_the_distances) {
    // A core stack that does not reach the centre leaves a genuine void. The
    // void must not attenuate, and must still occupy its place in the segment
    // list so the along-ray distances stay right.
    Material pb = make_Lead();
    Material fe = make_Iron();
    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                           Eigen::Matrix3d::Identity());
    sg.set_source_material(&pb);
    sg.add_core(&fe, 0.5);   //iron shell [1.5, 2], void ball r < 1.5

    const auto segs = sg.trace_source_segments({2.5, 0, 0}, {-1, 0, 0}, 662.0);
    BOOST_REQUIRE_EQUAL(segs.size(), 5u);
    BOOST_CHECK_EQUAL(segs[0].material, &pb);        //near lead wall, 0.5
    BOOST_CHECK_EQUAL(segs[1].material, &fe);        //near iron shell, 0.5
    BOOST_CHECK_EQUAL(segs[2].material, nullptr);    //void, 3.0
    BOOST_CHECK_EQUAL(segs[3].material, &fe);        //far iron shell, 0.5
    BOOST_CHECK_EQUAL(segs[4].material, &pb);        //far lead wall, 1.0
    BOOST_CHECK_CLOSE(segs[2].length, 3.0, 1e-9);
    BOOST_CHECK_CLOSE(traced_extent(sg, {2.5, 0, 0}, {-1, 0, 0}), 5.5, 1e-9);

    BOOST_CHECK_CLOSE(material_path(sg, &fe, {2.5, 0, 0}, {-1, 0, 0}), 1.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(nested_cylinder_is_not_a_pipe) {
    // A closed inner cavity (nested cylinders, InterSpec's model) attenuates
    // differently from a through-bore of the same radius: the material fills the
    // full radius beyond the cavity ends.
    Material pb = make_Lead();

    SourceGeometry pipe;   //bore runs the whole length
    pipe.configure_cylindrical(Eigen::Vector3d::Zero(), 3.0, 5.0,
                               Eigen::Matrix3d::Identity(), 2.0);
    pipe.set_source_material(&pb);

    SourceGeometry nested;  //cavity only over |z| < 1
    nested.configure_cylindrical(Eigen::Vector3d::Zero(), 3.0, 5.0,
                                 Eigen::Matrix3d::Identity(), 2.0, 1.0);
    nested.set_source_material(&pb);

    // Radially inward at |z| = 3, i.e. past the cavity's end: the nested case is
    // solid there, the pipe is not.
    const Eigen::Vector3d pos(2.5, 0.0, 3.0), dir(-1, 0, 0);
    BOOST_CHECK_CLOSE(material_path(nested, &pb, pos, dir), 5.5, 1e-9);
    BOOST_CHECK_CLOSE(material_path(pipe, &pb, pos, dir), 1.5, 1e-9);

    // Through the cavity itself the two agree.
    const Eigen::Vector3d mid(2.5, 0.0, 0.0);
    BOOST_CHECK_CLOSE(material_path(nested, &pb, mid, dir),
                      material_path(pipe, &pb, mid, dir), 1e-9);
}

BOOST_AUTO_TEST_CASE(box_core_of_same_material_is_a_solid_box) {
    Material pb = make_Lead();

    SourceGeometry solid;
    solid.configure_rectangular(Eigen::Vector3d::Zero(), Eigen::Vector3d(3, 3, 3),
                                Eigen::Matrix3d::Identity());
    solid.set_source_material(&pb);

    SourceGeometry cored;
    cored.configure_rectangular(Eigen::Vector3d::Zero(), Eigen::Vector3d(3, 3, 3),
                                Eigen::Matrix3d::Identity(),
                                Eigen::Vector3d(2, 2, 2));
    cored.set_source_material(&pb);
    cored.add_core(&pb, 2.0, 2.0, 2.0);

    for (const Eigen::Vector3d& dir : ray_fan()) {
        BOOST_CHECK_CLOSE(material_path(cored, &pb, {2.5, 0, 0}, dir),
                          material_path(solid, &pb, {2.5, 0, 0}, dir), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(shields_still_wrap_a_cored_source) {
    // Cores fill inward, shields grow outward, and both must appear - in the
    // right order - in one trace.
    Material pb = make_Lead();
    Material fe = make_Iron();
    Material cu = make_Copper();

    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                           Eigen::Matrix3d::Identity());
    sg.set_source_material(&pb);
    sg.add_core(&fe, 2.0);
    sg.add_shield(&cu, 1.0);

    const auto segs = sg.trace_source_segments({2.5, 0, 0}, {-1, 0, 0}, 662.0);
    BOOST_REQUIRE_EQUAL(segs.size(), 4u);
    BOOST_CHECK_EQUAL(segs[0].material, &pb);
    BOOST_CHECK_EQUAL(segs[1].material, &fe);
    BOOST_CHECK_EQUAL(segs[2].material, &pb);
    BOOST_CHECK_EQUAL(segs[3].material, &cu);
    BOOST_CHECK_CLOSE(segs[3].length, 1.0, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Segment-cap semantics (the Moliere walk asks for one segment)
// ============================================================

BOOST_AUTO_TEST_SUITE(SourceCoreSegmentCap)

BOOST_AUTO_TEST_CASE(one_segment_request_gets_the_whole_material_run) {
    // A caller asking for a single segment wants the distance to the next real
    // INTERFACE.  When the core is the same material as the shell there is no
    // interface, so the answer must be the whole chord - not the first interval.
    Material pb = make_Lead();
    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                           Eigen::Matrix3d::Identity());
    sg.set_source_material(&pb);
    sg.add_core(&pb, 2.0);

    std::vector<SourceGeometry::SourcePathSegment> segs;
    sg.trace_source_segments({2.5, 0, 0}, {-1, 0, 0}, 662.0, segs, 1);
    BOOST_REQUIRE_EQUAL(segs.size(), 1u);
    BOOST_CHECK_EQUAL(segs[0].material, &pb);
    BOOST_CHECK_CLOSE(segs[0].length, 5.5, 1e-9);

    // With a different core material the cap does bite, at the real interface.
    Material fe = make_Iron();
    SourceGeometry sg2;
    sg2.configure_spherical(Eigen::Vector3d::Zero(), 3.0, 2.0,
                            Eigen::Matrix3d::Identity());
    sg2.set_source_material(&pb);
    sg2.add_core(&fe, 2.0);

    sg2.trace_source_segments({2.5, 0, 0}, {-1, 0, 0}, 662.0, segs, 1);
    BOOST_REQUIRE_EQUAL(segs.size(), 1u);
    BOOST_CHECK_CLOSE(segs[0].length, 0.5, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SourceCoreTransport)

// The analytic suite above pins the GEOMETRY.  These pin the two TRANSPORT
// paths a core reaches that the geometry tests cannot see, both of which were
// silently wrong until an adversarial review found them.

BOOST_AUTO_TEST_CASE(fep_only_self_core_matches_the_solid_it_equals) {
    // compute_transmission_fep_only() walks the stack with a REUSED segment
    // buffer and re-traces after every Rayleigh turn.  trace_cored_segments()
    // appends, so a missing clear() made each turn re-process the segments
    // already walked: attenuation applied twice, and the exit point ran far
    // outside the source (a 2 cm source reported >100 cm), which then fed a
    // bogus air gap downstream.
    //
    // A core of the shell's own material is the same solid, so the two must
    // agree.  60 keV iron is chosen because Rayleigh is strong there
    // (mu_rs * L ~ 1.5), which is what drives the retrace.
    Material fe = make_Iron();

    SourceGeometry cored, solid;
    cored.configure_spherical(Eigen::Vector3d(0, 0, 0), 1.0, 0.5, Eigen::Matrix3d::Identity());
    cored.set_source_material(&fe);
    cored.add_core(&fe, 0.5);
    solid.configure_spherical(Eigen::Vector3d(0, 0, 0), 1.0, 0.0, Eigen::Matrix3d::Identity());
    solid.set_source_material(&fe);

    std::mt19937_64 rng_a(12345), rng_b(12345);
    double sum_a = 0.0, sum_b = 0.0, max_exit_a = 0.0;
    const int N = 4000;
    for (int i = 0; i < N; ++i) {
        const Eigen::Vector3d dir(0.0, 0.0, 1.0);
        // NOTE keV vs MeV: this entry point takes MeV.  Passing 60.0 here means
        //  60 MeV, where Rayleigh is negligible, no retrace ever happens and the
        //  test is vacuous -- which is exactly how the first version of it passed
        //  against the un-fixed code.
        auto a = cored.compute_transmission_fep_only(Eigen::Vector3d(0, 0, -1.0), dir, 0.060, rng_a);
        auto b = solid.compute_transmission_fep_only(Eigen::Vector3d(0, 0, -1.0), dir, 0.060, rng_b);
        sum_a += a.weight;
        sum_b += b.weight;
        max_exit_a = std::max(max_exit_a, a.exit_position.norm());
    }
    // The exit point must stay on the source, not run away down the ray.
    BOOST_CHECK_LT(max_exit_a, 1.01);
    // Same optical medium, same chord => same mean weight.  The tolerance is
    // loose only because the cored path re-traces the post-scatter leg, which
    // the aggregated legacy path does not model.
    BOOST_CHECK_CLOSE(sum_a / N, sum_b / N, 12.0);
}

BOOST_AUTO_TEST_CASE(an_electron_crossing_the_cavity_has_not_escaped) {
    // A partly-filled core leaves a null-material cavity segment.  The Moliere
    // source walk read that as "outside the source geometry" and declared the
    // electron escaped -- at a point deep INSIDE the source, with full residual
    // energy, skipping the far wall and every shield beyond it.
    Material soil = make_Soil(), fe = make_Iron();

    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d(0, 0, 0), 2.05, 2.0, Eigen::Matrix3d::Identity());
    sg.set_source_material(&soil);
    sg.add_core(&fe, 0.05);            // fills [1.95, 2.0]; r < 1.95 stays void
    sg.set_source_electron_transport(true);

    // Born in the near wall heading inward.  This geometry is deliberately thin
    // enough that escaping is physically possible, so the invariant is not HOW
    // OFTEN an electron gets out but WHERE: the only way out is the outer
    // surface at r = 2.05.  Exiting at the cavity wall (r = 1.95) is the bug.
    std::mt19937_64 rng(999);
    int escaped = 0, escaped_inside = 0;
    const int N = 200;
    for (int i = 0; i < N; ++i) {
        auto w = ElectronCsda::instance().walk_in_source_geometry(
            sg, soil, Eigen::Vector3d(0, 0, -2.02), Eigen::Vector3d(0, 0, 1),
            2000.0, rng);
        if (!w.escaped) continue;
        ++escaped;
        if (w.exit_position.norm() < 2.0) ++escaped_inside;
    }
    // Before the fix EVERY electron "escaped" from the cavity wall at r = 1.95
    // carrying ~1.9 MeV; now none does.
    BOOST_CHECK_EQUAL(escaped_inside, 0);
    BOOST_CHECK_GT(escaped, 0);   // the walk still reaches the outside at all
}

BOOST_AUTO_TEST_CASE(a_cavity_does_not_let_an_electron_skip_a_shield) {
    // The same defect, in the form that matters: with a real shield outside, an
    // electron that "escaped" at the cavity wall was handed to detector-side
    // transport having skipped the far wall AND the shield entirely.
    // The electron has to actually REACH the cavity for this to test anything:
    // born inside a thin core wall heading inward, so it enters the void at once.
    // (A first version started it in the outer shell, where 2 MeV stops well
    // before the cavity - it passed against the un-fixed code.)
    // Same geometry as the test above (which is known to get an electron as far
    // as the cavity), plus a shield.  Two earlier attempts were vacuous: one
    // started the electron where 2 MeV stops before reaching the cavity, the
    // other in a 0.01 cm sliver that exhausted the walk's step budget.  Both
    // passed against the un-fixed code, which is how they were caught.
    Material soil = make_Soil(), fe = make_Iron();

    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d(0, 0, 0), 2.05, 2.0, Eigen::Matrix3d::Identity());
    sg.set_source_material(&soil);
    sg.add_core(&fe, 0.05);            // fills [1.95, 2.0]; r < 1.95 stays void
    sg.add_shield(&fe, 0.5);           // 0.5 cm of iron outside, at [2.05, 2.55]
    sg.set_source_electron_transport(true);

    std::mt19937_64 rng(4242);
    int escaped_inside = 0;
    const int N = 200;
    for (int i = 0; i < N; ++i) {
        auto w = ElectronCsda::instance().walk_in_source_geometry(
            sg, soil, Eigen::Vector3d(0, 0, -2.02), Eigen::Vector3d(0, 0, 1),
            2000.0, rng);
        if (w.escaped && w.exit_position.norm() < 2.55) ++escaped_inside;
    }
    // Nothing may report escaping from anywhere inside the outer shield surface.
    // Before the fix these "escaped" at the cavity wall, r = 1.95, and were then
    // handed to detector-side transport having skipped the far wall AND the
    // 0.5 cm iron shield entirely.
    BOOST_CHECK_EQUAL(escaped_inside, 0);
}

BOOST_AUTO_TEST_CASE(a_solid_cylinder_is_not_secretly_hollow) {
    // configure_cylindrical() sets cyl_inner_half_length_ to the FULL half-length
    // for a solid cylinder, so testing the cavity triple's maxCoeff called a
    // solid cylinder hollow: it gained a phantom void layer, accepted add_core(),
    // and lost the closed-form and electron-containment fast paths for good.
    Material water = make_Water();
    SourceGeometry solid;
    solid.configure_cylindrical(Eigen::Vector3d(0, 0, 0), 3.0, 3.0,
                                Eigen::Matrix3d::Identity());
    solid.set_source_material(&water);
    BOOST_CHECK_EQUAL(solid.source_layer_index(), 0u);   // no leading void
    BOOST_CHECK(!solid.has_attenuating_interior());
    BOOST_CHECK_EQUAL(solid.layers().size(), 1u);
}

BOOST_AUTO_TEST_CASE(reconfiguring_to_a_shape_without_a_cavity_drops_the_cores) {
    // The gate that routes to the ordered march is cached.  Reconfiguring into a
    // shape that cannot hold cores used to leave it set, and the walker then
    // returned nothing for every trace -- silently dropping the shields.
    Material pb = make_Lead(), soil = make_Soil();
    SourceGeometry sg;
    sg.configure_spherical(Eigen::Vector3d(0, 0, 0), 3.0, 2.0, Eigen::Matrix3d::Identity());
    sg.set_source_material(&soil);
    sg.add_core(&pb, 2.0);
    BOOST_CHECK(sg.has_attenuating_interior());

    sg.configure_point(Eigen::Vector3d(0, 0, 0));
    sg.add_shield(&pb, 0.5);
    BOOST_CHECK(!sg.has_attenuating_interior());
    // The shield must still attenuate: exp(-mu * 0.5), not 1.
    const double t = sg.compute_transmission(Eigen::Vector3d(0, 0, 0),
                                             Eigen::Vector3d(0, 0, 1), 0.662);
    BOOST_CHECK_LT(t, 0.99);
    BOOST_CHECK_CLOSE(t, std::exp(-pb.mu_total(0.662) * 0.5), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
