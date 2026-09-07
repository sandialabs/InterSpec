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

#include <cmath>
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
