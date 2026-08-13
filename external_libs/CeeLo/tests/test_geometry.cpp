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

#define BOOST_TEST_MODULE GeometryTests
#include <boost/test/unit_test.hpp>

#include "geometry/Geometry.h"
#include "geometry/Cylinder.h"
#include "geometry/Box.h"
#include "geometry/RayTrace.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <optional>

using namespace ceelo;

constexpr double TOL = 1.0e-8;

// ============================================================
//  Ray-Cylinder Intersection Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(CylinderTests)

BOOST_AUTO_TEST_CASE(on_axis_ray_through_cylinder) {
    // Ray along z-axis through a cylinder centered on z-axis.
    // Cylinder: R=2.54 cm, z=[0, 5.08] (a 2" x 2" cylinder)
    Eigen::Vector3d origin(0.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_cylinder(origin, dir, 2.54, 0.0, 5.08);
    BOOST_REQUIRE(hit.has_value());
    BOOST_CHECK_CLOSE(hit->t_enter, 10.0, 1e-6);   // enters at z=0
    BOOST_CHECK_CLOSE(hit->t_exit, 15.08, 1e-6);    // exits at z=5.08
    BOOST_CHECK_CLOSE(hit->length(), 5.08, 1e-6);   // path = full length
}

BOOST_AUTO_TEST_CASE(off_axis_ray_through_cylinder) {
    // Ray parallel to z-axis but offset in x by 1.0 cm.
    // Cylinder: R=2.54, z=[0, 5.08]
    Eigen::Vector3d origin(1.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_cylinder(origin, dir, 2.54, 0.0, 5.08);
    BOOST_REQUIRE(hit.has_value());
    BOOST_CHECK_CLOSE(hit->t_enter, 10.0, 1e-6);
    BOOST_CHECK_CLOSE(hit->t_exit, 15.08, 1e-6);
    // Path length = 5.08 (same as on-axis for parallel ray)
    BOOST_CHECK_CLOSE(hit->length(), 5.08, 1e-6);
}

BOOST_AUTO_TEST_CASE(ray_misses_cylinder) {
    // Ray parallel to z-axis but outside the cylinder radius
    Eigen::Vector3d origin(3.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_cylinder(origin, dir, 2.54, 0.0, 5.08);
    BOOST_CHECK(!hit.has_value());
}

BOOST_AUTO_TEST_CASE(ray_misses_endcap) {
    // Ray hits the infinite cylinder but passes above the finite one
    Eigen::Vector3d origin(1.0, 0.0, 10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_cylinder(origin, dir, 2.54, 0.0, 5.08);
    BOOST_CHECK(!hit.has_value());
}

BOOST_AUTO_TEST_CASE(diagonal_ray_through_cylinder) {
    // Ray entering through the front face at an angle
    // Cylinder: R=5, z=[0, 10]
    // Ray from (-10, 0, -10) toward (1, 0, 1)/sqrt(2) — 45 degrees
    Eigen::Vector3d origin(-10.0, 0.0, -10.0);
    Eigen::Vector3d dir = Eigen::Vector3d(1.0, 0.0, 1.0).normalized();

    auto hit = intersect_cylinder(origin, dir, 5.0, 0.0, 10.0);
    BOOST_REQUIRE(hit.has_value());
    BOOST_CHECK(hit->valid());

    // The ray enters through the front face at z=0.
    // At z=0: t = 10/dir_z = 10/0.7071 = 14.142, x = -10 + 14.142*0.7071 = 0.0
    // That's inside radius 5, so it enters through the front face.
    double t_front = 10.0 / dir.z();
    BOOST_CHECK_CLOSE(hit->t_enter, t_front, 1e-4);

    // Should exit through the side (x² + y² = 25)
    // or through the back face at z=10.
    // At z=10: t = 20/dir_z = 28.28, x = -10 + 28.28*0.7071 = 10.0 > R=5
    // So it exits through the cylindrical surface.
    BOOST_CHECK(hit->t_exit > hit->t_enter);
}

BOOST_AUTO_TEST_CASE(perpendicular_ray_through_cylinder) {
    // Ray going in x-direction, perpendicular to cylinder axis
    // Cylinder: R=2.54, z=[0, 5.08]
    Eigen::Vector3d origin(-10.0, 0.0, 2.54);
    Eigen::Vector3d dir(1.0, 0.0, 0.0);

    auto hit = intersect_cylinder(origin, dir, 2.54, 0.0, 5.08);
    BOOST_REQUIRE(hit.has_value());

    // Enters at x = -2.54, exits at x = +2.54
    // t_enter = -2.54 - (-10) = 7.46, t_exit = 2.54 - (-10) = 12.54
    BOOST_CHECK_CLOSE(hit->t_enter, 7.46, 1e-4);
    BOOST_CHECK_CLOSE(hit->t_exit, 12.54, 1e-4);
    // Chord length through diameter = 2*R = 5.08
    BOOST_CHECK_CLOSE(hit->length(), 5.08, 1e-4);
}

BOOST_AUTO_TEST_CASE(ray_origin_inside_cylinder) {
    // Ray starting inside the cylinder
    Eigen::Vector3d origin(0.0, 0.0, 2.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_cylinder(origin, dir, 2.54, 0.0, 5.08);
    BOOST_REQUIRE(hit.has_value());
    // t_enter should be negative (behind the origin)
    BOOST_CHECK(hit->t_enter < 0.0);
    BOOST_CHECK_CLOSE(hit->t_exit, 3.08, 1e-6); // exits at z=5.08
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Ray-Box Intersection Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(BoxTests)

BOOST_AUTO_TEST_CASE(on_axis_ray_through_box) {
    // Ray along z-axis through a box centered on the z-axis
    // Box: [-2, 2] x [-3, 3] x [0, 5]
    Eigen::Vector3d origin(0.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_box(origin, dir, 2.0, 3.0, 0.0, 5.0);
    BOOST_REQUIRE(hit.has_value());
    BOOST_CHECK_CLOSE(hit->t_enter, 10.0, 1e-6);
    BOOST_CHECK_CLOSE(hit->t_exit, 15.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(diagonal_ray_through_box) {
    // Ray at 45 degrees through a cube: [-1,1]^2 x [0,2]
    Eigen::Vector3d origin(-5.0, 0.0, -5.0);
    Eigen::Vector3d dir = Eigen::Vector3d(1.0, 0.0, 1.0).normalized();

    auto hit = intersect_box(origin, dir, 1.0, 1.0, 0.0, 2.0);
    BOOST_REQUIRE(hit.has_value());
    BOOST_CHECK(hit->valid());
    BOOST_CHECK(hit->length() > 0.0);
}

BOOST_AUTO_TEST_CASE(ray_misses_box) {
    Eigen::Vector3d origin(5.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_box(origin, dir, 2.0, 3.0, 0.0, 5.0);
    BOOST_CHECK(!hit.has_value());
}

BOOST_AUTO_TEST_CASE(perpendicular_ray_through_box) {
    // Ray in x-direction through center of box
    // Box: [-2,2] x [-3,3] x [0,5], ray at y=0, z=2.5
    Eigen::Vector3d origin(-10.0, 0.0, 2.5);
    Eigen::Vector3d dir(1.0, 0.0, 0.0);

    auto hit = intersect_box(origin, dir, 2.0, 3.0, 0.0, 5.0);
    BOOST_REQUIRE(hit.has_value());
    // Enters at x=-2, exits at x=+2
    BOOST_CHECK_CLOSE(hit->t_enter, 8.0, 1e-6);
    BOOST_CHECK_CLOSE(hit->t_exit, 12.0, 1e-6);
    BOOST_CHECK_CLOSE(hit->length(), 4.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(ray_origin_inside_box) {
    Eigen::Vector3d origin(0.0, 0.0, 2.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_box(origin, dir, 2.0, 3.0, 0.0, 5.0);
    BOOST_REQUIRE(hit.has_value());
    BOOST_CHECK(hit->t_enter < 0.0);
    BOOST_CHECK_CLOSE(hit->t_exit, 3.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Bore Hole Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(BoreHoleTests)

BOOST_AUTO_TEST_CASE(on_axis_ray_through_bored_cylinder) {
    // Coaxial HPGe: R_outer=3.0, L=6.0, bore R=0.5, depth=4.0
    // Bore extends from z=2.0 to z=6.0 (from back face)
    // Ray along z-axis should pass through:
    //   [0, 2] — active crystal (in front of bore)
    //   [2, 6] — bore (no crystal)
    // So only 1 segment: [0, 2]
    Eigen::Vector3d origin(0.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    RayHit segs[2];
    int n = intersect_bored_cylinder(origin, dir,
                                      3.0, 0.5,  // outer R, bore R
                                      0.0, 6.0,  // z_min, z_max
                                      2.0, 6.0,  // bore z_start, z_end
                                      segs);

    // On-axis ray (x=0, y=0) is inside bore radius (0.5),
    // so it hits the bore from z=2 to z=6.
    // Active crystal only from z=0 to z=2.
    BOOST_REQUIRE_EQUAL(n, 1);
    BOOST_CHECK_CLOSE(segs[0].t_enter, 10.0, 1e-4);
    BOOST_CHECK_CLOSE(segs[0].t_exit, 12.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(off_axis_ray_misses_bore) {
    // Ray parallel to z-axis at x=1.5, well outside bore radius of 0.5
    // Should get a single full-length segment
    Eigen::Vector3d origin(1.5, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    RayHit segs[2];
    int n = intersect_bored_cylinder(origin, dir,
                                      3.0, 0.5,
                                      0.0, 6.0,
                                      2.0, 6.0,
                                      segs);

    BOOST_REQUIRE_EQUAL(n, 1);
    BOOST_CHECK_CLOSE(segs[0].t_enter, 10.0, 1e-4);
    BOOST_CHECK_CLOSE(segs[0].t_exit, 16.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(diagonal_ray_two_segments) {
    // Ray that enters the crystal, passes through the bore, exits the crystal.
    // Large bore: R_outer=5, bore R=2, bore depth=8 (bore from z=2 to z=10)
    // Crystal: z=[0, 10]
    // Ray in x-z plane at angle, starting outside
    Eigen::Vector3d origin(-8.0, 0.0, 5.0);  // At z=5, which is in the bore region
    Eigen::Vector3d dir(1.0, 0.0, 0.0);  // Perpendicular to z-axis

    RayHit segs[2];
    int n = intersect_bored_cylinder(origin, dir,
                                      5.0, 2.0,
                                      0.0, 10.0,
                                      2.0, 10.0,
                                      segs);

    // Ray enters outer cylinder at x=-5, bore at x=-2, exits bore at x=+2, exits outer at x=+5
    // Segment 1: x in [-5, -2] (active crystal)
    // Segment 2: x in [+2, +5] (active crystal)
    BOOST_REQUIRE_EQUAL(n, 2);

    // First segment: t from 3.0 to 6.0 (x=-5+3=-2 to x=-5+6=-2... wait)
    // origin_x = -8, enters outer at x = -5: t = (-5 - (-8))/1 = 3
    // enters bore at x = -2: t = (-2 - (-8))/1 = 6
    // exits bore at x = +2: t = (2 - (-8))/1 = 10
    // exits outer at x = +5: t = (5 - (-8))/1 = 13
    BOOST_CHECK_CLOSE(segs[0].t_enter, 3.0, 1e-4);
    BOOST_CHECK_CLOSE(segs[0].t_exit, 6.0, 1e-4);
    BOOST_CHECK_CLOSE(segs[1].t_enter, 10.0, 1e-4);
    BOOST_CHECK_CLOSE(segs[1].t_exit, 13.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(ray_in_front_of_bore) {
    // Ray perpendicular, at z=1.0 which is in front of the bore start
    // Bore from z=2 to z=6
    Eigen::Vector3d origin(-10.0, 0.0, 1.0);
    Eigen::Vector3d dir(1.0, 0.0, 0.0);

    RayHit segs[2];
    int n = intersect_bored_cylinder(origin, dir,
                                      3.0, 0.5,
                                      0.0, 6.0,
                                      2.0, 6.0,
                                      segs);

    // At z=1.0, which is before bore start (z=2), so no bore subtraction
    // Full chord through cylinder at z=1.0 and x,y passing through center
    BOOST_REQUIRE_EQUAL(n, 1);
    // Chord length at y=0 through cylinder of R=3: 2*3 = 6.0
    BOOST_CHECK_CLOSE(segs[0].length(), 6.0, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Full Geometry Trace Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(FullGeometryTests)

BOOST_AUTO_TEST_CASE(bare_cylinder_trace) {
    // Simple cylindrical detector with no attenuators or dead layer
    // 2" x 2" NaI: R=2.54 cm, L=5.08 cm
    auto nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {2.54, 5.08});

    Eigen::Vector3d origin(0.0, 0.0, -20.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto segments = geom.trace_ray(origin, dir);
    BOOST_REQUIRE_EQUAL(segments.size(), 1u);
    BOOST_CHECK(segments[0].is_scoring);
    BOOST_CHECK_CLOSE(segments[0].length(), 5.08, 1e-4);
    BOOST_CHECK_EQUAL(segments[0].material, &nai);
}

BOOST_AUTO_TEST_CASE(cylinder_with_dead_layer) {
    // HPGe with 0.5 mm front dead layer, 0.5 mm side dead layer
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &ge, {3.0, 6.0});
    geom.set_dead_layer(0.05, 0.05, 0.0); // 0.5 mm = 0.05 cm

    Eigen::Vector3d origin(0.0, 0.0, -20.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto segments = geom.trace_ray(origin, dir);

    // Should have: dead layer (front) + active crystal + (no back dead layer)
    // Dead layer: z=[0, 0.05], active: z=[0.05, 6.0]
    // Actually with dead layer, the segment structure is:
    //   1. Dead layer front: z in [0, 0.05]  — non-scoring
    //   2. Active crystal: z in [0.05, 6.0]  — scoring

    // Find scoring and non-scoring segments
    double scoring_length = 0.0;
    double non_scoring_length = 0.0;
    for (const auto& seg : segments) {
        if (seg.is_scoring)
            scoring_length += seg.length();
        else
            non_scoring_length += seg.length();
    }

    // Active length should be 6.0 - 0.05 = 5.95 cm
    BOOST_CHECK_CLOSE(scoring_length, 5.95, 1e-2);
    // Dead layer contributes 0.05 cm on front
    BOOST_CHECK_CLOSE(non_scoring_length, 0.05, 1e-2);
}

BOOST_AUTO_TEST_CASE(cylinder_with_attenuator) {
    // NaI with 1 mm Al housing
    auto nai = make_NaI();
    auto al = make_Aluminum();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {2.54, 5.08});
    geom.add_attenuator(&al, 0.1, 0.1, -0.1, 5.18); // 1mm Al all around

    Eigen::Vector3d origin(0.0, 0.0, -20.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto segments = geom.trace_ray(origin, dir);

    // Should have: Al (front) + NaI (scoring) + Al (back)
    double al_length = 0.0;
    double nai_length = 0.0;
    for (const auto& seg : segments) {
        if (seg.material == &al) al_length += seg.length();
        if (seg.material == &nai) nai_length += seg.length();
    }

    BOOST_CHECK_CLOSE(nai_length, 5.08, 1e-2);
    BOOST_CHECK_CLOSE(al_length, 0.2, 2.0); // 1mm front + 1mm back = 2mm
}

BOOST_AUTO_TEST_CASE(bare_box_trace) {
    // Simple box detector: 2x3x5 cm
    auto czt = make_CZT();
    Geometry geom;
    geom.set_detector(DetectorShape::Box, &czt, {1.0, 1.5, 5.0});

    Eigen::Vector3d origin(0.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto segments = geom.trace_ray(origin, dir);
    BOOST_REQUIRE_EQUAL(segments.size(), 1u);
    BOOST_CHECK(segments[0].is_scoring);
    BOOST_CHECK_CLOSE(segments[0].length(), 5.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(ray_misses_geometry) {
    auto nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {2.54, 5.08});

    // Ray that completely misses
    Eigen::Vector3d origin(10.0, 10.0, -20.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto segments = geom.trace_ray(origin, dir);
    BOOST_CHECK(segments.empty());
}

BOOST_AUTO_TEST_CASE(outer_bounding_radius_grows_with_attenuators) {
    auto nai = make_NaI();
    auto pb = make_Lead();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {2.54, 5.08});

    double r0 = geom.outer_bounding_radius();
    BOOST_CHECK_CLOSE(r0, 2.54, 1e-4);

    geom.add_attenuator(&pb, 0.5, 0.5, -0.5, 5.58);
    double r1 = geom.outer_bounding_radius();
    BOOST_CHECK_CLOSE(r1, 3.04, 1e-4); // 2.54 + 0.5

    // Dead layer also adds
    geom.set_dead_layer(0.05, 0.05);
    double r2 = geom.outer_bounding_radius();
    BOOST_CHECK_CLOSE(r2, 3.09, 1e-4); // 2.54 + 0.05 + 0.5
}

BOOST_AUTO_TEST_CASE(bore_hole_geometry_trace) {
    // Coaxial HPGe: R=3.0 cm, L=6.0 cm, bore R=0.5 cm, bore depth=4.0 cm
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &ge, {3.0, 6.0});
    geom.set_bore_hole(0.5, 4.0);

    // On-axis ray: should only traverse the non-bore part
    Eigen::Vector3d origin(0.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto segments = geom.trace_ray(origin, dir);

    double scoring_length = 0.0;
    for (const auto& seg : segments) {
        if (seg.is_scoring) scoring_length += seg.length();
    }

    // On-axis ray: bore from z=2 to z=6 (depth=4, back face at z=6)
    // Active crystal only from z=0 to z=2
    BOOST_CHECK_CLOSE(scoring_length, 2.0, 1e-2);
}

BOOST_AUTO_TEST_SUITE_END()
