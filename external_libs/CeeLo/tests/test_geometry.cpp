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
#include <random>
#include <vector>
#include <utility>

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
//  Bulletized (rounded front edge) Cylinder Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(BulletizedCylinderTests)

// ANGLE GEM35-70 crystal: R = 2.915 cm, L = 6.89 cm, bulletizing radius 8 mm.
constexpr double BR   = 2.915;
constexpr double BL   = 6.89;
constexpr double BRB  = 0.8;
constexpr double BRHO = BR - BRB;   // ring radius, 2.115
constexpr double BZC  = BRB;        // ring plane, front face at z = 0

FrontFillet gem_fillet() { return make_front_fillet(BR, 0.0, BRB); }

/// Membership test written straight from the geometric definition, with no
/// ray algebra in it: the oracle the intersection results are checked against.
bool bullet_contains(const Eigen::Vector3d& p) {
    const double rho = std::hypot(p.x(), p.y());
    if (rho > BR || p.z() < 0.0 || p.z() > BL) return false;
    if (rho <= BRHO || p.z() >= BZC) return true;
    const double dr = rho - BRHO, dz = p.z() - BZC;
    return dr * dr + dz * dz <= BRB * BRB;
}

/// Bisect the membership oracle for the boundary crossing between a parameter
/// known to be outside the solid and one known to be inside.
double bisect_surface(const Eigen::Vector3d& o, const Eigen::Vector3d& d,
                      double t_out, double t_in) {
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (t_out + t_in);
        if (mid == t_out || mid == t_in) break;
        if (bullet_contains(o + d * mid)) t_in = mid; else t_out = mid;
    }
    return 0.5 * (t_out + t_in);
}

/// Ray-torus quartic residual -- an algebraically independent check that a
/// corrected endpoint really lies on the fillet surface.  Normalised by the
/// largest term so the tolerance is scale-free.
double torus_quartic_residual(const Eigen::Vector3d& o, const Eigen::Vector3d& d,
                              double t) {
    const Eigen::Vector3d op(o.x(), o.y(), o.z() - BZC);
    const double beta  = op.dot(d);
    const double gamma = op.squaredNorm() + BRHO * BRHO - BRB * BRB;
    const double a = d.x() * d.x() + d.y() * d.y();
    const double b = 2.0 * (o.x() * d.x() + o.y() * d.y());
    const double c = o.x() * o.x() + o.y() * o.y();

    const double u = t * t + 2.0 * beta * t + gamma;
    const double v = 4.0 * BRHO * BRHO * (a * t * t + b * t + c);
    const double scale = std::max(1.0, std::max(std::abs(u * u), std::abs(v)));
    return std::abs(u * u - v) / scale;
}

BOOST_AUTO_TEST_CASE(bullet_axial_ray_inside_ring_radius_is_untouched) {
    // At rho < rho_c the ray meets the flat part of the front face, a surface
    // the sharp and bulletized solids share: the result must be *identical* to
    // the sharp cylinder, not merely close (no correction may be applied).
    Eigen::Vector3d origin(1.0, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto sharp  = intersect_cylinder(origin, dir, BR, 0.0, BL);
    auto bullet = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL,
                                                gem_fillet());
    BOOST_REQUIRE(sharp.has_value());
    BOOST_REQUIRE(bullet.has_value());
    BOOST_CHECK(bullet->t_enter == sharp->t_enter);
    BOOST_CHECK(bullet->t_exit  == sharp->t_exit);
    BOOST_CHECK_CLOSE(bullet->length(), BL, 1e-10);
}

BOOST_AUTO_TEST_CASE(bullet_axial_ray_through_fillet) {
    // rho between rho_c and R: the ray enters on the fillet, at
    // z = z_c - sqrt(r_b^2 - (rho - rho_c)^2).
    const double rho = 2.5;
    Eigen::Vector3d origin(rho, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    auto hit = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL,
                                             gem_fillet());
    BOOST_REQUIRE(hit.has_value());

    const double dr = rho - BRHO;
    const double z_entry = BZC - std::sqrt(BRB * BRB - dr * dr);
    BOOST_CHECK_CLOSE(hit->t_enter, 10.0 + z_entry, 1e-9);
    BOOST_CHECK_CLOSE(hit->t_exit, 10.0 + BL, 1e-12);
    BOOST_CHECK_SMALL(torus_quartic_residual(origin, dir, hit->t_enter), 1e-12);
}

BOOST_AUTO_TEST_CASE(bullet_radial_ray_below_ring_plane) {
    // Perpendicular to the axis at z < z_c: both ends land on the fillet, at
    // rho = rho_c + sqrt(r_b^2 - (z - z_c)^2).
    const double z = 0.4;
    Eigen::Vector3d origin(-10.0, 0.0, z);
    Eigen::Vector3d dir(1.0, 0.0, 0.0);

    auto hit = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL,
                                             gem_fillet());
    BOOST_REQUIRE(hit.has_value());

    const double dz = z - BZC;
    const double rho_s = BRHO + std::sqrt(BRB * BRB - dz * dz);
    BOOST_CHECK_CLOSE(hit->t_enter, 10.0 - rho_s, 1e-9);
    BOOST_CHECK_CLOSE(hit->t_exit, 10.0 + rho_s, 1e-9);
    BOOST_CHECK_SMALL(torus_quartic_residual(origin, dir, hit->t_enter), 1e-12);
    BOOST_CHECK_SMALL(torus_quartic_residual(origin, dir, hit->t_exit), 1e-12);
}

BOOST_AUTO_TEST_CASE(bullet_radial_ray_above_ring_plane_is_untouched) {
    // Above z_c the side wall is the shared surface: identical to sharp.
    Eigen::Vector3d origin(-10.0, 0.0, 1.5);
    Eigen::Vector3d dir(1.0, 0.0, 0.0);

    auto sharp  = intersect_cylinder(origin, dir, BR, 0.0, BL);
    auto bullet = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL,
                                                gem_fillet());
    BOOST_REQUIRE(sharp.has_value());
    BOOST_REQUIRE(bullet.has_value());
    BOOST_CHECK(bullet->t_enter == sharp->t_enter);
    BOOST_CHECK(bullet->t_exit  == sharp->t_exit);
}

BOOST_AUTO_TEST_CASE(bullet_tangent_ray_misses_and_offsets_flip_it) {
    // 45-degree ray in the x-z plane whose line is exactly r_b from the fillet
    // ring centre: the whole sharp chord lies in the corner zone, so this
    // exercises the convex-minimum branch.  Nudging the line inward makes it
    // clip the fillet; nudging it outward makes it miss cleanly.
    const double inv = 1.0 / std::sqrt(2.0);
    Eigen::Vector3d dir(inv, 0.0, inv);
    const double k_tangent = -(BRHO - BZC) - BRB * std::sqrt(2.0);

    auto shot = [&](double k) {
        return intersect_bulletized_cylinder(Eigen::Vector3d(0.0, 0.0, k), dir,
                                             BR, 0.0, BL, gem_fillet());
    };

    BOOST_CHECK(!shot(k_tangent - 1e-6).has_value());   // outside: clean miss
    auto inward = shot(k_tangent + 1e-3);               // inside: clips fillet
    BOOST_REQUIRE(inward.has_value());
    BOOST_CHECK(inward->length() > 0.0);
    BOOST_CHECK(inward->length() < 0.2);                // only a corner clip
}

BOOST_AUTO_TEST_CASE(bullet_grazing_corner_chord_matches_oracle) {
    // Both sharp endpoints sit in the removed corner and the chord never
    // leaves the corner zone.  Check the corrected interval against the
    // membership oracle and against the ray-torus quartic.
    const double inv = 1.0 / std::sqrt(2.0);
    Eigen::Vector3d dir(inv, 0.0, inv);
    const double k = -(BRHO - BZC) - BRB * std::sqrt(2.0) + 0.05;
    Eigen::Vector3d origin(0.0, 0.0, k);

    auto sharp  = intersect_cylinder(origin, dir, BR, 0.0, BL);
    auto bullet = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL,
                                                gem_fillet());
    BOOST_REQUIRE(sharp.has_value());
    BOOST_REQUIRE(bullet.has_value());

    // Corrected interval must be a sub-interval of the sharp one.
    BOOST_CHECK(bullet->t_enter >= sharp->t_enter);
    BOOST_CHECK(bullet->t_exit  <= sharp->t_exit);

    const double mid = 0.5 * (bullet->t_enter + bullet->t_exit);
    BOOST_REQUIRE(bullet_contains(origin + dir * mid));
    BOOST_CHECK_SMALL(
        bullet->t_enter - bisect_surface(origin, dir, sharp->t_enter, mid), 1e-11);
    BOOST_CHECK_SMALL(
        bullet->t_exit - bisect_surface(origin, dir, sharp->t_exit, mid), 1e-11);
    BOOST_CHECK_SMALL(torus_quartic_residual(origin, dir, bullet->t_enter), 1e-10);
    BOOST_CHECK_SMALL(torus_quartic_residual(origin, dir, bullet->t_exit), 1e-10);
}

BOOST_AUTO_TEST_CASE(bullet_zero_radius_is_bitwise_sharp) {
    // A zero bulletization radius must reproduce intersect_cylinder exactly,
    // ray for ray -- the regression gate for the common sharp-crystal case.
    const FrontFillet none = make_front_fillet(BR, 0.0, 0.0);
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    int hits = 0;
    for (int i = 0; i < 10000; ++i) {
        Eigen::Vector3d origin(6.0 * u(rng), 6.0 * u(rng), 5.0 * u(rng) - 3.0);
        Eigen::Vector3d dir(u(rng), u(rng), u(rng));
        if (dir.norm() < 1e-6) continue;
        dir.normalize();

        auto sharp  = intersect_cylinder(origin, dir, BR, 0.0, BL);
        auto bullet = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL, none);
        BOOST_REQUIRE_EQUAL(sharp.has_value(), bullet.has_value());
        if (sharp) {
            BOOST_REQUIRE(sharp->t_enter == bullet->t_enter);
            BOOST_REQUIRE(sharp->t_exit  == bullet->t_exit);
            ++hits;
        }
    }
    BOOST_CHECK_GT(hits, 100);  // the sweep actually exercised the solid
}

BOOST_AUTO_TEST_CASE(bullet_random_rays_match_oracle) {
    // Fuzz the general case against the membership oracle, including rays that
    // start inside the crystal (the electron walk re-traces from interior
    // vertices, so negative t_enter must be handled).
    const FrontFillet f = gem_fillet();
    std::mt19937_64 rng(987654321);
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    int checked = 0;
    for (int i = 0; i < 20000 && checked < 4000; ++i) {
        Eigen::Vector3d origin(3.2 * u(rng), 3.2 * u(rng), 3.6 * u(rng) + 3.4);
        Eigen::Vector3d dir(u(rng), u(rng), u(rng));
        if (dir.norm() < 1e-6) continue;
        dir.normalize();

        auto sharp  = intersect_cylinder(origin, dir, BR, 0.0, BL);
        auto bullet = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL, f);
        if (!sharp) { BOOST_REQUIRE(!bullet.has_value()); continue; }

        if (!bullet) {
            // A spurious miss is the dangerous failure mode -- the ray really
            // does cross the crystal and we would silently drop it -- so it
            // has to be asserted, not skipped. Sample the sharp chord densely
            // and require the solid to be empty along all of it.
            const double t0 = std::max(sharp->t_enter, 0.0);
            for (int k = 0; k <= 2000; ++k) {
                const double t = t0 + (sharp->t_exit - t0) * k / 2000.0;
                BOOST_REQUIRE_MESSAGE(!bullet_contains(origin + dir * t),
                    "returned no hit but t=" << t << " is inside the solid");
            }
            continue;
        }

        // Endpoints on the surface, interior actually inside.
        const double mid = 0.5 * (bullet->t_enter + bullet->t_exit);
        BOOST_REQUIRE(bullet_contains(origin + dir * mid));
        BOOST_REQUIRE(bullet->t_enter >= sharp->t_enter - 1e-12);
        BOOST_REQUIRE(bullet->t_exit  <= sharp->t_exit + 1e-12);
        BOOST_REQUIRE_SMALL(
            bullet->t_enter - bisect_surface(origin, dir, sharp->t_enter, mid), 1e-9);
        BOOST_REQUIRE_SMALL(
            bullet->t_exit - bisect_surface(origin, dir, sharp->t_exit, mid), 1e-9);
        ++checked;
    }
    BOOST_CHECK_GT(checked, 1000);
}

BOOST_AUTO_TEST_CASE(bullet_interior_origin_splits_full_chord) {
    // Firing forward and backward from an interior point must reconstruct the
    // full line chord.
    const FrontFillet f = gem_fillet();
    std::mt19937_64 rng(555);
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    int checked = 0;
    for (int i = 0; i < 4000 && checked < 400; ++i) {
        Eigen::Vector3d origin(2.6 * u(rng), 2.6 * u(rng), 0.5 * (u(rng) + 1.0) * BL);
        if (!bullet_contains(origin)) continue;
        Eigen::Vector3d dir(u(rng), u(rng), u(rng));
        if (dir.norm() < 1e-6) continue;
        dir.normalize();

        auto fwd = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL, f);
        auto bwd = intersect_bulletized_cylinder(origin, -dir, BR, 0.0, BL, f);
        BOOST_REQUIRE(fwd.has_value());
        BOOST_REQUIRE(bwd.has_value());
        // Origin inside => the interval straddles t = 0 in both directions.
        BOOST_REQUIRE_LT(fwd->t_enter, 1e-9);
        BOOST_REQUIRE_GT(fwd->t_exit, -1e-9);
        // Forward and backward reaches must add up to the full chord.
        BOOST_REQUIRE_SMALL(fwd->t_exit + bwd->t_exit
                            - (fwd->t_exit - fwd->t_enter), 1e-9);
        ++checked;
    }
    BOOST_CHECK_GT(checked, 100);
}

BOOST_AUTO_TEST_CASE(bullet_reverse_ray_hits_same_points) {
    // The corrected surface points must not depend on which way the ray runs.
    const FrontFillet f = gem_fillet();
    Eigen::Vector3d origin(2.4, 0.6, -4.0);
    Eigen::Vector3d dir(0.05, -0.02, 1.0);
    dir.normalize();

    auto fwd = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL, f);
    BOOST_REQUIRE(fwd.has_value());

    // Fire back along the same line from beyond the exit point.
    const Eigen::Vector3d far = origin + dir * (fwd->t_exit + 5.0);
    auto rev = intersect_bulletized_cylinder(far, -dir, BR, 0.0, BL, f);
    BOOST_REQUIRE(rev.has_value());

    const Eigen::Vector3d fwd_in  = origin + dir * fwd->t_enter;
    const Eigen::Vector3d rev_out = far - dir * rev->t_exit;
    BOOST_CHECK_SMALL((fwd_in - rev_out).norm(), 1e-10);
}

BOOST_AUTO_TEST_CASE(bullet_volume_matches_closed_form) {
    // Integrate the chord length the intersector returns over the front face
    // and compare with V = pi R^2 L - [2 pi rho_c r_b^2 (1 - pi/4) + pi r_b^3/3].
    const FrontFillet f = gem_fillet();
    const Eigen::Vector3d dir(0.0, 0.0, 1.0);

    const int N = 40000;
    double integral = 0.0;   // integral of chord(rho) * 2 pi rho drho
    for (int i = 0; i < N; ++i) {
        const double rho = BR * (i + 0.5) / N;
        auto hit = intersect_bulletized_cylinder(Eigen::Vector3d(rho, 0.0, -5.0),
                                                 dir, BR, 0.0, BL, f);
        const double chord = hit ? hit->length() : 0.0;
        integral += chord * 2.0 * M_PI * rho * (BR / N);
    }

    const double v_removed = 2.0 * M_PI * BRHO * BRB * BRB * (1.0 - M_PI / 4.0)
                           + M_PI * BRB * BRB * BRB / 3.0;
    BOOST_CHECK_CLOSE(v_removed, 2.361283, 0.01);  // sanity on the formula

    // Judge the fillet against the fillet, not against the whole crystal:
    // the corner is only 1.3 % of the volume, so normalising to the total
    // would let a fillet that is wrong by more than a percent pass.
    const double sharp_v = M_PI * BR * BR * BL;
    BOOST_CHECK_CLOSE(sharp_v - integral, v_removed, 0.5);   // 0.5 % of 2.36 cm^3
    BOOST_CHECK_CLOSE(integral, sharp_v - v_removed, 0.02);
}

BOOST_AUTO_TEST_CASE(bullet_axis_parallel_edge_cases) {
    const FrontFillet f = gem_fillet();
    const Eigen::Vector3d up(0.0, 0.0, 1.0);

    // Just inside the ring radius: untouched front face.
    auto a = intersect_bulletized_cylinder(Eigen::Vector3d(BRHO - 1e-9, 0.0, -1.0),
                                           up, BR, 0.0, BL, f);
    BOOST_REQUIRE(a.has_value());
    BOOST_CHECK_CLOSE(a->length(), BL, 1e-6);

    // Just outside it: the fillet is tangent to the front face at rho_c, so a
    // ray a nanometre beyond it still sees essentially the full length. This
    // and the case above straddle the predicate boundary in in_removed_corner.
    auto b = intersect_bulletized_cylinder(Eigen::Vector3d(BRHO + 1e-9, 0.0, -1.0),
                                           up, BR, 0.0, BL, f);
    BOOST_REQUIRE(b.has_value());
    BOOST_CHECK_CLOSE(b->length(), BL, 1e-6);

    // Hard against the side wall: the fillet removes almost the whole r_b
    // run-up, short by the sqrt term that vanishes only exactly at rho = R.
    const double rho_wall = BR - 1e-9;
    const double dr_wall = rho_wall - BRHO;
    auto c = intersect_bulletized_cylinder(Eigen::Vector3d(rho_wall, 0.0, -1.0),
                                           up, BR, 0.0, BL, f);
    BOOST_REQUIRE(c.has_value());
    BOOST_CHECK_CLOSE(c->length(),
                      BL - BRB + std::sqrt(BRB * BRB - dr_wall * dr_wall), 1e-8);

    // Origin already inside the z-slab and inside the solid.
    auto d = intersect_bulletized_cylinder(Eigen::Vector3d(1.0, 0.0, 3.0),
                                           up, BR, 0.0, BL, f);
    BOOST_REQUIRE(d.has_value());
    BOOST_CHECK_CLOSE(d->t_exit, BL - 3.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(rounded_bore_tip_lengthens_crystal_off_axis) {
    // GEM35-70 bore: r = 0.495, depth 5.54 from the back face.  A round-tipped
    // drill keeps the stated depth (apex at bore_z_start) but leaves material
    // in the corner ring, so off-axis rays meet the bore later by
    // r - sqrt(r^2 - rho^2).
    const double r_bore = 0.495, depth = 5.54;
    const double bore_z_start = BL - depth;
    const FrontFillet none = make_front_fillet(BR, 0.0, 0.0);
    const double rho = 0.3;
    Eigen::Vector3d origin(rho, 0.0, -10.0);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    RayHit flat[2], round[2];
    int n_flat = intersect_shaped_bored_cylinder(origin, dir, BR, r_bore, 0.0, BL,
                                                 bore_z_start, BL, none, false, flat);
    int n_round = intersect_shaped_bored_cylinder(origin, dir, BR, r_bore, 0.0, BL,
                                                  bore_z_start, BL, none, true, round);
    BOOST_REQUIRE_EQUAL(n_flat, 1);
    BOOST_REQUIRE_EQUAL(n_round, 1);

    const double extra = r_bore - std::sqrt(r_bore * r_bore - rho * rho);
    BOOST_CHECK_CLOSE(flat[0].length(), bore_z_start, 1e-9);
    BOOST_CHECK_CLOSE(round[0].length(), bore_z_start + extra, 1e-9);

    // On the axis the hemisphere apex sits at the flat bottom: no change.
    Eigen::Vector3d on_axis(0.0, 0.0, -10.0);
    intersect_shaped_bored_cylinder(on_axis, dir, BR, r_bore, 0.0, BL,
                                    bore_z_start, BL, none, false, flat);
    intersect_shaped_bored_cylinder(on_axis, dir, BR, r_bore, 0.0, BL,
                                    bore_z_start, BL, none, true, round);
    BOOST_CHECK_CLOSE(flat[0].length(), round[0].length(), 1e-9);
}

BOOST_AUTO_TEST_CASE(rounded_bore_tip_volume_matches_closed_form) {
    // Rounding the tip returns pi r^3 / 3 of germanium to the crystal.
    const double r_bore = 0.495, depth = 5.54;
    const double bore_z_start = BL - depth;
    const FrontFillet none = make_front_fillet(BR, 0.0, 0.0);
    const Eigen::Vector3d dir(0.0, 0.0, 1.0);

    const int N = 20000;
    double flat_v = 0.0, round_v = 0.0;
    for (int i = 0; i < N; ++i) {
        const double rho = BR * (i + 0.5) / N;
        const double w = 2.0 * M_PI * rho * (BR / N);
        RayHit segs[2];
        Eigen::Vector3d o(rho, 0.0, -5.0);

        int nf = intersect_shaped_bored_cylinder(o, dir, BR, r_bore, 0.0, BL,
                                                 bore_z_start, BL, none, false, segs);
        for (int s = 0; s < nf; ++s) flat_v += segs[s].length() * w;

        int nr = intersect_shaped_bored_cylinder(o, dir, BR, r_bore, 0.0, BL,
                                                 bore_z_start, BL, none, true, segs);
        for (int s = 0; s < nr; ++s) round_v += segs[s].length() * w;
    }

    BOOST_CHECK_CLOSE(flat_v, M_PI * BR * BR * BL - M_PI * r_bore * r_bore * depth,
                      0.02);
    BOOST_CHECK_CLOSE(round_v - flat_v, M_PI * r_bore * r_bore * r_bore / 3.0, 0.5);
}

BOOST_AUTO_TEST_CASE(bulletized_geometry_composes_with_bore_and_dead_layer) {
    // Full GEM35-70: bulletized front edge, coaxial bore with a rounded tip,
    // and a 0.7 mm dead layer.  The active volume's fillet is the crystal
    // fillet offset inward, which for equal front/side thicknesses shares the
    // same ring centre with radius r_b - t.
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{BR, BL});
    geom.set_bullet_radius(BRB);
    geom.set_bore_hole(0.495, 5.54, /*rounded_tip=*/true);
    geom.set_dead_layer(0.07, 0.07);

    const double t_dl = 0.07;
    const double rb_act = BRB - t_dl;

    // Axial ray inside the ring radius: front face and dead layer are flat
    // there, so the dead layer is exactly its own thickness.
    {
        std::vector<PathSegment> segs;
        geom.trace_ray(Eigen::Vector3d(1.0, 0.0, -10.0),
                       Eigen::Vector3d(0.0, 0.0, 1.0), segs);
        BOOST_REQUIRE(!segs.empty());
        BOOST_CHECK_CLOSE(segs.front().t_start, 10.0, 1e-9);
        BOOST_CHECK(!segs.front().is_scoring);
        BOOST_CHECK_CLOSE(segs.front().length(), t_dl, 1e-6);
    }

    // Axial ray through the fillet: entry into crystal and into active volume
    // are both on arcs about the shared centre (rho_c, z_c).
    {
        const double rho = 2.5;
        const double dr = rho - BRHO;
        const double z_outer = BZC - std::sqrt(BRB * BRB - dr * dr);
        const double z_active = BZC - std::sqrt(rb_act * rb_act - dr * dr);

        std::vector<PathSegment> segs;
        geom.trace_ray(Eigen::Vector3d(rho, 0.0, -10.0),
                       Eigen::Vector3d(0.0, 0.0, 1.0), segs);
        BOOST_REQUIRE_GE(segs.size(), 2u);
        BOOST_CHECK_CLOSE(segs[0].t_start, 10.0 + z_outer, 1e-8);
        BOOST_CHECK(!segs[0].is_scoring);
        BOOST_CHECK_CLOSE(segs[1].t_start, 10.0 + z_active, 1e-8);
        BOOST_CHECK(segs[1].is_scoring);
    }

    // Fan of rays: segments must stay sorted and non-overlapping, and the
    // dead-layer plus active length must equal the bulletized crystal chord.
    const FrontFillet outer_f = gem_fillet();
    for (int i = 0; i < 200; ++i) {
        const double ang = -0.9 + 1.8 * i / 199.0;
        Eigen::Vector3d origin(0.0, 0.0, -12.0);
        Eigen::Vector3d dir(std::sin(ang), 0.0, std::cos(ang));
        dir.normalize();

        std::vector<PathSegment> segs;
        geom.trace_ray(origin, dir, segs);
        if (segs.empty()) continue;

        double total = 0.0;
        for (size_t s = 0; s < segs.size(); ++s) {
            BOOST_REQUIRE_LE(segs[s].t_start, segs[s].t_end + 1e-12);
            if (s > 0) BOOST_REQUIRE_GE(segs[s].t_start, segs[s - 1].t_end - 1e-9);
            total += segs[s].length();
        }

        auto crystal = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL, outer_f);
        BOOST_REQUIRE(crystal.has_value());
        // The bore is the only material removed from the crystal chord.
        BOOST_CHECK_LE(total, crystal->length() + 1e-9);
    }
}

BOOST_AUTO_TEST_CASE(bulletized_geometry_survives_unequal_dead_layers) {
    // Unequal front/side dead layers use the documented fallback convention for
    // the active fillet radius.  The exact arc is a matter of convention, but
    // the trace must stay well formed: ordered, non-overlapping segments whose
    // active part never exceeds the crystal chord.  Includes the degenerate
    // case of a dead layer thicker than the fillet (active edge goes sharp).
    auto ge = make_HPGe();
    const FrontFillet outer_f = gem_fillet();

    for (auto dl : {std::make_pair(0.20, 0.05), std::make_pair(0.05, 0.20),
                    std::make_pair(1.00, 1.00)}) {
        Geometry geom;
        geom.set_detector(&ge, CylinderDims{BR, BL});
        geom.set_bullet_radius(BRB);
        geom.set_dead_layer(dl.first, dl.second);

        std::mt19937_64 rng(13579);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        for (int i = 0; i < 2000; ++i) {
            Eigen::Vector3d origin(8.0 * u(rng), 8.0 * u(rng), -6.0);
            Eigen::Vector3d dir(0.35 * u(rng), 0.35 * u(rng), 1.0);
            dir.normalize();

            std::vector<PathSegment> segs;
            geom.trace_ray(origin, dir, segs);
            for (size_t s = 0; s < segs.size(); ++s) {
                BOOST_REQUIRE_GE(segs[s].length(), -1e-12);
                if (s > 0) BOOST_REQUIRE_GE(segs[s].t_start, segs[s - 1].t_end - 1e-9);
            }
            auto crystal = intersect_bulletized_cylinder(origin, dir, BR, 0.0, BL,
                                                         outer_f);
            const double active = compute_active_path_length(segs);
            BOOST_REQUIRE_LE(active, (crystal ? crystal->length() : 0.0) + 1e-9);
        }
    }
}

BOOST_AUTO_TEST_CASE(bulletized_geometry_removes_only_corner_material) {
    // Against an otherwise identical sharp crystal, bulletization may only
    // shorten chords, and must shorten them strictly for corner-crossing rays.
    auto ge = make_HPGe();
    Geometry sharp, bullet;
    sharp.set_detector(&ge, CylinderDims{BR, BL});
    bullet.set_detector(&ge, CylinderDims{BR, BL});
    bullet.set_bullet_radius(BRB);

    std::mt19937_64 rng(24680);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    int shortened = 0, compared = 0;

    for (int i = 0; i < 5000; ++i) {
        Eigen::Vector3d origin(8.0 * u(rng), 8.0 * u(rng), -6.0);
        Eigen::Vector3d dir(0.3 * u(rng), 0.3 * u(rng), 1.0);
        dir.normalize();

        std::vector<PathSegment> a, b;
        sharp.trace_ray(origin, dir, a);
        bullet.trace_ray(origin, dir, b);
        const double la = compute_active_path_length(a);
        const double lb = compute_active_path_length(b);
        if (la <= 0.0) { BOOST_REQUIRE_SMALL(lb, 1e-12); continue; }

        BOOST_REQUIRE_LE(lb, la + 1e-12);
        ++compared;
        if (lb < la - 1e-9) ++shortened;
    }
    BOOST_CHECK_GT(compared, 500);
    BOOST_CHECK_GT(shortened, 20);   // corner rays really were clipped
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
    geom.set_detector(&nai, CylinderDims{2.54, 5.08});

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
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
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
    geom.set_detector(&nai, CylinderDims{2.54, 5.08});
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
    geom.set_detector(&czt, BoxDims{1.0, 1.5, 5.0});

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
    geom.set_detector(&nai, CylinderDims{2.54, 5.08});

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
    geom.set_detector(&nai, CylinderDims{2.54, 5.08});

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
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
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


// ============================================================
//  The bore is empty over its WHOLE length, including where it
//  passes through the dead layer.
//
//  Geometry::trace_cylinder_geometry used to build the dead layer as
//  (outer crystal - UNBORED active cylinder), which refills with germanium
//  every part of the bore that lies outside the active volume.  The bore
//  escapes the active volume in four ways:
//      z > dl_z_max    a back dead layer to drill through
//      z < dl_z_min    a bore deep enough to reach the front dead layer
//      rho > dl_r      a bore wider than the active radius
//      the removed front corner of a bulletized crystal
//  The tests below pin each one, plus the invariant that ties the tracer to
//  the polycone Geant4Export writes and to the volumes
//  src/DetectorGeometryDiagram.cpp reports.
//
//  These go through Geometry::trace_ray on purpose.  The older
//  rounded_bore_tip_volume_matches_closed_form calls
//  intersect_shaped_bored_cylinder directly, so it never reaches the
//  dead-layer block and structurally cannot see any of this.
// ============================================================

namespace {

struct Interval { double lo, hi; };

std::vector<Interval> segments_of(const std::vector<PathSegment>& segs, bool scoring) {
    std::vector<Interval> out;
    for (const auto& s : segs) {
        if (s.is_scoring == scoring) out.push_back({s.t_start, s.t_end});
    }
    return out;
}

double summed_length(const std::vector<PathSegment>& segs) {
    double s = 0.0;
    for (const auto& g : segs) s += g.length();
    return s;
}

/// Chord through the solid the crystal actually IS: the outer cylinder (filleted
/// or sharp) minus the bore.  Independent of the dead layer, and the same
/// description Geant4Export's polycone uses, so the tracer's active + dead
/// segments must reproduce it exactly.
double solid_chord(const Eigen::Vector3d& o, const Eigen::Vector3d& d,
                   double R, double L, double r_bore, double depth,
                   double bullet_radius, bool rounded) {
    const FrontFillet f = (bullet_radius > 0.0)
        ? make_front_fillet(R, 0.0, bullet_radius)
        : FrontFillet{0.0, 0.0, 0.0, 0.0};
    RayHit segs[2];
    const int n = intersect_shaped_bored_cylinder(o, d, R, r_bore, 0.0, L,
                                                  L - depth, L, f, rounded, segs);
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        if (segs[i].valid()) s += segs[i].length();
    }
    return s;
}

/// A seeded fan of rays: origins on a shell around the crystal plus some
/// interior ones, directions uniform on the sphere.  Deterministic.
std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>
ray_fan(unsigned seed, int n, double R, double L) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0), phi(0.0, 2.0 * M_PI);
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> out;
    out.reserve(n);
    for (int i = 0; i < n; ++i) {
        // Half the origins outside the crystal, half inside it.
        const double scale = (i % 2 == 0) ? 3.0 : 0.9;
        const double cz = u(rng), sz = std::sqrt(std::max(0.0, 1.0 - cz * cz));
        const double p = phi(rng);
        Eigen::Vector3d o(scale * R * sz * std::cos(p), scale * R * sz * std::sin(p),
                          0.5 * L + scale * L * cz);
        const double dcz = u(rng), dsz = std::sqrt(std::max(0.0, 1.0 - dcz * dcz));
        const double dp = phi(rng);
        Eigen::Vector3d d(dsz * std::cos(dp), dsz * std::sin(dp), dcz);
        out.emplace_back(o, d.normalized());
    }
    return out;
}

} // namespace


/// T7 -- the no-op proof.  Whenever the bore lies entirely inside the active
/// volume (no back dead layer, apex behind the front one, inside the active
/// radius, clear of the active fillet) the new bore subtraction is provably a
/// no-op, and "provably" here means BIT FOR BIT, not to some tolerance.  That
/// is what keeps the innermost transport loop of every ordinary coaxial HPGe --
/// and every stored golden response, all of which use back = 0 -- unchanged.
///
/// The pre-fix arithmetic is re-implemented below as the reference, so this
/// test is meaningful run against either version of the tracer.  It was written
/// first and confirmed green on the UNMODIFIED tracer.
BOOST_AUTO_TEST_CASE(bore_inside_active_volume_leaves_dead_layer_bit_identical) {
    auto ge = make_HPGe();
    const double R = 3.0, L = 6.0, r_bore = 0.5, depth = 4.0;
    const double t_f = 0.07, t_s = 0.07;

    for (double bullet : {0.0, 0.8}) {
        for (bool rounded : {false, true}) {
            Geometry geom;
            geom.set_detector(&ge, CylinderDims{R, L});
            if (bullet > 0.0) geom.set_bullet_radius(bullet);
            geom.set_bore_hole(r_bore, depth, rounded);
            geom.set_dead_layer(t_f, t_s, 0.0);   // back = 0: the bore cannot escape

            const double dl_r = R - t_s, dl_z0 = t_f, dl_z1 = L;
            const bool bulletized = (bullet > 0.0);
            const FrontFillet outer_f = bulletized
                ? make_front_fillet(R, 0.0, bullet) : FrontFillet{0.0, 0.0, 0.0, 0.0};
            const FrontFillet active_f = bulletized
                ? make_front_fillet(dl_r, dl_z0, std::max(0.0, bullet - std::max(t_f, t_s)))
                : FrontFillet{0.0, 0.0, 0.0, 0.0};

            for (const auto& ray : ray_fan(20260918u, 2000, R, L)) {
                const Eigen::Vector3d& o = ray.first;
                const Eigen::Vector3d& d = ray.second;

                // --- the pre-fix dead-layer arithmetic, verbatim ---
                std::vector<Interval> ref;
                auto outer_hit = bulletized
                    ? intersect_bulletized_cylinder(o, d, R, 0.0, L, outer_f)
                    : intersect_cylinder(o, d, R, 0.0, L);
                auto inner_hit = bulletized
                    ? intersect_bulletized_cylinder(o, d, dl_r, dl_z0, dl_z1, active_f)
                    : intersect_cylinder(o, d, dl_r, dl_z0, dl_z1);
                if (outer_hit && outer_hit->valid()) {
                    const double t0_out = std::max(outer_hit->t_enter, 0.0);
                    const double t1_out = outer_hit->t_exit;
                    if (inner_hit && inner_hit->valid()) {
                        const double t0_in = std::max(inner_hit->t_enter, 0.0);
                        const double t1_in = inner_hit->t_exit;
                        if (t0_out < t0_in - 1e-12) ref.push_back({t0_out, t0_in});
                        if (t1_in < t1_out - 1e-12) ref.push_back({t1_in, t1_out});
                    } else {
                        ref.push_back({t0_out, t1_out});
                    }
                }

                const std::vector<Interval> got =
                    segments_of(geom.trace_ray(o, d), /*scoring=*/false);

                BOOST_REQUIRE_EQUAL(got.size(), ref.size());
                for (size_t k = 0; k < ref.size(); ++k) {
                    // Deliberate == on doubles: the claim is bit identity.
                    BOOST_CHECK_EQUAL(got[k].lo, ref[k].lo);
                    BOOST_CHECK_EQUAL(got[k].hi, ref[k].hi);
                }
            }
        }
    }
}


/// T1 -- the plain case: a back dead layer the bore is drilled through.
BOOST_AUTO_TEST_CASE(back_dead_layer_does_not_plug_the_bore) {
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
    geom.set_bore_hole(0.5, 4.0);            // bore z in [2, 6]
    geom.set_dead_layer(0.05, 0.05, 0.5);    // active z in [0.05, 5.5]

    const auto segs = geom.trace_ray(Eigen::Vector3d(0.0, 0.0, -10.0),
                                     Eigen::Vector3d(0.0, 0.0, 1.0));

    // Front dead layer, then active crystal up to the bore apex.  Nothing
    // beyond z = 2 (t = 12): the rest of the chord is the hole.  Before the fix
    // there was a third segment, a spurious [15.5, 16] of dead germanium.
    BOOST_REQUIRE_EQUAL(segs.size(), 2u);
    BOOST_CHECK(!segs[0].is_scoring);
    BOOST_CHECK_CLOSE(segs[0].t_start, 10.0, 1e-9);
    BOOST_CHECK_CLOSE(segs[0].t_end, 10.05, 1e-9);
    BOOST_CHECK(segs[1].is_scoring);
    BOOST_CHECK_CLOSE(segs[1].t_start, 10.05, 1e-9);
    BOOST_CHECK_CLOSE(segs[1].t_end, 12.0, 1e-9);
    for (const auto& s : segs) BOOST_CHECK_LE(s.t_start, 12.0 + 1e-9);
}


/// T2 -- the control for T1: a ray that misses the bore must still see the
/// whole back dead layer.  Fails if someone "fixes" T1 by deleting it.
BOOST_AUTO_TEST_CASE(back_dead_layer_survives_clear_of_the_bore) {
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
    geom.set_bore_hole(0.5, 4.0);
    geom.set_dead_layer(0.05, 0.05, 0.5);

    const auto segs = geom.trace_ray(Eigen::Vector3d(0.8, 0.0, -10.0),
                                     Eigen::Vector3d(0.0, 0.0, 1.0));

    BOOST_REQUIRE_EQUAL(segs.size(), 3u);
    BOOST_CHECK(!segs[0].is_scoring);
    BOOST_CHECK(segs[1].is_scoring);
    BOOST_CHECK(!segs[2].is_scoring);
    BOOST_CHECK_CLOSE(segs[2].t_start, 15.5, 1e-9);
    BOOST_CHECK_CLOSE(segs[2].t_end, 16.0, 1e-9);
}


/// T3 -- a bore whose apex sits AHEAD of the active slab, so on the axis there
/// is no active crystal at all and only the front dead layer is material.
BOOST_AUTO_TEST_CASE(bore_reaching_past_the_front_dead_layer_is_empty) {
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
    geom.set_bore_hole(0.5, 5.9);            // apex at z = 0.1
    geom.set_dead_layer(0.2, 0.05, 0.3);     // active z in [0.2, 5.7]

    const auto segs = geom.trace_ray(Eigen::Vector3d(0.0, 0.0, -10.0),
                                     Eigen::Vector3d(0.0, 0.0, 1.0));

    // Only z in [0, 0.1] is germanium.  Pre-fix: two dead segments, 0.5 cm.
    BOOST_REQUIRE_EQUAL(segs.size(), 1u);
    BOOST_CHECK(!segs[0].is_scoring);
    BOOST_CHECK_CLOSE(segs[0].t_start, 10.0, 1e-9);
    BOOST_CHECK_CLOSE(segs[0].t_end, 10.1, 1e-9);
}


/// T4 -- the inner_hit == nullopt branch, which used to be the worst case: a ray
/// crossing the front dead layer ahead of a deep bore misses the active cylinder
/// entirely, and the WHOLE crystal chord -- hole included -- was scored as
/// germanium.
BOOST_AUTO_TEST_CASE(ray_missing_the_active_cylinder_still_sees_the_bore) {
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
    geom.set_bore_hole(0.5, 5.9);
    geom.set_dead_layer(0.2, 0.05, 0.3);

    // z = 0.15 is inside the front dead layer, ahead of the active slab.
    const auto segs = geom.trace_ray(Eigen::Vector3d(-10.0, 0.0, 0.15),
                                     Eigen::Vector3d(1.0, 0.0, 0.0));

    BOOST_REQUIRE_EQUAL(segs.size(), 2u);
    BOOST_CHECK(!segs[0].is_scoring);
    BOOST_CHECK(!segs[1].is_scoring);
    BOOST_CHECK_CLOSE(segs[0].t_start, 7.0, 1e-9);
    BOOST_CHECK_CLOSE(segs[0].t_end, 9.5, 1e-9);
    BOOST_CHECK_CLOSE(segs[1].t_start, 10.5, 1e-9);
    BOOST_CHECK_CLOSE(segs[1].t_end, 13.0, 1e-9);
}


/// T5 -- T4 with a round-tipped drill, so the hole's half-width at the ray's
/// depth comes off the tip hemisphere rather than the straight tube.
BOOST_AUTO_TEST_CASE(rounded_bore_tip_is_empty_in_the_front_dead_layer) {
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
    geom.set_bore_hole(0.5, 5.9, /*rounded_tip=*/true);  // apex 0.1, centre 0.6
    geom.set_dead_layer(0.2, 0.05, 0.3);

    const auto segs = geom.trace_ray(Eigen::Vector3d(-10.0, 0.0, 0.15),
                                     Eigen::Vector3d(1.0, 0.0, 0.0));

    // Ball centred (0,0,0.6) of radius 0.5; at z = 0.15 the half-width is
    // sqrt(0.5^2 - 0.45^2) = sqrt(0.0475).
    const double h = std::sqrt(0.0475);
    BOOST_REQUIRE_EQUAL(segs.size(), 2u);
    BOOST_CHECK_SMALL(segs[0].t_start - 7.0, 1e-12);
    BOOST_CHECK_SMALL(segs[0].t_end - (10.0 - h), 1e-12);
    BOOST_CHECK_SMALL(segs[1].t_start - (10.0 + h), 1e-12);
    BOOST_CHECK_SMALL(segs[1].t_end - 13.0, 1e-12);
}


/// T6 -- the radial escape: a bore wider than the ACTIVE radius, because the
/// side dead layer ate the crystal down past it.  GeometryDescriptor::problems()
/// rejects this on the load path (GeometryProblem::BoreInsideDeadLayer), but
/// Geometry's own API accepts it -- set_bore_hole only checks the OUTER radius
/// and set_dead_layer checks nothing -- so the tracer has to be right anyway.
BOOST_AUTO_TEST_CASE(bore_wider_than_the_active_radius_is_still_empty) {
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{3.0, 6.0});
    geom.set_dead_layer(0.0, 2.6, 0.0);      // active radius 0.4
    geom.set_bore_hole(0.5, 4.0);            // wider than that; bore z in [2, 6]

    // rho = 0.45 lies between the active radius and the bore radius.
    const auto segs = geom.trace_ray(Eigen::Vector3d(0.45, 0.0, -10.0),
                                     Eigen::Vector3d(0.0, 0.0, 1.0));

    BOOST_REQUIRE_EQUAL(segs.size(), 1u);
    BOOST_CHECK(!segs[0].is_scoring);
    BOOST_CHECK_CLOSE(segs[0].t_start, 10.0, 1e-9);
    BOOST_CHECK_CLOSE(segs[0].t_end, 12.0, 1e-9);   // pre-fix: 16.0
}


/// T11 -- a BLUNT round-tipped bore (depth < 2*bore_radius, which set_bore_hole
/// permits) must not eat material past its own mouth.  The capping ball would
/// otherwise reach beyond bore_z_end, where the union with the straight tube is
/// NOT convex -- the profile steps down from bore_radius to the ball's radius --
/// so the hull of the two intervals spans non-bore space.
///
/// Callers clamp the bore to the solid they cut, which hides this for every
/// geometry whose solid ends at the crystal back face.  A NEGATIVE back dead
/// layer is the exception: the active cylinder then extends past that face, the
/// gap falls inside it, and real active material was dropped (1.8% of the
/// overhang band before intersect_bore() started confining the ball).
BOOST_AUTO_TEST_CASE(blunt_rounded_bore_does_not_eat_past_its_mouth) {
    // The bore solid, stated directly rather than as an interval.
    const double r_bore = 0.5, L = 6.0, depth = 0.6;   // blunt: 0.6 < 2*0.5
    const double z_apex = L - depth, z_tip = z_apex + r_bore;   // ball top 6.4 > L
    // Tolerance is on the SURFACE, not along the ray: a grazing chord turns a
    // 1e-14 cm surface error into millimetres of arc length, so testing the
    // ray-parameter span would flag rounding as a defect.
    const double kSurf = 1e-9;
    auto in_bore = [&](const Eigen::Vector3d& p) {
        const double rho = std::hypot(p.x(), p.y());
        if (p.z() > L + kSurf || p.z() < z_apex - kSurf) return false;
        if (p.z() >= z_tip) return rho <= r_bore + kSurf;
        const double dz = p.z() - z_tip;
        return std::hypot(rho, dz) <= r_bore + kSurf;
    };

    // (a) the interval never covers anything that is not bore
    std::mt19937 rng(31337u);
    std::uniform_real_distribution<double> ut(0.0, 1.0);
    int checked = 0;
    for (const auto& ray : ray_fan(5150u, 14000, 3.0, L)) {
        const std::optional<RayHit> h = intersect_bore(ray.first, ray.second, r_bore,
                                                       z_apex, L, /*rounded_tip=*/true);
        if (!h || !h->valid()) continue;
        const double a = std::max(h->t_enter, 0.0), b = h->t_exit;
        for (int k = 0; k <= 60; ++k) {
            const Eigen::Vector3d p = ray.first + (a + (b - a) * k / 60.0) * ray.second;
            BOOST_REQUIRE_MESSAGE(in_bore(p),
                "intersect_bore covers a point outside the bore at z=" << p.z());
            ++checked;
        }
    }
    BOOST_CHECK_GT(checked, 1000);   // the fan must actually be hitting the bore

    // (b) and nothing downstream drops active material past the back face
    auto ge = make_HPGe();
    Geometry geom;
    geom.set_detector(&ge, CylinderDims{3.0, L});
    geom.set_bore_hole(r_bore, depth, /*rounded_tip=*/true);
    geom.set_dead_layer(0.07, 0.07, -0.4);   // NEGATIVE: active runs to z = 6.4
    const double dl_r = 3.0 - 0.07, dl_z1 = L + 0.4;

    int overhang = 0, dropped = 0;
    for (const auto& ray : ray_fan(6161u, 14000, 3.0, L)) {
        const auto segs = geom.trace_ray(ray.first, ray.second);
        for (int k = 0; k < 40; ++k) {
            const double t = 30.0 * ut(rng);
            const Eigen::Vector3d p = ray.first + t * ray.second;
            // the band above the crystal back face where the active cylinder still reaches
            if (p.z() <= L + 1e-6 || p.z() >= dl_z1 - 1e-6) continue;
            if (std::hypot(p.x(), p.y()) > dl_r - 1e-6 || in_bore(p)) continue;
            ++overhang;
            bool scoring = false;
            for (const auto& sg : segs) {
                if (t > sg.t_start && t < sg.t_end) { scoring = sg.is_scoring; break; }
            }
            if (!scoring) ++dropped;
        }
    }
    BOOST_CHECK_GT(overhang, 500);          // the probe has to reach the band
    BOOST_CHECK_EQUAL(dropped, 0);
}

/// T8 -- the invariant that ties all three descriptions of the crystal together:
/// active + dead is the chord of (outer solid - bore), whatever the ray.  That
/// is exactly what Geant4Export writes as a single polycone and what
/// src/DetectorGeometryDiagram.cpp reports as v_solid - v_bore.
BOOST_AUTO_TEST_CASE(active_plus_dead_is_the_crystal_minus_the_bore) {
    auto ge = make_HPGe();
    const double R = 3.0, L = 6.0, r_bore = 0.5, depth = 4.5;

    for (double bullet : {0.0, 0.8}) {
        for (bool rounded : {false, true}) {
            Geometry geom;
            geom.set_detector(&ge, CylinderDims{R, L});
            if (bullet > 0.0) geom.set_bullet_radius(bullet);
            geom.set_bore_hole(r_bore, depth, rounded);
            geom.set_dead_layer(0.07, 0.07, 0.3);   // the combination that used to break

            for (const auto& ray : ray_fan(4242u, 1500, R, L)) {
                const auto segs = geom.trace_ray(ray.first, ray.second);

                BOOST_CHECK_SMALL(summed_length(segs)
                                      - solid_chord(ray.first, ray.second, R, L,
                                                    r_bore, depth, bullet, rounded),
                                  1e-9);

                // Structure: sorted, non-overlapping, positive, and never more
                // than the three dead pieces the geometry can produce.
                int n_dead = 0;
                for (size_t i = 0; i < segs.size(); ++i) {
                    BOOST_CHECK_GE(segs[i].t_end, segs[i].t_start);
                    if (!segs[i].is_scoring) ++n_dead;
                    if (i) BOOST_CHECK_GE(segs[i].t_start, segs[i - 1].t_end - 1e-12);
                }
                BOOST_CHECK_LE(n_dead, 3);
            }
        }
    }
}


/// T9 -- absolute dead-layer volume against the closed form, by the same
/// parallel-beam quadrature rounded_bore_tip_volume_matches_closed_form uses,
/// but summing the NON-scoring segments of a full trace_ray.
BOOST_AUTO_TEST_CASE(traced_dead_layer_volume_matches_closed_form) {
    auto ge = make_HPGe();
    const double R = 3.0, L = 6.0, r_bore = 0.5, depth = 4.0;
    const double t_f = 0.05, t_s = 0.05, t_b = 0.5;

    Geometry geom;
    geom.set_detector(&ge, CylinderDims{R, L});
    geom.set_bore_hole(r_bore, depth);
    geom.set_dead_layer(t_f, t_s, t_b);

    // Midpoint quadrature error here is set by the radial discontinuities the
    // chord length has -- at the bore edge and at the active radius -- so this
    // needs more rings than rounded_bore_tip_volume_matches_closed_form, which
    // crosses only one.
    const int N = 200000;
    double dead_v = 0.0;
    std::vector<PathSegment> segs;
    for (int i = 0; i < N; ++i) {
        const double rho = R * (i + 0.5) / N;
        const double w = 2.0 * M_PI * rho * (R / N);
        geom.trace_ray(Eigen::Vector3d(rho, 0.0, -20.0), Eigen::Vector3d(0.0, 0.0, 1.0), segs);
        for (const auto& s : segs) {
            if (!s.is_scoring) dead_v += s.length() * w;
        }
    }

    // dead = solid - bore - active, with the bore counted over its whole length.
    const double dl_r = R - t_s, dl_z0 = t_f, dl_z1 = L - t_b;
    const double bore_z0 = L - depth;
    const double v_solid = M_PI * R * R * L;
    const double v_bore = M_PI * r_bore * r_bore * depth;
    const double v_bore_in_active =
        M_PI * r_bore * r_bore * std::max(0.0, std::min(dl_z1, L) - std::max(bore_z0, dl_z0));
    const double v_active = M_PI * dl_r * dl_r * (dl_z1 - dl_z0) - v_bore_in_active;
    const double v_dead = v_solid - v_bore - v_active;

    BOOST_CHECK_CLOSE(dead_v, v_dead, 0.02);

    // And record the size of the error the fix removed: the pre-fix tracer
    // reported the phantom plug pi*r^2*t_b on top, which this check rejects.
    const double v_dead_prefix = v_dead + M_PI * r_bore * r_bore * t_b;
    BOOST_CHECK_GT(std::fabs(dead_v - v_dead_prefix) / v_dead_prefix, 0.01);
}


/// T10 -- an independent oracle.  Classify points along the ray with direct
/// predicates instead of interval algebra, and require the segment list to
/// agree; immune to exactly the class of mistake being fixed.  Run on the
/// configuration where every code path meets: bulletized, round-tipped bore,
/// and a back dead layer.
BOOST_AUTO_TEST_CASE(segments_agree_with_a_direct_point_classifier) {
    auto ge = make_HPGe();
    const double R = 3.0, L = 6.0, bullet = 0.8;
    const double r_bore = 0.5, depth = 4.5;
    const double t_f = 0.07, t_s = 0.07, t_b = 0.3;

    Geometry geom;
    geom.set_detector(&ge, CylinderDims{R, L});
    geom.set_bullet_radius(bullet);
    geom.set_bore_hole(r_bore, depth, /*rounded_tip=*/true);
    geom.set_dead_layer(t_f, t_s, t_b);

    const double dl_r = R - t_s, dl_z0 = t_f, dl_z1 = L - t_b;
    const double bore_z0 = L - depth, z_tip = bore_z0 + r_bore;
    const FrontFillet outer_f = make_front_fillet(R, 0.0, bullet);
    const FrontFillet active_f =
        make_front_fillet(dl_r, dl_z0, std::max(0.0, bullet - std::max(t_f, t_s)));

    // Inside a (possibly filleted) cylinder: within the sharp cylinder and not
    // in the removed corner.
    auto in_cyl = [](const Eigen::Vector3d& p, double r, double z0, double z1,
                     const FrontFillet& f) {
        const double rho = std::hypot(p.x(), p.y());
        if (rho > r || p.z() < z0 || p.z() > z1) return false;
        if (f.r_b <= 0.0 || rho <= f.rho_c || p.z() >= f.z_c) return true;
        const double dr = rho - f.rho_c, dz = p.z() - f.z_c;
        return dr * dr + dz * dz <= f.r_b * f.r_b;
    };
    auto in_bore = [&](const Eigen::Vector3d& p) {
        const double rho = std::hypot(p.x(), p.y());
        if (p.z() > L || p.z() < bore_z0) return false;
        if (p.z() >= z_tip) return rho <= r_bore;
        const double dz = p.z() - z_tip;
        return rho * rho + dz * dz <= r_bore * r_bore;
    };
    // Distance-to-surface proxies, so samples that graze a boundary are skipped.
    auto near_surface = [&](const Eigen::Vector3d& p) {
        const double rho = std::hypot(p.x(), p.y());
        const double eps = 1e-7;
        return std::fabs(rho - R) < eps || std::fabs(rho - dl_r) < eps
            || std::fabs(rho - r_bore) < eps || std::fabs(p.z()) < eps
            || std::fabs(p.z() - L) < eps || std::fabs(p.z() - dl_z0) < eps
            || std::fabs(p.z() - dl_z1) < eps || std::fabs(p.z() - bore_z0) < eps
            || std::fabs(std::hypot(rho - outer_f.rho_c, p.z() - outer_f.z_c)
                         - outer_f.r_b) < eps
            || std::fabs(std::hypot(rho - active_f.rho_c, p.z() - active_f.z_c)
                         - active_f.r_b) < eps
            || std::fabs(std::hypot(rho, p.z() - z_tip) - r_bore) < eps;
    };

    std::mt19937 rng(777u);
    std::uniform_real_distribution<double> ut(0.0, 1.0);

    // Sample t only where the ray is inside the crystal's bounding cylinder.
    // Spraying t over the whole line puts almost every sample in empty space,
    // and the plug being tested for is a 0.24 cm3 disk -- it would be missed.
    int sampled_in_bore_beyond_active = 0;
    for (const auto& ray : ray_fan(909u, 2000, R, L)) {
        const auto segs = geom.trace_ray(ray.first, ray.second);
        auto span = intersect_cylinder(ray.first, ray.second, R, 0.0, L);
        if (!span || !span->valid()) continue;
        const double t_lo = std::max(span->t_enter, 0.0), t_hi = span->t_exit;
        for (int k = 0; k < 80; ++k) {
            const double t = t_lo + (t_hi - t_lo) * ut(rng);
            const Eigen::Vector3d p = ray.first + t * ray.second;
            if (near_surface(p)) continue;

            const bool solid = in_cyl(p, R, 0.0, L, outer_f) && !in_bore(p);
            const bool active = solid && in_cyl(p, dl_r, dl_z0, dl_z1, active_f);

            const PathSegment* found = nullptr;
            for (const auto& s : segs) {
                if (t > s.t_start && t < s.t_end) { found = &s; break; }
            }

            if (in_bore(p) && p.z() > dl_z1) ++sampled_in_bore_beyond_active;

            if (!solid) {
                BOOST_CHECK_MESSAGE(found == nullptr,
                    "point outside the crystal (or inside the bore) landed in a segment");
            } else {
                BOOST_REQUIRE_MESSAGE(found != nullptr, "solid point landed in no segment");
                BOOST_CHECK_EQUAL(found->is_scoring, active);
            }
        }
    }

    // The oracle is only worth anything if it actually visited the region the
    // fix is about -- the part of the bore past the active volume, where the
    // phantom plug used to be.  Guard against a future edit quietly sampling
    // its way around the bug.
    BOOST_CHECK_MESSAGE(sampled_in_bore_beyond_active > 20,
        "only " << sampled_in_bore_beyond_active << " samples landed in the bore"
        " beyond the active volume; the oracle is not exercising the fix");
}

BOOST_AUTO_TEST_SUITE_END()
