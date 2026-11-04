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

/// @file test_convex_separation.cpp
/// @brief convex_separation / convex_entry_distance against closed forms (and a brute-force check
///        where there is none), and GeometryDescriptor::Envelope::solid() against the envelope's
///        own exit distances.

#define BOOST_TEST_MODULE ConvexSeparationTests
#include <boost/test/unit_test.hpp>

#include "geometry/ConvexSeparation.h"
#include "io/DetectorResponse.h"
#include "materials/Material.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEps = 1.0e-7;

}  // namespace

BOOST_AUTO_TEST_CASE(separation_closed_forms) {
    // Sphere - sphere: centre distance less both radii; 0 once they overlap.
    BOOST_CHECK_SMALL(convex_separation(ConvexSolid::sphere({0, 0, 0}, 2.0),
                                        ConvexSolid::sphere({10, 0, 0}, 3.0)) - 5.0, kEps);
    BOOST_CHECK_SMALL(convex_separation(ConvexSolid::sphere({0, 0, 0}, 2.0),
                                        ConvexSolid::sphere({4, 0, 0}, 3.0)), kEps);

    // Parallel cylinders: radial and axial gaps combine in quadrature.
    const ConvexSolid a = ConvexSolid::cylinder({0, 0, 0}, Eigen::Vector3d::UnitZ(), 2.0, 5.0);
    const ConvexSolid b = ConvexSolid::cylinder({10, 0, 0}, Eigen::Vector3d::UnitZ(), 3.0, 5.0);
    BOOST_CHECK_SMALL(convex_separation(a, b) - 5.0, kEps);   // beside, z overlapping
    const ConvexSolid c = ConvexSolid::cylinder({10, 0, 16}, Eigen::Vector3d::UnitZ(), 3.0, 5.0);
    BOOST_CHECK_SMALL(convex_separation(a, c) - std::hypot(5.0, 6.0), kEps);   // off the corner

    // Box - cylinder, axis-aligned: rectangle-to-circle in x-y, combined with the z gap.
    const ConvexSolid box = ConvexSolid::box({0, 0, 0}, Eigen::Vector3d(2, 3, 4));
    const ConvexSolid cyl = ConvexSolid::cylinder({8, 0, 0}, Eigen::Vector3d::UnitZ(), 1.0, 2.0);
    BOOST_CHECK_SMALL(convex_separation(box, cyl) - 5.0, kEps);
    BOOST_CHECK_SMALL(convex_separation(ConvexSolid::point({0, 0, 9}), box) - 5.0, kEps);

    // The closest points lie on each solid, the separation apart.
    Eigen::Vector3d pa, pb;
    const double sep = convex_separation(box, cyl, &pa, &pb);
    BOOST_CHECK_SMALL((pa - pb).norm() - sep, kEps);
    BOOST_CHECK_SMALL((box.project(pa) - pa).norm(), kEps);
    BOOST_CHECK_SMALL((cyl.project(pb) - pb).norm(), kEps);
}

BOOST_AUTO_TEST_CASE(separation_perpendicular_cylinders) {
    // The side-on case has no closed form: brute force over one surface.
    const ConvexSolid a = ConvexSolid::cylinder({0, 0, 0}, Eigen::Vector3d::UnitZ(), 2.0, 6.0);
    const ConvexSolid b = ConvexSolid::cylinder({7, 1, 3}, Eigen::Vector3d::UnitX(), 1.5, 2.0);
    double brute = 1.0e30;
    for (int i = 0; i <= 120; ++i) {
        for (int j = 0; j <= 60; ++j) {
            const double phi = 2.0 * kPi * i / 120.0, z = -6.0 + 12.0 * j / 60.0;
            const Eigen::Vector3d p(2.0 * std::cos(phi), 2.0 * std::sin(phi), z);
            brute = std::min(brute, (b.project(p) - p).norm());
        }
    }
    const double sep = convex_separation(a, b);
    BOOST_CHECK_LE(sep, brute + 1.0e-9);
    BOOST_CHECK_GT(sep, brute - 0.02);   // the brute grid's own resolution
}

BOOST_AUTO_TEST_CASE(entry_distance_closed_forms) {
    constexpr double inf = std::numeric_limits<double>::infinity();
    const Eigen::Vector3d x = Eigen::Vector3d::UnitX();

    const ConvexSolid sphere = ConvexSolid::sphere({10, 0, 0}, 3.0);
    BOOST_CHECK_SMALL(convex_entry_distance(sphere, {0, 0, 0}, x) - 7.0, kEps);
    BOOST_CHECK_EQUAL(convex_entry_distance(sphere, {0, 0, 0}, -x), inf);         // pointing away
    BOOST_CHECK_EQUAL(convex_entry_distance(sphere, {10, 1, 0}, x), 0.0);         // starts inside
    BOOST_CHECK_EQUAL(convex_entry_distance(sphere, {0, 4, 0}, x), inf);          // passes by

    // A cylinder along z, entered through its side and through its end cap.
    const ConvexSolid cyl = ConvexSolid::cylinder({0, 0, 0}, Eigen::Vector3d::UnitZ(), 2.0, 5.0);
    BOOST_CHECK_SMALL(convex_entry_distance(cyl, {-10, 0, 1}, x) - 8.0, kEps);
    BOOST_CHECK_SMALL(convex_entry_distance(cyl, {0, 1, -9}, Eigen::Vector3d::UnitZ()) - 4.0, kEps);
    BOOST_CHECK_EQUAL(convex_entry_distance(cyl, {-10, 0, 6}, x), inf);           // above the cap

    // A box rotated 90 degrees about z: its x half-width now lies along y.
    const Eigen::Matrix3d rz = Eigen::AngleAxisd(0.5 * kPi, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const ConvexSolid box = ConvexSolid::box({0, 0, 0}, Eigen::Vector3d(1, 4, 1), rz);
    BOOST_CHECK_SMALL(convex_entry_distance(box, {-10, 0, 0}, x) - 6.0, kEps);
    BOOST_CHECK_SMALL(convex_entry_distance(box, {0, -10, 0}, Eigen::Vector3d::UnitY()) - 9.0, kEps);

    BOOST_CHECK_EQUAL(convex_entry_distance(ConvexSolid::point({5, 0, 0}), {0, 0, 0}, x), inf);
}

BOOST_AUTO_TEST_CASE(transformed_moves_rigidly) {
    // Rotating and shifting both solids together leaves their separation unchanged.
    const Eigen::Matrix3d R =
        Eigen::AngleAxisd(0.7, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
    const Eigen::Vector3d t(4.0, -2.0, 7.5);
    const ConvexSolid box = ConvexSolid::box({1, 0, 0}, Eigen::Vector3d(2, 3, 4));
    const ConvexSolid cyl = ConvexSolid::cylinder({9, 1, 2}, Eigen::Vector3d(0, 1, 1), 1.0, 2.0);
    BOOST_CHECK_SMALL(convex_separation(box.transformed(R, t), cyl.transformed(R, t))
                          - convex_separation(box, cyl), 1.0e-6);

    // A point maps through the transform as the solid does.
    const ConvexSolid moved = cyl.transformed(R, t);
    BOOST_CHECK_SMALL((moved.center - (R * cyl.center + t)).norm(), 1.0e-12);
    BOOST_CHECK_SMALL((moved.axis - R * cyl.axis).norm(), 1.0e-12);
}

BOOST_AUTO_TEST_CASE(envelope_solid_matches_exit_distances) {
    // An NaI crystal in an aluminium can with a back plate: the envelope's solid must agree with the
    // envelope's own exit distances, so a source's clearance reads the same either way.
    GeometryDescriptor gd;
    gd.set_dimensions(CylinderDims{2.54, 5.08});
    gd.materials = {MaterialSpec::from(make_NaI()), MaterialSpec::from(make_Aluminum())};
    gd.crystal_material_index = 0;
    LayerSpec can;
    can.material_index = 1;
    can.front_thickness_cm = 0.1;
    can.side_thickness_cm = 0.2;
    can.z_end_cm = 5.08;
    gd.layers.push_back(can);
    gd.set_layer_back_thickness_cm(0, 0.5);

    const GeometryDescriptor::Envelope env = gd.outer_envelope();
    const ConvexSolid solid = env.solid();
    BOOST_REQUIRE(solid.kind == ConvexSolid::Kind::Cylinder);

    for (const double theta_deg : {0.0, 30.0, 89.0, 91.0, 135.0, 180.0}) {
        const double th = theta_deg * kPi / 180.0;
        // Directions from the crystal-face origin; +z runs into the detector, so "in front" is -z.
        const Eigen::Vector3d u(std::sin(th), 0.0, -std::cos(th));
        const double exit = env.exit_distance_cm(u);
        for (const double clearance : {0.05, 1.0, 10.0}) {
            const Eigen::Vector3d p = (exit + clearance) * u;
            const double sep = convex_separation(ConvexSolid::point(p), solid);
            // The clearance along u is an upper bound: another part of the outside may be nearer.
            BOOST_CHECK_LE(sep, clearance + 1.0e-9);
            BOOST_CHECK_GT(sep, 0.0);
        }
    }

    // Straight in front and straight behind, the nearest outside is along the axis.
    const Eigen::Vector3d front(0.0, 0.0, env.z_min_cm - 3.0), back(0.0, 0.0, env.z_max_cm + 2.0);
    BOOST_CHECK_SMALL(convex_separation(ConvexSolid::point(front), solid) - 3.0, kEps);
    BOOST_CHECK_SMALL(convex_separation(ConvexSolid::point(back), solid) - 2.0, kEps);
    BOOST_CHECK_SMALL(convex_entry_distance(solid, front, Eigen::Vector3d::UnitZ()) - 3.0, kEps);
}
