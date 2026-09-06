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


// The crystal-dimension convention, pinned.
//
// CRYSTAL DIMENSION CONVENTION (geometry/Geometry.h) says transverse extents
// are HALVES and the axial extent is the FULL crystal length. Nothing crashes
// when that is misread -- the MC simply simulates a different detector than the
// one it names, and the result looks plausible. VirtualDepthFit did exactly
// that for a while: make_descriptor doubled dimensions[1] while set_detector
// read it as already-full, so one caller simulated a half-length crystal and
// the other advertised a double-length one.
//
// So these tests tie the three restatements of the convention together: the
// named structs, what Geometry actually traces, and what a DetectorDescriptor
// advertises.

#define BOOST_TEST_MODULE DimensionConventionTests
#include <boost/test/unit_test.hpp>

#include "geometry/Geometry.h"
#include "io/VirtualDepthFit.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <vector>

using namespace ceelo;

namespace {

// BOOST_CHECK_CLOSE takes a PERCENTAGE, not an epsilon; 1e-9% is ~1e-11
// relative, which every comparison here meets because the values are copied
// rather than computed.
constexpr double kTolPercent = 1.0e-9;

// A 3"x3" NaI, and a deliberately UNEQUAL box: three equal extents would let a
// permuted {half_x, half_y, full_length} pass every assertion below.
const CylinderDims kNaI3x3{3.81, 7.62};
const BoxDims kBox{0.5, 0.75, 1.2};

/// Axial extent of the crystal solid an on-axis ray passes through, measured by
/// tracing rather than by reading the members back -- this is the number the
/// physics actually sees.
double traced_crystal_length(const Geometry& g) {
    const std::vector<PathSegment> segs =
        g.trace_ray(Eigen::Vector3d(0.0, 0.0, -50.0), Eigen::Vector3d(0.0, 0.0, 1.0));
    double len = 0.0;
    for (const PathSegment& s : segs) {
        if (s.is_scoring)
            len += s.length();
    }
    return len;
}

} // namespace

// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(StructVectorRoundTrip)

BOOST_AUTO_TEST_CASE(cylinder_round_trips_through_the_serialized_vector) {
    const std::vector<double> v = to_dimensions_vector(kNaI3x3);
    BOOST_REQUIRE_EQUAL(v.size(), 2u);
    BOOST_CHECK_CLOSE(v[0], 3.81, kTolPercent);
    BOOST_CHECK_CLOSE(v[1], 7.62, kTolPercent);  // FULL length, not 3.81

    const CylinderDims back = cylinder_dims_from_vector(v);
    BOOST_CHECK_CLOSE(back.radius_cm, kNaI3x3.radius_cm, kTolPercent);
    BOOST_CHECK_CLOSE(back.full_length_cm, kNaI3x3.full_length_cm, kTolPercent);
}

BOOST_AUTO_TEST_CASE(box_round_trips_and_keeps_the_mixed_layout) {
    const std::vector<double> v = to_dimensions_vector(kBox);
    BOOST_REQUIRE_EQUAL(v.size(), 3u);
    // Two halves then a full length, IN THAT ORDER: the mixed layout is the
    // serialized one, and the three values differ so a permutation is caught.
    BOOST_CHECK_CLOSE(v[0], 0.5, kTolPercent);
    BOOST_CHECK_CLOSE(v[1], 0.75, kTolPercent);
    BOOST_CHECK_CLOSE(v[2], 1.2, kTolPercent);

    const BoxDims back = box_dims_from_vector(v);
    BOOST_CHECK_CLOSE(back.half_x_cm, kBox.half_x_cm, kTolPercent);
    BOOST_CHECK_CLOSE(back.half_y_cm, kBox.half_y_cm, kTolPercent);
    BOOST_CHECK_CLOSE(back.full_length_cm, kBox.full_length_cm, kTolPercent);
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(GeometryHonoursTheNames)

BOOST_AUTO_TEST_CASE(cylinder_traces_the_full_length_it_was_given) {
    const Material nai = make_NaI();
    Geometry g;
    g.set_detector(&nai, kNaI3x3);

    BOOST_CHECK_CLOSE(g.detector_radius(), kNaI3x3.radius_cm, kTolPercent);
    BOOST_CHECK_CLOSE(g.detector_length(), kNaI3x3.full_length_cm, kTolPercent);
    // The one that matters: a 3"x3" is 7.62 cm deep, not 3.81.
    BOOST_CHECK_CLOSE(traced_crystal_length(g), kNaI3x3.full_length_cm, 1e-6);
}

BOOST_AUTO_TEST_CASE(box_traces_the_full_length_and_keeps_halves_transverse) {
    const Material czt = make_CZT();
    Geometry g;
    g.set_detector(&czt, kBox);

    BOOST_CHECK_CLOSE(g.detector_half_x(), kBox.half_x_cm, kTolPercent);
    BOOST_CHECK_CLOSE(g.detector_half_y(), kBox.half_y_cm, kTolPercent);
    BOOST_CHECK_CLOSE(g.detector_length(), kBox.full_length_cm, kTolPercent);
    BOOST_CHECK_CLOSE(traced_crystal_length(g), kBox.full_length_cm, 1e-6);
}

BOOST_AUTO_TEST_CASE(the_vector_entry_point_agrees_with_the_typed_one) {
    const Material nai = make_NaI();
    Geometry typed, from_vector;
    typed.set_detector(&nai, kNaI3x3);
    from_vector.set_detector_from_dimensions_vector(DetectorShape::Cylinder, &nai,
                                                    to_dimensions_vector(kNaI3x3));
    BOOST_CHECK_CLOSE(typed.detector_radius(), from_vector.detector_radius(), kTolPercent);
    BOOST_CHECK_CLOSE(typed.detector_length(), from_vector.detector_length(), kTolPercent);

    const Material czt = make_CZT();
    Geometry btyped, bvector;
    btyped.set_detector(&czt, kBox);
    bvector.set_detector_from_dimensions_vector(DetectorShape::Box, &czt,
                                                to_dimensions_vector(kBox));
    BOOST_CHECK_CLOSE(btyped.detector_half_x(), bvector.detector_half_x(), kTolPercent);
    BOOST_CHECK_CLOSE(btyped.detector_half_y(), bvector.detector_half_y(), kTolPercent);
    BOOST_CHECK_CLOSE(btyped.detector_length(), bvector.detector_length(), kTolPercent);
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(DescriptorsAdvertiseWhatTheyTrace)

// This is the check that would have caught the VirtualDepthFit bug from either
// side: the length the descriptor advertises must be the length that gets traced.
BOOST_AUTO_TEST_CASE(vpd_descriptor_length_is_the_traced_length) {
    const Material nai = make_NaI();
    const DetectorDescriptor cyl = make_descriptor("NaI 3x3", nai, kNaI3x3);
    BOOST_CHECK_CLOSE(cyl.crystal_length_cm, kNaI3x3.full_length_cm, kTolPercent);
    BOOST_CHECK_CLOSE(cyl.crystal_radius_cm, kNaI3x3.radius_cm, kTolPercent);
    BOOST_CHECK_CLOSE(cyl.crystal_diameter_cm, 2.0 * kNaI3x3.radius_cm, kTolPercent);
    {
        Geometry g;
        g.set_detector(&nai, cyl.cylinder_dims());
        BOOST_CHECK_CLOSE(g.detector_length(), cyl.crystal_length_cm, kTolPercent);
        BOOST_CHECK_CLOSE(traced_crystal_length(g), cyl.crystal_length_cm, 1e-6);
    }

    const Material czt = make_CZT();
    const DetectorDescriptor box = make_descriptor("CZT slab", czt, kBox);
    BOOST_CHECK_CLOSE(box.crystal_length_cm, kBox.full_length_cm, kTolPercent);
    BOOST_CHECK_CLOSE(box.half_x_cm, kBox.half_x_cm, kTolPercent);
    BOOST_CHECK_CLOSE(box.half_y_cm, kBox.half_y_cm, kTolPercent);
    {
        Geometry g;
        g.set_detector(&czt, box.box_dims());
        BOOST_CHECK_CLOSE(g.detector_length(), box.crystal_length_cm, kTolPercent);
        BOOST_CHECK_CLOSE(traced_crystal_length(g), box.crystal_length_cm, 1e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END()
