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

#define BOOST_TEST_MODULE DetectorEtendueTests
#include <boost/test/unit_test.hpp>

// MC-free checks of the detector-side etendue line set (io/DetectorEtendue.h):
// the sampling measure, the reversed-integral identity against the per-point
// aperture kernel, determinism, and the stored-response probability accessor.

#include "geometry/Geometry.h"
#include "io/DetectorEtendue.h"
#include "io/DetectorResponse.h"
#include "io/EfficiencyTransfer.h"
#include "io/LowDiscrepancy.h"
#include "io/ResponseKernel.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;

const Material& mat_NaI() { static Material m = make_NaI();      return m; }
const Material& mat_Al()  { static Material m = make_Aluminum(); return m; }

// Bare 3"x3" NaI cylinder: clean analytic hull.
Geometry bare_nai_3x3() {
    Geometry g;
    g.set_detector(DetectorShape::Cylinder, &mat_NaI(), {3.81, 7.62});
    return g;
}

// A box crystal with a thin can, to exercise the five-face hull.
Geometry boxed_czt() {
    Geometry g;
    g.set_detector(DetectorShape::Box, &mat_NaI(), {1.0, 0.75, 1.5});
    g.add_attenuator(&mat_Al(), 0.05, 0.05, -0.05, 1.5);
    return g;
}

GeometryDescriptor nai3x3_descriptor() {
    GeometryDescriptor gd;
    gd.shape = DetectorShape::Cylinder;
    gd.dimensions_cm = {3.81, 7.62};
    gd.materials = {MaterialSpec::from(make_NaI()), MaterialSpec::from(make_Aluminum())};
    gd.crystal_material_index = 0;
    LayerSpec can;
    can.material_index = 1;
    can.front_thickness_cm = 0.05;
    can.side_thickness_cm = 0.05;
    can.z_end_cm = 7.62;
    gd.layers.push_back(can);
    return gd;
}

AnchorCurve make_anchor() {
    AnchorCurve a;
    a.energies_keV = {60.0, 100.0, 300.0, 662.0, 1000.0};
    a.eff = {2.0e-2, 1.3e-2, 5.0e-3, 3.0e-3, 2.2e-3};
    a.frac_sigma = {0.003, 0.003, 0.003, 0.003, 0.003};
    return a;
}

// Chord length of the (infinite) line through `o` along unit `d` inside the
// sphere (c, R); 0 when it misses.
double sphere_chord(const Eigen::Vector3d& o, const Eigen::Vector3d& d,
                    const Eigen::Vector3d& c, double R) {
    const Eigen::Vector3d b = o - c;
    const double t = -b.dot(d);
    const double disc = t * t - (b.dot(b) - R * R);
    return (disc > 0.0) ? 2.0 * std::sqrt(disc) : 0.0;
}

// The reversed integral for a TRANSPARENT uniform sphere source:
//   eps = sum_i omega_w_i * p_i * L_i / V
double line_set_sphere_eff(const EtendueLineSet& set, const Eigen::Vector3d& c, double R,
                           double energy_keV) {
    std::vector<double> p;
    line_interaction_probabilities(set.q, energy_keV, MuChoice::Total, p);
    const double V = 4.0 / 3.0 * kPi * R * R * R;
    double sum = 0.0;
    for (size_t i = 0; i < set.q.rays.size(); ++i)
        sum += set.q.rays[i].omega_w * p[i] * sphere_chord(set.origin[i], set.dir[i], c, R);
    return sum / V;
}

// The forward integral: the per-point aperture kernel averaged over Halton
// points uniform in the same sphere.
double point_kernel_sphere_eff(const Geometry& g, const Eigen::Vector3d& c, double R,
                               double energy_keV, int n_points, int n_rays) {
    double sum = 0.0;
    int used = 0;
    for (uint64_t i = 0; used < n_points; ++i) {
        const Eigen::Vector3d u(2.0 * halton(i, 2) - 1.0, 2.0 * halton(i, 3) - 1.0,
                                2.0 * halton(i, 5) - 1.0);
        if (u.squaredNorm() > 1.0) continue;
        const ApertureQuadrature q = make_aperture_quadrature(g, c + R * u, n_rays);
        sum += q.interaction_omega(energy_keV, MuChoice::Total);
        ++used;
    }
    return sum / n_points;
}

}  // namespace

// The sampling measure: with directions over the full sphere and the source
// side restricted to z < 0, the etendue is  pi * A_front + (pi/2) * A_side
// (front: the whole inward hemisphere; side wall: the half of its hemisphere
// that points toward z < 0).
BOOST_AUTO_TEST_CASE(etendue_measure_full_sphere) {
    const Geometry g = bare_nai_3x3();
    const double R = 3.81, L = 7.62;
    const EtendueLineSet set = build_etendue_lines(g, {0.0, 0.0, -1.0}, -1.0, 1 << 17);
    const double expect = kPi * (kPi * R * R) + 0.5 * kPi * (2.0 * kPi * R * L);
    BOOST_TEST_MESSAGE("etendue " << set.total_etendue << " vs analytic " << expect
                       << " (" << 100.0 * (set.total_etendue / expect - 1.0) << "%)");
    BOOST_CHECK_CLOSE(set.total_etendue, expect, 1.5 /*percent*/);

    // Every kept line: origin on the hull, unit direction into the crystal,
    // toward +z, sorted segments, an active chord, and a positive weight.
    BOOST_REQUIRE(!set.q.rays.empty());
    BOOST_REQUIRE_EQUAL(set.q.rays.size(), set.origin.size());
    BOOST_REQUIRE_EQUAL(set.q.rays.size(), set.dir.size());
    for (size_t i = 0; i < set.q.rays.size(); ++i) {
        const Eigen::Vector3d& o = set.origin[i];
        const double r = std::hypot(o.x(), o.y());
        BOOST_CHECK(r <= R + 1e-9);
        BOOST_CHECK(o.z() >= -1e-9 && o.z() <= L + 1e-9);
        BOOST_CHECK(std::abs(r - R) < 1e-9 || std::abs(o.z()) < 1e-9);   // side wall or front face
        BOOST_CHECK_CLOSE(set.dir[i].norm(), 1.0, 1e-9);
        BOOST_CHECK(set.dir[i].z() > 0.0);
        BOOST_CHECK(set.q.rays[i].omega_w > 0.0f);
        BOOST_CHECK(set.q.rays[i].active_len > 0.0f);
        BOOST_CHECK(!set.q.rays[i].segs.empty());
    }
    // Sum of the kept weights * 4pi is the etendue that survived the active-chord drop; for a
    // bare crystal every hull line has an active chord, so it is the whole measure.
    double kept = 0.0;
    for (const KernelRay& r : set.q.rays) kept += r.omega_w;
    BOOST_CHECK_CLOSE(kept * 4.0 * kPi, set.total_etendue, 1e-3);
}

// Reversed integral == forward integral for a transparent sphere source, on
// axis and off axis (side wall in play), at two energies.
BOOST_AUTO_TEST_CASE(transparent_sphere_source_identity) {
    const Geometry g = bare_nai_3x3();
    struct Case { Eigen::Vector3d c; double R; };
    const std::vector<Case> cases = {
        {{0.0, 0.0, -5.0}, 1.0},      // on axis, 5 cm
        {{3.0, 0.0, -3.0}, 1.5},      // off axis, near the rim
        {{0.0, 0.0, -1.5}, 1.2},      // contact
    };
    for (const Case& cs : cases) {
        const EtendueLineSet set = build_etendue_lines(g, cs.c, cs.R, 1 << 16);
        for (const double E : {122.0, 662.0}) {
            const double line = line_set_sphere_eff(set, cs.c, cs.R, E);
            const double point = point_kernel_sphere_eff(g, cs.c, cs.R, E, 512, 2048);
            BOOST_TEST_MESSAGE("sphere c=(" << cs.c.transpose() << ") R=" << cs.R << " E=" << E
                               << ": lines " << line << " vs points " << point << " ("
                               << 100.0 * (line / point - 1.0) << "%), " << set.q.rays.size()
                               << " lines kept of " << set.q.n_rays_total);
            BOOST_CHECK_CLOSE(line, point, 1.0 /*percent*/);
        }
    }
}

// A box crystal through the same identity (five-face hull, can layers).
BOOST_AUTO_TEST_CASE(box_crystal_identity) {
    const Geometry g = boxed_czt();
    const Eigen::Vector3d c(0.8, -0.4, -2.0);
    const double R = 0.8;
    const EtendueLineSet set = build_etendue_lines(g, c, R, 1 << 16);
    for (const double E : {122.0, 662.0}) {
        const double line = line_set_sphere_eff(set, c, R, E);
        const double point = point_kernel_sphere_eff(g, c, R, E, 512, 2048);
        BOOST_TEST_MESSAGE("box: E=" << E << ": lines " << line << " vs points " << point << " ("
                           << 100.0 * (line / point - 1.0) << "%)");
        BOOST_CHECK_CLOSE(line, point, 1.5 /*percent*/);
    }
}

// Determinism, and a disjoint index range gives an independent estimate of the
// same measure.
BOOST_AUTO_TEST_CASE(deterministic_and_offset_sets) {
    const Geometry g = bare_nai_3x3();
    const Eigen::Vector3d c(0.0, 0.0, -4.0);
    const EtendueLineSet a = build_etendue_lines(g, c, 2.0, 4096);
    const EtendueLineSet b = build_etendue_lines(g, c, 2.0, 4096);
    BOOST_REQUIRE_EQUAL(a.q.rays.size(), b.q.rays.size());
    for (size_t i = 0; i < a.q.rays.size(); ++i) {
        BOOST_CHECK_EQUAL(a.q.rays[i].omega_w, b.q.rays[i].omega_w);
        BOOST_CHECK_EQUAL(a.origin[i], b.origin[i]);
        BOOST_CHECK_EQUAL(a.dir[i], b.dir[i]);
    }
    const EtendueLineSet o = build_etendue_lines(g, c, 2.0, 4096, 4096);
    BOOST_CHECK(o.q.rays.size() > 0);
    BOOST_CHECK(o.origin.front() != a.origin.front());
    BOOST_CHECK_CLOSE(o.total_etendue, a.total_etendue, 3.0);
}

// The stored-response accessor is parallel to the rays and reproduces the
// compacted fep_ray_weights sum exactly.
BOOST_AUTO_TEST_CASE(stored_response_probabilities_match_ray_weights) {
    const GeometryDescriptor gd = nai3x3_descriptor();
    const Eigen::Vector3d ref(0.0, 0.0, -25.0);
    std::shared_ptr<DetectorResponse> resp = make_transfer_response(gd, make_anchor(), ref);
    BOOST_REQUIRE(resp);

    const EtendueLineSet set = build_etendue_lines(resp->geometry(), {0.0, 0.0, -3.0}, 2.0, 8192);
    for (const double E : {60.0, 300.0, 1000.0}) {
        std::vector<double> p;
        resp->fep_line_probabilities(E, set.q, p);
        BOOST_REQUIRE_EQUAL(p.size(), set.q.rays.size());
        double sum_p = 0.0;
        for (size_t i = 0; i < p.size(); ++i) {
            BOOST_CHECK(p[i] > 0.0 && p[i] <= 1.0);
            sum_p += set.q.rays[i].omega_w * p[i];
        }
        std::vector<double> w;
        std::vector<Eigen::Vector3d> dirs;
        resp->fep_ray_weights(E, set.q, w, dirs);
        double sum_w = 0.0;
        for (const double v : w) sum_w += v;
        BOOST_CHECK_CLOSE(sum_p, sum_w, 1e-9);

        // And the live-cross-section kernel agrees with the stored tables to table precision.
        std::vector<double> live;
        line_interaction_probabilities(set.q, E, MuChoice::Total, live);
        double sum_live = 0.0;
        for (size_t i = 0; i < live.size(); ++i) sum_live += set.q.rays[i].omega_w * live[i];
        BOOST_CHECK_CLOSE(sum_live, sum_p, 0.2);
    }
}
