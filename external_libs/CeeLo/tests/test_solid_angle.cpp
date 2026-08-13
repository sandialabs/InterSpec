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

#define BOOST_TEST_MODULE SolidAngleTests
#include <boost/test/unit_test.hpp>

#include "geometry/Geometry.h"
#include "io/SolidAngle.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <random>

using namespace ceelo;

namespace {
constexpr double kPi = 3.14159265358979323846;

/// Analytical solid angle fraction for an on-axis point source at distance d
/// from a disk of radius R (shared definition in io/SolidAngle.h):
///   Omega / 4pi = (1 - d / sqrt(d^2 + R^2)) / 2
double point_source_solid_angle_fraction(double d, double R) {
    return disk_solid_angle_fraction(d, R);
}

/// Monte Carlo geometric hit fraction: fire N isotropic rays from source_pos,
/// count fraction that hit the detector (trace_ray returns non-empty segments).
double mc_hit_fraction(const Geometry& geom,
                       const Eigen::Vector3d& source_pos,
                       uint64_t N, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    uint64_t n_hit = 0;

    for (uint64_t i = 0; i < N; ++i) {
        // Sample isotropic direction (full 4pi)
        double cos_theta = 2.0 * uniform(rng) - 1.0;
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        double phi = 2.0 * kPi * uniform(rng);

        Eigen::Vector3d dir(sin_theta * std::cos(phi),
                            sin_theta * std::sin(phi),
                            cos_theta);

        auto segments = geom.trace_ray(source_pos, dir);
        if (!segments.empty()) {
            ++n_hit;
        }
    }

    return static_cast<double>(n_hit) / static_cast<double>(N);
}
} // anonymous namespace


// ============================================================
//  Point Source Solid Angle Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(PointSourceSolidAngle)

BOOST_AUTO_TEST_CASE(on_axis_point_source_15cm) {
    // 3"x3" NaI bare detector (R=3.81 cm, L=7.62 cm)
    // Point source on axis at z=-15 cm
    //
    // Analytical: Omega/4pi = (1 - 15/sqrt(225 + 14.5161)) / 2
    //           = (1 - 15/15.4796) / 2
    //           = 0.01548
    //
    // The detector is a finite cylinder, not a disk, so the MC result will be
    // slightly larger than the disk formula (the cylinder side catches some
    // extra photons). The disk formula gives a lower bound.

    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    Eigen::Vector3d source_pos(0.0, 0.0, -15.0);

    const uint64_t N = 2000000;
    std::mt19937_64 rng(42);
    double mc_frac = mc_hit_fraction(geom, source_pos, N, rng);

    // Analytical disk approximation (lower bound)
    double analytic_disk = point_source_solid_angle_fraction(15.0, 3.81);

    // MC fraction should be close to (but >= ) the disk solid angle
    // For a 3" cylinder at 15 cm, the disk subtends ~1.55% and the full
    // cylinder subtends ~1.6% (side contribution is small at this distance).

    double mc_sigma = std::sqrt(mc_frac * (1.0 - mc_frac) / N);

    BOOST_TEST_MESSAGE("Point source 15cm on-axis:");
    BOOST_TEST_MESSAGE("  MC hit fraction:     " << mc_frac << " +/- " << mc_sigma);
    BOOST_TEST_MESSAGE("  Analytic disk Omega/4pi: " << analytic_disk);

    // MC should exceed the disk lower bound (within 3-sigma statistical noise)
    BOOST_CHECK_GT(mc_frac, analytic_disk - 3.0 * mc_sigma);

    // MC should not exceed the disk by more than ~10% (cylinder side is small)
    BOOST_CHECK_LT(mc_frac, analytic_disk * 1.15);

    // Overall agreement: within 10% relative
    double rel_diff = std::abs(mc_frac - analytic_disk) / analytic_disk;
    BOOST_CHECK_LT(rel_diff, 0.10);
}

BOOST_AUTO_TEST_CASE(on_axis_point_source_5cm_close) {
    // Closer source at 5 cm — larger solid angle, greater cylinder-side contribution.
    // Analytic disk: Omega/4pi = (1 - 5/sqrt(25+14.5161))/2 = (1 - 5/6.288)/2 = 0.1024

    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    Eigen::Vector3d source_pos(0.0, 0.0, -5.0);

    const uint64_t N = 1000000;
    std::mt19937_64 rng(123);
    double mc_frac = mc_hit_fraction(geom, source_pos, N, rng);

    double analytic_disk = point_source_solid_angle_fraction(5.0, 3.81);
    double mc_sigma = std::sqrt(mc_frac * (1.0 - mc_frac) / N);

    BOOST_TEST_MESSAGE("Point source 5cm on-axis:");
    BOOST_TEST_MESSAGE("  MC hit fraction:     " << mc_frac << " +/- " << mc_sigma);
    BOOST_TEST_MESSAGE("  Analytic disk Omega/4pi: " << analytic_disk);

    // MC must exceed disk lower bound
    BOOST_CHECK_GT(mc_frac, analytic_disk - 3.0 * mc_sigma);

    // At 5 cm the cylinder side subtends more, so allow up to 20% above disk
    BOOST_CHECK_LT(mc_frac, analytic_disk * 1.25);
}

BOOST_AUTO_TEST_CASE(off_axis_point_source_reduces_solid_angle) {
    // A point source moved off-axis should subtend a smaller solid angle
    // than the same source on-axis (for a cylindrical detector).

    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    const uint64_t N = 500000;

    // On-axis at 15 cm
    std::mt19937_64 rng_on(77);
    double frac_on = mc_hit_fraction(geom, {0.0, 0.0, -15.0}, N, rng_on);

    // Off-axis at 15 cm in z but shifted 5 cm in x
    std::mt19937_64 rng_off(78);
    double frac_off = mc_hit_fraction(geom, {5.0, 0.0, -15.0}, N, rng_off);

    BOOST_TEST_MESSAGE("On-axis  15cm hit fraction: " << frac_on);
    BOOST_TEST_MESSAGE("Off-axis 15cm hit fraction: " << frac_off);

    // Off-axis solid angle should be smaller
    BOOST_CHECK_LT(frac_off, frac_on);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Extended Source Solid Angle Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(ExtendedSourceSolidAngle)

BOOST_AUTO_TEST_CASE(cylindrical_source_geometric_hit_fraction) {
    // Cylindrical source (R=1.5 cm, half-length=2 cm) centered at z=-15 cm,
    // coaxial with a 3"x3" NaI detector.
    //
    // Sample source positions uniformly from the cylinder, fire isotropic rays,
    // count the fraction that hit the detector.
    //
    // Since the source is centered on axis, the average solid angle should be
    // close to (but slightly less than) the on-axis point source at z=-15,
    // because some source points are farther away and off-axis.
    //
    // This is a purely geometric test — no physics interaction, just ray tracing.

    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    const uint64_t N = 1000000;
    std::mt19937_64 rng(999);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    const double src_radius = 1.5;
    const double src_half_length = 2.0;
    const Eigen::Vector3d src_center(0.0, 0.0, -15.0);

    uint64_t n_hit = 0;

    for (uint64_t i = 0; i < N; ++i) {
        // Sample source position in cylinder
        double r = src_radius * std::sqrt(uniform(rng));
        double phi_src = 2.0 * kPi * uniform(rng);
        double z_local = src_half_length * (2.0 * uniform(rng) - 1.0);
        Eigen::Vector3d pos(r * std::cos(phi_src),
                            r * std::sin(phi_src),
                            src_center.z() + z_local);

        // Sample isotropic direction
        double cos_theta = 2.0 * uniform(rng) - 1.0;
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        double phi = 2.0 * kPi * uniform(rng);
        Eigen::Vector3d dir(sin_theta * std::cos(phi),
                            sin_theta * std::sin(phi),
                            cos_theta);

        auto segments = geom.trace_ray(pos, dir);
        if (!segments.empty()) {
            ++n_hit;
        }
    }

    double mc_frac = static_cast<double>(n_hit) / static_cast<double>(N);
    double mc_sigma = std::sqrt(mc_frac * (1.0 - mc_frac) / N);

    // Reference: on-axis point source at center of cylinder
    double analytic_center = point_source_solid_angle_fraction(15.0, 3.81);

    // The extended source averages over positions, some closer and some
    // farther, and some off-axis.  The average should be close to the
    // center point value.  The z-spread (-2 to +2 cm around -15 cm)
    // is small relative to the distance, so the effect is modest.
    // The radial offset (0 to 1.5 cm) also has a modest effect at 15 cm.

    BOOST_TEST_MESSAGE("Cylindrical source (R=1.5, h=4, z=-15):");
    BOOST_TEST_MESSAGE("  MC hit fraction:     " << mc_frac << " +/- " << mc_sigma);
    BOOST_TEST_MESSAGE("  Center point Omega/4pi: " << analytic_center);

    // MC should be positive and in a reasonable range
    BOOST_CHECK_GT(mc_frac, 0.005);  // Must detect something
    BOOST_CHECK_LT(mc_frac, 0.05);   // Can't be more than a few percent

    // Should be within ~20% of the center-point value
    // (averaged over the cylinder, slightly less due to 1/r^2 weighting)
    double rel_diff = std::abs(mc_frac - analytic_center) / analytic_center;
    BOOST_CHECK_LT(rel_diff, 0.20);
}

BOOST_AUTO_TEST_CASE(wider_source_at_fixed_z_reduces_average_solid_angle) {
    // A wider cylindrical source (larger radius, but SAME z-extent)
    // at a fixed center distance should have a lower average geometric
    // hit fraction, because off-axis points subtend a smaller solid angle.
    // Use a thin disk (half_length = 0.01 cm) to isolate the radial effect.

    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    const uint64_t N = 500000;
    const Eigen::Vector3d center(0.0, 0.0, -15.0);
    const double half_length = 0.01; // thin disk

    auto sample_cyl = [&](double radius, std::mt19937_64& rng) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        uint64_t n_hit = 0;
        for (uint64_t i = 0; i < N; ++i) {
            double r = radius * std::sqrt(uniform(rng));
            double phi_src = 2.0 * kPi * uniform(rng);
            double z_local = half_length * (2.0 * uniform(rng) - 1.0);
            Eigen::Vector3d pos(r * std::cos(phi_src),
                                r * std::sin(phi_src),
                                center.z() + z_local);

            double cos_theta = 2.0 * uniform(rng) - 1.0;
            double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
            double phi = 2.0 * kPi * uniform(rng);
            Eigen::Vector3d dir(sin_theta * std::cos(phi),
                                sin_theta * std::sin(phi),
                                cos_theta);

            if (!geom.trace_ray(pos, dir).empty()) ++n_hit;
        }
        return static_cast<double>(n_hit) / N;
    };

    // Narrow disk: R=0.1 cm (essentially a point)
    std::mt19937_64 rng_narrow(200);
    double frac_narrow = sample_cyl(0.1, rng_narrow);

    // Wide disk: R=8 cm (wider than detector, lots of off-axis points)
    std::mt19937_64 rng_wide(201);
    double frac_wide = sample_cyl(8.0, rng_wide);

    BOOST_TEST_MESSAGE("Narrow disk (R=0.1) hit fraction: " << frac_narrow);
    BOOST_TEST_MESSAGE("Wide disk   (R=8.0) hit fraction: " << frac_wide);

    // The narrow source should have a higher average hit fraction
    BOOST_CHECK_GT(frac_narrow, frac_wide);
}

BOOST_AUTO_TEST_SUITE_END()
