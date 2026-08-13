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

#include "geometry/Cylinder.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace ceelo {

namespace {

constexpr double kEpsilon = 1.0e-12;
constexpr double kInfinity = std::numeric_limits<double>::max();

/// Clip a ray interval [t_enter, t_exit] (currently intersecting an infinite cylinder)
/// against the z-slab [z_min, z_max].  Also checks endcap intersections.
/// Returns false if the clipped interval is empty.
bool clip_to_z_range(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double z_min, double z_max,
    double& t_enter, double& t_exit)
{
    double oz = origin.z();
    double dz = direction.z();

    if (std::abs(dz) < kEpsilon) {
        // Ray is parallel to endcaps
        if (oz < z_min || oz > z_max) {
            return false; // Outside the z-range entirely
        }
        // Inside z-range, t_enter/t_exit unchanged
        return t_enter < t_exit;
    }

    // Intersect with z=z_min and z=z_max planes
    double t_z_lo = (z_min - oz) / dz;
    double t_z_hi = (z_max - oz) / dz;
    if (t_z_lo > t_z_hi) std::swap(t_z_lo, t_z_hi);

    // Intersect the cylinder interval with the z-slab interval
    t_enter = std::max(t_enter, t_z_lo);
    t_exit  = std::min(t_exit,  t_z_hi);

    return t_enter < t_exit;
}

} // anonymous namespace

std::optional<RayHit> intersect_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double radius,
    double z_min,
    double z_max)
{
    // Solve the quadratic for the infinite cylinder x² + y² = R²
    // Ray: P(t) = O + t*D
    // (Ox + t*Dx)² + (Oy + t*Dy)² = R²
    // a*t² + b*t + c = 0
    double ox = origin.x(), oy = origin.y();
    double dx = direction.x(), dy = direction.y();

    double a = dx * dx + dy * dy;
    double b = 2.0 * (ox * dx + oy * dy);
    double c = ox * ox + oy * oy - radius * radius;

    double t_enter, t_exit;

    if (a < kEpsilon) {
        // Ray is parallel to the cylinder axis (along z)
        if (c > kEpsilon) {
            return std::nullopt; // Outside the cylinder, moving parallel
        }
        // Inside (or on surface of) the cylinder
        t_enter = -kInfinity;
        t_exit  =  kInfinity;
    } else {
        double discriminant = b * b - 4.0 * a * c;
        if (discriminant < 0.0) {
            return std::nullopt; // No intersection with infinite cylinder
        }

        double sqrt_disc = std::sqrt(discriminant);
        double inv_2a = 0.5 / a;
        t_enter = (-b - sqrt_disc) * inv_2a;
        t_exit  = (-b + sqrt_disc) * inv_2a;
    }

    // Clip to the finite z-range
    if (!clip_to_z_range(origin, direction, z_min, z_max, t_enter, t_exit)) {
        return std::nullopt;
    }

    // Must have some forward intersection
    if (t_exit <= 0.0) {
        return std::nullopt;
    }

    return RayHit{t_enter, t_exit};
}


int intersect_bored_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double outer_radius,
    double bore_radius,
    double z_min,
    double z_max,
    double bore_z_start,
    double bore_z_end,
    RayHit segments_out[2])
{
    // Step 1: Intersect ray with the outer cylinder
    auto outer_hit = intersect_cylinder(origin, direction, outer_radius, z_min, z_max);
    if (!outer_hit || !outer_hit->valid()) {
        return 0; // Ray misses the outer crystal entirely
    }

    // Step 2: Intersect ray with the bore cylinder (which goes from bore_z_start to bore_z_end)
    auto bore_hit = intersect_cylinder(origin, direction, bore_radius, bore_z_start, bore_z_end);

    if (!bore_hit || !bore_hit->valid()) {
        // Ray doesn't hit the bore — single continuous segment through crystal
        segments_out[0] = *outer_hit;
        return 1;
    }

    // Step 3: Subtract the bore from the outer crystal.
    // The outer crystal interval is [t_enter_out, t_exit_out].
    // The bore interval is [t_enter_bore, t_exit_bore].
    // Active segments = outer_interval minus bore_interval.

    double t_eo = std::max(outer_hit->t_enter, 0.0);
    double t_xo = outer_hit->t_exit;
    double t_eb = bore_hit->t_enter;
    double t_xb = bore_hit->t_exit;

    // Clamp bore interval to outer interval
    t_eb = std::max(t_eb, t_eo);
    t_xb = std::min(t_xb, t_xo);

    if (t_eb >= t_xb) {
        // Bore interval is empty within the outer crystal — single segment
        segments_out[0] = *outer_hit;
        return 1;
    }

    // We have bore [t_eb, t_xb] inside outer [t_eo, t_xo].
    // Active segments are: [t_eo, t_eb] and [t_xb, t_xo]
    int count = 0;

    if (t_eo < t_eb - kEpsilon) {
        segments_out[count].t_enter = outer_hit->t_enter;
        segments_out[count].t_exit  = t_eb;
        count++;
    }

    if (t_xb < t_xo - kEpsilon) {
        segments_out[count].t_enter = t_xb;
        segments_out[count].t_exit  = t_xo;
        count++;
    }

    return count;
}

} // namespace ceelo
