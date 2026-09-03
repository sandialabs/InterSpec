#pragma once
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

/// @file Cylinder.h
/// @brief Ray-finite-cylinder intersection.
///
/// A finite cylinder is defined by:
///   - Radius R, centered on the z-axis
///   - Axial extent from z_min to z_max
///
/// The intersection algorithm solves the quadratic for the infinite cylinder
/// (x² + y² = R²), then clips against the endcap planes (z = z_min, z = z_max).
///
/// Crystals may also be *bulletized*: the outer front edge is rounded by a
/// quarter-torus fillet instead of meeting the side wall at a sharp 90° corner
/// (common for HPGe).  See FrontFillet below.  Sharp cylinders never touch the
/// fillet code -- callers pass a fillet only when one is configured.

#include "geometry/Geometry.h"
#include <Eigen/Core>
#include <optional>

namespace ceelo {

/// Parameters of a rounded ("bulletized") outer front edge, precomputed from
/// the crystal radius, front-face z and the bulletization radius r_b.
///
/// The fillet is the quarter torus whose ring circle has radius `rho_c` and
/// lies in the plane z = `z_c`, with tube radius `r_b`; it is tangent to the
/// front face at rho = rho_c and to the side wall at z = z_c.  The bulletized
/// solid is the sharp cylinder minus the corner material outside that torus:
///
///     inside  <=>  inside sharp cylinder
///                  AND (rho <= rho_c OR z >= z_c
///                       OR (rho - rho_c)² + (z - z_c)² <= r_b²)
///
/// The solid is convex and a subset of the sharp cylinder, which is what makes
/// the endpoint correction in correct_hit_for_front_fillet() valid.
struct FrontFillet {
    double rho_c;   ///< Ring radius, R - r_b
    double z_c;     ///< Absolute z of the ring plane, z_min + r_b
    double r_b;     ///< Fillet (tube) radius
    double r_b_sq;  ///< r_b², precomputed
};

/// Build the fillet for a cylinder of the given radius and front-face z.
/// Requires 0 < bullet_radius < radius (rho_c > 0 keeps the surface smooth).
inline FrontFillet make_front_fillet(double radius, double z_min,
                                     double bullet_radius) {
    return FrontFillet{radius - bullet_radius, z_min + bullet_radius,
                       bullet_radius, bullet_radius * bullet_radius};
}

/// Compute ray intersection with a finite cylinder centered on the z-axis.
///
/// @param origin     Ray origin
/// @param direction  Ray direction (must be normalized)
/// @param radius     Cylinder radius
/// @param z_min      Bottom of cylinder (front face)
/// @param z_max      Top of cylinder (back face)
/// @return           Hit with t_enter and t_exit, or nullopt if no intersection
std::optional<RayHit> intersect_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double radius,
    double z_min,
    double z_max);

/// Compute ray intersection with a bore hole (subtracted cylinder).
/// Returns up to two active segments when the ray passes through the bore.
///
/// @param origin     Ray origin
/// @param direction  Ray direction (must be normalized)
/// @param outer_radius  Outer cylinder radius
/// @param bore_radius   Bore hole radius
/// @param z_min         Front face of outer cylinder
/// @param z_max         Back face of outer cylinder
/// @param bore_z_start  z where bore starts (= z_max - bore_depth)
/// @param bore_z_end    z where bore ends (= z_max, i.e., back face)
/// @param segments_out  Output: 0, 1, or 2 hit intervals in the active crystal
/// @return Number of active segments (0, 1, or 2)
int intersect_bored_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double outer_radius,
    double bore_radius,
    double z_min,
    double z_max,
    double bore_z_start,
    double bore_z_end,
    RayHit segments_out[2]);

/// Shrink a sharp-cylinder hit interval onto the bulletized solid.
///
/// `hit` must be the interval of the *same* ray against the sharp cylinder of
/// radius `fillet.rho_c + fillet.r_b` whose front face is at
/// `fillet.z_c - fillet.r_b`.  Because the bulletized solid is convex and a
/// subset of that cylinder, an endpoint lying on a surface the two solids
/// share (front face inside rho_c, side wall above z_c, back face) is already
/// exact; only an endpoint inside the removed corner needs to be moved onto
/// the fillet.  Endpoints are handled independently, so interior-origin rays
/// (electron walks, re-traced scatters) work unchanged: the test is on the hit
/// *point*, not on the sign of t.
///
/// @return false if the ray misses the bulletized solid entirely (it clipped
///         only removed corner material), in which case `hit` is unspecified.
bool correct_hit_for_front_fillet(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    const FrontFillet& fillet,
    RayHit& hit);

/// Ray intersection with a cylinder whose outer front edge is bulletized.
/// Equivalent to intersect_cylinder() followed by correct_hit_for_front_fillet().
std::optional<RayHit> intersect_bulletized_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double radius,
    double z_min,
    double z_max,
    const FrontFillet& fillet);

/// intersect_bored_cylinder() with optional rounded features.
///
/// `direction` must be normalized -- the rounded bore tip solves the ray-sphere
/// quadratic assuming a unit direction.
///
/// @param fillet         Front-edge fillet; pass r_b == 0 for a sharp edge.
/// @param rounded_bore_tip  If true the bore's closed end is a hemisphere of
///                       radius `bore_radius` whose apex sits at
///                       `bore_z_start` (a round-tipped drill: the stated bore
///                       depth is preserved), rather than a flat bottom.
int intersect_shaped_bored_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double outer_radius,
    double bore_radius,
    double z_min,
    double z_max,
    double bore_z_start,
    double bore_z_end,
    const FrontFillet& fillet,
    bool rounded_bore_tip,
    RayHit segments_out[2]);

} // namespace ceelo
