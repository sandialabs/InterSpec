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

#include "geometry/Geometry.h"
#include <Eigen/Core>
#include <optional>

namespace ceelo {

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

} // namespace ceelo
