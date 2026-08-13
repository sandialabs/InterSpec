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

/// @file Box.h
/// @brief Ray-AABB (axis-aligned bounding box) intersection using the slab method.
///
/// A box is defined by its center (always on the z-axis at x=0, y=0) and
/// half-extents (half_x, half_y) perpendicular to the z-axis, with
/// z-extent from z_min to z_max.

#include "geometry/Geometry.h"
#include <Eigen/Core>
#include <optional>

namespace ceelo {

/// Compute ray intersection with an axis-aligned box.
///
/// The box spans:
///   x: [-half_x, +half_x]
///   y: [-half_y, +half_y]
///   z: [z_min, z_max]
///
/// Uses the slab method: compute the intersection interval for each axis
/// independently, then intersect the three intervals.
///
/// @param origin     Ray origin
/// @param direction  Ray direction (must be normalized)
/// @param half_x     Half-width in x
/// @param half_y     Half-width in y
/// @param z_min      Front face z
/// @param z_max      Back face z
/// @return           Hit with t_enter and t_exit, or nullopt if no intersection
std::optional<RayHit> intersect_box(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double half_x,
    double half_y,
    double z_min,
    double z_max);

} // namespace ceelo
