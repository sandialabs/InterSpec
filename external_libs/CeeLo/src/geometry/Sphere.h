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

/// @file Sphere.h
/// @brief Ray-sphere intersection for spherical shielding shells around point sources.

#include "geometry/Geometry.h"

#include <Eigen/Core>
#include <optional>
#include <cmath>

namespace ceelo {

/// Compute ray intersection with a sphere centered at `center` with radius `radius`.
///
/// Solves |origin + t * direction - center|^2 = radius^2 for t.
/// Returns RayHit{t_enter, t_exit} or nullopt if no intersection.
///
/// @param origin     Ray origin
/// @param direction  Ray direction (must be normalized)
/// @param center     Sphere center
/// @param radius     Sphere radius
/// @return           Hit with t_enter and t_exit, or nullopt if no intersection
inline std::optional<RayHit> intersect_sphere(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    const Eigen::Vector3d& center,
    double radius)
{
    Eigen::Vector3d oc = origin - center;
    // direction is normalized, so a = 1
    double b = oc.dot(direction);
    double c = oc.squaredNorm() - radius * radius;
    double discriminant = b * b - c;

    if (discriminant < 0.0) return std::nullopt;

    double sqrt_disc = std::sqrt(discriminant);
    double t_enter = -b - sqrt_disc;
    double t_exit  = -b + sqrt_disc;

    if (t_exit < 0.0) return std::nullopt;

    return RayHit{t_enter, t_exit};
}

} // namespace ceelo
