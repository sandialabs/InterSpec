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

#include "geometry/Box.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace ceelo {

std::optional<RayHit> intersect_box(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double half_x,
    double half_y,
    double z_min,
    double z_max)
{
    // Slab method for AABB intersection.
    // For each axis i: compute t_min_i and t_max_i for the slab.
    // Overall: t_enter = max(t_min_x, t_min_y, t_min_z)
    //          t_exit  = min(t_max_x, t_max_y, t_max_z)
    // Hit if t_enter <= t_exit and t_exit > 0.

    // Box bounds
    double box_min[3] = {-half_x, -half_y, z_min};
    double box_max[3] = { half_x,  half_y, z_max};

    double t_enter = -std::numeric_limits<double>::max();
    double t_exit  =  std::numeric_limits<double>::max();

    for (int i = 0; i < 3; ++i) {
        double o = origin[i];
        double d = direction[i];

        if (std::abs(d) < 1.0e-12) {
            // Ray is parallel to this slab
            if (o < box_min[i] || o > box_max[i]) {
                return std::nullopt; // Outside the slab, no intersection
            }
            // Inside the slab for this axis, don't constrain t_enter/t_exit
            continue;
        }

        // Pre-compute 1/d for this axis
        double inv_d = 1.0 / d;

        double t1 = (box_min[i] - o) * inv_d;
        double t2 = (box_max[i] - o) * inv_d;

        if (t1 > t2) std::swap(t1, t2);

        t_enter = std::max(t_enter, t1);
        t_exit  = std::min(t_exit,  t2);

        if (t_enter > t_exit) {
            return std::nullopt; // No intersection
        }
    }

    // Must have forward intersection
    if (t_exit <= 0.0) {
        return std::nullopt;
    }

    return RayHit{t_enter, t_exit};
}

} // namespace ceelo
