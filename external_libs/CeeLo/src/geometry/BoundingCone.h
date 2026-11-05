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

/// @file BoundingCone.h
/// @brief Smallest cone from a point containing a z-aligned cylinder or box - used to aim
///        Monte-Carlo and aperture-quadrature rays at a detector from any position.

#include <Eigen/Core>

namespace ceelo {

/// Smallest cone from `src` (crystal frame, cm) containing a solid, z-aligned cylinder of radius
/// `R` (when `cylinder`) or box of half-widths `hx` x `hy`, spanning z in [z0, z1]: its axis and
/// cos(half-angle).  Returns false - the caller must sample the full sphere - when `src` is inside
/// that envelope, or no cone narrower than 90 degrees can be certified (see the safety note in
/// ResponseKernel.cpp).  make_aperture_quadrature uses it for the active crystal; the Monte-Carlo
/// cone sampler (EfficiencyCalculator::enable_general_cone_sampling) for the outer detector
/// envelope, where a source beside or behind the detector has no cone about the detector axis.
bool envelope_bounding_cone(bool cylinder, double R, double hx, double hy, double z0, double z1,
                            const Eigen::Vector3d& src, Eigen::Vector3d& axis_out,
                            double& cos_alpha_out);

}  // namespace ceelo
