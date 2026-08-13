#ifndef CEELO_IO_SOLID_ANGLE_H
#define CEELO_IO_SOLID_ANGLE_H
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

#include <cmath>

namespace ceelo {

/// On-axis disk solid-angle fraction Omega/4pi of a radius-R face seen from a
/// point at axial distance x in front of it:
///
///     Omega / 4pi = 0.5 * (1 - x / sqrt(x^2 + R^2))
///
/// For x <= 0 (point at or behind the face plane) the on-axis fraction
/// saturates at 0.5 (a full hemisphere). This is the single definition shared by
/// the VPD fitter, the efficiency-grid study, and the solid-angle unit test.
///
/// NOTE: this is deliberately NOT the importance-sampling cone construct inside
/// EfficiencyCalculator.cpp (`solid_angle_h_cm`) -- that is a different sampler
/// and must not be conflated with this geometric disk fraction.
inline double disk_solid_angle_fraction(double x, double R) {
    if (x <= 0.0) return 0.5;
    return 0.5 * (1.0 - x / std::sqrt(x * x + R * R));
}

} // namespace ceelo

#endif // CEELO_IO_SOLID_ANGLE_H
