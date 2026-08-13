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

/// @file RayTrace.h
/// @brief Full ray tracing through the nested detector geometry.
///
/// Traces a ray from an external point through the concentric shells
/// (attenuators, dead layer, active crystal), returning an ordered list
/// of path segments with their materials and scoring status.

#include "geometry/Geometry.h"
#include <Eigen/Core>
#include <vector>
#include <random>

namespace ceelo {

// RayTrace is implemented as methods on the Geometry class.
// This header exists to declare any additional helper functions
// that the ray trace implementation uses internally.

/// Compute the total attenuation (exp(-sum(mu*L))) for a ray passing through
/// a list of non-scoring path segments.
/// @param segments  Path segments from Geometry::trace_ray()
/// @param energy_MeV  Photon energy (needed for cross-section lookup)
/// @return  Transmission factor (0 to 1)
double compute_transmission(const std::vector<PathSegment>& segments,
                            double energy_MeV);

/// Compute the total path length through scoring (active) segments.
/// @param segments  Path segments from Geometry::trace_ray()
/// @return  Total active path length in cm
double compute_active_path_length(const std::vector<PathSegment>& segments);

/// Compute transmission through non-scoring segments excluding Rayleigh scattering.
/// Uses mu_pe + mu_cs + mu_pp (no mu_rs) for FEP-only mode, since Rayleigh
/// is elastic and doesn't remove photons from the full-energy peak.
/// @param segments  Path segments from Geometry::trace_ray()
/// @param energy_MeV  Photon energy (needed for cross-section lookup)
/// @return  Transmission factor (0 to 1)
double compute_transmission_no_rayleigh(const std::vector<PathSegment>& segments,
                                         double energy_MeV);

/// Result of stochastic Rayleigh transmission through non-scoring segments.
struct RayleighTransmissionResult {
    double weight;              ///< Product of exp(-mu_no_rayleigh * L) factors
    Eigen::Vector3d position;   ///< Entry point to scoring volume (or last position)
    Eigen::Vector3d direction;  ///< Final direction after any Rayleigh scatters
    bool hit_scoring;           ///< Whether ray still intersects scoring volume
};

/// Compute transmission through non-scoring segments with stochastic Rayleigh.
///
/// For each non-scoring segment:
///   1. Weight by exp(-mu_no_rs * L) (analytical for energy-changing processes)
///   2. Stochastic Rayleigh: sample number of scatters from Poisson(mu_rs * L),
///      apply direction changes, and re-trace from the scatter position.
///
/// @param start_position  Ray origin (source position)
/// @param start_direction Ray direction (normalized)
/// @param energy_keV      Photon energy in keV
/// @param geometry        Detector geometry for ray tracing
/// @param rng             Random number generator
/// @return Transmission result with weight, final position/direction, and hit flag
RayleighTransmissionResult compute_transmission_stochastic_rayleigh(
    const Eigen::Vector3d& start_position,
    const Eigen::Vector3d& start_direction,
    double energy_keV,
    const Geometry& geometry,
    std::mt19937_64& rng);

} // namespace ceelo
