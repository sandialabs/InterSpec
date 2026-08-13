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

/// @file TransportUtils.h
/// @brief Shared transport utility functions used by both PhotonTransport
///        and SourceGeometry for interaction sampling.

#include "materials/Material.h"
#include "transport/PhotonTransport.h"

#include <Eigen/Core>
#include <random>

namespace ceelo {

/// Select interaction type given macroscopic cross-sections and a random number.
InteractionType select_interaction(const MacroscopicXS& xs, double xi);

/// Sample interaction distance: s = -ln(xi) / Σ_total
double sample_interaction_distance(double mu_total, std::mt19937_64& rng);

/// Sample a direction uniformly over the full 4π sphere.
Eigen::Vector3d sample_isotropic_direction(std::mt19937_64& rng);

} // namespace ceelo
