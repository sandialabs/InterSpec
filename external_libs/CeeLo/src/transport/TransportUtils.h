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

struct FluorescenceData;
struct LFluorescenceData;

/// How a K-shell vacancy filled: a K x-ray of `line_keV` (with the rest of
/// the binding energy, `local_keV`, going to the L-shell Auger electrons), or
/// an Auger cascade (`line_keV` = 0, the whole binding energy local).
struct KRelaxation {
    double line_keV = 0.0;
    double local_keV = 0.0;
    /// EADL subshell whose electron filled the vacancy (EadlSubshell numbering),
    /// 0 for an Auger cascade. `k_line_family()` gives K-alpha / K-beta.
    uint8_t filling_subshell = 0;
};

/// Fill a K vacancy of the element `fl` describes: an x-ray with probability
/// omega_K, the line drawn by its relative intensity.  Draws the yield then,
/// for an x-ray, the line (two uniforms); the caller samples the direction.
KRelaxation sample_k_relaxation(const FluorescenceData& fl, std::mt19937_64& rng);

/// How an L-subshell vacancy filled: a retained L x-ray (`line_keV`, lines at or
/// above the 10 keV transport cut) or none (Auger and Coster-Kronig electrons and
/// sub-cut lines). `local_keV` is the starting subshell's binding energy less the
/// x-ray: everything that stays on the spot.
struct LRelaxation {
    double line_keV = 0.0;
    double local_keV = 0.0;
};

/// True when any L subshell of the element has a retained line; otherwise an L
/// vacancy relaxes locally and callers skip the L draws altogether.
bool has_l_lines(const LFluorescenceData* fl);

/// Fill a vacancy in L subshell `subshell` (0 = L1, 1 = L2, 2 = L3): radiate a
/// retained line with that subshell's yield, else move it by Coster-Kronig
/// transfer (f13 then f12 from L1, f23 from L2) and try again, else stop (Auger).
/// The cascade vacancy model draws in the same order.
LRelaxation sample_l_relaxation(const LFluorescenceData& fl, int subshell,
                                std::mt19937_64& rng);

/// The L subshell (0..2) a photoabsorption outside the K shell leaves its vacancy
/// in, -1 for M and beyond, from the element's tabulated L1/L2/L3 shares of the
/// photoelectric cross section at `energy_keV`; `f_K` is the K share at that
/// energy (0 below the K edge). One uniform draw. Returns -1 without drawing when
/// the element has no share table.
int sample_l_photoelectric_subshell(const LFluorescenceData& fl, double energy_keV,
                                    double f_K, std::mt19937_64& rng);

/// Probability-weighted mean energy and total relative intensity of an element's
/// K-alpha and K-beta lines (families from the lines' filling subshells).
struct KLineFamilyMeans {
    double alpha_keV = 0.0;
    double beta_keV = 0.0;
    double alpha_probability = 0.0;   ///< sum of K-alpha line probabilities
    double beta_probability = 0.0;    ///< sum of K-beta line probabilities
};
KLineFamilyMeans k_line_family_means(const FluorescenceData& fl);

} // namespace ceelo
