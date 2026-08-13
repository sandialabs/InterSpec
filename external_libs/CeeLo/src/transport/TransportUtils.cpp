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

#include "transport/TransportUtils.h"

#include <cmath>

namespace ceelo {

namespace {
constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
} // anonymous namespace

InteractionType select_interaction(const MacroscopicXS& xs, double xi) {
    double total = xs.mu_total();
    if (total <= 0.0) return InteractionType::Photoelectric; // fallback

    double threshold = xi * total;
    double cumulative = xs.mu_pe;
    if (threshold < cumulative) return InteractionType::Photoelectric;
    cumulative += xs.mu_cs;
    if (threshold < cumulative) return InteractionType::Compton;
    cumulative += xs.mu_rs;
    if (threshold < cumulative) return InteractionType::Rayleigh;
    if (xs.mu_pp > 0.0) return InteractionType::PairProduction;
    return InteractionType::Compton;  // float rounding fallback
}

double sample_interaction_distance(double mu_total, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double xi = uniform(rng);
    // Avoid log(0) — clamp xi away from 0
    if (xi < 1e-30) xi = 1e-30;
    return -std::log(xi) / mu_total;
}

Eigen::Vector3d sample_isotropic_direction(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double cos_theta = 2.0 * uniform(rng) - 1.0;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = kTwoPi * uniform(rng);
    return Eigen::Vector3d(
        sin_theta * std::cos(phi),
        sin_theta * std::sin(phi),
        cos_theta);
}

} // namespace ceelo
