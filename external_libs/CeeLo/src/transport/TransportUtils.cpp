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
#include "cross_sections/CrossSectionData.h"

#include <algorithm>
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

KRelaxation sample_k_relaxation(const FluorescenceData& fl, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    KRelaxation out;
    if (uniform(rng) >= fl.fluorescence_yield) {
        out.local_keV = fl.k_edge_keV;  // Auger electrons, sub-mm range
        return out;
    }
    const double xi_line = uniform(rng);
    double cum = 0.0;
    int chosen = 0;
    for (int li = 0; li < fl.num_lines; ++li) {
        cum += fl.line_probability[li];
        if (xi_line <= cum) {
            chosen = li;
            break;
        }
    }
    const float line = fl.line_energy_keV[chosen];
    out.line_keV = static_cast<double>(line);
    out.local_keV = fl.k_edge_keV - static_cast<double>(line);
    out.filling_subshell = fl.line_filling_subshell[chosen];
    return out;
}

bool has_l_lines(const LFluorescenceData* fl) {
    if (!fl) return false;
    for (const LSubshellFluor& sub : fl->sub)
        if (sub.num_lines > 0) return true;
    return false;
}

LRelaxation sample_l_relaxation(const LFluorescenceData& fl, int subshell,
                                std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    LRelaxation out;
    out.local_keV = fl.sub[subshell].edge_keV;
    const float* ck = fl.coster_kronig;  // {f12, f13, f23}
    int s = subshell;
    for (int guard = 0; guard < 3; ++guard) {
        const LSubshellFluor& sub = fl.sub[s];
        const double r = uniform(rng);
        if (r < sub.fluorescence_yield) {
            const double xi = uniform(rng);
            double cum = 0.0;
            int chosen = 0;
            for (int li = 0; li < sub.num_lines; ++li) {
                cum += sub.line_probability[li];
                if (xi <= cum) {
                    chosen = li;
                    break;
                }
            }
            if (sub.num_lines > 0) {
                out.line_keV = sub.line_energy_keV[chosen];
                out.local_keV -= out.line_keV;
            }
            return out;
        }
        const double w = sub.fluorescence_yield;
        if (s == 0) {
            if (r < w + ck[1]) { s = 2; continue; }          // L1 -> L3
            if (r < w + ck[1] + ck[0]) { s = 1; continue; }  // L1 -> L2
        } else if (s == 1) {
            if (r < w + ck[2]) { s = 2; continue; }          // L2 -> L3
        }
        return out;  // Auger, or an L3 vacancy that did not radiate
    }
    return out;
}

int sample_l_photoelectric_subshell(const LFluorescenceData& fl, double energy_keV,
                                    double f_K, std::mt19937_64& rng) {
    if (fl.num_pe_share_nodes == 0 || f_K >= 1.0) return -1;
    double share[3];
    l_photoelectric_shares(fl, energy_keV, share);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double r = uniform(rng) * (1.0 - f_K);
    double cum = 0.0;
    for (int s = 0; s < 3; ++s) {
        cum += share[s];
        if (r < cum) return s;
    }
    return -1;
}

KLineFamilyMeans k_line_family_means(const FluorescenceData& fl) {
    KLineFamilyMeans out;
    double sums[2] = {0.0, 0.0};
    double weights[2] = {0.0, 0.0};
    for (int li = 0; li < fl.num_lines; ++li) {
        const int family = static_cast<int>(k_line_family(fl.line_filling_subshell[li]));
        sums[family] += static_cast<double>(fl.line_energy_keV[li]) * fl.line_probability[li];
        weights[family] += fl.line_probability[li];
    }
    out.alpha_probability = weights[0];
    out.beta_probability = weights[1];
    out.alpha_keV = weights[0] > 0.0 ? sums[0] / weights[0] : 0.0;
    out.beta_keV = weights[1] > 0.0 ? sums[1] / weights[1] : 0.0;
    return out;
}

} // namespace ceelo
