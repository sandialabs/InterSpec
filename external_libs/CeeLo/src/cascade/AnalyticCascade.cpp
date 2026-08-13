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

#include "cascade/AnalyticCascade.h"

#include "cross_sections/CrossSectionData.h"

#include <cmath>

namespace ceelo {

// ===================== CachedEfficiencyProvider ==========================

const EfficiencyResult* CachedEfficiencyProvider::find(double E) const {
    if (cache_.empty()) return nullptr;
    const EfficiencyResult* best = nullptr;
    double best_d = tol_;
    auto hi = cache_.lower_bound(E);
    if (hi != cache_.end()) {
        const double d = std::abs(hi->first - E);
        if (d <= best_d) { best_d = d; best = &hi->second; }
    }
    if (hi != cache_.begin()) {
        auto lo = std::prev(hi);
        const double d = std::abs(lo->first - E);
        if (d <= best_d) { best_d = d; best = &lo->second; }
    }
    return best;
}
double CachedEfficiencyProvider::fep(double E) const {
    const EfficiencyResult* r = find(E);
    return r ? r->full_energy_peak_efficiency : 0.0;
}
double CachedEfficiencyProvider::total(double E) const {
    const EfficiencyResult* r = find(E);
    return r ? r->total_efficiency : 0.0;
}
double CachedEfficiencyProvider::fep_unc(double E) const {
    const EfficiencyResult* r = find(E);
    return r ? r->fep_uncertainty : 0.0;
}
double CachedEfficiencyProvider::total_unc(double E) const {
    const EfficiencyResult* r = find(E);
    return r ? r->total_uncertainty : 0.0;
}
bool CachedEfficiencyProvider::has(double E) const { return find(E) != nullptr; }

// ===================== x-ray line tables (detail) ==========================
// Defined here (not in AnalyticCascade_imp.hpp) so the template header stays
// free of the cross-section data dependency.

namespace detail {

std::vector<XrayLine> k_lines(int Z) {
    std::vector<XrayLine> out;
    const FluorescenceData* fl = CrossSectionData::instance().fluorescence(Z);
    if (!fl) return out;
    for (int i = 0; i < fl->num_lines; ++i)
        out.push_back({fl->line_energy_keV[i],
                       fl->fluorescence_yield * fl->line_probability[i]});
    return out;
}

std::vector<XrayLine> l_lines(int Z, int s) {
    std::vector<XrayLine> out;
    const LFluorescenceData* fl = CrossSectionData::instance().l_fluorescence(Z);
    if (!fl || s < 0 || s > 2) return out;
    const LSubshellFluor& sub = fl->sub[s];
    for (int i = 0; i < sub.num_lines; ++i)
        out.push_back({sub.line_energy_keV[i],
                       sub.fluorescence_yield * sub.line_probability[i]});
    return out;
}

}  // namespace detail

}  // namespace ceelo
