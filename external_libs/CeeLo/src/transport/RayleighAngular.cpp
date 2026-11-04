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

#include "transport/RayleighAngular.h"

#include "transport/GaussLegendre.h"
#include "cross_sections/CrossSectionData.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

namespace ceelo {

namespace {
// h*c in keV*Angstrom — same value ComptonScatter.cpp uses for the momentum
// transfer x = (E/hc) * sin(theta/2).
constexpr double kHcKeVAngstrom = 12.398;
} // namespace

double rayleigh_angular_pdf_unnorm(double energy_keV, int Z, double mu)
{
    mu = std::max(-1.0, std::min(1.0, mu));
    const double sin_half = std::sqrt(std::max(0.0, (1.0 - mu) * 0.5));
    const double x = (energy_keV / kHcKeVAngstrom) * sin_half;
    const double F = CrossSectionData::instance().form_factor_F(Z, x);
    return (1.0 + mu * mu) * F * F;
}

double rayleigh_angular_norm(double energy_keV, int Z)
{
    thread_local std::map<std::pair<int, double>, double> cache;
    const auto key = std::make_pair(Z, energy_keV);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    // 16 composite panels of 32-point GL in s = sin(theta/2), matching
    // compton_angular_norm(). F^2 is strongly forward-peaked, so the panels
    // concentrate resolution where it matters after the substitution.
    const double norm = gl::integrate32_composite(
        0.0, 1.0, 16, [&](double s) {
            return 4.0 * s *
                   rayleigh_angular_pdf_unnorm(energy_keV, Z, 1.0 - 2.0 * s * s);
        });

    if (cache.size() > 4096) cache.clear();  // unbounded-growth guard
    cache.emplace(key, norm);
    return norm;
}

} // namespace ceelo
