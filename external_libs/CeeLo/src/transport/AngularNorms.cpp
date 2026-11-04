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

#include "transport/AngularNorms.h"

#include "transport/ComptonScatter.h"
#include "transport/RayleighAngular.h"

#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

namespace ceelo {

namespace {

// Log-spaced grid spanning the library's useful photon range with margin.
//
// Grid density is set by the RAYLEIGH norm, not the Compton one. Compton's norm
// is smooth in log E and the cubic interpolation below reaches ~1e-6 with only
// 512 points. The Rayleigh norm inherits small KINKS from the tabulated form
// factor F(x,Z) -- integrating across the table's interpolation breakpoints
// leaves a piecewise-smooth function of energy, and cubic interpolation cannot
// improve on a kink, so only a denser grid helps (512 -> ~1.2e-3,
// 2048 -> well under 1e-4). Tables are built lazily per (process, Z) and cost a
// few tens of milliseconds each, so the denser grid is cheap.
constexpr double kEminKeV = 1.0;
constexpr double kEmaxKeV = 20000.0;
constexpr int kNPoints = 2048;

struct NormTable {
    // log_norm[i] at energy exp(ln_e_min + i*dln)
    std::vector<double> log_norm;
    double ln_e_min;
    double inv_dln;
};

double exact_norm(ScatterProcess p, double energy_keV, int Z)
{
    return (p == ScatterProcess::Compton) ? compton_angular_norm(energy_keV, Z)
                                          : rayleigh_angular_norm(energy_keV, Z);
}

std::unique_ptr<NormTable> build_table(ScatterProcess p, int Z)
{
    auto t = std::make_unique<NormTable>();
    t->ln_e_min = std::log(kEminKeV);
    const double dln = (std::log(kEmaxKeV) - t->ln_e_min) / (kNPoints - 1);
    t->inv_dln = 1.0 / dln;
    t->log_norm.resize(kNPoints);
    for (int i = 0; i < kNPoints; ++i) {
        const double e = std::exp(t->ln_e_min + i * dln);
        const double n = exact_norm(p, e, Z);
        // The norm is strictly positive; guard against an underflow to 0 at the
        // extreme low-energy end so the log is finite.
        t->log_norm[i] = std::log(n > 0.0 ? n : 1e-300);
    }
    return t;
}

using Key = std::pair<int, int>;  // (process, Z)

const NormTable& table_for(ScatterProcess p, int Z)
{
    // Process-wide store: built once, never mutated or erased afterwards, so
    // handing out const references is safe. std::map nodes are address-stable.
    static std::mutex mu;
    static std::map<Key, std::unique_ptr<NormTable>> store;

    const Key key{static_cast<int>(p), Z};

    // Thread-local memo of resolved pointers avoids taking the lock per call.
    thread_local std::map<Key, const NormTable*> local;
    auto lit = local.find(key);
    if (lit != local.end()) return *lit->second;

    const NormTable* ptr = nullptr;
    {
        std::lock_guard<std::mutex> lock(mu);
        auto it = store.find(key);
        if (it == store.end())
            it = store.emplace(key, build_table(p, Z)).first;
        ptr = it->second.get();
    }
    local.emplace(key, ptr);
    return *ptr;
}

} // namespace

double angular_norm_table_min_keV() { return kEminKeV; }
double angular_norm_table_max_keV() { return kEmaxKeV; }

double angular_norm_interp(ScatterProcess process, double energy_keV, int Z)
{
    if (energy_keV < kEminKeV || energy_keV > kEmaxKeV)
        return exact_norm(process, energy_keV, Z);

    const NormTable& t = table_for(process, Z);
    const double u = (std::log(energy_keV) - t.ln_e_min) * t.inv_dln;
    int i = static_cast<int>(u);
    if (i < 0) i = 0;
    if (i > kNPoints - 2) i = kNPoints - 2;
    const double f = u - i;

    // Four-point (cubic Catmull-Rom) interpolation of log(norm) vs log(E).
    // Linear interpolation on this grid is good to ~1.5e-5 for the Compton norm
    // but only ~1.6e-3 for Rayleigh, whose form factor makes the norm fall much
    // more steeply with energy. Cubic is O(h^4) and brings both well under
    // 1e-4 without paying for a denser table (which would be a per-process
    // build cost, since the tables are built lazily at first use).
    const int i0 = (i > 0) ? i - 1 : i;
    const int i3 = (i + 2 <= kNPoints - 1) ? i + 2 : i + 1;
    const double y0 = t.log_norm[i0], y1 = t.log_norm[i];
    const double y2 = t.log_norm[i + 1], y3 = t.log_norm[i3];
    const double a0 = y1;
    const double a1 = 0.5 * (y2 - y0);
    const double a2 = y0 - 2.5 * y1 + 2.0 * y2 - 0.5 * y3;
    const double a3 = 0.5 * (y3 - y0) + 1.5 * (y1 - y2);
    return std::exp(a0 + f * (a1 + f * (a2 + f * a3)));
}

} // namespace ceelo
