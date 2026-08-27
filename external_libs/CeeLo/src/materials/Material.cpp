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

#include "materials/Material.h"
#include "cross_sections/CrossSectionData.h"

#include <cassert>
#include <array>
#include <atomic>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace ceelo {

uint64_t Material::next_cache_id() {
    // Materials are built once and then queried millions of times, so one
    // relaxed increment per construction costs nothing measurable. Starts at 1
    // so a zero id can never be mistaken for a live one.
    static std::atomic<uint64_t> counter{1};
    return counter.fetch_add(1, std::memory_order_relaxed);
}

Material::Material(const Material& rhs)
    : name_(rhs.name_)
    , density_g_per_cm3_(rhs.density_g_per_cm3_)
    , composition_(rhs.composition_)
    , number_densities_(rhs.number_densities_)
    , cache_id_(next_cache_id()) {}

Material::Material(Material&& rhs) noexcept
    : name_(std::move(rhs.name_))
    , density_g_per_cm3_(rhs.density_g_per_cm3_)
    , composition_(std::move(rhs.composition_))
    , number_densities_(std::move(rhs.number_densities_))
    , cache_id_(next_cache_id()) {}

Material& Material::operator=(const Material& rhs) {
    if (this != &rhs) {
        name_ = rhs.name_;
        density_g_per_cm3_ = rhs.density_g_per_cm3_;
        composition_ = rhs.composition_;
        number_densities_ = rhs.number_densities_;
        cache_id_ = next_cache_id();   // the composition changed; the memo must not follow
    }
    return *this;
}

Material& Material::operator=(Material&& rhs) noexcept {
    if (this != &rhs) {
        name_ = std::move(rhs.name_);
        density_g_per_cm3_ = rhs.density_g_per_cm3_;
        composition_ = std::move(rhs.composition_);
        number_densities_ = std::move(rhs.number_densities_);
        cache_id_ = next_cache_id();
    }
    return *this;
}

Material::Material(std::string name, double density_g_per_cm3,
                   std::vector<MaterialComponent> composition)
    : name_(std::move(name))
    , density_g_per_cm3_(density_g_per_cm3)
    , composition_(std::move(composition))
    , cache_id_(next_cache_id())
{
    if (density_g_per_cm3_ <= 0.0) {
        throw std::invalid_argument("Material density must be positive");
    }
    if (composition_.empty()) {
        throw std::invalid_argument("Material must have at least one component");
    }

    // Validate mass fractions sum to approximately 1.0
    double sum = 0.0;
    for (const auto& c : composition_) {
        if (c.Z < 1 || c.Z > kMaxZ) {
            throw std::invalid_argument("Invalid atomic number Z=" + std::to_string(c.Z));
        }
        if (c.mass_fraction <= 0.0 || c.mass_fraction > 1.0) {
            throw std::invalid_argument("Invalid mass fraction for Z=" + std::to_string(c.Z));
        }
        sum += c.mass_fraction;
    }
    if (std::abs(sum - 1.0) > 0.01) {
        throw std::invalid_argument("Mass fractions sum to " + std::to_string(sum)
                                    + ", should be ~1.0");
    }

    // Pre-compute number densities: n_i = rho * N_A * w_i / A_i
    const auto& xs_data = CrossSectionData::instance();
    number_densities_.reserve(composition_.size());
    for (const auto& c : composition_) {
        double A_i = xs_data.atomic_weight(c.Z);
        double n_i = density_g_per_cm3_ * kAvogadro * c.mass_fraction / A_i;
        number_densities_.push_back(n_i);
    }
}

MacroscopicXS Material::macroscopic_xs(double energy_MeV) const {
    // Per-thread LRU cache keyed on (cache_id_, exact energy).
    // simulate_thread runs at a single fixed primary energy that only changes at
    // Compton vertices, so the same (material, energy) recurs across every transport
    // step / boundary crossing / segment.  macroscopic_xs is a pure function of
    // those inputs (number_densities_ are const after construction;
    // all_cross_sections reads immutable tables, no RNG), so an EXACT-match reuse
    // returns a value bit-for-bit identical to a recompute -- the FEP/total
    // results are unchanged.  The struct is returned BY VALUE (copied out of the
    // entry) because callers mutate their local copy (e.g. PhotonTransport zeroes
    // mu_rs/mu_pp).  Keyed on cache_id_, a per-object identity assigned by every
    // constructor and assignment: keying on `this` guarded only by the density
    // was not enough, because a destroyed Material and a new one built at the
    // same address with the same density but a DIFFERENT composition collided
    // and silently served the first one's cross-sections.  A worker thread's
    // cache lives only for one compute(), the main thread's may persist - which
    // is exactly the lifetime over which addresses get recycled.
    // Bounded LRU (kN entries): secondary energies (Compton/brems
    // continuum) produce many distinct doubles, so an unbounded map would grow
    // without bound -- kN covers the small per-thread working set of materials at
    // the locally-stable current energy.
    struct Entry { uint64_t id; double e; MacroscopicXS xs; };
    static constexpr int kN = 8;
    static thread_local std::array<Entry, kN> cache{};
    static thread_local int n_filled = 0;

    for (int i = 0; i < n_filled; ++i) {
        if (cache[i].id == cache_id_ && cache[i].e == energy_MeV) {
            MacroscopicXS hit = cache[i].xs;          // copy out, never a reference
            if (i != 0) {                              // move-to-front (keep hot)
                for (int j = i; j > 0; --j) cache[j] = cache[j - 1];
                cache[0] = Entry{cache_id_, energy_MeV, hit};
            }
            return hit;
        }
    }

    const auto& xs_data = CrossSectionData::instance();

    MacroscopicXS result{0.0, 0.0, 0.0, 0.0};

    for (size_t i = 0; i < composition_.size(); ++i) {
        int Z = composition_[i].Z;
        auto sigma = xs_data.all_cross_sections(Z, energy_MeV);

        // Σ_type = n_i * σ_type (barn) * 1e-24 (cm²/barn)
        double n = number_densities_[i] * kBarnToCm2;
        result.mu_pe += n * sigma.sigma_pe;
        result.mu_cs += n * sigma.sigma_cs;
        result.mu_rs += n * sigma.sigma_rs;
        result.mu_pp += n * sigma.sigma_pp;
    }

    // Insert at front (LRU eviction of the last slot).
    if (n_filled < kN) ++n_filled;
    for (int j = n_filled - 1; j > 0; --j) cache[j] = cache[j - 1];
    cache[0] = Entry{cache_id_, energy_MeV, result};

    return result;
}

double Material::mu_total(double energy_MeV) const {
    auto xs = macroscopic_xs(energy_MeV);
    return xs.mu_total();
}

double Material::mass_attenuation(double energy_MeV) const {
    return mu_total(energy_MeV) / density_g_per_cm3_;
}

int Material::select_element(double energy_MeV, int interaction_type, double xi) const {
    const auto& xs_data = CrossSectionData::instance();

    if (composition_.size() == 1) {
        return composition_[0].Z;
    }

    // Compute partial macroscopic cross-section contribution from each element
    // for the specified interaction type, then sample proportionally
    constexpr size_t kStackCap = 8;
    std::array<double, kStackCap> cumulative_stack{};
    std::vector<double> cumulative_heap;
    double* cumulative = cumulative_stack.data();
    if (composition_.size() > kStackCap) {
        cumulative_heap.resize(composition_.size());
        cumulative = cumulative_heap.data();
    }
    double running_sum = 0.0;

    for (size_t i = 0; i < composition_.size(); ++i) {
        int Z = composition_[i].Z;
        auto sigma = xs_data.all_cross_sections(Z, energy_MeV);

        double contribution = 0.0;
        switch (interaction_type) {
            case 0: contribution = number_densities_[i] * sigma.sigma_pe; break;
            case 1: contribution = number_densities_[i] * sigma.sigma_cs; break;
            case 2: contribution = number_densities_[i] * sigma.sigma_rs; break;
            case 3: contribution = number_densities_[i] * sigma.sigma_pp; break;
            default: assert(false && "Invalid interaction type");
        }
        running_sum += contribution;
        cumulative[i] = running_sum;
    }

    if (running_sum <= 0.0) {
        // Fallback: return first element
        return composition_[0].Z;
    }

    double threshold = xi * running_sum;
    for (size_t i = 0; i < composition_.size(); ++i) {
        if (threshold <= cumulative[i]) {
            return composition_[i].Z;
        }
    }

    // Rounding: return last element
    return composition_.back().Z;
}

double Material::number_density(size_t component_index) const {
    assert(component_index < number_densities_.size());
    return number_densities_[component_index];
}

} // namespace ceelo
