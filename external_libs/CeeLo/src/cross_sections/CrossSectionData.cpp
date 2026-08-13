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

#include "cross_sections/CrossSectionData.h"
#include "cross_sections/element_data.h"
#include "cross_sections/photon_epics_data.h"
#include "cross_sections/relaxation_epics_data.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace ceelo {

namespace {
// 10^x via exp2.  std::pow(10.0, x) routes through the general (slower) libm pow
// path; exp2(x*log2(10)) hits a dedicated fast routine.  NOT bit-identical to
// std::pow(10.0, x) -- the rounded log2(10) constant shifts the result by ~1e-16
// relative -- so this is a STATISTICALLY result-preserving substitution (validated
// by the test suite), used only on the cross-section interpolation leaf, never
// where bit reproducibility is required (e.g. the init-time Rayleigh CDF build
// keeps std::pow so its tables are unperturbed).
inline double pow10(double x) {
    constexpr double kLog2_10 = 3.321928094887362347870319429489; // log2(10)
    return std::exp2(x * kLog2_10);
}
} // namespace

const CrossSectionData& CrossSectionData::instance() {
    static const CrossSectionData inst;
    return inst;
}

const ElementData& CrossSectionData::element(int Z) const {
    assert(Z >= 1 && Z <= kMaxZ);
    return g_element_data[Z - 1];
}

double CrossSectionData::sample_rayleigh_x(int Z, double x_max,
                                           std::mt19937_64& rng) const {
    assert(Z >= 1 && Z <= kMaxZ);
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Below the grid (very low energy / pure forward): F≈const ⇒ G ∝ x^2,
    // so x = x_max * sqrt(U).
    const double lx = std::log10(x_max);
    if (lx <= kRayleighSamplingLogXMin) {
        return x_max * std::sqrt(uniform(rng));
    }

    const auto& cdf = g_rayleigh_sampling_cdf_q[Z - 1];
    constexpr int n = kRayleighSamplingNodes;
    constexpr double log_span =
        kRayleighSamplingLogXMax - kRayleighSamplingLogXMin;

    // Normalized G(x_max), retained in uint16 code units, and the grid index
    // just above x_max.
    double Gmax;
    int imax;
    if (lx >= kRayleighSamplingLogXMax) {
        Gmax = cdf[n - 1];
        imax = n - 1;
    } else {
        const double f = (lx - kRayleighSamplingLogXMin) * (n - 1) / log_span;
        int j = static_cast<int>(f);
        if (j < 0) j = 0;
        if (j >= n - 1) j = n - 2;
        const double t = f - j;
        Gmax = cdf[j] + t * (cdf[j + 1] - cdf[j]);
        imax = j + 1;
    }

    double target = uniform(rng) * Gmax;

    // The first CDF code represents the analytic [0, x_min] contribution.
    // Preserve its F≈constant, G∝x² sampling law rather than extrapolating the
    // first tabulated segment below x_min.
    if (cdf[0] > 0 && target <= cdf[0]) {
        return g_rayleigh_sampling_x[0] * std::sqrt(target / cdf[0]);
    }

    // Invert: largest index lo (< imax) with cdf[lo] <= target.
    int lo = 0, hi = imax;
    while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (cdf[mid] <= target) lo = mid; else hi = mid;
    }
    const double seg = static_cast<double>(cdf[lo + 1]) - cdf[lo];
    const double t = (seg > 0.0) ? (target - cdf[lo]) / seg : 0.0;
    const double x = g_rayleigh_sampling_x[lo]
        + t * (g_rayleigh_sampling_x[lo + 1] - g_rayleigh_sampling_x[lo]);
    return (x < x_max) ? x : x_max;
}

// --- Interpolation helpers ---

double CrossSectionData::interpolate_log_log(
    const PhotonProcessCurve& curve, double log_x)
{
    const uint16_t n = curve.size;
    if (n == 0) return 0.0;

    const auto grid = [&](size_t index) {
        const uint16_t pool_index =
            g_photon_process_grid_index[curve.data_offset + index];
        return static_cast<double>(g_photon_energy_pool[pool_index]);
    };

    const auto code_at = [&](size_t index) {
        return g_photon_process_log_value_q[curve.data_offset + index];
    };
    const auto decode_log = [&](size_t index) {
        const uint16_t code = code_at(index);
        return code == 0 ? -30.0
                         : static_cast<double>(curve.value_offset)
                           + static_cast<double>(code - 1) * curve.value_scale;
    };

    // Clamp to grid bounds
    if (log_x <= grid(0)) return code_at(0) ? pow10(decode_log(0)) : 0.0;
    if (log_x >= grid(n - 1))
        return code_at(n - 1) ? pow10(decode_log(n - 1)) : 0.0;

    size_t i_lo = 0, i_hi = n - 1;
    while (i_hi - i_lo > 1) {
        const size_t mid = (i_lo + i_hi) >> 1;
        if (grid(mid) <= log_x) i_lo = mid; else i_hi = mid;
    }

    // Log-log linear interpolation
    const double low_grid = grid(i_lo);
    double t = (log_x - low_grid) / (grid(i_hi) - low_grid);
    if (code_at(i_lo) == 0 || code_at(i_hi) == 0) {
        const double low_value = code_at(i_lo) ? pow10(decode_log(i_lo)) : 0.0;
        const double high_value = code_at(i_hi) ? pow10(decode_log(i_hi)) : 0.0;
        return low_value + t * (high_value - low_value);
    }
    const double low = decode_log(i_lo);
    const double log_result = low + t * (decode_log(i_hi) - low);
    return pow10(log_result);
}

double CrossSectionData::interpolate_rayleigh(
    int Z, const PhotonRayleighValues& values, double log_x)
{
    const uint16_t group = static_cast<uint16_t>(
        (Z - 1) / kRayleighXsElementsPerGroup
    );
    const uint16_t grid_begin = g_rayleigh_group_grid_offset[group];
    const uint16_t grid_end = g_rayleigh_group_grid_offset[group + 1];
    const uint16_t n = static_cast<uint16_t>(grid_end - grid_begin);
    const uint16_t lane = static_cast<uint16_t>(
        (Z - 1) % kRayleighXsElementsPerGroup
    );
    const uint16_t value_begin = static_cast<uint16_t>(
        g_rayleigh_group_value_offset[group] + lane * n
    );
    const float* grid = g_rayleigh_log_energy + grid_begin;
    const uint16_t* packed = g_rayleigh_log_value_q + value_begin;
    const auto decode_log = [&](size_t index) {
        const uint16_t code = packed[index];
        return code == 0 ? -30.0
                         : static_cast<double>(values.value_offset)
                           + static_cast<double>(code - 1) * values.value_scale;
    };

    if (log_x <= grid[0]) return packed[0] ? pow10(decode_log(0)) : 0.0;
    if (log_x >= grid[n - 1])
        return packed[n - 1] ? pow10(decode_log(n - 1)) : 0.0;

    const float log_x_f = static_cast<float>(log_x);
    const float* upper = std::upper_bound(grid, grid + n, log_x_f);
    size_t i_hi = static_cast<size_t>(upper - grid);
    if (i_hi == 0) i_hi = 1;
    const size_t i_lo = i_hi - 1;
    const double t = (log_x - static_cast<double>(grid[i_lo]))
                   / (static_cast<double>(grid[i_hi]) - grid[i_lo]);
    const double low = decode_log(i_lo);
    return pow10(low + t * (decode_log(i_hi) - low));
}

// Linear-in-value interpolation on a log-spaced grid. Used for S(x,Z), which
// rises smoothly from zero to Z.
double CrossSectionData::interpolate_log_lin(
    const PhotonAngularCurve& curve, double log_x)
{
    const uint16_t n = curve.size;
    if (n == 0) return 0.0;

    const float* log_grid = g_photon_angular_log_x + curve.data_offset;
    const uint16_t* packed_log_values =
        g_photon_angular_log_value_q + curve.data_offset;

    const auto decode = [=](size_t index) {
        const uint16_t code = packed_log_values[index];
        if (code == 0) return 0.0;
        const double stored = static_cast<double>(curve.value_offset)
                            + static_cast<double>(code - 1) * curve.value_scale;
        return pow10(stored);
    };

    // Clamp to grid bounds
    if (log_x <= log_grid[0]) return decode(0);
    if (log_x >= log_grid[n - 1]) return decode(n - 1);

    const float log_x_f = static_cast<float>(log_x);
    auto it = std::upper_bound(log_grid, log_grid + n, log_x_f);
    size_t i_hi = static_cast<size_t>(it - log_grid);
    if (i_hi == 0) i_hi = 1;
    size_t i_lo = i_hi - 1;

    double t = (log_x - static_cast<double>(log_grid[i_lo]))
             / (static_cast<double>(log_grid[i_hi]) - static_cast<double>(log_grid[i_lo]));
    const double low = decode(i_lo);
    return low + t * (decode(i_hi) - low);
}

// --- Individual cross-section queries ---

double CrossSectionData::sigma_photoelectric(int Z, double energy_MeV) const {
    assert(Z >= 1 && Z <= kMaxZ);
    const auto& curve = g_photon_epics_data[Z - 1].photoelectric;
    double log_E = std::log10(energy_MeV);
    return interpolate_log_log(curve, log_E);
}

double CrossSectionData::sigma_K_photoelectric(int Z, double energy_MeV) const {
    assert(Z >= 1 && Z <= kMaxZ);
    const auto& curve = g_photon_epics_data[Z - 1].k_photoelectric;
    if (curve.size == 0) return 0.0;
    const auto first_grid_index = g_photon_process_grid_index[curve.data_offset];
    if (energy_MeV < pow10(g_photon_energy_pool[first_grid_index])) return 0.0;
    double log_E = std::log10(energy_MeV);
    return interpolate_log_log(curve, log_E);
}

double CrossSectionData::sigma_compton(int Z, double energy_MeV) const {
    assert(Z >= 1 && Z <= kMaxZ);
    const auto& curve = g_photon_epics_data[Z - 1].compton;
    double log_E = std::log10(energy_MeV);
    return interpolate_log_log(curve, log_E);
}

double CrossSectionData::sigma_rayleigh(int Z, double energy_MeV) const {
    assert(Z >= 1 && Z <= kMaxZ);
    const auto& curve = g_photon_epics_data[Z - 1].rayleigh;
    double log_E = std::log10(energy_MeV);
    return interpolate_rayleigh(Z, curve, log_E);
}

double CrossSectionData::sigma_pair_production(int Z, double energy_MeV) const {
    assert(Z >= 1 && Z <= kMaxZ);
    if (energy_MeV < 1.022e-0) return 0.0; // Below pair production threshold
    const auto& curve = g_photon_epics_data[Z - 1].pair_production;
    double log_E = std::log10(energy_MeV);
    return interpolate_log_log(curve, log_E);
}

CrossSectionData::PartialCrossSections
CrossSectionData::all_cross_sections(int Z, double energy_MeV) const {
    assert(Z >= 1 && Z <= kMaxZ);

    // Per-thread (Z, energy) LRU cache.  all_cross_sections is a pure function of
    // (Z, energy) over the immutable compiled-in tables, so an EXACT-match reuse is
    // bit-for-bit identical to a recompute.  select_element() and the
    // macroscopic_xs() miss path each query every element of the current material
    // at the same (locally-stable) primary energy, so the small per-energy working
    // set fits; secondaries at many distinct energies are handled by LRU eviction.
    struct CacheEntry { int Z; double e; PartialCrossSections xs; };
    static constexpr int kCacheN = 16;
    static thread_local std::array<CacheEntry, kCacheN> cache{};
    static thread_local int n_filled = 0;
    for (int i = 0; i < n_filled; ++i) {
        if (cache[i].Z == Z && cache[i].e == energy_MeV) {
            PartialCrossSections hit = cache[i].xs;      // copy out, never a reference
            if (i != 0) {                                 // move-to-front (keep hot)
                for (int j = i; j > 0; --j) cache[j] = cache[j - 1];
                cache[0] = CacheEntry{Z, energy_MeV, hit};
            }
            return hit;
        }
    }

    // Compute each process against its own direct EPDL grid.
    const PartialCrossSections result = [&]() -> PartialCrossSections {
        double log_E = std::log10(energy_MeV);
        const auto& data = g_photon_epics_data[Z - 1];
        return PartialCrossSections{
            interpolate_log_log(data.photoelectric, log_E),
            interpolate_log_log(data.compton, log_E),
            interpolate_rayleigh(Z, data.rayleigh, log_E),
            energy_MeV >= 1.022
                ? interpolate_log_log(data.pair_production, log_E)
                : 0.0
        };
    }();

    // Insert at front (LRU eviction of the last slot).
    if (n_filled < kCacheN) ++n_filled;
    for (int j = n_filled - 1; j > 0; --j) cache[j] = cache[j - 1];
    cache[0] = CacheEntry{Z, energy_MeV, result};
    return result;
}

double CrossSectionData::scattering_function_S(int Z, double x) const {
    assert(Z >= 1 && Z <= kMaxZ);
    const auto& curve = g_photon_epics_data[Z - 1].scattering_function;
    if (curve.size == 0) return static_cast<double>(Z); // fallback: free-electron limit
    double log_x = std::log10(x);
    return interpolate_log_lin(curve, log_x);
}

const FluorescenceData* CrossSectionData::fluorescence(int Z) const {
    // Relaxation has a slightly wider domain than photon/electron transport:
    // radioactive-decay daughters can reach Z=99 (e.g. Np=93 from Am-241).
    if (Z < 1 || Z > kRelaxationMaxZ) return nullptr;
    const auto& data = g_epics_k_relaxation[Z];
    if (data.fluorescence_yield <= 0.0f || data.num_lines == 0) return nullptr;
    return &data;
}

const LFluorescenceData* CrossSectionData::l_fluorescence(int Z) const {
    if (Z < 1 || Z > kRelaxationMaxZ) return nullptr;
    const auto& data = g_epics_l_relaxation[Z];
    for (int shell = 0; shell < 3; ++shell) {
        if (data.sub[shell].fluorescence_yield > 0.0f && data.sub[shell].num_lines > 0)
            return &data;
    }
    return nullptr;
}

double CrossSectionData::sb_chi(int Z, double T_keV, double kappa) const {
    assert(Z >= 1 && Z <= kMaxZ);
    const auto& e = g_element_data[Z - 1];
    if (!e.sb_chi_quantized) return 1.0;  // fallback: flat chi => 1/k spectrum

    // Clamp kappa to table range
    if (kappa < kSB_kappa[0]) kappa = kSB_kappa[0];
    if (kappa > kSB_kappa[kSB_n_kappa - 1]) kappa = kSB_kappa[kSB_n_kappa - 1];

    double log_T = std::log10(T_keV);

    // Find energy bracket
    int iE = 0;
    for (int i = 1; i < kSB_n_energy; ++i) {
        if (kSB_log_E_keV[i] > log_T) break;
        iE = i;
    }
    if (iE >= kSB_n_energy - 1) iE = kSB_n_energy - 2;

    double tE = (log_T - kSB_log_E_keV[iE]) /
                (kSB_log_E_keV[iE + 1] - kSB_log_E_keV[iE]);
    tE = std::max(0.0, std::min(1.0, tE));

    // Find kappa bracket
    int iK = 0;
    for (int i = 1; i < kSB_n_kappa; ++i) {
        if (kSB_kappa[i] > kappa) break;
        iK = i;
    }
    if (iK >= kSB_n_kappa - 1) iK = kSB_n_kappa - 2;

    double tK = (kappa - kSB_kappa[iK]) / (kSB_kappa[iK + 1] - kSB_kappa[iK]);
    tK = std::max(0.0, std::min(1.0, tK));

    // Bilinear interpolation
    int stride = kSB_n_kappa;
    const double scale = kSB_chi_scale[Z - 1];
    double c00 = e.sb_chi_quantized[iE * stride + iK] * scale;
    double c01 = e.sb_chi_quantized[iE * stride + iK + 1] * scale;
    double c10 = e.sb_chi_quantized[(iE + 1) * stride + iK] * scale;
    double c11 = e.sb_chi_quantized[(iE + 1) * stride + iK + 1] * scale;

    double c0 = c00 + tK * (c01 - c00);
    double c1 = c10 + tK * (c11 - c10);
    return std::max(0.0, c0 + tE * (c1 - c0));
}

// --- Hoisted compound Seltzer-Berger chi (bit-identical split of sb_chi) ---
// The grids kSB_log_E_keV / kSB_kappa are global (Z-independent), so the energy
// and kappa brackets are computed once and shared across a material's elements.
// Each piece replicates sb_chi()'s exact arithmetic; for matching brackets,
// sb_chi_bracketed(Z, eb, kb) == sb_chi(Z, T_keV, kappa) to the last bit.

CrossSectionData::SBChiEBracket
CrossSectionData::sb_chi_energy_bracket(double T_keV) const {
    double log_T = std::log10(T_keV);
    int iE = 0;
    for (int i = 1; i < kSB_n_energy; ++i) {
        if (kSB_log_E_keV[i] > log_T) break;
        iE = i;
    }
    if (iE >= kSB_n_energy - 1) iE = kSB_n_energy - 2;
    double tE = (log_T - kSB_log_E_keV[iE]) /
                (kSB_log_E_keV[iE + 1] - kSB_log_E_keV[iE]);
    tE = std::max(0.0, std::min(1.0, tE));
    return {iE, tE};
}

CrossSectionData::SBChiKBracket
CrossSectionData::sb_chi_kappa_bracket(double kappa) const {
    if (kappa < kSB_kappa[0]) kappa = kSB_kappa[0];
    if (kappa > kSB_kappa[kSB_n_kappa - 1]) kappa = kSB_kappa[kSB_n_kappa - 1];
    int iK = 0;
    for (int i = 1; i < kSB_n_kappa; ++i) {
        if (kSB_kappa[i] > kappa) break;
        iK = i;
    }
    if (iK >= kSB_n_kappa - 1) iK = kSB_n_kappa - 2;
    double tK = (kappa - kSB_kappa[iK]) / (kSB_kappa[iK + 1] - kSB_kappa[iK]);
    tK = std::max(0.0, std::min(1.0, tK));
    return {iK, tK};
}

double CrossSectionData::sb_chi_bracketed(int Z, const SBChiEBracket& eb,
                                          const SBChiKBracket& kb) const {
    assert(Z >= 1 && Z <= kMaxZ);
    const auto& e = g_element_data[Z - 1];
    if (!e.sb_chi_quantized) return 1.0;  // fallback: flat chi => 1/k spectrum (matches sb_chi)

    const int stride = kSB_n_kappa;
    const double scale = kSB_chi_scale[Z - 1];
    double c00 = e.sb_chi_quantized[eb.iE * stride + kb.iK] * scale;
    double c01 = e.sb_chi_quantized[eb.iE * stride + kb.iK + 1] * scale;
    double c10 = e.sb_chi_quantized[(eb.iE + 1) * stride + kb.iK] * scale;
    double c11 = e.sb_chi_quantized[(eb.iE + 1) * stride + kb.iK + 1] * scale;

    double c0 = c00 + kb.tK * (c01 - c00);
    double c1 = c10 + kb.tK * (c11 - c10);
    return std::max(0.0, c0 + eb.tE * (c1 - c0));
}

double CrossSectionData::atomic_weight(int Z) const {
    assert(Z >= 1 && Z <= kMaxZ);
    return g_atomic_weights[Z - 1];
}

} // namespace ceelo
