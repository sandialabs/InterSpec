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

// Fully-analytic true-coincidence (cascade) summing.
//
// compute_cascade_analytic() returns per-peak summing factors from a nuclide's
// DecayCascade(s) and an injectable efficiency provider (eps_FEP(E), eps_tot(E)),
// WITHOUT running any Monte Carlo. It is the path InterSpec calls when
// integrating summing over a source volume with the parameterized detector
// response (DRF): for each sampled position the DRF supplies eps_FEP/eps_tot and
// this returns the per-peak factors, so the volume integral just accumulates.
//
// It implements the exact level-scheme joint-probability form (response function
// spec Eq. 10'): summing-OUT via a per-transition survival-factor DP over the
// daughter level scheme (exact joints -> mutually-exclusive transitions are not
// summed), including EC-feed and internal-conversion K/L vacancy x-rays as
// coincident members; and summing-IN over gamma-gamma, gamma+x-ray and
// x-ray+x-ray pairs (plus optional triples) with the exact joint emission
// probabilities. The FullRealization MC (EfficiencyCalculator::compute_cascade,
// CascadeMethod::FullRealization) is the ground-truth reference this matches.
//
// The math (survival DP, no-triple-subtraction summing-in weight, same-vacancy
// exclusivity, W(theta) g-factor) is documented at the call sites in
// AnalyticCascade_imp.hpp.
//
// GENERIC SCALAR TYPE: the provider and the computation are templated on the
// efficiency scalar type T so an automatic-differentiation dual number (e.g.
// ceres::Jet) can flow derivatives of the efficiencies through the summing
// factors. Emission probabilities, branching ratios and energies stay double —
// only efficiency-valued quantities are T. Requirements on T:
//   - constructible from double: T(0.0);
//   - arithmetic with itself and with double (+, -, *, /, +=, *=);
//   - comparison against double (>, <, >=, <=) on its scalar part;
//   - sqrt(T) findable via ADL or std:: (after `using std::sqrt;`).
// ceres::Jet satisfies all of these natively; T = double is the default path and
// `EfficiencyProvider` / `AnalyticPeakResult` / the double overload keep every
// pre-templating call site source-compatible.

#include "cascade/CascadeTypes.h"
#include "efficiency/EfficiencyCalculator.h"   // EfficiencyResult (for the cache provider)

#include <map>
#include <vector>

namespace ceelo {

/// Injectable source of detector efficiency at an arbitrary energy, templated on
/// the efficiency scalar type T (see the header comment). InterSpec plugs its
/// DRF here; the default CachedEfficiencyProvider reads a per-energy MC
/// efficiency cache. Implementations should be cheap and thread-safe for reads.
template <class T>
struct EfficiencyProviderT {
    virtual ~EfficiencyProviderT() = default;
    /// Full-energy-peak efficiency at E (keV).
    virtual T fep(double energy_keV) const = 0;
    /// Total (any-deposit) efficiency at E (keV).
    virtual T total(double energy_keV) const = 0;
    /// 1-sigma uncertainties (default 0 => no uncertainty propagated).
    virtual T fep_unc(double /*energy_keV*/) const { return T(0.0); }
    virtual T total_unc(double /*energy_keV*/) const { return T(0.0); }
    /// True if this provider can supply an efficiency for E within its match
    /// tolerance. A false return is reported through AnalyticPeakResult::
    /// unmatched_energies rather than silently treated as zero efficiency.
    virtual bool has(double /*energy_keV*/) const { return true; }
};

/// The plain-double provider interface (the pre-templating API).
using EfficiencyProvider = EfficiencyProviderT<double>;

/// Default provider over the existing std::map<double,EfficiencyResult> MC cache
/// (keyed by energy in keV). Nearest-neighbour lookup within `match_tol_keV`
/// (mirrors find_cached_efficiency in EfficiencyCalculator.cpp). Misses report
/// via has()==false and are surfaced in unmatched_energies.
class CachedEfficiencyProvider : public EfficiencyProvider {
public:
    CachedEfficiencyProvider(const std::map<double, EfficiencyResult>& cache,
                             double match_tol_keV = kCascadeEnergyMatchTolKeV)
        : cache_(cache), tol_(match_tol_keV) {}

    double fep(double energy_keV) const override;
    double total(double energy_keV) const override;
    double fep_unc(double energy_keV) const override;
    double total_unc(double energy_keV) const override;
    bool   has(double energy_keV) const override;

private:
    const EfficiencyResult* find(double energy_keV) const;
    const std::map<double, EfficiencyResult>& cache_;
    double tol_;
};

/// Options for compute_cascade_analytic().
struct AnalyticCascadeOptions {
    /// Energy-match tolerance the peak locator uses to find the peak's gamma in
    /// the cascades (not the provider's own lookup tolerance).
    double match_tol_keV = kCascadeEnergyMatchTolKeV;
    /// Enumerate triple-fed summing-in (three emissions into one window). Cheap;
    /// eps^3-small at far geometry, ~1% at contact.
    bool enumerate_triples = true;
    // NOTE: FullRealization's Kalpha->L secondary vacancy (a K vacancy relaxing
    // via Kalpha leaves an L vacancy that also radiates) is deliberately NOT
    // modeled here: for summing-OUT it is idempotent (the K x-ray already removed
    // the event, so the extra L deposit changes nothing) and for summing-IN it is
    // eps^2-small. It affects only the sum CONTINUUM, not the peak efficiencies —
    // Ba-133 (Cs Kalpha->L) already agrees to <1% without it.
    /// +/- half-width of the summing energy window (keV).
    double window_tol_keV = 1.5;
    /// Apply the gamma-gamma angular correlation W(theta) to correlated SUMMING-IN
    /// pairs, using the collinear limit g = W(0) = 1 + a2 + a4 (a2/a4 from the
    /// cascade links). This is geometry-independent and correct for coincidence
    /// FEP detection: both photons must point into the detector to both full-
    /// deposit, so their mutual angle is small and W(theta_ab) ~ W(0) (the FEP
    /// acceptance is on-axis-peaked; a uniform-cap average under-weights the small
    /// angles). Measured to recover Co-57 122+14->136 to ~0.2% (vs ~1.9% with
    /// g=1). Summing-OUT partners get g=1 (they need only ANY deposit, so the
    /// mutual angle is loosely constrained and the effect is negligible — Co-60
    /// summing-out agrees to <0.5% with g=1). No detector geometry is required.
    bool apply_angular_correlation = true;
};

/// Per-peak analytic summing result. Mirrors PeakCascadeResult where they
/// overlap, plus the decomposed C_out / C_in factors. Efficiency-valued fields
/// carry the provider's scalar type T; energies and the found flag stay double.
template <class T>
struct AnalyticPeakResultT {
    double energy_keV = 0.0;
    bool   found = false;                  ///< a matching gamma member was located
    T c_out = T(1.0);                      ///< summing-out factor Pi(1 - ...)
    T c_in = T(0.0);                       ///< summing-in term (absolute gain / eps_FEP(peak))
    T c_net = T(1.0);                      ///< c_out + c_in
    T eff_no_summing = T(0.0);             ///< eps_FEP(peak)
    T eff_with_summing = T(0.0);           ///< eps_FEP(peak) * c_net
    T eff_no_summing_unc = T(0.0);
    T eff_with_summing_unc = T(0.0);
    T summing_factor = T(1.0);             ///< eff_with_summing / eff_no_summing (= c_net)
    T summing_factor_unc = T(0.0);
    /// Emission energies requested from the provider that it could not supply
    /// (has()==false). Non-empty => an incomplete provider, treat as an error.
    std::vector<double> unmatched_energies;
    /// True only when every branch material to this peak has an accepted level
    /// graph.  False means the numeric result is retained for compatibility but
    /// is incomplete: an invalid primary uses the legacy pairwise summing-out
    /// approximation, or a plausible invalid-branch photon pair is omitted from
    /// analytic summing-in.  Use FullRealization for a numeric cross-check, while
    /// noting its invalid-branch sampler is itself only an approximate model.
    bool summing_model_complete = true;
};

/// The plain-double result (the pre-templating API).
using AnalyticPeakResult = AnalyticPeakResultT<double>;

/// Analytic per-peak cascade-summing factors. `cascades` are the correlated
/// decay branches (from the SandiaDecay adapter or hand-built); `peaks` the
/// photopeaks to report; `provider` supplies eps_FEP/eps_tot at any energy.
/// Defined in AnalyticCascade_imp.hpp (included below); T is deduced from the
/// provider, including through classes derived from EfficiencyProviderT<T>.
template <class T>
std::vector<AnalyticPeakResultT<T>> compute_cascade_analytic(
    const std::vector<DecayCascade>& cascades,
    const std::vector<PeakWindow>& peaks,
    const EfficiencyProviderT<T>& provider,
    const AnalyticCascadeOptions& options = {});

}  // namespace ceelo

#include "cascade/AnalyticCascade_imp.hpp"
