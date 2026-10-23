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

#include "cascade/LevelDag.h"
#include "cascade/CascadeTypes.h"
#include "efficiency/EfficiencyCalculator.h"   // EfficiencyResult (for the cache provider)

#include <map>
#include <mutex>
#include <memory>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <initializer_list>

namespace ceelo {

/// Injectable source of detector efficiency at an arbitrary energy, templated on
/// the efficiency scalar type T (see the header comment). InterSpec plugs its
/// DRF here; the default CachedEfficiencyProvider reads a per-energy MC
/// efficiency cache. Implementations should be cheap and thread-safe for reads,
/// and must be deterministic and side-effect free: compute_cascade_analytic
/// memoizes has()/fep()/total() per energy within a call, and may query
/// energies (and in an order) that a given window does not end up using.
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
    /// Also report pure coincidence-sum peaks. A window containing no emitted
    /// line normally comes back found == false with nothing computed; with this
    /// set it also gets `pure_sum_rate` - the absolute per-decay rate of photon
    /// pairs (and triples) summing into it, from exactly the joints, angular
    /// factor and coincident survival a line's summing-in uses - and the rate-
    /// weighted mean energy of those sums. Off by default; windows with a line
    /// are unaffected either way.
    bool report_pure_sums = false;
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
    /// With AnalyticCascadeOptions::report_pure_sums, for a window holding no
    /// emitted line (found == false): the absolute per-decay rate of coincidence
    /// sums into it (counts per parent decay, the efficiencies included), and
    /// their rate-weighted mean energy (0 when there are none).
    T pure_sum_rate = T(0.0);
    double pure_sum_energy_keV = 0.0;
};

/// The plain-double result (the pre-templating API).
using AnalyticPeakResult = AnalyticPeakResultT<double>;

namespace detail {

/// A single relaxation x-ray line: energy and its per-vacancy emission weight
/// (omega * line-branch), so a vacancy of probability p emits this line with
/// probability p*weight.
struct XrayLine { double energy; double weight; };

/// K fluorescence lines for element Z, each weight = omega_K * line_branch.
/// Defined in AnalyticCascade.cpp (needs the cross-section data tables).
std::vector<XrayLine> k_lines(int Z);

/// L fluorescence lines for element Z and subshell s (0=L1,1=L2,2=L3); each
/// weight = omega_Ls * line_branch. Coster-Kronig transfer is a second-order
/// redistribution ignored here (the analytic model uses the direct subshell
/// yield; FR samples the transfer, a <~0.1% effect on summing).
/// Defined in AnalyticCascade.cpp.
std::vector<XrayLine> l_lines(int Z, int s);

/// Absolute per-branch EC shell mass by daughter start level for a legacy
/// LevelScheme which has no explicit WeakOutcomeLaw. LevelDag::feed is
/// normalized conditional on entering the graph, so recover absolute starts
/// with entry_probability and restore the direct-to-ground terminal mass.
struct LegacyEcShellStarts {
    std::vector<double> k, l;
    double total_k = 0.0, total_l = 0.0;
};

/// One emission occurrence of a branch, for summing-in pair enumeration.
/// Kinds: gamma / IC-K / IC-L / EC-K / EC-L / annih-511. `src` keys the source
/// vacancy for exclusivity.
struct Occ {
    int trans;      // level-scheme transition, or -1 for EC / annih
    double energy;
    double factor;  // per-occurrence line/emit factor (applied after the joint)
    int src;        // exclusivity key within the branch
    bool is_ec;
    bool ec_isL;
    bool is_annih;
    int member;     // gamma/annih member index (for angular links); else -1
    int residual = -1;  // categorical unmatched transition, else -1
};

/// Everything about one decay branch that compute_cascade_analytic needs and
/// that does not depend on the efficiencies or the peak windows.
struct PreparedBranch {
    const DecayCascade* dc = nullptr;
    LevelDag dag;
    int Z = 0;
    std::vector<XrayLine> kl;            ///< k_lines(Z)
    std::vector<XrayLine> ll[3];         ///< l_lines(Z, s)
    std::vector<double> start_w;         ///< weak_unknown_start_weights(*dc, dag)
    LegacyEcShellStarts legacy_ec;       ///< empty when the weak outcome law is usable
    std::vector<int> member_residual;    ///< residual_transition_of(*dc, m) per member
    std::vector<Occ> occ;                ///< summing-in emission occurrences
    std::vector<int> ord;                ///< `occ` indices sorted by energy
    std::vector<double> es;              ///< `occ` energies, in `ord` order

    explicit PreparedBranch(const DecayCascade& c);
};

}  // namespace detail

/// The efficiency- and window-independent preparation of a set of decay
/// cascades for compute_cascade_analytic: level-scheme DAGs, x-ray line tables,
/// summing-in occurrence lists, and so on.  Build it once per cascade set and
/// reuse it across evaluations - it is logically immutable (its only mutable
/// state is an internally locked memo of pure functions), so it may be shared
/// between threads.  The lvalue constructor references `cascades`, which must
/// then outlive the model and not be modified; the other two keep the cascades
/// alive themselves.  Copies share the cascades and the memo.  Throws
/// std::invalid_argument for malformed member references.
class AnalyticCascadeModel {
public:
    explicit AnalyticCascadeModel(const std::vector<DecayCascade>& cascades);
    explicit AnalyticCascadeModel(std::vector<DecayCascade>&& cascades);
    explicit AnalyticCascadeModel(std::shared_ptr<const std::vector<DecayCascade>> cascades);

    const std::vector<DecayCascade>& cascades() const { return *cascades_; }
    const std::vector<detail::PreparedBranch>& branches() const { return branches_; }

    /// detail::weak_subset_probability() for branch `b`.  A pure function of
    /// the model, so it is memoized across evaluations (subsets of up to three
    /// transitions; larger ones are computed directly); thread-safe.
    double weak_subset_probability(int b, std::initializer_list<int> subset,
                                   const WeakOutcomeKind* kind = nullptr,
                                   int shell = 0) const;

private:
    void prepare();

    struct SharedMemo {
        explicit SharedMemo(std::size_t nbranch) : mutex(nbranch), wsp(nbranch) {}
        std::vector<std::mutex> mutex;                               ///< per branch
        std::vector<std::unordered_map<std::uint64_t, double>> wsp;  ///< per branch
    };

    std::shared_ptr<const std::vector<DecayCascade>> owned_;
    const std::vector<DecayCascade>* cascades_ = nullptr;
    std::vector<detail::PreparedBranch> branches_;
    std::shared_ptr<SharedMemo> memo_;
};

/// Analytic per-peak cascade-summing factors. `cascades` are the correlated
/// decay branches (from the SandiaDecay adapter or hand-built); `peaks` the
/// photopeaks to report; `provider` supplies eps_FEP/eps_tot at any energy.
/// Defined in AnalyticCascade_imp.hpp (included below); T is deduced from the
/// provider, including through classes derived from EfficiencyProviderT<T>.
/// Prepares a temporary AnalyticCascadeModel; callers evaluating the same
/// cascades repeatedly should build the model once and use the overload below.
template <class T>
std::vector<AnalyticPeakResultT<T>> compute_cascade_analytic(
    const std::vector<DecayCascade>& cascades,
    const std::vector<PeakWindow>& peaks,
    const EfficiencyProviderT<T>& provider,
    const AnalyticCascadeOptions& options = {});

/// The same, for a prepared model - identical results, without re-preparing.
template <class T>
std::vector<AnalyticPeakResultT<T>> compute_cascade_analytic(
    const AnalyticCascadeModel& model,
    const std::vector<PeakWindow>& peaks,
    const EfficiencyProviderT<T>& provider,
    const AnalyticCascadeOptions& options = {});

/// Every energy (keV, ascending) in [lo_keV, hi_keV] at which two emissions of
/// one decay branch can sum: each pair of summing-in occurrences (gamma lines
/// and the expanded vacancy x-ray lines) not from the same transition or
/// vacancy. Energies within `merge_keV` of each other are merged (keeping their
/// mean). No probabilities - a candidate list for pure-sum windows
/// (AnalyticCascadeOptions::report_pure_sums), whose rates say which matter.
inline std::vector<double> coincidence_sum_energies(const AnalyticCascadeModel& model,
                                             double lo_keV, double hi_keV,
                                             double merge_keV = 0.05);

}  // namespace ceelo

#include "cascade/AnalyticCascade_imp.hpp"
