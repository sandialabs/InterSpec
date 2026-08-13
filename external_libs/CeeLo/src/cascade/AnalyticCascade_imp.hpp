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

// Template implementation of compute_cascade_analytic (see AnalyticCascade.h for
// the API contract and the generic-scalar-type requirements on T). Included at
// the bottom of AnalyticCascade.h — do not include directly.
//
// Everything probability-shaped (branching, joints, emission fractions, x-ray
// line weights, angular g) is double; everything efficiency-shaped is T so a
// dual-number T carries d(eps)/d(param) through to the summing factors. The
// non-template x-ray line tables (detail::k_lines / detail::l_lines) are defined
// in AnalyticCascade.cpp so this header stays free of the cross-section data
// dependency.

#include "cascade/CascadeTypes.h"
#include "cascade/LevelDag.h"

#include <set>
#include <cmath>
#include <vector>
#include <algorithm>

namespace ceelo {

namespace detail {

// EC leaves one vacancy; when the L-subshell is unresolved FullRealization
// samples it by occupancy (L1:L2:L3 = 2:2:4). Match that mixture.
inline constexpr double kLOcc[3] = {0.25, 0.25, 0.5};
inline constexpr double kAnnih511 = 510.998950;  // annihilation photon energy (keV)

/// Angular-correlation g-factor for a coincident SUMMING-IN pair. Both photons
/// must be FEP-detected, i.e. both must point into the detector, so their mutual
/// angle theta_ab is small and W(theta_ab) ~ W(0) = 1 + a2 + a4 — the collinear
/// limit, independent of geometry. (A uniform-cap average under-weights the small
/// angles that coincident FEP detection actually favors: the FEP acceptance is
/// on-axis-peaked, so the two directions cluster near the axis; measured on Co-57
/// 122+14->136 the effect tracks W(0), not the cap average.) Used only for γ-γ
/// sum-in pairs; summing-OUT (coincident needs only ANY deposit, ε_tot, far less
/// on-axis-weighted) and x-ray partners (isotropic) get g = 1.
inline double angular_g_w0(double a2, double a4) { return 1.0 + a2 + a4; }

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

/// Efficiency lookups that record misses (energies the provider cannot supply)
/// for the caller to surface. Missing x-ray/low energies contribute 0.
template <class T>
struct EffLookup {
    const EfficiencyProviderT<T>& p;
    std::set<double>& unmatched;
    T tot(double E) const {
        if (p.has(E)) return p.total(E);
        unmatched.insert(E);
        return T(0.0);
    }
    T fep(double E) const {
        if (p.has(E)) return p.fep(E);
        unmatched.insert(E);
        return T(0.0);
    }
};

// ---- Per-branch precompute -------------------------------------------------

struct AnalyticBranch {
    const DecayCascade* dc = nullptr;
    LevelDag dag;
    int Z = 0;
    explicit AnalyticBranch(const DecayCascade& c)
        : dc(&c), dag(c.level_scheme),
          Z(c.daughter_Z ? c.daughter_Z : c.level_scheme.daughter_Z) {}
};

inline std::vector<double> weak_unknown_start_weights(const DecayCascade& dc,
                                                      const LevelDag& dag) {
    std::vector<double> w(dag.n, 0.0), known(dag.n, 0.0);
    double unknown_mass = 0.0;
    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
        if (o.fed_level >= 0 && o.fed_level < dag.n)
            known[o.fed_level] += o.selected_mass;
        else
            unknown_mass += o.selected_mass;
    }
    if (!(unknown_mass > 0.0)) return w;
    double residual = 0.0;
    for (int l = 0; l < dag.n; ++l) {
        w[l] = std::max(0.0, dc.level_scheme.levels[l].feeding - known[l]);
        residual += w[l];
    }
    const double scale = residual > unknown_mass ? unknown_mass / residual : 1.0;
    for (double& x : w) x *= scale / unknown_mass;
    return w;
}

inline double weak_outcome_reach(const DecayCascade& dc, const LevelDag& dag,
                                 const WeakOutcome& o,
                                 const std::vector<int>& subset) {
    if (subset.empty()) return 1.0;
    if (o.fed_level >= 0 && o.fed_level < dag.n)
        return dag.subset_from_level(o.fed_level, subset);
    const std::vector<double> w = weak_unknown_start_weights(dc, dag);
    double p = 0.0;
    for (int l = 0; l < dag.n; ++l)
        if (w[l] > 0.0) p += w[l] * dag.subset_from_level(l, subset);
    return p;
}

inline double weak_subset_probability(const DecayCascade& dc,
                                      const LevelDag& dag,
                                      const std::vector<int>& subset,
                                      const WeakOutcomeKind* only_kind = nullptr,
                                      int ec_shell = 0) {
    if (!dag.valid) return 0.0;
    if (!dc.weak_outcome_law.usable())
        return dc.level_scheme.entry_probability * dag.n_subset_joint(subset);
    double p = 0.0;
    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
        if (only_kind && o.kind != *only_kind) continue;
        double shell = 1.0;
        if (ec_shell == 1) shell = o.ec_K;
        else if (ec_shell == 2) shell = o.ec_L;
        const double reach = weak_outcome_reach(dc, dag, o, subset);
        p += o.selected_mass * shell * reach;
    }
    return std::clamp(p, 0.0, 1.0);
}

inline double weak_kind_mass(const DecayCascade& dc, WeakOutcomeKind kind) {
    double p = 0.0;
    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
        if (o.kind == kind) p += o.selected_mass;
    return std::clamp(p, 0.0, 1.0);
}

/// Absolute per-branch EC shell mass by daughter start level for a legacy
/// LevelScheme which has no explicit WeakOutcomeLaw. LevelDag::feed is
/// normalized conditional on entering the graph, so recover absolute starts
/// with entry_probability and restore the direct-to-ground terminal mass.
struct LegacyEcShellStarts {
    std::vector<double> k, l;
    double total_k = 0.0, total_l = 0.0;
};

inline LegacyEcShellStarts legacy_ec_shell_starts(const DecayCascade& dc,
                                                   const LevelDag& dag) {
    LegacyEcShellStarts result;
    if (!dag.valid) return result;
    result.k.assign(static_cast<std::size_t>(dag.n), 0.0);
    result.l.assign(static_cast<std::size_t>(dag.n), 0.0);
    const double ep = std::clamp(dc.level_scheme.entry_probability, 0.0, 1.0);
    for (int level = 0; level < dag.n; ++level) {
        double start = ep * dag.feed[static_cast<std::size_t>(level)];
        if (level == 0) start += 1.0 - ep;
        const CascadeLevel& lv = dc.level_scheme.levels[
            static_cast<std::size_t>(level)];
        result.k[static_cast<std::size_t>(level)] = start * lv.feed_ecK;
        result.l[static_cast<std::size_t>(level)] = start * lv.feed_ecL;
        result.total_k += result.k[static_cast<std::size_t>(level)];
        result.total_l += result.l[static_cast<std::size_t>(level)];
    }
    return result;
}

inline int residual_transition_of(const DecayCascade& dc, int member) {
    for (int r = 0; r < static_cast<int>(dc.residual_transitions.size()); ++r)
        if (dc.residual_transitions[static_cast<std::size_t>(r)].gamma_member == member)
            return r;
    return -1;
}

/// Deposit probability of transition t's summed x-ray channels into ANY window
/// (the sum used inside the survival factor and EC factor). `is_ic` selects the
/// per-transition IC yields; caller supplies the K/L line tables.
template <class T>
T xray_deposit_sum(double vac_prob, const std::vector<XrayLine>& lines,
                   const EffLookup<T>& eff) {
    if (vac_prob <= 0.0) return T(0.0);
    T s(0.0);
    for (const XrayLine& ln : lines) s += ln.weight * eff.tot(ln.energy);
    return vac_prob * s;
}

}  // namespace detail

// ===================== compute_cascade_analytic ==========================

template <class T>
std::vector<AnalyticPeakResultT<T>> compute_cascade_analytic(
    const std::vector<DecayCascade>& cascades,
    const std::vector<PeakWindow>& peaks,
    const EfficiencyProviderT<T>& provider,
    const AnalyticCascadeOptions& options) {

    using detail::XrayLine;
    using detail::kLOcc;
    using detail::kAnnih511;

    validate_cascade_member_references(cascades, "compute_cascade_analytic");

    std::vector<detail::AnalyticBranch> branches;
    branches.reserve(cascades.size());
    for (const auto& dc : cascades) branches.emplace_back(dc);

    const bool angular = options.apply_angular_correlation;

    // Angular a2/a4 for the coincidence link member_a -> member_b within branch b
    // (searched both directions); {0,0} if no correlated link.
    auto link_a2a4 = [&](int b, int ma, int mb) -> std::pair<double, double> {
        auto scan = [&](int from, int to) -> std::pair<double, double> {
            if (from < 0 || from >= (int)branches[b].dc->members.size())
                return {0.0, 0.0};
            for (const CoincidenceLink& lk : branches[b].dc->members[from].coincident)
                if (lk.partner == to && lk.has_correlation) return {lk.a2, lk.a4};
            return {0.0, 0.0};
        };
        auto r = scan(ma, mb);
        if (r.first == 0.0 && r.second == 0.0) r = scan(mb, ma);
        return r;
    };

    std::vector<AnalyticPeakResultT<T>> results;
    results.reserve(peaks.size());

    for (const PeakWindow& pk : peaks) {
        AnalyticPeakResultT<T> res;
        res.energy_keV = pk.energy_keV;
        const double win = pk.tolerance_keV > 0.0 ? pk.tolerance_keV
                                                  : options.window_tol_keV;
        std::set<double> unmatched;
        detail::EffLookup<T> eff{provider, unmatched};

        // --- Locate the window's emitted LINES (fitted-peak-area SF) ---
        // A window may contain more than one emitted line (e.g. Ba-133 79.6+81);
        // the SF is then defined on the summed peak AREA. Collect every emitted
        // line within the window, grouped by energy into "components". Each
        // component is branch-AVERAGED (the same gamma fed by e.g. a beta+ and an
        // EC branch has different coincident sets — only the beta+ carries the 511
        // pair). Per-decay window rates: R_no (no summing), R_out (after summing-
        // out), R_in (summing-in); SF = (R_out + R_in) / R_no.
        struct Prim { int b, m, t; double p_emit; };
        struct Component { double energy; std::vector<Prim> prims; double A = 0.0; };
        std::vector<Component> comps;
        for (int b = 0; b < (int)branches.size(); ++b) {
            const DecayCascade& dc = *branches[b].dc;
            const LevelDag& dag = branches[b].dag;
            for (int m = 0; m < (int)dc.members.size(); ++m) {
                const double E = dc.members[m].energy_keV;
                if (std::abs(E - pk.energy_keV) > win) continue;
                const int t = dag.valid ? dag.transition_of(m) : -1;
                // dag.pass() is conditional on the decay having entered the level
                // scheme (LevelDag normalizes feeding by its own sum), so an
                // absolute per-decay rate needs the scheme's entry probability.
                // This corrects the ABSOLUTE emission rate (Th-234's was 3.5x
                // high); it mostly cancels out of a reported summing FACTOR, where
                // it divides out between R_no and R_out/R_in, and survives only
                // where branches with different entry probabilities share a window
                // (measured: 0.000% on 16 nuclides except Ag-108 511 keV, -1.19%).
                double member_prob = dc.members[m].intensity;
                const int residual = t < 0
                    ? detail::residual_transition_of(dc, m) : -1;
                if (residual >= 0)
                    member_prob = dc.residual_transitions[
                        static_cast<std::size_t>(residual)].p_gamma;
                else if (t < 0 && dc.members[m].type == CascadeParticleType::Annih511 &&
                    dc.weak_outcome_law.usable())
                    member_prob = 2.0 * detail::weak_kind_mass(
                        dc, WeakOutcomeKind::Positron);  // two photons / positron
                else if (t < 0 &&
                         dc.members[m].type == CascadeParticleType::Annih511)
                    member_prob *= 2.0;
                const double p_emit =
                    (t >= 0) ? dc.branch_weight *
                                   detail::weak_subset_probability(dc, dag, {t}) *
                                   dag.ts[t].p_gamma
                             : dc.branch_weight * member_prob;
                if (p_emit <= 0.0) continue;
                Component* c = nullptr;
                for (auto& cc : comps)
                    if (std::abs(cc.energy - E) <= options.match_tol_keV) { c = &cc; break; }
                if (!c) { comps.push_back({E, {}, 0.0}); c = &comps.back(); }
                c->prims.push_back({b, m, t, p_emit});
                c->A += p_emit;
            }
        }
        if (comps.empty()) { res.found = false; results.push_back(res); continue; }
        res.found = true;
        // Rejected branches retain the legacy pairwise member fallback.  Its
        // coincidence links are not a coherent joint law database-wide, and
        // this analytic estimator does not enumerate their sum-fed pairs.  Keep
        // the historical numeric estimate, but make the affected peak explicit
        // rather than silently presenting it as complete.
        for (const Component& c : comps)
            for (const Prim& p : c.prims)
                if (!branches[p.b].dag.valid)
                    res.summing_model_complete = false;
        for (int b = 0; b < static_cast<int>(branches.size()); ++b) {
            const int excluded = !branches[b].dag.valid
                ? [&]() {
                    for (const Component& c : comps)
                        for (const Prim& p : c.prims)
                            if (p.b == b) return p.m;
                    return -1;
                  }()
                : -1;
            if (cascade_invalid_branch_can_feed_peak(
                    *branches[b].dc, pk.energy_keV, win, excluded))
                res.summing_model_complete = false;
        }

        // No-summing window rate and its emission-weighted mean FEP efficiency
        // (reduces to eps_FEP(peak) for a single line).
        T R_no(0.0);
        double A_total = 0.0;
        for (const auto& c : comps) { R_no += c.A * eff.fep(c.energy); A_total += c.A; }
        const T eff_no = (A_total > 0.0) ? T(R_no / A_total) : T(0.0);
        res.eff_no_summing = eff_no;
        res.eff_no_summing_unc = provider.fep_unc(pk.energy_keV);

        // ---------------- Summing-OUT (C_out) ----------------
        // C_out for one primary in branch b: the exact level-scheme survival DP
        // (or whole-scheme survival S_full for a non-scheme primary) times the
        // coincident annihilation-511 pair survival times any flat coincident
        // members. Averaged over branches below.
        auto cout_for_primary = [&](int b, int m, int t) -> T {
            const DecayCascade& dc = *branches[b].dc;
            const LevelDag& dag = branches[b].dag;
            const int Z = branches[b].Z;
            T c(1.0);
            if (dag.valid) {
                const std::vector<XrayLine> kl = detail::k_lines(Z);
                std::vector<std::vector<XrayLine>> ll = {detail::l_lines(Z, 0),
                                                         detail::l_lines(Z, 1),
                                                         detail::l_lines(Z, 2)};
                const int nts = (int)dag.ts.size();
                std::vector<T> f_t(nts, T(1.0));
                for (int i = 0; i < nts; ++i) {
                    const LevelDag::Tr& tr = dag.ts[i];
                    // gamma_member is -1 for an E0 transition, which has no photon
                    // member; only its conversion x-rays can deposit. Skip the
                    // efficiency lookup entirely rather than calling it at 0 keV --
                    // a provider records every miss, so a 0.0 would show up in
                    // AnalyticPeakResult::unmatched_energies and read as an
                    // incomplete provider.
                    const bool has_photon = tr.gamma_member >= 0
                        && static_cast<std::size_t>(tr.gamma_member)
                               < dc.members.size();
                    // Summing-OUT: the coincident gamma only needs ANY deposit
                    // (eps_tot, not on-axis-weighted), so its mutual angle with the
                    // FEP primary is loosely constrained and the W(theta) average is
                    // near 1 (measured negligible on Co-60). g = 1 here.
                    T dep = has_photon
                        ? tr.p_gamma * eff.tot(dc.members[tr.gamma_member].energy_keV)
                        : T(0.0);
                    dep += detail::xray_deposit_sum(tr.p_icK, kl, eff);
                    for (int s = 0; s < 3; ++s)
                        dep += detail::xray_deposit_sum(tr.p_icL[s], ll[s], eff);
                    f_t[i] = 1.0 - dep;
                }
                std::vector<T> f_EC(dag.n, T(1.0));
                T depK(0.0), depLmix(0.0);
                for (const XrayLine& ln : kl) depK += ln.weight * eff.tot(ln.energy);
                for (int s = 0; s < 3; ++s) {
                    T sL(0.0);
                    for (const XrayLine& ln : ll[s])
                        sL += ln.weight * eff.tot(ln.energy);
                    depLmix += kLOcc[s] * sL;
                }
                for (int l = 0; l < dag.n; ++l) {
                    const CascadeLevel& lv = dc.level_scheme.levels[l];
                    T dep = detail::xray_deposit_sum(lv.feed_ecK, kl, eff);
                    for (int s = 0; s < 3; ++s) {
                        T sL(0.0);
                        for (const XrayLine& ln : ll[s]) sL += ln.weight * eff.tot(ln.energy);
                        dep += lv.feed_ecL * kLOcc[s] * sL;
                    }
                    f_EC[l] = 1.0 - dep;
                }
                const std::vector<T> Dw = dag.survival_down(f_t);
                if (t >= 0) {
                    if (dc.weak_outcome_law.usable()) {
                        const double den = detail::weak_subset_probability(dc, dag, {t});
                        T num(0.0);
                        T pos_pair_survival = 1.0 - 2.0 * eff.tot(kAnnih511);
                        if (pos_pair_survival < 0.0) pos_pair_survival = T(0.0);
                        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                            if (o.fed_level >= 0 && o.fed_level < dag.n) {
                                const double reach = dag.pass_from_level(t, o.fed_level);
                                if (reach <= 0.0) continue;
                                T fs(1.0);
                                if (o.kind == WeakOutcomeKind::ElectronCapture)
                                    fs = fs - o.ec_K * depK - o.ec_L * depLmix;
                                else if (o.kind == WeakOutcomeKind::Positron)
                                    fs = fs * pos_pair_survival;
                                num += o.selected_mass * reach *
                                    dag.survival_out_from_level(t, o.fed_level,
                                                               f_t, Dw, fs);
                            } else {
                                T fs(1.0);
                                if (o.kind == WeakOutcomeKind::ElectronCapture)
                                    fs = fs - o.ec_K * depK - o.ec_L * depLmix;
                                else if (o.kind == WeakOutcomeKind::Positron)
                                    fs = fs * pos_pair_survival;
                                const std::vector<double> sw =
                                    detail::weak_unknown_start_weights(dc, dag);
                                for (int l = 0; l < dag.n; ++l) {
                                    const double reach = dag.pass_from_level(t, l);
                                    if (sw[l] <= 0.0 || reach <= 0.0) continue;
                                    num += o.selected_mass * sw[l] * reach *
                                        dag.survival_out_from_level(
                                            t, l, f_t, Dw, fs);
                                }
                            }
                        }
                        c *= (den > 0.0) ? T(num / den) : T(1.0);
                    } else {
                        const std::vector<T> Uw = dag.survival_up(f_t, f_EC);
                        c *= dag.survival_out(t, Dw, Uw);
                    }
                } else {
                    T s_full(0.0);  // whole-scheme survival (non-scheme primary)
                    for (int l = 0; l < dag.n; ++l)
                        s_full += dag.feed[l] * f_EC[l] * Dw[l];
                    if (dc.weak_outcome_law.usable()) {
                        T num(0.0);
                        double den = 0.0;
                        const bool primary_ann =
                            dc.members[m].type == CascadeParticleType::Annih511;
                        T pos_pair_survival = primary_ann
                            ? T(1.0)  // conditioned photon points at detector;
                                      // its back-to-back partner points away
                            : T(1.0) - 2.0 * eff.tot(kAnnih511);
                        if (pos_pair_survival < 0.0) pos_pair_survival = T(0.0);
                        const std::vector<double> sw =
                            detail::weak_unknown_start_weights(dc, dag);
                        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                            if (primary_ann && o.kind != WeakOutcomeKind::Positron)
                                continue;
                            den += o.selected_mass;
                            T so(0.0);
                            if (o.fed_level >= 0 && o.fed_level < dag.n) {
                                so = Dw[o.fed_level];
                            } else {
                                double assigned = 0.0;
                                for (int l = 0; l < dag.n; ++l) {
                                    so += sw[l] * Dw[l];
                                    assigned += sw[l];
                                }
                                so += std::max(0.0, 1.0 - assigned);
                            }
                            T fs(1.0);
                            if (o.kind == WeakOutcomeKind::ElectronCapture)
                                fs = fs - o.ec_K * depK - o.ec_L * depLmix;
                            else if (o.kind == WeakOutcomeKind::Positron)
                                fs = fs * pos_pair_survival;
                            num += o.selected_mass * so * fs;
                        }
                        c *= den > 0.0 ? T(num / den) : T(1.0);
                    } else {
                    // dag.feed is normalized, so s_full is the survival CONDITIONAL
                    // on the decay having entered the level scheme. A primary that
                    // is not itself a scheme transition (a 511, an x-ray, or a gamma
                    // with no level topology) is also emitted by decays that fed the
                    // daughter's ground state directly, and those have no scheme
                    // partner to sum with at all -- hence the (1 - ep) term. Without
                    // it this over-predicts summing-out by up to 1/ep and disagrees
                    // with the FullRealization walk, which does condition on entry.
                    //
                    // They are not partner-free, though: an electron capture makes
                    // its vacancy wherever it lands, so a ground-state-fed decay
                    // still carries the level-0 EC x-ray, and the walk fires it.
                    // f_EC[0] is that survival. dag.feed[0] is identically 0 (the
                    // ground state has no net out-flow), so s_full cannot carry it.
                    // Omitting it over-predicts C_out by up to (1-ep)*dep_EC0 --
                    // measured +6.7% (z = 9.2) on a synthetic Fe-daughter branch,
                    // and 614 of the database's schemes have (1-ep)*ec0 > 0.05.
                    const double ep = dc.level_scheme.entry_probability;
                    c *= (ep * s_full + (1.0 - ep) * f_EC[0]);
                    }
                }

                // Unmatched transitions are independent categorical residuals.
                // A residual primary is conditioned on its gamma alternative, so
                // its own conversion alternatives are impossible; every other
                // residual contributes its absolute gamma/x-ray deposit survival.
                const int primary_residual =
                    detail::residual_transition_of(dc, m);
                for (int rix = 0;
                     rix < static_cast<int>(dc.residual_transitions.size()); ++rix) {
                    if (rix == primary_residual) continue;
                    const ResidualTransition& r = dc.residual_transitions[
                        static_cast<std::size_t>(rix)];
                    T dep = (r.gamma_member >= 0
                             && static_cast<std::size_t>(r.gamma_member)
                                    < dc.members.size())
                        ? r.p_gamma * eff.tot(
                              dc.members[static_cast<std::size_t>(r.gamma_member)].energy_keV)
                        : T(0.0);
                    dep += detail::xray_deposit_sum(r.p_icK, kl, eff);
                    dep += detail::xray_deposit_sum(r.p_icL1, ll[0], eff);
                    dep += detail::xray_deposit_sum(r.p_icL2, ll[1], eff);
                    dep += detail::xray_deposit_sum(r.p_icL3, ll[2], eff);
                    c *= T(1.0) - dep;
                }
            }
            // Coincident annihilation 511 pairs (each fired member = TWO photons).
            // The two photons are BACK-TO-BACK, so for a single-sided detector at
            // most one can deposit: P(>=1 of the pair deposits) = eps_a + eps_b -
            // P(both) = 2*eps_tot(511), since P(both) = 0 (the partner points ~180
            // deg away from any convex detector subtending < 2pi). The independent
            // form 1-(1-eps)^2 = 2eps-eps^2 under-counts the removal by eps^2 (the
            // Na-22 1274 residual). For a 4pi/well detector both can deposit; that
            // regime would need the independent form, but the DRF is single-sided.
            const T et511 = eff.tot(kAnnih511);
            T surv_pair = 2.0 * et511;
            if (surv_pair > 1.0) surv_pair = T(1.0);
            for (int mm = 0; mm < (int)dc.members.size(); ++mm) {
                if (dc.members[mm].type != CascadeParticleType::Annih511) continue;
                if (dc.weak_outcome_law.usable())
                    continue;  // folded into the exact outcome mixture above
                if (mm == m)
                    c *= T(1.0);  // detected primary's partner points away
                else {
                    double pm = std::min(1.0, dc.members[mm].intensity);
                    if (dc.weak_outcome_law.usable()) {
                        if (t >= 0) {
                            const double den = detail::weak_subset_probability(
                                dc, dag, {t});
                            const WeakOutcomeKind pos = WeakOutcomeKind::Positron;
                            pm = den > 0.0 ? detail::weak_subset_probability(
                                dc, dag, {t}, &pos) / den : 0.0;
                        } else {
                            pm = detail::weak_kind_mass(
                                dc, WeakOutcomeKind::Positron);
                        }
                    }
                    c *= (1.0 - pm * surv_pair);
                }
            }
            // Flat coincident members when there is no level scheme (pairwise model).
            if (!dag.valid) {
                for (int mm = 0; mm < (int)dc.members.size(); ++mm) {
                    if (mm == m ||
                        dc.members[mm].type == CascadeParticleType::Annih511)
                        continue;
                    double pmate = dc.members[mm].intensity;
                    for (const CoincidenceLink& lk : dc.members[m].coincident)
                        if (lk.partner == mm) { pmate = lk.prob; break; }
                    c *= (1.0 - pmate * eff.tot(dc.members[mm].energy_keV));
                }
            }
            return c;
        };

        // Per-component summing-OUT, branch-averaged, then area-weighted into the
        // window rate R_out = Σ_comp A·eps_FEP(E_comp)·C_out_comp.
        T R_out(0.0);
        for (const auto& c : comps) {
            T num(0.0);
            for (const Prim& p : c.prims)
                num += p.p_emit * cout_for_primary(p.b, p.m, p.t);
            const T c_out_c = (c.A > 0.0) ? T(num / c.A) : T(1.0);
            R_out += c.A * eff.fep(c.energy) * c_out_c;
        }

        // ---------------- Summing-IN (C_in) ----------------
        // Absolute per-decay summing-in rate into the window
        //   R_in = Σ_pairs branch_weight_B · P_B(a ∧ b) · eps_FEP(a)·eps_FEP(b)·surv
        // (NO triple subtraction:
        // unlike the MC direct stream, the analytic estimator has none, so the
        // full joint P(a∧b) is used). Pairs are within one branch B (two members
        // from different branches are different decays and never coincide); B
        // ranges over all branches. Same-source (same transition / same vacancy)
        // pairs are mutually exclusive and skipped.
        T gain(0.0);

        // Unified emission-occurrence list per branch (built lazily per branch in
        // the pair loop). Occurrence kinds: gamma / IC-K / IC-L / EC-K / EC-L /
        // annih-511. `src` keys the source vacancy for exclusivity.
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

        for (int b = 0; b < (int)branches.size(); ++b) {
            const DecayCascade& dc = *branches[b].dc;
            const LevelDag& dag = branches[b].dag;
            const int Z = branches[b].Z;
            const int nts = dag.valid ? (int)dag.ts.size() : 0;
            const int ec_src = nts;  // shared EC exclusivity key (one capture/decay)
            const std::vector<XrayLine> kl = detail::k_lines(Z);
            const std::vector<std::vector<XrayLine>> ll = {detail::l_lines(Z, 0),
                                                           detail::l_lines(Z, 1),
                                                           detail::l_lines(Z, 2)};
            const detail::LegacyEcShellStarts legacy_ec =
                dc.weak_outcome_law.usable()
                    ? detail::LegacyEcShellStarts{}
                    : detail::legacy_ec_shell_starts(dc, dag);

            std::vector<Occ> occ;
            if (dag.valid) {
                for (int i = 0; i < nts; ++i) {
                    const LevelDag::Tr& tr = dag.ts[i];
                    if (tr.p_gamma > 0.0 && tr.gamma_member >= 0
                        && static_cast<std::size_t>(tr.gamma_member)
                               < dc.members.size())
                        occ.push_back({i, dc.members[tr.gamma_member].energy_keV,
                                       tr.p_gamma, i, false, false, false,
                                       tr.gamma_member});
                    if (tr.p_icK > 0.0)
                        for (const XrayLine& ln : kl)
                            occ.push_back({i, ln.energy, tr.p_icK * ln.weight, i,
                                           false, false, false, -1});
                    for (int s = 0; s < 3; ++s)
                        if (tr.p_icL[s] > 0.0)
                            for (const XrayLine& ln : ll[s])
                                occ.push_back({i, ln.energy, tr.p_icL[s] * ln.weight,
                                               i, false, false, false, -1});
                }
                // EC-feed vacancies (upstream of the whole cascade). Represented
                // once per branch; the joint uses up_reach_ec against the partner.
                double ecK_tot = 0.0, ecL_tot = 0.0;
                if (dc.weak_outcome_law.usable()) {
                    const WeakOutcomeKind ec = WeakOutcomeKind::ElectronCapture;
                    ecK_tot = detail::weak_subset_probability(dc, dag, {}, &ec, 1);
                    ecL_tot = detail::weak_subset_probability(dc, dag, {}, &ec, 2);
                } else {
                    ecK_tot = legacy_ec.total_k;
                    ecL_tot = legacy_ec.total_l;
                }
                if (ecK_tot > 0.0)
                    for (const XrayLine& ln : kl)
                        occ.push_back({-1, ln.energy, ln.weight, ec_src, true, false,
                                       false, -1});
                if (ecL_tot > 0.0)
                    for (int s = 0; s < 3; ++s)
                        for (const XrayLine& ln : ll[s])
                            occ.push_back({-1, ln.energy, kLOcc[s] * ln.weight, ec_src,
                                           true, true, false, -1});
            }
            // Categorical unmatched transitions coexist independently with the
            // matched graph. Their occurrence factors are already absolute per
            // branch firing (including the transition occurrence probability).
            for (int rix = 0;
                 rix < static_cast<int>(dc.residual_transitions.size()); ++rix) {
                const ResidualTransition& r = dc.residual_transitions[
                    static_cast<std::size_t>(rix)];
                const int src = nts + 1 + rix;
                if (r.p_gamma > 0.0 && r.gamma_member >= 0
                    && static_cast<std::size_t>(r.gamma_member)
                           < dc.members.size())
                    occ.push_back({-1,
                        dc.members[static_cast<std::size_t>(r.gamma_member)].energy_keV,
                        r.p_gamma, src, false, false, false, r.gamma_member, rix});
                if (r.p_icK > 0.0)
                    for (const XrayLine& ln : kl)
                        occ.push_back({-1, ln.energy, r.p_icK * ln.weight, src,
                                       false, false, false, -1, rix});
                const double pL[3] = {r.p_icL1, r.p_icL2, r.p_icL3};
                for (int s = 0; s < 3; ++s)
                    if (pL[s] > 0.0)
                        for (const XrayLine& ln : ll[s])
                            occ.push_back({-1, ln.energy, pL[s] * ln.weight, src,
                                           false, false, false, -1, rix});
            }
            // Annihilation 511 members (independent of the gamma path).
            for (int m = 0; m < (int)dc.members.size(); ++m) {
                if (dc.members[m].type != CascadeParticleType::Annih511) continue;
                occ.push_back({-1, dc.members[m].energy_keV,
                               2.0 * (dc.weak_outcome_law.usable() ? 1.0
                                      : dc.members[m].intensity),
                               nts + 1 + static_cast<int>(dc.residual_transitions.size()) + m,
                               false, false, true, m});
            }
            // A sum-in partner cannot be the primary line itself (an occurrence at
            // any component line is a peak line itself, already in R_no).
            auto is_peak_line = [&](const Occ& o) {
                for (const auto& c : comps)
                    if (std::abs(o.energy - c.energy) <= options.match_tol_keV) return true;
                return false;
            };

            // Unnormalised joint P_B(x ∧ y) for two occurrences of this branch.
            auto joint = [&](const Occ& x, const Occ& y) -> double {
                if (x.src == y.src) return 0.0;  // same transition / same vacancy
                // annih vs annih: two 511s of one positron are one physical event
                // (handled specially below); different-member 511s are rare -> skip.
                if (x.is_annih && y.is_annih) return 0.0;
                if (x.is_ec && y.is_ec) return 0.0;  // one capture => one vacancy
                auto ecmarg = [&](bool isL) {
                    if (dc.weak_outcome_law.usable()) {
                        const WeakOutcomeKind ec = WeakOutcomeKind::ElectronCapture;
                        return detail::weak_subset_probability(
                            dc, dag, {}, &ec, isL ? 2 : 1);
                    }
                    return isL ? legacy_ec.total_l : legacy_ec.total_k;
                };
                // Residual transition outcomes are independent of the matched
                // topology and weak mode. Their factors are absolute, whereas
                // DAG gamma/IC factors remain conditional on reaching `trans`.
                if (x.residual >= 0 || y.residual >= 0) {
                    const Occ& r = x.residual >= 0 ? x : y;
                    const Occ& o = x.residual >= 0 ? y : x;
                    if (o.residual >= 0)
                        return r.factor * o.factor;
                    if (o.trans >= 0)
                        return detail::weak_subset_probability(dc, dag, {o.trans})
                             * r.factor * o.factor;
                    if (o.is_ec)
                        return ecmarg(o.ec_isL) * r.factor * o.factor;
                    if (o.is_annih) {
                        const double pm = dc.weak_outcome_law.usable()
                            ? detail::weak_kind_mass(dc, WeakOutcomeKind::Positron)
                            : 1.0;  // legacy annihilation intensity is in factor
                        return pm * r.factor * o.factor;
                    }
                    return 0.0;
                }
                // Order so `a` is the transition-attached one when possible.
                const Occ& t = (x.trans >= 0) ? x : y;
                const Occ& o = (x.trans >= 0) ? y : x;
                if (t.trans < 0) {
                    // neither is a real transition: (ec/annih) x (ec/annih)
                    if ((x.is_ec && y.is_annih) || (x.is_annih && y.is_ec)) {
                        if (dc.weak_outcome_law.usable()) return 0.0;
                        return ecmarg(x.is_ec ? x.ec_isL : y.ec_isL) * x.factor * y.factor;
                    }
                    return 0.0;
                }
                if (o.trans >= 0)            // gamma/IC × gamma/IC
                    return detail::weak_subset_probability(
                               dc, dag, {t.trans, o.trans}) * t.factor * o.factor;
                if (o.is_ec) {               // gamma/IC × EC-feed vacancy
                    if (dc.weak_outcome_law.usable()) {
                        const WeakOutcomeKind ec = WeakOutcomeKind::ElectronCapture;
                        return detail::weak_subset_probability(
                                   dc, dag, {t.trans}, &ec, o.ec_isL ? 2 : 1) *
                               t.factor * o.factor;
                    }
                    return dc.level_scheme.entry_probability *
                           dag.up_reach_ec(t.trans, dc.level_scheme, o.ec_isL) *
                           t.factor * o.factor;
                }
                // gamma/IC × annih-511: exact beta+ mode joint when available.
                if (dc.weak_outcome_law.usable()) {
                    const WeakOutcomeKind pos = WeakOutcomeKind::Positron;
                    return detail::weak_subset_probability(
                               dc, dag, {t.trans}, &pos) * t.factor * o.factor;
                }
                return dc.level_scheme.entry_probability * dag.pass(t.trans) *
                       t.factor * o.factor;
            };

            // Per-shell "any-deposit" sums for the coincident-survival factor
            // below: Σ_lines (omega·branch)·eps_tot(line), constant per branch.
            T sumK(0.0);
            for (const XrayLine& ln : kl) sumK += ln.weight * eff.tot(ln.energy);
            T sumL[3] = {T(0.0), T(0.0), T(0.0)};
            T sumLmix(0.0);
            for (int s = 0; s < 3; ++s) {
                for (const XrayLine& ln : ll[s]) sumL[s] += ln.weight * eff.tot(ln.energy);
                sumLmix += kLOcc[s] * sumL[s];
            }
            const T et511 = eff.tot(kAnnih511);

            // Coincident-survival factor for a summing-in pair (x,y): the pair
            // already deposits E_x+E_y ≈ peak, so ANY OTHER coincident emission
            // depositing shifts the sum out of the window (summing-out applied to
            // the sum-fed channel). Product over the members co-emitted with the
            // pair, each conditioned on the pair via the level-scheme joint. This
            // is the near-field term the bare eps_FEP(x)·eps_FEP(y) misses (e.g.
            // Co-57 122+14→136 co-emits an EC Fe Kα that knocks the event out).
            auto pair_coincident_survival = [&](const Occ& x, const Occ& y) -> T {
                std::vector<int> defs;
                if (x.trans >= 0) defs.push_back(x.trans);
                if (y.trans >= 0) defs.push_back(y.trans);

                T residual_survival(1.0);
                for (int rix = 0;
                     rix < static_cast<int>(dc.residual_transitions.size()); ++rix) {
                    if (rix == x.residual || rix == y.residual) continue;
                    const ResidualTransition& r = dc.residual_transitions[
                        static_cast<std::size_t>(rix)];
                    T dep = (r.gamma_member >= 0
                             && static_cast<std::size_t>(r.gamma_member)
                                    < dc.members.size())
                        ? r.p_gamma * eff.tot(dc.members[
                              static_cast<std::size_t>(r.gamma_member)].energy_keV)
                        : T(0.0);
                    dep += r.p_icK * sumK;
                    dep += r.p_icL1 * sumL[0] + r.p_icL2 * sumL[1]
                         + r.p_icL3 * sumL[2];
                    residual_survival *= T(1.0) - dep;
                }
                std::vector<T> pair_f_t(static_cast<std::size_t>(nts), T(1.0));
                for (int m = 0; m < nts; ++m) {
                    const LevelDag::Tr& tr = dag.ts[m];
                    T dep = (tr.gamma_member >= 0
                             && static_cast<std::size_t>(tr.gamma_member)
                                    < dc.members.size())
                        ? tr.p_gamma * eff.tot(
                              dc.members[tr.gamma_member].energy_keV)
                        : T(0.0);
                    dep += tr.p_icK * sumK;
                    for (int s = 0; s < 3; ++s)
                        dep += tr.p_icL[s] * sumL[s];
                    pair_f_t[static_cast<std::size_t>(m)] = T(1.0) - dep;
                }
                if (defs.empty()) {
                    const std::vector<T> down = dag.survival_down(pair_f_t);
                    T depK(0.0), depL(0.0);
                    for (const XrayLine& ln : kl)
                        depK += ln.weight * eff.tot(ln.energy);
                    for (int s = 0; s < 3; ++s)
                        depL += kLOcc[s] * sumL[s];
                    T graph_survival(0.0);
                    const bool constrained_ec = x.is_ec || y.is_ec;
                    const bool constrained_pos = x.is_annih || y.is_annih;
                    const bool constrained_ecL = x.is_ec ? x.ec_isL : y.ec_isL;
                    if (dc.weak_outcome_law.usable()) {
                        const std::vector<double> sw =
                            detail::weak_unknown_start_weights(dc, dag);
                        double den = 0.0;
                        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                            if (constrained_ec
                                && o.kind != WeakOutcomeKind::ElectronCapture)
                                continue;
                            if (constrained_pos
                                && o.kind != WeakOutcomeKind::Positron)
                                continue;
                            const double mode_weight = constrained_ec
                                ? (constrained_ecL ? o.ec_L : o.ec_K) : 1.0;
                            const double outcome_weight =
                                o.selected_mass * mode_weight;
                            if (!(outcome_weight > 0.0)) continue;
                            T so(0.0);
                            if (o.fed_level >= 0 && o.fed_level < dag.n) {
                                so = down[o.fed_level];
                            } else {
                                double assigned = 0.0;
                                for (int l = 0; l < dag.n; ++l) {
                                    so += sw[l] * down[l];
                                    assigned += sw[l];
                                }
                                so += std::max(0.0, 1.0 - assigned);
                            }
                            if (!constrained_ec
                                && o.kind == WeakOutcomeKind::ElectronCapture)
                                so *= T(1.0) - o.ec_K * depK - o.ec_L * depL;
                            else if (!constrained_pos
                                     && o.kind == WeakOutcomeKind::Positron) {
                                T pair_dep = 2.0 * et511;
                                if (pair_dep > 1.0) pair_dep = T(1.0);
                                so *= T(1.0) - pair_dep;
                            }
                            graph_survival += outcome_weight * so;
                            den += outcome_weight;
                        }
                        if (den > 0.0) graph_survival *= 1.0 / den;
                    } else {
                        const double ep = dc.level_scheme.entry_probability;
                        if (constrained_ec) {
                            const std::vector<double>& starts = constrained_ecL
                                ? legacy_ec.l : legacy_ec.k;
                            const double den = constrained_ecL
                                ? legacy_ec.total_l : legacy_ec.total_k;
                            if (den > 0.0) {
                                for (int l = 0; l < dag.n; ++l)
                                    graph_survival += starts[
                                        static_cast<std::size_t>(l)] * down[l];
                                graph_survival *= 1.0 / den;
                            } else {
                                graph_survival = T(1.0);
                            }
                        } else {
                            for (int l = 0; l < dag.n; ++l) {
                                const CascadeLevel& lv = dc.level_scheme.levels[l];
                                const T fec = T(1.0) - lv.feed_ecK * depK
                                    - lv.feed_ecL * depL;
                                graph_survival += ep * dag.feed[l] * fec * down[l];
                            }
                            if (dag.n > 0) {
                                const CascadeLevel& g = dc.level_scheme.levels[0];
                                const T fec0 = T(1.0) - g.feed_ecK * depK
                                    - g.feed_ecL * depL;
                                graph_survival += (1.0 - ep) * fec0;
                            }
                        }
                    }
                    return residual_survival * graph_survival;
                }
                auto transition_survival = [&](const WeakOutcome* outcome,
                                               double base) -> T {
                    if (!(base > 0.0)) return T(1.0);
                    T joint_survival(0.0);
                    if (outcome && outcome->fed_level >= 0 &&
                        outcome->fed_level < dag.n) {
                        joint_survival = dag.constrained_survival_from_level(
                            outcome->fed_level, defs, pair_f_t);
                    } else if (outcome) {
                        const std::vector<double> sw =
                            detail::weak_unknown_start_weights(dc, dag);
                        for (int l = 0; l < dag.n; ++l)
                            if (sw[l] > 0.0)
                                joint_survival += sw[l] *
                                    dag.constrained_survival_from_level(
                                        l, defs, pair_f_t);
                    } else {
                        for (int l = 0; l < dag.n; ++l)
                            if (dag.feed[l] > 0.0)
                                joint_survival += dag.feed[l] *
                                    dag.constrained_survival_from_level(
                                        l, defs, pair_f_t);
                    }
                    return joint_survival / base;
                };

                if (dc.weak_outcome_law.usable()) {
                    T num(0.0);
                    double den = 0.0;
                    const bool constrained_ec = x.is_ec || y.is_ec;
                    const bool constrained_pos = x.is_annih || y.is_annih;
                    const bool ecL = x.is_ec ? x.ec_isL : y.ec_isL;
                    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                        if (constrained_ec && o.kind != WeakOutcomeKind::ElectronCapture)
                            continue;
                        if (constrained_pos && o.kind != WeakOutcomeKind::Positron)
                            continue;
                        const double base =
                            detail::weak_outcome_reach(dc, dag, o, defs);
                        double mode_weight = 1.0;
                        if (constrained_ec) mode_weight = ecL ? o.ec_L : o.ec_K;
                        const double w = o.selected_mass * base * mode_weight;
                        if (!(w > 0.0)) continue;
                        T so = transition_survival(&o, base);
                        if (!constrained_ec && !constrained_pos) {
                            if (o.kind == WeakOutcomeKind::ElectronCapture)
                                so *= 1.0 - o.ec_K * sumK - o.ec_L * sumLmix;
                            else if (o.kind == WeakOutcomeKind::Positron) {
                                T pair_dep = 2.0 * et511;
                                if (pair_dep > 1.0) pair_dep = T(1.0);
                                so *= 1.0 - pair_dep;
                            }
                        }
                        num += w * so;
                        den += w;
                    }
                    return residual_survival
                         * (den > 0.0 ? T(num / den) : T(1.0));
                }

                const double base = dag.n_subset_joint(defs);
                if (!(base > 0.0)) return T(1.0);
                const bool constrained_ec = x.is_ec || y.is_ec;
                const bool constrained_ecL = x.is_ec ? x.ec_isL : y.ec_isL;
                T surv(1.0);
                if (constrained_ec) {
                    const std::vector<double>& starts = constrained_ecL
                        ? legacy_ec.l : legacy_ec.k;
                    T num(0.0);
                    double den = 0.0;
                    for (int l = 0; l < dag.n; ++l) {
                        const double w = starts[static_cast<std::size_t>(l)];
                        if (!(w > 0.0)) continue;
                        const double reach = dag.subset_from_level(l, defs);
                        den += w * reach;
                        num += w * dag.constrained_survival_from_level(
                            l, defs, pair_f_t);
                    }
                    if (den > 0.0) surv = num / den;
                } else {
                    surv = transition_survival(nullptr, base);
                }
                // For a detected annihilation photon the back-to-back partner
                // points away in the single-sided geometry assumed here.
                if (!(x.is_ec || y.is_ec || x.is_annih || y.is_annih)) {
                    int a_hi = defs[0];
                    for (int t : defs)
                        if (dag.ts[t].from > dag.ts[a_hi].from) a_hi = t;
                    const double pass_hi = dag.pass(a_hi);
                    if (pass_hi > 0.0) {
                        const double pecK = dag.up_reach_ec(
                            a_hi, dc.level_scheme, false) / pass_hi;
                        const double pecL = dag.up_reach_ec(
                            a_hi, dc.level_scheme, true) / pass_hi;
                        surv *= (1.0 - pecK * sumK - pecL * sumLmix);
                    }
                }
                return residual_survival * surv;
            };

            // Pairs. Sort occurrences by energy; for each x, binary-search the
            // partner-energy window [peak-win - E_x, peak+win - E_x] so only pairs
            // that actually sum into the window are visited — O(n log n) instead
            // of O(n^2). Matters for dense schemes (Eu-152/Am-241) whose occurrence
            // lists include the ~25 expanded vacancy x-ray lines per transition.
            std::vector<int> ord(occ.size());
            for (int i = 0; i < (int)occ.size(); ++i) ord[i] = i;
            std::sort(ord.begin(), ord.end(),
                      [&](int p, int q) { return occ[p].energy < occ[q].energy; });
            std::vector<double> es(ord.size());
            for (int k = 0; k < (int)ord.size(); ++k) es[k] = occ[ord[k]].energy;
            const double plo = pk.energy_keV - win, phi = pk.energy_keV + win;
            for (int ia = 0; ia < (int)ord.size(); ++ia) {
                const Occ& x = occ[ord[ia]];
                if (is_peak_line(x)) continue;
                int jb = static_cast<int>(
                    std::lower_bound(es.begin() + ia + 1, es.end(), plo - x.energy)
                    - es.begin());
                const int je = static_cast<int>(
                    std::upper_bound(es.begin() + ia + 1, es.end(), phi - x.energy)
                    - es.begin());
                for (; jb < je; ++jb) {
                    const Occ& y = occ[ord[jb]];
                    if (is_peak_line(y)) continue;
                    const double J = joint(x, y);
                    if (J <= 0.0) continue;
                    double g = 1.0;
                    if (angular && x.member >= 0 && y.member >= 0 &&
                        !x.is_annih && !y.is_annih) {
                        auto a = link_a2a4(b, x.member, y.member);
                        if (a.first != 0.0 || a.second != 0.0)
                            g = detail::angular_g_w0(a.first, a.second);
                    }
                    const T surv = pair_coincident_survival(x, y);
                    gain += dc.branch_weight *
                            J * eff.fep(x.energy) *
                            eff.fep(y.energy) * g * surv;
                }
            }

            // No 511+511 -> 1022 term: this analytic path explicitly assumes a
            // convex, single-sided detector subtending <2pi.  The photons are
            // back-to-back, so if one enters that detector the other points away.
            // A well/4pi geometry requires an angular-response provider and is
            // outside this scalar-efficiency model's contract.

            // Triples (optional): three GAMMA lines on distinct transitions whose
            // energies sum into the window. Restricted to gammas (member >= 0):
            // triples involving expanded vacancy x-ray LINES (or EC/annih) are
            // eps^3-and-branch-small AND would make this O(n^3) over the ~25
            // x-ray lines per converting transition — seconds on dense schemes
            // (Eu-152, Am-241). Sorted + early-pruned on the running partial sum
            // so the inner loops collapse once the sum passes the window.
            if (options.enumerate_triples) {
                std::vector<int> real;
                for (int i = 0; i < (int)occ.size(); ++i)
                    if (occ[i].trans >= 0 && occ[i].member >= 0 && !is_peak_line(occ[i]))
                        real.push_back(i);
                std::sort(real.begin(), real.end(),
                          [&](int p, int q) { return occ[p].energy < occ[q].energy; });
                const double thi = pk.energy_keV + win, tlo = pk.energy_keV - win;
                for (int a = 0; a < (int)real.size(); ++a) {
                    const Occ& X = occ[real[a]];
                    if (X.energy >= thi) break;                    // sorted: rest larger
                    for (int b2 = a + 1; b2 < (int)real.size(); ++b2) {
                        const Occ& Y = occ[real[b2]];
                        const double xy = X.energy + Y.energy;
                        if (xy >= thi) break;                      // Z2 > 0 -> can't return
                        for (int c = b2 + 1; c < (int)real.size(); ++c) {
                            const Occ& Z2 = occ[real[c]];
                            const double s = xy + Z2.energy;
                            if (s > thi) break;                    // sorted in c: rest larger
                            if (s < tlo) continue;
                            if (X.src == Y.src || X.src == Z2.src || Y.src == Z2.src)
                                continue;
                            const double J = detail::weak_subset_probability(dc, dag,
                                {X.trans, Y.trans, Z2.trans}) *
                                X.factor * Y.factor * Z2.factor;
                            if (J <= 0.0) continue;
                            gain += dc.branch_weight * J *
                                    eff.fep(X.energy) *
                                    eff.fep(Y.energy) * eff.fep(Z2.energy);
                        }
                    }
                }
            }
        }

        // Fitted-peak-area summing factor: SF = (R_out + R_in) / R_no. c_out /
        // c_in are the window-rate-normalised summing-out / summing-in fractions
        // (for a single line these reduce to the per-line C_out and C_in).
        const T R_in = gain;  // absolute per-decay summing-in rate into window
        const bool has_rate = (R_no > 1e-30);
        const T sf = has_rate ? T((R_out + R_in) / R_no) : T(1.0);
        res.c_out = has_rate ? T(R_out / R_no) : T(1.0);
        res.c_in = has_rate ? T(R_in / R_no) : T(0.0);
        res.c_net = sf;
        res.summing_factor = sf;
        res.eff_with_summing = eff_no * sf;

        // First-order uncertainty: mean FEP uncertainty scaled by the summing-out
        // fraction (the summing-in gain's own uncertainty is second order here).
        const T var_fep = res.eff_no_summing_unc * res.eff_no_summing_unc;
        const T cc_var = res.c_out * res.c_out * var_fep;
        using std::sqrt;
        res.eff_with_summing_unc = (cc_var > 0.0) ? T(sqrt(cc_var)) : T(0.0);
        res.summing_factor_unc =
            (eff_no > 1e-15) ? T(res.eff_with_summing_unc / eff_no) : T(0.0);

        res.unmatched_energies.assign(unmatched.begin(), unmatched.end());
        results.push_back(std::move(res));
    }

    return results;
}

}  // namespace ceelo
