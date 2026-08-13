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

// Flattened daughter level scheme as a transition DAG, with the exact visit /
// reach dynamic program shared by the cascade level-scheme routines (the MC
// estimators cascade_level_pmate / cascade_sum_pair_channels /
// cascade_level_vacancies, and the analytic path compute_cascade_analytic).
// Mirrors the LevelWalk formulation validated in the cascade design study. All joint
// probabilities are Markov on the (descending) level order: a decay realization
// is exactly one descending path from the fed level to ground, so P(any subset
// of transitions passed) is closed-form.
//
// Extracted from the anonymous namespace of src/efficiency/EfficiencyCalculator.cpp
// so both the MC estimator and the analytic summing path can use it. Pure
// combinatorics: depends only on CascadeTypes.h (no cross-sections, no
// efficiency provider). The survival-factor DP below takes PRECOMPUTED per-
// transition / per-level factors so callers own the physics (ε lookups, x-ray
// line expansion) while this header owns the level-scheme bookkeeping.

#include "cascade/CascadeTypes.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ceelo {

/// Reject contradictory host-built cascade records before an estimator can
/// dereference a photon member.  A memberless E0 is represented by
/// gamma_member == -1 and p_gamma == 0; every positive photon outcome must name
/// an in-range Gamma member.  SandiaDecay-built cascades satisfy this by
/// construction, but DecayCascade is also a public aggregate and callers may
/// construct it directly.
inline void validate_cascade_member_references(
    const std::vector<DecayCascade>& cascades, const char* consumer)
{
    auto fail = [&](std::size_t branch, const std::string& where,
                    int member, double p_gamma) {
        std::ostringstream msg;
        msg << consumer << ": invalid gamma member in branch " << branch
            << ", " << where << " (gamma_member=" << member
            << ", p_gamma=" << p_gamma << ')';
        throw std::invalid_argument(msg.str());
    };
    auto valid_gamma = [](const DecayCascade& dc, int member) {
        return member >= 0
            && static_cast<std::size_t>(member) < dc.members.size()
            && dc.members[static_cast<std::size_t>(member)].type
                   == CascadeParticleType::Gamma;
    };

    for (std::size_t b = 0; b < cascades.size(); ++b) {
        const DecayCascade& dc = cascades[b];
        if (dc.level_scheme.valid) {
            for (std::size_t l = 0; l < dc.level_scheme.levels.size(); ++l) {
                const auto& out = dc.level_scheme.levels[l].out;
                for (std::size_t e = 0; e < out.size(); ++e) {
                    const LevelOutTransition& tr = out[e];
                    if (tr.gamma_member >= 0) {
                        if (!valid_gamma(dc, tr.gamma_member))
                            fail(b, "level " + std::to_string(l) + " edge "
                                        + std::to_string(e),
                                 tr.gamma_member, tr.p_gamma);
                    } else if (tr.gamma_member != -1 || tr.p_gamma > 0.0) {
                        fail(b, "level " + std::to_string(l) + " edge "
                                    + std::to_string(e),
                             tr.gamma_member, tr.p_gamma);
                    }
                }
            }
        }
        for (std::size_t r = 0; r < dc.residual_transitions.size(); ++r) {
            const ResidualTransition& tr = dc.residual_transitions[r];
            if (tr.gamma_member >= 0) {
                if (!valid_gamma(dc, tr.gamma_member))
                    fail(b, "residual transition " + std::to_string(r),
                         tr.gamma_member, tr.p_gamma);
            } else if (tr.gamma_member != -1 || tr.p_gamma > 0.0) {
                fail(b, "residual transition " + std::to_string(r),
                     tr.gamma_member, tr.p_gamma);
            }
        }
        for (std::size_t k = 0; k < dc.k_vacancies.size(); ++k) {
            const int member = dc.k_vacancies[k].gamma_member;
            if (member >= 0 && !valid_gamma(dc, member))
                fail(b, "K-vacancy source " + std::to_string(k), member, 0.0);
            if (member < -1)
                fail(b, "K-vacancy source " + std::to_string(k), member, 0.0);
        }
        for (std::size_t v = 0; v < dc.vacancy_groups.size(); ++v) {
            const VacancyGroup& group = dc.vacancy_groups[v];
            const int member = group.gamma_member;
            if (member >= 0 && !valid_gamma(dc, member))
                fail(b, "vacancy group " + std::to_string(v), member, 0.0);
            if (member < -1)
                fail(b, "vacancy group " + std::to_string(v), member, 0.0);
        }
    }
}

struct LevelDag {
    struct Tr {
        int from, to, gamma_member;
        double b, p_gamma, p_icK, p_icL[3];
        /// Transition energy. Carried here rather than recovered from
        /// `members[gamma_member]` because an E0 has no photon member
        /// (gamma_member == -1) yet still emits conversion electrons of this
        /// energy; looking it up through the member silently yielded 0 and
        /// dropped them.
        double gamma_keV;
        /// Below-pair E0 release approximated as a shell-unresolved conversion
        /// electron. It creates no modeled atomic vacancy.
        double p_icUnresolved;
        /// Above-pair E0 remainder carried explicitly but not transported.
        double p_unmodeled;
    };
    std::vector<Tr> ts;
    std::vector<double> feed;              ///< normalized direct feeding per level
    std::vector<double> V;                 ///< P(a realization visits level l)
    std::vector<std::vector<double>> R;    ///< R[i][l] = P(transition i passed | at level l)
    int n = 0;
    bool valid = false;

    explicit LevelDag(const LevelScheme& ls) {
        if (!ls.valid) return;
        n = static_cast<int>(ls.levels.size());

        // Flatten into transitions (to_level < from_level), normalized branching b.
        for (int l = 0; l < n; ++l) {
            double wsum = 0.0;
            for (const auto& t : ls.levels[l].out) wsum += t.weight;
            if (wsum <= 0.0) continue;
            for (const auto& t : ls.levels[l].out)
                if (t.to_level >= 0 && t.to_level < l)
                    ts.push_back({l, t.to_level, t.gamma_member, t.weight / wsum,
                                  t.p_gamma, t.p_icK,
                                  {t.p_icL1, t.p_icL2, t.p_icL3}, t.gamma_keV,
                                  t.p_icUnresolved, t.p_unmodeled});
        }

        // feed[l], then V[l] by propagating high -> low along the DAG.
        double fsum = 0.0;
        for (const auto& lv : ls.levels) fsum += lv.feeding;
        feed.assign(n, 0.0);
        for (int l = 0; l < n; ++l)
            feed[l] = (fsum > 0.0) ? ls.levels[l].feeding / fsum : 0.0;
        V = feed;
        for (int l = n - 1; l >= 0; --l)
            for (const auto& t : ts)
                if (t.from == l) V[t.to] += V[l] * t.b;

        // R[i][l] = reach(i)[l], the same low->high recurrence per target.
        R.assign(ts.size(), std::vector<double>(n, 0.0));
        for (int i = 0; i < static_cast<int>(ts.size()); ++i)
            for (int l = 1; l < n; ++l)
                for (int j = 0; j < static_cast<int>(ts.size()); ++j)
                    if (ts[j].from == l)
                        R[i][l] += ts[j].b * ((j == i) ? 1.0 : R[i][ts[j].to]);
        valid = true;
    }

    double pass(int i) const { return V[ts[i].from] * ts[i].b; }

    /// P(transition `i` is passed | the realization starts at `level`).  Unlike
    /// pass(), this does not average over or normalize LevelScheme::feeding.  It
    /// is the primitive needed by the explicit weak-decay outcome law, whose
    /// EC/positron records already identify the fed level (including non-zero
    /// terminal isomer levels).
    double pass_from_level(int i, int level) const {
        if (i < 0 || i >= static_cast<int>(ts.size()) ||
            level < 0 || level >= n)
            return 0.0;
        return R[i][level];
    }

    /// P(all transitions in `subset` are passed | start at `level`).  The fixed-
    /// start counterpart of n_subset_joint(); it deliberately contains no
    /// entry_probability or branch_weight, so callers apply those exactly once.
    double subset_from_level(int level, std::vector<int> subset) const {
        if (level < 0 || level >= n) return 0.0;
        if (subset.empty()) return 1.0;
        for (int t : subset)
            if (t < 0 || t >= static_cast<int>(ts.size())) return 0.0;
        std::sort(subset.begin(), subset.end(),
                  [this](int x, int y) { return ts[x].from > ts[y].from; });
        double p = R[subset[0]][level];
        for (std::size_t k = 1; k < subset.size() && p > 0.0; ++k)
            p *= R[subset[k]][ts[subset[k - 1]].to];
        return p;
    }

    /// Probability-weighted survival of every transition not in `required`,
    /// for a path starting at `level` and constrained to pass all required
    /// transitions. Unlike multiplying conditional per-transition marginals,
    /// this dynamic program preserves path exclusivity: two alternative exits
    /// from one level contribute a sum, never an impossible product cross-term.
    /// The return value is the JOINT numerator
    ///   P(required path) * E[product(f_other) | required path];
    /// divide by subset_from_level(level, required) for the conditional
    /// survival. Required transitions themselves get factor one because their
    /// deposits define the sum-fed channel rather than summing it out.
    template <class T>
    T constrained_survival_from_level(int level, std::vector<int> required,
                                      const std::vector<T>& f_t) const {
        if (level < 0 || level >= n || f_t.size() != ts.size()) return T(0.0);
        for (int t : required)
            if (t < 0 || t >= static_cast<int>(ts.size())) return T(0.0);
        std::sort(required.begin(), required.end(),
                  [this](int x, int y) { return ts[x].from > ts[y].from; });
        const int nr = static_cast<int>(required.size());
        std::vector<std::vector<T>> q(
            static_cast<std::size_t>(nr) + 1,
            std::vector<T>(static_cast<std::size_t>(n), T(0.0)));
        for (int l = 0; l < n; ++l) {
            bool any = false;
            for (const Tr& tr : ts)
                if (tr.from == l) { any = true; break; }
            if (!any) {
                q[static_cast<std::size_t>(nr)][static_cast<std::size_t>(l)] = T(1.0);
                continue;
            }
            for (int k = nr; k >= 0; --k) {
                T s(0.0);
                for (int j = 0; j < static_cast<int>(ts.size()); ++j) {
                    if (ts[j].from != l) continue;
                    const auto it = std::find(required.begin(), required.end(), j);
                    if (it != required.end()) {
                        const int rj = static_cast<int>(it - required.begin());
                        if (rj == k)
                            s += ts[j].b * q[static_cast<std::size_t>(k + 1)]
                                              [static_cast<std::size_t>(ts[j].to)];
                    } else {
                        s += ts[j].b * f_t[static_cast<std::size_t>(j)] *
                             q[static_cast<std::size_t>(k)]
                              [static_cast<std::size_t>(ts[j].to)];
                    }
                }
                q[static_cast<std::size_t>(k)][static_cast<std::size_t>(l)] = s;
            }
        }
        return q[0][static_cast<std::size_t>(level)];
    }

    /// Conditional survival of every transition other than primary `p`, for a
    /// realization known to start at `level` and to pass `p`.  `Dw` contains the
    /// already-computed survival below each level.  `start_factor` is the
    /// selected weak outcome's atomic-vacancy survival (one EC shell, or one for
    /// beta+/other).  This is the fixed-start analogue of survival_out().
    template <class T>
    T survival_out_from_level(int p, int level, const std::vector<T>& f_t,
                              const std::vector<T>& Dw,
                              const T& start_factor) const {
        if (p < 0 || p >= static_cast<int>(ts.size()) ||
            level < 0 || level >= n || R[p][level] <= 0.0)
            return T(1.0);
        std::vector<T> q(n, T(0.0));
        for (int l = 1; l < n; ++l)
            for (int j = 0; j < static_cast<int>(ts.size()); ++j)
                if (ts[j].from == l)
                    q[l] += ts[j].b * ((j == p) ? Dw[ts[j].to]
                                                   : f_t[j] * q[ts[j].to]);
        return start_factor * q[level] / R[p][level];
    }

    /// First transition emitting gamma member `m`, or -1 (each gamma maps to one).
    int transition_of(int m) const {
        for (int i = 0; i < static_cast<int>(ts.size()); ++i)
            if (ts[i].gamma_member == m) return i;
        return -1;
    }

    /// P(both transitions i and j are passed). 0 when they leave the same level
    /// (mutually exclusive) -- the exclusivity fix.
    double pair_joint(int i, int j) const {
        if (i == j) return pass(i);
        return pass(i) * R[j][ts[i].to] + pass(j) * R[i][ts[j].to];
    }

    /// P(all three transitions passed) in one downward walk. Sort by descending
    /// from-level (the only feasible pass order); same-level ties give 0 via R.
    double triple_joint(int i, int j, int k) const {
        int a = i, b = j, c = k;
        if (ts[b].from > ts[a].from) std::swap(a, b);
        if (ts[c].from > ts[b].from) std::swap(b, c);
        if (ts[b].from > ts[a].from) std::swap(a, b);
        return pass(a) * R[b][ts[a].to] * R[c][ts[b].to];
    }

    /// P(all transitions in `subset` passed on one downward walk), the N-subset
    /// generalization of pair_joint / triple_joint (Eq. 10' path enumeration).
    /// Sort by descending from-level (the only feasible pass order); any two
    /// leaving the same level give 0 via R (mutual exclusivity). Empty => 1;
    /// single => pass(); mirrors triple_joint for size 3. `subset` is copied so
    /// the caller's vector is untouched.
    double n_subset_joint(std::vector<int> subset) const {
        if (subset.empty()) return 1.0;
        std::sort(subset.begin(), subset.end(),
                  [this](int x, int y) { return ts[x].from > ts[y].from; });
        double p = pass(subset[0]);
        for (std::size_t k = 1; k < subset.size() && p > 0.0; ++k)
            p *= R[subset[k]][ts[subset[k - 1]].to];
        return p;
    }

    /// Reach-weighted EC-feed vacancy factor for transition i:
    /// Σ_L feed[L]·feed_ec{K,L}[L]·R[i][L]  — the (un-normalized) joint that the
    /// EC vacancy at the fed level is coincident with transition i being passed.
    /// Divide by pass(i) for the conditional probability. `ls` supplies the
    /// per-level EC yields (not stored on the DAG).
    double up_reach_ec(int i, const LevelScheme& ls, bool isL) const {
        double s = 0.0;
        for (int L = 0; L < n; ++L) {
            const double y = isL ? ls.levels[L].feed_ecL : ls.levels[L].feed_ecK;
            if (y != 0.0) s += feed[L] * y * R[i][L];
        }
        return s;
    }

    // ---- Survival-factor DP for analytic summing-OUT (Eq. 10') -------------
    //
    // Given per-transition survival factors f_t[i] = P(transition i contributes
    // no deposit to the primary window) and per-level EC survival factors
    // f_EC[l] = P(the EC-feed vacancy at level l contributes no deposit), the
    // summing-out factor for a primary transition p (A->B, emitting its gamma)
    // is  C_out(p) = Uw[A]·Dw[B] / V[A]  — the primary's own b·p_gamma cancels
    // against the "primary emitted" conditioning, and its factor f_p is excluded
    // (it emits the peak, not a partner). Dw/Uw depend only on f_t/f_EC (not on
    // which peak is primary), so compute them once per branch and reuse.
    //
    // Templated on the survival-factor scalar T (double, or an AD dual number
    // carrying d(eps)/d(param) — see AnalyticCascade.h); the DAG's own
    // probabilities (b, feed, V) stay double.

    /// Dw[l] = Σ over sub-paths l -> ground of Π(b·f_t): survival of everything
    /// strictly below level l. Terminal levels (no out-transitions, e.g. ground)
    /// = 1. `f_t` is indexed by transition.
    template <class T>
    std::vector<T> survival_down(const std::vector<T>& f_t) const {
        std::vector<T> Dw(n, T(1.0));
        for (int l = 1; l < n; ++l) {   // low -> high: children first
            T s(0.0);
            bool any = false;
            for (int i = 0; i < static_cast<int>(ts.size()); ++i)
                if (ts[i].from == l) { any = true; s += ts[i].b * f_t[i] * Dw[ts[i].to]; }
            if (any) Dw[l] = s;         // else terminal -> keep 1
        }
        return Dw;
    }

    /// Uw[l] = Σ over ways to arrive at level l from a fed level of
    /// feed·f_EC·Π(b·f_t): survival of everything strictly above level l,
    /// including the EC-feed vacancy at whatever level fed the path.
    template <class T>
    std::vector<T> survival_up(const std::vector<T>& f_t,
                               const std::vector<T>& f_EC) const {
        std::vector<T> Uw(n, T(0.0));
        for (int l = n - 1; l >= 0; --l) {  // high -> low: parents first
            T s = feed[l] * f_EC[l];
            for (int i = 0; i < static_cast<int>(ts.size()); ++i)
                if (ts[i].to == l) s += Uw[ts[i].from] * ts[i].b * f_t[i];
            Uw[l] = s;
        }
        return Uw;
    }

    /// C_out for primary transition `p` from precomputed Dw/Uw (see survival_up /
    /// survival_down). Returns 1.0 (no summing-out) if the primary level is
    /// unreachable.
    template <class T>
    T survival_out(int p, const std::vector<T>& Dw,
                   const std::vector<T>& Uw) const {
        const int A = ts[p].from, B = ts[p].to;
        if (V[A] <= 0.0) return T(1.0);
        return Uw[A] * Dw[B] / V[A];
    }
};

}  // namespace ceelo
