// CascadeSummingProgram.h - coincidence-summing pieces of the analytic cascade engine, repackaged for a caller
// that evaluates them many thousands of times (inside a fit's cost function).
//
// compute_cascade_analytic (AnalyticCascade.h) is the reference: it takes an efficiency provider, looks every
// energy up through a hash memo, locates each window's lines, and enumerates the summing-in pairs with their
// coincident survival.  That is right for a one-off evaluation, and too slow per fit iteration.  Here:
//
//  - SummingOutProgram evaluates the SAME summing-out factor (the exact level-scheme survival expectation) for
//    a fixed list of lines, as a function of the total efficiency at a fixed list of energies.  Every lookup is
//    resolved at construction, the per-branch x-ray line sums are formed once, and one upward pass per branch
//    serves all of its primaries, so an evaluation is a few passes over the transitions.  It also evaluates
//    the engine's coincident survival of chosen photon pairs (the summing-out of a sum peak).
//  - coincidence_sum_pairs lists the photon pairs of one decay whose energies add into an energy range, with the
//    probability the summing-in enumeration gives them - the input a caller needs to model coincidence-sum
//    peaks, or summing-in to a line, itself.
//
// Both are templated on the efficiency scalar (double, or an automatic-differentiation dual number), like the
// engine.

#ifndef CEELO_CASCADE_SUMMING_PROGRAM_H
#define CEELO_CASCADE_SUMMING_PROGRAM_H

#include "cascade/AnalyticCascade.h"

#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

namespace ceelo {

/// A photon pair of one decay branch whose energies add into a requested range.
struct CoincidenceSumPair {
    int branch = -1;           ///< index into AnalyticCascadeModel::branches()
    double energy_a = 0.0;     ///< keV; energy_a <= energy_b
    double energy_b = 0.0;
    /// Pairs emitted per parent decay: branch_weight x P_branch(a and b) x the angular factor W(0) (if asked for).
    double weight = 0.0;
    /// `weight` without the angular factor: the probability, per parent decay, that both are emitted.
    double joint = 0.0;
    /// Emission probability per parent decay of each member's line alone, within its branch (branch_weight x
    /// P_branch(occurrence)); joint / emission_a is P(b | a).
    double emission_a = 0.0;
    double emission_b = 0.0;
    /// Whether each member is a nuclear gamma (or annihilation photon), rather than a relaxation x-ray.
    bool a_is_gamma = false;
    bool b_is_gamma = false;
    /// The members' level-scheme transitions (LevelDag::ts index), or -1.
    int transition_a = -1;
    int transition_b = -1;
};

/// Every pair of summing-in occurrences (gammas, the K/L x-ray lines of converting transitions and of
/// electron-capture vacancies, 511) of one decay branch with lo_keV <= E_a + E_b <= hi_keV and weight > 0, with
/// exactly the joint probability and (if `apply_angular_correlation`) the collinear angular factor
/// compute_cascade_analytic gives the pair.  Pairs are not merged.  The engine's coincident-survival factor
/// (no third emission deposits) is efficiency dependent, and so is not part of `weight`.
inline std::vector<CoincidenceSumPair> coincidence_sum_pairs(const AnalyticCascadeModel& model,
                                                             const double lo_keV, const double hi_keV,
                                                             const bool apply_angular_correlation = true) {
    std::vector<CoincidenceSumPair> pairs;
    const std::vector<detail::PreparedBranch>& branches = model.branches();
    for (int b = 0; b < static_cast<int>(branches.size()); ++b) {
        const detail::PreparedBranch& pb = branches[static_cast<std::size_t>(b)];
        const DecayCascade& dc = *pb.dc;
        const std::vector<detail::Occ>& occ = pb.occ;
        const std::vector<int>& ord = pb.ord;
        const std::vector<double>& es = pb.es;
        const std::size_t nocc = ord.size();

        auto wsp = [&](std::initializer_list<int> subset, const WeakOutcomeKind* kind = nullptr,
                       int shell = 0) -> double {
            return model.weak_subset_probability(b, subset, kind, shell);
        };

        // P_branch(occurrence alone): the joint's single-occurrence analogue.
        auto marginal = [&](const detail::Occ& o) -> double {
            if (o.residual >= 0) return o.factor;
            if (o.trans >= 0) return wsp({o.trans}) * o.factor;
            if (o.is_ec) {
                if (dc.weak_outcome_law.usable()) {
                    const WeakOutcomeKind ec = WeakOutcomeKind::ElectronCapture;
                    return wsp({}, &ec, o.ec_isL ? 2 : 1) * o.factor;
                }
                return (o.ec_isL ? pb.legacy_ec.total_l : pb.legacy_ec.total_k) * o.factor;
            }
            if (o.is_annih) {
                const double pm = dc.weak_outcome_law.usable()
                    ? detail::weak_kind_mass(dc, WeakOutcomeKind::Positron) : 1.0;
                return pm * o.factor;
            }
            return 0.0;
        };

        for (std::size_t ia = 0; ia < nocc; ++ia) {
            const detail::Occ& x = occ[static_cast<std::size_t>(ord[ia])];
            if (2.0 * es[ia] > hi_keV) break;  // sorted: every later pair is above the range
            const auto jb = std::lower_bound(es.begin() + static_cast<std::ptrdiff_t>(ia) + 1, es.end(),
                                             lo_keV - x.energy);
            for (auto it = jb; (it != es.end()) && (*it <= hi_keV - x.energy); ++it) {
                const detail::Occ& y = occ[static_cast<std::size_t>(ord[static_cast<std::size_t>(it - es.begin())])];
                const double J = detail::occurrence_pair_joint(pb, x, y, wsp);
                if (J <= 0.0) continue;
                double g = 1.0;
                if (apply_angular_correlation && x.member >= 0 && y.member >= 0 && !x.is_annih && !y.is_annih) {
                    const std::pair<double, double> a = detail::member_link_a2a4(dc, x.member, y.member);
                    if (a.first != 0.0 || a.second != 0.0) g = detail::angular_g_w0(a.first, a.second);
                }
                CoincidenceSumPair p;
                p.branch = b;
                p.energy_a = x.energy;
                p.energy_b = y.energy;
                p.weight = dc.branch_weight * J * g;
                p.joint = dc.branch_weight * J;
                p.emission_a = dc.branch_weight * marginal(x);
                p.emission_b = dc.branch_weight * marginal(y);
                p.a_is_gamma = (x.member >= 0);
                p.b_is_gamma = (y.member >= 0);
                p.transition_a = x.trans;
                p.transition_b = y.trans;
                pairs.push_back(p);
            }
        }
    }
    return pairs;
}


/// The summing-out factors of a fixed list of lines, as a function of the total (any-deposit) efficiency at a
/// fixed list of energies - compute_cascade_analytic's summing-out, with nothing left to look up.
///
/// A "line" is every emitted gamma member within `tolerance_keV` of `energy_keV`, over all branches; its
/// factor is the emission-weighted mean of its members' factors.  (The engine weights the lines of a window by
/// their full-energy efficiencies as well; the two agree for a window of one line energy.)
///
/// Immutable after construction, so one program may be evaluated from several threads, each with its own
/// Workspace.
class SummingOutProgram {
public:
    struct Line {
        double energy_keV = 0.0;
        double tolerance_keV = 0.0;
    };

    /// Two level-scheme transitions of one branch, whose photons (the transitions' gammas, or the x-rays of
    /// their internal conversion) sum: CoincidenceSumPair::branch / transition_a / transition_b, both >= 0.
    struct Pair {
        int branch = -1;
        int transition_a = -1;
        int transition_b = -1;
    };

    template <class T>
    struct Workspace {
        std::vector<T> f_t, f_ec, down, up, start, residual, transfer;
    };

    SummingOutProgram(std::shared_ptr<const AnalyticCascadeModel> model, const std::vector<Line>& lines,
                      const std::vector<Pair>& pairs = {})
        : model_(std::move(model)), lines_(lines), pairs_(pairs) {
        if (!model_) throw std::invalid_argument("SummingOutProgram: null model");
        build();
    }

    /// Sorted, distinct energies (keV) whose total efficiency evaluate() reads.
    const std::vector<double>& energies() const { return energies_; }
    std::size_t num_lines() const { return lines_.size(); }
    std::size_t num_pairs() const { return pairs_.size(); }
    const Line& line(std::size_t i) const { return lines_[i]; }
    /// Emission probability of line `i` per parent decay (0 if no member matched it; its factor is then 1).
    double line_emission(std::size_t i) const { return line_emission_[i]; }
    /// False when a branch emitting line `i` has no accepted level graph, so that share of the factor is the
    /// engine's flat pairwise form.
    bool line_complete(std::size_t i) const { return line_complete_[i] != 0; }

    /// out[i] = summing-out factor of line i, with tot_eff[k] the total efficiency at energies()[k].  If
    /// `pair_survival` is given, pair_survival[j] = the coincident survival of pair j: the probability that
    /// nothing else emitted with the pair deposits energy (the engine's pair_coincident_survival).
    template <class T>
    void evaluate(const T* tot_eff, T* out, Workspace<T>& ws, T* pair_survival = nullptr) const {
        using std::size_t;
        const size_t nlines = lines_.size();
        for (size_t i = 0; i < nlines; ++i) out[i] = T(0.0);
        if (pair_survival)
            for (size_t j = 0; j < pairs_.size(); ++j) pair_survival[j] = T(1.0);

        for (const Branch& br : branches_) {
            const detail::PreparedBranch& pb = model_->branches()[static_cast<size_t>(br.model_branch)];
            const DecayCascade& dc = *pb.dc;
            const LevelDag& dag = pb.dag;

            T sumK(0.0), sumL[3] = {T(0.0), T(0.0), T(0.0)}, dep_lmix(0.0);
            for (const std::pair<int, double>& ln : br.k_lines) sumK += ln.second * tot_eff[ln.first];
            for (int s = 0; s < 3; ++s) {
                for (const std::pair<int, double>& ln : br.l_lines[s]) sumL[s] += ln.second * tot_eff[ln.first];
                dep_lmix += detail::kLOcc[s] * sumL[s];
            }
            const T et511 = tot_eff[idx511_];

            if (br.valid) {
                const int nts = static_cast<int>(dag.ts.size());
                ws.f_t.assign(static_cast<size_t>(nts), T(1.0));
                for (int i = 0; i < nts; ++i) {
                    const LevelDag::Tr& tr = dag.ts[static_cast<size_t>(i)];
                    const int e = br.trans_energy[static_cast<size_t>(i)];
                    T dep = (e >= 0) ? T(tr.p_gamma * tot_eff[e]) : T(0.0);
                    if (tr.p_icK > 0.0) dep += tr.p_icK * sumK;
                    for (int s = 0; s < 3; ++s)
                        if (tr.p_icL[s] > 0.0) dep += tr.p_icL[s] * sumL[s];
                    ws.f_t[static_cast<size_t>(i)] = 1.0 - dep;
                }
                down_pass(dag, ws.f_t, ws.down);

                const size_t nres = dc.residual_transitions.size();
                ws.residual.assign(nres, T(1.0));
                for (size_t rix = 0; rix < nres; ++rix) {
                    const ResidualTransition& r = dc.residual_transitions[rix];
                    const int e = br.residual_energy[rix];
                    T dep = (e >= 0) ? T(r.p_gamma * tot_eff[e]) : T(0.0);
                    if (r.p_icK > 0.0) dep += r.p_icK * sumK;
                    if (r.p_icL1 > 0.0) dep += r.p_icL1 * sumL[0];
                    if (r.p_icL2 > 0.0) dep += r.p_icL2 * sumL[1];
                    if (r.p_icL3 > 0.0) dep += r.p_icL3 * sumL[2];
                    ws.residual[rix] = T(1.0) - dep;
                }

                if (br.weak) {
                    // Survival of everything above each level, started from the weak-outcome law: each
                    // outcome's own atomic-vacancy (or positron-pair) survival times its mass, at its fed level.
                    if (br.any_scheme_primary) {
                        T pos_survival = 1.0 - 2.0 * et511;
                        if (pos_survival < 0.0) pos_survival = T(0.0);
                        ws.start.assign(static_cast<size_t>(dag.n), T(0.0));
                        for (int l = 0; l < dag.n; ++l) {
                            const StartCoef& c = br.start[static_cast<size_t>(l)];
                            ws.start[static_cast<size_t>(l)] = c.plain - c.ec_k * sumK - c.ec_l * dep_lmix
                                                             + c.positron * pos_survival;
                        }
                        up_pass(dag, ws.f_t, ws.start, ws.up);
                    }
                } else {
                    ws.f_ec.assign(static_cast<size_t>(dag.n), T(1.0));
                    for (int l = 0; l < dag.n; ++l) {
                        const CascadeLevel& lv = dc.level_scheme.levels[static_cast<size_t>(l)];
                        T dep = (lv.feed_ecK > 0.0) ? T(lv.feed_ecK * sumK) : T(0.0);
                        for (int s = 0; s < 3; ++s) dep += lv.feed_ecL * detail::kLOcc[s] * sumL[s];
                        ws.f_ec[static_cast<size_t>(l)] = 1.0 - dep;
                    }
                    if (br.any_scheme_primary) {
                        // LevelDag::survival_up, into the reused vectors.
                        ws.start.assign(static_cast<size_t>(dag.n), T(0.0));
                        for (int l = 0; l < dag.n; ++l)
                            ws.start[static_cast<size_t>(l)] = dag.feed[static_cast<size_t>(l)] * ws.f_ec[static_cast<size_t>(l)];
                        up_pass(dag, ws.f_t, ws.start, ws.up);
                    }
                }
            }

            T pair511 = 2.0 * et511;  // P(one photon of a back-to-back pair deposits)
            if (pair511 > 1.0) pair511 = T(1.0);

            // The whole-scheme survival seen by a primary that is not a level-scheme transition (an unplaced
            // gamma, a 511): the same for every such primary of the branch, bar whether it is a 511 itself.
            T off_scheme(1.0), off_scheme_annih(1.0);
            if (br.valid && br.any_off_scheme_primary) {
                if (br.weak) {
                    // Survival below a start level drawn from the unknown-start weights.
                    T unknown_start(0.0);
                    double assigned = 0.0;
                    for (int l = 0; l < dag.n; ++l) {
                        const double w = pb.start_w[static_cast<size_t>(l)];
                        if (w == 0.0) continue;
                        unknown_start += w * ws.down[static_cast<size_t>(l)];
                        assigned += w;
                    }
                    unknown_start += std::max(0.0, 1.0 - assigned);

                    T pos_survival = T(1.0 - 2.0 * et511);
                    if (pos_survival < 0.0) pos_survival = T(0.0);
                    T num(0.0), num_annih(0.0);
                    double den = 0.0, den_annih = 0.0;
                    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                        const T so = (o.fed_level >= 0 && o.fed_level < dag.n)
                            ? ws.down[static_cast<size_t>(o.fed_level)] : unknown_start;
                        if (o.kind == WeakOutcomeKind::Positron) {
                            num += o.selected_mass * so * pos_survival;
                            // The detected 511's back-to-back partner points away.
                            num_annih += o.selected_mass * so;
                            den_annih += o.selected_mass;
                        } else if (o.kind == WeakOutcomeKind::ElectronCapture) {
                            num += o.selected_mass * so * (T(1.0) - o.ec_K * sumK - o.ec_L * dep_lmix);
                        } else {
                            num += o.selected_mass * so;
                        }
                        den += o.selected_mass;
                    }
                    if (den > 0.0) off_scheme = num / den;
                    if (den_annih > 0.0) off_scheme_annih = num_annih / den_annih;
                } else {
                    T s_full(0.0);
                    for (int l = 0; l < dag.n; ++l)
                        s_full += dag.feed[static_cast<size_t>(l)] * ws.f_ec[static_cast<size_t>(l)]
                                * ws.down[static_cast<size_t>(l)];
                    const double ep = dc.level_scheme.entry_probability;
                    off_scheme = ep * s_full + (1.0 - ep) * ws.f_ec[0];
                    off_scheme_annih = off_scheme;
                }
            }

            // Every residual transition's survival; a primary that is one divides its own back out.
            T all_residual(1.0);
            if (br.valid)
                for (const T& r : ws.residual) all_residual *= r;

            if (pair_survival && br.valid && !br.pairs.empty()) {
                const T& residual_survival = all_residual;

                // transfer[k*n + l]: probability of going from level l to the k'th target level (a level one
                // of the pairs' transitions leaves) with every transition on the way undetected.  One pass
                // per target level serves every pair that uses it.
                const size_t nl = static_cast<size_t>(dag.n);
                ws.transfer.resize(br.pair_targets.size() * nl);
                for (size_t k = 0; k < br.pair_targets.size(); ++k) {
                    const int target = br.pair_targets[k];
                    T* const to_target = &ws.transfer[k * nl];
                    for (int l = 0; l < dag.n; ++l) {
                        if (l <= target) {
                            to_target[l] = T((l == target) ? 1.0 : 0.0);
                            continue;
                        }
                        T s(0.0);
                        for (int j = dag.out_begin[static_cast<size_t>(l)];
                             j < dag.out_begin[static_cast<size_t>(l) + 1]; ++j) {
                            const LevelDag::Tr& tr = dag.ts[static_cast<size_t>(j)];
                            s += tr.b * ws.f_t[static_cast<size_t>(j)] * to_target[tr.to];
                        }
                        to_target[l] = s;
                    }
                }

                for (const PairProg& pp : br.pairs) {
                    // Survival of every transition but the pair's two, for paths through both, from each
                    // level: down to the upper transition's level, that transition, down to the lower one's
                    // level, it, and on to the bottom (LevelDag::constrained_survival_table, factored so the
                    // passes over the level scheme are shared between pairs).
                    const bool ok = (pp.hi >= 0);
                    const T* const to_hi = ok ? &ws.transfer[static_cast<size_t>(pp.hi_target) * nl] : nullptr;
                    T tail(0.0);
                    if (ok) {
                        const LevelDag::Tr& hi = dag.ts[static_cast<size_t>(pp.hi)];
                        const LevelDag::Tr& lo = dag.ts[static_cast<size_t>(pp.lo)];
                        tail = hi.b * ws.transfer[static_cast<size_t>(pp.lo_target) * nl + static_cast<size_t>(hi.to)]
                               * lo.b * ws.down[static_cast<size_t>(lo.to)];
                    }
                    auto qat = [&](const int level) -> T {
                        return (ok && level >= 0 && level < dag.n) ? T(to_hi[level] * tail) : T(0.0);
                    };

                    T surv(1.0);
                    if (br.weak) {
                        T num(0.0), unknown_start_q(0.0);
                        bool have_unknown_start = false;
                        double den = 0.0;
                        for (const PairOutcome& po : pp.outcomes) {
                            const WeakOutcome& o = dc.weak_outcome_law.outcomes[static_cast<size_t>(po.outcome)];
                            T joint_survival(0.0);
                            if (o.fed_level >= 0 && o.fed_level < dag.n) {
                                joint_survival = qat(o.fed_level);
                            } else {
                                if (!have_unknown_start) {  // the same for every outcome of unknown start
                                    for (int l = 0; l < dag.n; ++l)
                                        if (pb.start_w[static_cast<size_t>(l)] > 0.0)
                                            unknown_start_q += pb.start_w[static_cast<size_t>(l)] * qat(l);
                                    have_unknown_start = true;
                                }
                                joint_survival = unknown_start_q;
                            }
                            T so = joint_survival / po.base;
                            if (o.kind == WeakOutcomeKind::ElectronCapture)
                                so *= 1.0 - o.ec_K * sumK - o.ec_L * dep_lmix;
                            else if (o.kind == WeakOutcomeKind::Positron)
                                so *= 1.0 - pair511;
                            num += po.weight * so;
                            den += po.weight;
                        }
                        surv = residual_survival * ((den > 0.0) ? T(num / den) : T(1.0));
                    } else if (pp.base > 0.0) {
                        T joint_survival(0.0);
                        for (int l = 0; l < dag.n; ++l)
                            if (dag.feed[static_cast<size_t>(l)] > 0.0)
                                joint_survival += dag.feed[static_cast<size_t>(l)] * qat(l);
                        surv = joint_survival / pp.base;
                        if (pp.has_ec) surv *= (1.0 - pp.pec_k * sumK - pp.pec_l * dep_lmix);
                        surv = residual_survival * surv;
                    }
                    pair_survival[static_cast<size_t>(pp.index)] = surv;
                }
            }

            for (const Prim& p : br.prims) {
                T c(1.0);
                if (br.valid) {
                    if (p.t >= 0) {
                        const LevelDag::Tr& tr = dag.ts[static_cast<size_t>(p.t)];
                        if (br.weak) {
                            // (up x b) is the survival-weighted probability of arriving through the primary;
                            // p.reach is the same with every survival set to 1.
                            c = (p.reach > 0.0)
                                ? T(ws.up[static_cast<size_t>(tr.from)] * tr.b * ws.down[static_cast<size_t>(tr.to)]
                                    / p.reach)
                                : T(1.0);
                        } else {
                            c = dag.survival_out(p.t, ws.down, ws.up);
                        }
                    } else {
                        c = p.is_annih ? off_scheme_annih : off_scheme;
                    }
                    if (p.residual < 0) {
                        c *= all_residual;
                    } else {
                        for (size_t rix = 0; rix < ws.residual.size(); ++rix)
                            if (static_cast<int>(rix) != p.residual) c *= ws.residual[rix];
                    }
                }

                if (!br.weak) {
                    for (const std::pair<int, double>& ann : br.annih) {
                        if (ann.first == p.m) continue;  // the detected photon's partner points away
                        c *= (1.0 - ann.second * pair511);
                    }
                }

                if (!br.valid) {
                    for (const std::pair<int, double>& mate : p.flat_mates)
                        c *= (1.0 - mate.second * tot_eff[mate.first]);
                }

                out[static_cast<size_t>(p.line)] += p.share * c;
            }
        }

        for (size_t i = 0; i < nlines; ++i)
            if (!(line_emission_[i] > 0.0)) out[i] = T(1.0);
    }

    /// For each line, the first-order loss per unit total efficiency at each energy: with every total efficiency
    /// zero except eps at energies()[k], the line's factor is 1 - coef[line][k] * eps + O(eps^2).  A measure of
    /// which partners matter (and which energies a caller may drop).
    std::vector<std::vector<double>> first_order_coefficients() const {
        const std::size_t nE = energies_.size(), nL = lines_.size();
        std::vector<std::vector<double>> coef(nL, std::vector<double>(nE, 0.0));
        std::vector<double> tot(nE, 0.0), out(nL, 1.0);
        Workspace<double> ws;
        const double eps = 1.0e-6;
        for (std::size_t k = 0; k < nE; ++k) {
            tot[k] = eps;
            evaluate(tot.data(), out.data(), ws);
            tot[k] = 0.0;
            for (std::size_t i = 0; i < nL; ++i) coef[i][k] = (1.0 - out[i]) / eps;
        }
        return coef;
    }

private:
    /// Coefficients of a level's start weight in (1, sumK, dep_lmix, positron survival).
    struct StartCoef { double plain = 0.0, ec_k = 0.0, ec_l = 0.0, positron = 0.0; };

    struct Prim {
        int line = -1;
        int m = -1, t = -1;
        double share = 0.0;        ///< p_emit / (the line's total emission)
        double reach = 0.0;        ///< weak law, t >= 0: P_branch(the primary's transition is passed)
        int residual = -1;         ///< the primary's own residual transition, or -1
        bool is_annih = false;
        std::vector<std::pair<int, double>> flat_mates;  ///< no level graph: (energy index, P(mate | primary))
    };

    /// A weak outcome that can produce both transitions of a pair.
    struct PairOutcome {
        int outcome = -1;
        double base = 0.0;    ///< P(both transitions passed | the outcome)
        double weight = 0.0;  ///< the outcome's mass x base
    };

    struct PairProg {
        int index = -1;       ///< position in the caller's pair list
        int ta = -1, tb = -1;
        int hi = -1, lo = -1;               ///< ta and tb, the one leaving the higher level first; -1 if invalid
        int hi_target = -1, lo_target = -1; ///< positions in Branch::pair_targets of the levels they leave
        std::vector<PairOutcome> outcomes;  ///< weak law
        double base = 0.0;                  ///< no weak law: P(both transitions passed)
        bool has_ec = false;                ///< no weak law: the EC-feed vacancy above the pair
        double pec_k = 0.0, pec_l = 0.0;
    };

    struct Branch {
        int model_branch = -1;
        bool valid = false, weak = false, any_scheme_primary = false, any_off_scheme_primary = false;
        std::vector<PairProg> pairs;
        std::vector<int> pair_targets;      ///< the distinct levels the pairs' transitions leave
        std::vector<std::pair<int, double>> k_lines, l_lines[3];  ///< (energy index, weight)
        std::vector<int> trans_energy;      ///< per transition: energy index of its gamma, or -1
        std::vector<int> residual_energy;   ///< per residual transition, likewise
        std::vector<StartCoef> start;       ///< per level (weak law)
        std::vector<std::pair<int, double>> annih;  ///< legacy 511 members: (member, min(1, intensity))
        std::vector<Prim> prims;
    };

    /// LevelDag::survival_down, into a vector that is reused.
    template <class T>
    static void down_pass(const LevelDag& dag, const std::vector<T>& f_t, std::vector<T>& down) {
        down.assign(static_cast<std::size_t>(dag.n), T(1.0));
        for (int l = 1; l < dag.n; ++l) {  // low -> high: children first
            const int begin = dag.out_begin[static_cast<std::size_t>(l)];
            const int end = dag.out_begin[static_cast<std::size_t>(l) + 1];
            if (begin == end) continue;  // terminal: stays 1
            T s(0.0);
            for (int i = begin; i < end; ++i) {
                const LevelDag::Tr& tr = dag.ts[static_cast<std::size_t>(i)];
                s += tr.b * f_t[static_cast<std::size_t>(i)] * down[static_cast<std::size_t>(tr.to)];
            }
            down[static_cast<std::size_t>(l)] = s;
        }
    }

    /// up[l] = start[l] + sum over transitions i arriving at l of up[from_i] * b_i * f_t[i]: LevelDag::survival_up
    /// with an arbitrary start vector.
    template <class T>
    static void up_pass(const LevelDag& dag, const std::vector<T>& f_t, const std::vector<T>& start,
                        std::vector<T>& up) {
        up.assign(static_cast<std::size_t>(dag.n), T(0.0));
        for (int l = dag.n - 1; l >= 0; --l) {
            T s = start[static_cast<std::size_t>(l)];
            for (int k = dag.in_begin[static_cast<std::size_t>(l)]; k < dag.in_begin[static_cast<std::size_t>(l) + 1]; ++k) {
                const int i = dag.in_tr[static_cast<std::size_t>(k)];
                const LevelDag::Tr& tr = dag.ts[static_cast<std::size_t>(i)];
                s += up[static_cast<std::size_t>(tr.from)] * tr.b * f_t[static_cast<std::size_t>(i)];
            }
            up[static_cast<std::size_t>(l)] = s;
        }
    }

    int energy_index(const double energy) const {
        const auto pos = std::lower_bound(energies_.begin(), energies_.end(), energy);
        if (pos == energies_.end() || *pos != energy)
            throw std::logic_error("SummingOutProgram: energy not indexed");
        return static_cast<int>(pos - energies_.begin());
    }

    void build() {
        using std::size_t;
        const std::vector<detail::PreparedBranch>& all = model_->branches();
        line_emission_.assign(lines_.size(), 0.0);
        line_complete_.assign(lines_.size(), 1);

        // Pass 1: the primaries of each branch (as compute_cascade_analytic locates a window's lines).
        struct Found { int b, m, t, line; double p_emit; };
        std::vector<Found> found;
        for (int b = 0; b < static_cast<int>(all.size()); ++b) {
            const detail::PreparedBranch& pb = all[static_cast<size_t>(b)];
            const DecayCascade& dc = *pb.dc;
            const LevelDag& dag = pb.dag;
            for (int m = 0; m < static_cast<int>(dc.members.size()); ++m) {
                const CascadeMember& mem = dc.members[static_cast<size_t>(m)];
                int line = -1;
                for (size_t i = 0; i < lines_.size(); ++i) {
                    if (std::abs(mem.energy_keV - lines_[i].energy_keV) <= lines_[i].tolerance_keV) {
                        line = static_cast<int>(i);
                        break;
                    }
                }
                if (line < 0) continue;
                const int t = dag.valid ? dag.transition_of(m) : -1;
                double member_prob = mem.intensity;
                const int residual = (t < 0) ? pb.member_residual[static_cast<size_t>(m)] : -1;
                if (residual >= 0)
                    member_prob = dc.residual_transitions[static_cast<size_t>(residual)].p_gamma;
                else if (t < 0 && mem.type == CascadeParticleType::Annih511 && dc.weak_outcome_law.usable())
                    member_prob = 2.0 * detail::weak_kind_mass(dc, WeakOutcomeKind::Positron);
                else if (t < 0 && mem.type == CascadeParticleType::Annih511)
                    member_prob *= 2.0;
                const double p_emit = (t >= 0)
                    ? dc.branch_weight * model_->weak_subset_probability(b, {t}) * dag.ts[static_cast<size_t>(t)].p_gamma
                    : dc.branch_weight * member_prob;
                if (p_emit <= 0.0) continue;
                found.push_back({b, m, t, line, p_emit});
                line_emission_[static_cast<size_t>(line)] += p_emit;
                if (!dag.valid) line_complete_[static_cast<size_t>(line)] = 0;
            }
        }

        // Pass 2: every energy the branches with a primary (or a pair) can deposit at.
        std::vector<char> used(all.size(), 0);
        for (const Found& f : found) used[static_cast<size_t>(f.b)] = 1;
        for (const Pair& pr : pairs_) {
            if (pr.branch < 0 || static_cast<size_t>(pr.branch) >= all.size())
                throw std::invalid_argument("SummingOutProgram: pair branch out of range");
            const LevelDag& dag = all[static_cast<size_t>(pr.branch)].dag;
            const int nts = dag.valid ? static_cast<int>(dag.ts.size()) : 0;
            if (pr.transition_a < 0 || pr.transition_a >= nts || pr.transition_b < 0 || pr.transition_b >= nts)
                throw std::invalid_argument("SummingOutProgram: pair transition out of range");
            used[static_cast<size_t>(pr.branch)] = 1;
        }
        std::vector<double> energies;
        energies.push_back(detail::kAnnih511);
        for (size_t b = 0; b < all.size(); ++b) {
            if (!used[b]) continue;
            const detail::PreparedBranch& pb = all[b];
            for (const CascadeMember& mem : pb.dc->members) energies.push_back(mem.energy_keV);
            for (const detail::XrayLine& ln : pb.kl) energies.push_back(ln.energy);
            for (int s = 0; s < 3; ++s)
                for (const detail::XrayLine& ln : pb.ll[s]) energies.push_back(ln.energy);
        }
        std::sort(energies.begin(), energies.end());
        energies.erase(std::unique(energies.begin(), energies.end()), energies.end());
        energies_ = std::move(energies);
        idx511_ = energy_index(detail::kAnnih511);

        // Pass 3: the per-branch tables.
        for (int b = 0; b < static_cast<int>(all.size()); ++b) {
            if (!used[static_cast<size_t>(b)]) continue;
            const detail::PreparedBranch& pb = all[static_cast<size_t>(b)];
            const DecayCascade& dc = *pb.dc;
            const LevelDag& dag = pb.dag;
            const size_t nmem = dc.members.size();

            Branch br;
            br.model_branch = b;
            br.valid = dag.valid;
            br.weak = dc.weak_outcome_law.usable();
            for (const detail::XrayLine& ln : pb.kl) br.k_lines.emplace_back(energy_index(ln.energy), ln.weight);
            for (int s = 0; s < 3; ++s)
                for (const detail::XrayLine& ln : pb.ll[s])
                    br.l_lines[s].emplace_back(energy_index(ln.energy), ln.weight);

            auto member_energy_index = [&](const int member) -> int {
                if (member < 0 || static_cast<size_t>(member) >= nmem) return -1;
                return energy_index(dc.members[static_cast<size_t>(member)].energy_keV);
            };

            if (dag.valid) {
                br.trans_energy.reserve(dag.ts.size());
                for (const LevelDag::Tr& tr : dag.ts) br.trans_energy.push_back(member_energy_index(tr.gamma_member));
                for (const ResidualTransition& r : dc.residual_transitions)
                    br.residual_energy.push_back(member_energy_index(r.gamma_member));

                if (br.weak) {
                    br.start.assign(static_cast<size_t>(dag.n), StartCoef());
                    auto add = [&](const WeakOutcome& o, const int level, const double mass) {
                        StartCoef& c = br.start[static_cast<size_t>(level)];
                        if (o.kind == WeakOutcomeKind::ElectronCapture) {
                            c.plain += mass;
                            c.ec_k += mass * o.ec_K;
                            c.ec_l += mass * o.ec_L;
                        } else if (o.kind == WeakOutcomeKind::Positron) {
                            c.positron += mass;
                        } else {
                            c.plain += mass;
                        }
                    };
                    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                        if (o.fed_level >= 0 && o.fed_level < dag.n) {
                            add(o, o.fed_level, o.selected_mass);
                        } else {
                            for (int l = 0; l < dag.n; ++l)
                                if (pb.start_w[static_cast<size_t>(l)] > 0.0)
                                    add(o, l, o.selected_mass * pb.start_w[static_cast<size_t>(l)]);
                        }
                    }
                }
            }

            if (!br.weak) {
                for (int mm = 0; mm < static_cast<int>(nmem); ++mm)
                    if (dc.members[static_cast<size_t>(mm)].type == CascadeParticleType::Annih511)
                        br.annih.emplace_back(mm, std::min(1.0, dc.members[static_cast<size_t>(mm)].intensity));
            }

            for (const Found& f : found) {
                if (f.b != b) continue;
                Prim p;
                p.line = f.line;
                p.m = f.m;
                p.t = f.t;
                p.share = f.p_emit / line_emission_[static_cast<size_t>(f.line)];
                p.residual = pb.member_residual[static_cast<size_t>(f.m)];
                p.is_annih = (dc.members[static_cast<size_t>(f.m)].type == CascadeParticleType::Annih511);
                if (dag.valid && f.t >= 0) {
                    br.any_scheme_primary = true;
                    if (br.weak) p.reach = model_->weak_subset_probability(b, {f.t});
                } else if (dag.valid) {
                    br.any_off_scheme_primary = true;
                }
                if (!dag.valid) {
                    const CascadeMember& primary = dc.members[static_cast<size_t>(f.m)];
                    for (int mm = 0; mm < static_cast<int>(nmem); ++mm) {
                        const CascadeMember& mate = dc.members[static_cast<size_t>(mm)];
                        if (mm == f.m || mate.type == CascadeParticleType::Annih511) continue;
                        double pmate = mate.intensity;
                        for (const CoincidenceLink& lk : primary.coincident)
                            if (lk.partner == mm) { pmate = lk.prob; break; }
                        p.flat_mates.emplace_back(energy_index(mate.energy_keV), pmate);
                    }
                    if (primary.self_pair_prob > 0.0 && primary.intensity > 0.0)
                        p.flat_mates.emplace_back(energy_index(primary.energy_keV),
                                                  std::min(1.0, 2.0 * primary.self_pair_prob / primary.intensity));
                }
                br.prims.push_back(std::move(p));
            }

            for (size_t j = 0; j < pairs_.size(); ++j) {
                const Pair& pr = pairs_[j];
                if (pr.branch != b) continue;
                PairProg pp;
                pp.index = static_cast<int>(j);
                pp.ta = pr.transition_a;
                pp.tb = pr.transition_b;
                const int num_trans = static_cast<int>(dag.ts.size());
                if (pp.ta >= 0 && pp.ta < num_trans && pp.tb >= 0 && pp.tb < num_trans) {
                    const bool b_first = dag.ts[static_cast<size_t>(pp.tb)].from > dag.ts[static_cast<size_t>(pp.ta)].from;
                    pp.hi = b_first ? pp.tb : pp.ta;
                    pp.lo = b_first ? pp.ta : pp.tb;
                    auto target_of = [&br](const int level) -> int {
                        for (size_t k = 0; k < br.pair_targets.size(); ++k)
                            if (br.pair_targets[k] == level) return static_cast<int>(k);
                        br.pair_targets.push_back(level);
                        return static_cast<int>(br.pair_targets.size()) - 1;
                    };
                    pp.hi_target = target_of(dag.ts[static_cast<size_t>(pp.hi)].from);
                    pp.lo_target = target_of(dag.ts[static_cast<size_t>(pp.lo)].from);
                }
                const std::vector<int> defs{pp.ta, pp.tb};
                if (br.weak) {
                    const std::vector<WeakOutcome>& outcomes = dc.weak_outcome_law.outcomes;
                    for (int o = 0; o < static_cast<int>(outcomes.size()); ++o) {
                        const double base = detail::weak_outcome_reach_sw(dag, outcomes[static_cast<size_t>(o)],
                                                                          defs, pb.start_w);
                        const double weight = outcomes[static_cast<size_t>(o)].selected_mass * base;
                        if (weight > 0.0) pp.outcomes.push_back({o, base, weight});
                    }
                } else {
                    pp.base = dag.n_subset_joint(defs);
                    const int a_hi = (dag.ts[static_cast<size_t>(pp.tb)].from > dag.ts[static_cast<size_t>(pp.ta)].from)
                                         ? pp.tb : pp.ta;
                    const double pass_hi = dag.pass(a_hi);
                    if (pass_hi > 0.0) {
                        pp.has_ec = true;
                        pp.pec_k = dag.up_reach_ec(a_hi, dc.level_scheme, false) / pass_hi;
                        pp.pec_l = dag.up_reach_ec(a_hi, dc.level_scheme, true) / pass_hi;
                    }
                }
                br.pairs.push_back(std::move(pp));
            }
            branches_.push_back(std::move(br));
        }
    }

    std::shared_ptr<const AnalyticCascadeModel> model_;
    std::vector<Line> lines_;
    std::vector<Pair> pairs_;
    std::vector<double> energies_;
    int idx511_ = 0;
    std::vector<double> line_emission_;
    std::vector<char> line_complete_;
    std::vector<Branch> branches_;
};

}  // namespace ceelo

#endif  // CEELO_CASCADE_SUMMING_PROGRAM_H
