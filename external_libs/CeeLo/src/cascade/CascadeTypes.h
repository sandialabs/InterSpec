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

// Data-agnostic core types for true-coincidence (cascade) summing.
//
// These describe a nuclide's correlated photon emissions in a form the MC engine
// can consume WITHOUT any dependency on a nuclear-data library. A host
// application (InterSpec) or the optional SandiaDecay adapter
// (src/cascade/SandiaDecayCascade.h, compiled only with -DCEELO_WITH_SANDIADECAY)
// populates these structs; the engine never parses decay data itself.
//
// See the cascade-summing design plan for the full rationale (the gold-standard
// single-history correlated-cascade MC reuses the existing per-history deposit
// accumulator and scores against multiple PeakWindows at once).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace ceelo {

/// Kind of a cascade emission member. All members are transported as photons;
/// the type only tags how a member participates in summing bookkeeping (e.g.
/// 511 keV photons come as a back-to-back pair; x-ray summing dominates for
/// electron-capture emitters). It does not change transport physics.
enum class CascadeParticleType : uint8_t {
    Gamma,    ///< nuclear de-excitation gamma
    Xray,     ///< atomic-relaxation / decay x-ray (K or L line)
    Annih511  ///< 511 keV annihilation photon (one of a back-to-back pair)
};

/// A pairwise coincidence link from one cascade member to a partner within the
/// same DecayCascade. Carries the conditional emission probability and the
/// optional gamma-gamma angular-correlation coefficients.
///
/// Aggregate-initializable as `{partner, prob}` (a2/a4 default to isotropic),
/// which mirrors the old `std::pair<uint16_t,double>` call sites.
struct CoincidenceLink {
    /// Index of the partner member within DecayCascade::members.
    uint16_t partner = 0;

    /// P(partner emitted | this member emitted).
    double prob = 0.0;

    /// Gamma-gamma angular correlation W(theta) = 1 + a2*P2(cos t) + a4*P4(cos t),
    /// theta measured between this member's and the partner's emission directions.
    /// Only meaningful when has_correlation is true; absent => isotropic relative
    /// emission. W integrates to 1 over the sphere (P2,P4 average to zero), so it
    /// introduces no estimator weight.
    double a2 = 0.0;
    double a4 = 0.0;
    bool has_correlation = false;
};

/// One potential photon emission within a single decay branch, together with its
/// pairwise coincidence links to other members of the same branch.
struct CascadeMember {
    double energy_keV = 0.0;
    CascadeParticleType type = CascadeParticleType::Gamma;

    /// Probability this member is emitted per decay of its branch (its intensity
    /// within the branch, in [0,1]). For an annihilation photon it is the
    /// per-photon emission probability of the beta+ branch.
    double intensity = 1.0;

    /// For an Annih511 member only: the beta+ endpoint (maximum) kinetic energy
    /// (keV) of the positron, intensity-weighted across the branch's positron
    /// spectra. The engine samples the positron KE from the allowed beta+
    /// spectrum up to this endpoint to range the positron (annihilation site /
    /// escape) and to sample in-flight annihilation. 0 => the legacy
    /// point-annihilation-at-the-vertex behavior (back-compat / no positron info).
    double positron_endpoint_keV = 0.0;

    /// Pairwise coincidence links to other members of the same branch. Mirrors
    /// SandiaDecay's RadParticle::coincidences. Used by the realization sampler
    /// to reconstruct correlated co-emission from pairwise fractions, and (when
    /// has_correlation is set) to sample the angular correlation W(theta).
    std::vector<CoincidenceLink> coincident;
};

/// A per-decay K-shell vacancy production channel (electron capture, or internal
/// conversion of a nuclear transition). When it fires (with probability `prob`),
/// the daughter atom emits exactly one K x-ray — a single line sampled from the
/// element's fluorescence data (omega_K + line branching) — or an Auger electron
/// (no photon). Emitting one line per vacancy enforces K-line exclusivity (a
/// single vacancy makes Kalpha OR Kbeta, never both) and the correct per-decay
/// x-ray multiplicity, which the flat per-line x-ray members cannot represent —
/// see the vacancy-level x-ray model in the engine.
struct KVacancySource {
    double prob = 0.0;  ///< P(vacancy produced per decay of this branch)
    bool is_L = false;  ///< false => K-shell vacancy (K x-ray); true => L-shell (L x-ray)
    /// For an L-shell vacancy (is_L), which subshell: 0 => L1, 1 => L2, 2 => L3.
    /// The engine emits from that subshell's fluorescence line set. -1 for a K
    /// vacancy, or for an L vacancy of unresolved subshell (sample by occupancy).
    int l_subshell = -1;

    /// For an internal-conversion vacancy, the index (within DecayCascade::members)
    /// of the gamma whose transition this vacancy converts. A transition either
    /// emits its gamma OR converts, so the vacancy may fire only when that gamma
    /// was NOT emitted in the realization; `prob` is then the conditional
    /// P(IC-K | gamma not emitted). -1 for an electron-capture vacancy, which is
    /// unconditional.
    int gamma_member = -1;
};

/// The weak-decay outcome selected for one firing of a DecayCascade branch.
/// Electron capture and positron emission are competing alternatives, never
/// independent Bernoulli trials. `Other` is the residual probability when the
/// tabulated EC + beta+ mass is below one.
enum class WeakOutcomeKind : uint8_t {
    Other,
    ElectronCapture,
    Positron
};

/// How much repair was required to turn the tabulated weak-outcome masses into
/// one categorical probability law.
enum class WeakOutcomeConfidence : uint8_t {
    AsTabulated,                 ///< raw EC + beta+ mass was at most one
    CommonScaleWithinTolerance, ///< 1 < raw sum <= 1.05; all masses scaled alike
    CommonScaleLowConfidence    ///< raw sum > 1.05; scaled, but data are suspect
};

/// One alternative in a branch-conditional weak-decay categorical law.
///
/// There is one entry for every positive-intensity SandiaDecay EC or positron
/// RadParticle, preserving distinct fed levels and positron endpoints. The
/// synthetic `Other` entry has raw_mass == selected_mass and no shell/level
/// metadata. `selected_mass` values sum to one within a WeakOutcomeLaw; the
/// DecayCascade::branch_weight remains external and must not be folded in.
struct WeakOutcome {
    WeakOutcomeKind kind = WeakOutcomeKind::Other;
    int fed_level = -1;          ///< EC destination, including terminal isomers
    double raw_mass = 0.0;       ///< tabulated branch-conditional intensity
    double selected_mass = 0.0;  ///< normalized categorical probability
    double ec_K = 0.0;           ///< P(K vacancy | this EC outcome)
    double ec_L = 0.0;           ///< P(unresolved L vacancy | this EC outcome)
    double positron_endpoint_keV = 0.0;
};

struct WeakOutcomeLaw {
    std::vector<WeakOutcome> outcomes;
    double raw_sum = 0.0;  ///< sum of EC + positron raw masses (excludes Other)
    double scale = 1.0;    ///< common factor applied to non-Other raw masses
    WeakOutcomeConfidence confidence = WeakOutcomeConfidence::AsTabulated;

    bool usable() const noexcept { return !outcomes.empty(); }
};

/// Kind of a mutually-exclusive fallback conversion group.
enum class VacancyGroupKind : uint8_t {
    InternalConversion, ///< sample only when gamma_member was not emitted
    ElectricMonopole    ///< memberless E0; group probabilities are unconditional
};

/// One categorical fallback conversion/vacancy source.
///
/// Exactly one of none/K/L1/L2/L3/outer/unmodelled is selected. For ordinary
/// internal conversion these probabilities are conditional on the gated gamma
/// member not being emitted. For a memberless E0 they are unconditional per
/// branch firing. `transition_energy_keV` is retained so an IC-electron consumer
/// can subtract the selected shell binding energy and transport the electron.
struct VacancyGroup {
    VacancyGroupKind kind = VacancyGroupKind::InternalConversion;
    int gamma_member = -1;
    double transition_energy_keV = 0.0;
    double p_none = 1.0;
    double p_K = 0.0;
    double p_L1 = 0.0;
    double p_L2 = 0.0;
    double p_L3 = 0.0;
    double p_outer = 0.0;
    /// Above-pair E0 outcome retained for probability accounting only. Selecting
    /// it emits and deposits nothing until internal-pair formation is modelled.
    double p_unmodeled = 0.0;
};

/// One SandiaDecay nuclear transition which could not be placed on the matched
/// daughter-level DAG of an otherwise valid branch.
///
/// These are absolute, branch-conditional categorical probabilities: exactly
/// one of none/gamma/K/L1/L2/L3/outer/unresolved/unmodelled is selected on each
/// branch firing.  Keeping the occurrence probability here (rather than making
/// the unmatched photon an independent member Bernoulli) preserves gamma/IC
/// exclusivity and gives every estimator the same residual law.  Different
/// residual transitions are independent because no topology is available to
/// establish a stronger correlation.
struct ResidualTransition {
    int gamma_member = -1;       ///< -1 for a memberless E0 transition
    double transition_energy_keV = 0.0;
    double p_none = 1.0;
    double p_gamma = 0.0;
    double p_icK = 0.0;
    double p_icL1 = 0.0;
    double p_icL2 = 0.0;
    double p_icL3 = 0.0;
    double p_icOuter = 0.0;
    double p_icUnresolved = 0.0;
    double p_unmodeled = 0.0;
};

/// One de-excitation transition leaving a level in the daughter level scheme,
/// for the level-path cascade realization.
struct LevelOutTransition {
    int to_level = -1;        ///< index (into LevelScheme::levels) of the daughter level
    int gamma_member = -1;    ///< index (into DecayCascade::members) of the gamma
    double gamma_keV = 0.0;
    /// Branching weight at the parent level = I_gamma * (1 + alpha_total) (the
    /// total transition rate, including conversion). Taken from SandiaDecay's
    /// gamma intensity (the trusted branching), with the G4 topology.
    double weight = 0.0;
    /// Per-occurrence outcome split: emit the gamma, convert in K, or convert in
    /// one of the L subshells (L1/L2/L3 resolved separately so the correct
    /// subshell fluorescence is emitted).
    double p_gamma = 1.0;
    double p_icK = 0.0;
    double p_icL1 = 0.0;
    double p_icL2 = 0.0;
    double p_icL3 = 0.0;
    /// Unresolved non-photon energy release modeled as one conversion electron
    /// with KE approximately gamma_keV and no shell vacancy. Populated only for
    /// a definitive ICC-sentinel, memberless E0 below the exact internal-pair
    /// threshold; its persisted K/L shell fractions identify the remainder as
    /// conversion in an unpersisted shell, though the electron KE still uses the
    /// full transition energy as an approximation to the unknown binding.
    double p_icUnresolved = 0.0;
    /// Explicitly unmodelled E0 outcome. Populated only above the internal-pair
    /// threshold, where the remainder after persisted conversion shells may be
    /// internal-pair formation. The level path still advances, but no photon,
    /// conversion electron, atomic vacancy, or deposited energy is invented.
    double p_unmodeled = 0.0;
};

/// One level in the daughter level scheme.
struct CascadeLevel {
    /// Direct feeding weight into this level (electron capture or beta), relative
    /// to the other levels; the realization picks the initial level proportional
    /// to this.
    double feeding = 0.0;
    /// K/L-shell vacancy probability of the electron capture that feeds this
    /// level (0 for beta feeding or none).
    double feed_ecK = 0.0;
    double feed_ecL = 0.0;
    std::vector<LevelOutTransition> out;  ///< de-excitation transitions
};

/// Daughter level scheme for the level-path cascade realization. When valid, the
/// engine walks it per decay (pick the fed level, then at each level pick a
/// transition by weight and sample gamma-emit / internal-conversion), which gives
/// the correct per-decay gamma multiplicity, single de-excitation path, and IC
/// vacancy count that the independent-per-gamma vacancy model approximates.
struct LevelScheme {
    int daughter_Z = 0;
    std::vector<CascadeLevel> levels;  ///< index 0 = ground state
    bool valid = false;
    /// Probability that a decay of this branch enters the level scheme at all,
    /// i.e. 1 - P(the daughter's ground/isomeric state is fed directly). Equal to
    /// the sum of `feeding` over levels, which is why consumers must NOT recover
    /// it from the normalized `feed` a LevelDag builds: LevelDag divides by that
    /// sum, so its `pass()` is conditional on having entered the scheme. An
    /// absolute per-decay emission probability for a metadata-free cascade is
    /// therefore branch_weight * entry_probability * dag.pass(t) * p_gamma.
    /// When DecayCascade::weak_outcome_law is usable, its categorical outcome
    /// masses and exact fed levels replace that aggregate primitive; consumers
    /// sum selected_mass * dag.pass_from_level(t, fed_level) instead.
    /// Below 1 whenever some decays go straight to the ground state -- Th-234
    /// feeds the Pa-234m isomer directly in ~71% of decays, so its scheme carries
    /// entry_probability = 0.286. 1.0 when unset, so schemes built by a host
    /// application behave exactly as before.
    ///
    /// Do NOT round values just below 1 up to 1 on the theory that they are data
    /// noise. They are not, at least not always: Na-22's 1 - 0.99943796 =
    /// 0.056204% reproduces ENSDF's 0.056% beta+ branch direct to the Ne-22
    /// ground state to four digits. The feeding sum resolves real sub-percent
    /// ground-state branches, and for beta/positron/alpha feeding it is the only
    /// route to them: RadParticle::fed_level is populated for 0 of 20081 beta and
    /// 0 of 2620 alpha products (it IS populated for 12062 of 12783 electron
    /// capture products, 6860 of them at level 0 -- which is what makes the
    /// per-level EC vacancy yields below possible).
    ///
    /// Some otherwise accepted schemes have a raw feeding sum between 1 and
    /// 1.05 because ENSDF intensities and matched ICCs come from different
    /// evaluations. The adapter now common-scales every transition and feeding
    /// in those graphs so the selected feeding sum is exactly one while all
    /// relative paths are preserved. BranchFeedingAudit retains the raw sum and
    /// scale; consumers see only the physically bounded graph.
    ///
    /// Nor should very small values be floored. They are the correct answer for a
    /// branch that almost always feeds the ground state directly: Cs-137's 5.30%
    /// beta branch to Ba-137 carries 1.1557e-4, reproducing the published
    /// 5.8e-6/decay yield of its 283.5 keV line exactly. Ignoring it emitted that
    /// line 8653x too often -- a spurious photopeak at 12% of the 662 keV area --
    /// and Fe-55's 126 keV line 7.7e8x too often. The cost is statistical: a
    /// FullRealization observable fed only through such a branch needs ~1/ep more
    /// histories. The Conditional estimator is unaffected, since it reports
    /// efficiencies conditioned on the primary having been emitted.
    double entry_probability = 1.0;
};

/// One decay branch: a set of mutually-correlated potential emissions plus the
/// branch's selection weight. A nuclide (decayed to a given age) is a list of
/// these — see the SandiaDecay adapter.
struct DecayCascade {
    std::vector<CascadeMember> members;

    /// Branch selection weight = transition branchRatio x parent activity
    /// fraction, so that summing branch_weight over all of a nuclide's branches
    /// (across its decay chain at the requested age) gives the expected number of
    /// branch firings per parent decay.
    double branch_weight = 1.0;

    /// Daughter element atomic number, for looking up the K x-ray lines emitted
    /// by `k_vacancies`. 0 => unset.
    int daughter_Z = 0;

    /// K-shell vacancy sources for the vacancy-level x-ray model. When non-empty,
    /// the engine generates correlated K x-rays per decay from these (one line per
    /// fired vacancy) and the branch's flat Xray members are omitted by the
    /// adapter; when empty, the legacy transition-group Xray members are used.
    std::vector<KVacancySource> k_vacancies;

    /// Exclusive EC/beta+ alternatives for this branch. Empty when the decay
    /// branch has neither EC nor positron products. Consumers select exactly one
    /// outcome; in particular, an EC vacancy can never coincide with this
    /// branch's annihilation pair.
    WeakOutcomeLaw weak_outcome_law;

    /// Mutually-exclusive fallback IC/E0 vacancy groups. The SandiaDecay
    /// adapter populates these instead of `k_vacancies`; the latter remains only
    /// as a source-compatibility bridge for host-built DecayCascades.
    std::vector<VacancyGroup> vacancy_groups;

    /// Categorical transitions not owned by an otherwise valid level DAG.
    /// Empty for an invalid scheme: rejected branches retain the pre-existing
    /// member + VacancyGroup fallback without a second emission path.
    std::vector<ResidualTransition> residual_transitions;

    /// Daughter level scheme. When valid, the engine uses the level-path
    /// realization (walking this scheme) instead of the pairwise-coincidence
    /// member realization + independent k_vacancies — exact for multi-level
    /// cascades. Gammas are still in `members` (referenced by gamma_member, for
    /// transport + the W(theta) coincidence links).
    LevelScheme level_scheme;
};

/// Bounded member marginal used by the invalid-branch forest and to decide
/// whether an invalid branch can materially affect a requested sum peak. An Annih511 member represents a
/// positron outcome (and hence a back-to-back photon pair), so its probability
/// is the selected positron mass when an explicit weak-outcome law is present.
inline double cascade_fallback_photon_probability(const DecayCascade& dc,
                                                   std::size_t member) {
    if (member >= dc.members.size()) return 0.0;
    const CascadeMember& m = dc.members[member];
    if (m.type != CascadeParticleType::Annih511
        || !dc.weak_outcome_law.usable())
        return std::clamp(m.intensity, 0.0, 1.0);
    double p = 0.0;
    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
        if (o.kind == WeakOutcomeKind::Positron) p += o.selected_mass;
    return std::clamp(p, 0.0, 1.0);
}

/// One node in the coherent fallback realization used when no accepted level
/// graph is available. Nodes are stored in parent-before-child order and form a
/// forest (at most one parent each), so one Bernoulli draw per member defines a
/// proper multivariate law while preserving every selected marginal.
struct CascadeFallbackNode {
    std::size_t member = 0;
    int parent_member = -1;
    double marginal = 0.0;
    double p_if_parent_emitted = 0.0;
    double p_if_parent_absent = 0.0;
    bool forced_annihilation = false;
};

/// Project the legacy directional pair links onto a deterministic, coherent
/// one-parent Bayesian forest. Members are ordered by decreasing bounded
/// marginal (index breaks ties). For each child, every link to an earlier member
/// implies J = I_from*P(to|from); non-finite inputs are discarded and J is
/// projected onto the Fréchet interval [max(0,Ip+Ic-1), min(Ip,Ic)]. The parent
/// with the strongest |J-Ip*Ic| is retained (lower member index breaks ties).
/// When both link directions exist, the stronger projected implication wins;
/// an exact tie retains parent->child before child->parent. This cannot honor all
/// mutually inconsistent links, but unlike the former sequential probability
/// boost it is one valid joint law and preserves every member marginal.
inline std::vector<CascadeFallbackNode> build_cascade_fallback_forest(
    const DecayCascade& dc) {
    const std::size_t n = dc.members.size();
    std::vector<double> marginal(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double p = cascade_fallback_photon_probability(dc, i);
        marginal[i] = std::isfinite(p) ? std::clamp(p, 0.0, 1.0) : 0.0;
    }
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (marginal[a] != marginal[b]) return marginal[a] > marginal[b];
        return a < b;
    });

    std::vector<CascadeFallbackNode> forest;
    forest.reserve(n);
    std::vector<char> earlier(n, 0);
    for (std::size_t child : order) {
        CascadeFallbackNode node;
        node.member = child;
        node.marginal = marginal[child];
        node.p_if_parent_emitted = node.marginal;
        node.p_if_parent_absent = node.marginal;
        node.forced_annihilation = dc.weak_outcome_law.usable()
            && dc.members[child].type == CascadeParticleType::Annih511;

        // A weak-law annihilation state is forced by the categorical decay-mode
        // draw. It may parent later nodes, but selecting a parent for it would be
        // misleading because its own Bernoulli conditional is not sampled.
        if (!node.forced_annihilation) {
            double best_strength = -1.0;
            std::size_t best_parent = n;
            double best_joint = 0.0;
            for (std::size_t parent = 0; parent < n; ++parent) {
                if (!earlier[parent]) continue;
                const double ip = marginal[parent], ic = marginal[child];
                const double lo = std::max(0.0, ip + ic - 1.0);
                const double hi = std::min(ip, ic);
                bool have_pair = false;
                double pair_strength = -1.0, pair_joint = 0.0;
                auto consider = [&](double from_marginal, double probability) {
                    if (!std::isfinite(probability)) return;
                    const double raw_joint = from_marginal
                        * std::clamp(probability, 0.0, 1.0);
                    const double joint = std::clamp(raw_joint, lo, hi);
                    const double strength = std::abs(joint - ip * ic);
                    // Strict comparison intentionally retains the first exact
                    // tie: parent->child links are visited before reverse links.
                    if (!have_pair || strength > pair_strength) {
                        have_pair = true;
                        pair_strength = strength;
                        pair_joint = joint;
                    }
                };
                for (const CoincidenceLink& link : dc.members[parent].coincident)
                    if (link.partner == child)
                        consider(ip, link.prob);
                for (const CoincidenceLink& link : dc.members[child].coincident)
                    if (link.partner == parent)
                        consider(ic, link.prob);
                if (!have_pair) continue;
                if (pair_strength > best_strength
                    || (pair_strength == best_strength && parent < best_parent)) {
                    best_strength = pair_strength;
                    best_parent = parent;
                    best_joint = pair_joint;
                }
            }
            if (best_parent < n) {
                node.parent_member = static_cast<int>(best_parent);
                const double ip = marginal[best_parent];
                node.p_if_parent_emitted = ip > 0.0
                    ? best_joint / ip : node.marginal;
                node.p_if_parent_absent = ip < 1.0
                    ? (node.marginal - best_joint) / (1.0 - ip)
                    : node.marginal;
                node.p_if_parent_emitted = std::clamp(
                    node.p_if_parent_emitted, 0.0, 1.0);
                node.p_if_parent_absent = std::clamp(
                    node.p_if_parent_absent, 0.0, 1.0);
            }
        }
        forest.push_back(node);
        earlier[child] = 1;
    }
    return forest;
}

/// Draw one coherent fallback realization. `forced_annihilation_state` is -1
/// without an exact weak law, otherwise 0/1 from its already-selected outcome.
/// A uniform is consumed for every member, including a forced Annih511 member,
/// so the fallback keeps a fixed one-draw-per-member RNG contract.
template <class UniformDraw>
inline void sample_cascade_fallback_forest(
    const std::vector<CascadeFallbackNode>& forest,
    int forced_annihilation_state, UniformDraw&& uniform,
    std::vector<char>& emitted) {
    emitted.assign(forest.size(), 0);
    for (const CascadeFallbackNode& node : forest) {
        const double draw = uniform();
        if (node.member >= emitted.size()) continue;
        if (node.forced_annihilation && forced_annihilation_state >= 0) {
            emitted[node.member] = forced_annihilation_state ? 1 : 0;
            continue;
        }
        double p = node.marginal;
        if (node.parent_member >= 0
            && static_cast<std::size_t>(node.parent_member) < emitted.size())
            p = emitted[static_cast<std::size_t>(node.parent_member)]
                ? node.p_if_parent_emitted : node.p_if_parent_absent;
        emitted[node.member] = draw < p ? 1 : 0;
    }
}

/// Conservative diagnostic for a gamma/x-ray/annihilation photon pair in an
/// invalid branch which could feed `peak_keV`.  This intentionally does not
/// invent a joint probability from the legacy pairwise coincidence metadata:
/// the database contains directionally inconsistent and marginal-infeasible
/// links, so no single exact joint law can satisfy all of them.  The caller uses
/// this predicate to surface an incomplete/approximate result instead.
inline bool cascade_invalid_branch_can_feed_peak(
    const DecayCascade& dc, double peak_keV, double tolerance_keV,
    int excluded_member = -1, double material_probability = 1e-9) {
    const double branch_mass = std::max(0.0, dc.branch_weight);
    if (dc.level_scheme.valid || !(branch_mass > 0.0)) return false;
    for (std::size_t i = 0; i < dc.members.size(); ++i) {
        if (static_cast<int>(i) == excluded_member) continue;
        const double pi = cascade_fallback_photon_probability(dc, i);
        if (!(branch_mass * pi > material_probability)) continue;
        // One Annih511 member represents both back-to-back photons of a single
        // positron outcome.  A surrounding/well detector can therefore receive
        // their 2E sum even though there is no second flat member to enumerate.
        if (dc.members[i].type == CascadeParticleType::Annih511
            && std::abs(2.0 * dc.members[i].energy_keV - peak_keV)
                   <= tolerance_keV)
            return true;
        for (std::size_t j = i + 1; j < dc.members.size(); ++j) {
            if (static_cast<int>(j) == excluded_member) continue;
            const double pj = cascade_fallback_photon_probability(dc, j);
            // min(pi,pj) is the conservative upper bound on any physical pair
            // joint.  This can flag a numerically tiny independent joint, but it
            // cannot silently clear a materially possible correlated pair.
            if (!(branch_mass * std::min(pi, pj) > material_probability))
                continue;
            const double sum = dc.members[i].energy_keV
                             + dc.members[j].energy_keV;
            if (std::isfinite(sum)
                && std::abs(sum - peak_keV) <= tolerance_keV)
                return true;
        }
    }
    return false;
}

/// A photopeak for which a summing-corrected efficiency is reported.
struct PeakWindow {
    double energy_keV = 0.0;
    double tolerance_keV = 1.5;  ///< +/- window (matches the engine's kFepTolerance)
};

} // namespace ceelo
