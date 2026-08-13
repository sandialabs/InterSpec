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

#include "cascade/SandiaDecayCascade.h"
#include "cascade/LevelDag.h"
#include "cross_sections/CrossSectionData.h"

#include "SandiaDecay.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace ceelo {
namespace cascade_adapter {

namespace {

using SandiaDecay::Nuclide;
using SandiaDecay::Transition;
using SandiaDecay::RadParticle;
using SandiaDecay::ProductType;
using SandiaDecay::NuclideActivityPair;
using SandiaDecay::NuclideMixture;

constexpr double kAnnih511keV = 510.998950; // electron rest mass energy (keV)
constexpr double kInternalPairThresholdKeV = 2.0 * kAnnih511keV;

// ---------------------------------------------------------------------------
// Internal-conversion data integrity.
//
// Two distinct defects reach this code through RadParticle::icc_total, and both
// destroy the level scheme's flux balance if taken at face value.
//
// (1) E0 SENTINEL. An electric-monopole (0+ -> 0+) transition cannot emit a
//     single gamma, so its conversion coefficient is formally infinite and
//     ENSDF / GEANT4 PhotonEvaporation encode that literally -- alpha up to
//     3.7e22, with the gamma's relative intensity set to 1e-20. Forming
//     I_gamma*(1+alpha) then inflates one level's out-flow by ~1e22.
//     Checked against all of PhotonEvaporation 6.1.2: among transitions with
//     alpha >= 1e4, every one above 100 keV is E0 (416 of them) and every
//     non-E0 one is at or below 100 keV -- the largest alpha above 100 keV that
//     is NOT E0 is 6130 (Hg-193, 101.25 keV). So the energy-qualified test has
//     no false positives, which matters: ~150 low-energy high-multipole
//     transitions carry genuinely huge PHYSICAL coefficients (Tc-99m 2.17 keV
//     E3 alpha=1.4e10, Sn-121m 6.29 keV M4 alpha=8.7e10) that must be kept.
//     Authoritative evaluated E0 provenance is carried independently as
//     RadParticle::e0_verified. GEANT4 multipole code 1 is used only to recognize
//     a low-energy transition that ALSO carries the >=1e4 sentinel. Code 1 by
//     itself is advisory: primary ENSDF leaves the M/CC blank or ambiguous for
//     some observed photons that G4 nevertheless labels 1 (Ac-232 373.3 keV,
//     Rb-78 189.8 keV).
//
// (2) INFEASIBLE TOTAL INTENSITY. A given nuclear transition occurs at most
//     once per decay, so I_gamma*(1+alpha) <= 1. Ambiguous GEANT4 <-> SandiaDecay
//     level matches violate this badly -- U-235's 19.59 keV line pairs
//     I_gamma = 0.61 with alpha = 114.7, i.e. 70.6 occurrences per decay, which
//     alone accounts for that branch's feeding sum of 72.3. This is the same
//     invariant, for the same reason, as the clamp build_radioactive_emissions
//     already applies below; it was simply never applied on the level-scheme
//     path.
constexpr double kE0AlphaSentinel = 1.0e4;
constexpr double kE0MinEnergyKeV = 100.0;

/// Largest level-scheme feeding sum accepted for a branch. The physical bound is
/// exactly 1 (a decay enters its daughter's level scheme at most once); the
/// margin absorbs the few-percent inconsistency present even in well-matched
/// schemes (Ba-133 1.020, Am-241 1.010) where ENSDF gamma intensities and GEANT4
/// conversion coefficients come from different evaluations.
constexpr double kMaxBranchFeed = 1.05;

/// GEANT4 PhotonEvaporation multipolarity code for an electric monopole.
constexpr int kMultipoleE0 = 1;

// Keep CeeLo source-compatible with SandiaDecay releases predating the
// enrichment duplicate marker.  A new library exposes the field and the branch
// below is compiled in; an older host simply has no records it can suppress.
template <class T, class = void>
struct HasSuppressedDuplicate : std::false_type {};
template <class T>
struct HasSuppressedDuplicate<T, std::void_t<decltype(
    std::declval<const T&>().suppressed_duplicate)>> : std::true_type {};

template <class T, class = void>
struct HasLevelFeeds : std::false_type {};
template <class T>
struct HasLevelFeeds<T, std::void_t<decltype(
    std::declval<const T&>().level_feeds)>> : std::true_type {};

template <class T, class = void>
struct HasE0Verified : std::false_type {};
template <class T>
struct HasE0Verified<T, std::void_t<decltype(
    std::declval<const T&>().e0_verified)>> : std::true_type {};

template <class T>
bool is_suppressed_duplicate_impl(const T& p)
{
    if constexpr (HasSuppressedDuplicate<T>::value)
        return p.suppressed_duplicate;
    return false;
}

bool is_suppressed_duplicate(const RadParticle& p)
{
    return is_suppressed_duplicate_impl(p);
}

template <class T>
bool is_e0_verified_impl(const T& p)
{
    if constexpr (HasE0Verified<T>::value)
        return p.e0_verified;
    return false;
}

template <class T>
bool exact_feeds_validate_partial_graph_impl(const T& tr,
                                             const LevelScheme& candidate)
{
    if constexpr (!HasLevelFeeds<T>::value) {
        return false;
    } else {
        if (tr.level_feeds.empty()) return false;
        double feed_sum = 0.0;
        std::vector<unsigned char> seen(candidate.levels.size(), 0);
        for (const auto& f : tr.level_feeds) {
            const double p = static_cast<double>(f.probability);
            if (f.level < 0 || f.level >= static_cast<int>(candidate.levels.size())
                || !(p > 0.0 && p <= 1.0)
                || seen[static_cast<std::size_t>(f.level)])
                return false;
            seen[static_cast<std::size_t>(f.level)] = 1;
            feed_sum += p;
        }
        if (std::abs(feed_sum - 1.0) > 2e-4) return false;

        // LevelDag deliberately ignores non-descending edges. That is a useful
        // defensive behavior for host-built graphs, but an enriched nuclear-
        // data graph must not become "compatible" merely because an invalid
        // candidate edge disappeared before replay.
        for (int level = 0; level < static_cast<int>(candidate.levels.size()); ++level)
            for (const LevelOutTransition& t : candidate.levels[
                     static_cast<std::size_t>(level)].out)
                if (t.to_level < 0 || t.to_level >= level)
                    return false;

        LevelScheme exact = candidate;
        exact.valid = true;
        const LevelDag dag(exact);
        if (!dag.valid) return false;
        std::vector<double> stored;
        for (int level = 0; level < static_cast<int>(exact.levels.size()); ++level) {
            double out_sum = 0.0;
            for (const LevelOutTransition& t : exact.levels[level].out)
                out_sum += t.weight;
            if (!(out_sum > 0.0)) continue;
            for (const LevelOutTransition& t : exact.levels[level].out)
                if (t.to_level >= 0 && t.to_level < level)
                    stored.push_back(t.weight);
        }
        if (stored.size() != dag.ts.size()) return false;

        for (int ti = 0; ti < static_cast<int>(dag.ts.size()); ++ti) {
            double predicted = 0.0;
            for (const auto& f : tr.level_feeds)
                predicted += static_cast<double>(f.probability)
                           * dag.pass_from_level(ti, f.level);
            const double observed = stored[static_cast<std::size_t>(ti)];
            if (std::max(predicted, observed) < 1e-5) continue;
            const double tolerance = std::max(
                5e-5, 0.05 * std::max(predicted, observed));
            if (std::abs(predicted - observed) > tolerance)
                return false;
        }
        return true;
    }
}

bool exact_feeds_validate_partial_graph(const Transition& tr,
                                        const LevelScheme& candidate)
{
    return exact_feeds_validate_partial_graph_impl(tr, candidate);
}

template <class T>
int max_exact_feed_level_impl(const T& tr)
{
    if constexpr (!HasLevelFeeds<T>::value) {
        return -1;
    } else {
        int result = -1;
        for (const auto& f : tr.level_feeds)
            result = std::max(result, f.level);
        return result;
    }
}

/// Identification of a definitive E0. Evaluated provenance takes precedence;
/// otherwise retain the conservative ICC-sentinel fallback. Above 100 keV the
/// >=1e4 sentinel signature is database-wide clean; below 100 keV code 1 is
/// additionally required so physical high-ICC E3/M4 transitions are preserved.
bool is_e0(const RadParticle& p)
{
    // Multipolarity code 1 is advisory unless the ICC also carries the E0
    // sentinel. GEANT4 labels several observed, finite-intensity photons as code
    // 1 despite ENSDF leaving their multipolarity/CC unresolved (Ac-232 373.3,
    // Rb-78 189.8); suppressing them on code alone discards measured lines.
    // The energy qualification preserves genuine low-energy, high-ICC physical
    // transitions, while code 1 identifies the rare low-energy sentinel E0.
    const double alpha = std::max(0.0, static_cast<double>(p.icc_total));
    const double energy_keV = static_cast<double>(p.energy);
    return is_e0_verified_impl(p)
        || (alpha >= kE0AlphaSentinel
            && (energy_keV > kE0MinEnergyKeV
                || p.multipole == kMultipoleE0));
}

/// The per-occurrence outcome split of one nuclear transition, with both data
/// defects above repaired. `weight` is the total transition intensity I_t; the
/// p_* are exclusive probabilities conditional on the transition occurring.
struct TransitionSplit {
    double weight = 0.0;
    double p_gamma = 1.0;
    double p_icK = 0.0, p_icL1 = 0.0, p_icL2 = 0.0, p_icL3 = 0.0;
    double p_icUnresolved = 0.0;
    double p_unmodeled = 0.0;
    bool e0_repaired = false;
    bool intensity_capped = false;
    bool split_renormalized = false;
};

/// The tabulated gamma intensity, bounded by what the transition can support.
///
/// A nuclear transition occurs at most once per decay, and emits its gamma on a
/// fraction 1/(1+alpha) of those occurrences, so
///     I_gamma  <=  1 / (1 + alpha)
/// is a hard physical ceiling -- true whichever of the two tabulated numbers is
/// the wrong one, so applying it needs no guess about the data's semantics.
///
/// This matters because the level-scheme path clamps the transition WEIGHT but
/// members are emitted at their own intensity, so a branch rejected by the flux
/// guard (exactly the branches whose data is known inconsistent) would otherwise
/// emit the offending line at its raw rate. U-235's 19.59 keV line pairs
/// I = 0.61 with alpha = 114.7 and was emitted 21x above the GEANT4/ENSDF scale,
/// feeding false sum peaks at 185.71+19.59 = 205.30 and 202.11+19.59 = 221.70
/// that land on real lines.
///
/// The bound degrades gracefully, which a semantic reinterpretation does not:
/// I-125's only Te-125 transition (I = 0.0668, alpha = 14.08, product 1.0073)
/// loses 0.7%, while U-235's 19.59 keV loses a factor of 70. That is the
/// distinction the raw product cannot make on its own.
double bounded_gamma_intensity(const RadParticle& p)
{
    const double inten = std::max(0.0, static_cast<double>(p.intensity));
    const double a = std::max(0.0, static_cast<double>(p.icc_total));
    return std::min(inten, 1.0 / (1.0 + a));
}

TransitionSplit transition_split(const RadParticle& p)
{
    const double inten = std::max(0.0, static_cast<double>(p.intensity));
    const double a = std::max(0.0, static_cast<double>(p.icc_total));
    const double aK = std::max(0.0, static_cast<double>(p.icc_K));
    const double aL1 = std::max(0.0, static_cast<double>(p.icc_L1));
    const double aL2 = std::max(0.0, static_cast<double>(p.icc_L2));
    const double aL3 = std::max(0.0, static_cast<double>(p.icc_L3));
    const double energy = static_cast<double>(p.energy);

    TransitionSplit s;
    if (is_e0(p)) {
        // No gamma, ever. The tabulated intensity is then the transition's
        // observed intensity rather than a gamma intensity, and alpha divides
        // out of the shell coefficients to leave the conversion-shell fractions
        // (GEANT4's icc_shell[] verbatim).
        //
        // Above the 1022-keV pair threshold an E0's unresolved remainder can be
        // INTERNAL PAIR FORMATION, which this engine does not model. Carry that
        // outcome explicitly so the level path advances without inventing a
        // photon, conversion electron, vacancy, or deposited energy.
        s.weight = inten;
        s.p_gamma = 0.0;
        if (a > 0.0) {
            s.p_icK = aK / a;
            s.p_icL1 = aL1 / a;
            s.p_icL2 = aL2 / a;
            s.p_icL3 = aL3 / a;
        }
        const double tracked = s.p_icK + s.p_icL1 + s.p_icL2 + s.p_icL3;
        if (energy < kInternalPairThresholdKeV) {
            s.p_icUnresolved = std::max(0.0, 1.0 - tracked);
            // This path is restricted to the definitive >=1e4 E0 sentinel. The
            // persisted K/L coefficients are shell fractions of that converted
            // transition; their below-pair remainder is an unpersisted shell.
        } else {
            s.p_unmodeled = std::max(0.0, 1.0 - tracked);
        }
        s.e0_repaired = true;
    } else {
        const double inv = 1.0 / (1.0 + a);
        s.weight = inten * (1.0 + a);
        s.p_gamma = inv;
        s.p_icK = aK * inv;
        s.p_icL1 = aL1 * inv;
        s.p_icL2 = aL2 * inv;
        s.p_icL3 = aL3 * inv;
    }

    // A transition occurs at most once per decay, so I_gamma*(1+alpha) <= 1.
    // 1052 transitions across the database violate this (1.6% of those carrying
    // level topology), by anything from 0.7% to a factor of 5e8.
    //
    // Clamping to the bound is deliberately the MINIMAL repair. The tempting
    // alternative -- reading the tabulated intensity as the TRANSITION's rate
    // rather than the gamma's, which the extreme cases invite (Ag-104m's 6.9 keV
    // isomeric transition lists intensity exactly 1.0 against alpha = 1e8) --
    // is catastrophic at the margin, and the margin is where most violations
    // are: 435 of the 1052 are within 2% of the bound. I-125's only Te-125
    // transition (I = 0.0668, alpha = 14.08 -> 1.0073) is a saturating
    // transition whose I_t genuinely is 1; reinterpreting it collapses the whole
    // branch 15x and puts the Te K x-rays 15x below GEANT4. Nothing but the size
    // of the violation distinguishes the two populations, and there is no
    // reference to validate a magnitude-based cut against, so the branch-level
    // flux guard rejects what the clamp cannot honestly repair.
    if (s.weight > 1.0) {
        s.weight = 1.0;
        s.intensity_capped = true;
    }

    // The outcomes are mutually exclusive, so they cannot sum past 1. Nothing
    // enforced this before; in the shipped data the violation was never material
    // (the worst pre-existing split is 1.00009, over 17 branches), but the E0
    // path above sets p_gamma = 0 and takes the shell fractions from a ratio, so
    // rounding there can push the sum over 1 and make later branches of the
    // engine's cumulative-u ladder unreachable.
    const double ic_sum = s.p_icK + s.p_icL1 + s.p_icL2 + s.p_icL3
                        + s.p_icUnresolved;
    // `p_unmodeled` already owns the above-pair E0 remainder. It is not room for
    // outer-shell conversion and must not be overwritten by a generic
    // `1-sum(tracked)` calculation in any consumer.
    const double ic_room = std::max(0.0, 1.0 - s.p_gamma - s.p_unmodeled);
    if (ic_sum > ic_room && ic_sum > 0.0) {
        const double f = ic_room / ic_sum;
        s.p_icK *= f;
        s.p_icL1 *= f;
        s.p_icL2 *= f;
        s.p_icL3 *= f;
        s.p_icUnresolved *= f;
        s.split_renormalized = true;
    }
    return s;
}

BetaShape beta_shape(SandiaDecay::ForbiddennessType forbiddenness)
{
    switch (forbiddenness) {
    case SandiaDecay::NoForbiddenness:
        return BetaShape::Allowed;
    case SandiaDecay::FirstForbidden:
        return BetaShape::FirstForbidden;
    case SandiaDecay::FirstUniqueForbidden:
        return BetaShape::FirstUniqueForbidden;
    case SandiaDecay::SecondForbidden:
        return BetaShape::SecondForbidden;
    case SandiaDecay::SecondUniqueForbidden:
        return BetaShape::SecondUniqueForbidden;
    case SandiaDecay::ThirdForbidden:
        return BetaShape::ThirdForbidden;
    case SandiaDecay::ThirdUniqueForbidden:
        return BetaShape::ThirdUniqueForbidden;
    case SandiaDecay::FourthForbidden:
        return BetaShape::FourthForbidden;
    }
    return BetaShape::Allowed;
}

template <class T>
bool append_exact_beta_level_feeds(const T& tr, double branch_ratio,
                                   int daughter_Z, int daughter_A,
                                   std::vector<BetaBranch>& out)
{
    if constexpr (!HasLevelFeeds<T>::value) {
        return false;
    } else {
        if (tr.mode != SandiaDecay::BetaDecay || tr.level_feeds.empty())
            return false;

        // SandiaDecay rejects these conditions while parsing. Keep the checks
        // here as well because callers may construct Transition objects by hand,
        // and silently falling back to the legacy aggregate beta record would
        // recreate the endpoint/shape error this exact feed law fixes.
        double sum = 0.0;
        std::vector<int> levels;
        levels.reserve(tr.level_feeds.size());
        for (const auto& feed : tr.level_feeds) {
            const double p = static_cast<double>(feed.probability);
            const double q = static_cast<double>(feed.q_keV);
            if (feed.level < 0 || !std::isfinite(p) || !(p > 0.0 && p <= 1.0)
                || !std::isfinite(q) || !(q > 0.0)
                || std::find(levels.begin(), levels.end(), feed.level) != levels.end())
                throw std::runtime_error("Malformed exact beta level-feed law");
            levels.push_back(feed.level);
            sum += p;
        }
        if (std::abs(sum - 1.0) > 2e-4)
            throw std::runtime_error(
                "Malformed exact beta level-feed law: probabilities do not sum to one");

        for (const auto& feed : tr.level_feeds) {
            out.push_back({static_cast<double>(feed.q_keV),
                           branch_ratio * static_cast<double>(feed.probability),
                           daughter_Z, daughter_A, false,
                           beta_shape(feed.forbiddenness)});
        }
        return true;
    }
}

WeakOutcomeLaw build_weak_outcome_law(const Transition& tr)
{
    WeakOutcomeLaw law;
    for (const RadParticle& p : tr.products) {
        if (is_suppressed_duplicate(p))
            continue;
        if (p.type != ProductType::CaptureElectronParticle
            && p.type != ProductType::PositronParticle)
            continue;
        const double mass = std::max(0.0, static_cast<double>(p.intensity));
        if (!(mass > 0.0))
            continue;

        WeakOutcome o;
        o.kind = (p.type == ProductType::CaptureElectronParticle)
            ? WeakOutcomeKind::ElectronCapture : WeakOutcomeKind::Positron;
        o.fed_level = p.fed_level;
        o.raw_mass = mass;
        if (o.kind == WeakOutcomeKind::ElectronCapture) {
            o.ec_K = std::max(0.0, static_cast<double>(p.ec_K));
            o.ec_L = std::max(0.0, static_cast<double>(p.ec_L));
            const double tracked = o.ec_K + o.ec_L;
            if (tracked > 1.0) {
                o.ec_K /= tracked;
                o.ec_L /= tracked;
            }
        } else {
            o.positron_endpoint_keV =
                std::max(0.0, static_cast<double>(p.energy));
        }
        law.raw_sum += mass;
        law.outcomes.push_back(o);
    }

    // A branch without EC or beta+ products has no weak categorical law. In
    // particular, do not retain every beta-/alpha branch merely by manufacturing
    // an all-Other law.
    if (law.outcomes.empty())
        return law;

    if (law.raw_sum > 1.0) {
        law.scale = 1.0 / law.raw_sum;
        law.confidence = (law.raw_sum <= kMaxBranchFeed)
            ? WeakOutcomeConfidence::CommonScaleWithinTolerance
            : WeakOutcomeConfidence::CommonScaleLowConfidence;
    }
    for (WeakOutcome& o : law.outcomes)
        o.selected_mass = o.raw_mass * law.scale;

    if (law.raw_sum < 1.0) {
        WeakOutcome other;
        other.kind = WeakOutcomeKind::Other;
        other.raw_mass = 1.0 - law.raw_sum;
        other.selected_mass = other.raw_mass;
        law.outcomes.push_back(other);
    }

    // Close the categorical law exactly despite floating-point roundoff. This
    // is a tiny correction to the last selected alternative, never an
    // independently sampled remainder.
    double selected = 0.0;
    for (const WeakOutcome& o : law.outcomes)
        selected += o.selected_mass;
    law.outcomes.back().selected_mass = std::max(
        0.0, law.outcomes.back().selected_mass + (1.0 - selected));
    return law;
}

void close_vacancy_group(VacancyGroup& g)
{
    g.p_K = std::max(0.0, g.p_K);
    g.p_L1 = std::max(0.0, g.p_L1);
    g.p_L2 = std::max(0.0, g.p_L2);
    g.p_L3 = std::max(0.0, g.p_L3);
    g.p_outer = std::max(0.0, g.p_outer);
    g.p_unmodeled = std::max(0.0, g.p_unmodeled);
    double active = g.p_K + g.p_L1 + g.p_L2 + g.p_L3 + g.p_outer
                  + g.p_unmodeled;
    if (active > 1.0) {
        const double f = 1.0 / active;
        g.p_K *= f;
        g.p_L1 *= f;
        g.p_L2 *= f;
        g.p_L3 *= f;
        g.p_outer *= f;
        g.p_unmodeled *= f;
        active = 1.0;
    }
    g.p_none = std::max(0.0, 1.0 - active);
}

void close_residual_transition(ResidualTransition& r)
{
    r.p_gamma = std::max(0.0, r.p_gamma);
    r.p_icK = std::max(0.0, r.p_icK);
    r.p_icL1 = std::max(0.0, r.p_icL1);
    r.p_icL2 = std::max(0.0, r.p_icL2);
    r.p_icL3 = std::max(0.0, r.p_icL3);
    r.p_icOuter = std::max(0.0, r.p_icOuter);
    r.p_icUnresolved = std::max(0.0, r.p_icUnresolved);
    r.p_unmodeled = std::max(0.0, r.p_unmodeled);
    double active = r.p_gamma + r.p_icK + r.p_icL1 + r.p_icL2
                  + r.p_icL3 + r.p_icOuter + r.p_icUnresolved
                  + r.p_unmodeled;
    if (active > 1.0) {
        const double f = 1.0 / active;
        r.p_gamma *= f;
        r.p_icK *= f;
        r.p_icL1 *= f;
        r.p_icL2 *= f;
        r.p_icL3 *= f;
        r.p_icOuter *= f;
        r.p_icUnresolved *= f;
        r.p_unmodeled *= f;
        active = 1.0;
    }
    r.p_none = std::max(0.0, 1.0 - active);
}

VacancyGroup make_e0_vacancy_group(const TransitionSplit& sp,
                                   double transition_energy_keV)
{
    VacancyGroup group;
    group.kind = VacancyGroupKind::ElectricMonopole;
    group.gamma_member = -1;
    group.transition_energy_keV = transition_energy_keV;
    group.p_K = sp.weight * sp.p_icK;
    group.p_L1 = sp.weight * sp.p_icL1;
    group.p_L2 = sp.weight * sp.p_icL2;
    group.p_L3 = sp.weight * sp.p_icL3;
    // VacancyGroup's compatibility field is named p_outer, but for a below-pair
    // E0 it means shell-UNRESOLVED conversion and produces only an approximate
    // electron. Above threshold the ambiguous remainder (which can be internal-
    // pair formation) is retained explicitly in p_unmodeled; p_none is only the
    // probability that this transition does not occur in the branch firing.
    group.p_outer = sp.weight * sp.p_icUnresolved;
    group.p_unmodeled = sp.weight * sp.p_unmodeled;
    close_vacancy_group(group);
    return group;
}

} // namespace

std::vector<DecayCascade> build_cascades(const Nuclide* parent,
                                         const CascadeOptions& opts)
{
    return build_cascades(parent, opts, nullptr);
}

std::vector<DecayCascade> build_cascades(
    const Nuclide* parent,
    const CascadeOptions& opts,
    std::vector<BranchFeedingAudit>* audit)
{
    std::vector<DecayCascade> out;
    if (!parent)
        return out;

    // Relative activity (decay rate) of each nuclide in the chain at the
    // requested age. Prompt equilibrium makes short-lived gamma-emitting
    // daughters present even at age 0 (e.g. Cs-137 -> Ba-137m).
    NuclideMixture mix;
    if (opts.prompt_equilibrium)
        mix.addNuclideInPromptEquilibrium(parent, 1.0);
    else
        mix.addNuclideByActivity(parent, 1.0);
    const std::vector<NuclideActivityPair> acts = mix.activity(opts.age_seconds);

    for (const NuclideActivityPair& na : acts) {
        const Nuclide* nuc = na.nuclide;
        const double act = static_cast<double>(na.activity);
        if (!nuc || act <= 0.0)
            continue;

        for (const Transition* tr : nuc->decaysToChildren) {
            if (!tr)
                continue;
            const double branch_weight = act * static_cast<double>(tr->branchRatio);
            if (branch_weight <= opts.min_branch_weight)
                continue;

            DecayCascade dc;
            dc.branch_weight = branch_weight;
            dc.weak_outcome_law = build_weak_outcome_law(*tr);

            BranchFeedingAudit ba;
            ba.parent = nuc->symbol;
            ba.child = tr->child ? tr->child->symbol : std::string("?");
            ba.branch_weight = branch_weight;
            for (const RadParticle& p : tr->products)
                if (is_suppressed_duplicate(p))
                    ++ba.n_suppressed_duplicates;

            // Decide whether to use the vacancy-level K x-ray model for this
            // branch: requires the daughter element and at least one K-vacancy
            // source (internal conversion of a gamma, or electron capture) with
            // the enriched ICC/EC data. Otherwise fall back to the transition-
            // group x-ray members.
            const int daughter_Z =
                tr->child ? static_cast<int>(tr->child->atomicNumber) : 0;
            bool any_vac = false;
            for (const RadParticle& p : tr->products) {
                if (is_suppressed_duplicate(p)) continue;
                if (p.intensity <= 0.0f) continue;
                if (p.type == ProductType::GammaParticle &&
                    (p.icc_K > 0.0f || p.icc_L1 > 0.0f || p.icc_L2 > 0.0f ||
                     p.icc_L3 > 0.0f)) any_vac = true;
                if (p.type == ProductType::CaptureElectronParticle &&
                    (p.ec_K > 0.0f || p.ec_L > 0.0f)) any_vac = true;
            }
            const bool use_vacancy =
                opts.vacancy_xray_model && daughter_Z > 0 && any_vac;

            // The vacancy model reproduces the daughter's K and L x-rays; pass 1
            // skips exactly those members (any others, e.g. M lines, are kept).
            const FluorescenceData* dfl =
                (use_vacancy ? CrossSectionData::instance().fluorescence(daughter_Z) : nullptr);
            const LFluorescenceData* dflL =
                (use_vacancy ? CrossSectionData::instance().l_fluorescence(daughter_Z) : nullptr);
            auto is_daughter_xray_line = [&](double e_keV) {
                if (dfl)
                    for (int li = 0; li < dfl->num_lines; ++li)
                        if (std::abs(e_keV - static_cast<double>(dfl->line_energy_keV[li])) < 1.5)
                            return true;
                if (dflL)
                    for (int s = 0; s < 3; ++s)
                        for (int li = 0; li < dflL->sub[s].num_lines; ++li)
                            if (std::abs(e_keV - static_cast<double>(dflL->sub[s].line_energy_keV[li])) < 1.0)
                                return true;
                return false;
            };

            // Map SandiaDecay product index -> index in dc.members for the kept
            // (photon) products; -1 for products we do not turn into members.
            const std::size_t n_prod = tr->products.size();
            std::vector<int> prod_to_member(n_prod, -1);
            std::vector<unsigned char> product_in_level_path(n_prod, 0);
            std::vector<std::size_t> gamma_members; // member indices of gammas

            // Pass 1: gammas and (optionally) x-rays become members. Skip
            // non-emitted (intensity <= 0) products and those below threshold.
            // When the vacancy model is active it supersedes the daughter's K
            // x-ray members (the model generates those); L / other x-rays are
            // kept as transition-group members.
            for (std::size_t i = 0; i < n_prod; ++i) {
                const RadParticle& p = tr->products[i];
                if (is_suppressed_duplicate(p))
                    continue;
                const bool is_gamma = (p.type == ProductType::GammaParticle);
                const bool is_xray  = (p.type == ProductType::XrayParticle);
                bool keep_xray = is_xray && opts.include_xrays;
                if (keep_xray && use_vacancy &&
                    is_daughter_xray_line(static_cast<double>(p.energy)))
                    keep_xray = false;  // daughter K/L x-ray -> generated by the vacancy model
                if (!(is_gamma || keep_xray))
                    continue;
                const double inten = static_cast<double>(p.intensity);
                if (inten <= 0.0 || inten < opts.min_intensity)
                    continue;
                // An E0 transition emits no photon at all, so it must not become
                // an emitting member. The level-path realization would suppress it
                // anyway (p_gamma = 0), but the pairwise/fallback path fires
                // members at their tabulated intensity and would invent a photon
                // the nucleus cannot produce -- Hg-182's 422 keV E0 at 64% per
                // decay, Pa-230's 634.9 keV at 29%. The transition is still
                // carried on the level-scheme path (with gamma_member = -1), and
                // the fallback path builds an unconditional vacancy source for it
                // below, so neither its flux nor its conversion products is lost.
                if (is_gamma && is_e0(p))
                    continue;

                CascadeMember m;
                m.energy_keV = static_cast<double>(p.energy);
                m.type = is_gamma ? CascadeParticleType::Gamma
                                  : CascadeParticleType::Xray;
                m.intensity = is_gamma ? bounded_gamma_intensity(p) : inten;
                prod_to_member[i] = static_cast<int>(dc.members.size());
                if (is_gamma)
                    gamma_members.push_back(dc.members.size());
                dc.members.push_back(std::move(m));
            }

            // Pass 2: copy gamma-gamma coincidences, remapping product indices to
            // member indices and dropping links to products we did not keep.
            for (std::size_t i = 0; i < n_prod; ++i) {
                const int mi = prod_to_member[i];
                if (mi < 0 || tr->products[i].type != ProductType::GammaParticle)
                    continue;
                for (const auto& c : tr->products[i].coincidences) {
                    const unsigned partner_prod = c.first;
                    const double frac = static_cast<double>(c.second);
                    if (partner_prod >= n_prod || frac <= 0.0)
                        continue;
                    const int partner_member = prod_to_member[partner_prod];
                    if (partner_member < 0)
                        continue;
                    // Carry the gamma-gamma angular-correlation coefficients from
                    // the enriched SandiaDecay data (0 => isotropic). They are
                    // only meaningful for an adjacent gamma pair; the enrichment
                    // tool leaves non-adjacent links at a2=a4=0.
                    CoincidenceLink link;
                    link.partner = static_cast<std::uint16_t>(partner_member);
                    link.prob = frac;
                    if (opts.angular_correlations &&
                        (c.a2 != 0.0f || c.a4 != 0.0f)) {
                        link.a2 = static_cast<double>(c.a2);
                        link.a4 = static_cast<double>(c.a4);
                        link.has_correlation = true;
                    }
                    dc.members[static_cast<std::size_t>(mi)].coincident.push_back(link);
                }
            }

            // Transition-group x-ray links (legacy model): each x-ray is
            // coincident with the branch's gammas. Used only when the vacancy
            // model is not active for this branch.
            if (opts.include_xrays && !use_vacancy) {
                for (CascadeMember& m : dc.members) {
                    if (m.type != CascadeParticleType::Xray)
                        continue;
                    for (std::size_t gi : gamma_members)
                        m.coincident.push_back({static_cast<std::uint16_t>(gi),
                                                dc.members[gi].intensity});
                }
            }

            // Vacancy-level K x-ray model. Prefer the LEVEL-PATH realization
            // (walk the daughter level scheme) when the full G4 topology is
            // present -- exact for multi-level cascades; otherwise fall back to
            // the independent per-gamma vacancy model.
            if (use_vacancy) {
                // Build the matched level graph and retain every unmatched nuclear
                // transition as an independent categorical residual.  Lack of
                // topology is not evidence that the matched graph is wrong (U-238
                // has a sound two-step graph plus several real high-energy lines);
                // Significant partial graphs additionally require independent,
                // complete RDM level feeds to reproduce their matched transition
                // occurrences; absent metadata fails closed to the old fallback.
                int max_lvl = 0;
                double matched_I = 0.0, total_I = 0.0;
                for (std::size_t i = 0; i < n_prod; ++i) {
                    const RadParticle& p = tr->products[i];
                    if (is_suppressed_duplicate(p))
                        continue;
                    if (p.type == ProductType::GammaParticle && prod_to_member[i] >= 0) {
                        total_I += static_cast<double>(p.intensity);
                        if (p.from_level >= 0 && p.to_level >= 0) {
                            matched_I += static_cast<double>(p.intensity);
                            max_lvl = std::max(max_lvl, std::max(p.from_level, p.to_level));
                        } else {
                            // Kept as an emitting member, but carries no topology, so
                            // its flux is absent from the level graph below.
                            ++ba.n_gammas_unmatched;
                            ba.unmatched_intensity += static_cast<double>(p.intensity);
                        }
                    }
                    // E0 transitions are not members (they emit no photon) but they
                    // still carry flux, so the level array must be sized to hold
                    // them or the build loop below indexes out of range.
                    if (p.type == ProductType::GammaParticle && prod_to_member[i] < 0
                        && p.from_level >= 0 && p.to_level >= 0
                        && is_e0(p))
                        max_lvl = std::max(max_lvl, std::max(p.from_level, p.to_level));
                    if (p.type == ProductType::CaptureElectronParticle && p.fed_level >= 0)
                        max_lvl = std::max(max_lvl, p.fed_level);

                    // Compute the physically relevant transition-occurrence
                    // coverage separately from the retained legacy raw gamma-
                    // intensity diagnostics above. Acceptance uses this mass.
                    // In particular, a memberless E0 can carry essentially all
                    // of a level's flux despite contributing no gamma intensity.
                    if (p.type == ProductType::GammaParticle) {
                        const TransitionSplit sp = transition_split(p);
                        const bool retained_member = prod_to_member[i] >= 0;
                        const bool retained_e0 = !retained_member && sp.e0_repaired
                            && sp.weight > 0.0 && sp.weight >= opts.min_intensity;
                        if (retained_member || retained_e0) {
                            ba.total_transition_occurrence_mass += sp.weight;
                            if (p.from_level >= 0 && p.to_level >= 0)
                                ba.matched_transition_occurrence_mass += sp.weight;
                        }
                    }
                }
                // Enrichment now supplies fed_level for positron outcomes too.
                // Size from the canonical weak law rather than only capture
                // products: a beta+ branch can terminate at a higher isomeric
                // level with no outgoing matched gamma, and that level is still
                // an exact categorical destination.
                for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
                    if (o.fed_level >= 0)
                        max_lvl = std::max(max_lvl, o.fed_level);
                max_lvl = std::max(max_lvl, max_exact_feed_level_impl(*tr));
                ba.matched_gamma_intensity = matched_I;
                ba.total_gamma_intensity = total_I;
                // Coverage is about nuclear-transition OCCURRENCES, not emitted
                // photons. Raw gamma intensity misses memberless E0 transitions
                // entirely and understates highly converted transitions. A graph
                // below 99% occurrence coverage is usable only when exact RDM
                // starts replay every matched transition and the uncorrelated
                // residual occurrence mass is at most 1% per branch.
                const bool materially_partial =
                    ba.total_transition_occurrence_mass > 0.0
                    && ba.matched_transition_occurrence_mass
                         < 0.99 * ba.total_transition_occurrence_mass;
                const double residual_occurrence = std::max(
                    0.0, ba.total_transition_occurrence_mass
                       - ba.matched_transition_occurrence_mass);
                bool built_levelpath = false;
                if (max_lvl >= 1) {
                    LevelScheme ls;
                    ls.daughter_Z = daughter_Z;
                    ls.levels.resize(static_cast<std::size_t>(max_lvl) + 1);
                    // Out-transitions per level, branching weight I_t = I_g*(1+alpha)
                    // from SandiaDecay's intensity (trusted) and G4 topology/ICC.
                    // transition_split() repairs definitive E0 records and the
                    // infeasible-intensity defects BEFORE the weight is formed, so
                    // the flux-conservation pass below sees only physical weights
                    // (a repair applied afterwards would be inert -- feeding is
                    // derived from these weights).
                    for (std::size_t i = 0; i < n_prod; ++i) {
                        const RadParticle& p = tr->products[i];
                        if (is_suppressed_duplicate(p)) continue;
                        if (p.type != ProductType::GammaParticle) continue;
                        const int mi = prod_to_member[i];
                        if (p.from_level < 0 || p.to_level < 0) continue;
                        const TransitionSplit sp = transition_split(p);
                        // mi < 0 is legal only for an E0 (no photon member); any
                        // other memberless gamma was filtered out deliberately and
                        // must not silently enter the flux graph. Apply pass 1's
                        // intensity gate to it as well, so a zero- or
                        // below-threshold E0 is dropped on the same terms as any
                        // other transition instead of entering at zero weight.
                        if (mi < 0 && !(sp.e0_repaired && sp.weight > 0.0
                                        && sp.weight >= opts.min_intensity))
                            continue;
                        // Every LevelDag consumer assumes a strictly descending
                        // level order and would silently discard this edge. Fail
                        // the entire candidate closed instead of returning a
                        // valid scheme with different topology from the audit.
                        if (p.to_level >= p.from_level) {
                            ++ba.n_invalid_topology_edges;
                            continue;
                        }
                        LevelOutTransition t;
                        t.to_level = p.to_level;
                        t.gamma_member = mi;
                        t.gamma_keV = static_cast<double>(p.energy);
                        t.weight = sp.weight;
                        t.p_gamma = sp.p_gamma;
                        t.p_icK = sp.p_icK;
                        t.p_icL1 = sp.p_icL1;
                        t.p_icL2 = sp.p_icL2;
                        t.p_icL3 = sp.p_icL3;
                        t.p_icUnresolved = sp.p_icUnresolved;
                        t.p_unmodeled = sp.p_unmodeled;
                        ls.levels[static_cast<std::size_t>(p.from_level)].out.push_back(t);
                        product_in_level_path[i] = 1;
                        if (sp.e0_repaired) {
                            ++ba.n_e0_repaired;
                            ++ba.n_graph_e0_repaired;
                        }
                        if (sp.intensity_capped) {
                            ++ba.n_intensity_capped;
                            ++ba.n_graph_intensity_capped;
                        }
                        if (sp.split_renormalized) {
                            ++ba.n_split_renormalized;
                            ++ba.n_graph_split_renormalized;
                        }
                    }
                    // Direct feeding by conservation: feeding(L) = out_flow(L) -
                    // in_flow(L). EC vacancy attached per fed level.
                    std::vector<double> out_flow(ls.levels.size(), 0.0),
                                        in_flow(ls.levels.size(), 0.0);
                    for (std::size_t L = 0; L < ls.levels.size(); ++L)
                        for (const auto& t : ls.levels[L].out) {
                            out_flow[L] += t.weight;
                            if (t.to_level >= 0 && t.to_level < static_cast<int>(in_flow.size()))
                                in_flow[static_cast<std::size_t>(t.to_level)] += t.weight;
                        }
                    double total_feed = 0.0, total_sink = 0.0;
                    int n_sinks = 0;
                    for (std::size_t L = 0; L < ls.levels.size(); ++L) {
                        const double net = out_flow[L] - in_flow[L];
                        ls.levels[L].feeding = std::max(0.0, net);
                        total_feed += ls.levels[L].feeding;
                        if (net < -1e-12) { ++n_sinks; total_sink += -net; }
                    }
                    const double raw_total_feed = total_feed;
                    ba.raw_total_feed = raw_total_feed;
                    if (raw_total_feed > 1.0
                        && raw_total_feed <= kMaxBranchFeed) {
                        // The graph is accepted, so do not leave it with an
                        // impossible absolute flux. Preserve every relative
                        // transition/feeding ratio with one common scale.
                        ba.feeding_scale = 1.0 / raw_total_feed;
                        for (CascadeLevel& level : ls.levels) {
                            level.feeding *= ba.feeding_scale;
                            for (LevelOutTransition& t : level.out)
                                t.weight *= ba.feeding_scale;
                        }
                        total_feed *= ba.feeding_scale;
                        total_sink *= ba.feeding_scale;
                    }
                    ba.total_feed = total_feed;
                    ba.total_sink = total_sink;
                    ba.n_sinks = n_sinks;
                    // feed_ec{K,L}[L] is P(the capture made a K/L vacancy | the decay
                    // landed on level L), so it is the capture-weighted shell fraction
                    // divided by the probability of landing there -- NOT ec_K itself.
                    // Assigning ec_K directly assumes every decay arriving at the level
                    // was a capture. That was harmless while only gamma-fed levels were
                    // reachable, but the ground state is also fed by betas and
                    // positrons, and it became reachable once the walk stopped entering
                    // the scheme unconditionally. Rb-82 -> Kr-82 (95% beta+, 5% EC)
                    // would otherwise report 0.76 K vacancies per decay against a true
                    // 0.040. Accumulating (rather than assigning) also fixes several
                    // capture branches feeding one level, where the last one used to win.
                    std::vector<double> ecK_num(ls.levels.size(), 0.0),
                                        ecL_num(ls.levels.size(), 0.0);
                    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
                        if (o.kind == WeakOutcomeKind::ElectronCapture
                            && o.fed_level >= 0
                            && o.fed_level < static_cast<int>(ls.levels.size())) {
                            const std::size_t li =
                                static_cast<std::size_t>(o.fed_level);
                            ecK_num[li] += o.selected_mass * o.ec_K;
                            ecL_num[li] += o.selected_mass * o.ec_L;
                        }
                    for (std::size_t L = 0; L < ls.levels.size(); ++L) {
                        // The ground state is where every decay that did not enter the
                        // scheme lands, so its landing probability is 1 - total_feed.
                        const double land = (L == 0) ? std::max(0.0, 1.0 - total_feed)
                                                     : ls.levels[L].feeding;
                        if (land <= 0.0) continue;
                        if (ecK_num[L] > 0.0)
                            ls.levels[L].feed_ecK = std::min(1.0, ecK_num[L] / land);
                        if (ecL_num[L] > 0.0)
                            ls.levels[L].feed_ecL = std::min(1.0, ecL_num[L] / land);
                    }
                    // A decay enters the level scheme at most once, so a RAW feeding sum
                    // above 1 means the graph is still carrying flux the decay data
                    // cannot support (a transition wired to the wrong levels, whose
                    // true branching we have no way to recover). Fall back to the
                    // independent-vacancy model rather than return a graph that
                    // looks usable and is not.
                    if (raw_total_feed > kMaxBranchFeed)
                        ba.rejected = true;
                    if (materially_partial) {
                        ba.partial_exact_feed_compatible =
                            exact_feeds_validate_partial_graph(*tr, ls);
                    }
                    // Exact starts validate only the matched graph. They contain
                    // no correlation information for topology-free residuals,
                    // so bound that independent approximation absolutely even
                    // when a many-step graph makes its relative coverage exceed
                    // 99%.
                    ba.partial_scheme_rejected = residual_occurrence > 0.01
                        || (materially_partial
                            && !ba.partial_exact_feed_compatible);
                    if (total_feed > 0.0 && !ba.rejected
                        && ba.n_invalid_topology_edges == 0
                        && !ba.partial_scheme_rejected) {
                        ls.entry_probability = std::min(1.0, total_feed);
                        ls.valid = true;
                        dc.level_scheme = std::move(ls);
                        dc.daughter_Z = daughter_Z;
                        built_levelpath = true;
                    }
                }
                ba.scheme_valid = built_levelpath;

                if (built_levelpath) {
                    // Unmatched Sandia gamma products still represent real nuclear
                    // transitions.  Sample each once as an absolute categorical
                    // law, independently of the matched graph because no topology
                    // exists with which to infer a stronger correlation.
                    for (std::size_t i = 0; i < n_prod; ++i) {
                        const RadParticle& p = tr->products[i];
                        if (is_suppressed_duplicate(p)
                            || p.type != ProductType::GammaParticle
                            || product_in_level_path[i])
                            continue;
                        const TransitionSplit sp = transition_split(p);
                        const int mi = prod_to_member[i];
                        if (mi < 0 && !(sp.e0_repaired && sp.weight > 0.0
                                        && sp.weight >= opts.min_intensity))
                            continue;
                        if (mi >= 0 && !(p.intensity > 0.0f
                                        && p.intensity >= opts.min_intensity))
                            continue;

                        ResidualTransition r;
                        r.gamma_member = mi;
                        r.transition_energy_keV = static_cast<double>(p.energy);
                        r.p_gamma = sp.weight * sp.p_gamma;
                        r.p_icK = sp.weight * sp.p_icK;
                        r.p_icL1 = sp.weight * sp.p_icL1;
                        r.p_icL2 = sp.weight * sp.p_icL2;
                        r.p_icL3 = sp.weight * sp.p_icL3;
                        r.p_icUnresolved = sp.weight * sp.p_icUnresolved;
                        r.p_unmodeled = sp.weight * sp.p_unmodeled;
                        const double conditional_used = sp.p_gamma + sp.p_icK
                            + sp.p_icL1 + sp.p_icL2 + sp.p_icL3
                            + sp.p_icUnresolved + sp.p_unmodeled;
                        r.p_icOuter = sp.weight
                            * std::max(0.0, 1.0 - conditional_used);
                        close_residual_transition(r);
                        if (r.p_none < 1.0) {
                            dc.residual_transitions.push_back(r);
                            ++ba.n_residual_transitions;
                            ba.residual_transition_occurrence_mass += 1.0 - r.p_none;
                            ba.residual_gamma_probability += r.p_gamma;
                        }
                        if (sp.e0_repaired) {
                            ++ba.n_e0_repaired;
                            ++ba.n_residual_e0_repaired;
                        }
                        if (sp.intensity_capped) {
                            ++ba.n_intensity_capped;
                            ++ba.n_residual_intensity_capped;
                        }
                        if (sp.split_renormalized) {
                            ++ba.n_split_renormalized;
                            ++ba.n_residual_split_renormalized;
                        }
                    }
                }

                if (!built_levelpath) {
                    // Fallback: one categorical conversion group per nuclear
                    // transition. Its shell alternatives are mutually exclusive;
                    // the adapter never emits independent K/L vacancy Bernoullis.
                    // EC alternatives already live in weak_outcome_law and are
                    // likewise selected categorically with positron emission.
                    dc.daughter_Z = daughter_Z;
                    for (std::size_t i = 0; i < n_prod; ++i) {
                        const RadParticle& p = tr->products[i];
                        if (is_suppressed_duplicate(p)) continue;
                        const double inten = static_cast<double>(p.intensity);
                        if (inten <= 0.0) continue;
                        if (p.type == ProductType::GammaParticle) {
                            const int mi = prod_to_member[i];
                            const TransitionSplit sp = transition_split(p);
                            VacancyGroup group;
                            group.transition_energy_keV =
                                static_cast<double>(p.energy);
                            if (sp.e0_repaired) {
                                dc.vacancy_groups.push_back(
                                    make_e0_vacancy_group(
                                        sp, static_cast<double>(p.energy)));
                                continue;
                            }
                            if (mi < 0) continue;
                            group.kind = VacancyGroupKind::InternalConversion;
                            group.gamma_member = mi;
                            const double p_gamma_abs =
                                std::max(0.0, dc.members[static_cast<std::size_t>(mi)].intensity);
                            const double gate_room = std::max(0.0, 1.0 - p_gamma_abs);
                            if (!(gate_room > 0.0))
                                continue;
                            group.p_K = sp.weight * sp.p_icK / gate_room;
                            group.p_L1 = sp.weight * sp.p_icL1 / gate_room;
                            group.p_L2 = sp.weight * sp.p_icL2 / gate_room;
                            group.p_L3 = sp.weight * sp.p_icL3 / gate_room;
                            const double tracked = sp.p_gamma + sp.p_icK
                                + sp.p_icL1 + sp.p_icL2 + sp.p_icL3;
                            group.p_outer = sp.weight
                                * std::max(0.0, 1.0 - tracked) / gate_room;
                            close_vacancy_group(group);
                            if (group.p_none < 1.0)
                                dc.vacancy_groups.push_back(group);
                        }
                    }
                }
            }

            // A definitive E0 can have no persisted K/L coefficients, so it does
            // not by itself activate use_vacancy. Still retain a memberless
            // categorical group (and its transition energy) for the explicitly
            // unmodelled internal-pair/remainder channel.
            if (opts.vacancy_xray_model && !use_vacancy) {
                dc.daughter_Z = daughter_Z;
                for (const RadParticle& p : tr->products) {
                    if (is_suppressed_duplicate(p))
                        continue;
                    if (p.type != ProductType::GammaParticle
                        || !(p.intensity > 0.0f))
                        continue;
                    const TransitionSplit sp = transition_split(p);
                    if (!sp.e0_repaired)
                        continue;
                    dc.vacancy_groups.push_back(make_e0_vacancy_group(
                        sp, static_cast<double>(p.energy)));
                }
            }

            // beta+ : one 511 keV annihilation member per positron intensity,
            // coincident with the branch gammas. The engine emits the back-to-
            // back pair from the annihilation point (Phase 3).
            if (opts.include_annihilation) {
                double positron_intensity = 0.0;
                double endpoint_weighted = 0.0;  // sum(selected mass * endpoint)
                for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
                    if (o.kind == WeakOutcomeKind::Positron) {
                        positron_intensity += o.selected_mass;
                        endpoint_weighted += o.selected_mass *
                            o.positron_endpoint_keV;
                    }
                if (positron_intensity > opts.min_intensity) {
                    CascadeMember m;
                    m.energy_keV = kAnnih511keV;
                    m.type = CascadeParticleType::Annih511;
                    m.intensity = positron_intensity;
                    // Intensity-weighted beta+ endpoint (SandiaDecay positron
                    // energy is the spectrum endpoint), for engine positron ranging.
                    m.positron_endpoint_keV = endpoint_weighted / positron_intensity;
                    for (std::size_t gi : gamma_members)
                        m.coincident.push_back({static_cast<std::uint16_t>(gi),
                                                dc.members[gi].intensity});
                    dc.members.push_back(std::move(m));
                }
            }

            if (audit)
                audit->push_back(ba);

            if (!dc.members.empty() || !dc.vacancy_groups.empty()
                || dc.weak_outcome_law.usable() || dc.level_scheme.valid)
                out.push_back(std::move(dc));
        }
    }

    return out;
}

std::vector<DecayCascade> build_cascades(const SandiaDecay::SandiaDecayDataBase& db,
                                         const std::string& nuclide_label,
                                         const CascadeOptions& opts)
{
    return build_cascades(db.nuclide(nuclide_label), opts, nullptr);
}

std::vector<DecayCascade> build_cascades(
    const SandiaDecay::SandiaDecayDataBase& db,
    const std::string& nuclide_label,
    const CascadeOptions& opts,
    std::vector<BranchFeedingAudit>* audit)
{
    return build_cascades(db.nuclide(nuclide_label), opts, audit);
}

RadioactiveEmissionSet build_radioactive_emissions(const Nuclide* parent)
{
    RadioactiveEmissionSet result;
    if (!parent)
        return result;
    result.nuclide = parent->symbol;

    for (const Transition* transition : parent->decaysToChildren) {
        if (!transition || !(transition->branchRatio > 0.0f))
            continue;
        const double branch_ratio =
            static_cast<double>(transition->branchRatio);
        const WeakOutcomeLaw weak_law = build_weak_outcome_law(*transition);
        const int daughter_Z = transition->child
            ? static_cast<int>(transition->child->atomicNumber) : 0;
        const int daughter_A = transition->child
            ? static_cast<int>(transition->child->massNumber) : 0;
        const FluorescenceData* k_data =
            daughter_Z >= 1 && daughter_Z <= kMaxZ
            ? CrossSectionData::instance().fluorescence(daughter_Z) : nullptr;
        const LFluorescenceData* l_data =
            daughter_Z >= 1 && daughter_Z <= kMaxZ
            ? CrossSectionData::instance().l_fluorescence(daughter_Z) : nullptr;
        const double k_binding =
            k_data ? static_cast<double>(k_data->k_edge_keV) : 0.0;
        const double l_binding =
            l_data ? static_cast<double>(l_data->l3_edge_keV) : 0.0;

        // A complete RDM beta-minus feed law carries one endpoint and shape per
        // daughter level. It supersedes SandiaDecay's legacy aggregate beta
        // products, whose lossy merge can duplicate endpoints and assign the
        // wrong intensity/forbiddenness to individual branches. Alpha level
        // feeds carry total decay Q rather than a beta endpoint and are ignored;
        // positron products retain the weak-outcome normalization above.
        const bool used_exact_beta_feeds = append_exact_beta_level_feeds(
            *transition, branch_ratio, daughter_Z, daughter_A,
            result.beta_branches);

        // Legacy SandiaDecay BetaParticle products are a lossy merge of level
        // feeds.  In 238 current fallback transitions their positive intensities
        // sum above one (occasionally by orders of magnitude), which would emit
        // more than one beta from a single nuclear transition.  Preserve their
        // endpoint/shape ratios but project the vector onto the physical
        // one-beta bound. Exact RDM level feeds above are authoritative and are
        // never touched by this fallback repair.
        double fallback_beta_sum = 0.0;
        if (!used_exact_beta_feeds) {
            for (const RadParticle& particle : transition->products)
                if (!is_suppressed_duplicate(particle)
                    && particle.type == ProductType::BetaParticle
                    && particle.intensity > 0.0f)
                    fallback_beta_sum +=
                        static_cast<double>(particle.intensity);
        }
        const double fallback_beta_scale = fallback_beta_sum > 1.0
            ? 1.0 / fallback_beta_sum : 1.0;

        for (const RadParticle& particle : transition->products) {
            if (is_suppressed_duplicate(particle))
                continue;
            const double particle_scale =
                particle.type == ProductType::PositronParticle
                    ? weak_law.scale
                    : particle.type == ProductType::BetaParticle
                    ? fallback_beta_scale : 1.0;
            const double intensity = branch_ratio * particle_scale *
                static_cast<double>(particle.intensity);
            if (!(intensity > 0.0))
                continue;

            if (particle.type == ProductType::BetaParticle) {
                if (used_exact_beta_feeds)
                    continue;
                result.beta_branches.push_back(
                    {static_cast<double>(particle.energy), intensity,
                     daughter_Z, daughter_A, false,
                     beta_shape(particle.forbiddenness)});
                continue;
            }
            if (particle.type == ProductType::PositronParticle) {
                result.beta_branches.push_back(
                    {static_cast<double>(particle.energy), intensity,
                     daughter_Z, daughter_A, true,
                     beta_shape(particle.forbiddenness)});
                continue;
            }
            if (particle.type != ProductType::GammaParticle)
                continue;

            const double gamma_energy = static_cast<double>(particle.energy);
            const TransitionSplit sp = transition_split(particle);
            const double gamma_yield = branch_ratio * sp.weight * sp.p_gamma;
            if (gamma_yield > 0.0)
                result.photons.push_back({gamma_energy, gamma_yield});

            // Reuse the same repaired occurrence weight and exclusive outcome
            // split as build_cascades(). This keeps the direct-source API from
            // bypassing the gamma ceiling on rejected branches (notably U-235's
            // 19.59-keV record).
            auto append_conversion =
                [&](double shell_probability, double binding,
                    ConversionShell shell) {
                const double electron_energy = gamma_energy - binding;
                const double electron_yield =
                    branch_ratio * sp.weight * shell_probability;
                if (electron_energy > 0.0 && electron_yield > 0.0) {
                    result.conversion_electrons.push_back(
                        {electron_energy, electron_yield, shell});
                }
            };
            append_conversion(sp.p_icK, k_binding, ConversionShell::K);
            append_conversion(sp.p_icL1, l_binding, ConversionShell::L1);
            append_conversion(sp.p_icL2, l_binding, ConversionShell::L2);
            append_conversion(sp.p_icL3, l_binding, ConversionShell::L3);
            // Storage compatibility: `Outer` here denotes a shell-unresolved
            // below-pair E0 release. Its KE is approximated by the transition
            // energy and no atomic vacancy is asserted.
            append_conversion(sp.p_icUnresolved, 0.0, ConversionShell::Outer);
            // `sp.p_unmodeled` is an explicit above-pair E0 outcome. Advancing
            // the nuclear level is its only represented effect; the direct API
            // intentionally emits no particle for it until pair formation is
            // modelled.
            // For non-E0 transitions the untracked split is outer-shell IC. The
            // below-pair E0 approximation was appended explicitly above; an E0
            // remainder above threshold may be internal-pair formation and stays
            // unmodelled rather than being mislabelled as a conversion electron.
            if (!sp.e0_repaired) {
                const double tracked = sp.p_gamma + sp.p_icK + sp.p_icL1
                    + sp.p_icL2 + sp.p_icL3;
                append_conversion(std::max(0.0, 1.0 - tracked), 0.0,
                                  ConversionShell::Outer);
            }
        }
    }
    // Split/duplicated parent transitions can remain mutually inconsistent even
    // after each fallback beta vector is physical. Do not silently common-scale
    // exact level-feed branches together with suspect legacy branches: without a
    // public provenance/audit channel there is no defensible choice of which
    // parent branch to alter. Fail closed instead. The tolerance covers only
    // float branch-ratio roundoff (~1e-7), not evaluated percent-level excesses.
    constexpr double kParentBetaYieldTolerance = 1e-6;
    const double total_beta_yield = result.beta_yield();
    if (total_beta_yield > 1.0 + kParentBetaYieldTolerance) {
        throw std::runtime_error(
            "Inconsistent beta/positron branch yield for " + result.nuclide
            + ": " + std::to_string(total_beta_yield));
    }
    if (total_beta_yield > 1.0) {
        // Numerical closure only (at most one part per million): unlike the
        // fail-closed path above, this does not conceal an evaluated branch
        // inconsistency and restores the strict one-beta physical invariant.
        const double scale = 1.0 / total_beta_yield;
        for (BetaBranch& branch : result.beta_branches)
            branch.yield_per_decay *= scale;
    }
    return result;
}

RadioactiveEmissionSet build_radioactive_emissions(
    const SandiaDecay::SandiaDecayDataBase& db,
    const std::string& nuclide_label)
{
    return build_radioactive_emissions(db.nuclide(nuclide_label));
}

} // namespace cascade_adapter
} // namespace ceelo
