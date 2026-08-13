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

#define BOOST_TEST_MODULE CascadeAdapterTests
#include <boost/test/unit_test.hpp>

#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace ceelo;
using namespace ceelo::cascade_adapter;

namespace {

// Load the decay database once (parses the ~31 MB sandia.decay.xml).
const SandiaDecay::SandiaDecayDataBase& db() {
    static SandiaDecay::SandiaDecayDataBase database(SANDIA_DECAY_XML_PATH);
    return database;
}

// First member within tol keV of E in a cascade, or nullptr.
const CascadeMember* find_member(const DecayCascade& dc, double E, double tol) {
    for (const auto& m : dc.members)
        if (std::abs(m.energy_keV - E) < tol)
            return &m;
    return nullptr;
}

} // namespace

BOOST_AUTO_TEST_CASE(database_loads) {
    BOOST_REQUIRE(db().initialized());
    BOOST_REQUIRE(db().nuclide("Co60") != nullptr);
}

BOOST_AUTO_TEST_CASE(co60_two_gamma_coincident) {
    // Co-60: 1173.2 + 1332.5 keV emitted in prompt coincidence (f ~ 1.0).
    const auto cascades = build_cascades(db(), "Co60");
    BOOST_REQUIRE(!cascades.empty());

    bool found_pair = false;
    for (const auto& dc : cascades) {
        const CascadeMember* g1 = find_member(dc, 1173.2, 2.0);
        const CascadeMember* g2 = find_member(dc, 1332.5, 2.0);
        if (!g1 || !g2)
            continue;
        found_pair = true;
        BOOST_REQUIRE(!g1->coincident.empty());

        bool links_to_1332 = false;
        for (const auto& c : g1->coincident) {
            const CascadeMember& partner = dc.members[c.partner];
            if (std::abs(partner.energy_keV - 1332.5) < 2.0) {
                links_to_1332 = true;
                BOOST_CHECK_GT(c.prob, 0.9);  // coincidence fraction ~ 1.0
                // Enriched data: the 1173-1332 (4(E2)2(E2)0) pair carries the
                // gamma-gamma angular correlation. G4's data gives the 1173 a
                // small E2/M3 admixture so a2 ~ 0.10 (pure-E2 literature 0.102).
                BOOST_CHECK(c.has_correlation);
                BOOST_CHECK_CLOSE(c.a2, 0.10, 10.0);   // within 10%
                BOOST_CHECK_CLOSE(c.a4, 0.0088, 25.0);
            }
        }
        BOOST_CHECK(links_to_1332);
    }
    BOOST_CHECK(found_pair);
}

BOOST_AUTO_TEST_CASE(cs137_single_gamma_no_gamma_coincidence) {
    // Cs-137 -> Ba-137m (prompt eq.) -> 661.7 keV single gamma: no coincident
    // gamma partner (it may carry x-ray links from the IC of that transition).
    const auto cascades = build_cascades(db(), "Cs137");
    BOOST_REQUIRE(!cascades.empty());

    const CascadeMember* g = nullptr;
    const DecayCascade* owner = nullptr;
    for (const auto& dc : cascades) {
        const CascadeMember* m = find_member(dc, 661.657, 1.0);
        if (m) { g = m; owner = &dc; }
    }
    BOOST_REQUIRE(g != nullptr);

    for (const auto& c : g->coincident) {
        const CascadeMember& partner = owner->members[c.partner];
        BOOST_CHECK(partner.type != CascadeParticleType::Gamma);
    }
}

BOOST_AUTO_TEST_CASE(ba133_has_gammas_and_vacancy_xray_model) {
    // Ba-133 (EC): rich gamma cascade + Cs K x-rays — a classic x-ray summer.
    // With the enriched ICC/EC data, the adapter uses the vacancy-level K x-ray
    // model: the Cs daughter Z + K-vacancy sources instead of flat Xray members.
    const auto cascades = build_cascades(db(), "Ba133");
    BOOST_REQUIRE(!cascades.empty());

    bool has_356 = false, has_gamma_coinc = false;
    bool has_levelpath = false, daughter_is_Cs = false;
    for (const auto& dc : cascades) {
        if (find_member(dc, 356.0, 2.0))
            has_356 = true;
        for (const auto& m : dc.members)
            if (m.type == CascadeParticleType::Gamma && !m.coincident.empty())
                has_gamma_coinc = true;
        // Ba-133 has the full G4 level topology -> the level-path realization.
        if (dc.level_scheme.valid) {
            has_levelpath = true;
            if (dc.daughter_Z == 55) daughter_is_Cs = true;  // Cs-133
            // A level must feed the cascade, and a level must produce a K vacancy.
            double feed = 0.0, ec = 0.0;
            for (const auto& lv : dc.level_scheme.levels) { feed += lv.feeding; ec += lv.feed_ecK; }
            BOOST_CHECK_GT(feed, 0.0);
            BOOST_CHECK_GT(ec, 0.0);
        }
    }
    BOOST_CHECK(has_356);
    BOOST_CHECK(has_gamma_coinc);
    BOOST_CHECK(has_levelpath);      // level-path vacancy x-ray model active
    BOOST_CHECK(daughter_is_Cs);     // Cs K x-rays
}

BOOST_AUTO_TEST_CASE(coincidence_partner_indices_in_range) {
    // Every coincidence partner index must be a valid member index, and the
    // fraction must be in (0, 1].
    for (const char* nuc : {"Co60", "Cs137", "Ba133", "Eu152", "Na22"}) {
        const auto cascades = build_cascades(db(), nuc);
        BOOST_REQUIRE_MESSAGE(!cascades.empty(), nuc);
        for (const auto& dc : cascades) {
            BOOST_CHECK_GT(dc.branch_weight, 0.0);
            for (const auto& m : dc.members) {
                for (const auto& c : m.coincident) {
                    BOOST_CHECK_LT(c.partner, dc.members.size());
                    BOOST_CHECK_GT(c.prob, 0.0);
                    BOOST_CHECK_LE(c.prob, 1.0 + 1e-6);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(na22_has_annihilation_member) {
    // Na-22 (beta+): a 511 keV annihilation member, coincident with the 1274 keV
    // gamma, must be produced.
    const auto cascades = build_cascades(db(), "Na22");
    BOOST_REQUIRE(!cascades.empty());

    bool has_511 = false, has_1274 = false;
    for (const auto& dc : cascades) {
        for (const auto& m : dc.members) {
            if (m.type == CascadeParticleType::Annih511)
                has_511 = true;
            if (m.type == CascadeParticleType::Gamma &&
                std::abs(m.energy_keV - 1274.5) < 3.0)
                has_1274 = true;
        }
    }
    BOOST_CHECK(has_511);
    BOOST_CHECK(has_1274);
}

BOOST_AUTO_TEST_CASE(uranium_chain_direct_electron_source_terms) {
    const auto u238 = build_radioactive_emissions(db(), "U238");
    const auto th234 = build_radioactive_emissions(db(), "Th234");
    const auto pa234m = build_radioactive_emissions(db(), "Pa234m");

    BOOST_CHECK_EQUAL(u238.nuclide, "U238");
    BOOST_CHECK_SMALL(u238.beta_yield(), 1e-14);
    BOOST_CHECK_GT(th234.beta_yield(), 0.9);
    BOOST_CHECK_GT(pa234m.beta_yield(), 0.9);

    double largest_pa_endpoint = 0.0;
    for (const auto& branch : pa234m.beta_branches)
        largest_pa_endpoint =
            std::max(largest_pa_endpoint, branch.endpoint_keV);
    BOOST_CHECK_GT(largest_pa_endpoint, 2000.0);

    for (const auto* emissions : {&u238, &th234, &pa234m}) {
        for (const auto& branch : emissions->beta_branches) {
            BOOST_CHECK_GT(branch.endpoint_keV, 0.0);
            BOOST_CHECK_GT(branch.yield_per_decay, 0.0);
            BOOST_CHECK_GT(branch.daughter_Z, 0);
            BOOST_CHECK_GT(branch.daughter_A, 0);
        }
        for (const auto& line : emissions->conversion_electrons) {
            BOOST_CHECK_GT(line.energy_keV, 0.0);
            BOOST_CHECK_GT(line.yield_per_decay, 0.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(exact_level_feeds_drive_direct_beta_branches) {
    struct Expected {
        const char* nuclide;
        std::size_t branches;
        double yield;
        int first_unique;
    };
    const Expected expected[] = {
        {"Th234",   4, 1.00000, 0},
        {"Pa234m", 23, 0.99840, 1},
        {"Pa234",  55, 1.00000, 3},
        {"Th231",  14, 1.00000, 1},
        {"Pb214",   6, 1.00000, 0},
        {"Bi214",  70, 0.99979, 1},
    };

    for (const Expected& x : expected) {
        const RadioactiveEmissionSet e =
            build_radioactive_emissions(db(), x.nuclide);
        BOOST_TEST_CONTEXT(x.nuclide) {
            BOOST_REQUIRE_EQUAL(e.beta_branches.size(), x.branches);
            double yield = 0.0;
            int first_unique = 0;
            std::set<double> endpoints;
            for (const BetaBranch& beta : e.beta_branches) {
                BOOST_CHECK(!beta.is_positron);
                BOOST_CHECK_GT(beta.endpoint_keV, 0.0);
                BOOST_CHECK_GT(beta.yield_per_decay, 0.0);
                endpoints.insert(beta.endpoint_keV);
                yield += beta.yield_per_decay;
                if (beta.shape == BetaShape::FirstUniqueForbidden)
                    ++first_unique;
            }
            // One exact RDM feed per daughter level: neither the lossy legacy
            // beta merge nor a second copy of it may survive alongside these.
            BOOST_CHECK_EQUAL(endpoints.size(), x.branches);
            BOOST_CHECK_SMALL(yield - x.yield, 5e-6);
            BOOST_CHECK_EQUAL(first_unique, x.first_unique);
        }
    }

    // SandiaDecay validates this on parse; the adapter also fails closed for a
    // host-constructed malformed exact law instead of reverting to the lossy
    // aggregate beta products.
    SandiaDecay::Transition* beta = nullptr;
    const SandiaDecay::Nuclide* th234 = db().nuclide("Th234");
    BOOST_REQUIRE(th234 != nullptr);
    for (const SandiaDecay::Transition* t : th234->decaysToChildren)
        if (t && t->mode == SandiaDecay::BetaDecay && !t->level_feeds.empty()) {
            beta = const_cast<SandiaDecay::Transition*>(t);
            break;
        }
    BOOST_REQUIRE(beta != nullptr);
    auto feeds = beta->level_feeds;
    beta->level_feeds.front().probability = 0.0f;
    BOOST_CHECK_THROW(build_radioactive_emissions(th234), std::runtime_error);
    beta->level_feeds = std::move(feeds);
}

BOOST_AUTO_TEST_CASE(fallback_direct_beta_vectors_are_bounded) {
    for (const char* label : {"Fr223", "Pb211"}) {
        const SandiaDecay::Nuclide* parent = db().nuclide(label);
        BOOST_REQUIRE(parent != nullptr);
        double expected = 0.0;
        double max_raw = 0.0;
        for (const SandiaDecay::Transition* t : parent->decaysToChildren) {
            if (!t || !t->level_feeds.empty()) continue;
            double raw = 0.0;
            for (const SandiaDecay::RadParticle& p : t->products)
                if (p.type == SandiaDecay::BetaParticle && p.intensity > 0.0f)
                    raw += static_cast<double>(p.intensity);
            max_raw = std::max(max_raw, raw);
            expected += static_cast<double>(t->branchRatio) * std::min(1.0, raw);
        }
        BOOST_TEST_CONTEXT(label) {
            BOOST_CHECK_GT(max_raw, 1.0);
            const RadioactiveEmissionSet e = build_radioactive_emissions(parent);
            double beta_minus = 0.0;
            for (const BetaBranch& b : e.beta_branches)
                if (!b.is_positron) beta_minus += b.yield_per_decay;
            BOOST_CHECK_CLOSE(beta_minus, expected, 1e-9);
            BOOST_CHECK_LE(e.beta_yield(), 1.0 + 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE(parent_beta_overflow_fails_closed_with_exact_throw_set) {
    const std::set<std::string> expected = {
        "Am244", "Ce148", "Hf181", "In131", "In131m2", "Lu179",
        "Nd152", "Pa235", "Pm155", "Sn129m", "Te136", "Tm151m"
    };
    std::set<std::string> actual;
    for (const SandiaDecay::Nuclide* parent : db().nuclides()) {
        if (!parent) continue;
        try {
            const RadioactiveEmissionSet e = build_radioactive_emissions(parent);
            BOOST_TEST_CONTEXT(parent->symbol)
                BOOST_CHECK_LE(e.beta_yield(), 1.0 + 1e-12);
        } catch (const std::runtime_error&) {
            actual.insert(parent->symbol);
        }
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(
        actual.begin(), actual.end(), expected.begin(), expected.end());

    // Duplicate an authoritative exact-feed transition by hand. The adapter
    // must not rescale those trusted feeds to conceal the parent inconsistency.
    SandiaDecay::Nuclide* th234 =
        const_cast<SandiaDecay::Nuclide*>(db().nuclide("Th234"));
    BOOST_REQUIRE(th234 != nullptr);
    const auto saved = th234->decaysToChildren;
    const SandiaDecay::Transition* exact = nullptr;
    for (const SandiaDecay::Transition* t : saved)
        if (t && t->mode == SandiaDecay::BetaDecay && !t->level_feeds.empty()) {
            exact = t;
            break;
        }
    BOOST_REQUIRE(exact != nullptr);
    th234->decaysToChildren.push_back(exact);
    const bool threw = [&]() {
        try { (void)build_radioactive_emissions(th234); }
        catch (const std::runtime_error&) { return true; }
        return false;
    }();
    th234->decaysToChildren = saved;
    BOOST_CHECK(threw);
    BOOST_CHECK_NO_THROW(build_radioactive_emissions(th234));
}

// Retargeted from U235 to Pa234 (2026-08-03): four unsupported U235 gamma
// records were corrected/removed against ENSDF + microcalorimeter evidence, so
// U235 -> Th231 now has a VALID level scheme (raw feed 2.373 -> 1.0013) and no
// longer builds a fallback forest at all. Pa234 -> U234 remains flux
// inconsistent (raw feed 3.19) and exercises the same code path.
BOOST_AUTO_TEST_CASE(fallback_forest_preserves_selected_member_marginals) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    const std::vector<DecayCascade> cascades = build_cascades(db(), "Pa234", opts);
    std::mt19937_64 rng(0x235f011ULL);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    constexpr int n = 300000;
    int checked = 0;
    for (const DecayCascade& dc : cascades) {
        if (dc.level_scheme.valid || dc.members.empty()) continue;
        const auto forest = build_cascade_fallback_forest(dc);
        std::vector<int> count(dc.members.size(), 0);
        std::vector<char> emitted;
        double positron_mass = 0.0;
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
            if (o.kind == WeakOutcomeKind::Positron)
                positron_mass += o.selected_mass;
        for (int k = 0; k < n; ++k) {
            const int forced = dc.weak_outcome_law.usable()
                ? (u(rng) < positron_mass ? 1 : 0) : -1;
            sample_cascade_fallback_forest(
                forest, forced, [&]() { return u(rng); }, emitted);
            for (std::size_t m = 0; m < emitted.size(); ++m)
                count[m] += emitted[m];
        }
        for (const CascadeFallbackNode& node : forest) {
            const double p = node.marginal;
            if (!(p >= 1e-4 && p <= 1.0 - 1e-4)) continue;
            const double observed = count[node.member] / static_cast<double>(n);
            const double sigma = std::sqrt(p * (1.0 - p) / n);
            const double energy = dc.members[node.member].energy_keV;
            if (std::abs(energy - 19.59) < 0.1
                || std::abs(energy - 202.11) < 0.1) {
                BOOST_TEST_MESSAGE("U235 fallback E=" << energy
                    << " expected=" << p << " observed=" << observed
                    << " count=" << count[node.member] << "/" << n
                    << " z=" << (observed - p) / sigma);
            }
            BOOST_TEST_CONTEXT("member " << node.member << " E="
                               << energy) {
                BOOST_CHECK_LT(std::abs(observed - p) / sigma, 5.0);
            }
            ++checked;
        }
    }
    BOOST_CHECK_GT(checked, 0);
}

// ---------------------------------------------------------------------------
// Level-scheme flux conservation.
//
// A decay enters its daughter's level scheme at most once, so the sum of direct
// feeding over levels is bounded by 1. Summing (out_flow - in_flow) over the
// closed level graph telescopes to zero, so the levels with net out-flow (the
// direct beta/EC/alpha feeding) must balance exactly against those with net
// in-flow. Nothing enforced either property before: an E0 transition's sentinel
// conversion coefficient (iccTotal up to 3.7e22) inflated one level's out-flow
// by ~1e22, and the max(0, .) in the feeding derivation silently destroyed the
// flux that arrived at the levels it fed.
//
// The feeding sum is BELOW 1 whenever some decays populate the daughter's ground
// or isomeric state directly, so 1.0 is an upper bound, not a target: Th-234
// feeds the Pa-234m isomer directly in ~71% of decays and correctly sums to
// ~0.29. Only daughters fed entirely through excited levels sum to 1.

namespace {

struct FluxCheck {
    double total_feed = 0.0, total_sink = 0.0;
    double worst_split = 0.0;   // max exclusive outcome sum over transitions
    int n_sinks = 0;
};

FluxCheck flux_check(const LevelScheme& ls) {
    FluxCheck f;
    const std::size_t n = ls.levels.size();
    std::vector<double> outf(n, 0.0), inf(n, 0.0);
    for (std::size_t L = 0; L < n; ++L)
        for (const auto& t : ls.levels[L].out) {
            outf[L] += t.weight;
            if (t.to_level >= 0 && static_cast<std::size_t>(t.to_level) < n)
                inf[static_cast<std::size_t>(t.to_level)] += t.weight;
            const double split = t.p_gamma + t.p_icK + t.p_icL1 + t.p_icL2
                               + t.p_icL3 + t.p_icUnresolved + t.p_unmodeled;
            f.worst_split = std::max(f.worst_split, split);
        }
    for (std::size_t L = 0; L < n; ++L) {
        const double net = outf[L] - inf[L];
        f.total_feed += std::max(0.0, net);
        if (net < -1e-12) { ++f.n_sinks; f.total_sink += -net; }
    }
    return f;
}

} // namespace

BOOST_AUTO_TEST_CASE(level_scheme_feeding_conserves_flux) {
    // The nuclides whose feeding sums the E0/ICC audit reported as damaged,
    // plus the three that were already sound (Co-60, Ba-133, Am-241).
    const char* nuclides[] = {"Co60", "Ba133", "Am241", "Th234", "Pa234m",
                              "Pa234", "U235", "U238", "Pu238", "Pu239",
                              "Pu240", "Bi212", "Bi214"};
    CascadeOptions opts;
    opts.prompt_equilibrium = false;   // isolate each named parent's own branches

    for (const char* nuc : nuclides) {
        BOOST_TEST_CONTEXT("nuclide " << nuc) {
            BOOST_REQUIRE(db().nuclide(nuc) != nullptr);
            for (const DecayCascade& dc : build_cascades(db(), nuc, opts)) {
                if (!dc.level_scheme.valid)
                    continue;
                const FluxCheck f = flux_check(dc.level_scheme);

                // A decay enters the scheme at most once.
                BOOST_CHECK_LE(f.total_feed, 1.05);

                // Sources balance sinks exactly (telescoping identity).
                BOOST_CHECK_SMALL(f.total_feed - f.total_sink,
                                  1e-9 * std::max(1.0, f.total_feed));

                // entry_probability is what consumers multiply back in to undo
                // LevelDag's normalization, so it must equal the feeding sum.
                BOOST_CHECK_CLOSE(dc.level_scheme.entry_probability,
                                  std::min(1.0, f.total_feed), 1e-6);
                BOOST_CHECK_GT(dc.level_scheme.entry_probability, 0.0);
                BOOST_CHECK_LE(dc.level_scheme.entry_probability, 1.0);

                // Per-occurrence outcomes are mutually exclusive.
                BOOST_CHECK_LE(f.worst_split, 1.0 + 1e-9);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(no_ground_feeding_nuclides_sum_to_one) {
    // Co-60, Ba-133 and Am-241 populate their daughter only through excited
    // levels, so for them the bound above is attained.
    for (const char* nuc : {"Co60", "Ba133", "Am241"}) {
        BOOST_TEST_CONTEXT("nuclide " << nuc) {
            CascadeOptions opts;
            opts.prompt_equilibrium = false;
            for (const DecayCascade& dc : build_cascades(db(), nuc, opts)) {
                if (!dc.level_scheme.valid || dc.branch_weight < 0.5)
                    continue;   // the dominant branch only
                BOOST_CHECK_CLOSE(dc.level_scheme.entry_probability, 1.0, 3.0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(whole_database_rejection_rate_is_stable) {
    // Sweeping only the ACCEPTED schemes for feeding <= 1.05 would be circular:
    // the adapter sets valid == true iff the feeding sum is in (0, 1.05], so such
    // a check can never fail and would not notice a data regression -- a refresh
    // that broke more transitions would simply reject more branches and still
    // pass. The quantity that actually moves is the REJECTION count, so freeze
    // that instead, together with the count of repairs the adapter had to apply.
    const std::vector<const SandiaDecay::Nuclide*>& all = db().nuclides();
    BOOST_REQUIRE_GT(all.size(), 1000u);

    CascadeOptions opts;
    opts.prompt_equilibrium = false;

    int n_valid = 0, n_raw_only = 0, n_partial_only = 0, n_both = 0;
    int n_invalid_topology_edges = 0;
    int n_e0 = 0, n_capped = 0;
    int n_accepted_e0 = 0, n_accepted_capped = 0;
    int n_graph_e0 = 0, n_residual_e0 = 0;
    int n_graph_capped = 0, n_residual_capped = 0;
    for (const SandiaDecay::Nuclide* nuc : all) {
        if (!nuc) continue;
        std::vector<BranchFeedingAudit> audit;
        const std::vector<DecayCascade> cs = build_cascades(nuc, opts, &audit);
        for (const BranchFeedingAudit& a : audit) {
            if (a.scheme_valid) {
                ++n_valid;
                n_accepted_e0 += a.n_e0_repaired;
                n_accepted_capped += a.n_intensity_capped;
            } else if (a.rejected && a.partial_scheme_rejected) {
                ++n_both;
            } else if (a.rejected) {
                ++n_raw_only;
            } else if (a.partial_scheme_rejected) {
                ++n_partial_only;
            }
            n_invalid_topology_edges += a.n_invalid_topology_edges;
            n_e0 += a.n_e0_repaired;
            n_capped += a.n_intensity_capped;
            n_graph_e0 += a.n_graph_e0_repaired;
            n_residual_e0 += a.n_residual_e0_repaired;
            n_graph_capped += a.n_graph_intensity_capped;
            n_residual_capped += a.n_residual_intensity_capped;
        }
        // Every accepted scheme still has to satisfy the invariant.
        for (const DecayCascade& dc : cs) {
            if (!dc.level_scheme.valid) continue;
            BOOST_REQUIRE_LE(flux_check(dc.level_scheme).total_feed, 1.05);
        }
    }
    BOOST_TEST_MESSAGE("valid=" << n_valid
                       << " raw_only=" << n_raw_only
                       << " partial_only=" << n_partial_only
                       << " raw_and_partial=" << n_both
                       << " invalid_topology_edges=" << n_invalid_topology_edges
                       << " e0_repaired=" << n_e0
                       << " (graph=" << n_graph_e0
                       << ", residual=" << n_residual_e0 << ")"
                       << " intensity_capped=" << n_capped
                       << " (graph=" << n_graph_capped
                       << ", residual=" << n_residual_capped << ")"
                       << " accepted_repairs=" << n_accepted_e0
                       << "/" << n_accepted_capped);

    // Exact branch-instance census for the reviewed regenerated XML pinned by
    // SANDIA_DECAY_XML_PATH. The old single `rejected` count is no longer a
    // meaningful baseline: a graph can independently fail raw flux and the
    // partial-topology policy. In the old code partial graphs were rejected
    // before their raw flux was inspected, so summing the two flags made the
    // apparent raw rejection count almost double. Keep the categories disjoint,
    // and freeze repair counts only for accepted graphs (repairs inspected on a
    // graph that later falls back do not affect engine output).
    // 2026-08-03, two data corrections, both raising valid at raw-only's expense:
    //   +28  U235  -> Th231, four unsupported gamma records corrected/removed
    //                against ENSDF NDS 185 (2022) and microcalorimeter data
    //   +26  Th231 -> Pa231, gamma set rebuilt from ENSDF (19 -> 47 records);
    //                its feed sum went 1.81827 -> 0.97581 and every fed level
    //                now agrees with the shipped evaluated feed law to <= 6%
    // Each branch appears once per decay chain that reaches it, hence the
    // multiplicity. No other nuclide moved.
    BOOST_CHECK_EQUAL(n_valid, 4558);
    BOOST_CHECK_EQUAL(n_raw_only, 763);
    BOOST_CHECK_EQUAL(n_partial_only, 1782);
    BOOST_CHECK_EQUAL(n_both, 1328);
    BOOST_CHECK_EQUAL(n_invalid_topology_edges, 0);
    BOOST_CHECK_EQUAL(n_accepted_e0, 323); // includes reviewed Eu152 615.4-keV E0
    BOOST_CHECK_EQUAL(n_accepted_capped, 801);
}

BOOST_AUTO_TEST_CASE(e0_transitions_emit_no_photon) {
    // A definitive ICC-sentinel E0 cannot emit a single gamma. Its conversion
    // coefficient is effectively infinite, so it must never become an emitting
    // member -- on the level-scheme path p_gamma = 0 would suppress it anyway,
    // but the pairwise/fallback path fires members at their tabulated intensity
    // and would invent the photon. Hg-182 -> Pt-178 lists its 422 keV E0 at 64%
    // per decay, so the failure mode is not subtle. The transition must still
    // carry its flux in the level graph, with gamma_member = -1.
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    int n_e0_transitions = 0;
    for (const char* nuc : {"Hg182", "Pa230", "Pa234", "Bi212", "Bi214", "Pu238"}) {
        const SandiaDecay::Nuclide* parent = db().nuclide(nuc);
        if (!parent) continue;
        BOOST_TEST_CONTEXT("nuclide " << nuc) {
            std::vector<BranchFeedingAudit> audit;
            const std::vector<DecayCascade> cs = build_cascades(parent, opts, &audit);
            int repaired = 0;
            for (const BranchFeedingAudit& a : audit) repaired += a.n_e0_repaired;
            if (repaired == 0) continue;

            for (const DecayCascade& dc : cs) {
                if (!dc.level_scheme.valid) continue;
                for (const CascadeLevel& lv : dc.level_scheme.levels)
                    for (const LevelOutTransition& t : lv.out) {
                        if (t.gamma_member >= 0) continue;
                        ++n_e0_transitions;
                        // A memberless transition is legal only as an E0: it emits
                        // nothing, but still carries flux and conversion vacancies.
                        BOOST_CHECK_EQUAL(t.p_gamma, 0.0);
                        BOOST_CHECK_GT(t.weight, 0.0);
                        // and no member may sit at its energy
                        for (const CascadeMember& m : dc.members)
                            if (m.type == CascadeParticleType::Gamma)
                                BOOST_CHECK_GT(std::abs(m.energy_keV - t.gamma_keV), 0.05);
                    }
            }
        }
    }
    BOOST_TEST_MESSAGE("checked " << n_e0_transitions << " memberless E0 transitions");
    BOOST_CHECK_GT(n_e0_transitions, 0);
}

BOOST_AUTO_TEST_CASE(gamma_member_intensity_respects_transition_bound) {
    // A transition occurs at most once per decay and emits its gamma on
    // 1/(1+alpha) of those, so I_gamma <= 1/(1+alpha) is a hard ceiling. Members
    // are emitted at their own intensity, so without this the branches the flux
    // guard rejects -- exactly the ones with known-inconsistent data -- emit the
    // offending line at its raw rate. U-235's 19.59 keV line pairs I = 0.61 with
    // alpha = 114.7 and was emitted 21x above the GEANT4/ENSDF scale, feeding
    // false sum peaks onto the 205.31 and 221.40 keV lines.
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    int n_bounded = 0;
    for (const char* nuc : {"U235", "I125", "Ba133", "Co60", "Am241", "Pa234"}) {
        const SandiaDecay::Nuclide* parent = db().nuclide(nuc);
        BOOST_REQUIRE(parent != nullptr);
        BOOST_TEST_CONTEXT("nuclide " << nuc) {
            for (const DecayCascade& dc : build_cascades(parent, opts)) {
                for (const SandiaDecay::Transition* tr : parent->decaysToChildren) {
                    if (!tr) continue;
                    for (const SandiaDecay::RadParticle& p : tr->products) {
                        if (p.type != SandiaDecay::ProductType::GammaParticle) continue;
                        const double a = std::max(0.0, double(p.icc_total));
                        const double bound = 1.0 / (1.0 + a);
                        for (const CascadeMember& m : dc.members) {
                            if (m.type != CascadeParticleType::Gamma) continue;
                            if (std::abs(m.energy_keV - double(p.energy)) > 1e-3) continue;
                            BOOST_CHECK_LE(m.intensity, bound + 1e-9);
                            if (double(p.intensity) > bound + 1e-9) ++n_bounded;
                        }
                    }
                }
            }
        }
    }
    BOOST_TEST_MESSAGE("members actually bounded by 1/(1+alpha): " << n_bounded);
    BOOST_CHECK_GT(n_bounded, 0);   // U-235's 19.59 keV at minimum
}

BOOST_AUTO_TEST_CASE(fallback_branches_keep_e0_conversion_vacancies) {
    // A definitive below-pair sentinel E0 emits no photon but is fully converted,
    // so it is a near-certain
    // vacancy source. On the fallback path there is no level scheme to carry it,
    // and skipping it (as the memberless representation first did) discarded
    // ~1.12 K vacancies per parent decay across the database. The vacancy is
    // unconditional -- there is no gamma to condition on -- hence gamma_member -1.
    // It is one exclusive group, not four independent shell Bernoullis.
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    for (const char* nuc : {"Hg182", "Pa230", "Co68m", "Y98"}) {
        const SandiaDecay::Nuclide* parent = db().nuclide(nuc);
        if (!parent) continue;
        BOOST_TEST_CONTEXT("nuclide " << nuc) {
            double kv = 0.0;
            for (const DecayCascade& dc : build_cascades(parent, opts)) {
                if (dc.level_scheme.valid) continue;
                BOOST_CHECK(dc.k_vacancies.empty());
                for (const VacancyGroup& v : dc.vacancy_groups) {
                    if (v.kind != VacancyGroupKind::ElectricMonopole) continue;
                    BOOST_CHECK_EQUAL(v.gamma_member, -1);
                    BOOST_CHECK_GT(v.transition_energy_keV, 0.0);
                    BOOST_CHECK_CLOSE(v.p_none + v.p_K + v.p_L1 + v.p_L2
                                          + v.p_L3 + v.p_outer + v.p_unmodeled,
                                      1.0, 1e-9);
                    kv += dc.branch_weight * v.p_K;
                    if (v.transition_energy_keV >= 2.0 * 510.998950)
                        BOOST_CHECK_EQUAL(v.p_outer, 0.0);
                }
            }
            BOOST_CHECK_GT(kv, 0.01);
        }
    }
}

BOOST_AUTO_TEST_CASE(ec_and_positron_probabilities_are_partitionable) {
    // Electron capture and beta+ compete for the same decay, so the realization
    // partitions one uniform over {positron, K capture, L capture, neither}. That
    // is only valid where the three fit inside one uniform; Rb-82 (95% beta+,
    // 5% EC) is the case that motivated it -- sampling them independently put a
    // spurious ~0.039 K-vacancy per decay in coincidence with the 511 pair.
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    const SandiaDecay::Nuclide* rb = db().nuclide("Rb82");
    BOOST_REQUIRE(rb != nullptr);
    double p_pos = 0.0, p_k = 0.0;
    for (const DecayCascade& dc : build_cascades(rb, opts)) {
        BOOST_REQUIRE(dc.weak_outcome_law.usable());
        double categorical_sum = 0.0;
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
            categorical_sum += o.selected_mass;
            if (o.kind == WeakOutcomeKind::Positron)
                p_pos += dc.branch_weight * o.selected_mass;
            else if (o.kind == WeakOutcomeKind::ElectronCapture)
                p_k += dc.branch_weight * o.selected_mass * o.ec_K;
        }
        BOOST_CHECK_CLOSE(categorical_sum, 1.0, 1e-9);
    }
    // Deterministic probability accounting: 0.0403469 K vacancies / decay.
    BOOST_CHECK_CLOSE(p_k, 0.0403469, 0.2);
    // One categorical draw chooses EC OR positron OR Other, so overlap is zero
    // by construction even though the marginals nearly saturate one.
    BOOST_CHECK_LE(p_pos + p_k, 1.0);
}

BOOST_AUTO_TEST_CASE(eu152_weak_products_are_transition_conditional) {
    // The source Eu-152 EC/beta+ vector is parent-absolute and closes to the
    // 0.721 Sm branch.  Enrichment canonicalizes it to a branch-conditional
    // law, so build_cascades must multiply by branch_weight exactly once while
    // retaining the small genuinely unrepresented outcome.
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    const SandiaDecay::Nuclide* eu = db().nuclide("Eu152");
    BOOST_REQUIRE(eu != nullptr);
    bool found_sm = false;
    double represented_parent_mass = 0.0;
    double ec_k_parent_mass = 0.0;
    double positron_parent_mass = 0.0;
    double other_parent_mass = 0.0;
    for (const DecayCascade& dc : build_cascades(eu, opts)) {
        if (dc.daughter_Z != 62 || std::abs(dc.branch_weight - 0.721) > 1e-5)
            continue;
        found_sm = true;
        BOOST_REQUIRE(dc.weak_outcome_law.usable());
        double categorical_sum = 0.0;
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
            categorical_sum += o.selected_mass;
            if (o.kind == WeakOutcomeKind::ElectronCapture) {
                represented_parent_mass += dc.branch_weight * o.selected_mass;
                ec_k_parent_mass +=
                    dc.branch_weight * o.selected_mass * o.ec_K;
            } else if (o.kind == WeakOutcomeKind::Positron) {
                represented_parent_mass += dc.branch_weight * o.selected_mass;
                positron_parent_mass += dc.branch_weight * o.selected_mass;
            } else {
                other_parent_mass += dc.branch_weight * o.selected_mass;
            }
        }
        BOOST_CHECK_CLOSE(categorical_sum, 1.0, 1e-9);
    }
    BOOST_REQUIRE(found_sm);
    BOOST_CHECK_CLOSE(represented_parent_mass, 0.720282, 2e-5);
    BOOST_CHECK_CLOSE(ec_k_parent_mass, 0.600803823, 2e-5);
    BOOST_CHECK_CLOSE(positron_parent_mass, 0.000138, 2e-4);
    BOOST_CHECK_SMALL(other_parent_mass - (0.721 - 0.720282), 5e-8);
}

BOOST_AUTO_TEST_CASE(eu152_615_e0_is_memberless_level2_to_ground) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    const SandiaDecay::Nuclide* eu = db().nuclide("Eu152");
    BOOST_REQUIRE(eu != nullptr);
    bool found_e0 = false;
    for (const DecayCascade& dc : build_cascades(eu, opts)) {
        if (dc.daughter_Z != 64 || !dc.level_scheme.valid) continue;
        for (std::size_t level_index = 0;
             level_index < dc.level_scheme.levels.size(); ++level_index)
            for (const LevelOutTransition& t :
                 dc.level_scheme.levels[level_index].out) {
                if (std::abs(t.gamma_keV - 615.4) > 0.1) continue;
                found_e0 = true;
                BOOST_CHECK_EQUAL(level_index, 2u);
                BOOST_CHECK_EQUAL(t.to_level, 0);
                BOOST_CHECK_EQUAL(t.gamma_member, -1);
                BOOST_CHECK_EQUAL(t.p_gamma, 0.0);
                BOOST_CHECK_CLOSE(t.p_icK, 0.877, 0.02);
                BOOST_CHECK_CLOSE(t.p_icUnresolved, 0.123, 0.2);
                BOOST_CHECK_EQUAL(t.p_unmodeled, 0.0);
                BOOST_CHECK_CLOSE(t.p_icK + t.p_icL1 + t.p_icL2 + t.p_icL3
                                      + t.p_icUnresolved + t.p_unmodeled,
                                  1.0, 1e-9);
            }
    }
    BOOST_CHECK(found_e0);

    double direct_615 = 0.0;
    for (const RadioactivePhotonLine& p :
         build_radioactive_emissions(eu).photons)
        if (std::abs(p.energy_keV - 615.4) < 0.1)
            direct_615 += p.yield_per_decay;
    BOOST_CHECK_SMALL(direct_615, 1e-14);
}

BOOST_AUTO_TEST_CASE(eu152_nocoinc_keeps_physical_ec_vacancy_yield) {
    std::string path = SANDIA_DECAY_XML_PATH;
    const std::string full_name = "sandia.decay.xml";
    const std::size_t at = path.rfind(full_name);
    BOOST_REQUIRE(at != std::string::npos);
    path.replace(at, full_name.size(), "sandia.decay.nocoinc.min.xml");
    const SandiaDecay::SandiaDecayDataBase nocoinc(path);
    const SandiaDecay::Nuclide* eu = nocoinc.nuclide("Eu152");
    BOOST_REQUIRE(eu != nullptr);
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    double ec_k_parent_mass = 0.0;
    for (const DecayCascade& dc : build_cascades(eu, opts))
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
            if (o.kind == WeakOutcomeKind::ElectronCapture)
                ec_k_parent_mass +=
                    dc.branch_weight * o.selected_mass * o.ec_K;
    // The minimized file rounds ecK to four significant digits.
    BOOST_CHECK_CLOSE(ec_k_parent_mass, 0.600803823, 0.01);
}

BOOST_AUTO_TEST_CASE(zr89_preserves_terminal_isomer_capture) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    const SandiaDecay::Nuclide* zr = db().nuclide("Zr89");
    BOOST_REQUIRE(zr != nullptr);
    bool found_level1_ec = false;
    double k_yield = 0.0;
    for (const DecayCascade& dc : build_cascades(zr, opts)) {
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
            if (o.kind != WeakOutcomeKind::ElectronCapture) continue;
            if (o.fed_level == 1) found_level1_ec = true;
            k_yield += dc.branch_weight * o.selected_mass * o.ec_K;
        }
    }
    BOOST_CHECK(found_level1_ec);
    BOOST_CHECK_CLOSE(k_yield, 0.67447, 0.2);
}

BOOST_AUTO_TEST_CASE(vacancy_only_ec_branch_is_retained) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    opts.include_xrays = false; // force the pure-EC branch to be memberless
    opts.min_intensity = 1e-6;  // omit Fe-55's vanishing 126-keV gamma
    const SandiaDecay::Nuclide* fe = db().nuclide("Fe55");
    BOOST_REQUIRE(fe != nullptr);
    bool found_memberless_ec = false;
    for (const DecayCascade& dc : build_cascades(fe, opts)) {
        if (!dc.members.empty()) continue;
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
            if (o.kind == WeakOutcomeKind::ElectronCapture && o.ec_K > 0.0)
                found_memberless_ec = true;
    }
    BOOST_CHECK(found_memberless_ec);
}

BOOST_AUTO_TEST_CASE(direct_u235_gamma_respects_transition_ceiling) {
    const RadioactiveEmissionSet emissions =
        build_radioactive_emissions(db(), "U235");
    double yield_19 = 0.0;
    for (const RadioactivePhotonLine& line : emissions.photons)
        if (std::abs(line.energy_keV - 19.59) < 0.1)
            yield_19 += line.yield_per_decay;
    BOOST_CHECK_GT(yield_19, 0.0);
    BOOST_CHECK_LE(yield_19, 1.0 / (1.0 + 114.7) + 1e-6);
    BOOST_CHECK_LT(yield_19, 0.01); // old direct API returned the raw 0.61
}

BOOST_AUTO_TEST_CASE(advisory_mult1_lines_remain_photons) {
    // G4 code 1 alone is not sufficient E0 provenance: ENSDF leaves these
    // observed lines' M/CC blank or ambiguous. They have no ICC sentinel and
    // must retain their tabulated photon yields.
    for (const auto& c : std::vector<std::pair<const char*, std::pair<double,double>>>{
             {"Ra232", {373.3, 0.62}}, {"Sr78", {189.8, 0.16}}}) {
        const RadioactiveEmissionSet e = build_radioactive_emissions(db(), c.first);
        double y = 0.0;
        for (const RadioactivePhotonLine& p : e.photons)
            if (std::abs(p.energy_keV - c.second.first) < 0.1)
                y += p.yield_per_decay;
        BOOST_TEST_CONTEXT(c.first)
            BOOST_CHECK_CLOSE(y, c.second.second, 0.2);
    }
}

BOOST_AUTO_TEST_CASE(definitive_below_pair_e0_keeps_unresolved_conversion) {
    // Nb-98's 734.6-keV transition has the definitive alpha=2e22 sentinel.
    // Only K is persisted for this light daughter; the remaining 9.82% is a
    // shell-unresolved conversion, not a photon and not four independent draws.
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    bool found = false;
    for (const DecayCascade& dc : build_cascades(db(), "Nb98", opts)) {
        if (!dc.level_scheme.valid) continue;
        for (const CascadeLevel& lv : dc.level_scheme.levels)
            for (const LevelOutTransition& t : lv.out) {
                if (std::abs(t.gamma_keV - 734.6) > 0.1) continue;
                found = true;
                BOOST_CHECK_EQUAL(t.gamma_member, -1);
                BOOST_CHECK_EQUAL(t.p_gamma, 0.0);
                BOOST_CHECK_CLOSE(t.p_icK, 0.9018, 0.02);
                BOOST_CHECK_CLOSE(t.p_icUnresolved, 0.0982, 0.2);
                BOOST_CHECK_CLOSE(t.p_icK + t.p_icL1 + t.p_icL2
                                      + t.p_icL3 + t.p_icUnresolved
                                      + t.p_unmodeled,
                                  1.0, 1e-9);
                BOOST_CHECK_EQUAL(t.p_unmodeled, 0.0);
            }
    }
    BOOST_CHECK(found);

    const RadioactiveEmissionSet direct =
        build_radioactive_emissions(db(), "Nb98");
    double photon = 0.0, unresolved = 0.0;
    for (const RadioactivePhotonLine& p : direct.photons)
        if (std::abs(p.energy_keV - 734.6) < 0.1) photon += p.yield_per_decay;
    for (const ConversionElectronLine& e : direct.conversion_electrons)
        if (e.shell == ConversionShell::Outer
            && std::abs(e.energy_keV - 734.6) < 0.1)
            unresolved += e.yield_per_decay;
    BOOST_CHECK_SMALL(photon, 1e-14);
    BOOST_CHECK_CLOSE(unresolved, 0.26 * 0.0982, 0.3);
}

BOOST_AUTO_TEST_CASE(verified_above_pair_e0_is_explicitly_unmodeled) {
    // N-16's evaluated 6048.2-keV 0+ -> 0+ is authoritative E0 provenance,
    // independently of G4's advisory code 1 and absent ICC sentinel. It cannot
    // emit a single photon. With no evaluated pair/conversion split, its one
    // conditional outcome is explicitly unmodelled and emits no substitute.
    const SandiaDecay::Nuclide* n16 = db().nuclide("N16");
    BOOST_REQUIRE(n16 != nullptr);
    bool verified_record = false;
    for (const SandiaDecay::Transition* tr : n16->decaysToChildren)
        if (tr)
            for (const SandiaDecay::RadParticle& p : tr->products)
                if (p.type == SandiaDecay::ProductType::GammaParticle
                    && std::abs(double(p.energy) - 6048.2) < 0.2) {
                    verified_record = true;
                    BOOST_CHECK(p.e0_verified);
                }
    BOOST_REQUIRE(verified_record);

    const RadioactiveEmissionSet e = build_radioactive_emissions(n16);
    double photon = 0.0, conversion = 0.0;
    for (const RadioactivePhotonLine& p : e.photons)
        if (std::abs(p.energy_keV - 6048.2) < 0.2) photon += p.yield_per_decay;
    for (const ConversionElectronLine& c : e.conversion_electrons)
        if (std::abs(c.energy_keV - 6048.2) < 0.2)
            conversion += c.yield_per_decay;
    BOOST_CHECK_SMALL(photon, 1e-14);
    BOOST_CHECK_SMALL(conversion, 1e-14);

    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    bool explicit_unmodeled = false;
    for (const DecayCascade& dc : build_cascades(n16, opts)) {
        for (const CascadeMember& m : dc.members)
            if (m.type == CascadeParticleType::Gamma)
                BOOST_CHECK_GT(std::abs(m.energy_keV - 6048.2), 0.2);
        for (const CascadeLevel& level : dc.level_scheme.levels)
            for (const LevelOutTransition& t : level.out)
                if (std::abs(t.gamma_keV - 6048.2) < 0.2) {
                    BOOST_CHECK_EQUAL(t.gamma_member, -1);
                    BOOST_CHECK_EQUAL(t.p_gamma, 0.0);
                    BOOST_CHECK_CLOSE(t.p_unmodeled, 1.0, 1e-9);
                    explicit_unmodeled = true;
                }
        for (const VacancyGroup& g : dc.vacancy_groups)
            if (g.kind == VacancyGroupKind::ElectricMonopole
                && std::abs(g.transition_energy_keV - 6048.2) < 0.2) {
                const double occurs = 1.0 - g.p_none;
                BOOST_REQUIRE_GT(occurs, 0.0);
                BOOST_CHECK_CLOSE(g.p_unmodeled / occurs, 1.0, 1e-9);
                BOOST_CHECK_EQUAL(g.p_K + g.p_L1 + g.p_L2 + g.p_L3
                                      + g.p_outer,
                                  0.0);
                explicit_unmodeled = true;
            }
    }
    BOOST_CHECK(explicit_unmodeled);

    // Disabling the vacancy model exercises the legacy member path: verified E0
    // suppression must not depend on building a LevelScheme or VacancyGroup.
    opts.vacancy_xray_model = false;
    for (const DecayCascade& dc : build_cascades(n16, opts))
        for (const CascadeMember& m : dc.members)
            if (m.type == CascadeParticleType::Gamma)
                BOOST_CHECK_GT(std::abs(m.energy_keV - 6048.2), 0.2);
}

BOOST_AUTO_TEST_CASE(above_pair_sentinel_remainder_is_not_conversion) {
    // Co-68m's 2511.9-keV transition is a definitive alpha=1.38e22 sentinel,
    // but above 2m_e the unresolved remainder can be internal-pair formation.
    // It must not be emitted as a photon OR mislabeled as an outer conversion.
    const RadioactiveEmissionSet e = build_radioactive_emissions(db(), "Co68m");
    double photon = 0.0, outer = 0.0;
    for (const RadioactivePhotonLine& p : e.photons)
        if (std::abs(p.energy_keV - 2511.9) < 0.2) photon += p.yield_per_decay;
    for (const ConversionElectronLine& c : e.conversion_electrons)
        if (c.shell == ConversionShell::Outer
            && std::abs(c.energy_keV - 2511.9) < 0.2)
            outer += c.yield_per_decay;
    BOOST_CHECK_SMALL(photon, 1e-14);
    BOOST_CHECK_SMALL(outer, 1e-14);
}

BOOST_AUTO_TEST_CASE(positron_consumers_use_selected_weak_mass) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    for (const char* name : {"Zr85m", "Rb82"}) {
        const SandiaDecay::Nuclide* parent = db().nuclide(name);
        BOOST_REQUIRE(parent != nullptr);
        double law_pos = 0.0, annih = 0.0;
        for (const DecayCascade& dc : build_cascades(parent, opts)) {
            for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
                if (o.kind == WeakOutcomeKind::Positron)
                    law_pos += dc.branch_weight * o.selected_mass;
            for (const CascadeMember& m : dc.members)
                if (m.type == CascadeParticleType::Annih511)
                    annih += dc.branch_weight * m.intensity;
        }
        double direct_pos = 0.0;
        for (const BetaBranch& beta :
             build_radioactive_emissions(parent).beta_branches)
            if (beta.is_positron) direct_pos += beta.yield_per_decay;
        BOOST_TEST_CONTEXT(name) {
            BOOST_CHECK_LE(law_pos, 1.0 + 1e-12);
            BOOST_CHECK_CLOSE(annih, law_pos, 1e-8);
            BOOST_CHECK_CLOSE(direct_pos, law_pos, 1e-5);
        }
    }
}

BOOST_AUTO_TEST_CASE(accepted_overfull_feed_is_common_scaled) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    std::vector<BranchFeedingAudit> audit;
    const auto cascades = build_cascades(db().nuclide("Ba133"), opts, &audit);
    bool found = false;
    for (const BranchFeedingAudit& a : audit) {
        if (!(a.raw_total_feed > 1.0 && a.raw_total_feed <= 1.05)) continue;
        found = true;
        BOOST_CHECK_CLOSE(a.feeding_scale, 1.0 / a.raw_total_feed, 1e-9);
        BOOST_CHECK_CLOSE(a.total_feed, 1.0, 1e-9);
        BOOST_CHECK(!a.rejected);
    }
    for (const DecayCascade& dc : cascades)
        if (dc.level_scheme.valid && dc.branch_weight > 0.5)
            BOOST_CHECK_CLOSE(flux_check(dc.level_scheme).total_feed, 1.0, 1e-9);
    BOOST_CHECK(found);
}

BOOST_AUTO_TEST_CASE(historical_two_argument_overloads_are_symbols) {
    using ParentFn = std::vector<DecayCascade> (*)(
        const SandiaDecay::Nuclide*, const CascadeOptions&);
    using DbFn = std::vector<DecayCascade> (*)(
        const SandiaDecay::SandiaDecayDataBase&, const std::string&,
        const CascadeOptions&);
    ParentFn parent_fn = static_cast<ParentFn>(&build_cascades);
    DbFn db_fn = static_cast<DbFn>(&build_cascades);
    BOOST_CHECK(parent_fn != nullptr);
    BOOST_CHECK(db_fn != nullptr);
}

BOOST_AUTO_TEST_CASE(categorical_laws_are_bounded_database_wide) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    int n_weak = 0, n_groups = 0, n_low_confidence = 0;
    for (const SandiaDecay::Nuclide* nuc : db().nuclides()) {
        if (!nuc) continue;
        for (const DecayCascade& dc : build_cascades(nuc, opts)) {
            if (dc.weak_outcome_law.usable()) {
                ++n_weak;
                double sum = 0.0;
                for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                    BOOST_REQUIRE_GE(o.selected_mass, 0.0);
                    BOOST_REQUIRE_LE(o.selected_mass, 1.0);
                    sum += o.selected_mass;
                    if (dc.level_scheme.valid && o.fed_level >= 0)
                        BOOST_REQUIRE_LT(static_cast<std::size_t>(o.fed_level),
                                         dc.level_scheme.levels.size());
                }
                BOOST_REQUIRE_SMALL(sum - 1.0, 2e-14);
                if (dc.weak_outcome_law.confidence
                    == WeakOutcomeConfidence::CommonScaleLowConfidence)
                    ++n_low_confidence;
            }
            for (const VacancyGroup& g : dc.vacancy_groups) {
                ++n_groups;
                const double sum = g.p_none + g.p_K + g.p_L1 + g.p_L2
                    + g.p_L3 + g.p_outer + g.p_unmodeled;
                BOOST_REQUIRE_SMALL(sum - 1.0, 2e-14);
                if (g.kind == VacancyGroupKind::InternalConversion) {
                    BOOST_REQUIRE_GE(g.gamma_member, 0);
                    BOOST_REQUIRE_LT(static_cast<std::size_t>(g.gamma_member),
                                     dc.members.size());
                } else {
                    BOOST_REQUIRE_EQUAL(g.gamma_member, -1);
                    if (g.transition_energy_keV >= 2.0 * 510.998950)
                        BOOST_REQUIRE_EQUAL(g.p_outer, 0.0);
                }
            }
        }
    }
    BOOST_CHECK_GT(n_weak, 0);
    BOOST_CHECK_GT(n_groups, 0);
    BOOST_CHECK_GT(n_low_confidence, 0);
}

BOOST_AUTO_TEST_CASE(valid_partial_schemes_retain_categorical_residuals) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;

    for (const char* name : {"Th234", "Pb214", "Bi214"}) {
        std::vector<BranchFeedingAudit> audit;
        const auto cascades = build_cascades(db().nuclide(name), opts, &audit);
        int residuals = 0;
        double gamma_probability = 0.0;
        for (const DecayCascade& dc : cascades) {
            if (!dc.level_scheme.valid) continue;
            std::vector<int> dag_members;
            for (const CascadeLevel& level : dc.level_scheme.levels)
                for (const LevelOutTransition& t : level.out)
                    if (t.gamma_member >= 0) dag_members.push_back(t.gamma_member);
            for (const ResidualTransition& r : dc.residual_transitions) {
                ++residuals;
                const double sum = r.p_none + r.p_gamma + r.p_icK + r.p_icL1
                    + r.p_icL2 + r.p_icL3 + r.p_icOuter
                    + r.p_icUnresolved + r.p_unmodeled;
                BOOST_REQUIRE_SMALL(sum - 1.0, 2e-14);
                BOOST_REQUIRE_GE(r.gamma_member, 0);
                BOOST_REQUIRE_LT(static_cast<std::size_t>(r.gamma_member),
                                 dc.members.size());
                BOOST_CHECK(std::find(dag_members.begin(), dag_members.end(),
                                      r.gamma_member) == dag_members.end());
                BOOST_CHECK_CLOSE(r.p_gamma,
                    dc.members[static_cast<std::size_t>(r.gamma_member)].intensity,
                    1e-9);
                gamma_probability += r.p_gamma;
            }
        }
        int audited_residuals = 0;
        double audited_gamma_probability = 0.0;
        for (const BranchFeedingAudit& a : audit) {
            audited_residuals += a.n_residual_transitions;
            audited_gamma_probability += a.residual_gamma_probability;
        }
        BOOST_TEST_CONTEXT(name) {
            BOOST_CHECK_GT(residuals, 0);
            BOOST_CHECK_GT(gamma_probability, 0.0);
            BOOST_CHECK_EQUAL(residuals, audited_residuals);
            BOOST_CHECK_CLOSE(gamma_probability, audited_gamma_probability, 1e-9);
        }
    }
}

BOOST_AUTO_TEST_CASE(u238_partial_graph_is_valid_and_residuals_coexist) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    std::vector<BranchFeedingAudit> audit;
    const auto cascades = build_cascades(db().nuclide("U238"), opts, &audit);
    bool found = false;
    for (const DecayCascade& dc : cascades) {
        if (!dc.level_scheme.valid || dc.residual_transitions.empty()) continue;
        int matched = 0;
        for (const CascadeLevel& level : dc.level_scheme.levels)
            matched += static_cast<int>(level.out.size());
        if (matched < 2) continue;
        found = true;
        BOOST_CHECK_GT(dc.level_scheme.entry_probability, 0.19);
        BOOST_CHECK_LT(dc.level_scheme.entry_probability, 0.22);
        BOOST_CHECK_GE(dc.residual_transitions.size(), 4u);
        double residual_gamma = 0.0;
        for (const ResidualTransition& r : dc.residual_transitions)
            residual_gamma += r.p_gamma;
        BOOST_CHECK_GT(residual_gamma, 1e-4);
    }
    BOOST_CHECK(found);
    bool audited_partial = false;
    for (const BranchFeedingAudit& a : audit)
        if (a.scheme_valid && a.n_residual_transitions >= 4) {
            audited_partial = true;
            BOOST_CHECK_GT(a.residual_gamma_probability, 1e-4);
        }
    BOOST_CHECK(audited_partial);
}

BOOST_AUTO_TEST_CASE(partial_graph_occurrence_gate_uses_exact_feeds) {
    auto transition_to = [](const SandiaDecay::Nuclide* parent,
                            const char* child) {
        if (!parent) return static_cast<SandiaDecay::Transition*>(nullptr);
        for (const SandiaDecay::Transition* t : parent->decaysToChildren)
            if (t && t->child && t->child->symbol == child)
                return const_cast<SandiaDecay::Transition*>(t);
        return static_cast<SandiaDecay::Transition*>(nullptr);
    };
    auto audit_for = [](const SandiaDecay::Nuclide* parent,
                        const char* child) {
        CascadeOptions opts;
        opts.prompt_equilibrium = false;
        std::vector<BranchFeedingAudit> audit;
        (void)build_cascades(parent, opts, &audit);
        for (const BranchFeedingAudit& a : audit)
            if (a.parent == parent->symbol && a.child == child)
                return a;
        return BranchFeedingAudit{};
    };

    // U-238 has several topology-free photon records, but they carry only
    // 0.0001365 occurrence per branch. Its occurrence coverage is above 99%, so
    // the graph does not need exact starts and remains valid if they are absent.
    SandiaDecay::Transition* u238 = transition_to(db().nuclide("U238"), "Th234");
    BOOST_REQUIRE(u238 != nullptr);
    auto u238_feeds = u238->level_feeds;
    u238->level_feeds.clear();
    const BranchFeedingAudit au = audit_for(db().nuclide("U238"), "Th234");
    BOOST_CHECK(au.scheme_valid);
    BOOST_CHECK(!au.partial_scheme_rejected);
    BOOST_CHECK_GE(au.matched_transition_occurrence_mass,
                   0.99 * au.total_transition_occurrence_mass);
    u238->level_feeds = std::move(u238_feeds);

    // Tl-206 -> Pb-206 has low relative occurrence coverage but only 0.0010
    // absolute residual occurrence. Complete compatible starts admit it; absent
    // or incompatible starts fail closed.
    SandiaDecay::Transition* tl206 = transition_to(db().nuclide("Tl206"), "Pb206");
    BOOST_REQUIRE(tl206 != nullptr);
    BOOST_REQUIRE(!tl206->level_feeds.empty());
    auto tl206_feeds = tl206->level_feeds;
    const BranchFeedingAudit ap0 = audit_for(db().nuclide("Tl206"), "Pb206");
    BOOST_CHECK(ap0.scheme_valid);
    BOOST_CHECK(ap0.partial_exact_feed_compatible);
    BOOST_CHECK_LT(ap0.matched_transition_occurrence_mass,
                   0.99 * ap0.total_transition_occurrence_mass);
    BOOST_CHECK_LE(ap0.total_transition_occurrence_mass
                       - ap0.matched_transition_occurrence_mass, 0.01);

    tl206->level_feeds.clear();
    const BranchFeedingAudit ap1 = audit_for(db().nuclide("Tl206"), "Pb206");
    BOOST_CHECK(!ap1.scheme_valid);
    BOOST_CHECK(ap1.partial_scheme_rejected);

    tl206->level_feeds = {
        SandiaDecay::LevelFeed(0, 1.0f, 0.0f,
                               SandiaDecay::NoForbiddenness)
    };
    const BranchFeedingAudit ap2 = audit_for(db().nuclide("Tl206"), "Pb206");
    BOOST_CHECK(!ap2.scheme_valid);
    BOOST_CHECK(ap2.partial_scheme_rejected);

    // Zero-mass feed records and a sum outside the parser's 2e-4 closure
    // tolerance are invalid even if a caller constructs Transition by hand.
    BOOST_REQUIRE_GE(tl206_feeds.size(), 2u);
    tl206->level_feeds = tl206_feeds;
    const double removed = tl206->level_feeds.front().probability;
    const double remainder = 1.0 - removed;
    BOOST_REQUIRE_GT(remainder, 0.0);
    tl206->level_feeds.front().probability = 0.0f;
    for (std::size_t i = 1; i < tl206->level_feeds.size(); ++i)
        tl206->level_feeds[i].probability = static_cast<float>(
            tl206->level_feeds[i].probability / remainder);
    const BranchFeedingAudit ap3 = audit_for(db().nuclide("Tl206"), "Pb206");
    BOOST_CHECK(ap3.partial_scheme_rejected);

    tl206->level_feeds = tl206_feeds;
    for (SandiaDecay::LevelFeed& f : tl206->level_feeds)
        f.probability *= 1.0003f;
    const BranchFeedingAudit ap4 = audit_for(db().nuclide("Tl206"), "Pb206");
    BOOST_CHECK(ap4.partial_scheme_rejected);

    // LevelDag ignores non-descending host edges defensively. The adapter must
    // reject such an enriched candidate explicitly rather than replaying a graph
    // from which the edge silently disappeared.
    tl206->level_feeds = tl206_feeds;
    SandiaDecay::RadParticle* edge = nullptr;
    for (SandiaDecay::RadParticle& p : tl206->products)
        if (p.type == SandiaDecay::GammaParticle && p.from_level > p.to_level
            && p.to_level >= 0) {
            edge = &p;
            break;
        }
    BOOST_REQUIRE(edge != nullptr);
    const int saved_to_level = edge->to_level;
    edge->to_level = edge->from_level;
    const BranchFeedingAudit ap5 = audit_for(db().nuclide("Tl206"), "Pb206");
    edge->to_level = saved_to_level;
    BOOST_CHECK_GT(ap5.n_invalid_topology_edges, 0);
    BOOST_CHECK(!ap5.scheme_valid);

    tl206->level_feeds = std::move(tl206_feeds);

    // The descending-edge invariant is global, not merely part of low-coverage
    // exact-feed replay. U-238 is above 99% occurrence coverage, so this mutation
    // verifies that an invalid candidate edge cannot bypass the exact-feed gate.
    SandiaDecay::RadParticle* uedge = nullptr;
    for (SandiaDecay::RadParticle& p : u238->products)
        if (p.type == SandiaDecay::GammaParticle && p.from_level > p.to_level
            && p.to_level >= 0) {
            uedge = &p;
            break;
        }
    BOOST_REQUIRE(uedge != nullptr);
    const int u_saved_to = uedge->to_level;
    uedge->to_level = uedge->from_level;
    const BranchFeedingAudit au_bad = audit_for(db().nuclide("U238"), "Th234");
    uedge->to_level = u_saved_to;
    BOOST_CHECK_GT(au_bad.n_invalid_topology_edges, 0);
    BOOST_CHECK(!au_bad.scheme_valid);

    // In-120's exact starts replay its small matched graph, but the unknown
    // residual occurrence is 0.242 per branch. Exact replay cannot recover the
    // missing correlations, so the independent-residual approximation is not
    // admitted at this magnitude.
    const BranchFeedingAudit ai = audit_for(db().nuclide("In120"), "Sn120");
    BOOST_CHECK(ai.partial_exact_feed_compatible);
    BOOST_CHECK(ai.partial_scheme_rejected);
    BOOST_CHECK(!ai.scheme_valid);
    BOOST_CHECK_GT(ai.total_transition_occurrence_mass
                       - ai.matched_transition_occurrence_mass, 0.01);

    // A large multi-step graph can exceed 99% relative coverage while still
    // leaving more than 1% absolute occurrence uncorrelated. Sb-128 is the
    // regression case: 4.70118 / 4.74718 is above 99%, but its 0.046 residual
    // occurrence must still reject the graph.
    const BranchFeedingAudit asb = audit_for(db().nuclide("Sb128"), "Te128");
    BOOST_CHECK_GE(asb.matched_transition_occurrence_mass,
                   0.99 * asb.total_transition_occurrence_mass);
    BOOST_CHECK_GT(asb.total_transition_occurrence_mass
                       - asb.matched_transition_occurrence_mass, 0.01);
    BOOST_CHECK(asb.partial_scheme_rejected);
    BOOST_CHECK(!asb.scheme_valid);
}

BOOST_AUTO_TEST_CASE(unmatched_memberless_e0_counts_toward_partial_gate) {
    const SandiaDecay::Nuclide* parent = db().nuclide("U238");
    BOOST_REQUIRE(parent != nullptr);
    SandiaDecay::Transition* transition = nullptr;
    for (const SandiaDecay::Transition* t : parent->decaysToChildren)
        if (t && t->child && t->child->symbol == "Th234") {
            transition = const_cast<SandiaDecay::Transition*>(t);
            break;
        }
    BOOST_REQUIRE(transition != nullptr);
    SandiaDecay::RadParticle* unmatched = nullptr;
    for (SandiaDecay::RadParticle& p : transition->products)
        if (p.type == SandiaDecay::GammaParticle
            && (p.from_level < 0 || p.to_level < 0)) {
            unmatched = &p;
            break;
        }
    BOOST_REQUIRE(unmatched != nullptr);
    const float saved_intensity = unmatched->intensity;
    const bool saved_e0 = unmatched->e0_verified;
    struct RestoreProduct {
        SandiaDecay::RadParticle* p;
        float intensity;
        bool e0;
        ~RestoreProduct() { p->intensity = intensity; p->e0_verified = e0; }
    } restore{unmatched, saved_intensity, saved_e0};

    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    auto branch_audit = [&]() {
        std::vector<BranchFeedingAudit> audit;
        (void)build_cascades(parent, opts, &audit);
        for (const BranchFeedingAudit& a : audit)
            if (a.parent == "U238" && a.child == "Th234")
                return a;
        return BranchFeedingAudit{};
    };
    const BranchFeedingAudit before = branch_audit();

    // A memberless E0 contributes no gamma intensity. Give the unmatched record
    // a 2% transition occurrence: occurrence-based coverage must nevertheless
    // engage the gate and the absolute residual bound must reject the graph.
    unmatched->intensity = 0.02f;
    unmatched->e0_verified = true;
    const BranchFeedingAudit after = branch_audit();
    BOOST_CHECK_GT(after.total_transition_occurrence_mass
                       - before.total_transition_occurrence_mass, 0.01);
    BOOST_CHECK_LE(after.total_gamma_intensity,
                   before.total_gamma_intensity + 1e-12);
    BOOST_CHECK_GT(after.total_transition_occurrence_mass
                       - after.matched_transition_occurrence_mass, 0.01);
    BOOST_CHECK(after.partial_scheme_rejected);
    BOOST_CHECK(!after.scheme_valid);
}

// U235 was removed from this list on 2026-08-03: its flux inconsistency was
// four specific unsupported gamma records (26.55, 38.90, 41.40 keV removed;
// 64.35 keV corrected 0.0065 -> 1.7e-5), not an irrecoverable topology problem.
// With those fixed the branch's feed sum closes to 1.0013. Pa234 and Pa234m
// remain inconsistent and are still guarded here.
BOOST_AUTO_TEST_CASE(uranium_flux_inconsistent_graphs_remain_invalid) {
    CascadeOptions opts;
    opts.prompt_equilibrium = false;
    for (const char* name : {"Pa234", "Pa234m"}) {
        std::vector<BranchFeedingAudit> audit;
        (void)build_cascades(db().nuclide(name), opts, &audit);
        bool rejected = false;
        for (const BranchFeedingAudit& a : audit)
            rejected = rejected || a.rejected;
        BOOST_TEST_CONTEXT(name) { BOOST_CHECK(rejected); }
    }
}
