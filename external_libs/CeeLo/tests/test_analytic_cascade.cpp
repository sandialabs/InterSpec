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

// Deterministic unit tests for the fully-analytic cascade-summing path
// (compute_cascade_analytic) and the shared LevelDag combinatorics. A mock
// EfficiencyProvider returns fixed efficiencies per energy so every summing
// factor is checkable by hand (exact arithmetic, NO Monte Carlo). Analytic-vs-
// FullRealization MC agreement is covered in test_analytic_cascade_nuclides
// (against a frozen eps + FR reference), so it is not duplicated here.

#define BOOST_TEST_MODULE AnalyticCascade
#include <boost/test/included/unit_test.hpp>

#include "cascade/AnalyticCascade.h"
#include "cascade/CascadeTypes.h"
#include "cascade/LevelDag.h"
#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"

#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <vector>

using namespace ceelo;

namespace {

// Provider returning a fixed FEP/total efficiency per energy (nearest match
// within tol), so the DP arithmetic is exact and independent of any MC.
class MockProvider : public EfficiencyProvider {
public:
    void set(double E, double fep, double total) { m_[E] = {fep, total}; }
    double fep(double E) const override { return get(E).first; }
    double total(double E) const override { return get(E).second; }
    bool has(double E) const override { return find(E) != nullptr; }
private:
    const std::pair<double, double>* find(double E) const {
        const std::pair<double, double>* best = nullptr;
        double bd = 0.2;
        for (const auto& kv : m_) {
            const double d = std::abs(kv.first - E);
            if (d <= bd) { bd = d; best = &kv.second; }
        }
        return best;
    }
    std::pair<double, double> get(double E) const {
        const auto* p = find(E);
        return p ? *p : std::pair<double, double>{0.0, 0.0};
    }
    std::map<double, std::pair<double, double>> m_;
};

// Co-57-like scheme: member 0 = 136.47 (direct-to-ground, primary), member 1 =
// 122.06 (via 14.41 level), member 2 = 14.41. EC-K feed on the top level.
DecayCascade make_co57_like(double feed_ecK) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 26;
    dc.members.resize(3);
    dc.members[0].energy_keV = 136.47;
    dc.members[1].energy_keV = 122.06;
    dc.members[2].energy_keV = 14.41;
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    LevelScheme& ls = dc.level_scheme;
    ls.daughter_Z = 26;
    ls.levels.resize(3);
    ls.levels[2].feeding = 1.0;
    ls.levels[2].feed_ecK = feed_ecK;
    ls.levels[2].out.push_back({0, 0, 136.47, 0.115, 0.91, 0, 0, 0, 0});
    ls.levels[2].out.push_back({1, 1, 122.06, 0.885, 0.91, 0, 0, 0, 0});
    ls.levels[1].out.push_back({0, 2, 14.41, 1.000, 0.10, 0, 0, 0, 0});
    ls.valid = true;
    return dc;
}

// Co-60-like sequential scheme: level2 --(1173,gm0)--> level1 --(1332,gm1)--> gnd.
DecayCascade make_co60_like(double pg = 0.9998) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 28;
    dc.members.resize(2);
    dc.members[0].energy_keV = 1173.2;
    dc.members[1].energy_keV = 1332.5;
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    LevelScheme& ls = dc.level_scheme;
    ls.daughter_Z = 28;
    ls.levels.resize(3);
    ls.levels[2].feeding = 1.0;
    ls.levels[2].out.push_back({1, 0, 1173.2, 1.0, pg, 0, 0, 0, 0});
    ls.levels[1].out.push_back({0, 1, 1332.5, 1.0, pg, 0, 0, 0, 0});
    ls.valid = true;
    return dc;
}

// Sequential 3-gamma scheme: level3 --(60,gm0)--> level2 --(40,gm1)--> level1
// --(100,gm2)--> gnd, per-transition p_gamma configurable.
DecayCascade make_seq3(double p0, double p1, double p2) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 0;
    dc.members.resize(3);
    dc.members[0].energy_keV = 60.0;
    dc.members[1].energy_keV = 40.0;
    dc.members[2].energy_keV = 100.0;
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    LevelScheme& ls = dc.level_scheme;
    ls.levels.resize(4);
    ls.levels[3].feeding = 1.0;
    ls.levels[3].out.push_back({2, 0, 60.0, 1.0, p0, 0, 0, 0, 0});
    ls.levels[2].out.push_back({1, 1, 40.0, 1.0, p1, 0, 0, 0, 0});
    ls.levels[1].out.push_back({0, 2, 100.0, 1.0, p2, 0, 0, 0, 0});
    ls.valid = true;
    return dc;
}

DecayCascade make_residual_like() {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members = {
        {100.0, CascadeParticleType::Gamma, 0.75},
        {50.0, CascadeParticleType::Gamma, 0.20},
        {150.0, CascadeParticleType::Gamma, 0.10}
    };
    dc.level_scheme.levels.resize(2);
    dc.level_scheme.levels[1].feeding = 1.0;
    dc.level_scheme.levels[1].out.push_back(
        {0, 0, 100.0, 1.0, 0.75, 0.25, 0, 0, 0});
    dc.level_scheme.valid = true;
    dc.residual_transitions.push_back(
        {1, 50.0, 0.40, 0.20, 0.30, 0.10, 0.0, 0.0, 0.0, 0.0, 0.0});
    dc.residual_transitions.push_back(
        {2, 150.0, 0.90, 0.10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    return dc;
}

DecayCascade make_single_gamma_branch(double energy, bool valid,
                                      double branch_weight = 1.0) {
    DecayCascade dc;
    dc.branch_weight = branch_weight;
    dc.members = {{energy, CascadeParticleType::Gamma, 1.0}};
    if (valid) {
        dc.level_scheme.levels.resize(2);
        dc.level_scheme.levels[1].feeding = 1.0;
        dc.level_scheme.levels[1].out.push_back(
            {0, 0, energy, 1.0, 1.0, 0, 0, 0, 0});
        dc.level_scheme.valid = true;
    }
    return dc;
}

}  // namespace

BOOST_AUTO_TEST_CASE(invalid_branch_completeness_is_peak_specific) {
    DecayCascade primary100 = make_single_gamma_branch(100.0, true);
    DecayCascade primary120 = make_single_gamma_branch(120.0, true);
    DecayCascade invalid_feeder;
    invalid_feeder.branch_weight = 0.5;
    invalid_feeder.members = {
        {40.0, CascadeParticleType::Gamma, 0.40},
        {60.0, CascadeParticleType::Gamma, 0.30}
    };
    invalid_feeder.members[0].coincident.push_back({1, 0.75});
    invalid_feeder.members[1].coincident.push_back({0, 1.00});
    DecayCascade invalid_primary = make_single_gamma_branch(200.0, false, 0.8);

    MockProvider eff;
    for (double e : {40.0, 60.0, 100.0, 120.0, 200.0})
        eff.set(e, 0.1, 0.2);
    const std::vector<DecayCascade> cascades = {
        primary100, primary120, invalid_feeder, invalid_primary};
    const auto result = compute_cascade_analytic(
        cascades, {{100.0, 0.1}, {120.0, 0.1}, {200.0, 0.1}}, eff);
    BOOST_REQUIRE_EQUAL(result.size(), 3u);
    BOOST_REQUIRE(result[0].found && result[1].found && result[2].found);
    BOOST_CHECK(!result[0].summing_model_complete);  // omitted invalid 40+60 pair
    BOOST_CHECK(result[1].summing_model_complete);   // rejected branch cannot feed 120
    BOOST_CHECK(!result[2].summing_model_complete);  // primary itself is rejected
}

BOOST_AUTO_TEST_CASE(malformed_gamma_member_references_fail_closed) {
    MockProvider eff;
    eff.set(100.0, 0.1, 0.2);

    auto malformed = [](int gamma_member, CascadeParticleType member_type) {
        DecayCascade dc;
        dc.branch_weight = 1.0;
        dc.members = {{100.0, member_type, 1.0}};
        dc.level_scheme.levels.resize(2);
        dc.level_scheme.levels[1].feeding = 1.0;
        LevelOutTransition tr;
        tr.to_level = 0;
        tr.gamma_member = gamma_member;
        tr.gamma_keV = 100.0;
        tr.weight = 1.0;
        tr.p_gamma = 0.5;
        dc.level_scheme.levels[1].out.push_back(tr);
        dc.level_scheme.valid = true;
        return dc;
    };

    const std::vector<DecayCascade> bad = {
        malformed(-1, CascadeParticleType::Gamma),
        malformed(1, CascadeParticleType::Gamma),
        malformed(0, CascadeParticleType::Xray)
    };
    for (const DecayCascade& dc : bad) {
        BOOST_CHECK_THROW(
            compute_cascade_analytic({dc}, {{100.0, 0.1}}, eff),
            std::invalid_argument);
        BOOST_CHECK_THROW(
            EfficiencyCalculator::cascade_sum_pair_channels(
                {dc}, 0, 0, 100.0, 0.1),
            std::invalid_argument);

        EfficiencyCalculator calc;
        calc.set_fep_window_keV(kTestFepWindowKeV);
        CascadeConfig cfg;
        cfg.cascades = {dc};
        cfg.peaks = {{100.0, 0.1}};
        cfg.num_events = 1;
        cfg.num_threads = 1;
        for (CascadeMethod method : {CascadeMethod::Conditional,
                                     CascadeMethod::FullRealization}) {
            cfg.method = method;
            BOOST_CHECK_THROW(calc.compute_cascade(cfg), std::invalid_argument);
        }
    }
}

BOOST_AUTO_TEST_CASE(invalid_pair_diagnostic_uses_positron_outcome_mass) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members = {
        {489.0, CascadeParticleType::Gamma, 1.0},
        {510.998950, CascadeParticleType::Annih511, 0.9}
    };
    dc.weak_outcome_law.outcomes = {
        {WeakOutcomeKind::Other, -1, 0.8, 0.8, 0.0, 0.0, 0.0},
        {WeakOutcomeKind::Positron, -1, 0.2, 0.0, 0.0, 0.0, 1000.0}
    };
    BOOST_CHECK(!cascade_invalid_branch_can_feed_peak(dc, 999.998950, 1e-6));
    BOOST_CHECK(!cascade_invalid_branch_can_feed_peak(dc, 1021.997900, 1e-6));
    dc.weak_outcome_law.outcomes[1].selected_mass = 0.2;
    BOOST_CHECK(cascade_invalid_branch_can_feed_peak(dc, 999.998950, 1e-6));
    BOOST_CHECK(cascade_invalid_branch_can_feed_peak(dc, 1021.997900, 1e-6));
}

BOOST_AUTO_TEST_CASE(fallback_forest_projects_links_and_preserves_marginals) {
    DecayCascade dc;
    dc.members = {
        {100.0, CascadeParticleType::Gamma, 0.8},
        {200.0, CascadeParticleType::Gamma, 0.4},
        {300.0, CascadeParticleType::Gamma,
         std::numeric_limits<double>::quiet_NaN()}
    };
    // Directional implications disagree: 0->1 gives raw J=.76 -> .40,
    // 1->0 gives raw J=.04 -> Frechet lower bound .20. The latter has the
    // stronger departure from independence (.12), so it deterministically wins.
    dc.members[0].coincident.push_back({1, 0.95});
    dc.members[1].coincident.push_back({0, 0.10});
    dc.members[0].coincident.push_back(
        {2, std::numeric_limits<double>::quiet_NaN()});
    const auto forest = build_cascade_fallback_forest(dc);
    BOOST_REQUIRE_EQUAL(forest.size(), 3u);
    BOOST_CHECK_EQUAL(forest[0].member, 0u);
    BOOST_CHECK_EQUAL(forest[1].member, 1u);
    BOOST_CHECK_EQUAL(forest[1].parent_member, 0);
    BOOST_CHECK_CLOSE(forest[1].p_if_parent_emitted, 0.25, 1e-10);
    BOOST_CHECK_CLOSE(forest[1].p_if_parent_absent, 1.0, 1e-10);
    BOOST_CHECK_EQUAL(forest[2].member, 2u);
    BOOST_CHECK_SMALL(forest[2].marginal, 1e-15);

    std::mt19937_64 rng(0x51a7eULL);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    constexpr int n = 400000;
    int c0 = 0, c1 = 0, c01 = 0;
    std::vector<char> emitted;
    for (int k = 0; k < n; ++k) {
        sample_cascade_fallback_forest(
            forest, -1, [&]() { return u(rng); }, emitted);
        c0 += emitted[0]; c1 += emitted[1];
        c01 += emitted[0] && emitted[1];
        BOOST_CHECK(!emitted[2]);
    }
    auto z = [&](int count, double p) {
        return (count / static_cast<double>(n) - p)
             / std::sqrt(p * (1.0 - p) / n);
    };
    BOOST_CHECK_LT(std::abs(z(c0, 0.8)), 5.0);
    BOOST_CHECK_LT(std::abs(z(c1, 0.4)), 5.0);
    BOOST_CHECK_LT(std::abs(z(c01, 0.2)), 5.0);
}

BOOST_AUTO_TEST_CASE(fallback_forest_parent_precedes_lower_index_child) {
    DecayCascade dc;
    dc.members = {
        {100.0, CascadeParticleType::Gamma, 0.4},  // child has lower member index
        {200.0, CascadeParticleType::Gamma, 0.9}
    };
    dc.members[1].coincident.push_back({0, 0.4 / 0.9});
    const auto forest = build_cascade_fallback_forest(dc);
    BOOST_REQUIRE_EQUAL(forest.size(), 2u);
    BOOST_CHECK_EQUAL(forest[0].member, 1u);
    BOOST_CHECK_EQUAL(forest[1].member, 0u);
    BOOST_CHECK_EQUAL(forest[1].parent_member, 1);
}

BOOST_AUTO_TEST_CASE(fallback_forest_weak_annihilation_forces_before_descendants) {
    DecayCascade dc;
    dc.members = {
        {510.998950, CascadeParticleType::Annih511, 0.1},
        {100.0, CascadeParticleType::Gamma, 0.4},
        {200.0, CascadeParticleType::Gamma, 0.9}
    };
    dc.weak_outcome_law.outcomes = {
        {WeakOutcomeKind::Other, -1, 0.2, 0.2, 0, 0, 0},
        {WeakOutcomeKind::Positron, -1, 0.8, 0.8, 0, 0, 1000.0}
    };
    // The 200-keV member is earlier than annihilation, but forced Annih511 must
    // not select it as a parent. Annihilation may in turn parent the later
    // 100-keV member, with J=.8*.5=.4.
    dc.members[2].coincident.push_back({0, 1.0});
    dc.members[0].coincident.push_back({1, 0.5});
    const auto forest = build_cascade_fallback_forest(dc);
    BOOST_REQUIRE_EQUAL(forest.size(), 3u);
    BOOST_CHECK_EQUAL(forest[0].member, 2u);
    BOOST_CHECK_EQUAL(forest[1].member, 0u);
    BOOST_CHECK(forest[1].forced_annihilation);
    BOOST_CHECK_EQUAL(forest[1].parent_member, -1);
    BOOST_CHECK_EQUAL(forest[2].member, 1u);
    BOOST_CHECK_EQUAL(forest[2].parent_member, 0);

    std::mt19937_64 rng(0xface511ULL);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    constexpr int n = 400000;
    int ca = 0, c1 = 0, c2 = 0, ca1 = 0;
    std::uint64_t draws = 0;
    std::vector<char> emitted;
    for (int k = 0; k < n; ++k) {
        const bool positron = u(rng) < 0.8;
        sample_cascade_fallback_forest(forest, positron ? 1 : 0, [&]() {
            ++draws;
            return u(rng);
        }, emitted);
        ca += emitted[0]; c1 += emitted[1]; c2 += emitted[2];
        ca1 += emitted[0] && emitted[1];
    }
    BOOST_CHECK_EQUAL(draws, static_cast<std::uint64_t>(n) * forest.size());
    auto z = [&](int count, double p) {
        return (count / static_cast<double>(n) - p)
             / std::sqrt(p * (1.0 - p) / n);
    };
    BOOST_CHECK_LT(std::abs(z(ca, 0.8)), 5.0);
    BOOST_CHECK_LT(std::abs(z(c1, 0.4)), 5.0);
    BOOST_CHECK_LT(std::abs(z(c2, 0.9)), 5.0);
    BOOST_CHECK_LT(std::abs(z(ca1, 0.4)), 5.0);
}

BOOST_AUTO_TEST_CASE(residual_transition_conditionals_are_exclusive) {
    const DecayCascade dc = make_residual_like();
    const auto matched_pmate = EfficiencyCalculator::cascade_level_pmate(dc, 0);
    BOOST_REQUIRE_EQUAL(matched_pmate.size(), dc.members.size());
    BOOST_CHECK_CLOSE(matched_pmate[1], 0.20, 1e-9);
    BOOST_CHECK_CLOSE(matched_pmate[2], 0.10, 1e-9);

    const auto residual_pmate = EfficiencyCalculator::cascade_level_pmate(dc, 1);
    BOOST_REQUIRE_EQUAL(residual_pmate.size(), dc.members.size());
    BOOST_CHECK_CLOSE(residual_pmate[0], 0.75, 1e-9);
    BOOST_CHECK_CLOSE(residual_pmate[2], 0.10, 1e-9);

    const auto matched_vac = EfficiencyCalculator::cascade_level_vacancies(dc, 0);
    int residual_group = -1;
    double pK = 0.0, pL = 0.0;
    for (const auto& v : matched_vac) {
        if (v.gamma_member != 1) continue;
        if (residual_group < 0) residual_group = v.exclusive_group;
        BOOST_CHECK_EQUAL(v.exclusive_group, residual_group);
        if (v.is_L) pL += v.prob;
        else if (v.produces_vacancy) pK += v.prob;
    }
    BOOST_CHECK_CLOSE(pK, 0.30, 1e-9);
    BOOST_CHECK_CLOSE(pL, 0.10, 1e-9);

    const auto own_vac = EfficiencyCalculator::cascade_level_vacancies(dc, 1);
    for (const auto& v : own_vac)
        BOOST_CHECK_NE(v.gamma_member, 1);
}

BOOST_AUTO_TEST_CASE(residual_transition_analytic_cout_and_sumin) {
    const DecayCascade dc = make_residual_like();
    MockProvider eff;
    eff.set(50.0, 0.20, 0.40);
    eff.set(100.0, 0.20, 0.30);
    // Keep the 150-keV primary from removing the 100+50 sum-fed event; the
    // analytic path deliberately counts the full pair joint before survival.
    eff.set(150.0, 0.20, 0.0);

    const auto matched = compute_cascade_analytic({dc}, {{100.0, 0.1}}, eff);
    BOOST_REQUIRE(matched[0].found);
    BOOST_CHECK_CLOSE(matched[0].c_out, 1.0 - 0.20 * 0.40, 1e-8);

    const auto residual = compute_cascade_analytic({dc}, {{50.0, 0.1}}, eff);
    BOOST_REQUIRE(residual[0].found);
    BOOST_CHECK_CLOSE(residual[0].c_out,
        (1.0 - 0.75 * 0.30) * (1.0 - 0.10 * 0.0), 1e-8);

    const auto sum = compute_cascade_analytic({dc}, {{150.0, 0.1}}, eff);
    BOOST_REQUIRE(sum[0].found);
    const double gain = 0.10 * sum[0].c_in * sum[0].eff_no_summing;
    BOOST_CHECK_CLOSE(gain, 0.75 * 0.20 * 0.20 * 0.20, 1e-8);

    const auto channels = EfficiencyCalculator::cascade_sum_pair_channels(
        {dc}, 0, 2, 150.0, 0.1);
    bool saw_residual_pair = false;
    for (const auto& ch : channels)
        if ((ch.member_a == 0 && ch.member_b == 1)
            || (ch.member_a == 1 && ch.member_b == 0)) {
            saw_residual_pair = true;
            // Conditional subtracts the direct primary overlap: residual 150
            // independently fires in 10% of otherwise sum-fed histories.
            BOOST_CHECK_CLOSE(ch.weight,
                0.75 * 0.20 * (1.0 - 0.10) / 0.10, 1e-8);
        }
    BOOST_CHECK(saw_residual_pair);
}

BOOST_AUTO_TEST_CASE(legacy_ec_ground_mass_is_absolute_in_residual_sumin) {
    const std::vector<detail::XrayLine> kl = detail::k_lines(26);
    BOOST_REQUIRE(!kl.empty());
    const double k_energy = kl.front().energy;
    double k_weight = 0.0;
    for (const detail::XrayLine& line : kl)
        if (std::abs(line.energy - k_energy) < 1e-9)
            k_weight += line.weight;

    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 26;
    const double residual_energy = 50.0;
    const double peak_energy = residual_energy + k_energy;
    dc.members = {
        {residual_energy, CascadeParticleType::Gamma, 1.0},
        {peak_energy, CascadeParticleType::Gamma, 0.10}
    };
    dc.level_scheme.daughter_Z = 26;
    dc.level_scheme.levels.resize(2);
    dc.level_scheme.levels[1].feeding = 1.0;
    dc.level_scheme.levels[0].feed_ecK = 1.0;
    dc.level_scheme.entry_probability = 0.20;
    dc.level_scheme.valid = true;
    dc.residual_transitions.push_back(
        {0, residual_energy, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    dc.residual_transitions.push_back(
        {1, peak_energy, 0.90, 0.10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    MockProvider eff;
    eff.set(residual_energy, 0.20, 0.0);
    eff.set(k_energy, 0.30, 0.0);
    eff.set(peak_energy, 0.10, 0.0);
    const auto result = compute_cascade_analytic(
        {dc}, {{peak_energy, 1e-6}}, eff);
    BOOST_REQUIRE(result[0].found);
    const double gain = 0.10 * result[0].c_in * result[0].eff_no_summing;
    // 80% of decays bypass the normalized DAG and land on ground. Their
    // ground-only K vacancy remains coincident with the independent residual.
    const double expected = 0.80 * k_weight * 0.20 * 0.30;
    BOOST_CHECK_CLOSE(gain, expected, 1e-8);
}

BOOST_AUTO_TEST_CASE(legacy_ec_pair_survival_is_shell_start_conditioned) {
    const std::vector<detail::XrayLine> kl = detail::k_lines(26);
    BOOST_REQUIRE(!kl.empty());
    const double k_energy = kl.front().energy;
    double k_weight = 0.0;
    for (const detail::XrayLine& line : kl)
        if (std::abs(line.energy - k_energy) < 1e-9)
            k_weight += line.weight;

    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 26;
    const double residual_energy = 50.0;
    const double peak_energy = residual_energy + k_energy;
    dc.members = {
        {residual_energy, CascadeParticleType::Gamma, 1.0},
        {peak_energy, CascadeParticleType::Gamma, 0.10},
        {20.0, CascadeParticleType::Gamma, 1.0}
    };
    dc.level_scheme.daughter_Z = 26;
    dc.level_scheme.levels.resize(3);
    dc.level_scheme.levels[1].feeding = 0.50;
    dc.level_scheme.levels[2].feeding = 0.50;
    dc.level_scheme.levels[2].feed_ecK = 1.0;
    dc.level_scheme.levels[2].out.push_back(
        {1, 2, 20.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0});
    dc.level_scheme.entry_probability = 1.0;
    dc.level_scheme.valid = true;
    dc.residual_transitions.push_back(
        {0, residual_energy, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    dc.residual_transitions.push_back(
        {1, peak_energy, 0.90, 0.10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    MockProvider eff;
    eff.set(residual_energy, 0.20, 0.0);
    eff.set(k_energy, 0.30, 0.0);
    eff.set(peak_energy, 0.10, 0.0);
    eff.set(20.0, 0.0, 0.40);
    const auto result = compute_cascade_analytic(
        {dc}, {{peak_energy, 1e-6}}, eff);
    BOOST_REQUIRE(result[0].found);
    const double gain = 0.10 * result[0].c_in * result[0].eff_no_summing;
    // Only the level-2 half of the starts carries a K vacancy. Conditioning on
    // observing that K x-ray therefore forces the 20-keV transition, whose 40%
    // total efficiency leaves a 0.60 pair-survival factor. Averaging the level-1
    // non-K starts would incorrectly give 0.80.
    const double expected = 0.50 * k_weight * 0.20 * 0.30 * 0.60;
    BOOST_CHECK_CLOSE(gain, expected, 1e-8);
}

// ============================= C_out (summing-out) =========================

// T1: Co-57 136 keV primary -> C_out = 1 exactly (122 and 14 are mutually
// exclusive with the direct 136 de-excitation, so nothing sums out).
BOOST_AUTO_TEST_CASE(t1_co57_exclusive_cout_is_unity) {
    MockProvider eff;
    eff.set(136.47, 0.05, 0.15);
    eff.set(122.06, 0.05, 0.10);
    eff.set(14.41, 0.05, 0.10);
    const auto r = compute_cascade_analytic({make_co57_like(0.0)},
                                            {{136.47, 1.5}}, eff);
    BOOST_REQUIRE_EQUAL(r.size(), 1u);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_CLOSE(r[0].c_out, 1.0, 1e-9);
    // The 122.06 + 14.41 pair sums into the 136.47 window (sum-fed summing-IN),
    // so c_in > 0: w_in = pair_joint(122,14)*pg*pg / (pass(136)*pg136)
    //   = 0.885*0.91*0.10 / (0.115*0.91) = 0.76957; gain = w_in*eps_fep(122)*eps_fep(14).
    const double w_in = 0.885 * 0.91 * 0.10 / (0.115 * 0.91);
    const double gain = w_in * 0.05 * 0.05;
    BOOST_CHECK_CLOSE(r[0].c_in * r[0].eff_no_summing, gain, 1e-4);
}

// T1 variant: Co-57 122 keV primary -> only the coincident 14.41 sums out,
// C_out = 1 - P(14.41 gamma|122)*eps_tot(14.41) = 1 - 0.10*0.10 = 0.99.
BOOST_AUTO_TEST_CASE(t1b_co57_122_cout) {
    MockProvider eff;
    eff.set(136.47, 0.05, 0.15);
    eff.set(122.06, 0.05, 0.10);
    eff.set(14.41, 0.05, 0.10);
    const auto r = compute_cascade_analytic({make_co57_like(0.0)},
                                            {{122.06, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_CLOSE(r[0].c_out, 0.99, 1e-6);
}

// T2: Co-60 1173 primary, sequential -> C_out = 1 - 0.9998*eps_tot(1332).
BOOST_AUTO_TEST_CASE(t2_co60_sequential_cout) {
    MockProvider eff;
    const double tau = 0.30;
    eff.set(1173.2, 0.04, 0.12);
    eff.set(1332.5, 0.04, tau);
    const auto r = compute_cascade_analytic({make_co60_like()}, {{1173.2, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_CLOSE(r[0].c_out, 1.0 - 0.9998 * tau, 1e-6);
}

// T3: 3-gamma sequential, primary = 60 (gm0). The two coincident gammas 40,100
// are independent given the single path, so C_out = (1-t40)(1-t100) EXACTLY
// (no spurious +t40*t100 cross term). p_gamma = 1 throughout.
BOOST_AUTO_TEST_CASE(t3_three_gamma_product_is_exact) {
    MockProvider eff;
    const double t40 = 0.20, t100 = 0.30;
    eff.set(60.0, 0.05, 0.10);
    eff.set(40.0, 0.05, t40);
    eff.set(100.0, 0.05, t100);
    const auto r = compute_cascade_analytic({make_seq3(1.0, 1.0, 1.0)},
                                            {{60.0, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_CLOSE(r[0].c_out, (1.0 - t40) * (1.0 - t100), 1e-9);
}

// ============================= C_in (summing-in) ===========================

// T4: the sum-fed pair (60,40) -> 100 window. The ANALYTIC weight does NOT
// subtract the triple overlap the MC direct stream removes: w_in =
// pair_joint*pg60*pg40 / (pass(100)*pg100) = 1*0.9*0.9/(1*0.5) = 1.62 (the MC
// channel weight for the same scheme is the SUBTRACTED 0.81). So the summing-in
// gain = 1.62 * eps_fep(60) * eps_fep(40).
BOOST_AUTO_TEST_CASE(t4_suminc_no_triple_subtraction) {
    MockProvider eff;
    const double f60 = 0.10, f40 = 0.12, f100 = 0.20;
    eff.set(60.0, f60, 0.3);
    eff.set(40.0, f40, 0.3);
    eff.set(100.0, f100, 0.0);  // eps_tot(100)=0 -> coincident-survival = 1 (isolate w_in)
    const auto r = compute_cascade_analytic({make_seq3(0.9, 0.9, 0.5)},
                                            {{100.0, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    const double expected_gain = 1.62 * f60 * f40;
    const double gain = r[0].c_in * r[0].eff_no_summing;  // c_in = gain/eff_no
    BOOST_CHECK_CLOSE(gain, expected_gain, 1e-6);
}

// T4c: coincident-survival factor on a sum-fed channel. In the sequential 60->40
// ->100 scheme, whenever 60+40 sum into the 100 window the 100 gamma is ALSO
// emitted (p_gamma=0.5) and if it deposits (eps_tot=0.4) the event is knocked out.
// So gain = 1.62*f60*f40 * (1 - 0.5*0.4). This is the near-field summing-out of
// the summing-in channel (Co-57 122+14->136 + coincident Fe Kalpha in real data).
BOOST_AUTO_TEST_CASE(t4c_suminc_coincident_survival) {
    MockProvider eff;
    const double f60 = 0.10, f40 = 0.12;
    eff.set(60.0, f60, 0.3);
    eff.set(40.0, f40, 0.3);
    eff.set(100.0, 0.20, 0.4);  // the coincident 100 gamma: p_gamma=0.5, eps_tot=0.4
    const auto r = compute_cascade_analytic({make_seq3(0.9, 0.9, 0.5)},
                                            {{100.0, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    const double expected_gain = 1.62 * f60 * f40 * (1.0 - 0.5 * 0.4);
    const double gain = r[0].c_in * r[0].eff_no_summing;
    BOOST_CHECK_CLOSE(gain, expected_gain, 1e-6);
}

// T4b: same-vacancy exclusivity. A single transition that converts in K produces
// TWO K x-ray occurrences (Kalpha, Kbeta); they come from the SAME vacancy and
// are mutually exclusive, so a pair of them must NEVER contribute to summing-in.
// Build a scheme with one converting transition whose two K lines could sum into
// a window; assert no gain from the same-vacancy pair.
BOOST_AUTO_TEST_CASE(t4b_same_vacancy_pair_excluded) {
    // Fe (Z=26): Kalpha ~6.40, Kbeta ~7.06 -> sum ~13.46 keV. One transition
    // converting in K (p_icK=1) is the only emission besides its (suppressed)
    // gamma. Peak at the Kalpha+Kbeta sum must get zero sum-in (same vacancy).
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 26;
    dc.members.resize(1);
    dc.members[0].energy_keV = 500.0;  // a gamma, mostly converted
    dc.members[0].type = CascadeParticleType::Gamma;
    LevelScheme& ls = dc.level_scheme;
    ls.daughter_Z = 26;
    ls.levels.resize(2);
    ls.levels[1].feeding = 1.0;
    ls.levels[1].out.push_back({0, 0, 500.0, 1.0, /*pg*/0.0, /*icK*/1.0, 0, 0, 0});
    ls.valid = true;

    MockProvider eff;
    eff.set(6.40, 0.5, 0.5);
    eff.set(7.06, 0.5, 0.5);
    eff.set(13.46, 0.3, 0.3);  // the (spurious) sum peak
    const auto r = compute_cascade_analytic({dc}, {{13.46, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found == false || r[0].c_in == 0.0);
    if (r[0].found) BOOST_CHECK_SMALL(r[0].c_in, 1e-12);
}

// ============================= T6 first-order ==============================

// T6: master cross-check. As eps -> 0, (1 - C_out)/eps must equal the sum of the
// per-coincident-member deposit coefficients that the exact MC combinatorics
// produce. On Co-60 that is P(1332|1173) = 0.9998 (cascade_level_pmate).
BOOST_AUTO_TEST_CASE(t6_first_order_matches_pmate) {
    const double eps = 1e-5;
    MockProvider eff;
    eff.set(1173.2, 0.04, eps);
    eff.set(1332.5, 0.04, eps);
    const auto r = compute_cascade_analytic({make_co60_like()}, {{1173.2, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    const double coeff = (1.0 - r[0].c_out) / eps;
    // exact pmate from the shared DP
    const auto pmate = EfficiencyCalculator::cascade_level_pmate(make_co60_like(), 0);
    BOOST_REQUIRE_GE(pmate.size(), 2u);
    BOOST_CHECK_CLOSE(coeff, pmate[1], 1e-3);
    BOOST_CHECK_CLOSE(coeff, 0.9998, 1e-3);
}

// ============================= T5 LevelDag golden ==========================

// T5: the shared LevelDag joints reproduce hand values, and n_subset_joint
// subsumes pair_joint / triple_joint.
BOOST_AUTO_TEST_CASE(t5_leveldag_joints_golden) {
    const DecayCascade co57 = make_co57_like(0.88);
    LevelDag dag(co57.level_scheme);
    BOOST_REQUIRE(dag.valid);
    const int t136 = dag.transition_of(0);
    const int t122 = dag.transition_of(1);
    const int t14 = dag.transition_of(2);
    // 136 is exclusive with 122 and 14 (same top level path split): pair_joint=0.
    BOOST_CHECK_SMALL(dag.pair_joint(t136, t122), 1e-12);
    BOOST_CHECK_SMALL(dag.pair_joint(t136, t14), 1e-12);
    // 122 and 14 are sequential: joint = pass(122)=0.885.
    BOOST_CHECK_CLOSE(dag.pair_joint(t122, t14), 0.885, 1e-6);
    // n_subset_joint subsumes pair / triple.
    BOOST_CHECK_CLOSE(dag.n_subset_joint({t122, t14}), dag.pair_joint(t122, t14), 1e-9);
    BOOST_CHECK_CLOSE(dag.n_subset_joint({t136, t122, t14}),
                      dag.triple_joint(t136, t122, t14), 1e-9);
    // triple with the exclusive 136 is 0.
    BOOST_CHECK_SMALL(dag.n_subset_joint({t136, t122, t14}), 1e-12);

    // Co-60 sequential: n_subset_joint({1173,1332}) = pass(1173) = 1.0.
    const DecayCascade co60 = make_co60_like();
    LevelDag d2(co60.level_scheme);
    BOOST_CHECK_CLOSE(d2.n_subset_joint({d2.transition_of(0), d2.transition_of(1)}),
                      1.0, 1e-9);
}

// T8: coincident annihilation-511 pair summing-out. A gamma primary coincident
// with a 511 PAIR (one Annih511 member = two back-to-back photons). For a single-
// sided detector at most ONE 511 can deposit, so P(>=1 deposits) = 2*eps_tot(511)
// (not the independent 2eps-eps^2). C_out = 1 - intensity * 2*eps_tot(511).
BOOST_AUTO_TEST_CASE(t8_annih511_backtoback_survival) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 10;
    dc.members.resize(2);
    dc.members[0].energy_keV = 1274.5;  // the gamma (primary)
    dc.members[0].type = CascadeParticleType::Gamma;
    dc.members[1].energy_keV = 511.0;   // annihilation pair member
    dc.members[1].type = CascadeParticleType::Annih511;
    dc.members[1].intensity = 1.0;      // always coincident
    LevelScheme& ls = dc.level_scheme;
    ls.levels.resize(2);
    ls.levels[1].feeding = 1.0;
    ls.levels[1].out.push_back({0, 0, 1274.5, 1.0, 1.0, 0, 0, 0, 0});
    ls.valid = true;

    MockProvider eff;
    const double et511 = 0.20;
    eff.set(1274.5, 0.03, 0.10);
    eff.set(511.0, 0.05, et511);
    const auto r = compute_cascade_analytic({dc}, {{1274.5, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_CLOSE(r[0].c_out, 1.0 - 1.0 * (2.0 * et511), 1e-6);  // = 0.60
}

// T9: overlapped window (fitted-peak-area SF). Two emitted lines 80.5 and 81.0
// fall in one 1.5 keV window; the no-summing efficiency is the emission-weighted
// mean FEP efficiency, eff_no = (A0·ε0 + A1·ε1)/(A0+A1). Here two levels feed
// ground directly (mutually exclusive per decay), no coincidences -> SF = 1.
BOOST_AUTO_TEST_CASE(t9_overlapped_window_area_weighted) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members.resize(2);
    dc.members[0].energy_keV = 80.5;
    dc.members[1].energy_keV = 81.0;
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    LevelScheme& ls = dc.level_scheme;
    ls.levels.resize(3);                       // 0 = ground
    ls.levels[1].feeding = 1.0;                // 80.5 level
    ls.levels[1].out.push_back({0, 0, 80.5, 1.0, 1.0, 0, 0, 0, 0});
    ls.levels[2].feeding = 0.5;                // 81.0 level
    ls.levels[2].out.push_back({0, 1, 81.0, 1.0, 1.0, 0, 0, 0, 0});
    ls.valid = true;

    MockProvider eff;
    eff.set(80.5, 0.10, 0.20);
    eff.set(81.0, 0.20, 0.30);
    const auto r = compute_cascade_analytic({dc}, {{81.0, 1.5}}, eff);
    BOOST_REQUIRE(r[0].found);
    // A(80.5) = 1.0/1.5, A(81.0) = 0.5/1.5 (per-decay feeding fractions).
    const double A0 = 1.0 / 1.5, A1 = 0.5 / 1.5;
    const double eff_no = (A0 * 0.10 + A1 * 0.20) / (A0 + A1);
    BOOST_CHECK_CLOSE(r[0].eff_no_summing, eff_no, 1e-6);
    BOOST_CHECK_CLOSE(r[0].summing_factor, 1.0, 1e-9);  // no coincidences
}

// ===================== T7: angular correlation (deterministic) =============

// The W(theta) angular correction multiplies a correlated SUMMING-IN pair by the
// collinear limit g = W(0) = 1 + a2 + a4 (both photons must be FEP-detected, so
// their mutual angle is small). On the Co-57 122+14->136 sum-fed channel with a
// correlated 122-14 link, c_in(on) = c_in(off) * (1 + a2 + a4); c_out is
// unaffected (136 is exclusive with 122/14, and summing-out gets g=1 anyway).
BOOST_AUTO_TEST_CASE(t7_angular_w0_suminc) {
    DecayCascade dc = make_co57_like(/*feed_ecK=*/0.0);
    const double a2 = 0.10, a4 = 0.02;
    dc.members[1].coincident.push_back({2, 1.0, a2, a4, true});  // 122 <-> 14
    dc.members[2].coincident.push_back({1, 1.0, a2, a4, true});

    MockProvider eff;
    eff.set(136.47, 0.05, 0.15);
    eff.set(122.06, 0.05, 0.10);
    eff.set(14.41, 0.05, 0.10);

    AnalyticCascadeOptions off; off.apply_angular_correlation = false;
    AnalyticCascadeOptions on;  on.apply_angular_correlation = true;
    const auto r_off = compute_cascade_analytic({dc}, {{136.47, 1.5}}, eff, off);
    const auto r_on  = compute_cascade_analytic({dc}, {{136.47, 1.5}}, eff, on);

    BOOST_CHECK_GT(r_off[0].c_in, 0.0);
    BOOST_CHECK_CLOSE(r_on[0].c_in, r_off[0].c_in * (1.0 + a2 + a4), 1e-6);
    BOOST_CHECK_CLOSE(r_on[0].c_out, r_off[0].c_out, 1e-9);
}

// Exact weak modes are selected before the daughter walk.  Here only EC feeds
// the gamma-emitting level; beta+ feeds ground.  Conditioning on the gamma must
// therefore make P(annihilation | gamma)=0, not the 80% unconditional beta+
// marginal used by the former independent-member model.
BOOST_AUTO_TEST_CASE(t10_weak_mode_posterior_removes_impossible_annihilation) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members.resize(2);
    dc.members[0] = {100.0, CascadeParticleType::Gamma, 0.2};
    dc.members[1] = {511.0, CascadeParticleType::Annih511, 0.8};
    dc.level_scheme.levels.resize(2);
    dc.level_scheme.levels[1].feeding = 0.2;
    dc.level_scheme.levels[1].out.push_back(
        {0, 0, 100.0, 1.0, 1.0, 0, 0, 0, 0});
    dc.level_scheme.entry_probability = 0.2;
    dc.level_scheme.valid = true;
    dc.weak_outcome_law.outcomes = {
        {WeakOutcomeKind::ElectronCapture, 1, 0.2, 0.2, 1.0, 0.0, 0.0},
        {WeakOutcomeKind::Positron,        0, 0.8, 0.8, 0.0, 0.0, 500.0}
    };
    dc.weak_outcome_law.raw_sum = 1.0;

    MockProvider eff;
    eff.set(100.0, 0.1, 0.2);
    eff.set(511.0, 0.1, 0.25);
    const auto r = compute_cascade_analytic({dc}, {{100.0, 0.2}}, eff);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_CLOSE(r[0].c_out, 1.0, 1e-9);

    const auto pm = EfficiencyCalculator::cascade_level_pmate(dc, 0);
    BOOST_REQUIRE_EQUAL(pm.size(), 2u);
    BOOST_CHECK_SMALL(pm[1], 1e-12);
    const auto vac = EfficiencyCalculator::cascade_level_vacancies(dc, 0);
    BOOST_REQUIRE_EQUAL(vac.size(), 1u);
    BOOST_CHECK_CLOSE(vac[0].prob, 1.0, 1e-9);
    BOOST_CHECK_EQUAL(vac[0].exclusive_group, 0);
}

// Conversely, when only beta+ reaches the gamma, conditioning must make the
// annihilation pair certain.  For a single-sided detector C_out=1-2 eps_tot511.
BOOST_AUTO_TEST_CASE(t11_weak_mode_posterior_makes_annihilation_certain) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members.resize(2);
    dc.members[0] = {100.0, CascadeParticleType::Gamma, 0.8};
    dc.members[1] = {511.0, CascadeParticleType::Annih511, 0.8};
    dc.level_scheme.levels.resize(2);
    dc.level_scheme.levels[1].feeding = 0.8;
    dc.level_scheme.levels[1].out.push_back(
        {0, 0, 100.0, 1.0, 1.0, 0, 0, 0, 0});
    dc.level_scheme.entry_probability = 0.8;
    dc.level_scheme.valid = true;
    dc.weak_outcome_law.outcomes = {
        {WeakOutcomeKind::ElectronCapture, 0, 0.2, 0.2, 1.0, 0.0, 0.0},
        {WeakOutcomeKind::Positron,        1, 0.8, 0.8, 0.0, 0.0, 500.0}
    };
    dc.weak_outcome_law.raw_sum = 1.0;

    MockProvider eff;
    eff.set(100.0, 0.1, 0.2);
    eff.set(511.0, 0.1, 0.25);
    const auto r = compute_cascade_analytic({dc}, {{100.0, 0.2}}, eff);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_CLOSE(r[0].c_out, 0.5, 1e-9);
    const auto pm = EfficiencyCalculator::cascade_level_pmate(dc, 0);
    BOOST_REQUIRE_EQUAL(pm.size(), 2u);
    BOOST_CHECK_CLOSE(pm[1], 1.0, 1e-9);
    BOOST_CHECK(EfficiencyCalculator::cascade_level_vacancies(dc, 0).empty());
}

// A selected weak outcome must condition every upstream nuclear partner and IC
// vacancy, not merely the annihilation member.  EC starts at the primary level;
// beta+ starts one level higher and must traverse a fully-converted transition.
BOOST_AUTO_TEST_CASE(t12_outcome_specific_upstream_conversion) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members = {{100.0, CascadeParticleType::Gamma, 1.0}};
    dc.level_scheme.levels.resize(3);
    dc.level_scheme.levels[1].feeding = 0.5;
    dc.level_scheme.levels[2].feeding = 0.5;
    dc.level_scheme.levels[1].out.push_back(
        {0, 0, 100.0, 1.0, 1.0, 0, 0, 0, 0});
    dc.level_scheme.levels[2].out.push_back(
        {1, -1, 40.0, 1.0, 0.0, 1.0, 0, 0, 0});
    dc.level_scheme.entry_probability = 1.0;
    dc.level_scheme.valid = true;
    dc.weak_outcome_law.outcomes = {
        {WeakOutcomeKind::ElectronCapture, 1, 0.5, 0.5, 0.8, 0.0, 0.0},
        {WeakOutcomeKind::Positron,        2, 0.5, 0.5, 0.0, 0.0, 500.0}
    };
    dc.weak_outcome_law.raw_sum = 1.0;

    const auto vac = EfficiencyCalculator::cascade_level_vacancies(dc, 0);
    bool found_ic = false;
    for (const auto& v : vac) {
        if (v.exclusive_group <= 0) continue;
        found_ic = true;
        BOOST_CHECK_CLOSE(v.prob, 0.5, 1e-9);
        BOOST_REQUIRE_EQUAL(v.weak_outcome_prob.size(), 2u);
        BOOST_CHECK_SMALL(v.weak_outcome_prob[0], 1e-12);
        BOOST_CHECK_CLOSE(v.weak_outcome_prob[1], 1.0, 1e-9);
    }
    BOOST_CHECK(found_ic);
}

// The supplemental pair stream represents pair gammas AND NOT the primary
// gamma.  Here the pair is certain, while the downstream primary transition is
// 50% gamma / 50% K conversion.  Its channel must never refire the excluded
// primary; instead the primary transition's K vacancy is certain in that stream.
BOOST_AUTO_TEST_CASE(t13_sum_pair_conditions_on_primary_not_emitted) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members = {
        {40.0,  CascadeParticleType::Gamma, 1.0},
        {60.0,  CascadeParticleType::Gamma, 1.0},
        {100.0, CascadeParticleType::Gamma, 0.5}
    };
    dc.level_scheme.levels.resize(4);
    dc.level_scheme.levels[3].feeding = 1.0;
    dc.level_scheme.levels[3].out.push_back(
        {2, 0, 40.0, 1.0, 1.0, 0, 0, 0, 0});
    dc.level_scheme.levels[2].out.push_back(
        {1, 1, 60.0, 1.0, 1.0, 0, 0, 0, 0});
    dc.level_scheme.levels[1].out.push_back(
        {0, 2, 100.0, 1.0, 0.5, 0.5, 0, 0, 0});
    dc.level_scheme.entry_probability = 1.0;
    dc.level_scheme.valid = true;
    dc.weak_outcome_law.outcomes = {
        {WeakOutcomeKind::Other, 3, 1.0, 1.0, 0, 0, 0}
    };
    dc.weak_outcome_law.raw_sum = 0.0;

    const auto channels = EfficiencyCalculator::cascade_sum_pair_channels(
        {dc}, 0, 2, 100.0, 0.1);
    BOOST_REQUIRE_EQUAL(channels.size(), 1u);
    const auto& ch = channels[0];
    BOOST_CHECK_CLOSE(ch.weight, 1.0, 1e-9);
    BOOST_REQUIRE_EQUAL(ch.partner_prob.size(), 3u);
    BOOST_CHECK_SMALL(ch.partner_prob[2], 1e-12);
    BOOST_REQUIRE_EQUAL(ch.weak_outcome_partner_prob.size(), 1u);
    BOOST_CHECK_SMALL(ch.weak_outcome_partner_prob[0][2], 1e-12);
    bool found_primary_ic = false;
    for (const auto& v : ch.vacancies) {
        if (std::abs(v.trans_keV - 100.0) > 1e-9 || v.is_L) continue;
        found_primary_ic = true;
        BOOST_CHECK_CLOSE(v.prob, 1.0, 1e-9);
        BOOST_REQUIRE_EQUAL(v.weak_outcome_prob.size(), 1u);
        BOOST_CHECK_CLOSE(v.weak_outcome_prob[0], 1.0, 1e-9);
    }
    BOOST_CHECK(found_primary_ic);
}

// An unknown fed level means unknown topology, not a missing weak decay.  The
// EC vacancy exists even when the residual outcome lands directly in a terminal
// state and never enters the level graph; only transition-containing joints use
// the residual graph weights.
BOOST_AUTO_TEST_CASE(t14_unknown_terminal_ec_vacancy_is_not_entry_gated) {
    DecayCascade dc;
    dc.level_scheme.levels.resize(2);
    dc.level_scheme.levels[1].feeding = 0.2;
    dc.level_scheme.levels[1].out.push_back(
        {0, 0, 100.0, 1.0, 1.0, 0, 0, 0, 0});
    dc.level_scheme.entry_probability = 0.2;
    dc.level_scheme.valid = true;
    dc.members = {{100.0, CascadeParticleType::Gamma, 0.2}};
    dc.weak_outcome_law.outcomes = {
        {WeakOutcomeKind::ElectronCapture, -1, 1.0, 1.0, 0.8, 0.1, 0.0}
    };
    const LevelDag dag(dc.level_scheme);
    const WeakOutcomeKind ec = WeakOutcomeKind::ElectronCapture;
    BOOST_CHECK_CLOSE(detail::weak_subset_probability(dc, dag, {}, &ec, 1),
                      0.8, 1e-9);
    BOOST_CHECK_CLOSE(detail::weak_subset_probability(dc, dag, {}, &ec, 2),
                      0.1, 1e-9);
    BOOST_CHECK_CLOSE(detail::weak_subset_probability(dc, dag, {0}, &ec),
                      0.2, 1e-9);
}

// A sum-fed pair followed by one of two mutually exclusive, certainly
// depositing branches has zero survival. Multiplying the two conditional
// marginals would invent (1-0.5)^2 = 0.25; the constrained-path DP must sum the
// alternatives and return exactly zero.
BOOST_AUTO_TEST_CASE(t15_pair_survival_preserves_path_exclusivity) {
    DecayCascade primary;
    primary.branch_weight = 1.0;
    primary.members = {{100.0, CascadeParticleType::Gamma, 1.0}};

    DecayCascade feeder;
    feeder.branch_weight = 1.0;
    feeder.members = {
        {40.0, CascadeParticleType::Gamma, 1.0},
        {60.0, CascadeParticleType::Gamma, 1.0},
        {10.0, CascadeParticleType::Gamma, 0.5},
        {20.0, CascadeParticleType::Gamma, 0.5}
    };
    feeder.level_scheme.levels.resize(4);
    feeder.level_scheme.levels[3].feeding = 1.0;
    feeder.level_scheme.levels[3].out.push_back(
        {2, 0, 40.0, 1.0, 1.0, 0, 0, 0, 0});
    feeder.level_scheme.levels[2].out.push_back(
        {1, 1, 60.0, 1.0, 1.0, 0, 0, 0, 0});
    feeder.level_scheme.levels[1].out.push_back(
        {0, 2, 10.0, 0.5, 1.0, 0, 0, 0, 0});
    feeder.level_scheme.levels[1].out.push_back(
        {0, 3, 20.0, 0.5, 1.0, 0, 0, 0, 0});
    feeder.level_scheme.entry_probability = 1.0;
    feeder.level_scheme.valid = true;

    MockProvider eff;
    eff.set(100.0, 0.1, 0.1);
    eff.set(40.0, 0.2, 0.2);
    eff.set(60.0, 0.2, 0.2);
    eff.set(10.0, 1.0, 1.0);
    eff.set(20.0, 1.0, 1.0);
    const auto r = compute_cascade_analytic(
        {primary, feeder}, {{100.0, 0.1}}, eff);
    BOOST_REQUIRE(r[0].found);
    BOOST_CHECK_SMALL(r[0].c_in, 1e-12);
}
