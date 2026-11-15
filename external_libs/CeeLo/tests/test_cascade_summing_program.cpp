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

// SummingOutProgram and coincidence_sum_pairs (cascade/CascadeSummingProgram.h) against the engine they
// repackage, compute_cascade_analytic, on real nuclides: every weak-decay kind (alpha, beta-, EC, beta+), aged
// chains, branches with and without an accepted level graph.  No Monte Carlo; both sides are deterministic.

#define BOOST_TEST_MODULE CascadeSummingProgram
#include <boost/test/unit_test.hpp>

#include "cascade/AnalyticCascade.h"
#include "cascade/SandiaDecayCascade.h"
#include "cascade/CascadeSummingProgram.h"
#include "SandiaDecay.h"

#include <map>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

using namespace ceelo;

namespace {

const SandiaDecay::SandiaDecayDataBase& db() {
    static SandiaDecay::SandiaDecayDataBase database(SANDIA_DECAY_XML_PATH);
    return database;
}

const double kYear = 365.25 * 24.0 * 3600.0;

struct Case { const char* nuclide; double age_years; };

// Long-lived chain parents at 20 y, as an isotopics fit ages them; calibration nuclides at 1 y.
const Case kCases[] = {
    {"U235", 20.0}, {"U238", 20.0}, {"U234", 20.0}, {"U232", 20.0}, {"Th232", 20.0}, {"Ra226", 20.0},
    {"Am241", 20.0}, {"Pu239", 20.0}, {"Eu152", 1.0}, {"Ba133", 1.0}, {"Co60", 1.0}, {"Na22", 1.0},
    {"Co57", 0.5}, {"Cs134", 1.0}, {"Y88", 0.3}, {"Eu154", 1.0}
};

std::shared_ptr<const AnalyticCascadeModel> make_model(const Case& c) {
    const SandiaDecay::Nuclide* const nuc = db().nuclide(c.nuclide);
    BOOST_REQUIRE(nuc);
    cascade_adapter::CascadeOptions opts;
    opts.age_seconds = c.age_years * kYear;
    opts.prompt_equilibrium = false;
    return std::make_shared<const AnalyticCascadeModel>(cascade_adapter::build_cascades(nuc, opts));
}

// A smooth, made-up detector: total efficiency up to ~0.3 (a close geometry), full-energy efficiency flat, so
// that a window's lines are weighted by emission alone - what SummingOutProgram does.
double total_eff(const double energy) {
    return 0.30 * std::exp(-energy / 700.0) * (1.0 - std::exp(-energy / 25.0)) + 0.02;
}

struct Provider final : EfficiencyProvider {
    double tot_scale = 1.0;
    double fep_value = 0.05;
    double fep(double) const override { return fep_value; }
    double total(double e) const override { return tot_scale * total_eff(e); }
};

// The distinct gamma (and annihilation) member energies emitted above `min_emission` per parent decay.
std::vector<double> line_energies(const AnalyticCascadeModel& model, const double min_emission) {
    std::map<double, double> emission;
    for (const DecayCascade& dc : model.cascades())
        for (const CascadeMember& m : dc.members)
            emission[m.energy_keV] += dc.branch_weight * m.intensity;
    std::vector<double> energies;
    for (const auto& e : emission)
        if (e.second >= min_emission && e.first >= 20.0) energies.push_back(e.first);
    return energies;
}

}  // namespace


BOOST_AUTO_TEST_CASE(summing_out_matches_engine) {
    for (const Case& c : kCases) {
        const std::shared_ptr<const AnalyticCascadeModel> model = make_model(c);
        const std::vector<double> energies = line_energies(*model, 1.0e-7);
        BOOST_REQUIRE_MESSAGE(!energies.empty(), c.nuclide);

        const double tol = 1.0e-4;  // each window holds one line energy
        std::vector<PeakWindow> windows;
        std::vector<SummingOutProgram::Line> lines;
        for (const double e : energies) {
            PeakWindow w;
            w.energy_keV = e;
            w.tolerance_keV = tol;
            windows.push_back(w);
            lines.push_back({e, tol});
        }

        const SummingOutProgram program(model, lines);
        SummingOutProgram::Workspace<double> ws;
        std::vector<double> tot(program.energies().size()), out(lines.size());

        for (const double scale : {0.0, 0.1, 1.0}) {
            Provider prov;
            prov.tot_scale = scale;
            const std::vector<AnalyticPeakResult> ref = compute_cascade_analytic(*model, windows, prov);

            for (std::size_t k = 0; k < tot.size(); ++k) tot[k] = prov.total(program.energies()[k]);
            program.evaluate(tot.data(), out.data(), ws);

            double worst = 0.0;
            std::size_t nfound = 0, nincomplete = 0;
            for (std::size_t i = 0; i < lines.size(); ++i) {
                if (!ref[i].found) {
                    BOOST_CHECK_EQUAL(program.line_emission(i), 0.0);
                    BOOST_CHECK_EQUAL(out[i], 1.0);
                    continue;
                }
                ++nfound;
                nincomplete += program.line_complete(i) ? 0 : 1;
                BOOST_CHECK_GT(program.line_emission(i), 0.0);
                const double diff = std::abs(out[i] - ref[i].c_out);
                worst = std::max(worst, diff);
                BOOST_CHECK_MESSAGE(diff < 1.0e-12, c.nuclide << " " << energies[i] << " keV, scale " << scale
                                    << ": program " << out[i] << " engine " << ref[i].c_out);
                if (scale == 0.0) BOOST_CHECK_SMALL(out[i] - 1.0, 1.0e-15);
            }
            BOOST_TEST_MESSAGE(c.nuclide << " scale " << scale << ": " << nfound << " lines (" << nincomplete
                               << " without a level graph), worst |diff| " << worst);
            BOOST_CHECK_GT(nfound, 0u);
        }
    }
}


BOOST_AUTO_TEST_CASE(first_order_coefficients_are_the_derivative) {
    const std::shared_ptr<const AnalyticCascadeModel> model = make_model({"Co60", 1.0});
    const SummingOutProgram program(model, {{1173.228, 0.01}, {1332.492, 0.01}});
    // Emission is per decay of the parent at age zero: after a year, 0.877 of Co-60's initial activity is left.
    BOOST_CHECK_CLOSE(program.line_emission(0), 0.9985 * std::exp(-std::log(2.0) / 5.2711), 0.1);
    BOOST_CHECK_CLOSE(program.line_emission(1), 0.9998 * std::exp(-std::log(2.0) / 5.2711), 0.1);
    const std::vector<std::vector<double>> coef = program.first_order_coefficients();
    // Each Co-60 line's partner is the other one, emitted in (nearly) every decay that emits it.
    const std::vector<double>& energies = program.energies();
    for (std::size_t k = 0; k < energies.size(); ++k) {
        if (std::abs(energies[k] - 1332.492) < 0.01) BOOST_CHECK_CLOSE(coef[0][k], 1.0, 0.5);
        if (std::abs(energies[k] - 1173.228) < 0.01) BOOST_CHECK_CLOSE(coef[1][k], 1.0, 0.5);
    }
}


BOOST_AUTO_TEST_CASE(sum_pairs_match_engine_pure_sum_rate) {
    // With unit full-energy efficiency and zero total efficiency, the engine's pure-sum rate of a window is the
    // sum of the weights of the pairs (and, with triples off, nothing else) adding into it.
    AnalyticCascadeOptions opts;
    opts.report_pure_sums = true;
    opts.enumerate_triples = false;
    Provider prov;
    prov.tot_scale = 0.0;
    prov.fep_value = 1.0;

    for (const Case& c : kCases) {
        const std::shared_ptr<const AnalyticCascadeModel> model = make_model(c);
        const double lo = 30.0, hi = 3000.0, half_window = 0.4;
        const std::vector<CoincidenceSumPair> pairs = coincidence_sum_pairs(*model, lo, hi);
        BOOST_REQUIRE_MESSAGE(!pairs.empty(), c.nuclide);

        // Windows at the sums of the strongest pairs.
        std::vector<CoincidenceSumPair> strongest = pairs;
        std::sort(strongest.begin(), strongest.end(),
                  [](const CoincidenceSumPair& a, const CoincidenceSumPair& b) { return a.weight > b.weight; });
        std::vector<PeakWindow> windows;
        for (const CoincidenceSumPair& p : strongest) {
            if (windows.size() >= 60) break;
            PeakWindow w;
            w.energy_keV = p.energy_a + p.energy_b;
            w.tolerance_keV = half_window;
            // The whole window must lie inside the range the pairs were listed for.
            if (w.energy_keV - half_window < lo || w.energy_keV + half_window > hi) continue;
            windows.push_back(w);
        }
        const std::vector<AnalyticPeakResult> ref = compute_cascade_analytic(*model, windows, prov, opts);

        std::size_t nchecked = 0;
        for (std::size_t i = 0; i < windows.size(); ++i) {
            if (ref[i].found) continue;  // holds an emitted line: not a pure-sum window
            double expected = 0.0;
            for (const CoincidenceSumPair& p : pairs) {
                const double sum = p.energy_a + p.energy_b;
                if (std::abs(sum - windows[i].energy_keV) <= windows[i].tolerance_keV) expected += p.weight;
            }
            ++nchecked;
            BOOST_CHECK_MESSAGE(std::abs(expected - ref[i].pure_sum_rate) <= 1.0e-12 * std::max(expected, 1.0e-30)
                                + 1.0e-300,
                                c.nuclide << " sum " << windows[i].energy_keV << " keV: pairs " << expected
                                << " engine " << ref[i].pure_sum_rate);
        }
        BOOST_TEST_MESSAGE(c.nuclide << ": " << pairs.size() << " pairs, " << nchecked << " pure-sum windows checked");
        BOOST_CHECK_GT(nchecked, 0u);
    }
}


BOOST_AUTO_TEST_CASE(pair_emissions_are_consistent) {
    // weight / emission is a conditional probability times W(0): bounded, and symmetric in which member is given.
    const std::shared_ptr<const AnalyticCascadeModel> model = make_model({"Ba133", 1.0});
    const std::vector<CoincidenceSumPair> pairs = coincidence_sum_pairs(*model, 100.0, 500.0, false);
    BOOST_REQUIRE(!pairs.empty());
    for (const CoincidenceSumPair& p : pairs) {
        BOOST_CHECK_LE(p.energy_a, p.energy_b);
        BOOST_CHECK_GT(p.emission_a, 0.0);
        BOOST_CHECK_GT(p.emission_b, 0.0);
        BOOST_CHECK_EQUAL(p.weight, p.joint);  // no angular factor asked for
        BOOST_CHECK_LE(p.joint, p.emission_a * (1.0 + 1.0e-9));
        BOOST_CHECK_LE(p.joint, p.emission_b * (1.0 + 1.0e-9));
    }
}


BOOST_AUTO_TEST_CASE(pair_survival_matches_engine) {
    // With real efficiencies the engine's pure-sum rate of a window is the sum over its pairs of
    // weight * eps_fep(a) * eps_fep(b) * (the pair's coincident survival); SummingOutProgram gives that survival
    // for pairs of level-scheme transitions.
    struct RealProvider final : EfficiencyProvider {
        double fep(double e) const override { return 0.25 * std::exp(-e / 500.0) + 0.005; }
        double total(double e) const override { return total_eff(e); }
    };
    const RealProvider prov;
    AnalyticCascadeOptions opts;
    opts.report_pure_sums = true;
    opts.enumerate_triples = false;

    for (const Case& c : kCases) {
        const std::shared_ptr<const AnalyticCascadeModel> model = make_model(c);
        const double lo = 60.0, hi = 3000.0, half_window = 0.4;
        const std::vector<CoincidenceSumPair> pairs = coincidence_sum_pairs(*model, lo, hi);

        std::vector<std::size_t> order(pairs.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        auto rate0 = [&](const CoincidenceSumPair& p) { return p.weight * prov.fep(p.energy_a) * prov.fep(p.energy_b); };
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return rate0(pairs[a]) > rate0(pairs[b]); });

        // Windows at the strongest pairs' sums, kept when every pair in them is transition x transition.
        std::vector<PeakWindow> windows;
        std::vector<std::vector<std::size_t>> window_pairs;
        std::vector<SummingOutProgram::Pair> requests;
        std::vector<int> request_of(pairs.size(), -1);
        for (const std::size_t i : order) {
            if (windows.size() >= 40) break;
            const double sum = pairs[i].energy_a + pairs[i].energy_b;
            if (sum - half_window < lo || sum + half_window > hi) continue;
            std::vector<std::size_t> in_window;
            bool all_transitions = true;
            for (std::size_t j = 0; j < pairs.size(); ++j) {
                if (std::abs(pairs[j].energy_a + pairs[j].energy_b - sum) > half_window) continue;
                in_window.push_back(j);
                all_transitions = all_transitions && (pairs[j].transition_a >= 0) && (pairs[j].transition_b >= 0);
            }
            if (!all_transitions) continue;
            for (const std::size_t j : in_window) {
                if (request_of[j] >= 0) continue;
                request_of[j] = static_cast<int>(requests.size());
                requests.push_back({pairs[j].branch, pairs[j].transition_a, pairs[j].transition_b});
            }
            PeakWindow w;
            w.energy_keV = sum;
            w.tolerance_keV = half_window;
            windows.push_back(w);
            window_pairs.push_back(std::move(in_window));
        }
        if (windows.empty()) {
            BOOST_TEST_MESSAGE(c.nuclide << ": no window of only transition pairs");
            continue;
        }

        const SummingOutProgram program(model, {}, requests);
        BOOST_REQUIRE_EQUAL(program.num_pairs(), requests.size());
        SummingOutProgram::Workspace<double> ws;
        std::vector<double> tot(program.energies().size()), survival(requests.size(), -1.0);
        for (std::size_t k = 0; k < tot.size(); ++k) tot[k] = prov.total(program.energies()[k]);
        program.evaluate<double>(tot.data(), nullptr, ws, survival.data());

        const std::vector<AnalyticPeakResult> ref = compute_cascade_analytic(*model, windows, prov, opts);
        double worst = 0.0, lowest = 1.0;
        std::size_t nchecked = 0;
        for (std::size_t w = 0; w < windows.size(); ++w) {
            if (ref[w].found) continue;
            double expected = 0.0;
            for (const std::size_t j : window_pairs[w]) {
                const double surv = survival[static_cast<std::size_t>(request_of[j])];
                BOOST_CHECK(surv > 0.0 && surv <= 1.0 + 1.0e-12);
                lowest = std::min(lowest, surv);
                expected += rate0(pairs[j]) * surv;
            }
            ++nchecked;
            const double rel = std::abs(expected - ref[w].pure_sum_rate) / std::max(ref[w].pure_sum_rate, 1.0e-300);
            worst = std::max(worst, rel);
            BOOST_CHECK_MESSAGE(rel < 1.0e-10, c.nuclide << " sum " << windows[w].energy_keV << " keV: program "
                                << expected << " engine " << ref[w].pure_sum_rate);
        }
        BOOST_TEST_MESSAGE(c.nuclide << ": " << nchecked << " windows, lowest pair survival " << lowest
                           << ", worst relative difference " << worst);
    }
}
