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

// Templating-refactor regression tests for compute_cascade_analytic.
//
// 1) stored_values_unchanged: the expected numbers below were captured from the
//    PRE-templating implementation (AnalyticCascade.cpp at commit dddbdba, plain
//    double throughout) with a temporary capture tool, printed at %.17g. The
//    templated T=double path must reproduce them EXACTLY (same operations, same
//    order) — proving the in-place templating changed no answers.
//    The synthetic schemes exercise every arithmetic path: level-scheme survival
//    DP, EC-K/L and IC-K/L x-ray channels, annihilation-511 pair handling,
//    pairwise fallback (no level scheme), angular-correlated summing-in pairs,
//    x-ray summing-in pairs, triples, multi-line windows, and multi-branch
//    averaging.
//
// 2) dual_number_derivatives: instantiate the template with a minimal forward-
//    mode dual number (value + d/dk for a provider whose efficiencies all scale
//    with k). Value parts must equal the double run exactly; derivative parts
//    must match central finite differences — the smoke test that an AD scalar
//    (e.g. ceres::Jet) flows correctly through the summing factors.

#define BOOST_TEST_MODULE AnalyticCascadeRegression
#include <boost/test/included/unit_test.hpp>

#include "cascade/AnalyticCascade.h"
#include "cascade/CascadeTypes.h"

#include <cmath>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

// ---------------------------------------------------------------------------
// Minimal forward-mode dual number: value `a`, derivative `b` (d/dk). Supports
// exactly the operations AnalyticCascade requires of a generic scalar T.
struct Dual {
    double a = 0.0;
    double b = 0.0;
    Dual() = default;
    Dual(double v) : a(v), b(0.0) {}
    Dual(double v, double d) : a(v), b(d) {}
    Dual& operator+=(const Dual& o) { a += o.a; b += o.b; return *this; }
    Dual& operator*=(const Dual& o) { b = a * o.b + b * o.a; a *= o.a; return *this; }
};
inline Dual operator+(const Dual& x, const Dual& y) { return {x.a + y.a, x.b + y.b}; }
inline Dual operator-(const Dual& x, const Dual& y) { return {x.a - y.a, x.b - y.b}; }
inline Dual operator*(const Dual& x, const Dual& y) {
    return {x.a * y.a, x.a * y.b + x.b * y.a};
}
inline Dual operator/(const Dual& x, const Dual& y) {
    return {x.a / y.a, (x.b * y.a - x.a * y.b) / (y.a * y.a)};
}
inline bool operator>(const Dual& x, double y) { return x.a > y; }
inline bool operator<(const Dual& x, double y) { return x.a < y; }
inline bool operator>=(const Dual& x, double y) { return x.a >= y; }
inline bool operator<=(const Dual& x, double y) { return x.a <= y; }
inline Dual sqrt(const Dual& x) {
    const double r = std::sqrt(x.a);
    return {r, x.b / (2.0 * r)};
}

// ---------------------------------------------------------------------------
// Smooth functional providers: defined for every energy (has() always true) so
// x-ray lines from the repo's fluorescence tables get consistent values. The
// double and Dual providers use identical value arithmetic; the Dual one seeds
// d/dk = the unscaled efficiency (eps = k * base  =>  d eps / dk = base).
double base_fep(double E) { return 0.03 * std::exp(-E / 800.0); }
double base_tot(double E) { return 0.20 * std::exp(-E / 1500.0); }

class FuncProvider : public EfficiencyProvider {
public:
    explicit FuncProvider(double scale = 1.0) : k_(scale) {}
    double fep(double E) const override { return k_ * base_fep(E); }
    double total(double E) const override { return k_ * base_tot(E); }
    double fep_unc(double E) const override { return 0.05 * fep(E); }
    bool has(double) const override { return true; }
private:
    double k_;
};

class DualProvider : public EfficiencyProviderT<Dual> {
public:
    explicit DualProvider(double scale = 1.0) : k_(scale) {}
    Dual fep(double E) const override { return {k_ * base_fep(E), base_fep(E)}; }
    Dual total(double E) const override { return {k_ * base_tot(E), base_tot(E)}; }
    Dual fep_unc(double E) const override {
        return {0.05 * k_ * base_fep(E), 0.05 * base_fep(E)};
    }
    bool has(double) const override { return true; }
private:
    double k_;
};

// ---------------------------------------------------------------------------
// Synthetic decay schemes (identical to the pre-templating capture tool).

// Co-60-like sequential scheme with an angular-correlated 1173<->1332 link.
DecayCascade make_co60_like() {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 28;
    dc.members.resize(2);
    dc.members[0].energy_keV = 1173.2;
    dc.members[1].energy_keV = 1332.5;
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    dc.members[0].coincident.push_back({1, 0.9998, 0.1020, 0.0091, true});
    LevelScheme& ls = dc.level_scheme;
    ls.daughter_Z = 28;
    ls.levels.resize(3);
    ls.levels[2].feeding = 1.0;
    ls.levels[2].out.push_back({1, 0, 1173.2, 1.0, 0.9998, 0, 0, 0, 0});
    ls.levels[1].out.push_back({0, 1, 1332.5, 1.0, 0.9998, 0, 0, 0, 0});
    ls.valid = true;
    return dc;
}

// Co-57-like EC + IC scheme: exercises EC-K/L x-rays, IC-K/L x-rays, exclusivity.
DecayCascade make_ec_ic() {
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
    ls.levels[2].feed_ecK = 0.86;
    ls.levels[2].feed_ecL = 0.09;
    ls.levels[2].out.push_back({0, 0, 136.47, 0.115, 0.86, 0.10, 0.02, 0.01, 0.005});
    ls.levels[2].out.push_back({1, 1, 122.06, 0.885, 0.86, 0.10, 0.02, 0.01, 0.005});
    ls.levels[1].out.push_back({0, 2, 14.41, 1.000, 0.10, 0.60, 0.20, 0.05, 0.02});
    ls.valid = true;
    return dc;
}

// Na-22-like beta+ scheme: gamma 1274.5 + back-to-back 511 pair member.
DecayCascade make_beta_plus() {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 10;
    dc.members.resize(2);
    dc.members[0].energy_keV = 1274.5;
    dc.members[0].type = CascadeParticleType::Gamma;
    dc.members[1].energy_keV = 510.998950;
    dc.members[1].type = CascadeParticleType::Annih511;
    dc.members[1].intensity = 0.899;  // positron-event probability; pair has 2 photons
    LevelScheme& ls = dc.level_scheme;
    ls.daughter_Z = 10;
    ls.levels.resize(2);
    ls.levels[1].feeding = 1.0;
    ls.levels[1].out.push_back({0, 0, 1274.5, 1.0, 0.9994, 0, 0, 0, 0});
    ls.valid = true;
    return dc;
}

// Sequential 3-gamma scheme with a 200 keV crossover, so the triple summing-in
// 60 + 40 + 100 = 200 lands on an emitted line (exercises the triples path).
DecayCascade make_seq3() {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 0;
    dc.members.resize(4);
    dc.members[0].energy_keV = 60.0;
    dc.members[1].energy_keV = 40.0;
    dc.members[2].energy_keV = 100.0;
    dc.members[3].energy_keV = 200.0;  // crossover level3 -> ground
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    LevelScheme& ls = dc.level_scheme;
    ls.levels.resize(4);
    ls.levels[3].feeding = 1.0;
    ls.levels[3].out.push_back({2, 0, 60.0, 0.95, 0.9, 0, 0, 0, 0});
    ls.levels[3].out.push_back({0, 3, 200.0, 0.05, 0.9, 0, 0, 0, 0});
    ls.levels[2].out.push_back({1, 1, 40.0, 1.0, 0.8, 0, 0, 0, 0});
    ls.levels[1].out.push_back({0, 2, 100.0, 1.0, 0.95, 0, 0, 0, 0});
    ls.valid = true;
    return dc;
}

// No-level-scheme pairwise fallback branch (dag invalid).
DecayCascade make_pairwise() {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 28;
    dc.members.resize(2);
    dc.members[0].energy_keV = 1173.2;
    dc.members[1].energy_keV = 1332.5;
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    dc.members[0].coincident.push_back({1, 0.9998});
    dc.members[1].coincident.push_back({0, 0.9985});
    return dc;
}

// Two lines inside one window (multi-line component test): 79.6 + 81.0.
DecayCascade make_two_line_window() {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.daughter_Z = 56;
    dc.members.resize(4);
    dc.members[0].energy_keV = 79.6;
    dc.members[1].energy_keV = 81.0;
    dc.members[2].energy_keV = 276.4;
    dc.members[3].energy_keV = 40.0;  // p_gamma = 0 transition (pure IC)
    for (auto& m : dc.members) m.type = CascadeParticleType::Gamma;
    LevelScheme& ls = dc.level_scheme;
    ls.daughter_Z = 56;
    ls.levels.resize(4);
    ls.levels[3].feeding = 1.0;
    ls.levels[3].out.push_back({2, 2, 276.4, 0.6, 0.9, 0.05, 0.01, 0.005, 0.002});
    ls.levels[3].out.push_back({1, 0, 79.6, 0.4, 0.7, 0.15, 0.05, 0.02, 0.01});
    ls.levels[2].out.push_back({1, 1, 81.0, 1.0, 0.75, 0.12, 0.04, 0.02, 0.01});
    ls.levels[1].out.push_back({0, 3, 40.0, 1.0, 0.0, 0.3, 0.1, 0.05, 0.02});
    ls.valid = true;
    return dc;
}

struct Case {
    std::string name;
    std::vector<DecayCascade> cascades;
    std::vector<PeakWindow> peaks;
};

std::vector<Case> make_cases() {
    std::vector<Case> cases;
    cases.push_back({"co60", {make_co60_like()},
                     {{1173.2, 1.5}, {1332.5, 1.5}, {2505.7, 1.5}}});
    cases.push_back({"ec_ic", {make_ec_ic()},
                     {{136.47, 1.5}, {122.06, 1.5}, {14.41, 1.5}}});
    cases.push_back({"beta_plus", {make_beta_plus()},
                     {{510.998950, 1.5}, {1274.5, 1.5}, {1785.5, 1.5}, {1021.9979, 1.5}}});
    cases.push_back({"seq3", {make_seq3()}, {{60.0, 1.5}, {100.0, 1.5}, {200.0, 1.5}}});
    cases.push_back({"pairwise", {make_pairwise()}, {{1173.2, 1.5}, {1332.5, 1.5}}});
    cases.push_back({"two_line", {make_two_line_window()}, {{80.0, 1.6}, {276.4, 1.5}}});
    cases.push_back({"multi_branch",
                     {make_co60_like(), make_pairwise(), make_beta_plus()},
                     {{1173.2, 1.5}, {511.0, 1.5}, {2505.7, 1.5}}});
    return cases;
}

struct Expected {
    const char* name;
    double energy;
    int found;
    double c_out, c_in, c_net;
    double eff_no, eff_with, eff_with_unc, sf_unc;
    std::size_t n_unmatched;
};

// Captured from the pre-templating implementation (see file header). %.17g.
// Two c_in values (ec_ic 14.41, seq3 200) were re-baselined by 1 ULP when the
// summing-IN pair loop was changed to a sorted O(n log n) energy-window search
// (the summed set is identical; only the float accumulation ORDER differs).
// The 511-primary rows were intentionally re-baselined for the annihilation
// geometry correction: once one back-to-back photon is conditioned into a
// convex single-sided detector, its partner points away and cannot sum it out.
// The beta-plus synthetic now stores the positron-event probability (0.899),
// not the two-photon line yield (1.798), matching CascadeMember's contract.
// The EC/IC and two-line x-ray rows were re-baselined in Aug 2026 when the
// relaxation fixtures moved from the Geant4 extraction to direct EPICS2023
// EADL. Peak-efficiency shifts are below 0.005%; every recorded deterministic
// field remains within the 1% migration gate (the largest is a small c_in term).
const Expected kExpected[] = {
    { "co60", 1173.200000, 1, 0.91774830369415894, 0, 0.91774830369415894, 0.0069219490296876587, 0.0063526069802532781, 0.00031763034901266394, 0.045887415184707951, 0 },
    { "co60", 1332.500000, 1, 0.90853247299122142, 0, 0.90853247299122142, 0.0056721735229304012, 0.0051533538380232856, 0.00025766769190116434, 0.045426623649561078, 0 },
    { "co60", 2505.700000, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
    { "ec_ic", 136.470000, 1, 0.94241371779491057, 0.02175746452822163, 0.96417118232313215, 0.025295079272586504, 0.024388786489207084, 0.0011919214849597615, 0.04712068588974553, 0 },
    { "ec_ic", 122.060000, 1, 0.88588278349663863, 0, 0.88588278349663863, 0.025754835133619532, 0.022815765036667893, 0.0011407882518333949, 0.044294139174831937, 0 },
    { "ec_ic", 14.410000, 1, 0.78667590393457154, 0.0040110572029001907, 0.79068696113747183, 0.029464462662663813, 0.023297166444290151, 0.0011589491399548746, 0.039333795196728583, 0 },
    { "beta_plus", 510.998950, 1, 0.91453982142315648, 0, 0.91453982142315648, 0.015838579338881949, 0.014485011520177593, 0.00072425057600887966, 0.045726991071157824, 0 },
    { "beta_plus", 1274.500000, 1, 0.74421779109340591, 0, 0.74421779109340591, 0.0060986801667498184, 0.0045387462822837144, 0.00022693731411418575, 0.0372108895546703, 0 },
    { "beta_plus", 1785.500000, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
    { "beta_plus", 1021.997900, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
    { "seq3", 60.000000, 1, 0.69415500346775605, 0, 0.69415500346775605, 0.027832304589856586, 0.019319933489087542, 0.000965996674454377, 0.0347077501733878, 0 },
    { "seq3", 100.000000, 1, 0.69821078222498911, 0.018695451929926568, 0.71690623415491561, 0.026474907077537864, 0.018980025932558993, 0.00092425327899708074, 0.034910539111249458, 0 },
    { "seq3", 200.000000, 1, 1, 0.012996000000000001, 1.012996, 0.023364023492142144, 0.023667662341446023, 0.0011682011746071072, 0.050000000000000003, 0 },
    { "pairwise", 1173.200000, 1, 0.91774830369415894, 0, 0.91774830369415894, 0.0069219490296876587, 0.0063526069802532781, 0.00031763034901266394, 0.045887415184707951, 0 },
    { "pairwise", 1332.500000, 1, 0.90865140456264715, 0, 0.90865140456264715, 0.0056721735229304012, 0.0051540284385337672, 0.0002577014219266884, 0.045432570228132366, 0 },
    { "two_line", 80.000000, 1, 0.85445141047467799, 0, 0.85445141047467799, 0.027129426205730972, 0.023180776486855521, 0.001159709412136637, 0.042747288620931224, 0 },
    { "two_line", 276.400000, 1, 0.79226887498114407, 0, 0.79226887498114407, 0.021235989954115724, 0.016824613870058142, 0.00084123069350290722, 0.03961344374905721, 0 },
    { "multi_branch", 1173.200000, 1, 0.91774830369415883, 0, 0.91774830369415883, 0.0069219490296876587, 0.0063526069802532772, 0.00031763034901266388, 0.045887415184707944, 0 },
    { "multi_branch", 511.000000, 1, 0.91453982142315648, 0, 0.91453982142315648, 0.015838579338881949, 0.014485011520177593, 0.00072424962543062238, 0.045726931054521421, 0 },
    { "multi_branch", 2505.700000, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
};

}  // namespace

// ===========================================================================

BOOST_AUTO_TEST_CASE(stored_values_unchanged) {
    const FuncProvider prov;
    const AnalyticCascadeOptions opts;  // defaults: triples + angular on

    std::size_t idx = 0;
    for (const Case& cs : make_cases()) {
        const std::vector<AnalyticPeakResult> rr =
            compute_cascade_analytic(cs.cascades, cs.peaks, prov, opts);
        BOOST_REQUIRE_EQUAL(rr.size(), cs.peaks.size());
        for (const AnalyticPeakResult& r : rr) {
            BOOST_REQUIRE_LT(idx, sizeof(kExpected) / sizeof(kExpected[0]));
            const Expected& e = kExpected[idx];
            BOOST_TEST_CONTEXT("case " << e.name << " @ " << e.energy << " keV") {
                BOOST_CHECK_EQUAL(cs.name, e.name);
                BOOST_CHECK_EQUAL(r.energy_keV, e.energy);
                BOOST_CHECK_EQUAL((int)r.found, e.found);
                BOOST_CHECK_EQUAL(r.c_out, e.c_out);
                BOOST_CHECK_EQUAL(r.c_in, e.c_in);
                BOOST_CHECK_EQUAL(r.c_net, e.c_net);
                BOOST_CHECK_EQUAL(r.eff_no_summing, e.eff_no);
                BOOST_CHECK_EQUAL(r.eff_with_summing, e.eff_with);
                BOOST_CHECK_EQUAL(r.eff_with_summing_unc, e.eff_with_unc);
                BOOST_CHECK_EQUAL(r.summing_factor_unc, e.sf_unc);
                BOOST_CHECK_EQUAL(r.unmatched_energies.size(), e.n_unmatched);
            }
            ++idx;
        }
    }
    BOOST_CHECK_EQUAL(idx, sizeof(kExpected) / sizeof(kExpected[0]));
}

// Dual-number instantiation: value parts identical to the double run; d/dk
// matches central finite differences over the provider scale k.
BOOST_AUTO_TEST_CASE(dual_number_derivatives) {
    const AnalyticCascadeOptions opts;
    const double h = 1e-6;

    for (const Case& cs : make_cases()) {
        const FuncProvider p0(1.0), pp(1.0 + h), pm(1.0 - h);
        const DualProvider pd(1.0);

        const std::vector<AnalyticPeakResult> r0 =
            compute_cascade_analytic(cs.cascades, cs.peaks, p0, opts);
        const std::vector<AnalyticPeakResult> rp =
            compute_cascade_analytic(cs.cascades, cs.peaks, pp, opts);
        const std::vector<AnalyticPeakResult> rm =
            compute_cascade_analytic(cs.cascades, cs.peaks, pm, opts);
        const std::vector<AnalyticPeakResultT<Dual>> rd =
            compute_cascade_analytic(cs.cascades, cs.peaks, pd, opts);

        BOOST_REQUIRE_EQUAL(rd.size(), r0.size());
        for (std::size_t i = 0; i < r0.size(); ++i) {
            BOOST_TEST_CONTEXT("case " << cs.name << " @ " << r0[i].energy_keV
                                       << " keV") {
                // Value parts: the double computation to within FP-contraction
                // noise (the compiler may fuse/reorder the double instantiation
                // differently from the Dual value lane; ~1 ulp).
                auto check_val = [](double dual_a, double dbl) {
                    if (dbl == 0.0)
                        BOOST_CHECK_EQUAL(dual_a, 0.0);
                    else
                        BOOST_CHECK_CLOSE(dual_a, dbl, 1e-10 /*percent*/);
                };
                check_val(rd[i].c_net.a, r0[i].c_net);
                check_val(rd[i].c_out.a, r0[i].c_out);
                check_val(rd[i].eff_with_summing.a, r0[i].eff_with_summing);

                // Derivatives vs central finite differences.
                const double fd_cnet = (rp[i].c_net - rm[i].c_net) / (2.0 * h);
                const double fd_eff =
                    (rp[i].eff_with_summing - rm[i].eff_with_summing) / (2.0 * h);
                const double tol_cnet =
                    1e-6 * std::max(std::abs(fd_cnet), 1e-6);
                const double tol_eff = 1e-6 * std::max(std::abs(fd_eff), 1e-9);
                BOOST_CHECK_SMALL(rd[i].c_net.b - fd_cnet, tol_cnet);
                BOOST_CHECK_SMALL(rd[i].eff_with_summing.b - fd_eff, tol_eff);
            }
        }
    }
}
