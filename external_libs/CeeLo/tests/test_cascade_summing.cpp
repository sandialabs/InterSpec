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

// Cascade true-coincidence summing regression gate: runs compute_cascade at LOW
// statistics and checks the engine's per-decay peak / sum-peak observables against
// GEANT4 references. Both the references AND this gate are driven from ONE CSV,
// tests/data/geant4_reference/cascade_summing_multi.csv, which is also consumed by
// profiling/compare_validation.py -- the numbers are no longer hardcoded here.
//
// The GEANT4 references were generated at HIGH statistics (8M decays for
// Co-60/Ba-133, 4M for the rest) with the GEANT4 11.4.0 harness
// (ceelo_g4val ... --histogram --correlated-gamma) in the SAME geometry this test
// builds: a 3"x3" NaI cylinder with the source distributed in a small solid Al
// cylinder (r = halfz = source_cm, from the CSV) 2 cm behind the front face.
// Critically, the G4 GPS sampling volume == the engine's source region -- a
// mismatch there silently biases the L x-ray self-absorption.
//
// Each observable is checked as an engine/G4 RATIO against a band (from the CSV)
// chosen to (a) comfortably pass the validated agreement plus low-stats scatter,
// and (b) catch a real regression (e.g. the level-path breaking would send Am-241
// 59.5+L back to ~3x). Every full row is gated via the FullRealization estimator;
// clean summing-OUT gamma lines additionally carry a "cond" row gated via the
// Conditional estimator (A_FR * k_cond/k_full). z-scores (engine vs G4, statistical
// only) are printed via BOOST_TEST_MESSAGE for diagnosis.

#define BOOST_TEST_MODULE CascadeSummingTests
#include <boost/test/unit_test.hpp>

// The alcyl compute_cascade run, the spectrum-window integrators, and the
// reference-CSV reader live in cascade_ref_common.h -- SHARED verbatim with the
// compare_validation.py producer (examples/cascade_observables.cpp) so the geometry
// and observable definitions have exactly one implementation.
#include "cascade_ref_common.h"

#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"
#include "cascade/AnalyticCascade.h"
#include "materials/Material.h"
#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

using namespace ceelo;
using namespace ceelo::cascade_adapter;
using namespace ceelo::cascade_ref;  // run_fullrealization / run_conditional / area / Obs / db

namespace {

// The reference CSV, parsed once and grouped by nuclide (single source of truth,
// shared with the producer via CASCADE_SUMMING_REF).
const std::map<std::string, NuclideRefs>& ref_by_nuclide() {
    static const std::map<std::string, NuclideRefs> g =
        group_reference(load_reference(CASCADE_SUMMING_REF));
    return g;
}

// Run the alcyl geometry for `nuc` and gate every reference observable against G4:
//   full rows -> FullRealization summed-spectrum area (ratio band + z-score)
//   cond rows -> Conditional per-decay area = A_FR(window) * k_cond/k_full, vs the
//                SAME G4 area (a slightly looser band absorbs the Conditional-vs-FR gap)
void check_nuclide_from_csv(const std::string& nuc, uint64_t nev_full,
                            uint64_t nev_cond) {
    const auto it = ref_by_nuclide().find(nuc);
    BOOST_REQUIRE_MESSAGE(it != ref_by_nuclide().end(),
                          nuc << ": not present in " << CASCADE_SUMMING_REF);
    const NuclideRefs& nr = it->second;

    std::vector<double> cond_peaks;
    for (const Obs& o : nr.cond) cond_peaks.push_back(o.peak_keV);

    const FullResult fr =
        run_fullrealization(nuc, nr.emax, nev_full, nr.source_cm, cond_peaks);

    // FullRealization rows: summed-spectrum window area vs the baked G4 area.
    for (const Obs& o : nr.full) {
        const double e = area(fr.spectrum, o.lo, o.hi);
        const double eu = area_unc(fr.spectrum_unc, o.lo, o.hi);
        const double ratio = (o.g4 > 0.0) ? e / o.g4 : 0.0;
        const double sig = std::sqrt(eu * eu + o.g4_sig * o.g4_sig);
        const double z = (sig > 0.0) ? (e - o.g4) / sig : 0.0;
        BOOST_TEST_MESSAGE(nuc << " full " << o.name << ": eng=" << e << " G4=" << o.g4
                           << " ratio=" << ratio << " (band " << o.rlo << "-" << o.rhi
                           << ") z=" << z);
        BOOST_CHECK_MESSAGE(ratio > o.rlo && ratio < o.rhi,
                            nuc << " full " << o.name << ": engine/G4 ratio " << ratio
                            << " outside [" << o.rlo << ", " << o.rhi << "] "
                            << "(eng=" << e << ", G4=" << o.g4 << ")");
    }

    // Conditional rows: A_FR(window) * k_cond/k_full vs the same G4 area.
    if (!nr.cond.empty()) {
        std::vector<double> kc, kcu;
        run_conditional(nuc, nev_cond, nr.source_cm, cond_peaks, kc, kcu);
        for (std::size_t i = 0; i < nr.cond.size(); ++i) {
            const Obs& o = nr.cond[i];
            const double afr = area(fr.spectrum, o.lo, o.hi);
            const double afru = area_unc(fr.spectrum_unc, o.lo, o.hi);
            const double kf = fr.k[i], kfu = fr.k_unc[i];
            const double kcc = kc[i], kccu = kcu[i];
            double e = 0.0, eu = 0.0;
            if (kf > 0.0 && kcc > 0.0 && afr > 0.0) {
                e = afr * kcc / kf;
                const double rel = (afru / afr) * (afru / afr)
                                 + (kccu / kcc) * (kccu / kcc)
                                 + (kfu / kf) * (kfu / kf);
                eu = e * std::sqrt(rel);
            }
            const double ratio = (o.g4 > 0.0) ? e / o.g4 : 0.0;
            const double sig = std::sqrt(eu * eu + o.g4_sig * o.g4_sig);
            const double z = (sig > 0.0) ? (e - o.g4) / sig : 0.0;
            BOOST_TEST_MESSAGE(nuc << " cond " << o.name << ": eng=" << e
                               << " (A_FR=" << afr << " k_cond=" << kcc
                               << " k_full=" << kf << ") G4=" << o.g4 << " ratio="
                               << ratio << " (band " << o.rlo << "-" << o.rhi
                               << ") z=" << z);
            BOOST_CHECK_MESSAGE(ratio > o.rlo && ratio < o.rhi,
                                nuc << " cond " << o.name << ": engine/G4 ratio "
                                << ratio << " outside [" << o.rlo << ", " << o.rhi << "]");
        }
    }
}

// --- Level-scheme fix: mutually-exclusive transitions are no longer summed ----
// The per-peak Conditional estimator used to give a peak's MUTUALLY-EXCLUSIVE
// same-level transitions a spurious "partner" probability (its no-pairwise-link
// fallback returned the partner's marginal intensity, assuming independence).
// The level-scheme visit/reach DP gives those partners p=0. These tests exercise
// that on real adapter-built decay data.
//
// NOTE: the Conditional estimator now captures gamma-gamma sum-peak-fed summing-IN
// (e.g. Co-57 122+14->136, via cascade_sum_pair_channels) and coincident vacancy
// K/L x-ray summing (cascade_level_vacancies), so it agrees with FullRealization to
// ~1-2% at contact and within statistics far away. Residual approximations remain
// (x-ray-FED feeding and 3+-gamma partner independence), so FullRealization stays
// the reference. The tests below isolate the exclusive-transition fix (spurious
// same-level summing removed) and the sum-peak-fed summing-IN.

// Marginal partner probability the OLD no-link fallback would have used for a
// same-branch member with no pairwise coincidence link: its marginal intensity.
double old_marginal_fallback(const DecayCascade& dc, std::size_t a, std::size_t b) {
    for (const auto& c : dc.members[a].coincident) if (c.partner == b) return -1.0;  // linked
    for (const auto& c : dc.members[b].coincident) if (c.partner == a) return -1.0;  // linked
    return dc.members[b].intensity;  // no link -> old code summed at marginal I
}

// Scan every ordered gamma pair that leaves the SAME level (mutually exclusive),
// each emitted by exactly one transition, and assert the level-scheme pmate
// zeroes the partner even though the old no-link fallback would have summed it at
// its marginal intensity. Also track the LARGEST such spurious partner removed
// and require it to exceed `min_big_partner`, proving a substantial over-sum
// (not just negligible ones) was eliminated for this nuclide.
void check_exclusive_pairs_zeroed(const std::string& nuc, double min_big_partner) {
    CascadeOptions opt;
    const auto casc = build_cascades(db(), nuc, opt);
    int demonstrated = 0;
    double max_removed = 0.0;
    for (const DecayCascade& dc : casc) {
        if (!dc.level_scheme.valid) continue;
        // Count transitions emitting each member (uniqueness guard, so pmate over
        // this member reflects only the same-level transition).
        std::vector<int> emit_count(dc.members.size(), 0);
        for (const auto& lvl : dc.level_scheme.levels)
            for (const auto& t : lvl.out)
                if (t.gamma_member >= 0
                    && t.gamma_member < static_cast<int>(emit_count.size()))
                    ++emit_count[t.gamma_member];

        for (const CascadeLevel& lvl : dc.level_scheme.levels) {
            std::vector<int> gms;  // gamma members leaving THIS level
            for (const auto& t : lvl.out)
                if (t.gamma_member >= 0 && t.p_gamma > 0.0
                    && emit_count[t.gamma_member] == 1)
                    gms.push_back(t.gamma_member);
            if (gms.size() < 2) continue;
            for (std::size_t ia = 0; ia < gms.size(); ++ia) {
                const std::size_t a = static_cast<std::size_t>(gms[ia]);
                const std::vector<double> pmate =
                    EfficiencyCalculator::cascade_level_pmate(dc, a);
                if (pmate.empty()) continue;
                for (std::size_t ib = 0; ib < gms.size(); ++ib) {
                    if (ib == ia) continue;
                    const std::size_t b = static_cast<std::size_t>(gms[ib]);
                    const double old_p = old_marginal_fallback(dc, a, b);
                    if (old_p <= 0.0) continue;  // linked (never for same-level) / zero
                    BOOST_CHECK_SMALL(pmate[b], 1e-9);  // fixed: exactly 0
                    if (old_p > max_removed) {
                        max_removed = old_p;
                        BOOST_TEST_MESSAGE(nuc << ": condition on "
                            << dc.members[a].energy_keV << " keV -> exclusive partner "
                            << dc.members[b].energy_keV << " keV old fallback p="
                            << old_p << " now pmate=" << pmate[b]);
                    }
                    ++demonstrated;
                }
            }
        }
    }
    BOOST_CHECK_MESSAGE(demonstrated > 0,
                        nuc << ": no mutually-exclusive same-level gamma pair found");
    BOOST_CHECK_MESSAGE(max_removed > min_big_partner,
                        nuc << ": largest spurious partner removed " << max_removed
                        << " <= " << min_big_partner << " (fix not meaningfully exercised)");
}

// Summing factor of `peak_keV` for `nuc` via the Conditional estimator, point
// source `z_cm` behind a 3"x3" NaI front face (no source material).
double conditional_summing_factor(const std::string& nuc, double peak_keV,
                                  uint64_t num_events, double z_cm,
                                  double& unc_out) {
    static Material nai = make_NaI();
    CascadeOptions opt;
    const auto casc = build_cascades(db(), nuc, opt);
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -z_cm));
    CascadeConfig cfg;
    cfg.cascades = casc;
    cfg.method = CascadeMethod::Conditional;
    cfg.num_events = num_events;
    cfg.num_threads = 4;
    cfg.peaks.push_back({peak_keV, 1.5});
    const CascadeResult res = calc.compute_cascade(cfg);
    BOOST_REQUIRE(!res.peaks.empty());
    BOOST_REQUIRE(res.peaks[0].found);
    unc_out = res.peaks[0].summing_factor_unc;
    return res.peaks[0].summing_factor;
}

// Summing factor of `peak_keV` for `nuc` via the FullRealization reference, same
// geometry as conditional_summing_factor (point source `z_cm` behind 3"x3" NaI).
double fullreal_summing_factor(const std::string& nuc, double peak_keV,
                               uint64_t num_events, double z_cm, double& unc_out) {
    static Material nai = make_NaI();
    CascadeOptions opt;
    const auto casc = build_cascades(db(), nuc, opt);
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -z_cm));
    CascadeConfig cfg;
    cfg.cascades = casc;
    cfg.method = CascadeMethod::FullRealization;
    cfg.num_events = num_events;
    cfg.num_threads = 4;
    cfg.peaks.push_back({peak_keV, 1.5});
    const CascadeResult res = calc.compute_cascade(cfg);
    BOOST_REQUIRE(!res.peaks.empty());
    BOOST_REQUIRE(res.peaks[0].found);
    unc_out = res.peaks[0].summing_factor_unc;
    return res.peaks[0].summing_factor;
}

}  // namespace

// ---- Level-scheme exclusivity fix on real decay data -------------------------
// Co-57, Ba-133 and Eu-152 all carry alternative de-excitations of a shared
// level. Verify the level-scheme pmate zeroes those spurious partners (which the
// old marginal fallback summed) on the actual adapter-built cascades.
BOOST_AUTO_TEST_CASE(co57_exclusive_transitions_zeroed) {
    // Conditioning on 136.47 keV, the old fallback summed the 122.06 keV
    // alternative de-excitation of the same level at p~0.856.
    check_exclusive_pairs_zeroed("Co57", /*min_big_partner=*/0.5);
}
BOOST_AUTO_TEST_CASE(ba133_exclusive_transitions_zeroed) {
    // Ba-133 has same-level alternatives summed at up to p~0.62 (356 keV).
    check_exclusive_pairs_zeroed("Ba133", /*min_big_partner=*/0.4);
}
BOOST_AUTO_TEST_CASE(eu152_exclusive_transitions_zeroed) {
    // Eu-152's dense scheme has same-level alternatives summed at up to p~0.46.
    check_exclusive_pairs_zeroed("Eu152", /*min_big_partner=*/0.3);
}

// ---- Co-57: the canonical sum-peak-fed case, at the engine level -------------
// 136.47 keV is the DIRECT-to-ground de-excitation of the 136.47 keV level; the
// 122.06 keV gamma is the ALTERNATIVE (via the 14.41 keV level) de-excitation of
// the SAME level -- mutually EXCLUSIVE with 136.47. But 122.06+14.41 = 136.47, so
// that path SUM-FEEDS the 136 window without emitting the 136 gamma: net summing-
// IN. Pre-fix history: the old fallback over-summed the exclusive 122 partner
// (SF~0.638); the exclusivity fix took it to ~1.0; the sum-peak-fed channels take
// it to ~1.18, matching the FullRealization reference (1.20) to ~2%. The residual
// is x-ray-FED feeding (122 + coincident Fe Kα x-rays), documented as out of scope.
BOOST_AUTO_TEST_CASE(co57_136_sum_peak_fed_summing_in) {
    double uc = 0.0, uf = 0.0;
    const double sc = conditional_summing_factor("Co57", 136.47, 2000000, 0.5, uc);
    const double sf = fullreal_summing_factor("Co57", 136.47, 4000000, 0.5, uf);
    const double z = (sf - sc) / std::sqrt(uc * uc + uf * uf);
    BOOST_TEST_MESSAGE("Co57 SF(136.47) @ 0.5 cm: Conditional " << sc << " +/- " << uc
                       << " vs FullRealization " << sf << " +/- " << uf
                       << " (z=" << z << ", pre-fix ~0.638/1.0)");
    // Net summing-IN is now captured (was ~1.0): SF well above 1.1, close to FR.
    BOOST_CHECK_GT(sc, 1.12);
    BOOST_CHECK_LT(sc, 1.26);
    // Agreement with the reference to ~3% (the x-ray-fed residual keeps it just
    // outside pure statistics at these high stats).
    BOOST_CHECK_LT(std::abs(sc - sf) / sf, 0.03);
}

// Far geometry: the sum-fed channel term is ~eps^2-suppressed, so the Conditional
// (cone-biased, high stats) must agree with the FullRealization reference within
// statistics -- no regression to the far-geometry accuracy the Conditional buys.
BOOST_AUTO_TEST_CASE(co57_136_far_geometry_matches_fullreal) {
    double uc = 0.0, uf = 0.0;
    const double sc = conditional_summing_factor("Co57", 136.47, 2000000, 10.0, uc);
    const double sf = fullreal_summing_factor("Co57", 136.47, 8000000, 10.0, uf);
    const double z = (sf - sc) / std::sqrt(uc * uc + uf * uf);
    BOOST_TEST_MESSAGE("Co57 SF(136.47) @ 10 cm: Conditional " << sc << " +/- " << uc
                       << " vs FullRealization " << sf << " +/- " << uf << " (z=" << z << ")");
    BOOST_CHECK_LT(std::abs(z), 4.0);   // within statistics
    BOOST_CHECK_GT(sc, 1.0);            // still net summing-IN
}

// Ba-133 356 keV: a net summing-OUT peak dominated by coincident gamma + K x-ray
// summing, with three small sum-fed pairs (81+276, 79.6+276, 302+53). The old
// Conditional (no vacancy x-rays) gave ~0.84; adding the coincident K-vacancy
// x-ray summing-OUT plus the sum-fed channels takes it to ~0.52, near the
// FullRealization reference (~0.49). Loose band: the residual is the multi-gamma
// partner-independence approximation + x-ray-fed feeding.
BOOST_AUTO_TEST_CASE(ba133_356_xray_and_sum_fed) {
    double uc = 0.0, uf = 0.0;
    const double sc = conditional_summing_factor("Ba133", 356.013, 2000000, 0.5, uc);
    const double sf = fullreal_summing_factor("Ba133", 356.013, 4000000, 0.5, uf);
    BOOST_TEST_MESSAGE("Ba133 SF(356) @ 0.5 cm: Conditional " << sc << " +/- " << uc
                       << " vs FullRealization " << sf << " +/- " << uf << " (was ~0.84)");
    // Captures the K-x-ray summing-OUT the old Conditional missed (was ~0.84).
    BOOST_CHECK_LT(sc, 0.60);
    BOOST_CHECK_GT(sc, 0.44);
    BOOST_CHECK_LT(std::abs(sc - sf) / sf, 0.10);  // within ~10% of the reference
}

// The G4 reference numbers, windows, ratio bands, per-nuclide source_cm and emax
// all live in tests/data/geant4_reference/cascade_summing_multi.csv (the single
// source of truth, also consumed by profiling/compare_validation.py). Each case
// below just names the nuclide + the low-stats event budget; check_nuclide_from_csv
// gates every full row (FullRealization area) and every cond row (Conditional
// per-decay area) for that nuclide.

// ---- Co-60: pure gamma-gamma summing (no x-rays) -----------------------------
// 1173 + 1332 cascade; the 2505 sum peak is pure summing-IN, the 1173/1332 peaks
// carry the summing-OUT suppression (also gated via the Conditional estimator).
BOOST_AUTO_TEST_CASE(co60_gamma_gamma_summing) {
    check_nuclide_from_csv("Co60", 2000000, 2000000);
}

// ---- Ba-133: gamma-gamma + K x-ray summing (K x-ray 31, 81/356 gammas, sum) --
BOOST_AUTO_TEST_CASE(ba133_kxray_summing) {
    check_nuclide_from_csv("Ba133", 2000000, 2000000);
}

// ---- Am-241: L x-ray summing (per-subshell + Coster-Kronig). The sum59+L band's
// upper bound guards the level path: an independent-vacancy regression jumps it ~3x.
BOOST_AUTO_TEST_CASE(am241_lxray_summing) {
    check_nuclide_from_csv("Am241", 4000000, 4000000);
}

// ---- I-125: lowest-energy K x-ray summing (EC -> Te; Te Kx ~27, 35.49+Kx) ----
BOOST_AUTO_TEST_CASE(i125_kxray_summing) {
    check_nuclide_from_csv("I125", 2000000, 2000000);
}

// ---- Eu-152: mixed EC + beta-, dense multi-gamma cascade (wider bands: the
// engine/G4 spread is SandiaDecay-vs-G4 emission-intensity, not summing). --------
BOOST_AUTO_TEST_CASE(eu152_mixed_mode_summing) {
    check_nuclide_from_csv("Eu152", 2000000, 2000000);
}

// ---- Na-22: beta+ 511 summing (source_cm=0.15 in the CSV so the positron
// annihilates near the source). P511 lower band guards the Annih511 emission. -----
BOOST_AUTO_TEST_CASE(na22_positron_511_summing) {
    check_nuclide_from_csv("Na22", 2000000, 2000000);
}

// ---- Co-57: G4-gated (8M ion decay). 122.06 & 136.47 are mutually exclusive
// (no 122+136 sum): P122 = summing-OUT, P136 = summing-IN (sum-fed by 122.06+14.41),
// P14 = summing-OUT; P122 also carries a Conditional row. Complements the
// co57_136_* Conditional-vs-FullRealization tests above (this one anchors to G4). ---
BOOST_AUTO_TEST_CASE(co57_gamma_summing_vs_g4) {
    check_nuclide_from_csv("Co57", 2000000, 2000000);
}

// ===== IC-electron summing (compute_cascade enable_ic_electrons) ==============
// These guard the internal-conversion electron deposition feature. They do NOT
// need a GEANT4 reference: they check (a) the DISTANCE GATE -- a distant source
// is unchanged, so the feature is safe to ship; and (b) that at CONTACT the
// coincident IC electrons produce the large, measured extra summing-OUT (the
// physics the photon-only path was missing). The magnitudes match the Phase-0
// study probe (the dev-only IC-electron probe). Run under BOTH estimators.
namespace {

// Summing factor of one peak of `nuc` on a bare 3"x3" NaI point source at
// `dist_cm`, with the IC-electron flag on/off, in the given method.
double ic_summing_factor(const std::string& nuc, double peak_keV, double dist_cm,
                         bool ic_on, CascadeMethod method, uint64_t events,
                         double& unc) {
    static Material nai = make_NaI();
    const auto casc = build_cascades(db(), nuc);
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -dist_cm));
    CascadeConfig cfg;
    cfg.cascades = casc;
    cfg.peaks.push_back(PeakWindow{peak_keV, 1.5});
    cfg.method = method;
    cfg.num_events = events;
    cfg.num_threads = 4;
    cfg.enable_ic_electrons = ic_on;
    const CascadeResult res = calc.compute_cascade(cfg);
    unc = res.peaks.empty() ? 0.0 : res.peaks[0].summing_factor_unc;
    return res.peaks.empty() ? 0.0 : res.peaks[0].summing_factor;
}

}  // namespace

// The distance gate: a source 10 cm from the face is far enough that the IC
// electrons (isotropic, short-range vs the solid angle) cannot reach the crystal,
// so the summing factor is statistically unchanged by the flag. This is the
// safety guarantee that lets the feature ship without perturbing distant sources.
BOOST_AUTO_TEST_CASE(ic_electron_distant_source_invariant) {
    for (CascadeMethod m : {CascadeMethod::FullRealization, CascadeMethod::Conditional}) {
        double uoff = 0.0, uon = 0.0;
        const double off = ic_summing_factor("Ba133", 356.01, 10.0, false, m, 3000000, uoff);
        const double on  = ic_summing_factor("Ba133", 356.01, 10.0, true,  m, 3000000, uon);
        const double sig = std::sqrt(uoff * uoff + uon * uon);
        const double z = sig > 0.0 ? (on - off) / sig : 0.0;
        BOOST_TEST_MESSAGE("Ba133 356 @10cm " << (m == CascadeMethod::Conditional ? "Cond" : "FR")
                           << ": SF_off=" << off << " SF_on=" << on
                           << " delta=" << (on - off) << " z=" << z);
        // Distant: flag must not move the summing factor by more than ~0.5%.
        BOOST_CHECK_MESSAGE(std::abs(on - off) < 0.006,
                            "Ba133 356 @10cm: IC flag changed SF by " << (on - off)
                            << " (expected ~0 -- distance gate failed)");
    }
}

// At CONTACT the coincident IC electrons deposit real energy on top of the peak
// gamma -> extra summing-OUT that the photon-only path misses. Measured (probe):
// Ba-133 356 SF drops from ~0.44 to ~0.31 (FR) / ~0.47 to ~0.36 (Cond). Guard the
// feature is active and physical (a large, correctly-signed drop) in both methods.
BOOST_AUTO_TEST_CASE(ic_electron_contact_extra_summing_out) {
    for (CascadeMethod m : {CascadeMethod::FullRealization, CascadeMethod::Conditional}) {
        double uoff = 0.0, uon = 0.0;
        const double off = ic_summing_factor("Ba133", 356.01, 0.1, false, m, 3000000, uoff);
        const double on  = ic_summing_factor("Ba133", 356.01, 0.1, true,  m, 3000000, uon);
        const double sig = std::sqrt(uoff * uoff + uon * uon);
        const double z = sig > 0.0 ? (on - off) / sig : 0.0;
        BOOST_TEST_MESSAGE("Ba133 356 @contact " << (m == CascadeMethod::Conditional ? "Cond" : "FR")
                           << ": SF_off=" << off << " SF_on=" << on
                           << " delta=" << (on - off) << " z=" << z);
        // The IC electrons must drop the summing factor substantially (>~5%),
        // and only downward (extra summing-out, never summing-in for this peak).
        BOOST_CHECK_MESSAGE(off - on > 0.05,
                            "Ba133 356 @contact: IC flag reduced SF by only " << (off - on)
                            << " (expected a large summing-out drop)");
        BOOST_CHECK_MESSAGE(on > 0.15 && on < off,
                            "Ba133 356 @contact: SF_on=" << on << " out of physical range");
    }
}

namespace {
// Per-decay summed-spectrum peak areas for a BARE point source (no source
// material) on a 3"x3" NaI at distance z, with the IC-electron flag on/off.
// Matches the geometry of the dev-only GEANT4 IC-electron cross-check (
// g4_findings.md): bare point, so the conversion electrons are unabsorbed and
// the channel is maximal/isolated.
std::vector<float> bare_point_spectrum(const std::string& nuc, double z_cm,
                                       bool ic_on, int emax, uint64_t events,
                                       double al_r_cm = 0.0) {
    static Material nai = make_NaI();
    static Material al = make_Aluminum();
    const auto casc = build_cascades(db(), nuc);
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    if (al_r_cm > 0.0) {  // solid Al cylinder source (self-attenuates; routes IC e- through the walk)
        calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -z_cm), al_r_cm, al_r_cm);
        calc.set_source_material(&al);
    } else {
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -z_cm));
    }
    CascadeConfig cfg;
    cfg.cascades = casc;
    cfg.method = CascadeMethod::FullRealization;
    cfg.num_events = events;
    cfg.num_threads = 4;
    cfg.enable_ic_electrons = ic_on;
    for (int e = 0; e <= emax; ++e) cfg.spectrum_bin_edges.push_back(static_cast<float>(e));
    return calc.compute_cascade(cfg).summed_spectrum;
}
}  // namespace

// GEANT4-anchored MAGNITUDE gate for the IC-electron channel. The GEANT4
// reference is a 1.5M-decay Ba-133 ion run (ceelo_g4val --histogram
// --correlated-gamma) on the SAME bare-point 2 cm geometry, in an AIR world so
// G4's conversion electrons air-degrade like the engine's (see g4_findings.md;
// use an air world + /process/had/rdm/thresholdForVeryLongDecayTime or the ion
// never decays). Validated result: flag-ON matches these G4 areas to <1% on the
// IC-sensitive peaks (Kx31/P356/sum356+81), where flag-OFF was +7-13% high. This
// locks that magnitude in: flag-ON must sit in the band AND be closer to G4 than
// flag-OFF on the IC-sensitive peaks.
BOOST_AUTO_TEST_CASE(ic_electron_g4_magnitude_ba133_2cm) {
    const std::vector<float> on  = bare_point_spectrum("Ba133", 2.0, true,  500, 3000000);
    const std::vector<float> off = bare_point_spectrum("Ba133", 2.0, false, 500, 3000000);
    struct Ref { const char* name; int lo, hi; double g4; double rlo, rhi; bool ic_sensitive; };
    // G4 (1.5M, 2 cm, AIR world): per-decay area, and the accepted ON/G4 band.
    const std::vector<Ref> refs = {
        {"Kx31",       29,  33, 1.460387e-1, 0.95, 1.04, true},   // ON 0.99
        {"P356",      353, 359, 5.337600e-2, 0.95, 1.06, true},   // ON 1.00
        {"sum356+81", 434, 440, 6.667333e-3, 0.88, 1.14, true},   // ON 1.01 (1% G4 stat)
        {"P81",        79,  83, 5.645200e-2, 0.97, 1.08, false},  // ON 1.03 -- IC-independent
    };
    for (const Ref& r : refs) {
        const double aon = area(on, r.lo, r.hi);
        const double aoff = area(off, r.lo, r.hi);
        BOOST_TEST_MESSAGE("Ba133 " << r.name << " @2cm bare/air: ON/G4=" << (aon / r.g4)
                           << " OFF/G4=" << (aoff / r.g4) << " (band " << r.rlo << "-" << r.rhi << ")");
        BOOST_CHECK_MESSAGE(aon / r.g4 > r.rlo && aon / r.g4 < r.rhi,
                            "Ba133 " << r.name << " @2cm: flag-ON/G4 " << (aon / r.g4)
                            << " outside [" << r.rlo << ", " << r.rhi << "]");
        // On the IC-sensitive peaks, turning the flag ON must move the engine
        // TOWARD G4 (flag-OFF over-predicts by missing the summing-out).
        if (r.ic_sensitive)
            BOOST_CHECK_MESSAGE(std::abs(aon - r.g4) < std::abs(aoff - r.g4),
                                "Ba133 " << r.name << " @2cm: flag-ON (" << aon
                                << ") not closer to G4 (" << r.g4 << ") than flag-OFF (" << aoff << ")");
    }
}

// Source-material branch (walk_in_source_geometry): a solid Al source cylinder
// ABSORBS the low-energy conversion electrons (a 45 keV e- range in Al ~ 0.3 mm ~
// the cylinder), so with source material the IC channel is largely contained --
// the flag-ON peak areas stay close to flag-OFF, unlike the bare-point case at
// the same distance where the effect is ~10x larger. GEANT4-validated to ~2-3%
// (Ba-133 2cm, 0.1 cm Al, air world; see g4_findings.md). This guards branch (a)
// against a regression that spuriously lets the electrons escape/deposit (which
// would drive the ON areas well below OFF). Self-contained (no G4 in the test).
BOOST_AUTO_TEST_CASE(ic_electron_source_material_contained) {
    const double z = 2.0, al_r = 0.1;
    const std::vector<float> on  = bare_point_spectrum("Ba133", z, true,  500, 3000000, al_r);
    const std::vector<float> off = bare_point_spectrum("Ba133", z, false, 500, 3000000, al_r);
    // Contrast: the SAME nuclide/distance with a BARE point (no absorbing material).
    const std::vector<float> bp_on  = bare_point_spectrum("Ba133", z, true,  500, 3000000);
    const std::vector<float> bp_off = bare_point_spectrum("Ba133", z, false, 500, 3000000);
    const double al_frac = std::abs(area(on, 353, 359) - area(off, 353, 359)) / area(off, 353, 359);
    const double bp_frac = std::abs(area(bp_on, 353, 359) - area(bp_off, 353, 359)) / area(bp_off, 353, 359);
    BOOST_TEST_MESSAGE("Ba133 P356 @2cm IC |on-off|/off: Al-source=" << al_frac
                       << " vs bare-point=" << bp_frac);
    // With source material the residual IC effect on P356 is small (electrons
    // contained) and much smaller than the bare-point effect.
    BOOST_CHECK_MESSAGE(al_frac < 0.02,
                        "Ba133 P356 @2cm Al source: IC changed the peak by " << al_frac
                        << " (>2% -- branch (a) may be over-escaping)");
    BOOST_CHECK_MESSAGE(al_frac < 0.4 * bp_frac,
                        "Ba133 P356 @2cm: Al-source IC effect " << al_frac
                        << " not much smaller than bare-point " << bp_frac
                        << " (source material should absorb the electrons)");
}

// The Conditional estimator's IC wiring (prim/independent/sum-pair vacancies +
// K-Auger) must stay consistent with the FullRealization reference. Both should
// show the same large IC summing-OUT on Ba-133 356 at contact. Guards the newer,
// more intricate Conditional path against silent drift from FR.
BOOST_AUTO_TEST_CASE(ic_electron_estimator_consistency) {
    double uf = 0.0, uc = 0.0;
    const double fr   = ic_summing_factor("Ba133", 356.01, 0.1, true,
                                          CascadeMethod::FullRealization, 3000000, uf);
    const double cond = ic_summing_factor("Ba133", 356.01, 0.1, true,
                                          CascadeMethod::Conditional, 3000000, uc);
    BOOST_TEST_MESSAGE("Ba133 356 @contact IC-on: FR SF=" << fr << " Conditional SF=" << cond);
    // Both must show the strong summing-out (SF well below 1) and agree to ~10%
    // (residual FR-vs-Conditional modeling differences are documented, ~1-2% at
    // this geometry; 0.10 is a loose guard against a wiring regression).
    BOOST_CHECK_MESSAGE(fr < 0.5 && cond < 0.5,
                        "Ba133 356 @contact IC-on: expected strong summing-out in both "
                        "(FR=" << fr << ", Cond=" << cond << ")");
    BOOST_CHECK_MESSAGE(std::abs(fr - cond) < 0.10,
                        "Ba133 356 @contact IC-on: FR (" << fr << ") and Conditional ("
                        << cond << ") disagree by more than 0.10");
}

BOOST_AUTO_TEST_CASE(residual_transition_full_conditional_analytic_agree) {
    DecayCascade dc;
    dc.branch_weight = 1.0;
    dc.members = {
        {500.0, CascadeParticleType::Gamma, 1.0},
        {200.0, CascadeParticleType::Gamma, 0.35}
    };
    dc.level_scheme.levels.resize(2);
    dc.level_scheme.levels[1].feeding = 1.0;
    dc.level_scheme.levels[1].out.push_back(
        {0, 0, 500.0, 1.0, 1.0, 0, 0, 0, 0});
    dc.level_scheme.valid = true;
    dc.residual_transitions.push_back(
        {1, 200.0, 0.65, 0.35, 0, 0, 0, 0, 0, 0, 0});

    static Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -2.0));

    CascadeConfig cfg;
    cfg.cascades = {dc};
    cfg.peaks = {{500.0, 1.5}};
    cfg.num_events = 400000;
    cfg.num_threads = 1;
    cfg.method = CascadeMethod::FullRealization;
    const auto fr = calc.compute_cascade(cfg);
    cfg.method = CascadeMethod::Conditional;
    const auto cond = calc.compute_cascade(cfg);
    BOOST_REQUIRE(fr.peaks[0].found);
    BOOST_REQUIRE(cond.peaks[0].found);

    std::map<double, EfficiencyResult> cache;
    cache[500.0] = calc.compute(500.0, 250000, 1);
    cache[200.0] = calc.compute(200.0, 250000, 1);
    const CachedEfficiencyProvider provider(cache);
    const auto an = compute_cascade_analytic({dc}, {{500.0, 1.5}}, provider);
    BOOST_REQUIRE(an[0].found);

    const double fc_sig = std::hypot(fr.peaks[0].summing_factor_unc,
                                     cond.peaks[0].summing_factor_unc);
    const double fa_sig = std::hypot(fr.peaks[0].summing_factor_unc,
                                     an[0].summing_factor_unc);
    const double z_fc = fc_sig > 0.0
        ? (cond.peaks[0].summing_factor - fr.peaks[0].summing_factor) / fc_sig : 0.0;
    const double z_fa = fa_sig > 0.0
        ? (an[0].summing_factor - fr.peaks[0].summing_factor) / fa_sig : 0.0;
    BOOST_TEST_MESSAGE("residual synthetic: Full=" << fr.peaks[0].summing_factor
        << " +/- " << fr.peaks[0].summing_factor_unc
        << " Conditional=" << cond.peaks[0].summing_factor
        << " +/- " << cond.peaks[0].summing_factor_unc << " z=" << z_fc
        << " Analytic=" << an[0].summing_factor
        << " +/- " << an[0].summing_factor_unc << " z=" << z_fa);
    BOOST_CHECK_LT(std::abs(z_fc), 4.0);
    BOOST_CHECK_LT(std::abs(z_fa), 4.0);
}

BOOST_AUTO_TEST_CASE(invalid_branch_completeness_flags_all_estimators) {
    auto single = [](double energy, bool valid, double weight) {
        DecayCascade dc;
        dc.branch_weight = weight;
        dc.members = {{energy, CascadeParticleType::Gamma, 1.0}};
        if (valid) {
            dc.level_scheme.levels.resize(2);
            dc.level_scheme.levels[1].feeding = 1.0;
            dc.level_scheme.levels[1].out.push_back(
                {0, 0, energy, 1.0, 1.0, 0, 0, 0, 0});
            dc.level_scheme.valid = true;
        }
        return dc;
    };
    DecayCascade invalid_feeder;
    invalid_feeder.branch_weight = 0.5;
    invalid_feeder.members = {
        {40.0, CascadeParticleType::Gamma, 0.40},
        {60.0, CascadeParticleType::Gamma, 0.30}
    };
    invalid_feeder.members[0].coincident.push_back({1, 0.75});
    invalid_feeder.members[1].coincident.push_back({0, 1.00});
    const std::vector<DecayCascade> cascades = {
        single(100.0, true, 1.0), single(120.0, true, 1.0),
        invalid_feeder, single(200.0, false, 0.8)};

    static Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -2.0));
    CascadeConfig cfg;
    cfg.cascades = cascades;
    cfg.peaks = {{100.0, 0.1}, {120.0, 0.1}, {200.0, 0.1}};
    cfg.num_events = 2000;
    cfg.num_threads = 1;

    for (CascadeMethod method : {CascadeMethod::Conditional,
                                 CascadeMethod::FullRealization}) {
        cfg.method = method;
        const auto result = calc.compute_cascade(cfg);
        BOOST_REQUIRE_EQUAL(result.peaks.size(), 3u);
        BOOST_REQUIRE(result.peaks[0].found && result.peaks[1].found
                      && result.peaks[2].found);
        BOOST_CHECK(!result.peaks[0].summing_model_complete);
        BOOST_CHECK(result.peaks[1].summing_model_complete);
        BOOST_CHECK(!result.peaks[2].summing_model_complete);
    }
}
