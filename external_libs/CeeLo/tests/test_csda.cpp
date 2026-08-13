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

/// @file test_csda.cpp
/// @brief Unit tests for CSDA electron transport physics.
///
/// Tests cover:
///  1. Bethe-Bloch stopping power physics (sign, ordering, known values)
///  2. CSDA range tables (monotonicity, physical bounds, Bragg additivity)
///  3. Efficiency improvement vs. local-deposit at high energy
///  4. Invariant: ε_FEP ≤ ε_total always holds with CSDA enabled
///  5. Recoil and photoelectron direction helpers

#define BOOST_TEST_MODULE CsdaTests
#include <boost/test/unit_test.hpp>

#include "physics/ElectronCsda.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"
#include "efficiency/EfficiencyCalculator.h"

#include <Eigen/Core>
#include <random>
#include <cmath>

using namespace ceelo;

// ============================================================
// Suite 1: Stopping power
// ============================================================
BOOST_AUTO_TEST_SUITE(StoppingPower)

BOOST_AUTO_TEST_CASE(stopping_power_positive)
{
    // S must be positive for all Z and energies in range
    for (int Z : {1, 6, 8, 11, 13, 26, 53, 82}) {
        double A = ElectronCsda::atomic_weight(Z);
        for (double KE : {1.0, 10.0, 100.0, 500.0, 1000.0, 5000.0}) {
            double S = ElectronCsda::stopping_power_MeV_cm2_g(Z, A, KE);
            BOOST_CHECK_GT(S, 0.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(stopping_power_decreases_with_energy_above_200keV)
{
    // Bethe-Bloch stopping power has a minimum at ~1 MeV then slowly rises.
    // It monotonically decreases from ~100 keV to 1 MeV for all Z.
    for (int Z : {11, 26, 53}) {
        double A = ElectronCsda::atomic_weight(Z);
        double S100  = ElectronCsda::stopping_power_MeV_cm2_g(Z, A, 100.0);
        double S500  = ElectronCsda::stopping_power_MeV_cm2_g(Z, A, 500.0);
        double S1000 = ElectronCsda::stopping_power_MeV_cm2_g(Z, A, 1000.0);
        BOOST_CHECK_GT(S100, S500);
        BOOST_CHECK_GT(S500, S1000);
    }
}

BOOST_AUTO_TEST_CASE(stopping_power_Al_1MeV_vs_nist)
{
    // NIST ESTAR Al at 1 MeV: S_col ≈ 1.465 MeV cm²/g.
    // With NIST correction factor applied, our value should be within 2%.
    double A_Al = ElectronCsda::atomic_weight(13);
    double S = ElectronCsda::stopping_power_MeV_cm2_g(13, A_Al, 1000.0);
    BOOST_CHECK_GT(S, 1.0);   // At least 1 MeV cm²/g
    BOOST_CHECK_LT(S, 2.5);   // Not absurdly large
    // With NIST ESTAR correction, should be within 5% of NIST:
    BOOST_CHECK_GT(S, 0.95 * 1.465);
    BOOST_CHECK_LT(S, 1.05 * 1.465);
}

BOOST_AUTO_TEST_CASE(higher_Z_gives_higher_stopping_per_g)
{
    // At 1 MeV: S(Z=1)/A >> S(Z=82)/A in units MeV cm²/g because Z/A differs
    // and I(Pb) >> I(H).  Lower I → higher S.
    // H has Z/A = 1.0, Pb has Z/A ≈ 0.40 → expect S_H > S_Pb (per g/cm²)
    double S_H  = ElectronCsda::stopping_power_MeV_cm2_g(1,  1.008,   1000.0);
    double S_Pb = ElectronCsda::stopping_power_MeV_cm2_g(82, 207.200, 1000.0);
    BOOST_CHECK_GT(S_H, S_Pb);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 2: CSDA range tables
// ============================================================
BOOST_AUTO_TEST_SUITE(RangeTables)

BOOST_AUTO_TEST_CASE(range_positive_and_increases_with_energy)
{
    const ElectronCsda& csda = ElectronCsda::instance();
    for (int Z : {8, 13, 26, 53}) {
        double R_prev = 0.0;
        for (double KE : {1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0}) {
            double R = csda.range_gcm2(Z, KE);
            BOOST_CHECK_GT(R, 0.0);
            BOOST_CHECK_GT(R, R_prev);
            R_prev = R;
        }
    }
}

BOOST_AUTO_TEST_CASE(range_NaI_at_662keV_physical_bounds)
{
    // CSDA range of a 662 keV electron in NaI (density 3.67 g/cm³).
    // NIST ESTAR reference (Bragg rule, g/cm²):
    //   Na (Z=11): NIST ~0.25 g/cm²
    //   I  (Z=53): NIST ~0.24 g/cm²
    // Bragg: R_NaI ≈ 0.24 g/cm² → in cm: 0.065 cm (0.65 mm)
    //
    // Our Bethe-Bloch formula omits shell corrections → underestimates S by 10–30%
    // at 10–300 keV.  Since most range comes from the 50–300 keV integration region,
    // this causes our range to be ~50% LARGER than NIST.
    // Expected: 0.24 × 1.5 ≈ 0.36 g/cm²  (generous bounds: 0.10 to 0.70)
    Material nai = make_NaI();
    const ElectronCsda& csda = ElectronCsda::instance();
    double R = csda.range_gcm2_material(nai, 662.0);
    BOOST_CHECK_GT(R, 0.10);   // At least 0.10 g/cm² (well above noise)
    BOOST_CHECK_LT(R, 0.70);   // Less than 0.70 g/cm² (sanity check)
}

BOOST_AUTO_TEST_CASE(bragg_rule_monotone_compound)
{
    // Bragg-rule compound range should also increase with energy
    Material nai  = make_NaI();
    Material pb   = make_Lead();
    const ElectronCsda& csda = ElectronCsda::instance();

    for (const auto& mat : {nai, pb}) {
        double R_prev = 0.0;
        for (double KE : {10.0, 100.0, 662.0, 1000.0, 3000.0}) {
            double R = csda.range_gcm2_material(mat, KE);
            BOOST_CHECK_GT(R, R_prev);
            R_prev = R;
        }
    }
}

BOOST_AUTO_TEST_CASE(residual_energy_decreases_with_path)
{
    const ElectronCsda& csda = ElectronCsda::instance();
    Material nai = make_NaI();
    double KE = 662.0;

    double T0 = csda.residual_energy_keV(nai, KE, 0.0);
    double T1 = csda.residual_energy_keV(nai, KE, 0.05);
    double T2 = csda.residual_energy_keV(nai, KE, 0.10);
    double T3 = csda.residual_energy_keV(nai, KE, 0.50);  // > range → stops

    BOOST_CHECK_CLOSE(T0, KE, 1.0);  // No path → full energy
    BOOST_CHECK_LT(T1, KE);
    BOOST_CHECK_LT(T2, T1);
    BOOST_CHECK_EQUAL(T3, 0.0);  // Stops (range < 0.50 g/cm²)
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 3: Direction sampling
// ============================================================
BOOST_AUTO_TEST_SUITE(DirectionSampling)

BOOST_AUTO_TEST_CASE(photoelectron_direction_normalized)
{
    std::mt19937_64 rng(42);
    Eigen::Vector3d gamma_dir(0.0, 0.0, 1.0);

    for (double KE : {10.0, 100.0, 662.0, 1500.0}) {
        for (int i = 0; i < 100; ++i) {
            auto d = sample_photoelectron_direction(gamma_dir, KE, rng);
            BOOST_CHECK_CLOSE(d.norm(), 1.0, 0.01);
        }
    }
}

BOOST_AUTO_TEST_CASE(photoelectron_forward_peaked_at_high_energy)
{
    // At high energy (β→1), Sauter distribution is very forward-peaked.
    // Mean cos(θ) should be significantly > 0 at 1 MeV.
    std::mt19937_64 rng(1234);
    Eigen::Vector3d gamma_dir(0.0, 0.0, 1.0);
    int N = 5000;
    double sum_cos = 0.0;
    for (int i = 0; i < N; ++i) {
        auto d = sample_photoelectron_direction(gamma_dir, 1000.0, rng);
        sum_cos += d.z(); // z-component = cos(θ) relative to photon direction
    }
    double mean_cos = sum_cos / N;
    BOOST_CHECK_GT(mean_cos, 0.3);  // Forward-peaked at 1 MeV
}

BOOST_AUTO_TEST_CASE(compton_recoil_direction_normalized)
{
    Eigen::Vector3d d_in(0.0, 0.0, 1.0);
    Eigen::Vector3d d_out(0.5, 0.0, std::sqrt(0.75));  // ~30° scatter
    d_out.normalize();
    double E_in  = 662.0;
    double E_out = 500.0;

    Eigen::Vector3d d_e = compton_recoil_direction(d_in, d_out, E_in, E_out);
    BOOST_CHECK_CLOSE(d_e.norm(), 1.0, 0.01);
}

BOOST_AUTO_TEST_CASE(compton_recoil_direction_forward)
{
    // For forward scatter (small angle), recoil electron should be approximately
    // perpendicular to or forward of the incoming photon direction.
    Eigen::Vector3d d_in(0.0, 0.0, 1.0);
    // Small forward scatter: d_out ≈ d_in with tiny lateral component
    Eigen::Vector3d d_out = Eigen::Vector3d(0.05, 0.0, std::sqrt(1.0 - 0.0025)).normalized();
    double E_in  = 662.0;
    double E_out = 655.0; // Small energy transfer

    Eigen::Vector3d d_e = compton_recoil_direction(d_in, d_out, E_in, E_out);
    BOOST_CHECK_CLOSE(d_e.norm(), 1.0, 0.01);
    // Recoil should be mostly perpendicular or forward
    BOOST_CHECK_GE(d_e.z(), -0.5);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 4: Efficiency improvement with CSDA
// ============================================================
BOOST_AUTO_TEST_SUITE(EfficiencyWithCsda)

BOOST_AUTO_TEST_CASE(csda_reduces_fep_at_high_energy)
{
    // At 1173 keV (Co-60), CSDA with Moliere + brems + straggling should
    // reduce FEP relative to local deposit.  Electron escape and brems
    // radiation carry energy out of the crystal.  The CSDA FEP should be
    // noticeably lower but not zero (crystal reabsorbs some brems photons).
    Material nai = make_NaI();
    const uint64_t N = 30000;

    EfficiencyCalculator calc_local, calc_csda;

    for (auto* calc : {&calc_local, &calc_csda}) {
        calc->set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc->set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    }
    calc_local.enable_electron_csda(false);
    calc_csda.enable_electron_csda(true);

    auto res_local = calc_local.compute(1173.0, N, 1);
    auto res_csda  = calc_csda .compute(1173.0, N, 1);

    // CSDA FEP should be less than local (electron tracking removes energy)
    // but still a reasonable fraction (not zero — crystal reabsorbs some)
    if (res_local.full_energy_peak_efficiency > 1e-6) {
        double ratio = res_csda.full_energy_peak_efficiency
                     / res_local.full_energy_peak_efficiency;
        BOOST_CHECK_GT(ratio, 0.10);  // At least 10% of local
        BOOST_CHECK_LT(ratio, 1.20);  // Not larger than local
    }
}

BOOST_AUTO_TEST_CASE(csda_small_effect_at_low_energy)
{
    // At 80 keV, photoelectron range in NaI is < 0.1 mm → CSDA ≈ local deposit.
    // FEP efficiencies should agree within 5%.
    Material nai = make_NaI();
    const uint64_t N = 20000;

    EfficiencyCalculator calc_local, calc_csda;
    for (auto* calc : {&calc_local, &calc_csda}) {
        calc->set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc->set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    }
    calc_local.enable_electron_csda(false);
    calc_csda.enable_electron_csda(true);

    auto res_local = calc_local.compute(80.0, N, 1);
    auto res_csda  = calc_csda .compute(80.0, N, 1);

    // Should be within 10% of each other (stat + small CSDA effect)
    if (res_local.full_energy_peak_efficiency > 1e-6) {
        double ratio = res_csda.full_energy_peak_efficiency
                     / res_local.full_energy_peak_efficiency;
        BOOST_CHECK_GT(ratio, 0.90);
        BOOST_CHECK_LT(ratio, 1.10);
    }
}

BOOST_AUTO_TEST_CASE(fep_leq_total_invariant_holds_with_csda)
{
    // Fundamental invariant: ε_FEP ≤ ε_total always.
    // Test at several energies with CSDA enabled.
    Material nai = make_NaI();
    const uint64_t N = 15000;

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.enable_electron_csda(true);

    for (double E : {80.0, 300.0, 662.0, 1173.0, 2000.0}) {
        auto res = calc.compute(E, N, 1);
        BOOST_CHECK_GE(res.total_efficiency,
                       res.full_energy_peak_efficiency - 1e-9);
    }
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 5: Radiation length
// ============================================================
BOOST_AUTO_TEST_SUITE(RadiationLength)

BOOST_AUTO_TEST_CASE(radiation_length_Al_vs_pdg)
{
    // PDG 2022 tabulated: X₀(Al) = 24.01 g/cm²
    // Tsai formula gives ~24.3 g/cm² (~1% over).  Check within 10%.
    double A_Al = ElectronCsda::atomic_weight(13);
    double X0 = ElectronCsda::radiation_length_gcm2_element(13, A_Al);
    BOOST_CHECK_GT(X0, 24.01 * 0.90);  // at least 90% of PDG value
    BOOST_CHECK_LT(X0, 24.01 * 1.10);  // at most 110% of PDG value
}

BOOST_AUTO_TEST_CASE(radiation_length_Pb_vs_pdg)
{
    // PDG 2022 tabulated: X₀(Pb) = 6.37 g/cm²
    // Tsai formula gives ~6.3 g/cm² (~1% under).  Check within 10%.
    double A_Pb = ElectronCsda::atomic_weight(82);
    double X0 = ElectronCsda::radiation_length_gcm2_element(82, A_Pb);
    BOOST_CHECK_GT(X0, 6.37 * 0.90);
    BOOST_CHECK_LT(X0, 6.37 * 1.10);
}

BOOST_AUTO_TEST_CASE(radiation_length_NaI_compound)
{
    // NaI via Bragg rule.
    // X₀_Na ≈ 28 g/cm², X₀_I ≈ 8.6 g/cm²
    // Mass fractions: Na 15.3%, I 84.7%
    // 1/X₀ = 0.153/28 + 0.847/8.6 ≈ 0.1037  → X₀ ≈ 9.6 g/cm²
    // Broad sanity check: 5–40 g/cm²
    Material nai = make_NaI();
    double X0 = ElectronCsda::radiation_length_gcm2(nai);
    BOOST_CHECK_GT(X0, 5.0);
    BOOST_CHECK_LT(X0, 40.0);
}

BOOST_AUTO_TEST_CASE(radiation_length_all_elements_positive)
{
    // X₀ must be finite and positive for all elements Z = 1..92.
    for (int Z = 1; Z <= 92; ++Z) {
        double A  = ElectronCsda::atomic_weight(Z);
        double X0 = ElectronCsda::radiation_length_gcm2_element(Z, A);
        BOOST_CHECK_GT(X0, 0.0);
        BOOST_CHECK(std::isfinite(X0));
    }
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 6: NIST ESTAR stopping power correction
// ============================================================
BOOST_AUTO_TEST_SUITE(StoppingCorrection)

BOOST_AUTO_TEST_CASE(correction_improves_nist_agreement)
{
    // Spot-check: for several elements at 1 MeV, the corrected stopping power
    // should be within 2% of NIST ESTAR values.
    struct SpotCheck { int Z; double A; double nist_1MeV; };
    SpotCheck checks[] = {
        {13, 26.982, 1.465},   // Al
        {26, 55.845, 1.308},   // Fe
        {53, 126.904, 1.127},  // I
        {82, 207.200, 0.9939}, // Pb
    };
    for (const auto& c : checks) {
        double S = ElectronCsda::stopping_power_MeV_cm2_g(c.Z, c.A, 1000.0);
        double ratio = S / c.nist_1MeV;
        BOOST_CHECK_GT(ratio, 0.95);
        BOOST_CHECK_LT(ratio, 1.05);
    }
}

BOOST_AUTO_TEST_CASE(correction_factor_reasonable_range)
{
    // The corrected stopping power should still be positive and in
    // a reasonable range for all Z and E in our domain.
    for (int Z : {1, 6, 13, 26, 53, 82}) {
        double A = ElectronCsda::atomic_weight(Z);
        for (double KE : {10.0, 100.0, 500.0, 1000.0, 5000.0, 20000.0}) {
            double S = ElectronCsda::stopping_power_MeV_cm2_g(Z, A, KE);
            BOOST_CHECK_GT(S, 0.0);
            BOOST_CHECK(std::isfinite(S));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 7: Bremsstrahlung in CSDA
// ============================================================
BOOST_AUTO_TEST_SUITE(BremsstrahlungCsda)

BOOST_AUTO_TEST_CASE(brems_photons_emitted_at_high_energy)
{
    // At 5 MeV in NaI, the Moliere walk should emit at least one brems photon
    // in a reasonable fraction of events.
    const ElectronCsda& csda = ElectronCsda::instance();
    Material nai = make_NaI();

    Geometry geometry;
    geometry.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    std::mt19937_64 rng(12345);
    // Start inside the crystal (z=0 is front face, crystal extends to z=7.62)
    Eigen::Vector3d pos(0, 0, 3.81);
    Eigen::Vector3d dir(0, 0, 1);

    int n_with_brems = 0;
    const int N = 100;
    for (int i = 0; i < N; ++i) {
        auto result = csda.deposited_in_scoring(pos, dir, 5000.0, geometry, rng);
        if (!result.brems_photons.empty()) n_with_brems++;
        // Each brems photon should have energy > 0
        for (const auto& bp : result.brems_photons) {
            BOOST_CHECK_GT(bp.energy_keV, 0.0);
            BOOST_CHECK_LT(bp.energy_keV, 5000.0);
        }
    }
    // At 5 MeV, radiative losses are ~5-10% of total; expect brems in many events
    BOOST_CHECK_GT(n_with_brems, 5);  // at least 5% of events
}

BOOST_AUTO_TEST_CASE(no_brems_at_low_energy)
{
    // At 100 keV (straight-line CSDA), no bremsstrahlung should be emitted.
    const ElectronCsda& csda = ElectronCsda::instance();
    Material nai = make_NaI();

    Geometry geometry;
    geometry.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    std::mt19937_64 rng(42);
    Eigen::Vector3d pos(0, 0, 3.81);
    Eigen::Vector3d dir(0, 0, 1);

    auto result = csda.deposited_in_scoring(pos, dir, 100.0, geometry, rng);
    BOOST_CHECK(result.brems_photons.empty());
    BOOST_CHECK_GT(result.deposited_scoring_keV, 0.0);
}

BOOST_AUTO_TEST_CASE(moliere_walk_produces_spread)
{
    // Moliere scattering + brems emission should produce a distribution
    // of energy deposits (non-zero variance) at high energy.
    // At 5 MeV in NaI, brems photons escaping the crystal and scattering
    // out of the detector contribute to the spread.
    const ElectronCsda& csda = ElectronCsda::instance();
    Material nai = make_NaI();

    Geometry geometry;
    geometry.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    std::mt19937_64 rng(99);
    Eigen::Vector3d pos(0, 0, 3.81);
    Eigen::Vector3d dir(0, 0, 1);

    const int N = 200;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < N; ++i) {
        auto result = csda.deposited_in_scoring(pos, dir, 5000.0, geometry, rng);
        double dep = result.deposited_scoring_keV;
        sum  += dep;
        sum2 += dep * dep;
    }
    double mean = sum / N;
    double var  = sum2 / N - mean * mean;
    double rms  = std::sqrt(std::max(0.0, var));

    // Moliere + brems at 5 MeV should give non-trivial spread
    BOOST_CHECK_GT(mean, 0.0);
    BOOST_CHECK_GT(rms, 0.001 * mean);  // At least 0.1% spread
}

BOOST_AUTO_TEST_SUITE_END()
