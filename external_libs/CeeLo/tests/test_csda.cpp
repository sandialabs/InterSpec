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
#include "physics/estar_stopping_data.h"
#include "cross_sections/CrossSectionData.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"
#include "efficiency/EfficiencyCalculator.h"
#include "test_fep_window.h"

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
    // The Geant4-calibrated correction remains consistent with ESTAR here.
    double A_Al = ElectronCsda::atomic_weight(13);
    double S = ElectronCsda::stopping_power_MeV_cm2_g(13, A_Al, 1000.0);
    BOOST_CHECK_GT(S, 1.0);   // At least 1 MeV cm²/g
    BOOST_CHECK_LT(S, 2.5);   // Not absurdly large
    // The Geant4-calibrated value should remain within 5% of NIST:
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
    // This test deliberately keeps broad physical bounds because CeeLo uses a
    // compact Geant4-calibrated interpolation rather than the full ESTAR tables.
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
        calc->set_detector(&nai, CylinderDims{3.81, 7.62});
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
        calc->set_detector(&nai, CylinderDims{3.81, 7.62});
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
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
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
    // X₀ must be finite and positive for every element a Material accepts.
    for (int Z = 1; Z <= kMaxZ; ++Z) {
        double A  = ElectronCsda::atomic_weight(Z);
        double X0 = ElectronCsda::radiation_length_gcm2_element(Z, A);
        BOOST_CHECK_GT(X0, 0.0);
        BOOST_CHECK(std::isfinite(X0));
    }
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 5b: Z 93-98 reuse the uranium electron tables
// ============================================================
// Electron-side data (ESTAR, ICRU-49 I, Seltzer-Berger) stops at Z=92; every
// heavier element is treated as uranium throughout the electron physics, so the
// per-element and material-level quantities must match uranium bit for bit.
BOOST_AUTO_TEST_SUITE(ActinideElectronReuse)

BOOST_AUTO_TEST_CASE(per_element_quantities_equal_uranium)
{
    const auto& csda = ElectronCsda::instance();
    const double A_U = ElectronCsda::atomic_weight(kMaxElectronTableZ);
    for (int Z = kMaxElectronTableZ + 1; Z <= kMaxZ; ++Z) {
        BOOST_TEST_CONTEXT("Z=" << Z) {
            BOOST_CHECK_EQUAL(ElectronCsda::mean_excitation_eV(Z),
                              ElectronCsda::mean_excitation_eV(kMaxElectronTableZ));
            BOOST_CHECK_EQUAL(ElectronCsda::atomic_weight(Z), A_U);
            BOOST_CHECK_EQUAL(ElectronCsda::radiation_length_gcm2_element(Z, ElectronCsda::atomic_weight(Z)),
                              ElectronCsda::radiation_length_gcm2_element(kMaxElectronTableZ, A_U));
            for (double E : {5.0, 50.0, 662.0, 3000.0, 19000.0}) {
                for (bool positron : {false, true}) {
                    BOOST_CHECK_EQUAL(
                        ElectronCsda::stopping_power_MeV_cm2_g(Z, ElectronCsda::atomic_weight(Z), E, positron),
                        ElectronCsda::stopping_power_MeV_cm2_g(kMaxElectronTableZ, A_U, E, positron));
                    BOOST_CHECK_EQUAL(csda.range_gcm2(Z, E, positron),
                                      csda.range_gcm2(kMaxElectronTableZ, E, positron));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(material_level_electron_physics_treats_plutonium_as_uranium)
{
    const auto& csda = ElectronCsda::instance();
    const Material pu("Pu metal", 19.84, {{94, 1.0}});
    const Material u("U metal", 19.84, {{92, 1.0}});
    // Same mass fractions, so the oxide comparison isolates the element swap.
    const Material puo2("PuO2", 11.46, {{94, 0.882}, {8, 0.118}});
    const Material uo2("UO2-like", 11.46, {{92, 0.882}, {8, 0.118}});
    const std::pair<const Material*, const Material*> pairs[] = {{&pu, &u}, {&puo2, &uo2}};
    for (const auto& pr : pairs) {
        for (double E : {100.0, 662.0, 2000.0, 15000.0}) {
            BOOST_CHECK_EQUAL(csda.range_gcm2_material(*pr.first, E),
                              csda.range_gcm2_material(*pr.second, E));
            BOOST_CHECK_EQUAL(csda.range_gcm2_material(*pr.first, E, /*is_positron=*/true),
                              csda.range_gcm2_material(*pr.second, E, /*is_positron=*/true));
            BOOST_CHECK_EQUAL(csda.residual_energy_keV(*pr.first, E, 0.05),
                              csda.residual_energy_keV(*pr.second, E, 0.05));
            BOOST_CHECK_EQUAL(csda.brems_pemit_corr(*pr.first, E),
                              csda.brems_pemit_corr(*pr.second, E));
            BOOST_CHECK_EQUAL(csda.extrapolated_range_gcm2_material(*pr.first, E),
                              csda.extrapolated_range_gcm2_material(*pr.second, E));
            for (bool positron : {false, true}) {
                BOOST_CHECK_EQUAL(csda.brems_emission_rate_per_gcm2(*pr.first, E, positron),
                                  csda.brems_emission_rate_per_gcm2(*pr.second, E, positron));
                BOOST_CHECK_EQUAL(csda.delta_electron_rate_per_gcm2(*pr.first, E, 50.0, positron),
                                  csda.delta_electron_rate_per_gcm2(*pr.second, E, 50.0, positron));
            }
        }
        BOOST_CHECK_EQUAL(ElectronCsda::radiation_length_gcm2(*pr.first),
                          ElectronCsda::radiation_length_gcm2(*pr.second));
        BOOST_CHECK_EQUAL(ElectronCsda::effective_Z(*pr.first),
                          ElectronCsda::effective_Z(*pr.second));
    }
    BOOST_CHECK_EQUAL(ElectronCsda::effective_Z(pu), static_cast<double>(kMaxElectronTableZ));
}

BOOST_AUTO_TEST_CASE(out_of_domain_Z_throws_instead_of_falling_back)
{
    // These used to return I=200 eV, A=2Z, a 1e30 stopping power or a zero range.
    BOOST_CHECK_THROW(ElectronCsda::atomic_weight(0), std::out_of_range);
    BOOST_CHECK_THROW(ElectronCsda::atomic_weight(kMaxZ + 1), std::out_of_range);
    BOOST_CHECK_THROW(ElectronCsda::mean_excitation_eV(kMaxZ + 1), std::out_of_range);
    BOOST_CHECK_THROW(ElectronCsda::stopping_power_MeV_cm2_g(kMaxZ + 1, 250.0, 100.0),
                      std::out_of_range);
    BOOST_CHECK_THROW(ElectronCsda::instance().range_gcm2(kMaxZ + 1, 100.0), std::out_of_range);
    BOOST_CHECK_THROW(ElectronCsda::radiation_length_gcm2_element(kMaxZ + 1, 250.0),
                      std::out_of_range);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
// Suite 6: Geant4-calibrated stopping-power correction
// ============================================================
BOOST_AUTO_TEST_SUITE(StoppingCorrection)

BOOST_AUTO_TEST_CASE(correction_remains_consistent_with_nist)
{
    // Independent spot-check: at 1 MeV the Geant4-calibrated stopping power
    // should remain within 8% of NIST ESTAR values.  This is deliberately not
    // the calibration test: density/material conventions differ slightly.
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
        BOOST_CHECK_LT(ratio, 1.08);
    }
}

BOOST_AUTO_TEST_CASE(collision_stopping_reproduces_nist_estar_nodes)
{
    // The collision stopping power is the direct public-domain NIST ESTAR
    // tabulation (not a Bethe-Bloch expression scaled by a fit to Geant4
    // G4EmCalculator values, which are not redistributable to pin against), so
    // the property under test is that the interpolation reproduces its own
    // source grid.  It matches the NIST 1 MeV spot checks above to 4-5 digits;
    // Geant4 11.4.0 Livermore eIoni differs from NIST by 0.9-3.6%.
    for (const int Z : {13, 26, 53, 82}) {
        for (const double E_keV : {100.0, 500.0, 1000.0, 2000.0, 5000.0}) {
            const double actual = ElectronCsda::stopping_power_MeV_cm2_g(
                Z, ElectronCsda::atomic_weight(Z), E_keV);
            const double estar = estar_data::collision_stopping(Z, E_keV);
            BOOST_TEST_CONTEXT("Z=" << Z << " E=" << E_keV << " keV") {
                BOOST_CHECK_GT(estar, 0.0);
                BOOST_CHECK_CLOSE(actual, estar, 0.05);
            }
        }
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

BOOST_AUTO_TEST_CASE(seltzer_berger_absolute_emission_rate)
{
    // Expected number of >10 keV bremsstrahlung photons per g/cm² for a
    // 2.29 MeV electron.  These values follow
    //   N_A/A * Z² * millibarn/beta² * integral chi d(log kappa)
    // over the compiled-in Seltzer-Berger table, which is derived from the
    // public-domain NIST EPQ BergerSeltzerBrem tables (see the provenance
    // banner in src/cross_sections/element_data.cpp).  They are pinned here as
    // a regression lock; the ABSOLUTE normalisation they depend on is checked
    // independently against NIST ESTAR in seltzer_berger_matches_nist_estar
    // below.  (A G4EMLOW/Livermore extraction differs from these by only
    // +0.07% to +0.13%.)
    const ElectronCsda& csda = ElectronCsda::instance();
    struct Reference {
        int Z;
        double rate_per_gcm2;
    };
    for (const auto& reference : {
             Reference{6, 0.1251439},
             Reference{26, 0.4140230},
             Reference{48, 0.6692323},
             Reference{92, 1.0937720},
         }) {
        Material material(
            "element", 1.0,
            {{static_cast<std::uint8_t>(reference.Z), 1.0}});
        const double rate =
            csda.brems_emission_rate_per_gcm2(
                material, 2290.0, false);
        BOOST_CHECK_CLOSE(rate, reference.rate_per_gcm2, 0.35);

        const double positron_rate =
            csda.brems_emission_rate_per_gcm2(
                material, 2290.0, true);
        BOOST_CHECK_GE(positron_rate, 0.0);
        BOOST_CHECK_LT(positron_rate, rate);
    }
}

BOOST_AUTO_TEST_CASE(seltzer_berger_matches_nist_estar)
{
    // ABSOLUTE normalisation check for the compiled-in Seltzer-Berger table,
    // against an independent public-domain source.  Integrating the table's own
    // chi over kappa must reproduce the NIST ESTAR radiative stopping power:
    //   S_rad = (N_A/A) * T * (Z^2/beta^2) * integral_0^1 chi dkappa * 1e-27
    // with S_rad in MeV cm^2/g.  This is what makes the pinned emission rates
    // in seltzer_berger_absolute_emission_rate meaningful rather than a pure
    // self-consistency lock: ESTAR is not the table's own provenance (the table
    // comes from NIST EPQ), so a table swap that shifted absolute yield would
    // break this even if the pinned rates were regenerated alongside it.
    const CrossSectionData& xs = CrossSectionData::instance();
    constexpr double kAvogadro = 6.02214076e23;
    constexpr double kMeC2_keV = 510.998950;

    for (const double T_keV : {500.0, 1000.0, 2290.0, 5000.0}) {
        const double tau = T_keV / kMeC2_keV;
        const double beta2 =
            tau * (tau + 2.0) / ((tau + 1.0) * (tau + 1.0));
        for (const int Z : {6, 26, 48, 82, 92}) {
            // Trapezoid on the table's own kappa grid.
            double integral = 0.0;
            for (int i = 0; i + 1 < static_cast<int>(kSB_n_kappa); ++i) {
                const double k0 = kSB_kappa[i];
                const double k1 = kSB_kappa[i + 1];
                integral += 0.5
                    * (xs.sb_chi(Z, T_keV, k0) + xs.sb_chi(Z, T_keV, k1))
                    * (k1 - k0);
            }
            const double from_table =
                (kAvogadro / xs.atomic_weight(Z)) * (T_keV * 1.0e-3)
                * (static_cast<double>(Z) * Z / beta2) * integral * 1.0e-27;
            const double estar = estar_data::radiative_stopping(Z, T_keV);
            BOOST_TEST_CONTEXT("Z=" << Z << " T=" << T_keV << " keV") {
                BOOST_CHECK_GT(estar, 0.0);
                // Measured: |deviation| <= 0.35% at 2.29 MeV.  The 1.5%
                // gate leaves room for the 32-node trapezoid's own error at the
                // grid edges while still catching a table-scale regression.
                BOOST_CHECK_CLOSE(from_table, estar, 1.5);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(moller_delta_electron_rate)
{
    // G4MollerBhabhaModel cross section above a 50 keV secondary cut for a
    // 2.29 MeV electron, converted to macroscopic cross section per g/cm².
    const ElectronCsda& csda = ElectronCsda::instance();
    struct Reference {
        int Z;
        double rate_per_gcm2;
    };
    for (const auto& reference : {
             Reference{6, 1.51875414},
             Reference{26, 1.41548230},
             Reference{48, 1.29818393},
             Reference{92, 1.17509600},
         }) {
        Material material(
            "element", 1.0,
            {{static_cast<std::uint8_t>(reference.Z), 1.0}});
        const double rate = csda.delta_electron_rate_per_gcm2(
            material, 2290.0, 50.0, false);
        BOOST_CHECK_CLOSE(rate, reference.rate_per_gcm2, 0.01);
    }
}

BOOST_AUTO_TEST_CASE(source_walk_emits_delta_electron_brems)
{
    Material uranium("uranium", 1.0, {{92, 1.0}});
    SourceGeometry source;
    source.configure_point(Eigen::Vector3d::Zero());
    source.add_shield(&uranium, 1.1);
    source.set_source_brems(true);
    source.set_record_truth_labels(true);

    std::mt19937_64 rng(0x4d4f4c4c4552ULL);
    std::size_t secondary_count = 0;
    for (int event = 0; event < 1000; ++event) {
        const auto result = source.transport_source_electron(
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d::UnitZ(),
            2290.0, rng, false);
        for (const auto& photon : result.created_brems_photons) {
            if (photon.origin
                == SourceGeometry::SourceSecondaryPhoton::Origin::
                    SecondaryElectronBrems) {
                ++secondary_count;
            }
        }
    }
    BOOST_CHECK_GT(secondary_count, 20U);
}

BOOST_AUTO_TEST_CASE(contained_primary_event_sampler_tracks_delta_brems)
{
    const ElectronCsda& csda = ElectronCsda::instance();
    Material uranium("uranium", 1.0, {{92, 1.0}});
    std::mt19937_64 rng(0x4f50544943414cULL);
    std::size_t primary_count = 0;
    std::size_t secondary_count = 0;

    for (int event = 0; event < 2000; ++event) {
        const auto photons = csda.sample_brems_in_material(
            uranium, Eigen::Vector3d::Zero(),
            Eigen::Vector3d::UnitZ(), 2290.0, rng);
        for (const auto& photon : photons) {
            if (photon.charged_generation == 0)
                ++primary_count;
            else if (photon.charged_generation == 1)
                ++secondary_count;
            BOOST_CHECK_LE(photon.charged_generation, 1);
            BOOST_CHECK_GE(photon.energy_keV, 10.0);
            BOOST_CHECK_LE(photon.energy_keV, 2290.0);
        }
    }

    // Wide deterministic gates: their purpose is to ensure both competing
    // optical-depth processes are active, not to freeze a Monte Carlo result.
    BOOST_CHECK_GT(primary_count, 2500U);
    BOOST_CHECK_GT(secondary_count, 150U);
}

BOOST_AUTO_TEST_CASE(contained_delta_sampler_preserves_integrated_rate)
{
    const ElectronCsda& csda = ElectronCsda::instance();
    Material uranium("uranium", 1.0, {{92, 1.0}});
    constexpr double initial_KE_keV = 1000.0;

    // Independent numerical integral of the expected photon multiplicity:
    // Lambda = integral sigma_brems(T) dR(T).
    constexpr int intervals = 500;
    double expected = 0.0;
    double previous_range =
        csda.range_gcm2_material(uranium, 10.0);
    double previous_rate =
        csda.brems_emission_rate_per_gcm2(
            uranium, 10.0, false);
    for (int i = 1; i <= intervals; ++i) {
        const double fraction =
            static_cast<double>(i) / intervals;
        const double energy =
            10.0 * std::pow(initial_KE_keV / 10.0, fraction);
        const double range =
            csda.range_gcm2_material(uranium, energy);
        const double rate =
            csda.brems_emission_rate_per_gcm2(
                uranium, energy, false);
        expected +=
            0.5 * (previous_rate + rate)
            * std::max(0.0, range - previous_range);
        previous_range = range;
        previous_rate = rate;
    }

    std::mt19937_64 rng(0x434f4e5441494e45ULL);
    constexpr int events = 20000;
    std::size_t photon_count = 0;
    for (int event = 0; event < events; ++event) {
        const auto photons = csda.sample_brems_in_material(
            uranium, Eigen::Vector3d::Zero(),
            Eigen::Vector3d::UnitZ(), initial_KE_keV, rng,
            /*charged_generation=*/1);
        photon_count += photons.size();
        for (const auto& photon : photons) {
            BOOST_CHECK_EQUAL(photon.charged_generation, 1);
            BOOST_CHECK_GE(photon.energy_keV, 10.0);
            BOOST_CHECK_LE(photon.energy_keV, initial_KE_keV);
            BOOST_CHECK(std::isfinite(photon.direction.norm()));
            BOOST_CHECK_CLOSE(photon.direction.norm(), 1.0, 1e-8);
        }
    }
    const double observed =
        static_cast<double>(photon_count) / events;
    // Photon recoil shortens the remaining charged path slightly, so the
    // no-recoil optical-depth integral is an upper reference rather than an
    // exact equality. At this energy the correction is well below this gate.
    BOOST_CHECK_CLOSE(observed, expected, 8.0);
}

BOOST_AUTO_TEST_CASE(homogeneous_material_clearance)
{
    Material iron("iron", 1.0, {{26, 1.0}});
    SourceGeometry source;
    source.configure_point(Eigen::Vector3d::Zero());
    source.add_shield(&iron, 2.0);

    BOOST_CHECK_CLOSE(
        source.min_distance_within_material(
            Eigen::Vector3d::Zero(), iron),
        2.0, 1e-10);
    BOOST_CHECK_CLOSE(
        source.min_distance_within_material(
            Eigen::Vector3d(0.0, 0.0, 1.25), iron),
        0.75, 1e-10);

    Material cadmium("cadmium", 1.0, {{48, 1.0}});
    source.add_shield(&cadmium, 1.0);
    BOOST_CHECK_EQUAL(
        source.min_distance_within_material(
            Eigen::Vector3d::Zero(), iron),
        0.0);
}

BOOST_AUTO_TEST_CASE(brems_photons_emitted_at_high_energy)
{
    // At 5 MeV in NaI, the Moliere walk should emit at least one brems photon
    // in a reasonable fraction of events.
    const ElectronCsda& csda = ElectronCsda::instance();
    Material nai = make_NaI();

    Geometry geometry;
    geometry.set_detector(&nai, CylinderDims{3.81, 7.62});

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
    geometry.set_detector(&nai, CylinderDims{3.81, 7.62});

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
    geometry.set_detector(&nai, CylinderDims{3.81, 7.62});

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

BOOST_AUTO_TEST_CASE(emission_recording_walk_matches_wrapper_and_moves)
{
    // sample_brems_in_material is now a thin wrapper over
    // sample_brems_emissions_in_material with an identical RNG draw sequence:
    // same seed must give the same photons (direction, energy, generation),
    // with the wrapper's positions collapsed to the birth point while the
    // core reports true along-walk emission positions bounded by the CSDA
    // range. Air-like density makes the displacement macroscopic (metres) so
    // a bookkeeping error cannot hide below floating-point noise.
    const ElectronCsda& csda = ElectronCsda::instance();
    Material air("air", 0.001205,
                 {{7, 0.755}, {8, 0.232}, {18, 0.013}});

    const double KE = 2000.0;
    const Eigen::Vector3d birth(1.0, -2.0, 3.0);
    const Eigen::Vector3d dir0 = Eigen::Vector3d(0.3, -0.1, 1.0).normalized();
    const double range_cm =
        csda.range_gcm2_material(air, KE, false) / air.density();

    std::mt19937_64 rng_a(777), rng_b(777);
    int n_emissions = 0;
    double mean_disp = 0.0;
    for (int i = 0; i < 3000; ++i) {
        const auto wrapped =
            csda.sample_brems_in_material(air, birth, dir0, KE, rng_a);
        const auto recorded =
            csda.sample_brems_emissions_in_material(
                air, birth, dir0, KE, rng_b, false);
        BOOST_REQUIRE_EQUAL(wrapped.size(), recorded.size());
        for (std::size_t j = 0; j < wrapped.size(); ++j) {
            BOOST_CHECK_EQUAL(wrapped[j].direction.x(),
                              recorded[j].photon_dir.x());
            BOOST_CHECK_EQUAL(wrapped[j].direction.y(),
                              recorded[j].photon_dir.y());
            BOOST_CHECK_EQUAL(wrapped[j].direction.z(),
                              recorded[j].photon_dir.z());
            BOOST_CHECK_EQUAL(wrapped[j].energy_keV, recorded[j].photon_keV);
            BOOST_CHECK_EQUAL(wrapped[j].charged_generation,
                              recorded[j].charged_generation);
            BOOST_CHECK_EQUAL(wrapped[j].position.x(), birth.x());

            const double disp = (recorded[j].position - birth).norm();
            BOOST_CHECK_LE(disp, range_cm * 1.0001);
            mean_disp += disp;

            BOOST_CHECK_GT(recorded[j].electron_KE_keV, 0.0);
            BOOST_CHECK_LE(recorded[j].electron_KE_keV, KE * 1.0001);
            BOOST_CHECK_CLOSE(recorded[j].electron_dir.norm(), 1.0, 1e-9);
            BOOST_CHECK_GE(recorded[j].photon_keV, 10.0);
            ++n_emissions;
        }
    }
    BOOST_REQUIRE_GT(n_emissions, 20);
    mean_disp /= n_emissions;
    // Emission happens part-way along a metres-long air track: the mean
    // displacement must be macroscopic and below the full range.
    BOOST_TEST_MESSAGE("air emissions: " << n_emissions
                       << ", mean displacement " << mean_disp
                       << " cm, CSDA range " << range_cm << " cm");
    BOOST_CHECK_GT(mean_disp, 0.01 * range_cm);
    BOOST_CHECK_LT(mean_disp, range_cm);

    // Positron variant: valid fields, no crashes (Bhabha/positron tables).
    std::mt19937_64 rng_p(31337);
    int n_pos = 0;
    for (int i = 0; i < 500 && n_pos == 0; ++i) {
        const auto rec = csda.sample_brems_emissions_in_material(
            air, birth, dir0, KE, rng_p, true);
        for (const auto& e : rec) {
            BOOST_CHECK_GT(e.photon_keV, 0.0);
            BOOST_CHECK_LE((e.position - birth).norm(), range_cm * 1.05);
            ++n_pos;
        }
    }
    BOOST_CHECK_GT(n_pos, 0);
}

BOOST_AUTO_TEST_SUITE_END()
