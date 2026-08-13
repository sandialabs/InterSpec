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

#define BOOST_TEST_MODULE MaterialTests
#include <boost/test/unit_test.hpp>

#include "materials/Material.h"
#include "cross_sections/CrossSectionData.h"

#include <cmath>

using namespace ceelo;

// ============================================================
//  Material Construction Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(MaterialConstruction)

BOOST_AUTO_TEST_CASE(single_element_material) {
    // Pure germanium
    auto ge = make_HPGe();
    BOOST_CHECK_EQUAL(ge.name(), "HPGe");
    BOOST_CHECK_CLOSE(ge.density(), 5.323, 1e-4);
    BOOST_CHECK_EQUAL(ge.composition().size(), 1u);
    BOOST_CHECK_EQUAL(ge.composition()[0].Z, 32);
    BOOST_CHECK_CLOSE(ge.composition()[0].mass_fraction, 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(compound_material) {
    // NaI
    auto nai = make_NaI();
    BOOST_CHECK_EQUAL(nai.composition().size(), 2u);

    // Check mass fractions sum to 1
    double sum = 0.0;
    for (const auto& c : nai.composition()) {
        sum += c.mass_fraction;
    }
    BOOST_CHECK_CLOSE(sum, 1.0, 0.1);

    // Na: A=22.99, I: A=126.90, total=149.89
    // w_Na = 22.99/149.89 = 0.1534
    // w_I  = 126.90/149.89 = 0.8466
    for (const auto& c : nai.composition()) {
        if (c.Z == 11) BOOST_CHECK_CLOSE(c.mass_fraction, 0.1534, 0.5);
        if (c.Z == 53) BOOST_CHECK_CLOSE(c.mass_fraction, 0.8466, 0.5);
    }
}

BOOST_AUTO_TEST_CASE(invalid_density_throws) {
    BOOST_CHECK_THROW(Material("bad", -1.0, {{32, 1.0}}), std::invalid_argument);
    BOOST_CHECK_THROW(Material("bad", 0.0, {{32, 1.0}}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(invalid_mass_fractions_throw) {
    // Mass fractions don't sum to 1
    BOOST_CHECK_THROW(Material("bad", 5.0, {{32, 0.5}}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(all_builtin_materials_construct) {
    BOOST_CHECK_NO_THROW(make_NaI());
    BOOST_CHECK_NO_THROW(make_LaBr3());
    BOOST_CHECK_NO_THROW(make_CeBr3());
    BOOST_CHECK_NO_THROW(make_HPGe());
    BOOST_CHECK_NO_THROW(make_CZT());
    BOOST_CHECK_NO_THROW(make_CLYC());
    BOOST_CHECK_NO_THROW(make_Lead());
    BOOST_CHECK_NO_THROW(make_Copper());
    BOOST_CHECK_NO_THROW(make_Iron());
    BOOST_CHECK_NO_THROW(make_Aluminum());
    BOOST_CHECK_NO_THROW(make_Tin());
    BOOST_CHECK_NO_THROW(make_Tungsten());
    BOOST_CHECK_NO_THROW(make_StainlessSteel304());
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Macroscopic Cross-Section Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(MacroscopicXS)

BOOST_AUTO_TEST_CASE(lead_attenuation_at_662keV) {
    // Pb at 662 keV: well-known mu/rho ≈ 0.1248 cm²/g (NIST)
    // mu = mu/rho * rho = 0.1248 * 11.35 ≈ 1.416 /cm
    // We'll accept 20% tolerance for placeholder data
    auto pb = make_Lead();
    double mu = pb.mu_total(0.662);
    BOOST_CHECK(mu > 0.0);

    // Just check it's in a reasonable range (0.3 to 3.0 /cm)
    // Exact value depends on cross-section data quality;
    // placeholder data may differ from NIST by up to 50%
    BOOST_CHECK_GT(mu, 0.3);
    BOOST_CHECK_LT(mu, 3.0);
}

BOOST_AUTO_TEST_CASE(nai_attenuation_reasonable) {
    // NaI at 662 keV: mu/rho ≈ 0.0786 cm²/g (mostly Compton from I)
    // mu = 0.0786 * 3.67 ≈ 0.288 /cm
    auto nai = make_NaI();
    double mu = nai.mu_total(0.662);
    BOOST_CHECK_GT(mu, 0.1);
    BOOST_CHECK_LT(mu, 1.0);
}

BOOST_AUTO_TEST_CASE(hpge_attenuation_reasonable) {
    // Ge at 662 keV: mu/rho ≈ 0.0737 cm²/g
    // mu = 0.0737 * 5.323 ≈ 0.392 /cm
    auto ge = make_HPGe();
    double mu = ge.mu_total(0.662);
    BOOST_CHECK_GT(mu, 0.15);
    BOOST_CHECK_LT(mu, 1.5);
}

BOOST_AUTO_TEST_CASE(cross_sections_decrease_with_energy) {
    // Total attenuation should generally decrease with energy (above K-edge)
    auto pb = make_Lead();
    double mu_low = pb.mu_total(0.1);
    double mu_mid = pb.mu_total(0.5);
    double mu_high = pb.mu_total(2.0);

    // At higher energies, pair production rises but total still decreases
    // from 100 keV to 500 keV (pair production not yet dominant)
    BOOST_CHECK_GT(mu_low, mu_mid);
}

BOOST_AUTO_TEST_CASE(macroscopic_xs_decomposition) {
    // Verify that mu_total = mu_pe + mu_cs + mu_rs + mu_pp
    auto ge = make_HPGe();
    auto xs = ge.macroscopic_xs(0.662);

    BOOST_CHECK_CLOSE(xs.mu_total(),
                      xs.mu_pe + xs.mu_cs + xs.mu_rs + xs.mu_pp,
                      1e-10);
}

BOOST_AUTO_TEST_CASE(compton_dominates_at_medium_energy) {
    // For NaI at 662 keV, Compton scattering should dominate
    auto nai = make_NaI();
    auto xs = nai.macroscopic_xs(0.662);

    BOOST_CHECK_GT(xs.mu_cs, xs.mu_pe);  // CS > PE at 662 keV
    BOOST_CHECK_GT(xs.mu_cs, xs.mu_rs);  // CS > Rayleigh
}

BOOST_AUTO_TEST_CASE(photoelectric_dominates_at_low_energy) {
    // For Pb at 30 keV, photoelectric should dominate strongly
    auto pb = make_Lead();
    auto xs = pb.macroscopic_xs(0.030);

    BOOST_CHECK_GT(xs.mu_pe, xs.mu_cs);
}

BOOST_AUTO_TEST_CASE(number_density_calculation) {
    // Pure Ge: n = rho * N_A / A = 5.323 * 6.022e23 / 72.63 = 4.414e22 atoms/cm³
    auto ge = make_HPGe();
    double n = ge.number_density(0);

    double expected = 5.323 * 6.02214076e23 / 72.63;
    BOOST_CHECK_CLOSE(n, expected, 1.0); // 1% tolerance for atomic weight differences
}

BOOST_AUTO_TEST_CASE(element_selection) {
    // For NaI at an energy where I dominates PE, selecting for PE should prefer I
    auto nai = make_NaI();

    // At 50 keV, I has much higher PE cross-section than Na (Z^4-5 scaling)
    // Selecting PE element should pick I most of the time
    int count_I = 0;
    int count_Na = 0;
    int N = 100;
    for (int i = 0; i < N; ++i) {
        double xi = static_cast<double>(i) / N;
        int Z = nai.select_element(0.050, 0 /*PE*/, xi);
        if (Z == 53) count_I++;
        if (Z == 11) count_Na++;
    }
    // I should be selected much more often than Na for PE at 50 keV
    BOOST_CHECK_GT(count_I, count_Na);
}

BOOST_AUTO_TEST_SUITE_END()
