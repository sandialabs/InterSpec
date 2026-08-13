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

#define BOOST_TEST_MODULE CrossSectionTests
#include <boost/test/unit_test.hpp>

#include "cross_sections/CrossSectionData.h"

#include <cmath>

using namespace ceelo;

// ============================================================
//  Basic Data Access Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(DataAccess)

BOOST_AUTO_TEST_CASE(singleton_returns_same_instance) {
    const auto& xs1 = CrossSectionData::instance();
    const auto& xs2 = CrossSectionData::instance();
    BOOST_CHECK_EQUAL(&xs1, &xs2);
}

BOOST_AUTO_TEST_CASE(element_data_available_for_key_elements) {
    const auto& xs = CrossSectionData::instance();

    // Key elements should have data
    int key_elements[] = {11, 13, 26, 29, 32, 48, 50, 52, 53, 55, 56, 57, 58, 74, 82};
    for (int Z : key_elements) {
        const auto& elem = xs.element(Z);
        BOOST_CHECK_EQUAL(elem.Z, Z);
        BOOST_REQUIRE(elem.sb_chi_quantized != nullptr);
        BOOST_CHECK_GT(elem.num_compton_shells, 0u);
        BOOST_REQUIRE(elem.shell_occupancy != nullptr);
        BOOST_REQUIRE(elem.shell_binding_keV != nullptr);
        BOOST_REQUIRE(elem.shell_J0 != nullptr);
    }
}

BOOST_AUTO_TEST_CASE(photon_tables_restore_historical_upper_clamp) {
    const auto& xs = CrossSectionData::instance();
    for (int Z : {1, 13, 32, 53, 82, 92}) {
        const auto at_limit = xs.all_cross_sections(Z, 10.0);
        const auto above_limit = xs.all_cross_sections(Z, 20.0);
        BOOST_CHECK_EQUAL(at_limit.sigma_pe, above_limit.sigma_pe);
        BOOST_CHECK_EQUAL(at_limit.sigma_cs, above_limit.sigma_cs);
        BOOST_CHECK_EQUAL(at_limit.sigma_rs, above_limit.sigma_rs);
        BOOST_CHECK_EQUAL(at_limit.sigma_pp, above_limit.sigma_pp);
    }
}

BOOST_AUTO_TEST_CASE(shared_rayleigh_grids_cover_every_element) {
    const auto& xs = CrossSectionData::instance();
    for (int Z = 1; Z <= 92; ++Z) {
        for (double energy_MeV : {0.01, 0.1, 1.0, 10.0}) {
            const double rayleigh = xs.sigma_rayleigh(Z, energy_MeV);
            BOOST_CHECK(std::isfinite(rayleigh));
            BOOST_CHECK_GT(rayleigh, 0.0);
            BOOST_CHECK_EQUAL(
                rayleigh, xs.all_cross_sections(Z, energy_MeV).sigma_rs
            );
        }
    }
}

BOOST_AUTO_TEST_CASE(atomic_weights_reasonable) {
    const auto& xs = CrossSectionData::instance();

    // Check some known atomic weights
    BOOST_CHECK_CLOSE(xs.atomic_weight(1), 1.008, 1.0);   // Hydrogen
    BOOST_CHECK_CLOSE(xs.atomic_weight(6), 12.011, 1.0);   // Carbon
    BOOST_CHECK_CLOSE(xs.atomic_weight(11), 22.990, 1.0);  // Sodium
    BOOST_CHECK_CLOSE(xs.atomic_weight(26), 55.845, 1.0);  // Iron
    BOOST_CHECK_CLOSE(xs.atomic_weight(32), 72.630, 1.0);  // Germanium
    BOOST_CHECK_CLOSE(xs.atomic_weight(53), 126.904, 1.0); // Iodine
    BOOST_CHECK_CLOSE(xs.atomic_weight(82), 207.2, 1.0);   // Lead
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Cross-Section Value Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(CrossSectionValues)

BOOST_AUTO_TEST_CASE(photoelectric_decreases_with_energy) {
    const auto& xs = CrossSectionData::instance();

    // For any element, PE should decrease with increasing energy (above edges)
    double pe_100 = xs.sigma_photoelectric(82, 0.100);  // Pb at 100 keV
    double pe_500 = xs.sigma_photoelectric(82, 0.500);  // Pb at 500 keV
    double pe_1000 = xs.sigma_photoelectric(82, 1.000); // Pb at 1 MeV

    BOOST_CHECK_GT(pe_100, pe_500);
    BOOST_CHECK_GT(pe_500, pe_1000);
}

BOOST_AUTO_TEST_CASE(compton_decreases_with_energy) {
    const auto& xs = CrossSectionData::instance();

    // Compton (KN) cross-section decreases with energy
    double cs_100 = xs.sigma_compton(32, 0.100);
    double cs_662 = xs.sigma_compton(32, 0.662);
    double cs_1000 = xs.sigma_compton(32, 1.000);

    BOOST_CHECK_GT(cs_100, cs_662);
    BOOST_CHECK_GT(cs_662, cs_1000);
}

BOOST_AUTO_TEST_CASE(pair_production_zero_below_threshold) {
    const auto& xs = CrossSectionData::instance();

    // PP should be zero below 1.022 MeV threshold
    double pp_500 = xs.sigma_pair_production(82, 0.500);
    BOOST_CHECK_LT(pp_500, 1e-10);

    double pp_1000 = xs.sigma_pair_production(82, 1.000);
    BOOST_CHECK_LT(pp_1000, 1e-10);
}

BOOST_AUTO_TEST_CASE(pair_production_increases_above_threshold) {
    const auto& xs = CrossSectionData::instance();

    // PP should increase with energy above threshold
    double pp_2 = xs.sigma_pair_production(82, 2.0);
    double pp_5 = xs.sigma_pair_production(82, 5.0);
    double pp_10 = xs.sigma_pair_production(82, 10.0);

    BOOST_CHECK_GT(pp_5, pp_2);
    BOOST_CHECK_GT(pp_10, pp_5);
}

BOOST_AUTO_TEST_CASE(cross_sections_scale_with_Z) {
    const auto& xs = CrossSectionData::instance();

    // At a given energy, PE scales roughly as Z^4-5
    // So Pb (Z=82) >> Ge (Z=32) >> Na (Z=11) for PE
    double pe_Na = xs.sigma_photoelectric(11, 0.100);
    double pe_Ge = xs.sigma_photoelectric(32, 0.100);
    double pe_Pb = xs.sigma_photoelectric(82, 0.100);

    BOOST_CHECK_GT(pe_Ge, pe_Na);
    BOOST_CHECK_GT(pe_Pb, pe_Ge);

    // Compton scales roughly as Z (number of electrons)
    double cs_Na = xs.sigma_compton(11, 0.662);
    double cs_Ge = xs.sigma_compton(32, 0.662);
    double cs_Pb = xs.sigma_compton(82, 0.662);

    BOOST_CHECK_GT(cs_Ge, cs_Na);
    BOOST_CHECK_GT(cs_Pb, cs_Ge);
}

BOOST_AUTO_TEST_CASE(all_cross_sections_returns_consistent_values) {
    const auto& xs = CrossSectionData::instance();

    // all_cross_sections should return the same values as individual queries
    double E = 0.662;
    int Z = 32;

    auto all = xs.all_cross_sections(Z, E);
    double pe = xs.sigma_photoelectric(Z, E);
    double cs = xs.sigma_compton(Z, E);
    double rs = xs.sigma_rayleigh(Z, E);
    double pp = xs.sigma_pair_production(Z, E);

    BOOST_CHECK_CLOSE(all.sigma_pe, pe, 0.01);
    BOOST_CHECK_CLOSE(all.sigma_cs, cs, 0.01);
    BOOST_CHECK_CLOSE(all.sigma_rs, rs, 0.01);
    BOOST_CHECK_SMALL(std::abs(all.sigma_pp - pp), 1e-10);
}

BOOST_AUTO_TEST_CASE(total_cross_section_is_sum_of_parts) {
    const auto& xs = CrossSectionData::instance();

    auto all = xs.all_cross_sections(82, 5.0);
    double total = all.sigma_pe + all.sigma_cs + all.sigma_rs + all.sigma_pp;
    BOOST_CHECK_CLOSE(all.sigma_total(), total, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Scattering Functions Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(ScatteringFunctions)

BOOST_AUTO_TEST_CASE(incoherent_scattering_function_limits) {
    const auto& xs = CrossSectionData::instance();

    // S(x,Z) should approach 0 at small x and Z at large x
    int Z = 32; // Germanium

    double S_small = xs.scattering_function_S(Z, 0.01);
    double S_large = xs.scattering_function_S(Z, 1.0e6);

    BOOST_CHECK_LT(S_small, static_cast<double>(Z));  // Less than Z at small x
    BOOST_CHECK_GT(S_large, 0.5 * Z);  // Approaches Z at large x
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Fluorescence Data Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(FluorescenceTests)

BOOST_AUTO_TEST_CASE(ge_has_fluorescence) {
    const auto& xs = CrossSectionData::instance();
    const auto* fl = xs.fluorescence(32);
    BOOST_REQUIRE(fl != nullptr);

    // Ge K-edge should be around 11.1 keV
    BOOST_CHECK_CLOSE(fl->k_edge_keV, 11.10, 5.0);

    // Fluorescence yield should be around 0.55
    BOOST_CHECK_GT(fl->fluorescence_yield, 0.3);
    BOOST_CHECK_LT(fl->fluorescence_yield, 0.8);

    // Should have at least 2 lines
    BOOST_CHECK_GE(fl->num_lines, 2);

    // Line energies should be below K-edge
    for (int i = 0; i < fl->num_lines; ++i) {
        BOOST_CHECK_GT(fl->line_energy_keV[i], 0.0);
        BOOST_CHECK_LT(fl->line_energy_keV[i], fl->k_edge_keV + 1.0);
    }

    // Probabilities should sum to approximately 1
    double prob_sum = 0.0;
    for (int i = 0; i < fl->num_lines; ++i) {
        prob_sum += fl->line_probability[i];
    }
    BOOST_CHECK_CLOSE(prob_sum, 1.0, 5.0);
}

BOOST_AUTO_TEST_CASE(pb_has_fluorescence) {
    const auto& xs = CrossSectionData::instance();
    const auto* fl = xs.fluorescence(82);
    BOOST_REQUIRE(fl != nullptr);

    // Pb K-edge around 88 keV
    BOOST_CHECK_CLOSE(fl->k_edge_keV, 88.0, 5.0);
    BOOST_CHECK_GT(fl->fluorescence_yield, 0.9);
}

BOOST_AUTO_TEST_CASE(np_decay_daughter_has_relaxation) {
    const auto& xs = CrossSectionData::instance();

    // Am-241 decays to Np-237 (Z=93). Relaxation therefore intentionally has a
    // wider domain than the Z<=92 transport cross-section tables.
    const auto* k = xs.fluorescence(93);
    const auto* l = xs.l_fluorescence(93);
    BOOST_REQUIRE(k != nullptr);
    BOOST_REQUIRE(l != nullptr);
    BOOST_CHECK_GT(k->fluorescence_yield, 0.9);
    BOOST_CHECK_GT(l->sub[1].fluorescence_yield, 0.4);
    BOOST_CHECK_GT(l->sub[2].fluorescence_yield, 0.4);
    BOOST_CHECK(xs.fluorescence(100) == nullptr);
    BOOST_CHECK(xs.l_fluorescence(100) == nullptr);
}

BOOST_AUTO_TEST_CASE(na_has_no_fluorescence) {
    const auto& xs = CrossSectionData::instance();
    const auto* fl = xs.fluorescence(11);
    // Na (Z=11) has very low fluorescence yield — data may or may not be present
    // but if present, yield should be very low
    if (fl != nullptr) {
        BOOST_CHECK_LT(fl->fluorescence_yield, 0.1);
    }
}

BOOST_AUTO_TEST_SUITE_END()
