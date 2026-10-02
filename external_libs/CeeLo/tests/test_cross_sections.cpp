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
#include "cross_sections/photon_epics_data.h"

#include <algorithm>
#include <cmath>
#include <vector>

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

BOOST_AUTO_TEST_CASE(photon_window_constants_match_tables) {
    // The published window is what InterSpec's adapter will refuse queries
    // outside of, so it must be exactly where the generated tables end. The top
    // node is a float32 log10 value (20 MeV decodes to 20.0000075 MeV), hence
    // the 1e-6 relative tolerance; 10 keV (log10 = -2) is exact.
    const double lo_MeV = kPhotonDataMinEnergy_keV * 1e-3;
    const double hi_MeV = kPhotonDataMaxEnergy_keV * 1e-3;
    const auto node_MeV = [](float log10_MeV) {
        return std::pow(10.0, static_cast<double>(log10_MeV));
    };
    double lowest_first = 1e300;
    for (int Z = 1; Z <= kMaxZ; ++Z) {
        const auto& d = g_photon_epics_data[Z - 1];
        const PhotonProcessCurve* curves[] = {
            &d.compton, &d.pair_production, &d.photoelectric, &d.k_photoelectric};
        for (const PhotonProcessCurve* c : curves) {
            BOOST_REQUIRE_GT(c->size, 1u);
            const float first =
                g_photon_energy_pool[g_photon_process_grid_index[c->data_offset]];
            const float last = g_photon_energy_pool[
                g_photon_process_grid_index[c->data_offset + c->size - 1]];
            BOOST_CHECK_SMALL(node_MeV(last) / hi_MeV - 1.0, 1e-6);
            BOOST_CHECK_GE(node_MeV(first), lo_MeV);
            lowest_first = std::min(lowest_first, node_MeV(first));
        }
        BOOST_CHECK_EQUAL(node_MeV(g_photon_energy_pool[
                              g_photon_process_grid_index[d.compton.data_offset]]), lo_MeV);
        BOOST_CHECK_EQUAL(node_MeV(g_photon_energy_pool[
                              g_photon_process_grid_index[d.photoelectric.data_offset]]), lo_MeV);
        const int group = (Z - 1) / kRayleighXsElementsPerGroup;
        const float r_first = g_rayleigh_log_energy[g_rayleigh_group_grid_offset[group]];
        const float r_last = g_rayleigh_log_energy[g_rayleigh_group_grid_offset[group + 1] - 1];
        BOOST_CHECK_EQUAL(node_MeV(r_first), lo_MeV);
        BOOST_CHECK_SMALL(node_MeV(r_last) / hi_MeV - 1.0, 1e-6);
    }
    BOOST_CHECK_EQUAL(lowest_first, lo_MeV);
}

BOOST_AUTO_TEST_CASE(photon_tables_clamp_only_above_the_window) {
    // Above the stored top node the accessors return the endpoint value; inside
    // the window (10-20 MeV is new) they interpolate.
    const auto& xs = CrossSectionData::instance();
    const double just_above = kPhotonDataMaxEnergy_keV * 1e-3 * (1.0 + 1e-6);
    for (int Z : {1, 13, 32, 53, 82, 92, 94, 98}) {
        const auto at_top = xs.all_cross_sections(Z, just_above);
        const auto beyond = xs.all_cross_sections(Z, 30.0);
        BOOST_CHECK_EQUAL(at_top.sigma_pe, beyond.sigma_pe);
        BOOST_CHECK_EQUAL(at_top.sigma_cs, beyond.sigma_cs);
        BOOST_CHECK_EQUAL(at_top.sigma_rs, beyond.sigma_rs);
        BOOST_CHECK_EQUAL(at_top.sigma_pp, beyond.sigma_pp);
        const auto at_10 = xs.all_cross_sections(Z, 10.0);
        const auto at_15 = xs.all_cross_sections(Z, 15.0);
        BOOST_CHECK_LT(at_15.sigma_cs, at_10.sigma_cs);   // Compton falls with E
        BOOST_CHECK_GT(at_15.sigma_pp, at_10.sigma_pp);   // pair rises with E
    }
}

BOOST_AUTO_TEST_CASE(shared_rayleigh_grids_cover_every_element) {
    const auto& xs = CrossSectionData::instance();
    for (int Z = 1; Z <= kMaxZ; ++Z) {
        for (double energy_MeV : {0.01, 0.1, 1.0, 10.0, 15.0}) {
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

BOOST_AUTO_TEST_CASE(actinide_atomic_weights_are_the_xraylib_conventions) {
    // Above uranium there is no standard atomic weight; xraylib 4.2.1 gives the
    // mass of a long-lived isotope. Pinned so a regeneration that changes the
    // convention (e.g. Pu 239.1 -> 244) is caught: it moves mu/rho by ~2%.
    const auto& xs = CrossSectionData::instance();
    BOOST_CHECK_EQUAL(xs.atomic_weight(93), 237.0);  // Np
    BOOST_CHECK_EQUAL(xs.atomic_weight(94), 239.1);  // Pu
    BOOST_CHECK_EQUAL(xs.atomic_weight(95), 243.0);  // Am
    BOOST_CHECK_EQUAL(xs.atomic_weight(96), 247.0);  // Cm
    BOOST_CHECK_EQUAL(xs.atomic_weight(97), 249.0);  // Bk
    BOOST_CHECK_EQUAL(xs.atomic_weight(98), 251.0);  // Cf
}

BOOST_AUTO_TEST_CASE(actinide_element_support_is_present) {
    const auto& xs = CrossSectionData::instance();
    for (int Z = 93; Z <= kMaxZ; ++Z) {
        const auto& e = xs.element(Z);
        BOOST_CHECK_EQUAL(e.Z, Z);
        BOOST_CHECK_GT(e.num_compton_shells, 0u);
        double occupancy = 0.0;
        for (int i = 0; i < e.num_compton_shells; ++i) occupancy += e.shell_occupancy[i];
        BOOST_CHECK_CLOSE(occupancy, static_cast<double>(Z), 0.01);
        // The bremsstrahlung pointer is uranium's table, not a duplicate.
        BOOST_CHECK(e.sb_chi_quantized == xs.element(kMaxElectronTableZ).sb_chi_quantized);
    }
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
//  Actinide photon data (Z 93-98) and the U-reusing electron tables
// ============================================================

BOOST_AUTO_TEST_SUITE(Actinides)

BOOST_AUTO_TEST_CASE(actinide_photon_data_present_and_positive) {
    const auto& xs = CrossSectionData::instance();
    for (int Z = 93; Z <= kMaxZ; ++Z) {
        for (double E : {0.015, 0.06, 0.2, 0.662, 1.332, 3.0, 10.0, 15.0, 20.0}) {
            const auto p = xs.all_cross_sections(Z, E);
            BOOST_CHECK(std::isfinite(p.sigma_total()));
            BOOST_CHECK_GT(p.sigma_pe, 0.0);
            BOOST_CHECK_GT(p.sigma_cs, 0.0);
            BOOST_CHECK_GT(p.sigma_rs, 0.0);
            if (E > 1.03) BOOST_CHECK_GT(p.sigma_pp, 0.0);
            else          BOOST_CHECK_EQUAL(p.sigma_pp, 0.0);
        }
        BOOST_CHECK_GT(xs.sigma_K_photoelectric(Z, 0.2), 0.0);
        BOOST_CHECK_GT(xs.scattering_function_S(Z, 1.0), 0.0);
        BOOST_CHECK_CLOSE(xs.scattering_function_S(Z, 1.0e6), static_cast<double>(Z), 1.0);
    }
}

BOOST_AUTO_TEST_CASE(actinide_edges_match_eadl_binding_energies) {
    // The EPDL photoelectric curve retains both sides of each absorption edge as
    // adjacent grid nodes; the K and L3 ones must sit at the EADL binding
    // energies used for fluorescence (same EPICS2023 evaluation).
    const auto& xs = CrossSectionData::instance();
    for (int Z = 93; Z <= kMaxZ; ++Z) {
        const auto& c = g_photon_epics_data[Z - 1].photoelectric;
        std::vector<double> edges_keV;
        for (uint16_t i = 1; i < c.size; ++i) {
            const float a = g_photon_energy_pool[g_photon_process_grid_index[c.data_offset + i - 1]];
            const float b = g_photon_energy_pool[g_photon_process_grid_index[c.data_offset + i]];
            if (static_cast<double>(b) - a < 1e-5)
                edges_keV.push_back(1e3 * std::pow(10.0, static_cast<double>(b)));
        }
        const auto* k = xs.fluorescence(Z);
        const auto* l = xs.l_fluorescence(Z);
        BOOST_REQUIRE(k != nullptr);
        BOOST_REQUIRE(l != nullptr);
        const auto nearest = [&](double target) {
            double best = 1e300;
            for (double e : edges_keV)
                if (std::abs(e - target) < std::abs(best - target)) best = e;
            return best;
        };
        BOOST_CHECK_CLOSE(nearest(k->k_edge_keV), static_cast<double>(k->k_edge_keV), 0.05);
        BOOST_CHECK_CLOSE(nearest(l->l3_edge_keV), static_cast<double>(l->l3_edge_keV), 0.05);
        // K photoelectric starts at the K edge.
        BOOST_CHECK_EQUAL(xs.sigma_K_photoelectric(Z, 0.999e-3 * k->k_edge_keV), 0.0);
        BOOST_CHECK_GT(xs.sigma_K_photoelectric(Z, 1.001e-3 * k->k_edge_keV), 0.0);
    }
}

BOOST_AUTO_TEST_CASE(attenuation_is_smooth_in_Z_through_the_actinides) {
    // Away from absorption edges the per-atom total cross section varies smoothly
    // with Z from Th to Cf: neighbouring ratios sigma(Z)/sigma(Z-1) may differ by
    // at most 1.5% (the largest change measured on generation, 2026-09-23, is
    // 0.42% at 662 keV, from EPDL's coherent term).
    const auto& xs = CrossSectionData::instance();
    for (double E : {0.06, 0.2, 0.662, 1.0, 3.0, 10.0, 15.0}) {
        std::vector<double> ratio;
        for (int Z = 91; Z <= kMaxZ; ++Z)
            ratio.push_back(xs.all_cross_sections(Z, E).sigma_total()
                            / xs.all_cross_sections(Z - 1, E).sigma_total());
        for (size_t i = 0; i < ratio.size(); ++i) {
            BOOST_CHECK_GT(ratio[i], 1.0);  // more electrons, more attenuation
            if (i > 0) BOOST_CHECK_SMALL(ratio[i] / ratio[i - 1] - 1.0, 0.015);
        }
    }
}

BOOST_AUTO_TEST_CASE(attenuation_agrees_with_nist_xcom) {
    // NIST XCOM (an evaluation independent of EPDL), photoelectric + Compton +
    // pair (XCOM's total without coherent), barns/atom, fetched 2026-09-23. On
    // generation these agreed to <= 0.27% for Z 92-98 and to <= 0.42% for every
    // element compared at 12-20 MeV; the pins guard the new elements and the new
    // 10-20 MeV range. Coherent scattering is left out on purpose: EPDL's coherent
    // cross section includes anomalous scattering and XCOM's does not, so they
    // differ by up to tens of percent just below high-Z K edges (see DESIGN.md).
    struct Pin { int Z; double E_MeV; double xcom_b; };
    const Pin pins[] = {
        {94, 0.100, 739.5}, {94, 0.662, 50.99}, {94, 1.000, 31.23}, {94, 15.0, 24.21},
        {92, 0.662, 48.00}, {92, 15.0, 23.42}, {82, 15.0, 19.46},
    };
    const auto& xs = CrossSectionData::instance();
    for (const Pin& p : pins) {
        BOOST_TEST_CONTEXT("Z=" << p.Z << " E=" << p.E_MeV << " MeV") {
            const auto c = xs.all_cross_sections(p.Z, p.E_MeV);
            BOOST_CHECK_CLOSE(c.sigma_pe + c.sigma_cs + c.sigma_pp, p.xcom_b, 0.5);
        }
    }
}

BOOST_AUTO_TEST_CASE(bremsstrahlung_tables_reuse_uranium_above_92) {
    // Electron-side tables stop at Z=92; Np..Cf read uranium's, bit for bit.
    const auto& xs = CrossSectionData::instance();
    for (int Z = kMaxElectronTableZ + 1; Z <= kMaxZ; ++Z) {
        BOOST_CHECK_EQUAL(electron_table_z(Z), kMaxElectronTableZ);
        for (double T : {15.0, 250.0, 1000.0, 7500.0, 19000.0})
            for (double kappa : {0.001, 0.1, 0.5, 0.9, 1.0}) {
                BOOST_CHECK_EQUAL(xs.sb_chi(Z, T, kappa), xs.sb_chi(kMaxElectronTableZ, T, kappa));
                const auto eb = xs.sb_chi_energy_bracket(T);
                const auto kb = xs.sb_chi_kappa_bracket(kappa);
                BOOST_CHECK_EQUAL(xs.sb_chi_bracketed(Z, eb, kb),
                                  xs.sb_chi_bracketed(kMaxElectronTableZ, eb, kb));
            }
    }
    for (int Z = 1; Z <= kMaxElectronTableZ; ++Z) BOOST_CHECK_EQUAL(electron_table_z(Z), Z);
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

    // Am-241 decays to Np-237 (Z=93). Relaxation covers decay daughters through
    // Z=99, one beyond the Z<=98 photon tables.
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
