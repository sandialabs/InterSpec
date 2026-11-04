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

/// @file test_spectrum_fidelity.cpp
/// @brief Escaping positrons, source and Compton K fluorescence, the geometry
///        safety distance and the Detailed electron walk.
///
/// Each test checks a physical bound or invariant rather than a frozen number;
/// the GEANT4 comparisons live in the full-spectrum study.

#define BOOST_TEST_MODULE SpectrumFidelityTests
#include <boost/test/unit_test.hpp>

#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"
#include "physics/ElectronCsda.h"
#include "cross_sections/CrossSectionData.h"
#include "transport/PhotonTransport.h"
#include "transport/TransportUtils.h"

#include <Eigen/Core>
#include <cmath>
#include <random>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kMeC2 = 510.998950;

TransportConfig config_with(ElectronTransport walk) {
    TransportConfig c;
    c.electron_transport = walk;
    return c;
}

Eigen::Vector3d random_unit(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double ct = 2.0 * u(rng) - 1.0, st = std::sqrt(1.0 - ct * ct), ph = 6.283185307 * u(rng);
    return {st * std::cos(ph), st * std::sin(ph), ct};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Positrons)

BOOST_AUTO_TEST_CASE(both_walks_flag_an_electron_leaving_the_geometry)
{
    const Material ge = make_HPGe();
    Geometry g;
    g.set_detector(&ge, BoxDims{2.0, 2.0, 2.0});
    std::mt19937_64 rng(11);
    for (auto walk : {ElectronTransport::Fast, ElectronTransport::Detailed}) {
        TransportConfig c = config_with(walk);
        c.disable_brems = true;
        // 1 µm inside the front face, heading straight out at 500 keV.
        const auto out = walk_electron({0.0, 0.0, 1e-4}, {0.0, 0.0, -1.0}, 500.0, g, c, rng);
        BOOST_CHECK(out.escaped);
        BOOST_CHECK_GT(out.exit_KE_keV, 400.0);
        BOOST_CHECK_CLOSE(out.deposited_scoring_keV + out.exit_KE_keV + out.escaped_delta_keV,
                          500.0, 1e-6);
        // Deep inside: stops, nothing escapes.
        const auto in = walk_electron({0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, 500.0, g, c, rng);
        BOOST_CHECK(!in.escaped);
        BOOST_CHECK_CLOSE(in.deposited_scoring_keV, 500.0, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(an_escaped_positron_emits_no_annihilation_photons)
{
    // A bare 1 cm CZT cube at 2614.5 keV: positrons from pairs near its faces
    // leave. Such a history can deposit at most E - 2 m_e c^2.
    const Material czt = make_CZT();
    Geometry g;
    g.set_detector(&czt, BoxDims{0.5, 0.5, 1.0});
    TransportConfig c;
    c.record_interactions = true;
    std::mt19937_64 rng(5);
    const double E = 2614.511;
    int escaped = 0;
    for (int i = 0; i < 200000; ++i) {
        const auto r = transport_photon({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, E, g, c, rng);
        if (!r.pair_history.occurred || !r.pair_history.positron_escaped) continue;
        ++escaped;
        BOOST_CHECK_EQUAL(r.pair_history.annihilation[0].emitted_energy_keV, 0.0);
        BOOST_CHECK_EQUAL(r.pair_history.annihilation[1].emitted_energy_keV, 0.0);
        BOOST_CHECK_LE(r.energy_deposited_scoring, E - 2.0 * kMeC2 + 1e-6);
    }
    BOOST_CHECK_GT(escaped, 50);
}

BOOST_AUTO_TEST_CASE(annihilation_pair_kinematics)
{
    std::mt19937_64 rng(3);
    const auto rest = sample_annihilation_pair(false, 0.0, rng);
    BOOST_CHECK_CLOSE(rest.energy_keV[0], kMeC2, 1e-9);
    BOOST_CHECK_CLOSE(rest.energy_keV[1], kMeC2, 1e-9);
    BOOST_CHECK_CLOSE(rest.direction[0].dot(rest.direction[1]), -1.0, 1e-9);
    const auto flight = sample_annihilation_pair(true, 300.0, rng);
    BOOST_CHECK_CLOSE(flight.energy_keV[0] + flight.energy_keV[1], 2.0 * kMeC2 + 300.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(source_positron_transport_is_on_by_default)
{
    SourceGeometry s;
    BOOST_CHECK(s.source_positron_transport_enabled());
    BOOST_CHECK(s.source_fluorescence_enabled());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Fluorescence)

BOOST_AUTO_TEST_CASE(a_lead_shell_emits_its_k_xrays)
{
    const Material pb = make_Lead();
    std::mt19937_64 rng(17);
    for (bool on : {true, false}) {
        SourceGeometry s;
        s.configure_point(Eigen::Vector3d(0.0, 0.0, -10.0));
        s.add_shield(&pb, 0.2);
        s.set_source_fluorescence(on);
        int k_lines = 0, fluorescence = 0;
        for (int i = 0; i < 20000; ++i) {
            const auto r = s.transport_source_photon(
                Eigen::Vector3d(0.0, 0.0, -10.0), random_unit(rng), 200.0, rng);
            for (const auto& sec : r.secondaries) {
                if (sec.origin != SourceGeometry::SourceSecondaryPhoton::Origin::Fluorescence)
                    continue;
                ++fluorescence;
                if (sec.energy_keV > 72.0 && sec.energy_keV < 89.0) ++k_lines;
            }
        }
        if (on) {
            BOOST_CHECK_GT(k_lines, 200);
            BOOST_CHECK_GT(k_lines, 0.8 * fluorescence);  // unscattered Pb K lines dominate
        } else {
            BOOST_CHECK_EQUAL(fluorescence, 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(k_lines_carry_their_filling_subshell)
{
    // Line energies are EADL binding-energy differences: Pb K-alpha1 (filled
    // from L3) at the measured 74.97 keV, not EADL's tabulated 75.25 keV.
    const FluorescenceData* pb = CrossSectionData::instance().fluorescence(82);
    BOOST_REQUIRE(pb && pb->num_lines > 0);
    BOOST_CHECK_EQUAL(pb->line_filling_subshell[0], kEadlL3);
    BOOST_CHECK_CLOSE(pb->line_energy_keV[0], 74.970, 0.01);
    for (int i = 0; i < pb->num_lines; ++i)
        BOOST_CHECK_LT(pb->line_energy_keV[i], pb->k_edge_keV);
    const KLineFamilyMeans means = k_line_family_means(*pb);
    BOOST_CHECK(means.alpha_keV > 72.0 && means.alpha_keV < 75.5);
    BOOST_CHECK(means.beta_keV > 84.0 && means.beta_keV < 88.0);
    BOOST_CHECK_CLOSE(means.alpha_probability + means.beta_probability, 1.0, 1e-3);
}

BOOST_AUTO_TEST_CASE(l_relaxation_conserves_the_subshell_energy)
{
    const LFluorescenceData* pb = CrossSectionData::instance().l_fluorescence(82);
    BOOST_REQUIRE(has_l_lines(pb));
    BOOST_CHECK(!has_l_lines(CrossSectionData::instance().l_fluorescence(32)));  // Ge
    std::mt19937_64 rng(5);
    for (int s = 0; s < 3; ++s) {
        int radiated = 0;
        for (int i = 0; i < 20000; ++i) {
            const LRelaxation r = sample_l_relaxation(*pb, s, rng);
            BOOST_REQUIRE_CLOSE(r.line_keV + r.local_keV, pb->sub[s].edge_keV, 1e-4);
            if (r.line_keV > 0.0) {
                ++radiated;
                BOOST_REQUIRE(r.line_keV >= 10.0 && r.line_keV < pb->sub[0].edge_keV);
            }
        }
        // Pb L-shell fluorescence yields are 0.1-0.4 (plus Coster-Kronig feed).
        BOOST_CHECK_GT(radiated, 1000);
    }
}

BOOST_AUTO_TEST_CASE(l_photoelectric_shares_respect_the_l_edges)
{
    const LFluorescenceData* pb = CrossSectionData::instance().l_fluorescence(82);
    BOOST_REQUIRE(pb && pb->num_pe_share_nodes > 0);
    std::mt19937_64 rng(9);
    // 14 keV lies between Pb's L3 (13.04) and L2 (15.21) edges, 15.5 keV
    // between L2 and L1 (15.87); 50 keV is below the K edge.
    int counts[3][4] = {};
    const double energies[3] = {14.0, 15.5, 50.0};
    for (int e = 0; e < 3; ++e)
        for (int i = 0; i < 20000; ++i)
            ++counts[e][sample_l_photoelectric_subshell(*pb, energies[e], 0.0, rng) + 1];
    BOOST_CHECK_EQUAL(counts[0][1] + counts[0][2], 0);  // no L1/L2 below their edges
    BOOST_CHECK_GT(counts[0][3], 10000);                // mostly L3
    BOOST_CHECK_EQUAL(counts[1][1], 0);                 // no L1 below its edge
    BOOST_CHECK_GT(counts[1][2], 2000);
    for (int sub = 1; sub <= 3; ++sub)
        BOOST_CHECK_GT(counts[2][sub], 1000);           // all three above every L edge
}

BOOST_AUTO_TEST_CASE(a_thin_lead_foil_emits_l_xrays_below_its_k_edge)
{
    const Material pb = make_Lead();
    std::mt19937_64 rng(31);
    for (bool on : {true, false}) {
        SourceGeometry s;
        s.configure_point(Eigen::Vector3d(0.0, 0.0, -10.0));
        s.add_shield(&pb, 0.002);
        s.set_source_fluorescence(on);
        int l_lines = 0, fluorescence = 0;
        for (int i = 0; i < 20000; ++i) {
            const auto r = s.transport_source_photon(
                Eigen::Vector3d(0.0, 0.0, -10.0), random_unit(rng), 50.0, rng);
            for (const auto& sec : r.secondaries) {
                if (sec.origin != SourceGeometry::SourceSecondaryPhoton::Origin::Fluorescence)
                    continue;
                ++fluorescence;
                if (sec.energy_keV > 9.0 && sec.energy_keV < 16.0) ++l_lines;
            }
        }
        if (on) BOOST_CHECK_GT(l_lines, 100);
        else    BOOST_CHECK_EQUAL(fluorescence, 0);
    }
}

BOOST_AUTO_TEST_CASE(compton_k_vacancy_relaxes_only_when_enabled)
{
    const Material nai = make_NaI();
    Geometry g;
    g.set_detector(&nai, CylinderDims{3.81, 7.62});
    std::mt19937_64 rng(23);
    for (bool on : {true, false}) {
        TransportConfig c;
        c.record_interactions = true;
        c.enable_compton_fluorescence = on;
        int from_compton = 0, iodine = 0;
        for (int i = 0; i < 20000; ++i) {
            const auto r = transport_photon({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 661.657, g, c, rng);
            BOOST_REQUIRE_LE(r.energy_deposited_scoring, 661.657 + 1e-6);
            for (const auto& f : r.fluorescence_records)
                if (f.from_compton) {
                    ++from_compton;
                    BOOST_CHECK(f.element_Z == 53 || f.element_Z == 11);
                    if (f.element_Z == 53) {
                        ++iodine;
                        BOOST_CHECK(f.fluorescence_energy_keV > 28.0 && f.fluorescence_energy_keV < 33.2);
                    }
                }
        }
        if (on) BOOST_CHECK_GT(iodine, 50);
        else    BOOST_CHECK_EQUAL(from_compton, 0);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DetailedWalk)

BOOST_AUTO_TEST_CASE(safety_distance_is_a_lower_bound)
{
    // A bulletized, round-bored HPGe with dead layers and two housing layers.
    const Material ge = make_HPGe(), al = make_Aluminum();
    Geometry g;
    g.set_detector(&ge, CylinderDims{3.095, 5.755});
    g.set_bore_hole(0.45, 4.555, /*rounded_tip=*/true);
    g.set_bullet_radius(0.5);
    g.set_dead_layer(0.0739, 0.07, 0.01);
    g.add_attenuator(&al, 0.1, 0.15, 0.0, 5.755);
    g.add_attenuator(&al, 0.1, 0.1, -0.5, 6.0);
    std::mt19937_64 rng(29);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<PathSegment> segs;
    int checked = 0;
    for (int i = 0; i < 4000; ++i) {
        const Eigen::Vector3d p(-3.5 + 7.0 * u(rng), -3.5 + 7.0 * u(rng), -0.4 + 6.5 * u(rng));
        const double safety = g.safety_distance(p);
        BOOST_REQUIRE_GE(safety, 0.0);
        for (int k = 0; k < 20; ++k) {
            g.trace_ray(p, random_unit(rng), segs);
            // Distance along this ray to the first boundary of any region.
            double d = 1e30;
            for (const auto& s : segs) {
                if (s.t_start > 1e-9) d = std::min(d, s.t_start);
                if (s.t_end > 1e-9) d = std::min(d, s.t_end);
            }
            if (d < 1e30) {
                BOOST_CHECK_LE(safety, d + 1e-9);
                ++checked;
            }
        }
    }
    BOOST_CHECK_GT(checked, 10000);
}

BOOST_AUTO_TEST_CASE(contained_and_surface_electrons_conserve_energy)
{
    const Material ge = make_HPGe();
    Geometry g;
    g.set_detector(&ge, BoxDims{5.0, 5.0, 4.0});
    TransportConfig c = config_with(ElectronTransport::Detailed);
    c.disable_brems = true;
    std::mt19937_64 rng(31);
    int back = 0;
    const int n = 20000;
    for (int i = 0; i < n; ++i) {
        const auto r = walk_electron({0.0, 0.0, 1e-7}, {0.0, 0.0, 1.0}, 50.0, g, c, rng);
        BOOST_REQUIRE_CLOSE(r.deposited_scoring_keV + (r.escaped ? r.exit_KE_keV : 0.0)
                                + r.escaped_delta_keV, 50.0, 1e-6);
        if (r.escaped && r.exit_direction.z() < 0.0) ++back;
    }
    // Normal-incidence backscatter of 50 keV electrons from germanium is ~0.3
    // (Tabata; GEANT4 option4 ~0.33). The Fast walk gives zero below 200 keV.
    const double eta = double(back) / n;
    BOOST_CHECK_GT(eta, 0.24);
    BOOST_CHECK_LT(eta, 0.40);
    const auto deep = walk_electron({0.0, 0.0, 2.0}, {1.0, 0.0, 0.0}, 300.0, g, c, rng);
    BOOST_CHECK(!deep.escaped);
    BOOST_CHECK_CLOSE(deep.deposited_scoring_keV, 300.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(impact_ionization_conserves_energy_and_matches_the_csda_yield)
{
    const Material ge = make_HPGe();
    Geometry g;
    g.set_detector(&ge, BoxDims{5.0, 5.0, 4.0});
    const auto& csda = ElectronCsda::instance();
    BOOST_CHECK_EQUAL(csda.impact_xray_yield(ge, 10.0), 0.0);  // below Ge's K edge
    // Just below a threshold that falls between the yield table's log-grid
    // nodes (Yb, 61.335 keV) the yield is exactly zero, and just above it is not.
    const Material ytterbium("Yb", 6.9, {{70, 1.0}});
    BOOST_CHECK_EQUAL(csda.impact_xray_yield(ytterbium, 61.0), 0.0);
    BOOST_CHECK_GT(csda.impact_xray_yield(ytterbium, 64.0), 0.0);
    const double y300 = csda.impact_xray_yield(ge, 300.0);
    BOOST_CHECK(y300 > 0.01 && y300 < 0.3);
    BOOST_CHECK_GT(csda.impact_xray_rate_per_gcm2(ge, 40.0),
                   csda.impact_xray_rate_per_gcm2(ge, 12.0));

    std::mt19937_64 rng(37);
    for (bool on : {true, false}) {
        TransportConfig c = config_with(ElectronTransport::Detailed);
        c.disable_brems = true;
        c.enable_impact_ionization = on;
        const int n = 20000;
        int surface_xrays = 0, deep_xrays = 0;
        for (int i = 0; i < n; ++i) {
            for (bool deep : {false, true}) {
                const auto r = deep
                    ? walk_electron({0.0, 0.0, 2.0}, {1.0, 0.0, 0.0}, 300.0, g, c, rng)
                    : walk_electron({0.0, 0.0, 1e-7}, {0.0, 0.0, 1.0}, 300.0, g, c, rng);
                double xray_keV = 0.0;
                for (const auto& bp : r.brems_photons) {
                    BOOST_REQUIRE_EQUAL(bp.impact_Z, 32);
                    BOOST_REQUIRE(bp.energy_keV > 9.8 && bp.energy_keV < 11.11);
                    xray_keV += bp.energy_keV;
                    ++(deep ? deep_xrays : surface_xrays);
                }
                BOOST_REQUIRE_CLOSE(r.deposited_scoring_keV + (r.escaped ? r.exit_KE_keV : 0.0)
                                        + r.escaped_delta_keV + xray_keV - r.unpaid_impact_keV,
                                    300.0, 1e-6);
            }
        }
        if (!on) {
            BOOST_CHECK_EQUAL(surface_xrays + deep_xrays, 0);
            continue;
        }
        // A contained electron emits Poisson(yield) x-rays in place; the
        // walked surface electrons, flight by flight, fall short only by the
        // share of the slowing-down that backscattered electrons take outside.
        const double deep_mean = double(deep_xrays) / n;
        BOOST_CHECK_LT(std::abs(deep_mean - y300), 4.0 * std::sqrt(y300 / n));
        const double surface_mean = double(surface_xrays) / n;
        BOOST_CHECK(surface_mean > 0.6 * y300 && surface_mean < y300 + 4.0 * std::sqrt(y300 / n));
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SourceImpactIonization)

BOOST_AUTO_TEST_CASE(source_electrons_emit_impact_xrays_only_when_enabled)
{
    const Material pb = make_Lead();
    std::mt19937_64 rng(41);
    for (bool on : {true, false}) {
        SourceGeometry source;
        source.configure_point(Eigen::Vector3d::Zero());
        source.add_shield(&pb, 0.05);   // thinner than a 1 MeV electron's range
        source.set_source_impact_ionization(on);
        int k_lines = 0, impact = 0;
        for (int i = 0; i < 20000; ++i) {
            const auto r = source.transport_source_electron(
                Eigen::Vector3d::Zero(), random_unit(rng), 1000.0, rng);
            for (const auto& photon : r.photons) {
                if (photon.origin != SourceGeometry::SourceSecondaryPhoton::Origin::ImpactIonization)
                    continue;
                ++impact;
                if (photon.energy_keV > 72.0 && photon.energy_keV < 89.0) ++k_lines;
            }
        }
        if (on) {
            BOOST_CHECK_GT(impact, 20);
            BOOST_CHECK_GT(k_lines, 0.8 * impact);
        } else {
            BOOST_CHECK_EQUAL(impact, 0);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
