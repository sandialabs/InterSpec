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
/// @brief Escaping electrons and positrons, source and Compton K fluorescence.
///
/// Each test checks a physical bound or invariant rather than a frozen number;
/// the GEANT4 comparisons live in the full-spectrum study.

#define BOOST_TEST_MODULE SpectrumFidelityTests
#include <boost/test/unit_test.hpp>

#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"
#include "physics/ElectronCsda.h"
#include "transport/PhotonTransport.h"

#include <Eigen/Core>
#include <cmath>
#include <random>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kMeC2 = 510.998950;

Eigen::Vector3d random_unit(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double ct = 2.0 * u(rng) - 1.0, st = std::sqrt(1.0 - ct * ct), ph = 6.283185307 * u(rng);
    return {st * std::cos(ph), st * std::sin(ph), ct};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Positrons)

BOOST_AUTO_TEST_CASE(the_walk_flags_a_particle_leaving_the_geometry)
{
    const Material ge = make_HPGe();
    Geometry g;
    g.set_detector(&ge, BoxDims{2.0, 2.0, 2.0});
    std::mt19937_64 rng(11);
    TransportConfig c;
    c.disable_brems = true;
    for (bool positron : {false, true}) {
        // 1 µm inside the front face, heading straight out at 500 keV.
        const auto out = walk_electron({0.0, 0.0, 1e-4}, {0.0, 0.0, -1.0}, 500.0, g, c, rng,
                                       positron);
        BOOST_CHECK(out.escaped);
        BOOST_CHECK_GT(out.exit_KE_keV, 400.0);
        if (!out.annihilated_in_flight)
            BOOST_CHECK_CLOSE(out.deposited_scoring_keV + out.exit_KE_keV, 500.0, 1e-6);
        // Deep inside: stops, nothing escapes.
        const auto in = walk_electron({0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, 500.0, g, c, rng,
                                      positron);
        BOOST_CHECK(!in.escaped);
    }
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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Fluorescence)

BOOST_AUTO_TEST_CASE(source_fluorescence_is_on_by_default)
{
    SourceGeometry s;
    BOOST_CHECK(s.source_fluorescence_enabled());
}

BOOST_AUTO_TEST_CASE(a_lead_shell_emits_its_k_xrays)
{
    // At 200 keV nothing else leaves 2 mm of Pb at 72-89 keV: a Compton photon
    // keeps at least 112 keV and the photoelectron is below the brems floor.
    const Material pb = make_Lead();
    std::mt19937_64 rng(17);
    for (bool on : {true, false}) {
        SourceGeometry s;
        s.configure_point(Eigen::Vector3d(0.0, 0.0, -10.0));
        s.add_shield(&pb, 0.2);
        s.set_source_fluorescence(on);
        int k_lines = 0;
        for (int i = 0; i < 20000; ++i) {
            const auto r = s.transport_source_photon(
                Eigen::Vector3d(0.0, 0.0, -10.0), random_unit(rng), 200.0, rng);
            for (const auto& sec : r.secondaries)
                if (sec.energy_keV > 72.0 && sec.energy_keV < 89.0) ++k_lines;
        }
        if (on) BOOST_CHECK_GT(k_lines, 200);
        else    BOOST_CHECK_EQUAL(k_lines, 0);
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
