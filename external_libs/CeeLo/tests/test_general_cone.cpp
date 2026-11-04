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

/// @file test_general_cone.cpp
/// @brief EfficiencyCalculator::enable_general_cone_sampling: sources beside the detector get a
///        cone with its own axis instead of the full sphere - the same answer, far cheaper - and
///        nothing else moves.

#define BOOST_TEST_MODULE GeneralConeTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "geometry/BoundingCone.h"
#include "materials/Material.h"
#include "test_fep_window.h"

#include <cmath>
#include <string>

using namespace ceelo;

namespace {

const Material& mat_NaI() { static Material m = make_NaI();      return m; }
const Material& mat_Ge()  { static Material m = make_HPGe();     return m; }
const Material& mat_Al()  { static Material m = make_Aluminum(); return m; }

void nai3x3(EfficiencyCalculator& calc) {
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&mat_NaI(), CylinderDims{3.81, 7.62});
    calc.add_attenuator(&mat_Al(), 0.05, 0.05, 0.0, 7.62);
}

/// Bored, dead-layered coax in an Al can with a 2 mm back plate (z_end past the crystal).
void hpge_bored(EfficiencyCalculator& calc) {
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&mat_Ge(), CylinderDims{2.915, 6.89});
    calc.set_bore_hole(0.5, 4.0);
    calc.set_dead_layer(0.07, 0.07, 0.07);
    calc.add_attenuator(&mat_Al(), 0.15, 0.15, -0.5, 6.89 + 0.2);
}

void nai_box(EfficiencyCalculator& calc) {
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&mat_NaI(), BoxDims{2.5, 5.0, 20.0});
}

SimulationConfig cfg(const double energy_keV, const uint64_t events, const uint64_t seed) {
    SimulationConfig c;
    c.energy_keV = energy_keV;
    c.termination.max_events = events;
    c.termination.min_events = events;
    c.num_threads = 4;
    c.batch_size = 20000;
    c.seed = seed;
    return c;
}

/// |a - b| in combined sigmas.
double pull(const double a, const double sa, const double b, const double sb) {
    return std::abs(a - b) / std::sqrt(sa * sa + sb * sb);
}

/// Runs the same point source with the general cone off (isotropic beside the detector) and on,
/// and checks both efficiencies agree; returns the history ratio for the same precision.
void check_point(void (*setup)(EfficiencyCalculator&), const Eigen::Vector3d& src,
                 const double energy_keV, const std::string& label) {
    EfficiencyCalculator iso, gen;
    setup(iso);
    setup(gen);
    iso.set_point_source(src);
    gen.set_point_source(src);
    gen.enable_general_cone_sampling(true);

    const EfficiencyResult r_iso = iso.compute(cfg(energy_keV, 4000000, 11));
    const EfficiencyResult r_gen = gen.compute(cfg(energy_keV, 400000, 12));

    const double p_fep = pull(r_iso.full_energy_peak_efficiency, r_iso.fep_uncertainty,
                              r_gen.full_energy_peak_efficiency, r_gen.fep_uncertainty);
    const double p_tot = pull(r_iso.total_efficiency, r_iso.total_uncertainty,
                              r_gen.total_efficiency, r_gen.total_uncertainty);
    BOOST_TEST_MESSAGE(label << " @" << energy_keV << " keV: FEP iso "
                       << r_iso.full_energy_peak_efficiency << " +- " << r_iso.fep_uncertainty
                       << " vs coned " << r_gen.full_energy_peak_efficiency << " +- "
                       << r_gen.fep_uncertainty << " (" << p_fep << " sigma); total "
                       << r_iso.total_efficiency << " vs " << r_gen.total_efficiency
                       << " (" << p_tot << " sigma)");
    BOOST_CHECK_MESSAGE(p_fep < 4.0, label << ": FEP differs by " << p_fep << " sigma");
    BOOST_CHECK_MESSAGE(p_tot < 4.0, label << ": total differs by " << p_tot << " sigma");

    // The point of it: 10x fewer histories and still a tighter estimate.
    BOOST_CHECK_MESSAGE(r_gen.fep_uncertainty / r_gen.full_energy_peak_efficiency
                            < r_iso.fep_uncertainty / r_iso.full_energy_peak_efficiency,
                        label << ": the coned run is not more precise");
}

}  // namespace


BOOST_AUTO_TEST_CASE(envelope_cone_contains_the_envelope) {
    // Sources beside and behind a cylinder: every envelope rim point is inside the cone.
    const double R = 4.0, z0 = -0.5, z1 = 8.0;
    for (const Eigen::Vector3d src : {Eigen::Vector3d(20.0, 0.0, 3.0), Eigen::Vector3d(10.0, 5.0, 12.0),
                                      Eigen::Vector3d(0.0, 0.0, 30.0), Eigen::Vector3d(6.0, 0.0, 8.0)}) {
        Eigen::Vector3d axis;
        double cos_alpha = -1.0;
        BOOST_REQUIRE(envelope_bounding_cone(true, R, 0.0, 0.0, z0, z1, src, axis, cos_alpha));
        BOOST_CHECK(cos_alpha > 0.0);
        for (int i = 0; i < 360; ++i) {
            const double phi = 2.0 * 3.14159265358979323846 * i / 360.0;
            for (const double z : {z0, z1}) {
                const Eigen::Vector3d p(R * std::cos(phi), R * std::sin(phi), z);
                BOOST_CHECK_GE(axis.dot((p - src).normalized()), cos_alpha - 1e-12);
            }
        }
    }

    // Inside the envelope: no cone.
    Eigen::Vector3d axis;
    double cos_alpha = 0.0;
    BOOST_CHECK(!envelope_bounding_cone(true, R, 0.0, 0.0, z0, z1, Eigen::Vector3d(1.0, 0.0, 2.0),
                                        axis, cos_alpha));
}


BOOST_AUTO_TEST_CASE(front_source_is_unchanged) {
    // In front, the legacy cone exists, so the general option must not move a single bit.
    for (const Eigen::Vector3d src : {Eigen::Vector3d(0.0, 0.0, -25.0), Eigen::Vector3d(8.0, 0.0, -5.0)}) {
        EfficiencyCalculator a, b;
        nai3x3(a);
        nai3x3(b);
        a.set_point_source(src);
        b.set_point_source(src);
        b.enable_general_cone_sampling(true);
        SimulationConfig c = cfg(662.0, 50000, 7);
        c.num_threads = 1;
        const EfficiencyResult ra = a.compute(c);
        const EfficiencyResult rb = b.compute(c);
        BOOST_CHECK_EQUAL(ra.full_energy_peak_efficiency, rb.full_energy_peak_efficiency);
        BOOST_CHECK_EQUAL(ra.total_efficiency, rb.total_efficiency);
    }
}


BOOST_AUTO_TEST_CASE(point_beside_matches_isotropic) {
    check_point(&nai3x3, Eigen::Vector3d(25.0, 0.0, 3.8), 662.0, "NaI 3x3 beside, mid-crystal");
    check_point(&nai3x3, Eigen::Vector3d(20.0, 0.0, -0.02), 122.0, "NaI 3x3 beside, at the face plane");
    check_point(&hpge_bored, Eigen::Vector3d(15.0, 5.0, 5.0), 1332.0, "bored HPGe beside, toward the back");
    check_point(&nai_box, Eigen::Vector3d(0.0, 30.0, 10.0), 662.0, "NaI log beside its long face");
}


BOOST_AUTO_TEST_CASE(extended_beside_matches_isotropic) {
    // A water-free (vacuum-host) cylinder beside the detector: per-position cones, each with its
    // own axis, against the full sphere.
    EfficiencyCalculator iso, gen;
    nai3x3(iso);
    nai3x3(gen);
    for (EfficiencyCalculator* calc : {&iso, &gen})
        calc->set_cylindrical_source(Eigen::Vector3d(20.0, 0.0, 4.0), 2.0, 4.0);
    gen.enable_general_cone_sampling(true);

    const EfficiencyResult r_iso = iso.compute(cfg(662.0, 4000000, 21));
    const EfficiencyResult r_gen = gen.compute(cfg(662.0, 400000, 22));
    const double p_fep = pull(r_iso.full_energy_peak_efficiency, r_iso.fep_uncertainty,
                              r_gen.full_energy_peak_efficiency, r_gen.fep_uncertainty);
    const double p_tot = pull(r_iso.total_efficiency, r_iso.total_uncertainty,
                              r_gen.total_efficiency, r_gen.total_uncertainty);
    BOOST_TEST_MESSAGE("cylinder beside: FEP " << r_iso.full_energy_peak_efficiency << " vs "
                       << r_gen.full_energy_peak_efficiency << " (" << p_fep << " sigma), total "
                       << r_iso.total_efficiency << " vs " << r_gen.total_efficiency << " ("
                       << p_tot << " sigma)");
    BOOST_CHECK_MESSAGE(p_fep < 4.0, "extended FEP differs by " << p_fep << " sigma");
    BOOST_CHECK_MESSAGE(p_tot < 4.0, "extended total differs by " << p_tot << " sigma");
}
