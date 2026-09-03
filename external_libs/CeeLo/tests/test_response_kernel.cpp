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

#define BOOST_TEST_MODULE ResponseKernelTests
#include <boost/test/unit_test.hpp>

// MC-free geometry-kernel self-checks, promoted from the July 2026
// parameterization campaign (its MC truth bank).
// All tests are analytic / quadrature-only -- no Monte Carlo, so they run in
// well under a second.

#include "geometry/Geometry.h"
#include "io/ResponseKernel.h"
#include "test_fep_window.h"
#include "io/SolidAngle.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;

const Material& mat_NaI() { static Material m = make_NaI();      return m; }
const Material& mat_Al()  { static Material m = make_Aluminum(); return m; }
const Material& mat_W()   { static Material m = make_Tungsten(); return m; }

// Bare 3"x3" NaI cylinder, no can, for clean analytic limits.
Geometry bare_nai_3x3() {
    Geometry g;
    g.set_detector(DetectorShape::Cylinder, &mat_NaI(), {3.81, 7.62});
    return g;
}

}  // namespace

// Far-field on-axis: the active solid angle converges to the disk formula
// (side-wall visibility ~0 at 300 cm). Campaign-measured agreement: 0.003%.
BOOST_AUTO_TEST_CASE(far_field_disk_limit) {
    const Geometry g = bare_nai_3x3();
    const double d = 300.0;
    const ApertureQuadrature q = make_aperture_quadrature(g, {0, 0, -d}, 8192);
    const double disk = disk_solid_angle_fraction(d, 3.81);
    BOOST_CHECK_CLOSE(q.omega_frac_active, disk, 0.5 /*percent*/);
}

// Quadrature convergence: N=4096 vs N=16384 at a close point, 45 deg off-axis.
BOOST_AUTO_TEST_CASE(n_ray_convergence) {
    const Geometry g = bare_nai_3x3();
    const Eigen::Vector3d src = source_position(2.0, std::cos(45.0 * kPi / 180.0));
    const ApertureQuadrature q1 = make_aperture_quadrature(g, src, 4096);
    const ApertureQuadrature q2 = make_aperture_quadrature(g, src, 16384);
    BOOST_CHECK_CLOSE(q1.omega_frac_active, q2.omega_frac_active, 0.3);
    BOOST_CHECK_CLOSE(q1.interaction_omega(662.0, MuChoice::Total),
                      q2.interaction_omega(662.0, MuChoice::Total), 0.3);
}

// Saturation limit: at low energy (mu ~ 10+/cm in NaI) every geometric hit
// interacts, so interaction_omega -> active solid angle. Grazing chords keep
// the ratio slightly below 1.
BOOST_AUTO_TEST_CASE(low_energy_saturation) {
    const Geometry g = bare_nai_3x3();
    const ApertureQuadrature q = make_aperture_quadrature(g, {0, 0, -10.0}, 8192);
    const double ratio =
        q.interaction_omega(40.0, MuChoice::Total) / q.omega_frac_active;
    BOOST_CHECK_LE(ratio, 1.0 + 1e-9);
    BOOST_CHECK_GE(ratio, 0.95);
}

// Ordering invariants that hold at any energy:
//   interaction_omega(NoRayleigh) <= interaction_omega(Total) <=
//   omega_frac_active <= cone_omega_frac,
// and the transmitted envelope is bounded by the active solid angle for a
// bare crystal (no passive layers => transmitted == active).
BOOST_AUTO_TEST_CASE(kernel_ordering_invariants) {
    const Geometry g = bare_nai_3x3();
    const ApertureQuadrature q = make_aperture_quadrature(g, {2.0, 0, -6.0}, 8192);
    for (const double E : {40.0, 121.8, 662.0, 2614.0}) {
        const double io_tot = q.interaction_omega(E, MuChoice::Total);
        const double io_nors = q.interaction_omega(E, MuChoice::NoRayleigh);
        BOOST_CHECK_LE(io_nors, io_tot * (1.0 + 1e-12));
        BOOST_CHECK_LE(io_tot, q.omega_frac_active + 1e-12);
    }
    BOOST_CHECK_LE(q.omega_frac_active, q.cone_omega_frac + 1e-12);
    BOOST_CHECK_CLOSE(q.transmitted_omega(662.0), q.omega_frac_active, 1e-6);
}

// A trivial (T_src == 1) source-transmission functor must reproduce
// interaction_omega exactly.
BOOST_AUTO_TEST_CASE(tsrc_identity) {
    const Geometry g = bare_nai_3x3();
    const ApertureQuadrature q = make_aperture_quadrature(g, {0, 0, -5.0}, 4096);
    const double base = q.interaction_omega(200.0, MuChoice::Total);
    const double with_tsrc = q.interaction_omega_with_tsrc(
        200.0, MuChoice::Total, [](const Eigen::Vector3d&) { return 1.0; });
    BOOST_CHECK_CLOSE(base, with_tsrc, 1e-9);
}

// Attenuated can: adding an Al can must strictly reduce the transmitted and
// interaction-weighted solid angles, more so at low energy.
BOOST_AUTO_TEST_CASE(can_attenuation_monotonic) {
    const Geometry bare = bare_nai_3x3();
    Geometry canned;
    canned.set_detector(DetectorShape::Cylinder, &mat_NaI(), {3.81, 7.62});
    canned.add_attenuator(&mat_Al(), 0.05, 0.05, 0.0, 7.62);

    const Eigen::Vector3d src(0, 0, -10.0);
    const ApertureQuadrature qb = make_aperture_quadrature(bare, src, 8192);
    const ApertureQuadrature qc = make_aperture_quadrature(canned, src, 8192);

    const double b60 = qb.interaction_omega(60.0, MuChoice::Total);
    const double c60 = qc.interaction_omega(60.0, MuChoice::Total);
    const double b662 = qb.interaction_omega(662.0, MuChoice::Total);
    const double c662 = qc.interaction_omega(662.0, MuChoice::Total);
    BOOST_CHECK_LT(c60, b60);
    BOOST_CHECK_LT(c662, b662);
    // Fractional loss larger at 60 keV than at 662 keV.
    BOOST_CHECK_LT(c60 / b60, c662 / b662);
}

// Collimator hole fraction: on-axis view down a W collimator bore stays open
// (s ~ 1); a deep off-axis view is shadowed (s small). Uses the campaign Z8
// geometry (3"x3" NaI + 1.5 cm W side collimator extending 5 cm past the face).
BOOST_AUTO_TEST_CASE(collimator_hole_fraction) {
    Geometry g;
    g.set_detector(DetectorShape::Cylinder, &mat_NaI(), {3.81, 7.62});
    g.add_attenuator(&mat_Al(), 0.05, 0.05, 0.0, 7.62);
    g.add_collimator(&mat_W(), 1.5, -5.0, 7.62);

    const ApertureQuadrature on_axis =
        make_aperture_quadrature(g, {0, 0, -30.0}, 8192);
    BOOST_CHECK_GT(on_axis.hole_fraction(662.0), 0.5);

    const Eigen::Vector3d side = source_position(30.0, std::cos(75.0 * kPi / 180.0));
    const ApertureQuadrature off_axis = make_aperture_quadrature(g, side, 8192);
    BOOST_CHECK_LT(off_axis.hole_fraction(200.0), 0.3);
}

// Klein-Nishina in-window fraction: large forward-cone fraction at low E,
// small at high E, monotonically decreasing, bounded in (0, 1).
BOOST_AUTO_TEST_CASE(kn_in_window_fraction_behavior) {
    const double f60 = kn_in_window_fraction(60.0, kTestFepWindowKeV);
    const double f200 = kn_in_window_fraction(200.0, kTestFepWindowKeV);
    const double f662 = kn_in_window_fraction(662.0, kTestFepWindowKeV);
    const double f2614 = kn_in_window_fraction(2614.0, kTestFepWindowKeV);
    BOOST_CHECK_GT(f60, f200);
    BOOST_CHECK_GT(f200, f662);
    BOOST_CHECK_GT(f662, f2614);
    BOOST_CHECK_GT(f60, 0.15);  // wide window cone at 60 keV (theta <= ~37 deg)
    BOOST_CHECK_LT(f2614, 0.02);
    BOOST_CHECK_GT(f2614, 0.0);
}

// Survival removal mu: mu_rem = mu_tot - mu_RS - f_win*mu_CS, strictly below
// mu_total and above the pure-absorption floor mu_tot - mu_RS - mu_CS.
BOOST_AUTO_TEST_CASE(survival_removal_mu_bounds) {
    for (const double E : {60.0, 200.0, 662.0}) {
        const MacroscopicXS xs = mat_Al().macroscopic_xs(E * 1e-3);
        const double f_win = kn_in_window_fraction(E, kTestFepWindowKeV);
        const double mu_rem = fep_survival_removal_mu(xs, f_win);
        BOOST_CHECK_LT(mu_rem, xs.mu_total());
        BOOST_CHECK_GT(mu_rem, xs.mu_total() - xs.mu_rs - xs.mu_cs - 1e-12);
    }
}

// Depth-aware survival credit (stage E3 A2): g(tau_src) = exp(-tau/tau_c) in
// (0,1], g(0)=1, strictly decreasing; and the depth-aware removal mu reduces
// EXACTLY to fep_survival_removal_mu at tau_src=0 (thin-source identity) while
// converging to the pure mu_tot - mu_RS floor as tau_src -> infinity.
BOOST_AUTO_TEST_CASE(fep_depth_survival_credit_behavior) {
    BOOST_CHECK_EQUAL(fep_depth_survival_credit(0.0), 1.0);
    double prev = fep_depth_survival_credit(0.0);
    for (const double tau : {0.5, 1.0, 2.0, 5.0, 11.0, 30.0}) {
        const double g = fep_depth_survival_credit(tau);
        BOOST_CHECK_GT(g, 0.0);
        BOOST_CHECK_LT(g, prev);          // strictly decreasing
        prev = g;
    }
    BOOST_CHECK_LT(fep_depth_survival_credit(60.0), 1e-4);  // deep => credit gone
}

BOOST_AUTO_TEST_CASE(fep_survival_removal_mu_depth_identity_and_bounds) {
    const Material fe = make_Iron();
    const Material h2o = make_Water();
    for (const double E : {60.0, 122.0, 662.0}) {
        for (const Material* m : {&fe, &h2o}) {
            const MacroscopicXS xs = m->macroscopic_xs(E * 1e-3);
            const double f_win = kn_in_window_fraction(E, kTestFepWindowKeV, *m);
            const double flat = fep_survival_removal_mu(xs, f_win);
            // tau_src = 0 must be bit-identical to the flat credit.
            BOOST_CHECK_EQUAL(fep_survival_removal_mu_depth(xs, f_win, 0.0), flat);
            // Growing tau_src removes the Compton credit => mu_rem rises toward
            // the mu_tot - mu_RS floor, staying between flat and that floor.
            const double floor = xs.mu_total() - xs.mu_rs;
            double prev = flat;
            for (const double tau : {1.0, 3.0, 8.0, 25.0}) {
                const double md = fep_survival_removal_mu_depth(xs, f_win, tau);
                BOOST_CHECK_GE(md, prev - 1e-15);   // non-decreasing in tau
                BOOST_CHECK_LE(md, floor + 1e-12);
                prev = md;
            }
            BOOST_CHECK_CLOSE(fep_survival_removal_mu_depth(xs, f_win, 1e4),
                              floor, 1e-6);
        }
    }
}

// S(x,Z)-corrected in-window fraction (stage E1): bound electrons suppress
// forward scatter, so the material-aware fraction is a proper fraction BELOW
// the free-electron value, ordered by Z at low energy.
BOOST_AUTO_TEST_CASE(kn_in_window_material_below_free_electron) {
    const Material fe = make_Iron();
    const Material pb = make_Lead();
    const Material h2o = make_Water();
    for (const double E : {60.0, 122.0, 344.0, 662.0, 1332.0, 2614.0}) {
        const double f_free = kn_in_window_fraction(E, kTestFepWindowKeV);
        for (const Material* m : {&h2o, &fe, &pb}) {
            const double f_mat = kn_in_window_fraction(E, kTestFepWindowKeV, *m);
            BOOST_CHECK_GT(f_mat, 0.0);
            BOOST_CHECK_LT(f_mat, 1.0);
            BOOST_CHECK_LT(f_mat, f_free);
        }
    }
    // Suppression strengthens with Z (measured @60 keV: water 0.89, Fe 0.77,
    // Pb 0.63 of the free-electron fraction).
    const double f_free60 = kn_in_window_fraction(60.0, kTestFepWindowKeV);
    const double r_h2o = kn_in_window_fraction(60.0, kTestFepWindowKeV, h2o) / f_free60;
    const double r_fe = kn_in_window_fraction(60.0, kTestFepWindowKeV, fe) / f_free60;
    const double r_pb = kn_in_window_fraction(60.0, kTestFepWindowKeV, pb) / f_free60;
    BOOST_CHECK_GT(r_h2o, r_fe);
    BOOST_CHECK_GT(r_fe, r_pb);
    BOOST_CHECK_LT(r_fe, 0.85);  // material correction is real at 60 keV
    BOOST_CHECK_LT(r_pb, 0.72);
}

// High-energy convergence: the window-edge momentum transfer is energy-
// independent (~1.6 A^-1) so the f ratio does NOT approach 1, but f itself
// -> 0, so the ABSOLUTE correction and its effect on the survival removal
// mu vanish: mu_rem(material-aware) within 0.1% of free-electron >= 1332 keV.
BOOST_AUTO_TEST_CASE(kn_in_window_material_high_E_convergence) {
    const Material fe = make_Iron();
    const Material pb = make_Lead();
    for (const double E : {1332.0, 2614.0}) {
        const double f_free = kn_in_window_fraction(E, kTestFepWindowKeV);
        for (const Material* m : {&fe, &pb}) {
            const double f_mat = kn_in_window_fraction(E, kTestFepWindowKeV, *m);
            BOOST_CHECK_LT(f_free - f_mat, 1.0e-3);  // absolute gap vanishes
            const MacroscopicXS xs = m->macroscopic_xs(E * 1e-3);
            const double mr_free = fep_survival_removal_mu(xs, f_free);
            const double mr_mat = fep_survival_removal_mu(xs, f_mat);
            BOOST_CHECK_CLOSE(mr_free, mr_mat, 0.1 /*percent*/);
        }
    }
}

// Window-size behavior: monotonically non-decreasing in win_keV and -> the
// full incoherent fraction (1.0) as the window swallows the whole spectrum.
BOOST_AUTO_TEST_CASE(kn_in_window_material_window_monotonic) {
    const Material fe = make_Iron();
    for (const double E : {60.0, 344.0}) {
        double prev = 0.0;
        for (const double win : {0.5, 1.5, 3.0, 6.0, 12.0}) {
            const double f = kn_in_window_fraction(E, win, fe);
            BOOST_CHECK_GE(f, prev - 1e-12);
            BOOST_CHECK_LE(f, 1.0 + 1e-12);
            prev = f;
        }
        // Window wider than the maximum Compton energy loss => every
        // incoherent scatter stays in-window.
        BOOST_CHECK_CLOSE(kn_in_window_fraction(E, E, fe), 1.0, 1e-6);
    }
}
