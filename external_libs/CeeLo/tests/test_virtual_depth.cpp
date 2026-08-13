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

#define BOOST_TEST_MODULE VirtualDepthTests
#include <boost/test/unit_test.hpp>

#include "io/SolidAngle.h"
#include "io/VirtualDepthFit.h"
#include "materials/Material.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace ceelo;

namespace {
const std::vector<double> kLadder = {1, 2, 3, 5, 8, 12, 20, 30};
const double kR = 3.81; // NaI 3x3 crystal radius

// Synthesize an exact ideal efficiency scan for a known (intrinsic, delta).
std::vector<double> ideal_scan(double intrinsic, double delta, double R,
                               const std::vector<double>& d) {
    std::vector<double> eps;
    for (double x : d) eps.push_back(intrinsic * disk_solid_angle_fraction(x + delta, R));
    return eps;
}
} // namespace

// ---- Test 1: synthetic delta recovery (no MC) ----
BOOST_AUTO_TEST_SUITE(SyntheticRecovery)

BOOST_AUTO_TEST_CASE(recovers_injected_delta_exact) {
    const double intrinsic = 1.234e-2;
    const double delta_true = 0.80;
    const auto eps = ideal_scan(intrinsic, delta_true, kR, kLadder);
    const std::vector<double> unc(kLadder.size(), 0.0);

    VpdPoint pt = fit_delta_one_energy(kR, kLadder, eps, unc, 6.0, 0.02);
    // Recovered within one grid step.
    BOOST_CHECK_SMALL(pt.delta_cm - delta_true, 0.021);
    BOOST_CHECK_CLOSE(pt.intrinsic_mean, intrinsic, 1e-3);
    BOOST_CHECK_SMALL(pt.max_abs_res, 1e-4); // an exact model has ~zero residual
}

BOOST_AUTO_TEST_CASE(recovers_with_noise) {
    const double intrinsic = 5.0e-3;
    const double delta_true = 1.20;
    auto eps = ideal_scan(intrinsic, delta_true, kR, kLadder);
    std::vector<double> unc(kLadder.size(), 0.0);

    std::mt19937_64 rng(42);
    std::normal_distribution<double> noise(0.0, 0.01); // 1% relative jitter
    for (size_t i = 0; i < eps.size(); ++i) {
        eps[i] *= (1.0 + noise(rng));
        unc[i] = 0.01 * eps[i];
    }

    VpdPoint pt = fit_delta_one_energy(kR, kLadder, eps, unc, 6.0, 0.02);
    // With 1% noise the argmin should still land within a few grid steps.
    BOOST_CHECK_SMALL(pt.delta_cm - delta_true, 0.15);
    BOOST_CHECK(pt.cv_residual < 0.05); // small spread => good fit
}

BOOST_AUTO_TEST_SUITE_END()

// ---- Test 2: exp-of-log fit + file round-trip ----
BOOST_AUTO_TEST_SUITE(ExpLogRoundTrip)

BOOST_AUTO_TEST_CASE(coeffs_roundtrip) {
    const std::vector<double> coeffs = {-0.6, 0.5, -0.04}; // delta=exp(p0+p1 u+p2 u^2)
    const std::vector<double> energies = {122, 344, 662, 1173, 1332};
    std::vector<double> deltas;
    for (double E : energies) deltas.push_back(eval_exp_log(coeffs, E));

    auto refit = fit_exp_log_coeffs(energies, deltas, 3);
    BOOST_REQUIRE_EQUAL(refit.size(), coeffs.size());
    for (size_t j = 0; j < coeffs.size(); ++j)
        BOOST_CHECK_SMALL(refit[j] - coeffs[j], 1e-6);

    // Re-evaluation reproduces the deltas.
    for (size_t k = 0; k < energies.size(); ++k)
        BOOST_CHECK_CLOSE(eval_exp_log(refit, energies[k]), deltas[k], 1e-4);
}

BOOST_AUTO_TEST_CASE(file_save_load_roundtrip) {
    // Build a small VpdFit by hand and round-trip it through the text format.
    VpdFit fit;
    fit.detector = make_descriptor("Test NaI", DetectorShape::Cylinder,
                                   make_NaI(), {kR, 7.62});
    const std::vector<double> energies = {122, 662, 1332};
    fit.coeffs = {-0.6, 0.5, -0.04};
    for (double E : energies) {
        VpdPoint p;
        p.energy_keV = E;
        p.delta_cm = eval_exp_log(fit.coeffs, E);
        p.intrinsic_mean = 1e-2;
        p.cv_residual = 0.003;
        p.max_abs_res = 0.004;
        p.dist_cm = {1, 5, 10};
        p.eps = ideal_scan(p.intrinsic_mean, p.delta_cm, kR, p.dist_cm);
        p.eps_unc = {1e-5, 2e-6, 8e-7};
        fit.points.push_back(p);
        fit.log_fit_residuals.push_back(0.0);
    }
    fit.valid_e_min_keV = 122;
    fit.valid_e_max_keV = 1332;
    fit.rms_log_residual = 1e-9;

    const std::string fname = "test_vpd_roundtrip.txt";
    BOOST_REQUIRE(save_vpd(fit, fname));
    VpdFit loaded = load_vpd(fname);

    BOOST_CHECK_EQUAL(loaded.detector.name, "Test NaI");
    BOOST_CHECK_EQUAL(loaded.detector.material_name, "NaI");
    BOOST_CHECK(loaded.detector.shape == DetectorShape::Cylinder);
    BOOST_CHECK_CLOSE(loaded.detector.crystal_radius_cm, kR, 1e-4);
    BOOST_CHECK_CLOSE(loaded.valid_e_min_keV, 122.0, 1e-4);
    BOOST_CHECK_CLOSE(loaded.valid_e_max_keV, 1332.0, 1e-4);

    BOOST_REQUIRE_EQUAL(loaded.coeffs.size(), fit.coeffs.size());
    for (size_t j = 0; j < fit.coeffs.size(); ++j)
        BOOST_CHECK_SMALL(loaded.coeffs[j] - fit.coeffs[j], 1e-6);

    BOOST_REQUIRE_EQUAL(loaded.points.size(), fit.points.size());
    for (size_t k = 0; k < fit.points.size(); ++k) {
        BOOST_CHECK_CLOSE(loaded.points[k].energy_keV, fit.points[k].energy_keV, 1e-4);
        BOOST_CHECK_CLOSE(loaded.points[k].delta_cm, fit.points[k].delta_cm, 1e-3);
        BOOST_REQUIRE_EQUAL(loaded.points[k].dist_cm.size(), fit.points[k].dist_cm.size());
        for (size_t i = 0; i < fit.points[k].eps.size(); ++i)
            BOOST_CHECK_CLOSE(loaded.points[k].eps[i], fit.points[k].eps[i], 1e-3);
    }

    // delta_at evaluates the fitted model and clamps outside the valid range.
    BOOST_CHECK_CLOSE(loaded.delta_at(662.0), eval_exp_log(fit.coeffs, 662.0), 1e-4);
    BOOST_CHECK_CLOSE(loaded.delta_at(50.0), eval_exp_log(fit.coeffs, 122.0), 1e-4); // clamped low
    BOOST_CHECK_CLOSE(loaded.delta_at(5000.0), eval_exp_log(fit.coeffs, 1332.0), 1e-4); // clamped high
}

BOOST_AUTO_TEST_SUITE_END()

// ---- Test 3: degeneracy guard ----
BOOST_AUTO_TEST_SUITE(Degeneracy)

BOOST_AUTO_TEST_CASE(flat_efficiency_no_nan) {
    // A nearly distance-flat efficiency (tiny R, far ladder) makes the CV
    // objective flat in delta. The fit must not NaN and should report a finite
    // (and, here, small) CV without crashing.
    const double tinyR = 0.01;
    const std::vector<double> far = {20, 25, 30, 40, 50};
    std::vector<double> eps(far.size(), 1e-4); // identical -> perfectly flat
    std::vector<double> unc(far.size(), 1e-6);

    VpdPoint pt = fit_delta_one_energy(tinyR, far, eps, unc, 6.0, 0.02);
    BOOST_CHECK(std::isfinite(pt.delta_cm));
    BOOST_CHECK(std::isfinite(pt.intrinsic_mean));
    BOOST_CHECK(std::isfinite(pt.cv_residual));
}

BOOST_AUTO_TEST_SUITE_END()

// ---- Test 4: physics sanity (small-stat MC) ----
BOOST_AUTO_TEST_SUITE(PhysicsSanity)

BOOST_AUTO_TEST_CASE(delta_increases_with_energy) {
    // delta(E) should grow with energy: interactions move deeper into the crystal.
    // Keep events low -- this is the only MC-bearing test.
    Material nai = make_NaI();
    DetectorDescriptor desc = make_descriptor("NaI 3x3", DetectorShape::Cylinder,
                                              nai, {kR, 7.62});
    VpdFitConfig cfg;
    cfg.energies_keV = {122.0, 1332.0};
    cfg.distance_ladder_cm = {1, 2, 3, 5, 8, 12};
    cfg.events_per_point = 80000;
    cfg.num_threads = 0;
    cfg.poly_order = 2;

    VpdFit fit = fit_virtual_depth(desc, nai, {kR, 7.62}, cfg);
    BOOST_REQUIRE_EQUAL(fit.points.size(), 2u);
    const double d122 = fit.points[0].delta_cm;
    const double d1332 = fit.points[1].delta_cm;
    BOOST_TEST_MESSAGE("delta(122)=" << d122 << " cm, delta(1332)=" << d1332 << " cm");
    // Margin well above the CV grid step (0.02) and statistical noise.
    BOOST_CHECK(d1332 > d122 + 0.2);
}

BOOST_AUTO_TEST_SUITE_END()
