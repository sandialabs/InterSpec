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

#define BOOST_TEST_MODULE BiasingConsistencyTests
#include <boost/test/unit_test.hpp>

/// @file test_biasing_consistency.cpp
/// @brief Gate A consistency tests for event-biasing variance reduction:
/// biased and analog estimators must agree within statistics (5 sigma),
/// plus analytic unit checks of the forced-collision weight.

#include "efficiency/EfficiencyCalculator.h"
#include "geometry/Geometry.h"
#include "materials/Material.h"
#include "transport/ComptonScatter.h"
#include "transport/PhotonTransport.h"

#include <cmath>
#include <random>

using namespace ceelo;

namespace {

/// SimulationConfig targeting 1% FEP relative precision (fast Gate A runs).
SimulationConfig precision_config(double energy_keV, double target = 0.01) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.max_events = 50000000;
    config.termination.min_events = 5000;
    config.termination.max_wall_seconds = 30.0;
    config.num_threads = 0;
    config.batch_size = 10000;
    return config;
}

/// 5-sigma z-test between two independent estimates.
void check_consistent(double a, double sig_a, double b, double sig_b,
                      const char* what) {
    double sigma = std::sqrt(sig_a * sig_a + sig_b * sig_b);
    double diff = std::abs(a - b);
    BOOST_TEST_INFO(what << ": a=" << a << " +- " << sig_a
                         << ", b=" << b << " +- " << sig_b
                         << ", z=" << (sigma > 0 ? diff / sigma : 0.0));
    BOOST_CHECK_LT(diff, 5.0 * sigma + 1e-9);
}

/// Run forced-vs-analog Gate A on a configured calculator at one energy.
/// The calculator is configured by the caller; this runs both modes.
void gate_a_forced(EfficiencyCalculator& calc, double energy_keV,
                   const char* what) {
    BiasingConfig analog;
    calc.set_biasing(analog);
    auto res_analog = calc.compute(precision_config(energy_keV));

    BiasingConfig forced;
    forced.force_detector_interaction = true;
    calc.set_biasing(forced);
    auto res_forced = calc.compute(precision_config(energy_keV));

    check_consistent(res_analog.full_energy_peak_efficiency,
                     res_analog.fep_uncertainty,
                     res_forced.full_energy_peak_efficiency,
                     res_forced.fep_uncertainty, what);
    check_consistent(res_analog.total_efficiency,
                     res_analog.total_uncertainty,
                     res_forced.total_efficiency,
                     res_forced.total_uncertainty, what);
}

} // anonymous namespace

// ============================================================
//  Forced first interaction: analytic unit checks
// ============================================================

BOOST_AUTO_TEST_SUITE(ForcedCollisionUnit)

BOOST_AUTO_TEST_CASE(forced_weight_matches_analytic_on_axis) {
    // On-axis photon through a bare 3"x3" NaI cylinder: every history sees
    // the same chord (7.62 cm), so the forced-collision weight must equal
    // 1 - exp(-mu_total * L) exactly, event by event.
    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    const double energy_keV = 662.0;
    MacroscopicXS xs = nai.macroscopic_xs(energy_keV * 1e-3);
    // PP masked below 1022 keV in transport; Rayleigh enabled by default.
    const double mu = xs.mu_pe + xs.mu_cs + xs.mu_rs;
    const double expected_w = 1.0 - std::exp(-mu * 7.62);

    TransportConfig config;
    config.force_first_interaction = true;

    std::mt19937_64 rng(12345);
    for (int i = 0; i < 50; ++i) {
        auto res = transport_photon(
            Eigen::Vector3d(0.0, 0.0, -5.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            energy_keV, geom, config, rng);
        BOOST_CHECK_CLOSE(res.weight, expected_w, 1e-6);
        // Forced events always interact somewhere on the path
        BOOST_CHECK(res.any_interaction_in_scoring ||
                    res.energy_deposited_total > 0.0 || res.escaped);
        BOOST_CHECK(res.any_interaction_in_scoring);
    }
}

BOOST_AUTO_TEST_CASE(forced_weight_is_one_when_ray_misses) {
    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    TransportConfig config;
    config.force_first_interaction = true;

    std::mt19937_64 rng(999);
    // Photon heading away from the detector: no material on the ray
    auto res = transport_photon(
        Eigen::Vector3d(0.0, 0.0, -5.0), Eigen::Vector3d(0.0, 0.0, -1.0),
        662.0, geom, config, rng);
    BOOST_CHECK_CLOSE(res.weight, 1.0, 1e-12);
    BOOST_CHECK(res.escaped);
    BOOST_CHECK(!res.any_interaction_in_scoring);
}

BOOST_AUTO_TEST_CASE(forced_mean_deposit_matches_analog) {
    // Distributional check of the forced first-interaction depth: the
    // weighted mean scoring deposit over forced histories must match the
    // analog mean deposit (estimator of the same physical quantity).
    Material nai = make_NaI();
    Geometry geom;
    geom.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

    const double energy_keV = 2614.0;
    const int n = 200000;

    TransportConfig analog;
    TransportConfig forced;
    forced.force_first_interaction = true;

    std::mt19937_64 rng_a(2024), rng_f(4048);
    double sum_a = 0.0, sum_a2 = 0.0, sum_f = 0.0, sum_f2 = 0.0;
    for (int i = 0; i < n; ++i) {
        auto ra = transport_photon(
            Eigen::Vector3d(0.3, -0.2, -5.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            energy_keV, geom, analog, rng_a);
        sum_a += ra.energy_deposited_scoring;
        sum_a2 += ra.energy_deposited_scoring * ra.energy_deposited_scoring;

        auto rf = transport_photon(
            Eigen::Vector3d(0.3, -0.2, -5.0), Eigen::Vector3d(0.0, 0.0, 1.0),
            energy_keV, geom, forced, rng_f);
        double score = rf.weight * rf.energy_deposited_scoring;
        sum_f += score;
        sum_f2 += score * score;
    }
    double mean_a = sum_a / n, mean_f = sum_f / n;
    double sig_a = std::sqrt((sum_a2 / n - mean_a * mean_a) / n);
    double sig_f = std::sqrt((sum_f2 / n - mean_f * mean_f) / n);
    check_consistent(mean_a, sig_a, mean_f, sig_f, "mean scoring deposit");
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Gate A: forced vs analog efficiency consistency (5 sigma)
// ============================================================

BOOST_AUTO_TEST_SUITE(ForcedCollisionGateA)

BOOST_AUTO_TEST_CASE(nai_bare_662) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    gate_a_forced(calc, 662.0, "config1 662 keV");
}

BOOST_AUTO_TEST_CASE(nai_bare_2614) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    gate_a_forced(calc, 2614.0, "config1 2614 keV");
}

BOOST_AUTO_TEST_CASE(czt_box_662) {
    Material czt = make_CZT();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Box, &czt, {0.5, 0.5, 0.5});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
    gate_a_forced(calc, 662.0, "config5 662 keV");
}

BOOST_AUTO_TEST_CASE(nai_al_pb_attenuators_662) {
    // Multi-layer attenuators exercise the multi-segment tau ladder
    Material nai = make_NaI();
    Material al = make_Aluminum();
    Material pb = make_Lead();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.add_attenuator(&al, 0.1, 0.1, 0.0, 7.62);
    calc.add_attenuator(&pb, 0.2, 0.2, 0.0, 7.62);
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -15.0));
    gate_a_forced(calc, 662.0, "config7 662 keV");
}

BOOST_AUTO_TEST_CASE(fep_only_mode_2614) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.enable_fep_only_mode(true);

    BiasingConfig analog;
    calc.set_biasing(analog);
    auto res_analog = calc.compute(precision_config(2614.0));

    BiasingConfig forced;
    forced.force_detector_interaction = true;
    calc.set_biasing(forced);
    auto res_forced = calc.compute(precision_config(2614.0));

    check_consistent(res_analog.full_energy_peak_efficiency,
                     res_analog.fep_uncertainty,
                     res_forced.full_energy_peak_efficiency,
                     res_forced.fep_uncertainty, "fep-only 2614 keV");
}

BOOST_AUTO_TEST_CASE(fe_source_shield_662) {
    // Source-effects (isotropic full-mode) path: forced primary initial
    // detector transport for a non-Marinelli shielded point source
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);
    gate_a_forced(calc, 662.0, "config11 662 keV");
}

BOOST_AUTO_TEST_CASE(fe_source_shield_2000) {
    // Above the pair-production threshold: shield PP produces 511 keV
    // annihilation-gamma secondaries, so forced runs exercise the per-event
    // force gating (events with secondaries transport the primary analog).
    // Guards estimator blow-ups; the 4M-event validation sweep is the
    // fine-grained bias check.
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);
    gate_a_forced(calc, 2000.0, "config11 2000 keV (PP-secondary gating)");
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Gate A: mixture angular biasing vs isotropic (5 sigma)
// ============================================================

namespace {

/// Run mixture-vs-isotropic Gate A on a configured source-effects calculator.
void gate_a_mixture(EfficiencyCalculator& calc, double energy_keV,
                    double alpha, const char* what) {
    BiasingConfig iso;
    calc.set_biasing(iso);
    auto res_iso = calc.compute(precision_config(energy_keV));

    BiasingConfig mix;
    mix.mixture_cone_alpha = alpha;
    calc.set_biasing(mix);
    auto res_mix = calc.compute(precision_config(energy_keV));

    check_consistent(res_iso.full_energy_peak_efficiency,
                     res_iso.fep_uncertainty,
                     res_mix.full_energy_peak_efficiency,
                     res_mix.fep_uncertainty, what);
    check_consistent(res_iso.total_efficiency,
                     res_iso.total_uncertainty,
                     res_mix.total_efficiency,
                     res_mix.total_uncertainty, what);
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(MixtureBiasingGateA)

BOOST_AUTO_TEST_CASE(fe_source_shield_multi_energy) {
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);
    gate_a_mixture(calc, 100.0, 0.2, "config11 100 keV mixture");
    gate_a_mixture(calc, 662.0, 0.2, "config11 662 keV mixture");
    gate_a_mixture(calc, 2614.0, 0.2, "config11 2614 keV mixture");
}

BOOST_AUTO_TEST_CASE(steel_box_cellulose_662) {
    // Shrunken config-12-style box source for test speed
    Material nai = make_NaI();
    Material ss = make_StainlessSteel304();
    Material cellulose = make_Cellulose();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(
        Eigen::Vector3d(0.0, 0.0, -12.0),
        Eigen::Vector3d(2.5, 2.5, 4.0),
        Eigen::Matrix3d::Identity());
    calc.set_source_material(&cellulose);
    calc.add_source_shield(&ss, 0.2);
    gate_a_mixture(calc, 662.0, 0.2, "box source 662 keV mixture");
}

BOOST_AUTO_TEST_CASE(mixture_alpha_one_equals_isotropic) {
    // alpha = 1 makes the mixture pdf exactly isotropic; every weight is 1
    // and results must agree with the analog isotropic path.
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);
    gate_a_mixture(calc, 662.0, 1.0, "config11 662 keV alpha=1");
}

BOOST_AUTO_TEST_CASE(auto_enable_policy) {
    // Pin the compute_effective_biasing() heuristics: biasing must be ON
    // where it measurably wins and OFF where it measurably loses.
    SimulationConfig fixed_n;          // no precision targets: both matter
    fixed_n.termination.max_events = 1000;
    SimulationConfig fep_only_target;  // only FEP precision targeted
    fep_only_target.termination.target_fep_rel_precision = 0.005;
    fep_only_target.termination.max_events = 1000;
    SimulationConfig total_only_target;
    total_only_target.termination.target_total_rel_precision = 0.005;
    total_only_target.termination.max_events = 1000;

    // Thick NaI (config-1-like): analog at 662 keV (tau ~ 2, forcing loses
    // FEP FOM); forcing allowed when only total efficiency is targeted.
    {
        Material nai = make_NaI();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
        fixed_n.energy_keV = 662.0;
        auto b = calc.compute_effective_biasing(fixed_n);
        BOOST_CHECK(!b.force_detector_interaction);
        BOOST_CHECK_EQUAL(b.mixture_cone_alpha, 0.0);  // no source effects
        total_only_target.energy_keV = 662.0;
        b = calc.compute_effective_biasing(total_only_target);
        BOOST_CHECK(b.force_detector_interaction);
    }

    // Thin CZT (config-5-like): forced at 662 keV (tau ~ 0.4), analog at
    // 60 keV (tau >> 1, photoelectric-dominated).
    {
        Material czt = make_CZT();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Box, &czt, {0.5, 0.5, 0.5});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
        fixed_n.energy_keV = 662.0;
        BOOST_CHECK(calc.compute_effective_biasing(fixed_n).force_detector_interaction);
        fixed_n.energy_keV = 60.0;
        BOOST_CHECK(!calc.compute_effective_biasing(fixed_n).force_detector_interaction);
    }

    // Source-effect config (config-11-like): two-stream estimator with
    // target-dependent f and Compton-bias gamma (mixture alpha is no
    // longer auto-enabled); manual set_biasing() overrides everything.
    {
        Material nai = make_NaI();
        Material fe = make_Iron();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
        calc.add_source_shield(&fe, 0.5);
        fixed_n.energy_keV = 662.0;
        auto b = calc.compute_effective_biasing(fixed_n);
        BOOST_CHECK(b.two_stream);
        BOOST_CHECK_EQUAL(b.two_stream_direct_fraction, 0.25);
        BOOST_CHECK_EQUAL(b.compton_cone_fraction, 0.3);
        BOOST_CHECK_EQUAL(b.mixture_cone_alpha, 0.0);
        fep_only_target.energy_keV = 662.0;
        b = calc.compute_effective_biasing(fep_only_target);
        BOOST_CHECK(b.two_stream);
        BOOST_CHECK_EQUAL(b.two_stream_direct_fraction, 0.5);
        BOOST_CHECK_EQUAL(b.compton_cone_fraction, 0.0);

        calc.set_biasing(BiasingConfig{});  // manual analog override
        b = calc.compute_effective_biasing(fixed_n);
        BOOST_CHECK(!b.force_detector_interaction);
        BOOST_CHECK(!b.two_stream);
        BOOST_CHECK_EQUAL(b.mixture_cone_alpha, 0.0);
        BOOST_CHECK_EQUAL(b.compton_cone_fraction, 0.0);
    }

    // Close geometry (water puck ~1 cm from the face): the detector
    // subtends omega/4pi >= 0.15 from the source center, where the
    // two-stream split measurably loses — must stay analog.
    {
        Material nai = make_NaI();
        Material water = make_Water();
        Material pe = make_Polyethylene();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -1.6), 2.5, 0.5);
        calc.set_source_material(&water);
        calc.add_source_shield(&pe, 0.1);
        fixed_n.energy_keV = 662.0;
        auto b = calc.compute_effective_biasing(fixed_n);
        BOOST_CHECK(!b.two_stream);
        BOOST_CHECK_EQUAL(b.compton_cone_fraction, 0.0);
    }

    // Marinelli (config-8-like): always fully analog under auto-enable.
    {
        Material nai = make_NaI();
        Material al = make_Aluminum();
        Material water = make_Water();
        Material pe = make_Polyethylene();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
        calc.set_marinelli_beaker(4.3, 6.0, 7.5, 4.0, 0.5, &water, &pe, 0.2);
        fixed_n.energy_keV = 2614.0;
        auto b = calc.compute_effective_biasing(fixed_n);
        BOOST_CHECK(!b.force_detector_interaction);
        BOOST_CHECK_EQUAL(b.mixture_cone_alpha, 0.0);
        BOOST_CHECK(!b.two_stream);
        BOOST_CHECK_EQUAL(b.compton_cone_fraction, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(marinelli_with_mixture_is_unbiased) {
    // Marinelli positions surround the crystal: the per-position cone is
    // degenerate for most positions and the mixture must reduce to a
    // near-no-op while staying unbiased.
    Material nai = make_NaI();
    Material al = make_Aluminum();
    Material water = make_Water();
    Material pe = make_Polyethylene();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
    calc.set_marinelli_beaker(4.3, 6.0, 7.5, 4.0, 0.5, &water, &pe, 0.2);
    gate_a_mixture(calc, 662.0, 0.2, "marinelli 662 keV mixture");
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Two-stream stratified estimator: unit checks + Gate A
// ============================================================

namespace {

/// Run two-stream-vs-analog Gate A on a configured source-effects calculator.
void gate_a_two_stream(EfficiencyCalculator& calc, double energy_keV,
                       double f, bool with_forcing, const char* what) {
    BiasingConfig analog;
    calc.set_biasing(analog);
    auto res_a = calc.compute(precision_config(energy_keV));

    BiasingConfig ts;
    ts.two_stream = true;
    ts.two_stream_direct_fraction = f;
    ts.force_detector_interaction = with_forcing;
    calc.set_biasing(ts);
    auto res_t = calc.compute(precision_config(energy_keV));

    check_consistent(res_a.full_energy_peak_efficiency,
                     res_a.fep_uncertainty,
                     res_t.full_energy_peak_efficiency,
                     res_t.fep_uncertainty, what);
    check_consistent(res_a.total_efficiency,
                     res_a.total_uncertainty,
                     res_t.total_efficiency,
                     res_t.total_uncertainty, what);
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(TwoStreamUnit)

BOOST_AUTO_TEST_CASE(no_interaction_probability_matches_analytic) {
    // The direct-stream transmission weight must equal exp(-sum mu_tot*l)
    // with the full (unmasked) mu_total, for each source shape.
    const double e_keV = 662.0;
    const double e_MeV = e_keV * 1e-3;

    // Point source + 0.5 cm Fe spherical shell: path = 0.5 cm any direction.
    {
        Material nai = make_NaI();
        Material fe = make_Iron();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
        calc.add_source_shield(&fe, 0.5);
        const double mu = fe.macroscopic_xs(e_MeV).mu_total();
        double path = 0.0;
        double T = calc.source_geometry().no_interaction_probability(
            Eigen::Vector3d(0.0, 0.0, -10.0),
            Eigen::Vector3d(0.0, 0.0, 1.0), e_keV, &path);
        BOOST_CHECK_CLOSE(T, std::exp(-mu * 0.5), 1e-6);
        BOOST_CHECK_CLOSE(path, 0.5, 1e-6);
    }

    // Rectangular cellulose box + 0.2 cm SS304, ray along +z from center:
    // 4.0 cm cellulose + 0.2 cm steel.
    {
        Material nai = make_NaI();
        Material ss = make_StainlessSteel304();
        Material cellulose = make_Cellulose();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -12.0),
                                    Eigen::Vector3d(2.5, 2.5, 4.0),
                                    Eigen::Matrix3d::Identity());
        calc.set_source_material(&cellulose);
        calc.add_source_shield(&ss, 0.2);
        const double mu_c = cellulose.macroscopic_xs(e_MeV).mu_total();
        const double mu_s = ss.macroscopic_xs(e_MeV).mu_total();
        double path = 0.0;
        double T = calc.source_geometry().no_interaction_probability(
            Eigen::Vector3d(0.0, 0.0, -12.0),
            Eigen::Vector3d(0.0, 0.0, 1.0), e_keV, &path);
        BOOST_CHECK_CLOSE(T, std::exp(-mu_c * 4.0 - mu_s * 0.2), 1e-6);
        BOOST_CHECK_CLOSE(path, 4.2, 1e-6);
    }

    // Cylindrical water source, ray along the axis from center:
    // half_length of water (no shield).
    {
        Material nai = make_NaI();
        Material water = make_Water();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0),
                                    2.0, 1.5);
        calc.set_source_material(&water);
        const double mu_w = water.macroscopic_xs(e_MeV).mu_total();
        double path = 0.0;
        double T = calc.source_geometry().no_interaction_probability(
            Eigen::Vector3d(0.0, 0.0, -10.0),
            Eigen::Vector3d(0.0, 0.0, 1.0), e_keV, &path);
        BOOST_CHECK_CLOSE(T, std::exp(-mu_w * 1.5), 1e-6);
        BOOST_CHECK_CLOSE(path, 1.5, 1e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TwoStreamGateA)

BOOST_AUTO_TEST_CASE(fe_source_shield_662) {
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);
    gate_a_two_stream(calc, 662.0, 0.5, false, "two-stream cfg11 662 keV");
    gate_a_two_stream(calc, 662.0, 0.85, false,
                      "two-stream cfg11 662 keV f=0.85");
}

BOOST_AUTO_TEST_CASE(fe_source_shield_2000_pp_secondaries) {
    // Above the PP threshold: shield PP 511-keV secondaries and escaped
    // electrons exercise the interacted-classification (all secondary-
    // producing events must land in the scatter stream) and the per-event
    // forcing gate.
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);
    calc.enable_source_electron_transport(true);
    gate_a_two_stream(calc, 2000.0, 0.5, true,
                      "two-stream cfg11 2000 keV (PP + forcing)");
}

BOOST_AUTO_TEST_CASE(steel_box_cellulose_662) {
    Material nai = make_NaI();
    Material ss = make_StainlessSteel304();
    Material cellulose = make_Cellulose();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -12.0),
                                Eigen::Vector3d(2.5, 2.5, 4.0),
                                Eigen::Matrix3d::Identity());
    calc.set_source_material(&cellulose);
    calc.add_source_shield(&ss, 0.2);
    gate_a_two_stream(calc, 662.0, 0.5, false, "two-stream box 662 keV");
}

BOOST_AUTO_TEST_CASE(low_energy_fep_s_channel) {
    // At low energy, small-angle Compton in the source loses < 1.5 keV and
    // the degraded primary still lands inside the FEP window: the scatter
    // stream must carry that FEP contribution (measured ~24% of FEP for a
    // cellulose box at 59 keV). A two-stream estimator that assumed
    // FEP_s = 0 would fail this gate.
    Material nai = make_NaI();
    Material ss = make_StainlessSteel304();
    Material cellulose = make_Cellulose();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -12.0),
                                Eigen::Vector3d(2.5, 2.5, 4.0),
                                Eigen::Matrix3d::Identity());
    calc.set_source_material(&cellulose);
    calc.add_source_shield(&ss, 0.1);
    gate_a_two_stream(calc, 59.0, 0.5, false, "two-stream box 59 keV");
}

BOOST_AUTO_TEST_CASE(close_water_puck_122) {
    // Close geometry: near-degenerate per-position cones (direct stream
    // falls back to isotropic for some positions) + recoil electrons.
    Material nai = make_NaI();
    Material water = make_Water();
    Material pe = make_Polyethylene();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -1.6), 2.5, 0.5);
    calc.set_source_material(&water);
    calc.add_source_shield(&pe, 0.1);
    calc.enable_source_electron_transport(true);
    gate_a_two_stream(calc, 122.0, 0.5, false, "two-stream puck 122 keV");
}

BOOST_AUTO_TEST_CASE(absorbed_primary_electron_channel_662) {
    // Absorbed-primary escaped-electron channel: the primary Compton
    // recoil electron escapes the source toward the crystal, then the
    // degraded photon is absorbed in the source geometry with no surviving
    // secondary photon. Before June 2026 an early-out dropped these events
    // entirely. Assert the channel is populated, carries no FEP (electron
    // KE < E0 by construction), and that analog and two-stream agree on
    // its contribution within 5 sigma.
    Material nai = make_NaI();
    Material water = make_Water();
    Material pe = make_Polyethylene();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -1.6), 2.5, 0.5);
    calc.set_source_material(&water);
    calc.add_source_shield(&pe, 0.1);
    calc.enable_source_electron_transport(true);

    const double e = 662.0;
    const uint64_t n = 1500000;

    BiasingConfig analog;
    calc.set_biasing(analog);
    auto res_a = calc.compute(e, n);

    BiasingConfig ts;
    ts.two_stream = true;
    ts.two_stream_direct_fraction = 0.25;
    calc.set_biasing(ts);
    auto res_t = calc.compute(e, n);

    const auto& da = res_a.src_diag;
    const auto& dt = res_t.src_diag;
    BOOST_CHECK_GT(da.n_e_only, 0u);
    BOOST_CHECK_GT(dt.n_e_only, 0u);
    // No FEP from this channel here: a single recoil carries < 478 keV,
    // and multi-electron sums are capped at E0 minus the photon energy
    // remaining at absorption (>~5 keV) minus CSDA losses, i.e. < 657 keV
    // — kinematically outside the 662 +- 1.5 keV window. (In thick high-Z
    // shells multi-electron continuum coincidences CAN land in the FEP
    // window; this geometry cannot.)
    BOOST_CHECK_EQUAL(da.e_only_fep_w, 0.0);
    BOOST_CHECK_EQUAL(dt.e_only_fep_w, 0.0);

    const double Na = static_cast<double>(res_a.num_events_simulated);
    const double Nt = static_cast<double>(res_t.num_events_simulated);
    check_consistent(da.e_only_any_w / Na, std::sqrt(da.e_only_any_w2) / Na,
                     dt.e_only_any_w / Nt, std::sqrt(dt.e_only_any_w2) / Nt,
                     "e-only channel eps analog vs two-stream");
    check_consistent(res_a.total_efficiency, res_a.total_uncertainty,
                     res_t.total_efficiency, res_t.total_uncertainty,
                     "total with e-only channel analog vs two-stream");
}

BOOST_AUTO_TEST_CASE(off_axis_fe_662) {
    // Off-axis source: cone axis correctness in the direct stream.
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    const double d = 15.0, th = 45.0 * M_PI / 180.0;
    calc.set_point_source(
        Eigen::Vector3d(d * std::sin(th), 0.0, -d * std::cos(th)));
    calc.add_source_shield(&fe, 0.5);
    gate_a_two_stream(calc, 662.0, 0.5, false, "two-stream off-axis 662 keV");
}

BOOST_AUTO_TEST_CASE(spectrum_consistency_662) {
    // Gate-C style: the two-stream pulse-height spectrum must agree with
    // the analog spectrum per coarse bin. Conservative per-bin sigma:
    // analog sqrt(p(1-p)/N); two-stream sqrt(w_max * p / N) with
    // w_max = 1/(1-f) = 2 bounding the scatter-stream weight (direct
    // weights are far below 1).
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);

    const double e = 662.0;
    std::vector<float> edges;
    for (float b = 0.0f; b <= 672.0f; b += 42.0f) edges.push_back(b);

    const uint64_t n = 600000;
    BiasingConfig analog;
    calc.set_biasing(analog);
    auto res_a = calc.compute(e, n, 0, edges);

    BiasingConfig ts;
    ts.two_stream = true;
    ts.two_stream_direct_fraction = 0.5;
    calc.set_biasing(ts);
    auto res_t = calc.compute(e, n, 0, edges);

    const double Na = static_cast<double>(res_a.num_events_simulated);
    const double Nt = static_cast<double>(res_t.num_events_simulated);
    BOOST_REQUIRE_EQUAL(res_a.pulse_height_distribution.size(),
                        res_t.pulse_height_distribution.size());
    for (size_t i = 0; i < res_a.pulse_height_distribution.size(); ++i) {
        double pa = res_a.pulse_height_distribution[i];
        double pt = res_t.pulse_height_distribution[i];
        if (pa * Na < 100.0) continue;  // skip low-count bins
        double sig_a = std::sqrt(pa * (1.0 - pa) / Na);
        double sig_t = std::sqrt(2.0 * pt / Nt);
        double z = std::abs(pa - pt) / std::sqrt(sig_a * sig_a + sig_t * sig_t);
        BOOST_TEST_INFO("bin " << i << ": p_analog=" << pa
                               << " p_two_stream=" << pt << " z=" << z);
        BOOST_CHECK_LT(z, 5.0);
    }
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Compton-angle mixture biasing: unit checks + Gate A
// ============================================================

namespace {

/// Fine Simpson integration of the unnormalized Compton angular pdf,
/// reference for the composite Gauss-Legendre norm. Integrated in
/// s = sin(theta/2) (mu = 1 - 2s^2, dmu = -4s ds) so that S(x,Z) is
/// piecewise smooth in the integration variable (no sqrt cusp at mu = 1).
double compton_norm_simpson(double energy_keV, int Z, int n = 40000) {
    auto g = [&](double s) {
        return 4.0 * s *
               compton_angular_pdf_unnorm(energy_keV, Z, 1.0 - 2.0 * s * s);
    };
    double h = 1.0 / n;
    double sum = g(0.0) + g(1.0);
    for (int i = 1; i < n; ++i) {
        sum += g(i * h) * ((i % 2) ? 4.0 : 2.0);
    }
    return sum * h / 3.0;
}

/// Run two-stream + Compton-bias vs analog Gate A.
void gate_a_compton_bias(EfficiencyCalculator& calc, double energy_keV,
                         double gamma, const char* what) {
    BiasingConfig analog;
    calc.set_biasing(analog);
    auto res_a = calc.compute(precision_config(energy_keV));

    BiasingConfig b;
    b.two_stream = true;
    b.two_stream_direct_fraction = 0.25;
    b.compton_cone_fraction = gamma;
    calc.set_biasing(b);
    auto res_b = calc.compute(precision_config(energy_keV));

    check_consistent(res_a.full_energy_peak_efficiency,
                     res_a.fep_uncertainty,
                     res_b.full_energy_peak_efficiency,
                     res_b.fep_uncertainty, what);
    check_consistent(res_a.total_efficiency,
                     res_a.total_uncertainty,
                     res_b.total_efficiency,
                     res_b.total_uncertainty, what);
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(ComptonBiasUnit)

BOOST_AUTO_TEST_CASE(angular_norm_matches_fine_integration) {
    // Any norm error is a direct multiplicative bias on the IS weight, so
    // the composite-GL norm must match a fine reference integration.
    for (auto [Z, E] : {std::pair<int, double>{26, 100.0}, {26, 662.0},
                        {26, 2614.0}, {82, 300.0}, {82, 3000.0}, {8, 59.0},
                        {0, 662.0}}) {
        double gl = compton_angular_norm(E, Z);
        double ref = compton_norm_simpson(E, Z);
        BOOST_TEST_INFO("Z=" << Z << " E=" << E);
        BOOST_CHECK_CLOSE(gl, ref, 1e-3);  // 1e-5 relative
    }
}

BOOST_AUTO_TEST_CASE(angular_pdf_matches_sampler) {
    // The pdf used in the IS weight must be the density the analog sampler
    // actually draws from: bin-wise z-test of sampled mu against
    // pdf/norm integrated per bin (midpoint x bin width is adequate at 40
    // bins for a 5-sigma gate).
    const int Z = 26;
    const double E = 662.0;
    const int n = 400000, nbins = 40;
    std::mt19937_64 rng(777);
    std::vector<int> counts(nbins, 0);
    for (int i = 0; i < n; ++i) {
        double mu = sample_compton_cos_theta(E, Z, rng);
        int b = std::min(nbins - 1,
                         static_cast<int>((mu + 1.0) * 0.5 * nbins));
        counts[b]++;
    }
    const double norm = compton_angular_norm(E, Z);
    const double dx = 2.0 / nbins;
    for (int b = 0; b < nbins; ++b) {
        double mu_mid = -1.0 + (b + 0.5) * dx;
        double p = compton_angular_pdf_unnorm(E, Z, mu_mid) * dx / norm;
        double expect = n * p;
        if (expect < 50.0) continue;
        double sigma = std::sqrt(expect);
        double z = std::abs(counts[b] - expect) / sigma;
        BOOST_TEST_INFO("bin " << b << " mu=" << mu_mid << " counts="
                               << counts[b] << " expect=" << expect);
        BOOST_CHECK_LT(z, 5.0);
    }
}

BOOST_AUTO_TEST_CASE(finish_at_angle_matches_kinematics) {
    // finish_compton_at_angle at a given angle must reproduce the analog
    // Compton line (doppler off) and honor the direction override.
    std::mt19937_64 rng(5);
    const double E = 662.0;
    const double mu = 0.3;
    Eigen::Vector3d dir(0.0, 0.0, 1.0);
    Eigen::Vector3d w(0.6, 0.0, 0.8);  // |w|=1, w.z = 0.8 (any unit vector)
    auto r = finish_compton_at_angle(E, dir, mu, rng, 26, false, &w);
    BOOST_CHECK_CLOSE(r.scattered_energy_keV,
                      compton_scattered_energy(E, mu), 1e-9);
    BOOST_CHECK_SMALL((r.new_direction - w).norm(), 1e-12);
    BOOST_CHECK_CLOSE(r.scattered_energy_keV + r.deposited_energy_keV, E,
                      1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ComptonBiasGateA)

BOOST_AUTO_TEST_CASE(fe_source_shield_662) {
    Material nai = make_NaI();
    Material fe = make_Iron();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&fe, 0.5);
    gate_a_compton_bias(calc, 662.0, 0.3, "compton-bias cfg11 662 gamma=0.3");
    gate_a_compton_bias(calc, 662.0, 0.5, "compton-bias cfg11 662 gamma=0.5");
}

BOOST_AUTO_TEST_CASE(thick_pb_shield_662) {
    // The target regime: total efficiency dominated by shield-scattered
    // photons (2 cm Pb, scattered class ~47% of total).
    Material nai = make_NaI();
    Material pb = make_Lead();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    calc.add_source_shield(&pb, 2.0);
    calc.enable_source_electron_transport(true);
    gate_a_compton_bias(calc, 662.0, 0.5, "compton-bias 2cm Pb 662");
}

BOOST_AUTO_TEST_CASE(steel_box_cellulose_2614_pp) {
    // Above the PP threshold with an extended source: biased Compton
    // vertices coexist with annihilation/brems secondaries (which must
    // stay analog) and the whole-event weight composition.
    Material nai = make_NaI();
    Material ss = make_StainlessSteel304();
    Material cellulose = make_Cellulose();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -12.0),
                                Eigen::Vector3d(2.5, 2.5, 4.0),
                                Eigen::Matrix3d::Identity());
    calc.set_source_material(&cellulose);
    calc.add_source_shield(&ss, 0.2);
    calc.enable_source_electron_transport(true);
    gate_a_compton_bias(calc, 2614.0, 0.3, "compton-bias box 2614");
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Compton Doppler broadening: FEP/total efficiency invariance
// ============================================================

namespace {

/// Doppler ON vs OFF must give the same FEP and total efficiency: the
/// subshell binding energy is deposited locally, so the per-event total
/// deposit is identical, only the *spectral shape* near the Compton edge
/// changes.  5-sigma consistency (the combined ~1.4% statistical band is far
/// wider than any genuine sub-1% FEP shift, so this checks for accidental
/// efficiency bias, not the physical edge change).
void gate_a_doppler(EfficiencyCalculator& calc, double energy_keV,
                    const char* what) {
    BiasingConfig analog;
    calc.set_biasing(analog);

    calc.enable_doppler_broadening(false);
    auto res_off = calc.compute(precision_config(energy_keV));

    calc.enable_doppler_broadening(true);
    auto res_on = calc.compute(precision_config(energy_keV));

    check_consistent(res_off.full_energy_peak_efficiency,
                     res_off.fep_uncertainty,
                     res_on.full_energy_peak_efficiency,
                     res_on.fep_uncertainty, what);
    check_consistent(res_off.total_efficiency,
                     res_off.total_uncertainty,
                     res_on.total_efficiency,
                     res_on.total_uncertainty, what);
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(DopplerConsistencyGateA)

BOOST_AUTO_TEST_CASE(nai_bare_662) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    gate_a_doppler(calc, 662.0, "doppler config1 662 keV");
}

BOOST_AUTO_TEST_CASE(nai_bare_2614) {
    Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    gate_a_doppler(calc, 2614.0, "doppler config1 2614 keV");
}

BOOST_AUTO_TEST_SUITE_END()
