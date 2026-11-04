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

#define BOOST_TEST_MODULE EnvScatterTest
#include <boost/test/included/unit_test.hpp>

#include "geometry/EnvironmentGeometry.h"
#include "physics/BremsAngular.h"
#include "transport/AngularNorms.h"
#include "transport/ComptonScatter.h"
#include "transport/FluenceTally.h"
#include "transport/GaussLegendre.h"
#include "transport/NextEventEstimator.h"
#include "transport/RayleighAngular.h"
#include "io/EnvScatterKernel.h"
#include "materials/Material.h"

#include <cstdio>

#include <array>
#include <cmath>
#include <random>
#include <utility>
#include <vector>

using namespace ceelo;

namespace {
constexpr double kTwoPi = 6.283185307179586;
}

// ---------------------------------------------------------------------------
// Rayleigh angular evaluator (net-new: only a sampler existed before)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(rayleigh_pdf_is_positive_and_forward_peaked)
{
    // F(x,Z) falls with momentum transfer, so the density must be strongly
    // forward-peaked and monotonically decreasing away from mu = 1 (until the
    // (1+mu^2) polarization factor lifts the backscatter tail slightly).
    const double E = 100.0;
    const int Z = 8;  // oxygen, the dominant air/soil/concrete scatterer
    const double f_fwd = rayleigh_angular_pdf_unnorm(E, Z, 0.999);
    const double f_mid = rayleigh_angular_pdf_unnorm(E, Z, 0.0);
    const double f_bwd = rayleigh_angular_pdf_unnorm(E, Z, -1.0);
    BOOST_CHECK_GT(f_fwd, 0.0);
    BOOST_CHECK_GT(f_fwd, f_mid);
    BOOST_CHECK_GT(f_mid, f_bwd);
}

BOOST_AUTO_TEST_CASE(rayleigh_norm_matches_direct_quadrature)
{
    // The norm helper uses the s = sin(theta/2) substitution. Cross-check it
    // against a brute-force fine Simpson integration in mu, which needs a very
    // large number of points precisely BECAUSE of the cusp the substitution
    // removes -- so agreement here validates the substitution itself.
    for (double E : {60.0, 662.0}) {
        for (int Z : {7, 14}) {
            const int N = 200001;  // odd, for Simpson
            const double h = 2.0 / (N - 1);
            double sum = 0.0;
            for (int i = 0; i < N; ++i) {
                const double mu = -1.0 + i * h;
                const double w = (i == 0 || i == N - 1) ? 1.0
                               : ((i % 2) ? 4.0 : 2.0);
                sum += w * rayleigh_angular_pdf_unnorm(E, Z, mu);
            }
            const double brute = sum * h / 3.0;
            const double fast = rayleigh_angular_norm(E, Z);
            BOOST_CHECK_CLOSE(fast, brute, 0.5);  // 0.5% -- Simpson is the weak one
        }
    }
}

BOOST_AUTO_TEST_CASE(angular_densities_integrate_to_unity_on_the_sphere)
{
    // p(Omega) = pdf(mu) / (2*pi*norm) must integrate to 1 over 4*pi. This is
    // the identity the next-event estimator relies on: the cone integral of a
    // normalized density cannot exceed 1, which is what bounds the estimator.
    for (double E : {59.5, 662.0, 2614.0}) {
        for (int Z : {7, 8, 18}) {
            const double nc = compton_angular_norm(E, Z);
            const double integral_c = gl::integrate_mu_in_s(
                -1.0, 1.0, 16,
                [&](double mu) { return compton_angular_pdf_unnorm(E, Z, mu); });
            BOOST_CHECK_CLOSE(integral_c, nc, 1e-3);

            const double nr = rayleigh_angular_norm(E, Z);
            const double integral_r = gl::integrate_mu_in_s(
                -1.0, 1.0, 16,
                [&](double mu) { return rayleigh_angular_pdf_unnorm(E, Z, mu); });
            BOOST_CHECK_CLOSE(integral_r, nr, 1e-3);
        }
    }
}

// ---------------------------------------------------------------------------
// 2BS bremsstrahlung angular evaluator (physics/BremsAngular.{h,cpp})
//
// These three closure tests must hold BEFORE the evaluator is wired into the
// next-event estimator: the realized sampler density carries the envelope's
// 1/(1+y)^2 factor and the accept step's clamps, and an evaluator that drops
// either is evaluating a different distribution than the sampler draws.
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(brems_2bs_norm_matches_direct_quadrature)
{
    // The norm uses the envelope-map q-substitution, under which the
    // (1+y)^-2 factor cancels exactly against dy/dq. Cross-check against a
    // brute-force fine Simpson integration in mu (which has to resolve the
    // forward peak of width ~2/y_max the substitution removes), across the
    // kinematic range the estimator will use -- including y_max ~ 150 at
    // 2.5 MeV -- and the air / soil / iron effective Z values.
    for (double KE : {150.0, 500.0, 1000.0, 2500.0}) {
        for (double kf : {0.0, 0.3, 0.9}) {
            const double k = (kf == 0.0) ? 30.0 : kf * KE;
            for (double Z : {7.3, 11.6, 26.0}) {
                const double y_max = brems_2bs_ymax(KE);
                const int N = 200001;  // odd, for Simpson
                const double h = 2.0 / (N - 1);
                double sum = 0.0;
                for (int i = 0; i < N; ++i) {
                    const double mu = -1.0 + i * h;
                    const double w = (i == 0 || i == N - 1) ? 1.0
                                   : ((i % 2) ? 4.0 : 2.0);
                    const double y = y_max * (1.0 - mu) * 0.5;
                    sum += w * brems_2bs_pdf_y_unnorm(KE, k, Z, y);
                }
                const double brute = sum * h / 3.0 * (y_max * 0.5);
                const double fast = brems_2bs_norm(KE, k, Z);
                BOOST_CHECK_CLOSE(fast, brute, 0.5);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(brems_2bs_pdf_integrates_to_unity_on_sphere)
{
    // p(Omega) must integrate to 1 over 4*pi -- the identity that bounds the
    // estimator's brems emission score, exactly as for Compton/Rayleigh above.
    for (double KE : {150.0, 662.0, 2500.0}) {
        for (double kf : {0.1, 0.5, 0.9}) {
            const double k = kf * KE;
            for (double Z : {7.3, 11.6}) {
                const double integral = gl::integrate_mu_in_s(
                    -1.0, 1.0, 16,
                    [&](double mu) {
                        return kTwoPi * brems_2bs_pdf_sr(KE, k, Z, mu);
                    });
                BOOST_CHECK_CLOSE(integral, 1.0, 0.02);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(brems_2bs_sampler_matches_pdf_histogram)
{
    // Deterministic evaluator vs the actual rejection sampler, binned in
    // forward-weighted cos(theta) bands with binomial-sigma discipline (the
    // test class that exposes a mis-implemented rejection sampler). This
    // is the arbiter for the (1+y)^-2 envelope factor and the ceiling clamp:
    // an evaluator missing either fails these bands by far more than 5 sigma.
    struct Config { double KE, k, Z; };
    const Config configs[] = {
        {1000.0, 100.0, 11.6},   // MeV electron, soft photon, soil
        {300.0, 150.0, 7.3},     // mid-energy, hard photon, air
    };
    const double edges[] = {-1.0, 0.0, 0.9, 0.99, 0.999, 1.0};
    const int n_bands = 5;
    const int n_samples = 500000;

    for (const Config& c : configs) {
        // Expected band fractions from the evaluator.
        double expected[n_bands];
        double expected_sum = 0.0;
        for (int b = 0; b < n_bands; ++b) {
            expected[b] = gl::integrate_mu_in_s(
                edges[b], edges[b + 1], 8,
                [&](double mu) {
                    return kTwoPi * brems_2bs_pdf_sr(c.KE, c.k, c.Z, mu);
                });
            expected_sum += expected[b];
        }
        BOOST_CHECK_CLOSE(expected_sum, 1.0, 0.05);

        std::mt19937_64 rng(20260805);
        const Eigen::Vector3d e_dir = Eigen::Vector3d::UnitZ();
        int counts[n_bands] = {};
        for (int i = 0; i < n_samples; ++i) {
            const Eigen::Vector3d d = sample_brems_direction_2bs(
                e_dir, c.KE, c.k, c.Z, rng);
            const double mu = d.z();
            for (int b = 0; b < n_bands; ++b) {
                if (mu < edges[b + 1] || (b == n_bands - 1 && mu <= 1.0)) {
                    ++counts[b];
                    break;
                }
            }
        }
        for (int b = 0; b < n_bands; ++b) {
            const double obs = double(counts[b]) / n_samples;
            const double sigma = std::sqrt(
                std::max(expected[b] * (1.0 - expected[b]), 1e-12)
                / n_samples);
            BOOST_TEST_MESSAGE("KE=" << c.KE << " k=" << c.k
                               << " band [" << edges[b] << "," << edges[b + 1]
                               << "): obs " << obs << " exp " << expected[b]
                               << " sigma " << sigma);
            BOOST_CHECK_LT(std::abs(obs - expected[b]), 5.0 * sigma + 2e-4);
        }
    }
}

// ---------------------------------------------------------------------------
// Interpolated norm tables (the hot-path mitigation)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(interpolated_norms_match_exact_quadrature)
{
    // The estimator scores at every vertex of a multiple-scatter walk, so the
    // energy differs at essentially every call and the exact-key cache in
    // ComptonScatter never hits. These tables replace the quadrature; they must
    // agree with it to well inside the kernel's accuracy budget.
    const std::vector<int> zs{1, 7, 8, 14, 18, 20, 26};
    double worst_c = 0.0, worst_r = 0.0;
    for (int Z : zs) {
        for (double E = 31.0; E < 2900.0; E *= 1.37) {
            const double ec = compton_angular_norm(E, Z);
            const double ic = angular_norm_interp(ScatterProcess::Compton, E, Z);
            worst_c = std::max(worst_c, std::abs(ic - ec) / ec);

            const double er = rayleigh_angular_norm(E, Z);
            const double ir = angular_norm_interp(ScatterProcess::Rayleigh, E, Z);
            worst_r = std::max(worst_r, std::abs(ir - er) / er);
        }
    }
    BOOST_TEST_MESSAGE("worst relative norm interpolation error: Compton "
                       << worst_c << ", Rayleigh " << worst_r);
    BOOST_CHECK_LT(worst_c, 1e-5);
    // Rayleigh's floor is set by kinks inherited from the tabulated form factor,
    // not by grid curvature, so it does not respond to densification the way
    // Compton does (512 -> 1.2e-3, 2048 -> 2.2e-4). It does not need to: this
    // error propagates only into the Rayleigh part of the differential
    // scattering coefficient, which is at most ~half of the scatter at 60 keV in
    // air and far less elsewhere, so the worst-case effect on the scored fluence
    // is ~1e-4 -- roughly a thousand times below the ~10-20% accuracy this
    // kernel targets. Spending 4x more table-build time to chase it would buy
    // nothing physical.
    BOOST_CHECK_LT(worst_r, 5e-4);
}

// ---------------------------------------------------------------------------
// Analytic air transmission modes (the cheap, always-on path)
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(collimated_air_mode_removes_more_than_open_mode)
{
    // The two analytic modes differ only in whether Rayleigh counts as removal.
    // Open geometry: coherent scatter out of the line is compensated by scatter
    // in from the surrounding air, so it must NOT be removed -- on the air
    // benchmark an apparent 9% "attenuation gap" at 59 keV / 50 m is exactly
    // the Rayleigh mu/rho, a cone-biasing artifact. A collimated path blocks
    // the compensating in-scatter, so there Rayleigh IS removal.
    const Material air = make_Air();
    const double d = 5000.0;  // 50 m, where the difference is visible
    for (double E : {59.5, 662.0}) {
        const MacroscopicXS xs = air.macroscopic_xs(E * 1e-3);
        const double t_open = std::exp(-(xs.mu_pe + xs.mu_cs + xs.mu_pp) * d);
        const double t_coll = std::exp(-xs.mu_total() * d);
        BOOST_CHECK_LT(t_coll, t_open);  // collimated always removes more
        // At 59.5 keV Rayleigh is a large share of the air cross-section, so the
        // two modes must differ by several percent -- this is the regime where
        // choosing the wrong one shows up.
        if (E < 100.0) BOOST_CHECK_GT((t_open - t_coll) / t_open, 0.02);
    }
}

// ---------------------------------------------------------------------------
// EnvironmentGeometry tracing
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(open_air_traces_a_single_segment)
{
    const Material air = make_Air();
    EnvironmentGeometry env(air);
    std::vector<EnvSegment> segs;
    env.trace(Eigen::Vector3d(0, 0, -100), Eigen::Vector3d(0, 0, 1), segs, 100.0);
    BOOST_REQUIRE_EQUAL(segs.size(), 1u);
    BOOST_CHECK_EQUAL(segs[0].slab_index, -1);
    BOOST_CHECK_CLOSE(segs[0].t_end_cm - segs[0].t_start_cm, 100.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(optical_depth_in_air_matches_closed_form)
{
    const Material air = make_Air();
    EnvironmentGeometry env(air);
    const double d = 5000.0;  // 50 m
    const double tau = env.optical_depth(Eigen::Vector3d(0, 0, -d),
                                         Eigen::Vector3d(0, 0, 0), 0.662);
    BOOST_CHECK_CLOSE(tau, air.mu_total(0.662) * d, 1e-6);
}

BOOST_AUTO_TEST_CASE(reusable_trace_scratch_is_exactly_equivalent)
{
    const Material air = make_Air();
    const Material concrete = make_Concrete();
    const Eigen::Vector3d src(0.0, 0.0, -300.0);
    const EnvironmentGeometry env = EnvironmentGeometry::room(
        air, concrete, src, 100.0, 100.0, 300.0, 250.0,
        150.0, 150.0, 30.0);

    const Eigen::Vector3d pos(5.0, 7.0, -120.0);
    const Eigen::Vector3d dir = Eigen::Vector3d(0.31, -0.44, 0.83).normalized();
    std::vector<EnvSegment> baseline;
    std::vector<EnvSegment> reused;
    EnvTraceScratch scratch;
    env.trace(pos, dir, baseline, 900.0);
    env.trace(pos, dir, reused, scratch, 900.0);
    BOOST_REQUIRE_EQUAL(reused.size(), baseline.size());
    for (std::size_t i = 0; i < baseline.size(); ++i) {
        BOOST_CHECK_EQUAL(reused[i].material, baseline[i].material);
        BOOST_CHECK_EQUAL(reused[i].slab_index, baseline[i].slab_index);
        BOOST_CHECK_EQUAL(reused[i].t_start_cm, baseline[i].t_start_cm);
        BOOST_CHECK_EQUAL(reused[i].t_end_cm, baseline[i].t_end_cm);
    }

    const std::size_t first_capacity = scratch.crossings.capacity();
    env.trace(pos, dir, reused, scratch, 900.0);
    BOOST_CHECK_EQUAL(scratch.crossings.capacity(), first_capacity);
}

BOOST_AUTO_TEST_CASE(paired_detector_leg_attenuation_matches_two_traces)
{
    const Material air = make_Air();
    const Material concrete = make_Concrete();
    const Eigen::Vector3d src(0.0, 0.0, -300.0);
    const EnvironmentGeometry env = EnvironmentGeometry::room(
        air, concrete, src, 100.0, 100.0, 300.0, 250.0,
        150.0, 150.0, 30.0);
    const std::array<std::pair<Eigen::Vector3d, Eigen::Vector3d>, 4> chords{{
        {Eigen::Vector3d(0.0, 0.0, -300.0), Eigen::Vector3d::Zero()},
        {Eigen::Vector3d(0.0, -130.0, -20.0), Eigen::Vector3d::Zero()},
        {Eigen::Vector3d(260.0, 0.0, -40.0), Eigen::Vector3d::Zero()},
        {Eigen::Vector3d(-10.0, 180.0, 30.0), Eigen::Vector3d(5.0, 2.0, 0.0)}
    }};
    EnvTraceScratch scratch;
    for (const auto& chord : chords) {
        for (const auto& energies :
             {std::pair<double, double>{0.122, 0.122},
              std::pair<double, double>{0.180, 0.662},
              std::pair<double, double>{0.511, 2.614}}) {
            const std::array<double, 2> paired = env.transmission_pair(
                chord.first, chord.second, energies.first, energies.second,
                scratch);
            const double first = env.transmission(
                chord.first, chord.second, energies.first);
            const double second = env.transmission(
                chord.first, chord.second, energies.second);
            BOOST_CHECK_CLOSE(paired[0], first, 1e-12);
            BOOST_CHECK_CLOSE(paired[1], second, 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE(ground_slab_is_entered_and_exited)
{
    const Material air = make_Air();
    const Material soil = make_Soil();
    EnvironmentGeometry env(air);
    // Ground plane 100 cm below the origin, normal +y (up).
    env.add_slab(Eigen::Vector3d(0, -100, 0), Eigen::Vector3d(0, 1, 0), 50.0, soil);

    std::vector<EnvSegment> segs;
    // Fire straight down from the origin: 100 cm of air, then 50 cm of soil,
    // then out the back of the slab into void.
    env.trace(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, -1, 0), segs, 1000.0);
    BOOST_REQUIRE_EQUAL(segs.size(), 2u);
    BOOST_CHECK_EQUAL(segs[0].slab_index, -1);
    BOOST_CHECK_CLOSE(segs[0].t_end_cm, 100.0, 1e-6);
    BOOST_CHECK_EQUAL(segs[1].slab_index, 0);
    BOOST_CHECK_CLOSE(segs[1].t_end_cm - segs[1].t_start_cm, 50.0, 1e-6);

    // Optical depth of that chord must be the two-material sum.
    const double tau = env.optical_depth(Eigen::Vector3d(0, 0, 0),
                                         Eigen::Vector3d(0, -150, 0), 0.662);
    const double expect = air.mu_total(0.662) * 100.0 + soil.mu_total(0.662) * 50.0;
    BOOST_CHECK_CLOSE(tau, expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(field_builder_places_the_plane_at_both_heights)
{
    const Material air = make_Air();
    const Material soil = make_Soil();
    // Source 5 m away, source 30 cm off the ground, detector 100 cm up.
    const double d = 500.0, h_src = 30.0, h_det = 100.0;
    const Eigen::Vector3d src(0, 0, -d);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, soil, src, h_src, h_det);
    BOOST_REQUIRE_EQUAL(env.num_slabs(), 1u);

    const EnvSlab& s = env.slabs()[0];
    const double height_det = (Eigen::Vector3d::Zero() - s.point).dot(s.normal);
    const double height_src = (src - s.point).dot(s.normal);
    BOOST_CHECK_CLOSE(height_det, h_det, 1e-6);
    BOOST_CHECK_CLOSE(height_src, h_src, 1e-6);
}

BOOST_AUTO_TEST_CASE(room_builds_six_surfaces_at_the_right_distances)
{
    const Material air = make_Air();
    const Material concrete = make_Concrete();
    const Eigen::Vector3d src(0, 0, -300.0);
    const double h = 100.0, height = 300.0, half_w = 250.0;
    EnvironmentGeometry env = EnvironmentGeometry::room(
        air, concrete, src, h, h, height, half_w, 150.0, 150.0);

    BOOST_REQUIRE_EQUAL(env.num_slabs(), 6u);

    // Every surface must face the detector at the origin: the detector has to
    // be on the air side of all six planes, or the room is inside out.
    for (const EnvSlab& s : env.slabs()) {
        const double d_det = (Eigen::Vector3d::Zero() - s.point).dot(s.normal);
        BOOST_CHECK_GT(d_det, 0.0);
        const double d_src = (src - s.point).dot(s.normal);
        BOOST_CHECK_GT(d_src, 0.0);
    }
    // Floor and ceiling should sit at the requested heights from the detector.
    const double floor_h = (Eigen::Vector3d::Zero() - env.slabs()[0].point)
                               .dot(env.slabs()[0].normal);
    const double ceil_h = (Eigen::Vector3d::Zero() - env.slabs()[1].point)
                              .dot(env.slabs()[1].normal);
    BOOST_CHECK_CLOSE(floor_h, h, 1e-6);
    BOOST_CHECK_CLOSE(ceil_h, height - h, 1e-6);

    // The detector must classify as being in the air, not inside a wall.
    BOOST_CHECK_EQUAL(env.classify(Eigen::Vector3d::Zero()).slab_index, -1);
    BOOST_CHECK_EQUAL(env.classify(src).slab_index, -1);
}

BOOST_AUTO_TEST_CASE(field_builds_a_perpendicular_end_wall)
{
    // |h_src - h_det| == distance is the PERPENDICULAR-wall case (the plane
    // normal lies along the source->detector axis), not a degenerate one. A room
    // assembled by superposing per-surface kernels needs its end walls this way,
    // so the builder must accept them rather than silently return no slab.
    const Material air = make_Air();
    const Material concrete = make_Concrete();
    const double d = 200.0, back = 100.0;
    const Eigen::Vector3d src(0, 0, -d);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, concrete, src, back + d, back, 60.0);
    BOOST_REQUIRE_EQUAL(env.num_slabs(), 1u);

    const EnvSlab& s = env.slabs()[0];
    BOOST_CHECK_CLOSE((Eigen::Vector3d::Zero() - s.point).dot(s.normal), back, 1e-6);
    BOOST_CHECK_CLOSE((src - s.point).dot(s.normal), back + d, 1e-6);
    // The normal must lie along the axis (pointing from the wall toward the
    // source, i.e. -z here), since that is what "perpendicular wall" means.
    BOOST_CHECK_SMALL(std::abs(std::abs(s.normal.z()) - 1.0), 1e-9);
}

BOOST_AUTO_TEST_CASE(transverse_extent_bounds_the_air)
{
    const Material air = make_Air();
    const Eigen::Vector3d src(0, 0, -100), det(0, 0, 0);
    EnvironmentGeometry env =
        EnvironmentGeometry::open_air(air, src, det, /*extent*/ 50.0, /*margin*/ 10.0);

    // A ray fired sideways from the midpoint must leave the air after 50 cm.
    std::vector<EnvSegment> segs;
    env.trace(Eigen::Vector3d(0, 0, -50), Eigen::Vector3d(1, 0, 0), segs, 1000.0);
    BOOST_REQUIRE_EQUAL(segs.size(), 1u);
    BOOST_CHECK_CLOSE(segs[0].t_end_cm, 50.0, 1e-6);
}

// ---------------------------------------------------------------------------
// FluenceTally bookkeeping
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tally_correlated_scores_inflate_the_variance)
{
    // Two scores of 0.5 within ONE history must give the same mean as a single
    // score of 1.0, but the variance must reflect the history total (1.0), not
    // the individual scores. This is the per-history pre-summation contract.
    FluenceAxes ax;
    FluenceTally a(ax), b(ax);
    for (int i = 0; i < 100; ++i) {
        a.add(662.0, 0.5, 0.5, ChainClass::AirOnly);
        a.add(662.0, 0.5, 0.5, ChainClass::AirOnly);
        a.end_history();
        b.add(662.0, 0.5, 1.0, ChainClass::AirOnly);
        b.end_history();
    }
    const int ie = ax.energy_bin(662.0), imu = ax.mu_bin(0.5);
    BOOST_CHECK_CLOSE(a.mean(ie, imu), b.mean(ie, imu), 1e-9);
    BOOST_CHECK_CLOSE(a.sigma(ie, imu), b.sigma(ie, imu), 1e-9);
    // Zero variance: every history scored exactly the same total.
    BOOST_CHECK_SMALL(a.sigma(ie, imu), 1e-12);
}

BOOST_AUTO_TEST_CASE(tally_counts_empty_histories_in_the_denominator)
{
    FluenceAxes ax;
    FluenceTally t(ax);
    t.add(662.0, 0.0, 1.0, ChainClass::AirOnly);
    t.end_history();
    for (int i = 0; i < 9; ++i) t.end_history();  // nine histories scored nothing
    BOOST_CHECK_EQUAL(t.num_histories(), 10u);
    BOOST_CHECK_CLOSE(t.total_mean(), 0.1, 1e-9);
}

BOOST_AUTO_TEST_CASE(tally_reports_out_of_range_weight_instead_of_dropping_it)
{
    FluenceAxes ax;
    FluenceTally t(ax);
    t.add(5.0, 0.0, 2.5, ChainClass::AirOnly);  // below e_min
    t.end_history();
    BOOST_CHECK_CLOSE(t.dropped_weight(), 2.5, 1e-9);
    BOOST_CHECK_SMALL(t.total_mean(), 1e-12);
}

BOOST_AUTO_TEST_CASE(fit_components_keep_per_history_sufficient_statistics)
{
    FluenceAxes ax;
    ax.e_min_keV = 10.0;
    ax.e_max_keV = 1000.0;
    ax.n_energy = 12;
    ax.n_mu = 4;
    FluenceTally t(ax, {{662.0, LineOrigin::ElasticE0},
                        {510.998950, LineOrigin::Annihilation}});

    // Multiple scores into one component/cell must be summed BEFORE S2.
    t.add(100.0, 0.25, 0.25, ChainClass::AirOnly,
          FluxOrigin::Primary, PostCreationOrder::One);
    t.add(100.0, 0.25, 0.75, ChainClass::AirOnly,
          FluxOrigin::Primary, PostCreationOrder::One);
    t.add(200.0, -0.25, 2.0, ChainClass::OneBounce,
          FluxOrigin::Primary, PostCreationOrder::TwoPlus);
    t.add_line(1, 0.5, 3.0, ChainClass::OneBounce,
               FluxOrigin::Annihilation, PostCreationOrder::Zero);
    t.add_line(1, -0.5, 4.0, ChainClass::MultiSameSlab,
               FluxOrigin::Annihilation, PostCreationOrder::One);
    t.add(5.0, 0.0, 2.0, ChainClass::AirOnly);  // same-history drop
    t.add(5.0, 0.0, 3.0, ChainClass::AirOnly);
    t.end_history(7);
    t.end_history(7);  // an independent all-zero history

    BOOST_CHECK_EQUAL(t.num_histories(), 2u);
    BOOST_CHECK_EQUAL(t.replica_histories()[7], 2u);
    BOOST_CHECK_CLOSE(t.total_mean(), 5.0, 1e-9);  // (1+2+3+4)/2
    BOOST_CHECK_CLOSE(
        t.mean_component(ax.energy_bin(100.0), ax.mu_bin(0.25),
                         FluxOrigin::Primary, PostCreationOrder::One),
        0.5, 1e-9);
    BOOST_CHECK_CLOSE(
        t.sigma_component(ax.energy_bin(100.0), ax.mu_bin(0.25),
                          FluxOrigin::Primary, PostCreationOrder::One),
        0.5, 1e-9);
    BOOST_CHECK_CLOSE(
        t.line_mean_component(1, ax.mu_bin(0.5), FluxOrigin::Annihilation,
                              PostCreationOrder::Zero),
        1.5, 1e-9);

    const int c = flux_component_index(FluxOrigin::Primary,
                                       PostCreationOrder::One);
    const std::size_t component_cell = t.component_offsets()[c] +
        static_cast<std::size_t>(ax.index(ax.energy_bin(100.0),
                                          ax.mu_bin(0.25)));
    BOOST_CHECK_CLOSE(t.component_sum_w()[component_cell], 1.0, 1e-9);
    BOOST_CHECK_CLOSE(t.component_sum_w2()[component_cell], 1.0, 1e-9);
    BOOST_REQUIRE_EQUAL(t.replica_component_sum_w()[7].size(),
                        t.component_sum_w().size());
    BOOST_CHECK_CLOSE(t.replica_component_sum_w()[7][component_cell],
                      1.0f, 1e-5);
    const std::size_t total_cell = static_cast<std::size_t>(
        ax.index(ax.energy_bin(100.0), ax.mu_bin(0.25)));
    BOOST_REQUIRE_EQUAL(t.replica_continuum_sum_w()[7].size(),
                        t.continuum_sum_w().size());
    BOOST_CHECK_CLOSE(t.replica_continuum_sum_w()[7][total_cell],
                      1.0f, 1e-5);
    const std::size_t line_cell = static_cast<std::size_t>(
        1 * ax.n_mu + ax.mu_bin(0.5));
    BOOST_REQUIRE_EQUAL(t.replica_line_sum_w()[7].size(),
                        t.line_sum_w().size());
    BOOST_CHECK_CLOSE(t.replica_line_sum_w()[7][line_cell], 3.0f, 1e-5);

    const int class_cell =
        ax.index(ax.energy_bin(100.0), ax.mu_bin(0.25)) *
            static_cast<int>(ChainClass::Count) +
        static_cast<int>(ChainClass::AirOnly);
    BOOST_CHECK_CLOSE(t.class_sum_w()[class_cell], 1.0, 1e-9);
    BOOST_CHECK_CLOSE(t.class_sum_w2()[class_cell], 1.0, 1e-9);
    BOOST_CHECK_CLOSE(t.dropped_weight(), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(t.dropped_weight_sum_w2(), 25.0, 1e-9);
    BOOST_CHECK_CLOSE(t.max_history_total(), 10.0, 1e-9);

    double component_total = 0.0;
    for (int origin = 0; origin < static_cast<int>(FluxOrigin::Count); ++origin)
        for (int order = 0;
             order < static_cast<int>(PostCreationOrder::Count); ++order)
            component_total += t.total_mean_component(
                static_cast<FluxOrigin>(origin),
                static_cast<PostCreationOrder>(order));
    BOOST_CHECK_CLOSE(component_total, t.total_mean(), 1e-9);
}

BOOST_AUTO_TEST_CASE(fit_component_axes_are_native_and_merge_checked)
{
    const FluenceAxes common =
        make_registered_landmark_axes(2614.0, 12, 4, 10.0);
    std::vector<FluenceAxes> native =
        make_fit_component_axes(2614.0, 12, 4, 10.0);
    native[static_cast<int>(FluxOrigin::Annihilation)] =
        make_registered_landmark_axes(510.998950, 8, 4, 10.0);
    native[static_cast<int>(FluxOrigin::Brems)].n_energy = 10;

    FluenceTally a(common, env_line_channels(2614.0), native);
    FluenceTally b(common, env_line_channels(2614.0), native);
    a.add(100.0, 0.0, 1.0, ChainClass::AirOnly,
          FluxOrigin::Primary, PostCreationOrder::One);
    a.add(300.0, 0.0, 2.0, ChainClass::AirOnly,
          FluxOrigin::Annihilation, PostCreationOrder::One);
    a.add(50.0, 0.0, 3.0, ChainClass::AirOnly,
          FluxOrigin::Brems, PostCreationOrder::Zero);
    a.end_history(2);
    b.end_history(3);
    BOOST_CHECK_NO_THROW(a.merge(b));
    BOOST_CHECK_EQUAL(a.replica_histories()[2], 1u);
    BOOST_CHECK_EQUAL(a.replica_histories()[3], 1u);

    const FluenceAxes& ann = a.component_axes(FluxOrigin::Annihilation);
    const FluenceAxes& brems = a.component_axes(FluxOrigin::Brems);
    BOOST_CHECK_CLOSE(ann.e_max_keV, 510.998950, 1e-9);
    BOOST_CHECK_EQUAL(ann.n_energy, 8);
    BOOST_CHECK(!brems.has_edges());
    BOOST_CHECK_EQUAL(brems.n_energy, 10);
    BOOST_CHECK_CLOSE(
        a.mean_component(ann.energy_bin(300.0), ann.mu_bin(0.0),
                         FluxOrigin::Annihilation, PostCreationOrder::One),
        1.0, 1e-9);
    for (double w : a.component_dropped_weight())
        BOOST_CHECK_SMALL(w, 1e-12);

    std::vector<FluenceAxes> wrong = native;
    wrong[static_cast<int>(FluxOrigin::Brems)].n_energy = 11;
    FluenceTally mismatched(common, env_line_channels(2614.0), wrong);
    BOOST_CHECK_THROW(a.merge(mismatched), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(coarse_stopping_vector_has_exact_outer_moments)
{
    FluenceAxes ax;
    ax.e_min_keV = 10.0;
    ax.e_max_keV = 662.0;
    ax.n_energy = 12;
    ax.n_mu = 4;
    FluenceTally t(ax, {{662.0, LineOrigin::ElasticE0}});
    t.add_line(0, 0.5, 2.0, ChainClass::AirOnly,
               FluxOrigin::Primary, PostCreationOrder::One);
    t.end_history(4);
    t.end_history(4);  // zero history contributes to N, not S1/S2

    const auto& s1 = t.coarse_sum_w();
    const auto& s2 = t.coarse_outer_sum_w();
    BOOST_CHECK_CLOSE(s1[11], 2.0, 1e-9);       // exact E0 -> final band
    BOOST_CHECK_CLOSE(s1[12], 1.0, 1e-9);       // 2 P1(0.5)
    BOOST_CHECK_CLOSE(s1[13], -0.25, 1e-9);     // 2 P2(0.5)
    BOOST_CHECK_CLOSE(s1[14], -0.875, 1e-9);    // 2 P3(0.5)
    BOOST_CHECK_CLOSE(s1[15], -0.578125, 1e-9); // 2 P4(0.5)
    BOOST_CHECK_CLOSE(s2[11 * FluenceTally::kCoarseDimension + 11],
                      4.0, 1e-9);
    BOOST_CHECK_CLOSE(s2[11 * FluenceTally::kCoarseDimension + 12],
                      2.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(default_replica_assignment_is_balanced)
{
    FluenceTally t;
    for (int i = 0; i < FluenceTally::kNumLogicalReplicas; ++i)
        t.end_history();
    for (uint64_t n : t.replica_histories()) BOOST_CHECK_EQUAL(n, 1u);
    BOOST_CHECK_THROW(t.end_history(FluenceTally::kNumLogicalReplicas),
                      std::out_of_range);
}

BOOST_AUTO_TEST_CASE(transport_diagnostics_count_analog_events_and_vertex_caps)
{
    const Material air = make_Air();
    const EnvironmentGeometry env(air);
    EnvScatterConfig cfg;
    cfg.axes = make_registered_landmark_axes(122.0, 12, 4, 10.0);
    cfg.e_floor_keV = 10.0;
    cfg.max_vertices = 1;
    cfg.include_annihilation = false;
    cfg.include_brems = false;
    const Eigen::Vector3d source(0.0, 0.0, -100.0);
    const NextEventEstimator estimator(env, cfg, source);
    EnvScatterResult result(cfg.axes, 122.0);
    constexpr std::uint64_t histories = 128;
    estimator.run(122.0, histories, 0x8a55ULL, result);

    const auto& counts = result.diagnostics.interaction_counts;
    const std::uint64_t photon_interactions =
        counts[static_cast<std::size_t>(EnvTransportCounter::Photoelectric)] +
        counts[static_cast<std::size_t>(EnvTransportCounter::Compton)] +
        counts[static_cast<std::size_t>(EnvTransportCounter::Rayleigh)] +
        counts[static_cast<std::size_t>(EnvTransportCounter::PairProduction)];
    BOOST_CHECK_EQUAL(photon_interactions, histories);
    BOOST_CHECK_EQUAL(
        result.diagnostics.max_vertex_terminations +
            counts[static_cast<std::size_t>(EnvTransportCounter::Photoelectric)],
        histories);
    BOOST_CHECK_EQUAL(
        counts[static_cast<std::size_t>(EnvTransportCounter::BremsEmission)],
        0U);
}

// ---------------------------------------------------------------------------
// Surface flux projection (the deliverable) and kernel table I/O
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(surface_flux_projects_onto_the_right_faces)
{
    // mu_in is the cosine between the photon's direction of travel and +z.
    // Head-on (mu = +1) must land entirely on the front face; mu = -1 entirely
    // on the back; mu = 0 entirely on the side. Getting this convention backwards
    // would silently swap face and back in every result.
    DetectorEnvelope det;  // 3"x3"
    const double face_area = M_PI * det.radius_cm * det.radius_cm;
    const double side_area = 2.0 * det.radius_cm * det.length_cm;

    const SurfaceFlux head_on = direct_surface_flux(1.0, +1.0, det);
    BOOST_CHECK_CLOSE(head_on.face, face_area, 1e-9);
    BOOST_CHECK_SMALL(head_on.back, 1e-12);
    BOOST_CHECK_SMALL(head_on.side, 1e-12);

    const SurfaceFlux from_behind = direct_surface_flux(1.0, -1.0, det);
    BOOST_CHECK_CLOSE(from_behind.back, face_area, 1e-9);
    BOOST_CHECK_SMALL(from_behind.face, 1e-12);

    const SurfaceFlux broadside = direct_surface_flux(1.0, 0.0, det);
    BOOST_CHECK_CLOSE(broadside.side, side_area, 1e-9);
    BOOST_CHECK_SMALL(broadside.face, 1e-12);
    BOOST_CHECK_SMALL(broadside.back, 1e-12);
}

BOOST_AUTO_TEST_CASE(kernel_table_round_trips)
{
    EnvKernelTable t;
    t.kind = "soil";
    t.air_name = "Air";
    t.air_density = 0.001205;
    t.slab_name = "Soil";
    t.slab_density = 1.5;
    t.scoring_radius_cm = 4.0;
    t.e0_grid_keV = {59.5, 662.0};
    t.distance_grid_cm = {100.0, 500.0};
    t.h_src_grid_cm = {30.0};
    t.h_det_grid_cm = {100.0};
    t.provenance = "unit test";

    EnvKernelNode n;
    n.e0_keV = 662.0;
    n.distance_cm = 500.0;
    n.h_src_cm = 30.0;
    n.h_det_cm = 100.0;
    n.direct_transmission = 0.9547;
    n.direct_fluence = 3.04e-7;
    n.histories = 123456;
    n.max_rel_sigma = 0.012;
    n.multi_slab_fraction = 0.34;
    n.s_hat.assign(t.axes.n_bins(), 0.0);
    n.s_hat_slab.assign(t.axes.n_bins(), 0.0);
    n.s_hat[5] = 1.25e-3;
    n.s_hat_slab[5] = 1.05e-3;
    EnvKernelNode::Line L;
    L.energy_keV = 511.0;
    L.origin = 1;
    L.s_hat.assign(t.axes.n_mu, 0.0);
    L.s_hat_slab.assign(t.axes.n_mu, 0.0);
    L.s_hat[3] = 7.5e-4;
    L.s_hat_slab[3] = 6.5e-4;
    n.lines.push_back(L);
    t.nodes.push_back(n);

    const std::string path = "/tmp/ceelo_env_kernel_roundtrip.txt";
    BOOST_REQUIRE(save_env_kernel(t, path));

    EnvKernelTable back;
    BOOST_REQUIRE(load_env_kernel(back, path));
    BOOST_CHECK_EQUAL(back.kind, t.kind);
    BOOST_REQUIRE_EQUAL(back.nodes.size(), 1u);
    BOOST_CHECK_CLOSE(back.nodes[0].e0_keV, 662.0, 1e-9);
    BOOST_CHECK_CLOSE(back.nodes[0].direct_transmission, 0.9547, 1e-6);
    BOOST_CHECK_EQUAL(back.nodes[0].histories, 123456u);
    BOOST_CHECK_CLOSE(back.nodes[0].multi_slab_fraction, 0.34, 1e-6);
    BOOST_REQUIRE_EQUAL(back.nodes[0].s_hat.size(), n.s_hat.size());
    BOOST_CHECK_CLOSE(back.nodes[0].s_hat[5], 1.25e-3, 1e-6);
    BOOST_CHECK_CLOSE(back.nodes[0].s_hat_slab[5], 1.05e-3, 1e-6);
    BOOST_REQUIRE_EQUAL(back.nodes[0].lines.size(), 1u);
    BOOST_CHECK_CLOSE(back.nodes[0].lines[0].energy_keV, 511.0, 1e-9);
    BOOST_CHECK_EQUAL(int(back.nodes[0].lines[0].origin), 1);
    BOOST_CHECK_CLOSE(back.nodes[0].lines[0].s_hat[3], 7.5e-4, 1e-6);
    BOOST_CHECK_CLOSE(back.nodes[0].lines[0].s_hat_slab[3], 6.5e-4, 1e-6);
    BOOST_CHECK_EQUAL(back.e0_grid_keV.size(), 2u);
    std::remove(path.c_str());

    // Binary float32 form: same content round-trips (to float32 precision on
    // payload, exact on metadata), the truncated above-E0 tail reconstitutes
    // as exact zeros, and the generic load path auto-detects the magic.
    // s_hat[5] sits in energy bin 0 (index = ie * n_mu + imu), far below the
    // node's E0 bin, so truncation keeps it while the file drops the empty
    // above-E0 rows.
    const std::string bpath = "/tmp/ceelo_env_kernel_roundtrip.bin";
    BOOST_REQUIRE(save_env_kernel_binary(t, bpath));

    EnvKernelTable bin;
    BOOST_REQUIRE(load_env_kernel(bin, bpath));   // auto-detect, not _binary
    BOOST_CHECK_EQUAL(bin.kind, t.kind);
    BOOST_CHECK_EQUAL(bin.provenance, t.provenance);
    BOOST_REQUIRE_EQUAL(bin.nodes.size(), 1u);
    BOOST_CHECK_CLOSE(bin.nodes[0].direct_transmission, 0.9547, 1e-9);
    BOOST_CHECK_EQUAL(bin.nodes[0].histories, 123456u);
    BOOST_REQUIRE_EQUAL(bin.nodes[0].s_hat.size(), n.s_hat.size());
    BOOST_CHECK_CLOSE(bin.nodes[0].s_hat[5], 1.25e-3, 1e-4);      // float32
    BOOST_CHECK_CLOSE(bin.nodes[0].s_hat_slab[5], 1.05e-3, 1e-4);
    BOOST_CHECK_EQUAL(bin.nodes[0].s_hat.back(), 0.0);            // rebuilt tail
    BOOST_REQUIRE_EQUAL(bin.nodes[0].lines.size(), 1u);
    BOOST_CHECK_CLOSE(bin.nodes[0].lines[0].s_hat[3], 7.5e-4, 1e-4);
    BOOST_CHECK_EQUAL(int(bin.nodes[0].lines[0].origin), 1);

    // The binary file must actually be smaller than the text one would be for
    // a full node -- and the E0 truncation must have dropped the empty tail:
    // for E0 = 662 on the default 30-3000/48 axis, bins above ~697 keV go.
    {
        std::ifstream fb(bpath, std::ios::binary | std::ios::ate);
        const auto bin_size = fb.tellg();
        const int ie0 = t.axes.energy_bin(662.0);
        const std::size_t kept = static_cast<std::size_t>(ie0 + 1) * t.axes.n_mu;
        BOOST_CHECK_LT(kept, n.s_hat.size());
        // 2 payload arrays * kept float32 + small fixed overhead bounds it.
        BOOST_CHECK_LT(static_cast<std::size_t>(bin_size),
                       2u * kept * 4u + 1024u);
    }
    std::remove(bpath.c_str());

    // A node whose above-E0 tail is NOT zero (a broken invariant) must fall
    // back to full storage rather than silently truncate.
    EnvKernelTable t2 = t;
    t2.nodes[0].s_hat.back() = 1e-9;
    BOOST_REQUIRE(save_env_kernel_binary(t2, bpath));
    EnvKernelTable bin2;
    BOOST_REQUIRE(load_env_kernel_binary(bin2, bpath));
    BOOST_CHECK_CLOSE(bin2.nodes[0].s_hat.back(), 1e-9, 1e-4);
    std::remove(bpath.c_str());
}

BOOST_AUTO_TEST_CASE(kernel_shards_merge_and_report_completeness)
{
    // A distributed bake returns one file per shard; reassembly must refuse
    // anything that would silently corrupt a write-once product.
    EnvKernelTable base;
    base.kind = "soil";
    base.air_name = "Air";
    base.air_density = 0.001205;
    base.slab_name = "Soil";
    base.slab_density = 1.5;
    base.scoring_radius_cm = 15.0;
    base.e0_grid_keV = {662.0};
    base.distance_grid_cm = {100.0, 500.0};
    base.h_src_grid_cm = {100.0};
    base.h_det_grid_cm = {100.0};

    auto make_node = [&](double d) {
        EnvKernelNode n;
        n.e0_keV = 662.0;
        n.distance_cm = d;
        n.h_src_cm = 100.0;
        n.h_det_cm = 100.0;
        n.energy_edges = make_landmark_axes(662.0, 96, base.axes.n_mu).edges;
        n.s_hat.assign((n.energy_edges.size() - 1) * base.axes.n_mu, 0.0);
        n.s_hat[7] = 1e-3;
        return n;
    };

    EnvKernelTable sh_a = base, sh_b = base;
    sh_a.nodes.push_back(make_node(100.0));
    sh_b.nodes.push_back(make_node(500.0));

    EnvKernelTable merged;
    BOOST_REQUIRE(merge_env_kernel_shard(merged, sh_a));
    BOOST_REQUIRE(merge_env_kernel_shard(merged, sh_b));
    BOOST_CHECK_EQUAL(merged.nodes.size(), 2u);

    std::string missing;
    BOOST_CHECK(env_kernel_is_complete(merged, &missing));
    BOOST_CHECK(missing.empty());

    // A duplicate node means overlapping shard assignment -- some other node
    // is therefore missing, so averaging it away would hide a hole.
    BOOST_CHECK(!merge_env_kernel_shard(merged, sh_b));

    // A shard from a different grid is refused outright.
    EnvKernelTable wrong = base;
    wrong.distance_grid_cm = {100.0, 700.0};
    wrong.nodes.push_back(make_node(700.0));
    BOOST_CHECK(!merge_env_kernel_shard(merged, wrong));

    // An incomplete table names what is absent.
    EnvKernelTable partial;
    BOOST_REQUIRE(merge_env_kernel_shard(partial, sh_a));
    BOOST_CHECK(!env_kernel_is_complete(partial, &missing));
    BOOST_CHECK(!missing.empty());

    // Round-trip through the binary container preserves per-node edges.
    const std::string p = "/tmp/ceelo_env_shard.bin";
    BOOST_REQUIRE(save_env_kernel_binary(merged, p));
    EnvKernelTable back;
    BOOST_REQUIRE(load_env_kernel(back, p));
    BOOST_REQUIRE_EQUAL(back.nodes.size(), 2u);
    BOOST_REQUIRE_EQUAL(back.nodes[0].energy_edges.size(),
                        merged.nodes[0].energy_edges.size());
    BOOST_CHECK_CLOSE(back.nodes[0].energy_edges.back(), 662.0, 1e-9);
    BOOST_CHECK_CLOSE(back.nodes[0].s_hat[7], 1e-3, 1e-4);
    std::remove(p.c_str());
}

// ---------------------------------------------------------------------------
// Landmark-blocked energy axis
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(landmark_axes_pin_the_physical_edges)
{
    // The whole point is that no bin STRADDLES a spectral landmark: the
    // backscatter edge, 511 keV, and E0 must each be an exact bin boundary,
    // and the axis must stop at E0 (nothing scatters up).
    for (double e0 : {122.0, 662.0, 1332.0, 2614.0}) {
        const FluenceAxes ax = make_landmark_axes(e0, 192, 32, 30.0);
        BOOST_REQUIRE(ax.has_edges());
        BOOST_REQUIRE_EQUAL(static_cast<int>(ax.edges.size()), ax.n_energy + 1);

        // Ascending, spanning [30, E0] exactly.
        for (int i = 0; i + 1 < static_cast<int>(ax.edges.size()); ++i)
            BOOST_CHECK_LT(ax.edges[i], ax.edges[i + 1]);
        BOOST_CHECK_CLOSE(ax.edges.front(), 30.0, 1e-9);
        BOOST_CHECK_CLOSE(ax.edges.back(), e0, 1e-9);

        auto is_edge = [&](double e) {
            for (double x : ax.edges)
                if (std::abs(x - e) < 1e-9 * std::max(1.0, e)) return true;
            return false;
        };
        const double e_bs = e0 / (1.0 + 2.0 * e0 / 510.998950);
        BOOST_CHECK_MESSAGE(is_edge(e_bs),
                            "backscatter edge " << e_bs << " not pinned for E0=" << e0);
        if (e0 > 511.0 * 1.05) BOOST_CHECK(is_edge(510.998950));

        // Binning is consistent with the edges, and E0 is in range while
        // anything above it is dropped.
        BOOST_CHECK_EQUAL(ax.energy_bin(e0 * 1.001), -1);
        BOOST_CHECK_EQUAL(ax.energy_bin(29.0), -1);
        const int ib = ax.energy_bin(0.5 * (ax.edges[3] + ax.edges[4]));
        BOOST_CHECK_EQUAL(ib, 3);
        // A photon exactly AT the backscatter edge lands in the block above
        // it, and the bin below it ends there -- i.e. the edge separates.
        const int i_bs = ax.energy_bin(e_bs);
        BOOST_CHECK_CLOSE(ax.energy_low(i_bs), e_bs, 1e-6);

        // Bin count is close to the requested target.
        BOOST_CHECK_GT(ax.n_energy, 150);
        BOOST_CHECK_LT(ax.n_energy, 240);
    }
}

BOOST_AUTO_TEST_CASE(registered_landmark_axes_put_bin_k_at_the_same_place)
{
    // THE PROPERTY THAT MAKES E0 AN INTERPOLABLE AXIS. Every node gets the
    // same number of channels spanning its OWN [e_cut, E0], split equally
    // between the two blocks, so bin k sits at an identical fractional
    // position u for every E0 and the backscatter edge is always the middle
    // bin. Interpolating bin k against bin k across E0 is then automatically
    // feature-registered -- no coordinate transform at query time.
    //
    // Contrast make_landmark_axes, which pins the same EDGES but allocates
    // bins by log-width, so bin k drifts with E0 (180/12 vs 139/53 at 192
    // bins) -- fine within one node, wrong for interpolating across nodes.
    const std::vector<double> e0s{59.0, 122.0, 662.0, 1332.0, 2614.0};
    std::vector<double> u_at_bin20, u_at_bin150;

    for (double e0 : e0s) {
        const FluenceAxes ax = make_registered_landmark_axes(e0, 192, 32, 10.0);
        BOOST_REQUIRE(ax.has_edges());
        BOOST_CHECK_EQUAL(ax.n_energy, 192);
        BOOST_CHECK_CLOSE(ax.edges.front(), 10.0, 1e-9);
        BOOST_CHECK_CLOSE(ax.edges.back(), e0, 1e-9);

        // The backscatter edge is EXACTLY the middle edge, for every E0.
        const double e_bs = e0 / (1.0 + 2.0 * e0 / 510.998950);
        BOOST_CHECK_CLOSE(ax.edges[96], e_bs, 1e-9);

        auto u_of_bin = [&](int k) {
            const double em = std::sqrt(ax.energy_low(k) * ax.energy_high(k));
            return (em <= e_bs)
                ? 0.5 * std::log(em / 10.0) / std::log(e_bs / 10.0)
                : 0.5 + 0.5 * std::log(em / e_bs) / std::log(e0 / e_bs);
        };
        u_at_bin20.push_back(u_of_bin(20));
        u_at_bin150.push_back(u_of_bin(150));
    }

    // Same bin index -> same u at every E0 (that is the whole point).
    for (std::size_t i = 1; i < e0s.size(); ++i) {
        BOOST_CHECK_CLOSE(u_at_bin20[i], u_at_bin20[0], 1e-6);
        BOOST_CHECK_CLOSE(u_at_bin150[i], u_at_bin150[0], 1e-6);
    }

    // A 59 keV source gets as many channels as a 2614 keV one -- which is NOT
    // true on the fixed absolute grid (~15 vs ~47 over [10 keV, E0]).
    const FluenceAxes lo = make_registered_landmark_axes(59.0, 192, 32, 10.0);
    const FluenceAxes hi = make_registered_landmark_axes(2614.0, 192, 32, 10.0);
    BOOST_CHECK_EQUAL(lo.n_energy, hi.n_energy);
    FluenceAxes fixed;
    fixed.e_min_keV = 10.0;
    fixed.e_max_keV = 3000.0;
    fixed.n_energy = 48;
    const int n_lo = fixed.energy_bin(59.0) + 1;
    const int n_hi = fixed.energy_bin(2614.0) + 1;
    BOOST_CHECK_LT(n_lo * 2, n_hi);   // the imbalance the registered axis removes
}

BOOST_AUTO_TEST_CASE(landmark_axes_tally_and_merge_consistently)
{
    // A tally on a landmark axis must score, sum and merge exactly like a
    // log-uniform one -- and merging ACROSS different axes must throw rather
    // than silently produce a mixed-grid table (the shardable-bake hazard).
    const FluenceAxes ax = make_landmark_axes(662.0, 96, 16, 30.0);
    FluenceTally a(ax, {{662.0, LineOrigin::ElasticE0}});
    FluenceTally b(ax, {{662.0, LineOrigin::ElasticE0}});

    a.add(100.0, 0.5, 2.0, ChainClass::OneBounce);
    a.add(400.0, -0.5, 3.0, ChainClass::AirOnly, FluxOrigin::Brems);
    a.end_history();
    b.add(100.0, 0.5, 4.0, ChainClass::OneBounce);
    b.end_history();

    a.merge(b);
    BOOST_CHECK_CLOSE(a.total_mean(), (2.0 + 3.0 + 4.0) / 2.0, 1e-9);
    BOOST_CHECK_CLOSE(a.total_mean_origin(FluxOrigin::Brems), 1.5, 1e-9);

    // Out-of-range (above E0) is dropped, not silently binned.
    FluenceTally c(ax, {});
    c.add(700.0, 0.0, 1.0, ChainClass::AirOnly);
    c.end_history();
    BOOST_CHECK_EQUAL(c.total_mean(), 0.0);
    BOOST_CHECK_CLOSE(c.dropped_weight(), 1.0, 1e-9);

    FluenceTally other(make_landmark_axes(1332.0, 96, 16, 30.0), {});
    BOOST_CHECK_THROW(other.merge(c), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Discrete line channels
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(line_channels_route_without_changing_totals)
{
    // A line score and a grid score must both land in total_mean(): "total
    // scattered fluence" has to be independent of how the flux is split between
    // lines and continuum, or every existing comparison against it would shift
    // when the split was introduced.
    FluenceAxes ax;
    std::vector<LineDef> lines{{662.0, LineOrigin::ElasticE0},
                               {510.998950, LineOrigin::Annihilation}};
    FluenceTally t(ax, lines);
    BOOST_REQUIRE_EQUAL(t.n_lines(), 2);

    t.add(300.0, 0.5, 2.0, ChainClass::AirOnly);
    t.add_line(0, 0.5, 3.0, ChainClass::OneBounce);
    t.add_line(1, -0.5, 5.0, ChainClass::OneBounce);
    t.end_history();

    BOOST_CHECK_CLOSE(t.total_mean(), 10.0, 1e-9);           // 2 + 3 + 5
    BOOST_CHECK_CLOSE(t.line_total_mean(0), 3.0, 1e-9);
    BOOST_CHECK_CLOSE(t.line_total_mean(1), 5.0, 1e-9);
    // Line scores must NOT leak into the continuum grid.
    double grid = 0.0;
    for (int ie = 0; ie < ax.n_energy; ++ie)
        for (int j = 0; j < ax.n_mu; ++j) grid += t.mean(ie, j);
    BOOST_CHECK_CLOSE(grid, 2.0, 1e-9);
    // Angle and class resolution survive.
    BOOST_CHECK_CLOSE(t.line_mean(1, ax.mu_bin(-0.5)), 5.0, 1e-9);
    BOOST_CHECK_CLOSE(t.line_mean_class(1, ax.mu_bin(-0.5),
                                        ChainClass::OneBounce), 5.0, 1e-9);
    BOOST_CHECK_SMALL(t.line_mean_class(1, ax.mu_bin(-0.5),
                                        ChainClass::AirOnly), 1e-12);
}

BOOST_AUTO_TEST_CASE(line_channel_variance_uses_per_history_totals)
{
    // Same contract as the grid: two correlated scores in ONE history must give
    // the same mean AND the same sigma as one score of their sum.
    FluenceAxes ax;
    std::vector<LineDef> lines{{662.0, LineOrigin::ElasticE0}};
    FluenceTally a(ax, lines), b(ax, lines);
    for (int i = 0; i < 100; ++i) {
        a.add_line(0, 0.5, 0.5, ChainClass::AirOnly);
        a.add_line(0, 0.5, 0.5, ChainClass::AirOnly);
        a.end_history();
        b.add_line(0, 0.5, 1.0, ChainClass::AirOnly);
        b.end_history();
    }
    BOOST_CHECK_CLOSE(a.line_total_mean(0), b.line_total_mean(0), 1e-9);
    BOOST_CHECK_CLOSE(a.line_total_sigma(0), b.line_total_sigma(0), 1e-9);
    BOOST_CHECK_SMALL(a.line_total_sigma(0), 1e-12);
}

BOOST_AUTO_TEST_CASE(line_channel_merge_rejects_mismatched_sets)
{
    FluenceAxes ax;
    FluenceTally a(ax, {{662.0, LineOrigin::ElasticE0}});
    FluenceTally b(ax, {{1332.0, LineOrigin::ElasticE0}});
    FluenceTally c(ax, {{662.0, LineOrigin::ElasticE0}});
    c.add_line(0, 0.0, 1.0, ChainClass::AirOnly);
    c.end_history();
    BOOST_CHECK_NO_THROW(a.merge(c));
    BOOST_CHECK_THROW(a.merge(b), std::invalid_argument);   // different energy
}

BOOST_AUTO_TEST_CASE(annihilation_line_exists_only_above_pair_threshold)
{
    BOOST_CHECK_EQUAL(env_line_channels(662.0).size(), 1u);
    BOOST_CHECK_EQUAL(env_line_channels(2614.0).size(), 2u);
    BOOST_CHECK(env_line_channels(2614.0)[1].origin == LineOrigin::Annihilation);
    // A result built for one energy must not be silently run at another.
    EnvScatterConfig cfg;
    const Material air = make_Air();
    EnvironmentGeometry env(air);
    NextEventEstimator est(env, cfg, Eigen::Vector3d(0, 0, -100));
    EnvScatterResult wrong(cfg.axes, 662.0);
    BOOST_CHECK_THROW(est.run(1332.0, 10, 1, wrong), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Closure: the differential coefficient must integrate to the material's mu
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(differential_scatter_integrates_to_material_mu)
{
    // The estimator weights each vertex by mu_{s,e} * p_e(mu_s). If that product
    // did not integrate over the sphere to the SAME mu_cs / mu_rs the analog walk
    // samples its collision rate from, the estimator and the walk would disagree
    // by a constant factor -- a silent scale error invisible to any single-run
    // sanity check. This pins the normalization.
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Material concrete = make_Concrete();

    for (const Material* m : {&air, &soil, &concrete}) {
        for (double E : {59.5, 662.0, 2614.0}) {
            const MacroscopicXS xs = m->macroscopic_xs(E * 1e-3);
            // Integrate dMu/dOmega over 4pi: 2*pi * integral over mu.
            const double int_c = kTwoPi * gl::integrate_mu_in_s(
                -1.0, 1.0, 16, [&](double mu) {
                    return differential_scatter(*m, E, mu, true).compton;
                });
            const double int_r = kTwoPi * gl::integrate_mu_in_s(
                -1.0, 1.0, 16, [&](double mu) {
                    return differential_scatter(*m, E, mu, true).rayleigh;
                });
            BOOST_TEST_CONTEXT(m->name() << " @ " << E << " keV") {
                BOOST_CHECK_CLOSE(int_c, xs.mu_cs, 0.05);
                BOOST_CHECK_CLOSE(int_r, xs.mu_rs, 0.05);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// THE UNBIASEDNESS GATE: next-event estimator vs independent analog truth
// ---------------------------------------------------------------------------
//
// These compare the estimator against a fully independent analog truth that
// needs tens of millions of histories to converge, so each takes minutes. They
// are the most important tests in this file and MUST be run whenever the
// estimator changes:
//     cmake -B build -S . -DCEELO_RUN_MC_TESTS=ON && ctest -R env_scatter
// They are gated out of the default ctest run only because of their runtime.
#ifdef CEELO_RUN_MC_TESTS

BOOST_AUTO_TEST_CASE(estimator_matches_brute_force_single_scatter_only)
{
    // Restricting both estimators to ONE vertex isolates the per-vertex scoring
    // algebra (cone solid angle, chord, angular density, attenuation leg) from
    // the multiple-scatter bookkeeping. If this passes and the full multi-scatter
    // case does not, the fault is in the walk, not the score.
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Eigen::Vector3d src(0, 0, -40.0);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, soil, src, 25.0, 25.0, 60.0);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 5.0;
    cfg.max_vertices = 1;  // estimator: score at the first collision only

    // The brute needs ONE MORE vertex budget than the estimator to cover the
    // same physics: the estimator scores AT vertex k, whereas the brute tallies
    // the flight LEAVING vertex k, which it only reaches on the following
    // iteration.
    EnvScatterConfig brute_cfg = cfg;
    brute_cfg.max_vertices = cfg.max_vertices + 1;

    NextEventEstimator est(env, cfg, src);
    NextEventEstimator est_brute(env, brute_cfg, src);
    EnvScatterResult ne(cfg.axes, 662.0), br(cfg.axes, 662.0);
    est.run(662.0, 400000, 31337, ne);
    est_brute.run_brute(662.0, 40000000, 27182, br);

    const double m_ne = ne.tally.total_mean(), s_ne = ne.tally.total_sigma();
    const double m_br = br.tally.total_mean(), s_br = br.tally.total_sigma();
    const double z = (m_ne - m_br) / std::sqrt(s_ne * s_ne + s_br * s_br);
    BOOST_TEST_MESSAGE("single-scatter only: next-event " << m_ne << " +/- " << s_ne
                       << ", brute " << m_br << " +/- " << s_br
                       << "  (ratio " << m_ne / m_br << ", z = " << z << ")");
    BOOST_CHECK_LT(std::abs(z), 4.0);
}

BOOST_AUTO_TEST_CASE(estimator_matches_brute_force_in_air)
{
    // A close geometry where brute-force analog actually converges. The two
    // estimators share only the geometry and the physics kernels -- the scoring
    // paths are completely different (deterministic cone score vs. tallying real
    // track length inside the sphere), so agreement is a genuine gate on the
    // bounded-estimator algebra.
    const Material air = make_Air();
    const Eigen::Vector3d src(0, 0, -25.0);
    EnvironmentGeometry env =
        EnvironmentGeometry::open_air(air, src, Eigen::Vector3d::Zero(), 60.0, 30.0);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 4.0;
    cfg.e_floor_keV = 30.0;

    // The brute tallies the flight LEAVING a vertex, so it needs one more
    // vertex of budget than the estimator to cover the same physics.
    EnvScatterConfig brute_cfg = cfg;
    brute_cfg.max_vertices = cfg.max_vertices + 1;

    NextEventEstimator est(env, cfg, src);
    NextEventEstimator est_brute(env, brute_cfg, src);

    EnvScatterResult ne(cfg.axes, 662.0), br(cfg.axes, 662.0);
    est.run(662.0, 3000000, 12345, ne);
    est_brute.run_brute(662.0, 40000000, 6789, br);

    const double m_ne = ne.tally.total_mean(), s_ne = ne.tally.total_sigma();
    const double m_br = br.tally.total_mean(), s_br = br.tally.total_sigma();
    const double sigma = std::sqrt(s_ne * s_ne + s_br * s_br);
    const double z = (m_ne - m_br) / sigma;

    BOOST_TEST_MESSAGE("scattered fluence per source photon: next-event "
                       << m_ne << " +/- " << s_ne << ", brute " << m_br
                       << " +/- " << s_br << "  (z = " << z << ")");
    BOOST_CHECK_GT(m_ne, 0.0);
    BOOST_CHECK_GT(m_br, 0.0);
    // The z-test is the gate. A fixed-percent band is NOT appropriate here: open
    // air is exactly the regime where the estimator's variance is largest (most
    // primaries fly off and contribute nothing, while rare near-detector
    // vertices contribute a lot -- the very problem importance-sampled vertices
    // are meant to fix), so a band tighter than the noise would just be flaky.
    // The loose band below is a coarse order-of-magnitude guard only.
    BOOST_CHECK_LT(std::abs(z), 4.0);
    BOOST_CHECK_CLOSE(m_ne, m_br, 20.0);
}

BOOST_AUTO_TEST_CASE(estimator_matches_brute_force_over_ground)
{
    // Same gate with a dense slab in play: this exercises the in-slab vertex
    // path, where the vertex->detector leg must attenuate through the remaining
    // soil thickness before the air. That leg is what makes albedo scoring
    // correct, and it is not exercised by the air-only case at all.
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Eigen::Vector3d src(0, 0, -40.0);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, soil, src, 25.0, 25.0, 60.0);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 5.0;

    EnvScatterConfig brute_cfg = cfg;
    brute_cfg.max_vertices = cfg.max_vertices + 1;

    NextEventEstimator est(env, cfg, src);
    NextEventEstimator est_brute(env, brute_cfg, src);
    EnvScatterResult ne(cfg.axes, 662.0), br(cfg.axes, 662.0);
    est.run(662.0, 1000000, 999, ne);
    est_brute.run_brute(662.0, 40000000, 4242, br);

    const double m_ne = ne.tally.total_mean(), s_ne = ne.tally.total_sigma();
    const double m_br = br.tally.total_mean(), s_br = br.tally.total_sigma();
    const double sigma = std::sqrt(s_ne * s_ne + s_br * s_br);
    const double z = (m_ne - m_br) / sigma;

    BOOST_TEST_MESSAGE("ground geometry, scattered fluence: next-event "
                       << m_ne << " +/- " << s_ne << ", brute " << m_br
                       << " +/- " << s_br << "  (z = " << z << ")");
    // Note on the threshold: with a small scoring sphere the ANALOG side is
    // heavy-tailed (most histories score zero, a few score a lot), so its sample
    // variance systematically understates the true variance and z runs slightly
    // hot. The definitive gate is the high-acceptance case below, where the
    // analog truth converges tightly and the two agree to 0.01%.
    BOOST_CHECK_LT(std::abs(z), 4.0);
    BOOST_CHECK_CLOSE(m_ne, m_br, 10.0);

    // The ground must actually contribute: some chains have to touch the slab.
    const double slab_part = ne.tally.total_mean_class(ChainClass::OneBounce) +
                             ne.tally.total_mean_class(ChainClass::MultiSameSlab) +
                             ne.tally.total_mean_class(ChainClass::MultiDistinct);
    BOOST_TEST_MESSAGE("slab-touching fraction: " << slab_part / m_ne);
    BOOST_CHECK_GT(slab_part, 0.0);
}

BOOST_AUTO_TEST_CASE(estimator_matches_brute_force_in_a_room)
{
    // THE MULTI-SLAB GATE. Every other gate here uses open air or a SINGLE soil
    // slab, so the multi-slab tracing, the wall-to-wall next-event legs, and the
    // ChainClass::MultiDistinct labelling were entirely ungated -- and the
    // load-bearing "18-34% of room scatter has bounced twice or more" result
    // rides on exactly those paths. A room is a cheap place to gate: acceptance
    // is high (scattered flux exceeds direct), so the analog truth converges
    // fast with a generous scoring sphere.
    const Material air = make_Air();
    const Material concrete = make_Concrete();
    const Eigen::Vector3d src(0, 0, -200.0);
    EnvironmentGeometry env = EnvironmentGeometry::room(
        air, concrete, src, 100.0, 100.0,
        /*height*/ 250.0, /*half_width*/ 150.0, /*back*/ 100.0, /*front*/ 100.0);
    BOOST_REQUIRE_EQUAL(env.num_slabs(), 6u);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 15.0;

    EnvScatterConfig brute_cfg = cfg;
    brute_cfg.max_vertices = cfg.max_vertices + 1;

    NextEventEstimator est(env, cfg, src);
    NextEventEstimator est_brute(env, brute_cfg, src);
    EnvScatterResult ne(cfg.axes, 662.0), br(cfg.axes, 662.0);
    est.run(662.0, 300000, 20260803, ne);
    est_brute.run_brute(662.0, 8000000, 30260803, br);

    const double m_ne = ne.tally.total_mean(), s_ne = ne.tally.total_sigma();
    const double m_br = br.tally.total_mean(), s_br = br.tally.total_sigma();
    const double z = (m_ne - m_br) / std::sqrt(s_ne * s_ne + s_br * s_br);
    BOOST_TEST_MESSAGE("ROOM: next-event " << m_ne << " +/- " << s_ne
                       << ", brute " << m_br << " +/- " << s_br
                       << "  (ratio " << m_ne / m_br << ", z = " << z << ")");
    BOOST_CHECK_LT(std::abs(z), 4.0);
    BOOST_CHECK_CLOSE(m_ne, m_br, 5.0);

    // The multi-bounce component must be substantial AND must agree between the
    // two estimators -- that is the specific claim being gated, so checking only
    // the total would let a mislabelled chain class pass.
    const double mb_ne = ne.tally.total_mean_class(ChainClass::MultiDistinct) / m_ne;
    const double mb_br = br.tally.total_mean_class(ChainClass::MultiDistinct) / m_br;
    BOOST_TEST_MESSAGE("ROOM multi-slab fraction: next-event " << mb_ne
                       << ", brute " << mb_br);
    BOOST_CHECK_GT(mb_ne, 0.10);
    BOOST_CHECK_SMALL(std::abs(mb_ne - mb_br), 0.05);
}

BOOST_AUTO_TEST_CASE(estimator_matches_brute_force_with_annihilation)
{
    // THE ANNIHILATION GATE. Until the 511s were walked on BOTH sides this
    // channel was covered by nothing: the estimator scored the pair while
    // run_brute simply terminated at pair production, so the two were not
    // comparable above 1022 keV and no internal check could see the difference.
    // GEANT4 measured that gap as +8.3% at 2614 keV.
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Eigen::Vector3d src(0, 0, -30.0);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, soil, src, 20.0, 20.0, 40.0);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 15.0;
    EnvScatterConfig brute_cfg = cfg;
    brute_cfg.max_vertices = cfg.max_vertices + 1;

    NextEventEstimator est(env, cfg, src);
    NextEventEstimator est_brute(env, brute_cfg, src);
    EnvScatterResult ne(cfg.axes, 2614.0), br(cfg.axes, 2614.0);
    est.run(2614.0, 400000, 511511, ne);
    est_brute.run_brute(2614.0, 30000000, 611611, br);

    const double m_ne = ne.tally.total_mean(), s_ne = ne.tally.total_sigma();
    const double m_br = br.tally.total_mean(), s_br = br.tally.total_sigma();
    const double z = (m_ne - m_br) / std::sqrt(s_ne * s_ne + s_br * s_br);
    BOOST_TEST_MESSAGE("ANNIHILATION (2614 keV): next-event " << m_ne
                       << " +/- " << s_ne << ", brute " << m_br << " +/- " << s_br
                       << "  (ratio " << m_ne / m_br << ", z = " << z << ")");
    BOOST_CHECK_LT(std::abs(z), 4.0);
    BOOST_CHECK_CLOSE(m_ne, m_br, 5.0);

    // The channel must be populated, with BOTH pieces present: the uncollided
    // 511 line and the continuum from walking those photons.
    const double line511 = ne.tally.line_total_mean(1);
    BOOST_TEST_MESSAGE("  511 line = " << line511 << " ("
                       << 100.0 * line511 / m_ne << "% of scattered)");
    BOOST_CHECK_GT(line511, 0.0);
    BOOST_CHECK_GT(m_ne - line511, 0.0);
}

BOOST_AUTO_TEST_CASE(estimator_matches_brute_force_with_brems)
{
    // THE BREMS GATE. Same lesson as the annihilation gate above: the
    // estimator scores each emission next-event through the 2BS density AND
    // walks the analog photon, so run_brute must grow the same channel
    // (electrons formed at Compton/PE/pair vertices, photons pushed
    // scattered-from-birth) or the internal gate is blind to exactly the
    // term being added. GEANT4's origin decomposition measures brems at
    // 4.26% of scattered flux at 2614 keV in a soil field.
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Eigen::Vector3d src(0, 0, -30.0);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, soil, src, 20.0, 20.0, 40.0);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 15.0;
    EnvScatterConfig brute_cfg = cfg;
    brute_cfg.max_vertices = cfg.max_vertices + 1;

    NextEventEstimator est(env, cfg, src);
    NextEventEstimator est_brute(env, brute_cfg, src);
    EnvScatterResult ne(cfg.axes, 2614.0), br(cfg.axes, 2614.0);
    est.run(2614.0, 400000, 20260805, ne);
    est_brute.run_brute(2614.0, 30000000, 50806020, br);

    const double m_ne = ne.tally.total_mean(), s_ne = ne.tally.total_sigma();
    const double m_br = br.tally.total_mean(), s_br = br.tally.total_sigma();
    const double z = (m_ne - m_br) / std::sqrt(s_ne * s_ne + s_br * s_br);
    BOOST_TEST_MESSAGE("BREMS (2614 keV): next-event " << m_ne
                       << " +/- " << s_ne << ", brute " << m_br << " +/- " << s_br
                       << "  (ratio " << m_ne / m_br << ", z = " << z << ")");
    BOOST_CHECK_LT(std::abs(z), 4.0);
    BOOST_CHECK_CLOSE(m_ne, m_br, 5.0);

    // The channel itself, not just the total it hides in: both sides must
    // populate it, and they must agree. A z-gate, not a %-gate -- the brute
    // channel is a ~4% slice of an already-noisy analog tally.
    const double c_ne = ne.tally.total_mean_origin(FluxOrigin::Brems);
    const double cs_ne = ne.tally.total_sigma_origin(FluxOrigin::Brems);
    const double c_br = br.tally.total_mean_origin(FluxOrigin::Brems);
    const double cs_br = br.tally.total_sigma_origin(FluxOrigin::Brems);
    const double zc =
        (c_ne - c_br) / std::sqrt(cs_ne * cs_ne + cs_br * cs_br);
    BOOST_TEST_MESSAGE("  brems channel: next-event " << c_ne << " +/- " << cs_ne
                       << " (" << 100.0 * c_ne / m_ne << "% of scattered), brute "
                       << c_br << " +/- " << cs_br << ", z = " << zc);
    BOOST_CHECK_GT(c_ne, 0.0);
    BOOST_CHECK_GT(c_br, 0.0);
    BOOST_CHECK_LT(std::abs(zc), 4.0);

    // Shape spot-check: G4 says brems is 100% of the 30-33 keV bin. Weak
    // gate (> 50%) at these statistics; the real number comes from the
    // origin-tagged G4 comparison.
    const int ie_soft = cfg.axes.energy_bin(31.0);
    double soft_all = 0.0, soft_brems = 0.0;
    for (int imu = 0; imu < cfg.axes.n_mu; ++imu) {
        soft_all += ne.tally.mean(ie_soft, imu);
        soft_brems += ne.tally.mean_origin(ie_soft, imu, FluxOrigin::Brems);
    }
    BOOST_TEST_MESSAGE("  30-33 keV bin: brems share "
                       << (soft_all > 0.0 ? 100.0 * soft_brems / soft_all : 0.0)
                       << "%");
    if (soft_all > 0.0) BOOST_CHECK_GT(soft_brems / soft_all, 0.5);
}

BOOST_AUTO_TEST_CASE(estimator_matches_brute_force_high_acceptance)
{
    // A deliberately high-acceptance geometry: a large scoring sphere close to a
    // dense scatterer, so the analog truth converges to well under 1% and any
    // systematic in the estimator is exposed rather than buried in the tail
    // variance that plagues small-sphere/long-distance comparisons.
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Eigen::Vector3d src(0, 0, -30.0);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, soil, src, 20.0, 20.0, 40.0);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 15.0;

    EnvScatterConfig brute_cfg = cfg;
    brute_cfg.max_vertices = cfg.max_vertices + 1;  // see the single-scatter case

    NextEventEstimator est(env, cfg, src);
    NextEventEstimator est_brute(env, brute_cfg, src);
    EnvScatterResult ne(cfg.axes, 662.0), br(cfg.axes, 662.0);
    est.run(662.0, 400000, 5150, ne);
    est_brute.run_brute(662.0, 40000000, 8080, br);

    const double m_ne = ne.tally.total_mean(), s_ne = ne.tally.total_sigma();
    const double m_br = br.tally.total_mean(), s_br = br.tally.total_sigma();
    const double z = (m_ne - m_br) / std::sqrt(s_ne * s_ne + s_br * s_br);
    BOOST_TEST_MESSAGE("high acceptance: next-event " << m_ne << " +/- " << s_ne
                       << ", brute " << m_br << " +/- " << s_br
                       << "  (ratio " << m_ne / m_br << ", z = " << z << ")");
    BOOST_CHECK_LT(std::abs(z), 4.0);
    BOOST_CHECK_CLOSE(m_ne, m_br, 3.0);
}

#endif  // CEELO_RUN_MC_TESTS

BOOST_AUTO_TEST_CASE(origin_decomposition_partitions_the_total)
{
    // Fast (ungated) invariants of the FluxOrigin tag on a short real run:
    //   * the origins PARTITION the total exactly (identity, not statistics);
    //   * the Brems channel is populated when include_brems is on, and is
    //     exactly zero when off;
    //   * the ElasticE0 line channel is populated by the PRIMARY chain --
    //     the regression for the positional-aggregate-init bug that left
    //     elastic_line = -1 on the primary and silently emptied line 0
    //     while every total stayed correct.
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Eigen::Vector3d src(0, 0, -30.0);
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, soil, src, 20.0, 20.0, 40.0);

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = 15.0;

    for (double e0 : {662.0, 2614.0}) {
        NextEventEstimator est(env, cfg, src);
        EnvScatterResult ne(cfg.axes, e0);
        est.run(e0, 20000, 424242, ne);

        const double total = ne.tally.total_mean();
        const double parts =
            ne.tally.total_mean_origin(FluxOrigin::Primary) +
            ne.tally.total_mean_origin(FluxOrigin::Annihilation) +
            ne.tally.total_mean_origin(FluxOrigin::Brems);
        BOOST_REQUIRE_GT(total, 0.0);
        BOOST_CHECK_CLOSE(parts, total, 1e-9);

        BOOST_CHECK_GT(ne.tally.total_mean_origin(FluxOrigin::Brems), 0.0);
        // The fit taxonomy is creator-local: the analytic primary direct term
        // is absent, first and multiple collision products are separated, and
        // a newly created brems photon contributes at order zero.
        BOOST_CHECK_EQUAL(
            ne.tally.total_mean_component(FluxOrigin::Primary,
                                          PostCreationOrder::Zero),
            0.0);
        BOOST_CHECK_GT(
            ne.tally.total_mean_component(FluxOrigin::Primary,
                                          PostCreationOrder::One),
            0.0);
        BOOST_CHECK_GT(
            ne.tally.total_mean_component(FluxOrigin::Primary,
                                          PostCreationOrder::TwoPlus),
            0.0);
        BOOST_CHECK_GT(
            ne.tally.total_mean_component(FluxOrigin::Brems,
                                          PostCreationOrder::Zero),
            0.0);
        if (e0 >= 1022.0)
            BOOST_CHECK_GT(
                ne.tally.total_mean_origin(FluxOrigin::Annihilation), 0.0);
        if (e0 >= 1022.0)
            BOOST_CHECK_GT(
                ne.tally.total_mean_component(FluxOrigin::Annihilation,
                                              PostCreationOrder::Zero),
                0.0);

        // ElasticE0 (line 0) must be fed by the primary's elastic chain.
        BOOST_CHECK_GT(ne.tally.line_total_mean(0), 0.0);
    }

    // Brems off: the channel is identically zero and the rest still runs.
    EnvScatterConfig off = cfg;
    off.include_brems = false;
    NextEventEstimator est_off(env, off, src);
    EnvScatterResult ne_off(off.axes, 662.0);
    est_off.run(662.0, 20000, 424242, ne_off);
    BOOST_CHECK_EQUAL(ne_off.tally.total_mean_origin(FluxOrigin::Brems), 0.0);
    BOOST_CHECK_GT(ne_off.tally.total_mean(), 0.0);
}

BOOST_AUTO_TEST_CASE(direct_transmission_is_analytic_and_separate)
{
    // The peak term must be the analytic uncollided transmission, and must NOT
    // appear in the scattered tally -- they are mutually exclusive outcomes of
    // one primary.
    const Material air = make_Air();
    const Eigen::Vector3d src(0, 0, -1000.0);
    EnvironmentGeometry env(air);
    EnvScatterConfig cfg;
    NextEventEstimator est(env, cfg, src);

    const double t = est.direct_transmission(662.0);
    BOOST_CHECK_CLOSE(t, std::exp(-air.mu_total(0.662) * 1000.0), 1e-6);
    BOOST_CHECK_LT(t, 1.0);
    BOOST_CHECK_GT(t, 0.0);
}
