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

#define BOOST_TEST_MODULE DetectorResponseIoTests
#include <boost/test/unit_test.hpp>

// Storage/evaluation/XML tests for the parameterized detector response --
// all MC-free (fabricated tables on real geometries), so they run fast.

#include "geometry/Geometry.h"
#include "io/DetectorResponse.h"
#include "io/Pchip.h"
#include "io/ResponseKernel.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Build a minimal-but-complete response for a 3"x3" NaI + 0.5 mm Al can with
// a CONSTANT eta (so eps_fep == eta0 * K exactly and kernel effects separate
// from table effects in the checks below).
std::shared_ptr<DetectorResponse> make_synthetic_nai(double eta0 = 0.5,
                                                     bool with_grounding = false) {
    auto r = std::make_shared<DetectorResponse>();

    MaterialSpec nai = MaterialSpec::from(make_NaI());
    MaterialSpec al = MaterialSpec::from(make_Aluminum());
    r->descriptor.shape = DetectorShape::Cylinder;
    r->descriptor.dimensions_cm = {3.81, 7.62};
    r->descriptor.crystal_material_index = 0;
    r->descriptor.materials = {nai, al};
    LayerSpec can;
    can.material_index = 1;
    can.front_thickness_cm = 0.05;
    can.side_thickness_cm = 0.05;
    can.z_start_cm = 0.0;
    can.z_end_cm = 7.62;
    r->descriptor.layers.push_back(can);

    r->mu_tables.push_back(MuTable::sample(make_NaI(), 0));
    r->mu_tables.push_back(MuTable::sample(make_Aluminum(), 1));

    EtaTable& t = r->eta_fep;
    t.energies_keV = {35.0, 33.17 * 0.999, 33.17 * 1.001, 60.0, 150.0, 400.0,
                      1000.0, 3000.0};
    std::sort(t.energies_keV.begin(), t.energies_keV.end());
    t.cos_thetas = {0.02, 0.3, 0.6, 0.85, 1.0};
    t.edges_keV = {33.17};
    const size_t n = t.energies_keV.size() * t.cos_thetas.size();
    t.ln_eta.assign(n, std::log(eta0));
    t.frac_sigma.assign(n, 0.003);

    r->tot_eff.tier = TotEffTier::BCurve;
    r->tot_eff.b_energies_keV = {30.0, 100.0, 300.0, 1000.0, 3000.0};
    r->tot_eff.ln_b.assign(5, std::log(1.02));

    r->provenance.profile = ResponseProfile::General;
    r->provenance.valid_e_min_keV = 33.0;
    r->provenance.valid_e_max_keV = 3000.0;
    r->provenance.kernel_n_rays = 4096;
    r->provenance.detector_name = "synthetic NaI 3x3";

    if (with_grounding) {
        GroundingBlock& g = r->grounding;
        g.knot_ln_energies = {std::log(60.0), std::log(400.0), std::log(2000.0)};
        g.ln_k = {std::log(1.05), std::log(1.02), std::log(0.99)};
        g.cov = {4e-4, 1e-4, 0.0,
                 1e-4, 2.5e-4, 5e-5,
                 0.0, 5e-5, 9e-4};
        GroundingPoint p;
        p.energy_keV = 661.7;
        p.measured_eff = 1.1e-3;
        p.model_eff = 1.08e-3;
        p.frac_stat_sigma = 0.01;
        p.frac_cert_sigma = 0.03;
        p.source_key = "Cs137/cert1";
        p.distance_cm = 50.0;
        g.points.push_back(p);
    }

    r->finalize();
    return r;
}

}  // namespace

// --- Pchip vs scipy ---------------------------------------------------------

BOOST_AUTO_TEST_CASE(pchip_matches_scipy) {
    // Reference values from scipy.interpolate.PchipInterpolator (see the
    // generation snippet in the git history / CP-A2 notes).
    const std::vector<double> x{0.5, 1.0, 1.7, 2.1, 3.0, 4.6};
    const std::vector<double> y{1.2, 0.8, 0.95, 1.4, 1.38, 0.2};
    const Pchip f(x, y);
    const std::vector<double> q{0.5, 0.6, 0.9, 1.05, 1.3, 1.69, 1.9, 2.0,
                                2.5, 2.9, 3.5, 4.0, 4.59};
    const std::vector<double> ref{
        1.2, 1.080152380952381, 0.82203809523809523, 0.80091404502000141,
        0.8308461590616314, 0.94618068004610467, 1.1941860465116274,
        1.3368822674418601, 1.3955767878981793, 1.3838194550099014,
        1.2253447025903179, 0.84792549796007877, 0.21193907681715418};
    for (size_t i = 0; i < q.size(); ++i)
        BOOST_CHECK_CLOSE(f(q[i]), ref[i], 1e-9);

    // Log-log monotone-decreasing curve (the response's usual shape).
    const std::vector<double> x2{std::log(30.0), std::log(60.0), std::log(150.0),
                                 std::log(400.0), std::log(1000.0), std::log(3000.0)};
    const std::vector<double> y2{std::log(0.9), std::log(0.7), std::log(0.4),
                                 std::log(0.2), std::log(0.09), std::log(0.03)};
    const Pchip f2(x2, y2);
    const std::vector<double> q2{std::log(45.0), std::log(100.0),
                                 std::log(662.0), std::log(2614.0)};
    const std::vector<double> ref2{-0.23632633615678272, -0.6470127793176178,
                                   -2.0320600197215275, -3.3603763239033047};
    for (size_t i = 0; i < q2.size(); ++i)
        BOOST_CHECK_CLOSE(f2(q2[i]), ref2[i], 1e-9);

    // Clamping (no extrapolation).
    BOOST_CHECK_EQUAL(f(0.0), y.front());
    BOOST_CHECK_EQUAL(f(9.9), y.back());
}

BOOST_AUTO_TEST_CASE(pchip_monotone_no_overshoot) {
    // Monotone data must stay monotone (the PCHIP guarantee).
    const std::vector<double> x{1, 2, 3, 4, 5, 6};
    const std::vector<double> y{0.0, 0.1, 0.11, 3.0, 3.05, 3.1};
    const Pchip f(x, y);
    double prev = -1e9;
    for (double q = 1.0; q <= 6.0; q += 0.01) {
        const double v = f(q);
        BOOST_REQUIRE_GE(v + 1e-12, prev);
        prev = v;
    }
    BOOST_CHECK_LE(f(4.5), 3.1);  // no overshoot past the data range
}

// --- K-edge segmentation ----------------------------------------------------

BOOST_AUTO_TEST_CASE(eta_table_never_bridges_edge) {
    // eta drops 10% across the 33.17 keV NaI K-edge; the interpolant must
    // reproduce the discontinuity, not smooth it.
    auto r = make_synthetic_nai();
    EtaTable t = r->eta_fep;  // copy structure, then install the jump
    const double edge = 33.17;
    for (size_t e = 0; e < t.energies_keV.size(); ++e)
        for (size_t c = 0; c < t.cos_thetas.size(); ++c)
            t.ln_eta[t.index(e, c, 0)] =
                std::log(t.energies_keV[e] < edge ? 0.5 : 0.45);
    t.finalize();

    bool clamped = false;
    const double below = std::exp(t.eval_ln(edge * 0.9995, 1.0, 0.0, clamped));
    const double above = std::exp(t.eval_ln(edge * 1.0005, 1.0, 0.0, clamped));
    BOOST_CHECK_CLOSE(below, 0.5, 0.2);
    BOOST_CHECK_CLOSE(above, 0.45, 0.2);
}

// --- synthetic response consistency ----------------------------------------

BOOST_AUTO_TEST_CASE(eps_fep_equals_eta_times_kernel) {
    auto r = make_synthetic_nai(0.5);
    const Eigen::Vector3d pos = source_position(20.0, 1.0);
    const ApertureQuadrature q = r->make_quadrature(pos);
    const double K = r->kernel_K(662.0, q, MuChoice::Total);
    const EffResult res = r->eps_fep_at(662.0, pos, q);
    BOOST_CHECK_EQUAL(static_cast<int>(res.flag), static_cast<int>(ResponseFlag::Ok));
    BOOST_CHECK_CLOSE(res.value, 0.5 * K, 1e-6);
    BOOST_CHECK_GT(res.sigma, 0.0);

    // eps_total (BCurve tier) = b(E) * K_noRS.
    const double K_nors = r->kernel_K(662.0, q, MuChoice::NoRayleigh);
    const EffResult tot = r->eps_total_at(662.0, pos, q);
    BOOST_CHECK_CLOSE(tot.value, 1.02 * K_nors, 1e-6);
    BOOST_CHECK_GT(tot.value, res.value);
}

BOOST_AUTO_TEST_CASE(stored_mu_tables_match_live_xs) {
    // The stored-mu kernel path must reproduce the live-XS kernel to well
    // under the node precision (log-log interpolation on ~100 pts/decade).
    auto r = make_synthetic_nai();
    const Eigen::Vector3d pos = source_position(10.0, std::cos(30.0 * kPi / 180.0));
    const ApertureQuadrature q = r->make_quadrature(pos);
    for (const double E : {40.0, 60.0, 150.0, 662.0, 1332.0, 2614.0}) {
        const double stored = r->kernel_K(E, q, MuChoice::Total);
        const double live = q.interaction_omega(E, MuChoice::Total);
        BOOST_CHECK_CLOSE(stored, live, 0.1);
    }
    // Transmitted envelope too (collimator gate path).
    BOOST_CHECK_CLOSE(r->kernel_transmitted(120.0, q),
                      q.transmitted_omega(120.0), 0.1);
}

BOOST_AUTO_TEST_CASE(tsrc_folds_into_kernel) {
    auto r = make_synthetic_nai(0.5);
    const Eigen::Vector3d pos = source_position(5.0, 1.0);
    const ApertureQuadrature q = r->make_quadrature(pos);
    const std::function<double(const Eigen::Vector3d&)> half =
        [](const Eigen::Vector3d&) { return 0.5; };
    const EffResult full = r->eps_fep_at(662.0, pos, q);
    const EffResult attn = r->eps_fep_element(662.0, pos, q, half);
    BOOST_CHECK_CLOSE(attn.value, 0.5 * full.value, 1e-9);
}

// Build-up seam (stage E3 A1): with a stub model returning 1.3 and a supplied
// ShieldContext, eps_total scales by exactly 1.3 and its sigma inflates; with
// no ShieldContext (or no model) the value/sigma are byte-identical to the
// pre-seam path.  The FEP path never sees the seam.
BOOST_AUTO_TEST_CASE(buildup_seam_stub) {
    auto r = make_synthetic_nai(0.5);
    const Eigen::Vector3d pos = source_position(20.0, 1.0);
    const ApertureQuadrature q = r->make_quadrature(pos);

    const EffResult base = r->eps_total_at(662.0, pos, q);          // sc = null
    // Null-ShieldContext path is byte-identical even after a model is installed.
    r->buildup_model = [](double, const ShieldContext&) { return 1.3; };
    const EffResult still_base = r->eps_total_at(662.0, pos, q);    // sc = null
    BOOST_CHECK_EQUAL(still_base.value, base.value);
    BOOST_CHECK_EQUAL(still_base.sigma, base.sigma);

    // With a ShieldContext, eps_total scales by exactly 1.3; sigma inflates.
    ShieldContext sc;
    sc.areal_density_g_cm2 = 5.0;
    sc.eff_atomic_number = 26.0;
    const EffResult boosted = r->eps_total_at(662.0, pos, q, &sc);
    BOOST_CHECK_CLOSE(boosted.value, 1.3 * base.value, 1e-9);
    BOOST_CHECK_GT(boosted.sigma / boosted.value, base.sigma / base.value);
    // Build-up is ratio-only (>= 1): a model returning < 1 is clamped to 1.
    r->buildup_model = [](double, const ShieldContext&) { return 0.4; };
    const EffResult clamped = r->eps_total_at(662.0, pos, q, &sc);
    BOOST_CHECK_CLOSE(clamped.value, base.value, 1e-9);
    // FEP is never touched by the seam.
    const EffResult fep = r->eps_fep_at(662.0, pos, q);
    BOOST_CHECK_GT(fep.value, 0.0);
}

// --- flags ------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(flag_out_of_range_clamped) {
    auto r = make_synthetic_nai();
    const EffResult lo = r->eps_fep(10.0, 0.0, 0.0, 25.0);
    BOOST_CHECK(lo.flag == ResponseFlag::OutOfRangeClamped);
    const EffResult hi = r->eps_fep(8000.0, 0.0, 0.0, 25.0);
    BOOST_CHECK(hi.flag == ResponseFlag::OutOfRangeClamped);
    const EffResult ok = r->eps_fep(662.0, 0.0, 0.0, 25.0);
    BOOST_CHECK(ok.flag == ResponseFlag::Ok);
}

BOOST_AUTO_TEST_CASE(flag_near_field_unmodeled_for_farfield_profile) {
    auto r = make_synthetic_nai();
    // Far-field profile: no near model, validity floor at 2a.
    r->provenance.profile = ResponseProfile::FarField;
    r->provenance.min_distance_cm = 2.0 * r->transverse_half_extent();
    const EffResult close = r->eps_fep(662.0, 0.0, 0.0, 2.0);
    BOOST_CHECK(close.flag == ResponseFlag::NearFieldUnmodeled);
    const EffResult far = r->eps_fep(662.0, 0.0, 0.0, 100.0);
    BOOST_CHECK(far.flag == ResponseFlag::Ok);
    // The unmodeled-near sigma must be inflated relative to the far query.
    BOOST_CHECK_GT(close.sigma / close.value, far.sigma / far.value);
}

BOOST_AUTO_TEST_CASE(flag_behind_plane_needs_mc) {
    auto r = make_synthetic_nai();
    const EffResult behind = r->eps_fep(662.0, 2.5 /*~143 deg*/, 0.0, 30.0);
    BOOST_CHECK(behind.flag == ResponseFlag::NeedsMc);
}

BOOST_AUTO_TEST_CASE(flag_collimator_shadow) {
    // Z8-style W collimator; a steep side view is shadow-dominated.
    auto r = std::make_shared<DetectorResponse>();
    r->descriptor.shape = DetectorShape::Cylinder;
    r->descriptor.dimensions_cm = {3.81, 7.62};
    r->descriptor.crystal_material_index = 0;
    r->descriptor.materials = {MaterialSpec::from(make_NaI()),
                               MaterialSpec::from(make_Aluminum()),
                               MaterialSpec::from(make_Tungsten())};
    LayerSpec can;
    can.material_index = 1;
    can.front_thickness_cm = 0.05;
    can.side_thickness_cm = 0.05;
    can.z_end_cm = 7.62;
    r->descriptor.layers.push_back(can);
    CollimatorSpec col;
    col.material_index = 2;
    col.side_thickness_cm = 1.5;
    col.z_start_cm = -5.0;
    col.z_end_cm = 7.62;
    r->descriptor.collimator = col;
    r->mu_tables.push_back(MuTable::sample(make_NaI(), 0));
    r->mu_tables.push_back(MuTable::sample(make_Aluminum(), 1));
    r->mu_tables.push_back(MuTable::sample(make_Tungsten(), 2));
    r->eta_fep.energies_keV = {40.0, 200.0, 662.0, 3000.0};
    r->eta_fep.cos_thetas = {0.02, 0.5, 1.0};
    r->eta_fep.ln_eta.assign(12, std::log(0.5));
    r->eta_fep.frac_sigma.assign(12, 0.003);
    r->provenance.kernel_n_rays = 4096;
    r->finalize();

    const EffResult open = r->eps_fep(662.0, 0.0, 0.0, 30.0);
    BOOST_CHECK(open.flag == ResponseFlag::Ok);
    const EffResult deep = r->eps_fep(200.0, 75.0 * kPi / 180.0, 0.0, 30.0);
    BOOST_CHECK(deep.flag == ResponseFlag::Shadowed ||
                deep.flag == ResponseFlag::NeedsMc);
    BOOST_CHECK_GT(deep.sigma / std::max(deep.value, 1e-300),
                   open.sigma / open.value);
}

// --- grounding --------------------------------------------------------------

BOOST_AUTO_TEST_CASE(grounding_k_and_covariance) {
    auto r = make_synthetic_nai(0.5, /*with_grounding=*/true);
    auto plain = make_synthetic_nai(0.5, /*with_grounding=*/false);

    // k at a knot is exact; between knots piecewise-linear in ln E.
    const EffResult at_knot = r->eps_fep(400.0, 0.0, 0.0, 50.0);
    const EffResult base = plain->eps_fep(400.0, 0.0, 0.0, 50.0);
    BOOST_CHECK_CLOSE(at_knot.value / base.value, 1.02, 1e-6);

    // Covariance: exact at knots, symmetric, strongly correlated nearby.
    BOOST_CHECK_CLOSE(r->grounding.var_ln_k(400.0), 2.5e-4, 1e-6);
    BOOST_CHECK_CLOSE(r->grounding.cov_ln_k(60.0, 400.0), 1e-4, 1e-6);
    BOOST_CHECK_CLOSE(r->grounding.cov_ln_k(60.0, 400.0),
                      r->grounding.cov_ln_k(400.0, 60.0), 1e-9);

    // Outside the calibrated range: clamped + flagged.
    const EffResult outside = r->eps_fep(2900.0, 0.0, 0.0, 50.0);
    BOOST_CHECK(outside.flag == ResponseFlag::OutOfRangeClamped);

    // Multi-energy covariance: PSD-ish sanity + near-energy correlation.
    const std::vector<double> Es{121.8, 344.3, 661.7, 1408.0};
    const std::vector<double> C = r->frac_covariance(Es, 0.0, 50.0);
    const size_t n = Es.size();
    for (size_t i = 0; i < n; ++i) {
        BOOST_CHECK_GT(C[i * n + i], 0.0);
        for (size_t j = 0; j < n; ++j) {
            BOOST_CHECK_CLOSE(C[i * n + j], C[j * n + i], 1e-9);
            const double rho = C[i * n + j] /
                std::sqrt(C[i * n + i] * C[j * n + j]);
            BOOST_CHECK_LE(std::fabs(rho), 1.0 + 1e-9);
            if (i != j) BOOST_CHECK_GT(rho, 0.3);  // strongly correlated
        }
    }
}

BOOST_AUTO_TEST_CASE(sigma_transfer_shape) {
    const SigmaTransferModel m;
    // Far on-axis small; off-axis grows; low-E grazing large; near grows.
    const double far_on = m.eval(20.0, 1.0, 662.0);
    const double far_off = m.eval(20.0, std::cos(75.0 * kPi / 180.0), 662.0);
    const double low_graze = m.eval(20.0, 0.05, 45.0);
    const double contact = m.eval(1.0, 1.0, 662.0);
    BOOST_CHECK_LT(far_on, 0.01);
    BOOST_CHECK_GT(far_off, far_on);
    BOOST_CHECK_GT(low_graze, 0.15);
    BOOST_CHECK_GT(contact, 0.05);
}

// --- XML round trip ---------------------------------------------------------

BOOST_AUTO_TEST_CASE(xml_round_trip_bit_stable) {
    auto r = make_synthetic_nai(0.5, /*with_grounding=*/true);
    // Exercise the optional blocks too: a small 2 x 3 x 4 near-field table.
    NearFieldModel& nf = r->near_field;
    nf.energies_keV = {60.0, 662.0};
    nf.cos_thetas = {0.02, 0.5, 1.0};
    nf.dists_cm = {2.0, 5.0, 12.0, 20.0};
    nf.ln_n.assign(2 * 3 * 4, 0.0);
    for (size_t i = 0; i < nf.ln_n.size(); ++i)
        nf.ln_n[i] = 0.12 - 0.01 * double(i);   // arbitrary smooth values
    nf.frac_sigma.assign(nf.ln_n.size(), 0.005);
    nf.break_cos_thetas = {0.02, 1.0};
    nf.break_d_cm = {8.0, 6.0, 7.0, 5.0};
    nf.finalize();   // r was finalized before these fields were set

    const std::string xml1 = r->to_xml_string();
    std::shared_ptr<DetectorResponse> r2 = DetectorResponse::from_xml_string(xml1);
    const std::string xml2 = r2->to_xml_string();
    BOOST_REQUIRE_EQUAL(xml1, xml2);
    BOOST_CHECK_EQUAL(r->content_hash(), r2->content_hash());

    // And the reloaded object evaluates identically.
    for (const double E : {40.0, 121.8, 662.0, 2614.0}) {
        const EffResult a = r->eps_fep(E, 0.3, 0.1, 12.0);
        const EffResult b = r2->eps_fep(E, 0.3, 0.1, 12.0);
        BOOST_CHECK_CLOSE(a.value, b.value, 1e-12);
        BOOST_CHECK_CLOSE(a.sigma, b.sigma, 1e-12);
        BOOST_CHECK(a.flag == b.flag);
    }
}

BOOST_AUTO_TEST_CASE(near_field_table_eval) {
    // ln N = f(ct, ln d) linear in both coords and equal across energies, so
    // PCHIP is exact at nodes AND midpoints (monotone cubic reproduces linear
    // data), and the energy blend is a no-op.
    auto f = [](double ct, double d) { return 0.3 * ct + 0.1 * std::log(d); };

    NearFieldModel nf;
    nf.energies_keV = {60.0, 600.0};
    nf.cos_thetas = {0.1, 0.4, 0.7, 1.0};
    nf.dists_cm = {1.0, 3.0, 9.0};
    const size_t nc = nf.cos_thetas.size();
    const size_t nd = nf.dists_cm.size();
    nf.ln_n.assign(nf.energies_keV.size() * nc * nd, 0.0);
    nf.frac_sigma.assign(nf.ln_n.size(), 0.0);
    for (size_t e = 0; e < nf.energies_keV.size(); ++e)
        for (size_t c = 0; c < nc; ++c)
            for (size_t d = 0; d < nd; ++d) {
                nf.ln_n[nf.index(e, c, d)] = f(nf.cos_thetas[c], nf.dists_cm[d]);
                nf.frac_sigma[nf.index(e, c, d)] =
                    0.001 * double(c + 1) * double(d + 1);
            }

    // (e) ln_boost before finalize() throws.
    BOOST_CHECK_THROW(nf.ln_boost(100.0, 0.5, 5.0), std::exception);

    nf.finalize();

    // (a) exact at nodes.
    for (size_t c = 0; c < nc; ++c)
        for (size_t d = 0; d < nd; ++d)
            BOOST_CHECK_CLOSE(nf.ln_boost(300.0, nf.cos_thetas[c], nf.dists_cm[d]),
                              f(nf.cos_thetas[c], nf.dists_cm[d]), 1e-9);

    // (b) exact at an interior midpoint (linear data).
    BOOST_CHECK_CLOSE(nf.ln_boost(300.0, 0.55, 2.0), f(0.55, 2.0), 1e-6);

    // (c) clamping outside the grid.
    BOOST_CHECK_CLOSE(nf.ln_boost(300.0, 0.5, 0.2),      // d below front -> front
                      f(0.5, nf.dists_cm.front()), 1e-6);
    BOOST_CHECK_CLOSE(nf.ln_boost(300.0, 0.5, 100.0),    // d above back  -> back
                      f(0.5, nf.dists_cm.back()), 1e-6);
    BOOST_CHECK_CLOSE(nf.ln_boost(300.0, 0.01, 3.0),     // ct below front -> front
                      f(nf.cos_thetas.front(), 3.0), 1e-6);

    // (d) node_frac_sigma exact at a node.
    BOOST_CHECK_CLOSE(nf.node_frac_sigma(60.0, nf.cos_thetas[2], nf.dists_cm[1]),
                      0.001 * 3.0 * 2.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(xml_rejects_future_version) {
    auto r = make_synthetic_nai();
    std::string xml = r->to_xml_string();
    const size_t pos = xml.find("version=\"1\"");
    BOOST_REQUIRE(pos != std::string::npos);
    xml.replace(pos, 11, "version=\"99\"");
    BOOST_CHECK_THROW(DetectorResponse::from_xml_string(xml), std::exception);
}

// --- Accuracy certificate (additive metadata) -------------------------------

namespace {
AccuracyCertificate make_synthetic_certificate() {
    AccuracyCertificate c;
    c.converged = true;
    c.iterations = 0;
    c.cpu_seconds = 12.3456789;
    c.probe_seed_base = 1;
    c.fep_median = 0.0071;
    c.fep_p95 = 0.0138;
    c.fep_max = 0.0301;
    c.tot_median = 0.0092;
    c.tot_p95 = 0.0164;
    for (int i = 0; i < 5; ++i) {
        AccuracyCertificate::Row r;
        r.E_keV = 60.0 + 500.0 * i;
        r.d_cm = 0.5 + 3.0 * i;
        r.cos_theta = 0.2 + 0.15 * i;
        r.phi_deg = 0.0;
        r.mc = 1.0e-3 * (i + 1);
        r.mc_sig = 3.0e-6 * (i + 1);
        r.model = 1.01e-3 * (i + 1);
        r.model_sig = 1.4e-5 * (i + 1);
        r.tag = 0;
        r.pass = (i != 2);
        c.rows.push_back(r);
    }
    return c;
}
}  // namespace

BOOST_AUTO_TEST_CASE(certificate_round_trip) {
    auto r = make_synthetic_nai();
    r->certificate = make_synthetic_certificate();

    const std::string xml1 = r->to_xml_string();
    BOOST_CHECK(xml1.find("<Certificate") != std::string::npos);

    std::shared_ptr<DetectorResponse> r2 = DetectorResponse::from_xml_string(xml1);
    const AccuracyCertificate& c2 = r2->certificate;
    BOOST_REQUIRE(!c2.empty());
    BOOST_CHECK_EQUAL(c2.converged, true);
    BOOST_CHECK_EQUAL(c2.iterations, 0);
    BOOST_CHECK_CLOSE(c2.cpu_seconds, 12.3456789, 1e-12);
    BOOST_CHECK_EQUAL(c2.probe_seed_base, 1u);
    BOOST_CHECK_CLOSE(c2.fep_p95, 0.0138, 1e-12);
    BOOST_CHECK_CLOSE(c2.tot_p95, 0.0164, 1e-12);
    BOOST_REQUIRE_EQUAL(c2.rows.size(), 5u);
    for (size_t i = 0; i < c2.rows.size(); ++i) {
        const AccuracyCertificate::Row& a = r->certificate.rows[i];
        const AccuracyCertificate::Row& b = c2.rows[i];
        BOOST_CHECK_CLOSE(a.E_keV, b.E_keV, 1e-12);
        BOOST_CHECK_CLOSE(a.d_cm, b.d_cm, 1e-12);
        BOOST_CHECK_CLOSE(a.cos_theta, b.cos_theta, 1e-12);
        BOOST_CHECK_CLOSE(a.mc, b.mc, 1e-12);
        BOOST_CHECK_CLOSE(a.model, b.model, 1e-12);
        BOOST_CHECK_CLOSE(a.model_sig, b.model_sig, 1e-12);
        BOOST_CHECK_EQUAL(a.tag, b.tag);
        BOOST_CHECK_EQUAL(a.pass, b.pass);
    }

    // save/load/save is string-identical.
    const std::string xml2 = r2->to_xml_string();
    BOOST_CHECK_EQUAL(xml1, xml2);
}

BOOST_AUTO_TEST_CASE(certificate_absence_tolerated) {
    // A response WITHOUT a certificate has no <Certificate> element and loads
    // to an empty certificate (back-compat with pre-D-a XML).
    auto r = make_synthetic_nai();
    const std::string xml = r->to_xml_string();
    BOOST_CHECK(xml.find("<Certificate") == std::string::npos);
    std::shared_ptr<DetectorResponse> r2 = DetectorResponse::from_xml_string(xml);
    BOOST_CHECK(r2->certificate.empty());
}

BOOST_AUTO_TEST_CASE(certificate_invariant_content_hash) {
    // The certificate is metadata: adding one must NOT change content_hash,
    // and a with/without pair must hash identically.
    auto bare = make_synthetic_nai();
    const uint64_t h_bare = bare->content_hash();

    auto certd = make_synthetic_nai();
    certd->certificate = make_synthetic_certificate();
    const uint64_t h_certd = certd->content_hash();
    BOOST_CHECK_EQUAL(h_bare, h_certd);

    // And a reload of the certificate-bearing XML still hashes the same.
    std::shared_ptr<DetectorResponse> reload =
        DetectorResponse::from_xml_string(certd->to_xml_string());
    BOOST_CHECK(!reload->certificate.empty());
    BOOST_CHECK_EQUAL(reload->content_hash(), h_bare);
}
