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

#define BOOST_TEST_MODULE EfficiencyTransferTests
#include <boost/test/unit_test.hpp>

// EFFTRAN-style efficiency-transfer self-checks. All MC-free (analytic kernel
// ratios + table round-trips), so the suite runs in well under a second. The
// MC-vs-transfer consistency check is guarded behind CEELO_RUN_MC_TESTS.

#include "geometry/Geometry.h"
#include "io/DetectorResponse.h"
#include "io/EfficiencyTransfer.h"
#include "io/ResponseKernel.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;

const Material& mat_NaI() { static Material m = make_NaI();      return m; }
const Material& mat_Ge()  { static Material m = make_HPGe();     return m; }
const Material& mat_Al()  { static Material m = make_Aluminum(); return m; }

// Bare 3"x3" NaI cylinder for the standalone (live-mu) transfer tests.
Geometry bare_nai_3x3() {
    Geometry g;
    g.set_detector(DetectorShape::Cylinder, &mat_NaI(), {3.81, 7.62});
    return g;
}

// HPGe with a Ge dead layer + Al endcap -- passive layers whose Compton
// scatter-in the recapture term credits.
Geometry hpge_with_passive() {
    Geometry g;
    g.set_detector(DetectorShape::Cylinder, &mat_Ge(), {3.4, 5.5});
    g.set_dead_layer(0.07, 0.07, 0.0);
    g.add_attenuator(&mat_Al(), 0.15, 0.10, -0.5, 6.0);
    return g;
}

// A NaI 3"x3" + 0.5 mm Al can descriptor for the DetectorResponse producer.
GeometryDescriptor nai3x3_descriptor() {
    GeometryDescriptor gd;
    gd.shape = DetectorShape::Cylinder;
    gd.dimensions_cm = {3.81, 7.62};
    gd.materials = {MaterialSpec::from(make_NaI()),
                    MaterialSpec::from(make_Aluminum())};
    gd.crystal_material_index = 0;
    LayerSpec can;
    can.material_index = 1;
    can.front_thickness_cm = 0.05;
    can.side_thickness_cm = 0.05;
    can.z_end_cm = 7.62;
    gd.layers.push_back(can);
    return gd;
}

// A plausible (arbitrary but positive, decreasing) anchor curve. The transfer
// identities do not depend on physical correctness of these values.
AnchorCurve make_anchor(double scale = 1.0) {
    AnchorCurve a;
    a.energies_keV = {60.0, 100.0, 300.0, 662.0, 1000.0};
    a.eff = {2.0e-2, 1.3e-2, 5.0e-3, 3.0e-3, 2.2e-3};
    for (double& e : a.eff) e *= scale;
    a.frac_sigma = {0.003, 0.003, 0.003, 0.003, 0.003};
    return a;
}

}  // namespace

// The transfer to the SAME reference position returns the anchor exactly:
// eps(E_a, ref) = eta_ref(E_a) * K(E_a, ref) = (eff/K) * K = eff.
BOOST_AUTO_TEST_CASE(same_position_identity) {
    const Geometry g = bare_nai_3x3();
    const Eigen::Vector3d ref(0.0, 0.0, -40.0);
    const AnchorCurve a = make_anchor();
    EfficiencyTransfer et(g, ref, a, MuChoice::Total, 8192);
    const ApertureQuadrature q_ref = et.make_target_quadrature(ref);
    for (size_t i = 0; i < a.energies_keV.size(); ++i)
        BOOST_CHECK_CLOSE(et.eps_at(a.energies_keV[i], q_ref), a.eff[i], 1e-4);
}

// The transfer factor is exactly the kernel ratio, and efficiency falls with
// distance at a fixed angle.
BOOST_AUTO_TEST_CASE(transfer_is_kernel_ratio) {
    const Geometry g = bare_nai_3x3();
    const Eigen::Vector3d ref(0.0, 0.0, -40.0);
    EfficiencyTransfer et(g, ref, make_anchor(), MuChoice::Total, 8192);

    const Eigen::Vector3d t1(0.0, 0.0, -20.0);
    const Eigen::Vector3d t2(0.0, 0.0, -60.0);
    const ApertureQuadrature q1 = et.make_target_quadrature(t1);
    const ApertureQuadrature q2 = et.make_target_quadrature(t2);
    const double E = 662.0;
    const double ratio_eps = et.eps_at(E, q1) / et.eps_at(E, q2);
    const double ratio_K = q1.interaction_omega(E, MuChoice::Total) /
                           q2.interaction_omega(E, MuChoice::Total);
    BOOST_CHECK_CLOSE(ratio_eps, ratio_K, 1e-6);
    BOOST_CHECK_GT(et.eps_at(E, q1), et.eps_at(E, q2));  // nearer = higher
}

// A unit source-transmission functor must reproduce the plain transfer.
BOOST_AUTO_TEST_CASE(tsrc_identity) {
    const Geometry g = bare_nai_3x3();
    const Eigen::Vector3d ref(0.0, 0.0, -40.0);
    EfficiencyTransfer et(g, ref, make_anchor(), MuChoice::Total, 4096);
    const ApertureQuadrature q = et.make_target_quadrature({0.0, 0.0, -25.0});
    const double E = 300.0;
    const double plain = et.eps_at(E, q);
    const double with_unit =
        et.eps_at_with_tsrc(E, q, [](const Eigen::Vector3d&) { return 1.0; });
    BOOST_CHECK_CLOSE(plain, with_unit, 1e-9);
}

// Behind the crystal-face plane (cos_theta < 0) the position overload refuses.
BOOST_AUTO_TEST_CASE(behind_plane_refused) {
    const Geometry g = bare_nai_3x3();
    const Eigen::Vector3d ref(0.0, 0.0, -40.0);
    EfficiencyTransfer et(g, ref, make_anchor(), MuChoice::Total, 2048);
    const double v = et.eps_at(662.0, Eigen::Vector3d(0.0, 0.0, +20.0));
    BOOST_CHECK(std::isnan(v));
}

// make_transfer_response reproduces the anchor at the reference geometry (the
// eta table passes through its nodes), and off-axis inflates sigma via the
// attached SigmaTransferModel.
BOOST_AUTO_TEST_CASE(response_reproduces_anchor_and_inflates_offaxis) {
    const GeometryDescriptor gd = nai3x3_descriptor();
    const AnchorCurve a = make_anchor();
    const double a_ext =
        3.81 + 0.05;                    // outer radius incl. can (transverse a)
    const double d_ref = 10.0 * a_ext;  // far field, on axis
    const Eigen::Vector3d ref(0.0, 0.0, -d_ref);

    TransferResponseOptions opts;
    opts.detector_name = "nai3x3-transfer";
    auto resp = make_transfer_response(gd, a, ref, nullptr, opts);

    // At each anchor energy, on-axis far, eps_fep == anchor eff.
    for (size_t i = 0; i < a.energies_keV.size(); ++i) {
        const EffResult r = resp->eps_fep(a.energies_keV[i], 0.0, 0.0, d_ref);
        BOOST_CHECK_CLOSE(r.value, a.eff[i], 0.2);
    }

    // Off-axis sigma > on-axis sigma at the same distance (model_transfer).
    const EffResult on = resp->eps_fep(662.0, 0.0, 0.0, d_ref);
    const EffResult off =
        resp->eps_fep(662.0, 60.0 * kPi / 180.0, 0.0, d_ref);
    BOOST_CHECK(resp->model_transfer.has_value());
    BOOST_CHECK_GT(off.sigma / off.value, on.sigma / on.value);
}

// XML round-trips bit-stable and preserves model_transfer + evaluation.
BOOST_AUTO_TEST_CASE(response_xml_roundtrip) {
    const GeometryDescriptor gd = nai3x3_descriptor();
    const Eigen::Vector3d ref(0.0, 0.0, -40.0);
    auto resp = make_transfer_response(gd, make_anchor(), ref);

    const std::string xml = resp->to_xml_string();
    auto resp2 = DetectorResponse::from_xml_string(xml);
    BOOST_CHECK(resp2->model_transfer.has_value());
    BOOST_CHECK_EQUAL(resp2->to_xml_string(), xml);  // stable

    const EffResult r1 = resp->eps_fep(500.0, 20.0 * kPi / 180.0, 0.0, 30.0);
    const EffResult r2 = resp2->eps_fep(500.0, 20.0 * kPi / 180.0, 0.0, 30.0);
    BOOST_CHECK_CLOSE(r1.value, r2.value, 1e-9);
    BOOST_CHECK_CLOSE(r1.sigma, r2.sigma, 1e-9);
}

// A total anchor above the FEP anchor transfers to eps_total >= eps_fep.
BOOST_AUTO_TEST_CASE(total_ge_fep) {
    const GeometryDescriptor gd = nai3x3_descriptor();
    const AnchorCurve fep = make_anchor();
    AnchorCurve tot = make_anchor(2.0);  // total = 2x fep at the anchor
    const Eigen::Vector3d ref(0.0, 0.0, -40.0);
    auto resp = make_transfer_response(gd, fep, ref, &tot);

    for (double d : {15.0, 40.0, 80.0}) {
        const EffResult f = resp->eps_fep(662.0, 0.0, 0.0, d);
        const EffResult t = resp->eps_total(662.0, 0.0, 0.0, d);
        BOOST_CHECK_GE(t.value, f.value * (1.0 - 1e-6));
    }
}

// Transferring a response built with make_transfer_response is self-consistent
// with a hand-computed kernel-ratio transfer using the SAME stored mu tables.
BOOST_AUTO_TEST_CASE(response_self_consistent_kernel_ratio) {
    const GeometryDescriptor gd = nai3x3_descriptor();
    const AnchorCurve a = make_anchor();
    const Eigen::Vector3d ref(0.0, 0.0, -40.0);
    auto resp = make_transfer_response(gd, a, ref);

    const ApertureQuadrature q_ref = resp->make_quadrature(ref);
    const Eigen::Vector3d tgt(0.0, 0.0, -70.0);
    const ApertureQuadrature q_tgt = resp->make_quadrature(tgt);
    const double E = 662.0;  // an anchor energy
    // Hand transfer: eta = eff/K_ref, eps = eta * K_tgt (stored mu on both).
    const double K_ref = resp->kernel_K(E, q_ref, MuChoice::Total);
    const double K_tgt = resp->kernel_K(E, q_tgt, MuChoice::Total);
    const double eps_manual = (a.eff[3] / K_ref) * K_tgt;  // a.eff[3] == E 662
    const EffResult r = resp->eps_fep_at(E, tgt, q_tgt);
    BOOST_CHECK_CLOSE(r.value, eps_manual, 0.05);
}

// The scatter-in recapture raises the TOTAL kernel (credits passive-layer
// Compton) but leaves the FEP kernel bit-identical, and is a no-op on a bare
// crystal (no passive segments).
BOOST_AUTO_TEST_CASE(recapture_raises_total_not_fep) {
    const Geometry g = hpge_with_passive();
    const Eigen::Vector3d src = source_position(2.0 * 3.5, 1.0);  // near-field
    const ApertureQuadrature q = make_aperture_quadrature(g, src, 16384);
    const double E = 662.0;
    const double tot0 = q.interaction_omega(E, MuChoice::NoRayleigh, 0.0);
    const double tot8 = q.interaction_omega(E, MuChoice::NoRayleigh, 0.8);
    BOOST_CHECK_GT(tot8, tot0);                       // scatter-in credit
    BOOST_CHECK_LE(tot8, q.omega_frac_active + 1e-9); // still bounded
    // FEP kernel ignores the recapture.
    BOOST_CHECK_EQUAL(q.interaction_omega(E, MuChoice::Total, 0.8),
                      q.interaction_omega(E, MuChoice::Total, 0.0));
}

BOOST_AUTO_TEST_CASE(recapture_noop_on_bare_crystal) {
    const Geometry g = bare_nai_3x3();
    const ApertureQuadrature q = make_aperture_quadrature(g, {0, 0, -8.0}, 8192);
    BOOST_CHECK_EQUAL(q.interaction_omega(662.0, MuChoice::NoRayleigh, 0.8),
                      q.interaction_omega(662.0, MuChoice::NoRayleigh, 0.0));
}

// make_transfer_response with a total anchor + recapture still reproduces the
// total anchor at the reference geometry (the folded-kernel ratio is 1 there),
// and the recapture round-trips through XML.
BOOST_AUTO_TEST_CASE(recapture_response_reproduces_total_anchor) {
    const GeometryDescriptor gd = nai3x3_descriptor();
    const AnchorCurve fep = make_anchor();
    AnchorCurve tot = make_anchor(2.0);
    const double d_ref = 10.0 * (3.81 + 0.05);
    const Eigen::Vector3d ref(0.0, 0.0, -d_ref);
    TransferResponseOptions opts;
    opts.scatter_in_recapture = 0.8;
    auto resp = make_transfer_response(gd, fep, ref, &tot, opts);
    BOOST_CHECK_EQUAL(resp->scatter_in_recapture, 0.8);
    for (size_t i = 0; i < tot.energies_keV.size(); ++i) {
        const EffResult r = resp->eps_total(tot.energies_keV[i], 0.0, 0.0, d_ref);
        BOOST_CHECK_CLOSE(r.value, tot.eff[i], 0.3);
    }
    auto resp2 = DetectorResponse::from_xml_string(resp->to_xml_string());
    BOOST_CHECK_EQUAL(resp2->scatter_in_recapture, 0.8);
    BOOST_CHECK_EQUAL(resp2->to_xml_string(), resp->to_xml_string());
}

#ifdef CEELO_RUN_MC_TESTS
#include "efficiency/EfficiencyCalculator.h"
#include "io/ResponseGenerator.h"

// Coarse MC anchor on-axis, transfer to two far-field targets, compare to
// direct MC. |z| <= 3 is the acceptance bar (transfer ~1% far-field).
BOOST_AUTO_TEST_CASE(mc_consistency) {
    const GeometryDescriptor gd = nai3x3_descriptor();
    std::vector<std::unique_ptr<Material>> owned;
    const Geometry geom = gd.build_geometry(owned);

    EfficiencyCalculator calc;
    ResponseGenerator::configure_calculator(calc, gd, owned);

    const double a_ext = 3.81 + 0.05;
    const double d_ref = 10.0 * a_ext;
    const std::vector<double> Es{100.0, 662.0};
    AnchorCurve anchor;
    for (double E : Es) {
        calc.set_point_source({0.0, 0.0, -d_ref});
        SimulationConfig cfg;
        cfg.energy_keV = E;
        cfg.termination.target_fep_rel_precision = 0.01;
        cfg.termination.max_events = 4000000;
        cfg.seed = 12345;
        const EfficiencyResult r = calc.compute(cfg);
        anchor.energies_keV.push_back(E);
        anchor.eff.push_back(r.full_energy_peak_efficiency);
        anchor.frac_sigma.push_back(r.fep_uncertainty /
                                    r.full_energy_peak_efficiency);
    }

    EfficiencyTransfer et(geom, {0.0, 0.0, -d_ref}, anchor, MuChoice::Total);

    for (double d : {6.0 * a_ext, 20.0 * a_ext}) {
        const Eigen::Vector3d tgt(0.0, 0.0, -d);
        for (size_t i = 0; i < Es.size(); ++i) {
            const double eps_t = et.eps_at(Es[i], tgt);
            calc.set_point_source(tgt);
            SimulationConfig cfg;
            cfg.energy_keV = Es[i];
            cfg.termination.target_fep_rel_precision = 0.01;
            cfg.termination.max_events = 4000000;
            cfg.seed = 999;
            const EfficiencyResult r = calc.compute(cfg);
            const double z = (eps_t - r.full_energy_peak_efficiency) /
                             std::max(r.fep_uncertainty, 1e-30);
            BOOST_CHECK_LT(std::fabs(z), 3.0);
        }
    }
}
#endif  // CEELO_RUN_MC_TESTS
