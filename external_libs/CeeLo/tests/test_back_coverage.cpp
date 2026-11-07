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

/// @file test_back_coverage.cpp
/// @brief GenerationOptions::back_coverage: a response characterised beside and behind the detector
///        answers there like direct Monte Carlo, keeps every front value of a response without it,
///        and is stamped version 4.

#define BOOST_TEST_MODULE BackCoverageTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "io/ResponseKernel.h"
#include "materials/Material.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;

GeometryDescriptor small_nai() {
    GeometryDescriptor gd;
    gd.set_dimensions(CylinderDims{2.54, 5.08});
    gd.materials = {MaterialSpec::from(make_NaI()), MaterialSpec::from(make_Aluminum())};
    gd.crystal_material_index = 0;
    LayerSpec can;
    can.material_index = 1;
    can.front_thickness_cm = 0.05;
    can.side_thickness_cm = 0.05;
    can.z_end_cm = 5.08;
    gd.layers.push_back(can);
    gd.set_layer_back_thickness_cm(0, 0.1);   // a closed can: back structure described
    return gd;
}

GenerationOptions quick(const GenerationOptions::BackCoverage back) {
    GenerationOptions o;
    o.profile = ResponseProfile::FarField;
    o.node_fep_precision = 0.03;
    o.n_energy_scan = 10;
    o.n_energy_nodes = 8;
    o.n_shape_energies = 2;
    o.n_cos_theta_scan = 5;
    o.n_cos_theta_nodes = 5;
    o.e_min_keV = 100.0;
    o.e_max_keV = 1500.0;
    o.num_threads = 1;   // bit-exact across runs (multi-threaded batches race)
    o.base_seed = 77;
    o.back_coverage = back;
    return o;
}

// The General profile, so side & back coverage also measures the back near-field tables (stage 5).
//  The side & back nodes (stages 4 and 5) run at 0.5% here (1% after back_precision_scale), so the
//  far-field back value the near field is relative to is not 6% noise; everything else at 3%.
GenerationOptions general(const bool back_near, const double body_depth_cm = 0.0) {
    GenerationOptions o = quick(GenerationOptions::BackCoverage::Coarse);
    o.profile = ResponseProfile::General;
    o.back_near_field = back_near;
    o.back_body_depth_cm = body_depth_cm;
    o.node_precision = [](const uint32_t stage, double) { return (stage >= 4) ? 0.005 : 0.03; };
    return o;
}

// Direct Monte Carlo FEP at `pos` (general cone), the truth the tables are checked against.
EfficiencyResult direct_mc(const GeometryDescriptor& gd, const Eigen::Vector3d& pos,
                           const double energy_keV, const uint64_t seed) {
    std::vector<std::unique_ptr<Material>> owned;
    EfficiencyCalculator calc;
    ResponseGenerator::configure_calculator(calc, gd, owned);
    calc.enable_general_cone_sampling(true);
    calc.set_point_source(pos);
    SimulationConfig cfg;
    cfg.energy_keV = energy_keV;
    cfg.termination.target_fep_rel_precision = 0.005;
    cfg.termination.max_events = 40000000;
    cfg.termination.min_events = 20000;
    cfg.num_threads = 3;
    cfg.seed = seed;
    return calc.compute(cfg);
}

// `r` without its back near-field tables (and the body they assume), via an XML copy.
std::shared_ptr<DetectorResponse> without_back_near(const DetectorResponse& r) {
    std::shared_ptr<DetectorResponse> out = DetectorResponse::from_xml_string(r.to_xml_string());
    out->back_near_field = NearFieldModel{};
    out->tot_eff.back_near_field = NearFieldModel{};
    out->provenance.back_body_depth_cm = 0.0;
    out->finalize();
    return out;
}

// A node of a back near-field table: along the node's direction, at its clearance past the body.
Eigen::Vector3d node_position(const DetectorResponse& r, const NearFieldModel& m, const size_t ic,
                              const size_t id) {
    const double a = r.transverse_half_extent();
    return r.back_near_body().position_at_clearance(source_position(1.0, m.cos_thetas[ic]),
                                                    m.dists_cm[id] - a);
}

}  // namespace


BOOST_AUTO_TEST_CASE(back_coverage_front_unchanged_back_answered) {
    const GeometryDescriptor gd = small_nai();

    const std::shared_ptr<DetectorResponse> none = ResponseGenerator::generate(gd, quick(GenerationOptions::BackCoverage::None));
    const std::shared_ptr<DetectorResponse> covered = ResponseGenerator::generate(gd, quick(GenerationOptions::BackCoverage::Coarse));
    BOOST_REQUIRE(none && covered);

    // The choices are recorded, so a host can offer the same run again.
    BOOST_CHECK(none->provenance.back_coverage == BackCoverage::None);
    BOOST_CHECK(covered->provenance.back_coverage == BackCoverage::Coarse);
    BOOST_CHECK(covered->provenance.precision_profile == PrecisionProfile::Uniform);
    BOOST_CHECK(!covered->provenance.n_anchor_angles);

    // The plan counts the back nodes.
    BOOST_CHECK_EQUAL(ResponseGenerator::plan_nodes(gd, quick(GenerationOptions::BackCoverage::Coarse)).n_back(),
                      2 * 6);

    // Covered: the eta table reaches behind the face plane.
    BOOST_REQUIRE(covered->eta_fep.cos_thetas.front() < -0.99);
    BOOST_CHECK(!(none->eta_fep.cos_thetas.front() < 0.0));
    // The far-field profile measures no near field, in front or behind.
    BOOST_CHECK(covered->back_near_field.empty() && covered->tot_eff.back_near_field.empty());
    // Every stored back value is a measurement or held from one - never a log of zero counts; a held
    //  node carries a 100% sigma on top of the on-axis one, so at most sqrt(2).
    for (size_t i = 0; i < covered->eta_fep.ln_eta.size(); ++i) {
        BOOST_CHECK(std::isfinite(covered->eta_fep.ln_eta[i]) && (covered->eta_fep.ln_eta[i] > -20.0));
        BOOST_CHECK(covered->eta_fep.frac_sigma[i] <= std::sqrt(2.0));
    }

    // Every front value is what the uncovered response gives - same nodes, same seeds, and the
    //  front interpolated on its own - with its uncertainty and flag; also below the grazing node
    //  (cos theta 0.02, about 88.85 deg), where both clamp to it.
    for (const double theta_deg : {0.0, 30.0, 60.0, 85.0, 88.0, 89.5, 90.0}) {
        for (const double d : {5.0, 30.0}) {
            const Eigen::Vector3d pos = source_position(d, std::cos(theta_deg * kPi / 180.0));
            for (const double e : {120.0, 662.0, 1332.0}) {
                const EffResult a = none->eps_fep_at(e, pos);
                const EffResult b = covered->eps_fep_at(e, pos);
                BOOST_CHECK_EQUAL(a.value, b.value);
                BOOST_CHECK_EQUAL(a.sigma, b.sigma);
                BOOST_CHECK(a.flag == b.flag);
            }
        }
    }

    // Behind: a back node position against direct Monte Carlo (general cone), and no refuse flag.
    {
        const double a = gd.transverse_half_extent();
        const double d_far = std::max(10.0 * a, 10.0);
        const Eigen::Vector3d pos = source_position(d_far, std::cos(135.0 * kPi / 180.0));
        std::vector<std::unique_ptr<Material>> owned;
        EfficiencyCalculator calc;
        ResponseGenerator::configure_calculator(calc, gd, owned);
        calc.enable_general_cone_sampling(true);
        calc.set_point_source(pos);
        SimulationConfig cfg;
        cfg.energy_keV = 662.0;
        cfg.termination.target_fep_rel_precision = 0.005;
        cfg.termination.max_events = 40000000;
        cfg.termination.min_events = 20000;
        cfg.num_threads = 3;
        cfg.seed = 99;
        const EfficiencyResult mc = calc.compute(cfg);
        const EffResult m = covered->eps_fep_at(662.0, pos);
        const double rel = m.value / mc.full_energy_peak_efficiency - 1.0;
        BOOST_TEST_MESSAGE("behind (135 deg, far field) @ 662 keV: model/MC - 1 = " << 100.0 * rel << "%");
        BOOST_CHECK_LT(std::abs(rel), 0.05);
        BOOST_CHECK(m.flag != ResponseFlag::NeedsMc);
        BOOST_CHECK(m.flag != ResponseFlag::OutOfRangeClamped);
    }

    // Stamped v4, and it round-trips.
    const std::string xml = covered->to_xml_string();
    BOOST_CHECK(xml.find("version=\"4\"") != std::string::npos);
    BOOST_CHECK(none->to_xml_string().find("version=\"4\"") == std::string::npos);
    const std::shared_ptr<DetectorResponse> back = DetectorResponse::from_xml_string(xml);
    BOOST_REQUIRE(back);
    const Eigen::Vector3d behind = source_position(40.0, std::cos(150.0 * kPi / 180.0));
    BOOST_CHECK_EQUAL(back->eps_fep_at(662.0, behind).value, covered->eps_fep_at(662.0, behind).value);
    BOOST_CHECK_EQUAL(back->content_hash(), covered->content_hash());
}


// The back near-field tables (General / Contact profile with side & back coverage): measured on their
//  own stage, so every other table is what a run without them gives; they reproduce their own nodes
//  on the coordinate the evaluator uses; and close beside the detector they bring the response to
//  direct Monte Carlo.
BOOST_AUTO_TEST_CASE(back_near_field_tables) {
    const GeometryDescriptor gd = small_nai();
    const double a = gd.transverse_half_extent();

    const std::shared_ptr<DetectorResponse> with = ResponseGenerator::generate(gd, general(true));
    const std::shared_ptr<DetectorResponse> without = ResponseGenerator::generate(gd, general(false));
    BOOST_REQUIRE(with && without);

    // Built only when asked for, and only with a near-field profile.
    BOOST_REQUIRE(!with->back_near_field.empty());
    BOOST_REQUIRE(!with->tot_eff.back_near_field.empty());
    BOOST_CHECK(without->back_near_field.empty() && without->tot_eff.back_near_field.empty());
    const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, general(true));
    BOOST_CHECK_EQUAL(plan.n_back_near(), 2 * 6 * 4);
    BOOST_CHECK_EQUAL(ResponseGenerator::plan_nodes(gd, general(false)).n_back_near(), 0);
    BOOST_CHECK_EQUAL(ResponseGenerator::plan_nodes(gd, quick(GenerationOptions::BackCoverage::Coarse)).n_back_near(), 0);
    BOOST_CHECK_EQUAL(ResponseGenerator::estimated_node_count(gd, general(true)),
                      ResponseGenerator::estimated_node_count(gd, general(false)) + 2 * 6 * 4);

    // Behind the face plane only: 6 angles, 4 clearances and the ln N = 0 anchor row.
    const NearFieldModel& bn = with->back_near_field;
    BOOST_REQUIRE_EQUAL(bn.cos_thetas.size(), 6u);
    BOOST_REQUIRE_EQUAL(bn.dists_cm.size(), 5u);
    BOOST_CHECK(bn.cos_thetas.back() < 0.0);
    BOOST_CHECK_CLOSE(bn.dists_cm.front(), 1.05 * a, 1e-9);
    for (size_t ie = 0; ie < bn.energies_keV.size(); ++ie)
        for (size_t ic = 0; ic < bn.cos_thetas.size(); ++ic) {
            BOOST_CHECK_EQUAL(bn.ln_n[bn.index(ie, ic, bn.dists_cm.size() - 1)], 0.0);
            for (size_t id = 0; id < bn.dists_cm.size(); ++id) {
                BOOST_CHECK(std::isfinite(bn.ln_n[bn.index(ie, ic, id)]));
                BOOST_CHECK(bn.frac_sigma[bn.index(ie, ic, id)] <= 1.0);
            }
        }

    // Their own stage: every other table, node for node, is the run's without them.
    BOOST_CHECK(with->eta_fep.cos_thetas == without->eta_fep.cos_thetas);
    BOOST_CHECK(with->eta_fep.ln_eta == without->eta_fep.ln_eta);
    BOOST_CHECK(with->eta_fep.frac_sigma == without->eta_fep.frac_sigma);
    BOOST_CHECK(with->near_field.ln_n == without->near_field.ln_n);
    BOOST_CHECK(with->near_field.frac_sigma == without->near_field.frac_sigma);
    BOOST_CHECK(with->tot_eff.tier == without->tot_eff.tier);
    BOOST_CHECK(with->tot_eff.eta_tot.ln_eta == without->tot_eff.eta_tot.ln_eta);
    BOOST_CHECK(with->tot_eff.near_field.ln_n == without->tot_eff.near_field.ln_n);

    // In front of the face plane nothing reads them.
    for (const double theta_deg : {0.0, 45.0, 80.0, 88.0, 90.0}) {
        for (const double d : {3.0, 8.0, 30.0}) {
            const Eigen::Vector3d pos = source_position(d, std::cos(theta_deg * kPi / 180.0));
            for (const double e : {120.0, 662.0, 1332.0}) {
                const EffResult x = with->eps_fep_at(e, pos), y = without->eps_fep_at(e, pos);
                BOOST_CHECK_EQUAL(x.value, y.value);
                BOOST_CHECK_EQUAL(x.sigma, y.sigma);
                BOOST_CHECK(x.flag == y.flag);
                const EffResult xt = with->eps_total_at(e, pos), yt = without->eps_total_at(e, pos);
                BOOST_CHECK_EQUAL(xt.value, yt.value);
                BOOST_CHECK_EQUAL(xt.sigma, yt.sigma);
            }
        }
    }

    // Each measured node is reproduced on the evaluator's own coordinate: model / (eta * K) = N.
    for (size_t ic = 0; ic < bn.cos_thetas.size(); ++ic) {
        for (size_t id = 0; id + 1 < bn.dists_cm.size(); ++id) {
            const Eigen::Vector3d pos = node_position(*with, bn, ic, id);
            BOOST_CHECK_GT(pos.z(), 0.0);
            const ApertureQuadrature q = with->make_quadrature(pos);
            for (size_t ie = 0; ie < bn.energies_keV.size(); ++ie) {
                const size_t idx = bn.index(ie, ic, id);
                if (bn.frac_sigma[idx] >= 1.0)
                    continue;   // held, not measured
                const double e = bn.energies_keV[ie];
                bool clamped = false;
                const double eta = std::exp(with->eta_fep.eval_ln(e, bn.cos_thetas[ic], 0.0, clamped));
                const double ratio = with->eps_fep_at(e, pos, q).value
                                     / (eta * with->kernel_K(e, q, MuChoice::Total));
                BOOST_CHECK_CLOSE(ratio, std::exp(bn.ln_n[idx]), 1e-6);
                BOOST_CHECK(with->eps_fep_at(e, pos, q).flag != ResponseFlag::NearFieldUnmodeled);
            }
        }
    }

    // Past the anchor row the tables are out of their gate - the value is the one without them.
    {
        const Eigen::Vector3d far = with->back_near_body().position_at_clearance(
            source_position(1.0, -1.0), bn.dists_cm.back() - a + 1.0);
        for (const double e : {120.0, 662.0}) {
            BOOST_CHECK_EQUAL(with->eps_fep_at(e, far).value, without->eps_fep_at(e, far).value);
            BOOST_CHECK_EQUAL(with->eps_total_at(e, far).value, without->eps_total_at(e, far).value);
        }
    }

    // Between the face plane and the first back node (93 deg) the grazing column holds.
    {
        const Eigen::Vector3d pos = source_position(6.0, std::cos(91.0 * kPi / 180.0));
        const EffResult r = with->eps_fep_at(662.0, pos);
        BOOST_CHECK(std::isfinite(r.value) && (r.value > 0.0));
    }

    // Off-node, close beside the front of the detector (between the 93 and 100 deg columns and the
    //  0.05a and 0.3a rows) at a shape energy, where the kernel alone reads ~5% high: the table takes
    //  the response to direct MC.
    {
        const double e = bn.energies_keV.back();
        const Eigen::Vector3d pos = with->back_near_body().position_at_clearance(
            source_position(1.0, std::cos(96.5 * kPi / 180.0)), 0.15 * a);
        const EfficiencyResult mc = direct_mc(gd, pos, e, 1234);
        const double rel_with = with->eps_fep_at(e, pos).value / mc.full_energy_peak_efficiency - 1.0;
        const double rel_without = without->eps_fep_at(e, pos).value / mc.full_energy_peak_efficiency - 1.0;
        BOOST_TEST_MESSAGE("beside (96.5 deg, 0.15a) @ " << e << " keV: with table " << 100.0 * rel_with
                           << "%, without " << 100.0 * rel_without << "%");
        BOOST_CHECK_LT(std::abs(rel_with), 0.03);
        BOOST_CHECK_LT(std::abs(rel_with), std::abs(rel_without));
    }

    // Stored only when present; the copy evaluates the same.
    const std::string xml = with->to_xml_string();
    BOOST_CHECK(xml.find("<BackNearField") != std::string::npos);
    BOOST_CHECK(without->to_xml_string().find("<BackNearField") == std::string::npos);
    const std::shared_ptr<DetectorResponse> copy = DetectorResponse::from_xml_string(xml);
    BOOST_CHECK_EQUAL(copy->content_hash(), with->content_hash());
    const Eigen::Vector3d beside = source_position(4.0, std::cos(120.0 * kPi / 180.0));
    BOOST_CHECK_EQUAL(copy->eps_fep_at(662.0, beside).value, with->eps_fep_at(662.0, beside).value);
    BOOST_CHECK_EQUAL(copy->eps_total_at(662.0, beside).value, with->eps_total_at(662.0, beside).value);
}


// The assumed detector body (GenerationOptions::back_body_depth_cm): no back near node sits inside
//  it, and a source inside it - outside the described detector, where a real one keeps its cryostat
//  or electronics - is answered exactly as by a response without the tables.
BOOST_AUTO_TEST_CASE(back_near_field_body) {
    const GeometryDescriptor gd = small_nai();
    const double a = gd.transverse_half_extent();
    const double depth = 15.0;

    const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, general(true, depth));
    BOOST_CHECK_EQUAL(plan.back_body_depth_cm, depth);
    BOOST_CHECK_EQUAL(plan.n_back_near(), 2 * 6 * 4);

    // A detector so small that its own far-field back node falls inside 15 cm keeps its housing as
    //  described (no MC needed to tell).
    {
        GeometryDescriptor tiny;
        tiny.set_dimensions(CylinderDims{0.5, 1.0});
        tiny.materials = {MaterialSpec::from(make_NaI())};
        tiny.crystal_material_index = 0;
        BOOST_CHECK_EQUAL(ResponseGenerator::plan_nodes(tiny, general(true, depth)).back_body_depth_cm, 0.0);
    }

    // A crystal whose back rim falls between the angle nodes gets a column there (6.5 x 5 cm: the rim
    //  is at 147 deg, between 135 and 155); the 1 x 2 in NaI's rim (153 deg) is within 4 of 155.
    {
        GeometryDescriptor hpge;
        hpge.set_dimensions(CylinderDims{3.25, 5.0});
        hpge.materials = {MaterialSpec::from(make_NaI())};
        hpge.crystal_material_index = 0;
        BOOST_CHECK_EQUAL(ResponseGenerator::plan_nodes(hpge, general(true)).n_back_near_positions, 7 * 4);
        BOOST_CHECK_EQUAL(ResponseGenerator::plan_nodes(gd, general(true)).n_back_near_positions, 6 * 4);
    }

    const std::shared_ptr<DetectorResponse> r = ResponseGenerator::generate(gd, general(true, depth));
    BOOST_REQUIRE(r && !r->back_near_field.empty());
    BOOST_CHECK_EQUAL(r->provenance.back_body_depth_cm, depth);
    BOOST_CHECK_EQUAL(r->back_near_body().z_max_cm, depth);
    const std::shared_ptr<DetectorResponse> stripped = without_back_near(*r);

    // Directly behind, the nodes start at the body's back face.
    const NearFieldModel& bn = r->back_near_field;
    for (size_t id = 0; id + 1 < bn.dists_cm.size(); ++id) {
        const Eigen::Vector3d pos = node_position(*r, bn, 0, id);   // 180 deg
        BOOST_CHECK_CLOSE(pos.z(), depth + bn.dists_cm[id] - a, 1e-6);
    }

    // Inside the body but outside the described detector (its back is at 5.18 cm): no table.
    for (const Eigen::Vector3d& pos : {Eigen::Vector3d(0.0, 0.0, 8.0), Eigen::Vector3d(1.0, 0.5, 12.0)}) {
        for (const double e : {120.0, 662.0}) {
            const EffResult x = r->eps_fep_at(e, pos), y = stripped->eps_fep_at(e, pos);
            BOOST_CHECK_EQUAL(x.value, y.value);
            BOOST_CHECK_EQUAL(x.sigma, y.sigma);
            BOOST_CHECK(x.flag == y.flag);
            BOOST_CHECK_EQUAL(r->eps_total_at(e, pos).value, stripped->eps_total_at(e, pos).value);
        }
    }

    // Beside the body (radially outside it, behind the crystal) the table applies.
    const Eigen::Vector3d beside_body(a + 1.0, 0.0, 10.0);
    BOOST_CHECK(r->eps_fep_at(662.0, beside_body).flag != ResponseFlag::NearFieldUnmodeled);

    // The body depth is stored with the tables and survives the round trip.
    const std::shared_ptr<DetectorResponse> copy = DetectorResponse::from_xml_string(r->to_xml_string());
    BOOST_CHECK_EQUAL(copy->provenance.back_body_depth_cm, depth);
    BOOST_CHECK_EQUAL(copy->content_hash(), r->content_hash());
    BOOST_CHECK_EQUAL(copy->eps_fep_at(662.0, beside_body).value, r->eps_fep_at(662.0, beside_body).value);
}
