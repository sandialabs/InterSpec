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

#define BOOST_TEST_MODULE NodePlanTests
#include <boost/test/unit_test.hpp>

// ResponseGenerator::plan_nodes() / backbone_scan_energies() (pure, no MC) and
// the per-node GenerationOptions::node_progress callback (a tiny real MC run:
// far-field profile, a handful of nodes at a few thousand events each).

#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "materials/Material.h"

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

using namespace ceelo;

namespace {

GeometryDescriptor small_nai() {
    GeometryDescriptor gd;
    gd.set_dimensions(CylinderDims{1.27, 2.54});  // NaI 1"x1"
    gd.materials = {MaterialSpec::from(make_NaI()),
                    MaterialSpec::from(make_Aluminum())};
    gd.crystal_material_index = 0;
    LayerSpec can;
    can.material_index = 1;
    can.front_thickness_cm = 0.05;
    can.side_thickness_cm = 0.05;
    can.z_end_cm = 2.54;
    gd.layers.push_back(can);
    return gd;
}

GeometryDescriptor small_box() {
    GeometryDescriptor gd;
    gd.set_dimensions(BoxDims{1.0, 1.5, 2.0});
    gd.symmetry = ResponseSymmetry::Quadrant;
    gd.materials = {MaterialSpec::from(make_NaI())};
    gd.crystal_material_index = 0;
    return gd;
}

GenerationOptions tiny_options() {
    GenerationOptions o;
    o.profile = ResponseProfile::FarField;
    o.node_fep_precision = 0.05;
    o.min_events_per_node = 1000;
    o.max_events_per_node = 4000;
    o.max_cpu_seconds_per_node = 2.0;
    o.n_energy_scan = 6;
    o.n_energy_nodes = 4;
    o.n_cos_theta_scan = 3;
    o.n_cos_theta_nodes = 3;
    o.n_shape_energies = 3;
    o.base_seed = 11;
    o.detector_name = "node plan test";
    return o;
}

}  // namespace

BOOST_AUTO_TEST_CASE(plan_matches_estimate_and_stage_structure) {
    const GeometryDescriptor gd = small_nai();

    {   // General profile (the defaults)
        const GenerationOptions o;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, o);
        BOOST_CHECK_EQUAL(plan.total(), ResponseGenerator::estimated_node_count(gd, o));
        BOOST_CHECK(plan.backbone_energies_keV ==
                    ResponseGenerator::backbone_scan_energies(gd, o));
        BOOST_CHECK(plan.n_backbone() >= o.n_energy_scan);  // grid + K-edge flanks
        BOOST_CHECK_EQUAL(static_cast<int>(plan.shape_energies_keV.size()),
                          o.n_shape_energies);
        BOOST_CHECK_EQUAL(plan.n_cos_theta, o.n_cos_theta_scan);
        BOOST_CHECK_EQUAL(plan.n_phi, 1);
        BOOST_CHECK_EQUAL(plan.n_near_positions, 9 * 8);
        BOOST_CHECK(!plan.transfer_anchors);
        BOOST_CHECK_EQUAL(plan.total(),
                          plan.n_backbone() + o.n_shape_energies * o.n_cos_theta_scan +
                              o.n_shape_energies * 72);

        // Backbone grid: strictly increasing, spanning the requested range (the
        // K-edge flanks sit 0.1% either side of an edge, always inside it).
        for (size_t i = 1; i < plan.backbone_energies_keV.size(); ++i)
            BOOST_CHECK(plan.backbone_energies_keV[i] > plan.backbone_energies_keV[i - 1]);
        BOOST_CHECK_CLOSE(plan.backbone_energies_keV.front(), o.e_min_keV, 1e-9);
        BOOST_CHECK_CLOSE(plan.backbone_energies_keV.back(), o.e_max_keV, 1e-9);
        for (size_t i = 1; i < plan.shape_energies_keV.size(); ++i)
            BOOST_CHECK(plan.shape_energies_keV[i] > plan.shape_energies_keV[i - 1]);
    }

    {   // Contact profile: one more MC distance row
        GenerationOptions o;
        o.profile = ResponseProfile::Contact;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, o);
        BOOST_CHECK_EQUAL(plan.n_near_positions, 9 * 9);
        BOOST_CHECK_EQUAL(plan.total(), ResponseGenerator::estimated_node_count(gd, o));
    }

    {   // Far-field profile: no near-field tensor at all
        GenerationOptions o;
        o.profile = ResponseProfile::FarField;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, o);
        BOOST_CHECK_EQUAL(plan.n_near_positions, 0);
        BOOST_CHECK_EQUAL(plan.n_near(), 0);
        BOOST_CHECK(plan.n_angular() > 0);
        BOOST_CHECK_EQUAL(plan.total(), ResponseGenerator::estimated_node_count(gd, o));
    }

    {   // Flat transfer: the backbone only
        GenerationOptions o;
        o.transfer_mode = true;
        o.n_anchor_angles = 1;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, o);
        BOOST_CHECK(plan.transfer_anchors);
        BOOST_CHECK_EQUAL(plan.n_cos_theta, 0);
        BOOST_CHECK_EQUAL(plan.n_angular(), 0);
        BOOST_CHECK_EQUAL(plan.n_near(), 0);
        BOOST_CHECK_EQUAL(plan.total(), plan.n_backbone());
        BOOST_CHECK_EQUAL(plan.total(), ResponseGenerator::estimated_node_count(gd, o));
    }

    {   // Transfer with angle anchors: 3 forced angles per shape energy
        GenerationOptions o;
        o.transfer_mode = true;
        o.n_anchor_angles = 3;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, o);
        BOOST_CHECK(plan.transfer_anchors);
        BOOST_CHECK_EQUAL(plan.n_cos_theta, 3);
        BOOST_CHECK_EQUAL(plan.n_near(), 0);
        BOOST_CHECK_EQUAL(plan.total(), plan.n_backbone() + 3 * o.n_shape_energies);
        BOOST_CHECK_EQUAL(plan.total(), ResponseGenerator::estimated_node_count(gd, o));
    }

    {   // Boxes scan azimuths too
        const GenerationOptions o;
        const ResponseGenerator::NodePlan plan =
            ResponseGenerator::plan_nodes(small_box(), o);
        BOOST_CHECK_EQUAL(plan.n_phi, o.n_phi_nodes);
        BOOST_CHECK_EQUAL(plan.n_angular(),
                          o.n_shape_energies * o.n_cos_theta_scan * o.n_phi_nodes);
        BOOST_CHECK_EQUAL(plan.total(),
                          ResponseGenerator::estimated_node_count(small_box(), o));
    }
}

BOOST_AUTO_TEST_CASE(node_progress_reports_every_node) {
    const GeometryDescriptor gd = small_nai();
    GenerationOptions opts = tiny_options();
    const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, opts);

    std::vector<NodeProgress> seen;
    opts.node_progress = [&](const NodeProgress& p) { seen.push_back(p); };
    int legacy_calls = 0;
    opts.progress = [&](double, const std::string&) { ++legacy_calls; };

    const std::shared_ptr<DetectorResponse> resp = ResponseGenerator::generate(gd, opts);
    BOOST_REQUIRE(resp);
    BOOST_REQUIRE(!seen.empty());

    // Every MC node reported exactly once, counting from 1, never more than the
    // plan (shape energies that snap to the same backbone node merge), and the
    // legacy callback fired at least as often (it also announces "Done").
    BOOST_CHECK(static_cast<int>(seen.size()) <= plan.total());
    BOOST_CHECK(static_cast<int>(seen.size()) >= plan.n_backbone());
    BOOST_CHECK(legacy_calls >= static_cast<int>(seen.size()));

    for (size_t i = 0; i < seen.size(); ++i) {
        BOOST_CHECK_EQUAL(seen[i].nodes_done, static_cast<int>(i) + 1);
        BOOST_CHECK_EQUAL(seen[i].nodes_total, plan.total());
        BOOST_CHECK(seen[i].stage >= 1u);
        BOOST_CHECK(seen[i].stage <= 3u);
        if (i > 0)
            BOOST_CHECK(seen[i].stage >= seen[i - 1].stage);
        BOOST_CHECK(seen[i].events > 0u);
        BOOST_CHECK(seen[i].energy_keV > 0.0);
        BOOST_CHECK(seen[i].cpu_s >= 0.0);
        BOOST_CHECK(seen[i].wall_s >= 0.0);
    }

    // The first n_backbone nodes are the backbone grid, in order, at stage 1.
    for (int i = 0; i < plan.n_backbone(); ++i) {
        BOOST_CHECK_EQUAL(seen[static_cast<size_t>(i)].stage, 1u);
        BOOST_CHECK_CLOSE(seen[static_cast<size_t>(i)].energy_keV,
                          plan.backbone_energies_keV[static_cast<size_t>(i)], 1e-9);
    }
}

// Degenerate grid options. No current caller supplies any of these -- InterSpec and
// make_golden_response both leave every field below at its default -- which is exactly
// why they are worth pinning: nothing else in the tree would notice a regression. Each
// one used to fail silently rather than loudly: a NaN grid handed to std::sort, a node
// count the run never matched, or a "successful" response that was quietly angle-flat.

BOOST_AUTO_TEST_CASE(degenerate_options_give_a_sane_plan) {
    const GeometryDescriptor gd = small_nai();

    {   // A one-point energy scan is just the low end of the range -- never 0/0.
        GenerationOptions o;
        o.n_energy_scan = 1;
        const std::vector<double> scan = ResponseGenerator::backbone_scan_energies(gd, o);
        BOOST_REQUIRE(!scan.empty());
        for (const double e : scan) BOOST_CHECK(std::isfinite(e));
        for (size_t i = 1; i < scan.size(); ++i) BOOST_CHECK(scan[i] > scan[i - 1]);
    }

    {   // The default grid stays finite and ascending (the sort below it has no
        // defined behaviour once a NaN reaches it).
        const GenerationOptions o;
        const std::vector<double> scan = ResponseGenerator::backbone_scan_energies(gd, o);
        BOOST_REQUIRE(!scan.empty());
        for (const double e : scan) BOOST_CHECK(std::isfinite(e));
        for (size_t i = 1; i < scan.size(); ++i) BOOST_CHECK(scan[i] > scan[i - 1]);
    }

    // Too few scan angles: the plan must never report a negative node count, and
    // must never claim more angular nodes than the run would execute (the scan
    // loop runs zero times for n <= 0).
    for (const int n : {1, 0, -4}) {
        GenerationOptions o;
        o.n_cos_theta_scan = n;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, o);
        BOOST_CHECK(plan.n_cos_theta >= 0);
        BOOST_CHECK(plan.n_angular() >= 0);
        BOOST_CHECK(plan.n_backbone() > 0);
        BOOST_CHECK(plan.total() >= plan.n_backbone());
        BOOST_CHECK_EQUAL(plan.n_cos_theta, std::max(0, n));
        BOOST_CHECK_EQUAL(plan.total(), ResponseGenerator::estimated_node_count(gd, o));
    }

    // A non-positive or inverted energy range: e_max/e_min is inf, so every grid
    // point past the first is 0 * inf = NaN. Rejected rather than sorted.
    std::vector<GenerationOptions> bad(4);
    bad[0].e_min_keV = 0.0;
    bad[1].e_min_keV = -5.0;
    bad[2].e_min_keV = std::numeric_limits<double>::quiet_NaN();
    bad[3].e_max_keV = bad[3].e_min_keV;
    for (const GenerationOptions& o : bad) {
        BOOST_CHECK_THROW(ResponseGenerator::backbone_scan_energies(gd, o), std::runtime_error);
        BOOST_CHECK_THROW(ResponseGenerator::plan_nodes(gd, o), std::runtime_error);
        BOOST_CHECK_THROW(ResponseGenerator::estimated_node_count(gd, o), std::runtime_error);
    }
}

// generate() refuses the same inputs, and does so before any MC runs. The rejected
// cases are therefore free; the one accepted case below is a real (tiny) backbone run.
BOOST_AUTO_TEST_CASE(degenerate_options_rejected_by_generate) {
    const GeometryDescriptor gd = small_nai();

    // Fewer than two scan angles. One point sits at the grazing cutoff, which is
    // the very point shape_at() normalizes against, so the response would come
    // back angle-flat at the on-axis value; zero points index an empty scan_ct.
    for (const int n : {1, 0, -4}) {
        GenerationOptions o = tiny_options();
        o.n_cos_theta_scan = n;
        BOOST_CHECK_THROW(ResponseGenerator::generate(gd, o), std::runtime_error);
    }

    // ... but the angle-flat transfer variant never runs an angular scan, so it is still
    // allowed to leave n_cos_theta_scan alone. This one is not free: it runs the stage-1
    // backbone MC (tiny_options, so ~0.1 s).
    {
        GenerationOptions o = tiny_options();
        o.transfer_mode = true;
        o.n_anchor_angles = 1;
        o.n_cos_theta_scan = 1;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(gd, o);
        BOOST_CHECK_EQUAL(plan.n_cos_theta, 0);
        BOOST_CHECK_NO_THROW(ResponseGenerator::generate(gd, o));
    }

    std::vector<GenerationOptions> bad(4);
    for (GenerationOptions& o : bad) o = tiny_options();
    bad[0].e_min_keV = 0.0;
    bad[1].e_min_keV = -5.0;
    bad[2].e_min_keV = std::numeric_limits<double>::quiet_NaN();
    bad[3].e_max_keV = bad[3].e_min_keV;
    for (const GenerationOptions& o : bad)
        BOOST_CHECK_THROW(ResponseGenerator::generate(gd, o), std::runtime_error);

    // The probe banks build their own log-spaced energies instead of calling
    // backbone_scan_energies, so they carry the same range check separately. Left
    // unguarded, MC at a NaN energy never returns (the rejection sampler loops).
    for (const GenerationOptions& o : bad) {
        BOOST_CHECK_THROW(ResponseGenerator::probe_bank(gd, o, 3, 7000), std::runtime_error);
        const DetectorResponse empty_resp;
        BOOST_CHECK_THROW(
            ResponseGenerator::plan_structured_probes(gd, o, empty_resp,
                                                      probe_family_bit(ProbeFamily::Random),
                                                      2, 7000),
            std::runtime_error);
    }
}

// The rest of the grid options. Each of these used to abort the process outright
// (an empty vector reached a `n - 1` size_t underflow), except cos_theta_min, which
// reached Pchip as a non-ascending axis and blamed Pchip rather than the option.
BOOST_AUTO_TEST_CASE(other_degenerate_options_rejected) {
    const GeometryDescriptor gd = small_nai();

    for (const int n : {0, -3}) {
        GenerationOptions o = tiny_options();
        o.n_energy_scan = n;
        BOOST_CHECK_THROW(ResponseGenerator::backbone_scan_energies(gd, o), std::runtime_error);
        BOOST_CHECK_THROW(ResponseGenerator::plan_nodes(gd, o), std::runtime_error);
        BOOST_CHECK_THROW(ResponseGenerator::generate(gd, o), std::runtime_error);
    }

    for (const int n : {0, -3}) {
        GenerationOptions o = tiny_options();
        o.n_shape_energies = n;
        BOOST_CHECK_THROW(ResponseGenerator::generate(gd, o), std::runtime_error);
    }

    for (const double ct : {1.0, 1.5}) {
        GenerationOptions o = tiny_options();
        o.cos_theta_min = ct;
        BOOST_CHECK_THROW(ResponseGenerator::generate(gd, o), std::runtime_error);
    }

    // A box with no azimuths is clamped to one, not rejected: plan_nodes() already
    // promises max(1, n_phi_nodes), so the run has to agree with it.
    for (const int n : {0, -3}) {
        GenerationOptions o = tiny_options();
        o.n_phi_nodes = n;
        const ResponseGenerator::NodePlan plan = ResponseGenerator::plan_nodes(small_box(), o);
        BOOST_CHECK_EQUAL(plan.n_phi, 1);
        const std::shared_ptr<DetectorResponse> resp =
            ResponseGenerator::generate(small_box(), o);
        BOOST_REQUIRE(resp);
        BOOST_CHECK_EQUAL(plan.total(), ResponseGenerator::estimated_node_count(small_box(), o));
    }

    // One scan energy remains legal all the way through a real run.
    {
        GenerationOptions o = tiny_options();
        o.n_energy_scan = 1;
        BOOST_CHECK_NO_THROW(ResponseGenerator::generate(gd, o));
    }
}
