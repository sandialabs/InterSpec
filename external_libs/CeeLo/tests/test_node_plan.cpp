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
#include <string>
#include <vector>

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
