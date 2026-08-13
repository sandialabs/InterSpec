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

#define BOOST_TEST_MODULE StructuredProbeBankTests
#include <boost/test/unit_test.hpp>

// Pure (no-MC) tests of the D-b structured probe family planner: correct tagged
// family counts, midpoints landing inside the node gaps, K-edge segments never
// bridged. plan_structured_probes reads only the response's node grids, so a
// hand-built response is enough -- no MC, no golden file.

#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "materials/Material.h"

#include <algorithm>
#include <cmath>

using namespace ceelo;

namespace {

GeometryDescriptor nai_descriptor() {
    GeometryDescriptor gd;
    gd.shape = DetectorShape::Cylinder;
    gd.dimensions_cm = {3.81, 7.62};  // NaI 3"x3"
    gd.materials = {MaterialSpec::from(make_NaI()),
                    MaterialSpec::from(make_Aluminum())};
    gd.crystal_material_index = 0;
    return gd;
}

// Build a response with explicit node grids (no MC, no finalize needed: the
// planner only reads coordinate vectors + edges).
std::shared_ptr<DetectorResponse> make_response(
    const std::vector<double>& energies, const std::vector<double>& cts,
    const std::vector<double>& edges, const std::vector<double>& near_E,
    const std::vector<double>& near_ct, const std::vector<double>& near_d) {
    auto resp = std::make_shared<DetectorResponse>();
    resp->descriptor = nai_descriptor();
    EtaTable& t = resp->eta_fep;
    t.energies_keV = energies;
    t.cos_thetas = cts;
    t.edges_keV = edges;
    t.ln_eta.assign(energies.size() * cts.size(), 0.0);
    t.frac_sigma.assign(t.ln_eta.size(), 0.003);
    if (!near_E.empty()) {
        NearFieldModel& nf = resp->near_field;
        nf.energies_keV = near_E;
        nf.cos_thetas = near_ct;
        nf.dists_cm = near_d;
        nf.ln_n.assign(near_E.size() * near_ct.size() * near_d.size(), 0.0);
        nf.frac_sigma.assign(nf.ln_n.size(), 0.003);
    }
    return resp;
}

int count_family(const std::vector<ResponseGenerator::StructuredProbe>& v,
                 ProbeFamily f) {
    return static_cast<int>(std::count_if(
        v.begin(), v.end(),
        [f](const ResponseGenerator::StructuredProbe& p) { return p.family == f; }));
}

}  // namespace

BOOST_AUTO_TEST_CASE(family_counts_and_coordinates) {
    const std::vector<double> E = {35, 60, 122, 344, 662, 1332, 3000};  // 6 gaps
    const std::vector<double> CT = {0.02, 0.3, 0.6, 1.0};               // 3 gaps
    const std::vector<double> NE = {35, 122, 662};                      // 2 gaps
    const std::vector<double> NC = {0.02, 0.3, 1.0};                    // 2 gaps
    const std::vector<double> ND = {1.0, 3.0, 10.0, 40.0};              // 3 gaps
    std::shared_ptr<DetectorResponse> resp = make_response(E, CT, {}, NE, NC, ND);

    GenerationOptions opts;
    opts.node_fep_precision = 0.003;
    const int n = 5;
    const std::vector<ResponseGenerator::StructuredProbe> plan =
        ResponseGenerator::plan_structured_probes(
            resp->descriptor, opts, *resp, kAllProbeFamilies, n, 20000);

    // Counts: capped by both n and the number of gaps.
    BOOST_CHECK_EQUAL(count_family(plan, ProbeFamily::EnergyGap), 5);  // min(5,6)
    BOOST_CHECK_EQUAL(count_family(plan, ProbeFamily::CtGap), 3);      // min(5,3)
    BOOST_CHECK_EQUAL(count_family(plan, ProbeFamily::Random), 5);
    BOOST_CHECK_GT(count_family(plan, ProbeFamily::ShapeEGap), 0);
    BOOST_CHECK_GT(count_family(plan, ProbeFamily::NearCell), 0);

    const double a = resp->descriptor.transverse_half_extent();
    const double d_far = std::max(opts.far_distance_a * a, 10.0);

    for (const ResponseGenerator::StructuredProbe& p : plan) {
        if (p.family == ProbeFamily::EnergyGap) {
            // On-axis, far, node precision, and E is a geometric-mean midpoint
            // strictly inside some adjacent-node gap.
            BOOST_CHECK_CLOSE(p.cos_theta, 1.0, 1e-9);
            BOOST_CHECK_CLOSE(p.d_cm, d_far, 1e-9);
            BOOST_CHECK_LE(p.node_precision, 0.003 + 1e-12);
            bool in_a_gap = false;
            for (size_t i = 0; i + 1 < E.size(); ++i)
                if (p.energy_keV > E[i] && p.energy_keV < E[i + 1]) in_a_gap = true;
            BOOST_CHECK(in_a_gap);
        }
        if (p.family == ProbeFamily::CtGap) {
            BOOST_CHECK_CLOSE(p.d_cm, d_far, 1e-9);
            bool in_a_gap = false;
            for (size_t i = 0; i + 1 < CT.size(); ++i)
                if (p.cos_theta > CT[i] && p.cos_theta < CT[i + 1]) in_a_gap = true;
            BOOST_CHECK(in_a_gap);
        }
    }
}

BOOST_AUTO_TEST_CASE(energygap_never_bridges_a_k_edge) {
    // Nodes 80 and 96 straddle a K-edge at 88 keV: that gap must NOT yield an
    // EnergyGap probe (the interpolant is segmented at the edge).
    const std::vector<double> E = {35, 60, 80, 96, 344, 662};
    const std::vector<double> CT = {0.02, 0.5, 1.0};
    const std::vector<double> edges = {88.0};
    std::shared_ptr<DetectorResponse> resp =
        make_response(E, CT, edges, {}, {}, {});

    GenerationOptions opts;
    opts.node_fep_precision = 0.003;
    const std::vector<ResponseGenerator::StructuredProbe> plan =
        ResponseGenerator::plan_structured_probes(
            resp->descriptor, opts, *resp,
            probe_family_bit(ProbeFamily::EnergyGap), 20, 20000);

    for (const ResponseGenerator::StructuredProbe& p : plan) {
        BOOST_REQUIRE(p.family == ProbeFamily::EnergyGap);
        // No midpoint may sit inside the edge-bridging (80, 96) gap.
        BOOST_CHECK(!(p.energy_keV > 80.0 && p.energy_keV < 96.0));
    }
    // 5 adjacent-node pairs, minus the edge-bridging (80,96) gap -> 4 probes.
    BOOST_CHECK_EQUAL(plan.size(), 4u);
}

BOOST_AUTO_TEST_CASE(family_mask_selects_only_requested) {
    const std::vector<double> E = {35, 122, 662, 3000};
    const std::vector<double> CT = {0.02, 0.5, 1.0};
    std::shared_ptr<DetectorResponse> resp = make_response(E, CT, {}, {}, {}, {});

    GenerationOptions opts;
    const std::vector<ResponseGenerator::StructuredProbe> plan =
        ResponseGenerator::plan_structured_probes(
            resp->descriptor, opts, *resp,
            probe_family_bit(ProbeFamily::CtGap), 4, 20000);
    BOOST_CHECK_GT(plan.size(), 0u);
    for (const ResponseGenerator::StructuredProbe& p : plan)
        BOOST_CHECK(p.family == ProbeFamily::CtGap);
}
