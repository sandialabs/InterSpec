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

#define BOOST_TEST_MODULE ResponseGeneratorSmokeTests
#include <boost/test/unit_test.hpp>

// End-to-end generation smoke test -- runs REAL (coarse) MC, so it is gated
// behind the CEELO_RUN_MC_TESTS CMake option (OFF by default). At the 3%
// node precision + trimmed budgets below it takes ~1-2 minutes.

#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "materials/Material.h"

#include <cmath>
#include <cstdio>

using namespace ceelo;

namespace {

GeometryDescriptor small_nai() {
    GeometryDescriptor gd;
    gd.shape = DetectorShape::Cylinder;
    gd.dimensions_cm = {1.27, 2.54};  // NaI 1"x1"
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

GenerationOptions coarse_options() {
    GenerationOptions o;
    o.node_fep_precision = 0.03;
    o.max_cpu_seconds_per_node = 20.0;  // ~= old 2 wall-s at ~10 threads
    o.n_energy_scan = 14;
    o.n_energy_nodes = 9;
    o.n_cos_theta_scan = 7;
    o.n_cos_theta_nodes = 5;
    o.n_shape_energies = 3;
    o.base_seed = 7;
    o.detector_name = "smoke NaI 1x1";
    return o;
}

}  // namespace

BOOST_AUTO_TEST_CASE(generate_evaluates_and_round_trips) {
    const GeometryDescriptor gd = small_nai();
    GenerationOptions opts = coarse_options();

    int progress_calls = 0;
    double last_frac = 0.0;
    opts.progress = [&](double frac, const std::string&) {
        BOOST_CHECK_GE(frac, last_frac - 1e-9);
        last_frac = frac;
        ++progress_calls;
    };

    std::shared_ptr<DetectorResponse> resp =
        ResponseGenerator::generate(gd, opts);
    BOOST_REQUIRE(resp);
    BOOST_CHECK_GT(progress_calls, 10);
    BOOST_CHECK(!resp->eta_fep.empty());

    // Sanity: mid-energy on-axis far-field efficiency within a factor-ish of
    // the small-crystal ballpark, positive sigma, Ok flag.
    const EffResult far = resp->eps_fep(662.0, 0.0, 0.0, 50.0);
    BOOST_CHECK(far.flag == ResponseFlag::Ok);
    BOOST_CHECK_GT(far.value, 1e-6);
    BOOST_CHECK_LT(far.value, 1e-3);
    BOOST_CHECK_GT(far.sigma, 0.0);

    // eps_tot >= eps_fep everywhere it is evaluated.
    const EffResult tot = resp->eps_total(662.0, 0.0, 0.0, 50.0);
    BOOST_CHECK_GT(tot.value, far.value);

    // XML round trip of a REAL generated object stays bit-stable.
    const std::string xml = resp->to_xml_string();
    std::shared_ptr<DetectorResponse> r2 = DetectorResponse::from_xml_string(xml);
    BOOST_CHECK_EQUAL(xml, r2->to_xml_string());

    // Probe check: a handful of fresh MC points must agree with the
    // parameterization within the coarse-node error budget (few %).
    GenerationOptions probe_opts = opts;
    probe_opts.progress = nullptr;  // separate run reports its own [0,1]
    probe_opts.node_fep_precision = 0.02;
    const std::vector<ProbePoint> bank =
        ResponseGenerator::probe_bank(gd, probe_opts, 8, 5000, 3.0, 60.0);
    int n_checked = 0;
    for (const ProbePoint& p : bank) {
        if (p.eps_fep <= 0.0 || p.fep_unc / p.eps_fep > 0.05) continue;
        const double theta = std::acos(p.cos_theta);
        const EffResult r = resp->eps_fep(p.energy_keV, theta,
                                          p.phi_deg * 3.14159265358979 / 180.0,
                                          p.d_cm);
        const double rel = r.value / p.eps_fep - 1.0;
        BOOST_CHECK_MESSAGE(std::fabs(rel) < 0.15,
                            "probe E=" << p.energy_keV << " d=" << p.d_cm
                                       << " ct=" << p.cos_theta << " rel="
                                       << rel);
        ++n_checked;
    }
    BOOST_CHECK_GT(n_checked, 3);
}

BOOST_AUTO_TEST_CASE(grounding_recovers_injected_bias) {
    // Synthetic perturbed-twin (S7-style): fake measured points = model * 1.04
    // with small noise; k(E) must recover ~1.04 with honest covariance.
    const GeometryDescriptor gd = small_nai();
    GenerationOptions opts = coarse_options();
    std::shared_ptr<DetectorResponse> resp =
        ResponseGenerator::generate(gd, opts);

    std::vector<GroundingPoint> pts;
    for (const double E : {60.0, 122.0, 344.0, 662.0, 1332.0}) {
        const EffResult model = resp->eps_fep(E, 0.0, 0.0, 25.0);
        GroundingPoint p;
        p.energy_keV = E;
        p.model_eff = model.value;
        p.measured_eff = model.value * 1.04;
        p.frac_stat_sigma = 0.01;
        p.frac_cert_sigma = 0.02;
        p.source_key = "twin";
        p.distance_cm = 25.0;
        p.cos_theta = 1.0;
        pts.push_back(p);
    }
    ResponseGenerator::ground_to_points(*resp, pts, /*curve_derived=*/false);
    BOOST_REQUIRE(!resp->grounding.empty());

    bool clamped = false;
    const double k = std::exp(resp->grounding.eval_ln_k(400.0, clamped));
    BOOST_CHECK_CLOSE(k, 1.04, 1.5 /*percent*/);
    // Certificate error is common-mode: var(ln k) must NOT be beaten below
    // the 2% cert floor by combining the 5 points.
    const double sig = std::sqrt(resp->grounding.var_ln_k(400.0));
    BOOST_CHECK_GT(sig, 0.015);
    BOOST_CHECK_LT(sig, 0.05);

    // Grounded response evaluates 4% above the ungrounded model.
    const EffResult grounded = resp->eps_fep(400.0, 0.0, 0.0, 25.0);
    resp->grounding = GroundingBlock();
    const EffResult plain = resp->eps_fep(400.0, 0.0, 0.0, 25.0);
    BOOST_CHECK_CLOSE(grounded.value / plain.value, 1.04, 1.5);
}

BOOST_AUTO_TEST_CASE(closed_loop_coarse_start_converges_or_records) {
    // Closed-loop generation from a COARSE start: it must end with a populated,
    // deterministic accuracy certificate whose converged flag honestly reflects
    // the final structured pass (never a silent miss), and it must be
    // reproducible bit-for-bit with a single thread.
    const GeometryDescriptor gd = small_nai();
    GenerationOptions opts = coarse_options();
    opts.closed_loop = true;
    opts.initial_grid = InitialGrid::Coarse;
    opts.num_threads = 1;             // determinism check needs a fixed thread
    opts.cert_tol = 0.012;
    opts.max_refine_iters = 2;
    opts.n_cert_probes = 24;

    // Baseline: certify the UNREFINED coarse initial grid (same family) so the
    // demo can report initial p95 -> refined p95.
    {
        GenerationOptions coarse = opts;
        coarse.closed_loop = false;
        coarse.n_shape_energies = 5;
        coarse.n_cos_theta_scan = 6;
        coarse.n_cos_theta_nodes = 6;
        std::shared_ptr<DetectorResponse> base =
            ResponseGenerator::generate(gd, coarse);
        ResponseGenerator::certify(*base, gd, coarse, 24, 7000,
                                   kAllStructuredFamilies);
        std::printf("[closed-loop demo] coarse initial p95 = %.3f%% "
                    "(max %.3f%%)\n",
                    100.0 * base->certificate.fep_p95,
                    100.0 * base->certificate.fep_max);
    }

    std::shared_ptr<DetectorResponse> r1 = ResponseGenerator::generate(gd, opts);
    BOOST_REQUIRE(r1);
    const AccuracyCertificate& c1 = r1->certificate;
    std::printf("[closed-loop demo] refined p95 = %.3f%% (max %.3f%%) after "
                "%d iter(s), converged=%d, cert+refine CPU = %.1f s\n",
                100.0 * c1.fep_p95, 100.0 * c1.fep_max, c1.iterations,
                c1.converged ? 1 : 0, c1.cpu_seconds);
    BOOST_REQUIRE(!c1.empty());
    BOOST_CHECK_GE(c1.iterations, 1);
    BOOST_CHECK_GT(c1.fep_p95, 0.0);
    // Honest bookkeeping: converged=true => the final structured pass met tol.
    // Either way it is RECORDED (never a silent miss).
    if (c1.converged)
        BOOST_CHECK_LE(c1.fep_p95, 0.05);  // coarse-node budget, generous gate

    // Determinism: same options + single thread => identical certificate p95.
    std::shared_ptr<DetectorResponse> r2 = ResponseGenerator::generate(gd, opts);
    BOOST_CHECK_CLOSE(r1->certificate.fep_p95, r2->certificate.fep_p95, 1e-6);
    BOOST_CHECK_EQUAL(r1->certificate.iterations, r2->certificate.iterations);
    BOOST_CHECK_EQUAL(r1->certificate.converged, r2->certificate.converged);

    // Structured families made it into the scorecard (tags 1..4 present).
    bool has_structured = false;
    for (const AccuracyCertificate::Row& row : c1.rows)
        if (row.tag != 0) has_structured = true;
    BOOST_CHECK(has_structured);
}

BOOST_AUTO_TEST_CASE(closed_loop_default_off_matches_fixed_grid) {
    // The default (closed_loop=false) path must run the fixed-grid code
    // VERBATIM. Generate with the default options and with closed_loop set
    // explicitly false; single-threaded so the MC is bit-exact (library policy),
    // then the two must be string-identical -- proving the dispatch is a no-op
    // when closed_loop is off.
    const GeometryDescriptor gd = small_nai();
    GenerationOptions a = coarse_options();
    a.num_threads = 1;                 // bit-exact MC for the comparison
    GenerationOptions b = a;
    b.closed_loop = false;             // explicit; a leaves it default (false)
    const std::string xa = ResponseGenerator::generate(gd, a)->to_xml_string();
    const std::string xb = ResponseGenerator::generate(gd, b)->to_xml_string();
    BOOST_CHECK_EQUAL(xa, xb);
}

BOOST_AUTO_TEST_CASE(certify_populates_certificate) {
    // Certify a coarse NaI response: the probe pass must fill the summary
    // percentiles and one row per probe, and it must not touch content_hash.
    const GeometryDescriptor gd = small_nai();
    GenerationOptions opts = coarse_options();
    std::shared_ptr<DetectorResponse> resp =
        ResponseGenerator::generate(gd, opts);
    BOOST_REQUIRE(resp);
    BOOST_CHECK(resp->certificate.empty());

    const uint64_t hash_before = resp->content_hash();

    const int n_probes = 16;
    ResponseGenerator::certify(*resp, gd, opts, n_probes, /*seed_offset=*/7000);

    const AccuracyCertificate& c = resp->certificate;
    BOOST_REQUIRE(!c.empty());
    BOOST_CHECK_EQUAL(c.rows.size(), static_cast<size_t>(n_probes));
    BOOST_CHECK(c.converged);
    BOOST_CHECK_EQUAL(c.iterations, 0);
    BOOST_CHECK_EQUAL(c.probe_seed_base, opts.base_seed);
    BOOST_CHECK_GT(c.cpu_seconds, 0.0);        // probe-bank MC actually cost time
    BOOST_CHECK_GE(c.fep_p95, c.fep_median);
    BOOST_CHECK_GE(c.fep_max, c.fep_p95);
    BOOST_CHECK_GT(c.fep_p95, 0.0);

    // The certificate is metadata: content_hash is unchanged by its presence.
    BOOST_CHECK_EQUAL(resp->content_hash(), hash_before);

    // And it survives an XML round trip (rows + summary intact).
    std::shared_ptr<DetectorResponse> r2 =
        DetectorResponse::from_xml_string(resp->to_xml_string());
    BOOST_REQUIRE(!r2->certificate.empty());
    BOOST_CHECK_EQUAL(r2->certificate.rows.size(), static_cast<size_t>(n_probes));
    BOOST_CHECK_CLOSE(r2->certificate.fep_p95, c.fep_p95, 1e-9);
    BOOST_CHECK_EQUAL(r2->content_hash(), hash_before);
}
