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
