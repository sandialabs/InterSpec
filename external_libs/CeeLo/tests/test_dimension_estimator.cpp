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

#define BOOST_TEST_MODULE DimensionEstimatorTests
#include <boost/test/unit_test.hpp>

#include "efficiency/SourceDimensionEstimator.h"
#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include <cmath>

using namespace ceelo;

namespace {

/// Helper: set up a standard 3"x3" NaI detector.
Geometry make_nai_3x3() {
    Material nai = make_NaI();
    Geometry geo;
    geo.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    return geo;
}

} // anonymous namespace

// ============================================================
//  Depth Estimation Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(DepthEstimation)

BOOST_AUTO_TEST_CASE(uniform_depth_sanity_662keV) {
    // 662 keV in soil: mu_total ~ 0.11 cm^-1 (rough estimate)
    // 95% coverage depth: -ln(0.05)/mu ~ 27 cm
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Uniform;

    auto result = estimate_source_dimensions(config, geo);

    // Depth should be in the ballpark of ~20-35 cm for 95% coverage
    BOOST_CHECK_GT(result.recommended_depth_cm, 10.0);
    BOOST_CHECK_LT(result.recommended_depth_cm, 60.0);

    // mu_medium should be positive and reasonable
    BOOST_CHECK_GT(result.mu_medium_per_cm, 0.05);
    BOOST_CHECK_LT(result.mu_medium_per_cm, 0.30);

    // Coverage should be at or above target
    BOOST_CHECK_GE(result.depth_coverage_fraction, 0.949);
}

BOOST_AUTO_TEST_CASE(exponential_depth_sanity_662keV) {
    // Exponential with lambda=1 cm: L_eff = 1/(1+mu) ~ 0.9 cm
    // 95% coverage: -L_eff * ln(0.05) ~ 2.7 cm
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Exponential;
    config.relaxation_length_cm = 1.0;

    auto result = estimate_source_dimensions(config, geo);

    // Depth should be much smaller than uniform case
    BOOST_CHECK_GT(result.recommended_depth_cm, 1.0);
    BOOST_CHECK_LT(result.recommended_depth_cm, 10.0);

    // Effective depth scale ~ 1/(1+mu)
    BOOST_CHECK_GT(result.effective_depth_scale_cm, 0.5);
    BOOST_CHECK_LT(result.effective_depth_scale_cm, 1.5);
}

BOOST_AUTO_TEST_CASE(monotonicity_coverage_vs_depth) {
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.depth_distribution = DepthDistribution::Uniform;

    config.coverage_fraction = 0.90;
    auto res90 = estimate_source_dimensions(config, geo);

    config.coverage_fraction = 0.95;
    auto res95 = estimate_source_dimensions(config, geo);

    config.coverage_fraction = 0.99;
    auto res99 = estimate_source_dimensions(config, geo);

    // Higher coverage should require greater depth
    BOOST_CHECK_LT(res90.recommended_depth_cm, res95.recommended_depth_cm);
    BOOST_CHECK_LT(res95.recommended_depth_cm, res99.recommended_depth_cm);
}

BOOST_AUTO_TEST_CASE(exponential_depth_leq_uniform_depth) {
    // Exponential concentrates activity near surface, so required depth
    // for same coverage should be smaller.
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.coverage_fraction = 0.95;

    config.depth_distribution = DepthDistribution::Uniform;
    auto res_uni = estimate_source_dimensions(config, geo);

    config.depth_distribution = DepthDistribution::Exponential;
    config.relaxation_length_cm = 2.0;
    auto res_exp = estimate_source_dimensions(config, geo);

    BOOST_CHECK_LT(res_exp.recommended_depth_cm, res_uni.recommended_depth_cm);
}

BOOST_AUTO_TEST_CASE(lower_energy_smaller_depth) {
    // Lower energy photons are more attenuated, so less depth is needed
    // to capture the same fraction of signal.
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Uniform;

    config.energy_keV = 100.0;
    auto res_low = estimate_source_dimensions(config, geo);

    config.energy_keV = 662.0;
    auto res_mid = estimate_source_dimensions(config, geo);

    config.energy_keV = 2614.0;
    auto res_high = estimate_source_dimensions(config, geo);

    // Lower energy → higher mu → less depth needed
    BOOST_CHECK_LT(res_low.recommended_depth_cm, res_mid.recommended_depth_cm);
    BOOST_CHECK_LT(res_mid.recommended_depth_cm, res_high.recommended_depth_cm);
}

BOOST_AUTO_TEST_CASE(thin_medium_caps_depth) {
    // If medium is thin (e.g., 1 mm surface contamination), depth is capped
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 0.1;  // 1 mm
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Uniform;

    auto result = estimate_source_dimensions(config, geo);

    // Depth should be capped at medium thickness
    BOOST_CHECK_LE(result.recommended_depth_cm, 0.1 + 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Lateral Extent Estimation Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(LateralEstimation)

BOOST_AUTO_TEST_CASE(lateral_radius_sanity) {
    // 3"x3" NaI at 10 cm distance: the lateral radius for 95% coverage
    // should be on the order of 30-100 cm (solid angle dominated).
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Uniform;

    auto result = estimate_source_dimensions(config, geo);

    BOOST_CHECK_GT(result.recommended_radius_cm, 10.0);
    BOOST_CHECK_LT(result.recommended_radius_cm, 300.0);

    // half_x and half_y should equal radius (conservative circular)
    BOOST_CHECK_CLOSE(result.recommended_half_x_cm,
                       result.recommended_radius_cm, 1e-6);
    BOOST_CHECK_CLOSE(result.recommended_half_y_cm,
                       result.recommended_radius_cm, 1e-6);
}

BOOST_AUTO_TEST_CASE(larger_distance_requires_larger_radius) {
    // Further from detector should require larger lateral extent
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Uniform;

    config.distance_cm = 10.0;
    auto res_near = estimate_source_dimensions(config, geo);

    config.distance_cm = 100.0;
    auto res_far = estimate_source_dimensions(config, geo);

    BOOST_CHECK_GT(res_far.recommended_radius_cm,
                    res_near.recommended_radius_cm);
}

BOOST_AUTO_TEST_CASE(monotonicity_coverage_vs_radius) {
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.depth_distribution = DepthDistribution::Uniform;

    config.coverage_fraction = 0.90;
    auto res90 = estimate_source_dimensions(config, geo);

    config.coverage_fraction = 0.95;
    auto res95 = estimate_source_dimensions(config, geo);

    config.coverage_fraction = 0.99;
    auto res99 = estimate_source_dimensions(config, geo);

    BOOST_CHECK_LT(res90.recommended_radius_cm, res95.recommended_radius_cm);
    BOOST_CHECK_LT(res95.recommended_radius_cm, res99.recommended_radius_cm);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Source Center Position Tests
// ============================================================

BOOST_AUTO_TEST_SUITE(CenterPosition)

BOOST_AUTO_TEST_CASE(center_position_correct) {
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 200.0;
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Exponential;
    config.relaxation_length_cm = 2.0;

    auto result = estimate_source_dimensions(config, geo);

    // Center should be on-axis
    BOOST_CHECK_CLOSE(result.center.x(), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(result.center.y(), 0.0, 1e-10);

    // Center z = -(distance + depth/2)
    double expected_z = -(config.distance_cm + result.recommended_depth_cm / 2.0);
    BOOST_CHECK_CLOSE(result.center.z(), expected_z, 1e-6);

    // Verify the source surface is at -distance (nearest to detector)
    double surface_z = result.center.z() + result.recommended_depth_cm / 2.0;
    BOOST_CHECK_CLOSE(surface_z, -config.distance_cm, 1e-6);
}

BOOST_AUTO_TEST_CASE(box_detector_works) {
    // Verify the estimator works with a box detector too
    Material czt = make_CZT();
    Geometry geo;
    geo.set_detector(DetectorShape::Box, &czt, {0.5, 0.5, 0.5});

    Material soil = make_Soil();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 5.0;
    config.medium = &soil;
    config.medium_thickness_cm = 100.0;
    config.coverage_fraction = 0.95;
    config.depth_distribution = DepthDistribution::Uniform;

    auto result = estimate_source_dimensions(config, geo);

    BOOST_CHECK_GT(result.recommended_depth_cm, 0.0);
    BOOST_CHECK_GT(result.recommended_radius_cm, 0.0);
    BOOST_CHECK_GT(result.mu_medium_per_cm, 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Analytical Reference Value Tests
// ============================================================
//
// The depth and lateral formulas have exact closed-form solutions
// in limiting cases. These tests verify the estimator matches
// those analytical values to high precision.
//
// Depth formula (for D_max >> L_eff, i.e. infinite-medium limit):
//   d = -L_eff * ln(1 - coverage)
//   where L_eff = 1/mu (uniform) or 1/(1/lambda + mu) (exponential)
//
// Lateral formula (for mu_air = 0, i.e. no air attenuation):
//   F(R) = 1 - h / sqrt(R^2 + h^2)
//   => R = h * sqrt(1/(1-f)^2 - 1)
//   This is the exact integral of r*h/(r^2+h^2)^(3/2) dr from 0 to R,
//   normalized by the integral from 0 to infinity (which equals 1).

BOOST_AUTO_TEST_SUITE(AnalyticalReferenceValues)

BOOST_AUTO_TEST_CASE(depth_uniform_exact_infinite_medium) {
    // For D_max >> L_eff = 1/mu, the depth formula simplifies to:
    //   d = -(1/mu) * ln(1 - coverage)
    // We use a large D_max (1000 cm) so exp(-D_max * mu) ≈ 0.
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    double energy_MeV = 0.662;
    double mu = soil.mu_total(energy_MeV);
    double L_eff = 1.0 / mu;

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 1000.0;  // >> L_eff, so infinite-medium limit
    config.depth_distribution = DepthDistribution::Uniform;

    for (double coverage : {0.90, 0.95, 0.99}) {
        config.coverage_fraction = coverage;
        auto result = estimate_source_dimensions(config, geo);

        double expected_depth = -L_eff * std::log(1.0 - coverage);

        // Should match to within 0.01% (numerical precision)
        BOOST_CHECK_CLOSE(result.recommended_depth_cm, expected_depth, 0.01);
    }
}

BOOST_AUTO_TEST_CASE(depth_exponential_exact_infinite_medium) {
    // For D_max >> L_eff = 1/(1/lambda + mu):
    //   d = -L_eff * ln(1 - coverage)
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    double energy_MeV = 0.662;
    double mu = soil.mu_total(energy_MeV);

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 1000.0;
    config.depth_distribution = DepthDistribution::Exponential;

    // Test multiple relaxation lengths
    for (double lambda : {0.5, 1.0, 3.0, 10.0}) {
        config.relaxation_length_cm = lambda;
        double L_eff = 1.0 / (1.0 / lambda + mu);

        for (double coverage : {0.90, 0.95, 0.99}) {
            config.coverage_fraction = coverage;
            auto result = estimate_source_dimensions(config, geo);

            double expected_depth = -L_eff * std::log(1.0 - coverage);

            BOOST_CHECK_CLOSE(result.recommended_depth_cm, expected_depth, 0.01);

            // Also verify L_eff is reported correctly
            BOOST_CHECK_CLOSE(result.effective_depth_scale_cm, L_eff, 0.01);
        }
    }
}

BOOST_AUTO_TEST_CASE(depth_specific_reference_values) {
    // Verify specific depth values that can be computed by hand.
    // For coverage f: d = L_eff * ln(1/(1-f))
    //   90%: d = L_eff * ln(10)     ≈ 2.302585 * L_eff
    //   95%: d = L_eff * ln(20)     ≈ 2.995732 * L_eff
    //   99%: d = L_eff * ln(100)    ≈ 4.605170 * L_eff
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    double energy_MeV = 0.662;
    double mu = soil.mu_total(energy_MeV);
    double L_eff = 1.0 / mu;

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 1000.0;
    config.depth_distribution = DepthDistribution::Uniform;

    // 90% coverage: d/L_eff = ln(10)
    config.coverage_fraction = 0.90;
    auto res90 = estimate_source_dimensions(config, geo);
    BOOST_CHECK_CLOSE(res90.recommended_depth_cm / L_eff, std::log(10.0), 0.01);

    // 95% coverage: d/L_eff = ln(20)
    config.coverage_fraction = 0.95;
    auto res95 = estimate_source_dimensions(config, geo);
    BOOST_CHECK_CLOSE(res95.recommended_depth_cm / L_eff, std::log(20.0), 0.01);

    // 99% coverage: d/L_eff = ln(100)
    config.coverage_fraction = 0.99;
    auto res99 = estimate_source_dimensions(config, geo);
    BOOST_CHECK_CLOSE(res99.recommended_depth_cm / L_eff, std::log(100.0), 0.01);
}

BOOST_AUTO_TEST_CASE(depth_finite_medium_exact) {
    // For finite D_max, the exact formula is:
    //   d = -L_eff * ln(1 - coverage * (1 - exp(-D_max/L_eff)))
    // Test with moderate D_max where the correction matters.
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    double energy_MeV = 0.662;
    double mu = soil.mu_total(energy_MeV);
    double L_eff = 1.0 / mu;
    double D_max = 3.0 * L_eff;  // Moderate thickness: exp(-3) ≈ 0.05

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = D_max;
    config.depth_distribution = DepthDistribution::Uniform;
    config.coverage_fraction = 0.95;

    auto result = estimate_source_dimensions(config, geo);

    double exp_neg_ratio = std::exp(-D_max / L_eff);
    double expected_depth = -L_eff * std::log(1.0 - 0.95 * (1.0 - exp_neg_ratio));

    BOOST_CHECK_CLOSE(result.recommended_depth_cm, expected_depth, 0.01);
}

BOOST_AUTO_TEST_CASE(lateral_radius_bounded_by_analytical_no_air) {
    // Without air attenuation, the exact lateral coverage function is:
    //   F(R) = 1 - h / sqrt(R^2 + h^2)
    // Solving for R at coverage f:
    //   R = h * sqrt(1/(1-f)^2 - 1)
    //
    // The estimator includes air attenuation (mu_air ≈ 8.5e-5 /cm at 662 keV),
    // which suppresses distant contributions. Therefore the estimator's R should
    // always be LESS than the air-free analytical value (air helps converge the
    // integral, so less lateral extent is needed for the same coverage).
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.medium = &soil;
    config.medium_thickness_cm = 1000.0;
    config.depth_distribution = DepthDistribution::Uniform;

    for (double h : {1.0, 5.0, 10.0, 50.0}) {
        config.distance_cm = h;
        for (double f : {0.90, 0.95}) {
            config.coverage_fraction = f;
            auto result = estimate_source_dimensions(config, geo);

            double exact_no_air = h * std::sqrt(1.0 / ((1.0 - f) * (1.0 - f)) - 1.0);

            // Estimator radius must be strictly less than air-free value
            BOOST_CHECK_LT(result.recommended_radius_cm, exact_no_air);

            // But should be at least 50% of it (air correction is modest at 662 keV)
            BOOST_CHECK_GT(result.recommended_radius_cm, exact_no_air * 0.5);
        }
    }
}

BOOST_AUTO_TEST_CASE(lateral_scaling_with_distance) {
    // For a fixed coverage fraction and negligible air, R should scale
    // linearly with h. Verify R(2h)/R(h) ≈ 2. Use small distances where
    // air attenuation is negligible (mu_air ≈ 8.5e-5 /cm at 662 keV,
    // so at R ≈ 20 cm, exp(-mu*R) ≈ 0.998).
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.energy_keV = 662.0;
    config.medium = &soil;
    config.medium_thickness_cm = 1000.0;
    config.depth_distribution = DepthDistribution::Uniform;
    config.coverage_fraction = 0.90;

    config.distance_cm = 1.0;
    auto res1 = estimate_source_dimensions(config, geo);

    config.distance_cm = 2.0;
    auto res2 = estimate_source_dimensions(config, geo);

    // R should scale roughly linearly with distance at small h.
    // The ~2% deviation from exact 2.0 comes from numerical discretization
    // (trapezoidal with 4000 steps) and residual air correction.
    double ratio = res2.recommended_radius_cm / res1.recommended_radius_cm;
    BOOST_CHECK_CLOSE(ratio, 2.0, 3.0);  // within 3% of exact linear scaling
}

BOOST_AUTO_TEST_CASE(lateral_air_correction_increases_with_energy_decrease) {
    // At lower energies, mu_air is larger, so the air correction is bigger,
    // meaning the lateral radius should be smaller. Compare 662 keV vs 100 keV.
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.distance_cm = 50.0;
    config.medium = &soil;
    config.medium_thickness_cm = 1000.0;
    config.depth_distribution = DepthDistribution::Uniform;
    config.coverage_fraction = 0.95;

    config.energy_keV = 662.0;
    auto res_high = estimate_source_dimensions(config, geo);

    config.energy_keV = 100.0;
    auto res_low = estimate_source_dimensions(config, geo);

    // Lower energy → more air attenuation → smaller lateral radius
    BOOST_CHECK_LT(res_low.recommended_radius_cm, res_high.recommended_radius_cm);

    // Both should still be less than the air-free analytical value
    double h = config.distance_cm;
    double exact_no_air = h * std::sqrt(1.0 / (0.05 * 0.05) - 1.0);
    BOOST_CHECK_LT(res_high.recommended_radius_cm, exact_no_air);
    BOOST_CHECK_LT(res_low.recommended_radius_cm, exact_no_air);
}

BOOST_AUTO_TEST_CASE(effective_depth_scale_values) {
    // Verify L_eff = 1/(1/lambda + mu) for exponential
    // and L_eff = 1/mu for uniform, using the actual mu from the material.
    Material soil = make_Soil();
    auto geo = make_nai_3x3();

    DimensionEstimatorConfig config;
    config.distance_cm = 10.0;
    config.medium = &soil;
    config.medium_thickness_cm = 1000.0;
    config.coverage_fraction = 0.95;

    // Test at multiple energies
    for (double E_keV : {100.0, 662.0, 1332.0, 2614.0}) {
        config.energy_keV = E_keV;
        double mu = soil.mu_total(E_keV * 1e-3);

        // Uniform: L_eff = 1/mu
        config.depth_distribution = DepthDistribution::Uniform;
        auto res_uni = estimate_source_dimensions(config, geo);
        BOOST_CHECK_CLOSE(res_uni.effective_depth_scale_cm, 1.0 / mu, 0.01);
        BOOST_CHECK_CLOSE(res_uni.mu_medium_per_cm, mu, 0.01);

        // Exponential: L_eff = 1/(1/lambda + mu)
        config.depth_distribution = DepthDistribution::Exponential;
        for (double lambda : {0.5, 2.0, 10.0}) {
            config.relaxation_length_cm = lambda;
            auto res_exp = estimate_source_dimensions(config, geo);
            double expected_L = 1.0 / (1.0 / lambda + mu);
            BOOST_CHECK_CLOSE(res_exp.effective_depth_scale_cm, expected_L, 0.01);

            // L_eff should always be less than both lambda and 1/mu
            BOOST_CHECK_LT(res_exp.effective_depth_scale_cm, lambda + 1e-10);
            BOOST_CHECK_LT(res_exp.effective_depth_scale_cm, 1.0 / mu + 1e-10);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
