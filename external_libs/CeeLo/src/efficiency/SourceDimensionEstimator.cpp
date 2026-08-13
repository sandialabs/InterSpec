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

#include "efficiency/SourceDimensionEstimator.h"
#include "materials/Material.h"

#include <cmath>
#include <cassert>
#include <algorithm>

namespace ceelo {

namespace {

/// Estimate the required source depth to capture the target coverage fraction.
///
/// For exponential distribution, the signal from depth d is proportional to:
///   exp(-d/lambda) * exp(-mu*d) = exp(-d/L_eff)
/// where L_eff = 1/(1/lambda + mu).
///
/// For uniform distribution, the signal from depth d is:
///   exp(-mu*d)
/// which is the same formula with L_eff = 1/mu.
///
/// The coverage fraction from [0, d] out of [0, D_max] is:
///   F(d) = (1 - exp(-d/L_eff)) / (1 - exp(-D_max/L_eff))
///
/// Invert to find d for a target F:
///   d = -L_eff * ln(1 - F * (1 - exp(-D_max/L_eff)))
double estimate_depth(double mu_medium, double D_max, double coverage,
                      DepthDistribution dist, double relaxation_length)
{
    assert(coverage > 0.0 && coverage < 1.0);
    assert(D_max > 0.0);

    double L_eff;
    if (dist == DepthDistribution::Exponential) {
        assert(relaxation_length > 0.0);
        L_eff = 1.0 / (1.0 / relaxation_length + mu_medium);
    } else {
        if (mu_medium < 1e-10) {
            // Very low attenuation — need the full depth for coverage
            return D_max * coverage;
        }
        L_eff = 1.0 / mu_medium;
    }

    double ratio = D_max / L_eff;
    double exp_neg_ratio = (ratio > 700.0) ? 0.0 : std::exp(-ratio);

    // F(D_max) = 1.0 by definition (normalized). The actual coverage of the
    // full medium thickness is 1.0, so we just need to find the depth d where
    // F(d) = coverage.
    double arg = 1.0 - coverage * (1.0 - exp_neg_ratio);
    if (arg <= 0.0) {
        // Coverage requires the full depth (or more)
        return D_max;
    }

    double d = -L_eff * std::log(arg);
    return std::min(d, D_max);
}

/// Estimate the required lateral radius to capture the target coverage fraction.
///
/// The signal contribution from an annular ring at lateral radius r is:
///   dC/dr = r * A_det * distance / (r^2 + distance^2)^(3/2) * exp(-mu_air * R)
/// where R = sqrt(r^2 + distance^2) is the total path length.
///
/// We integrate numerically from 0 to r_max and find the radius where the
/// cumulative integral reaches the target fraction of the total.
double estimate_lateral_radius(double distance, double mu_air, double coverage)
{
    assert(distance > 0.0);
    assert(coverage > 0.0 && coverage < 1.0);

    // Integrate to a large radius: the integrand falls off as ~1/r^2 for large r,
    // so the integral converges. Use r_max large enough that the remaining tail
    // contribution is negligible.
    double r_max = std::max(50.0 * distance, 500.0);

    // Use 4000 steps for good accuracy
    const int N = 4000;
    double dr = r_max / N;

    // Trapezoidal integration
    // Integrand: r * distance / (r^2 + distance^2)^(3/2) * exp(-mu_air * sqrt(r^2 + distance^2))
    // At r=0: integrand = 0 (due to the r factor)
    double cumulative = 0.0;
    double prev_f = 0.0;  // f(0) = 0

    double target_radius = r_max;
    bool found = false;

    // Pre-compute total integral first to get the target value
    double total_integral = 0.0;
    for (int i = 1; i <= N; ++i) {
        double r = i * dr;
        double R2 = r * r + distance * distance;
        double R = std::sqrt(R2);
        double f = r * distance / (R2 * R) * std::exp(-mu_air * R);
        total_integral += 0.5 * (prev_f + f) * dr;
        prev_f = f;
    }

    double target = coverage * total_integral;

    // Now integrate again to find where cumulative reaches the target
    prev_f = 0.0;
    cumulative = 0.0;
    for (int i = 1; i <= N; ++i) {
        double r = i * dr;
        double R2 = r * r + distance * distance;
        double R = std::sqrt(R2);
        double f = r * distance / (R2 * R) * std::exp(-mu_air * R);
        cumulative += 0.5 * (prev_f + f) * dr;
        prev_f = f;

        if (!found && cumulative >= target) {
            // Linear interpolation for refinement
            double prev_cumulative = cumulative - 0.5 * (prev_f + f) * dr;
            double frac = (target - prev_cumulative) / (cumulative - prev_cumulative);
            target_radius = (i - 1 + frac) * dr;
            found = true;
            break;
        }
    }

    return target_radius;
}

} // anonymous namespace

DimensionEstimatorResult estimate_source_dimensions(
    const DimensionEstimatorConfig& config,
    const Geometry& /*detector*/)
{
    assert(config.medium != nullptr);
    assert(config.energy_keV > 0.0);
    assert(config.distance_cm > 0.0);
    assert(config.coverage_fraction > 0.0 && config.coverage_fraction < 1.0);
    assert(config.medium_thickness_cm > 0.0);

    DimensionEstimatorResult result;

    double energy_MeV = config.energy_keV * 1e-3;
    double mu_medium = config.medium->mu_total(energy_MeV);
    result.mu_medium_per_cm = mu_medium;

    // Effective depth scale
    if (config.depth_distribution == DepthDistribution::Exponential) {
        assert(config.relaxation_length_cm > 0.0);
        result.effective_depth_scale_cm = 1.0 / (1.0 / config.relaxation_length_cm + mu_medium);
    } else {
        result.effective_depth_scale_cm = (mu_medium > 1e-10) ? (1.0 / mu_medium) : config.medium_thickness_cm;
    }

    // --- Depth estimation ---
    double depth = estimate_depth(mu_medium, config.medium_thickness_cm,
                                  config.coverage_fraction,
                                  config.depth_distribution,
                                  config.relaxation_length_cm);
    result.recommended_depth_cm = depth;

    // Compute actual depth coverage achieved
    double L_eff = result.effective_depth_scale_cm;
    double ratio_d = depth / L_eff;
    double ratio_D = config.medium_thickness_cm / L_eff;
    if (ratio_D > 700.0) {
        result.depth_coverage_fraction = (ratio_d > 700.0) ? 1.0 : (1.0 - std::exp(-ratio_d));
    } else {
        result.depth_coverage_fraction =
            (1.0 - std::exp(-ratio_d)) / (1.0 - std::exp(-ratio_D));
    }

    // --- Lateral extent estimation ---
    // The lateral integrand is r * Omega(r) * T_air(r), where
    // Omega(r) ~ A_det * cos(alpha) / R^2. Since A_det is constant and cancels
    // in the coverage ratio, we only need the geometric factor:
    //   r * distance / (r^2 + distance^2)^(3/2) * exp(-mu_air * R)

    // Air attenuation coefficient
    Material air = make_Air();
    double mu_air = air.mu_total(energy_MeV);

    double lateral_radius = estimate_lateral_radius(
        config.distance_cm, mu_air, config.coverage_fraction);

    result.recommended_radius_cm = lateral_radius;
    result.recommended_half_x_cm = lateral_radius;
    result.recommended_half_y_cm = lateral_radius;

    // Compute actual lateral coverage
    // Re-integrate to the recommended radius and compare with total
    // (this is already done implicitly by the algorithm, so set directly)
    result.lateral_coverage_fraction = config.coverage_fraction;

    // --- Source center position ---
    // Detector front face at z=0, source surface at z = -distance,
    // source extends from -distance to -(distance + depth).
    // Center is at the midpoint: z = -(distance + depth/2).
    result.center = Eigen::Vector3d(0.0, 0.0, -(config.distance_cm + depth / 2.0));

    return result;
}

} // namespace ceelo
