#pragma once
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

/// @file SourceDimensionEstimator.h
/// @brief Estimates minimum source dimensions to capture a target fraction
///        of the FEP signal from an extended (soil contamination) source.
///
/// This is a fast analytical/numerical pre-computation (no MC). Given a detector
/// geometry, source medium, photon energy, and depth distribution, it estimates
/// the minimum source depth and lateral extent needed so that a finite source
/// volume captures at least X% of the FEP signal from an infinitely large source.
///
/// Typical usage:
/// @code
///   DimensionEstimatorConfig cfg;
///   cfg.energy_keV = 662.0;
///   cfg.distance_cm = 100.0;       // 1 m above ground
///   cfg.medium = &soil;
///   cfg.medium_thickness_cm = 200.0;
///   cfg.coverage_fraction = 0.95;
///   cfg.depth_distribution = DepthDistribution::Exponential;
///   cfg.relaxation_length_cm = 3.0;
///
///   auto result = estimate_source_dimensions(cfg, calc.geometry());
///
///   // Use result directly:
///   calc.set_cylindrical_source(result.center,
///                               result.recommended_radius_cm,
///                               result.recommended_depth_cm / 2.0);
/// @endcode

#include "efficiency/EfficiencyCalculator.h"  // for DepthDistribution
#include "geometry/Geometry.h"
#include "materials/Material.h"

#include <Eigen/Core>

namespace ceelo {

/// Configuration for source dimension estimation.
struct DimensionEstimatorConfig {
    double energy_keV = 662.0;            ///< Photon energy
    double distance_cm = 10.0;            ///< Distance from detector front face to source surface
    const Material* medium = nullptr;     ///< Source medium material (e.g., soil)
    double medium_thickness_cm = 100.0;   ///< Total available medium thickness (cm)
    double coverage_fraction = 0.95;      ///< Target coverage (e.g., 0.95 for 95%)
    DepthDistribution depth_distribution = DepthDistribution::Uniform;
    double relaxation_length_cm = 1.0;    ///< Relaxation length (cm), only for Exponential
};

/// Result of source dimension estimation.
/// Provides ready-to-use source configuration parameters.
struct DimensionEstimatorResult {
    /// Source center position in the detector frame (cm).
    /// Placed on-axis at z = -(distance + depth/2).
    Eigen::Vector3d center{0, 0, 0};

    /// Recommended source depth (total, not half). Use depth/2 as half_length or half_z.
    double recommended_depth_cm = 0.0;

    /// Recommended source radius (for cylindrical source).
    double recommended_radius_cm = 0.0;

    /// Recommended half-widths (for rectangular source). Set equal (square) by default.
    double recommended_half_x_cm = 0.0;
    double recommended_half_y_cm = 0.0;

    // --- Diagnostics ---

    /// Actual depth coverage achieved (may exceed target due to medium thickness limit).
    double depth_coverage_fraction = 0.0;

    /// Actual lateral coverage achieved.
    double lateral_coverage_fraction = 0.0;

    /// Total attenuation coefficient of the medium at the given energy (1/cm).
    double mu_medium_per_cm = 0.0;

    /// Effective signal depth scale (cm):
    /// 1/(1/lambda + mu) for exponential, 1/mu for uniform.
    double effective_depth_scale_cm = 0.0;
};

/// Estimate minimum source dimensions (depth + lateral extent) to capture
/// the target fraction of FEP signal from an infinitely large source.
///
/// The depth estimate uses the analytical integral of the depth-weighted
/// attenuation function. The lateral estimate uses numerical integration
/// of the solid-angle-weighted contribution from annular rings.
///
/// @param config   Estimator configuration (energy, distance, medium, coverage)
/// @param detector The configured detector geometry (for face area / bounding box)
/// @return         Recommended source dimensions and diagnostic info
DimensionEstimatorResult estimate_source_dimensions(
    const DimensionEstimatorConfig& config,
    const Geometry& detector);

} // namespace ceelo
