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

/// @file Geometry.h
/// @brief Detector geometry definitions and ray-tracing interface.
///
/// The geometry consists of:
///   1. A detector crystal (cylinder or box), optionally with a bore hole
///   2. A dead layer (innermost non-scoring shell, same material as detector)
///   3. N attenuator layers (concentric shells around the detector)
///
/// All geometry is centered on the z-axis with the detector front face at z=0.
/// The detector extends from z=0 to z=L.

#include <Eigen/Core>
#include <vector>
#include <cstdint>
#include <optional>
#include <limits>
#include <string>

namespace ceelo {

class Material; // Forward declaration

/// Shape type for the detector
enum class DetectorShape : uint8_t {
    Cylinder,
    Box
};

/// Result of a ray intersection with a single geometric primitive.
struct RayHit {
    double t_enter;  ///< Parameter along ray where it enters the shape
    double t_exit;   ///< Parameter along ray where it exits the shape

    bool valid() const {
        return t_enter < t_exit
            && t_exit > 0.0;
    }

    double length() const {
        return t_exit - std::max(t_enter, 0.0);
    }
};

/// A segment of a ray through a specific material region.
struct PathSegment {
    double t_start;         ///< Start parameter along ray (from ray origin)
    double t_end;           ///< End parameter along ray
    const Material* material;  ///< Material in this segment (nullptr = vacuum)
    bool is_scoring;        ///< True if this is the active detector volume

    double length() const { return t_end - t_start; }
};

/// Configuration for a bore hole (coaxial HPGe).
struct BoreHoleConfig {
    double radius;    ///< Bore radius in cm
    double depth;     ///< Bore depth from back face in cm
    bool rounded_tip = false;  ///< Hemispherical closed end (round-tipped drill)
                               ///< instead of a flat bottom; the stated depth is
                               ///< kept, so the hemisphere's apex sits at the
                               ///< same z the flat bottom would.
};

/// Configuration for the dead layer.
struct DeadLayerConfig {
    double front;     ///< Front face dead layer thickness in cm
    double side;      ///< Side dead layer thickness in cm
    double back;      ///< Back dead layer thickness in cm
};

/// Configuration for an attenuator shell.
struct AttenuatorConfig {
    const Material* material;     ///< Attenuator material
    double front_thickness;  ///< Thickness on front face (cm)
    double side_thickness;   ///< Thickness on sides (cm)
    double z_start;          ///< Axial start of this attenuator
    double z_end;            ///< Axial end of this attenuator
    bool side_only = false;  ///< Tube geometry: side walls only, no front/back disks.
                             ///< Used for collimators that extend past the crystal.
};

/// The complete detector geometry: crystal + dead layer + attenuator shells.
class Geometry {
public:
    Geometry();

    /// Set the detector crystal shape and dimensions.
    /// For Cylinder: dimensions = {radius, length}
    /// For Box: dimensions = {half_width_x, half_width_y, length}
    ///
    /// Clears any bore hole, dead layer and bulletization radius: those are
    /// sized against the crystal, so redefining it invalidates them. Call this
    /// first, then declare the rest.
    void set_detector(DetectorShape shape, const Material* material,
                      const std::vector<double>& dimensions);

    /// Set bore hole for coaxial HPGe.
    /// Bore extends from the back face inward along the z-axis.
    /// With `rounded_tip` the closed end is a hemisphere of the bore radius
    /// (a round-tipped drill) rather than a flat bottom; the bore depth is
    /// measured to the apex either way.
    void set_bore_hole(double bore_radius, double bore_depth,
                       bool rounded_tip = false);

    /// Round ("bulletize") the outer front edge of a cylindrical crystal with
    /// a quarter-torus fillet of the given radius, as HPGe crystals usually
    /// are.  Pass 0 (the default) for a sharp 90-degree edge -- that path is
    /// bit-for-bit the same as before this feature existed.
    ///
    /// Composes with a bore hole and with a dead layer: the active volume gets
    /// the inward-offset fillet (see trace_cylinder_geometry).
    void set_bullet_radius(double bullet_radius);

    /// Set dead layer thicknesses.
    /// The dead layer uses the same material as the detector crystal.
    void set_dead_layer(double front, double side, double back = 0.0);

    /// Add an attenuator layer (call from innermost to outermost).
    void add_attenuator(const Material* material,
                        double front_thickness, double side_thickness,
                        double z_start, double z_end);

    /// Add a collimator (tube-shaped, side walls only, no front/back disks).
    /// The inner bore extends the full z range, creating an open tube.
    /// Used for W collimators that extend past the crystal front face.
    void add_collimator(const Material* material, double side_thickness,
                        double z_start, double z_end);

    /// Trace a ray through the complete geometry (attenuators + dead layer + crystal).
    /// Returns a list of path segments from the ray origin, ordered by distance.
    /// Each segment identifies the material and whether it's scoring (active crystal).
    ///
    /// @param origin  Ray origin (typically the source position)
    /// @param direction  Ray direction (must be normalized)
    /// @return  Ordered list of path segments through materials
    std::vector<PathSegment> trace_ray(const Eigen::Vector3d& origin,
                                       const Eigen::Vector3d& direction) const;

    /// Allocation-aware overload that writes into a caller-provided vector.
    /// Reuses vector capacity across repeated calls in tight transport loops.
    void trace_ray(const Eigen::Vector3d& origin,
                   const Eigen::Vector3d& direction,
                   std::vector<PathSegment>& segments) const;

    /// Quick test: does the ray hit the outermost boundary of the geometry?
    /// Useful for cone-sampling checks.
    bool ray_hits_outer_boundary(const Eigen::Vector3d& origin,
                                 const Eigen::Vector3d& direction) const;

    /// Get the bounding radius of the outermost shell (for cone sampling).
    /// For a cylindrical detector, this is the outer radius including all attenuators.
    /// For a box detector, this is the half-diagonal of the outermost shell.
    double outer_bounding_radius() const;

    /// Get the full axial extent [z_min, z_max] of the outermost shell.
    std::pair<double, double> outer_z_extent() const;

    // Accessors
    DetectorShape shape() const { return shape_; }
    const Material* detector_material() const { return detector_material_; }
    bool has_bore_hole() const { return bore_hole_.has_value(); }
    bool has_dead_layer() const { return dead_layer_.has_value(); }
    bool has_bullet_radius() const { return bullet_radius_ > 0.0; }
    double bullet_radius() const { return bullet_radius_; }
    size_t num_attenuators() const { return attenuators_.size(); }

    // Accessors for geometry export
    const std::optional<BoreHoleConfig>& bore_hole() const { return bore_hole_; }
    const std::optional<DeadLayerConfig>& dead_layer() const { return dead_layer_; }
    const std::vector<AttenuatorConfig>& attenuators() const { return attenuators_; }

    // Detector dimensions
    double detector_radius() const; ///< Only for cylindrical
    double detector_half_x() const; ///< Only for box
    double detector_half_y() const; ///< Only for box
    double detector_length() const;

private:
    DetectorShape shape_ = DetectorShape::Cylinder;
    const Material* detector_material_ = nullptr;

    // Cylinder: radius_, length_
    // Box: half_x_, half_y_, length_
    double radius_ = 0.0;
    double half_x_ = 0.0;
    double half_y_ = 0.0;
    double length_ = 0.0;

    /// Rounded front outer edge; 0 = sharp corner (the common case, and the
    /// only value the pre-bulletization code paths ever see).
    double bullet_radius_ = 0.0;

    std::optional<BoreHoleConfig> bore_hole_;
    std::optional<DeadLayerConfig> dead_layer_;
    std::vector<AttenuatorConfig> attenuators_;

    // Internal helpers for building the path segments
    void trace_cylinder_geometry(const Eigen::Vector3d& origin,
                                 const Eigen::Vector3d& direction,
                                 std::vector<PathSegment>& segments) const;
    void trace_box_geometry(const Eigen::Vector3d& origin,
                            const Eigen::Vector3d& direction,
                            std::vector<PathSegment>& segments) const;
};

} // namespace ceelo
