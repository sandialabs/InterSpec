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

/// @file EnvironmentGeometry.h
/// @brief The environment tier: the air between source and detector, plus the
///        material slabs (ground, floor, walls, ceiling) that surround them.
///
/// CeeLo has two geometry systems today: `SourceGeometry` (the radioactive
/// volume and its shields) and `Geometry` (the crystal and its attenuators).
/// Between them sits open space that is not modeled at all. `EnvironmentGeometry`
/// fills that gap: a homogeneous air region bounded by an optional extent, plus
/// N half-space material slabs.
///
/// One representation covers every scenario of interest:
///   - 0 slabs              open air (a source and detector in free space)
///   - 1 slab               a field measurement over ground/soil
///   - 6 slabs              a room (floor, ceiling, four walls)
/// An infinite ground half-space is simply the degenerate one-slab room, so the
/// ground-plane and full-room cases share all of the tracing and scoring code.
///
/// IMPORTANT — this is an OFFLINE truth/bake tier, not a per-event transport
/// tier. It is used by the next-event fluence estimator to generate kernels; it
/// is deliberately NOT wired into `EfficiencyCalculator::compute()`. Routing a
/// full analog environment through the per-event efficiency loop would add cost
/// and an FEP-coupling hazard to the main path for a peak effect already
/// covered by the analytic `AirAttenuation` modes.
///
/// Frame: the detector frame used throughout CeeLo — crystal front face at
/// z = 0, crystal extending toward +z, source at negative z (see
/// `io/ResponseKernel.h::source_position`).

#include "materials/Material.h"

#include <Eigen/Core>
#include <array>
#include <vector>
#include <cstddef>

namespace ceelo {

/// A semi-infinite (or finite-thickness) material slab bounding the air region.
///
/// The slab occupies the region behind its plane: points p with
/// `s = (p - point).dot(normal)` in `[-thickness_cm, 0)`. The `normal` points
/// INTO the air, so positive `s` is the air side.
struct EnvSlab {
    Eigen::Vector3d point;    ///< A point on the inner (air-facing) plane.
    Eigen::Vector3d normal;   ///< Unit normal, pointing into the air region.
    double thickness_cm;      ///< Material depth measured along -normal.
    Material material;        ///< Slab material (held by value; owns lifetime).
};

/// Optional bound on the air region: a cylinder about an axis. This is the
/// "transverse extent", a real, independent kernel axis (measured ~20% effect
/// on the collected scatter at fixed distance) — the
/// surrounding air halo that scatters photons back toward the detector. An
/// unbounded (open-field) air region is the default.
struct EnvBounds {
    bool enabled = false;
    Eigen::Vector3d axis_point{0.0, 0.0, 0.0};
    Eigen::Vector3d axis_dir{0.0, 0.0, 1.0};  ///< Unit.
    double radius_cm = 0.0;                    ///< Transverse extent.
    double t_min_cm = 0.0;                     ///< Along axis from axis_point.
    double t_max_cm = 0.0;
};

/// One homogeneous stretch of a traced ray.
struct EnvSegment {
    double t_start_cm;
    double t_end_cm;
    const Material* material;  ///< nullptr => void (the ray has left the problem).
    int slab_index;            ///< -1 for air, else index into slabs().
};

/// What region a point lies in.
struct EnvPointClass {
    const Material* material;  ///< nullptr => void.
    int slab_index;            ///< -1 for air, else the slab it is inside.
};

/// Reusable storage for ray-boundary crossings.  Transport owns one of these
/// per worker so tracing a history and its detector legs does not allocate a
/// temporary vector at every collision.
struct EnvTraceScratch {
    std::vector<double> crossings;
};

class EnvironmentGeometry {
public:
    /// Open air with no slabs and no extent limit (free space).
    explicit EnvironmentGeometry(const Material& air);

    /// Open air bounded by a cylinder of radius `transverse_extent_cm` about
    /// the source->detector axis, running from behind the source to behind the
    /// detector. `source` and `detector` are in the detector frame.
    static EnvironmentGeometry open_air(const Material& air,
                                        const Eigen::Vector3d& source,
                                        const Eigen::Vector3d& detector,
                                        double transverse_extent_cm,
                                        double longitudinal_margin_cm = 0.0);

    /// A field measurement: source and detector at heights `h_src_cm` and
    /// `h_det_cm` above a semi-infinite ground half-space.
    ///
    /// The ground plane is inferred from the two positions and their heights:
    /// it is the unique plane at perpendicular distance h_det from the detector
    /// and h_src from the source (the two heights may differ — a source on the
    /// ground viewed by a detector at waist height is the common case). The
    /// source is assumed to be at `source` in the detector frame with the
    /// detector at the origin.
    static EnvironmentGeometry field(const Material& air,
                                     const Material& ground,
                                     const Eigen::Vector3d& source,
                                     double h_src_cm,
                                     double h_det_cm,
                                     double ground_thickness_cm = 200.0);

    /// A rectangular room: floor, ceiling and four walls, all of one material.
    ///
    /// The room is built in a "survey" frame whose axes are: `u` along the
    /// source->detector line (horizontal), `v` horizontal and perpendicular to
    /// it, and `w` vertical. The detector sits at the origin at height
    /// `h_det_cm`; the source is at `source`, at height `h_src_cm`. The room
    /// extends `back_cm` behind the detector and `front_cm` beyond the source
    /// along u, `half_width_cm` either side along v, and up to `height_cm`.
    ///
    /// A ground-only field measurement is the same construction with the
    /// ceiling and walls removed, which is why both share every line of the
    /// tracing and scoring code.
    static EnvironmentGeometry room(const Material& air,
                                    const Material& surfaces,
                                    const Eigen::Vector3d& source,
                                    double h_src_cm,
                                    double h_det_cm,
                                    double height_cm,
                                    double half_width_cm,
                                    double back_cm,
                                    double front_cm,
                                    double wall_thickness_cm = 30.0);

    /// Add a slab. `normal` need not be normalized.
    void add_slab(const Eigen::Vector3d& point,
                  const Eigen::Vector3d& normal,
                  double thickness_cm,
                  const Material& material);

    /// Bound the air region (see EnvBounds).
    void set_bounds(const EnvBounds& bounds) { bounds_ = bounds; }
    const EnvBounds& bounds() const { return bounds_; }

    /// Classify a point: which region (air, a slab, or void) contains it.
    ///
    /// Corner rule: when a point lies behind more than one slab plane (a corner
    /// wedge, e.g. where a wall meets the floor), it is assigned to the slab it
    /// has penetrated LEAST — the surface it is nearest to. Such wedges are a
    /// measure-~0 population and the adjacent materials are usually identical,
    /// so this is a documented convention rather than a modeled effect.
    EnvPointClass classify(const Eigen::Vector3d& p) const;

    /// Trace a ray from `pos` along unit `dir`, appending homogeneous segments
    /// to `out` (cleared first) until the ray leaves the problem (enters void)
    /// or `max_distance_cm` is reached. Void is not appended.
    void trace(const Eigen::Vector3d& pos,
               const Eigen::Vector3d& dir,
               std::vector<EnvSegment>& out,
               double max_distance_cm = 1.0e9) const;

    /// Allocation-free overload using caller-owned crossing storage.
    void trace(const Eigen::Vector3d& pos,
               const Eigen::Vector3d& dir,
               std::vector<EnvSegment>& out,
               EnvTraceScratch& scratch,
               double max_distance_cm = 1.0e9) const;

    /// Optical depth tau = sum(mu_total * length) along the straight chord from
    /// `from` to `to`. This is the vertex->detector attenuation leg of the
    /// next-event estimator; it must cross every region uniformly, which is what
    /// makes scoring correct for a vertex buried inside a wall or the ground.
    double optical_depth(const Eigen::Vector3d& from,
                         const Eigen::Vector3d& to,
                         double energy_MeV) const;

    /// Allocation-free optical-depth overload.
    double optical_depth(const Eigen::Vector3d& from,
                         const Eigen::Vector3d& to,
                         double energy_MeV,
                         EnvTraceScratch& scratch) const;

    /// exp(-optical_depth(...)) — the attenuated arrival probability.
    double transmission(const Eigen::Vector3d& from,
                        const Eigen::Vector3d& to,
                        double energy_MeV) const;

    /// Allocation-free transmission overload.
    double transmission(const Eigen::Vector3d& from,
                        const Eigen::Vector3d& to,
                        double energy_MeV,
                        EnvTraceScratch& scratch) const;

    /// Transmissions for two energies using one exact geometry trace.  The two
    /// optical-depth sums retain the same segment order as independent calls;
    /// only boundary discovery and classification are shared.
    std::array<double, 2> transmission_pair(
        const Eigen::Vector3d& from,
        const Eigen::Vector3d& to,
        double first_energy_MeV,
        double second_energy_MeV,
        EnvTraceScratch& scratch) const;

    const Material& air() const { return air_; }
    const std::vector<EnvSlab>& slabs() const { return slabs_; }
    std::size_t num_slabs() const { return slabs_.size(); }

private:
    /// Collect and sort candidate boundary crossing distances along a ray.
    void gather_crossings(const Eigen::Vector3d& pos,
                          const Eigen::Vector3d& dir,
                          double max_distance_cm,
                          std::vector<double>& ts) const;

    Material air_;
    std::vector<EnvSlab> slabs_;
    EnvBounds bounds_;
};

} // namespace ceelo
