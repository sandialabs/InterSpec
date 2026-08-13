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

/// @file SourceGeometry.h
/// @brief Source self-attenuation and source-side shielding.
///
/// Models source material (fills the source volume) and concentric shielding
/// shells around the source. For point sources, shielding is spherical shells
/// (direction-independent path). For extended sources (cylindrical/rectangular),
/// shielding follows the source shape.

#include "materials/Material.h"

#include <Eigen/Core>
#include <vector>
#include <optional>
#include <random>
#include <cstddef>
#include <cstdint>

namespace ceelo {

/// A single shielding layer around the source.
/// Thicknesses (cm) are per local axis of the source shape:
///   - Point / Marinelli: uniform — tx == ty == tz; use scalar_thickness().
///   - Cylindrical: tx == ty = radial thickness, tz = end-cap thickness
///     (same on both end caps).
///   - Rectangular: tx, ty, tz = thickness on both +/- faces of each axis.
/// At least one component is > 0; individual components may be 0.
struct SourceShieldLayer {
    const Material* material;
    double tx;
    double ty;
    double tz;

    bool is_uniform() const { return tx == ty && ty == tz; }
    double scalar_thickness() const { return tx; }  ///< Only valid when is_uniform()
};

/// Source geometry: source material + concentric shielding shells.
class SourceGeometry {
public:
    SourceGeometry() = default;

    /// Set the material filling the source volume (e.g., soil, water).
    void set_source_material(const Material* mat) { source_material_ = mat; }
    const Material* source_material() const { return source_material_; }

    /// Add a uniform shielding layer (innermost first). Valid for all shapes.
    void add_shield(const Material* mat, double thickness);

    /// Add a cylindrical-source shielding layer with independent radial and
    /// end-cap thicknesses (cm). The end thickness applies to both end caps.
    /// One (but not both) may be zero, e.g. (t_radial, 0) = side wall only.
    /// Only valid for cylindrical sources.
    void add_shield(const Material* mat, double t_radial, double t_end);

    /// Add a rectangular-source shielding layer with independent x/y/z
    /// thicknesses (cm), applied on both +/- faces of each axis. One or two
    /// may be zero (but not all three). Only valid for rectangular sources.
    void add_shield(const Material* mat, double t_x, double t_y, double t_z);

    const std::vector<SourceShieldLayer>& shields() const { return shields_; }

    /// Configure for a point source.
    void configure_point(const Eigen::Vector3d& position);

    /// Configure for a cylindrical extended source.
    /// @param inner_radius  Inner (bore) radius for a hollow/annular cylinder
    ///   (tube, pipe, ring). 0 = solid cylinder. The active material occupies the
    ///   annulus [inner_radius, radius]; the central bore is an inactive,
    ///   non-attenuating void.
    void configure_cylindrical(const Eigen::Vector3d& center, double radius,
                               double half_length, const Eigen::Matrix3d& rotation,
                               double inner_radius = 0.0);

    /// Configure for a rectangular extended source.
    /// @param inner_half_dims  Inner void half-dimensions for a hollow box
    ///   shell (crate, container wall). All-zero = solid box. The active
    ///   material occupies the outer box minus the inner box (both centered,
    ///   same rotation); the inner box is an inactive, non-attenuating void.
    ///   Must satisfy 0 <= inner < outer componentwise, or be all zero.
    void configure_rectangular(const Eigen::Vector3d& center,
                               const Eigen::Vector3d& half_dims,
                               const Eigen::Matrix3d& rotation,
                               const Eigen::Vector3d& inner_half_dims
                                   = Eigen::Vector3d::Zero());

    /// Configure for a spherical extended source.
    /// @param inner_radius  Inner (void) radius for a hollow spherical shell.
    ///   0 = solid ball. The active material occupies the shell
    ///   [inner_radius, outer_radius]; the central void is non-attenuating.
    /// @param rotation  Stored for API symmetry; physically irrelevant for a
    ///   sphere (the volume is rotation-invariant).
    void configure_spherical(const Eigen::Vector3d& center, double inner_radius,
                             double outer_radius, const Eigen::Matrix3d& rotation);

    /// Configure for a Marinelli beaker source.
    /// All z-coordinates are absolute (in detector frame), pre-computed by
    /// EfficiencyCalculator from user-facing parameters relative to z_det_min.
    /// @param well_inner_radius  Inner radius of re-entrant well (cm)
    /// @param outer_radius       Outer radius of beaker body (cm)
    /// @param z_bk               Beaker well opening z (= z_det_min - endcap_to_beaker)
    /// @param z_we               Well end z (= z_det_min + well_depth)
    /// @param z_bot              Beaker bottom z (= z_bk - fill_height)
    void configure_marinelli(double well_inner_radius, double outer_radius,
                             double z_bk, double z_we, double z_bot);

    /// Compute transmission from an emission point through source material + shields.
    /// @param position   Photon emission point (detector frame)
    /// @param direction  Photon direction (normalized)
    /// @param energy_MeV Photon energy
    /// @return Transmission factor (0 to 1)
    double compute_transmission(const Eigen::Vector3d& position,
                                const Eigen::Vector3d& direction,
                                double energy_MeV) const;

    /// Point-source shortcut: direction-independent transmission (spherical shells).
    /// Does not include source material (point sources have no volume).
    double point_source_transmission(double energy_MeV) const;

    /// Outermost extent radius (for overlap checking).
    double outermost_extent_radius() const;

    /// Result of FEP-only source transmission with stochastic Rayleigh.
    struct SourceTransmissionResult {
        double weight;                  ///< Transmission weight (product of exp(-mu_no_rs * L))
        Eigen::Vector3d direction;      ///< May be modified by Rayleigh scatters
        Eigen::Vector3d exit_position;  ///< 3D point where photon exits outermost shield
    };

    /// Compute transmission through source material + shields with stochastic Rayleigh.
    /// Same as compute_transmission() but treats Rayleigh elastically: weights by
    /// exp(-mu_no_rs * path) and stochastically applies Rayleigh direction changes.
    SourceTransmissionResult compute_transmission_fep_only(
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& direction,
        double energy_MeV,
        std::mt19937_64& rng) const;

    /// A secondary photon produced in source material (e.g., 511 keV from PP).
    struct SourceSecondaryPhoton {
        Eigen::Vector3d position;   ///< Exit position from source material
        Eigen::Vector3d direction;  ///< Exit direction
        double energy_keV;          ///< Photon energy
    };

    /// A Compton recoil electron that escaped source material toward the detector.
    /// Produced when source_electron_enabled_ is true and the electron passes
    /// geometric + range filtering (energy > 50 keV, aimed at detector, CSDA range
    /// sufficient to traverse remaining source + wall material).
    struct SourceElectron {
        Eigen::Vector3d position;   ///< Exit position (just past source geometry boundary)
        Eigen::Vector3d direction;  ///< Electron direction
        double energy_keV;          ///< Residual KE at exit (after CSDA through source + wall)
    };

    /// Result of full MC transport through source material + shielding.
    struct SourceFullTransportResult {
        bool survived;              ///< Photon not absorbed
        Eigen::Vector3d position;   ///< Exit position
        Eigen::Vector3d direction;  ///< Exit direction (may have changed)
        double energy_keV;          ///< Exit energy (may have decreased from Compton)

        /// Secondary photons that exited the source geometry: 511 keV
        /// annihilation gammas from PP, and bremsstrahlung from electrons
        /// stopping in source material/shielding.
        std::vector<SourceSecondaryPhoton> secondaries;

        /// Compton recoil electrons that escaped source material toward the detector.
        /// Only populated when source electron transport is enabled.
        std::vector<SourceElectron> source_electrons;

        /// True if the photon underwent ANY interaction (Compton, PE, PP, or
        /// Rayleigh) in the source geometry. False means the photon exited
        /// with its original direction and energy and produced no
        /// secondaries/electrons (the "unscattered" event class).
        bool interacted = false;

        /// Product of angular importance-sampling weight factors from
        /// Compton-angle biasing (see ComptonBiasConfig). 1.0 when biasing
        /// is off. Must be multiplied into the EVENT weight by the caller
        /// (one weight per event, applied once at the tally).
        double bias_weight = 1.0;
    };

    /// Compton-angle mixture biasing for the primary photon's first few
    /// Compton vertices in the source geometry: the scatter direction is
    /// sampled from
    ///   q(w) = (1-gamma) * p_analog(w) + gamma * 1[w in detector cone]/Omega_d
    /// where p_analog is the KN x S(x,Z) angular density on the sphere and
    /// the cone subtends the detector bounding sphere from the scatter
    /// point. The event weight gains p/q <= 1/(1-gamma) per biased vertex.
    /// Unbiased (path-space importance sampling); feeds the scattered-in
    /// component of total efficiency. Applied only to the primary photon
    /// (never to annihilation/brems recursion) and never for Marinelli.
    /// Requires set_detector_bounds().
    struct ComptonBiasConfig {
        double cone_fraction = 0.0;   ///< gamma in [0, 0.9]; 0 = disabled
        /// Bias only the first K Compton vertices. K = 1 keeps every
        /// angular-norm evaluation at the (constant) primary energy, where
        /// it is cached — later vertices have continuously varying energies
        /// whose per-call norm quadrature would dominate the vertex cost.
        int max_biased_vertices = 1;
    };
    void set_compton_bias(const ComptonBiasConfig& cfg) { compton_bias_ = cfg; }
    const ComptonBiasConfig& compton_bias() const { return compton_bias_; }

    /// Probability that a photon emitted at `position` along `direction`
    /// crosses the entire source geometry with ZERO interactions:
    ///   T = exp(-sum_seg mu_total * length)
    /// over trace_source_segments(), with the same unmasked mu_total the
    /// analog sampler in transport_source_photon() uses (including
    /// Rayleigh), so T is exactly P(interacted == false) for that sampler.
    /// Used by the two-stream estimator's direct stream.
    ///
    /// @param path_cm  Optional out: total geometric path length (cm)
    ///                 through the source geometry along the ray (for air
    ///                 attenuation bookkeeping).
    double no_interaction_probability(const Eigen::Vector3d& position,
                                      const Eigen::Vector3d& direction,
                                      double energy_keV,
                                      double* path_cm = nullptr) const;

    /// Full analog MC transport through source material + shielding.
    /// Tracks the photon through source material and concentric shielding shells,
    /// handling PE absorption, Compton scattering (with energy loss), and Rayleigh
    /// scattering (direction change only). After any scatter, re-traces from the
    /// new position/direction.
    ///
    /// @param position   Photon emission point (detector frame)
    /// @param direction  Photon direction (normalized)
    /// @param energy_keV Photon energy in keV
    /// @param rng        Random number generator
    /// @return           Transport result with survived flag, exit position/direction/energy
    SourceFullTransportResult transport_source_photon(
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& direction,
        double energy_keV,
        std::mt19937_64& rng) const;

    /// Result of checking whether a photon from outside can re-enter the
    /// Marinelli water volume (e.g., after exiting the crystal).
    struct MarinelliReentryInfo {
        bool can_reenter = false;
        double distance_to_water = 0.0;       ///< Distance from pos to water entry
        double wall_path = 0.0;               ///< Path through beaker wall (0 if no wall, e.g., well opening)
        Eigen::Vector3d water_entry_pos{0,0,0}; ///< Position at water entry (just inside water)
    };

    /// Check if a photon from outside the Marinelli water can re-enter.
    /// Computes the nearest intersection with the water boundary, accounting
    /// for beaker wall traversal. Used for environmental re-scattering loop.
    MarinelliReentryInfo compute_marinelli_reentry(
        const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) const;

    /// Set detector bounding geometry for source electron filtering.
    /// Called by EfficiencyCalculator after detector configuration.
    void set_detector_bounds(double det_radius, double det_z_min, double det_z_max);

    /// Enable/disable source electron transport (Compton recoil electrons
    /// from source material tracked through source geometry toward detector).
    void set_source_electron_transport(bool enable) { source_electron_enabled_ = enable; }

    /// Set minimum electron KE threshold (keV) for source electron tracking.
    void set_source_electron_threshold(double keV) { source_electron_threshold_keV_ = keV; }

    /// Enable/disable bremsstrahlung from electrons (Compton recoil,
    /// photoelectron, pair e+/e-) slowing down in source material/shielding.
    /// Brems photons are transported through the remaining source geometry and
    /// survivors join SourceFullTransportResult::secondaries. Default on
    /// (physics, not variance reduction); the off switch is for A/B
    /// quantification.
    void set_source_brems(bool enable) { source_brems_enabled_ = enable; }
    bool source_brems_enabled() const { return source_brems_enabled_; }

    /// Enable/disable geometric pre-check (electron ray must hit detector bounding cylinder).
    void set_source_electron_geom_check(bool enable) { source_electron_geom_check_ = enable; }

    /// Whether any source material or shields are configured.
    bool has_source_effects() const {
        return source_material_ != nullptr || !shields_.empty();
    }

    bool is_configured() const { return configured_; }

    // Shape query for GDML export.
    enum class Shape { Point, Cylindrical, Rectangular, Marinelli, Sphere };
    Shape shape() const { return shape_; }

    // Accessors for GDML export.
    const Eigen::Vector3d& point_position() const { return point_position_; }
    const Eigen::Vector3d& cyl_center() const { return cyl_center_; }
    double cyl_radius() const { return cyl_radius_; }
    double cyl_inner_radius() const { return cyl_inner_r_; }
    double cyl_half_length() const { return cyl_half_length_; }
    const Eigen::Vector3d& rect_center() const { return rect_center_; }
    const Eigen::Vector3d& rect_half_dims() const { return rect_half_dims_; }
    const Eigen::Vector3d& rect_inner_half_dims() const { return rect_inner_half_dims_; }
    const Eigen::Vector3d& sphere_center() const { return sphere_center_; }
    double sphere_inner_radius() const { return sphere_inner_r_; }
    double sphere_radius() const { return sphere_radius_; }

    // Marinelli accessors (absolute z-coordinates in detector frame).
    double marinelli_well_inner_radius() const { return marinelli_well_r_; }
    double marinelli_outer_radius() const { return marinelli_outer_r_; }
    double marinelli_z_bk() const { return marinelli_z_bk_; }   ///< Beaker well opening z
    double marinelli_z_we() const { return marinelli_z_we_; }   ///< Well end z
    double marinelli_z_bot() const { return marinelli_z_bot_; } ///< Beaker bottom z

private:
    const Material* source_material_ = nullptr;
    std::vector<SourceShieldLayer> shields_;

    bool configured_ = false;

    // Bremsstrahlung from electrons stopping in source material/shielding
    bool source_brems_enabled_ = true;

    // Source electron transport
    bool source_electron_enabled_ = false;
    double source_electron_threshold_keV_ = 50.0;
    bool source_electron_geom_check_ = false;
    double det_bound_radius_ = 0.0;
    double det_bound_z_min_ = 0.0;
    double det_bound_z_max_ = 0.0;

    // Compton-angle mixture biasing (primary photon only)
    ComptonBiasConfig compton_bias_;

    /// transport_source_photon() body. bias_budget = number of Compton
    /// vertices that may still be angle-biased; the public entry point
    /// passes compton_bias_.max_biased_vertices (0 when disabled or
    /// Marinelli), recursive calls for annihilation gammas and brems
    /// photons always pass 0.
    SourceFullTransportResult transport_source_photon_impl(
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& direction,
        double energy_keV,
        std::mt19937_64& rng,
        int bias_budget) const;

    Shape shape_ = Shape::Point;

    // Point source
    Eigen::Vector3d point_position_{0, 0, 0};

    // Cylindrical source
    Eigen::Vector3d cyl_center_{0, 0, 0};
    double cyl_radius_ = 0.0;
    double cyl_inner_r_ = 0.0;   ///< Inner bore radius (0 = solid); annular when > 0
    double cyl_half_length_ = 0.0;
    Eigen::Matrix3d cyl_rotation_ = Eigen::Matrix3d::Identity();

    // Rectangular source
    Eigen::Vector3d rect_center_{0, 0, 0};
    Eigen::Vector3d rect_half_dims_{0, 0, 0};
    Eigen::Vector3d rect_inner_half_dims_{0, 0, 0};  ///< Inner void half-dims (all-zero = solid); box shell when > 0
    Eigen::Matrix3d rect_rotation_ = Eigen::Matrix3d::Identity();

    // Spherical source
    Eigen::Vector3d sphere_center_{0, 0, 0};
    double sphere_inner_r_ = 0.0;  ///< Inner void radius (0 = solid ball); shell when > 0
    double sphere_radius_ = 0.0;   ///< Outer radius
    Eigen::Matrix3d sphere_rotation_ = Eigen::Matrix3d::Identity();

    // Marinelli beaker (absolute z-coordinates in detector frame)
    double marinelli_well_r_ = 0.0;   ///< Inner radius of re-entrant well
    double marinelli_outer_r_ = 0.0;  ///< Outer radius of beaker body
    double marinelli_z_bk_ = 0.0;     ///< Beaker well opening z
    double marinelli_z_we_ = 0.0;     ///< Well end z (= z_det_min + well_depth)
    double marinelli_z_bot_ = 0.0;    ///< Beaker bottom z (= z_bk - fill_height)

    /// Test if a point is inside the Marinelli L-shaped sample volume.
    bool is_inside_marinelli(const Eigen::Vector3d& pos) const;

    /// Which surface a photon exits through in the Marinelli volume.
    enum class MarinelliExitSurface {
        OuterCylinder,  ///< Side wall (r = outer_r)
        WellWall,       ///< Well inner wall (r = well_r)
        Bottom,         ///< Bottom plane (z = z_bot)
        Top,            ///< Top of ring (z = z_we)
        WellOpening,    ///< Open top (z = z_bk, r < well_r) — no beaker wall
        None            ///< Fallback (shouldn't happen)
    };

    /// Exit info from Marinelli boundary: distance and which surface.
    struct MarinelliExitInfo {
        double distance;
        MarinelliExitSurface surface;
    };

    /// Distance from a point inside the Marinelli volume to exit boundary.
    double marinelli_boundary_distance(const Eigen::Vector3d& pos,
                                       const Eigen::Vector3d& dir) const;

    /// Distance + exit surface from a point inside the Marinelli volume.
    MarinelliExitInfo marinelli_exit_info(const Eigen::Vector3d& pos,
                                          const Eigen::Vector3d& dir) const;

    /// Compute angle correction factor for beaker wall at exit surface.
    /// Returns 0.0 if photon exits through well opening (no wall) or is outside sample.
    /// Otherwise returns 1/|cos(incidence)| capped at 20 (grazing incidence limit).
    double marinelli_wall_factor(const Eigen::Vector3d& position,
                                  const Eigen::Vector3d& direction) const;

    /// Path through source material from emission point to boundary.
    double source_material_path(const Eigen::Vector3d& position,
                                const Eigen::Vector3d& direction) const;

public:
    /// Exact (or conservatively under-estimated) distance from an interior
    /// point to the nearest point of the source geometry's OUTER boundary
    /// (source material + shield envelope). Used to cheaply rule out charged-
    /// particle escape: a particle whose CSDA range is shorter than this distance
    /// cannot reach the boundary in any direction (straight-line displacement
    /// <= path length <= range). Returns a value <= the true distance, so a
    /// `range < min_distance_to_boundary` test never wrongly concludes "cannot
    /// escape". Public so the cascade β⁺ positron containment fast-path
    /// (EfficiencyCalculator) can reuse it.
    ///
    /// @param include_shields  When true, measures to the outer shield surface
    ///   (the source-geometry boundary). When false, measures to the
    ///   source-material surface only (the material/shield interface) — used to
    ///   test whether an electron born in the source material stays within it.
    ///   Returns 0.0 for Point (no material volume) when false.
    double min_distance_to_boundary(const Eigen::Vector3d& position,
                                    bool include_shields = true) const;

    /// A segment through a single material layer (for source transport).
    struct SourcePathSegment {
        const Material* material;
        double length;  ///< Path length through this segment (cm)
    };

    /// Build path segments through source material + shielding from a given
    /// position/direction. Used by transport_source_photon() and by
    /// ElectronCsda::walk_in_source_geometry().
    std::vector<SourcePathSegment> trace_source_segments(
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& direction,
        double energy_keV) const;

    /// Reusable-output overload of trace_source_segments: fills `out` (cleared
    /// first) instead of heap-allocating a fresh vector per call. `max_segments`
    /// stops the trace once that many segments have been recorded — the Molière
    /// source-escape walk (ElectronCsda::walk_in_source_geometry) passes 1, since
    /// it needs only the first (current material + distance to next interface)
    /// segment and must not pay to trace downstream shells every substep.
    /// Bit-identical to the by-value form: the first `max_segments` segments are
    /// produced by exactly the same arithmetic in the same order (this function
    /// makes no RNG draws).
    void trace_source_segments(
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& direction,
        double energy_keV,
        std::vector<SourcePathSegment>& out,
        std::size_t max_segments = SIZE_MAX) const;

private:

};

} // namespace ceelo
