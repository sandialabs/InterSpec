#ifndef CEELO_GEOMETRY_GEOMETRY_DESCRIPTOR_H
#define CEELO_GEOMETRY_GEOMETRY_DESCRIPTOR_H
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

/// @file GeometryDescriptor.h
/// @brief The storable detector geometry: crystal, fillet, bore, dead layer,
/// housing layers, collimator and the materials they index. Rebuilds the
/// ray-trace Geometry (build_geometry) and declares it onto any Geometry-like
/// sink (apply_detector_side). Always built; the XML codec lives with
/// io/DetectorResponse's.

#include "geometry/ConvexSeparation.h"
#include "geometry/Geometry.h"
#include "materials/Material.h"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace ceelo {

/// Full material definition so the stored response is self-contained
/// (user-defined materials round-trip; nothing resolves by name at load).
struct MaterialSpec {
    std::string name;
    double density_g_per_cm3 = 0.0;
    std::vector<MaterialComponent> composition;  ///< {Z, mass_fraction}

    static MaterialSpec from(const Material& m);
    Material to_material() const { return Material(name, density_g_per_cm3, composition); }
};

/// One concentric cup-shaped layer (endcap, window, housing, wrap...).
/// Multiple layers are supported, ordered innermost -> outermost; each maps
/// to one Geometry::add_attenuator call, which builds it as a cup over the
/// envelope inside it (the crystal, then each inner layer): front and side
/// walls, plus a BACK PLATE wherever z_end reaches past that envelope's back.
/// So a layer's back thickness is `z_end_cm` minus the inner envelope's z_max
/// (see GeometryDescriptor::layer_back_thickness_cm); a vacuum layer with a
/// back thickness is a gap behind the crystal, and lets the outer walls run on
/// past it.  A layer ending at the crystal's back (z_end = L) has no back.
struct LayerSpec {
    int material_index = -1;      ///< into GeometryDescriptor::materials
    double front_thickness_cm = 0.0;
    double side_thickness_cm = 0.0;
    double z_start_cm = 0.0;      ///< axial extent (crystal front face = 0)
    double z_end_cm = 0.0;
};

/// Side-only collimator tube (maps to Geometry::add_collimator).
struct CollimatorSpec {
    int material_index = -1;
    double side_thickness_cm = 0.0;
    double z_start_cm = 0.0;      ///< typically negative (extends past the face)
    double z_end_cm = 0.0;
};

enum class ResponseSymmetry : uint8_t {
    Axial,     ///< cylinders: no phi axis
    Quadrant   ///< boxes: phi axis over one quadrant (reflect by symmetry)
};

/// Which plane a caller's `dist_cm` is measured FROM. Distances and the kernel
/// must share one origin -- a mismatch silently corrupts everything (master
/// plan sec 3.2 warning) rather than failing.
///
/// Both data below are planes perpendicular to the axis, and BOTH are measured
/// to the source, so a larger value always means further away:
///
///   CrystalFace  -- the front face of the CRYSTAL SOLID, z = 0 in the crystal
///                   frame. NOT the front of the active volume: with a dead
///                   layer the active volume starts a further t_front behind
///                   this plane (see Geometry::set_dead_layer).
///   EndcapFront  -- the outermost front surface of the attenuator stack, i.e.
///                   the face a user can physically touch. This is InterSpec's
///                   user-facing convention. It sits
///                   endcap_front_offset_cm() in FRONT of CrystalFace.
///
/// The conversion between them is GeometryDescriptor::endcap_front_offset_cm(),
/// applied in DetectorResponse::query_position(). Nothing else should be doing
/// that arithmetic by hand.
enum class ReferencePoint : uint8_t {
    CrystalFace,
    EndcapFront
};

/// Ways a GeometryDescriptor violates a Geometry/RayTrace precondition.
///
/// Those preconditions are asserts, so they vanish in a release build and a
/// violating descriptor traces silent garbage instead of failing. Every
/// descriptor must therefore be run through GeometryDescriptor::problems()
/// before it reaches a Geometry.
enum class GeometryProblem : uint8_t {
    DimensionsMissing,     ///< dimensions_cm too short for the shape
    BulletOnNonCylinder,   ///< fillet requested on a box
    BulletNotFinite,       ///< NaN / inf / negative
    BulletTooWide,         ///< bullet_radius_cm >= crystal radius (rho_c <= 0)
    BulletTooLong,         ///< bullet_radius_cm >= crystal length
    BulletNoDeadLayerRoom, ///< dead layer leaves no active fillet
    DeadLayerTooThick,     ///< dead layer consumes the whole crystal
    BoreOnNonCylinder,
    BoreNotFinite,
    BoreTooWide,           ///< bore radius >= crystal radius
    BoreTooDeep,           ///< bore depth >= crystal length
    BoreTipTooBlunt,       ///< rounded_tip && radius > depth
    BoreOutsideFillet,     ///< bore_fits() fails against the fillet
    BoreInsideDeadLayer    ///< bore radius >= (crystal radius - side dead layer)
};

/// Short English description of `p`; for developer-facing messages (InterSpec
/// maps the enum onto its own localized strings).
const char* to_string(GeometryProblem p);

/// The storable geometry: shape, dimensions, bore, dead layer, layers,
/// collimator, symmetry, reference point + the material table they index.
struct GeometryDescriptor {
    DetectorShape shape = DetectorShape::Cylinder;
    /// The serialized crystal dimensions, in cm: Cylinder {radius, FULL length};
    /// Box {half_x, half_y, FULL length}.  See CRYSTAL DIMENSION CONVENTION in
    /// geometry/Geometry.h -- and prefer set_dimensions() / cylinder_dims() /
    /// box_dims() below, which name the fields, over indexing this by hand.
    ///
    /// This layout is written into every response file; it does not change.
    std::vector<double> dimensions_cm;
    /// Cylinder only: quarter-torus fillet radius on the outer FRONT edge
    /// ("bulletization", ANGLE's `bulletizingRadius`), cm. 0 = a sharp
    /// 90-degree edge -- the default, and bit-for-bit the pre-feature trace.
    /// Must be < radius and < length; see problems().
    double bullet_radius_cm = 0.0;
    int crystal_material_index = -1;
    std::optional<BoreHoleConfig> bore;          ///< coax HPGe finger
    std::optional<DeadLayerConfig> dead_layer;   ///< cm, crystal material
    std::vector<LayerSpec> layers;               ///< innermost -> outermost
    std::optional<CollimatorSpec> collimator;
    ResponseSymmetry symmetry = ResponseSymmetry::Axial;
    ReferencePoint reference_point = ReferencePoint::CrystalFace;
    std::vector<MaterialSpec> materials;

    /// Set `shape` and `dimensions_cm` together, so a Box can never end up
    /// declared with two numbers.  See CRYSTAL DIMENSION CONVENTION.
    void set_dimensions(const CylinderDims& dims);
    void set_dimensions(const BoxDims& dims);

    /// Named read-back of `dimensions_cm`.  Asserts the shape matches and the
    /// vector is long enough; consult problems() first for untrusted input.
    CylinderDims cylinder_dims() const;
    BoxDims box_dims() const;

    /// Build the ray-trace geometry. The returned Geometry references the
    /// Material instances appended to `owned` -- keep them alive as long as
    /// the Geometry (DetectorResponse holds both).
    Geometry build_geometry(std::vector<std::unique_ptr<Material>>& owned) const;

    /// Transverse half-extent `a` -- the near/far regime SCALE PARAMETER
    /// (cylinder: outer radius; box: half-diagonal), not a physical dimension.
    ///
    /// It intentionally includes the side dead layer, even though the dead
    /// layer is internal to the crystal (Geometry::set_dead_layer), so it runs
    /// slightly large. That is harmless and must stay: `a` is computed from the
    /// descriptor at BOTH generation time (ResponseGenerator picks its far
    /// reference distance as far_distance_a * a) and query time (the
    /// near-regime test is d < near_regime_a * a). Changing the value would
    /// desynchronise every already-generated response from the code that made
    /// it. Contrast endcap_front_offset_cm(), which is query-time only and so
    /// could be -- and was -- corrected in place.
    ///
    /// If you need a true physical outer radius, sum the layer side
    /// thicknesses onto dimensions_cm yourself; do not use this.
    double transverse_half_extent() const;

    /// |z_min| of the outermost front surface: the (positive) offset from the
    /// endcap front to the CRYSTAL SOLID's face (z = 0), i.e. the summed front
    /// thicknesses of the attenuator shells. The dead layer is deliberately
    /// excluded -- it is carved out of the inside of the crystal, so it sets
    /// where the ACTIVE volume starts, not where the crystal face is.
    double endcap_front_offset_cm() const;

    /// The outer surface of the detector as described - crystal, layers and collimator - in the
    /// crystal frame (cm): the solid a source must stay outside of.  Unlike
    /// transverse_half_extent(), the dead layer (inside the crystal) adds nothing.
    struct Envelope {
        DetectorShape shape = DetectorShape::Cylinder;
        double radius_cm = 0.0;     ///< Cylinder: outer radius
        double half_x_cm = 0.0;     ///< Box: outer half-widths
        double half_y_cm = 0.0;
        double z_min_cm = 0.0;      ///< front-most outer surface (<= 0)
        double z_max_cm = 0.0;      ///< back-most outer surface (>= crystal length)

        /// Distance from the crystal-face origin, along the unit direction `u`, to where the ray
        /// leaves this envelope (cm).  The origin is always inside it.  A source at distance d
        /// along `u` has a CLEARANCE of d - exit_distance_cm(u) beyond the detector's outside.
        double exit_distance_cm(const Eigen::Vector3d& u) const;

        /// The point along the unit direction `u` from the crystal-face origin whose distance to
        /// this envelope - to its nearest point, as solid().project() finds it - is
        /// `clearance_cm` (>= 0).  Past the exit point that distance only grows (the envelope is
        /// convex), so it is found by bisection.  Beside a long detector it is closer to the
        /// envelope than exit_distance_cm(u) + clearance_cm, which measures along the ray.
        Eigen::Vector3d position_at_clearance(const Eigen::Vector3d& u, double clearance_cm) const;

        /// This envelope as a solid in the crystal frame (cm), for convex_separation() from a
        /// source and convex_entry_distance() along a ray; place it elsewhere with
        /// ConvexSolid::transformed().
        ConvexSolid solid() const;
    };
    Envelope outer_envelope() const;

    /// Back plate thickness of layer `k` (cm): how far its z_end reaches past the envelope inside
    /// it - 0 when the layer ends at, or before, that envelope's back.
    double layer_back_thickness_cm(size_t k) const;

    /// Gives layer `k` a back plate `t_cm` thick (>= 0) by moving its z_end, and moves every outer
    /// layer's z_end (and the collimator's, if it reached the old back) by the same amount, so
    /// their own back plates are unchanged.
    void set_layer_back_thickness_cm(size_t k, double t_cm);

    /// How far the described structure reaches behind the crystal's back face (cm): the sum of the
    /// layers' back thicknesses, vacuum gaps included.  0 when nothing behind the crystal is
    /// described - the response then knows nothing of what really sits there.
    double back_structure_extent_cm() const;

    /// Crystal K-edges within (e_min, e_max) keV -- the mandatory segment
    /// breaks for any ln-eta(E) interpolation (both flanks become nodes).
    std::vector<double> crystal_k_edges(double e_min_keV,
                                        double e_max_keV) const;

    /// Every Geometry/RayTrace precondition this descriptor violates; empty
    /// means it is safe to build. See GeometryProblem.
    std::vector<GeometryProblem> problems() const;

    /// Standalone XML for JUST the geometry, so a host can store a detector's
    /// shape before (or without) any response having been generated for it.
    ///
    /// The payload is the same <Detector> element DetectorResponse's own codec
    /// writes -- shared code, so the two cannot drift -- wrapped in a
    /// <CeeLoGeometry> root. Round-trips through from_xml_string().
    std::string to_xml_string() const;
    static GeometryDescriptor from_xml_string(const std::string& xml);

    /// Plain-text form of every field above, "CEELO_GEOMETRY 1" then one
    /// keyword line per item (see GeometryDescriptor.cpp). Always built and
    /// dependency-free, unlike the XML pair, which needs rapidxml and is
    /// defined only with CEELO_BUILD_RESPONSE. Round-trips through
    /// from_text(), which throws std::runtime_error on malformed input and,
    /// like from_xml_string(), does not run problems().
    std::string to_text() const;
    static GeometryDescriptor from_text(const std::string& text);
};

/// Declare a descriptor's detector side - crystal, fillet, bore, dead layer, attenuator
/// layers, collimator - onto `sink`.
///
/// ONE definition, used by all three paths that build a detector from a descriptor:
/// GeometryDescriptor::build_geometry (the query-time kernel), ResponseGenerator's per-node
/// setup, and ResponseGenerator::configure_calculator.  `Geometry` and `EfficiencyCalculator`
/// expose the same setter signatures, so both are valid sinks.
///
/// Keep it that way.  The eta table is measured through a calculator configured here while
/// the query-time kernel K traces a Geometry configured here; if the two ever declare
/// different solids the response is internally inconsistent, and the generator's own probe
/// banks cannot detect it because they route through the same path.
///
/// `mat` maps a descriptor material index to an instantiated Material the caller owns and
/// keeps alive for as long as the sink is used.
template <class Sink, class MatFn>
void apply_detector_side(Sink& sink, const GeometryDescriptor& gd, MatFn mat) {
    sink.set_detector_from_dimensions_vector(gd.shape, mat(gd.crystal_material_index),
                                             gd.dimensions_cm);
    // set_detector() clears the fillet/bore/dead layer, so declare them after it; fillet
    // first, so bore_fits() sees the final crystal profile.
    if (gd.bullet_radius_cm > 0.0) sink.set_bullet_radius(gd.bullet_radius_cm);
    if (gd.bore)
        sink.set_bore_hole(gd.bore->radius, gd.bore->depth, gd.bore->rounded_tip);
    if (gd.dead_layer)
        sink.set_dead_layer(gd.dead_layer->front, gd.dead_layer->side, gd.dead_layer->back);
    for (const LayerSpec& l : gd.layers)
        sink.add_attenuator(mat(l.material_index), l.front_thickness_cm,
                            l.side_thickness_cm, l.z_start_cm, l.z_end_cm);
    if (gd.collimator)
        sink.add_collimator(mat(gd.collimator->material_index),
                            gd.collimator->side_thickness_cm, gd.collimator->z_start_cm,
                            gd.collimator->z_end_cm);
}

}  // namespace ceelo

#endif  // CEELO_GEOMETRY_GEOMETRY_DESCRIPTOR_H
