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

#include "geometry/GeometryDescriptor.h"

#include "cross_sections/CrossSectionData.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ceelo {

MaterialSpec MaterialSpec::from(const Material& m) {
    MaterialSpec s;
    s.name = m.name();
    s.density_g_per_cm3 = m.density();
    s.composition = m.composition();
    return s;
}

void GeometryDescriptor::set_dimensions(const CylinderDims& dims) {
    shape = DetectorShape::Cylinder;
    dimensions_cm = to_dimensions_vector(dims);
}

void GeometryDescriptor::set_dimensions(const BoxDims& dims) {
    shape = DetectorShape::Box;
    dimensions_cm = to_dimensions_vector(dims);
}

CylinderDims GeometryDescriptor::cylinder_dims() const {
    assert(shape == DetectorShape::Cylinder);
    return cylinder_dims_from_vector(dimensions_cm);
}

BoxDims GeometryDescriptor::box_dims() const {
    assert(shape == DetectorShape::Box);
    return box_dims_from_vector(dimensions_cm);
}

Geometry GeometryDescriptor::build_geometry(
    std::vector<std::unique_ptr<Material>>& owned) const {
    if (crystal_material_index < 0 ||
        crystal_material_index >= static_cast<int>(materials.size()))
        throw std::runtime_error("GeometryDescriptor: bad crystal material index");

    // Validate BEFORE rebuilding `owned`: finalize() assigns geometry_ only on
    // success, so throwing after the clear would leave a previously-finalized
    // response holding Material* into the vector we just freed. Geometry's own
    // preconditions are asserts, so in a release build this is the only thing
    // standing between a bad descriptor and a silently garbage trace.
    const std::vector<GeometryProblem> probs = problems();
    if (!probs.empty())
        throw std::runtime_error(std::string("GeometryDescriptor: ")
                                 + to_string(probs.front()));

    owned.clear();
    owned.reserve(materials.size());
    for (const MaterialSpec& spec : materials)
        owned.push_back(std::make_unique<Material>(spec.to_material()));

    auto mat_at = [&](int idx) -> const Material* {
        if (idx < 0 || idx >= static_cast<int>(owned.size()))
            throw std::runtime_error("GeometryDescriptor: bad material index");
        return owned[static_cast<size_t>(idx)].get();
    };

    Geometry g;
    apply_detector_side(g, *this, mat_at);
    return g;
}

const char* to_string(GeometryProblem p) {
    switch (p) {
        case GeometryProblem::DimensionsMissing:
            return "the crystal dimensions are missing (need radius+length, or 3 half-widths)";
        case GeometryProblem::DeadLayerTooThick:
            return "the dead layer consumes the whole crystal";
        case GeometryProblem::BulletOnNonCylinder:
            return "a front-edge fillet needs a cylindrical crystal";
        case GeometryProblem::BulletNotFinite:
            return "the front-edge fillet radius is not a finite, non-negative length";
        case GeometryProblem::BulletTooWide:
            return "the front-edge fillet radius is not less than the crystal radius";
        case GeometryProblem::BulletTooLong:
            return "the front-edge fillet radius is not less than the crystal length";
        case GeometryProblem::BulletNoDeadLayerRoom:
            return "the dead layer leaves no room for the front-edge fillet";
        case GeometryProblem::BoreOnNonCylinder:
            return "a bore hole needs a cylindrical crystal";
        case GeometryProblem::BoreNotFinite:
            return "the bore radius and depth are not finite, positive lengths";
        case GeometryProblem::BoreTooWide:
            return "the bore radius is not less than the crystal radius";
        case GeometryProblem::BoreTooDeep:
            return "the bore depth is not less than the crystal length";
        case GeometryProblem::BoreTipTooBlunt:
            return "a rounded bore tip needs the bore at least as deep as its radius";
        case GeometryProblem::BoreOutsideFillet:
            return "the bore ends where the filleted crystal is narrower than the bore";
        case GeometryProblem::BoreInsideDeadLayer:
            return "the bore radius is not less than the active (post dead layer) radius";
    }
    assert(0 && "unhandled GeometryProblem");
    return "invalid geometry";
}

std::vector<GeometryProblem> GeometryDescriptor::problems() const {
    std::vector<GeometryProblem> out;

    const bool cyl = (shape == DetectorShape::Cylinder);
    const double R = (cyl && dimensions_cm.size() > 0) ? dimensions_cm[0] : 0.0;
    const double L = cyl ? (dimensions_cm.size() > 1 ? dimensions_cm[1] : 0.0)
                         : (dimensions_cm.size() > 2 ? dimensions_cm[2] : 0.0);
    const double df = dead_layer ? dead_layer->front : 0.0;
    const double ds = dead_layer ? dead_layer->side : 0.0;
    const double db = dead_layer ? dead_layer->back : 0.0;
    const double rb = bullet_radius_cm;

    // --- crystal itself: set_detector_from_dimensions_vector() indexes
    //     dimensions_cm unchecked ---
    // Its only guard is assert(size() >= 2 / >= 3), so a short vector (a saved
    // response with a missing or truncated <Dimensions>) would read past the
    // end of the vector in a release build.
    const size_t need = cyl ? 2u : 3u;
    if (dimensions_cm.size() < need) {
        out.push_back(GeometryProblem::DimensionsMissing);
        return out;   // every check below reads R / L
    }

    // A dead layer at least as thick as the crystal leaves the tracer working
    // on a negative-extent solid. Checked here rather than only inside the
    // fillet branch, since it is wrong with or without a fillet. NOTE the
    // transverse extent is shape-dependent -- R is meaningless for a box, and
    // using it here would reject every box outright.
    const double half_extent = cyl ? R
                                   : std::min(dimensions_cm[0], dimensions_cm[1]);
    if (!std::isfinite(half_extent) || !std::isfinite(L)
        || half_extent <= 0.0 || L <= 0.0
        || (half_extent - ds) <= 0.0 || (L - df - db) <= 0.0)
        out.push_back(GeometryProblem::DeadLayerTooThick);

    // --- front-edge fillet: Geometry::set_bullet_radius asserts ---
    // Anything but an exact 0 is a request for a fillet, so negatives and NaN
    // get reported rather than silently reading as "sharp edge".
    if (rb != 0.0) {
        if (!cyl) {
            out.push_back(GeometryProblem::BulletOnNonCylinder);
        } else if (!std::isfinite(rb) || rb < 0.0) {
            out.push_back(GeometryProblem::BulletNotFinite);
        } else {
            if (rb >= R) out.push_back(GeometryProblem::BulletTooWide);
            if (rb >= L) out.push_back(GeometryProblem::BulletTooLong);
            // The active volume gets the fillet offset inward by the dead
            // layer; trace_cylinder_geometry asserts it still has room.
            const double dl_r = R - ds;
            const double rba = std::max(0.0, rb - std::max(df, ds));
            if (!(dl_r > 0.0 && rba < dl_r && (df + rba) <= (L - db)))
                out.push_back(GeometryProblem::BulletNoDeadLayerRoom);
        }
    }

    // --- bore: Geometry::set_bore_hole asserts, plus the active-radius one
    // that RayTrace relies on but no assert covers ---
    if (bore) {
        const double br = bore->radius;
        const double bd = bore->depth;
        if (!cyl) {
            out.push_back(GeometryProblem::BoreOnNonCylinder);
        } else if (!std::isfinite(br) || !std::isfinite(bd)
                   || br <= 0.0 || bd <= 0.0) {
            out.push_back(GeometryProblem::BoreNotFinite);
        } else {
            if (br >= R) out.push_back(GeometryProblem::BoreTooWide);
            if (bd >= L) out.push_back(GeometryProblem::BoreTooDeep);
            if (bore->rounded_tip && br > bd)
                out.push_back(GeometryProblem::BoreTipTooBlunt);
            // Only meaningful once the bore is narrower than the crystal:
            // bore_fits() returns false for br >= R too, and blaming that on
            // the fillet would make relaxGeometryFeatures() drop a perfectly
            // good fillet for a problem it cannot fix.
            if (br < R && rb > 0.0 && std::isfinite(rb)) {
                // The outer solid, and then the dead-layer-offset active solid
                // -- which carries a tighter fillet (radius rb - max(df,ds)
                // about a centre at z = df), so it can be violated when the
                // outer one is not. trace_cylinder_geometry cuts the bore out
                // of that one.
                const double act_r = R - ds;
                const double act_len = L - df - db;
                const double act_rb = std::max(0.0, rb - std::max(df, ds));
                if (!bore_fits(br, bd, R, L, rb)
                    || (act_r > 0.0 && act_len > 0.0
                        && !bore_fits(br, bd - db, act_r, act_len, act_rb)))
                    out.push_back(GeometryProblem::BoreOutsideFillet);
            }
            // trace_cylinder_geometry cuts the bore out of the ACTIVE cylinder
            // (radius R - side dead layer), which set_bore_hole never checks.
            // No assert covers this: such a bore simply yields zero active
            // volume. Rejecting is a deliberate policy choice -- a response
            // that silently computes zero efficiency is worse than one that
            // refuses to load -- so unlike the fillet it is NOT relaxable.
            if (br >= (R - ds))
                out.push_back(GeometryProblem::BoreInsideDeadLayer);
        }
    }

    return out;
}

double GeometryDescriptor::transverse_half_extent() const {
    double r = 0.0;
    if (shape == DetectorShape::Cylinder) {
        r = dimensions_cm.size() > 0 ? dimensions_cm[0] : 0.0;
    } else {
        const double hx = dimensions_cm.size() > 0 ? dimensions_cm[0] : 0.0;
        const double hy = dimensions_cm.size() > 1 ? dimensions_cm[1] : 0.0;
        r = std::hypot(hx, hy);
    }
    double side = 0.0;
    if (dead_layer) side += dead_layer->side;
    for (const LayerSpec& l : layers) side += l.side_thickness_cm;
    if (collimator) side += collimator->side_thickness_cm;
    return r + side;
}

double GeometryDescriptor::endcap_front_offset_cm() const {
    double front = 0.0;
    // ONLY the attenuator shells: they are what stands between the endcap front
    // and the crystal. The dead layer must NOT be added -- it is carved out of
    // the INSIDE of the crystal solid (trace_cylinder_geometry insets the
    // active volume within [0, L]; the solid still starts at z = 0), so it does
    // not push the crystal face away from the endcap. Adding it made
    // query_position() place every EndcapFront-referenced source one dead-layer
    // thickness too far away -- ~0.7 mm on a typical HPGe, worth ~3% of the
    // efficiency at contact geometry, and nothing at all in the far field.
    for (const LayerSpec& l : layers) front += l.front_thickness_cm;
    // A collimator extending past the face does NOT define the endcap front
    // (it is an open tube); the reference plane is the front surface stack.
    return front;
}

namespace {
/// z_max of the envelope inside layer `k`: the crystal's back, grown by each inner layer's
/// back plate (RayTrace's cup construction: shell z_max = max(inner z_max, z_end)).
double inner_back_z(const GeometryDescriptor& gd, const size_t k) {
    double z = (gd.dimensions_cm.empty() ? 0.0 : gd.dimensions_cm.back());
    for (size_t i = 0; (i < k) && (i < gd.layers.size()); ++i)
        z = std::max(z, gd.layers[i].z_end_cm);
    return z;
}

/// A back plate thinner than this is rounding residue (from shifting z_end values in place), not
/// structure.
constexpr double kMinBackPlateCm = 1.0e-9;

double back_plate(const double t_cm) {
    return (t_cm > kMinBackPlateCm) ? t_cm : 0.0;
}
}  // namespace

GeometryDescriptor::Envelope GeometryDescriptor::outer_envelope() const {
    Envelope e;
    e.shape = shape;
    const bool cyl = (shape == DetectorShape::Cylinder);
    double r  = cyl ? (dimensions_cm.size() > 0 ? dimensions_cm[0] : 0.0) : 0.0;
    double hx = cyl ? 0.0 : (dimensions_cm.size() > 0 ? dimensions_cm[0] : 0.0);
    double hy = cyl ? 0.0 : (dimensions_cm.size() > 1 ? dimensions_cm[1] : 0.0);
    double z_min = 0.0;
    double z_max = dimensions_cm.empty() ? 0.0 : dimensions_cm.back();
    for (const LayerSpec& l : layers) {
        r += l.side_thickness_cm;
        hx += l.side_thickness_cm;
        hy += l.side_thickness_cm;
        z_min = std::min(z_min - l.front_thickness_cm, l.z_start_cm);
        z_max = std::max(z_max, l.z_end_cm);
    }
    if (collimator) {
        r += collimator->side_thickness_cm;
        hx += collimator->side_thickness_cm;
        hy += collimator->side_thickness_cm;
        z_min = std::min(z_min, collimator->z_start_cm);
        z_max = std::max(z_max, collimator->z_end_cm);
    }
    e.radius_cm = r;
    e.half_x_cm = hx;
    e.half_y_cm = hy;
    e.z_min_cm = z_min;
    e.z_max_cm = z_max;
    return e;
}

double GeometryDescriptor::Envelope::exit_distance_cm(const Eigen::Vector3d& u) const {
    constexpr double inf = std::numeric_limits<double>::infinity();
    double t = inf;
    if (u.z() > 0.0) t = std::min(t, z_max_cm / u.z());
    if (u.z() < 0.0) t = std::min(t, z_min_cm / u.z());
    if (shape == DetectorShape::Cylinder) {
        const double uxy = std::hypot(u.x(), u.y());
        if (uxy > 0.0) t = std::min(t, radius_cm / uxy);
    } else {
        if (u.x() != 0.0) t = std::min(t, half_x_cm / std::abs(u.x()));
        if (u.y() != 0.0) t = std::min(t, half_y_cm / std::abs(u.y()));
    }
    return t;
}

ConvexSolid GeometryDescriptor::Envelope::solid() const {
    const Eigen::Vector3d c(0.0, 0.0, 0.5 * (z_min_cm + z_max_cm));
    const double half_length = 0.5 * (z_max_cm - z_min_cm);
    if (shape == DetectorShape::Cylinder)
        return ConvexSolid::cylinder(c, Eigen::Vector3d::UnitZ(), radius_cm, half_length);
    return ConvexSolid::box(c, Eigen::Vector3d(half_x_cm, half_y_cm, half_length));
}

double GeometryDescriptor::layer_back_thickness_cm(const size_t k) const {
    if (k >= layers.size())
        return 0.0;
    return back_plate(layers[k].z_end_cm - inner_back_z(*this, k));
}

void GeometryDescriptor::set_layer_back_thickness_cm(const size_t k, const double t_cm) {
    if (k >= layers.size())
        throw std::out_of_range("set_layer_back_thickness_cm: no such layer");
    if (!(t_cm >= 0.0))
        throw std::invalid_argument("set_layer_back_thickness_cm: thickness must be >= 0");

    const double old_back = std::max(inner_back_z(*this, k), layers[k].z_end_cm);
    const double new_end = inner_back_z(*this, k) + t_cm;
    const double delta = new_end - old_back;
    layers[k].z_end_cm = new_end;
    for (size_t i = k + 1; i < layers.size(); ++i)
        layers[i].z_end_cm += delta;
    if (collimator && (collimator->z_end_cm >= old_back))
        collimator->z_end_cm += delta;
}

double GeometryDescriptor::back_structure_extent_cm() const {
    const double crystal_back = dimensions_cm.empty() ? 0.0 : dimensions_cm.back();
    return back_plate(inner_back_z(*this, layers.size()) - crystal_back);
}

std::vector<double> GeometryDescriptor::crystal_k_edges(double e_min_keV,
                                                        double e_max_keV) const {
    std::vector<double> edges;
    if (crystal_material_index < 0 ||
        static_cast<size_t>(crystal_material_index) >= materials.size())
        return edges;
    const CrossSectionData& db = CrossSectionData::instance();
    const MaterialSpec& crystal =
        materials[static_cast<size_t>(crystal_material_index)];
    for (const MaterialComponent& c : crystal.composition) {
        if (const FluorescenceData* f = db.fluorescence(c.Z)) {
            const double e = f->k_edge_keV;
            if (e > e_min_keV * 1.02 && e < e_max_keV * 0.98) edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end(),
                            [](double a, double b) { return b - a < 0.01; }),
                edges.end());
    return edges;
}

// ---------------------------------------------------------------------------
// Plain-text codec
//
//   CEELO_GEOMETRY 1
//   shape cylinder|box
//   dimensions <cm...>                       radius length | half_x half_y length
//   bullet_radius <cm>
//   crystal_material <index>
//   bore <radius> <depth> <rounded 0|1>      optional
//   dead_layer <front> <side> <back>         optional
//   layer <material> <front> <side> <z_start> <z_end>    zero or more, innermost first
//   collimator <material> <side> <z_start> <z_end>       optional
//   symmetry axial|quadrant
//   reference_point crystal_face|endcap_front
//   material <index> <density g/cm3> <n> <Z> <mass fraction> ... <name>
//                                            indices 0, 1, ... in order; the
//                                            name is the rest of the line
//   end
//
// Numbers are written with %.17g, so a round trip is exact.
// ---------------------------------------------------------------------------

namespace {

std::string number(double v) {
    char buffer[40];
    std::snprintf(buffer, sizeof buffer, "%.17g", v);
    return buffer;
}

[[noreturn]] void text_error(int line, const std::string& what) {
    throw std::runtime_error("GeometryDescriptor::from_text: line "
                             + std::to_string(line) + ": " + what);
}

}  // namespace

std::string GeometryDescriptor::to_text() const {
    std::ostringstream out;
    out << "CEELO_GEOMETRY 1\n";
    out << "shape " << (shape == DetectorShape::Box ? "box" : "cylinder") << '\n';
    out << "dimensions";
    for (double v : dimensions_cm) out << ' ' << number(v);
    out << '\n';
    out << "bullet_radius " << number(bullet_radius_cm) << '\n';
    out << "crystal_material " << crystal_material_index << '\n';
    if (bore)
        out << "bore " << number(bore->radius) << ' ' << number(bore->depth) << ' '
            << (bore->rounded_tip ? 1 : 0) << '\n';
    if (dead_layer)
        out << "dead_layer " << number(dead_layer->front) << ' '
            << number(dead_layer->side) << ' ' << number(dead_layer->back) << '\n';
    for (const LayerSpec& l : layers)
        out << "layer " << l.material_index << ' ' << number(l.front_thickness_cm) << ' '
            << number(l.side_thickness_cm) << ' ' << number(l.z_start_cm) << ' '
            << number(l.z_end_cm) << '\n';
    if (collimator)
        out << "collimator " << collimator->material_index << ' '
            << number(collimator->side_thickness_cm) << ' '
            << number(collimator->z_start_cm) << ' ' << number(collimator->z_end_cm) << '\n';
    out << "symmetry "
        << (symmetry == ResponseSymmetry::Quadrant ? "quadrant" : "axial") << '\n';
    out << "reference_point "
        << (reference_point == ReferencePoint::EndcapFront ? "endcap_front" : "crystal_face")
        << '\n';
    for (size_t i = 0; i < materials.size(); ++i) {
        const MaterialSpec& m = materials[i];
        out << "material " << i << ' ' << number(m.density_g_per_cm3) << ' '
            << m.composition.size();
        for (const MaterialComponent& c : m.composition)
            out << ' ' << static_cast<int>(c.Z) << ' ' << number(c.mass_fraction);
        out << ' ' << m.name << '\n';
    }
    out << "end\n";
    return out.str();
}

GeometryDescriptor GeometryDescriptor::from_text(const std::string& text) {
    GeometryDescriptor gd;
    std::istringstream in(text);
    std::string line;
    int number_of_line = 0;
    bool header = false, ended = false, have_shape = false, have_dimensions = false;
    bool have_crystal = false;
    while (std::getline(in, line)) {
        ++number_of_line;
        std::istringstream fields(line);
        std::string key;
        if (!(fields >> key) || key[0] == '#') continue;
        if (ended) text_error(number_of_line, "content after 'end'");
        if (!header) {
            int version = 0;
            if (key != "CEELO_GEOMETRY" || !(fields >> version) || version != 1)
                text_error(number_of_line, "expected 'CEELO_GEOMETRY 1'");
            header = true;
            continue;
        }
        // Explicit captures: MSVC fails to implicitly capture names used next to a
        // fold expression in a generic lambda (C2065 on number_of_line / key).
        auto need = [&fields, &number_of_line, &key](auto&... values) {
            if (!((fields >> values) && ...))
                text_error(number_of_line, "malformed '" + key + "'");
        };
        if (key == "shape") {
            std::string s;
            need(s);
            if (s == "cylinder") gd.shape = DetectorShape::Cylinder;
            else if (s == "box") gd.shape = DetectorShape::Box;
            else text_error(number_of_line, "unknown shape '" + s + "'");
            have_shape = true;
        } else if (key == "dimensions") {
            double v;
            while (fields >> v) gd.dimensions_cm.push_back(v);
            if (gd.dimensions_cm.empty()) text_error(number_of_line, "no dimensions");
            have_dimensions = true;
        } else if (key == "bullet_radius") {
            need(gd.bullet_radius_cm);
        } else if (key == "crystal_material") {
            need(gd.crystal_material_index);
            have_crystal = true;
        } else if (key == "bore") {
            BoreHoleConfig b{0.0, 0.0};
            int rounded = 0;
            need(b.radius, b.depth, rounded);
            b.rounded_tip = rounded != 0;
            gd.bore = b;
        } else if (key == "dead_layer") {
            DeadLayerConfig d{0.0, 0.0, 0.0};
            need(d.front, d.side, d.back);
            gd.dead_layer = d;
        } else if (key == "layer") {
            LayerSpec l;
            need(l.material_index, l.front_thickness_cm, l.side_thickness_cm,
                 l.z_start_cm, l.z_end_cm);
            gd.layers.push_back(l);
        } else if (key == "collimator") {
            CollimatorSpec c;
            need(c.material_index, c.side_thickness_cm, c.z_start_cm, c.z_end_cm);
            gd.collimator = c;
        } else if (key == "symmetry") {
            std::string s;
            need(s);
            if (s == "axial") gd.symmetry = ResponseSymmetry::Axial;
            else if (s == "quadrant") gd.symmetry = ResponseSymmetry::Quadrant;
            else text_error(number_of_line, "unknown symmetry '" + s + "'");
        } else if (key == "reference_point") {
            std::string s;
            need(s);
            if (s == "crystal_face") gd.reference_point = ReferencePoint::CrystalFace;
            else if (s == "endcap_front") gd.reference_point = ReferencePoint::EndcapFront;
            else text_error(number_of_line, "unknown reference point '" + s + "'");
        } else if (key == "material") {
            size_t index = 0, count = 0;
            MaterialSpec m;
            need(index, m.density_g_per_cm3, count);
            if (index != gd.materials.size())
                text_error(number_of_line, "materials must be numbered 0, 1, ... in order");
            for (size_t k = 0; k < count; ++k) {
                int z = 0;
                double fraction = 0.0;
                need(z, fraction);
                if (z < 1 || z > 255)
                    text_error(number_of_line, "atomic number out of range");
                m.composition.push_back({static_cast<uint8_t>(z), fraction});
            }
            std::getline(fields >> std::ws, m.name);
            gd.materials.push_back(m);
        } else if (key == "end") {
            ended = true;
        } else {
            text_error(number_of_line, "unknown keyword '" + key + "'");
        }
    }
    if (!header) text_error(number_of_line, "missing 'CEELO_GEOMETRY 1'");
    if (!ended) text_error(number_of_line, "missing 'end'");
    if (!have_shape || !have_dimensions || !have_crystal)
        text_error(number_of_line, "'shape', 'dimensions' and 'crystal_material' are required");
    return gd;
}

}  // namespace ceelo
