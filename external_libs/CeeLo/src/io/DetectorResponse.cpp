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

/// @file DetectorResponse.cpp
/// @brief Storable MC-parameterized detector response (see DetectorResponse.h).

#include "io/DetectorResponse.h"

#include "cross_sections/CrossSectionData.h"

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"  // third-party header
#endif
// rapidxml is NOT vendored: the host application's copy is used (located by
// CEELO_RAPIDXML_DIR in CMake, e.g. InterSpec's SpecUtils/3rdparty) for ODR safety.
#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ceelo {

// v2 added the crystal fillet radius and the rounded bore tip. A pre-v2 reader
// would not merely lose them, it would MISREAD the file -- silently tracing a
// sharp-edged, flat-bored crystal (worth up to ~30% efficiency at the worst
// near-field angles) and computing a different content_hash than the writer.
//
// So the version is written PER FILE, not per build: serialize_xml() emits "1"
// whenever every v2 field is at its default, and "2" only when a fillet or a
// rounded tip is actually present. Old builds therefore keep reading every file
// they can read correctly, and hard-fail (via the range test in from_xml)
// on exactly the ones they would get wrong. Existing responses keep their bytes
// and their hash. Follow this pattern for the next additive field.
const int DetectorResponse::sm_xmlSerializationVersion = 2;

namespace {

constexpr double kLnClampEps = 1e-30;

// ---- small helpers ---------------------------------------------------------

std::string fmt_double(double v) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return buf;
}

std::string join_doubles(const std::vector<double>& v) {
    std::string out;
    out.reserve(v.size() * 20);
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out += ' ';
        out += fmt_double(v[i]);
    }
    return out;
}

std::vector<double> parse_doubles(const char* txt) {
    std::vector<double> out;
    if (!txt) return out;
    const char* p = txt;
    char* end = nullptr;
    while (*p) {
        const double v = std::strtod(p, &end);
        if (end == p) break;
        out.push_back(v);
        p = end;
    }
    return out;
}

double clamp(double v, double lo, double hi) {
    return std::min(std::max(v, lo), hi);
}

// Linear interpolation helper on ascending grid; clamps outside.
// Returns interpolation weight/bracket via out params.
void bracket(const std::vector<double>& x, double q, size_t& i0, size_t& i1,
             double& w) {
    if (x.size() < 2 || q <= x.front()) { i0 = i1 = 0; w = 0.0; return; }
    if (q >= x.back()) { i0 = i1 = x.size() - 1; w = 0.0; return; }
    i1 = static_cast<size_t>(std::upper_bound(x.begin(), x.end(), q) - x.begin());
    i0 = i1 - 1;
    w = (q - x[i0]) / (x[i1] - x[i0]);
}

// Fold an azimuth (deg) into [0, span] using the reflective symmetry of the
// box quadrant/octant table (period 2*span).
double fold_phi_deg(double phi_deg, double span_deg) {
    if (span_deg <= 0.0) return 0.0;
    double p = std::fmod(std::fabs(phi_deg), 2.0 * span_deg);
    if (p > span_deg) p = 2.0 * span_deg - p;
    return p;
}

uint64_t fnv1a_64(const std::string& s) {
    uint64_t h = 1469598103934665603ULL;
    for (const char c : s) {
        h ^= static_cast<unsigned char>(c);
        h *= 1099511628211ULL;
    }
    return h;
}

using XmlNode = rapidxml::xml_node<char>;
using XmlDoc = rapidxml::xml_document<char>;

XmlNode* append_node(XmlDoc& doc, XmlNode* parent, const char* name) {
    XmlNode* n = doc.allocate_node(rapidxml::node_element,
                                   doc.allocate_string(name));
    parent->append_node(n);
    return n;
}

void append_attrib(XmlDoc& doc, XmlNode* node, const char* name,
                   const std::string& value) {
    node->append_attribute(doc.allocate_attribute(
        doc.allocate_string(name), doc.allocate_string(value.c_str())));
}

XmlNode* append_value_node(XmlDoc& doc, XmlNode* parent, const char* name,
                           const std::string& value) {
    XmlNode* n = append_node(doc, parent, name);
    n->value(doc.allocate_string(value.c_str()));
    return n;
}

const char* attrib_value(const XmlNode* node, const char* name) {
    const rapidxml::xml_attribute<char>* a = node->first_attribute(name);
    return a ? a->value() : nullptr;
}

double attrib_double(const XmlNode* node, const char* name, double fallback) {
    const char* v = attrib_value(node, name);
    return v ? std::strtod(v, nullptr) : fallback;
}

bool attrib_bool(const XmlNode* node, const char* name, bool fallback) {
    const char* v = attrib_value(node, name);
    if (!v) return fallback;
    return (v[0] == '1' || v[0] == 't' || v[0] == 'T'
            || v[0] == 'y' || v[0] == 'Y');
}

const char* child_value(const XmlNode* node, const char* name) {
    const XmlNode* c = node->first_node(name);
    return c ? c->value() : nullptr;
}

}  // namespace

// ---------------------------------------------------------------------------
// enum <-> string
// ---------------------------------------------------------------------------

const char* to_string(ResponseFlag f) {
    switch (f) {
        case ResponseFlag::Ok:                 return "ok";
        case ResponseFlag::OutOfRangeClamped:  return "out_of_range_clamped";
        case ResponseFlag::NearFieldUnmodeled: return "near_field_unmodeled";
        case ResponseFlag::Shadowed:           return "shadowed";
        case ResponseFlag::NeedsMc:            return "needs_mc";
    }
    return "?";
}

const char* to_string(ResponseProfile p) {
    switch (p) {
        case ResponseProfile::FarField: return "far-field";
        case ResponseProfile::General:  return "general";
        case ResponseProfile::Contact:  return "contact";
    }
    return "?";
}

namespace {
ResponseProfile profile_from_string(const char* s) {
    if (s) {
        if (std::strcmp(s, "far-field") == 0) return ResponseProfile::FarField;
        if (std::strcmp(s, "contact") == 0) return ResponseProfile::Contact;
    }
    return ResponseProfile::General;
}
}  // namespace

// ---------------------------------------------------------------------------
// MaterialSpec / GeometryDescriptor
// ---------------------------------------------------------------------------

MaterialSpec MaterialSpec::from(const Material& m) {
    MaterialSpec s;
    s.name = m.name();
    s.density_g_per_cm3 = m.density();
    s.composition = m.composition();
    return s;
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
    g.set_detector(shape, mat_at(crystal_material_index), dimensions_cm);
    // set_detector() clears the fillet/bore/dead layer, so declare them after
    // it; fillet first, so bore_fits() sees the final crystal profile.
    if (bullet_radius_cm > 0.0) g.set_bullet_radius(bullet_radius_cm);
    if (bore) g.set_bore_hole(bore->radius, bore->depth, bore->rounded_tip);
    if (dead_layer)
        g.set_dead_layer(dead_layer->front, dead_layer->side, dead_layer->back);
    for (const LayerSpec& l : layers)
        g.add_attenuator(mat_at(l.material_index), l.front_thickness_cm,
                         l.side_thickness_cm, l.z_start_cm, l.z_end_cm);
    if (collimator)
        g.add_collimator(mat_at(collimator->material_index),
                         collimator->side_thickness_cm, collimator->z_start_cm,
                         collimator->z_end_cm);
    return g;
}

const char* to_string(ProductionMethod m) {
    switch (m) {
        case ProductionMethod::FullMc:          return "full_mc";
        case ProductionMethod::QuickMcTransfer: return "quick_mc_transfer";
        case ProductionMethod::CurveTransfer:   return "curve_transfer";
    }
    return "full_mc";
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

    // --- crystal itself: set_detector() indexes dimensions_cm unchecked ---
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
// MuTable
// ---------------------------------------------------------------------------

MacroscopicXS MuTable::eval(double e_keV) const {
    MacroscopicXS xs{0.0, 0.0, 0.0, 0.0};
    const size_t n = energy_keV.size();
    if (n == 0) return xs;
    size_t i0, i1;
    double w;
    bracket(energy_keV, e_keV, i0, i1, w);

    auto interp = [&](const std::vector<double>& mu) -> double {
        const double v0 = mu[i0], v1 = mu[i1];
        if (i0 == i1) return v0;
        if (v0 > 0.0 && v1 > 0.0) {
            // log-log linear (grid is log-spaced; w is linear in E -- use
            // ln-E weight for correctness).
            const double lw = (std::log(e_keV) - std::log(energy_keV[i0])) /
                              (std::log(energy_keV[i1]) - std::log(energy_keV[i0]));
            return std::exp((1.0 - lw) * std::log(v0) + lw * std::log(v1));
        }
        return (1.0 - w) * v0 + w * v1;  // zero-crossing (e.g. PP onset)
    };
    xs.mu_pe = interp(mu_pe);
    xs.mu_cs = interp(mu_cs);
    xs.mu_rs = interp(mu_rs);
    xs.mu_pp = interp(mu_pp);
    return xs;
}

MuTable MuTable::sample(const Material& mat, int material_index,
                        double e_min_keV, double e_max_keV, int n_per_decade) {
    MuTable t;
    t.material_index = material_index;

    std::vector<double> grid;
    const double decades = std::log10(e_max_keV / e_min_keV);
    const int n = std::max(2, static_cast<int>(std::ceil(decades * n_per_decade)));
    for (int i = 0; i <= n; ++i)
        grid.push_back(e_min_keV * std::pow(e_max_keV / e_min_keV,
                                            static_cast<double>(i) / n));

    // Absorption-edge flank pairs so interpolation never bridges an edge.
    const CrossSectionData& db = CrossSectionData::instance();
    for (const MaterialComponent& c : mat.composition()) {
        std::vector<double> edges;
        if (const FluorescenceData* f = db.fluorescence(c.Z))
            edges.push_back(f->k_edge_keV);
        if (const LFluorescenceData* l = db.l_fluorescence(c.Z))
            edges.push_back(l->l3_edge_keV);
        for (const double e : edges) {
            if (e > e_min_keV * 1.001 && e < e_max_keV * 0.999) {
                grid.push_back(e * (1.0 - 4e-4));
                grid.push_back(e * (1.0 + 4e-4));
            }
        }
    }
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end(),
                           [](double a, double b) { return b - a < a * 1e-7; }),
               grid.end());

    t.energy_keV = grid;
    t.mu_pe.reserve(grid.size());
    for (const double e : grid) {
        const MacroscopicXS xs = mat.macroscopic_xs(e * 1e-3);
        t.mu_pe.push_back(xs.mu_pe);
        t.mu_cs.push_back(xs.mu_cs);
        t.mu_rs.push_back(xs.mu_rs);
        t.mu_pp.push_back(xs.mu_pp);
    }
    return t;
}

// ---------------------------------------------------------------------------
// EtaTable
// ---------------------------------------------------------------------------

namespace {
// Segment index of E among ascending edges: number of edges strictly below E.
size_t segment_of(double e_keV, const std::vector<double>& edges) {
    size_t s = 0;
    for (const double e : edges)
        if (e_keV > e) ++s;
    return s;
}
}  // namespace

double EtaTable::SegCurve::eval(double lnE, bool& clamped) const {
    if (segs.empty()) { clamped = true; return 0.0; }
    // Find the segment whose range contains lnE (ranges tile the node span).
    size_t best = 0;
    for (size_t i = 0; i < segs.size(); ++i) {
        if (lnE >= seg_lo[i] - 1e-12 && lnE <= seg_hi[i] + 1e-12) { best = i; break; }
        // otherwise remember nearest by range distance
        const double d = std::min(std::fabs(lnE - seg_lo[i]),
                                  std::fabs(lnE - seg_hi[i]));
        const double db = std::min(std::fabs(lnE - seg_lo[best]),
                                   std::fabs(lnE - seg_hi[best]));
        if (d < db) best = i;
    }
    const double q = clamp(lnE, seg_lo[best], seg_hi[best]);
    if (q != lnE) clamped = true;
    return segs[best](q);
}

void EtaTable::finalize() {
    curves_.clear();
    if (empty()) return;
    const size_t ne = energies_keV.size();
    const size_t nc = cos_thetas.size();
    const size_t np = phis_deg.empty() ? 1 : phis_deg.size();
    if (ln_eta.size() != ne * nc * np)
        throw std::runtime_error("EtaTable: ln_eta size mismatch");
    if (frac_sigma.size() != ln_eta.size())
        throw std::runtime_error("EtaTable: frac_sigma size mismatch");

    // Group energy nodes into K-edge segments (both flanks are nodes).
    std::vector<size_t> seg_id(ne);
    for (size_t e = 0; e < ne; ++e)
        seg_id[e] = segment_of(energies_keV[e], edges_keV);

    curves_.resize(nc * np);
    for (size_t c = 0; c < nc; ++c) {
        for (size_t p = 0; p < np; ++p) {
            SegCurve& curve = curves_[c * np + p];
            size_t e0 = 0;
            while (e0 < ne) {
                size_t e1 = e0;
                while (e1 + 1 < ne && seg_id[e1 + 1] == seg_id[e0]) ++e1;
                std::vector<double> xs, ys;
                for (size_t e = e0; e <= e1; ++e) {
                    xs.push_back(std::log(energies_keV[e]));
                    ys.push_back(ln_eta[index(e, c, p)]);
                }
                if (xs.size() == 1) {  // lone node: constant stub segment
                    xs.push_back(xs[0] + 1e-9);
                    ys.push_back(ys[0]);
                }
                curve.seg_lo.push_back(xs.front());
                curve.seg_hi.push_back(xs.back());
                curve.segs.emplace_back(std::move(xs), std::move(ys));
                e0 = e1 + 1;
            }
        }
    }
}

double EtaTable::eval_ln(double energy_keV_q, double cos_theta, double phi_deg,
                         bool& clamped) const {
    if (curves_.empty())
        throw std::runtime_error("EtaTable::eval_ln before finalize()");
    const size_t nc = cos_thetas.size();
    const size_t np = phis_deg.empty() ? 1 : phis_deg.size();
    const double lnE = std::log(std::max(energy_keV_q, kLnClampEps));

    const double ct = clamp(cos_theta, cos_thetas.front(), cos_thetas.back());
    if (ct != cos_theta) clamped = true;

    double pq = 0.0;
    if (np > 1) {
        pq = fold_phi_deg(phi_deg, phis_deg.back());
        pq = clamp(pq, phis_deg.front(), phis_deg.back());
    }

    // Per phi node: energy eval per ct node, then PCHIP across cos-theta.
    auto eval_at_phi = [&](size_t p) -> double {
        std::vector<double> y(nc);
        for (size_t c = 0; c < nc; ++c)
            y[c] = curves_[c * np + p].eval(lnE, clamped);
        if (nc == 1) return y[0];
        const Pchip ang(cos_thetas, y);
        return ang(ct);
    };

    if (np == 1) return eval_at_phi(0);

    size_t p0, p1;
    double w;
    bracket(phis_deg, pq, p0, p1, w);
    const double v0 = eval_at_phi(p0);
    if (p0 == p1) return v0;
    return (1.0 - w) * v0 + w * eval_at_phi(p1);
}

double EtaTable::node_frac_sigma(double energy_keV_q, double cos_theta,
                                 double phi_deg) const {
    if (empty()) return 0.0;
    const size_t np = phis_deg.empty() ? 1 : phis_deg.size();

    // Nearest phi node.
    size_t p = 0;
    if (np > 1) {
        const double pq = fold_phi_deg(phi_deg, phis_deg.back());
        double best = 1e300;
        for (size_t i = 0; i < np; ++i) {
            const double d = std::fabs(phis_deg[i] - pq);
            if (d < best) { best = d; p = i; }
        }
    }

    // Bilinear on sigma^2 in (lnE, ct).
    std::vector<double> lnEs(energies_keV.size());
    for (size_t i = 0; i < lnEs.size(); ++i) lnEs[i] = std::log(energies_keV[i]);
    size_t e0, e1, c0, c1;
    double we, wc;
    bracket(lnEs, std::log(std::max(energy_keV_q, kLnClampEps)), e0, e1, we);
    bracket(cos_thetas, cos_theta, c0, c1, wc);

    auto s2 = [&](size_t e, size_t c) {
        const double s = frac_sigma[index(e, c, p)];
        return s * s;
    };
    const double v = (1.0 - we) * ((1.0 - wc) * s2(e0, c0) + wc * s2(e0, c1)) +
                     we * ((1.0 - wc) * s2(e1, c0) + wc * s2(e1, c1));
    return std::sqrt(std::max(0.0, v));
}

// ---------------------------------------------------------------------------
// NearFieldModel
// ---------------------------------------------------------------------------

void NearFieldModel::finalize() {
    d_curves_.clear();
    if (empty()) return;
    const size_t ne = energies_keV.size();
    const size_t nc = cos_thetas.size();
    const size_t nd = dists_cm.size();
    if (nc < 2 || nd < 2)
        throw std::runtime_error(
            "NearFieldModel: needs >= 2 cos_theta and >= 2 distance nodes");
    if (ln_n.size() != ne * nc * nd)
        throw std::runtime_error("NearFieldModel: ln_n size mismatch");
    if (frac_sigma.size() != ln_n.size())
        throw std::runtime_error("NearFieldModel: frac_sigma size mismatch");

    std::vector<double> ln_d(nd);
    for (size_t i = 0; i < nd; ++i) ln_d[i] = std::log(dists_cm[i]);

    // One monotone ln(d) curve per (E, cos_theta) node.
    d_curves_.reserve(ne * nc);
    for (size_t e = 0; e < ne; ++e) {
        for (size_t c = 0; c < nc; ++c) {
            std::vector<double> y(nd);
            for (size_t d = 0; d < nd; ++d) y[d] = ln_n[index(e, c, d)];
            d_curves_.emplace_back(ln_d, y);
        }
    }
}

double NearFieldModel::ln_boost(double energy_keV_q, double cos_theta,
                                double d_cm) const {
    if (empty()) return 0.0;
    if (d_curves_.empty())
        throw std::runtime_error("NearFieldModel::ln_boost before finalize()");
    const size_t nc = cos_thetas.size();
    const double ln_d = std::log(std::max(d_cm, 1e-6));
    const double ct = clamp(cos_theta, cos_thetas.front(), cos_thetas.back());

    std::vector<double> lnEs(energies_keV.size());
    for (size_t i = 0; i < lnEs.size(); ++i) lnEs[i] = std::log(energies_keV[i]);
    size_t e0, e1;
    double w;
    bracket(lnEs, std::log(std::max(energy_keV_q, kLnClampEps)), e0, e1, w);

    // At one E-node: eval each cos_theta node's ln(d) curve, then PCHIP in
    // cos_theta (which clamps ct outside the grid).
    auto at = [&](size_t e) -> double {
        std::vector<double> y(nc);
        for (size_t c = 0; c < nc; ++c) y[c] = d_curves_[e * nc + c](ln_d);
        const Pchip ang(cos_thetas, y);
        return ang(ct);
    };
    const double v0 = at(e0);
    if (e0 == e1) return v0;
    return (1.0 - w) * v0 + w * at(e1);
}

double NearFieldModel::node_frac_sigma(double energy_keV_q, double cos_theta,
                                       double d_cm) const {
    if (empty()) return 0.0;

    // Trilinear on sigma^2 in (lnE, cos_theta, ln d).
    std::vector<double> lnEs(energies_keV.size());
    for (size_t i = 0; i < lnEs.size(); ++i) lnEs[i] = std::log(energies_keV[i]);
    std::vector<double> lnDs(dists_cm.size());
    for (size_t i = 0; i < lnDs.size(); ++i) lnDs[i] = std::log(dists_cm[i]);

    size_t e0, e1, c0, c1, d0, d1;
    double we, wc, wd;
    bracket(lnEs, std::log(std::max(energy_keV_q, kLnClampEps)), e0, e1, we);
    bracket(cos_thetas, cos_theta, c0, c1, wc);
    bracket(lnDs, std::log(std::max(d_cm, 1e-6)), d0, d1, wd);

    auto s2 = [&](size_t e, size_t c, size_t d) {
        const double s = frac_sigma[index(e, c, d)];
        return s * s;
    };
    auto at_e = [&](size_t e) {
        const double v0 =
            (1.0 - wc) * s2(e, c0, d0) + wc * s2(e, c1, d0);
        const double v1 =
            (1.0 - wc) * s2(e, c0, d1) + wc * s2(e, c1, d1);
        return (1.0 - wd) * v0 + wd * v1;
    };
    const double v = (1.0 - we) * at_e(e0) + we * at_e(e1);
    return std::sqrt(std::max(0.0, v));
}

double NearFieldModel::breakpoint_d_cm(double energy_keV_q,
                                       double cos_theta) const {
    if (break_cos_thetas.empty() || break_d_cm.empty()) return 0.0;
    const size_t nc = break_cos_thetas.size();
    if (break_d_cm.size() != energies_keV.size() * nc) return 0.0;

    std::vector<double> lnEs(energies_keV.size());
    for (size_t i = 0; i < lnEs.size(); ++i) lnEs[i] = std::log(energies_keV[i]);
    size_t e0, e1, c0, c1;
    double we, wc;
    bracket(lnEs, std::log(std::max(energy_keV_q, kLnClampEps)), e0, e1, we);
    bracket(break_cos_thetas, cos_theta, c0, c1, wc);

    auto at = [&](size_t e, size_t c) { return break_d_cm[e * nc + c]; };
    return (1.0 - we) * ((1.0 - wc) * at(e0, c0) + wc * at(e0, c1)) +
           we * ((1.0 - wc) * at(e1, c0) + wc * at(e1, c1));
}

// ---------------------------------------------------------------------------
// TotEffPayload
// ---------------------------------------------------------------------------

void TotEffPayload::finalize() {
    if (tier == TotEffTier::BCurve) {
        if (b_energies_keV.size() < 2)
            throw std::runtime_error("TotEffPayload: BCurve needs >= 2 nodes");
        std::vector<double> x(b_energies_keV.size());
        for (size_t i = 0; i < x.size(); ++i) x[i] = std::log(b_energies_keV[i]);
        b_curve_ = Pchip(std::move(x), ln_b);
    }
    if (tier == TotEffTier::EtaTotTable)
        eta_tot.finalize();
}

double TotEffPayload::ln_b_at(double energy_keV) const {
    if (!b_curve_.valid()) return 0.0;
    return b_curve_(std::log(std::max(energy_keV, kLnClampEps)));
}

// ---------------------------------------------------------------------------
// Grounding
// ---------------------------------------------------------------------------

double SigmaTransferModel::eval(double d_over_a, double cos_theta,
                                double energy_keV) const {
    const double s2 = std::max(0.0, 1.0 - cos_theta * cos_theta);
    // Low-E weight: 0 at/above mid_e_ref, 1 at/below low_e_ref (ln ramp).
    double w = 0.0;
    if (energy_keV < mid_e_ref_keV) {
        w = (std::log(mid_e_ref_keV) - std::log(std::max(energy_keV, 1.0))) /
            (std::log(mid_e_ref_keV) - std::log(low_e_ref_keV));
        w = clamp(w, 0.0, 1.0);
    }
    const double off = s2 * (offaxis_mid + offaxis_low_e * w * w);
    double near = 0.0;
    if (d_over_a < near_gate_a) {
        const double t = clamp((near_gate_a - d_over_a) / (near_gate_a - 1.0),
                               0.0, 1.0);
        near = near_contact * t;
    }
    return std::sqrt(far_onaxis * far_onaxis + off * off + near * near);
}

double GroundingBlock::eval_ln_k(double energy_keV, bool& clamped) const {
    if (empty()) return 0.0;
    const double lnE = std::log(std::max(energy_keV, kLnClampEps));
    if (lnE < knot_ln_energies.front() || lnE > knot_ln_energies.back())
        clamped = true;
    const std::vector<double> B = hat_basis(lnE, knot_ln_energies);
    double v = 0.0;
    for (size_t i = 0; i < B.size(); ++i) v += B[i] * ln_k[i];
    return v;
}

double GroundingBlock::cov_ln_k(double e1_keV, double e2_keV) const {
    if (empty() || cov.empty()) return 0.0;
    const size_t n = knot_ln_energies.size();
    if (cov.size() != n * n) return 0.0;
    const std::vector<double> B1 =
        hat_basis(std::log(std::max(e1_keV, kLnClampEps)), knot_ln_energies);
    const std::vector<double> B2 =
        hat_basis(std::log(std::max(e2_keV, kLnClampEps)), knot_ln_energies);
    double v = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (B1[i] == 0.0) continue;
        for (size_t j = 0; j < n; ++j)
            v += B1[i] * cov[i * n + j] * B2[j];
    }
    return v;
}

double GroundingBlock::var_ln_k(double energy_keV) const {
    return cov_ln_k(energy_keV, energy_keV);
}

// ---------------------------------------------------------------------------
// DetectorResponse -- assembly / geometry
// ---------------------------------------------------------------------------

void DetectorResponse::finalize() {
    geometry_ = descriptor.build_geometry(owned_materials_);
    geometry_built_ = true;

    // Map every material the ray tracer can hand back to its mu table.
    mat_to_mu_.clear();
    for (size_t mi = 0; mi < owned_materials_.size(); ++mi) {
        size_t table = SIZE_MAX;
        for (size_t t = 0; t < mu_tables.size(); ++t)
            if (mu_tables[t].material_index == static_cast<int>(mi))
                table = t;
        if (table == SIZE_MAX)
            throw std::runtime_error(
                "DetectorResponse: missing mu table for material '" +
                descriptor.materials[mi].name + "'");
        mat_to_mu_.push_back({owned_materials_[mi].get(), table});
    }

    eta_fep.finalize();
    tot_eff.finalize();
    near_field.finalize();  // tolerates empty (far-field profile)
}

Eigen::Vector3d DetectorResponse::query_position(double theta_rad,
                                                 double phi_rad,
                                                 double dist_cm) const {
    Eigen::Vector3d pos =
        source_position(dist_cm, std::cos(theta_rad), phi_rad);
    if (descriptor.reference_point == ReferencePoint::EndcapFront)
        pos.z() -= descriptor.endcap_front_offset_cm();
    return pos;
}

ApertureQuadrature DetectorResponse::make_quadrature(
    const Eigen::Vector3d& src_cm) const {
    return make_aperture_quadrature(geometry_, src_cm,
                                    provenance.kernel_n_rays);
}

void DetectorResponse::kernel_ray_weights_impl(
    double energy_keV, const ApertureQuadrature& q, MuChoice mu, double recap,
    std::vector<double>& w_out, std::vector<Eigen::Vector3d>& dirs_out) const {
    // Same recurrence as ApertureQuadrature::interaction_omega, but with the
    // STORED generation-time mu tables (self-containment guarantee).
    struct MuT { double tot, nors, cs; };
    std::vector<std::pair<const Material*, MuT>> cache;
    auto mu_of = [&](const Material* m) -> MuT {
        for (const auto& e : cache)
            if (e.first == m) return e.second;
        size_t table = SIZE_MAX;
        for (const auto& mm : mat_to_mu_)
            if (mm.first == m) { table = mm.second; break; }
        if (table == SIZE_MAX)
            throw std::runtime_error("DetectorResponse::kernel_K: unknown material");
        const MacroscopicXS xs = mu_tables[table].eval(energy_keV);
        const MuT v{xs.mu_total(), xs.mu_pe + xs.mu_cs + xs.mu_pp, xs.mu_cs};
        cache.push_back({m, v});
        return v;
    };

    w_out.clear();
    dirs_out.clear();
    w_out.reserve(q.rays.size());
    dirs_out.reserve(q.rays.size());

    for (const KernelRay& r : q.rays) {
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;
        double p_int = 0.0;
        for (const RaySegment& s : r.segs) {
            const MuT m = mu_of(s.material);
            if (s.is_scoring) {
                const double mu_star = (mu == MuChoice::Total) ? m.tot : m.nors;
                p_int += std::exp(-tau_before) *
                         (1.0 - std::exp(-mu_star * s.length));
                tau_before += m.tot * s.length;
            } else {
                tau_before += (m.tot - recap * m.cs) * s.length;
            }
        }
        if (p_int <= 0.0) continue;
        w_out.push_back(r.omega_w * p_int);
        dirs_out.push_back(r.dir.cast<double>().eval());
    }
}

double DetectorResponse::kernel_K(
    double energy_keV, const ApertureQuadrature& q, MuChoice mu,
    const std::function<double(const Eigen::Vector3d&)>* t_src) const {
    // Dead-layer / endcap scatter-in credit (total kernel only); see
    // ApertureQuadrature::interaction_omega. Uses the stored recapture so the
    // anchor and target kernels are folded with the same coefficient.
    const double recap = (mu == MuChoice::NoRayleigh)
        ? std::max(0.0, std::min(1.0, scatter_in_recapture)) : 0.0;

    // Deliberately the inline loop rather than a call to kernel_ray_weights_impl: this runs per
    //  volumetric element per fit iteration, and routing it through the vector-filling form would
    //  add two heap allocations and a per-ray Vector3d materialization even when t_src is null.
    //  The risk that the two bodies drift is covered by
    //  tests/test_efficiency_transfer.cpp::ray_weight_decomposition_matches_full_query, which pins
    //  prefactor * sum(weights) == the full query exactly.
    struct MuT { double tot, nors, cs; };
    std::vector<std::pair<const Material*, MuT>> cache;
    auto mu_of = [&](const Material* m) -> MuT {
        for (const auto& e : cache)
            if (e.first == m) return e.second;
        size_t table = SIZE_MAX;
        for (const auto& mm : mat_to_mu_)
            if (mm.first == m) { table = mm.second; break; }
        if (table == SIZE_MAX)
            throw std::runtime_error("DetectorResponse::kernel_K: unknown material");
        const MacroscopicXS xs = mu_tables[table].eval(energy_keV);
        const MuT v{xs.mu_total(), xs.mu_pe + xs.mu_cs + xs.mu_pp, xs.mu_cs};
        cache.push_back({m, v});
        return v;
    };

    double total = 0.0;
    for (const KernelRay& r : q.rays) {
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;
        double p_int = 0.0;
        for (const RaySegment& s : r.segs) {
            const MuT m = mu_of(s.material);
            if (s.is_scoring) {
                const double mu_star = (mu == MuChoice::Total) ? m.tot : m.nors;
                p_int += std::exp(-tau_before) *
                         (1.0 - std::exp(-mu_star * s.length));
                tau_before += m.tot * s.length;
            } else {
                tau_before += (m.tot - recap * m.cs) * s.length;
            }
        }
        if (p_int <= 0.0) continue;
        double w = r.omega_w * p_int;
        if (t_src) w *= (*t_src)(r.dir.cast<double>().eval());
        total += w;
    }
    return total;
}

void DetectorResponse::fep_ray_weights(
    double energy_keV, const ApertureQuadrature& q,
    std::vector<double>& w_out, std::vector<Eigen::Vector3d>& dirs_out) const {
    // FEP always uses the total mu on the scoring segments, and takes no scatter-in credit: a
    // degraded Compton photon cannot land in the full-energy peak.
    kernel_ray_weights_impl(energy_keV, q, MuChoice::Total, 0.0, w_out, dirs_out);
}

void DetectorResponse::total_ray_weights(
    double energy_keV, const ApertureQuadrature& q,
    std::vector<double>& w_out, std::vector<Eigen::Vector3d>& dirs_out) const {
    // Mirrors eps_total_impl's tier dispatch - the EtaTotTable tier folds a measured eta_tot over
    // the FULL-mu kernel, the other tiers use the Rayleigh-free one.
    const MuChoice mu = (tot_eff.tier == TotEffTier::EtaTotTable) ? MuChoice::Total
                                                                  : MuChoice::NoRayleigh;
    const double recap = (mu == MuChoice::NoRayleigh)
        ? std::max(0.0, std::min(1.0, scatter_in_recapture)) : 0.0;
    kernel_ray_weights_impl(energy_keV, q, mu, recap, w_out, dirs_out);
}

// ---------------------------------------------------------------------------
// DetectorResponse -- evaluation
// ---------------------------------------------------------------------------

struct DetectorResponse::EvalCommon {
    double d_cm = 0.0;          // from crystal-face origin
    double cos_theta = 1.0;
    double phi_deg = 0.0;
    double a_cm = 0.0;          // transverse half-extent
    bool near_regime = false;
    ResponseFlag flag = ResponseFlag::Ok;
    double extra_sigma2 = 0.0;  // shadow / behind-plane inflation
};

namespace {
void raise_flag(ResponseFlag& current, ResponseFlag f) {
    if (static_cast<int>(f) > static_cast<int>(current)) current = f;
}
}  // namespace

DetectorResponse::EvalCommon DetectorResponse::common_eval(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q) const {
    EvalCommon ec;
    ec.d_cm = src_cm.norm();
    ec.cos_theta = (ec.d_cm > 0.0) ? (-src_cm.z() / ec.d_cm) : 1.0;
    ec.phi_deg = std::atan2(src_cm.y(), src_cm.x()) * 180.0 / 3.14159265358979323846;
    ec.a_cm = descriptor.transverse_half_extent();
    ec.near_regime = ec.d_cm < floors.near_regime_a * ec.a_cm;

    // Energy validity from provenance (falls back to the table range).
    const double e_lo = provenance.valid_e_min_keV > 0.0
        ? provenance.valid_e_min_keV
        : (eta_fep.empty() ? 0.0 : eta_fep.energies_keV.front());
    const double e_hi = provenance.valid_e_max_keV > 0.0
        ? provenance.valid_e_max_keV
        : (eta_fep.empty() ? 0.0 : eta_fep.energies_keV.back());
    if (e_hi > 0.0 && (energy_keV < e_lo || energy_keV > e_hi))
        raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);

    // Behind the face plane: unmodeled (theta > 90; Marinelli-class needs
    // dedicated nodes -- deferred). Value is a clamped guess.
    if (ec.cos_theta < 0.0) {
        raise_flag(ec.flag, ResponseFlag::NeedsMc);
        ec.extra_sigma2 += 0.30 * 0.30;
    }

    // EFFTRAN transfer envelope: inflate sigma where the (angle-flat) eta of a
    // transfer response cannot carry the true eta(E,theta) / near-field
    // residual. Applied UNCONDITIONALLY (independent of grounding), and to both
    // FEP and total via this shared path.
    if (model_transfer) {
        const double st = model_transfer->eval(
            ec.a_cm > 0.0 ? ec.d_cm / ec.a_cm : 1e6, ec.cos_theta, energy_keV);
        ec.extra_sigma2 += st * st;
    }

    // Collimator shadow gate (spec sec 4.5): s = transmitted/geometric.
    if (descriptor.collimator && q.omega_frac_active > 0.0) {
        const double s = kernel_transmitted(energy_keV, q) / q.omega_frac_active;
        if (s < 0.05) {
            raise_flag(ec.flag, ResponseFlag::NeedsMc);
            ec.extra_sigma2 += 1.0;  // sigma ~ 100%
        } else if (s < 0.3) {
            raise_flag(ec.flag, ResponseFlag::Shadowed);
            const double sh = 0.5 * (0.3 - s) / 0.3;
            ec.extra_sigma2 += sh * sh;
        }
    }
    return ec;
}

double DetectorResponse::kernel_transmitted(double energy_keV,
                                            const ApertureQuadrature& q) const {
    std::vector<std::pair<const Material*, double>> cache;
    auto mu_of = [&](const Material* m) -> double {
        for (const auto& e : cache)
            if (e.first == m) return e.second;
        size_t table = SIZE_MAX;
        for (const auto& mm : mat_to_mu_)
            if (mm.first == m) { table = mm.second; break; }
        if (table == SIZE_MAX)
            throw std::runtime_error("DetectorResponse: unknown material");
        const double v = mu_tables[table].eval(energy_keV).mu_total();
        cache.push_back({m, v});
        return v;
    };
    double total = 0.0;
    for (const KernelRay& r : q.rays) {
        if (r.active_len <= 0.0f) continue;
        double tau = 0.0;
        for (const RaySegment& s : r.segs) {
            if (s.is_scoring) break;
            tau += mu_of(s.material) * s.length;
        }
        total += r.omega_w * std::exp(-tau);
    }
    return total;
}

EffResult DetectorResponse::fep_prefactor(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q) const {
    // Deliberately the same body as eps_fep_impl with the `* K` dropped, so a host that assembles
    // K itself still gets the near-field gate, the grounding, the sigma budget and the flag.
    EvalCommon ec = common_eval(energy_keV, src_cm, q);

    bool clamped = false;
    const double ln_eta =
        eta_fep.eval_ln(energy_keV, ec.cos_theta, ec.phi_deg, clamped);
    if (clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);

    double ln_N = 0.0;
    const double d_break = near_field.breakpoint_d_cm(energy_keV, ec.cos_theta);
    const double d_gate = std::max(d_break, provenance.min_distance_cm);
    if (ec.d_cm < d_gate) {
        if (!near_field.empty()) {
            ln_N = near_field.ln_boost(energy_keV, ec.cos_theta, ec.d_cm);
            const double nf_sig =
                near_field.node_frac_sigma(energy_keV, ec.cos_theta, ec.d_cm);
            ec.extra_sigma2 += nf_sig * nf_sig;
        } else {
            raise_flag(ec.flag, ResponseFlag::NearFieldUnmodeled);
            ec.extra_sigma2 += 0.05 * 0.05;
        }
    }

    double ln_k = 0.0, ground_var = 0.0;
    if (!grounding.empty()) {
        bool k_clamped = false;
        ln_k = grounding.eval_ln_k(energy_keV, k_clamped);
        if (k_clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
        const double st = grounding.transfer.eval(
            ec.a_cm > 0.0 ? ec.d_cm / ec.a_cm : 1e6, ec.cos_theta, energy_keV);
        ground_var = grounding.var_ln_k(energy_keV) + st * st;
    }

    EffResult res;
    res.value = std::exp(ln_eta + ln_N + ln_k);
    const double node_sig =
        eta_fep.node_frac_sigma(energy_keV, ec.cos_theta, ec.phi_deg);
    const double floor = ec.near_regime ? floors.fep_near : floors.fep_far;
    const double frac2 = node_sig * node_sig + floor * floor + ground_var +
                         ec.extra_sigma2;
    res.sigma = res.value * std::sqrt(frac2);
    res.flag = ec.flag;
    return res;
}

EffResult DetectorResponse::total_prefactor(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q) const {
    // Mirrors eps_total_impl's tier dispatch with the kernel factored out.  The build-up seam is
    // NOT applied here: it needs a ShieldContext the host supplies, and folding it in silently
    // would double-count against a host that applies its own scatter augment.
    EvalCommon ec = common_eval(energy_keV, src_cm, q);

    double value = 1.0, node_sig = 0.0;
    switch (tot_eff.tier) {
        case TotEffTier::KernelExact:
            break;                                   //bare kernel; multiplier is 1
        case TotEffTier::BCurve:
            value = std::exp(tot_eff.ln_b_at(energy_keV));
            break;
        case TotEffTier::EtaTotTable: {
            bool clamped = false;
            const double ln_eta = tot_eff.eta_tot.eval_ln(
                energy_keV, ec.cos_theta, ec.phi_deg, clamped);
            if (clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
            double ln_k = 0.0;
            if (!grounding.empty()) {
                bool k_clamped = false;
                ln_k = grounding.eval_ln_k(energy_keV, k_clamped);
                if (k_clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
            }
            value = std::exp(ln_eta + ln_k);
            node_sig = tot_eff.eta_tot.node_frac_sigma(energy_keV, ec.cos_theta,
                                                       ec.phi_deg);
            break;
        }
    }

    EffResult res;
    res.value = value;
    const double floor = ec.near_regime ? floors.tot_near : floors.tot_far;
    const double frac2 = node_sig * node_sig + floor * floor + ec.extra_sigma2;
    res.sigma = res.value * std::sqrt(frac2);
    res.flag = ec.flag;
    return res;
}

EffResult DetectorResponse::eps_fep_impl(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q,
    const std::function<double(const Eigen::Vector3d&)>* t_src) const {
    EvalCommon ec = common_eval(energy_keV, src_cm, q);

    bool clamped = false;
    const double ln_eta =
        eta_fep.eval_ln(energy_keV, ec.cos_theta, ec.phi_deg, clamped);
    if (clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);

    const double K = kernel_K(energy_keV, q, MuChoice::Total, t_src);

    // Near field: below the measured breakpoint apply N (general/contact
    // profiles) or flag as unmodeled (far-field profile).
    double ln_N = 0.0;
    const double d_break = near_field.breakpoint_d_cm(energy_keV, ec.cos_theta);
    const double d_gate = std::max(d_break, provenance.min_distance_cm);
    if (ec.d_cm < d_gate) {
        if (!near_field.empty()) {
            ln_N = near_field.ln_boost(energy_keV, ec.cos_theta, ec.d_cm);
            const double nf_sig =
                near_field.node_frac_sigma(energy_keV, ec.cos_theta, ec.d_cm);
            ec.extra_sigma2 += nf_sig * nf_sig;
        } else {
            raise_flag(ec.flag, ResponseFlag::NearFieldUnmodeled);
            ec.extra_sigma2 += 0.05 * 0.05;  // kernel-only near error (S1)
        }
    }

    // Grounding k(E) + transfer inflation.
    double ln_k = 0.0, ground_var = 0.0;
    if (!grounding.empty()) {
        bool k_clamped = false;
        ln_k = grounding.eval_ln_k(energy_keV, k_clamped);
        if (k_clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
        const double st = grounding.transfer.eval(
            ec.a_cm > 0.0 ? ec.d_cm / ec.a_cm : 1e6, ec.cos_theta, energy_keV);
        ground_var = grounding.var_ln_k(energy_keV) + st * st;
    }

    EffResult res;
    res.value = std::exp(ln_eta + ln_N + ln_k) * K;
    const double node_sig =
        eta_fep.node_frac_sigma(energy_keV, ec.cos_theta, ec.phi_deg);
    const double floor = ec.near_regime ? floors.fep_near : floors.fep_far;
    const double frac2 = node_sig * node_sig + floor * floor + ground_var +
                         ec.extra_sigma2;
    res.sigma = res.value * std::sqrt(frac2);
    res.flag = ec.flag;
    return res;
}

EffResult DetectorResponse::eps_total_impl(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q,
    const std::function<double(const Eigen::Vector3d&)>* t_src,
    const ShieldContext* sc) const {
    EvalCommon ec = common_eval(energy_keV, src_cm, q);

    double value = 0.0;
    double node_sig = 0.0;
    switch (tot_eff.tier) {
        case TotEffTier::KernelExact:
            value = kernel_K(energy_keV, q, MuChoice::NoRayleigh, t_src);
            break;
        case TotEffTier::BCurve:
            value = std::exp(tot_eff.ln_b_at(energy_keV)) *
                    kernel_K(energy_keV, q, MuChoice::NoRayleigh, t_src);
            break;
        case TotEffTier::EtaTotTable: {
            bool clamped = false;
            const double ln_eta = tot_eff.eta_tot.eval_ln(
                energy_keV, ec.cos_theta, ec.phi_deg, clamped);
            if (clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
            double ln_k = 0.0;
            if (!grounding.empty()) {
                bool k_clamped = false;
                ln_k = grounding.eval_ln_k(energy_keV, k_clamped);
                if (k_clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
            }
            value = std::exp(ln_eta + ln_k) *
                    kernel_K(energy_keV, q, MuChoice::Total, t_src);
            node_sig = tot_eff.eta_tot.node_frac_sigma(energy_keV, ec.cos_theta,
                                                       ec.phi_deg);
            break;
        }
    }

    EffResult res;
    res.value = value;
    const double floor = ec.near_regime ? floors.tot_near : floors.tot_far;
    double frac2 = node_sig * node_sig + floor * floor + ec.extra_sigma2;

    // Build-up seam (Stage E3 A1): only when a ShieldContext is supplied AND a
    // model is installed.  Ratio-only (>= 1), scales eps_total up; the build-up
    // sigma floor is added in quadrature.  Null sc => this whole block is
    // skipped, so the value/sigma are byte-identical to the pre-seam path.
    if (sc && buildup_model) {
        const double B = std::max(1.0, buildup_model(energy_keV, *sc));
        res.value = value * B;
        frac2 += kBuildupSigmaFloor * kBuildupSigmaFloor;
    }
    res.sigma = res.value * std::sqrt(frac2);
    res.flag = ec.flag;
    return res;
}

// --- public wrappers ---

EffResult DetectorResponse::eps_fep(double energy_keV, double theta_rad,
                                    double phi_rad, double dist_cm) const {
    const Eigen::Vector3d pos = query_position(theta_rad, phi_rad, dist_cm);
    const ApertureQuadrature q = make_quadrature(pos);
    return eps_fep_impl(energy_keV, pos, q, nullptr);
}

EffResult DetectorResponse::eps_total(double energy_keV, double theta_rad,
                                      double phi_rad, double dist_cm) const {
    const Eigen::Vector3d pos = query_position(theta_rad, phi_rad, dist_cm);
    const ApertureQuadrature q = make_quadrature(pos);
    return eps_total_impl(energy_keV, pos, q, nullptr);
}

EffResult DetectorResponse::eps_fep_at(double energy_keV,
                                       const Eigen::Vector3d& src_cm) const {
    const ApertureQuadrature q = make_quadrature(src_cm);
    return eps_fep_impl(energy_keV, src_cm, q, nullptr);
}

EffResult DetectorResponse::eps_total_at(double energy_keV,
                                         const Eigen::Vector3d& src_cm,
                                         const ShieldContext* sc) const {
    const ApertureQuadrature q = make_quadrature(src_cm);
    return eps_total_impl(energy_keV, src_cm, q, nullptr, sc);
}

EffResult DetectorResponse::eps_fep_at(double energy_keV,
                                       const Eigen::Vector3d& src_cm,
                                       const ApertureQuadrature& q) const {
    return eps_fep_impl(energy_keV, src_cm, q, nullptr);
}

EffResult DetectorResponse::eps_total_at(double energy_keV,
                                         const Eigen::Vector3d& src_cm,
                                         const ApertureQuadrature& q,
                                         const ShieldContext* sc) const {
    return eps_total_impl(energy_keV, src_cm, q, nullptr, sc);
}

EffResult DetectorResponse::eps_fep_element(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q,
    const std::function<double(const Eigen::Vector3d&)>& t_src) const {
    return eps_fep_impl(energy_keV, src_cm, q, &t_src);
}

EffResult DetectorResponse::eps_total_element(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q,
    const std::function<double(const Eigen::Vector3d&)>& t_src,
    const ShieldContext* sc) const {
    return eps_total_impl(energy_keV, src_cm, q, &t_src, sc);
}

std::vector<double> DetectorResponse::frac_covariance(
    const std::vector<double>& energies_keV, double theta_rad,
    double dist_cm) const {
    const size_t n = energies_keV.size();
    std::vector<double> C(n * n, 0.0);
    const double ct = std::cos(theta_rad);
    const double a = descriptor.transverse_half_extent();
    const double d_over_a = a > 0.0 ? dist_cm / a : 1e6;
    const bool near = dist_cm < floors.near_regime_a * a;
    const double floor = near ? floors.fep_near : floors.fep_far;

    std::vector<double> node(n, 0.0), transfer(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        node[i] = eta_fep.node_frac_sigma(energies_keV[i], ct, 0.0);
        if (!grounding.empty())
            transfer[i] =
                grounding.transfer.eval(d_over_a, ct, energies_keV[i]);
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            double v = floor * floor + transfer[i] * transfer[j];
            if (!grounding.empty())
                v += grounding.cov_ln_k(energies_keV[i], energies_keV[j]);
            if (i == j) v += node[i] * node[i];
            C[i * n + j] = C[j * n + i] = v;
        }
    }
    return C;
}

// ---------------------------------------------------------------------------
// XML
// ---------------------------------------------------------------------------

namespace {

void eta_to_xml(XmlDoc& doc, XmlNode* parent, const char* name,
                const EtaTable& t) {
    XmlNode* n = append_node(doc, parent, name);
    append_value_node(doc, n, "Energies", join_doubles(t.energies_keV));
    append_value_node(doc, n, "CosThetas", join_doubles(t.cos_thetas));
    if (!t.phis_deg.empty())
        append_value_node(doc, n, "Phis", join_doubles(t.phis_deg));
    if (!t.edges_keV.empty())
        append_value_node(doc, n, "Edges", join_doubles(t.edges_keV));
    append_value_node(doc, n, "LnEta", join_doubles(t.ln_eta));
    append_value_node(doc, n, "FracSigma", join_doubles(t.frac_sigma));
}

void eta_from_xml(const XmlNode* n, EtaTable& t) {
    t.energies_keV = parse_doubles(child_value(n, "Energies"));
    t.cos_thetas = parse_doubles(child_value(n, "CosThetas"));
    t.phis_deg = parse_doubles(child_value(n, "Phis"));
    t.edges_keV = parse_doubles(child_value(n, "Edges"));
    t.ln_eta = parse_doubles(child_value(n, "LnEta"));
    t.frac_sigma = parse_doubles(child_value(n, "FracSigma"));
}

/// Writes the <Detector> element - the whole storable geometry, materials
/// included - under `parent`.  Shared by DetectorResponse::serialize_xml and
/// GeometryDescriptor::to_xml_string, so a geometry stored on its own and the
/// geometry inside a stored response cannot drift apart.  Factoring it out left
/// the emitted bytes unchanged, so every existing response keeps its
/// content_hash.
void write_detector_node(XmlDoc& doc, XmlNode* parent,
                         const GeometryDescriptor& descriptor) {
    XmlNode* d = append_node(doc, parent, "Detector");
    append_attrib(doc, d, "shape",
                  descriptor.shape == DetectorShape::Cylinder ? "cylinder"
                                                              : "box");
    append_attrib(doc, d, "crystalMaterial",
                  std::to_string(descriptor.crystal_material_index));
    append_attrib(doc, d, "referencePoint",
                  descriptor.reference_point == ReferencePoint::EndcapFront
                      ? "endcap_front" : "crystal_face");
    append_attrib(doc, d, "symmetry",
                  descriptor.symmetry == ResponseSymmetry::Quadrant
                      ? "quadrant" : "axial");
    append_value_node(doc, d, "Dimensions",
                      join_doubles(descriptor.dimensions_cm));
    // Both of the following are written only when non-default, so files
    // for sharp-edged, flat-bored crystals stay byte-identical (and keep
    // their content_hash) across this feature being added.
    if (descriptor.bullet_radius_cm > 0.0)
        append_attrib(doc, d, "bulletRadius",
                      fmt_double(descriptor.bullet_radius_cm));
    if (descriptor.bore) {
        XmlNode* b = append_node(doc, d, "Bore");
        append_attrib(doc, b, "radius", fmt_double(descriptor.bore->radius));
        append_attrib(doc, b, "depth", fmt_double(descriptor.bore->depth));
        if (descriptor.bore->rounded_tip)
            append_attrib(doc, b, "roundedTip", "1");
    }
    if (descriptor.dead_layer) {
        XmlNode* dl = append_node(doc, d, "DeadLayer");
        append_attrib(doc, dl, "front", fmt_double(descriptor.dead_layer->front));
        append_attrib(doc, dl, "side", fmt_double(descriptor.dead_layer->side));
        append_attrib(doc, dl, "back", fmt_double(descriptor.dead_layer->back));
    }
    for (const LayerSpec& l : descriptor.layers) {
        XmlNode* ln = append_node(doc, d, "Layer");
        append_attrib(doc, ln, "material", std::to_string(l.material_index));
        append_attrib(doc, ln, "frontThickness", fmt_double(l.front_thickness_cm));
        append_attrib(doc, ln, "sideThickness", fmt_double(l.side_thickness_cm));
        append_attrib(doc, ln, "zStart", fmt_double(l.z_start_cm));
        append_attrib(doc, ln, "zEnd", fmt_double(l.z_end_cm));
    }
    if (descriptor.collimator) {
        XmlNode* c = append_node(doc, d, "Collimator");
        append_attrib(doc, c, "material",
                      std::to_string(descriptor.collimator->material_index));
        append_attrib(doc, c, "sideThickness",
                      fmt_double(descriptor.collimator->side_thickness_cm));
        append_attrib(doc, c, "zStart", fmt_double(descriptor.collimator->z_start_cm));
        append_attrib(doc, c, "zEnd", fmt_double(descriptor.collimator->z_end_cm));
    }
    XmlNode* mats = append_node(doc, d, "Materials");
    for (size_t i = 0; i < descriptor.materials.size(); ++i) {
        const MaterialSpec& m = descriptor.materials[i];
        XmlNode* mn = append_node(doc, mats, "Material");
        append_attrib(doc, mn, "index", std::to_string(i));
        append_attrib(doc, mn, "name", m.name);
        append_attrib(doc, mn, "density", fmt_double(m.density_g_per_cm3));
        for (const MaterialComponent& c : m.composition) {
            XmlNode* el = append_node(doc, mn, "El");
            append_attrib(doc, el, "Z", std::to_string(int(c.Z)));
            append_attrib(doc, el, "frac", fmt_double(c.mass_fraction));
        }
    }

}

/// Reads a <Detector> element written by write_detector_node.
GeometryDescriptor read_detector_node(const XmlNode* d) {
    GeometryDescriptor gd;
    const char* shape = attrib_value(d, "shape");
    gd.shape = (shape && std::strcmp(shape, "box") == 0)
                   ? DetectorShape::Box : DetectorShape::Cylinder;
    gd.crystal_material_index =
        static_cast<int>(attrib_double(d, "crystalMaterial", -1));
    const char* rp = attrib_value(d, "referencePoint");
    gd.reference_point = (rp && std::strcmp(rp, "endcap_front") == 0)
                             ? ReferencePoint::EndcapFront
                             : ReferencePoint::CrystalFace;
    const char* sym = attrib_value(d, "symmetry");
    gd.symmetry = (sym && std::strcmp(sym, "quadrant") == 0)
                      ? ResponseSymmetry::Quadrant : ResponseSymmetry::Axial;
    gd.dimensions_cm = parse_doubles(child_value(d, "Dimensions"));
    // Absent in files written before the fillet/rounded-tip feature; the
    // defaults are what those crystals actually were.
    gd.bullet_radius_cm = attrib_double(d, "bulletRadius", 0.0);
    if (const XmlNode* b = d->first_node("Bore")) {
        BoreHoleConfig bc;
        bc.radius = attrib_double(b, "radius", 0.0);
        bc.depth = attrib_double(b, "depth", 0.0);
        bc.rounded_tip = attrib_bool(b, "roundedTip", false);
        gd.bore = bc;
    }
    if (const XmlNode* dl = d->first_node("DeadLayer")) {
        DeadLayerConfig dc;
        dc.front = attrib_double(dl, "front", 0.0);
        dc.side = attrib_double(dl, "side", 0.0);
        dc.back = attrib_double(dl, "back", 0.0);
        gd.dead_layer = dc;
    }
    for (const XmlNode* ln = d->first_node("Layer"); ln;
         ln = ln->next_sibling("Layer")) {
        LayerSpec l;
        l.material_index = static_cast<int>(attrib_double(ln, "material", -1));
        l.front_thickness_cm = attrib_double(ln, "frontThickness", 0.0);
        l.side_thickness_cm = attrib_double(ln, "sideThickness", 0.0);
        l.z_start_cm = attrib_double(ln, "zStart", 0.0);
        l.z_end_cm = attrib_double(ln, "zEnd", 0.0);
        gd.layers.push_back(l);
    }
    if (const XmlNode* c = d->first_node("Collimator")) {
        CollimatorSpec cs;
        cs.material_index = static_cast<int>(attrib_double(c, "material", -1));
        cs.side_thickness_cm = attrib_double(c, "sideThickness", 0.0);
        cs.z_start_cm = attrib_double(c, "zStart", 0.0);
        cs.z_end_cm = attrib_double(c, "zEnd", 0.0);
        gd.collimator = cs;
    }
    if (const XmlNode* mats = d->first_node("Materials")) {
        for (const XmlNode* mn = mats->first_node("Material"); mn;
             mn = mn->next_sibling("Material")) {
            MaterialSpec m;
            if (const char* v = attrib_value(mn, "name")) m.name = v;
            m.density_g_per_cm3 = attrib_double(mn, "density", 0.0);
            for (const XmlNode* el = mn->first_node("El"); el;
                 el = el->next_sibling("El")) {
                MaterialComponent c;
                c.Z = static_cast<uint8_t>(attrib_double(el, "Z", 0));
                c.mass_fraction = attrib_double(el, "frac", 0.0);
                m.composition.push_back(c);
            }
            gd.materials.push_back(std::move(m));
        }
    }

    return gd;
}

}  // namespace

std::string GeometryDescriptor::to_xml_string() const {
    XmlDoc doc;
    XmlNode* root = doc.allocate_node(rapidxml::node_element, "CeeLoGeometry");
    doc.append_node(root);
    write_detector_node(doc, root, *this);

    std::string out;
    rapidxml::print(std::back_inserter(out), doc, 0);
    return out;
}


GeometryDescriptor GeometryDescriptor::from_xml_string(const std::string& xml) {
    std::vector<char> buf(xml.begin(), xml.end());
    buf.push_back('\0');
    XmlDoc doc;
    doc.parse<rapidxml::parse_trim_whitespace>(buf.data());

    const XmlNode* root = doc.first_node("CeeLoGeometry");
    if (!root) throw std::runtime_error("CeeLoGeometry: missing root node");

    const XmlNode* d = root->first_node("Detector");
    if (!d) throw std::runtime_error("CeeLoGeometry: missing Detector node");

    return read_detector_node(d);
}


std::string DetectorResponse::to_xml_string() const {
    return serialize_xml(/*include_certificate=*/true);
}

std::string DetectorResponse::serialize_xml(bool include_certificate) const {
    XmlDoc doc;
    XmlNode* root = doc.allocate_node(rapidxml::node_element, "CeeLoResponse");
    doc.append_node(root);
    {
        // Lowest version that can read this file correctly -- see the comment on
        // sm_xmlSerializationVersion. Only a file that actually carries a v2
        // field is stamped v2, so responses that predate the fillet keep their
        // exact bytes (and content_hash) and stay loadable by older builds.
        const bool needs_v2 = (descriptor.bullet_radius_cm > 0.0)
                              || (descriptor.bore && descriptor.bore->rounded_tip);
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%i", needs_v2 ? 2 : 1);
        append_attrib(doc, root, "version", buf);
    }

    // Provenance
    {
        XmlNode* p = append_node(doc, root, "Provenance");
        append_attrib(doc, p, "ceeloVersion", provenance.ceelo_version);
        append_attrib(doc, p, "createdUtc", provenance.created_utc);
        append_attrib(doc, p, "profile", to_string(provenance.profile));
        append_attrib(doc, p, "nodeFepPrecision",
                      fmt_double(provenance.node_fep_precision));
        append_attrib(doc, p, "fepWindowKeV",
                      fmt_double(provenance.fep_window_keV));
        append_attrib(doc, p, "generationSeed",
                      std::to_string(provenance.generation_seed));
        append_attrib(doc, p, "kernelNRays",
                      std::to_string(provenance.kernel_n_rays));
        append_attrib(doc, p, "detectorName", provenance.detector_name);
        XmlNode* ve = append_node(doc, p, "ValidEnergy");
        append_attrib(doc, ve, "minKeV", fmt_double(provenance.valid_e_min_keV));
        append_attrib(doc, ve, "maxKeV", fmt_double(provenance.valid_e_max_keV));
        append_value_node(doc, p, "MinDistanceCm",
                          fmt_double(provenance.min_distance_cm));
        append_attrib(doc, p, "method", to_string(provenance.method));
    }

    // Detector geometry
    write_detector_node(doc, root, descriptor);

    // Mu tables
    {
        XmlNode* ms = append_node(doc, root, "MuTables");
        for (const MuTable& t : mu_tables) {
            XmlNode* mn = append_node(doc, ms, "MuTable");
            append_attrib(doc, mn, "material", std::to_string(t.material_index));
            append_value_node(doc, mn, "Energies", join_doubles(t.energy_keV));
            append_value_node(doc, mn, "MuPE", join_doubles(t.mu_pe));
            append_value_node(doc, mn, "MuCS", join_doubles(t.mu_cs));
            append_value_node(doc, mn, "MuRS", join_doubles(t.mu_rs));
            append_value_node(doc, mn, "MuPP", join_doubles(t.mu_pp));
        }
    }

    if (!eta_fep.empty())
        eta_to_xml(doc, root, "EtaFep", eta_fep);

    if (!near_field.empty()) {
        XmlNode* nf = append_node(doc, root, "NearField");
        append_value_node(doc, nf, "Energies",
                          join_doubles(near_field.energies_keV));
        append_value_node(doc, nf, "CosThetas",
                          join_doubles(near_field.cos_thetas));
        append_value_node(doc, nf, "DistsCm",
                          join_doubles(near_field.dists_cm));
        // ln N flattened [e][c][d] energy-major: index = (e*nc + c)*nd + d.
        append_value_node(doc, nf, "LnN", join_doubles(near_field.ln_n));
        append_value_node(doc, nf, "FracSigma",
                          join_doubles(near_field.frac_sigma));
        append_value_node(doc, nf, "BreakCosThetas",
                          join_doubles(near_field.break_cos_thetas));
        append_value_node(doc, nf, "BreakDistCm",
                          join_doubles(near_field.break_d_cm));
    }

    {
        XmlNode* te = append_node(doc, root, "TotalEfficiency");
        const char* tier = tot_eff.tier == TotEffTier::KernelExact ? "kernel"
                           : tot_eff.tier == TotEffTier::BCurve    ? "bcurve"
                                                                   : "etatable";
        append_attrib(doc, te, "tier", tier);
        if (scatter_in_recapture != 0.0)
            append_attrib(doc, te, "scatterInRecapture",
                          fmt_double(scatter_in_recapture));
        if (tot_eff.tier == TotEffTier::BCurve) {
            append_value_node(doc, te, "BEnergies",
                              join_doubles(tot_eff.b_energies_keV));
            append_value_node(doc, te, "LnB", join_doubles(tot_eff.ln_b));
        }
        if (tot_eff.tier == TotEffTier::EtaTotTable)
            eta_to_xml(doc, te, "EtaTot", tot_eff.eta_tot);
    }

    if (!grounding.empty()) {
        XmlNode* g = append_node(doc, root, "Grounding");
        append_attrib(doc, g, "curveDerived", grounding.curve_derived ? "1" : "0");
        append_value_node(doc, g, "KnotLnEnergies",
                          join_doubles(grounding.knot_ln_energies));
        append_value_node(doc, g, "LnK", join_doubles(grounding.ln_k));
        append_value_node(doc, g, "Cov", join_doubles(grounding.cov));
        {
            const SigmaTransferModel& t = grounding.transfer;
            XmlNode* tn = append_node(doc, g, "Transfer");
            append_attrib(doc, tn, "farOnAxis", fmt_double(t.far_onaxis));
            append_attrib(doc, tn, "offAxisMid", fmt_double(t.offaxis_mid));
            append_attrib(doc, tn, "offAxisLowE", fmt_double(t.offaxis_low_e));
            append_attrib(doc, tn, "lowERefKeV", fmt_double(t.low_e_ref_keV));
            append_attrib(doc, tn, "midERefKeV", fmt_double(t.mid_e_ref_keV));
            append_attrib(doc, tn, "nearContact", fmt_double(t.near_contact));
            append_attrib(doc, tn, "nearGateA", fmt_double(t.near_gate_a));
        }
        if (!grounding.points.empty()) {
            XmlNode* pts = append_node(doc, g, "Points");
            for (const GroundingPoint& p : grounding.points) {
                XmlNode* pn = append_node(doc, pts, "Pt");
                append_attrib(doc, pn, "E", fmt_double(p.energy_keV));
                append_attrib(doc, pn, "eff", fmt_double(p.measured_eff));
                append_attrib(doc, pn, "model", fmt_double(p.model_eff));
                append_attrib(doc, pn, "statSig", fmt_double(p.frac_stat_sigma));
                append_attrib(doc, pn, "certSig", fmt_double(p.frac_cert_sigma));
                append_attrib(doc, pn, "src", p.source_key);
                append_attrib(doc, pn, "d", fmt_double(p.distance_cm));
                append_attrib(doc, pn, "ct", fmt_double(p.cos_theta));
                append_attrib(doc, pn, "phi", fmt_double(p.phi_deg));
            }
        }
    }

    {
        XmlNode* f = append_node(doc, root, "SigmaFloors");
        append_attrib(doc, f, "fepFar", fmt_double(floors.fep_far));
        append_attrib(doc, f, "fepNear", fmt_double(floors.fep_near));
        append_attrib(doc, f, "totFar", fmt_double(floors.tot_far));
        append_attrib(doc, f, "totNear", fmt_double(floors.tot_near));
        append_attrib(doc, f, "nearRegimeA", fmt_double(floors.near_regime_a));
    }

    // Unconditional transfer-sigma envelope (EFFTRAN transfer responses).
    if (model_transfer) {
        const SigmaTransferModel& t = *model_transfer;
        XmlNode* tn = append_node(doc, root, "ModelTransfer");
        append_attrib(doc, tn, "farOnAxis", fmt_double(t.far_onaxis));
        append_attrib(doc, tn, "offAxisMid", fmt_double(t.offaxis_mid));
        append_attrib(doc, tn, "offAxisLowE", fmt_double(t.offaxis_low_e));
        append_attrib(doc, tn, "lowERefKeV", fmt_double(t.low_e_ref_keV));
        append_attrib(doc, tn, "midERefKeV", fmt_double(t.mid_e_ref_keV));
        append_attrib(doc, tn, "nearContact", fmt_double(t.near_contact));
        append_attrib(doc, tn, "nearGateA", fmt_double(t.near_gate_a));
    }

    // Accuracy certificate -- additive METADATA, excluded from content_hash
    // (serialize_xml(false)). Old readers ignore the unknown <Certificate>
    // element; new readers tolerate its absence (certificate.empty()). The
    // Rows payload is a flat, fixed-width (10 columns) double stream so the
    // existing join_doubles/parse_doubles helpers give string-identical
    // save/load/save.
    if (include_certificate && !certificate.empty()) {
        const AccuracyCertificate& c = certificate;
        XmlNode* cn = append_node(doc, root, "Certificate");
        append_attrib(doc, cn, "converged", c.converged ? "1" : "0");
        append_attrib(doc, cn, "iterations", std::to_string(c.iterations));
        append_attrib(doc, cn, "cpuSeconds", fmt_double(c.cpu_seconds));
        append_attrib(doc, cn, "probeSeed", std::to_string(c.probe_seed_base));
        XmlNode* sn = append_node(doc, cn, "Summary");
        append_attrib(doc, sn, "fepMed", fmt_double(c.fep_median));
        append_attrib(doc, sn, "fepP95", fmt_double(c.fep_p95));
        append_attrib(doc, sn, "fepMax", fmt_double(c.fep_max));
        append_attrib(doc, sn, "totMed", fmt_double(c.tot_median));
        append_attrib(doc, sn, "totP95", fmt_double(c.tot_p95));
        std::vector<double> flat;
        flat.reserve(c.rows.size() * 10);
        for (const AccuracyCertificate::Row& r : c.rows) {
            flat.push_back(r.E_keV);      flat.push_back(r.d_cm);
            flat.push_back(r.cos_theta);  flat.push_back(r.phi_deg);
            flat.push_back(r.mc);         flat.push_back(r.mc_sig);
            flat.push_back(r.model);      flat.push_back(r.model_sig);
            flat.push_back(static_cast<double>(r.tag));
            flat.push_back(r.pass ? 1.0 : 0.0);
        }
        append_value_node(doc, cn, "Rows", join_doubles(flat));
    }

    std::string out;
    rapidxml::print(std::back_inserter(out), doc, 0);
    return out;
}

std::shared_ptr<DetectorResponse> DetectorResponse::from_xml_string(
    const std::string& xml) {
    std::vector<char> buf(xml.begin(), xml.end());
    buf.push_back('\0');
    XmlDoc doc;
    doc.parse<rapidxml::parse_trim_whitespace>(buf.data());

    const XmlNode* root = doc.first_node("CeeLoResponse");
    if (!root)
        throw std::runtime_error("CeeLoResponse: missing root node");
    const char* ver = attrib_value(root, "version");
    const int version = ver ? std::atoi(ver) : -1;
    if (version < 1 || version > sm_xmlSerializationVersion)
        throw std::runtime_error("CeeLoResponse: unsupported version");

    auto resp = std::make_shared<DetectorResponse>();

    if (const XmlNode* p = root->first_node("Provenance")) {
        ResponseProvenance& pr = resp->provenance;
        if (const char* v = attrib_value(p, "ceeloVersion")) pr.ceelo_version = v;
        if (const char* v = attrib_value(p, "createdUtc")) pr.created_utc = v;
        pr.profile = profile_from_string(attrib_value(p, "profile"));
        pr.node_fep_precision = attrib_double(p, "nodeFepPrecision", 0.003);
        pr.fep_window_keV =
            attrib_double(p, "fepWindowKeV", kDefaultFepWindowKeV);
        if (const char* v = attrib_value(p, "generationSeed"))
            pr.generation_seed = std::strtoull(v, nullptr, 10);
        pr.kernel_n_rays =
            static_cast<int>(attrib_double(p, "kernelNRays", 2048));
        if (const char* v = attrib_value(p, "detectorName")) pr.detector_name = v;
        if (const XmlNode* ve = p->first_node("ValidEnergy")) {
            pr.valid_e_min_keV = attrib_double(ve, "minKeV", 0.0);
            pr.valid_e_max_keV = attrib_double(ve, "maxKeV", 0.0);
        }
        if (const char* v = child_value(p, "MinDistanceCm"))
            pr.min_distance_cm = std::strtod(v, nullptr);
        if (const char* v = attrib_value(p, "method")) {
            if (std::strcmp(v, "quick_mc_transfer") == 0)
                pr.method = ProductionMethod::QuickMcTransfer;
            else if (std::strcmp(v, "curve_transfer") == 0)
                pr.method = ProductionMethod::CurveTransfer;
            else
                pr.method = ProductionMethod::FullMc;
        }
    }

    const XmlNode* d = root->first_node("Detector");
    if (!d) throw std::runtime_error("CeeLoResponse: missing Detector node");
    resp->descriptor = read_detector_node(d);

    if (const XmlNode* ms = root->first_node("MuTables")) {
        for (const XmlNode* mn = ms->first_node("MuTable"); mn;
             mn = mn->next_sibling("MuTable")) {
            MuTable t;
            t.material_index = static_cast<int>(attrib_double(mn, "material", -1));
            t.energy_keV = parse_doubles(child_value(mn, "Energies"));
            t.mu_pe = parse_doubles(child_value(mn, "MuPE"));
            t.mu_cs = parse_doubles(child_value(mn, "MuCS"));
            t.mu_rs = parse_doubles(child_value(mn, "MuRS"));
            t.mu_pp = parse_doubles(child_value(mn, "MuPP"));
            resp->mu_tables.push_back(std::move(t));
        }
    }

    if (const XmlNode* n = root->first_node("EtaFep"))
        eta_from_xml(n, resp->eta_fep);

    if (const XmlNode* nf = root->first_node("NearField")) {
        NearFieldModel& m = resp->near_field;
        m.energies_keV = parse_doubles(child_value(nf, "Energies"));
        m.cos_thetas = parse_doubles(child_value(nf, "CosThetas"));
        m.dists_cm = parse_doubles(child_value(nf, "DistsCm"));
        m.ln_n = parse_doubles(child_value(nf, "LnN"));
        m.frac_sigma = parse_doubles(child_value(nf, "FracSigma"));
        if (m.ln_n.size() !=
            m.energies_keV.size() * m.cos_thetas.size() * m.dists_cm.size())
            throw std::runtime_error("CeeLoResponse: NearField LnN size");
        if (m.frac_sigma.size() != m.ln_n.size())
            throw std::runtime_error("CeeLoResponse: NearField FracSigma size");
        m.break_cos_thetas = parse_doubles(child_value(nf, "BreakCosThetas"));
        m.break_d_cm = parse_doubles(child_value(nf, "BreakDistCm"));
    }

    if (const XmlNode* te = root->first_node("TotalEfficiency")) {
        const char* tier = attrib_value(te, "tier");
        if (tier && std::strcmp(tier, "bcurve") == 0)
            resp->tot_eff.tier = TotEffTier::BCurve;
        else if (tier && std::strcmp(tier, "etatable") == 0)
            resp->tot_eff.tier = TotEffTier::EtaTotTable;
        else
            resp->tot_eff.tier = TotEffTier::KernelExact;
        resp->scatter_in_recapture = attrib_double(te, "scatterInRecapture", 0.0);
        resp->tot_eff.b_energies_keV = parse_doubles(child_value(te, "BEnergies"));
        resp->tot_eff.ln_b = parse_doubles(child_value(te, "LnB"));
        if (const XmlNode* et = te->first_node("EtaTot"))
            eta_from_xml(et, resp->tot_eff.eta_tot);
    }

    if (const XmlNode* g = root->first_node("Grounding")) {
        GroundingBlock& gb = resp->grounding;
        const char* cd = attrib_value(g, "curveDerived");
        gb.curve_derived = cd && std::strcmp(cd, "1") == 0;
        gb.knot_ln_energies = parse_doubles(child_value(g, "KnotLnEnergies"));
        gb.ln_k = parse_doubles(child_value(g, "LnK"));
        gb.cov = parse_doubles(child_value(g, "Cov"));
        if (const XmlNode* tn = g->first_node("Transfer")) {
            SigmaTransferModel& t = gb.transfer;
            t.far_onaxis = attrib_double(tn, "farOnAxis", t.far_onaxis);
            t.offaxis_mid = attrib_double(tn, "offAxisMid", t.offaxis_mid);
            t.offaxis_low_e = attrib_double(tn, "offAxisLowE", t.offaxis_low_e);
            t.low_e_ref_keV = attrib_double(tn, "lowERefKeV", t.low_e_ref_keV);
            t.mid_e_ref_keV = attrib_double(tn, "midERefKeV", t.mid_e_ref_keV);
            t.near_contact = attrib_double(tn, "nearContact", t.near_contact);
            t.near_gate_a = attrib_double(tn, "nearGateA", t.near_gate_a);
        }
        if (const XmlNode* pts = g->first_node("Points")) {
            for (const XmlNode* pn = pts->first_node("Pt"); pn;
                 pn = pn->next_sibling("Pt")) {
                GroundingPoint p;
                p.energy_keV = attrib_double(pn, "E", 0.0);
                p.measured_eff = attrib_double(pn, "eff", 0.0);
                p.model_eff = attrib_double(pn, "model", 0.0);
                p.frac_stat_sigma = attrib_double(pn, "statSig", 0.0);
                p.frac_cert_sigma = attrib_double(pn, "certSig", 0.0);
                if (const char* v = attrib_value(pn, "src")) p.source_key = v;
                p.distance_cm = attrib_double(pn, "d", 0.0);
                p.cos_theta = attrib_double(pn, "ct", 1.0);
                p.phi_deg = attrib_double(pn, "phi", 0.0);
                gb.points.push_back(std::move(p));
            }
        }
    }

    if (const XmlNode* f = root->first_node("SigmaFloors")) {
        SigmaFloors& fl = resp->floors;
        fl.fep_far = attrib_double(f, "fepFar", fl.fep_far);
        fl.fep_near = attrib_double(f, "fepNear", fl.fep_near);
        fl.tot_far = attrib_double(f, "totFar", fl.tot_far);
        fl.tot_near = attrib_double(f, "totNear", fl.tot_near);
        fl.near_regime_a = attrib_double(f, "nearRegimeA", fl.near_regime_a);
    }

    if (const XmlNode* tn = root->first_node("ModelTransfer")) {
        SigmaTransferModel t;
        t.far_onaxis = attrib_double(tn, "farOnAxis", t.far_onaxis);
        t.offaxis_mid = attrib_double(tn, "offAxisMid", t.offaxis_mid);
        t.offaxis_low_e = attrib_double(tn, "offAxisLowE", t.offaxis_low_e);
        t.low_e_ref_keV = attrib_double(tn, "lowERefKeV", t.low_e_ref_keV);
        t.mid_e_ref_keV = attrib_double(tn, "midERefKeV", t.mid_e_ref_keV);
        t.near_contact = attrib_double(tn, "nearContact", t.near_contact);
        t.near_gate_a = attrib_double(tn, "nearGateA", t.near_gate_a);
        resp->model_transfer = t;
    }

    // Accuracy certificate (optional; absence tolerated -> empty()).
    if (const XmlNode* cn = root->first_node("Certificate")) {
        AccuracyCertificate& c = resp->certificate;
        const char* cv = attrib_value(cn, "converged");
        c.converged = cv && std::strcmp(cv, "1") == 0;
        c.iterations = static_cast<int>(attrib_double(cn, "iterations", 0.0));
        c.cpu_seconds = attrib_double(cn, "cpuSeconds", 0.0);
        if (const char* v = attrib_value(cn, "probeSeed"))
            c.probe_seed_base = std::strtoull(v, nullptr, 10);
        if (const XmlNode* sn = cn->first_node("Summary")) {
            c.fep_median = attrib_double(sn, "fepMed", 0.0);
            c.fep_p95 = attrib_double(sn, "fepP95", 0.0);
            c.fep_max = attrib_double(sn, "fepMax", 0.0);
            c.tot_median = attrib_double(sn, "totMed", 0.0);
            c.tot_p95 = attrib_double(sn, "totP95", 0.0);
        }
        const std::vector<double> flat = parse_doubles(child_value(cn, "Rows"));
        for (size_t i = 0; i + 10 <= flat.size(); i += 10) {
            AccuracyCertificate::Row r;
            r.E_keV = flat[i];      r.d_cm = flat[i + 1];
            r.cos_theta = flat[i + 2]; r.phi_deg = flat[i + 3];
            r.mc = flat[i + 4];     r.mc_sig = flat[i + 5];
            r.model = flat[i + 6];  r.model_sig = flat[i + 7];
            r.tag = static_cast<uint8_t>(flat[i + 8]);
            r.pass = flat[i + 9] != 0.0;
            c.rows.push_back(r);
        }
    }

    resp->finalize();
    return resp;
}

uint64_t DetectorResponse::content_hash() const {
    // The certificate is metadata ABOUT the response, not content -- hash the
    // certificate-free payload so presence/absence never changes the hash.
    return fnv1a_64(serialize_xml(/*include_certificate=*/false));
}

} // namespace ceelo
