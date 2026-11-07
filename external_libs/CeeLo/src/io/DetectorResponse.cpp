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
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

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
// and their hash (content_hash hashes an optional field only when it is present).
// Follow this pattern for the next additive field.
//
// v3 followed the pattern for TotEffTier::NotCharacterized: a pre-v3 reader maps
// the unknown tier string onto KernelExact and would serve a bare kernel as a
// verified total efficiency, so only responses that use that tier are stamped
// v3.
//
// v4 = side & back coverage: eta tables with nodes below cos theta = 0.  A pre-v4 reader would
// interpolate one PCHIP across the front and back node sets (EtaTable::eval_ln keeps them apart),
// changing values in FRONT of the face too, so only covered responses are stamped v4.
const int DetectorResponse::sm_xmlSerializationVersion = 4;

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

/// 64-bit FNV-1a over the bytes of a response's stored content; see DetectorResponse::content_hash.
/// Each variable-length field is preceded by its length, and an optional part contributes nothing
/// at all when absent and a #tag naming it ahead of its content when present - so a field added
/// later, absent from every existing response, leaves their hashes unchanged.  Numbers are hashed in
/// native byte order - little-endian on every platform built for.
class ContentHasher {
public:
    /// Names the optional part whose content follows; unique within the response.
    void tag(const char* name) { text(name); }
    void flag(bool b) {
        const unsigned char c = b ? 1 : 0;
        bytes(&c, 1);
    }
    void integer(int64_t v) { bytes(&v, sizeof(v)); }
    void number(double v) {
        // A NaN's payload is not stored (the text form is just "nan"), so it must not count.
        if (std::isnan(v)) v = std::numeric_limits<double>::quiet_NaN();
        bytes(&v, sizeof(v));
    }
    void numbers(const std::vector<double>& v) {
        integer(static_cast<int64_t>(v.size()));
        for (const double d : v) number(d);
    }
    void text(const std::string& s) {
        integer(static_cast<int64_t>(s.size()));
        bytes(s.data(), s.size());
    }
    uint64_t value() const { return h_; }

private:
    void bytes(const void* p, size_t n) {
        const unsigned char* c = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) {
            h_ ^= c[i];
            h_ *= 1099511628211ULL;
        }
    }
    uint64_t h_ = 14695981039346656037ULL;
};

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

const char* to_string(PrecisionProfile p) {
    switch (p) {
        case PrecisionProfile::Uniform:   return "uniform";
        case PrecisionProfile::RelaxMild: return "relax_mild";
    }
    return "?";
}

const char* to_string(BackCoverage b) {
    switch (b) {
        case BackCoverage::None:   return "none";
        case BackCoverage::Coarse: return "coarse";
        case BackCoverage::Fine:   return "fine";
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

// The recorded generation choices are optional: absent or unrecognised reads as not recorded.
std::optional<PrecisionProfile> precision_profile_from_string(const char* s) {
    for (const PrecisionProfile p : {PrecisionProfile::Uniform, PrecisionProfile::RelaxMild})
        if (s && (std::strcmp(s, to_string(p)) == 0)) return p;
    return std::nullopt;
}

std::optional<BackCoverage> back_coverage_from_string(const char* s) {
    for (const BackCoverage b : {BackCoverage::None, BackCoverage::Coarse, BackCoverage::Fine})
        if (s && (std::strcmp(s, to_string(b)) == 0)) return b;
    return std::nullopt;
}
}  // namespace

const char* to_string(ProductionMethod m) {
    switch (m) {
        case ProductionMethod::FullMc:          return "full_mc";
        case ProductionMethod::QuickMcTransfer: return "quick_mc_transfer";
        case ProductionMethod::CurveTransfer:   return "curve_transfer";
    }
    return "full_mc";
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
    // Every reader takes cos_thetas.front() / .back() once the table has energies.
    if (nc == 0)
        throw std::runtime_error("EtaTable: energies but no cos-theta nodes");
    if (ln_eta.size() != ne * nc * np)
        throw std::runtime_error("EtaTable: ln_eta size mismatch");
    if (frac_sigma.size() != ln_eta.size())
        throw std::runtime_error("EtaTable: frac_sigma size mismatch");

    // Normalize the edge list before anything reads it positionally.  `edges_keV`
    //  is a public member and the deserializer takes it verbatim from the XML, so
    //  an ascending, de-duplicated, in-range list cannot be assumed - while
    //  `segment_of` only counts edges below E (order-agnostic), the stub-reach
    //  test below indexes edges[s-1]/edges[s] as the segment's BOUNDS, which is
    //  only meaningful when sorted.  Unsorted input could otherwise invert those
    //  bounds and silently restore the fabricated-discontinuity bug this merge
    //  exists to prevent.  Edges at or outside the node span are dropped: every
    //  node lies on one side of such an edge, so it cannot separate any nodes,
    //  and keeping it would only corrupt a neighbouring segment's bounds.
    std::sort(edges_keV.begin(), edges_keV.end());
    edges_keV.erase(std::unique(edges_keV.begin(), edges_keV.end(),
                                [](double a, double b) { return b - a < 0.01; }),
                    edges_keV.end());
    edges_keV.erase(std::remove_if(edges_keV.begin(), edges_keV.end(),
                                   [this](double e) {
                                       return !(e > energies_keV.front())
                                              || !(e < energies_keV.back());
                                   }),
                    edges_keV.end());

    // Group energy nodes into K-edge segments (both flanks are nodes).  The
    //  runs depend only on the energy axis, not on (c, p), so compute them once.
    std::vector<size_t> seg_id(ne);
    for (size_t e = 0; e < ne; ++e)
        seg_id[e] = segment_of(energies_keV[e], edges_keV);

    std::vector<std::pair<size_t, size_t>> runs;   // inclusive [first, last]
    for (size_t e0 = 0; e0 < ne;) {
        size_t e1 = e0;
        while (e1 + 1 < ne && seg_id[e1 + 1] == seg_id[e0]) ++e1;
        runs.emplace_back(e0, e1);
        e0 = e1 + 1;
    }

    // Fold a single-node run into a neighbour (following preferred, else
    //  preceding) WHEN its lone node cannot stand in for its whole segment.
    //  A one-node run becomes a constant stub below, and `SegCurve::eval`
    //  clamps to the nearest segment, so that node's value is returned across
    //  the segment's entire width.  Whether that is faithful or fabricated
    //  depends on how far it has to reach:
    //
    //   - A lone node that IS a flank sits ~0.1% from the edge it brackets, and
    //     pinning the edge value from one side is exactly what a flank is for.
    //     The stub is then correct to within the flank width, and the clamp
    //     reproduces the discontinuity - keep it.
    //   - A lone node far from its segment's bounds is a fabrication.  Seen in
    //     the wild: a 10/12 keV node pair with a retained 11.107 keV Ge K-edge
    //     between them left 10 keV alone in its segment, so eval froze at the
    //     10 keV value below the geometric crossover and at the 12 keV value
    //     above it - a hard 95x jump at 10.95 keV that is nowhere in the data
    //     (the two source pages were bit-identical).
    //
    //  `kStubReach` is 10x the 1e-3 relative flank half-width producers use, so
    //  a somewhat wider flank convention still keeps its stub.  A lone node
    //  bounded by two edges cannot pin both, and merging it does bridge an
    //  edge; a monotone interpolant spanning the edge is a weaker guarantee
    //  than segmentation, but it is continuous and bounded by its node values,
    //  which freezing flat and inventing a jump is not.  A genuinely
    //  single-node table (one run) always keeps the stub: constant is then the
    //  only defensible answer.
    const double kStubReach = 1.0e-2;
    for (size_t i = 0; runs.size() > 1 && i < runs.size();) {
        if (runs[i].second > runs[i].first) { ++i; continue; }

        // How far this node must be extrapolated to cover its segment.
        const double e_lone = energies_keV[runs[i].first];
        const size_t s = seg_id[runs[i].first];
        const double lo = (s == 0) ? energies_keV.front() : edges_keV[s - 1];
        const double hi = (s >= edges_keV.size()) ? energies_keV.back()
                                                  : edges_keV[s];
        const double reach = std::max(e_lone / std::max(lo, 1e-300),
                                      std::max(hi, 0.0) / std::max(e_lone, 1e-300))
                             - 1.0;
        if (reach <= kStubReach) { ++i; continue; }

        if (i + 1 < runs.size())
            runs[i + 1].first = runs[i].first;
        else
            runs[i - 1].second = runs[i].second;
        runs.erase(runs.begin() + static_cast<long>(i));
    }

    curves_.resize(nc * np);
    for (size_t c = 0; c < nc; ++c) {
        for (size_t p = 0; p < np; ++p) {
            SegCurve& curve = curves_[c * np + p];
            curve.segs.reserve(runs.size());
            curve.seg_lo.reserve(runs.size());
            curve.seg_hi.reserve(runs.size());
            for (const std::pair<size_t, size_t>& run : runs) {
                std::vector<double> xs, ys;
                xs.reserve(run.second - run.first + 2);
                ys.reserve(run.second - run.first + 2);
                for (size_t e = run.first; e <= run.second; ++e) {
                    xs.push_back(std::log(energies_keV[e]));
                    ys.push_back(ln_eta[index(e, c, p)]);
                }
                if (xs.size() == 1) {  // lone node table: constant stub segment
                    xs.push_back(xs[0] + 1e-9);
                    ys.push_back(ys[0]);
                }
                curve.seg_lo.push_back(xs.front());
                curve.seg_hi.push_back(xs.back());
                curve.segs.emplace_back(std::move(xs), std::move(ys));
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

    // Side & back coverage (nodes below cos theta = 0): the front nodes and the back nodes are
    //  interpolated as separate sets, so the front of a covered table evaluates exactly as a table
    //  without coverage does - including below its grazing node, where both clamp to it.  Behind
    //  the face plane, the gap from the last back node up to the plane is bridged linearly to that
    //  clamped grazing value.
    size_t c_lo = 0, c_hi = nc;
    bool gap = false;
    if ((nc > 1) && (cos_thetas.front() < 0.0)) {
        size_t k0 = 0;
        while ((k0 < nc) && (cos_thetas[k0] < 0.0))
            ++k0;
        if ((k0 < nc) && (ct >= cos_thetas[k0])) {
            c_lo = k0;
        } else if ((k0 > 0) && (ct <= cos_thetas[k0 - 1])) {
            c_hi = k0;
        } else if ((k0 < nc) && (ct >= 0.0)) {
            c_lo = k0;
            c_hi = k0 + 1;
            clamped = true;
        } else if ((k0 > 0) && (k0 < nc)) {
            gap = true;
            c_lo = k0 - 1;
            c_hi = k0 + 1;
        }
    }

    // Per phi node: energy eval per ct node, then PCHIP across cos-theta.
    auto eval_at_phi = [&](size_t p) -> double {
        const size_t n = c_hi - c_lo;
        std::vector<double> y(n);
        for (size_t c = c_lo; c < c_hi; ++c)
            y[c - c_lo] = curves_[c * np + p].eval(lnE, clamped);
        if (n == 1) return y[0];
        if (gap) {
            const double x0 = cos_thetas[c_lo];   // < 0; the grazing value holds from 0 up
            const double f = (ct - x0) / (0.0 - x0);
            return (1.0 - f) * y[0] + f * y[1];
        }
        if (n == nc) {
            const Pchip ang(cos_thetas, y);
            return ang(ct);
        }
        const std::vector<double> xs(cos_thetas.begin() + static_cast<std::ptrdiff_t>(c_lo),
                                     cos_thetas.begin() + static_cast<std::ptrdiff_t>(c_hi));
        const Pchip ang(xs, y);
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

bool EtaTable::is_angle_flat() const {
    const size_t nc = cos_thetas.size();
    const size_t np = phis_deg.empty() ? 1 : phis_deg.size();
    if (nc < 2) return true;
    for (size_t e = 0; e < energies_keV.size(); ++e)
        for (size_t p = 0; p < np; ++p)
            for (size_t c = 1; c < nc; ++c)
                if (ln_eta[index(e, c, p)] != ln_eta[index(e, 0, p)])
                    return false;
    return true;
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

    // In front of the face plane a covered table brackets only its front nodes, as a table without
    //  coverage does (clamping below the grazing node).
    size_t k0 = 0;
    if (cos_theta >= 0.0) {
        while ((k0 < cos_thetas.size()) && (cos_thetas[k0] < 0.0)) ++k0;
        if (k0 == cos_thetas.size()) k0 = 0;
    }
    if (k0 == 0) {
        bracket(cos_thetas, cos_theta, c0, c1, wc);
    } else {
        const std::vector<double> front(cos_thetas.begin() + static_cast<std::ptrdiff_t>(k0),
                                        cos_thetas.end());
        bracket(front, cos_theta, c0, c1, wc);
        c0 += k0;
        c1 += k0;
    }

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

    // In front of the face plane a covered table brackets only its front nodes, as a table without
    //  coverage does (clamping below the grazing node).
    size_t k0 = 0;
    if (cos_theta >= 0.0) {
        while ((k0 < cos_thetas.size()) && (cos_thetas[k0] < 0.0)) ++k0;
        if (k0 == cos_thetas.size()) k0 = 0;
    }
    if (k0 == 0) {
        bracket(cos_thetas, cos_theta, c0, c1, wc);
    } else {
        const std::vector<double> front(cos_thetas.begin() + static_cast<std::ptrdiff_t>(k0),
                                        cos_thetas.end());
        bracket(front, cos_theta, c0, c1, wc);
        c0 += k0;
        c1 += k0;
    }
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
    near_field.finalize();  // tolerates empty (no near-field total)
}

double TotEffPayload::ln_b_at(double energy_keV) const {
    if (!b_curve_.valid()) return 0.0;
    return b_curve_(std::log(std::max(energy_keV, kLnClampEps)));
}

// ---------------------------------------------------------------------------
// Grounding
// ---------------------------------------------------------------------------

bool crystal_is_cdte_class(const GeometryDescriptor& gd) {
    if (gd.crystal_material_index < 0 ||
        static_cast<size_t>(gd.crystal_material_index) >= gd.materials.size())
        return false;
    double cd_te = 0.0;
    for (const MaterialComponent& c : gd.materials[static_cast<size_t>(
             gd.crystal_material_index)].composition) {
        if (c.Z == 48 || c.Z == 52) cd_te += c.mass_fraction;
    }
    return cd_te > 0.5;
}

SigmaTransferModel::Components SigmaTransferModel::components(
    double d_over_a, double cos_theta, double energy_keV) const {
    Components c;
    c.far_onaxis = far_onaxis;
    // Behind the face plane the crystal is seen from the side or back - the regime the off-axis
    //  envelope saturates toward - so sin^2 stays at its 90-degree value rather than falling
    //  back to 0 as if the source were on axis again.
    const double s2 = (cos_theta < 0.0) ? 1.0 : std::max(0.0, 1.0 - cos_theta * cos_theta);
    // Low-E weight: 0 at/above mid_e_ref, 1 at/below low_e_ref (ln ramp).
    double w = 0.0;
    if (energy_keV < mid_e_ref_keV) {
        w = (std::log(mid_e_ref_keV) - std::log(std::max(energy_keV, 1.0))) /
            (std::log(mid_e_ref_keV) - std::log(low_e_ref_keV));
        w = clamp(w, 0.0, 1.0);
    }
    // Saturating in sin^2(theta): the measured residual flattens past ~30 degrees
    // because beyond there you are looking at the crystal from the side and the
    // distribution of path lengths stops changing much with further angle.
    // offaxis_s2_half <= 0 selects the legacy unbounded form (old stored files).
    const double amp = offaxis_mid + offaxis_low_e * w * w;
    c.offaxis = (offaxis_s2_half > 0.0) ? (amp * s2 / (s2 + offaxis_s2_half))
                                        : (amp * s2);
    if (d_over_a < near_gate_a) {
        const double t = clamp((near_gate_a - d_over_a) / (near_gate_a - 1.0),
                               0.0, 1.0);
        c.near_field = near_contact * t;
    }
    return c;
}

double SigmaTransferModel::eval(double d_over_a, double cos_theta,
                                double energy_keV) const {
    const Components c = components(d_over_a, cos_theta, energy_keV);
    return std::sqrt(c.far_onaxis * c.far_onaxis + c.offaxis * c.offaxis +
                     c.near_field * c.near_field);
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

namespace {
/// The STORED (generation-time) cross sections of every material `q` crosses, indexed by
/// RaySegment::material_slot - looked up once per call rather than once per ray segment.
std::vector<MacroscopicXS> slot_cross_sections(
    double energy_keV, const ApertureQuadrature& q,
    const std::vector<std::pair<const Material*, size_t>>& mat_to_mu,
    const std::vector<MuTable>& mu_tables, const char* caller) {
    std::vector<MacroscopicXS> xs;
    xs.reserve(q.materials.size());
    for (const Material* m : q.materials) {
        size_t table = SIZE_MAX;
        for (const auto& mm : mat_to_mu)
            if (mm.first == m) { table = mm.second; break; }
        if (table == SIZE_MAX)
            throw std::runtime_error(std::string(caller) + ": unknown material");
        xs.push_back(mu_tables[table].eval(energy_keV));
    }
    return xs;
}
}  // namespace

void DetectorResponse::kernel_ray_weights_impl(
    double energy_keV, const ApertureQuadrature& q, MuChoice mu, double recap,
    std::vector<double>& w_out, std::vector<Eigen::Vector3d>& dirs_out) const {
    // Same recurrence as ApertureQuadrature::interaction_omega, but with the
    // STORED generation-time mu tables (self-containment guarantee).
    struct MuT { double tot, nors, cs; };
    std::vector<MuT> mus;  // by RaySegment::material_slot
    for (const MacroscopicXS& xs : slot_cross_sections(energy_keV, q, mat_to_mu_, mu_tables,
                                                       "DetectorResponse::kernel_K"))
        mus.push_back({xs.mu_total(), xs.mu_pe + xs.mu_cs + xs.mu_pp, xs.mu_cs});

    w_out.clear();
    dirs_out.clear();
    w_out.reserve(q.rays.size());
    dirs_out.reserve(q.rays.size());

    for (const KernelRay& r : q.rays) {
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;
        double p_int = 0.0;
        for (const RaySegment& s : q.segs(r)) {
            const MuT& m = mus[s.material_slot];
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
    std::vector<MuT> mus;  // by RaySegment::material_slot
    for (const MacroscopicXS& xs : slot_cross_sections(energy_keV, q, mat_to_mu_, mu_tables,
                                                       "DetectorResponse::kernel_K"))
        mus.push_back({xs.mu_total(), xs.mu_pe + xs.mu_cs + xs.mu_pp, xs.mu_cs});

    double total = 0.0;
    for (const KernelRay& r : q.rays) {
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;
        double p_int = 0.0;
        for (const RaySegment& s : q.segs(r)) {
            const MuT& m = mus[s.material_slot];
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

void DetectorResponse::fep_line_probabilities(
    double energy_keV, const ApertureQuadrature& q, std::vector<double>& p_out) const {
    // Same mu lookup and chord model as kernel_ray_weights_impl (FEP: total mu on the scoring
    // segments, no scatter-in credit), but one entry per ray, in order.
    std::vector<double> mu_tots;  // by RaySegment::material_slot
    for (const MacroscopicXS& xs : slot_cross_sections(energy_keV, q, mat_to_mu_, mu_tables,
                                                       "DetectorResponse::fep_line_probabilities"))
        mu_tots.push_back(xs.mu_total());

    p_out.assign(q.rays.size(), 0.0);
    for (size_t i = 0; i < q.rays.size(); ++i) {
        const KernelRay& r = q.rays[i];
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;
        double p_int = 0.0;
        for (const RaySegment& s : q.segs(r)) {
            const double mu_tot = mu_tots[s.material_slot];
            if (s.is_scoring)
                p_int += std::exp(-tau_before) * (1.0 - std::exp(-mu_tot * s.length));
            tau_before += mu_tot * s.length;
        }
        p_out[i] = p_int;
    }
}

void DetectorResponse::total_ray_weights(
    double energy_keV, const ApertureQuadrature& q,
    std::vector<double>& w_out, std::vector<Eigen::Vector3d>& dirs_out) const {
    // Mirrors eps_total_impl's tier dispatch - the EtaTotTable tier folds a measured eta_tot over
    // the FULL-mu kernel, the other tiers use the Rayleigh-free one.
    if (!tot_eff.characterized()) {
        // No total to differentiate through; empty weights make K == 0, matching
        // eps_total_impl's refusal instead of implying a bare-kernel total.
        w_out.clear();
        dirs_out.clear();
        return;
    }
    const MuChoice mu = (tot_eff.tier == TotEffTier::EtaTotTable) ? MuChoice::Total
                                                                  : MuChoice::NoRayleigh;
    const double recap = (mu == MuChoice::NoRayleigh)
        ? std::max(0.0, std::min(1.0, scatter_in_recapture)) : 0.0;
    kernel_ray_weights_impl(energy_keV, q, mu, recap, w_out, dirs_out);
}

// ---------------------------------------------------------------------------
// DetectorResponse -- evaluation
// ---------------------------------------------------------------------------

/// The fractional 1-sigma budget of one efficiency query, split by what it is made of: the
/// data-derived variances (the MC's / the user's own statistics - trustable) and the model
/// envelopes (model_sigma constants - how wrong the model may be where it is not measured).
/// Every envelope is a common mode across energies at the query geometry, so frac_covariance
/// builds its off-diagonals from exactly the terms the diagonal is made of.  Fixed size: this
/// sits in the per-element hot path of the volumetric integrals.
struct DetectorResponse::SigmaBudget {
    // Data-derived
    double node2 = 0.0;     // MC node variance (eta + near-field tables); independent per energy
    double var_ln_k = 0.0;  // anchor / grounding fit variance; off-diagonals via cov_ln_k
    // Model envelopes, each fully correlated across energies
    enum Term : size_t { Floor, BehindPlane, Shadow, NearUnmodeled, Buildup,
                         ModelFar, ModelOff, ModelNear, GroundFar, GroundOff, GroundNear,
                         BehindExtrapolation, BackStructureApprox,
                         NumTerms };
    std::array<double, NumTerms> model{};

    // `var_ln_k` comes from a stored covariance matrix that nothing validates as positive
    // semi-definite (it is read verbatim from XML, or handed over by set_anchor_covariance), so
    // clamp rather than let a bad file turn into a NaN sigma halfway through somebody's fit.
    double data2() const { return node2 + std::max(0.0, var_ln_k); }
    double model2() const {
        double s = 0.0;
        for (const double m : model) s += m * m;
        return s;
    }
    double frac2() const { return data2() + model2(); }
};

struct DetectorResponse::EvalCommon {
    double d_cm = 0.0;          // from crystal-face origin
    double cos_theta = 1.0;
    double phi_deg = 0.0;
    double a_cm = 0.0;          // transverse half-extent
    // The distance the near-field gates measure: d_cm in front of the face plane; behind it,
    //  a + the clearance past the detector's outside - from the face-centre origin a source
    //  touching the side wall or the back is d >~ L away, which says nothing about how near it is.
    double d_near_cm = 0.0;
    bool near_regime = false;
    ResponseFlag flag = ResponseFlag::Ok;
    SigmaBudget budget;         // geometry-only terms filled here; the rest per quantity
};

namespace {
void raise_flag(ResponseFlag& current, ResponseFlag f) {
    if (static_cast<int>(f) > static_cast<int>(current)) current = f;
}

/// The cos-theta an eta table is evaluated at: the query's, except behind the face plane on an
/// angle-flat table, which is evaluated at its own grazing node - the same value, without a clamp
/// flag that would say something was lost.
double eta_cos_theta(const EtaTable& t, const double cos_theta) {
    if ((cos_theta < 0.0) && !t.empty() && t.is_angle_flat())
        return t.cos_thetas.front();
    return cos_theta;
}

/// Behind the face plane, whether a table that stops at its grazing node is being extrapolated -
/// which common_eval, judging coverage by the FEP table, does not charge a TOTAL table for (an
/// imported grid has side & back FEP nodes, but its total may come from a front-only donor).
bool extrapolated_behind(const EtaTable& t, const double cos_theta) {
    return (cos_theta < 0.0) && !t.empty() && !(t.cos_thetas.front() < 0.0) && !t.is_angle_flat();
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
    ec.d_near_cm = ec.d_cm;
    bool behind_back_plane = false;
    if (ec.cos_theta < 0.0) {
        // How near the source comes to the detector's outside - the nearest point of it, not the
        //  exit along the ray from the face centre, which overstates it beside a long detector.
        const ConvexSolid outside = descriptor.outer_envelope().solid();
        const double clearance = (src_cm - outside.project(src_cm)).norm();
        ec.d_near_cm = ec.a_cm + clearance;
        behind_back_plane = (!descriptor.dimensions_cm.empty())
                            && (src_cm.z() > descriptor.dimensions_cm.back());
    }
    ec.near_regime = ec.d_near_cm < floors.near_regime_a * ec.a_cm;

    // Energy validity from provenance (falls back to the table range).
    const double e_lo = provenance.valid_e_min_keV > 0.0
        ? provenance.valid_e_min_keV
        : (eta_fep.empty() ? 0.0 : eta_fep.energies_keV.front());
    const double e_hi = provenance.valid_e_max_keV > 0.0
        ? provenance.valid_e_max_keV
        : (eta_fep.empty() ? 0.0 : eta_fep.energies_keV.back());
    if (e_hi > 0.0 && (energy_keV < e_lo || energy_keV > e_hi))
        raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);

    // Behind the face plane (theta > 90 degrees).  The kernel K is ray-traced from the real
    //  position through the geometry as described, so what is in question is the residual tables
    //  and the description itself:
    //   - behind the crystal's back plane with no back structure described: the response knows
    //     nothing of what sits there - refuse-grade, as before;
    //   - angle-flat tables (curve / flat transfers) lose nothing to the clamp: the transfer
    //     envelope (saturated in sin^2 below) carries it;
    //   - tables that stop at the grazing node: extrapolated, with a measured envelope;
    //   - behind a described back: what approximating that back costs.
    SigmaBudget& b = ec.budget;
    bool behind_uncovered = false;
    if (ec.cos_theta < 0.0) {
        const bool has_back = descriptor.back_structure_extent_cm() > 0.0;
        const bool covered = !eta_fep.empty() && (eta_fep.cos_thetas.front() < 0.0);
        behind_uncovered = !covered;
        if (behind_back_plane && !has_back) {
            raise_flag(ec.flag, ResponseFlag::NeedsMc);
            b.model[SigmaBudget::BehindPlane] = model_sigma::behind_plane;
        } else {
            if (!covered && !eta_fep.is_angle_flat()) {
                raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
                b.model[SigmaBudget::BehindExtrapolation] = model_sigma::behind_extrapolation;
            }
            if (behind_back_plane)
                b.model[SigmaBudget::BackStructureApprox] = model_sigma::back_structure_approx;
        }
    }

    // EFFTRAN transfer envelope: inflate sigma where the (angle-flat) eta of a
    // transfer response cannot carry the true eta(E,theta) / near-field
    // residual. Applied UNCONDITIONALLY (independent of grounding), and to both
    // FEP and total via this shared path.
    if (model_transfer) {
        const SigmaTransferModel::Components c = model_transfer->components(
            ec.a_cm > 0.0 ? ec.d_near_cm / ec.a_cm : 1e6, ec.cos_theta, energy_keV);
        b.model[SigmaBudget::ModelFar] = c.far_onaxis;
        b.model[SigmaBudget::ModelOff] = c.offaxis;
        b.model[SigmaBudget::ModelNear] = c.near_field;
    }

    // Collimator shadow gate (spec sec 4.5): s = transmitted/geometric.  In front of the face, for
    //  a collimator.  Behind it, the same gate on what the direct path crosses: once that blocks
    //  most of the path, photons scattered in at full energy - which the removal-only kernel does
    //  not count - are much of what arrives, and vary steeply with energy and direction.  Against
    //  MC around two HPGe and two NaI detectors (2026-11):
    //   - no side & back nodes: behind the crystal's back plane the model read ~-11% at s 0.1-0.3,
    //     ~-33% at s 1e-4..1e-2 and ~-100% below, within ~5% above 0.3; beside the crystal,
    //     through thin side walls, it stayed within ~7% down to s 0.01 (the grazing node's eta
    //     carries that scatter) - so gated behind the back plane, or behind a collimator wall;
    //   - with side & back nodes: within ~5-8% above s 0.3, but 20-70% off below it, beside the
    //     crystal as well as behind it (the nodes' scatter-dominated values interpolate badly
    //     between them) - so gated everywhere past the face plane.
    const bool shadow_gate = (ec.cos_theta >= 0.0)
        ? descriptor.collimator.has_value()
        : (!behind_uncovered || behind_back_plane || descriptor.collimator.has_value());
    if (shadow_gate && q.omega_frac_active > 0.0) {
        const double s = kernel_transmitted(energy_keV, q) / q.omega_frac_active;
        if (s < model_sigma::shadow_refuse_s) {
            raise_flag(ec.flag, ResponseFlag::NeedsMc);
            b.model[SigmaBudget::Shadow] = model_sigma::shadow_refuse;
        } else if (s < model_sigma::shadow_ramp_s) {
            raise_flag(ec.flag, ResponseFlag::Shadowed);
            b.model[SigmaBudget::Shadow] = model_sigma::shadow_ramp_max *
                (model_sigma::shadow_ramp_s - s) / model_sigma::shadow_ramp_s;
        }
    }
    return ec;
}

double DetectorResponse::kernel_transmitted(double energy_keV,
                                            const ApertureQuadrature& q) const {
    std::vector<double> mu_tots;  // by RaySegment::material_slot
    for (const MacroscopicXS& xs : slot_cross_sections(energy_keV, q, mat_to_mu_, mu_tables,
                                                       "DetectorResponse::kernel_transmitted"))
        mu_tots.push_back(xs.mu_total());
    double total = 0.0;
    for (const KernelRay& r : q.rays) {
        if (r.active_len <= 0.0f) continue;
        double tau = 0.0;
        for (const RaySegment& s : q.segs(r)) {
            if (s.is_scoring) break;
            tau += mu_tots[s.material_slot] * s.length;
        }
        total += r.omega_w * std::exp(-tau);
    }
    return total;
}

void DetectorResponse::fep_budget(double energy_keV, EvalCommon& ec,
                                  double& ln_N, double& ln_k) const {
    SigmaBudget& b = ec.budget;

    // Near field: below the measured breakpoint apply N (general/contact profiles) or flag as
    // unmodeled (far-field profile).  The gate is the near-field model's own breakpoint (energy
    // and angle dependent) or the provenance floor - NOT the floor regime `near_regime`.
    ln_N = 0.0;
    // Behind the face plane the near-field table has no nodes: the boost is the kernel's alone,
    //  and the gate measures how near the source is to the detector's outside (d_near_cm).
    const bool behind = (ec.cos_theta < 0.0);
    const double d_break = near_field.breakpoint_d_cm(energy_keV, ec.cos_theta);
    double d_gate = std::max(d_break, provenance.min_distance_cm);
    // With no near-field table the gate would be provenance.min_distance_cm, which the
    // generator sets to 2a - but the boost the response cannot represent is still
    // 1.75-2.07% out at 2.4-3a (measured on the goldens' ln_n), a regime where this term
    // switched off and nothing else covered it.  Widen it to the floor's own near regime.
    if (near_field.empty())
        d_gate = std::max(d_gate, floors.near_regime_a * ec.a_cm);
    if (ec.d_near_cm < d_gate) {
        if (!near_field.empty() && !behind) {
            ln_N = near_field.ln_boost(energy_keV, ec.cos_theta, ec.d_cm);
            const double nf_sig =
                near_field.node_frac_sigma(energy_keV, ec.cos_theta, ec.d_cm);
            b.node2 += nf_sig * nf_sig;
        } else {
            raise_flag(ec.flag, ResponseFlag::NearFieldUnmodeled);
            // Charge this ONCE.  A curve transfer has no near-field table AND carries a
            // model_transfer, so it used to pay near_unmodeled here and the transfer's
            // near ramp in common_eval for the same physics - the ray-traced kernel
            // cannot know the near-field boost, which is one limitation, not two.  The
            // 2026-09 corpus measured both at ~4% and they are the same measurement.
            // Where there is no model_transfer to carry it, this term still applies.
            if (!model_transfer)
                b.model[SigmaBudget::NearUnmodeled] = model_sigma::near_unmodeled;
        }
    }

    // Grounding k(E): its fit covariance (data) and its transfer envelope (model).
    ln_k = 0.0;
    if (!grounding.empty()) {
        bool k_clamped = false;
        ln_k = grounding.eval_ln_k(energy_keV, k_clamped);
        if (k_clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);
        b.var_ln_k = grounding.var_ln_k(energy_keV);
        const SigmaTransferModel::Components c = grounding.transfer.components(
            ec.a_cm > 0.0 ? ec.d_near_cm / ec.a_cm : 1e6, ec.cos_theta, energy_keV);
        b.model[SigmaBudget::GroundFar] = c.far_onaxis;
        b.model[SigmaBudget::GroundOff] = c.offaxis;
        b.model[SigmaBudget::GroundNear] = c.near_field;
    }

    const double node_sig =
        eta_fep.node_frac_sigma(energy_keV, ec.cos_theta, ec.phi_deg);
    b.node2 += node_sig * node_sig;
    b.model[SigmaBudget::Floor] = ec.near_regime ? floors.fep_near : floors.fep_far;
}

double DetectorResponse::total_near_ln(double energy_keV, EvalCommon& ec) const {
    // Gated on the table's own breakpoint, as the FEP's is; the outermost distance node is an
    // ln N = 0 anchor, so the correction fades to nothing there.
    const NearFieldModel& nf = tot_eff.near_field;
    if (nf.empty() || (ec.cos_theta < 0.0)
        || (ec.d_cm >= nf.breakpoint_d_cm(energy_keV, ec.cos_theta)))
        return 0.0;
    const double sig = nf.node_frac_sigma(energy_keV, ec.cos_theta, ec.d_cm);
    ec.budget.node2 += sig * sig;
    return nf.ln_boost(energy_keV, ec.cos_theta, ec.d_cm);
}

EffResult DetectorResponse::fep_prefactor(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q) const {
    EvalCommon ec = common_eval(energy_keV, src_cm, q);

    bool clamped = false;
    const double ln_eta =
        eta_fep.eval_ln(energy_keV, eta_cos_theta(eta_fep, ec.cos_theta), ec.phi_deg, clamped);
    if (clamped) raise_flag(ec.flag, ResponseFlag::OutOfRangeClamped);

    double ln_N = 0.0, ln_k = 0.0;
    fep_budget(energy_keV, ec, ln_N, ln_k);

    EffResult res;
    res.value = std::exp(ln_eta + ln_N + ln_k);
    const double model2 = ec.budget.model2();
    res.sigma = res.value * std::sqrt(ec.budget.data2() + model2);
    res.sigma_model = res.value * std::sqrt(model2);
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

    // Mirrors eps_total_impl's refusal: no total data, no prefactor.
    if (!tot_eff.characterized()) {
        EffResult res;
        res.value = 0.0;
        res.sigma = 0.0;
        res.flag = ResponseFlag::NeedsMc;
        return res;
    }

    double value = 1.0, node_sig = 0.0;
    switch (tot_eff.tier) {
        case TotEffTier::NotCharacterized:
            break;                                   //handled above
        case TotEffTier::KernelExact:
            break;                                   //bare kernel; multiplier is 1
        case TotEffTier::BCurve:
            value = std::exp(tot_eff.ln_b_at(energy_keV));
            break;
        case TotEffTier::EtaTotTable: {
            bool clamped = false;
            const double ln_eta = tot_eff.eta_tot.eval_ln(
                energy_keV, eta_cos_theta(tot_eff.eta_tot, ec.cos_theta), ec.phi_deg, clamped);
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
            if (extrapolated_behind(tot_eff.eta_tot, ec.cos_theta))
                ec.budget.model[SigmaBudget::BehindExtrapolation] = model_sigma::behind_extrapolation;
            break;
        }
    }

    EffResult res;
    ec.budget.node2 = node_sig * node_sig;
    res.value = value * std::exp(total_near_ln(energy_keV, ec));
    ec.budget.model[SigmaBudget::Floor] = ec.near_regime ? floors.tot_near : floors.tot_far;
    const double model2 = ec.budget.model2();
    res.sigma = res.value * std::sqrt(ec.budget.data2() + model2);
    res.sigma_model = res.value * std::sqrt(model2);
    res.flag = ec.flag;
    return res;
}

EffResult DetectorResponse::eps_fep_impl(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q,
    const std::function<double(const Eigen::Vector3d&)>* t_src) const {
    // One body: everything but the kernel is the prefactor, so the near-field gate, the
    // grounding, the sigma budget and the flag cannot drift between the two entry points.
    EffResult res = fep_prefactor(energy_keV, src_cm, q);
    const double K = kernel_K(energy_keV, q, MuChoice::Total, t_src);
    res.value *= K;
    res.sigma *= K;
    res.sigma_model *= K;
    return res;
}

EffResult DetectorResponse::eps_total_impl(
    double energy_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q,
    const std::function<double(const Eigen::Vector3d&)>* t_src,
    const ShieldContext* sc) const {
    EvalCommon ec = common_eval(energy_keV, src_cm, q);

    // No total-efficiency data: refuse rather than serve the bare kernel as if
    // it were a verified total (see TotEffTier::NotCharacterized).
    if (!tot_eff.characterized()) {
        EffResult res;
        res.value = 0.0;
        res.sigma = 0.0;
        res.flag = ResponseFlag::NeedsMc;
        return res;
    }

    double value = 0.0;
    double node_sig = 0.0;
    switch (tot_eff.tier) {
        case TotEffTier::NotCharacterized:
            break;                                   //handled above
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
                energy_keV, eta_cos_theta(tot_eff.eta_tot, ec.cos_theta), ec.phi_deg, clamped);
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
            if (extrapolated_behind(tot_eff.eta_tot, ec.cos_theta))
                ec.budget.model[SigmaBudget::BehindExtrapolation] = model_sigma::behind_extrapolation;
            break;
        }
    }

    ec.budget.node2 = node_sig * node_sig;
    value *= std::exp(total_near_ln(energy_keV, ec));

    EffResult res;
    res.value = value;
    ec.budget.model[SigmaBudget::Floor] = ec.near_regime ? floors.tot_near : floors.tot_far;

    // Build-up seam (Stage E3 A1): only when a ShieldContext is supplied AND a
    // model is installed.  Ratio-only (>= 1), scales eps_total up; the build-up
    // sigma floor is added in quadrature.  Null sc => this whole block is
    // skipped, so the value/sigma are byte-identical to the pre-seam path.
    if (sc && buildup_model) {
        const double B = std::max(1.0, buildup_model(energy_keV, *sc));
        res.value = value * B;
        ec.budget.model[SigmaBudget::Buildup] = model_sigma::buildup_floor;
    }
    const double model2 = ec.budget.model2();
    res.sigma = res.value * std::sqrt(ec.budget.data2() + model2);
    res.sigma_model = res.value * std::sqrt(model2);
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
    const std::vector<double>& energies_keV, const Eigen::Vector3d& src_cm,
    const ApertureQuadrature& q, std::vector<double>* model_part) const {
    const size_t n = energies_keV.size();

    // The same budget eps_fep_at would report at each energy, through the same code.
    std::vector<SigmaBudget> b(n);
    for (size_t i = 0; i < n; ++i) {
        EvalCommon ec = common_eval(energies_keV[i], src_cm, q);
        double ln_N = 0.0, ln_k = 0.0;
        fep_budget(energies_keV[i], ec, ln_N, ln_k);
        b[i] = ec.budget;
    }

    std::vector<double> C(n * n, 0.0);
    if (model_part) model_part->assign(n * n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            // Every envelope is one rank-one block: m_t[i] * m_t[j].
            double m = 0.0;
            for (size_t t = 0; t < SigmaBudget::NumTerms; ++t)
                m += b[i].model[t] * b[j].model[t];
            double v = m;
            if (i == j)
                v = b[i].frac2();   // the diagonal IS the per-query budget
            else if (!grounding.empty())
                v += grounding.cov_ln_k(energies_keV[i], energies_keV[j]);
            C[i * n + j] = C[j * n + i] = v;
            if (model_part) (*model_part)[i * n + j] = (*model_part)[j * n + i] = m;
        }
    }
    return C;
}

std::vector<double> DetectorResponse::frac_covariance(
    const std::vector<double>& energies_keV, const Eigen::Vector3d& src_cm,
    std::vector<double>* model_part) const {
    const ApertureQuadrature q = make_quadrature(src_cm);
    return frac_covariance(energies_keV, src_cm, q, model_part);
}

std::vector<double> DetectorResponse::frac_covariance(
    const std::vector<double>& energies_keV, double theta_rad, double phi_rad,
    double dist_cm, std::vector<double>* model_part) const {
    return frac_covariance(energies_keV, query_position(theta_rad, phi_rad, dist_cm),
                           model_part);
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

// The FEP's near-field table and the total's share one layout.
void near_field_to_xml(XmlDoc& doc, XmlNode* parent, const NearFieldModel& m) {
    XmlNode* nf = append_node(doc, parent, "NearField");
    append_value_node(doc, nf, "Energies", join_doubles(m.energies_keV));
    append_value_node(doc, nf, "CosThetas", join_doubles(m.cos_thetas));
    append_value_node(doc, nf, "DistsCm", join_doubles(m.dists_cm));
    // ln N flattened [e][c][d] energy-major: index = (e*nc + c)*nd + d.
    append_value_node(doc, nf, "LnN", join_doubles(m.ln_n));
    append_value_node(doc, nf, "FracSigma", join_doubles(m.frac_sigma));
    append_value_node(doc, nf, "BreakCosThetas", join_doubles(m.break_cos_thetas));
    append_value_node(doc, nf, "BreakDistCm", join_doubles(m.break_d_cm));
}

void near_field_from_xml(const XmlNode* nf, NearFieldModel& m) {
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

// The content_hash counterparts of eta_to_xml, near_field_to_xml and write_detector_node: they hash
// exactly what those write, under the same conditions (see DetectorResponse::content_hash).
void hash_eta(ContentHasher& h, const EtaTable& t) {
    h.numbers(t.energies_keV);
    h.numbers(t.cos_thetas);
    if (!t.phis_deg.empty()) {
        h.tag("Phis");
        h.numbers(t.phis_deg);
    }
    if (!t.edges_keV.empty()) {
        h.tag("Edges");
        h.numbers(t.edges_keV);
    }
    h.numbers(t.ln_eta);
    h.numbers(t.frac_sigma);
}

void hash_near_field(ContentHasher& h, const NearFieldModel& m) {
    h.numbers(m.energies_keV);
    h.numbers(m.cos_thetas);
    h.numbers(m.dists_cm);
    h.numbers(m.ln_n);
    h.numbers(m.frac_sigma);
    h.numbers(m.break_cos_thetas);
    h.numbers(m.break_d_cm);
}

void hash_detector(ContentHasher& h, const GeometryDescriptor& descriptor) {
    h.flag(descriptor.shape == DetectorShape::Cylinder);
    h.integer(descriptor.crystal_material_index);
    h.flag(descriptor.reference_point == ReferencePoint::EndcapFront);
    h.flag(descriptor.symmetry == ResponseSymmetry::Quadrant);
    h.numbers(descriptor.dimensions_cm);

    if (descriptor.bullet_radius_cm > 0.0) {
        h.tag("bulletRadius");
        h.number(descriptor.bullet_radius_cm);
    }

    if (descriptor.bore) {
        h.tag("Bore");
        h.number(descriptor.bore->radius);
        h.number(descriptor.bore->depth);
        if (descriptor.bore->rounded_tip) h.tag("roundedTip");
    }

    if (descriptor.dead_layer) {
        h.tag("DeadLayer");
        h.number(descriptor.dead_layer->front);
        h.number(descriptor.dead_layer->side);
        h.number(descriptor.dead_layer->back);
    }

    h.integer(static_cast<int64_t>(descriptor.layers.size()));
    for (const LayerSpec& l : descriptor.layers) {
        h.integer(l.material_index);
        h.number(l.front_thickness_cm);
        h.number(l.side_thickness_cm);
        h.number(l.z_start_cm);
        h.number(l.z_end_cm);
    }

    if (descriptor.collimator) {
        h.tag("Collimator");
        h.integer(descriptor.collimator->material_index);
        h.number(descriptor.collimator->side_thickness_cm);
        h.number(descriptor.collimator->z_start_cm);
        h.number(descriptor.collimator->z_end_cm);
    }

    h.integer(static_cast<int64_t>(descriptor.materials.size()));
    for (const MaterialSpec& m : descriptor.materials) {
        h.text(m.name);
        h.number(m.density_g_per_cm3);
        h.integer(static_cast<int64_t>(m.composition.size()));
        for (const MaterialComponent& c : m.composition) {
            h.integer(int(c.Z));
            h.number(c.mass_fraction);
        }
    }
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
        // v3 = the eps_tot "not characterized" tier. A v2 reader maps an unknown
        // tier string onto KernelExact, which is precisely the silent misread
        // this tier exists to prevent (it would serve a bare kernel as a
        // verified total), so such a file must hard-fail on older builds rather
        // than load wrongly. Only responses that actually use the tier are
        // stamped v3, so every existing file keeps its bytes and content_hash.
        const bool needs_v3 = !tot_eff.characterized();
        // v4 = side & back coverage (see sm_xmlSerializationVersion).
        const bool needs_v4 = (!eta_fep.cos_thetas.empty() && (eta_fep.cos_thetas.front() < 0.0))
                              || (!tot_eff.eta_tot.empty() && !tot_eff.eta_tot.cos_thetas.empty()
                                  && (tot_eff.eta_tot.cos_thetas.front() < 0.0));
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%i", needs_v4 ? 4 : (needs_v3 ? 3 : (needs_v2 ? 2 : 1)));
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
        if (provenance.precision_profile)
            append_attrib(doc, p, "precisionProfile", to_string(*provenance.precision_profile));
        if (provenance.back_coverage)
            append_attrib(doc, p, "backCoverage", to_string(*provenance.back_coverage));
        if (provenance.n_anchor_angles)
            append_attrib(doc, p, "anchorAngles", std::to_string(*provenance.n_anchor_angles));
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

    if (!near_field.empty())
        near_field_to_xml(doc, root, near_field);

    {
        XmlNode* te = append_node(doc, root, "TotalEfficiency");
        const char* tier = tot_eff.tier == TotEffTier::KernelExact ? "kernel"
                           : tot_eff.tier == TotEffTier::BCurve    ? "bcurve"
                           : tot_eff.tier == TotEffTier::EtaTotTable ? "etatable"
                                                                   : "none";
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
        // Only when present, so every response without one keeps its bytes and content_hash.
        if (!tot_eff.near_field.empty())
            near_field_to_xml(doc, te, tot_eff.near_field);
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
            append_attrib(doc, tn, "offAxisS2Half", fmt_double(t.offaxis_s2_half));
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
            append_attrib(doc, tn, "offAxisS2Half", fmt_double(t.offaxis_s2_half));
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
        // Total efficiency rides in its OWN 4-column stream rather than
        // widening "Rows": a reader of an older file parses Rows at its
        // historical stride of 10 and simply finds no RowsTot, instead of
        // silently mis-striding every row.
        std::vector<double> flat_tot;
        flat_tot.reserve(c.rows.size() * 4);
        for (const AccuracyCertificate::Row& r : c.rows) {
            flat_tot.push_back(r.mc_tot);    flat_tot.push_back(r.mc_tot_sig);
            flat_tot.push_back(r.model_tot); flat_tot.push_back(r.model_tot_sig);
        }
        append_value_node(doc, cn, "RowsTot", join_doubles(flat_tot));
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
        pr.precision_profile = precision_profile_from_string(attrib_value(p, "precisionProfile"));
        pr.back_coverage = back_coverage_from_string(attrib_value(p, "backCoverage"));
        if (const char* v = attrib_value(p, "anchorAngles")) {
            char* end = nullptr;
            const long n = std::strtol(v, &end, 10);
            if ((end != v) && (*end == '\0') && (n > 0) && (n < 1000))
                pr.n_anchor_angles = static_cast<int>(n);
        }
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

    if (const XmlNode* nf = root->first_node("NearField"))
        near_field_from_xml(nf, resp->near_field);

    if (const XmlNode* te = root->first_node("TotalEfficiency")) {
        const char* tier = attrib_value(te, "tier");
        if (tier && std::strcmp(tier, "bcurve") == 0)
            resp->tot_eff.tier = TotEffTier::BCurve;
        else if (tier && std::strcmp(tier, "etatable") == 0)
            resp->tot_eff.tier = TotEffTier::EtaTotTable;
        else if (tier && std::strcmp(tier, "none") == 0)
            resp->tot_eff.tier = TotEffTier::NotCharacterized;
        else
            resp->tot_eff.tier = TotEffTier::KernelExact;
        resp->scatter_in_recapture = attrib_double(te, "scatterInRecapture", 0.0);
        resp->tot_eff.b_energies_keV = parse_doubles(child_value(te, "BEnergies"));
        resp->tot_eff.ln_b = parse_doubles(child_value(te, "LnB"));
        if (const XmlNode* et = te->first_node("EtaTot"))
            eta_from_xml(et, resp->tot_eff.eta_tot);
        if (const XmlNode* nf = te->first_node("NearField"))
            near_field_from_xml(nf, resp->tot_eff.near_field);
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
            // Absent => 0 => the legacy sin^2 form, so old files are unchanged.
            t.offaxis_s2_half = attrib_double(tn, "offAxisS2Half", 0.0);
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
            // Absent => 0 => the legacy sin^2 form, so old files are unchanged.
            t.offaxis_s2_half = attrib_double(tn, "offAxisS2Half", 0.0);
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
        // Optional; absent in certificates written before the total-efficiency
        // columns existed, which then keep their zeros.
        const std::vector<double> flat_tot = parse_doubles(child_value(cn, "RowsTot"));
        if (flat_tot.size() >= 4 * c.rows.size()) {
            for (size_t k = 0; k < c.rows.size(); ++k) {
                c.rows[k].mc_tot = flat_tot[4 * k];
                c.rows[k].mc_tot_sig = flat_tot[4 * k + 1];
                c.rows[k].model_tot = flat_tot[4 * k + 2];
                c.rows[k].model_tot_sig = flat_tot[4 * k + 3];
            }
        }
    }

    resp->finalize();
    return resp;
}

uint64_t DetectorResponse::content_hash() const {
    // Everything serialize_xml(false) stores, in its order and under its conditions - so a response
    // hashes the same as its save/load round trip - but hashed from the members: formatting the
    // ~10k doubles of a response as text cost about a millisecond, and hosts hash on every edit.
    // Left out, as there: the certificate (metadata ABOUT the response, so presence/absence never
    // changes the hash) and the version stamp (derived from the content).
    //
    // A field added to serialize_xml must be added here too.  If it is optional (written only when
    // present or non-default), hash it the same way the optional parts below are: nothing when
    // absent, `h.tag(...)` and its content when present - so every existing response keeps its hash.
    // test_detector_response_io's content_hash_tracks_every_serialized_field fails if a field that
    // serialize_xml writes is missing here (for any field one of its test responses populates).
    ContentHasher h;

    h.text(provenance.ceelo_version);
    h.text(provenance.created_utc);
    h.text(to_string(provenance.profile));
    h.number(provenance.node_fep_precision);
    if (provenance.precision_profile) {
        h.tag("precisionProfile");
        h.text(to_string(*provenance.precision_profile));
    }
    if (provenance.back_coverage) {
        h.tag("backCoverage");
        h.text(to_string(*provenance.back_coverage));
    }
    if (provenance.n_anchor_angles) {
        h.tag("anchorAngles");
        h.integer(*provenance.n_anchor_angles);
    }
    h.number(provenance.fep_window_keV);
    h.integer(static_cast<int64_t>(provenance.generation_seed));
    h.integer(provenance.kernel_n_rays);
    h.text(provenance.detector_name);
    h.number(provenance.valid_e_min_keV);
    h.number(provenance.valid_e_max_keV);
    h.number(provenance.min_distance_cm);
    h.text(to_string(provenance.method));

    hash_detector(h, descriptor);

    h.integer(static_cast<int64_t>(mu_tables.size()));
    for (const MuTable& t : mu_tables) {
        h.integer(t.material_index);
        h.numbers(t.energy_keV);
        h.numbers(t.mu_pe);
        h.numbers(t.mu_cs);
        h.numbers(t.mu_rs);
        h.numbers(t.mu_pp);
    }

    if (!eta_fep.empty()) {
        h.tag("EtaFep");
        hash_eta(h, eta_fep);
    }

    if (!near_field.empty()) {
        h.tag("NearField");
        hash_near_field(h, near_field);
    }

    h.integer(static_cast<int64_t>(tot_eff.tier));
    if (scatter_in_recapture != 0.0) {
        h.tag("scatterInRecapture");
        h.number(scatter_in_recapture);
    }
    if (tot_eff.tier == TotEffTier::BCurve) {
        h.numbers(tot_eff.b_energies_keV);
        h.numbers(tot_eff.ln_b);
    }
    if (tot_eff.tier == TotEffTier::EtaTotTable) hash_eta(h, tot_eff.eta_tot);
    if (!tot_eff.near_field.empty()) {
        h.tag("TotalNearField");
        hash_near_field(h, tot_eff.near_field);
    }

    if (!grounding.empty()) {
        h.tag("Grounding");
        h.flag(grounding.curve_derived);
        h.numbers(grounding.knot_ln_energies);
        h.numbers(grounding.ln_k);
        h.numbers(grounding.cov);
        const SigmaTransferModel& t = grounding.transfer;
        for (const double v : {t.far_onaxis, t.offaxis_mid, t.offaxis_low_e, t.low_e_ref_keV,
                               t.mid_e_ref_keV, t.near_contact, t.near_gate_a, t.offaxis_s2_half})
            h.number(v);
        h.integer(static_cast<int64_t>(grounding.points.size()));
        for (const GroundingPoint& p : grounding.points) {
            h.number(p.energy_keV);
            h.number(p.measured_eff);
            h.number(p.model_eff);
            h.number(p.frac_stat_sigma);
            h.number(p.frac_cert_sigma);
            h.text(p.source_key);
            h.number(p.distance_cm);
            h.number(p.cos_theta);
            h.number(p.phi_deg);
        }
    }

    for (const double v : {floors.fep_far, floors.fep_near, floors.tot_far, floors.tot_near,
                           floors.near_regime_a})
        h.number(v);

    if (model_transfer) {
        h.tag("ModelTransfer");
        const SigmaTransferModel& t = *model_transfer;
        for (const double v : {t.far_onaxis, t.offaxis_mid, t.offaxis_low_e, t.low_e_ref_keV,
                               t.mid_e_ref_keV, t.near_contact, t.near_gate_a, t.offaxis_s2_half})
            h.number(v);
    }

    return h.value();
}

} // namespace ceelo
