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

#include "io/EnvScatterKernel.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace ceelo {

namespace {

constexpr double kPi = 3.141592653589793;

/// Projected areas of the cylinder envelope for a parallel beam arriving at
/// incidence cosine `mu` (photon direction of travel vs +z).
struct ProjAreas {
    double face = 0.0, side = 0.0, back = 0.0;
};

ProjAreas projected_areas(double mu, const DetectorEnvelope& det)
{
    ProjAreas a;
    const double cap = kPi * det.radius_cm * det.radius_cm * std::abs(mu);
    if (mu > 0.0) a.face = cap;       // travelling +z: enters the front face
    else if (mu < 0.0) a.back = cap;  // travelling -z: enters the rear face
    a.side = 2.0 * det.radius_cm * det.length_cm *
             std::sqrt(std::max(0.0, 1.0 - mu * mu));
    return a;
}

/// Bin-AVERAGED projection weights, in closed form.
///
/// Evaluating the weights at the bin midpoint is not good enough. |mu| is linear
/// so the midpoint is exact for it, but sqrt(1-mu^2) is concave and the midpoint
/// OVERESTIMATES the side weight by 6.8% in the outermost |mu| bins (8 bins) --
/// precisely where forward-going flux concentrates, and the face-vs-side split is
/// a headline deliverable. The exact averages over [a,b] are
///     <|mu|>          = (b^2 - a^2) / (2 (b-a))            (a, b same sign)
///     <sqrt(1-mu^2)>  = [ mu sqrt(1-mu^2) + asin(mu) ]_a^b / (2 (b-a))
struct BinWeights {
    double abs_mu = 0.0;   ///< <|mu|> over the bin
    double sin_th = 0.0;   ///< <sqrt(1-mu^2)> over the bin
    double mu_mid = 0.0;   ///< midpoint, used only to pick face vs back
};

BinWeights bin_weights(const FluenceAxes& ax, int j)
{
    const double a = ax.mu_low(j), b = ax.mu_high(j);
    const double w = b - a;
    BinWeights r;
    r.mu_mid = 0.5 * (a + b);
    if (w <= 0.0) return r;

    // <|mu|>: integrate |mu| exactly, handling a bin that straddles zero (it
    // cannot with an even n_mu, but do not silently depend on that).
    auto half_sq = [](double x) { return 0.5 * x * x; };
    double int_abs;
    if (a >= 0.0)      int_abs = half_sq(b) - half_sq(a);
    else if (b <= 0.0) int_abs = half_sq(a) - half_sq(b);
    else               int_abs = half_sq(a) + half_sq(b);   // straddles 0
    r.abs_mu = int_abs / w;

    auto anti = [](double x) {
        x = std::max(-1.0, std::min(1.0, x));
        return 0.5 * (x * std::sqrt(std::max(0.0, 1.0 - x * x)) + std::asin(x));
    };
    r.sin_th = (anti(b) - anti(a)) / w;
    return r;
}

/// Bin-averaged projected areas (see bin_weights for why averaging matters).
ProjAreas binned_areas(const FluenceAxes& ax, int j, const DetectorEnvelope& det)
{
    const BinWeights w = bin_weights(ax, j);
    ProjAreas a;
    const double cap = kPi * det.radius_cm * det.radius_cm * w.abs_mu;
    if (w.mu_mid > 0.0) a.face = cap;
    else if (w.mu_mid < 0.0) a.back = cap;
    a.side = 2.0 * det.radius_cm * det.length_cm * w.sin_th;
    return a;
}

} // namespace

SurfaceFlux surface_flux(const FluenceTally& tally, const DetectorEnvelope& det)
{
    const FluenceAxes& ax = tally.axes();
    SurfaceFlux f;
    for (int ie = 0; ie < ax.n_energy; ++ie) {
        for (int j = 0; j < ax.n_mu; ++j) {
            const double phi = tally.mean(ie, j);
            if (phi == 0.0) continue;
            const ProjAreas a = binned_areas(ax, j, det);
            f.face += phi * a.face;
            f.side += phi * a.side;
            f.back += phi * a.back;
        }
    }
    // Lines are part of the incident flux too.
    for (int l = 0; l < tally.n_lines(); ++l) {
        const SurfaceFlux lf = line_surface_flux(tally, det, l);
        f.face += lf.face; f.side += lf.side; f.back += lf.back;
    }
    return f;
}

SurfaceFlux line_surface_flux(const FluenceTally& tally,
                              const DetectorEnvelope& det, int line)
{
    const FluenceAxes& ax = tally.axes();
    SurfaceFlux f;
    if (line < 0 || line >= tally.n_lines()) return f;
    for (int j = 0; j < ax.n_mu; ++j) {
        const double phi = tally.line_mean(line, j);
        if (phi == 0.0) continue;
        const ProjAreas a = binned_areas(ax, j, det);
        f.face += phi * a.face;
        f.side += phi * a.side;
        f.back += phi * a.back;
    }
    return f;
}

SurfaceFlux surface_flux_in_bin(const FluenceTally& tally,
                                const DetectorEnvelope& det, int energy_bin)
{
    const FluenceAxes& ax = tally.axes();
    SurfaceFlux f;
    if (energy_bin < 0 || energy_bin >= ax.n_energy) return f;
    for (int j = 0; j < ax.n_mu; ++j) {
        const double phi = tally.mean(energy_bin, j);
        if (phi == 0.0) continue;
        const ProjAreas a = binned_areas(ax, j, det);
        f.face += phi * a.face;
        f.side += phi * a.side;
        f.back += phi * a.back;
    }
    return f;
}

SurfaceFlux surface_flux_sigma(const FluenceTally& tally,
                               const DetectorEnvelope& det)
{
    const FluenceAxes& ax = tally.axes();
    double v_face = 0.0, v_side = 0.0, v_back = 0.0;
    for (int ie = 0; ie < ax.n_energy; ++ie) {
        for (int j = 0; j < ax.n_mu; ++j) {
            const double s = tally.sigma(ie, j);
            if (s == 0.0) continue;
            const ProjAreas a = binned_areas(ax, j, det);
            v_face += (s * a.face) * (s * a.face);
            v_side += (s * a.side) * (s * a.side);
            v_back += (s * a.back) * (s * a.back);
        }
    }
    return {std::sqrt(v_face), std::sqrt(v_side), std::sqrt(v_back)};
}

SurfaceFlux direct_surface_flux(double direct_fluence, double mu_in,
                                const DetectorEnvelope& det)
{
    const ProjAreas a = projected_areas(mu_in, det);
    return {direct_fluence * a.face, direct_fluence * a.side,
            direct_fluence * a.back};
}

// ---------------------------------------------------------------------------
// Table I/O
// ---------------------------------------------------------------------------

namespace {

void write_grid(std::ostream& os, const char* name,
                const std::vector<double>& g)
{
    os << name << " " << g.size();
    for (double v : g) os << " " << v;
    os << "\n";
}

bool read_grid(std::istream& is, std::vector<double>& g)
{
    std::size_t n = 0;
    if (!(is >> n)) return false;
    g.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        if (!(is >> g[i])) return false;
    return true;
}

} // namespace

bool save_env_kernel(const EnvKernelTable& t, const std::string& path)
{
    std::ofstream os(path);
    if (!os) return false;
    os << std::setprecision(10);

    os << "# CeeLo environment-scatter kernel\n"
       << "# format: 1\n"
       << "# kind: " << t.kind << "\n"
       << "# quantity: scattered fluence at the detector, per source photon,\n"
       << "#   NORMALIZED to the uncollided direct fluence T(E0)/(4 pi d^2).\n"
       << "#   Dimensionless. The peak transmission T is analytic, not tabulated.\n"
       << "# NOTE: no detector efficiency is folded in -- this is incident flux\n"
       << "#   only, reusable across detectors.\n"
       << "# air: " << t.air_name << " rho=" << t.air_density << " g/cm3\n";
    if (!t.slab_name.empty())
        os << "# slab: " << t.slab_name << " rho=" << t.slab_density << " g/cm3\n";
    os << "# scoring_radius_cm: " << t.scoring_radius_cm << "\n";
    if (!t.provenance.empty()) os << "# provenance: " << t.provenance << "\n";

    os << "AXES " << t.axes.e_min_keV << " " << t.axes.e_max_keV << " "
       << t.axes.n_energy << " " << t.axes.n_mu << "\n";
    write_grid(os, "E0GRID", t.e0_grid_keV);
    write_grid(os, "DGRID", t.distance_grid_cm);
    write_grid(os, "HSRCGRID", t.h_src_grid_cm);
    write_grid(os, "HDETGRID", t.h_det_grid_cm);
    write_grid(os, "EXTENTGRID", t.extent_grid_cm);

    os << "NODES " << t.nodes.size() << "\n";
    for (const EnvKernelNode& n : t.nodes) {
        os << "NODE " << n.e0_keV << " " << n.distance_cm << " " << n.h_src_cm
           << " " << n.h_det_cm << " " << n.transverse_extent_cm << " "
           << n.direct_transmission << " " << n.direct_fluence << " "
           << n.histories << " " << n.max_rel_sigma << " "
           << n.multi_slab_fraction << "\n";
        os << "SHAT " << n.s_hat.size();
        for (double v : n.s_hat) os << " " << v;
        os << "\n";
        os << "SHATSLAB " << n.s_hat_slab.size();
        for (double v : n.s_hat_slab) os << " " << v;
        os << "\n";
        os << "LINES " << n.lines.size() << "\n";
        for (const EnvKernelNode::Line& L : n.lines) {
            os << "LINE " << L.energy_keV << " " << int(L.origin) << "\n";
            os << "LSHAT " << L.s_hat.size();
            for (double v : L.s_hat) os << " " << v;
            os << "\n";
            os << "LSHATSLAB " << L.s_hat_slab.size();
            for (double v : L.s_hat_slab) os << " " << v;
            os << "\n";
        }
    }
    return static_cast<bool>(os);
}

// ---------------------------------------------------------------------------
// Binary float32 form (production container)
// ---------------------------------------------------------------------------

namespace {

constexpr char kBinMagic[4] = {'C', 'E', 'E', 'K'};
/// Binary kernel format version. Every node carries its energy edge array
/// (empty = the table-level log-uniform axis; otherwise its landmark axis).
/// Only this version is read.
constexpr uint32_t kBinVersion = 1;

void put_u32(std::ostream& os, uint32_t v) { os.write(reinterpret_cast<const char*>(&v), 4); }
void put_u64(std::ostream& os, uint64_t v) { os.write(reinterpret_cast<const char*>(&v), 8); }
void put_f64(std::ostream& os, double v) { os.write(reinterpret_cast<const char*>(&v), 8); }
void put_str(std::ostream& os, const std::string& s)
{
    put_u32(os, static_cast<uint32_t>(s.size()));
    os.write(s.data(), static_cast<std::streamsize>(s.size()));
}
void put_grid(std::ostream& os, const std::vector<double>& g)
{
    put_u32(os, static_cast<uint32_t>(g.size()));
    for (double v : g) put_f64(os, v);
}
/// Payload array as float32, truncated to `keep` leading entries.
void put_f32_array(std::ostream& os, const std::vector<double>& a,
                   std::size_t keep)
{
    keep = std::min(keep, a.size());
    put_u32(os, static_cast<uint32_t>(keep));
    put_u32(os, static_cast<uint32_t>(a.size()));
    for (std::size_t i = 0; i < keep; ++i) {
        const float f = static_cast<float>(a[i]);
        os.write(reinterpret_cast<const char*>(&f), 4);
    }
}

bool get_u32(std::istream& is, uint32_t& v) { return static_cast<bool>(is.read(reinterpret_cast<char*>(&v), 4)); }
bool get_u64(std::istream& is, uint64_t& v) { return static_cast<bool>(is.read(reinterpret_cast<char*>(&v), 8)); }
bool get_f64(std::istream& is, double& v) { return static_cast<bool>(is.read(reinterpret_cast<char*>(&v), 8)); }
bool get_str(std::istream& is, std::string& s)
{
    uint32_t n = 0;
    if (!get_u32(is, n)) return false;
    s.resize(n);
    return n == 0 || static_cast<bool>(is.read(&s[0], n));
}
bool get_grid(std::istream& is, std::vector<double>& g)
{
    uint32_t n = 0;
    if (!get_u32(is, n)) return false;
    g.resize(n);
    for (uint32_t i = 0; i < n; ++i)
        if (!get_f64(is, g[i])) return false;
    return true;
}
bool get_f32_array(std::istream& is, std::vector<double>& a)
{
    uint32_t keep = 0, full = 0;
    if (!get_u32(is, keep) || !get_u32(is, full) || keep > full) return false;
    a.assign(full, 0.0);   // truncated tail reconstitutes as exact zeros
    for (uint32_t i = 0; i < keep; ++i) {
        float f = 0.0f;
        if (!is.read(reinterpret_cast<char*>(&f), 4)) return false;
        a[i] = static_cast<double>(f);
    }
    return true;
}

/// Number of leading s_hat entries to keep for a node: everything up to and
/// including the energy bin containing E0 (row-major [ie * n_mu + imu]), IF
/// the dropped tail is actually zero -- verified, not assumed.
///
/// A landmark axis already stops at E0, so nothing is dropped there; the
/// truncation only earns its keep on a log-uniform-to-3-MeV axis.
std::size_t node_keep_entries(const FluenceAxes& ax, const EnvKernelNode& n,
                              const std::vector<double>& a)
{
    const int n_e = n.n_energy(ax);
    const int n_mu = ax.n_mu;
    if (a.size() != static_cast<std::size_t>(n_e) * n_mu) return a.size();

    int ie0 = n_e - 1;
    if (n.energy_edges.empty()) {
        const int b = ax.energy_bin(n.e0_keV);
        if (b >= 0) ie0 = b;
    }
    const std::size_t keep = static_cast<std::size_t>(ie0 + 1) * n_mu;
    for (std::size_t i = keep; i < a.size(); ++i)
        if (a[i] != 0.0) return a.size();  // invariant broken: store full
    return keep;
}

bool same_grid(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i]) return false;
    return true;
}

bool same_node_coords(const EnvKernelNode& a, const EnvKernelNode& b)
{
    return a.e0_keV == b.e0_keV && a.distance_cm == b.distance_cm &&
           a.h_src_cm == b.h_src_cm && a.h_det_cm == b.h_det_cm &&
           a.transverse_extent_cm == b.transverse_extent_cm;
}

} // namespace

bool save_env_kernel_binary(const EnvKernelTable& t, const std::string& path)
{
    std::ofstream os(path, std::ios::binary);
    if (!os) return false;

    os.write(kBinMagic, 4);
    put_u32(os, kBinVersion);
    put_str(os, t.kind);
    put_str(os, t.air_name);
    put_f64(os, t.air_density);
    put_str(os, t.slab_name);
    put_f64(os, t.slab_density);
    put_f64(os, t.scoring_radius_cm);
    put_str(os, t.provenance);
    put_f64(os, t.axes.e_min_keV);
    put_f64(os, t.axes.e_max_keV);
    put_u32(os, static_cast<uint32_t>(t.axes.n_energy));
    put_u32(os, static_cast<uint32_t>(t.axes.n_mu));
    put_grid(os, t.e0_grid_keV);
    put_grid(os, t.distance_grid_cm);
    put_grid(os, t.h_src_grid_cm);
    put_grid(os, t.h_det_grid_cm);
    put_grid(os, t.extent_grid_cm);

    put_u32(os, static_cast<uint32_t>(t.nodes.size()));
    for (const EnvKernelNode& n : t.nodes) {
        put_f64(os, n.e0_keV);
        put_f64(os, n.distance_cm);
        put_f64(os, n.h_src_cm);
        put_f64(os, n.h_det_cm);
        put_f64(os, n.transverse_extent_cm);
        put_f64(os, n.direct_transmission);
        put_f64(os, n.direct_fluence);
        put_u64(os, n.histories);
        put_f64(os, n.max_rel_sigma);
        put_f64(os, n.multi_slab_fraction);
        put_grid(os, n.energy_edges);   // per-node landmark axis (may be empty)
        put_f32_array(os, n.s_hat, node_keep_entries(t.axes, n, n.s_hat));
        put_f32_array(os, n.s_hat_slab,
                      node_keep_entries(t.axes, n, n.s_hat_slab));
        put_u32(os, static_cast<uint32_t>(n.lines.size()));
        for (const EnvKernelNode::Line& L : n.lines) {
            put_f64(os, L.energy_keV);
            put_u32(os, L.origin);
            put_f32_array(os, L.s_hat, L.s_hat.size());
            put_f32_array(os, L.s_hat_slab, L.s_hat_slab.size());
        }
    }
    return static_cast<bool>(os);
}

bool load_env_kernel_binary(EnvKernelTable& t, const std::string& path)
{
    std::ifstream is(path, std::ios::binary);
    if (!is) return false;

    char magic[4] = {};
    uint32_t version = 0;
    if (!is.read(magic, 4) || std::memcmp(magic, kBinMagic, 4) != 0)
        return false;
    if (!get_u32(is, version) || version != kBinVersion)
        return false;

    t = EnvKernelTable{};
    uint32_t u = 0;
    if (!get_str(is, t.kind) || !get_str(is, t.air_name) ||
        !get_f64(is, t.air_density) || !get_str(is, t.slab_name) ||
        !get_f64(is, t.slab_density) || !get_f64(is, t.scoring_radius_cm) ||
        !get_str(is, t.provenance) || !get_f64(is, t.axes.e_min_keV) ||
        !get_f64(is, t.axes.e_max_keV) || !get_u32(is, u))
        return false;
    t.axes.n_energy = static_cast<int>(u);
    if (!get_u32(is, u)) return false;
    t.axes.n_mu = static_cast<int>(u);
    if (!get_grid(is, t.e0_grid_keV) || !get_grid(is, t.distance_grid_cm) ||
        !get_grid(is, t.h_src_grid_cm) || !get_grid(is, t.h_det_grid_cm) ||
        !get_grid(is, t.extent_grid_cm))
        return false;

    uint32_t n_nodes = 0;
    if (!get_u32(is, n_nodes)) return false;
    t.nodes.resize(n_nodes);
    for (EnvKernelNode& n : t.nodes) {
        if (!get_f64(is, n.e0_keV) || !get_f64(is, n.distance_cm) ||
            !get_f64(is, n.h_src_cm) || !get_f64(is, n.h_det_cm) ||
            !get_f64(is, n.transverse_extent_cm) ||
            !get_f64(is, n.direct_transmission) ||
            !get_f64(is, n.direct_fluence) || !get_u64(is, n.histories) ||
            !get_f64(is, n.max_rel_sigma) ||
            !get_f64(is, n.multi_slab_fraction))
            return false;
        if (!get_grid(is, n.energy_edges)) return false;
        if (!get_f32_array(is, n.s_hat) || !get_f32_array(is, n.s_hat_slab))
            return false;
        uint32_t n_lines = 0;
        if (!get_u32(is, n_lines)) return false;
        n.lines.resize(n_lines);
        for (EnvKernelNode::Line& L : n.lines) {
            uint32_t origin = 0;
            if (!get_f64(is, L.energy_keV) || !get_u32(is, origin) ||
                !get_f32_array(is, L.s_hat) ||
                !get_f32_array(is, L.s_hat_slab))
                return false;
            L.origin = static_cast<uint8_t>(origin);
        }
    }
    return true;
}

bool merge_env_kernel_shard(EnvKernelTable& t, const EnvKernelTable& s)
{
    // An empty destination adopts the shard wholesale (first shard in).
    if (t.nodes.empty() && t.kind.empty()) {
        t = s;
        return true;
    }
    if (t.kind != s.kind || t.air_name != s.air_name ||
        t.slab_name != s.slab_name || t.air_density != s.air_density ||
        t.slab_density != s.slab_density ||
        t.scoring_radius_cm != s.scoring_radius_cm)
        return false;
    if (!same_grid(t.e0_grid_keV, s.e0_grid_keV) ||
        !same_grid(t.distance_grid_cm, s.distance_grid_cm) ||
        !same_grid(t.h_src_grid_cm, s.h_src_grid_cm) ||
        !same_grid(t.h_det_grid_cm, s.h_det_grid_cm) ||
        !same_grid(t.extent_grid_cm, s.extent_grid_cm))
        return false;
    // Log-uniform tables must share the table-level energy axis; a
    // landmark table carries its axis per node, so only n_mu has to match.
    if (t.axes.n_mu != s.axes.n_mu) return false;
    const bool landmark = !s.nodes.empty() && !s.nodes.front().energy_edges.empty();
    if (!landmark && !t.axes.same_axis(s.axes)) return false;

    // Duplicates mean overlapping shard assignments -- refuse rather than
    // average, because the duplicate implies some OTHER node was never baked.
    for (const EnvKernelNode& ns : s.nodes)
        for (const EnvKernelNode& nt : t.nodes)
            if (same_node_coords(ns, nt)) return false;

    t.nodes.insert(t.nodes.end(), s.nodes.begin(), s.nodes.end());
    return true;
}

bool env_kernel_is_complete(const EnvKernelTable& t, std::string* missing)
{
    auto axis_or_zero = [](const std::vector<double>& g) {
        return g.empty() ? std::vector<double>{0.0} : g;
    };
    const auto e0s = axis_or_zero(t.e0_grid_keV);
    const auto ds = axis_or_zero(t.distance_grid_cm);
    const auto hs = axis_or_zero(t.h_src_grid_cm);
    const auto hd = axis_or_zero(t.h_det_grid_cm);
    const auto ex = axis_or_zero(t.extent_grid_cm);

    std::ostringstream miss;
    int n_missing = 0;
    for (double e0 : e0s)
    for (double d : ds)
    for (double a : hs)
    for (double b : hd)
    for (double x : ex) {
        // Degenerate field geometries are legitimately absent: field()
        // produces no slab when |h_src - h_det| >= d, a different geometry
        // family that must not be baked into a slab table.
        if (!t.slab_name.empty() && std::abs(a - b) >= d) continue;
        bool found = false;
        for (const EnvKernelNode& n : t.nodes) {
            if (n.e0_keV == e0 && n.distance_cm == d && n.h_src_cm == a &&
                n.h_det_cm == b && n.transverse_extent_cm == x) {
                found = true;
                break;
            }
        }
        if (!found) {
            ++n_missing;
            if (n_missing <= 20)
                miss << "  E0=" << e0 << " d=" << d << " h_src=" << a
                     << " h_det=" << b << " extent=" << x << "\n";
        }
    }
    if (missing) {
        *missing = miss.str();
        if (n_missing > 20)
            *missing += "  ... and " + std::to_string(n_missing - 20) + " more\n";
    }
    return n_missing == 0;
}

bool load_env_kernel(EnvKernelTable& t, const std::string& path)
{
    // Auto-detect the binary container by its magic so consumers have one
    // load path for either representation.
    {
        std::ifstream probe(path, std::ios::binary);
        char magic[4] = {};
        if (probe && probe.read(magic, 4) &&
            std::memcmp(magic, kBinMagic, 4) == 0)
            return load_env_kernel_binary(t, path);
    }

    std::ifstream is(path);
    if (!is) return false;

    t = EnvKernelTable{};
    std::string line;
    while (is.peek() == '#') {
        std::getline(is, line);
        // Recover the identifying fields from the header comments so a loaded
        // table can be checked against what a consumer expects.
        auto after = [&](const char* key) -> std::string {
            const std::string k = std::string("# ") + key + ": ";
            return (line.rfind(k, 0) == 0) ? line.substr(k.size()) : std::string();
        };
        std::string v;
        if (!(v = after("kind")).empty()) t.kind = v;
        else if (!(v = after("provenance")).empty()) t.provenance = v;
        else if (!(v = after("format")).empty() && v != "1") {
            // Only the current text format is read.
            return false;
        }
    }

    std::string tok;
    while (is >> tok) {
        if (tok == "AXES") {
            if (!(is >> t.axes.e_min_keV >> t.axes.e_max_keV >>
                  t.axes.n_energy >> t.axes.n_mu)) return false;
        } else if (tok == "E0GRID") {
            if (!read_grid(is, t.e0_grid_keV)) return false;
        } else if (tok == "DGRID") {
            if (!read_grid(is, t.distance_grid_cm)) return false;
        } else if (tok == "HSRCGRID") {
            if (!read_grid(is, t.h_src_grid_cm)) return false;
        } else if (tok == "HDETGRID") {
            if (!read_grid(is, t.h_det_grid_cm)) return false;
        } else if (tok == "EXTENTGRID") {
            if (!read_grid(is, t.extent_grid_cm)) return false;
        } else if (tok == "NODES") {
            std::size_t n = 0;
            if (!(is >> n)) return false;
            t.nodes.reserve(n);
        } else if (tok == "NODE") {
            EnvKernelNode n;
            if (!(is >> n.e0_keV >> n.distance_cm >> n.h_src_cm >> n.h_det_cm >>
                  n.transverse_extent_cm >> n.direct_transmission >>
                  n.direct_fluence >> n.histories >> n.max_rel_sigma >>
                  n.multi_slab_fraction)) return false;
            t.nodes.push_back(std::move(n));
        } else if (tok == "LINES") {
            std::size_t n = 0;
            if (!(is >> n)) return false;
            if (!t.nodes.empty()) t.nodes.back().lines.reserve(n);
        } else if (tok == "LINE") {
            if (t.nodes.empty()) return false;
            EnvKernelNode::Line L;
            int origin = 0;
            if (!(is >> L.energy_keV >> origin)) return false;
            L.origin = static_cast<uint8_t>(origin);
            t.nodes.back().lines.push_back(std::move(L));
        } else if (tok == "LSHAT" || tok == "LSHATSLAB") {
            if (t.nodes.empty() || t.nodes.back().lines.empty()) return false;
            std::size_t m = 0;
            if (!(is >> m)) return false;
            std::vector<double>& dst = (tok == "LSHAT")
                ? t.nodes.back().lines.back().s_hat
                : t.nodes.back().lines.back().s_hat_slab;
            dst.resize(m);
            for (std::size_t i = 0; i < m; ++i)
                if (!(is >> dst[i])) return false;
        } else if (tok == "SHAT" || tok == "SHATSLAB") {
            if (t.nodes.empty()) return false;
            std::size_t m = 0;
            if (!(is >> m)) return false;
            std::vector<double>& dst =
                (tok == "SHAT") ? t.nodes.back().s_hat : t.nodes.back().s_hat_slab;
            dst.resize(m);
            for (std::size_t i = 0; i < m; ++i)
                if (!(is >> dst[i])) return false;
        }
    }
    return true;
}

} // namespace ceelo
