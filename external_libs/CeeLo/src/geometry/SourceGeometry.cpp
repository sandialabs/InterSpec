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

#include "geometry/SourceGeometry.h"
#include "geometry/Cylinder.h"
#include "geometry/Box.h"
#include "geometry/Sphere.h"
#include "physics/ElectronCsda.h"
#include "transport/ComptonScatter.h"
#include "transport/TransportUtils.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <random>
#include <utility>

namespace ceelo {

namespace {

/// Path length of a ray through a spherical shell [r_in, r_out] centered at
/// `center`. Sums the near-wall and far-wall crossings (a ray born in the shell
/// can cross the inner void and re-enter the far wall). r_in = 0 → solid ball.
/// Mirrors the inner+outer sphere logic the point-source shield loop has always
/// used, factored so the spherical *source material* and shields share it.
double spherical_shell_path(const Eigen::Vector3d& center,
                            const Eigen::Vector3d& position,
                            const Eigen::Vector3d& direction,
                            double r_in, double r_out) {
    auto outer_hit = intersect_sphere(position, direction, center, r_out);
    if (!outer_hit || !outer_hit->valid()) return 0.0;

    double t0_out = std::max(outer_hit->t_enter, 0.0);
    double t1_out = outer_hit->t_exit;

    if (r_in > 1e-10) {
        auto inner_hit = intersect_sphere(position, direction, center, r_in);
        if (inner_hit && inner_hit->valid()) {
            double path = 0.0;
            double t1_in = inner_hit->t_exit;
            if (t1_in < t1_out) path = t1_out - t1_in;
            double t0_in = std::max(inner_hit->t_enter, 0.0);
            if (t0_out < t0_in) path += t0_in - t0_out;
            return path;
        }
    }
    return t1_out - t0_out;
}

/// Path length of a ray (in the cylinder-local frame) through the region between
/// an inner and outer coaxial cylinder. Inner cylinder = (r_in, half-length
/// hz_in); outer = (r_out, hz_out). r_in = 0 → solid. Mirrors the inner+outer
/// tube logic the cylindrical shield loop uses; reused for the annular
/// (hollow-bore) source material where hz_in == hz_out.
double cyl_shell_path(const Eigen::Vector3d& local_pos,
                      const Eigen::Vector3d& local_dir,
                      double r_in, double hz_in,
                      double r_out, double hz_out) {
    auto outer_hit = intersect_cylinder(local_pos, local_dir, r_out, -hz_out, hz_out);
    if (!outer_hit || !outer_hit->valid()) return 0.0;

    double t0_out = std::max(outer_hit->t_enter, 0.0);
    double t1_out = outer_hit->t_exit;

    if (r_in > 1e-10) {
        auto inner_hit = intersect_cylinder(local_pos, local_dir, r_in, -hz_in, hz_in);
        if (inner_hit && inner_hit->valid()) {
            double path = 0.0;
            double t1_in = inner_hit->t_exit;
            if (t1_in < t1_out) path = t1_out - t1_in;
            double t0_in = std::max(inner_hit->t_enter, 0.0);
            if (t0_out < t0_in) path += t0_in - t0_out;
            return path;
        }
    }
    return t1_out - t0_out;
}

/// Path length of a ray (in the box-local frame) through the shell between an
/// inner and outer centered box. Inner all-zero → solid. Mirrors cyl_shell_path
/// (and the inner+outer box logic the rectangular shield loop uses); used for
/// the hollow-box (crate/container-wall) source material.
double box_shell_path(const Eigen::Vector3d& local_pos,
                      const Eigen::Vector3d& local_dir,
                      const Eigen::Vector3d& inner_hd,
                      const Eigen::Vector3d& outer_hd) {
    auto outer_hit = intersect_box(local_pos, local_dir,
                                   outer_hd.x(), outer_hd.y(),
                                   -outer_hd.z(), outer_hd.z());
    if (!outer_hit || !outer_hit->valid()) return 0.0;

    double t0_out = std::max(outer_hit->t_enter, 0.0);
    double t1_out = outer_hit->t_exit;

    if (inner_hd.minCoeff() > 1e-10) {
        auto inner_hit = intersect_box(local_pos, local_dir,
                                       inner_hd.x(), inner_hd.y(),
                                       -inner_hd.z(), inner_hd.z());
        if (inner_hit && inner_hit->valid()) {
            double path = 0.0;
            double t1_in = inner_hit->t_exit;
            if (t1_in < t1_out) path = t1_out - t1_in;
            double t0_in = std::max(inner_hit->t_enter, 0.0);
            if (t0_out < t0_in) path += t0_in - t0_out;
            return path;
        }
    }
    return t1_out - t0_out;
}

bool inside_surface(SourceGeometry::Shape shape, const SourceLayer& s,
                    const Eigen::Vector3d& p) {
    // Points land on a boundary constantly here (every interval endpoint is one),
    // so the test is deliberately inclusive; midpoints are what it is asked about.
    constexpr double kEps = 1e-12;
    switch (shape) {
    case SourceGeometry::Shape::Sphere:
        return p.norm() <= s.dims[0] + kEps;
    case SourceGeometry::Shape::Cylindrical:
        return std::hypot(p.x(), p.y()) <= s.dims[0] + kEps
               && std::fabs(p.z()) <= s.dims[1] + kEps;
    case SourceGeometry::Shape::Rectangular:
        return std::fabs(p.x()) <= s.dims[0] + kEps
               && std::fabs(p.y()) <= s.dims[1] + kEps
               && std::fabs(p.z()) <= s.dims[2] + kEps;
    default:
        return false;
    }
}

void surface_crossings(SourceGeometry::Shape shape, const SourceLayer& s,
                       const Eigen::Vector3d& p, const Eigen::Vector3d& d,
                       std::vector<double>& out) {
    std::optional<RayHit> hit;
    switch (shape) {
    case SourceGeometry::Shape::Sphere:
        hit = intersect_sphere(p, d, Eigen::Vector3d::Zero(), s.dims[0]);
        break;
    case SourceGeometry::Shape::Cylindrical:
        hit = intersect_cylinder(p, d, s.dims[0], -s.dims[1], s.dims[1]);
        break;
    case SourceGeometry::Shape::Rectangular:
        hit = intersect_box(p, d, s.dims[0], s.dims[1], -s.dims[2], s.dims[2]);
        break;
    default:
        return;
    }
    if (!hit || !hit->valid()) return;
    if (hit->t_enter > 0.0) out.push_back(hit->t_enter);
    if (hit->t_exit > 0.0) out.push_back(hit->t_exit);
}

/// Ordered ray march through a concentric stack given innermost-layer-first.
///
/// Unlike the per-layer shell-path helpers above - which return one AGGREGATE
/// length per layer, merging the near-wall and far-wall crossings of a hollow
/// shell into a single number - this produces segments in true traversal order.
/// That distinction only matters once the inner region attenuates: a merged
/// near+far length still gives the right exp(-mu L), but transport_source_photon()
/// locates interactions by accumulating segment lengths along the ray, so a
/// merged segment would place a far-wall interaction short by the whole core
/// chord.  Hence: cores => this walker; no cores => the original loops, untouched.
///
/// Method: collect every forward crossing of every surface, sort, and label each
/// interval by the innermost surface containing its midpoint.  O(n log n) in the
/// number of surfaces and correct for any nested convex stack, including rays
/// that miss the inner surfaces entirely.
void trace_concentric(SourceGeometry::Shape shape,
                      const std::vector<SourceLayer>& surfs,
                      const Eigen::Vector3d& local_pos,
                      const Eigen::Vector3d& local_dir,
                      std::vector<SourceGeometry::SourcePathSegment>& out,
                      std::size_t max_segments) {
    std::vector<double> ts;
    ts.reserve(2 * surfs.size() + 1);
    ts.push_back(0.0);
    for (const SourceLayer& s : surfs)
        surface_crossings(shape, s, local_pos, local_dir, ts);

    std::sort(ts.begin(), ts.end());

    for (std::size_t i = 0; i + 1 < ts.size(); ++i) {
        const double t0 = ts[i], t1 = ts[i + 1];
        const double len = t1 - t0;
        if (len <= 1e-10) continue;

        const Eigen::Vector3d mid = local_pos + local_dir * (0.5 * (t0 + t1));

        const Material* mat = nullptr;
        bool in_assembly = false;
        for (const SourceLayer& s : surfs) {
            if (inside_surface(shape, s, mid)) {
                mat = s.material;
                in_assembly = true;
                break;  //surfs is innermost-first, so the first hit is the region
            }
        }
        if (!in_assembly) break;  //left the outermost surface; nothing beyond

        // Merge with the previous segment when the material repeats, so a stack
        // whose neighbours share a material reads the way the layer loops do.
        // Merging is checked BEFORE the cap: a caller asking for one segment
        // (the Moliere walk) wants the full run of its current material, not the
        // first interval of it.
        if (!out.empty() && out.back().material == mat) {
            out.back().length += len;
            continue;
        }
        if (out.size() >= max_segments) return;
        out.push_back({mat, len});
    }
}

} // namespace

Eigen::Vector3d SourceGeometry::layer_delta(const SourceLayer& l) const {
    // The per-shape reading of a thickness triple, in the same meaning `dims`
    //  has: Sphere {radius}; Cylinder {radial, end}; Box {x, y, z}.
    switch (shape_) {
    case Shape::Cylindrical: return Eigen::Vector3d(l.tx, l.tz, 0.0);
    case Shape::Rectangular: return Eigen::Vector3d(l.tx, l.ty, l.tz);
    default:                 return Eigen::Vector3d(l.tx, 0.0, 0.0);  //sphere
    }
}

Eigen::Vector3d SourceGeometry::shape_inner_dims() const {
    switch (shape_) {
    case Shape::Sphere:      return Eigen::Vector3d(sphere_inner_r_, 0.0, 0.0);
    case Shape::Cylindrical:
        // cyl_inner_half_length_ is set to the full half-length for a SOLID
        //  cylinder (that is what "through-bore" degenerates to), so testing the
        //  triple's maxCoeff would call a solid cylinder hollow: it would gain a
        //  phantom void layer, accept add_core(), and lose the closed-form and
        //  electron-containment fast paths for good.  The cavity exists only
        //  when there is a bore.
        if (cyl_inner_r_ <= 1e-10) return Eigen::Vector3d::Zero();
        return Eigen::Vector3d(cyl_inner_r_, cyl_inner_half_length_, 0.0);
    case Shape::Rectangular: return rect_inner_half_dims_;
    default:                 return Eigen::Vector3d::Zero();  //Point, Marinelli
    }
}

Eigen::Vector3d SourceGeometry::shape_outer_dims() const {
    switch (shape_) {
    case Shape::Sphere:      return Eigen::Vector3d(sphere_radius_, 0.0, 0.0);
    case Shape::Cylindrical: return Eigen::Vector3d(cyl_radius_, cyl_half_length_, 0.0);
    case Shape::Rectangular: return rect_half_dims_;
    default:                 return Eigen::Vector3d::Zero();  //Point, Marinelli
    }
}

void SourceGeometry::rebuild_layer_dims() {
    // Split the current stack back into its three runs.  n_core_layers_ plus
    //  has_center_void_ describe the layout exactly, so this survives any
    //  ordering of configure_*/set_source_material/add_shield/add_core calls.
    const std::size_t first_core = has_center_void_ ? 1u : 0u;
    const std::size_t src_at = first_core + n_core_layers_;
    assert(src_at < layers_.size());

    std::vector<SourceLayer> cores(layers_.begin() + first_core,
                                   layers_.begin() + src_at);
    std::vector<SourceLayer> shields(layers_.begin() + src_at + 1, layers_.end());
    SourceLayer src = layers_[src_at];

    const Eigen::Vector3d inner = shape_inner_dims();
    const Eigen::Vector3d outer = shape_outer_dims();

    // Cores are stored innermost-first, but they FILL inward: the last one is
    //  the layer against the source's inner surface.  Walk them outermost-first
    //  (i.e. backwards), taking each one's thickness off the cavity as we go.
    Eigen::Vector3d d = inner;
    for (std::size_t i = cores.size(); i-- > 0; ) {
        cores[i].dims = d;
        d = (d - layer_delta(cores[i])).cwiseMax(0.0);
    }

    layers_.clear();
    layers_.reserve(cores.size() + shields.size() + 2);

    // Whatever the cores did not reach stays a genuine void.  Keeping it in the
    //  list (null material, which every consumer already skips) is what makes
    //  the along-ray distances right for a partly-filled core.
    has_center_void_ = (d.maxCoeff() > 1e-10);
    if (has_center_void_)
        layers_.push_back({nullptr, 0.0, 0.0, 0.0, d});

    for (const SourceLayer& c : cores) layers_.push_back(c);

    source_layer_index_ = layers_.size();
    src.material = source_material_;
    src.dims = outer;
    layers_.push_back(src);

    Eigen::Vector3d o = outer;
    for (SourceLayer sh : shields) {
        o += layer_delta(sh);
        sh.dims = o;
        layers_.push_back(sh);
    }

    // Derive the gate from what the stack actually holds.  It is READ once per
    //  Moliere substep so it must stay a cached bool, but deriving it here means
    //  it cannot drift out of step with the layers the way a write-once flag did.
    has_attenuating_interior_ = false;
    for (std::size_t i = 0; i < source_layer_index_; ++i) {
        if (layers_[i].material) { has_attenuating_interior_ = true; break; }
    }
}

void SourceGeometry::init_layer_stack() {
    // Seed the stack the first time a shape is configured; a re-configure keeps
    //  whatever layers were already added and just re-sizes them.
    if (layers_.empty()) {
        n_core_layers_ = 0;
        has_center_void_ = false;
        source_layer_index_ = 0;
        layers_.push_back({source_material_, 0.0, 0.0, 0.0,
                           Eigen::Vector3d::Zero()});
    }

    // Point and Marinelli have no hollow interior, so they cannot carry cores.
    //  Reconfiguring INTO one of those shapes after add_core() would otherwise
    //  leave core layers behind and the gate set, routing every trace into a
    //  walker that only handles the concentric shapes - an assert in a debug
    //  build, and in release a silent "no segments", which dropped the shields.
    if ((shape_ == Shape::Point || shape_ == Shape::Marinelli)
        && (n_core_layers_ > 0 || has_center_void_)) {
        const std::size_t first = has_center_void_ ? 1u : 0u;
        layers_.erase(layers_.begin(), layers_.begin() + first + n_core_layers_);
        n_core_layers_ = 0;
        has_center_void_ = false;
        source_layer_index_ = 0;
    }
    rebuild_layer_dims();
}

// A shield grows the stack outward; a core fills the cavity inward.  Which side
//  a layer lands on is the ONLY difference between the two.
void SourceGeometry::add_shield(const Material* mat, double thickness) {
    assert(mat && thickness > 0.0);
    push_shield_layer({mat, thickness, thickness, thickness, {}});
}

void SourceGeometry::add_shield(const Material* mat, double t_radial, double t_end) {
    assert(mat && t_radial >= 0.0 && t_end >= 0.0 && (t_radial > 0.0 || t_end > 0.0));
    assert(!configured_ || shape_ == Shape::Cylindrical);
    push_shield_layer({mat, t_radial, t_radial, t_end, {}});
}

void SourceGeometry::add_shield(const Material* mat, double t_x, double t_y, double t_z) {
    assert(mat && t_x >= 0.0 && t_y >= 0.0 && t_z >= 0.0
           && (t_x > 0.0 || t_y > 0.0 || t_z > 0.0));
    assert(!configured_ || shape_ == Shape::Rectangular);
    push_shield_layer({mat, t_x, t_y, t_z, {}});
}

void SourceGeometry::add_core(const Material* mat, double thickness) {
    assert(mat && thickness > 0.0);
    push_core_layer({mat, thickness, thickness, thickness, {}});
}

void SourceGeometry::add_core(const Material* mat, double t_radial, double t_end) {
    assert(mat && t_radial >= 0.0 && t_end >= 0.0 && (t_radial > 0.0 || t_end > 0.0));
    assert(!configured_ || shape_ == Shape::Cylindrical);
    push_core_layer({mat, t_radial, t_radial, t_end, {}});
}

void SourceGeometry::add_core(const Material* mat, double t_x, double t_y, double t_z) {
    assert(mat && t_x >= 0.0 && t_y >= 0.0 && t_z >= 0.0
           && (t_x > 0.0 || t_y > 0.0 || t_z > 0.0));
    assert(!configured_ || shape_ == Shape::Rectangular);
    push_core_layer({mat, t_x, t_y, t_z, {}});
}

void SourceGeometry::push_shield_layer(SourceLayer layer) {
    if (layers_.empty()) {
        // add_shield() before configure_*(): seed the emitting layer so the
        //  shells have something to accumulate from.  A point source's shells
        //  legitimately start at r = 0, which is what a zero-dims seed gives.
        n_core_layers_ = 0;
        has_center_void_ = false;
        source_layer_index_ = 0;
        layers_.push_back({source_material_, 0.0, 0.0, 0.0,
                           Eigen::Vector3d::Zero()});
    }
    layers_.push_back(layer);
    rebuild_layer_dims();
}

void SourceGeometry::push_core_layer(SourceLayer layer) {
    assert(!layers_.empty() && "add_core() requires a configured source");
    assert(shape_ == Shape::Sphere || shape_ == Shape::Cylindrical
           || shape_ == Shape::Rectangular);
    assert(shape_inner_dims().maxCoeff() > 1e-10
           && "add_core() requires a hollow source");

    // Cores are given OUTERMOST FIRST, so each new one goes immediately inside
    //  the previous one - i.e. at the innermost end of the core run, which is
    //  where the stack stores it.  rebuild_layer_dims() then re-derives every
    //  boundary from the cavity inward.
    const std::size_t first_core = has_center_void_ ? 1u : 0u;
    layers_.insert(layers_.begin() + first_core, layer);
    ++n_core_layers_;
    rebuild_layer_dims();   // re-derives has_attenuating_interior_
}

void SourceGeometry::configure_point(const Eigen::Vector3d& position) {
    shape_ = Shape::Point;
    point_position_ = position;
    // A point source has no volume: the emitting layer is a zero-dims seed and
    //  the shells accumulate outward from r = 0, exactly as before.
    init_layer_stack();
    configured_ = true;
}

void SourceGeometry::configure_cylindrical(const Eigen::Vector3d& center, double radius,
                                            double half_length, const Eigen::Matrix3d& rotation,
                                            double inner_radius,
                                            double inner_half_length) {
    assert(inner_radius >= 0.0 && inner_radius < radius);
    shape_ = Shape::Cylindrical;
    cyl_center_ = center;
    cyl_radius_ = radius;
    cyl_inner_r_ = inner_radius;
    // Negative = through-bore (a pipe): the hollow region runs the full length.
    cyl_inner_half_length_ = (inner_half_length < 0.0) ? half_length : inner_half_length;
    assert(cyl_inner_half_length_ >= 0.0 && cyl_inner_half_length_ <= half_length);
    cyl_half_length_ = half_length;
    cyl_rotation_ = rotation;
    init_layer_stack();
    configured_ = true;
}

void SourceGeometry::configure_spherical(const Eigen::Vector3d& center, double outer_radius,
                                          double inner_radius, const Eigen::Matrix3d& rotation) {
    assert(inner_radius >= 0.0 && inner_radius < outer_radius);
    shape_ = Shape::Sphere;
    sphere_center_ = center;
    sphere_inner_r_ = inner_radius;
    sphere_radius_ = outer_radius;
    sphere_rotation_ = rotation;
    init_layer_stack();
    configured_ = true;
}

void SourceGeometry::configure_rectangular(const Eigen::Vector3d& center,
                                            const Eigen::Vector3d& half_dims,
                                            const Eigen::Matrix3d& rotation,
                                            const Eigen::Vector3d& inner_half_dims) {
    // Inner void must be all-zero (solid) or strictly inside the outer box.
    assert((inner_half_dims.array() == 0.0).all()
           || ((inner_half_dims.array() >= 0.0).all()
               && (inner_half_dims.array() < half_dims.array()).all()));
    shape_ = Shape::Rectangular;
    rect_center_ = center;
    rect_half_dims_ = half_dims;
    rect_inner_half_dims_ = inner_half_dims;
    rect_rotation_ = rotation;
    init_layer_stack();
    configured_ = true;
}

void SourceGeometry::configure_marinelli(double well_inner_radius, double outer_radius,
                                          double z_bk, double z_we, double z_bot) {
    shape_ = Shape::Marinelli;
    marinelli_well_r_ = well_inner_radius;
    marinelli_outer_r_ = outer_radius;
    marinelli_z_bk_ = z_bk;
    marinelli_z_we_ = z_we;
    marinelli_z_bot_ = z_bot;
    // Marinelli is an L-shaped fill on absolute z-planes, not a concentric
    //  stack, so its layers carry thicknesses only and `dims` stays zero; the
    //  beaker-wall paths read scalar_thickness() and the absolute z-planes.
    init_layer_stack();
    configured_ = true;
}

void SourceGeometry::set_detector_bounds(double det_radius, double det_z_min, double det_z_max) {
    det_bound_radius_ = det_radius;
    det_bound_z_min_ = det_z_min;
    det_bound_z_max_ = det_z_max;
}

bool SourceGeometry::is_inside_marinelli(const Eigen::Vector3d& pos) const {
    double r2 = pos.x() * pos.x() + pos.y() * pos.y();
    double z = pos.z();
    double or2 = marinelli_outer_r_ * marinelli_outer_r_;
    double wr2 = marinelli_well_r_ * marinelli_well_r_;

    if (r2 > or2) return false;
    if (z < marinelli_z_bot_ || z > marinelli_z_we_) return false;

    // In the disk region (z < z_bk): all r ∈ [0, outer_r] are valid
    if (z <= marinelli_z_bk_) return true;

    // In the ring region (z > z_bk): only r ∈ [well_r, outer_r]
    return r2 >= wr2;
}

SourceGeometry::MarinelliExitInfo SourceGeometry::marinelli_exit_info(
    const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) const
{
    // Find nearest exit distance from the L-shaped Marinelli volume,
    // tracking which boundary surface is hit.
    // Boundaries:
    //   1. Outer cylinder: r = outer_r, z ∈ [z_bot, z_we]
    //   2. Well inner wall: r = well_r, z ∈ [z_bk, z_we]
    //   3. Bottom plane: z = z_bot
    //   4. Top of ring: z = z_we (r >= well_r)
    //   5. Well opening: z = z_bk, r < well_r (open — no beaker wall)

    constexpr double kEps = 1e-10;
    double t_min = 1e30;
    MarinelliExitSurface surface = MarinelliExitSurface::None;

    double px = pos.x(), py = pos.y(), pz = pos.z();
    double dx = dir.x(), dy = dir.y(), dz = dir.z();
    double r2 = px * px + py * py;

    // --- Outer cylinder (exit from inside) ---
    {
        double a = dx * dx + dy * dy;
        double b = 2.0 * (px * dx + py * dy);
        double c = r2 - marinelli_outer_r_ * marinelli_outer_r_;
        double disc = b * b - 4.0 * a * c;
        if (disc > 0.0 && a > kEps) {
            double sqrt_disc = std::sqrt(disc);
            double t1 = (-b + sqrt_disc) / (2.0 * a);
            if (t1 > kEps) {
                double z_hit = pz + t1 * dz;
                if (z_hit >= marinelli_z_bot_ - kEps && z_hit <= marinelli_z_we_ + kEps) {
                    if (t1 < t_min) {
                        t_min = t1;
                        surface = MarinelliExitSurface::OuterCylinder;
                    }
                }
            }
        }
    }

    // --- Well inner wall (r = well_r, z ∈ [z_bk, z_we]) ---
    {
        double a = dx * dx + dy * dy;
        double b = 2.0 * (px * dx + py * dy);
        double c = r2 - marinelli_well_r_ * marinelli_well_r_;
        double disc = b * b - 4.0 * a * c;
        if (disc > 0.0 && a > kEps) {
            double sqrt_disc = std::sqrt(disc);
            double t1 = (-b - sqrt_disc) / (2.0 * a);
            double t2 = (-b + sqrt_disc) / (2.0 * a);

            for (double t : {t1, t2}) {
                if (t > kEps) {
                    double z_hit = pz + t * dz;
                    if (z_hit >= marinelli_z_bk_ - kEps && z_hit <= marinelli_z_we_ + kEps) {
                        if (t < t_min) {
                            t_min = t;
                            surface = MarinelliExitSurface::WellWall;
                        }
                        break;
                    }
                }
            }
        }
    }

    // --- Bottom plane: z = z_bot ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_bot_ - pz) / dz;
        if (t > kEps) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            if (r2_hit <= marinelli_outer_r_ * marinelli_outer_r_ + kEps) {
                if (t < t_min) {
                    t_min = t;
                    surface = MarinelliExitSurface::Bottom;
                }
            }
        }
    }

    // --- Top plane: z = z_we (ring region, r >= well_r) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_we_ - pz) / dz;
        if (t > kEps) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            double wr2 = marinelli_well_r_ * marinelli_well_r_;
            if (r2_hit >= wr2 - kEps && r2_hit <= marinelli_outer_r_ * marinelli_outer_r_ + kEps) {
                if (t < t_min) {
                    t_min = t;
                    surface = MarinelliExitSurface::Top;
                }
            }
        }
    }

    // --- Well opening: z = z_bk, r < well_r (open face toward detector) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_bk_ - pz) / dz;
        if (t > kEps) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            double wr2 = marinelli_well_r_ * marinelli_well_r_;
            if (r2_hit < wr2 + kEps) {
                if (t < t_min) {
                    t_min = t;
                    surface = MarinelliExitSurface::WellOpening;
                }
            }
        }
    }

    if (t_min < 1e29) {
        return {t_min, surface};
    }
    return {0.0, MarinelliExitSurface::None};
}

double SourceGeometry::marinelli_boundary_distance(const Eigen::Vector3d& pos,
                                                    const Eigen::Vector3d& dir) const {
    return marinelli_exit_info(pos, dir).distance;
}

double SourceGeometry::marinelli_wall_factor(const Eigen::Vector3d& position,
                                              const Eigen::Vector3d& direction) const {
    if (!is_inside_marinelli(position)) return 0.0;

    auto exit_info = marinelli_exit_info(position, direction);

    // Well opening has no beaker wall
    if (exit_info.surface == MarinelliExitSurface::WellOpening ||
        exit_info.surface == MarinelliExitSurface::None) {
        return 0.0;
    }

    // Compute cosine of incidence angle at exit surface
    double cos_inc = 0.0;
    Eigen::Vector3d exit_pos = position + direction * exit_info.distance;

    switch (exit_info.surface) {
    case MarinelliExitSurface::OuterCylinder:
    case MarinelliExitSurface::WellWall: {
        // Cylindrical surface: normal is radial in xy-plane
        double ex = exit_pos.x(), ey = exit_pos.y();
        double r_exit = std::sqrt(ex * ex + ey * ey);
        if (r_exit > 1e-10) {
            cos_inc = std::abs(direction.x() * ex / r_exit +
                               direction.y() * ey / r_exit);
        }
        break;
    }
    case MarinelliExitSurface::Bottom:
    case MarinelliExitSurface::Top:
        // Flat surface: normal is along z
        cos_inc = std::abs(direction.z());
        break;
    default:
        break;
    }

    // Cap to prevent extreme paths at grazing incidence (max 20x thickness)
    cos_inc = std::max(cos_inc, 0.05);

    return 1.0 / cos_inc;
}

SourceGeometry::MarinelliReentryInfo SourceGeometry::compute_marinelli_reentry(
    const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) const
{
    // Check if a photon from OUTSIDE the Marinelli water can re-enter.
    // Three entry paths:
    //   1. Through well wall → ring water (outward from well cavity)
    //   2. Through well opening → disk water (downward, no wall)
    //   3. Through ring top → ring water (downward from above z_we, no wall)

    constexpr double kEps = 1e-9;
    double t_best = 1e30;
    double wall_path_best = 0.0;

    double px = pos.x(), py = pos.y(), pz = pos.z();
    double dx = dir.x(), dy = dir.y(), dz = dir.z();
    double or2 = marinelli_outer_r_ * marinelli_outer_r_;

    // Beaker wall thickness (from shields — first layer)
    double wall_t = 0.0;
    if (!shields().empty()) wall_t = shields()[0].scalar_thickness();

    // --- 1. Ring entry through well wall (r = marinelli_well_r_) ---
    // Photon crosses the well outer surface (= water inner boundary) at r = well_r
    {
        double a = dx * dx + dy * dy;
        double b = 2.0 * (px * dx + py * dy);
        double r2_cur = px * px + py * py;
        double c = r2_cur - marinelli_well_r_ * marinelli_well_r_;
        double disc = b * b - 4.0 * a * c;
        if (disc > 0.0 && a > kEps) {
            double sqrt_disc = std::sqrt(disc);
            // We want the FIRST positive intersection (entering the cylinder from inside)
            double t1 = (-b - sqrt_disc) / (2.0 * a);
            double t2 = (-b + sqrt_disc) / (2.0 * a);
            // From inside the well (r < well_r), the exit is the + root
            for (double t : {t1, t2}) {
                if (t > kEps) {
                    double z_hit = pz + t * dz;
                    if (z_hit >= marinelli_z_bk_ - kEps && z_hit <= marinelli_z_we_ + kEps) {
                        // Check that the hit point is within the ring (r_hit ≈ well_r, which is ≥ well_r)
                        if (t < t_best) {
                            t_best = t;
                            // Compute wall path: distance through the cylindrical shell
                            // from r = well_r - wall_t to r = well_r
                            if (wall_t > 0.0) {
                                double r_inner = marinelli_well_r_ - wall_t;
                                double c_inner = r2_cur - r_inner * r_inner;
                                double disc_inner = b * b - 4.0 * a * c_inner;
                                if (disc_inner > 0.0) {
                                    double t_inner = (-b + std::sqrt(disc_inner)) / (2.0 * a);
                                    if (t_inner > kEps && t_inner < t) {
                                        wall_path_best = t - t_inner;
                                    } else {
                                        // Already inside wall
                                        wall_path_best = t;
                                    }
                                } else {
                                    wall_path_best = wall_t; // fallback
                                }
                            } else {
                                wall_path_best = 0.0;
                            }
                        }
                        break;
                    }
                }
            }
        }
    }

    // --- 2. Disk entry through well opening (z = z_bk, r < well_r, no wall) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_bk_ - pz) / dz;
        if (t > kEps && t < t_best) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            // Well opening: r < well_r enters disk water (no beaker wall)
            if (r2_hit < marinelli_well_r_ * marinelli_well_r_) {
                // Also verify it's within the outer radius
                if (r2_hit <= or2 + kEps) {
                    t_best = t;
                    wall_path_best = 0.0; // No wall at well opening
                }
            }
        }
    }

    // --- 3. Ring top entry (z = z_we, r ∈ [well_r, outer_r], no wall) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_we_ - pz) / dz;
        if (t > kEps && t < t_best) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            double wr2 = marinelli_well_r_ * marinelli_well_r_;
            if (r2_hit >= wr2 - kEps && r2_hit <= or2 + kEps) {
                t_best = t;
                wall_path_best = 0.0; // No wall at ring top
            }
        }
    }

    if (t_best < 1e29) {
        // Advance slightly past the boundary into the water
        Eigen::Vector3d entry = pos + dir * (t_best + 1e-8);
        return {true, t_best, wall_path_best, entry};
    }

    return {false, 0.0, 0.0, {0, 0, 0}};
}

double SourceGeometry::source_material_path(const Eigen::Vector3d& position,
                                             const Eigen::Vector3d& direction) const {
    if (!source_material_) return 0.0;

    if (shape_ == Shape::Point) {
        return 0.0; // Point sources have no volume
    }

    if (shape_ == Shape::Cylindrical) {
        // Transform to source local frame
        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);
        Eigen::Vector3d local_dir = cyl_rotation_ * direction;

        // Annulus when a bore is present (subtract the inactive inner region);
        // cyl_inner_r_ == 0 reduces to the solid-cylinder exit path.
        // cyl_inner_half_length_ == cyl_half_length_ for a through-bore, which
        // is what configure_cylindrical() gives unless told otherwise.
        return cyl_shell_path(local_pos, local_dir,
                              cyl_inner_r_, cyl_inner_half_length_,
                              cyl_radius_, cyl_half_length_);
    }

    if (shape_ == Shape::Sphere) {
        // Shell when a void center is present; inner_r == 0 → solid ball.
        return spherical_shell_path(sphere_center_, position, direction,
                                    sphere_inner_r_, sphere_radius_);
    }

    if (shape_ == Shape::Rectangular) {
        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);
        Eigen::Vector3d local_dir = rect_rotation_ * direction;

        // Hollow box shell: subtract ray overlap with the non-attenuating
        // inner void (a ray born in the shell can cross the void and re-enter
        // the far wall). Solid boxes keep the original exit-path expression
        // below, unchanged.
        if (rect_inner_half_dims_.minCoeff() > 1e-10) {
            return box_shell_path(local_pos, local_dir,
                                  rect_inner_half_dims_, rect_half_dims_);
        }

        auto hit = intersect_box(local_pos, local_dir,
                                 rect_half_dims_.x(), rect_half_dims_.y(),
                                 -rect_half_dims_.z(), rect_half_dims_.z());
        if (!hit || !hit->valid()) return 0.0;

        return std::max(hit->t_exit, 0.0);
    }

    if (shape_ == Shape::Marinelli) {
        if (!is_inside_marinelli(position)) return 0.0;
        return marinelli_boundary_distance(position, direction);
    }

    return 0.0;
}

double SourceGeometry::min_distance_to_boundary(
    const Eigen::Vector3d& position, bool include_shields) const {
    // Envelope = source-material extent, plus (when include_shields) the
    // summed shield thicknesses per axis. For a convex envelope the nearest
    // boundary point from an interior point is the perpendicular foot on the
    // closest face, so the per-axis min below is exact (never an over-estimate).
    switch (shape_) {
    case Shape::Point: {
        // Point source: no material volume — electrons originate in the
        // spherical shield shells, which fill outward from the point.
        if (!include_shields) return 0.0;
        double r_out = 0.0;
        for (const auto& layer : shields()) r_out += layer.scalar_thickness();
        double r = (position - point_position_).norm();
        return std::max(r_out - r, 0.0);
    }
    case Shape::Cylindrical: {
        double r_out = cyl_radius_;
        double l_out = cyl_half_length_;
        if (include_shields)
            for (const auto& layer : shields()) { r_out += layer.tx; l_out += layer.tz; }
        Eigen::Vector3d local = cyl_rotation_ * (position - cyl_center_);
        double rho = std::hypot(local.x(), local.y());
        double d = std::min(r_out - rho, l_out - std::abs(local.z()));
        // Hollow bore: the inner wall is also a material boundary when we
        // measure to the source-material surface (include_shields == false) -
        // whether the inside is a void or a core.  A closed cavity
        // (cyl_inner_half_length_ < cyl_half_length_) adds its end faces, on the
        // same conservative footing as the box below.
        if (!include_shields && cyl_inner_r_ > 1e-10) {
            d = std::min(d, rho - cyl_inner_r_);
            const double beyond = std::abs(local.z()) - cyl_inner_half_length_;
            if (beyond > 0.0)
                d = std::min(d, beyond);
        }
        return std::max(d, 0.0);
    }
    case Shape::Sphere: {
        double r_out = sphere_radius_;
        if (include_shields)
            for (const auto& layer : shields()) r_out += layer.scalar_thickness();
        double r = (position - sphere_center_).norm();
        double d = r_out - r;
        if (!include_shields && sphere_inner_r_ > 1e-10)
            d = std::min(d, r - sphere_inner_r_);
        return std::max(d, 0.0);
    }
    case Shape::Rectangular: {
        Eigen::Vector3d h = rect_half_dims_;
        if (include_shields)
            for (const auto& layer : shields()) {
                h.x() += layer.tx; h.y() += layer.ty; h.z() += layer.tz;
            }
        Eigen::Vector3d local = rect_rotation_ * (position - rect_center_);
        double dx = h.x() - std::abs(local.x());
        double dy = h.y() - std::abs(local.y());
        double dz = h.z() - std::abs(local.z());
        double d = std::min({dx, dy, dz});
        // Hollow shell: the inner-void wall is also a material boundary when
        // we measure to the source-material surface (include_shields == false).
        // For each axis where the point lies beyond the inner box, the
        // perpendicular foot on that inner face is a candidate. Conservative
        // (<= true distance) since the shell is non-convex — same caveat as
        // the cylinder bore.
        if (!include_shields && rect_inner_half_dims_.minCoeff() > 1e-10) {
            for (int i = 0; i < 3; ++i) {
                const double beyond = std::abs(local[i]) - rect_inner_half_dims_[i];
                if (beyond > 0.0)
                    d = std::min(d, beyond);
            }
        }
        return std::max(d, 0.0);
    }
    case Shape::Marinelli: {
        // Conservative (≤ true) distance from an interior water point to the
        // nearest water boundary (beaker wall or re-entrant well void) of the
        // L-shaped fill.  Under-estimate is safe for the containment test; the
        // beaker walls ARE the boundary, so include_shields does not extend it.
        // Water = disk (z∈[z_bot,z_bk], r∈[0,outer]) ∪ ring (z∈[z_bk,z_we],
        // r∈[well,outer]).
        const double rho = std::hypot(position.x(), position.y());
        const double z   = position.z();
        double d = marinelli_outer_r_ - rho;            // outer beaker wall
        d = std::min(d, z - marinelli_z_bot_);          // beaker bottom
        if (z >= marinelli_z_bk_) {                      // ring region
            d = std::min(d, rho - marinelli_well_r_);   // well inner wall (void)
            d = std::min(d, marinelli_z_we_ - z);       // ring top
        } else if (rho < marinelli_well_r_) {           // disk under the well
            d = std::min(d, marinelli_z_bk_ - z);       // well floor (void above)
        }
        return std::max(d, 0.0);
    }
    default:
        return 0.0;
    }
}

double SourceGeometry::min_distance_within_material(
    const Eigen::Vector3d& position,
    const Material& material) const
{
    // If every finite source region has the same Material object, the nearest
    // material boundary is the outer source+shield envelope. Pointer identity
    // is intentional: distinct Material objects may differ in density or
    // composition even when their names happen to match.
    const Material* homogeneous_material = nullptr;
    bool homogeneous = true;
    auto include_region = [&](const Material* candidate) {
        if (!candidate || !homogeneous)
            return;
        if (!homogeneous_material)
            homogeneous_material = candidate;
        else if (homogeneous_material != candidate)
            homogeneous = false;
    };

    include_region(source_material_);
    for (const auto& layer : shields())
        include_region(layer.material);
    // A cored source always has an interface at the source/core wall - and
    //  possibly a void behind the cores - so it is never homogeneous; the same
    //  rule as the single-material electron-containment test.
    if (has_attenuating_interior_)
        homogeneous = false;

    if (homogeneous && homogeneous_material == &material)
        return min_distance_to_boundary(position, /*include_shields=*/true);

    // In a heterogeneous geometry the source-material surface is still a
    // cheap, conservative boundary. Shield-born particles fall back to the
    // full walk because not every supported shape has a general nearest-layer
    // distance calculation.
    if (source_material_ == &material)
        return min_distance_to_boundary(position, /*include_shields=*/false);

    return 0.0;
}

double SourceGeometry::compute_transmission(const Eigen::Vector3d& position,
                                             const Eigen::Vector3d& direction,
                                             double energy_MeV) const {
    double total_mu_L = 0.0;
    const auto attenuation_mu = [&](const Material* material) {
        const MacroscopicXS xs = material->macroscopic_xs(energy_MeV);
        return rayleigh_enabled_
            ? xs.mu_total()
            : xs.mu_pe + xs.mu_cs + xs.mu_pp;
    };

    // With cores the concentric stack is walked once, in order, and this is just
    //  the sum over the segments it returns; without them, the original per-layer
    //  loops below run unchanged.  (Attenuation itself does not care about order,
    //  but having one geometry description rather than two does.)
    if (has_attenuating_interior_) {
        std::vector<SourcePathSegment> segments;
        trace_cored_segments(position, direction, segments, SIZE_MAX);
        for (const SourcePathSegment& seg : segments) {
            if (seg.material && seg.length > 0.0)
                total_mu_L += attenuation_mu(seg.material) * seg.length;
        }
        return std::exp(-total_mu_L);
    }

    // Source material attenuation
    if (source_material_) {
        double path = source_material_path(position, direction);
        if (path > 0.0) {
            total_mu_L += attenuation_mu(source_material_) * path;
        }
    }

    // Shielding layers
    if (shape_ == Shape::Point) {
        // Spherical shells: path = thickness (direction-independent)
        for (const auto& layer : shields()) {
            assert(layer.is_uniform());
            total_mu_L += attenuation_mu(layer.material)
                * layer.scalar_thickness();
        }
    } else if (shape_ == Shape::Cylindrical) {
        // Cylindrical shells: compute path through each shell
        double inner_r = cyl_radius_;
        double inner_half_z = cyl_half_length_;

        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);
        Eigen::Vector3d local_dir = cyl_rotation_ * direction;

        for (const auto& layer : shields()) {
            assert(layer.tx == layer.ty);
            double outer_r = inner_r + layer.tx;
            double outer_half_z = inner_half_z + layer.tz;

            auto outer_hit = intersect_cylinder(local_pos, local_dir,
                                                outer_r, -outer_half_z, outer_half_z);
            auto inner_hit = intersect_cylinder(local_pos, local_dir,
                                                inner_r, -inner_half_z, inner_half_z);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;

                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    // Path from inner exit to outer exit
                    if (t1_in < t1_out) {
                        path = t1_out - t1_in;
                    }
                    // Also front segment if applicable
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) {
                        path += t0_in - t0_out;
                    }
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 0.0) {
                total_mu_L += attenuation_mu(layer.material) * path;
            }

            inner_r = outer_r;
            inner_half_z = outer_half_z;
        }
    } else if (shape_ == Shape::Rectangular) {
        double inner_hx = rect_half_dims_.x();
        double inner_hy = rect_half_dims_.y();
        double inner_hz = rect_half_dims_.z();

        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);
        Eigen::Vector3d local_dir = rect_rotation_ * direction;

        for (const auto& layer : shields()) {
            double outer_hx = inner_hx + layer.tx;
            double outer_hy = inner_hy + layer.ty;
            double outer_hz = inner_hz + layer.tz;

            auto outer_hit = intersect_box(local_pos, local_dir,
                                           outer_hx, outer_hy, -outer_hz, outer_hz);
            auto inner_hit = intersect_box(local_pos, local_dir,
                                           inner_hx, inner_hy, -inner_hz, inner_hz);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;

                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) {
                        path = t1_out - t1_in;
                    }
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) {
                        path += t0_in - t0_out;
                    }
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 0.0) {
                total_mu_L += attenuation_mu(layer.material) * path;
            }

            inner_hx = outer_hx;
            inner_hy = outer_hy;
            inner_hz = outer_hz;
        }
    } else if (shape_ == Shape::Sphere) {
        // Spherical shells: ray-trace each shell (the emission point is offset
        // from the center, so the path is direction-dependent — unlike a point
        // source emitting from the exact center).
        double inner_r = sphere_radius_;
        for (const auto& layer : shields()) {
            assert(layer.is_uniform());
            double outer_r = inner_r + layer.scalar_thickness();
            double path = spherical_shell_path(sphere_center_, position, direction,
                                               inner_r, outer_r);
            if (path > 0.0)
                total_mu_L += attenuation_mu(layer.material) * path;
            inner_r = outer_r;
        }
    } else if (shape_ == Shape::Marinelli) {
        // Beaker wall: path depends on exit surface. Well opening has no wall.
        double wall_factor = marinelli_wall_factor(position, direction);
        if (wall_factor > 0.0) {
            for (const auto& layer : shields()) {
                assert(layer.is_uniform());
                total_mu_L += attenuation_mu(layer.material)
                    * layer.scalar_thickness() * wall_factor;
            }
        }
    }

    return std::exp(-total_mu_L);
}

double SourceGeometry::point_source_transmission(double energy_MeV) const {
    double total_mu_L = 0.0;
    for (const auto& layer : shields()) {
        assert(layer.is_uniform());
        const MacroscopicXS xs =
            layer.material->macroscopic_xs(energy_MeV);
        const double mu = rayleigh_enabled_
            ? xs.mu_total()
            : xs.mu_pe + xs.mu_cs + xs.mu_pp;
        total_mu_L += mu * layer.scalar_thickness();
    }
    return std::exp(-total_mu_L);
}

SourceGeometry::SourceTransmissionResult SourceGeometry::compute_transmission_fep_only(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_MeV,
    std::mt19937_64& rng) const
{
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double energy_keV = energy_MeV * 1e3;
    double weight = 1.0;
    Eigen::Vector3d dir = direction;
    double total_path = 0.0;  // accumulate path for exit_position

    // Helper: process a single material layer (path length known)
    auto process_layer = [&](const Material* mat, double path) {
        if (!mat || path <= 0.0) return;
        total_path += path;
        MacroscopicXS xs = mat->macroscopic_xs(energy_MeV);
        double mu_no_rs = xs.mu_pe + xs.mu_cs + xs.mu_pp;
        double mu_rs = rayleigh_enabled_ ? xs.mu_rs : 0.0;

        weight *= std::exp(-mu_no_rs * path);

        // Stochastic Rayleigh
        double lambda_rs = mu_rs * path;
        int n_scatter = 0;
        if (lambda_rs > 0.01) {
            std::poisson_distribution<int> poisson(lambda_rs);
            n_scatter = poisson(rng);
        } else if (lambda_rs > 1e-10) {
            n_scatter = (uniform(rng) < lambda_rs) ? 1 : 0;
        }

        if (n_scatter > 0) {
            int Z_rs = mat->select_element(energy_MeV, 2 /*RS*/, uniform(rng));
            for (int i = 0; i < n_scatter; ++i) {
                double cos_theta = sample_rayleigh_cos_theta(Z_rs, energy_keV, rng);
                double phi = kTwoPi * uniform(rng);
                dir = rotate_direction(dir, cos_theta, phi);
            }
        }
    };

    // With cores, one ordered walk replaces the per-layer loops (see
    //  compute_transmission).  process_layer can turn the photon (Rayleigh), and
    //  every layer past a turn must be traced along the NEW direction - which is
    //  what the per-shape loops below do by re-deriving the direction each
    //  iteration.  Here that means re-tracing the rest of the stack from the
    //  turn point.  Voids advance the position without being processed.
    if (has_attenuating_interior_) {
        constexpr int kMaxRetrace = 32;
        Eigen::Vector3d pos = position;
        std::vector<SourcePathSegment> segments;
        for (int retrace = 0; retrace < kMaxRetrace; ++retrace) {
            const Eigen::Vector3d dir_before = dir;
            trace_cored_segments(pos, dir_before, segments, SIZE_MAX);

            bool turned = false;
            for (const SourcePathSegment& seg : segments) {
                process_layer(seg.material, seg.length);
                pos += dir_before * seg.length;
                if (dir != dir_before) { turned = true; break; }
            }
            if (!turned) break;
        }
        return {weight, dir, pos};
    }

    // Source material attenuation
    if (source_material_) {
        double path = source_material_path(position, direction);
        process_layer(source_material_, path);
    }

    // Shielding layers
    if (shape_ == Shape::Point) {
        for (const auto& layer : shields()) {
            assert(layer.is_uniform());
            process_layer(layer.material, layer.scalar_thickness());
        }
    } else if (shape_ == Shape::Cylindrical) {
        double inner_r = cyl_radius_;
        double inner_half_z = cyl_half_length_;
        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);

        for (const auto& layer : shields()) {
            // Recompute local_dir each iteration: process_layer may change dir
            // via Rayleigh scattering in an earlier layer.
            Eigen::Vector3d local_dir = cyl_rotation_ * dir;

            assert(layer.tx == layer.ty);
            double outer_r = inner_r + layer.tx;
            double outer_half_z = inner_half_z + layer.tz;

            auto outer_hit = intersect_cylinder(local_pos, local_dir,
                                                outer_r, -outer_half_z, outer_half_z);
            auto inner_hit = intersect_cylinder(local_pos, local_dir,
                                                inner_r, -inner_half_z, inner_half_z);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            process_layer(layer.material, path);

            inner_r = outer_r;
            inner_half_z = outer_half_z;
        }
    } else if (shape_ == Shape::Rectangular) {
        double inner_hx = rect_half_dims_.x();
        double inner_hy = rect_half_dims_.y();
        double inner_hz = rect_half_dims_.z();
        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);

        for (const auto& layer : shields()) {
            // Recompute local_dir each iteration: process_layer may change dir
            // via Rayleigh scattering in an earlier layer.
            Eigen::Vector3d local_dir = rect_rotation_ * dir;

            double outer_hx = inner_hx + layer.tx;
            double outer_hy = inner_hy + layer.ty;
            double outer_hz = inner_hz + layer.tz;

            auto outer_hit = intersect_box(local_pos, local_dir,
                                           outer_hx, outer_hy, -outer_hz, outer_hz);
            auto inner_hit = intersect_box(local_pos, local_dir,
                                           inner_hx, inner_hy, -inner_hz, inner_hz);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            process_layer(layer.material, path);

            inner_hx = outer_hx;
            inner_hy = outer_hy;
            inner_hz = outer_hz;
        }
    } else if (shape_ == Shape::Sphere) {
        // Spherical shells: ray-trace each shell with the current direction
        // (Rayleigh in an earlier layer may have changed `dir`).
        double inner_r = sphere_radius_;
        for (const auto& layer : shields()) {
            assert(layer.is_uniform());
            double outer_r = inner_r + layer.scalar_thickness();
            double path = spherical_shell_path(sphere_center_, position, dir,
                                               inner_r, outer_r);
            process_layer(layer.material, path);
            inner_r = outer_r;
        }
    } else if (shape_ == Shape::Marinelli) {
        // Beaker wall: path depends on exit surface. Well opening has no wall.
        double wall_factor = marinelli_wall_factor(position, direction);
        if (wall_factor > 0.0) {
            for (const auto& layer : shields()) {
                assert(layer.is_uniform());
                process_layer(layer.material, layer.scalar_thickness() * wall_factor);
            }
        }
    }

    return {weight, dir, position + dir * total_path};
}

std::vector<SourceGeometry::SourcePathSegment> SourceGeometry::trace_source_segments(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV) const
{
    std::vector<SourcePathSegment> segments;
    trace_source_segments(position, direction, energy_keV, segments);
    return segments;
}

void SourceGeometry::trace_source_segments(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double /*energy_keV*/,
    std::vector<SourcePathSegment>& segments,
    std::size_t max_segments) const
{
    segments.clear();
    segments.reserve(1 + shields().size());

    // An attenuating core makes segment ORDER load-bearing, so the whole stack
    //  goes through the ordered walker instead of the per-layer loops below.
    //  There is no attenuating interior for any source that predates
    //  add_core(), so the loops
    //  below - and every number they have ever produced - are untouched.
    if (has_attenuating_interior_) {
        trace_cored_segments(position, direction, segments, max_segments);
        return;
    }

    // Source material (extended sources only)
    if (source_material_ && shape_ != Shape::Point) {
        double path = source_material_path(position, direction);
        if (path > 1e-10) {
            segments.push_back({source_material_, path});
            if (segments.size() >= max_segments) return;
        }
    }

    // Shielding layers
    if (shape_ == Shape::Point || shape_ == Shape::Sphere) {
        // Concentric spherical shells. Point: shells grow from the point
        // (inner_r starts at 0). Sphere: shells wrap the source ball (inner_r
        // starts at the source outer radius). Both use the shared shell-path
        // helper (the emission point may be offset from the center for a
        // spherical volume source, so the path is direction-dependent).
        const Eigen::Vector3d& center =
            (shape_ == Shape::Point) ? point_position_ : sphere_center_;
        double inner_r = (shape_ == Shape::Point) ? 0.0 : sphere_radius_;
        for (const auto& layer : shields()) {
            assert(layer.is_uniform());
            double outer_r = inner_r + layer.scalar_thickness();
            double path = spherical_shell_path(center, position, direction,
                                               inner_r, outer_r);
            if (path > 1e-10) {
                segments.push_back({layer.material, path});
                if (segments.size() >= max_segments) return;
            }
            inner_r = outer_r;
        }
    } else if (shape_ == Shape::Cylindrical) {
        double inner_r = cyl_radius_;
        double inner_half_z = cyl_half_length_;
        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);
        Eigen::Vector3d local_dir = cyl_rotation_ * direction;

        for (const auto& layer : shields()) {
            assert(layer.tx == layer.ty);
            double outer_r = inner_r + layer.tx;
            double outer_half_z = inner_half_z + layer.tz;

            auto outer_hit = intersect_cylinder(local_pos, local_dir,
                                                outer_r, -outer_half_z, outer_half_z);
            auto inner_hit = intersect_cylinder(local_pos, local_dir,
                                                inner_r, -inner_half_z, inner_half_z);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 1e-10) {
                segments.push_back({layer.material, path});
                if (segments.size() >= max_segments) return;
            }

            inner_r = outer_r;
            inner_half_z = outer_half_z;
        }
    } else if (shape_ == Shape::Rectangular) {
        double inner_hx = rect_half_dims_.x();
        double inner_hy = rect_half_dims_.y();
        double inner_hz = rect_half_dims_.z();
        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);
        Eigen::Vector3d local_dir = rect_rotation_ * direction;

        for (const auto& layer : shields()) {
            double outer_hx = inner_hx + layer.tx;
            double outer_hy = inner_hy + layer.ty;
            double outer_hz = inner_hz + layer.tz;

            auto outer_hit = intersect_box(local_pos, local_dir,
                                           outer_hx, outer_hy, -outer_hz, outer_hz);
            auto inner_hit = intersect_box(local_pos, local_dir,
                                           inner_hx, inner_hy, -inner_hz, inner_hz);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 1e-10) {
                segments.push_back({layer.material, path});
                if (segments.size() >= max_segments) return;
            }

            inner_hx = outer_hx;
            inner_hy = outer_hy;
            inner_hz = outer_hz;
        }
    } else if (shape_ == Shape::Marinelli) {
        // Beaker wall: path depends on exit surface. Well opening has no wall.
        // Only add wall if photon is inside sample volume.
        double wall_factor = marinelli_wall_factor(position, direction);
        if (wall_factor > 0.0) {
            for (const auto& layer : shields()) {
                assert(layer.is_uniform());
                double path = layer.scalar_thickness() * wall_factor;
                if (path > 1e-10) {
                    segments.push_back({layer.material, path});
                    if (segments.size() >= max_segments) return;
                }
            }
        }
    }
}

void SourceGeometry::trace_cored_segments(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    std::vector<SourcePathSegment>& segments,
    std::size_t max_segments) const
{
    // Clear, like trace_source_segments() does.  compute_transmission_fep_only()
    //  reuses ONE buffer across its retrace iterations, so appending here made
    //  every Rayleigh turn re-process the segments already walked: the
    //  attenuation was applied twice and `pos` ran away past the source (a 2 cm
    //  source reported exit positions over 100 cm out, feeding a bogus air gap).
    segments.clear();

    // Cores only make sense for an extended source with a hollow interior.
    assert(shape_ == Shape::Sphere || shape_ == Shape::Cylindrical
           || shape_ == Shape::Rectangular);

    Eigen::Vector3d local_pos, local_dir;
    switch (shape_) {
    case Shape::Sphere:
        // A sphere is rotation invariant, so the local frame is a translation.
        local_pos = position - sphere_center_;
        local_dir = direction;
        break;
    case Shape::Cylindrical:
        local_pos = cyl_rotation_ * (position - cyl_center_);
        local_dir = cyl_rotation_ * direction;
        break;
    case Shape::Rectangular:
        local_pos = rect_rotation_ * (position - rect_center_);
        local_dir = rect_rotation_ * direction;
        break;
    default:
        return;
    }

    // layers_ IS the stack the march wants: innermost-first, each entry
    //  carrying the outer boundary of the region its material fills, with an
    //  unfilled cavity kept as a null-material entry so the along-ray distances
    //  stay right.  It used to be stitched together here from an outermost-first
    //  core list, the source material and an innermost-first shield list, on
    //  every single ray; maintaining it in the mutators instead is what lets the
    //  GDML export emit the same geometry this traces.
    trace_concentric(shape_, layers_, local_pos, local_dir, segments, max_segments);
}

double SourceGeometry::no_interaction_probability(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV,
    double* path_cm) const
{
    const auto segments = trace_source_segments(position, direction, energy_keV);
    const double energy_MeV = energy_keV * 1e-3;
    double tau = 0.0;
    double path = 0.0;
    for (const auto& seg : segments) {
        if (!seg.material || seg.length <= 0.0) continue;
        path += seg.length;
        const MacroscopicXS xs =
            seg.material->macroscopic_xs(energy_MeV);
        const double mu = rayleigh_enabled_
            ? xs.mu_total()
            : xs.mu_pe + xs.mu_cs + xs.mu_pp;
        tau += mu * seg.length;
    }
    if (path_cm) *path_cm = path;
    return std::exp(-tau);
}

SourceGeometry::SourceFullTransportResult SourceGeometry::transport_source_photon(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV,
    std::mt19937_64& rng) const
{
    // Compton-angle biasing applies only to the primary photon (this public
    // entry point); annihilation/brems recursion always runs analog. Never
    // biased for Marinelli (re-entry transports share this entry point and
    // the auto-policy never enables it there).
    const int budget = (compton_bias_.cone_fraction > 0.0 &&
                        shape_ != Shape::Marinelli)
                           ? compton_bias_.max_biased_vertices
                           : 0;
    return transport_source_photon_impl(position, direction, energy_keV, rng,
                                        budget);
}

SourceGeometry::SourcePrimaryElectronResult
SourceGeometry::transport_source_electron(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double kinetic_energy_keV,
    std::mt19937_64& rng,
    bool is_positron,
    SourceSecondaryPhoton::Origin brems_origin) const
{
    SourcePrimaryElectronResult out;
    out.stop_position = position;
    if (!(kinetic_energy_keV > 0.0) || !std::isfinite(kinetic_energy_keV))
        return out;

    Eigen::Vector3d dir = direction;
    const double norm = dir.norm();
    if (!(norm > 0.0) || !std::isfinite(norm))
        return out;
    dir /= norm;

    std::vector<SourcePathSegment> birth_segments;
    trace_source_segments(position, dir, kinetic_energy_keV, birth_segments, 1);
    if (birth_segments.empty() || !birth_segments.front().material) {
        out.escaped = true;
        out.escaped_electron = {position, dir, kinetic_energy_keV};
        return out;
    }

    const ElectronSourceWalkResult walk =
        ElectronCsda::instance().walk_in_source_geometry(
            *this, *birth_segments.front().material, position, dir,
            kinetic_energy_keV, rng, is_positron);
    out.escaped = walk.escaped;
    out.traversed_gcm2 = walk.traversed_gcm2;
    out.stop_position = walk.stop_position;
    out.annihilated_in_flight = walk.annihilated_in_flight;
    out.residual_KE_at_annih_keV = walk.residual_KE_at_annih_keV;
    if (walk.escaped) {
        out.escaped_electron = {
            walk.exit_position, walk.exit_direction, walk.exit_KE_keV,
            is_positron};
    }

    auto transport_and_append =
        [&](const Eigen::Vector3d& photon_position,
            const Eigen::Vector3d& photon_direction,
            double photon_energy_keV,
            SourceSecondaryPhoton::Origin origin,
            bool pair_history) {
        auto transported = transport_source_photon_impl(
            photon_position, photon_direction, photon_energy_keV, rng,
            /*bias_budget=*/0);
        if (transported.survived) {
            out.photons.push_back(
                {transported.position, transported.direction,
                 transported.energy_keV, origin,
                 transported.truth.compton_interactions,
                 transported.truth.rayleigh_interactions,
                 pair_history || transported.truth.pair_history});
        }
        for (auto& photon : transported.secondaries)
            out.photons.push_back(std::move(photon));
    };

    for (const BremsPhoton& photon : walk.brems_photons) {
        if (photon.impact_Z > 0) {
            transport_and_append(
                photon.position, photon.direction, photon.energy_keV,
                SourceSecondaryPhoton::Origin::ImpactIonization, false);
            continue;
        }
        if (!source_brems_enabled_) continue;
        {
            const SourceSecondaryPhoton::Origin photon_origin =
                photon.charged_generation == 0
                    ? brems_origin
                    : SourceSecondaryPhoton::Origin::
                        SecondaryElectronBrems;
            if (record_truth_labels_) {
                out.created_brems_photons.push_back(
                    {photon.position, photon.direction, photon.energy_keV,
                     photon_origin, 0, 0, false});
            }
            transport_and_append(
                photon.position, photon.direction, photon.energy_keV,
                photon_origin, false);
        }
    }

    // A positron that escapes the source carries its annihilation photons with
    // it; they are not part of the source-escape photon measure.
    if (!is_positron || walk.escaped)
        return out;

    constexpr double kElectronMassKeV = 510.998950;
    double energy1 = kElectronMassKeV;
    double energy2 = kElectronMassKeV;
    Eigen::Vector3d direction1 = sample_isotropic_direction(rng);
    Eigen::Vector3d direction2 = -direction1;
    auto origin = SourceSecondaryPhoton::Origin::PrimaryPositronAnnihilation;

    if (walk.annihilated_in_flight) {
        const double total_energy =
            2.0 * kElectronMassKeV + walk.residual_KE_at_annih_keV;
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        const double fraction =
            0.5 + (uniform(rng) - 0.5)
                * walk.residual_KE_at_annih_keV / total_energy;
        energy1 = fraction * total_energy;
        energy2 = total_energy - energy1;
        direction1 = sample_isotropic_direction(rng);
        direction2 = sample_isotropic_direction(rng);
        origin = SourceSecondaryPhoton::Origin::
            PrimaryPositronAnnihilationInFlight;
    }

    transport_and_append(
        walk.stop_position, direction1, energy1, origin, true);
    transport_and_append(
        walk.stop_position, direction2, energy2, origin, true);
    return out;
}

SourceGeometry::SourceFullTransportResult SourceGeometry::transport_source_photon_impl(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV,
    std::mt19937_64& rng,
    int bias_budget) const
{
    constexpr int kMaxSourceScatters = 50;
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    Eigen::Vector3d pos = position;
    Eigen::Vector3d dir = direction;
    double energy = energy_keV;

    // Accumulate source electrons and secondary photons (annihilation gammas,
    // bremsstrahlung) across multiple Compton scatters
    std::vector<SourceElectron> result_electrons;
    std::vector<SourceSecondaryPhoton> result_secondaries;
    SourceFullTransportResult::TruthLabels truth;

    // True once any interaction (Compton/PE/PP/Rayleigh) has occurred,
    // including scatters in earlier re-trace iterations.
    bool any_interaction = false;

    // Accumulated Compton-angle importance-sampling weight (1.0 = analog).
    double w_bias = 1.0;

    // Transport a brems photon through the remaining source geometry and keep
    // the survivor (and any nested secondaries) as event secondaries.
    auto transport_brems_photon =
        [&](const BremsPhoton& bp,
            SourceSecondaryPhoton::Origin brems_origin) {
        const SourceSecondaryPhoton::Origin photon_origin =
            bp.impact_Z > 0
                ? SourceSecondaryPhoton::Origin::ImpactIonization
                : bp.charged_generation == 0
                ? brems_origin
                : SourceSecondaryPhoton::Origin::
                    SecondaryElectronBrems;
        auto br = transport_source_photon_impl(
            bp.position, bp.direction, bp.energy_keV, rng, /*bias_budget=*/0);
        if (br.survived) {
            result_secondaries.push_back(
                {br.position, br.direction, br.energy_keV,
                 photon_origin,
                 br.truth.compton_interactions,
                 br.truth.rayleigh_interactions,
                 br.truth.pair_history});
        }
        for (auto& sec : br.secondaries)
            result_secondaries.push_back(std::move(sec));
        for (auto& e : br.source_electrons)
            result_electrons.push_back(std::move(e));
    };

    // Transport a relaxation x-ray (K or L line) isotropically from `at`
    // through the remaining source geometry. Lines below 1 keV stay local.
    auto transport_relaxation_xray = [&](const Eigen::Vector3d& at, double x_keV) {
        if (x_keV < 1.0) return;
        auto x = transport_source_photon_impl(
            at, sample_isotropic_direction(rng), x_keV, rng, /*bias_budget=*/0);
        if (x.survived) {
            result_secondaries.push_back(
                {x.position, x.direction, x.energy_keV,
                 SourceSecondaryPhoton::Origin::Fluorescence,
                 x.truth.compton_interactions,
                 x.truth.rayleigh_interactions,
                 x.truth.pair_history});
        }
        for (auto& sec : x.secondaries)
            result_secondaries.push_back(std::move(sec));
        for (auto& e : x.source_electrons)
            result_electrons.push_back(std::move(e));
    };

    // Relax a K vacancy of element Z: the K x-ray (or Auger), then the L
    // vacancy a K-alpha leaves when the element has retained L lines.
    auto relax_k_vacancy = [&](const Eigen::Vector3d& at, const FluorescenceData& fl,
                               const LFluorescenceData* fl_L) {
        const KRelaxation relax = sample_k_relaxation(fl, rng);
        LRelaxation l_relax;
        if (relax.line_keV > 0.0 && has_l_lines(fl_L)
            && relax.filling_subshell >= kEadlL1 && relax.filling_subshell <= kEadlL3)
            l_relax = sample_l_relaxation(*fl_L, relax.filling_subshell - kEadlL1, rng);
        transport_relaxation_xray(at, relax.line_keV);
        transport_relaxation_xray(at, l_relax.line_keV);
    };

    // Detect whether the source geometry is a single material (no internal
    // interfaces). Computed once per primary.  A cored source always has an
    // interface at the source/core wall - and possibly a void behind the cores -
    // so it never takes the single-material electron-containment fast path
    // below.  Erring toward `false` only costs the slower, exact route.
    const Material* first_mat = source_material_;
    bool single_material = !has_attenuating_interior_;
    for (const auto& layer : shields()) {
        if (!layer.material) continue;
        if (!first_mat) first_mat = layer.material;
        else if (layer.material != first_mat) { single_material = false; break; }
    }

    // Handle an electron (Compton recoil, photoelectron, pair e+/e-) born in
    // the source geometry at e_pos in material `mat`.
    //
    // ALL shapes (incl. Marinelli as of the regime-aware escape work): Moliere
    // condensed-history walk through the source segments — emits bremsstrahlung
    // along the way and detects escape through the outer boundary (the direction
    // diffusion is essential in dense shields; straight-line CSDA overestimates
    // escape several-fold).  Marinelli previously used the legacy straight-line
    // try_source_electron for escape; it now shares this walk + the regime-aware
    // exit-state skin-escape gate (light water → fully analog), so cfg 8 joins
    // the gated/validated set.
    //
    // Perf: the geometry-tracing Moliere walk dominates the cost of the
    // secondary channels — it re-traces the source segments at every substep
    // (profiling: trace_source_segments is the second-hottest leaf). Most
    // electrons are born deep enough that they stop without leaving their
    // birth material. For those we (a) cannot escape and (b) can sample brems
    // with the geometry-free ElectronCsda::sample_brems_in_material (the same
    // SB physics, creation-point approximation — exact when the range is small
    // vs the medium, which is precisely the contained case; this is the
    // sampler already used for Marinelli). The walk is reserved for electrons
    // that can reach a material interface or escape.
    //
    // Containment test: the electron stays within material `mat` if its CSDA
    // range in `mat` is below the distance to the nearest surface of the `mat`
    // region. `d_stay` always measures to the boundary of the birth material's
    // own region, and containment means the electron never leaves `mat`, so the
    // matching reach bound is `mat`'s own range/density — exact, not merely
    // conservative (using any other material's range/density would be a
    // mismatch, since the g/cm^2 CSDA range is itself Z/A-dependent). The
    // per-config escape-energy floor (set_source_electron_threshold, default
    // 50 keV) is the other knob.
    auto process_electron = [&](const Eigen::Vector3d& e_pos,
                                const Eigen::Vector3d& e_dir,
                                double T_keV, const Material& mat,
                                SourceSecondaryPhoton::Origin brems_origin,
                                bool is_positron = false)
        -> ElectronSourceWalkResult {
        ElectronSourceWalkResult no_walk;
        no_walk.stop_position = e_pos;
        // An electron that is not walked emits the impact-ionization K x-rays
        // of its whole slowing-down at its birth point.
        auto impact_in_place = [&]() {
            if (!source_impact_ionization_enabled_) return;
            const auto& csda = ElectronCsda::instance();
            const double ke = std::min(T_keV, ElectronCsda::impact_ionization_max_keV());
            std::vector<ElectronCsda::ImpactXray> xrays;
            csda.sample_impact_xrays(mat, ke, csda.impact_xray_yield(mat, ke), true, rng, xrays);
            for (const auto& x : xrays) {
                BremsPhoton bp{e_pos, sample_isotropic_direction(rng), x.line_keV};
                bp.impact_Z = static_cast<std::uint8_t>(x.Z);
                bp.impact_filling_subshell = x.filling_subshell;
                transport_brems_photon(bp, brems_origin);
            }
        };
        const bool want_brems = source_brems_enabled_ &&
            T_keV > kMoliereBremsThreshold_keV;
        const bool want_escape = source_electron_enabled_ &&
            T_keV > source_electron_threshold_keV_;
        const bool want_positron_transport =
            is_positron && source_positron_transport_enabled_;
        if (!want_brems && !want_escape && !want_positron_transport) {
            impact_in_place();
            return no_walk;
        }

        // The electron's reach (cm) in its birth material `mat`: the straight-
        // line displacement cannot exceed the CSDA path length. Exact for the
        // containment test below, which only ever asks whether the electron
        // stays inside the `mat` region.
        double range_cm = 0.0;
        if (mat.density() > 0.0) {
            range_cm =
                ElectronCsda::instance().range_gcm2_material(
                    mat, T_keV, want_positron_transport)
                / mat.density();
        }

        // Distance the electron is guaranteed to remain inside material `mat`:
        // the outer boundary for a single-material geometry, otherwise the
        // source-material surface for a source-material electron. Shield-born
        // electrons in a multi-material geometry take the full walk (0 here).
        double d_stay = 0.0;
        if (single_material)
            d_stay = min_distance_to_boundary(e_pos, /*include_shields=*/true);
        else if (&mat == source_material_)
            d_stay = min_distance_to_boundary(e_pos, /*include_shields=*/false);

        if (!want_positron_transport && range_cm < d_stay) {
            // Contained in `mat`: cannot escape, cannot cross an interface.
            if (want_brems) {
                auto brems = ElectronCsda::instance().sample_brems_in_material(
                    mat, e_pos, e_dir, T_keV, rng);
                for (const auto& bp : brems)
                    transport_brems_photon(bp, brems_origin);
            }
            impact_in_place();
            return no_walk;  // no escaped electron possible
        }

        auto walk = ElectronCsda::instance().walk_in_source_geometry(
            *this, mat, e_pos, e_dir, T_keV, rng,
            want_positron_transport);

        for (const auto& bp : walk.brems_photons)
            if (source_brems_enabled_ || bp.impact_Z > 0)
                transport_brems_photon(bp, brems_origin);

        // A slow positron still annihilates wherever it stops, so the
        // kinetic-energy floor applies to electrons only.
        if (source_electron_enabled_ && walk.escaped &&
            (want_positron_transport ||
             walk.exit_KE_keV > source_electron_threshold_keV_)) {
            result_electrons.push_back(
                {walk.exit_position, walk.exit_direction, walk.exit_KE_keV,
                 want_positron_transport});
        }
        return walk;
    };

    std::vector<SourcePathSegment> segments;  // reused across retrace iters (E2)
    for (int retrace = 0; retrace < kMaxSourceScatters; ++retrace) {
        trace_source_segments(pos, dir, energy, segments);

        if (segments.empty()) {
            return {true, pos, dir, energy, std::move(result_secondaries),
                    std::move(result_electrons), any_interaction, w_bias, truth};
        }

        bool scattered = false;
        double distance_along_ray = 0.0;

        for (const auto& seg : segments) {
            // A null material is a void the trace kept in the list so the
            //  along-ray distance stays right (a partly-filled core); cross it.
            if (!seg.material) {
                distance_along_ray += seg.length;
                continue;
            }

            double energy_MeV = energy * 1e-3;
            MacroscopicXS xs = seg.material->macroscopic_xs(energy_MeV);
            if (!rayleigh_enabled_)
                xs.mu_rs = 0.0;
            double mu_total = xs.mu_total();

            if (mu_total <= 0.0) {
                distance_along_ray += seg.length;
                continue;
            }

            double s = sample_interaction_distance(mu_total, rng);

            if (s >= seg.length) {
                // No interaction in this segment — advance past it
                distance_along_ray += seg.length;
                continue;
            }

            // Interaction occurs at distance s within this segment
            pos += dir * (distance_along_ray + s);

            InteractionType type = select_interaction(xs, uniform(rng));
            any_interaction = true;

            switch (type) {
            case InteractionType::Photoelectric: {
                if (record_truth_labels_)
                    ++truth.photoelectric_interactions;
                // A K-shell vacancy relaxes as in the crystal: the K x-ray
                // is transported through the remaining source geometry, the
                // rest of the binding energy stays local, and the
                // photoelectron carries E - E_K.  Lines below 1 keV (low-Z
                // K shells) are absorbed on the spot.
                double pe_KE = energy;
                if (source_fluorescence_enabled_) {
                    const auto& xs_data = CrossSectionData::instance();
                    const int Z_pe = seg.material->select_element(
                        energy_MeV, 0 /*PE*/, uniform(rng));
                    const FluorescenceData* fl = xs_data.fluorescence(Z_pe);
                    const LFluorescenceData* fl_L = xs_data.l_fluorescence(Z_pe);
                    double f_K = 0.0;
                    bool k_shell = false;
                    if (fl && energy > fl->k_edge_keV) {
                        const double sigma_K =
                            xs_data.sigma_K_photoelectric(Z_pe, energy_MeV);
                        const double sigma_tot =
                            xs_data.sigma_photoelectric(Z_pe, energy_MeV);
                        f_K = (sigma_tot > 0.0)
                            ? std::min(sigma_K / sigma_tot, 1.0) : 1.0;
                        k_shell = uniform(rng) < f_K;
                    }
                    if (k_shell) {
                        pe_KE = energy - fl->k_edge_keV;
                        relax_k_vacancy(pos, *fl, fl_L);
                    } else if (has_l_lines(fl_L)) {
                        // Outer-shell absorption: an L vacancy (by the tabulated
                        // subshell shares) may radiate a retained L line.
                        const int s_l = sample_l_photoelectric_subshell(
                            *fl_L, energy, f_K, rng);
                        if (s_l >= 0 && energy > fl_L->sub[s_l].edge_keV) {
                            pe_KE = energy - fl_L->sub[s_l].edge_keV;
                            transport_relaxation_xray(
                                pos, sample_l_relaxation(*fl_L, s_l, rng).line_keV);
                        }
                    }
                }
                // Below the 200 keV brems floor the photoelectron neither
                // radiates nor (above the escape threshold) gets walked; with
                // impact ionization it still emits its K x-rays in place.
                if (pe_KE > kMoliereBremsThreshold_keV
                    || source_impact_ionization_enabled_) {
                    Eigen::Vector3d e_dir =
                        sample_photoelectron_direction(dir, pe_KE, rng);
                    process_electron(
                        pos, e_dir, pe_KE, *seg.material,
                        SourceSecondaryPhoton::Origin::PhotoelectronBrems);
                }
                return {false, pos, dir, energy, std::move(result_secondaries),
                        std::move(result_electrons), true, w_bias, truth};
            }

            case InteractionType::PairProduction: {
                if (record_truth_labels_) {
                    ++truth.pair_production_interactions;
                    truth.pair_history = true;
                }
                // Generate two back-to-back annihilation gammas (electron rest mass energy)
                // and transport each through the remaining source material.
                constexpr double kElectronMassKeV = 510.998950;

                // The e+/e- pair stopping in the shield: bremsstrahlung and
                // possible escape. Pair kinetic energy split sampled
                // uniformly; both particles start (approximately) along the
                // photon direction.
                double pair_T = energy - 2.0 * kElectronMassKeV;
                ElectronSourceWalkResult positron_walk;
                positron_walk.stop_position = pos;
                if (pair_T > 0.0) {
                    double u = uniform(rng);
                    process_electron(
                        pos, dir, u * pair_T, *seg.material,
                        SourceSecondaryPhoton::Origin::PairElectronBrems,
                        false);
                    positron_walk = process_electron(
                        pos, dir, (1.0 - u) * pair_T, *seg.material,
                        SourceSecondaryPhoton::Origin::PairElectronBrems,
                        true);
                }

                SourceFullTransportResult pp_result{false, pos, dir, energy,
                                                    std::move(result_secondaries),
                                                    std::move(result_electrons),
                                                    true, w_bias, truth};
                // With positron transport off, the 511-keV pair is emitted at
                // the pair vertex. Full-spectrum campaigns track the positron to its
                // stop/annihilation point, omit source annihilation photons
                // when it escapes, and include the in-flight continuum.
                if (source_positron_transport_enabled_
                    && positron_walk.escaped) {
                    return pp_result;
                }
                const Eigen::Vector3d annihilation_position =
                    source_positron_transport_enabled_
                    ? positron_walk.stop_position : pos;
                double annihilation_energy_1 = kElectronMassKeV;
                double annihilation_energy_2 = kElectronMassKeV;
                Eigen::Vector3d annihilation_direction_1 =
                    sample_isotropic_direction(rng);
                Eigen::Vector3d annihilation_direction_2 =
                    -annihilation_direction_1;
                if (source_positron_transport_enabled_
                    && positron_walk.annihilated_in_flight) {
                    const double total_energy =
                        2.0 * kElectronMassKeV
                        + positron_walk.residual_KE_at_annih_keV;
                    const double fraction =
                        0.5 + (uniform(rng) - 0.5)
                        * (positron_walk.residual_KE_at_annih_keV
                           / total_energy);
                    annihilation_energy_1 = fraction * total_energy;
                    annihilation_energy_2 =
                        total_energy - annihilation_energy_1;
                    annihilation_direction_1 =
                        sample_isotropic_direction(rng);
                    annihilation_direction_2 =
                        sample_isotropic_direction(rng);
                }
                const std::array<std::pair<Eigen::Vector3d, double>, 2>
                    annihilation_photons{{
                        {annihilation_direction_1, annihilation_energy_1},
                        {annihilation_direction_2, annihilation_energy_2},
                    }};
                const auto annihilation_origin =
                    source_positron_transport_enabled_
                    && positron_walk.annihilated_in_flight
                    ? SourceSecondaryPhoton::Origin::PairAnnihilationInFlight
                    : SourceSecondaryPhoton::Origin::PairAnnihilation;
                for (const auto& annihilation : annihilation_photons) {
                    auto ann = transport_source_photon_impl(
                        annihilation_position, annihilation.first,
                        annihilation.second, rng, /*bias_budget=*/0);
                    if (ann.survived) {
                        pp_result.secondaries.push_back(
                            {ann.position, ann.direction, ann.energy_keV,
                             annihilation_origin,
                             ann.truth.compton_interactions,
                             ann.truth.rayleigh_interactions,
                             true});
                    }
                    // Propagate secondaries and source electrons from recursive calls
                    for (auto& sec : ann.secondaries)
                        pp_result.secondaries.push_back(std::move(sec));
                    for (auto& e : ann.source_electrons)
                        pp_result.source_electrons.push_back(std::move(e));
                }
                return pp_result;
            }

            case InteractionType::Compton: {
                if (record_truth_labels_)
                    ++truth.compton_interactions;
                int Z_cs = seg.material->select_element(energy_MeV, 1 /*CS*/, uniform(rng));

                // Compton-angle mixture biasing (primary's first K vertices):
                // sample the outgoing direction from
                //   q(w) = (1-gamma) p(w) + gamma * 1[w in cone]/Omega_d,
                // p(w) = f(mu)/(2pi*N) the analog KN x S density on the
                // sphere, cone subtending the detector bounding sphere from
                // this vertex; multiply the event weight by p/q <= 1/(1-gamma).
                // Energy/recoil at the sampled angle go through the exact
                // analog kinematics (finish_compton_at_angle).
                ComptonScatterResult scatter;
                bool biased_vertex = false;
                if (bias_budget > 0 && det_bound_radius_ > 0.0) {
                    const double gamma = compton_bias_.cone_fraction;
                    const Eigen::Vector3d det_c(
                        0.0, 0.0, 0.5 * (det_bound_z_min_ + det_bound_z_max_));
                    const double half_h =
                        0.5 * (det_bound_z_max_ - det_bound_z_min_);
                    const double Rb = 1.05 * std::sqrt(
                        det_bound_radius_ * det_bound_radius_ + half_h * half_h);
                    Eigen::Vector3d axis = det_c - pos;
                    const double dist = axis.norm();
                    if (dist > Rb) {
                        axis /= dist;
                        const double sin_half = Rb / dist;
                        const double cos_half =
                            std::sqrt(1.0 - sin_half * sin_half);
                        const double omega_d = kTwoPi * (1.0 - cos_half);

                        Eigen::Vector3d w_dir;
                        if (uniform(rng) < gamma) {
                            double ct = 1.0 - uniform(rng) * (1.0 - cos_half);
                            double ph = kTwoPi * uniform(rng);
                            w_dir = rotate_direction(axis, ct, ph);
                        } else {
                            double mu_a =
                                sample_compton_cos_theta(energy, Z_cs, rng);
                            double ph = kTwoPi * uniform(rng);
                            w_dir = rotate_direction(dir, mu_a, ph);
                        }
                        const double mu = std::max(
                            -1.0, std::min(1.0, w_dir.dot(dir)));
                        const double p_sphere =
                            compton_angular_pdf_unnorm(energy, Z_cs, mu) /
                            (kTwoPi * compton_angular_norm(energy, Z_cs));
                        const double c_sphere =
                            (w_dir.dot(axis) >= cos_half) ? 1.0 / omega_d : 0.0;
                        const double q =
                            (1.0 - gamma) * p_sphere + gamma * c_sphere;
                        if (q > 0.0) {
                            w_bias *= p_sphere / q;
                            scatter = finish_compton_at_angle(
                                energy, dir, mu, rng, Z_cs,
                                compton_doppler_enabled_, &w_dir);
                            biased_vertex = true;
                            --bias_budget;
                        }
                    }
                }
                if (!biased_vertex) {
                    scatter = sample_compton_scatter(
                        energy, dir, rng, Z_cs,
                        compton_doppler_enabled_);
                }

                // A Doppler-sampled collision leaves a vacancy in a known
                // subshell (scatter.shell_index, xraylib order K, L1, L2, L3):
                // relax K and retained-line L vacancies as after
                // photoabsorption. The recoil electron then carries
                // E - E' - U rather than E - E'.
                double recoil_T = energy - scatter.scattered_energy_keV;
                if (source_fluorescence_enabled_ && scatter.binding_deposit_keV > 0.0) {
                    recoil_T = scatter.deposited_energy_keV;
                    const auto& xs_data = CrossSectionData::instance();
                    const double U = scatter.binding_deposit_keV;
                    const FluorescenceData* fl = xs_data.fluorescence(Z_cs);
                    const LFluorescenceData* fl_L = xs_data.l_fluorescence(Z_cs);
                    if (scatter.shell_index == 0 && fl
                        && std::abs(U - fl->k_edge_keV) < 0.05 * fl->k_edge_keV) {
                        relax_k_vacancy(pos, *fl, fl_L);
                    } else if (scatter.shell_index >= 1 && scatter.shell_index <= 3
                               && has_l_lines(fl_L)) {
                        const int s_l = scatter.shell_index - 1;
                        const double edge = fl_L->sub[s_l].edge_keV;
                        if (std::abs(U - edge) < 0.05 * edge)
                            transport_relaxation_xray(
                                pos, sample_l_relaxation(*fl_L, s_l, rng).line_keV);
                    }
                }
                if (recoil_T > std::min(kMoliereBremsThreshold_keV,
                                        source_electron_threshold_keV_)
                    || source_impact_ionization_enabled_) {
                    Eigen::Vector3d e_dir = compton_recoil_direction(
                        dir, scatter.new_direction,
                        energy, scatter.scattered_energy_keV);
                    process_electron(
                        pos, e_dir, recoil_T, *seg.material,
                        SourceSecondaryPhoton::Origin::ComptonElectronBrems);
                }

                energy = scatter.scattered_energy_keV;
                dir = scatter.new_direction;
                if (energy < 1.0) {
                    return {false, pos, dir, energy, std::move(result_secondaries),
                            std::move(result_electrons), true, w_bias, truth};
                }
                scattered = true;
                break;
            }

            case InteractionType::Rayleigh: {
                if (record_truth_labels_)
                    ++truth.rayleigh_interactions;
                int Z_rs = seg.material->select_element(energy_MeV, 2 /*RS*/, uniform(rng));
                double cos_theta = sample_rayleigh_cos_theta(Z_rs, energy, rng);
                double phi = kTwoPi * uniform(rng);
                dir = rotate_direction(dir, cos_theta, phi);
                scattered = true;
                break;
            }
            }

            if (scattered) break;  // Re-trace from new position/direction
        }

        if (!scattered) {
            // Walked all segments without interaction — photon exited
            pos += dir * distance_along_ray;
            return {true, pos, dir, energy, std::move(result_secondaries),
                    std::move(result_electrons), any_interaction, w_bias, truth};
        }
    }

    // Hit max scatters — kill the photon
    return {false, pos, dir, energy, std::move(result_secondaries),
            std::move(result_electrons), true, w_bias, truth};
}

double SourceGeometry::outermost_extent_radius() const {
    double r = 0.0;

    if (shape_ == Shape::Point) {
        // Total shield thickness
        for (const auto& layer : shields()) {
            r += layer.scalar_thickness();
        }
    } else if (shape_ == Shape::Cylindrical) {
        // Bounding-sphere radius: shields grow radially by tx and axially by tz.
        double rr = cyl_radius_;
        double hz = cyl_half_length_;
        for (const auto& layer : shields()) {
            rr += layer.tx;
            hz += layer.tz;
        }
        r = std::sqrt(rr * rr + hz * hz);
    } else if (shape_ == Shape::Rectangular) {
        // Bounding-sphere radius: half-diagonal of the outermost box.
        Eigen::Vector3d hd = rect_half_dims_;
        for (const auto& layer : shields()) {
            hd += Eigen::Vector3d(layer.tx, layer.ty, layer.tz);
        }
        r = hd.norm();
    } else if (shape_ == Shape::Sphere) {
        r = sphere_radius_;
        for (const auto& layer : shields()) {
            r += layer.scalar_thickness();
        }
    } else if (shape_ == Shape::Marinelli) {
        r = marinelli_outer_r_;
        for (const auto& layer : shields()) {
            r += layer.scalar_thickness();
        }
    }

    return r;
}

} // namespace ceelo
