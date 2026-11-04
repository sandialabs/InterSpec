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

#include "geometry/EnvironmentGeometry.h"

#include <Eigen/Geometry>   // for .cross(); Eigen/Core alone is not enough
#include <algorithm>
#include <cmath>
#include <limits>

namespace ceelo {

namespace {

/// Rays exactly parallel to a plane produce no crossing; this guards the
/// division. Also used to nudge midpoint sampling off a boundary.
constexpr double kParallelEps = 1e-12;

} // namespace

EnvironmentGeometry::EnvironmentGeometry(const Material& air)
    : air_(air) {}

EnvironmentGeometry EnvironmentGeometry::open_air(
    const Material& air,
    const Eigen::Vector3d& source,
    const Eigen::Vector3d& detector,
    double transverse_extent_cm,
    double longitudinal_margin_cm)
{
    EnvironmentGeometry env(air);
    const Eigen::Vector3d delta = detector - source;
    const double len = delta.norm();
    if (len <= 0.0 || transverse_extent_cm <= 0.0) return env;

    const double margin = (longitudinal_margin_cm > 0.0) ? longitudinal_margin_cm
                                                         : transverse_extent_cm;
    EnvBounds b;
    b.enabled = true;
    b.axis_point = source;
    b.axis_dir = delta / len;
    b.radius_cm = transverse_extent_cm;
    b.t_min_cm = -margin;
    b.t_max_cm = len + margin;
    env.set_bounds(b);
    return env;
}

EnvironmentGeometry EnvironmentGeometry::field(
    const Material& air,
    const Material& ground,
    const Eigen::Vector3d& source,
    double h_src_cm,
    double h_det_cm,
    double ground_thickness_cm)
{
    EnvironmentGeometry env(air);

    // The detector sits at the origin; the source at `source`. We need the plane
    // whose perpendicular distance is h_det from the origin and h_src from the
    // source. Write the plane as n.p + c = 0 with |n| = 1 and n pointing up (into
    // the air). Then the signed height of a point p is n.p + c, so:
    //     c = h_det                       (detector at the origin)
    //     n.source + c = h_src  =>  n.source = h_src - h_det
    //
    // Any n satisfying that constraint is a valid ground plane; the physical one
    // is the one that stays "as horizontal as possible" — i.e. n lies in the
    // plane spanned by the source direction and the detector axis, tilted just
    // enough to meet the height difference.
    const double len = source.norm();
    if (len <= 0.0) {
        // Degenerate: source at the detector. Fall back to a plane below the
        // detector axis, normal along +y (the frame's "up" for a nominal
        // horizontally-aimed detector).
        env.add_slab(Eigen::Vector3d(0.0, -h_det_cm, 0.0),
                     Eigen::Vector3d(0.0, 1.0, 0.0), ground_thickness_cm, ground);
        return env;
    }

    const Eigen::Vector3d u = source / len;          // detector -> source
    const double dh = h_src_cm - h_det_cm;           // required n.source
    // Decompose n = a*u + b*w, where w is a unit vector perpendicular to u
    // chosen as the "up" direction. n.source = a*len = dh  =>  a = dh/len.
    const double a = dh / len;

    // Pick w perpendicular to u, preferring the frame's y axis as "up" so a
    // level geometry (h_src == h_det) gives an exactly horizontal ground plane.
    Eigen::Vector3d up(0.0, 1.0, 0.0);
    if (std::abs(up.dot(u)) > 0.99) up = Eigen::Vector3d(1.0, 0.0, 0.0);
    Eigen::Vector3d w = up - u * up.dot(u);
    const double wn = w.norm();
    if (wn < kParallelEps) return env;  // cannot construct; leave slab-free
    w /= wn;

    const double b2 = 1.0 - a * a;
    if (b2 < -1e-12) {
        // |h_src - h_det| EXCEEDS the separation: geometrically impossible for
        // any plane. Leave the environment slab-free rather than silently
        // producing a bogus one.
        return env;
    }
    // b2 == 0 is NOT degenerate -- it is the perpendicular-wall case, where the
    // plane normal lies along the source->detector axis (an end wall of a room,
    // for which h_src - h_det is exactly the separation). Superposing per-surface
    // kernels for a room requires it, so clamp rather than reject.
    const Eigen::Vector3d n =
        (u * a + w * std::sqrt(std::max(0.0, b2))).normalized();

    // Plane point: the foot of the perpendicular from the detector (origin).
    const Eigen::Vector3d plane_point = -n * h_det_cm;
    env.add_slab(plane_point, n, ground_thickness_cm, ground);
    return env;
}

EnvironmentGeometry EnvironmentGeometry::room(
    const Material& air,
    const Material& surfaces,
    const Eigen::Vector3d& source,
    double h_src_cm,
    double h_det_cm,
    double height_cm,
    double half_width_cm,
    double back_cm,
    double front_cm,
    double wall_thickness_cm)
{
    // Start from the field construction: it already solves for the floor plane
    // that sits at h_det below the detector and h_src below the source, which is
    // the only non-trivial part of the geometry (the two heights may differ).
    EnvironmentGeometry env =
        EnvironmentGeometry::field(air, surfaces, source, h_src_cm, h_det_cm,
                                   wall_thickness_cm);
    if (env.num_slabs() != 1) return env;  // degenerate; field() already bailed

    // Survey frame: w = up (the floor normal), u = horizontal along the
    // source->detector line, v = the remaining horizontal direction.
    const Eigen::Vector3d w = env.slabs()[0].normal;
    const Eigen::Vector3d floor_pt = env.slabs()[0].point;

    Eigen::Vector3d u = source - source.dot(w) * w;  // horizontal part of the
    if (u.norm() < 1e-9) {                            // source direction
        u = Eigen::Vector3d(1, 0, 0) - w * w.dot(Eigen::Vector3d(1, 0, 0));
    }
    u.normalize();
    const Eigen::Vector3d v = w.cross(u).normalized();

    // Ceiling: height_cm above the floor, normal pointing back down into the air.
    env.add_slab(floor_pt + w * height_cm, -w, wall_thickness_cm, surfaces);

    // Side walls at +/- half_width along v.
    env.add_slab(floor_pt + v * half_width_cm, -v, wall_thickness_cm, surfaces);
    env.add_slab(floor_pt - v * half_width_cm, v, wall_thickness_cm, surfaces);

    // End walls. u points from the detector toward the source, so the "front"
    // wall is beyond the source and the "back" wall is behind the detector.
    const double u_det = -floor_pt.dot(u);            // detector's u coordinate
    const double u_src = (source - floor_pt).dot(u);  // source's u coordinate
    env.add_slab(floor_pt + u * (u_src + front_cm), -u, wall_thickness_cm, surfaces);
    env.add_slab(floor_pt + u * (u_det - back_cm), u, wall_thickness_cm, surfaces);

    return env;
}

void EnvironmentGeometry::add_slab(const Eigen::Vector3d& point,
                                   const Eigen::Vector3d& normal,
                                   double thickness_cm,
                                   const Material& material)
{
    slabs_.push_back(EnvSlab{point, normal.normalized(), thickness_cm, material});
}

EnvPointClass EnvironmentGeometry::classify(const Eigen::Vector3d& p) const
{
    if (bounds_.enabled) {
        const Eigen::Vector3d rel = p - bounds_.axis_point;
        const double t = rel.dot(bounds_.axis_dir);
        if (t < bounds_.t_min_cm || t > bounds_.t_max_cm)
            return {nullptr, -1};
        const double r2 = (rel - bounds_.axis_dir * t).squaredNorm();
        if (r2 > bounds_.radius_cm * bounds_.radius_cm)
            return {nullptr, -1};
    }

    // Behind any slab plane? Pick the slab penetrated least (largest s among
    // the negatives) — the surface the point is nearest to. See classify() doc.
    int best = -1;
    double best_s = -std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < slabs_.size(); ++k) {
        const double s = (p - slabs_[k].point).dot(slabs_[k].normal);
        if (s >= 0.0) continue;              // air side of this plane
        if (s > best_s) { best_s = s; best = static_cast<int>(k); }
    }
    if (best < 0) return {&air_, -1};        // air side of every plane

    // Inside the slab material only within its thickness; beyond that the ray
    // has passed clean through and left the problem.
    if (-best_s <= slabs_[best].thickness_cm)
        return {&slabs_[best].material, best};
    return {nullptr, -1};
}

void EnvironmentGeometry::gather_crossings(const Eigen::Vector3d& pos,
                                           const Eigen::Vector3d& dir,
                                           double max_distance_cm,
                                           std::vector<double>& ts) const
{
    ts.clear();
    ts.push_back(0.0);

    auto push = [&](double t) {
        if (t > 0.0 && t < max_distance_cm) ts.push_back(t);
    };

    // Slab planes: the inner (air-facing) plane and the far face at depth.
    for (const EnvSlab& s : slabs_) {
        const double denom = dir.dot(s.normal);
        if (std::abs(denom) < kParallelEps) continue;
        const double s0 = (pos - s.point).dot(s.normal);
        push(-s0 / denom);                              // inner plane, s = 0
        push((-s.thickness_cm - s0) / denom);           // far face, s = -thickness
    }

    // Bounding cylinder: the two end caps and the lateral surface.
    if (bounds_.enabled) {
        const Eigen::Vector3d rel = pos - bounds_.axis_point;
        const double vd = dir.dot(bounds_.axis_dir);
        const double pd = rel.dot(bounds_.axis_dir);
        if (std::abs(vd) > kParallelEps) {
            push((bounds_.t_min_cm - pd) / vd);
            push((bounds_.t_max_cm - pd) / vd);
        }
        // Radial: |rel_perp + t*dir_perp|^2 = R^2
        const Eigen::Vector3d rp = rel - bounds_.axis_dir * pd;
        const Eigen::Vector3d dp = dir - bounds_.axis_dir * vd;
        const double A = dp.squaredNorm();
        if (A > kParallelEps) {
            const double B = 2.0 * rp.dot(dp);
            const double C = rp.squaredNorm() - bounds_.radius_cm * bounds_.radius_cm;
            const double disc = B * B - 4.0 * A * C;
            if (disc >= 0.0) {
                const double sq = std::sqrt(disc);
                push((-B - sq) / (2.0 * A));
                push((-B + sq) / (2.0 * A));
            }
        }
    }

    ts.push_back(max_distance_cm);
    std::sort(ts.begin(), ts.end());
}

void EnvironmentGeometry::trace(const Eigen::Vector3d& pos,
                                const Eigen::Vector3d& dir,
                                std::vector<EnvSegment>& out,
                                double max_distance_cm) const
{
    EnvTraceScratch scratch;
    trace(pos, dir, out, scratch, max_distance_cm);
}

void EnvironmentGeometry::trace(const Eigen::Vector3d& pos,
                                const Eigen::Vector3d& dir,
                                std::vector<EnvSegment>& out,
                                EnvTraceScratch& scratch,
                                double max_distance_cm) const
{
    out.clear();
    std::vector<double>& ts = scratch.crossings;
    gather_crossings(pos, dir, max_distance_cm, ts);

    for (std::size_t i = 0; i + 1 < ts.size(); ++i) {
        const double t0 = ts[i], t1 = ts[i + 1];
        if (t1 - t0 <= kParallelEps) continue;
        // Classify at the midpoint: unambiguous, and immune to landing exactly
        // on a boundary.
        const EnvPointClass c = classify(pos + dir * (0.5 * (t0 + t1)));
        if (c.material == nullptr) break;  // entered void: the ray is gone
        if (!out.empty() && out.back().slab_index == c.slab_index &&
            out.back().material == c.material) {
            out.back().t_end_cm = t1;      // merge adjacent identical regions
        } else {
            out.push_back({t0, t1, c.material, c.slab_index});
        }
    }
}

double EnvironmentGeometry::optical_depth(const Eigen::Vector3d& from,
                                          const Eigen::Vector3d& to,
                                          double energy_MeV) const
{
    EnvTraceScratch scratch;
    return optical_depth(from, to, energy_MeV, scratch);
}

double EnvironmentGeometry::optical_depth(const Eigen::Vector3d& from,
                                          const Eigen::Vector3d& to,
                                          double energy_MeV,
                                          EnvTraceScratch& scratch) const
{
    const Eigen::Vector3d delta = to - from;
    const double dist = delta.norm();
    if (dist <= 0.0) return 0.0;
    const Eigen::Vector3d dir = delta / dist;

    std::vector<double>& ts = scratch.crossings;
    gather_crossings(from, dir, dist, ts);

    // Cache mu by material pointer: a chord crosses at most a handful of
    // distinct materials (air plus one or two slabs), and macroscopic_xs()
    // re-interpolates the cross-section tables on every call.
    const Material* cached_mat = nullptr;
    double cached_mu = 0.0;

    double tau = 0.0;
    for (std::size_t i = 0; i + 1 < ts.size(); ++i) {
        const double t0 = ts[i], t1 = ts[i + 1];
        if (t1 - t0 <= kParallelEps) continue;
        const EnvPointClass c = classify(from + dir * (0.5 * (t0 + t1)));
        if (c.material == nullptr) {
            // The chord leaves the modeled world. Anything beyond contributes
            // no attenuation (there is no material out there).
            continue;
        }
        if (c.material != cached_mat) {
            cached_mat = c.material;
            cached_mu = c.material->mu_total(energy_MeV);
        }
        tau += cached_mu * (t1 - t0);
    }
    return tau;
}

double EnvironmentGeometry::transmission(const Eigen::Vector3d& from,
                                         const Eigen::Vector3d& to,
                                         double energy_MeV) const
{
    return std::exp(-optical_depth(from, to, energy_MeV));
}

double EnvironmentGeometry::transmission(const Eigen::Vector3d& from,
                                         const Eigen::Vector3d& to,
                                         double energy_MeV,
                                         EnvTraceScratch& scratch) const
{
    return std::exp(-optical_depth(from, to, energy_MeV, scratch));
}

std::array<double, 2> EnvironmentGeometry::transmission_pair(
    const Eigen::Vector3d& from,
    const Eigen::Vector3d& to,
    double first_energy_MeV,
    double second_energy_MeV,
    EnvTraceScratch& scratch) const
{
    const Eigen::Vector3d delta = to - from;
    const double dist = delta.norm();
    if (dist <= 0.0) return {{1.0, 1.0}};
    const Eigen::Vector3d dir = delta / dist;

    std::vector<double>& ts = scratch.crossings;
    gather_crossings(from, dir, dist, ts);

    const Material* cached_mat = nullptr;
    double cached_first_mu = 0.0;
    double cached_second_mu = 0.0;
    double first_tau = 0.0;
    double second_tau = 0.0;
    for (std::size_t i = 0; i + 1 < ts.size(); ++i) {
        const double t0 = ts[i], t1 = ts[i + 1];
        if (t1 - t0 <= kParallelEps) continue;
        const EnvPointClass c = classify(from + dir * (0.5 * (t0 + t1)));
        if (c.material == nullptr) continue;
        if (c.material != cached_mat) {
            cached_mat = c.material;
            cached_first_mu = c.material->mu_total(first_energy_MeV);
            cached_second_mu = first_energy_MeV == second_energy_MeV
                                   ? cached_first_mu
                                   : c.material->mu_total(second_energy_MeV);
        }
        const double length = t1 - t0;
        first_tau += cached_first_mu * length;
        second_tau += cached_second_mu * length;
    }
    return {{std::exp(-first_tau), std::exp(-second_tau)}};
}

} // namespace ceelo
