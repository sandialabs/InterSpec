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

#include "geometry/ConvexSeparation.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace ceelo {

ConvexSolid ConvexSolid::point(const Eigen::Vector3d& c) {
    ConvexSolid s;
    s.kind = Kind::Point;
    s.center = c;
    return s;
}

ConvexSolid ConvexSolid::sphere(const Eigen::Vector3d& c, const double r) {
    ConvexSolid s;
    s.kind = Kind::Sphere;
    s.center = c;
    s.radius = std::max(0.0, r);
    return s;
}

ConvexSolid ConvexSolid::cylinder(const Eigen::Vector3d& c, const Eigen::Vector3d& axis,
                                  const double r, const double half_length) {
    ConvexSolid s;
    s.kind = Kind::Cylinder;
    s.center = c;
    s.axis = axis.normalized();
    s.radius = std::max(0.0, r);
    s.half_length = std::max(0.0, half_length);
    return s;
}

ConvexSolid ConvexSolid::box(const Eigen::Vector3d& c, const Eigen::Vector3d& half_extents,
                             const Eigen::Matrix3d& rotation) {
    ConvexSolid s;
    s.kind = Kind::Box;
    s.center = c;
    s.half_extents = half_extents.cwiseMax(0.0);
    s.rotation = rotation;
    return s;
}

Eigen::Vector3d ConvexSolid::project(const Eigen::Vector3d& p) const {
    const Eigen::Vector3d rel = p - center;
    switch (kind) {
        case Kind::Point:
            return center;

        case Kind::Sphere: {
            const double n = rel.norm();
            return (n <= radius) ? p : Eigen::Vector3d(center + rel * (radius / n));
        }

        case Kind::Cylinder: {
            // A solid cylinder's nearest point clamps the axial coordinate and the radial
            // distance independently.
            const double t = rel.dot(axis);
            const Eigen::Vector3d radial = rel - t * axis;
            const double rn = radial.norm();
            const double tc = std::clamp(t, -half_length, half_length);
            const Eigen::Vector3d rc =
                (rn <= radius) ? radial : Eigen::Vector3d(radial * (radius / rn));
            return center + tc * axis + rc;
        }

        case Kind::Box: {
            const Eigen::Vector3d local = rotation.transpose() * rel;
            const Eigen::Vector3d clamped = local.cwiseMax(-half_extents).cwiseMin(half_extents);
            return center + rotation * clamped;
        }
    }
    return p;
}

ConvexSolid ConvexSolid::transformed(const Eigen::Matrix3d& R, const Eigen::Vector3d& t) const {
    ConvexSolid s = *this;
    s.center = R * center + t;
    s.axis = R * axis;
    s.rotation = R * rotation;
    return s;
}

namespace {

/// Alternating projections from b's point nearest a's centre.  Each pass can only shorten the
/// pair's distance, so the loop ends on whichever of convergence or the pass budget comes first.
double alternate(const ConvexSolid& a, const ConvexSolid& b, Eigen::Vector3d& pa,
                 Eigen::Vector3d& pb, const double tol, const int max_iter) {
    pb = b.project(a.center);
    pa = a.project(pb);
    for (int i = 0; i < max_iter; ++i) {
        const Eigen::Vector3d nb = b.project(pa);
        const Eigen::Vector3d na = a.project(nb);
        const double moved = (nb - pb).norm() + (na - pa).norm();
        pa = na;
        pb = nb;
        if (moved < tol) break;
    }
    return (pa - pb).norm();
}

}  // namespace

double convex_separation(const ConvexSolid& a, const ConvexSolid& b, Eigen::Vector3d* closest_a,
                         Eigen::Vector3d* closest_b, const double tol, const int max_iter) {
    // Where the solids are (nearly) tangent - a curved face against a flat one - the iteration
    //  creeps, and how far it gets depends on which solid it starts from: start from each and keep
    //  the closer pair.
    Eigen::Vector3d pa, pb, qa, qb;
    const double d_ab = alternate(a, b, pa, pb, tol, max_iter);
    const double d_ba = alternate(b, a, qb, qa, tol, max_iter);
    if (d_ba < d_ab) {
        pa = qa;
        pb = qb;
    }
    if (closest_a) *closest_a = pa;
    if (closest_b) *closest_b = pb;
    return std::min(d_ab, d_ba);
}

namespace {

/// Intersects [lo, hi] with the parameter interval where `a + t*b` lies in [-h, h]; false if empty.
bool clip_slab(const double a, const double b, const double h, double& lo, double& hi) {
    if (std::fabs(b) < 1.0e-300) return (std::fabs(a) <= h);
    double t0 = (-h - a) / b, t1 = (h - a) / b;
    if (t0 > t1) std::swap(t0, t1);
    lo = std::max(lo, t0);
    hi = std::min(hi, t1);
    return (lo <= hi);
}

}  // namespace

double convex_entry_distance(const ConvexSolid& solid, const Eigen::Vector3d& origin,
                             const Eigen::Vector3d& dir) {
    constexpr double inf = std::numeric_limits<double>::infinity();
    const Eigen::Vector3d o = origin - solid.center;
    double lo = -inf, hi = inf;

    switch (solid.kind) {
        case ConvexSolid::Kind::Point:
            return inf;

        case ConvexSolid::Kind::Sphere: {
            // |o + t d|^2 = R^2
            const double b = o.dot(dir), c = o.squaredNorm() - solid.radius * solid.radius;
            const double disc = b * b - c;
            if (disc < 0.0) return inf;
            lo = -b - std::sqrt(disc);
            hi = -b + std::sqrt(disc);
            break;
        }

        case ConvexSolid::Kind::Cylinder: {
            const double oa = o.dot(solid.axis), da = dir.dot(solid.axis);
            if (!clip_slab(oa, da, solid.half_length, lo, hi)) return inf;
            const Eigen::Vector3d orad = o - oa * solid.axis, drad = dir - da * solid.axis;
            const double A = drad.squaredNorm(), B = orad.dot(drad),
                         C = orad.squaredNorm() - solid.radius * solid.radius;
            if (A < 1.0e-300) {
                if (C > 0.0) return inf;
            } else {
                const double disc = B * B - A * C;
                if (disc < 0.0) return inf;
                const double sq = std::sqrt(disc);
                lo = std::max(lo, (-B - sq) / A);
                hi = std::min(hi, (-B + sq) / A);
            }
            break;
        }

        case ConvexSolid::Kind::Box: {
            const Eigen::Vector3d ol = solid.rotation.transpose() * o,
                                  dl = solid.rotation.transpose() * dir;
            for (int i = 0; i < 3; ++i)
                if (!clip_slab(ol[i], dl[i], solid.half_extents[i], lo, hi)) return inf;
            break;
        }
    }

    if (!(lo <= hi) || (hi < 0.0)) return inf;
    return std::max(0.0, lo);
}

}  // namespace ceelo
