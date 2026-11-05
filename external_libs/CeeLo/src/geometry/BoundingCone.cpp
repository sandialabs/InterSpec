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

/// @file BoundingCone.cpp
/// @brief See BoundingCone.h.

#include "geometry/BoundingCone.h"

#include <Eigen/Geometry>  // .cross()

#include <algorithm>
#include <cmath>
#include <vector>

namespace ceelo {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Smallest cos(angle) between `axis` and the direction from the source to any point of the
/// circle of radius R centred at (0,0,z_k). Exact for an IN-PLANE axis (see the caller).
///
/// With the source rotated into the x-z plane at (s_r, 0, s_z) and the axis in that same plane
/// as (a_x, 0, a_z), writing u = cos(phi) over the circle gives
///     cos(angle)(u) = (A u + B) / sqrt(P - Q u)
/// whose only interior stationary point solves A*(P - Q u) + (Q/2)*(A u + B) = 0. So the extremum
/// is that point or an endpoint -- three evaluations, no iteration.
double min_cos_to_circle(double s_r, double s_z, double a_x, double a_z, double R, double z_k) {
    const double dz = z_k - s_z;
    const double A = R * a_x;
    const double B = -a_x * s_r + a_z * dz;
    const double P = R * R + s_r * s_r + dz * dz;
    const double Q = 2.0 * R * s_r;

    const auto f = [&](double u) {
        const double D = P - Q * u;
        if (D <= 1e-300) return -1.0;   // source sits ON the circle
        return (A * u + B) / std::sqrt(D);
    };

    double best = std::min(f(-1.0), f(1.0));
    if (std::abs(A) > 0.0 && std::abs(Q) > 0.0) {
        const double u_star = (2.0 * A * P + Q * B) / (A * Q);
        if (u_star > -1.0 && u_star < 1.0) best = std::min(best, f(u_star));
    }
    return best;
}

/// Smallest cos(angle) between `axis` and the direction from `src` to any point of the box
/// [-hx,hx] x [-hy,hy] x [z0,z1], evaluated at the 8 corners.
///
/// VALID ONLY WHEN THE RESULT IS >= 0, and the caller must enforce that. The corner argument needs
/// angular distance to `axis` to be convex along geodesics, which holds only inside the hemisphere
/// about `axis`: the set {x : a.x >= c|x|} is a convex (second-order) cone for c >= 0 and NOT
/// convex for c < 0. With c < 0 the extreme over the body is not attained at an extreme point --
/// e.g. for a source just off the middle of a face, the direction to the nearest face point is not
/// in the span of the corner directions at all -- and the returned value is an UNDER-estimate of
/// the true half-angle, i.e. a cone that silently clips the crystal.
double min_cos_to_box(const Eigen::Vector3d& src, const Eigen::Vector3d& axis,
                      double hx, double hy, double z0, double z1) {
    double best = 1.0;
    for (int i = 0; i < 8; ++i) {
        const Eigen::Vector3d p((i & 1) ? hx : -hx, (i & 2) ? hy : -hy, (i & 4) ? z1 : z0);
        const Eigen::Vector3d d = p - src;
        const double n = d.norm();
        if (n <= 1e-12) return -1.0;    // source at a corner
        best = std::min(best, axis.dot(d) / n);
    }
    return best;
}

}  // namespace

/// Search notes for envelope_bounding_cone.
///
/// SAFETY: both extremum routines above are exact only where the resulting cos is >= 0, so the
/// search here accepts a candidate axis ONLY while its computed worst-case cos stays positive, and
/// the function refuses (returns false) rather than hand back a cone it cannot vouch for. That is
/// not a corner case to wave at:
///   - For a box the plane-restricted search below cannot always reach the optimal axis (the
///     inward normal of a +/-y face is out of that plane when the source is offset in both x and
///     y), and an unguarded search settled on a NEGATIVE cos that the corner evaluator then
///     flattered -- measured up to 43% of the active solid angle sampled away, silently.
///   - For a cylinder, `min_cos_to_circle` sees only the in-plane axis components, so an
///     out-of-plane axis is evaluated as its in-plane projection: worst_cos scales by
///     sqrt(1 - a_t^2), which SHRINKS a positive value (so such a step is rejected) but GROWS a
///     negative one, letting an out-of-plane step always "improve" and run away. Hence both the
///     positivity guard and the in-plane-only refinement for cylinders.
/// A source outside a convex body always admits a cone under 90 degrees (separating hyperplane),
/// so a positive cos is not a restriction in exact arithmetic -- but it is not robust to round-off
/// on its own: on the back-face plane sin(pi) evaluates to 1.2e-16 rather than 0, which is enough
/// to make the true 90-degree optimum come out a hair negative and trip the runaway above.
bool envelope_bounding_cone(const bool cylinder, const double R, const double hx, const double hy,
                            const double z0, const double z1, const Eigen::Vector3d& src,
                            Eigen::Vector3d& axis_out, double& cos_alpha_out) {
    // Refuse anything not comfortably inside the evaluators' valid domain. A cone at 90 degrees is
    // worth at most a 2x yield gain over the full sphere, so declining near there costs ~nothing.
    constexpr double kMinCos = 1e-6;

    if (!(z1 > z0)) return false;

    const double s_r = std::hypot(src.x(), src.y());

    // Inside the envelope => the body subtends more than any cone can hold.
    if (src.z() > z0 && src.z() < z1) {
        if (cylinder ? (s_r < R) : (std::abs(src.x()) < hx && std::abs(src.y()) < hy))
            return false;
    }

    const Eigen::Vector3d e_r = (s_r > 1e-9)
        ? Eigen::Vector3d(src.x() / s_r, src.y() / s_r, 0.0) : Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d e_z(0.0, 0.0, 1.0);

    const auto worst_cos = [&](const Eigen::Vector3d& axis) {
        if (cylinder) {
            const double a_x = axis.dot(e_r), a_z = axis.z();
            return std::min(min_cos_to_circle(s_r, src.z(), a_x, a_z, R, z0),
                            min_cos_to_circle(s_r, src.z(), a_x, a_z, R, z1));
        }
        return min_cos_to_box(src, axis, hx, hy, z0, z1);
    };

    double best_cos = -2.0;
    double best_psi = 0.0;
    Eigen::Vector3d best_axis = e_z;

    // A cylinder is axisymmetric, so its optimal axis lies in the plane spanned by the detector
    // axis and the source's radial offset: scanning that plane suffices.
    for (int i = 0; i < 720; ++i) {
        const double psi = 2.0 * kPi * i / 720.0;
        const Eigen::Vector3d axis = std::cos(psi) * e_z + std::sin(psi) * e_r;
        const double c = worst_cos(axis);
        if (c > best_cos) { best_cos = c; best_axis = axis; best_psi = psi; }
    }

    if (cylinder) {
        // In-plane refinement ONLY -- see the safety note above.
        double step = 2.0 * kPi / 720.0;
        for (int iter = 0; iter < 60 && step > 1e-12; ++iter) {
            bool improved = false;
            for (const double sgn : {+1.0, -1.0}) {
                const double psi = best_psi + sgn * step;
                const Eigen::Vector3d cand = std::cos(psi) * e_z + std::sin(psi) * e_r;
                const double c = worst_cos(cand);
                if (c > best_cos) { best_cos = c; best_axis = cand; best_psi = psi; improved = true; }
            }
            if (!improved) step *= 0.5;
        }
    } else {
        // A box is not axisymmetric, so the in-plane scan above is only one seed among several.
        // The optimum can sit anywhere on the sphere -- for a source just off the middle of a
        // large face it is essentially that face's NORMAL, which is neither in the scan plane nor
        // along any corner direction -- so seed globally and cheaply. A few thousand corner
        // evaluations cost nothing next to the n_rays full-geometry ray traces that follow.
        const double zc = 0.5 * (z0 + z1);
        std::vector<Eigen::Vector3d> seeds = {
            {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0},      // face normals
            {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0},
            {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0},
        };
        for (const Eigen::Vector3d& p : {                    // face centres and corners
                 Eigen::Vector3d{0.0, 0.0, zc},
                 { hx, 0.0, zc}, {-hx, 0.0, zc}, {0.0,  hy, zc}, {0.0, -hy, zc},
                 {0.0, 0.0, z0}, {0.0, 0.0, z1},
                 { hx,  hy, z0}, { hx, -hy, z0}, {-hx,  hy, z0}, {-hx, -hy, z0},
                 { hx,  hy, z1}, { hx, -hy, z1}, {-hx,  hy, z1}, {-hx, -hy, z1}}) {
            const Eigen::Vector3d d = p - src;
            if (d.norm() > 1e-12) seeds.push_back(d.normalized());
        }
        constexpr int kBoxSeedDirs = 1024;                   // global Fibonacci sweep
        const double golden = kPi * (3.0 - std::sqrt(5.0));
        for (int i = 0; i < kBoxSeedDirs; ++i) {
            const double ct = 1.0 - 2.0 * (i + 0.5) / kBoxSeedDirs;
            const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
            const double ph = golden * i;
            seeds.emplace_back(st * std::cos(ph), st * std::sin(ph), ct);
        }
        for (const Eigen::Vector3d& cand : seeds) {
            const double c = worst_cos(cand);
            if (c > best_cos) { best_cos = c; best_axis = cand; }
        }

        double step = 0.25;
        for (int iter = 0; iter < 200 && step > 1e-12; ++iter) {
            const Eigen::Vector3d t1 = best_axis.unitOrthogonal();
            const Eigen::Vector3d t2 = best_axis.cross(t1).normalized();
            bool improved = false;
            for (const Eigen::Vector3d& t : {t1, t2}) {
                for (const double sgn : {+1.0, -1.0}) {
                    const Eigen::Vector3d cand = (best_axis + sgn * step * t).normalized();
                    const double c = worst_cos(cand);
                    if (c > best_cos) { best_cos = c; best_axis = cand; improved = true; }
                }
            }
            if (!improved) step *= 0.5;
        }
    }

    // Only vouch for a cone the evaluators were entitled to measure.
    if (!(best_cos > kMinCos)) return false;

    // Pad by a hair against round-off in the exact extremum; a cone a whisker too wide costs a
    // little efficiency, one a whisker too narrow silently loses paths.
    const double cos_alpha = best_cos - 1e-9;
    if (!(cos_alpha > 0.0) || 0.5 * (1.0 - cos_alpha) >= 1.0) return false;

    axis_out = best_axis;
    cos_alpha_out = cos_alpha;
    return true;
}

}  // namespace ceelo
