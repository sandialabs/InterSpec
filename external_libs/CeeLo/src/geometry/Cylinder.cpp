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

#include "geometry/Cylinder.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace ceelo {

namespace {

constexpr double kEpsilon = 1.0e-12;
constexpr double kInfinity = std::numeric_limits<double>::max();

/// Clip a ray interval [t_enter, t_exit] (currently intersecting an infinite cylinder)
/// against the z-slab [z_min, z_max].  Also checks endcap intersections.
/// Returns false if the clipped interval is empty.
bool clip_to_z_range(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double z_min, double z_max,
    double& t_enter, double& t_exit)
{
    double oz = origin.z();
    double dz = direction.z();

    if (std::abs(dz) < kEpsilon) {
        // Ray is parallel to endcaps
        if (oz < z_min || oz > z_max) {
            return false; // Outside the z-range entirely
        }
        // Inside z-range, t_enter/t_exit unchanged
        return t_enter < t_exit;
    }

    // Intersect with z=z_min and z=z_max planes
    double t_z_lo = (z_min - oz) / dz;
    double t_z_hi = (z_max - oz) / dz;
    if (t_z_lo > t_z_hi) std::swap(t_z_lo, t_z_hi);

    // Intersect the cylinder interval with the z-slab interval
    t_enter = std::max(t_enter, t_z_lo);
    t_exit  = std::min(t_exit,  t_z_hi);

    return t_enter < t_exit;
}

} // anonymous namespace

std::optional<RayHit> intersect_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double radius,
    double z_min,
    double z_max)
{
    // Solve the quadratic for the infinite cylinder x² + y² = R²
    // Ray: P(t) = O + t*D
    // (Ox + t*Dx)² + (Oy + t*Dy)² = R²
    // a*t² + b*t + c = 0
    double ox = origin.x(), oy = origin.y();
    double dx = direction.x(), dy = direction.y();

    double a = dx * dx + dy * dy;
    double b = 2.0 * (ox * dx + oy * dy);
    double c = ox * ox + oy * oy - radius * radius;

    double t_enter, t_exit;

    if (a < kEpsilon) {
        // Ray is parallel to the cylinder axis (along z)
        if (c > kEpsilon) {
            return std::nullopt; // Outside the cylinder, moving parallel
        }
        // Inside (or on surface of) the cylinder
        t_enter = -kInfinity;
        t_exit  =  kInfinity;
    } else {
        double discriminant = b * b - 4.0 * a * c;
        if (discriminant < 0.0) {
            return std::nullopt; // No intersection with infinite cylinder
        }

        double sqrt_disc = std::sqrt(discriminant);
        double inv_2a = 0.5 / a;
        t_enter = (-b - sqrt_disc) * inv_2a;
        t_exit  = (-b + sqrt_disc) * inv_2a;
    }

    // Clip to the finite z-range
    if (!clip_to_z_range(origin, direction, z_min, z_max, t_enter, t_exit)) {
        return std::nullopt;
    }

    // Must have some forward intersection
    if (t_exit <= 0.0) {
        return std::nullopt;
    }

    return RayHit{t_enter, t_exit};
}


namespace {

/// Subtract the bore interval from the outer-crystal interval, yielding the
/// 0-2 active segments.  Factored out of intersect_bored_cylinder() so the
/// bulletized and rounded-tip variants reuse the exact same arithmetic.
int subtract_bore_interval(const RayHit& outer_hit,
                           const std::optional<RayHit>& bore_hit,
                           RayHit segments_out[2])
{
    if (!bore_hit || !bore_hit->valid()) {
        // Ray doesn't hit the bore — single continuous segment through crystal
        segments_out[0] = outer_hit;
        return 1;
    }

    // Subtract the bore from the outer crystal.
    // The outer crystal interval is [t_enter_out, t_exit_out].
    // The bore interval is [t_enter_bore, t_exit_bore].
    // Active segments = outer_interval minus bore_interval.

    double t_eo = std::max(outer_hit.t_enter, 0.0);
    double t_xo = outer_hit.t_exit;
    double t_eb = bore_hit->t_enter;
    double t_xb = bore_hit->t_exit;

    // Clamp bore interval to outer interval
    t_eb = std::max(t_eb, t_eo);
    t_xb = std::min(t_xb, t_xo);

    if (t_eb >= t_xb) {
        // Bore interval is empty within the outer crystal — single segment
        segments_out[0] = outer_hit;
        return 1;
    }

    // We have bore [t_eb, t_xb] inside outer [t_eo, t_xo].
    // Active segments are: [t_eo, t_eb] and [t_xb, t_xo]
    int count = 0;

    if (t_eo < t_eb - kEpsilon) {
        segments_out[count].t_enter = outer_hit.t_enter;
        segments_out[count].t_exit  = t_eb;
        count++;
    }

    if (t_xb < t_xo - kEpsilon) {
        segments_out[count].t_enter = t_xb;
        segments_out[count].t_exit  = t_xo;
        count++;
    }

    return count;
}

} // anonymous namespace

int intersect_bored_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double outer_radius,
    double bore_radius,
    double z_min,
    double z_max,
    double bore_z_start,
    double bore_z_end,
    RayHit segments_out[2])
{
    // Step 1: Intersect ray with the outer cylinder
    auto outer_hit = intersect_cylinder(origin, direction, outer_radius, z_min, z_max);
    if (!outer_hit || !outer_hit->valid()) {
        return 0; // Ray misses the outer crystal entirely
    }

    // Step 2: Intersect ray with the bore cylinder (which goes from bore_z_start to bore_z_end)
    auto bore_hit = intersect_cylinder(origin, direction, bore_radius, bore_z_start, bore_z_end);

    // Step 3: Subtract the bore from the outer crystal.
    return subtract_bore_interval(*outer_hit, bore_hit, segments_out);
}


// ---- Bulletized (rounded) front outer edge ----
//
// Notation used throughout this section, for a ray P(t) = O + t*D with D
// normalized, and rho(t) = sqrt(x(t)² + y(t)²):
//
//     q(t)  = a t² + b t + c   with a, b, c the usual infinite-cylinder
//                              quadratic coefficients, so rho = sqrt(q)
//     G(t)  = (rho(t) - rho_c)² + (z(t) - z_c)² - r_b²
//
// G is the signed "outside-ness" of the fillet torus: G < 0 inside the tube,
// G > 0 outside.  Setting G = 0 and clearing the square root gives the usual
// ray-torus quartic
//
//     t⁴ + 4β t³ + (4β² + 2γ - 4ρ_c² a) t² + (4βγ - 4ρ_c² b) t + γ² - 4ρ_c² c
//
// with β = O'·D, γ = O'·O' + ρ_c² - r_b², O' = O - (0,0,z_c).  We deliberately
// do NOT solve that quartic in closed form: its coefficients cancel
// catastrophically for exactly the rays we care about (grazing corner chords,
// i.e. near-double roots), and origins sit tens of cm from a sub-cm feature.
//
// Instead we exploit two structural facts:
//
//   1. G is CONVEX wherever rho >= rho_c.  rho(t) is a norm of an affine
//      function (convex), x -> (x - rho_c)² is convex and non-decreasing for
//      x >= rho_c, and (z - z_c)² is convex.  The corner zone requires
//      rho > rho_c, so G is convex on every corner-zone sub-segment.
//   2. Brackets are guaranteed.  Every point where the ray leaves the corner
//      zone while still inside the sharp cylinder lies on rho = rho_c with
//      z in [z_c - r_b, z_c], or on z = z_c with rho <= rho_c + r_b; both are
//      within r_b of the ring centre, i.e. G <= 0 there.  Those crossings are
//      closed form (a linear and a quadratic equation), so we can always
//      produce a sign-changing bracket, and convexity makes the root in it
//      unique.
//
// The result is a safeguarded Newton iteration that converges in a handful of
// steps and cannot land on the wrong root.

namespace {

/// G(t) evaluated at a point.
inline double fillet_G(const Eigen::Vector3d& p, const FrontFillet& f) {
    const double rho = std::sqrt(p.x() * p.x() + p.y() * p.y());
    const double dr  = rho - f.rho_c;
    const double dz  = p.z() - f.z_c;
    return dr * dr + dz * dz - f.r_b_sq;
}

inline double fillet_G_at(const Eigen::Vector3d& origin,
                          const Eigen::Vector3d& direction,
                          const FrontFillet& f, double t) {
    return fillet_G(origin + direction * t, f);
}

/// dG/dt.  q'(t) = 2(a t) + b; d/dt (rho - rho_c)² = q' (1 - rho_c/rho).
/// rho >= rho_c > 0 on every interval this is used, so there is no division
/// hazard.
inline double fillet_dG_at(const Eigen::Vector3d& origin,
                           const Eigen::Vector3d& direction,
                           const FrontFillet& f, double t) {
    const Eigen::Vector3d p = origin + direction * t;
    const double rho = std::sqrt(p.x() * p.x() + p.y() * p.y());
    const double dq  = 2.0 * (p.x() * direction.x() + p.y() * direction.y());
    const double drho_term = (rho > 0.0) ? dq * (1.0 - f.rho_c / rho) : 0.0;
    return drho_term + 2.0 * direction.z() * (p.z() - f.z_c);
}

/// True when `p` is in the material the bulletization removes: inside the
/// corner zone {rho > rho_c, z < z_c} and outside the fillet torus.  All three
/// tests are strict, so a point exactly on a surface the sharp and bulletized
/// solids share is left untouched (and therefore stays exact).
inline bool in_removed_corner(const Eigen::Vector3d& p, const FrontFillet& f) {
    if (p.z() >= f.z_c) return false;
    const double rho_sq = p.x() * p.x() + p.y() * p.y();
    if (rho_sq <= f.rho_c * f.rho_c) return false;
    return fillet_G(p, f) > 0.0;
}

/// First t > `t_from` (capped at `t_cap`) at which the ray leaves the corner
/// zone, by rising to z = z_c or by falling to rho = rho_c.
double corner_exit_forward(const Eigen::Vector3d& origin,
                           const Eigen::Vector3d& direction,
                           const FrontFillet& f, double t_from, double t_cap) {
    double t_stop = t_cap;

    if (direction.z() > kEpsilon) {
        const double tz = (f.z_c - origin.z()) / direction.z();
        if (tz > t_from && tz < t_stop) t_stop = tz;
    }

    const double a = direction.x() * direction.x() + direction.y() * direction.y();
    if (a > kEpsilon) {
        const double b = 2.0 * (origin.x() * direction.x() + origin.y() * direction.y());
        const double c = origin.x() * origin.x() + origin.y() * origin.y()
                       - f.rho_c * f.rho_c;
        const double disc = b * b - 4.0 * a * c;
        if (disc > 0.0) {
            // Entry into {rho <= rho_c}; the ray is outside it at t_from.
            const double s0 = (-b - std::sqrt(disc)) * (0.5 / a);
            if (s0 > t_from && s0 < t_stop) t_stop = s0;
        }
    }

    return t_stop;
}

/// Mirror of corner_exit_forward() walking backwards from `t_from` (floored
/// at `t_floor`).
double corner_exit_backward(const Eigen::Vector3d& origin,
                            const Eigen::Vector3d& direction,
                            const FrontFillet& f, double t_from, double t_floor) {
    double t_stop = t_floor;

    if (direction.z() < -kEpsilon) {
        const double tz = (f.z_c - origin.z()) / direction.z();
        if (tz < t_from && tz > t_stop) t_stop = tz;
    }

    const double a = direction.x() * direction.x() + direction.y() * direction.y();
    if (a > kEpsilon) {
        const double b = 2.0 * (origin.x() * direction.x() + origin.y() * direction.y());
        const double c = origin.x() * origin.x() + origin.y() * origin.y()
                       - f.rho_c * f.rho_c;
        const double disc = b * b - 4.0 * a * c;
        if (disc > 0.0) {
            const double s1 = (-b + std::sqrt(disc)) * (0.5 / a);
            if (s1 < t_from && s1 > t_stop) t_stop = s1;
        }
    }

    return t_stop;
}

/// Root of G on a bracket with G(t_out) > 0 >= G(t_in).  Safeguarded Newton:
/// Newton where it stays in the bracket, bisection otherwise.
double solve_fillet_root(const Eigen::Vector3d& origin,
                         const Eigen::Vector3d& direction,
                         const FrontFillet& f, double t_out, double t_in) {
    double x_pos = t_out;  // G > 0
    double x_neg = t_in;   // G <= 0
    double t = 0.5 * (x_pos + x_neg);

    for (int iter = 0; iter < 80; ++iter) {
        const double g = fillet_G_at(origin, direction, f, t);
        if (g > 0.0) x_pos = t; else x_neg = t;

        const double lo = std::min(x_neg, x_pos);
        const double hi = std::max(x_neg, x_pos);

        const double dg = fillet_dG_at(origin, direction, f, t);
        double t_next = (dg != 0.0) ? (t - g / dg) : 0.5 * (lo + hi);
        if (!(t_next > lo && t_next < hi)) {
            t_next = 0.5 * (lo + hi);
        }

        const double step = std::abs(t_next - t);
        t = t_next;
        if (step <= 1.0e-14 * (1.0 + std::abs(t))) break;
    }

    return t;
}

/// Minimiser of the convex G on [lo, hi], found by bisecting its monotone
/// derivative.  Only reached for chords that stay inside the corner zone from
/// end to end -- a grazing corner clip, or a complete miss.
double minimize_fillet_G(const Eigen::Vector3d& origin,
                         const Eigen::Vector3d& direction,
                         const FrontFillet& f, double lo, double hi) {
    if (fillet_dG_at(origin, direction, f, lo) >= 0.0) return lo;
    if (fillet_dG_at(origin, direction, f, hi) <= 0.0) return hi;

    for (int iter = 0; iter < 80; ++iter) {
        const double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) break;  // exhausted double precision
        if (fillet_dG_at(origin, direction, f, mid) > 0.0) hi = mid; else lo = mid;
    }

    return 0.5 * (lo + hi);
}

} // anonymous namespace

bool correct_hit_for_front_fillet(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    const FrontFillet& fillet,
    RayHit& hit)
{
    const bool fix_enter =
        in_removed_corner(origin + direction * hit.t_enter, fillet);
    const bool fix_exit =
        in_removed_corner(origin + direction * hit.t_exit, fillet);

    if (!fix_enter && !fix_exit) {
        return true;  // both endpoints lie on shared surfaces — already exact
    }

    // Each endpoint needing correction gets its own bracket end inside the
    // solid, taken at the point where the ray leaves the corner zone.
    double in_enter = hit.t_enter, in_exit = hit.t_exit;
    bool have_enter = false, have_exit = false;

    if (fix_enter) {
        const double t = corner_exit_forward(origin, direction, fillet,
                                             hit.t_enter, hit.t_exit);
        if (fillet_G_at(origin, direction, fillet, t) <= 0.0) {
            in_enter = t;
            have_enter = true;
        }
    }
    if (fix_exit) {
        const double t = corner_exit_backward(origin, direction, fillet,
                                              hit.t_exit, hit.t_enter);
        if (fillet_G_at(origin, direction, fillet, t) <= 0.0) {
            in_exit = t;
            have_exit = true;
        }
    }

    if ((fix_enter && !have_enter) || (fix_exit && !have_exit)) {
        // The chord never leaves the corner zone: it either clips the fillet
        // as a grazing chord or misses the solid altogether.
        const double t_min = minimize_fillet_G(origin, direction, fillet,
                                               hit.t_enter, hit.t_exit);
        if (fillet_G_at(origin, direction, fillet, t_min) >= 0.0) {
            return false;  // passes through removed corner material only
        }
        if (fix_enter && !have_enter) { in_enter = t_min; have_enter = true; }
        if (fix_exit  && !have_exit)  { in_exit  = t_min; have_exit  = true; }
    }

    if (fix_enter) {
        hit.t_enter = solve_fillet_root(origin, direction, fillet,
                                        hit.t_enter, in_enter);
    }
    if (fix_exit) {
        hit.t_exit = solve_fillet_root(origin, direction, fillet,
                                       hit.t_exit, in_exit);
    }

    return hit.t_enter < hit.t_exit;
}


std::optional<RayHit> intersect_bulletized_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double radius,
    double z_min,
    double z_max,
    const FrontFillet& fillet)
{
    auto hit = intersect_cylinder(origin, direction, radius, z_min, z_max);
    if (!hit) return std::nullopt;

    if (!correct_hit_for_front_fillet(origin, direction, fillet, *hit)) {
        return std::nullopt;
    }
    if (hit->t_exit <= 0.0) {
        return std::nullopt;
    }

    return hit;
}


namespace {

/// Ray interval through a bore with a hemispherical ("round-tipped drill")
/// closed end: the straight part from z_tip upward, unioned with the ball
/// tangent to it.  The union is convex, and the ball's upper half lies inside
/// the straight part, so the interval is simply the hull of the two pieces.
std::optional<RayHit> rounded_bore_interval(const Eigen::Vector3d& origin,
                                            const Eigen::Vector3d& direction,
                                            double bore_radius,
                                            double bore_z_start,
                                            double bore_z_end)
{
    const double z_tip = bore_z_start + bore_radius;

    auto straight = intersect_cylinder(origin, direction, bore_radius,
                                       z_tip, bore_z_end);

    // Ball at the tip.  Rolled out by hand rather than via intersect_sphere()
    // so an interval lying entirely behind the origin still contributes its
    // endpoints to the hull.
    const double oz = origin.z() - z_tip;
    const double b  = origin.x() * direction.x() + origin.y() * direction.y()
                    + oz * direction.z();
    const double c  = origin.x() * origin.x() + origin.y() * origin.y()
                    + oz * oz - bore_radius * bore_radius;
    const double disc = b * b - c;

    std::optional<RayHit> ball;
    if (disc >= 0.0) {
        const double sqrt_disc = std::sqrt(disc);
        ball = RayHit{-b - sqrt_disc, -b + sqrt_disc};
    }

    if (!straight) return ball;
    if (!ball)     return straight;

    return RayHit{std::min(straight->t_enter, ball->t_enter),
                  std::max(straight->t_exit,  ball->t_exit)};
}

} // anonymous namespace


int intersect_shaped_bored_cylinder(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    double outer_radius,
    double bore_radius,
    double z_min,
    double z_max,
    double bore_z_start,
    double bore_z_end,
    const FrontFillet& fillet,
    bool rounded_bore_tip,
    RayHit segments_out[2])
{
    auto outer_hit = intersect_cylinder(origin, direction, outer_radius, z_min, z_max);
    if (!outer_hit || !outer_hit->valid()) {
        return 0; // Ray misses the outer crystal entirely
    }

    if (fillet.r_b > 0.0) {
        if (!correct_hit_for_front_fillet(origin, direction, fillet, *outer_hit)
            || !outer_hit->valid()) {
            return 0;
        }
    }

    auto bore_hit = rounded_bore_tip
        ? rounded_bore_interval(origin, direction, bore_radius,
                                bore_z_start, bore_z_end)
        : intersect_cylinder(origin, direction, bore_radius,
                             bore_z_start, bore_z_end);

    return subtract_bore_interval(*outer_hit, bore_hit, segments_out);
}

} // namespace ceelo
