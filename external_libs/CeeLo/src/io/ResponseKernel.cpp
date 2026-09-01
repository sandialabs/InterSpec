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

/// @file ResponseKernel.cpp
/// @brief MC-free aperture-quadrature geometry kernel (see ResponseKernel.h).

#include "io/ResponseKernel.h"

#include "cross_sections/CrossSectionData.h"

#include <Eigen/Geometry>  // .cross()

#include <algorithm>
#include <vector>

namespace ceelo {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kHcKeVAngstrom = 12.398;  // hc, keV*Angstrom (x = E/hc sin(t/2))

// --- Smallest cone containing the active crystal ------------------------------------------
//
// SAFETY MODEL: for any candidate axis the routines below compute the EXACT extreme angle over
// the body, so the resulting cone always contains it. The axis search only affects how tight the
// cone is, never whether it is valid -- a poor search costs efficiency, it cannot lose paths.
// (This is not hypothetical: a mean-direction-plus-bisector heuristic tried first produced cones
// wider than the bounding sphere's, which is impossible for a body inside that sphere.)

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

/// A cone from `src` containing the active crystal envelope. Returns false when none is returned
/// safely, in which case the caller must sample the full sphere.
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
bool active_bounding_cone(const Geometry& geom, const Eigen::Vector3d& src,
                          Eigen::Vector3d& axis_out, double& cos_alpha_out) {
    // Refuse anything not comfortably inside the evaluators' valid domain. A cone at 90 degrees is
    // worth at most a 2x yield gain over the full sphere, so declining near there costs ~nothing.
    constexpr double kMinCos = 1e-6;

    // The scoring envelope. Normally this is just the crystal: a dead layer eats INWARD and the
    // bore hole and bulletization only remove material, so radius x [0, L] bounds the active
    // volume. But `set_dead_layer` accepts NEGATIVE thicknesses (nothing rejects them, and they
    // survive an XML round-trip), and the tracer then scores
    //     active_r = R - side,  active_z = [front, L - back]
    // (RayTrace.cpp, trace_cylinder_geometry) -- which for a negative component lies OUTSIDE the
    // crystal. The whole cone rests on "active is inside the envelope", so derive the envelope
    // from those same expressions rather than assuming the signs. Measured before this was
    // handled: a -0.5 cm front dead layer silently deleted 1.66% of the active solid angle at an
    // ordinary far-field query, with no error anywhere.
    double z0 = 0.0;
    double z1 = geom.detector_length();
    double grow = 0.0;                          // outward growth of the transverse extent
    if (geom.has_dead_layer()) {
        const DeadLayerConfig& dl = *geom.dead_layer();
        z0 = std::min(z0, dl.front);
        z1 = std::max(z1, geom.detector_length() - dl.back);
        grow = std::max(0.0, -dl.side);
    }
    if (!(z1 > z0)) return false;

    const bool cylinder = (geom.shape() == DetectorShape::Cylinder);
    const double R  = cylinder ? (geom.detector_radius() + grow) : 0.0;
    const double hx = cylinder ? 0.0 : (geom.detector_half_x() + grow);
    const double hy = cylinder ? 0.0 : (geom.detector_half_y() + grow);

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
}

ApertureQuadrature make_aperture_quadrature(const Geometry& geom,
                                            const Eigen::Vector3d& src,
                                            int n_rays) {
    ApertureQuadrature q;

    // Smallest cone containing the ACTIVE CRYSTAL (see header). Only rays with an active chord
    // are ever consumed, so the endcap and attenuators do not need to be inside the cone.
    double cos_alpha = -1.0;  // full sphere
    Eigen::Vector3d axis(0, 0, 1);
    if (!active_bounding_cone(geom, src, axis, cos_alpha)) {
        cos_alpha = -1.0;
        axis = Eigen::Vector3d(0, 0, 1);
    }
    q.cone_omega_frac = 0.5 * (1.0 - cos_alpha);
    q.cone_axis = axis;

    // Orthonormal frame around the cone axis.
    const Eigen::Vector3d u = std::abs(axis.z()) < 0.9
        ? axis.cross(Eigen::Vector3d(0, 0, 1)).normalized()
        : axis.cross(Eigen::Vector3d(1, 0, 0)).normalized();
    const Eigen::Vector3d v = axis.cross(u);

    const double golden = kPi * (3.0 - std::sqrt(5.0));  // Fibonacci spiral step
    const double w_per_ray = q.cone_omega_frac / n_rays; // equal-solid-angle rays
    q.n_rays_total = n_rays;

    double sum_chord = 0.0, sum_cos = 0.0;
    std::vector<PathSegment> segs;
    for (int i = 0; i < n_rays; ++i) {
        const double f = (i + 0.5) / n_rays;
        const double ct = 1.0 - f * (1.0 - cos_alpha);
        const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
        const double ph = golden * i;
        const Eigen::Vector3d dir =
            (ct * axis + st * (std::cos(ph) * u + std::sin(ph) * v)).normalized();

        segs.clear();
        geom.trace_ray(src, dir, segs);
        if (segs.empty()) continue;

        // trace_ray does not guarantee global ordering across shells; sort.
        std::sort(segs.begin(), segs.end(),
                  [](const PathSegment& a, const PathSegment& b) {
                      return a.t_start < b.t_start;
                  });

        KernelRay kr;
        kr.omega_w = static_cast<float>(w_per_ray);
        kr.cos_incidence = static_cast<float>(std::abs(dir.z()));
        kr.dir = dir.cast<float>();
        double active = 0.0;
        for (const auto& s : segs) {
            const double len = s.length();
            if (len <= 1e-12) continue;
            if (s.is_scoring) active += len;
            if (s.material != nullptr)
                kr.segs.push_back({s.material, static_cast<float>(len),
                                   s.is_scoring});
        }
        kr.active_len = static_cast<float>(active);
        if (kr.segs.empty()) continue;

        if (active > 0.0) {
            // Accumulate with the SAME (float) weight the ray carries, not the double w_per_ray.
            // transmitted_omega() and interaction_omega() sum kr.omega_w, so mixing the two
            // precisions here left omega_frac_active disagreeing with them by ~float epsilon --
            // enough to break the bare-detector identity transmitted_omega == omega_frac_active,
            // which is exact physics, at a tolerance that had been passing only by luck.
            const double w = kr.omega_w;
            q.omega_frac_active += w;
            sum_chord += w * active;
            sum_cos += w * kr.cos_incidence;
        }
        q.rays.push_back(std::move(kr));
    }
    if (q.omega_frac_active > 0.0) {
        q.mean_chord_cm = sum_chord / q.omega_frac_active;
        q.mean_cos_incidence = sum_cos / q.omega_frac_active;
    }
    return q;
}

double ApertureQuadrature::interaction_omega(
    double energy_keV, MuChoice mu, double passive_compton_recapture) const {
    const double E_MeV = energy_keV * 1e-3;
    // Scatter-in credit only for the TOTAL kernel (see header): a degraded
    // Compton photon can still deposit but never lands in the FEP peak.
    const double recap = (mu == MuChoice::NoRayleigh)
        ? std::max(0.0, std::min(1.0, passive_compton_recapture)) : 0.0;
    double total = 0.0;
    // Per-call cache of (mu_total, mu_nors, mu_cs) by material pointer.
    struct Mu { double tot, nors, cs; };
    std::vector<std::pair<const Material*, Mu>> mu_cache;
    auto mu_of = [&](const Material* m) -> Mu {
        for (const auto& e : mu_cache)
            if (e.first == m) return e.second;
        const MacroscopicXS xs = m->macroscopic_xs(E_MeV);
        const Mu val{xs.mu_total(), xs.mu_pe + xs.mu_cs + xs.mu_pp, xs.mu_cs};
        mu_cache.push_back({m, val});
        return val;
    };
    for (const auto& r : rays) {
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;  // optical depth of everything traversed so far
        double p_int = 0.0;
        for (const auto& s : r.segs) {
            const Mu m = mu_of(s.material);
            if (s.is_scoring) {
                const double mu_star = (mu == MuChoice::Total) ? m.tot : m.nors;
                p_int += std::exp(-tau_before) * (1.0 - std::exp(-mu_star * s.length));
                tau_before += m.tot * s.length;
            } else {
                // Passive segment: credit forward-Compton scatter-in (total only).
                tau_before += (m.tot - recap * m.cs) * s.length;
            }
        }
        total += r.omega_w * p_int;
    }
    return total;
}

double ApertureQuadrature::transmitted_omega(double energy_keV) const {
    const double E_MeV = energy_keV * 1e-3;
    double total = 0.0;
    std::vector<std::pair<const Material*, double>> mu_cache;
    auto mu_of = [&](const Material* m) -> double {
        for (const auto& e : mu_cache)
            if (e.first == m) return e.second;
        const double v = m->mu_total(E_MeV);
        mu_cache.push_back({m, v});
        return v;
    };
    for (const auto& r : rays) {
        if (r.active_len <= 0.0f) continue;
        double tau = 0.0;
        for (const auto& s : r.segs) {
            if (s.is_scoring) break;  // transmission up to first active entry
            tau += mu_of(s.material) * s.length;
        }
        total += r.omega_w * std::exp(-tau);
    }
    return total;
}

double kn_in_window_fraction(double E_keV, double win_keV) {
    const double k = E_keV / 511.0;
    auto dsdmu = [&](double mu) {
        const double r = 1.0 / (1.0 + k * (1.0 - mu));  // E'/E
        return r * r * (r + 1.0 / r - (1.0 - mu * mu));
    };
    // Window edge: E' = E/(1+k(1-mu)) >= E - win  =>  1-mu <= win/(k(E-win)).
    const double one_minus_mu_max =
        (E_keV > win_keV) ? win_keV / (k * (E_keV - win_keV)) : 2.0;
    const double mu_lo = std::max(-1.0, 1.0 - one_minus_mu_max);
    auto integrate = [&](double a, double b) {
        const int n = 200;
        double s = 0.0;
        for (int i = 0; i < n; ++i) {
            const double x0 = a + (b - a) * i / n, x1 = a + (b - a) * (i + 1) / n;
            s += (dsdmu(x0) + 4.0 * dsdmu(0.5 * (x0 + x1)) + dsdmu(x1)) *
                 (x1 - x0) / 6.0;
        }
        return s;
    };
    return integrate(mu_lo, 1.0) / integrate(-1.0, 1.0);
}

double kn_in_window_fraction(double E_keV, double win_keV, const Material& mat) {
    const CrossSectionData& xsd = CrossSectionData::instance();
    // Per-element ATOM weights n_i ~ w_mass_i / A_i: the macroscopic incoherent
    // angular distribution is sum_i n_i (dsigma_KN/dmu) S(x, Z_i). (Identical
    // to electron-fraction weights w_i Z_i/A_i times the normalized S/Z_i; the
    // overall scale cancels in the ratio.)
    struct ElemWeight { int Z; double n; };
    std::vector<ElemWeight> elems;
    for (const MaterialComponent& c : mat.composition())
        elems.push_back({c.Z, c.mass_fraction / xsd.atomic_weight(c.Z)});
    const double k = E_keV / 511.0;
    auto dsdmu = [&](double mu) {
        const double r = 1.0 / (1.0 + k * (1.0 - mu));  // E'/E
        const double kn = r * r * (r + 1.0 / r - (1.0 - mu * mu));
        // x = sin(theta/2)/lambda = (E/hc) sin(theta/2), inverse Angstrom --
        // same convention as the transport sampler (ComptonScatter.cpp).
        const double sin_half = std::sqrt(std::max(0.0, (1.0 - mu) * 0.5));
        const double x = (E_keV / kHcKeVAngstrom) * sin_half;
        double sw = 0.0;
        for (const ElemWeight& e : elems)
            sw += e.n * xsd.scattering_function_S(e.Z, x);
        return kn * sw;
    };
    // Window edge: E' = E/(1+k(1-mu)) >= E - win  =>  1-mu <= win/(k(E-win)).
    const double one_minus_mu_max =
        (E_keV > win_keV) ? win_keV / (k * (E_keV - win_keV)) : 2.0;
    const double mu_lo = std::max(-1.0, 1.0 - one_minus_mu_max);
    auto integrate = [&](double a, double b) {
        const int n = 200;
        double s = 0.0;
        for (int i = 0; i < n; ++i) {
            const double x0 = a + (b - a) * i / n, x1 = a + (b - a) * (i + 1) / n;
            s += (dsdmu(x0) + 4.0 * dsdmu(0.5 * (x0 + x1)) + dsdmu(x1)) *
                 (x1 - x0) / 6.0;
        }
        return s;
    };
    const double denom = integrate(-1.0, 1.0);
    return (denom > 0.0) ? integrate(mu_lo, 1.0) / denom : 0.0;
}

} // namespace ceelo
