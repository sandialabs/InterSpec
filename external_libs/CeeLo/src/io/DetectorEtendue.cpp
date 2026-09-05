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

#include "io/DetectorEtendue.h"
#include "io/LowDiscrepancy.h"
#include "materials/Material.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace ceelo {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// One face of the crystal's convex hull.
struct HullFace {
    enum class Kind { FrontDisc, FrontRect, CylSide, BoxSide };
    Kind kind;
    double area;            ///< cm^2
    Eigen::Vector3d normal; ///< outward unit normal (CylSide: per point, unused here)
    int axis = 0;           ///< BoxSide: 0 = x, 1 = y
    double sign = 1.0;      ///< BoxSide: +1 / -1
    double prob = 0.0;      ///< sampling probability
};

/// Point on a face from two unit-square coordinates, plus its outward normal.
void face_point(const HullFace& f, const Geometry& geom, double u, double v,
                Eigen::Vector3d& p, Eigen::Vector3d& n) {
    const double L = geom.detector_length();
    switch (f.kind) {
        case HullFace::Kind::FrontDisc: {
            const double R = geom.detector_radius();
            const double r = R * std::sqrt(u);
            const double ph = 2.0 * kPi * v;
            p = Eigen::Vector3d(r * std::cos(ph), r * std::sin(ph), 0.0);
            n = Eigen::Vector3d(0.0, 0.0, -1.0);
            return;
        }
        case HullFace::Kind::FrontRect: {
            p = Eigen::Vector3d((2.0 * u - 1.0) * geom.detector_half_x(),
                                (2.0 * v - 1.0) * geom.detector_half_y(), 0.0);
            n = Eigen::Vector3d(0.0, 0.0, -1.0);
            return;
        }
        case HullFace::Kind::CylSide: {
            const double R = geom.detector_radius();
            const double ph = 2.0 * kPi * u;
            const double c = std::cos(ph), s = std::sin(ph);
            p = Eigen::Vector3d(R * c, R * s, v * L);
            n = Eigen::Vector3d(c, s, 0.0);
            return;
        }
        case HullFace::Kind::BoxSide: {
            const double hx = geom.detector_half_x(), hy = geom.detector_half_y();
            if (f.axis == 0) {
                p = Eigen::Vector3d(f.sign * hx, (2.0 * u - 1.0) * hy, v * L);
                n = Eigen::Vector3d(f.sign, 0.0, 0.0);
            } else {
                p = Eigen::Vector3d((2.0 * u - 1.0) * hx, f.sign * hy, v * L);
                n = Eigen::Vector3d(0.0, f.sign, 0.0);
            }
            return;
        }
    }
}

/// The hull faces with their sampling probabilities.  Allocation is by
/// projected area toward the target direction `t` (from the hull centre), with
/// a floor of 5% of each face's true area so no face is starved: a side wall
/// still sees a far on-axis source at grazing incidence, and the per-line
/// |w.n| weight keeps the estimate unbiased whatever the allocation.
std::vector<HullFace> hull_faces(const Geometry& geom, const Eigen::Vector3d& t,
                                 bool have_direction) {
    std::vector<HullFace> faces;
    const double L = geom.detector_length();
    if (geom.shape() == DetectorShape::Cylinder) {
        const double R = geom.detector_radius();
        faces.push_back({HullFace::Kind::FrontDisc, kPi * R * R, {0.0, 0.0, -1.0}});
        faces.push_back({HullFace::Kind::CylSide, 2.0 * kPi * R * L, {0.0, 0.0, 0.0}});
    } else {
        const double hx = geom.detector_half_x(), hy = geom.detector_half_y();
        faces.push_back({HullFace::Kind::FrontRect, 4.0 * hx * hy, {0.0, 0.0, -1.0}});
        for (const double sgn : {+1.0, -1.0}) {
            HullFace fx{HullFace::Kind::BoxSide, 2.0 * hy * L, {sgn, 0.0, 0.0}, 0, sgn};
            HullFace fy{HullFace::Kind::BoxSide, 2.0 * hx * L, {0.0, sgn, 0.0}, 1, sgn};
            faces.push_back(fx);
            faces.push_back(fy);
        }
    }

    double norm = 0.0;
    for (HullFace& f : faces) {
        double proj = f.area;
        if (have_direction) {
            if (f.kind == HullFace::Kind::CylSide) {
                // Projected area of a cylinder's side wall: 2 R L sin(theta).
                const double st = std::sqrt(std::max(0.0, 1.0 - t.z() * t.z()));
                proj = 2.0 * geom.detector_radius() * L * st;
            } else {
                proj = f.area * std::max(0.0, t.dot(f.normal));
            }
        }
        f.prob = std::max(proj, 0.05 * f.area);
        norm += f.prob;
    }
    for (HullFace& f : faces) f.prob /= norm;
    return faces;
}

/// Orthonormal (e1, e2) perpendicular to a unit `a`.
void tangent_frame(const Eigen::Vector3d& a, Eigen::Vector3d& e1, Eigen::Vector3d& e2) {
    e1 = (std::abs(a.z()) < 0.9) ? a.cross(Eigen::Vector3d(0.0, 0.0, 1.0)).normalized()
                                 : a.cross(Eigen::Vector3d(1.0, 0.0, 0.0)).normalized();
    e2 = a.cross(e1);
}

}  // namespace

void sample_hull_points(const Geometry& geom, int n, const Eigen::Vector3d& toward,
                        bool have_direction, uint64_t index_offset,
                        std::vector<HullPoint>& out) {
    if (n <= 0) throw std::invalid_argument("sample_hull_points: n must be > 0");
    if (geom.detector_length() <= 0.0)
        throw std::invalid_argument("sample_hull_points: geometry has no crystal");

    const std::vector<HullFace> faces = hull_faces(geom, toward, have_direction);
    out.clear();
    out.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const uint64_t idx = index_offset + static_cast<uint64_t>(i);
        const double u_face = halton(idx, 2);
        const double u1 = halton(idx, 3);
        const double u2 = halton(idx, 5);

        size_t fi = 0;
        double acc = 0.0;
        for (; fi + 1 < faces.size(); ++fi) {
            acc += faces[fi].prob;
            if (u_face < acc) break;
        }
        const HullFace& face = faces[fi];

        HullPoint hp;
        face_point(face, geom, u1, u2, hp.point, hp.normal);
        hp.area_weight = face.area / face.prob;
        out.push_back(hp);
    }
}

bool append_etendue_line(EtendueLineSet& set, const Geometry& geom,
                         const Eigen::Vector3d& x, const Eigen::Vector3d& w_out,
                         double etendue_cm2_sr) {
    assert(w_out.z() < 0.0);
    assert(etendue_cm2_sr > 0.0);

    // The trace must start outside the outermost shell, on the source side.
    const std::pair<double, double> zext = geom.outer_z_extent();
    const double back = 2.0 * (geom.outer_bounding_radius()
                               + std::max(std::abs(zext.first), std::abs(zext.second)))
                        + x.norm() + 1.0;

    const Eigen::Vector3d dir_in = -w_out;
    const Eigen::Vector3d far_origin = x + w_out * back;
    std::vector<PathSegment> segs = geom.trace_ray(far_origin, dir_in);
    if (segs.empty()) return false;
    std::sort(segs.begin(), segs.end(),
              [](const PathSegment& a, const PathSegment& b) { return a.t_start < b.t_start; });

    KernelRay kr;
    kr.omega_w = static_cast<float>(etendue_cm2_sr / (4.0 * kPi));
    kr.cos_incidence = static_cast<float>(std::abs(dir_in.z()));
    kr.dir = dir_in.cast<float>();
    double active = 0.0;
    for (const PathSegment& s : segs) {
        const double len = s.length();
        if (len <= 1e-12) continue;
        if (s.is_scoring) active += len;
        if (s.material != nullptr)
            kr.segs.push_back({s.material, static_cast<float>(len), s.is_scoring});
    }
    kr.active_len = static_cast<float>(active);
    if (!(active > 0.0) || kr.segs.empty()) return false;   // can never contribute

    // Accumulate with the (float) weight the ray carries, as make_aperture_quadrature does.
    const double wf = kr.omega_w;
    const double prev = set.q.omega_frac_active;
    set.q.omega_frac_active += wf;
    // Running weighted means of the chord and incidence cosine.
    if (set.q.omega_frac_active > 0.0) {
        set.q.mean_chord_cm = (set.q.mean_chord_cm * prev + wf * active) / set.q.omega_frac_active;
        set.q.mean_cos_incidence =
            (set.q.mean_cos_incidence * prev + wf * kr.cos_incidence) / set.q.omega_frac_active;
    }

    set.q.rays.push_back(std::move(kr));
    set.origin.push_back(x);
    set.dir.push_back(dir_in);
    return true;
}

EtendueLineSet build_etendue_lines(const Geometry& geom,
                                   const Eigen::Vector3d& target_centre_cm,
                                   double target_radius_cm,
                                   int n_lines,
                                   uint64_t index_offset) {
    if (n_lines <= 0) throw std::invalid_argument("build_etendue_lines: n_lines must be > 0");

    EtendueLineSet set;
    set.target_centre = target_centre_cm;
    set.target_radius = target_radius_cm;
    set.q.cone_omega_frac = 0.0;
    set.q.n_rays_total = n_lines;
    set.q.rays.reserve(static_cast<size_t>(n_lines));
    set.origin.reserve(static_cast<size_t>(n_lines));
    set.dir.reserve(static_cast<size_t>(n_lines));

    // Face allocation from the hull centre toward the target (plain areas when
    // the target sphere contains the centre, or no target was given).
    const Eigen::Vector3d hull_centre(0.0, 0.0, 0.5 * geom.detector_length());
    const Eigen::Vector3d to_target = target_centre_cm - hull_centre;
    const bool have_dir = (target_radius_cm > 0.0) && (to_target.norm() > target_radius_cm);
    const Eigen::Vector3d t = have_dir ? to_target.normalized() : Eigen::Vector3d(0.0, 0.0, -1.0);

    std::vector<HullPoint> hull;
    sample_hull_points(geom, n_lines, t, have_dir, index_offset, hull);

    const double inv_n = 1.0 / static_cast<double>(n_lines);
    for (int i = 0; i < n_lines; ++i) {
        const uint64_t idx = index_offset + static_cast<uint64_t>(i);
        const double u3 = halton(idx, 7);
        const double u4 = halton(idx, 11);
        const Eigen::Vector3d& x = hull[static_cast<size_t>(i)].point;
        const Eigen::Vector3d& n = hull[static_cast<size_t>(i)].normal;

        // Direction cone toward the target sphere as seen from x.  `w` is the
        // OUTWARD direction (from the hull point toward the source side).
        double omega = 4.0 * kPi;
        Eigen::Vector3d w;
        const Eigen::Vector3d dvec = target_centre_cm - x;
        const double D = dvec.norm();
        if ((target_radius_cm > 0.0) && (D > target_radius_cm)) {
            const double ratio = target_radius_cm / D;
            const double cos_a = std::sqrt(std::max(0.0, 1.0 - ratio * ratio));
            omega = 2.0 * kPi * (1.0 - cos_a);
            const Eigen::Vector3d axis = dvec / D;
            Eigen::Vector3d e1, e2;
            tangent_frame(axis, e1, e2);
            const double ct = 1.0 - u3 * (1.0 - cos_a);
            const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
            const double ph = 2.0 * kPi * u4;
            w = (ct * axis + st * (std::cos(ph) * e1 + std::sin(ph) * e2)).normalized();
        } else {
            const double ct = 1.0 - 2.0 * u3;
            const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
            const double ph = 2.0 * kPi * u4;
            w = Eigen::Vector3d(st * std::cos(ph), st * std::sin(ph), ct);
        }

        // Leaving the hull, and toward the source side of the face plane;
        // otherwise this line is counted elsewhere (or nowhere) - weight zero,
        // and NOT renormalised away.
        const double cos_n = w.dot(n);
        if (!(cos_n > 0.0) || !(w.z() < 0.0)) continue;

        const double etendue = hull[static_cast<size_t>(i)].area_weight * omega * cos_n * inv_n;
        set.total_etendue += etendue;
        append_etendue_line(set, geom, x, w, etendue);
    }
    return set;
}

void line_interaction_probabilities(const ApertureQuadrature& q,
                                    double energy_keV, MuChoice mu,
                                    std::vector<double>& p_out,
                                    double passive_compton_recapture) {
    const double E_MeV = energy_keV * 1e-3;
    const double recap = (mu == MuChoice::NoRayleigh)
        ? std::max(0.0, std::min(1.0, passive_compton_recapture)) : 0.0;
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

    p_out.assign(q.rays.size(), 0.0);
    for (size_t i = 0; i < q.rays.size(); ++i) {
        const KernelRay& r = q.rays[i];
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;
        double p_int = 0.0;
        for (const RaySegment& s : r.segs) {
            const Mu m = mu_of(s.material);
            if (s.is_scoring) {
                const double mu_star = (mu == MuChoice::Total) ? m.tot : m.nors;
                p_int += std::exp(-tau_before) * (1.0 - std::exp(-mu_star * s.length));
                tau_before += m.tot * s.length;
            } else {
                tau_before += (m.tot - recap * m.cs) * s.length;
            }
        }
        p_out[i] = p_int;
    }
}

}  // namespace ceelo
