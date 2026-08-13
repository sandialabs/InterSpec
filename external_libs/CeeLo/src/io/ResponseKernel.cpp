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

namespace ceelo {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kHcKeVAngstrom = 12.398;  // hc, keV*Angstrom (x = E/hc sin(t/2))
}

ApertureQuadrature make_aperture_quadrature(const Geometry& geom,
                                            const Eigen::Vector3d& src,
                                            int n_rays) {
    ApertureQuadrature q;

    // Bounding sphere of the outermost shell.
    const double r_out = geom.outer_bounding_radius();
    const auto [z_min, z_max] = geom.outer_z_extent();
    const Eigen::Vector3d center(0.0, 0.0, 0.5 * (z_min + z_max));
    const double half_z = 0.5 * (z_max - z_min);
    const double r_sphere = std::sqrt(r_out * r_out + half_z * half_z);

    const Eigen::Vector3d to_c = center - src;
    const double dist = to_c.norm();

    double cos_alpha = -1.0;  // full sphere
    Eigen::Vector3d axis(0, 0, 1);
    if (dist > r_sphere * 1.0001) {
        const double sin_a = r_sphere / dist;
        cos_alpha = std::sqrt(std::max(0.0, 1.0 - sin_a * sin_a));
        axis = to_c / dist;
    }
    q.cone_omega_frac = 0.5 * (1.0 - cos_alpha);

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
            q.omega_frac_active += w_per_ray;
            sum_chord += w_per_ray * active;
            sum_cos += w_per_ray * kr.cos_incidence;
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
