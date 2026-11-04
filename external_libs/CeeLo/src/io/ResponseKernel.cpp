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
#include "geometry/BoundingCone.h"

#include "cross_sections/CrossSectionData.h"
#include "transport/ComptonScatter.h"   // sample_rayleigh_cos_theta

#include <Eigen/Geometry>  // .cross()

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <stdexcept>
#include <utility>
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

/// A cone from `src` containing the active crystal envelope (envelope_bounding_cone, in
/// geometry/BoundingCone.cpp, does the search and carries its safety notes). Returns false when
/// none is returned safely, in which case the caller must sample the full sphere.
bool active_bounding_cone(const Geometry& geom, const Eigen::Vector3d& src,
                          Eigen::Vector3d& axis_out, double& cos_alpha_out) {
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

    return envelope_bounding_cone(cylinder, R, hx, hy, z0, z1, src, axis_out, cos_alpha_out);
}
}

void ApertureQuadrature::add_segment(KernelRay& ray, const Material* material, float length,
                                     bool is_scoring) {
    if (ray.seg_count == 0)
        ray.seg_begin = static_cast<uint32_t>(segments.size());
    assert(ray.seg_begin + ray.seg_count == segments.size());

    size_t slot = 0;
    while ((slot < materials.size()) && (materials[slot] != material))
        ++slot;
    if (slot == materials.size()) {
        if (slot > std::numeric_limits<uint8_t>::max())
            throw std::runtime_error("ApertureQuadrature: too many distinct materials");
        materials.push_back(material);
    }

    segments.push_back({material, length, is_scoring, static_cast<uint8_t>(slot)});
    ++ray.seg_count;
}

void ApertureQuadrature::discard_segments(KernelRay& ray) {
    if (ray.seg_count == 0)
        return;  // nothing added - and seg_begin may never have been set
    assert(ray.seg_begin + ray.seg_count == segments.size());
    segments.resize(segments.size() - ray.seg_count);
    ray.seg_count = 0;
}

void ApertureQuadrature::append_rays(const ApertureQuadrature& part) {
    rays.reserve(rays.size() + part.rays.size());
    segments.reserve(segments.size() + part.segments.size());
    for (KernelRay ray : part.rays) {
        const RaySegmentRange from = part.segs(ray);
        ray.seg_begin = static_cast<uint32_t>(segments.size());
        ray.seg_count = 0;
        for (const RaySegment& s : from)
            add_segment(ray, s.material, s.length, s.is_scoring);
        rays.push_back(ray);
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
    q.rays.reserve(n_rays);
    q.segments.reserve(4 * static_cast<size_t>(n_rays));

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

        // Geometry::trace_ray returns the segments in traversal order (both shape tracers sort).
        assert(std::is_sorted(segs.begin(), segs.end(),
                              [](const PathSegment& a, const PathSegment& b) {
                                  return a.t_start < b.t_start;
                              }));

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
                q.add_segment(kr, s.material, static_cast<float>(len), s.is_scoring);
        }
        kr.active_len = static_cast<float>(active);
        if (kr.seg_count == 0) continue;

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
    // (mu_total, mu_nors, mu_cs) of each material, by RaySegment::material_slot.
    struct Mu { double tot, nors, cs; };
    std::vector<Mu> mus;
    mus.reserve(materials.size());
    for (const Material* m : materials) {
        const MacroscopicXS xs = m->macroscopic_xs(E_MeV);
        mus.push_back({xs.mu_total(), xs.mu_pe + xs.mu_cs + xs.mu_pp, xs.mu_cs});
    }
    for (const auto& r : rays) {
        if (r.active_len <= 0.0f) continue;
        double tau_before = 0.0;  // optical depth of everything traversed so far
        double p_int = 0.0;
        for (const auto& s : segs(r)) {
            const Mu& m = mus[s.material_slot];
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
    std::vector<double> mus;  // mu_total by RaySegment::material_slot
    mus.reserve(materials.size());
    for (const Material* m : materials)
        mus.push_back(m->mu_total(E_MeV));
    for (const auto& r : rays) {
        if (r.active_len <= 0.0f) continue;
        double tau = 0.0;
        for (const auto& s : segs(r)) {
            if (s.is_scoring) break;  // transmission up to first active entry
            tau += mus[s.material_slot] * s.length;
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

// ---- Rayleigh deflection loss ---------------------------------------------------------------

namespace {

/// Histogram of a = sec(theta) - 1 over the Rayleigh angular law of one (Z, E): forward samples
/// binned in log(a), the backward (theta >= 90 deg) fraction kept separately.  Weights sum to 1.
struct RayleighDeflectionHist {
    std::vector<double> a;     ///< bin centre, sec(theta) - 1
    std::vector<double> w;     ///< bin weight
    double backward = 0.0;     ///< fraction scattered through >= 90 deg
};

const RayleighDeflectionHist& rayleigh_deflection_hist(int Z, double E_keV) {
    static std::mutex mutex;
    static std::map<std::pair<int, long long>, RayleighDeflectionHist> cache;

    const std::pair<int, long long> key(Z, std::llround(E_keV * 1000.0));
    std::lock_guard<std::mutex> lock(mutex);
    const auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    constexpr int kSamples = 50000;
    constexpr int kBins = 128;
    constexpr double kLogAMin = -4.0, kLogAMax = 3.0;   // a in [1e-4, 1e3]; below => ~no loss
    constexpr double kLogSpan = kLogAMax - kLogAMin;

    RayleighDeflectionHist h;
    h.a.resize(kBins);
    h.w.assign(kBins, 0.0);
    for (int b = 0; b < kBins; ++b)
        h.a[b] = std::pow(10.0, kLogAMin + (b + 0.5) * kLogSpan / kBins);

    // Fixed seed: the histogram (and so every downstream coefficient) is reproducible.
    std::mt19937_64 rng(0x5EED2026ULL ^ (static_cast<uint64_t>(Z) << 32)
                        ^ static_cast<uint64_t>(key.second));
    const double dw = 1.0 / kSamples;
    for (int i = 0; i < kSamples; ++i) {
        const double ct = sample_rayleigh_cos_theta(Z, E_keV, rng);
        if (ct <= 0.0) { h.backward += dw; continue; }
        const double a = std::max(1.0 / ct - 1.0, 1e-300);
        int b = static_cast<int>(std::floor((std::log10(a) - kLogAMin) / kLogSpan * kBins));
        b = std::max(0, std::min(kBins - 1, b));
        h.w[b] += dw;
    }
    return cache.emplace(key, std::move(h)).first->second;
}

}  // namespace

double rayleigh_deflection_loss_fraction(double E_keV, double tau_nr, const Material& mat) {
    if (!(tau_nr > 0.0) || !(E_keV > 0.0)) return 0.0;

    const CrossSectionData& xsd = CrossSectionData::instance();
    const double E_MeV = E_keV * 1e-3;

    // Per-element weights n_i * sigma_R,i: each element's share of the Rayleigh interactions.
    double num = 0.0, den = 0.0;
    for (const MaterialComponent& c : mat.composition()) {
        const double w = c.mass_fraction / xsd.atomic_weight(c.Z) * xsd.sigma_rayleigh(c.Z, E_MeV);
        if (!(w > 0.0)) continue;
        const RayleighDeflectionHist& h = rayleigh_deflection_hist(c.Z, E_keV);
        // For a scatter uniformly placed along the layer, the surviving fraction of the
        //  lengthened remaining path averages to (1 - e^{-aT}) / (aT).
        double lost = h.backward;
        for (size_t b = 0; b < h.a.size(); ++b) {
            const double x = h.a[b] * tau_nr;
            const double surv = (x < 1e-8) ? 1.0 : (-std::expm1(-x) / x);
            lost += h.w[b] * (1.0 - surv);
        }
        num += w * lost;
        den += w;
    }
    return (den > 0.0) ? num / den : 0.0;
}

} // namespace ceelo
