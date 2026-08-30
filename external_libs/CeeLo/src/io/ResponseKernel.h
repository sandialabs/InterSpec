#ifndef CEELO_IO_RESPONSE_KERNEL_H
#define CEELO_IO_RESPONSE_KERNEL_H
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

/// @file ResponseKernel.h
/// @brief MC-free geometry kernel for the parameterized detector response.
///
/// Deterministic aperture quadrature over the exact ray-traced detector
/// geometry: from a source point, equal-solid-angle rays (Fibonacci spiral)
/// are cast into the cone circumscribing the geometry bounding sphere and
/// traced through every shell (attenuators, dead layer, active crystal).
/// The resulting per-ray material path lists let any energy-dependent
/// quantity be evaluated analytically -- no Monte Carlo:
///
///   - omega_frac_active: fractional solid angle of the active crystal
///   - interaction_omega(E, mu*): the interaction-weighted effective solid
///     angle  K(E) = sum_i w_i sum_{scoring segs} exp(-tau_before(E)) *
///     (1 - exp(-mu*(E) l_seg))  -- the distance/angle kernel of the
///     parameterized response (response_function_spec.md Eq. 1)
///   - transmitted_omega(E): passive-layer transmission envelope (no chord
///     saturation); used for the collimator hole-fraction validity gate
///   - interaction_omega_with_tsrc: K with a per-ray source-transmission
///     factor folded in (extended/shielded sources, spec Eq. 5/8)
///
/// Validated in the July 2026 parameterization campaign (S1: far-field flat
/// at MC noise with no fitted parameter; 1%-validity down to d ~ 1.1-1.8 a).
/// Promoted from the detector-response parameterization campaign.
///
/// Conventions: crystal front face at z = 0, crystal along +z; source at
/// (d sin(t) cos(p), d sin(t) sin(p), -d cos(t)) for distance d measured to
/// the crystal-face origin, polar angle t from the detector axis, azimuth p.

#include "physics/FepWindow.h"
#include "geometry/Geometry.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace ceelo {

/// One material stretch along a quadrature ray, in traversal order.
struct RaySegment {
    const Material* material;  ///< never nullptr (vacuum stretches are skipped)
    float length;              ///< cm
    bool is_scoring;           ///< active crystal
};

/// One quadrature ray that intersects the geometry's outer boundary.
struct KernelRay {
    float omega_w;             ///< this ray's solid-angle weight, fraction of 4pi
    float active_len;          ///< total active-crystal chord length (cm)
    float cos_incidence;       ///< |cos| of angle vs detector axis (+z)
    Eigen::Vector3f dir;       ///< ray direction (for per-ray source attenuation)
    std::vector<RaySegment> segs;  ///< ordered material stretches (few entries)
};

/// Which mu* the scoring (active-crystal) segments use in interaction_omega.
enum class MuChoice : uint8_t {
    Total,       ///< mu_pe + mu_cs + mu_rs + mu_pp  (FEP kernel K)
    NoRayleigh   ///< mu_pe + mu_cs + mu_pp  (eps_tot kernel; Rayleigh is elastic)
};

/// Deterministic aperture quadrature from one source point (see file docs).
/// Built by make_aperture_quadrature(); all evaluations are MC-free.
struct ApertureQuadrature {
    double cone_omega_frac = 0.0;   ///< sampling-cone Omega/4pi (1.0 = full sphere)
    int n_rays_total = 0;           ///< rays sampled (incl. misses)
    std::vector<KernelRay> rays;    ///< rays with any geometry intersection
    double omega_frac_active = 0.0; ///< Omega/4pi of the ACTIVE crystal (chord > 0)
    double mean_chord_cm = 0.0;     ///< Omega-weighted mean active chord
    double mean_cos_incidence = 0.0;///< Omega-weighted (flux) mean |cos| vs axis;
                                    ///< the near-field model's theta_eff coordinate

    /// Interaction-weighted effective solid angle (spec Eq. 1); also the
    /// chord-model estimate of "any interaction in the active crystal":
    ///   sum_i w_i * sum_{scoring segs} T_before(E) * (1 - exp(-mu*(E) l_seg))
    /// T_before uses mu_total of every preceding segment (any interaction in a
    /// passive layer removes/degrades the photon); mu* on the scoring segment
    /// is picked by MuChoice.
    ///
    /// `passive_compton_recapture` r in [0,1] credits DEAD-LAYER / ENDCAP
    /// scatter-in for the TOTAL kernel only (mu == NoRayleigh): a fraction r of
    /// the Compton interactions in each NON-scoring segment forward-scatter into
    /// the active crystal and still deposit, so they are removed from the
    /// pre-scoring attenuation -- T_before uses (mu_total - r*mu_cs) on passive
    /// segments instead of mu_total. r acts on the per-ray passive PATH length,
    /// so the credit grows with oblique/near-field incidence (longer passive
    /// chords) and vanishes far-field -- exactly the near-field total systematic
    /// it corrects. r = 0 (default) is bit-identical to the removal-only kernel;
    /// the FEP kernel (mu == Total) ignores r (a degraded Compton photon cannot
    /// land in the peak).
    double interaction_omega(double energy_keV, MuChoice mu,
                             double passive_compton_recapture = 0.0) const;

    /// Geometric envelope solid angle with passive-layer transmission only (no
    /// chord saturation): sum_i w_i * T_passive_before_first_scoring(E).
    double transmitted_omega(double energy_keV) const;

    /// Collimator validity gate: fraction of the active geometric solid angle
    /// that survives the passive layers, s = transmitted_omega / omega_active.
    /// s < ~0.3 means the response is dominated by collimator scatter/leak the
    /// removal-only kernel cannot represent (warn / inflate sigma); s < ~0.05
    /// should refuse or fall back to bespoke MC (spec sec 4.5).
    double hole_fraction(double energy_keV) const {
        return (omega_frac_active > 0.0)
            ? transmitted_omega(energy_keV) / omega_frac_active : 0.0;
    }

    /// interaction_omega with an additional per-ray transmission factor
    /// T_src(dir) folded in BEFORE the detector layers -- e.g. the source-
    /// geometry survival probability for extended/shielded sources (the
    /// per-element extended-source kernel, spec Eq. 5; use
    /// fep_survival_removal_mu for the FEP transmission).
    template <typename TsrcFn>
    double interaction_omega_with_tsrc(double energy_keV, MuChoice mu,
                                       TsrcFn&& t_src,
                                       double passive_compton_recapture = 0.0) const;
};

/// Build the quadrature from a source point. n_rays is the ray count in the
/// sampling cone (Fibonacci spiral; deterministic). The cone subtends the
/// geometry bounding sphere; if the source is inside it, the full sphere is
/// used. Cost ~ n_rays ray traces (sub-ms at 1k-8k rays); cacheable per
/// position, and in fit loops the position set is fixed so kernels are
/// computed once.
ApertureQuadrature make_aperture_quadrature(const Geometry& geom,
                                            const Eigen::Vector3d& src,
                                            int n_rays = 8192);

/// Source position for (d, cos_theta, phi); d measured from the origin
/// (center of the crystal front face).
inline Eigen::Vector3d source_position(double d_cm, double cos_theta,
                                       double phi_rad = 0.0) {
    const double st = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    return Eigen::Vector3d(d_cm * st * std::cos(phi_rad),
                           d_cm * st * std::sin(phi_rad),
                           -d_cm * cos_theta);
}

/// Fraction of Klein-Nishina Compton scatters whose energy loss stays within
/// the +-win_keV FEP window (forward cone; theta <= 37 deg at 60 keV).
/// Free-electron KN only -- S(x,Z) suppresses forward scatter, so this
/// slightly overestimates the surviving fraction (measured +1..2% mid-E
/// overshoot, S9b). Simpson integration over cos(theta).
double kn_in_window_fraction(double E_keV, double win_keV = kDefaultFepWindowKeV);

/// Material-aware in-window fraction: weights the Klein-Nishina integrand by
/// the incoherent scattering function of the material's elements,
/// sum_i n_i S(x, Z_i) with atom weights n_i ~ w_mass_i / A_i and momentum
/// transfer x = (E/hc) sin(theta/2) in inverse Angstrom (the same convention
/// as CrossSectionData::scattering_function_S and the transport sampler).
/// Normalized by the identical S-weighted integral over the full mu range,
/// so it remains the proper fraction of INCOHERENT scatters staying
/// in-window -- consistent with multiplying the tabulated (bound) mu_CS in
/// fep_survival_removal_mu. S(x,Z) suppresses forward scatter, so this is
/// BELOW the free-electron fraction, most strongly at high Z (measured
/// ratios @60 keV: water 0.89, Fe 0.77, Pb 0.63). The RATIO does not
/// approach 1 at high E -- the window-edge momentum transfer
/// x_edge ~ sqrt(m_e c^2 win/2)/hc ~ 1.6 A^-1 is energy-independent -- but
/// the absolute correction vanishes (f_win -> 0), so mu_rem converges to
/// the free-electron value to <0.1% at >= 1332 keV. Use when the
/// source/shield material is known; hoist per (E, material) out of
/// per-element/per-ray loops.
double kn_in_window_fraction(double E_keV, double win_keV, const Material& mat);

/// Survival removal cross-section for per-element extended-source FEP
/// transmission (spec Eq. 8):
///     mu_rem = mu_total - mu_Rayleigh - f_win(E) * mu_Compton
/// Rayleigh is elastic (cannot remove a photon from the FEP window) and the
/// KN in-window fraction captures forward Compton scatters whose energy loss
/// stays in-window. Measured (S9b): every tested shielded/self-attenuating
/// case within +-2.4% FEP across 60-2614 keV (plain mu_total gave -8..-16%
/// at 60 keV). For eps_tot transmission use mu_total (no survival credit).
/// Pass f_win = kn_in_window_fraction(E) (hoist out of per-element loops).
inline double fep_survival_removal_mu(const MacroscopicXS& xs, double f_win) {
    return xs.mu_total() - xs.mu_rs - f_win * xs.mu_cs;
}

/// Depth-aware forward-Compton survival credit g(tau_src), in (0, 1].
/// fep_survival_removal_mu credits ALL forward-in-window Compton as surviving,
/// but the +-win_keV energy window admits scattering angles (37 deg at 60 keV)
/// far wider than a distant detector's acceptance; after tau_src source optical
/// depths a forward-Compton "survivor" has been deflected out of acceptance
/// (and multiply-scattered survivors lose the forward correlation), so the
/// credit must decay with the accumulated source-path optical depth.  Modeled
/// as g(tau) = exp(-tau/tau_c) with g(0) = 1 (so the depth-aware removal mu
/// reduces EXACTLY to fep_survival_removal_mu at the surface).  tau_c fitted to
/// the deep-shield steel-pipe contents truth (Stage E3 A2): tau_c = 5 mfp best
/// removes the forward-Compton over-credit on the r<=1-mfp/few-mfp pipes
/// (60 keV +9..+15% -> ~1%) while leaving the >=122 keV rows within a few %.
/// It does NOT (and physically cannot) close the deepest 9-mfp / 60 keV case
/// (+195% -> +128%): that residual is a scatter-dominated regime the removal
/// kernel cannot represent -- callers should inflate sigma for large tau_src
/// (see kBuildupSigmaPerMfp) rather than trust the point estimate there.
constexpr double kFepDepthTauC = 5.0;   ///< source-path mfp scale for g(tau)
inline double fep_depth_survival_credit(double tau_src) {
    return std::exp(-std::max(0.0, tau_src) / kFepDepthTauC);
}

/// Depth-aware fep_survival_removal_mu: the in-window Compton credit is scaled
/// by fep_depth_survival_credit(tau_src).  tau_src is the FEP photon's
/// accumulated source-path optical depth (mu_total-weighted path length, in
/// mfp).  Reduces EXACTLY to fep_survival_removal_mu(xs, f_win) at tau_src = 0
/// (thin-shell / self-attenuating skin sources, where the credit is valid),
/// so existing thin-source callers are unaffected.  Prefer this over the flat
/// fep_survival_removal_mu for extended/shielded sources where a ray can
/// accumulate several mfp before leaving the source.
inline double fep_survival_removal_mu_depth(const MacroscopicXS& xs,
                                            double f_win, double tau_src) {
    return xs.mu_total() - xs.mu_rs
           - f_win * fep_depth_survival_credit(tau_src) * xs.mu_cs;
}

// ---- template implementation ----------------------------------------------

template <typename TsrcFn>
double ApertureQuadrature::interaction_omega_with_tsrc(
    double energy_keV, MuChoice mu, TsrcFn&& t_src,
    double passive_compton_recapture) const {
    const double E_MeV = energy_keV * 1e-3;
    // Scatter-in credit only for the TOTAL kernel (a degraded Compton photon
    // can still deposit energy, but never lands in the FEP peak).
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
        double tau_before = 0.0;
        double p_int = 0.0;
        for (const auto& s : r.segs) {
            const Mu m = mu_of(s.material);
            if (s.is_scoring) {
                const double mu_star = (mu == MuChoice::Total) ? m.tot : m.nors;
                p_int += std::exp(-tau_before) * (1.0 - std::exp(-mu_star * s.length));
                tau_before += m.tot * s.length;
            } else {
                // Passive (non-scoring) segment: credit forward-Compton
                // scatter-in for the total kernel.
                tau_before += (m.tot - recap * m.cs) * s.length;
            }
        }
        if (p_int > 0.0)
            total += r.omega_w * p_int * t_src(r.dir.template cast<double>().eval());
    }
    return total;
}

} // namespace ceelo

#endif // CEELO_IO_RESPONSE_KERNEL_H
