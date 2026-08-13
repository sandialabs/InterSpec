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

#include "transport/ComptonScatter.h"
#include "cross_sections/CrossSectionData.h"

#include <algorithm>
#include <cmath>
#include <cassert>
#include <map>
#include <utility>

namespace ceelo {

namespace {

/// Electron rest mass energy in keV
constexpr double kElectronMassKeV = 510.998950;

/// Thomson cross section in barns
constexpr double kThomsonXS = 0.665245873e0; // 8π/3 * r_e² in barn

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

/// hc in keV·Å — for momentum transfer parameter x = (E/hc) * sin(θ/2)
constexpr double kHcKeVAngstrom = 12.398;

/// Fine-structure constant: 1 atomic unit of momentum = kAlphaFS * m_e * c.
/// Converts Compton-profile momenta (atomic units) to the dimensionless
/// t = p_z/(m_e c) used in the impulse-approximation kinematics.
constexpr double kAlphaFS = 7.2973525693e-3;

} // anonymous namespace

double compton_scattered_energy(double energy_keV, double cos_theta) {
    double alpha = energy_keV / kElectronMassKeV;
    return energy_keV / (1.0 + alpha * (1.0 - cos_theta));
}

// ---------------------------------------------------------------------------
// Compton Doppler broadening (impulse approximation, PENELOPE-style analytic
// one-parameter subshell profiles).  The angle is sampled from KN × S(x,Z)/Z
// as before; the profile sample only redistributes the scattered energy at
// fixed angle — PENELOPE's factorization of the IA double-differential cross
// section (Ribberfors), so there is no double counting with S(x,Z).
// ---------------------------------------------------------------------------

double compton_profile_cdf(double J0, double p_au) {
    // F(p) = ½ exp{½[1 − (1 − 2 J0 p)²]}        p ≤ 0
    //      = 1 − ½ exp{½[1 − (1 + 2 J0 p)²]}    p ≥ 0
    if (p_au <= 0.0) {
        double u = 1.0 - 2.0 * J0 * p_au;
        return 0.5 * std::exp(0.5 * (1.0 - u * u));
    }
    double u = 1.0 + 2.0 * J0 * p_au;
    return 1.0 - 0.5 * std::exp(0.5 * (1.0 - u * u));
}

double compton_profile_invcdf(double J0, double xi) {
    xi = std::max(1e-12, std::min(1.0 - 1e-12, xi));
    if (xi <= 0.5) {
        return (1.0 - std::sqrt(1.0 - 2.0 * std::log(2.0 * xi))) / (2.0 * J0);
    }
    return (std::sqrt(1.0 - 2.0 * std::log(2.0 * (1.0 - xi))) - 1.0) / (2.0 * J0);
}

double doppler_scattered_energy(double energy_keV, double cos_theta, double t) {
    // E′(t) = E [ (A − t²c) + t √((A−c)² + (1−t²)(1−c²)) ] / (A² − t²)
    // with A = 1 + (E/m_ec²)(1−c).  The radical's prefactor is the SIGNED t:
    // this selects the branch continuous and monotone in t with E′(0) = E_C.
    //
    // t may exceed +1 (deeply bound, relativistic inner-shell electrons); the
    // formula stays valid for t up to t_max < A.  On the negative side t < −1
    // can turn the discriminant negative, so clamp there — it only affects a
    // ~1e-6 probability tail of the deepest profiles.
    const double c = cos_theta;
    const double A = 1.0 + (energy_keV / kElectronMassKeV) * (1.0 - c);
    t = std::max(-0.99995, t);
    const double Ac = A - c;
    const double disc = Ac * Ac + (1.0 - t * t) * (1.0 - c * c);
    const double denom = A * A - t * t;
    if (disc < 0.0 || denom < 1e-12) {
        // Defensive: outside the formula's validity (caller's t_max truncation
        // should prevent this).  Cap at the incident energy; the caller clamps
        // to E − U.
        return energy_keV;
    }
    return energy_keV * (A - t * t * c + t * std::sqrt(disc)) / denom;
}

double doppler_t_max(double energy_keV, double cos_theta, double U_keV) {
    // The t for which E′ = E − U at this angle:
    // t = [ E·E′(1−c)/m_ec² − (E − E′) ] / √(E² + E′² − 2 E E′ c)
    const double Ep = energy_keV - U_keV;
    const double c = cos_theta;
    const double denom2 = energy_keV * energy_keV + Ep * Ep
                          - 2.0 * energy_keV * Ep * c;
    if (denom2 <= 0.0) return 0.0;  // only possible for U=0, c=1 (no scatter)
    return (energy_keV * Ep * (1.0 - c) / kElectronMassKeV - U_keV)
           / std::sqrt(denom2);
}

double klein_nishina_total(double energy_keV) {
    double alpha = energy_keV / kElectronMassKeV;

    if (alpha < 1e-4) {
        // Low-energy limit: Thomson cross section
        return kThomsonXS;
    }

    // Full Klein-Nishina integrated cross section per electron (barns)
    // σ_KN = (3σ_T / 4) * { (1+α)/α³ * [2α(1+α)/(1+2α) - ln(1+2α)] + ln(1+2α)/(2α) - (1+3α)/(1+2α)² }
    double a = alpha;
    double a2 = 1.0 + 2.0 * a;
    double ln_a2 = std::log(a2);
    double a3 = a * a * a;

    double sigma = 0.75 * kThomsonXS * (
        (1.0 + a) / a3 * (2.0 * a * (1.0 + a) / a2 - ln_a2)
        + ln_a2 / (2.0 * a)
        - (1.0 + 3.0 * a) / (a2 * a2)
    );

    return sigma;
}

double sample_rayleigh_cos_theta(int Z, double energy_keV, std::mt19937_64& rng) {
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const CrossSectionData& xsdata = CrossSectionData::instance();

    // Sample momentum transfer x ∝ F(x,Z)^2 by inverse transform (no rejection),
    // then accept with the polarization factor (1+cos^2θ)/2.  This targets the same
    // PDF as the former uniform-cosθ rejection — dσ/dcosθ ∝ (1+cos^2θ) F(x)^2 —
    // but acceptance here is >= 0.5 (vs. the old loop's low acceptance for the
    // forward-peaked form factor), using the offline-generated sampling CDF.
    double x_max = energy_keV / kHcKeVAngstrom;   // x at the backscatter limit
    double cos_theta;
    do {
        double x = xsdata.sample_rayleigh_x(Z, x_max, rng);
        double r = x / x_max;                     // r = sin(θ/2) = sqrt((1-cosθ)/2)
        cos_theta = 1.0 - 2.0 * r * r;
        if (cos_theta < -1.0) cos_theta = -1.0;
        else if (cos_theta > 1.0) cos_theta = 1.0;
    } while (uniform(rng) > 0.5 * (1.0 + cos_theta * cos_theta));

    return cos_theta;
}

Eigen::Vector3d rotate_direction(const Eigen::Vector3d& direction,
                                  double cos_theta, double phi) {
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    double ux = direction.x();
    double uy = direction.y();
    double uz = direction.z();

    double uz2 = uz * uz;

    Eigen::Vector3d new_dir;

    if (std::abs(uz) > 0.99999) {
        // Direction is nearly along z-axis — use simplified rotation
        double sign = (uz > 0.0) ? 1.0 : -1.0;
        new_dir.x() = sin_theta * cos_phi;
        new_dir.y() = sign * sin_theta * sin_phi;
        new_dir.z() = sign * cos_theta;
    } else {
        // General rotation using Rodrigues-like formula
        double sin_theta_0 = std::sqrt(1.0 - uz2);
        double inv_sin_theta_0 = 1.0 / sin_theta_0;

        new_dir.x() = ux * cos_theta
                     + (ux * uz * cos_phi - uy * sin_phi) * sin_theta * inv_sin_theta_0;
        new_dir.y() = uy * cos_theta
                     + (uy * uz * cos_phi + ux * sin_phi) * sin_theta * inv_sin_theta_0;
        new_dir.z() = uz * cos_theta
                     - sin_theta_0 * cos_phi * sin_theta;
    }

    // Re-normalize to prevent drift
    new_dir.normalize();
    return new_dir;
}

/// Sample cos_theta from Klein-Nishina using Butcher & Messel's method.
/// This is the inner sampling without bound-electron correction.
///
/// Reference: Butcher & Messel, Nucl. Phys. 20 (1960) 15.
///
/// PROVENANCE: independently implemented here from the published 1960 method;
/// NOT derived from, ported from, or copied out of any other code base.  The
/// same published algorithm is used by GEANT4 (G4LowEPComptonModel,
/// G4KleinNishinaCompton), EGS, PENELOPE and MCNP, so the algebraic structure
/// -- and hence some variable naming (epsilon0, alpha1/alpha2, the rejection
/// function) -- necessarily resembles theirs.  GEANT4 is named throughout this
/// file only as the reference implementation we validate against.
///
/// The energy ratio epsilon = E'/E is sampled from the KN differential
/// cross section using composition + rejection. Two proposal distributions
/// cover the [epsilon0, 1] range:
///   Branch 1: epsilon = exp(-alpha1 * U)  (density ∝ 1/epsilon)
///   Branch 2: epsilon² uniform on [epsilon0², 1]  (density ∝ epsilon)
/// Both are rejected with g(epsilon) = 1 - epsilon*sin²θ/(1+epsilon²).
static double sample_kn_cos_theta(double alpha,
                                   std::uniform_real_distribution<double>& uniform,
                                   std::mt19937_64& rng)
{
    double cos_theta;

    if (alpha < 1e-4) {
        // Thomson limit: dσ/dΩ ∝ (1 + cos²θ)
        while (true) {
            double xi1 = uniform(rng);
            double xi2 = uniform(rng);
            cos_theta = 2.0 * xi1 - 1.0;
            double g = 0.5 * (1.0 + cos_theta * cos_theta);
            if (xi2 <= g) break;
        }
    } else {
        // Butcher-Messel method for Klein-Nishina
        double epsilon0 = 1.0 / (1.0 + 2.0 * alpha);
        double epsilon0_sq = epsilon0 * epsilon0;
        double alpha1 = -std::log(epsilon0);    // = ln(1 + 2α)
        double alpha2 = 0.5 * (1.0 - epsilon0_sq);

        while (true) {
            double epsilon, epsilon_sq;

            if (uniform(rng) * (alpha1 + alpha2) < alpha1) {
                // Branch 1: sample epsilon from density ∝ 1/epsilon
                epsilon = std::exp(-alpha1 * uniform(rng));
                epsilon_sq = epsilon * epsilon;
            } else {
                // Branch 2: sample epsilon² uniformly on [epsilon0², 1]
                epsilon_sq = epsilon0_sq + (1.0 - epsilon0_sq) * uniform(rng);
                epsilon = std::sqrt(epsilon_sq);
            }

            // cos_theta from Compton kinematics: epsilon = 1/(1 + alpha*(1-cosθ))
            double one_minus_cos = (1.0 - epsilon) / (epsilon * alpha);
            double sin_sq = one_minus_cos * (2.0 - one_minus_cos);

            // Rejection: accept with probability 1 - epsilon*sin²θ/(1+epsilon²)
            double g_reject = 1.0 - epsilon * sin_sq / (1.0 + epsilon_sq);
            if (uniform(rng) <= g_reject) {
                cos_theta = 1.0 - one_minus_cos;
                break;
            }
        }

        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
    }

    return cos_theta;
}

double sample_compton_cos_theta(double energy_keV, int Z, std::mt19937_64& rng)
{
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double alpha = energy_keV / kElectronMassKeV;
    double cos_theta;

    if (Z > 0) {
        // Bound-electron correction: sample from KN then reject with S(x,Z)/Z
        const CrossSectionData& xsdata = CrossSectionData::instance();
        const double Z_d = static_cast<double>(Z);

        while (true) {
            cos_theta = sample_kn_cos_theta(alpha, uniform, rng);

            // Momentum transfer: x = (E/hc) * sin(θ/2)  [Å⁻¹]
            double sin_half = std::sqrt(std::max(0.0, (1.0 - cos_theta) * 0.5));
            double x = (energy_keV / kHcKeVAngstrom) * sin_half;
            double S = xsdata.scattering_function_S(Z, x);
            if (uniform(rng) <= S / Z_d) break;  // accept
        }
    } else {
        // Free-electron Klein-Nishina (no binding correction)
        cos_theta = sample_kn_cos_theta(alpha, uniform, rng);
    }
    return cos_theta;
}

double compton_angular_pdf_unnorm(double energy_keV, int Z, double mu)
{
    mu = std::max(-1.0, std::min(1.0, mu));
    const double alpha = energy_keV / kElectronMassKeV;
    const double eps = 1.0 / (1.0 + alpha * (1.0 - mu));
    // dσ_KN/dμ ∝ ε² (ε + 1/ε − sin²θ)
    const double sin_sq = 1.0 - mu * mu;
    double f = eps * eps * (eps + 1.0 / eps - sin_sq);
    if (Z > 0) {
        const double sin_half = std::sqrt(std::max(0.0, (1.0 - mu) * 0.5));
        const double x = (energy_keV / kHcKeVAngstrom) * sin_half;
        f *= CrossSectionData::instance().scattering_function_S(Z, x);
    }
    return f;
}

namespace {

/// 32-point Gauss-Legendre nodes/weights on [-1, 1] (positive half; the
/// nodes are symmetric).
///
/// These are the standard quadrature abscissae/weights -- roots of P_32 and the
/// associated weights -- tabulated in Abramowitz & Stegun, "Handbook of
/// Mathematical Functions" (NBS, 1964), Table 25.4. Mathematical constants in
/// the public domain, computable from scratch to any precision.
constexpr int kGL32 = 16;
constexpr double kGL32_x[kGL32] = {
    0.0483076656877383162, 0.1444719615827964934, 0.2392873622521370745,
    0.3318686022821276497, 0.4213512761306353454, 0.5068999089322293900,
    0.5877157572407623290, 0.6630442669302152010, 0.7321821187402896804,
    0.7944837959679424070, 0.8493676137325699701, 0.8963211557660521240,
    0.9349060759377396892, 0.9647622555875064308, 0.9856115115452683354,
    0.9972638618494815635};
constexpr double kGL32_w[kGL32] = {
    0.0965400885147278006, 0.0956387200792748594, 0.0938443990808045654,
    0.0911738786957638847, 0.0876520930044038111, 0.0833119242269467552,
    0.0781938957870703065, 0.0723457941088485062, 0.0658222227763618468,
    0.0586840934785355471, 0.0509980592623761762, 0.0428358980222266807,
    0.0342738629130214331, 0.0253920653092620595, 0.0162743947309056706,
    0.0070186100094700966};

// Integrand for the norm in the substituted variable s = sin(theta/2):
// mu = 1 - 2 s^2, dmu = -4 s ds, so
//   integral_{-1}^{1} f(mu) dmu = integral_0^1 f(1 - 2 s^2) * 4 s ds.
// The substitution makes the momentum transfer x = (E/hc) * s LINEAR in the
// integration variable: it removes the sqrt(1-mu) cusp that S(x,Z) imprints
// on the mu-space integrand at mu -> 1 (which otherwise limits any
// quadrature to ~1e-4 at MeV energies).
double norm_integrand_s(double energy_keV, int Z, double s) {
    const double mu = 1.0 - 2.0 * s * s;
    return 4.0 * s * compton_angular_pdf_unnorm(energy_keV, Z, mu);
}

double gl32_integrate_s(double a, double b, double energy_keV, int Z) {
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    double sum = 0.0;
    for (int i = 0; i < kGL32; ++i) {
        sum += kGL32_w[i] *
               (norm_integrand_s(energy_keV, Z, c - h * kGL32_x[i]) +
                norm_integrand_s(energy_keV, Z, c + h * kGL32_x[i]));
    }
    return sum * h;
}

} // anonymous namespace

double compton_angular_norm(double energy_keV, int Z)
{
    // Exact-key thread-local cache: biasing is restricted to the first
    // Compton vertex, so every call is at the (constant) primary energy and
    // hits the cache after the first event per thread.
    thread_local std::map<std::pair<int, double>, double> cache;
    const auto key = std::make_pair(Z, energy_keV);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    // Composite 32-point GL in s = sin(theta/2), 16 uniform panels: the
    // remaining nonsmoothness is only the S(x,Z) interpolation kinks;
    // ~1e-6 relative (unit-tested against fine Simpson integration).
    double norm = 0.0;
    const double h = 1.0 / 16.0;
    for (int k = 0; k < 16; ++k) {
        norm += gl32_integrate_s(k * h, (k + 1) * h, energy_keV, Z);
    }

    if (cache.size() > 4096) cache.clear();  // unbounded-growth guard
    cache.emplace(key, norm);
    return norm;
}

ComptonScatterResult sample_compton_scatter(
    double energy_keV,
    const Eigen::Vector3d& direction,
    std::mt19937_64& rng,
    int Z,
    bool doppler)
{
    // Angle sampling, then energy/direction at that angle. RNG draw order
    // is unchanged from the original monolithic implementation.
    double cos_theta = sample_compton_cos_theta(energy_keV, Z, rng);
    return finish_compton_at_angle(energy_keV, direction, cos_theta, rng, Z,
                                   doppler, nullptr);
}

ComptonScatterResult finish_compton_at_angle(
    double energy_keV,
    const Eigen::Vector3d& direction,
    double cos_theta,
    std::mt19937_64& rng,
    int Z,
    bool doppler,
    const Eigen::Vector3d* new_dir_override)
{
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Compute scattered energy: free-electron Compton line, optionally
    // Doppler-broadened by the bound electron's momentum.
    double scattered_energy = compton_scattered_energy(energy_keV, cos_theta);
    double binding_keV = 0.0;

    const ElementData* elem = (doppler && Z > 0)
        ? &CrossSectionData::instance().element(Z) : nullptr;
    if (elem != nullptr && elem->num_compton_shells > 0) {
        // Select a subshell among those ionizable at this energy (U_i < E),
        // with probability proportional to occupancy.
        double occ_eligible = 0.0;
        for (int i = 0; i < elem->num_compton_shells; ++i) {
            if (elem->shell_binding_keV[i] < energy_keV)
                occ_eligible += elem->shell_occupancy[i];
        }
        int shell = -1;
        if (occ_eligible > 0.0) {
            double target = uniform(rng) * occ_eligible;
            double cum = 0.0;
            for (int i = 0; i < elem->num_compton_shells; ++i) {
                if (elem->shell_binding_keV[i] >= energy_keV) continue;
                cum += elem->shell_occupancy[i];
                shell = i;
                if (target <= cum) break;
            }
        }

        if (shell >= 0) {
            const double U = elem->shell_binding_keV[shell];
            const double J0 = elem->shell_J0[shell];
            const double t_max = doppler_t_max(energy_keV, cos_theta, U);
            const double F_max = compton_profile_cdf(J0, t_max / kAlphaFS);

            if (F_max > 1e-12) {
                double xi = uniform(rng) * F_max;
                double t = compton_profile_invcdf(J0, xi) * kAlphaFS;
                t = std::min(t, t_max);
                scattered_energy =
                    doppler_scattered_energy(energy_keV, cos_theta, t);
            } else {
                // Degenerate truncation (E barely above U): free line, capped.
                scattered_energy = std::min(scattered_energy, energy_keV - U);
            }
            scattered_energy =
                std::max(0.0, std::min(scattered_energy, energy_keV - U));
            binding_keV = U;
        }
        // No eligible shell (E below every binding energy): free-line fallback.
    }

    double deposited = std::max(0.0, energy_keV - scattered_energy - binding_keV);

    Eigen::Vector3d new_dir;
    if (new_dir_override) {
        // Direction supplied by the caller (angular importance sampling);
        // no azimuth draw.
        new_dir = *new_dir_override;
    } else {
        // Sample azimuthal angle uniformly
        double phi = kTwoPi * uniform(rng);
        new_dir = rotate_direction(direction, cos_theta, phi);
    }

    return {scattered_energy, deposited, cos_theta, new_dir, binding_keV};
}

} // namespace ceelo
