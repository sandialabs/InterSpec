#pragma once
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

/// @file ComptonScatter.h
/// @brief Compton scattering sampling using Butcher-Messel method for
///        Klein-Nishina distribution, with optional bound-electron correction.
///
/// Given the incident photon energy, samples:
///   - Scattering angle (cos_theta)
///   - Scattered photon energy
///   - New photon direction (requires the incident direction)

#include <Eigen/Core>
#include <random>

namespace ceelo {

/// Result of a Compton scatter event.
struct ComptonScatterResult {
    double scattered_energy_keV;  ///< Energy of the scattered photon
    double deposited_energy_keV;  ///< Recoil electron kinetic energy
    double cos_theta;             ///< Cosine of the scattering angle
    Eigen::Vector3d new_direction; ///< New photon direction (normalized)
    double binding_deposit_keV = 0.0; ///< Subshell binding energy U_i, deposited
                                      ///< locally (Doppler sampling only)
};

/// Sample a Compton scatter from the Klein-Nishina distribution.
///
/// Uses the Butcher-Messel method (same as GEANT4 G4LowEPComptonModel).
/// The scattering angle is sampled from the Klein-Nishina differential cross section:
///   dσ/dΩ = (r_e²/2) × (E'/E)² × (E/E' + E'/E - sin²θ)
/// where E' = E / (1 + (E/m_e c²)(1 - cosθ))
///
/// If Z > 0, applies bound-electron correction using the incoherent scattering
/// function S(x,Z). The effective differential cross section becomes:
///   dσ/dΩ = (dσ_KN/dΩ) × S(x,Z) / Z
/// where x = (E/hc) × sin(θ/2) is the momentum transfer parameter.
///
/// If doppler is true (and Z > 0), the scattered photon energy at the sampled
/// angle is additionally Doppler-broadened by the bound electron's momentum
/// (impulse approximation with PENELOPE-style analytic one-parameter subshell
/// profiles).  The subshell binding energy U_i is returned in
/// binding_deposit_keV for local deposition (no Compton-vacancy relaxation),
/// so E' + KE + U_i == E exactly.
///
/// @param energy_keV  Incident photon energy in keV
/// @param direction   Current photon direction (normalized)
/// @param rng         Random number generator
/// @param Z           Atomic number for bound-electron correction (0 = free electron)
/// @param doppler     Enable Doppler broadening of the scattered energy
/// @return            Scatter result with new energy, direction, and deposited energy
ComptonScatterResult sample_compton_scatter(
    double energy_keV,
    const Eigen::Vector3d& direction,
    std::mt19937_64& rng,
    int Z = 0,
    bool doppler = false);

/// @name Compton-profile helpers (exposed for unit testing)
/// PENELOPE-style analytic per-electron Compton profile
///   J(p) = J0 (1 + 2 J0 |p|) exp{ ½[1 − (1 + 2 J0 |p|)²] },  p in atomic units.
///@{

/// CDF of the analytic profile: F(p), with F(0) = 1/2.
double compton_profile_cdf(double J0, double p_au);

/// Inverse CDF: p such that F(p) = xi, xi in (0,1).
double compton_profile_invcdf(double J0, double xi);

/// Doppler-shifted scattered photon energy (impulse approximation).
/// t = p_z/(m_e c) is the projection of the electron's pre-collision momentum
/// on the scattering vector; t = 0 recovers the free-electron Compton line.
double doppler_scattered_energy(double energy_keV, double cos_theta, double t);

/// Maximum t allowed by kinematics for subshell binding energy U_keV
/// (the t for which E' = E - U at the given angle).
double doppler_t_max(double energy_keV, double cos_theta, double U_keV);
///@}

/// Sample only the Compton scattering cosine from KN x S(x,Z)/Z (the same
/// rejection chain sample_compton_scatter() uses, with identical RNG draw
/// order). Z = 0 gives the free-electron Klein-Nishina angle.
double sample_compton_cos_theta(double energy_keV, int Z, std::mt19937_64& rng);

/// Complete a Compton scatter at a GIVEN angle: Doppler-broadened scattered
/// energy (if doppler && Z>0), recoil energy, and new direction. This is the
/// post-angle half of sample_compton_scatter(); with new_dir_override ==
/// nullptr it draws the azimuth uniformly (same RNG order as
/// sample_compton_scatter). A non-null new_dir_override supplies the new
/// direction explicitly (no azimuth draw) — used by angular importance
/// sampling, where the direction was sampled on the sphere and cos_theta =
/// new_dir.dot(direction).
ComptonScatterResult finish_compton_at_angle(
    double energy_keV,
    const Eigen::Vector3d& direction,
    double cos_theta,
    std::mt19937_64& rng,
    int Z,
    bool doppler = false,
    const Eigen::Vector3d* new_dir_override = nullptr);

/// Unnormalized Compton angular density in mu = cos(theta):
///   f(mu) = (dsigma_KN/dmu)(mu; E) * S(x(mu,E), Z)      (S omitted for Z=0)
/// This is exactly the density sample_compton_cos_theta() draws from, up to
/// the normalization compton_angular_norm(). Used to weight angular
/// importance sampling of source-geometry Compton scatters.
double compton_angular_pdf_unnorm(double energy_keV, int Z, double mu);

/// Integral of compton_angular_pdf_unnorm over mu in [-1, 1], by composite
/// Gauss-Legendre quadrature (two 32-point panels split inside the forward
/// peak). Relative accuracy ~1e-6 (unit-tested); cached per (Z, energy) in a
/// thread-local map (exact-key, so repeated calls at the primary energy are
/// free).
double compton_angular_norm(double energy_keV, int Z);

/// Compute the Klein-Nishina cross section per electron at the given energy.
/// Returns the total integrated cross section in barns.
double klein_nishina_total(double energy_keV);

/// Given a scattering cosine, compute the scattered photon energy.
/// E' = E / (1 + alpha * (1 - cos_theta))
/// where alpha = E / (m_e c²)
double compton_scattered_energy(double energy_keV, double cos_theta);

/// Sample cos(theta) for Rayleigh (coherent) scattering using form-factor
/// rejection sampling.  dσ/dΩ ∝ F²(x,Z) × (1+cos²θ)/2.
///
/// @param Z          Atomic number of the scattering element
/// @param energy_keV Photon energy in keV
/// @param rng        Random number generator
/// @return           cos(theta) of the scattered photon
double sample_rayleigh_cos_theta(int Z, double energy_keV, std::mt19937_64& rng);

/// Rotate a direction vector by polar angle theta and uniform azimuthal angle phi.
/// Used to apply scattering angles to the photon direction.
Eigen::Vector3d rotate_direction(const Eigen::Vector3d& direction,
                                  double cos_theta, double phi);

} // namespace ceelo
