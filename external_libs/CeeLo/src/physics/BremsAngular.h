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

/// @file BremsAngular.h
/// @brief 2BS bremsstrahlung angular generator (sampler + matching evaluator).
///
/// The 2BS generator of Bielajew, Mohan, and Chui (PIRS-0203) -- the same
/// angular model G4LivermoreBremsstrahlungModel uses via G4Generator2BS.
/// Photon polar angle about the electron direction is expressed through
///     y in [0, y_max],   y_max = 2 beta (1 + beta) gamma^2,
///     cos(theta) = 1 - 2 y / y_max          (constant Jacobian dy/dmu = y_max/2).
///
/// The SAMPLER draws y from the analytically-invertible envelope density
/// h(y) = [(1 + y_max)/y_max] / (1 + y)^2 and accepts with probability
/// f(y)/M, where f is the 2BS rejection function and M = max(f(0), f(y_max))
/// (a flat ceiling; negative f is rejected outright). The realized density is
/// therefore
///     p(y)  PROPORTIONAL TO  min(max(f(y), 0), M) / (1 + y)^2
/// -- the 1/(1+y)^2 is part of the physical density, not a rejection
/// convenience. Any evaluator must carry it, and must apply the same clamps,
/// or it is evaluating a different distribution than the sampler draws
/// (without the (1+y)^-2 factor, wide angles are over-weighted by ~1e3-1e4 at
/// MeV energies).
///
/// The normalization exploits an exact cancellation: substituting the
/// envelope map y(q) = q y_max / (1 + y_max (1 - q)) into
/// N = integral of clamp(f(y)) / (1+y)^2 dy gives
///     N = [y_max / (1 + y_max)] * integral_0^1 clamp(f(y(q))) dq,
/// a bounded smooth integrand (single kink where f crosses zero), evaluated
/// with composite Gauss-Legendre. Sampler and evaluator share ONE rejection
/// function in BremsAngular.cpp so they cannot drift apart.
///
/// `atomic_number` is a double: callers pass the material effective Z
/// (ElectronCsda::effective_Z), matching how the sampler has always been used.

#include <Eigen/Core>

#include <cstdint>
#include <random>

namespace ceelo {

/// y_max = 2 beta (1 + beta) gamma^2 for the given electron kinetic energy.
/// Returns 0 for non-positive kinetic energy.
double brems_2bs_ymax(double electron_kinetic_energy_keV);

/// Sampler-faithful unnormalized density in y:
///     min(max(f(y), 0), M) / (1 + y)^2,   M = max(f(0), f(y_max)).
/// Zero outside [0, y_max] or for kinematically invalid (KE, k).
double brems_2bs_pdf_y_unnorm(double electron_kinetic_energy_keV,
                              double photon_energy_keV,
                              double atomic_number,
                              double y);

/// Normalization integral of brems_2bs_pdf_y_unnorm over [0, y_max].
/// Computed per call (emission events are rare; no table).
double brems_2bs_norm(double electron_kinetic_energy_keV,
                      double photon_energy_keV,
                      double atomic_number);

/// Normalized per-steradian density of the emitted photon direction about the
/// electron direction, at polar cosine `mu`:
///     p(Omega) = [pdf_y(y(mu)) * (y_max / 2) / N] / (2 pi),
/// with y(mu) = y_max (1 - mu) / 2. Integrates to 1 over the sphere.
/// Returns 0 when the kinematics are invalid (the sampler's fallback there is
/// a delta at the electron direction, which no density can represent; callers
/// in the next-event estimator never hit that regime).
double brems_2bs_pdf_sr(double electron_kinetic_energy_keV,
                        double photon_energy_keV,
                        double atomic_number,
                        double mu);

/// Sample the photon direction from the same 2BS angular generator used by
/// G4LivermoreBremsstrahlungModel (Bielajew, Mohan, and Chui, PIRS-0203).
/// Hoisted verbatim from ElectronCsda.cpp so the evaluator above shares its
/// rejection function; the RNG draw sequence is unchanged.
Eigen::Vector3d sample_brems_direction_2bs(
    const Eigen::Vector3d& electron_direction,
    double electron_kinetic_energy_keV,
    double photon_energy_keV,
    double atomic_number,
    std::mt19937_64& rng);

} // namespace ceelo
