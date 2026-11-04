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

/// @file RayleighAngular.h
/// @brief Deterministic EVALUATION of the Rayleigh angular density.
///
/// `ComptonScatter.h` provides both a sampler (`sample_rayleigh_cos_theta`) and,
/// for Compton, an evaluator + normalization pair
/// (`compton_angular_pdf_unnorm` / `compton_angular_norm`). Rayleigh had only
/// the sampler. A next-event estimator cannot sample — it must evaluate the
/// differential cross-section in the one direction that points at the detector —
/// so this file supplies the missing Rayleigh half, mirroring the Compton
/// pattern exactly.

namespace ceelo {

/// Unnormalized Rayleigh angular density in mu = cos(theta):
///   f(mu) = (1 + mu^2) * F^2(x(mu,E), Z),   x = (E/hc) * sin(theta/2)
/// This is exactly the density `sample_rayleigh_cos_theta()` draws from, up to
/// the normalization `rayleigh_angular_norm()`.
double rayleigh_angular_pdf_unnorm(double energy_keV, int Z, double mu);

/// Integral of `rayleigh_angular_pdf_unnorm` over mu in [-1, 1].
///
/// Uses the same composite 32-point Gauss-Legendre quadrature in
/// s = sin(theta/2) that `compton_angular_norm()` uses. The substitution
/// mu = 1 - 2s^2 makes the momentum transfer x = (E/hc)*s linear in the
/// integration variable, removing the sqrt(1-mu) cusp the form factor imprints
/// at mu -> 1. That matters MORE here than for Compton: F^2 is more sharply
/// forward-peaked than S, so the cusp is stronger.
double rayleigh_angular_norm(double energy_keV, int Z);

} // namespace ceelo
