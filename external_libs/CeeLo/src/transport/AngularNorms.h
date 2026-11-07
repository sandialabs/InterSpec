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

/// @file AngularNorms.h
/// @brief Interpolated angular-normalization tables for Compton and Rayleigh.
///
/// WHY THIS EXISTS. `compton_angular_norm()` and `rayleigh_angular_norm()` each
/// run a 512-point composite Gauss-Legendre quadrature and memoize the result
/// under an EXACT (Z, energy) key. That is the right design for the existing
/// caller — Compton-angle mixture biasing is restricted to the FIRST vertex, so
/// every call is at the constant primary energy and hits the cache after the
/// first event (this is precisely why `ComptonBiasConfig::max_biased_vertices`
/// defaults to 1).
///
/// The next-event estimator breaks that assumption. It scores at EVERY vertex of
/// a multiple-scatter walk, so the energy is different at essentially every call
/// and the exact-key cache never hits — the quadrature would move into the hot
/// path and dominate the run.
///
/// These tables replace the per-call quadrature with a lazily-built, cubic
/// log-log interpolated table per (process, Z), shared process-wide. Measured
/// worst-case relative error against the exact quadrature (tests/
/// test_env_scatter.cpp): ~2e-7 for Compton and ~2e-4 for Rayleigh, whose norm
/// inherits kinks from the tabulated form factor and so is grid-limited rather
/// than curvature-limited. Both are far below the ~10-20% accuracy target of the
/// scatter continuum itself.

namespace ceelo {

enum class ScatterProcess { Compton, Rayleigh };

/// Interpolated integral of the angular density over mu in [-1, 1].
///
/// Equivalent to `compton_angular_norm(energy_keV, Z)` /
/// `rayleigh_angular_norm(energy_keV, Z)` to <= 1e-4 relative inside the
/// tabulated range, and falls through to the exact quadrature outside it.
double angular_norm_interp(ScatterProcess process, double energy_keV, int Z);

/// Tabulated energy range (keV). Outside this the exact quadrature is used.
double angular_norm_table_min_keV();
double angular_norm_table_max_keV();

} // namespace ceelo
