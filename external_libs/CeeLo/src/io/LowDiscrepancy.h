#ifndef CEELO_IO_LOW_DISCREPANCY_H
#define CEELO_IO_LOW_DISCREPANCY_H
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
/// @file LowDiscrepancy.h
/// @brief Deterministic low-discrepancy sequences shared by the response
/// generator (probe banks) and the detector-side etendue line sampler.

#include <cstdint>

namespace ceelo {

/// Halton low-discrepancy value: index i (>= 0), prime base. The i-th point of
/// a d-dimensional Halton set is (halton(i,2), halton(i,3), halton(i,5), ...).
/// Deterministic; consecutive indices fill the unit interval evenly, so a
/// contiguous index range is a usable quadrature set at any size.
inline double halton(uint64_t i, uint32_t base) {
    double f = 1.0, r = 0.0;
    for (uint64_t n = i + 1; n > 0; n /= base) {
        f /= base;
        r += f * static_cast<double>(n % base);
    }
    return r;
}

}  // namespace ceelo

#endif  // CEELO_IO_LOW_DISCREPANCY_H
