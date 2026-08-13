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

/// @file element_data.h
/// @brief Extern declarations for compiled-in element support data.
///
/// Generated from locked public-domain NIST EPQ bremsstrahlung tables and
/// BSD-licensed xraylib Compton-profile support. Photon interactions and atomic
/// relaxation are generated separately from EPICS2023 EPDL/EADL.
///
/// IMPORTANT: This header is included by CrossSectionData.h — do NOT include
/// any project headers here to avoid circular dependencies. The struct
/// the ElementData definition comes from CrossSectionData.h.

#include <cstdint>

namespace ceelo {

// Forward declarations — the full definitions are in CrossSectionData.h,
// which includes this header. We rely on the include order:
//   CrossSectionData.h  ->  defines ElementData
//   CrossSectionData.cpp ->  #include "element_data.h"
// So by the time element_data.h is processed, the structs are already defined.

/// Cross-section data for all 92 elements (Z=1..92).
/// Indexed by Z-1 (i.e., g_element_data[0] is hydrogen, g_element_data[91] is uranium).
/// Photon processes use their own compact tables and are queried through
/// CrossSectionData rather than exposed as raw pointers here.
extern const ElementData g_element_data[92];

/// Standard atomic weights for all 92 elements (g/mol).
/// Indexed by Z-1.
extern const double g_atomic_weights[92];

/// Seltzer-Berger bremsstrahlung spectral shape data (shared grids).
/// Defined in element_data.cpp (auto-generated).
extern const uint16_t kSB_n_kappa;       ///< Number of k/T fraction grid points (32)
extern const uint16_t kSB_n_energy;      ///< Number of electron energy grid points (27)
extern const float kSB_kappa[];          ///< k/T fraction values, ascending [kSB_n_kappa]
extern const float kSB_log_E_keV[];      ///< log10(electron KE / keV), ascending [kSB_n_energy]
extern const float kSB_chi_scale[kMaxZ]; ///< Per-element uint16 decode scales

} // namespace ceelo
