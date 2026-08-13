/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC.
 This library is distributed under the GNU Lesser General Public License,
 version 2.1 or (at your option) any later version.
 */

/* Direct EPICS2023 EADL atomic-relaxation data. See sources.lock.json. */
#pragma once

#include "cross_sections/CrossSectionData.h"

namespace ceelo {

inline constexpr int kRelaxationMaxZ = 99;
extern const FluorescenceData g_epics_k_relaxation[kRelaxationMaxZ + 1];
extern const LFluorescenceData g_epics_l_relaxation[kRelaxationMaxZ + 1];

} // namespace ceelo
