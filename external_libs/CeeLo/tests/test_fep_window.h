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

/// @file test_fep_window.h
/// @brief The full-energy window these tests are baselined against.

#include "physics/FepWindow.h"

namespace ceelo {

/// Every numeric expectation in this suite was established while the engine
/// scored full-energy with a 1.5 keV half-window, and the committed GEANT4
/// references in data/geant4_reference/ were generated at that window too.
///
/// The engine default is now narrower (kDefaultFepWindowKeV).  Rather than
/// re-baseline several hundred assertions against a value that is still being
/// chosen, the tests pin the window they were written for - so a failure means
/// the physics moved, not that the window did.  Re-baselining is a deliberate
/// job for whenever the default settles, and it has to move the GEANT4
/// references with it.
constexpr double kTestFepWindowKeV = 1.5;

} // namespace ceelo
