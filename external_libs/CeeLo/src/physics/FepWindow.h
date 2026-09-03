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

/// @file FepWindow.h
/// @brief The one definition of "full energy": the photopeak acceptance window.

namespace ceelo {

/// Half-width of the full-energy-peak window, in keV.  A scoring event is
/// full-energy when |E_deposited - E_source| < this.
///
/// ONE definition, shared by everything that has to agree on what "in the peak"
/// means: the MC's FEP tally, the transport layer's FEP-only early kill, the
/// cascade peak windows, and the analytic in-window Compton fraction
/// (kn_in_window_fraction).  Agreement is not cosmetic -- scoring a window
/// wider than the early kill assumes silently under-counts FEP, and a model
/// whose in-window credit assumes a different window than the MC scored is
/// being validated against the wrong truth.
///
/// 0.75 keV suits HPGe.  Lower-resolution detectors want a wider window, and
/// the peak-dependent callers already do that themselves -- CascadeSummingCalc
/// sets each window from the fitted peak's own sigma.  Widening this default
/// for the MC is deliberately NOT wired up yet.
constexpr double kDefaultFepWindowKeV = 0.75;
} // namespace ceelo
