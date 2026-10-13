#ifndef MassAttenuationTool_h
#define MassAttenuationTool_h
/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
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

#include "InterSpec_config.h"

/** Photon mass attenuation coefficients of the elements.

 A thin forward to CeeLo's compiled-in EPICS2023 (EPDL) photon cross sections
 (external_libs/CeeLo/src/cross_sections/CrossSectionData.h), converted from barns/atom to mass attenuation with
 CeeLo's (xraylib 4.2.1) atomic weights, so InterSpec's attenuation is the same as CeeLo's Monte Carlo transport.
 The data are immutable and compiled in: thread safe, and no file I/O.  CeeLo's photoelectric + Compton + pair agree
 with NIST XCOM to within 0.3% for Z 92-98 from 15 keV to 10 MeV, and 0.5% from 12 to 20 MeV for the elements compared
 (see external_libs/CeeLo/DESIGN.md).

 These replaced the old data/em_xs_data tables (1 keV - 100 MeV) in Oct 2026.  Measured against them, per
 gram, for Z = 1-98 on 4000 log-spaced energies from 10 keV to 20 MeV:
   - Total (photoelectric + Compton + pair), away from absorption edges: median 0.03%, p95 0.16%, p99 1.3%.  Water,
     air, concrete, Fe, Pb, and U at common gamma lines from 60 keV to 15 MeV (all away from edges) moved by under 0.1%.
     The ~1% outliers are atomic weights: Pm -1.3% and Tc -1.0% (xraylib's 147 and 99), and Tm +1.0% (xraylib lists Tm
     at Er's 167.27, not 168.93).
   - Absorption edges: em_xs_data interpolated across them, so it stayed >1% off for a median 1% (worst 3.9%, Nd) in
     energy above each K edge - e.g., Pb's 88.01 keV edge ramped over 88.2-90.6 keV, and U's sat at ~116.1 instead of
     115.61 keV.  CeeLo's edges are sharp at the EADL energies; within those bands mu/rho changed by up to ~5.6x.
   - Rayleigh (not part of the total): median 2.9%, up to ~23% just below heavy-element K edges.
   - CPU (Apple M1 Pro, Release): ~80 ns per total-coefficient call vs ~91 ns for em_xs_data (random Z and energy;
     74 vs 81 ns cycling the same 240 element/energy pairs, like a fit does); fractional-AN calls are ~355 vs ~335 ns.
     There is no first-use cost; em_xs_data took 7-22 ms to read all 98 element files.
 */
namespace MassAttenuation
{
  static const int sm_min_xs_atomic_number = 1;
  static const int sm_max_xs_atomic_number = 98;

  /** The energy range, in keV, that the cross-section data covers; energies outside of this range are evaluated at the
   nearest end of it.
   */
  static constexpr double sm_min_xs_energy_keV = 10.0;
  static constexpr double sm_max_xs_energy_keV = 20000.0;

  enum class GammaEmProcces : int
  {
    ComptonScatter,
    RayleighScatter,
    PhotoElectric,
    PairProduction,
    NumGammaEmProcces
  };//enum GammaEmProcces


  /** Gives the mass attenuation coefficient of an element for photoelectric + Compton + pair production (coherent, i.e.
   Rayleigh, scattering is not included), in units of PhysicalUnits; that is, divide by
   (PhysicalUnits::cm2 / PhysicalUnits::g) to get cm2/g.

   Energies outside [sm_min_xs_energy_keV, sm_max_xs_energy_keV] are clamped into that range (CeeLo's own behavior), so
   below 10 keV the attenuation is under-estimated; e.g., air at 5 keV gives 4.9 instead of 39.7 cm2/g.

   Throws std::runtime_error if atomic_number is outside [sm_min_xs_atomic_number, sm_max_xs_atomic_number], or
   energy is not finite.  (The em_xs_data implementation instead returned 0, i.e. no attenuation, below ~1 keV.)

   \param atomic_number Atomic number ranging from 1 to 98, inclusive
   \param energy Energy, in keV
   */
  double massAttenuationCoefficientElement( const int atomic_number, const double energy );

  /** Similar to the other #massAttenuationCoefficientElement function, but instead
   * only for a specific sub-proccess.
   */
  double massAttenuationCoefficientElement( const int atomic_number, const double energy, MassAttenuation::GammaEmProcces process );

  /** Similar to #massAttenuationCoefficientElement, but for fractional atomic
   * numbers: a C1-continuous (Catmull-Rom cubic in log(mu) vs atomic number)
   * interpolation across the neighboring integer atomic numbers - see
   * #mass_atten_coef_frac_an in MassAttenuationTool_imp.hpp, which this
   * forwards to, and which can also preserve ceres::Jet derivatives.
   * The cubic is not monotone: where an absorption edge falls between the
   * neighboring elements at this energy (e.g., Z 91-92 at 115.7 keV), it can
   * overshoot the bracketing elements' values by up to ~14%.
   */
  double massAttenuationCoefficientFracAN( const double atomic_number, const double energy );
}//namespace MassAttenuation


#endif  //MassAttenuationTool_h
