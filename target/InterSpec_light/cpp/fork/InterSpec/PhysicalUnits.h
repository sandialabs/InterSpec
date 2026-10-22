#ifndef PhysicalUnits_h
#define PhysicalUnits_h
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

/* InterSpec-light: the small subset of InterSpec's PhysicalUnits the forked code uses. */

#include <string>

namespace PhysicalUnits
{
  static const double MeV = 1.0E3;
  static const double  eV = 1.e-6*MeV;
  static const double keV = 1.e-3*MeV;

  /** Prints value and uncertainty to the given number of significant figures,
   e.g. (1.23457, 0.0123, 3) -> "1.23 ± 0.0123".
   */
  std::string printValueWithUncertainty( double value, double uncert, size_t nsigfig );
}//namespace PhysicalUnits

#endif //PhysicalUnits_h
