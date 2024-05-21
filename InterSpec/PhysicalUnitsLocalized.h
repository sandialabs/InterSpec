//Note these units are NOT CLHEP compatible!

#ifndef PhysicalUnitsLocalized_h
#define PhysicalUnitsLocalized_h 1
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

#include <vector>
#include <string>
#include <utility>

#include "InterSpec/PhysicalUnits.h"

namespace Wt
{
  class WMessageResourceBundle;
}

/** This file patches writing units to/from strings, to take into account the users current locale.
*/
namespace PhysicalUnitsLocalized
{
  double stringToTimeDuration( std::string str, double second_def = PhysicalUnits::second );
  double stringToTimeDurationPossibleHalfLife( const std::string &str,
                                               const double halflife,
                                                double second_def = PhysicalUnits::second );
  std::string printToBestTimeUnits( double time,
                                    int maxNpostDecimal = 2,
                                    double second_definition = PhysicalUnits::second );
  
  
  // Some definitions to allow unit testing
  double stringToTimeDurationPossibleHalfLife( const std::string &str,
                                               const double halflife,
                                                double second_def = PhysicalUnits::second,
                                              Wt::WMessageResourceBundle &bundle );
  std::string printToBestTimeUnits( double time,
                                    int maxNpostDecimal = 2,
                                    double second_definition = PhysicalUnits::second,
                                   Wt::WMessageResourceBundle &bundle );
 
}//namespace PhysicalUnits
#endif  //PhysicalUnitsLocalized_h
