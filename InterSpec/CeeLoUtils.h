#ifndef CeeLoUtils_h
#define CeeLoUtils_h
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

#include "io/DetectorResponse.h"

struct Material;

/** Small shared helpers for mapping InterSpec objects onto the CeeLo
 Monte-Carlo library's types.  Used by the detector-geometry MC UI
 (DetectorGeometryInput / MakeMcResponseForDrf), the fixed-geometry response
 builder, and the cascade-summing truth-generation tests.
 */
namespace CeeLoUtils
{
  /** Converts an InterSpec MaterialDB material to a self-contained CeeLo
   material spec (per-element mass fractions; nuclide fractions folded into
   their element).  Throws for elements CeeLo has no data for (Z > 92).
   */
  ceelo::MaterialSpec to_ceelo_material( const Material &mat );
}//namespace CeeLoUtils

#endif //CeeLoUtils_h
