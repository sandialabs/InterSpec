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

#include <map>
#include <string>
#include <stdexcept>

#include "SandiaDecay.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"

using namespace std;


namespace CeeLoUtils
{

ceelo::MaterialSpec to_ceelo_material( const Material &mat )
{
  ceelo::MaterialSpec spec;
  spec.name = mat.name;
  spec.density_g_per_cm3 = mat.density * PhysicalUnits::cm3 / PhysicalUnits::g;

  map<int,double> frac_by_z;
  for( const Material::ElementFractionPair &efp : mat.elements )
    frac_by_z[efp.first->atomicNumber] += efp.second;
  for( const Material::NuclideFractionPair &nfp : mat.nuclides )
    frac_by_z[nfp.first->atomicNumber] += nfp.second;

  double sum = 0.0;
  for( const auto &zf : frac_by_z )
    sum += zf.second;
  if( sum <= 0.0 )
    throw runtime_error( "Material '" + mat.name + "' has no composition." );

  for( const auto &zf : frac_by_z )
  {
    if( (zf.first < 1) || (zf.first > 92) )
      throw runtime_error( "Material '" + mat.name + "' contains an element"
                           " (Z=" + std::to_string(zf.first) + ") the Monte"
                           " Carlo has no cross-section data for." );
    ceelo::MaterialComponent c;
    c.Z = static_cast<uint8_t>( zf.first );
    c.mass_fraction = zf.second / sum;
    spec.composition.push_back( c );
  }

  return spec;
}//to_ceelo_material(...)

}//namespace CeeLoUtils
