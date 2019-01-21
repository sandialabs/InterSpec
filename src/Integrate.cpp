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

#include "InterSpec/Integrate.h"

#include <vector>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "Cuba-3.0/cuba.h"


using namespace std;


namespace Integrate
{
void CuhreIntegrate( const int ndim,
                       Integrate::Integrand integrand,
                       void *userdata,
                       const double epsrel,
                       const double epsabs,
                       const unsigned int flags,
                       const size_t mineval,
                       const size_t maxeval,
                       int &pnregions,
                       int &pneval,
                       int &pfail,
                       double &integral,
                       double &error,
                       double &prob )
{
  const int key = 0;
  
  const int ncomp = 1;
  Cuhre( ndim, ncomp, integrand, userdata,
         epsrel, epsabs, flags, mineval, maxeval, key,
         &pnregions, &pneval, &pfail,
         &integral, &error, &prob );
}//CuhreIntegrate(...)

}//namespace Integrate
