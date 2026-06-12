#ifndef MassAttenuationTool_imp_hpp
#define MassAttenuationTool_imp_hpp
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

#include <cmath>
#include <algorithm>
#include <type_traits>

#include "InterSpec/MassAttenuationTool.h"

namespace ceres
{
  /* dummy namespace for when using this file for only doubles, and ceres.h hasnt been included */
}

namespace MassAttenuation
{

/** Mass attenuation coefficient for a fractional atomic number.

 C1-continuous in atomic number: Catmull-Rom cubic interpolation of log(mu)
 across the integer-AN nodes {floor(an)-1 ... floor(an)+2}, with one-sided
 tangents at the Z=1 and Z=98 boundaries.  Exactly equals
 massAttenuationCoefficientElement(an,energy) at integer atomic numbers, and
 unlike linear interpolation, both the value and the derivative with respect
 to atomic number are continuous across integer atomic numbers - which
 matters when an optimizer is fitting atomic number using gradients.

 `T` may be `double`, or a `ceres::Jet<>` (in which case the derivative with
 respect to atomic number is preserved); when using Jets, include ceres.h
 before this header.

 Outside [1,98] the value is clamped to the boundary element (derivative zero).

 Returns mass attenuation coefficient in PhysicalUnits (divide by
 PhysicalUnits::cm2/PhysicalUnits::g for cm2/g); throws the same as
 #massAttenuationCoefficientElement for invalid energies.
 */
template<typename T>
T mass_atten_coef_frac_an( const T &atomic_number, const float energy )
{
  using namespace std;
  using namespace ceres;

  double an_scalar;
  if constexpr ( std::is_same_v<T, double> )
    an_scalar = atomic_number;
  else
    an_scalar = atomic_number.a;

  if( an_scalar <= sm_min_xs_atomic_number )
    return T( static_cast<double>( massAttenuationCoefficientElement( sm_min_xs_atomic_number, energy ) ) );

  if( an_scalar >= sm_max_xs_atomic_number )
    return T( static_cast<double>( massAttenuationCoefficientElement( sm_max_xs_atomic_number, energy ) ) );

  const int lower_an = std::clamp( static_cast<int>( std::floor(an_scalar) ),
                                   sm_min_xs_atomic_number, sm_max_xs_atomic_number - 1 );
  const int upper_an = lower_an + 1;

  // Interpolate in log(mu): mu varies roughly geometrically with Z, and log
  //  keeps the interpolation positive by construction.
  const double log_mu_1 = std::log( massAttenuationCoefficientElement( lower_an, energy ) );
  const double log_mu_2 = std::log( massAttenuationCoefficientElement( upper_an, energy ) );

  // Catmull-Rom tangents; one-sided at the table boundaries
  double m1, m2;
  if( lower_an > sm_min_xs_atomic_number )
    m1 = 0.5*(log_mu_2 - std::log( massAttenuationCoefficientElement( lower_an - 1, energy ) ));
  else
    m1 = log_mu_2 - log_mu_1;

  if( upper_an < sm_max_xs_atomic_number )
    m2 = 0.5*(std::log( massAttenuationCoefficientElement( upper_an + 1, energy ) ) - log_mu_1);
  else
    m2 = log_mu_2 - log_mu_1;

  const T t = atomic_number - static_cast<double>( lower_an );  //preserves any Jet derivative
  const T t2 = t*t;
  const T t3 = t2*t;

  // Cubic Hermite basis
  const T h00 = 2.0*t3 - 3.0*t2 + 1.0;
  const T h10 = t3 - 2.0*t2 + t;
  const T h01 = -2.0*t3 + 3.0*t2;
  const T h11 = t3 - t2;

  const T log_mu = h00*log_mu_1 + h10*m1 + h01*log_mu_2 + h11*m2;

  return exp( log_mu );
}//T mass_atten_coef_frac_an( const T &atomic_number, const float energy )

}//namespace MassAttenuation

#endif //MassAttenuationTool_imp_hpp
