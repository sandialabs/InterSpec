#ifndef LightMath_h
#define LightMath_h
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
#include <cstdint>
#include <numbers>
#include <utility>
#include <stdexcept>

/** The few Boost.Math pieces the forked code used, so the light build does not need Boost.

 `erf_inv`, `erfc_inv`, `normal_quantile`, and `bisect` are ports of Boost 1.84's double-precision
 implementations (same algorithms, constants, and order of operations, so the same results), and
 throw the same standard exception types for the same bad input.  The constants are the same
 correctly-rounded values as Boost's.
 */
namespace LightMath
{
  constexpr double pi = std::numbers::pi;
  constexpr double e = std::numbers::e;
  constexpr double ln_ten = std::numbers::ln10;
  constexpr double root_two = std::numbers::sqrt2;
  constexpr double root_pi = 1.772453850905516027298167483341145182797549456;
  constexpr double root_two_pi = 2.506628274631000502415765284811045253006986740;
  constexpr double root_half_pi = 1.253314137315500251207882642405522626503493370;
  constexpr double one_div_root_two = 0.707106781186547524400844362104849039284835938;

  /** Inverse error function; throws std::domain_error outside [-1,1], and std::overflow_error at +-1. */
  double erf_inv( const double z );

  /** Inverse complementary error function; throws std::domain_error outside [0,2], and
   std::overflow_error at 0 and 2. */
  double erfc_inv( const double z );

  /** Quantile of the standard normal distribution (boost::math::normal_distribution<double>(0,1)). */
  double normal_quantile( const double p );

  /** Quantile of a Poisson distribution, rounded "outwards" like Boost's default discrete-quantile
   policy: for p < 0.5 the largest k with P(X <= k) <= p, otherwise the smallest k with
   P(X <= k) >= p; 0 if p <= P(X = 0).  Sums the probabilities, so is for means up to ~1000.
   */
  double poisson_quantile( const double mean, const double p );

  /** Port of boost::math::tools::bisect: bisects [min, max], where f changes sign, until `tol(a, b)`
   is true, or `max_iter` evaluations; returns the final bracket, and sets `max_iter` to the number
   of evaluations used.  Throws std::runtime_error if min >= max, or f does not change sign.

   From boost/math/tools/roots.hpp: (C) Copyright John Maddock 2006.  Distributed under the Boost
   Software License, Version 1.0 (see the license text in LightMath.cpp).
   */
  template<class F, class Tol>
  std::pair<double,double> bisect( F f, double min, double max, Tol tol, std::uintmax_t &max_iter )
  {
    double fmin = f( min );
    const double fmax = f( max );
    if( fmin == 0.0 )
    {
      max_iter = 2;
      return { min, min };
    }
    if( fmax == 0.0 )
    {
      max_iter = 2;
      return { max, max };
    }

    if( min >= max )
      throw std::runtime_error( "Arguments in wrong order in bisect" );
    if( fmin * fmax >= 0.0 )
      throw std::runtime_error( "No change of sign in bisect: either there is no root to find, or there are multiple roots in the interval" );

    const auto sign = []( const double z ) -> int { return (z == 0.0) ? 0 : (std::signbit( z ) ? -1 : 1); };

    // Three function evaluations so far
    std::uintmax_t count = max_iter;
    count = (count < 3) ? 0 : (count - 3);

    while( count && !tol( min, max ) )
    {
      const double mid = (min + max) / 2;
      const double fmid = f( mid );
      if( (mid == max) || (mid == min) )
        break;
      if( fmid == 0.0 )
      {
        min = max = mid;
        break;
      }else if( sign( fmid ) * sign( fmin ) < 0 )
      {
        max = mid;
      }else
      {
        min = mid;
        fmin = fmid;
      }
      --count;
    }//while( count && !tol( min, max ) )

    max_iter -= count;
    return { min, max };
  }//bisect(...)
}//namespace LightMath

#endif //LightMath_h
