#ifndef DetectorPeakResponse_h
#define DetectorPeakResponse_h
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

/* InterSpec-light: reduced DetectorPeakResponse - only the FWHM functional forms the forked
 peak-fitting code uses.  The light version never has a DRF loaded, so `isValid()` and
 `hasResolutionInfo()` are always false.
 */

#include <cmath>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <type_traits>

#include "InterSpec/PhysicalUnits.h"

class DetectorPeakResponse
{
public:
  enum ResolutionFnctForm
  {
    kGadrasResolutionFcn, //See peakResolutionFWHM() implementation

    /** FWHM = sqrt( Sum_i{A_i*pow(x/1000,i)} ); */
    kSqrtPolynomial,

    /** FWHM = `sqrt(A0 + A1*E + A2/E)` */
    kSqrtEnergyPlusInverse,

    /** FWHM = `A0 + A1*sqrt(E)` */
    kConstantPlusSqrtEnergy,

    kNumResolutionFnctForm
  };//enum ResolutionFnctForm

  bool isValid() const { return false; }
  bool hasResolutionInfo() const { return false; }

  float peakResolutionSigma( const float ) const
  {
    throw std::runtime_error( "DetectorPeakResponse: no resolution info in InterSpec-light" );
  }

  static float peakResolutionFWHM( float energy, ResolutionFnctForm fcnFrm,
                                   const std::vector<float> &pars )
  {
    return static_cast<float>( peakResolutionFWHM( static_cast<double>(energy), fcnFrm,
                                                   pars.data(), pars.size() ) );
  }

  static float peakResolutionSigma( const float energy, ResolutionFnctForm fcnFrm,
                                    const std::vector<float> &pars )
  {
    return peakResolutionFWHM( energy, fcnFrm, pars ) / 2.35482f;
  }

  /** FWHM functional forms, for double or `ceres::Jet<>` T; energy clamped to >= 10 keV. */
  template<typename T, typename ParT>
  static T peakResolutionFWHM( T energy, const ResolutionFnctForm fcnFrm,
                               const ParT * const pars, const size_t num_pars );

  template<typename T>
  static T positive_c1_continuation( const T &value, const double floor );

  template<typename T>
  static T upper_c1_continuation( const T &value, const T &upper );
};//class DetectorPeakResponse


namespace DetectorPeakResponseImp
{
  /** The value part of a double or a ceres::Jet<>. */
  template<typename T>
  inline double scalar_value( const T &v )
  {
    if constexpr( std::is_arithmetic_v<T> )
      return static_cast<double>( v );
    else
      return static_cast<double>( v.a );
  }
}//namespace DetectorPeakResponseImp


template<typename T>
T DetectorPeakResponse::positive_c1_continuation( const T &value, const double floor )
{
  if( DetectorPeakResponseImp::scalar_value(value) >= floor )
    return value;

  // Equals `value` in both value and derivative at `floor`, stays strictly positive for every
  //  finite trial, and approaches zero monotonically as the raw argument tends to -infinity.
  return T(floor*floor) / (T(2.0*floor) - value);
}//positive_c1_continuation(...)


template<typename T>
T DetectorPeakResponse::upper_c1_continuation( const T &value, const T &upper )
{
  const T join = T(0.5) * upper;
  if( DetectorPeakResponseImp::scalar_value(value) <= 0.5*DetectorPeakResponseImp::scalar_value(upper) )
    return value;

  // Exact below half of `upper`, C1 at the join, and monotonically approaches (but never
  //  reaches) `upper`.  Physical values are far below the join; this only regularizes
  //  otherwise-invalid optimizer trials.
  const T distance = upper - join;
  return upper - distance*distance / (value - join + distance);
}//upper_c1_continuation(...)


template<typename T, typename ParT>
T DetectorPeakResponse::peakResolutionFWHM( T energy, const ResolutionFnctForm fcnFrm,
                                            const ParT * const pars, const size_t num_pars )
{
  using std::log;
  using std::pow;
  using std::sqrt;

  // Below ~10 keV the forms are physically ill-defined (e.g. kSqrtEnergyPlusInverse divides by
  //  energy) and carry no useful information for our use cases; treat as 10 keV.
  if( DetectorPeakResponseImp::scalar_value(energy) < 10.0*PhysicalUnits::keV )
    energy = T( 10.0*PhysicalUnits::keV );

  switch( fcnFrm )
  {
    case kGadrasResolutionFcn:
    {
      if( num_pars != 3 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionFWHM(): pars not defined" );

      // Straight-forward translation of the GADRAS Fortran GetFWHM ("form C", shared with
      //  PeakDists::gadras_fwhm).  The sign of `a` selects the low-energy branch: a >= 0 adds
      //  a linear "FWHM at zero energy" offset in quadrature; a < 0 instead bends the power law.
      const T a = T( pars[0] );   // resolution offset ("FWHM @ 0")
      const T b = T( pars[1] );   // resolution @ 661 (percent)
      const T c = T( pars[2] );   // resolution power

      if( DetectorPeakResponseImp::scalar_value(energy) > 661.0 )
        return 6.61 * b * pow( energy/661.0, c );

      if( DetectorPeakResponseImp::scalar_value(a) >= 0.0 )
      {
        // a >= 0 here, so fabs(a) == a; (661 - energy) >= 0 since energy <= 661 in this branch.
        T zero_limit = a * (661.0 - energy) / 661.0;
        if( DetectorPeakResponseImp::scalar_value(zero_limit) < 0.0 )
          zero_limit = T( 0.0 );
        const T fwhm = 6.61 * b * pow( energy/661.0, c );
        return sqrt( zero_limit*zero_limit + fwhm*fwhm );
      }//if( a >= 0.0 )

      T e_clamped = energy;
      if( DetectorPeakResponseImp::scalar_value(e_clamped) < 30.0 )
        e_clamped = T( 30.0 );
      const T p = pow( c, T(1.0)/log(1.0 - a) );
      return 6.61 * b * pow( e_clamped/661.0, p );
    }//case kGadrasResolutionFcn:

    case kSqrtEnergyPlusInverse:
    {
      if( num_pars != 3 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionFWHM(): pars not defined" );
      energy /= PhysicalUnits::keV;

      const T width_squared = T(pars[0]) + T(pars[1])*energy + T(pars[2])/energy;
      return sqrt( positive_c1_continuation(width_squared, 1.0e-12) );
    }//case kSqrtEnergyPlusInverse:

    case kConstantPlusSqrtEnergy:
    {
      if( num_pars != 2 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionFWHM(): pars not defined" );
      energy /= PhysicalUnits::keV;

      return T(pars[0]) + T(pars[1])*sqrt(energy);
    }//case kConstantPlusSqrtEnergy:

    case kSqrtPolynomial:
    {
      if( num_pars < 1 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionFWHM(): pars not defined" );

      energy /= PhysicalUnits::MeV;

      // Horner's rule - more stable than summing powers
      T val = T( pars[num_pars - 1] );
      for( int i = static_cast<int>(num_pars) - 2; i >= 0; i -= 1 )
        val = val * energy + T( pars[i] );

      return sqrt( positive_c1_continuation(val, 1.0e-12) );
    }//case kSqrtPolynomial:

    case kNumResolutionFnctForm:
      throw std::runtime_error( "DetectorPeakResponse::peakResolutionFWHM(): Resolution not defined" );
      break;
  }//switch( fcnFrm )

  assert( 0 );
  return T( 0.0 );
}//peakResolutionFWHM(...)

#endif //DetectorPeakResponse_h
