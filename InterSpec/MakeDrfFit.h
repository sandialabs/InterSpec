#ifndef MakeDrfFit_h
#define MakeDrfFit_h
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

#include <deque>
#include <vector>
#include <memory>

#include "InterSpec/DetectorPeakResponse.h"

class PeakDef;
class Measurement;

namespace MakeDrfFit
{
  /* If result passed into the function is of the proper size, then it will be
    used as the starting parameters for the fit, otherwise some default values
    will be chosen to start the fit with.
    The provided Measurement should be same one peaks where fit for.
    Currently for kSqrtPolynomial, only the first two coefficients (A1 and A2)
    are fit for.
   
    @param sqrtEqnOrder Only used if fnctnlForm==DetectorPeakResponse::kSqrtPolynomial
   
    @returns chi2 (not divided by dof).
   
    Throws exception on error with a kinda explanatory message.
  */
  double performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                             const size_t num_gamma_channels,
                             const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                             const int sqrtEqnOrder,
                             std::vector<float> &result,
                             std::vector<float> &uncerts );
  
  /** Fits the FWHM as a sqrt( Sum_i {A_i *pow(x,i)} )
   */
  double fit_fwhm_least_linear_squares( const std::deque< std::shared_ptr<const PeakDef> > &peaks,
                                       const int order,
                                       std::vector<float> &coeffs,
                                       std::vector<float> &coeff_uncerts );
  
  double peak_width_chi2( double predicted_sigma, const PeakDef &peak );
  
  
  
  struct DetEffDataPoint{ float energy, efficiency, efficiency_uncert; };
  
  /** Performs fit of data efficiencies to:
      eff(x) = exp(A + B*log(x) + C*log(x)^2 + ...)
   
      Only barely kinda superficially seems to be working; much work and testing
      to be done (20190512)
   
      Detects if DetEffDataPoint::energy is in keV or MeV by testing if largest
      energy is greater than 30; if so, in keV, else MeV.  Returned coefficients
      will be in the same energy units as DetEffDataPoint::energy.
   
      return value is chi2/dof.
   
      Throws exception on error with a kinda explanatory message.
   */
  double performEfficiencyFit( const std::vector<DetEffDataPoint> data,
                             const int fcnOrder,
                             std::vector<float> &result,
                             std::vector<float> &uncerts );
  
}//namespace MakeDrfFit

#endif  //MakeDrfFit_h
