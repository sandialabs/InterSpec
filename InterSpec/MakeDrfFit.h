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

namespace MakeDrfFit
{
  /* If result passed into the function is of the proper size, then it will be
    used as the starting parameters for the fit, otherwise some default values
    will be chosen to start the fit with.
    The provided Measurement should be same one peaks where fit for.
   
    @param sqrtEqnOrder Only used if fnctnlForm==DetectorPeakResponse::kSqrtPolynomial
   
    @returns chi2 (not divided by dof).
   
    Throws exception on error with a kinda explanatory message.
  */
  double performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                             const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                             const bool highResolution,
                             const int sqrtEqnOrder,
                             std::vector<float> &result,
                             std::vector<float> &uncerts );
  
  /** Fits the FWHM as a sqrt( Sum_i {A_i *pow(x,i)} )
   
   @param peaks The peaks to use for the fit
   @param num_fit_coefficients The number of coefficients to fit _not_ the fit order (which would be num energy
          dependent coefficients).  Must be larger than zero
   @param include_inv_term If true, fits sqrt(A + B*x + C/x) instead of sqrt(A + B*x + C*x*x + D*x*x + ... )
          Right now, if true, assumes will fit for three coefficients.
   @param[out] coeffs The resulting fit coefficients.
   @param[out] coeffs The resulting fit coefficient uncertainties.
   @returns the chi2 of the fit
   
   Throws exception if fitting for more parameters than there are peaks, or if number of coefficients isnt at least 1, or the fit fails.
   */
  double fit_sqrt_poly_fwhm_lls( const std::deque< std::shared_ptr<const PeakDef> > &peaks,
                                       const int num_fit_coefficients,
                                       const bool include_inv_term,
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
