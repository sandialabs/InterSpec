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
    Throws exception on error with a kinda explanitory message.
  */
  void performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                             const size_t num_gamma_channels,
                             const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                             std::vector<float> &result,
                             std::vector<float> &uncerts );
  
  
  double peak_width_chi2( double predicted_sigma, const PeakDef &peak );
}//namespace MakeDrfFit

#endif  //MakeDrfFit_h
