#ifndef EnergyCal_h
#define EnergyCal_h
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
#include <utility>

// Forward declarations.
namespace SpecUtils{ enum class EnergyCalType : int; }


namespace EnergyCal
{

struct RecalPeakInfo
{
  double peakMean;
  double peakMeanUncert;
  double peakMeanBinNumber;
  double photopeakEnergy;
};//struct RecalPeakInfo



/** Fits polynomial energy equation based on peaks with assigned
 nuclide-photopeaks.
 
 @param peakinfos Relevant information collected from peaks.
 @param fitfor Whether to fit for each coefficient.  This vector must be sized
        for the number of coefficients in the calabration equation (e.g.,
        calibration order).  A true value for an index indicates fit for the
        coefficient.  If any value is false (i.e. dont fit for that parameter)
        then the 'coefs' parameter must be exactly the same size as this vector
        and the value at the corresponding index will be used for that parameter.
 @param dev_pairs The non-linear deviation pairs.
 @param coefs Provides the fit coeficients for the energy calibration. Also,
        if any parameter is not being fit for than this vector also provides
        its (fixed) value, and this vector must be same size as the 'fitfor'
        vector.
 @param coefs_uncert The uncertainties on the fit parameters.
 @returns Chi2 of the found solution.
 
 Throws exception on error.
 */
double fit_energy_cal_poly( const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                            const std::vector<bool> fitfor,
                            const std::vector<std::pair<float,float>> &dev_pairs,
                            std::vector<float> &coefs,
                            std::vector<float> &coefs_uncert );


/** Uses Minuit2 to fit for calibration coefficients.
 
 
 Throws exception on error.
 
 @param peakinfos Relevant information collected from peaks.
 @param fitfor Whether to fit for each coefficient.  This vector must be sized
        for the number of coefficients in the calabration equation (e.g.,
        calibration order).  A true value for an index indicates fit for the
        coefficient.  If any value is false (i.e. dont fit for that parameter)
        then the 'coefs' parameter must be exactly the same size as this vector
        and the value at the corresponding index will be used for that parameter.
 @param nchannels The number of gamma channels in the spectrum the calibration is being evaluated
        for.  Must be at least 7.
 @param eqnType The type of energy calibration being fit for; must be either polynomial or FRF.
 @param startingCoefs The starting value of energy coefficients - must be same size as 'fitfor'.
 @param dev_pairs The non-linear deviation pairs.
 @param coefs Provides the fit coeficients for the energy calibration.
 @param coefs_uncert The uncertainties on the fit parameters.
 @param warning_msg If warnings were encountered (like it isnt a great fit), they will be placed
        here.  Will be empty if all looks good.
 @returns Chi2 of the found solution.
 
 \deprecated Please use #fit_energy_cal_poly.  This function is around while #fit_energy_cal_poly
             continues to be tested.
 */
double fit_energy_cal_iterative( const std::vector<EnergyCal::RecalPeakInfo> &peakInfo,
                                 const size_t nchannels,
                                 const SpecUtils::EnergyCalType eqnType,
                                 const std::vector<bool> fitfor,
                                 std::vector<float> &startingCoefs,
                                 const std::vector<std::pair<float,float>> &devpair,
                                 std::vector<float> &coefs,
                                 std::vector<float> &coefs_uncert,
                                 std::string &warning_msg );

}//namespace EnergyCal

#endif //EnergyCal_h
