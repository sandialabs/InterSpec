#ifndef PeakFitLM_h
#define PeakFitLM_h
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
#include <memory>
#include <utility>


// Forward declarations
class PeakDef;
class DetectorPeakResponse;
namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils


namespace PeakFitLM
{
void fit_peak_for_user_click_LM( std::vector< std::shared_ptr<const PeakDef> > &results,
                                double &chi2Dof,
                                const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                std::vector< std::shared_ptr<const PeakDef> > coFitPeaks,
                                const double mean0, const double sigma0,
                                const double area0,
                                const float roiLowerEnergy,
                                const float roiUpperEnergy );
}//namespace PeakFitLM


#endif //PeakFitLM_h
