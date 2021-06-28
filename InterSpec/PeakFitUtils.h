#ifndef PeakFitUtils_h
#define PeakFitUtils_h
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

#include <memory>

// Forward declarations
namespace SpecUtils
{
  class Measurement;
}


namespace PeakFitUtils
{

/** Tries to guess if the passed in spectrum is from a high-resolution system (i.e., HPGe), or a lower resolution system.
 
 Primarily intended to determine defaults for peak-fitting.
 
 Currently looks at number of channels and/or average channel width, but may look at actual spectrum features for ambiguous cases
 in the future.
 

 */
bool is_high_res( const std::shared_ptr<const SpecUtils::Measurement> &meas );

}//namespace PeakFitUtils

#endif //PeakFitUtils_h
