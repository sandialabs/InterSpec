#ifndef BatchPeak_h
#define BatchPeak_h
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
#include <memory>

// Forward declarations
class PeakDef;
class SpecMeas;

namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils


namespace BatchPeak
{
  std::shared_ptr<std::deque<std::shared_ptr<const PeakDef>>>
  fit_peaks_from_exemplar( const std::shared_ptr<const SpecUtils::Measurement> &data,
                          std::vector<PeakDef> input_peaks,
                          std::vector<PeakDef> orig_peaks );
  
  // implement from fit_template_peaks
  // Steps:
  //  - Add self-atten sources testing to Act/Shield; add some tests
  //  - Add trace source tests to Act/Shield
  //  - Run all avaliable Act/Shield tests
  //  - Make structs that hold the information of ShieldingSelect's; have the ShieldingSelect
  //    put its information into these structs, and able to set ShieldingSelect from tehse structs;
  //    have the XML go to/from these structs.
  //  - Create struct to represent Nuclide fitting; e.g., if to fit activity, or age, or nuclides
  //    age, etc.
  //  - Create struct to represent user options (atten for air, etc)
  //  - ShieldingSourceDisplay::shieldingFitnessFcn a static function that accepts only non-widget
  //    input (e.g., the above)
  //  - Test nothing broke.
  
}//namespace BatchPeak

#endif //BatchPeak_h
