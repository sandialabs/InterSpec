#ifndef AnalystChecks_h
#define AnalystChecks_h
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

#include <string>
#include <vector>
#include <optional>
#include <memory>

#include "SpecUtils/SpecFile.h"

// Forward declarations
class PeakDef;
class InterSpec;

namespace AnalystChecks
{
  /** A simple "hello world" type function for testing purposes.
   
   Returns a greeting message indicating the AnalystChecks module is working.
   
   @return A string containing a hello world message.
   */
  InterSpec_API std::string hello_world();
  
  /** Options for peak detection analysis. */
  struct DetectedPeaksOptions {
    SpecUtils::SpectrumType specType;
    std::optional<std::string> userSession;
  };
  
  /** Results of peak detection analysis. */
  struct DetectedPeakStatus {
    std::string userSession;
    std::vector<std::shared_ptr<const PeakDef>> peaks;
  };
  
  /** Perform automated peak detection on the specified spectrum.
   
   This function combines user-defined peaks with automatically detected peaks,
   avoiding duplicates by checking peak overlap.
   
   @param options Specifies which spectrum to analyze and optional user session
   @param interspec Pointer to the InterSpec session containing the spectrum
   @return DetectedPeakStatus containing the found peaks and session info
   @throws std::runtime_error if spectrum is not available or other errors occur
   */
  InterSpec_API DetectedPeakStatus detected_peaks(const DetectedPeaksOptions& options, InterSpec* interspec);
  
  
  struct FitPeakOptions {
    bool addToUsersPeaks;
    double energy;
    SpecUtils::SpectrumType specType;
    std::optional<std::string> source;
    std::optional<std::string> userSession;
  };
  
  struct FitPeakStatus {
    std::string userSession;
    std::shared_ptr<const PeakDef> fitPeak;
    std::vector<std::shared_ptr<const PeakDef>> peaksInRoi;
  };
  
  InterSpec_API FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec );
  
  struct GetUserPeakOptions {
    SpecUtils::SpectrumType specType;
    std::optional<std::string> userSession;
  };
  
  struct GetUserPeakStatus {
    std::string userSession;
    std::vector<std::shared_ptr<const PeakDef>> peaks;
  };
  
  InterSpec_API GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec );
  
  InterSpec_API std::vector<float> get_characteristic_gammas( const std::string &nuclide );
  
} // namespace AnalystChecks

#endif // AnalystChecks_h 
