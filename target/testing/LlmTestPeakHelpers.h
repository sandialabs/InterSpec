#ifndef LlmTestPeakHelpers_h
#define LlmTestPeakHelpers_h

/* InterSpec: an application to analyze spectral gamma radiation data.

 Shared helpers for the LLM tool unit-test suites.
 */

#include <string>
#include <optional>
#include <stdexcept>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AnalystChecks.h"

class InterSpec;

namespace LlmTestHelpers
{
  /** Drive the `add_analysis_peak` tool synchronously from a unit test.

   `add_analysis_peak` (and `add_analysis_peaks_for_source`) are async tools: their asyncExecutor
   posts the result back through WServer to the GUI thread, which a unit test's blocked thread can't
   pump - so calling them through the synchronous `ToolRegistry::executeTool()` throws "is async".
   This helper drives the identical synchronous fitting core (`AnalystChecks::fit_user_peak`) and
   returns a minimal result mirroring the tool's JSON contract: `{"fitPeakEnergy": <mean>}` on
   success, or `{"error": <message>}` on failure - which is all the add_analysis_peak assertions in
   the LLM test suites inspect.  The peak is actually added to the PeakModel, so downstream
   get_peaks / edit_analysis_peak / activity-fit assertions still hold.
   */
  inline nlohmann::json add_analysis_peak_sync( const nlohmann::json &params, InterSpec *interspec )
  {
    try
    {
      AnalystChecks::FitPeakOptions options;
      options.energy = params.at("energy").get<double>();
      options.doNotAddToAnalysisPeaks = params.value("DoNotAddToAnalysisPeaks", false);

      const std::string specTypeStr = params.value("specType", std::string("Foreground"));
      if( specTypeStr.empty() || specTypeStr == "Foreground" )
        options.specType = SpecUtils::SpectrumType::Foreground;
      else if( specTypeStr == "Background" )
        options.specType = SpecUtils::SpectrumType::Background;
      else if( specTypeStr == "Secondary" )
        options.specType = SpecUtils::SpectrumType::SecondForeground;
      else
        throw std::runtime_error( "Invalid spectrum type: " + specTypeStr );

      // Mirror the tool's from_json source resolution order.
      options.source = std::nullopt;
      for( const char *key : {"source", "nuclide", "element", "xray", "x-ray", "reaction"} )
      {
        if( params.contains(key) )
        {
          options.source = params[key].get<std::string>();
          break;
        }
      }

      const AnalystChecks::FitPeakStatus status = AnalystChecks::fit_user_peak( options, interspec );
      if( !status.fitPeak )
        return nlohmann::json{ {"error", "No peak was fit"} };

      nlohmann::json result;
      result["fitPeakEnergy"] = status.fitPeak->mean();
      return result;
    }catch( const std::exception &e )
    {
      return nlohmann::json{ {"error", std::string(e.what())} };
    }
  }//add_analysis_peak_sync
}//namespace LlmTestHelpers

#endif //LlmTestPeakHelpers_h
