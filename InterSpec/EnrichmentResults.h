#ifndef EnrichmentResults_h
#define EnrichmentResults_h

#include "InterSpec_config.h"

#include <ctime>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace Farm
{

struct EnrichmentResults
{
  std::string analysis_program;  // "FRAM" or "RelActCalcAuto"
  std::string program_version;
  std::time_t analysis_time = 0;

  struct NuclideResult
  {
    std::string nuclide;
    double mass_fraction = 0.0;
    double mass_fraction_rsd_percent = -1.0;      // -1 if not available
    double mass_fraction_systematic_percent = -1.0;
    double total_power_watts_per_kg = -1.0;       // watts/kg, -1 if not applicable
  };//struct NuclideResult

  std::vector<NuclideResult> nuclide_results;

  struct MassRatio
  {
    std::string numerator_nuclide;
    std::string denominator_nuclide;
    double mass_ratio = 0.0;
    double mass_ratio_rsd = -1.0;
    double mass_ratio_systematic = -1.0;
  };//struct MassRatio

  std::vector<MassRatio> mass_ratios;

  double days_since_separation = -1.0;      // -1 for not applicable
  double days_since_separation_uncert = -1.0;

  std::vector<std::string> warnings;

  nlohmann::json toJson() const;
  static EnrichmentResults fromJson( const std::string &json_str );
};//struct EnrichmentResults

} // namespace Farm

#endif // EnrichmentResults_h
