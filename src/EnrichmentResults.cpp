#include "InterSpec/EnrichmentResults.h"

#include <nlohmann/json.hpp>

using namespace std;

namespace Farm
{

nlohmann::json EnrichmentResults::toJson() const
{
  nlohmann::json j;
  j["analysis_program"] = analysis_program;
  j["program_version"] = program_version;
  j["analysis_time"] = analysis_time;

  nlohmann::json nucs = nlohmann::json::array();
  for( const NuclideResult &nr : nuclide_results )
  {
    nlohmann::json nuc_obj;
    nuc_obj["nuclide"] = nr.nuclide;
    nuc_obj["mass_fraction"] = nr.mass_fraction;
    nuc_obj["mass_fraction_rsd_percent"] = nr.mass_fraction_rsd_percent;
    nuc_obj["mass_fraction_systematic_percent"] = nr.mass_fraction_systematic_percent;
    nuc_obj["total_power_watts_per_kg"] = nr.total_power_watts_per_kg;
    nucs.push_back( nuc_obj );
  }
  j["nuclide_results"] = nucs;

  nlohmann::json ratios = nlohmann::json::array();
  for( const MassRatio &mr : mass_ratios )
  {
    nlohmann::json ratio_obj;
    ratio_obj["numerator_nuclide"] = mr.numerator_nuclide;
    ratio_obj["denominator_nuclide"] = mr.denominator_nuclide;
    ratio_obj["mass_ratio"] = mr.mass_ratio;
    ratio_obj["mass_ratio_rsd"] = mr.mass_ratio_rsd;
    ratio_obj["mass_ratio_systematic"] = mr.mass_ratio_systematic;
    ratios.push_back( ratio_obj );
  }
  j["mass_ratios"] = ratios;

  j["days_since_separation"] = days_since_separation;
  j["days_since_separation_uncert"] = days_since_separation_uncert;
  j["warnings"] = warnings;

  return j;
}//toJson()


EnrichmentResults EnrichmentResults::fromJson( const std::string &json_str )
{
  EnrichmentResults result;

  if( json_str.empty() )
    return result;

  try
  {
    const nlohmann::json j = nlohmann::json::parse( json_str );

    result.analysis_program = j.value( "analysis_program", string() );
    result.program_version = j.value( "program_version", string() );
    result.analysis_time = j.value( "analysis_time", std::time_t(0) );

    if( j.contains("nuclide_results") && j["nuclide_results"].is_array() )
    {
      for( const auto &nuc_obj : j["nuclide_results"] )
      {
        NuclideResult nr;
        nr.nuclide = nuc_obj.value( "nuclide", string() );
        nr.mass_fraction = nuc_obj.value( "mass_fraction", 0.0 );
        nr.mass_fraction_rsd_percent = nuc_obj.value( "mass_fraction_rsd_percent", -1.0 );
        nr.mass_fraction_systematic_percent = nuc_obj.value( "mass_fraction_systematic_percent", -1.0 );
        nr.total_power_watts_per_kg = nuc_obj.value( "total_power_watts_per_kg", -1.0 );
        result.nuclide_results.push_back( nr );
      }
    }

    if( j.contains("mass_ratios") && j["mass_ratios"].is_array() )
    {
      for( const auto &ratio_obj : j["mass_ratios"] )
      {
        MassRatio mr;
        mr.numerator_nuclide = ratio_obj.value( "numerator_nuclide", string() );
        mr.denominator_nuclide = ratio_obj.value( "denominator_nuclide", string() );
        mr.mass_ratio = ratio_obj.value( "mass_ratio", 0.0 );
        mr.mass_ratio_rsd = ratio_obj.value( "mass_ratio_rsd", -1.0 );
        mr.mass_ratio_systematic = ratio_obj.value( "mass_ratio_systematic", -1.0 );
        result.mass_ratios.push_back( mr );
      }
    }

    result.days_since_separation = j.value( "days_since_separation", -1.0 );
    result.days_since_separation_uncert = j.value( "days_since_separation_uncert", -1.0 );

    if( j.contains("warnings") && j["warnings"].is_array() )
      result.warnings = j["warnings"].get<vector<string>>();
  }
  catch( const std::exception & )
  {
    // Return empty result on parse error
  }

  return result;
}//fromJson(...)

} // namespace Farm
