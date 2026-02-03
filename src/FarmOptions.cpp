#include "InterSpec/FarmOptions.h"

#include <iostream>
#include <stdexcept>

#include <nlohmann/json.hpp>

using namespace std;

namespace Farm
{

const int FarmOptions::sm_version;

std::string FarmOptions::toJson() const
{
  nlohmann::json j;
  j["version"] = sm_version;
  j["enable_farm_analysis"] = enable_farm_analysis;
  j["enable_gadras_rid"] = enable_gadras_rid;
  j["gadras_exe_path"] = gadras_exe_path;
  j["synthesize_background_if_missing"] = synthesize_background_if_missing;
  j["enable_relact_isotopics"] = enable_relact_isotopics;
  j["enable_fram_isotopics"] = enable_fram_isotopics;
  j["fram_exe_path"] = fram_exe_path;
  j["fram_output_path"] = fram_output_path;
  j["write_fertilized_n42"] = write_fertilized_n42;

  return j.dump();
}//toJson()


FarmOptions FarmOptions::fromJson( const std::string &json_str )
{
  FarmOptions opts;

  if( json_str.empty() )
    return opts;

  try
  {
    const nlohmann::json j = nlohmann::json::parse( json_str );

    // Version check for future compatibility
    const int version = j.value( "version", 0 );
    (void)version; // Currently unused, but available for migrations

    opts.enable_farm_analysis = j.value( "enable_farm_analysis", false );
    opts.enable_gadras_rid = j.value( "enable_gadras_rid", false );
    opts.gadras_exe_path = j.value( "gadras_exe_path", string() );
    opts.synthesize_background_if_missing = j.value( "synthesize_background_if_missing", true );
    opts.enable_relact_isotopics = j.value( "enable_relact_isotopics", false );
    opts.enable_fram_isotopics = j.value( "enable_fram_isotopics", false );
    opts.fram_exe_path = j.value( "fram_exe_path", string() );
    opts.fram_output_path = j.value( "fram_output_path", string() );
    opts.write_fertilized_n42 = j.value( "write_fertilized_n42", false );
  }
  catch( const std::exception &e )
  {
    // Return default options on parse error
    std::cerr << "Failed to parse FARM options" << std::endl;
  }

  return opts;
}//fromJson(...)

} // namespace Farm
