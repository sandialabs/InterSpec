#ifndef FarmOptions_h
#define FarmOptions_h

#include "InterSpec_config.h"

#include <string>
#include <vector>

namespace Farm
{

struct FarmOptions
{
  // Master enable - if false, no FARM analysis is performed
  bool enable_farm_analysis = false;

  // GADRAS Remote RID options
  bool enable_gadras_rid = false;
  std::string gadras_exe_path;  // Path to "Gadras Full Spectrum Isotope ID" EXE
  bool synthesize_background_if_missing = true;

  // Isotopics via RelActCalcAuto
  bool enable_relact_isotopics = false;
  // Uses "data/rel_act/HPGe U (120-1001 keV).xml" or "data/rel_act/HPGe Pu (120-780 keV).xml"

  // Isotopics via FRAM EXE
  bool enable_fram_isotopics = false;
  std::string fram_exe_path;
  std::string fram_output_path;  // Directory FRAM writes output to

  // Output options
  bool write_fertilized_n42 = false;

  bool operator==( const FarmOptions &rhs ) const
  {
    return enable_farm_analysis == rhs.enable_farm_analysis
        && enable_gadras_rid == rhs.enable_gadras_rid
        && gadras_exe_path == rhs.gadras_exe_path
        && synthesize_background_if_missing == rhs.synthesize_background_if_missing
        && enable_relact_isotopics == rhs.enable_relact_isotopics
        && enable_fram_isotopics == rhs.enable_fram_isotopics
        && fram_exe_path == rhs.fram_exe_path
        && fram_output_path == rhs.fram_output_path
        && write_fertilized_n42 == rhs.write_fertilized_n42;
  }

  bool operator!=( const FarmOptions &rhs ) const { return !(*this == rhs); }

  // Serialization
  std::string toJson() const;
  static FarmOptions fromJson( const std::string &json_str );

  static const int sm_version = 1;
};

} // namespace Farm

#endif // FarmOptions_h
