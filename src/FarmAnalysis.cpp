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

#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <boost/process.hpp>
#if( defined(_WIN32) )
#include <boost/process/windows.hpp>
#endif
#include <boost/filesystem.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/FarmAnalysis.h"
#include "InterSpec/EnrichmentResults.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/SpecFileQueryDbCache.h"

#include <nlohmann/json.hpp>


using namespace std;

namespace Farm
{

std::vector<std::shared_ptr<const PeakDef>> perform_peak_search(
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::shared_ptr<const DetectorPeakResponse> &drf,
    const bool is_hpge,
    SpecFileInfoToQuery &info )
{
  info.farm_peaks_json.clear();

  if( !foreground || foreground->num_gamma_channels() < 64 )
    return std::vector<std::shared_ptr<const PeakDef>>{};

  try
  {
    std::shared_ptr<std::deque<std::shared_ptr<const PeakDef>>> existing_peaks
        = std::make_shared<std::deque<std::shared_ptr<const PeakDef>>>();

    // Use single-threaded for batch processing stability
    std::vector<std::shared_ptr<const PeakDef>> peaks =
        ExperimentalAutomatedPeakSearch::search_for_peaks(
            foreground,
            drf,
            existing_peaks,
            true,   // singleThreaded - important for batch processing
            is_hpge );

    // Serialize peaks to JSON array
    nlohmann::json peaks_array = nlohmann::json::array();
    for( const std::shared_ptr<const PeakDef> &peak : peaks )
    {
      if( !peak )
        continue;

      nlohmann::json peak_obj;
      peak_obj["mean"] = peak->mean();
      peak_obj["fwhm"] = peak->fwhm();
      peak_obj["amplitude"] = peak->amplitude();
      peak_obj["area"] = peak->peakArea();
      peak_obj["area_uncert"] = peak->peakAreaUncert();
      peak_obj["chi2dof"] = peak->chi2dof();
      peak_obj["roi_lower"] = peak->lowerX();
      peak_obj["roi_upper"] = peak->upperX();

      if( peak->parentNuclide() )
        peak_obj["nuclide"] = peak->parentNuclide()->symbol;

      peaks_array.push_back( peak_obj );
    }//for( peak )

    info.farm_peaks_json = peaks_array.dump();

#if( SpecUtils_ENABLE_D3_CHART )
    // Full ROI-grouped JSON (same format the D3 chart uses); includes continuum
    // coefficients, skew params, colors, etc.
    info.farm_peaks_full_json = PeakDef::peak_json( peaks, foreground,
                                                    Wt::WColor(0,51,255), 255 );
#endif

    return peaks;
  }
  catch( const std::exception &e )
  {
    // Store error in JSON for debugging
    nlohmann::json error_obj;
    error_obj["error"] = e.what();
    info.farm_peaks_json = error_obj.dump();
  }

  return std::vector<std::shared_ptr<const PeakDef>>{};
}//perform_peak_search(...)


  /*
std::string detector_type_to_gadras_drf( const SpecUtils::DetectorType type )
{
   // Avaiable DRF names are:
   //[ "auto", "1x1/BGO Side", "1x1/CsI Side", "1x1/LaCl3", "1x1/NaI Front", "1x1/NaI Side", "3x3/NaI InCorner",
   //"3x3/NaI LowScat", "3x3/NaI MidScat", "3x3/NaI OnGround",
   //"Atomex-AT6102",
   //"D3S", "Detective", "Detective-EX", "Detective-EX100", "Detective-EX200", "Detective-Micro", "Detective-Micro/Variant-LowEfficiency", "Detective-X",
   //"Falcon 5000", "FieldSpec", "Fulcrum40h", "GR130", "GR135", "GR135Plus",
   //"IdentiFINDER-LaBr3", "IdentiFINDER-N","IdentiFINDER-NG", "IdentiFINDER-NGH", "IdentiFINDER-R300", "IdentiFINDER-R425", "IdentiFINDER-R500-NaI",
   //"InSpector 1000 LaBr3", "InSpector 1000 NaI",
   //"Interceptor", "Kromek-D5", "Kromek-GR1-CZT", "MKC-A03", "Mirion PDS-100", "NaI 2x4x16", "Polimaster PM1704-GN",
   //"RIIDEyeX-GN1", "RadEagle", "RadEye", "RadPack", "RadSeeker-NaI", "Radiacode-102", "Radseeker-LaBr3",
   //"Raider", "Ranger", "Raysid", "SAM-935", "SAM-945", "SAM-950GN-N30", "SAM-Eagle-LaBr3", "SAM-Eagle-NaI-3x3",
   //"SpiR-ID/LaBr3", "SpiR-ID/NaI",
   //"Thermo ARIS Portal", "Verifinder"]
   
  switch( type )
  {
    case SpecUtils::DetectorType::DetectiveEx:
      return "Detective-EX";
    case SpecUtils::DetectorType::DetectiveEx100:
      return "Detective-EX100";
    case SpecUtils::DetectorType::DetectiveEx200:
      return "Detective-EX200";
    case SpecUtils::DetectorType::IdentiFinder:
    case SpecUtils::DetectorType::IdentiFinderNG:
      return "identiFINDER-NGH";
    case SpecUtils::DetectorType::IdentiFinderLaBr3:
      return "identiFINDER-LaBr";
    case SpecUtils::DetectorType::Sam940:
      return "SAM-940";
    case SpecUtils::DetectorType::Sam945:
      return "SAM-945";
    case SpecUtils::DetectorType::Rsi701:
    case SpecUtils::DetectorType::Rsi705:
      return "NaI 2x4x16";
  ...
    default:
      return "";
  }
}//detector_type_to_gadras_drf(...)
*/
  

std::string run_gadras_full_spectrum_id_analysis(
    const std::string &exe_path,
    std::shared_ptr<SpecUtils::SpecFile> spec_file,
    const SpecUtils::DetectorType detector_type,
    const bool synthesize_background )
{
  if( exe_path.empty() || !spec_file || (detector_type == SpecUtils::DetectorType::Unknown) )
    return "";

  // Build command-line arguments
  // Pattern from RemoteRid.cpp:1128-1237
  std::vector<std::string> arguments;
  arguments.push_back( "--mode=command-line" );
  arguments.push_back( "--out-format=json" );
  //arguments.push_back( "--drf" );
  //const std::string drf_name = detector_type_to_gadras_drf( detector_type );
  //arguments.push_back( drf_name );

  if( spec_file->num_measurements() < 2 )
  {
    if( synthesize_background )
      arguments.push_back( "--synthesize-background=1" );
    else
      return "";
  }

  // Create temporary N42 file
  const std::string tmpfilename = SpecUtils::temp_file_name(
      "farm_gadras_", SpecUtils::temp_dir() );

  try
  {
    // Write spectrum to temp file
    {
      std::ofstream tmpfile( tmpfilename.c_str(), std::ios::out | std::ios::binary );
      if( !tmpfile.is_open() )
        throw std::runtime_error( "Failed to create temp file for GADRAS" );
      spec_file->write_2012_N42( tmpfile );
    }

    arguments.push_back( tmpfilename );

    // Execute GADRAS using boost::process
    // Pattern from RemoteRid.cpp:106-190 run_external_command()
    namespace bp = boost::process;

    const boost::filesystem::path exe_parent = boost::filesystem::path(exe_path).parent_path();
    bp::ipstream proc_stdout, proc_stderr;

#ifdef _WIN32
    bp::child c( exe_path, bp::args(arguments), bp::start_dir(exe_parent),
                 bp::std_out > proc_stdout, bp::std_err > proc_stderr,
                 bp::windows::create_no_window );
#else
    bp::child c( exe_path, bp::args(arguments), bp::start_dir(exe_parent),
                 bp::std_out > proc_stdout, bp::std_err > proc_stderr );
#endif

    c.wait();

    std::string output( std::istreambuf_iterator<char>(proc_stdout), {} );
    std::string error( std::istreambuf_iterator<char>(proc_stderr), {} );

    const int result = c.exit_code();

    // Clean up temp file
    SpecUtils::remove_file( tmpfilename );

    if( (result != EXIT_SUCCESS) && (!error.empty() || output.empty()) )
    {
      nlohmann::json err_json;
      err_json["error"] = error.empty() ? "GADRAS returned non-zero exit code" : error;
      err_json["exit_code"] = result;
      return err_json.dump();
    }

    SpecUtils::trim( output );
    return output;  // JSON from GADRAS
  }catch( const std::exception &e )
  {
    SpecUtils::remove_file( tmpfilename );

    nlohmann::json err_json;
    err_json["error"] = e.what();
    return err_json.dump();
  }
}//run_gadras_full_spectrum_id_analysis(...)


bool should_do_uranium_isotopics( const std::string &peaks_json )
{
  if( peaks_json.empty() )
    return false;

  try
  {
    const nlohmann::json peaks = nlohmann::json::parse( peaks_json );
    if( !peaks.is_array() )
      return false;

    bool has_185_peak = false;
    bool has_205_or_144_peak = false;

    for( const auto &peak : peaks )
    {
      const double mean = peak.value( "mean", 0.0 );
      const double fwhm = peak.value( "fwhm", 2.0 );  // Default ~2 keV
      const double tolerance = 0.75 * fwhm;

      // Check for 185 keV (U-235 signature)
      if( std::fabs( mean - 185.7 ) < tolerance )
        has_185_peak = true;

      // Check for 205.3 keV or 143.8 keV (U-235 signatures)
      if( std::fabs( mean - 205.3 ) < tolerance || std::fabs( mean - 143.8 ) < tolerance )
        has_205_or_144_peak = true;
    }

    return has_185_peak && has_205_or_144_peak;
  }
  catch( ... )
  {
    return false;
  }
}//should_do_uranium_isotopics(...)


bool should_do_plutonium_isotopics( const std::string &peaks_json )
{
  if( peaks_json.empty() )
    return false;

  try
  {
    const nlohmann::json peaks = nlohmann::json::parse( peaks_json );
    if( !peaks.is_array() )
      return false;

    bool has_375_peak = false;
    bool has_secondary_peak = false;

    for( const auto &peak : peaks )
    {
      const double mean = peak.value( "mean", 0.0 );
      const double fwhm = peak.value( "fwhm", 2.0 );
      const double tolerance = 0.75 * fwhm;

      // Check for 375 keV (Pu-239 signature)
      if( std::fabs( mean - 375.1 ) < tolerance )
        has_375_peak = true;

      // Check for secondary Pu peaks: 129.3, 413.7, or 662 keV
      if( std::fabs( mean - 129.3 ) < tolerance ||
          std::fabs( mean - 413.7 ) < tolerance ||
          std::fabs( mean - 662.0 ) < tolerance )
        has_secondary_peak = true;
    }

    return has_375_peak && has_secondary_peak;
  }
  catch( ... )
  {
    return false;
  }
}//should_do_plutonium_isotopics(...)


EnrichmentResults run_relact_isotopics(
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::shared_ptr<const SpecUtils::Measurement> &background,
    const std::vector<std::shared_ptr<const PeakDef>> &peaks,
    const std::string &options_xml_path,
    const std::shared_ptr<const DetectorPeakResponse> &drf,
    MaterialDB *materialDB )
{
  EnrichmentResults result;
  result.analysis_program = "RelActCalcAuto";
  result.analysis_time = std::time( nullptr );

  if( !foreground )
  {
    result.warnings.push_back( "No foreground spectrum provided" );
    return result;
  }

  try
  {
    std::ifstream xml_file( options_xml_path );
    if( !xml_file.is_open() )
    {
      result.warnings.push_back( "Could not open isotopics config: " + options_xml_path );
      return result;
    }

    const std::string xml_content( (std::istreambuf_iterator<char>(xml_file)),
                                    std::istreambuf_iterator<char>() );

    rapidxml::xml_document<char> doc;
    std::vector<char> xml_copy( xml_content.begin(), xml_content.end() );
    xml_copy.push_back( '\0' );
    doc.parse<rapidxml::parse_trim_whitespace>( xml_copy.data() );

    rapidxml::xml_node<char> *base_node = doc.first_node( "RelActCalcAuto" );
    if( !base_node )
      base_node = doc.first_node();

    rapidxml::xml_node<char> *opts_node = base_node ? base_node->first_node("Options") : nullptr;
    if( !opts_node )
    {
      result.warnings.push_back( "Invalid isotopics XML config format: " + options_xml_path );
      return result;
    }

    RelActCalcAuto::Options options;
    options.fromXml( opts_node, materialDB );

    const RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
        options, foreground, background, drf, peaks, nullptr );

    if( solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      result.warnings.push_back( "RelActCalcAuto failed: " + solution.m_error_message );
      return result;
    }

    // Extract mass enrichment fractions from solution
    for( size_t curve_idx = 0; curve_idx < solution.m_rel_activities.size(); ++curve_idx )
    {
      for( const RelActCalcAuto::NuclideRelAct &nuc_act : solution.m_rel_activities[curve_idx] )
      {
        const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( nuc_act.source );
        if( !nuc )
          continue;

        const std::pair<double,std::optional<double>> mass_frac_pair =
            solution.mass_enrichment_fraction( nuc, curve_idx );

        EnrichmentResults::NuclideResult nr;
        nr.nuclide = nuc->symbol;
        nr.mass_fraction = mass_frac_pair.first;
        if( mass_frac_pair.second.has_value() && mass_frac_pair.first > 0.0 )
          nr.mass_fraction_rsd_percent = 100.0 * (*mass_frac_pair.second) / mass_frac_pair.first;

        result.nuclide_results.push_back( nr );
      }
    }

    // Pu-242 by correlation correction (only populated for Pu options)
    for( size_t i = 0; i < solution.m_corrected_pu.size(); ++i )
    {
      if( solution.m_corrected_pu[i] )
      {
        const RelActCalc::Pu242ByCorrelationOutput<double> &pu_corr = *solution.m_corrected_pu[i];
        EnrichmentResults::NuclideResult nr;
        nr.nuclide = "Pu242";
        nr.mass_fraction = pu_corr.pu242_mass_frac;
        result.nuclide_results.push_back( nr );
      }
    }

    result.warnings = solution.m_warnings;
    result.program_version = "RelActCalcAuto";
  }
  catch( const std::exception &e )
  {
    result.warnings.push_back( std::string("Exception: ") + e.what() );
  }

  return result;
}//run_relact_isotopics(...)


EnrichmentResults run_fram_isotopics(
    const std::string &fram_exe_path,
    const std::string &fram_output_path,
    std::shared_ptr<const SpecUtils::SpecFile> foreground,
    std::shared_ptr<const SpecUtils::SpecFile> background,
    const bool is_uranium,
    const bool is_plutonium )
{
  EnrichmentResults result;
  result.analysis_program = "FRAM";
  result.analysis_time = std::time( nullptr );

  // TODO: Implement FRAM executable calling
  // The user will fill in the details of the call and output parsing
  // Placeholder structure:

  if( fram_exe_path.empty() )
  {
    result.warnings.push_back( "FRAM executable path not specified" );
    return result;
  }

  if( !SpecUtils::is_file( fram_exe_path ) )
  {
    result.warnings.push_back( "FRAM executable not found: " + fram_exe_path );
    return result;
  }

  assert( foreground && (foreground->num_measurements() == 1) );
  if( !foreground || (foreground->num_measurements() != 1) )
  {
    result.warnings.push_back( "FRAM computation input not expected number of records." );
    return result;
  }

  assert( !background || (background->num_measurements() == 1) );
  if( background && (background->num_measurements() != 1) );
  {
    result.warnings.push_back( "FRAM computation input background not expected number of records." );
    return result;
  }


  // Write a temp N42 for FRAM input
  const SpecUtils::SaveSpectrumAsType output_format = SpecUtils::SaveSpectrumAsType::SpcBinaryInt;
  const std::string fram_fore_tmp = SpecUtils::temp_file_name( "farm_fram_foreground_", SpecUtils::temp_dir() )
                                    + string(".") + SpecUtils::suggestedNameEnding( output_format );

  try
  {
    foreground->write_to_file( fram_fore_tmp, output_format );
  }catch( std::exception &e )
  {
    result.warnings.push_back( "FRAM computation error: failed to write temporary foreground file: " + string(e.what()) );
    return result;
  }


  std::string fram_back_tmp;
  if( background )
  {
    fram_back_tmp = SpecUtils::temp_file_name( "farm_fram_foreground_", SpecUtils::temp_dir() )
                    + string(".") + SpecUtils::suggestedNameEnding( output_format );
    try
    {
      background->write_to_file( fram_back_tmp, output_format );
    }catch( std::exception &e )
    {
      result.warnings.push_back( "FRAM computation error: failed to write temporary background file: " + string(e.what()) );
      return result;
    }
  }//if( background )

  assert( foreground->measurements().size() == 1 );
  std::shared_ptr<const SpecUtils::Measurement> fore_spec = foreground->measurements()[0];
  assert( fore_spec );

  double energy_offset = 0.0, energy_gain = -9999.9;
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = fore_spec->energy_calibration();
  if( !cal || !cal->valid() )
  {
    result.warnings.push_back( "FRAM computation error: energy calibration is invalid" );
    return result;
  }

  switch( cal->type() )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      assert( cal->coefficients().size() >= 2 );
      energy_offset = cal->coefficients()[0];
      energy_gain = cal->coefficients()[1];
      break;

    case SpecUtils::EnergyCalType::FullRangeFraction:
    {
      const std::vector<float> poly = SpecUtils::fullrangefraction_coef_to_polynomial(cal->coefficients(), cal->num_channels());
      assert( poly.size() >= 2 );
      energy_offset = poly[0];
      energy_gain = poly[1];
      break;
    }

    case SpecUtils::EnergyCalType::LowerChannelEdge:
      energy_offset = cal->lower_energy();
      energy_gain = (cal->upper_energy() - cal->lower_energy()) / cal->num_channels();
      break;

    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
      break;
  }//switch( cal->type() )

  // Path to foreground file: fram_fore_tmp
  // Path to background file: fram_back_tmp (empty if no background)
  // Energy offset: energy_offset
  // Energy gain: energy_gain
  // Type of isotopics: is_uranium, is_plutonium (both may be true)

  result.warnings.push_back( "FRAM integration not yet implemented - user to fill in details" );

  // To execute the FRAM exe, a rough sketch is:
  /*
   std::vector<std::string> arguments;
   arguments.push_back( "--offset=" + std::to_string(energy_gain) );
   arguments.push_back( "--gain=" + std::to_string(energy_gain) );
   arguments.push_back( "--foreground='" + fram_fore_tmp + "'" );
   if( !fram_back_tmp.empty() )
    arguments.push_back( "--background='" + fram_back_tmp + "'" );
   //...

   namespace bp = boost::process;

   const boost::filesystem::path exe_parent = boost::filesystem::path(fram_exe_path).parent_path();
   bp::ipstream proc_stdout, proc_stderr;

#ifdef _WIN32
   bp::child c( fram_exe_path, bp::args(arguments), bp::start_dir(exe_parent),
                bp::std_out > proc_stdout, bp::std_err > proc_stderr,
                bp::windows::create_no_window );
#else
   bp::child c( fram_exe_path, bp::args(arguments), bp::start_dir(exe_parent),
                bp::std_out > proc_stdout, bp::std_err > proc_stderr );
#endif

   c.wait();

   std::string output( std::istreambuf_iterator<char>(proc_stdout), {} );
   std::string error( std::istreambuf_iterator<char>(proc_stderr), {} );

   const int result_code = c.exit_code();

   if( (result_code != EXIT_SUCCESS) && (!error.empty() || output.empty()) )
   {
     nlohmann::json err_json;
     err_json["error"] = error.empty() ? "FRAM returned non-zero exit code" : error;
     err_json["exit_code"] = result_code;
     return err_json.dump();
   }

   // TODO: Populate result.nuclide_results from parsed output
   result = ...
   */


  //Note: NuclideResult::nuclide should be in format "U235", "Pu239", etc, for consistent search results

  // Clean up temp file
  assert( SpecUtils::is_file( fram_fore_tmp ) );
  if( !fram_back_tmp.empty() )
    SpecUtils::remove_file( fram_back_tmp );

  return result;
}//run_fram_isotopics(...)


void write_fertilized_n42(
    const std::string &original_file_path,
    const SpecUtils::SpecFile &meas,
    const SpecFileInfoToQuery &info )
{
  // Build JSON remark containing all FARM results
  nlohmann::json farm_remark;
  farm_remark["farm_version"] = 1;
  farm_remark["original_file"] = original_file_path;
  farm_remark["analysis_time"] = std::time( nullptr );

  // Add simplified peaks summary if present
  if( !info.farm_peaks_json.empty() && info.farm_peaks_json != "null" )
  {
    try
    {
      farm_remark["peaks"] = nlohmann::json::parse( info.farm_peaks_json );
    }catch( ... )
    {
      assert( 0 );
    }
  }

  // Add full ROI-grouped peak JSON (continuum coefficients, skew, etc.)
  if( !info.farm_peaks_full_json.empty() && info.farm_peaks_full_json != "null" )
  {
    try
    {
      farm_remark["peaks_full"] = nlohmann::json::parse( info.farm_peaks_full_json );
    }catch( ... )
    {
      assert( 0 );
    }
  }

  // Add GADRAS results if present
  if( !info.gadras_rid_json.empty() && info.gadras_rid_json != "null" )
  {
    try
    {
      farm_remark["gadras_rid"] = nlohmann::json::parse( info.gadras_rid_json );
    }catch( ... )
    {
      assert( 0 );
    }
  }

  // Add isotopics results if present
  if( !info.isotopics_result_json.empty() && info.isotopics_result_json != "null" )
  {
    try
    {
      farm_remark["computed_isotopics"] = nlohmann::json::parse( info.isotopics_result_json );
    }
    catch( ... ) {}
  }

  // Add statistical moments
  nlohmann::json stats;
  stats["counts_mean"] = info.spectrum_mean;
  stats["counts_variance"] = info.spectrum_variance;
  stats["counts_standard_deviation"] = sqrt(info.spectrum_variance);
  stats["counts_skewness"] = info.spectrum_skewness;
  stats["counts_kurtosis"] = info.spectrum_kurtosis;
  farm_remark["statistics"] = stats;

  // Foreground summary fields
  if( info.farm_min_channel_with_data >= 0 )
    farm_remark["min_channel_with_data"] = info.farm_min_channel_with_data;
  if( info.farm_max_channel_with_data >= 0 )
    farm_remark["max_channel_with_data"] = info.farm_max_channel_with_data;
  farm_remark["total_gamma_counts"]  = info.farm_foreground_total_gamma_counts;
  farm_remark["num_gamma_channels"]  = info.farm_foreground_num_gamma_channels;
  farm_remark["min_gamma_count"]     = info.farm_foreground_min_gamma_count;
  farm_remark["max_gamma_count"]     = info.farm_foreground_max_gamma_count;
  farm_remark["has_neutrons"]        = info.farm_foreground_has_neutrons;
  if( info.farm_foreground_has_neutrons )
  {
    farm_remark["neutron_count"]     = info.farm_foreground_neutron_count;
    farm_remark["neutron_live_time"] = info.farm_foreground_neutron_live_time;
  }

  // Energy calibration object (omitted entirely when cal was invalid)
  if( !info.farm_energy_cal_json.empty() )
  {
    try
    {
      farm_remark["energy_calibration"] = nlohmann::json::parse( info.farm_energy_cal_json );
    }catch( ... )
    {
      assert( 0 );
    }
  }

  // Instrument RID/NuclID results from the original SpecFile
  {
    const std::shared_ptr<const SpecUtils::DetectorAnalysis> det_ana = meas.detectors_analysis();
    if( det_ana && !det_ana->is_empty() )
    {
      nlohmann::json ana_json;

      if( !det_ana->remarks_.empty() )
        ana_json["remarks"] = det_ana->remarks_;
      if( !det_ana->algorithm_name_.empty() )
        ana_json["algorithm_name"] = det_ana->algorithm_name_;
      if( !det_ana->algorithm_component_versions_.empty() )
      {
        nlohmann::json versions = nlohmann::json::array();
        for( const std::pair<std::string,std::string> &v : det_ana->algorithm_component_versions_ )
        {
          nlohmann::json entry;
          entry["component"] = v.first;
          entry["version"]   = v.second;
          versions.push_back( entry );
        }
        ana_json["algorithm_component_versions"] = versions;
      }
      if( !det_ana->algorithm_creator_.empty() )
        ana_json["algorithm_creator"] = det_ana->algorithm_creator_;
      if( !det_ana->algorithm_description_.empty() )
        ana_json["algorithm_description"] = det_ana->algorithm_description_;
      if( det_ana->analysis_start_time_ != SpecUtils::time_point_t{} )
      {
        const auto epoch_us = std::chrono::duration_cast<std::chrono::microseconds>(
            det_ana->analysis_start_time_.time_since_epoch() ).count();
        ana_json["analysis_start_time_epoch_us"] = epoch_us;
      }
      if( det_ana->analysis_computation_duration_ > 0.0f )
        ana_json["analysis_computation_duration_sec"] = det_ana->analysis_computation_duration_;
      if( !det_ana->algorithm_result_description_.empty() )
        ana_json["algorithm_result_description"] = det_ana->algorithm_result_description_;

      if( !det_ana->results_.empty() )
      {
        nlohmann::json results_arr = nlohmann::json::array();
        for( const SpecUtils::DetectorAnalysisResult &r : det_ana->results_ )
        {
          if( r.isEmpty() )
            continue;
          nlohmann::json res;
          if( !r.nuclide_.empty() )
            res["nuclide"] = r.nuclide_;
          if( !r.nuclide_type_.empty() )
            res["nuclide_type"] = r.nuclide_type_;
          if( r.activity_ >= 0.0f )
            res["activity_bq"] = r.activity_;
          if( !r.id_confidence_.empty() )
            res["id_confidence"] = r.id_confidence_;
          if( r.distance_ >= 0.0f )
            res["distance_mm"] = r.distance_;
          if( r.dose_rate_ >= 0.0f )
            res["dose_rate_uSv_per_hr"] = r.dose_rate_;
          if( r.real_time_ >= 0.0f )
            res["real_time_sec"] = r.real_time_;
          if( !r.detector_.empty() )
            res["detector"] = r.detector_;
          if( !r.remark_.empty() )
            res["remark"] = r.remark_;
          results_arr.push_back( res );
        }
        if( !results_arr.empty() )
          ana_json["results"] = results_arr;
      }

      farm_remark["instrument_rid_nucid"] = ana_json;
    }
  }

  // Create copy of SpecFile and add remark
  SpecUtils::SpecFile output_meas = meas;

  std::vector<std::string> remarks = output_meas.remarks();
  
  nlohmann::json final_remark_json = nlohmann::json::object();
  final_remark_json["InterSpec"] = farm_remark;
  remarks.push_back( "FARM_REMARK: " + final_remark_json.dump() );
  
  
  nlohmann::json orig_uuid_json = nlohmann::json::object();
  orig_uuid_json["OriginalFileUuid"] = meas.uuid();
  remarks.push_back( "FARM_REMARK: " + orig_uuid_json.dump() );
  
  output_meas.set_remarks( remarks );

  // Determine output path
  const std::string parent_dir = SpecUtils::parent_path( original_file_path );
  const std::string filename = SpecUtils::filename( original_file_path );
  const std::string output_path = SpecUtils::append_path(
      parent_dir, filename + ".farm.fertilized.n42" );

  output_meas.set_uuid( "" ); // Make sure UUID will get updated.
  
  
  // Write output file
  std::ofstream out( output_path, std::ios::binary );
  if( out.is_open() )
    output_meas.write_2012_N42( out );
}//write_fertilized_n42(...)

} // namespace Farm
