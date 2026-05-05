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

#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "rapidxml/rapidxml.hpp"

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"
#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/BatchRelActAuto.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActAutoReport.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;

namespace BatchRelActAuto
{

const char *to_str( const ResultCode code )
{
  switch( code )
  {
    case ResultCode::Success:                              return "Success";
    case ResultCode::NoExemplar:                           return "NoExemplar";
    case ResultCode::CouldntOpenExemplar:                  return "CouldntOpenExemplar";
    case ResultCode::ExemplarMissingRelActState:           return "ExemplarMissingRelActState";
    case ResultCode::CouldntOpenStateOverride:             return "CouldntOpenStateOverride";
    case ResultCode::CouldntOpenInputFile:                 return "CouldntOpenInputFile";
    case ResultCode::CouldntOpenBackgroundFile:            return "CouldntOpenBackgroundFile";
    case ResultCode::ForegroundSampleNumberUnderSpecified: return "ForegroundSampleNumberUnderSpecified";
    case ResultCode::BackgroundSampleNumberUnderSpecified: return "BackgroundSampleNumberUnderSpecified";
    case ResultCode::FwhmMethodNeedsDrfButNoneAvailable:   return "FwhmMethodNeedsDrfButNoneAvailable";
    case ResultCode::SolveFailedToSetup:                   return "SolveFailedToSetup";
    case ResultCode::SolveFailedToSolve:                   return "SolveFailedToSolve";
    case ResultCode::SolveUserCanceled:                    return "SolveUserCanceled";
    case ResultCode::SolveThrewException:                  return "SolveThrewException";
    case ResultCode::UnknownStatus:                        return "UnknownStatus";
  }//switch( code )

  return "InvalidResultCode";
}//const char *to_str( const ResultCode code )


shared_ptr<RelActCalcAuto::RelActAutoGuiState>
  load_state_from_xml_file( const std::string &filename )
{
  if( filename.empty() )
    throw runtime_error( "load_state_from_xml_file: empty filename" );

  if( !SpecUtils::is_file(filename) )
    throw runtime_error( "load_state_from_xml_file: file does not exist: '" + filename + "'" );

  vector<char> data;
  try
  {
    SpecUtils::load_file_data( filename.c_str(), data );
  }catch( std::exception &e )
  {
    throw runtime_error( "load_state_from_xml_file: error reading '" + filename
                        + "': " + e.what() );
  }

  if( data.empty() )
    throw runtime_error( "load_state_from_xml_file: empty file '" + filename + "'" );

  // load_file_data null-terminates already; rapidxml needs this for in-place parse
  rapidxml::xml_document<char> doc;
  try
  {
    doc.parse<rapidxml::parse_trim_whitespace>( &data[0] );
  }catch( std::exception &e )
  {
    throw runtime_error( "load_state_from_xml_file: invalid XML in '" + filename
                        + "': " + e.what() );
  }

  const rapidxml::xml_node<char> *base_node = doc.first_node( "RelActCalcAuto" );
  if( !base_node )
    throw runtime_error( "load_state_from_xml_file: '" + filename
                        + "' has no <RelActCalcAuto> root element" );

  shared_ptr<RelActCalcAuto::RelActAutoGuiState> state
        = make_shared<RelActCalcAuto::RelActAutoGuiState>();

  try
  {
    state->deSerialize( base_node );
  }catch( std::exception &e )
  {
    throw runtime_error( "load_state_from_xml_file: failed to deserialize '" + filename
                        + "': " + e.what() );
  }

  return state;
}//load_state_from_xml_file(...)


// Pick the single foreground or background sample number from a SpecMeas, the same
// way `BatchActivity::fit_activities_in_file` does it.  Throws if the choice is
// ambiguous.
static int find_single_sample( const shared_ptr<const SpecMeas> &meas,
                               const SpecUtils::SourceType wanted )
{
  if( !meas )
    throw runtime_error( "find_single_sample: null SpecMeas" );

  const vector<string> &detectors = meas->detector_names();
  const set<int> &sample_nums = meas->sample_numbers();

  set<int> foreground_samples, background_samples;
  for( const int sample : sample_nums )
  {
    bool classified = false;
    for( const string &det : detectors )
    {
      const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det );
      if( !m )
        continue;

      switch( m->source_type() )
      {
        case SpecUtils::SourceType::IntrinsicActivity:
        case SpecUtils::SourceType::Calibration:
          break;
        case SpecUtils::SourceType::Background:
          classified = true;
          background_samples.insert( sample );
          break;
        case SpecUtils::SourceType::Foreground:
        case SpecUtils::SourceType::Unknown:
          classified = true;
          foreground_samples.insert( sample );
          break;
      }//switch

      if( classified )
        break;
    }//for( det )
  }//for( sample )

  const set<int> &samples = (wanted == SpecUtils::SourceType::Foreground) ? foreground_samples
                                                                          : background_samples;

  // If a file has only a single sample marked as background but we want foreground,
  // use it (mirrors BatchActivity behaviour).
  if( samples.empty()
      && (wanted == SpecUtils::SourceType::Foreground)
      && (background_samples.size() == 1) )
    return *begin(background_samples);

  if( samples.size() != 1 )
    throw runtime_error( "Sample number to use could not be uniquely identified" );

  return *begin(samples);
}//find_single_sample(...)


// Sum (or pass through) the measurements of a given sample-number set into a
// single Measurement.  Returns nullptr on failure.
static shared_ptr<const SpecUtils::Measurement>
    extract_measurement( const shared_ptr<const SpecMeas> &meas,
                         const set<int> &sample_nums )
{
  if( !meas || sample_nums.empty() )
    return nullptr;

  if( sample_nums.size() == 1 )
  {
    const vector<shared_ptr<const SpecUtils::Measurement>> sample_meass
                                              = meas->sample_measurements( *begin(sample_nums) );
    if( sample_meass.size() == 1 )
      return sample_meass[0];
  }

  try
  {
    return meas->sum_measurements( sample_nums, meas->detector_names(), nullptr );
  }catch( std::exception & )
  {
    return nullptr;
  }
}//extract_measurement(...)


static bool fwhm_method_requires_drf( const RelActCalcAuto::FwhmEstimationMethod m )
{
  return (m == RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency)
         || (m == RelActCalcAuto::FwhmEstimationMethod::StartingFromDetectorEfficiency);
}


// Build the on-disk output filename for a rendered template.  Mirrors
// BatchInfoLog::suggested_output_report_filename's name-cleanup logic, but
// keeps the isotopics path independent of BatchPeak's option struct.
//   - filename:   input spectrum filename (empty for summary reports)
//   - tmplt:      shorthand ("html","txt","json") or template path
//   - is_summary: true → use "rel_eff_summary"; false → use "rel_eff" as base
//   - output_dir: target directory for the file
static std::string isotopics_output_filename( const std::string &filename,
                                              const std::string &tmplt,
                                              const bool is_summary,
                                              const std::string &output_dir )
{
  std::string outname = SpecUtils::filename( filename );
  const std::string file_ext = SpecUtils::file_extension( outname );
  if( !file_ext.empty() )
    outname = outname.substr( 0, outname.size() - file_ext.size() );

  std::string tmplt_name = SpecUtils::filename( tmplt );
  std::string tmplt_ext = SpecUtils::file_extension( tmplt_name );

  const bool shorthand = SpecUtils::iequals_ascii( tmplt, "html" )
                      || SpecUtils::iequals_ascii( tmplt, "txt" )
                      || SpecUtils::iequals_ascii( tmplt, "text" )
                      || SpecUtils::iequals_ascii( tmplt, "json" )
                      || SpecUtils::iequals_ascii( tmplt, "html-summary" )
                      || SpecUtils::iequals_ascii( tmplt, "summary" );
  if( shorthand )
  {
    tmplt_name = is_summary ? "rel_eff_summary" : "rel_eff";
    if( SpecUtils::iequals_ascii(tmplt, "html-summary")
        || SpecUtils::iequals_ascii(tmplt, "summary") )
      tmplt_ext = ".html";
    else if( SpecUtils::iequals_ascii(tmplt, "text") )
      tmplt_ext = ".txt";
    else
    {
      tmplt_ext = "." + tmplt;
      SpecUtils::to_lower_ascii( tmplt_ext );
    }
  }
  else
  {
    // Strip "_tmplt" / "_template" markers from filename to get a clean base
    size_t pos = SpecUtils::ifind_substr_ascii( tmplt_name, "tmplt" );
    if( pos == std::string::npos )
      pos = SpecUtils::ifind_substr_ascii( tmplt_name, "template" );
    if( pos != std::string::npos )
      tmplt_name = tmplt_name.substr( 0, pos );
    if( SpecUtils::iends_with(tmplt_name, "_")
        || SpecUtils::iends_with(tmplt_name, ".")
        || SpecUtils::iends_with(tmplt_name, "-") )
      tmplt_name = tmplt_name.substr( 0, tmplt_name.size() - 1 );

    if( tmplt_ext.empty()
        || SpecUtils::iequals_ascii(tmplt_ext, "tmplt")
        || SpecUtils::iequals_ascii(tmplt_ext, "template") )
      tmplt_ext = SpecUtils::file_extension( tmplt_name );
    if( tmplt_ext.empty() )
      tmplt_ext = ".txt";
  }

  outname += (outname.empty() ? "" : "_") + tmplt_name + tmplt_ext;
  return SpecUtils::append_path( output_dir, outname );
}//isotopics_output_filename(...)


Result run_on_file( const std::string &exemplar_filename,
                    std::set<int> exemplar_sample_nums,
                    std::shared_ptr<const SpecMeas> cached_exemplar,
                    const std::string &filename,
                    std::shared_ptr<SpecMeas> cached_file,
                    const Options &options )
{
  Result result;
  result.m_options = options;
  result.m_filename = filename;
  result.m_exemplar_sample_numbers = exemplar_sample_nums;
  result.m_result_code = ResultCode::UnknownStatus;

  // ---- 1. Load exemplar ------------------------------------------------------
  if( !cached_exemplar )
  {
    if( exemplar_filename.empty() && !options.state_override )
    {
      result.m_error_msg = "No exemplar file or state-override XML supplied.";
      result.m_result_code = ResultCode::NoExemplar;
      return result;
    }

    if( !exemplar_filename.empty() )
    {
      shared_ptr<SpecMeas> tmp = make_shared<SpecMeas>();
      if( !tmp->load_file( exemplar_filename, SpecUtils::ParserType::Auto ) )
      {
        result.m_error_msg = "Could not load exemplar '" + exemplar_filename + "'.";
        result.m_result_code = ResultCode::CouldntOpenExemplar;
        return result;
      }
      cached_exemplar = tmp;
    }
  }//if( !cached_exemplar )

  result.m_exemplar_file = cached_exemplar;

  // ---- 2. Resolve RelActAutoGuiState ----------------------------------------
  shared_ptr<RelActCalcAuto::RelActAutoGuiState> state;

  if( options.state_override )
  {
    // Make a mutable copy so we can apply scalar overrides on top.
    state = make_shared<RelActCalcAuto::RelActAutoGuiState>( *options.state_override );
  }else if( cached_exemplar )
  {
    unique_ptr<RelActCalcAuto::RelActAutoGuiState> tmp = cached_exemplar->getRelActAutoGuiState();
    if( !tmp )
    {
      result.m_error_msg = "Exemplar '" + exemplar_filename
                           + "' does not contain a stored <RelActCalcAuto> state.";
      result.m_result_code = ResultCode::ExemplarMissingRelActState;
      return result;
    }
    state = make_shared<RelActCalcAuto::RelActAutoGuiState>( std::move(*tmp) );
  }
  else
  {
    result.m_error_msg = "No exemplar and no state-override supplied.";
    result.m_result_code = ResultCode::NoExemplar;
    return result;
  }

  assert( state );

  // ---- 3. Apply scalar overrides --------------------------------------------
  if( options.energy_cal_type )
    state->options.energy_cal_type = *options.energy_cal_type;
  if( options.fwhm_form )
    state->options.fwhm_form = *options.fwhm_form;
  if( options.skew_type )
    state->options.skew_type = *options.skew_type;
  if( options.background_subtract )
    state->background_subtract = *options.background_subtract;

  // ---- 4. Load foreground spectrum ------------------------------------------
  shared_ptr<SpecMeas> specfile = cached_file;
  if( !specfile )
  {
    specfile = make_shared<SpecMeas>();
    if( !specfile->load_file( filename, SpecUtils::ParserType::Auto ) )
    {
      result.m_error_msg = "Could not load foreground '" + filename + "'.";
      result.m_result_code = ResultCode::CouldntOpenInputFile;
      return result;
    }
  }
  result.m_foreground_file = specfile;

  shared_ptr<const SpecUtils::Measurement> foreground;
  set<int> foreground_sample_numbers;
  try
  {
    const int sample_num = find_single_sample( specfile, SpecUtils::SourceType::Foreground );
    foreground_sample_numbers.insert( sample_num );
    foreground = extract_measurement( specfile, foreground_sample_numbers );
    if( !foreground )
      throw runtime_error( "Missing foreground measurement." );
  }catch( std::exception &e )
  {
    result.m_error_msg = string("Could not determine foreground: ") + e.what();
    result.m_result_code = ResultCode::ForegroundSampleNumberUnderSpecified;
    return result;
  }
  result.m_foreground_sample_numbers = foreground_sample_numbers;

  // ---- 5. Resolve background --------------------------------------------------
  // Precedence (decisions confirmed with user):
  //   a) explicit --back-sub-file (or cached_background_subtract_spec) → use it
  //      and force state->background_subtract = true
  //   b) options.background_subtract == false → force off
  //   c) state->background_subtract == true → reuse exemplar's stored background
  //      (look up state->options.background_filename / sample_numbers)
  //   d) otherwise → no background
  shared_ptr<const SpecUtils::Measurement> background;
  shared_ptr<const SpecMeas> background_file;
  set<int> background_sample_numbers;

  const bool have_explicit_back = !options.background_subtract_file.empty()
                                  || (options.cached_background_subtract_spec != nullptr);

  if( have_explicit_back )
  {
    state->background_subtract = true;

    shared_ptr<SpecMeas> backfile = options.cached_background_subtract_spec;
    if( !backfile )
    {
      if( options.background_subtract_file == filename )
      {
        backfile = specfile;
      }else if( options.background_subtract_file == exemplar_filename && cached_exemplar )
      {
        // Keep null backfile sentinel — we'll source the measurement from cached_exemplar below.
      }else
      {
        backfile = make_shared<SpecMeas>();
        if( !backfile->load_file( options.background_subtract_file, SpecUtils::ParserType::Auto ) )
        {
          result.m_error_msg = "Could not load background '" + options.background_subtract_file + "'.";
          result.m_result_code = ResultCode::CouldntOpenBackgroundFile;
          return result;
        }
      }
    }//if( !backfile )

    const shared_ptr<const SpecMeas> source =
        backfile ? std::const_pointer_cast<const SpecMeas>(backfile)
                 : cached_exemplar;
    if( source )
    {
      try
      {
        if( !options.background_subtract_samples.empty() )
        {
          background_sample_numbers = options.background_subtract_samples;
        }else
        {
          const int n = find_single_sample( source, SpecUtils::SourceType::Background );
          background_sample_numbers.insert( n );
        }
        background = extract_measurement( source, background_sample_numbers );
        if( !background )
          throw runtime_error( "could not extract background measurement" );
      }catch( std::exception &e )
      {
        result.m_error_msg = string("Could not determine background: ") + e.what();
        result.m_result_code = ResultCode::BackgroundSampleNumberUnderSpecified;
        return result;
      }
      background_file = source;
    }

    // Update state's bookkeeping fields so reports show the right filenames.
    state->options.background_filename = options.background_subtract_file;
    state->options.background_sample_numbers = background_sample_numbers;
  }
  else if( options.background_subtract && !*options.background_subtract )
  {
    state->background_subtract = false;
  }
  else if( state->background_subtract )
  {
    // Use exemplar's stored background.  First, try to find it inside the
    // exemplar SpecMeas (typical case — the exemplar N42 includes both).
    shared_ptr<const SpecMeas> source = cached_exemplar;
    set<int> bk_samples = state->options.background_sample_numbers;

    auto try_extract = [&]() -> bool {
      if( !source )
        return false;
      try
      {
        if( bk_samples.empty() )
        {
          const int n = find_single_sample( source, SpecUtils::SourceType::Background );
          bk_samples.insert( n );
        }
        background = extract_measurement( source, bk_samples );
      }catch( std::exception & )
      {
        background = nullptr;
      }
      return !!background;
    };

    bool ok = try_extract();
    if( !ok && !state->options.background_filename.empty() )
    {
      // Background lives in a separate file referenced by the exemplar.
      shared_ptr<SpecMeas> tmp = make_shared<SpecMeas>();
      if( tmp->load_file( state->options.background_filename, SpecUtils::ParserType::Auto ) )
      {
        source = tmp;
        bk_samples = state->options.background_sample_numbers;
        ok = try_extract();
      }
    }

    if( !ok )
    {
      result.m_error_msg = "Exemplar requested background subtraction, but the background"
                           " measurement could not be located.";
      result.m_result_code = ResultCode::BackgroundSampleNumberUnderSpecified;
      return result;
    }

    background_file = source;
    background_sample_numbers = bk_samples;
  }
  result.m_background_file = background_file;
  result.m_background_sample_numbers = background_sample_numbers;

  // ---- 6. Resolve DRF and validate ------------------------------------------
  shared_ptr<const DetectorPeakResponse> drf = options.drf_override;
  if( !drf && cached_exemplar )
    drf = cached_exemplar->detector();

  if( !drf && fwhm_method_requires_drf( state->options.fwhm_estimation_method )  )
  {
    result.m_error_msg = "The exemplar / state's FwhmEstimationMethod is "
                          + string( RelActCalcAuto::to_str(state->options.fwhm_estimation_method) )
                         + ", which requires a DRF, but none was provided.";
    result.m_result_code = ResultCode::FwhmMethodNeedsDrfButNoneAvailable;
    return result;
  }

  // ---- 7. Solve --------------------------------------------------------------
  RelActCalcAuto::RelActAutoSolution solution;
  try
  {
    solution = RelActCalcAuto::solve( state->options, foreground, background, drf,
                                      /*all_peaks=*/{}, /*cancel_calc=*/nullptr );
  }catch( std::exception &e )
  {
    result.m_error_msg = string("RelActCalcAuto::solve threw: ") + e.what();
    result.m_result_code = ResultCode::SolveThrewException;
    return result;
  }

  // Move warnings out so the per-file reporter sees them.
  result.m_warnings = solution.m_warnings;

  switch( solution.m_status )
  {
    case RelActCalcAuto::RelActAutoSolution::Status::Success:
      result.m_result_code = ResultCode::Success;
      break;
    case RelActCalcAuto::RelActAutoSolution::Status::NotInitiated:
      result.m_result_code = ResultCode::SolveFailedToSetup;
      result.m_error_msg = solution.m_error_message;
      break;
    case RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem:
      result.m_result_code = ResultCode::SolveFailedToSetup;
      result.m_error_msg = solution.m_error_message;
      break;
    case RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem:
      result.m_result_code = ResultCode::SolveFailedToSolve;
      result.m_error_msg = solution.m_error_message;
      break;
    case RelActCalcAuto::RelActAutoSolution::Status::UserCanceled:
      result.m_result_code = ResultCode::SolveUserCanceled;
      result.m_error_msg = solution.m_error_message;
      break;
  }//switch

  result.m_solution = std::move(solution);
  return result;
}//run_on_file(...)


void run_in_files( const std::string &exemplar_filename,
                   std::shared_ptr<const SpecMeas> cached_exemplar,
                   const std::set<int> &exemplar_sample_nums,
                   const std::vector<std::string> &files,
                   std::vector<std::shared_ptr<SpecMeas>> cached_files,
                   const Options &options,
                   Summary *summary )
{
  vector<string> warnings;

  if( files.empty() )
    throw runtime_error( "No input files specified." );

  if( !cached_files.empty() && (cached_files.size() != files.size()) )
    throw runtime_error( "If you specify cached files, you must specify the same number of files." );

  if( !options.output_dir.empty() && !SpecUtils::is_directory(options.output_dir) )
    throw runtime_error( "Output directory ('" + options.output_dir + "'), is not a directory." );

  if( options.write_n42_with_results && options.output_dir.empty() )
    throw runtime_error( "If you specify to write N42 files with results, you must specify an output directory." );

  // ---- Set up Inja env (templates resolved against include_dir → writable → bundled)
  // We keep the un-resolved `options.template_include_dir` (could be "", "none", "default",
  // or a real path) and pass it through to `render_template` -- it needs to distinguish the
  // user-explicit case from the CLI default to know whether to validate template paths.
  const string user_include_dir = options.template_include_dir;
  const string tmplt_dir = RelActAutoReport::template_include_dir( user_include_dir );
  inja::Environment env = RelActAutoReport::get_default_inja_env( tmplt_dir );

  nlohmann::json summary_json;
  summary_json["ExemplarFile"] = exemplar_filename;
  if( !exemplar_sample_nums.empty() )
    summary_json["ExemplarSampleNumbers"]
        = vector<int>{ begin(exemplar_sample_nums), end(exemplar_sample_nums) };
  summary_json["InputFiles"] = files;
  summary_json["Files"] = nlohmann::json::array();

  if( summary )
  {
    summary->options = options;
    summary->exemplar_filename = exemplar_filename;
    summary->exemplar = cached_exemplar;
    summary->exemplar_sample_nums = exemplar_sample_nums;
  }

  for( size_t i = 0; i < files.size(); ++i )
  {
    const string &fname = files[i];
    const shared_ptr<SpecMeas> cached = cached_files.empty() ? nullptr : cached_files[i];

    Result file_result = run_on_file( exemplar_filename, exemplar_sample_nums,
                                      cached_exemplar, fname, cached, options );

    // Hard-stop errors that aren't per-file: bad exemplar / bad background.
    if( (file_result.m_result_code == ResultCode::CouldntOpenExemplar)
        || (file_result.m_result_code == ResultCode::CouldntOpenStateOverride) )
      throw runtime_error( file_result.m_error_msg );

    if( !cached_exemplar )
      cached_exemplar = file_result.m_exemplar_file;

    for( const string &w : file_result.m_warnings )
      warnings.push_back( "File '" + SpecUtils::filename(fname) + "': " + w );

    // Build per-file JSON via RelActAutoReport, then add bookkeeping for templates.
    nlohmann::json file_json = RelActAutoReport::solution_to_json( file_result.m_solution );
    file_json["Filepath"] = fname;
    file_json["Filename"] = SpecUtils::filename(fname);
    file_json["ParentDir"] = SpecUtils::parent_path(fname);
    file_json["ResultCode"] = to_str( file_result.m_result_code );
    file_json["ResultCodeInt"] = static_cast<int>( file_result.m_result_code );
    file_json["Success"] = (file_result.m_result_code == ResultCode::Success);
    file_json["HasErrorMessage"] = !file_result.m_error_msg.empty();
    if( !file_result.m_error_msg.empty() )
      file_json["ErrorMessage"] = file_result.m_error_msg;
    file_json["HasWarnings"] = !file_result.m_warnings.empty();
    if( !file_result.m_warnings.empty() )
      file_json["Warnings"] = file_result.m_warnings;

    const bool success = (file_result.m_result_code == ResultCode::Success);
    if( success )
      cout << "Success analyzing '" << fname << "'!" << endl;
    else
      cout << "Failure analyzing '" << fname << "': " << file_result.m_error_msg << endl;

    // Per-file reports - rendered with assets present in `file_json`.  We
    // strip them after rendering and lift `assets` to `summary_json["assets"]`
    // once (mirrors BatchActivity's pattern), so `summary_json["Files"][i]`
    // doesn't carry N copies of hundreds of KB of JS/CSS.
    if( summary )
    {
      summary->file_results.push_back( file_result );
      summary->file_json.push_back( file_json.dump() );
      summary->file_reports.push_back( vector<string>() );
    }

    for( const string &tmplt : options.report_templates )
    {
      try
      {
        const string rpt = RelActAutoReport::render_template( env, file_json, tmplt, user_include_dir );

        if( summary )
          summary->file_reports.back().push_back( rpt );

        if( options.to_stdout && !SpecUtils::iequals_ascii(tmplt, "html") )
          cout << "\n\n" << rpt << endl << endl;

        if( !options.output_dir.empty() )
        {
          const string out_file = isotopics_output_filename( fname, tmplt,
                                                             /*is_summary=*/false,
                                                             options.output_dir );

          if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
          {
            warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                                " See the '--overwrite-output-files' option to force writing." );
          }else
          {
#ifdef _WIN32
            const std::wstring wout = SpecUtils::convert_from_utf8_to_utf16(out_file);
            std::ofstream output( wout.c_str(), ios::binary | ios::out );
#else
            std::ofstream output( out_file.c_str(), ios::binary | ios::out );
#endif
            if( !output )
              warnings.push_back( "Failed to open report output '" + out_file + "', for writing." );
            else
              output.write( rpt.c_str(), rpt.size() );
          }
        }
      }catch( inja::InjaError &e )
      {
        const string msg = "Error templating results (" + e.type + ": line "
                          + std::to_string(e.location.line) + ", column "
                          + std::to_string(e.location.column) + " of '" + tmplt + "'): "
                          + e.message + ". While processing '" + fname + "'.";
        cerr << msg << endl;
        warnings.push_back( msg );
      }catch( std::exception &e )
      {
        cerr << "Error templating results: " << e.what() << endl;
        warnings.push_back( string("Error templating results: ") + e.what() );
      }
    }//for each per-file template

    // Lift `assets` to summary_json["assets"] once, then strip from per-file
    // JSON so `summary_json["Files"][i]` doesn't carry duplicate hundreds-of-KB
    // blobs.  Mirrors BatchActivity::fit_activities_in_files.
    if( !summary_json.contains("assets")
        && file_json.contains("assets")
        && !file_json["assets"].empty() )
    {
      summary_json["assets"] = file_json["assets"];
    }
    file_json.erase( "assets" );

    summary_json["Files"].push_back( file_json );
  }//for each input file

  // Summary-level counts and metadata used by the bundled multi-file summary
  // template.
  size_t num_succeeded = 0, num_failed = 0;
  for( const auto &f : summary_json["Files"] )
  {
    if( f.value("Success", false) )
      ++num_succeeded;
    else
      ++num_failed;
  }
  summary_json["NumFiles"] = summary_json["Files"].size();
  summary_json["NumSucceeded"] = num_succeeded;
  summary_json["NumFailed"] = num_failed;

  for( const string &w : warnings )
    summary_json["Warnings"].push_back( w );

  // ---- Summary reports ------------------------------------------------------
  // Bundled-shorthand mapping for the summary path: "html" / "" map to the
  // multi-file `html-summary` template, which walks `summary_json["Files"]`.
  // Custom templates (filenames / absolute paths) flow through unchanged, as
  // does "json".
  for( const string &tmplt_in : options.summary_report_templates )
  {
    string tmplt = tmplt_in;
    if( tmplt.empty() || SpecUtils::iequals_ascii(tmplt, "html") )
      tmplt = "html-summary";

    try
    {
      const string rpt = RelActAutoReport::render_template( env, summary_json, tmplt, user_include_dir );

      if( summary )
        summary->summary_reports.push_back( rpt );

      if( options.to_stdout && !SpecUtils::iequals_ascii(tmplt, "html") )
        cout << "\n\n" << rpt << endl << endl;

      if( !options.output_dir.empty() )
      {
        const string out_file = isotopics_output_filename( "", tmplt,
                                                            /*is_summary=*/true,
                                                            options.output_dir );

        if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
        {
          warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                              " See the '--overwrite-output-files' option to force writing." );
        }else
        {
#ifdef _WIN32
          const std::wstring wout = SpecUtils::convert_from_utf8_to_utf16(out_file);
          std::ofstream output( wout.c_str(), ios::binary | ios::out );
#else
          std::ofstream output( out_file.c_str(), ios::binary | ios::out );
#endif
          if( !output )
            warnings.push_back( "Failed to open summary report output, '" + out_file + "'" );
          else
            output.write( rpt.c_str(), rpt.size() );
        }
      }
    }catch( inja::InjaError &e )
    {
      const string msg = "Error templating summary output (" + e.type + ": line "
                        + std::to_string(e.location.line) + ", column "
                        + std::to_string(e.location.column) + " of '" + tmplt + "'): "
                        + e.message + ".";
      cerr << msg << endl;
      warnings.push_back( msg );
    }catch( std::exception &e )
    {
      warnings.push_back( string("Error making summary output: ") + e.what() );
    }
  }//for each summary template

  if( !warnings.empty() )
    cerr << endl;
  for( const string &w : warnings )
    cerr << w << endl;

  if( summary )
  {
    summary->summary_json = summary_json.dump();
    summary->warnings = warnings;
  }
}//run_in_files(...)

}//namespace BatchRelActAuto
