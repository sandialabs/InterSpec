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
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <filesystem>
#include <system_error>

#include <Wt/WRandom.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/PeakModel.h"
#include "target/mcp_http/McpAgentTools.h"
#include "InterSpec/ExportSpecFile.h"
#include "target/mcp_http/McpSessionManager.h"

using namespace std;
using json = nlohmann::json;
namespace fs = std::filesystem;


namespace
{
  /** Reads a text file from the LLM config directory, or else the static data directory. */
  std::string read_config_or_data_file( const std::string &relpath )
  {
    vector<string> candidates;
    const string config_dir = LlmConfig::configDirectory();
    if( !config_dir.empty() )
      candidates.push_back( SpecUtils::append_path( config_dir, relpath ) );
    try
    {
      candidates.push_back( SpecUtils::append_path( InterSpec::staticDataDirectory(), relpath ) );
    }catch( std::exception & )
    {
    }

    for( const string &path : candidates )
    {
      if( !SpecUtils::is_file( path ) )
        continue;
      std::ifstream input( path.c_str(), ios::in | ios::binary );
      return string( std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>() );
    }

    return string();
  }//read_config_or_data_file(...)


  /** Writes `data` to `<session>/<sub_dir>/` as a uniquely named file: first to a temporary name,
   then renamed into place, so a reader never sees a partial file.  Returns the path relative to the
   session directory.
   */
  std::string write_session_file( const std::string &handle, const std::string &sub_dir,
                                  const std::string &filename, std::istream &in,
                                  const uint64_t max_bytes, uint64_t &num_bytes )
  {
    num_bytes = 0;
    const shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
    if( !manager )
      throw runtime_error( "Not running as an MCP agent server." );

    const McpSessionManager::AgentSettings &settings = manager->settings();
    const uint64_t used = manager->sessionBytesUsed( handle );
    if( used >= settings.maxSessionBytes )
      throw runtime_error( "The session's storage limit (" + std::to_string(settings.maxSessionBytes)
                           + " bytes) has been reached." );
    const uint64_t limit = std::min( max_bytes, settings.maxSessionBytes - used );

    const string dir = manager->resolveSessionPath( handle, sub_dir );
    const string tmp_path = SpecUtils::append_path( dir, ".partial-" + Wt::WRandom::generateId(16) );

    {
      std::ofstream output( tmp_path.c_str(), ios::out | ios::binary | ios::trunc );
      if( !output.is_open() )
        throw runtime_error( "Could not create the file." );

      char buffer[64*1024];
      while( in )
      {
        in.read( buffer, sizeof(buffer) );
        const std::streamsize nread = in.gcount();
        if( nread <= 0 )
          break;
        num_bytes += static_cast<uint64_t>( nread );
        if( num_bytes > limit )
        {
          output.close();
          std::error_code ec;
          fs::remove( tmp_path, ec );
          throw runtime_error( "The file is larger than allowed (limit " + std::to_string(limit)
                               + " bytes, given the session's remaining storage)." );
        }
        output.write( buffer, nread );
      }//while( in )

      if( !output.good() )
      {
        output.close();
        std::error_code ec;
        fs::remove( tmp_path, ec );
        throw runtime_error( "Error writing the file." );
      }
    }

    const string name = McpSessionManager::uniqueFilename( dir, filename );
    std::error_code ec;
    fs::rename( tmp_path, SpecUtils::append_path( dir, name ), ec );
    if( ec )
    {
      fs::remove( tmp_path, ec );
      throw runtime_error( "Could not save the file." );
    }

    return sub_dir + "/" + name;
  }//write_session_file(...)


  /** The user-session argument added to the session tools' schemas. */
  json user_session_property()
  {
    return json{
      {"type", "string"},
      {"description", "Required: the userSession value returned by create_session."}
    };
  }
}//namespace


namespace McpAgentTools
{

const std::vector<ToolDef> &toolDefinitions()
{
  static const std::vector<ToolDef> sm_tools = [](){
    vector<ToolDef> tools;

    {
      ToolDef tool;
      tool.name = "create_session";
      tool.description =
        "Call this first.  Creates a new InterSpec analysis session, and returns its 'userSession'"
        " handle, which every other session tool requires.  Each session has its own spectra,"
        " analysis, and assistant conversation, plus a directory for exchanging files: put input files"
        " in its \"input/\" directory - either by writing into session_dir (if the exchange directory"
        " is shared with you), or by HTTP PUT of the file to <upload_url>/<filename> (with the same"
        " bearer token), or with upload_file - then load them with load_spectrum_file using paths like"
        " \"input/<filename>\".  Results are written to \"output/\", each with a download_url.  Call"
        " close_session when done; a session also ends after a period of inactivity, or at its"
        " maximum lifetime (expires_at), and its files are then deleted unless the server keeps them.";
      tool.schema = json{ {"type", "object"}, {"properties", json::object()} };
      tools.push_back( tool );
    }

    {
      ToolDef tool;
      tool.name = "close_session";
      tool.description =
        "Ends a session created with create_session; its URLs stop working, and (unless the server is"
        " configured to keep them) its files are deleted - so download or export_analysis anything you"
        " need first.";
      tool.schema = json{
        {"type", "object"},
        {"properties", { {"userSession", user_session_property()} }},
        {"required", json::array({"userSession"})}
      };
      tool.needs_session = true;
      tools.push_back( tool );
    }

    {
      ToolDef tool;
      tool.name = "upload_file";
      tool.description =
        "Saves a (base64 encoded) file into the session's \"input/\" directory, returning its path to"
        " use with load_spectrum_file, or with assistant_submit_prompt's images.  For large files,"
        " prefer an HTTP PUT to the session's upload_url.";
      tool.schema = json{
        {"type", "object"},
        {"properties", {
          {"userSession", user_session_property()},
          {"filename", { {"type", "string"}, {"description", "Name for the file; characters other than letters, digits, '.', '_', and '-' become '_', and a number is added if the name is taken."} }},
          {"base64Data", { {"type", "string"}, {"description", "The file contents, base64 encoded."} }}
        }},
        {"required", json::array({"userSession", "filename", "base64Data"})}
      };
      tool.needs_session = true;
      tools.push_back( tool );
    }

    {
      ToolDef tool;
      tool.name = "list_session_files";
      tool.description =
        "Lists the files in the session's \"input/\" and \"output/\" directories, with a download_url"
        " for each output file.";
      tool.schema = json{
        {"type", "object"},
        {"properties", { {"userSession", user_session_property()} }},
        {"required", json::array({"userSession"})}
      };
      tool.needs_session = true;
      tools.push_back( tool );
    }

    {
      ToolDef tool;
      tool.name = "export_analysis";
      tool.description =
        "Writes the session's work to its \"output/\" directory, returning each file's path and"
        " download_url.  \"n42\" is the foreground spectrum file with all of the analysis - peaks,"
        " the activity/shielding fit, isotopics setup, detector response, energy calibration, and the"
        " assistant conversation (including images sent with prompts) - plus the background if it is"
        " from a different file.  Loading that N42 into a new session with load_spectrum_file and"
        " conversation=\"resume\" continues the work.  \"peaks_csv\" is the peak table, and"
        " \"spectrum_html\" a standalone interactive plot.  Not available while an assistant turn is"
        " running.";
      tool.schema = json{
        {"type", "object"},
        {"properties", {
          {"userSession", user_session_property()},
          {"items", {
            {"type", "array"},
            {"items", { {"type", "string"}, {"enum", json::array({"n42", "peaks_csv", "spectrum_html"})} }},
            {"description", "What to export; default [\"n42\"]."}
          }},
          {"basename", { {"type", "string"}, {"description", "Optional start of the file names; default is from the foreground file name."} }}
        }},
        {"required", json::array({"userSession"})}
      };
      tool.needs_session = true;
      tool.runs_in_session = true;
      tools.push_back( tool );
    }

    {
      ToolDef tool;
      tool.name = "get_usage_guide";
      tool.description =
        "Returns the guide to using this server: the session workflow, exchanging files, working with"
        " the assistant, exporting and resuming, and how to present results.";
      tool.schema = json{ {"type", "object"}, {"properties", json::object()} };
      tools.push_back( tool );
    }

    return tools;
  }();

  return sm_tools;
}//toolDefinitions()


const ToolDef *toolDefinition( const std::string &name )
{
  for( const ToolDef &tool : toolDefinitions() )
  {
    if( tool.name == name )
      return &tool;
  }
  return nullptr;
}//toolDefinition(...)


std::string serverInstructions()
{
  return read_config_or_data_file( "mcp_agent_instructions.md" );
}


std::string usageGuide()
{
  return read_config_or_data_file( SpecUtils::append_path( SpecUtils::append_path( "skills", "interspec-mcp-agent" ), "SKILL.md" ) );
}


std::string saveInputFile( const std::string &handle, const std::string &filename,
                           std::istream &in, uint64_t &num_bytes )
{
  const shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( !manager )
    throw runtime_error( "Not running as an MCP agent server." );

  return write_session_file( handle, "input", filename, in, manager->settings().maxUploadBytes, num_bytes );
}//saveInputFile(...)


std::string saveOutputFile( const std::string &handle, const std::string &filename,
                            const std::string &data )
{
  const shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( !manager )
    throw runtime_error( "Not running as an MCP agent server." );

  std::istringstream in( data );
  uint64_t num_bytes = 0;
  return write_session_file( handle, "output", filename, in, manager->settings().maxSessionBytes, num_bytes );
}//saveOutputFile(...)


nlohmann::json listSessionFiles( const std::string &handle )
{
  const shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( !manager )
    throw runtime_error( "Not running as an MCP agent server." );

  json answer;
  for( const char * const sub_dir : { "input", "output" } )
  {
    json files = json::array();
    const string dir = manager->resolveSessionPath( handle, sub_dir );
    std::error_code ec;
    for( fs::directory_iterator iter( dir, ec ), iter_end; !ec && (iter != iter_end); iter.increment(ec) )
    {
      const string name = iter->path().filename().string();
      std::error_code status_ec;
      if( name.empty() || (name[0] == '.') || !fs::is_regular_file( fs::symlink_status( iter->path(), status_ec ) ) )
        continue;
      files.push_back( json{ {"path", string(sub_dir) + "/" + name},
                             {"bytes", static_cast<uint64_t>( iter->file_size( status_ec ) )} } );
    }//for( loop over directory )
    answer[sub_dir] = std::move( files );
  }//for( input and output )

  return answer;
}//listSessionFiles(...)


nlohmann::json exportAnalysis( InterSpec *viewer, const std::string &handle, const nlohmann::json &args )
{
  if( !viewer )
    throw runtime_error( "No InterSpec session." );

  const shared_ptr<SpecMeas> foreground = viewer->measurment( SpecUtils::SpectrumType::Foreground );
  if( !foreground )
    throw runtime_error( "No foreground spectrum is loaded; load one with load_spectrum_file first." );

  vector<string> items;
  if( args.contains( "items" ) && args["items"].is_array() )
  {
    for( const json &item : args["items"] )
    {
      if( item.is_string() )
        items.push_back( item.get<string>() );
    }
  }
  if( items.empty() )
    items.push_back( "n42" );

  for( const string &item : items )
  {
    if( (item != "n42") && (item != "peaks_csv") && (item != "spectrum_html") )
      throw runtime_error( "Unknown export item '" + item + "'; use n42, peaks_csv, or spectrum_html." );
  }

  string basename = args.value( "basename", string() );
  if( basename.empty() )
  {
    basename = SpecUtils::filename( foreground->filename() );
    const size_t dot = basename.find_last_of( '.' );
    if( (dot != string::npos) && (dot > 0) )
      basename = basename.substr( 0, dot );
  }
  basename = McpSessionManager::sanitizeFilename( basename.empty() ? string("analysis") : basename );

  // Bring the state open tools hold into the SpecMeas, so the N42 has everything.
  viewer->syncAllToolStateToForegroundSpecMeas();

  json files = json::array();
  const auto add_file = [&files,&handle]( const string &item, const string &filename, const string &data ){
    const string path = saveOutputFile( handle, filename, data );
    files.push_back( json{ {"path", path}, {"item", item}, {"bytes", static_cast<uint64_t>(data.size())} } );
  };

  for( const string &item : items )
  {
    if( item == "n42" )
    {
      // The live SpecMeas is written, not a copy: copies renumber samples (which would break the
      //  saved conversation's association with them).
      std::ostringstream n42;
      if( !foreground->write_2012_N42( n42 ) )
        throw runtime_error( "Failed to write the N42 file." );
      add_file( item, basename + "_analysis.n42", n42.str() );

      const shared_ptr<SpecMeas> background = viewer->measurment( SpecUtils::SpectrumType::Background );
      if( background && (background != foreground) )
      {
        std::ostringstream back_n42;
        if( !background->write_2012_N42( back_n42 ) )
          throw runtime_error( "Failed to write the background N42 file." );
        add_file( item, basename + "_background.n42", back_n42.str() );
        files.back()["note"] = "The background spectrum, from a different file; when resuming, load it"
                               " with specType Background after the analysis N42.";
      }
    }else if( item == "peaks_csv" )
    {
      const shared_ptr<PeakModel> peak_model = viewer->peakModel();
      if( !peak_model )
        throw runtime_error( "No peaks are available." );
      std::ostringstream csv;
      peak_model->writePeaksCsv( csv );
      add_file( item, basename + "_peaks.csv", csv.str() );
    }else if( item == "spectrum_html" )
    {
#if( SpecUtils_ENABLE_D3_CHART )
      std::ostringstream html;
      write_displayed_spectra_d3_html( html, viewer );
      add_file( item, basename + "_spectrum.html", html.str() );
#else
      throw runtime_error( "HTML export is not available in this build." );
#endif
    }
  }//for( const string &item : items )

  return json{ {"files", files} };
}//exportAnalysis(...)

}//namespace McpAgentTools
