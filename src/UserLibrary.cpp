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

#include <map>
#include <set>
#include <ctime>
#include <mutex>
#include <atomic>
#include <random>
#include <regex>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <filesystem>

#include "rapidxml/rapidxml.hpp"

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"
#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "io/EnvFitArchive.h"

#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfLibrary.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/UserLibrary.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/ReferenceLinePredef.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ShieldingSourceDisplay.h"

#if( USE_REL_ACT_TOOL )
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActAutoReport.h"
#endif

#if( USE_LLM_INTERFACE )
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmKnowledgeBase.h"
#endif

using namespace std;
using nlohmann::json;


namespace
{
  /** Folders deeper than this, below a type's folder, are not searched. */
  const int ns_max_depth = 12;

  /** Files larger than this are not read as manifests or the history file. */
  const size_t ns_max_json_bytes = 16*1024*1024;

  std::filesystem::path to_fs_path( const string &utf8 )
  {
#ifdef _WIN32
    return std::filesystem::path( SpecUtils::convert_from_utf8_to_utf16( utf8 ) );
#else
    return std::filesystem::path( utf8 );
#endif
  }//to_fs_path(...)


  string from_fs_path( const std::filesystem::path &path )
  {
#ifdef _WIN32
    return SpecUtils::convert_from_utf16_to_utf8( path.wstring() );
#else
    return path.string();
#endif
  }//from_fs_path(...)


  bool is_dot_name( const string &path )
  {
    const string name = SpecUtils::filename( path );
    return name.empty() || (name[0] == '.');
  }


  string normalize( const string &path )
  {
    return SpecUtils::lexically_normalize_path( path );
  }


  /** Reads a whole file; false if it can not be read or is larger than `max_bytes`. */
  bool read_file( const string &path, string &contents, const size_t max_bytes )
  {
    contents.clear();
    std::ifstream input( to_fs_path( path ), ios::in | ios::binary );
    if( !input )
      return false;

    input.seekg( 0, ios::end );
    const std::streamoff len = input.tellg();
    if( (len < 0) || (static_cast<size_t>(len) > max_bytes) )
      return false;
    input.seekg( 0, ios::beg );

    contents.resize( static_cast<size_t>(len) );
    if( len && !input.read( &contents[0], len ) )
      return false;
    return true;
  }//read_file(...)


  /** A file's modification time and size, to tell whether a cached reading is still current. */
  struct Stamp
  {
    bool exists = false;
    int64_t mtime = 0;
    uint64_t size = 0;

    bool operator==( const Stamp &rhs ) const
    {
      return (exists == rhs.exists) && (mtime == rhs.mtime) && (size == rhs.size);
    }
  };//struct Stamp


  Stamp stamp_of( const string &path )
  {
    Stamp s;
    std::error_code ec;
    const std::filesystem::path p = to_fs_path( path );
    const std::filesystem::file_time_type t = std::filesystem::last_write_time( p, ec );
    if( ec )
      return s;
    // A folder has no size; its time changes as files are added or removed (or it is replaced).
    const uintmax_t size = std::filesystem::is_directory( p, ec ) ? 0 : std::filesystem::file_size( p, ec );
    if( ec )
      return s;
    s.exists = true;
    s.mtime = static_cast<int64_t>( t.time_since_epoch().count() );
    s.size = static_cast<uint64_t>( size );
    return s;
  }//stamp_of(...)


  /** Files below `root` (recursively), skipping names starting with a dot; `exts` lower case with a
   dot, or empty for every file.
   */
  vector<string> list_files( const string &root, const vector<string> &exts, const bool recursive )
  {
    vector<string> answer;
    std::function<void(const string &, int)> walk;
    walk = [&]( const string &dir, const int depth ){
      if( depth > ns_max_depth )
        return;

      for( const string &f : SpecUtils::ls_files_in_directory( dir ) )
      {
        if( is_dot_name( f ) )
          continue;
        if( !exts.empty() )
        {
          const string ext = SpecUtils::file_extension( f );
          bool match = false;
          for( const string &e : exts )
            match = match || SpecUtils::iequals_ascii( ext, e );
          if( !match )
            continue;
        }
        answer.push_back( f );
      }//for( loop over files )

      if( recursive )
      {
        for( const string &sub : SpecUtils::ls_directories_in_directory( dir ) )
        {
          if( !is_dot_name( sub ) )
            walk( sub, depth + 1 );
        }
      }
    };//walk

    if( SpecUtils::is_directory( root ) )
      walk( root, 0 );
    std::sort( begin(answer), end(answer) );
    return answer;
  }//list_files(...)


  string file_stem( const string &path )
  {
    const string name = SpecUtils::filename( path );
    const string ext = SpecUtils::file_extension( name );
    return name.substr( 0, name.size() - ext.size() );
  }


  /** Cuts UTF-8 text to at most `max_bytes`, at a character boundary. */
  string truncate_utf8( string text, const size_t max_bytes )
  {
    if( text.size() <= max_bytes )
      return text;
    size_t len = max_bytes;
    while( (len > 0) && ((static_cast<unsigned char>(text[len]) & 0xC0) == 0x80) )
      --len;
    text.resize( len );
    return text;
  }//truncate_utf8(...)


  // -----------------------------------------------------------------------------------------------
  //  Detector efficiency functions, through DrfLibrary (so the same files are detectors in both).
  // -----------------------------------------------------------------------------------------------

  /** Folders in the user data folder that are not searched for hand-placed GADRAS detectors/tables. */
  vector<string> user_top_excludes( const string &user_dir )
  {
    vector<string> answer;
    try
    {
      answer.push_back( UserLibrary::packagesDirectory() );
    }catch( std::exception & )
    {
    }

    for( const UserLibrary::AssetType &type : UserLibrary::assetTypes() )
    {
      vector<string> folders = type.legacyFolders;
      folders.push_back( type.folder );
      for( const string &folder : folders )
        answer.push_back( SpecUtils::append_path( user_dir, folder ) );  //includes "drfs" - its own root
    }
    return answer;
  }//user_top_excludes(...)


  DrfLibrary::Location drf_location( const string &root, const bool user_top )
  {
    DrfLibrary::Location loc;
    loc.id = loc.path = root;
    loc.isFolder = !SpecUtils::is_file( root );
    try
    {
      const string static_dir = normalize( InterSpec::staticDataDirectory() );
      loc.builtIn = !static_dir.empty() && SpecUtils::starts_with( normalize( root ), static_dir.c_str() );
    }catch( std::exception & )
    {
    }

    if( user_top )
    {
      loc.formats = DrfLibrary::FormatGadras | DrfLibrary::FormatTable;
      loc.excludePaths = user_top_excludes( root );
    }
    return loc;
  }//drf_location(...)


  DrfLibrary::Snapshot drf_scan( const DrfLibrary::Location &loc )
  {
    DrfLibrary::ScanRequest request;
    request.locations.push_back( loc );
    return DrfLibrary::scan( request );
  }


  vector<UserLibrary::FoundItem> drf_find_items( const string &root, const bool user_top )
  {
    const DrfLibrary::Snapshot snapshot = drf_scan( drf_location( root, user_top ) );

    map<string,UserLibrary::FoundItem> by_path;
    for( const DrfLibrary::Entry &e : snapshot.entries )
    {
      if( e.path.empty() )
        continue;
      UserLibrary::FoundItem &item = by_path[e.path];
      item.path = e.path;
      item.isFolder = (e.kind == DrfLibrary::EntryKind::GadrasFolder);
      if( !e.companionPath.empty()
          && (std::find( begin(item.support), end(item.support), e.companionPath ) == end(item.support)) )
        item.support.push_back( e.companionPath );
    }//for( loop over entries )

    // Detector-looking files that could not be read are still items (so they are flagged, not lost).
    for( const DrfLibrary::UnreadableItem &bad : snapshot.unreadable )
    {
      if( by_path.count( bad.path ) )
        continue;
      UserLibrary::FoundItem &item = by_path[bad.path];
      item.path = bad.path;
      item.isFolder = bad.isFolder;
    }

    vector<UserLibrary::FoundItem> answer;
    for( const pair<const string,UserLibrary::FoundItem> &p : by_path )
      answer.push_back( p.second );
    return answer;
  }//drf_find_items(...)


  /** The detector names in one item; `error` set if none could be read. */
  vector<string> drf_item_names( const UserLibrary::FoundItem &item, string &error )
  {
    DrfLibrary::Location loc = drf_location( item.path, false );
    if( item.isFolder )
      loc.formats = DrfLibrary::FormatGadras;
    const DrfLibrary::Snapshot snapshot = drf_scan( loc );

    vector<string> names;
    for( const DrfLibrary::Entry &e : snapshot.entries )
    {
      if( e.summary && (normalize( e.path ) == normalize( item.path ))
          && (std::find( begin(names), end(names), e.summary->name ) == end(names)) )
        names.push_back( e.summary->name );
    }

    if( names.empty() )
    {
      error = snapshot.unreadable.empty() ? string() : snapshot.unreadable.front().error;
      if( error.empty() )
        error = "No detector efficiency could be read.";
    }
    return names;
  }//drf_item_names(...)


  UserLibrary::Validation drf_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;

    // An InterSpec DRF export written by a newer InterSpec reads as "not a detector", so check its
    //  version first.
    if( !item.isFolder && SpecUtils::iends_with( item.path, ".xml" ) )
    {
      string contents;
      if( read_file( item.path, contents, 64*1024*1024 ) )
      {
        try
        {
          rapidxml::xml_document<char> doc;
          doc.parse<rapidxml::parse_non_destructive>( &contents[0] );
          const rapidxml::xml_node<char> *node = doc.first_node( "DetectorPeakResponse" );
          for( const rapidxml::xml_node<char> *n = doc.first_node(); !node && n; n = n->next_sibling() )
            node = n->first_node( "DetectorPeakResponse" );
          if( node )
          {
            DetectorPeakResponse drf;
            drf.fromXml( node );
          }
        }catch( XmlUtils::FormatTooNew &e )
        {
          answer.status = UserLibrary::Status::NeedsNewer;
          answer.message = e.what();
          return answer;
        }catch( std::exception & )
        {
          // Left to DrfLibrary to say why
        }
      }//if( read_file(...) )
    }//if( an XML file )

    string error;
    answer.names = drf_item_names( item, error );
    answer.status = answer.names.empty() ? UserLibrary::Status::Unreadable : UserLibrary::Status::Ok;
    answer.message = error;
    return answer;
  }//drf_validate(...)


  // -----------------------------------------------------------------------------------------------
  //  Isotopics-by-nuclides presets
  // -----------------------------------------------------------------------------------------------

  vector<UserLibrary::FoundItem> xml_file_items( const string &root )
  {
    vector<UserLibrary::FoundItem> answer;
    for( const string &f : list_files( root, { ".xml" }, true ) )
    {
      UserLibrary::FoundItem item;
      item.path = f;
      answer.push_back( item );
    }
    return answer;
  }//xml_file_items(...)


  UserLibrary::Validation rel_act_preset_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    answer.names.push_back( file_stem( item.path ) );

    string contents;
    if( !read_file( item.path, contents, 16*1024*1024 ) )
    {
      answer.message = "Could not read the file.";
      return answer;
    }

    try
    {
      rapidxml::xml_document<char> doc;
      doc.parse<rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace>( &contents[0] );
      const rapidxml::xml_node<char> *base = doc.first_node( "RelActCalcAuto" );
      if( !base )
        throw runtime_error( "Not an isotopics-by-nuclides preset (no RelActCalcAuto element)." );

#if( USE_REL_ACT_TOOL )
      RelActCalcAuto::RelActAutoGuiState state;
      state.deSerialize( base );
#endif
      answer.status = UserLibrary::Status::Ok;
    }catch( XmlUtils::FormatTooNew &e )
    {
      answer.status = UserLibrary::Status::NeedsNewer;
      answer.message = e.what();
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }

    return answer;
  }//rel_act_preset_validate(...)


  // -----------------------------------------------------------------------------------------------
  //  Activity/Shielding Fit setups ("<name>.act.shield.xml", as the tool exports them)
  // -----------------------------------------------------------------------------------------------

  /** The name a setup is listed under: its file name without ".act.shield.xml" (or ".xml"). */
  string act_shield_setup_name( const string &path )
  {
    string name = SpecUtils::filename( path );
    if( SpecUtils::iends_with( name, ".act.shield.xml" ) )
      return name.substr( 0, name.size() - 15 );
    return file_stem( path );
  }


  /** XML files whose top element is an Activity/Shielding Fit setup ("ShieldingSourceFit"). */
  vector<UserLibrary::FoundItem> act_shield_setup_items( const string &root )
  {
    vector<UserLibrary::FoundItem> answer;
    for( const string &f : list_files( root, { ".xml" }, true ) )
    {
      // A cheap look at the start of the file, rather than parsing every one.
      std::ifstream input( to_fs_path( f ), ios::in | ios::binary );
      char start[1024] = { '\0' };
      input.read( start, sizeof(start) - 1 );
      if( string( start, static_cast<size_t>( std::max<std::streamsize>( 0, input.gcount() ) ) ).find( "<ShieldingSourceFit" ) == string::npos )
        continue;

      UserLibrary::FoundItem item;
      item.path = f;
      answer.push_back( item );
    }
    return answer;
  }//act_shield_setup_items(...)


  UserLibrary::Validation act_shield_setup_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    answer.names.push_back( act_shield_setup_name( item.path ) );

    string contents;
    if( !read_file( item.path, contents, 64*1024*1024 ) )
    {
      answer.message = "Could not read the file.";
      return answer;
    }

    try
    {
      rapidxml::xml_document<char> doc;
      doc.parse<rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace>( &contents[0] );
      const rapidxml::xml_node<char> *base = doc.first_node( "ShieldingSourceFit" );
      if( !base )
        throw runtime_error( "Not an Activity/Shielding Fit setup (no ShieldingSourceFit element)." );

      ShieldingSourceDisplay::ShieldingSourceDisplayState state;
      state.deSerialize( base );
      answer.status = UserLibrary::Status::Ok;
    }catch( XmlUtils::FormatTooNew &e )
    {
      answer.status = UserLibrary::Status::NeedsNewer;
      answer.message = e.what();
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//act_shield_setup_validate(...)


  // -----------------------------------------------------------------------------------------------
  //  Reference-line sources and shielding materials: files merged by name over InterSpec's own,
  //  read when InterSpec starts.  A root may be a single file (the older add_ref_line.xml or
  //  MaterialDataBase.txt), or a folder of them.
  // -----------------------------------------------------------------------------------------------

  vector<UserLibrary::FoundItem> file_or_folder_items( const string &root, const char *ext )
  {
    vector<UserLibrary::FoundItem> answer;
    if( SpecUtils::is_file( root ) )
    {
      UserLibrary::FoundItem item;
      item.path = root;
      answer.push_back( item );
      return answer;
    }
    for( const string &f : list_files( root, { ext }, true ) )
    {
      UserLibrary::FoundItem item;
      item.path = f;
      answer.push_back( item );
    }
    return answer;
  }//file_or_folder_items(...)


  UserLibrary::Validation ref_lines_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    try
    {
      map<string,ReferenceLinePredef::NucMix> mixes;
      map<string,ReferenceLinePredef::CustomSrcLines> lines;
      ReferenceLinePredef::load_ref_line_file( item.path, mixes, lines, nullptr );
      // The loader keys sources by a normalized name; list them as they are written.
      for( const auto &m : mixes )
        answer.names.push_back( m.second.m_name.empty() ? m.first : m.second.m_name );
      for( const auto &l : lines )
        answer.names.push_back( l.second.m_name.empty() ? l.first : l.second.m_name );
      if( answer.names.empty() )
        throw runtime_error( "No reference-line sources in the file." );
      answer.status = UserLibrary::Status::Ok;
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//ref_lines_validate(...)


  UserLibrary::Validation materials_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    try
    {
      answer.names = MaterialDB::namesInFile( item.path );
      if( answer.names.empty() )
        throw runtime_error( "No materials in the file." );
      answer.status = UserLibrary::Status::Ok;
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//materials_validate(...)


  /** A check-source library (the "Source.lib" format Make Detector Response reads); names are its sources. */
  UserLibrary::Validation source_lib_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    try
    {
      {
        std::ifstream input( to_fs_path( item.path ), ios::in | ios::binary );
        if( !input || !SrcLibLineInfo::is_candidate_src_lib( input ) )
          throw runtime_error( "Not a check-source library: each line should be like"
                               " \"Cs137_ABC123 3.7E+04 01-Jan-2020 optional comments\"." );
      }
      for( const SrcLibLineInfo &src : SrcLibLineInfo::sources_in_lib( item.path ) )
        answer.names.push_back( src.m_source_name );
      if( answer.names.empty() )
        throw runtime_error( "No sources in the file." );
      answer.status = UserLibrary::Status::Ok;
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//source_lib_validate(...)


  /** A file of Nuclide Search categories (the "NuclideSearchCatagories.xml" format); names are the
   categories' (i18n keys for InterSpec's own).
   */
  UserLibrary::Validation nuc_search_category_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    try
    {
      for( const IsotopeSearchByEnergy::NucSearchCategory &cat : IsotopeSearchByEnergy::categories_in_file( item.path ) )
        answer.names.push_back( cat.key() );
      if( answer.names.empty() )
        throw runtime_error( "No categories in the file." );
      answer.status = UserLibrary::Status::Ok;
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//nuc_search_category_validate(...)


#if( USE_LLM_INTERFACE )
  /** An LLM agent definition file; names are the agents it defines (only InterSpec's agents exist). */
  UserLibrary::Validation llm_agent_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    try
    {
      for( const LlmConfig::AgentConfig &agent : LlmConfig::loadAgentsFromFile( item.path ) )
        answer.names.push_back( agent.name );
      if( answer.names.empty() )
        throw runtime_error( "No agents in the file." );
      answer.status = UserLibrary::Status::Ok;
    }catch( XmlUtils::FormatTooNew &e )
    {
      answer.status = UserLibrary::Status::NeedsNewer;
      answer.message = e.what();
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//llm_agent_validate(...)
#endif


  // -----------------------------------------------------------------------------------------------
  //  Colour themes (the JSON the Color Themes window exports)
  // -----------------------------------------------------------------------------------------------

  vector<UserLibrary::FoundItem> color_theme_items( const string &root )
  {
    vector<UserLibrary::FoundItem> answer;
    for( const string &f : list_files( root, { ".json" }, true ) )
    {
      UserLibrary::FoundItem item;
      item.path = f;
      answer.push_back( item );
    }
    return answer;
  }//color_theme_items(...)


  UserLibrary::Validation color_theme_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    string contents;
    if( !read_file( item.path, contents, 4*1024*1024 ) )
    {
      answer.message = "Could not read the file.";
      return answer;
    }

    try
    {
      ColorTheme theme;
      ColorTheme::fromJson( contents, theme );
      const string name = theme.theme_name.toUTF8();
      answer.names.push_back( name.empty() ? file_stem( item.path ) : name );
      answer.status = UserLibrary::Status::Ok;
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//color_theme_validate(...)


  // -----------------------------------------------------------------------------------------------
  //  Report templates (Inja)
  // -----------------------------------------------------------------------------------------------

  bool is_template_file_name( const string &path )
  {
    const string name = SpecUtils::filename( path );
    return SpecUtils::icontains( name, "tmplt." ) || SpecUtils::icontains( name, "template." )
           || SpecUtils::icontains( name, "inja." );
  }


  /** The files a template includes, extends, or imports (recursively), resolved in its folder. */
  vector<string> template_dependencies( const UserLibrary::FoundItem &item )
  {
    static const std::regex include_re( "\\{%-?\\s*(?:include|extends|import)\\s+[\"']([^\"']+)[\"']" );

    vector<string> answer;
    set<string> seen{ normalize( item.path ) };
    vector<string> to_check{ item.path };
    while( !to_check.empty() && (answer.size() < 100) )
    {
      const string path = to_check.back();
      to_check.pop_back();

      string contents;
      if( !read_file( path, contents, 8*1024*1024 ) )
        continue;

      const string dir = SpecUtils::parent_path( path );
      for( std::sregex_iterator it( begin(contents), end(contents), include_re ), it_end; it != it_end; ++it )
      {
        const string target = (*it)[1].str();
        if( target.empty() || SpecUtils::is_absolute_path( target ) )
          continue;
        const string full = normalize( SpecUtils::append_path( dir, target ) );
        if( !SpecUtils::is_file( full ) || !seen.insert( full ).second )
          continue;
        // Only files beside the template (or below it): a template may not pull in arbitrary files.
        if( UserLibrary::relativePath( SpecUtils::parent_path( item.path ), full ).empty() )
          continue;
        answer.push_back( full );
        to_check.push_back( full );
      }//for( loop over include statements )
    }//while( !to_check.empty() )

    return answer;
  }//template_dependencies(...)


  /** Templates are the files directly in the folder with "tmplt."/"template."/"inja." in their name
   (as the tools have always offered them); the files they include come along as support files.
   */
  vector<UserLibrary::FoundItem> template_items( const string &root )
  {
    vector<UserLibrary::FoundItem> answer;
    for( const string &f : list_files( root, {}, false ) )
    {
      if( !is_template_file_name( f ) )
        continue;
      UserLibrary::FoundItem item;
      item.path = f;
      item.support = template_dependencies( item );
      answer.push_back( item );
    }
    return answer;
  }//template_items(...)


  UserLibrary::Validation template_validate( const UserLibrary::FoundItem &item, const bool isotopics )
  {
    UserLibrary::Validation answer;
    answer.names.push_back( SpecUtils::filename( item.path ) );

#if( defined(_WIN32) )
    const char path_sep = '\\';
#else
    const char path_sep = '/';
#endif
    const string dir = SpecUtils::parent_path( item.path ) + path_sep;

    try
    {
#if( USE_REL_ACT_TOOL )
      inja::Environment env = isotopics ? RelActAutoReport::get_default_inja_env( dir )
                                        : BatchInfoLog::get_default_inja_env( dir );
#else
      inja::Environment env = BatchInfoLog::get_default_inja_env( dir );
#endif
      // A template in a package is rendered with includes only from its own folder; check it that way.
      bool in_packages = false;
      try
      {
        in_packages = UserLibrary::available()
                      && UserLibrary::isInPackages( item.path );
      }catch( std::exception & )
      {
      }

      if( in_packages )
      {
        const vector<string> outside = UserLibrary::templateIncludesOutsideFolder( item.path );
        if( !outside.empty() )
          throw runtime_error( "Includes '" + outside.front() + "', which is not beside it." );
        UserLibrary::confineTemplateIncludes( env, item.path );
      }

      string text;
      if( !read_file( item.path, text, 16*1024*1024 ) )
        throw runtime_error( "Could not read the file." );
      env.parse( text );
      answer.status = UserLibrary::Status::Ok;
    }catch( std::exception &e )
    {
      answer.status = UserLibrary::Status::Unreadable;
      answer.message = e.what();
    }
    return answer;
  }//template_validate(...)


  // -----------------------------------------------------------------------------------------------
  //  Reference spectra: each folder directly in `reference_spectra/` (a set, e.g. one detector's)
  //  is one item, as is each spectrum file directly in it.
  // -----------------------------------------------------------------------------------------------

  bool is_readme( const string &path )
  {
    const string name = SpecUtils::filename( path );
    return SpecUtils::iequals_ascii( name, "readme.txt" ) || SpecUtils::iequals_ascii( name, "readme.xml" );
  }


  vector<UserLibrary::FoundItem> ref_spectra_items( const string &root )
  {
    vector<UserLibrary::FoundItem> answer;
    if( !SpecUtils::is_directory( root ) )
      return answer;

    vector<string> dirs = SpecUtils::ls_directories_in_directory( root );
    std::sort( begin(dirs), end(dirs) );
    for( const string &d : dirs )
    {
      if( is_dot_name( d ) || list_files( d, {}, true ).empty() )
        continue;
      UserLibrary::FoundItem item;
      item.path = d;
      item.isFolder = true;
      answer.push_back( item );
    }

    for( const string &f : list_files( root, {}, false ) )
    {
      if( is_readme( f ) )
        continue;
      UserLibrary::FoundItem item;
      item.path = f;
      answer.push_back( item );
    }
    return answer;
  }//ref_spectra_items(...)


  UserLibrary::Validation ref_spectra_validate( const UserLibrary::FoundItem &item )
  {
    UserLibrary::Validation answer;
    answer.names.push_back( SpecUtils::filename( item.path ) );

    vector<string> files;
    if( item.isFolder )
    {
      for( const string &f : list_files( item.path, {}, true ) )
        if( !is_readme( f ) )
          files.push_back( f );
    }else
    {
      files.push_back( item.path );
    }

    // A set is readable if a spectrum in it reads (trying a few - reading them all can be slow).
    for( size_t i = 0; (i < files.size()) && (i < 5); ++i )
    {
      SpecUtils::SpecFile spec;
      if( spec.load_file( files[i], SpecUtils::ParserType::Auto, files[i] ) )
      {
        answer.status = UserLibrary::Status::Ok;
        return answer;
      }
    }

    answer.status = UserLibrary::Status::Unreadable;
    answer.message = files.empty() ? "No files." : "No spectrum file could be read.";
    return answer;
  }//ref_spectra_validate(...)


  // -----------------------------------------------------------------------------------------------
  //  LLM knowledge skills: a folder with a SKILL.md (and optional `references/`).
  // -----------------------------------------------------------------------------------------------

  vector<UserLibrary::FoundItem> skill_items( const string &root )
  {
    vector<string> folders;
    for( const string &f : list_files( root, { ".md" }, true ) )
    {
      if( SpecUtils::filename( f ) != "SKILL.md" )
        continue;
      const string folder = SpecUtils::parent_path( f );
      if( normalize( folder ) != normalize( root ) )
        folders.push_back( folder );
    }
    std::sort( begin(folders), end(folders) );

    // A SKILL.md inside another skill's folder is part of that skill.
    vector<UserLibrary::FoundItem> answer;
    for( const string &folder : folders )
    {
      bool nested = false;
      for( const UserLibrary::FoundItem &outer : answer )
        nested = nested || !UserLibrary::relativePath( outer.path, folder ).empty();
      if( nested )
        continue;

      UserLibrary::FoundItem item;
      item.path = folder;
      item.isFolder = true;
      answer.push_back( item );
    }
    return answer;
  }//skill_items(...)


  /** The canonical skill names in a skill folder; `error` set if there are none. */
  vector<string> skill_names( const UserLibrary::FoundItem &item, string &error )
  {
    vector<string> answer;
#if( USE_LLM_INTERFACE )
    for( const LlmKnowledgeBase::SkillInfo &skill : LlmKnowledgeBase::loadSkillsFromKnowledgeDir( item.path ) )
      answer.push_back( skill.canonical_name );
#else
    // Without the LLM interface the name only labels the item: the front-matter `name:`, or the folder.
    string contents;
    if( read_file( SpecUtils::append_path( item.path, "SKILL.md" ), contents, 4*1024*1024 ) )
    {
      std::smatch m;
      static const std::regex name_re( "(?:^|\\n)name:\\s*([^\\r\\n]+)" );
      string name = std::regex_search( contents, m, name_re ) ? SpecUtils::trim_copy( m[1].str() )
                                                              : SpecUtils::filename( item.path );
      SpecUtils::to_lower_ascii( name );
      answer.push_back( name );
    }
#endif
    if( answer.empty() )
      error = "No readable SKILL.md (it needs a front-matter block with a name and description).";
    return answer;
  }//skill_names(...)


  // -----------------------------------------------------------------------------------------------
  //  Caches and the history file
  // -----------------------------------------------------------------------------------------------

  std::atomic<uint64_t> ns_generation( 1 );

  /** Guards the caches below. */
  std::mutex ns_cache_mutex;

  struct HistoryCache
  {
    string path;
    Stamp stamp;
    shared_ptr<const json> value;  ///< shared, not copied: name lookups read it for every name
    bool loaded = false;
  };
  HistoryCache ns_history_cache;

  struct ManifestCache
  {
    Stamp stamp;
    shared_ptr<const UserLibrary::Manifest> manifest;
  };
  map<string,ManifestCache> ns_manifest_cache;

  struct ValidationCache
  {
    vector<Stamp> stamps;
    UserLibrary::Validation validation;
  };
  map<string,ValidationCache> ns_validation_cache;  //key: "<type>|<path>"

  /** Serializes #modifyHistory within the process. */
  std::recursive_mutex ns_modify_mutex;


  /** Reads and parses the history file; an empty object if missing.  Throws on invalid JSON. */
  json read_history_file( const string &path )
  {
    string contents;
    if( !SpecUtils::is_file( path ) )
      return json::object();
    if( !read_file( path, contents, ns_max_json_bytes ) )
      throw runtime_error( "Could not read '" + path + "'" );
    json value = json::parse( contents );
    if( !value.is_object() )
      throw runtime_error( "'" + path + "' is not a JSON object" );
    return value;
  }//read_history_file(...)


  /** The history file's contents, re-read only when the file changes; never null. */
  shared_ptr<const json> cached_history()
  {
    if( !UserLibrary::available() )
      return make_shared<const json>( json::object() );

    const string path = UserLibrary::historyPath();
    const Stamp stamp = stamp_of( path );

    std::lock_guard<std::mutex> lock( ns_cache_mutex );
    if( ns_history_cache.loaded && ns_history_cache.value
        && (ns_history_cache.path == path) && (ns_history_cache.stamp == stamp) )
      return ns_history_cache.value;

    json value = json::object();
    try
    {
      value = read_history_file( path );
    }catch( std::exception &e )
    {
      cerr << "UserLibrary: could not read history file: " << e.what() << endl;
    }

    ns_history_cache.path = path;
    ns_history_cache.stamp = stamp;
    ns_history_cache.value = make_shared<const json>( std::move( value ) );
    ns_history_cache.loaded = true;
    return ns_history_cache.value;
  }//cached_history()


  shared_ptr<const UserLibrary::Manifest> cached_manifest( const string &path )
  {
    const Stamp stamp = stamp_of( path );
    if( !stamp.exists )
      return nullptr;

    {
      std::lock_guard<std::mutex> lock( ns_cache_mutex );
      const map<string,ManifestCache>::const_iterator pos = ns_manifest_cache.find( path );
      if( (pos != end(ns_manifest_cache)) && (pos->second.stamp == stamp) )
        return pos->second.manifest;
    }

    shared_ptr<const UserLibrary::Manifest> manifest;
    string contents;
    try
    {
      if( !read_file( path, contents, ns_max_json_bytes ) )
        throw runtime_error( "could not read" );
      manifest = make_shared<const UserLibrary::Manifest>( UserLibrary::Manifest::fromJson( contents ) );
    }catch( std::exception &e )
    {
      cerr << "UserLibrary: invalid manifest '" << path << "': " << e.what() << endl;
    }

    std::lock_guard<std::mutex> lock( ns_cache_mutex );
    ns_manifest_cache[path] = ManifestCache{ stamp, manifest };
    return manifest;
  }//cached_manifest(...)


  string json_str( const json &j, const char *key )
  {
    const json::const_iterator pos = j.find( key );
    if( pos == j.end() )
      return string();
    if( pos->is_string() )
      return pos->get<string>();
    if( pos->is_number() )
      return pos->dump();
    return string();
  }//json_str(...)


  /** The disabled entries (item_id or path) of a package, from the history. */
  set<string> package_disabled( const json &history, const string &package_id )
  {
    set<string> answer;
    const json::const_iterator packages = history.find( "packages" );
    if( (packages == history.end()) || !packages->is_object() )
      return answer;
    const json::const_iterator pkg = packages->find( package_id );
    if( (pkg == packages->end()) || !pkg->is_object() )
      return answer;
    const json::const_iterator disabled = pkg->find( "disabled" );
    if( (disabled != pkg->end()) && disabled->is_array() )
    {
      for( const json &d : *disabled )
        if( d.is_string() )
          answer.insert( d.get<string>() );
    }
    return answer;
  }//package_disabled(...)


  /** The user's own disabled items, as (type, relative path). */
  set<pair<string,string>> user_disabled( const json &history )
  {
    set<pair<string,string>> answer;
    const json::const_iterator disabled = history.find( "disabled_user" );
    if( (disabled == history.end()) || !disabled->is_array() )
      return answer;
    for( const json &d : *disabled )
    {
      if( d.is_object() )
        answer.insert( { json_str( d, "type" ), json_str( d, "path" ) } );
    }
    return answer;
  }//user_disabled(...)


  /** The InterSpec docroot folder (holding `InterSpec_resources`): the static data folder's parent. */
  string app_root_dir()
  {
    const string static_dir = InterSpec::staticDataDirectory().empty() ? string("./data")
                                                                       : InterSpec::staticDataDirectory();
    return normalize( SpecUtils::append_path( static_dir, ".." ) );
  }
  /** UserLibrary::resolve, with the recorded choices passed in (so listing many names reads them once). */
  std::optional<UserLibrary::Item> resolve_among( const std::string &typeId, const std::string &name,
                                                  const std::vector<UserLibrary::Item> &candidates,
                                                  const std::map<std::string,UserLibrary::Choice> &choices )
  {
    const UserLibrary::AssetType * const type = UserLibrary::assetType( typeId );
    if( !type )
      return std::nullopt;

    vector<const UserLibrary::Item *> matching;
    for( const UserLibrary::Item &item : candidates )
    {
      if( item.disabled || item.needsNewer || (item.typeId != typeId) )
        continue;
      for( const string &n : item.names )
      {
        if( UserLibrary::sameName( *type, n, name ) )
        {
          matching.push_back( &item );
          break;
        }
      }
    }//for( loop over candidates )

    if( matching.empty() )
      return std::nullopt;

    const string prefix = typeId + ":";
    for( const pair<const string,UserLibrary::Choice> &c : choices )
    {
      if( !SpecUtils::starts_with( c.first, prefix.c_str() )
          || !UserLibrary::sameName( *type, c.first.substr( prefix.size() ), name ) )
        continue;

      for( const UserLibrary::Item *item : matching )
      {
        if( item->source.key() == c.second.source )
          return *item;
      }
    }//for( loop over choices )

    return *matching.front();
  }//resolve_among(...)
}//namespace


namespace UserLibrary
{
const char * const ManifestFileName = "interspec-package.json";


bool available()
{
#if( USE_USER_LIBRARY )
  try
  {
    const string user_dir = InterSpec::writableDataDirectory();
    if( user_dir.empty() )
      return false;

    string static_dir;
    try
    {
      static_dir = InterSpec::staticDataDirectory();
    }catch( std::exception & )
    {
    }

    return static_dir.empty() || (normalize( user_dir ) != normalize( static_dir ));
  }catch( std::exception & )
  {
  }
#endif
  return false;
}//bool available()


string packagesDirectory()
{
#if( USE_USER_LIBRARY )
  if( !available() )
    throw runtime_error( "No user library in this process." );
  return SpecUtils::append_path( InterSpec::writableDataDirectory(), "packages" );
#else
  throw runtime_error( "No user library in this build." );
#endif
}//string packagesDirectory()


string historyPath()
{
#if( USE_USER_LIBRARY )
  if( !available() )
    throw runtime_error( "No user library in this process." );
  return SpecUtils::append_path( InterSpec::writableDataDirectory(), "user_library.json" );
#else
  throw runtime_error( "No user library in this build." );
#endif
}//string historyPath()


const std::vector<AssetType> &assetTypes()
{
  static const vector<AssetType> types = [](){
    vector<AssetType> answer;

    {
      AssetType t;
      t.id = "drf";
      t.titleKey = "ul-type-drf";
      t.folder = "drfs";
      t.builtInPaths = { "GenericGadrasDetectors", "common_drfs.tsv" };
      t.kind = Kind::List;
      t.findItems = []( const string &root ){ return drf_find_items( root, false ); };
      t.validate = &drf_validate;
      t.names = []( const FoundItem &item ){ string error; return drf_item_names( item, error ); };
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "rel_act_preset";
      t.titleKey = "ul-type-rel-act-preset";
      t.folder = "rel_act";
      t.builtInPaths = { "rel_act" };
      t.caseInsensitiveNames = true;
      t.findItems = &xml_file_items;
      t.validate = &rel_act_preset_validate;
      t.names = []( const FoundItem &item ){ return vector<string>{ file_stem( item.path ) }; };
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "act_shield_setup";
      t.titleKey = "ul-type-act-shield-setup";
      t.folder = "act_shield_setups";
      t.caseInsensitiveNames = true;
      t.findItems = &act_shield_setup_items;
      t.validate = &act_shield_setup_validate;
      t.names = []( const FoundItem &item ){ return vector<string>{ act_shield_setup_name( item.path ) }; };
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "act_shield_template";
      t.titleKey = "ul-type-act-shield-template";
      t.folder = "report_templates/act_shield";
      t.legacyFolders = { "act_shield_fit_templates" };
      t.caseInsensitiveNames = true;
      t.findItems = &template_items;
      t.validate = []( const FoundItem &item ){ return template_validate( item, false ); };
      t.names = []( const FoundItem &item ){ return vector<string>{ SpecUtils::filename( item.path ) }; };
      t.dependencies = &template_dependencies;
      answer.push_back( t );
    }

    // Batch analysis report templates: what each renders differs (the JSON for one file, or a summary
    //  of all files, of a peak fit, activity/shielding fit, or isotopics), so each kind has its folder.
    //  Per-file Activity/Shielding reports use "act_shield_template" (the same JSON as the tool's log).
    for( const char * const kind : { "peak_fit_template", "peak_fit_summary_template", "act_shield_summary_template" } )
    {
      AssetType t;
      t.id = kind;
      t.caseInsensitiveNames = true;
      t.findItems = &template_items;
      t.validate = []( const FoundItem &item ){ return template_validate( item, false ); };
      t.names = []( const FoundItem &item ){ return vector<string>{ SpecUtils::filename( item.path ) }; };
      t.dependencies = &template_dependencies;
      if( t.id == "peak_fit_template" )
      {
        t.titleKey = "ul-type-peak-fit-template";
        t.folder = "report_templates/peak_fit";
      }else if( t.id == "peak_fit_summary_template" )
      {
        t.titleKey = "ul-type-peak-fit-summary-template";
        t.folder = "report_templates/peak_fit_summary";
      }else
      {
        t.titleKey = "ul-type-act-shield-summary-template";
        t.folder = "report_templates/act_shield_summary";
      }
      answer.push_back( t );
    }//for( loop over batch template kinds )

    {
      AssetType t;
      t.id = "isotopics_template";
      t.titleKey = "ul-type-isotopics-template";
      t.folder = "report_templates/isotopics_auto";
      t.legacyFolders = { "IsotopicsByNuclidesReportTmplts" };
      t.builtInPaths = { "InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts" };
      t.caseInsensitiveNames = true;
      t.findItems = &template_items;
      t.validate = []( const FoundItem &item ){ return template_validate( item, true ); };
      t.names = []( const FoundItem &item ){ return vector<string>{ SpecUtils::filename( item.path ) }; };
      t.dependencies = &template_dependencies;
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "isotopics_summary_template";
      t.titleKey = "ul-type-isotopics-summary-template";
      t.folder = "report_templates/isotopics_auto_summary";
      t.caseInsensitiveNames = true;
      t.findItems = &template_items;
      t.validate = []( const FoundItem &item ){ return template_validate( item, true ); };
      t.names = []( const FoundItem &item ){ return vector<string>{ SpecUtils::filename( item.path ) }; };
      t.dependencies = &template_dependencies;
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "ref_lines";
      t.titleKey = "ul-type-ref-lines";
      // As the reference-line loader matches them: ignoring case and spaces.
      t.nameKey = []( const string &name ){
        string key = SpecUtils::trim_copy( name );
        SpecUtils::ireplace_all( key, " ", "" );
        return SpecUtils::to_lower_ascii_copy( key );
      };
      t.folder = "ref_lines";
      t.legacyUserFiles = { "add_ref_line.xml" };
      t.builtInPaths = { "add_ref_line.xml" };
      t.kind = Kind::MergeByName;
      t.apply = ApplyMode::Restart;
      t.findItems = []( const string &root ){ return file_or_folder_items( root, ".xml" ); };
      t.validate = &ref_lines_validate;
      t.names = []( const FoundItem &item ){ return ref_lines_validate( item ).names; };
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "materials";
      t.titleKey = "ul-type-materials";
      t.folder = "materials";
      // As the material database matches them: ignoring case, spaces, '-', and '_'.
      t.nameKey = []( const string &name ){
        string key = SpecUtils::to_lower_ascii_copy( name );
        SpecUtils::erase_any_character( key, "_- \t" );
        return key;
      };
      t.legacyUserFiles = { "MaterialDataBase.txt" };
      t.builtInPaths = { "MaterialDataBase.txt" };
      t.kind = Kind::MergeByName;
      t.apply = ApplyMode::Restart;
      t.findItems = []( const string &root ){ return file_or_folder_items( root, ".txt" ); };
      t.validate = &materials_validate;
      t.names = []( const FoundItem &item ){ return materials_validate( item ).names; };
      answer.push_back( t );
    }

    {
      // The Nuclide Search tool's categories; re-read when the library changes (in open sessions too).
      AssetType t;
      t.id = "nuc_search_category";
      t.titleKey = "ul-type-nuc-search-category";
      t.folder = "nuclide_search_categories";
      t.legacyUserFiles = { "NuclideSearchCatagories.xml" };
      t.builtInPaths = { "NuclideSearchCatagories.xml" };
      t.kind = Kind::MergeByName;
      t.findItems = []( const string &root ){ return file_or_folder_items( root, ".xml" ); };
      t.validate = &nuc_search_category_validate;
      t.names = []( const FoundItem &item ){ return nuc_search_category_validate( item ).names; };
      answer.push_back( t );
    }

    {
      // Make Detector Response's check-source libraries.  It also still reads "Source.lib" files from
      //  anywhere else in the user data folder (see SrcLibLineInfo::sources_in_all_libs).
      AssetType t;
      t.id = "source_lib";
      t.titleKey = "ul-type-source-lib";
      t.folder = "source_libs";
      t.legacyUserFiles = { "Source.lib" };
      t.findItems = []( const string &root ){ return file_or_folder_items( root, ".lib" ); };
      t.validate = &source_lib_validate;
      t.names = []( const FoundItem &item ){ return source_lib_validate( item ).names; };
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "color_theme";
      t.titleKey = "ul-type-color-theme";
      t.folder = "color_themes";
      t.findItems = &color_theme_items;
      t.validate = &color_theme_validate;
      t.names = []( const FoundItem &item ){
        const Validation v = color_theme_validate( item );
        return v.names.empty() ? vector<string>{ file_stem( item.path ) } : v.names;
      };
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "ref_spectra";
      t.titleKey = "ul-type-ref-spectra";
      t.folder = "reference_spectra";
      t.builtInPaths = { "reference_spectra" };
      t.caseInsensitiveNames = true;
      t.findItems = &ref_spectra_items;
      t.validate = &ref_spectra_validate;
      t.names = []( const FoundItem &item ){ return vector<string>{ SpecUtils::filename( item.path ) }; };
      answer.push_back( t );
    }

    {
      AssetType t;
      t.id = "llm_skill";
      t.titleKey = "ul-type-llm-skill";
      t.folder = "llm_knowledge";
      t.legacyFolders = { "gamma_spec_skills/llm_knowledge" };
      // The skills repository checkout first: it overrides same-named skills in "llm_knowledge" (as
      //  InterSpec has always done).
      t.builtInPaths = { "gamma_spec_skills/llm_knowledge", "llm_knowledge" };
      t.kind = Kind::MergeByName;
      t.sensitive = true;
      t.findItems = &skill_items;
      t.validate = []( const FoundItem &item ){
        Validation v;
        v.names = skill_names( item, v.message );
        v.status = v.names.empty() ? Status::Unreadable : Status::Ok;
        return v;
      };
      t.names = []( const FoundItem &item ){ string error; return skill_names( item, error ); };
      answer.push_back( t );
    }

#if( USE_LLM_INTERFACE )
    {
      // The AI assistant's agents: each agent's definition comes from one file (LlmConfig::agentFiles).
      AssetType t;
      t.id = "llm_agent";
      t.titleKey = "ul-type-llm-agent";
      t.folder = "llm_agents";
      for( const AgentType type : definedAgentTypes() )
      {
        const string filename = "llm_agent_" + agentTypeToString( type ) + ".xml";
        t.legacyUserFiles.push_back( filename );
        t.builtInPaths.push_back( filename );
      }
      t.kind = Kind::MergeByName;
      t.apply = ApplyMode::Restart;
      t.sensitive = true;
      t.findItems = []( const string &root ){ return file_or_folder_items( root, ".xml" ); };
      t.validate = &llm_agent_validate;
      t.names = []( const FoundItem &item ){ return llm_agent_validate( item ).names; };
      answer.push_back( t );
    }
#endif

    return answer;
  }();

  return types;
}//const std::vector<AssetType> &assetTypes()


const AssetType *assetType( const std::string &typeId )
{
  for( const AssetType &t : assetTypes() )
  {
    if( t.id == typeId )
      return &t;
  }
  return nullptr;
}//assetType(...)


const AssetType *typeForRelativePath( const std::string &relPath, std::string &folder )
{
  const AssetType *answer = nullptr;
  folder.clear();
  for( const AssetType &t : assetTypes() )
  {
    vector<string> folders = t.legacyFolders;
    folders.push_back( t.folder );
    for( const string &f : folders )
    {
      // Longest match wins ("report_templates/act_shield" over a hypothetical "report_templates").
      if( (f.size() > folder.size()) && (relPath.size() > f.size())
          && SpecUtils::istarts_with( relPath, f ) && (relPath[f.size()] == '/') )
      {
        answer = &t;
        folder = relPath.substr( 0, f.size() );
      }
    }//for( loop over folders )
  }//for( loop over types )
  return answer;
}//typeForRelativePath(...)


std::string Source::key() const
{
  switch( kind )
  {
    case SourceKind::BuiltIn: return "builtin";
    case SourceKind::User:    return "user";
    case SourceKind::Package: return "package:" + packageId;
  }
  return "user";
}//Source::key()


bool Source::operator==( const Source &rhs ) const
{
  return (kind == rhs.kind) && (packageId == rhs.packageId);
}


std::vector<Root> searchRoots( const std::string &typeId )
{
  vector<Root> answer;
  const AssetType * const type = assetType( typeId );
  if( !type )
    return answer;

  vector<string> folders{ type->folder };
  folders.insert( end(folders), begin(type->legacyFolders), end(type->legacyFolders) );

#if( USE_USER_LIBRARY )
  if( available() )
  {
    const string user_dir = InterSpec::writableDataDirectory();
    for( const string &folder : folders )
    {
      const string path = SpecUtils::append_path( user_dir, folder );
      if( SpecUtils::is_directory( path ) )
      {
        Root root;
        root.source.kind = SourceKind::User;
        root.path = path;
        root.base = user_dir;
        answer.push_back( root );
      }
    }//for( loop over folders )

    if( type->id == "drf" )
    {
      Root root;
      root.source.kind = SourceKind::User;
      root.path = root.base = user_dir;
      root.userDataTop = true;
      answer.push_back( root );
    }

    for( const string &file : type->legacyUserFiles )
    {
      const string path = SpecUtils::append_path( user_dir, file );
      if( SpecUtils::is_file( path ) )
      {
        Root root;
        root.source.kind = SourceKind::User;
        root.path = path;
        root.base = user_dir;
        answer.push_back( root );
      }
    }//for( loop over older single files )

    for( const Package &pkg : installedPackages() )
    {
      for( const string &folder : folders )
      {
        const string path = SpecUtils::append_path( pkg.path, folder );
        if( SpecUtils::is_directory( path ) )
        {
          Root root;
          root.source.kind = SourceKind::Package;
          root.source.packageId = pkg.id;
          root.source.label = pkg.manifest.name;
          root.path = path;
          root.base = pkg.path;
          answer.push_back( root );
        }
      }//for( loop over folders )
    }//for( loop over packages )
  }//if( available() )
#endif

  string static_dir;
  try
  {
    static_dir = InterSpec::staticDataDirectory();
  }catch( std::exception & )
  {
  }

  if( !static_dir.empty() )
  {
    for( const string &builtin : type->builtInPaths )
    {
      const bool in_docroot = SpecUtils::starts_with( builtin, "InterSpec_resources/" );
      const string base = in_docroot ? app_root_dir() : static_dir;
      const string path = SpecUtils::append_path( base, builtin );
      if( SpecUtils::is_directory( path ) || SpecUtils::is_file( path ) )
      {
        Root root;
        root.source.kind = SourceKind::BuiltIn;
        root.path = path;
        root.base = base;
        answer.push_back( root );
      }
    }//for( loop over built-in paths )
  }//if( !static_dir.empty() )

  return answer;
}//searchRoots(...)


Validation validateItem( const std::string &typeId, const FoundItem &item )
{
  const AssetType * const type = assetType( typeId );
  if( !type || !type->validate )
  {
    Validation v;
    v.message = "Unknown type '" + typeId + "'.";
    return v;
  }

  // A folder's own time changes when files are added or removed, and package folders are replaced
  //  as a whole rather than edited, so it stands in for the files in it.
  vector<Stamp> stamps{ stamp_of( item.path ) };
  for( const string &s : item.support )
    stamps.push_back( stamp_of( s ) );

  const string key = typeId + "|" + normalize( item.path );
  {
    std::lock_guard<std::mutex> lock( ns_cache_mutex );
    const map<string,ValidationCache>::const_iterator pos = ns_validation_cache.find( key );
    if( (pos != end(ns_validation_cache)) && (pos->second.stamps == stamps) )
      return pos->second.validation;
  }

  const Validation answer = type->validate( item );

  std::lock_guard<std::mutex> lock( ns_cache_mutex );
  ns_validation_cache[key] = ValidationCache{ stamps, answer };
  return answer;
}//validateItem(...)


std::string nameKey( const AssetType &type, const std::string &name )
{
  if( type.nameKey )
    return type.nameKey( name );
  return type.caseInsensitiveNames ? SpecUtils::to_lower_ascii_copy( name ) : name;
}


bool sameName( const AssetType &type, const std::string &lhs, const std::string &rhs )
{
  return nameKey( type, lhs ) == nameKey( type, rhs );
}


std::vector<Item> items( const std::string &typeId, const bool includeDisabled )
{
  vector<Item> answer;
  const AssetType * const type = assetType( typeId );
  if( !type || !type->findItems )
    return answer;

  const shared_ptr<const json> history_owner = cached_history();
  const json &history = *history_owner;
  const set<pair<string,string>> disabled_user = user_disabled( history );
  const vector<Package> packages = available() ? installedPackages() : vector<Package>();

  set<string> seen_paths;
  for( const Root &root : searchRoots( typeId ) )
  {
    const Package *pkg = nullptr;
    for( const Package &p : packages )
      pkg = (!pkg && (root.source.kind == SourceKind::Package) && (p.id == root.source.packageId)) ? &p : pkg;

    const vector<FoundItem> found = (type->id == "drf") ? drf_find_items( root.path, root.userDataTop )
                                                        : type->findItems( root.path );
    for( const FoundItem &f : found )
    {
      if( !seen_paths.insert( normalize( f.path ) ).second )
        continue;

      Item item;
      item.typeId = typeId;
      item.source = root.source;
      item.found = f;
      item.relPath = relativePath( root.base, f.path );
      if( item.relPath.empty() )
        item.relPath = SpecUtils::filename( f.path );

      if( pkg )
      {
        const ManifestItem * const mi = pkg->manifest.itemForPath( item.relPath );
        if( mi )
          item.itemId = mi->itemId;
        item.disabled = (!item.itemId.empty() && pkg->disabled.count( item.itemId ))
                        || pkg->disabled.count( item.relPath );

        // Written for a newer InterSpec: installed, but unused until InterSpec is upgraded.
        item.needsNewer = pkg->manifest.needsNewerInterSpec()
                          || (validateItem( typeId, f ).status == Status::NeedsNewer);
      }else if( root.source.kind == SourceKind::User )
      {
        item.disabled = disabled_user.count( { typeId, item.relPath } ) > 0;
      }

      if( (item.disabled || item.needsNewer) && !includeDisabled )
        continue;

      if( type->names )
        item.names = type->names( f );
      answer.push_back( item );
    }//for( loop over found items )
  }//for( loop over roots )

  // Which items answer lookups by name.
  const map<string,Choice> current_choices = choices();
  for( Item &item : answer )
  {
    if( item.disabled || item.needsNewer || item.names.empty() )
    {
      item.inUse = !item.disabled && !item.needsNewer;
      continue;
    }

    item.inUse = false;
    for( const string &name : item.names )
    {
      const optional<Item> chosen = resolve_among( typeId, name, answer, current_choices );
      item.inUse = item.inUse || (chosen && (normalize( chosen->found.path ) == normalize( item.found.path )));
    }
  }//for( Item &item : answer )

  return answer;
}//items(...)


std::optional<Item> resolve( const std::string &typeId, const std::string &name )
{
  return resolve( typeId, name, items( typeId, false ) );
}


std::optional<Item> resolve( const std::string &typeId, const std::string &name,
                             const std::vector<Item> &candidates )
{
  return resolve_among( typeId, name, candidates, choices() );
}//resolve(...)


// -------------------------------------------------------------------------------------------------
//  Manifest
// -------------------------------------------------------------------------------------------------

Manifest Manifest::fromJson( const std::string &text )
{
  return fromJsonObject( json::parse( text ) );
}


Manifest Manifest::fromJsonObject( const nlohmann::json &j )
{
  if( !j.is_object() )
    throw runtime_error( "The manifest is not a JSON object." );

  Manifest m;
  const json::const_iterator fv = j.find( "format_version" );
  if( (fv == j.end()) || !fv->is_number_integer() )
    throw runtime_error( "The manifest has no format_version." );
  m.formatVersion = fv->get<int>();

  // Untrusted text: keep each field to a sane length.
  const size_t short_len = 512, long_len = 64*1024;
  m.id = truncate_utf8( json_str( j, "id" ), 128 );
  m.name = truncate_utf8( json_str( j, "name" ), short_len );
  m.version = truncate_utf8( json_str( j, "version" ), 64 );
  if( m.id.empty() )
    throw runtime_error( "The manifest has no id." );
  if( m.name.empty() )
    throw runtime_error( "The manifest has no name." );
  if( m.version.empty() )
    throw runtime_error( "The manifest has no version." );

  m.created = truncate_utf8( json_str( j, "created" ), 64 );
  m.createdWith = truncate_utf8( json_str( j, "created_with" ), short_len );
  m.minInterSpecVersion = truncate_utf8( json_str( j, "min_interspec_version" ), 64 );
  m.author = truncate_utf8( json_str( j, "author" ), short_len );
  m.organization = truncate_utf8( json_str( j, "organization" ), short_len );
  m.contact = truncate_utf8( json_str( j, "contact" ), short_len );
  m.description = truncate_utf8( json_str( j, "description" ), 2*short_len );
  m.notes = truncate_utf8( json_str( j, "notes" ), long_len );
  m.changes = truncate_utf8( json_str( j, "changes" ), long_len );
  m.handling = truncate_utf8( json_str( j, "handling" ), short_len );
  m.reviewBy = truncate_utf8( json_str( j, "review_by" ), 64 );

  // Bounded, as manifests are untrusted and are kept in the history file.
  const size_t max_previous_versions = 200, max_items = 10000;

  const json::const_iterator prev = j.find( "previous_versions" );
  if( (prev != j.end()) && prev->is_array() )
  {
    for( const json &p : *prev )
    {
      if( !p.is_object() )
        continue;
      if( m.previousVersions.size() >= max_previous_versions )
        break;
      ManifestPreviousVersion v;
      v.version = truncate_utf8( json_str( p, "version" ), 64 );
      v.created = truncate_utf8( json_str( p, "created" ), 64 );
      v.changes = truncate_utf8( json_str( p, "changes" ), long_len );
      m.previousVersions.push_back( v );
    }
  }//if( previous_versions )

  const json::const_iterator items = j.find( "items" );
  if( (items != j.end()) && items->is_array() )
  {
    if( items->size() > max_items )
      throw runtime_error( "The manifest lists too many items." );
    for( const json &i : *items )
    {
      if( !i.is_object() )
        throw runtime_error( "A manifest item is not a JSON object." );

      ManifestItem item;
      item.itemId = truncate_utf8( json_str( i, "item_id" ), 128 );
      item.path = json_str( i, "path" );
      SpecUtils::ireplace_all( item.path, "\\", "/" );
      while( !item.path.empty() && (item.path.back() == '/') )
        item.path.pop_back();
      item.type = json_str( i, "type" );
      if( item.path.empty() || item.type.empty() )
        throw runtime_error( "A manifest item has no path or type." );

      // A path inside the package: relative, without "." or ".." components or a drive letter.
      vector<string> components;
      SpecUtils::split( components, item.path, "/" );
      bool path_ok = (item.path.front() != '/') && (item.path.find( ':' ) == string::npos)
                     && (item.path.find( "//" ) == string::npos)
                     && (item.path.find( '\0' ) == string::npos) && (item.path.size() <= 1024);
      for( const string &c : components )
        path_ok = path_ok && !c.empty() && (c != ".") && (c != "..");
      if( !path_ok )
        throw runtime_error( "A manifest item's path is not a path inside the package." );
      item.sha256 = json_str( i, "sha256" );
      item.notes = truncate_utf8( json_str( i, "notes" ), long_len );

      const json::const_iterator derived = i.find( "derived_from" );
      if( (derived != i.end()) && derived->is_object() )
      {
        item.derivedFromName = truncate_utf8( json_str( *derived, "name" ), short_len );
        item.derivedFromId = truncate_utf8( json_str( *derived, "id" ), 128 );
        item.derivedFromVersion = truncate_utf8( json_str( *derived, "version" ), 64 );
      }

      const json::const_iterator override_pos = i.find( "intended_override" );
      item.intendedOverride = (override_pos != i.end()) && override_pos->is_boolean() && override_pos->get<bool>();
      m.items.push_back( item );
    }//for( loop over items )
  }//if( items )

  return m;
}//Manifest::fromJsonObject(...)


nlohmann::json Manifest::toJsonObject() const
{
  json j = json::object();
  j["format_version"] = formatVersion;
  j["id"] = id;
  j["name"] = name;
  j["version"] = version;
  j["created"] = created;
  j["created_with"] = createdWith;
  j["min_interspec_version"] = minInterSpecVersion;
  j["author"] = author;
  j["organization"] = organization;
  j["contact"] = contact;
  j["description"] = description;
  j["notes"] = notes;
  j["changes"] = changes;
  j["handling"] = handling;
  j["review_by"] = reviewBy;

  if( !previousVersions.empty() )
  {
    json prev = json::array();
    for( const ManifestPreviousVersion &v : previousVersions )
      prev.push_back( { {"version", v.version}, {"created", v.created}, {"changes", v.changes} } );
    j["previous_versions"] = prev;
  }

  json item_array = json::array();
  for( const ManifestItem &item : items )
  {
    json i = json::object();
    i["item_id"] = item.itemId;
    i["path"] = item.path;
    i["type"] = item.type;
    i["sha256"] = item.sha256;
    i["notes"] = item.notes;
    if( item.derivedFromId.empty() && item.derivedFromName.empty() )
      i["derived_from"] = nullptr;
    else
      i["derived_from"] = { {"name", item.derivedFromName}, {"id", item.derivedFromId},
                            {"version", item.derivedFromVersion} };
    i["intended_override"] = item.intendedOverride;
    item_array.push_back( i );
  }
  j["items"] = item_array;

  return j;
}//Manifest::toJsonObject()


std::string Manifest::toJson() const
{
  return toJsonObject().dump( 2 );
}


const ManifestItem *Manifest::itemForPath( const std::string &relPath ) const
{
  for( const ManifestItem &item : items )
  {
    if( item.path == relPath )
      return &item;
  }
  return nullptr;
}//itemForPath(...)


bool Manifest::needsNewerInterSpec() const
{
  if( formatVersion > sm_formatVersion )
    return true;
  return !minInterSpecVersion.empty() && (compareVersions( minInterSpecVersion, interSpecVersion() ) > 0);
}//needsNewerInterSpec()


bool Manifest::reviewDue() const
{
  if( (reviewBy.size() < 10) || (reviewBy[4] != '-') || (reviewBy[7] != '-') )
    return false;
  return reviewBy.substr( 0, 10 ) < nowIsoUtc().substr( 0, 10 );
}//reviewDue()


// -------------------------------------------------------------------------------------------------
//  Installed and created packages
// -------------------------------------------------------------------------------------------------

std::vector<Package> installedPackages()
{
  vector<Package> answer;
  if( !available() )
    return answer;

  const string pkgs_dir = packagesDirectory();
  if( !SpecUtils::is_directory( pkgs_dir ) )
    return answer;

  const shared_ptr<const json> history_owner = cached_history();
  const json &history = *history_owner;
  const json::const_iterator hist_pkgs = history.find( "packages" );
  const bool have_hist_pkgs = (hist_pkgs != history.end()) && hist_pkgs->is_object();

  set<string> ids_seen;
  vector<string> dirs = SpecUtils::ls_directories_in_directory( pkgs_dir );
  std::sort( begin(dirs), end(dirs) );
  for( const string &dir : dirs )
  {
    if( is_dot_name( dir ) )
      continue;

    const shared_ptr<const Manifest> manifest = cached_manifest( SpecUtils::append_path( dir, ManifestFileName ) );
    if( !manifest )
      continue;

    Package pkg;
    pkg.id = manifest->id;
    pkg.folder = SpecUtils::filename( dir );
    pkg.path = dir;
    pkg.manifest = *manifest;
    pkg.recorded = false;

    const json::const_iterator rec = have_hist_pkgs ? hist_pkgs->find( pkg.id ) : json::const_iterator();
    if( have_hist_pkgs && (rec != hist_pkgs->end()) && rec->is_object() && (json_str( *rec, "folder" ) == pkg.folder) )
    {
      pkg.recorded = true;
      pkg.firstInstalled = json_str( *rec, "first_installed" );

      const json::const_iterator installs = rec->find( "installs" );
      if( (installs != rec->end()) && installs->is_array() )
      {
        for( const json &inst : *installs )
        {
          if( !inst.is_object() )
            continue;
          InstalledVersion v;
          v.installed = json_str( inst, "installed" );
          v.archiveSha256 = json_str( inst, "archive_sha256" );
          try
          {
            const json::const_iterator m = inst.find( "manifest" );
            if( m != inst.end() )
              v.manifest = Manifest::fromJsonObject( *m );
          }catch( std::exception & )
          {
          }
          pkg.installs.push_back( v );
        }//for( loop over installs )
      }//if( installs )

      const json::const_iterator files = rec->find( "files" );
      if( (files != rec->end()) && files->is_object() )
      {
        for( json::const_iterator f = files->begin(); f != files->end(); ++f )
          if( f->is_string() )
            pkg.files[f.key()] = f->get<string>();
      }

    }//if( recorded in the history )

    // Disabled items count for a package added by hand (not recorded) too.
    pkg.disabled = package_disabled( history, pkg.id );

    // Two folders claiming the same package: the recorded one (or the first) wins.
    if( ids_seen.count( pkg.id ) )
    {
      if( !pkg.recorded )
        continue;
      answer.erase( std::remove_if( begin(answer), end(answer), [&pkg]( const Package &p ){
        return p.id == pkg.id;
      } ), end(answer) );
    }
    ids_seen.insert( pkg.id );
    answer.push_back( pkg );
  }//for( loop over package folders )

  // Earliest installed first; folders not in the history last.
  std::stable_sort( begin(answer), end(answer), []( const Package &lhs, const Package &rhs ){
    if( lhs.recorded != rhs.recorded )
      return lhs.recorded;
    if( lhs.firstInstalled != rhs.firstInstalled )
      return lhs.firstInstalled < rhs.firstInstalled;
    return lhs.folder < rhs.folder;
  } );

  return answer;
}//installedPackages()


std::optional<Package> installedPackage( const std::string &packageId )
{
  for( const Package &pkg : installedPackages() )
  {
    if( pkg.id == packageId )
      return pkg;
  }
  return std::nullopt;
}//installedPackage(...)


bool isInPackages( const std::string &path )
{
  if( !available() || path.empty() )
    return false;

  const string pkgs = packagesDirectory();
  if( !relativePath( pkgs, path ).empty() )
    return true;

  // The same folder reached another way (a symbolic link in the user data path, e.g. macOS's /var ->
  //  /private/var, or one path made canonical and the other not): compare the real paths.
  std::error_code pkgs_ec, path_ec;
  const std::filesystem::path real_pkgs = std::filesystem::weakly_canonical( to_fs_path( pkgs ), pkgs_ec );
  const std::filesystem::path real_path = std::filesystem::weakly_canonical( to_fs_path( path ), path_ec );
  return !pkgs_ec && !path_ec && !relativePath( from_fs_path( real_pkgs ), from_fs_path( real_path ) ).empty();
}//isInPackages(...)


void deleteUserItem( const Item &item )
{
#if( USE_USER_LIBRARY )
  if( !available() )
    throw runtime_error( "No user library in this process." );
  if( item.source.kind != SourceKind::User )
    throw runtime_error( "Only your own items can be deleted here." );

  const string user_dir = InterSpec::writableDataDirectory();
  const string path = normalize( item.found.path );
  if( relativePath( user_dir, path ).empty() )
    throw runtime_error( "'" + item.found.path + "' is not in the user data folder." );

  // Never a folder that is, or holds, a type's folder or the packages.
  const auto is_or_holds = [&path]( const string &other ) -> bool {
    return (normalize( other ) == path) || !relativePath( path, other ).empty();
  };
  vector<string> protected_paths{ packagesDirectory(), historyPath() };
  for( const AssetType &type : assetTypes() )
  {
    protected_paths.push_back( SpecUtils::append_path( user_dir, type.folder ) );
    for( const string &folder : type.legacyFolders )
      protected_paths.push_back( SpecUtils::append_path( user_dir, folder ) );
  }
  for( const string &p : protected_paths )
  {
    if( is_or_holds( p ) )
      throw runtime_error( "'" + item.found.path + "' holds other data; delete it outside InterSpec if you are sure." );
  }

  std::error_code ec;
  if( !item.found.isFolder )
  {
    std::filesystem::remove( to_fs_path( path ), ec );
  }else if( item.typeId == "drf" )
  {
    for( const std::filesystem::directory_entry &entry : std::filesystem::directory_iterator( to_fs_path( path ), ec ) )
    {
      std::error_code file_ec;
      if( entry.is_regular_file( file_ec ) )
        std::filesystem::remove( entry.path(), file_ec );
      if( file_ec && !ec )
        ec = file_ec;
    }
    std::error_code dir_ec;
    if( !ec && std::filesystem::is_empty( to_fs_path( path ), dir_ec ) )
      std::filesystem::remove( to_fs_path( path ), dir_ec );
  }else
  {
    std::filesystem::remove_all( to_fs_path( path ), ec );
  }

  notifyChanged();
  if( ec )
    throw runtime_error( ec.message() );
#else
  throw runtime_error( "No user library in this build." );
#endif
}//deleteUserItem(...)


std::optional<Package> packageForPath( const std::string &path )
{
  if( path.empty() )
    return std::nullopt;
  for( const Package &pkg : installedPackages() )
  {
    if( !relativePath( pkg.path, path ).empty() )
      return pkg;
  }
  return std::nullopt;
}//packageForPath(...)


std::map<std::string,std::vector<CreatedVersion>> createdPackages()
{
  map<string,vector<CreatedVersion>> answer;
  const shared_ptr<const json> history_owner = cached_history();
  const json &history = *history_owner;
  const json::const_iterator created = history.find( "created" );
  if( (created == history.end()) || !created->is_object() )
    return answer;

  for( json::const_iterator pkg = created->begin(); pkg != created->end(); ++pkg )
  {
    if( !pkg->is_array() )
      continue;
    for( const json &v : *pkg )
    {
      if( !v.is_object() )
        continue;
      try
      {
        CreatedVersion cv;
        cv.created = json_str( v, "created" );
        cv.manifest = Manifest::fromJsonObject( v.at( "manifest" ) );
        const json::const_iterator sources = v.find( "sources" );
        if( (sources != v.end()) && sources->is_object() )
        {
          for( json::const_iterator s = sources->begin(); s != sources->end(); ++s )
            if( s->is_string() )
              cv.sources[s.key()] = s->get<string>();
        }
        answer[pkg.key()].push_back( cv );
      }catch( std::exception &e )
      {
        cerr << "UserLibrary: skipping a created-package record: " << e.what() << endl;
      }
    }//for( loop over versions )
  }//for( loop over created packages )

  return answer;
}//createdPackages()


std::vector<LooseImport> looseImports()
{
  vector<LooseImport> answer;
  const shared_ptr<const json> history_owner = cached_history();
  const json &history = *history_owner;
  const json::const_iterator loose = history.find( "loose_imports" );
  if( (loose == history.end()) || !loose->is_array() )
    return answer;

  for( const json &l : *loose )
  {
    if( !l.is_object() )
      continue;
    LooseImport li;
    li.relPath = json_str( l, "path" );
    li.originalName = json_str( l, "original_name" );
    li.imported = json_str( l, "imported" );
    li.sha256 = json_str( l, "sha256" );
    answer.push_back( li );
  }
  return answer;
}//looseImports()


void recordLooseImport( const std::string &absPath, const std::string &originalName )
{
  if( !available() )
    return;

  try
  {
    FoundItem item;
    item.path = absPath;
    item.isFolder = SpecUtils::is_directory( absPath );
    const string sha = sha256OfItem( item );
#if( USE_USER_LIBRARY )
    const string rel = relativePath( InterSpec::writableDataDirectory(), absPath );
#else
    const string rel;
#endif
    if( rel.empty() )
      return;

    modifyHistory( [&]( json &history ){
      json &loose = history["loose_imports"];
      if( !loose.is_array() )
        loose = json::array();
      loose.push_back( { {"path", rel}, {"original_name", originalName},
                         {"imported", nowIsoUtc()}, {"sha256", sha} } );
    } );
  }catch( std::exception &e )
  {
    cerr << "UserLibrary::recordLooseImport(): " << e.what() << endl;
  }
}//recordLooseImport(...)


Creator lastCreator()
{
  Creator answer;
  const shared_ptr<const json> history_owner = cached_history();
  const json &history = *history_owner;
  const json::const_iterator c = history.find( "last_creator" );
  if( (c != history.end()) && c->is_object() )
  {
    answer.author = json_str( *c, "author" );
    answer.organization = json_str( *c, "organization" );
    answer.contact = json_str( *c, "contact" );
  }
  return answer;
}//lastCreator()


// -------------------------------------------------------------------------------------------------
//  Choices and disabled items
// -------------------------------------------------------------------------------------------------

std::map<std::string,Choice> choices()
{
  map<string,Choice> answer;
  const shared_ptr<const json> history_owner = cached_history();
  const json &history = *history_owner;
  const json::const_iterator c = history.find( "choices" );
  if( (c == history.end()) || !c->is_object() )
    return answer;

  for( json::const_iterator it = c->begin(); it != c->end(); ++it )
  {
    if( !it->is_object() )
      continue;
    Choice choice;
    choice.source = json_str( *it, "source" );
    choice.sha256 = json_str( *it, "sha256" );
    choice.chosen = json_str( *it, "chosen" );
    answer[it.key()] = choice;
  }
  return answer;
}//choices()


void setChoice( const std::string &typeId, const std::string &name, const Item &item )
{
  string sha;
  try
  {
    sha = sha256OfItem( item.found );
  }catch( std::exception & )
  {
  }

  const AssetType * const type = assetType( typeId );
  modifyHistory( [&]( json &history ){
    json &c = history["choices"];
    if( !c.is_object() )
      c = json::object();

    // One choice per name (a case-insensitive type may have it recorded with other capitalization).
    if( type )
    {
      vector<string> to_remove;
      const string prefix = typeId + ":";
      for( json::iterator it = c.begin(); it != c.end(); ++it )
      {
        if( SpecUtils::starts_with( it.key(), prefix.c_str() )
            && sameName( *type, it.key().substr( prefix.size() ), name ) )
          to_remove.push_back( it.key() );
      }
      for( const string &k : to_remove )
        c.erase( k );
    }//if( type )

    c[typeId + ":" + name] = { {"source", item.source.key()}, {"sha256", sha}, {"chosen", nowIsoUtc()} };
  } );
}//setChoice(...)


void clearChoice( const std::string &typeId, const std::string &name )
{
  const AssetType * const type = assetType( typeId );
  modifyHistory( [&]( json &history ){
    json &c = history["choices"];
    if( !c.is_object() )
    {
      c = json::object();
      return;
    }

    vector<string> to_remove;
    const string prefix = typeId + ":";
    for( json::iterator it = c.begin(); it != c.end(); ++it )
    {
      if( SpecUtils::starts_with( it.key(), prefix.c_str() )
          && (type ? sameName( *type, it.key().substr( prefix.size() ), name ) : (it.key() == prefix + name)) )
        to_remove.push_back( it.key() );
    }
    for( const string &k : to_remove )
      c.erase( k );
  } );
}//clearChoice(...)


void setDisabled( const Item &item, const bool disabled )
{
  if( item.source.kind == SourceKind::BuiltIn )
    throw runtime_error( "InterSpec's own items can not be disabled." );

  modifyHistory( [&]( json &history ){
    if( item.source.kind == SourceKind::Package )
    {
      json &pkgs = history["packages"];
      if( !pkgs.is_object() )
        pkgs = json::object();
      json &pkg = pkgs[item.source.packageId];
      if( !pkg.is_object() )
        pkg = json::object();
      json &list = pkg["disabled"];
      if( !list.is_array() )
        list = json::array();

      const string key = item.itemId.empty() ? item.relPath : item.itemId;
      json updated = json::array();
      for( const json &d : list )
        if( !d.is_string() || ((d.get<string>() != key) && (d.get<string>() != item.relPath)) )
          updated.push_back( d );
      if( disabled )
        updated.push_back( key );
      list = updated;
    }else
    {
      json &list = history["disabled_user"];
      if( !list.is_array() )
        list = json::array();
      json updated = json::array();
      for( const json &d : list )
      {
        if( !d.is_object() || (json_str( d, "type" ) != item.typeId) || (json_str( d, "path" ) != item.relPath) )
          updated.push_back( d );
      }
      if( disabled )
        updated.push_back( { {"type", item.typeId}, {"path", item.relPath} } );
      list = updated;
    }
  } );
}//setDisabled(...)


std::set<std::string> disabledPaths( const std::string &typeId )
{
  set<string> answer;
  if( !available() )
    return answer;

  const shared_ptr<const json> history_owner = cached_history();
  const json &history = *history_owner;
#if( USE_USER_LIBRARY )
  for( const pair<string,string> &d : user_disabled( history ) )
  {
    if( d.first == typeId )
      answer.insert( normalize( SpecUtils::append_path( InterSpec::writableDataDirectory(), d.second ) ) );
  }
#endif

  for( const Package &pkg : installedPackages() )
  {
    for( const string &d : pkg.disabled )
    {
      string rel = d;
      for( const ManifestItem &mi : pkg.manifest.items )
      {
        if( mi.itemId == d )
          rel = mi.path;
      }

      string folder;
      const AssetType * const type = typeForRelativePath( rel, folder );
      if( type && (type->id == typeId) )
        answer.insert( normalize( SpecUtils::append_path( pkg.path, rel ) ) );
    }//for( loop over disabled entries )
  }//for( loop over packages )

  return answer;
}//disabledPaths(...)


// -------------------------------------------------------------------------------------------------
//  History file
// -------------------------------------------------------------------------------------------------

void modifyHistory( const std::function<void( nlohmann::json &history )> &fcn )
{
  if( !available() )
    throw runtime_error( "No user library in this process." );

  std::lock_guard<std::recursive_mutex> modify_lock( ns_modify_mutex );

  const string path = historyPath();
  json history = json::object();
  if( SpecUtils::is_file( path ) )
  {
    // A file we can't read (locked, too large, an I/O error) is left alone, and nothing is changed;
    //  only one that is not valid JSON is set aside (kept for the user, or a developer, to look at).
    string contents;
    if( !read_file( path, contents, ns_max_json_bytes ) )
      throw runtime_error( "Could not read '" + path + "' (or it is over "
                           + std::to_string( ns_max_json_bytes / (1024*1024) ) + " MB)." );
    try
    {
      history = json::parse( contents );
      if( !history.is_object() )
        throw runtime_error( "not a JSON object" );
    }catch( std::exception &e )
    {
      const string backup = path + ".unreadable-" + std::to_string( std::time(nullptr) );
      cerr << "UserLibrary: history file is not valid (" << e.what() << "); moving it to '" << backup << "'" << endl;
      std::error_code ec;
      std::filesystem::rename( to_fs_path( path ), to_fs_path( backup ), ec );
      if( ec )
        throw runtime_error( "History file '" + path + "' is not valid, and could not be set aside: " + ec.message() );
      history = json::object();
    }
  }//if( SpecUtils::is_file( path ) )

  fcn( history );
  history["format_version"] = 1;

  const string tmp_path = path + ".tmp-" + newUuid().substr( 0, 8 );
  {
    std::ofstream output( to_fs_path( tmp_path ), ios::out | ios::binary | ios::trunc );
    if( !output )
      throw runtime_error( "Could not write '" + tmp_path + "'" );
    output << history.dump( 2 );
    output.flush();
    if( !output )
    {
      output.close();
      SpecUtils::remove_file( tmp_path );
      throw runtime_error( "Failed writing '" + tmp_path + "'" );
    }
  }

  std::error_code ec;
  std::filesystem::rename( to_fs_path( tmp_path ), to_fs_path( path ), ec );
  if( ec )
  {
    SpecUtils::remove_file( tmp_path );
    throw runtime_error( "Could not replace '" + path + "': " + ec.message() );
  }

  {
    std::lock_guard<std::mutex> lock( ns_cache_mutex );
    ns_history_cache.path = path;
    ns_history_cache.stamp = stamp_of( path );
    ns_history_cache.value = make_shared<const json>( history );
    ns_history_cache.loaded = true;
  }

  notifyChanged();
}//modifyHistory(...)


nlohmann::json readHistory()
{
  return *cached_history();
}


uint64_t generation()
{
  return ns_generation.load();
}


void notifyChanged()
{
  {
    std::lock_guard<std::mutex> lock( ns_cache_mutex );
    ns_manifest_cache.clear();
  }
  ++ns_generation;
}//notifyChanged()


void cleanupStaleFolders()
{
  if( !available() )
    return;

  try
  {
    const string dir = packagesDirectory();
    if( !SpecUtils::is_directory( dir ) )
      return;

    const shared_ptr<const json> history_owner = cached_history();
    const json &history = *history_owner;
    const json::const_iterator packages = history.find( "packages" );

    for( const string &sub : SpecUtils::ls_directories_in_directory( dir ) )
    {
      const string name = SpecUtils::filename( sub );
      const bool incoming = SpecUtils::starts_with( name, ".incoming-" );
      const bool old = SpecUtils::starts_with( name, ".old-" );
      if( !incoming && !old )
        continue;

      // Left alone for an hour: another InterSpec using this folder may be mid-install.
      std::error_code ec;
      const std::filesystem::file_time_type modified = std::filesystem::last_write_time( to_fs_path( sub ), ec );
      if( ec || ((std::filesystem::file_time_type::clock::now() - modified) < std::chrono::hours( 1 )) )
        continue;

      // An installed version set aside during an update that never finished: put it back, if its
      //  package's folder is missing.
      if( old && (packages != history.end()) && packages->is_object() )
      {
        try
        {
          string contents;
          const string manifest_path = SpecUtils::append_path( sub, ManifestFileName );
          if( read_file( manifest_path, contents, ns_max_json_bytes ) )
          {
            const Manifest m = Manifest::fromJson( contents );
            const json::const_iterator rec = packages->find( m.id );
            const string folder = ((rec != packages->end()) && rec->is_object()) ? rec->value( "folder", string() ) : string();
            const string folder_path = folder.empty() ? string() : SpecUtils::append_path( dir, folder );
            if( !folder.empty() && (folder.find_first_of( "/\\" ) == string::npos) && !SpecUtils::is_directory( folder_path ) )
            {
              std::filesystem::rename( to_fs_path( sub ), to_fs_path( folder_path ), ec );
              if( !ec )
              {
                cerr << "UserLibrary: restored '" << folder_path << "' from an unfinished update." << endl;
                continue;
              }
            }
          }//if( it has a manifest )
        }catch( std::exception &e )
        {
          cerr << "UserLibrary: could not check '" << sub << "': " << e.what() << endl;
        }
      }//if( old )

      ec.clear();
      std::filesystem::remove_all( to_fs_path( sub ), ec );
      if( ec )
        cerr << "UserLibrary: could not remove '" << sub << "': " << ec.message() << endl;
    }//for( loop over sub-folders )
  }catch( std::exception &e )
  {
    cerr << "UserLibrary::cleanupStaleFolders(): " << e.what() << endl;
  }
}//cleanupStaleFolders()


// -------------------------------------------------------------------------------------------------
//  Utilities
// -------------------------------------------------------------------------------------------------

std::string sha256OfData( const void *data, const size_t size )
{
  return ceelo::env_fit_sha256_hex( data, size );
}


std::string sha256OfFile( const std::string &path )
{
  string contents;
  if( !read_file( path, contents, std::numeric_limits<size_t>::max() ) )
    throw runtime_error( "Could not read '" + path + "'" );
  return sha256OfData( contents.data(), contents.size() );
}//sha256OfFile(...)


std::string sha256OfItem( const FoundItem &item )
{
  if( !item.isFolder )
    return sha256OfFile( item.path );

  vector<string> lines;
  for( const string &f : list_files( item.path, {}, true ) )
    lines.push_back( relativePath( item.path, f ) + "\t" + sha256OfFile( f ) + "\n" );
  std::sort( begin(lines), end(lines) );

  string all;
  for( const string &l : lines )
    all += l;
  return sha256OfData( all.data(), all.size() );
}//sha256OfItem(...)


int compareVersions( const std::string &lhs, const std::string &rhs, bool *numeric )
{
  auto split_numbers = []( const string &v, vector<unsigned long long> &parts ) -> bool {
    parts.clear();
    vector<string> fields;
    SpecUtils::split( fields, SpecUtils::trim_copy( v ), "." );
    if( fields.empty() )
      return false;
    for( const string &f : fields )
    {
      if( f.empty() || (f.size() > 18) || !std::all_of( begin(f), end(f), []( char c ){ return (c >= '0') && (c <= '9'); } ) )
        return false;
      parts.push_back( std::stoull( f ) );
    }
    return true;
  };//split_numbers

  vector<unsigned long long> a, b;
  const bool is_numeric = split_numbers( lhs, a ) && split_numbers( rhs, b );
  if( numeric )
    *numeric = is_numeric;

  if( !is_numeric )
    return (lhs < rhs) ? -1 : ((lhs == rhs) ? 0 : 1);

  const size_t n = std::max( a.size(), b.size() );
  a.resize( n, 0 );
  b.resize( n, 0 );
  for( size_t i = 0; i < n; ++i )
  {
    if( a[i] != b[i] )
      return (a[i] < b[i]) ? -1 : 1;
  }
  return 0;
}//compareVersions(...)


std::string interSpecVersion()
{
#ifdef InterSpec_VERSION
  return InterSpec_VERSION;
#else
  return "0";
#endif
}


std::string nowIsoUtc()
{
  const std::time_t t = std::time( nullptr );
  std::tm tm{};
#ifdef _WIN32
  gmtime_s( &tm, &t );
#else
  gmtime_r( &t, &tm );
#endif
  char buffer[32] = { '\0' };
  std::strftime( buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &tm );
  return buffer;
}//nowIsoUtc()


std::string newUuid()
{
  static std::mutex m;
  static std::mt19937_64 gen = [](){
    std::random_device rd;
    const uint64_t seed = (static_cast<uint64_t>(rd()) << 32) ^ rd()
               ^ static_cast<uint64_t>( std::chrono::high_resolution_clock::now().time_since_epoch().count() );
    return std::mt19937_64( seed );
  }();

  uint64_t hi, lo;
  {
    std::lock_guard<std::mutex> lock( m );
    hi = gen();
    lo = gen();
  }

  hi = (hi & 0xFFFFFFFFFFFF0FFFULL) | 0x0000000000004000ULL;  // version 4
  lo = (lo & 0x3FFFFFFFFFFFFFFFULL) | 0x8000000000000000ULL;  // RFC 4122 variant

  char buffer[40];
  snprintf( buffer, sizeof(buffer), "%08x-%04x-%04x-%04x-%012llx",
            static_cast<unsigned>( hi >> 32 ), static_cast<unsigned>( (hi >> 16) & 0xFFFF ),
            static_cast<unsigned>( hi & 0xFFFF ), static_cast<unsigned>( lo >> 48 ),
            static_cast<unsigned long long>( lo & 0xFFFFFFFFFFFFULL ) );
  return buffer;
}//newUuid()


std::string collisionKey( const std::string &utf8 )
{
  // Base letters of U+00C0 to U+017F (Latin-1 Supplement and Latin Extended-A), so a precomposed
  //  letter and its decomposed form (macOS file names are often NFD) collide.  '\0' keeps the char.
  static const char latin1[] = "aaaaaaaceeeeiiiidnooooo\0ouuuuyts" "aaaaaaaceeeeiiiidnooooo\0ouuuuyty";
  static const char latin_ext_a[] = "aaaaaaccccccccddddeeeeeeeeeegggggggghhhhiiiiiiiiiiiijjkkklllllll"
                                    "lllnnnnnnnnnoooooooorrrrrrssssssssttttttuuuuuuuuuuuuwwyyyzzzzzzs";

  string answer;
  answer.reserve( utf8.size() );
  for( size_t i = 0; i < utf8.size(); )
  {
    const unsigned char c = static_cast<unsigned char>( utf8[i] );
    if( c < 0x80 )
    {
      answer.push_back( static_cast<char>( std::tolower( c ) ) );
      i += 1;
      continue;
    }

    size_t len = 1;
    uint32_t cp = 0;
    if( ((c & 0xE0) == 0xC0) && ((i + 1) < utf8.size()) )
    {
      len = 2;
      cp = ((c & 0x1F) << 6) | (static_cast<unsigned char>(utf8[i+1]) & 0x3F);
    }else if( ((c & 0xF0) == 0xE0) && ((i + 2) < utf8.size()) )
    {
      len = 3;
    }else if( ((c & 0xF8) == 0xF0) && ((i + 3) < utf8.size()) )
    {
      len = 4;
    }

    if( (cp >= 0x300) && (cp <= 0x36F) )
    {
      // combining mark: dropped
    }else if( (cp >= 0xC0) && (cp <= 0xFF) && latin1[cp - 0xC0] )
    {
      answer.push_back( latin1[cp - 0xC0] );
    }else if( (cp >= 0x100) && (cp <= 0x17F) && ((cp - 0x100) < (sizeof(latin_ext_a) - 1)) )
    {
      answer.push_back( latin_ext_a[cp - 0x100] );
    }else
    {
      answer.append( utf8, i, len );
    }
    i += len;
  }//for( loop over bytes )

  return answer;
}//collisionKey(...)


bool isSecretPath( const std::string &relPath )
{
  string lower = relPath;
  SpecUtils::to_lower_ascii( lower );
  SpecUtils::ireplace_all( lower, "\\", "/" );

  vector<string> parts;
  SpecUtils::split( parts, lower, "/" );
  if( parts.empty() )
    return false;

  // Folders that only hold logs, test data, or other packages.
  const char * const secret_folders[] = { "llm_logs", "llm_benchmarks", "test_datasets", "packages" };
  for( size_t i = 0; (i + 1) < parts.size(); ++i )
  {
    for( const char * const f : secret_folders )
      if( parts[i] == f )
        return true;
  }

  // API keys and endpoints, the user database, and user library's own bookkeeping.
  const string &name = parts.back();
  const char * const secret_files[] = { "llm_config.xml", "llm_tools_config.xml", "arcgis_key.txt",
                                        "prev_spec_uri.xml", "user_library.json" };
  for( const char * const f : secret_files )
    if( name == f )
      return true;

  const char * const secret_exts[] = { ".db", ".sqlite", ".sqlite3", ".db-journal", ".db-wal", ".db-shm" };
  for( const char * const e : secret_exts )
    if( SpecUtils::iends_with( name, e ) )
      return true;

  return SpecUtils::starts_with( name, "interspecuserdata" );
}//isSecretPath(...)


std::vector<std::string> templateIncludesOutsideFolder( const std::string &templatePath )
{
  static const std::regex include_re( "\\{%-?\\s*(?:include|extends|import)\\s+[\"']([^\"']+)[\"']" );

  UserLibrary::FoundItem item;
  item.path = templatePath;
  const string dir = SpecUtils::parent_path( templatePath );
  const vector<string> deps = template_dependencies( item );

  set<string> allowed;
  for( const string &dep : deps )
    allowed.insert( relativePath( dir, dep ) );

  vector<string> answer;
  vector<string> to_check{ templatePath };
  to_check.insert( end(to_check), begin(deps), end(deps) );
  for( const string &f : to_check )
  {
    string contents;
    if( !read_file( f, contents, 16*1024*1024 ) )
      continue;
    for( std::sregex_iterator it( begin(contents), end(contents), include_re ), it_end; it != it_end; ++it )
    {
      const string target = (*it)[1].str();
      // InterSpec's built-in partials are registered by name, not read from files.
      if( !SpecUtils::starts_with( target, "default-" ) && !allowed.count( target )
          && (std::find( begin(answer), end(answer), target ) == end(answer)) )
        answer.push_back( target );
    }
  }//for( loop over the template and its includes )

  return answer;
}//templateIncludesOutsideFolder(...)


void confineTemplateIncludes( inja::Environment &env, const std::string &templatePath )
{
  env.set_search_included_templates_in_files( false );

  UserLibrary::FoundItem item;
  item.path = templatePath;
  const string dir = SpecUtils::parent_path( templatePath );
  for( const string &dep : template_dependencies( item ) )
  {
    string text;
    if( !read_file( dep, text, 16*1024*1024 ) )
      throw runtime_error( "Could not read '" + dep + "'" );
    env.include_template( relativePath( dir, dep ), env.parse( text ) );
  }
}//confineTemplateIncludes(...)


std::string relativePath( const std::string &base, const std::string &path )
{
  string b = normalize( base ), p = normalize( path );
  SpecUtils::ireplace_all( b, "\\", "/" );
  SpecUtils::ireplace_all( p, "\\", "/" );
  while( !b.empty() && (b.back() == '/') )
    b.pop_back();

  if( b.empty() || (p.size() <= (b.size() + 1)) || (p[b.size()] != '/') || !SpecUtils::starts_with( p, b.c_str() ) )
    return string();

  string answer = p.substr( b.size() + 1 );
  while( !answer.empty() && (answer.back() == '/') )
    answer.pop_back();
  return answer;
}//relativePath(...)


std::string sanitizeFileName( const std::string &name, const std::string &fallback )
{
  string answer;
  for( const char c : name )
  {
    const unsigned char uc = static_cast<unsigned char>( c );
    if( (uc < 0x20) || (string( "\\/:?\"<>|*" ).find( c ) != string::npos) )
      answer.push_back( '_' );
    else
      answer.push_back( c );
  }

  answer = truncate_utf8( SpecUtils::trim_copy( answer ), 60 );
  while( !answer.empty() && ((answer.back() == '.') || (answer.back() == ' ')) )
    answer.pop_back();
  while( !answer.empty() && (answer.front() == '.') )
    answer.erase( 0, 1 );

  // Names Windows reserves for devices (with or without an extension).
  string stem = answer.substr( 0, answer.find( '.' ) );
  SpecUtils::to_lower_ascii( stem );
  const bool reserved = (stem == "con") || (stem == "prn") || (stem == "aux") || (stem == "nul")
                        || (((stem.size() == 4) && (SpecUtils::starts_with( stem, "com" ) || SpecUtils::starts_with( stem, "lpt" )))
                            && (stem[3] >= '0') && (stem[3] <= '9'));
  if( reserved )
    answer = "_" + answer;

  return answer.empty() ? fallback : answer;
}//sanitizeFileName(...)
}//namespace UserLibrary
