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

/** Tests of the user library (InterSpec/UserLibrary.h) and its packages
 (InterSpec/UserLibraryPackage.h): archive entry checks, manifest parsing, version comparison,
 validation statuses (including a format too new), secrets refused, and a
 create -> install -> update -> uninstall round trip between two user data folders, with the
 recorded choices and disabled items.
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#define BOOST_TEST_MODULE TestUserLibrary
#include <boost/test/included/unit_test.hpp>

#include <set>
#include <chrono>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <filesystem>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"
#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/ZipArchive.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ColorThemeWindow.h"
#include "InterSpec/DrfLibrary.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/UserLibrary.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/UserLibraryPackage.h"
#include "InterSpec/DetectorPeakResponse.h"

#if( USE_LLM_INTERFACE )
#include "InterSpec/LlmConfig.h"
#endif

using namespace std;

namespace
{
string g_data_dir;

void set_dirs()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;
  s_have_set = true;

  const int argc = boost::unit_test::framework::master_test_suite().argc;
  char ** const argv = boost::unit_test::framework::master_test_suite().argv;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      g_data_dir = arg.substr( 10 );
  }
  SpecUtils::ireplace_all( g_data_dir, "%20", " " );

  if( g_data_dir.empty() )
  {
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        g_data_dir = d;
        break;
      }
    }
  }//if( g_data_dir.empty() )

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory( SpecUtils::append_path(g_data_dir, "rel_act") ),
                         "No data dir found; pass --datadir=..." );
  InterSpec::setStaticDataDirectory( g_data_dir );
  BOOST_REQUIRE( DecayDataBaseServer::database() );
}//set_dirs()


std::filesystem::path fs_path( const string &p )
{
#ifdef _WIN32
  return std::filesystem::path( SpecUtils::convert_from_utf8_to_utf16( p ) );
#else
  return std::filesystem::path( p );
#endif
}


/** A fresh, empty temporary directory. */
string make_temp_dir( const string &tag )
{
  const string dir = SpecUtils::append_path( SpecUtils::temp_dir(),
                                             "UserLibTest_" + tag + "_" + UserLibrary::newUuid().substr( 0, 8 ) );
  std::filesystem::remove_all( fs_path( dir ) );
  std::filesystem::create_directories( fs_path( dir ) );
  return dir;
}


void write_file( const string &path, const string &contents )
{
  std::filesystem::create_directories( fs_path( SpecUtils::parent_path( path ) ) );
  std::ofstream out( fs_path( path ), ios::out | ios::binary | ios::trunc );
  out << contents;
  BOOST_REQUIRE( out.good() );
}


string read_file( const string &path )
{
  std::ifstream in( fs_path( path ), ios::in | ios::binary );
  std::ostringstream strm;
  strm << in.rdbuf();
  return strm.str();
}


/** A built-in isotopics preset's text, optionally with a note changed (so contents differ). */
string preset_text( const string &note = "" )
{
  string text = read_file( SpecUtils::append_path( SpecUtils::append_path( g_data_dir, "rel_act" ),
                                                   "HPGe U (120-1001 keV).xml" ) );
  BOOST_REQUIRE( text.size() > 100 );
  if( !note.empty() )
  {
    const size_t pos = text.find( "<Description>" );
    BOOST_REQUIRE( pos != string::npos );
    text.insert( pos + 13, note + " " );
  }
  return text;
}


void copy_folder( const string &from, const string &to )
{
  std::filesystem::create_directories( fs_path( to ) );
  std::filesystem::copy( fs_path( from ), fs_path( to ), std::filesystem::copy_options::recursive );
}


void use_user_dir( const string &dir )
{
  InterSpec::setWritableDataDirectory( dir );
  UserLibrary::notifyChanged();
  BOOST_REQUIRE( UserLibrary::available() );
}


/** Writes a zip of the given (name, contents) entries. */
string make_zip( const string &dir, const string &filename, const vector<pair<string,string>> &entries )
{
  const string path = SpecUtils::append_path( dir, filename );
  std::ofstream out( fs_path( path ), ios::out | ios::binary | ios::trunc );
  ZipArchive::ZipWriter zip( out );
  for( const pair<string,string> &e : entries )
    zip.add_file( e.first, e.second );
  zip.finish();
  out.close();
  return path;
}


string create_package( const UserLibraryPackage::CreateRequest &request, const string &path )
{
  std::ofstream out( fs_path( path ), ios::out | ios::binary | ios::trunc );
  const UserLibrary::Manifest m = UserLibraryPackage::create( request, out );
  out.close();
  return m.id;
}


vector<UserLibrary::Item> user_items( const string &type )
{
  vector<UserLibrary::Item> answer;
  for( const UserLibrary::Item &item : UserLibrary::items( type ) )
    if( item.source.kind == UserLibrary::SourceKind::User )
      answer.push_back( item );
  return answer;
}


const UserLibraryPackage::IncomingItem *incoming( const UserLibraryPackage::Inspection &insp, const string &path )
{
  for( const UserLibraryPackage::IncomingItem &item : insp.items )
    if( item.relPath == path )
      return &item;
  return nullptr;
}
}//namespace


BOOST_AUTO_TEST_CASE( entry_names )
{
  const UserLibraryPackage::Limits limits;
  using UserLibraryPackage::checkEntryName;

  BOOST_CHECK_EQUAL( checkEntryName( "drfs/a.xml", limits ), "drfs/a.xml" );
  BOOST_CHECK_EQUAL( checkEntryName( "drfs\\b.tsv", limits ), "drfs/b.tsv" );
  BOOST_CHECK_EQUAL( checkEntryName( "rel_act/sub/", limits ), "rel_act/sub" );
  BOOST_CHECK_EQUAL( checkEntryName( "rel_act/D\xC3\xA9tecteur.xml", limits ), "rel_act/D\xC3\xA9tecteur.xml" );

  const char * const bad[] = {
    "/etc/passwd", "\\windows\\system32", "C:/x.txt", "c:x.txt", "../x", "drfs/../../x", "drfs/./x",
    "a//b", "con.txt", "drfs/aux", "drfs/COM1.xml", "drfs/x:y", "drfs/trailing.", "drfs/space ",
    "drfs/a\x01" "b", "drfs/bad\xC3", "drfs/a?b", ""
  };
  for( const char * const name : bad )
    BOOST_CHECK_THROW( checkEntryName( name, limits ), std::exception );

  BOOST_CHECK_THROW( checkEntryName( "drfs/" + string( 120, 'a' ), limits ), std::exception );
  BOOST_CHECK_THROW( checkEntryName( string( 100, 'a' ) + "/" + string( 70, 'b' ), limits ), std::exception );
}//BOOST_AUTO_TEST_CASE( entry_names )


BOOST_AUTO_TEST_CASE( collision_keys_and_secrets )
{
  // NFC "é" vs NFD "e" + U+0301, and capitalization, all collide.
  BOOST_CHECK_EQUAL( UserLibrary::collisionKey( "D\xC3\xA9tecteur" ), UserLibrary::collisionKey( "De\xCC\x81tecteur" ) );
  BOOST_CHECK_EQUAL( UserLibrary::collisionKey( "D\xC3\xA9tecteur" ), UserLibrary::collisionKey( "DETECTEUR" ) );
  BOOST_CHECK( UserLibrary::collisionKey( "abc" ) != UserLibrary::collisionKey( "abd" ) );
  BOOST_CHECK_EQUAL( UserLibrary::collisionKey( "\xC5\x81\xC3\xB3" ), "lo" );  // "Łó"

  BOOST_CHECK( UserLibrary::isSecretPath( "llm_config.xml" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "drfs/LLM_CONFIG.xml" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "llm_tools_config.xml" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "arcgis_key.txt" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "InterSpecUserData.db" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "prev_spec_uri.xml" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "llm_logs/a.txt" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "llm_benchmarks/x/y.xml" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "test_datasets/x.n42" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "packages/other/interspec-package.json" ) );
  BOOST_CHECK( UserLibrary::isSecretPath( "user_library.json" ) );
  BOOST_CHECK( !UserLibrary::isSecretPath( "drfs/a.xml" ) );
  BOOST_CHECK( !UserLibrary::isSecretPath( "rel_act/llm_logs.xml" ) );

  BOOST_CHECK_EQUAL( UserLibrary::sanitizeFileName( "Team: HPGe/Kit?", "x" ), "Team_ HPGe_Kit_" );
  BOOST_CHECK_EQUAL( UserLibrary::sanitizeFileName( "...", "fallback" ), "fallback" );
  BOOST_CHECK_EQUAL( UserLibrary::sanitizeFileName( "CON", "x" ), "_CON" );
}//BOOST_AUTO_TEST_CASE( collision_keys_and_secrets )


BOOST_AUTO_TEST_CASE( versions )
{
  bool numeric = false;
  BOOST_CHECK( UserLibrary::compareVersions( "1.10", "1.9", &numeric ) > 0 );
  BOOST_CHECK( numeric );
  BOOST_CHECK_EQUAL( UserLibrary::compareVersions( "2", "2.0" ), 0 );
  BOOST_CHECK( UserLibrary::compareVersions( "1.0.14", "1.0.15" ) < 0 );
  BOOST_CHECK( UserLibrary::compareVersions( "3", "12" ) < 0 );
  BOOST_CHECK( UserLibrary::compareVersions( "beta", "alpha", &numeric ) > 0 );
  BOOST_CHECK( !numeric );

  BOOST_CHECK_EQUAL( UserLibraryPackage::nextVersion( "3" ), "4" );
  BOOST_CHECK_EQUAL( UserLibraryPackage::nextVersion( "1.9" ), "1.10" );
  BOOST_CHECK_EQUAL( UserLibraryPackage::nextVersion( "v2" ), "v2.1" );
  BOOST_CHECK_EQUAL( UserLibraryPackage::nextVersion( "" ), "1" );
}//BOOST_AUTO_TEST_CASE( versions )


BOOST_AUTO_TEST_CASE( manifest_parsing )
{
  UserLibrary::Manifest m;
  m.id = UserLibrary::newUuid();
  m.name = "Field Team Kit";
  m.version = "3";
  m.author = "A. Person";
  m.handling = "Internal use only";
  m.previousVersions.push_back( { "2", "2026-01-01T00:00:00Z", "Added the LaBr detector" } );
  UserLibrary::ManifestItem item;
  item.itemId = UserLibrary::newUuid();
  item.path = "drfs/Det.xml";
  item.type = "drf";
  item.sha256 = "abc";
  item.intendedOverride = true;
  item.derivedFromName = "Other";
  item.derivedFromId = "other-id";
  item.derivedFromVersion = "7";
  m.items.push_back( item );

  const UserLibrary::Manifest back = UserLibrary::Manifest::fromJson( m.toJson() );
  BOOST_CHECK_EQUAL( back.id, m.id );
  BOOST_CHECK_EQUAL( back.name, m.name );
  BOOST_CHECK_EQUAL( back.version, "3" );
  BOOST_CHECK_EQUAL( back.handling, m.handling );
  BOOST_REQUIRE_EQUAL( back.previousVersions.size(), 1 );
  BOOST_CHECK_EQUAL( back.previousVersions[0].changes, "Added the LaBr detector" );
  BOOST_REQUIRE_EQUAL( back.items.size(), 1 );
  BOOST_CHECK_EQUAL( back.items[0].path, "drfs/Det.xml" );
  BOOST_CHECK( back.items[0].intendedOverride );
  BOOST_CHECK_EQUAL( back.items[0].derivedFromId, "other-id" );
  BOOST_CHECK( !back.needsNewerInterSpec() );

  // Versions as numbers, trailing slashes, and unknown fields are fine.
  const UserLibrary::Manifest loose = UserLibrary::Manifest::fromJson(
    R"({"format_version":1,"id":"x","name":"n","version":2,"unknown":{"a":1},
        "items":[{"path":"rel_act/a.xml/","type":"rel_act_preset"}]})" );
  BOOST_CHECK_EQUAL( loose.version, "2" );
  BOOST_CHECK_EQUAL( loose.items.at(0).path, "rel_act/a.xml" );

  BOOST_CHECK_THROW( UserLibrary::Manifest::fromJson( R"({"id":"x","name":"n","version":"1"})" ), std::exception );
  BOOST_CHECK_THROW( UserLibrary::Manifest::fromJson( R"({"format_version":1,"name":"n","version":"1"})" ), std::exception );
  BOOST_CHECK_THROW( UserLibrary::Manifest::fromJson( R"({"format_version":1,"id":"x","name":"n","version":"1","items":[{"type":"drf"}]})" ), std::exception );
  BOOST_CHECK_THROW( UserLibrary::Manifest::fromJson( "not json" ), std::exception );

  UserLibrary::Manifest newer = back;
  newer.formatVersion = 2;
  BOOST_CHECK( newer.needsNewerInterSpec() );
  newer.formatVersion = 1;
  newer.minInterSpecVersion = "999.0";
  BOOST_CHECK( newer.needsNewerInterSpec() );
  newer.minInterSpecVersion = "1.0.0";
  BOOST_CHECK( !newer.needsNewerInterSpec() );
}//BOOST_AUTO_TEST_CASE( manifest_parsing )


BOOST_AUTO_TEST_CASE( zip_writer_round_trip )
{
  const string dir = make_temp_dir( "zip" );
  const string big( 100000, 'x' );
  const string path = make_zip( dir, "a.zip", { {"a.txt", "hello"}, {"sub/b.bin", big}, {"empty.txt", ""} } );

  std::ifstream in( fs_path( path ), ios::in | ios::binary );
  const vector<shared_ptr<const ZipArchive::ZipFileHeader>> entries = ZipArchive::list_zip_entries( in );
  BOOST_REQUIRE_EQUAL( entries.size(), 3 );
  BOOST_CHECK_EQUAL( entries[1]->filename, "sub/b.bin" );
  BOOST_CHECK( !entries[1]->is_symlink() );

  ZipArchive::ExtractionLimits limits;
  limits.max_uncompressed_bytes = 1024*1024;
  limits.max_compression_ratio = 2000;
  std::ostringstream out;
  ZipArchive::read_file_from_zip( in, entries[1], out, limits );
  BOOST_CHECK( out.str() == big );

  std::filesystem::remove_all( fs_path( dir ) );
}//BOOST_AUTO_TEST_CASE( zip_writer_round_trip )


BOOST_AUTO_TEST_CASE( hostile_archives )
{
  set_dirs();
  const string user = make_temp_dir( "hostile_user" );
  const string work = make_temp_dir( "hostile_work" );
  use_user_dir( user );

  const string manifest = R"({"format_version":1,"id":"hostile","name":"Hostile","version":"1"})";

  // Escaping the package folder.
  {
    const string pkg = make_zip( work, "escape.ispkg", { {"interspec-package.json", manifest}, {"../evil.txt", "x"} } );
    const shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "escape.ispkg" );
    BOOST_CHECK( !insp->error.empty() );
    BOOST_CHECK( insp->stagingDir.empty() );
    BOOST_CHECK( !SpecUtils::is_file( SpecUtils::append_path( UserLibrary::packagesDirectory(), "evil.txt" ) ) );
  }

  {
    const string pkg = make_zip( work, "abs.ispkg", { {"interspec-package.json", manifest}, {"/tmp/evil.txt", "x"} } );
    BOOST_CHECK( !UserLibraryPackage::inspect( pkg, "abs.ispkg" )->error.empty() );
  }

  // Names differing only by case.
  {
    const string pkg = make_zip( work, "case.ispkg", { {"interspec-package.json", manifest},
                                 {"rel_act/A.xml", preset_text()}, {"rel_act/a.xml", preset_text()} } );
    BOOST_CHECK( !UserLibraryPackage::inspect( pkg, "case.ispkg" )->error.empty() );
  }

  // Limits on count and size.
  {
    UserLibraryPackage::Limits limits;
    limits.maxEntries = 2;
    const string pkg = make_zip( work, "count.ispkg", { {"interspec-package.json", manifest},
                                 {"rel_act/a.xml", "x"}, {"rel_act/b.xml", "y"} } );
    BOOST_CHECK( !UserLibraryPackage::inspect( pkg, "count.ispkg", limits )->error.empty() );

    limits = UserLibraryPackage::Limits();
    limits.maxTotalBytes = 1000;
    const string big = make_zip( work, "big.ispkg", { {"interspec-package.json", manifest},
                                 {"rel_act/a.xml", string( 5000, 'z' )} } );
    BOOST_CHECK( !UserLibraryPackage::inspect( big, "big.ispkg", limits )->error.empty() );
  }

  // A symbolic link entry: patch the central directory of an ordinary archive to say so.
  {
    const string pkg = make_zip( work, "link.ispkg", { {"interspec-package.json", manifest},
                                 {"rel_act/link.xml", "/etc/passwd"} } );
    string bytes = read_file( pkg );
    size_t pos = bytes.rfind( string( "\x50\x4b\x01\x02", 4 ) );  // last central directory header
    BOOST_REQUIRE( pos != string::npos );
    const uint32_t attrs = 0120777u << 16;
    memcpy( &bytes[pos + 38], &attrs, 4 );
    write_file( pkg, bytes );
    const shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "link.ispkg" );
    BOOST_CHECK( insp->error.find( "symbolic link" ) != string::npos );
  }

  // Secrets are refused, unknown files ignored; neither is installed.
  {
    const string pkg = make_zip( work, "secrets.ispkg", { {"interspec-package.json", manifest},
                                 {"llm_config.xml", "<key>sk-123</key>"},
                                 {"InterSpecUserData.db", "db"},
                                 {"notes/readme.txt", "hello"},
                                 {"rel_act/Good.xml", preset_text()} } );
    shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "secrets.ispkg" );
    BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
    size_t num_refused = 0, num_ignored = 0;
    for( const UserLibraryPackage::SkippedFile &s : insp->skipped )
      (s.refused ? num_refused : num_ignored) += 1;
    BOOST_CHECK_EQUAL( num_refused, 2 );
    BOOST_CHECK_EQUAL( num_ignored, 1 );
    BOOST_REQUIRE_EQUAL( insp->items.size(), 1 );
    BOOST_CHECK_EQUAL( insp->items[0].typeId, "rel_act_preset" );

    const UserLibrary::Package installed = UserLibraryPackage::install( *insp, {} );
    BOOST_CHECK( !SpecUtils::is_file( SpecUtils::append_path( installed.path, "llm_config.xml" ) ) );
    BOOST_CHECK( !SpecUtils::is_file( SpecUtils::append_path( installed.path, "InterSpecUserData.db" ) ) );
    BOOST_CHECK( !SpecUtils::is_file( SpecUtils::append_path( installed.path, "notes/readme.txt" ) ) );
    BOOST_CHECK( SpecUtils::is_file( SpecUtils::append_path( installed.path, "rel_act/Good.xml" ) ) );
    UserLibraryPackage::uninstall( installed.id );
  }

  // Nothing left behind in the packages folder.
  for( const string &d : SpecUtils::ls_directories_in_directory( UserLibrary::packagesDirectory() ) )
    BOOST_ERROR( "Left behind: " << d );

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( hostile_archives )


BOOST_AUTO_TEST_CASE( validation_statuses )
{
  set_dirs();
  const string user = make_temp_dir( "validation" );
  use_user_dir( user );

  const string rel_act = SpecUtils::append_path( user, "rel_act" );
  write_file( SpecUtils::append_path( rel_act, "Fine.xml" ), preset_text( "fine" ) );
  write_file( SpecUtils::append_path( rel_act, "Garbage.xml" ), "<RelActCalcAuto version=\"0\"><nope" );

  // Written by a newer InterSpec: the top-level version, or a nested element's version.
  string newer_top = preset_text();
  SpecUtils::ireplace_all( newer_top, "<RelActCalcAuto version=\"0\">", "<RelActCalcAuto version=\"1\">" );
  write_file( SpecUtils::append_path( rel_act, "NewerTop.xml" ), newer_top );
  string newer_options = preset_text();
  SpecUtils::ireplace_all( newer_options, "<Options version=\"7\"", "<Options version=\"99\"" );
  write_file( SpecUtils::append_path( rel_act, "NewerOptions.xml" ), newer_options );

  map<string,UserLibrary::Status> status;
  for( const UserLibrary::Item &item : user_items( "rel_act_preset" ) )
    status[item.names.at(0)] = UserLibrary::validateItem( "rel_act_preset", item.found ).status;

  BOOST_CHECK( status["Fine"] == UserLibrary::Status::Ok );
  BOOST_CHECK( status["Garbage"] == UserLibrary::Status::Unreadable );
  BOOST_CHECK( status["NewerTop"] == UserLibrary::Status::NeedsNewer );
  BOOST_CHECK( status["NewerOptions"] == UserLibrary::Status::NeedsNewer );

  // Templates are parsed with their includes; a missing include is unreadable.
  const string tmplts = SpecUtils::append_path( user, "report_templates/act_shield" );
  write_file( SpecUtils::append_path( tmplts, "good.tmplt.html" ), "Hi {% include \"part.html\" %}" );
  write_file( SpecUtils::append_path( tmplts, "part.html" ), "there {{ 1 }}" );
  write_file( SpecUtils::append_path( tmplts, "bad.tmplt.html" ), "Hi {% include \"missing.html\" %}" );
  const vector<UserLibrary::Item> templates = user_items( "act_shield_template" );
  BOOST_REQUIRE_EQUAL( templates.size(), 2 );
  for( const UserLibrary::Item &item : templates )
  {
    const UserLibrary::Validation v = UserLibrary::validateItem( item.typeId, item.found );
    if( item.names.at(0) == "good.tmplt.html" )
    {
      BOOST_CHECK_MESSAGE( v.status == UserLibrary::Status::Ok, v.message );
      BOOST_REQUIRE_EQUAL( item.found.support.size(), 1 );
      BOOST_CHECK_EQUAL( SpecUtils::filename( item.found.support[0] ), "part.html" );
    }else
    {
      BOOST_CHECK( v.status == UserLibrary::Status::Unreadable );
    }
  }

  // The old template folder is still read.
  write_file( SpecUtils::append_path( user, "act_shield_fit_templates/old.tmplt.txt" ), "old" );
  BOOST_CHECK_EQUAL( user_items( "act_shield_template" ).size(), 3 );

  std::filesystem::remove_all( fs_path( user ) );
}//BOOST_AUTO_TEST_CASE( validation_statuses )


BOOST_AUTO_TEST_CASE( create_install_update_uninstall )
{
  set_dirs();
  const string author_dir = make_temp_dir( "author" );
  const string recipient_dir = make_temp_dir( "recipient" );
  const string work = make_temp_dir( "work" );

  // The author's own files: a preset, a GADRAS detector, and a template with an include.
  use_user_dir( author_dir );
  write_file( SpecUtils::append_path( author_dir, "rel_act/Team Preset.xml" ), preset_text( "v1" ) );
  copy_folder( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/HPGe 40%" ),
               SpecUtils::append_path( author_dir, "drfs/Team HPGe" ) );
  write_file( SpecUtils::append_path( author_dir, "act_shield_fit_templates/team.tmplt.html" ),
              "Report {% include \"team_part.html\" %}" );
  write_file( SpecUtils::append_path( author_dir, "act_shield_fit_templates/team_part.html" ), "part" );

  UserLibraryPackage::CreateRequest request;
  request.manifest.name = "Team Kit";
  request.manifest.version = "1";
  request.manifest.author = "Author A";
  request.manifest.organization = "Org";
  request.manifest.changes = "First version";
  for( const string type : { "rel_act_preset", "drf", "act_shield_template" } )
  {
    for( const UserLibrary::Item &item : user_items( type ) )
    {
      UserLibraryPackage::CreateItem ci;
      ci.item = item;
      ci.notes = "Notes for " + item.relPath;
      request.items.push_back( ci );
    }
  }
  BOOST_REQUIRE_EQUAL( request.items.size(), 3 );

  for( const UserLibraryPackage::PreflightIssue &issue : UserLibraryPackage::preflight( request ) )
    BOOST_ERROR( "Preflight: " << issue.message );

  const string v1_path = SpecUtils::append_path( work, "Team Kit v1.ispkg" );
  const string package_id = create_package( request, v1_path );
  BOOST_CHECK( UserLibraryPackage::looksLikePackage( "Team Kit v1.ispkg", v1_path ) );

  // The author's machine records the package as created, not installed.
  BOOST_CHECK( UserLibrary::installedPackages().empty() );
  BOOST_REQUIRE_EQUAL( UserLibrary::createdPackages().count( package_id ), 1 );
  BOOST_CHECK_EQUAL( UserLibrary::lastCreator().author, "Author A" );

  // A preflight with a duplicate name is refused.
  {
    UserLibraryPackage::CreateRequest dup = request;
    dup.items.push_back( dup.items.front() );
    bool blocking = false;
    for( const UserLibraryPackage::PreflightIssue &issue : UserLibraryPackage::preflight( dup ) )
      blocking = blocking || issue.blocking;
    BOOST_CHECK( blocking );
  }

  // The recipient already has their own "Team Preset", and installs the package.
  use_user_dir( recipient_dir );
  write_file( SpecUtils::append_path( recipient_dir, "rel_act/Team Preset.xml" ), preset_text( "recipient's own" ) );

  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( v1_path, "Team Kit v1.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_CHECK( insp->hadManifest );
  BOOST_CHECK( insp->relation == UserLibraryPackage::Relation::NewPackage );
  BOOST_CHECK_EQUAL( insp->manifest.name, "Team Kit" );
  BOOST_REQUIRE_EQUAL( insp->items.size(), 3 );
  BOOST_CHECK( insp->skipped.empty() );

  const UserLibraryPackage::IncomingItem *preset = incoming( *insp, "rel_act/Team Preset.xml" );
  BOOST_REQUIRE( preset );
  BOOST_CHECK( preset->state == UserLibraryPackage::ItemState::Conflict );
  BOOST_REQUIRE_EQUAL( preset->conflicts.size(), 1 );
  BOOST_CHECK( preset->conflicts[0].existing.source.kind == UserLibrary::SourceKind::User );
  BOOST_CHECK( !preset->conflicts[0].useIncomingByDefault );
  BOOST_CHECK_EQUAL( preset->notes, "Notes for rel_act/Team Preset.xml" );

  const UserLibraryPackage::IncomingItem *drf = incoming( *insp, "drfs/Team HPGe" );
  BOOST_REQUIRE( drf );
  BOOST_CHECK( drf->found.isFolder );
  BOOST_CHECK_MESSAGE( drf->validation.status == UserLibrary::Status::Ok, drf->validation.message );
  BOOST_CHECK( drf->state == UserLibraryPackage::ItemState::New );

  const UserLibraryPackage::IncomingItem *tmplt = incoming( *insp, "report_templates/act_shield/team.tmplt.html" );
  BOOST_REQUIRE( tmplt );
  BOOST_CHECK_MESSAGE( tmplt->validation.status == UserLibrary::Status::Ok, tmplt->validation.message );

  // Keep the recipient's own preset (the default).
  const UserLibrary::Package installed = UserLibraryPackage::install( *insp, {} );
  BOOST_CHECK( insp->stagingDir.empty() );
  BOOST_CHECK( installed.recorded );
  BOOST_CHECK_EQUAL( installed.folder, "Team Kit" );
  BOOST_CHECK( SpecUtils::is_file( SpecUtils::append_path( installed.path, "report_templates/act_shield/team_part.html" ) ) );

  // Both presets are listed; the user's own answers lookups.
  {
    const vector<UserLibrary::Item> presets = UserLibrary::items( "rel_act_preset" );
    size_t num_team = 0;
    for( const UserLibrary::Item &item : presets )
      num_team += (item.names.at(0) == "Team Preset");
    BOOST_CHECK_EQUAL( num_team, 2 );

    const optional<UserLibrary::Item> chosen = UserLibrary::resolve( "rel_act_preset", "team preset" );
    BOOST_REQUIRE( chosen );
    BOOST_CHECK( chosen->source.kind == UserLibrary::SourceKind::User );
  }

  // The detector is in the DRF library, under the package's name.
  {
    bool found = false;
    for( const DrfLibrary::Location &loc : DrfLibrary::standardLocations() )
      found = found || ((loc.role == DrfLibrary::LocationRole::Package) && (loc.label == "Team Kit"));
    BOOST_CHECK( found );
    const optional<UserLibrary::Item> det = UserLibrary::resolve( "drf", "Team HPGe" );
    BOOST_REQUIRE( det );
    BOOST_CHECK( det->source.kind == UserLibrary::SourceKind::Package );
  }

  // Choose the package's preset instead.
  {
    vector<UserLibrary::Item> presets = UserLibrary::items( "rel_act_preset" );
    const UserLibrary::Item *pkg_preset = nullptr;
    for( const UserLibrary::Item &item : presets )
      if( (item.source.kind == UserLibrary::SourceKind::Package) && (item.names.at(0) == "Team Preset") )
        pkg_preset = &item;
    BOOST_REQUIRE( pkg_preset );
    UserLibrary::setChoice( "rel_act_preset", "Team Preset", *pkg_preset );
    const optional<UserLibrary::Item> chosen = UserLibrary::resolve( "rel_act_preset", "Team Preset" );
    BOOST_REQUIRE( chosen );
    BOOST_CHECK( chosen->source.kind == UserLibrary::SourceKind::Package );
  }

  // Version 2 from the author: preset changed, template removed, a new preset added.
  use_user_dir( author_dir );
  write_file( SpecUtils::append_path( author_dir, "rel_act/Team Preset.xml" ), preset_text( "v2" ) );
  write_file( SpecUtils::append_path( author_dir, "rel_act/Second.xml" ), preset_text( "second" ) );
  {
    UserLibraryPackage::CreateRequest v2;
    v2.manifest = UserLibrary::createdPackages().at( package_id ).back().manifest;
    v2.manifest.version = UserLibraryPackage::nextVersion( v2.manifest.version );
    v2.manifest.changes = "Changed the preset; added Second; dropped the template";
    v2.manifest.previousVersions.clear();
    for( const string type : { "rel_act_preset", "drf" } )
    {
      for( const UserLibrary::Item &item : user_items( type ) )
      {
        UserLibraryPackage::CreateItem ci;
        ci.item = item;
        v2.items.push_back( ci );
      }
    }
    BOOST_REQUIRE_EQUAL( v2.items.size(), 3 );
    const string id2 = create_package( v2, SpecUtils::append_path( work, "Team Kit v2.ispkg" ) );
    BOOST_CHECK_EQUAL( id2, package_id );
    BOOST_CHECK_EQUAL( UserLibrary::createdPackages().at( package_id ).size(), 2 );
  }

  use_user_dir( recipient_dir );
  insp = UserLibraryPackage::inspect( SpecUtils::append_path( work, "Team Kit v2.ispkg" ), "Team Kit v2.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_CHECK( insp->relation == UserLibraryPackage::Relation::Update );
  BOOST_REQUIRE( insp->installed );
  BOOST_CHECK_EQUAL( insp->added.size(), 1 );
  BOOST_CHECK_EQUAL( insp->changed.size(), 1 );
  BOOST_REQUIRE_EQUAL( insp->removed.size(), 1 );
  BOOST_CHECK_EQUAL( insp->removed[0], "report_templates/act_shield/team.tmplt.html" );
  BOOST_CHECK( insp->locallyModified.empty() );
  BOOST_CHECK( !insp->authorChanged );
  BOOST_REQUIRE_EQUAL( insp->changesSince.size(), 1 );
  BOOST_CHECK_EQUAL( insp->changesSince[0].version, "2" );

  preset = incoming( *insp, "rel_act/Team Preset.xml" );
  BOOST_REQUIRE( preset );
  BOOST_CHECK( preset->state == UserLibraryPackage::ItemState::Update );
  BOOST_REQUIRE_EQUAL( preset->conflicts.size(), 1 );
  BOOST_CHECK( preset->conflicts[0].alreadyChosen );
  BOOST_CHECK( preset->conflicts[0].useIncomingByDefault );

  drf = incoming( *insp, "drfs/Team HPGe" );
  BOOST_REQUIRE( drf );
  BOOST_CHECK( drf->state == UserLibraryPackage::ItemState::Unchanged );

  const UserLibrary::Package updated = UserLibraryPackage::install( *insp, {} );
  BOOST_CHECK_EQUAL( updated.folder, "Team Kit" );
  BOOST_CHECK_EQUAL( updated.manifest.version, "2" );
  BOOST_CHECK_EQUAL( updated.installs.size(), 2 );
  BOOST_CHECK( !SpecUtils::is_file( SpecUtils::append_path( updated.path, "report_templates/act_shield/team.tmplt.html" ) ) );

  // The choice survived the update.
  {
    const optional<UserLibrary::Item> chosen = UserLibrary::resolve( "rel_act_preset", "Team Preset" );
    BOOST_REQUIRE( chosen );
    BOOST_CHECK( chosen->source.kind == UserLibrary::SourceKind::Package );
    BOOST_CHECK( read_file( chosen->found.path ).find( "v2 " ) != string::npos );
  }

  // Disabling the package's preset hands lookups back to the user's own; the flag survives an update.
  {
    vector<UserLibrary::Item> presets = UserLibrary::items( "rel_act_preset" );
    for( const UserLibrary::Item &item : presets )
    {
      if( (item.source.kind == UserLibrary::SourceKind::Package) && (item.names.at(0) == "Second") )
        UserLibrary::setDisabled( item, true );
    }
    BOOST_CHECK( !UserLibrary::resolve( "rel_act_preset", "Second" ) );
    size_t num_disabled = 0;
    for( const UserLibrary::Item &item : UserLibrary::items( "rel_act_preset", true ) )
      num_disabled += item.disabled;
    BOOST_CHECK_EQUAL( num_disabled, 1 );

    // Re-install the same version: still disabled.
    insp = UserLibraryPackage::inspect( SpecUtils::append_path( work, "Team Kit v2.ispkg" ), "Team Kit v2.ispkg" );
    BOOST_REQUIRE( insp->error.empty() );
    BOOST_CHECK( insp->relation == UserLibraryPackage::Relation::Reinstall );
    UserLibraryPackage::install( *insp, {} );
    BOOST_CHECK( !UserLibrary::resolve( "rel_act_preset", "Second" ) );
  }

  // A file edited inside the installed package is noticed by the next update.
  write_file( SpecUtils::append_path( updated.path, "rel_act/Second.xml" ), preset_text( "edited" ) );
  insp = UserLibraryPackage::inspect( SpecUtils::append_path( work, "Team Kit v2.ispkg" ), "Team Kit v2.ispkg" );
  BOOST_REQUIRE( insp->error.empty() );
  BOOST_CHECK_EQUAL( insp->locallyModified.size(), 1 );
  UserLibraryPackage::discard( *insp );
  BOOST_CHECK( insp->stagingDir.empty() );

  // Uninstalling brings back the user's own preset, and forgets the choices that named the package.
  UserLibraryPackage::uninstall( package_id );
  BOOST_CHECK( UserLibrary::installedPackages().empty() );
  BOOST_CHECK( !SpecUtils::is_directory( updated.path ) );
  {
    const optional<UserLibrary::Item> chosen = UserLibrary::resolve( "rel_act_preset", "Team Preset" );
    BOOST_REQUIRE( chosen );
    BOOST_CHECK( chosen->source.kind == UserLibrary::SourceKind::User );
    BOOST_CHECK( UserLibrary::choices().empty() );
    BOOST_CHECK( !UserLibrary::resolve( "drf", "Team HPGe" ) );
  }

  std::filesystem::remove_all( fs_path( author_dir ) );
  std::filesystem::remove_all( fs_path( recipient_dir ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( create_install_update_uninstall )


BOOST_AUTO_TEST_CASE( package_without_manifest_and_newer_items )
{
  set_dirs();
  const string user = make_temp_dir( "nomanifest" );
  const string work = make_temp_dir( "nomanifest_work" );
  use_user_dir( user );

  // A plain zip, laid out like the user data folder, but without a manifest: not recognized as a
  //  package by name, but installable as one when given to inspect().
  const string zip = make_zip( work, "Shared Presets.zip", { {"rel_act/Handmade.xml", preset_text( "hand" )} } );
  BOOST_CHECK( !UserLibraryPackage::looksLikePackage( "Shared Presets.zip", zip ) );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( zip, "Shared Presets.zip" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_CHECK( !insp->hadManifest );
  BOOST_CHECK_EQUAL( insp->manifest.name, "Shared Presets" );
  BOOST_REQUIRE_EQUAL( insp->manifest.items.size(), 1 );
  const UserLibrary::Package pkg = UserLibraryPackage::install( *insp, {} );
  BOOST_CHECK( SpecUtils::is_file( SpecUtils::append_path( pkg.path, UserLibrary::ManifestFileName ) ) );
  BOOST_REQUIRE( UserLibrary::resolve( "rel_act_preset", "Handmade" ) );

  // A package for a newer InterSpec: installed, but its items are inactive.
  const string newer_manifest = R"({"format_version":1,"id":"newer-pkg","name":"Newer","version":"1",
                                    "min_interspec_version":"999.0"})";
  const string newer = make_zip( work, "Newer.ispkg", { {"interspec-package.json", newer_manifest},
                                 {"rel_act/FromTheFuture.xml", preset_text( "future" )} } );
  insp = UserLibraryPackage::inspect( newer, "Newer.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_REQUIRE_EQUAL( insp->items.size(), 1 );
  BOOST_CHECK( insp->items[0].state == UserLibraryPackage::ItemState::NeedsNewer );
  UserLibraryPackage::install( *insp, {} );
  BOOST_CHECK( !UserLibrary::resolve( "rel_act_preset", "FromTheFuture" ) );
  bool listed_inactive = false;
  for( const UserLibrary::Item &item : UserLibrary::items( "rel_act_preset", true ) )
    listed_inactive = listed_inactive || (item.needsNewer && (item.names.at(0) == "FromTheFuture"));
  BOOST_CHECK( listed_inactive );

  // A package format this InterSpec does not know is refused outright.
  const string future_format = make_zip( work, "Format2.ispkg", { {"interspec-package.json",
    R"({"format_version":2,"id":"f2","name":"F2","version":"1"})"}, {"rel_act/X.xml", preset_text()} } );
  BOOST_CHECK( !UserLibraryPackage::inspect( future_format, "Format2.ispkg" )->error.empty() );

  // Leftover staging folders are cleaned up at startup - once over an hour old (another InterSpec
  //  may be using a newer one).  An installed version set aside by an update that never finished is
  //  put back if its package's folder is gone.
  const auto age = []( const string &path ){
    std::filesystem::last_write_time( fs_path( path ), std::filesystem::file_time_type::clock::now() - std::chrono::hours( 2 ) );
  };
  const string stale = SpecUtils::append_path( UserLibrary::packagesDirectory(), ".incoming-deadbeef" );
  const string fresh = SpecUtils::append_path( UserLibrary::packagesDirectory(), ".incoming-feedface" );
  std::filesystem::create_directories( fs_path( stale ) );
  std::filesystem::create_directories( fs_path( fresh ) );
  age( stale );

  const optional<UserLibrary::Package> future_pkg = UserLibrary::installedPackage( "newer-pkg" );
  BOOST_REQUIRE( future_pkg );
  const string set_aside = SpecUtils::append_path( UserLibrary::packagesDirectory(), ".old-0badf00d" );
  std::filesystem::rename( fs_path( future_pkg->path ), fs_path( set_aside ) );
  age( set_aside );

  UserLibrary::cleanupStaleFolders();
  BOOST_CHECK( !SpecUtils::is_directory( stale ) );
  BOOST_CHECK( SpecUtils::is_directory( fresh ) );
  BOOST_CHECK( !SpecUtils::is_directory( set_aside ) );
  BOOST_CHECK( SpecUtils::is_directory( future_pkg->path ) );

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( package_without_manifest_and_newer_items )


BOOST_AUTO_TEST_CASE( template_includes_confined )
{
  set_dirs();
  const string user = make_temp_dir( "confine" );
  const string work = make_temp_dir( "confine_work" );
  use_user_dir( user );

  write_file( SpecUtils::append_path( user, "secret.txt" ), "TOP SECRET" );

  // A package template reaching outside its folder is flagged, and can not read the file when rendered.
  const string manifest = R"({"format_version":1,"id":"confine-test","name":"Confine","version":"1"})";
  const string pkg = make_zip( work, "confine.ispkg", {
    {"interspec-package.json", manifest},
    {"report_templates/act_shield/evil.tmplt.html", "x {% include \"../../../secret.txt\" %}"},
    {"report_templates/act_shield/good.tmplt.html", "A {% include \"part.html\" %}"},
    {"report_templates/act_shield/part.html", "B"} } );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "confine.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  const UserLibraryPackage::IncomingItem *evil = incoming( *insp, "report_templates/act_shield/evil.tmplt.html" );
  const UserLibraryPackage::IncomingItem *good = incoming( *insp, "report_templates/act_shield/good.tmplt.html" );
  BOOST_REQUIRE( evil && good );
  BOOST_CHECK( evil->state == UserLibraryPackage::ItemState::Unreadable );
  BOOST_CHECK_MESSAGE( good->validation.status == UserLibrary::Status::Ok, good->validation.message );

  const UserLibrary::Package installed = UserLibraryPackage::install( *insp, {} );
  const string dir = SpecUtils::append_path( installed.path, "report_templates/act_shield" );
  {
    inja::Environment env = BatchInfoLog::get_default_inja_env( dir + "/" );
    UserLibrary::confineTemplateIncludes( env, SpecUtils::append_path( dir, "good.tmplt.html" ) );
    BOOST_CHECK_EQUAL( env.render_file( "good.tmplt.html", nlohmann::json::object() ), "A B" );
  }
  {
    inja::Environment env = BatchInfoLog::get_default_inja_env( dir + "/" );
    UserLibrary::confineTemplateIncludes( env, SpecUtils::append_path( dir, "evil.tmplt.html" ) );
    string out;
    try
    {
      out = env.render_file( "evil.tmplt.html", nlohmann::json::object() );
    }catch( std::exception & )
    {
    }
    BOOST_CHECK( out.find( "TOP SECRET" ) == string::npos );
  }

  // The user's own template that includes a file outside its folder can not be packaged.
  write_file( SpecUtils::append_path( user, "report_templates/shared.html" ), "shared" );
  write_file( SpecUtils::append_path( user, "report_templates/act_shield/mine.tmplt.html" ), "{% include \"../shared.html\" %}" );
  UserLibraryPackage::CreateRequest request;
  request.manifest.name = "Mine";
  request.manifest.version = "1";
  for( const UserLibrary::Item &item : user_items( "act_shield_template" ) )
  {
    UserLibraryPackage::CreateItem ci;
    ci.item = item;
    request.items.push_back( ci );
  }
  BOOST_REQUIRE_EQUAL( request.items.size(), 1 );
  bool blocking = false;
  for( const UserLibraryPackage::PreflightIssue &issue : UserLibraryPackage::preflight( request ) )
    blocking = blocking || issue.blocking;
  BOOST_CHECK( blocking );

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( template_includes_confined )


BOOST_AUTO_TEST_CASE( create_from_mixed_sources )
{
  set_dirs();
  const string work = make_temp_dir( "mixed_work" );
  const string maker1 = make_temp_dir( "mixed_maker1" ), maker2 = make_temp_dir( "mixed_maker2" );
  const string user = make_temp_dir( "mixed_user" ), recipient = make_temp_dir( "mixed_recipient" );

  auto make_package = [&]( const string &dir, const string &name, const string &file ) -> string {
    use_user_dir( dir );
    UserLibraryPackage::CreateRequest request;
    request.manifest.name = name;
    request.manifest.version = "1";
    for( const string type : { "drf", "rel_act_preset", "act_shield_template" } )
    {
      for( const UserLibrary::Item &item : user_items( type ) )
      {
        UserLibraryPackage::CreateItem ci;
        ci.item = item;
        request.items.push_back( ci );
      }
    }
    const string path = SpecUtils::append_path( work, file );
    create_package( request, path );
    return path;
  };

  // Two packages, from two other people.
  copy_folder( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/NaI 3x3" ),
               SpecUtils::append_path( maker1, "drfs/Kit1 NaI" ) );
  write_file( SpecUtils::append_path( maker1, "report_templates/act_shield/kit1.tmplt.html" ), "Kit1 {% include \"kit1_part.html\" %}" );
  write_file( SpecUtils::append_path( maker1, "report_templates/act_shield/kit1_part.html" ), "part" );
  const string pkg1 = make_package( maker1, "Kit One", "kit1.ispkg" );

  write_file( SpecUtils::append_path( maker2, "rel_act/Kit2 Preset.xml" ), preset_text( "kit2" ) );
  write_file( SpecUtils::append_path( maker2, "rel_act/Kit2 Other.xml" ), preset_text( "kit2 other" ) );
  const string pkg2 = make_package( maker2, "Kit Two", "kit2.ispkg" );

  // This user installs both, and has their own files in the various places InterSpec reads.
  use_user_dir( user );
  for( const string &pkg : { pkg1, pkg2 } )
  {
    shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, SpecUtils::filename( pkg ) );
    BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
    UserLibraryPackage::install( *insp, {} );
  }
  copy_folder( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/HPGe 40%" ),
               SpecUtils::append_path( user, "drfs/My HPGe" ) );
  copy_folder( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/LaBr 10%" ),
               SpecUtils::append_path( user, "Old Detectors/My LaBr" ) );  //a GADRAS folder elsewhere in the user data folder
  write_file( SpecUtils::append_path( user, "rel_act/Uranium/My Preset.xml" ), preset_text( "mine" ) );
  write_file( SpecUtils::append_path( user, "rel_act/Not Shared.xml" ), preset_text( "not shared" ) );
  write_file( SpecUtils::append_path( user, "act_shield_fit_templates/mine.tmplt.html" ), "Mine" );  //the older template folder

  // Pick a subset: some of their own files, plus one item from each package.
  const set<string> wanted{ "My HPGe", "My LaBr", "My Preset", "mine.tmplt.html", "Kit1 NaI", "Kit2 Preset" };
  UserLibraryPackage::CreateRequest request;
  request.manifest.name = "Picked";
  request.manifest.version = "1";
  for( const string type : { "drf", "rel_act_preset", "act_shield_template" } )
  {
    for( const UserLibrary::Item &item : UserLibrary::items( type ) )
    {
      if( (item.source.kind == UserLibrary::SourceKind::BuiltIn) || item.names.empty() || !wanted.count( item.names.front() ) )
        continue;
      UserLibraryPackage::CreateItem ci;
      ci.item = item;
      request.items.push_back( ci );
    }
  }
  BOOST_REQUIRE_EQUAL( request.items.size(), wanted.size() );
  for( const UserLibraryPackage::PreflightIssue &issue : UserLibraryPackage::preflight( request ) )
    BOOST_ERROR( "Preflight: " << issue.message );

  const string picked = SpecUtils::append_path( work, "picked.ispkg" );
  std::ofstream out( fs_path( picked ), ios::out | ios::binary | ios::trunc );
  const UserLibrary::Manifest m = UserLibraryPackage::create( request, out );
  out.close();

  map<string,const UserLibrary::ManifestItem *> by_path;
  for( const UserLibrary::ManifestItem &mi : m.items )
    by_path[mi.path] = &mi;
  for( const char * const path : { "drfs/My HPGe", "drfs/Old Detectors/My LaBr", "rel_act/Uranium/My Preset.xml",
                                   "report_templates/act_shield/mine.tmplt.html", "drfs/Kit1 NaI", "rel_act/Kit2 Preset.xml" } )
    BOOST_CHECK_MESSAGE( by_path.count( path ), "No '" << path << "' in the package" );

  // Items from the packages say where they came from; the user's own do not.
  BOOST_REQUIRE( by_path.count( "drfs/Kit1 NaI" ) && by_path.count( "rel_act/Kit2 Preset.xml" ) && by_path.count( "drfs/My HPGe" ) );
  BOOST_CHECK_EQUAL( by_path["drfs/Kit1 NaI"]->derivedFromName, "Kit One" );
  BOOST_CHECK_EQUAL( by_path["rel_act/Kit2 Preset.xml"]->derivedFromName, "Kit Two" );
  BOOST_CHECK( by_path["drfs/My HPGe"]->derivedFromId.empty() );

  // Someone else installs it: every item reads, nothing is skipped.
  use_user_dir( recipient );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( picked, "picked.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_CHECK_EQUAL( insp->items.size(), wanted.size() );
  BOOST_CHECK( insp->skipped.empty() );
  for( const UserLibraryPackage::IncomingItem &item : insp->items )
    BOOST_CHECK_MESSAGE( item.validation.status == UserLibrary::Status::Ok, item.relPath << ": " << item.validation.message );
  UserLibraryPackage::discard( *insp );

  for( const string &d : { work, maker1, maker2, user, recipient } )
    std::filesystem::remove_all( fs_path( d ) );
}//BOOST_AUTO_TEST_CASE( create_from_mixed_sources )


BOOST_AUTO_TEST_CASE( create_with_database_and_outside_items )
{
  set_dirs();
  const string user = make_temp_dir( "outside_user" ), recipient = make_temp_dir( "outside_recipient" );
  const string outside = make_temp_dir( "outside_folder" ), work = make_temp_dir( "outside_work" );
  use_user_dir( user );

  // A detector efficiency "saved in InterSpec" (in memory here, as if loaded from the database).
  auto db_drf = make_shared<DetectorPeakResponse>( "Saved Det", "" );
  db_drf->fromGadrasDirectory( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/HPGe 20%" ), false );
  BOOST_REQUIRE( db_drf->isValid() );
  db_drf->setName( "Saved Det" );

  // A GADRAS detector in a folder outside the user data folder (as added in the Locations page), and a
  //  reference-spectra set elsewhere.
  const string added_folder = SpecUtils::append_path( outside, "Team Detectors" );
  copy_folder( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/NaI 2x2" ),
               SpecUtils::append_path( added_folder, "NaI Field" ) );
  const string spec_dir = SpecUtils::append_path( g_data_dir, "reference_spectra/Common_Field_Nuclides/IdentiFINDER-R500-NaI" );
  const string ref_set = SpecUtils::append_path( outside, "My Ref Set" );
  std::filesystem::create_directories( fs_path( ref_set ) );
  vector<string> spec_files = SpecUtils::ls_files_in_directory( spec_dir );
  std::sort( begin(spec_files), end(spec_files) );
  for( size_t i = 0; (i < spec_files.size()) && (i < 2); ++i )
    std::filesystem::copy_file( fs_path( spec_files[i] ),
                                fs_path( SpecUtils::append_path( ref_set, SpecUtils::filename( spec_files[i] ) ) ) );

  UserLibraryPackage::CreateRequest request;
  request.manifest.name = "Outside";
  request.manifest.version = "1";

  UserLibraryPackage::CreateItem db_item;
  db_item.item.typeId = "drf";
  db_item.item.source.kind = UserLibrary::SourceKind::User;
  db_item.item.names.push_back( "Saved Det" );
  db_item.databaseDrf = db_drf;
  db_item.sourceKey = "database:1234";
  request.items.push_back( db_item );

  UserLibraryPackage::CreateItem det_item;
  det_item.item.typeId = "drf";
  det_item.item.source.kind = UserLibrary::SourceKind::User;
  det_item.item.found.path = SpecUtils::append_path( added_folder, "NaI Field" );
  det_item.item.found.isFolder = true;
  det_item.item.names.push_back( "NaI Field" );
  det_item.root = SpecUtils::parent_path( added_folder );
  det_item.sourceKey = "location:" + det_item.item.found.path;
  request.items.push_back( det_item );

  UserLibraryPackage::CreateItem ref_item;
  ref_item.item.typeId = "ref_spectra";
  ref_item.item.source.kind = UserLibrary::SourceKind::User;
  ref_item.item.found.path = ref_set;
  ref_item.item.found.isFolder = true;
  ref_item.item.names.push_back( "My Ref Set" );
  ref_item.root = SpecUtils::parent_path( ref_set );
  ref_item.sourceKey = "location:" + ref_set;
  request.items.push_back( ref_item );

  BOOST_CHECK_EQUAL( UserLibraryPackage::packagePathFor( db_item ), "drfs/Saved Det.drf.xml" );
  BOOST_CHECK_EQUAL( UserLibraryPackage::packagePathFor( det_item ), "drfs/Team Detectors/NaI Field" );
  BOOST_CHECK_EQUAL( UserLibraryPackage::packagePathFor( ref_item ), "reference_spectra/My Ref Set" );
  for( const UserLibraryPackage::PreflightIssue &issue : UserLibraryPackage::preflight( request ) )
    BOOST_ERROR( "Preflight: " << issue.message );

  const string pkg = SpecUtils::append_path( work, "outside.ispkg" );
  std::ofstream out( fs_path( pkg ), ios::out | ios::binary | ios::trunc );
  UserLibraryPackage::create( request, out );
  out.close();

  // The record of what was made keeps where each came from, for the package's next version.
  const vector<UserLibrary::CreatedVersion> made = UserLibrary::createdPackages().begin()->second;
  BOOST_REQUIRE_EQUAL( made.size(), 1 );
  set<string> sources;
  for( const pair<const string,string> &src : made.back().sources )
    sources.insert( src.second );
  BOOST_CHECK( sources.count( "database:1234" ) );
  BOOST_CHECK( sources.count( "location:" + ref_set ) );

  // Someone else installs it: the database detector arrives as a readable .drf.xml.
  use_user_dir( recipient );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "outside.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_REQUIRE_EQUAL( insp->items.size(), 3 );
  const UserLibraryPackage::IncomingItem *saved = incoming( *insp, "drfs/Saved Det.drf.xml" );
  BOOST_REQUIRE( saved );
  BOOST_CHECK_MESSAGE( saved->validation.status == UserLibrary::Status::Ok, saved->validation.message );
  BOOST_REQUIRE_EQUAL( saved->validation.names.size(), 1 );
  BOOST_CHECK_EQUAL( saved->validation.names[0], "Saved Det" );
  for( const UserLibraryPackage::IncomingItem &item : insp->items )
    BOOST_CHECK_MESSAGE( item.validation.status == UserLibrary::Status::Ok, item.relPath << ": " << item.validation.message );
  UserLibraryPackage::discard( *insp );

  for( const string &d : { user, recipient, outside, work } )
    std::filesystem::remove_all( fs_path( d ) );
}//BOOST_AUTO_TEST_CASE( create_with_database_and_outside_items )


BOOST_AUTO_TEST_CASE( act_shield_setups )
{
  set_dirs();
  const string user = make_temp_dir( "actshield_user" );
  use_user_dir( user );

  // A real setup, as the Activity/Shielding Fit tool writes it (embedded in a test spectrum file).
  const string n42 = read_file( SpecUtils::append_path( g_data_dir, "../target/testing/test_data/ActShieldBackDecay/exemplar.n42" ) );
  const size_t start = n42.find( "<ShieldingSourceFit" );
  const size_t end = n42.find( "</ShieldingSourceFit>", start );
  BOOST_REQUIRE( (start != string::npos) && (end != string::npos) );
  const string setup = n42.substr( start, end + 21 - start );

  const string dir = SpecUtils::append_path( user, "act_shield_setups" );
  write_file( SpecUtils::append_path( dir, "Steel Drum.act.shield.xml" ), setup );
  string newer = setup;
  SpecUtils::ireplace_all( newer, "<ShieldingSourceFit version=\"0.", "<ShieldingSourceFit version=\"9." );
  write_file( SpecUtils::append_path( dir, "From The Future.act.shield.xml" ), newer );
  write_file( SpecUtils::append_path( dir, "notes.xml" ), "<Notes>not a setup</Notes>" );

  map<string,UserLibrary::Status> status;
  for( const UserLibrary::Item &item : user_items( "act_shield_setup" ) )
    status[item.names.at(0)] = UserLibrary::validateItem( "act_shield_setup", item.found ).status;
  BOOST_CHECK_EQUAL( status.size(), 2 );  //notes.xml is not a setup
  BOOST_CHECK( status["Steel Drum"] == UserLibrary::Status::Ok );
  BOOST_CHECK( status["From The Future"] == UserLibrary::Status::NeedsNewer );

  // Packaged as "act_shield_setups/<name>.act.shield.xml".
  for( const UserLibrary::Item &item : user_items( "act_shield_setup" ) )
  {
    if( item.names.at(0) == "Steel Drum" )
      BOOST_CHECK_EQUAL( UserLibraryPackage::packagePathFor( item ), "act_shield_setups/Steel Drum.act.shield.xml" );
  }

  std::filesystem::remove_all( fs_path( user ) );
}//BOOST_AUTO_TEST_CASE( act_shield_setups )


BOOST_AUTO_TEST_CASE( batch_report_templates )
{
  set_dirs();
  const string user = make_temp_dir( "batch_user" ), work = make_temp_dir( "batch_work" );
  use_user_dir( user );
  write_file( SpecUtils::append_path( user, "secret.txt" ), "TOP SECRET" );

  // A package with batch templates of each kind; one peak-fit template uses an include beside it and
  //  one of InterSpec's own partials, another tries to read a file outside its folder.
  const string manifest = R"({"format_version":1,"id":"batch-tmplts","name":"Batch Templates","version":"1"})";
  const string pkg = make_zip( work, "batch.ispkg", {
    {"interspec-package.json", manifest},
    {"report_templates/peak_fit/peaks.tmplt.txt", "Peaks: {% include \"part.txt\" %}"},
    {"report_templates/peak_fit/part.txt", "[part]"},
    {"report_templates/peak_fit/evil.tmplt.txt", "{% include \"../../../../secret.txt\" %}"},
    {"report_templates/peak_fit_summary/sum.tmplt.csv", "summary"},
    {"report_templates/act_shield_summary/act_sum.tmplt.html", "act summary"},
    {"report_templates/isotopics_auto_summary/iso_sum.tmplt.html", "iso summary"} } );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "batch.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  map<string,string> type_of;
  for( const UserLibraryPackage::IncomingItem &item : insp->items )
    type_of[SpecUtils::filename( item.relPath )] = item.typeId;
  BOOST_CHECK_EQUAL( type_of["peaks.tmplt.txt"], "peak_fit_template" );
  BOOST_CHECK_EQUAL( type_of["sum.tmplt.csv"], "peak_fit_summary_template" );
  BOOST_CHECK_EQUAL( type_of["act_sum.tmplt.html"], "act_shield_summary_template" );
  BOOST_CHECK_EQUAL( type_of["iso_sum.tmplt.html"], "isotopics_summary_template" );
  const UserLibraryPackage::IncomingItem *evil = incoming( *insp, "report_templates/peak_fit/evil.tmplt.txt" );
  BOOST_REQUIRE( evil );
  BOOST_CHECK( evil->state == UserLibraryPackage::ItemState::Unreadable );

  const UserLibrary::Package installed = UserLibraryPackage::install( *insp, {} );
  BOOST_CHECK_EQUAL( UserLibrary::items( "peak_fit_template" ).size(), 2 );

  // Rendered by the batch code: the include beside it works; the outside file is not read.
  BatchPeak::BatchPeakFitOptions options;
  inja::Environment env = BatchInfoLog::get_default_inja_env( "" );
  const string dir = SpecUtils::append_path( installed.path, "report_templates/peak_fit" );
  const string rendered = BatchInfoLog::render_template( SpecUtils::append_path( dir, "peaks.tmplt.txt" ), env,
                              BatchInfoLog::TemplateRenderType::PeakFitIndividual, options, nlohmann::json::object() );
  BOOST_CHECK_EQUAL( rendered, "Peaks: [part]" );
  string evil_out;
  try
  {
    evil_out = BatchInfoLog::render_template( SpecUtils::append_path( dir, "evil.tmplt.txt" ), env,
                              BatchInfoLog::TemplateRenderType::PeakFitIndividual, options, nlohmann::json::object() );
  }catch( std::exception & )
  {
  }
  BOOST_CHECK( evil_out.find( "TOP SECRET" ) == string::npos );

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( batch_report_templates )


BOOST_AUTO_TEST_CASE( color_themes )
{
  set_dirs();
  const string user = make_temp_dir( "themes_user" ), work = make_temp_dir( "themes_work" );
  use_user_dir( user );

  ColorTheme theme;
  theme.theme_name = "Team Night";
  const string json = ColorTheme::toJson( theme );

  // A package with a theme file, and a package made with a theme "from the database" (in memory).
  const string manifest = R"({"format_version":1,"id":"themes","name":"Themes","version":"1"})";
  const string pkg = make_zip( work, "themes.ispkg", { {"interspec-package.json", manifest},
                                                       {"color_themes/team_night.json", json},
                                                       {"color_themes/broken.json", "{ not json"} } );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "themes.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  const UserLibraryPackage::IncomingItem *night = incoming( *insp, "color_themes/team_night.json" );
  const UserLibraryPackage::IncomingItem *broken = incoming( *insp, "color_themes/broken.json" );
  BOOST_REQUIRE( night && broken );
  BOOST_CHECK_EQUAL( night->typeId, "color_theme" );
  BOOST_REQUIRE_EQUAL( night->validation.names.size(), 1 );
  BOOST_CHECK_EQUAL( night->validation.names[0], "Team Night" );
  BOOST_CHECK( broken->state == UserLibraryPackage::ItemState::Unreadable );
  const UserLibrary::Package installed = UserLibraryPackage::install( *insp, {} );

  // The theme is used from its file, relative to the user data folder; disabling it hides it.
  const string rel = "packages/" + installed.folder + "/color_themes/team_night.json";
  shared_ptr<ColorTheme> from_file = ColorThemeWindow::themeFromFile( rel );
  BOOST_REQUIRE( from_file );
  BOOST_CHECK( from_file->theme_name == theme.theme_name );
  for( const UserLibrary::Item &item : UserLibrary::items( "color_theme" ) )
  {
    if( item.found.path.find( "team_night" ) != string::npos )
      UserLibrary::setDisabled( item, true );
  }
  BOOST_CHECK( !ColorThemeWindow::themeFromFile( rel ) );

  // A theme from the database goes into a package as generated JSON.
  UserLibraryPackage::CreateRequest request;
  request.manifest.name = "Mine";
  request.manifest.version = "1";
  UserLibraryPackage::CreateItem ci;
  ci.item.typeId = "color_theme";
  ci.item.source.kind = UserLibrary::SourceKind::User;
  ci.item.names.push_back( "My Theme" );
  ci.generatedData = make_shared<const string>( json );
  ci.generatedName = "My Theme.json";
  ci.sourceKey = "database-theme:7";
  request.items.push_back( ci );
  BOOST_CHECK_EQUAL( UserLibraryPackage::packagePathFor( ci ), "color_themes/My Theme.json" );
  for( const UserLibraryPackage::PreflightIssue &issue : UserLibraryPackage::preflight( request ) )
    BOOST_ERROR( "Preflight: " << issue.message );
  std::ostringstream out;
  const UserLibrary::Manifest m = UserLibraryPackage::create( request, out );
  BOOST_REQUIRE_EQUAL( m.items.size(), 1 );
  BOOST_CHECK_EQUAL( m.items[0].path, "color_themes/My Theme.json" );

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( color_themes )


BOOST_AUTO_TEST_CASE( ref_lines_and_materials )
{
  set_dirs();
  const string user = make_temp_dir( "fragments_user" ), work = make_temp_dir( "fragments_work" );
  use_user_dir( user );

  // The older single files in the user data folder are still read, as one item each.
  write_file( SpecUtils::append_path( user, "add_ref_line.xml" ),
              "<RefLineDefinitions><SourceLines name=\"Team Lines\"><Line>100 1.0</Line></SourceLines></RefLineDefinitions>" );
  // A material file holding just a couple of materials (no "Air" entry, as a full copy would have).
  write_file( SpecUtils::append_path( user, "materials/team.txt" ),
              "!Team materials\nTeam Steel\nA steel\n1 0 7.9\nFe\n100\n!\nTeam Lead\nSome lead\n1 0 11.35\nPb\n100\n" );

  map<string,vector<string>> names;
  for( const string type : { "ref_lines", "materials" } )
  {
    for( const UserLibrary::Item &item : user_items( type ) )
    {
      const UserLibrary::Validation v = UserLibrary::validateItem( type, item.found );
      BOOST_CHECK_MESSAGE( v.status == UserLibrary::Status::Ok, item.relPath << ": " << v.message );
      names[type].insert( end(names[type]), begin(v.names), end(v.names) );
    }
  }
  BOOST_REQUIRE_EQUAL( names["ref_lines"].size(), 1 );
  BOOST_CHECK_EQUAL( names["ref_lines"][0], "Team Lines" );
  BOOST_REQUIRE_EQUAL( names["materials"].size(), 2 );
  BOOST_CHECK( (names["materials"][0] == "Team Steel") || (names["materials"][1] == "Team Steel") );

  // A package defining the same reference-line source: the user's own answers by default; a choice
  //  switches it.  (InterSpec reads these when it starts.)
  const string manifest = R"({"format_version":1,"id":"lines-pkg","name":"Lines","version":"1"})";
  const string pkg = make_zip( work, "lines.ispkg", { {"interspec-package.json", manifest},
    {"ref_lines/team.xml", "<RefLineDefinitions><SourceLines name=\"Team Lines\"><Line>200 1.0</Line></SourceLines></RefLineDefinitions>"} } );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "lines.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_REQUIRE_EQUAL( insp->items.size(), 1 );
  BOOST_CHECK( insp->items[0].state == UserLibraryPackage::ItemState::Conflict );
  UserLibraryPackage::install( *insp, {} );

  optional<UserLibrary::Item> chosen = UserLibrary::resolve( "ref_lines", "Team Lines" );
  BOOST_REQUIRE( chosen );
  BOOST_CHECK( chosen->source.kind == UserLibrary::SourceKind::User );
  for( const UserLibrary::Item &item : UserLibrary::items( "ref_lines" ) )
  {
    if( item.source.kind == UserLibrary::SourceKind::Package )
      UserLibrary::setChoice( "ref_lines", "Team Lines", item );
  }
  chosen = UserLibrary::resolve( "ref_lines", "Team Lines" );
  BOOST_REQUIRE( chosen );
  BOOST_CHECK( chosen->source.kind == UserLibrary::SourceKind::Package );

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( ref_lines_and_materials )


BOOST_AUTO_TEST_CASE( source_libs_and_agents )
{
  set_dirs();
  const string user = make_temp_dir( "srclib_user" ), work = make_temp_dir( "srclib_work" );
  use_user_dir( user );

  // Check-source libraries: the older "Source.lib" in the user data folder, and `source_libs/*.lib`.
  write_file( SpecUtils::append_path( user, "Source.lib" ), "Cs137_ABC123 3.7E+04 01-Jan-2020 the old one\n" );
  write_file( SpecUtils::append_path( user, "source_libs/team.lib" ), "Co60_T1 1.0E+05 01-Jan-2021\n" );
  write_file( SpecUtils::append_path( user, "source_libs/notes.lib" ), "just some notes\n" );
  map<string,UserLibrary::Status> statuses;
  for( const UserLibrary::Item &item : user_items( "source_lib" ) )
    statuses[SpecUtils::filename( item.found.path )] = UserLibrary::validateItem( "source_lib", item.found ).status;
  BOOST_CHECK( statuses["Source.lib"] == UserLibrary::Status::Ok );
  BOOST_CHECK( statuses["team.lib"] == UserLibrary::Status::Ok );
  BOOST_CHECK( statuses["notes.lib"] == UserLibrary::Status::Unreadable );

  // A package's library named "Source.lib": found once, through the user library, and not at all when
  //  disabled (the older search through the user data folder skips packages).
  const string manifest = R"({"format_version":1,"id":"srclib-pkg","name":"Sources","version":"1"})";
  const string pkg = make_zip( work, "sources.ispkg", { {"interspec-package.json", manifest},
                               {"source_libs/Source.lib", "Ba133_P1 2.0E+05 01-Jan-2022\n"} } );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "sources.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_REQUIRE_EQUAL( insp->items.size(), 1 );
  UserLibraryPackage::install( *insp, {} );

  const auto source_names = [](){
    vector<string> names;
    for( const SrcLibLineInfo &src : SrcLibLineInfo::sources_in_all_libs() )
      names.push_back( src.m_source_name );
    return names;
  };
  // Which package a file is in (e.g. for a detector's review reminder).
  for( const UserLibrary::Item &item : UserLibrary::items( "source_lib" ) )
  {
    const optional<UserLibrary::Package> owner = UserLibrary::packageForPath( item.found.path );
    BOOST_CHECK_EQUAL( !!owner, (item.source.kind == UserLibrary::SourceKind::Package) );
    if( owner )
      BOOST_CHECK_EQUAL( owner->id, "srclib-pkg" );
  }
  BOOST_CHECK( !UserLibrary::packageForPath( UserLibrary::packagesDirectory() ) );

  // Review dates: only a well-formed date before today is due.
  UserLibrary::Manifest review;
  BOOST_CHECK( !review.reviewDue() );
  review.reviewBy = "2001-02-03";
  BOOST_CHECK( review.reviewDue() );
  review.reviewBy = "9999-12-31";
  BOOST_CHECK( !review.reviewDue() );
  review.reviewBy = "next spring";
  BOOST_CHECK( !review.reviewDue() );

  vector<string> names = source_names();
  BOOST_CHECK_EQUAL( std::count( begin(names), end(names), "Ba133_P1" ), 1 );
  BOOST_CHECK_EQUAL( std::count( begin(names), end(names), "Co60_T1" ), 1 );
  BOOST_CHECK_EQUAL( std::count( begin(names), end(names), "Cs137_ABC123" ), 1 );
  for( const UserLibrary::Item &item : UserLibrary::items( "source_lib" ) )
  {
    if( item.source.kind == UserLibrary::SourceKind::Package )
      UserLibrary::setDisabled( item, true );
  }
  names = source_names();
  BOOST_CHECK_EQUAL( std::count( begin(names), end(names), "Ba133_P1" ), 0 );
  BOOST_CHECK_EQUAL( std::count( begin(names), end(names), "Co60_T1" ), 1 );

#if( USE_LLM_INTERFACE )
  // LLM agents: each agent comes from one file - the user's own over InterSpec's; a package's only if
  //  chosen, and then with InterSpec's as the fallback.
  const auto agent_xml = []( const string &name, const string &desc ){
    return "<AgentDefinitions version=\"0\"><Agent name=\"" + name + "\"><Description>" + desc
           + "</Description><SystemPrompt>Be careful.</SystemPrompt></Agent></AgentDefinitions>";
  };
  write_file( SpecUtils::append_path( user, "llm_agents/main.xml" ), agent_xml( "MainAgent", "Team main agent" ) );
  const vector<UserLibrary::Item> agents = user_items( "llm_agent" );
  BOOST_REQUIRE_EQUAL( agents.size(), 1 );
  const UserLibrary::Validation v = UserLibrary::validateItem( "llm_agent", agents[0].found );
  BOOST_CHECK( v.status == UserLibrary::Status::Ok );
  BOOST_REQUIRE_EQUAL( v.names.size(), 1 );
  BOOST_CHECK_EQUAL( v.names[0], "MainAgent" );

  vector<string> files = LlmConfig::agentFiles( AgentType::MainAgent );
  BOOST_REQUIRE_EQUAL( files.size(), 1 );
  BOOST_CHECK_EQUAL( files[0], agents[0].found.path );
  files = LlmConfig::agentFiles( AgentType::NuclideId );
  BOOST_REQUIRE_EQUAL( files.size(), 1 );
  BOOST_CHECK_EQUAL( SpecUtils::filename( files[0] ), "llm_agent_NuclideId.xml" );
  const string builtin_nuc_id = files[0];

  // An agent with a name InterSpec doesn't have, or a too-new file, isn't installed for use.
  const string agent_manifest = R"({"format_version":1,"id":"agent-pkg","name":"Agents","version":"1"})";
  const string agent_pkg = make_zip( work, "agents.ispkg", { {"interspec-package.json", agent_manifest},
    {"llm_agents/nuc.xml", agent_xml( "NuclideId", "Team nuclide ID" )},
    {"llm_agents/odd.xml", agent_xml( "NotAnAgent", "?" )},
    {"llm_agents/newer.xml", "<AgentDefinitions version=\"99\"><Agent name=\"InSitu\"/></AgentDefinitions>"} } );
  insp = UserLibraryPackage::inspect( agent_pkg, "agents.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  map<string,UserLibraryPackage::ItemState> states;
  for( const UserLibraryPackage::IncomingItem &item : insp->items )
    states[SpecUtils::filename( item.found.path )] = item.state;
  BOOST_CHECK( states["nuc.xml"] == UserLibraryPackage::ItemState::Conflict );
  BOOST_CHECK( states["odd.xml"] == UserLibraryPackage::ItemState::Unreadable );
  BOOST_CHECK( states["newer.xml"] == UserLibraryPackage::ItemState::NeedsNewer );
  UserLibraryPackage::install( *insp, {} );

  files = LlmConfig::agentFiles( AgentType::NuclideId );
  BOOST_REQUIRE_EQUAL( files.size(), 1 );
  BOOST_CHECK_EQUAL( files[0], builtin_nuc_id );  //the package's is not used unless chosen
  for( const UserLibrary::Item &item : UserLibrary::items( "llm_agent" ) )
  {
    if( (item.source.kind == UserLibrary::SourceKind::Package) && (SpecUtils::filename( item.found.path ) == "nuc.xml") )
      UserLibrary::setChoice( "llm_agent", "NuclideId", item );
  }
  files = LlmConfig::agentFiles( AgentType::NuclideId );
  BOOST_REQUIRE_EQUAL( files.size(), 2 );
  BOOST_CHECK_EQUAL( SpecUtils::filename( files[0] ), "nuc.xml" );
  BOOST_CHECK_EQUAL( files[1], builtin_nuc_id );
  files = LlmConfig::agentFiles( AgentType::InSitu );
  BOOST_REQUIRE_EQUAL( files.size(), 1 );
  BOOST_CHECK_EQUAL( SpecUtils::filename( files[0] ), "llm_agent_InSitu.xml" );  //the too-new one is not used
#endif

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( source_libs_and_agents )


BOOST_AUTO_TEST_CASE( nuclide_search_categories )
{
  typedef IsotopeSearchByEnergy::NucSearchCategory Category;

  set_dirs();
  const string user = make_temp_dir( "nsc_user" ), work = make_temp_dir( "nsc_work" );
  use_user_dir( user );

  const auto category = []( const string &name, const string &nuclide ){
    return "<NucSearchCategory><Name>" + name + "</Name><Description>d</Description>"
           "<AllowNucGammas>1</AllowNucGammas><Nuclides><Nuclide>" + nuclide + "</Nuclide></Nuclides>"
           "</NucSearchCategory>";
  };
  const auto categories_xml = []( const string &inner ){ return "<NucSearchCategories>" + inner + "</NucSearchCategories>"; };

  // The older single file in the user data folder (which InterSpec never actually read before), and
  //  a file in the folder.
  write_file( SpecUtils::append_path( user, "NuclideSearchCatagories.xml" ), categories_xml( category( "My Favs", "Co60" ) ) );
  write_file( SpecUtils::append_path( user, "nuclide_search_categories/team.xml" ),
              categories_xml( category( "Team Sources", "Cs137" ) + category( "isbe-category-norm", "K40" ) ) );
  for( const UserLibrary::Item &item : user_items( "nuc_search_category" ) )
  {
    const UserLibrary::Validation v = UserLibrary::validateItem( "nuc_search_category", item.found );
    BOOST_CHECK_MESSAGE( v.status == UserLibrary::Status::Ok, item.relPath << ": " << v.message );
    BOOST_CHECK( !v.names.empty() );
  }

  const auto load = [](){
    vector<Category> cats;
    IsotopeSearchByEnergy::init_category_info( cats );
    return cats;
  };
  const auto index_of = []( const vector<Category> &cats, const string &key ) -> int {
    for( size_t i = 0; i < cats.size(); ++i )
    {
      if( cats[i].key() == key )
        return static_cast<int>( i );
    }
    return -1;
  };
  const auto only_nuclide = []( const Category &cat ) -> string {
    return (cat.m_specific_nuclides.size() == 1) ? cat.m_specific_nuclides[0]->symbol : string();
  };

  vector<Category> cats = load();
  const int medical = index_of( cats, "isbe-category-medical" );
  BOOST_REQUIRE( medical > 0 );
  BOOST_REQUIRE( index_of( cats, "My Favs" ) > medical );    //after InterSpec's own
  BOOST_REQUIRE( index_of( cats, "Team Sources" ) > medical );
  BOOST_CHECK_EQUAL( only_nuclide( cats[index_of( cats, "Team Sources" )] ), "Cs137" );
  BOOST_CHECK( only_nuclide( cats[medical] ).empty() );
  const int norm = index_of( cats, "isbe-category-norm" );
  BOOST_REQUIRE( norm > 0 );
  BOOST_CHECK_EQUAL( only_nuclide( cats[norm] ), "K40" );  //the user's own replaces InterSpec's
  const size_t num_before = cats.size();

  // A package with a "Medical" of its own, a "Team Sources" of its own, and a new category: the names
  //  in use are conflicts, kept as they were unless chosen.
  const string manifest = R"({"format_version":1,"id":"nsc-pkg","name":"Search","version":"1"})";
  const string pkg = make_zip( work, "search.ispkg", { {"interspec-package.json", manifest},
    {"nuclide_search_categories/pkg.xml", categories_xml( category( "isbe-category-medical", "Tc99m" )
                                                          + category( "Team Sources", "Ba133" )
                                                          + category( "isbe-category-norm", "Th232" )
                                                          + category( "Pkg Only", "Am241" ) )} } );
  shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "search.ispkg" );
  BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
  BOOST_REQUIRE_EQUAL( insp->items.size(), 1 );
  BOOST_CHECK( insp->items[0].state == UserLibraryPackage::ItemState::Conflict );
  BOOST_CHECK_EQUAL( insp->items[0].conflicts.size(), 3 );  //one per name: against the one in use
  UserLibraryPackage::install( *insp, {} );

  cats = load();
  BOOST_CHECK_EQUAL( cats.size(), num_before + 1 );
  // Kept: the user's own, not InterSpec's (which also has the name).
  BOOST_CHECK_EQUAL( only_nuclide( cats[index_of( cats, "isbe-category-norm" )] ), "K40" );
  BOOST_REQUIRE( index_of( cats, "Pkg Only" ) > medical );
  BOOST_CHECK_EQUAL( only_nuclide( cats[index_of( cats, "Team Sources" )] ), "Cs137" );
  BOOST_CHECK( only_nuclide( cats[medical] ).empty() );

  // Choosing the package's "Medical": it replaces InterSpec's, in the same place.
  for( const UserLibrary::Item &item : UserLibrary::items( "nuc_search_category" ) )
  {
    if( item.source.kind == UserLibrary::SourceKind::Package )
      UserLibrary::setChoice( "nuc_search_category", "isbe-category-medical", item );
  }
  cats = load();
  BOOST_REQUIRE_EQUAL( index_of( cats, "isbe-category-medical" ), medical );
  BOOST_CHECK_EQUAL( only_nuclide( cats[medical] ), "Tc99m" );
  BOOST_CHECK_EQUAL( cats.size(), num_before + 1 );

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( nuclide_search_categories )


BOOST_AUTO_TEST_CASE( review_fixes )
{
  set_dirs();
  const string user = make_temp_dir( "fixes_user" ), work = make_temp_dir( "fixes_work" );
  use_user_dir( user );

  // Deleting a user's GADRAS detector removes just its files: not detectors in folders below it, and
  //  never the user data folder (here, a Detector.dat at its top), a type's folder, or the packages.
  {
    const string det = SpecUtils::append_path( user, "my_dets/Det A" );
    const string nested = SpecUtils::append_path( det, "Det B" );
    for( const string &d : { det, nested } )
    {
      write_file( SpecUtils::append_path( d, "Detector.dat" ), "x" );
      write_file( SpecUtils::append_path( d, "Efficiency.csv" ), "x" );
    }
    UserLibrary::Item item;
    item.typeId = "drf";
    item.source.kind = UserLibrary::SourceKind::User;
    item.found.path = det;
    item.found.isFolder = true;
    UserLibrary::deleteUserItem( item );
    BOOST_CHECK( !SpecUtils::is_file( SpecUtils::append_path( det, "Detector.dat" ) ) );
    BOOST_CHECK( SpecUtils::is_file( SpecUtils::append_path( nested, "Detector.dat" ) ) );

    item.found.path = user;  //a Detector.dat at the top of the user data folder
    BOOST_CHECK_THROW( UserLibrary::deleteUserItem( item ), std::exception );
    item.found.path = SpecUtils::append_path( user, "drfs" );
    BOOST_CHECK_THROW( UserLibrary::deleteUserItem( item ), std::exception );
    item.found.path = SpecUtils::parent_path( user );
    BOOST_CHECK_THROW( UserLibrary::deleteUserItem( item ), std::exception );
    item.source.kind = UserLibrary::SourceKind::Package;
    item.found.path = nested;
    BOOST_CHECK_THROW( UserLibrary::deleteUserItem( item ), std::exception );
    BOOST_CHECK( SpecUtils::is_directory( user ) );
  }

  // Colours from shared files can't carry script into CSS or JavaScript.
  BOOST_CHECK( ColorTheme::colorFromText( "red\"}]});alert(1);//" ).isDefault() );
  BOOST_CHECK( ColorTheme::colorFromText( "url(javascript:x)" ).isDefault() );
  BOOST_CHECK( !ColorTheme::colorFromText( "#aabbcc" ).isDefault() );
  BOOST_CHECK( !ColorTheme::colorFromText( "rgba(1, 2, 3, 0.5)" ).isDefault() );

  // Manifest item paths must be inside the package.
  for( const string bad : { "../x.xml", "/etc/x", "rel_act/../../x", "C:/x", "rel_act//x" } )
  {
    const string manifest = R"({"format_version":1,"id":"p","name":"P","version":"1","items":[{"path":")"
                            + bad + R"(","type":"rel_act_preset"}]})";
    BOOST_CHECK_THROW( UserLibrary::Manifest::fromJson( manifest ), std::exception );
  }

  // Hidden and very deep archive entries are not unpacked (just listed as ignored).
  {
    string deep = "rel_act";
    for( int i = 0; i < 20; ++i )
      deep += "/d";
    const string pkg = make_zip( work, "hidden.ispkg", {
      {"interspec-package.json", R"({"format_version":1,"id":"hidden-pkg","name":"Hidden","version":"1"})"},
      {"llm_knowledge/.sneaky/SKILL.md", "---\nname: sneaky\ndescription: x\n---\nx"},
      {deep + "/x.xml", "x"} } );
    shared_ptr<UserLibraryPackage::Inspection> insp = UserLibraryPackage::inspect( pkg, "hidden.ispkg" );
    BOOST_REQUIRE_MESSAGE( insp->error.empty(), insp->error );
    BOOST_CHECK_EQUAL( insp->skipped.size(), 2 );
    BOOST_CHECK( !SpecUtils::is_file( SpecUtils::append_path( insp->stagingDir, "llm_knowledge/.sneaky/SKILL.md" ) ) );
    UserLibraryPackage::discard( *insp );
  }

  // Materials and reference-line sources are the same name as their loaders see them.
  {
    const UserLibrary::AssetType * const materials = UserLibrary::assetType( "materials" );
    const UserLibrary::AssetType * const ref_lines = UserLibrary::assetType( "ref_lines" );
    BOOST_REQUIRE( materials && ref_lines );
    BOOST_CHECK( UserLibrary::sameName( *materials, "Stainless-Steel", "stainless steel" ) );
    BOOST_CHECK( UserLibrary::sameName( *ref_lines, "Co 60 Mix", "co60mix" ) );
    BOOST_CHECK( !UserLibrary::sameName( *ref_lines, "Co60 Mix", "Co57 Mix" ) );
  }

  // A history file that is not valid JSON is set aside (and a new one started); one that can't be
  //  read is left alone, and nothing is changed.
  {
    const string history = UserLibrary::historyPath();
    write_file( history, "{ not json" );
    UserLibrary::Item item;
    item.typeId = "rel_act_preset";
    item.source.kind = UserLibrary::SourceKind::User;
    item.found.path = SpecUtils::append_path( user, "rel_act/A.xml" );
    UserLibrary::setChoice( "rel_act_preset", "A", item );
    BOOST_CHECK( UserLibrary::choices().count( "rel_act_preset:A" ) );
    bool set_aside = false;
    for( const string &f : SpecUtils::ls_files_in_directory( user ) )
      set_aside = set_aside || SpecUtils::icontains( SpecUtils::filename( f ), ".unreadable-" );
    BOOST_CHECK( set_aside );

#ifndef _WIN32
    // Can't be read (no permission): left as it is, with the choice made above still in it.
    std::filesystem::permissions( fs_path( history ), std::filesystem::perms::none );
    BOOST_CHECK_THROW( UserLibrary::setChoice( "rel_act_preset", "B", item ), std::exception );
    std::filesystem::permissions( fs_path( history ), std::filesystem::perms::owner_read | std::filesystem::perms::owner_write );
    BOOST_CHECK( UserLibrary::choices().count( "rel_act_preset:A" ) );
    BOOST_CHECK( !UserLibrary::choices().count( "rel_act_preset:B" ) );
#endif
  }

  std::filesystem::remove_all( fs_path( user ) );
  std::filesystem::remove_all( fs_path( work ) );
}//BOOST_AUTO_TEST_CASE( review_fixes )
