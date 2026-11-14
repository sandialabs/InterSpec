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

/* Tests of `DrfLibrary` - the data behind the Detector Efficiency dialog's Library: the summary of
 each detector, the process-wide cache (reused unless a file or database row changed), merging the
 user's database rows into the detectors found on disk, and the rows the list shows for a filter.
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#define BOOST_TEST_MODULE TestDrfLibrary
#include <boost/test/included/unit_test.hpp>

#include <map>
#include <cmath>
#include <ctime>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <filesystem>

#include <Wt/Dbo/Dbo.h>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfLibrary.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace DrfLibrary;

namespace
{
string g_data_dir, g_test_file_dir;

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
    else if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }
  SpecUtils::ireplace_all( g_data_dir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

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

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory( SpecUtils::append_path(g_data_dir, "GenericGadrasDetectors") ),
                         "No data dir found; pass --datadir=..." );
  InterSpec::setStaticDataDirectory( g_data_dir );

  if( g_test_file_dir.empty() )
    g_test_file_dir = SpecUtils::append_path( g_data_dir, "../target/testing/test_data" );
}//set_dirs()


string generic_dir( const string &name )
{
  return SpecUtils::append_path( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" ), name );
}


shared_ptr<DetectorPeakResponse> load_generic( const string &name )
{
  auto drf = make_shared<DetectorPeakResponse>( name, "" );
  drf->fromGadrasDirectory( generic_dir( name ), false );
  return drf;
}


/** A fresh, empty temporary directory. */
string make_temp_dir( const string &tag )
{
  const string base = SpecUtils::temp_dir();
  const string dir = SpecUtils::append_path( base, "DrfLibraryTest_" + tag + "_" + std::to_string( std::time(nullptr) )
                                                   + "_" + std::to_string( rand() ) );
  std::filesystem::remove_all( dir );
  std::filesystem::create_directories( dir );
  return dir;
}


Location folder_location( const string &path, const bool built_in )
{
  Location loc;
  loc.id = loc.path = path;
  loc.isFolder = true;
  loc.builtIn = built_in;
  loc.isAuto = built_in;
  return loc;
}


Location file_location( const string &path, const bool built_in )
{
  Location loc = folder_location( path, built_in );
  loc.isFolder = false;
  return loc;
}


const Entry *find_entry( const vector<Entry> &entries, const string &name )
{
  for( const Entry &e : entries )
    if( e.summary && (e.summary->name == name) )
      return &e;
  return nullptr;
}


/** A synthetic entry for the merge/row tests. */
Entry synthetic( const string &id, const string &name, const uint64_t hash,
                 const uint64_t parent, const Origin origin, const int64_t last_used,
                 const bool geometry = false,
                 const DetectorPeakResponse::DrfSource source = DetectorPeakResponse::DrfSource::UnknownDrfSource,
                 const vector<uint64_t> &aliases = {} )
{
  auto s = make_shared<Summary>();
  s->name = name;
  s->hash = hash;
  s->source = source;
  s->sourceAliasHashes = aliases;
  s->parentHash = parent;
  s->hasEfficiency = true;
  s->hasGeometry = geometry;

  Entry e;
  e.id = id;
  e.origin = origin;
  e.summary = s;
  e.lastUsedUtc = last_used;
  e.knownHashes.push_back( hash );
  e.knownHashes.insert( end(e.knownHashes), begin(aliases), end(aliases) );
  return e;
}


/** The entry names of the detector rows, with "+" before nested copies; "[...]" for headers and
 "(N)"/"(N open)" for copies toggles - a compact picture of the list to compare against.
 */
vector<string> describe( const vector<Entry> &entries, const vector<Row> &rows )
{
  vector<string> answer;
  for( const Row &r : rows )
  {
    switch( r.type )
    {
      case RowType::GroupHeader:
      {
        const char *names[] = { "[Suggested]", "[Recent]", "[Everything]", "[All]" };
        answer.push_back( names[static_cast<int>(r.group)] );
        break;
      }
      case RowType::Detector:
        answer.push_back( string(r.nested ? "+" : "") + (r.copyOfHidden ? "copy:" : "")
                          + entries[r.entry].summary->name );
        break;
      case RowType::CopiesToggle:
        answer.push_back( "(" + std::to_string(r.numCopies) + (r.expanded ? " open)" : ")") );
        break;
    }
  }
  return answer;
}


void check_rows( const vector<string> &got, const vector<string> &expected )
{
  string got_str, exp_str;
  for( const string &s : got )
    got_str += s + ", ";
  for( const string &s : expected )
    exp_str += s + ", ";
  BOOST_CHECK_MESSAGE( got == expected, "Rows were {" << got_str << "} expected {" << exp_str << "}" );
}
}//namespace


BOOST_AUTO_TEST_CASE( ieee325_relative_efficiency )
{
  set_dirs();

  // A 3x3 NaI is the reference of the definition, so its value should be close to 100%; the generic
  //  GADRAS curves are not measured detectors, so allow a fair margin.
  const shared_ptr<DetectorPeakResponse> nai3x3 = load_generic( "NaI 3x3" );
  const double nai = nai3x3->relativeEfficiencyIeee325();
  BOOST_TEST_MESSAGE( "NaI 3x3 IEEE 325 rel. eff. = " << nai << "% (flat-disk model: "
                      << 100.0*nai3x3->flatDiskEfficiency( 1332.49f, 25.0*PhysicalUnits::cm )/1.2E-3 << "%)" );
  BOOST_CHECK_MESSAGE( (nai > 80.0) && (nai < 130.0), "NaI 3x3 rel. eff. " << nai << "% not near 100%" );

  const shared_ptr<DetectorPeakResponse> hpge40 = load_generic( "HPGe 40%" );
  const double hpge = hpge40->relativeEfficiencyIeee325();
  BOOST_TEST_MESSAGE( "HPGe 40% IEEE 325 rel. eff. = " << hpge << "% (flat-disk model: "
                      << 100.0*hpge40->flatDiskEfficiency( 1332.49f, 25.0*PhysicalUnits::cm )/1.2E-3 << "%)" );
  BOOST_CHECK_MESSAGE( (hpge > 25.0) && (hpge < 60.0), "HPGe 40% rel. eff. " << hpge << "% not near 40%" );

  // Bigger crystal, bigger number.
  const shared_ptr<DetectorPeakResponse> nai2x2 = load_generic( "NaI 2x2" );
  BOOST_CHECK_LT( nai2x2->relativeEfficiencyIeee325(), nai );

  // Not meaningful for a fixed-geometry DRF, or one with no efficiency.
  const shared_ptr<DetectorPeakResponse> fixed = nai3x3->reinterpretAsFixedGeom(
                                         DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  BOOST_REQUIRE( fixed );
  BOOST_CHECK_LT( fixed->relativeEfficiencyIeee325(), 0.0 );
  BOOST_CHECK_LT( DetectorPeakResponse().relativeEfficiencyIeee325(), 0.0 );
}//ieee325_relative_efficiency


BOOST_AUTO_TEST_CASE( hash_with_other_source )
{
  set_dirs();

  using DrfSource = DetectorPeakResponse::DrfSource;

  const shared_ptr<DetectorPeakResponse> a = load_generic( "NaI 3x3" );
  BOOST_CHECK_EQUAL( a->hashValueWithSource( a->drfSource() ), a->hashValue() );

  auto b = make_shared<DetectorPeakResponse>( *a );
  const DrfSource other = (a->drfSource() == DrfSource::DefaultGadrasDrf) ? DrfSource::UserAddedGadrasDrf
                                                                          : DrfSource::DefaultGadrasDrf;
  b->setDrfSource( other );
  BOOST_CHECK_NE( b->hashValue(), a->hashValue() );
  BOOST_CHECK_EQUAL( a->hashValueWithSource( other ), b->hashValue() );
  BOOST_CHECK_EQUAL( b->hashValueWithSource( a->drfSource() ), a->hashValue() );
}//hash_with_other_source


BOOST_AUTO_TEST_CASE( crystal_name_from_text )
{
  set_dirs();

  BOOST_CHECK_EQUAL( CeeLoUtils::crystalNameFromText( "Generic 3x3 NaI detector" ), "NaI" );
  BOOST_CHECK_EQUAL( CeeLoUtils::crystalNameFromText( "LaBr3 1.5x1.5" ), "LaBr3" );
  BOOST_CHECK_EQUAL( CeeLoUtils::crystalNameFromText( "60% HPGe coax" ), "HPGe" );
  BOOST_CHECK_EQUAL( CeeLoUtils::crystalNameFromText( "Since yesterday" ), "" );  //not "Si"
}//crystal_name_from_text


BOOST_AUTO_TEST_CASE( summaries_of_built_in_detectors )
{
  set_dirs();
  clearCache();

  ScanRequest request;
  Location generic = folder_location( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" ), true );
  generic.formats = FormatGadras;
  request.locations.push_back( generic );
  request.locations.push_back( file_location( SpecUtils::append_path( g_data_dir, "common_drfs.tsv" ), true ) );

  const Snapshot snap = scan( request );
  BOOST_REQUIRE_EQUAL( snap.locations.size(), 2 );
  BOOST_CHECK( snap.locations[0].found && snap.locations[1].found );
  BOOST_CHECK_GE( snap.locations[0].numDetectors, 10 );
  BOOST_CHECK_GE( snap.locations[1].numDetectors, 20 );
  BOOST_CHECK_EQUAL( snap.entries.size(), snap.locations[0].numDetectors + snap.locations[1].numDetectors );

  for( const Entry &e : snap.entries )
  {
    BOOST_REQUIRE( e.summary );
    BOOST_CHECK( e.origin == Origin::BuiltIn );
    BOOST_CHECK( !e.pending && !e.hidden );
    const DetectorPeakResponse::DrfSource expected = (e.kind == EntryKind::GadrasFolder)
                                        ? DetectorPeakResponse::DrfSource::DefaultGadrasDrf
                                        : DetectorPeakResponse::DrfSource::DefaultRelativeEfficiencyDrf;
    BOOST_CHECK_MESSAGE( e.summary->source == expected, "Wrong source for " << e.summary->name );
  }

  const Entry *nai = find_entry( snap.entries, "NaI 3x3" );
  BOOST_REQUIRE( nai );
  const Summary &s = *nai->summary;
  BOOST_CHECK( nai->kind == EntryKind::GadrasFolder );
  BOOST_CHECK_EQUAL( s.crystal.material, "NaI" );
  BOOST_CHECK( s.hasEfficiency && s.hasFwhm && s.hasGeometry );
  BOOST_CHECK( s.range == Range::NearFieldOk );  // GADRAS detectors get a curve-transfer response
  BOOST_CHECK_GT( s.relEffPercent, 80.0 );
  BOOST_CHECK( !s.nameStatesRelEff );
  BOOST_CHECK_GT( s.fwhmAt662, 20.0f );  //~7% at 662 keV

  const Entry *hpge = find_entry( snap.entries, "HPGe 40%" );
  BOOST_REQUIRE( hpge );
  BOOST_CHECK( hpge->summary->nameStatesRelEff );
  BOOST_CHECK_EQUAL( hpge->summary->crystal.material, "HPGe" );
  BOOST_CHECK_LT( hpge->summary->fwhmAt662, 3.0f );

  // Table rows: far-field coefficient curves.
  size_t num_rows = 0;
  for( const Entry &e : snap.entries )
  {
    if( e.kind != EntryKind::TableRow )
      continue;
    ++num_rows;
    BOOST_CHECK( e.summary->hasEfficiency );
    BOOST_CHECK( e.summary->range != Range::FixedGeometry );
    if( e.summary->efficiencyForm == DetectorPeakResponse::kExpOfLogPowerSeries )
      BOOST_CHECK_GT( e.summary->efficiencyParamCount, 0 );
    BOOST_CHECK_EQUAL( e.summary->sourceAliasHashes.size(), 1 );
    BOOST_CHECK( SpecUtils::istarts_with( e.id, "row:" ) );
  }
  BOOST_CHECK_GE( num_rows, 20 );

  // Every entry's full DRF loads, and is the detector summarized.
  for( const Entry *e : { nai, hpge, &snap.entries.back() } )
  {
    shared_ptr<const DetectorPeakResponse> drf;
    BOOST_REQUIRE_NO_THROW( drf = loadDrf( *e, nullptr ) );
    BOOST_REQUIRE( drf );
    BOOST_CHECK_EQUAL( drf->hashValue(), e->summary->hash );
  }
}//summaries_of_built_in_detectors


BOOST_AUTO_TEST_CASE( cache_reused_until_files_change )
{
  set_dirs();
  clearCache();

  // A user folder: two GADRAS detectors and a relative-efficiency table.
  const string dir = make_temp_dir( "cache" );
  for( const char * const name : { "NaI 3x3", "LaBr 10%" } )
  {
    const string dest = SpecUtils::append_path( dir, name );
    std::filesystem::create_directories( dest );
    for( const char * const file : { "Detector.dat", "Efficiency.csv" } )
      std::filesystem::copy_file( SpecUtils::append_path( generic_dir(name), file ),
                                  SpecUtils::append_path( dest, file ) );
  }

  const string table = SpecUtils::append_path( dir, "site_effs.tsv" );
  {
    ifstream input( SpecUtils::append_path( g_data_dir, "common_drfs.tsv" ).c_str() );
    ofstream output( table.c_str() );
    string line;
    size_t num_data_lines = 0;
    while( std::getline( input, line ) && (num_data_lines < 3) )
    {
      output << line << "\n";
      if( !line.empty() && (line[0] != '#') )
        ++num_data_lines;
    }
  }

  ScanRequest request;
  request.locations.push_back( folder_location( dir, false ) );

  const size_t computed_start = numSummariesComputed();
  const Snapshot first = scan( request );
  const size_t num_detectors = first.entries.size();
  BOOST_CHECK_GE( num_detectors, 4 );  // 2 GADRAS + the table's rows (3 lines, one may be a header)
  BOOST_CHECK_EQUAL( numSummariesComputed() - computed_start, num_detectors );

  for( const Entry &e : first.entries )
  {
    BOOST_REQUIRE( e.summary );
    BOOST_CHECK( e.origin == ((e.kind == EntryKind::GadrasFolder) ? Origin::Gadras : Origin::RelEffTable) );
    BOOST_CHECK( e.summary->source == ((e.kind == EntryKind::GadrasFolder)
                                       ? DetectorPeakResponse::DrfSource::UserAddedGadrasDrf
                                       : DetectorPeakResponse::DrfSource::UserAddedRelativeEfficiencyDrf) );
  }

  // Nothing changed: nothing is read again.
  const Snapshot second = scan( request );
  BOOST_CHECK_EQUAL( second.entries.size(), num_detectors );
  BOOST_CHECK_EQUAL( numSummariesComputed() - computed_start, num_detectors );

  // Change one detector's Efficiency.csv: only that detector is read again.
  {
    ofstream output( SpecUtils::append_path( SpecUtils::append_path( dir, "LaBr 10%" ), "Efficiency.csv" ).c_str(),
                     ios::app );
    output << "\n";
  }
  const Snapshot third = scan( request );
  BOOST_CHECK_EQUAL( third.entries.size(), num_detectors );
  BOOST_CHECK_EQUAL( numSummariesComputed() - computed_start, num_detectors + 1 );

  // A new detector folder appears; it is the only one read.
  {
    const string dest = SpecUtils::append_path( dir, "NaI 2x2" );
    std::filesystem::create_directories( dest );
    for( const char * const file : { "Detector.dat", "Efficiency.csv" } )
      std::filesystem::copy_file( SpecUtils::append_path( generic_dir("NaI 2x2"), file ),
                                  SpecUtils::append_path( dest, file ) );
  }
  const Snapshot fourth = scan( request );
  BOOST_CHECK_EQUAL( fourth.entries.size(), num_detectors + 1 );
  BOOST_CHECK_EQUAL( numSummariesComputed() - computed_start, num_detectors + 2 );

  // A folder holding only a Detector.dat is listed as needing characterization.
  std::filesystem::remove( SpecUtils::append_path( SpecUtils::append_path( dir, "NaI 2x2" ), "Efficiency.csv" ) );
  const Snapshot fifth = scan( request );
  const Entry *geom_only = nullptr;
  for( const Entry &e : fifth.entries )
    if( (e.kind == EntryKind::GadrasFolder) && SpecUtils::iends_with( e.path, "NaI 2x2" ) )
      geom_only = &e;
  BOOST_REQUIRE( geom_only && geom_only->summary );
  BOOST_CHECK( geom_only->summary->range == Range::NeedsCharacterization );
  BOOST_CHECK( !geom_only->summary->hasEfficiency );

  // A pending pass is reported when something needs reading.
  invalidatePath( dir );
  bool saw_progress = false;
  request.progress = [&]( const Snapshot &partial ){
    saw_progress = true;
    bool any_pending = false;
    for( const Entry &e : partial.entries )
      any_pending = any_pending || e.pending;
    BOOST_CHECK( any_pending );
  };
  const Snapshot sixth = scan( request );
  BOOST_CHECK( saw_progress );
  for( const Entry &e : sixth.entries )
    BOOST_CHECK( !e.pending );

  // Cancelling returns promptly with no entries.
  auto cancel = make_shared<std::atomic<bool>>( true );
  request.progress = nullptr;
  request.cancel = cancel;
  invalidatePath( dir );
  BOOST_CHECK( scan( request ).entries.empty() );

  std::filesystem::remove_all( dir );
}//cache_reused_until_files_change


BOOST_AUTO_TEST_CASE( isocs_folder_layout )
{
  set_dirs();
  clearCache();

  // Genie 2000 keeps .par grids in DCG sub-folders, and every detector's record in one
  //  PARAMETERS\DETECTOR.txt; .ecc results can be anywhere below isocs\data.
  const string det_eff = SpecUtils::append_path( g_test_file_dir, "det_eff" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path( det_eff, "mock.par" ) ),
                         "No test files found; pass --testfiledir=..." );

  const string root = make_temp_dir( "isocs" );
  const string data = SpecUtils::append_path( SpecUtils::append_path( root, "isocs" ), "data" );
  const string params = SpecUtils::append_path( data, "PARAMETERS" );
  const string dcg = SpecUtils::append_path( data, "DCG" );
  const string det1 = SpecUtils::append_path( dcg, "DET1" ), det2 = SpecUtils::append_path( dcg, "DET2" );
  const string geom = SpecUtils::append_path( SpecUtils::append_path( data, "GEOMETRY" ), "In-Situ" );
  for( const string &d : { params, det1, det2, geom } )
    std::filesystem::create_directories( d );

  const string grid = SpecUtils::append_path( det1, "mock.par" );
  std::filesystem::copy_file( SpecUtils::append_path( det_eff, "mock.par" ), grid );
  // A grid no record names: not listed, as guessing a record would show the wrong crystal.
  std::filesystem::copy_file( SpecUtils::append_path( det_eff, "mock.par" ),
                              SpecUtils::append_path( det2, "orphan.par" ) );
  std::filesystem::copy_file( SpecUtils::append_path( det_eff, "Detective-X_in-situ.ecc" ),
                              SpecUtils::append_path( geom, "Detective-X_in-situ.ecc" ) );

  // The mock record plus another detector's, so the record must be picked by the grid's name.
  const string detector_txt = SpecUtils::append_path( params, "DETECTOR.TXT" );
  {
    ifstream input( SpecUtils::append_path( det_eff, "mock_DETECTOR.txt" ).c_str(), ios::binary );
    ofstream output( detector_txt.c_str(), ios::binary );
    output << input.rdbuf()
           << "\r\nOtherDet,50.0,40.0,0,60.0,120.0,3.0,3.0,26,other.par,4, #\r\n"
              "ge,0.7,5.35, #\r\nal,0.5,2.7, #\r\n";
  }

  Location isocs;
  isocs.id = "ISOCS";
  isocs.path = data;
  isocs.isAuto = true;
  isocs.formats = FormatOther;

  ScanRequest request;
  request.locations.push_back( isocs );

  const Snapshot snap = scan( request );
  BOOST_REQUIRE_EQUAL( snap.locations.size(), 1u );
  BOOST_CHECK( snap.locations[0].found );
  BOOST_CHECK_EQUAL( snap.locations[0].numUnreadable, 1u );  // orphan.par
  BOOST_REQUIRE_EQUAL( snap.entries.size(), 2u );            // mock.par and the .ecc

  const Entry *par_entry = nullptr;
  for( const Entry &e : snap.entries )
  {
    BOOST_CHECK( e.origin == Origin::Isocs );
    if( e.path == grid )
      par_entry = &e;
  }
  BOOST_REQUIRE( par_entry && par_entry->summary );
  BOOST_CHECK_EQUAL( par_entry->companionPath, detector_txt );
  BOOST_CHECK( par_entry->summary->hasGeometry );
  BOOST_CHECK_CLOSE( par_entry->summary->crystal.diameter_cm, 5.0, 1.0E-3 );  // MockAngleCoax, not OtherDet

  // The grid is re-read when the DETECTOR.txt it uses changes.
  const size_t computed = numSummariesComputed();
  scan( request );
  BOOST_CHECK_EQUAL( numSummariesComputed(), computed );
  {
    ofstream output( detector_txt.c_str(), ios::binary | ios::app );
    output << "\r\n";
  }
  scan( request );
  BOOST_CHECK_EQUAL( numSummariesComputed(), computed + 1 );

  // A location of isocs\data\DCG still finds the PARAMETERS folder beside it.
  Location dcg_only = isocs;
  dcg_only.id = dcg_only.path = dcg;
  request.locations = { dcg_only };
  const Snapshot from_dcg = scan( request );
  BOOST_REQUIRE_EQUAL( from_dcg.entries.size(), 1u );
  BOOST_CHECK_EQUAL( from_dcg.entries[0].companionPath, detector_txt );

  std::filesystem::remove_all( root );
}//isocs_folder_layout


BOOST_AUTO_TEST_CASE( database_rows_merge_and_cache )
{
  set_dirs();
  clearCache();

  const string dir = make_temp_dir( "db" );
  DataBaseUtils::setPreferenceDatabaseFile( SpecUtils::append_path( dir, "test.db" ) );
  auto sql = make_shared<DataBaseUtils::DbSession>();

  long long user_id = -1;
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    Wt::Dbo::ptr<InterSpecUser> user = sql->session()->add(
                       std::make_unique<InterSpecUser>( "drf-library-test", InterSpecUser::Desktop ) );
    transaction.commit();
    user_id = user.id();
  }
  BOOST_REQUIRE_GE( user_id, 0 );

  auto add_row = [&]( const DetectorPeakResponse &drf ) -> long long {
    DataBaseUtils::DbTransaction transaction( *sql );
    auto row = std::make_unique<DetectorPeakResponse>( drf );
    row->m_user = static_cast<int>( user_id );
    Wt::Dbo::ptr<DetectorPeakResponse> ptr = sql->session()->add( std::move(row) );
    transaction.commit();
    return ptr.id();
  };

  // The built-in NaI 3x3, as the old GADRAS tab tagged it (UserAddedGadrasDrf) - must still fold into
  //  the built-in entry, which is now tagged DefaultGadrasDrf.
  const shared_ptr<DetectorPeakResponse> nai = load_generic( "NaI 3x3" );
  auto used_nai = make_shared<DetectorPeakResponse>( *nai );
  used_nai->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedGadrasDrf );
  used_nai->updateLastUsedTimeToNow();
  const long long nai_row = add_row( *used_nai );

  // A customized copy of it.
  auto copy = make_shared<DetectorPeakResponse>( *used_nai );
  copy->setName( "NaI 3x3, my FWHM" );
  copy->setParentHashValue( used_nai->hashValue() );
  add_row( *copy );

  // A usage record of the shipped NaI 2x2 from before InterSpec's modeling of it changed (stand-in:
  //  a different description, so its hash matches nothing) - folds into the shipped entry by name.
  auto stale = make_shared<DetectorPeakResponse>( *load_generic( "NaI 2x2" ) );
  stale->setDescription( "as modeled by an older InterSpec" );
  stale->setDrfSource( DetectorPeakResponse::DrfSource::DefaultGadrasDrf );
  const long long stale_row = add_row( *stale );

  // A formula detector.
  auto formula = make_shared<DetectorPeakResponse>( "Bench HPGe", "Entered as a formula" );
  formula->setIntrinsicEfficiencyFormula( "exp(-1.2 - 0.6*log(x/1000))", 6.0*PhysicalUnits::cm,
                                          PhysicalUnits::keV, 0.0f, 0.0f,
                                          DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  formula->setDrfSource( DetectorPeakResponse::DrfSource::UserSpecifiedFormulaDrf );
  const long long formula_row = add_row( *formula );

  ScanRequest request;
  Location generic = folder_location( SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" ), true );
  generic.formats = FormatGadras;
  request.locations.push_back( generic );
  request.sql = sql;
  request.dbUserId = user_id;

  const size_t computed_start = numSummariesComputed();
  const Snapshot snap = scan( request );

  const Entry *nai_entry = find_entry( snap.entries, "NaI 3x3" );
  BOOST_REQUIRE( nai_entry );
  BOOST_CHECK( nai_entry->kind == EntryKind::GadrasFolder );
  BOOST_CHECK_EQUAL( nai_entry->dbId, nai_row );
  BOOST_CHECK_GT( nai_entry->lastUsedUtc, 0 );
  BOOST_CHECK_EQUAL( nai_entry->knownHashes.size(), 2 );

  size_t num_nai = 0;
  for( const Entry &e : snap.entries )
    num_nai += (e.summary && (e.summary->name == "NaI 3x3"));
  BOOST_CHECK_EQUAL( num_nai, 1 );

  const Entry *nai2x2 = find_entry( snap.entries, "NaI 2x2" );
  BOOST_REQUIRE( nai2x2 );
  BOOST_CHECK( nai2x2->kind == EntryKind::GadrasFolder );
  BOOST_CHECK_EQUAL( nai2x2->dbId, stale_row );

  const Entry *copy_entry = find_entry( snap.entries, "NaI 3x3, my FWHM" );
  BOOST_REQUIRE( copy_entry );
  BOOST_CHECK( copy_entry->kind == EntryKind::DbRow );
  BOOST_CHECK( copy_entry->origin == Origin::UserMade );
  const int parent = parentIndex( snap.entries, static_cast<size_t>( copy_entry - &snap.entries[0] ) );
  BOOST_REQUIRE_GE( parent, 0 );
  BOOST_CHECK_EQUAL( snap.entries[static_cast<size_t>(parent)].id, nai_entry->id );

  const Entry *formula_entry = find_entry( snap.entries, "Bench HPGe" );
  BOOST_REQUIRE( formula_entry );
  BOOST_CHECK( formula_entry->origin == Origin::UserMade );
  BOOST_CHECK_EQUAL( formula_entry->dbId, formula_row );

  // The database row loads.
  shared_ptr<const DetectorPeakResponse> loaded;
  BOOST_REQUIRE_NO_THROW( loaded = loadDrf( *formula_entry, sql ) );
  BOOST_REQUIRE( loaded );
  BOOST_CHECK_EQUAL( loaded->hashValue(), formula->hashValue() );

  // Scanning again reads nothing.
  const size_t computed_after_first = numSummariesComputed();
  BOOST_CHECK_GT( computed_after_first, computed_start );
  const Snapshot again = scan( request );
  BOOST_CHECK_EQUAL( again.entries.size(), snap.entries.size() );
  BOOST_CHECK_EQUAL( numSummariesComputed(), computed_after_first );

  // Using a detector updates only its last-used time - no re-read - and it shows.
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    sql->session()->execute( "UPDATE DetectorPeakResponse SET m_lastUsedUtc = ? WHERE id = ?" )
                    .bind( static_cast<long long>(4102444800LL) ).bind( formula_row );
    transaction.commit();
  }
  const Snapshot used = scan( request );
  BOOST_CHECK_EQUAL( numSummariesComputed(), computed_after_first );
  BOOST_CHECK_EQUAL( find_entry( used.entries, "Bench HPGe" )->lastUsedUtc, 4102444800LL );

  // Changing a row's content (here, its appearance-free description) re-reads just that row.
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    Wt::Dbo::ptr<DetectorPeakResponse> row = sql->session()->find<DetectorPeakResponse>()
                                                .where( "id = ?" ).bind( formula_row ).resultValue();
    row.modify()->setDescription( "Now with a longer description than before" );
    transaction.commit();
  }
  const Snapshot changed = scan( request );
  BOOST_CHECK_EQUAL( numSummariesComputed(), computed_after_first + 1 );
  BOOST_CHECK_EQUAL( find_entry( changed.entries, "Bench HPGe" )->summary->description,
                     "Now with a longer description than before" );

  // A suggestion matching a listed detector marks it; one matching nothing is added.
  ScanRequest with_suggestions = request;
  with_suggestions.suggestions.push_back( Suggestion{ nai, SuggestDefaultForType } );
  auto from_file = make_shared<DetectorPeakResponse>( *formula );
  from_file->setName( "Factory efficiency" );
  with_suggestions.suggestions.push_back( Suggestion{ from_file, SuggestInSpectrumFile } );
  const Snapshot suggested = scan( with_suggestions );
  BOOST_CHECK_EQUAL( find_entry( suggested.entries, "NaI 3x3" )->suggestReasons, SuggestDefaultForType );
  const Entry *file_entry = find_entry( suggested.entries, "Factory efficiency" );
  BOOST_REQUIRE( file_entry );
  BOOST_CHECK( file_entry->kind == EntryKind::InMemory );
  BOOST_CHECK( file_entry->origin == Origin::ThisSpectrum );
  BOOST_CHECK_EQUAL( loadDrf( *file_entry, nullptr ).get(), from_file.get() );

  sql.reset();
}//database_rows_merge_and_cache


BOOST_AUTO_TEST_CASE( merge_entries )
{
  using DrfSource = DetectorPeakResponse::DrfSource;
  vector<Entry> located{
    synthetic( "gadras:a", "NaI 3x3", 1, 0, Origin::BuiltIn, 0, true, DrfSource::DefaultGadrasDrf, {11} ),
    synthetic( "gadras:b", "HPGe 40%", 2, 0, Origin::BuiltIn, 0, true, DrfSource::DefaultGadrasDrf, {4} )
  };
  vector<Entry> db{
    synthetic( "db:3", "HPGe 40%", 4, 0, Origin::UserMade, 900, false, DrfSource::UserAddedGadrasDrf ), // B, other source
    synthetic( "db:4", "HPGe 40%", 5, 0, Origin::UserMade, 100, false, DrfSource::DefaultGadrasDrf ),   // B, stale
    synthetic( "db:5", "Bench", 6, 0, Origin::UserMade, 50, false, DrfSource::UserSpecifiedFormulaDrf ),
    synthetic( "db:6", "HPGe 40%", 7, 0, Origin::UserMade, 10, false, DrfSource::UserImportedGadrasDrf ) // not B
  };
  db[0].dbId = 3;
  db[1].dbId = 4;
  db[2].dbId = 5;
  db[3].dbId = 6;
  for( Entry &e : located )
    e.kind = EntryKind::GadrasFolder;
  for( Entry &e : db )
    e.kind = EntryKind::DbRow;

  const vector<Entry> merged = mergeEntries( located, db, {} );
  BOOST_REQUIRE_EQUAL( merged.size(), 4 );
  BOOST_CHECK_EQUAL( merged[1].dbId, 3 );            // the newest matching row
  BOOST_CHECK_EQUAL( merged[1].lastUsedUtc, 900 );
  BOOST_CHECK_EQUAL( merged[1].knownHashes.size(), 3 ); // 2, alias 4, stale 5
  BOOST_CHECK_EQUAL( merged[2].summary->name, "Bench" );
  BOOST_CHECK( merged[2].origin == Origin::UserMade );
  BOOST_CHECK_EQUAL( merged[3].dbId, 6 );            // an import is never folded by name
  BOOST_CHECK( merged[3].origin == Origin::UserMade );
}//merge_entries


BOOST_AUTO_TEST_CASE( rows_groups_nesting_and_filters )
{
  vector<Entry> e{
    synthetic( "a", "NaI 3x3", 1, 0, Origin::BuiltIn, 0, true ),          // 0
    synthetic( "b", "HPGe 40%", 2, 0, Origin::BuiltIn, 900, true ),       // 1
    synthetic( "c", "NaI 3x3, my FWHM", 3, 1, Origin::UserMade, 500 ),    // 2: copy of a
    synthetic( "f", "Bench", 6, 0, Origin::UserMade, 0 ),                 // 3
    synthetic( "g", "NaI 3x3, jig", 7, 3, Origin::UserMade, 0 ),          // 4: copy of c
    synthetic( "h", "Zeta", 8, 0, Origin::Isocs, 0, true ),               // 5
    synthetic( "s", "Factory", 9, 0, Origin::ThisSpectrum, 0 )            // 6
  };
  e[6].suggestReasons = SuggestInSpectrumFile;

  Filter all;
  const map<string,bool> none;

  // Default: suggested, recent (family last-used: a's copy was used), then the rest by name.
  //  Both copies (including the copy of a copy) nest under "NaI 3x3", collapsed.
  check_rows( describe( e, buildRows( e, all, none, "", "" ) ),
              { "[Suggested]", "Factory",
                "[Recent]", "HPGe 40%", "NaI 3x3", "(2)",
                "[Everything]", "Bench", "Zeta" } );

  // Picking a copy (or it being in use) opens its family.
  check_rows( describe( e, buildRows( e, all, none, "c", "" ) ),
              { "[Suggested]", "Factory",
                "[Recent]", "HPGe 40%", "NaI 3x3", "(2 open)", "+NaI 3x3, jig", "+NaI 3x3, my FWHM",
                "[Everything]", "Bench", "Zeta" } );
  check_rows( describe( e, buildRows( e, all, none, "", "g" ) ),
              { "[Suggested]", "Factory",
                "[Recent]", "HPGe 40%", "NaI 3x3", "(2 open)", "+NaI 3x3, jig", "+NaI 3x3, my FWHM",
                "[Everything]", "Bench", "Zeta" } );

  // ... unless the user closed it.
  check_rows( describe( e, buildRows( e, all, { {"a", false} }, "c", "" ) ),
              { "[Suggested]", "Factory",
                "[Recent]", "HPGe 40%", "NaI 3x3", "(2)",
                "[Everything]", "Bench", "Zeta" } );

  // "Made by you": the originals are filtered out, so the copies are listed on their own.
  Filter mine;
  mine.allOrigins = false;
  mine.origin = Origin::UserMade;
  check_rows( describe( e, buildRows( e, mine, none, "", "" ) ),
              { "[Recent]", "copy:NaI 3x3, my FWHM",
                "[Everything]", "Bench", "copy:NaI 3x3, jig" } );

  // "With geometry" combines with the chip.
  Filter geom;
  geom.withGeometry = true;
  check_rows( describe( e, buildRows( e, geom, none, "", "" ) ),
              { "[Recent]", "HPGe 40%", "NaI 3x3", "(2)",
                "[Everything]", "Zeta" } );
  geom.allOrigins = false;
  geom.origin = Origin::Isocs;
  check_rows( describe( e, buildRows( e, geom, none, "", "" ) ), { "[All]", "Zeta" } );

  // Search: every word must match, case-insensitively; extra text (e.g. a translated label) counts.
  Filter search;
  search.query = "nai FWHM";
  check_rows( describe( e, buildRows( e, search, none, "", "" ) ), { "[Recent]", "copy:NaI 3x3, my FWHM" } );
  search.query = "isocs";
  check_rows( describe( e, buildRows( e, search, none, "", "" ) ), {} );
  const ExtraSearchText labels = []( const Entry &entry ) -> string {
    return (entry.origin == Origin::Isocs) ? "ISOCS" : "";
  };
  check_rows( describe( e, buildRows( e, search, none, "", "", labels ) ), { "[All]", "Zeta" } );

  // "Recently used" is capped; the overflow is listed with everything else.
  Filter capped;
  capped.maxRecent = 1;
  check_rows( describe( e, buildRows( e, capped, none, "", "" ) ),
              { "[Suggested]", "Factory",
                "[Recent]", "HPGe 40%",
                "[Everything]", "Bench", "NaI 3x3", "(2)", "Zeta" } );

  // Hidden entries (a hidden location) are left out, and their copies are listed on their own.
  vector<Entry> hidden = e;
  hidden[0].hidden = true;
  check_rows( describe( hidden, buildRows( hidden, all, none, "", "" ) ),
              { "[Suggested]", "Factory",
                "[Recent]", "HPGe 40%", "copy:NaI 3x3, my FWHM", "(1)",
                "[Everything]", "Bench", "Zeta" } );
}//rows_groups_nesting_and_filters



BOOST_AUTO_TEST_CASE( location_settings_round_trip )
{
  set_dirs();

  const string dir = make_temp_dir( "locs" );
  const string table = SpecUtils::append_path( dir, "site.tsv" );
  { ofstream out( table.c_str() ); out << "#nothing\n"; }

  // Entries as the preference stores them; blanks and repeats are tolerated
  const string pref = dir + ";  !" + table + " ;;gadras:C:\\GADRAS\\Old;@ISOCS=D:\\GENIE2K\\isocs\\data;!@ANGLE;"
                      + dir + ";!@";
  const LocationSettings settings = parseLocationSettings( pref );

  BOOST_REQUIRE_EQUAL( settings.userLocations.size(), 3u );
  BOOST_CHECK_EQUAL( settings.userLocations[0].path, dir );
  BOOST_CHECK( !settings.userLocations[0].hidden );
  BOOST_CHECK_EQUAL( settings.userLocations[0].formats, uint32_t(FormatAll) );
  BOOST_CHECK_EQUAL( settings.userLocations[1].path, table );
  BOOST_CHECK( settings.userLocations[1].hidden );
  BOOST_CHECK_EQUAL( settings.userLocations[2].path, "C:\\GADRAS\\Old" );
  BOOST_CHECK_EQUAL( settings.userLocations[2].formats, uint32_t(FormatGadras) );
  BOOST_CHECK( settings.hiddenStandard == set<string>{ "ANGLE" } );
  BOOST_REQUIRE_EQUAL( settings.standardPaths.size(), 1u );
  BOOST_CHECK_EQUAL( settings.standardPaths.at("ISOCS"), "D:\\GENIE2K\\isocs\\data" );

  // Writing and reading back gives the same settings
  const string written = serializeLocationSettings( settings );
  BOOST_CHECK( parseLocationSettings( written ) == settings );
  BOOST_CHECK( parseLocationSettings( "" ) == LocationSettings() );

  // Applied to standard locations: a hidden one and a located one, then the user's own; a user
  //  location that is a standard one is not listed twice.
  vector<Location> standard( 3 );
  standard[0].id = standard[0].path = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" );
  standard[1].id = "ISOCS";
  standard[1].path = "C:\\GENIE2K\\isocs\\data";
  standard[2].id = "ANGLE";
  standard[2].path = "C:\\Users\\x\\Documents\\ANGLE";
  LocationSettings with_dup = settings;
  Location dup;
  dup.id = dup.path = standard[0].path;
  with_dup.userLocations.push_back( dup );

  const vector<Location> applied = applyLocationSettings( standard, with_dup );
  BOOST_REQUIRE_EQUAL( applied.size(), 3u + 3u );
  BOOST_CHECK( !applied[0].hidden );
  BOOST_CHECK_EQUAL( applied[1].path, "D:\\GENIE2K\\isocs\\data" );
  BOOST_CHECK( applied[2].hidden );
  BOOST_CHECK( applied[3].isFolder );    //dir
  BOOST_CHECK( !applied[4].isFolder );   //site.tsv
  BOOST_CHECK( applied[4].role == LocationRole::User );

  std::filesystem::remove_all( dir );
}//location_settings_round_trip
