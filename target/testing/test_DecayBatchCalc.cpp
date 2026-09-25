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
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE test_DecayBatchCalc_suite
#include <boost/test/included/unit_test.hpp>

#include <map>
#include <set>
#include <cmath>
#include <regex>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayBatchCalc.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace DecayBatchCalc;

// Data directory globals (set from the command line).
std::string g_data_dir = "";
std::string g_test_data_dir = "";

namespace
{
  /** Reads the data directories from the command line and points the decay database at them.
   Runs once, on first use, and must happen before the first DecayDataBaseServer::database() call.
   */
  void init_data_dirs()
  {
    static bool s_have_init = false;
    if( s_have_init )
      return;
    s_have_init = true;

    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( arg.find("--datadir=") == 0 )
        g_data_dir = arg.substr( 10 );
      else if( arg.find("--testfiledir=") == 0 )
        g_test_data_dir = arg.substr( 14 );
    }

    SpecUtils::ireplace_all( g_data_dir, "%20", " " );
    SpecUtils::ireplace_all( g_test_data_dir, "%20", " " );

    // Search around a little for the data directories, if they were not given.
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

    if( g_test_data_dir.empty() )
    {
      for( const char * const d : { "target/testing/test_data", "../target/testing/test_data",
                                    "../../target/testing/test_data", "test_data", "../test_data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path(d, "BatchDecay/probeResults.csv") ) )
        {
          g_test_data_dir = d;
          break;
        }
      }
    }//if( g_test_data_dir.empty() )

    BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Could not find the data directory (--datadir=...)" );
    BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(),
                          "Could not find the test data directory (--testfiledir=...)" );

    const string decay_xml = SpecUtils::append_path( g_data_dir, "sandia.decay.xml" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(decay_xml), "No sandia.decay.xml at " + decay_xml );

    BOOST_REQUIRE_NO_THROW( DecayDataBaseServer::setDecayXmlFile( decay_xml ) );
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
  }//init_data_dirs()


  string batch_decay_dir()
  {
    init_data_dirs();
    return SpecUtils::append_path( g_test_data_dir, "BatchDecay" );
  }

  string read_file( const string &path )
  {
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(path), "Missing fixture: " + path );
    ifstream input( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE_MESSAGE( input.is_open(), "Could not open: " + path );
    stringstream contents;
    contents << input.rdbuf();
    return contents.str();
  }//read_file(...)


  /** The motivating multi-location file: 35 probe locations x 40 nuclides. */
  vector<BatchNuclide> probe_results()
  {
    return parse_csv( read_file( SpecUtils::append_path( batch_decay_dir(), "probeResults.csv" ) ) );
  }


  /** Decay options for a plain single-step run of `time_str`. */
  BatchDecayOptions options_for( const string &time_str )
  {
    BatchDecayOptions opts;
    opts.time_span = PhysicalUnits::stringToTimeDuration( time_str );
    opts.time_span_str = time_str;
    opts.num_steps = 1;
    opts.include_activity = true;
    return opts;
  }//options_for(...)


  /** Index of a column by header name; -1 when absent. */
  int column_index( const BatchDecayResult &result, const string &name )
  {
    for( size_t i = 0; i < result.column_headers.size(); ++i )
    {
      if( result.column_headers[i] == name )
        return static_cast<int>( i );
    }
    return -1;
  }//column_index(...)


  /** Rows of one location, as (nuclide symbol -> value) in row order. */
  vector<pair<string,double>> location_rows( const BatchDecayResult &result, const string &location )
  {
    const int name_col = column_index( result, "Probe name" );
    const int prod_col = column_index( result, "Product" );
    const int val_col = column_index( result, "Value" );
    BOOST_REQUIRE( (name_col >= 0) && (prod_col >= 0) && (val_col >= 0) );

    vector<pair<string,double>> answer;
    for( const vector<string> &row : result.rows )
    {
      if( row[name_col] != location )
        continue;

      const string &product = row[prod_col];
      const string::size_type sp = product.find( ' ' );
      answer.emplace_back( (sp == string::npos) ? product : product.substr(0,sp),
                           std::stod( row[val_col] ) );
    }//for( each row )

    return answer;
  }//location_rows(...)


  double value_for( const vector<pair<string,double>> &rows, const string &symbol )
  {
    for( const pair<string,double> &nv : rows )
    {
      if( nv.first == symbol )
        return nv.second;
    }
    BOOST_FAIL( "No row for nuclide " + symbol );
    return 0.0;
  }//value_for(...)


  /** The value of the first result row whose label starts with `symbol + " "`, at `step`. */
  double wide_value( const BatchDecayResult &result, const string &symbol, const size_t step )
  {
    for( const vector<string> &row : result.rows )
    {
      if( SpecUtils::starts_with( row.front(), (symbol + " ").c_str() ) )
      {
        BOOST_REQUIRE( (step + 1) < row.size() );
        return std::stod( row[step + 1] );
      }
    }
    BOOST_FAIL( "No row for " + symbol );
    return 0.0;
  }//wide_value(...)


  /** A synthetic multi-location CSV, for the key-rule fallbacks. */
  string make_csv( const vector<string> &data_rows )
  {
    string csv = "Probe name,Product,Latitude,Longitude,Value,Unit\r\n";
    for( const string &row : data_rows )
      csv += row + "\r\n";
    return csv;
  }//make_csv(...)
}//namespace


// ---------------------------------------------------------------------------------------------
// 1. parse_csv on the real multi-location file.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( parse_multi_location )
{
  const vector<BatchNuclide> inputs = probe_results();

  BOOST_REQUIRE_EQUAL( inputs.size(), 1400u );

  // 35 locations, each with the same 40 nuclides, no nuclide twice in a location.
  map<string,set<string>> nuclides_by_location;
  map<string,vector<size_t>> rows_by_location;
  for( size_t i = 0; i < inputs.size(); ++i )
  {
    const BatchNuclide &in = inputs[i];
    BOOST_REQUIRE_MESSAGE( !in.location.empty(), "Input " + std::to_string(i) + " has no location" );
    BOOST_REQUIRE( in.nuclide );
    nuclides_by_location[in.location].insert( in.nuclide->symbol );
    rows_by_location[in.location].push_back( i );
  }

  BOOST_CHECK_EQUAL( nuclides_by_location.size(), 35u );
  for( const pair<const string,set<string>> &nv : nuclides_by_location )
  {
    BOOST_CHECK_MESSAGE( nv.second.size() == 40u,
                        nv.first + " has " + std::to_string(nv.second.size()) + " distinct nuclides,"
                        " expected 40 (a duplicate would mean two locations were merged)" );
    BOOST_CHECK_EQUAL( rows_by_location[nv.first].size(), 40u );
  }

  // The row-index suffix is stripped, so "TeamA-01_7" groups with "TeamA-01".
  BOOST_CHECK( nuclides_by_location.count("TeamA-01") );
  BOOST_CHECK( nuclides_by_location.count("TeamG-05") );
  for( const pair<const string,set<string>> &nv : nuclides_by_location )
    BOOST_CHECK_MESSAGE( nv.first.find('_') == string::npos,
                        "Location '" + nv.first + "' still carries a row-index suffix" );

  // One lat/lon per location, echoed through for the output.
  for( const pair<const string,vector<size_t>> &nv : rows_by_location )
  {
    const vector<pair<string,string>> &first = inputs[nv.second.front()].extra_columns;
    BOOST_REQUIRE_EQUAL( first.size(), 2u );
    BOOST_CHECK_EQUAL( first[0].first, string("Latitude") );
    BOOST_CHECK_EQUAL( first[1].first, string("Longitude") );

    for( const size_t i : nv.second )
      BOOST_CHECK_MESSAGE( inputs[i].extra_columns == first,
                          nv.first + " has inconsistent Latitude/Longitude between its rows" );
  }

  // The unit is split into an activity unit and the areal label, and the trailing product text is kept.
  BOOST_CHECK_EQUAL( inputs.front().activity_unit, string("uCi") );
  BOOST_CHECK_EQUAL( inputs.front().unit_label, string("/m2") );
  BOOST_CHECK_EQUAL( inputs.front().product_suffix, string(" Deposition at 28 hrs") );

  // Activity is the value in its stated unit.
  const double expected = 3.815873E-10 * 1.0E-6 * PhysicalUnits::curie;
  BOOST_CHECK_CLOSE( inputs.front().activity, expected, 1.0E-6 );
}


// ---------------------------------------------------------------------------------------------
// 1b. The location-key rule, and the two guards that protect it.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( location_key_rule )
{
  init_data_dirs();

  // Numeric suffix is a row index.
  BOOST_CHECK_EQUAL( location_key("TeamA-01_39"), string("TeamA-01") );
  BOOST_CHECK_EQUAL( location_key("Grid_7"), string("Grid") );

  // A mixed suffix is stripped only when the head looks like a structured ID *and* the suffix carries
  //  a digit, i.e. reads as an index rather than as part of a place name.
  BOOST_CHECK_EQUAL( location_key("TeamA-01_run3"), string("TeamA-01") );
  BOOST_CHECK_EQUAL( location_key("A-1_2"), string("A-1") );

  // Anything else is left exactly as given, so distinct names are never merged.  In particular a
  //  wholly alphabetic suffix names the location, even on a structured-looking head - "Site-A_North"
  //  and "Site-A_South" are two places, not two rows of one.
  BOOST_CHECK_EQUAL( location_key("Site-A_North"), string("Site-A_North") );
  BOOST_CHECK_EQUAL( location_key("Site-A_South"), string("Site-A_South") );
  BOOST_CHECK_EQUAL( location_key("Site_North"), string("Site_North") );
  BOOST_CHECK_EQUAL( location_key("Site_South"), string("Site_South") );
  BOOST_CHECK_EQUAL( location_key("plain"), string("plain") );
  BOOST_CHECK_EQUAL( location_key("_leading"), string("_leading") );
  BOOST_CHECK_EQUAL( location_key("trailing_"), string("trailing_") );

  // Fallback 1: stripping "Site_1"/"Site_2" would merge two rows of the same nuclide into one
  //  location, which cannot be right - so the raw names must be kept.
  {
    const vector<BatchNuclide> inputs = parse_csv( make_csv( {
      "Site_1,CS-137 stuff,1.0,2.0,5.0,uCi",
      "Site_1,BA-133 stuff,1.0,2.0,5.0,uCi",
      "Site_2,CS-137 stuff,3.0,4.0,7.0,uCi",
      "Site_2,BA-133 stuff,3.0,4.0,7.0,uCi",
    } ) );

    set<string> locations;
    for( const BatchNuclide &in : inputs )
      locations.insert( in.location );

    BOOST_CHECK_EQUAL( locations.size(), 2u );
    BOOST_CHECK( locations.count("Site_1") );
    BOOST_CHECK( locations.count("Site_2") );
  }

  // Fallback 1, the other trigger: no nuclide repeats under the merge, but the passthrough cells
  //  disagree, so the rows cannot be one physical location either.
  {
    const vector<BatchNuclide> inputs = parse_csv( make_csv( {
      "Plot-1_a,CS-137 stuff,1.0,2.0,5.0,uCi",
      "Plot-1_a,BA-133 stuff,1.0,2.0,5.0,uCi",
      "Plot-1_b,CO-60 stuff,9.0,9.0,5.0,uCi",
      "Plot-1_b,NA-22 stuff,9.0,9.0,5.0,uCi",
    } ) );

    set<string> locations;
    for( const BatchNuclide &in : inputs )
      locations.insert( in.location );

    BOOST_CHECK_MESSAGE( locations.size() == 2u,
                        "Rows with differing Latitude/Longitude must not be merged into one location" );
    BOOST_CHECK( locations.count("Plot-1_a") );
    BOOST_CHECK( locations.count("Plot-1_b") );
  }

  // Fallback 2: no location ends up with more than one nuclide, so there is no grouping information
  //  in the file and it comes back ungrouped (rather than one group per row).
  {
    const vector<BatchNuclide> inputs = parse_csv( make_csv( {
      "Alpha,CS-137 stuff,1.0,2.0,5.0,uCi",
      "Bravo,BA-133 stuff,3.0,4.0,7.0,uCi",
      "Charlie,CO-60 stuff,5.0,6.0,9.0,uCi",
    } ) );

    BOOST_REQUIRE_EQUAL( inputs.size(), 3u );
    for( const BatchNuclide &in : inputs )
      BOOST_CHECK_MESSAGE( in.location.empty(),
                          "One-nuclide-per-name input should be ungrouped, got location '"
                          + in.location + "'" );
  }

  // The two fallbacks compose: when the raw names that fallback 1 restores hold one nuclide each,
  //  fallback 2 then applies and the file is ungrouped - rather than becoming one group per row.
  {
    const vector<BatchNuclide> inputs = parse_csv( make_csv( {
      "Plot-1_1,CS-137 stuff,1.0,2.0,5.0,uCi",
      "Plot-1_2,BA-133 stuff,9.0,9.0,5.0,uCi",
      "Plot-1_3,CO-60 stuff,1.0,2.0,5.0,uCi",
    } ) );

    BOOST_REQUIRE_EQUAL( inputs.size(), 3u );
    for( const BatchNuclide &in : inputs )
      BOOST_CHECK_MESSAGE( in.location.empty(),
                          "Should be ungrouped once each restored name holds a single nuclide, got '"
                          + in.location + "'" );
  }
}


// ---------------------------------------------------------------------------------------------
// 2. Grouped DecayBatchCalc::decay(): one block per location, same nuclide list in (A, Z, isomer) order.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( grouped_decay_shape )
{
  const vector<BatchNuclide> inputs = probe_results();
  const BatchDecayResult result = DecayBatchCalc::decay( inputs, options_for("48h") );

  // The output mirrors the input's columns.
  const vector<string> expected_header{ "Probe name", "Product", "Latitude", "Longitude",
                                        "Value", "Unit" };
  BOOST_CHECK( result.column_headers == expected_header );

  // 35 locations, each listing the union of all inputs' descendant closures.
  vector<string> location_order;
  map<string,size_t> rows_per_location;
  for( const vector<string> &row : result.rows )
  {
    BOOST_REQUIRE_EQUAL( row.size(), expected_header.size() );
    if( !rows_per_location.count( row[0] ) )
      location_order.push_back( row[0] );
    rows_per_location[row[0]] += 1;
  }

  BOOST_CHECK_EQUAL( location_order.size(), 35u );
  BOOST_CHECK_EQUAL( location_order.front(), string("TeamA-01") );
  BOOST_CHECK_EQUAL( location_order.back(), string("TeamG-05") );

  const size_t nuclides_per_location = rows_per_location[location_order.front()];
  BOOST_CHECK_MESSAGE( nuclides_per_location == 103u,
                      "Expected 103 nuclides per location (union of all inputs' progeny), got "
                      + std::to_string(nuclides_per_location) );

  for( const pair<const string,size_t> &nv : rows_per_location )
    BOOST_CHECK_MESSAGE( nv.second == nuclides_per_location,
                        nv.first + " lists " + std::to_string(nv.second) + " nuclides, but "
                        + location_order.front() + " lists "
                        + std::to_string(nuclides_per_location) + "; every block must list the same"
                        " nuclides so the blocks are comparable" );

  BOOST_CHECK_EQUAL( result.rows.size(), 35u * 103u );

  // The rows are location-major, and within a location in (mass number, atomic number, isomer) order.
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const vector<pair<string,double>> first_block = location_rows( result, "TeamA-01" );
  BOOST_REQUIRE_EQUAL( first_block.size(), 103u );

  for( size_t i = 1; i < first_block.size(); ++i )
  {
    const SandiaDecay::Nuclide * const prev = db->nuclide( first_block[i-1].first );
    const SandiaDecay::Nuclide * const cur = db->nuclide( first_block[i].first );
    BOOST_REQUIRE( prev && cur );
    BOOST_CHECK_MESSAGE( SandiaDecay::Nuclide::lessThanForOrdering( prev, cur ),
                        first_block[i-1].first + " should not come before " + first_block[i].first );
  }

  // Stable nuclides are listed (as an activity, hence zero), matching the legacy output.
  size_t num_stable = 0;
  for( const pair<string,double> &nv : first_block )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( nv.first );
    BOOST_REQUIRE( nuc );
    if( nuc->isStable() )
    {
      ++num_stable;
      BOOST_CHECK_MESSAGE( nv.second == 0.0,
                          "Stable " + nv.first + " should have zero activity, got "
                          + std::to_string(nv.second) );
    }
  }
  BOOST_CHECK_EQUAL( num_stable, 17u );

  // Every location is one sample set: no location's value may depend on another's.  Decaying a single
  //  location on its own must give exactly what it gives within the full file.
  vector<BatchNuclide> just_one;
  for( const BatchNuclide &in : inputs )
  {
    if( in.location == "TeamB-03" )
      just_one.push_back( in );
  }
  BOOST_REQUIRE_EQUAL( just_one.size(), 40u );

  const BatchDecayResult alone = DecayBatchCalc::decay( just_one, options_for("48h") );
  const vector<pair<string,double>> in_full = location_rows( result, "TeamB-03" );
  const vector<pair<string,double>> on_own = location_rows( alone, "TeamB-03" );

  // On its own the union of progeny is over that location's inputs only, so compare the nuclides it
  //  does report - the values must be identical, i.e. nothing leaked between locations.
  BOOST_REQUIRE( !on_own.empty() );
  for( const pair<string,double> &nv : on_own )
  {
    BOOST_CHECK_MESSAGE( fabs( value_for(in_full, nv.first) - nv.second ) <= 0.0,
                        "TeamB-03's " + nv.first + " differs when decayed alone ("
                        + std::to_string(nv.second) + ") vs. within the whole file ("
                        + std::to_string(value_for(in_full, nv.first))
                        + "); values must never mix between locations" );
  }
}


// ---------------------------------------------------------------------------------------------
// 3. The CSV format the legacy tool produced.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( grouped_csv_format )
{
  const BatchDecayResult result = DecayBatchCalc::decay( probe_results(), options_for("48h") );
  const string csv = result_to_csv( result );

  // Header, then one row per (location, nuclide).
  BOOST_CHECK( SpecUtils::starts_with( csv, "Probe name,Product,Latitude,Longitude,Value,Unit\r\n" ) );

  // CRLF on every line, including the last, and no bare LF anywhere.
  BOOST_REQUIRE( csv.size() > 2 );
  BOOST_CHECK_EQUAL( csv.substr( csv.size() - 2 ), string("\r\n") );
  for( string::size_type i = 0; i < csv.size(); ++i )
  {
    if( csv[i] == '\n' )
      BOOST_REQUIRE_MESSAGE( (i > 0) && (csv[i-1] == '\r'), "Bare LF at offset " + std::to_string(i) );
  }

  vector<string> lines;
  SpecUtils::split( lines, csv, "\r\n" );
  BOOST_REQUIRE_EQUAL( lines.size(), 1u + result.rows.size() );

  // "<Symbol> <input trailing text>+<time>", with no space before the '+'.
  BOOST_CHECK_EQUAL( result.rows.front()[1], string("Sr89 Deposition at 28 hrs+48h") );

  // Values are %.6E, and a zero prints as 0.000000E+00.
  const std::regex sci_re( "^-?[0-9]\\.[0-9]{6}E[+-][0-9]{2}$" );
  size_t num_zero = 0;
  for( const vector<string> &row : result.rows )
  {
    BOOST_REQUIRE_MESSAGE( std::regex_match( row[4], sci_re ),
                          "Value '" + row[4] + "' is not in %.6E form" );
    if( std::stod( row[4] ) == 0.0 )
    {
      ++num_zero;
      BOOST_CHECK_EQUAL( row[4], string("0.000000E+00") );
    }

    // The input's own unit is carried through; the user's Ci/Bq preference does not apply here.
    BOOST_CHECK_EQUAL( row[5], string("uCi/m2") );
  }
  BOOST_CHECK_MESSAGE( num_zero > 0u, "Expected the stable nuclides to print as exact zeros" );
}


// ---------------------------------------------------------------------------------------------
// 4. The lower-activity cut.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( min_activity_cut )
{
  const vector<BatchNuclide> inputs = probe_results();

  BatchDecayOptions opts = options_for( "48h" );
  const BatchDecayResult uncut = DecayBatchCalc::decay( inputs, opts );

  opts.min_activity = 1.0E-20;
  const BatchDecayResult cut = DecayBatchCalc::decay( inputs, opts );

  BOOST_CHECK_MESSAGE( cut.rows.size() < uncut.rows.size(), "The cut dropped nothing" );

  // Nothing below the cut survives, and nothing above it was dropped.
  set<string> kept_keys;
  for( const vector<string> &row : cut.rows )
  {
    BOOST_REQUIRE( std::stod( row[4] ) >= 1.0E-20 );
    kept_keys.insert( row[0] + "|" + row[1] );
  }

  for( const vector<string> &row : uncut.rows )
  {
    if( std::stod( row[4] ) >= 1.0E-20 )
      BOOST_REQUIRE_MESSAGE( kept_keys.count( row[0] + "|" + row[1] ),
                            "The cut dropped " + row[0] + " " + row[1] + " = " + row[4]
                            + ", which is above it" );
  }

  // A location whose every row falls below the cut disappears entirely.
  set<string> uncut_locations, cut_locations;
  for( const vector<string> &row : uncut.rows )
    uncut_locations.insert( row[0] );
  for( const vector<string> &row : cut.rows )
    cut_locations.insert( row[0] );

  BOOST_CHECK_EQUAL( uncut_locations.size(), 35u );
  BOOST_CHECK_MESSAGE( cut_locations.size() < uncut_locations.size(),
                      "Expected the all-zero locations to drop out under the cut" );

  // In the wide (ungrouped) format a row spans the time steps, so it is dropped only when every step
  //  is below the cut - which keeps the table rectangular.
  vector<BatchNuclide> ungrouped;
  for( const BatchNuclide &in : inputs )
  {
    if( in.location == "TeamA-01" )
    {
      BatchNuclide copy = in;
      copy.location.clear();
      ungrouped.push_back( copy );
    }
  }
  BOOST_REQUIRE_EQUAL( ungrouped.size(), 40u );

  BatchDecayOptions wide_opts = options_for( "48h" );
  wide_opts.num_steps = 4;
  wide_opts.min_activity = 1.0E-6;

  const BatchDecayResult wide = DecayBatchCalc::decay( ungrouped, wide_opts );
  BOOST_REQUIRE( !wide.rows.empty() );
  for( const vector<string> &row : wide.rows )
  {
    BOOST_CHECK_EQUAL( row.size(), wide.column_headers.size() );

    bool any_above = false;
    for( size_t i = 1; i < row.size(); ++i )
      any_above |= (std::stod( row[i] ) >= wide_opts.min_activity);
    BOOST_CHECK_MESSAGE( any_above, "Row '" + row[0] + "' is below the cut at every step" );
  }
}


// ---------------------------------------------------------------------------------------------
// 5. The physics: a co-decayed location reaches the expected equilibrium ratios.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( grouped_decay_equilibrium_ratios )
{
  // The reference output is a different dataset (its coordinates are disjoint from the input's), so
  //  check the physics through ratios that secular equilibrium fixes, not by diffing values.
  const BatchDecayResult result = DecayBatchCalc::decay( probe_results(), options_for("48h") );
  const vector<pair<string,double>> rows = location_rows( result, "TeamA-01" );

  // Ba137m is in equilibrium with Cs137 at its branching ratio.
  const double cs137 = value_for( rows, "Cs137" );
  const double ba137m = value_for( rows, "Ba137m" );
  BOOST_REQUIRE( cs137 > 0.0 );
  BOOST_CHECK_CLOSE( ba137m / cs137, 0.9470, 0.5 );

  // Rh106 (30 s) is in equilibrium with Ru106 (372 d) one-to-one.
  const double ru106 = value_for( rows, "Ru106" );
  const double rh106 = value_for( rows, "Rh106" );
  BOOST_REQUIRE( ru106 > 0.0 );
  BOOST_CHECK_CLOSE( rh106 / ru106, 1.0, 0.5 );

  // Rh103m from Ru103, at the 99.75% branch that feeds it.
  const double ru103 = value_for( rows, "Ru103" );
  const double rh103m = value_for( rows, "Rh103m" );
  BOOST_REQUIRE( ru103 > 0.0 );
  BOOST_CHECK_CLOSE( rh103m / ru103, 0.98926, 1.0 );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  // The input activity of one nuclide at this location, in the file's own uCi.
  const vector<BatchNuclide> inputs = probe_results();
  auto input_activity = [&inputs]( const SandiaDecay::Nuclide * const nuc ) -> double {
    for( const BatchNuclide &in : inputs )
    {
      if( (in.location == "TeamA-01") && (in.nuclide == nuc) )
        return in.activity / (1.0E-6 * PhysicalUnits::curie);
    }
    return 0.0;
  };//input_activity

  const double t = 48.0 * PhysicalUnits::hour;

  // A nuclide with no ancestor among this location's inputs is simple exponential decay.  (Cs137 is
  //  such an input; Am241 is not, since Pu241 is also an input and feeds it - see below.)
  const SandiaDecay::Nuclide * const cs137_nuc = db->nuclide( "Cs137" );
  BOOST_REQUIRE( cs137_nuc );
  const double cs137_input = input_activity( cs137_nuc );
  BOOST_REQUIRE( cs137_input > 0.0 );
  BOOST_CHECK_CLOSE( value_for( rows, "Cs137" ),
                    cs137_input * exp( -cs137_nuc->decayConstant() * t ), 0.01 );

  // Co-decaying a location means an input that feeds another input contributes to it: the location's
  //  Pu241 grows into its Am241, so Am241 must end up *above* its own bare decay.  This is exactly the
  //  mixing that must happen within a location, and never between locations.
  const SandiaDecay::Nuclide * const am241 = db->nuclide( "Am241" );
  const SandiaDecay::Nuclide * const pu241 = db->nuclide( "Pu241" );
  BOOST_REQUIRE( am241 && pu241 );

  const double am241_input = input_activity( am241 );
  const double pu241_input = input_activity( pu241 );
  BOOST_REQUIRE( (am241_input > 0.0) && (pu241_input > 0.0) );

  const double am241_bare = am241_input * exp( -am241->decayConstant() * t );
  const double am241_got = value_for( rows, "Am241" );
  BOOST_CHECK_MESSAGE( am241_got > am241_bare,
                      "Am241 (" + std::to_string(am241_got) + ") should exceed its own bare decay ("
                      + std::to_string(am241_bare) + ") thanks to in-growth from the location's Pu241" );

  // And it must equal bare decay plus exactly what the Pu241 grows in.
  SandiaDecay::NuclideMixture pu_only;
  pu_only.addNuclideByActivity( pu241, pu241_input * 1.0E-6 * PhysicalUnits::curie );
  const double from_pu241 = pu_only.activity( t, am241 ) / (1.0E-6 * PhysicalUnits::curie);
  BOOST_CHECK_CLOSE( am241_got, am241_bare + from_pu241, 0.01 );
}


// ---------------------------------------------------------------------------------------------
// 6. Negative decay times.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( negative_time_single_nuclide )
{
  init_data_dirs();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const co60 = db->nuclide( "Co60" );
  BOOST_REQUIRE( co60 );

  // A nuclide with no ancestor in the set must reduce to A_past = A_now / exp(-lambda*|t|).
  BatchNuclide in;
  in.nuclide = co60;
  in.nuclide_str = "Co60";
  in.activity = 1.0 * PhysicalUnits::curie;

  BatchDecayOptions opts = options_for( "-1y" );
  BOOST_REQUIRE( opts.time_span < 0.0 );

  const BatchDecayResult result = DecayBatchCalc::decay( vector<BatchNuclide>{in}, opts );
  BOOST_REQUIRE_EQUAL( result.rows.size(), 1u );
  BOOST_REQUIRE_EQUAL( result.rows.front().size(), 2u );

  // Step 0 is the measurement itself; a single-step run evaluates at time_span, i.e. the past state.
  const double past = std::stod( result.rows.front()[1] );
  const double expected = 1.0 / exp( -co60->decayConstant() * 1.0 * PhysicalUnits::year );
  BOOST_CHECK_CLOSE( past, expected, 0.01 );
  BOOST_CHECK_MESSAGE( past > 1.0, "Looking back, a decaying nuclide must have had more activity" );

  // Zero is the one time span that has no meaning.
  BatchDecayOptions zero_opts = opts;
  zero_opts.time_span = 0.0;
  BOOST_CHECK_THROW( DecayBatchCalc::decay( vector<BatchNuclide>{in}, zero_opts ), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( negative_time_coupled_nuclides )
{
  init_data_dirs();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const cs137 = db->nuclide( "Cs137" );
  const SandiaDecay::Nuclide * const ba137m = db->nuclide( "Ba137m" );
  BOOST_REQUIRE( cs137 && ba137m );

  // Cs137 and its Ba137m progeny are coupled: part of today's Ba137m grew in from the Cs137, and over
  //  48 h (1129 Ba137m half-lives) none of the original Ba137m is left - so its own past activity is
  //  not recoverable at all, rather than being A_now/exp(-lambda*t) (which would divide by zero).
  vector<BatchNuclide> inputs( 2 );
  inputs[0].nuclide = cs137;
  inputs[0].nuclide_str = "Cs137";
  inputs[0].activity = 1.0 * PhysicalUnits::curie;
  inputs[1].nuclide = ba137m;
  inputs[1].nuclide_str = "Ba137m";
  inputs[1].activity = 0.9470 * PhysicalUnits::curie;   // in equilibrium, so the input is consistent

  BatchDecayOptions opts = options_for( "-48h" );
  opts.mix_input = true;       // one mixture, so they are solved together

  BatchDecayResult result = DecayBatchCalc::decay( inputs, opts );

  // Neither value may be infinite or NaN.
  for( const vector<string> &row : result.rows )
  {
    for( size_t i = 1; i < row.size(); ++i )
    {
      const double v = std::stod( row[i] );
      BOOST_REQUIRE_MESSAGE( !std::isnan(v) && !std::isinf(v),
                            "Got " + row[i] + " for '" + row[0] + "'" );
    }
  }

  // Ba137m's unrecoverable past activity is called out, not silently reported as a number.
  BOOST_CHECK_MESSAGE( result.warnings.find("Ba137m") != string::npos,
                      "Expected a warning about Ba137m's past activity; warnings were: "
                      + result.warnings );

  // A consistent input must not also be reported as inconsistent.
  BOOST_CHECK_MESSAGE( result.warnings.find("not self-consistent") == string::npos,
                      "An equilibrium input should need no inconsistency warning; got: "
                      + result.warnings );

  // The recovered past state, decayed forward again, must reproduce the measured Cs137.
  const double past_cs137 = std::stod( result.rows.front()[1] );
  BOOST_CHECK_EQUAL( SpecUtils::starts_with( result.rows.front()[0], "Cs137" ), true );

  SandiaDecay::NuclideMixture check;
  check.addNuclideByActivity( cs137, past_cs137 * PhysicalUnits::curie );
  const double forward = check.activity( 48.0 * PhysicalUnits::hour, cs137 ) / PhysicalUnits::curie;
  BOOST_CHECK_CLOSE( forward, 1.0, 0.01 );

  // An input that equilibrium forbids cannot have come from any past state, and must say so.
  inputs[1].activity = 0.890 * PhysicalUnits::curie;   // the real file's ratio, ~6% off equilibrium
  result = DecayBatchCalc::decay( inputs, opts );
  BOOST_CHECK_MESSAGE( result.warnings.find("not self-consistent") != string::npos,
                      "An input inconsistent with equilibrium should be flagged; got: "
                      + result.warnings );
}


BOOST_AUTO_TEST_CASE( negative_time_too_far_back )
{
  // 48 h is ~167 half-lives of the file's Pr144, and no parent of it is among the 40 inputs, so its
  //  past activity is not recoverable at all - recovering it would mean scaling the measurement up by
  //  1E+50.  That has to be an error the user can act on, not a number.
  const vector<BatchNuclide> inputs = probe_results();

  try
  {
    DecayBatchCalc::decay( inputs, options_for("-48h") );
    BOOST_FAIL( "Looking back 48 h should have been refused (Pr144 cannot be recovered)" );
  }catch( std::exception &e )
  {
    const string msg = e.what();
    BOOST_CHECK_MESSAGE( msg.find("Pr144") != string::npos,
                        "The error should name the nuclide at fault; got: " + msg );
    BOOST_CHECK_MESSAGE( msg.find("shorter time") != string::npos,
                        "The error should tell the user what to do; got: " + msg );
  }

  // A nuclide *with* a parent among the inputs is not an error, however far back: Ba137m is over a
  //  thousand half-lives in 48 h, but its Cs137 parent accounts for the measurement.
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  vector<BatchNuclide> pair_inputs( 2 );
  pair_inputs[0].nuclide = db->nuclide( "Cs137" );
  pair_inputs[0].nuclide_str = "Cs137";
  pair_inputs[0].activity = 1.0 * PhysicalUnits::curie;
  pair_inputs[1].nuclide = db->nuclide( "Ba137m" );
  pair_inputs[1].nuclide_str = "Ba137m";
  pair_inputs[1].activity = 0.9470 * PhysicalUnits::curie;
  BOOST_REQUIRE( pair_inputs[0].nuclide && pair_inputs[1].nuclide );

  BatchDecayOptions pair_opts = options_for( "-48h" );
  pair_opts.mix_input = true;
  BOOST_CHECK_NO_THROW( DecayBatchCalc::decay( pair_inputs, pair_opts ) );
}


BOOST_AUTO_TEST_CASE( negative_time_grouped )
{
  // The real multi-location file back-decayed: every location is solved on its own, and the same
  //  warning must not be repeated once per location.  12 h is within reach of every input that has no
  //  parent among the others (the tightest is Pr144, at 13.5 h - see negative_time_too_far_back).
  BatchDecayOptions opts = options_for( "-12h" );
  const BatchDecayResult result = DecayBatchCalc::decay( probe_results(), opts );

  BOOST_CHECK_EQUAL( result.rows.size(), 35u * 103u );

  for( const vector<string> &row : result.rows )
  {
    const double v = std::stod( row[4] );
    BOOST_REQUIRE_MESSAGE( !std::isnan(v) && !std::isinf(v),
                          "Got " + row[4] + " for " + row[0] + " " + row[1] );
    BOOST_REQUIRE_MESSAGE( v >= 0.0, "Negative activity " + row[4] + " for " + row[0] + " " + row[1] );
  }

  // A negative time supplies its own sign, rather than reading "...+-12h".
  BOOST_CHECK_EQUAL( result.rows.front()[1], string("Sr89 Deposition at 28 hrs-12h") );

  // This file's parent/progeny ratios are not those of a single 12-hour-old deposition, so the solve
  //  must say the inputs are not self-consistent rather than silently returning a state that does not
  //  decay forward to them.
  BOOST_REQUIRE_MESSAGE( result.warnings.find("not self-consistent") != string::npos,
                        "Expected an inconsistency warning for the real file; got: "
                        + result.warnings );

  // The warnings quote activities in the table's own unit, not becquerel.
  BOOST_CHECK_MESSAGE( result.warnings.find(" Bq") == string::npos,
                      "Warnings should quote the input's unit, not Bq; got: " + result.warnings );
  BOOST_CHECK_MESSAGE( result.warnings.find("uCi/m2") != string::npos,
                      "Expected the input's unit in the warnings; got: " + result.warnings );

  // One message per affected nuclide, however many of the 35 locations hit it.
  size_t num_unrecoverable_msgs = 0;
  string::size_type pos = result.warnings.find( "cannot be determined" );
  while( pos != string::npos )
  {
    ++num_unrecoverable_msgs;
    pos = result.warnings.find( "cannot be determined", pos + 1 );
  }
  BOOST_CHECK_MESSAGE( num_unrecoverable_msgs <= 1u,
                      "The unrecoverable-nuclide warning should be collapsed to one message, saw "
                      + std::to_string(num_unrecoverable_msgs) );
}


BOOST_AUTO_TEST_CASE( negative_time_duplicate_nuclides )
{
  init_data_dirs();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const cs137 = db->nuclide( "Cs137" );
  BOOST_REQUIRE( cs137 );

  // Two rows of one nuclide add, in a mixture; neither may be treated as having produced the other.
  //  Both orders, since the smaller-second case used to produce a spurious inconsistency warning.
  const double survive = exp( -cs137->decayConstant() * PhysicalUnits::year );
  for( const pair<double,double> &acts : vector<pair<double,double>>{ {1.0, 2.0}, {2.0, 1.0} } )
  {
    vector<BatchNuclide> inputs( 2 );
    for( size_t i = 0; i < 2; ++i )
    {
      inputs[i].nuclide = cs137;
      inputs[i].nuclide_str = "Cs137";
    }
    inputs[0].activity = acts.first * PhysicalUnits::curie;
    inputs[1].activity = acts.second * PhysicalUnits::curie;

    BatchDecayOptions opts = options_for( "-1y" );
    opts.mix_input = true;
    opts.num_steps = 2;

    const BatchDecayResult result = DecayBatchCalc::decay( inputs, opts );

    BOOST_CHECK_CLOSE( wide_value( result, "Cs137", 0 ), 3.0, 0.01 );            // now
    BOOST_CHECK_CLOSE( wide_value( result, "Cs137", 1 ), 3.0 / survive, 0.01 );  // a year ago
    BOOST_CHECK_MESSAGE( result.warnings.find("not self-consistent") == string::npos,
                        "Repeated rows of one nuclide are consistent; got: " + result.warnings );
  }

  // With a parent among the inputs, its in-growth into a repeated nuclide must be taken off once in
  //  total, not once per row.  A day of 1 Ci Sr90 grows in ~0.23 Ci of Y90, so 0.2 + 0.2 Ci is also a
  //  case where each row alone is below the in-growth, though their sum is not.
  const SandiaDecay::Nuclide * const sr90 = db->nuclide( "Sr90" );
  const SandiaDecay::Nuclide * const y90 = db->nuclide( "Y90" );
  BOOST_REQUIRE( sr90 && y90 );
  for( const pair<double,double> &acts : vector<pair<double,double>>{ {0.5, 0.6}, {0.6, 0.5}, {0.2, 0.2} } )
  {
    vector<BatchNuclide> inputs( 3 );
    inputs[0].nuclide = sr90;
    inputs[0].nuclide_str = "Sr90";
    inputs[0].activity = 1.0 * PhysicalUnits::curie;
    inputs[1].nuclide = inputs[2].nuclide = y90;
    inputs[1].nuclide_str = inputs[2].nuclide_str = "Y90";
    inputs[1].activity = acts.first * PhysicalUnits::curie;
    inputs[2].activity = acts.second * PhysicalUnits::curie;

    BatchDecayOptions opts = options_for( "-1d" );
    opts.mix_input = true;
    opts.num_steps = 2;

    const BatchDecayResult result = DecayBatchCalc::decay( inputs, opts );
    BOOST_CHECK_CLOSE( wide_value( result, "Sr90", 0 ), 1.0, 0.01 );
    BOOST_CHECK_CLOSE( wide_value( result, "Y90", 0 ), acts.first + acts.second, 0.01 );
    BOOST_CHECK_MESSAGE( result.warnings.empty(), "Expected no warnings; got: " + result.warnings );
  }
}


BOOST_AUTO_TEST_CASE( long_initial_ages )
{
  init_data_dirs();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const ba140 = db->nuclide( "Ba140" );
  const SandiaDecay::Nuclide * const mo99 = db->nuclide( "Mo99" );
  BOOST_REQUIRE( ba140 && mo99 );

  // 2 y is 57 half-lives of Ba140.  Looking back used to trip SandiaDecay's precision check in
  //  addAgedNuclideByNumAtoms (at ~45 half-lives), though decaying the same input forwards works.
  BatchNuclide aged;
  aged.nuclide = ba140;
  aged.nuclide_str = "Ba140";
  aged.activity = 1.0 * PhysicalUnits::curie;
  aged.age = 2.0 * PhysicalUnits::year;

  const double survive = exp( -ba140->decayConstant() * PhysicalUnits::day );
  for( const bool mix : { false, true } )
  {
    BatchDecayOptions opts = options_for( "-1d" );
    opts.mix_input = mix;
    opts.num_steps = 2;

    BatchDecayResult result;
    BOOST_REQUIRE_NO_THROW( result = DecayBatchCalc::decay( vector<BatchNuclide>{ aged }, opts ) );
    BOOST_CHECK_CLOSE( wide_value( result, "Ba140", 0 ), 1.0, 0.01 );
    BOOST_CHECK_CLOSE( wide_value( result, "Ba140", 1 ), 1.0 / survive, 0.01 );
  }

  // An age so long that none of the nuclide could have survived it (10 y is ~1300 half-lives of Mo99)
  //  cannot be computed at all, and must be refused with a message saying why, not give NaN or inf.
  BatchNuclide ancient = aged;
  ancient.nuclide = mo99;
  ancient.nuclide_str = "Mo99";
  ancient.age = 10.0 * PhysicalUnits::year;

  for( const char * const time : { "1d", "-1d" } )
  {
    try
    {
      DecayBatchCalc::decay( vector<BatchNuclide>{ ancient }, options_for( time ) );
      BOOST_ERROR( string("An initial age of 1300 half-lives should be refused, decaying ") + time );
    }catch( std::exception &e )
    {
      const string msg = e.what();
      BOOST_CHECK_MESSAGE( (msg.find("Mo99") != string::npos) && (msg.find("half-lives") != string::npos),
                          "The error should name the nuclide and say why; got: " + msg );
    }
  }
}


BOOST_AUTO_TEST_CASE( negative_time_aged_inputs )
{
  init_data_dirs();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const pu241 = db->nuclide( "Pu241" );
  const SandiaDecay::Nuclide * const am241 = db->nuclide( "Am241" );
  BOOST_REQUIRE( pu241 && am241 );

  const double ci = PhysicalUnits::curie;
  const double year = PhysicalUnits::year;

  // A 10-year-old Pu241 sample measured along with its Am241, which is what 10 years of in-growth
  //  gives plus a little of its own.  Looking back 2 years the Pu241 was 8 years old; the solve must
  //  account for the Am241 that aged sample already held, rather than count it twice.
  SandiaDecay::NuclideMixture aged10;
  aged10.addAgedNuclideByActivity( pu241, 1.0 * ci, 10.0 * year );
  const double am_ingrown = aged10.activity( 0.0, am241 ) / ci;
  const double am_measured = am_ingrown + 0.01;

  vector<BatchNuclide> inputs( 2 );
  inputs[0].nuclide = pu241;
  inputs[0].nuclide_str = "Pu241";
  inputs[0].activity = 1.0 * ci;
  inputs[0].age = 10.0 * year;
  inputs[1].nuclide = am241;
  inputs[1].nuclide_str = "Am241";
  inputs[1].activity = am_measured * ci;

  BatchDecayOptions opts = options_for( "-2y" );
  opts.mix_input = true;
  opts.num_steps = 2;

  const BatchDecayResult result = DecayBatchCalc::decay( inputs, opts );

  // Step 0 is the measurement, which the recovered past state must decay forward to.
  BOOST_CHECK_CLOSE( wide_value( result, "Pu241", 0 ), 1.0, 0.01 );
  BOOST_CHECK_CLOSE( wide_value( result, "Am241", 0 ), am_measured, 0.01 );
  BOOST_CHECK_MESSAGE( result.warnings.empty(), "Expected no warnings; got: " + result.warnings );

  // The past state, built independently: an 8-year-old Pu241 sample plus whatever Am241 of its own
  //  decays to the part of the measurement the Pu241 does not account for.
  const double pu_past = 1.0 / exp( -pu241->decayConstant() * 2.0 * year );
  SandiaDecay::NuclideMixture aged8;
  aged8.addAgedNuclideByActivity( pu241, pu_past * ci, 8.0 * year );
  const double am_from_pu = aged8.activity( 2.0 * year, am241 ) / ci;
  const double am_own_past = (am_measured - am_from_pu) / exp( -am241->decayConstant() * 2.0 * year );
  const double am_past = aged8.activity( 0.0, am241 ) / ci + am_own_past;

  BOOST_CHECK_CLOSE( wide_value( result, "Pu241", 1 ), pu_past, 0.01 );
  BOOST_CHECK_CLOSE( wide_value( result, "Am241", 1 ), am_past, 0.01 );

  // A stated age shorter than the look-back: the sample did not exist then, so it is taken as fresh
  //  at the past time (no Am241 yet), and a note says so.
  vector<BatchNuclide> young( 1, inputs[0] );
  young[0].age = 1.0 * year;
  BatchDecayOptions young_opts = options_for( "-2y" );
  young_opts.show_progeny = true;
  young_opts.num_steps = 2;

  const BatchDecayResult young_result = DecayBatchCalc::decay( young, young_opts );
  BOOST_CHECK_MESSAGE( young_result.warnings.find( "initial age of Pu241" ) != string::npos,
                      "Expected a note about Pu241's short age; got: " + young_result.warnings );
  BOOST_CHECK_EQUAL( wide_value( young_result, "Am241", 1 ), 0.0 );
  BOOST_CHECK( wide_value( young_result, "Am241", 0 ) > 0.0 );
}


namespace
{
  /** Every nuclide's value in one step of a wide result, keyed by the label's leading symbol. */
  map<string,double> wide_step( const BatchDecayResult &result, const size_t step )
  {
    map<string,double> answer;
    for( const vector<string> &row : result.rows )
    {
      BOOST_REQUIRE( (step + 1) < row.size() );
      const string &label = row.front();
      answer[label.substr( 0, label.find(' ') )] += std::stod( row[step + 1] );
    }
    return answer;
  }//wide_step(...)


  /** Seeds a fresh mixture with every nuclide of a reported past state (activities in any one unit),
   decays it forward by `dt`, and returns the activity of each of `symbols` then, in that same unit.
   The past state holds every descendant, so no ages are needed to reproduce the present. */
  map<string,double> decay_forward( const map<string,double> &past, const double dt,
                                    const vector<string> &symbols )
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE( db );

    SandiaDecay::NuclideMixture mix;
    for( const pair<const string,double> &nv : past )
    {
      const SandiaDecay::Nuclide * const nuc = db->nuclide( nv.first );
      BOOST_REQUIRE_MESSAGE( nuc, "Unknown nuclide '" + nv.first + "' in the result" );
      if( (nv.second > 0.0) && !nuc->isStable() )
        mix.addNuclideByActivity( nuc, nv.second );
    }

    map<string,double> answer;
    for( const string &sym : symbols )
      answer[sym] = mix.activity( dt, db->nuclide( sym ) );
    return answer;
  }//decay_forward(...)
}//namespace


// The reported past state must decay forward to the measurement, when an input is part progeny of
//  another input and part its own (so the parent alone does not explain it).  Checked independently
//  of the tool: its whole past column seeds a fresh NuclideMixture, which is decayed to T=0.
BOOST_AUTO_TEST_CASE( negative_time_round_trip )
{
  init_data_dirs();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const double ci = PhysicalUnits::curie;
  const double year = PhysicalUnits::year;

  struct Input { string symbol; double activity_ci; double age; };
  struct Case { string name; vector<Input> inputs; string time; };

  // Ba140 aged 2 y is 57 half-lives old; La140 is a little above its transient equilibrium (1.15).
  const Case long_age = { "Ba140 aged 2 y + La140", { {"Ba140", 1.0, 2.0 * year}, {"La140", 1.2, 0.0} }, "-1d" };

  // A repeated nuclide, one row aged (so with Po210 grown in) and one fresh, below its parent.
  const Case repeated = { "Ra226 + aged and fresh Pb210 + Po210",
                          { {"Ra226", 1.0, 0.0}, {"Pb210", 0.2, 5.0 * year}, {"Pb210", 0.3, 0.0},
                            {"Po210", 0.8, 0.0} }, "-1y" };

  // Measured progeny above what the parent's in-growth gives, so each has a past amount of its own.
  //  Aged Ra226 already holds ~0.8 Ci of Pb210 and Po210, which the look-back must not count twice.
  SandiaDecay::NuclideMixture ra_aged;
  ra_aged.addAgedNuclideByActivity( db->nuclide("Ra226"), 1.0 * ci, 50.0 * year );
  const double pb_in_aged = ra_aged.activity( 0.0, db->nuclide("Pb210") ) / ci;
  const double po_in_aged = ra_aged.activity( 0.0, db->nuclide("Po210") ) / ci;

  const vector<Case> cases = {
    { "Ra226 chain, fresh", { {"Ra226", 1.0, 0.0}, {"Pb210", 0.3, 0.0}, {"Po210", 0.5, 0.0} }, "-1y" },
    { "Ra226 chain, aged 50 y", { {"Ra226", 1.0, 50.0 * year}, {"Pb210", pb_in_aged + 0.2, 0.0},
                                  {"Po210", po_in_aged + 0.3, 0.0} }, "-1y" },
    { "Pu241 aged 10 y + Am241", { {"Pu241", 1.0, 10.0 * year}, {"Am241", 0.05, 0.0} }, "-5y" },
    { "Sr90 + excess Y90", { {"Sr90", 1.0, 0.0}, {"Y90", 1.5, 0.0} }, "-2d" },
    { "U238 + U234 + Th230, aged", { {"U238", 1.0, 1000.0 * year}, {"U234", 1.2, 0.0},
                                     {"Th230", 0.1, 0.0} }, "-100y" },
    long_age,
    repeated
  };

  for( const Case &c : cases )
  {
    BOOST_TEST_MESSAGE( "Round trip: " << c.name );

    vector<BatchNuclide> inputs;
    vector<string> symbols;
    map<string,double> measured;   // rows of one nuclide add
    for( const Input &in : c.inputs )
    {
      BatchNuclide bn;
      bn.nuclide = db->nuclide( in.symbol );
      BOOST_REQUIRE( bn.nuclide );
      bn.nuclide_str = in.symbol;
      bn.activity = in.activity_ci * ci;
      bn.age = in.age;
      inputs.push_back( bn );
      symbols.push_back( in.symbol );
      measured[in.symbol] += in.activity_ci;
    }

    BatchDecayOptions opts = options_for( c.time );
    opts.mix_input = true;
    opts.num_steps = 2;

    const BatchDecayResult result = DecayBatchCalc::decay( inputs, opts );
    BOOST_CHECK_MESSAGE( result.warnings.empty(), c.name + ": expected no warnings; got: " + result.warnings );

    const map<string,double> now = decay_forward( wide_step( result, 1 ), -opts.time_span, symbols );
    for( const pair<const string,double> &nv : measured )
    {
      BOOST_CHECK_MESSAGE( fabs( now.at(nv.first) - nv.second ) <= 1.0E-5 * nv.second,
                          c.name + ": " + nv.first + " decays forward to "
                          + std::to_string( now.at(nv.first) ) + " Ci, but was measured at "
                          + std::to_string( nv.second ) + " Ci" );
    }
  }//for( each case )

  // Progeny *under* what its parent explains: no non-negative past amount of it fits, so its own share
  //  is zero and a warning says the inputs are inconsistent.  The parent must still round-trip, and
  //  the progeny comes back at the parent's in-growth - above the measurement, never below.
  {
    vector<BatchNuclide> inputs( 2 );
    inputs[0].nuclide = db->nuclide( "Ra226" );
    inputs[0].nuclide_str = "Ra226";
    inputs[0].activity = 1.0 * ci;
    inputs[0].age = 50.0 * year;
    inputs[1].nuclide = db->nuclide( "Pb210" );
    inputs[1].nuclide_str = "Pb210";
    inputs[1].activity = 0.5 * pb_in_aged * ci;

    BatchDecayOptions opts = options_for( "-1y" );
    opts.mix_input = true;
    opts.num_steps = 2;

    const BatchDecayResult result = DecayBatchCalc::decay( inputs, opts );
    BOOST_CHECK_MESSAGE( (result.warnings.find("not self-consistent") != string::npos)
                         && (result.warnings.find("Pb210") != string::npos),
                        "Expected an inconsistency warning naming Pb210; got: " + result.warnings );

    const map<string,double> now = decay_forward( wide_step( result, 1 ), 1.0 * year, {"Ra226", "Pb210"} );
    BOOST_CHECK_CLOSE( now.at("Ra226"), 1.0, 1.0E-3 );
    BOOST_CHECK_CLOSE( now.at("Pb210"), pb_in_aged, 0.1 );
  }

  // Grouped (multi-location) input goes through the same solve, one location at a time.
  {
    const vector<BatchNuclide> inputs = parse_csv(
      "Probe name,Product,Value,Unit\n"
      "SiteA_1,Ra226 soil,1,uCi\n"
      "SiteA_2,Pb210 soil,0.3,uCi\n"
      "SiteA_3,Po210 soil,0.5,uCi\n"
      "SiteB_1,Sr90 soil,1,uCi\n"
      "SiteB_2,Y90 soil,1.5,uCi\n" );
    BOOST_REQUIRE_EQUAL( inputs.size(), 5u );

    BatchDecayOptions opts = options_for( "-1d" );
    const BatchDecayResult result = DecayBatchCalc::decay( inputs, opts );
    BOOST_CHECK_MESSAGE( result.warnings.empty(), "Grouped: expected no warnings; got: " + result.warnings );

    for( const BatchNuclide &in : inputs )
    {
      map<string,double> past;
      for( const pair<string,double> &nv : location_rows( result, in.location ) )
        past[nv.first] += nv.second;

      const string sym = in.nuclide->symbol;
      const double now = decay_forward( past, 1.0 * PhysicalUnits::day, {sym} ).at( sym );
      BOOST_CHECK_CLOSE( now, in.activity / (1.0E-6 * ci), 1.0E-3 );
    }
  }
}


BOOST_AUTO_TEST_CASE( negative_time_unmixed_hint )
{
  init_data_dirs();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  vector<BatchNuclide> inputs( 2 );
  inputs[0].nuclide = db->nuclide( "Cs137" );
  inputs[0].nuclide_str = "Cs137";
  inputs[0].activity = 1.0 * PhysicalUnits::curie;
  inputs[1].nuclide = db->nuclide( "Ba137m" );
  inputs[1].nuclide_str = "Ba137m";
  inputs[1].activity = 0.9470 * PhysicalUnits::curie;
  BOOST_REQUIRE( inputs[0].nuclide && inputs[1].nuclide );

  // Unmixed, the Ba137m is solved on its own, so its parent being an input cannot help - but the
  //  message must say that mixing would, rather than claim there is no parent among the inputs.
  try
  {
    DecayBatchCalc::decay( inputs, options_for("-48h") );
    BOOST_FAIL( "Looking back 48 h at unmixed Ba137m should have been refused" );
  }catch( std::exception &e )
  {
    const string msg = e.what();
    BOOST_CHECK_MESSAGE( (msg.find("Mix inputs") != string::npos) && (msg.find("Cs137") != string::npos),
                        "The error should point at Mix inputs and the parent; got: " + msg );
  }
}


// ---------------------------------------------------------------------------------------------
// A byte-order mark (spreadsheet "CSV UTF-8" exports) must not break either format.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( parse_csv_byte_order_mark )
{
  init_data_dirs();
  const string bom = "\xEF\xBB\xBF";

  const vector<BatchNuclide> simple = parse_csv( bom + "I131, 5.9 uCi\r\nCs137, 1 mCi\r\n" );
  BOOST_REQUIRE_EQUAL( simple.size(), 2u );
  BOOST_REQUIRE( simple[0].nuclide );
  BOOST_CHECK_EQUAL( simple[0].nuclide->symbol, string("I131") );

  const vector<BatchNuclide> keyed = parse_csv( bom + "Product,Value,Unit\r\nCs137 stuff,5,uCi\r\n" );
  BOOST_REQUIRE_EQUAL( keyed.size(), 1u );
  BOOST_REQUIRE( keyed[0].nuclide );
  BOOST_CHECK_EQUAL( keyed[0].nuclide->symbol, string("Cs137") );
  BOOST_CHECK_CLOSE( keyed[0].activity, 5.0E-6 * PhysicalUnits::curie, 1.0E-6 );
}


// A value that is not a measurable activity must be refused, in every format, rather than silently
//  decayed as zero (and a file holding one is then not claimed when dropped).
BOOST_AUTO_TEST_CASE( parse_csv_rejects_bad_activities )
{
  init_data_dirs();

  for( const string bad : { "inf", "-inf", "nan", "-5", "1e999", "5abc" } )
  {
    const vector<string> texts = {
      "Cs137, " + bad + "\n",
      "Product,Value,Unit\nCs137," + bad + ",uCi\n",
      "Product,Value\nCs137," + bad + "\n",
      "Probe name,Product,Value,Unit\nA_1,Cs137," + bad + ",uCi\nA_2,Co60,1,uCi\n",
      "Probe name,Product,Value\nA_1,Cs137," + bad + "\nA_2,Co60,1\n"
    };

    for( const string &text : texts )
    {
      BOOST_CHECK_MESSAGE( !is_candidate_file( text, true ), "Claimed: " + text );
      BOOST_CHECK_THROW( parse_csv( text ), std::runtime_error );
    }
  }//for( each bad value )

  for( const string bad : { "-5 uCi", "nan uCi", "inf Bq" } )
    BOOST_CHECK_THROW( parse_csv( "Cs137, " + bad ), std::runtime_error );

  // With its own Unit column, a keyed format's Value is just the number.
  BOOST_CHECK_THROW( parse_csv( "Product,Value\nCs137,5 uCi\n" ), std::runtime_error );

  // Zero is a measurement ("none detected").
  BOOST_CHECK_EQUAL( parse_csv( "Cs137, 0" ).at(0).activity, 0.0 );
  BOOST_CHECK_EQUAL( parse_csv( "Product,Value,Unit\nCs137,0,uCi\n" ).at(0).activity, 0.0 );

  // Nor may other callers slip one past decay().
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  for( const double bad : { -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() } )
  {
    BatchNuclide in;
    in.nuclide = db->nuclide( "Cs137" );
    in.nuclide_str = "Cs137";
    in.activity = bad * PhysicalUnits::curie;
    BOOST_REQUIRE( in.nuclide );

    for( const bool mix : { false, true } )
    {
      BatchDecayOptions opts = options_for( "1y" );
      opts.mix_input = mix;
      BOOST_CHECK_THROW( DecayBatchCalc::decay( vector<BatchNuclide>{ in }, opts ), std::runtime_error );
    }
  }
}


// Spreadsheets quote a cell that holds the delimiter, and R's write.csv quotes every text cell.
BOOST_AUTO_TEST_CASE( parse_csv_quoted_cells )
{
  init_data_dirs();

  // A quoted Product holding a comma, with no Unit column: split at that comma, Latitude would
  //  silently be read as the Value.
  const vector<BatchNuclide> excel = parse_csv(
    "Probe name,Product,Latitude,Longitude,Value\n"
    "TeamA-01_1,\"Cs137 fallout, day 2\",42.5,-88.1,5\n"
    "TeamA-01_2,\"I131 fallout, day 2\",42.5,-88.1,7\n" );
  BOOST_REQUIRE_EQUAL( excel.size(), 2u );
  BOOST_CHECK_EQUAL( excel[0].location, string("TeamA-01") );
  BOOST_CHECK_EQUAL( excel[1].location, string("TeamA-01") );
  BOOST_CHECK_EQUAL( excel[0].product_suffix, string(" fallout, day 2") );
  BOOST_CHECK_CLOSE( excel[0].activity, 5.0 * PhysicalUnits::becquerel, 1.0E-9 );
  BOOST_CHECK_CLOSE( excel[1].activity, 7.0 * PhysicalUnits::becquerel, 1.0E-9 );
  BOOST_REQUIRE_EQUAL( excel[0].extra_columns.size(), 2u );
  BOOST_CHECK_EQUAL( excel[0].extra_columns[0].second, string("42.5") );
  BOOST_CHECK_EQUAL( excel[0].extra_columns[1].second, string("-88.1") );

  // What the tool itself writes (it quotes that Product) must read back the same.
  const BatchDecayResult result = DecayBatchCalc::decay( excel, options_for( "1d" ) );
  const vector<BatchNuclide> reread = parse_csv( result_to_csv( result ) );
  BOOST_REQUIRE_EQUAL( reread.size(), result.rows.size() );
  BOOST_CHECK_EQUAL( reread[0].product_suffix, string(" fallout, day 2+1d") );
  BOOST_CHECK( reread[0].extra_columns == excel[0].extra_columns );

  // Every cell quoted, as R writes them; also recognized when dropped.
  const string r_style = "\"Probe name\",\"Product\",\"Value\",\"Unit\"\n"
                         "\"A_1\",\"Cs137 soil\",5,\"uCi\"\n"
                         "\"A_2\",\"Co60 soil\",2.5,\"uCi\"\n";
  const vector<BatchNuclide> r = parse_csv( r_style );
  BOOST_REQUIRE_EQUAL( r.size(), 2u );
  BOOST_CHECK_EQUAL( r[0].location, string("A") );
  BOOST_CHECK_CLOSE( r[1].activity, 2.5E-6 * PhysicalUnits::curie, 1.0E-9 );
  BOOST_CHECK( is_candidate_file( r_style, true ) );

  // The simple format quoted, and quotes in a tab-separated file.
  const vector<BatchNuclide> simple = parse_csv( "\"Cs137\",\"5 uCi\"\n\"Co60\",3\n" );
  BOOST_REQUIRE_EQUAL( simple.size(), 2u );
  BOOST_CHECK_CLOSE( simple[0].activity, 5.0E-6 * PhysicalUnits::curie, 1.0E-9 );

  const vector<BatchNuclide> tsv = parse_csv( "Product\tValue\tUnit\n\"Cs137 from site A, north\"\t5\tuCi\n" );
  BOOST_REQUIRE_EQUAL( tsv.size(), 1u );
  BOOST_CHECK_CLOSE( tsv[0].activity, 5.0E-6 * PhysicalUnits::curie, 1.0E-9 );

  // A backslash is just text (by default the tokenizer would take it as an escape, and throw).
  BOOST_CHECK_EQUAL( parse_csv( "Product,Value,Unit,Notes\nCs137,5,uCi,C:\\data\\run1.csv\n" ).size(), 1u );
}


// Tab-separated input (e.g. pasted from a spreadsheet) must parse the same as its CSV equivalent.
BOOST_AUTO_TEST_CASE( parse_tab_separated )
{
  init_data_dirs();

  const auto same = []( const vector<BatchNuclide> &a, const vector<BatchNuclide> &b ){
    BOOST_REQUIRE_EQUAL( a.size(), b.size() );
    for( size_t i = 0; i < a.size(); ++i )
    {
      BOOST_CHECK( a[i].nuclide == b[i].nuclide );
      BOOST_CHECK_EQUAL( a[i].activity, b[i].activity );
      BOOST_CHECK_EQUAL( a[i].unit_label, b[i].unit_label );
      BOOST_CHECK_EQUAL( a[i].location, b[i].location );
      BOOST_CHECK_EQUAL( a[i].product_suffix, b[i].product_suffix );
      BOOST_CHECK_EQUAL( a[i].activity_unit, b[i].activity_unit );
      BOOST_CHECK( a[i].extra_columns == b[i].extra_columns );
      BOOST_CHECK( a[i].fixed_column_names == b[i].fixed_column_names );
    }
  };

  for( const char *name : { "batch_decay_simple.csv", "batch_decay_product_value_unit.csv", "probeResults.csv" } )
  {
    const string csv = read_file( SpecUtils::append_path( batch_decay_dir(), name ) );
    string tsv = csv;
    std::replace( begin(tsv), end(tsv), ',', '\t' );

    BOOST_TEST_MESSAGE( "Comparing CSV and TSV parse of " << name );
    same( parse_csv( csv ), parse_csv( tsv ) );
    BOOST_CHECK( is_candidate_file( tsv, true ) );
  }

  // An empty cell must not shift the later columns, in either format.
  for( const string delim : { string(","), string("\t") } )
  {
    const string text = "Product" + delim + "Value" + delim + "Unit" + delim + "Notes\n"
                        "Cs137" + delim + "5" + delim + delim + "no unit, so Bq\n";
    const vector<BatchNuclide> keyed = parse_csv( text );
    BOOST_REQUIRE_EQUAL( keyed.size(), 1u );
    BOOST_CHECK_CLOSE( keyed[0].activity, 5.0 * PhysicalUnits::becquerel, 1.0E-9 );
  }

  // A comma inside a tab-separated cell is just text.
  const vector<BatchNuclide> commas =
    parse_csv( "Site\tProduct\tValue\tUnit\tNotes\n"
               "TeamA-01_1\tCs137 fallout, 2 days\t2\tuCi\tnorth, near road\n"
               "TeamA-01_2\tI131 fallout, 2 days\t3\tuCi\tnorth, near road\n" );
  BOOST_REQUIRE_EQUAL( commas.size(), 2u );
  BOOST_CHECK_EQUAL( commas[0].location, commas[1].location );
  BOOST_CHECK_EQUAL( commas[0].product_suffix, string(" fallout, 2 days") );
  BOOST_CHECK_CLOSE( commas[1].activity, 3.0E-6 * PhysicalUnits::curie, 1.0E-9 );
  BOOST_REQUIRE_EQUAL( commas[0].extra_columns.size(), 1u );
  BOOST_CHECK_EQUAL( commas[0].extra_columns[0].second, string("north, near road") );
}


// ---------------------------------------------------------------------------------------------
// Recognizing batch-decay input dropped on the main window (see SpecMeasManager's classifier).
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( candidate_file_detection )
{
  init_data_dirs();

  for( const char * const name : { "batch_decay_simple.csv", "batch_decay_product_value_unit.csv",
                                   "probeResults.csv" } )
  {
    const string path = SpecUtils::append_path( batch_decay_dir(), name );
    const string contents = read_file( path );

    // The whole file, and the 1024-byte start the drop classifier sees.
    BOOST_CHECK_MESSAGE( is_candidate_file( contents, true ), string(name) + " not recognized" );
    BOOST_CHECK_MESSAGE( is_candidate_file( contents.substr( 0, 1024 ), (contents.size() <= 1024) ),
                        string(name) + " not recognized from its first 1024 bytes" );

    // A drop only reaches the batch-decay check after the spectrum parse fails, so SpecUtils must not
    //  take these for spectra.
    SpecUtils::SpecFile spec;
    BOOST_CHECK_MESSAGE( !spec.load_file( path, SpecUtils::ParserType::Auto, "csv" ),
                        string(name) + " was parsed as a spectrum file" );
  }//for( each fixture )

  // A hand-written list can be a single line - but only a complete one.
  BOOST_CHECK( is_candidate_file( "Cs137, 1 uCi", true ) );
  BOOST_CHECK( !is_candidate_file( "Cs137, 1 uCi", false ) );
  BOOST_CHECK( is_candidate_file( "Cs137, 1 uCi\nCo60, 2", false ) );

  // A bare "nuclide, number" list could as well be nuclides and energies, so a drop only claims a
  //  simple list that gives a unit somewhere - though the tool's own upload still takes it.
  BOOST_CHECK( !is_candidate_file( "Cs137, 5\nCo60, 3\n", true ) );
  BOOST_CHECK_EQUAL( parse_csv( "Cs137, 5\nCo60, 3\n" ).size(), 2u );
  BOOST_CHECK( is_candidate_file( "Product,Value\nCs137,5\n", true ) );   // the header says what it is

  const string gammas = read_file( SpecUtils::append_path( g_data_dir, "CharacteristicGammas.txt" ) );
  BOOST_CHECK( !is_candidate_file( gammas, true ) );
  BOOST_CHECK( !is_candidate_file( gammas.substr( 0, 1024 ), false ) );

  // Nothing to decay, so nothing to open the tool for.
  BOOST_CHECK( !is_candidate_file( "B10, 5 uCi\n", true ) );

  // Other text must not be claimed.
  BOOST_CHECK( !is_candidate_file( "Channel,Counts\n0,5\n1,7\n2,9\n", true ) );
  BOOST_CHECK( !is_candidate_file( "Energy (keV), Counts\n10.1, 5\n", true ) );
  BOOST_CHECK( !is_candidate_file( "Some notes about Cs137, which is 1 uCi\n", true ) );
  BOOST_CHECK( !is_candidate_file( string( "Cs137, 1 uCi\n\0\n", 15 ), true ) );
  BOOST_CHECK( !is_candidate_file( "", true ) );

  vector<string> others;
  for( const char * const dir : { "det_eff", "gadras_detectors", "ceelo_drf", "cascade_truth" } )
  {
    for( const string &path : SpecUtils::recursive_ls( SpecUtils::append_path( g_test_data_dir, dir ), ".csv" ) )
      others.push_back( path );
  }
  others.push_back( SpecUtils::append_path( g_test_data_dir, "AnalystTests/beta_check_cases.csv" ) );
  others.push_back( SpecUtils::append_path( g_test_data_dir, "manual_rel_eff/source.txt" ) );
  BOOST_CHECK( others.size() > 20u );

  for( const string &path : others )
  {
    const string contents = read_file( path );
    BOOST_CHECK_MESSAGE( !is_candidate_file( contents, true ), "Claimed " + path );
    BOOST_CHECK_MESSAGE( !is_candidate_file( contents.substr( 0, 1024 ), (contents.size() <= 1024) ),
                        "Claimed the start of " + path );
  }
}


// ---------------------------------------------------------------------------------------------
// The two simpler CSV formats must be untouched by all of the above.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( simple_formats_unchanged )
{
  const vector<BatchNuclide> simple =
    parse_csv( read_file( SpecUtils::append_path( batch_decay_dir(), "batch_decay_simple.csv" ) ) );
  BOOST_REQUIRE( !simple.empty() );
  for( const BatchNuclide &in : simple )
  {
    BOOST_CHECK( in.nuclide );
    BOOST_CHECK_MESSAGE( in.location.empty(), "A 'nuclide, activity' file has no locations" );
  }

  const vector<BatchNuclide> keyed =
    parse_csv( read_file( SpecUtils::append_path( batch_decay_dir(),
                                                "batch_decay_product_value_unit.csv" ) ) );
  BOOST_REQUIRE( !keyed.empty() );
  for( const BatchNuclide &in : keyed )
  {
    BOOST_CHECK( in.nuclide );
    BOOST_CHECK_MESSAGE( in.location.empty(),
                        "A Product/Value/Unit file with no leading name column has no locations" );
  }

  // Both still produce the wide format: a label column plus one column per step.
  BatchDecayOptions opts = options_for( "1y" );
  opts.num_steps = 3;

  const BatchDecayResult result = DecayBatchCalc::decay( simple, opts );
  BOOST_CHECK_EQUAL( result.column_headers.size(), 1u + opts.num_steps );
  BOOST_CHECK_EQUAL( result.column_headers.front(), string("Nuclide") );
  for( const vector<string> &row : result.rows )
    BOOST_CHECK_EQUAL( row.size(), 1u + opts.num_steps );
}


// ---------------------------------------------------------------------------------------------
// 12. The activities parse_csv hands out must survive being rendered as text and read back.
//
// DecayBatchCalcWidget shows each input in a row's edit box and re-parses that text when it computes,
// so any rendering that loses digits changes the measurement.  `printToBestActivityUnits` counts
// *decimal places*, not significant figures, which zeroed 152 of these 1400 rows; the widget uses
// "%.9G" in the file's own unit instead (see `activity_text_in_unit`).  This pins the property that
// fix depends on, at the granularity the GUI needs it.
// ---------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( activity_text_round_trip )
{
  const vector<BatchNuclide> inputs =
    parse_csv( read_file( SpecUtils::append_path( batch_decay_dir(), "probeResults.csv" ) ) );
  BOOST_REQUIRE_EQUAL( inputs.size(), 1400u );

  size_t num_checked = 0;
  for( const BatchNuclide &in : inputs )
  {
    BOOST_REQUIRE( !in.activity_unit.empty() );
    const double unit = PhysicalUnits::stringToActivity( "1" + in.activity_unit );
    BOOST_REQUIRE( unit > 0.0 );

    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.9G", in.activity / unit );

    const double read_back = PhysicalUnits::stringToActivity( string(buffer) + in.activity_unit );

    if( in.activity == 0.0 )
    {
      BOOST_CHECK_EQUAL( read_back, 0.0 );
    }else
    {
      // A non-zero measurement must never come back as zero, and must keep its value to the 7
      //  significant figures the file supplies.
      BOOST_REQUIRE_MESSAGE( read_back > 0.0,
                            in.nuclide_str + " activity " + std::to_string(in.activity)
                            + " became zero when rendered as '" + buffer + "'" );
      BOOST_CHECK_CLOSE( read_back, in.activity, 1.0E-5 );
    }

    ++num_checked;
  }//for( each input )

  BOOST_CHECK_EQUAL( num_checked, 1400u );
}


// Data directories are found lazily by init_data_dirs(), on first use by a test case; doing it from a
//  BOOST_GLOBAL_FIXTURE instead aborts the whole run on any BOOST_REQUIRE there.
