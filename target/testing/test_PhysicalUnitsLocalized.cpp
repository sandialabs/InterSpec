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
#include <map>
#include <tuple>
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>

#include <Wt/WMessageResourceBundle>


//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testPhysicalUnits_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include <boost/regex.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PhysicalUnitsLocalized.h"


using namespace std;
using namespace boost::unit_test;

using Wt::WMessageResourceBundle;
using PhysicalUnitsLocalized::printToBestTimeUnits;
using PhysicalUnitsLocalized::stringToTimeDuration;
using PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife;
using PhysicalUnits::second;
using PhysicalUnits::minute;
using PhysicalUnits::hour;
using PhysicalUnits::day;
using PhysicalUnits::year;
using PhysicalUnits::picosecond;
using PhysicalUnits::nanosecond;
using PhysicalUnits::microsecond;
using PhysicalUnits::millisecond;
using PhysicalUnits::month;
using PhysicalUnits::curie;
using PhysicalUnits::becquerel;

string g_app_text_dir;

void set_app_text_dir()
{
  // We only need to initialize things once
  static bool s_have_set = false;
  if( s_have_set )
    return;
  
  s_have_set = true;
  
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;
  
  string datadir;
  
  for( int i = 1; i < argc; ++i )
  {
    cout << "Arg " << i << ": '" << argv[i] << "'" << endl;
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( g_app_text_dir, "%20", " " );
  
  // Search around a little for the app_text directory, if it wasnt specified
  const string possible_paths[] = {
    SpecUtils::append_path(datadir,".."), 
    SpecUtils::append_path(datadir,"../.."),
    ".", "../", "../../", "../../../", "/Users/wcjohns/rad_ana/InterSpec/"
  };
  
  for( const auto &d : possible_paths )
  {
    if( SpecUtils::is_file( SpecUtils::append_path(d, "InterSpec_resources/app_text/InterSpec_fr.xml") ) )
    {
      g_app_text_dir = SpecUtils::append_path(d, "InterSpec_resources/app_text/");
      break;
    }
  }//for( loop over candidate dirs )
  
  
  BOOST_REQUIRE_MESSAGE( !g_app_text_dir.empty(), "Error finding app_text dir" );
}//void set_app_text_dir()



BOOST_AUTO_TEST_CASE( printToBestTimeUnitsLocalized ) {
  set_app_text_dir();
  
  WMessageResourceBundle bundle;
  
  std::string path = SpecUtils::append_path(g_app_text_dir, "InterSpec_fr");
  
  const bool loadInMemory = true;
  bundle.use( path, loadInMemory );
  
  bundle.refresh();
  
  const set<string> keys = bundle.keys( Wt::WMessageResourceBundle::Default );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-label-seconds-short" ), "Missing seconds short label" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-label-minutes-short" ), "Missing seconds minutes label" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-label-hours-short" ), "Missing seconds hours label" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-label-days-short" ), "Missing seconds days label" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-label-years-short" ), "Missing seconds years label" );
  
  BOOST_CHECK_EQUAL( printToBestTimeUnits(1.2E-6*PhysicalUnits::second, 2, PhysicalUnits::second, bundle), "1.20 us" );
  BOOST_CHECK_EQUAL( printToBestTimeUnits(1.2*PhysicalUnits::second, 2, PhysicalUnits::second, bundle), "1.20 s" );
  BOOST_CHECK_EQUAL( printToBestTimeUnits(0.0015*PhysicalUnits::second, 2, PhysicalUnits::second, bundle), "1.50 ms" );
  BOOST_CHECK_EQUAL( printToBestTimeUnits(1*PhysicalUnits::year, 2, PhysicalUnits::second, bundle), "1.00 a" );
  BOOST_CHECK_EQUAL( printToBestTimeUnits(3.2*60*PhysicalUnits::second, 2, PhysicalUnits::second, bundle), "3.20 m" );
  BOOST_CHECK_EQUAL( printToBestTimeUnits(3.51*3600*PhysicalUnits::second, 2, PhysicalUnits::second, bundle), "3.51 h" );
  BOOST_CHECK_EQUAL( printToBestTimeUnits(3.51*3600*24*PhysicalUnits::second, 2, PhysicalUnits::second, bundle), "3.51 j" );
  
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.2*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.1E-6*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.09E-9*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.09E-12*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.3E-4*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.4E-2*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 3600*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 2*3600*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 2*2*3600*PhysicalUnits::second );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.2*PhysicalUnits::year );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.2E8*PhysicalUnits::year );
  PhysicalUnitsLocalized::printToBestTimeUnits( 1.2E3*PhysicalUnits::year );
}//BOOST_AUTO_TEST_CASE( printToBestTimeUnits )


BOOST_AUTO_TEST_CASE( testStringToTimeDurationLocalized ) {
  set_app_text_dir();
  
  WMessageResourceBundle bundle;
  
  std::string path = SpecUtils::append_path(g_app_text_dir, "InterSpec_fr");
  
  const bool loadInMemory = true;
  bundle.use( path, loadInMemory );
  
  bundle.refresh();
  
  const set<string> keys = bundle.keys( Wt::WMessageResourceBundle::Default );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-labels-second" ), "Missing seconds localization" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-labels-minute" ), "Missing minutes localization" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-labels-hours" ), "Missing hours localization" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-labels-days" ), "Missing days localization" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-labels-year" ), "Missing years localization" );
  BOOST_REQUIRE_MESSAGE( keys.count( "units-labels-half-lives" ), "Missing half-lives localization" );
  
  
  map<string,double> validCases;
  
  // Make sure English will always work
  validCases["5y"] = 5*year;
  validCases["5 y"] = 5*year;
  validCases["+5 y"] = 5*year;
  validCases["5 year"] = 5*year;
  validCases["5.year"] = 5*year;
  validCases["5 year 2d 3h 4s"] = 5*year + 2*day + +3*hour + 4*second;
  validCases["3s 5min"] = 3*second + 5*minute;
  validCases[".5s .25m"] = 0.5*second + 0.25*minute;
  validCases[".5 second .25 m"] = 0.5*second + 0.25*minute;
  validCases[" .3ps "] = 0.3*picosecond;
  validCases["8 picosecond"] = 8*picosecond;
  validCases["3 minute"] = 3*minute;
  validCases[".9999s"] = 0.9999*second;
  validCases["10.s"] = 10.0*second;
  validCases["10. s"] = 10.0*second;
  validCases["-0.9999s"] = -0.9999*second;
  validCases["1s - 0.1s"] = 0.9*second;
  validCases["1s - .1s"] = 0.9*second;
  validCases["1s - 1s"] = 0.0*second;
  validCases["-.9999s"] = -0.9999*second;
  validCases["1.0E-6s"] = 1*microsecond;
  validCases["1.0E-6s 0.1E-5s"] = 2*microsecond;
  validCases["88.3E+2s"] = (88.3E+2) * second;
  validCases["23:59:59.100"] = 23*hour + 59*minute + 59.1*second;
  validCases["23:59:59.100 + 10m"] = 23*hour + 69*minute + 59.1*second;
  validCases["10m 23:59:59.100"] = 23*hour + 69*minute + 59.1*second;
  validCases["10m 23:59:59.100"] = 23*hour + 69*minute + 59.1*second;
  validCases["0:1:59.9 + 0:00:1.1"] = 0*hour + 1*minute + 61*second;
  validCases["0:1:59.9 + 0:00:1.1 5s 8m"] = 0*hour + 9*minute + 66*second;
  validCases["-0:1:15.5 "] = -(1*minute + 15.5*second);
  validCases["1s - 0.9s"] = 0.1*second;
  validCases["-0.9s + 1s"] = 0.1*second;
  validCases["-0.9s  1s"] = 0.1*second;
  validCases["0s"] = 0.0;
  validCases["0"] = 0.0;
  validCases["PT1S"] = 1.0*second;
  validCases["1s - PT1S"] = 0.0;
  validCases["1y P1Y1D"] = 2*year + 1*day;
  validCases["3s -P1Y2M2W3DT10H30M10.1S 10y"] = 3.0*second - (1.0*year + 2*month + 2*7*day + 3*day + 10*hour + 30*minute + 10.1*second) + 10.0*year;
  validCases["PT10.1S 1y"] = 10.1*second + 1*year;
  validCases["P1.1YT0S"] = 1.1*year;
  validCases["-P1.1YT0S"] = -1.1*year;
  validCases["PT3H5M2.343S"] = 3*hour + 5*minute + 2.343*second;
  validCases["PT3m5h2s"] = 3*minute + 5*hour + 2*second;
  
  
  // French time periods.  Both pure French, and French mixed with English should work
  validCases["+18.1a"] = 18.1*PhysicalUnits::year;
  validCases["1.2E-3 year 2m -5s 14.5 ms +18.1y"] = 1.2E-3*PhysicalUnits::year + 2*60*second - 5*second + 0.0145*second + 18.1*year;
  validCases["1.2 an 3jours 5.1 heure"] = 1.2*PhysicalUnits::year + 3*PhysicalUnits::day + 5.1*PhysicalUnits::hour;
  validCases["2.11années"] = 2.11*PhysicalUnits::year;
  validCases["1h 3ms"] = 1*PhysicalUnits::hour + 3*millisecond;
  validCases["1.1E-3jour 1   heures -3.2min"] = 1.1E-3*day + 1*hour - 3.2*minute;
  validCases["1.1E-3 sec 3ms"] = 1.1E-3*second + 3*millisecond;
  validCases["1.1E-3années"] = 1.1E-3*year;


  for( map<string,double>::const_iterator iter = validCases.begin();
       iter != validCases.end(); ++iter )
  {
    try
    {
      const double expected = iter->second;
      const double parsed = stringToTimeDurationPossibleHalfLife( iter->first, -1.0, second, bundle );
      BOOST_CHECK_CLOSE( parsed, expected, 1.0E-6 * fabs(expected) );
      BOOST_CHECK_CLOSE_FRACTION( parsed, expected, 1.0E-6 );
    }catch( std::exception &e )
    {
      BOOST_CHECK_MESSAGE( false, "Failed to parse '" << iter->first << "' using stringToTimeDuration, recieved: " << e.what() << "." );
    }
  }
}


BOOST_AUTO_TEST_CASE( testStringToTimeDurationPossibleHalfLifeLocalized ) {

  set_app_text_dir();
  
  WMessageResourceBundle bundle;
  
  std::string path = SpecUtils::append_path(g_app_text_dir, "InterSpec_fr");
  
  const bool loadInMemory = true;
  bundle.use( path, loadInMemory );
  
  bundle.refresh();
  
  
  typedef std::tuple<string,double,double> StrHlVal;
  vector<StrHlVal> cases;
  
  cases.push_back( StrHlVal("5hl", 0.289*millisecond, 5*0.289*millisecond ) );
  cases.push_back( StrHlVal(".5hl", 8*year, 4*year ) );
  cases.push_back( StrHlVal("+0.5hl", 8*year, 4*year ) );
  cases.push_back( StrHlVal("-.5hl", 8*year, -4*year ) );
  cases.push_back( StrHlVal("-0.5hl", 8*year, -4*year ) );
  cases.push_back( StrHlVal("5 half lives", 1.3*second, 5*1.3*second) );
  cases.push_back( StrHlVal("5 half lives 8min", 1.3*second, 5*1.3*second + 8*minute) );
  cases.push_back( StrHlVal("8min 5hl", 1.3*second, 5*1.3*second + 8*minute) );
  cases.push_back( StrHlVal("-5hl + 8min", 1.3*second, -5*1.3*second + 8*minute) );
  cases.push_back( StrHlVal("8y \t5m 3.4s ", 1.3*second, 8*year + 5*minute + 3.4*second) );
  
  cases.push_back( StrHlVal("5 demi vie", 1.3*second, 5*1.3*second) );
  cases.push_back( StrHlVal("5dv", 1.3*second, 5*1.3*second) );
  cases.push_back( StrHlVal("5demi-vies", 1.3*second, 5*1.3*second) );
  cases.push_back( StrHlVal("5demi vie", 1.3*second, 5*1.3*second) );
  cases.push_back( StrHlVal("5demi vies", 1.3*second, 5*1.3*second) );
  cases.push_back( StrHlVal("1.1E3dv", 1.3*second, 1.1E3*1.3*second) );
  cases.push_back( StrHlVal("5dv 0.53jour", 1.3*second, 5*1.3*second + 0.53*24*3600*second) );
  
  for( const StrHlVal &val : cases  )
  {
    const string input = get<0>(val);
    const double hl = get<1>(val);
    const double expected = get<2>(val);
    
    try
    {
      const double parsed = stringToTimeDurationPossibleHalfLife( input, hl, second, bundle );
      BOOST_CHECK_CLOSE( parsed, expected, 1.0E-6 * fabs(expected) );
      BOOST_CHECK_CLOSE_FRACTION( parsed, expected, 1.0E-6 );
    }catch( std::exception &e )
    {
      BOOST_CHECK_MESSAGE( false, "Failed to parse '" << input << "' using localized stringToTimeDurationPossibleHalfLife, recieved: " << e.what() << "." );
    }
  }
}

