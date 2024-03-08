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
#include <random>
#include <string>
#include <vector>
#include <cstdlib>


//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testPhysicalUnits_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include <boost/regex.hpp>


#include "InterSpec/PhysicalUnits.h"


using namespace std;
using namespace boost::unit_test;

using PhysicalUnits::stringToDistance;
using PhysicalUnits::stringToTimeDuration;
using PhysicalUnits::stringToTimeDurationPossibleHalfLife;
using PhysicalUnits::stringToActivity;
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


BOOST_AUTO_TEST_CASE( testStringToActivity ) {
  map<string,double> validCases;
  validCases["1.0 mci"] = 0.001*curie;
  validCases["1.0 mCi"] = 0.001*curie;
  validCases["128 bq"] = 128*becquerel;
  validCases["128 Mbq"] = 128*becquerel*1.0E6;
  validCases["128 becquerel"] = 128*becquerel;
  validCases["128 mbecquerel"] = 0.128*becquerel;
  validCases["128 millibecquerel"] = 0.128*becquerel;
  validCases["128 milli-becquerel"] = 0.128*becquerel;
  validCases["128 milli becquerel"] = 0.128*becquerel;
  validCases["128 milli bq"] = 0.128*becquerel;
  validCases["1.0 curie"] = curie;
  validCases["1.0 mcurie"] = 0.001*curie;
  validCases["1.0 m-curie"] = 0.001*curie;
  validCases["1.0 milli-curie"] = 0.001*curie;
  validCases["1.0 milli curie"] = 0.001*curie;
  validCases["128 mbq"] = 0.128*becquerel;
  validCases["128 ci"] = 128*curie;
  validCases["128. ci"] = 128*curie;
  validCases["128.ci"] = 128*curie;
  validCases["128.0ci"] = 128*curie;
  validCases["128.0mci"] = 0.128*curie;
  validCases["128.0Mci"] = 128*curie*1.0E6;
  validCases["128mbq"] = 0.128*becquerel;
  validCases["0.128mbq"] = 0.000128*becquerel;
  validCases["+0.128mbq"] = 0.000128*becquerel;
  validCases["+0.128mbq"] = 0.000128*becquerel;
  validCases["+0.128m-bq"] = 0.000128*becquerel;
  validCases["+0.128\t mbq"] = 0.000128*becquerel;
  validCases["+0.128\t milli-bq"] = 0.000128*becquerel;
  validCases["1.27e-06ci"] = 1.27e-06 * curie;
  validCases["1.27E+03mci"] = 1.27*curie;


  //need to make this more complete 
  
  for( map<string,double>::const_iterator iter = validCases.begin();
       iter != validCases.end(); ++iter )
  {
    try
    {
      const double expected = iter->second;
      const double parsed = stringToActivity( iter->first );
      BOOST_CHECK_CLOSE( parsed, expected, 1.0E-6 * fabs(expected) );
      BOOST_CHECK_CLOSE_FRACTION( parsed, expected, 1.0E-6 );
    }catch( std::exception &e )
    {
      BOOST_CHECK_MESSAGE( false, "Failed to parse '" << iter->first << "' using stringToActivity, received: " << e.what() << "." );
    }
  }
 
}

BOOST_AUTO_TEST_CASE( testStringToActivityShouldFail ) {
  const char * invalidCases[] =
    {
      "a33", "m33", "3", "5 micro-gram", "5 ugram", "5 microgram"
    };
    //Need to make this list more complete!
  const size_t ncases = sizeof(invalidCases) / sizeof(invalidCases[0]);


  for( size_t i = 0; i < ncases; ++i )
  {
    //BOOST_CHECK_THROW( stringToTimeDuration(*iter), std::exception );
    try
    {
      stringToActivity( invalidCases[i] );
      BOOST_CHECK_MESSAGE( false, "Did not fail to parse '" << invalidCases[i] << "' using stringToActivity, as expected." );
    }catch( std::exception &e )
    {
      BOOST_CHECK( true );
    }
  }
}


BOOST_AUTO_TEST_CASE( testStringToDistance ) {
  using PhysicalUnits::cm;
  using PhysicalUnits::m;
  const double inch = 2.54*PhysicalUnits::cm;
  const double foot = 12*inch;
  
  BOOST_CHECK_EQUAL( stringToDistance(" +0 "), 0.0 );
  BOOST_CHECK_EQUAL( stringToDistance(" 0"), 0.0 );
  BOOST_CHECK_EQUAL( stringToDistance("0"), 0.0 );
  BOOST_CHECK_EQUAL( stringToDistance("+0 "), 0.0 );
  BOOST_CHECK_EQUAL( stringToDistance("-0"), 0.0 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("+8.3E-02ft .9cm"), 8.3E-02*foot + 0.9*cm, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("3m"), 3.0*m, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("+3.2m"), 3.2*m, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("3'4\""), 3*foot + 4*inch, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("3'4.5\""), 3*foot + 4.5*inch, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("3.2'4.5\""), 3.2*foot + 4.5*inch, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("8ft 9cm"), 8*foot + 9*cm, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("3in 4ft"), 4*foot + 3*inch, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("3 in + 4  ft"), 4*foot + 3*inch, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("1.2E-10m"), 1.2E-10 * m, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance(".1m"), 0.1*m, 1.0E-6 );
  BOOST_CHECK_CLOSE( stringToDistance("0m"), 0.0, 1.0E-12 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("1m 2ft 3' 18cm"), 1*m + 2*foot + 3*foot + 18*cm, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("1000000.1 '"), 1000000.1*foot, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance(".3m"), 0.3*m, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToDistance("8m + 9cm"), 8*m + 9*cm, 1.0E-6 );

  const string shouldFailDist[] ={
    "8", "+8", "8m 3", "m8", "e.3m", "+-0.3m", "+-0.3m", "", "*1m", "1m ?",
    "1m +- 3cm", "1m \\xC2\\xB1 3cm", "8a1m", "0 4m"
  };
  
  for( const string &val : shouldFailDist )
  {
    BOOST_CHECK_THROW( stringToDistance( val ), std::exception );
  }
}//BOOST_AUTO_TEST_CASE( testStringToDistance )


BOOST_AUTO_TEST_CASE( testStringToTimeDuration ) {
  map<string,double> validCases;

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
  


  for( map<string,double>::const_iterator iter = validCases.begin();
       iter != validCases.end(); ++iter )
  {
    try
    {
      const double expected = iter->second;
      const double parsed = stringToTimeDuration( iter->first );
      BOOST_CHECK_CLOSE( parsed, expected, 1.0E-6 * fabs(expected) );
      BOOST_CHECK_CLOSE_FRACTION( parsed, expected, 1.0E-6 );
    }catch( std::exception &e )
    {
      BOOST_CHECK_MESSAGE( false, "Failed to parse '" << iter->first << "' using stringToTimeDuration, recieved: " << e.what() << "." );
    }
  }

  BOOST_CHECK_CLOSE_FRACTION( stringToTimeDuration( "1s", 1.0 ), 1.0, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToTimeDuration( "1m", 1.0E-12 ), 60*1.0E-12, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToTimeDuration( "1m", 1.0E+12 ), 60*1.0E+12, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToTimeDuration( "1s 3h", 88.323 ), (3*3600 + 1)*88.323, 1.0E-6 );
  BOOST_CHECK_CLOSE_FRACTION( stringToTimeDuration( "0", 99.9 ), 0.0, 1.0E-6 );
  BOOST_CHECK_THROW( stringToTimeDuration( "1m", -1 ), std::exception );
  BOOST_CHECK_THROW( stringToTimeDuration( "1m", -0 ), std::exception );
  BOOST_CHECK_THROW( stringToTimeDuration( "1m", +0 ), std::exception );
}


BOOST_AUTO_TEST_CASE( testStringToTimeDurationShouldFail ) {
  const char * invalidCases[] =
    {
      "a33", "m33", "33m 5", "3 19m", "3 3", "3 9 m sec", "minute", "second",
      "1 secd", "1.0.0 s", ".0E-6ss", "+5 + y", "23:59:59.100 10",
      "23:59:59.100m", "23:59:59.100 m", "a23:59:59.100", "5 23:59:59.100",
      "3hl", "2.3 hl", ".01 hl", "3halflives", "3 halflives", "3half lives",
      "P", "P5YT", "P-5Y" "PT-3h"
    };
  const size_t ncases = sizeof(invalidCases) / sizeof(invalidCases[0]);


  for( size_t i = 0; i < ncases; ++i )
  {
    //BOOST_CHECK_THROW( stringToTimeDuration(*iter), std::exception );
    try
    {
      stringToTimeDuration( invalidCases[i] );
      BOOST_CHECK_MESSAGE( false, "Did not fail to parse '" << invalidCases[i] << "' using stringToTimeDuration, as expected." );
    }catch( std::exception &e )
    {
      BOOST_CHECK( true );
    }
  }
}


BOOST_AUTO_TEST_CASE( testStringToTimeDurationPossibleHalfLife ) {

   //double stringToTimeDurationPossibleHalfLife( std::string str, double halflife, double sec_def = becquerel );

   typedef std::tuple<string,double,double> StrHlVal;
   vector<StrHlVal> cases;
   cases.push_back( StrHlVal("5hl", 0.289*millisecond, 5*0.289*millisecond ) );
   cases.push_back( StrHlVal(".5hl", 8*year, 4*year ) );
   cases.push_back( StrHlVal("+.5hl", 8*year, 4*year ) );
   cases.push_back( StrHlVal("+0.5hl", 8*year, 4*year ) );
   cases.push_back( StrHlVal("-.5hl", 8*year, -4*year ) );
   cases.push_back( StrHlVal("-0.5hl", 8*year, -4*year ) );
   cases.push_back( StrHlVal("5 half lives", 1.3*second, 5*1.3*second) );
   cases.push_back( StrHlVal("5 half lives 8min", 1.3*second, 5*1.3*second + 8*minute) );
   cases.push_back( StrHlVal("8min 5hl", 1.3*second, 5*1.3*second + 8*minute) );
   cases.push_back( StrHlVal("-5hl + 8min", 1.3*second, -5*1.3*second + 8*minute) );
   cases.push_back( StrHlVal("5hl + 2hl", 1.3*second, 7*1.3*second) );
   cases.push_back( StrHlVal("5hl  2.3hl", 1.3*second, 7.3*1.3*second) );
   cases.push_back( StrHlVal("5hl  .3hl", 1.3*second, 5.3*1.3*second) );
   cases.push_back( StrHlVal(".5hl  .3hl", 1.3*second, 0.8*1.3*second) );
   cases.push_back( StrHlVal("8y \t5m 3.4s ", 1.3*second, 8*year + 5*minute + 3.4*second) );
   


   for( const StrHlVal &val : cases  )
   {
     const string input = get<0>(val);
     const double hl = get<1>(val);
     const double expected = get<2>(val);

     try
     {
       const double parsed = stringToTimeDurationPossibleHalfLife( input, hl );
       BOOST_CHECK_CLOSE( parsed, expected, 1.0E-6 * fabs(expected) );
       BOOST_CHECK_CLOSE_FRACTION( parsed, expected, 1.0E-6 );
     }catch( std::exception &e )
     {
       BOOST_CHECK_MESSAGE( false, "Failed to parse '" << input << "' using stringToTimeDurationPossibleHalfLife, recieved: " << e.what() << "." );
     }
   }
}


BOOST_AUTO_TEST_CASE( userInputDistanceRegexTest ) {

  boost::regex distanceRegex( PhysicalUnits::sm_distanceRegex,
                           boost::regex::ECMAScript | boost::regex::icase );
  boost::regex distanceUnitOptionalRegex( PhysicalUnits::sm_distanceUnitOptionalRegex,
                           boost::regex::ECMAScript | boost::regex::icase );
  boost::regex uncertRegex( PhysicalUnits::sm_distanceUncertaintyRegex,
                           boost::regex::ECMAScript | boost::regex::icase );
  //sm_distanceUncertaintyUnitsOptionalRegex
     
  const string shouldMatchDist[] = {
    "3m", "3'4\"", "3in 4ft", "8m + 9cm", "8ft 9cm", "8ft 4m", "8ft4m", 
    "8ft+4m", "8 ft+4 m", "8\tft\t+\t4.3\tm", "1.2E-10m", "1.2E-10 m", 
    "1.2e-10 m", "1.2e+10feet", "1.2e+1ft 3.1cm", "3.1cm 1.2e+1ft", ".3m",
    ".3E-1m", "1000000.1 '", "1m", ".1m", "1.m", "0", "1m 2ft 3' 18cm",
    "1m + 2ft +3' +18cm", "+1m+2ft+3'+18cm"
  };

  const string shouldFailDist[] ={
    "8", "+8", "8m 3", "m8", "e.3m", "+-0.3m", "+-0.3m", "", "*1m", "1m ?", 
    "1m +- 3cm", "1m \\xC2\\xB1 3cm", "8a1m" 
  };

  const string shouldFailDistOpt[] = {
    "8a", "mm", "6.mmmm", "1 1", ".1 1", "1 .1", "1m .1"
  };

  const string noUnitShouldMatch[] = {
    "0.2", ".1", ".1E-6", "3.2", " 3.2 ", " 3.2", "3.2 ", "1.1E-6", "1.E-6",
    "0.E-6", "0.3e+6", ".3e+6"
  };
  
  const string uncertStringMatch[] = {
    "1m (\xC2\xB1 3cm)", "1 (\xC2\xB1 3) m",
    "1m (\xC2\xB1 .3cm)", ".1 (\xC2\xB1 .3) m",
    "1 +- 3 m", "0.93m +- 1cm", "0.93m+-1cm",
    "0.93+-1cm", "0.93\xC2\xB1 1.3 cm", 
    "0.93 \xC2\xB1 1.3 cm", ".93 \xC2\xB1 .3 cm"
  };
  
  const string uncertStringMatchShouldFail[] = {
    "1m\xC2\xB1", "1m (\xC2\xB1 .3cm) m", "1 ft (\xC2\xB1 .3)", "..1 ft",
    "1.a ft", "1 + - 3 ft", "1 - + 3 ft", "1 -/+ 3 ft"
  };
  
  for( const string str : shouldMatchDist )
  {
    BOOST_CHECK_MESSAGE( boost::regex_match( str, distanceRegex ), 
                         "sm_distanceRegex failed to match '" << str << "'" );
    BOOST_CHECK_MESSAGE( boost::regex_match( str, distanceUnitOptionalRegex ), 
                         "sm_distanceUnitOptionalRegex failed to match '" << str << "'" );
  }

  for( const string str : noUnitShouldMatch )
  {
    BOOST_CHECK_MESSAGE( boost::regex_match( str, distanceUnitOptionalRegex ), 
                         "sm_distanceUnitOptionalRegex failed to match '" << str << "'" );
  }

  for( const string str : shouldFailDist )
  {
    BOOST_CHECK_MESSAGE( !boost::regex_match( str, distanceRegex ), 
                         "sm_distanceRegex erroneously matched '" << str << "'" );
  }
  
  for( const string str : shouldFailDistOpt )
  {
    BOOST_CHECK_MESSAGE( !boost::regex_match( str, distanceUnitOptionalRegex ), 
                         "sm_distanceUnitOptionalRegex erroneously matched '" << str << "'" );
  }
  
   for( const string uncert : uncertStringMatch )
   {
      BOOST_CHECK_MESSAGE( boost::regex_match( uncert, uncertRegex ), 
                           "sm_distanceUncertaintyRegex failed to match '"
                            << (uncert) << "'" );
   }
   
   for( const string uncert : uncertStringMatchShouldFail )
   {
      BOOST_CHECK_MESSAGE( !boost::regex_match( uncert, uncertRegex ), 
                           "sm_distanceUncertaintyRegex erroneously matched '"
                            << (uncert) << "'" );
   }
  
}


BOOST_AUTO_TEST_CASE( userInputDurationRegexTest ) {
  boost::regex testregex( PhysicalUnits::sm_timeDurationRegex,
                           boost::regex::ECMAScript | boost::regex::icase );
  boost::regex testregexHLOpt( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex,
                           boost::regex::ECMAScript | boost::regex::icase );
  
  //could probably add even more formats to the checks.
  const string shouldMatchNoHl[] = {
    "3h", "4m 3 hour", "8s", "01:22:13.2", "8s 01:22:13.2", 
    "4m3h 01:22:13.2 8s", "01:22:13", "01:22:1.2", "+01:22:1.2", 
    "4m3h 01:22:13 8s", "4m3h 01:22:13.1 8s",
    "+4m + 3h", "01:22:13.1", 
    "+4m + 3h + 01:22:13.1 + 8s", "+32423:22:1.4334"
  };
  
  const string shouldMatchHl[] = {
    "3h", ".2hl", "+9.3 halflife", "0.1 halflives", "1half-life", 
    "9.9\thalf-lives", "8 half lives", "9 half life", "3h 4hl",
    ".3h 4hl", " 4hl + .3h ", " .1hl + .3h "
  };

  const string shouldFailNoHl[] = {
    "3hl", "8 01:22:13.2", "*8 01:22:13.2", "8m - 01:22:13.2",
    "3m 2hl", "3m + 2half lives", ".3hl", "2 2m", ".1 .99s"
    "01:22:13.2 hl", "8m - 01:22:13.2 hl"
  };
  
  const string shouldFailHl[] = {
    "3", "hl", "2:22 m", "-1m", "-01:22:33.1", "1h - 01:22:33.1",
    "1 01:22:33.1", "1 1hl"
  };

  for( const string str : shouldMatchNoHl ) 
  {
    BOOST_CHECK_MESSAGE( boost::regex_match( str, testregex ), 
                         "sm_timeDurationRegex failed to match '" << str << "'" );
    BOOST_CHECK_MESSAGE( boost::regex_match( str, testregexHLOpt ), 
                         "sm_timeDurationHalfLiveOptionalRegex failed to match '" << str << "'" );
  }
  
  for( const string str : shouldMatchHl )
  {
    BOOST_CHECK_MESSAGE( boost::regex_match( str, testregexHLOpt ), 
                         "sm_timeDurationHalfLiveOptionalRegex failed to match '" << str << "'" );
  }

  for( const string str : shouldFailNoHl )
  {
    BOOST_CHECK_MESSAGE( !boost::regex_match( str, testregex ), 
                         "sm_timeDurationRegex erroneously matched '" << str << "'" );
  }
  
  for( const string str : shouldFailHl )
  {
    BOOST_CHECK_MESSAGE( !boost::regex_match( str, testregex ), 
                         "sm_timeDurationRegex erroneously matched '" << str << "'" );
    BOOST_CHECK_MESSAGE( !boost::regex_match( str, testregexHLOpt ), 
                         "sm_timeDurationHalfLiveOptionalRegex erroneously matched '" << str << "'" );
  }
}


BOOST_AUTO_TEST_CASE( userInputActivityRegexTest ) {

  boost::regex testregex( PhysicalUnits::sm_activityRegex,
                           boost::regex::ECMAScript | boost::regex::icase );
  boost::regex testregexUnitOpt( PhysicalUnits::sm_activityUnitOptionalRegex,
                           boost::regex::ECMAScript | boost::regex::icase );

  const string shouldMatch[] = {
    "5bq", ".5bq", "1 bq", "1.2 becquerel", "0.2 ci", "+1cu", ".1curie", "1c",
    "5mbq", "5kbq", "5gbq",
    "5ubq", "5pbq", "5nbq", "5millibq", "5microbq", "5picobq", "5nanobq", 
    "5kilobq", "5megabq", "5terra-bq",
    "5.0u-bq", "5.p-bq", "5.3n-bq", "5.1milli-bq", "5micro-bq", "5pico-bq", 
    "5 nano bq", "5.2kilo bq", ".5 mega bq", "0.35 terra bq",
    "5.0 u bq", "5. p bq", "5.3 n bq", "5.1milli bq", "5micro ci", "5 pico bq", 
    "5 nano  bq", "5.2 kilo\tbq", ".5 mega\t bq", "0.35 terra \tbq",
    "5 bq", "+5 bq", "+.5 bq", "+0.5 bq",
    "1.4E-6 C", "1.4E-6 bq", "1.4E-6 mbq"
  };

  const string shouldFailUnits[] ={
    "90a", "bq 90", "bq.1", "1 1bq", "1 .1bq", "1 .1 bq", "1bq - 0.1 bq", "1",
    "5 bq 2 c", "1cu 3", "1cu bq", "1cu mbq", "1.4E-6", "1.4E-6a"
  };

  const string shouldFailUnitsOptional[] ={
    "90a", "bq 90", "bq.1", "1 1bq", "1 .1bq", "1 .1 bq", "1bq - 0.1 bq",
    "5 bq 2 c", "1cu 3", "1cu bq", "1cu mbq", "1.4E-6a"
  };
  
  const string noUnitShouldMatch[] = {
    ".2", "1.3", "1.4E-6", ".2", ".2E+10"
  };

  for( const string str : shouldMatch )
  {
    BOOST_CHECK_MESSAGE( boost::regex_match( str, testregex ), 
                         "sm_activityRegex failed to match '" << str << "'" );
    BOOST_CHECK_MESSAGE( boost::regex_match( str, testregexUnitOpt ), 
                         "sm_activityUnitOptionalRegex failed to match '" << str << "'" );
  }
  
  for( const string str : noUnitShouldMatch )
  {
    BOOST_CHECK_MESSAGE( boost::regex_match( str, testregexUnitOpt ), 
                         "sm_activityUnitOptionalRegex failed to match '" << str << "'" );
  }

  for( const string str : shouldFailUnits )
  {
    BOOST_CHECK_MESSAGE( !boost::regex_match( str, testregex ), 
                         "sm_activityRegex erroneously matched '" << str << "'" );
  }
  
  for( const string str : shouldFailUnitsOptional )
  {
    BOOST_CHECK_MESSAGE( !boost::regex_match( str, testregexUnitOpt ),
                        "sm_activityUnitOptionalRegex erroneously matched '" << str << "'" );
  }
}


BOOST_AUTO_TEST_CASE( test_stringToMass )
{
  using PhysicalUnits::gram;
  
  const double ounce = 28.3495*gram;
  const double pound = 453.592*gram;
  const double stone = 6350.29*gram;
  const double grain = 0.0647989*gram;
  
  const pair<string,double> valid_values[] = {
    {"1.2 gram", 1.2*gram},
    {"3g", 3*gram},
    {"5 kg", 5000*gram},
    {"5 mg", 0.005*gram},
    {"5 micro-gram", 5E-6*gram},
    {"5 ugram", 5E-6*gram},
    {"5 microgram", 5E-6*gram},
    {"5 micro gram", 5E-6*gram},
    {"5 kilogram", 5000*gram},
    {"5.0E0 kilo-gram", 5000*gram},
    {"2lb", 2*pound},
    {"2.34 pounds", 2.34*pound},
    {".1 oz", 0.1*ounce},
    {"0.1 oz", 0.1*ounce},
    {"+0.1 oz", 0.1*ounce},
    {"+ 0.1 oz", 0.1*ounce},
    {"2lb 3oz", 2*pound + 3*ounce},
    {"2.1lb 3E-01oz", 2.1*pound + 0.3*ounce},
    {"0oz", 0},
    {"2kg - 1g", 2000*gram - gram},
    {"-1.3E-01 kg + 10lbs", -130*gram + 10*pound},
    {"-5kg", -5000*gram}
  };
  
  
  auto my_check_close = []( const double parsed, const double expected ) -> bool {
    return (fabs(parsed - expected) <= (1.0E-6*fabs(expected)));
  };
  
  for( const auto &val : valid_values )
  {
    const double expected = val.second;
    double parsed;
    BOOST_REQUIRE_NO_THROW( parsed = PhysicalUnits::stringToMass( val.first ) );
    BOOST_CHECK_MESSAGE( my_check_close(parsed,expected), "Mass string '" << val.first
                        << "' should have value " << expected << ", but stringToMass returned "
                        << parsed << "\n");
  }//for( const auto &val : valid_values )
  
  
  const string invalid_values[] = {
    "13.2 gams",
    "13.2",
    "oz 13.2",
    "oz 1",
    "1.a gram",
    "1g 3",
    ""
  };
  
  for( const string &val : invalid_values )
  {
    BOOST_CHECK_THROW( PhysicalUnits::stringToMass(val), std::exception );
  }
}//BOOST_AUTO_TEST_CASE( test_stringToMass )


BOOST_AUTO_TEST_CASE( test_printToBestSpecificActivityUnits ) {
//Note 20200127: test cases are not exaustive, and output subject to change if string formatting 
//  of printToBestSpecificActivityUnits changes (which is likely)
  using namespace std;
  using namespace PhysicalUnits;
  const bool useCurry = true;
  const double actUnit = useCurry ? PhysicalUnits::curie : PhysicalUnits::becquerel;
  const double massUnit = PhysicalUnits::gram;
  const int nSigFigs = 3;

  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-18*actUnit/massUnit,nSigFigs,useCurry) == "1.23e-18 Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-16*actUnit/massUnit,nSigFigs,useCurry) == "1.23e-16 Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-15*actUnit/massUnit,nSigFigs,useCurry) == "1.23e-15 Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-14*actUnit/massUnit,nSigFigs,useCurry) == "12.3 pCi/kg" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-12*actUnit/massUnit,nSigFigs,useCurry) == "1.23 pCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-11*actUnit/massUnit,nSigFigs,useCurry) == "12.3 pCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-10*actUnit/massUnit,nSigFigs,useCurry) == "123 nCi/kg" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-9*actUnit/massUnit,nSigFigs,useCurry) == "1.23 nCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-8*actUnit/massUnit,nSigFigs,useCurry) == "12.3 nCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-7*actUnit/massUnit,nSigFigs,useCurry) == "123 uCi/kg" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-6*actUnit/massUnit,nSigFigs,useCurry) == "1.23 uCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-5*actUnit/massUnit,nSigFigs,useCurry) == "12.3 uCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-4*actUnit/massUnit,nSigFigs,useCurry) == "123 mCi/kg" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-3*actUnit/massUnit,nSigFigs,useCurry) == "1.23 mCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-2*actUnit/massUnit,nSigFigs,useCurry) == "12.3 mCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E-1*actUnit/massUnit,nSigFigs,useCurry) == "123 mCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*actUnit/massUnit,nSigFigs,useCurry) == "1.23 Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E1*actUnit/massUnit,nSigFigs,useCurry) == "12.3 Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E2*actUnit/massUnit,nSigFigs,useCurry) == "123 Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E3*actUnit/massUnit,nSigFigs,useCurry) == "1.23 kCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E4*actUnit/massUnit,nSigFigs,useCurry) == "12.3 kCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E5*actUnit/massUnit,nSigFigs,useCurry) == "123 kCi/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E6*actUnit/massUnit,nSigFigs,useCurry) == "1.23 mega-Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E7*actUnit/massUnit,nSigFigs,useCurry) == "12.3 mega-Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E8*actUnit/massUnit,nSigFigs,useCurry) == "123 mega-Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E9*actUnit/massUnit,nSigFigs,useCurry) == "1.23 giga-Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E10*actUnit/massUnit,nSigFigs,useCurry) == "12.3 giga-Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E12*actUnit/massUnit,nSigFigs,useCurry) == "1.23 tera-Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E14*actUnit/massUnit,nSigFigs,useCurry) == "123 tera-Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E15*actUnit/massUnit,nSigFigs,useCurry) == "1.23e+15 Ci/g" );
  BOOST_CHECK( printToBestSpecificActivityUnits(1.23456*1.0E16*actUnit/massUnit,nSigFigs,useCurry) == "1.23e+16 Ci/g" );
}//printToBestSpecificActivityUnits

