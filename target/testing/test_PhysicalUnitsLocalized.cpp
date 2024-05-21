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


#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PhysicalUnitsLocalized.h"


using namespace std;
using namespace boost::unit_test;

using PhysicalUnitsLocalized::stringToDistance;
using PhysicalUnitsLocalized::stringToTimeDuration;



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
