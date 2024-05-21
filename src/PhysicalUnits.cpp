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

#include <string>
#include <cctype>

#include <boost/regex.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/PhysicalUnits.h"

using namespace std;

namespace PhysicalUnits
{
#define MU_CHARACTER_1 "\xCE\xBC"
#define MU_CHARACTER_2 "\xC2\xB5"
#define DIAERESIS_O  "\xC3\xB6"

#define DECIMAL_REGEX "\\s*[\\+\\-]?\\s*((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?\\s*"
#define POS_DECIMAL_REGEX "\\s*\\+?\\s*((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?\\s*"
#define DIST_UNITS_REGEX "(meter|cm|km|mm|um|nm|m|ft|feet|'|inches|inch|in|\")"
#define METRIC_PREFIX_UNITS "m|M|k|g|G|t|T|u|" MU_CHARACTER_1 "|" MU_CHARACTER_2 "|p|n|f|milli|micro|pico|nano|femto|kilo|mega|giga|terra"

#define PLUS_MINUS_REGEX "(\\xC2?\\xB1|\\+\\-\\s*|\\-\\+\\s*)"
#define TIME_UNIT_REGEX "(years|year|yr|y|days|day|d|hrs|hours|hour|h|minute|minutes|min|m|second|seconds|sec|s|ms|microseconds|us|nanosecond|nanoseconds|ns)"
#define HALF_LIFE_REGEX "(hl|halflife|halflives|half-life|half-lives|half lives|half lifes|half live|halflive|half life)"
#define ACTIVITY_UNIT_REGEX "(bq|becquerel|ci|cu|curie|c)"
#define ABSORBED_DOSE_UNIT_REGEX "(gray|Gy|gy|rad|erg|erg\\/g|erg per gram)"
#define EQUIVALENT_DOSE_UNIT_REGEX "(sievert|Sv|rem|roentgen|r" DIAERESIS_O "entgen)"

#define DURATION_REGEX "(\\+?\\d+:\\d\\d:\\d+(\\.\\d+)?)"
#define DURATION_REGEX_POS_NEG "([\\+\\-]?\\d+:\\d\\d:\\d+(\\.\\d+)?)"
#define ISO_8601_DURATION_REGEX "[\\-+]?[Pp](?!$)(\\d+(?:\\.\\d+)?[Yy])?(\\d+(?:\\.\\d+)?[Mm])?(\\d+(?:\\.\\d+)?[Ww])?(\\d+(?:\\.\\d+)?[Dd])?([Tt](?=\\d)(\\d+(?:\\.\\d+)?[Hh])?(\\d+(?:\\.\\d+)?[Mm])?(\\d+(?:\\.\\d+)?[Ss])?)?"
  
const char * const sm_distanceRegex
    = "^(((\\s*[+\\-]?\\s*0\\s*)|((" POS_DECIMAL_REGEX ")\\s*" DIST_UNITS_REGEX ")+\\s*))+$";
const char * const sm_distanceUnitOptionalRegex = "^("
"(((\\s*" DECIMAL_REGEX ")\\s*" DIST_UNITS_REGEX ")+\\s*)"
"|(\\s*" DECIMAL_REGEX "\\s*)"
")$";


const char * const sm_distanceUncertaintyRegex
   = "^\\s*((" POS_DECIMAL_REGEX "(\\(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*\\s*\\)\\s*)?" DIST_UNITS_REGEX ")"
     "|(" POS_DECIMAL_REGEX "(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*)?" DIST_UNITS_REGEX ")"
     "|(" POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX ")?)"
     "|(" POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "(\\(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "\\s*\\))?)"
     "|(" POS_DECIMAL_REGEX DIST_UNITS_REGEX "\\s*(\\(" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "+\\s*\\)\\s*)?)"
     "|(" POS_DECIMAL_REGEX DIST_UNITS_REGEX "\\s*(" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "+\\s*)?)"
      ")+\\s*$";
  
  
const char * const sm_distanceUncertaintyUnitsOptionalRegex
  = "^\\s*((" POS_DECIMAL_REGEX "(\\(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*\\s*\\)\\s*)?" DIST_UNITS_REGEX "?)"
  "|(" POS_DECIMAL_REGEX "(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*)?" DIST_UNITS_REGEX ")"
  "|(" POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "?(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "?)?)"
  "|(" POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "?(\\(\\s*" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "?\\s*\\))?)"
  "|(" POS_DECIMAL_REGEX DIST_UNITS_REGEX "?\\s*(" PLUS_MINUS_REGEX POS_DECIMAL_REGEX "\\s*" DIST_UNITS_REGEX "?\\s*)?)"
  ")+\\s*$";
  
  //Non escaped version of regex is: "^(\s*\+?\d+(\.\d*)?(?:[Ee][+-]?\d+)?\s*(?:m|M|k|g|G|t|T|u|p|n|f|milli|micro|pico|nano|femto|kilo|mega|giga|terra)?[\s-]*([Bb][Qq]|[Bb][Ee][Cc][Qq][Uu][Ee][Rr][Ee][Ll]|[Cc][Ii]|[Cc][Uu]|[Cc][Uu][Rr][Ii][Ee]))?\s*$"
const char * const sm_activityRegex = "^(\\s*" POS_DECIMAL_REGEX "\\s*(?:(?:" METRIC_PREFIX_UNITS ")?[\\s-]*" ACTIVITY_UNIT_REGEX "))?\\s*$";
                                       

const char * const sm_absorbedDoseRegex = "^(\\s*" POS_DECIMAL_REGEX "\\s*(?:(?:" METRIC_PREFIX_UNITS ")?[\\s-]*" ABSORBED_DOSE_UNIT_REGEX "))?\\s*$";

const char * const sm_equivalentDoseRegex = "^(\\s*" POS_DECIMAL_REGEX "\\s*(?:(?:" METRIC_PREFIX_UNITS ")?[\\s-]*" EQUIVALENT_DOSE_UNIT_REGEX "))?\\s*$";

const char * const sm_activityUnitOptionalRegex
            = "^(\\s*" POS_DECIMAL_REGEX "\\s*"
              "(?:(?:" METRIC_PREFIX_UNITS ")?"
              "[\\s-]*"
              "([Bb][Qq]|[Bb][Ee][Cc][Qq][Uu][Ee][Rr][Ee][Ll]|[Cc][Ii]|[Cc]|[Cc][Uu]|[Cc][Uu][Rr][Ii][Ee]))?)?"
              "\\s*$";
  
const char * const sm_timeDurationRegex
  = "^(\\s*\\+?\\s*((" POS_DECIMAL_REGEX "\\s*" TIME_UNIT_REGEX "\\s*)|(" DURATION_REGEX "\\s*)\\s*)|(" ISO_8601_DURATION_REGEX "\\s*))+$";
  
// Note: (at least as of Wt 3.7.1) the RegEx validator expects the first match group to match the
//       entire string, not just a subset - so for example with an input string of "2.5 hl", we need
//       to keep this from matching with "2.5 h" (which given the half-live optional component is
//       last in the regex, will happen), so we force the the regex to have start at beginning of
//       input string, and go to the end of the string, which I guess we should do for any of the
//       WRegExpValidator's.
const char * const sm_timeDurationHalfLiveOptionalRegex
= "^(\\s*\\+?\\s*((" POS_DECIMAL_REGEX "\\s*" TIME_UNIT_REGEX "\\s*)"
  "|(" DURATION_REGEX "\\s*)"
  "|(" POS_DECIMAL_REGEX "\\s*" HALF_LIFE_REGEX "\\s*)"
  "|(0+(\\.0*)?)" //Allow the user to enter "0", or "0.000", with no time units
  "))+$";

const char * const sm_timeDurationHalfLiveOptionalPosOrNegRegex
   = "^\\s*((" DECIMAL_REGEX "\\s*" TIME_UNIT_REGEX "\\s*)"
     "|(" DURATION_REGEX_POS_NEG "\\s*)"
     "|(" DECIMAL_REGEX "\\s*" HALF_LIFE_REGEX "\\s*))+$";



const char * const sm_positiveDecimalRegex = POS_DECIMAL_REGEX;
  
const char * const sm_decimalRegex = DECIMAL_REGEX;

const UnitNameValuePairV sm_lengthUnitNameValues{
  {"nm", nm},
  {"um", um},
  {"mm", mm},
  {"cm", cm},
  {"m", meter},
  {"km", 1.0E3*meter}
};
  
const UnitNameValuePairV sm_lengthUnitHtmlNameValues{
  {"nm", nm},
  {"&mu;m", um},
  {"mm", mm},
  {"cm", cm},
  {"m", meter},
  {"km", 1.0E3*meter}
};
    
  
const UnitNameValuePairV sm_activityUnitNameValues{
  {"uBq", 1.0E-6*bq},
  {"mBq", 1.0E-3*bq},
  {"bq", bq},
  {"kBq", kBq},
  {"MBq", MBq},
  {"GBq", GBq},
  {"TBq", TBq},
  {"fCi", fCi},
  {"pCi", pCi},
  {"nCi", nCi},
  {"microCi", microCi},
  {"mCi", mCi},
  {"ci", ci}
};

const UnitNameValuePairV sm_activityUnitHtmlNameValues{
  {"&mu;Bq", 1.0E-6*bq},
  {"mBq", 0.001*bq},
  {"bq", bq},
  {"kBq", kBq},
  {"MBq", MBq},
  {"GBq", GBq},
  {"TBq", TBq},
  {"fCi", fCi},
  {"pCi", pCi},
  {"nCi", nCi},
  {"&mu;Ci", microCi},
  {"mCi", mCi},
  {"ci", ci}
};

const UnitNameValuePairV sm_doseRateUnitHtmlNameValues{
  {"&mu;R/hr", (rem*1.0E-6)/hour},
  {"mR/hr", (rem*1.0E-3)/hour},
  {"R/hr", rem/hour},
  {"&mu;Sv/hr", (sievert*1.0E-6)/hour},
  {"mSv/hr", (sievert*1.0E-3)/hour},
  {"Sv/hr", sievert/hour}
};

namespace
{
  double metrix_prefix_value( const string &prefix )
  {
    if( prefix == "" )
      return 1.0;
    if( prefix == "u" || prefix == "micro" || prefix == MU_CHARACTER_1 || prefix == MU_CHARACTER_2 )
      return 1.0E-6;
    if( prefix == "m" || prefix == "milli" )
      return 1.0E-3;
    if( prefix == "k" || prefix == "kilo" )
      return 1.0E3;
    if( prefix == "M" || prefix == "mega" )
      return 1.0E6;
    if( prefix == "g" || prefix == "giga" )
      return 1.0E9;
    if( prefix == "t" || prefix == "terra" )
      return 1.0E12;
    if( prefix == "n" || prefix == "nano" )
      return 1.0E-9;
    if( prefix == "p" || prefix == "pico" )
      return 1.0E-12;
    
    throw std::runtime_error( "Unrecognized prefix: '" + prefix + "'" );
  }//double metrix_prefix_value( const string &prefix )

  /** Function that takes in a metric prefix string, and returns its numeric value. */
  double prefix_str_to_val( const string &letters )
  {
    if( letters.empty() )
      return 1.0;
    
    if( SpecUtils::istarts_with(letters, "f" ) )
      return 1.0E-15;
    
    if( SpecUtils::istarts_with(letters, "p" ) )
      return 1.0E-12;
    
    if( SpecUtils::istarts_with(letters, "n" ) )
      return 1.0E-9;
    
    if( SpecUtils::istarts_with(letters, "u" )
             || SpecUtils::istarts_with(letters, "micro" )
             || SpecUtils::starts_with(letters, MU_CHARACTER_1 )
             || SpecUtils::starts_with(letters, MU_CHARACTER_2 ) )
      return 1.0E-6;
    
    if( SpecUtils::starts_with(letters, "m" )
             || SpecUtils::istarts_with(letters, "milli" ) )
      return 1.0E-3;
    
    if( SpecUtils::istarts_with(letters, "k" ) )
      return 1.0E+3;
    
    if( SpecUtils::starts_with(letters, "M" )
             || SpecUtils::istarts_with(letters, "mega" ) )
      return 1.0E+6;
    
    if( SpecUtils::starts_with(letters, "G" ) )
      return 1.0E+9;
    
    if( SpecUtils::starts_with(letters, "T" ) )
      return 1.0E+12;
    
    throw std::runtime_error( "'" + letters + "' is not a recognized prefix unit" );
  }//double prefix_str_to_val( const string &letters )
}//namespace



std::string printToBestLengthUnits( double length, int maxNpostDecimal,
                                    const double cm_definition )
{
  using namespace std;
  length *= cm / cm_definition;

  char formatflag[32], buffer[64];
  snprintf(formatflag, sizeof(formatflag), "%%.%if %%s", maxNpostDecimal );
  
  if( length < (1000.0 * nm) )
    snprintf(buffer, sizeof(buffer), formatflag, (length/nm), "nm" );
  else if( length < (1000.0 * um) )
    snprintf(buffer, sizeof(buffer), formatflag, (length/um), "um" );
  else if( length < (10.0 * mm) )
    snprintf(buffer, sizeof(buffer), formatflag, (length/mm), "mm" );
  else if( length < (100.0 * cm) )
    snprintf(buffer, sizeof(buffer), formatflag, (length/cm), "cm" );
  else if( length < (1000.0 * m) )
    snprintf(buffer, sizeof(buffer), formatflag, (length/m), "m" );
  else
    snprintf(buffer, sizeof(buffer), formatflag, (length/(1000.0*m)), "km" );

  return buffer;
}//std::string printToBestLengthUnits( ... )

std::wstring printToBestLengthUnits( double length, double uncert,
                                     int maxNpostDecimal,
                                     const double cm_definition )
{
  using namespace std;
  length *= cm / cm_definition;
  uncert *= cm / cm_definition;
  
  double unit = 1000.0*m;
  wstring unitLabel = L" km";
  
  if( length < (1000.0 * nm) )
  {
    unit = nm;
    unitLabel = L" nm";
  }else if( length < (1000.0 * um) )
  {
    unit = um;
    unitLabel = L" um";
  }else if( length < (10.0 * mm) )
  {
    unit = mm;
    unitLabel = L" mm";
  }else if( length < (100.0 * cm) )
  {
    unit = cm;
    unitLabel = L" cm";
  }else if( length < (1000.0 * m) )
  {
    unit = m;
    unitLabel = L" m";
  }
  
  wstringstream txt;
  txt << setprecision(maxNpostDecimal) << fixed << (length/unit) << L" (\x00B1";
  if( ((uncert > 100.0*length) || (uncert < length/100.0)) && uncert!=0.0 )
    txt << setprecision(maxNpostDecimal) << scientific << (uncert/unit) << ")";
  else
    txt << setprecision(maxNpostDecimal) << fixed << (uncert/unit) << ")";
  txt << unitLabel;

  return txt.str();
}//std::string printToBestLengthUnits( ... )
  
  
std::string printToBestActivityUnits( double activity,
                                      int maxNpostDecimal,
                                      bool useCurries,
                                      double bq_definition )
{
  using namespace std;
  activity *= becquerel / bq_definition;

  char formatflag[32], buffer[64];
  const char *unitstr = useCurries ? "Ci" : "Bq";
  
  snprintf(formatflag, sizeof(formatflag), "%%.%if %%s%s", maxNpostDecimal, unitstr );
  

  if( useCurries )
    activity /= curie;
  else
    activity /= becquerel;

  if( activity == 0 )
    snprintf(buffer, sizeof(buffer), "0 %s", unitstr );
  else if( activity < 1.0E-12 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E15), "f" );
  else if( activity < 1.0E-9 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E12), "p" );
  else if( activity < 1.0E-6 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E9), "n" );
  else if( activity < 1.0E-3 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E6), "u" );
  else if( activity < 1.0 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E3), "m" );
  else if( activity < 1.0E3 )
    snprintf(buffer, sizeof(buffer), formatflag, activity, "" );
  else if( activity < 1.0E6 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E-3), "k" );
  else if( activity < 1.0E9 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E-6), "M" );
  else if( activity < 1.0E12 )
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E-9), "G" );
  else
    snprintf(buffer, sizeof(buffer), formatflag, (activity*1.0E-12), "T" );
  
  return buffer;
}//std::string printToBestActivityUnits(...)


std::string printToBestTimeUnits( double time,
                                  int maxDecimal,
                                  double second_definition )
{
  using namespace std;
  time *= second / second_definition;

  char formatflag[32], buffer[64];
  snprintf(formatflag, sizeof(formatflag), "%%.%if %%s", maxDecimal );
  
  if( fabs(time) < second*1.0E-3 )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time*1.0E6/second), "us" );
  }else if( fabs(time) < second*0.1 )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time*1.0E3/second), "ms");
  }else if( fabs(time) < minute )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/second), "s");
  }else if( fabs(time) > 1.0E3*year )
  {
    snprintf(formatflag, sizeof(formatflag), "%%.%ig y", maxDecimal );
    snprintf(buffer, sizeof(buffer), formatflag, (time/year) );
  }else if( fabs(time) > 0.99*year )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/year), "y");
  }else if( fabs(time) > 3.0*day )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/day), "d");
  }else if( fabs(time) > 3.0*hour )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/hour), "h");
  }else if( fabs(time) > 3.0*minute )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/minute), "m");
  }else if( fabs(time) > 1.0*second )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/second), "s");
  }else if( fabs(time) > 1.0*millisecond )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/millisecond), "ms");
  }else if( fabs(time) > 1.0*microsecond )
  {
//#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
    const unsigned char utf8mus[] = { 0xCE, 0xBC, 0x73, 0 };
    snprintf(buffer, sizeof(buffer), formatflag, (time/microsecond), (char *)utf8mus );
//#else
//    snprintf(buffer, sizeof(buffer), formatflag, (time/microsecond), "&mu;s" );
//#endif
  }else if( fabs(time) > 1.0*picosecond )
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/picosecond), "ps");
  }else
  {
    snprintf(buffer, sizeof(buffer), formatflag, (time/second), "s");
  }
  /*
  else
  {
    //Less than 4 days we'll return time in the format 'hh:mm:ss.ss'
    const int nhour = static_cast<int>( floor( time / hour ) );
    time -= nhour*hour;
    const int nminute = static_cast<int>( floor( time / minute ) );
    time -= nminute*minute;
    const int nsec = static_cast<int>( floor( time / second ) );
    time -= nsec*second;

    //Need to round the fractional seconds to nearest number of ticks
    typedef boost::posix_time::time_duration::tick_type tick_type;
    const double mult = std::pow(10.0,static_cast<double>(maxDecimal));
    tick_type ticks = boost::posix_time::time_duration::ticks_per_second();
    const int count = static_cast<int>( floor(time*mult*ticks+0.5) / mult );
  
    boost::posix_time::time_duration duration( nhour, nminute, nsec, count );
    string durstr = boost::posix_time::to_simple_string( duration );
  
    //Get rid of unwanted decimals (they should be all zeros from rounding above)
    const size_t pos = durstr.find_last_of( '.' );
    if( pos != string::npos )
      durstr = durstr.substr( 0, min(size_t(pos+maxDecimal+1), durstr.size()) );
    
    snprintf(buffer, sizeof(buffer), "%s", durstr.c_str() );
  }
   */
  
  return buffer;
}//printToBestTimeUnits(...)

  
std::string printToBestMassUnits( const double mass,
                                   const int maxDecimal,
                                   const double gram_definition )
{
  const double mass_grams = mass / gram_definition;
  
  if( IsInf(mass_grams) || IsNan(mass_grams) )
    return "Nan";
  
  char formatflag[32], massstr[64];
  snprintf(formatflag, sizeof(formatflag), "%%.%if %%s", maxDecimal );
  
  if( mass_grams < 1.0E-9 )
    snprintf( massstr, sizeof(massstr), formatflag, (1.0E12 * mass_grams), "pg" );
  else if( mass_grams < 1.0E-6 )
    snprintf( massstr, sizeof(massstr), formatflag, (1.0E9 * mass_grams), "ng" );
  else if( mass_grams < 1.0E-3 )
    snprintf( massstr, sizeof(massstr), formatflag, (1.0E6 * mass_grams), "ug" );
  else if( mass_grams < 1.0E0 )
    snprintf( massstr, sizeof(massstr), formatflag, (1.0E3 * mass_grams), "mg" );
  else if( mass_grams < 1.0E3 )
    snprintf( massstr, sizeof(massstr), formatflag, (1.0E0 * mass_grams), "g" );
  else// if( mass_grams < 1.0E6 )
    snprintf( massstr, sizeof(massstr), formatflag, (1.0E-3 * mass_grams), "kg" );
  //      else
  //      snprintf( massstr, sizeof(massstr), "%.4g g", mass_grams );

  return massstr;
}//std::string printToBestMassUnits(...)
  
  
  
string printToBestEquivalentDoseRateUnits( double dose, const int maxNpostDecimal,
                             const bool useSievert, const double sievertPerHourDefinition )
{
  dose *= (sievert/hour) / sievertPerHourDefinition;
  
  char formatflag[32], buffer[64];
  const char *unitstr = useSievert ? "sv/hr" : "rem/hr";
  
  snprintf(formatflag, sizeof(formatflag), "%%.%if %%s%s", maxNpostDecimal, unitstr );
  
  
  if( useSievert )
    dose /= (sievert/hour);
  else
    dose /= (rem/hour);
  
  if( dose < 1.0E-6 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E9), "n" );
  else if( dose < 1.0E-3 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E6), "u" );
  else if( dose < 1.0 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E3), "m" );
  else if( dose < 1.0E3 )
    snprintf(buffer, sizeof(buffer), formatflag, dose, "" );
  else if( dose < 1.0E6 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E-3), "k" );
  else
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E-6), "M" );
  
  return buffer;
}//string printToBestEquivalentDoseRateUnits(...)


std::string printToBestEquivalentDoseUnits( double dose, const int maxNpostDecimal,
                                            const bool use_sievert,
                                            const double sievert_definition )
{
  dose *= sievert / sievert_definition;
  
  char formatflag[32], buffer[64];
  const char *unitstr = use_sievert ? "sv" : "rem";
  
  snprintf(formatflag, sizeof(formatflag), "%%.%if %%s%s", maxNpostDecimal, unitstr );
  
  dose /= (use_sievert ? sievert : rem);
  
  if( dose < 1.0E-6 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E9), "n" );
  else if( dose < 1.0E-3 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E6), "u" );
  else if( dose < 1.0 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E3), "m" );
  else if( dose < 1.0E3 )
    snprintf(buffer, sizeof(buffer), formatflag, dose, "" );
  else if( dose < 1.0E6 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E-3), "k" );
  else
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E-6), "M" );
  
  return buffer;
}//printToBestEquivalentDoseUnits(...)


std::string printToBestAbsorbedDoseUnits( double dose, const int maxNpostDecimal,
                                          const bool use_gray, const double gray_definition )
{
  dose *= gray / gray_definition;
  
  char formatflag[32], buffer[64];
  const char *unitstr = use_gray ? "Gy" : "rad";
  
  snprintf(formatflag, sizeof(formatflag), "%%.%if %%s%s", maxNpostDecimal, unitstr );
  
  dose /= (use_gray ? gray : rad);
  
  if( dose < 1.0E-6 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E9), "n" );
  else if( dose < 1.0E-3 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E6), "u" );
  else if( dose < 1.0 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E3), "m" );
  else if( dose < 1.0E3 )
    snprintf(buffer, sizeof(buffer), formatflag, dose, "" );
  else if( dose < 1.0E6 )
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E-3), "k" );
  else
    snprintf(buffer, sizeof(buffer), formatflag, (dose*1.0E-6), "M" );
  
  return buffer;
}//printToBestEquivalentDoseUnits(...)



std::string printToBestSpecificActivityUnits( double actPerMass,
                                               const int numSignificantFigures,
                                               const bool useCurie )
{
  //This is something barely more than a minimum viable function for printing
  //  out specific activity.  Issues remaining that should be checked at least
  //  are:
  //  - The activity units are modified first before mass units (e.g., kCi/g
  //    output instead of Ci/kg).  Is this the convention people normally use?
  //    Or is there a standardized convention for this?
  //  - Rounding has only kinda been tested on macOS, and not Linux or Windows
  //  - The mega, giga, and tera activity prefixes are currently fully spelled
  //    out, and all other prefixes us just one letter abreviations (k, m, u, n,
  //    p) - this should probably be made consistent.
  //  - Could use some better tests implemented.
  
  char formatflag[64], buffer[64];
  string numerator_unitstr = useCurie ? "Ci" : "Bq";
  string denominator_unitstr = "g";
  
  actPerMass /= (useCurie ? curie : becquerel);
  actPerMass *= gram;
  
  //Get the exponent
  double order = std::floor( std::log10(actPerMass) );
  
  //Define some ranges on power of 10 of the input were we will try to adjust
  //  to reasonable units.  Outside this range we'll use scientific notation and
  //  Ci/g or Bq/g.
  const int min_order_adjust = -15;
  const int max_order_adjust = 15;
  
  double value = actPerMass;
  
  if( order <= min_order_adjust )
  {
    //Do nothing, we will put value in scientific notation in Ci/g
  }else if( order < -10 )
  {
    numerator_unitstr = "p" + numerator_unitstr;
    value /= 1.0E-12;
  }else if( order < -7 )
  {
    numerator_unitstr = "n" + numerator_unitstr;
    value /= 1.0E-9;
  }else if( order < -4 )
  {
    numerator_unitstr = "u" + numerator_unitstr;
    value /= 1.0E-6;
  }else if( order < 0 )
  {
    numerator_unitstr = "m" + numerator_unitstr;
    value /= 1.0E-3;
  }else if( order < 3 )
  {
    
  }else if( order < 6 )
  {
    numerator_unitstr = "k" + numerator_unitstr;
    value /= 1.0E3;
  }else if( order < 9 )
  {
    numerator_unitstr = "mega-" + numerator_unitstr;
    value /= 1.0E6;
  }else if( order < 12 )
  {
    numerator_unitstr = "giga-" + numerator_unitstr;
    value /= 1.0E9;
  }else if( order < max_order_adjust )
  {
    numerator_unitstr = "tera-" + numerator_unitstr;
    value /= 1.0E12;
  }
  
  order = std::floor( std::log10(value) );
  
  if( order <= min_order_adjust )
  {
    //Do nothing, we will put value in scientific notation in Ci/g
  }else if( order < 0 )
  {
    denominator_unitstr = "kg";
    value /= 1.0E-3;
  }else if( order < 3 )
  {
    
  }else if( order < 6 )
  {
    denominator_unitstr = "mg";
    value /= 1.0E3;
  }else if( order < 9 )
  {
    denominator_unitstr = "ug";
    value /= 1.0E6;
  }else if( order < 12 )
  {
    denominator_unitstr = "ng";
    value /= 1.0E9;
  }else if( order < max_order_adjust )
  {
    denominator_unitstr = "pg";
    value /= 1.0E12;
  }
  
  if( order > min_order_adjust && order < max_order_adjust )
  {
    order = std::ceil( std::fabs(std::log10(value)) );
    const int ndecimals = std::max( numSignificantFigures - static_cast<int>(order), 0 );
    snprintf( formatflag, sizeof(formatflag), "%%.%if %s/%s",
              ndecimals, numerator_unitstr.c_str(), denominator_unitstr.c_str() );
  }else
  {
    value = actPerMass;
    numerator_unitstr = useCurie ? "Ci" : "Bq";
    denominator_unitstr = "g";
    snprintf( formatflag, sizeof(formatflag), "%%.%ig %s/%s",
             numSignificantFigures, numerator_unitstr.c_str(), denominator_unitstr.c_str() );
  }
  
  snprintf( buffer, sizeof(buffer), formatflag, value );
  
  return buffer;
}//printToBestSpecificActivityUnits(...)
  
  
double stringToTimeDuration( std::string str, double second_def )
{
  return stringToTimeDurationPossibleHalfLife( str, -1.0, second_def );
}//double stringToTimeDuration( std::string str, double second_def )
  

double stringToTimeDurationPossibleHalfLife( std::string str,
                                             const double hl,
                                             double second_def )
{
  double time_dur = 0.0;

  if( second_def <= 0.0 )
    throw runtime_error( "Second definition must be larger than zero" );
  
  SpecUtils::trim( str );
  SpecUtils::to_lower_ascii( str );
  
  while( str.find("+ ") != string::npos )
    SpecUtils::ireplace_all(str, "+ ", "+" );
  while( str.find("- ") != string::npos )
    SpecUtils::ireplace_all(str, "- ", "-" );

  const size_t colonpos = str.find(":");
  if( colonpos != std::string::npos )
  {
    size_t startpos, endpos;
    for( startpos = colonpos; startpos > 0; --startpos )
      if( isspace(str[startpos-1]) )
        break;
    for( endpos = colonpos; endpos < str.size(); ++endpos )
      if( isspace(str[endpos]) )
        break;
    
    
    string durstr = str.substr( startpos, endpos-startpos );
    str = str.substr(0,startpos) + str.substr(endpos);
    time_dur += SpecUtils::delimited_duration_string_to_seconds( durstr );
    
    //SpecUtils::trim( str );
    //boost::posix_time::time_duration dur;
    //dur = boost::posix_time::duration_from_string( durstr );
    //if( dur.is_special() )
    //  throw runtime_error( "stringToTimeDuration(...): couldnt translate '"
    //                        + durstr + "' to a format like '23:59:59.000'" );
    //time_dur = (1.0E-6 * dur.total_microseconds() * second_def);
    if( str.length() )
      time_dur += stringToTimeDurationPossibleHalfLife( str, hl, second_def );
    
    return time_dur;
  }//if( str.find(":") != std::string::npos )

  
  //if( (hl<=0.0) && (str.size()>2)
  //    && (str[0]=='p' || ((str[0]=='+' || str[0]=='-') && str[0]=='p') ) )
  if( str.find('p') != string::npos )
  {
    //extended format of ISO 8601, i.e., PnYnMnDTnHnMnS,
    //  ex 'P1Y2M2W3DT10H30M10.1S', '-P1347M', 'PT1M'
    
    boost::smatch matches;
    boost::regex expression( ISO_8601_DURATION_REGEX );
    if( boost::regex_search( str, matches, expression ) )
    {
      assert( matches.size() == 9 );
      //for( size_t i = 0; i < matches.size(); ++i )
      //  cout << "Match[" << i << "] = '" << matches[i] << "'" << endl;
      //For: 'P1Y2M2W3DT10H30M10.1S'
      //matches[0] = 'P1Y2M2W3DT10H30M10.1S'
      //matches[1] = '1y'
      //matches[2] = '2m'
      //matches[3] = '2w'
      //matches[4] = '3d'
      //matches[5] = 'T10H30M10.1S'
      //matches[6] = '10h'
      //matches[7] = '30m'
      //matches[8] = '10.1s'
      double iso_dur = 0.0;
      const int nfields = std::min( 9, static_cast<int>(matches.size()) );
      for( int i = 1; i < nfields; ++i )
      {
        string valstr = matches[i].str();
        if( i==5 || (valstr.size() <= 1) )
          continue;
        
        valstr = valstr.substr(0,valstr.size()-1);
        double val = std::stod( valstr );  //could throw, if our regex is messed up
        
        switch( i )
        {
          case 1: val *= PhysicalUnits::year;    break;
          case 2: val *= PhysicalUnits::month;   break;
          case 3: val *= 7.0*PhysicalUnits::day; break;
          case 4: val *= PhysicalUnits::day;     break;
          case 6: val *= PhysicalUnits::hour;    break;
          case 7: val *= PhysicalUnits::minute;  break;
          case 8: val *= PhysicalUnits::second;  break;
        }
        
        iso_dur += val;
      }//for( size_t i = 1; i < 10; ++i )
      
      const string fullmatch = matches[0].str();
      assert( fullmatch.size() > 0 );
      
      if( fullmatch[0] == '-' )
        iso_dur *= -1.0;
      
      if( fullmatch.size() > 1 ) //protect against regex just matching "p"
      {
        time_dur += iso_dur;
      
        if( str.size() == fullmatch.size() )
        return time_dur;
        
        const size_t startpos = str.find( fullmatch );
        if( startpos != string::npos )  //should always be true, but jic
        {
          str.erase( begin(str)+startpos, begin(str)+startpos+fullmatch.size() );
          SpecUtils::trim( str );
          if( !str.empty() )
            time_dur += stringToTimeDurationPossibleHalfLife( str, hl, second_def );
        }
        
        return time_dur;
      }//if( fullmatch.size() > 1 )
    }//if( boost::regex_match( str, matches, expression ) )
    //If not a required format, keep trying.
  }//if( str.size()>2 && (str[0]=='p' || ((str[0]=='+' || str[0]=='-') && str[0]=='p') )  )
  
  

  const string time_regex ="\\s*((\\+|\\-)?((\\d+(\\.\\d*)?)|(\\.\\d*))(?:[Ee][+\\-]?\\d+)?)"
                            "\\s*("
                            "ps|picosecond|pico\\-*sec|pico\\-*second"
                            "|ns|nano\\-*second|nano\\-*sec"
                            "|us|micro\\-*sec|micro\\-*second"
                            "|ms|milli\\-*sec|milli\\-*second"
                            "|second|seconds|sec|s|minute|minutes|min|m"
                            "|hours|hour|hrs|h|days|day|d|year|years|yr|yrs|y"
                            "|hl|halflife|halflives|half\\-life|half\\-lives|half\\slives|half\\slife"
                            ")"
                            "\\s*((\\+|\\-)?(((\\d+(\\.\\d*)?)|(\\.\\d*)).*))?";
  
  boost::smatch matches;
  boost::regex expression( time_regex );
  if( boost::regex_match( str, matches, expression ) )
  {
    //for( size_t i = 0; i < matches.size(); ++i )
      //cerr << "stringToTimeDuration Match " << i << " is " << string( matches[i].first, matches[i].second ) << endl;
    
    string number = string( matches[1].first, matches[1].second );
    string letters = string( matches[7].first, matches[7].second );
    SpecUtils::trim( number );
    SpecUtils::trim( letters );
    SpecUtils::to_lower_ascii( letters );
    
    double value;
    if( !(stringstream(number) >> value) )
      throw std::runtime_error( "Unable to convert '" + number + "' to a number" );
    
    double unit = 0.0;
    if( SpecUtils::starts_with(letters, "ps" )
       || SpecUtils::starts_with(letters, "pico" ) )
      unit = second * 1.0E-12;
    else if( SpecUtils::starts_with(letters, "ns" )
        || SpecUtils::starts_with(letters, "nano" ) )
      unit = second * 1.0E-9;
    else if( SpecUtils::starts_with(letters, "us" )
        || SpecUtils::starts_with(letters, "micro" )
        || SpecUtils::starts_with(letters, MU_CHARACTER_1 )
        || SpecUtils::starts_with(letters, MU_CHARACTER_2 ) )
      unit = second * 1.0E-6;
    else if( SpecUtils::starts_with(letters, "ms" )
        || SpecUtils::starts_with(letters, "mS" )
        || SpecUtils::starts_with(letters, "milli" ) )
      unit = second * 1.0E-3;
    else if( SpecUtils::starts_with(letters, "s" ) )
      unit = second;
    else if( (letters=="m") || SpecUtils::contains(letters, "min" )  )
      unit = minute;
    else if( letters=="h"
             || SpecUtils::contains(letters, "hour" )
             || SpecUtils::contains(letters, "hrs" ) )
      unit = hour;
    else if( letters=="d" || SpecUtils::contains(letters, "day" ) )
      unit = day;
    else if( letters=="y"
             || SpecUtils::contains(letters, "yr" )
             || SpecUtils::contains(letters, "year" ) )
      unit = year;
    else if( SpecUtils::starts_with(letters, "h" ) && (hl > 0.0) )  //the regex hopefully makes sure
      unit = hl;
    else
    {
      std::string msg = "cant identify label '";
      msg += letters + "' as a unit of time";
      //std::cerr << msg << endl;
      throw runtime_error( msg );
    }

    time_dur += value*unit*(second_def/second);

    if( matches[8].length() > 0 )
      time_dur += stringToTimeDurationPossibleHalfLife( matches[8], hl, second_def );
  }else
  {
    //If we are here, it may be the case the user didnt enter any units
    //  which is okay only if they mean for zero time, so well check
    //  for that case, and if not throw an exception.
    const string zero_regex ="(?|\\+|\\-)?0+(?|\\.0*)*";
    boost::smatch zeromatches;
    boost::regex xeroexpression( zero_regex );
    if( boost::regex_match( str, zeromatches, xeroexpression ) )
      time_dur = 0.0;
    else
      throw runtime_error( "stringToTimeDuration(...): could not translate '"
                                + str + "' to a time");
  }//if( match regex ) / else

  return time_dur;
}//double stringToTimeDuration( std::string str, double second_def )
  
  
//stringToActivity(...): Takes a user string, and attempts to interpret
//  if as an activity.  Ex. '52.3 bq', '88 MCi', etc
//  throws std::runtime_error on failure
double stringToActivity( std::string str, double bq_def )
{
  SpecUtils::trim( str );
  boost::smatch mtch;
  boost::regex expr( "(" POS_DECIMAL_REGEX ")" "\\s*([a-zA-Z \\-]+)" );
  
  /*
   // 20231113: A shorter and probably more robust implementation of this function is commented
   //           out here - it passes unit-tests, but I didnt enable yet incase unit-tests dont fully
   //           cover possibilities.
  boost::smatch test_mtch;
  boost::regex test_expr( "(" POS_DECIMAL_REGEX ")" "\\s*((" METRIC_PREFIX_UNITS ")*[\\s\\-]*(bq|ci|bec.*|cu.*))" );
  if( !boost::regex_match( str, test_mtch, test_expr ) )
    throw std::runtime_error( std::string("'") + str + "' is not appropriate activity string" );
  
  string number = SpecUtils::trim_copy( string( test_mtch[1].first, test_mtch[1].second ) );
  string prefix = SpecUtils::trim_copy( string( test_mtch[7].first, test_mtch[7].second ) );
  string units = SpecUtils::trim_copy( string( test_mtch[8].first, test_mtch[8].second ) );
    
  const double prefix_val = prefix_str_to_val( prefix );
    
  const bool hasb = (SpecUtils::icontains(units, "bq" ) || SpecUtils::icontains(units, "bec" ));
  const bool hasc = (SpecUtils::icontains(units, "ci" ) || SpecUtils::icontains(units, "cu" ));
    
  if( hasb == hasc )
    throw runtime_error( "didnt contain unit specification (bq or ci) in '" + units + "'" );
    
  return std::stod(number) * prefix_val * (hasb ? becquerel : curie) * bq_def / becquerel;
  */

  if( boost::regex_match( str, mtch, expr ) )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    if( mtch.size() != 7 )
    {
      char buffer[512];
      snprintf( buffer, sizeof(buffer),
               "stringToActivity(..) regex yields %i matches instead of 6",
                static_cast<int>(mtch.size()) );
      log_developer_error( __func__, buffer );
    }
#endif
    
    assert( mtch.size() == 7 );
    
    //cout << "mtch.size()=" << mtch.size() << endl;
    //for( size_t i = 0; i < mtch.size(); ++i )
    //  cout << "\t" << string( mtch[i].first, mtch[i].second ) << endl;
    //cout << endl;
    
    string number = string( mtch[1].first, mtch[1].second );
    string letters = string( mtch[6].first, mtch[6].second );
    
    const double value = std::stod( number );
    double unit = 0.0;
    
    if( SpecUtils::istarts_with(letters, "f" ) )
      unit = 1.0E-15;
    else if( SpecUtils::istarts_with(letters, "p" ) )
      unit = 1.0E-12;
    else if( SpecUtils::istarts_with(letters, "n" ) )
      unit = 1.0E-9;
    else if( SpecUtils::istarts_with(letters, "u" )
             || SpecUtils::istarts_with(letters, "micro" )
             || SpecUtils::starts_with(letters, MU_CHARACTER_1 )
             || SpecUtils::starts_with(letters, MU_CHARACTER_2 ) )
      unit = 1.0E-6;
    else if( SpecUtils::starts_with(letters, "m" )
             || SpecUtils::istarts_with(letters, "milli" ) )
      unit = 1.0E-3;
    else if( ((letters.size() == 1) && (letters[0]=='b' || letters[0]=='B' || letters[0]=='c' || letters[0]=='C'))
            || SpecUtils::istarts_with(letters, "bq")
            || SpecUtils::istarts_with(letters, "becq")
            || SpecUtils::istarts_with(letters, "ci")
            || SpecUtils::istarts_with(letters, "cu") )
      unit = 1.0;
    else if( SpecUtils::istarts_with(letters, "k" ) )
      unit = 1.0E+3;
    else if( SpecUtils::starts_with(letters, "M" )
             || SpecUtils::istarts_with(letters, "mega" ) )
      unit = 1.0E+6;
    else if( SpecUtils::starts_with(letters, "G" ) )
      unit = 1.0E+9;
    else if( SpecUtils::starts_with(letters, "T" ) )
      unit = 1.0E+12;
    else
      throw std::runtime_error( "'" + letters + "' is not a recognized activity unit" );
    
    //There is no need for a prefix
    //else
    //{
      //std::string msg = "didnt contain valid prefix specification in '";
      //msg += letters + "'";
      //std::cerr << msg << endl;
      //throw runtime_error( msg );
    //}

    const bool hasb = (SpecUtils::icontains(letters, "bq" )
                      || SpecUtils::icontains(letters, "bec" ));
    const bool hasc = (SpecUtils::icontains(letters, "ci" )
                       || SpecUtils::icontains(letters, "curi" )
                       || SpecUtils::icontains(letters, "cu" ) );

    if( hasb && !hasc )
     unit *= becquerel;
    else if( hasc && !hasb )
     unit *= curie;
    else
    {
      std::string msg = "didnt contain unit specification (bq or ci) in '";
      msg += letters + "'";
      throw runtime_error( msg );
    }

    return value * unit * bq_def / becquerel;
  }else
  {
    throw std::runtime_error( std::string("'") + str
                              + "' is not appropriate activity string" );
  } //if( boost::regex_match( substr, mtch, expr ) ) / else

  return 0.0;
}//double stringToActivity( std::string str, double bq_def )



double stringToDistance( std::string str, double cm_definition )
{
  double distance = 0.0;

  string::size_type openParan = str.find( '(' );
  while( openParan != string::npos )
  {
    const string::size_type closeParan = str.find( ')', openParan );
    if( closeParan == string::npos )
      break;
    str = str.substr( 0, openParan ) + str.substr( closeParan + 1 );
    openParan = str.find( '(', openParan );
  }//while( openParan != string::npos )
  
  const char * const regex_str
  = "\\s*(([+\\-]?\\s*0\\s*)|(((" POS_DECIMAL_REGEX ")\\s*(" DIST_UNITS_REGEX "))(\\s*" POS_DECIMAL_REGEX ".+)*\\s*))";
  
  boost::smatch matches;
  boost::regex expression( regex_str, boost::regex::ECMAScript|boost::regex::icase );

  if( !boost::regex_match( str, matches, expression ) )
  {
    char msg[128];
    snprintf(msg, sizeof(msg), "'%s' is an invalid distance", str.c_str() );
    //cerr << endl << msg << endl;
    throw std::runtime_error( msg );
  }//if( we dont have a match )

  //cout << "\nstringToDistance - str='" << str << "':\n";
  //for( size_t i = 0; i < matches.size(); ++i )
  //  cerr << "\tMatch " << i << " is '" << string( matches[i].first, matches[i].second ) << "'" << endl;
  
  string floatstr( matches[5].first, matches[5].second );
  string unitstr( matches[10].first, matches[10].second );
  const string leftover = SpecUtils::trim_copy( string(matches[12].first, matches[12].second) );
  SpecUtils::trim( floatstr );
  SpecUtils::ireplace_all(floatstr, " ", "");
  SpecUtils::ireplace_all(floatstr, "\t", "");
  SpecUtils::trim( unitstr );
  SpecUtils::to_lower_ascii( unitstr );
  
  // We will allow the input string of "0", with no units
  if( (matches[2].first != matches[2].second) && unitstr.empty() && leftover.empty() )
  {
    assert( matches[4].second == matches[4].first );
    assert( matches[5].second == matches[5].first );
    assert( matches[6].second == matches[6].first );
    //.... all the way through 16
    return 0.0;
  }

  double unitval = 0.0;
  if( unitstr == "meter" )       unitval = m;
  else if( unitstr == "cm" )     unitval = cm;
  else if( unitstr == "km" )     unitval = 1000.0*m;
  else if( unitstr == "mm" )     unitval = mm;
  else if( unitstr == "um" )     unitval = 0.001*mm;
  else if( unitstr == "nm" )     unitval = 0.000001*mm;
  else if( unitstr == "m" )      unitval = m;
  else if( unitstr == "ft" )     unitval = 12.0*2.54*cm;
  else if( unitstr == "'" )      unitval = 12.0*2.54*cm;
  else if( unitstr == "feet" )   unitval = 12.0*2.54*cm;
  else if( unitstr == "in" )     unitval = 2.54*cm;
  else if( unitstr == "\"" )     unitval = 2.54*cm;
  else if( unitstr == "i" )      unitval = 2.54*cm;
  else if( unitstr == "in" )     unitval = 2.54*cm;
  else if( unitstr == "inch" )   unitval = 2.54*cm;
  else if( unitstr == "inchs" )  unitval = 2.54*cm;
  else if( unitstr == "inches" ) unitval = 2.54*cm;
  else if( unitstr == "attoparsec" )  unitval = 3.086*cm;
  else if( unitstr == "attoparsecs" ) unitval = 3.086*cm;
  
//  else assert(0);

  unitval *= (cm_definition / cm);

  distance = std::stod( floatstr );  //shouldnt ever throw, right?
  distance *= unitval;

  //if there are characters past what we needed for a match, maybe they are
  //  another distance string, lets try to add them on
  if( leftover.length() > 1 )
  {
    try
    {
      distance += stringToDistance( leftover, cm_definition );
    }catch(...)
    {}
  }//if( string(matches[4]).length() )

//  cerr << str << " = " << distance/cm_definition << " cm, or "
//       << distance/cm_definition/2.54 << " inches" << endl;

  return distance;
}//double stringToDistance( std::string str )


double stringToAbsorbedDose( const std::string &str, const double gray_definition )
{
  /*
  //A poor-persons test for this function:
  std::cout << "Doses:" << std::endl;
  const std::vector<std::string> dosestr = {"1 rad", "1.123 gray", "1.1E-3 Gy", "1.1krad", "1urad",
    "10.2uGy", "10.1gy", "10.2uGy", "1.4gy 20 rad", "8\xCE\xBCgy", "8.1 \xCE\xBCgy", "8\xCE\xBCrad",
    "1.2 gray", "1.23gy", "100 urad", "100 mrad", "100 milli-rad", "100 millirad", "13 erg/g", "1 gray 1 Gy"
  };
  for( const std::string str : dosestr )
  {
    const double val = PhysicalUnits::stringToAbsorbedDose( str );
    std::cout << "\t\"" << str << "\" --> " << PhysicalUnits::printToBestAbsorbedDoseUnits(val,4,true)
              << " and " << PhysicalUnits::printToBestAbsorbedDoseUnits(val,4,false) << std::endl;
  }
  */
  
  const string regex_str ="\\s*\\+?(\\d+(\\.\\d*)?(?:[Ee][+\\-]?\\d+)?)"
                            "\\s*(" METRIC_PREFIX_UNITS ")*?"
                            "\\s*\\-*\\s*"
                            ABSORBED_DOSE_UNIT_REGEX
                            "\\s*(\\d.+)?";
  
  boost::smatch matches;
  boost::regex expression( regex_str, boost::regex::ECMAScript|boost::regex::icase );

  if( !boost::regex_match( str, matches, expression ) )
  {
    char msg[128];
    snprintf(msg, sizeof(msg), "'%s' is not an absorbed dose", str.c_str() );
    //cerr << endl << msg << endl;
    throw std::runtime_error( msg );
  }//if( we dont have a match )

  //cout << endl << endl;
  //for( size_t i = 0; i < matches.size(); ++i )
  //  cerr << "stringToDistance Match " << i << " is " << string( matches[i].first, matches[i].second ) << endl;
  
  string floatstr( matches[1].first, matches[1].second );
  string prefix( matches[3].first, matches[3].second );
  string unitstr( matches[4].first, matches[4].second );
  string trailingstr( matches[5].first, matches[5].second );
  SpecUtils::trim( floatstr );  //These trims may not be necassary, but whatever
  SpecUtils::trim( prefix );
  SpecUtils::trim( unitstr );
  SpecUtils::to_lower_ascii( unitstr );
  
  //cout << "\n\"" << str << "\" : {number->" << floatstr << ", prefix->" << prefix
  //     << ", unit->" << unitstr << ", trailing->" << trailingstr << "}" << endl;
  
  double unitval = 0.0;
  if( unitstr == "gray" || unitstr == "gy" )     unitval = gray;
  else if( unitstr == "rad" )                    unitval = rad;
  else if( SpecUtils::contains(unitstr, "erg") ) unitval = gray * 1.0E-04;
  else
    throw runtime_error( "Unexpeced absorbed dose unit: '" + unitstr + "'" );

  double dose = unitval * std::stod( floatstr );  //shouldnt ever throw, right?
  dose *= metrix_prefix_value( prefix );
  dose *= (gray_definition / gray);

  //if there are characters past what we needed for a match, maybe they are
  //  another distance string, lets try to add them on
  if( trailingstr.length() > 2 )
  {
    try
    {
      dose += stringToAbsorbedDose( trailingstr, gray_definition );
    }catch(...)
    {}
  }//if( string(matches[4]).length() )

  return dose;
}//double stringToAbsorbedDose( std::string str, double gray_definition );


double stringToEquivalentDose( const std::string &str, const double sievert_definition )
{
  /*
   //A poor-persons test for this function:
   std::cout << "\n\nEquivalent Doses:" << std::endl;
   const std::vector<std::string> equivdosestr = { "1.2 sv", "1.23 urem", "100 kilo-rem", "1 gigasv",
     "1gsv", "1.2mrem", "1.2Mrem"
   };
   for( const std::string str : equivdosestr )
   {
     const double val = PhysicalUnits::stringToEquivalentDose( str );
     std::cout << "\t\"" << str << "\" --> " << PhysicalUnits::printToBestEquivalentDoseUnits(val,4,true)
               << " and " << PhysicalUnits::printToBestEquivalentDoseUnits(val,4,false) << std::endl;
   }
   */
  const string regex_str = "\\s*\\+?(\\d+(\\.\\d*)?(?:[Ee][+\\-]?\\d+)?)"
                             "\\s*(" METRIC_PREFIX_UNITS ")*?"
                             "\\s*\\-*\\s*"
                             EQUIVALENT_DOSE_UNIT_REGEX
                             "\\s*(\\d.+)?";
  
  boost::smatch matches;
  boost::regex expression( regex_str, boost::regex::ECMAScript|boost::regex::icase );

  if( !boost::regex_match( str, matches, expression ) )
  {
    char msg[128];
    snprintf(msg, sizeof(msg), "'%s' is not an equivalent dose", str.c_str() );
    //cerr << endl << msg << endl;
    throw std::runtime_error( msg );
  }//if( we dont have a match )

  //cout << endl << endl;
  //for( size_t i = 0; i < matches.size(); ++i )
  //  cerr << "stringToDistance Match " << i << " is " << string( matches[i].first, matches[i].second ) << endl;
  
  string floatstr( matches[1].first, matches[1].second );
  string prefix( matches[3].first, matches[3].second );
  string unitstr( matches[4].first, matches[4].second );
  string trailingstr( matches[5].first, matches[5].second );
  SpecUtils::trim( floatstr );  //These trims may not be necassary, but whatever
  SpecUtils::trim( prefix );
  SpecUtils::trim( unitstr );
  SpecUtils::to_lower_ascii( unitstr );
  
  //cout << "\n\"" << str << "\" : {number->" << floatstr << ", prefix->" << prefix
  //     << ", unit->" << unitstr << ", trailing->" << trailingstr << "}" << endl;
  
  double unitval = 0.0;
  if( unitstr == "sievert" || unitstr == "sv" )
    unitval = sievert;
  else if( unitstr == "rem" || SpecUtils::contains(unitstr,"entgen") )
    unitval = rem;
  else
    throw runtime_error( "Unexpeced equivalent dose unit: '" + unitstr + "'" );

  double dose = unitval * std::stod( floatstr );  //shouldnt ever throw, right?
  dose *= metrix_prefix_value( prefix );
  dose *= (sievert_definition / sievert);

  //if there are characters past what we needed for a match, maybe they are
  //  another distance string, lets try to add them on
  if( trailingstr.length() > 2 )
  {
    try
    {
      dose += stringToEquivalentDose( trailingstr, sievert_definition );
    }catch(...)
    {}
  }//if( string(matches[4]).length() )

  return dose;
}//double stringToEquivalentDose( std::string str, double gray_definition )


double stringToMass( const std::string &str, const double gram_def )
{
  const string regex_str =
  "(" DECIMAL_REGEX ")\\s*"
  "(" METRIC_PREFIX_UNITS ")*?"
  "\\s*\\-*\\s*"
  "(gram|grams|g|ounces|ounce|oz|grain|grains|pounds|pound|lbs|lb|stone|stones)\\s*"
  "(" DECIMAL_REGEX ".+)?";
  
  boost::smatch matches;
  boost::regex expression( regex_str, boost::regex::ECMAScript|boost::regex::icase );

  if( !boost::regex_match( str, matches, expression ) )
  {
    char msg[128];
    snprintf(msg, sizeof(msg), "'%s' is not a mass", str.c_str() );
    //cerr << endl << msg << endl;
    throw std::runtime_error( msg );
  }//if( we dont have a match )
    
  string floatstr( matches[1].first, matches[1].second );
  string prefix( matches[6].first, matches[6].second );
  string unitstr( matches[7].first, matches[7].second );
  string trailingstr( matches[8].first, matches[8].second );
  
  SpecUtils::ireplace_all(floatstr, " ", ""); //The string could be like "+ 2.1", which will cause stod to throw
  SpecUtils::trim( floatstr );  //These trims may not be necessary, but whatever
  SpecUtils::trim( prefix );
  SpecUtils::trim( unitstr );
  SpecUtils::to_lower_ascii( unitstr );
  SpecUtils::trim( trailingstr );

  double unitval = 0.0;
  if( unitstr == "gram" || unitstr == "grams" || unitstr == "g" )
    unitval = gram;
  else if( unitstr == "ounce" || unitstr == "ounces" || unitstr == "oz" )
    unitval = 28.3495*gram;
  else if( unitstr == "pounds" || unitstr == "pound" || unitstr == "lbs" || unitstr == "lb" )
    unitval = 453.592*gram;
  else if( unitstr == "stone" || unitstr == "stones" )
    unitval = 6350.29*gram;
  else if( unitstr == "grain" || unitstr == "grains" )
    unitval = 0.0647989*gram;
  else
    throw runtime_error( "Unexpected mass unit: '" + unitstr + "'" );

  double mass = unitval * std::stod( floatstr );  //shouldnt ever throw, right?
  
  mass *= metrix_prefix_value( prefix );
  mass *= (gram_def / gram);
  
  if( !trailingstr.empty() )
    mass += stringToMass( trailingstr, gram_def );

  return mass;
}//double stringToMass( const std::string &str, const double gram_def )


double stringToEnergy( const std::string &str, const double keV_def )
{
  const string regex_str = "(" DECIMAL_REGEX ")"
                           "\\s*(k|M|killo|mega)?"
                           "\\s*\\-*\\s*"
                           "(eV|electron-volts|electronvolts)"
                           "\\s*(" DECIMAL_REGEX ".+)?";

  boost::smatch matches;
  boost::regex expression( regex_str, boost::regex::ECMAScript | boost::regex::icase );

  if( !boost::regex_match( str, matches, expression ) )
  {
    char msg[128];
    snprintf( msg, sizeof( msg ), "'%s' is not an energy", str.c_str() );
    throw std::runtime_error( msg );
  }//if( we dont have a match )

  //cout << "matches.size() = " << matches.size() << endl;
  //for( size_t i = 0; i < matches.size(); ++i )
  //  cout << "matches[" << i << "]:\t'" << string( matches[i].first, matches[i].second ) << "'" << endl;
  //matches[0] : '511 keV -11 keV'
  //matches[1] : '511 '
  //matches[2] : '511'
  //matches[3] : '511'
  //matches[4] : ''
  //matches[5] : ''
  //matches[6] : 'k'
  //matches[7] : 'eV'
  //matches[8] : '-11 keV'
    
  assert( matches.size() == 13 );

  string floatstr( matches[1].first, matches[1].second );
  string prefix( matches[6].first, matches[6].second );
  SpecUtils::trim( floatstr );
  SpecUtils::trim( prefix );
  
#if( !defined(NDEBUG) )
  string unitstr( matches[7].first, matches[7].second );
  SpecUtils::trim( unitstr );
  SpecUtils::to_lower_ascii( unitstr );
  assert( unitstr == "ev" );
#endif

  double energy = std::stod( floatstr );  //shouldnt ever throw, right?
  energy *= metrix_prefix_value( prefix ) * 0.001;
  energy *= (keV_def / keV);

  string trailingstr( matches[8].first, matches[8].second );
  SpecUtils::trim( trailingstr );
  if( !trailingstr.empty() )
  {
    //remove spaces so like "2614 keV - 511keV" will work (the "- 511keV" doesnt work)
    SpecUtils::ireplace_all( trailingstr, " ", "" ); 
    energy += stringToEnergy( trailingstr, keV_def );
  }

  return energy;
}//stringToEnergy(...)


std::string printValueWithUncertainty( double value, double uncert, size_t unsigfig )
{
  // We will explicitly round value/uncert, 
  //cout << "Value=" << value << ", uncert=" << uncert << ", nsigfig=" << nsigfig << endl;
  const int nsigfig = std::max(1,static_cast<int>(unsigfig));
  
  //Get the exponent
  double valorder = std::floor( std::log10(value) );  //1.2345E-06 will give value -6
  //cout << "\tvalorder=" << valorder << endl;
  double normalizer = std::pow( 10.0, -valorder );
  //cout << "\tnormalizer=" << normalizer << endl;
  double roundedval = value * normalizer;  //value will now be like 1.2345
  //cout << "\troundedval_0=" << roundedval << endl;
  roundedval = std::floor(roundedval * std::pow(10.0,nsigfig-1) + 0.5) / std::pow(10.0,nsigfig-1);
  //cout << "\troundedval_1=" << roundedval << endl;
  roundedval /= normalizer;
  //cout << "\troundedval_2=" << roundedval << endl;
  
  // \TODO: Currently, if value uncertainty and uncerainty "overlap" with more than one digit, then
  //        we will only print uncertainty out to the same decimal order as the value; if they only
  //        "overlap" by one digit, or not at all, then we will print out one more decimal point
  //        than value
  double uncertorder = std::floor( std::log10(uncert) );
  const double numoverlap = nsigfig - (valorder - uncertorder);
  //cout << "\tnumoverlap=" << numoverlap << endl;
  //cout << "\t(valorder - uncertorder)=" << (valorder - uncertorder) << ", (nsigfig - 1)=" << (nsigfig - 1) << endl;
  double uncertnorm = (numoverlap < 2) ? (normalizer * 10) : normalizer;
  
  //cout << "\tuncertnorm=" << uncertnorm << endl;
  double rounduncert = std::fabs( uncert ) * uncertnorm;
  //cout << "\trounduncert_1=" << rounduncert << endl;
  rounduncert = std::floor(rounduncert * std::pow(10.0,nsigfig-1) + 0.5) / std::pow(10.0,nsigfig-1);
  //cout << "\trounduncert_2=" << rounduncert << endl;
  rounduncert /= uncertnorm;
  
  char buffer[64] = {'\0'};
  snprintf( buffer, sizeof(buffer), "%.*g \xC2\xB1 %.*g", nsigfig, roundedval, nsigfig, rounduncert );
  
  // Incase "%.*g" isnt supprted somewhere, could instead do
  //char formatflag[64] = {'\0'}, buffer[64] = {'\0'};
  //snprintf( formatflag, sizeof(formatflag), "%%.%ig \xC2\xB1 %%.%ig", nsigfig, nsigfig );
  //snprintf( buffer, sizeof(buffer), formatflag, roundedval, rounduncert );
  
  //cout << "\tAnswer='" << buffer << "'" << endl << endl;
  //cout << "(" << value << ", " << uncert << ", " << nsigfig << ") -> '" << buffer << "'" << endl;
  
  return buffer;
}//printValueWithUncertainty(...)


std::string printCompact( const double value, const size_t sig_figs )
{
  return SpecUtils::printCompact( value, sig_figs );
}


const UnitNameValuePair &bestActivityUnit( const double activity,
                                           bool useCurries )
{
  assert( sm_activityUnitNameValues.size() == 13 );
  
  UnitNameValuePairV::const_iterator begin, end, iter;
  
  if( useCurries )
  {
    begin = sm_activityUnitNameValues.begin() + 7;
    end = sm_activityUnitNameValues.end();
    
    if( activity == 0.0 )
      return sm_activityUnitNameValues[12];  // return just "ci", no prefix
  }else
  {
    begin = sm_activityUnitNameValues.begin();
    end = begin + 7;
    
    if( activity == 0.0 )
      return sm_activityUnitNameValues[2];  // return just "bq", no prefix
  }//if( useCurries ) / else


  for( iter = begin; iter != end; ++iter )
  {
    if( activity <= 100.0*iter->second )
      return *iter;
  }//for( iter = begin; begin != end; ++begin )

  return *(end-1);
}//const UnitNameValuePair &bestActivityUnit( const double activity )

  
const UnitNameValuePair &bestDoseUnitHtml( const double activity,
                                            bool useRem )
{
  assert( sm_doseRateUnitHtmlNameValues.size() == 6 );
  
  UnitNameValuePairV::const_iterator begin, end, iter;
  
  if( useRem )
  {
    begin = sm_doseRateUnitHtmlNameValues.begin();
    end = begin + 3;
  }else
  {
    begin = sm_doseRateUnitHtmlNameValues.begin() + 3;
    end = sm_doseRateUnitHtmlNameValues.end();
  }//if( useRem ) / else
  
  
  for( iter = begin; iter != end; ++iter )
  {
    if( activity <= 100.0*iter->second )
      return *iter;
  }//for( iter = begin; begin != end; ++begin )
  
  return *(end-1);
}//const UnitNameValuePair &bestDoseUnitHtml(...)
  
  
const UnitNameValuePair &bestActivityUnitHtml( const double activity,
                                                bool useCurries )
{
  const UnitNameValuePair &a = bestActivityUnit( activity, useCurries );
  
  const auto begin = sm_activityUnitNameValues.begin();
  const auto end = sm_activityUnitNameValues.end();
  
  assert( sm_activityUnitHtmlNameValues.size() == sm_activityUnitNameValues.size() );
  assert( sm_activityUnitHtmlNameValues.front().second == sm_activityUnitNameValues.front().second );
  assert( sm_activityUnitHtmlNameValues.back().second == sm_activityUnitNameValues.back().second );
  assert( sm_activityUnitHtmlNameValues[sm_activityUnitHtmlNameValues.size()/2].second
                  == sm_activityUnitNameValues[sm_activityUnitNameValues.size()/2].second );
  
  for( auto iter = begin; iter != end; ++iter )
  {
    if( a.second == iter->second )
      return sm_activityUnitHtmlNameValues[iter-begin];
  }
  
  return sm_activityUnitHtmlNameValues[0];  //shouldnt ever happen
}//const UnitNameValuePair &bestActivityUnitHtml(...)


const UnitNameValuePair &bestLengthUnit( const double length )
{
  if( length == 0.0 )
  {
    assert( sm_lengthUnitNameValues[sm_lengthUnitNameValues.size() - 2].first == "m" );
    return sm_lengthUnitNameValues[sm_lengthUnitNameValues.size() - 2];  // return just "m", no prefix
  }

  for( size_t index = 0; (index + 1) < sm_lengthUnitNameValues.size(); ++index )
  {
    const UnitNameValuePair &next_unit = sm_lengthUnitNameValues[index + 1];
    if( length <= next_unit.second )
      return sm_lengthUnitNameValues[index];
  }//for( iter = begin; begin != end; ++begin )

  return sm_lengthUnitNameValues.back();
}//const UnitNameValuePair &bestLengthUnit( const double length )
  
  
const UnitNameValuePair &bestLengthUnitHtml( const double length )
{
  if( length == 0.0 )
  {
    assert( sm_lengthUnitHtmlNameValues[sm_lengthUnitHtmlNameValues.size() - 2].first == "m" );
    return sm_lengthUnitHtmlNameValues[sm_lengthUnitHtmlNameValues.size() - 2];  // return just "m", no prefix
  }

  for( size_t index = 0; (index + 1) < sm_lengthUnitHtmlNameValues.size(); ++index )
  {
    const UnitNameValuePair &next_unit = sm_lengthUnitHtmlNameValues[index + 1];
    if( length <= next_unit.second )
      return sm_lengthUnitHtmlNameValues[index];
  }//for( iter = begin; begin != end; ++begin )

  return sm_lengthUnitHtmlNameValues.back();
}//const UnitNameValuePair &bestLengthUnitHtml( const double length )
  
  
UnitNameValuePair bestTimeUnit( double t )
{
  t = fabs( t );
  
  if( t <= 100.0*picosecond )
    return UnitNameValuePair("pico-seconds", picosecond );
  if( t <= 100.0*nanosecond )
    return UnitNameValuePair("nano-seconds", nanosecond );
  if( t <= 100.0*microsecond )
    return UnitNameValuePair("micro-seconds", microsecond );
  if( t <= 100.0*millisecond )
    return UnitNameValuePair("milli-seconds", millisecond );
  if( t <= 5.0*hour )
    return UnitNameValuePair("seconds", second );
  if( t <= 1.0*year )
    return UnitNameValuePair("days", day );

  return UnitNameValuePair("years", year );
}//UnitNameValuePair bestTimeUnit( const double time )

}//namespace PhysicalUnits
