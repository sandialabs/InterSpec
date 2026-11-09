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

#include <cmath>
#include <cfloat>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayActivityCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace std;


namespace
{
  bool isNucKey( const string &key )
  {
    return (key == "iso") || (key == "nuc") || (key == "isotope") || (key == "nuclide");
  }

  bool isActKey( const string &key )
  {
    return (key == "act") || (key == "activity");
  }

  bool isAgeKey( const string &key )
  {
    return (key == "age") || (key == "initialage");
  }

  /** Drops trailing zeros (and then a trailing decimal point) from the number at the start of
   `value`; ex., "1.500 mCi" -> "1.5 mCi".  Numbers with an exponent are left alone.
   */
  string trimTrailingZeros( const string &value )
  {
    const string::size_type space = value.find( ' ' );
    string number = value.substr( 0, space );
    if( (number.find( '.' ) != string::npos) && (number.find_first_of( "eE" ) == string::npos) )
    {
      while( !number.empty() && (number.back() == '0') )
        number.pop_back();
      if( !number.empty() && (number.back() == '.') )
        number.pop_back();
    }

    return number + ((space == string::npos) ? string() : value.substr( space ));
  }//trimTrailingZeros(...)


  bool isTrueValue( const string &value )
  {
    return (value == "1") || SpecUtils::iequals_ascii( value, "true" )
           || SpecUtils::iequals_ascii( value, "yes" );
  }

  /** Parses a time; URLs are written with English units, but a user may have typed in their own. */
  double parseTime( const string &txt, const double halfLife )
  {
    try
    {
      return PhysicalUnits::stringToTimeDurationPossibleHalfLife( txt, halfLife );
    }catch( std::exception & )
    {
    }

    return PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( txt, halfLife );
  }//parseTime(...)


  /** Returns the name of the display unit with the same value as the unit written `unit` (ex.,
   "microCi" for "uCi"), or empty if there is not one.
   */
  string unitNameForValue( const string &unit )
  {
    double value = -1.0;
    try
    {
      value = PhysicalUnits::stringToActivity( "1 " + unit );
    }catch( std::exception & )
    {
      return "";
    }

    for( const PhysicalUnits::UnitNameValuePair &u : PhysicalUnits::sm_activityUnitNameValues )
    {
      if( fabs(u.second - value) < 1.0E-6 * u.second )
        return u.first;
    }

    return "";
  }//unitNameForValue(...)
}//namespace


namespace DecayActivityCalc
{
const char * const sm_arbitraryDisplayUnit = "arbitrary";


ToolState fromAppUrl( const string &path, string query, const bool defaultUseBq,
                      vector<Wt::WString> &errors )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Decay database not available." );

  ToolState state;

  if( path.empty() || SpecUtils::iequals_ascii(path, "chart") )
    state.tab = Tab::Chart;
  else if( SpecUtils::iequals_ascii(path, "chain") )
    state.tab = Tab::Chain;
  else if( SpecUtils::iequals_ascii(path, "calc") )
    state.tab = Tab::Calc;
  else
    errors.push_back( Wt::WString::tr("dad-url-err-path").arg(path) );

  if( !query.empty() && (query[0] == '?') )
    query = query.substr( 1 );

  vector<string> query_args;
  SpecUtils::split( query_args, query, "&" );

  vector<pair<string,string>> field_values;
  for( string comp : query_args )
  {
    SpecUtils::trim( comp );
    if( comp.empty() )
      continue;

    string::size_type pos = comp.find( '=' );
    if( pos == string::npos )
      pos = comp.find( ':' );

    string key = (pos == string::npos) ? comp : comp.substr( 0, pos );
    string value = (pos == string::npos) ? string() : comp.substr( pos + 1 );
    SpecUtils::trim( key );
    SpecUtils::trim( value );
    SpecUtils::to_lower_ascii( key );

    field_values.emplace_back( key, value );
  }//for( string comp : query_args )

  // With a single nuclide the ordering of the arguments should not matter, so move the nuclide to
  //  the front, letting its activity and age come before it.
  const size_t num_nucs = std::count_if( begin(field_values), end(field_values),
                            []( const pair<string,string> &kv ){ return isNucKey(kv.first); } );
  if( num_nucs == 1 )
  {
    const vector<pair<string,string>>::iterator nuc_pos
          = std::find_if( begin(field_values), end(field_values),
                          []( const pair<string,string> &kv ){ return isNucKey(kv.first); } );
    std::rotate( begin(field_values), nuc_pos, nuc_pos + 1 );
  }//if( num_nucs == 1 )

  // First the settings that are not specific to a nuclide
  bool useBq = defaultUseBq;
  bool haveActUnits = false, haveDispUnits = false;
  for( const pair<string,string> &kv : field_values )
  {
    const string &key = kv.first;
    const string &value = kv.second;

    if( key == "actunits" )
    {
      if( SpecUtils::iequals_ascii(value, "becquerel") || SpecUtils::iequals_ascii(value, "bq") )
      {
        useBq = haveActUnits = true;
      }else if( SpecUtils::iequals_ascii(value, "curie") || SpecUtils::iequals_ascii(value, "ci") )
      {
        useBq = false;
        haveActUnits = true;
      }else
      {
        errors.push_back( Wt::WString::tr("dad-url-err-act-units").arg(value) );
      }
    }else if( key == "dispunits" )
    {
      haveDispUnits = true;
      // Unit names differ only by case (mBq vs MBq), so require an exact match, unless the
      //  case-insensitive match is unambiguous.
      state.displayUnit.clear();
      int num_icase_matches = 0;
      string icase_match;
      for( const PhysicalUnits::UnitNameValuePair &unit : PhysicalUnits::sm_activityUnitNameValues )
      {
        if( unit.first == value )
        {
          state.displayUnit = unit.first;
          break;
        }

        if( SpecUtils::iequals_ascii(unit.first, value) )
        {
          num_icase_matches += 1;
          icase_match = unit.first;
        }
      }//for( loop over activity units )

      if( state.displayUnit.empty() && SpecUtils::iequals_ascii(value, sm_arbitraryDisplayUnit) )
        state.displayUnit = sm_arbitraryDisplayUnit;
      else if( state.displayUnit.empty() && (num_icase_matches == 1) )
        state.displayUnit = icase_match;
      else if( state.displayUnit.empty() && (num_icase_matches == 0) )
        state.displayUnit = unitNameForValue( value ); // ex., "uCi" -> "microCi"

      if( state.displayUnit.empty() )
        errors.push_back( Wt::WString::tr("dad-url-err-disp-units").arg(value) );
    }else if( (key == "time") || (key == "timespan") )
    {
      state.timeSpan = value;
    }else if( (key == "grid") || (key == "gridlines") )
    {
      state.gridLines = isTrueValue( value );
    }else if( key == "logy" )
    {
      state.logY = isTrueValue( value );
    }else if( key == "charttype" )
    {
      if( SpecUtils::iequals_ascii(value, "act") || SpecUtils::iequals_ascii(value, "activity") )
        state.chartType = ChartType::Activity;
      else if( SpecUtils::iequals_ascii(value, "gamma") || SpecUtils::iequals_ascii(value, "gammas") )
        state.chartType = ChartType::Gamma;
      else if( SpecUtils::iequals_ascii(value, "beta") || SpecUtils::iequals_ascii(value, "betas") )
        state.chartType = ChartType::Beta;
      else if( SpecUtils::iequals_ascii(value, "alpha") || SpecUtils::iequals_ascii(value, "alphas") )
        state.chartType = ChartType::Alpha;
      else
        errors.push_back( Wt::WString::tr("dad-url-err-chart-type").arg(value) );
    }else if( !isNucKey(key) && !isActKey(key) && !isAgeKey(key) )
    {
      errors.push_back( Wt::WString::tr("dad-url-err-unknown-key").arg(key) );
    }
  }//for( const pair<string,string> &kv : field_values )

  // "actunits" alone picks the chart units, as it used to be the only way to say them.
  if( state.displayUnit.empty() && haveActUnits )
    state.displayUnit = useBq ? "MBq" : "mCi";

  // URLs written before "dispunits" existed always gave "act" in curies, with the users choice of
  //  units only in "actunits" - so for them, "actunits" decides.
  const bool legacyUnits = haveActUnits && !haveDispUnits;

  // Now the nuclides, each followed by its activity and age
  for( size_t i = 0; i < field_values.size(); ++i )
  {
    if( !isNucKey(field_values[i].first) )
      continue;

    const string &nucstr = field_values[i].second;
    const SandiaDecay::Nuclide * const nuc = db->nuclide( nucstr );
    if( !nuc )
    {
      errors.push_back( Wt::WString::tr("dad-url-err-bad-nuc").arg(nucstr) );
      continue;
    }

    if( nuc->isStable() )
    {
      errors.push_back( Wt::WString::tr("dad-url-err-stable-nuc").arg(nucstr) );
      continue;
    }

    string activity_str, age_str;
    bool have_act = false, have_age = false;
    for( size_t j = i + 1; (j < field_values.size()) && !isNucKey(field_values[j].first); ++j )
    {
      if( isActKey(field_values[j].first) )
      {
        have_act = true;
        activity_str = field_values[j].second;
      }else if( isAgeKey(field_values[j].first) )
      {
        have_age = true;
        age_str = field_values[j].second;
      }
    }//for( loop over the arguments following this nuclide )

    SourceNuclide src;
    src.nuclide = nuc;
    src.activity = 1.0 * PhysicalUnits::mCi;
    src.useCurie = !useBq;

    if( have_act )
    {
      double activity = -1.0;
      try
      {
        activity = PhysicalUnits::stringToActivity( activity_str );
      }catch( std::exception & )
      {
      }

      if( (activity > 0.0) && !IsInf(activity) && !IsNan(activity) )
      {
        const std::optional<bool> actIsCurie = activityStringIsCurie( activity_str );
        src.activity = activity;
        src.useCurie = legacyUnits ? !useBq : actIsCurie.value_or( !useBq );

        // Keep the string as written, unless it is in the other units (then it is reformatted).
        if( !actIsCurie.has_value() || (*actIsCurie == src.useCurie) )
          src.activityStr = activity_str;
      }else
      {
        errors.push_back( Wt::WString::tr("dad-url-err-bad-act").arg(activity_str)
                         .arg(nuc->symbol)
                         .arg(PhysicalUnits::printToBestActivityUnits(src.activity, 2, src.useCurie)) );
      }
    }//if( have_act )

    normalizeActivityString( src );

    if( have_age )
    {
      double age = -1.0;
      try
      {
        age = parseTime( age_str, nuc->halfLife );
      }catch( std::exception & )
      {
      }

      if( (age >= 0.0) && !IsInf(age) && !IsNan(age) )
        src.age = age;
      else
        errors.push_back( Wt::WString::tr("dad-url-err-bad-age").arg(age_str).arg(nuc->symbol) );
    }//if( have_age )

    state.nuclides.push_back( src );
  }//for( size_t i = 0; i < field_values.size(); ++i )

  // Check the time span parses, now we know the first nuclides half-life
  if( !state.timeSpan.empty() )
  {
    const double halfLife = state.nuclides.empty() ? -1.0 : state.nuclides.front().nuclide->halfLife;
    double span = -1.0;
    try
    {
      span = parseTime( state.timeSpan, halfLife );
    }catch( std::exception & )
    {
    }

    if( (span <= 0.0) || IsInf(span) || IsNan(span) )
    {
      errors.push_back( Wt::WString::tr("dad-url-err-bad-time-span").arg(state.timeSpan) );
      state.timeSpan.clear();
    }
  }//if( !state.timeSpan.empty() )

  return state;
}//ToolState fromAppUrl(...)


string toAppUrl( const ToolState &state )
{
  string path;
  switch( state.tab )
  {
    case Tab::Chart: path = "chart"; break;
    case Tab::Chain: path = "chain"; break;
    case Tab::Calc:  path = "calc";  break;
  }//switch( state.tab )

  string query = "grid=" + string(state.gridLines ? "1" : "0");
  query += "&logy=" + string(state.logY ? "1" : "0");

  // "actunits" is still written, for versions of InterSpec before "dispunits" existed.
  const bool arbitrary = SpecUtils::iequals_ascii( state.displayUnit, sm_arbitraryDisplayUnit );
  if( !state.displayUnit.empty() && !arbitrary )
    query += "&actunits=" + string(SpecUtils::icontains(state.displayUnit, "bq") ? "bq" : "ci");
  else if( !state.nuclides.empty() )
    query += "&actunits=" + string(state.nuclides.front().useCurie ? "ci" : "bq");

  if( !state.displayUnit.empty() )
    query += "&dispunits=" + (arbitrary ? string(sm_arbitraryDisplayUnit) : state.displayUnit);

  switch( state.chartType )
  {
    case ChartType::Activity: query += "&charttype=activity"; break;
    case ChartType::Gamma:    query += "&charttype=gamma";    break;
    case ChartType::Beta:     query += "&charttype=beta";     break;
    case ChartType::Alpha:    query += "&charttype=alpha";    break;
  }//switch( state.chartType )

  if( !state.timeSpan.empty() )
    query += "&timespan=" + state.timeSpan;

  for( const SourceNuclide &src : state.nuclides )
  {
    if( !src.nuclide )
      continue;

    SourceNuclide normalized = src;
    normalizeActivityString( normalized );

    query += "&nuc=" + src.nuclide->symbol;
    query += "&act=" + SpecUtils::trim_copy( normalized.activityStr );
    if( src.age > FLT_EPSILON )
      query += "&initialage=" + printTimeForUrl( src.age );
  }//for( const SourceNuclide &src : state.nuclides )

  return path + "?" + query;
}//string toAppUrl( const ToolState &state )


std::optional<bool> activityStringIsCurie( const std::string &activity )
{
  // Look only at the units, after the number (whose exponent may have an 'e' in it).
  const string::size_type last_digit = activity.find_last_of( "0123456789" );
  string units = (last_digit == string::npos) ? activity : activity.substr( last_digit + 1 );
  SpecUtils::trim( units );
  SpecUtils::to_lower_ascii( units );

  // Check for becquerels first, so the 'c' in a "micro" prefix is not taken as curies.
  if( SpecUtils::contains(units, "bq") || SpecUtils::contains(units, "bec") )
    return false;

  if( SpecUtils::contains(units, "ci") || SpecUtils::contains(units, "cu")
     || (!units.empty() && (units.back() == 'c')) )
    return true;

  return std::nullopt;
}//activityStringIsCurie(...)


void normalizeActivityString( SourceNuclide &src )
{
  if( !src.activityStr.empty() )
  {
    try
    {
      const double parsed = PhysicalUnits::stringToActivity( src.activityStr );
      if( fabs(parsed - src.activity) <= 0.001 * fabs(src.activity) )
        return;
    }catch( std::exception & )
    {
    }
  }//if( !src.activityStr.empty() )

  // Three decimals keeps the round-trip to within 0.1%; drop trailing zeros so "1.000 mCi" does
  //  not count as four significant figures.
  src.activityStr = trimTrailingZeros( PhysicalUnits::printToBestActivityUnits( src.activity, 3, src.useCurie ) );
}//void normalizeActivityString( SourceNuclide &src )


string printTimeForUrl( const double time )
{
  return trimTrailingZeros( PhysicalUnits::printToBestTimeUnits( time, 6 ) );
}//string printTimeForUrl( const double time )


double findTimeForActivityFrac( const SandiaDecay::NuclideMixture &mixture, const double frac,
                                const double searchStart, const double searchEnd )
{
  if( (frac <= 0.0) || (frac > 1.0) || (searchStart > searchEnd) )
    throw runtime_error( "findTimeForActivityFrac(...): Invalid frac=" + std::to_string(frac)
                         + ", or time range (" + std::to_string(searchStart) + ", "
                         + std::to_string(searchEnd) + ")" );

  const double epsilon = 0.01 * frac;

  const double initialActivity = mixture.totalActivity( 0.0 );
  if( initialActivity <= 0.0 )
    return 0.0;

  double maxHalfLife = 0.0;
  for( const SandiaDecay::NuclideTimeEvolution &evo : mixture.decayedToNuclidesEvolutions() )
  {
    const double halfLife = evo.nuclide->halfLife;
    if( !IsInf(halfLife) && !IsNan(halfLife) )
      maxHalfLife = std::max( maxHalfLife, halfLife );
  }//for( loop over nuclides )

  const bool ranged = ((searchStart >= 0.0) && (searchEnd > 0.0));
  const double dt = ranged ? (searchEnd - searchStart) / 10.0 : maxHalfLife / 2.0;
  const double startT = ranged ? searchStart : 0.0;

  if( dt <= 0.0 )
    return startT;

  int iteration = 0;
  double time = startT, lastTime = startT;

  while( iteration++ < 5000000 )
  {
    const double currentFrac = mixture.totalActivity( time ) / initialActivity;

    if( fabs(currentFrac - frac) < epsilon )
      return time;

    if( currentFrac <= frac )
    {
      // The start of the range is already below the wanted fraction; shouldnt happen, but the
      //  start is the best answer we have.
      if( time == startT )
        return time;

      // Search the step that crossed the wanted fraction
      return findTimeForActivityFrac( mixture, frac, lastTime, time );
    }//if( currentFrac <= frac )

    lastTime = time;
    time += dt;
  }//while( not yet decayed to the wanted fraction )

  cerr << "DecayActivityCalc::findTimeForActivityFrac(...): couldnt find the time; searched "
       << (time - startT) / SandiaDecay::second << " seconds, at a dt="
       << dt / SandiaDecay::second << endl;

  return -1.0;
}//double findTimeForActivityFrac(...)


double particleRate( const SandiaDecay::Nuclide &nuc, const double activity,
                     const SandiaDecay::ProductType type )
{
  double rate = 0.0;
  for( const SandiaDecay::Transition *trans : nuc.decaysToChildren )
  {
    for( const SandiaDecay::RadParticle &particle : trans->products )
    {
      // Pure E0 transitions emit no photon; NuclideMixture::decayParticle() skips them too.
      if( (particle.type != type)
         || ((type == SandiaDecay::GammaParticle) && particle.e0_verified) )
        continue;

      rate += trans->branchRatio * particle.intensity;
    }//for( loop over particles )
  }//for( loop over decays )

  return rate * activity / SandiaDecay::becquerel;
}//double particleRate(...)

}//namespace DecayActivityCalc
