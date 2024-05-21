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
#include <array>
#include <regex>
#include <utility>

#include <Wt/WLocale>
#include <Wt/WApplication>
#include <Wt/WMessageResourceBundle>

#include "SpecUtils/StringAlgo.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace std;

namespace PhysicalUnitsLocalized
{
#define PARTIAL_DECIMAL_REGEX "((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?"
  
  double stringToTimeDuration( std::string str, double second_def )
  {
    return PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( str, -1.0, second_def );
  }
  
  
  double stringToTimeDurationPossibleHalfLife( const string &str, const double hl, double sec_def )
  {
    Wt::WApplication *app = wApp;
    const string locale_name = app ? app->locale().name() : string();
    if( locale_name.empty() || SpecUtils::istarts_with(locale_name, "en") )
      return PhysicalUnits::stringToTimeDurationPossibleHalfLife( str, hl, sec_def );
  
    map<string,vector<string>> to_units_from;
    
    Wt::WMessageResourceBundle &bundle = app->messageResourceBundle();
    auto get_translations = [&]( const string &map_key, const string &val_key, const bool allow_space ){
      string csv_list;
      if( bundle.resolveKey( val_key, csv_list ) )
      {
        if( allow_space )
        {
          SpecUtils::split( to_units_from[map_key], csv_list, ", " );
        }else
        {
          SpecUtils::split( to_units_from[map_key], csv_list, "," );
          for( auto &str : to_units_from[map_key] )
            SpecUtils::trim( str );
        }
      }else
      {
        cerr << "Warning: Failed to resolve string key '" << val_key << "' for locale " << locale_name << endl;
      }
    };//get_translations
    
    get_translations( "seconds", "units-labels-second", true );
    get_translations( "minutes", "units-labels-minute", true );
    get_translations( "hours", "units-labels-hours", true );
    get_translations( "days", "units-labels-days", true );
    get_translations( "years", "units-labels-year", true );
    //"half lives" may have space in it, so require comma delimiter
    get_translations( "halflives", "units-labels-half-lives", false );

    
    // We'll split every pair of number {number,unit}, and do a simple replace.
    //  Its simple, but something.
    
    // "[^\\s\\d\\+-]" means anything but space number , or +-.  We want something like \\w, but
    //  \\w wont work for non-ascii characters forming a word
    std::regex regex( "(?:^|\\s|,|\\+|-)((" PARTIAL_DECIMAL_REGEX ")\\s*([^\\s\\d\\+-]+))" );
    std::sregex_iterator iter( begin(str), end(str), regex );
    const std::sregex_iterator str_end;
    
    //cout << "Input '" << str << "' -->" << endl;
    string localized_str;
    string::const_iterator last_pos = begin(str);
    while( iter != str_end )
    {
      const std::smatch &match = *iter;
      
      localized_str += string(last_pos, match[2].first);
      last_pos = match[7].second;
      
      //for( size_t i = 0; i < match.size(); ++i )
      //  cout << "   Match " << i << ": '" << match[i]<< "'"<< endl;
      
      const std::string number = match[2].str();
      const std::string unit = match[7].str();

      bool found_translation = false;
      string en_unit = unit;
      for( const auto iter : to_units_from )
      {
        if( iter.first == "seconds" )
        {
          // We will allow seconds to be prefixed with u, n, p, m
          for( const string &val : iter.second )
          {
            if( SpecUtils::iends_with(unit, val) )
            {
              string prefix = unit.substr(0, unit.size() - val.size());
              if( !prefix.empty() && (prefix.back() == '-') )
                prefix = prefix.substr(0,prefix.size()-1);
              SpecUtils::to_lower_ascii( prefix ); //TODO: this call is totally locale un-aware
              
              if( (prefix == "p") || (prefix == "pico") 
                 || (prefix == "n") || (prefix == "nano")
                 || (prefix == "u") || (prefix == "μ") || (prefix == "micro")
                 || (prefix == "m") || (prefix == "milli") )
              {
                found_translation = true;
                en_unit = prefix + "s";
                break;
              }
            }//if( SpecUtils::iends_with(val, unit) )
          }//for( const string &val : iter.second )
          
          if( found_translation )
            break;
        }//if( iter.first == "seconds" )
        
        
        for( const string &val : iter.second )
        {
          found_translation = SpecUtils::iequals_ascii( val, unit ); //TODO: this call is totally locale un-aware
          if( found_translation )
          {
            en_unit = iter.first;
            break;
          }
        }//for( const string &val : iter.second )
        
        if( found_translation )
          break;
      }//for( const auto iter : to_units_from )
      
      
      //std::cout << "Number: " << number << ", Unit: " << unit << std::endl;
      
      localized_str += (localized_str.empty() ? "" : " ") + number + en_unit;

      ++iter;
    }//while( iter != str_end )
    
    localized_str += string( last_pos, end(str) );
    
    // cout << "Localized string is: '" << localized_str << "'" << endl << endl;
    
    return 0;
  }//double stringToTimeDurationPossibleHalfLife( std::string str, const double halflife, double second_def )
  
  
  std::string printToBestTimeUnits( double time, int maxNpostDecimal, double sec_def  )
  {
    Wt::WApplication *app = wApp;
    const string locale_name = app ? app->locale().name() : string();
    if( locale_name.empty() || SpecUtils::istarts_with(locale_name, "en") )
      return PhysicalUnits::printToBestTimeUnits( time, maxNpostDecimal, sec_def );
    
    const string en_val = PhysicalUnits::printToBestTimeUnits( time, maxNpostDecimal, sec_def );
    
    
    Wt::WMessageResourceBundle &bundle = app->messageResourceBundle();

    std::regex regex( "(?:^|\\s|,|\\+|-)((" PARTIAL_DECIMAL_REGEX "\\s*)([^\\s\\d\\+-]+))" );
    std::sregex_iterator iter( begin(en_val), end(en_val), regex );
    const std::sregex_iterator str_end;
    
    string localized_str;
    string::const_iterator last_pos = begin(en_val);
    while( iter != str_end )
    {
      const std::smatch &match = *iter;
      
      //for( size_t i = 0; i < match.size(); ++i )
      //  cout << "   Match " << i << ": '" << match[i]<< "'"<< endl;
      
      localized_str += string(last_pos, match[2].first);
      last_pos = match[7].second;
      
      const std::string number = match[2].str();
      localized_str += number;
      
      const std::string unit = match[7].str();
      
      bool found = false;
      string localized_unit;
      if( unit == "y" )
        found = bundle.resolveKey( "units-label-years-short", localized_unit);
      else if( unit == "d" )
        found = bundle.resolveKey( "units-label-days-short", localized_unit);
      else if( unit == "h" )
        found = bundle.resolveKey( "units-label-hours-short", localized_unit);
      else if( unit == "m" )
        found = bundle.resolveKey( "units-label-minutes-short", localized_unit);
      else if( SpecUtils::iends_with(unit, "s") )
      {
        found = bundle.resolveKey( "units-label-seconds-short", localized_unit);
        if( found )
        {
          assert( (unit == "s") || (unit == "ms") || (unit == "ps") || (unit == "us") || (unit == "\xCE\xBCs") );
          localized_unit = unit.substr(0,unit.size()-1) + localized_unit;
        }
      }else
      {
        assert( 0 );
      }
      assert( found );
      localized_str += found ? localized_unit : unit;
      
      ++iter;
    }//while( iter != str_end )
  
    //cout << "Converted time duration '" << en_val << "' to '" << localized_str << "'" << endl;
    
    return localized_str;
  }//std::string printToBestTimeUnits( double time, int maxNpostDecimal, double sec_def  )

}//namespace PhysicalUnits
