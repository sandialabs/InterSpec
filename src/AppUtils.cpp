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

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"

using namespace std;

namespace AppUtils
{
  void split_uri( string url, string &host, string &path, string &query, string &fragment )
  {
    // But also see rfc3986 for a regex that might work, or be the start of making one that
    //  will work, to parse the URL.
    
    host = path = query = fragment = "";
    
    //Get rid of (optional) leading "interspec://", so URL will look like: 'drf/specify?v=1&n=MyName&"My other Par"'
    const string interspec_scheme = "interspec://";
    const string spectrum_scheme = "RADDATA://";
    if( SpecUtils::istarts_with(url, interspec_scheme) )
      url = url.substr(interspec_scheme.size());
    else if( SpecUtils::istarts_with(url, spectrum_scheme) )
      url = url.substr(spectrum_scheme.size());
    
    
    const string::size_type q_pos = url.find( '?' );
    //string::size_type q_end_pos = q_pos + 1;
    
    const string::size_type hash_pos = url.find( '#' );
    const string::size_type path_end = std::min( hash_pos, q_pos );
    
    /*
    if( q_pos == string::npos )
    {
      // Since '?' isnt a QR alphanumeric code also allow the URL encoded value of it, "%3F"
     q_pos = url.find( "%3F" );
      if( q_pos == string::npos )
     q_pos = url.find( "%3f" );
      
      if( q_pos != string::npos )
        q_end_pos = q_pos + 3;
    }//
     */
      
    //if( q_pos == string::npos )
    //  throw runtime_error( "App URL did not contain the host/path component that specifies the intent"
    //                       " of the URL (e.g., what tool to use, or what info is contained in the URL)." );
    
    // the q_pos+1 should be safe, even if q_pos is last character in string.
    //if( url.find(q_end_pos, '?') != string::npos )
    //  throw runtime_error( "App URL contained more than one '?' character, which isnt allowed" );

    const string host_path = url.substr( 0, path_end );
    const string::size_type host_end = host_path.find( '/' );
    
    // If "host/path" is only one compnent (i.e., no slash), then populate host,
    //  and have path be empty
    host = (host_end == string::npos) ? host_path : host_path.substr(0,host_end);
    path = (host_end == string::npos) ? string("") : host_path.substr(host_end+1);
    
    if( (q_pos != string::npos) && (q_pos < hash_pos) )
      query = url.substr( q_pos+1, (hash_pos==string::npos) ? string::npos : (hash_pos - q_pos) );
      
    if( hash_pos != string::npos )
      fragment = query.substr( hash_pos + 1 );
  }//void split_uri( string url, string &host, string &path, string &query )

  
  std::map<std::string,std::string> query_str_key_values( const std::string &query_str )
  {
    map<string,string> parts;
    vector<string> components;
    SpecUtils::split( components, query_str, "&" );
    
    for( string comp : components )
    {
      SpecUtils::trim( comp );
      if( comp.empty() )
        continue;
      
      auto pos = comp.find('=');
      if( pos == string::npos )
        pos = comp.find(':');
      
      if( pos == string::npos )
        continue;
      
      string key = comp.substr(0,pos);
      string value = comp.substr(pos+1);
      SpecUtils::trim(key);
      SpecUtils::to_upper_ascii(key);
      
      if( key.empty() )
        throw runtime_error( "fromAppUrl: query portion '" + comp + "' has empty name" );
      
      if( parts.count(key) )
        throw runtime_error( "fromAppUrl: query portion contains duplicate key '" + key + "'" );
      
      parts[key] = value;
    }//for( string comp : components )
    
    return parts;
  }//std::map<std::string,std::string> split_query_str( const std::string &query )
  
  /*
  vector<pair<string,string>> query_key_values( const string &query_str )
  {
    vector<pair<string,string>> parts;
    vector<string> components;
    SpecUtils::split( components, query_str, "&" );
    
    for( string comp : components )
    {
      SpecUtils::trim( comp );
      if( comp.empty() )
        continue;
      
      auto pos = comp.find('=');
      if( pos == string::npos )
        pos = comp.find(':');
      
      string key = comp.substr(0,pos);
      string value = (pos == string::npos) ? string("") : comp.substr(pos+1);
      SpecUtils::trim(key);
      SpecUtils::to_upper_ascii(key);
      
      parts.emplace_back( key, value );
    }//for( string comp : components )
    
    return parts;
  }//vector<pair<string,string>> query_key_values( const string &query );
   */
}//namespace AppUtils
