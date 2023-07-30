#ifndef AppUtils_h
#define AppUtils_h
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

#include <map>
#include <string>
#include <vector>
#include <utility>


namespace AppUtils
{
  /** Splits an app url URI up into components.
   
   This is intended for handling app-urls only, and is definetly not a general solution.
   
   if you pass in a URI of "interspec://decay/chain?nuclide=U238&activity=3uCi..."
   then host=="decay", path=="chain", and query=="nuclide=U238&activity=3uCi..."
   
   You should url-decode the uri, before passing it in, because, for example the '?' character
   is not in the allowed QR code ascii set
   */
  void split_uri( std::string uri, std::string &host, std::string &path,
                 std::string &query, std::string &fragment );
  
  
  /** Splits the query string (i.e., everything that comes after the first '?' character), into key-value pairs.
   
   E.g., "p=1235&val=Agoo asd&dud=p&a:9"  will return { {"A","9"}, {"DUD","p"}, {"P","1235"}, {"VAL","Agoo asd"} }
   
   All keys will be uppercased; values will not be altered.
   The '&' character must be used to seperate key-value pairs, and either the '=' or ':' characters can be used
   to seperate key and value.  If a field (e.g., delimited by '&' characters), does not have a '=' or ':', then it will
   be skipped over.
   
   If a key value is empty (i.e., either '=' or ':' is first character in field), or the same key is specified multiple times
   an exception will be thrown (this isnt standard - but specific to how we use URIs in InterSpec).
   */
  std::map<std::string,std::string> query_str_key_values( const std::string &query );
  
  /** Similar to #query_str_key_values, but keeps key-value pairs in original order, allows duplicates, and empty values. */
  // std::vector<std::pair<std::string,std::string>> query_key_values( const std::string &query );
  
  
#if( USE_BATCH_TOOLS || BUILD_AS_LOCAL_SERVER )
  /** Returns the terminal character width */
  unsigned terminal_width();
#endif
  
#ifdef _WIN32
/** Get command line arguments encoded as UTF-8.
    This function just leaks the memory
 
 Note that environment variables are not in UTF-8, we could account for this
 similar to:
 wchar_t *wenvstrings = GetEnvironmentStringsW();
 ...
 */
void getUtf8Args( int &argc, char **&argv );
#endif
}//namespace AppUtils

#endif //AppUtils_h
