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

#include <string>
#include <exception>

#include "nlohmann/json.hpp"

#include "Session.h"

#if defined(__EMSCRIPTEN__)
#include <emscripten/emscripten.h>
#define LIGHT_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define LIGHT_EXPORT
#endif

/** The single entry point from JavaScript (and light_cli).

 Takes `{"method": "...", "params": {...}}`, and returns the JSON result, or `{"error": "..."}`.
 The returned pointer is valid until the next call.
 */
extern "C" LIGHT_EXPORT const char *light_call( const char *request )
{
  static Session session;
  static std::string result;

  try
  {
    const nlohmann::json req = nlohmann::json::parse( request ? request : "" );
    const std::string method = req.at( "method" ).get<std::string>();
    const nlohmann::json params = req.value( "params", nlohmann::json::object() );
    result = session.call( method, params ).dump();
  }catch( std::exception &e )
  {
    result = nlohmann::json( { {"error", e.what()} } ).dump();
  }catch( ... )
  {
    result = "{\"error\":\"Unknown error\"}";
  }

  return result.c_str();
}//light_call(...)
