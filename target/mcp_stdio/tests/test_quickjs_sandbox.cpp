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

/** Tests of the QuickJS sandbox for the assistant's `run_javascript` tool (QuickJsSandbox), run
 directly: results, console capture, errors, limits, and that no host access is available.
 */

#include <string>
#include <variant>
#include <iostream>

#include "QuickJsSandbox.h"

using namespace std;
using json = nlohmann::json;

namespace
{
  int ns_failures = 0;

  void check( const string &name, const bool ok, const json &detail = json() )
  {
    cout << (ok ? "PASS: " : "FAIL: ") << name;
    if( !ok && !detail.is_null() )
      cout << " - " << detail.dump().substr( 0, 400 );
    cout << endl;
    if( !ok )
      ++ns_failures;
  }

  /** Runs `code` (with `data`, and `timeout_ms` if positive), and returns the tool result. */
  json run( const string &code, const json &data = json(), const int timeout_ms = 0 )
  {
    json params{ {"code", code} };
    if( !data.is_null() )
      params["data"] = data;
    if( timeout_ms > 0 )
      params["timeoutMs"] = timeout_ms;

    const variant<json, string> result = QuickJsSandbox::run( params );
    if( const string *err = get_if<string>( &result ) )
      return json{ {"callError", *err} };
    return get<json>( result );
  }
}//namespace


int main()
{
  json r = run( "1 + 2" );
  check( "expression", (r.value( "result", json() ) == 3) && (r.value( "resultType", "" ) == "number"), r );

  r = run( "var total = 0;\nfor( var i = 0; i < data.values.length; ++i )\n  total += data.values[i];\nreturn { total: total };",
           json{ {"values", {1, 2, 3.5}} } );
  check( "statements ending in return, with data", r.value( "result", json() ) == json{ {"total", 6.5} }, r );

  // (Statements give no value without a `return`, as in the browser sandbox.)
  r = run( "console.log('a', 1, {b: 2}); console.info('c'); console.warn('w'); console.error('e'); return 7;" );
  check( "console capture",
         (r["stdout"] == json::array({ "a 1 {\"b\":2}", "c" })) && (r["stderr"] == json::array({ "w", "e" }))
         && (r.value( "result", json() ) == 7), r );

  r = run( "var x = 1;\nthrow new RangeError('bad value');" );
  check( "thrown error",
         (r.value( "errorName", "" ) == "RangeError")
         && (r.value( "error", "" ).find( "bad value" ) != string::npos) && r.contains( "errorLine" ), r );

  r = run( "var = ;" );
  check( "syntax error", r.value( "errorName", "" ) == "SyntaxError", r );

  r = run( "while( true ){}", json(), 300 );
  check( "infinite loop times out", (r.value( "timedOut", false ) == true) && (r.value( "executionTimeMs", 0 ) < 2000), r );

  r = run( "var a = []; while( true ) a.push( new Array(100000).fill(1.5) );", json(), 20000 );
  check( "memory bomb is stopped", r.contains( "error" ) && (r.value( "timedOut", true ) == false), r );

  r = run( "function f(n){ return f(n + 1) + 1; } f(0)" );
  check( "unbounded recursion is stopped", r.contains( "error" ), r );

  r = run( "[NaN, Infinity, -Infinity, 1]" );
  check( "NaN and Infinity become strings", r.value( "result", json() ) == json::array({ "NaN", "Infinity", "-Infinity", 1 }), r );

  r = run( "12345678901234567890n" );
  check( "bigint", (r.value( "result", json() ) == "12345678901234567890") && (r.value( "resultType", "" ) == "bigint"), r );

  r = run( "(function hello(){ return 1; })" );
  check( "function result", (r.value( "resultType", "" ) == "function")
                             && r.value( "result", json() ).is_string(), r );

  r = run( "undefined" );
  check( "undefined result", (r.value( "resultType", "" ) == "undefined") && r["result"].is_null(), r );

  r = run( "[typeof require, typeof std, typeof os, typeof fetch, typeof XMLHttpRequest, typeof setTimeout, typeof process]" );
  check( "no host access", r.value( "result", json() ) == json::array({ "undefined", "undefined", "undefined", "undefined",
                                                                        "undefined", "undefined", "undefined" }), r );

  r = run( "" );
  check( "missing code is an error", r.contains( "callError" ), r );

  r = run( "globalThis.leaked = 5; 1" );
  r = run( "typeof leaked" );
  check( "each run is fresh", r.value( "result", json() ) == "undefined", r );

  cout << (ns_failures ? (to_string( ns_failures ) + " FAILED") : string( "ALL PASSED" )) << endl;
  return ns_failures ? 1 : 0;
}//main()
