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

#include "QuickJsSandbox.h"

#include <chrono>
#include <string>
#include <variant>

#include <quickjs.h>

#include "InterSpec/LlmJsSandboxBridge.h"

using namespace std;
using json = nlohmann::json;


namespace
{
  /** Runs `code` as the browser sandbox's iframe does (see the `SANDBOX_HTML` bootstrap in
   InterSpec_resources/LlmJsSandbox.js), and returns its payload as a JSON string.  Differences:
   `data` is passed as JSON text, and as QuickJS's `Error.stack` holds only the call stack,
   `errorMessage` is the error followed by its stack, and `errorLine` comes from the stack.
   */
  const char * const ns_bootstrap = R"JS(
(function( code, dataJson ){
  var __stdout = [], __stderr = [];
  function __fmt( args ){
    try {
      return Array.prototype.map.call( args, function(a){
        if( typeof a === 'string' ) return a;
        try { return JSON.stringify(a); } catch(e){ return String(a); }
      }).join(' ');
    } catch(e){ return ''; }
  }
  globalThis.console = {
    log:   function(){ __stdout.push(__fmt(arguments)); },
    info:  function(){ __stdout.push(__fmt(arguments)); },
    debug: function(){ __stdout.push(__fmt(arguments)); },
    warn:  function(){ __stderr.push(__fmt(arguments)); },
    error: function(){ __stderr.push(__fmt(arguments)); }
  };

  var result; var resultType = 'undefined';
  var errorMessage = null, errorLine = null, errorName = null;
  try {
    var data = JSON.parse( dataJson );
    // An expression gives its value; otherwise the code is a function body (ending in a `return`).
    var fn;
    try {
      fn = new Function( 'data', 'return (' + code + ');' );
    } catch( exprErr ) {
      fn = new Function( 'data', code );
    }
    result = fn( data );
    resultType = (result === null) ? 'null' : (typeof result);
  } catch( e ){
    errorMessage = String(e) + ((e && e.stack) ? ('\n' + String(e.stack)) : '');
    errorName = (e && e.name) ? String(e.name) : null;
    if( e && (typeof e.lineNumber === 'number') )
      errorLine = e.lineNumber;
    else if( e && e.stack ){
      var m = /:(\d+):\d+/.exec( String(e.stack) );
      if( m ) errorLine = parseInt( m[1], 10 );
    }
  }

  var resultJson = null, resultString = null;
  if( !errorMessage ){
    try {
      if( resultType === 'undefined' ){
        resultJson = null;
      } else if( resultType === 'function' || resultType === 'symbol' ){
        resultString = String(result);
      } else {
        resultJson = JSON.stringify( result, function(k, v){
          if( typeof v === 'number' && !isFinite(v) ) return String(v);
          if( typeof v === 'bigint' ) return v.toString();
          return v;
        });
        if( resultJson === undefined ){
          resultJson = null;
          resultString = String(result);
        }
      }
    } catch( e ){
      resultJson = null;
      try { resultString = String(result); } catch( e2 ){ resultString = '[unserializable]'; }
    }
  }

  return JSON.stringify({ resultJson: resultJson, resultString: resultString, resultType: resultType,
                          stdout: __stdout, stderr: __stderr, errorMessage: errorMessage,
                          errorName: errorName, errorLine: errorLine });
})
)JS";


  struct Deadline
  {
    std::chrono::steady_clock::time_point end;
    bool passed = false;
  };

  /** Interrupts the run (with an uncatchable error) once past its deadline. */
  int interrupt_handler( JSRuntime *, void *opaque )
  {
    Deadline * const deadline = static_cast<Deadline *>( opaque );
    if( std::chrono::steady_clock::now() >= deadline->end )
      deadline->passed = true;
    return deadline->passed ? 1 : 0;
  }


  /** The pending exception as text. */
  std::string exception_text( JSContext *ctx )
  {
    const JSValue exception = JS_GetException( ctx );
    std::string text = "unknown error";
    if( const char * const str = JS_ToCString( ctx, exception ) )
    {
      text = str;
      JS_FreeCString( ctx, str );
    }
    JS_FreeValue( ctx, exception );
    return text;
  }


  /** Frees the runtime and context however run_code() returns. */
  struct QuickJs
  {
    JSRuntime *rt = nullptr;
    JSContext *ctx = nullptr;
    ~QuickJs()
    {
      if( ctx )
        JS_FreeContext( ctx );
      if( rt )
        JS_FreeRuntime( rt );
    }
  };


  /** Runs `code`, and returns the sandbox payload (see LlmJsSandboxBridge::toolResult). */
  json run_code( const std::string &code, const std::string &data_json, const int timeout_ms )
  {
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    json payload = json::object();

    Deadline deadline;
    deadline.end = start + std::chrono::milliseconds( timeout_ms );

    QuickJs js;
    js.rt = JS_NewRuntime();
    if( !js.rt )
      throw std::runtime_error( "Could not create the JavaScript runtime." );
    JS_SetMemoryLimit( js.rt, QuickJsSandbox::sm_memory_limit_bytes );
    JS_SetMaxStackSize( js.rt, QuickJsSandbox::sm_max_stack_bytes );
    JS_SetInterruptHandler( js.rt, &interrupt_handler, &deadline );

    js.ctx = JS_NewContext( js.rt );
    if( !js.ctx )
      throw std::runtime_error( "Could not create the JavaScript context." );

    const std::string bootstrap = ns_bootstrap;
    const JSValue run_fn = JS_Eval( js.ctx, bootstrap.c_str(), bootstrap.size(), "<sandbox>", JS_EVAL_TYPE_GLOBAL );
    if( JS_IsException( run_fn ) )
      throw std::runtime_error( "The JavaScript sandbox failed to start: " + exception_text( js.ctx ) );

    JSValue args[2] = { JS_NewStringLen( js.ctx, code.data(), code.size() ),
                        JS_NewStringLen( js.ctx, data_json.data(), data_json.size() ) };
    const JSValue ret = JS_Call( js.ctx, run_fn, JS_UNDEFINED, 2, args );
    JS_FreeValue( js.ctx, args[0] );
    JS_FreeValue( js.ctx, args[1] );
    JS_FreeValue( js.ctx, run_fn );

    if( JS_IsException( ret ) )
    {
      // Errors in the code are caught by the bootstrap; what escapes is uncatchable: the time limit,
      //  or running out of memory or stack.
      const std::string message = exception_text( js.ctx );
      payload["timedOut"] = deadline.passed;
      if( !deadline.passed )
        payload["errorMessage"] = message;
    }else
    {
      std::string text;
      if( const char * const str = JS_ToCString( js.ctx, ret ) )
      {
        text = str;
        JS_FreeCString( js.ctx, str );
      }
      JS_FreeValue( js.ctx, ret );

      try
      {
        payload = json::parse( text );
      }catch( std::exception &e )
      {
        payload = json::object();
        payload["errorMessage"] = std::string( "Could not read the sandbox result: " ) + e.what();
      }
      payload["timedOut"] = false;
    }//if( exception ) / else

    payload["executionTimeMs"] = std::chrono::duration_cast<std::chrono::milliseconds>(
                                               std::chrono::steady_clock::now() - start ).count();
    return payload;
  }//run_code(...)
}//namespace


namespace QuickJsSandbox
{

std::variant<json, std::string> run( const json &params )
{
  // The same parameter handling as LlmJsSandboxBridge::run().
  std::string code;
  if( params.contains("code") && params["code"].is_string() )
    code = params["code"].get<std::string>();
  if( code.empty() )
    return std::string( "'code' parameter is required and must be a non-empty string" );

  const json data = params.contains("data") ? params["data"] : json( nullptr );

  int timeout_ms = 5000;
  if( params.contains("timeoutMs") && params["timeoutMs"].is_number_integer() )
    timeout_ms = params["timeoutMs"].get<int>();
  if( timeout_ms < 1 )
    timeout_ms = 5000;
  if( timeout_ms > 30000 )
    timeout_ms = 30000;

  try
  {
    return LlmJsSandboxBridge::toolResult( run_code( code, data.dump(), timeout_ms ) );
  }catch( std::exception &e )
  {
    return std::string( e.what() );
  }
}//run(...)

}//namespace QuickJsSandbox
