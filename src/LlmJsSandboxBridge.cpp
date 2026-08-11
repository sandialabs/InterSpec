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

#if( USE_LLM_INTERFACE )

#include <string>
#include <iostream>
#include <exception>

#include <Wt/WApplication>
#include <Wt/WStringStream>

#include "InterSpec/LlmJsSandboxBridge.h"

using namespace std;
using namespace Wt;

using json = nlohmann::json;


LlmJsSandboxBridge::LlmJsSandboxBridge( WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_resultSignal( this, "jsSandboxResult", false ),
    m_pending(),
    m_nextRequestId( 1 )
{
  // Off-screen, zero-size container that hosts the sandboxed iframe the JS
  // helper attaches to us.  We use off-screen positioning rather than
  // display:none so timer throttling rules (which can apply to fully
  // hidden frames in some browsers) do not bite synchronous JS execution.
  setPositionScheme( PositionScheme::Absolute );
  setOffsets( -9999, Wt::Left | Wt::Top );
  setWidth( 1 );
  setHeight( 1 );

  WApplication *app = WApplication::instance();
  if( app )
    app->require( "InterSpec_resources/LlmJsSandbox.js" );

  m_resultSignal.connect( this, &LlmJsSandboxBridge::handleJsResult );

  // Initialize the JS-side singleton with our widget id so it knows
  // where to attach the iframe and where to emit results to.
  WStringStream init_js;
  init_js << "console.log('Initing InterSpec.LlmJsSandbox'); InterSpec.LlmJsSandbox.init('" << id() << "');";
  doJavaScript( init_js.str() );
}//LlmJsSandboxBridge::LlmJsSandboxBridge( ... )


LlmJsSandboxBridge::~LlmJsSandboxBridge()
{
  // Notify any in-flight callbacks that the bridge is going away so
  // higher-level async tool plumbing does not hang waiting forever.
  for( auto &entry : m_pending )
  {
    try { entry.second( std::string("LLM JS sandbox shut down before request completed") ); }
    catch( std::exception &e ) { cerr << "LlmJsSandboxBridge dtor callback threw: " << e.what() << endl; }
  }
  m_pending.clear();
}


void LlmJsSandboxBridge::run( const json &params, AsyncCallback callback )
{
  if( !callback )
    return;

  // Extract and validate parameters.  We are forgiving here because the
  // LLM produced these and may have minor schema slips; on hard errors
  // we fail the callback synchronously rather than dispatching to JS.
  std::string code;
  if( params.contains("code") && params["code"].is_string() )
    code = params["code"].get<std::string>();

  if( code.empty() )
  {
    callback( std::string("'code' parameter is required and must be a non-empty string") );
    return;
  }

  json data = json(nullptr);
  if( params.contains("data") )
    data = params["data"];

  int timeoutMs = 5000;
  if( params.contains("timeoutMs") && params["timeoutMs"].is_number_integer() )
    timeoutMs = params["timeoutMs"].get<int>();
  if( timeoutMs < 1 )
    timeoutMs = 5000;
  if( timeoutMs > 30000 )
    timeoutMs = 30000;

  const std::string requestId = "r" + std::to_string( ++m_nextRequestId );
  m_pending[requestId] = std::move(callback);

  // Build the JS dispatch.  We rely on nlohmann::json::dump() to produce
  // a JS-safe literal for `code` (a string) and `data` (arbitrary JSON);
  // both are valid JS expressions.
  const std::string codeLiteral = json(code).dump();
  const std::string dataLiteral = data.dump();
  const std::string idLiteral   = json(requestId).dump();

  WStringStream js;
  js << "console.log('Calling InterSpec.LlmJsSandbox.run'); InterSpec.LlmJsSandbox.run("
     << idLiteral << ", "
     << codeLiteral << ", "
     << dataLiteral << ", "
     << timeoutMs << ");";

  doJavaScript( js.str() );
}//void LlmJsSandboxBridge::run( ... )


void LlmJsSandboxBridge::handleJsResult( std::string payload )
{
  json incoming;
  try
  {
    incoming = json::parse( payload );
  }catch( std::exception &e )
  {
    cerr << "LlmJsSandboxBridge: failed to parse JS result payload: " << e.what() << endl;
    return;
  }

  const std::string requestId = incoming.value( "requestId", std::string() );
  if( requestId.empty() )
    return;

  auto it = m_pending.find( requestId );
  if( it == m_pending.end() )
    return;  // Stale (e.g. arrived after our local timeout fired)

  AsyncCallback cb = std::move( it->second );
  m_pending.erase( it );

  // Reshape the JS-side payload into the documented tool result schema.
  json response = json::object();

  const bool timedOut = incoming.value( "timedOut", false );
  response["timedOut"]        = timedOut;
  response["executionTimeMs"] = incoming.value( "executionTimeMs", 0 );
  response["stdout"]          = incoming.value( "stdout", json::array() );
  response["stderr"]          = incoming.value( "stderr", json::array() );

  const bool hasErrorMessage
    = incoming.contains("errorMessage") && !incoming["errorMessage"].is_null();

  if( timedOut || hasErrorMessage )
  {
    response["error"] = hasErrorMessage
      ? incoming["errorMessage"].get<std::string>()
      : std::string("Execution timed out");
    if( incoming.contains("errorName") && !incoming["errorName"].is_null() )
      response["errorName"] = incoming["errorName"];
    if( incoming.contains("errorLine") && !incoming["errorLine"].is_null() )
      response["errorLine"] = incoming["errorLine"];
  }else
  {
    const std::string resultType = incoming.value( "resultType", std::string("undefined") );
    response["resultType"] = resultType;

    // Prefer the JSON-serialized form; fall back to the String(value)
    // representation for non-serializable values like functions/symbols.
    if( incoming.contains("resultJson") && !incoming["resultJson"].is_null() )
    {
      try {
        response["result"] = json::parse( incoming["resultJson"].get<std::string>() );
      } catch( std::exception &e ) {
        response["result"] = incoming["resultJson"];
        response["resultParseWarning"] = std::string("Could not parse JSON result: ") + e.what();
      }
    }else if( incoming.contains("resultString") && !incoming["resultString"].is_null() )
    {
      response["result"] = incoming["resultString"];
    }else
    {
      response["result"] = nullptr;
    }
  }//if( error or timeout ) / else

  try
  {
    cb( response );
  }catch( std::exception &e )
  {
    cerr << "LlmJsSandboxBridge: result callback threw: " << e.what() << endl;
  }
}//void LlmJsSandboxBridge::handleJsResult( ... )

#endif // USE_LLM_INTERFACE
