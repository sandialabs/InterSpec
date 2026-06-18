#include "InterSpec_config.h"
#include "InterSpec/LlmMcpResource.h"

#include <map>
#include <future>
#include <string>
#include <variant>
#include <iostream>

#include <Wt/WServer>
#include <Wt/WResource>
#include <Wt/WIOService>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/Http/ResponseContinuation>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmToolRegistry.h"


static_assert( USE_LLM_INTERFACE, "This file should not be compiled unless USE_LLM_INTERFACE is enabled" );
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) || defined(ANDROID) )
static_assert( 0, "MCP server can not be enabled for web or mobile builds" );
#endif


using json = nlohmann::json;


/** Build MCP result content array from a tool result JSON.

 If the result contains an "_image" key (convention from image-producing tools),
 the content array will include both a text block and an MCP image block.
 Otherwise, the result is wrapped as a single text content block.
 */
static json formatMcpToolResultContent( json &result )
{
  json contentArray = json::array();

  if( result.contains( "_image" ) )
  {
    const json &img = result["_image"];
    const std::string base64Data = img.value( "base64Data", "" );
    const std::string mimeType = img.value( "mimeType", "" );

    // Remove _image before dumping the text portion
    result.erase( "_image" );

    json textBlock;
    textBlock["type"] = "text";
    textBlock["text"] = result.dump();
    contentArray.push_back( textBlock );

    if( !base64Data.empty() )
    {
      json imageBlock;
      imageBlock["type"] = "image";
      imageBlock["data"] = base64Data;
      imageBlock["mimeType"] = mimeType;
      contentArray.push_back( imageBlock );
    }
  }
  else
  {
    json textBlock;
    textBlock["type"] = "text";
    textBlock["text"] = result.dump();
    contentArray.push_back( textBlock );
  }

  return contentArray;
}//formatMcpToolResultContent(...)


// --- LlmMcpResource Implementation ---

LlmMcpResource::LlmMcpResource( const std::shared_ptr<const LlmConfig> &config )
  : Wt::WResource(),
    llm_config_( config ),
    tool_registry_( config ? std::make_unique<LlmTools::ToolRegistry>( *config )
                           : std::unique_ptr<LlmTools::ToolRegistry>() )
#if( MCP_ENABLE_AUTH )
    , secret_token_( config ? config->mcpServer.bearerToken : std::string() )
#endif
{
  if( !llm_config_ )
    throw std::logic_error( "LlmMcpResource must be initialized with a valid config" );

  if( !config->mcpServer.enabled )
    throw std::logic_error( "Cannot initialize LlmMcpResource if MCP isnt enabled in the config" );

  m_alive = std::make_shared<std::atomic<bool>>( true );

  register_default_tools();
}


LlmMcpResource::~LlmMcpResource()
{
  m_alive->store( false );
  beingDeleted();
}


void LlmMcpResource::handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response )
{
  // Handle SSE continuation re-entries for async tool calls before anything else.
  // CORS headers were already set on the initial request that started this SSE stream.
  if( request.continuation() )
  {
    handle_sse_continuation( request, response );
    return;
  }

  // CORS headers required for some MCP clients
  response.addHeader( "Access-Control-Allow-Origin", "*" );
  response.addHeader( "Access-Control-Allow-Methods", "GET, POST, OPTIONS" );
  response.addHeader( "Access-Control-Allow-Headers", "Content-Type, Authorization" );

  // Handle CORS preflight
  if( request.method() == "OPTIONS" )
  {
    response.setStatus( 204 );
    return;
  }

  try
  {
    const std::string &path = request.pathInfo();

    // Handle root path — this is what MCP clients (e.g., LM Studio) send to
    if( path.empty() || (path == "/") )
    {
      if( request.method() == "GET" )
      {
        const std::string accept = request.headerValue( "Accept" );
        if( accept.find( "text/event-stream" ) != std::string::npos )
        {
          // SSE connection request
          response.setMimeType( "text/event-stream" );
          response.addHeader( "Cache-Control", "no-cache" );
          response.addHeader( "Connection", "keep-alive" );
          response.out() << "data: {\"jsonrpc\":\"2.0\",\"method\":\"notifications/initialized\",\"params\":{}}\n\n";
          response.out().flush();
          return;
        }

        // Regular GET — return auth discovery info
        handle_get_auth( request, response );
        return;
      }

      if( request.method() == "POST" )
      {
#if( MCP_ENABLE_AUTH )
        if( !validate_request( request, response ) )
          return;
#endif
        handle_jsonrpc( request, response );
        return;
      }
    }//if( root path )

    // GET on /get_auth for MCP discovery
    if( path == "/get_auth" )
    {
      handle_get_auth( request, response );
      return;
    }

    // All other paths: 404
    response.setStatus( 404 );
    response.setMimeType( "application/json" );
    response.out() << json{{"error", "Not Found"}, {"path", path}}.dump();
  }catch( const std::exception &e )
  {
    response.setStatus( 500 );
    response.setMimeType( "application/json" );
    response.out() << json{{"error", e.what()}}.dump();
  }
}//void handleRequest(...)


#if( MCP_ENABLE_AUTH )
bool LlmMcpResource::validate_request( const Wt::Http::Request &request, Wt::Http::Response &response )
{
  if( secret_token_.empty() )
    return true;

  const std::string auth_header = request.headerValue( "Authorization" );
  if( auth_header.empty() )
  {
    response.setStatus( 401 );
    response.setMimeType( "application/json" );
    response.out() << json{{"error", "Missing Authorization header"}}.dump();
    return false;
  }

  const std::string prefix = "Bearer ";
  if( auth_header.rfind( prefix, 0 ) != 0 )
  {
    response.setStatus( 401 );
    response.setMimeType( "application/json" );
    response.out() << json{{"error", "Invalid Authorization header format"}}.dump();
    return false;
  }

  const std::string token = auth_header.substr( prefix.length() );
  if( token != secret_token_ )
  {
    response.setStatus( 401 );
    response.setMimeType( "application/json" );
    response.out() << json{{"error", "Invalid token"}}.dump();
    return false;
  }

  return true;
}//validate_request(...)
#endif // MCP_ENABLE_AUTH


void LlmMcpResource::handle_get_auth( const Wt::Http::Request &request, Wt::Http::Response &response )
{
  response.setMimeType( "application/json" );
#if( MCP_ENABLE_AUTH )
  if( secret_token_.empty() )
    response.out() << json{{"type", "none"}}.dump();
  else
    response.out() << json{{"type", "bearer"}}.dump();
#else
  response.out() << json{{"type", "none"}}.dump();
#endif
}//handle_get_auth(...)


void LlmMcpResource::handle_jsonrpc( const Wt::Http::Request &request, Wt::Http::Response &response )
{
  if( request.method() != "POST" )
  {
    response.setStatus( 405 );
    response.setMimeType( "application/json" );
    response.out() << json{{"error", "Method must be POST"}}.dump();
    return;
  }

  json id = nullptr;  // Track request id for error responses

  try
  {
    std::string request_body;
    request.in().seekg( 0, std::ios::end );
    const size_t length = request.in().tellg();
    request.in().seekg( 0, std::ios::beg );
    request_body.resize( length );
    request.in().read( &request_body[0], length );

    const json mcp_request = json::parse( request_body );

    if( !mcp_request.contains( "method" ) || !mcp_request["method"].is_string() )
    {
      response.setMimeType( "application/json" );
      response.out() << json{
        {"jsonrpc", "2.0"},
        {"id", nullptr},
        {"error", {
          {"code", -32600},
          {"message", "Invalid Request: missing 'method' field"}
        }}
      }.dump();
      return;
    }

    const std::string method = mcp_request["method"];
    id = mcp_request.value( "id", json( nullptr ) );


    // --- initialize ---
    if( method == "initialize" )
    {
      response.setMimeType( "application/json" );
      response.out() << json{
        {"jsonrpc", "2.0"},
        {"id", id},
        {"result", {
          {"protocolVersion", "2024-11-05"},
          {"capabilities", {
            {"tools", json::object()}
          }},
          {"serverInfo", {
            {"name", "InterSpec MCP Server"},
            {"version", "1.0.0"}
          }}
        }}
      }.dump();
      return;
    }


    // --- notifications/initialized ---
    if( method == "notifications/initialized" )
    {
      // Notification — no response required by JSON-RPC spec
      response.setStatus( 204 );
      return;
    }


    // --- tools/list ---
    if( method == "tools/list" )
    {
      json tool_definitions = json::array();
      for( const auto &[name, info] : registered_tools_ )
      {
        tool_definitions.push_back( {
          {"name", info.name},
          {"description", info.description},
          {"inputSchema", info.parameters_schema}
        } );
      }

      response.setMimeType( "application/json" );
      response.out() << json{
        {"jsonrpc", "2.0"},
        {"id", id},
        {"result", {{"tools", tool_definitions}}}
      }.dump();
      return;
    }


    // --- tools/call ---
    if( method == "tools/call" )
    {
      if( !mcp_request.contains( "params" ) )
      {
        response.setMimeType( "application/json" );
        response.out() << json{
          {"jsonrpc", "2.0"},
          {"id", id},
          {"error", {
            {"code", -32602},
            {"message", "Missing 'params' for tools/call"}
          }}
        }.dump();
        return;
      }

      const json &call_params = mcp_request["params"];
      const std::string tool_name = call_params.at( "name" );
      json arguments = call_params.value( "arguments", json::object() );

      // Check tool exists in our registered set
      if( registered_tools_.find( tool_name ) == registered_tools_.end() )
      {
        response.setMimeType( "application/json" );
        response.out() << json{
          {"jsonrpc", "2.0"},
          {"id", id},
          {"error", {
            {"code", -32602},
            {"message", "Unknown tool: " + tool_name}
          }}
        }.dump();
        return;
      }

      // Extract and remove userSession from arguments before passing to tool
      std::string user_session_id;
      if( arguments.contains( "userSession" ) && arguments["userSession"].is_string() )
        user_session_id = arguments["userSession"].get<std::string>();
      arguments.erase( "userSession" );

      // Check if this is an async tool and client supports SSE streaming.
      // If so, use the SSE continuation path to send keepalive messages while
      // the tool runs, preventing client-side timeouts.
      const LlmTools::SharedTool * const tool = tool_registry_->getTool( tool_name );
      const std::string accept = request.headerValue( "Accept" );
      const bool client_accepts_sse = (accept.find( "text/event-stream" ) != std::string::npos);

      if( tool && tool->isAsync() && client_accepts_sse )
      {
        // --- SSE async path: stream keepalives, then final result ---

        // Find the target InterSpecApp (same logic as executeWithSession)
        const std::set<InterSpecApp *> instances = InterSpecApp::runningInstances();
        InterSpecApp *targetApp = nullptr;

        if( user_session_id.empty() )
        {
          std::chrono::steady_clock::time_point earliest_start{};
          for( InterSpecApp *app : instances )
          {
            Wt::WApplication::UpdateLock lock( app );
            if( lock && (!targetApp || (earliest_start > app->startTime())) )
            {
              targetApp = app;
              earliest_start = app->startTime();
            }
          }
        }else
        {
          for( InterSpecApp *app : instances )
          {
            Wt::WApplication::UpdateLock lock( app );
            if( lock && ((app->sessionId() == user_session_id)
                          || (app->externalToken() == user_session_id)) )
            {
              targetApp = app;
              break;
            }
          }
        }

        if( !targetApp )
        {
          throw std::runtime_error( "No active InterSpec session found"
            + (user_session_id.empty() ? std::string()
               : (" for session: " + user_session_id)) );
        }

        // Set up shared state for the continuation chain
        std::shared_ptr<AsyncSseState> state = std::make_shared<AsyncSseState>();
        state->request_id = id;
        state->prom = std::make_shared<std::promise<std::variant<json, std::string>>>();
        state->fut = std::make_shared<std::future<std::variant<json, std::string>>>(
          state->prom->get_future()
        );

        {
          Wt::WApplication::UpdateLock lock( targetApp );
          if( !lock )
            throw std::runtime_error( "Failed to acquire session lock for async tool" );

          InterSpec *viewer = targetApp->viewer();
          if( !viewer )
            throw std::runtime_error( "No InterSpec viewer in session" );

          state->actual_session_id = targetApp->sessionId();

          // Start async operation. The callback captures the shared state and
          // wakes the SSE continuation when the tool finishes.
          LlmMcpResource *self = this;
          tool->asyncExecutor( arguments, viewer, nullptr, nullptr,
            [state, self]( std::variant<json, std::string> result )
            {
              state->prom->set_value( std::move( result ) );
              state->tool_complete.store( true );

              // Wake any waiting SSE continuations so they can send the result.
              // This is safe because the callback runs on the GUI thread via
              // WServer::post(), and the resource is alive during normal operation.
              self->haveMoreData();
            }
          );

          // Flush queued JS (e.g., from doJavaScript) to the browser so the
          // async operation can actually begin on the client side.
          targetApp->triggerUpdate();
        } // UpdateLock released

        // Begin SSE response stream
        response.setMimeType( "text/event-stream" );
        response.addHeader( "Cache-Control", "no-cache" );
        response.out() << ": keepalive\n\n";

        Wt::Http::ResponseContinuation *cont = response.createContinuation();
        cont->setData( state );
        cont->waitForMoreData();

        schedule_sse_keepalive();

        return;
      }//if( async tool with SSE-capable client )

      // --- Standard (blocking) path for sync tools or non-SSE clients ---
      std::string actual_session_id;
      const json result_content = executeWithSession( user_session_id, tool_name,
                                                      arguments, actual_session_id );

      // Inject userSession into result
      json final_result = result_content;
      if( final_result.is_object() )
      {
        final_result["userSession"] = actual_session_id;
      }else
      {
        json wrapped = json::object();
        wrapped["result"] = final_result;
        wrapped["userSession"] = actual_session_id;
        final_result = std::move( wrapped );
      }

      response.setMimeType( "application/json" );
      response.out() << json{
        {"jsonrpc", "2.0"},
        {"id", id},
        {"result", {
          {"content", formatMcpToolResultContent( final_result )}
        }}
      }.dump();
      return;
    }


    // --- Unknown method ---
    response.setMimeType( "application/json" );
    response.out() << json{
      {"jsonrpc", "2.0"},
      {"id", id},
      {"error", {
        {"code", -32601},
        {"message", "Method not found: " + method}
      }}
    }.dump();

  }catch( const json::exception &e )
  {
    // JSON parse or access error
    response.setMimeType( "application/json" );
    response.out() << json{
      {"jsonrpc", "2.0"},
      {"id", id},
      {"error", {
        {"code", -32700},
        {"message", std::string( "Parse error: " ) + e.what()}
      }}
    }.dump();
  }catch( const std::exception &e )
  {
    // Tool execution or other error
    response.setMimeType( "application/json" );
    response.out() << json{
      {"jsonrpc", "2.0"},
      {"id", id},
      {"error", {
        {"code", -32603},
        {"message", e.what()}
      }}
    }.dump();
  }
}//handle_jsonrpc(...)


void LlmMcpResource::register_default_tools()
{
  assert( tool_registry_ && llm_config_ );

  if( !tool_registry_ )
    throw std::logic_error( "LlmMcpResource: should not have null tool registry" );
  if( !llm_config_ )
    throw std::logic_error( "LlmMcpResource: should not have null llm config" );

  // Get tools filtered for MCP (excludes invoke_*, deep research, agent-only tools)
  for( const auto &[name, sharedTool] : tool_registry_->getToolsForMcp() )
  {
    ToolInfo info;
    info.name = sharedTool.name;
    info.description = sharedTool.description;
    info.parameters_schema = sharedTool.parameters_schema;

    // Add userSession parameter to schema
    if( info.parameters_schema.is_object() )
    {
      if( !info.parameters_schema.contains( "properties" ) )
        info.parameters_schema["properties"] = json::object();

      if( !info.parameters_schema["properties"].contains( "userSession" ) )
      {
        info.parameters_schema["properties"]["userSession"] = {
          {"type", "string"},
          {"description", "Optional: the user session identifier. If not specified, will use most recent session."}
        };
      }
    }

    registered_tools_[name] = std::move( info );
  }//for( loop over MCP-filtered tools )

#if( PERFORM_DEVELOPER_CHECKS )
  std::cout << "MCP Resource registered " << registered_tools_.size()
            << " tools from shared registry" << std::endl;
#endif
}//register_default_tools()


json LlmMcpResource::executeWithSession( const std::string &userSessionId,
                                          const std::string &toolName,
                                          const json &params,
                                          std::string &actualSessionId )
{
  const LlmTools::SharedTool * const tool = tool_registry_->getTool( toolName );
  if( !tool )
    throw std::runtime_error( "Tool not found: " + toolName );

  // Find the target InterSpecApp
  const std::set<InterSpecApp *> instances = InterSpecApp::runningInstances();
  InterSpecApp *targetApp = nullptr;

  if( userSessionId.empty() )
  {
    // No session specified — find the most recently started instance
    std::chrono::steady_clock::time_point earliest_start{};
    for( InterSpecApp *app : instances )
    {
      Wt::WApplication::UpdateLock lock( app );
      if( lock && (!targetApp || (earliest_start > app->startTime())) )
      {
        targetApp = app;
        earliest_start = app->startTime();
      }
    }
  }else
  {
    // Match by session ID or external token
    for( InterSpecApp *app : instances )
    {
      Wt::WApplication::UpdateLock lock( app );
      if( lock && ((app->sessionId() == userSessionId)
                    || (app->externalToken() == userSessionId)) )
      {
        targetApp = app;
        break;
      }
    }
  }

  if( !targetApp )
  {
    throw std::runtime_error( "No active InterSpec session found"
      + (userSessionId.empty() ? std::string() : (" for session: " + userSessionId)) );
  }


  if( tool->isAsync() )
  {
    // Async path: acquire lock for Stage A (data capture), then release before blocking.
    // The async callback is delivered via WServer::post() which re-acquires the lock,
    // so we must NOT hold it while waiting — otherwise deadlock.
    std::promise<std::variant<json, std::string>> prom;
    std::future<std::variant<json, std::string>> fut = prom.get_future();

    {
      Wt::WApplication::UpdateLock lock( targetApp );
      if( !lock )
        throw std::runtime_error( "Failed to acquire session lock for async tool" );

      InterSpec *viewer = targetApp->viewer();
      if( !viewer )
        throw std::runtime_error( "No InterSpec viewer in session" );

      actualSessionId = targetApp->sessionId();

      // Start async operation — Stage A runs here with lock held (wApp is set).
      // The callback will be called on the GUI thread (via WServer::post) after
      // Stages B (background computation) and C (result post-processing) complete.
      tool->asyncExecutor( params, viewer, nullptr, nullptr,
        [&prom]( std::variant<json, std::string> result )
        {
          prom.set_value( std::move( result ) );
        }
      );

      // Flush queued JS (e.g., from doJavaScript) to the browser so the
      // async operation can actually begin on the client side.
      targetApp->triggerUpdate();
    } // UpdateLock released — Stages B and C can now proceed

    // Block until the async operation completes
    std::variant<json, std::string> result_or_error = fut.get();

    if( const std::string *err = std::get_if<std::string>( &result_or_error ) )
      throw std::runtime_error( *err );

    return std::get<json>( std::move( result_or_error ) );
  }
  else
  {
    // Sync path: hold lock throughout execution
    Wt::WApplication::UpdateLock lock( targetApp );
    if( !lock )
      throw std::runtime_error( "Failed to acquire session lock" );

    InterSpec *viewer = targetApp->viewer();
    if( !viewer )
      throw std::runtime_error( "No InterSpec viewer in session" );

    actualSessionId = targetApp->sessionId();

    const json result = tool_registry_->executeTool( toolName, params, viewer );

    targetApp->triggerUpdate();

    return result;
  }
}//executeWithSession(...)


void LlmMcpResource::handle_sse_continuation( const Wt::Http::Request &request,
                                               Wt::Http::Response &response )
{
  Wt::Http::ResponseContinuation *continuation = request.continuation();
  assert( continuation );
  if( !continuation )
    return;

  std::shared_ptr<AsyncSseState> state;
  try
  {
    state = boost::any_cast<std::shared_ptr<AsyncSseState>>( continuation->data() );
  }catch( const boost::bad_any_cast & )
  {
    assert( 0 );
    return;
  }

  assert( state );
  if( !state )
    return;

  if( !state->tool_complete.load() )
  {
    // Tool still running: send a keepalive comment and re-suspend.
    response.out() << ": keepalive\n\n";

    Wt::Http::ResponseContinuation *cont = response.createContinuation();
    cont->setData( state );
    cont->waitForMoreData();

    schedule_sse_keepalive();
    return;
  }

  // Tool is complete: send the final JSON-RPC response as an SSE event.
  std::variant<json, std::string> result_or_error = state->fut->get();

  json rpc_response;

  if( const std::string *err = std::get_if<std::string>( &result_or_error ) )
  {
    rpc_response = json{
      {"jsonrpc", "2.0"},
      {"id", state->request_id},
      {"error", {
        {"code", -32603},
        {"message", *err}
      }}
    };
  }else
  {
    json result_content = std::get<json>( std::move( result_or_error ) );

    // Inject userSession into result (same logic as the blocking path)
    json final_result = result_content;
    if( final_result.is_object() )
    {
      final_result["userSession"] = state->actual_session_id;
    }else
    {
      json wrapped = json::object();
      wrapped["result"] = final_result;
      wrapped["userSession"] = state->actual_session_id;
      final_result = std::move( wrapped );
    }

    rpc_response = json{
      {"jsonrpc", "2.0"},
      {"id", state->request_id},
      {"result", {
        {"content", formatMcpToolResultContent( final_result )}
      }}
    };
  }

  // Write as SSE event; no continuation created means Wt sends ResponseDone.
  response.out() << "event: message\ndata: " << rpc_response.dump() << "\n\n";
}//handle_sse_continuation(...)


void LlmMcpResource::schedule_sse_keepalive()
{
  Wt::WServer *server = Wt::WServer::instance();
  if( !server )
    return;

  std::weak_ptr<std::atomic<bool>> weak_alive = m_alive;

  server->ioService().schedule( 10000, [this, weak_alive](){
    const std::shared_ptr<std::atomic<bool>> alive = weak_alive.lock();
    if( alive && alive->load() )
      this->haveMoreData();
  } );
}//schedule_sse_keepalive()
