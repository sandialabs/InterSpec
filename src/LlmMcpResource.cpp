#include "InterSpec_config.h"
#include "InterSpec/LlmMcpResource.h"

#include <set>
#include <mutex>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <optional>
#include <algorithm>
#include <functional>

#include <Wt/Utils.h>
#include <Wt/WServer.h>
#include <Wt/WResource.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>
#include <Wt/Http/Request.h>
#include <Wt/Http/Response.h>
#include <Wt/Http/ResponseContinuation.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmMcpDispatcher.h"

#if( BUILD_AS_MCP_HTTP_SERVER )
#include "target/mcp_http/McpAgentTools.h"
#include "target/mcp_http/McpSessionManager.h"
#endif


static_assert( USE_LLM_INTERFACE, "This file should not be compiled unless USE_LLM_INTERFACE is enabled" );
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) || defined(ANDROID) )
static_assert( 0, "MCP server can not be enabled for web or mobile builds" );
#endif


using namespace std;
using json = nlohmann::json;


namespace
{
  std::string dump_message( const json &message )
  {
    return LlmMcpDispatcher::dumpMessage( message );
  }

  json rpc_result( const json &id, json tool_result )
  {
    return LlmMcpDispatcher::rpcResult( id, std::move(tool_result) );
  }

  json rpc_error( const json &id, const int code, const std::string &message )
  {
    return LlmMcpDispatcher::rpcError( id, code, message );
  }

  /** Compares in time independent of where the strings differ. */
  bool constant_time_equals( const std::string &given, const std::string &secret )
  {
    unsigned char diff = (given.size() == secret.size()) ? 0 : 1;
    for( size_t i = 0; i < secret.size(); ++i )
    {
      const char g = given.empty() ? '\0' : given[i % given.size()];
      diff |= static_cast<unsigned char>( g ^ secret[i] );
    }
    return (diff == 0);
  }
}//namespace


// --- LlmMcpResource Implementation ---

LlmMcpResource::LlmMcpResource( const std::shared_ptr<const LlmConfig> &config )
  : Wt::WResource(),
    llm_config_( config ),
    dispatcher_{}
#if( MCP_ENABLE_AUTH )
    , secret_token_( config ? config->mcpServer.bearerToken : std::string() )
#endif
{
  if( !llm_config_ )
    throw std::logic_error( "LlmMcpResource must be initialized with a valid config" );

  if( !config->mcpServer.enabled )
    throw std::logic_error( "Cannot initialize LlmMcpResource if MCP isnt enabled in the config" );

  LlmMcpDispatcher::Options options;
  options.session_argument = true;
#if( BUILD_AS_MCP_HTTP_SERVER )
  if( !McpSessionManager::instance() )
    throw std::logic_error( "The MCP session manager must be started before the MCP resource." );

  options.session_argument_required = true;
  options.session_argument_description = "Required: the userSession value returned by create_session.";
  options.instructions = McpAgentTools::serverInstructions();
#else
  options.session_argument_description = "Optional: the user session identifier. If not specified, will use most recent session.";
#endif

  dispatcher_ = std::make_unique<LlmMcpDispatcher>( config, options, &LlmMcpResource::resolve_target_session );

#if( BUILD_AS_MCP_HTTP_SERVER )
  std::vector<LlmMcpDispatcher::ToolInfo> agent_tools;
  for( const McpAgentTools::ToolDef &tool : McpAgentTools::toolDefinitions() )
    agent_tools.push_back( LlmMcpDispatcher::ToolInfo{ tool.name, tool.description, tool.schema } );

  dispatcher_->addTools( agent_tools,
    [this]( const json &id, const std::string &tool_name, const std::string &user_session,
            json &arguments, const LlmMcpDispatcher::RequestContext &context ){
      return handle_agent_tool( id, tool_name, user_session, arguments, context.base_url );
    } );
#endif
}

LlmMcpResource::~LlmMcpResource()
{
  beingDeleted();
}


void LlmMcpResource::handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response )
{
  // Re-entry of a request waiting on a tool, or an assistant turn.
  if( Wt::Http::ResponseContinuation * const continuation = request.continuation() )
  {
    std::shared_ptr<LlmMcpDispatcher::PendingReply> call;
    try
    {
      call = Wt::cpp17::any_cast<std::shared_ptr<LlmMcpDispatcher::PendingReply>>( continuation->data() );
    }catch( const std::exception & )
    {
      assert( 0 );
    }

    if( call )
      continue_deferred( response, call );
    return;
  }//if( request.continuation() )

  // CORS headers required for some (browser-based) MCP clients; an agent server is only for
  //  back-ends, so a browser page should not be able to talk to it.
#if( !BUILD_AS_MCP_HTTP_SERVER )
  response.addHeader( "Access-Control-Allow-Origin", "*" );
  response.addHeader( "Access-Control-Allow-Methods", "GET, POST, OPTIONS" );
  response.addHeader( "Access-Control-Allow-Headers", "Content-Type, Authorization" );
#endif

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

#if( BUILD_AS_MCP_HTTP_SERVER )
    if( path == "/healthz" )
    {
      handle_health( response );
      return;
    }

    if( SpecUtils::starts_with( path, "/files/" ) )
    {
#if( MCP_ENABLE_AUTH )
      if( !validate_request( request, response ) )
        return;
#endif
      handle_file_request( path.substr( 7 ), request, response );
      return;
    }
#endif //BUILD_AS_MCP_HTTP_SERVER

    // All other paths: 404
    response.setStatus( 404 );
    response.setMimeType( "application/json" );
    response.out() << json{{"error", "Not Found"}}.dump();
  }catch( const std::exception &e )
  {
    response.setStatus( 500 );
    response.setMimeType( "application/json" );
    response.out() << dump_message( json{{"error", e.what()}} );
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
  if( !constant_time_equals( token, secret_token_ ) )
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

  json message;
  try
  {
    std::string request_body;
    request.in().seekg( 0, std::ios::end );
    const size_t length = request.in().tellg();
    request.in().seekg( 0, std::ios::beg );
    request_body.resize( length );
    request.in().read( &request_body[0], length );

    message = json::parse( request_body );
  }catch( const std::exception &e )
  {
    response.setMimeType( "application/json" );
    response.out() << dump_message( rpc_error( nullptr, -32700, std::string( "Parse error: " ) + e.what() ) );
    return;
  }

  LlmMcpDispatcher::RequestContext context;
  context.can_stream = (request.headerValue( "Accept" ).find( "text/event-stream" ) != std::string::npos);
#if( BUILD_AS_MCP_HTTP_SERVER )
  context.base_url = public_base_url( request );
#endif

  const LlmMcpDispatcher::Outcome outcome = dispatcher_->handle( message, context );

  if( outcome.pending )
  {
    start_deferred( response, outcome.pending );
  }else if( outcome.reply )
  {
    response.setMimeType( "application/json" );
    response.out() << dump_message( *outcome.reply );
  }else
  {
    // A notification gets no JSON-RPC response; MCP's HTTP transport answers it with 202.
    response.setStatus( 202 );
  }
}//handle_jsonrpc(...)


#if( BUILD_AS_MCP_HTTP_SERVER )
std::optional<LlmMcpDispatcher::Outcome> LlmMcpResource::handle_agent_tool( const json &id,
                                                                        const std::string &tool_name,
                                                                        const std::string &user_session,
                                                                        const json &arguments,
                                                                        const std::string &base_url )
{
  // Tools run in the session (export_analysis) are handled by the dispatcher, with the others.
  const McpAgentTools::ToolDef * const def = McpAgentTools::toolDefinition( tool_name );
  if( !def || def->runs_in_session )
    return std::nullopt;

  const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( !manager )
    throw std::runtime_error( "Not running as an MCP agent server." );

  if( def->needs_session && (tool_name != "close_session") )
  {
    // Checks the session is live, with a helpful message if not.
    resolve_target_session( user_session );
    manager->touch( user_session );
  }

  const auto reply_now = []( json message ) -> std::optional<LlmMcpDispatcher::Outcome> {
    LlmMcpDispatcher::Outcome outcome;
    outcome.reply = std::move( message );
    return outcome;
  };

  if( tool_name == "create_session" )
  {
    const std::string handle = manager->createSession(); // throws if at the session limit

    const std::shared_ptr<LlmMcpDispatcher::PendingReply> call = std::make_shared<LlmMcpDispatcher::PendingReply>();
    call->request_id = id;
    call->deadline = std::chrono::steady_clock::now() + std::chrono::seconds( 45 );
    call->version = [](){ return uint64_t(0); };
    call->listen = [manager, handle]( const void *, uint64_t, std::function<void()> fcn ){
      manager->addSessionListener( handle, std::move(fcn) );
    };
    call->poll = [manager, handle, id, base_url]( const bool deadline_passed ) -> std::optional<json> {
      const std::optional<McpSessionManager::SessionInfo> info = manager->session( handle );
      if( !info || (info->state == McpSessionManager::State::Ended) )
        return rpc_error( id, -32603, "The session could not be started"
                          + (info ? (": " + info->end_reason) : std::string()) + "." );

      if( info->state == McpSessionManager::State::Starting )
      {
        if( !deadline_passed )
          return std::nullopt;
        manager->endSession( handle, "the session did not start in time" );
        return rpc_error( id, -32603, "The session did not start in time; please try again." );
      }

      const McpSessionManager::AgentSettings &settings = manager->settings();
      const std::time_t expires = std::chrono::system_clock::to_time_t( info->expires );
      char expires_str[32] = { '\0' };
      std::strftime( expires_str, sizeof(expires_str), "%Y-%m-%dT%H:%M:%SZ", std::gmtime( &expires ) );

      return rpc_result( id, json{
        {"userSession", handle},
        {"session_dir", info->dir_name},
        {"upload_url", base_url + "/mcp-api/files/" + info->upload_token},
        {"expires_at", expires_str},
        {"limits", {
          {"max_upload_bytes", settings.maxUploadBytes},
          {"max_session_bytes", settings.maxSessionBytes},
          {"max_image_bytes", settings.maxImageBytes},
          {"idle_timeout_seconds", settings.sessionIdleTimeoutSeconds}
        }},
        {"next_steps", "Put spectrum files in input/ (PUT <upload_url>/<filename>, upload_file, or the shared"
                       " session_dir), load them with load_spectrum_file (path \"input/<filename>\"), then ask"
                       " the assistant with assistant_submit_prompt."}
      } );
    };

    LlmMcpDispatcher::Outcome outcome;
    outcome.pending = call;
    return outcome;
  }//if( tool_name == "create_session" )

  if( tool_name == "close_session" )
  {
    const std::optional<McpSessionManager::SessionInfo> info = user_session.empty() ? std::nullopt
                                                                                     : manager->session( user_session );
    if( !info )
      throw std::runtime_error( "Unknown userSession." );

    if( info->state == McpSessionManager::State::Ended )
    {
      return reply_now( rpc_result( id, json{ {"closed", false}, {"userSession", user_session},
                                       {"message", "The session had already ended (" + info->end_reason + ")."} } ) );
    }

    manager->endSession( user_session, "closed by the client" );
    return reply_now( rpc_result( id, json{ {"closed", true}, {"userSession", user_session},
                                     {"files_kept", manager->settings().keepSessionFiles} } ) );
  }

  if( tool_name == "upload_file" )
  {
    const std::string filename = arguments.value( "filename", std::string() );
    const std::string base64 = arguments.value( "base64Data", std::string() );
    if( filename.empty() || base64.empty() )
      throw std::runtime_error( "upload_file needs 'filename' and 'base64Data'." );

    std::istringstream in( Wt::Utils::base64Decode( base64 ) );
    uint64_t num_bytes = 0;
    const std::string path = McpAgentTools::saveInputFile( user_session, filename, in, num_bytes );
    return reply_now( rpc_result( id, json{ {"path", path}, {"bytes", num_bytes} } ) );
  }//if( tool_name == "upload_file" )

  if( tool_name == "list_session_files" )
  {
    json files = McpAgentTools::listSessionFiles( user_session );
    for( json &file : files["output"] )
    {
      const std::string relpath = file.value( "path", std::string() );
      file["download_url"] = base_url + "/mcp-api/files/" + manager->addDownload( user_session, relpath )
                             + "/" + SpecUtils::filename( relpath );
    }
    files["userSession"] = user_session;
    return reply_now( rpc_result( id, files ) );
  }//if( tool_name == "list_session_files" )

  if( tool_name == "get_usage_guide" )
  {
    const std::string guide = McpAgentTools::usageGuide();
    if( guide.empty() )
      throw std::runtime_error( "The usage guide is not available on this server." );
    return reply_now( rpc_result( id, json{ {"guide", guide} } ) );
  }

  throw std::runtime_error( "Unhandled tool: " + tool_name );
}//handle_agent_tool(...)


void LlmMcpResource::handle_file_request( const std::string &path_info, const Wt::Http::Request &request,
                                          Wt::Http::Response &response )
{
  const auto write_error = [&response]( const int status, const std::string &message ){
    response.setStatus( status );
    response.setMimeType( "application/json" );
    response.out() << dump_message( json{ {"error", message} } );
  };

  const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( !manager )
    return write_error( 404, "Not Found" );

  const size_t slash = path_info.find( '/' );
  const std::string token = path_info.substr( 0, slash );
  const std::string filename = (slash == std::string::npos) ? std::string() : path_info.substr( slash + 1 );

  const std::optional<McpSessionManager::FileToken> file_token = manager->lookupFileToken( token );
  if( token.empty() || !file_token )
    return write_error( 404, "Unknown or expired file URL." );

  const std::optional<McpSessionManager::SessionInfo> info = manager->session( file_token->handle );
  if( !info || (info->state == McpSessionManager::State::Ended) )
    return write_error( 404, "The session has ended." );

  manager->touch( file_token->handle );

  try
  {
    if( file_token->is_upload )
    {
      if( request.tooLarge() )
        return write_error( 413, "The upload is larger than the server allows." );

      if( request.method() == "PUT" )
      {
        if( filename.empty() )
          return write_error( 400, "PUT the file to <upload_url>/<filename>." );

        uint64_t num_bytes = 0;
        const std::string path = McpAgentTools::saveInputFile( file_token->handle, filename, request.in(), num_bytes );
        response.setStatus( 201 );
        response.setMimeType( "application/json" );
        response.out() << json{ {"path", path}, {"bytes", num_bytes} }.dump();
        return;
      }//if( PUT )

      if( request.method() == "POST" )
      {
        const Wt::Http::UploadedFileMap &uploads = request.uploadedFiles();
        if( uploads.empty() )
          return write_error( 400, "POST the file(s) as multipart/form-data, or PUT to <upload_url>/<filename>." );

        json files = json::array();
        for( const auto &field_file : uploads )
        {
          const Wt::Http::UploadedFile &upload = field_file.second;
          std::ifstream input( upload.spoolFileName().c_str(), std::ios::in | std::ios::binary );
          const std::string name = upload.clientFileName().empty() ? (filename.empty() ? field_file.first : filename)
                                                                   : upload.clientFileName();
          uint64_t num_bytes = 0;
          const std::string path = McpAgentTools::saveInputFile( file_token->handle, name, input, num_bytes );
          files.push_back( json{ {"path", path}, {"bytes", num_bytes} } );
        }

        response.setStatus( 201 );
        response.setMimeType( "application/json" );
        response.out() << json{ {"files", files} }.dump();
        return;
      }//if( POST )

      return write_error( 405, "Use PUT (or a multipart POST) to upload." );
    }//if( file_token->is_upload )

    // A download
    if( (request.method() != "GET") && (request.method() != "HEAD") )
      return write_error( 405, "Use GET to download." );

    const std::string path = manager->resolveSessionPath( file_token->handle, file_token->relpath );
    if( !SpecUtils::is_file( path ) )
      return write_error( 404, "The file no longer exists." );

    const std::string name = SpecUtils::filename( path );
    std::string ext = SpecUtils::file_extension( name );
    SpecUtils::to_lower_ascii( ext );

    std::string mime = "application/octet-stream";
    if( (ext == ".n42") || (ext == ".xml") )
      mime = "application/xml";
    else if( ext == ".csv" )
      mime = "text/csv";
    else if( ext == ".html" )
      mime = "text/html";
    else if( ext == ".png" )
      mime = "image/png";
    else if( (ext == ".jpg") || (ext == ".jpeg") )
      mime = "image/jpeg";
    else if( ext == ".svg" )
      mime = "image/svg+xml";
    else if( ext == ".json" )
      mime = "application/json";

    response.setMimeType( mime );
    response.addHeader( "Content-Disposition", "attachment; filename=\"" + name + "\"" );
    response.addHeader( "X-Content-Type-Options", "nosniff" );
    response.setContentLength( static_cast<::uint64_t>( SpecUtils::file_size( path ) ) );

    if( request.method() == "HEAD" )
      return;

    std::ifstream input( path.c_str(), std::ios::in | std::ios::binary );
    char buffer[64*1024];
    while( input )
    {
      input.read( buffer, sizeof(buffer) );
      if( input.gcount() > 0 )
        response.out().write( buffer, input.gcount() );
    }
  }catch( const std::exception &e )
  {
    write_error( 400, e.what() );
  }
}//handle_file_request(...)


void LlmMcpResource::handle_health( Wt::Http::Response &response )
{
  const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  response.setMimeType( "application/json" );
  if( !manager )
  {
    response.setStatus( 503 );
    response.out() << json{ {"status", "unavailable"} }.dump();
    return;
  }

  const McpSessionManager::HealthInfo health = manager->health();
  response.out() << json{
    {"status", "ok"},
    {"sessions", health.live_sessions},
    {"starting_sessions", health.starting_sessions},
    {"max_sessions", health.max_sessions},
    {"browser_running", health.browser_running}
  }.dump();
}//handle_health(...)


std::string LlmMcpResource::public_base_url( const Wt::Http::Request &request ) const
{
  const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  std::string base = manager ? manager->settings().publicBaseUrl : std::string();
  if( base.empty() )
  {
    const std::string host = request.headerValue( "Host" );
    base = "http://" + (host.empty() ? (request.serverName() + ":" + request.serverPort()) : host);
  }

  while( !base.empty() && (base.back() == '/') )
    base.pop_back();

  return base;
}//public_base_url(...)
#endif //BUILD_AS_MCP_HTTP_SERVER


std::string LlmMcpResource::resolve_target_session( const std::string &user_session )
{
#if( BUILD_AS_MCP_HTTP_SERVER )
  {
    if( user_session.empty() )
      throw std::runtime_error( "This tool works on a session: call create_session first, and pass the"
                                " userSession value it returns." );

    const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
    const std::optional<McpSessionManager::SessionInfo> info = manager ? manager->session( user_session ) : std::nullopt;
    if( !info )
      throw std::runtime_error( "Unknown userSession; sessions are created with create_session (an old"
                                " session may have expired)." );

    switch( info->state )
    {
      case McpSessionManager::State::Ended:
        throw std::runtime_error( "This session has ended (" + info->end_reason + "); call create_session"
                                  " to start a new one." );

      case McpSessionManager::State::Starting:
        throw std::runtime_error( "This session is still starting; try again in a few seconds." );

      case McpSessionManager::State::Live:
        break;
    }//switch( info->state )

    if( info->wt_session_id.empty() )
      throw std::runtime_error( "This session is still starting; try again in a few seconds." );

    return info->wt_session_id;
  }
#else
  // Find the target InterSpecApp: by session id or external token, or else the most recent one.
  const std::set<InterSpecApp *> instances = InterSpecApp::runningInstances();
  std::string session_id;
  std::chrono::steady_clock::time_point latest_start{};

  for( InterSpecApp *app : instances )
  {
    Wt::WApplication::UpdateLock lock( app );
    if( !lock )
      continue;

    if( user_session.empty() )
    {
      if( session_id.empty() || (app->startTime() > latest_start) )
      {
        session_id = app->sessionId();
        latest_start = app->startTime();
      }
    }else if( (app->sessionId() == user_session) || (app->externalToken() == user_session) )
    {
      return app->sessionId();
    }
  }//for( InterSpecApp *app : instances )

  if( session_id.empty() )
    throw std::runtime_error( "No active InterSpec session found"
      + (user_session.empty() ? std::string() : (" for session: " + user_session)) );

  return session_id;
#endif //BUILD_AS_MCP_HTTP_SERVER / else
}//resolve_target_session(...)


void LlmMcpResource::start_deferred( Wt::Http::Response &response, const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call )
{
  const uint64_t version = call->version();

  if( const std::optional<json> message = call->poll( false ) )
  {
    response.setMimeType( "application/json" );
    response.out() << dump_message( *message );
    return;
  }

  if( call->streaming )
  {
    response.setMimeType( "text/event-stream" );
    response.addHeader( "Cache-Control", "no-cache" );
    response.out() << ": keepalive\n\n";
  }else
  {
    response.setMimeType( "application/json" );
  }

  suspend( response, call, version );
}//start_deferred(...)


void LlmMcpResource::continue_deferred( Wt::Http::Response &response, const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call )
{
  // Without a continuation Wt completes the response, so on an error the client must still be sent
  //  a JSON-RPC reply, rather than an empty body (or an SSE stream with no result).
  try
  {
    continue_deferred_impl( response, call );
  }catch( const std::exception &e )
  {
    const json message = rpc_error( call->request_id, -32603, std::string("Internal error: ") + e.what() );
    if( call->streaming )
      response.out() << "event: message\ndata: " << dump_message( message ) << "\n\n";
    else
      response.out() << dump_message( message );
  }
}//continue_deferred(...)


void LlmMcpResource::continue_deferred_impl( Wt::Http::Response &response,
                                             const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call )
{
  const uint64_t version = call->version();

  if( call->streaming )
  {
    for( const json &notification : call->takeProgressNotifications() )
      response.out() << "event: message\ndata: " << dump_message( notification ) << "\n\n";
  }

  const bool deadline_passed = (std::chrono::steady_clock::now() >= call->deadline);
  const std::optional<json> message = call->poll( deadline_passed );
  if( message )
  {
    // No new continuation, so Wt completes the response.
    if( call->streaming )
      response.out() << "event: message\ndata: " << dump_message( *message ) << "\n\n";
    else
      response.out() << dump_message( *message );
    return;
  }

  assert( !deadline_passed );
  if( call->streaming )
    response.out() << ": keepalive\n\n";

  suspend( response, call, version );
}//continue_deferred_impl(...)


void LlmMcpResource::suspend( Wt::Http::Response &response, const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call,
                              const uint64_t seen_version )
{
  Wt::Http::ResponseContinuation * const continuation = response.createContinuation();
  continuation->setData( call );
  continuation->waitForMoreData();

  // Wake only this request; the continuation is held weakly, as it goes away with the request.
  const std::weak_ptr<Wt::Http::ResponseContinuation> weak_continuation = continuation->shared_from_this();
  const std::function<void()> wake = [weak_continuation](){
    const std::shared_ptr<Wt::Http::ResponseContinuation> cont = weak_continuation.lock();
    if( cont )
      cont->haveMoreData();
  };

  // (A later registration replaces this listener, so only the timer below can go stale.)
  call->listen( call.get(), seen_version, wake );

  const uint64_t generation = ++call->wake_generation;
  const std::weak_ptr<LlmMcpDispatcher::PendingReply> weak_call = call;
  const std::function<void()> timer_wake = [weak_call, generation, wake](){
    const std::shared_ptr<LlmMcpDispatcher::PendingReply> deferred = weak_call.lock();
    if( deferred && (deferred->wake_generation == generation) )
      wake();
  };

  // And at the deadline, or sooner: in time to send a keepalive (SSE), and in any case at least every
  //  LlmMcpDispatcher::sm_keepalive_ms - Wt's io_service does not finish stopping until every pending timer has
  //  fired, so a long deadline must not become a long timer (a woken call that is not done just
  //  suspends again).
  const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  long long delay_ms = std::chrono::duration_cast<std::chrono::milliseconds>( call->deadline - now ).count();
  delay_ms = std::min<long long>( delay_ms, LlmMcpDispatcher::sm_keepalive_ms );
  delay_ms = std::max<long long>( delay_ms, LlmMcpDispatcher::sm_min_wake_ms );

  Wt::WServer * const server = Wt::WServer::instance();
  if( server )
    server->ioService().schedule( std::chrono::milliseconds( delay_ms ), timer_wake );
}//suspend(...)
