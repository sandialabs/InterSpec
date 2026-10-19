#include "InterSpec_config.h"
#include "InterSpec/LlmMcpResource.h"

#include <map>
#include <deque>
#include <mutex>
#include <atomic>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <variant>
#include <iostream>
#include <optional>
#include <algorithm>
#include <functional>

#include <Wt/Utils.h>
#include <Wt/WServer.h>
#include <Wt/WRandom.h>
#include <Wt/WResource.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include <Wt/Http/Request.h>
#include <Wt/Http/Response.h>
#include <Wt/Http/ResponseContinuation.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmConversationHistory.h"
#include "InterSpec/LlmToolRegistry.h"

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
  /** Floor for how soon a suspended request is woken by its timer, so a nearly-expired deadline
   can not spin.
   */
  const int sm_min_wake_ms = 100;

  /** Deadline for a synchronous tool call (and agent tools run in the session).  Such a call runs
   on the session's thread and can not be cancelled, so giving up early would only orphan the work
   (and the calls queued behind it); some legitimately take many minutes (e.g., isotopics). */
  const int sm_default_call_timeout_ms = 3600000;

  /** How long a finished assistant turn is kept for `assistant_get_result`. */
  const std::chrono::minutes sm_turn_retention( 60 );

  /** Serializes a message to send; invalid UTF-8 (e.g., from a file's text) is replaced, rather than
   throwing, so the client always gets its reply.
   */
  std::string dump_message( const json &message )
  {
    return message.dump( -1, ' ', false, json::error_handler_t::replace );
  }


  /** Build MCP result content array from a tool result JSON.

   If the result contains an "_image" key (convention from image-producing tools),
   the content array will include both a text block and an MCP image block.
   Otherwise, the result is wrapped as a single text content block.
   */
  json formatMcpToolResultContent( json &result )
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
      textBlock["text"] = dump_message( result );
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
      textBlock["text"] = dump_message( result );
      contentArray.push_back( textBlock );
    }

    return contentArray;
  }//formatMcpToolResultContent(...)


  json rpc_result( const json &id, json tool_result )
  {
    if( !tool_result.is_object() )
      tool_result = json{ {"result", std::move(tool_result)} };
    return json{ {"jsonrpc", "2.0"}, {"id", id}, {"result", { {"content", formatMcpToolResultContent( tool_result )} }} };
  }


  json rpc_error( const json &id, const int code, const std::string &message )
  {
    return json{ {"jsonrpc", "2.0"}, {"id", id}, {"error", { {"code", code}, {"message", message} }} };
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


  /** Lets a suspended request be woken when what it waits for may have changed.  A version number
   avoids lost wake-ups: a listener registered against a version that is already stale is called
   right away.
   */
  class ChangeNotifier
  {
  public:
    uint64_t version() const
    {
      std::lock_guard<std::mutex> lock( m_mutex );
      return m_version;
    }

    void notify()
    {
      vector<std::function<void()>> to_call;
      {
        std::lock_guard<std::mutex> lock( m_mutex );
        ++m_version;
        for( auto &key_fcn : m_listeners )
          to_call.push_back( std::move(key_fcn.second) );
        m_listeners.clear();
      }
      for( const std::function<void()> &fcn : to_call )
        fcn();
    }//notify()

    /** Calls `fcn`, once, on the first notification after `seen_version` (right away, if there
     already was one).  A later registration with the same `key` replaces this one.
     */
    void listen( const void *key, const uint64_t seen_version, std::function<void()> fcn )
    {
      {
        std::lock_guard<std::mutex> lock( m_mutex );
        if( m_version == seen_version )
        {
          m_listeners[key] = std::move( fcn );
          return;
        }
      }
      fcn();
    }//listen(...)

  private:
    mutable std::mutex m_mutex;
    uint64_t m_version = 0;
    std::map<const void *, std::function<void()>> m_listeners;
  };//class ChangeNotifier


  /** The eventual result of a tool call (or assistant turn) run in a session; completed from the
   session's thread.
   */
  struct PendingCall
  {
    ChangeNotifier changes;

    mutable std::mutex mutex;
    bool done = false;
    std::variant<json, std::string> result;
    std::deque<json> progress;

    /** The most recent progress item (kept when `progress` is drained), for polling clients. */
    json last_progress;

    /** For assistant turns: set by assistant_cancel. */
    bool cancel_requested = false;

    /** Called once, when the call finishes; set before the call is dispatched. */
    std::function<void()> on_finish;

    void finish( std::variant<json, std::string> value )
    {
      std::function<void()> finished_fcn;
      {
        std::lock_guard<std::mutex> lock( mutex );
        if( done )
          return;
        done = true;
        result = std::move( value );
        finished_fcn.swap( on_finish );
      }
      if( finished_fcn )
        finished_fcn();
      changes.notify();
    }

    void addProgress( const json &item )
    {
      {
        std::lock_guard<std::mutex> lock( mutex );
        if( done )
          return;
        progress.push_back( item );
        last_progress = item;
        while( progress.size() > 200 )  // no one is draining (a non-streaming client); drop the oldest
          progress.pop_front();
      }
      changes.notify();
    }

    std::vector<json> drainProgress()
    {
      std::lock_guard<std::mutex> lock( mutex );
      std::vector<json> items( begin(progress), end(progress) );
      progress.clear();
      return items;
    }

    bool isDone() const
    {
      std::lock_guard<std::mutex> lock( mutex );
      return done;
    }

    json lastProgress() const
    {
      std::lock_guard<std::mutex> lock( mutex );
      return last_progress;
    }

    /** Only meaningful once `isDone()`. */
    std::variant<json, std::string> getResult() const
    {
      std::lock_guard<std::mutex> lock( mutex );
      return result;
    }

    bool cancelRequested() const
    {
      std::lock_guard<std::mutex> lock( mutex );
      return cancel_requested;
    }

    void requestCancel()
    {
      std::lock_guard<std::mutex> lock( mutex );
      cancel_requested = true;
    }

    void clearCancelRequest()
    {
      std::lock_guard<std::mutex> lock( mutex );
      cancel_requested = false;
    }
  };//struct PendingCall


  /** Posts `fcn` to run in the session (with the session lock held, and wApp set), or calls
   `gone` if the session no longer exists.
   */
  void post_to_session( const std::string &wt_session_id, std::function<void()> fcn,
                        std::function<void()> gone )
  {
    Wt::WServer * const server = Wt::WServer::instance();
    if( !server )
    {
      gone();
      return;
    }

    // A delayed schedule() rather than post(): post() goes through one asio strand shared by the
    //  whole server, and runs the function on it when the session is idle - so a long tool call in
    //  one session would hold up every other session's posted work until it finished.
    server->schedule( std::chrono::milliseconds(1), wt_session_id, [fcn](){
      fcn();
      if( wApp )
        wApp->triggerUpdate();
    }, gone );
  }//post_to_session(...)
}//namespace


namespace
{
  /** True for tools that only read state, so may be called while an assistant turn is running;
   other tools are refused then, as they could change what the assistant is working on.
   */
  bool allowed_while_turn_running( const std::string &tool_name )
  {
    if( SpecUtils::starts_with( tool_name, "get_" ) || SpecUtils::starts_with( tool_name, "list_" ) )
      return true;

    // Nuclear data, and calculations that do not touch the session's analysis.
    static const std::set<std::string> sm_read_only = {
      "sources_with_gammas_near_energy", "source_info", "nuclide_decay_chain", "decay_calculator",
      "source_photons", "search_sources_by_energy", "currie_mda_calc", "calculate_dose",
      "photopeak_detection_efficiency", "available_detector_efficiency_functions",
      "detector_efficiency_function_info", "compton_scatter_peak_check", "escape_peak_check",
      "sum_peak_check", "beta_continuum_check", "activity_fit_one_off",
      // Turn handling itself, and the agent server's session tools.
      "assistant_submit_prompt", "invoke_agent", "assistant_get_result", "assistant_cancel",
      "upload_file", "list_session_files", "close_session"
    };

    return sm_read_only.count( tool_name ) > 0;
  }//allowed_while_turn_running(...)


  /** A short description of a progress item from `LlmToolGui::buildProgressFromTurn`. */
  std::string progress_message( const json &item )
  {
    if( !item.is_object() )
      return "";

    const std::string type = item.value( "type", std::string() );
    std::string names;
    const char * const list_key = item.contains( "results" ) ? "results" : "tools";
    if( item.contains( list_key ) && item[list_key].is_array() )
    {
      for( const json &tool : item[list_key] )
      {
        const std::string name = tool.is_string() ? tool.get<std::string>()
                                 : (tool.is_object() ? tool.value( "name", std::string() ) : std::string());
        if( !name.empty() )
          names += (names.empty() ? "" : ", ") + name;
      }
    }//if( item has tool names )

    if( type == "tool_call" )
      return "Calling " + (names.empty() ? std::string("tools") : names);
    if( type == "tool_result" )
      return "Finished " + (names.empty() ? std::string("tools") : names);
    if( type == "assistant_response" )
      return "Writing the response";
    return type;
  }//progress_message(...)


  /** The reply for a turn still running, with its latest step for clients that poll. */
  json running_turn_status( const PendingCall &pending, const std::string &turn_id, const std::string &message )
  {
    json status{ {"status", "running"}, {"turn_id", turn_id}, {"message", message} };
    const std::string latest = progress_message( pending.lastProgress() );
    if( !latest.empty() )
      status["latest_step"] = latest;
    return status;
  }//running_turn_status(...)
}//namespace


/** The assistant turns started through MCP, so they can be polled and cancelled. */
class LlmMcpResource::TurnStore
{
public:
  struct Turn
  {
    std::string id;

    /** The session the turn runs in: the agent-server handle, or else the Wt session id. */
    std::string session_key;

    std::shared_ptr<PendingCall> call;
    std::chrono::steady_clock::time_point created;
  };//struct Turn

  std::shared_ptr<PendingCall> create( const std::string &turn_id, const std::string &session_key )
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    prune();
    Turn turn;
    turn.id = turn_id;
    turn.session_key = session_key;
    turn.call = std::make_shared<PendingCall>();
    turn.created = std::chrono::steady_clock::now();
    m_turns[turn_id] = turn;
    return turn.call;
  }

  std::optional<Turn> find( const std::string &turn_id ) const
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    const auto pos = m_turns.find( turn_id );
    if( pos == end(m_turns) )
      return std::nullopt;
    return pos->second;
  }

private:
  /** Forgets old finished turns; caller must hold m_mutex. */
  void prune()
  {
    const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    for( auto iter = begin(m_turns); iter != end(m_turns); )
    {
      const bool old = ((now - iter->second.created) > sm_turn_retention);
      if( old && (iter->second.call->isDone() || (m_turns.size() > 1000)) )
        iter = m_turns.erase( iter );
      else
        ++iter;
    }
  }//prune()

  mutable std::mutex m_mutex;
  std::map<std::string, Turn> m_turns;
};//class LlmMcpResource::TurnStore


/** Conversation histories for MCP tool calls in sessions without an LLM assistant (e.g., a desktop
 session whose assistant was never opened), so tools that keep state across calls - peak and
 energy-calibration checkpoints - work there too.  Keyed by Wt session id; entries unused for a day
 are dropped.
 */
class LlmMcpResource::HistoryStore
{
public:
  std::shared_ptr<LlmConversationHistory> get( const std::string &wt_session_id )
  {
    const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lock( m_mutex );
    for( auto iter = begin(m_histories); iter != end(m_histories); )
    {
      if( (now - iter->second.second) > std::chrono::hours(24) )
        iter = m_histories.erase( iter );
      else
        ++iter;
    }

    std::pair<std::shared_ptr<LlmConversationHistory>,std::chrono::steady_clock::time_point> &entry
                                                                    = m_histories[wt_session_id];
    if( !entry.first )
      entry.first = std::make_shared<LlmConversationHistory>();
    entry.second = now;
    return entry.first;
  }//get(...)

private:
  std::mutex m_mutex;
  std::map<std::string,std::pair<std::shared_ptr<LlmConversationHistory>,std::chrono::steady_clock::time_point>> m_histories;
};//class LlmMcpResource::HistoryStore


/** A request that is suspended until its result is ready (passed through
 `ResponseContinuation::setData()`).
 */
struct LlmMcpResource::DeferredCall
{
  json request_id;

  /** If true, the response is an SSE stream: keepalives and progress, then the result. */
  bool sse = false;

  std::chrono::steady_clock::time_point deadline;

  /** The current version of what is waited on (see ChangeNotifier). */
  std::function<uint64_t()> version;

  /** Registers a one-shot wake-up for after `seen_version`. */
  std::function<void( const void *key, uint64_t seen_version, std::function<void()> )> listen;

  /** Returns the complete JSON-RPC message to send, or nullopt if not ready; once the deadline has
   passed (`deadline_passed`) it must return a message.
   */
  std::function<std::optional<json>( bool deadline_passed )> poll;

  /** For SSE: progress items (see `LlmToolGui::buildProgressFromTurn`) to report before the result. */
  std::function<std::vector<json>()> drain_progress;

  /** The client's `_meta.progressToken`, or null.  Progress notifications are only sent when the
   client gave one (as MCP requires); each wake-up then sends at least a heartbeat, so clients that
   time out calls with no activity keep waiting.
   */
  json progress_token;

  /** Number of progress notifications sent; MCP requires `progress` to increase. */
  uint64_t progress_sent = 0;

  /** What is being waited on, e.g. "Assistant turn turn-abc", for progress messages. */
  std::string what;

  /** Incremented by each suspend(); a wake-up timer from an earlier suspend is then stale, and does
   nothing (the continuation is reused, so it would otherwise wake the request early).
   */
  std::atomic<uint64_t> wake_generation{ 0 };
};//struct LlmMcpResource::DeferredCall


// --- LlmMcpResource Implementation ---

LlmMcpResource::LlmMcpResource( const std::shared_ptr<const LlmConfig> &config )
  : Wt::WResource(),
    llm_config_( config ),
    tool_registry_( config ? std::make_shared<LlmTools::ToolRegistry>( *config )
                           : std::shared_ptr<LlmTools::ToolRegistry>() ),
    turns_( std::make_shared<TurnStore>() ),
    histories_( std::make_shared<HistoryStore>() )
#if( MCP_ENABLE_AUTH )
    , secret_token_( config ? config->mcpServer.bearerToken : std::string() )
#endif
{
  if( !llm_config_ )
    throw std::logic_error( "LlmMcpResource must be initialized with a valid config" );

  if( !config->mcpServer.enabled )
    throw std::logic_error( "Cannot initialize LlmMcpResource if MCP isnt enabled in the config" );

#if( BUILD_AS_MCP_HTTP_SERVER )
  if( !McpSessionManager::instance() )
    throw std::logic_error( "The MCP session manager must be started before the MCP resource." );

  server_instructions_ = McpAgentTools::serverInstructions();
#endif

  register_default_tools();
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
    std::shared_ptr<DeferredCall> call;
    try
    {
      call = Wt::cpp17::any_cast<std::shared_ptr<DeferredCall>>( continuation->data() );
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
      response.out() << dump_message( rpc_error( nullptr, -32600, "Invalid Request: missing 'method' field" ) );
      return;
    }

    const std::string method = mcp_request["method"];
    id = mcp_request.value( "id", json( nullptr ) );


    // --- initialize ---
    if( method == "initialize" )
    {
      json result = {
        {"protocolVersion", "2024-11-05"},
        {"capabilities", {
          {"tools", json::object()}
        }},
        {"serverInfo", {
          {"name", "InterSpec MCP Server"},
          {"version", "1.0.0"}
        }}
      };

      if( !server_instructions_.empty() )
        result["instructions"] = server_instructions_;

      response.setMimeType( "application/json" );
      response.out() << json{ {"jsonrpc", "2.0"}, {"id", id}, {"result", result} }.dump();
      return;
    }


    // --- notifications (notifications/initialized, notifications/cancelled, ...) ---
    if( SpecUtils::starts_with( method, "notifications/" ) || !mcp_request.contains( "id" ) )
    {
      // A notification gets no JSON-RPC response; MCP's HTTP transport answers it with 202.
      response.setStatus( 202 );
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
        response.out() << dump_message( rpc_error( id, -32602, "Missing 'params' for tools/call" ) );
        return;
      }

      const json &call_params = mcp_request["params"];
      const std::string tool_name = call_params.at( "name" );
      json arguments = call_params.value( "arguments", json::object() );
      if( !arguments.is_object() )
        arguments = json::object();

      json progress_token;
      if( call_params.contains( "_meta" ) && call_params["_meta"].is_object()
          && call_params["_meta"].contains( "progressToken" ) )
        progress_token = call_params["_meta"]["progressToken"];

      handle_tool_call( id, tool_name, std::move(arguments), progress_token, request, response );
      return;
    }//if( method == "tools/call" )


    // --- Unknown method ---
    response.setMimeType( "application/json" );
    response.out() << dump_message( rpc_error( id, -32601, "Method not found: " + method ) );
  }catch( const json::exception &e )
  {
    // JSON parse or access error
    response.setMimeType( "application/json" );
    response.out() << dump_message( rpc_error( id, -32700, std::string( "Parse error: " ) + e.what() ) );
  }catch( const std::exception &e )
  {
    // Tool execution or other error
    response.setMimeType( "application/json" );
    response.out() << dump_message( rpc_error( id, -32603, e.what() ) );
  }
}//handle_jsonrpc(...)


void LlmMcpResource::handle_tool_call( const json &id, const std::string &tool_name, json arguments,
                                       const json &progress_token,
                                       const Wt::Http::Request &request, Wt::Http::Response &response )
{
  if( registered_tools_.find( tool_name ) == registered_tools_.end() )
  {
    response.setMimeType( "application/json" );
    response.out() << dump_message( rpc_error( id, -32602, "Unknown tool: " + tool_name ) );
    return;
  }

  // Extract and remove userSession from arguments before passing to tool
  std::string user_session;
  if( arguments.contains( "userSession" ) && arguments["userSession"].is_string() )
    user_session = arguments["userSession"].get<std::string>();
  arguments.erase( "userSession" );

  // Set by us for turns; never taken from a client.
  arguments.erase( "_mcp_turn_id" );

#if( BUILD_AS_MCP_HTTP_SERVER )
  const McpAgentTools::ToolDef * const agent_tool = McpAgentTools::toolDefinition( tool_name );
  if( agent_tool && !agent_tool->runs_in_session )
  {
    handle_agent_tool( id, tool_name, user_session, std::move(arguments), request, response );
    return;
  }
#endif

  const std::string wt_session_id = resolve_target_session( user_session ); // throws with a message for the client

#if( BUILD_AS_MCP_HTTP_SERVER )
  const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( manager )
    manager->touch( user_session );

  // The session a turn belongs to, for later assistant_get_result / assistant_cancel calls.
  const std::string &session_key = user_session;
#else
  const std::string &session_key = wt_session_id;
#endif

  if( tool_name == "assistant_get_result" )
  {
    handle_turn_result( id, session_key, arguments, progress_token, request, response );
    return;
  }

  const LlmTools::SharedTool * const tool = tool_registry_->getTool( tool_name );
  bool known_tool = (tool || (tool_name == "assistant_cancel"));
#if( BUILD_AS_MCP_HTTP_SERVER )
  known_tool = (known_tool || agent_tool);
#endif
  if( !known_tool )
    throw std::runtime_error( "Tool not found: " + tool_name );

  const std::string accept = request.headerValue( "Accept" );
  const bool client_accepts_sse = (accept.find( "text/event-stream" ) != std::string::npos);

  const bool is_turn = ((tool_name == "assistant_submit_prompt") || (tool_name == "invoke_agent"));

  // How long a non-streaming client waits for a turn, before being told it is still running.
  int wait_ms = sm_default_turn_wait_s * 1000;
  if( arguments.contains( "wait_seconds" ) && arguments["wait_seconds"].is_number() )
    wait_ms = static_cast<int>( 1000.0 * std::clamp( arguments["wait_seconds"].get<double>(), 0.0,
                                                      static_cast<double>(sm_max_turn_wait_s) ) );
  arguments.erase( "wait_seconds" );

  std::string turn_id;
  std::shared_ptr<PendingCall> pending, cancel_target;

  if( is_turn )
  {
    turn_id = "turn-" + Wt::WRandom::generateId( 16 );
    pending = turns_->create( turn_id, session_key );
    arguments["_mcp_turn_id"] = turn_id;
  }else if( tool_name == "assistant_cancel" )
  {
    const std::string target_turn = arguments.value( "turn_id", std::string() );
    const std::optional<TurnStore::Turn> turn = turns_->find( target_turn );
    if( !turn || (turn->session_key != session_key) )
      throw std::runtime_error( "Unknown turn_id '" + target_turn + "' for this session." );
    // Set now, so a client woken by the cancelled turn finishing sees it was cancelled; cleared
    //  below if there was nothing to cancel.
    turn->call->requestCancel();
    cancel_target = turn->call;
    pending = std::make_shared<PendingCall>();
  }else
  {
    pending = std::make_shared<PendingCall>();
  }//if( is_turn ) / else

#if( BUILD_AS_MCP_HTTP_SERVER )
  // A session with a call in progress is not idle, however long the call takes.
  if( manager )
  {
    manager->callStarted( user_session );
    const std::weak_ptr<McpSessionManager> weak_manager = manager;
    pending->on_finish = [weak_manager, user_session](){
      const std::shared_ptr<McpSessionManager> mgr = weak_manager.lock();
      if( mgr )
        mgr->callFinished( user_session );
    };
  }//if( manager )
#endif

  // Run it in the session.
  const std::shared_ptr<const LlmTools::ToolRegistry> registry = tool_registry_;
  const std::shared_ptr<HistoryStore> histories = histories_;

  post_to_session( wt_session_id,
    [registry, histories, tool_name, arguments, pending, cancel_target, user_session](){
      InterSpecApp * const app = dynamic_cast<InterSpecApp *>( Wt::WApplication::instance() );
      InterSpec * const viewer = app ? app->viewer() : nullptr;
      if( !viewer )
      {
        pending->finish( std::string( "No InterSpec viewer in the session." ) );
        return;
      }

      try
      {
        // A tool that changes the analysis could disrupt (or, by changing the foreground, cancel)
        //  a running assistant turn.
        LlmToolGui * const gui = viewer->currentLlmTool();
        if( gui && gui->isTurnInProgress() && !allowed_while_turn_running( tool_name ) )
          throw std::runtime_error( "The assistant is working on a turn in this session; wait for it to"
                                    " finish (assistant_get_result), or cancel it (assistant_cancel), before"
                                    " calling '" + tool_name + "'." );

#if( BUILD_AS_MCP_HTTP_SERVER )
        if( tool_name == "export_analysis" )
        {
          pending->finish( McpAgentTools::exportAnalysis( viewer, user_session, arguments ) );
        }else
#endif
        if( tool_name == "assistant_cancel" )
        {
          const bool cancelled = gui && gui->cancelMcpPrompt( arguments.value( "turn_id", std::string() ) );
          if( !cancelled && cancel_target )
            cancel_target->clearCancelRequest();
          pending->finish( json{ {"cancelled", cancelled},
                                 {"message", cancelled ? "The turn was cancelled."
                                                       : "The turn was not queued or running (it may have already finished)."} } );
        }else
        {
          const LlmTools::SharedTool * const tool = registry->getTool( tool_name );
          if( !tool )
            throw std::runtime_error( "Tool not found: " + tool_name );

          if( tool->isAsync() )
          {
            const LlmTools::SharedTool::AsyncCallback done = [pending]( std::variant<json, std::string> result ){
              pending->finish( std::move(result) );
            };

            if( tool->isStreaming() )
              tool->asyncStreamingExecutor( arguments, viewer, nullptr, nullptr, done,
                                            [pending]( const json &item ){ pending->addProgress( item ); } );
            else
              tool->asyncExecutor( arguments, viewer, nullptr, nullptr, done );
          }else
          {
            // Tools that keep state across calls (peak and energy-calibration checkpoints) use the
            //  conversation history; share the session assistant's, so both see the same checkpoints.
            LlmInterface * const llm = gui ? gui->llmInterface() : nullptr;
            std::shared_ptr<LlmConversationHistory> history = llm ? llm->getHistory() : nullptr;
            if( !history )
              history = histories->get( Wt::WApplication::instance()->sessionId() );
            pending->finish( registry->executeTool( tool_name, arguments, viewer, nullptr, history.get() ) );
          }
        }//if( tool_name == ... ) / else
      }catch( const std::exception &e )
      {
        pending->finish( std::string( e.what() ) );
      }
    },
    [pending](){
      pending->finish( std::string( "The session is no longer running." ) );
    } );

  // Then wait for it.
  const int tool_timeout_ms = (tool && tool->isAsync()) ? LlmTools::effective_async_timeout_ms( *tool )
                                                         : sm_default_call_timeout_ms;

  const std::shared_ptr<DeferredCall> call = std::make_shared<DeferredCall>();
  call->request_id = id;
  call->sse = client_accepts_sse && (is_turn || (tool && tool->isAsync()));
  // A turn returns its turn_id within wait_seconds even when streaming: clients commonly give up on
  //  a call after a few minutes, and without the turn_id they could neither poll nor cancel it.
  call->deadline = std::chrono::steady_clock::now()
                   + std::chrono::milliseconds( is_turn ? wait_ms : tool_timeout_ms );
  call->progress_token = progress_token;
  call->what = is_turn ? ("Assistant turn " + turn_id) : ("Tool " + tool_name);
  call->version = [pending](){ return pending->changes.version(); };
  call->listen = [pending]( const void *key, uint64_t seen, std::function<void()> fcn ){
    pending->changes.listen( key, seen, std::move(fcn) );
  };
  call->drain_progress = [pending](){ return pending->drainProgress(); };

#if( BUILD_AS_MCP_HTTP_SERVER )
  const std::string base_url = public_base_url( request );
#endif

  call->poll = [=]( const bool deadline_passed ) -> std::optional<json> {
    if( !pending->isDone() )
    {
      if( !deadline_passed )
        return std::nullopt;

      if( is_turn )
        return rpc_result( id, running_turn_status( *pending, turn_id,
          "The assistant is still working; call assistant_get_result with this turn_id to wait for the result." ) );

      return rpc_error( id, -32603, LlmTools::async_timeout_error_message( tool_name, tool_timeout_ms ) );
    }//if( !pending->isDone() )

    std::variant<json, std::string> result_or_error = pending->getResult();

    if( const std::string * const err = std::get_if<std::string>( &result_or_error ) )
    {
      if( is_turn )
        return rpc_result( id, json{ {"status", pending->cancelRequested() ? "cancelled" : "error"},
                                     {"turn_id", turn_id}, {"error", *err} } );
      return rpc_error( id, -32603, *err );
    }

    json result = std::get<json>( std::move(result_or_error) );
    if( !result.is_object() )
      result = json{ {"result", std::move(result)} };

    // Charts for the in-app assistant's display (inline SVG, chart data) - tens of KB that an MCP
    //  client has no use for; the assistant strips them before results reach its LLM, too.
    result.erase( "_visuals" );

    if( is_turn )
    {
      result["status"] = "done";
      result["turn_id"] = turn_id;
    }

#if( BUILD_AS_MCP_HTTP_SERVER )
    {
      const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
      try
      {
        // Images are also saved as files, so a client can hand them on.
        if( manager && result.contains( "_image" ) && result["_image"].is_object() )
        {
          const std::string mime = result["_image"].value( "mimeType", "image/png" );
          const std::string data = Wt::Utils::base64Decode( result["_image"].value( "base64Data", "" ) );
          if( !data.empty() )
          {
            const std::string ext = (mime == "image/svg+xml") ? ".svg" : ((mime == "image/jpeg") ? ".jpg" : ".png");
            const std::string relpath = McpAgentTools::saveOutputFile( user_session, tool_name + ext, data );
            const std::string token = manager->addDownload( user_session, relpath );
            result["image_file"] = json{ {"path", relpath},
              {"download_url", base_url + "/mcp-api/files/" + token + "/" + SpecUtils::filename(relpath)} };
          }
        }//if( an image )

        if( manager && result.contains( "files" ) && result["files"].is_array() )
        {
          for( json &file : result["files"] )
          {
            const std::string relpath = file.value( "path", std::string() );
            if( !relpath.empty() && SpecUtils::starts_with( relpath, "output/" ) )
              file["download_url"] = base_url + "/mcp-api/files/" + manager->addDownload( user_session, relpath )
                                     + "/" + SpecUtils::filename( relpath );
          }
        }//if( exported files )
      }catch( const std::exception &e )
      {
        result["file_error"] = e.what();
      }

      result["userSession"] = user_session;
    }
#else
    result["userSession"] = wt_session_id;
#endif

    return rpc_result( id, std::move(result) );
  };//call->poll

  start_deferred( response, call );
}//handle_tool_call(...)


void LlmMcpResource::handle_turn_result( const json &id, const std::string &session_key, const json &arguments,
                                         const json &progress_token,
                                         const Wt::Http::Request &request, Wt::Http::Response &response )
{
  const std::string turn_id = arguments.value( "turn_id", std::string() );
  const std::optional<TurnStore::Turn> turn = turns_->find( turn_id );
  if( !turn || (turn->session_key != session_key) )
    throw std::runtime_error( "Unknown turn_id '" + turn_id + "' for this session (results are kept for "
                              + std::to_string(sm_turn_retention.count()) + " minutes)." );

  int wait_ms = sm_default_turn_wait_s * 1000;
  if( arguments.contains( "wait_seconds" ) && arguments["wait_seconds"].is_number() )
    wait_ms = static_cast<int>( 1000.0 * std::clamp( arguments["wait_seconds"].get<double>(), 0.0,
                                                      static_cast<double>(sm_max_turn_wait_s) ) );

  const std::shared_ptr<PendingCall> pending = turn->call;
  const std::shared_ptr<DeferredCall> call = std::make_shared<DeferredCall>();
  call->request_id = id;
  call->deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds( wait_ms );
  call->sse = (request.headerValue( "Accept" ).find( "text/event-stream" ) != std::string::npos);
  call->progress_token = progress_token;
  call->what = "Assistant turn " + turn_id;
  call->drain_progress = [pending](){ return pending->drainProgress(); };
  call->version = [pending](){ return pending->changes.version(); };
  call->listen = [pending]( const void *key, uint64_t seen, std::function<void()> fcn ){
    pending->changes.listen( key, seen, std::move(fcn) );
  };
  call->poll = [pending, id, turn_id, session_key]( const bool deadline_passed ) -> std::optional<json> {
    if( !pending->isDone() )
    {
      if( !deadline_passed )
        return std::nullopt;
      return rpc_result( id, running_turn_status( *pending, turn_id,
                                                  "Still working; call assistant_get_result again." ) );
    }

    std::variant<json, std::string> result_or_error = pending->getResult();
    if( const std::string * const err = std::get_if<std::string>( &result_or_error ) )
      return rpc_result( id, json{ {"status", pending->cancelRequested() ? "cancelled" : "error"},
                                   {"turn_id", turn_id}, {"error", *err} } );

    json result = std::get<json>( std::move(result_or_error) );
    if( !result.is_object() )
      result = json{ {"result", std::move(result)} };
    result["status"] = "done";
    result["turn_id"] = turn_id;
    result["userSession"] = session_key;
    return rpc_result( id, std::move(result) );
  };

  start_deferred( response, call );
}//handle_turn_result(...)


#if( BUILD_AS_MCP_HTTP_SERVER )
void LlmMcpResource::handle_agent_tool( const json &id, const std::string &tool_name,
                                        const std::string &user_session, json arguments,
                                        const Wt::Http::Request &request, Wt::Http::Response &response )
{
  const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( !manager )
    throw std::runtime_error( "Not running as an MCP agent server." );

  const std::string base_url = public_base_url( request );

  const McpAgentTools::ToolDef * const def = McpAgentTools::toolDefinition( tool_name );
  if( def && def->needs_session && (tool_name != "close_session") )
  {
    // Checks the session is live, with a helpful message if not.
    resolve_target_session( user_session );
    manager->touch( user_session );
  }

  const auto write_now = [&response]( const json &message ){
    response.setMimeType( "application/json" );
    response.out() << dump_message( message );
  };

  if( tool_name == "create_session" )
  {
    const std::string handle = manager->createSession(); // throws if at the session limit

    const std::shared_ptr<DeferredCall> call = std::make_shared<DeferredCall>();
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

    start_deferred( response, call );
    return;
  }//if( tool_name == "create_session" )

  if( tool_name == "close_session" )
  {
    const std::optional<McpSessionManager::SessionInfo> info = user_session.empty() ? std::nullopt
                                                                                     : manager->session( user_session );
    if( !info )
      throw std::runtime_error( "Unknown userSession." );

    if( info->state == McpSessionManager::State::Ended )
    {
      write_now( rpc_result( id, json{ {"closed", false}, {"userSession", user_session},
                                       {"message", "The session had already ended (" + info->end_reason + ")."} } ) );
      return;
    }

    manager->endSession( user_session, "closed by the client" );
    write_now( rpc_result( id, json{ {"closed", true}, {"userSession", user_session},
                                     {"files_kept", manager->settings().keepSessionFiles} } ) );
    return;
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
    write_now( rpc_result( id, json{ {"path", path}, {"bytes", num_bytes} } ) );
    return;
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
    write_now( rpc_result( id, files ) );
    return;
  }//if( tool_name == "list_session_files" )

  if( tool_name == "get_usage_guide" )
  {
    const std::string guide = McpAgentTools::usageGuide();
    if( guide.empty() )
      throw std::runtime_error( "The usage guide is not available on this server." );
    write_now( rpc_result( id, json{ {"guide", guide} } ) );
    return;
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


std::string LlmMcpResource::resolve_target_session( const std::string &user_session ) const
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


void LlmMcpResource::start_deferred( Wt::Http::Response &response, const std::shared_ptr<DeferredCall> &call )
{
  const uint64_t version = call->version();

  if( const std::optional<json> message = call->poll( false ) )
  {
    response.setMimeType( "application/json" );
    response.out() << dump_message( *message );
    return;
  }

  if( call->sse )
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


void LlmMcpResource::continue_deferred( Wt::Http::Response &response, const std::shared_ptr<DeferredCall> &call )
{
  // Without a continuation Wt completes the response, so on an error the client must still be sent
  //  a JSON-RPC reply, rather than an empty body (or an SSE stream with no result).
  try
  {
    continue_deferred_impl( response, call );
  }catch( const std::exception &e )
  {
    const json message = rpc_error( call->request_id, -32603, std::string("Internal error: ") + e.what() );
    if( call->sse )
      response.out() << "event: message\ndata: " << dump_message( message ) << "\n\n";
    else
      response.out() << dump_message( message );
  }
}//continue_deferred(...)


void LlmMcpResource::continue_deferred_impl( Wt::Http::Response &response,
                                             const std::shared_ptr<DeferredCall> &call )
{
  const uint64_t version = call->version();

  if( call->sse && call->drain_progress )
  {
    std::vector<json> items = call->drain_progress();
    if( !call->progress_token.is_null() )
    {
      if( items.empty() )
        items.push_back( json() );   // a heartbeat

      for( const json &item : items )
      {
        std::string message = progress_message( item );
        message = call->what + ": " + (message.empty() ? std::string("still working") : message);
        const json notification{ {"jsonrpc", "2.0"}, {"method", "notifications/progress"},
          {"params", { {"progressToken", call->progress_token},
                       {"progress", ++call->progress_sent},
                       {"message", message} }} };
        response.out() << "event: message\ndata: " << dump_message( notification ) << "\n\n";
      }
    }//if( the client asked for progress )
  }//if( call->sse && call->drain_progress )

  const bool deadline_passed = (std::chrono::steady_clock::now() >= call->deadline);
  const std::optional<json> message = call->poll( deadline_passed );
  if( message )
  {
    // No new continuation, so Wt completes the response.
    if( call->sse )
      response.out() << "event: message\ndata: " << dump_message( *message ) << "\n\n";
    else
      response.out() << dump_message( *message );
    return;
  }

  assert( !deadline_passed );
  if( call->sse )
    response.out() << ": keepalive\n\n";

  suspend( response, call, version );
}//continue_deferred_impl(...)


void LlmMcpResource::suspend( Wt::Http::Response &response, const std::shared_ptr<DeferredCall> &call,
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
  const std::weak_ptr<DeferredCall> weak_call = call;
  const std::function<void()> timer_wake = [weak_call, generation, wake](){
    const std::shared_ptr<DeferredCall> deferred = weak_call.lock();
    if( deferred && (deferred->wake_generation == generation) )
      wake();
  };

  // And at the deadline, or sooner: in time to send a keepalive (SSE), and in any case at least every
  //  sm_sse_keepalive_ms - Wt's io_service does not finish stopping until every pending timer has
  //  fired, so a long deadline must not become a long timer (a woken call that is not done just
  //  suspends again).
  const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
  long long delay_ms = std::chrono::duration_cast<std::chrono::milliseconds>( call->deadline - now ).count();
  delay_ms = std::min<long long>( delay_ms, sm_sse_keepalive_ms );
  delay_ms = std::max<long long>( delay_ms, sm_min_wake_ms );

  Wt::WServer * const server = Wt::WServer::instance();
  if( server )
    server->ioService().schedule( std::chrono::milliseconds( delay_ms ), timer_wake );
}//suspend(...)


void LlmMcpResource::register_default_tools()
{
  assert( tool_registry_ && llm_config_ );

  if( !tool_registry_ )
    throw std::logic_error( "LlmMcpResource: should not have null tool registry" );
  if( !llm_config_ )
    throw std::logic_error( "LlmMcpResource: should not have null llm config" );

#if( BUILD_AS_MCP_HTTP_SERVER )
  const char * const user_session_desc = "Required: the userSession value returned by create_session.";
#else
  const char * const user_session_desc = "Optional: the user session identifier. If not specified, will use most recent session.";
#endif

  // Every session tool takes the session to act on.
  const auto add_session_param = [user_session_desc]( json &schema ){
    if( !schema.is_object() )
      return;

    if( !schema.contains( "properties" ) )
      schema["properties"] = json::object();

    if( !schema["properties"].contains( "userSession" ) )
      schema["properties"]["userSession"] = { {"type", "string"}, {"description", user_session_desc} };

#if( BUILD_AS_MCP_HTTP_SERVER )
    json &required = schema["required"];
    if( !required.is_array() )
      required = json::array();
    if( std::find( begin(required), end(required), json("userSession") ) == end(required) )
      required.push_back( "userSession" );
#endif
  };//add_session_param

  // Get tools filtered for MCP (excludes invoke_*, KnowledgeBase, agent-only tools)
  for( const auto &[name, sharedTool] : tool_registry_->getToolsForMcp() )
  {
    ToolInfo info;
    info.name = sharedTool->name;
    // Use the MCP-surface rendering of the description (rendered at registry build); fall back to the
    // agent-surface description if unset.  For marker-free descriptions these are identical.
    info.description = sharedTool->mcpDescription.empty() ? sharedTool->description
                                                          : sharedTool->mcpDescription;
    info.parameters_schema = sharedTool->parameters_schema;
    add_session_param( info.parameters_schema );

    registered_tools_[name] = std::move( info );
  }//for( loop over MCP-filtered tools )

  // Polling and cancelling the assistant turns that the prompt and agent tools start, whenever one
  //  of those is offered (see `<McpServer><AllowAssistantPrompt>` / `<AllowAgentInvocation>`).
  if( registered_tools_.count( "assistant_submit_prompt" ) || registered_tools_.count( "invoke_agent" ) )
  {
    const json turn_id_property{ {"type", "string"}, {"description", "The turn_id returned by assistant_submit_prompt (or invoke_agent)."} };

    ToolInfo result_tool;
    result_tool.name = "assistant_get_result";
    result_tool.description =
      "Waits for an assistant turn started by assistant_submit_prompt (or invoke_agent) that"
      " returned status \"running\", and returns its result: status \"done\" with the response, or"
      " \"error\" / \"cancelled\", or still \"running\" if it did not finish within wait_seconds -"
      " in which case call this again.";
    result_tool.parameters_schema = json{
      {"type", "object"},
      {"properties", {
        {"turn_id", turn_id_property},
        {"wait_seconds", { {"type", "number"}, {"description", "How long to wait for the turn to finish (default 30, at most 50)."} }}
      }},
      {"required", json::array({"turn_id"})}
    };
    add_session_param( result_tool.parameters_schema );
    registered_tools_[result_tool.name] = std::move( result_tool );

    ToolInfo cancel_tool;
    cancel_tool.name = "assistant_cancel";
    cancel_tool.description =
      "Cancels an assistant turn: a queued prompt is dropped, and a running turn is stopped."
      "  The turn then reports status \"cancelled\" from assistant_get_result.";
    cancel_tool.parameters_schema = json{
      {"type", "object"},
      {"properties", { {"turn_id", turn_id_property} }},
      {"required", json::array({"turn_id"})}
    };
    add_session_param( cancel_tool.parameters_schema );
    registered_tools_[cancel_tool.name] = std::move( cancel_tool );
  }//if( an assistant turn tool is offered )

#if( BUILD_AS_MCP_HTTP_SERVER )
  for( const McpAgentTools::ToolDef &tool : McpAgentTools::toolDefinitions() )
    registered_tools_[tool.name] = ToolInfo{ tool.name, tool.description, tool.schema };
#endif

#if( PERFORM_DEVELOPER_CHECKS )
  std::cout << "MCP Resource registered " << registered_tools_.size()
            << " tools from shared registry" << std::endl;
#endif
}//register_default_tools()
