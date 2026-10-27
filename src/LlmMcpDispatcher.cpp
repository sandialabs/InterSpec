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

#include <set>
#include <map>
#include <deque>
#include <mutex>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <variant>
#include <iostream>
#include <optional>
#include <algorithm>
#include <functional>

#include <boost/asio/steady_timer.hpp>

#include <Wt/Utils.h>
#include <Wt/WServer.h>
#include <Wt/WRandom.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmMcpDispatcher.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/LlmConversationHistory.h"

#if( BUILD_AS_MCP_HTTP_SERVER )
#include "target/mcp_http/McpAgentTools.h"
#include "target/mcp_http/McpSessionManager.h"
#endif

static_assert( USE_MCP_SERVER, "This file should not be compiled unless USE_MCP_SERVER is enabled" );


using namespace std;
using json = nlohmann::json;


namespace
{
  /** Deadline for a synchronous tool call (and agent tools run in the session).  Such a call runs
   on the session's thread and can not be cancelled, so giving up early would only orphan the work
   (and the calls queued behind it); some legitimately take many minutes (e.g., isotopics). */
  const int sm_default_call_timeout_ms = 3600000;

  /** How long a finished assistant turn is kept for `assistant_get_result`. */
  const std::chrono::minutes sm_turn_retention( 60 );


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
      textBlock["text"] = LlmMcpDispatcher::dumpMessage( result );
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
      textBlock["text"] = LlmMcpDispatcher::dumpMessage( result );
      contentArray.push_back( textBlock );
    }

    return contentArray;
  }//formatMcpToolResultContent(...)


  /** Lets a waiting reply be woken when what it waits for may have changed.  A version number
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


  /** Wires a `PendingReply` to wait on `pending`. */
  void wait_on( LlmMcpDispatcher::PendingReply &reply, const std::shared_ptr<PendingCall> &pending )
  {
    reply.version = [pending](){ return pending->changes.version(); };
    reply.listen = [pending]( const void *key, uint64_t seen, std::function<void()> fcn ){
      pending->changes.listen( key, seen, std::move(fcn) );
    };
    reply.drain_progress = [pending](){ return pending->drainProgress(); };
  }//wait_on(...)


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


  /** The `wait_seconds` argument of turn tools, in milliseconds. */
  int turn_wait_ms( const json &arguments, const int default_s, const int max_s )
  {
    int wait_ms = default_s * 1000;
    if( arguments.contains( "wait_seconds" ) && arguments["wait_seconds"].is_number() )
      wait_ms = static_cast<int>( 1000.0 * std::clamp( arguments["wait_seconds"].get<double>(), 0.0,
                                                        static_cast<double>(max_s) ) );
    return wait_ms;
  }//turn_wait_ms(...)


  /** Drives a `PendingReply` for `LlmMcpDispatcher::deliverWhenReady`: sends its progress
   notifications, and then its reply.  Kept alive by its pending timer wait.

   It uses its own timer, rather than `WIOService::schedule()`, so the wait can be cancelled once the
   reply has been sent (or at shutdown) - Wt's IO service does not finish stopping until every
   pending timer has fired.
   */
  struct Delivery : public std::enable_shared_from_this<Delivery>
  {
    std::shared_ptr<LlmMcpDispatcher::PendingReply> reply;
    std::function<void( const json & )> send;

    std::mutex mutex;
    bool finished = false;
    std::unique_ptr<boost::asio::steady_timer> timer;

    /** Stops waiting, without sending anything more. */
    void cancel()
    {
      std::lock_guard<std::mutex> lock( mutex );
      finished = true;
      if( timer )
        timer->cancel();
    }

    void step()
    {
      std::lock_guard<std::mutex> lock( mutex );
      if( finished )
        return;

      try
      {
        Wt::WServer * const server = Wt::WServer::instance();
        const uint64_t version = reply->version();

        if( reply->streaming )
        {
          for( const json &notification : reply->takeProgressNotifications() )
            send( notification );
        }

        const bool deadline_passed = (std::chrono::steady_clock::now() >= reply->deadline);
        const std::optional<json> message = reply->poll( deadline_passed || !server );
        if( message )
        {
          finished = true;
          if( timer )
            timer->cancel();
          send( *message );
          return;
        }

        // Wait again: woken by a change (rescheduled, so a listener called right away does not
        //  re-enter this function), or by the timer - at the deadline, or sooner to send a heartbeat.
        const std::weak_ptr<Delivery> weak_self = shared_from_this();
        reply->listen( this, version, [weak_self](){
          Wt::WServer * const server = Wt::WServer::instance();
          if( server )
            server->ioService().post( [weak_self](){
              if( const std::shared_ptr<Delivery> self = weak_self.lock() )
                self->step();
            } );
        } );

        const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        long long delay_ms = std::chrono::duration_cast<std::chrono::milliseconds>( reply->deadline - now ).count();
        delay_ms = std::min<long long>( delay_ms, LlmMcpDispatcher::sm_keepalive_ms );
        delay_ms = std::max<long long>( delay_ms, LlmMcpDispatcher::sm_min_wake_ms );

        if( !timer )
          timer = std::make_unique<boost::asio::steady_timer>( server->ioService() );

        // (Re-arming cancels an earlier wait, whose handler then sees operation_aborted.)
        timer->expires_after( std::chrono::milliseconds( delay_ms ) );
        const std::shared_ptr<Delivery> self = shared_from_this();
        timer->async_wait( [self]( const boost::system::error_code &ec ){
          if( !ec )
            self->step();
        } );
      }catch( const std::exception &e )
      {
        finished = true;
        if( timer )
          timer->cancel();
        send( LlmMcpDispatcher::rpcError( reply->request_id, -32603, std::string("Internal error: ") + e.what() ) );
      }
    }//void step()
  };//struct Delivery


  /** The deliveries not yet finished, so they can be cancelled at shutdown. */
  std::mutex ns_deliveries_mutex;
  std::vector<std::weak_ptr<Delivery>> ns_deliveries;
}//namespace


/** The assistant turns started through MCP, so they can be polled and cancelled. */
class LlmMcpDispatcher::TurnStore
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
};//class LlmMcpDispatcher::TurnStore


/** Conversation histories for MCP tool calls in sessions without an LLM assistant (e.g., a desktop
 session whose assistant was never opened), so tools that keep state across calls - peak and
 energy-calibration checkpoints - work there too.  Keyed by Wt session id; entries unused for a day
 are dropped.
 */
class LlmMcpDispatcher::HistoryStore
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
};//class LlmMcpDispatcher::HistoryStore


std::vector<json> LlmMcpDispatcher::PendingReply::takeProgressNotifications()
{
  std::vector<json> notifications;
  if( !streaming || !drain_progress )
    return notifications;

  std::vector<json> items = drain_progress();
  if( progress_token.is_null() )
    return notifications;

  if( items.empty() )
    items.push_back( json() );   // a heartbeat

  for( const json &item : items )
  {
    std::string message = progress_message( item );
    message = what + ": " + (message.empty() ? std::string("still working") : message);
    notifications.push_back( json{ {"jsonrpc", "2.0"}, {"method", "notifications/progress"},
      {"params", { {"progressToken", progress_token},
                   {"progress", ++progress_sent},
                   {"message", message} }} } );
  }

  return notifications;
}//takeProgressNotifications()


std::string LlmMcpDispatcher::dumpMessage( const json &message )
{
  return message.dump( -1, ' ', false, json::error_handler_t::replace );
}


json LlmMcpDispatcher::rpcResult( const json &id, json tool_result )
{
  if( !tool_result.is_object() )
    tool_result = json{ {"result", std::move(tool_result)} };
  return json{ {"jsonrpc", "2.0"}, {"id", id}, {"result", { {"content", formatMcpToolResultContent( tool_result )} }} };
}


json LlmMcpDispatcher::rpcError( const json &id, const int code, const std::string &message )
{
  return json{ {"jsonrpc", "2.0"}, {"id", id}, {"error", { {"code", code}, {"message", message} }} };
}


json LlmMcpDispatcher::toolErrorResult( const json &id, const std::string &message )
{
  const json content = json::array( { json{ {"type", "text"}, {"text", message} } } );
  return json{ {"jsonrpc", "2.0"}, {"id", id}, {"result", { {"content", content}, {"isError", true} }} };
}


const std::vector<std::string> &LlmMcpDispatcher::supportedProtocolVersions()
{
  // The revisions with an `initialize` handshake; they differ for us only in optional features.
  //  (2026-07-28 replaced the handshake with per-request metadata, and is not implemented; clients
  //  that also speak the older revisions fall back to `initialize` when `server/discover` fails.)
  static const std::vector<std::string> sm_versions{ "2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05" };
  return sm_versions;
}


void LlmMcpDispatcher::deliverWhenReady( const std::shared_ptr<PendingReply> &pending,
                                      std::function<void( const json & )> send )
{
  assert( pending && send );
  const std::shared_ptr<Delivery> delivery = std::make_shared<Delivery>();
  delivery->reply = pending;
  delivery->send = std::move( send );

  {
    std::lock_guard<std::mutex> lock( ns_deliveries_mutex );
    ns_deliveries.erase( std::remove_if( begin(ns_deliveries), end(ns_deliveries),
                                         []( const std::weak_ptr<Delivery> &d ){ return d.expired(); } ),
                         end(ns_deliveries) );
    ns_deliveries.push_back( delivery );
  }

  delivery->step();
}//deliverWhenReady(...)


void LlmMcpDispatcher::cancelDeliveries()
{
  std::vector<std::weak_ptr<Delivery>> deliveries;
  {
    std::lock_guard<std::mutex> lock( ns_deliveries_mutex );
    deliveries.swap( ns_deliveries );
  }

  for( const std::weak_ptr<Delivery> &weak_delivery : deliveries )
  {
    if( const std::shared_ptr<Delivery> delivery = weak_delivery.lock() )
      delivery->cancel();
  }
}//cancelDeliveries()


LlmMcpDispatcher::LlmMcpDispatcher( const std::shared_ptr<const LlmConfig> &config, Options options,
                              SessionResolver resolver )
  : m_config( config ),
    m_options( std::move(options) ),
    m_resolver( std::move(resolver) ),
    m_registry( config ? std::make_shared<LlmTools::ToolRegistry>( *config )
                       : std::shared_ptr<LlmTools::ToolRegistry>() ),
    m_turns( std::make_shared<TurnStore>() ),
    m_histories( std::make_shared<HistoryStore>() ),
    m_tools{},
    m_interceptor{}
{
  if( !m_config )
    throw std::logic_error( "LlmMcpDispatcher must be initialized with a valid config" );
  if( !m_resolver )
    throw std::logic_error( "LlmMcpDispatcher must be given a session resolver" );

  register_default_tools();
}//LlmMcpDispatcher constructor


LlmMcpDispatcher::~LlmMcpDispatcher()
{
}


void LlmMcpDispatcher::addTools( const std::vector<ToolInfo> &tools, ToolInterceptor interceptor )
{
  for( const ToolInfo &tool : tools )
    m_tools[tool.name] = tool;
  m_interceptor = std::move( interceptor );
}//addTools(...)


LlmMcpDispatcher::Outcome LlmMcpDispatcher::handle( const json &message, const RequestContext &context )
{
  Outcome outcome;
  json id = nullptr;  // Track request id for error responses

  try
  {
    if( !message.is_object() || !message.contains( "method" ) || !message["method"].is_string() )
    {
      outcome.reply = rpcError( nullptr, -32600, "Invalid Request: missing 'method' field" );
      return outcome;
    }

    const std::string method = message["method"];
    id = message.value( "id", json( nullptr ) );

    // --- initialize ---
    if( method == "initialize" )
    {
      // Reply with the client's version if we support it, and otherwise our latest.
      const std::vector<std::string> &versions = supportedProtocolVersions();
      std::string version = versions.front();
      if( message.contains( "params" ) && message["params"].is_object() )
      {
        const std::string requested = message["params"].value( "protocolVersion", std::string() );
        if( std::find( begin(versions), end(versions), requested ) != end(versions) )
          version = requested;
      }

      json result = {
        {"protocolVersion", version},
        {"capabilities", {
          {"tools", json::object()}
        }},
        {"serverInfo", {
          {"name", "InterSpec MCP Server"},
          {"version", InterSpec_VERSION}
        }}
      };

      if( !m_options.instructions.empty() )
        result["instructions"] = m_options.instructions;

      outcome.reply = json{ {"jsonrpc", "2.0"}, {"id", id}, {"result", result} };
      return outcome;
    }//if( method == "initialize" )

    // --- notifications (notifications/initialized, notifications/cancelled, ...): no reply ---
    if( SpecUtils::starts_with( method, "notifications/" ) || !message.contains( "id" ) )
      return outcome;

    // --- ping ---
    if( method == "ping" )
    {
      outcome.reply = json{ {"jsonrpc", "2.0"}, {"id", id}, {"result", json::object()} };
      return outcome;
    }

    // --- tools/list ---
    if( method == "tools/list" )
    {
      json tool_definitions = json::array();
      for( const auto &[name, info] : m_tools )
      {
        tool_definitions.push_back( {
          {"name", info.name},
          {"description", info.description},
          {"inputSchema", info.parameters_schema}
        } );
      }

      outcome.reply = json{ {"jsonrpc", "2.0"}, {"id", id}, {"result", {{"tools", tool_definitions}}} };
      return outcome;
    }//if( method == "tools/list" )

    // --- tools/call ---
    if( method == "tools/call" )
    {
      if( !message.contains( "params" ) )
      {
        outcome.reply = rpcError( id, -32602, "Missing 'params' for tools/call" );
        return outcome;
      }

      const json &call_params = message["params"];
      if( !call_params.is_object() || !call_params.contains( "name" ) || !call_params["name"].is_string() )
      {
        outcome.reply = rpcError( id, -32602, "Missing or invalid 'name' for tools/call" );
        return outcome;
      }
      const std::string tool_name = call_params["name"].get<std::string>();
      json arguments = call_params.value( "arguments", json::object() );
      if( !arguments.is_object() )
        arguments = json::object();

      json progress_token;
      if( call_params.contains( "_meta" ) && call_params["_meta"].is_object()
          && call_params["_meta"].contains( "progressToken" ) )
        progress_token = call_params["_meta"]["progressToken"];

      try
      {
        return handle_tool_call( id, tool_name, std::move(arguments), progress_token, context );
      }catch( const std::exception &e )
      {
        // The call could not be started (no such session or turn, a bad argument, ...).
        outcome = Outcome{};
        outcome.reply = toolErrorResult( id, e.what() );
        return outcome;
      }
    }//if( method == "tools/call" )

    // --- Unknown method ---
    outcome.reply = rpcError( id, -32601, "Method not found: " + method );
  }catch( const json::exception &e )
  {
    // JSON access error
    outcome = Outcome{};
    outcome.reply = rpcError( id, -32700, std::string( "Parse error: " ) + e.what() );
  }catch( const std::exception &e )
  {
    // Tool execution or other error
    outcome = Outcome{};
    outcome.reply = rpcError( id, -32603, e.what() );
  }

  return outcome;
}//handle(...)


LlmMcpDispatcher::Outcome LlmMcpDispatcher::handle_tool_call( const json &id, const std::string &tool_name,
                                                        json arguments, const json &progress_token,
                                                        const RequestContext &context )
{
  Outcome outcome;

  if( m_tools.find( tool_name ) == m_tools.end() )
  {
    outcome.reply = rpcError( id, -32602, "Unknown tool: " + tool_name );
    return outcome;
  }

  // Extract and remove userSession from arguments before passing to tool
  std::string user_session;
  if( m_options.session_argument && arguments.contains( "userSession" ) && arguments["userSession"].is_string() )
    user_session = arguments["userSession"].get<std::string>();
  arguments.erase( "userSession" );

  // Set by us for turns; never taken from a client.
  arguments.erase( "_mcp_turn_id" );

  if( m_interceptor )
  {
    std::optional<Outcome> intercepted = m_interceptor( id, tool_name, user_session, arguments, context );
    if( intercepted )
      return std::move( *intercepted );
  }

  const std::string wt_session_id = m_resolver( user_session ); // throws with a message for the client

#if( BUILD_AS_MCP_HTTP_SERVER )
  const McpAgentTools::ToolDef * const agent_tool = McpAgentTools::toolDefinition( tool_name );

  const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
  if( manager )
    manager->touch( user_session );

  // The session a turn belongs to, for later assistant_get_result / assistant_cancel calls.
  const std::string &session_key = user_session;
#else
  const std::string &session_key = wt_session_id;
#endif

  if( tool_name == "assistant_get_result" )
    return handle_turn_result( id, session_key, arguments, progress_token, context );

  const LlmTools::SharedTool * const tool = m_registry->getTool( tool_name );
  bool known_tool = (tool || (tool_name == "assistant_cancel"));
#if( BUILD_AS_MCP_HTTP_SERVER )
  known_tool = (known_tool || agent_tool);
#endif
  if( !known_tool )
    throw std::runtime_error( "Tool not found: " + tool_name );

  const bool is_turn = ((tool_name == "assistant_submit_prompt") || (tool_name == "invoke_agent"));

  // How long a non-streaming client waits for a turn, before being told it is still running.
  const int wait_ms = turn_wait_ms( arguments, sm_default_turn_wait_s, sm_max_turn_wait_s );
  arguments.erase( "wait_seconds" );

  std::string turn_id;
  std::shared_ptr<PendingCall> pending, cancel_target;

  if( is_turn )
  {
    turn_id = "turn-" + Wt::WRandom::generateId( 16 );
    pending = m_turns->create( turn_id, session_key );
    arguments["_mcp_turn_id"] = turn_id;
  }else if( tool_name == "assistant_cancel" )
  {
    const std::string target_turn = arguments.value( "turn_id", std::string() );
    const std::optional<TurnStore::Turn> turn = m_turns->find( target_turn );
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
  const std::shared_ptr<const LlmTools::ToolRegistry> registry = m_registry;
  const std::shared_ptr<HistoryStore> histories = m_histories;

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

  const std::shared_ptr<PendingReply> reply = std::make_shared<PendingReply>();
  reply->request_id = id;
  reply->streaming = context.can_stream && (is_turn || (tool && tool->isAsync()));
  // A turn returns its turn_id within wait_seconds even when streaming: clients commonly give up on
  //  a call after a few minutes, and without the turn_id they could neither poll nor cancel it.
  reply->deadline = std::chrono::steady_clock::now()
                    + std::chrono::milliseconds( is_turn ? wait_ms : tool_timeout_ms );
  reply->progress_token = progress_token;
  reply->what = is_turn ? ("Assistant turn " + turn_id) : ("Tool " + tool_name);
  wait_on( *reply, pending );

  const bool session_argument = m_options.session_argument;
#if( BUILD_AS_MCP_HTTP_SERVER )
  const std::string base_url = context.base_url;
#endif

  reply->poll = [=]( const bool deadline_passed ) -> std::optional<json> {
    if( !pending->isDone() )
    {
      if( !deadline_passed )
        return std::nullopt;

      if( is_turn )
        return rpcResult( id, running_turn_status( *pending, turn_id,
          "The assistant is still working; call assistant_get_result with this turn_id to wait for the result." ) );

      return toolErrorResult( id, LlmTools::async_timeout_error_message( tool_name, tool_timeout_ms ) );
    }//if( !pending->isDone() )

    std::variant<json, std::string> result_or_error = pending->getResult();

    if( const std::string * const err = std::get_if<std::string>( &result_or_error ) )
    {
      if( is_turn )
        return rpcResult( id, json{ {"status", pending->cancelRequested() ? "cancelled" : "error"},
                                    {"turn_id", turn_id}, {"error", *err} } );
      return toolErrorResult( id, *err );
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
            const std::string ext = (mime == "image/jpeg") ? ".jpg" : ".png";
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
    if( session_argument )
      result["userSession"] = wt_session_id;
#endif

    return rpcResult( id, std::move(result) );
  };//reply->poll

  outcome.pending = reply;
  return outcome;
}//handle_tool_call(...)


LlmMcpDispatcher::Outcome LlmMcpDispatcher::handle_turn_result( const json &id, const std::string &session_key,
                                                          const json &arguments, const json &progress_token,
                                                          const RequestContext &context )
{
  const std::string turn_id = arguments.value( "turn_id", std::string() );
  const std::optional<TurnStore::Turn> turn = m_turns->find( turn_id );
  if( !turn || (turn->session_key != session_key) )
    throw std::runtime_error( "Unknown turn_id '" + turn_id + "' for this session (results are kept for "
                              + std::to_string(sm_turn_retention.count()) + " minutes)." );

  const int wait_ms = turn_wait_ms( arguments, sm_default_turn_wait_s, sm_max_turn_wait_s );
  const bool session_argument = m_options.session_argument;

  const std::shared_ptr<PendingCall> pending = turn->call;
  const std::shared_ptr<PendingReply> reply = std::make_shared<PendingReply>();
  reply->request_id = id;
  reply->deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds( wait_ms );
  reply->streaming = context.can_stream;
  reply->progress_token = progress_token;
  reply->what = "Assistant turn " + turn_id;
  wait_on( *reply, pending );

  reply->poll = [pending, id, turn_id, session_key, session_argument]( const bool deadline_passed ) -> std::optional<json> {
    if( !pending->isDone() )
    {
      if( !deadline_passed )
        return std::nullopt;
      return rpcResult( id, running_turn_status( *pending, turn_id,
                                                 "Still working; call assistant_get_result again." ) );
    }

    std::variant<json, std::string> result_or_error = pending->getResult();
    if( const std::string * const err = std::get_if<std::string>( &result_or_error ) )
      return rpcResult( id, json{ {"status", pending->cancelRequested() ? "cancelled" : "error"},
                                  {"turn_id", turn_id}, {"error", *err} } );

    json result = std::get<json>( std::move(result_or_error) );
    if( !result.is_object() )
      result = json{ {"result", std::move(result)} };
    result["status"] = "done";
    result["turn_id"] = turn_id;
    if( session_argument )
      result["userSession"] = session_key;
    return rpcResult( id, std::move(result) );
  };

  Outcome outcome;
  outcome.pending = reply;
  return outcome;
}//handle_turn_result(...)


void LlmMcpDispatcher::add_session_param( json &schema ) const
{
  if( !m_options.session_argument || !schema.is_object() )
    return;

  if( !schema.contains( "properties" ) )
    schema["properties"] = json::object();

  if( !schema["properties"].contains( "userSession" ) )
    schema["properties"]["userSession"] = { {"type", "string"},
                                            {"description", m_options.session_argument_description} };

  if( m_options.session_argument_required )
  {
    json &required = schema["required"];
    if( !required.is_array() )
      required = json::array();
    if( std::find( begin(required), end(required), json("userSession") ) == end(required) )
      required.push_back( "userSession" );
  }
}//add_session_param(...)


void LlmMcpDispatcher::register_default_tools()
{
  assert( m_registry && m_config );

  // Tools filtered for MCP (excludes invoke_*, KnowledgeBase, agent-only tools)
  for( const auto &[name, sharedTool] : m_registry->getToolsForMcp() )
  {
    if( sharedTool->requiresBrowser && !m_options.browser_available )
      continue;

    ToolInfo info;
    info.name = sharedTool->name;
    // Use the MCP-surface rendering of the description (rendered at registry build); fall back to the
    // agent-surface description if unset.  For marker-free descriptions these are identical.
    info.description = sharedTool->mcpDescription.empty() ? sharedTool->description
                                                          : sharedTool->mcpDescription;
    info.parameters_schema = sharedTool->mcpParametersSchema.is_null() ? sharedTool->parameters_schema
                                                                       : sharedTool->mcpParametersSchema;
    add_session_param( info.parameters_schema );

    m_tools[name] = std::move( info );
  }//for( loop over MCP-filtered tools )

  // Polling and cancelling the assistant turns that the prompt and agent tools start, whenever one
  //  of those is offered (see `<McpServer><AllowAssistantPrompt>` / `<AllowAgentInvocation>`).
  if( m_tools.count( "assistant_submit_prompt" ) || m_tools.count( "invoke_agent" ) )
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
    m_tools[result_tool.name] = std::move( result_tool );

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
    m_tools[cancel_tool.name] = std::move( cancel_tool );
  }//if( an assistant turn tool is offered )

#if( PERFORM_DEVELOPER_CHECKS )
  std::cerr << "MCP dispatcher registered " << m_tools.size() << " tools from shared registry" << std::endl;
#endif
}//register_default_tools()
