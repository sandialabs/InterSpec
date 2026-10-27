#ifndef LLM_MCP_DISPATCHER_H
#define LLM_MCP_DISPATCHER_H
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
#include <mutex>
#include <atomic>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <optional>
#include <functional>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_MCP_SERVER, "You should not include this header unless USE_MCP_SERVER is enabled" );

// Forward declarations
class LlmConfig;

namespace LlmTools
{
  class ToolRegistry;
}


/** The transport-independent part of InterSpec's MCP server: takes one JSON-RPC message, and gives
 back the reply - right away, or later through a `PendingReply` the transport waits on.

 It handles `initialize` (negotiating the protocol version), `ping`, `tools/list` and `tools/call`.
 Tools run in an InterSpec session, scheduled into it with a delayed `WServer::schedule` (not
 `WServer::post`, whose server-wide strand would let one session's long tool call hold up every
 other session's work), so `handle()` never blocks on a tool, and must not be called with a
 session's lock held.

 Assistant prompts (`assistant_submit_prompt`, `invoke_agent`) become "turns" with an id; the call
 waits up to `wait_seconds` (reporting progress meanwhile, if the transport can), and the client then
 polls `assistant_get_result`, or cancels with `assistant_cancel`.

 Used by `LlmMcpResource` (HTTP, at `/mcp-api`), and by the stdio server (target/mcp_stdio).

 This implements the MCP revisions with an `initialize` handshake (2024-11-05 through 2025-11-25).
 */
class LlmMcpDispatcher
{
public:
  /** Shared state of the assistant turns; see the .cpp. */
  class TurnStore;

  /** Per-session conversation histories, for tools that keep state (checkpoints) when the session
   has no LLM assistant to share one with; see the .cpp.
   */
  class HistoryStore;

  /** How the server presents itself; fixed for the life of the dispatcher. */
  struct Options
  {
    /** If true, every session tool takes a `userSession` argument naming the session to act on
     (passed to the `SessionResolver`), and results say which session they came from.  Off for a
     server with only one session (stdio).
     */
    bool session_argument = true;

    /** If true, `userSession` is a required argument (the multi-session agent server). */
    bool session_argument_required = false;

    /** Description of the `userSession` argument. */
    std::string session_argument_description;

    /** The `instructions` returned by `initialize`, if not empty. */
    std::string instructions;

    /** If false, the sessions have no browser attached (stdio), so tools that need one
     (`SharedTool::requiresBrowser`) are not offered.
     */
    bool browser_available = true;
  };//struct Options

  /** What only the transport knows about a request. */
  struct RequestContext
  {
    /** If true, progress notifications can be sent while a call runs (an SSE response, or stdio). */
    bool can_stream = false;

    /** Base of URLs given to clients (no trailing slash), e.g. "http://host:8078"; agent server only. */
    std::string base_url;
  };//struct RequestContext

  /** Returns the Wt session id to run a session tool in, given the `userSession` argument (empty if
   not given).  Called on the transport's thread, not in a session.
   @throws std::runtime_error with a message for the client.
   */
  using SessionResolver = std::function<std::string( const std::string &user_session )>;


  /** A reply that is not ready yet: the transport waits for it, by registering with `listen` (and
   with a timer, for the deadline and heartbeats), and calling `poll` when woken.
   */
  struct PendingReply
  {
    nlohmann::json request_id;

    /** If true, progress is reported while waiting (`takeProgressNotifications`), and an HTTP reply
     is an SSE stream.
     */
    bool streaming = false;

    std::chrono::steady_clock::time_point deadline;

    /** The current version of what is waited on; see `listen`. */
    std::function<uint64_t()> version;

    /** Registers a one-shot wake-up for the first change after `seen_version` (called right away if
     there already was one); a later registration with the same `key` replaces an earlier one.
     */
    std::function<void( const void *key, uint64_t seen_version, std::function<void()> )> listen;

    /** Returns the complete JSON-RPC message to send, or nullopt if not ready; once the deadline has
     passed (`deadline_passed`) it must return a message.
     */
    std::function<std::optional<nlohmann::json>( bool deadline_passed )> poll;

    /** Progress items (see `LlmToolGui::buildProgressFromTurn`) to report before the result. */
    std::function<std::vector<nlohmann::json>()> drain_progress;

    /** The client's `_meta.progressToken`, or null.  Progress notifications are only sent when the
     client gave one (as MCP requires).
     */
    nlohmann::json progress_token;

    /** Number of progress notifications sent; MCP requires `progress` to increase. */
    uint64_t progress_sent = 0;

    /** What is being waited on, e.g. "Assistant turn turn-abc", for progress messages. */
    std::string what;

    /** Incremented each time the transport starts a new wait; a wake-up timer from an earlier wait is
     then stale, and should do nothing.
     */
    std::atomic<uint64_t> wake_generation{ 0 };

    /** The `notifications/progress` messages to send now, for a streaming reply whose client gave a
     progress token: one per progress item since the last call, or else a heartbeat, so clients that
     time out calls with no activity keep waiting.  Empty otherwise.
     */
    std::vector<nlohmann::json> takeProgressNotifications();
  };//struct PendingReply


  /** The result of handling a message: a reply to send now, or one to wait for, or (for a
   notification) neither.
   */
  struct Outcome
  {
    std::optional<nlohmann::json> reply;
    std::shared_ptr<PendingReply> pending;
  };//struct Outcome


  /** Lets a transport handle some tools itself (the agent server's session management and file
   tools); called for each `tools/call`, with the `userSession` argument already removed from
   `arguments`.  Returns nullopt to let the dispatcher handle the call.
   @throws std::exception with a message for the client.
   */
  using ToolInterceptor = std::function<std::optional<Outcome>( const nlohmann::json &id,
                                                                const std::string &tool_name,
                                                                const std::string &user_session,
                                                                nlohmann::json &arguments,
                                                                const RequestContext &context )>;

  /** A tool offered to clients, besides those of the tool registry. */
  struct ToolInfo
  {
    std::string name;
    std::string description;
    nlohmann::json parameters_schema;
  };//struct ToolInfo


  /** @throws std::exception if `config` is null. */
  LlmMcpDispatcher( const std::shared_ptr<const LlmConfig> &config, Options options,
                 SessionResolver resolver );

  ~LlmMcpDispatcher();

  /** Offers tools handled by the `ToolInterceptor` (their schemas are given as-is). */
  void addTools( const std::vector<ToolInfo> &tools, ToolInterceptor interceptor );

  /** Handles one JSON-RPC message; never throws (errors become JSON-RPC error replies). */
  Outcome handle( const nlohmann::json &message, const RequestContext &context );

  /** Serializes a message to send, on a single line; invalid UTF-8 (e.g., from a file's text) is
   replaced rather than throwing, so the client always gets its reply.
   */
  static std::string dumpMessage( const nlohmann::json &message );

  /** A `tools/call` reply carrying `tool_result` (an "_image" in it becomes an MCP image block). */
  static nlohmann::json rpcResult( const nlohmann::json &id, nlohmann::json tool_result );

  static nlohmann::json rpcError( const nlohmann::json &id, int code, const std::string &message );

  /** A `tools/call` reply for a tool that failed: `message` as text, with `isError` set.  MCP reports
   tool failures this way (rather than as JSON-RPC errors, which are for protocol problems), so the
   model sees the message and can correct its call.
   */
  static nlohmann::json toolErrorResult( const nlohmann::json &id, const std::string &message );

  /** For message-based transports (stdio): sends `pending`'s progress notifications as they come,
   and then its reply, through `send` (called from Wt IO threads, one call at a time per reply).
   Returns right away.
   */
  static void deliverWhenReady( const std::shared_ptr<PendingReply> &pending,
                                std::function<void( const nlohmann::json & )> send );

  /** Stops every `deliverWhenReady` still waiting, without sending anything more; for shutdown,
   before the Wt server is destroyed.
   */
  static void cancelDeliveries();

  /** The MCP revisions supported, newest first. */
  static const std::vector<std::string> &supportedProtocolVersions();

  /** Interval between heartbeats while waiting on a reply; short enough to keep clients and
   intermediaries from timing the call out.
   */
  static constexpr int sm_keepalive_ms = 10000;

  /** Floor for how soon a wait is woken by its timer, so a nearly-expired deadline can not spin. */
  static constexpr int sm_min_wake_ms = 100;

private:
  Outcome handle_tool_call( const nlohmann::json &id, const std::string &tool_name, nlohmann::json arguments,
                            const nlohmann::json &progress_token, const RequestContext &context );

  /** `assistant_get_result`: waits for the turn (of the session `session_key`) to finish, up to
   `wait_seconds`.
   */
  Outcome handle_turn_result( const nlohmann::json &id, const std::string &session_key,
                              const nlohmann::json &arguments, const nlohmann::json &progress_token,
                              const RequestContext &context );

  /** Populates `m_tools` from the tool registry, filtered for MCP. */
  void register_default_tools();

  /** Adds the `userSession` argument to a tool schema, if `Options::session_argument`. */
  void add_session_param( nlohmann::json &schema ) const;

  /** Default, and maximum, of the `wait_seconds` argument of assistant turns. */
  static constexpr int sm_default_turn_wait_s = 30;
  static constexpr int sm_max_turn_wait_s = 50;

  const std::shared_ptr<const LlmConfig> m_config;
  const Options m_options;
  const SessionResolver m_resolver;

  /** Shared with the work scheduled into sessions, which may outlive this dispatcher. */
  const std::shared_ptr<const LlmTools::ToolRegistry> m_registry;

  const std::shared_ptr<TurnStore> m_turns;
  const std::shared_ptr<HistoryStore> m_histories;

  /** The tools offered, keyed by name. */
  std::map<std::string, ToolInfo> m_tools;

  ToolInterceptor m_interceptor;
};//class LlmMcpDispatcher

#endif //LLM_MCP_DISPATCHER_H
