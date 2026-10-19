#ifndef LLM_MCP_RESOURCE_H
#define LLM_MCP_RESOURCE_H
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
#include <memory>
#include <string>

#include <Wt/WResource.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) || defined(ANDROID) )
  static_assert( 0, "MCP server can not be enabled for web or mobile builds" );
#endif

// Forward declarations
class LlmConfig;

namespace LlmTools
{
  class ToolRegistry;
}

namespace Wt
{
  namespace Http
  {
    class Request;
    class Response;
  }
}


/** The MCP server: JSON-RPC 2.0 over HTTP at `/mcp-api`, exposing the tools of
 `LlmTools::ToolRegistry` that are available to MCP clients.

 Each tool call runs in one InterSpec session (named by the `userSession` argument), scheduled into
 it with a delayed `WServer::schedule` (not `WServer::post`, whose server-wide strand would let one
 session's long tool call hold up every other session).  The HTTP request waits for the result with a
 `ResponseContinuation`, so no HTTP worker thread is blocked while a tool (or a whole assistant turn)
 runs.  Clients that accept `text/event-stream` get async tools' progress streamed as SSE events.

 Assistant prompts (`assistant_submit_prompt`, `invoke_agent`) become "turns" with an id; the call
 waits up to `wait_seconds` (streaming MCP progress meanwhile, if the client asked for it), and the
 client then polls `assistant_get_result`, or cancels with `assistant_cancel`.

 Built as the multi-session agent server (`BUILD_AS_MCP_HTTP_SERVER`; see target/mcp_http), this
 also serves:
 - the tools of `McpAgentTools` (`create_session`, `close_session`, file exchange, exporting, ...),
   and `userSession` is required by every session tool, naming a session from `create_session`;
 - `PUT|POST /mcp-api/files/<upload token>[/<filename>]` to upload into a sessions `input/`;
 - `GET /mcp-api/files/<download token>/<filename>` to download a sessions output file;
 - `GET /mcp-api/healthz`, unauthenticated, with session counts only.
 All but the health check need the bearer token, and no CORS headers are sent, as only a back-end
 should be talking to an agent server.
 */
class LlmMcpResource : public Wt::WResource
{
public:
  /** May throw exception. */
  LlmMcpResource( const std::shared_ptr<const LlmConfig> &config );

  virtual ~LlmMcpResource();

  /** Shared state of the assistant turns started through MCP; see the .cpp. */
  class TurnStore;

  /** Whats being waited on by a suspended request; see the .cpp. */
  struct DeferredCall;

  /** Per-session conversation histories for tools that keep state (checkpoints) when the session
   has no LLM assistant to share one with; see the .cpp. */
  class HistoryStore;

protected:
  void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response ) override;

private:
  /** Tool metadata registered for MCP clients (name, description, schema). */
  struct ToolInfo
  {
    std::string name;
    std::string description;
    nlohmann::json parameters_schema;  // Includes the added "userSession" field
  };

  /** Populates `registered_tools_` from the shared tool registry, filtered for MCP, plus (for an
   agent server) the `McpAgentTools`.
   */
  void register_default_tools();

#if( MCP_ENABLE_AUTH )
  /** Validates the bearer token in the request (compared in constant time).
   If invalid, sets response status to 401 and writes a JSON error body.
   @returns true if the validation passed (or no token is configured).
   */
  bool validate_request( const Wt::Http::Request &request, Wt::Http::Response &response );
#endif

  /** Returns auth type info for MCP discovery (GET on root). */
  void handle_get_auth( const Wt::Http::Request &request, Wt::Http::Response &response );

  /** Handles all JSON-RPC 2.0 MCP methods (initialize, tools/list, tools/call, etc). */
  void handle_jsonrpc( const Wt::Http::Request &request, Wt::Http::Response &response );

  /** Handles a `tools/call`; writes the response, or suspends the request until it is ready.
   `progress_token` is the request's `_meta.progressToken` (null if none).
   */
  void handle_tool_call( const nlohmann::json &id, const std::string &tool_name, nlohmann::json arguments,
                         const nlohmann::json &progress_token,
                         const Wt::Http::Request &request, Wt::Http::Response &response );

  /** Handles `assistant_get_result`: waits for the turn (of the session `session_key`) to finish,
   up to `wait_seconds`.
   */
  void handle_turn_result( const nlohmann::json &id, const std::string &session_key,
                           const nlohmann::json &arguments, const nlohmann::json &progress_token,
                           const Wt::Http::Request &request, Wt::Http::Response &response );

#if( BUILD_AS_MCP_HTTP_SERVER )
  /** Handles the agent-server tools that are not run inside a session. */
  void handle_agent_tool( const nlohmann::json &id, const std::string &tool_name,
                          const std::string &user_session, nlohmann::json arguments,
                          const Wt::Http::Request &request, Wt::Http::Response &response );

  /** Agent server: `/files/...` uploads and downloads. */
  void handle_file_request( const std::string &path_info, const Wt::Http::Request &request,
                            Wt::Http::Response &response );

  /** Agent server: the unauthenticated health check. */
  void handle_health( Wt::Http::Response &response );

  /** Base for the upload/download URLs given to clients (no trailing slash), e.g. "http://host:8078". */
  std::string public_base_url( const Wt::Http::Request &request ) const;
#endif //BUILD_AS_MCP_HTTP_SERVER

  /** The Wt session id to run a session tool in: for an agent server, the live session named by
   `user_session` (required); otherwise the session with that id or external token, or the most
   recently started one if `user_session` is empty.

   @throws std::runtime_error with a message for the client.
   */
  std::string resolve_target_session( const std::string &user_session ) const;

  /** Writes the response now, if `call` is ready; otherwise suspends the request until it is. */
  void start_deferred( Wt::Http::Response &response, const std::shared_ptr<DeferredCall> &call );

  /** Re-entry of a suspended request. */
  void continue_deferred( Wt::Http::Response &response, const std::shared_ptr<DeferredCall> &call );
  void continue_deferred_impl( Wt::Http::Response &response, const std::shared_ptr<DeferredCall> &call );

  /** Suspends the request, to be woken when `call` changes, at its deadline, or (for SSE) to send a
   keepalive.
   */
  void suspend( Wt::Http::Response &response, const std::shared_ptr<DeferredCall> &call,
                const uint64_t seen_version );

  /** Interval between `:keepalive` comments written to a waiting SSE stream; short enough to keep
   clients and intermediaries from timing the connection out.
   */
  static constexpr int sm_sse_keepalive_ms = 10000;

  /** Default, and maximum, of the `wait_seconds` argument of assistant turns. */
  static constexpr int sm_default_turn_wait_s = 30;
  static constexpr int sm_max_turn_wait_s = 50;

  const std::shared_ptr<const LlmConfig> llm_config_;

  /** Shared with the work posted to sessions, which may outlive this resource. */
  const std::shared_ptr<const LlmTools::ToolRegistry> tool_registry_;

  const std::shared_ptr<TurnStore> turns_;

  const std::shared_ptr<HistoryStore> histories_;

  /** Tool metadata map, keyed by tool name. Schema has "userSession" field added. */
  std::map<std::string, ToolInfo> registered_tools_;

  /** For the agent server, the `instructions` returned by `initialize`. */
  std::string server_instructions_;

#if( MCP_ENABLE_AUTH )
  /** Bearer token required to use the MCP server. If empty, no authorization is required. */
  const std::string secret_token_;
#endif
};//class LlmMcpResource

#endif //LLM_MCP_RESOURCE_H
