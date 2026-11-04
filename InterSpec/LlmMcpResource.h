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

#include <memory>
#include <string>
#include <optional>

#include <Wt/WResource.h>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "InterSpec/LlmMcpDispatcher.h"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) || defined(ANDROID) )
  static_assert( 0, "MCP server can not be enabled for web or mobile builds" );
#endif

// Forward declarations
class LlmConfig;

namespace Wt
{
  namespace Http
  {
    class Request;
    class Response;
  }
}


/** The MCP server over HTTP: JSON-RPC 2.0 at `/mcp-api`, with the messages themselves handled by
 `LlmMcpDispatcher`; this class adds the bearer-token check, CORS, and the waiting.

 A reply that is not ready (a tool still running in its session, or an assistant turn) is waited for
 with a `ResponseContinuation`, so no HTTP worker thread is blocked meanwhile.  Clients that accept
 `text/event-stream` get async tools' progress streamed as SSE events.

 Each tool call runs in one InterSpec session, named by the `userSession` argument (by default, the
 most recently started session).

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

protected:
  void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response ) override;

private:
#if( MCP_ENABLE_AUTH )
  /** Validates the bearer token in the request (compared in constant time).
   If invalid, sets response status to 401 and writes a JSON error body.
   @returns true if the validation passed (or no token is configured).
   */
  bool validate_request( const Wt::Http::Request &request, Wt::Http::Response &response );
#endif

  /** Returns auth type info for MCP discovery (GET on root). */
  void handle_get_auth( const Wt::Http::Request &request, Wt::Http::Response &response );

  /** Handles a POSTed JSON-RPC message. */
  void handle_jsonrpc( const Wt::Http::Request &request, Wt::Http::Response &response );

#if( BUILD_AS_MCP_HTTP_SERVER )
  /** The agent-server tools that are not run inside a session (the dispatcher's tool interceptor);
   returns nullopt for other tools.
   */
  std::optional<LlmMcpDispatcher::Outcome> handle_agent_tool( const nlohmann::json &id, const std::string &tool_name,
                                                           const std::string &user_session,
                                                           const nlohmann::json &arguments,
                                                           const std::string &base_url );

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
  static std::string resolve_target_session( const std::string &user_session );

  /** Writes the response now, if `call` is ready; otherwise suspends the request until it is. */
  void start_deferred( Wt::Http::Response &response, const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call );

  /** Re-entry of a suspended request. */
  void continue_deferred( Wt::Http::Response &response, const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call );
  void continue_deferred_impl( Wt::Http::Response &response,
                               const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call );

  /** Suspends the request, to be woken when `call` changes, at its deadline, or (for SSE) to send a
   keepalive.
   */
  void suspend( Wt::Http::Response &response, const std::shared_ptr<LlmMcpDispatcher::PendingReply> &call,
                const uint64_t seen_version );

  const std::shared_ptr<const LlmConfig> llm_config_;

  std::unique_ptr<LlmMcpDispatcher> dispatcher_;

#if( MCP_ENABLE_AUTH )
  /** Bearer token required to use the MCP server. If empty, no authorization is required. */
  const std::string secret_token_;
#endif
};//class LlmMcpResource

#endif //LLM_MCP_RESOURCE_H
