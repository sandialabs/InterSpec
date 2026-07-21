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
#include <atomic>
#include <future>
#include <memory>
#include <string>
#include <variant>

#include <Wt/WResource>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) || defined(ANDROID) )
  static_assert( 0, "MCP server can not be enabled for web or mobile builds" );
#endif

// Forward declarations
class InterSpec;
class InterSpecApp;
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


class LlmMcpResource : public Wt::WResource
{
public:
  /** May throw exception. */
  LlmMcpResource( const std::shared_ptr<const LlmConfig> &config );

  virtual ~LlmMcpResource();

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

  /** Populates `registered_tools_` from the shared tool registry, filtered for MCP. */
  void register_default_tools();

#if( MCP_ENABLE_AUTH )
  /** Validates the bearer token in the request.
   If invalid, sets response status to 401 and writes a JSON error body.
   @returns true if the validation passed (or no token is configured).
   */
  bool validate_request( const Wt::Http::Request &request, Wt::Http::Response &response );
#endif

  /** Returns auth type info for MCP discovery (GET on root). */
  void handle_get_auth( const Wt::Http::Request &request, Wt::Http::Response &response );

  /** Handles all JSON-RPC 2.0 MCP methods (initialize, tools/list, tools/call, etc). */
  void handle_jsonrpc( const Wt::Http::Request &request, Wt::Http::Response &response );

  /** Execute a tool within a WApplication session context, with proper locking.
   For sync tools, holds the WApplication::UpdateLock throughout execution.
   For async tools, releases the lock before blocking on the result future,
   to avoid deadlock with the async callback which re-acquires the lock via WServer::post.

   @param userSessionId  Session ID or external token to find. Empty = most recent session.
   @param toolName       Name of the tool to execute.
   @param params         JSON parameters for the tool (without userSession field).
   @param[out] actualSessionId  Set to the actual session ID used.
   @returns JSON result from the tool execution.
   @throws std::runtime_error on failure.
   */
  nlohmann::json executeWithSession( const std::string &userSessionId,
                                     const std::string &toolName,
                                     const nlohmann::json &params,
                                     std::string &actualSessionId );

  /** State carried across SSE continuation callbacks for async tool calls.
   Created when an async tool call starts on an SSE-capable client, and passed
   through `ResponseContinuation::setData()` to each resumption of `handleRequest`.
   */
  struct AsyncSseState
  {
    nlohmann::json request_id;
    std::string actual_session_id;
    std::shared_ptr<std::promise<std::variant<nlohmann::json, std::string>>> prom;
    std::shared_ptr<std::future<std::variant<nlohmann::json, std::string>>> fut;
    std::atomic<bool> tool_complete{false};
  };

  /** Handles continuation re-entries for SSE streaming of async tool results.
   Writes a keepalive comment if the tool is still running, or the final
   JSON-RPC response if the tool has completed.
   */
  void handle_sse_continuation( const Wt::Http::Request &request,
                                Wt::Http::Response &response );

  /** Schedules a single keepalive timer (~10s) that calls `haveMoreData()` to
   wake any waiting SSE continuations, so they can write a `:keepalive` comment
   and prevent client-side timeouts.
   */
  void schedule_sse_keepalive();

  const std::shared_ptr<const LlmConfig> llm_config_;
  const std::unique_ptr<const LlmTools::ToolRegistry> tool_registry_;

  /** Tool metadata map, keyed by tool name. Schema has "userSession" field added. */
  std::map<std::string, ToolInfo> registered_tools_;

  /** Guard for the keepalive timer lambda: set to false in destructor so that any
   pending timer callbacks do not call `haveMoreData()` on a destroyed resource.
   */
  std::shared_ptr<std::atomic<bool>> m_alive;

#if( MCP_ENABLE_AUTH )
  /** Bearer token required to use the MCP server. If empty, no authorization is required. */
  const std::string secret_token_;
#endif
};//class LlmMcpResource

#endif //LLM_MCP_RESOURCE_H
