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
#include <string>
#include <functional>

#include <Wt/WResource>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) || defined(ANDROID) )
  static_assert( 0, "MCP server can not be enabled for web or mobile builds" );
#endif

// Forward declarations
class InterSpec;

class LlmConfig;
namespace LlmTools {
  class ToolRegistry;
}

namespace Wt {
  namespace Http {
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
    void handleRequest(const Wt::Http::Request& request, Wt::Http::Response& response) override;

private:
  // A struct to hold all information about a single tool, including how to execute it.
  struct Tool {
    std::string name;
    std::string description;
    nlohmann::json parameters_schema;
    // The executor takes the parameters as JSON and returns the result as JSON.
    // NOTE: For shared tool registry compatibility, the executor should be compatible with
    // std::function<nlohmann::json(const nlohmann::json&, InterSpec*)> when possible
    std::function<nlohmann::json( nlohmann::json )> executor;
  };
  
  // Method to register a new tool with the server.
  void register_tool(Tool tool);
  
  void register_default_tools();
  
#if( MCP_ENABLE_AUTH )
  // This function is only called if authentication is enabled at compile time.
  //  If invalid request, will set response status to 401, and content to a short message.
  //  @returns true if the validation passed.
  bool validate_request( const Wt::Http::Request& request, Wt::Http::Response &response );
#endif

    void handle_get_auth(const Wt::Http::Request& request, Wt::Http::Response& response);

    // Dynamically builds the tool list from the registered tools.
    void handle_list_tools(const Wt::Http::Request& request, Wt::Http::Response& response);

    // Dispatches the call to the correct registered tool's executor.
    void handle_call_tool(const Wt::Http::Request& request, Wt::Http::Response& response);
    
    // Find InterSpec instance by session ID for session-agnostic tool execution
    // If actualSessionId is not nullptr, it will be set to the actual session ID found
    InterSpec* findInterSpecBySessionId(const std::string& sessionId, std::string* actualSessionId = nullptr);

  const std::shared_ptr<const LlmConfig> llm_config_;
  const std::unique_ptr<const LlmTools::ToolRegistry> tool_registry_;
  
  /** This class adds a "userSession" field to the parameters JSON schema, so we will make a copy. */
  std::map<std::string, Tool> registered_tools_;
  
#if( MCP_ENABLE_AUTH )
  /** Bearer token required to use the MCP server.
   If empty, no authorization is required.
   */
  const std::string secret_token_;
#endif
};//class McpResource

#endif //LLM_MCP_RESOURCE_H
