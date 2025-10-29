#ifndef LLM_CONFIG_H
#define LLM_CONFIG_H
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
#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

/** Agent types for LLM sub-agents */
enum class AgentType : int
{
  MainAgent,
  NuclideId,
  ActivityFit
};//enum class AgentType

/** Convert AgentType to string name */
std::string agentTypeToString( AgentType type );

/** Convert string name to AgentType (throws if invalid) */
AgentType stringToAgentType( const std::string &name );

/** Configuration settings for LLM interface and MCP server.

 This class handles loading and saving LLM configuration from XML files.
 It follows the InterSpec pattern of looking first in writableDataDirectory()
 then falling back to staticDataDirectory().
 */
class LlmConfig
{
  static const int sm_xmlSerializationVersion = 0;

public:

  /** Configuration for a specific agent (MainAgent, NuclideId, ActivityFit, etc.) */
  struct AgentConfig
  {
    static const int sm_xmlSerializationVersion = 0;

    AgentType type;           // Agent type enum
    std::string name;         // Agent name string (e.g., "MainAgent", "NuclideId")
    std::string systemPrompt; // System prompt for this agent
  };//struct AgentConfig

  /** Configuration for a tool with role-specific descriptions */
  struct ToolConfig
  {
    static const int sm_xmlSerializationVersion = 0;

    std::string name;                              // Tool name
    std::string defaultDescription;                 // Default description for all agents
    std::map<AgentType, std::string> roleDescriptions; // Role-specific descriptions (AgentType -> description)
    std::vector<AgentType> availableForAgents;     // List of agent types that can use this tool (empty = all)
    nlohmann::json parametersSchema;               // JSON schema for tool parameters (validated at parse time)
  };//struct ToolConfig

  /** LLM API settings - this is for the "LLM Assistant" tool. */
  struct LlmApi
  {
    static const int sm_xmlSerializationVersion = 0;

    bool enabled = false;
    std::string apiEndpoint;    // Example: "https://api.openai.com/v1/chat/completions"
    std::string bearerToken;    //
    std::string model;          // Example: "gpt-4"
    int maxTokens = 0;          // Example: 4000
    int contextLengthLimit = 0; // Example: 128000;
  };//struct LlmApi


  /** MCP server settings */
  struct McpServer
  {
    static const int sm_xmlSerializationVersion = 0;
    static const std::string sm_invalid_bearer_token;

    bool enabled = false;

#if( MCP_ENABLE_AUTH )
    /** The bearer token that requestors must supply.

     `sm_invalid_bearer_token` is a canary, in principle it should never be set to this value if `enabled` is true,
     but just to be sure,  then you should refuse to create a server if set to this value.
     */
    std::string bearerToken = McpServer::sm_invalid_bearer_token;
#endif
  };//struct McpServer
  
public:
  LlmApi llmApi;
  McpServer mcpServer;

  std::vector<AgentConfig> agents;  // Agent-specific configurations
  std::vector<ToolConfig> tools;    // Tool configurations with role-specific descriptions
  
  /** Load configuration from XML files with fallback logic.
   
   Tries to load from:
   1. InterSpec::writableDataDirectory() + "/" + {llm_config.xml", "llm_agents.xml", "llm_tools_config.xml"}
   2. InterSpec::staticDataDirectory() + "/" + {llm_config.xml", "llm_agents.xml", "llm_tools_config.xml"}
   
   If a user config file exists, but is invalid, throws exception.
   
   Users can specify one or none of these files.
   
   @return LlmConfig instance with loaded settings.
   */
  static std::shared_ptr<LlmConfig> load();
  
  /** Load LLM API and MCP server configuration from XML file.

   @param llmConfigPath Path to the llm_config.xml file
   @return pair of LlmApi and McpServer settings loaded from XML
   @throws std::runtime_error if file cannot be parsed
   */
  static std::pair<LlmApi, McpServer> loadApiAndMcpConfigs( const std::string &llmConfigPath );
  
  /** Save configuration to a specific XML file.

   @param config The configuration to save
   @param filename Path where to save the XML configuration file
   @return true if save was successful, false otherwise

   Suggest saving to InterSpec::writableDataDirectory() + "/llm_config.xml".
   */
  static bool saveToFile( const LlmConfig &config, const std::string &filename );

  /** Load agent configurations from a specific XML file.

   @param agentsConfigPath Path to the llm_agents.xml file
   @return vector of AgentConfig instances loaded from XML
   @throws std::runtime_error if file cannot be parsed or doesn't exist
   */
  static std::vector<AgentConfig> loadAgentsFromFile( const std::string &agentsConfigPath );

  /** Load tool configurations from a specific XML file.

   @param toolsConfigPath Path to the llm_tools_config.xml file
   @return vector of ToolConfig instances loaded from XML
   @throws std::runtime_error if file cannot be parsed or doesn't exist
   */
  static std::vector<ToolConfig> loadToolConfigsFromFile( const std::string &toolsConfigPath );
};

#endif // LLM_CONFIG_H
