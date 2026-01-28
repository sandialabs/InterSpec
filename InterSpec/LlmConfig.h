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
#include <optional>

#include <nlohmann/json.hpp>

// Forward declarations
namespace rapidxml
{
  template<class Ch> class xml_node;
}//namespace rapidxml

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

/** Agent types for LLM sub-agents */
enum class AgentType : int
{
  MainAgent,
  NuclideId,
  ActivityFit,
  Isotopics
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
/** Optional state machine for guiding agent workflows.

 State machines are completely optional - agents without one continue to work normally.
 For agents with state machines, they provide:
 - Structured workflow guidance
 - Validation of state transitions (soft enforcement - warns but doesn't block)
 - Dynamic prompt injection with state-specific guidance
 */
class AgentStateMachine
{
public:
  /** Definition of a single state in the workflow */
  struct StateDefinition
  {
    std::string name;                           // State name (e.g., "ANALYZE_REQUEST")
    std::string description;                    // Human-readable description
    std::vector<std::string> allowed_transitions; // States this can transition to
    std::vector<std::string> suggested_tools;    // Tools typically used in this state
    std::string prompt_guidance;                // Guidance text for the agent
    bool is_final = false;                      // True if this is a terminal state
  };//struct StateDefinition

private:
  std::string m_initial_state;                      // Starting state name
  std::string m_current_state;                      // Current state in this conversation
  std::map<std::string, StateDefinition> m_states;  // All state definitions

public:
  AgentStateMachine();

  /** Create a copy of this state machine with its own independent state.

   The copy shares the same state definitions but has its own current state tracker.
   */
  std::shared_ptr<AgentStateMachine> copy() const;

  // Load from XML <StateMachine> node
  void fromXml( const rapidxml::xml_node<char> *state_machine_node );

  // State queries
  const std::string& getCurrentState() const { return m_current_state; }
  const std::string& getInitialState() const { return m_initial_state; }
  const StateDefinition& getStateDefinition( const std::string &state_name ) const;
  bool hasState( const std::string &state_name ) const;
  bool isFinalState( const std::string &state_name ) const;

  // Transition validation
  bool canTransitionTo( const std::string &new_state ) const;

  // State updates
  void transitionTo( const std::string &new_state );
  void reset();

  // Guidance
  std::string getPromptGuidanceForCurrentState() const;
  std::vector<std::string> getAllowedTransitions() const;
};//class AgentStateMachine


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
    std::string description;  // Description of what this agent does (used for tool invocation description)
    std::string systemPrompt; // System prompt for this agent
    std::vector<AgentType> availableForAgents; // List of agent types that can invoke this agent (empty = all agents can invoke)
    std::shared_ptr<AgentStateMachine> state_machine; // Optional state machine for workflow guidance (nullptr if not used)
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
    std::optional<double> temperature;  // Optional: valid range 0.0-2.0
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
