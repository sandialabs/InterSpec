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
#include <variant>
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
  Isotopics,
  DeepResearch,
  EnergyCalibration,
  SpectrumTimeHistory,
  InSitu,          // Specialized agent for in-situ gamma spectrometry surface/soil contamination analysis
  McpServer  // Pseudo-agent for MCP-only tool filtering; not a real agent
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
    std::string ephemeral_message_txt;          // Optional per-state ephemeral message text
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

  /** Access all state definitions (state name -> definition); used for config validation. */
  const std::map<std::string, StateDefinition>& allStates() const { return m_states; }

  // Transition validation
  bool canTransitionTo( const std::string &new_state ) const;

  // State updates
  void transitionTo( const std::string &new_state );
  void reset();

  // Guidance
  std::string getPromptGuidanceForCurrentState() const;
  std::vector<std::string> getAllowedTransitions() const;
  std::string getEphemeralMessageTxtForCurrentState() const;

  /** Generate a compact summary of all states and transitions for LLM context.

   Produces a markdown-formatted overview of the full state machine graph,
   auto-generated from the XML definitions so it stays in sync automatically.
   */
  std::string getFullStateMapSummary() const;
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
    static const int sm_xmlSerializationVersion = 1;

    /** Reasoning effort level for OpenAI-style reasoning models */
    enum class ReasoningEffort
    {
      low,
      medium,
      high
    };//enum class ReasoningEffort

    /** How verbose/explicit the rendered instruction text (system prompts, tool descriptions,
     state-machine guidance) should be for a given model.  Weaker models generally benefit from more
     explicit, step-by-step instructions (Verbose); frontier models do well with terse prompts.

     Purely an input to the Inja instruction-template rendering (see LlmPromptTemplate); it is never
     sent on the wire.  Exposed to templates as the `verbosity` variable ("terse"/"normal"/"verbose").
     */
    enum class InstructionVerbosity
    {
      Terse,
      Normal,
      Verbose
    };//enum class InstructionVerbosity

    /** The wire-format / protocol an API provider speaks.

     - OpenAiChat: the legacy OpenAI Chat Completions API (`/v1/chat/completions`), and the many
       OpenAI-compatible providers (OpenRouter, Google OpenAI-compat, Ask Sage, local servers, ...).
     - OpenAiResponses: the newer OpenAI Responses API (`/v1/responses`).
     - Anthropic: the native Anthropic Messages API (`/v1/messages`).
     */
    enum class ApiFormat
    {
      OpenAiChat,
      OpenAiResponses,
      Anthropic
    };//enum class ApiFormat

    /** Map an `apiFormat` attribute string to an ApiFormat; throws std::runtime_error if unknown. */
    static ApiFormat parseApiFormat( const std::string &str );

    /** The canonical attribute string for an ApiFormat (returns a static literal). */
    static const char *apiFormatToString( ApiFormat fmt );

    /** Best-effort detection of the wire format from an endpoint URL, used when a provider does
     not explicitly specify `apiFormat`.  Contains "api.anthropic.com" or path ends with
     "/messages" -> Anthropic; ends with "/responses" -> OpenAiResponses; otherwise OpenAiChat.
     */
    static ApiFormat detectApiFormat( const std::string &endpoint );

    /** Per-model configuration within an ApiProvider. */
    struct ModelInfo
    {
      std::string name;           // e.g., "gpt-5-mini"

      /** Whether this model supports image/vision input.
       When true, the `get_spectrum_image` tool will be available to LLM agents, and image content
       will be included in tool results sent to the LLM API.
       When false, the tool is hidden from agents (but still available via MCP for external clients).
       */
      bool supportsImages = false;

      /** Reasoning configuration for LLM requests.
       Can be either:
       - bool: If true, adds `reasoning: {enabled: true}` to requests (OpenRouter style)
       - ReasoningEffort: Adds `reasoning_effort: "low"|"medium"|"high"` to requests (OpenAI o1/o3 style)
       Default is false (no reasoning).
       */
      std::variant<bool, ReasoningEffort> reasoning = false;

      /** The number of reasoning tokens to use.
       For OpenAI models, it corresponds to the `"max_completion_tokens"` field; for other models
       its the `"max_tokens"` in the JSON sent to the LLM.
       A value of 0 means use the API default.
       */
      int maxTokens = 0;

      /** Maximum context window size for this model.  A value of 0 means use the API default. */
      int contextLengthLimit = 0;

      /** Optional model temperature; valid range 0.0-2.0.  Unset means use the API default. */
      std::optional<double> temperature;

      /** When true and a state machine is active in a non-final state, sets
       `tool_choice: "required"` in the API request to force the model to make
       a tool call.  Default is false (uses `tool_choice: "auto"`).
       */
      bool requireToolInStateMachine = false;

      /** How verbose the rendered instruction text should be for this model.  Default Normal.
       Only affects instruction-template rendering (LlmPromptTemplate); never sent on the wire.
       */
      InstructionVerbosity instructionVerbosity = InstructionVerbosity::Normal;
    };//struct ModelInfo

    /** A single API provider (endpoint + token) with one or more models. */
    struct ApiProvider
    {
      std::string apiEndpoint;  // e.g., "https://api.openai.com/v1/chat/completions"
      std::string bearerToken;

      /** The wire-format this provider speaks.  An empty optional means "auto-detect from the
       endpoint URL" (see `detectApiFormat`), which keeps existing configs (that predate this
       attribute) working as OpenAiChat.  Parsed from the optional `apiFormat` attribute on the
       `<ApiProvider>` XML element.  (Added as a purely additive, default-preserving field, so the
       `LlmApi` serialization version was intentionally NOT bumped.)
       */
      std::optional<ApiFormat> apiFormat;

      std::vector<ModelInfo> models;
      size_t activeModelIndex = 0;  // Index of the active model; resolved at parse time
    };//struct ApiProvider

    bool enabled = false;

    /** A HTTP enpoint to post questions to, to do deeper research.
     At the moment, this is a seperate service that is not distributed, so it likely will not be of use to anyone besides the primary InterSpec developer.

     This field will be a value like: "http://localhost:8000/query", which if you fill it out, will cause the
     DeepResearch sub-agent to have access to the `query_deep_research_endpoint` tool-call, and post a message
     formatted like:
     ```
     {"messages": [{"content": "what is the primary use of Ba133?"}], "corpora": ["example_corpra"]}
     ```

     If you leave this field blank, the DeepResearch sub-agent still remains available, but without
     the remote endpoint query tool-call.
     */
    std::string deep_research_url;

    /** System prompt used for context compaction (summarization of older conversation history).

     If empty, a default prompt is used.  This prompt is sent as the system message when the LLM
     is asked to summarize older conversation history to manage context length.
     */
    std::string compactionSystemPrompt;

    /** Optional debug logging file path or stream name.

     If empty, no debug logging will be done.
     If "stdout" or "stderr", debug messages will be logged to standard streams.
     Otherwise, treated as a file path where debug logs will be written.

     Debug logging includes detailed information about LLM interactions, state transitions,
     tool calls, and internal processing steps useful for development and troubleshooting.
     */
    std::string debug_file;

    std::vector<ApiProvider> providers;
    size_t activeProviderIndex = 0;  // Index of the active provider; resolved at parse time

    // Convenience accessors for the active provider and model
    const ApiProvider &activeProvider() const { return providers.at( activeProviderIndex ); }
    const ModelInfo &activeModel() const { return activeProvider().models.at( activeProvider().activeModelIndex ); }

    const std::string &apiEndpoint() const { return activeProvider().apiEndpoint; }
    const std::string &bearerToken() const { return activeProvider().bearerToken; }
    const std::string &model() const { return activeModel().name; }
    bool supportsImages() const { return activeModel().supportsImages; }
    const std::variant<bool, ReasoningEffort> &reasoning() const { return activeModel().reasoning; }
    int maxTokens() const { return activeModel().maxTokens; }
    int contextLengthLimit() const { return activeModel().contextLengthLimit; }
    std::optional<double> temperature() const { return activeModel().temperature; }
    bool requireToolInStateMachine() const { return activeModel().requireToolInStateMachine; }
    InstructionVerbosity instructionVerbosity() const { return activeModel().instructionVerbosity; }

    /** The resolved wire-format for the active provider: the explicitly-configured `apiFormat`,
     or, if unset, auto-detected from the active provider's endpoint URL.
     */
    ApiFormat apiFormat() const
    {
      const ApiProvider &p = activeProvider();
      return p.apiFormat.value_or( detectApiFormat( p.apiEndpoint ) );
    }
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
   1. InterSpec::writableDataDirectory() + "/" + {filename}
   2. InterSpec::staticDataDirectory() + "/" + {filename}

   Config files loaded:
   - "llm_config.xml" - API and MCP server settings
   - "llm_tools_config.xml" - Tool definitions
   - "llm_agent_{AgentName}.xml" - Per-agent configuration files (e.g., llm_agent_NuclideId.xml)

   Each agent is loaded from its own file, allowing users to override individual agents
   by placing a modified copy in the writable data directory.

   If a user config file exists, but is invalid, throws exception.

   @return LlmConfig instance with loaded settings.
   */
  static std::shared_ptr<LlmConfig> load();
  
  /** Load LLM API and MCP server configuration from XML file.

   @param llmConfigPath Path to the llm_config.xml file
   @return pair of LlmApi and McpServer settings loaded from XML
   @throws std::runtime_error if file cannot be parsed
   */
  static std::pair<LlmApi, McpServer> loadApiAndMcpConfigs( const std::string &llmConfigPath );
  
  /** Serialize a configuration to an XML string.

   Produces exactly the bytes `saveToFile` writes, including documentation comments describing each
   field (so a GUI "preview" matches the saved file).  Only `llmApi` and `mcpServer` are serialized;
   `agents` and `tools` live in their own files and are not written here.

   @param config The configuration to serialize.
   @return The XML document as a UTF-8 string.
   */
  static std::string toXmlString( const LlmConfig &config );

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
