#ifndef LLM_CONVERSATION_HISTORY_H
#define LLM_CONVERSATION_HISTORY_H
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
#include <vector>
#include <string>
#include <chrono>
#include <optional>

#include <Wt/WSignal>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

enum class AgentType : int;
struct LlmInteraction;

// Forward declarations
namespace rapidxml {
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}

/** Represents a single tool call request and its result.

 When the LLM requests multiple tool calls in a single response, each tool
 call is stored as an LlmToolCall within the toolCalls vector.
 */
struct LlmToolCall
{
  enum class CallStatus {
    Pending, Success, Error
  };
  
  CallStatus status;
  std::string toolName;        // Name of the tool to call
  std::string invocationId;    // Unique ID for this tool invocation
  nlohmann::json toolParameters; // Parameters to pass to the tool
  std::string content;         // Result content from executing the tool
  std::optional<std::chrono::milliseconds> executionDuration; // How long the tool took to execute

  // For sub-agent tool calls (invoke_*), this holds the sub-agent's conversation
  std::shared_ptr<LlmInteraction> sub_agent_conversation;

  LlmToolCall( const std::string &name, const std::string &id, const nlohmann::json &params )
    : status(CallStatus::Pending), toolName(name), invocationId(id), toolParameters(params)
  {
  }
};//struct LlmToolCall


/** Base class representing a turn in an LLM interaction.

 This is a linguistics "turn" - a single exchange unit in the interaction.
 Derived classes represent specific types of turns: final assistant (LLM) responses, tool requests, tool results, and errors.
 */
class LlmInteractionTurn
{
public:
  /** Type enum kept for backward compatibility and type identification */
  enum class Type
  {
    InitialRequest,   // Initial user or system message that starts a conversation
    FinalLlmResponse, // LLM response
    ToolCall,         // LLM requesting to call a tool
    ToolResult,       // Result from tool execution
    Error,            // Error message
    AutoReply         // Automatic prompt to continue when LLM has reasoning but no content/tool calls
  };

protected:
  Type m_type;
  std::chrono::system_clock::time_point m_timestamp;
  std::weak_ptr<LlmInteraction> m_conversation;
  std::string m_thinkingContent;  // Raw thinking content from LLM (e.g., <think>...</think>)
  std::optional<std::chrono::milliseconds> m_callDuration;  // Duration of API call or tool execution
  std::string m_rawContent;  // Raw content: JSON sent to LLM for tool results, or JSON received from LLM
  bool m_exclude_from_history;  // If true, exclude this turn from conversation history sent to LLM

  LlmInteractionTurn( Type t, const std::shared_ptr<LlmInteraction> &convo )
    : m_type(t), m_timestamp(std::chrono::system_clock::now()), m_conversation(convo), m_exclude_from_history(false)
  {
  }

public:
  virtual ~LlmInteractionTurn() = default;

  Type type() const { return m_type; }
  const std::chrono::system_clock::time_point& timestamp() const { return m_timestamp; }
  std::weak_ptr<LlmInteraction> conversation() const { return m_conversation; }
  const std::string& thinkingContent() const { return m_thinkingContent; }
  void setThinkingContent( const std::string &thinking ) { m_thinkingContent = thinking; }
  const std::optional<std::chrono::milliseconds>& callDuration() const { return m_callDuration; }
  void setCallDuration( std::chrono::milliseconds duration ) { m_callDuration = duration; }
  const std::string& rawContent() const { return m_rawContent; }
  void setRawContent( const std::string &raw ) { m_rawContent = raw; }
  bool excludeFromHistory() const { return m_exclude_from_history; }
  void setExcludeFromHistory( bool exclude ) { m_exclude_from_history = exclude; }
};//class LlmInteractionTurn


/** Represents a final response from the LLM (assistant message) - there will be no further messages/exchanges/tool-calls for thei `LlmInteraction`. */
class LlmInteractionFinalResponse : public LlmInteractionTurn
{
protected:
  std::string m_content;

public:
  LlmInteractionFinalResponse( const std::string &content, const std::shared_ptr<LlmInteraction> &convo )
    : LlmInteractionTurn(Type::FinalLlmResponse, convo), m_content(content)
  {
  }

  const std::string& content() const { return m_content; }
  void setContent( const std::string &content ) { m_content = content; }
};//class LlmInteractionFinalResponse


/** Represents a request from the LLM to call one or more tools. */
class LlmToolRequest : public LlmInteractionTurn
{
protected:
  std::vector<LlmToolCall> m_toolCalls;

public:
  LlmToolRequest( const std::shared_ptr<LlmInteraction> &convo )
    : LlmInteractionTurn(Type::ToolCall, convo)
  {
  }

  const std::vector<LlmToolCall>& toolCalls() const { return m_toolCalls; }
  std::vector<LlmToolCall>& toolCalls() { return m_toolCalls; }
  void setToolCalls( std::vector<LlmToolCall> &&calls ) { m_toolCalls = std::move(calls); }
};//class LlmToolRequest


/** Represents the results of tool execution(s). */
class LlmToolResults : public LlmInteractionTurn
{
protected:
  std::vector<LlmToolCall> m_toolCalls;

public:
  LlmToolResults( const std::shared_ptr<LlmInteraction> &convo )
    : LlmInteractionTurn(Type::ToolResult, convo)
  {
  }

  const std::vector<LlmToolCall>& toolCalls() const { return m_toolCalls; }
  std::vector<LlmToolCall>& toolCalls() { return m_toolCalls; }
  void setToolCalls( std::vector<LlmToolCall> &&calls ) { m_toolCalls = std::move(calls); }
};//class LlmToolResults


/** Represents an error that occurred during conversation processing. */
class LlmInteractionError : public LlmInteractionTurn
{
public:
  enum class ErrorType
  {
    Unknown,           // Parsing errors, unexpected errors
    Timeout,          // timeout_error from JavaScript
    Network,          // network_error from JavaScript
    LlmApi,           // LLM API errors (rate limit, invalid request, etc)
    JsonParse         // JSON parsing failures
  };

protected:
  std::string m_errorMessage;
  ErrorType m_errorType;
  bool m_retryAttempted;  // True if JavaScript automatic retry was attempted

public:
  LlmInteractionError( const std::string &errorMsg,
                       const std::shared_ptr<LlmInteraction> &convo,
                       ErrorType type = ErrorType::Unknown )
    : LlmInteractionTurn(Type::Error, convo),
      m_errorMessage(errorMsg),
      m_errorType(type),
      m_retryAttempted(false)
  {
  }

  const std::string& errorMessage() const { return m_errorMessage; }
  void setErrorMessage( const std::string &msg ) { m_errorMessage = msg; }

  ErrorType errorType() const { return m_errorType; }
  void setErrorType( ErrorType type ) { m_errorType = type; }

  bool retryAttempted() const { return m_retryAttempted; }
  void setRetryAttempted( bool attempted ) { m_retryAttempted = attempted; }

  /** Returns true if this error supports automatic retry (timeout/network errors) */
  bool supportsAutomaticRetry() const
  {
    return (m_errorType == ErrorType::Timeout || m_errorType == ErrorType::Network);
  }
};//class LlmInteractionError


/** Represents an automatic prompt sent to the LLM to continue when it has reasoning but no content or tool calls.

 This is used when the LLM returns only reasoning/thinking content but no actual response or tool calls,
 suggesting it needs prompting to continue with its planned actions.
 */
class LlmInteractionAutoReply : public LlmInteractionTurn
{
protected:
  std::string m_content;

public:
  LlmInteractionAutoReply( const std::string &content, const std::shared_ptr<LlmInteraction> &convo )
    : LlmInteractionTurn(Type::AutoReply, convo), m_content(content)
  {
  }

  const std::string& content() const { return m_content; }
  void setContent( const std::string &content ) { m_content = content; }
};//class LlmInteractionAutoReply


/** Represents the initial request that starts a conversation (user or system message).

 This replaces the old approach of storing the initial request in LlmInteraction::content.
 Now the initial request is stored as the first turn in the responses array.
 */
class LlmInteractionInitialRequest : public LlmInteractionTurn
{
public:
  enum class RequestType
  {
    System,  // System prompt or system-generated message
    User     // User question or request
  };

protected:
  RequestType m_requestType;
  std::string m_content;

public:
  LlmInteractionInitialRequest( RequestType reqType, const std::string &content,
                                const std::shared_ptr<LlmInteraction> &convo )
    : LlmInteractionTurn(Type::InitialRequest, convo), m_requestType(reqType), m_content(content)
  {
  }

  RequestType requestType() const { return m_requestType; }
  void setRequestType( RequestType type ) { m_requestType = type; }

  const std::string& content() const { return m_content; }
  void setContent( const std::string &content ) { m_content = content; }
};//class LlmInteractionInitialRequest


/** Represents a user asking the LLM a question or to perform a task.

 This is typically a user message or system message that initiates a conversation.
 */
struct LlmInteraction
{
  /** Marks this conversation as either the user asking the question, or as InterSpec doing its own thing - which isnt implemented yet. */
  enum class Type {
    System,      // System prompt or system-generated message
    User         // User question or request
  };
  
  Type type;
  AgentType agent_type;

  /** The time this conversation was started. */
  std::chrono::system_clock::time_point timestamp;
  
  /** ID for the entire conversation thread */
  std::string conversationId;

  /** Time when this conversation finished (completed or errored). Empty if still in progress. */
  std::optional<std::chrono::system_clock::time_point> finishTime;

  // Token usage tracking from LLM API responses
  // Using optional<size_t> because some models/APIs don't provide token usage information,
  // making it clear when this data is unavailable rather than defaulting to 0
  std::optional<size_t> promptTokens;      // Tokens used in the input/prompt
  std::optional<size_t> completionTokens;  // Tokens generated in the response
  std::optional<size_t> totalTokens;       // Total tokens used (prompt + completion)
  
  /** The conversation turns sent to, or recieved from the LLM.
   
   The first entry is the initial prompt sent to the LLM to start the conversation, with following entries being
   assistant responses, tool calls, and tool results.
   */
  std::vector<std::shared_ptr<LlmInteractionTurn>> responses;

  /** Function called when the conversation with the LLM has ended.
   For the main agent, this will be to update the GUI.
   For sub-agents, this will be to fill-in the agent summary and send the chat history back to the LLM
   */
  std::function<void(void)> conversation_completion_handler;

  /** Signal emitted when a new response is added to this conversation.
   Emits the newly added response.
   */
  Wt::Signal<std::shared_ptr<LlmInteractionTurn>> responseAdded;

  /** Signal emitted when a sub-agent conversation completes.
   Emits the parent response (containing the sub-agent) and the sub-agent conversation.
   This is emitted on the parent conversation, not the sub-agent conversation.
   */
  Wt::Signal<std::shared_ptr<LlmInteractionTurn>, std::shared_ptr<LlmInteraction>> subAgentFinished;

  /** Signal emitted when this conversation finishes (no more tool calls to process).
   Can be used by the GUI to re-enable input.
   */
  Wt::Signal<> conversationFinished;

  /** Factory method to create a new LlmInteraction with an initial message.

   This ensures the InitialRequest turn is properly added to the responses array.

   @param t The conversation type (User or System)
   @param initialMessage The initial message content
   @param a The agent type
   @return A shared_ptr to the newly created LlmInteraction
   */
  static std::shared_ptr<LlmInteraction> create( Type t, const std::string& initialMessage, AgentType a );

  /** Factory method for creating an empty LlmInteraction for deserialization.

   This is used when loading conversations from XML where fields will be populated separately.
   Should only be used by LlmConversationHistory deserialization code.

   @return A shared_ptr to an empty LlmInteraction with default values
   */
  static std::shared_ptr<LlmInteraction> createEmpty();

private:
  LlmInteraction(Type t, AgentType a )
  : type(t), agent_type(a), timestamp(std::chrono::system_clock::now()), conversationId{}
  {
  }

  friend class LlmConversationHistory; // Allow LlmConversationHistory to call createEmpty()
};



/** Manages conversation history for LLM interactions.
 
 This class handles:
 - Storing messages in chronological order
 - Serializing to/from XML for persistence with SpecMeas
 - Context length management (future)
 - Conversation summarization (future)
 */
class LlmConversationHistory {
public:
  /** Construct empty conversation history */
  LlmConversationHistory();
  
  /** Add a user message - this is a message that the user typed in, and will cause a new `LlmInteraction` to be created, with all the
   LLM responses, and tool-call results being stored in that newly created object.  The newly created object will be added to `m_conversations` and returned.
   
   This function does not send anything to the LLM - it should be called when you send this new message to the LLM.
   
   Note: the returned `LlmInteraction::conversationId` will be a newly generated conversation identifier
   */
  std::shared_ptr<LlmInteraction> addUserMessageToMainConversation( const std::string &message );
  

  /** Add an assistant response with thinking content as a follow-up to the last message */
  std::shared_ptr<LlmInteractionFinalResponse> addAssistantMessageWithThinking(const std::string &message,
                                       const std::string &thinkingContent,
                                       const std::string &rawContent,
                                       std::shared_ptr<LlmInteraction> conversation );
  
  /** Add a system message */
  std::shared_ptr<LlmInteraction> addSystemMessageToMainConversation( const std::string &message );
  
  /** Add tool call requests (potentially multiple) as a follow-up to the last message.

   This creates a single LlmToolRequest containing all the tool calls.
   */
  std::shared_ptr<LlmToolRequest> addToolCalls( std::vector<LlmToolCall> &&toolCalls,
                                               const std::string &rawResponseContent,
                     const std::shared_ptr<LlmInteraction> &convo );

  /** Add tool call results (potentially multiple) as a follow-up to the last message.

   This creates a single LlmToolResults containing all the results.
   The jsonSentToLlm parameter contains the JSON that was sent back to the LLM.
   */
  std::shared_ptr<LlmToolResults> addToolResults( std::vector<LlmToolCall> &&toolResults,
                       const std::string &jsonSentToLlm,
                       const std::shared_ptr<LlmInteraction> &convo );
  
  /** Add an error message as a follow-up to the last message */
  std::shared_ptr<LlmInteractionError> addErrorMessage( const std::string &errorMessage,
                       const std::string &rawResponseContent,
                       const std::shared_ptr<LlmInteraction> &convo,
                       LlmInteractionError::ErrorType errorType = LlmInteractionError::ErrorType::Unknown );

  /** Add an automatic prompt message to continue the conversation */
  std::shared_ptr<LlmInteractionAutoReply> addAutoReplyMessage( const std::string &promptMessage,
                       const std::shared_ptr<LlmInteraction> &convo );


  /** Add token usage to a specific conversation by conversation ID (accumulates across API calls) */
  static void addTokenUsage( std::shared_ptr<LlmInteraction> conversation,
                     std::optional<int> promptTokens,
                     std::optional<int> completionTokens,
                     std::optional<int> totalTokens);
  
  /** Find a conversation start by conversation ID */
  std::shared_ptr<LlmInteraction> findConversationByConversationId(const std::string& conversationId);

  /** Get all conversation starts (const) */
  const std::vector<std::shared_ptr<LlmInteraction>>& getConversations() const;

  /** Get all conversation starts (mutable) */
  std::vector<std::shared_ptr<LlmInteraction>>& getConversations();
  
  /** Clear all messages */
  void clear();
  
  /** Check if there are any messages */
  bool isEmpty() const;
  
  /** Get number of messages */
  size_t size() const;
  
  /** Convert messages to OpenAI API format for sending to LLM */
  nlohmann::json toApiFormat() const;
  
  /** Adds the given conversation to the messages array - to send to the LLM.
   @param conv The conversation you want to add to the messages
   @param messages The "messages" array that you will sent to the LLM.  Just be a `json::array()`, or exception will be thrown..
   */
  static void addConversationToLlmApiHistory( const LlmInteraction &conv, nlohmann::json &messages );
  
  /** Serialize to XML for saving with SpecMeas */
  void toXml(rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) const;
  
  /** Load from XML when loading SpecMeas */
  void fromXml(const rapidxml::xml_node<char>* node);
  
  /** Static function to serialize conversations to XML */
  static void toXml(const std::vector<std::shared_ptr<LlmInteraction>>& conversations,
                    rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc);

  /** Static function to deserialize conversations from XML */
  static void fromXml(const rapidxml::xml_node<char>* node,
                      std::vector<std::shared_ptr<LlmInteraction>>& conversations);

private:
  std::vector<std::shared_ptr<LlmInteraction>> m_conversations;
  
  /** Convert conversation start type to string for XML */
  static std::string conversationTypeToString(LlmInteraction::Type type);
  
  /** Convert string to conversation start type from XML */
  static LlmInteraction::Type stringToConversationType(const std::string& str);
  
  /** Convert response type to string for XML */
  static std::string responseTypeToString(LlmInteractionTurn::Type type);

  /** Convert string to response type from XML */
  static LlmInteractionTurn::Type stringToResponseType(const std::string& str);

  /** Convert tool call status to string for XML */
  static std::string callStatusToString(LlmToolCall::CallStatus status);

  /** Convert string to tool call status from XML */
  static LlmToolCall::CallStatus stringToCallStatus(const std::string& str);
};

#endif // LLM_CONVERSATION_HISTORY_H
