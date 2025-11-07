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

#include <vector>
#include <string>
#include <chrono>
#include <optional>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

enum class AgentType : int;
struct LlmConversationStart;

// Forward declarations
namespace rapidxml {
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}

/** Represents a response within a conversation thread.

 This can be an assistant response, tool call, or tool result that follows
 up on a conversation start.
 */
struct LlmConversationResponse {
  enum class Type {
    Assistant,   // LLM response
    ToolCall,    // LLM requesting to call a tool
    ToolResult,  // Result from tool execution
    Error        // Error message
  };

  Type type;
  std::string content;
  std::string thinkingContent;  // Raw thinking content from LLM (e.g., <think>...</think>)
  std::chrono::system_clock::time_point timestamp;

  // For tool calls/results
  std::string toolName;
  std::string invocationId;    // ID for specific tool invocation
  nlohmann::json toolParameters;

  // A pointer back to the conversation that owns this response - not currently used, but maybe useful in the future
  std::weak_ptr<LlmConversationStart> conversation;

  std::shared_ptr<LlmConversationStart> sub_agent_conversation;
  
  LlmConversationResponse(Type t, const std::string& c, const std::shared_ptr<LlmConversationStart> &convo)
    : type(t), content(c), timestamp(std::chrono::system_clock::now()), conversation(convo) {}
};

/** Represents the start of a conversation thread.
 
 This is typically a user message or system message that initiates a conversation.
 */
struct LlmConversationStart
{
  /** Marks this conversation as either the user asking the question, or as InterSpec doing its own thing - which isnt implemented yet. */
  enum class Type {
    System,      // System prompt or system-generated message
    User         // User question or request
  };
  
  Type type;
  AgentType agent_type;
  std::string content;
  std::chrono::system_clock::time_point timestamp;
  std::string conversationId;  // ID for the entire conversation thread
  
  // Token usage tracking from LLM API responses
  // Using optional<size_t> because some models/APIs don't provide token usage information,
  // making it clear when this data is unavailable rather than defaulting to 0
  std::optional<size_t> promptTokens;      // Tokens used in the input/prompt
  std::optional<size_t> completionTokens;  // Tokens generated in the response
  std::optional<size_t> totalTokens;       // Total tokens used (prompt + completion)
  
  // Nested follow-up responses (assistant responses, tool calls, tool results)
  std::vector<LlmConversationResponse> responses;
  
  /** Function called when the conversation with the LLM has ended.
   For the main agent, this will be to update the GUI.
   For sub-agents, this will be to fill-in the agent summary and send the chat history back to the LLM
   */
  std::function<void(void)> conversation_completion_handler;
  
  LlmConversationStart(Type t, const std::string& c, AgentType a )
    : type(t), agent_type(a), content(c), timestamp(std::chrono::system_clock::now()) {}
};

// Legacy alias for backward compatibility during transition
using LlmMessage = LlmConversationStart;

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
  
  /** Add a user message - this is a message that the user typed in, and will cause a new `LlmConversationStart` to be created, with all the
   LLM responses, and tool-call results being stored in that newly created object.  The newly created object will be added to `m_conversations` and returned.
   
   This function does not send anything to the LLM - it should be called when you send this new message to the LLM.
   
   Note: the returned `LlmConversationStart::conversationId` will be a newly generated conversation identifier
   */
  std::shared_ptr<LlmConversationStart> addUserMessageToMainConversation( const std::string &message  );
  

  /** Add an assistant response with thinking content as a follow-up to the last message */
  void addAssistantMessageWithThinking(const std::string &message,
                                       const std::string &thinkingContent,
                                       std::shared_ptr<LlmConversationStart> conversation );
  
  /** Add a system message */
  std::shared_ptr<LlmConversationStart> addSystemMessageToMainConversation( const std::string &message );
  
  /** Add a tool call request as a follow-up to the last message */
  void addToolCall(const std::string &toolName,
                   const std::string &invocationId,
                   const nlohmann::json &parameters,
                   const std::shared_ptr<LlmConversationStart> &convo );

  /** Add a tool call result as a follow-up to the last message */
  void addToolResult( const std::string &invocationId,
                      const nlohmann::json &result,
                      const std::shared_ptr<LlmConversationStart> &convo );
  
  /** Add an error message as a follow-up to the last message */
  void addErrorMessage( const std::string& errorMessage, const std::shared_ptr<LlmConversationStart> &convo );
  
  
  /** Add token usage to a specific conversation by conversation ID (accumulates across API calls) */
  static void addTokenUsage( std::shared_ptr<LlmConversationStart> conversation,
                     std::optional<int> promptTokens,
                     std::optional<int> completionTokens,
                     std::optional<int> totalTokens);
  
  /** Find a conversation start by conversation ID */
  std::shared_ptr<LlmConversationStart> findConversationByConversationId(const std::string& conversationId);

  /** Get all conversation starts (const) */
  const std::vector<std::shared_ptr<LlmConversationStart>>& getConversations() const;

  /** Get all conversation starts (mutable) */
  std::vector<std::shared_ptr<LlmConversationStart>>& getConversations();
  
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
  static void addConversationToLlmApiHistory( const LlmConversationStart &conv, nlohmann::json &messages );
  
  /** Serialize to XML for saving with SpecMeas */
  void toXml(rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) const;
  
  /** Load from XML when loading SpecMeas */
  void fromXml(const rapidxml::xml_node<char>* node);
  
  /** Static function to serialize conversations to XML */
  static void toXml(const std::vector<std::shared_ptr<LlmConversationStart>>& conversations,
                    rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc);

  /** Static function to deserialize conversations from XML */
  static void fromXml(const rapidxml::xml_node<char>* node,
                      std::vector<std::shared_ptr<LlmConversationStart>>& conversations);

private:
  std::vector<std::shared_ptr<LlmConversationStart>> m_conversations;
  
  /** Convert conversation start type to string for XML */
  static std::string conversationTypeToString(LlmConversationStart::Type type);
  
  /** Convert string to conversation start type from XML */
  static LlmConversationStart::Type stringToConversationType(const std::string& str);
  
  /** Convert response type to string for XML */
  static std::string responseTypeToString(LlmConversationResponse::Type type);
  
  /** Convert string to response type from XML */
  static LlmConversationResponse::Type stringToResponseType(const std::string& str);
};

#endif // LLM_CONVERSATION_HISTORY_H
