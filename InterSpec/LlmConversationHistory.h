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
  
  LlmConversationResponse(Type t, const std::string& c) 
    : type(t), content(c), timestamp(std::chrono::system_clock::now()) {}
};

/** Represents the start of a conversation thread.
 
 This is typically a user message or system message that initiates a conversation.
 */
struct LlmConversationStart {
  enum class Type {
    System,      // System prompt or system-generated message
    User         // User question or request
  };
  
  Type type;
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
  
  LlmConversationStart(Type t, const std::string& c) 
    : type(t), content(c), timestamp(std::chrono::system_clock::now()) {}
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
  
  /** Add a user message */
  void addUserMessage(const std::string& message, const std::string& conversationId = "");
  
  /** Add an assistant response as a follow-up to the last message */
  void addAssistantMessage(const std::string& message, const std::string& conversationId = "");
  
  /** Add an assistant response with thinking content as a follow-up to the last message */
  void addAssistantMessageWithThinking(const std::string& message, const std::string& thinkingContent, 
                                     const std::string& conversationId = "");
  
  /** Add a system message */
  void addSystemMessage(const std::string& message, const std::string& conversationId = "");
  
  /** Add a tool call request as a follow-up to the last message */
  void addToolCall(const std::string& toolName, const std::string& conversationId, 
                   const std::string& invocationId, const nlohmann::json& parameters);
  
  /** Add a tool call result as a follow-up to the last message */
  void addToolResult(const std::string& conversationId, const std::string& invocationId, 
                     const nlohmann::json& result);
  
  /** Add an error message as a follow-up to the last message */
  void addErrorMessage(const std::string& errorMessage, const std::string& conversationId = "");
  
  /** Add a follow-up response to a specific conversation by conversation ID */
  void addFollowUpResponse(const std::string& conversationId, const LlmConversationResponse& response);
  
  /** Add token usage to a specific conversation by conversation ID (accumulates across API calls) */
  void addTokenUsage(const std::string& conversationId, 
                     std::optional<int> promptTokens,
                     std::optional<int> completionTokens,
                     std::optional<int> totalTokens);
  
  /** Find a conversation start by conversation ID */
  LlmConversationStart* findConversationByConversationId(const std::string& conversationId);
  
  /** Get all conversation starts */
  const std::vector<LlmConversationStart>& getConversations() const;
  
  /** Set conversations from shared pointer */
  void setConversations(std::shared_ptr<std::vector<LlmConversationStart>> conversations);
  
  /** Get conversations as shared pointer */
  std::shared_ptr<std::vector<LlmConversationStart>> getConversationsPtr() const;
  
  /** Clear all messages */
  void clear();
  
  /** Check if there are any messages */
  bool isEmpty() const;
  
  /** Get number of messages */
  size_t size() const;
  
  /** Convert messages to OpenAI API format for sending to LLM */
  nlohmann::json toApiFormat() const;
  
  /** Serialize to XML for saving with SpecMeas */
  void toXml(rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) const;
  
  /** Load from XML when loading SpecMeas */
  void fromXml(const rapidxml::xml_node<char>* node);
  
  /** Static function to serialize conversations to XML */
  static void toXml(const std::vector<LlmConversationStart>& conversations, 
                    rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc);
  
  /** Static function to deserialize conversations from XML */
  static void fromXml(const rapidxml::xml_node<char>* node, 
                      std::vector<LlmConversationStart>& conversations);
  
private:
  std::shared_ptr<std::vector<LlmConversationStart>> m_conversations;
  
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
