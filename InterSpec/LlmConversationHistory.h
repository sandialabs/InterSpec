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

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

// Forward declarations
namespace rapidxml {
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}

/** Represents a single message in an LLM conversation.
 
 This can be a user message, assistant response, system message, or tool call/result.
 */
struct LlmMessage {
  enum class Type {
    System,      // System prompt or system-generated message
    User,        // User question or request
    Assistant,   // LLM response
    ToolCall,    // LLM requesting to call a tool
    ToolResult   // Result from tool execution
  };
  
  Type type;
  std::string content;
  std::chrono::system_clock::time_point timestamp;
  
  // For tool calls/results
  std::string toolName;
  std::string toolCallId;
  nlohmann::json toolParameters;
  
  LlmMessage(Type t, const std::string& c) 
    : type(t), content(c), timestamp(std::chrono::system_clock::now()) {}
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
  
  /** Add a user message */
  void addUserMessage(const std::string& message, const std::string& toolCallId = "");
  
  /** Add an assistant response */
  void addAssistantMessage(const std::string& message, const std::string& toolCallId = "");
  
  /** Add a system message */
  void addSystemMessage(const std::string& message, const std::string& toolCallId = "");
  
  /** Add a tool call request */
  void addToolCall(const std::string& toolName, const std::string& callId, 
                   const nlohmann::json& parameters);
  
  /** Add a tool call result */
  void addToolResult(const std::string& callId, const nlohmann::json& result);
  
  /** Get all messages */
  const std::vector<LlmMessage>& getMessages() const;
  
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
  
private:
  std::vector<LlmMessage> m_messages;
  
  /** Convert message type to string for XML */
  static std::string typeToString(LlmMessage::Type type);
  
  /** Convert string to message type from XML */
  static LlmMessage::Type stringToType(const std::string& str);
};

#endif // LLM_CONVERSATION_HISTORY_H