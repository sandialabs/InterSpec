#include "InterSpec_config.h"
#include "InterSpec/LlmConversationHistory.h"

#if( USE_LLM_INTERFACE )

#include <iostream>
#include <iomanip>
#include <sstream>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

using namespace std;
using json = nlohmann::json;

LlmConversationHistory::LlmConversationHistory() = default;

void LlmConversationHistory::addUserMessage(const std::string& message, const std::string& toolCallId) {
  LlmMessage msg(LlmMessage::Type::User, message);
  if (!toolCallId.empty()) {
    msg.toolCallId = toolCallId;
  }
  m_messages.push_back(msg);
}

void LlmConversationHistory::addAssistantMessage(const std::string& message, const std::string& toolCallId) {
  LlmMessage msg(LlmMessage::Type::Assistant, message);
  if (!toolCallId.empty()) {
    msg.toolCallId = toolCallId;
  }
  m_messages.push_back(msg);
}

void LlmConversationHistory::addSystemMessage(const std::string& message, const std::string& toolCallId) {
  LlmMessage msg(LlmMessage::Type::System, message);
  if (!toolCallId.empty()) {
    msg.toolCallId = toolCallId;
  }
  m_messages.push_back(msg);
}

void LlmConversationHistory::addToolCall(const std::string& toolName, const std::string& callId, 
                                       const nlohmann::json& parameters) {
  LlmMessage msg(LlmMessage::Type::ToolCall, "");
  msg.toolName = toolName;
  msg.toolCallId = callId;
  msg.toolParameters = parameters;
  m_messages.push_back(msg);
}

void LlmConversationHistory::addToolResult(const std::string& callId, const nlohmann::json& result) {
  LlmMessage msg(LlmMessage::Type::ToolResult, result.dump());
  msg.toolCallId = callId;
  m_messages.push_back(msg);
}

const std::vector<LlmMessage>& LlmConversationHistory::getMessages() const {
  return m_messages;
}

void LlmConversationHistory::clear() {
  m_messages.clear();
}

bool LlmConversationHistory::isEmpty() const {
  return m_messages.empty();
}

size_t LlmConversationHistory::size() const {
  return m_messages.size();
}

nlohmann::json LlmConversationHistory::toApiFormat() const {
  json messages = json::array();
  
  for (const auto& msg : m_messages) {
    json apiMsg;
    
    switch (msg.type) {
      case LlmMessage::Type::System:
        apiMsg["role"] = "system";
        apiMsg["content"] = msg.content;
        break;
        
      case LlmMessage::Type::User:
        apiMsg["role"] = "user";
        apiMsg["content"] = msg.content;
        break;
        
      case LlmMessage::Type::Assistant:
        apiMsg["role"] = "assistant";
        apiMsg["content"] = msg.content;
        break;
        
      case LlmMessage::Type::ToolCall:
        {
          apiMsg["role"] = "assistant";
          apiMsg["tool_calls"] = json::array();
          json toolCall;
          toolCall["id"] = msg.toolCallId;
          toolCall["type"] = "function";
          toolCall["function"]["name"] = msg.toolName;
          toolCall["function"]["arguments"] = msg.toolParameters.dump();
          apiMsg["tool_calls"].push_back(toolCall);
          break;
        }
        
      case LlmMessage::Type::ToolResult:
        apiMsg["role"] = "tool";
        apiMsg["tool_call_id"] = msg.toolCallId;
        apiMsg["content"] = msg.content;
        break;
    }
    
    messages.push_back(apiMsg);
  }
  
  return messages;
}

void LlmConversationHistory::toXml(rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) const {
  rapidxml::xml_node<char>* historyNode = doc->allocate_node(rapidxml::node_element, "LlmHistory");
  parent->append_node(historyNode);
  
  for (const auto& msg : m_messages) {
    rapidxml::xml_node<char>* msgNode = doc->allocate_node(rapidxml::node_element, "Message");
    historyNode->append_node(msgNode);
    
    // Add type attribute
    msgNode->append_attribute(doc->allocate_attribute("type", 
                             doc->allocate_string(typeToString(msg.type).c_str())));
    
    // Add timestamp attribute
    auto timeT = chrono::system_clock::to_time_t(msg.timestamp);
    string timeStr = std::to_string(timeT);
    msgNode->append_attribute(doc->allocate_attribute("timestamp", 
                             doc->allocate_string(timeStr.c_str())));
    
    // Add content
    if (!msg.content.empty()) {
      rapidxml::xml_node<char>* contentNode = doc->allocate_node(rapidxml::node_element, "Content", 
                                              doc->allocate_string(msg.content.c_str()));
      msgNode->append_node(contentNode);
    }
    
    // Add tool-specific fields
    if (msg.type == LlmMessage::Type::ToolCall || msg.type == LlmMessage::Type::ToolResult) {
      if (!msg.toolName.empty()) {
        rapidxml::xml_node<char>* toolNameNode = doc->allocate_node(rapidxml::node_element, "ToolName", 
                                                 doc->allocate_string(msg.toolName.c_str()));
        msgNode->append_node(toolNameNode);
      }
      if (!msg.toolCallId.empty()) {
        rapidxml::xml_node<char>* callIdNode = doc->allocate_node(rapidxml::node_element, "ToolCallId", 
                                               doc->allocate_string(msg.toolCallId.c_str()));
        msgNode->append_node(callIdNode);
      }
      if (!msg.toolParameters.empty()) {
        rapidxml::xml_node<char>* paramsNode = doc->allocate_node(rapidxml::node_element, "ToolParameters", 
                                               doc->allocate_string(msg.toolParameters.dump().c_str()));
        msgNode->append_node(paramsNode);
      }
    }
  }
}

void LlmConversationHistory::fromXml(const rapidxml::xml_node<char>* node) {
  m_messages.clear();
  
  if (!node) return;
  
  for (rapidxml::xml_node<char>* msgNode = node->first_node("Message"); 
       msgNode; msgNode = msgNode->next_sibling("Message")) {
    
    LlmMessage msg(LlmMessage::Type::User, ""); // Default, will be overridden
    
    // Read type
    if (rapidxml::xml_attribute<char>* typeAttr = msgNode->first_attribute("type")) {
      msg.type = stringToType(typeAttr->value());
    }
    
    // Read timestamp
    if (rapidxml::xml_attribute<char>* timeAttr = msgNode->first_attribute("timestamp")) {
      auto timeT = static_cast<time_t>(std::stoll(timeAttr->value()));
      msg.timestamp = chrono::system_clock::from_time_t(timeT);
    }
    
    // Read content
    if (rapidxml::xml_node<char>* contentNode = msgNode->first_node("Content")) {
      msg.content = contentNode->value();
    }
    
    // Read tool fields
    if (rapidxml::xml_node<char>* toolNameNode = msgNode->first_node("ToolName")) {
      msg.toolName = toolNameNode->value();
    }
    if (rapidxml::xml_node<char>* callIdNode = msgNode->first_node("ToolCallId")) {
      msg.toolCallId = callIdNode->value();
    }
    if (rapidxml::xml_node<char>* paramsNode = msgNode->first_node("ToolParameters")) {
      try {
        msg.toolParameters = json::parse(paramsNode->value());
      } catch (const std::exception& e) {
        cout << "Failed to parse tool parameters: " << e.what() << endl;
      }
    }
    
    m_messages.push_back(msg);
  }
  
  cout << "Loaded " << m_messages.size() << " messages from XML" << endl;
}

std::string LlmConversationHistory::typeToString(LlmMessage::Type type) {
  switch (type) {
    case LlmMessage::Type::System: return "system";
    case LlmMessage::Type::User: return "user";
    case LlmMessage::Type::Assistant: return "assistant";
    case LlmMessage::Type::ToolCall: return "tool_call";
    case LlmMessage::Type::ToolResult: return "tool_result";
    default: return "unknown";
  }
}

LlmMessage::Type LlmConversationHistory::stringToType(const std::string& str) {
  if (str == "system") return LlmMessage::Type::System;
  if (str == "user") return LlmMessage::Type::User;
  if (str == "assistant") return LlmMessage::Type::Assistant;
  if (str == "tool_call") return LlmMessage::Type::ToolCall;
  if (str == "tool_result") return LlmMessage::Type::ToolResult;
  return LlmMessage::Type::User; // Default fallback
}

#endif // USE_LLM_INTERFACE