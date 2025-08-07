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

LlmConversationHistory::LlmConversationHistory() 
  : m_conversations(std::make_shared<std::vector<LlmConversationStart>>()) {}

void LlmConversationHistory::addUserMessage(const std::string& message, const std::string& conversationId) {
  LlmConversationStart conv(LlmConversationStart::Type::User, message);
  if (!conversationId.empty()) {
    conv.conversationId = conversationId;
  }
  m_conversations->push_back(conv);
}

void LlmConversationHistory::addAssistantMessage(const std::string& message, const std::string& conversationId) {
  // If we have a conversation ID, try to add this as a follow-up to an existing conversation
  if (!conversationId.empty()) {
    LlmConversationStart* parentConv = findConversationByConversationId(conversationId);
    if (parentConv) {
      LlmConversationResponse response(LlmConversationResponse::Type::Assistant, message);
      parentConv->responses.push_back(response);
      return;
    }
  }
  
  // Otherwise, add as a new top-level conversation (this shouldn't happen in normal flow)
  LlmConversationStart conv(LlmConversationStart::Type::User, ""); // Empty user message as placeholder
  conv.conversationId = conversationId;
  LlmConversationResponse response(LlmConversationResponse::Type::Assistant, message);
  conv.responses.push_back(response);
  m_conversations->push_back(conv);
}

void LlmConversationHistory::addSystemMessage(const std::string& message, const std::string& conversationId) {
  LlmConversationStart conv(LlmConversationStart::Type::System, message);
  if (!conversationId.empty()) {
    conv.conversationId = conversationId;
  }
  m_conversations->push_back(conv);
}

void LlmConversationHistory::addAssistantMessageWithThinking(const std::string& message, const std::string& thinkingContent, 
                                                           const std::string& conversationId) {
  // If we have a conversation ID, try to add this as a follow-up to an existing conversation
  if (!conversationId.empty()) {
    LlmConversationStart* parentConv = findConversationByConversationId(conversationId);
    if (parentConv) {
      LlmConversationResponse response(LlmConversationResponse::Type::Assistant, message);
      response.thinkingContent = thinkingContent;
      parentConv->responses.push_back(response);
      return;
    }
  }
  
  // Otherwise, add as a new top-level conversation (this shouldn't happen in normal flow)
  LlmConversationStart conv(LlmConversationStart::Type::User, ""); // Empty user message as placeholder
  conv.conversationId = conversationId;
  LlmConversationResponse response(LlmConversationResponse::Type::Assistant, message);
  response.thinkingContent = thinkingContent;
  conv.responses.push_back(response);
  m_conversations->push_back(conv);
}

void LlmConversationHistory::addToolCall(const std::string& toolName, const std::string& conversationId, 
                                       const std::string& invocationId, const nlohmann::json& parameters) {
  // If we have a conversation ID, try to add this as a follow-up to an existing conversation
  if (!conversationId.empty()) {
    LlmConversationStart* parentConv = findConversationByConversationId(conversationId);
    if (parentConv) {
      LlmConversationResponse response(LlmConversationResponse::Type::ToolCall, "");
      response.toolName = toolName;
      response.invocationId = invocationId;
      response.toolParameters = parameters;
      parentConv->responses.push_back(response);
      return;
    }
  }
  
  // Otherwise, add as a new top-level conversation (this shouldn't happen in normal flow)
  LlmConversationStart conv(LlmConversationStart::Type::User, ""); // Empty user message as placeholder
  conv.conversationId = conversationId;
  LlmConversationResponse response(LlmConversationResponse::Type::ToolCall, "");
  response.toolName = toolName;
  response.invocationId = invocationId;
  response.toolParameters = parameters;
  conv.responses.push_back(response);
  m_conversations->push_back(conv);
}

void LlmConversationHistory::addToolResult(const std::string& conversationId, const std::string& invocationId, 
                                         const nlohmann::json& result) {
  // If we have a conversation ID, try to add this as a follow-up to an existing conversation
  if (!conversationId.empty()) {
    LlmConversationStart* parentConv = findConversationByConversationId(conversationId);
    if (parentConv) {
      LlmConversationResponse response(LlmConversationResponse::Type::ToolResult, result.dump());
      response.invocationId = invocationId;
      parentConv->responses.push_back(response);
      return;
    }
  }
  
  // Otherwise, add as a new top-level conversation (this shouldn't happen in normal flow)
  LlmConversationStart conv(LlmConversationStart::Type::User, ""); // Empty user message as placeholder
  conv.conversationId = conversationId;
  LlmConversationResponse response(LlmConversationResponse::Type::ToolResult, result.dump());
  response.invocationId = invocationId;
  conv.responses.push_back(response);
  m_conversations->push_back(conv);
}

void LlmConversationHistory::addErrorMessage(const std::string& errorMessage, const std::string& conversationId) {
  // If we have a conversation ID, try to add this as a follow-up to an existing conversation
  if (!conversationId.empty()) {
    LlmConversationStart* parentConv = findConversationByConversationId(conversationId);
    if (parentConv) {
      LlmConversationResponse response(LlmConversationResponse::Type::Error, errorMessage);
      parentConv->responses.push_back(response);
      return;
    }
  }
  
  // If no conversation ID or conversation not found, add as a new top-level conversation
  LlmConversationStart conv(LlmConversationStart::Type::User, ""); // Empty user message as placeholder
  conv.conversationId = conversationId;
  LlmConversationResponse response(LlmConversationResponse::Type::Error, errorMessage);
  conv.responses.push_back(response);
  m_conversations->push_back(conv);
}

void LlmConversationHistory::addFollowUpResponse(const std::string& conversationId, const LlmConversationResponse& response) {
  LlmConversationStart* parentConv = findConversationByConversationId(conversationId);
  if (parentConv) {
    parentConv->responses.push_back(response);
  }
}

LlmConversationStart* LlmConversationHistory::findConversationByConversationId(const std::string& conversationId) {
  for (auto& conv : *m_conversations) {
    if (conv.conversationId == conversationId) {
      return &conv;
    }
  }
  return nullptr;
}

const std::vector<LlmConversationStart>& LlmConversationHistory::getConversations() const {
  return *m_conversations;
}

void LlmConversationHistory::setConversations(std::shared_ptr<std::vector<LlmConversationStart>> conversations) {
  m_conversations = conversations;
}

std::shared_ptr<std::vector<LlmConversationStart>> LlmConversationHistory::getConversationsPtr() const {
  return m_conversations;
}

void LlmConversationHistory::clear() {
  m_conversations->clear();
}

bool LlmConversationHistory::isEmpty() const {
  return m_conversations->empty();
}

size_t LlmConversationHistory::size() const {
  return m_conversations->size();
}

nlohmann::json LlmConversationHistory::toApiFormat() const {
  json messages = json::array();
  
  for (const auto& conv : *m_conversations) {
    // Add the conversation start message
    json apiMsg;
    
    switch (conv.type) {
      case LlmConversationStart::Type::System:
        apiMsg["role"] = "system";
        apiMsg["content"] = conv.content;
        break;
        
      case LlmConversationStart::Type::User:
        apiMsg["role"] = "user";
        apiMsg["content"] = conv.content;
        break;
    }
    
    messages.push_back(apiMsg);
    
    // Add all responses in chronological order
    for (const auto& response : conv.responses) {
      json responseMsg;
      
      switch (response.type) {
        case LlmConversationResponse::Type::Assistant:
          responseMsg["role"] = "assistant";
          responseMsg["content"] = response.content;
          break;
          
        case LlmConversationResponse::Type::ToolCall:
          {
            responseMsg["role"] = "assistant";
            responseMsg["tool_calls"] = json::array();
            json toolCall;
            toolCall["id"] = conv.conversationId + ":" + response.invocationId;
            toolCall["type"] = "function";
            toolCall["function"]["name"] = response.toolName;
            toolCall["function"]["arguments"] = response.toolParameters.dump();
            responseMsg["tool_calls"].push_back(toolCall);
            break;
          }
          
        case LlmConversationResponse::Type::ToolResult:
          responseMsg["role"] = "tool";
          responseMsg["tool_call_id"] = conv.conversationId + ":" + response.invocationId;
          responseMsg["content"] = response.content;
          break;
          
        case LlmConversationResponse::Type::Error:
          responseMsg["role"] = "assistant";
          responseMsg["content"] = "Error: " + response.content;
          break;
      }
      
      messages.push_back(responseMsg);
    }
  }
  
  return messages;
}

void LlmConversationHistory::toXml(rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) const {
  toXml(*m_conversations, parent, doc);
}

void LlmConversationHistory::toXml(const std::vector<LlmConversationStart>& conversations, 
                                   rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) {
  rapidxml::xml_node<char>* historyNode = doc->allocate_node(rapidxml::node_element, "LlmHistory");
  parent->append_node(historyNode);
  
  cout << "Serializing " << conversations.size() << " conversations to XML" << endl;
  
  for (const auto& conv : conversations) {
    rapidxml::xml_node<char>* convNode = doc->allocate_node(rapidxml::node_element, "Conversation");
    historyNode->append_node(convNode);
    
    // Add type attribute
    convNode->append_attribute(doc->allocate_attribute("type", 
                             doc->allocate_string(conversationTypeToString(conv.type).c_str())));
    
    // Add timestamp attribute
    auto timeT = chrono::system_clock::to_time_t(conv.timestamp);
    string timeStr = std::to_string(timeT);
    convNode->append_attribute(doc->allocate_attribute("timestamp", 
                             doc->allocate_string(timeStr.c_str())));
    
    // Add conversation ID
    if (!conv.conversationId.empty()) {
      convNode->append_attribute(doc->allocate_attribute("conversationId", 
                               doc->allocate_string(conv.conversationId.c_str())));
    }
    
    // Add content
    if (!conv.content.empty()) {
      rapidxml::xml_node<char>* contentNode = doc->allocate_node(rapidxml::node_element, "Content", 
                                              doc->allocate_string(conv.content.c_str()));
      convNode->append_node(contentNode);
    }
    
    // Add responses
    if (!conv.responses.empty()) {
      rapidxml::xml_node<char>* responsesNode = doc->allocate_node(rapidxml::node_element, "Responses");
      convNode->append_node(responsesNode);
      
      for (const auto& response : conv.responses) {
        rapidxml::xml_node<char>* responseNode = doc->allocate_node(rapidxml::node_element, "Response");
        responsesNode->append_node(responseNode);
        
        // Add type attribute
        responseNode->append_attribute(doc->allocate_attribute("type", 
                                   doc->allocate_string(responseTypeToString(response.type).c_str())));
        
        // Add timestamp attribute
        auto responseTimeT = chrono::system_clock::to_time_t(response.timestamp);
        string responseTimeStr = std::to_string(responseTimeT);
        responseNode->append_attribute(doc->allocate_attribute("timestamp", 
                                   doc->allocate_string(responseTimeStr.c_str())));
        
        // Add content
        if (!response.content.empty()) {
          rapidxml::xml_node<char>* responseContentNode = doc->allocate_node(rapidxml::node_element, "Content", 
                                                  doc->allocate_string(response.content.c_str()));
          responseNode->append_node(responseContentNode);
        }
        
        // Add thinking content
        if (!response.thinkingContent.empty()) {
          rapidxml::xml_node<char>* thinkingContentNode = doc->allocate_node(rapidxml::node_element, "ThinkingContent", 
                                                  doc->allocate_string(response.thinkingContent.c_str()));
          responseNode->append_node(thinkingContentNode);
        }
        
        // Add tool-specific fields for responses
        if (response.type == LlmConversationResponse::Type::ToolCall || response.type == LlmConversationResponse::Type::ToolResult) {
          if (!response.toolName.empty()) {
            rapidxml::xml_node<char>* responseToolNameNode = doc->allocate_node(rapidxml::node_element, "ToolName", 
                                                     doc->allocate_string(response.toolName.c_str()));
            responseNode->append_node(responseToolNameNode);
          }
          if (!response.invocationId.empty()) {
            rapidxml::xml_node<char>* responseInvocationIdNode = doc->allocate_node(rapidxml::node_element, "InvocationId", 
                                                   doc->allocate_string(response.invocationId.c_str()));
            responseNode->append_node(responseInvocationIdNode);
          }
          if (!response.toolParameters.empty()) {
            rapidxml::xml_node<char>* responseParamsNode = doc->allocate_node(rapidxml::node_element, "ToolParameters", 
                                                   doc->allocate_string(response.toolParameters.dump().c_str()));
            responseNode->append_node(responseParamsNode);
          }
        }
      }
    }
  }
}

void LlmConversationHistory::fromXml(const rapidxml::xml_node<char>* node) {
  fromXml(node, *m_conversations);
}

void LlmConversationHistory::fromXml(const rapidxml::xml_node<char>* node, std::vector<LlmConversationStart>& conversations) {
  conversations.clear();
  
  if (!node) {
    cout << "fromXml: No node provided" << endl;
    return;
  }
  
  cout << "fromXml: Looking for conversations in node: " << node->name() << endl;
  
  int convCount = 0;
  for (rapidxml::xml_node<char>* convNode = node->first_node("Conversation"); 
       convNode; convNode = convNode->next_sibling("Conversation")) {
    convCount++;
    
    LlmConversationStart conv(LlmConversationStart::Type::User, ""); // Default, will be overridden
    
    // Read type
    if (rapidxml::xml_attribute<char>* typeAttr = convNode->first_attribute("type")) {
      conv.type = stringToConversationType(typeAttr->value());
    }
    
    // Read timestamp
    if (rapidxml::xml_attribute<char>* timeAttr = convNode->first_attribute("timestamp")) {
      auto timeT = static_cast<time_t>(std::stoll(timeAttr->value()));
      conv.timestamp = chrono::system_clock::from_time_t(timeT);
    }
    
    // Read conversation ID
    if (rapidxml::xml_attribute<char>* convIdAttr = convNode->first_attribute("conversationId")) {
      conv.conversationId = convIdAttr->value();
    }
    
    // Read content
    if (rapidxml::xml_node<char>* contentNode = convNode->first_node("Content")) {
      conv.content = contentNode->value();
    }
    
    // Read responses
    if (rapidxml::xml_node<char>* responsesNode = convNode->first_node("Responses")) {
      for (rapidxml::xml_node<char>* responseNode = responsesNode->first_node("Response"); 
           responseNode; responseNode = responseNode->next_sibling("Response")) {
        
        LlmConversationResponse response(LlmConversationResponse::Type::Assistant, ""); // Default, will be overridden
        
        // Read response type
        if (rapidxml::xml_attribute<char>* responseTypeAttr = responseNode->first_attribute("type")) {
          response.type = stringToResponseType(responseTypeAttr->value());
        }
        
        // Read response timestamp
        if (rapidxml::xml_attribute<char>* responseTimeAttr = responseNode->first_attribute("timestamp")) {
          auto responseTimeT = static_cast<time_t>(std::stoll(responseTimeAttr->value()));
          response.timestamp = chrono::system_clock::from_time_t(responseTimeT);
        }
        
        // Read response content
        if (rapidxml::xml_node<char>* responseContentNode = responseNode->first_node("Content")) {
          response.content = responseContentNode->value();
        }
        
        // Read thinking content
        if (rapidxml::xml_node<char>* thinkingContentNode = responseNode->first_node("ThinkingContent")) {
          response.thinkingContent = thinkingContentNode->value();
        }
        
        // Read response tool fields
        if (rapidxml::xml_node<char>* responseToolNameNode = responseNode->first_node("ToolName")) {
          response.toolName = responseToolNameNode->value();
        }
        if (rapidxml::xml_node<char>* responseInvocationIdNode = responseNode->first_node("InvocationId")) {
          response.invocationId = responseInvocationIdNode->value();
        }
        if (rapidxml::xml_node<char>* responseParamsNode = responseNode->first_node("ToolParameters")) {
          try {
            response.toolParameters = json::parse(responseParamsNode->value());
          } catch (const std::exception& e) {
            cout << "Failed to parse response tool parameters: " << e.what() << endl;
          }
        }
        
        conv.responses.push_back(response);
      }
    }
    
    conversations.push_back(conv);
  }
  
  cout << "fromXml: Found " << convCount << " conversation nodes, loaded " << conversations.size() << " conversations from XML" << endl;
}

std::string LlmConversationHistory::conversationTypeToString(LlmConversationStart::Type type) {
  switch (type) {
    case LlmConversationStart::Type::System: return "system";
    case LlmConversationStart::Type::User: return "user";
    default: return "unknown";
  }
}

LlmConversationStart::Type LlmConversationHistory::stringToConversationType(const std::string& str) {
  if (str == "system") return LlmConversationStart::Type::System;
  if (str == "user") return LlmConversationStart::Type::User;
  return LlmConversationStart::Type::User; // Default fallback
}

std::string LlmConversationHistory::responseTypeToString(LlmConversationResponse::Type type) {
  switch (type) {
    case LlmConversationResponse::Type::Assistant: return "assistant";
    case LlmConversationResponse::Type::ToolCall: return "tool_call";
    case LlmConversationResponse::Type::ToolResult: return "tool_result";
    case LlmConversationResponse::Type::Error: return "error";
    default: return "unknown";
  }
}

LlmConversationResponse::Type LlmConversationHistory::stringToResponseType(const std::string& str) {
  if (str == "assistant") return LlmConversationResponse::Type::Assistant;
  if (str == "tool_call") return LlmConversationResponse::Type::ToolCall;
  if (str == "tool_result") return LlmConversationResponse::Type::ToolResult;
  if (str == "error") return LlmConversationResponse::Type::Error;
  return LlmConversationResponse::Type::Assistant; // Default fallback
}

#endif // USE_LLM_INTERFACE