#include "InterSpec_config.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/LlmConfig.h"
#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/LlmConversationHistory.h"


static_assert( USE_LLM_INTERFACE, "You should not be compiling this file without USE_LLM_INTERFACE enabled" );

using namespace std;
using json = nlohmann::json;

LlmConversationHistory::LlmConversationHistory()
  : m_conversations()
{
}

std::shared_ptr<LlmConversationStart> LlmConversationHistory::addUserMessageToMainConversation( const string &message )
{
  auto conv = make_shared<LlmConversationStart>( LlmConversationStart::Type::User, message, AgentType::MainAgent );
  conv->conversationId = "conv_" + std::to_string(chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now().time_since_epoch()).count());

  m_conversations.push_back(conv);
  
  return conv;
}


std::shared_ptr<LlmConversationStart> LlmConversationHistory::addSystemMessageToMainConversation( const std::string &message )
{
  auto conv = std::make_shared<LlmConversationStart>(LlmConversationStart::Type::System, message, AgentType::MainAgent );
  conv->conversationId = "conv_" + std::to_string(chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now().time_since_epoch()).count());
  
  m_conversations.push_back(conv);
  
  return conv;
}

void LlmConversationHistory::addAssistantMessageWithThinking(const std::string& message,
                                                             const std::string& thinkingContent,
                                                             std::shared_ptr<LlmConversationStart> conversation )
{
  assert( conversation );
  if( !conversation )
  {
    cerr << "LlmConversationHistory::addAssistantMessageWithThinking: null conversation - not adding LLM response" << endl;
    return;
  }
  
  // Add message as a follow-up to an existing conversation
  LlmConversationResponse response( LlmConversationResponse::Type::Assistant, message, conversation );
  response.thinkingContent = thinkingContent;
  conversation->responses.push_back(response);
}//void addAssistantMessageWithThinking(...)


void LlmConversationHistory::addToolCall(const std::string& toolName,
                                         const std::string& invocationId,
                                         const nlohmann::json& parameters,
                                         const std::shared_ptr<LlmConversationStart> &convo )
{
  // If we have a conversation ID, try to add this as a follow-up to an existing conversation
  assert( convo );
  if( !convo )
  {
    cerr << "LlmConversationHistory::addToolCall: invalid LlmConversationStart passed in." << endl;
    throw runtime_error( "LlmConversationHistory::addToolCall: invalid LlmConversationStart passed in." );
    return;
  }
  
  LlmConversationResponse response( LlmConversationResponse::Type::ToolCall, "", convo );
  response.toolName = toolName;
  response.invocationId = invocationId;
  response.toolParameters = parameters;
  convo->responses.push_back(response);
}

void LlmConversationHistory::addToolResult( const std::string &invocationId,
                                           const nlohmann::json &result,
                                           const std::shared_ptr<LlmConversationStart> &convo )
{
  assert( convo );
  if( !convo )
  {
    cerr << "LlmConversationHistory::addToolResult - got null conversation pointer" << endl;
    return;
  }
  
  // Add this as a follow-up to an existing conversation
  LlmConversationResponse response(LlmConversationResponse::Type::ToolResult, result.dump(), convo );
  response.invocationId = invocationId;
  convo->responses.push_back(response);
}

void LlmConversationHistory::addErrorMessage( const std::string &errorMessage,
                                             const std::shared_ptr<LlmConversationStart> &convo )
{
  assert( convo );
  if( !convo )
  {
    cerr << "LlmConversationHistory::addErrorMessage - got null conversation pointer: errorMessage='" << errorMessage << "'" << endl;
    return;
  }
  
  // Add this as a follow-up to an existing conversation
  convo->responses.emplace_back(LlmConversationResponse::Type::Error, errorMessage, convo);
}


void LlmConversationHistory::addTokenUsage( std::shared_ptr<LlmConversationStart> conversation,
                                           std::optional<int> promptTokens,
                                           std::optional<int> completionTokens,
                                           std::optional<int> totalTokens )
{
  assert( conversation );
  if( !conversation )
    return;
  
  // Accumulate token usage across API calls within this conversation
  if (promptTokens.has_value() && (promptTokens.value() > 0) )
  {
    if( conversation->promptTokens.has_value() )
      conversation->promptTokens = conversation->promptTokens.value() + promptTokens.value();
    else
      conversation->promptTokens = static_cast<size_t>( promptTokens.value() );
  }
  
  if( completionTokens.has_value() && (completionTokens.value() > 0) )
  {
    if (conversation->completionTokens.has_value())
      conversation->completionTokens = conversation->completionTokens.value() + completionTokens.value();
    else
      conversation->completionTokens = static_cast<size_t>( completionTokens.value() );
  }
  
  if( totalTokens.has_value() && (totalTokens.value() > 0) )
  {
    if (conversation->totalTokens.has_value())
      conversation->totalTokens = conversation->totalTokens.value() + totalTokens.value();
    else
      conversation->totalTokens = static_cast<size_t>( totalTokens.value() );
  }
}//void addTokenUsage(...)


std::shared_ptr<LlmConversationStart> LlmConversationHistory::findConversationByConversationId(const std::string& conversationId)
{
  for( shared_ptr<LlmConversationStart> &conv : m_conversations)
  {
    if (conv->conversationId == conversationId)
      return conv;
  }
  return nullptr;
}

const std::vector<std::shared_ptr<LlmConversationStart>>& LlmConversationHistory::getConversations() const {
  return m_conversations;
}

std::vector<std::shared_ptr<LlmConversationStart>>& LlmConversationHistory::getConversations() {
  return m_conversations;
}

void LlmConversationHistory::clear() {
  m_conversations.clear();
}

bool LlmConversationHistory::isEmpty() const {
  return m_conversations.empty();
}

size_t LlmConversationHistory::size() const {
  return m_conversations.size();
}


void LlmConversationHistory::addConversationToLlmApiHistory( const LlmConversationStart &conv, nlohmann::json &messages )
{
  assert( messages.is_array() );
  if( !messages.is_array() )
    throw logic_error( "addConversationToLlmApiHistory: messages must be an array." );
  
  json apiMsg;
  
  switch( conv.type )
  {
    case LlmConversationStart::Type::System:
      apiMsg["role"] = "system";
      apiMsg["content"] = conv.content;
      break;
      
    case LlmConversationStart::Type::User:
      apiMsg["role"] = "user";
      apiMsg["content"] = conv.content;
      break;
  }//switch( conv.type )
  
  messages.push_back(apiMsg);
  
  // Add all responses in chronological order
  for( const LlmConversationResponse &response : conv.responses )
  {
    json responseMsg;
    
    switch( response.type )
    {
      case LlmConversationResponse::Type::Assistant:
        responseMsg["role"] = "assistant";
        responseMsg["content"] = response.content;
        break;
        
      case LlmConversationResponse::Type::ToolCall:
      {
        responseMsg["role"] = "assistant";
        responseMsg["tool_calls"] = json::array();
        json toolCall;
        // Use just the invocationId to keep within OpenAI's 40-character limit
        toolCall["id"] = response.invocationId;
        toolCall["type"] = "function";
        toolCall["function"]["name"] = response.toolName;
        toolCall["function"]["arguments"] = response.toolParameters.dump();
        responseMsg["tool_calls"].push_back(toolCall);
        break;
      }
        
      case LlmConversationResponse::Type::ToolResult:
        responseMsg["role"] = "tool";
        responseMsg["tool_call_id"] = response.invocationId;
        responseMsg["content"] = response.content;
        break;
        
      case LlmConversationResponse::Type::Error:
        responseMsg["role"] = "assistant";
        responseMsg["content"] = "Error: " + response.content;
        break;
    }//switch( response.type )
    
    messages.push_back(responseMsg);
  }//for( const LlmConversationResponse &response : conv.responses )
}//nlohmann::json addConversationToLlmApiHistory(...)


nlohmann::json LlmConversationHistory::toApiFormat() const {
  json messages = json::array();

  for (const shared_ptr<LlmConversationStart> &conv : m_conversations)
  {
    assert( conv );
    addConversationToLlmApiHistory( *conv, messages );
  }
  
  return messages;
}

void LlmConversationHistory::toXml(rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) const {
  toXml(m_conversations, parent, doc);
}

void LlmConversationHistory::toXml( const vector<shared_ptr<LlmConversationStart>> &conversations,
                                   rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc)
{
  rapidxml::xml_node<char>* historyNode = doc->allocate_node(rapidxml::node_element, "LlmHistory");
  parent->append_node(historyNode);

  cout << "Serializing " << conversations.size() << " conversations to XML" << endl;

  for( const shared_ptr<LlmConversationStart> &conv : conversations )
  {
    rapidxml::xml_node<char>* convNode = doc->allocate_node(rapidxml::node_element, "Conversation");
    historyNode->append_node(convNode);

    // Add type attribute
    const string convoType = conversationTypeToString(conv->type);
    XmlUtils::append_attrib(convNode, "type", convoType );
    
    // Add agent type
    XmlUtils::append_attrib( convNode, "agentType", agentTypeToString(conv->agent_type) );
    
    // Add timestamp attribute
    auto timeT = chrono::system_clock::to_time_t(conv->timestamp);
    string timeStr = std::to_string(timeT);
    XmlUtils::append_attrib(convNode, "timestamp", timeStr );
    

    // Add conversation ID
    if( !conv->conversationId.empty() )
      XmlUtils::append_attrib(convNode, "conversationId", conv->conversationId );

    // Add content
    if( !conv->content.empty() )
      XmlUtils::append_string_node( convNode, "Content", conv->content );

    // Add responses
    rapidxml::xml_node<char> *responsesNode = nullptr; //Will create at first message
    for( const LlmConversationResponse &response : conv->responses )
    {
      if( !responsesNode ) // create "Responses" at the first message, so we dont create it if responses is empty
      {
        responsesNode = doc->allocate_node(rapidxml::node_element, "Responses");
        convNode->append_node(responsesNode);
      }
      
      rapidxml::xml_node<char>* responseNode = doc->allocate_node(rapidxml::node_element, "Response");
      responsesNode->append_node(responseNode);
      
      // Add type attribute
      XmlUtils::append_attrib(responseNode, "type", responseTypeToString(response.type) );
      
      // Add timestamp attribute
      const auto responseTimeT = chrono::system_clock::to_time_t(response.timestamp);
      const string responseTimeStr = std::to_string(responseTimeT);
      XmlUtils::append_attrib(responseNode, "timestamp", responseTimeStr );
      
      // Add content
      if( !response.content.empty() )
        XmlUtils::append_string_node( responseNode, "Content", response.content );
      
      // Add thinking content
      if( !response.thinkingContent.empty() )
        XmlUtils::append_string_node( responseNode, "ThinkingContent", response.thinkingContent );
      
      // Add tool-specific fields for responses
      if( (response.type == LlmConversationResponse::Type::ToolCall)
         || (response.type == LlmConversationResponse::Type::ToolResult) )
      {
        if( !response.toolName.empty() )
          XmlUtils::append_string_node( responseNode, "ToolName", response.toolName );
        
        if( !response.invocationId.empty() )
          XmlUtils::append_string_node( responseNode, "InvocationId", response.invocationId );
        
        if( !response.toolParameters.empty() )
          XmlUtils::append_string_node( responseNode, "ToolParameters", response.toolParameters.dump() );
      }//if( response is ToolCall or ToolResult )
    }//for( const LlmConversationResponse &response : conv->responses )
  }//for( const shared_ptr<LlmConversationStart> &conv : conversations )
}//void toXml( const vector<shared_ptr<LlmConversationStart>> &conversations ... )

void LlmConversationHistory::fromXml( const rapidxml::xml_node<char> *node )
{
  fromXml(node, m_conversations);
}

void LlmConversationHistory::fromXml( const rapidxml::xml_node<char> *node, std::vector<std::shared_ptr<LlmConversationStart>> &conversations )
{
  // TODO: this function needs a bit of refactoring (LlmConversationStart should get its own fromXml function) and we should be more exacting in requing the various fields and throw exceptions if they are not there, or invalid
  conversations.clear();

  if( !node )
  {
    cout << "fromXml: No node provided" << endl;
    return;
  }

  cout << "fromXml: Looking for conversations in node: " << node->name() << endl;

  int convCount = 0;
  XML_FOREACH_CHILD(convNode, node, "Conversation")
  {
    convCount++;

    auto conv = std::make_shared<LlmConversationStart>( LlmConversationStart::Type::User, "", AgentType::MainAgent ); // Default, will be overridden

    // Read type
    if( rapidxml::xml_attribute<char>* typeAttr = XML_FIRST_ATTRIB(convNode, "type") )
      conv->type = stringToConversationType( SpecUtils::xml_value_str(typeAttr) );

    // Read agent type
    if( rapidxml::xml_attribute<char> *agentTypeAttrib = XML_FIRST_ATTRIB(convNode,"agentType") )
    {
      const string agent_type = SpecUtils::xml_value_str(agentTypeAttrib);
      conv->agent_type = stringToAgentType(agent_type);
    }
    
    // Read timestamp
    if( rapidxml::xml_attribute<char> *timeAttr = XML_FIRST_ATTRIB(convNode, "timestamp") )
    {
      const string time_str = SpecUtils::xml_value_str(timeAttr);
      auto timeT = static_cast<time_t>( std::stoll( time_str.c_str() ) );
      conv->timestamp = chrono::system_clock::from_time_t(timeT);
    }

    // Read conversation ID
    if( rapidxml::xml_attribute<char>* convIdAttr = XML_FIRST_ATTRIB(convNode, "conversationId") )
      conv->conversationId = SpecUtils::xml_value_str(convIdAttr);

    // Read content
    if( rapidxml::xml_node<char> *contentNode = XML_FIRST_NODE(convNode,"Content") )
      conv->content = SpecUtils::xml_value_str(contentNode);
    
    // Read responses
    if( rapidxml::xml_node<char> *responsesNode = XML_FIRST_NODE(convNode, "Responses") )
    {
      XML_FOREACH_CHILD( responseNode, responsesNode, "Response" )
      {
        LlmConversationResponse response( LlmConversationResponse::Type::Assistant, "", conv ); // Default, will be overridden
        
        // Read response type
        if( rapidxml::xml_attribute<char>* responseTypeAttr = XML_FIRST_ATTRIB(responseNode,"type") )
          response.type = stringToResponseType( SpecUtils::xml_value_str(responseTypeAttr) );
        
        // Read response timestamp
        if( rapidxml::xml_attribute<char> *responseTimeAttr = XML_FIRST_ATTRIB(responseNode,"timestamp") )
        {
          const string time_str = SpecUtils::xml_value_str(responseTimeAttr);
          auto responseTimeT = static_cast<time_t>( std::stoll(time_str.c_str()) );
          response.timestamp = chrono::system_clock::from_time_t(responseTimeT);
        }
        
        // Read response content
        if( rapidxml::xml_node<char>* responseContentNode = XML_FIRST_NODE(responseNode,"Content") )
          response.content = SpecUtils::xml_value_str(responseContentNode );
        
        // Read thinking content
        if( rapidxml::xml_node<char>* thinkingContentNode = XML_FIRST_NODE(responseNode,"ThinkingContent") )
          response.thinkingContent = SpecUtils::xml_value_str(thinkingContentNode);

        // Read response tool fields
        if( rapidxml::xml_node<char>* responseToolNameNode = XML_FIRST_NODE(responseNode, "ToolName") )
          response.toolName = SpecUtils::xml_value_str(responseToolNameNode);
        
        if( rapidxml::xml_node<char>* responseInvocationIdNode = XML_FIRST_NODE(responseNode, "InvocationId") )
          response.invocationId = SpecUtils::xml_value_str(responseInvocationIdNode);
        
        if( rapidxml::xml_node<char> *responseParamsNode = XML_FIRST_NODE(responseNode, "ToolParameters") )
        {
          try
          {
            response.toolParameters = json::parse(responseParamsNode->value());
          }catch( const std::exception &e )
          {
            cout << "Failed to parse response tool parameters: " << e.what() << endl;
          }
        }
        
        conv->responses.push_back(response);
      }//XML_FOREACH_CHILD( responseNode, responsesNode, "Response" )
    }//if( rapidxml::xml_node<char> *responsesNode = XML_FIRST_NODE(convNode, "Responses") )

    conversations.push_back(conv);
  }//XML_FOREACH_CHILD(convNode, node, "Conversation")
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
