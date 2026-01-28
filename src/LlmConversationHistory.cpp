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

std::shared_ptr<LlmInteraction> LlmInteraction::create( Type t, const std::string& initialMessage, AgentType a )
{
  // Create the interaction with empty content field (legacy field for backward compatibility)
  std::shared_ptr<LlmInteraction> conv = std::shared_ptr<LlmInteraction>( new LlmInteraction(t, a) );

  // Add InitialRequest as the first response if message is not empty
  if( !initialMessage.empty() )
  {
    const LlmInteractionInitialRequest::RequestType reqType = (t == Type::System)
                                                                ? LlmInteractionInitialRequest::RequestType::System
                                                                : LlmInteractionInitialRequest::RequestType::User;
    auto initialReq = std::make_shared<LlmInteractionInitialRequest>(
      reqType,
      initialMessage,
      conv
    );
    conv->responses.push_back( initialReq );
  }

  return conv;
}

std::shared_ptr<LlmInteraction> LlmInteraction::createEmpty()
{
  return std::shared_ptr<LlmInteraction>( new LlmInteraction(Type::User, AgentType::MainAgent) );
}

std::shared_ptr<LlmInteraction> LlmConversationHistory::addUserMessageToMainConversation( const string &message )
{
  auto conv = LlmInteraction::create( LlmInteraction::Type::User, message, AgentType::MainAgent );
  conv->conversationId = "conv_" + std::to_string(chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now().time_since_epoch()).count());

  m_conversations.push_back(conv);

  return conv;
}


std::shared_ptr<LlmInteraction> LlmConversationHistory::addSystemMessageToMainConversation( const std::string &message )
{
  auto conv = LlmInteraction::create( LlmInteraction::Type::System, message, AgentType::MainAgent );
  conv->conversationId = "conv_" + std::to_string(chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now().time_since_epoch()).count());

  m_conversations.push_back(conv);

  return conv;
}

std::shared_ptr<LlmInteractionFinalResponse> LlmConversationHistory::addAssistantMessageWithThinking(const std::string& message,
                                                             const std::string& thinkingContent,
                                                             const std::string &rawContent,
                                                             std::shared_ptr<LlmInteraction> conversation )
{
  assert( conversation );
  if( !conversation )
  {
    cerr << "LlmConversationHistory::addAssistantMessageWithThinking: null conversation - not adding LLM response" << endl;
    return nullptr;
  }
  
  // Add message as a follow-up to an existing conversation
  auto response = std::make_shared<LlmInteractionFinalResponse>( message, conversation );
  response->setThinkingContent( thinkingContent );
  response->setRawContent( rawContent );
  conversation->responses.push_back(response);
  conversation->responseAdded.emit( response );
  
  return response;
}//void addAssistantMessageWithThinking(...)


shared_ptr<LlmToolRequest> LlmConversationHistory::addToolCalls( std::vector<LlmToolCall> &&toolCalls,
                                                                const std::string &rawResponseContent,
                                           const std::shared_ptr<LlmInteraction> &convo )
{
  assert( convo );
  if( !convo )
  {
    cerr << "LlmConversationHistory::addToolCalls: invalid LlmInteraction passed in." << endl;
    throw runtime_error( "LlmConversationHistory::addToolCalls: invalid LlmInteraction passed in." );
  }

  if( toolCalls.empty() )
  {
    cerr << "LlmConversationHistory::addToolCalls: called with empty toolCalls vector." << endl;
    return nullptr;
  }

  auto response = std::make_shared<LlmToolRequest>( convo );
  response->setToolCalls( std::move(toolCalls) );
  response->setRawContent( rawResponseContent );
  convo->responses.push_back(response);
  convo->responseAdded.emit( response );
  
  return response;
}//void addToolCalls(...)

std::shared_ptr<LlmToolResults> LlmConversationHistory::addToolResults( std::vector<LlmToolCall> &&toolResults,
                                             const std::string &jsonSentToLlm,
                                             const std::shared_ptr<LlmInteraction> &convo )
{
  assert( convo );
  if( !convo )
  {
    cerr << "LlmConversationHistory::addToolResults - got null conversation pointer" << endl;
    return nullptr;
  }

  if( toolResults.empty() )
  {
    cerr << "LlmConversationHistory::addToolResults - called with empty toolResults vector." << endl;
    return nullptr;
  }

  // Add this as a follow-up to an existing conversation
  auto response = std::make_shared<LlmToolResults>( convo );
  response->setToolCalls( std::move(toolResults) );
  response->setRawContent( jsonSentToLlm );
  convo->responses.push_back(response);
  convo->responseAdded.emit( response );
  
  return response;
}//void addToolResults(...)

std::shared_ptr<LlmInteractionError> LlmConversationHistory::addErrorMessage( const std::string &errorMessage,
                                             const std::string &rawResponseContent,
                                             const std::shared_ptr<LlmInteraction> &convo,
                                             LlmInteractionError::ErrorType errorType )
{
  assert( convo );
  if( !convo )
  {
    cerr << "LlmConversationHistory::addErrorMessage - got null conversation pointer: errorMessage='" << errorMessage << "'" << endl;
    return nullptr;
  }

  // Add this as a follow-up to an existing conversation
  auto response = std::make_shared<LlmInteractionError>( errorMessage, convo, errorType );
  response->setRawContent( rawResponseContent );
  convo->responses.push_back( response );
  convo->responseAdded.emit( response );

  return response;
}

std::shared_ptr<LlmInteractionAutoReply> LlmConversationHistory::addAutoReplyMessage( const std::string &promptMessage,
                                             const std::shared_ptr<LlmInteraction> &convo )
{
  assert( convo );
  if( !convo )
  {
    cerr << "LlmConversationHistory::addAutoReplyMessage - got null conversation pointer" << endl;
    return nullptr;
  }

  // Add this as a follow-up to an existing conversation
  auto response = std::make_shared<LlmInteractionAutoReply>( promptMessage, convo );
  convo->responses.push_back( response );
  convo->responseAdded.emit( response );

  return response;
}


void LlmConversationHistory::addTokenUsage( std::shared_ptr<LlmInteraction> conversation,
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


std::shared_ptr<LlmInteraction> LlmConversationHistory::findConversationByConversationId(const std::string& conversationId)
{
  for( shared_ptr<LlmInteraction> &conv : m_conversations)
  {
    if (conv->conversationId == conversationId)
      return conv;
  }
  return nullptr;
}

const std::vector<std::shared_ptr<LlmInteraction>>& LlmConversationHistory::getConversations() const {
  return m_conversations;
}

std::vector<std::shared_ptr<LlmInteraction>>& LlmConversationHistory::getConversations() {
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


void LlmConversationHistory::addConversationToLlmApiHistory( const LlmInteraction &conv, nlohmann::json &messages )
{
  assert( messages.is_array() );
  if( !messages.is_array() )
    throw logic_error( "addConversationToLlmApiHistory: messages must be an array." );

  // Add all responses in chronological order (including the initial request which is now part of responses)
  for( const std::shared_ptr<LlmInteractionTurn> &response : conv.responses )
  {
    // Skip responses marked to exclude from history
    if( response->excludeFromHistory() )
    {
      cout << "Skipping response excluded from history (type: " << static_cast<int>(response->type()) << ")" << endl;
      continue;
    }

    json responseMsg;

    switch( response->type() )
    {
      case LlmInteractionTurn::Type::InitialRequest:
      {
        const LlmInteractionInitialRequest *initialReq = dynamic_cast<const LlmInteractionInitialRequest*>(response.get());
        assert( initialReq );
        if( !initialReq )
        {
          cerr << "addConversationToLlmApiHistory: InitialRequest cast failed" << endl;
          break;
        }
        responseMsg["role"] = (initialReq->requestType() == LlmInteractionInitialRequest::RequestType::System) ? "system" : "user";
        responseMsg["content"] = initialReq->content();
        break;
      }

      case LlmInteractionTurn::Type::FinalLlmResponse:
      {
        const LlmInteractionFinalResponse *llmResp = dynamic_cast<const LlmInteractionFinalResponse*>(response.get());
        assert( llmResp );
        if( !llmResp )
        {
          cerr << "addConversationToLlmApiHistory: Assistant response cast failed" << endl;
          break;
        }
        responseMsg["role"] = "assistant";
        responseMsg["content"] = llmResp->content();
        break;
      }

      case LlmInteractionTurn::Type::ToolCall:
      {
        const LlmToolRequest *toolReq = dynamic_cast<const LlmToolRequest*>(response.get());
        assert( toolReq );
        if( !toolReq )
        {
          cerr << "addConversationToLlmApiHistory: ToolCall response cast failed" << endl;
          break;
        }
        responseMsg["role"] = "assistant";
        responseMsg["content"] = nlohmann::json::value_t::null; //Some LLM providers need this, notably asksage.
        responseMsg["tool_calls"] = json::array();

        // Iterate through all tool calls in this response (batched tool calls)
        for( const LlmToolCall &toolCall : toolReq->toolCalls() )
        {
          json toolCallJson;
          // Use just the invocationId to keep within OpenAI's 40-character limit
          toolCallJson["id"] = toolCall.invocationId;
          toolCallJson["type"] = "function";
          toolCallJson["function"]["name"] = toolCall.toolName;
          toolCallJson["function"]["arguments"] = toolCall.toolParameters.dump();
          responseMsg["tool_calls"].push_back(toolCallJson);
        }
        break;
      }

      case LlmInteractionTurn::Type::ToolResult:
      {
        const LlmToolResults *toolResult = dynamic_cast<const LlmToolResults*>(response.get());
        assert( toolResult );
        if( !toolResult )
        {
          cerr << "addConversationToLlmApiHistory: ToolResult response cast failed" << endl;
          break;
        }
        // Each tool result needs its own message with role="tool"
        // We need to add multiple messages for batched tool results
        for( const LlmToolCall &toolRes : toolResult->toolCalls() )
        {
          json toolResultMsg;
          toolResultMsg["role"] = "tool";
          toolResultMsg["tool_call_id"] = toolRes.invocationId;
          toolResultMsg["content"] = toolRes.content;

          // The Ask Sage tool-call *might* expect a different return format...
          //toolResultMsg["content"] = toolRes.toolName;
          //toolResultMsg["tool_response"] = toolRes.content;

          messages.push_back(toolResultMsg);
        }
        // Skip the normal push_back at the end since we've already added the messages
        continue;
      }

      case LlmInteractionTurn::Type::Error:
      {
        const LlmInteractionError *errorResp = dynamic_cast<const LlmInteractionError*>(response.get());
        assert( errorResp );
        if( !errorResp )
        {
          cerr << "addConversationToLlmApiHistory: Error response cast failed" << endl;
          break;
        }

        responseMsg["role"] = "assistant";
        responseMsg["content"] = "Error: " + errorResp->errorMessage();
        break;
      }

      case LlmInteractionTurn::Type::AutoReply:
      {
        const LlmInteractionAutoReply *autoReply = dynamic_cast<const LlmInteractionAutoReply*>(response.get());
        assert( autoReply );
        if( !autoReply )
        {
          cerr << "addConversationToLlmApiHistory: AutoReply response cast failed" << endl;
          break;
        }
        responseMsg["role"] = "user";
        responseMsg["content"] = autoReply->content();
        break;
      }
    }//switch( response->type() )

    messages.push_back(responseMsg);
  }//for( const std::shared_ptr<LlmInteractionTurn> &response : conv.responses )
}//nlohmann::json addConversationToLlmApiHistory(...)


nlohmann::json LlmConversationHistory::toApiFormat() const {
  json messages = json::array();

  for (const shared_ptr<LlmInteraction> &conv : m_conversations)
  {
    assert( conv );
    addConversationToLlmApiHistory( *conv, messages );
  }
  
  return messages;
}

void LlmConversationHistory::toXml(rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc) const {
  toXml(m_conversations, parent, doc);
}

void LlmConversationHistory::toXml( const vector<shared_ptr<LlmInteraction>> &conversations,
                                   rapidxml::xml_node<char>* parent, rapidxml::xml_document<char>* doc)
{
  rapidxml::xml_node<char>* historyNode = doc->allocate_node(rapidxml::node_element, "LlmHistory");
  parent->append_node(historyNode);

  cout << "Serializing " << conversations.size() << " conversations to XML" << endl;

  for( const shared_ptr<LlmInteraction> &conv : conversations )
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

    // Add finish time
    if( conv->finishTime.has_value() )
    {
      const auto finishTimeT = chrono::system_clock::to_time_t( conv->finishTime.value() );
      const string finishTimeStr = std::to_string( finishTimeT );
      XmlUtils::append_attrib(convNode, "finishTime", finishTimeStr );
    }

    // Add responses
    rapidxml::xml_node<char> *responsesNode = nullptr; //Will create at first message
    for( const std::shared_ptr<LlmInteractionTurn> &response : conv->responses )
    {
      if( !responsesNode ) // create "Responses" at the first message, so we dont create it if responses is empty
      {
        responsesNode = doc->allocate_node(rapidxml::node_element, "Responses");
        convNode->append_node(responsesNode);
      }

      rapidxml::xml_node<char>* responseNode = doc->allocate_node(rapidxml::node_element, "Response");
      responsesNode->append_node(responseNode);

      // Add type attribute
      XmlUtils::append_attrib(responseNode, "type", responseTypeToString(response->type()) );

      // Add timestamp attribute
      const auto responseTimeT = chrono::system_clock::to_time_t(response->timestamp());
      const string responseTimeStr = std::to_string(responseTimeT);
      XmlUtils::append_attrib(responseNode, "timestamp", responseTimeStr );

      // Add excludeFromHistory flag if true
      if( response->excludeFromHistory() )
        XmlUtils::append_attrib(responseNode, "excludeFromHistory", "true" );

      // Add thinking content (common to all types)
      if( !response->thinkingContent().empty() )
        XmlUtils::append_string_node( responseNode, "ThinkingContent", response->thinkingContent() );

      // Add raw content (JSON sent to LLM or received from LLM)
      if( !response->rawContent().empty() )
        XmlUtils::append_string_node( responseNode, "RawContent", response->rawContent() );

      // Add call duration if available
      if( response->callDuration().has_value() )
        XmlUtils::append_attrib( responseNode, "callDurationMs", std::to_string(response->callDuration().value().count()) );

      // Handle type-specific fields using dynamic_cast
      switch( response->type() )
      {
        case LlmInteractionTurn::Type::InitialRequest:
        {
          const LlmInteractionInitialRequest *initialReq = dynamic_cast<const LlmInteractionInitialRequest*>(response.get());
          if( initialReq )
          {
            if( !initialReq->content().empty() )
              XmlUtils::append_string_node( responseNode, "Content", initialReq->content() );
            // Add request type attribute (system or user)
            const string reqTypeStr = (initialReq->requestType() == LlmInteractionInitialRequest::RequestType::System) ? "system" : "user";
            XmlUtils::append_attrib( responseNode, "requestType", reqTypeStr );
          }
          break;
        }

        case LlmInteractionTurn::Type::FinalLlmResponse:
        {
          const LlmInteractionFinalResponse *llmResp = dynamic_cast<const LlmInteractionFinalResponse*>(response.get());
          if( llmResp && !llmResp->content().empty() )
            XmlUtils::append_string_node( responseNode, "Content", llmResp->content() );
          break;
        }

        case LlmInteractionTurn::Type::ToolCall:
        case LlmInteractionTurn::Type::ToolResult:
        {
          const LlmToolCall *toolCalls = nullptr;
          size_t numToolCalls = 0;

          if( response->type() == LlmInteractionTurn::Type::ToolCall )
          {
            const LlmToolRequest *toolReq = dynamic_cast<const LlmToolRequest*>(response.get());
            if( toolReq && !toolReq->toolCalls().empty() )
            {
              toolCalls = toolReq->toolCalls().data();
              numToolCalls = toolReq->toolCalls().size();
            }
          }else
          {
            const LlmToolResults *toolRes = dynamic_cast<const LlmToolResults*>(response.get());
            if( toolRes && !toolRes->toolCalls().empty() )
            {
              toolCalls = toolRes->toolCalls().data();
              numToolCalls = toolRes->toolCalls().size();
            }
          }

          if( toolCalls && numToolCalls > 0 )
          {
            rapidxml::xml_node<char> *toolCallsNode = doc->allocate_node(rapidxml::node_element, "ToolCalls");
            responseNode->append_node(toolCallsNode);

            for( size_t i = 0; i < numToolCalls; ++i )
            {
              const LlmToolCall &toolCall = toolCalls[i];
              rapidxml::xml_node<char> *toolCallNode = doc->allocate_node(rapidxml::node_element, "ToolCall");
              toolCallsNode->append_node(toolCallNode);

              // Add status attribute
              XmlUtils::append_attrib( toolCallNode, "status", callStatusToString(toolCall.status) );

              if( !toolCall.toolName.empty() )
                XmlUtils::append_string_node( toolCallNode, "ToolName", toolCall.toolName );

              if( !toolCall.invocationId.empty() )
                XmlUtils::append_string_node( toolCallNode, "InvocationId", toolCall.invocationId );

              if( !toolCall.toolParameters.empty() )
                XmlUtils::append_string_node( toolCallNode, "ToolParameters", toolCall.toolParameters.dump() );

              if( !toolCall.content.empty() )
                XmlUtils::append_string_node( toolCallNode, "ToolContent", toolCall.content );

              if( toolCall.executionDuration.has_value() )
                XmlUtils::append_attrib( toolCallNode, "executionDurationMs", std::to_string(toolCall.executionDuration.value().count()) );
            }//for( loop over toolCalls )
          }//if( toolCalls && numToolCalls > 0 )
          break;
        }

        case LlmInteractionTurn::Type::Error:
        {
          const LlmInteractionError *errorResp = dynamic_cast<const LlmInteractionError*>(response.get());
          if( errorResp && !errorResp->errorMessage().empty() )
            XmlUtils::append_string_node( responseNode, "ErrorMessage", errorResp->errorMessage() );
          break;
        }

        case LlmInteractionTurn::Type::AutoReply:
        {
          const LlmInteractionAutoReply *autoReply = dynamic_cast<const LlmInteractionAutoReply*>(response.get());
          if( autoReply && !autoReply->content().empty() )
            XmlUtils::append_string_node( responseNode, "Content", autoReply->content() );
          break;
        }
      }//switch( response->type() )
    }//for( const std::shared_ptr<LlmInteractionTurn> &response : conv->responses )
  }//for( const shared_ptr<LlmInteraction> &conv : conversations )
}//void toXml( const vector<shared_ptr<LlmInteraction>> &conversations ... )

void LlmConversationHistory::fromXml( const rapidxml::xml_node<char> *node )
{
  fromXml(node, m_conversations);
}

void LlmConversationHistory::fromXml( const rapidxml::xml_node<char> *node, std::vector<std::shared_ptr<LlmInteraction>> &conversations )
{
  // TODO: this function needs a bit of refactoring (LlmInteraction should get its own fromXml function) and we should be more exacting in requing the various fields and throw exceptions if they are not there, or invalid
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

    auto conv = LlmInteraction::createEmpty(); // Create empty conversation for deserialization

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

    // Read finish time
    if( rapidxml::xml_attribute<char>* finishTimeAttr = XML_FIRST_ATTRIB(convNode, "finishTime") )
    {
      const string finishTime_str = SpecUtils::xml_value_str(finishTimeAttr);
      const auto finishTimeT = static_cast<time_t>( std::stoll( finishTime_str.c_str() ) );
      conv->finishTime = chrono::system_clock::from_time_t(finishTimeT);
    }else
    {
      // If loaded from XML without finishTime, mark as completed with current time for backward compatibility
      conv->finishTime = std::nullopt;
    }

    // Read responses
    if( rapidxml::xml_node<char> *responsesNode = XML_FIRST_NODE(convNode, "Responses") )
    {
      XML_FOREACH_CHILD( responseNode, responsesNode, "Response" )
      {
        // Read response type first to determine which derived class to create
        LlmInteractionTurn::Type responseType = LlmInteractionTurn::Type::FinalLlmResponse; // Default
        if( rapidxml::xml_attribute<char>* responseTypeAttr = XML_FIRST_ATTRIB(responseNode,"type") )
          responseType = stringToResponseType( SpecUtils::xml_value_str(responseTypeAttr) );

        // Create the appropriate derived class based on type
        std::shared_ptr<LlmInteractionTurn> response;

        switch( responseType )
        {
          case LlmInteractionTurn::Type::InitialRequest:
          {
            // Read request type attribute (system or user)
            LlmInteractionInitialRequest::RequestType reqType = LlmInteractionInitialRequest::RequestType::User;
            if( rapidxml::xml_attribute<char>* reqTypeAttr = XML_FIRST_ATTRIB(responseNode, "requestType") )
            {
              const string reqTypeStr = SpecUtils::xml_value_str(reqTypeAttr);
              reqType = (reqTypeStr == "system") ? LlmInteractionInitialRequest::RequestType::System
                                                 : LlmInteractionInitialRequest::RequestType::User;
            }
            response = std::make_shared<LlmInteractionInitialRequest>( reqType, "", conv );
            break;
          }
          case LlmInteractionTurn::Type::FinalLlmResponse:
            response = std::make_shared<LlmInteractionFinalResponse>( "", conv );
            break;
          case LlmInteractionTurn::Type::ToolCall:
            response = std::make_shared<LlmToolRequest>( conv );
            break;
          case LlmInteractionTurn::Type::ToolResult:
            response = std::make_shared<LlmToolResults>( conv );
            break;
          case LlmInteractionTurn::Type::Error:
            response = std::make_shared<LlmInteractionError>( "", conv );
            break;
          case LlmInteractionTurn::Type::AutoReply:
            response = std::make_shared<LlmInteractionAutoReply>( "", conv );
            break;
        }

        // Read response timestamp (common to all types, but protected, so we can't set it directly - will use current time from constructor)
        // TODO: Add a setter or make timestamp accessible for deserialization

        // Read thinking content (common to all types)
        if( rapidxml::xml_node<char>* thinkingContentNode = XML_FIRST_NODE(responseNode,"ThinkingContent") )
          response->setThinkingContent( SpecUtils::xml_value_str(thinkingContentNode) );

        // Read raw content (JSON sent to LLM or received from LLM) - also check for old "JsonSentToLlm" name for backward compatibility
        if( rapidxml::xml_node<char> *rawContentNode = XML_FIRST_NODE(responseNode, "RawContent") )
          response->setRawContent( SpecUtils::xml_value_str(rawContentNode) );
        else if( rapidxml::xml_node<char> *jsonSentNode = XML_FIRST_NODE(responseNode, "JsonSentToLlm") )
          response->setRawContent( SpecUtils::xml_value_str(jsonSentNode) );

        // Read call duration - also check for old "apiCallDurationMs" name for backward compatibility
        if( rapidxml::xml_attribute<char> *durationAttr = XML_FIRST_ATTRIB(responseNode, "callDurationMs") )
        {
          const string duration_str = SpecUtils::xml_value_str(durationAttr);
          response->setCallDuration( std::chrono::milliseconds( std::stoll(duration_str) ) );
        }
        else if( rapidxml::xml_attribute<char> *durationAttr = XML_FIRST_ATTRIB(responseNode, "apiCallDurationMs") )
        {
          const string duration_str = SpecUtils::xml_value_str(durationAttr);
          response->setCallDuration( std::chrono::milliseconds( std::stoll(duration_str) ) );
        }

        // Read excludeFromHistory flag
        if( rapidxml::xml_attribute<char> *excludeAttr = XML_FIRST_ATTRIB(responseNode, "excludeFromHistory") )
        {
          const string exclude_str = SpecUtils::xml_value_str(excludeAttr);
          response->setExcludeFromHistory( exclude_str == "true" || exclude_str == "1" );
        }

        // Read type-specific fields
        switch( responseType )
        {
          case LlmInteractionTurn::Type::InitialRequest:
          {
            LlmInteractionInitialRequest *initialReq = dynamic_cast<LlmInteractionInitialRequest*>(response.get());
            if( initialReq )
            {
              if( rapidxml::xml_node<char>* contentNode = XML_FIRST_NODE(responseNode,"Content") )
                initialReq->setContent( SpecUtils::xml_value_str(contentNode) );
            }
            break;
          }

          case LlmInteractionTurn::Type::FinalLlmResponse:
          {
            LlmInteractionFinalResponse *llmResp = dynamic_cast<LlmInteractionFinalResponse*>(response.get());
            if( llmResp )
            {
              if( rapidxml::xml_node<char>* responseContentNode = XML_FIRST_NODE(responseNode,"Content") )
                llmResp->setContent( SpecUtils::xml_value_str(responseContentNode) );
            }
            break;
          }

          case LlmInteractionTurn::Type::ToolCall:
          case LlmInteractionTurn::Type::ToolResult:
          {
            // Read tool calls (for both ToolCall and ToolResult types)
            if( rapidxml::xml_node<char> *toolCallsNode = XML_FIRST_NODE(responseNode, "ToolCalls") )
            {
              std::vector<LlmToolCall> toolCalls;

              XML_FOREACH_CHILD( toolCallNode, toolCallsNode, "ToolCall" )
              {
                string toolName, invocationId, toolContent;
                nlohmann::json toolParameters;
                std::optional<std::chrono::milliseconds> executionDuration;
                LlmToolCall::CallStatus status = LlmToolCall::CallStatus::Success; // Default for backward compatibility

                // Read status attribute
                if( rapidxml::xml_attribute<char> *statusAttr = XML_FIRST_ATTRIB(toolCallNode, "status") )
                {
                  const string status_str = SpecUtils::xml_value_str(statusAttr);
                  status = stringToCallStatus(status_str);
                }

                if( rapidxml::xml_node<char> *toolNameNode = XML_FIRST_NODE(toolCallNode, "ToolName") )
                  toolName = SpecUtils::xml_value_str(toolNameNode);

                if( rapidxml::xml_node<char> *invocationIdNode = XML_FIRST_NODE(toolCallNode, "InvocationId") )
                  invocationId = SpecUtils::xml_value_str(invocationIdNode);

                if( rapidxml::xml_node<char> *paramsNode = XML_FIRST_NODE(toolCallNode, "ToolParameters") )
                {
                  try
                  {
                    toolParameters = json::parse(paramsNode->value());
                  }catch( const std::exception &e )
                  {
                    cout << "Failed to parse tool parameters: " << e.what() << endl;
                  }
                }

                if( rapidxml::xml_node<char> *contentNode = XML_FIRST_NODE(toolCallNode, "ToolContent") )
                  toolContent = SpecUtils::xml_value_str(contentNode);

                if( rapidxml::xml_attribute<char> *execDurationAttr = XML_FIRST_ATTRIB(toolCallNode, "executionDurationMs") )
                {
                  const string duration_str = SpecUtils::xml_value_str(execDurationAttr);
                  executionDuration = std::chrono::milliseconds( std::stoll(duration_str) );
                }

                LlmToolCall toolCall( toolName, invocationId, toolParameters );
                toolCall.status = status;
                toolCall.content = toolContent;
                toolCall.executionDuration = executionDuration;
                toolCalls.push_back(toolCall);
              }//XML_FOREACH_CHILD( toolCallNode, toolCallsNode, "ToolCall" )

              if( responseType == LlmInteractionTurn::Type::ToolCall )
              {
                LlmToolRequest *toolReq = dynamic_cast<LlmToolRequest*>(response.get());
                if( toolReq )
                  toolReq->setToolCalls( std::move(toolCalls) );
              }else
              {
                LlmToolResults *toolRes = dynamic_cast<LlmToolResults*>(response.get());
                if( toolRes )
                  toolRes->setToolCalls( std::move(toolCalls) );
              }
            }//if( rapidxml::xml_node<char> *toolCallsNode ... )
            break;
          }

          case LlmInteractionTurn::Type::Error:
          {
            LlmInteractionError *errorResp = dynamic_cast<LlmInteractionError*>(response.get());
            if( errorResp )
            {
              // Check for new "ErrorMessage" node or old "Content" node for backward compatibility
              if( rapidxml::xml_node<char>* errorMsgNode = XML_FIRST_NODE(responseNode,"ErrorMessage") )
                errorResp->setErrorMessage( SpecUtils::xml_value_str(errorMsgNode) );
              else if( rapidxml::xml_node<char>* responseContentNode = XML_FIRST_NODE(responseNode,"Content") )
                errorResp->setErrorMessage( SpecUtils::xml_value_str(responseContentNode) );
            }
            break;
          }

          case LlmInteractionTurn::Type::AutoReply:
          {
            LlmInteractionAutoReply *autoReply = dynamic_cast<LlmInteractionAutoReply*>(response.get());
            if( autoReply )
            {
              if( rapidxml::xml_node<char>* responseContentNode = XML_FIRST_NODE(responseNode,"Content") )
                autoReply->setContent( SpecUtils::xml_value_str(responseContentNode) );
            }
            break;
          }
        }//switch( responseType )

        conv->responses.push_back(response);
      }//XML_FOREACH_CHILD( responseNode, responsesNode, "Response" )
    }//if( rapidxml::xml_node<char> *responsesNode = XML_FIRST_NODE(convNode, "Responses") )

    conversations.push_back(conv);
  }//XML_FOREACH_CHILD(convNode, node, "Conversation")
}

std::string LlmConversationHistory::conversationTypeToString(LlmInteraction::Type type) {
  switch (type) {
    case LlmInteraction::Type::System: return "system";
    case LlmInteraction::Type::User: return "user";
    default: return "unknown";
  }
}

LlmInteraction::Type LlmConversationHistory::stringToConversationType(const std::string& str) {
  if (str == "system") return LlmInteraction::Type::System;
  if (str == "user") return LlmInteraction::Type::User;
  return LlmInteraction::Type::User; // Default fallback
}

std::string LlmConversationHistory::responseTypeToString(LlmInteractionTurn::Type type) {
  switch (type) {
    case LlmInteractionTurn::Type::InitialRequest: return "initial_request";
    case LlmInteractionTurn::Type::FinalLlmResponse: return "final_response";
    case LlmInteractionTurn::Type::ToolCall: return "tool_call";
    case LlmInteractionTurn::Type::ToolResult: return "tool_result";
    case LlmInteractionTurn::Type::Error: return "error";
    case LlmInteractionTurn::Type::AutoReply: return "auto_reply";
    default: return "unknown";
  }
}

LlmInteractionTurn::Type LlmConversationHistory::stringToResponseType(const std::string& str) {
  if (str == "initial_request") return LlmInteractionTurn::Type::InitialRequest;
  if (str == "final_response") return LlmInteractionTurn::Type::FinalLlmResponse;
  if (str == "tool_call") return LlmInteractionTurn::Type::ToolCall;
  if (str == "tool_result") return LlmInteractionTurn::Type::ToolResult;
  if (str == "error") return LlmInteractionTurn::Type::Error;
  if (str == "auto_reply") return LlmInteractionTurn::Type::AutoReply;
  return LlmInteractionTurn::Type::FinalLlmResponse; // Default fallback
}

std::string LlmConversationHistory::callStatusToString(LlmToolCall::CallStatus status) {
  switch (status) {
    case LlmToolCall::CallStatus::Pending: return "pending";
    case LlmToolCall::CallStatus::Success: return "success";
    case LlmToolCall::CallStatus::Error: return "error";
    default: return "unknown";
  }
}

LlmToolCall::CallStatus LlmConversationHistory::stringToCallStatus(const std::string& str) {
  if (str == "pending") return LlmToolCall::CallStatus::Pending;
  if (str == "success") return LlmToolCall::CallStatus::Success;
  if (str == "error") return LlmToolCall::CallStatus::Error;
  return LlmToolCall::CallStatus::Success; // Default fallback
}
