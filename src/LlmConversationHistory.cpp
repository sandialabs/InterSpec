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

std::shared_ptr<LlmInteraction> LlmInteraction::createAgent( const Type t, std::string initialMessage, const AgentType a, const LlmConfig * const config )
{
  // Create the interaction with empty content field (legacy field for backward compatibility)
  std::shared_ptr<LlmInteraction> conv = std::shared_ptr<LlmInteraction>( new LlmInteraction(t, a) );

  if( config )
  {
    for( const LlmConfig::AgentConfig &agent : config->agents )
    {
      if( agent.type == a )
      {
        if( agent.state_machine )
          conv->state_machine = agent.state_machine->copy();
        break;
      }
    }
  }//if( config )

  if( conv->state_machine )
  {
    const string &currentState = conv->state_machine->getCurrentState();
    const vector<string> allowedTransitions = conv->state_machine->getAllowedTransitions();
    const string guidance = conv->state_machine->getPromptGuidanceForCurrentState();

    string stateContent = "You are currently in workflow state '" + currentState
      + "'.\n\n";

    string transitionList;
    for( size_t i = 0; i < allowedTransitions.size(); ++i )
      transitionList += string((i > 0) ? ", " : "") + allowedTransitions[i];

    if( !transitionList.empty() )
      stateContent += "Allowed transitions from this state: " + transitionList + ".\n\n";

    if( !guidance.empty() )
      stateContent += "Guidance for current state: " + guidance;

    stateContent += "\n\nPlease continue making tool calls and reasoning for this current state, and then explicitly call the `set_workflow_state` tool-call, with an allowed transition state, to transition to the next state.";

    if( !initialMessage.empty() )
      initialMessage += "\n\n";

    initialMessage += stateContent;
  }

  // Add InitialRequest as the first response if message is not empty
  if( !initialMessage.empty() )
  {
    const LlmInteractionInitialRequest::RequestType reqType = (t == Type::System)
                                                                ? LlmInteractionInitialRequest::RequestType::System
                                                                : LlmInteractionInitialRequest::RequestType::User;
    auto initialReq = std::make_shared<LlmInteractionInitialRequest>( reqType, initialMessage, conv );
    conv->responses.push_back( initialReq );
  }

  return conv;
}


shared_ptr<LlmInteraction> LlmInteraction::createUserInteractionForMainConvo( const std::string &initialMessage )
{
  return createAgent( LlmInteraction::Type::User, initialMessage, AgentType::MainAgent, nullptr );
}


shared_ptr<LlmInteraction> LlmInteraction::createSystemInteractionForMainConvo( const std::string &initialMessage )
{
  return createAgent( LlmInteraction::Type::System, initialMessage, AgentType::MainAgent, nullptr );
}


std::shared_ptr<LlmInteraction> LlmInteraction::createEmpty()
{
  return std::shared_ptr<LlmInteraction>( new LlmInteraction(Type::User, AgentType::MainAgent) );
}


std::shared_ptr<LlmInteraction> LlmInteraction::shallowClone() const
{
  // Create a new interaction with the same type/agent, then copy scalar fields.
  // The responses vector is copied so the clone has an independent container -
  // new turns appended by LlmInterface will not appear in the original, and vice versa.
  // Turn objects themselves are shared (their content is immutable after creation).
  // Signals, handlers, and runtime-only tracking fields are intentionally not copied.
  std::shared_ptr<LlmInteraction> clone(
    new LlmInteraction( type, agent_type ) );

  clone->conversationId   = conversationId;
  clone->timestamp        = timestamp;
  clone->finishTime       = finishTime;
  clone->promptTokens     = promptTokens;
  clone->completionTokens = completionTokens;
  clone->totalTokens      = totalTokens;
  clone->responses        = responses; // shared_ptr copies; turns are immutable read-only objects
  if( state_machine )
    clone->state_machine = state_machine->copy();

  return clone;
}

std::shared_ptr<LlmInteraction> LlmConversationHistory::addUserMessageToMainConversation( const string &message )
{
  auto conv = LlmInteraction::createUserInteractionForMainConvo( message );
  conv->conversationId = "conv_" + std::to_string(chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now().time_since_epoch()).count());

  m_conversations.push_back(conv);

  return conv;
}


std::shared_ptr<LlmInteraction> LlmConversationHistory::addUserMessageToMainConversation( const string &message,
                                                                                           vector<LlmToolCall::ImageContent> images )
{
  shared_ptr<LlmInteraction> conv = addUserMessageToMainConversation( message );

  if( !images.empty() && conv && !conv->responses.empty() )
  {
    LlmInteractionInitialRequest *initialReq
      = dynamic_cast<LlmInteractionInitialRequest *>( conv->responses.front().get() );
    assert( initialReq );
    if( initialReq )
      initialReq->setImageContent( std::move( images ) );
  }

  return conv;
}


std::shared_ptr<LlmInteraction> LlmConversationHistory::addSystemMessageToMainConversation( const std::string &message )
{
  auto conv = LlmInteraction::createSystemInteractionForMainConvo( message );
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

        if( !initialReq->imageContent().empty()
            && (initialReq->requestType() == LlmInteractionInitialRequest::RequestType::User) )
        {
          // Multi-part content array: text + image(s)
          json contentArray = json::array();

          json textBlock;
          textBlock["type"] = "text";
          textBlock["text"] = initialReq->content();
          contentArray.push_back( textBlock );

          for( const LlmToolCall::ImageContent &img : initialReq->imageContent() )
          {
            json imageBlock;
            imageBlock["type"] = "image_url";
            imageBlock["image_url"] = {
              {"url", "data:" + img.mimeType + ";base64," + img.base64Data}
            };
            contentArray.push_back( imageBlock );
          }

          responseMsg["content"] = contentArray;
        }
        else
        {
          responseMsg["content"] = initialReq->content();
        }
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

        const std::string &thinkingSig = llmResp->thinkingSignature();
        const std::string &thinkingText = llmResp->thinkingContent();

        if( !thinkingSig.empty() && !thinkingText.empty() )
        {
          // Claude extended thinking - emit content block array
          json contentArray = json::array();

          json thinkingBlock;
          thinkingBlock["type"] = "thinking";
          thinkingBlock["thinking"] = thinkingText;
          thinkingBlock["signature"] = thinkingSig;
          contentArray.push_back( thinkingBlock );

          if( !llmResp->content().empty() )
          {
            json textBlock;
            textBlock["type"] = "text";
            textBlock["text"] = llmResp->content();
            contentArray.push_back( textBlock );
          }

          responseMsg["content"] = contentArray;
        }
        else
        {
          // Standard string format
          responseMsg["content"] = llmResp->content();
          // Note: reasoning_details is NOT included for FinalLlmResponse - it's only needed
          // during tool call sequences, not for final responses
        }
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

        const std::string &thinkingSig = toolReq->thinkingSignature();
        const std::string &thinkingText = toolReq->thinkingContent();
        const std::string &reasoningText = toolReq->reasoningContent();

        if( !thinkingSig.empty() && !thinkingText.empty() )
        {
          // Claude extended thinking with tool use - emit content block array
          json contentArray = json::array();

          json thinkingBlock;
          thinkingBlock["type"] = "thinking";
          thinkingBlock["thinking"] = thinkingText;
          thinkingBlock["signature"] = thinkingSig;
          contentArray.push_back( thinkingBlock );

          // Add tool_use blocks to the content array (Claude format)
          for( const LlmToolCall &toolCall : toolReq->toolCalls() )
          {
            json toolUseBlock;
            toolUseBlock["type"] = "tool_use";
            toolUseBlock["id"] = toolCall.invocationId;
            toolUseBlock["name"] = toolCall.toolName;
            toolUseBlock["input"] = toolCall.toolParameters;
            contentArray.push_back( toolUseBlock );
          }

          responseMsg["content"] = contentArray;
        }else if( !reasoningText.empty() )
        {
          // Some models MUST include reasoning_content during tool calls
          // Check if this tool call is followed by a tool result (same turn)
          // If not, we shouldn't include reasoning_content (it's a new turn)
          // (this logic has not actually been tested using one of these models)
          bool isFollowedByToolResult = false;
          auto it = std::find( conv.responses.begin(), conv.responses.end(), response );
          if( it != conv.responses.end() )
          {
            ++it;
            if( it != conv.responses.end() )
              isFollowedByToolResult = ((*it)->type() == LlmInteractionTurn::Type::ToolResult);
          }

          responseMsg["content"] = nlohmann::json::value_t::null;
          if( isFollowedByToolResult )
            responseMsg["reasoning_content"] = reasoningText;

          responseMsg["tool_calls"] = json::array();
          for( const LlmToolCall &toolCall : toolReq->toolCalls() )
          {
            json toolCallJson;
            toolCallJson["id"] = toolCall.invocationId;
            toolCallJson["type"] = "function";
            toolCallJson["function"]["name"] = toolCall.toolName;
            toolCallJson["function"]["arguments"] = toolCall.toolParameters.dump();
            if( !toolCall.extra_content.is_null() )
              toolCallJson["extra_content"] = toolCall.extra_content;
            responseMsg["tool_calls"].push_back( toolCallJson );
          }
        }else
        {
          // Standard OpenAI format with tool_calls array
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
            if( !toolCall.extra_content.is_null() )
              toolCallJson["extra_content"] = toolCall.extra_content;
            responseMsg["tool_calls"].push_back( toolCallJson );
          }

          // OpenRouter: Add reasoning_details if present, but ONLY when the conversation
          // ends with a ToolResult (i.e., we're sending tool results back to the LLM).
          // It should only be included on the most recent tool request, not on past ones.
          const std::string &reasoningDetailsJson = toolReq->reasoningDetails();
          if( !reasoningDetailsJson.empty() && !conv.responses.empty()
            && (conv.responses.back()->type() == LlmInteractionTurn::Type::ToolResult) )
          {
            // Check if this ToolCall immediately precedes the final ToolResult(s)
            auto it = std::find( conv.responses.begin(), conv.responses.end(), response );
            if( (it != conv.responses.end()) && (++it != conv.responses.end())
              && ((*it)->type() == LlmInteractionTurn::Type::ToolResult) )
            {
              try
              {
                responseMsg["reasoning_details"] = json::parse( reasoningDetailsJson );
              }catch( const std::exception &e )
              {
                // If parsing fails, skip reasoning_details
                cerr << "\n\nFailed to parse 'reasoning_details' - so not including in response - this shouldnt happen." << endl;
              }
            }
          }
        }
        break;
      }//case LlmInteractionTurn::Type::ToolCall:

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

          if( toolRes.imageContent.has_value() )
          {
            // Format as OpenAI content array with text + image
            json contentArray = json::array();

            json textBlock;
            textBlock["type"] = "text";
            textBlock["text"] = toolRes.content;
            contentArray.push_back( textBlock );

            json imageBlock;
            imageBlock["type"] = "image_url";
            imageBlock["image_url"] = {
              {"url", "data:" + toolRes.imageContent->mimeType + ";base64," + toolRes.imageContent->base64Data}
            };
            contentArray.push_back( imageBlock );

            toolResultMsg["content"] = contentArray;
          }
          else
          {
            toolResultMsg["content"] = toolRes.content;
          }

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

      case LlmInteractionTurn::Type::ConversationSummary:
      {
        const LlmConversationSummary *summary = dynamic_cast<const LlmConversationSummary *>( response.get() );
        assert( summary );
        if( !summary )
        {
          cerr << "addConversationToLlmApiHistory: ConversationSummary cast failed" << endl;
          break;
        }
        responseMsg["role"] = "user";
        responseMsg["content"] = "[System: Summary of previous conversation history]\n\n"
          + summary->summary();
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

      // Add thinking signature (Claude extended thinking - must be passed back unmodified)
      if( !response->thinkingSignature().empty() )
        XmlUtils::append_string_node( responseNode, "ThinkingSignature", response->thinkingSignature() );

      // Add reasoning content for models that require it MUST be passed back
      if( !response->reasoningContent().empty() )
        XmlUtils::append_string_node( responseNode, "ReasoningContent", response->reasoningContent() );

      // Add reasoning details JSON array (OpenRouter format - must be passed back unmodified)
      if( !response->reasoningDetails().empty() )
        XmlUtils::append_string_node( responseNode, "ReasoningDetails", response->reasoningDetails() );

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

              if( !toolCall.toolParameters.empty() )
                XmlUtils::append_string_node( toolCallNode, "ExtraContent", toolCall.extra_content.dump() );

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

        case LlmInteractionTurn::Type::ConversationSummary:
        {
          const LlmConversationSummary *summary = dynamic_cast<const LlmConversationSummary *>( response.get() );
          if( summary )
          {
            if( !summary->summary().empty() )
              XmlUtils::append_string_node( responseNode, "Content", summary->summary() );
            XmlUtils::append_attrib( responseNode, "summarizedCount", std::to_string( summary->summarizedCount() ) );

            const auto earliestT = chrono::system_clock::to_time_t( summary->earliestTimestamp() );
            XmlUtils::append_attrib( responseNode, "earliestTimestamp", std::to_string( earliestT ) );

            const auto latestT = chrono::system_clock::to_time_t( summary->latestTimestamp() );
            XmlUtils::append_attrib( responseNode, "latestTimestamp", std::to_string( latestT ) );

            // Recursively serialize the archived original conversations
            if( !summary->originalConversations().empty() )
              toXml( summary->originalConversations(), responseNode, doc );
          }
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

  // The XML structure is: <LlmConversationHistory> -> <LlmHistory> -> <Conversation>...
  // So we need to navigate to the LlmHistory child first, if present.
  const rapidxml::xml_node<char> *historyNode = XML_FIRST_NODE( node, "LlmHistory" );
  if( !historyNode )
  {
    // Fall back to looking for Conversation directly under node, for backwards compatibility
    historyNode = node;
  }

  int convCount = 0;
  XML_FOREACH_CHILD(convNode, historyNode, "Conversation")
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
          case LlmInteractionTurn::Type::ConversationSummary:
          {
            size_t summarizedCount = 0;
            std::chrono::system_clock::time_point earliest = std::chrono::system_clock::now();
            std::chrono::system_clock::time_point latest = earliest;

            if( rapidxml::xml_attribute<char> *countAttr = XML_FIRST_ATTRIB(responseNode, "summarizedCount") )
              summarizedCount = static_cast<size_t>( std::stoull( SpecUtils::xml_value_str(countAttr) ) );
            if( rapidxml::xml_attribute<char> *earlyAttr = XML_FIRST_ATTRIB(responseNode, "earliestTimestamp") )
              earliest = chrono::system_clock::from_time_t( static_cast<time_t>( std::stoll( SpecUtils::xml_value_str(earlyAttr) ) ) );
            if( rapidxml::xml_attribute<char> *lateAttr = XML_FIRST_ATTRIB(responseNode, "latestTimestamp") )
              latest = chrono::system_clock::from_time_t( static_cast<time_t>( std::stoll( SpecUtils::xml_value_str(lateAttr) ) ) );

            response = std::make_shared<LlmConversationSummary>( "", summarizedCount, earliest, latest, conv );
            break;
          }
        }

        // Read response timestamp (common to all types, but protected, so we can't set it directly - will use current time from constructor)
        // TODO: Add a setter or make timestamp accessible for deserialization

        // Read thinking content (common to all types)
        if( rapidxml::xml_node<char>* thinkingContentNode = XML_FIRST_NODE(responseNode,"ThinkingContent") )
          response->setThinkingContent( SpecUtils::xml_value_str(thinkingContentNode) );

        // Read thinking signature (Claude extended thinking - must be passed back unmodified)
        if( rapidxml::xml_node<char>* thinkingSigNode = XML_FIRST_NODE(responseNode, "ThinkingSignature") )
          response->setThinkingSignature( SpecUtils::xml_value_str(thinkingSigNode) );

        // Read reasoning content for models that require it must be passed back
        if( rapidxml::xml_node<char>* reasoningNode = XML_FIRST_NODE(responseNode, "ReasoningContent") )
          response->setReasoningContent( SpecUtils::xml_value_str(reasoningNode) );

        // Read reasoning details JSON array (OpenRouter format - must be passed back unmodified)
        if( rapidxml::xml_node<char>* reasoningDetailsNode = XML_FIRST_NODE(responseNode, "ReasoningDetails") )
          response->setReasoningDetails( SpecUtils::xml_value_str(reasoningDetailsNode) );

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
                nlohmann::json toolParameters, extraContent;
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

                if( rapidxml::xml_node<char> *extraContentNode = XML_FIRST_NODE(toolCallNode, "ExtraContent") )
                {
                  try
                  {
                    extraContent = json::parse(extraContentNode->value());
                  }catch( const std::exception &e )
                  {
                    cout << "Failed to parse extra_content: " << e.what() << endl;
                  }
                }

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
                toolCall.extra_content = extraContent;
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

          case LlmInteractionTurn::Type::ConversationSummary:
          {
            LlmConversationSummary *summary = dynamic_cast<LlmConversationSummary *>( response.get() );
            if( summary )
            {
              if( rapidxml::xml_node<char> *contentNode = XML_FIRST_NODE(responseNode, "Content") )
                summary->setSummary( SpecUtils::xml_value_str(contentNode) );

              // Restore archived original conversations if present
              if( rapidxml::xml_node<char> *archivedNode = XML_FIRST_NODE(responseNode, "LlmHistory") )
              {
                vector<shared_ptr<LlmInteraction>> originals;
                fromXml( archivedNode, originals );
                summary->setOriginalConversations( std::move( originals ) );
              }
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
    case LlmInteractionTurn::Type::ConversationSummary: return "conversation_summary";
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
  if (str == "conversation_summary") return LlmInteractionTurn::Type::ConversationSummary;
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


nlohmann::json LlmConversationHistory::buildMessagesForMainAgent(
  const std::shared_ptr<LlmInteraction> &convo,
  const std::string &systemPrompt ) const
{
  assert( convo );
  if( !convo )
    throw logic_error( "buildMessagesForMainAgent: null conversation passed in." );

  json messages = json::array();

  if( !systemPrompt.empty() )
  {
    json systemMsg;
    systemMsg["role"] = "system";
    systemMsg["content"] = systemPrompt;
    messages.push_back( systemMsg );
  }

  for( const shared_ptr<LlmInteraction> &previous_conversation : m_conversations )
  {
    assert( previous_conversation );
    addConversationToLlmApiHistory( *previous_conversation, messages );

    if( convo == previous_conversation )
      break;
  }

  return messages;
}//json buildMessagesForMainAgent(...)


bool LlmConversationHistory::hasSummary() const
{
  if( m_conversations.empty() )
    return false;

  const shared_ptr<LlmInteraction> &first = m_conversations.front();
  if( !first || first->responses.empty() )
    return false;

  return (first->responses.front()->type() == LlmInteractionTurn::Type::ConversationSummary);
}//bool hasSummary()


size_t LlmConversationHistory::estimateTokenCount() const
{
  // Rough heuristic: ~4 characters per token for English text
  static const double chars_per_token = 4.0;

  size_t total_chars = 0;

  for( const shared_ptr<LlmInteraction> &conv : m_conversations )
  {
    if( !conv )
      continue;

    for( const shared_ptr<LlmInteractionTurn> &turn : conv->responses )
    {
      if( !turn || turn->excludeFromHistory() )
        continue;

      switch( turn->type() )
      {
        case LlmInteractionTurn::Type::InitialRequest:
        {
          const LlmInteractionInitialRequest *req = dynamic_cast<const LlmInteractionInitialRequest *>( turn.get() );
          if( req )
            total_chars += req->content().size();
          break;
        }

        case LlmInteractionTurn::Type::FinalLlmResponse:
        {
          const LlmInteractionFinalResponse *resp = dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
          if( resp )
            total_chars += resp->content().size();
          break;
        }

        case LlmInteractionTurn::Type::ToolCall:
        {
          const LlmToolRequest *toolReq = dynamic_cast<const LlmToolRequest *>( turn.get() );
          if( toolReq )
          {
            for( const LlmToolCall &tc : toolReq->toolCalls() )
              total_chars += tc.toolName.size() + tc.toolParameters.dump().size();
          }
          break;
        }

        case LlmInteractionTurn::Type::ToolResult:
        {
          const LlmToolResults *toolRes = dynamic_cast<const LlmToolResults *>( turn.get() );
          if( toolRes )
          {
            for( const LlmToolCall &tc : toolRes->toolCalls() )
              total_chars += tc.content.size();
          }
          break;
        }

        case LlmInteractionTurn::Type::Error:
        {
          const LlmInteractionError *err = dynamic_cast<const LlmInteractionError *>( turn.get() );
          if( err )
            total_chars += err->errorMessage().size();
          break;
        }

        case LlmInteractionTurn::Type::AutoReply:
        {
          const LlmInteractionAutoReply *reply = dynamic_cast<const LlmInteractionAutoReply *>( turn.get() );
          if( reply )
            total_chars += reply->content().size();
          break;
        }

        case LlmInteractionTurn::Type::ConversationSummary:
        {
          const LlmConversationSummary *summary = dynamic_cast<const LlmConversationSummary *>( turn.get() );
          if( summary )
            total_chars += summary->summary().size();
          break;
        }
      }//switch( turn->type() )

      // Include thinking content in the estimate
      total_chars += turn->thinkingContent().size();
    }//for( turns )
  }//for( conversations )

  return static_cast<size_t>( total_chars / chars_per_token );
}//size_t estimateTokenCount()


bool LlmConversationHistory::shouldSummarize( const size_t contextLengthLimit ) const
{
  if( contextLengthLimit == 0 )
    return false;

  // Need at least a few non-summary conversations before summarizing
  const size_t startIdx = hasSummary() ? 1 : 0;
  const size_t nonSummaryCount = m_conversations.size() - startIdx;
  if( nonSummaryCount <= 6 )
    return false;

  // Check actual prompt token count from the most recent completed conversation.
  // Walk backwards to find the most recent conversation with token data.
  for( auto it = m_conversations.rbegin(); it != m_conversations.rend(); ++it )
  {
    const shared_ptr<LlmInteraction> &conv = *it;
    if( conv && conv->promptTokens.has_value() )
    {
      const size_t promptTokens = conv->promptTokens.value();
      const size_t threshold = static_cast<size_t>( contextLengthLimit * 0.85 );
      return (promptTokens >= threshold);
    }
  }//for( walk backwards )

  // Fallback: estimate from character counts
  const size_t estimatedTokens = estimateTokenCount();
  const size_t threshold = static_cast<size_t>( contextLengthLimit * 0.85 );
  return (estimatedTokens >= threshold);
}//bool shouldSummarize(...)


std::string LlmConversationHistory::buildSummarizationPrompt( const size_t keepRecentCount ) const
{
  const size_t startIdx = hasSummary() ? 1 : 0;
  const size_t nonSummaryCount = m_conversations.size() - startIdx;

  if( nonSummaryCount <= keepRecentCount )
    return "";

  const size_t endIdx = m_conversations.size() - keepRecentCount;

  string prompt = "Please provide a concise summary of the following conversation history. "
    "Capture: (1) what the user asked, (2) key tools used and their outcomes, "
    "(3) conclusions reached, (4) any important state changes (peaks added/removed, "
    "calibrations changed, nuclides identified). "
    "The summary should be self-contained so you can understand the full context "
    "without the original messages.\n\n";

  // Include any existing summary as context
  if( hasSummary() && !m_conversations.empty() )
  {
    const shared_ptr<LlmInteraction> &summaryConvo = m_conversations.front();
    if( summaryConvo && !summaryConvo->responses.empty() )
    {
      const LlmConversationSummary *prevSummary =
        dynamic_cast<const LlmConversationSummary *>( summaryConvo->responses.front().get() );
      if( prevSummary && !prevSummary->summary().empty() )
      {
        prompt += "Previous summary of even earlier conversations:\n"
          + prevSummary->summary() + "\n\n";
      }
    }
  }//if( hasSummary() )

  prompt += "Conversations to summarize:\n\n";

  for( size_t i = startIdx; i < endIdx; ++i )
  {
    const shared_ptr<LlmInteraction> &conv = m_conversations[i];
    if( !conv )
      continue;

    prompt += "--- Exchange " + to_string( i - startIdx + 1 ) + " ---\n";

    for( const shared_ptr<LlmInteractionTurn> &turn : conv->responses )
    {
      if( !turn )
        continue;

      switch( turn->type() )
      {
        case LlmInteractionTurn::Type::InitialRequest:
        {
          const LlmInteractionInitialRequest *req =
            dynamic_cast<const LlmInteractionInitialRequest *>( turn.get() );
          if( req )
          {
            const string &content = req->content();
            const string typeStr = (req->requestType() == LlmInteractionInitialRequest::RequestType::System)
              ? "System" : "User";
            // Truncate long messages to keep the summarization prompt manageable
            if( content.size() > 500 )
              prompt += typeStr + ": " + content.substr( 0, 500 ) + "...\n";
            else
              prompt += typeStr + ": " + content + "\n";
          }
          break;
        }

        case LlmInteractionTurn::Type::ToolCall:
        {
          const LlmToolRequest *toolReq = dynamic_cast<const LlmToolRequest *>( turn.get() );
          if( toolReq )
          {
            prompt += "Tools called: ";
            for( size_t j = 0; j < toolReq->toolCalls().size(); ++j )
            {
              if( j > 0 )
                prompt += ", ";
              prompt += toolReq->toolCalls()[j].toolName;
            }
            prompt += "\n";
          }
          break;
        }

        case LlmInteractionTurn::Type::ToolResult:
        {
          const LlmToolResults *toolRes = dynamic_cast<const LlmToolResults *>( turn.get() );
          if( toolRes )
          {
            for( const LlmToolCall &tc : toolRes->toolCalls() )
            {
              const string statusStr = (tc.status == LlmToolCall::CallStatus::Success) ? "Success" : "Error";
              const string &content = tc.content;
              if( content.size() > 300 )
                prompt += "  " + tc.toolName + " (" + statusStr + "): " + content.substr( 0, 300 ) + "...\n";
              else
                prompt += "  " + tc.toolName + " (" + statusStr + "): " + content + "\n";
            }
          }
          break;
        }

        case LlmInteractionTurn::Type::FinalLlmResponse:
        {
          const LlmInteractionFinalResponse *resp =
            dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
          if( resp )
          {
            const string &content = resp->content();
            if( content.size() > 500 )
              prompt += "Assistant: " + content.substr( 0, 500 ) + "...\n";
            else
              prompt += "Assistant: " + content + "\n";
          }
          break;
        }

        case LlmInteractionTurn::Type::Error:
        {
          const LlmInteractionError *err = dynamic_cast<const LlmInteractionError *>( turn.get() );
          if( err )
            prompt += "Error: " + err->errorMessage() + "\n";
          break;
        }

        case LlmInteractionTurn::Type::AutoReply:
        case LlmInteractionTurn::Type::ConversationSummary:
          break;
      }//switch( turn->type() )
    }//for( turns )

    prompt += "\n";
  }//for( conversations to summarize )

  return prompt;
}//string buildSummarizationPrompt(...)


void LlmConversationHistory::applySummary( const std::string &summaryText, const size_t keepRecentCount )
{
  const size_t startIdx = hasSummary() ? 1 : 0;
  const size_t nonSummaryCount = m_conversations.size() - startIdx;

  if( nonSummaryCount <= keepRecentCount )
    return;

  const size_t endIdx = m_conversations.size() - keepRecentCount;

  // Compute timestamp range of summarized conversations
  const auto earliest = m_conversations[startIdx]->timestamp;
  const auto latest = m_conversations[endIdx - 1]->timestamp;

  // Count total summarized, including any previously summarized count
  size_t totalSummarizedCount = endIdx - startIdx;
  if( hasSummary() )
  {
    const LlmConversationSummary *prevSummary =
      dynamic_cast<const LlmConversationSummary *>( m_conversations.front()->responses.front().get() );
    if( prevSummary )
      totalSummarizedCount += prevSummary->summarizedCount();
  }

  // Create summary LlmInteraction
  shared_ptr<LlmInteraction> summaryConvo = LlmInteraction::createEmpty();
  summaryConvo->type = LlmInteraction::Type::System;
  summaryConvo->agent_type = AgentType::MainAgent;
  summaryConvo->conversationId = "summary_" + to_string(
    chrono::duration_cast<chrono::milliseconds>(
      chrono::system_clock::now().time_since_epoch() ).count() );
  summaryConvo->timestamp = earliest;
  summaryConvo->finishTime = chrono::system_clock::now();

  // Collect the original conversations being summarized (excluding any previous summary placeholder)
  vector<shared_ptr<LlmInteraction>> originals;
  originals.reserve( endIdx - startIdx );

  // If there was a previous summary, pull its archived originals forward
  if( hasSummary() )
  {
    LlmConversationSummary *prevSummary =
      dynamic_cast<LlmConversationSummary *>( m_conversations.front()->responses.front().get() );
    if( prevSummary )
    {
      vector<shared_ptr<LlmInteraction>> prevOriginals = prevSummary->takeOriginalConversations();
      originals.insert( originals.end(),
        std::make_move_iterator( prevOriginals.begin() ),
        std::make_move_iterator( prevOriginals.end() ) );
    }
  }

  for( size_t i = startIdx; i < endIdx; ++i )
    originals.push_back( m_conversations[i] );

  shared_ptr<LlmConversationSummary> summaryTurn = make_shared<LlmConversationSummary>(
    summaryText, totalSummarizedCount, earliest, latest, summaryConvo );
  summaryTurn->setOriginalConversations( std::move( originals ) );
  summaryConvo->responses.push_back( summaryTurn );

  // Erase conversations [0..endIdx) and insert summary at position 0
  m_conversations.erase( m_conversations.begin(), m_conversations.begin() + static_cast<ptrdiff_t>( endIdx ) );
  m_conversations.insert( m_conversations.begin(), summaryConvo );
}//void applySummary(...)
