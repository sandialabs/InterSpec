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

#include <regex>
#include <chrono>
#include <memory>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#if( PERFORM_DEVELOPER_CHECKS )
#include <mutex>
#endif

#include <Wt/WApplication>
#include <Wt/WResource>
#include <Wt/WWebWidget>
#include <Wt/Http/Request>
#include <Wt/Http/Response>

#ifdef USE_JS_BRIDGE_FOR_LLM
#include <Wt/WJavaScriptSlot>
#endif

#ifdef USE_WT_HTTP_FOR_LLM
#include <Wt/Http/Client>
#include <Wt/Http/Message>
#endif

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/LlmConversationHistory.h"


#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
#include "SpecUtils/DateTime.h"
#endif
#include "SpecUtils/StringAlgo.h"

static_assert( USE_LLM_INTERFACE, "This file shouldnt be compiled unless USE_LLM_INTERFACE is true" );

using namespace std;
using json = nlohmann::json;

LlmInterface::LlmInterface( InterSpec* interspec, const std::shared_ptr<const LlmConfig> &config )
  : m_interspec(interspec),
    m_config( config ),
    m_tool_registry( config ? make_shared<LlmTools::ToolRegistry>(*config) : shared_ptr<LlmTools::ToolRegistry>() ),
    m_history(std::make_shared<LlmConversationHistory>()),
#ifdef USE_JS_BRIDGE_FOR_LLM
    m_responseSignal(std::make_unique<Wt::JSignal<std::string, int>>(interspec, "llmResponse")),
#endif
#ifdef USE_WT_HTTP_FOR_LLM
    m_httpClient(std::make_unique<Wt::Http::Client>(interspec)),
    m_currentRequestId(0),
#endif
    m_nextRequestId(1)
{
  if( !m_interspec )
    throw std::runtime_error("InterSpec instance cannot be null");
  
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );
  
  if( !m_tool_registry )
    throw std::logic_error( "LlmInterface: couldnt create tool registry from config" );
  
  // Get session ID through WApplication
  string sessionId = "unknown";
  if (Wt::WApplication* app = Wt::WApplication::instance()) {
    sessionId = app->sessionId();
  }
  cout << "LlmInterface created for session: " << sessionId << endl;
  
#ifdef USE_JS_BRIDGE_FOR_LLM
  // Connect the JavaScript response signal
  m_responseSignal->connect(this, &LlmInterface::handleJavaScriptResponse);
  
  setupJavaScriptBridge();
#endif

#ifdef USE_WT_HTTP_FOR_LLM
  // Configure HTTP client
  m_httpClient->setTimeout(120); // 120 second timeout
  m_httpClient->setMaximumResponseSize(1024 * 1024); // 1MB max response
  
  // Connect the done signal to our response handler
  m_httpClient->done().connect(this, &LlmInterface::handleWtHttpClientResponse);
#endif
}

LlmInterface::~LlmInterface() {
  cout << "LlmInterface destroyed" << endl;
}


std::shared_ptr<const LlmTools::ToolRegistry> LlmInterface::toolRegistry()
{
  return m_tool_registry;
}


void LlmInterface::sendUserMessage( const std::string &message )
{
  if( !isConfigured() )
    throw std::logic_error( "LLM interface is not properly configured" );
  
  cout << "User message: " << message << endl;
  
  // Add to history FIRST to ensure it's there before any async responses
  shared_ptr<LlmConversationStart> convo = m_history->addUserMessageToMainConversation( message );
  assert( convo );
  cout << "Starting new conversation with ID: " << convo->conversationId << endl;
  cout << "Added user message to history. History now has " << m_history->getConversations().size() << " conversations" << endl;
  
  // Debug: Show what's actually in history
  cout << "Current history contents:" << endl;
  for (size_t i = 0; i < m_history->getConversations().size(); ++i) {
    const auto& conv = m_history->getConversations()[i];
    string typeStr;
    switch (conv->type) {
      case LlmConversationStart::Type::User: typeStr = "user"; break;
      case LlmConversationStart::Type::System: typeStr = "system"; break;
    }
    cout << "  " << i << ". " << typeStr << ": " << conv->content.substr(0, 50) << "..." << " (responses: " << conv->responses.size() << ")" << endl;
  }
  
  // Build API request (buildMessagesArray will include the message from history + current message)
  json requestJson = buildMessagesArray( convo );
  
  // Make tracked API call
  int requestId = makeTrackedApiCall( requestJson, convo );
  
  cout << "Sent user message with request ID: " << requestId << endl;
}

void LlmInterface::sendSystemMessage( const std::string &message )
{
  if( !isConfigured() )
    throw std::runtime_error("LLM interface is not properly configured");
  
  cout << "System message: " << message << endl;
  
  std::shared_ptr<LlmConversationStart> convo = m_history->addSystemMessageToMainConversation( message );
  
  // Build API request (marked as system-generated)
  json requestJson = buildMessagesArray( convo );
  
  // Make tracked API call
  int requestId = makeTrackedApiCall( requestJson, convo );
  cout << "Sent system message with request ID: " << requestId << endl;
}


std::shared_ptr<LlmConversationHistory> LlmInterface::getHistory() const {
  return m_history;
}

void LlmInterface::setHistory(std::shared_ptr<LlmConversationHistory> history) {
  m_history = history ? history : std::make_shared<LlmConversationHistory>();
}

bool LlmInterface::isConfigured() const {
  return (m_config && m_config->llmApi.enabled && !m_config->llmApi.apiEndpoint.empty());
}


Wt::Signal<>& LlmInterface::responseReceived() {
  return m_responseReceived;
}

bool LlmInterface::isRequestPending(int requestId) const {
  return m_pendingRequests.find(requestId) != m_pendingRequests.end();
}


/** Manually serialize JSON request with specific field ordering for LLM provider caching.

 Some LLM providers cache parts of requests between calls to reduce costs, but only if
 the beginning of the request string exactly matches previous requests. This function
 ensures that cacheable fields (tools, tool_choice, model, max_completion_tokens, max_tokens)
 appear first and in a reliable order.

 @param requestJson The JSON object to serialize
 @return JSON string with fields in cache-friendly order
 */
static std::string serializeRequestForCaching( const nlohmann::json &requestJson )
{
  std::string result = "{";
  bool first = true;

  // Helper lambda to add a field if it exists
  auto addField = [&]( const std::string &key )
  {
    if( !requestJson.contains(key) )
      return;

    if( !first )
      result += ",";
    first = false;

    result += "\"" + key + "\":";
    result += requestJson[key].dump();
  };

  // Add fields in priority order for caching (these should stay stable across requests)
  addField( "model" );
  addField( "max_completion_tokens" );
  addField( "max_tokens" );
  addField( "tools" );
  addField( "tool_choice" );

  // Add remaining fields (messages, etc.)
  for( auto it = requestJson.begin(); it != requestJson.end(); ++it )
  {
    const std::string &key = it.key();

    // Skip fields we already added
    if( key == "model" || key == "max_completion_tokens" || key == "max_tokens" ||
        key == "tools" || key == "tool_choice" )
      continue;

    if( !first )
      result += ",";
    first = false;

    result += "\"" + key + "\":";
    result += it.value().dump();
  }

  result += "}";
  return result;
}//serializeRequestForCaching(...)


void LlmInterface::makeApiCallWithId(const nlohmann::json& requestJson, int requestId)
{
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );

#if( PERFORM_DEVELOPER_CHECKS )
  // Track request string prefix matching for cache verification
  static std::mutex s_cache_check_mutex;
  static std::string s_previous_request_str;
#endif

  cout << "=== Making LLM API Call with ID " << requestId << " ===" << endl;
  cout << "Endpoint: " << m_config->llmApi.apiEndpoint << endl;
  cout << "Request JSON:" << endl;
  nlohmann::json debugJson = requestJson;
  // Truncate system prompt to first 100 characters
  if( debugJson.contains("messages") && debugJson["messages"].is_array() && debugJson["messages"].size() && debugJson["messages"][0].is_object() )
  {
    if( debugJson["messages"][0].contains("content") && debugJson["messages"][0]["content"].is_string() )
    {
      std::string systemPrompt = debugJson["messages"][0]["content"].get<std::string>();
      if( systemPrompt.length() > 100 )
        debugJson["messages"][0]["content"] = systemPrompt.substr(0, 100) + "...";
    }
  }

  // Truncate all message contents to 100 characters for readability
  if( debugJson.contains("messages") && debugJson["messages"].is_array() )
  {
    for( auto &msg : debugJson["messages"] )
    {
      if( msg.is_object() && msg.contains("content") && msg["content"].is_string() )
      {
        std::string content = msg["content"].get<std::string>();
        if( content.length() > 100 )
        {
          msg["content"] = content.substr(0, 100) + "...";
        }
      }
    }
  }

  if( debugJson.contains("tools") )
  {
    const size_t ntools = debugJson["tools"].size();
    debugJson.erase("tools");
    debugJson["tools"] = std::to_string(ntools) + " tools ommited";
  }
  if( debugJson.contains("tool_choice") )
  {
    const size_t ntools = debugJson["tool_choice"].size();
    debugJson.erase("tool_choice");
    debugJson["tool_choice"] = std::to_string(ntools) + " tool_choice ommited";
  }
  cout << debugJson.dump(2) << endl;
  cout << "=========================" << endl;

#ifdef USE_JS_BRIDGE_FOR_LLM
  // Serialize request with consistent field ordering for LLM provider caching
  const std::string requestStr = serializeRequestForCaching( requestJson );

  /*
#if( PERFORM_DEVELOPER_CHECKS )
  // Developer check: Compare this request with the previous one to verify caching optimization
  {
    std::lock_guard<std::mutex> lock( s_cache_check_mutex );

    if( !s_previous_request_str.empty() )
    {
      // Find how many characters match from the beginning
      const size_t min_len = std::min( requestStr.length(), s_previous_request_str.length() );
      size_t matching_chars = 0;
      for( size_t i = 0; i < min_len; ++i )
      {
        if( requestStr[i] != s_previous_request_str[i] )
          break;
        ++matching_chars;
      }

      const double match_percent = (matching_chars * 100.0) / std::max( requestStr.length(), s_previous_request_str.length() );

      cout << "=== LLM Request Cache Analysis ===" << endl;
      cout << "Current request length: " << requestStr.length() << " chars" << endl;
      cout << "Previous request length: " << s_previous_request_str.length() << " chars" << endl;
      cout << "Matching prefix: " << matching_chars << " chars (" << std::fixed << std::setprecision(1) << match_percent << "%)" << endl;

      if( matching_chars > 0 )
      {
        // Show the matching prefix (truncated for readability)
        const size_t preview_len = std::min( matching_chars, size_t(150) );
        cout << "Matching prefix: \"" << requestStr.substr(0, preview_len);
        if( preview_len < matching_chars )
          cout << "...\" (+" << (matching_chars - preview_len) << " more chars)";
        else
          cout << "\"";
        cout << endl;
      }

      // Show where the difference starts
      if( matching_chars < min_len )
      {
        const size_t diff_preview_len = std::min( size_t(50), min_len - matching_chars );
        cout << "First difference at position " << matching_chars << ":" << endl;
        cout << "  Current:  \"" << requestStr.substr(matching_chars, diff_preview_len) << "...\"" << endl;
        cout << "  Previous: \"" << s_previous_request_str.substr(matching_chars, diff_preview_len) << "...\"" << endl;
      }

      cout << "===================================" << endl;
    }else
    {
      cout << "=== First LLM request in session (no cache comparison) ===" << endl;
    }

    // Store current request for next comparison
    s_previous_request_str = requestStr;
  }
#endif // PERFORM_DEVELOPER_CHECKS
   */

  // Use Wt::WWebWidget::jsStringLiteral to properly escape the JSON string for JavaScript
  const std::string jsLiteralRequestStr = Wt::WWebWidget::jsStringLiteral( requestStr );

  string jsCall =
    "var llmRequestData = " + jsLiteralRequestStr + ";\n"
    "window.llmHttpRequest('" + m_config->llmApi.apiEndpoint + "', llmRequestData, '" +
    m_config->llmApi.bearerToken + "', " + std::to_string(requestId) + ");";

  // Execute JavaScript to make the HTTP request
  auto app = Wt::WApplication::instance();
  if( !app )
    throw runtime_error( "Error: No WApplication instance available for JavaScript bridge" );

  app->doJavaScript(jsCall);
  cout << "JavaScript HTTP request with ID " << requestId << " initiated..." << endl;

#endif

#ifdef USE_WT_HTTP_FOR_LLM
  // Use Wt HTTP client to make the request
  cout << "Using Wt HTTP client with ID " << requestId << "..." << endl;

  try {
    // Create HTTP message with JSON body
    Wt::Http::Message request;

    // Set Content-Type header
    request.setHeader("Content-Type", "application/json");

    // Add Authorization header if bearer token is provided
    if (!m_config->llmApi.bearerToken.empty()) {
      request.setHeader("Authorization", "Bearer " + m_config->llmApi.bearerToken);
    }

    // Set request body to JSON with consistent field ordering for caching
    const std::string requestBody = serializeRequestForCaching( requestJson );
    request.addBodyText(requestBody);

    cout << "Request body length: " << requestBody.length() << " bytes" << endl;

    // Store request ID for response correlation (we'll use member variable for now)
    // TODO: Better approach would be to store this in a map with client instance
    m_currentRequestId = requestId;

    // Make the HTTP POST request
    bool success = m_httpClient->post(m_config->llmApi.apiEndpoint, request);

    if (success) {
      cout << "Wt HTTP POST request with ID " << requestId << " initiated successfully" << endl;
    } else {
      cout << "Failed to initiate Wt HTTP request" << endl;
      handleWtHttpError(requestId, "Failed to initiate HTTP request");
    }

  } catch (const std::exception& e) {
    cout << "Exception in Wt HTTP request: " << e.what() << endl;
    handleWtHttpError(requestId, std::string("Exception: ") + e.what());
  }
#endif
}

void LlmInterface::handleApiResponse( const std::string &response,
                                     const std::shared_ptr<LlmConversationStart> &conversation,
                                     const int requestId )
{
  assert( conversation );
  if( !conversation )
  {
    cerr << "LlmInterface::handleApiResponse: recieved null conversation for response: " << response << endl << endl;
    return;
  }
  
  size_t number_tool_calls = 0;
  
  try
  {
    json responseJson = json::parse(response);
    
    // Parse and accumulate token usage information if available
    if( responseJson.contains("usage") )
    {
      const auto &usage = responseJson["usage"];
      
      std::optional<int> promptTokens, completionTokens, totalTokens;
      if (usage.contains("prompt_tokens") && usage["prompt_tokens"].is_number())
        promptTokens = usage["prompt_tokens"].get<int>();
      if (usage.contains("completion_tokens") && usage["completion_tokens"].is_number())
        completionTokens = usage["completion_tokens"].get<int>();
      if (usage.contains("total_tokens") && usage["total_tokens"].is_number())
        totalTokens = usage["total_tokens"].get<int>();
      
      // Accumulate token usage for this conversation
      m_history->addTokenUsage( conversation, promptTokens, completionTokens, totalTokens );
      
      if( completionTokens.has_value() )
      {
        cout << "=== Token Usage This Call ===" << endl;
        cout << "Prompt tokens: " << (promptTokens.has_value() ? std::to_string(promptTokens.value()) : "N/A") << endl;
        cout << "Completion tokens: " << completionTokens.value() << endl;
        cout << "Total tokens: " << (totalTokens.has_value() ? std::to_string(totalTokens.value()) : "N/A") << endl;
        cout << "=============================" << endl;
      }
    }//if( responseJson.contains("usage") )
    
    
    if( responseJson.contains("choices") && !responseJson["choices"].empty() )
    {
      const json &choice = responseJson["choices"][0];
      
      if( choice.contains("message") )
      {
        const json &message = choice["message"];
        string role = message.value("role", "");
        string content;
        if( message.contains("content") && message["content"].is_string() )
          content = message["content"];

        if( role == "assistant" )
        {
          // Extract thinking content and clean content
          auto [cleanContent, thinkingContent] = extractThinkingAndContent(content);

          // Handle structured tool calls first (OpenAI format)
          if (message.contains("tool_calls"))
          {
            number_tool_calls += executeToolCallsAndSendResults( message["tool_calls"], conversation, requestId );
          }else
          {
            // Parse content for text-based tool requests (use cleaned content)
            number_tool_calls += parseContentForToolCallsAndSendResults( cleanContent, conversation, requestId );
          }
          
          // Add assistant message to history with thinking content and current agent name - only if there were no
          //  tool calls - if there were tool calls `executeToolCallsAndSendResults` will add these to the
          //  history (although we will lose the thinking content, but whatever)
          //  TODO: add thinking content to tool call LlmConversationResponse's
          if( number_tool_calls == 0 )
            m_history->addAssistantMessageWithThinking( cleanContent, thinkingContent, conversation );
        }//if( role == "assistant" )
      }//if( choice.contains("message") )
    }//if( responseJson.contains("choices") && !responseJson["choices"].empty() )
  }catch( const std::exception &e )
  {
    cout << "Error parsing LLM response: " << e.what()
    << "\n\tresponse="
    << response << endl << endl;
  }
  
  if( !number_tool_calls && conversation->conversation_completion_handler )
    conversation->conversation_completion_handler();
  
  // Only emit signal if there are no pending requests (i.e., this is the final response)
  if( m_pendingRequests.empty() )
  {
    cout << "No pending requests, emitting response received signal" << endl;
    m_responseReceived.emit();
  }else
  {
    cout << "Still have " << m_pendingRequests.size() << " pending requests, not emitting signal yet" << endl;
  }
}//void handleApiResponse(...)


size_t LlmInterface::executeToolCallsAndSendResults( const nlohmann::json &toolCalls,
                                    const std::shared_ptr<LlmConversationStart> &convo,
                                    const int parentRequestId )
{
  assert( wApp ); //Consistency requires we have the WApplication::UpdateLock
  
  assert( convo );
  if( !convo )
  {
    cerr << "LlmInterface::executeToolCallsAndSendResults: recieved null conversation for tool calls: " << toolCalls.dump(2) << endl << endl;
    return 0;
  }
  
  cout << "Executing " << toolCalls.size() << " tool calls" << endl;

  // Track executed tool calls for follow-up
  std::vector<std::string> executedToolCallIds;
  vector<int> subAgentRequestIds;

  for( const nlohmann::json &toolCall : toolCalls )
  {
    const string callId = toolCall.value("id", "");

    try
    {
      const string toolName = toolCall["function"]["name"];
      cout << "--- about to parse tool call arguments for toolName=" << toolName << endl;
      json arguments = json::parse(toolCall["function"]["arguments"].get<string>());
      cout << "--- done parsing tool call arguments for toolName=" << toolName << endl;

      cout << "Calling tool: " << toolName << " with ID: " << callId << endl;

      // Add tool call to history with conversation ID, invocation ID, and current agent name
      m_history->addToolCall( toolName, callId, arguments, convo );

      // Check if this is an invoke_* tool (sub-agent invocation) - handle specially
    
      if( SpecUtils::istarts_with(toolName, "invoke_") && (toolName.length() > 7) )
      {
        // Extract agent name from tool name (e.g., "invoke_NuclideId" -> "NuclideId")
        const string agentName = toolName.substr(7);
        const AgentType agent_type = stringToAgentType(agentName);
        assert( agent_type != AgentType::MainAgent );
        if( agent_type == AgentType::MainAgent )
          throw runtime_error( "Can not invoke a agent of type AgentType::MainAgent" );
        
        const string context = arguments.at("context").get<string>();
        const string task = arguments.at("task").get<string>();

        cout << "Detected sub-agent invocation: " << agentName << " via tool: " << toolName << endl;

        // Invoke sub-agent (returns request ID)
        const string combinedMessage = "Context:\n" + context + "\n\nTask:\n" + task + "\n\nWhen you are done, provide a summary of your reasoning, as well as how sure you are, or alternate possibilities to investigate.";
        shared_ptr<LlmConversationStart> sub_agent_convo = make_shared<LlmConversationStart>( convo->type, combinedMessage, agent_type );
        
        // TODO: make sure we dont have two sub-agent calls for the same type of agent - among other problems, this would make `sub_agent_convo->conversationId` non-unique
        const auto current_ticks = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        sub_agent_convo->conversationId = "subagent_" + agentName + "_" + std::to_string( current_ticks );
        
        const std::weak_ptr<LlmConversationStart> &parent_convo_wk = convo;
        const std::weak_ptr<LlmConversationStart> sub_agent_convo_wk = sub_agent_convo;
        
        const int subAgentRequestId = invokeSubAgent( sub_agent_convo );
        subAgentRequestIds.push_back( subAgentRequestId );
        
        
        {// Begin placeholder result that will be updated when sub-agent completes
          json result;
          result["status"] = "pending";
          result["message"] = "Sub-agent " + agentName + " is processing (will update when complete)";
          result["requestId"] = subAgentRequestId;
          
          LlmConversationResponse response(LlmConversationResponse::Type::ToolResult, result.dump(), convo );
          response.invocationId = callId;
          response.sub_agent_conversation = sub_agent_convo;
          
          convo->responses.push_back( std::move(response) );
        }// End placeholder result that will be updated when sub-agent completes
        
        // We will define a completion handler that will get called when the sub-agent conversation is complete.
        LlmInterface * self = this;
        sub_agent_convo->conversation_completion_handler = [parent_convo_wk, parentRequestId, subAgentRequestId, self, sub_agent_convo_wk, callId](){
          shared_ptr<LlmConversationStart> parent_conv = parent_convo_wk.lock();
          assert( parent_conv );
          if( !parent_conv )
          {
            cerr << "Sub-agent completed with parent conversation no longer alive - ignoring results." << endl;
            return;
          }
          
          shared_ptr<LlmConversationStart> sub_agent_convo = sub_agent_convo_wk.lock();
          assert( sub_agent_convo );
          if( !sub_agent_convo )
          {
            cerr << "Sub-agent cconversation no longer alive - ignoring results." << endl;
            return;
          }
          
          InterSpec *viewer = InterSpec::instance();
          if( !viewer )
          {
            cerr << "Sub-agent completed for dead session - ignoring results." << endl;
            return;
          }
          
          // TODO: maybe `LlmInterface` should use the Wt object life and/or boost signal/slot mechanism to protect against LlmInterface lifetime - and not this hack of a system of relying on checking stuff through the GUI
          LlmToolGui * const llm_gui = viewer->currentLlmTool();
          if( !llm_gui )
          {
            cerr << "Sub-agent completed and no llm_gui avaiable - ignoring results." << endl;
            return;
          }
          
          LlmInterface * const interface = llm_gui->llmInterface();
          if( interface != self )
          {
            cerr << "Sub-agent completed and with different LlmInterface now present - ignoring results." << endl;
            return;
          }
          
          // Here is where we fill in `parent_conv` with the result, and send back to the LLM
          cout << "Sub-agent complete (no tool calls left) - will extracting summary and continue main conversation." << endl;
          
          assert( !parent_conv->responses.empty() );
          LlmConversationResponse *response = nullptr;
          for( size_t rspns_index = 0; !response && (rspns_index < parent_conv->responses.size()); ++rspns_index )
          {
            if( (parent_conv->responses[rspns_index].invocationId == callId)
               && (parent_conv->responses[rspns_index].type == LlmConversationResponse::Type::ToolResult) )
            {
              response = &(parent_conv->responses[rspns_index]);
              assert( response->content.find( "(will update when complete)" ) != string::npos );
            }
          }//for( LlmConversationResponse &rspns : parent_conv->responses )
          
          assert( response );
          
          if( response )
          {
            assert( !sub_agent_convo->responses.empty() );
            
            json updatedResult;
            
            string subAgentSummary, subAgentThinkingContent;
            if( sub_agent_convo->responses.empty() )
            {
              updatedResult["status"] = "failed";
              updatedResult["summary"] = "Failed to get agent summary - because of internal logic error.";
              updatedResult["agentName"] = agentTypeToString( sub_agent_convo->agent_type );
            }else
            {
              string last_response = sub_agent_convo->responses.back().content;
              cout << "--- about to parse sub_agent_convo->responses='" << last_response << "'" << endl;
              SpecUtils::trim( last_response );
              nlohmann::json responseJson = nlohmann::json::object();
              try
              {
                if( !last_response.empty() )
                  responseJson = nlohmann::json::parse( last_response );
                cout << "--- Done parsing sub_agent_convo->responses" << endl;
                
                if( responseJson.contains("choices") && !responseJson["choices"].empty() )
                {
                  const json &choice = responseJson["choices"][0];
                  if( choice.contains("message") && choice["message"].contains("content") )
                  {
                    const string content = choice["message"]["content"].get<string>();
                    // Extract clean content (strip thinking tags if present)
                    auto [cleanContent, thinkingContent] = extractThinkingAndContent(content);
                    subAgentSummary = cleanContent;
                    subAgentThinkingContent = thinkingContent;
                  }
                }
                
                cout << "--- Done extracting JSON subAgentSummary='" << subAgentSummary << "'\n\n" << endl;
              }catch( std::exception &e )
              {
                // `LlmConversationResponse::content` may not / likely is not, JSON
                auto [cleanContent, thinkingContent] = extractThinkingAndContent(last_response);
                subAgentSummary = cleanContent;
                subAgentThinkingContent = thinkingContent;
                cout << "--- Done extracting non-JSON subAgentSummary='" << subAgentSummary << "'\n\n" << endl;
              }
              
              
              cout << "Sub-agent summary: " << subAgentSummary.substr(0, 100) << (subAgentSummary.length() > 100 ? "..." : "") << endl;
              cout << "subAgentThinkingContent='" << subAgentThinkingContent << "'" << endl;
              
              updatedResult["status"] = "completed";
              updatedResult["summary"] = subAgentSummary;
              updatedResult["agentName"] = agentTypeToString( sub_agent_convo->agent_type );
            }//if( sub_agent_convo->responses.empty() )
            
            response->content = updatedResult.dump();
            response->thinkingContent = subAgentThinkingContent;
          }else
          {
            cerr << "Failed to find tool-call response - shouldnt happen!!!" << endl;
          }//if( response )
          
        
          const auto defered_pos = interface->m_deferredToolResults.find(parentRequestId);
          assert( defered_pos != end(interface->m_deferredToolResults) );
          if( defered_pos == end(interface->m_deferredToolResults) )
          {
            cerr << "Sub-agent completed with no defered request being found..." << endl;
            return;
          }
          
          
          DeferredToolResult &deffered_result = defered_pos->second;
          vector<int> &subAgentToolCallIds = deffered_result.subAgentToolCallIds;
          const auto pos = std::find( begin(subAgentToolCallIds), end(subAgentToolCallIds), subAgentRequestId );
          assert( pos != end(subAgentToolCallIds) );
          subAgentToolCallIds.erase( pos );
          if( subAgentToolCallIds.empty() )
          {
            // No more agent calls is pending, lets cleanup the results
            interface->m_deferredToolResults.erase( defered_pos );
            
            cout << "=== Sending parent conversation with sub-agent results back to LLM!" << endl;
            interface->sendToolResultsToLLM( parent_conv );
          }
        };//sub_agent_convo->conversation_completion_handler

        cout << "Sub-agent invocation deferred - will resume main agent when complete" << endl;
      }else
      {
        // Normal tool execution
        json result = m_tool_registry->executeTool(toolName, arguments, m_interspec);

        LlmConversationResponse response(LlmConversationResponse::Type::ToolResult, result.dump(), convo );
        response.invocationId = callId;
        convo->responses.push_back( response );
      }
    }catch( const std::exception &e )
    {
      cout << "Tool execution error: " << e.what() << endl;

      json result;
      result["error"] = "Tool call failed: " + string(e.what());
      m_history->addToolResult( callId, result, convo );
    }//try / catch

    // Track this tool call for follow-up
    executedToolCallIds.push_back(callId);
  }

  // If we have a sub-agent invocation, defer sending results until sub-agent completes
  if( !subAgentRequestIds.empty() )
  {
    cout << "Deferring tool results - for sub-agent to complete" << endl;

    DeferredToolResult deferred;
    deferred.conversationId = convo->conversationId;
    deferred.toolCallIds = executedToolCallIds;
    deferred.subAgentToolCallIds = subAgentRequestIds; // The invoke_sub_agent call should be last

    m_deferredToolResults[parentRequestId] = deferred;

    // Don't send results yet - wait for sub-agent
  }else
  {
    // No sub-agent - send results immediately
    if (!executedToolCallIds.empty()) {
      cout << "Sending tool results back to LLM for " << executedToolCallIds.size() << " executed tools" << endl;
      sendToolResultsToLLM( convo );
    }
  }
  
  return toolCalls.size();
}

size_t LlmInterface::parseContentForToolCallsAndSendResults( const std::string &content, const std::shared_ptr<LlmConversationStart> &convo, const int requestId )
{
  cout << "Parsing content for text-based tool calls..." << endl;
  
  nlohmann::json tool_calls = nlohmann::json::array();
  // We will add entries to this array that look like:
  // { "id": "call_abc123", "type": "function", "function": { "name": "get_gamma_spectrum", "arguments": "{\"energy_range\": [0, 3000], \"detector\": \"HPGe\"}" }
  
  // Look for various tool call patterns
  std::vector<std::regex> toolCallPatterns = {
    // Pattern: [TOOL_REQUEST] {"name": "tool_name", "arguments": {...}} [END_TOOL_REQUEST]
    std::regex("\\[TOOL_REQUEST\\]\\s*([^\\[]+)\\s*\\[END_TOOL_REQUEST\\]"),
    // Pattern: [TOOL_REQUEST] tool_name {"param": "value"} [/TOOL_REQUEST]
    std::regex("\\[TOOL_REQUEST\\]\\s*(\\w+)\\s*(\\{[^}]*\\})\\s*\\[/TOOL_REQUEST\\]"),
    // Pattern: <tool_call name="tool_name">{"param": "value"}</tool_call>
    std::regex("<tool_call\\s+name=\"([^\"]+)\">([^<]*)</tool_call>"),
    // Pattern: <tool_call>{"name": "tool_name", "arguments": {...}}</tool_call>
    std::regex("<tool_call>\\s*([^<]+)\\s*</tool_call>"),
    // Pattern: Tool: tool_name Arguments: {"param": "value"}
    std::regex("Tool:\\s*(\\w+)\\s*Arguments:\\s*(\\{[^}]*\\})"),
    // Pattern: detected_peaks({"specType": "Foreground"})
    std::regex("(\\w+)\\s*\\(\\s*(\\{[^}]*\\})\\s*\\)")
  };
  
  std::smatch match;
  for( const auto &pattern : toolCallPatterns )
  {
    std::sregex_iterator iter(content.begin(), content.end(), pattern);
    std::sregex_iterator end;
    
    for( ; iter != end; ++iter )
    {
      const std::smatch& match = *iter;
      
      if( match.size() < 2 )
        continue;
      
      // Generate a unique invocation ID for this tool call
      const auto epoch_ticks = chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now().time_since_epoch()).count();
      const string invocationId = "text_call_" + std::to_string( epoch_ticks );
      
      string toolName;
      json arguments = json::object();
      
      try
      {
        if( match.size() == 2 )
        {
          // Single capture group - JSON object with name and arguments fields
          string toolCallJson = match[1].str();
          cout << "Found JSON tool call: " << toolCallJson << " -- and about to parse" << endl;
          
          json toolCallObj = json::parse(toolCallJson);
          cout << "Done parsing tool call JSON" << endl;
          
          if( toolCallObj.contains("name") )
            toolName = toolCallObj["name"];
          else
            throw std::runtime_error( "Invalid tool call JSON format - missing name" );
          
          if( toolCallObj.contains("arguments") )
            arguments = toolCallObj["arguments"];
        }else
        {
          assert( match.size() >= 3 );
          
          // Two capture groups - tool name and arguments separately
          toolName = match[1].str();
          string argumentsStr = match[2].str();
          
          // Parse arguments JSON
          if( !argumentsStr.empty() && (argumentsStr != "{}") )
          {
            cout << "About to parse alone tool call JSON" << endl;
            
            arguments = json::parse(argumentsStr);
            
            cout << "Done parsing alone tool call JSON" << endl;
          }
        }
        
        cout << "Found text-based tool call: " << toolName << " with args: " << arguments.dump() << endl;
      }catch( const std::exception &e )
      {
        cout << "Error parsing text-based tool call: " << e.what() << endl;
        if( toolName.empty() )
          toolName = match[0];
      }//try / catch
      
      
      // Add an entry to `tool_calls` that look like:
      // { "id": "call_abc123", "type": "function", "function": { "name": "get_gamma_spectrum", "arguments": "{...}" }
      // (note, adding function call, even if we got an exception extracting the tool call - we will let
      //  `executeToolCallsAndSendResults(...)` deal with errors).
      nlohmann::json function_def = nlohmann::json::object();
      function_def["name"] = toolName;
      function_def["arguments"] = std::move(arguments);
      
      nlohmann::json call_def = nlohmann::json::object();
      call_def["id"] = invocationId;
      call_def["type"] = "function";
      call_def["function"] = std::move(function_def);
      
      tool_calls.push_back( std::move(call_def) );
    }//for( ; iter != end; ++iter )
  }//for( const auto &pattern : toolCallPatterns )
  
  cout << "=== Complete extracting " << tool_calls.size() << " text-based tool calls ===" << endl;
  
  const size_t num_calls = executeToolCallsAndSendResults( tool_calls, convo, requestId );

  return num_calls;
}//parseContentForToolCallsAndSendResults



std::string LlmInterface::stripThinkingContent(const std::string& content) {
  std::string result = content;
  
  // Remove <think>...</think> blocks (case insensitive, supports nested and multiline)
  // Use [\s\S] instead of . to match any character including newlines
  std::regex thinkRegex("<think[^>]*>[\\s\\S]*?</think>", 
    std::regex_constants::icase | std::regex_constants::ECMAScript);
  
  // Keep removing think blocks until no more are found (handles nested cases)
  std::string prevResult;
  do {
    prevResult = result;
    result = std::regex_replace(result, thinkRegex, "");
  } while (result != prevResult);
  
  // Clean up extra whitespace that may be left after removing think blocks
  // Replace multiple newlines with at most two newlines
  std::regex multiNewlineRegex("\n\n\n+");
  result = std::regex_replace(result, multiNewlineRegex, "\n\n");
  
  // Trim leading and trailing whitespace
  result = std::regex_replace(result, std::regex("^\\s+"), "");
  result = std::regex_replace(result, std::regex("\\s+$"), "");
  
  return result;
}

std::pair<std::string, std::string> LlmInterface::extractThinkingAndContent(const std::string& content) {
  std::string cleanContent = content;
  std::string thinkingContent;
  
  // Extract <think>...</think> blocks (case insensitive, supports nested and multiline)
  std::regex thinkRegex("<think[^>]*>([\\s\\S]*?)</think>", 
    std::regex_constants::icase | std::regex_constants::ECMAScript);
  
  std::smatch match;
  std::string::const_iterator searchStart(content.cbegin());
  
  // Collect all thinking content from properly formatted <think>...</think> blocks
  while (std::regex_search(searchStart, content.cend(), match, thinkRegex)) {
    if (!thinkingContent.empty()) {
      thinkingContent += "\n";
    }
    thinkingContent += match[1].str();
    searchStart = match.suffix().first;
  }
  
  // Check for orphaned </think> tags (closing tag without opening tag)
  std::regex orphanedCloseRegex("</think>", std::regex_constants::icase);
  std::smatch orphanedMatch;
  if (std::regex_search(content, orphanedMatch, orphanedCloseRegex)) {
    // Find the position of the first </think> tag
    size_t closePos = orphanedMatch.position();
    
    // Extract all text before the </think> tag as thinking content
    std::string orphanedThinking = content.substr(0, closePos);
    
    // Only add if we found some content and it's not just whitespace
    if (!orphanedThinking.empty() && orphanedThinking.find_first_not_of("\\s\\t\\n\\r") != std::string::npos) {
      if (!thinkingContent.empty()) {
        thinkingContent += "\n";
      }
      thinkingContent += orphanedThinking;
      
      // Remove the orphaned thinking content and </think> tag from clean content
      cleanContent = content.substr(closePos + 7); // 7 is length of "</think>"
    }
  } else {
    // No orphaned tags, use normal stripping
    cleanContent = stripThinkingContent(content);
  }
  
  return {cleanContent, thinkingContent};
}

std::string LlmInterface::getSystemPromptForAgent( const AgentType agentType ) const
{
  if( !m_config || !m_config->llmApi.enabled )
    return "";

  // Search for agent in config
  for( const LlmConfig::AgentConfig &agent : m_config->agents )
  {
    if( agent.type == agentType )
      return agent.systemPrompt;
  }

  assert( 0 );
  throw std::runtime_error( "Failed to find agent config for agent " + agentTypeToString(agentType) );
  
  return "";
}//getSystemPromptForAgent(...)


int LlmInterface::invokeSubAgent( std::shared_ptr<LlmConversationStart> sub_agent_convo )
{
  assert( sub_agent_convo );
  if( !sub_agent_convo )
    throw std::logic_error( "LlmInterface::invokeSubAgent called with null conversation" );
  
  assert( sub_agent_convo->agent_type != AgentType::MainAgent );
  assert( !sub_agent_convo->content.empty() );
  assert( sub_agent_convo->responses.empty() );
  
  cout << "=== Invoking sub-agent: " << agentTypeToString(sub_agent_convo->agent_type) << " (async with pause) ===" << endl;
  cout << "Context/Task: " << sub_agent_convo->content.substr(0, 100) << (sub_agent_convo->content.length() > 100 ? "..." : "") << endl;
  cout << "Invoke tool call ID: " << sub_agent_convo->conversationId << endl;

  // Build API request for this specific agent - since the combined message is already in the
  const json requestJson = buildMessagesArray( sub_agent_convo );
  
  // Make tracked API call
  const int requestId = m_nextRequestId++;

  // Create pending request with sub-agent info
  PendingRequest pending;
  pending.requestId = requestId;
  pending.conversation = sub_agent_convo;
  pending.isSubAgentRequest = true;

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  pending.requestJson = requestJson;
#endif

  m_pendingRequests[requestId] = pending;

  cout << "Main agent will pause until sub-agent completes" << endl;

  // Make the actual API call
  makeApiCallWithId(requestJson, requestId);

  return requestId;
}//invokeSubAgent(...)


nlohmann::json LlmInterface::buildMessagesArray( const std::shared_ptr<LlmConversationStart> &convo )
{
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );
  
  assert( convo );
  if( !convo )
    throw std::logic_error( "LlmInterface::buildMessagesArray: null conversation history passed in." );
  
  json request;
  request["model"] = m_config->llmApi.model;
  
  // Use max_completion_tokens for newer OpenAI models, max_tokens for others
  string modelName = m_config->llmApi.model;
  if (modelName.find("gpt-4") != string::npos || modelName.find("gpt-3.5") != string::npos || 
      modelName.find("o1") != string::npos || modelName.find("gpt-5") != string::npos) {
    request["max_completion_tokens"] = m_config->llmApi.maxTokens;
  } else {
    request["max_tokens"] = m_config->llmApi.maxTokens;
  }
  
  json messages = json::array();

  // Add system prompt (use MainAgent's prompt)
  const string systemPrompt = getSystemPromptForAgent( convo->agent_type );
  if( !systemPrompt.empty() )
  {
    json systemMsg;
    systemMsg["role"] = "system";
    systemMsg["content"] = systemPrompt;
    messages.push_back(systemMsg);
  }
  
  // Add conversation history
  if( convo->agent_type == AgentType::MainAgent )
  {
    const vector<shared_ptr<LlmConversationStart>> &conversations = m_history->getConversations();
    if( !conversations.empty() )
    {
      cout << "=== Including " << conversations.size() << " history messages in request ===" << endl;
      
      for( const shared_ptr<LlmConversationStart> &previous_conversation : conversations )
      {
        assert( previous_conversation );
        LlmConversationHistory::addConversationToLlmApiHistory( *previous_conversation, messages );
        
        if( convo == previous_conversation )
          break;
      }
    }else
    {
      cout << "=== No history to include ===" << endl;
    }
  }else
  {
    // For sub-agents, we don't include full chat history of the user - just the task context
    LlmConversationHistory::addConversationToLlmApiHistory( *convo, messages );
  }
  
  request["messages"] = messages;
  
  // Add tools
  json tools = json::array();
  
  const map<string, LlmTools::SharedTool> agentTools = m_tool_registry->getToolsForAgent(convo->agent_type);
  for( const auto &[name, tool] : agentTools )
  {
    json toolDef;
    toolDef["type"] = "function";
    toolDef["function"]["name"] = tool.name;
    toolDef["function"]["description"] = m_tool_registry->getDescriptionForAgent(tool.name, convo->agent_type);
    
    // Functions used to include "userSession" argument for the MCP server (now thats added in by the MCP server) - but we'll make sure of this here for the moment until we verify they have been totally removed.
    nlohmann::json par_schema = tool.parameters_schema;
    assert( !par_schema.contains("properties") || !par_schema["properties"].contains("userSession") );
    if( par_schema.contains("properties") && par_schema["properties"].contains("userSession") )
      par_schema["properties"].erase( "userSession" );
    
    toolDef["function"]["parameters"] = par_schema;
    tools.push_back(toolDef);
  }//for( const auto &[name, tool] : m_tool_registry->getTools() )
  
  assert( !tools.empty() );
  if( !tools.empty() )
  {
    request["tools"] = tools;
    request["tool_choice"] = "auto";
  }
  
  return request;
}//nlohmann::json LlmInterface::buildMessagesArray( convo )


#ifdef USE_JS_BRIDGE_FOR_LLM
void LlmInterface::setupJavaScriptBridge() {
  auto app = Wt::WApplication::instance();
  if (!app) {
    cout << "Warning: No WApplication instance for JavaScript bridge setup" << endl;
    return;
  }
  
  cout << "Setting up JavaScript bridge for HTTP requests..." << endl;
  
  // Set up the JavaScript function to handle HTTP requests
  string jsCode = R"(
    window.llmHttpRequest = function(endpoint, requestJsonString, bearerToken, requestId) {
      var startTime = Date.now();
      var requestObj = null;

      try {
        requestObj = JSON.parse(requestJsonString);
      } catch(e) {
        console.error('Failed to parse request JSON:', e);
      }

      // Log request diagnostics
      console.log('=== LLM Request ID', requestId, '===');
      console.log('Endpoint:', endpoint);
      console.log('Submit time:', new Date().toISOString());
      console.log('Request size:', requestJsonString.length, 'bytes');

      if (requestObj) {
        console.log('Model:', requestObj.model || 'not specified');
        console.log('Messages count:', requestObj.messages ? requestObj.messages.length : 0);
        console.log('Tools count:', requestObj.tools ? requestObj.tools.length : 0);

        // Calculate total conversation tokens (rough estimate)
        var totalChars = 0;
        if (requestObj.messages) {
          requestObj.messages.forEach(function(msg) {
            if (msg.content && typeof msg.content === 'string') {
              totalChars += msg.content.length;
            }
          });
        }
        console.log('Estimated input chars:', totalChars);

        // Log truncated request JSON
        var jsonStr = JSON.stringify(requestObj);
        if (jsonStr.length > 80) {
          var truncated = jsonStr.substring(0, 70) + '...' + jsonStr.substring(jsonStr.length - 10);
          console.log('Request JSON (truncated):', truncated);
        } else {
          console.log('Request JSON:', jsonStr);
        }
        // Uncomment to see full JSON:
        // console.log('Request JSON (full):', JSON.stringify(requestObj, null, 2));
      } else {
        var truncated = requestJsonString.length > 80
          ? requestJsonString.substring(0, 70) + '...' + requestJsonString.substring(requestJsonString.length - 10)
          : requestJsonString;
        console.log('Request JSON (raw, truncated):', truncated);
      }

      var headers = {
        'Content-Type': 'application/json'
      };

      if (bearerToken && bearerToken.trim() !== '') {
        headers['Authorization'] = 'Bearer ' + bearerToken;
      }

      // Create AbortController for timeout handling
      var controller = new AbortController();
      var timeoutMs = 120000; // 2 minutes
      var timeoutId = setTimeout(function() {
        var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
        console.error('=== LLM Request TIMEOUT (ID ' + requestId + ') ===');
        console.error('Elapsed time: ' + elapsed + 's');
        console.error('Timeout limit: ' + (timeoutMs / 1000) + 's');
        console.error('Endpoint: ' + endpoint);
        controller.abort();
      }, timeoutMs);

      console.log('Timeout set for:', (timeoutMs / 1000) + 's');
      console.log('Initiating fetch...');

      fetch(endpoint, {
        method: 'POST',
        headers: headers,
        body: requestJsonString,
        signal: controller.signal
      })
      .then(function(response) {
        var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
        console.log('=== LLM Response received (ID ' + requestId + ') ===');
        console.log('Status:', response.status, response.statusText);
        console.log('Elapsed time:', elapsed + 's');
        return response.text();
      })
      .then(function(responseText) {
        var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
        console.log('Response parsed, size:', responseText.length, 'bytes');
        console.log('Total round-trip time:', elapsed + 's');

        // Clear the timeout since we got a response
        clearTimeout(timeoutId);

        // Call back to C++ with response and request ID as separate parameters
        if (window.llmResponseCallback) {
          window.llmResponseCallback(responseText, requestId);
        }
      })
      .catch(function(error) {
        var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
        console.error('=== LLM Request ERROR (ID ' + requestId + ') ===');
        console.error('Error type:', error.name);
        console.error('Error message:', error.message);
        console.error('Elapsed time:', elapsed + 's');
        console.error('Is timeout?', error.name === 'AbortError');

        // Clear the timeout since we got an error
        clearTimeout(timeoutId);

        var errorResponse = JSON.stringify({
          error: {
            message: 'HTTP request failed: ' + error.message,
            type: error.name === 'AbortError' ? 'timeout_error' : 'network_error',
            elapsed_seconds: parseFloat(elapsed),
            request_id: requestId
          }
        });
        if (window.llmResponseCallback) {
          window.llmResponseCallback(errorResponse, requestId);
        }
      });
    };
    
    console.log('LLM JavaScript bridge initialized');
  )";
  
  app->doJavaScript(jsCode);
  
  // Set up the response callback using JSignal to emit signal to C++
  string callbackJs = 
    "window.llmResponseCallback = function(response, requestId) { "
    //"console.log('Emitting signal to C++ with response length:', response.length, 'requestId:', requestId); "
    "" + m_responseSignal->createCall("response", "requestId") + ";"
    "};";
  
  app->doJavaScript(callbackJs);
  
  cout << "JavaScript bridge setup complete" << endl;
}

void LlmInterface::handleJavaScriptResponse(std::string response, int requestId)
{
  cout << "\n\n=== Recieved requestId=" << requestId << " ===" << endl;
  SpecUtils::trim( response );
  
  std::string responsePreview = response;
  SpecUtils::ireplace_all( responsePreview, "\n", " ");
  SpecUtils::ireplace_all( responsePreview, "\r", "");
  if( responsePreview.length() > 300 )
    responsePreview = responsePreview.substr(0, 300) + "...";
  cout << "Response: " << responsePreview << endl;

  try
  {
    // Find and remove the pending request
    PendingRequest pendingRequest;
    if( m_pendingRequests.find(requestId) != m_pendingRequests.end() )
    {
      pendingRequest = std::move(m_pendingRequests[requestId]);
      m_pendingRequests.erase(requestId);
    }else
    {
      cerr << "Got response that didnt have pending request: requestId=" << requestId << ", response='" << response << "'" << endl;
      assert( 0 );
      return;
    }
    
    std::shared_ptr<LlmConversationStart> convo = pendingRequest.conversation.lock();
    if( !convo )
    {
      cerr << "For JavaScript response, found original request, but conversation is nullptr, so ending this conversation." << endl;
      assert( convo );
      
      return;
    }//if( !convo )
    
    // Check for errors first
    cout << "--- about to parse response to json --- " << endl;
    json responseJson = json::parse(response);
    cout << "--- done parsing response to json --- " << endl;


    if( responseJson.contains("error") && !responseJson["error"].is_null() )
    {
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
      string now_str = SpecUtils::to_iso_string( now );
      const string::size_type period_pos = now_str.find('.');
      if( period_pos != string::npos )
        now_str = now_str.substr(0,period_pos);
        
      string debug_name = "llm_request_with_error_id" + std::to_string(requestId) + "_" + now_str + ".json";
      string debug_result = "llm_result_with_error_id" + std::to_string(requestId) + "_" + now_str + ".json";
#ifdef _WIN32
      const std::wstring wdebug_name = SpecUtils::convert_from_utf8_to_utf16(debug_name);
      std::ofstream output_request_json( wdebug_name.c_str(), ios::binary | ios::out );
      const std::wstring wdebug_result = SpecUtils::convert_from_utf8_to_utf16(debug_result);
      std::ofstream output_result_json( wdebug_result.c_str(), ios::binary | ios::out );
#else
      std::ofstream output_request_json( debug_name.c_str(), ios::binary | ios::out);
      std::ofstream output_result_json( debug_result.c_str(), ios::binary | ios::out );
#endif
      if( output_request_json )
        output_request_json << pendingRequest.requestJson.dump(2);
      if( output_result_json )
        output_result_json << response;
      cout << "\nLLM request error: wrote request input and output to '"
      << debug_name << "', and '" << debug_result << "', respectively."
      << endl << endl;
#endif //#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      
      string errorMsg = "LLM API Error: " + responseJson["error"].dump(2);
      cout << errorMsg << endl;
      
      // Add error to conversation history
      if (m_history)
        m_history->addErrorMessage( errorMsg, convo );
      
      // Signal that a response was received (even if it's an error)
      // Only emit if no pending requests (this is the final response)
      if (m_pendingRequests.empty()) {
        cout << "Error response - no pending requests, emitting signal" << endl;
        m_responseReceived.emit();
      } else {
        cout << "Error response - still have " << m_pendingRequests.size() << " pending requests, not emitting signal yet" << endl;
      }
      return;
    }
    
    if( pendingRequest.isSubAgentRequest )
      cout << "=== Processing sub-agent response for: " << agentTypeToString(convo->agent_type) << " ===" << endl;
    
    handleApiResponse( response, convo, requestId );
  }catch( const json::parse_error &e )
  {
    cout << "Failed to parse LLM response as JSON: " << e.what() << endl;
    cout << "Raw response: " << response << endl;
  }catch( const std::exception &e )
  {
    cout << "Error processing LLM response: " << e.what() << endl;
  }
}
#endif // USE_JS_BRIDGE_FOR_LLM

int LlmInterface::makeTrackedApiCall( const nlohmann::json& requestJson,
                                      std::shared_ptr<LlmConversationStart> convo )
{
  assert( convo );
  
  const int requestId = m_nextRequestId++;  

  // Store the pending request
  PendingRequest pending;
  pending.requestId = requestId;
  pending.conversation = convo;
  pending.isSubAgentRequest = (convo->agent_type != AgentType::MainAgent);
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  pending.requestJson = requestJson;
#endif

  m_pendingRequests[requestId] = pending;
  
  // Make the call with request ID tracking
  makeApiCallWithId(requestJson, requestId);
  
  return requestId;
}

int LlmInterface::sendToolResultsToLLM( std::shared_ptr<LlmConversationStart> convo )
{
  // The history already contains the tool calls and results, so we just need to
  // build a new request with the current conversation state
  json followupRequest = buildMessagesArray( convo ); // System-generated followup
  
  // Make tracked API call to get LLM's response to the tool results
  const int requestId = makeTrackedApiCall( followupRequest, convo );
  
  return requestId;
}

#ifdef USE_WT_HTTP_FOR_LLM
void LlmInterface::handleWtHttpClientResponse(boost::system::error_code err, const Wt::Http::Message& response) {
  cout << "=== Wt HTTP Client Response Received ===" << endl;
  
  // Get the request ID from our stored current request
  int requestId = m_currentRequestId;
  cout << "Processing response for request ID: " << requestId << endl;
  
  if (err) {
    cout << "HTTP error: " << err.message() << endl;
    handleWtHttpError(requestId, err.message());
    return;
  }
  
  try {
    cout << "Response status: " << response.status() << endl;
    
    // Check for successful HTTP status codes
    if (response.status() < 200 || response.status() >= 300) {
      string errorMsg = "HTTP error status: " + std::to_string(response.status());
      cout << errorMsg << endl;
      handleWtHttpError(requestId, errorMsg);
      return;
    }
    
    // Get response body
    string responseBody = response.body();
    cout << "Response length: " << responseBody.length() << " characters" << endl;
    
    // Find and remove the pending request
    PendingRequest pendingRequest;
    bool foundPending = false;
    if (m_pendingRequests.find(requestId) != m_pendingRequests.end()) {
      pendingRequest = m_pendingRequests[requestId];
      foundPending = true;
      m_pendingRequests.erase(requestId);
      cout << "Matched to pending request: " << (pendingRequest.isToolResultFollowup ? "tool followup" : "original") << endl;
    } else {
      cout << "Warning: No pending request found for ID " << requestId << endl;
    }
    
    // Process the LLM response using unified handler
    handleApiResponse(responseBody);
    
  } catch (const std::exception& e) {
    cout << "Exception processing Wt HTTP response: " << e.what() << endl;
    handleWtHttpError(requestId, std::string("Exception: ") + e.what());
  }
}

void LlmInterface::handleWtHttpResponse(int requestId, const std::string& responseBody, int statusCode) {
  cout << "=== Wt HTTP Response Received for ID " << requestId << " ===" << endl;
  cout << "Response status: " << statusCode << endl;
  
  try {
    cout << "Response length: " << responseBody.length() << " characters" << endl;
    
    // Find and remove the pending request
    PendingRequest pendingRequest;
    bool foundPending = false;
    if (m_pendingRequests.find(requestId) != m_pendingRequests.end()) {
      pendingRequest = m_pendingRequests[requestId];
      foundPending = true;
      m_pendingRequests.erase(requestId);
      cout << "Matched to pending request: " << (pendingRequest.isToolResultFollowup ? "tool followup" : "original") << endl;
    } else {
      cout << "Warning: No pending request found for ID " << requestId << endl;
    }
    
    // Check HTTP status
    if (statusCode != 200) {
      cout << "HTTP Error: Status " << statusCode << endl;
      return;
    }
    
    // Process the LLM response
    handleApiResponse(responseBody);
    
  } catch (const std::exception& e) {
    cout << "Error processing Wt HTTP response: " << e.what() << endl;
  }
}

void LlmInterface::handleWtHttpError(int requestId, const std::string& error) {
  cout << "=== Wt HTTP Error for ID " << requestId << " ===" << endl;
  cout << "Error: " << error << endl;
  
  // Find and remove the pending request
  if (m_pendingRequests.find(requestId) != m_pendingRequests.end()) {
    m_pendingRequests.erase(requestId);
    cout << "Removed pending request for ID " << requestId << endl;
  }
  
  // Add error to conversation history
  if (m_history) {
    m_history->addErrorMessage("HTTP Error: " + error, m_currentConversationId);
  }
  
  // Signal that a response was received (even if it's an error)
  // Only emit if no pending requests (this is the final response)
  if (m_pendingRequests.empty()) {
    cout << "HTTP Error - no pending requests, emitting signal" << endl;
    m_responseReceived.emit();
  } else {
    cout << "HTTP Error - still have " << m_pendingRequests.size() << " pending requests, not emitting signal yet" << endl;
  }
}
#endif // USE_WT_HTTP_FOR_LLM


