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
#include <iostream>
#include <stdexcept>

#include <Wt/WApplication>
#include <Wt/WResource>
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
    m_nextRequestId(1),
    m_currentConversationId(""),
    m_currentAgentName("MainAgent")
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

void LlmInterface::sendUserMessage(const std::string& message) {
  if (!isConfigured()) {
    throw std::runtime_error("LLM interface is not properly configured");
  }
  
  cout << "User message: " << message << endl;
  
  // Generate a new conversation ID for each user message
  m_currentConversationId = "conv_" + std::to_string(chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now().time_since_epoch()).count());
  cout << "Starting new conversation with ID: " << m_currentConversationId << endl;
  
  // Add to history FIRST to ensure it's there before any async responses
  cout << "About to add user message to history..." << endl;
  m_history->addUserMessage(message, m_currentConversationId);
  cout << "Added user message to history. History now has " << m_history->getConversations().size() << " conversations" << endl;
  
  // Debug: Show what's actually in history
  cout << "Current history contents:" << endl;
  for (size_t i = 0; i < m_history->getConversations().size(); ++i) {
    const auto& conv = m_history->getConversations()[i];
    string typeStr;
    switch (conv.type) {
      case LlmConversationStart::Type::User: typeStr = "user"; break;
      case LlmConversationStart::Type::System: typeStr = "system"; break;
    }
    cout << "  " << i << ". " << typeStr << ": " << conv.content.substr(0, 50) << "..." << " (responses: " << conv.responses.size() << ")" << endl;
  }
  
  // Build API request (buildMessagesArray will include the message from history + current message)
  json requestJson = buildMessagesArray(message, false);
  
  // Make tracked API call
  int requestId = makeTrackedApiCall(requestJson, message, false);
  
  cout << "Sent user message with request ID: " << requestId << endl;
}

void LlmInterface::sendSystemMessage(const std::string& message) {
  if (!isConfigured()) {
    throw std::runtime_error("LLM interface is not properly configured");
  }
  
  cout << "System message: " << message << endl;
  
  // Build API request (marked as system-generated)
  json requestJson = buildMessagesArray(message, true);
  
  // Make tracked API call
  int requestId = makeTrackedApiCall(requestJson, message, false);
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

namespace {
  // Map from WApplication pointer to current LlmInterface instance
  // This allows tools to access the active LlmInterface during execution
  std::map<Wt::WApplication*, LlmInterface*> s_currentInstances;
}

LlmInterface* LlmInterface::getCurrentInstance()
{
  auto app = Wt::WApplication::instance();
  if( !app )
    return nullptr;

  const auto iter = s_currentInstances.find(app);
  if( iter != s_currentInstances.end() )
    return iter->second;

  return nullptr;
}//getCurrentInstance()


void LlmInterface::setCurrentInstance( LlmInterface *instance )
{
  auto app = Wt::WApplication::instance();
  if( !app )
    return;

  if( instance )
    s_currentInstances[app] = instance;
  else
    s_currentInstances.erase(app);
}//setCurrentInstance(...)


Wt::Signal<>& LlmInterface::responseReceived() {
  return m_responseReceived;
}

bool LlmInterface::isRequestPending(int requestId) const {
  return m_pendingRequests.find(requestId) != m_pendingRequests.end();
}

void LlmInterface::makeApiCall(const nlohmann::json& requestJson) {
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );
  
  cout << "=== Making LLM API Call ===" << endl;
  cout << "Endpoint: " << m_config->llmApi.apiEndpoint << endl;
  cout << "Request JSON:" << endl;
  nlohmann::json debugJson = requestJson;
  if( debugJson.contains("messages") && debugJson["messages"].is_array() && debugJson["messages"].size() && debugJson["messages"][0].is_object() )
    debugJson["messages"][0]["content"] = "...system prompt not shown..."; //Erase System prompt
  if( debugJson.contains("tools") )
    debugJson.erase("tools");
  if( debugJson.contains("tool_choice") )
    debugJson.erase("tool_choice");
  cout << debugJson.dump(2) << endl;
  cout << "=========================" << endl;
  
#ifdef USE_JS_BRIDGE_FOR_LLM
  // Use JavaScript bridge to make the HTTP request
  string requestStr = requestJson.dump();
  
  string jsCall = 
    "var llmRequestData = " + requestStr + ";\n"
    "window.llmHttpRequest('" + m_config->llmApi.apiEndpoint + "', JSON.stringify(llmRequestData), '" + 
    m_config->llmApi.bearerToken + "');";
  
  // Execute JavaScript to make the HTTP request
  auto app = Wt::WApplication::instance();
  if (app) {
    app->doJavaScript(jsCall);
    cout << "JavaScript HTTP request initiated..." << endl;
  } else {
    cout << "Error: No WApplication instance available for JavaScript bridge" << endl;
    
    // Fallback to simulation for testing
    string simulatedResponse = R"({
      "choices": [{
        "message": {
          "role": "assistant",
          "content": "This is a simulated response. JavaScript bridge not available."
        }
      }]
    })";
    handleApiResponse(simulatedResponse);
  }
#endif

#ifdef USE_WT_HTTP_FOR_LLM
  // Use Wt HTTP client to make the request
  cout << "Using Wt HTTP client for request..." << endl;
  
  // For now, fall back to simulation since we need request ID tracking
  // This should be called via makeApiCallWithId instead
  string simulatedResponse = R"({
    "choices": [{
      "message": {
        "role": "assistant",
        "content": "This is a simulated response. Direct makeApiCall not supported with Wt HTTP - use makeApiCallWithId."
      }
    }]
  })";
  handleApiResponse(simulatedResponse);
#endif
}



void LlmInterface::makeApiCallWithId(const nlohmann::json& requestJson, int requestId)
{
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );
  
  cout << "=== Making LLM API Call with ID " << requestId << " ===" << endl;
  cout << "Endpoint: " << m_config->llmApi.apiEndpoint << endl;
  cout << "Request JSON:" << endl;
  nlohmann::json debugJson = requestJson;
  if( debugJson.contains("messages") && debugJson["messages"].is_array() && debugJson["messages"].size() && debugJson["messages"][0].is_object() )
    debugJson["messages"][0]["content"] = "...system prompt not shown..."; //Erase System prompt
  if( debugJson.contains("tools") )
    debugJson.erase("tools");
  if( debugJson.contains("tool_choice") )
    debugJson.erase("tool_choice");
  cout << debugJson.dump(2) << endl;
  cout << "=========================" << endl;
  
#ifdef USE_JS_BRIDGE_FOR_LLM
  // Use JavaScript bridge to make the HTTP request with request ID
  string requestStr = requestJson.dump();
  
  string jsCall = 
    "var llmRequestData = " + requestStr + ";\n"
    "window.llmHttpRequest('" + m_config->llmApi.apiEndpoint + "', JSON.stringify(llmRequestData), '" + 
    m_config->llmApi.bearerToken + "', " + std::to_string(requestId) + ");";
  
  // Execute JavaScript to make the HTTP request
  auto app = Wt::WApplication::instance();
  if (app) {
    app->doJavaScript(jsCall);
    cout << "JavaScript HTTP request with ID " << requestId << " initiated..." << endl;
  } else {
    cout << "Error: No WApplication instance available for JavaScript bridge" << endl;
    
    // Fallback to simulation for testing
    string simulatedResponse = R"({
      "choices": [{
        "message": {
          "role": "assistant",
          "content": "This is a simulated response. JavaScript bridge not available."
        }
      }]
    })";
    handleJavaScriptResponse(simulatedResponse, requestId);
  }
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
    
    // Set request body to JSON
    string requestBody = requestJson.dump();
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

void LlmInterface::handleApiResponse(const std::string& response) {
  try {
    json responseJson = json::parse(response);
    
    // Parse and accumulate token usage information if available
    if (responseJson.contains("usage") && m_history && !m_currentConversationId.empty()) {
      const auto& usage = responseJson["usage"];
      
      std::optional<int> promptTokens, completionTokens, totalTokens;
      if (usage.contains("prompt_tokens") && usage["prompt_tokens"].is_number())
        promptTokens = usage["prompt_tokens"].get<int>();
      if (usage.contains("completion_tokens") && usage["completion_tokens"].is_number())
        completionTokens = usage["completion_tokens"].get<int>();
      if (usage.contains("total_tokens") && usage["total_tokens"].is_number())
        totalTokens = usage["total_tokens"].get<int>();
      
      // Accumulate token usage for this conversation
      m_history->addTokenUsage(m_currentConversationId, promptTokens, completionTokens, totalTokens);
      
      if (completionTokens.has_value()) {
        cout << "=== Token Usage This Call ===" << endl;
        cout << "Prompt tokens: " << (promptTokens.has_value() ? std::to_string(promptTokens.value()) : "N/A") << endl;
        cout << "Completion tokens: " << completionTokens.value() << endl;
        cout << "Total tokens: " << (totalTokens.has_value() ? std::to_string(totalTokens.value()) : "N/A") << endl;
        cout << "=============================" << endl;
      }
    }
    
    if (responseJson.contains("choices") && !responseJson["choices"].empty()) {
      json choice = responseJson["choices"][0];
      if (choice.contains("message")) {
        json message = choice["message"];
        string role = message.value("role", "");
        string content;
        if( message.contains("content") && message["content"].is_string() )
          content = message["content"];

        if (role == "assistant") {
          // Extract thinking content and clean content
          auto [cleanContent, thinkingContent] = extractThinkingAndContent(content);

          // Add assistant message to history with thinking content and current agent name
          m_history->addAssistantMessageWithThinking(cleanContent, thinkingContent, m_currentConversationId, m_currentAgentName);

          // Handle structured tool calls first (OpenAI format)
          if (message.contains("tool_calls")) {
            executeToolCalls(message["tool_calls"]);
          } else {
            // Parse content for text-based tool requests (use cleaned content)
            parseContentForToolCalls(cleanContent);
          }
        }
      }
    }
    
  } catch (const std::exception& e) {
    cout << "Error parsing LLM response: " << e.what()
    << "\n\tresponse="
    << response << endl << endl;
  }
  
  // Only emit signal if there are no pending requests (i.e., this is the final response)
  if (m_pendingRequests.empty()) {
    cout << "No pending requests, emitting response received signal" << endl;
    m_responseReceived.emit();
  } else {
    cout << "Still have " << m_pendingRequests.size() << " pending requests, not emitting signal yet" << endl;
  }
}

void LlmInterface::executeToolCalls( const nlohmann::json &toolCalls )
{
  cout << "Executing " << toolCalls.size() << " tool calls" << endl;

  // Track executed tool calls for follow-up
  std::vector<std::string> executedToolCallIds;
  bool hasSubAgentInvocation = false;
  int subAgentRequestId = -1;

  for( const nlohmann::json &toolCall : toolCalls )
  {
    const string callId = toolCall.value("id", "");

    try
    {
      const string toolName = toolCall["function"]["name"];
      json arguments = json::parse(toolCall["function"]["arguments"].get<string>());

      cout << "Calling tool: " << toolName << " with ID: " << callId << endl;

      // Add tool call to history with conversation ID, invocation ID, and current agent name
      m_history->addToolCall(toolName, m_currentConversationId, callId, arguments, m_currentAgentName);

      // Set current instance so tools can access this LlmInterface (e.g., for sub-agent invocation)
      setCurrentInstance(this);

      // Check if this is an invoke_* tool (sub-agent invocation) - handle specially
    
      if( SpecUtils::istarts_with(toolName, "invoke_") && (toolName.length() > 7) )
      {
        // Extract agent name from tool name (e.g., "invoke_NuclideId" -> "NuclideId")
        const string agentName = toolName.substr(7);
        const string context = arguments.at("context").get<string>();
        const string task = arguments.at("task").get<string>();

        cout << "Detected sub-agent invocation: " << agentName << " via tool: " << toolName << endl;

        // Invoke sub-agent (returns request ID)
        subAgentRequestId = invokeSubAgent(agentName, context, task, callId);

        hasSubAgentInvocation = true;

        // Add placeholder result that will be updated when sub-agent completes
        json result;
        result["status"] = "pending";
        result["message"] = "Sub-agent " + agentName + " is processing (will update when complete)";
        result["requestId"] = subAgentRequestId;
        m_history->addToolResult(m_currentConversationId, callId, result, m_currentAgentName);

        cout << "Sub-agent invocation deferred - will resume main agent when complete" << endl;
      }else
      {
        // Normal tool execution
        json result = m_tool_registry->executeTool(toolName, arguments, m_interspec);

        // Add result to history with conversation ID, invocation ID, and current agent name
        m_history->addToolResult(m_currentConversationId, callId, result, m_currentAgentName);
      }

      // Clear current instance
      setCurrentInstance(nullptr);
    }catch( const std::exception &e )
    {
      cout << "Tool execution error: " << e.what() << endl;

      json result;
      result["error"] = "Tool call failed: " + string(e.what());
      m_history->addToolResult(m_currentConversationId, callId, result, m_currentAgentName);
    }//try / catch

    // Track this tool call for follow-up
    executedToolCallIds.push_back(callId);
  }

  // If we have a sub-agent invocation, defer sending results until sub-agent completes
  if( hasSubAgentInvocation )
  {
    cout << "Deferring tool results - sub-agent must complete first (requestId: " << subAgentRequestId << ")" << endl;

    DeferredToolResult deferred;
    deferred.conversationId = m_currentConversationId;
    deferred.toolCallIds = executedToolCallIds;
    deferred.subAgentToolCallId = executedToolCallIds.back(); // The invoke_sub_agent call should be last

    m_deferredToolResults[subAgentRequestId] = deferred;

    // Don't send results yet - wait for sub-agent
  }else
  {
    // No sub-agent - send results immediately
    if (!executedToolCallIds.empty()) {
      cout << "Sending tool results back to LLM for " << executedToolCallIds.size() << " executed tools" << endl;
      sendToolResultsToLLM(executedToolCallIds);
    }
  }
}

void LlmInterface::parseContentForToolCalls(const std::string& content) {
  cout << "Parsing content for text-based tool calls..." << endl;
  
  // Track executed tool calls for follow-up
  std::vector<std::string> executedToolCallIds;
  
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
  for (const auto& pattern : toolCallPatterns) {
    std::sregex_iterator iter(content.begin(), content.end(), pattern);
    std::sregex_iterator end;
    
    for (; iter != end; ++iter) {
      const std::smatch& match = *iter;
      
      // Generate a unique invocation ID for this tool call
      string invocationId = "text_call_" + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                                                                                                std::chrono::system_clock::now().time_since_epoch()).count());
      
      if (match.size() >= 2) {
        try {
          string toolName;
          json arguments;
          
          if (match.size() == 2) {
            // Single capture group - JSON object with name and arguments fields
            string toolCallJson = match[1].str();
            cout << "Found JSON tool call: " << toolCallJson << endl;
            
            json toolCallObj = json::parse(toolCallJson);
            if (toolCallObj.contains("name") && toolCallObj.contains("arguments")) {
              toolName = toolCallObj["name"];
              arguments = toolCallObj["arguments"];
            } else {
              string name = toolCallObj.contains("name") ? toolCallObj["name"] : "missing_tool_name";
              string args = toolCallObj.contains("arguments") ? toolCallObj["arguments"] : "missing_arguments";
              
              m_history->addToolCall( name, m_currentConversationId, invocationId, args);
              throw std::runtime_error( "Invalid tool call JSON format - missing name or arguments" );
            }
          } else if (match.size() >= 3) {
            // Two capture groups - tool name and arguments separately
            toolName = match[1].str();
            string argumentsStr = match[2].str();
            
            // Parse arguments JSON
            if (argumentsStr.empty() || argumentsStr == "{}") {
              arguments = json::object();
            } else {
              arguments = json::parse(argumentsStr);
            }
          }
          
          cout << "Found text-based tool call: " << toolName << " with args: " << arguments.dump() << endl;


          // Add tool call to history with conversation ID, invocation ID, and current agent name
          m_history->addToolCall(toolName, m_currentConversationId, invocationId, arguments, m_currentAgentName);
          cout << "Added tool call to history. History now has " << m_history->getConversations().size() << " conversations" << endl;

          // Set current instance so tools can access this LlmInterface (e.g., for sub-agent invocation)
          setCurrentInstance(this);

          // Execute the tool
          json result = m_tool_registry->executeTool(toolName, arguments, m_interspec);

          // Clear current instance
          setCurrentInstance(nullptr);

          // Add result to history with conversation ID, invocation ID, and current agent name
          m_history->addToolResult(m_currentConversationId, invocationId, result, m_currentAgentName);
          cout << "Added tool result to history. History now has " << m_history->getConversations().size() << " conversations" << endl;

          // Track this tool call for follow-up
          executedToolCallIds.push_back(invocationId);

          cout << "Tool " << toolName << " executed successfully" << endl;

        } catch (const std::exception& e) {
          cout << "Error parsing/executing text-based tool call: " << e.what() << endl;

          json result;
          result["error"] = "Tool call failed: " + string(e.what());
          m_history->addToolResult(m_currentConversationId, invocationId, result, m_currentAgentName);
        }
      }
    }
  }
  
  // If we executed any tools, send the results back to the LLM for processing
  if (!executedToolCallIds.empty()) {
    cout << "Sending tool results back to LLM for " << executedToolCallIds.size() << " executed tools" << endl;
    sendToolResultsToLLM(executedToolCallIds);
  }
  
  cout << "=== Response Processing Complete ===" << endl;
}

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

std::string LlmInterface::getSystemPromptForAgent( const std::string &agentName ) const
{
  if( !m_config || !m_config->llmApi.enabled )
    return "";

  // Search for agent in config
  for( const LlmConfig::AgentConfig &agent : m_config->agents )
  {
    if( agent.name == agentName )
      return agent.systemPrompt;
  }

  return "";
}//getSystemPromptForAgent(...)


nlohmann::json LlmInterface::buildMessagesArrayForAgent( const std::string &agentName, const std::string &userMessage, bool isSystemGenerated )
{
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );

  json request;
  request["model"] = m_config->llmApi.model;

  // Use max_completion_tokens for newer OpenAI models, max_tokens for others
  const string modelName = m_config->llmApi.model;
  if( modelName.find("gpt-4") != string::npos || modelName.find("gpt-3.5") != string::npos ||
      modelName.find("o1") != string::npos || modelName.find("gpt-5") != string::npos )
  {
    request["max_completion_tokens"] = m_config->llmApi.maxTokens;
  }else
  {
    request["max_tokens"] = m_config->llmApi.maxTokens;
  }

  json messages = json::array();

  // Add system prompt for this specific agent
  const string systemPrompt = getSystemPromptForAgent(agentName);
  if( !systemPrompt.empty() )
  {
    json systemMsg;
    systemMsg["role"] = "system";
    systemMsg["content"] = systemPrompt;
    messages.push_back(systemMsg);
  }

  // For sub-agents, we typically don't include full history - just the task context
  // Main agent uses buildMessagesArray() which includes history

  // Add current user message (the task/context for sub-agent)
  if( !userMessage.empty() )
  {
    json userMsg;
    userMsg["role"] = "user";
    userMsg["content"] = userMessage;
    messages.push_back(userMsg);
  }

  request["messages"] = messages;

  // Add tools available for this agent
  const AgentType agentType = stringToAgentType(agentName);
  const std::map<std::string, LlmTools::SharedTool> agentTools = m_tool_registry->getToolsForAgent(agentType);
  json tools = json::array();

  for( const auto &[toolName, tool] : agentTools )
  {
    json toolDef;
    toolDef["type"] = "function";
    toolDef["function"]["name"] = tool.name;

    // Use role-specific description if available
    const string description = m_tool_registry->getDescriptionForAgent(toolName, agentType);
    toolDef["function"]["description"] = description;

    // Most functions include a "userSession" argument for the MCP server - but we dont need that here since we are in a Wt session
    nlohmann::json par_schema = tool.parameters_schema;
    if( par_schema.contains("properties") && par_schema["properties"].contains("userSession") )
      par_schema["properties"].erase( "userSession" );

    toolDef["function"]["parameters"] = par_schema;
    tools.push_back(toolDef);
  }

  if( !tools.empty() )
  {
    request["tools"] = tools;
    request["tool_choice"] = "auto";
  }

  return request;
}//buildMessagesArrayForAgent(...)


int LlmInterface::invokeSubAgent( const std::string &agentName, const std::string &context, const std::string &task,
                                  const std::string &invokeToolCallId )
{
  cout << "=== Invoking sub-agent: " << agentName << " (async with pause) ===" << endl;
  cout << "Context: " << context.substr(0, 100) << (context.length() > 100 ? "..." : "") << endl;
  cout << "Task: " << task << endl;
  cout << "Invoke tool call ID: " << invokeToolCallId << endl;

  // Build a combined message with context and task
  const string combinedMessage = "Context:\n" + context + "\n\nTask:\n" + task;

  // Build API request for this specific agent
  const json requestJson = buildMessagesArrayForAgent(agentName, combinedMessage, false);

  // Save current state (will be restored when sub-agent completes)
  const string savedAgentName = m_currentAgentName;
  const string savedConversationId = m_currentConversationId;

  // Set current agent and create a new conversation ID for sub-agent
  m_currentAgentName = agentName;
  m_currentConversationId = "subagent_" + agentName + "_" + std::to_string(
    chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count()
  );

  // Add sub-agent context message to history
  m_history->addUserMessage(combinedMessage, m_currentConversationId);

  // Make tracked API call
  const int requestId = m_nextRequestId++;

  // Create pending request with sub-agent info
  PendingRequest pending;
  pending.requestId = requestId;
  pending.originalUserMessage = combinedMessage;
  pending.isToolResultFollowup = false;
  pending.isSubAgentRequest = true;
  pending.subAgentName = agentName;
  pending.savedAgentName = savedAgentName;
  pending.savedConversationId = savedConversationId;
  // No callback needed - we'll handle resumption directly

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  pending.requestJson = requestJson;
#endif

  m_pendingRequests[requestId] = pending;

  cout << "Sub-agent " << agentName << " request ID: " << requestId << " (conversation: " << m_currentConversationId << ")" << endl;
  cout << "Main agent will pause until sub-agent completes" << endl;

  // Make the actual API call
  makeApiCallWithId(requestJson, requestId);

  // Note: m_currentAgentName will be restored when the response is received
  // For now, keep it set to the sub-agent so any immediate processing is correctly attributed

  return requestId;
}//invokeSubAgent(...)


nlohmann::json LlmInterface::buildMessagesArray(const std::string& userMessage, bool isSystemGenerated)
{
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );
  
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
  const string systemPrompt = getSystemPromptForAgent("MainAgent");
  if( !systemPrompt.empty() )
  {
    json systemMsg;
    systemMsg["role"] = "system";
    systemMsg["content"] = systemPrompt;
    messages.push_back(systemMsg);
  }
  

  
          // Add conversation history
        if (!m_history->isEmpty()) {
          json historyMessages = m_history->toApiFormat();
          cout << "=== Including " << historyMessages.size() << " history messages in request ===" << endl;
          
          for (size_t i = 0; i < historyMessages.size(); ++i) {
            const auto& msg = historyMessages[i];
            
            cout << "  " << i << ". " << msg["role"].get<string>() << ": " 
                 << (msg.contains("content") ? msg["content"].get<string>().substr(0, 50) + "..." : "tool_call") << endl;
            messages.push_back(msg);
          }
          cout << "=== End history messages ===" << endl;
        } else {
          cout << "=== No history to include ===" << endl;
        }
  
      // Add current user message (only if not empty, not system-generated, and not already in history)
    if (!userMessage.empty() && !isSystemGenerated) {
      // Check if this message is already the last message in history to prevent duplication
      const auto& historyConversations = m_history->getConversations();
              bool isAlreadyLastMessage = false;
        if (!historyConversations.empty()) {
          const auto& lastConversation = historyConversations.back();
          isAlreadyLastMessage = (lastConversation.type == LlmConversationStart::Type::User && 
                                lastConversation.content == userMessage);
        }
      
      if (!isAlreadyLastMessage) {
        json userMsg;
        userMsg["role"] = "user";
        userMsg["content"] = userMessage;
        messages.push_back(userMsg);
      }
    }
  
  request["messages"] = messages;
  
  // Add tools
  json tools = json::array();
  for (const auto& [name, tool] : m_tool_registry->getTools()) {
    json toolDef;
    toolDef["type"] = "function";
    toolDef["function"]["name"] = tool.name;
    toolDef["function"]["description"] = tool.description;
    
    // Most functions include a "userSession" argument for the MCP server - but we dont need that here since we are in a Wt session
    nlohmann::json par_schema = tool.parameters_schema;
    if( par_schema.contains("properties") && par_schema["properties"].contains("userSession") )
      par_schema["properties"].erase( "userSession" );
    
    toolDef["function"]["parameters"] = par_schema;
    tools.push_back(toolDef);
  }
  
  if (!tools.empty()) {
    request["tools"] = tools;
    request["tool_choice"] = "auto";
  }
  
  return request;
}

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
      //console.log('LLM HTTP Request to:', endpoint, 'For requestID', requestId);
      //console.log('Request data:', requestJsonString);
      
      var headers = {
        'Content-Type': 'application/json'
      };
      
      if (bearerToken && bearerToken.trim() !== '') {
        headers['Authorization'] = 'Bearer ' + bearerToken;
      }
      
      // Create AbortController for timeout handling
      var controller = new AbortController();
      var timeoutId = setTimeout(function() {
        console.log('LLM Request timeout after 2 minutes');
        controller.abort();
      }, 120000); // 2 minutes = 120 seconds
      
      fetch(endpoint, {
        method: 'POST',
        headers: headers,
        body: requestJsonString,
        signal: controller.signal
      })
      .then(function(response) {
        //console.log('LLM Response status:', response.status);
        return response.text();
      })
      .then(function(responseText) {
        //console.log('LLM Response:', responseText);
        //console.log( 'Got LLM Response text', responseText ); 
        
        // Clear the timeout since we got a response
        clearTimeout(timeoutId);
        
        // Call back to C++ with response and request ID as separate parameters
        if (window.llmResponseCallback) {
          window.llmResponseCallback(responseText, requestId);
        }
      })
      .catch(function(error) {
        console.error('LLM Request error:', error);
        
        // Clear the timeout since we got an error
        clearTimeout(timeoutId);
        
        var errorResponse = JSON.stringify({
          error: {
            message: 'HTTP request failed: ' + error.message,
            type: error.name === 'AbortError' ? 'timeout_error' : 'network_error'
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
  try
  {
    // Find and remove the pending request
    PendingRequest pendingRequest;
    bool foundPending = false;
    if( m_pendingRequests.find(requestId) != m_pendingRequests.end() )
    {
      pendingRequest = std::move(m_pendingRequests[requestId]);
      foundPending = true;
      m_pendingRequests.erase(requestId);
    }
    
    // Check for errors first
    json responseJson = json::parse(response);
    if( responseJson.contains("error") && !responseJson["error"].is_null() )
    {
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      if( foundPending )
      {
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
      }//if( foundPending )
#endif //#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      
      string errorMsg = "LLM API Error: " + responseJson["error"].dump(2);
      cout << errorMsg << endl;
      
      // Add error to conversation history
      if (m_history) {
        m_history->addErrorMessage(errorMsg, m_currentConversationId);
      }
      
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
    
    // Check if this is a sub-agent request
    if( foundPending && pendingRequest.isSubAgentRequest )
    {
      cout << "=== Processing sub-agent response for: " << pendingRequest.subAgentName << " ===" << endl;

      // Process the response normally first (this adds messages to history with correct agent name)
      handleApiResponse(response);

      // Extract the final assistant message as the sub-agent summary
      string subAgentSummary;
      if( responseJson.contains("choices") && !responseJson["choices"].empty() )
      {
        const json &choice = responseJson["choices"][0];
        if( choice.contains("message") && choice["message"].contains("content") )
        {
          const string content = choice["message"]["content"].get<string>();
          // Extract clean content (strip thinking tags if present)
          auto [cleanContent, thinkingContent] = extractThinkingAndContent(content);
          subAgentSummary = cleanContent;
        }
      }

      cout << "Sub-agent summary: " << subAgentSummary.substr(0, 100) << (subAgentSummary.length() > 100 ? "..." : "") << endl;

      // Restore the saved agent name and conversation ID
      m_currentAgentName = pendingRequest.savedAgentName;
      m_currentConversationId = pendingRequest.savedConversationId;

      cout << "Restored agent: " << m_currentAgentName << ", conversation: " << m_currentConversationId << endl;

      // Check if there are deferred tool results waiting for this sub-agent
      const auto deferredIter = m_deferredToolResults.find(requestId);
      if( deferredIter != m_deferredToolResults.end() )
      {
        const DeferredToolResult &deferred = deferredIter->second;

        cout << "=== Resuming main agent with sub-agent summary ===" << endl;

        // Update the invoke_sub_agent tool result with the actual summary
        // Find the tool result in history and update it
        if( m_history )
        {
          std::vector<LlmConversationStart> &conversations = m_history->getConversations();
          for( LlmConversationStart &conv : conversations )
          {
            if( conv.conversationId == deferred.conversationId )
            {
              for( LlmConversationResponse &resp : conv.responses )
              {
                if( resp.type == LlmConversationResponse::Type::ToolResult &&
                   resp.invocationId == deferred.subAgentToolCallId )
                {
                  // Update the tool result with the real summary
                  json updatedResult;
                  updatedResult["status"] = "completed";
                  updatedResult["summary"] = subAgentSummary;
                  updatedResult["agentName"] = pendingRequest.subAgentName;
                  resp.content = updatedResult.dump();

                  cout << "Updated invoke_sub_agent tool result with summary" << endl;
                  break;
                }
              }
              break;
            }
          }
        }

        // Now send all the deferred tool results back to the main agent LLM
        cout << "Sending " << deferred.toolCallIds.size() << " deferred tool results back to main agent" << endl;
        sendToolResultsToLLM(deferred.toolCallIds);

        // Clean up deferred results
        m_deferredToolResults.erase(deferredIter);
      }else
      {
        cout << "Warning: No deferred tool results found for sub-agent request " << requestId << endl;
      }

      cout << "=== Sub-agent processing complete ===" << endl;
    }else
    {
      // Normal (non-sub-agent) request - use the enhanced handleApiResponse method for consistent processing
      handleApiResponse(response);
    }

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

int LlmInterface::makeTrackedApiCall(const nlohmann::json& requestJson, const std::string& originalMessage, bool isToolFollowup) {
  int requestId = m_nextRequestId++;
  
  // Store the pending request
  PendingRequest pending;
  pending.requestId = requestId;
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  pending.requestJson = requestJson;
#endif
  pending.originalUserMessage = originalMessage;
  pending.isToolResultFollowup = isToolFollowup;
  m_pendingRequests[requestId] = pending;
  

  
  // Make the call with request ID tracking
  makeApiCallWithId(requestJson, requestId);
  
  return requestId;
}

void LlmInterface::sendToolResultsToLLM(const std::vector<std::string>& toolCallIds) {
  // The history already contains the tool calls and results, so we just need to 
  // build a new request with the current conversation state
  json followupRequest = buildMessagesArray("", true); // System-generated followup
  
  // Make tracked API call to get LLM's response to the tool results
  int requestId = makeTrackedApiCall(followupRequest, "", true);
  
  // Update the pending request to track which tool calls this is following up on
  if (m_pendingRequests.find(requestId) != m_pendingRequests.end()) {
    m_pendingRequests[requestId].toolCallIds = toolCallIds;
  }
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


