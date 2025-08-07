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

#include <iostream>
#include <memory>
#include <stdexcept>
#include <regex>
#include <chrono>

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

#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmConversationHistory.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/InterSpec.h"
#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

static_assert( USE_LLM_INTERFACE, "This file shouldnt be compiled unless USE_LLM_INTERFACE is true" );

using namespace std;
using json = nlohmann::json;



LlmInterface::LlmInterface(InterSpec* interspec) 
  : m_interspec(interspec),
    m_config(std::make_unique<LlmConfig>(LlmConfig::load())),
    m_history(std::make_shared<LlmConversationHistory>()),
#ifdef USE_JS_BRIDGE_FOR_LLM
    m_responseSignal(std::make_unique<Wt::JSignal<std::string, int>>(interspec, "llmResponse")),
#endif
#ifdef USE_WT_HTTP_FOR_LLM
    m_httpClient(std::make_unique<Wt::Http::Client>(interspec)),
    m_currentRequestId(0),
#endif
    m_nextRequestId(1),
    m_currentConversationId("")
{
  if (!m_interspec)
    throw std::runtime_error("InterSpec instance cannot be null");
  
  // Get session ID through WApplication
  string sessionId = "unknown";
  if (Wt::WApplication* app = Wt::WApplication::instance()) {
    sessionId = app->sessionId();
  }
  cout << "LlmInterface created for session: " << sessionId << endl;
  
  // Register default tools
  LlmTools::ToolRegistry::instance().registerDefaultTools();
  
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

void LlmInterface::testConnection() {
  cout << "=== LLM Interface Test ===" << endl;
  cout << "Config API endpoint: " << m_config->llmApi.apiEndpoint << endl;
  cout << "Config model: " << m_config->llmApi.model << endl;
  cout << "Bearer token configured: " << (!m_config->llmApi.bearerToken.empty() ? "Yes" : "No") << endl;
  cout << "Is configured: " << (isConfigured() ? "Yes" : "No") << endl;
  
  // Test tool registry
  cout << "Available tools: " << LlmTools::ToolRegistry::instance().getTools().size() << endl;
  for (const auto& [name, tool] : LlmTools::ToolRegistry::instance().getTools()) {
    cout << "  - " << name << ": " << tool.description << endl;
  }
  
  // Test a simple tool call
  try {
    json testParams;
    json result = LlmTools::ToolRegistry::instance().executeTool("test_tool", testParams, m_interspec);
    cout << "Test tool result: " << result.dump(2) << endl;
  } catch (const std::exception& e) {
    cout << "Test tool error: " << e.what() << endl;
  }
  
  // Test with a real API call if configured
  if (isConfigured()) {
    cout << "Making test API call to LLM..." << endl;
    // Use sendUserMessage to properly add to history
    sendUserMessage("What peaks do you see in my foreground spectrum? Please analyze the detected peaks.");
  } else {
    cout << "LLM not configured - skipping API call test" << endl;
    // Create a test API request
    json testRequest = buildMessagesArray("Hello, this is a test message", false);
    cout << "Test API request would be:" << endl;
    cout << testRequest.dump(2) << endl;
  }
  

  
  cout << "=========================" << endl;
}

void LlmInterface::reloadConfig() {
  m_config = std::make_unique<LlmConfig>(LlmConfig::load());
  cout << "LLM configuration reloaded" << endl;
}

const LlmConfig& LlmInterface::getConfig() const {
  return *m_config;
}

std::shared_ptr<LlmConversationHistory> LlmInterface::getHistory() const {
  return m_history;
}

void LlmInterface::setHistory(std::shared_ptr<LlmConversationHistory> history) {
  m_history = history ? history : std::make_shared<LlmConversationHistory>();
}

bool LlmInterface::isConfigured() const {
  return !m_config->llmApi.apiEndpoint.empty();
}

Wt::Signal<>& LlmInterface::responseReceived() {
  return m_responseReceived;
}

bool LlmInterface::isRequestPending(int requestId) const {
  return m_pendingRequests.find(requestId) != m_pendingRequests.end();
}

void LlmInterface::makeApiCall(const nlohmann::json& requestJson) {
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



void LlmInterface::makeApiCallWithId(const nlohmann::json& requestJson, int requestId) {
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
    
    if (responseJson.contains("choices") && !responseJson["choices"].empty()) {
      json choice = responseJson["choices"][0];
      if (choice.contains("message")) {
        json message = choice["message"];
        string role = message.value("role", "");
        string content = message.value("content", "");
        
                if (role == "assistant") {
          // Extract thinking content and clean content
          auto [cleanContent, thinkingContent] = extractThinkingAndContent(content);
                  
                  cout
                  << "=== Start Cleaned Response Content ===" << endl
                  << cleanContent
                  << "\n=== End Cleaned Response Content   ===" << endl
                  << endl;

          // Add assistant message to history with thinking content
          m_history->addAssistantMessageWithThinking(cleanContent, thinkingContent, m_currentConversationId);
          
          // Handle structured tool calls first (OpenAI format)
          if (message.contains("tool_calls")) {
            cout << "Found structured tool_calls" << endl;
            executeToolCalls(message["tool_calls"]);
          } else {
            // Parse content for text-based tool requests (use cleaned content)
            parseContentForToolCalls(cleanContent);
          }
        }
      }
    }
    
  } catch (const std::exception& e) {
    cout << "Error parsing LLM response: " << e.what() << endl;
  }
  
  // Only emit signal if there are no pending requests (i.e., this is the final response)
  if (m_pendingRequests.empty()) {
    cout << "No pending requests, emitting response received signal" << endl;
    m_responseReceived.emit();
  } else {
    cout << "Still have " << m_pendingRequests.size() << " pending requests, not emitting signal yet" << endl;
  }
}

void LlmInterface::executeToolCalls(const nlohmann::json& toolCalls) {
  cout << "Executing " << toolCalls.size() << " tool calls" << endl;
  
  // Track executed tool calls for follow-up
  std::vector<std::string> executedToolCallIds;
  
  for (const auto& toolCall : toolCalls) {
    const string callId = toolCall.value("id", "");
    
    try {
      const string toolName = toolCall["function"]["name"];
      json arguments = json::parse(toolCall["function"]["arguments"].get<string>());
      
      cout << "Calling tool: " << toolName << " with ID: " << callId << endl;
      
      // Add tool call to history with conversation ID and invocation ID
      m_history->addToolCall(toolName, m_currentConversationId, callId, arguments);
      
      // Execute the tool
      json result = LlmTools::ToolRegistry::instance().executeTool(toolName, arguments, m_interspec);
      
      // Add result to history with conversation ID and invocation ID
      m_history->addToolResult(m_currentConversationId, callId, result);
    } catch (const std::exception& e) {
      cout << "Tool execution error: " << e.what() << endl;
      
      json result;
      result["error"] = "Tool call failed: " + string(e.what());
      m_history->addToolResult(m_currentConversationId, callId, result);
    }
    
    // Track this tool call for follow-up
    executedToolCallIds.push_back(callId);
  }
  
  // If we executed any tools, send the results back to the LLM for processing
  if (!executedToolCallIds.empty()) {
    cout << "Sending tool results back to LLM for " << executedToolCallIds.size() << " executed tools" << endl;
    sendToolResultsToLLM(executedToolCallIds);
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
          
          
          // Add tool call to history with conversation ID and invocation ID
          m_history->addToolCall(toolName, m_currentConversationId, invocationId, arguments);
          cout << "Added tool call to history. History now has " << m_history->getConversations().size() << " conversations" << endl;
          
          // Execute the tool
          json result = LlmTools::ToolRegistry::instance().executeTool(toolName, arguments, m_interspec);
          
          // Add result to history with conversation ID and invocation ID
          m_history->addToolResult(m_currentConversationId, invocationId, result);
          cout << "Added tool result to history. History now has " << m_history->getConversations().size() << " conversations" << endl;
          
          // Track this tool call for follow-up
          executedToolCallIds.push_back(invocationId);
          
          cout << "Tool " << toolName << " executed successfully" << endl;
          
        } catch (const std::exception& e) {
          cout << "Error parsing/executing text-based tool call: " << e.what() << endl;
          
          json result;
          result["error"] = "Tool call failed: " + string(e.what());
          m_history->addToolResult(m_currentConversationId, invocationId, result);
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

nlohmann::json LlmInterface::buildMessagesArray(const std::string& userMessage, bool isSystemGenerated) {
  json request;
  request["model"] = m_config->llmApi.model;
  request["max_tokens"] = m_config->llmApi.maxTokens;
  
  json messages = json::array();
  
  // Add system prompt
  json systemMsg;
  systemMsg["role"] = "system";
  systemMsg["content"] = m_config->llmApi.systemPrompt;
  messages.push_back(systemMsg);
  

  
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
  for (const auto& [name, tool] : LlmTools::ToolRegistry::instance().getTools()) {
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
      console.log('LLM HTTP Request to:', endpoint, 'For requestID', requestId);
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
        console.log('LLM Response status:', response.status);
        return response.text();
      })
      .then(function(responseText) {
        //console.log('LLM Response:', responseText);
        console.log( 'Got LLM Response text' ); 
        
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
    "console.log('Emitting signal to C++ with response length:', response.length, 'requestId:', requestId); "
    "" + m_responseSignal->createCall("response", "requestId") + ";"
    "};";
  
  app->doJavaScript(callbackJs);
  
  cout << "JavaScript bridge setup complete" << endl;
}

void LlmInterface::handleJavaScriptResponse(std::string response, int requestId) {
  try {
    // Find and remove the pending request
    PendingRequest pendingRequest;
    bool foundPending = false;
    if (m_pendingRequests.find(requestId) != m_pendingRequests.end()) {
      pendingRequest = m_pendingRequests[requestId];
      foundPending = true;
      m_pendingRequests.erase(requestId);
    }
    
    // Check for errors first
    json responseJson = json::parse(response);
    if (responseJson.contains("error")) {
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
    
    // Use the enhanced handleApiResponse method for consistent processing
    handleApiResponse(response);
    
  } catch (const json::parse_error& e) {
    cout << "Failed to parse LLM response as JSON: " << e.what() << endl;
    cout << "Raw response: " << response << endl;
  } catch (const std::exception& e) {
    cout << "Error processing LLM response: " << e.what() << endl;
  }
}
#endif // USE_JS_BRIDGE_FOR_LLM

int LlmInterface::makeTrackedApiCall(const nlohmann::json& requestJson, const std::string& originalMessage, bool isToolFollowup) {
  int requestId = m_nextRequestId++;
  
  // Store the pending request
  PendingRequest pending;
  pending.requestId = requestId;
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


