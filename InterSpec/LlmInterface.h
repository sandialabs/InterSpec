#ifndef LlmInterface_h
#define LlmInterface_h
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

#include <string>
#include <memory>
#include <functional>

#include <Wt/WContainerWidget>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

// Choose LLM HTTP implementation method
//#define USE_WT_HTTP_FOR_LLM  // Use Wt::Http::Client for HTTPS requests (needs further debugging)
#define USE_JS_BRIDGE_FOR_LLM   // Use JavaScript bridge for HTTPS requests (default)

// Ensure only one implementation is selected
#if defined(USE_WT_HTTP_FOR_LLM) && defined(USE_JS_BRIDGE_FOR_LLM)
#error "Cannot define both USE_WT_HTTP_FOR_LLM and USE_JS_BRIDGE_FOR_LLM"
#endif

#if !defined(USE_WT_HTTP_FOR_LLM) && !defined(USE_JS_BRIDGE_FOR_LLM)
#error "Must define either USE_WT_HTTP_FOR_LLM or USE_JS_BRIDGE_FOR_LLM"
#endif

#ifdef USE_WT_HTTP_FOR_LLM
#include <boost/system/error_code.hpp>
#include <Wt/Http/Client>
#include <Wt/Http/Message>
#endif


// Forward declarations
class InterSpec;
class LlmConfig;
class LlmConversationHistory;

namespace Wt {
  class WResource;

#ifdef USE_JS_BRIDGE_FOR_LLM
  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6> 
  class JSignal;
#endif
}

/** Main LLM interface class that handles communication with OpenAI-compatible API endpoints.
 
 This class manages:
 - API calls to LLM endpoints using JavaScript bridge for HTTPS
 - Tool calling and execution
 - Conversation history management  
 - Integration with InterSpec session
 */
class LlmInterface {
public:
  /** Construct LLM interface for the given InterSpec instance.
   @param interspec The InterSpec instance this interface belongs to
   */
  explicit LlmInterface(InterSpec* interspec);
  
  /** Destructor - implementation in .cpp to handle incomplete types */
  ~LlmInterface();
  
  /** Send a user message to the LLM.
   @param message The user's message/question
   */
  void sendUserMessage(const std::string& message);
  
  /** Send a system-generated message to the LLM.
   
   This is used when the application automatically asks the LLM for assistance,
   rather than the user explicitly asking a question.
   
   @param message The system-generated query
   */
  void sendSystemMessage(const std::string& message);
  
  /** Test the LLM connection with a simple request.
   
   This will print request/response to stdout for debugging.
   */
  void testConnection();
  
  /** Reload configuration from XML files */
  void reloadConfig();
  
  /** Get current configuration */
  const LlmConfig& getConfig() const;
  
  /** Get conversation history (may be null if no history yet) */
  std::shared_ptr<LlmConversationHistory> getHistory() const;
  
  /** Set conversation history (typically loaded from SpecMeas) */
  void setHistory(std::shared_ptr<LlmConversationHistory> history);
  
  /** Check if LLM interface is properly configured and ready to use */
  bool isConfigured() const;
  
#ifdef USE_JS_BRIDGE_FOR_LLM
  /** JavaScript callback to handle LLM response */
  void handleJavaScriptResponse(std::string response, int requestId);
#endif

private:
  InterSpec* m_interspec;
  std::unique_ptr<LlmConfig> m_config;
  std::shared_ptr<LlmConversationHistory> m_history;
  
#ifdef USE_JS_BRIDGE_FOR_LLM
  std::unique_ptr<Wt::JSignal<std::string, int>> m_responseSignal; // For JavaScript bridge (response, requestId)
#endif

#ifdef USE_WT_HTTP_FOR_LLM
  std::unique_ptr<Wt::Http::Client> m_httpClient; // Wt HTTP client for native C++ requests
  int m_currentRequestId; // Current request ID for Wt HTTP client correlation
#endif
  
  // Request tracking
  int m_nextRequestId;
  std::string m_currentToolCallId; // Track the current tool call ID for message association
  struct PendingRequest {
    int requestId;
    std::string originalUserMessage;
    bool isToolResultFollowup;
    std::vector<std::string> toolCallIds; // For tracking which tool calls this request contains
  };
  std::map<int, PendingRequest> m_pendingRequests;
  
  /** Make an API call to the LLM endpoint */
  void makeApiCall(const nlohmann::json& requestJson);
  
  /** Make an API call with request ID tracking */
  void makeApiCallWithId(const nlohmann::json& requestJson, int requestId);
  
  /** Make an API call with request tracking */
  int makeTrackedApiCall(const nlohmann::json& requestJson, const std::string& originalMessage = "", bool isToolFollowup = false);
  
  /** Send tool results back to LLM for processing */
  void sendToolResultsToLLM(const std::vector<std::string>& toolCallIds);
  
  /** Handle response from LLM API */
  void handleApiResponse(const std::string& response);
  
  /** Execute tool calls requested by the LLM */
  void executeToolCalls(const nlohmann::json& toolCalls);
  
  /** Parse text content for tool call requests (for models that don't support structured tool calls) */
  void parseContentForToolCalls(const std::string& content);
  
  /** Strip <think>...</think> content from LLM responses */
  std::string stripThinkingContent(const std::string& content);
  
  /** Build the messages array for the API request including history and system prompt */
  nlohmann::json buildMessagesArray(const std::string& userMessage, bool isSystemGenerated = false);
  
#ifdef USE_JS_BRIDGE_FOR_LLM
  /** Set up the JavaScript bridge for making HTTPS requests */
  void setupJavaScriptBridge();
#endif

#ifdef USE_WT_HTTP_FOR_LLM
  /** Handle HTTP response from Wt::Http::Client */
  void handleWtHttpResponse(int requestId, const std::string& responseBody, int statusCode = 200);
  
  /** Handle HTTP error from Wt::Http::Client */
  void handleWtHttpError(int requestId, const std::string& error);
  
  /** Handle HTTP client done signal (boost::system::error_code, Http::Message) */
  void handleWtHttpClientResponse(boost::system::error_code err, const Wt::Http::Message& response);
#endif
};

// Legacy namespace for backward compatibility
namespace LlmInterfaceNS
{
  /** A simple "hello world" type function for testing purposes.
   
   Returns a greeting message indicating the LlmInterface module is working.
   
   @return A string containing a hello world message.
   */
  InterSpec_API std::string hello_world();
  
} // namespace LlmInterfaceNS

#endif // LlmInterface_h 
