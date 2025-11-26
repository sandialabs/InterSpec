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
struct LlmInteraction;
class LlmConversationHistory;

enum class AgentType : int;

namespace LlmTools {
  class ToolRegistry;
}

class LlmToolRequest;
class LlmToolResults;

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
class LlmInterface : public Wt::Signals::trackable
{
public:
  /** Construct LLM interface for the given InterSpec instance.
   @param interspec The InterSpec instance this interface belongs to
   */
  explicit LlmInterface( InterSpec *interspec, const std::shared_ptr<const LlmConfig> &config );
  
  /** Destructor - implementation in .cpp to handle incomplete types */
  ~LlmInterface();
  
  /** Returns the tool registry.
   
   Will be a valid pointer.
   */
  std::shared_ptr<const LlmTools::ToolRegistry> toolRegistry();
  
  /** Send a user message to the LLM.
   @param message The user's message/question
   
   @returns The `LlmInteraction` created for this message.
   */
  std::shared_ptr<LlmInteraction> sendUserMessage(const std::string& message);
  
  /** Send a system-generated message to the LLM.
   
   This is used when the application automatically asks the LLM for assistance,
   rather than the user explicitly asking a question.
   
   @param message The system-generated query
   */
  void sendSystemMessage(const std::string& message);
  
  /** Test chat history recording and reconstruction without calling an actual LLM.
   
   This simulates various conversation scenarios including tool calls to verify
   that history is being properly recorded and can be reconstructed correctly.
   */
  
  /** Get conversation history (may be null if no history yet) */
  std::shared_ptr<LlmConversationHistory> getHistory() const;
  
  /** Set conversation history (typically loaded from SpecMeas) */
  void setHistory(std::shared_ptr<LlmConversationHistory> history);
  
  /** Check if LLM interface is properly configured and ready to use */
  bool isConfigured() const;
  
  /** Signal emitted when a new response is received from the LLM */
  Wt::Signal<>& conversationFinished();

  /** Signal emitted when an error response is received from the LLM */
  Wt::Signal<>& responseError();
  
  /** Check if a specific request ID is still pending */
  bool isRequestPending(int requestId) const;

  /** Invoke a sub-agent to handle a specific task (async).
   @param sub_agent_convo The conversation to send to the LLM to start the sub-agent
   @return Request ID for the sub-agent invocation
   */
  int invokeSubAgent( std::shared_ptr<LlmInteraction> sub_agent_convo );

  /** Build the messages array for a conversation.

   This is exposed publicly to support retry functionality from LlmToolGui.

   @param convo The conversation to build messages for
   @return JSON object ready to send to LLM API
   */
  nlohmann::json buildMessagesArray( const std::shared_ptr<LlmInteraction> &convo );

  /** Make a tracked API call with the given JSON request.

   This is exposed publicly to support retry functionality from LlmToolGui.

   @param requestJson The JSON request to send
   @param convo The conversation this request belongs to
   @return Pair of (requestId, request content as string)
   */
  std::pair<int,std::string> makeTrackedApiCall( const nlohmann::json& requestJson,
                         std::shared_ptr<LlmInteraction> convo );

  /** Send tool results back to LLM for processing.

   This is exposed publicly to support retry functionality from LlmToolGui.

   @param convo The conversation to send tool results for
   @return Request ID for the API call
   */
  int sendToolResultsToLLM( std::shared_ptr<LlmInteraction> convo );

#ifdef USE_JS_BRIDGE_FOR_LLM
  /** JavaScript callback to handle LLM response */
  void handleJavaScriptResponse(std::string response, int requestId);
#endif

private:
  void emitConversationFinished(){ conversationFinished().emit(); }
  
  InterSpec* m_interspec;
  const std::shared_ptr<const LlmConfig> m_config;
  const std::shared_ptr<const LlmTools::ToolRegistry> m_tool_registry;
  std::shared_ptr<LlmConversationHistory> m_history;
  
  Wt::Signal<> m_conversationFinished; // Signal emitted when succesful final response from LLM is recieved.
  Wt::Signal<> m_responseError;    // Signal emitted when error responses are received
  
#ifdef USE_JS_BRIDGE_FOR_LLM
  std::unique_ptr<Wt::JSignal<std::string, int>> m_responseSignal; // For JavaScript bridge (response, requestId)
#endif

#ifdef USE_WT_HTTP_FOR_LLM
  std::unique_ptr<Wt::Http::Client> m_httpClient; // Wt HTTP client for native C++ requests
  int m_currentRequestId; // Current request ID for Wt HTTP client correlation
#endif
  
  // Request tracking
  int m_nextRequestId;

  struct PendingRequest
  {
    int requestId;
    std::weak_ptr<LlmInteraction> conversation;

    // Sub-agent support
    bool isSubAgentRequest = false;       // True if this request is from a sub-agent (we can probably get rid of this variable)

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
    nlohmann::json requestJson;
#endif
  };//struct PendingRequest
  
  std::map<int, PendingRequest> m_pendingRequests;

  // Deferred tool results (for sub-agent invocations that need to pause main agent)
  struct DeferredToolResult {
    std::string conversationId;
    std::vector<std::string> toolCallIds;  // Tool call IDs to send back when sub-agent completes (does not include sub-agent calls)
    std::vector<int> subAgentToolCallIds;  // The invoke_sub_agent tool call ID to update with summary
  };
  std::map<int, DeferredToolResult> m_deferredToolResults; // Key is sub-agent requestId
  
  /** Make an API call with request ID tracking

   @returns The request body content (e.g., the JSON as a string).
   */
  std::string makeApiCallWithId(const nlohmann::json& requestJson, int requestId);

  /** Handle response from LLM API */
  void handleApiResponse( const std::string &response, const std::shared_ptr<LlmInteraction> &convo, const int requestId );
  
  /** Execute tool calls requested by the LLM, and sends the LLM back a response with the results.
   
   Returns the number of tool calls processed.
   */
  std::pair<std::shared_ptr<LlmToolRequest>, std::shared_ptr<LlmToolResults>>
  executeToolCallsAndSendResults( const nlohmann::json& toolCalls,
                                         const std::shared_ptr<LlmInteraction> &convo,
                                         const int requestId,
                                         const std::string &rawResponseContent,
                                         std::optional<size_t> promptTokens = std::nullopt,
                                         std::optional<size_t> completionTokens = std::nullopt );
  
  /** Parse text content for tool call requests (for models that don't support structured tool calls)
   
   Returns the number of tool calls processed.
   */
  std::pair<std::shared_ptr<LlmToolRequest>, std::shared_ptr<LlmToolResults>>
  parseContentForToolCallsAndSendResults( const std::string &content,
                                                const std::shared_ptr<LlmInteraction> &convo,
                                                const int requestId,
                                                const std::string &rawResponseContent );
  
  /** Strip <think>...</think> content from LLM responses */
  static std::string stripThinkingContent(const std::string& content);
  
  /** Extract thinking content and clean content from LLM responses */
  static std::pair<std::string, std::string> extractThinkingAndContent(const std::string& content);

  /** Get the system prompt for a specific agent from config
   */
  std::string getSystemPromptForAgent( const AgentType agentType ) const;
  
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


#if( PERFORM_DEVELOPER_CHECKS )
// Test namespace to expose static functions for unit testing
namespace LlmInterfaceTests
{
  /** Exposed for unit testing only - parse JSON with lenient error handling */
  nlohmann::json lenientlyParseJson( const std::string &jsonStr );

  /** Exposed for unit testing only - sanitize JSON string before parsing */
  std::string sanitizeJsonString( const std::string &jsonStr );

  /** Exposed for unit testing only - repair structurally incomplete JSON */
  std::string repairIncompleteJson( const std::string &jsonStr,
                                    std::string *repairLog = nullptr );
}//namespace LlmInterfaceTests
#endif


#endif // LlmInterface_h 
