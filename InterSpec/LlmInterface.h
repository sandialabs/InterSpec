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

#include <map>
#include <chrono>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include <functional>

#include <Wt/WContainerWidget>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "InterSpec/LlmConversationHistory.h"

#if( USE_NATIVE_HTTP_CLIENT )
#include "InterSpec/NativeHttpClient.h"
#endif

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );


// Forward declarations
class InterSpec;
class LlmConfig;
class LlmApiProtocol;
struct LlmInteraction;
class LlmConversationHistory;

enum class AgentType : int;

namespace LlmTools {
  class ToolRegistry;
}

class LlmToolRequest;
class LlmToolResults;

namespace Wt {
  class WTimer;
  class WResource;

  template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
  class JSignal;
}

/** Main LLM interface class that handles communication with OpenAI-compatible API endpoints.
 
 This class manages:
 - API calls to LLM endpoints using JavaScript bridge for HTTPS
 - Tool calling and execution
 - Conversation history management  
 - Integration with InterSpec session
 */
class LlmInterface : public Wt::Signals::trackable,
                     public std::enable_shared_from_this<LlmInterface>
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
   @param images Optional image content to include with the message (e.g., from drag-and-drop)

   @returns The `LlmInteraction` created for this message.
   */
  std::shared_ptr<LlmInteraction> sendUserMessage( const std::string &message,
                                                    std::vector<LlmToolCall::ImageContent> images = {} );
  
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
  
  /** Get the LLM configuration (for building API requests externally). */
  std::shared_ptr<const LlmConfig> config() const;

  /** Check if LLM interface is properly configured and ready to use */
  bool isConfigured() const;

  /** Reset the interface with a new config, clearing all conversation state.

   Keeps the existing JavaScript bridge and JSignal intact (avoids signal
   re-registration issues), but swaps the config, recreates the tool registry,
   and clears all conversation history and pending requests.

   @param config The new LLM configuration to use
   @throws std::logic_error if config is null or LLM API is not enabled
   */
  void resetWithConfig( const std::shared_ptr<const LlmConfig> &config );

  /** How a config change relates to the active provider, used to decide whether the existing
   conversation can be carried over to the new config (see applyConfigPreservingHistory). */
  enum class ConfigChange
  {
    SameModel,     ///< Same wire format and same active model name: fully compatible.
    ModelChanged,  ///< Same wire format, different active model: history is kept but model-specific
                   ///<   reasoning blocks must not be replayed (strip them).
    FormatChanged  ///< Different wire format: history cannot be reliably carried over.
  };//enum class ConfigChange

  /** Classify how `newCfg` differs from `oldCfg` for the active provider/model.  Compares the
   resolved wire format first (any difference -> FormatChanged), then the active model name.
   Defaults to FormatChanged on any error (the safe, "reset" choice). */
  static ConfigChange classifyConfigChange( const LlmConfig &oldCfg, const LlmConfig &newCfg );

  /** Apply a new config while KEEPING the existing conversation history.

   Swaps the config, tool registry, protocol, and debug logging (like resetWithConfig), but does
   not clear `m_history`/`m_currentConversation`.  Any in-flight request is finalized.  When
   `stripReasoning` is true, the model-specific reasoning state (thinking signatures,
   reasoning_content, reasoning_details) is cleared from every stored turn so it is not replayed to
   a model that did not produce it; the human-readable thinking text is kept for display.

   @param config The new LLM configuration to use (same wire format as the current one).
   @param stripReasoning Clear stored reasoning blocks from history (use when the model changed).
   @throws std::logic_error if config is null or LLM API is not enabled
   */
  void applyConfigPreservingHistory( const std::shared_ptr<const LlmConfig> &config, bool stripReasoning );
  
  /** Signal emitted when a new response is received from the LLM */
  Wt::Signal<>& conversationFinished();

  /** Signal emitted when an error response is received from the LLM */
  Wt::Signal<>& responseError();
  
  /** Check if a specific request ID is still pending */
  bool isRequestPending(int requestId) const;

  /** Get the conversation currently being processed (for tool execution context).

   Returns nullptr if no conversation is currently being processed.
   This is used by tools that need access to conversation-specific state like state machines.
   */
  std::shared_ptr<LlmInteraction> getCurrentConversation() const;

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

  /** JavaScript callback to handle LLM response.

   Also the entry point used by the native HTTP transport, which accumulates the whole response
   body and then hands it over here - so both transports converge on one response path.
   */
  void handleJavaScriptResponse(std::string response, int requestId);

  /** When true, tool calls will not be executed; instead an error result is
   returned to the LLM for each tool call.  Useful for follow-up conversations
   where we want the LLM to respond with text only.
   Defaults to false.
   */
  void setBlockToolCalls( bool block );

  /** Manually trigger context compaction (summarization of older exchanges).

   Unlike the automatic trigger in sendUserMessage(), this does not queue a user
   message to be sent after summarization completes.

   @return The visible summarization LlmInteraction (for creating a display widget),
           or nullptr if compaction could not be triggered (e.g., not enough history).
   */
  std::shared_ptr<LlmInteraction> triggerManualCompaction();

  /** Cancel all in-flight requests for this interface.

   Aborts the underlying browser fetch(es), finalizes every pending conversation with a "cancelled"
   error turn, clears the pending/deferred/summarization bookkeeping, and emits conversationFinished
   so the GUI re-enables input.  Safe to call when nothing is pending (no-op).
   */
  void cancelAll();

  /** Returns true if any request is currently in flight (pending API call or deferred tool/sub-agent
   results), i.e. cancelAll() would do something / a Stop button should be enabled.
   */
  bool hasActiveRequests() const;

  /** Direct all debug logging to the specified file path.
   Pass an empty string to disable file logging.
   Opens the file in append mode; logs a warning to stderr if the file cannot be opened.
   */
  void setDebugFile( const std::string &filePath );

private:

  
  InterSpec* m_interspec;
  bool m_block_tool_calls; // If true, tool calls return an error instead of being executed
  std::shared_ptr<const LlmConfig> m_config;
  std::shared_ptr<const LlmTools::ToolRegistry> m_tool_registry;
  std::shared_ptr<LlmConversationHistory> m_history;

  /** The wire-format translator for the active provider (OpenAI Chat / OpenAI Responses /
   Anthropic).  Built from the config's resolved apiFormat() in the constructor and resetWithConfig().
   */
  std::unique_ptr<const LlmApiProtocol> m_protocol;

  /** Instruction text pre-rendered for the active model by prepareModelInstructions() (at construction
   and on every config adoption), then read - never re-rendered - during per-turn request assembly.
   See LlmPromptTemplate.
   */
  std::map<AgentType, std::string> m_renderedSystemPrompt;                          // agent -> rendered base system prompt
  std::map<std::pair<AgentType,std::string>, std::string> m_renderedStateGuidance;  // (agent,state) -> rendered PromptGuidance
  std::map<std::pair<AgentType,std::string>, std::string> m_renderedStateEphemeral; // (agent,state) -> rendered ephemeral text
  std::string m_renderedCompaction;                                                 // rendered compaction prompt ("" -> use default)

  // Debug logging support
  std::ostream* m_debug_stream;              // Pointer to debug output stream (nullptr if no logging)
  std::unique_ptr<std::ofstream> m_debug_file; // File stream if logging to a file
  
  Wt::Signal<> m_conversationFinished; // Signal emitted when succesful final response from LLM is recieved.
  Wt::Signal<> m_responseError;    // Signal emitted when error responses are received

  /** Backstop against a hung conversation.  A single-shot timer (re)armed while work is outstanding
   (a request in flight, a response being rendered, or async tools running); if no terminal signal
   (conversationFinished/responseError) fires within sm_watchdog_timeout_ms, onWatchdogTimeout()
   forces the in-flight conversation(s) to fail and emits responseError so the GUI/benchmark cannot
   hang.  Every terminal emit goes through emitConversationFinished()/emitResponseError(), which
   disarm it.  Must exceed the JS-side fetch timeout (300 s) so normal stalls report through the
   existing error path first.
   */
  std::unique_ptr<Wt::WTimer> m_watchdogTimer;
  static constexpr int sm_watchdog_timeout_ms = 600000; // 10 minutes

  const std::string m_instanceId; // Unique ID for JavaScript bridge routing of responses to this instance
  std::unique_ptr<Wt::JSignal<std::string, int>> m_responseSignal; // For JavaScript bridge (response, requestId)

  // Request tracking
  int m_nextRequestId;

  struct PendingRequest
  {
    int requestId;
    std::weak_ptr<LlmInteraction> conversation;

    // Sub-agent support
    bool isSubAgentRequest = false;       // True if this request is from a sub-agent (we can probably get rid of this variable)

    // Summarization support
    bool isSummarizationRequest = false;  // True if this is a context summarization request

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
    nlohmann::json requestJson;
#endif
  };//struct PendingRequest
  
  std::map<int, PendingRequest> m_pendingRequests;

#if( USE_NATIVE_HTTP_CLIENT )
  /** Whether this request should go out through the native C++ HTTP client rather than the
   browser's fetch(), given the active provider's configured transport and whether a usable
   native backend was actually compiled in.  Logs once if native was asked for but is
   unavailable, then falls back to the browser. */
  bool useNativeHttpBackend() const;

  /** Issue `requestStr` to `endpoint` via the native transport.

   Accumulates the whole response body (the native path is deliberately non-streaming, matching
   the browser path), then delivers it to handleJavaScriptResponse() on the session thread.
   Mirrors the browser path's rate-limit and network retries, and its normalization of non-JSON
   HTTP failures, so that behaviour does not depend on which transport carried the request.

   @param attempt Zero for the first try; incremented by the internal retry logic.
   @param rateLimitRetry How many rate-limit retries have already been made.
   */
  void startNativeRequest( const std::string &endpoint,
                          const std::vector<std::pair<std::string,std::string>> &headers,
                          const std::string &requestStr,
                          int requestId,
                          int attempt,
                          int rateLimitRetry );

  /** Re-issue a native request after `delayMs`, for a rate-limit or network retry.  The timer is
   owned per request id so that a retry still pending at teardown cannot fire into a destroyed
   interface. */
  void scheduleNativeRetry( const std::string &endpoint,
                           const std::vector<std::pair<std::string,std::string>> &headers,
                           const std::string &requestStr,
                           int requestId,
                           int attempt,
                           int rateLimitRetry,
                           int delayMs );

  /** Abandon every in-flight native request and pending retry.  Silent by design: whoever asked
   for the cancel has already finalized the affected conversations. */
  void cancelNativeRequests();

  /** Live native requests, keyed by request id; erased when a request finishes or is cancelled.
   Destroying an entry cancels its transfer silently, which is what teardown wants. */
  std::map<int, std::unique_ptr<NativeHttp::Call>> m_nativeCalls;

  /** Timers backing the native path's retry delays.  Kept per request id so a retry that is
   pending when the interface is torn down does not fire into a destroyed object. */
  std::map<int, std::unique_ptr<Wt::WTimer>> m_nativeRetryTimers;
#endif //USE_NATIVE_HTTP_CLIENT

  // Track current conversation being processed (for tool execution context)
  std::weak_ptr<LlmInteraction> m_currentConversation;

  /** A tool round whose results cannot be sent back to the LLM yet, because one or more sub-agent
   invocations and/or async tool calls from that round are still outstanding.  The round finalizes - the
   entry is erased and the results are sent - when the last of them resolves, whichever that is.

   Async calls are tracked by invocation id rather than by a bare count so that resolution is
   *identity-based and exactly-once*: a callback that arrives twice, or arrives after the deadline sweep
   already reported the call as timed out, is recognized and dropped instead of double-resolving; and a
   callback that is skipped cannot silently pin the entry (the sweep reaps it by id).  An orphaned entry
   here blocks the terminal conversation signal, so this bookkeeping must be airtight.
   */
  struct DeferredToolResult {
    std::string conversationId;
    /** The conversation this round belongs to.  Held as well as the id because sub-agent conversations
     are NOT in LlmConversationHistory's top-level list (they hang off LlmToolCall::sub_agent_conversation),
     so findConversationByConversationId() cannot find them - resolving a sub-agent's round by id alone
     would silently never send its results.  The id remains the fallback if the weak_ptr has expired.
     */
    std::weak_ptr<LlmInteraction> conversation;
    std::vector<int> subAgentToolCallIds;  // The invoke_sub_agent tool call ID to update with summary

    /** An async tool call dispatched from this round and not yet resolved. */
    struct PendingAsyncCall {
      std::string toolName;
      int timeoutMs = 0;                               // the deadline that was applied, for messages
      std::chrono::steady_clock::time_point deadline;  // when onDeferredSweep() force-fails this call
    };
    std::map<std::string,PendingAsyncCall> pendingAsyncCalls;  // key is the tool-call invocation id
  };
  std::map<int, DeferredToolResult> m_deferredToolResults; // Key is the requestId the tool round came from

  // The deadline for an async tool call is not owned here: it comes from
  // LlmTools::effective_async_timeout_ms(), i.e. the tool's own LlmTools::SharedTool::asyncTimeoutMs
  // falling back to LlmTools::sm_default_async_timeout_ms, so that every consumer of asyncExecutor
  // (this class' deferred-round sweep, LlmMcpResource) applies the same deadline to a given tool.
  // That default is well under sm_watchdog_timeout_ms - static_assert'ed where it is applied - so a
  // stuck async tool is reported as a tool error, and the conversation continues, long before the
  // watchdog would fail the whole turn.

  /** How often m_deferredSweepTimer checks for overdue async calls. */
  static constexpr int sm_deferred_sweep_interval_ms = 5000;

  /** Repeating timer, running while any tool round has an outstanding async call, that force-resolves
   async tool calls which blew their deadline.  This - not the watchdog - is what makes a round with a
   stuck ASYNC TOOL finalize: the watchdog is an *inactivity* timer that is re-armed by every request,
   response and async step, so it never fires while the model keeps working around a stuck entry.

   NOTE: this covers `pendingAsyncCalls` only.  A round's `subAgentToolCallIds` leg has no deadline -
   a sub-agent turn is legitimately open-ended - so it is still resolved solely by the sub-agent's
   conversation_completion_handler (or, failing that, the watchdog).  A sub-agent whose own request
   ends in an API error therefore still pins its parent's round, since the error path pauses that
   conversation for the user to retry rather than completing it.
   */
  std::unique_ptr<Wt::WTimer> m_deferredSweepTimer;

  /** True if `callId` is still an outstanding async call of the tool round `parentRequestId`.

   Used to drop a *late* async result: once the deadline sweep has reported the call as an error, that
   error has already been sent to the model, so applying the late result would leave the stored history
   disagreeing with what was actually transmitted (and invalidate the provider's cached prefix).
   */
  bool isAsyncCallOutstanding( const int parentRequestId, const std::string &callId ) const;

  /** The single resolution point for an async tool call; idempotent.

   No-ops (with a debug log) when the round or the call id is not outstanding, so a duplicate or late
   callback cannot double-resolve.  Otherwise drops the call from the round's bookkeeping and, if that
   was the last outstanding sub-agent/async call, finalizes the round.
   */
  void resolveAsyncToolCall( const int parentRequestId, const std::string &callId );

  /** Sends the round's tool results back to the LLM and erases the entry, if nothing is outstanding.

   Shared by all three drivers (async callback, sub-agent completion, deadline sweep) so "who finished
   last" does not matter.  If the conversation no longer exists there is nothing to send, so the
   terminal signal is re-evaluated instead - the turn is never left stranded.
   */
  void finalizeDeferredIfReady( std::map<int, DeferredToolResult>::iterator pos );

  /** The conversation a deferred round belongs to, or nullptr if it no longer exists.

   Must be used instead of a bare findConversationByConversationId(): sub-agent conversations are not
   in the history's top-level list, so an id-only lookup silently fails for them.
   */
  std::shared_ptr<LlmInteraction> conversationForRound( const DeferredToolResult &round ) const;

  /** Emits conversationFinished if no requests and no deferred tool rounds remain.

   The success-path emit gate in handleApiResponse() only runs when an API response arrives; this lets
   the paths that reap a deferred entry outside of a response still terminate the turn.
   */
  void maybeEmitTerminalSignal();

  /** Starts/stops m_deferredSweepTimer to match whether any tool round is outstanding. */
  void updateDeferredSweepTimer();

  /** Force-fails async tool calls past their deadline, then re-evaluates the terminal signal. */
  void onDeferredSweep();

  /** Finds the Pending placeholder result for `callId` in `convo`, or nullptr.

   `owner` (optional) receives the LlmToolResults turn holding it, needed to emit asyncToolCompleted().
   */
  static LlmToolCall *findPendingToolResult( const std::shared_ptr<LlmInteraction> &convo,
                                             const std::string &callId,
                                             LlmToolResults **owner );

  /** Writes a short human-readable description of every outstanding deferred round (request id,
   conversation, and each pending async call's id/tool/age) to the debug stream.  Turns a recurrence of
   the orphaned-deferred-result hang into a one-line diagnosis.
   */
  void logOutstandingDeferred( const char * const context ) const;

#if( PERFORM_DEVELOPER_CHECKS )
  /** Asserts the invariants of m_deferredToolResults; call after every mutation of the map. */
  void assertDeferredInvariants() const;
#else
  void assertDeferredInvariants() const {}
#endif

  // Context summarization support.  A queue (not a single slot) so a second user message arriving
  // while summarization is in flight is not silently dropped.
  std::vector<std::weak_ptr<LlmInteraction>> m_summarizationPendingConvos;

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
                                         const std::string &thinkingContent,
                                         const std::string &thinkingSignature,
                                         const std::string &reasoningContent,
                                         const std::string &reasoningDetails,
                                         std::optional<size_t> promptTokens = std::nullopt,
                                         std::optional<size_t> completionTokens = std::nullopt );
  
  /** Parse text content for tool call requests (for models that don't support structured tool calls)

   Returns the number of tool calls processed.
   */
  std::pair<std::shared_ptr<LlmToolRequest>, std::shared_ptr<LlmToolResults>>
  parseContentForToolCallsAndSendResults( const std::string &content,
                                                const std::shared_ptr<LlmInteraction> &convo,
                                                const int requestId,
                                                const std::string &rawResponseContent,
                                                const std::string &thinkingContent,
                                                const std::string &thinkingSignature,
                                                const std::string &reasoningContent,
                                                const std::string &reasoningDetails );
  
  /** Strip <think>...</think> content from LLM responses */
  static std::string stripThinkingContent(const std::string& content);
  
  /** Extract thinking content and clean content from LLM responses */
  static std::pair<std::string, std::string> extractThinkingAndContent(const std::string& content);

  /** Get the system prompt for a specific agent from config
   */
  std::string getSystemPromptForAgent( const AgentType agentType ) const;

  /** Render all model-dependent instruction text (agent system prompts, per-state guidance/ephemeral
   text, and the compaction prompt) for the active model and cache it in the m_rendered* members.

   Called from the constructor and whenever a new config is adopted (resetWithConfig /
   applyConfigPreservingHistory) - i.e. exactly when the active model can change, NOT per request.
   render() is a no-op for marker-free text, so non-templated instructions are unchanged.
   */
  void prepareModelInstructions();

  /** Initialize state machine for a conversation if the agent has one defined.

   Creates a fresh copy of the agent's state machine and resets it to the initial state.

   @param convo The conversation to initialize state machine for
   */
  //void initializeStateMachineForConversation( std::shared_ptr<LlmInteraction> convo ) const;

  /** Build and send a compaction (summarization) request to the LLM.

   Creates a visible system conversation, builds a minimal request (no tools),
   sends it, and marks it as a summarization request.

   @param summarizationPrompt The prompt asking the LLM to summarize
   @param displayMessage The message shown to the user in the widget
   @return The visible summarization LlmInteraction
   */
  std::shared_ptr<LlmInteraction> sendCompactionRequest( const std::string &summarizationPrompt,
                                                          const std::string &displayMessage );

  /** Extract raw content verification from sendToolResultsToLLM into a helper. */
  void storeRawContentForToolResults( const std::shared_ptr<LlmInteraction> &convo,
                                      const nlohmann::json &followupRequest );

  /** Set up the JavaScript bridge for making HTTPS requests */
  void setupJavaScriptBridge();

  /** Finalize every still-pending conversation (set finishTime and record an error turn) and clear
   the pending/deferred bookkeeping.  Used on destruction and resetWithConfig() so in-flight
   conversations are not left perpetually "in progress" (which also stalls benchmark runs waiting on
   completion).  Does not abort the underlying browser fetch; late responses for cleared request IDs
   are ignored by handleJavaScriptResponse().

   @param reason Error message recorded on each finalized conversation.
   @param recordErrorTurn If true, append a (display-only) error turn to each conversation, which
          emits responseAdded to the GUI.  Pass false from the destructor to avoid emitting signals
          while this object is being torn down.
   */
  void failInFlightConversations( const std::string &reason, bool recordErrorTurn = true );

  /** Clear the model-specific reasoning state (thinking signatures, reasoning_content,
   reasoning_details) from every stored turn, so it is not replayed to a model that did not produce
   it.  The human-readable thinking text (thinkingContent) is kept for display. */
  void stripReasoningFromHistory();

  /** (Re)point m_debug_stream/m_debug_file from m_config->llmApi.debug_file. */
  void initDebugLoggingFromConfig();

  /** Send every conversation queued in m_summarizationPendingConvos (and clear the queue).  Called
   when a compaction request completes or fails so queued user messages are not stranded.
   */
  void flushSummarizationQueue();

  /** Watchdog helpers (see m_watchdogTimer).  armWatchdog() (re)starts the single-shot timer;
   disarmWatchdog() stops it; onWatchdogTimeout() force-fails outstanding conversations and emits
   responseError. */
  void armWatchdog();
  void disarmWatchdog();
  void onWatchdogTimeout();

  /** Sole emit points for the interface-level terminal signals; each disarms the watchdog first so
   there is exactly one disarm site per terminal transition. */
  void emitConversationFinished();
  void emitResponseError();
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
