#ifndef LlmToolGui_h
#define LlmToolGui_h
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

#if( USE_LLM_INTERFACE )

#include <set>
#include <map>
#include <deque>
#include <string>
#include <vector>
#include <memory>
#include <variant>
#include <functional>

#include <Wt/WSignal>
#include <Wt/WLineEdit>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

//Forward declarations
class SpecMeas;
class InterSpec;
class LlmConfig;
class LlmInterface;
class LlmConfigWindow;
class LlmBenchmarkRunner;
class LlmJsSandboxBridge;
class LlmInteractionDisplay;
class LlmConversationHistory;
struct LlmInteraction;
class LlmInteractionTurn;

namespace Wt
{
  class WPushButton;
  class WGridLayout;
}

/** Widget for interacting with Large Language Models.
 * 
 * This widget provides a chat-like interface for users to interact with LLMs.
 * It displays conversation history in a scrollable text area and provides
 * an input field for new messages.
 */
class LlmToolGui : public Wt::WContainerWidget
{
public:
  /** Constructor.
   * @param viewer The InterSpec instance this tool belongs to
   * @param parent The parent widget (optional)
   */
  LlmToolGui(InterSpec *viewer, Wt::WContainerWidget *parent = nullptr);
  ~LlmToolGui();
  
  /** Returns if the user has added a `llm_config.xml` file, and enabled the `LlmApi`, and this was read in when the
   InterSpec server was started.
   */
  static bool llmToolIsConfigured();
  
  LlmInterface *llmInterface();
  
  /** Focus the input text field */
  void focusInput();

  /** Clear the conversation history display */
  void clearHistory();

  /** Handle successful response received from LLM interface */
  void handleConversationFinished();

  /** Handle error response received from LLM interface */
  void handleResponseError();
  
  /** Get the current conversation history for saving to SpecMeas */
  std::shared_ptr<std::vector<std::shared_ptr<LlmInteraction>>> getConversationHistory() const;

  /** Set the conversation history from SpecMeas */
  void setConversationHistory(const std::shared_ptr<std::vector<std::shared_ptr<LlmInteraction>>>& history);
  
  /** Clear the current conversation history */
  void clearConversationHistory();

  /** Submit a message as if the user typed it.
   Sets the text on the input field, then triggers send — identical to the user pressing Enter.
   Used by LlmBenchmarkRunner to drive automated evaluation.
   */
  void submitMessageAsUser( const std::string &message );

  /** Callback used by the MCP `assistant_submit_prompt` tool: invoked exactly once, on the session
   event loop, with either a result JSON (final assistant text + compact tool-call trace) or an
   error string. */
  using McpPromptCallback = std::function<void( std::variant<nlohmann::json, std::string> )>;

  /** Queue a prompt to be submitted into the conversation as if the user typed it, then invoke
   `callback` once the resulting turn (including all tool calls / sub-agents) finishes.
   Prompts are drained FIFO; only one runs at a time and never while a turn is already in progress.
   Resolves `callback` immediately with an error while a benchmark is running or the assistant is
   not configured.  Safe to call from the session event loop (e.g. an MCP async tool's Stage A,
   under the session UpdateLock).
   */
  void queuePromptForMcp( const std::string &prompt, McpPromptCallback callback );

  /** Disconnect the spectrum-changed handler (used by LlmBenchmarkRunner during SpectrumSequence). */
  void disconnectSpectrumChangedForBenchmark();

  /** Reconnect the spectrum-changed handler after a benchmark SpectrumSequence. */
  void reconnectSpectrumChangedForBenchmark();

  /** Stage an image to be sent with the next user message.
   Called from UploadedImgDisplay when user clicks "Send to AI assistant".
   Shows a thumbnail preview near the input area.
   */
  void stageImage( const std::string &base64Data, const std::string &mimeType,
                   const std::string &displayName, int widthPx, int heightPx );

  /** Returns true if the LLM model supports images and the GUI is available. */
  bool canAcceptImages() const;

  /** Returns true if a benchmark is currently running. */
  bool isBenchmarkRunning() const;

  /** Returns true if the LLM conversation should be preserved during a foreground spectrum change.
   Only true during benchmark SpectrumSequence steps where the conversation must survive mid-sequence
   spectrum loads.  For normal problem transitions, this returns false so the conversation is cleared.
   */
  bool shouldPreserveConversation() const;

  /** Returns the sandboxed-JavaScript bridge used by the `run_javascript` tool, or
   nullptr if it could not be constructed.  Owned by this widget. */
  LlmJsSandboxBridge *jsSandboxBridge();

protected:
  /** Handle user pressing Enter in the input field */
  void handleInputSubmit();

  /** Handle user clicking the send button */
  void handleSendButton();

  /** Process the user's input message */
  void sendMessage(const std::string& message);

  /** Handle retry request from error display */
  //void handleRetry( std::shared_ptr<LlmInteraction> interaction );

  /** Export entire conversation as JSON */
  void exportConversationJson();

  /** Show confirmation dialog and clear conversation */
  void handleClearConversation();

  /** Re-read LLM config from disk, clear GUI, and create a new LlmInterface */
  void handleResetLlmConfig();

  /** Manually trigger context compaction (summarization of older conversation history) */
  void handleCompactConversation();
  


private:
  InterSpec *m_viewer;              ///< The InterSpec instance
  std::shared_ptr<LlmInterface> m_llmInterface;  ///< The LLM interface for sending messages

  /** The single swappable content container (child of this widget's outer layout).
   Holds either the configured chat UI or the "not configured" panel; recreated when the
   configured/unconfigured state changes. */
  Wt::WContainerWidget *m_root;

  /** The open provider-settings window (nullptr if none), so it can be cleaned up. */
  LlmConfigWindow *m_configWindow;

  Wt::WContainerWidget *m_conversationContainer;  ///< Container holding LlmInteractionDisplay widgets
  Wt::WLineEdit *m_inputEdit;             ///< Input field for user messages
  Wt::WPushButton *m_sendButton;          ///< Button to send messages
  Wt::WPushButton *m_menuIcon;            ///< Menu icon for conversation-level actions
  Wt::WGridLayout *m_layout;              ///< Main layout manager

  bool m_isRequestPending;                 ///< Whether a request is currently pending

  /** Stored connection for spectrum change signal, so it can be disconnected/reconnected. */
  boost::signals2::connection m_spectrumChangedConnection;

  /** Benchmark runner (created on demand, nullable). */
  LlmBenchmarkRunner *m_benchmarkRunner;

  /** Hidden child container that bridges to the sandboxed JS iframe used
   by the `run_javascript` LLM tool.  May be null if construction failed. */
  LlmJsSandboxBridge *m_jsSandboxBridge;

  /** A staged image waiting to be sent with the next user message. */
  struct StagedImage
  {
    std::string base64Data;
    std::string mimeType;
    std::string displayName;
    int widthPx = 0;
    int heightPx = 0;
  };//struct StagedImage

  std::vector<StagedImage> m_stagedImages;

  /** Container for thumbnail previews of staged images, shown above the input field. */
  Wt::WContainerWidget *m_imagePreviewContainer;

  /** Remove a staged image by index and update the preview strip. */
  void removeStagedImage( size_t index );

  /** Remove all staged images and hide the preview strip. */
  void clearStagedImages();

  /** Initialize the configured chat UI (conversation, input, menu) into m_root. */
  void initializeUI();

  /** Build the full chat UI (creates the LlmInterface, hooks signals, calls initializeUI). */
  void buildConfiguredUi( const std::shared_ptr<const LlmConfig> &config );

  /** Build the "not configured" panel (message + a button to open the settings window). */
  void buildUnconfiguredUi();

  /** Tear down m_root and recreate an empty one in the outer layout; nulls child pointers. */
  void resetRoot();

  /** Open (or focus) the LLM provider settings window. */
  void openConfigWindow();

  /** Whether there is any (non-empty) conversation history right now.  Used by the settings window
   to decide whether to warn before applying a config change. */
  bool hasConversationHistory() const;

  /** Re-read config after the settings window saved, swap cache, and refresh this widget.
   Preserves the conversation when the change is compatible (same wire format); only does a full
   reset when the wire format changed (see LlmInterface::classifyConfigChange). */
  void handleConfigSaved();

  /** Enable or disable the input field and send button */
  void setInputEnabled(bool enabled);
  
  /** Handle spectrum changes - cancel pending requests if foreground spectrum changes.
   For passthrough/search-mode data, foreground sample changes are expected (the LLM is changing them)
   so we do not cancel the pending request.
   */
  void handleSpectrumChanged( SpecUtils::SpectrumType specType,
                              std::shared_ptr<SpecMeas> meas,
                              std::set<int> samples,
                              std::vector<std::string> detectors );
  
  /** Cancel the current pending request */
  void cancelCurrentRequest();

  /** Start running a benchmark from the given XML file path. */
  void handleStartBenchmark( const std::string &xmlPath );

  /** A prompt queued by the MCP `assistant_submit_prompt` tool, with its completion callback. */
  struct PendingMcpPrompt
  {
    std::string prompt;
    McpPromptCallback callback;
  };//struct PendingMcpPrompt

  /** FIFO of MCP-submitted prompts awaiting their turn (front = next to run). */
  std::deque<PendingMcpPrompt> m_mcpPromptQueue;

  /** Completion callback for the MCP prompt whose turn is currently in flight; empty when the
   active turn was not MCP-initiated, or when idle. */
  McpPromptCallback m_activeMcpCallback;

  /** If idle and no MCP prompt is in flight, submit the next queued prompt. */
  void drainMcpPromptQueue();

  /** Resolve the in-flight MCP prompt's callback (success → result JSON; else error string), then
   drain the next queued prompt. */
  void resolveActiveMcpPrompt( bool success );

  /** Build the JSON result (final assistant text + compact tool-call trace) for an MCP prompt. */
  nlohmann::json buildMcpPromptResult() const;
};

#endif // USE_LLM_INTERFACE
#endif // LlmToolGui_h
