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

#include <string>
#include <map>

#include <Wt/WSignal>
#include <Wt/WLineEdit>
#include <Wt/WContainerWidget>

//Forward declarations
class InterSpec;
class LlmInterface;
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

  /** Handle response received from LLM interface */
  void handleResponseReceived();
  
  /** Get the current conversation history for saving to SpecMeas */
  std::shared_ptr<std::vector<std::shared_ptr<LlmInteraction>>> getConversationHistory() const;

  /** Set the conversation history from SpecMeas */
  void setConversationHistory(const std::shared_ptr<std::vector<std::shared_ptr<LlmInteraction>>>& history);
  
  /** Clear the current conversation history */
  void clearConversationHistory();

protected:
  /** Handle user pressing Enter in the input field */
  void handleInputSubmit();

  /** Handle user clicking the send button */
  void handleSendButton();

  /** Process the user's input message */
  void sendMessage(const std::string& message);

  /** Handle retry request from error display */
  void handleRetry( std::shared_ptr<LlmInteraction> interaction );

  /** Export entire conversation as JSON */
  void exportConversationJson();

  /** Show confirmation dialog and clear conversation */
  void handleClearConversation();
  


private:
  InterSpec *m_viewer;              ///< The InterSpec instance
  std::unique_ptr<LlmInterface> m_llmInterface;  ///< The LLM interface for sending messages

  Wt::WContainerWidget *m_conversationContainer;  ///< Container holding LlmInteractionDisplay widgets
  Wt::WLineEdit *m_inputEdit;             ///< Input field for user messages
  Wt::WPushButton *m_sendButton;          ///< Button to send messages
  Wt::WPushButton *m_menuIcon;            ///< Menu icon for conversation-level actions
  Wt::WGridLayout *m_layout;              ///< Main layout manager

  bool m_isRequestPending;                 ///< Whether a request is currently pending
  
  /** Initialize the UI layout and widgets */
  void initializeUI();
  
  /** Enable or disable the input field and send button */
  void setInputEnabled(bool enabled);
  
  /** Handle spectrum changes - cancel pending requests if foreground spectrum changes */
  void handleSpectrumChanged();
  
  /** Cancel the current pending request */
  void cancelCurrentRequest();
};

#endif // USE_LLM_INTERFACE
#endif // LlmToolGui_h
