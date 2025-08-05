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
#include <Wt/WTextArea>
#include <Wt/WLineEdit>
#include <Wt/WContainerWidget>

//Forward declarations
class InterSpec;
class LlmInterface;
class LlmConversationHistory;
struct LlmConversationStart;
struct LlmConversationResponse;

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
  
  /** Focus the input text field */
  void focusInput();
  
  /** Clear the conversation history display */
  void clearHistory();
  
  /** Refresh the conversation display from current history */
  void refreshDisplay();
  
  /** Handle response received from LLM interface */
  void handleResponseReceived();

protected:
  /** Handle user pressing Enter in the input field */
  void handleInputSubmit();
  
  /** Handle user clicking the send button */
  void handleSendButton();
  
  /** Process the user's input message */
  void sendMessage(const std::string& message);
  
  /** Update the display when new messages are received */
  void updateConversationDisplay();
  
  /** Format a single conversation for display */
  std::string formatMessage(const LlmConversationStart& conversation, int requestId = -1);
  
  /** Group conversations by request ID for display */
  std::map<int, std::vector<const LlmConversationStart*>> groupConversationsByRequest();
  


private:
  InterSpec *m_viewer;              ///< The InterSpec instance
  LlmInterface *m_llmInterface;     ///< The LLM interface for sending messages
  
  Wt::WTextArea *m_conversationDisplay;  ///< Text area showing conversation history
  Wt::WLineEdit *m_inputEdit;             ///< Input field for user messages
  Wt::WPushButton *m_sendButton;          ///< Button to send messages
  Wt::WGridLayout *m_layout;              ///< Main layout manager
  
  int m_nextRequestId;                    ///< Counter for tracking request IDs
  
  bool m_isRequestPending;                 ///< Whether a request is currently pending
  int m_currentRequestId;                  ///< ID of the current pending request
  
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
