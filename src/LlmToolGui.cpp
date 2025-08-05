#include "InterSpec_config.h"

#if( USE_LLM_INTERFACE )

#include "InterSpec/LlmToolGui.h"

#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>

#include <Wt/WText>
#include <Wt/WLineEdit>
#include <Wt/WTextArea>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WStringUtil>
#include <Wt/WSignal>

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace std;
using namespace Wt;

LlmToolGui::LlmToolGui(InterSpec *viewer, WContainerWidget *parent)
  : WContainerWidget(parent),
    m_viewer(viewer),
    m_llmInterface(nullptr),
    m_conversationDisplay(nullptr),
    m_inputEdit(nullptr),
    m_sendButton(nullptr),
    m_layout(nullptr),
    m_nextRequestId(1),
    m_isRequestPending(false),
    m_currentRequestId(-1)
{
  if (!m_viewer) {
    throw std::runtime_error("InterSpec instance cannot be null");
  }
  
  // Get the LLM interface from InterSpec
  m_llmInterface = m_viewer->llmInterface();
  if (!m_llmInterface) {
    throw std::runtime_error("LLM interface is not available");
  }
  
  wApp->useStyleSheet( "InterSpec_resources/LlmToolGui.css" );
  
  // Connect to LLM interface response signal
  m_llmInterface->responseReceived().connect(this, &LlmToolGui::handleResponseReceived);
  
  // Connect to spectrum change signal to cancel pending requests when foreground spectrum changes
  m_viewer->displayedSpectrumChanged().connect(this, &LlmToolGui::handleSpectrumChanged);
  
  initializeUI();
  
  // Initial display refresh
  refreshDisplay();
}

void LlmToolGui::initializeUI()
{
  addStyleClass("LlmToolGui");
  
  // Create the main layout
  m_layout = new WGridLayout();
  setLayout(m_layout);
  m_layout->setContentsMargins(5, 5, 5, 5);
  m_layout->setVerticalSpacing(5);
  m_layout->setHorizontalSpacing(5);
  
  // Create conversation display area (top part)
  m_conversationDisplay = new WTextArea();
  m_conversationDisplay->addStyleClass("llm-conversation-display");
  m_conversationDisplay->setReadOnly(true);
  m_conversationDisplay->setPlaceholderText("LLM conversation will appear here...");
  
  // Add conversation display to layout (row 0, spans 2 columns)
  m_layout->addWidget(m_conversationDisplay, 0, 0, 1, 2);
  m_layout->setRowStretch(0, 1); // Make conversation display expand
  
  // Create input line edit (bottom left)
  m_inputEdit = new WLineEdit();
  m_inputEdit->addStyleClass("llm-input-edit");
  m_inputEdit->setPlaceholderText("Type your question here...");
  m_inputEdit->setAutoComplete(false);
  
#if( BUILD_AS_OSX_APP || IOS )
  m_inputEdit->setAttributeValue("autocorrect", "off");
  m_inputEdit->setAttributeValue("spellcheck", "off");
#endif
  
  // Add input to layout (row 1, column 0)
  m_layout->addWidget(m_inputEdit, 1, 0);
  m_layout->setColumnStretch(0, 1); // Make input expand horizontally
  
  // Create send button (bottom right)
  m_sendButton = new WPushButton("Send");
  m_sendButton->addStyleClass("llm-send-button");
  
  // Add send button to layout (row 1, column 1)
  m_layout->addWidget(m_sendButton, 1, 1);
  

  
  // Connect signals
  m_inputEdit->enterPressed().connect(this, &LlmToolGui::handleInputSubmit);
  m_sendButton->clicked().connect(this, &LlmToolGui::handleSendButton);
  
  // Basic styling - can be enhanced with CSS later
  //m_conversationDisplay->setMinimumSize(WLength::Auto, WLength(200, WLength::Pixel));
  //m_inputEdit->setMinimumSize(WLength(200, WLength::Pixel), WLength::Auto);
  //m_sendButton->setMinimumSize(WLength(60, WLength::Pixel), WLength::Auto);
}


void LlmToolGui::focusInput()
{
  if (m_inputEdit) {
    m_inputEdit->setFocus(true);
  }
}

void LlmToolGui::clearHistory()
{
  if (m_conversationDisplay) {
    m_conversationDisplay->setText("");
  }
  
  // Also clear the LLM interface history if needed
  if (m_llmInterface) {
    auto history = m_llmInterface->getHistory();
    if (history) {
      history->clear();
    }
  }
}

void LlmToolGui::handleInputSubmit()
{
  handleSendButton(); // Same logic as clicking send button
}

void LlmToolGui::handleSendButton()
{
  if (!m_inputEdit) {
    return;
  }
  
  string message = m_inputEdit->text().toUTF8();
  if (message.empty()) {
    return;
  }
  
  // Clear the input field
  m_inputEdit->setText("");
  
  // Send the message
  sendMessage(message);
}

void LlmToolGui::sendMessage(const std::string& message)
{
  if (!m_llmInterface) {
    cout << "Error: LLM interface not available" << endl;
    return;
  }
  
  // Don't allow sending if a request is already pending
  if (m_isRequestPending) {
    cout << "Warning: Request already pending, ignoring new message" << endl;
    return;
  }
  
  try {
    // Disable input while request is pending
    setInputEnabled(false);
    m_isRequestPending = true;
    
    // Send message to LLM
    m_llmInterface->sendUserMessage(message);
    
    // Note: We can't directly track the request ID from sendUserMessage
    // We'll use a simple approach where we assume the request is pending
    // until we get a response or error
    
    cout << "Sent message to LLM" << endl;
    
    // Update display immediately to show user message
    updateConversationDisplay();
    
  } catch (const std::exception& e) {
    cout << "Error sending message to LLM: " << e.what() << endl;
    
    // Re-enable input on error
    setInputEnabled(true);
    m_isRequestPending = false;
    m_currentRequestId = -1;
    
    // Show error in conversation display
    string errorText = m_conversationDisplay->text().toUTF8();
    errorText += "\n[ERROR] Failed to send message: " + string(e.what()) + "\n";
    m_conversationDisplay->setText(WString::fromUTF8(errorText));
  }
}

void LlmToolGui::updateConversationDisplay()
{
  refreshDisplay();
}

void LlmToolGui::refreshDisplay()
{
  if( !m_conversationDisplay) {
    return;
  }
  
  if (!m_llmInterface ){
    m_conversationDisplay->setText("No conversation yet.");
    return;
  }
  
  // Check if we have a pending request
  if (m_isRequestPending) {
    // Check if any requests are still pending in the LLM interface
    // Since we can't track specific request IDs, we'll use a simpler approach
    // and re-enable input when we get a response (handled by responseReceived signal)
  }
  
  auto history = m_llmInterface->getHistory();
  if (!history || history->isEmpty()) {
    m_conversationDisplay->setText("No conversation history yet. Ask a question to get started!");
    return;
  }
  
  // Group messages by request for better display
  auto groupedConversations = groupConversationsByRequest();
  
  stringstream displayText;
  displayText << "=== LLM Conversation History ===\n\n";
  
  if (groupedConversations.empty()) {
    // Fallback: show all messages in chronological order
    const auto& conversations = history->getConversations();
    for (size_t i = 0; i < conversations.size(); ++i) {
      displayText << formatMessage(conversations[i]) << "\n";
    }
  } else {
    // Show conversations grouped by request ID
    for (const auto& [requestId, conversations] : groupedConversations) {
      if (requestId > 0) {
        displayText << "--- Request " << requestId << " ---\n";
      }
      
      for (const auto* conversation : conversations) {
        displayText << formatMessage(*conversation, requestId) << "\n";
      }
      displayText << "\n";
    }
  }
  
  m_conversationDisplay->setText(WString::fromUTF8(displayText.str()));
  
  // Scroll to bottom
  m_conversationDisplay->doJavaScript(
                                      m_conversationDisplay->jsRef() + ".scrollTop = " +
                                      m_conversationDisplay->jsRef() + ".scrollHeight;"
                                      );
}

void LlmToolGui::handleResponseReceived()
{
  // Re-enable input when we receive a response
  if (m_isRequestPending) {
    cout << "Response received, re-enabling input" << endl;
    setInputEnabled(true);
    m_isRequestPending = false;
    m_currentRequestId = -1;
  }
  
  // Refresh the display to show the new response
  refreshDisplay();
}

std::string LlmToolGui::formatMessage(const LlmConversationStart& conversation, int requestId)
{
  stringstream formatted;
  
  // Format timestamp
  auto time_t = std::chrono::system_clock::to_time_t(conversation.timestamp);
  formatted << "[" << std::put_time(std::localtime(&time_t), "%H:%M:%S") << "] ";
  
  // Format conversation start type and content
  switch (conversation.type) {
    case LlmConversationStart::Type::User:
      formatted << "YOU: " << conversation.content;
      break;
      
    case LlmConversationStart::Type::System:
      formatted << "SYSTEM: " << conversation.content;
      break;
      
    default:
      formatted << "UNKNOWN: " << conversation.content;
      break;
  }
  
  // // Add conversation ID if present
  //if (!conversation.conversationId.empty()) {
  //  formatted << " (Conv: " << conversation.conversationId << ")";
  //}
  
  // Add responses if any
  if (!conversation.responses.empty()) {
    formatted << "\n  Responses:";
    for (const auto& response : conversation.responses) {
      formatted << "\n    ";
      
      // Format response timestamp
      auto responseTime_t = std::chrono::system_clock::to_time_t(response.timestamp);
      formatted << "[" << std::put_time(std::localtime(&responseTime_t), "%H:%M:%S") << "] ";
      
      // Format response type and content
      switch (response.type) {
        case LlmConversationResponse::Type::Assistant:
          formatted << "LLM: " << response.content;
          break;
          
        case LlmConversationResponse::Type::ToolCall:
          formatted << "TOOL_CALL[" << response.toolName << "]: " << response.toolParameters.dump();
          break;
          
        case LlmConversationResponse::Type::ToolResult:
          formatted << "TOOL_RESULT: ";
          if( response.content.size() < 128 )
            formatted << response.content;
          else
            formatted << response.content.substr(0,125) + "...";
          break;
          
        case LlmConversationResponse::Type::Error:
          formatted << "ERROR: " << response.content;
          break;
          
        default:
          formatted << "UNKNOWN: " << response.content;
          break;
      }
    }
  }
  
  return formatted.str();
}

std::map<int, std::vector<const LlmConversationStart*>> LlmToolGui::groupConversationsByRequest()
{
  std::map<int, std::vector<const LlmConversationStart*>> groups;
  
  if (!m_llmInterface) {
    return groups;
  }
  
  auto history = m_llmInterface->getHistory();
  if (!history) {
    return groups;
  }
  
  const auto& conversations = history->getConversations();
  
  // Group conversations by conversation ID
  std::map<std::string, std::vector<const LlmConversationStart*>> conversationGroups;
  
  for (const auto& conversation : conversations) {
    std::string conversationId = conversation.conversationId.empty() ? "default" : conversation.conversationId;
    conversationGroups[conversationId].push_back(&conversation);
  }
  
  // Convert conversation groups to numbered groups for display
  int groupNumber = 1;
  for (const auto& [conversationId, conversationGroup] : conversationGroups) {
    groups[groupNumber] = conversationGroup;
    groupNumber++;
  }
  
  return groups;
}

void LlmToolGui::setInputEnabled(bool enabled)
{
  if (m_inputEdit) {
    m_inputEdit->setEnabled(enabled);
  }
  if (m_sendButton) {
    m_sendButton->setEnabled(enabled);
  }
  
  // Update visual appearance
  if (enabled) {
    m_inputEdit->removeStyleClass("disabled");
    m_sendButton->removeStyleClass("disabled");
  } else {
    m_inputEdit->addStyleClass("disabled");
    m_sendButton->addStyleClass("disabled");
  }
}

void LlmToolGui::handleSpectrumChanged()
{
  // If there's a pending request, cancel it and re-enable input
  if (m_isRequestPending) {
    cout << "Foreground spectrum changed during LLM request - cancelling request" << endl;
    cancelCurrentRequest();
  }
}

void LlmToolGui::cancelCurrentRequest()
{
  if (!m_isRequestPending || m_currentRequestId == -1) {
    return;
  }
  
  cout << "Cancelling LLM request ID: " << m_currentRequestId << endl;
  
  // Clear the pending request state
  m_isRequestPending = false;
  m_currentRequestId = -1;
  
  // Re-enable the input
  setInputEnabled(true);
  
  // Note: We can't actually cancel the HTTP request from C++, but the response will be ignored
  // when it comes back since we've cleared the pending request state
}

#endif // USE_LLM_INTERFACE
