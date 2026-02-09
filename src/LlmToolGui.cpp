#include "InterSpec_config.h"

#if( USE_LLM_INTERFACE )

#include "InterSpec/LlmToolGui.h"

#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>

#include <Wt/WText>
#include <Wt/WResource>
#include <Wt/WLineEdit>
#include <Wt/WScrollArea>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WStringUtil>
#include <Wt/WMemoryResource>
#include <Wt/WContainerWidget>

#include "SpecUtils/DateTime.h"

#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/LlmInteractionDisplay.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace std;
using namespace Wt;

LlmToolGui::LlmToolGui(InterSpec *viewer, WContainerWidget *parent)
  : WContainerWidget(parent),
    m_viewer(viewer),
    m_llmInterface(nullptr),
    m_conversationContainer(nullptr),
    m_inputEdit(nullptr),
    m_sendButton(nullptr),
    m_menuIcon(nullptr),
    m_layout(nullptr),
    m_isRequestPending(false)
{
  if( !m_viewer )
    throw std::runtime_error("InterSpec instance cannot be null");
  
  wApp->useStyleSheet( "InterSpec_resources/LlmToolGui.css" );
  
  std::shared_ptr<const LlmConfig> config;
  
  try
  {
    config = InterSpecServer::llm_config();
  }catch( std::exception &e )
  {
    cerr << "LLM config reading in failed: " << e.what() << endl;
  }
  
  if( !config || !config->llmApi.enabled || config->llmApi.apiEndpoint.empty() )
    throw std::runtime_error("LLM API not enabled");
  
  // Create our own LLM interface instance
  try
  {
    m_llmInterface = std::make_shared<LlmInterface>(m_viewer, config);
    
    // Connect to LLM interface response signals
    m_llmInterface->conversationFinished().connect(this, &LlmToolGui::handleConversationFinished);
    m_llmInterface->responseError().connect(this, &LlmToolGui::handleResponseError);

    // Connect to spectrum change signal to cancel pending requests when foreground spectrum changes
    m_viewer->displayedSpectrumChanged().connect(this, &LlmToolGui::handleSpectrumChanged);

    initializeUI();
  }catch( std::exception &e )
  {
    m_llmInterface.reset();
    WText *err_msg = new WText( "Error initializing LLM assistant: " + string(e.what()), this );
  }
}

LlmToolGui::~LlmToolGui()
{
  
}

bool LlmToolGui::llmToolIsConfigured()
{
  try
  {
    shared_ptr<const LlmConfig> config = InterSpecServer::llm_config();
    return (config && config->llmApi.enabled && !config->llmApi.apiEndpoint.empty());
  }catch( std::exception &e )
  {
    cerr << "LLM config invalid: " << e.what() << endl;
  }
  
  return false;
}


LlmInterface *LlmToolGui::llmInterface()
{
  return m_llmInterface.get();
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

  // Create wrapper container for scroll area and menu icon
  WContainerWidget *scrollWrapper = new WContainerWidget();
  scrollWrapper->setStyleClass( "LlmScrollWrapper" );

  // Create menu icon for conversation-level actions (positioned in upper-right via CSS)
  m_menuIcon = new WPushButton( scrollWrapper );
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn LlmMenuIcon" );
  m_menuIcon->clicked().preventPropagation();

  // Create popup menu
  PopupDivMenu *menu = new PopupDivMenu( m_menuIcon, PopupDivMenu::TransientMenu );

  // Add "Export Conversation as JSON" option
  PopupDivMenuItem *exportJsonItem = menu->addMenuItem( "Export Conversation as JSON" );
  exportJsonItem->triggered().connect( this, &LlmToolGui::exportConversationJson );

  // Add "Clear Conversation" option
  PopupDivMenuItem *clearItem = menu->addMenuItem( "Clear Conversation" );
  clearItem->triggered().connect( this, &LlmToolGui::handleClearConversation );

  // Add "Reset LLM Config" option - re-reads config files and creates a new LlmInterface
  PopupDivMenuItem *resetConfigItem = menu->addMenuItem( "Reset LLM Config" );
  resetConfigItem->triggered().connect( this, &LlmToolGui::handleResetLlmConfig );

  // Create scrollable container for conversation displays
  WScrollArea *scrollArea = new WScrollArea( scrollWrapper );
  scrollArea->addStyleClass("LlmConversationScrollArea");

  m_conversationContainer = new WContainerWidget();
  m_conversationContainer->addStyleClass("LlmConversationContainer");
  scrollArea->setWidget(m_conversationContainer);

  // Add scroll wrapper to layout (row 0, spans 2 columns)
  m_layout->addWidget(scrollWrapper, 0, 0, 1, 2);
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
}


void LlmToolGui::focusInput()
{
  if( m_inputEdit )
    m_inputEdit->setFocus(true);
}

void LlmToolGui::clearHistory()
{
  // Clear all interaction displays from the container
  if( m_conversationContainer )
    m_conversationContainer->clear();

  // Also clear the LLM interface history
  shared_ptr<LlmConversationHistory> history = m_llmInterface ? m_llmInterface->getHistory() : nullptr;
  if( history )
    history->clear();
}

void LlmToolGui::handleInputSubmit()
{
  handleSendButton(); // Same logic as clicking send button
}

void LlmToolGui::handleSendButton()
{
  if( !m_inputEdit )
    return;
  
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

    // Send message to LLM - this creates a new LlmInteraction in the history
    shared_ptr<LlmInteraction> convo = m_llmInterface->sendUserMessage(message);

    // Get the conversation history to find the new interaction
    if( convo )
    {
      // Create display widget for this interaction
      LlmInteractionDisplay *interactionDisplay =
      new LlmInteractionDisplay( convo, m_llmInterface, 0, m_conversationContainer );

      // Scroll to bottom
      const string js = "setTimeout(function() {"
      "  var scrollArea = " + m_conversationContainer->jsRef() + ".parentElement;"
      "  if (scrollArea) scrollArea.scrollTop = scrollArea.scrollHeight;"
      "}, 100);";
      doJavaScript( js );
    }//if( convo )

    cout << "Sent message to LLM" << endl;
  }catch( const std::exception &e )
  {
    cerr << "Error sending message to LLM: " << e.what() << endl;

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
    const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    const string now_str = SpecUtils::to_iso_string( now );
    const string debug_name = "llm_request_with_error_id_" + now_str + ".json";
#ifdef _WIN32
    const std::wstring wdebug_name = SpecUtils::convert_from_utf8_to_utf16(debug_name);
    std::ofstream output_request_json( wdebug_name.c_str(), ios::binary | ios::out );
#else
    std::ofstream output_request_json( debug_name.c_str(), ios::binary | ios::out);
#endif
    output_request_json << message;
    cerr << "Wrote request string to '" << debug_name << "'" << endl;
#endif
    
    // Re-enable input on error
    setInputEnabled(true);
    m_isRequestPending = false;

    // Show error in a dialog
    SimpleDialog *errorDialog = new SimpleDialog( "Error" );
    WText *errorText = new WText( "Failed to send message: " + string(e.what()), errorDialog->contents() );
    WPushButton *okBtn = errorDialog->addButton( "OK" );
    okBtn->clicked().connect( errorDialog, &SimpleDialog::accept );
  }
}

/*
void LlmToolGui::handleRetry( shared_ptr<LlmInteraction> interaction )
{
  if( !interaction || !m_llmInterface )
    return;

  cout << "\n=== Retry requested for conversation: " << interaction->conversationId << " ===" << endl;

  // Validate the conversation has an error as the last response
  if( interaction->responses.empty() )
  {
    cerr << "Cannot retry - conversation has no responses" << endl;
    return;
  }

  std::shared_ptr<LlmInteractionTurn> lastResponse = interaction->responses.back();
  if( lastResponse->type() != LlmInteractionTurn::Type::Error )
  {
    cerr << "Cannot retry - last response is not an error (type: "
         << static_cast<int>(lastResponse->type()) << ")" << endl;
    return;
  }

  std::shared_ptr<LlmInteractionError> error =
    std::dynamic_pointer_cast<LlmInteractionError>( lastResponse );

  if( !error )
  {
    cerr << "Cannot retry - failed to cast to LlmInteractionError" << endl;
    return;
  }

  cout << "Retrying error (type: " << static_cast<int>(error->errorType()) << ")" << endl;

  // Remove the error from conversation history
  // This allows the conversation to continue as if the error never happened
  interaction->responses.pop_back();

  // Determine what to retry by looking at second-to-last response
  // (since we just removed the error)

  if( interaction->responses.empty() )
  {
    // The error was the only response - retry the initial user message
    cout << "Retrying initial user message: " << interaction->content.substr(0, 50) << "..." << endl;

    // Disable input during retry
    setInputEnabled( false );
    m_isRequestPending = true;

    // Rebuild and send the request
    try
    {
      // Build API request with the conversation history (which now excludes the error)
      const nlohmann::json requestJson = m_llmInterface->buildMessagesArray( interaction );

      // Make tracked API call
      std::pair<int,std::string> request_id_content =
        m_llmInterface->makeTrackedApiCall( requestJson, interaction );

      cout << "Retry request sent with ID: " << request_id_content.first << endl;
    }
    catch( const std::exception &e )
    {
      cerr << "Error retrying request: " << e.what() << endl;

      // Re-add error so user can see what happened
      interaction->responses.push_back( error );

      // Re-enable input
      setInputEnabled( true );
      m_isRequestPending = false;
    }
  }
  else
  {
    // There were responses before the error - check what they were
    std::shared_ptr<LlmInteractionTurn> lastNonErrorResponse = interaction->responses.back();

    if( lastNonErrorResponse->type() == LlmInteractionTurn::Type::ToolResult )
    {
      // The error happened after tool results were sent
      // Retry by re-sending tool results to LLM
      cout << "Retrying after tool results" << endl;

      // Disable input during retry
      setInputEnabled( false );
      m_isRequestPending = true;

      try
      {
        const int requestId = m_llmInterface->sendToolResultsToLLM( interaction );
        cout << "Retry request sent with ID: " << requestId << endl;
      }
      catch( const std::exception &e )
      {
        cerr << "Error retrying tool results: " << e.what() << endl;

        // Re-add error
        interaction->responses.push_back( error );

        // Re-enable input
        setInputEnabled( true );
        m_isRequestPending = false;
      }
    }
    else if( lastNonErrorResponse->type() == LlmInteractionTurn::Type::AutoReply )
    {
      // The error happened after an auto-reply
      // Retry by re-sending the request with auto-reply
      cout << "Retrying after auto-reply" << endl;

      // Disable input during retry
      setInputEnabled( false );
      m_isRequestPending = true;

      try
      {
        // Build API request with conversation including auto-reply
        const nlohmann::json requestJson = m_llmInterface->buildMessagesArray( interaction );

        // Make tracked API call
        std::pair<int,std::string> request_id_content =
          m_llmInterface->makeTrackedApiCall( requestJson, interaction );

        cout << "Retry request sent with ID: " << request_id_content.first << endl;
      }
      catch( const std::exception &e )
      {
        cerr << "Error retrying request: " << e.what() << endl;

        // Re-add error
        interaction->responses.push_back( error );

        // Re-enable input
        setInputEnabled( true );
        m_isRequestPending = false;
      }
    }
    else
    {
      // Unexpected state - generic retry
      cerr << "Warning: unexpected conversation state for retry (last response type: "
           << static_cast<int>(lastNonErrorResponse->type()) << ")" << endl;

      // Disable input during retry
      setInputEnabled( false );
      m_isRequestPending = true;

      try
      {
        // Build API request with current conversation state
        const nlohmann::json requestJson = m_llmInterface->buildMessagesArray( interaction );

        // Make tracked API call
        std::pair<int,std::string> request_id_content =
          m_llmInterface->makeTrackedApiCall( requestJson, interaction );

        cout << "Retry request sent with ID: " << request_id_content.first << endl;
      }
      catch( const std::exception &e )
      {
        cerr << "Error retrying request: " << e.what() << endl;

        // Re-add error
        interaction->responses.push_back( error );

        // Re-enable input
        setInputEnabled( true );
        m_isRequestPending = false;
      }
    }
  }

  cout << "=== Retry handling complete ===" << endl;
}
*/


void LlmToolGui::handleConversationFinished()
{
  // Re-enable input when we receive a successful response
  // (errors are handled by handleResponseError() which keeps input disabled)
  m_isRequestPending = false;
  cout << "Conversation finished, re-enabling input" << endl;

  setInputEnabled( true );
  
  // The display is automatically updated via signals in LlmInteractionDisplay
  // No need to manually refresh
}

void LlmToolGui::handleResponseError()
{
  // Keep input disabled when we receive an error response
  // User must click "Retry" or "Continue Anyway" to proceed
  if( m_isRequestPending )
  {
    m_isRequestPending = false;
    cout << "Error response received, keeping input disabled" << endl;
    // Input stays disabled - user must use retry/continue buttons
  }

  // The display is automatically updated via signals in LlmInteractionDisplay
  // No need to manually refresh
}

void LlmToolGui::exportConversationJson()
{
  if( !m_llmInterface )
    return;

  shared_ptr<LlmConversationHistory> history = m_llmInterface->getHistory();
  if( !history || history->isEmpty() )
  {
    SimpleDialog *dialog = new SimpleDialog( "No Conversation" );
    WText *msg = new WText( "There is no conversation history to export.", dialog->contents() );
    WPushButton *okBtn = dialog->addButton( "OK" );
    okBtn->clicked().connect( dialog, &SimpleDialog::accept );
    return;
  }

  // Convert conversation history to JSON
  const nlohmann::json conversationsJson = history->toApiFormat();

  const string jsonStr = conversationsJson.dump(2);  // Pretty print with 2-space indent

  // Create memory resource for download
  WMemoryResource *resource = new WMemoryResource( "application/json" );
  resource->setData( reinterpret_cast<const unsigned char*>(jsonStr.data()), jsonStr.size() );
  resource->suggestFileName( "llm_conversation.json" );
  resource->setDispositionType( WResource::Attachment );

  // Trigger download
  WApplication::instance()->redirect( resource->url() );

  // Resource will auto-delete after download
}

void LlmToolGui::handleClearConversation()
{
  // Show confirmation dialog
  SimpleDialog *dialog = new SimpleDialog( "Clear Conversation" );
  WText *msg = new WText( "Are you sure you want to clear the entire conversation history? This cannot be undone.",
                         dialog->contents() );

  WPushButton *yesBtn = dialog->addButton( "Yes, Clear" );
  WPushButton *noBtn = dialog->addButton( "Cancel" );

  yesBtn->clicked().connect( std::bind( [this, dialog]() {
    clearHistory();
    dialog->accept();
  }));

  noBtn->clicked().connect( dialog, &SimpleDialog::accept );
}


void LlmToolGui::handleResetLlmConfig()
{
  // Show confirmation dialog before resetting
  SimpleDialog *dialog = new SimpleDialog( "Reset LLM Config" );
  new WText( "This will re-read the LLM configuration files and clear the current conversation. Continue?",
            dialog->contents() );

  WPushButton *yesBtn = dialog->addButton( "Yes, Reset" );
  WPushButton *noBtn = dialog->addButton( "Cancel" );

  yesBtn->clicked().connect( std::bind( [this, dialog]() {
    dialog->accept();

    // Try to load the config from disk; if it fails, show error and bail out
    shared_ptr<LlmConfig> config;

    try
    {
      config = LlmConfig::load();
    }catch( std::exception &e )
    {
      SimpleDialog *errDialog = new SimpleDialog( "Error Loading LLM Config" );
      new WText( e.what(), errDialog->contents() );
      errDialog->addButton( "Okay" );
      return;
    }

    if( !config || !config->llmApi.enabled || config->llmApi.apiEndpoint.empty() )
    {
      SimpleDialog *errDialog = new SimpleDialog( "Error Loading LLM Config" );
      new WText( "LLM API is not enabled in the configuration.", errDialog->contents() );
      errDialog->addButton( "Okay" );
      return;
    }

    // Update the global cached config
    InterSpecServer::set_llm_config( config );

    // Cancel any in-flight request and clear all conversation UI
    cancelCurrentRequest();
    clearHistory();

    // Reset the existing LlmInterface with the new config (preserves JS bridge and JSignal)
    try
    {
      m_llmInterface->resetWithConfig( config );
    }catch( std::exception &e )
    {
      SimpleDialog *errDialog = new SimpleDialog( "Error Resetting LLM" );
      new WText( e.what(), errDialog->contents() );
      errDialog->addButton( "Okay" );
      return;
    }

    setInputEnabled( true );
  }));

  noBtn->clicked().connect( dialog, &SimpleDialog::accept );
}


void LlmToolGui::setInputEnabled(bool enabled)
{
  if( m_inputEdit )
    m_inputEdit->setEnabled(enabled);
  if( m_sendButton )
    m_sendButton->setEnabled(enabled);
  
  // Update visual appearance
  if( m_inputEdit && m_sendButton )
  {
    if( enabled )
    {
      m_inputEdit->removeStyleClass("disabled");
      m_sendButton->removeStyleClass("disabled");
    }else
    {
      m_inputEdit->addStyleClass("disabled");
      m_sendButton->addStyleClass("disabled");
    }
  }//if( m_inputEdit && m_sendButton )
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
  if (!m_isRequestPending) {
    return;
  }

  cout << "Cancelling LLM request" << endl;

  // Clear the pending request state
  m_isRequestPending = false;

  // Re-enable the input
  setInputEnabled(true);

  // Note: We can't actually cancel the HTTP request from C++, but the response will be ignored
  // when it comes back since we've cleared the pending request state
}

std::shared_ptr<std::vector<std::shared_ptr<LlmInteraction>>> LlmToolGui::getConversationHistory() const
{
  if (!m_llmInterface) {
    return nullptr;
  }

  shared_ptr<LlmConversationHistory> history = m_llmInterface->getHistory();
  if (!history) {
    return nullptr;
  }

  const auto& conversations = history->getConversations();
  if (conversations.empty()) {
    return nullptr;
  }

  // Return a copy of the conversations vector for storage
  return std::make_shared<std::vector<std::shared_ptr<LlmInteraction>>>(conversations.begin(), conversations.end());
}

void LlmToolGui::setConversationHistory(const std::shared_ptr<std::vector<std::shared_ptr<LlmInteraction>>>& history)
{
  if (!m_llmInterface) {
    return;
  }

  shared_ptr<LlmConversationHistory> llmHistory = m_llmInterface ? m_llmInterface->getHistory() : nullptr;
  if( !llmHistory )
    return;

  // Clear current display
  if( m_conversationContainer )
    m_conversationContainer->clear();

  if( history && !history->empty() )
  {
    // Clear current history and add conversations from saved history
    llmHistory->clear();
    for( const std::shared_ptr<LlmInteraction> &conv : *history )
    {
      llmHistory->getConversations().push_back(conv);

      // Create display widget for this interaction
      LlmInteractionDisplay *interactionDisplay =
        new LlmInteractionDisplay( conv, m_llmInterface, 0, m_conversationContainer );
    }
  }else
  {
    // Clear the history if no valid history provided
    llmHistory->clear();
  }
}

void LlmToolGui::clearConversationHistory()
{
  if( !m_llmInterface )
    return;

  shared_ptr<LlmConversationHistory> history = m_llmInterface ? m_llmInterface->getHistory() : nullptr;
  if( history )
    history->clear();

  // Clear the display
  if( m_conversationContainer )
    m_conversationContainer->clear();
}

#endif // USE_LLM_INTERFACE
