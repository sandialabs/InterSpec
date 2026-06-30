#include "InterSpec_config.h"

#if( USE_LLM_INTERFACE )

#include "InterSpec/LlmToolGui.h"

#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WResource>
#include <Wt/WLineEdit>
#include <Wt/WScrollArea>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WStringUtil>
#include <Wt/Utils>
#include <Wt/WMemoryResource>
#include <Wt/WContainerWidget>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
#include "SpecUtils/StringAlgo.h" //For SpecUtils::convert_from_utf8_to_utf16
#endif

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmConfigWindow.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/LlmBenchmarkRunner.h"
#include "InterSpec/LlmJsSandboxBridge.h"
#include "InterSpec/LlmInteractionDisplay.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace std;
using namespace Wt;

LlmToolGui::LlmToolGui(InterSpec *viewer, WContainerWidget *parent)
  : WContainerWidget(parent),
    m_viewer(viewer),
    m_llmInterface(nullptr),
    m_root(nullptr),
    m_configWindow(nullptr),
    m_conversationContainer(nullptr),
    m_inputEdit(nullptr),
    m_sendButton(nullptr),
    m_menuIcon(nullptr),
    m_layout(nullptr),
    m_isRequestPending(false),
    m_imagePreviewContainer(nullptr),
    m_benchmarkRunner( nullptr ),
    m_jsSandboxBridge( nullptr )
{
  if( !m_viewer )
    throw std::runtime_error("InterSpec instance cannot be null");

  wApp->useStyleSheet( "InterSpec_resources/LlmToolGui.css" );
  // The "not configured" panel and the settings window share styles + strings from these bundles.
  wApp->useStyleSheet( "InterSpec_resources/LlmConfigWindow.css" );
  m_viewer->useMessageResourceBundle( "LlmConfigWindow" );

  // A single-cell outer layout holds m_root, which we swap between the configured chat UI and
  // the "not configured" panel.
  WGridLayout *outer = new WGridLayout();
  outer->setContentsMargins( 0, 0, 0, 0 );
  outer->setVerticalSpacing( 0 );
  outer->setHorizontalSpacing( 0 );
  setLayout( outer );
  m_root = new WContainerWidget();
  outer->addWidget( m_root, 0, 0 );
  outer->setRowStretch( 0, 1 );
  outer->setColumnStretch( 0, 1 );

  std::shared_ptr<const LlmConfig> config;
  try
  {
    config = InterSpecServer::llm_config();
  }catch( std::exception &e )
  {
    cerr << "LLM config reading in failed: " << e.what() << endl;
  }

  bool configured = false;
  try
  {
    configured = config && config->llmApi.enabled && !config->llmApi.apiEndpoint().empty();
  }catch( std::exception & )
  {
    configured = false; // e.g. a config with no providers
  }

  if( configured )
    buildConfiguredUi( config );
  else
    buildUnconfiguredUi();
}

LlmToolGui::~LlmToolGui()
{
  // If the settings window is still open, delete it now so its onSaved callback (which captures
  // this) can never fire against a destroyed widget.
  if( m_configWindow )
  {
    delete m_configWindow;
    m_configWindow = nullptr;
  }
}


void LlmToolGui::buildConfiguredUi( const std::shared_ptr<const LlmConfig> &config )
{
  try
  {
    m_llmInterface = std::make_shared<LlmInterface>(m_viewer, config);

    // Connect to LLM interface response signals
    m_llmInterface->conversationFinished().connect(this, &LlmToolGui::handleConversationFinished);
    m_llmInterface->responseError().connect(this, &LlmToolGui::handleResponseError);

    // Connect to spectrum change signal to cancel pending requests when foreground spectrum changes
    m_spectrumChangedConnection =
      m_viewer->displayedSpectrumChanged().connect(this, &LlmToolGui::handleSpectrumChanged);

    initializeUI();
  }catch( std::exception &e )
  {
    m_llmInterface.reset();
    new WText( "Error initializing LLM assistant: " + string(e.what()), m_root );
  }
}


void LlmToolGui::buildUnconfiguredUi()
{
  m_root->addStyleClass( "LlmToolGui" );

  WContainerWidget *panel = new WContainerWidget( m_root );
  panel->addStyleClass( "LlmNotConfigured" );

  WText *msg = new WText( WString::tr("lcw-not-configured"), panel );
  msg->addStyleClass( "LlmNotConfiguredMsg" );

  WPushButton *btn = new WPushButton( WString::tr("lcw-configure-btn"), panel );
  btn->clicked().connect( this, &LlmToolGui::openConfigWindow );
}


void LlmToolGui::resetRoot()
{
  WGridLayout *outer = dynamic_cast<WGridLayout *>( layout() );

  // Tear down the configured-UI state.  The spectrum-changed slot must be disconnected here (not
  // just reassigned in buildConfiguredUi) - boost::signals2 reassignment drops the handle without
  // disconnecting, which would leave an orphaned slot firing into a rebuilt widget.
  m_spectrumChangedConnection.disconnect();

  // Staged images live in this object (not under m_root); drop them so they don't leak into the
  // next conversation after a reconfigure.
  m_stagedImages.clear();

  if( m_root )
  {
    delete m_root;   // detaches from the layout and deletes all child widgets
    m_root = nullptr;
  }

  // Child widgets owned by the old root are now gone; clear dangling pointers.
  m_layout = nullptr;
  m_conversationContainer = nullptr;
  m_inputEdit = nullptr;
  m_sendButton = nullptr;
  m_menuIcon = nullptr;
  m_imagePreviewContainer = nullptr;
  m_jsSandboxBridge = nullptr;

  m_root = new WContainerWidget();
  if( outer )
  {
    outer->addWidget( m_root, 0, 0 );
    outer->setRowStretch( 0, 1 );
    outer->setColumnStretch( 0, 1 );
  }
}


void LlmToolGui::openConfigWindow()
{
  if( m_configWindow )
    return;  // already open

  m_configWindow = new LlmConfigWindow( m_viewer,
                                        std::bind( &LlmToolGui::handleConfigSaved, this ),
                                        std::bind( &LlmToolGui::hasConversationHistory, this ) );
  m_configWindow->finished().connect( std::bind( [this](){ m_configWindow = nullptr; } ) );
}


bool LlmToolGui::hasConversationHistory() const
{
  if( !m_llmInterface )
    return false;
  const std::shared_ptr<LlmConversationHistory> history = m_llmInterface->getHistory();
  return history && !history->isEmpty();
}


void LlmToolGui::handleConfigSaved()
{
  // Re-read the full config from disk (also loads agents/tools), swap the server-wide cache, then
  // refresh this widget to reflect the new state.
  std::shared_ptr<LlmConfig> config;
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

  InterSpecServer::set_llm_config( config );

  bool configured = false;
  try
  {
    configured = config && config->llmApi.enabled && !config->llmApi.apiEndpoint().empty();
  }catch( std::exception & )
  {
    configured = false;
  }

  if( configured )
  {
    if( m_llmInterface )
    {
      // Decide whether the conversation can be carried over to the new config.  Classify against the
      // interface's current (old) config before applying.
      const std::shared_ptr<const LlmConfig> oldConfig = m_llmInterface->config();
      LlmInterface::ConfigChange change = LlmInterface::ConfigChange::FormatChanged;
      if( oldConfig )
        change = LlmInterface::classifyConfigChange( *oldConfig, *config );

      cancelCurrentRequest();
      try
      {
        if( change == LlmInterface::ConfigChange::FormatChanged )
        {
          // The wire format changed - the conversation can't be reliably carried over, so reset it.
          clearHistory();
          m_llmInterface->resetWithConfig( config );
        }else
        {
          // Same wire format: keep the conversation.  If the active model changed, strip the
          // model-specific reasoning blocks so they aren't replayed to a different model.
          const bool stripReasoning = (change == LlmInterface::ConfigChange::ModelChanged);
          m_llmInterface->applyConfigPreservingHistory( config, stripReasoning );
        }
      }catch( std::exception &e )
      {
        SimpleDialog *errDialog = new SimpleDialog( "Error Updating LLM" );
        new WText( e.what(), errDialog->contents() );
        errDialog->addButton( "Okay" );
        return;
      }
      setInputEnabled( true );
    }else
    {
      // Was unconfigured: build the full chat UI now.
      resetRoot();
      buildConfiguredUi( config );
    }
  }else
  {
    // User disabled the LLM (or left it unconfigured): show the not-configured panel.
    m_llmInterface.reset();
    resetRoot();
    buildUnconfiguredUi();
  }
}

bool LlmToolGui::llmToolIsConfigured()
{
  try
  {
    shared_ptr<const LlmConfig> config = InterSpecServer::llm_config();
    return (config && config->llmApi.enabled && !config->llmApi.apiEndpoint().empty());
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
  m_root->addStyleClass("LlmToolGui");

  // Create the main layout
  m_layout = new WGridLayout();
  m_root->setLayout(m_layout);
  m_layout->setContentsMargins(5, 5, 5, 5);
  m_layout->setVerticalSpacing(5);
  m_layout->setHorizontalSpacing(5);

  // Create wrapper container for scroll area and menu icon
  WContainerWidget *scrollWrapper = new WContainerWidget();
  scrollWrapper->setStyleClass( "LlmScrollWrapper" );


  // Sandboxed JS bridge for the `run_javascript` LLM tool.  Hidden child
  // of this widget; its iframe is positioned off-screen by the JS helper.
  try
  {
    m_jsSandboxBridge = new LlmJsSandboxBridge( scrollWrapper );
  }catch( std::exception &e )
  {
    cerr << "Failed to construct LlmJsSandboxBridge: " << e.what() << endl;
    m_jsSandboxBridge = nullptr;
  }

  // Create menu icon for conversation-level actions (positioned in upper-right via CSS)
  m_menuIcon = new WPushButton( scrollWrapper );
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn LlmMenuIcon" );
  // Note: deliberately do NOT preventPropagation() here (doing this lets menu open, but not close when you click the button).

  // Create popup menu
  PopupDivMenu *menu = new PopupDivMenu( m_menuIcon, PopupDivMenu::TransientMenu );

  // Add "Export Conversation as JSON" option
  PopupDivMenuItem *exportJsonItem = menu->addMenuItem( "Export Conversation as JSON" );
  exportJsonItem->triggered().connect( this, &LlmToolGui::exportConversationJson );

  // Add "Clear Conversation" option
  PopupDivMenuItem *clearItem = menu->addMenuItem( "Clear Conversation" );
  clearItem->triggered().connect( this, &LlmToolGui::handleClearConversation );

  // Add "Compact Conversation" option - manually trigger context summarization
  PopupDivMenuItem *compactItem = menu->addMenuItem( "Compact Conversation" );
  compactItem->triggered().connect( this, &LlmToolGui::handleCompactConversation );

  // Add "Configure LLM Provider…" option - opens the graphical provider/model settings window
  PopupDivMenuItem *configItem = menu->addMenuItem( WString::tr("lcw-configure-menu") );
  configItem->triggered().connect( this, &LlmToolGui::openConfigWindow );

  // Add "Reset LLM Config" option - re-reads config files and creates a new LlmInterface
  PopupDivMenuItem *resetConfigItem = menu->addMenuItem( "Reset LLM Config" );
  resetConfigItem->triggered().connect( this, &LlmToolGui::handleResetLlmConfig );

  // Add benchmark submenu - scan for *_llm_benchmark.xml files in test_datasets/
  {
    const string dataDir = InterSpec::staticDataDirectory();
    // Go up from data/ to the repo root, then into test_datasets/
    const string repoRoot = SpecUtils::parent_path( dataDir );
    const string testDir = SpecUtils::append_path( repoRoot, "test_datasets" );

    const vector<pair<string,string>> benchmarkFiles = LlmBenchmarkRunner::findBenchmarkFiles( testDir );

    if( !benchmarkFiles.empty() )
    {
      menu->addSeparator();

      PopupDivMenu *benchmarkSubMenu = menu->addPopupMenuItem( "Benchmarks / Quizzes" );

      for( const pair<string,string> &benchFile : benchmarkFiles )
      {
        const string displayName = benchFile.first;
        const string filePath = benchFile.second;

        PopupDivMenuItem *item = benchmarkSubMenu->addMenuItem( displayName );
        item->triggered().connect(
          std::bind( [this]( const string &path ){
            handleStartBenchmark( path );
          }, filePath )
        );
      }
    }
  }

  // Create scrollable container for conversation displays
  WScrollArea *scrollArea = new WScrollArea( scrollWrapper );
  scrollArea->addStyleClass("LlmConversationScrollArea");

  m_conversationContainer = new WContainerWidget();
  m_conversationContainer->addStyleClass("LlmConversationContainer");
  scrollArea->setWidget(m_conversationContainer);

  // Add scroll wrapper to layout (row 0, spans 2 columns)
  m_layout->addWidget(scrollWrapper, 0, 0, 1, 2);
  m_layout->setRowStretch(0, 1); // Make conversation display expand

  // Create a vertical container for the image preview strip + input edit
  WContainerWidget *inputColumn = new WContainerWidget();
  inputColumn->addStyleClass( "LlmInputColumn" );

  // Image preview strip (hidden by default)
  m_imagePreviewContainer = new WContainerWidget( inputColumn );
  m_imagePreviewContainer->addStyleClass( "LlmStagedImagePreview" );
  m_imagePreviewContainer->hide();

  // Create input line edit
  m_inputEdit = new WLineEdit( inputColumn );
  m_inputEdit->addStyleClass("llm-input-edit");
  m_inputEdit->setPlaceholderText( WString( "Type your question here..." ) );
  m_inputEdit->setAutoComplete(false);

#if( BUILD_AS_OSX_APP || IOS )
  m_inputEdit->setAttributeValue("autocorrect", "off");
  m_inputEdit->setAttributeValue("spellcheck", "off");
#endif

  // Add input column to layout (row 1, column 0)
  m_layout->addWidget(inputColumn, 1, 0);
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


bool LlmToolGui::canAcceptImages() const
{
  try
  {
    const shared_ptr<const LlmConfig> config = InterSpecServer::llm_config();
    return config && config->llmApi.supportsImages();
  }catch( ... )
  {
    return false;
  }
}//bool canAcceptImages() const


void LlmToolGui::stageImage( const string &base64Data, const string &mimeType,
                              const string &displayName, int widthPx, int heightPx )
{
  StagedImage staged;
  staged.base64Data = base64Data;
  staged.mimeType = mimeType;
  staged.displayName = displayName;
  staged.widthPx = widthPx;
  staged.heightPx = heightPx;
  m_stagedImages.push_back( std::move( staged ) );

  if( !m_imagePreviewContainer )
    return;

  // Create a thumbnail entry for this image
  const size_t index = m_stagedImages.size() - 1;

  WContainerWidget *entry = new WContainerWidget( m_imagePreviewContainer );
  entry->addStyleClass( "LlmStagedImageEntry" );

  // Create thumbnail using WMemoryResource (same pattern as UploadedImgDisplay)
  const string decoded = Wt::Utils::base64Decode( base64Data );
  vector<unsigned char> bytes( decoded.begin(), decoded.end() );
  WMemoryResource *thumbResource = new WMemoryResource( mimeType, entry );
  thumbResource->setData( bytes );

  WImage *thumb = new WImage( WLink( thumbResource ), entry );
  thumb->addStyleClass( "LlmStagedImageThumb" );

  // Filename label
  WText *nameLabel = new WText( WString::fromUTF8( displayName ), entry );
  nameLabel->addStyleClass( "LlmStagedImageName" );

  // Remove button
  WPushButton *removeBtn = new WPushButton( entry );
  removeBtn->setStyleClass( "LlmStagedImageRemove" );
  removeBtn->setText( WString::fromUTF8( "\xc3\x97" ) ); // multiplication sign as X
  removeBtn->clicked().connect( std::bind( [this]( size_t idx ){ removeStagedImage( idx ); }, index ) );

  m_imagePreviewContainer->show();
}//void stageImage(...)


void LlmToolGui::removeStagedImage( const size_t index )
{
  if( index >= m_stagedImages.size() )
    return;

  m_stagedImages.erase( m_stagedImages.begin() + index );

  // Rebuild the preview strip (simpler than managing individual widget indices)
  if( m_imagePreviewContainer )
  {
    m_imagePreviewContainer->clear();

    if( m_stagedImages.empty() )
    {
      m_imagePreviewContainer->hide();
      return;
    }

    // Re-add all remaining thumbnails
    for( size_t i = 0; i < m_stagedImages.size(); ++i )
    {
      const StagedImage &img = m_stagedImages[i];

      WContainerWidget *entry = new WContainerWidget( m_imagePreviewContainer );
      entry->addStyleClass( "LlmStagedImageEntry" );

      const string decoded = Wt::Utils::base64Decode( img.base64Data );
      vector<unsigned char> bytes( decoded.begin(), decoded.end() );
      WMemoryResource *thumbResource = new WMemoryResource( img.mimeType, entry );
      thumbResource->setData( bytes );

      WImage *thumb = new WImage( WLink( thumbResource ), entry );
      thumb->addStyleClass( "LlmStagedImageThumb" );

      WText *nameLabel = new WText( WString::fromUTF8( img.displayName ), entry );
      nameLabel->addStyleClass( "LlmStagedImageName" );

      WPushButton *removeBtn = new WPushButton( entry );
      removeBtn->setStyleClass( "LlmStagedImageRemove" );
      removeBtn->setText( WString::fromUTF8( "\xc3\x97" ) );
      removeBtn->clicked().connect( std::bind( [this]( size_t idx ){ removeStagedImage( idx ); }, i ) );
    }
  }
}//void removeStagedImage(...)


void LlmToolGui::clearStagedImages()
{
  m_stagedImages.clear();
  if( m_imagePreviewContainer )
  {
    m_imagePreviewContainer->clear();
    m_imagePreviewContainer->hide();
  }
}//void clearStagedImages()


void LlmToolGui::handleInputSubmit()
{
  handleSendButton(); // Same logic as clicking send button
}

void LlmToolGui::handleSendButton()
{
  // While a request is in flight the Send button acts as a Stop button.
  if( m_isRequestPending )
  {
    if( m_llmInterface )
      m_llmInterface->cancelAll();
    return;
  }

  if( !m_inputEdit )
    return;

  const string message = m_inputEdit->text().toUTF8();
  if( message.empty() && m_stagedImages.empty() )
    return;

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

    // Collect staged images
    vector<LlmToolCall::ImageContent> images;
    for( StagedImage &staged : m_stagedImages )
    {
      LlmToolCall::ImageContent img;
      img.base64Data = std::move( staged.base64Data );
      img.mimeType = std::move( staged.mimeType );
      img.widthPx = staged.widthPx;
      img.heightPx = staged.heightPx;
      images.push_back( std::move( img ) );
    }
    clearStagedImages();

    // Send message to LLM - this creates a new LlmInteraction in the history
    shared_ptr<LlmInteraction> convo = m_llmInterface->sendUserMessage( message, std::move( images ) );

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

    //cout << "Sent message to LLM" << endl;
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
  //cout << "Conversation finished, re-enabling input" << endl;

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

// Set to 1 to export the full API request (system prompt, tools, model params, all messages)
// matching exactly what is sent to the LLM. Set to 0 for the legacy bare messages array.
#define LLM_EXPORT_FULL_API_JSON 1

#if( LLM_EXPORT_FULL_API_JSON )
  // Build the exact same JSON payload that is sent to the LLM API - this includes
  // the system prompt, tool definitions, tool_choice, model parameters, all messages,
  // and the ephemeral state machine reminder (if applicable).
  // Prefer the currently-active conversation; fall back to the most recent from history.
  shared_ptr<LlmInteraction> currentConvo = m_llmInterface->getCurrentConversation();
  if( !currentConvo )
  {
    const vector<shared_ptr<LlmInteraction>> &convos = history->getConversations();
    if( !convos.empty() )
      currentConvo = convos.back();
  }

  if( !currentConvo )
  {
    SimpleDialog *dialog = new SimpleDialog( "No Active Conversation" );
    WText *msg = new WText( "There is no active conversation to export.", dialog->contents() );
    WPushButton *okBtn = dialog->addButton( "OK" );
    okBtn->clicked().connect( dialog, &SimpleDialog::accept );
    return;
  }

  nlohmann::json requestJson;
  try
  {
    requestJson = m_llmInterface->buildMessagesArray( currentConvo );
  }
  catch( std::exception &e )
  {
    SimpleDialog *dialog = new SimpleDialog( "Export Failed" );
    WText *msg = new WText( string("Failed to build export JSON: ") + e.what(), dialog->contents() );
    WPushButton *okBtn = dialog->addButton( "OK" );
    okBtn->clicked().connect( dialog, &SimpleDialog::accept );
    return;
  }

  const string jsonStr = requestJson.dump(2);  // Pretty print with 2-space indent
#else
  // Legacy format: bare messages array only (no system prompt, tools, or model params).
  const nlohmann::json conversationsJson = history->toApiFormat();
  const string jsonStr = conversationsJson.dump(2);  // Pretty print with 2-space indent
#endif  //LLM_EXPORT_FULL_API_JSON

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

    if( !config || !config->llmApi.enabled || config->llmApi.apiEndpoint().empty() )
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


void LlmToolGui::handleCompactConversation()
{
  if( !m_llmInterface )
    return;

  if( m_isRequestPending )
  {
    SimpleDialog *dialog = new SimpleDialog( "Cannot Compact" );
    new WText( "Please wait for the current request to finish before compacting.", dialog->contents() );
    dialog->addButton( "OK" );
    return;
  }

  try
  {
    shared_ptr<LlmInteraction> summarizationConvo = m_llmInterface->triggerManualCompaction();

    if( !summarizationConvo )
    {
      SimpleDialog *dialog = new SimpleDialog( "Cannot Compact" );
      new WText( "Not enough conversation history to compact.", dialog->contents() );
      dialog->addButton( "OK" );
      return;
    }

    // Create display widget for the summarization conversation
    new LlmInteractionDisplay( summarizationConvo, m_llmInterface, 0, m_conversationContainer );

    // Scroll to bottom
    const string js = "setTimeout(function() {"
    "  var scrollArea = " + m_conversationContainer->jsRef() + ".parentElement;"
    "  if (scrollArea) scrollArea.scrollTop = scrollArea.scrollHeight;"
    "}, 100);";
    doJavaScript( js );

    // Disable input while compaction is in-flight
    setInputEnabled( false );
    m_isRequestPending = true;
  }catch( const std::exception &e )
  {
    cerr << "Error triggering conversation compaction: " << e.what() << endl;

    SimpleDialog *dialog = new SimpleDialog( "Compaction Error" );
    new WText( string("Failed to compact conversation: ") + e.what(), dialog->contents() );
    dialog->addButton( "OK" );
  }
}


void LlmToolGui::setInputEnabled(bool enabled)
{
  if( m_inputEdit )
  {
    m_inputEdit->setEnabled( enabled );
    if( enabled )
      m_inputEdit->setPlaceholderText( WString( "Type your question here..." ) );
    else
      m_inputEdit->setPlaceholderText( WString( "Working..." ) );
  }

  // While working, keep the button enabled but relabel it "Stop" so the user can cancel the request
  // (handleSendButton() routes to cancelAll() when a request is pending).
  if( m_sendButton )
  {
    m_sendButton->setEnabled( true );
    m_sendButton->setText( enabled ? "Send" : "Stop" );
    if( enabled )
      m_sendButton->removeStyleClass( "llm-stop-button" );
    else
      m_sendButton->addStyleClass( "llm-stop-button" );
  }

  // Update visual appearance.  Only the input is greyed out while working; the button stays active
  // (as "Stop"), so it must not get the "disabled" class.
  if( m_inputEdit )
  {
    if( enabled )
      m_inputEdit->removeStyleClass("disabled");
    else
      m_inputEdit->addStyleClass("disabled");
  }
  if( m_sendButton )
    m_sendButton->removeStyleClass("disabled");
}

void LlmToolGui::handleSpectrumChanged( SpecUtils::SpectrumType specType,
                                        std::shared_ptr<SpecMeas> meas,
                                        std::set<int> /*samples*/,
                                        std::vector<std::string> /*detectors*/ )
{
  // For passthrough/search-mode data, foreground sample changes are expected
  //  (the LLM SpectrumTimeHistory agent changes them), so do not cancel the pending request.
  if( meas && meas->passthrough() )
    return;

  // If there's a pending request, cancel it and re-enable input
  if( m_isRequestPending )
  {
    cout << "Foreground spectrum changed during LLM request - cancelling request" << endl;
    cancelCurrentRequest();
  }

  setInputEnabled(true);
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
      if( !conv )
        continue;

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


void LlmToolGui::submitMessageAsUser( const string &message )
{
  if( !m_inputEdit )
    return;

  // During benchmarks, force-clear any pending/disabled state from previous errors
  if( isBenchmarkRunning() )
  {
    m_isRequestPending = false;
    setInputEnabled( true );
  }

  m_inputEdit->setText( WString::fromUTF8( message ) );
  handleSendButton();
}


void LlmToolGui::disconnectSpectrumChangedForBenchmark()
{
  m_spectrumChangedConnection.disconnect();
}


void LlmToolGui::reconnectSpectrumChangedForBenchmark()
{
  m_spectrumChangedConnection =
    m_viewer->displayedSpectrumChanged().connect( this, &LlmToolGui::handleSpectrumChanged );
}


bool LlmToolGui::isBenchmarkRunning() const
{
  return m_benchmarkRunner && m_benchmarkRunner->isRunning();
}


bool LlmToolGui::shouldPreserveConversation() const
{
  // During benchmark SpectrumSequence steps, the spectrum-changed handler is disconnected
  //  so that conversation history survives mid-sequence spectrum loads.
  //  When disconnected, the connection object is "empty" / not connected.
  return isBenchmarkRunning() && !m_spectrumChangedConnection.connected();
}


void LlmToolGui::handleStartBenchmark( const string &xmlPath )
{
  if( m_benchmarkRunner && m_benchmarkRunner->isRunning() )
  {
    SimpleDialog *dialog = new SimpleDialog( "Benchmark Running" );
    new WText( "A benchmark is already running. Please wait for it to finish.", dialog->contents() );
    dialog->addButton( "OK" );
    return;
  }

  // Clean up any previous runner
  if( m_benchmarkRunner )
  {
    delete m_benchmarkRunner;
    m_benchmarkRunner = nullptr;
  }

  cout << "[LlmBenchmark] Starting benchmark: " << xmlPath << endl;

  m_benchmarkRunner = new LlmBenchmarkRunner( m_viewer, this );
  m_benchmarkRunner->startBenchmark( xmlPath );
}


LlmJsSandboxBridge *LlmToolGui::jsSandboxBridge()
{
  return m_jsSandboxBridge;
}


#endif // USE_LLM_INTERFACE
