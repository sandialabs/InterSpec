#include "InterSpec_config.h"

#if( USE_LLM_INTERFACE )

#include "InterSpec/LlmToolGui.h"

#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <algorithm>

#include <boost/bind/bind.hpp>

#include <Wt/WText.h>
#include <Wt/WImage.h>
#include <Wt/WAnchor.h>
#include <Wt/WRandom.h>
#include <Wt/WTextArea.h>
#include <Wt/WResource.h>
#include <Wt/WPushButton.h>
#include <Wt/WVBoxLayout.h>
#include <Wt/WApplication.h>
#include <Wt/WStringUtil.h>
#include <Wt/Utils.h>
#include <Wt/WMemoryResource.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WJavaScriptPreamble.h>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
#include "SpecUtils/StringAlgo.h" //For SpecUtils::convert_from_utf8_to_utf16
#endif

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmConfigWindow.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/LlmBenchmarkRunner.h"
#include "InterSpec/LlmJsSandboxBridge.h"
#include "InterSpec/LlmInteractionDisplay.h"
#include "InterSpec/LlmConversationHistory.h"

#if( BUILD_AS_MCP_HTTP_SERVER )
#include "target/mcp_http/McpSessionManager.h"
#endif

using namespace std;
using namespace Wt;

namespace
{
/** Enable/disable the submit button from the current input state.
 Leaves the button enabled while it is acting as a "Stop" button (a request is in flight);
 otherwise enables it iff there is non-empty text or a staged image (data-llm-has-image="1").
 */
WT_DECLARE_WT_MEMBER
(LlmEvalSend, Wt::JavaScriptFunction, "LlmEvalSend",
 function( ta, btn )
{
  if( !ta || !btn )
    return;
  if( btn.classList.contains('llm-stop-button') ){ btn.disabled = false; return; }
  var hasText = (ta.value.trim().length > 0);
  var hasImg = (ta.getAttribute('data-llm-has-image') === '1');
  btn.disabled = !(hasText || hasImg);
}
);

/** Wire the input text area: auto-grow, Enter-to-submit (Shift/Ctrl/Meta+Enter = newline),
 submit-button enable/disable, and clipboard-image paste.
   ta          - the input textarea element
   btn         - the submit/stop button element
   scrollArea  - the scrollable conversation container (overflow:auto)
   container   - the inner conversation container (whose bottom padding we adjust)
   wrapper     - the relatively-positioned wrapper (used to cap the grow height)
   inputArea   - the floating input area (its height sets the container's bottom padding)
   acceptImages- whether pasted images should be captured
   pasteCb     - function(base64, mimeType, widthPx, heightPx) to stage a pasted image
   enterNewline- if true, Enter just inserts a newline (phones, whose on-screen keyboards have no
                 Shift+Enter); the submit button is then the only way to send.
 */
WT_DECLARE_WT_MEMBER
(LlmSetupInput, Wt::JavaScriptFunction, "LlmSetupInput",
 function( ta, btn, scrollArea, container, wrapper, inputArea, acceptImages, pasteCb, enterNewline )
{
  if( !ta )
    return;

  function grow(){
    // Hidden (e.g., tab not yet shown): scrollHeight would be 0 and pin the height to 0px, so keep
    //  the current height; the ResizeObserver below re-runs this once it becomes visible.
    if( !ta.offsetWidth )
      return;
    ta.style.height = 'auto';
    var max = Math.max( 60, Math.round( (wrapper ? wrapper.clientHeight : 300) * 0.4 ) );
    var h = Math.min( ta.scrollHeight, max );
    ta.style.height = h + 'px';
    ta.style.overflowY = (ta.scrollHeight > max) ? 'auto' : 'hidden';
    // Just enough to clear the floating input (its 8px bottom inset + height) plus a 4px gap, so the
    // newest content sits right above the input instead of behind a tall empty pad.
    if( container && inputArea )
      container.style.paddingBottom = (inputArea.offsetHeight + 12) + 'px';
    if( scrollArea && scrollArea._llmAutoScroll )
      scrollArea.scrollTop = scrollArea.scrollHeight;
  }
  ta._llmGrow = grow;

  ta.addEventListener( 'input', function(){
    grow();
    Wt.WT.LlmEvalSend( ta, btn );
  } );

  ta.addEventListener( 'keydown', function( e ){
    if( !enterNewline && e.key === 'Enter' && !e.shiftKey && !e.ctrlKey && !e.metaKey && !e.altKey ){
      e.preventDefault();
      if( btn && !btn.disabled )
        btn.click();
    }
  } );

  if( acceptImages ){
    ta.addEventListener( 'paste', function( e ){
      var items = (e.clipboardData && e.clipboardData.items) ? e.clipboardData.items : [];
      for( var i = 0; i < items.length; ++i ){
        if( items[i].type && items[i].type.indexOf('image/') === 0 ){
          var blob = items[i].getAsFile();
          if( !blob )
            continue;
          e.preventDefault();
          var reader = new FileReader();
          reader.onload = function( ev ){
            var img = new Image();
            img.onload = function(){
              var maxSize = 2048;
              var w = img.naturalWidth, h = img.naturalHeight;
              if( w > maxSize || h > maxSize ){
                var scale = maxSize / Math.max(w, h);
                w = Math.round(w*scale); h = Math.round(h*scale);
              }
              var canvas = document.createElement('canvas');
              canvas.width = w; canvas.height = h;
              canvas.getContext('2d').drawImage( img, 0, 0, w, h );
              var outMime = (blob.type === 'image/png') ? 'image/png' : 'image/jpeg';
              var quality = (outMime === 'image/jpeg') ? 0.85 : undefined;
              var dataUrl = canvas.toDataURL( outMime, quality );
              var base64 = dataUrl.split(',')[1];
              if( pasteCb ) pasteCb( base64, outMime, w, h );
            };
            img.src = ev.target.result;
          };
          reader.readAsDataURL( blob );
          return;
        }
      }
    } );
  }

  // Re-size when the tool is shown or resized (wrapper), or when the stylesheet finishes loading or
  //  the image strip toggles (inputArea); the setup call below can run before either has settled.
  //  grow() is idempotent, so the height change it makes does not cause further notifications.
  if( window.ResizeObserver ){
    var ro = new ResizeObserver( function(){ grow(); } );
    if( wrapper )
      ro.observe( wrapper );
    if( inputArea )
      ro.observe( inputArea );
  }

  // Initialize sizing and button state.
  grow();
  Wt.WT.LlmEvalSend( ta, btn );
}
);

/** Set up auto-scroll on the conversation scroll area: stay pinned to the bottom as content
 streams in, unless the user has scrolled up; scrolling back to the bottom re-enables it. */
WT_DECLARE_WT_MEMBER
(LlmSetupAutoScroll, Wt::JavaScriptFunction, "LlmSetupAutoScroll",
 function( scrollArea, container )
{
  if( !scrollArea || !container )
    return;
  scrollArea._llmAutoScroll = true;
  scrollArea.addEventListener( 'scroll', function(){
    var atBottom = (scrollArea.scrollHeight - scrollArea.scrollTop - scrollArea.clientHeight) < 24;
    scrollArea._llmAutoScroll = atBottom;
  } );
  var obs = new MutationObserver( function(){
    if( scrollArea._llmAutoScroll )
      scrollArea.scrollTop = scrollArea.scrollHeight;
  } );
  obs.observe( container, { childList:true, subtree:true, characterData:true } );
}
);

  /** When the settings window is opened from within a (non-modal) window, the click that got us there
   raises that window to the front after the settings window is shown - which on phones (where both
   are full-screen) completely hides the settings window.  So, like SimpleDialog does, raise the
   settings window again on a delay.
   */
  void bring_to_front_after_click( AuxWindow *window )
  {
    window->doJavaScript( "setTimeout(function(){Wt.WT.AuxWindowBringToFront('" + window->id() + "');},5);" );
  }
}//namespace


LlmToolGui::LlmToolGui(InterSpec *viewer)
  : WContainerWidget(),
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
    m_logFilePath(),
    m_logWriteFailed(false),
    m_imagePreviewContainer(nullptr),
    m_inputBanner(nullptr),
    m_pasteImageSignal( this, "LlmPasteImage", false ),
    m_benchmarkRunner( nullptr ),
    m_jsSandboxBridge( nullptr ),
    m_menuIconHome( nullptr ),
    m_phoneHeaderTools(),
    m_phoneDockButton( nullptr ),
    m_phoneCanDock( false )
{
  if( !m_viewer )
    throw std::runtime_error("InterSpec instance cannot be null");

  // Phone styling is keyed off this class (on us, rather than m_root, so it survives UI rebuilds).
  if( m_viewer->isPhone() )
    addStyleClass( "LlmPhone" );

  wApp->useStyleSheet( "InterSpec_resources/LlmToolGui.css" );
  // The "not configured" panel and the settings window share styles + strings from these bundles.
  wApp->useStyleSheet( "InterSpec_resources/LlmConfigWindow.css" );
  m_viewer->useMessageResourceBundle( "LlmConfigWindow" );

  // Get marked/KaTeX loading now, so they are ready before the first response comes back.
  LlmMarkdown::load_resources( wApp );

  createRoot();

  // Connect the clipboard-image paste signal once; it persists across UI rebuilds (the input widgets
  // that emit it are recreated by initializeUI(), but this JSignal belongs to the tool).
  m_pasteImageSignal.connect( this, &LlmToolGui::handlePastedImage );

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
    configured = config && config->llmApi.isConfigured();
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
  // Final flush of the log, to capture anything since the last completed turn.
  writeConversationLog();

  // On phones our dock button and menu icon live in our windows header; don't leave them there
  //  (they would call back into this destroyed widget).  A no-op if the header is already gone.
  if( m_phoneHeaderTools )
  {
    const std::unique_ptr<WWidget> doomed = m_phoneHeaderTools->removeFromParent();
  }

  // Fail any outstanding MCP `assistant_submit_prompt` callbacks so a blocked MCP request thread
  //  doesn't hang forever waiting on a turn that can no longer complete.
  if( m_activeMcpCallback )
  {
    McpPromptCallback cb = std::move( m_activeMcpCallback );
    m_activeMcpCallback = McpPromptCallback();
    cb( string( "The assistant session was torn down before the prompt completed." ) );
  }

  for( PendingMcpPrompt &pending : m_mcpPromptQueue )
  {
    if( pending.callback )
      pending.callback( string( "The assistant session was torn down before the prompt ran." ) );
  }
  m_mcpPromptQueue.clear();

  // If the settings window is still open, delete it now so its onSaved callback (which captures
  // this) can never fire against a destroyed widget.
  if( m_configWindow )
  {
    // AuxWindow::make() gave ownership to wApp, so hand it back rather than `delete`ing a
    //  widget the app owns (that would double-free); the observing_ptr clears itself.
    AuxWindow::deleteAuxWindow( m_configWindow.get() );
    assert( !m_configWindow );
  }
}


void LlmToolGui::buildConfiguredUi( const std::shared_ptr<const LlmConfig> &config )
{
  try
  {
    m_llmInterface = std::make_shared<LlmInterface>(m_viewer, config);
    m_logFilePath.clear();

    // Connect to LLM interface response signals
    m_llmInterface->conversationFinished().connect(this, &LlmToolGui::handleConversationFinished);
    m_llmInterface->responseError().connect(this, &LlmToolGui::handleResponseError);
    m_llmInterface->activityChanged().connect( this, [this]( const LlmInterface::ActivityState state ){
      if( m_sendButton )
        m_sendButton->toggleStyleClass( "llm-tools-running",
                                        (state == LlmInterface::ActivityState::RunningTools) );

      // A request not sent through sendMessage() - e.g., Retry on an error - is still work the user
      //  may want to stop; the terminal signals return the input to idle as usual.  The user started
      //  it, so is there to answer any error it hits.
      if( (state != LlmInterface::ActivityState::Idle) && !m_isRequestPending )
      {
        m_isRequestPending = true;
        setInputEnabled( false );
        m_llmInterface->setUnattended( false );
      }
    } );

    // Connect to spectrum change signal to cancel pending requests when foreground spectrum changes
    m_spectrumChangedConnection =
      m_viewer->displayedSpectrumChanged().connect(this, &LlmToolGui::handleSpectrumChanged);

    initializeUI();
  }catch( std::exception &e )
  {
    m_llmInterface.reset();
    m_root->addNew<WText>( "Error initializing LLM assistant: " + string(e.what()) );
  }
}


void LlmToolGui::buildUnconfiguredUi()
{
  m_root->addStyleClass( "LlmToolGui" );

  WContainerWidget *panel = m_root->addNew<WContainerWidget>();
  panel->addStyleClass( "LlmNotConfigured" );

  WText *msg = panel->addNew<WText>( WString::tr("lcw-not-configured") );
  msg->addStyleClass( "LlmNotConfiguredMsg" );

  WPushButton *btn = panel->addNew<WPushButton>( WString::tr("lcw-configure-btn") );
  btn->clicked().connect( this, &LlmToolGui::openConfigWindow );

  if( m_viewer->isPhone() )
    addPhoneExpandButton( panel );
}


void LlmToolGui::showWelcomeMessage()
{
  if( m_welcomeMessage || !m_conversationContainer )
    return;

  m_welcomeMessage = m_conversationContainer->addNew<WContainerWidget>();
  m_welcomeMessage->addStyleClass( "LlmWelcome" );

  WText *title = m_welcomeMessage->addNew<WText>( WString::tr("lcw-welcome-title") );
  title->addStyleClass( "LlmWelcomeTitle" );

  const char * const summaryKey = m_viewer->isPhone() ? "lcw-welcome-summary-phone" : "lcw-welcome-summary";
  WText *summary = m_welcomeMessage->addNew<WText>( WString::tr(summaryKey) );
  summary->addStyleClass( "LlmWelcomeSummary" );

  // Empty WLink => an <a> with no href, so the browser does nothing on click while the
  // server-side clicked() still fires and opens the help window.
  WAnchor *learnMore = m_welcomeMessage->addNew<WAnchor>( WLink() );
  learnMore->setText( WString::tr("lcw-welcome-learn-more") );
  learnMore->addStyleClass( "LlmWelcomeLearnMore" );
  learnMore->clicked().connect( this, [](){ HelpSystem::createHelpWindow( "llm-assistant" ); } );
}//showWelcomeMessage()


void LlmToolGui::updateWelcomeMessage()
{
  if( !m_conversationContainer )   // unconfigured / torn down
    return;

  const int realMsgs = m_conversationContainer->count() - (m_welcomeMessage ? 1 : 0);
  if( realMsgs > 0 )
  {
    if( m_welcomeMessage )
    {
      m_conversationContainer->removeWidget( m_welcomeMessage );
      m_welcomeMessage = nullptr;
    }
  }else if( !m_welcomeMessage )
  {
    showWelcomeMessage();
  }
}//updateWelcomeMessage()


void LlmToolGui::clearConversationDisplays()
{
  if( m_conversationContainer )
    m_conversationContainer->clear();  // deletes the welcome panel and all interaction displays
  m_welcomeMessage = nullptr;          // pointer was just invalidated by clear()
}//clearConversationDisplays()


void LlmToolGui::createRoot()
{
  // A single-item outer layout holds m_root, which we swap between the configured chat UI and
  // the "not configured" panel.  A box layout, rather than a WGridLayout, as Wt4 implements box
  // layouts with CSS flex, while grid layouts use its JavaScript layout engine - whose bookkeeping
  // is left with stale entries when we are moved out of the tool-tab strip on phones (see
  // InterSpec::setLlmToolDocked), causing a client-side exception.
  //
  // When swapping UIs, the whole outer layout is replaced, rather than just m_root or m_root's
  // contents: Wt's flex layout only sizes things right when a layout and its contents are rendered
  // together.  Adding m_root to an already-rendered layout left its item without `overflow: hidden`
  // (so the chat UI grew past the tab), and setting a layout on an already-rendered m_root made it
  // `height: fit-content` (so it did not fill the tab) - either way the input field was not at the
  // bottom of the tab after the user first configured a provider.
  m_root = nullptr;
  WVBoxLayout *outer = setLayout( std::make_unique<WVBoxLayout>() );  // deletes any previous layout
  outer->setContentsMargins( 0, 0, 0, 0 );
  outer->setSpacing( 0 );
  m_root = outer->addWidget( std::make_unique<WContainerWidget>(), 1 );
}//void createRoot()


void LlmToolGui::resetRoot()
{
  // Tear down the configured-UI state.  The spectrum-changed slot must be disconnected here (not
  // just reassigned in buildConfiguredUi) - boost::signals2 reassignment drops the handle without
  // disconnecting, which would leave an orphaned slot firing into a rebuilt widget.
  m_spectrumChangedConnection.disconnect();

  // Staged images live in this object (not under m_root); drop them so they don't leak into the
  // next conversation after a reconfigure.
  m_stagedImages.clear();

  // On phones the menu icon lives in the window header, so would not go away with m_root.
  removePhoneMenuIcon();

  // Replace the whole outer layout (which deletes m_root and all its children).
  createRoot();

  // Child widgets owned by the old root are now gone; clear dangling pointers.
  m_layout = nullptr;
  m_menuIconHome = nullptr;
  m_conversationContainer = nullptr;
  m_welcomeMessage = nullptr;
  m_inputEdit = nullptr;
  m_sendButton = nullptr;
  m_menuIcon = nullptr;
  m_imagePreviewContainer = nullptr;
  m_inputBanner = nullptr;
  m_jsSandboxBridge = nullptr;
}


void LlmToolGui::openConfigWindow()
{
  if( m_configWindow )
    return;  // already open

  m_configWindow = AuxWindow::make<LlmConfigWindow>( m_viewer,
                                        std::bind( &LlmToolGui::handleConfigSaved, this ),
                                        std::bind( &LlmToolGui::hasConversationHistory, this ) );
  m_configWindow->finished().connect( this, [this](){ AuxWindow::deleteAuxWindow( m_configWindow.get() ); } );
  bring_to_front_after_click( m_configWindow.get() );
}


void LlmToolGui::refitPhoneWindows()
{
  if( m_configWindow )
    m_configWindow->refitPhoneFullScreen();
}//void refitPhoneWindows()


void LlmToolGui::setPhoneWindowHeader( WContainerWidget *header )
{
  // Take down the controls in any previous windows header, bringing the menu icon back to the
  //  conversation area.
  if( m_phoneHeaderTools )
  {
    if( m_menuIcon && m_menuIconHome && (m_menuIcon->parent() == m_phoneHeaderTools.get()) )
    {
      m_menuIcon->removeStyleClass( "LlmHeaderMenuIcon" );
      m_menuIconHome->addWidget( m_menuIcon->removeFromParent() );
    }

    const std::unique_ptr<WWidget> doomed = m_phoneHeaderTools->removeFromParent();
    m_phoneDockButton = nullptr;
  }//if( m_phoneHeaderTools )

  // While docked, CSS hides the menu icon, and shows the expand button.
  toggleStyleClass( "LlmDocked", !header );

  if( !header )
    return;

  // Floated right, as a unit, so the dock button stays left of the menu icon even when the icon is
  //  recreated (see initializeUI()).
  m_phoneHeaderTools = header->addNew<WContainerWidget>();
  m_phoneHeaderTools->addStyleClass( "LlmHeaderTools" );

  m_phoneDockButton = m_phoneHeaderTools->addNew<WPushButton>();
  m_phoneDockButton->addStyleClass( "LlmDockBtn InvertInDark" );
  m_phoneDockButton->setToolTip( WString::tr("lcw-phone-dock-tt") );
  m_phoneDockButton->setHidden( !m_phoneCanDock );
  m_phoneDockButton->clicked().connect( this, [this](){ m_viewer->setLlmToolDocked( true ); } );

  movePhoneMenuIcon();
}//void setPhoneWindowHeader( WContainerWidget *header )


void LlmToolGui::setPhoneCanDock( const bool canDock )
{
  m_phoneCanDock = canDock;
  if( m_phoneHeaderTools && m_phoneDockButton )
    m_phoneDockButton->setHidden( !canDock );
}//void setPhoneCanDock( const bool canDock )


void LlmToolGui::addPhoneExpandButton( WContainerWidget *parent )
{
  WPushButton *expand = parent->addNew<WPushButton>();
  expand->addStyleClass( "LlmExpandBtn InvertInDark" );
  expand->setIcon( "InterSpec_resources/images/llm_expand.svg" );
  expand->setToolTip( WString::tr("lcw-phone-expand-tt") );
  expand->clicked().connect( this, [this](){ m_viewer->setLlmToolDocked( false ); } );
}//void addPhoneExpandButton( WContainerWidget *parent )


void LlmToolGui::movePhoneMenuIcon()
{
  if( !m_phoneHeaderTools || !m_menuIcon || (m_menuIcon->parent() == m_phoneHeaderTools.get()) )
    return;

  // The icons popup menu is owned by, and connected to, the button itself, so it comes along.
  m_menuIcon->addStyleClass( "LlmHeaderMenuIcon" );
  m_phoneHeaderTools->addWidget( m_menuIcon->removeFromParent() );
}//void movePhoneMenuIcon()


void LlmToolGui::removePhoneMenuIcon()
{
  if( m_menuIcon && m_phoneHeaderTools && (m_menuIcon->parent() == m_phoneHeaderTools.get()) )
  {
    const std::unique_ptr<WWidget> doomed = m_menuIcon->removeFromParent();
    m_menuIcon = nullptr;
  }
}//void removePhoneMenuIcon()


void LlmToolGui::openConfigWindowToImport( const std::string &configFilePath )
{
  // Replace any window that is already open so the imported config is what the user sees.
  if( m_configWindow )
  {
    AuxWindow::deleteAuxWindow( m_configWindow.get() );  // deletes the window; does not emit finished()
    assert( !m_configWindow );   //observing_ptr auto-clears
  }

  m_configWindow = AuxWindow::make<LlmConfigWindow>( m_viewer,
                                        std::bind( &LlmToolGui::handleConfigSaved, this ),
                                        std::bind( &LlmToolGui::hasConversationHistory, this ),
                                        configFilePath );
  m_configWindow->finished().connect( this, [this](){ AuxWindow::deleteAuxWindow( m_configWindow.get() ); } );
  bring_to_front_after_click( m_configWindow.get() );
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
    SimpleDialog *errDialog = SimpleDialog::make<SimpleDialog>( "Error Loading LLM Config" );
    errDialog->contents()->addNew<WText>( e.what() );
    errDialog->addButton( "Okay", WidgetUtils::ButtonRole::Affirm );
    return;
  }

  InterSpecServer::set_llm_config( config );

  bool configured = false;
  try
  {
    configured = config && config->llmApi.isConfigured();
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
        SimpleDialog *errDialog = SimpleDialog::make<SimpleDialog>( "Error Updating LLM" );
        errDialog->contents()->addNew<WText>( e.what() );
        errDialog->addButton( "Okay", WidgetUtils::ButtonRole::Affirm );
        return;
      }
      setInputEnabled( true );
      updateInputBanner();
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
    return (config && config->llmApi.isConfigured());
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

  // Create the main layout (phones are narrow enough that we keep the margins small)
  m_layout = m_root->setLayout( std::make_unique<WVBoxLayout>() );  // flex; see createRoot()
  const int margin = m_viewer->isPhone() ? 2 : 5;
  m_layout->setContentsMargins( margin, margin, margin, margin );
  m_layout->setSpacing( 5 );

  // Create wrapper container for scroll area and menu icon
  auto scrollWrapperOwned = std::make_unique<WContainerWidget>();
  WContainerWidget *scrollWrapper = scrollWrapperOwned.get();
  scrollWrapper->setStyleClass( "LlmScrollWrapper" );


  // Sandboxed JS bridge for the `run_javascript` LLM tool.  Hidden child
  // of this widget; its iframe is positioned off-screen by the JS helper.
  try
  {
    m_jsSandboxBridge = scrollWrapper->addNew<LlmJsSandboxBridge>();
  }catch( std::exception &e )
  {
    cerr << "Failed to construct LlmJsSandboxBridge: " << e.what() << endl;
    m_jsSandboxBridge = nullptr;
  }

  // Create menu icon for conversation-level actions (positioned in upper-right via CSS)
  m_menuIcon = scrollWrapper->addNew<WPushButton>();
  m_menuIconHome = scrollWrapper;
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn LlmMenuIcon" );
  // Note: deliberately do NOT preventPropagation() here (doing this lets menu open, but not close when you click the button).

  // Create popup menu
  PopupDivMenu *menu = makePopupMenu( m_menuIcon );

  // On phones we leave out the items that are only of use when developing the LLM tools (raw JSON,
  //  re-reading the config XML files, and benchmarks).
  const bool phone = m_viewer->isPhone();
  if( phone )
    menu->addPhoneBackItem( nullptr );

  // Add "Export Conversation as JSON" option
  if( !phone )
  {
    PopupDivMenuItem *exportJsonItem = menu->addMenuItem( "Export Conversation as JSON" );
    exportJsonItem->triggered().connect( this, &LlmToolGui::exportConversationJson );
  }

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
  if( !phone )
  {
    PopupDivMenuItem *resetConfigItem = menu->addMenuItem( "Reset LLM Config" );
    resetConfigItem->triggered().connect( this, &LlmToolGui::handleResetLlmConfig );
  }

  // Add "Help" option - opens the standard InterSpec help window on the LLM Assistant topic
  PopupDivMenuItem *helpItem = menu->addMenuItem( WString::tr("lcw-help-menu") );
  helpItem->triggered().connect( this, [](){ HelpSystem::createHelpWindow( "llm-assistant" ); } );

  // Add benchmark submenu - scan for *_llm_benchmark.xml files in test_datasets/
  if( !phone )
  {
    const vector<pair<string,string>> benchmarkFiles = LlmBenchmarkRunner::findAvailableBenchmarks();

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

  // On phones, the icon goes in the windows header bar, if we have been given it.
  movePhoneMenuIcon();

  // Create scrollable container for conversation displays.  Wt4 removed WScrollArea; the
  //  replacement is a plain container with overflow:auto.
  WContainerWidget *scrollArea = scrollWrapper->addNew<WContainerWidget>();
  scrollArea->addStyleClass("LlmConversationScrollArea");
  scrollArea->setOverflow( Overflow::Auto );

  m_conversationContainer = scrollArea->addNew<WContainerWidget>();
  m_conversationContainer->addStyleClass("LlmConversationContainer");

  // Add scroll wrapper to layout (fills the whole tool; the input area floats over it)
  m_layout->addWidget( std::move(scrollWrapperOwned), 1 );  // Make conversation display expand

  // Floating input area, laid over the bottom of the scroll wrapper (CSS position:absolute).  It is
  // a horizontal flex row: the input column (image previews + growing text area) on the left, and
  // the icon submit/stop button pinned to the right.
  // Its client-side behaviour (see setupInputJavaScript()) is set up whenever its element is created.
  LlmJsRenderedDiv *inputArea = scrollWrapper->addNew<LlmJsRenderedDiv>();
  inputArea->addStyleClass( "LlmInputArea" );

  // Phones: while docked in the tools tab strip, an "expand to full-screen" button leads the row.
  if( m_viewer->isPhone() )
    addPhoneExpandButton( inputArea );

  WContainerWidget *inputColumn = inputArea->addNew<WContainerWidget>();
  inputColumn->addStyleClass( "LlmInputColumn" );

  // Optional per-provider banner, directly above the text area (shown by updateInputBanner())
  m_inputBanner = inputColumn->addNew<WText>();
  m_inputBanner->setTextFormat( TextFormat::Plain );
  m_inputBanner->addStyleClass( "LlmInputBanner" );
  m_inputBanner->hide();

  // Image preview strip (hidden by default)
  m_imagePreviewContainer = inputColumn->addNew<WContainerWidget>();
  m_imagePreviewContainer->addStyleClass( "LlmStagedImagePreview" );
  m_imagePreviewContainer->hide();

  // Auto-growing text area (replaces the old single-line edit)
  m_inputEdit = inputColumn->addNew<WTextArea>();
  m_inputEdit->addStyleClass("llm-input-edit");
  m_inputEdit->setRows( 1 );
  m_inputEdit->setPlaceholderText( WString( "Type your question here..." ) );
  m_inputEdit->setAttributeValue( "autocomplete", "off" );

#if( BUILD_AS_OSX_APP || IOS )
  m_inputEdit->setAttributeValue("autocorrect", "off");
  m_inputEdit->setAttributeValue("spellcheck", "off");
#endif

  // Icon submit/stop button, pinned to the right of the text area.  The button is transparent/
  // borderless; the CSS inverts just the icon glyph for the dark theme (not the whole button, which
  // would also invert the ring colors).  The ring is an animated "working" outline around the stop
  // icon, shown only while the button is a Stop button (CSS); its color shows whether the LLM or
  // tools are being waited on (the "llm-tools-running" class).
  m_sendButton = inputArea->addNew<WPushButton>();
  m_sendButton->addStyleClass("llm-send-button");
  m_sendButton->setIcon( "InterSpec_resources/images/llm_submit_arrow.svg" );
  m_sendButton->setTextFormat( TextFormat::XHTML );
  m_sendButton->setText( "<svg class=\"LlmWorkingRing\" width=\"28\" height=\"28\" viewBox=\"0 0 28 28\">"
                         "<rect x=\"1\" y=\"1\" width=\"26\" height=\"26\" rx=\"6\" ry=\"6\" /></svg>" );

  // Connect signals (Enter-to-submit is handled client-side in setupInputJavaScript)
  m_sendButton->clicked().connect(this, &LlmToolGui::handleSendButton);

  // Wire client-side behavior: auto-grow, Enter submit, enable/disable, paste, and auto-scroll.
  setupInputJavaScript();

  updateInputBanner();

  // Show the getting-started panel when the conversation area is still empty.  Any subsequent
  // restore of saved history (setConversationHistory) reconciles this away.
  updateWelcomeMessage();
}


void LlmToolGui::updateInputBanner()
{
  if( !m_inputBanner )
    return;

  string banner;
  try
  {
    const shared_ptr<const LlmConfig> config = m_llmInterface ? m_llmInterface->config() : nullptr;
    if( config )
      banner = config->llmApi.activeProvider().banner;
  }catch( std::exception & )
  {
    // No active provider - leave the banner hidden
  }

  // The CSS clamps long text to two lines, so the tooltip carries the full text.
  m_inputBanner->setText( WString::fromUTF8(banner) );
  m_inputBanner->setToolTip( WString::fromUTF8(banner) );
  m_inputBanner->setHidden( banner.empty() );
}//void updateInputBanner()


void LlmToolGui::focusInput()
{
  if( m_inputEdit )
    m_inputEdit->setFocus(true);
}


void LlmToolGui::setupInputJavaScript()
{
  if( !m_inputEdit || !m_sendButton || !m_conversationContainer )
    return;

  LOAD_JAVASCRIPT(wApp, "LlmToolGui.cpp", "LlmToolGui", wtjsLlmEvalSend);
  LOAD_JAVASCRIPT(wApp, "LlmToolGui.cpp", "LlmToolGui", wtjsLlmSetupInput);
  LOAD_JAVASCRIPT(wApp, "LlmToolGui.cpp", "LlmToolGui", wtjsLlmSetupAutoScroll);

  const bool acceptImages = canAcceptImages();

  // When images are supported, hand the setup a JS callback that emits the paste signal back to us.
  const string pasteCb = acceptImages
    ? ("function(b64,mime,w,h){" + m_pasteImageSignal.createCall( {"b64","mime","w","h"} ) + ";}")
    : string("null");

  ostringstream js;
  js << "(function(){"
     << "var ta=" << m_inputEdit->jsRef() << ", btn=" << m_sendButton->jsRef()
     << ", container=" << m_conversationContainer->jsRef() << ";"
     << "if(!ta||!btn||!container) return;"
     << "var scrollArea=container.parentElement;"
     << "var wrapper=scrollArea?scrollArea.parentElement:null;"
     << "var inputArea=btn.parentElement;"
     << "Wt.WT.LlmSetupAutoScroll(scrollArea, container);"
     << "Wt.WT.LlmSetupInput(ta, btn, scrollArea, container, wrapper, inputArea, "
     << (acceptImages ? "true" : "false") << ", " << pasteCb << ", "
     << (m_viewer->isPhone() ? "true" : "false") << ");"
     << "})();";

  // Run by the input area whenever its element is (re-)created - the event listeners set up here
  //  go with the old elements when the tool is moved between the tool tabs and a window.
  LlmJsRenderedDiv * const inputArea = dynamic_cast<LlmJsRenderedDiv *>( m_sendButton->parent() );
  assert( inputArea );
  if( inputArea )
    inputArea->setRenderJs( js.str() );
  else
    doJavaScript( js.str() );
}//setupInputJavaScript()


void LlmToolGui::updateSendButtonEnabled()
{
  if( !m_inputEdit || !m_sendButton )
    return;

  const bool hasImg = !m_stagedImages.empty();

  ostringstream js;
  js << "(function(){var ta=" << m_inputEdit->jsRef() << ", btn=" << m_sendButton->jsRef() << ";"
     << "if(!ta||!btn) return;"
     << "ta.setAttribute('data-llm-has-image','" << (hasImg ? "1" : "0") << "');"
     << "if(Wt.WT.LlmEvalSend) Wt.WT.LlmEvalSend(ta, btn);"
     << "})();";
  doJavaScript( js.str() );
}//updateSendButtonEnabled()


void LlmToolGui::resetInputHeight()
{
  // Only call this after `m_inputEdit->setText("")`.
  if( !m_inputEdit )
    return;

  // Queued on the textarea itself so it runs after Wt applies the cleared value; otherwise grow()
  //  can measure the old (long) text and leave the area expanded.  Clearing `value` here too makes
  //  this independent of Wt's update ordering.
  ostringstream js;
  js << "(function(){var ta=" << m_inputEdit->jsRef() << ";"
     << "if(!ta) return;"
     << "ta.value='';"
     << "ta.style.height='';"
     << "if(ta._llmGrow) ta._llmGrow();"
     << "})();";
  m_inputEdit->doJavaScript( js.str() );
}//resetInputHeight()


void LlmToolGui::handlePastedImage( const string &base64Data, const string &mimeType,
                                    int widthPx, int heightPx )
{
  if( base64Data.empty() )
    return;

  stageImage( base64Data, mimeType, "Pasted image", widthPx, heightPx );
}//handlePastedImage(...)

void LlmToolGui::clearHistory()
{
  // Clear all interaction displays from the container, then re-show the getting-started panel
  clearConversationDisplays();
  updateWelcomeMessage();

  // Also clear the LLM interface history
  shared_ptr<LlmConversationHistory> history = m_llmInterface ? m_llmInterface->getHistory() : nullptr;
  if( history )
    history->clear();

  m_logFilePath.clear(); // The next conversation gets its own log file
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
                              const string &displayName, int widthPx, int heightPx,
                              const string &caption )
{
  StagedImage staged;
  staged.base64Data = base64Data;
  staged.mimeType = mimeType;
  staged.displayName = displayName;
  staged.caption = caption;
  staged.widthPx = widthPx;
  staged.heightPx = heightPx;
  m_stagedImages.push_back( std::move( staged ) );

  if( !m_imagePreviewContainer )
    return;

  // Create a thumbnail entry for this image
  const size_t index = m_stagedImages.size() - 1;

  WContainerWidget *entry = m_imagePreviewContainer->addNew<WContainerWidget>();
  entry->addStyleClass( "LlmStagedImageEntry" );

  // Create thumbnail using WMemoryResource (same pattern as UploadedImgDisplay)
  const string decoded = Wt::Utils::base64Decode( base64Data );
  vector<unsigned char> bytes( decoded.begin(), decoded.end() );
  auto thumbResource = std::make_shared<WMemoryResource>( mimeType );
  thumbResource->setData( bytes );

  WImage *thumb = entry->addNew<WImage>( WLink( thumbResource ) );
  thumb->addStyleClass( "LlmStagedImageThumb" );

  // Filename label
  WText *nameLabel = entry->addNew<WText>( WString::fromUTF8( displayName ) );
  nameLabel->addStyleClass( "LlmStagedImageName" );

  // Remove button
  WPushButton *removeBtn = entry->addNew<WPushButton>();
  removeBtn->setStyleClass( "LlmStagedImageRemove" );
  removeBtn->setText( WString::fromUTF8( "\xc3\x97" ) ); // multiplication sign as X
  removeBtn->clicked().connect( std::bind( [this]( size_t idx ){ removeStagedImage( idx ); }, index ) );

  m_imagePreviewContainer->show();

  // A staged image is enough to allow sending even with no typed text.
  updateSendButtonEnabled();
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
      updateSendButtonEnabled();
      return;
    }

    // Re-add all remaining thumbnails
    for( size_t i = 0; i < m_stagedImages.size(); ++i )
    {
      const StagedImage &img = m_stagedImages[i];

      WContainerWidget *entry = m_imagePreviewContainer->addNew<WContainerWidget>();
      entry->addStyleClass( "LlmStagedImageEntry" );

      const string decoded = Wt::Utils::base64Decode( img.base64Data );
      vector<unsigned char> bytes( decoded.begin(), decoded.end() );
      auto thumbResource = std::make_shared<WMemoryResource>( img.mimeType );
      thumbResource->setData( bytes );

      WImage *thumb = entry->addNew<WImage>( WLink( thumbResource ) );
      thumb->addStyleClass( "LlmStagedImageThumb" );

      WText *nameLabel = entry->addNew<WText>( WString::fromUTF8( img.displayName ) );
      nameLabel->addStyleClass( "LlmStagedImageName" );

      WPushButton *removeBtn = entry->addNew<WPushButton>();
      removeBtn->setStyleClass( "LlmStagedImageRemove" );
      removeBtn->setText( WString::fromUTF8( "\xc3\x97" ) );
      removeBtn->clicked().connect( std::bind( [this]( size_t idx ){ removeStagedImage( idx ); }, i ) );
    }
  }

  updateSendButtonEnabled();
}//void removeStagedImage(...)


void LlmToolGui::clearStagedImages()
{
  m_stagedImages.clear();
  if( m_imagePreviewContainer )
  {
    m_imagePreviewContainer->clear();
    m_imagePreviewContainer->hide();
  }

  updateSendButtonEnabled();
}//void clearStagedImages()


void LlmToolGui::handleInputSubmit()
{
  handleSendButton(); // Same logic as clicking send button
}

void LlmToolGui::handleSendButton()
{
  // While a request is in flight the Send button acts as a Stop button.  (cancelCurrentRequest() also
  //  returns the input to idle if nothing turns out to be outstanding.)
  if( m_isRequestPending )
  {
    cancelCurrentRequest();
    return;
  }

  if( !m_inputEdit )
    return;

  const string message = m_inputEdit->text().toUTF8();
  if( message.empty() && m_stagedImages.empty() )
    return;

  // Clear the input field and shrink it back to its minimal one-row height.  (The button's
  // enabled/Stop state is set authoritatively by sendMessage() -> setInputEnabled(false) below.)
  m_inputEdit->setText("");
  resetInputHeight();

  // Send the message
  sendMessage(message);
}

bool LlmToolGui::sendMessage(const std::string& message)
{
  if (!m_llmInterface) {
    cout << "Error: LLM interface not available" << endl;
    return false;
  }

  // Don't allow sending if a request is already pending
  if (m_isRequestPending) {
    cout << "Warning: Request already pending, ignoring new message" << endl;
    return false;
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
      img.caption = std::move( staged.caption );
      img.widthPx = staged.widthPx;
      img.heightPx = staged.heightPx;
      images.push_back( std::move( img ) );
    }
    clearStagedImages();

    // An MCP client or benchmark cant answer a sub-agent's error panel, so have the agent that invoked
    //  it carry on instead.
    m_llmInterface->setUnattended( static_cast<bool>( m_activeMcpCallback ) || isBenchmarkRunning() );

    // Send message to LLM - this creates a new LlmInteraction in the history
    shared_ptr<LlmInteraction> convo = m_llmInterface->sendUserMessage( message, std::move( images ) );

    // Get the conversation history to find the new interaction
    if( convo )
    {
      // Create display widget for this interaction
      LlmInteractionDisplay *interactionDisplay =
        m_conversationContainer->addNew<LlmInteractionDisplay>( convo, m_llmInterface, 0 );

      // Remove the getting-started panel now that there is a real interaction
      updateWelcomeMessage();

      // Sending resumes auto-scroll and jumps to the bottom.
      const string js = "setTimeout(function() {"
      "  var scrollArea = " + m_conversationContainer->jsRef() + ".parentElement;"
      "  if (scrollArea){ scrollArea._llmAutoScroll = true; scrollArea.scrollTop = scrollArea.scrollHeight; }"
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

    // Show error in a dialog - except to a benchmark, which handles (and logs) the failure itself,
    //  and may be running with no one to dismiss the dialog.
    if( !isBenchmarkRunning() )
    {
      SimpleDialog *errorDialog = SimpleDialog::make<SimpleDialog>( "Error" );
      errorDialog->contents()->addNew<WText>( "Failed to send message: " + string(e.what()) );
      WPushButton *okBtn = errorDialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
      okBtn->clicked().connect( errorDialog, &SimpleDialog::accept );
    }

    return false;
  }

  return true;
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

  writeConversationLog();

  // Resolve any in-flight MCP `assistant_submit_prompt` call, then start the next queued prompt.
  resolveActiveMcpPrompt( true );

  // The display is automatically updated via signals in LlmInteractionDisplay
  // No need to manually refresh
}

void LlmToolGui::handleResponseError()
{
  // The request is over, so return the input to its idle state: this swaps the button back to the
  // submit arrow (dropping the now-meaningless "Stop" icon), restores the normal placeholder, and
  // re-enables typing.  The errored turn's own Retry / Continue Anyway buttons remain available.
  if( m_isRequestPending )
  {
    m_isRequestPending = false;
    cout << "Error response received, re-enabling input" << endl;
    setInputEnabled( true );
  }

  writeConversationLog();

  // If this errored turn was an MCP-submitted prompt, resolve its callback with an error and keep
  //  draining the queue so automation doesn't get stuck on the human retry flow.
  resolveActiveMcpPrompt( false );

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
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "No Conversation" );
    WText *msg = dialog->contents()->addNew<WText>( "There is no conversation history to export." );
    WPushButton *okBtn = dialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
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
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "No Active Conversation" );
    WText *msg = dialog->contents()->addNew<WText>( "There is no active conversation to export." );
    WPushButton *okBtn = dialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
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
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Export Failed" );
    WText *msg = dialog->contents()->addNew<WText>( string("Failed to build export JSON: ") + e.what() );
    WPushButton *okBtn = dialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
    okBtn->clicked().connect( dialog, &SimpleDialog::accept );
    return;
  }

  const string jsonStr = requestJson.dump(2);  // Pretty print with 2-space indent
#else
  // Legacy format: bare messages array only (no system prompt, tools, or model params).
  const nlohmann::json conversationsJson = history->toApiFormat();
  const string jsonStr = conversationsJson.dump(2);  // Pretty print with 2-space indent
#endif  //LLM_EXPORT_FULL_API_JSON

  // Present the download as an anchor the user clicks, rather than redirecting the page to the
  //  resource - a redirect just navigates away from the app in the native (WKWebView) builds.
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Export Conversation" );

  // In Wt4 a WLink holds a shared_ptr to the resource, so the link itself keeps it alive for as
  //  long as the download anchor exists.
  auto resource = std::make_shared<WMemoryResource>( "application/json" );
  resource->setData( reinterpret_cast<const unsigned char*>(jsonStr.data()), jsonStr.size() );
  resource->suggestFileName( "llm_conversation.json" );
  resource->setDispositionType( ContentDisposition::Attachment );

  WLink dlLink( resource );
  dlLink.setTarget( LinkTarget::NewWindow );
  WAnchor *download = dialog->contents()->addNew<WAnchor>( dlLink );
  download->setStyleClass( "LinkBtn DownloadLink" );
  download->setText( "Download conversation JSON" );

  WPushButton *closeBtn = dialog->addButton( "Close", WidgetUtils::ButtonRole::Dismiss );
  closeBtn->clicked().connect( dialog, &SimpleDialog::accept );
}

void LlmToolGui::handleClearConversation()
{
  // Show confirmation dialog
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Clear Conversation" );
  WText *msg = dialog->contents()->addNew<WText>( "Are you sure you want to clear the entire conversation history? This cannot be undone." );

  WPushButton *yesBtn = dialog->addButton( "Yes, Clear", WidgetUtils::ButtonRole::Destructive );
  WPushButton *noBtn = dialog->addButton( "Cancel", WidgetUtils::ButtonRole::Dismiss );

  yesBtn->clicked().connect( std::bind( [this, dialog]() {
    clearHistory();
    dialog->accept();
  }));

  noBtn->clicked().connect( dialog, &SimpleDialog::accept );
}


void LlmToolGui::handleResetLlmConfig()
{
  // Show confirmation dialog before resetting
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Reset LLM Config" );
  dialog->contents()->addNew<WText>( "This will re-read the LLM configuration files and clear the current conversation. Continue?" );

  WPushButton *yesBtn = dialog->addButton( "Yes, Reset", WidgetUtils::ButtonRole::Destructive );
  WPushButton *noBtn = dialog->addButton( "Cancel", WidgetUtils::ButtonRole::Dismiss );

  yesBtn->clicked().connect( std::bind( [this, dialog]() {
    dialog->accept();

    // Try to load the config from disk; if it fails, show error and bail out
    shared_ptr<LlmConfig> config;

    try
    {
      config = LlmConfig::load();
    }catch( std::exception &e )
    {
      SimpleDialog *errDialog = SimpleDialog::make<SimpleDialog>( "Error Loading LLM Config" );
      errDialog->contents()->addNew<WText>( e.what() );
      errDialog->addButton( "Okay", WidgetUtils::ButtonRole::Affirm );
      return;
    }

    if( !config || !config->llmApi.isConfigured() )
    {
      SimpleDialog *errDialog = SimpleDialog::make<SimpleDialog>( "Error Loading LLM Config" );
      errDialog->contents()->addNew<WText>( "LLM API is not enabled in the configuration." );
      errDialog->addButton( "Okay", WidgetUtils::ButtonRole::Affirm );
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
      SimpleDialog *errDialog = SimpleDialog::make<SimpleDialog>( "Error Resetting LLM" );
      errDialog->contents()->addNew<WText>( e.what() );
      errDialog->addButton( "Okay", WidgetUtils::ButtonRole::Affirm );
      return;
    }

    setInputEnabled( true );
    updateInputBanner();
  }));

  noBtn->clicked().connect( dialog, &SimpleDialog::accept );
}


void LlmToolGui::handleCompactConversation()
{
  if( !m_llmInterface )
    return;

  if( m_isRequestPending )
  {
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Cannot Compact" );
    dialog->contents()->addNew<WText>( "Please wait for the current request to finish before compacting." );
    dialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
    return;
  }

  try
  {
    shared_ptr<LlmInteraction> summarizationConvo = m_llmInterface->triggerManualCompaction();

    if( !summarizationConvo )
    {
      SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Cannot Compact" );
      dialog->contents()->addNew<WText>( "Not enough conversation history to compact." );
      dialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
      return;
    }

    // Create display widget for the summarization conversation
    m_conversationContainer->addNew<LlmInteractionDisplay>( summarizationConvo, m_llmInterface, 0 );

    // Remove the getting-started panel now that there is a real interaction
    updateWelcomeMessage();

    // Resume auto-scroll and jump to the bottom.
    const string js = "setTimeout(function() {"
    "  var scrollArea = " + m_conversationContainer->jsRef() + ".parentElement;"
    "  if (scrollArea){ scrollArea._llmAutoScroll = true; scrollArea.scrollTop = scrollArea.scrollHeight; }"
    "}, 100);";
    doJavaScript( js );

    // Disable input while compaction is in-flight
    setInputEnabled( false );
    m_isRequestPending = true;
  }catch( const std::exception &e )
  {
    cerr << "Error triggering conversation compaction: " << e.what() << endl;

    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Compaction Error" );
    dialog->contents()->addNew<WText>( string("Failed to compact conversation: ") + e.what() );
    dialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
  }
}


void LlmToolGui::writeConversationLog()
{
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  try
  {
    if( !m_llmInterface )
      return;

    const shared_ptr<const LlmConfig> config = m_llmInterface->config();
    if( !config || !config->llmApi.logConversations )
      return;

    const shared_ptr<LlmConversationHistory> history = m_llmInterface->getHistory();
    if( !history || history->isEmpty() )
      return;

    const vector<shared_ptr<LlmInteraction>> &convos = history->getConversations();
    assert( !convos.empty() && convos.front() );
    // After compaction, front() is the summary placeholder, whose timestamp is that of the earliest
    //  summarized message - so this is always when the conversation started.
    const chrono::system_clock::time_point start_time = convos.front()->timestamp;

    if( m_logFilePath.empty() )
    {
#if( BUILD_AS_MCP_HTTP_SERVER )
      // A session of an MCP agent server logs into its own directory, which goes away with it.
      const std::shared_ptr<McpSessionManager> manager = McpSessionManager::instance();
      const std::optional<McpSessionManager::SessionInfo> session
                    = manager ? manager->session( McpSessionManager::currentHandle() ) : std::nullopt;
      if( !session.has_value() || session->dir.empty() )
        throw runtime_error( "no session directory to log to" );
      const string log_dir = SpecUtils::append_path( session->dir, "logs" );
#else
      const string base_dir = InterSpec::writableDataDirectory(); //throws if not set
      const string log_dir = SpecUtils::append_path( base_dir, "llm_logs" );
#endif
      if( SpecUtils::create_directory( log_dir ) == 0 ) //-1 already existed, 1 created
        throw runtime_error( "could not create directory '" + log_dir + "'" );

      // A random id, rather than the session id, as one session may have many conversations,
      //  and we dont want session ids on disk.
      string path;
      for( size_t attempt = 0; path.empty() || SpecUtils::is_file( path ); ++attempt )
      {
        if( attempt > 10 )
          throw runtime_error( "could not find unused log file name" );
        path = SpecUtils::append_path( log_dir,
                       LlmConversationHistory::logFileName( start_time, Wt::WRandom::generateId(6) ) );
      }
      m_logFilePath = path;
    }//if( m_logFilePath.empty() )

    vector<pair<string,string>> attribs;
    const char * const iso_fmt = "%Y-%m-%dT%H:%M:%S%z"; //local time, with UTC offset
    attribs.emplace_back( "started", LlmConversationHistory::formatLocalTime( start_time, iso_fmt ) );
    attribs.emplace_back( "lastWritten",
                          LlmConversationHistory::formatLocalTime( chrono::system_clock::now(), iso_fmt ) );
#ifdef InterSpec_VERSION
    attribs.emplace_back( "InterSpecVersion", InterSpec_VERSION );
#endif
    try
    {
      attribs.emplace_back( "model", config->llmApi.model() );
      attribs.emplace_back( "apiFormat", LlmConfig::LlmApi::apiFormatToString( config->llmApi.apiFormat() ) );
    }catch( std::exception & )
    {
      // No active provider/model - not important for the log
    }

    history->writeLogFile( m_logFilePath, attribs );
    m_logWriteFailed = false;
  }catch( std::exception &e )
  {
    if( !m_logWriteFailed )
      cerr << "Failed to write LLM conversation log: " << e.what() << endl;
    m_logWriteFailed = true;
  }
#endif
}//void writeConversationLog()


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

  // While working, keep the button clickable but swap its icon to a "Stop" square so the user can
  // cancel the request (handleSendButton() routes to cancelAll() when a request is pending).
  if( m_sendButton )
  {
    m_sendButton->setEnabled( true );
    if( enabled )
    {
      m_sendButton->setIcon( "InterSpec_resources/images/llm_submit_arrow.svg" );
      m_sendButton->removeStyleClass( "llm-stop-button" );
    }else
    {
      m_sendButton->setIcon( "InterSpec_resources/images/llm_stop.svg" );
      m_sendButton->addStyleClass( "llm-stop-button" );
    }
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

  // Sync the button's client-side clickability.  When returning to idle, re-apply the empty/non-empty
  // enable state.  While working, force the Stop button enabled directly (don't route through
  // LlmEvalSend, whose result would depend on the "llm-stop-button" class DOM diff having already been
  // applied client-side) so the user can always cancel an in-flight request.
  if( enabled )
  {
    updateSendButtonEnabled();
  }else if( m_sendButton )
  {
    doJavaScript( "(function(){var btn=" + m_sendButton->jsRef() + ";if(btn) btn.disabled=false;})();" );
  }
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

  if( m_activeMcpCallback )
    m_activeMcpCancelled = true;

  // Stop the turn itself (in-flight requests and deferred tool rounds), or it would keep calling
  //  tools and share the history with the next turn; cancelAll() emits conversationFinished, so
  //  handleConversationFinished() re-enables input and resolves any MCP prompt (as cancelled).
  if( m_llmInterface && m_llmInterface->hasActiveRequests() )
  {
    m_llmInterface->cancelAll();
    return;
  }

  m_isRequestPending = false;
  setInputEnabled(true);
  writeConversationLog();

  // If an MCP prompt was in flight, resolve it so the SSE stream ends and the queue can drain to
  //  the next entry.
  resolveActiveMcpPrompt( false );
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

  m_logFilePath.clear(); // A restored history is logged to a new file at its next turn end

  // Clear current display
  clearConversationDisplays();

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
        m_conversationContainer->addNew<LlmInteractionDisplay>( conv, m_llmInterface, 0 );
    }
  }else
  {
    // Clear the history if no valid history provided
    llmHistory->clear();
  }

  // Re-show the getting-started panel if nothing was restored; remove it otherwise
  updateWelcomeMessage();
}

void LlmToolGui::clearConversationHistory()
{
  if( !m_llmInterface )
    return;

  shared_ptr<LlmConversationHistory> history = m_llmInterface ? m_llmInterface->getHistory() : nullptr;
  if( history )
    history->clear();

  m_logFilePath.clear();

  // Clear the display, then re-show the getting-started panel
  clearConversationDisplays();
  updateWelcomeMessage();
}


bool LlmToolGui::submitMessageAsUser( const string &message )
{
  if( !m_inputEdit )
    return false;

  // During benchmarks, force-clear any pending/disabled state from previous errors
  if( isBenchmarkRunning() )
  {
    m_isRequestPending = false;
    setInputEnabled( true );
  }

  // As handleSendButton() does with what the user typed, but reporting whether it was sent.
  if( m_isRequestPending || (message.empty() && m_stagedImages.empty()) )
    return false;

  m_inputEdit->setText( "" );
  resetInputHeight();

  return sendMessage( message );
}


void LlmToolGui::queuePromptForMcp( const string &prompt, McpPromptCallback callback,
                                    McpProgressCallback progressCallback,
                                    vector<LlmToolCall::ImageContent> images,
                                    const string &id )
{
  if( !callback )
    return;

  if( prompt.empty() )
  {
    callback( string( "Empty prompt: nothing to submit." ) );
    return;
  }

  if( !m_llmInterface )
  {
    callback( string( "The LLM assistant is not configured for this session." ) );
    return;
  }

  if( isBenchmarkRunning() )
  {
    callback( string( "A benchmark is currently running; cannot submit a prompt right now." ) );
    return;
  }

  if( !images.empty() && !canAcceptImages() )
  {
    callback( string( "The configured model does not accept images." ) );
    return;
  }

  m_mcpPromptQueue.push_back( PendingMcpPrompt{ prompt, std::move( callback ), std::move( progressCallback ),
                                                std::move( images ), id } );
  drainMcpPromptQueue();
}


bool LlmToolGui::cancelMcpPrompt( const string &id )
{
  if( id.empty() )
    return false;

  for( auto iter = begin(m_mcpPromptQueue); iter != end(m_mcpPromptQueue); ++iter )
  {
    if( iter->id == id )
    {
      McpPromptCallback cb = std::move( iter->callback );
      m_mcpPromptQueue.erase( iter );
      if( cb )
        cb( string( "The prompt was cancelled before it ran." ) );
      return true;
    }
  }//for( loop over queued prompts )

  if( m_activeMcpCallback && (m_activeMcpId == id) )
  {
    m_activeMcpCancelled = true;
    cancelCurrentRequest();  //resolves the active prompt, then drains the queue

    // (The queue may have started the next prompt, so check the id.)
    if( m_activeMcpCallback && (m_activeMcpId == id) )
      resolveActiveMcpPrompt( false );
    return true;
  }

  return false;
}//bool cancelMcpPrompt( const string &id )


bool LlmToolGui::isTurnInProgress() const
{
  return m_isRequestPending || static_cast<bool>(m_activeMcpCallback) || !m_mcpPromptQueue.empty();
}


void LlmToolGui::drainMcpPromptQueue()
{
  // Only dispatch when no MCP request is already in flight and the conversation is idle.
  if( m_activeMcpCallback || m_isRequestPending || m_mcpPromptQueue.empty() )
    return;

  // A benchmark may have started after these prompts were queued; fail them rather than interfere
  //  with the benchmark or hang waiting for it (queuePromptForMcp already refuses at enqueue time).
  if( isBenchmarkRunning() )
  {
    for( PendingMcpPrompt &pending : m_mcpPromptQueue )
    {
      if( pending.callback )
        pending.callback( string( "A benchmark started; queued prompt was dropped." ) );
    }
    m_mcpPromptQueue.clear();
    return;
  }

  PendingMcpPrompt next = std::move( m_mcpPromptQueue.front() );
  m_mcpPromptQueue.pop_front();
  m_activeMcpCallback = std::move( next.callback );
  m_activeProgressCallback = std::move( next.progressCallback );
  m_activeMcpId = next.id;
  m_activeMcpCancelled = false;

  // Staged right before sending, since sending takes every staged image along.
  for( const LlmToolCall::ImageContent &img : next.images )
    stageImage( img.base64Data, img.mimeType, img.caption.empty() ? string("image") : img.caption,
                img.widthPx, img.heightPx, img.caption );

  submitMessageAsUser( next.prompt );

  // If the turn did not actually start (e.g. input unavailable / validation refused it),
  //  fail fast so the waiting MCP call doesn't hang, and move on to the next queued prompt.
  if( !m_isRequestPending )
  {
    resolveActiveMcpPrompt( false );
  }else if( m_activeProgressCallback && m_llmInterface )
  {
    // Connect to the conversation's responseAdded signal for streaming progress
    const std::shared_ptr<LlmConversationHistory> history = m_llmInterface->getHistory();
    if( history )
    {
      const std::vector<std::shared_ptr<LlmInteraction>> &conversations = history->getConversations();
      if( !conversations.empty() )
        connectMcpProgressToConversation( conversations.back() );
    }
  }
}


void LlmToolGui::resolveActiveMcpPrompt( bool success )
{
  disconnectMcpProgress();

  if( m_activeMcpCallback )
  {
    McpPromptCallback cb = std::move( m_activeMcpCallback );
    m_activeMcpCallback = McpPromptCallback();
    const bool cancelled = m_activeMcpCancelled;
    m_activeMcpId.clear();
    m_activeMcpCancelled = false;

    // (A cancelled turn ends through conversationFinished, i.e. "success", so check cancelled first.)
    if( cancelled )
      cb( string( "The assistant's turn was cancelled." ) );
    else if( success )
      cb( buildMcpPromptResult() );
    else
      cb( string( "The assistant's turn ended with an error before producing a final response." ) );
  }

  drainMcpPromptQueue();
}


nlohmann::json LlmToolGui::buildMcpPromptResult() const
{
  nlohmann::json result;
  result["response"] = "";
  result["tool_calls"] = nlohmann::json::array();

  if( !m_llmInterface )
    return result;

  const shared_ptr<LlmConversationHistory> history = m_llmInterface->getHistory();
  if( !history )
    return result;

  const vector<shared_ptr<LlmInteraction>> &conversations = history->getConversations();
  if( conversations.empty() )
    return result;

  const shared_ptr<LlmInteraction> &convo = conversations.back();
  if( !convo )
    return result;

  // Final assistant text: the last FinalLlmResponse in the most recent conversation.
  for( int i = static_cast<int>( convo->responses.size() ) - 1; i >= 0; --i )
  {
    const shared_ptr<LlmInteractionTurn> &turn = convo->responses[i];
    if( turn && (turn->type() == LlmInteractionTurn::Type::FinalLlmResponse) )
    {
      const LlmInteractionFinalResponse *finalResp
                       = dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
      if( finalResp )
      {
        result["response"] = finalResp->content();
        break;
      }
    }
  }

  // Compact tool-call trace: the names of tools requested this conversation, in order.
  for( const shared_ptr<LlmInteractionTurn> &turn : convo->responses )
  {
    if( !turn || (turn->type() != LlmInteractionTurn::Type::ToolCall) )
      continue;

    const LlmToolRequest *req = dynamic_cast<const LlmToolRequest *>( turn.get() );
    if( !req )
      continue;

    for( const LlmToolCall &call : req->toolCalls() )
      result["tool_calls"].push_back( call.toolName );
  }

  // Questions the assistant asked the GUI's user (ask_user_question), with their answers, so an
  //  MCP client can surface them.  (An MCP HTTP server has no such user, nor the tool.)
  nlohmann::json questions = nlohmann::json::array();
  for( const shared_ptr<LlmInteractionTurn> &turn : convo->responses )
  {
    if( !turn || (turn->type() != LlmInteractionTurn::Type::ToolResult) )
      continue;

    const LlmToolResults *results = dynamic_cast<const LlmToolResults *>( turn.get() );
    if( !results )
      continue;

    for( const LlmToolCall &call : results->toolCalls() )
    {
      if( call.toolName != "ask_user_question" )
        continue;
      nlohmann::json entry;
      entry["question"] = call.toolParameters.value( "question", nlohmann::json("") );
      entry["answer"] = call.content;
      questions.push_back( std::move(entry) );
    }
  }//for( loop over turns )

  if( !questions.empty() )
    result["questions_asked"] = std::move( questions );

  return result;
}


nlohmann::json LlmToolGui::buildMcpAgentResult( const shared_ptr<LlmInteraction> &convo ) const
{
  nlohmann::json result;
  result["response"] = "";
  result["tool_calls"] = nlohmann::json::array();

  if( !convo )
    return result;

  for( int i = static_cast<int>( convo->responses.size() ) - 1; i >= 0; --i )
  {
    const shared_ptr<LlmInteractionTurn> &turn = convo->responses[i];
    if( turn && (turn->type() == LlmInteractionTurn::Type::FinalLlmResponse) )
    {
      const LlmInteractionFinalResponse *finalResp
                       = dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
      if( finalResp )
      {
        result["response"] = finalResp->content();
        break;
      }
    }
  }

  for( const shared_ptr<LlmInteractionTurn> &turn : convo->responses )
  {
    if( !turn || (turn->type() != LlmInteractionTurn::Type::ToolCall) )
      continue;

    const LlmToolRequest *req = dynamic_cast<const LlmToolRequest *>( turn.get() );
    if( !req )
      continue;

    for( const LlmToolCall &call : req->toolCalls() )
      result["tool_calls"].push_back( call.toolName );
  }

  return result;
}


nlohmann::json LlmToolGui::buildProgressFromTurn( const shared_ptr<LlmInteractionTurn> &turn ) const
{
  nlohmann::json progress;
  if( !turn )
    return progress;

  switch( turn->type() )
  {
    case LlmInteractionTurn::Type::FinalLlmResponse:
    {
      progress["type"] = "assistant_response";
      const LlmInteractionFinalResponse *resp
        = dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
      if( resp )
        progress["content"] = resp->content();
      break;
    }

    case LlmInteractionTurn::Type::ToolCall:
    {
      progress["type"] = "tool_call";
      const LlmToolRequest *req = dynamic_cast<const LlmToolRequest *>( turn.get() );
      if( req )
      {
        nlohmann::json tools = nlohmann::json::array();
        for( const LlmToolCall &call : req->toolCalls() )
          tools.push_back( call.toolName );
        progress["tools"] = std::move( tools );
      }
      break;
    }

    case LlmInteractionTurn::Type::ToolResult:
    {
      progress["type"] = "tool_result";
      const LlmToolResults *results = dynamic_cast<const LlmToolResults *>( turn.get() );
      if( results )
      {
        nlohmann::json tools = nlohmann::json::array();
        for( const LlmToolCall &call : results->toolCalls() )
        {
          nlohmann::json entry;
          entry["name"] = call.toolName;
          if( call.status == LlmToolCall::CallStatus::Success )
          {
            // Include a truncated version of the result for progress reporting
            const string &content = call.content;
            if( content.size() > 500 )
              entry["content"] = content.substr( 0, 500 ) + "...";
            else
              entry["content"] = content;
          }else if( call.status == LlmToolCall::CallStatus::Error )
          {
            entry["error"] = call.content;
          }

          if( call.imageContent.has_value() )
          {
            entry["image"] = nlohmann::json{
              {"mimeType", call.imageContent->mimeType},
              {"data", call.imageContent->base64Data}
            };
          }

          tools.push_back( std::move( entry ) );
        }
        progress["results"] = std::move( tools );
      }
      break;
    }

    case LlmInteractionTurn::Type::Error:
    {
      progress["type"] = "error";
      const LlmInteractionError *err = dynamic_cast<const LlmInteractionError *>( turn.get() );
      if( err )
        progress["content"] = err->errorMessage();
      break;
    }

    default:
      progress["type"] = "other";
      break;
  }//switch( turn->type() )

  return progress;
}


void LlmToolGui::connectMcpProgressToConversation( const shared_ptr<LlmInteraction> &convo )
{
  if( !convo || !m_activeProgressCallback )
    return;

  m_mcpResponseAddedConn = convo->responseAdded.connect(
    [this]( shared_ptr<LlmInteractionTurn> turn )
    {
      if( !m_activeProgressCallback )
        return;

      nlohmann::json progress = buildProgressFromTurn( turn );
      if( !progress.empty() )
        m_activeProgressCallback( progress );
    }
  );
}


void LlmToolGui::disconnectMcpProgress()
{
  m_mcpResponseAddedConn.disconnect();
  m_activeProgressCallback = McpProgressCallback();
}


void LlmToolGui::invokeAgentForMcp( AgentType agentType, const string &prompt,
                                    McpPromptCallback callback, McpProgressCallback progressCallback,
                                    const string &id )
{
  if( !callback )
    return;

  if( prompt.empty() )
  {
    callback( string( "Empty prompt: nothing to submit." ) );
    return;
  }

  if( !m_llmInterface )
  {
    callback( string( "The LLM assistant is not configured for this session." ) );
    return;
  }

  if( agentType == AgentType::MainAgent || agentType == AgentType::McpServer )
  {
    callback( string( "Cannot directly invoke MainAgent or McpServer agent types." ) );
    return;
  }

  // Checked here too (not just when the main agent calls the invoke_* tool), so we dont spend a main
  //  agent LLM round-trip only to have the invocation refused.
  const std::optional<std::string> prereq_error = m_llmInterface->subAgentPrerequisiteError( agentType );
  if( prereq_error )
  {
    callback( *prereq_error );
    return;
  }

  if( isBenchmarkRunning() )
  {
    callback( string( "A benchmark is currently running; cannot invoke agent right now." ) );
    return;
  }

  // Build a prompt that instructs the assistant to invoke the specific agent tool.
  // This goes through the normal UI conversation flow, so the user sees the agent in the chat.
  const string agentName = agentTypeToString( agentType );
  const string wrappedPrompt = "Please use the invoke_" + agentName
    + " tool to handle the following request:\n\n" + prompt;

  queuePromptForMcp( wrappedPrompt, std::move( callback ), std::move( progressCallback ), {}, id );
}


void LlmToolGui::disconnectSpectrumChangedForBenchmark()
{
  m_spectrumChangedConnection.disconnect();
}


void LlmToolGui::reconnectSpectrumChangedForBenchmark()
{
  // Assigning over a live connection drops its handle without disconnecting it (see resetRoot()),
  //  which would leave a second copy of the slot that disconnectSpectrumChangedForBenchmark() can't reach.
  m_spectrumChangedConnection.disconnect();
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
  return isBenchmarkRunning() && !m_spectrumChangedConnection.isConnected();
}


void LlmToolGui::handleStartBenchmark( const string &xmlPath )
{
  if( m_benchmarkRunner && m_benchmarkRunner->isRunning() )
  {
    SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( "Benchmark Running" );
    dialog->contents()->addNew<WText>( "A benchmark is already running. Please wait for it to finish." );
    dialog->addButton( "OK", WidgetUtils::ButtonRole::Affirm );
    return;
  }

  startBenchmark( xmlPath, BenchmarkRunOptions() );
}


bool LlmToolGui::startBenchmark( const string &xmlPath, const BenchmarkRunOptions &options,
                                 std::function<void(const LlmBenchmarkRunner &)> finished )
{
  if( isBenchmarkRunning() || !m_llmInterface )
    return false;

  // Clean up any previous runner
  m_benchmarkRunner.reset();

  cout << "[LlmBenchmark] Starting benchmark: " << xmlPath << endl;

  m_benchmarkRunner = std::make_unique<LlmBenchmarkRunner>( m_viewer, this, options );
  if( finished )
  {
    LlmBenchmarkRunner * const runner = m_benchmarkRunner.get();
    runner->benchmarkFinished().connect( [runner, finished](){ finished( *runner ); } );
  }

  m_benchmarkRunner->startBenchmark( xmlPath );
  return true;
}//startBenchmark(...)


void LlmToolGui::cancelBenchmark()
{
  if( !isBenchmarkRunning() )
    return;

  m_benchmarkRunner->cancelBenchmark();

  // Stop the in-flight request too, so its reply does not arrive into a cancelled run.
  if( m_llmInterface )
    m_llmInterface->cancelAll();
}//cancelBenchmark()


LlmJsSandboxBridge *LlmToolGui::jsSandboxBridge()
{
  return m_jsSandboxBridge;
}


#endif // USE_LLM_INTERFACE
