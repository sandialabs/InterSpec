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

#include <chrono>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WPanel.h>
#include <Wt/WAnchor.h>
#include <Wt/WServer.h>
#include <Wt/WTextArea.h>
#include <Wt/WResource.h>
#include <Wt/WPushButton.h>
#include <Wt/WGridLayout.h>
#include <Wt/WApplication.h>
#include <Wt/Utils.h>
#include <Wt/WImage.h>
#include <Wt/WMemoryResource.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WJavaScriptPreamble.h>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/PopupDiv.h"
#include "SpecUtils/SpecFile.h"

#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmInteractionDisplay.h"
#include "InterSpec/LlmSubAgentFollowup.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace Wt;
using namespace std;

namespace
{
WT_DECLARE_WT_MEMBER
(LlmCopyTextToClipboard, Wt::JavaScriptFunction, "LlmCopyTextToClipboard",
  function( sender, event, textAreaIdOrText )
{
  // `textAreaIdOrText` is either the id of an element to copy the .value of, or (when no such
  //  element exists) the literal text to copy - the latter lets callers copy a string directly.
  var el = document.getElementById(textAreaIdOrText);
  var text = el ? el.value : textAreaIdOrText;
  if( text === null || text === undefined ) {
    console.warn('LlmCopyTextToClipboard: nothing to copy for:', textAreaIdOrText);
    return;
  }

  if( !navigator.clipboard )
  {
    let tempArea = document.createElement("textarea");
    tempArea.value = text;
    tempArea.style.position = "fixed";
    document.body.appendChild(tempArea);
    tempArea.focus();
    tempArea.select();

    let successful = false;
    try
    {
      successful = document.execCommand('copy');
    }catch( ex )
    {
      console.warn('Failed to copy using execCommand:', ex);
    }

    document.body.removeChild( tempArea );

    if( successful )
    {
      console.log('Copied text to clipboard via appending textarea to DOM');
    }else
    {
      console.warn('Failed to copy to clipboard using textarea');
    }

    return successful;
  }//if( !navigator.clipboard )

  navigator.clipboard.writeText(text).then(function() {
    console.log('Copied text to clipboard using async method.');
  }, function(err) {
    console.warn('Failed copying text to clipboard using async method:', err);
  });
}
);
}//namespace


namespace
{
  /** Format timestamp as HH:MM:SS */
  string formatTime( const chrono::system_clock::time_point &tp )
  {
    const time_t tt = chrono::system_clock::to_time_t( tp );
    const tm *local_tm = localtime( &tt );

    ostringstream oss;
    oss << setfill('0')
        << setw(2) << local_tm->tm_hour << ":"
        << setw(2) << local_tm->tm_min << ":"
        << setw(2) << local_tm->tm_sec;

    return oss.str();
  }//formatTime(...)


  /** Format duration in milliseconds as human-readable string */
  string formatDuration( const chrono::milliseconds &duration )
  {
    const long long ms = duration.count();

    if( ms < 1000 )
      return to_string(ms) + "ms";

    const long long seconds = ms / 1000;
    if( seconds < 60 )
      return to_string(seconds) + "s";

    const long long minutes = seconds / 60;
    const long long sec_remainder = seconds % 60;

    ostringstream oss;
    oss << minutes << "m " << sec_remainder << "s";

    return oss.str();
  }//formatDuration(...)


  /** Truncate string to max length, adding ellipsis if needed */
  string truncateString( const string &str, const size_t maxLen )
  {
    if( str.length() <= maxLen )
      return str;

    return str.substr(0, maxLen - 3) + "...";
  }//truncateString(...)


  /** Get first line of a multi-line string */
  string getFirstLine( const string &str )
  {
    const size_t pos = str.find('\n');
    if( pos == string::npos )
      return str;

    return str.substr(0, pos);
  }//getFirstLine(...)


  /** Format JSON string with indentation for display */
  string formatJson( const nlohmann::json &j )
  {
    try
    {
      return j.dump(2);  // Indent with 2 spaces
    }catch( std::exception & )
    {
      return j.dump();  // Fallback to compact format
    }
  }//formatJson(...)


  /** Format JSON string with indentation for display, from string */
  string formatJsonString( const string &jsonStr )
  {
    try
    {
      const nlohmann::json j = nlohmann::json::parse( jsonStr );
      return formatJson( j );
    }catch( std::exception & )
    {
      return jsonStr;  // Not valid JSON, return as-is
    }
  }//formatJsonString(...)


  /** Whether a turn's rawContent() is a request sent *to* the LLM (true) or a response received
   *from* the LLM (false).  A turn stores exactly one raw JSON blob; its direction is fixed by type:
   InitialRequest/ToolResult hold the request body sent to the API; AutoReply is an automatic
   prompt sent to the LLM; everything else (FinalLlmResponse/ToolCall/Error) holds the API response.
   */
  bool isRequestDirectionTurn( LlmInteractionTurn::Type type )
  {
    switch( type )
    {
      case LlmInteractionTurn::Type::InitialRequest:
      case LlmInteractionTurn::Type::ToolResult:
      case LlmInteractionTurn::Type::AutoReply:
      case LlmInteractionTurn::Type::ConversationSummary:
        return true;

      case LlmInteractionTurn::Type::FinalLlmResponse:
      case LlmInteractionTurn::Type::ToolCall:
      case LlmInteractionTurn::Type::Error:
        return false;
    }//switch( type )

    return false;
  }//isRequestDirectionTurn(...)


  /** Returns rawContent() of the nearest turn adjacent to `turn` in the same interaction, in the
   given direction (`dir` = -1 to search backwards, +1 to search forwards), whose direction matches
   `wantRequestDir`.  Used to pair a response turn with the request that produced it (search
   backwards for a request-direction turn), or a request turn with its response (search forwards for
   a response-direction turn).  Returns empty string if the interaction is gone, `turn` is not found
   (e.g. the AutoReply wrapper, which is not stored in `responses`), or no matching turn exists.
   */
  string adjacentTurnRawJson( const shared_ptr<LlmInteractionTurn> &turn,
                              const int dir, const bool wantRequestDir )
  {
    const shared_ptr<LlmInteraction> convo = turn ? turn->conversation().lock() : nullptr;
    if( !convo )
      return string();

    const vector<shared_ptr<LlmInteractionTurn>> &turns = convo->responses;
    size_t index = turns.size();
    for( size_t i = 0; i < turns.size(); ++i )
    {
      if( turns[i] == turn )
      {
        index = i;
        break;
      }
    }

    if( index >= turns.size() )
      return string();  // `turn` not found in the interaction (e.g. a display-only wrapper)

    for( long long i = static_cast<long long>(index) + dir;
         (i >= 0) && (i < static_cast<long long>(turns.size()));
         i += dir )
    {
      const shared_ptr<LlmInteractionTurn> &other = turns[static_cast<size_t>(i)];
      if( other && (isRequestDirectionTurn( other->type() ) == wantRequestDir) )
        return other->rawContent();
    }

    return string();
  }//adjacentTurnRawJson(...)

}//namespace



//
// LlmInteractionTurnDisplay implementation
//

LlmInteractionTurnDisplay::LlmInteractionTurnDisplay( shared_ptr<LlmInteractionTurn> turn )
  : WPanel(),
    m_turn( turn ),
    m_menuIcon( nullptr ),
    m_menu( nullptr ),
    m_bodyContainer( nullptr ),
    m_bodyCreated( false )
{
  addStyleClass( "LlmInteractionTurnDisplay" );

  // Create empty body container - content will be created on construction for now
  auto bodyContainerOwned = std::make_unique<WContainerWidget>();
  m_bodyContainer = bodyContainerOwned.get();
  setCentralWidget( std::move(bodyContainerOwned) );

  // Note: Lazy loading will be added later
  // For now, we create content immediately in derived constructors
}//LlmInteractionTurnDisplay constructor


LlmInteractionTurnDisplay::~LlmInteractionTurnDisplay()
{
  // Nothing to do for m_menu: it is owned by m_menuIcon (WObject::addChild), so it is destroyed
  //  with the button.  Deleting it here would be a double-free.  (Under Wt3 the menu really was
  //  app-global-owned and did need an explicit delete.)
}//~LlmInteractionTurnDisplay()


shared_ptr<LlmInteractionTurn> LlmInteractionTurnDisplay::turn() const
{
  return m_turn;
}//turn()


void LlmInteractionTurnDisplay::handleExpansion( bool expanded )
{
  // Note: This will be connected to WPanel::expanded() signal later
  // For now, content is created immediately in derived constructors
  if( expanded && !m_bodyCreated )
  {
    createBodyContent();
    m_bodyCreated = true;
  }
}//handleExpansion(...)


void LlmInteractionTurnDisplay::createMenuIcon()
{
  // Create menu button in title bar.  The menu itself is built lazily on first
  // click (showMenu), so turns whose menu is never opened cost nothing - which is
  // every turn during headless/benchmark runs.
  m_menuIcon = titleBarWidget()->addNew<WPushButton>();
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn" );
  m_menuIcon->clicked().preventPropagation();
  m_menuIcon->clicked().connect( this, &LlmInteractionTurnDisplay::showMenu );
}//createMenuIcon()


void LlmInteractionTurnDisplay::showMenu()
{
  if( !m_menu )
  {
    // Owned by m_menuIcon, so it dies with this panel.  Not makePopupMenu(): that also wires
    //  m_menuIcon->clicked() to popup(), which would double up with the connection that got us
    //  here (and we want the items built lazily, on first click).
    m_menu = makeButtonOwnedPopupMenu( m_menuIcon );
    addMenuItems( m_menu );
  }
  m_menu->popup( m_menuIcon, Wt::Orientation::Vertical );
  m_menu->doJavaScript( "Wt.WT.BringAboveDialogs('" + m_menu->id() + "');" );
}//showMenu()


void LlmInteractionTurnDisplay::addMenuItems( PopupDivMenu *menu )
{
  // Each turn stores exactly one raw JSON blob (m_turn->rawContent()) - either the request sent to
  // the LLM or the response received from it, depending on the turn type (see isRequestDirectionTurn).
  // The opposite direction is not stored on this turn, but can be recovered from the adjacent turn:
  // a response turn's request is the nearest preceding request-direction turn's payload, and a
  // request turn's response is the following response-direction turn's payload.
  const bool thisIsRequest = isRequestDirectionTurn( m_turn->type() );

  // Show Request JSON
  PopupDivMenuItem *showRequestJson = menu->addMenuItem( "Show Request JSON" );
  showRequestJson->triggered().connect( std::bind( [this,thisIsRequest]() {
    const string json = thisIsRequest ? m_turn->rawContent()
                                       : adjacentTurnRawJson( m_turn, -1, true );
    if( json.empty() )
      showJsonDialog( "Request JSON", "{\"note\": \"Request JSON not available\"}", false );
    else
      showJsonDialog( "Request JSON", json, true );
  }));

  // Show Response JSON
  PopupDivMenuItem *showResponseJson = menu->addMenuItem( "Show Response JSON" );
  showResponseJson->triggered().connect( std::bind( [this,thisIsRequest]() {
    const string json = thisIsRequest ? adjacentTurnRawJson( m_turn, +1, false )
                                       : m_turn->rawContent();
    if( json.empty() )
      showJsonDialog( "Response JSON", "{\"note\": \"Response JSON not available\"}", false );
    else
      showJsonDialog( "Response JSON", json, true );
  }));
}//addMenuItems()


void LlmInteractionTurnDisplay::showJsonDialog( const WString &title,
                                                const string &jsonStr,
                                                bool allowDownload )
{
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( title );

  // Format and display JSON
  const string formattedJson = formatJsonString( jsonStr );

  WTextArea *jsonArea = dialog->contents()->addNew<WTextArea>();
  jsonArea->setReadOnly( true );
  jsonArea->setText( formattedJson );
  jsonArea->setRows( 20 );
  jsonArea->setColumns( 80 );
  jsonArea->addStyleClass( "LlmJsonDisplay" );

  // Add buttons
  LOAD_JAVASCRIPT(wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmCopyTextToClipboard );

  WPushButton *copyBtn = dialog->addButton( "Copy to Clipboard" );
  const string js = "function(s,e){ Wt.WT.LlmCopyTextToClipboard(s,e,'" + jsonArea->id() + "'); }";
  copyBtn->clicked().connect( js );

  if( allowDownload )
  {
    // Create memory resource for download; parented to the dialog so it gets cleaned up with it.
    auto resource = std::make_shared<WMemoryResource>( "application/json" );
    resource->setData( reinterpret_cast<const unsigned char*>(formattedJson.data()),
                      static_cast<int>(formattedJson.size()) );

    // Generate filename with timestamp from the interaction turn
    const time_t tt = chrono::system_clock::to_time_t( m_turn->timestamp() );
    const tm *local_tm = localtime( &tt );
    char timeBuffer[64];
    strftime( timeBuffer, sizeof(timeBuffer), "%Y%m%d_%H%M%S", local_tm );
    const string filename = string("llm_interaction_") + timeBuffer + ".json";
    resource->suggestFileName( filename );
    resource->setDispositionType( ContentDisposition::Attachment );

    // On macOS/iOS only clicks on real anchors trigger a download; the `window.open(...)` a
    //  WPushButton link emits is silently dropped by the WKWebView.
    // Note: we deliberately do not use `SimpleDialog::addButton(...)` here - it closes the dialog
    //  on click, which would destroy `resource` (its child) before the browser has fetched it.
#if( BUILD_AS_OSX_APP || IOS )
    WLink dlLink( resource );
    dlLink.setTarget( Wt::LinkTarget::NewWindow );
    WAnchor *downloadBtn = dialog->footer()->addNew<WAnchor>( dlLink );
    downloadBtn->setStyleClass( "simple-dialog-btn DownloadLink" );
    downloadBtn->setText( "Download" );
#else
    WPushButton *downloadBtn = dialog->footer()->addNew<WPushButton>( "Download" );
    downloadBtn->setIcon( "InterSpec_resources/images/download_small.svg" );
    { WLink lnk( resource ); lnk.setTarget( Wt::LinkTarget::NewWindow ); downloadBtn->setLink( lnk ); }
#endif
  }

  WPushButton *closeBtn = dialog->addButton( "Close" );
  closeBtn->clicked().connect( dialog, &SimpleDialog::accept );
}//showJsonDialog(...)


// No longer needed - clipboard copy now happens directly in client-side JavaScript
// void LlmInteractionTurnDisplay::copyToClipboard( const string &text )
// {
//   // Load the JavaScript function if not already loaded
//   LOAD_JAVASCRIPT(wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmCopyTextToClipboard );
//
//   // Call the JavaScript function directly with the text
//   const string escapedText = WString(text).jsStringLiteral();
//   const string js = "Wt.WT.LlmCopyTextToClipboard(null, null, " + escapedText + ");";
//   doJavaScript( js );
// }//copyToClipboard(...)


//
// LlmInteractionFinalResponseDisplay implementation
//

LlmInteractionFinalResponseDisplay::LlmInteractionFinalResponseDisplay(
  shared_ptr<LlmInteractionFinalResponse> response )
  : LlmInteractionTurnDisplay( response ),
    m_response( response )
{
  addStyleClass( "LlmFinalResponseDisplay" );

  // Set title
  setTitle( getTitleText() );

  // NOT collapsed by default
  setCollapsible( true );
  setCollapsed( false );

  // Create title bar with menu icon (we have to do this AFTER calling setTitle, or the button will disapear)
  createMenuIcon();
  
  // Create content immediately
  createBodyContent();
  m_bodyCreated = true;
}//LlmInteractionFinalResponseDisplay constructor


LlmInteractionFinalResponseDisplay::~LlmInteractionFinalResponseDisplay()
{
}//~LlmInteractionFinalResponseDisplay()


void LlmInteractionFinalResponseDisplay::addMenuItems( PopupDivMenu *menu )
{
  // Base items (Request/Response JSON) first, then a convenience "Copy Response Text" that puts the
  // human-readable assistant text (not the raw JSON) on the clipboard.
  LlmInteractionTurnDisplay::addMenuItems( menu );

  PopupDivMenuItem *copyText = menu->addMenuItem( "Copy Response Text" );
  copyText->triggered().connect( std::bind( [this]() {
    LOAD_JAVASCRIPT(wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmCopyTextToClipboard );
    // Passing the literal text (not an element id) - LlmCopyTextToClipboard falls back to copying it
    //  directly when no element with that "id" exists.
    const string js = "Wt.WT.LlmCopyTextToClipboard(null,null,"
                      + WString(m_response->content()).jsStringLiteral() + ");";
    wApp->doJavaScript( js );
  }));
}//addMenuItems()


WString LlmInteractionFinalResponseDisplay::getTitleText() const
{
  const string content = m_response->content();
  const string firstLine = getFirstLine( content );
  const string truncated = truncateString( firstLine, 80 );

  // Escape model-derived text: the title WText renders as XHTMLText, and WString::arg() does not
  // HTML-escape, so raw '&'/'<'/'>' or stray entities would break Wt's XHTML parse.
  return WString( "Response: {1}" ).arg( Wt::Utils::htmlEncode( truncated ) );
}//getTitleText()


void LlmInteractionFinalResponseDisplay::createBodyContent()
{
  // Display the full response content.  Pass TextFormat::Plain in the constructor: the 3-arg WText ctor sets
  // the format BEFORE the text, so the raw model content is never run through Wt's XHTML parser (which
  // would log "Error reading XHTML string" on a bare '&' or '<' before silently falling back).
  WText *contentText = m_bodyContainer->addNew<WText>( m_response->content(), TextFormat::Plain );
  contentText->addStyleClass( "LlmResponseContent" );

  // Add metadata section
  WContainerWidget *metadataDiv = m_bodyContainer->addNew<WContainerWidget>();
  metadataDiv->addStyleClass( "LlmMetadata" );

  // Timestamp
  const string timeStr = formatTime( m_response->timestamp() );
  WText *timeText = metadataDiv->addNew<WText>( "Time: " + timeStr );
  timeText->addStyleClass( "LlmTimestamp" );

  // Duration
  if( m_response->callDuration().has_value() )
  {
    const string durationStr = formatDuration( m_response->callDuration().value() );
    WText *durationText = metadataDiv->addNew<WText>( " | Duration: " + durationStr );
    durationText->addStyleClass( "LlmDuration" );
  }

  // Token usage (if available from parent conversation)
  if( const shared_ptr<LlmInteraction> convo = m_response->conversation().lock() )
  {
    if( convo->promptTokens.has_value() )
    {
      string tokenStr = " | Tokens: " + to_string(convo->promptTokens.value())
                      + " in + " + to_string(convo->completionTokens.value_or(0))
                      + " out = " + to_string(convo->totalTokens.value_or(0));
      if( convo->cachedTokens.has_value() )
        tokenStr += " (" + to_string(convo->cachedTokens.value()) + " cached)";
      WText *tokenText = metadataDiv->addNew<WText>( tokenStr );
      tokenText->addStyleClass( "LlmTokenUsage" );
    }
  }

  // Thinking content (if available)
  if( !m_response->thinkingContent().empty() )
  {
    WContainerWidget *thinkingDiv = m_bodyContainer->addNew<WContainerWidget>();
    thinkingDiv->addStyleClass( "LlmThinkingContent" );

    WText *thinkingLabel = thinkingDiv->addNew<WText>( "Thinking:" );
    thinkingLabel->addStyleClass( "LlmThinkingLabel" );

    // TextFormat::Plain via the 3-arg ctor (format before text) so raw thinking content skips the XHTML parser.
    WText *thinkingText = thinkingDiv->addNew<WText>( m_response->thinkingContent(), TextFormat::Plain );
  }
}//createBodyContent()


//
// LlmToolRequestDisplay implementation
//

LlmToolRequestDisplay::LlmToolRequestDisplay( shared_ptr<LlmToolRequest> request )
  : LlmInteractionTurnDisplay( request ),
    m_request( request )
{
  addStyleClass( "LlmToolRequestDisplay" );

  // Set title
  setTitle( getTitleText() );

  // Collapsed by default
  setCollapsible( true );
  setCollapsed( true );

  // Create title bar with menu icon (we have to do this AFTER calling setTitle, or the button will disapear)
  createMenuIcon();
  
  // Create content immediately (lazy loading to be added later)
  createBodyContent();
  m_bodyCreated = true;
}//LlmToolRequestDisplay constructor


LlmToolRequestDisplay::~LlmToolRequestDisplay()
{
}//~LlmToolRequestDisplay()


WString LlmToolRequestDisplay::getTitleText() const
{
  const vector<LlmToolCall> &calls = m_request->toolCalls();

  if( calls.empty() )
    return "Tool call requests: (none)";

  // Build comma-separated list of tool names
  ostringstream oss;
  oss << "Tool call requests: ";
  for( size_t i = 0; i < calls.size(); ++i )
  {
    if( i > 0 )
      oss << ", ";
    oss << calls[i].toolName;
  }

  // Escape model-derived tool names before they reach the XHTMLText panel title.
  const string truncated = truncateString( oss.str(), 100 );
  return WString( Wt::Utils::htmlEncode( truncated ) );
}//getTitleText()


void LlmToolRequestDisplay::createBodyContent()
{
  const vector<LlmToolCall> &calls = m_request->toolCalls();

  if( calls.empty() )
  {
    WText *emptyText = m_bodyContainer->addNew<WText>( "No tool calls" );
    emptyText->addStyleClass( "LlmEmptyMessage" );
    return;
  }

  // Display each tool call
  for( const LlmToolCall &call : calls )
  {
    WContainerWidget *callDiv = m_bodyContainer->addNew<WContainerWidget>();
    callDiv->addStyleClass( "LlmToolCallItem" );

    // Tool name - escape the model-derived name but keep the intentional <b> markup.
    WText *nameText = callDiv->addNew<WText>( "<b>" + Wt::Utils::htmlEncode( call.toolName ) + "</b>" );
    nameText->setTextFormat( TextFormat::UnsafeXHTML );
    nameText->addStyleClass( "LlmToolName" );

    // Invocation ID - TextFormat::Plain so a stray '&'/'<' in an id can't hit the XHTML parser.
    WText *idText = callDiv->addNew<WText>( " (ID: " + call.invocationId + ")", TextFormat::Plain );
    idText->addStyleClass( "LlmInvocationId" );

    // Parameters
    if( !call.toolParameters.empty() )
    {
      WText *paramsLabel = callDiv->addNew<WText>( "Parameters:" );
      paramsLabel->addStyleClass( "LlmParamsLabel" );

      WTextArea *paramsArea = callDiv->addNew<WTextArea>();
      paramsArea->setReadOnly( true );
      paramsArea->setText( formatJson(call.toolParameters) );
      paramsArea->setRows( 5 );
      paramsArea->addStyleClass( "LlmToolParams" );
    }
  }
}//createBodyContent()


//
// LlmToolResultsDisplay implementation
//

LlmToolResultsDisplay::LlmToolResultsDisplay( shared_ptr<LlmToolResults> results,
                                              std::weak_ptr<LlmInterface> llmInterface,
                                              int nestingLevel )
  : LlmInteractionTurnDisplay( results ),
    m_results( results ),
    m_llmInterface( llmInterface ),
    m_nestingLevel( nestingLevel )
{
  addStyleClass( "LlmToolResultsDisplay" );

  // Set title
  setTitle( getTitleText() );

  // Collapsed by default, but expand if any tool call has visual content
  setCollapsible( true );

  bool hasVisuals = false;
  for( const LlmToolCall &call : m_results->toolCalls() )
  {
    if( !call.visualContent.empty() || call.imageContent.has_value() )
    {
      hasVisuals = true;
      break;
    }
  }

  setCollapsed( !hasVisuals );

  // Create title bar with menu icon (we have to do this AFTER calling setTitle, or the button will disapear)
  createMenuIcon();
  
  // Create content immediately (lazy loading to be added later)
  createBodyContent();
  m_bodyCreated = true;
}//LlmToolResultsDisplay constructor


LlmToolResultsDisplay::~LlmToolResultsDisplay()
{
}//~LlmToolResultsDisplay()


WString LlmToolResultsDisplay::getTitleText() const
{
  const vector<LlmToolCall> &calls = m_results->toolCalls();

  if( calls.empty() )
    return "Tool results: (none)";

  // Check if any errors occurred
  bool hasErrors = false;
  for( const LlmToolCall &call : calls )
  {
    if( call.status == LlmToolCall::CallStatus::Error )
    {
      hasErrors = true;
      break;
    }
  }

  // Build comma-separated list of tool names
  ostringstream oss;
  oss << "Tool results";
  if( hasErrors )
    oss << " [HAS ERRORS]";
  oss << ": ";
  for( size_t i = 0; i < calls.size(); ++i )
  {
    if( i > 0 )
      oss << ", ";
    oss << calls[i].toolName;
  }

  // Escape model-derived tool names before they reach the XHTMLText panel title.
  const string truncated = truncateString( oss.str(), 100 );
  return WString( Wt::Utils::htmlEncode( truncated ) );
}//getTitleText()


void LlmToolResultsDisplay::renderToolCallResult( const LlmToolCall &call,
                                                    WContainerWidget *resultDiv )
{
  // Error styling
  resultDiv->removeStyleClass( "LlmToolResultError" );
  if( call.status == LlmToolCall::CallStatus::Error )
    resultDiv->addStyleClass( "LlmToolResultError" );

  // Tool name with status indicator - escape the model-derived name but keep the intentional markup.
  ostringstream nameHtml;
  nameHtml << "<b>" << Wt::Utils::htmlEncode( call.toolName ) << "</b>";
  if( call.status == LlmToolCall::CallStatus::Error )
    nameHtml << " <span class='LlmToolStatusError'>[ERROR]</span>";
  else if( call.status == LlmToolCall::CallStatus::Pending && !call.sub_agent_conversation )
    nameHtml << " <span class='LlmToolStatusPending'>[PENDING]</span>";

  WText *nameText = resultDiv->addNew<WText>( nameHtml.str() );
  nameText->setTextFormat( TextFormat::UnsafeXHTML );
  nameText->addStyleClass( "LlmToolName" );

  // Invocation ID - TextFormat::Plain so a stray '&'/'<' in an id can't hit the XHTML parser.
  WText *idText = resultDiv->addNew<WText>( " (ID: " + call.invocationId + ")", TextFormat::Plain );
  idText->addStyleClass( "LlmInvocationId" );

  // Execution duration
  if( call.executionDuration.has_value() )
  {
    const string durationStr = formatDuration( call.executionDuration.value() );
    WText *durationText = resultDiv->addNew<WText>( " - took " + durationStr );
    durationText->addStyleClass( "LlmDuration" );
  }

  // Result content - show if not pending, available, and (no sub-agent OR sub-agent completed)
  if( !call.content.empty()
      && call.status != LlmToolCall::CallStatus::Pending
      && (!call.sub_agent_conversation || call.sub_agent_conversation->finishTime.has_value()) )
  {
    WText *resultLabel = resultDiv->addNew<WText>( "Result:" );
    resultLabel->addStyleClass( "LlmResultLabel" );

    WTextArea *resultArea = resultDiv->addNew<WTextArea>();
    resultArea->setReadOnly( true );
    resultArea->setText( call.content );
    resultArea->setRows( 5 );
    resultArea->addStyleClass( "LlmToolResult" );
  }

  // Show image preview if the tool result included an image sent to the LLM
  if( call.imageContent.has_value() )
  {
    const LlmToolCall::ImageContent &img = call.imageContent.value();
    const string decoded = Wt::Utils::base64Decode( img.base64Data );
    vector<unsigned char> bytes( decoded.begin(), decoded.end() );
    auto imgResource = std::make_shared<WMemoryResource>( img.mimeType );
    imgResource->setData( bytes );

    WImage *thumb = resultDiv->addNew<WImage>( WLink( imgResource ) );
    thumb->addStyleClass( "LlmToolResultImageThumb" );
  }

  // Render visual content (display-only, not sent to LLM)
  for( const LlmToolCall::VisualContent &visual : call.visualContent )
    renderVisualContent( visual, resultDiv );
}//renderToolCallResult(...)


void LlmToolResultsDisplay::renderVisualContent( const LlmToolCall::VisualContent &visual,
                                                   WContainerWidget *parentDiv )
{
  WContainerWidget *vizDiv = parentDiv->addNew<WContainerWidget>();
  vizDiv->addStyleClass( "LlmToolVisual" );

  if( !visual.title.empty() )
  {
    // TextFormat::Plain via the 3-arg ctor (format before text) keeps model-derived text out of the XHTML parser.
    WText *titleText = vizDiv->addNew<WText>( WString::fromUTF8( visual.title ), Wt::TextFormat::Plain );
    titleText->addStyleClass( "LlmToolVisualTitle" );
  }

  if( !visual.description.empty() )
  {
    WText *descText = vizDiv->addNew<WText>( WString::fromUTF8( visual.description ), Wt::TextFormat::Plain );
    descText->addStyleClass( "LlmToolVisualDesc" );
  }

  WContainerWidget *chartDiv = vizDiv->addNew<WContainerWidget>();
  chartDiv->addStyleClass( "LlmToolVisualChart" );

  const int defaultHeight = [&visual]() -> int {
    switch( visual.type )
    {
      case LlmToolCall::VisualContent::VisualType::ChiSquaredPlot:     return 250;
      case LlmToolCall::VisualContent::VisualType::ShieldingDiagram2D: return 300;
      case LlmToolCall::VisualContent::VisualType::SpectrumChart:      return 250;
      case LlmToolCall::VisualContent::VisualType::Svg:                return 200;
      case LlmToolCall::VisualContent::VisualType::Html:               return 400;
    }
    return 250;
  }();

  const int height = (visual.suggestedHeightPx > 0) ? visual.suggestedHeightPx : defaultHeight;
  chartDiv->setHeight( WLength( height ) );

  // Width: use suggested if provided, otherwise CSS class default (500px, max-width:75vw).
  if( visual.suggestedWidthPx > 0 )
    chartDiv->setWidth( WLength( visual.suggestedWidthPx ) );

  switch( visual.type )
  {
    case LlmToolCall::VisualContent::VisualType::ChiSquaredPlot:
    {
      wApp->require( "InterSpec_resources/ShieldingSourceFitPlot.js" );
      const string chartId = chartDiv->id();

      // Create the chart and set data, then attach ResizeObserver
      const string js =
        "var el=document.getElementById('" + chartId + "');"
        "if(el){"
          "var c=new ShieldingSourceFitPlot(el, {margins:{top:5,right:10,bottom:30,left:45}});"
          "el._chartObj=c;"
          "c.setData(" + visual.data + ");"
          "new ResizeObserver(function(){"
            "if(el._chartObj) el._chartObj.handleResize();"
          "}).observe(el);"
        "}";
      chartDiv->doJavaScript( js );
      break;
    }

    case LlmToolCall::VisualContent::VisualType::ShieldingDiagram2D:
    {
      wApp->require( "InterSpec_resources/Shielding2DView.js" );
      const string chartId = chartDiv->id();
      chartDiv->doJavaScript(
        "new Shielding2DView('" + chartId + "', " + visual.data + ");" );
      break;
    }

    case LlmToolCall::VisualContent::VisualType::SpectrumChart:
    {
      // Button-only - remove the chart div and put button directly in vizDiv
      vizDiv->removeWidget( chartDiv );
      chartDiv = nullptr;

      const string buttonLabel = visual.title.empty()
        ? "View Spectrum" : ("View " + visual.title);

      WPushButton *viewBtn = vizDiv->addNew<WPushButton>( WString::fromUTF8( buttonLabel ) );
      viewBtn->addStyleClass( "LlmToolVisualReportBtn" );

      // Capture by value - measurement shared_ptr is lightweight (shares channel data)
      const shared_ptr<const SpecUtils::Measurement> specMeas = visual.spectrumMeasurement;
      const string peakJson = visual.peakJson;
      const double xMin = visual.xRangeMin;
      const double xMax = visual.xRangeMax;
      const string dlgTitle = visual.title.empty() ? "Spectrum" : visual.title;

      viewBtn->clicked().connect( std::bind( [specMeas, peakJson, xMin, xMax, dlgTitle](){
        if( !specMeas )
          return;

        AuxWindow *window = AuxWindow::make( WString::fromUTF8( dlgTitle ),
          AuxWindowProperties::EnableResize | AuxWindowProperties::DisableCollapse );
        window->rejectWhenEscapePressed();
        window->setClosable( true );

        InterSpec *viewer = InterSpec::instance();
        const int dlgWidth = viewer ? static_cast<int>( 0.8 * viewer->renderedWidth() ) : 700;
        const int dlgHeight = viewer ? static_cast<int>( 0.8 * viewer->renderedHeight() ) : 500;
        window->resize( WLength( dlgWidth ), WLength( dlgHeight ) );

        WContainerWidget *contents = window->contents();
        contents->setOverflow( Overflow::Hidden );

        D3SpectrumDisplayDiv *chart = contents->addNew<D3SpectrumDisplayDiv>();
        chart->setCompactAxis( true );
        chart->setWidth( WLength( 100, WLength::Unit::Percentage ) );
        chart->setHeight( WLength( 100, WLength::Unit::Percentage ) );
        chart->setData( specMeas, true );
        chart->setRoiData( peakJson, SpecUtils::SpectrumType::Foreground );
        if( xMax > xMin )
          chart->setXAxisRange( xMin, xMax );

        WPushButton *closeBtn = window->addCloseButtonToFooter();
        closeBtn->clicked().connect( window, [window](){ AuxWindow::deleteAuxWindow( window ); } );

        window->centerWindow();
        window->show();

        window->finished().connect( window, [window](){ AuxWindow::deleteAuxWindow( window ); } );
      } ) );

      break;
    }

    case LlmToolCall::VisualContent::VisualType::Svg:
    {
      WText *svgWidget = chartDiv->addNew<WText>( WString::fromUTF8( visual.data ) );
      svgWidget->setTextFormat( TextFormat::UnsafeXHTML );
      break;
    }

    case LlmToolCall::VisualContent::VisualType::Html:
    {
      // Button-only - remove the chart div and put button directly in vizDiv
      vizDiv->removeWidget( chartDiv );
      chartDiv = nullptr;

      const string buttonLabel = visual.title.empty()
        ? "View Report" : ("View " + visual.title);

      WPushButton *viewBtn = vizDiv->addNew<WPushButton>( WString::fromUTF8( buttonLabel ) );
      viewBtn->addStyleClass( "LlmToolVisualReportBtn" );

      // Capture the HTML data by value for the click handler
      const string htmlData = visual.data;
      const string dlgTitle = visual.title.empty() ? "Report" : visual.title;

      viewBtn->clicked().connect( std::bind( [htmlData, dlgTitle](){
        AuxWindow *window = AuxWindow::make( WString::fromUTF8( dlgTitle ),
          AuxWindowProperties::EnableResize | AuxWindowProperties::DisableCollapse );
        window->rejectWhenEscapePressed();
        window->setClosable( true );

        // Size to 80% of InterSpec rendered area
        InterSpec *viewer = InterSpec::instance();
        const int dlgWidth = viewer ? static_cast<int>( 0.8 * viewer->renderedWidth() ) : 700;
        const int dlgHeight = viewer ? static_cast<int>( 0.8 * viewer->renderedHeight() ) : 500;
        window->resize( WLength( dlgWidth ), WLength( dlgHeight ) );

        WContainerWidget *contents = window->contents();
        contents->setOverflow( Overflow::Hidden );

        // Create iframe via JS after the dialog is rendered
        WContainerWidget *iframeDiv = contents->addNew<WContainerWidget>();
        iframeDiv->setWidth( WLength( 100, WLength::Unit::Percentage ) );
        iframeDiv->setHeight( WLength( 100, WLength::Unit::Percentage ) );

        // Escape the HTML for embedding in JS string literal
        string escapedHtml = htmlData;
        for( size_t pos = 0; (pos = escapedHtml.find( '\\', pos )) != string::npos; pos += 2 )
          escapedHtml.insert( pos, "\\" );
        for( size_t pos = 0; (pos = escapedHtml.find( '\'', pos )) != string::npos; pos += 2 )
          escapedHtml.insert( pos, "\\" );
        for( size_t pos = 0; (pos = escapedHtml.find( '\n', pos )) != string::npos; )
          escapedHtml.replace( pos, 1, "\\n" );
        for( size_t pos = 0; (pos = escapedHtml.find( '\r', pos )) != string::npos; )
          escapedHtml.replace( pos, 1, "\\r" );

        const string iframeId = iframeDiv->id();
        const string js =
          "var el=document.getElementById('" + iframeId + "');"
          "if(el){"
            "var iframe=document.createElement('iframe');"
            "iframe.style.width='100%';"
            "iframe.style.height='100%';"
            "iframe.style.border='none';"
            "iframe.setAttribute('sandbox','allow-scripts');"
            "iframe.srcdoc='" + escapedHtml + "';"
            "el.appendChild(iframe);"
          "}";
        iframeDiv->doJavaScript( js );

        // Add close button to the footer
        WPushButton *closeBtn = window->addCloseButtonToFooter();
        closeBtn->clicked().connect( window, [window](){ AuxWindow::deleteAuxWindow( window ); } );

        window->centerWindow();
        window->show();

        window->finished().connect( window, [window](){ AuxWindow::deleteAuxWindow( window ); } );
      } ) );

      break;
    }
  }//switch( visual.type )
}//renderVisualContent(...)


void LlmToolResultsDisplay::createBodyContent()
{
  const vector<LlmToolCall> &calls = m_results->toolCalls();

  if( calls.empty() )
  {
    WText *emptyText = m_bodyContainer->addNew<WText>( "No tool results" );
    emptyText->addStyleClass( "LlmEmptyMessage" );
    return;
  }

  // Display each tool result
  for( const LlmToolCall &call : calls )
  {
    WContainerWidget *resultDiv = m_bodyContainer->addNew<WContainerWidget>();
    resultDiv->addStyleClass( "LlmToolResultItem" );

    // Render the common parts (name, status, duration, result content)
    renderToolCallResult( call, resultDiv );

    // Sub-agent conversation (if present)
    if( call.sub_agent_conversation )
    {
        WText *subAgentLabel = resultDiv->addNew<WText>( "Sub-agent conversation:" );
        subAgentLabel->addStyleClass( "LlmSubAgentLabel" );

      /*
        // Display input arguments (context and task)
        if( !call.toolParameters.is_null() )
        {
          // Extract context field
          if( call.toolParameters.contains("context") && call.toolParameters["context"].is_string() )
          {
            const string context = call.toolParameters["context"].get<string>();
            if( !context.empty() )
            {
              WText *contextLabel = resultDiv->addNew<WText>( "Context:" );
              contextLabel->addStyleClass( "LlmParamsLabel" );

              WTextArea *contextArea = resultDiv->addNew<WTextArea>();
              contextArea->setReadOnly( true );
              contextArea->setText( context );
              contextArea->setRows( 3 );
              contextArea->addStyleClass( "LlmToolParams" );
            }
          }

          // Extract task field
          if( call.toolParameters.contains("task") && call.toolParameters["task"].is_string() )
          {
            const string task = call.toolParameters["task"].get<string>();
            if( !task.empty() )
            {
              WText *taskLabel = resultDiv->addNew<WText>( "Task:" );
              taskLabel->addStyleClass( "LlmParamsLabel" );

              WTextArea *taskArea = resultDiv->addNew<WTextArea>();
              taskArea->setReadOnly( true );
              taskArea->setText( task );
              taskArea->setRows( 3 );
              taskArea->addStyleClass( "LlmToolParams" );
            }
          }
        }
       */

        // Create nested LlmInteractionDisplay for sub-agent
        LlmInteractionDisplay *subAgentDisplay = resultDiv->addNew<LlmInteractionDisplay>(
                    call.sub_agent_conversation, m_llmInterface, m_nestingLevel + 1 );
        subAgentDisplay->addStyleClass( "LlmSubAgentNesting" );

        m_subAgentDisplays.push_back( subAgentDisplay );

        // Store the result div for later when sub-agent finishes
        m_subAgentResultDivs[call.sub_agent_conversation] = resultDiv;

        // Connect to sub-agent finished signal to update results
        {
          const auto turn = m_turn;
          const auto conv = call.sub_agent_conversation;
          call.sub_agent_conversation->conversationFinished.connect(
            this, [this,turn,conv](){ handleSubAgentFinished( turn, conv ); } );
        }

      // Make sure we expand the parent widget so we can see the sub-agent conversation going on.
      if( isCollapsed() )
        setCollapsed( false );
    }//if( call.sub_agent_conversation )

    // For pending async tools (no sub-agent), store div for later update
    if( call.status == LlmToolCall::CallStatus::Pending && !call.sub_agent_conversation )
    {
      m_asyncToolResultDivs[call.invocationId] = resultDiv;
    }
  }

  // Connect to async tool completion signal if there are any pending async tools
  if( !m_asyncToolResultDivs.empty() )
  {
    m_results->asyncToolCompleted().connect( this, [this]( const std::string &a1 ){
      handleAsyncToolCompleted( a1 );
    } );
  }
}//createBodyContent()


void LlmToolResultsDisplay::addMenuItems( PopupDivMenu *menu )
{
  // Base items first, then the tool-results-specific item
  LlmInteractionTurnDisplay::addMenuItems( menu );

  // Add "Show Tool Results JSON" menu item
  PopupDivMenuItem *showToolResultsJson = menu->addMenuItem( "Show Tool Results JSON" );
  showToolResultsJson->triggered().connect( std::bind( [this]() {
    // Build JSON for all tool results
    const vector<LlmToolCall> &calls = m_results->toolCalls();

    nlohmann::json resultsJson = nlohmann::json::array();
    for( const LlmToolCall &call : calls )
    {
      nlohmann::json resultObj = nlohmann::json::object();

      // Add status
      string statusStr;
      switch( call.status )
      {
        case LlmToolCall::CallStatus::Pending: statusStr = "pending"; break;
        case LlmToolCall::CallStatus::Success: statusStr = "success"; break;
        case LlmToolCall::CallStatus::Error: statusStr = "error"; break;
      }
      resultObj["status"] = statusStr;

      resultObj["tool_name"] = call.toolName;
      resultObj["invocation_id"] = call.invocationId;
      resultObj["content"] = call.content;

      if( !call.toolParameters.is_null() )
      {
        resultObj["parameters"] = call.toolParameters;
      }

      if( call.executionDuration.has_value() )
      {
        resultObj["execution_duration_ms"] = call.executionDuration.value().count();
      }

      if( call.sub_agent_conversation )
      {
        resultObj["has_sub_agent"] = true;
      }

      resultsJson.push_back( resultObj );
    }

    const string jsonStr = resultsJson.dump( 2 );  // Pretty print with 2-space indent
    showJsonDialog( "Tool Results JSON", jsonStr, true );
  }));
}//addMenuItems()


void LlmToolResultsDisplay::handleSubAgentFinished( shared_ptr<LlmInteractionTurn> turn,
                                                    shared_ptr<LlmInteraction> subAgent )
{
  // Sub-agent finished - create the result text area if content is available
  // Find the corresponding tool call in m_results
  for( const LlmToolCall &call : m_results->toolCalls() )
  {
    if( call.sub_agent_conversation == subAgent && !call.content.empty() )
    {
      // Find the result div for this sub-agent
      const auto divIter = m_subAgentResultDivs.find( subAgent );
      if( divIter == m_subAgentResultDivs.end() )
      {
        cerr << "LlmToolResultsDisplay::handleSubAgentFinished(): Could not find result div for sub-agent" << endl;
        return;
      }

      WContainerWidget *resultDiv = divIter->second;

      // Create the result label and text area
      WText *resultLabel = resultDiv->addNew<WText>( "Result:" );
      resultLabel->addStyleClass( "LlmResultLabel" );

      WTextArea *resultArea = resultDiv->addNew<WTextArea>();
      resultArea->setReadOnly( true );
      resultArea->setText( call.content );
      resultArea->setRows( 5 );
      resultArea->addStyleClass( "LlmToolResult" );

      break;
    }
  }
}//handleSubAgentFinished(...)


void LlmToolResultsDisplay::handleAsyncToolCompleted( const string &invocationId )
{
  const auto divIter = m_asyncToolResultDivs.find( invocationId );
  if( divIter == m_asyncToolResultDivs.end() )
    return;

  // Find the updated tool call data
  for( const LlmToolCall &call : m_results->toolCalls() )
  {
    if( call.invocationId != invocationId )
      continue;

    WContainerWidget *resultDiv = divIter->second;

    // Clear and re-render using the shared helper
    resultDiv->clear();
    renderToolCallResult( call, resultDiv );

    // Expand the panel if the completed tool call has visual content
    if( !call.visualContent.empty() && isCollapsed() )
      setCollapsed( false );

    // Remove from tracking map since it's now complete
    m_asyncToolResultDivs.erase( divIter );
    break;
  }
}//handleAsyncToolCompleted(...)


//
// LlmInteractionErrorDisplay implementation
//

LlmInteractionErrorDisplay::LlmInteractionErrorDisplay( shared_ptr<LlmInteractionError> error,
                                                        std::weak_ptr<LlmInterface> llmInterface )
  : LlmInteractionTurnDisplay( error ),
    m_error( error ),
    m_llmInterface( llmInterface ),
    m_retryBtn( nullptr ),
    m_continueBtn( nullptr ),
    m_helpText( nullptr )
{
  addStyleClass( "LlmErrorDisplay" );

  // Set title
  setTitle( getTitleText() );

  // NOT collapsed by default
  setCollapsible( true );
  setCollapsed( false );

  // Create title bar with menu icon (we have to do this AFTER calling setTitle, or the button will disapear)
  createMenuIcon();
  
  // Create content immediately
  createBodyContent();
  m_bodyCreated = true;
}//LlmInteractionErrorDisplay constructor


LlmInteractionErrorDisplay::~LlmInteractionErrorDisplay()
{
}//~LlmInteractionErrorDisplay()



WString LlmInteractionErrorDisplay::getTitleText() const
{
  const string errorMsg = m_error->errorMessage();
  const string firstLine = getFirstLine( errorMsg );
  const string truncated = truncateString( firstLine, 80 );

  // Add error type indicator to title
  string prefix;
  switch( m_error->errorType() )
  {
    case LlmInteractionError::ErrorType::Timeout:
      prefix = "Timeout Error";
      break;
    case LlmInteractionError::ErrorType::Network:
      prefix = "Network Error";
      break;
    case LlmInteractionError::ErrorType::LlmApi:
      prefix = "LLM API Error";
      break;
    case LlmInteractionError::ErrorType::JsonParse:
      prefix = "Parse Error";
      break;
    default:
      prefix = "Error";
  }

  // Escape the model-derived message portion; prefix is a fixed literal.
  return WString( "{1}: {2}" ).arg( prefix ).arg( Wt::Utils::htmlEncode( truncated ) );
}//getTitleText()


void LlmInteractionErrorDisplay::createBodyContent()
{
  // Display error message
  WText *errorText = m_bodyContainer->addNew<WText>( m_error->errorMessage(), Wt::TextFormat::Plain );
  errorText->addStyleClass( "LlmErrorMessage" );
  cout << "=======Setting error body text to:\n" << m_error->errorMessage() << "=======" << endl;

  // Add retry status if automatic retry was attempted
  if( m_error->supportsAutomaticRetry() && m_error->retryAttempted() )
  {
    WText *retryInfo = m_bodyContainer->addNew<WText>( "Automatic retry was attempted" );
    retryInfo->addStyleClass( "LlmErrorRetryInfo" );
  }

  // Add retry button (always shown for all errors)
  m_retryBtn = m_bodyContainer->addNew<WPushButton>( "Retry" );
  m_retryBtn->addStyleClass( "LlmRetryButton" );
  m_retryBtn->clicked().connect( this, &LlmInteractionErrorDisplay::handleRetry );

  // Add "Continue Anyway" button (always shown for all errors)
  m_continueBtn = m_bodyContainer->addNew<WPushButton>( "Continue Anyway" );
  m_continueBtn->addStyleClass( "LlmContinueButton" );
  m_continueBtn->clicked().connect( this, &LlmInteractionErrorDisplay::handleContinueAnyway );

  // Add help text
  m_helpText = m_bodyContainer->addNew<WText>( "Click Retry to attempt this request again, or Continue Anyway to start a new conversation." );
  m_helpText->addStyleClass( "LlmErrorHelpText" );
}//createBodyContent()


void LlmInteractionErrorDisplay::handleRetry()
{
  // Disable both buttons to prevent double-clicks
  if( m_retryBtn )
    m_retryBtn->setEnabled( false );
  if( m_continueBtn )
  {
    m_continueBtn->setEnabled( false );
    m_continueBtn->hide();  // Hide the button that wasn't clicked
  }

  // Update help text to indicate retry action was taken
  if( m_helpText )
    m_helpText->setText( "Retrying request..." );

  // Lock the weak_ptr to get access to LlmInterface
  std::shared_ptr<LlmInterface> llmInterface = m_llmInterface.lock();
  if( !llmInterface )
  {
    cerr << "Cannot retry - LlmInterface no longer available" << endl;
    if( m_helpText )
      m_helpText->setText( "Error: LLM interface no longer available" );
    return;
  }

  // Get the conversation
  std::shared_ptr<LlmInteraction> convo = m_error->conversation().lock();
  if( !convo )
  {
    cerr << "Cannot retry - conversation no longer available" << endl;
    if( m_helpText )
      m_helpText->setText( "Error: Conversation no longer available" );
    return;
  }

  // Validate the conversation has this error as the last response
  if( convo->responses.empty() )
  {
    cerr << "Cannot retry - conversation has no responses" << endl;
    if( m_helpText )
      m_helpText->setText( "Error: No responses to retry" );
    return;
  }

  std::shared_ptr<LlmInteractionTurn> lastResponse = convo->responses.back();
  if( lastResponse != m_error )
  {
    cerr << "Cannot retry - error is not the last response" << endl;
    if( m_helpText )
      m_helpText->setText( "Error: Can only retry the most recent error" );
    return;
  }

  cout << "\n=== Retrying request for conversation: " << convo->conversationId << " ===" << endl;
  cout << "Error type: " << static_cast<int>(m_error->errorType()) << endl;

  try
  {
    // Remove the error from conversation history
    // This allows the conversation to continue as if the error never happened
    convo->responses.pop_back();

    // Build API request with the conversation history (which now excludes the error)
    const nlohmann::json requestJson = llmInterface->buildMessagesArray( convo );

    // Make tracked API call
    std::pair<int,std::string> request_id_content =
      llmInterface->makeTrackedApiCall( requestJson, convo );

    cout << "Retry request sent with ID: " << request_id_content.first << endl;

    // The response will be handled by the normal LlmInterface::handleApiResponse flow
    // which will add new response turns to the conversation and emit signals
  }
  catch( const std::exception &e )
  {
    cerr << "Error during retry: " << e.what() << endl;

    // Re-add the error to the conversation since retry failed
    convo->responses.push_back( m_error );

    // Update help text
    if( m_helpText )
      m_helpText->setText( std::string("Retry failed: ") + e.what() );

    // Re-enable the retry button so user can try again
    if( m_retryBtn )
      m_retryBtn->setEnabled( true );
    if( m_continueBtn )
    {
      m_continueBtn->setEnabled( true );
      m_continueBtn->show();
    }
  }
}//handleRetry()


void LlmInteractionErrorDisplay::handleContinueAnyway()
{
  // Disable both buttons to prevent double-clicks
  if( m_continueBtn )
    m_continueBtn->setEnabled( false );
  if( m_retryBtn )
  {
    m_retryBtn->setEnabled( false );
    m_retryBtn->hide();  // Hide the button that wasn't clicked
  }

  // Update help text to indicate continue action was taken
  if( m_helpText )
    m_helpText->setText( "Continuing with new conversation. The error has been excluded from history." );

  // Mark the error to exclude from history sent to LLM
  m_error->setExcludeFromHistory( true );

  // Find and mark the request that led to this error to also exclude from history
  // Search backwards from the error to find the closest request we sent to the LLM:
  // - InitialRequest (user/system message)
  // - AutoReply (automatic continuation)
  // - ToolResults (results from tool execution)
  std::shared_ptr<LlmInteraction> convo = m_error->conversation().lock();
  if( convo && !convo->responses.empty() )
  {
    // Find the error in the responses array
    auto errorIt = std::find( convo->responses.begin(), convo->responses.end(), m_error );
    if( errorIt != convo->responses.end() )
    {
      // Search backwards from the error to find the request that triggered it
      auto it = errorIt;
      while( it != convo->responses.begin() )
      {
        --it;
        const LlmInteractionTurn::Type turnType = (*it)->type();

        // Check if this is a request type we send to the LLM
        if( turnType == LlmInteractionTurn::Type::InitialRequest ||
            turnType == LlmInteractionTurn::Type::AutoReply ||
            turnType == LlmInteractionTurn::Type::ToolResult )
        {
          (*it)->setExcludeFromHistory( true );
          cout << "Marked error and preceding request (type: " << static_cast<int>(turnType)
               << ") to exclude from history" << endl;

          // A ToolResult answers a preceding assistant ToolCall (tool_use block).  Excluding the
          // results while leaving the request would create a dangling tool_use, which the Anthropic
          // API rejects - so also exclude the nearest preceding ToolCall turn.
          if( turnType == LlmInteractionTurn::Type::ToolResult )
          {
            for( auto callIt = it; callIt != convo->responses.begin(); )
            {
              --callIt;
              if( (*callIt)->type() == LlmInteractionTurn::Type::ToolCall )
              {
                (*callIt)->setExcludeFromHistory( true );
                break;
              }
            }
          }//if( excluded turn was a ToolResult )
          break;
        }
      }
    }

    // Mark conversation as finished
    convo->finishTime = std::chrono::system_clock::now();
    if( convo->conversation_completion_handler )
      convo->conversation_completion_handler();
    convo->conversationFinished.emit();
  }

  // Call the continue callback to re-enable input
  //blah blah blah

  // The conversation is already marked as finished in handleContinueAnyway()
  // Trigger the status update and call parent's continue callback
  // blah blah blah
  //parent->updateStatus();
}//handleContinueAnyway()


//
// LlmInteractionDisplay implementation
//

int LlmInteractionDisplay::s_nextTimerId = 0;

LlmInteractionDisplay::LlmInteractionDisplay( shared_ptr<LlmInteraction> interaction,
                                              std::weak_ptr<LlmInterface> llmInterface,
                                              int nestingLevel )
  : WPanel(),
    m_interaction( interaction ),
    m_llmInterface( llmInterface ),
    m_turnContainer( nullptr ),
    m_statusText( nullptr ),
    m_timerText( nullptr ),
    m_menuIcon( nullptr ),
    m_menu( nullptr ),
    m_nestingLevel( nestingLevel ),
    m_isFinished( false )
{
  assert( m_interaction );

  addStyleClass( "LlmInteractionDisplay" );

  // Generate unique timer ID
  m_timerId = "llmTimer_" + to_string( s_nextTimerId++ );

  // Set collapsible to create title bar
  setCollapsible( true );
  setCollapsed( false );

  // Set title with status - create status text widget
  m_statusText = titleBarWidget()->addNew<WText>();
  m_statusText->setTextFormat( TextFormat::Plain );  // Plain text only
  m_statusText->addStyleClass( "LlmStatusIndicator" );

  // Create separate timer widget with the timer ID
  m_timerText = titleBarWidget()->addNew<WText>();
  m_timerText->setId( m_timerId );
  m_timerText->setInline( true );
  m_timerText->setText( "0s" );
  m_timerText->addStyleClass( "LlmStatusIndicator" );

  updateStatus();

  // Create menu icon in title bar
  createMenuIcon();

  // Create main body container
  WContainerWidget *bodyContainer = setCentralWidget( std::make_unique<WContainerWidget>() );

  // Display the user's question/content if non-empty
  // Get content from the InitialRequest turn if it exists
  std::string questionContent;
  const LlmInteractionInitialRequest *initialReq = nullptr;
  if( !interaction->responses.empty() )
  {
    const std::shared_ptr<LlmInteractionTurn> &firstTurn = interaction->responses.front();
    if( firstTurn->type() == LlmInteractionTurn::Type::InitialRequest )
    {
      initialReq = dynamic_cast<const LlmInteractionInitialRequest*>(firstTurn.get());
      if( initialReq )
        questionContent = initialReq->content();
    }
  }

  if( !questionContent.empty() || (initialReq && !initialReq->imageContent().empty()) )
  {
    WContainerWidget *questionDiv = bodyContainer->addNew<WContainerWidget>();
    questionDiv->addStyleClass( "LlmUserQuestion" );

    WText *questionLabel = questionDiv->addNew<WText>( "Question:" );
    questionLabel->addStyleClass( "LlmQuestionLabel" );

    if( !questionContent.empty() )
    {
      // TextFormat::Plain via the 3-arg ctor (format before text) so raw user text skips the XHTML parser.
      WText *questionText = questionDiv->addNew<WText>( questionContent, TextFormat::Plain );
      questionText->addStyleClass( "LlmQuestionContent" );
    }

    // Show image thumbnails if the user attached images
    if( initialReq && !initialReq->imageContent().empty() )
    {
      WContainerWidget *imageStrip = questionDiv->addNew<WContainerWidget>();
      imageStrip->addStyleClass( "LlmUserQuestionImages" );

      for( const LlmToolCall::ImageContent &img : initialReq->imageContent() )
      {
        const string decoded = Wt::Utils::base64Decode( img.base64Data );
        vector<unsigned char> bytes( decoded.begin(), decoded.end() );
        auto imgResource = std::make_shared<WMemoryResource>( img.mimeType );
        imgResource->setData( bytes );

        WImage *thumb = imageStrip->addNew<WImage>( WLink( imgResource ) );
        thumb->addStyleClass( "LlmUserQuestionImageThumb" );
      }
    }
  }

  // Create container for turns
  m_turnContainer = bodyContainer->addNew<WContainerWidget>();

  // Add summary section after the turns
  m_summaryDiv = bodyContainer->addNew<WContainerWidget>();
  m_summaryDiv->addStyleClass( "LlmInteractionSummary" );

  // Add any existing responses
  for( const shared_ptr<LlmInteractionTurn> &turn : interaction->responses )
  {
    handleResponseAdded( turn );
  }

  // Connect to signals
  interaction->responseAdded.connect( this, &LlmInteractionDisplay::handleResponseAdded );
  interaction->subAgentFinished.connect( this, &LlmInteractionDisplay::handleSubAgentFinished );
  interaction->conversationFinished.connect( this, &LlmInteractionDisplay::handleConversationFinished );

  // Check if conversation is already finished
  bool hasError = false;
  bool hasFinalResponse = false;
  for( const shared_ptr<LlmInteractionTurn> &turn : interaction->responses )
  {
    if( turn->type() == LlmInteractionTurn::Type::Error )
      hasError = true;
    if( turn->type() == LlmInteractionTurn::Type::FinalLlmResponse )
      hasFinalResponse = true;
  }

  if( hasError || hasFinalResponse )
  {
    m_isFinished = true;
    updateSummary();
  }else
  {
    // Start status timer
    startStatusTimer();
  }

  // Apply nesting indentation
  if( m_nestingLevel > 0 )
  {
    addStyleClass( "LlmSubAgentNesting LlmNestLevel" + to_string(m_nestingLevel) );
  }
}//LlmInteractionDisplay constructor


LlmInteractionDisplay::~LlmInteractionDisplay()
{
  stopStatusTimer();

  // See ~LlmInteractionTurnDisplay(): m_menu is owned by m_menuIcon, not by us.
}//~LlmInteractionDisplay()


shared_ptr<LlmInteraction> LlmInteractionDisplay::interaction() const
{
  return m_interaction;
}//interaction()


void LlmInteractionDisplay::handleResponseAdded( shared_ptr<LlmInteractionTurn> turn )
{
  // Create appropriate display widget for this turn.  Rendering must never throw back into the
  // conversation engine: this signal is emitted synchronously from LlmInterface's turn processing,
  // and an escaped exception here would skip the engine's terminal-signal emit and hang the
  // conversation (GUI input stays disabled; the benchmark runner blocks forever).  So catch any
  // display failure, fall back to a plain-text placeholder, and keep the finished-state bookkeeping.
  LlmInteractionTurnDisplay *turnDisplay = nullptr;

  try
  {
  switch( turn->type() )
  {
    case LlmInteractionTurn::Type::InitialRequest:
    {
      // Initial request is displayed separately in the conversation header/summary
      // We don't need a separate turn display for it
      // Just update status and return early
      updateStatus();
      return;
    }

    case LlmInteractionTurn::Type::FinalLlmResponse:
    {
      shared_ptr<LlmInteractionFinalResponse> response =
        dynamic_pointer_cast<LlmInteractionFinalResponse>( turn );
      turnDisplay = m_turnContainer->addNew<LlmInteractionFinalResponseDisplay>( response );

      // Final response means conversation is finished
      if( !m_isFinished )
      {
        m_isFinished = true;
        stopStatusTimer();
        updateSummary();
      }
      break;
    }

    case LlmInteractionTurn::Type::ToolCall:
    {
      shared_ptr<LlmToolRequest> request =
        dynamic_pointer_cast<LlmToolRequest>( turn );
      turnDisplay = m_turnContainer->addNew<LlmToolRequestDisplay>( request );
      break;
    }

    case LlmInteractionTurn::Type::ToolResult:
    {
      shared_ptr<LlmToolResults> results =
        dynamic_pointer_cast<LlmToolResults>( turn );
      turnDisplay = m_turnContainer->addNew<LlmToolResultsDisplay>( results, m_llmInterface, m_nestingLevel );
      break;
    }

    case LlmInteractionTurn::Type::Error:
    {
      shared_ptr<LlmInteractionError> error =
        dynamic_pointer_cast<LlmInteractionError>( turn );
      LlmInteractionErrorDisplay *errorDisplay =
            m_turnContainer->addNew<LlmInteractionErrorDisplay>( error, m_llmInterface );

      turnDisplay = errorDisplay;

      // Error means conversation is finished (at least for now)
      if( !m_isFinished )
      {
        m_isFinished = true;
        stopStatusTimer();
        updateSummary();
      }
      break;
    }

    case LlmInteractionTurn::Type::AutoReply:
    {
      // AutoReply is an automatic prompt to the LLM to continue
      // We'll display it similar to a final response, showing the prompt message
      shared_ptr<LlmInteractionAutoReply> autoReply =
        dynamic_pointer_cast<LlmInteractionAutoReply>( turn );

      // Reuse the FinalResponse display class since AutoReply also has content
      // We create a temporary FinalResponse wrapper to use the existing display code
      shared_ptr<LlmInteractionFinalResponse> tempResponse =
        make_shared<LlmInteractionFinalResponse>( autoReply->content(), turn->conversation().lock() );
      tempResponse->setThinkingContent( autoReply->thinkingContent() );
      tempResponse->setRawContent( autoReply->rawContent() );

      turnDisplay = m_turnContainer->addNew<LlmInteractionFinalResponseDisplay>( tempResponse );

      // AutoReply doesn't mean conversation is finished - we're waiting for LLM response
      break;
    }
  }//switch( turn->type() )
  }catch( const std::exception &e )
  {
    cerr << "LlmInteractionDisplay::handleResponseAdded: failed to render turn: " << e.what() << endl;

    WText *fallback = m_turnContainer->addNew<WText>();
    fallback->setTextFormat( Wt::TextFormat::Plain );
    fallback->setText( string("[Turn could not be displayed: ") + e.what() + "]" );

    // A FinalLlmResponse or Error turn still terminates the conversation - keep that bookkeeping
    // consistent even though its widget failed to build.
    const LlmInteractionTurn::Type type = turn->type();
    if( !m_isFinished
        && ((type == LlmInteractionTurn::Type::FinalLlmResponse)
            || (type == LlmInteractionTurn::Type::Error)) )
    {
      m_isFinished = true;
      stopStatusTimer();
      updateSummary();
    }
  }//try / catch

  // Update status
  updateStatus();
}//handleResponseAdded(...)


void LlmInteractionDisplay::handleSubAgentFinished( shared_ptr<LlmInteractionTurn> turn,
                                                    shared_ptr<LlmInteraction> subAgent )
{
  // Sub-agent finished - the display should already be updated via its own signals
  updateStatus();
}//handleSubAgentFinished(...)


void LlmInteractionDisplay::handleConversationFinished()
{
  m_isFinished = true;
  stopStatusTimer();
  updateSummary();
  updateStatus();
  
  
}//handleConversationFinished()


void LlmInteractionDisplay::updateStatus()
{
  const WString titleWithStatus = getTitleWithStatus();
  //cerr << "LlmInteractionDisplay::updateStatus(): Setting title to: '"
  //     << titleWithStatus.toUTF8() << "', m_isFinished=" << m_isFinished << endl;
  m_statusText->setText( titleWithStatus );

  // Show/hide timer based on finished status
  if( m_timerText )
  {
    m_timerText->setHidden( false ); // Make sure it's visible when running
    
    if( m_isFinished && m_interaction )
    {
      auto endTime = m_interaction->timestamp;
      if( m_interaction->finishTime.has_value() )
      {
        endTime = m_interaction->finishTime.value();
      }else
      {
        for( const shared_ptr<LlmInteractionTurn> &turn : m_interaction->responses )
        {
          if( turn->timestamp() > endTime )
            endTime = turn->timestamp();
        }
      }
      
      const auto duration = chrono::duration_cast<chrono::milliseconds>( endTime - m_interaction->timestamp );
      const string durationStr = formatDuration( duration );
      m_timerText->setText( "(total time: " + durationStr + ")" );
    }//if( m_isFinished && m_interaction )
  }
}//updateStatus()


void LlmInteractionDisplay::updateSummary()
{
  assert( m_summaryDiv );
  if( !m_summaryDiv )
    return;

  m_summaryDiv->clear();

  // Start time
  const string startTimeStr = formatTime( m_interaction->timestamp );
  WText *startTime = m_summaryDiv->addNew<WText>( "Start Time: " + startTimeStr );
  startTime->addStyleClass( "LlmSummaryItem" );

  // Total duration (if conversation is finished)
  if( m_isFinished )
  {
    // Calculate total duration from start to last response
    auto endTime = m_interaction->timestamp;
    for( const shared_ptr<LlmInteractionTurn> &turn : m_interaction->responses )
    {
      if( turn->timestamp() > endTime )
        endTime = turn->timestamp();
    }

    const auto duration = chrono::duration_cast<chrono::milliseconds>( endTime - m_interaction->timestamp );
    const string durationStr = formatDuration( duration );

    WText *totalDuration = m_summaryDiv->addNew<WText>( " | Total Duration: " + durationStr );
    totalDuration->addStyleClass( "LlmSummaryItem" );
  }

  // Token usage.  Note the input count is cumulative across all of this conversation's API calls
  // (the growing history is re-sent each tool-call round), which is why it can be large; the cache
  // breakdown below shows how much of that input was re-processed cheaply from the provider's cache.
  if( m_interaction->totalTokens.has_value() && m_interaction->totalTokens.value() > 0 )
  {
    const size_t promptTok = m_interaction->promptTokens.value_or(0);
    const size_t completionTok = m_interaction->completionTokens.value_or(0);

    const string tokenStr = " | Tokens: " + to_string(promptTok) + " in + "
                         + to_string(completionTok) + " out = "
                         + to_string(m_interaction->totalTokens.value()) + " total";
    WText *tokens = m_summaryDiv->addNew<WText>( tokenStr );
    tokens->addStyleClass( "LlmSummaryItem" );

    // Cache usage: how much of the (cumulative) input was served from / written to cache.
    if( m_interaction->cachedTokens.has_value() )
    {
      const size_t cachedTok = m_interaction->cachedTokens.value();
      string cacheStr = " | Cache: " + to_string(cachedTok) + " read";
      if( promptTok > 0 )
      {
        const int pct = static_cast<int>( (100.0 * cachedTok / promptTok) + 0.5 );
        cacheStr += " (" + to_string(pct) + "% of input)";
      }
      const size_t writtenTok = m_interaction->cacheCreationTokens.value_or(0);
      if( writtenTok > 0 )
        cacheStr += ", " + to_string(writtenTok) + " written";

      WText *cache = m_summaryDiv->addNew<WText>( cacheStr );
      cache->addStyleClass( "LlmSummaryItem" );
    }
  }
}//updateSummary()


WString LlmInteractionDisplay::getTitleWithStatus() const
{
  ostringstream oss;

  // User question (truncated) - get from InitialRequest turn if available
  std::string questionContent;
  if( !m_interaction->responses.empty() )
  {
    const std::shared_ptr<LlmInteractionTurn> &firstTurn = m_interaction->responses.front();
    if( firstTurn->type() == LlmInteractionTurn::Type::InitialRequest )
    {
      const LlmInteractionInitialRequest *initialReq = dynamic_cast<const LlmInteractionInitialRequest*>(firstTurn.get());
      if( initialReq )
        questionContent = initialReq->content();
    }
  }

  const string question = truncateString( questionContent, 60 );
  oss << question;

  // Status
  if( m_isFinished )
  {
    // Show completion status without timer
    oss << " - Complete";
  }else
  {
    // Check if performing tool calls
    bool hasToolRequest = false;
    for( const shared_ptr<LlmInteractionTurn> &turn : m_interaction->responses )
    {
      if( turn->type() == LlmInteractionTurn::Type::ToolCall )
      {
        hasToolRequest = true;
        break;
      }
    }

    if( hasToolRequest )
    {
      // Check if we have results for all tool calls
      bool hasToolResults = false;
      for( const shared_ptr<LlmInteractionTurn> &turn : m_interaction->responses )
      {
        if( turn->type() == LlmInteractionTurn::Type::ToolResult )
        {
          hasToolResults = true;
          break;
        }
      }

      if( !hasToolResults )
      {
        oss << " - Performing tool calls...";
      }else
      {
        oss << " - Pending LLM response...";
      }
    }else
    {
      oss << " - Pending LLM response...";
    }
  }

  return WString( oss.str() );
}//getTitleWithStatus()


void LlmInteractionDisplay::startStatusTimer()
{
  // Start JavaScript timer to update the timer display
  const auto now = chrono::system_clock::now();
  const long long startMs = chrono::duration_cast<chrono::milliseconds>(
    now.time_since_epoch() ).count();

  //cerr << "LlmInteractionDisplay::startStatusTimer(): Starting timer with ID='" << m_timerId
  //     << "', startMs=" << startMs << endl;

  ostringstream js;
  js << "console.log('Starting timer for ID: " << m_timerId << ", startMs: " << startMs << "');"
     << "window.llmTimers = window.llmTimers || {};"
     << "window.llmTimers['" << m_timerId << "'] = {"
     << "  startTime: " << startMs << ","
     << "  interval: setInterval(function() {"
     << "    var elem = document.getElementById('" << m_timerId << "');"
     //<< "    console.log('Timer tick for " << m_timerId << ":', elem ? 'elem found' : 'elem NOT found');"
     << "    if (!elem) {"
     << "      console.warn('Timer " << m_timerId << ": element not found, stopping timer');"
     << "      clearInterval(window.llmTimers['" << m_timerId << "'].interval);"
     << "      delete window.llmTimers['" << m_timerId << "'];"
     << "      return;"
     << "    }"
     << "    var elapsed = Date.now() - " << startMs << ";"
     << "    var seconds = Math.floor(elapsed / 1000);"
     //<< "    console.log('Timer " << m_timerId << ": elapsed=' + elapsed + 'ms, seconds=' + seconds + 's');"
     << "    elem.innerHTML = seconds + 's';"
     << "  }, 1000)"
     << "};"
     //<< "console.log('Timer started for " << m_timerId << ":', window.llmTimers['" << m_timerId << "']);"
  ;

  doJavaScript( js.str() );
}//startStatusTimer()


void LlmInteractionDisplay::stopStatusTimer()
{
  // Calculate final duration
  
  //conversation->finishTime
  std::chrono::system_clock::time_point end_time;
  if( m_interaction->finishTime.has_value() )
    end_time = m_interaction->finishTime.value();
  else
    end_time = chrono::system_clock::now();
  
  const auto duration = chrono::duration_cast<chrono::milliseconds>( end_time - m_interaction->timestamp );
  const string durationStr = formatDuration( duration );

  //cerr << "LlmInteractionDisplay::stopStatusTimer(): Stopping timer ID='" << m_timerId
  //     << "', final duration=" << durationStr << endl;

  // Stop JavaScript timer and set final time
  ostringstream js;
  js << "console.log('Stopping timer for ID: " << m_timerId << ", final duration: " << durationStr << "');"
     << "if (window.llmTimers && window.llmTimers['" << m_timerId << "']) {"
     << "  clearInterval(window.llmTimers['" << m_timerId << "'].interval);"
     << "  delete window.llmTimers['" << m_timerId << "'];"
     << "  var elem = document.getElementById('" << m_timerId << "');"
     //<< "  console.log('Timer " << m_timerId << ": elem for final update:', elem);"
     << "  if (elem) elem.innerHTML = '" << durationStr << ")';"  // Include closing parenthesis
     << "}";

  doJavaScript( js.str() );

  // Also update the widget directly
  if( m_timerText )
    m_timerText->setText( "(total time: " + durationStr + ")" );
}//stopStatusTimer()


void LlmInteractionDisplay::createMenuIcon()
{
  // Create menu button in title bar; the menu is built lazily on first click
  // (showMenu) so interactions whose menu is never opened cost nothing.
  m_menuIcon = titleBarWidget()->addNew<WPushButton>();
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn" );
  m_menuIcon->clicked().preventPropagation();
  m_menuIcon->clicked().connect( this, &LlmInteractionDisplay::showMenu );
}//createMenuIcon()


void LlmInteractionDisplay::showMenu()
{
  if( m_menu )
  {
    m_menu->popup( m_menuIcon, Wt::Orientation::Vertical );
    m_menu->doJavaScript( "Wt.WT.BringAboveDialogs('" + m_menu->id() + "');" );
    return;
  }

  // See LlmInteractionTurnDisplay::showMenu() for why this isnt makePopupMenu().
  PopupDivMenu *menu = makeButtonOwnedPopupMenu( m_menuIcon );
  m_menu = menu;

  // Add "Show Initial Request" option
  PopupDivMenuItem *showInitialRequest = menu->addMenuItem( "Show Initial Request" );
  showInitialRequest->triggered().connect( std::bind( [this]() {
    string requestContent;
    if( !m_interaction->responses.empty() )
    {
      const std::shared_ptr<LlmInteractionTurn> &firstTurn = m_interaction->responses.front();
      if( firstTurn && firstTurn->type() == LlmInteractionTurn::Type::InitialRequest )
        requestContent = firstTurn->rawContent();
    }

    if( requestContent.empty() )
    {
      showJsonDialog( "Initial Request", "{\"note\": \"Initial request content not available\"}", false );
    }
    else
    {
      showJsonDialog( "Initial Request", requestContent, true );
    }
  }));

  // Add "Export Interaction JSON" option
  PopupDivMenuItem *exportJson = menu->addMenuItem( "Export Interaction JSON" );
  exportJson->triggered().connect( std::bind( [this]() {
    // Export just this interaction to JSON
    try
    {
      // Create a JSON array with just this interaction's messages
      nlohmann::json messages = nlohmann::json::array();
      LlmConversationHistory::addConversationToLlmApiHistory( *m_interaction, messages );

      const string jsonStr = messages.dump( 2 );  // Pretty print with 2-space indent
      showJsonDialog( "Interaction JSON", jsonStr, true );
    }
    catch( const exception &e )
    {
      showJsonDialog( "Error", string("{\"error\": \"") + e.what() + "\"}", false );
    }
  }));

  // Add "Ask Followup Question" option for debugging
  PopupDivMenuItem *askFollowup = menu->addMenuItem( "Ask Followup Question" );
  askFollowup->triggered().connect( std::bind( [this](){
    if( !m_interaction )
      return;

    // LlmSubAgentFollowup self-deletes when done
    try
    {
      // Owned by wApp (Wt4 has no parent-taking WObject ctor); it hands itself back via
      //  removeChild() in startDeleteSelf() when the user closes its dialog.
      LlmSubAgentFollowup *followup
            = wApp->addChild( std::make_unique<LlmSubAgentFollowup>( m_interaction ) );
      followup->showDialog();
    }catch( const std::exception &e )
    {
      cerr << "LlmSubAgentFollowup: failed to create: " << e.what() << endl;
    }
  } ) );

  // Note: "Copy Question" menu item removed - clipboard copy only works via showJsonDialog

  m_menu->popup( m_menuIcon );
}//showMenu()


void LlmInteractionDisplay::showJsonDialog( const WString &title,
                                            const string &jsonStr,
                                            bool allowDownload )
{
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>( title );

  // Create text area to display JSON
  WTextArea *jsonArea = dialog->contents()->addNew<WTextArea>();
  jsonArea->setReadOnly( true );
  jsonArea->setText( jsonStr );  // jsonStr is already formatted
  jsonArea->setRows( 20 );
  jsonArea->addStyleClass( "LlmJsonDisplay" );

  // Add download button if allowed
  if( allowDownload )
  {
    // Create a memory resource for the JSON; parented to the dialog so it is cleaned up with it.
    auto resource = std::make_shared<WMemoryResource>( "application/json" );
    resource->setData( reinterpret_cast<const unsigned char*>(jsonStr.c_str()),
                      static_cast<int>(jsonStr.length()) );

    // Generate filename from title with timestamp from the interaction
    const time_t tt = chrono::system_clock::to_time_t( m_interaction->timestamp );
    const tm *local_tm = localtime( &tt );
    char timeBuffer[64];
    strftime( timeBuffer, sizeof(timeBuffer), "%Y%m%d_%H%M%S", local_tm );

    string filename = "interaction";
    if( !title.empty() )
    {
      string titleStr = title.toUTF8();
      // Replace spaces with underscores
      for( char &c : titleStr )
      {
        if( c == ' ' )
          c = '_';
      }
      filename = titleStr;
    }
    filename = filename + "_" + timeBuffer + ".json";

    resource->suggestFileName( filename );
    resource->setDispositionType( ContentDisposition::Attachment );

    // On macOS/iOS only clicks on real anchors trigger a download; the `window.open(...)` a
    //  WPushButton link emits is silently dropped by the WKWebView.
#if( BUILD_AS_OSX_APP || IOS )
    WLink dlLink( resource );
    dlLink.setTarget( Wt::LinkTarget::NewWindow );
    WAnchor *downloadBtn = dialog->footer()->addNew<WAnchor>( dlLink );
    downloadBtn->setStyleClass( "simple-dialog-btn DownloadLink" );
    downloadBtn->setText( "Download" );
#else
    WPushButton *downloadBtn = dialog->footer()->addNew<WPushButton>( "Download" );
    downloadBtn->setIcon( "InterSpec_resources/images/download_small.svg" );
    { WLink lnk( resource ); lnk.setTarget( Wt::LinkTarget::NewWindow ); downloadBtn->setLink( lnk ); }
#endif
  }

  // Add copy button
  LOAD_JAVASCRIPT(wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmCopyTextToClipboard );

  WPushButton *copyBtn = dialog->footer()->addNew<WPushButton>( "Copy" );
  const string copyJs = "function(s,e){ Wt.WT.LlmCopyTextToClipboard(s,e,'" + jsonArea->id() + "'); }";
  copyBtn->clicked().connect( copyJs );

  // Add close button
  WPushButton *closeBtn = dialog->footer()->addNew<WPushButton>( "Close" );
  closeBtn->clicked().connect( dialog, &WDialog::accept );

  dialog->show();
}//showJsonDialog(...)


// No longer needed - clipboard copy now happens directly in client-side JavaScript
// void LlmInteractionDisplay::copyToClipboard( const string &text )
// {
//   // Load the JavaScript function if not already loaded
//   LOAD_JAVASCRIPT(wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmCopyTextToClipboard );
//
//   // Call the JavaScript function directly with the text
//   const string escapedText = WString(text).jsStringLiteral();
//   const string js = "Wt.WT.LlmCopyTextToClipboard(null, null, " + escapedText + ");";
//   doJavaScript( js );
// }//copyToClipboard(...)
