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

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WPanel>
#include <Wt/WServer>
#include <Wt/WTextArea>
#include <Wt/WResource>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WMemoryResource>
#include <Wt/WContainerWidget>
#include <Wt/WJavaScriptPreamble>

#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmInteractionDisplay.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace Wt;
using namespace std;

namespace
{
WT_DECLARE_WT_MEMBER
(LlmCopyTextToClipboard, Wt::JavaScriptFunction, "LlmCopyTextToClipboard",
  function( sender, event, textAreaId )
{
  var textArea = document.getElementById(textAreaId);
  if( !textArea ) {
    console.warn('Could not find text area with id:', textAreaId);
    return;
  }

  var text = textArea.value;

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

}//namespace



//
// LlmInteractionTurnDisplay implementation
//

LlmInteractionTurnDisplay::LlmInteractionTurnDisplay( shared_ptr<LlmInteractionTurn> turn,
                                                      WContainerWidget *parent )
  : WPanel( parent ),
    m_turn( turn ),
    m_menuIcon( nullptr ),
    m_bodyContainer( nullptr ),
    m_bodyCreated( false )
{
  addStyleClass( "LlmInteractionTurnDisplay" );

  // Create empty body container - content will be created on construction for now
  m_bodyContainer = new WContainerWidget();
  setCentralWidget( m_bodyContainer );

  // Note: Lazy loading will be added later
  // For now, we create content immediately in derived constructors
}//LlmInteractionTurnDisplay constructor


LlmInteractionTurnDisplay::~LlmInteractionTurnDisplay()
{
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
  // Create menu button in title bar
  m_menuIcon = new WPushButton( titleBarWidget() );
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn" );
  m_menuIcon->clicked().preventPropagation();

  // Create popup menu
  PopupDivMenu *menu = new PopupDivMenu( m_menuIcon, PopupDivMenu::TransientMenu );

  // Add "Show Full Content" menu item for all turn types
  PopupDivMenuItem *showFullContent = menu->addMenuItem( "Show Full Content" );
  showFullContent->triggered().connect( std::bind( [this]() {
    const string rawContent = m_turn->rawContent();
    if( rawContent.empty() )
    {
      showJsonDialog( "Full Content", "{\"note\": \"Content not available\"}", false );
    }else
    {
      showJsonDialog( "Full Content", rawContent, true );
    }
  }));

  // Add "Show Request JSON" option
  PopupDivMenuItem *showRequestJson = menu->addMenuItem( "Show Request JSON" );
  showRequestJson->triggered().connect( std::bind( [this]() {
    const string rawContent = m_turn->rawContent();
    if( rawContent.empty() )
    {
      // Try to get from parent conversation
      if( const shared_ptr<LlmInteraction> convo = m_turn->conversation().lock() )
      {
        // TODO: Build request JSON from conversation context
        showJsonDialog( "Request JSON", "{\"note\": \"Request JSON not available\"}", false );
      }else
      {
        showJsonDialog( "Request JSON", "{\"note\": \"Request JSON not available\"}", false );
      }
    }else
    {
      showJsonDialog( "Request JSON", rawContent, true );
    }
  }));

  // Add "Show Response JSON" option
  PopupDivMenuItem *showResponseJson = menu->addMenuItem( "Show Response JSON" );
  showResponseJson->triggered().connect( std::bind( [this]() {
    const string rawContent = m_turn->rawContent();
    if( rawContent.empty() )
    {
      showJsonDialog( "Response JSON", "{\"note\": \"Response JSON not available\"}", false );
    }else
    {
      showJsonDialog( "Response JSON", rawContent, true );
    }
  }));
}//createMenuIcon()


void LlmInteractionTurnDisplay::showJsonDialog( const WString &title,
                                                const string &jsonStr,
                                                bool allowDownload )
{
  SimpleDialog *dialog = new SimpleDialog( title );

  // Format and display JSON
  const string formattedJson = formatJsonString( jsonStr );

  WTextArea *jsonArea = new WTextArea( dialog->contents() );
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
    // Create memory resource for download
    WMemoryResource *resource = new WMemoryResource( "application/json" );
    resource->setData( reinterpret_cast<const unsigned char*>(formattedJson.data()),
                      static_cast<int>(formattedJson.size()) );

    // Generate filename with timestamp from the interaction turn
    const time_t tt = chrono::system_clock::to_time_t( m_turn->timestamp() );
    const tm *local_tm = localtime( &tt );
    char timeBuffer[64];
    strftime( timeBuffer, sizeof(timeBuffer), "%Y%m%d_%H%M%S", local_tm );
    const string filename = string("llm_interaction_") + timeBuffer + ".json";
    resource->suggestFileName( filename );
    resource->setDispositionType( WResource::Attachment );

    WPushButton *downloadBtn = dialog->addButton( "Download" );
    downloadBtn->setIcon( "InterSpec_resources/images/download_small.svg" );
    downloadBtn->setLink( WLink(resource) );
    downloadBtn->setLinkTarget( Wt::TargetNewWindow );
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
  shared_ptr<LlmInteractionFinalResponse> response,
  WContainerWidget *parent )
  : LlmInteractionTurnDisplay( response, parent ),
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


WString LlmInteractionFinalResponseDisplay::getTitleText() const
{
  const string content = m_response->content();
  const string firstLine = getFirstLine( content );
  const string truncated = truncateString( firstLine, 80 );

  return WString( "Response: {1}" ).arg( truncated );
}//getTitleText()


void LlmInteractionFinalResponseDisplay::createBodyContent()
{
  // Display the full response content
  WText *contentText = new WText( m_bodyContainer );
  contentText->setText( m_response->content() );
  contentText->setTextFormat( PlainText );
  contentText->addStyleClass( "LlmResponseContent" );

  // Add metadata section
  WContainerWidget *metadataDiv = new WContainerWidget( m_bodyContainer );
  metadataDiv->addStyleClass( "LlmMetadata" );

  // Timestamp
  const string timeStr = formatTime( m_response->timestamp() );
  WText *timeText = new WText( "Time: " + timeStr, metadataDiv );
  timeText->addStyleClass( "LlmTimestamp" );

  // Duration
  if( m_response->callDuration().has_value() )
  {
    const string durationStr = formatDuration( m_response->callDuration().value() );
    WText *durationText = new WText( " | Duration: " + durationStr, metadataDiv );
    durationText->addStyleClass( "LlmDuration" );
  }

  // Token usage (if available from parent conversation)
  if( const shared_ptr<LlmInteraction> convo = m_response->conversation().lock() )
  {
    if( convo->promptTokens.has_value() )
    {
      const string tokenStr = " | Tokens: " + to_string(convo->promptTokens.value())
                           + " + " + to_string(convo->completionTokens.value_or(0))
                           + " = " + to_string(convo->totalTokens.value_or(0));
      WText *tokenText = new WText( tokenStr, metadataDiv );
      tokenText->addStyleClass( "LlmTokenUsage" );
    }
  }

  // Thinking content (if available)
  if( !m_response->thinkingContent().empty() )
  {
    WContainerWidget *thinkingDiv = new WContainerWidget( m_bodyContainer );
    thinkingDiv->addStyleClass( "LlmThinkingContent" );

    WText *thinkingLabel = new WText( "Thinking:", thinkingDiv );
    thinkingLabel->addStyleClass( "LlmThinkingLabel" );

    WText *thinkingText = new WText( m_response->thinkingContent(), thinkingDiv );
    thinkingText->setTextFormat( PlainText );
  }
}//createBodyContent()


//
// LlmToolRequestDisplay implementation
//

LlmToolRequestDisplay::LlmToolRequestDisplay( shared_ptr<LlmToolRequest> request,
                                              WContainerWidget *parent )
  : LlmInteractionTurnDisplay( request, parent ),
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

  const string truncated = truncateString( oss.str(), 100 );
  return WString( truncated );
}//getTitleText()


void LlmToolRequestDisplay::createBodyContent()
{
  const vector<LlmToolCall> &calls = m_request->toolCalls();

  if( calls.empty() )
  {
    WText *emptyText = new WText( "No tool calls", m_bodyContainer );
    emptyText->addStyleClass( "LlmEmptyMessage" );
    return;
  }

  // Display each tool call
  for( const LlmToolCall &call : calls )
  {
    WContainerWidget *callDiv = new WContainerWidget( m_bodyContainer );
    callDiv->addStyleClass( "LlmToolCallItem" );

    // Tool name
    WText *nameText = new WText( "<b>" + call.toolName + "</b>", callDiv );
    nameText->setTextFormat( TextFormat::XHTMLUnsafeText );
    nameText->addStyleClass( "LlmToolName" );

    // Invocation ID
    WText *idText = new WText( " (ID: " + call.invocationId + ")", callDiv );
    idText->addStyleClass( "LlmInvocationId" );

    // Parameters
    if( !call.toolParameters.empty() )
    {
      WText *paramsLabel = new WText( "Parameters:", callDiv );
      paramsLabel->addStyleClass( "LlmParamsLabel" );

      WTextArea *paramsArea = new WTextArea( callDiv );
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
                                              int nestingLevel,
                                              WContainerWidget *parent )
  : LlmInteractionTurnDisplay( results, parent ),
    m_results( results ),
    m_nestingLevel( nestingLevel )
{
  addStyleClass( "LlmToolResultsDisplay" );

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

  const string truncated = truncateString( oss.str(), 100 );
  return WString( truncated );
}//getTitleText()


void LlmToolResultsDisplay::createBodyContent()
{
  const vector<LlmToolCall> &calls = m_results->toolCalls();

  if( calls.empty() )
  {
    WText *emptyText = new WText( "No tool results", m_bodyContainer );
    emptyText->addStyleClass( "LlmEmptyMessage" );
    return;
  }

  // Display each tool result
  for( const LlmToolCall &call : calls )
  {
    WContainerWidget *resultDiv = new WContainerWidget( m_bodyContainer );
    resultDiv->addStyleClass( "LlmToolResultItem" );

    // Add error styling if tool call had an error
    if( call.status == LlmToolCall::CallStatus::Error )
      resultDiv->addStyleClass( "LlmToolResultError" );

    // Tool name with status indicator
    ostringstream nameHtml;
    nameHtml << "<b>" << call.toolName << "</b>";
    if( call.status == LlmToolCall::CallStatus::Error )
      nameHtml << " <span class='LlmToolStatusError'>[ERROR]</span>";
    else if( call.status == LlmToolCall::CallStatus::Pending  && !call.sub_agent_conversation )
      nameHtml << " <span class='LlmToolStatusPending'>[PENDING]</span>";

    WText *nameText = new WText( nameHtml.str(), resultDiv );
    nameText->setTextFormat( TextFormat::XHTMLUnsafeText );
    nameText->addStyleClass( "LlmToolName" );

    // Invocation ID
    WText *idText = new WText( " (ID: " + call.invocationId + ")", resultDiv );
    idText->addStyleClass( "LlmInvocationId" );

    // Execution duration
    if( call.executionDuration.has_value() )
    {
      const string durationStr = formatDuration( call.executionDuration.value() );
      WText *durationText = new WText( " - took " + durationStr, resultDiv );
      durationText->addStyleClass( "LlmDuration" );
    }

    // Sub-agent conversation (if present)
    if( call.sub_agent_conversation )
    {
        WText *subAgentLabel = new WText( "Sub-agent conversation:", resultDiv );
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
              WText *contextLabel = new WText( "Context:", resultDiv );
              contextLabel->addStyleClass( "LlmParamsLabel" );

              WTextArea *contextArea = new WTextArea( resultDiv );
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
              WText *taskLabel = new WText( "Task:", resultDiv );
              taskLabel->addStyleClass( "LlmParamsLabel" );

              WTextArea *taskArea = new WTextArea( resultDiv );
              taskArea->setReadOnly( true );
              taskArea->setText( task );
              taskArea->setRows( 3 );
              taskArea->addStyleClass( "LlmToolParams" );
            }
          }
        }
       */

        // Create nested LlmInteractionDisplay for sub-agent
        LlmInteractionDisplay *subAgentDisplay =
          new LlmInteractionDisplay( call.sub_agent_conversation, m_nestingLevel + 1, resultDiv );
        subAgentDisplay->addStyleClass( "LlmSubAgentNesting" );

        m_subAgentDisplays.push_back( subAgentDisplay );

        // Store the result div for later when sub-agent finishes
        m_subAgentResultDivs[call.sub_agent_conversation] = resultDiv;

        // Connect to sub-agent finished signal to update results
        call.sub_agent_conversation->conversationFinished.connect(
          boost::bind( &LlmToolResultsDisplay::handleSubAgentFinished,
                    this, m_turn, call.sub_agent_conversation ) );

      // Make sure we expand the parent widget so we can see the sub-agent conversation going on.
      if( isCollapsed() )
        setCollapsed( false );
    }//if( call.sub_agent_conversation )

    // Result content - show if available and (no sub-agent OR sub-agent has completed)
    if( !call.content.empty() && (!call.sub_agent_conversation || call.sub_agent_conversation->finishTime.has_value()) )
    {
      WText *resultLabel = new WText( "Result:", resultDiv );
      resultLabel->addStyleClass( "LlmResultLabel" );

      WTextArea *resultArea = new WTextArea( resultDiv );
      resultArea->setReadOnly( true );
      resultArea->setText( call.content );
      resultArea->setRows( 5 );
      resultArea->addStyleClass( "LlmToolResult" );
    }
  }
}//createBodyContent()


void LlmToolResultsDisplay::createMenuIcon()
{
  // Call base class implementation first
  LlmInteractionTurnDisplay::createMenuIcon();

  // Find the popup menu that was created by the base class
  if( !m_menuIcon )
    return;
  
  

  // Get the popup menu from the button's children
  PopupDivMenu *menu = dynamic_cast<PopupDivMenu *>( m_menuIcon->menu() );
  assert( menu );
  if( !menu )
    return;

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
}//createMenuIcon()


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
      WText *resultLabel = new WText( "Result:", resultDiv );
      resultLabel->addStyleClass( "LlmResultLabel" );

      WTextArea *resultArea = new WTextArea( resultDiv );
      resultArea->setReadOnly( true );
      resultArea->setText( call.content );
      resultArea->setRows( 5 );
      resultArea->addStyleClass( "LlmToolResult" );

      break;
    }
  }
}//handleSubAgentFinished(...)


//
// LlmInteractionErrorDisplay implementation
//

LlmInteractionErrorDisplay::LlmInteractionErrorDisplay( shared_ptr<LlmInteractionError> error,
                                                        WContainerWidget *parent )
  : LlmInteractionTurnDisplay( error, parent ),
    m_error( error ),
    m_retryCallback( nullptr )
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


void LlmInteractionErrorDisplay::setRetryCallback( function<void()> callback )
{
  m_retryCallback = callback;
}//setRetryCallback(...)


WString LlmInteractionErrorDisplay::getTitleText() const
{
  const string errorMsg = m_error->errorMessage();
  const string firstLine = getFirstLine( errorMsg );
  const string truncated = truncateString( firstLine, 80 );

  return WString( "Error: {1}" ).arg( truncated );
}//getTitleText()


void LlmInteractionErrorDisplay::createBodyContent()
{
  // Display error message
  WText *errorText = new WText( m_bodyContainer );
  errorText->setText( m_error->errorMessage() );
  errorText->setTextFormat( PlainText );
  errorText->addStyleClass( "LlmErrorMessage" );

  // Add retry button
  WPushButton *retryBtn = new WPushButton( "Retry", m_bodyContainer );
  retryBtn->addStyleClass( "LlmRetryButton" );
  retryBtn->clicked().connect( this, &LlmInteractionErrorDisplay::handleRetry );
}//createBodyContent()


void LlmInteractionErrorDisplay::handleRetry()
{
  if( m_retryCallback )
    m_retryCallback();
}//handleRetry()


//
// LlmInteractionDisplay implementation
//

int LlmInteractionDisplay::s_nextTimerId = 0;

LlmInteractionDisplay::LlmInteractionDisplay( shared_ptr<LlmInteraction> interaction,
                                              int nestingLevel,
                                              WContainerWidget *parent )
  : WPanel( parent ),
    m_interaction( interaction ),
    m_turnContainer( nullptr ),
    m_statusText( nullptr ),
    m_timerText( nullptr ),
    m_menuIcon( nullptr ),
    m_retryCallback( nullptr ),
    m_nestingLevel( nestingLevel ),
    m_isFinished( false )
{
  addStyleClass( "LlmInteractionDisplay" );

  // Generate unique timer ID
  m_timerId = "llmTimer_" + to_string( s_nextTimerId++ );

  // Set collapsible to create title bar
  setCollapsible( true );
  setCollapsed( false );

  // Set title with status - create status text widget
  m_statusText = new WText( titleBarWidget() );
  m_statusText->setTextFormat( PlainText );  // Plain text only
  m_statusText->addStyleClass( "LlmStatusIndicator" );

  // Create separate timer widget with the timer ID
  m_timerText = new WText( titleBarWidget() );
  m_timerText->setId( m_timerId );
  m_timerText->setInline( true );
  m_timerText->setText( "0s" );
  m_timerText->addStyleClass( "LlmStatusIndicator" );

  updateStatus();

  // Create menu icon in title bar
  createMenuIcon();

  // Create main body container
  WContainerWidget *bodyContainer = new WContainerWidget();
  setCentralWidget( bodyContainer );

  // Display the user's question/content if non-empty
  if( !interaction->content.empty() )
  {
    WContainerWidget *questionDiv = new WContainerWidget( bodyContainer );
    questionDiv->addStyleClass( "LlmUserQuestion" );

    WText *questionLabel = new WText( "Question:", questionDiv );
    questionLabel->addStyleClass( "LlmQuestionLabel" );

    WText *questionText = new WText( interaction->content, questionDiv );
    questionText->setTextFormat( PlainText );
    questionText->addStyleClass( "LlmQuestionContent" );
  }

  // Create container for turns
  m_turnContainer = new WContainerWidget( bodyContainer );

  // Add summary section after the turns
  m_summaryDiv = new WContainerWidget( bodyContainer );
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
}//~LlmInteractionDisplay()


shared_ptr<LlmInteraction> LlmInteractionDisplay::interaction() const
{
  return m_interaction;
}//interaction()


void LlmInteractionDisplay::setRetryCallback( function<void(shared_ptr<LlmInteraction>)> callback )
{
  m_retryCallback = callback;
}//setRetryCallback(...)


void LlmInteractionDisplay::handleResponseAdded( shared_ptr<LlmInteractionTurn> turn )
{
  // Create appropriate display widget for this turn
  LlmInteractionTurnDisplay *turnDisplay = nullptr;

  switch( turn->type() )
  {
    case LlmInteractionTurn::Type::FinalLlmResponse:
    {
      shared_ptr<LlmInteractionFinalResponse> response =
        dynamic_pointer_cast<LlmInteractionFinalResponse>( turn );
      turnDisplay = new LlmInteractionFinalResponseDisplay( response, m_turnContainer );

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
      turnDisplay = new LlmToolRequestDisplay( request, m_turnContainer );
      break;
    }

    case LlmInteractionTurn::Type::ToolResult:
    {
      shared_ptr<LlmToolResults> results =
        dynamic_pointer_cast<LlmToolResults>( turn );
      turnDisplay = new LlmToolResultsDisplay( results, m_nestingLevel, m_turnContainer );
      break;
    }

    case LlmInteractionTurn::Type::Error:
    {
      shared_ptr<LlmInteractionError> error =
        dynamic_pointer_cast<LlmInteractionError>( turn );
      LlmInteractionErrorDisplay *errorDisplay =
        new LlmInteractionErrorDisplay( error, m_turnContainer );

      // Set retry callback
      errorDisplay->setRetryCallback( [this]() {
        if( m_retryCallback )
          m_retryCallback( m_interaction );
      });

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

      turnDisplay = new LlmInteractionFinalResponseDisplay( tempResponse, m_turnContainer );

      // AutoReply doesn't mean conversation is finished - we're waiting for LLM response
      break;
    }
  }//switch( turn->type() )

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
  WText *startTime = new WText( "Start Time: " + startTimeStr, m_summaryDiv );
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

    WText *totalDuration = new WText( " | Total Duration: " + durationStr, m_summaryDiv );
    totalDuration->addStyleClass( "LlmSummaryItem" );
  }

  // Token usage
  if( m_interaction->totalTokens.has_value() && m_interaction->totalTokens.value() > 0 )
  {
    const string tokenStr = " | Tokens: " + to_string(m_interaction->promptTokens.value_or(0))
                         + " + " + to_string(m_interaction->completionTokens.value_or(0))
                         + " = " + to_string(m_interaction->totalTokens.value());
    WText *tokens = new WText( tokenStr, m_summaryDiv );
    tokens->addStyleClass( "LlmSummaryItem" );
  }
}//updateSummary()


WString LlmInteractionDisplay::getTitleWithStatus() const
{
  ostringstream oss;

  // User question (truncated)
  const string question = truncateString( m_interaction->content, 60 );
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

  cerr << "LlmInteractionDisplay::startStatusTimer(): Starting timer with ID='" << m_timerId
       << "', startMs=" << startMs << endl;

  ostringstream js;
  js << "console.log('Starting timer for ID: " << m_timerId << ", startMs: " << startMs << "');"
     << "window.llmTimers = window.llmTimers || {};"
     << "window.llmTimers['" << m_timerId << "'] = {"
     << "  startTime: " << startMs << ","
     << "  interval: setInterval(function() {"
     << "    var elem = document.getElementById('" << m_timerId << "');"
     << "    console.log('Timer tick for " << m_timerId << ":', elem ? 'elem found' : 'elem NOT found');"
     << "    if (!elem) {"
     << "      console.warn('Timer " << m_timerId << ": element not found, stopping timer');"
     << "      clearInterval(window.llmTimers['" << m_timerId << "'].interval);"
     << "      delete window.llmTimers['" << m_timerId << "'];"
     << "      return;"
     << "    }"
     << "    var elapsed = Date.now() - " << startMs << ";"
     << "    var seconds = Math.floor(elapsed / 1000);"
     << "    console.log('Timer " << m_timerId << ": elapsed=' + elapsed + 'ms, seconds=' + seconds + 's');"
     << "    elem.innerHTML = seconds + 's';"
     << "  }, 1000)"
     << "};"
     << "console.log('Timer started for " << m_timerId << ":', window.llmTimers['" << m_timerId << "']);";

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
  // Create menu button in title bar
  m_menuIcon = new WPushButton( titleBarWidget() );
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn" );
  m_menuIcon->clicked().preventPropagation();
  
  // Create popup menu
  PopupDivMenu *menu = new PopupDivMenu( m_menuIcon, PopupDivMenu::TransientMenu );

  // Add "Show Initial Request" option
  PopupDivMenuItem *showInitialRequest = menu->addMenuItem( "Show Initial Request" );
  showInitialRequest->triggered().connect( std::bind( [this]() {
    const string &requestContent = m_interaction->initialRequestContent;
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

  // Note: "Copy Question" menu item removed - clipboard copy only works via showJsonDialog
}//createMenuIcon()


void LlmInteractionDisplay::showJsonDialog( const WString &title,
                                            const string &jsonStr,
                                            bool allowDownload )
{
  SimpleDialog *dialog = new SimpleDialog( title );

  // Create text area to display JSON
  WTextArea *jsonArea = new WTextArea( dialog->contents() );
  jsonArea->setReadOnly( true );
  jsonArea->setText( jsonStr );  // jsonStr is already formatted
  jsonArea->setRows( 20 );
  jsonArea->addStyleClass( "LlmJsonDisplay" );

  // Add download button if allowed
  if( allowDownload )
  {
    // Create a memory resource for the JSON
    WMemoryResource *resource = new WMemoryResource( "application/json" );
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
    resource->setDispositionType( WResource::Attachment );

    WPushButton *downloadBtn = new WPushButton( "Download", dialog->footer() );
    downloadBtn->setIcon( "InterSpec_resources/images/download_small.svg" );
    downloadBtn->setLink( WLink(resource) );
    downloadBtn->setLinkTarget( Wt::TargetNewWindow );
  }

  // Add copy button
  LOAD_JAVASCRIPT(wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmCopyTextToClipboard );

  WPushButton *copyBtn = new WPushButton( "Copy", dialog->footer() );
  const string copyJs = "function(s,e){ Wt.WT.LlmCopyTextToClipboard(s,e,'" + jsonArea->id() + "'); }";
  copyBtn->clicked().connect( copyJs );

  // Add close button
  WPushButton *closeBtn = new WPushButton( "Close", dialog->footer() );
  closeBtn->clicked().connect( dialog, &WDialog::accept );

  dialog->show();
  dialog->finished().connect( std::bind( [dialog]() {
    delete dialog;
  }));
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
