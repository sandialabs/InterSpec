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

#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WPanel.h>
#include <Wt/WTheme.h>
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
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/LlmInteractionDisplay.h"
#include "InterSpec/LlmSubAgentFollowup.h"
#include "InterSpec/LlmConversationHistory.h"

using namespace Wt;
using namespace std;

/** Leaves each tool calls invocation ID (e.g., "(ID: toolu_bdrk_01Dfm4fer8...)") out of the display
 with debug options on; set to 0 to show them.  They are also in the "Show Tool Results JSON" dialog.
 */
#define NO_TOOL_CALL_ID_TEXT 1

namespace
{
  /** Phones get a pared-down conversation display; see LlmInteractionDisplay::m_phoneLayout. */
  bool is_phone_layout()
  {
    const InterSpec * const viewer = InterSpec::instance();
    return viewer && viewer->isPhone();
  }

  /** If the user wants to see the details only of use when developing the assistant (the
   "LlmShowDebugOptions" preference); never on phones.  Only needed for things CSS can't hide by the
   "LlmDebugOnly" class (i.e., popup menu items).
   */
  bool show_debug_options()
  {
    InterSpec * const viewer = InterSpec::instance();
    if( !viewer || viewer->isPhone() )
      return false;

    try
    {
      return UserPreferences::preferenceValue<bool>( "LlmShowDebugOptions", viewer );
    }catch( std::exception & )
    {
    }
    return false;
  }//bool show_debug_options()

  /** Marks the title text of `panel` so CSS can limit it to two lines with an ellipsis - a list of tool
   names can otherwise overflow the title bar, especially in the narrow LLM side panel.
   WPanel gives its title text no style class; it is the first WText in the title bar (after the
   collapse icon, and before our status text).
   */
  void mark_panel_title( WPanel *panel )
  {
    WContainerWidget * const titleBar = panel->titleBarWidget();
    for( int i = 0; titleBar && (i < titleBar->count()); ++i )
    {
      WText * const text = dynamic_cast<WText *>( titleBar->widget(i) );
      if( text && !text->hasStyleClass( "LlmTurnStatus" ) )
      {
        text->addStyleClass( "LlmTurnTitle" );
        break;
      }
    }//for( loop over title bar widgets )
  }//void mark_panel_title( WPanel *panel )

  /** If `call` was cancelled by the user (e.g., the agent it invoked was cancelled), rather than failed;
   see LlmInterface::markOutstandingCallsCancelled(). */
  bool is_cancelled( const LlmToolCall &call )
  {
    if( call.status != LlmToolCall::CallStatus::Error )
      return false;
    const nlohmann::json result = nlohmann::json::parse( call.content, nullptr, false );
    return result.is_object() && (result.value( "status", string() ) == "cancelled");
  }//bool is_cancelled( const LlmToolCall &call )

  /** If `call` failed (other than by being cancelled by the user). */
  bool is_failed( const LlmToolCall &call )
  {
    return (call.status == LlmToolCall::CallStatus::Error) && !is_cancelled( call );
  }

  /** The status marker of a tool call (e.g., "[ERROR]"), as XHTML; empty if it succeeded. */
  string tool_status_html( const LlmToolCall &call )
  {
    if( is_cancelled( call ) )
      return "<span class='LlmToolStatusPending'>[CANCELLED]</span>";
    if( call.status == LlmToolCall::CallStatus::Error )
      return "<span class='LlmToolStatusError'>[ERROR]</span>";
    if( call.status == LlmToolCall::CallStatus::Pending && !call.sub_agent_conversation )
      return "<span class='LlmToolStatusPending'>[PENDING]</span>";
    return "";
  }//string tool_status_html( const LlmToolCall &call )


  /** The arguments the LLM gave a tool call.  If the call could not be run, LlmInterface replaced them
   with {"_raw": <the arguments as sent>, "_execution_error": ...}, so they are parsed back out.
   */
  nlohmann::json tool_call_args( const LlmToolCall &call )
  {
    const nlohmann::json &params = call.toolParameters;
    const nlohmann::json::const_iterator raw_iter = params.is_object() ? params.find( "_raw" ) : params.end();
    if( (raw_iter == params.end()) || !raw_iter->is_string() )
      return params;

    const string &raw = raw_iter->get_ref<const string &>();
    if( raw.empty() )
      return nlohmann::json::object();

    const nlohmann::json parsed = nlohmann::json::parse( raw, nullptr, false );
    return parsed.is_discarded() ? nlohmann::json( raw ) : parsed;
  }//nlohmann::json tool_call_args( const LlmToolCall &call )


  /** Compact JSON with a space after each ':' and ',', e.g. `{"energy": 2103.65, "refit": true}`;
   invalid UTF-8 from the model is replaced, rather than throwing. */
  string compact_json( const nlohmann::json &j )
  {
    if( !j.is_structured() )
      return j.dump( -1, ' ', false, nlohmann::json::error_handler_t::replace );

    string out = j.is_object() ? "{" : "[";
    for( nlohmann::json::const_iterator iter = j.begin(); iter != j.end(); ++iter )
    {
      if( iter != j.begin() )
        out += ", ";
      if( j.is_object() )
        out += compact_json( nlohmann::json( iter.key() ) ) + ": ";
      out += compact_json( iter.value() );
    }
    return out + (j.is_object() ? "}" : "]");
  }//string compact_json( const nlohmann::json &j )


  /** compact_json() of `j`, shortened to about `max_bytes`. */
  string short_json( const nlohmann::json &j, const size_t max_bytes )
  {
    string str = compact_json( j );
    if( str.size() <= max_bytes )
      return str;
    SpecUtils::utf8_limit_str_size( str, max_bytes );
    return str + (j.is_string() ? "...\"" : "...");
  }//string short_json(...)


  /** The names of the arguments of a tool to show first, as they tell the most about a call: for some
   tools, as listed here; and then the tool's required arguments, in the order its definition
   (`schema`, its JSON-schema parameters) lists them.  May include duplicates.
   */
  vector<string> leading_arg_names( const string &toolName, const nlohmann::json &schema )
  {
    static const map<string, vector<string>> s_tool_orders = {
      { "add_analysis_peak", { "source", "energy" } },
      { "edit_analysis_peak", { "energy", "action", "source", "userLabel" } },
      { "set_peak_shape", { "energy", "skewType", "continuumType" } },
      { "get_peaks", { "filter", "lowerEnergy", "upperEnergy" } },
      { "currie_mda_calc", { "nuclide", "energy", "distance" } },
      { "load_spectrum_file", { "path", "specType" } },
      { "load_detector_efficiency_function", { "identifier", "drfName", "detectorName" } },
      { "modify_shielding_source_config", { "operation", "nuclide", "material", "geometry" } },
      { "activity_fit_one_off", { "sources", "peak_energies" } },
      { "peak_based_relative_efficiency", { "sources", "peak_energies" } },
      { "get_spectrum_image", { "energy_range" } }
    };

    vector<string> names;
    const map<string, vector<string>>::const_iterator order_iter = s_tool_orders.find( toolName );
    if( order_iter != s_tool_orders.end() )
      names = order_iter->second;

    const nlohmann::json::const_iterator required_iter = schema.is_object() ? schema.find( "required" ) : schema.end();
    if( (required_iter != schema.end()) && required_iter->is_array() )
    {
      for( const nlohmann::json &name : *required_iter )
      {
        if( name.is_string() )
          names.push_back( name.get<string>() );
      }
    }//if( the tool has required arguments )

    return names;
  }//vector<string> leading_arg_names(...)


  /** A tool call's arguments as (name, value) pairs, in the order to show them: those named in
   `leading` first, then the rest (alphabetically).  Empty if `params` is not a JSON object.
   */
  vector<pair<string, const nlohmann::json *>> ordered_args( const nlohmann::json &params,
                                                              const vector<string> &leading )
  {
    vector<pair<string, const nlohmann::json *>> args;
    if( !params.is_object() )
      return args;

    const auto add = [&args]( const string &name, const nlohmann::json &value ){
      for( const pair<string, const nlohmann::json *> &arg : args )
      {
        if( arg.first == name )
          return;
      }
      args.emplace_back( name, &value );
    };

    for( const string &name : leading )
    {
      const nlohmann::json::const_iterator iter = params.find( name );
      if( iter != params.end() )
        add( name, *iter );
    }

    for( nlohmann::json::const_iterator iter = params.begin(); iter != params.end(); ++iter )
      add( iter.key(), iter.value() );

    return args;
  }//ordered_args(...)


  /** A tool call's arguments as XHTML for its one-line header, `{"energy": 2103.65, "refit": true}`:
   each of `args` (with long names and values shortened) is in a "LlmArg" span, and there is a
   hidden ",..." "LlmArgMore" span before the closing brace, so LlmFitToolArgs can show only the
   arguments that fit.  Empty if there are no arguments.
   */
  string tool_args_html( const nlohmann::json &params, const vector<pair<string, const nlohmann::json *>> &args )
  {
    // E.g., arguments that were not valid JSON; just left to the CSS ellipsis.
    if( !params.is_null() && !params.is_object() )
      return Wt::Utils::htmlEncode( short_json( params, 80 ) );

    if( args.empty() )
      return "";

    string html = "{";
    for( size_t i = 0; i < args.size(); ++i )
    {
      const string item = short_json( nlohmann::json( args[i].first ), 40 ) + ": " + short_json( *args[i].second, 40 );
      html += "<span class=\"LlmArg\">" + Wt::Utils::htmlEncode( (i ? ", " : "") + item ) + "</span>";
    }

    return html + "<span class=\"LlmArgMore\">,...</span><span class=\"LlmArgClose\">}</span>";
  }//string tool_args_html(...)


  /** All of a tool call's arguments, in full, as compact JSON in the order of `args`; for the
   debug-options header, which wraps it onto as many lines as needed.  Empty if there are none.
   */
  string tool_args_full( const nlohmann::json &params, const vector<pair<string, const nlohmann::json *>> &args )
  {
    if( !params.is_null() && !params.is_object() )
      return compact_json( params );

    if( args.empty() )
      return "";

    string text = "{";
    for( size_t i = 0; i < args.size(); ++i )
      text += (i ? ", " : "") + compact_json( nlohmann::json( args[i].first ) ) + ": " + compact_json( *args[i].second );
    return text + "}";
  }//string tool_args_full(...)

  /** The name of a tool for a title: a sub-agent invocation ("invoke_X") is shown as "X agent". */
  string tool_display_name( const string &toolName, bool &isAgent )
  {
    const string prefix = "invoke_";
    isAgent = (toolName.size() > prefix.size()) && (toolName.compare( 0, prefix.size(), prefix ) == 0);
    if( !isAgent )
      return toolName;
    return WString::tr( "lcw-round-agent" ).arg( WString::fromUTF8( toolName.substr( prefix.size() ) ) ).toUTF8();
  }//string tool_display_name( const string &toolName, bool &isAgent )
}//namespace

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


WT_DECLARE_WT_MEMBER
(LlmRenderMarkdown, Wt::JavaScriptFunction, "LlmRenderMarkdown",
  function( elemId, markdown )
{
  // InterSpec_resources/LlmMarkdown.js is pulled in with WApplication::require(), so it may not have
  //  finished loading when this runs.  This function is inlined into the JavaScript stream (not
  //  fetched), so it is always defined; it waits for LlmMarkdown, then shows plain text if the
  //  libraries never arrive.
  var tries = 0;
  var attempt = function(){
    if( typeof LlmMarkdown !== 'undefined' ){
      LlmMarkdown.render( elemId, markdown );
      return;
    }

    tries += 1;
    if( tries > 100 ){  // ~5 seconds
      console.warn( 'LlmRenderMarkdown: LlmMarkdown.js never loaded; showing plain text.' );
      var el = document.getElementById( elemId );
      if( el ){
        // Stash the source too: if LlmMarkdown.js does turn up later, "Show Raw Text" reads it from
        //  here, and would otherwise blank the response.
        el.llmMarkdownSource = markdown;
        el.textContent = markdown;
        el.classList.add( 'LlmMarkdownRaw' );
      }
      return;
    }

    setTimeout( attempt, 50 );
  };

  attempt();
}
);


/** Shows a once-a-second counter in `el`, counting up from `elapsedMs`: `fmt` with its "{1}" replaced
 by the whole seconds elapsed; see LlmElapsedText.
 */
WT_DECLARE_WT_MEMBER
(LlmElapsedStart, Wt::JavaScriptFunction, "LlmElapsedStart",
  function( el, fmt, elapsedMs )
{
  if( !el )
    return;
  if( el._llmTick )
    clearInterval( el._llmTick );

  var t0 = Date.now() - elapsedMs;
  var update = function(){
    // The element was removed (e.g., the conversation was cleared) - so stop.
    if( !el.isConnected ){
      clearInterval( el._llmTick );
      el._llmTick = null;
      return;
    }
    el.textContent = fmt.replace( "{1}", Math.floor( (Date.now() - t0) / 1000 ) );
  };
  el._llmTick = setInterval( update, 1000 );
  update();
}
);


/** Stops the LlmElapsedStart counter in `el`, and shows `text` - Wt sends a WText's new text before
 the JS queued in the same response, so a start earlier in that response would otherwise overwrite it.
 */
WT_DECLARE_WT_MEMBER
(LlmElapsedStop, Wt::JavaScriptFunction, "LlmElapsedStop",
  function( el, text )
{
  if( !el )
    return;
  if( el._llmTick ){
    clearInterval( el._llmTick );
    el._llmTick = null;
  }
  el.textContent = text;
}
);


/** Shows as many whole arguments as fit on a tool call's one-line header (see
 LlmToolResultsDisplay::renderCallHeader()), ending in ",...}" if some don't.  Done when the header is
 first shown (it is often made in a collapsed panel), and again if its width changes in the following
 second; if it later gets narrower, CSS cuts it with an ellipsis.
 */
WT_DECLARE_WT_MEMBER
(LlmFitToolArgs, Wt::JavaScriptFunction, "LlmFitToolArgs",
  function( header )
{
  var el = header ? header.querySelector( '.LlmToolArgs' ) : null;
  var more = el ? el.querySelector( '.LlmArgMore' ) : null;
  var close = el ? el.querySelector( '.LlmArgClose' ) : null;
  if( !more || !close )
    return;

  var tries = 0, fitWidth = -1;
  var fit = function(){
    // LlmToolGui.css may not have been applied yet, the first time the assistant is shown.
    if( getComputedStyle( header ).display !== 'flex' ){
      if( ++tries < 50 )
        setTimeout( fit, 100 );
      return;
    }

    // Lay out every argument, and the ",...", then see where each ends; if they don't all fit, the
    //  arguments element is shrunk to the width available to it.
    var args = el.querySelectorAll( '.LlmArg' );
    if( !args.length )
      return;
    for( var i = 0; i < args.length; ++i )
      args[i].style.display = 'inline';
    more.style.display = 'inline';

    // The unrounded width (clientWidth is rounded, so could overstate it), less a pixel to spare, as
    //  overflowing by even a fraction of one has the CSS ellipsis cut off the ",...}".
    var elRect = el.getBoundingClientRect();
    var left = elRect.left;
    var avail = elRect.width - 1;
    var moreWidth = more.getBoundingClientRect().width;
    var closeWidth = close.getBoundingClientRect().width;
    var end = function( i ){ return args[i].getBoundingClientRect().right - left; };

    // All of them, or as many as fit followed by ",...}" - but at least the first (cut by the ellipsis).
    var nshow = args.length;
    if( end( nshow - 1 ) + closeWidth > avail ){
      nshow = 1;
      while( (nshow < args.length) && ((end( nshow ) + moreWidth + closeWidth) <= avail) )
        ++nshow;
    }

    for( var i = nshow; i < args.length; ++i )
      args[i].style.display = 'none';
    more.style.display = (nshow < args.length) ? 'inline' : 'none';
    fitWidth = header.clientWidth;
  };

  // Wt may still be laying out the page (e.g., just after the assistant moved to the side panel), so
  //  fit again if the width soon changes.
  var fitAndRecheck = function(){
    fit();
    var recheck = function(){
      if( (fitWidth >= 0) && (header.clientWidth !== fitWidth) )
        fit();
    };
    setTimeout( recheck, 200 );
    setTimeout( recheck, 1000 );
  };

  if( el.getClientRects().length ){
    fitAndRecheck();
  }else if( window.IntersectionObserver ){
    // Not displayed yet (e.g., a collapsed panel, a hidden tab, or debug options are on, which show
    //  all the arguments instead), so fit it once it first is.
    var observer = new IntersectionObserver( function( entries ){
      for( var i = 0; i < entries.length; ++i ){
        if( entries[i].isIntersecting ){
          observer.disconnect();
          fitAndRecheck();
          return;
        }
      }
    } );
    observer.observe( el );
  }
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


  /** Formats a duration like "4 ms", "13.1 s", or "6 m 14 s". */
  WString formatDuration( const chrono::milliseconds &duration )
  {
    const long long ms = duration.count();
    if( ms < 1000 )
      return WString::tr( "lcw-time-ms" ).arg( ms );

    if( ms < 60000 )
    {
      ostringstream oss;
      oss << fixed << setprecision(1) << (0.001 * ms);
      return WString::tr( "lcw-time-s" ).arg( oss.str() );
    }

    return WString::tr( "lcw-time-min-s" ).arg( ms / 60000 ).arg( (ms % 60000) / 1000 );
  }//formatDuration(...)


  /** Formats a USD cost like "Cost: $0.00423" - three significant figures below $1 (per-request
   costs are often fractions of a cent), and cents at $1 or more. */
  WString formatCost( const double usd )
  {
    int decimals = 2;
    if( (usd > 0.0) && (usd < 1.0) )
      decimals = std::min( 8, std::max( 2, 2 - static_cast<int>( std::floor( std::log10( usd ) ) ) ) );

    ostringstream oss;
    oss << fixed << setprecision( decimals ) << usd;
    return WString::tr( "lcw-summary-cost" ).arg( oss.str() );
  }//formatCost(...)


  /** How long a finished conversation took: to its finishTime, or else to its latest turn. */
  chrono::milliseconds totalDuration( const LlmInteraction &convo )
  {
    chrono::system_clock::time_point endTime = convo.timestamp;
    if( convo.finishTime.has_value() )
    {
      endTime = convo.finishTime.value();
    }else if( convo.interruption.has_value() )
    {
      endTime = std::max( endTime, convo.interruption->time );
    }else
    {
      for( const shared_ptr<LlmInteractionTurn> &turn : convo.responses )
        endTime = std::max( endTime, turn->timestamp() );
    }

    return std::max( chrono::milliseconds(0),
                     chrono::duration_cast<chrono::milliseconds>( endTime - convo.timestamp ) );
  }//totalDuration(...)


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


  /** Appends `j` to `out` for people to read: indented, but with any object or array that fits within
   about 100 characters kept on one line (so a table of numbers is not one number per line). */
  void append_readable_json( const nlohmann::json &j, const size_t indent, string &out )
  {
    const string compact = compact_json( j );
    if( !j.is_structured() || j.empty() || ((indent + compact.size()) <= 100) )
    {
      out += compact;
      return;
    }

    out += j.is_object() ? "{\n" : "[\n";
    for( nlohmann::json::const_iterator iter = j.begin(); iter != j.end(); ++iter )
    {
      out.append( indent + 2, ' ' );
      if( j.is_object() )
        out += compact_json( nlohmann::json( iter.key() ) ) + ": ";
      append_readable_json( iter.value(), indent + 2, out );
      out += (std::next( iter ) == j.end()) ? "\n" : ",\n";
    }
    out.append( indent, ' ' );
    out += j.is_object() ? "}" : "]";
  }//append_readable_json(...)


  /** Format JSON for display; see append_readable_json(). */
  string formatJson( const nlohmann::json &j )
  {
    string out;
    append_readable_json( j, 0, out );
    return out;
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


  /** Shows a tool's result in `area`: formatted for reading, if JSON, with rows for up to 15 lines. */
  void set_tool_result_text( WTextArea *area, const string &content )
  {
    const string text = formatJsonString( content );
    const int nlines = 1 + static_cast<int>( std::count( begin(text), end(text), '\n' ) );
    area->setText( WString::fromUTF8( text ) );
    area->setRows( std::min( std::max( nlines, 2 ), 15 ) );
  }//set_tool_result_text(...)


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


void LlmJsRenderedDiv::setRenderJs( const std::string &js )
{
  m_renderJs = js;
  if( isRendered() && !js.empty() )
    doJavaScript( js );
}//setRenderJs(...)


void LlmJsRenderedDiv::setStateJs( const std::string &js )
{
  m_stateJs = js;
  if( isRendered() && !js.empty() )
    doJavaScript( js );
}//setStateJs(...)


void LlmJsRenderedDiv::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  // A full render means the DOM element is being (re-)created, empty.
  if( flags.test( RenderFlag::Full ) )
  {
    if( !m_renderJs.empty() )
      doJavaScript( m_renderJs );
    if( !m_stateJs.empty() )
      doJavaScript( m_stateJs );
  }

  WContainerWidget::render( flags );
}//render(...)


LlmElapsedText::LlmElapsedText()
  : WText(),
    m_counting( false ),
    m_format(),
    m_start()
{
  setTextFormat( TextFormat::Plain );
}//LlmElapsedText constructor


void LlmElapsedText::start( const WString &label, const chrono::milliseconds alreadyElapsed )
{
  m_counting = true;
  m_start = chrono::steady_clock::now() - alreadyElapsed;

  // The counter's text, with its "{1}" left for the seconds (Wt fills arguments in order, so the
  //  "{1}" of the seconds' format survives being put in as argument 2).
  const WString secs = WString::tr( "lcw-time-s" );
  m_format = (label.empty() ? secs : WString::tr( "lcw-status-counting" ).arg( label ).arg( secs )).toUTF8();

  // What the client shows until its first tick (and what a full render starts from).
  setText( WString::fromUTF8( m_format ).arg( chrono::duration_cast<chrono::seconds>( alreadyElapsed ).count() ) );

  // A tooltip about a previous final text (e.g., before a retry) no longer applies.
  setToolTip( WString() );

  if( isRendered() )
    startClientCounter();
}//start(...)


void LlmElapsedText::setFinalText( const WString &text )
{
  // Run the stop JS even if no longer counting: it also sets the text, and Wt runs JavaScript after
  //  applying the DOM changes, so a stop issued earlier in this same response (e.g., finalizing the
  //  counter just before showing "Stopped") would otherwise overwrite this text.
  m_counting = false;
  if( isRendered() )
  {
    LOAD_JAVASCRIPT( wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmElapsedStop );
    doJavaScript( "Wt.WT.LlmElapsedStop(" + jsRef() + "," + text.jsStringLiteral() + ");" );
  }

  setText( text );
}//setFinalText(...)


void LlmElapsedText::startClientCounter()
{
  LOAD_JAVASCRIPT( wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmElapsedStart );

  const long long elapsedMs = chrono::duration_cast<chrono::milliseconds>(
                                                   chrono::steady_clock::now() - m_start ).count();
  doJavaScript( "Wt.WT.LlmElapsedStart(" + jsRef() + "," + WString::fromUTF8( m_format ).jsStringLiteral()
                + "," + to_string(elapsedMs) + ");" );
}//startClientCounter()


void LlmElapsedText::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  // A full render means the DOM element is being (re-)created, without a counter running.
  if( m_counting && flags.test( RenderFlag::Full ) )
    startClientCounter();

  WText::render( flags );
}//render(...)


namespace LlmMarkdown
{
void load_resources( WApplication *app )
{
  if( !app )
    return;

  // Order matters: the marked-katex bridge resolves the `katex` global inside its UMD factory, so
  //  KaTeX must be loaded before it.  Wt de-duplicates these by URI, so repeat calls cost nothing.
  app->useStyleSheet( "InterSpec_resources/assets/js/katex-0.18.5/katex.min.css" );
  app->require( "InterSpec_resources/assets/js/katex-0.18.5/katex.min.js" );
  app->require( "InterSpec_resources/assets/js/marked-18.0.11/marked.umd.js" );
  app->require( "InterSpec_resources/assets/js/marked-katex-extension-5.1.12/index.umd.js" );
  app->require( "InterSpec_resources/LlmMarkdown.js" );
}//load_resources(...)


WContainerWidget *add_text( WContainerWidget *parent, const string &markdown,
                            const string &style_class )
{
  assert( parent );
  if( !parent )
    return nullptr;

  load_resources( wApp );

  LlmJsRenderedDiv *div = parent->addNew<LlmJsRenderedDiv>();
  div->addStyleClass( "LlmMarkdownBody" );
  if( !style_class.empty() )
    div->addStyleClass( style_class );

  // LlmRenderMarkdown waits for LlmMarkdown.js, which in turn waits for marked/KaTeX, and either one
  //  falls back to showing the plain text - so the content can not be lost here.
  LOAD_JAVASCRIPT( wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmRenderMarkdown );
  div->setRenderJs( "Wt.WT.LlmRenderMarkdown('" + div->id() + "',"
                    + WString(markdown).jsStringLiteral() + ");" );

  return div;
}//add_text(...)


void set_raw_mode( WWidget *widget, const bool raw )
{
  assert( widget );
  if( !widget )
    return;

  // Guarded so a failed script load leaves the (already plain) text alone rather than throwing.
  const string js = "if(typeof LlmMarkdown!=='undefined')"
                    " LlmMarkdown.setRawMode('" + widget->id() + "',"
                    + string(raw ? "true" : "false") + ");";

  // Remembered on the widget (made by add_text), so a re-created element gets the same mode.
  LlmJsRenderedDiv * const div = dynamic_cast<LlmJsRenderedDiv *>( widget );
  if( div )
    div->setStateJs( js );
  else
    widget->doJavaScript( js );
}//set_raw_mode(...)
}//namespace LlmMarkdown



//
// LlmInteractionTurnDisplay implementation
//

LlmInteractionTurnDisplay::LlmInteractionTurnDisplay( shared_ptr<LlmInteractionTurn> turn,
                                                      LlmToolRoundDisplay *round )
  : WPanel(),
    m_turn( turn ),
    m_menuIcon( nullptr ),
    m_menu( nullptr ),
    m_bodyContainer( nullptr ),
    m_bodyCreated( false ),
    m_statusText( nullptr ),
    m_round( round ),
    m_debugMenuItems{},
    m_phoneLayout( is_phone_layout() )
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


LlmElapsedText *LlmInteractionTurnDisplay::statusText()
{
  if( m_round )
    return m_round->statusText();

  if( m_statusText )
    return m_statusText;

  WContainerWidget * const titleBar = titleBarWidget();
  if( !titleBar )
    return nullptr;  // e.g., the compact phone tool-results display

  // Right-aligned, just before the menu icon (if there is one).
  std::unique_ptr<LlmElapsedText> status = std::make_unique<LlmElapsedText>();
  m_statusText = status.get();
  m_statusText->addStyleClass( "LlmTurnStatus" );
  if( m_menuIcon && (titleBar->indexOf( m_menuIcon ) >= 0) )
    titleBar->insertBefore( std::move(status), m_menuIcon );
  else
    titleBar->addWidget( std::move(status) );

  return m_statusText;
}//statusText()


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


void LlmInteractionTurnDisplay::setTurnTitle( const WString &title )
{
  setTitle( title );
  mark_panel_title( this );
}//void setTurnTitle( const WString &title )


void LlmInteractionTurnDisplay::createMenuIcon()
{
  // Create menu button in title bar.  The menu itself is built lazily on first
  // click (showMenu), so turns whose menu is never opened cost nothing - which is
  // every turn during headless/benchmark runs.
  m_menuIcon = titleBarWidget()->addNew<WPushButton>();
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn" );
  if( !menuAlwaysShown() )
    m_menuIcon->addStyleClass( "LlmDebugOnly" );
  m_menuIcon->clicked().preventPropagation();
  m_menuIcon->clicked().connect( this, &LlmInteractionTurnDisplay::showMenu );
}//createMenuIcon()


bool LlmInteractionTurnDisplay::menuAlwaysShown() const
{
  return false;
}//bool menuAlwaysShown() const


void LlmInteractionTurnDisplay::updateMenuItemsShown()
{
  const bool showDebug = show_debug_options();
  for( WMenuItem * const item : m_debugMenuItems )
    item->setHidden( !showDebug );
}//void updateMenuItemsShown()


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

  updateMenuItemsShown();

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

  m_debugMenuItems.push_back( showRequestJson );
  m_debugMenuItems.push_back( showResponseJson );
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

  WPushButton *copyBtn = dialog->addButton( "Copy to Clipboard", WidgetUtils::ButtonRole::Neutral );
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

  WPushButton *closeBtn = dialog->addButton( "Close", WidgetUtils::ButtonRole::Dismiss );
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
    m_response( response ),
    m_contentDiv( nullptr ),
    m_showRawText( false ),
    m_rawTextItem( nullptr )
{
  addStyleClass( "LlmFinalResponseDisplay" );

  // Set title
  setTurnTitle( getTitleText() );

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


bool LlmInteractionFinalResponseDisplay::menuAlwaysShown() const
{
  return true;  // "Copy Response Text" is for everyone
}//bool menuAlwaysShown() const


void LlmInteractionFinalResponseDisplay::addMenuItems( PopupDivMenu *menu )
{
  // Base items (Request/Response JSON) first, then a convenience "Copy Response Text" that puts the
  // human-readable assistant text (not the raw JSON) on the clipboard.  Phones get only the latter,
  // and only "Copy Response Text" is shown on desktop unless debug options are on.
  if( m_phoneLayout )
    menu->addPhoneBackItem( nullptr );
  else
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

  if( m_phoneLayout )
    return;

  // Escape hatch to see the Markdown/LaTeX the model actually wrote (and for when rendering it
  //  obscures something).  "Copy Response Text" above always copies the raw source, either way.
  PopupDivMenuItem *rawText = menu->addMenuItem( "Show Raw Text" );
  rawText->triggered().connect( std::bind( [this,rawText](){
    if( !m_contentDiv )
      return;

    m_showRawText = !m_showRawText;
    LlmMarkdown::set_raw_mode( m_contentDiv, m_showRawText );
    rawText->setText( m_showRawText ? "Show Rendered Text" : "Show Raw Text" );
  }));
  m_debugMenuItems.push_back( rawText );
  m_rawTextItem = rawText;
}//addMenuItems()


void LlmInteractionFinalResponseDisplay::updateMenuItemsShown()
{
  LlmInteractionTurnDisplay::updateMenuItemsShown();

  // Showing the raw text must stay undoable, even if debug options were turned off since.
  if( m_rawTextItem && m_showRawText )
    m_rawTextItem->setHidden( false );
}//void updateMenuItemsShown()


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
  // Display the full response content, rendered as Markdown + LaTeX client-side.  The raw model text
  // deliberately never goes through Wt's XHTML parser (a bare '&' or '<' makes it log
  // "Error reading XHTML string" before silently falling back); LlmMarkdown.js escapes everything it
  // does not explicitly allow.
  m_contentDiv = LlmMarkdown::add_text( m_bodyContainer, m_response->content(), "LlmResponseContent" );

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
    const WString durationStr = formatDuration( m_response->callDuration().value() );
    WText *durationText = metadataDiv->addNew<WText>( WString::tr( "lcw-response-duration" ).arg( durationStr ) );
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

    if( convo->cost.has_value() )
    {
      WText *costText = metadataDiv->addNew<WText>( " | " + formatCost( convo->cost.value() ) );
      costText->addStyleClass( "LlmTokenUsage" );
    }
  }

  // Thinking content (if available)
  if( !m_response->thinkingContent().empty() )
  {
    WContainerWidget *thinkingDiv = m_bodyContainer->addNew<WContainerWidget>();
    thinkingDiv->addStyleClass( "LlmThinkingContent" );

    WText *thinkingLabel = thinkingDiv->addNew<WText>( "Thinking:" );
    thinkingLabel->addStyleClass( "LlmThinkingLabel" );

    // Thinking content is also Markdown/LaTeX from the model, so render it the same way.
    LlmMarkdown::add_text( thinkingDiv, m_response->thinkingContent(), "" );
  }
}//createBodyContent()


//
// LlmToolRequestDisplay implementation
//

LlmToolRequestDisplay::LlmToolRequestDisplay( shared_ptr<LlmToolRequest> request,
                                              LlmToolRoundDisplay *round )
  : LlmInteractionTurnDisplay( request, round ),
    m_request( request )
{
  addStyleClass( "LlmToolRequestDisplay" );

  // In a tool round, the results display shows each tools parameters (with debug options on), so
  //  this panel is just for its menu of the raw JSON, and only shown with debug options on.
  if( m_round )
    addStyleClass( "LlmDebugOnly" );

  // Set title
  setTurnTitle( getTitleText() );

  // Collapsed by default
  setCollapsible( true );
  setCollapsed( true );

  // Create title bar with menu icon (we have to do this AFTER calling setTitle, or the button will disapear)
  if( !m_phoneLayout )
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

#if( !NO_TOOL_CALL_ID_TEXT )
    // Invocation ID - TextFormat::Plain so a stray '&'/'<' in an id can't hit the XHTML parser.
    WText *idText = callDiv->addNew<WText>( " (ID: " + call.invocationId + ")", TextFormat::Plain );
    idText->addStyleClass( "LlmInvocationId LlmDebugOnly" );
#endif

    // The parameters are shown by the tool round's results display, on each call's header line.
  }
}//createBodyContent()


//
// LlmToolResultsDisplay implementation
//

LlmToolResultsDisplay::LlmToolResultsDisplay( shared_ptr<LlmToolResults> results,
                                              std::weak_ptr<LlmInterface> llmInterface,
                                              int nestingLevel,
                                              WContainerWidget *phoneVisuals,
                                              LlmToolRoundDisplay *round )
  : LlmInteractionTurnDisplay( results, round ),
    m_results( results ),
    m_llmInterface( llmInterface ),
    m_nestingLevel( nestingLevel ),
    m_phoneVisuals( phoneVisuals )
{
  addStyleClass( "LlmToolResultsDisplay" );

  if( m_phoneVisuals )
  {
    // Compact phone entry in the interactions tool-summary row: no title bar (so no collapsing or
    //  menu), just the per-tool lines.
    addStyleClass( "LlmPhoneCompact" );
    createBodyContent();
    m_bodyCreated = true;
    return;
  }//if( m_phoneVisuals )

  // Set title
  setTurnTitle( getTitleText() );

  if( m_round )
  {
    // The round panel does the collapsing (and has the title the user sees), so our title bar is
    //  only for the debug menu; we must not be collapsible, or our body could be hidden while the
    //  title bar to expand it is.
    titleBarWidget()->addStyleClass( "LlmDebugOnly" );
  }else
  {
    // Collapsed by default, but expand if any tool call has visual content
    setCollapsible( true );
    setCollapsed( !wantsExpanded() );
  }//if( m_round ) / else

  // Create title bar with menu icon (we have to do this AFTER calling setTitle, or the button will disapear)
  createMenuIcon();
  
  // Create content immediately (lazy loading to be added later)
  createBodyContent();
  m_bodyCreated = true;
}//LlmToolResultsDisplay constructor


LlmToolResultsDisplay::~LlmToolResultsDisplay()
{
}//~LlmToolResultsDisplay()


bool LlmToolResultsDisplay::wantsExpanded() const
{
  for( const LlmToolCall &call : m_results->toolCalls() )
  {
    if( !call.visualContent.empty() || call.imageContent.has_value() || call.sub_agent_conversation )
      return true;
  }
  return false;
}//bool wantsExpanded() const


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
  // Error styling (not for a call the user cancelled)
  resultDiv->toggleStyleClass( "LlmToolResultError", is_failed( call ) );

  // On desktop, a lone sub-agents conversation follows, in a panel already titled "X agent", so its
  //  tool name and duration line is only for debugging - unless it failed (see also
  //  handleSubAgentFinished()), or the round has other calls it must be told apart from.
  const bool debugOnlyHeader = !m_phoneVisuals && call.sub_agent_conversation && !is_failed( call )
                               && (m_results->toolCalls().size() == 1);

  // The one line of tool name, arguments, status, and duration
  LlmJsRenderedDiv * const header = resultDiv->addNew<LlmJsRenderedDiv>();
  header->addStyleClass( "LlmToolCallHeader" );
  if( debugOnlyHeader )
    header->addStyleClass( "LlmDebugOnly" );
  if( call.sub_agent_conversation )
    m_subAgentHeaders[call.sub_agent_conversation] = header;
  renderCallHeader( call, header );

  // Result content - show if not pending, available, and (no sub-agent OR sub-agent completed).
  //  The raw result is for developing the tools, so phones skip it, and it is debug-only otherwise.
  if( !m_phoneVisuals
      && !call.content.empty()
      && call.status != LlmToolCall::CallStatus::Pending
      && (!call.sub_agent_conversation || call.sub_agent_conversation->finishTime.has_value()) )
  {
    WContainerWidget * const rawResultDiv = resultDiv->addNew<WContainerWidget>();
    rawResultDiv->addStyleClass( "LlmDebugOnly" );

    WText * const resultLabel = rawResultDiv->addNew<WText>( "Result:" );
    resultLabel->addStyleClass( "LlmResultLabel" );

    WTextArea * const resultArea = rawResultDiv->addNew<WTextArea>();
    resultArea->setReadOnly( true );
    set_tool_result_text( resultArea, call.content );
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
    renderVisualContent( visual, m_phoneVisuals ? m_phoneVisuals : resultDiv );
}//renderToolCallResult(...)


void LlmToolResultsDisplay::renderCallHeader( const LlmToolCall &call, LlmJsRenderedDiv *header )
{
  // The duration comes first, so that with debug options on, where the header wraps, it can float to
  //  the right of the first line; otherwise CSS (flex `order`) shows it after the arguments and status.
  if( call.executionDuration.has_value() )
  {
    const WString durationStr = formatDuration( call.executionDuration.value() );
    WText * const durationText = header->addNew<WText>( WString::tr( "lcw-tool-duration" ).arg( durationStr ) );
    durationText->addStyleClass( "LlmDuration" );
  }

  // The model-derived name is escaped, but the intentional <b> markup kept.
  WText * const nameText = header->addNew<WText>( "<b>" + Wt::Utils::htmlEncode( call.toolName ) + "</b>",
                                                  TextFormat::UnsafeXHTML );
  nameText->addStyleClass( "LlmToolName" );

  // The arguments: normally on one line, as many as fit (client-side) with all of them in the
  //  tooltip, and with debug options on, all of them in full, wrapping as needed.  Phones skip them,
  //  and a sub-agents arguments are the task it was given, which its conversation shows as its question.
  const nlohmann::json callArgs = (m_phoneVisuals || call.sub_agent_conversation) ? nlohmann::json()
                                                                                  : tool_call_args( call );

  // The tool's own required arguments are generally the most telling; see leading_arg_names().
  const shared_ptr<LlmInterface> llm = m_llmInterface.lock();
  const shared_ptr<const LlmTools::ToolRegistry> registry = llm ? llm->toolRegistry() : nullptr;
  const LlmTools::SharedTool * const tool = registry ? registry->getTool( call.toolName ) : nullptr;
  const vector<pair<string, const nlohmann::json *>> args
                = ordered_args( callArgs, leading_arg_names( call.toolName, tool ? tool->parameters_schema : nlohmann::json() ) );

  const string argsHtml = tool_args_html( callArgs, args );
  if( !argsHtml.empty() )
  {
    WText * const argsText = header->addNew<WText>( argsHtml, TextFormat::UnsafeXHTML );
    argsText->addStyleClass( "LlmToolArgs" );

    // One argument per line
    string tooltip = args.empty() ? compact_json( callArgs ) : string();
    for( const pair<string, const nlohmann::json *> &arg : args )
      tooltip += (tooltip.empty() ? "" : "\n") + arg.first + ": " + compact_json( *arg.second );
    if( tooltip.size() > 2000 )
    {
      SpecUtils::utf8_limit_str_size( tooltip, 2000 );
      tooltip += "...";
    }
    argsText->setToolTip( WString::fromUTF8( tooltip ), TextFormat::Plain );

    WText * const fullArgsText = header->addNew<WText>( WString::fromUTF8( tool_args_full( callArgs, args ) ),
                                                        TextFormat::Plain );
    fullArgsText->addStyleClass( "LlmToolArgsFull LlmDebugOnly" );

    // Run again whenever the DOM element is re-created, e.g., moving the assistant to a window, where
    //  the width will differ.
    LOAD_JAVASCRIPT( wApp, "LlmInteractionDisplay.cpp", "LlmInteractionDisplay", wtjsLlmFitToolArgs );
    header->setRenderJs( "Wt.WT.LlmFitToolArgs(" + header->jsRef() + ");" );
  }//if( !argsHtml.empty() )

  const string status = tool_status_html( call );
  if( !status.empty() )
    header->addNew<WText>( status, TextFormat::UnsafeXHTML );

#if( !NO_TOOL_CALL_ID_TEXT )
  // Invocation ID - TextFormat::Plain so a stray '&'/'<' in an id can't hit the XHTML parser.
  if( !m_phoneVisuals )
  {
    WText * const idText = header->addNew<WText>( "(ID: " + call.invocationId + ")", TextFormat::Plain );
    idText->addStyleClass( "LlmInvocationId LlmDebugOnly" );
  }
#endif
}//renderCallHeader(...)


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

  // Charts are drawn client-side, so must be redrawn if the element is re-created.
  LlmJsRenderedDiv *chartDiv = vizDiv->addNew<LlmJsRenderedDiv>();
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
      chartDiv->setRenderJs( js );
      break;
    }

    case LlmToolCall::VisualContent::VisualType::ShieldingDiagram2D:
    {
      wApp->require( "InterSpec_resources/Shielding2DView.js" );
      const string chartId = chartDiv->id();
      // Redrawn on resize (e.g., when its element is re-created before it has its final width)
      chartDiv->setRenderJs(
        "var el=document.getElementById('" + chartId + "');"
        "if(el){"
          "el._chartObj=new Shielding2DView('" + chartId + "', " + visual.data + ");"
          "new ResizeObserver(function(){"
            "if(el._chartObj) el._chartObj.handleResize();"
          "}).observe(el);"
        "}" );
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
        if( !window->isPhone() )  // full-screen on phones
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

        // Size to 80% of InterSpec rendered area (it is full-screen on phones)
        InterSpec *viewer = InterSpec::instance();
        const int dlgWidth = viewer ? static_cast<int>( 0.8 * viewer->renderedWidth() ) : 700;
        const int dlgHeight = viewer ? static_cast<int>( 0.8 * viewer->renderedHeight() ) : 500;
        if( !window->isPhone() )
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
      // On desktop, the sub-agent is in a round panel titled "X agent", so needs no label.
      if( m_phoneVisuals )
      {
        WText * const subAgentLabel = resultDiv->addNew<WText>( "Sub-agent conversation:" );
        subAgentLabel->addStyleClass( "LlmSubAgentLabel" );
      }

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

        // Create nested LlmInteractionDisplay for sub-agent (it styles itself as nested)
        LlmInteractionDisplay *subAgentDisplay = resultDiv->addNew<LlmInteractionDisplay>(
                    call.sub_agent_conversation, m_llmInterface, m_nestingLevel + 1, m_phoneVisuals );

        // While it runs, trim the margins below it, so its newest turn isn't scrolled out of view.
        if( !call.sub_agent_conversation->finishTime.has_value() )
          resultDiv->addStyleClass( "LlmSubAgentRunning" );

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

      // Make sure we expand the parent widget so we can see the sub-agent conversation going on (a
      //  round expands itself, once we are added to it).
      if( !m_round && isCollapsed() )
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
  m_debugMenuItems.push_back( showToolResultsJson );
}//addMenuItems()


void LlmToolResultsDisplay::handleSubAgentFinished( shared_ptr<LlmInteractionTurn> turn,
                                                    shared_ptr<LlmInteraction> subAgent )
{
  // However the sub-agent ended, it is no longer running, so restore its normal spacing.
  const auto runningDivIter = m_subAgentResultDivs.find( subAgent );
  if( runningDivIter != m_subAgentResultDivs.end() )
    runningDivIter->second->removeStyleClass( "LlmSubAgentRunning" );

  // The sub-agent may have failed, which the rounds title shows.
  if( m_round )
    m_round->updateTitle();

  // If it failed, mark its result item, and always show its header line, with the error marker (it
  //  was rendered while the sub-agent was running, so may be debug-only).
  for( const LlmToolCall &call : m_results->toolCalls() )
  {
    if( (call.sub_agent_conversation != subAgent) || !is_failed( call ) )
      continue;

    if( runningDivIter != m_subAgentResultDivs.end() )
      runningDivIter->second->addStyleClass( "LlmToolResultError" );

    const auto headerIter = m_subAgentHeaders.find( subAgent );
    if( headerIter != m_subAgentHeaders.end() )
    {
      LlmJsRenderedDiv * const header = headerIter->second;
      header->clear();
      renderCallHeader( call, header );
      header->removeStyleClass( "LlmDebugOnly" );
    }
    break;
  }//for( loop over calls to find a failed sub-agent )

  // The compact phone display does not show raw tool results.
  if( m_phoneVisuals )
    return;

  // Sub-agent finished - create the result text area if content is available
  // Find the corresponding tool call in m_results.  A sub-agent ended early (cancel, timeout) still
  //  has its "processing" placeholder, so skip it while it is Pending.
  for( const LlmToolCall &call : m_results->toolCalls() )
  {
    if( call.sub_agent_conversation == subAgent && !call.content.empty()
        && (call.status != LlmToolCall::CallStatus::Pending) )
    {
      // Find the result div for this sub-agent
      const auto divIter = m_subAgentResultDivs.find( subAgent );
      if( divIter == m_subAgentResultDivs.end() )
      {
        cerr << "LlmToolResultsDisplay::handleSubAgentFinished(): Could not find result div for sub-agent" << endl;
        return;
      }

      WContainerWidget *resultDiv = divIter->second;

      // Create the result label and text area; only for debugging, as the sub-agents final response
      //  above is what it returned.
      WContainerWidget * const rawResultDiv = resultDiv->addNew<WContainerWidget>();
      rawResultDiv->addStyleClass( "LlmDebugOnly" );

      WText * const resultLabel = rawResultDiv->addNew<WText>( "Result:" );
      resultLabel->addStyleClass( "LlmResultLabel" );

      WTextArea * const resultArea = rawResultDiv->addNew<WTextArea>();
      resultArea->setReadOnly( true );
      set_tool_result_text( resultArea, call.content );
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
    if( m_round )
    {
      m_round->updateTitle();  // it may have failed
      if( !call.visualContent.empty() )
        m_round->expandIfWanted();
    }else if( !call.visualContent.empty() && isCollapsed() )
    {
      setCollapsed( false );
    }

    // Remove from tracking map since it's now complete
    m_asyncToolResultDivs.erase( divIter );
    break;
  }
}//handleAsyncToolCompleted(...)


//
// LlmToolRoundDisplay implementation
//

LlmToolRoundDisplay::LlmToolRoundDisplay( shared_ptr<LlmToolRequest> request )
  : WPanel(),
    m_request( request ),
    m_results( nullptr ),
    m_resultsDisplay( nullptr ),
    m_body( nullptr ),
    m_status( nullptr )
{
  // Styled as a turn panel, which it stands in for.
  addStyleClass( "LlmInteractionTurnDisplay LlmToolRoundDisplay" );

  // The title must be set before the status text is added, so it comes first in the title bar.
  updateTitle();
  setCollapsible( true );

  m_status = titleBarWidget()->addNew<LlmElapsedText>();
  m_status->addStyleClass( "LlmTurnStatus" );

  m_body = setCentralWidget( std::make_unique<WContainerWidget>() );
  m_body->addStyleClass( "LlmToolRoundBody" );

  setCollapsed( true );  // Until results with something to see are added (see expandIfWanted())
}//LlmToolRoundDisplay constructor


LlmToolRoundDisplay::~LlmToolRoundDisplay()
{
}//~LlmToolRoundDisplay()


LlmToolRequestDisplay *LlmToolRoundDisplay::addRequestDisplay()
{
  if( !m_request )
    return nullptr;
  return m_body->addNew<LlmToolRequestDisplay>( m_request, this );
}//LlmToolRequestDisplay *addRequestDisplay()


LlmToolResultsDisplay *LlmToolRoundDisplay::addResultsDisplay( shared_ptr<LlmToolResults> results,
                                                               std::weak_ptr<LlmInterface> llmInterface,
                                                               int nestingLevel )
{
  assert( !m_results && !m_resultsDisplay );

  m_results = results;
  m_resultsDisplay = m_body->addNew<LlmToolResultsDisplay>( results, llmInterface, nestingLevel,
                                                            nullptr, this );
  updateTitle();
  expandIfWanted();

  return m_resultsDisplay;
}//LlmToolResultsDisplay *addResultsDisplay(...)


LlmElapsedText *LlmToolRoundDisplay::statusText()
{
  return m_status;
}//LlmElapsedText *statusText()


bool LlmToolRoundDisplay::hasTurn( const shared_ptr<LlmInteractionTurn> &turn ) const
{
  return turn && ((turn == m_request) || (turn == m_results));
}//bool hasTurn(...)


void LlmToolRoundDisplay::showFinalDurations()
{
  // The tools time is the request turns, and the LLMs time (to respond to the results) the results
  //  turns; either may be missing (e.g., the round is still going, or the LLM call failed).
  const std::optional<chrono::milliseconds> tools = m_request ? m_request->callDuration() : std::nullopt;
  const std::optional<chrono::milliseconds> llm = m_results ? m_results->callDuration() : std::nullopt;

  WString text, tooltip;
  if( tools.has_value() && llm.has_value() )
  {
    text = WString::tr( "lcw-duration-round" ).arg( formatDuration( llm.value() ) )
                                               .arg( formatDuration( tools.value() ) );
    tooltip = WString::tr( "lcw-duration-round-tt" );
  }else if( llm.has_value() )
  {
    text = WString::tr( "lcw-duration-llm" ).arg( formatDuration( llm.value() ) );
    tooltip = WString::tr( "lcw-duration-llm-tt" );
  }else if( tools.has_value() )
  {
    text = WString::tr( "lcw-duration-tools" ).arg( formatDuration( tools.value() ) );
    tooltip = WString::tr( "lcw-duration-tools-tt" );
  }

  m_status->setFinalText( text );
  m_status->setToolTip( tooltip );
}//void showFinalDurations()


void LlmToolRoundDisplay::updateTitle()
{
  // The tools are listed from the request, if we have it, as that is what the LLM asked for.
  const vector<LlmToolCall> * const calls = m_request ? &m_request->toolCalls()
                                        : (m_results ? &m_results->toolCalls() : nullptr);

  bool allAgents = true;
  string names;
  for( size_t i = 0; calls && (i < calls->size()); ++i )
  {
    bool isAgent = false;
    names += (i ? ", " : "") + tool_display_name( (*calls)[i].toolName, isAgent );
    allAgents = (allAgents && isAgent);
  }

  bool hasErrors = false, hasCancelled = false;
  for( size_t i = 0; m_results && (i < m_results->toolCalls().size()); ++i )
  {
    hasErrors = (hasErrors || is_failed( m_results->toolCalls()[i] ));
    hasCancelled = (hasCancelled || is_cancelled( m_results->toolCalls()[i] ));
  }

  const bool single = (calls && (calls->size() == 1));
  const char *key = nullptr;  // null for just the agent names
  if( names.empty() )
    key = "lcw-round-no-tools";
  else if( allAgents && hasErrors )
    key = single ? "lcw-round-agent-failed" : "lcw-round-agents-errors";
  else if( allAgents && hasCancelled )
    key = single ? "lcw-round-agent-cancelled" : "lcw-round-agents-cancelled";
  else if( !allAgents )
    key = hasErrors ? "lcw-round-tools-errors" : (hasCancelled ? "lcw-round-tools-cancelled" : "lcw-round-tools");

  const WString title = key ? WString::tr( key ).arg( WString::fromUTF8( names ) ) : WString::fromUTF8( names );

  // Escape model-derived tool names before they reach the XHTMLText panel title.
  const string truncated = truncateString( title.toUTF8(), 100 );
  setTitle( WString::fromUTF8( Wt::Utils::htmlEncode( truncated ) ) );
  mark_panel_title( this );
  toggleStyleClass( "LlmToolRoundError", hasErrors );
}//void updateTitle()


void LlmToolRoundDisplay::expandIfWanted()
{
  if( m_resultsDisplay && m_resultsDisplay->wantsExpanded() && isCollapsed() )
    setCollapsed( false );
}//void expandIfWanted()


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
  setTurnTitle( getTitleText() );

  // NOT collapsed by default
  setCollapsible( true );
  setCollapsed( false );

  // Create title bar with menu icon (we have to do this AFTER calling setTitle, or the button will
  //  disapear).  Phones dont get the menu, as it only has the raw request/response JSON.
  if( !m_phoneLayout )
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

  // A sub-agent's error pauses it (and the agent that invoked it), so it can be continued, or
  //  cancelled to have the invoking agent carry on without it; for the main agent, the failed step
  //  can be skipped.
  const shared_ptr<LlmInteraction> convo = m_error->conversation().lock();
  const bool subAgent = convo && (convo->agent_type != AgentType::MainAgent);

  m_retryBtn = m_bodyContainer->addNew<WPushButton>( WString::tr( "lcw-error-retry" ) );
  m_retryBtn->addStyleClass( "LlmRetryButton" );
  m_retryBtn->clicked().connect( this, &LlmInteractionErrorDisplay::handleRetry );

  m_continueBtn = m_bodyContainer->addNew<WPushButton>(
                                 WString::tr( subAgent ? "lcw-error-cancel-agent" : "lcw-error-skip" ) );
  m_continueBtn->addStyleClass( "LlmContinueButton" );
  m_continueBtn->clicked().connect( this, &LlmInteractionErrorDisplay::handleContinueAnyway );

  m_helpText = m_bodyContainer->addNew<WText>( WString::tr( subAgent ? "lcw-error-help-agent"
                                                                     : "lcw-error-help-main" ) );
  m_helpText->addStyleClass( "LlmErrorHelpText" );

  // If the sub-agent is ended before anything is chosen here (e.g., a new question is asked), say so.
  if( subAgent )
  {
    convo->conversationFinished.connect( this, [this](){
      if( m_retryBtn && m_retryBtn->isEnabled() && m_continueBtn && m_continueBtn->isEnabled() )
        refuseIfEndedSubAgent();
    } );
    if( convo->finishTime.has_value() )
      refuseIfEndedSubAgent();
  }//if( subAgent )
}//createBodyContent()


bool LlmInteractionErrorDisplay::refuseIfEndedSubAgent()
{
  const shared_ptr<LlmInteraction> convo = m_error->conversation().lock();
  if( !convo || (convo->agent_type == AgentType::MainAgent) || !convo->finishTime.has_value() )
    return false;

  if( m_retryBtn )
    m_retryBtn->setEnabled( false );
  if( m_continueBtn )
    m_continueBtn->setEnabled( false );
  if( m_helpText )
    m_helpText->setText( WString::tr( "lcw-subagent-ended" ) );

  return true;
}//refuseIfEndedSubAgent()


bool LlmInteractionErrorDisplay::refuseIfSubAgentNotAwaited()
{
  const shared_ptr<LlmInteraction> convo = m_error->conversation().lock();
  const shared_ptr<LlmInterface> llm = m_llmInterface.lock();
  if( !convo || !llm || (convo->agent_type == AgentType::MainAgent) || llm->resumeInterruptedFor( convo ) )
    return false;

  if( m_retryBtn )
    m_retryBtn->setEnabled( false );
  if( m_continueBtn )
    m_continueBtn->setEnabled( false );
  if( m_helpText )
    m_helpText->setText( WString::tr( "lcw-subagent-ended" ) );

  return true;
}//refuseIfSubAgentNotAwaited()


void LlmInteractionErrorDisplay::handleRetry()
{
  if( refuseIfEndedSubAgent() || refuseIfSubAgentNotAwaited() )
    return;

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

  // The conversation is running again (e.g., after a watchdog timeout, which set this), so that a later
  //  cancel or timeout ends it properly.
  const std::optional<std::chrono::system_clock::time_point> previousFinishTime = convo->finishTime;

  // For a sub-agent, the tool rounds waiting on it (stopped while it was paused) go back in flight.
  const bool resumedRounds = llmInterface->resumeStoppedRounds( convo );

  try
  {
    // Remove the error from conversation history
    // This allows the conversation to continue as if the error never happened
    convo->responses.pop_back();
    convo->finishTime.reset();

    // Build API request with the conversation history (which now excludes the error)
    const nlohmann::json requestJson = llmInterface->buildMessagesArray( convo );

    // Make tracked API call
    std::pair<int,std::string> request_id_content =
      llmInterface->makeTrackedApiCall( requestJson, convo );

    cout << "Retry request sent with ID: " << request_id_content.first << endl;

    if( m_helpText )
      m_helpText->setText( WString::tr( "lcw-error-retried" ) );

    // The response will be handled by the normal LlmInterface::handleApiResponse flow
    // which will add new response turns to the conversation and emit signals
  }
  catch( const std::exception &e )
  {
    cerr << "Error during retry: " << e.what() << endl;

    // Re-add the error to the conversation since retry failed
    convo->responses.push_back( m_error );
    convo->finishTime = previousFinishTime;
    if( resumedRounds )
      llmInterface->stopIdleRounds( true );

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
  if( refuseIfEndedSubAgent() || refuseIfSubAgentNotAwaited() )
    return;

  // Disable both buttons to prevent double-clicks
  if( m_continueBtn )
    m_continueBtn->setEnabled( false );
  if( m_retryBtn )
  {
    m_retryBtn->setEnabled( false );
    m_retryBtn->hide();  // Hide the button that wasn't clicked
  }

  // Mark the error to exclude from history sent to LLM, and that the user moved on from it (for a
  //  sub-agent, the agent that invoked it is then told the user cancelled it - see the sub-agent
  //  completion handler).
  m_error->setExcludeFromHistory( true );
  m_error->setSkippedByUser( true );

  std::shared_ptr<LlmInteraction> convo = m_error->conversation().lock();
  const bool subAgent = convo && (convo->agent_type != AgentType::MainAgent);

  // Update help text to indicate continue action was taken
  if( m_helpText )
    m_helpText->setText( WString::tr( subAgent ? "lcw-error-agent-cancelled" : "lcw-error-skipped" ) );

  // For a sub-agent, the tool rounds waiting on it (stopped while it was paused) go back in flight,
  //  so ending it, below, lets the agent that invoked it continue.
  if( subAgent )
  {
    if( const std::shared_ptr<LlmInterface> llmInterface = m_llmInterface.lock() )
      llmInterface->resumeStoppedRounds( convo );
  }

  // Find and mark the request that led to this error to also exclude from history
  // Search backwards from the error to find the closest request we sent to the LLM:
  // - InitialRequest (user/system message)
  // - AutoReply (automatic continuation)
  // - ToolResults (results from tool execution)
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

    // Mark conversation as finished.  The completion handler (a sub-agent's resolves its parent's
    //  tool round) is cleared before it runs, so nothing can run it a second time.
    convo->finishTime = std::chrono::system_clock::now();
    std::function<void(void)> handler;
    handler.swap( convo->conversation_completion_handler );
    if( handler )
      handler();
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

LlmInteractionDisplay::LlmInteractionDisplay( shared_ptr<LlmInteraction> interaction,
                                              std::weak_ptr<LlmInterface> llmInterface,
                                              int nestingLevel,
                                              WContainerWidget *phoneVisuals )
  : WPanel(),
    m_interaction( interaction ),
    m_llmInterface( llmInterface ),
    m_turnContainer( nullptr ),
    m_statusText( nullptr ),
    m_timerText( nullptr ),
    m_menuIcon( nullptr ),
    m_menu( nullptr ),
    m_nestingLevel( nestingLevel ),
    m_isFinished( false ),
    m_interruptedPanel( nullptr ),
    m_isStopped( false ),
    m_isLive( false ),
    m_activeStatus( nullptr ),
    m_activeStatusTurn( nullptr ),
    m_waitingOnTools( false ),
    m_openToolRound( nullptr ),
    m_toolRounds{},
    m_lastRequestTurn( nullptr ),
    m_lastRequestDisplay( nullptr ),
    m_phoneLayout( is_phone_layout() ),
    m_phoneToolSummary( nullptr ),
    m_phoneToolList( nullptr ),
    m_phoneVisuals( m_phoneLayout ? phoneVisuals : nullptr ),
    m_phoneOpenToolRequests( 0 ),
    m_phoneToolErrors( false ),
    m_phoneStatus( nullptr )
{
  assert( m_interaction );

  addStyleClass( "LlmInteractionDisplay" );

  // On desktop, a sub-agents conversation is shown in its parents tool-round panel (titled
  //  "X agent"), which does the collapsing; so it is not shown as a panel of its own, and its title
  //  bar (with the menu) is only shown with debug options on.  It can not be collapsed, or it could
  //  be hidden while its title bar is.
  const bool flatSubAgent = (m_nestingLevel > 0) && !m_phoneLayout;
  if( flatSubAgent )
  {
    setTitleBar( true );
  }else
  {
    // Set collapsible to create title bar
    setCollapsible( true );
    setCollapsed( false );
  }//if( flatSubAgent ) / else

  // Wt4's WPanel only applies the themes title-bar style (".titlebar") from setTitle(), which we
  //  dont call - our title is made of our own status and timer texts - so apply it ourselves.
  scheduleThemeStyleApply( wApp->theme(), titleBarWidget(), Wt::PanelTitleBar );
  if( flatSubAgent )
    titleBarWidget()->addStyleClass( "LlmDebugOnly" );

  // Set title with status - create status text widget
  m_statusText = titleBarWidget()->addNew<WText>();
  m_statusText->setTextFormat( TextFormat::Plain );  // Plain text only
  m_statusText->addStyleClass( "LlmStatusIndicator" );

  // Separate widget for the total-time counter
  m_timerText = titleBarWidget()->addNew<LlmElapsedText>();
  m_timerText->setText( WString::tr( "lcw-time-s" ).arg( 0 ) );
  m_timerText->addStyleClass( "LlmStatusIndicator" );

  updateStatus();

  // Create menu icon in title bar; phones dont get it, as it only has developer-oriented items.
  if( !m_phoneLayout )
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

      // A sub-agents "question" is the whole task prompt it was given, so CSS limits it to a few
      //  lines; clicking it shows it all, or limits it again.
      if( flatSubAgent )
      {
        questionText->addStyleClass( "LlmSubAgentPrompt" );
        questionText->clicked().connect( "function(o,e){ o.classList.toggle('LlmExpanded'); }" );
      }
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
  m_turnContainer->addStyleClass( "LlmTurnContainer" );

  // Add summary section after the turns; it is empty, so hidden, until updateSummary() fills it - so
  //  it doesn't add blank space below the newest turn while the conversation is running.
  m_summaryDiv = bodyContainer->addNew<WContainerWidget>();
  m_summaryDiv->addStyleClass( "LlmInteractionSummary" );
  m_summaryDiv->setHidden( true );

  // Add any existing responses
  for( const shared_ptr<LlmInteractionTurn> &turn : interaction->responses )
  {
    handleResponseAdded( turn );
  }
  m_isLive = true;

  // Connect to signals
  interaction->responseAdded.connect( this, &LlmInteractionDisplay::handleResponseAdded );
  interaction->subAgentFinished.connect( this, &LlmInteractionDisplay::handleSubAgentFinished );
  interaction->conversationFinished.connect( this, &LlmInteractionDisplay::handleConversationFinished );
  interaction->llmRequestSent.connect( this, &LlmInteractionDisplay::handleRequestSent );
  interaction->stoppedChanged.connect( this, &LlmInteractionDisplay::handleStoppedChanged );
  interaction->interruptionChanged.connect( this, &LlmInteractionDisplay::handleInterruptionChanged );

  // Check if conversation is already finished; a set finishTime with no final or error turn means it
  //  was ended without one (e.g., a cancelled sub-agent, or the interface was closed mid-request).
  bool hasError = false;
  bool hasFinalResponse = false;
  for( const shared_ptr<LlmInteractionTurn> &turn : interaction->responses )
  {
    if( turn->type() == LlmInteractionTurn::Type::Error )
      hasError = true;
    if( turn->type() == LlmInteractionTurn::Type::FinalLlmResponse )
      hasFinalResponse = true;
  }

  if( hasError || hasFinalResponse || m_interaction->finishTime.has_value() )
  {
    m_isFinished = true;
    updateSummary();
    updateStatus();
  }else
  {
    // Count up the total time since the conversation started - but only if it has work in flight; a
    //  conversation restored unfinished has nothing that would end it, so it just shows "0s".
    const shared_ptr<LlmInterface> llm = m_llmInterface.lock();
    if( llm && llm->isConversationActive( m_interaction ) )
    {
      const chrono::system_clock::duration sinceStart = chrono::system_clock::now() - m_interaction->timestamp;
      m_timerText->start( WString(), std::max( chrono::milliseconds(0),
                                    chrono::duration_cast<chrono::milliseconds>( sinceStart ) ) );

      // The turns were replayed without live counters, so start the one for what is in flight now
      //  (e.g., this display was rebuilt from a restored history).
      resumeActiveStatus( *llm );
    }else if( llm && llm->isConversationStopped( m_interaction ) )
    {
      handleStoppedChanged( true );
    }else if( m_interaction->interruption.has_value() )
    {
      // Cut off by the assistant, or InterSpec, closing, or another spectrum being shown; show how long
      //  it ran, and offer to resume it.
      stopStatusTimer();
      showInterruptedPanel();
      updateStatus();
    }
  }

  // Apply nesting indentation
  if( flatSubAgent )
  {
    addStyleClass( "LlmSubAgentFlat" );
  }else if( m_nestingLevel > 0 )
  {
    addStyleClass( "LlmSubAgentNesting LlmNestLevel" + to_string(m_nestingLevel) );

    // On phones, a sub-agent is already inside its parents (collapsed) tool-summary row; start it
    //  collapsed too, as its "question" is the full task prompt the sub-agent was given.
    if( m_phoneLayout )
      setCollapsed( true );
  }
}//LlmInteractionDisplay constructor


LlmInteractionDisplay::~LlmInteractionDisplay()
{
  // Nothing to do for the JS counters: they stop themselves once their element is gone.
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

  // A tool round ends with its results; anything else (e.g., an Error after the tool calls) means
  //  the open round won't get any.
  if( turn->type() != LlmInteractionTurn::Type::ToolResult )
    m_openToolRound = nullptr;

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
      finalizeActiveStatus();  // The LLM answered

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

      finalizeActiveStatus();  // The LLM answered

      // A round whose tools have all run already has its duration set, so this only counts while
      //  async tools or sub-agents are still running.
      const bool countRunningTools = m_isLive && !turn->callDuration().has_value();

      if( m_phoneLayout )
      {
        // No request panel on phones; just note what is running in the summary row title.
        ensurePhoneToolSummary();
        m_phoneOpenToolRequests += 1;
        if( request )
        {
          for( const LlmToolCall &call : request->toolCalls() )
          {
            if( std::find( begin(m_phoneToolNames), end(m_phoneToolNames), call.toolName ) == end(m_phoneToolNames) )
              m_phoneToolNames.push_back( call.toolName );
          }
        }//if( request )
        updatePhoneToolSummaryTitle();

        if( countRunningTools )
          startActiveStatus( m_phoneStatus, true, turn );
        break;
      }//if( m_phoneLayout )

      // The request and the results that follow it are shown together, in one panel.
      assert( request );
      if( !request )
        break;

      LlmToolRoundDisplay * const round = m_turnContainer->addNew<LlmToolRoundDisplay>( request );
      m_toolRounds.push_back( round );
      m_openToolRound = round;

      turnDisplay = round->addRequestDisplay();

      if( countRunningTools )
        startActiveStatus( round->statusText(), true, turn );
      else if( turn->callDuration().has_value() )
        showFinalDuration( round->statusText(), turn );
      break;
    }

    case LlmInteractionTurn::Type::ToolResult:
    {
      shared_ptr<LlmToolResults> results =
        dynamic_pointer_cast<LlmToolResults>( turn );

      if( m_phoneLayout )
      {
        ensurePhoneToolSummary();
        m_phoneOpenToolRequests = std::max( 0, m_phoneOpenToolRequests - 1 );
        if( results )
        {
          for( const LlmToolCall &call : results->toolCalls() )
          {
            if( call.status == LlmToolCall::CallStatus::Error )
              m_phoneToolErrors = true;
            if( std::find( begin(m_phoneToolNames), end(m_phoneToolNames), call.toolName ) == end(m_phoneToolNames) )
              m_phoneToolNames.push_back( call.toolName );
          }
        }//if( results )

        turnDisplay = m_phoneToolList->addNew<LlmToolResultsDisplay>( results, m_llmInterface,
                                                                     m_nestingLevel, m_phoneVisuals );
        updatePhoneToolSummaryTitle();

        m_lastRequestTurn = turn;
        m_lastRequestDisplay = nullptr;  // phones use m_phoneStatus
        break;
      }//if( m_phoneLayout )

      // Into the round of the request these are the results of (there should always be one).
      LlmToolRoundDisplay *round = m_openToolRound;
      m_openToolRound = nullptr;
      if( !round )
      {
        round = m_turnContainer->addNew<LlmToolRoundDisplay>( nullptr );
        m_toolRounds.push_back( round );
      }

      turnDisplay = round->addResultsDisplay( results, m_llmInterface, m_nestingLevel );

      // When the results are sent to the LLM, this panel (i.e., its round) shows the "Pending LLM
      //  response" counter.
      m_lastRequestTurn = turn;
      m_lastRequestDisplay = turnDisplay;
      if( turn->callDuration().has_value() )
        showFinalDuration( turnDisplay->statusText(), turn );
      break;
    }

    case LlmInteractionTurn::Type::Error:
    {
      finalizeActiveStatus();  // The request failed or was cancelled

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

      if( m_phoneLayout )
      {
        // The auto-reply is an internal nudge to the LLM, so on phones just note it in the summary.
        ensurePhoneToolSummary();
        WText *note = m_phoneToolList->addNew<WText>( WString::tr("lcw-phone-auto-reply") );
        note->addStyleClass( "LlmPhoneAutoReply" );

        m_lastRequestTurn = turn;
        m_lastRequestDisplay = nullptr;
        break;
      }//if( m_phoneLayout )

      // Reuse the FinalResponse display class since AutoReply also has content
      // We create a temporary FinalResponse wrapper to use the existing display code
      shared_ptr<LlmInteractionFinalResponse> tempResponse =
        make_shared<LlmInteractionFinalResponse>( autoReply->content(), turn->conversation().lock() );
      tempResponse->setThinkingContent( autoReply->thinkingContent() );
      tempResponse->setRawContent( autoReply->rawContent() );

      turnDisplay = m_turnContainer->addNew<LlmInteractionFinalResponseDisplay>( tempResponse );

      // AutoReply doesn't mean conversation is finished - we're waiting for LLM response.  Keyed on
      //  `turn` (not the wrapper), as that is what the request is sent with.
      m_lastRequestTurn = turn;
      m_lastRequestDisplay = turnDisplay;
      if( turn->callDuration().has_value() )
        showFinalDuration( turnDisplay->statusText(), turn );
      break;
    }

    case LlmInteractionTurn::Type::ConversationSummary:
      // Stands in for the older conversations it summarizes, for the LLM; nothing to show.
      break;
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


void LlmInteractionDisplay::ensurePhoneToolSummary()
{
  assert( m_phoneLayout );
  if( m_phoneToolSummary )
    return;

  m_phoneToolSummary = m_turnContainer->addNew<WPanel>();
  m_phoneToolSummary->addStyleClass( "LlmPhoneToolSummary" );
  m_phoneToolSummary->setTitle( WString() );
  m_phoneToolSummary->setCollapsible( true );
  m_phoneToolSummary->setCollapsed( true );
  m_phoneToolList = m_phoneToolSummary->setCentralWidget( std::make_unique<WContainerWidget>() );
  m_phoneToolList->addStyleClass( "LlmPhoneToolList" );

  m_phoneStatus = m_phoneToolSummary->titleBarWidget()->addNew<LlmElapsedText>();
  m_phoneStatus->addStyleClass( "LlmTurnStatus" );

  // Visuals go right after the summary row (and so before the final response).  Sub-agents were
  //  handed their top-level interactions container instead, so their visuals are not hidden inside
  //  this (collapsed) row.
  if( !m_phoneVisuals )
  {
    m_phoneVisuals = m_turnContainer->addNew<WContainerWidget>();
    m_phoneVisuals->addStyleClass( "LlmPhoneVisuals" );
  }
}//ensurePhoneToolSummary()


void LlmInteractionDisplay::updatePhoneToolSummaryTitle()
{
  if( !m_phoneToolSummary )
    return;

  string names;
  for( const string &name : m_phoneToolNames )
    names += (names.empty() ? "" : ", ") + name;

  const char * const key = (m_phoneOpenToolRequests > 0) ? "lcw-phone-tools-running"
                          : (m_phoneToolErrors ? "lcw-phone-tools-used-errors" : "lcw-phone-tools-used");
  const WString title = WString::tr( key ).arg( WString::fromUTF8( names ) );

  // The panel title is XHTML, so escape the model-derived tool names.
  m_phoneToolSummary->setTitle( WString::fromUTF8( Wt::Utils::htmlEncode( title.toUTF8() ) ) );
}//updatePhoneToolSummaryTitle()


void LlmInteractionDisplay::handleSubAgentFinished( shared_ptr<LlmInteractionTurn> turn,
                                                    shared_ptr<LlmInteraction> subAgent )
{
  // Sub-agent finished - the display should already be updated via its own signals
  updateStatus();
}//handleSubAgentFinished(...)


void LlmInteractionDisplay::handleConversationFinished()
{
  // Also emitted for a conversation ended without a final or error turn (e.g., a cancelled
  //  sub-agent), so nothing it shows may still be counting.
  finalizeActiveStatus();

  m_isFinished = true;
  stopStatusTimer();
  updateSummary();
  updateStatus();
}//handleConversationFinished()


void LlmInteractionDisplay::handleStoppedChanged( const bool stopped )
{
  if( stopped == m_isStopped )
    return;
  m_isStopped = stopped;

  const chrono::milliseconds sinceStart = std::max( chrono::milliseconds(0),
          chrono::duration_cast<chrono::milliseconds>( chrono::system_clock::now() - m_interaction->timestamp ) );

  // The round that is stopped is the newest one (on phones, the tool-summary row).
  LlmElapsedText * const roundStatus = m_phoneLayout ? m_phoneStatus
                                     : (m_toolRounds.empty() ? nullptr : m_toolRounds.back()->statusText());

  if( stopped )
  {
    // Nothing runs for it now: it waits on the user to Retry, or Cancel, the paused agent.
    finalizeActiveStatus();
    if( roundStatus )
    {
      roundStatus->setFinalText( WString::tr( "lcw-status-stopped" ) );
      roundStatus->setToolTip( WString::tr( "lcw-status-stopped-tt" ) );
    }
    m_timerText->setFinalText( WString::tr( "lcw-total-time" ).arg( formatDuration( sinceStart ) ) );
  }else
  {
    if( roundStatus )
      roundStatus->setToolTip( WString() );
    m_timerText->start( WString(), sinceStart );
    if( const shared_ptr<LlmInterface> llm = m_llmInterface.lock() )
      resumeActiveStatus( *llm );
  }//if( stopped ) / else

  updateStatus();
}//handleStoppedChanged(...)


void LlmInteractionDisplay::showInterruptedPanel()
{
  // A sub-agent is resumed along with the conversation that invoked it, so has no buttons of its own.
  const shared_ptr<LlmInterface> llm = m_llmInterface.lock();
  if( (m_nestingLevel > 0) || !llm || !llm->canResumeInterrupted( m_interaction ) )
    return;

  // One hidden when the conversation was resumed is replaced, to give this interruption's time and cause.
  if( m_interruptedPanel )
  {
    if( !m_interruptedPanel->isHidden() )
      return;
    m_interruptedPanel->removeFromParent();
    m_interruptedPanel = nullptr;
  }

  const LlmInteraction::Interruption &interruption = m_interaction->interruption.value();

  m_interruptedPanel = m_turnContainer->addNew<WContainerWidget>();
  m_interruptedPanel->addStyleClass( "LlmInterruptedPanel" );

  const char *message_key = "lcw-interrupted-session";
  switch( interruption.cause )
  {
    case LlmInteraction::Interruption::Cause::SessionEnded:    message_key = "lcw-interrupted-session";  break;
    case LlmInteraction::Interruption::Cause::AssistantClosed: message_key = "lcw-interrupted-closed";   break;
    case LlmInteraction::Interruption::Cause::SpectrumChanged: message_key = "lcw-interrupted-spectrum"; break;
  }//switch( interruption.cause )

  WText *message = m_interruptedPanel->addNew<WText>( WString::tr( message_key ).arg( formatTime( interruption.time ) ) );
  message->addStyleClass( "LlmInterruptedMessage" );

  WContainerWidget *buttons = m_interruptedPanel->addNew<WContainerWidget>();
  buttons->addStyleClass( "LlmInterruptedButtons" );

  // The buttons are hidden with the panel, by handleInterruptionChanged(), from within their own click.
  WPushButton *resumeBtn = buttons->addNew<WPushButton>( WString::tr( "lcw-interrupted-resume" ) );
  resumeBtn->addStyleClass( "LlmResumeButton" );
  resumeBtn->clicked().connect( this, [this](){
    if( const shared_ptr<LlmInterface> llm = m_llmInterface.lock() )
      llm->resumeInterrupted( m_interaction );
  } );

  WPushButton *dismissBtn = buttons->addNew<WPushButton>( WString::tr( "lcw-interrupted-dismiss" ) );
  dismissBtn->addStyleClass( "LlmContinueButton" );
  dismissBtn->clicked().connect( this, [this](){
    if( const shared_ptr<LlmInterface> llm = m_llmInterface.lock() )
      llm->abandonInterrupted( m_interaction );
  } );

  WText *help = m_interruptedPanel->addNew<WText>( WString::tr( "lcw-interrupted-help" ) );
  help->addStyleClass( "LlmErrorHelpText" );
}//showInterruptedPanel()


void LlmInteractionDisplay::handleInterruptionChanged()
{
  if( m_interaction->interruption.has_value() )
  {
    // Cut off while shown (e.g., the detectors summed changed): nothing runs for it now.
    m_isStopped = false;
    finalizeActiveStatus();
    stopStatusTimer();
    showInterruptedPanel();
    updateStatus();
    return;
  }//if( interrupted )

  if( m_interruptedPanel )
    m_interruptedPanel->hide();

  // Resumed: its total time counts again, from where it was started.
  const shared_ptr<LlmInterface> llm = m_llmInterface.lock();
  if( !m_isFinished && llm && llm->isConversationActive( m_interaction ) )
  {
    const chrono::system_clock::duration sinceStart = chrono::system_clock::now() - m_interaction->timestamp;
    m_timerText->start( WString(), std::max( chrono::milliseconds(0),
                                  chrono::duration_cast<chrono::milliseconds>( sinceStart ) ) );
  }

  updateStatus();
}//handleInterruptionChanged()


void LlmInteractionDisplay::handleRequestSent( shared_ptr<LlmInteractionTurn> turn )
{
  // A new request ends whatever was being waited on before it (e.g., a tool round).
  finalizeActiveStatus();

  // A request for a finished conversation - a Retry of its error, or a re-prompt after a response
  //  that did not end it (an auto-reply) - means it is running again.
  if( m_isFinished )
  {
    m_isFinished = false;
    const chrono::system_clock::duration sinceStart = chrono::system_clock::now() - m_interaction->timestamp;
    m_timerText->start( WString(), std::max( chrono::milliseconds(0),
                                  chrono::duration_cast<chrono::milliseconds>( sinceStart ) ) );
    updateStatus();
    updateSummary();
  }//if( m_isFinished )

  startActiveStatus( requestStatusText( turn ), false, turn );
}//handleRequestSent(...)


LlmElapsedText *LlmInteractionDisplay::requestStatusText( const shared_ptr<LlmInteractionTurn> &turn )
{
  // The InitialRequest has no panel; the conversation title already shows that wait.
  if( !turn || (turn != m_lastRequestTurn) )
    return nullptr;

  if( m_phoneLayout )
    return m_phoneStatus;

  return m_lastRequestDisplay ? m_lastRequestDisplay->statusText() : nullptr;
}//requestStatusText(...)


void LlmInteractionDisplay::resumeActiveStatus( const LlmInterface &llm )
{
  const std::optional<LlmInterface::InFlightWait> wait = llm.inFlightWait( m_interaction );
  if( !wait || !wait->turn )
    return;

  LlmElapsedText *status = nullptr;
  if( !wait->runningTools )
  {
    status = requestStatusText( wait->turn );
  }else if( m_phoneLayout )
  {
    status = m_phoneStatus;
  }else
  {
    // The running tools are of a tool-call turn, which is shown in a tool round.
    LlmToolRoundDisplay * const round = toolRound( wait->turn );
    status = round ? round->statusText() : nullptr;
  }//if( request ) / else

  const chrono::milliseconds elapsed = chrono::duration_cast<chrono::milliseconds>(
                                                        chrono::steady_clock::now() - wait->since );
  startActiveStatus( status, wait->runningTools, wait->turn, std::max( chrono::milliseconds(0), elapsed ) );
}//resumeActiveStatus(...)


void LlmInteractionDisplay::handleActivityIdle()
{
  // Nothing is in flight any more, so the turn being waited on cant still be running.  (Not the
  //  total-time counter, which is only stopped when the conversation finishes - e.g., an error leaves
  //  it paused, not finished.)
  finalizeActiveStatus();
}//handleActivityIdle()


void LlmInteractionDisplay::updateIdleListener()
{
  const bool counting = (m_activeStatus != nullptr);
  if( counting == m_idleConnection.isConnected() )
    return;

  if( !counting )
  {
    m_idleConnection.disconnect();
    return;
  }

  const shared_ptr<LlmInterface> llm = m_llmInterface.lock();
  if( !llm )
    return;

  m_idleConnection = llm->activityChanged().connect( this, [this]( const LlmInterface::ActivityState state ){
    if( state == LlmInterface::ActivityState::Idle )
      handleActivityIdle();
  } );
}//updateIdleListener()


void LlmInteractionDisplay::startActiveStatus( LlmElapsedText *status, const bool runningTools,
                                               shared_ptr<LlmInteractionTurn> turn,
                                               const chrono::milliseconds alreadyElapsed )
{
  finalizeActiveStatus();

  // The title says what is being waited on, even if there is no counter for it.
  m_waitingOnTools = runningTools;
  updateStatus();

  if( !status )
    return;

  m_activeStatus = status;
  m_activeStatusTurn = turn;
  status->start( WString::tr( runningTools ? "lcw-status-running-tools" : "lcw-status-pending-llm" ),
                 alreadyElapsed );
  updateIdleListener();
}//startActiveStatus(...)


void LlmInteractionDisplay::finalizeActiveStatus()
{
  if( m_waitingOnTools )
  {
    m_waitingOnTools = false;
    updateStatus();
  }

  if( !m_activeStatus )
    return;

  LlmElapsedText * const status = m_activeStatus;
  const shared_ptr<LlmInteractionTurn> turn = m_activeStatusTurn;
  m_activeStatus = nullptr;
  m_activeStatusTurn.reset();

  // The phone summary row is shared by all the calls, so it doesnt show any one calls duration.
  if( status == m_phoneStatus )
    status->setFinalText( WString() );
  else
    showFinalDuration( status, turn );

  updateIdleListener();
}//finalizeActiveStatus()


LlmToolRoundDisplay *LlmInteractionDisplay::toolRound( const shared_ptr<LlmInteractionTurn> &turn ) const
{
  for( LlmToolRoundDisplay * const round : m_toolRounds )
  {
    if( round->hasTurn( turn ) )
      return round;
  }
  return nullptr;
}//LlmToolRoundDisplay *toolRound(...)


void LlmInteractionDisplay::showFinalDuration( LlmElapsedText *status,
                                               const shared_ptr<LlmInteractionTurn> &turn )
{
  if( !status )
    return;

  // A tool rounds status shows both the tools and the LLMs time.
  LlmToolRoundDisplay * const round = toolRound( turn );
  if( round && (status == round->statusText()) )
  {
    round->showFinalDurations();
    return;
  }

  if( !turn || !turn->callDuration().has_value() )
  {
    status->setFinalText( WString() );
    status->setToolTip( WString() );
    return;
  }

  // A ToolCall turn's duration is the time to run its tools here; any other turn's is the time the
  //  LLM took to respond to it.  They look alike, so say which.
  const bool isTools = (turn->type() == LlmInteractionTurn::Type::ToolCall);
  status->setFinalText( WString::tr( isTools ? "lcw-duration-tools" : "lcw-duration-llm" )
                          .arg( formatDuration( turn->callDuration().value() ) ) );
  status->setToolTip( WString::tr( isTools ? "lcw-duration-tools-tt" : "lcw-duration-llm-tt" ) );
}//showFinalDuration(...)


void LlmInteractionDisplay::updateStatus()
{
  WString titleWithStatus = getTitleWithStatus();

  // On phones the title bar is a single line, so dont let newlines in the question (which a plain-
  //  text WText renders as <br/>) break it up - e.g., a sub-agents "Context: ...\nTask: ..." prompt.
  if( m_phoneLayout )
  {
    string title = titleWithStatus.toUTF8();
    std::replace( begin(title), end(title), '\n', ' ' );
    std::replace( begin(title), end(title), '\r', ' ' );
    titleWithStatus = WString::fromUTF8( title );
  }

  //cerr << "LlmInteractionDisplay::updateStatus(): Setting title to: '"
  //     << titleWithStatus.toUTF8() << "', m_isFinished=" << m_isFinished << endl;
  m_statusText->setText( titleWithStatus );

  if( m_isFinished )
    stopStatusTimer();
}//updateStatus()


void LlmInteractionDisplay::updateSummary()
{
  assert( m_summaryDiv );
  if( !m_summaryDiv )
    return;

  m_summaryDiv->clear();
  m_summaryDiv->setHidden( false );

  // Start time
  const string startTimeStr = formatTime( m_interaction->timestamp );
  WText *startTime = m_summaryDiv->addNew<WText>( WString::tr( "lcw-summary-start" ).arg( startTimeStr ) );
  startTime->addStyleClass( "LlmSummaryItem" );

  // Total duration (if conversation is finished)
  if( m_isFinished )
  {
    const WString duration = WString::tr( "lcw-summary-duration" ).arg( formatDuration( totalDuration( *m_interaction ) ) );
    WText *durationText = m_summaryDiv->addNew<WText>( " | " + duration );
    durationText->addStyleClass( "LlmSummaryItem" );
  }

  // Token usage.  Note the input count is cumulative across all of this conversation's API calls
  // (the growing history is re-sent each tool-call round), which is why it can be large; the cache
  // breakdown below shows how much of that input was re-processed cheaply from the provider's cache.
  if( m_interaction->totalTokens.has_value() && m_interaction->totalTokens.value() > 0 )
  {
    const size_t promptTok = m_interaction->promptTokens.value_or(0);
    const size_t completionTok = m_interaction->completionTokens.value_or(0);

    const WString tokenStr = WString::tr( "lcw-summary-tokens" ).arg( to_string(promptTok) )
                         .arg( to_string(completionTok) ).arg( to_string(m_interaction->totalTokens.value()) );
    WText *tokens = m_summaryDiv->addNew<WText>( " | " + tokenStr );
    tokens->addStyleClass( "LlmSummaryItem" );

    // Cache usage: how much of the (cumulative) input was served from / written to cache.
    if( m_interaction->cachedTokens.has_value() )
    {
      const size_t cachedTok = m_interaction->cachedTokens.value();
      WString cacheStr = WString::tr( "lcw-summary-cache" ).arg( to_string(cachedTok) );
      if( promptTok > 0 )
      {
        const int pct = static_cast<int>( (100.0 * cachedTok / promptTok) + 0.5 );
        cacheStr = WString::tr( "lcw-summary-cache-pct" ).arg( to_string(cachedTok) ).arg( pct );
      }
      const size_t writtenTok = m_interaction->cacheCreationTokens.value_or(0);
      if( writtenTok > 0 )
        cacheStr = WString::tr( "lcw-summary-cache-written" ).arg( cacheStr ).arg( to_string(writtenTok) );

      WText *cache = m_summaryDiv->addNew<WText>( " | " + cacheStr );
      cache->addStyleClass( "LlmSummaryItem" );
    }
  }

  // Cost - only shown when the provider reports it (e.g. OpenRouter).
  if( m_interaction->cost.has_value() )
  {
    WText *cost = m_summaryDiv->addNew<WText>( " | " + formatCost( m_interaction->cost.value() ) );
    cost->addStyleClass( "LlmSummaryItem" );
  }
}//updateSummary()


WString LlmInteractionDisplay::getTitleWithStatus() const
{
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

  const WString question = WString::fromUTF8( truncateString( questionContent, 60 ) );

  // Status: stopped while waiting on a paused agent, or ended on an error still waiting on a Retry or Skip
  const LlmInteractionError * const lastError = m_interaction->responses.empty() ? nullptr
                     : dynamic_cast<const LlmInteractionError *>( m_interaction->responses.back().get() );
  const bool awaitingErrorChoice = m_isFinished && lastError && !lastError->skippedByUser();
  if( !m_isFinished && m_interaction->interruption.has_value() )
    return WString::tr( "lcw-title-interrupted" ).arg( question );

  // Also one ended before it got an answer (e.g., an interrupted conversation that was not resumed).
  if( m_isStopped || awaitingErrorChoice || (m_isFinished && m_interaction->isAwaitingWork()) )
    return WString::tr( "lcw-title-stopped" ).arg( question );

  if( m_isFinished )
    return WString::tr( "lcw-title-complete" ).arg( question );

  // Tools (sub-agents, async tools) still running; a round whose tools have all run has its results
  //  sent right away, so the LLM is what is being waited on.
  if( m_waitingOnTools )
    return WString::tr( "lcw-title-running-tools" ).arg( question );

  return WString::tr( "lcw-title-pending-llm" ).arg( question );
}//getTitleWithStatus()


void LlmInteractionDisplay::stopStatusTimer()
{
  if( m_timerText )
    m_timerText->setFinalText( WString::tr( "lcw-total-time" ).arg( formatDuration( totalDuration( *m_interaction ) ) ) );
}//stopStatusTimer()


void LlmInteractionDisplay::createMenuIcon()
{
  // Create menu button in title bar; the menu is built lazily on first click
  // (showMenu) so interactions whose menu is never opened cost nothing.
  // All of its items are for developing the assistant, so it is only shown with debug options on.
  m_menuIcon = titleBarWidget()->addNew<WPushButton>();
  m_menuIcon->setStyleClass( "RoundMenuIcon InvertInDark dropdown-toggle Wt-btn LlmDebugOnly" );
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
  WidgetUtils::applyButtonRole( closeBtn, WidgetUtils::ButtonRole::Dismiss );
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
