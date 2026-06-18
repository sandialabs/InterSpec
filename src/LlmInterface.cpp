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

#include <regex>
#include <stack>
#include <atomic>
#include <chrono>
#include <memory>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#if( PERFORM_DEVELOPER_CHECKS )
#include <mutex>
#endif

#include <Wt/WApplication>
#include <Wt/WResource>
#include <Wt/WWebWidget>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WJavaScriptSlot>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include "InterSpec/InterSpec.h"

#include "SpecUtils/SpecFile.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/LlmConversationHistory.h"


#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
#include <fstream>

#include "SpecUtils/DateTime.h"
#endif
#include "SpecUtils/StringAlgo.h"

static_assert( USE_LLM_INTERFACE, "This file shouldnt be compiled unless USE_LLM_INTERFACE is true" );

using namespace std;
using json = nlohmann::json;

// Max consecutive auto-replies when sub-agent is not in a final state (resets on tool use)
static constexpr size_t sm_max_consecutive_nonfinal_autoreplies = 4;

// Max total non-final-state auto-replies per sub-agent invocation.
//  The open-source models can need significant nudging during nulcide ID of complex spectra
static constexpr size_t sm_max_total_nonfinal_autoreplies = 100;

// Forward declarations
static std::string sanitizeUtf8( const std::string &str );
static std::string sanitizeToolName( const std::string &toolName );


/** Generate a tool call ID compatible with OpenAI's format: "call_" + 24 alphanumeric chars.
 
 This is used when parsing text-based tool calls from models that don't natively support
 the OpenAI tool call format, to create IDs that match the expected format.
 */
static std::string generateToolCallId()
{
  static const char alphanum[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  
  // Use combination of timestamp and counter for uniqueness
  const auto epoch_ticks = std::chrono::duration_cast<std::chrono::microseconds>(
    std::chrono::system_clock::now().time_since_epoch() ).count();
  
  std::string id = "call_";
  
  // Use epoch as seed, combined with a static counter for same-microsecond calls
  static std::atomic<uint32_t> counter{0};
  uint64_t seed = static_cast<uint64_t>( epoch_ticks ) ^ (static_cast<uint64_t>( counter.fetch_add(1) ) << 48);
  
  // Generate 24 alphanumeric characters using a simple LCG
  for( int i = 0; i < 24; ++i )
  {
    id += alphanum[seed % 62];
    seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
  }
  
  return id;
}//std::string generateToolCallId()


#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
/** Write cache debug information to files for inspection when cache hit rate is low.
 * Creates timestamped files with current and previous request JSONs.
 */

/*
static void writeCacheDebugFiles( const std::string &currentRequest,
                                   const std::string &previousRequest,
                                   const double matchPercent )
{
  try
  {
    const std::string timestamp = SpecUtils::to_common_string( std::chrono::system_clock::now(), true );
    // Convert timestamp to safe filename format (replace colons and spaces)
    std::string safeTimestamp = timestamp;
    std::replace( safeTimestamp.begin(), safeTimestamp.end(), ':', '-' );
    std::replace( safeTimestamp.begin(), safeTimestamp.end(), ' ', '_' );
    
    const std::string basePath = "/tmp/llm_cache_debug_" + safeTimestamp;
    const std::string currentPath = basePath + "_current.json";
    const std::string previousPath = basePath + "_previous.json";
    const std::string summaryPath = basePath + "_summary.txt";
    
    // Write current request
    std::ofstream currentFile( currentPath );
    if( currentFile.is_open() )
    {
      currentFile << currentRequest;
      currentFile.close();
    }
    
    // Write previous request
    std::ofstream previousFile( previousPath );
    if( previousFile.is_open() )
    {
      previousFile << previousRequest;
      previousFile.close();
    }
    
    // Write summary
    std::ofstream summaryFile( summaryPath );
    if( summaryFile.is_open() )
    {
      summaryFile << "Cache hit prediction: " << std::fixed << std::setprecision(1) 
                  << matchPercent << "%" << std::endl;
      summaryFile << "Current request length: " << currentRequest.length() << " chars" << std::endl;
      summaryFile << "Previous request length: " << previousRequest.length() << " chars" << std::endl;
      summaryFile << std::endl;
      summaryFile << "Files written:" << std::endl;
      summaryFile << "  Current:  " << currentPath << std::endl;
      summaryFile << "  Previous: " << previousPath << std::endl;
      summaryFile << std::endl;
      summaryFile << "To diff these files, run:" << std::endl;
      summaryFile << "  diff -u \"" << previousPath << "\" \"" << currentPath << "\"" << std::endl;
      summaryFile << "Or for a side-by-side view:" << std::endl;
      summaryFile << "  sdiff -w 200 \"" << previousPath << "\" \"" << currentPath << "\" | less" << std::endl;
      summaryFile.close();
      
      cout << "Cache debug files written to: " << basePath << "_*.{json,txt}" << endl;
      cout << "Summary: " << summaryPath << endl;
    }
  }
  catch( const std::exception &e )
  {
    cerr << "Error writing cache debug files: " << e.what() << endl;
  }
}
 */
#endif // PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER

LlmInterface::LlmInterface( InterSpec* interspec, const std::shared_ptr<const LlmConfig> &config )
  : Wt::Signals::trackable(),
    m_interspec(interspec),
    m_block_tool_calls(false),
    m_config( config ),
    m_tool_registry( config ? make_shared<LlmTools::ToolRegistry>(*config) : shared_ptr<LlmTools::ToolRegistry>() ),
    m_history(std::make_shared<LlmConversationHistory>()),
    m_debug_stream( nullptr ),
    m_debug_file( nullptr ),
    m_instanceId( std::to_string( []{ static std::atomic<int> s_instance_count{0}; return s_instance_count.fetch_add(1); }() ) ),
    m_responseSignal(std::make_unique<Wt::JSignal<std::string, int>>(interspec, "llmResponse_" + m_instanceId)),
    m_nextRequestId(1)
{
  if( !m_interspec )
    throw std::runtime_error("InterSpec instance cannot be null");
  
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );
  
  if( !m_tool_registry )
    throw std::logic_error( "LlmInterface: couldnt create tool registry from config" );
  
  // Initialize debug logging if configured
  if( !m_config->llmApi.debug_file.empty() )
  {
    const string debug_dest = m_config->llmApi.debug_file;
    if( debug_dest == "stdout" )
    {
      m_debug_stream = &std::cout;
    }else if( debug_dest == "stderr" )
    {
      m_debug_stream = &std::cerr;
    }else
    {
      // Open file for writing
      m_debug_file = std::make_unique<std::ofstream>( debug_dest, std::ios::out | std::ios::app );
      if( m_debug_file->is_open() )
      {
        m_debug_stream = m_debug_file.get();
      }else
      {
        std::cerr << "Warning: Failed to open debug log file: " << debug_dest << std::endl;
      }
    }
  }
  
  // Get session ID through WApplication
  string sessionId = "unknown";
  if (Wt::WApplication* app = Wt::WApplication::instance()) {
    sessionId = app->sessionId();
  }
  //cout << "LlmInterface created for session: " << sessionId << endl;
  
  // Connect the JavaScript response signal
  m_responseSignal->connect(this, &LlmInterface::handleJavaScriptResponse);
  
  setupJavaScriptBridge();
}

LlmInterface::~LlmInterface()
{
  // Unregister this instance's callback from the JavaScript registry
  if( Wt::WApplication *app = Wt::WApplication::instance() )
    app->doJavaScript( "if(window.llmResponseCallbacks) delete window.llmResponseCallbacks['" + m_instanceId + "'];" );

  cout << "LlmInterface destroyed (instanceId=" << m_instanceId << ")" << endl;
}


void LlmInterface::setBlockToolCalls( bool block )
{
  m_block_tool_calls = block;
}


shared_ptr<LlmInteraction> LlmInterface::triggerManualCompaction()
{
  if( !m_history )
    return nullptr;

  // Need at least a few conversations to compact
  const vector<shared_ptr<LlmInteraction>> &convos = m_history->getConversations();
  const size_t numNonSummary = m_history->hasSummary() ? (convos.size() - 1) : convos.size();
  if( numNonSummary < 2 )
    return nullptr;

  const string summarizationPrompt = m_history->buildSummarizationPrompt();
  if( summarizationPrompt.empty() )
    return nullptr;

  if( m_debug_stream )
    (*m_debug_stream) << "=== Manual context compaction triggered ===" << endl;

  return sendCompactionRequest( summarizationPrompt,
    "[Context Compaction] Summarizing conversation history..." );
}


shared_ptr<LlmInteraction> LlmInterface::sendCompactionRequest( const string &summarizationPrompt,
                                                                  const string &displayMessage )
{
  assert( m_history && m_config );

  // Create a visible system conversation so the user sees the compaction in the widget
  shared_ptr<LlmInteraction> summarizationConvo
    = m_history->addSystemMessageToMainConversation( displayMessage );

  // Build a minimal request (no tools, just the summarization prompt)
  json sumRequest;
  sumRequest["model"] = m_config->llmApi.model();

  if( m_config->llmApi.maxTokens() > 0 )
  {
    const string &modelName = m_config->llmApi.model();
    if( (modelName.find("gpt-") != string::npos)
       || (modelName.find("o1") != string::npos)
       || (modelName.find("o3") != string::npos)
       || (modelName.find("o4") != string::npos)
       || (modelName.find("gemini") != string::npos) )
    {
      sumRequest["max_completion_tokens"] = m_config->llmApi.maxTokens();
    }else
    {
      sumRequest["max_tokens"] = m_config->llmApi.maxTokens();
    }
  }//if( m_config->llmApi.maxTokens() > 0 )

  if( m_config->llmApi.temperature().has_value() )
    sumRequest["temperature"] = m_config->llmApi.temperature().value();

  static const char * const sm_defaultCompactionPrompt =
    "You are a specialized assistant that summarizes conversation history for a nuclear"
    " spectroscopy application (InterSpec). Produce a concise summary that preserves:"
    "\n- Key user questions and the conclusions reached"
    "\n- Important analysis results (nuclide identifications, activity values, shielding estimates)"
    "\n- Tool calls that modified application state (peaks added/removed, sources assigned,"
    " calibration changes) and their outcomes"
    "\n- Any unresolved issues or follow-up items"
    "\nOmit tool call details, intermediate reasoning, and error retries."
    " Focus on what was asked, what was found, and what changed.";

  const string &compactionPrompt = m_config->llmApi.compactionSystemPrompt.empty()
    ? string( sm_defaultCompactionPrompt )
    : m_config->llmApi.compactionSystemPrompt;

  json sumMessages = json::array();

  json systemMsg;
  systemMsg["role"] = "system";
  systemMsg["content"] = compactionPrompt;
  sumMessages.push_back( systemMsg );

  json userMsg;
  userMsg["role"] = "user";
  userMsg["content"] = summarizationPrompt;
  sumMessages.push_back( userMsg );

  sumRequest["messages"] = sumMessages;

  // Send the summarization request and mark it
  std::pair<int,string> sumReqResult = makeTrackedApiCall( sumRequest, summarizationConvo );

  {
    std::map<int, PendingRequest>::iterator it = m_pendingRequests.find( sumReqResult.first );
    if( it != m_pendingRequests.end() )
      it->second.isSummarizationRequest = true;
  }

  if( m_debug_stream )
    (*m_debug_stream) << "Compaction request sent with ID: " << sumReqResult.first << endl;

  return summarizationConvo;
}


void LlmInterface::setDebugFile( const std::string &filePath )
{
  m_debug_stream = nullptr;
  m_debug_file.reset();

  if( filePath.empty() )
    return;

  m_debug_file = std::make_unique<std::ofstream>( filePath, std::ios::out | std::ios::app );
  if( m_debug_file->is_open() )
  {
    m_debug_stream = m_debug_file.get();
  }else
  {
    std::cerr << "LlmInterface::setDebugFile: failed to open '" << filePath << "'" << std::endl;
    m_debug_file.reset();
  }
}


void LlmInterface::resetWithConfig( const std::shared_ptr<const LlmConfig> &config )
{
  if( !config || !config->llmApi.enabled )
    throw std::logic_error( "LlmInterface::resetWithConfig: config is null or LLM API not enabled" );

  // Create the tool registry first - if this throws, we haven't modified any state yet
  std::shared_ptr<const LlmTools::ToolRegistry> tool_registry
    = make_shared<LlmTools::ToolRegistry>( *config );

  m_config = config;
  m_tool_registry = tool_registry;
  m_history = make_shared<LlmConversationHistory>();
  m_pendingRequests.clear();
  m_deferredToolResults.clear();
  m_currentConversation.reset();
  m_summarizationPendingConvo.reset();
  m_serializedCache.clear();
  m_nextRequestId = 1;

  // Re-initialize debug logging
  m_debug_stream = nullptr;
  m_debug_file.reset();

  if( !m_config->llmApi.debug_file.empty() )
  {
    const string debug_dest = m_config->llmApi.debug_file;
    if( debug_dest == "stdout" )
    {
      m_debug_stream = &std::cout;
    }else if( debug_dest == "stderr" )
    {
      m_debug_stream = &std::cerr;
    }else
    {
      // Open file for writing
      m_debug_file = std::make_unique<std::ofstream>( debug_dest, std::ios::out | std::ios::app );
      if( m_debug_file->is_open() )
      {
        m_debug_stream = m_debug_file.get();
      }else
      {
        std::cerr << "Warning: Failed to open debug log file: " << debug_dest << std::endl;
      }
    }
  }
}


std::shared_ptr<const LlmTools::ToolRegistry> LlmInterface::toolRegistry()
{
  return m_tool_registry;
}


shared_ptr<LlmInteraction> LlmInterface::sendUserMessage( const std::string &message,
                                                          vector<LlmToolCall::ImageContent> images )
{
  if( !isConfigured() )
    throw std::logic_error( "LLM interface is not properly configured" );

  if( m_debug_stream )
  {
    (*m_debug_stream) << "User message: " << message << endl;
    if( !images.empty() )
      (*m_debug_stream) << "  (with " << images.size() << " attached image(s))" << endl;
  }

  // Add to history FIRST to ensure it's there before any async responses
  shared_ptr<LlmInteraction> convo = images.empty()
    ? m_history->addUserMessageToMainConversation( message )
    : m_history->addUserMessageToMainConversation( message, std::move( images ) );
  assert( convo );
  if( m_debug_stream )
    (*m_debug_stream) << "Starting new conversation with ID: " << convo->conversationId << endl;

  // Initialize state machine if the agent has one
  //initializeStateMachineForConversation( convo );

  // Debug: Show what's actually in history
  if( m_debug_stream )
  {
    (*m_debug_stream) << "Added user message to history. History now has " << m_history->getConversations().size() << " conversations" << endl;
    (*m_debug_stream) << "Current history contents:" << endl;
    
    for (size_t i = 0; i < m_history->getConversations().size(); ++i) {
      const auto& conv = m_history->getConversations()[i];
      string typeStr;
      switch (conv->type) {
        case LlmInteraction::Type::User: typeStr = "user"; break;
        case LlmInteraction::Type::System: typeStr = "system"; break;
      }
      
      // Get content from InitialRequest
      string contentPreview;
      if( !conv->responses.empty() && conv->responses.front()->type() == LlmInteractionTurn::Type::InitialRequest )
      {
        const LlmInteractionInitialRequest *initialReq = dynamic_cast<const LlmInteractionInitialRequest*>(conv->responses.front().get());
        if( initialReq )
          contentPreview = initialReq->content();
      }
      
      const string preview = contentPreview.length() > 50 ? contentPreview.substr(0, 50) + "..." : contentPreview;
      (*m_debug_stream) << "  " << i << ". " << typeStr << ": " << preview << " (responses: " << conv->responses.size() << ")" << endl;
    }
  }//if( m_debug_stream )
  
  // Check if context summarization is needed before building the request
  const int contextLimit = m_config ? m_config->llmApi.contextLengthLimit() : 0;
  if( (contextLimit > 0) && m_history->shouldSummarize( static_cast<size_t>( contextLimit ) ) )
  {
    if( m_debug_stream )
      (*m_debug_stream) << "=== Context summarization triggered ===" << endl;

    const string summarizationPrompt = m_history->buildSummarizationPrompt();

    if( !summarizationPrompt.empty() )
    {
      // Store the user's real conversation as pending - it will be sent after summarization completes
      m_summarizationPendingConvo = convo;

      sendCompactionRequest( summarizationPrompt,
        "[Context Compaction] Summarizing older conversation history to manage context length..." );

      return convo;
    }
  }//if( context summarization needed )

  // Build API request (buildMessagesArray will include the message from history + current message)
  json requestJson = buildMessagesArray( convo );

  // Make tracked API call
  std::pair<int,std::string> request_id_content = makeTrackedApiCall( requestJson, convo );

  if( m_debug_stream )
    (*m_debug_stream) << "Sent user message with request ID: " << request_id_content.first << endl;

  // Store raw JSON request in the InitialRequest turn's rawContent field
  if( !convo->responses.empty() )
  {
    std::shared_ptr<LlmInteractionTurn> firstTurn = convo->responses.front();
    if( firstTurn && firstTurn->type() == LlmInteractionTurn::Type::InitialRequest )
      firstTurn->setRawContent( std::move(request_id_content.second) );
  }

  return convo;
}

void LlmInterface::sendSystemMessage( const std::string &message )
{
  if( !isConfigured() )
    throw std::runtime_error("LLM interface is not properly configured");
  
  if( m_debug_stream )
    (*m_debug_stream) << "System message: " << message << endl;

  std::shared_ptr<LlmInteraction> convo = m_history->addSystemMessageToMainConversation( message );

  // Initialize state machine if the agent has one
  //initializeStateMachineForConversation( convo );

  // Build API request (marked as system-generated)
  json requestJson = buildMessagesArray( convo );
  
  // Make tracked API call
  std::pair<int,std::string> request_id_content = makeTrackedApiCall( requestJson, convo );
  if( m_debug_stream )
    (*m_debug_stream) << "Sent system message with request ID: " << request_id_content.first << endl;

  // Store raw JSON request in the InitialRequest turn's rawContent field
  if( !convo->responses.empty() )
  {
    std::shared_ptr<LlmInteractionTurn> firstTurn = convo->responses.front();
    if( firstTurn && firstTurn->type() == LlmInteractionTurn::Type::InitialRequest )
      firstTurn->setRawContent( std::move(request_id_content.second) );
  }
}


std::shared_ptr<LlmConversationHistory> LlmInterface::getHistory() const {
  return m_history;
}

void LlmInterface::setHistory(std::shared_ptr<LlmConversationHistory> history) {
  m_history = history ? history : std::make_shared<LlmConversationHistory>();
  m_serializedCache.clear();
}

std::shared_ptr<const LlmConfig> LlmInterface::config() const
{
  return m_config;
}


bool LlmInterface::isConfigured() const {
  return (m_config && m_config->llmApi.enabled && !m_config->llmApi.apiEndpoint().empty());
}


Wt::Signal<>& LlmInterface::conversationFinished() {
  return m_conversationFinished;
}

Wt::Signal<>& LlmInterface::responseError() {
  return m_responseError;
}

bool LlmInterface::isRequestPending(int requestId) const {
  return m_pendingRequests.find(requestId) != m_pendingRequests.end();
}


/** Sanitize a tool name by removing angle brackets and any characters after them.
 
 Some LLM models return tool names with extra characters in angle brackets, like
 "get_identified_sources<|channel|>" when it should be "get_identified_sources".
 This function removes the angle bracket and everything after it, then trims whitespace.
 
 @param toolName The potentially malformed tool name
 @return Sanitized tool name
 */
static std::string sanitizeToolName( const std::string &toolName )
{
  std::string result = toolName;
  
  // Find the first '<' character
  const size_t pos = result.find( '<' );
  if( pos != std::string::npos )
  {
    // Remove '<' and everything after it
    result = result.substr( 0, pos );
  }
  
  // Trim whitespace
  SpecUtils::trim( result );
  
  return result;
}


/** Resolve a tool name case-insensitively against the tool registry.

 Some LLMs return tool names in uppercase (e.g., GET_UNIDENTIFIED_PEAKS) when the actual
 registered tool name is lowercase (get_unidentified_peaks).  This function tries to find
 the correct registered name.

 @param candidateName The tool name to resolve (may be any case)
 @param registry The tool registry to look up against
 @return The matching registered tool name if found, or the original candidateName unchanged
 */
static std::string resolveToolNameCaseInsensitive( const std::string &candidateName,
                                                    const LlmTools::ToolRegistry &registry )
{
  // Fast path: exact match
  if( registry.getTool( candidateName ) )
    return candidateName;

  // Try simple lowercase (covers the common UPPERCASE -> lowercase case)
  const std::string lowered = SpecUtils::to_lower_ascii_copy( candidateName );
  if( registry.getTool( lowered ) )
    return lowered;

  // Full case-insensitive search as fallback
  for( const auto &entry : registry.getTools() )
  {
    if( SpecUtils::iequals_ascii( entry.first, candidateName ) )
      return entry.first;
  }

  return candidateName;
}//resolveToolNameCaseInsensitive


/** Extract JSON from a markdown code fence if one is present.

 LLMs frequently wrap their JSON output in triple-backtick blocks, e.g.:
   ```json
   { ... }
   ```
 This function strips the fence and returns the inner content.  If no
 recognisable fence is found, the original string is returned unchanged.
 **/
static std::string extractFromMarkdown( const std::string &s )
{
  // Look for ```json or ``` followed by content and a closing ```.
  // We try the labelled variant first, then unlabelled.
  static const std::string patterns[] = { "```json", "```JSON", "```" };
  for( const std::string &fence : patterns )
  {
    const size_t start = s.find( fence );
    if( start == std::string::npos )
      continue;

    const size_t contentStart = start + fence.size();
    // Skip the newline immediately after the opening fence
    size_t bodyStart = contentStart;
    while( bodyStart < s.size() && (s[bodyStart] == '\r' || s[bodyStart] == '\n') )
      ++bodyStart;

    const size_t closePos = s.find( "```", contentStart );
    if( closePos == std::string::npos )
      continue;

    // For the unlabelled fence we need the body to look like JSON
    std::string body = s.substr( bodyStart, closePos - bodyStart );
    // Trim trailing whitespace from body
    const size_t lastNonSpace = body.find_last_not_of( " \t\r\n" );
    if( lastNonSpace != std::string::npos )
      body = body.substr( 0, lastNonSpace + 1 );

    if( body.empty() )
      continue;

    // For the unlabelled fence, require the body to start with { or [
    if( fence == "```" && body[0] != '{' && body[0] != '[' )
      continue;

    return body;
  }

  return s;
}//extractFromMarkdown(...)


/** Insert missing commas between object fields and array elements.

 Some LLMs omit the comma separator between fields, e.g.:
   { "a": 1
     "b": 2 }
 This function does a line-by-line pass and inserts a trailing comma on lines
 that end with a complete JSON value when the next non-empty line is neither
 a closing delimiter nor already preceded by one.

 Only operates on multi-line input; single-line JSON is returned unchanged
 because the regex approach needed there risks corrupting string values.
 **/
static std::string fixMissingCommas( const std::string &s )
{
  // Only attempt for multi-line strings
  if( s.find('\n') == std::string::npos )
    return s;

  // Split into lines
  std::vector<std::string> lines;
  {
    std::istringstream stream( s );
    std::string line;
    while( std::getline( stream, line ) )
      lines.push_back( line );
  }

  int depth = 0;
  std::vector<std::string> result;
  result.reserve( lines.size() );

  for( size_t li = 0; li < lines.size(); ++li )
  {
    const std::string &line = lines[li];

    // Process the line to track depth and find where it ends structurally
    std::string processed;
    bool inStr = false;
    bool esc = false;
    for( const char c : line )
    {
      if( esc ) { esc = false; processed += c; continue; }
      if( c == '\\' && inStr ) { esc = true; processed += c; continue; }
      if( c == '"' ) { inStr = !inStr; processed += c; continue; }
      if( !inStr )
      {
        if( c == '{' || c == '[' ) ++depth;
        else if( c == '}' || c == ']' ) --depth;
      }
      processed += c;
    }

    // Determine trimmed content
    const size_t trimStart = processed.find_first_not_of( " \t\r\n" );
    const size_t trimEnd   = processed.find_last_not_of( " \t\r\n" );
    const std::string trimmed = (trimStart == std::string::npos)
                                ? "" : processed.substr( trimStart, trimEnd - trimStart + 1 );

    // Lines that already have a comma, colon, or open-delimiter don't need one
    if( trimmed.empty()
        || trimmed.back() == ','
        || trimmed.back() == ':'
        || trimmed.back() == '{'
        || trimmed.back() == '[' )
    {
      result.push_back( processed );
      continue;
    }

    bool needsComma = false;
    if( depth > 0 && li + 1 < lines.size() )
    {
      // Find next non-empty line
      size_t nextLi = li + 1;
      while( nextLi < lines.size() )
      {
        const size_t ns = lines[nextLi].find_first_not_of( " \t\r\n" );
        if( ns != std::string::npos )
          break;
        ++nextLi;
      }

      if( nextLi < lines.size() )
      {
        const size_t ns = lines[nextLi].find_first_not_of( " \t\r\n" );
        const char nextFirst = (ns != std::string::npos) ? lines[nextLi][ns] : '\0';
        if( nextFirst != '\0' && nextFirst != '}' && nextFirst != ']' )
        {
          // Current line ends with a complete value?
          const char last = trimmed.back();
          const bool endsWithValue = (last == '"' || last == '}' || last == ']'
                                      || (last >= '0' && last <= '9'));
          // Also check for boolean/null endings
          const bool endsWithKeyword = (trimmed.size() >= 4
                                        && (trimmed.substr( trimmed.size() - 4 ) == "true"
                                            || trimmed.substr( trimmed.size() - 4 ) == "null"))
                                       || (trimmed.size() >= 5
                                           && trimmed.substr( trimmed.size() - 5 ) == "false");
          if( endsWithValue || endsWithKeyword )
            needsComma = true;
        }
      }
    }

    if( needsComma )
      result.push_back( processed + "," );
    else
      result.push_back( processed );
  }//for( lines )

  std::string out;
  out.reserve( s.size() + lines.size() );
  for( size_t i = 0; i < result.size(); ++i )
  {
    if( i > 0 )
      out += '\n';
    out += result[i];
  }
  return out;
}//fixMissingCommas(...)


/** Sanitize a JSON string to make it more likely to parse successfully.

 This function attempts to fix common JSON formatting issues:
 - Trims leading and trailing whitespace
 - Removes BOM (Byte Order Mark) if present
 - Replaces invalid UTF-8 sequences
 - Strips null bytes

 @param jsonStr The potentially malformed JSON string
 @return Sanitized JSON string
 */
static std::string sanitizeJsonString( const std::string &jsonStr )
{
  std::string result = jsonStr;

  // Trim leading and trailing whitespace
  SpecUtils::trim( result );

  // Remove UTF-8 BOM if present
  if( result.size() >= 3 &&
      static_cast<unsigned char>(result[0]) == 0xEF &&
      static_cast<unsigned char>(result[1]) == 0xBB &&
      static_cast<unsigned char>(result[2]) == 0xBF )
  {
    result = result.substr(3);
  }

  // Remove null bytes which can cause parse failures
  result.erase( std::remove(result.begin(), result.end(), '\0'), result.end() );

  // Extract JSON from markdown code fences (e.g. ```json ... ```) that LLMs emit.
  // Only do this if the input doesn't already look like valid JSON - otherwise we might
  // incorrectly match backticks that appear inside JSON string values.
  if( !result.empty() && result[0] != '{' && result[0] != '[' )
  {
    result = extractFromMarkdown( result );
    SpecUtils::trim( result );  // re-trim after extraction
  }

  // If the string doesn't start with a valid JSON value character, strip leading garbage
  // until we reach one.  This handles LLM artefacts like "<tool_call>{...}" or other
  // non-JSON prefixes that some models emit before the actual payload.
  // Valid JSON root characters: '{', '[', '"', digit, '-', 't'(rue), 'f'(alse), 'n'(ull).
  // For the search we only anchor on '{', '[', '"' since digits/letters are too likely
  // to appear inside the garbage prefix itself (e.g. "<foo1>{...}").
  if( !result.empty() )
  {
    const char first = result[0];
    const bool isValidJsonStart = (first == '{' || first == '[' || first == '"'
                                   || first == '-' || (first >= '0' && first <= '9')
                                   || first == 't' || first == 'f' || first == 'n');
    if( !isValidJsonStart )
    {
      const size_t jsonStart = result.find_first_of( "{[\"" );
      if( jsonStart != std::string::npos )
        result = result.substr( jsonStart );
      // If no recognisable JSON start is found, leave result as-is and let parsing fail.
    }
  }

  // Replace any invalid UTF-8 sequences
  // This is a simple implementation that handles common cases
  std::string cleaned;
  cleaned.reserve( result.size() );

  for( size_t i = 0; i < result.size(); ++i )
  {
    unsigned char c = static_cast<unsigned char>(result[i]);

    // Valid single-byte UTF-8 (ASCII)
    if( c < 0x80 )
    {
      cleaned += c;
      continue;
    }

    // Check for valid multi-byte UTF-8 sequences
    size_t seqLen = 0;
    if( (c & 0xE0) == 0xC0 ) seqLen = 2;      // 110xxxxx
    else if( (c & 0xF0) == 0xE0 ) seqLen = 3; // 1110xxxx
    else if( (c & 0xF8) == 0xF0 ) seqLen = 4; // 11110xxx

    if( seqLen > 0 && i + seqLen <= result.size() )
    {
      // Verify continuation bytes (10xxxxxx)
      bool validSeq = true;
      for( size_t j = 1; j < seqLen; ++j )
      {
        unsigned char cont = static_cast<unsigned char>(result[i + j]);
        if( (cont & 0xC0) != 0x80 )
        {
          validSeq = false;
          break;
        }
      }

      if( validSeq )
      {
        // Copy the valid UTF-8 sequence
        for( size_t j = 0; j < seqLen; ++j )
          cleaned += result[i + j];
        i += seqLen - 1;
        continue;
      }
    }

    // Invalid UTF-8 sequence - skip this byte and continue
    // (nlohmann::json will handle it or we'll get a parse error)
  }

  return cleaned;
}//sanitizeJsonString(...)


/** Attempt to repair structurally incomplete JSON.

 This function tries to make incomplete JSON parseable by:
 - Inserting missing commas between object fields / array elements
 - Closing unclosed string literals
 - Adding missing closing braces/brackets

 The function performs a single-pass scan to track JSON structure state and
 appends missing closing delimiters at the end. It handles:
 - Escape sequences (backslash handling)
 - Nested braces and brackets
 - String delimiters

 @param jsonStr The potentially incomplete JSON string (already sanitized)
 @param repairLog Optional output parameter describing repairs made
 @return Repaired JSON string (may still be unparseable if issues are complex)
 */
static std::string repairIncompleteJson( const std::string &jsonStr,
                                         std::string *repairLog )
{
  if( jsonStr.empty() )
    return jsonStr;

  std::string result = jsonStr;

  // Fix missing commas between object fields / array elements before the
  // closing-delimiter pass, since the latter doesn't care about commas.
  result = fixMissingCommas( result );

  // Track JSON structure state
  bool inString = false;
  bool escaped = false;
  std::stack<char> structureStack;  // Track { and [ characters

  // Scan through the string to track state
  for( size_t i = 0; i < jsonStr.length(); ++i )
  {
    const char c = jsonStr[i];

    if( escaped )
    {
      // Previous character was backslash, so this character is escaped
      escaped = false;
      continue;
    }

    if( c == '\\' )
    {
      // Mark next character as escaped
      escaped = true;
      continue;
    }

    if( c == '"' )
    {
      // Toggle string state (only if not escaped)
      inString = !inString;
      continue;
    }

    // Only track structure characters when not inside a string
    if( !inString )
    {
      if( c == '{' || c == '[' )
      {
        structureStack.push( c );
      }
      else if( c == '}' )
      {
        // Pop matching opening brace if present
        if( !structureStack.empty() && structureStack.top() == '{' )
          structureStack.pop();
      }
      else if( c == ']' )
      {
        // Pop matching opening bracket if present
        if( !structureStack.empty() && structureStack.top() == '[' )
          structureStack.pop();
      }
    }
  }

  // Now append missing closing delimiters
  std::vector<std::string> repairs;

  // Close unclosed string first
  if( inString )
  {
    result += '"';
    repairs.push_back( "Added missing closing quote" );
  }

  // Close unclosed structures (in reverse order - LIFO)
  int braceCount = 0;
  int bracketCount = 0;

  while( !structureStack.empty() )
  {
    const char opener = structureStack.top();
    structureStack.pop();

    if( opener == '{' )
    {
      result += '}';
      braceCount++;
    }
    else if( opener == '[' )
    {
      result += ']';
      bracketCount++;
    }
  }

  if( braceCount > 0 )
  {
    repairs.push_back( std::to_string(braceCount) + " closing brace(s)" );
  }

  if( bracketCount > 0 )
  {
    repairs.push_back( std::to_string(bracketCount) + " closing bracket(s)" );
  }

  // Build repair log if requested
  if( repairLog && !repairs.empty() )
  {
    *repairLog = "Added: ";
    for( size_t i = 0; i < repairs.size(); ++i )
    {
      if( i > 0 )
        *repairLog += ", ";
      *repairLog += repairs[i];
    }
  }

  return result;
}//repairIncompleteJson(...)


/** Parse JSON string with lenient preprocessing.

 This function sanitizes the input string to handle common formatting issues
 (whitespace, BOM, invalid UTF-8, null bytes) before parsing. This makes parsing
 more robust when dealing with potentially malformed JSON from external sources.

 @param jsonStr The JSON string to parse
 @return Parsed JSON object
 @throws nlohmann::json::parse_error if parsing fails after sanitization
 */
static nlohmann::json parseLenientJson( const std::string &jsonStr )
{
  const std::string sanitized = sanitizeJsonString( jsonStr );
  return nlohmann::json::parse( sanitized );
}//parseLenientJson(...)


static nlohmann::json lenientlyParseJson( const std::string &jsonStr )
{
  // Phase 0: try as original string - this will catch most input
  try
  {
    return nlohmann::json::parse( jsonStr );
  }catch( std::exception & )
  {
    
  }
  
  // Phase 1: Sanitize the string
  const std::string sanitized = sanitizeJsonString( jsonStr );

  // Phase 2: Try parsing sanitized JSON (fast path)
  try
  {
    return nlohmann::json::parse( sanitized );
  }catch( const nlohmann::json::parse_error &e )
  {
    // Phase 3: If the trimmed string starts with '{' or '[', strip any trailing
    // junk after the matching closing bracket (e.g. " <|call|>" appended by some LLMs).
    // Only attempt this when the top-level structure is an object or array so we
    // don't accidentally truncate a plain-string value.
    const size_t firstNonSpace = sanitized.find_first_not_of( " \t\r\n" );
    if( (firstNonSpace != std::string::npos)
        && ((sanitized[firstNonSpace] == '{') || (sanitized[firstNonSpace] == '[')) )
    {
      const char openChar  = sanitized[firstNonSpace];
      const char closeChar = (openChar == '{') ? '}' : ']';
      size_t depth = 0;
      bool inString = false;
      bool escaped = false;
      size_t closePos = std::string::npos;
      for( size_t i = firstNonSpace; i < sanitized.size(); ++i )
      {
        const char c = sanitized[i];
        if( escaped )
        {
          escaped = false;
          continue;
        }
        if( inString )
        {
          if( c == '\\' )
            escaped = true;
          else if( c == '"' )
            inString = false;
          continue;
        }
        if( c == '"' )
        {
          inString = true;
        }
        else if( c == openChar )
        {
          ++depth;
        }
        else if( c == closeChar )
        {
          --depth;
          if( depth == 0 )
          {
            closePos = i;
            break;
          }
        }
      } //for( loop over sanitized chars )

      if( closePos != std::string::npos && closePos + 1 < sanitized.size() )
      {
        // There is trailing content after the closing bracket - try truncating it.
        const std::string truncated = sanitized.substr( firstNonSpace, closePos - firstNonSpace + 1 );
        try
        {
          return nlohmann::json::parse( truncated );
        }
        catch( ... )
        {
          // Truncation didn't help; fall through to repair logic.
        }
      }
    } //if top-level object or array

    // Phase 5: Attempt structural repair
    std::string repairLog;
    const std::string repaired = repairIncompleteJson( sanitized, &repairLog );

#if( PERFORM_DEVELOPER_CHECKS )
    if( !repairLog.empty() )
    {
      const string msg = "=== JSON Repair Attempted ===\n"
      "Original: " + sanitized.substr( 0, 200 )
      + string(sanitized.length() > 200 ? "...\n" : "\n")
      + "Repairs: " + repairLog + "\n"
      "Repaired: " + repaired.substr( 0, 200 )
      + string(repaired.length() > 200 ? "...\n" : "\n")
      + "=============================";
      
      cout << msg << endl;
      //if( m_debug_stream )
      //  (*m_debug_stream) << msg << endl;
    }
#endif

    // Phase 6: Try parsing repaired JSON
    try
    {
      return nlohmann::json::parse( repaired );
    }
    catch( const nlohmann::json::parse_error &e2 )
    {
      // Repair didn't help - throw original error for better diagnostics
      throw e;
    }
  }
}



/** Manually serialize JSON request with specific field ordering for LLM provider caching.

 Some LLM providers cache parts of requests between calls to reduce costs, but only if
 the beginning of the request string exactly matches previous requests. This function
 ensures that cacheable fields (tools, tool_choice, model, max_completion_tokens, max_tokens)
 appear first and in a reliable order.

 @param requestJson The JSON object to serialize
 @return JSON string with fields in cache-friendly order
 */
/** Helper function to recursively sanitize UTF-8 in JSON values */
static nlohmann::json sanitizeJsonUtf8( const nlohmann::json &value )
{
  if( value.is_string() )
  {
    // Sanitize string values to ensure valid UTF-8
    return sanitizeUtf8( value.get<std::string>() );
  }
  else if( value.is_array() )
  {
    nlohmann::json sanitizedArray = nlohmann::json::array();
    for( const auto &item : value )
      sanitizedArray.push_back( sanitizeJsonUtf8(item) );
    return sanitizedArray;
  }
  else if( value.is_object() )
  {
    nlohmann::json sanitizedObj = nlohmann::json::object();
    for( auto it = value.begin(); it != value.end(); ++it )
      sanitizedObj[it.key()] = sanitizeJsonUtf8( it.value() );
    return sanitizedObj;
  }
  else
  {
    // For numbers, booleans, null - return as-is
    return value;
  }
}

static std::string serializeRequestForCaching( const nlohmann::json &requestJson )
{
  try
  {
    std::string result = "{";
    bool first = true;

    // Helper lambda to add a field if it exists
    auto addField = [&]( const std::string &key )
    {
      if( !requestJson.contains(key) )
        return;

      if( !first )
        result += ",";
      first = false;

      result += "\"" + key + "\":";

      // Sanitize UTF-8 in the value before dumping to prevent JSON encoding errors
      nlohmann::json sanitizedValue = sanitizeJsonUtf8( requestJson[key] );
      result += sanitizedValue.dump();
    };

    // Add fields in priority order for caching (these should stay stable across requests)
    addField( "model" );
    addField( "temperature" );
    addField( "max_completion_tokens" );
    addField( "max_tokens" );
    addField( "reasoning" );
    addField( "reasoning_effort" );
    addField( "tools" );
    addField( "tool_choice" );

    // Add remaining fields (messages, etc.)
    for( auto it = requestJson.begin(); it != requestJson.end(); ++it )
    {
      const std::string &key = it.key();

      // Skip fields we already added in the priority order above
      if( key == "model" || key == "temperature"
         || key == "max_completion_tokens" || key == "max_tokens"
         || key == "reasoning" || key == "reasoning_effort"
         || key == "tools" || key == "tool_choice" )
        continue;

      if( !first )
        result += ",";
      first = false;

      result += "\"" + key + "\":";

      // Sanitize UTF-8 in the value before dumping to prevent JSON encoding errors
      nlohmann::json sanitizedValue = sanitizeJsonUtf8( it.value() );
      result += sanitizedValue.dump();
    }

    result += "}";
    return result;
  }
  catch( const std::exception &e )
  {
    std::cerr << "ERROR in serializeRequestForCaching: " << e.what() << std::endl;
    std::cerr << "  Attempting to return fallback serialization" << std::endl;

    // Fallback: try to return a basic error structure that can be serialized
    nlohmann::json errorJson;
    errorJson["error"] = "Failed to serialize request";
    errorJson["details"] = e.what();
    return errorJson.dump();
  }
}//serializeRequestForCaching(...)


std::string LlmInterface::makeApiCallWithId(const nlohmann::json& requestJson, int requestId)
{
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );

#if( PERFORM_DEVELOPER_CHECKS )
  // Track request string prefix matching for cache verification
  static std::mutex s_cache_check_mutex;
  static std::string s_previous_request_str;
  static nlohmann::json s_previous_request_json;
#endif

  if( m_debug_stream )
  {
    (*m_debug_stream) << "=== Making LLM API Call with ID " << requestId << " ===" << endl;
    (*m_debug_stream) << "Endpoint: " << m_config->llmApi.apiEndpoint() << endl;
    (*m_debug_stream) << "Request JSON:" << endl;
    
    // Wrap debug output in try-catch to prevent debug code from breaking actual functionality
    try
    {
      nlohmann::json debugJson = requestJson;
      
      // Truncate system prompt to first 100 characters, then sanitize UTF-8
      if( debugJson.contains("messages") && debugJson["messages"].is_array() && debugJson["messages"].size() && debugJson["messages"][0].is_object() )
      {
        if( debugJson["messages"][0].contains("content") && debugJson["messages"][0]["content"].is_string() )
        {
          std::string systemPrompt = debugJson["messages"][0]["content"].get<std::string>();
          if( systemPrompt.length() > 100 )
            systemPrompt = systemPrompt.substr(0, 100) + "...";
          systemPrompt = sanitizeUtf8( systemPrompt ); // Sanitize after substring to fix any broken UTF-8
          debugJson["messages"][0]["content"] = systemPrompt;
        }
      }
      
      // Truncate all message contents to 100 characters, then sanitize UTF-8
      if( debugJson.contains("messages") && debugJson["messages"].is_array() )
      {
        for( auto &msg : debugJson["messages"] )
        {
          if( msg.is_object() && msg.contains("content") && msg["content"].is_string() )
          {
            std::string content = msg["content"].get<std::string>();
            if( content.length() > 100 )
              content = content.substr(0, 100) + "...";
            content = sanitizeUtf8( content ); // Sanitize after substring to fix any broken UTF-8
            msg["content"] = content;
          }
        }
      }
      
      if( debugJson.contains("tools") )
      {
        const size_t ntools = debugJson["tools"].size();
        debugJson.erase("tools");
        debugJson["tools"] = std::to_string(ntools) + " tools ommited";
      }
      if( debugJson.contains("tool_choice") )
      {
        const size_t ntools = debugJson["tool_choice"].size();
        debugJson.erase("tool_choice");
        debugJson["tool_choice"] = std::to_string(ntools) + " tool_choice ommited";
      }
      
      (*m_debug_stream) << debugJson.dump(2) << endl;
    }catch( const std::exception &e )
    {
      (*m_debug_stream) << "Error generating debug output (this is non-fatal): " << e.what() << endl;
      (*m_debug_stream) << "  Request will still be sent to LLM API" << endl;
    }
    
    (*m_debug_stream) << "=========================" << endl;
  }//if( m_debug_stream )

  // Serialize request - use incremental caching for MainAgent, full serialization otherwise.
  // We check the messages array to determine if this looks like a MainAgent request (has history),
  // but for sub-agents or summarization requests, we use the stable full serialization.
  // The m_serializedCache is only valid for MainAgent conversations with growing message histories.
  const std::string requestStr = serializeRequestForCaching( requestJson );

#if( PERFORM_DEVELOPER_CHECKS )
  // Developer check: Compare this request with the previous one to predict cache hit rate
  //    Good to check occasionally, but we dont need to clutter up the debug output all the time with it.
  /*
  {
    std::lock_guard<std::mutex> lock( s_cache_check_mutex );

    if( !s_previous_request_str.empty() )
    {
      // String-level prefix comparison
      const size_t min_len = std::min( requestStr.length(), s_previous_request_str.length() );
      size_t matching_chars = 0;
      for( size_t i = 0; i < min_len; ++i )
      {
        if( requestStr[i] != s_previous_request_str[i] )
          break;
        ++matching_chars;
      }

      const double match_percent = (matching_chars * 100.0)
                                 / std::max( requestStr.length(), s_previous_request_str.length() );

      cout << "=== LLM Request Cache Prediction ===" << endl;
      cout << "Current request: " << requestStr.length() << " chars, "
           << "Previous: " << s_previous_request_str.length() << " chars" << endl;
      cout << "Matching prefix: " << matching_chars << " chars ("
           << std::fixed << std::setprecision(1) << match_percent << "%)" << endl;

      // Show where the string-level difference starts
      if( matching_chars < min_len )
      {
        const size_t ctx = std::min( size_t(200), min_len - matching_chars );
        cout << "First difference at char " << matching_chars << ":" << endl;
        cout << "  Current:  \"" << sanitizeUtf8( requestStr.substr( matching_chars, ctx ) ) << "\"" << endl;
        cout << "  Previous: \"" << sanitizeUtf8( s_previous_request_str.substr( matching_chars, ctx ) ) << "\"" << endl;
      }

      // Message-level comparison to identify which message diverges
      if( requestJson.contains("messages") && s_previous_request_json.contains("messages") )
      {
        const auto &curMsgs = requestJson["messages"];
        const auto &prevMsgs = s_previous_request_json["messages"];
        const size_t minMsgs = std::min( curMsgs.size(), prevMsgs.size() );

        for( size_t i = 0; i < minMsgs; ++i )
        {
          if( curMsgs[i].dump() != prevMsgs[i].dump() )
          {
            const string curRole = curMsgs[i].value( "role", "?" );
            const string prevRole = prevMsgs[i].value( "role", "?" );
            cout << "First message diff at index " << i << " of " << minMsgs
                 << " (cur role=\"" << curRole << "\", prev role=\"" << prevRole << "\"):" << endl;

            // Show truncated content of diverging messages
            string curContent = curMsgs[i].dump();
            string prevContent = prevMsgs[i].dump();
            if( curContent.length() > 300 )
              curContent = curContent.substr( 0, 300 ) + "...";
            if( prevContent.length() > 300 )
              prevContent = prevContent.substr( 0, 300 ) + "...";
            cout << "  Current:  " << sanitizeUtf8( curContent ) << endl;
            cout << "  Previous: " << sanitizeUtf8( prevContent ) << endl;
            break;
          }
        }

        if( curMsgs.size() != prevMsgs.size() )
        {
          cout << "Message count: current=" << curMsgs.size()
               << ", previous=" << prevMsgs.size() << endl;
        }
      }
      
#if( BUILD_AS_LOCAL_SERVER )
      // Write debug files if cache hit rate is below threshold
      const double debug_threshold = 40.0; // Write files if cache hit < 40%
      if( match_percent < debug_threshold )
      {
        cout << "Cache hit rate below " << debug_threshold << "% - writing debug files..." << endl;
        writeCacheDebugFiles( requestStr, s_previous_request_str, match_percent );
      }
#endif

      cout << "====================================" << endl;
    }else
    {
      cout << "=== First LLM request (no cache comparison) ===" << endl;
    }

    s_previous_request_str = requestStr;
    s_previous_request_json = requestJson;
  }
   */
#endif // PERFORM_DEVELOPER_CHECKS

  // Use Wt::WWebWidget::jsStringLiteral to properly escape the JSON string for JavaScript
  const std::string jsLiteralRequestStr = Wt::WWebWidget::jsStringLiteral( requestStr );

  string jsCall =
    "var llmRequestData = " + jsLiteralRequestStr + ";\n"
    "window.llmHttpRequest('" + m_config->llmApi.apiEndpoint() + "', llmRequestData, '" +
    m_config->llmApi.bearerToken() + "', " + std::to_string(requestId) + ", '" + m_instanceId + "');";

  // Execute JavaScript to make the HTTP request
  auto app = Wt::WApplication::instance();
  if( !app )
    throw runtime_error( "Error: No WApplication instance available for JavaScript bridge" );

  app->doJavaScript(jsCall);
  
  if( m_debug_stream )
    (*m_debug_stream) << "JavaScript HTTP request with ID " << requestId << " initiated..." << endl;

  return requestStr;
}

void LlmInterface::handleApiResponse( const std::string &response,
                                     const std::shared_ptr<LlmInteraction> &conversation,
                                     const int requestId )
{
  assert( conversation );
  if( !conversation )
  {
    cerr << "LlmInterface::handleApiResponse: received null conversation for response: " << response << endl << endl;
    return;
  }
  
  size_t number_tool_calls = 0;
  shared_ptr<LlmInteractionTurn> request_interaction, response_interaction;
  try
  {
    json responseJson = parseLenientJson( response );
    
    // Parse and accumulate token usage information if available
    std::optional<size_t> promptTokens, completionTokens;
    if( responseJson.contains("usage") )
    {
      const auto &usage = responseJson["usage"];

      std::optional<int> promptTokensInt, completionTokensInt, totalTokens;
      if (usage.contains("prompt_tokens") && usage["prompt_tokens"].is_number())
      {
        promptTokensInt = usage["prompt_tokens"].get<int>();
        promptTokens = static_cast<size_t>(promptTokensInt.value());
      }
      if (usage.contains("completion_tokens") && usage["completion_tokens"].is_number())
      {
        completionTokensInt = usage["completion_tokens"].get<int>();
        completionTokens = static_cast<size_t>(completionTokensInt.value());
      }
      if (usage.contains("total_tokens") && usage["total_tokens"].is_number())
        totalTokens = usage["total_tokens"].get<int>();

      // Accumulate token usage for this conversation
      m_history->addTokenUsage( conversation, promptTokensInt, completionTokensInt, totalTokens );

      if( m_debug_stream )
      {
        const bool has_completion = completionTokens.has_value();
        const bool has_cache_detals = (usage.contains("prompt_tokens_details")
                                       && usage["prompt_tokens_details"].is_object()
                                       && usage["prompt_tokens_details"].contains("cached_tokens")
                                       && usage["prompt_tokens_details"].is_number());
        
        if( has_completion || has_cache_detals )
          (*m_debug_stream) << "=== Token Usage This Call ===" << endl;
        
        if( has_completion )
        {
          (*m_debug_stream) << "Prompt tokens: " << (promptTokens.has_value() ? std::to_string(promptTokens.value()) : "N/A") << endl;
          (*m_debug_stream) << "Completion tokens: " << completionTokens.value() << endl;
          (*m_debug_stream) << "Total tokens: " << (totalTokens.has_value() ? std::to_string(totalTokens.value()) : "N/A") << endl;
        }
        
        // Log prompt cache hit info (OpenRouter and OpenAI both use usage.prompt_tokens_details.cached_tokens)
        if( has_cache_detals )
        {
          const auto &details = usage["prompt_tokens_details"];
          const int cachedTokens = details["cached_tokens"].get<int>();
          (*m_debug_stream) << "Cached tokens: " << cachedTokens;
          if( promptTokens.has_value() && (promptTokens.value() > 0) )
          {
            const double pct = 100.0 * cachedTokens / static_cast<double>( promptTokens.value() );
            (*m_debug_stream) << " (" << static_cast<int>( pct ) << "% of prompt)";
          }
        }//if( has_cache_detals )
        
        if( has_completion || has_cache_detals )
          (*m_debug_stream) << "=============================" << endl;
      }//if( m_debug_stream )
    }//if( responseJson.contains("usage") )
    
    
    string content, reasoning, thinkingContent, thinkingSignature, reasoningContent;
    if( responseJson.contains("choices") && !responseJson["choices"].empty() )
    {
      const json &choice = responseJson["choices"][0];

      if( choice.contains("message") )
      {
        const json &message = choice["message"];
        const string role = message.value("role", "");

        // Handle content field - can be string or array (Claude content blocks)
        if( message.contains("content") )
        {
          if( message["content"].is_string() )
          {
            content = message["content"].get<string>();
          }
          else if( message["content"].is_array() )
          {
            // Claude content block array format with thinking blocks
            for( const auto &block : message["content"] )
            {
              const string blockType = block.value("type", "");
              if( blockType == "thinking" )
              {
                thinkingContent = block.value("thinking", "");
                thinkingSignature = block.value("signature", "");
              }
              else if( blockType == "text" )
              {
                content += block.value("text", "");
              }
              // tool_use blocks are handled separately via tool_calls
            }
          }
        }

        // reasoning_content (some models) at message level (same level as content)
        if( message.contains("reasoning_content") && message["reasoning_content"].is_string() )
        {
          reasoningContent = message["reasoning_content"].get<string>();
          if( thinkingContent.empty() )
            thinkingContent = reasoningContent;  // Also store in thinkingContent for consistency
        }

        // OpenAI o1/o3: reasoning field (older Chat Completions format)
        if( message.contains("reasoning") && message["reasoning"].is_string() )
        {
          reasoning = message["reasoning"].get<string>();
          if( thinkingContent.empty() )
            thinkingContent = reasoning;
        }

        // OpenRouter: reasoning_details array (must be passed back unmodified)
        string reasoningDetails;
        if( message.contains("reasoning_details") && message["reasoning_details"].is_array() )
        {
          // Store the entire array as JSON string to preserve exact structure
          reasoningDetails = message["reasoning_details"].dump();

          // Also extract text content for display purposes
          if( thinkingContent.empty() )
          {
            for( const auto &detail : message["reasoning_details"] )
            {
              if( !detail.is_object() )
                continue;

              string detailType;
              if( detail.contains("type") && detail["type"].is_string() )
                detailType = detail["type"].get<string>();

              if( detailType == "reasoning.text" )
              {
                if( detail.contains("text") && detail["text"].is_string() )
                {
                  const string text = detail["text"].get<string>();
                  if( !text.empty() )
                  {
                    if( !thinkingContent.empty() )
                      thinkingContent += "\n";
                    thinkingContent += text;
                  }
                }
              }
              else if( detailType == "reasoning.summary" )
              {
                if( detail.contains("summary") && detail["summary"].is_string() )
                {
                  const string summary = detail["summary"].get<string>();
                  if( !summary.empty() )
                  {
                    if( !thinkingContent.empty() )
                      thinkingContent += "\n";
                    thinkingContent += summary;
                  }
                }
              }
            }
          }
        }

        if( role == "assistant" )
        {
          // If we haven't already extracted thinking content from content blocks,
          // extract <think>...</think> tags from string content
          if( thinkingContent.empty() && !content.empty() )
          {
            auto [cleanContent, extractedThinking] = extractThinkingAndContent(content);
            content = cleanContent;
            thinkingContent = extractedThinking;
          }

          // Handle structured tool calls first (OpenAI format)
          pair<shared_ptr<LlmToolRequest>,shared_ptr<LlmToolResults>> call_interactions;
          if (message.contains("tool_calls"))
          {
            call_interactions = executeToolCallsAndSendResults( message["tool_calls"], conversation, requestId, response,
              thinkingContent, thinkingSignature, reasoningContent, reasoningDetails, promptTokens, completionTokens );
          }else
          {
            // Parse content for text-based tool requests (use cleaned content)
            call_interactions = parseContentForToolCallsAndSendResults( content, conversation, requestId, response,
              thinkingContent, thinkingSignature, reasoningContent, reasoningDetails );
          }

          request_interaction = call_interactions.first;
          response_interaction = call_interactions.second;

          if( call_interactions.first )
          {
            number_tool_calls += call_interactions.first->toolCalls().size();
            conversation->consecutiveNonFinalAutoReplies = 0;
          }

          // Add assistant message to history with thinking content - only if there were no tool calls
          if( number_tool_calls == 0 )
          {
            request_interaction = m_history->addAssistantMessageWithThinking( content, thinkingContent, response, conversation );

            // Store additional reasoning fields
            if( request_interaction )
            {
              request_interaction->setThinkingSignature( thinkingSignature );
              request_interaction->setReasoningContent( reasoningContent );
              request_interaction->setReasoningDetails( reasoningDetails );
            }
          }
        }//if( role == "assistant" )
      }//if( choice.contains("message") )
    }

    // Lets check for Ask Sage style tool-calls to the `/server/query` endpoint
    if( (number_tool_calls == 0)
       && responseJson.contains("tool_calls")
       && !responseJson["tool_calls"].is_null()
       && !responseJson["tool_calls"].empty() )
    {
      pair<shared_ptr<LlmToolRequest>,shared_ptr<LlmToolResults>> call_interactions;
      call_interactions = executeToolCallsAndSendResults( responseJson["tool_calls"], conversation, requestId, response,
        "", "", "", "", promptTokens, completionTokens );
      if( call_interactions.first )
      {
        number_tool_calls += call_interactions.first->toolCalls().size();
        conversation->consecutiveNonFinalAutoReplies = 0;
      }
    }//if( responseJson.contains("choices") && !responseJson["choices"].empty() )

    // Lets check for Ask Sage style message to the `/server/query` endpoint
    if( content.empty()
       && responseJson.contains("message")
       && !responseJson["message"].is_null()
       && !responseJson["message"].empty() )
    {
      content = responseJson["message"];
      auto [cleanContent, extractedThinking] = extractThinkingAndContent(content);
      content = cleanContent;
      if( thinkingContent.empty() )
        thinkingContent = extractedThinking;

      if( number_tool_calls == 0 )
        m_history->addAssistantMessageWithThinking( content, thinkingContent, response, conversation );
    }

    // Sometimes a model might need a little prompting to continue....
    SpecUtils::trim( content );
    if( (number_tool_calls == 0) && content.empty() && !thinkingContent.empty() )
    {
      // Check if the last response was already an AutoReply
      bool lastWasAutoReply = false;
      if( !conversation->responses.empty() )
      {
        const std::shared_ptr<LlmInteractionTurn> &lastResponse = conversation->responses.back();
        lastWasAutoReply = (lastResponse->type() == LlmInteractionTurn::Type::AutoReply);
      }

      if( lastWasAutoReply )
      {
        // We already sent an auto-reply prompt and the LLM still didn't respond properly.
        // Just end the conversation instead of looping.
        cout << "LLM provided reasoning but no content/tool calls after auto-reply prompt. Ending conversation." << endl;
        if( m_debug_stream )
          (*m_debug_stream) << "LLM provided reasoning but no content/tool calls after auto-reply prompt. Ending conversation." << endl;
      }else
      {
        // Add an auto-reply message to prompt the LLM to continue
        const string autoReplyContent =
        "Please continue the analysis - and retry any failed tool calls (please re-try liberally) and make any additional tool calls that would be helpful."
        "Then applying careful reasoning at each step, and continuing to to working on the problem until you have"
        " arrived at an answer, so you can provide a summary.";
        cout << "LLM provided reasoning but no content/tool calls. Sending auto-reply prompt." << endl;

        shared_ptr<LlmInteractionAutoReply> autoReply = m_history->addAutoReplyMessage( autoReplyContent, conversation );

        // Build API request with the auto-reply included
        const json requestJson = buildMessagesArray( conversation );

        // Make tracked API call to continue the conversation
        std::pair<int,std::string> request_id_content = makeTrackedApiCall( requestJson, conversation );

        if( m_debug_stream )
          (*m_debug_stream) << "Sent auto-reply prompt with request ID: " << request_id_content.first << endl;

        // Don't finish the conversation yet - we're waiting for the LLM to respond to our prompt
        return;
      }
    }//if( (number_tool_calls == 0) && content.empty() && !reasoning.empty() )

    // For sub-agents with state machines: if the LLM responded with text but no tool calls,
    // and the current state is NOT a final state, auto-reply to prompt it to continue.
    if( (number_tool_calls == 0)
       && !content.empty()
       && (conversation->agent_type != AgentType::MainAgent)
       && conversation->state_machine
       && !conversation->state_machine->isFinalState( conversation->state_machine->getCurrentState() )
       && !conversation->requestedNonFinalSummary )
    {
      const bool consecutiveLimitReached
        = (conversation->consecutiveNonFinalAutoReplies >= sm_max_consecutive_nonfinal_autoreplies);
      const bool totalLimitReached
        = (conversation->totalNonFinalAutoReplies >= sm_max_total_nonfinal_autoreplies);

      if( consecutiveLimitReached || totalLimitReached )
      {
        // Limits exhausted - ask for a summary of progress before finishing
        if( m_debug_stream )
          (*m_debug_stream) << "Sub-agent " << agentTypeToString( conversation->agent_type )
             << " in non-final state '" << conversation->state_machine->getCurrentState()
             << "' but auto-reply limit reached (consecutive="
             << conversation->consecutiveNonFinalAutoReplies
             << ", total=" << conversation->totalNonFinalAutoReplies
             << "). Requesting summary before ending." << endl;

        conversation->requestedNonFinalSummary = true;

        const string summaryRequest =
          "You have not completed the full workflow, but no further retries are available."
          " Please provide a concise summary of: (1) what work was completed and any changes made,"
          " (2) the current state of the analysis, and (3) any remaining steps that were not completed."
          " This summary will be provided to the coordinating agent.";

        shared_ptr<LlmInteractionAutoReply> autoReply
          = m_history->addAutoReplyMessage( summaryRequest, conversation );

        const json requestJson = buildMessagesArray( conversation );
        std::pair<int,std::string> request_id_content
          = makeTrackedApiCall( requestJson, conversation );

        if( m_debug_stream )
          (*m_debug_stream) << "Sent non-final summary request with request ID: "
             << request_id_content.first << endl;

        return;
      }else
      {
        // Send a re-prompt to nudge the LLM to continue
        const string &currentState = conversation->state_machine->getCurrentState();
        const vector<string> allowedTransitions = conversation->state_machine->getAllowedTransitions();
        const string guidance = conversation->state_machine->getPromptGuidanceForCurrentState();

        string transitionList;
        for( size_t i = 0; i < allowedTransitions.size(); ++i )
        {
          if( i > 0 )
            transitionList += ", ";
          transitionList += allowedTransitions[i];
        }

        // Scan response content and thinking/reasoning text for tool names
        // the model may have intended to call
        const map<string, LlmTools::SharedTool> agentTools
          = m_tool_registry->getToolsForAgent( conversation->agent_type );
        vector<string> mentionedTools;
        for( const auto &[name, tool] : agentTools )
        {
          if( (content.find( name ) != string::npos)
             || (!thinkingContent.empty() && (thinkingContent.find( name ) != string::npos))
             || (!reasoningContent.empty() && (reasoningContent.find( name ) != string::npos)) )
          {
            mentionedTools.push_back( name );
          }
        }

        string autoReplyContent;
        if( !mentionedTools.empty() )
        {
          autoReplyContent = "A tool call is expected, but your previous response contained"
            " only text. It looks like you were considering calling: ";
          for( size_t i = 0; i < mentionedTools.size(); ++i )
          {
            if( i > 0 )
              autoReplyContent += ", ";
            autoReplyContent += "`" + mentionedTools[i] + "`";
          }
          autoReplyContent += ". Please make that tool call now if it is what you intended.";
        }else
        {
          autoReplyContent = "You are currently in workflow state '" + currentState
            + "', which is NOT a final state. You must continue working - please make"
            " the necessary tool calls to complete the current phase, then call"
            " `set_workflow_state` to transition to the next state.";
        }

        if( !transitionList.empty() )
          autoReplyContent += " Allowed transitions from this state: " + transitionList + ".";

        if( !guidance.empty() )
          autoReplyContent += " Guidance for current state: " + guidance;

        autoReplyContent += " Do not provide a final summary until you have reached a"
          " final state in the workflow.";

        autoReplyContent += "\n\nPlease continue making tool calls and reasoning for this current state, or explicitly call the `set_workflow_state` tool-call, with an allowed transition state, to transition to the next state.";

        if( m_debug_stream )
          (*m_debug_stream) << "Sub-agent " << agentTypeToString( conversation->agent_type )
             << " in non-final state '" << currentState
             << "'. Sending non-final-state auto-reply (consecutive="
             << (conversation->consecutiveNonFinalAutoReplies + 1)
             << ", total=" << (conversation->totalNonFinalAutoReplies + 1)
             << ")." << endl;

        conversation->consecutiveNonFinalAutoReplies++;
        conversation->totalNonFinalAutoReplies++;

        shared_ptr<LlmInteractionAutoReply> autoReply
          = m_history->addAutoReplyMessage( autoReplyContent, conversation );

        const json requestJson = buildMessagesArray( conversation );
        std::pair<int,std::string> request_id_content
          = makeTrackedApiCall( requestJson, conversation );

        if( m_debug_stream )
          (*m_debug_stream) << "Sent non-final-state auto-reply with request ID: "
             << request_id_content.first << endl;

        return;
      }
    }//if( sub-agent, no tool calls, non-final state )
  }catch( const std::exception &e )
  {
    cerr << "Error parsing LLM response: " << e.what() << "\n\tresponse=" << response << endl << endl;
    if( m_debug_stream )
      (*m_debug_stream) << "Error parsing LLM response: " << e.what() << "\n\tresponse=" << response << endl << endl;
  }

  if( !number_tool_calls )
  {
    conversation->finishTime = std::chrono::system_clock::now();
    if( conversation->conversation_completion_handler )
      conversation->conversation_completion_handler();
    conversation->conversationFinished.emit();
  }
  
  // Only emit signal if there are no pending requests AND no deferred async tools still running.
  // This is the SUCCESS path - emit m_conversationFinished (errors emit m_responseError)
  if( m_pendingRequests.empty() && m_deferredToolResults.empty() )
  {
    if( m_debug_stream )
      (*m_debug_stream) << "No pending requests or deferred tools, emitting response received signal" << endl;
    m_conversationFinished.emit();
  }else
  {
    if( m_debug_stream )
      (*m_debug_stream) << "Still have " << m_pendingRequests.size() << " pending requests"
        << " and " << m_deferredToolResults.size() << " deferred tool results, not emitting signal yet" << endl;
  }
}//void handleApiResponse(...)


std::pair<std::shared_ptr<LlmToolRequest>, std::shared_ptr<LlmToolResults>>
LlmInterface::executeToolCallsAndSendResults( const nlohmann::json &toolCalls,
                                    const std::shared_ptr<LlmInteraction> &convo,
                                    const int parentRequestId,
                                    const std::string &rawResponseContent,
                                    const std::string &thinkingContent,
                                    const std::string &thinkingSignature,
                                    const std::string &reasoningContent,
                                    const std::string &reasoningDetails,
                                    std::optional<size_t> promptTokens,
                                    std::optional<size_t> completionTokens )
{
  assert( wApp ); //Consistency requires we have the WApplication::UpdateLock
  
  assert( convo );
  if( !convo )
  {
    cerr << "LlmInterface::executeToolCallsAndSendResults: recieved null conversation for tool calls: " << toolCalls.dump(2) << endl << endl;
    return {nullptr, nullptr};
  }

  // Set current conversation for tool execution context
  m_currentConversation = convo;

  if( m_debug_stream )
    (*m_debug_stream) << "Executing " << toolCalls.size() << " tool calls" << endl;

  // Track executed tool calls for follow-up
  std::vector<std::string> executedToolCallIds;
  vector<int> subAgentRequestIds;
  size_t asyncToolCount = 0;  // Number of async tool calls dispatched to background threads

  // Collect all tool call requests for batching
  std::vector<LlmToolCall> toolCallRequests;
  std::vector<LlmToolCall> toolCallResults;

  for( const nlohmann::json &toolCall : toolCalls )
  {
    const string callId = toolCall.value("id", "");
    string toolName;
    string rawArguments;  // Keep raw string before parsing

    try
    {
      if( !toolCall.contains("function") )
        throw runtime_error( "No 'function' definition in the JSON" );

      const nlohmann::json &function_def = toolCall["function"];

      if( !function_def.contains("name") )
        throw runtime_error( "No tool 'name' given definition in the JSON" );
      toolName = function_def["name"];
      
      // Sanitize tool name to remove any spurious characters like angle brackets
      // (some LLM models return things like "tool_name<|channel|>")
      toolName = sanitizeToolName( toolName );

      // Get raw arguments string (don't parse yet - we'll parse based on tool type)
      if( function_def.contains("arguments") )
        rawArguments = function_def["arguments"].get<string>();

      if( m_debug_stream )
        (*m_debug_stream) << "Calling tool: " << toolName << " with ID: " << callId << endl;

      // Check if this is an invoke_* tool (sub-agent invocation) - handle specially
      if( SpecUtils::istarts_with(toolName, "invoke_") && (toolName.length() > 7) )
      {
        // ===== SUB-AGENT INVOCATION BRANCH =====
        // For sub-agents: Extract context/task robustly, don't require perfect JSON

        string context, task;
        json arguments;  // Will be constructed after extraction

        // Try to parse as JSON and extract context/task
        try
        {
          if( m_debug_stream )
            (*m_debug_stream) << "--- about to parse sub-agent arguments (lenient)" << endl;
          
          try
          {
            arguments = lenientlyParseJson( rawArguments );
          }catch( ...)
          {
            if( m_debug_stream )
              (*m_debug_stream) << "--- Failed to parse sub-agent arguments (lenient) of:\n```" << rawArguments << "\n```" << endl;
            throw;
          }

          if( m_debug_stream )
            (*m_debug_stream) << "--- done parsing sub-agent arguments" << endl;
          
          // Try to extract context and task fields
          if( arguments.contains("context") && arguments["context"].is_string() )
            context = arguments["context"].get<string>();

          if( arguments.contains("task") && arguments["task"].is_string() )
            task = arguments["task"].get<string>();

          // If we didn't get both fields, try to extract any string values from malformed JSON
          if( context.empty() && task.empty() && arguments.is_object() )
          {
            cerr << "Warning: Sub-agent arguments missing context/task fields, extracting string values" << endl;
            vector<string> stringValues;
            for( auto it = arguments.begin(); it != arguments.end(); ++it )
            {
              if( it.value().is_string() )
                stringValues.push_back( it.value().get<string>() );
            }

            // Assume first string is context, second is task
            if( stringValues.size() >= 1 )
              context = stringValues[0];
            if( stringValues.size() >= 2 )
              task = stringValues[1];
            else if( stringValues.size() == 1 )
              task = stringValues[0];  // Use single string as task
          }
        }
        catch( const std::exception &e )
        {
          // Complete JSON parsing failure - use raw string as task
          cerr << "Failed to parse sub-agent arguments as JSON: " << e.what() << endl;
          cerr << "Using raw arguments string as task" << endl;
          task = rawArguments;
        }

        // If we still don't have any content, use raw string as task
        if( context.empty() && task.empty() )
        {
          cerr << "Warning: Sub-agent invocation with empty context and task, using raw arguments" << endl;
          task = rawArguments;
        }

        // Reconstruct arguments JSON for consistency
        arguments = json::object();
        arguments["context"] = context;
        arguments["task"] = task;

        // Add to requests BEFORE creating sub-agent (so it's always in sync with results)
        toolCallRequests.emplace_back( toolName, callId, arguments );

        // Some LLM providers require the "extra_content" that accompanied the tool-call
        // request to be echoed back with the tool result.  Sub-agent invocations are just
        // tool calls from the LLMs perspective, so preserve it here exactly like the normal
        // tool branch does, so it gets serialized back when sending results to the LLM.
        if( toolCall.contains("extra_content") )
          toolCallRequests.back().extra_content = toolCall["extra_content"];

        // Extract agent name from tool name (e.g., "invoke_NuclideId" -> "NuclideId")
        const string agentName = toolName.substr(7);
        const AgentType agent_type = stringToAgentType(agentName);
        assert( agent_type != AgentType::MainAgent );
        if( agent_type == AgentType::MainAgent )
          throw runtime_error( "Can not invoke a agent of type AgentType::MainAgent" );

        if( m_debug_stream )
        {
          (*m_debug_stream) << "Detected sub-agent invocation: " << agentName << " (context: " << context.substr(0,50)
                            << (context.length() > 50 ? "..." : "")
                            << ", task: " << task.substr(0,50) << (task.length() > 50 ? "..." : "") << ")" << endl;
        }

        // Invoke sub-agent (returns request ID)
        const string combinedMessage = "Context:\n" + context + "\n\nTask:\n" + task + "\n\nWhen you are done, provide a summary of your reasoning, as well as how sure you are, or alternate possibilities to investigate.";
        shared_ptr<LlmInteraction> sub_agent_convo = LlmInteraction::createAgent( convo->type, combinedMessage, agent_type, m_config.get() );
        // finishTime will be set when the conversation completes
        
        // TODO: make sure we dont have two sub-agent calls for the same type of agent - among other problems, this would make `sub_agent_convo->conversationId` non-unique
        const auto current_ticks = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
        sub_agent_convo->conversationId = "subagent_" + agentName + "_" + std::to_string( current_ticks );
        
        const std::weak_ptr<LlmInteraction> &parent_convo_wk = convo;
        const std::weak_ptr<LlmInteraction> sub_agent_convo_wk = sub_agent_convo;
        
        const int subAgentRequestId = invokeSubAgent( sub_agent_convo );
        subAgentRequestIds.push_back( subAgentRequestId );

        {// Begin placeholder result that will be updated when sub-agent completes
          json result;
          result["status"] = "pending";
          result["message"] = "Sub-agent " + agentName + " is processing (will update when complete)";
          result["requestId"] = subAgentRequestId;

          // Create tool result entry with placeholder
          LlmToolCall toolResult( toolName, callId, arguments );
          toolResult.status = LlmToolCall::CallStatus::Pending;
          toolResult.content = result.dump();
          toolResult.sub_agent_conversation = sub_agent_convo;

          toolCallResults.push_back( std::move(toolResult) );
        }// End placeholder result that will be updated when sub-agent completes
        
        // We will define a completion handler that will get called when the sub-agent conversation is complete.
        LlmInterface * self = this;
        sub_agent_convo->conversation_completion_handler = [parent_convo_wk, parentRequestId, subAgentRequestId, self, sub_agent_convo_wk, callId](){
          shared_ptr<LlmInteraction> parent_conv = parent_convo_wk.lock();
          assert( parent_conv );
          if( !parent_conv )
          {
            cerr << "Sub-agent completed with parent conversation no longer alive - ignoring results." << endl;
            return;
          }
          
          shared_ptr<LlmInteraction> sub_agent_convo = sub_agent_convo_wk.lock();
          assert( sub_agent_convo );
          if( !sub_agent_convo )
          {
            cerr << "Sub-agent cconversation no longer alive - ignoring results." << endl;
            return;
          }
          
          InterSpec *viewer = InterSpec::instance();
          if( !viewer )
          {
            cerr << "Sub-agent completed for dead session - ignoring results." << endl;
            return;
          }
          
          // TODO: maybe `LlmInterface` should use the Wt object life and/or boost signal/slot mechanism to protect against LlmInterface lifetime - and not this hack of a system of relying on checking stuff through the GUI
          LlmToolGui * const llm_gui = viewer->currentLlmTool();
          if( !llm_gui )
          {
            cerr << "Sub-agent completed and no llm_gui avaiable - ignoring results." << endl;
            return;
          }
          
          LlmInterface * const interface = llm_gui->llmInterface();
          if( interface != self )
          {
            cerr << "Sub-agent completed and with different LlmInterface now present - ignoring results." << endl;
            return;
          }
          
          // Here is where we fill in `parent_conv` with the result, and send back to the LLM
          if( interface->m_debug_stream )
            (*interface->m_debug_stream) << "Sub-agent complete (no tool calls left) - will extracting summary and continue main conversation." << endl;

          assert( !parent_conv->responses.empty() );
          std::shared_ptr<LlmInteractionTurn> response = nullptr;
          LlmToolCall *toolResult = nullptr;

          // Find the ToolResult response containing this tool call
          for( size_t rspns_index = 0; !toolResult && (rspns_index < parent_conv->responses.size()); ++rspns_index )
          {
            if( parent_conv->responses[rspns_index]->type() == LlmInteractionTurn::Type::ToolResult )
            {
              response = parent_conv->responses[rspns_index];
              // Search within the toolCalls vector for the matching invocationId - need to cast to access toolCalls
              LlmToolResults *toolResResp = dynamic_cast<LlmToolResults*>(response.get());
              if( toolResResp )
              {
                for( LlmToolCall &tc : toolResResp->toolCalls() )
                {
                  if( tc.invocationId == callId )
                  {
                    toolResult = &tc;
                    assert( tc.content.find( "(will update when complete)" ) != string::npos );
                    break;
                  }
                }
              }
            }
          }//for( loop over responses )

          assert( toolResult );

          if( toolResult )
          {
            assert( !sub_agent_convo->responses.empty() );
            
            json updatedResult;
            
            string subAgentSummary, subAgentThinkingContent;
            if( sub_agent_convo->responses.empty() )
            {
              updatedResult["status"] = "failed";
              updatedResult["summary"] = "Failed to get agent summary - because of internal logic error.";
              updatedResult["agentName"] = agentTypeToString( sub_agent_convo->agent_type );
            }else
            {
              string last_response;
              // Extract content from the last response (should be an Assistant response)
              const LlmInteractionFinalResponse *llmResp = dynamic_cast<const LlmInteractionFinalResponse*>(sub_agent_convo->responses.back().get());
              if( llmResp )
                last_response = llmResp->content();
              if( interface->m_debug_stream )
                (*interface->m_debug_stream) << "--- about to parse sub_agent_convo->responses='" << last_response << "'" << endl;

              nlohmann::json responseJson = nlohmann::json::object();
              try
              {
                if( !last_response.empty() )
                  responseJson = parseLenientJson( last_response );
                if( interface->m_debug_stream )
                  (*interface->m_debug_stream) << "--- Done parsing sub_agent_convo->responses" << endl;
                
                if( responseJson.contains("choices") && !responseJson["choices"].empty() )
                {
                  const json &choice = responseJson["choices"][0];
                  if( choice.contains("message") && choice["message"].contains("content") )
                  {
                    const string content = choice["message"]["content"].get<string>();
                    // Extract clean content (strip thinking tags if present)
                    auto [cleanContent, thinkingContent] = extractThinkingAndContent(content);
                    subAgentSummary = cleanContent;
                    subAgentThinkingContent = thinkingContent;
                  }
                }
                
                if( interface->m_debug_stream )
                  (*interface->m_debug_stream) << "--- Done extracting JSON subAgentSummary='" << subAgentSummary << "'\n\n" << endl;
              }catch( std::exception &e )
              {
                // `LlmConversationResponse::content` may not / likely is not, JSON
                auto [cleanContent, thinkingContent] = extractThinkingAndContent(last_response);
                subAgentSummary = cleanContent;
                subAgentThinkingContent = thinkingContent;
                if( interface->m_debug_stream )
                  (*interface->m_debug_stream) << "--- Done extracting non-JSON subAgentSummary='" << subAgentSummary << "'\n\n" << endl;
              }
              
              
              // Sanitize UTF-8 before adding to JSON (nlohmann::json requires valid UTF-8)
              subAgentSummary = sanitizeUtf8( subAgentSummary );
              subAgentThinkingContent = sanitizeUtf8( subAgentThinkingContent );

              if( sub_agent_convo->agent_type == AgentType::DeepResearch )
              {
                std::string trimmed_summary = subAgentSummary;
                SpecUtils::trim( trimmed_summary );
                if( trimmed_summary.empty() )
                  subAgentSummary = "no further information found";
              }

              if( interface->m_debug_stream )
              {
                (*interface->m_debug_stream) << "Sub-agent summary: " << subAgentSummary.substr(0, 100) << (subAgentSummary.length() > 100 ? "..." : "") << endl;
                (*interface->m_debug_stream) << "subAgentThinkingContent='" << subAgentThinkingContent << "'" << endl;
              }
              
              updatedResult["status"] = "completed";
              updatedResult["summary"] = subAgentSummary;
              updatedResult["agentName"] = agentTypeToString( sub_agent_convo->agent_type );
            }//if( sub_agent_convo->responses.empty() )

            toolResult->status = LlmToolCall::CallStatus::Success;
            
            // Update the tool result content
            toolResult->content = updatedResult.dump();

            // Store sub-agent conversation reference in the tool call request
            toolResult->sub_agent_conversation = sub_agent_convo;

            // Store thinking content in the parent response
            if( response )
              response->setThinkingContent( subAgentThinkingContent );

            // Aggregate sub-agent token usage to parent conversation
            if( sub_agent_convo->promptTokens.has_value() ||
                sub_agent_convo->completionTokens.has_value() ||
                sub_agent_convo->totalTokens.has_value() )
            {
              if( interface->m_debug_stream )
              {
                (*interface->m_debug_stream) << "=== Aggregating sub-agent token usage to parent ===" << endl;
                (*interface->m_debug_stream) << "Sub-agent prompt tokens: " << (sub_agent_convo->promptTokens.has_value() ? std::to_string(sub_agent_convo->promptTokens.value()) : "N/A") << endl;
                (*interface->m_debug_stream) << "Sub-agent completion tokens: " << (sub_agent_convo->completionTokens.has_value() ? std::to_string(sub_agent_convo->completionTokens.value()) : "N/A") << endl;
                (*interface->m_debug_stream) << "Sub-agent total tokens: " << (sub_agent_convo->totalTokens.has_value() ? std::to_string(sub_agent_convo->totalTokens.value()) : "N/A") << endl;
              }

              // Add sub-agent's tokens to parent conversation
              if( sub_agent_convo->promptTokens.has_value() )
              {
                if( parent_conv->promptTokens.has_value() )
                  parent_conv->promptTokens = parent_conv->promptTokens.value() + sub_agent_convo->promptTokens.value();
                else
                  parent_conv->promptTokens = sub_agent_convo->promptTokens.value();
              }

              if( sub_agent_convo->completionTokens.has_value() )
              {
                if( parent_conv->completionTokens.has_value() )
                  parent_conv->completionTokens = parent_conv->completionTokens.value() + sub_agent_convo->completionTokens.value();
                else
                  parent_conv->completionTokens = sub_agent_convo->completionTokens.value();
              }

              if( sub_agent_convo->totalTokens.has_value() )
              {
                if( parent_conv->totalTokens.has_value() )
                  parent_conv->totalTokens = parent_conv->totalTokens.value() + sub_agent_convo->totalTokens.value();
                else
                  parent_conv->totalTokens = sub_agent_convo->totalTokens.value();
              }

              if( interface->m_debug_stream )
              {
                (*interface->m_debug_stream) << "Parent conversation now has:" << endl;
                (*interface->m_debug_stream) << "  Prompt tokens: " << (parent_conv->promptTokens.has_value() ? std::to_string(parent_conv->promptTokens.value()) : "N/A") << endl;
                (*interface->m_debug_stream) << "  Completion tokens: " << (parent_conv->completionTokens.has_value() ? std::to_string(parent_conv->completionTokens.value()) : "N/A") << endl;
                (*interface->m_debug_stream) << "  Total tokens: " << (parent_conv->totalTokens.has_value() ? std::to_string(parent_conv->totalTokens.value()) : "N/A") << endl;
                (*interface->m_debug_stream) << "===================================================" << endl;
              }
            }

            // Emit signal that this sub-agent has finished
            parent_conv->subAgentFinished.emit( response, sub_agent_convo );
          }else
          {
            cerr << "Failed to find tool-call result - shouldnt happen!!!" << endl;
            if( interface->m_debug_stream )
              (*interface->m_debug_stream) << "Failed to find tool-call result - shouldnt happen!!!" << endl;
          }//if( toolResult )


          const auto defered_pos = interface->m_deferredToolResults.find(parentRequestId);
          assert( defered_pos != end(interface->m_deferredToolResults) );
          if( defered_pos == end(interface->m_deferredToolResults) )
          {
            if( interface->m_debug_stream )
              (*interface->m_debug_stream) << "Sub-agent completed with no defered request being found..." << endl;
            return;
          }
          
          
          DeferredToolResult &deffered_result = defered_pos->second;
          vector<int> &subAgentToolCallIds = deffered_result.subAgentToolCallIds;
          const auto pos = std::find( begin(subAgentToolCallIds), end(subAgentToolCallIds), subAgentRequestId );
          assert( pos != end(subAgentToolCallIds) );
          subAgentToolCallIds.erase( pos );
          if( subAgentToolCallIds.empty() && (deffered_result.pendingAsyncToolCount == 0) )
          {
            // No more agent calls or async tools pending, lets cleanup the results
            interface->m_deferredToolResults.erase( defered_pos );
            if( interface->m_debug_stream )
              (*interface->m_debug_stream) << "=== Sending parent conversation with sub-agent results back to LLM!" << endl;
            interface->sendToolResultsToLLM( parent_conv );
          }
        };//sub_agent_convo->conversation_completion_handler

        if( m_debug_stream )
          (*m_debug_stream) << "Sub-agent invocation deferred - will resume main agent when complete" << endl;
      }else
      {
        // ===== NORMAL TOOL EXECUTION BRANCH =====
        // For normal tools: Parse arguments as JSON (strict)

        // Add to requests BEFORE parsing arguments, so the request entry always
        // exists even if lenientlyParseJson throws (keeping requests/results in sync)
        toolCallRequests.emplace_back( toolName, callId, json::object() );

        if( m_debug_stream )
          (*m_debug_stream) << "--- about to parse normal tool arguments for toolName=" << toolName << endl;
        
        json arguments;
        if( !rawArguments.empty() )
          arguments = lenientlyParseJson( rawArguments );
        
        if( m_debug_stream )
          (*m_debug_stream) << "--- done parsing normal tool arguments for toolName=" << toolName << endl;

        // Update the stored request with the parsed arguments
        toolCallRequests.back().toolParameters = arguments;

        if( toolCall.contains("extra_content") )
          toolCallRequests.back().extra_content = toolCall["extra_content"];

        // Create tool result entry
        LlmToolCall toolResult( toolName, callId, arguments );

        if( m_block_tool_calls )
        {
          // Tool calls are blocked; return an error result without executing the tool
          toolResult.status = LlmToolCall::CallStatus::Error;
          toolResult.content = R"({"status": false, "error": "Tool calls are disabled now, and for the future for this conversation - please answer using only the tool calls that were previously made."})";
          toolCallResults.push_back( std::move(toolResult) );
        }else
        {
          const LlmTools::SharedTool *tool = m_tool_registry->getTool( toolName );

          if( tool && tool->isAsync() )
          {
            // === ASYNC TOOL PATH ===
            // Tool runs expensive computation in background; callback fires on GUI thread.
            toolResult.status = LlmToolCall::CallStatus::Pending;
            toolResult.content = R"({"status": "pending"})";
            toolCallResults.push_back( std::move(toolResult) );

            asyncToolCount++;

            const string capturedCallId = callId;
            const string capturedToolName = toolName;
            const int capturedParentRequestId = parentRequestId;
            weak_ptr<LlmInteraction> weakConvo = convo;
            LlmInterface *self = this;

            const auto exec_start = std::chrono::steady_clock::now();

            if( m_debug_stream )
              (*m_debug_stream) << "Dispatching async tool: " << toolName << " (callId=" << callId << ")" << endl;

            try
            {
            tool->asyncExecutor( arguments, m_interspec, convo, m_history.get(),
              [self, weakConvo, capturedCallId, capturedToolName, capturedParentRequestId, exec_start]
              ( std::variant<nlohmann::json, std::string> result_or_error )
              {
                // This callback runs on the GUI thread (posted by the async function)
                shared_ptr<LlmInteraction> convo = weakConvo.lock();
                if( !convo )
                  return;

                InterSpec *viewer = InterSpec::instance();
                if( !wApp || !viewer )
                  return;

                const auto exec_end = std::chrono::steady_clock::now();
                const auto exec_duration = std::chrono::duration_cast<std::chrono::milliseconds>( exec_end - exec_start );

                if( self->m_debug_stream )
                  (*self->m_debug_stream) << "Async tool callback for: " << capturedToolName
                    << " (callId=" << capturedCallId << ", duration=" << exec_duration.count() << "ms)" << endl;

                // Find the Pending tool result in the conversation and update it
                LlmToolCall *pendingResult = nullptr;
                LlmToolResults *owningToolResults = nullptr;
                for( shared_ptr<LlmInteractionTurn> &resp : convo->responses )
                {
                  if( resp->type() != LlmInteractionTurn::Type::ToolResult )
                    continue;

                  LlmToolResults *toolResResp = dynamic_cast<LlmToolResults *>( resp.get() );
                  if( !toolResResp )
                    continue;

                  for( LlmToolCall &tc : toolResResp->toolCalls() )
                  {
                    if( tc.invocationId == capturedCallId )
                    {
                      pendingResult = &tc;
                      owningToolResults = toolResResp;
                      break;
                    }
                  }
                  if( pendingResult )
                    break;
                }//for( each response )

                if( !pendingResult )
                {
                  if( self->m_debug_stream )
                    (*self->m_debug_stream) << "Async callback: could not find pending result for callId="
                      << capturedCallId << " - discarding" << endl;
                  return;
                }

                if( const string *error = std::get_if<string>( &result_or_error ) )
                {
                  json errJson;
                  errJson["error"] = "Tool call failed: " + *error;
                  errJson["success"] = false;
                  pendingResult->status = LlmToolCall::CallStatus::Error;
                  pendingResult->content = errJson.dump();
                }else
                {
                  json &result = std::get<json>( result_or_error );

                  // Normalize result to JSON object
                  if( result.type() != nlohmann::detail::value_t::object )
                  {
                    json result_copy = result;
                    result = json::object();
                    result["value"] = std::move( result_copy );
                  }

                  if( !result.contains("success") && !result.contains("Success") )
                    result["success"] = true;

                  pendingResult->status = LlmToolCall::CallStatus::Success;

                  // Extract image content if present
                  if( result.contains( "_image" ) )
                  {
                    const json &img = result["_image"];
                    LlmToolCall::ImageContent ic;
                    ic.base64Data = img.value( "base64Data", "" );
                    ic.mimeType = img.value( "mimeType", "" );
                    ic.widthPx = img.value( "widthPx", 0 );
                    ic.heightPx = img.value( "heightPx", 0 );
                    if( !ic.base64Data.empty() )
                      pendingResult->imageContent = std::move( ic );
                    result.erase( "_image" );
                  }

                  // Extract visual content if present (display-only, not sent to LLM)
                  if( result.contains( "_visuals" ) )
                  {
                    const json &visuals = result["_visuals"];
                    if( visuals.is_array() )
                    {
                      for( const json &vis : visuals )
                      {
                        LlmToolCall::VisualContent vc;
                        const std::string typeStr = vis.value( "type", "" );
                        if( typeStr == "chi_squared_plot" )
                          vc.type = LlmToolCall::VisualContent::VisualType::ChiSquaredPlot;
                        else if( typeStr == "shielding_diagram_2d" )
                          vc.type = LlmToolCall::VisualContent::VisualType::ShieldingDiagram2D;
                        else if( typeStr == "spectrum_chart" )
                          vc.type = LlmToolCall::VisualContent::VisualType::SpectrumChart;
                        else if( typeStr == "svg" )
                          vc.type = LlmToolCall::VisualContent::VisualType::Svg;
                        else if( typeStr == "html" )
                          vc.type = LlmToolCall::VisualContent::VisualType::Html;
                        else
                          continue;

                        // For string types (HTML, SVG), get as string; for JSON objects, dump to string
                        if( vis.contains( "data" ) )
                        {
                          if( vis["data"].is_string() )
                            vc.data = vis["data"].get<std::string>();
                          else
                            vc.data = vis["data"].dump();
                        }
                        vc.title = vis.value( "title", "" );
                        vc.description = vis.value( "description", "" );
                        vc.suggestedWidthPx = vis.value( "widthPx", 0 );
                        vc.suggestedHeightPx = vis.value( "heightPx", 0 );

                        // For SpectrumChart: extract measurement pointer, peak JSON, and energy range
                        if( vc.type == LlmToolCall::VisualContent::VisualType::SpectrumChart )
                        {
                          vc.peakJson = vis.value( "peakJson", "" );
                          vc.xRangeMin = vis.value( "xRangeMin", 0.0 );
                          vc.xRangeMax = vis.value( "xRangeMax", 0.0 );
                          InterSpec *viewer = InterSpec::instance();
                          if( viewer )
                          {
                            std::shared_ptr<const SpecUtils::Measurement> dispMeas
                              = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
                            if( dispMeas )
                              vc.spectrumMeasurement = std::make_shared<SpecUtils::Measurement>( *dispMeas );
                          }
                        }

                        pendingResult->visualContent.push_back( std::move( vc ) );
                      }//for( const json &vis : visuals )
                    }//if( visuals.is_array() )

                    result.erase( "_visuals" );
                  }//if( result.contains( "_visuals" ) )

                  pendingResult->content = result.dump();
                }

                pendingResult->executionDuration = exec_duration;

                // Notify the display widget that this tool result has been updated
                if( owningToolResults )
                  owningToolResults->asyncToolCompleted().emit( capturedCallId );

                // Decrement pending count; send results when all async tools are done
                auto defPos = self->m_deferredToolResults.find( capturedParentRequestId );
                if( defPos == self->m_deferredToolResults.end() )
                {
                  if( self->m_debug_stream )
                    (*self->m_debug_stream) << "Async callback: no deferred entry for parentRequestId="
                      << capturedParentRequestId << " - discarding" << endl;
                  return;
                }

                assert( defPos->second.pendingAsyncToolCount > 0 );
                defPos->second.pendingAsyncToolCount--;

                if( self->m_debug_stream )
                  (*self->m_debug_stream) << "Async tool complete: " << capturedToolName
                    << ", remaining async=" << defPos->second.pendingAsyncToolCount
                    << ", remaining sub-agents=" << defPos->second.subAgentToolCallIds.size() << endl;

                if( defPos->second.pendingAsyncToolCount == 0
                    && defPos->second.subAgentToolCallIds.empty() )
                {
                  self->m_deferredToolResults.erase( defPos );

                  if( self->m_debug_stream )
                    (*self->m_debug_stream) << "All async/sub-agent tools complete - sending results to LLM" << endl;

                  self->sendToolResultsToLLM( convo );
                }
              }
            );//asyncExecutor callback
            }catch( const std::exception &e )
            {
              // The asyncExecutor threw synchronously (e.g., during argument parsing)
              // before establishing the async callback.  asyncToolCount was already
              // incremented and a Pending result was pushed, but the callback will
              // never fire, so we must handle the error here.
              asyncToolCount--;

              if( m_debug_stream )
                (*m_debug_stream) << "Async tool threw synchronously: " << e.what() << endl;

              // Update the already-pushed Pending result to Error
              LlmToolCall &pendingEntry = toolCallResults.back();
              assert( pendingEntry.invocationId == callId );
              nlohmann::json errJson;
              errJson["error"] = "Tool call failed: " + std::string( e.what() );
              errJson["success"] = false;
              pendingEntry.status = LlmToolCall::CallStatus::Error;
              pendingEntry.content = errJson.dump();
            }//try / catch for asyncExecutor
          }else
          {
            // === SYNCHRONOUS TOOL PATH ===
            const auto exec_start = std::chrono::steady_clock::now();
            json result = m_tool_registry->executeTool(toolName, arguments, m_interspec, convo, m_history.get());
            const auto exec_end = std::chrono::steady_clock::now();
            const auto exec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(exec_end - exec_start);

            switch( result.type() )
            {
              case nlohmann::detail::value_t::object:
                // We're good, nothing to do here
                break;


              case nlohmann::detail::value_t::null:
              case nlohmann::detail::value_t::discarded:
                result = nlohmann::json::object();
                assert( 0 ); //for development to help ID tool calls where this is the case
                break;

              case nlohmann::detail::value_t::boolean:
              case nlohmann::detail::value_t::array:
              case nlohmann::detail::value_t::string:
              case nlohmann::detail::value_t::number_integer:
              case nlohmann::detail::value_t::number_unsigned:
              case nlohmann::detail::value_t::number_float:
              case nlohmann::detail::value_t::binary:
              {
                auto result_copy = result;
                result = nlohmann::json::object();
                result["value"] = std::move(result_copy);
                assert( 0 ); //for development to help ID tool calls where this is the case
                break;
              }
            }//switch( result.type() )

            assert( result.type() == nlohmann::detail::value_t::object );

            // Add in an indicator, if it isnt already present, that the tool call was succesful - some LLMs or providers
            //  relly on at least `result["status"] = true;`, or they will keep doing the tool calls multiple times
            if( !result.contains("success") && !result.contains("Success") )
              result["success"] = true;

            if( result.contains("status") || result.contains("Status") )
            {
              if( m_debug_stream )
                (*m_debug_stream) << "Tool call result contains 'status'" << endl;
            }

            toolResult.status = LlmToolCall::CallStatus::Success;

            // Extract image content if present (convention: tools put image data in "_image" key)
            if( result.contains( "_image" ) )
            {
              const json &img = result["_image"];
              LlmToolCall::ImageContent ic;
              ic.base64Data = img.value( "base64Data", "" );
              ic.mimeType = img.value( "mimeType", "" );
              ic.widthPx = img.value( "widthPx", 0 );
              ic.heightPx = img.value( "heightPx", 0 );
              if( !ic.base64Data.empty() )
                toolResult.imageContent = std::move( ic );
              result.erase( "_image" );
            }

            // Extract visual content if present (display-only, not sent to LLM)
            if( result.contains( "_visuals" ) )
            {
              const json &visuals = result["_visuals"];
              if( visuals.is_array() )
              {
                for( const json &vis : visuals )
                {
                  LlmToolCall::VisualContent vc;
                  const std::string typeStr = vis.value( "type", "" );
                  if( typeStr == "chi_squared_plot" )
                    vc.type = LlmToolCall::VisualContent::VisualType::ChiSquaredPlot;
                  else if( typeStr == "shielding_diagram_2d" )
                    vc.type = LlmToolCall::VisualContent::VisualType::ShieldingDiagram2D;
                  else if( typeStr == "spectrum_chart" )
                    vc.type = LlmToolCall::VisualContent::VisualType::SpectrumChart;
                  else if( typeStr == "svg" )
                    vc.type = LlmToolCall::VisualContent::VisualType::Svg;
                  else if( typeStr == "html" )
                    vc.type = LlmToolCall::VisualContent::VisualType::Html;
                  else
                    continue;

                  // For string types (HTML, SVG), get as string; for JSON objects, dump to string
                        if( vis.contains( "data" ) )
                        {
                          if( vis["data"].is_string() )
                            vc.data = vis["data"].get<std::string>();
                          else
                            vc.data = vis["data"].dump();
                        }
                  vc.title = vis.value( "title", "" );
                  vc.description = vis.value( "description", "" );
                  vc.suggestedWidthPx = vis.value( "widthPx", 0 );
                  vc.suggestedHeightPx = vis.value( "heightPx", 0 );

                  // For SpectrumChart: extract measurement pointer, peak JSON, and energy range
                  if( vc.type == LlmToolCall::VisualContent::VisualType::SpectrumChart )
                  {
                    vc.peakJson = vis.value( "peakJson", "" );
                    vc.xRangeMin = vis.value( "xRangeMin", 0.0 );
                    vc.xRangeMax = vis.value( "xRangeMax", 0.0 );
                    InterSpec *viewer = InterSpec::instance();
                    if( viewer )
                    {
                      std::shared_ptr<const SpecUtils::Measurement> dispMeas
                        = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
                      if( dispMeas )
                        vc.spectrumMeasurement = std::make_shared<SpecUtils::Measurement>( *dispMeas );
                    }
                  }

                  toolResult.visualContent.push_back( std::move( vc ) );
                }//for( const json &vis : visuals )
              }//if( visuals.is_array() )

              result.erase( "_visuals" );
            }//if( result.contains( "_visuals" ) )

            toolResult.content = result.dump();
            toolResult.executionDuration = exec_duration;
            toolCallResults.push_back( std::move(toolResult) );
          }//if( async tool ) / else ( synchronous tool )
        }//if( m_block_tool_calls ) / else
      }
    }catch( const json::parse_error &e )
    {
      if( m_debug_stream )
        (*m_debug_stream) << "Tool execution JSON parse error: " << e.what() << endl;
      cerr << "Tool execution JSON parse error: " << e.what() << endl;
      

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      {
        const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
        const string now_str = SpecUtils::to_iso_string( now );
        const string debug_name = "llm_toolcall_json_error_" + now_str + ".json";
#ifdef _WIN32
        const std::wstring wdebug_name = SpecUtils::convert_from_utf8_to_utf16(debug_name);
        std::ofstream output_request_json( wdebug_name.c_str(), ios::binary | ios::out );
#else
        std::ofstream output_request_json( debug_name.c_str(), ios::binary | ios::out);
#endif
        output_request_json << rawResponseContent;

        cerr << "Wrote message content to '" << debug_name << "'" << endl;
      }
#endif

      // Create placeholder arguments with error info for tracking
      json arguments = json::object();
      arguments["_raw"] = rawArguments;
      arguments["_parse_error"] = e.what();

      // Note: toolCallRequests was already populated before the parse attempt for normal tools,
      // or before the sub-agent branch for invoke_* tools - so no need to add it here.

      json result;
      result["error"] = "Tool call failed due to JSON parsing: " + string(e.what());

      result["success"] = false;
      
      // Create error tool result entry
      LlmToolCall toolResult( toolName, callId, arguments );
      toolResult.status = LlmToolCall::CallStatus::Error;
      toolResult.content = result.dump();
      toolCallResults.push_back( std::move(toolResult) );
    }catch( const json::exception &e )
    {
      //Will catch: json::exception, json::invalid_iterator, json::type_error, json::out_of_range, and json::other_error here
      cerr << "Tool execution JSON error: " << e.what() << endl;
      if( m_debug_stream )
        (*m_debug_stream) << "Tool execution JSON error: " << e.what() << endl;
      
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      {
        const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
        const string now_str = SpecUtils::to_iso_string( now );
        const string debug_name = "llm_toolcall_json_error_" + now_str + ".json";
#ifdef _WIN32
        const std::wstring wdebug_name = SpecUtils::convert_from_utf8_to_utf16(debug_name);
        std::ofstream output_request_json( wdebug_name.c_str(), ios::binary | ios::out );
#else
        std::ofstream output_request_json( debug_name.c_str(), ios::binary | ios::out);
#endif
        output_request_json << rawResponseContent;
        
        cerr << "Wrote message content to '" << debug_name << "'" << endl;
      }
#endif
      
      // Create placeholder arguments with error info for tracking
      json arguments = json::object();
      arguments["_raw"] = rawArguments;
      arguments["_json_error"] = e.what();
      
      // Note: Do NOT add to toolCallRequests here - it was already added before execution
      // (line 1104 for sub-agents, line 1384 for normal tools)
      
      json result;
      result["error"] = "Tool call failed due to JSON parsing: " + string(e.what());
      
      result["success"] = false;
      
      // Create error tool result entry
      LlmToolCall toolResult( toolName, callId, arguments );
      toolResult.status = LlmToolCall::CallStatus::Error;
      toolResult.content = result.dump();
      toolCallResults.push_back( std::move(toolResult) );
    }catch( const std::exception &e )
    {
      const string msg = e.what();
      if( m_debug_stream )
        (*m_debug_stream) << "Tool execution error: " << msg << endl;

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      {
        const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
        const string now_str = SpecUtils::to_iso_string( now );
        const string debug_name = "llm_toolcall_error_" + now_str + ".json";
#ifdef _WIN32
        const std::wstring wdebug_name = SpecUtils::convert_from_utf8_to_utf16(debug_name);
        std::ofstream output_request_json( wdebug_name.c_str(), ios::binary | ios::out );
#else
        std::ofstream output_request_json( debug_name.c_str(), ios::binary | ios::out);
#endif
        output_request_json << rawResponseContent;

        cerr << "Wrote message content to '" << debug_name << "'" << endl;
      }
#endif

      // Create placeholder arguments with error info for tracking
      json arguments = json::object();
      arguments["_raw"] = rawArguments;
      arguments["_execution_error"] = e.what();

      // Note: Do NOT add to toolCallRequests here - it was already added before execution
      // (line 1104 for sub-agents, line 1384 for normal tools)

      json result;
      result["error"] = "Tool call failed: " + string(e.what());

      result["success"] = false;
      
      // Create error tool result entry
      LlmToolCall toolResult( toolName, callId, arguments );
      toolResult.status = LlmToolCall::CallStatus::Error;
      toolResult.content = result.dump();
      toolCallResults.push_back( std::move(toolResult) );
    }//try / catch

    // Track this tool call for follow-up
    executedToolCallIds.push_back(callId);
  }

  // Add all tool calls to history as a single batched entry
  shared_ptr<LlmToolRequest> tool_request;
  if( !toolCallRequests.empty() )
  {
    tool_request = m_history->addToolCalls( std::move(toolCallRequests), rawResponseContent, convo );

    // Set reasoning details on the tool request BEFORE sendToolResultsToLLM is called
    // These must be set now so they can be included when building the messages array
    if( tool_request )
    {
      tool_request->setThinkingContent( thinkingContent );
      tool_request->setThinkingSignature( thinkingSignature );
      tool_request->setReasoningContent( reasoningContent );
      tool_request->setReasoningDetails( reasoningDetails );
    }
  }

  // Add all tool results to history as a single batched entry (if we have any)
  // For sub-agents, results will be added now as placeholders and updated later
  std::shared_ptr<LlmToolResults> tool_result;
  if( !toolCallResults.empty() )
  {
    // We'll populate jsonSentToLlm when we actually send the results
    // Note: sub_agent_conversation is already stored in each ToolCallRequest for sub-agent invocations
    tool_result = m_history->addToolResults( std::move(toolCallResults), "", convo );
  }

  // If we have sub-agent invocations or async tool calls, defer sending results
  if( !subAgentRequestIds.empty() || (asyncToolCount > 0) )
  {
    if( m_debug_stream )
    {
      (*m_debug_stream) << "Deferring tool results - "
        << subAgentRequestIds.size() << " sub-agent(s), "
        << asyncToolCount << " async tool(s) pending" << endl;
    }

    DeferredToolResult deferred;
    deferred.conversationId = convo->conversationId;
    deferred.toolCallIds = executedToolCallIds;
    deferred.subAgentToolCallIds = subAgentRequestIds;
    deferred.pendingAsyncToolCount = asyncToolCount;

    m_deferredToolResults[parentRequestId] = deferred;

    // Don't send results yet - wait for sub-agent and/or async tools
  }else
  {
    // No sub-agent and no async tools - send results immediately
    if (!executedToolCallIds.empty()) {
      if( m_debug_stream )
        (*m_debug_stream) << "Sending tool results back to LLM for " << executedToolCallIds.size() << " executed tools" << endl;
      sendToolResultsToLLM( convo );
    }
  }

  // Clear current conversation after tool execution
  m_currentConversation.reset();

  return { tool_request, tool_result };
}

std::pair<std::shared_ptr<LlmToolRequest>, std::shared_ptr<LlmToolResults>>
LlmInterface::parseContentForToolCallsAndSendResults( const std::string &content,
                                                            const std::shared_ptr<LlmInteraction> &convo,
                                                            const int requestId,
                                                            const std::string &rawResponseContent,
                                                            const std::string &thinkingContent,
                                                            const std::string &thinkingSignature,
                                                            const std::string &reasoningContent,
                                                            const std::string &reasoningDetails )
{
  if( m_debug_stream )
    (*m_debug_stream) << "Parsing content for text-based tool calls..." << endl;
  
  nlohmann::json tool_calls = nlohmann::json::array();
  // We will add entries to this array that look like:
  // { "id": "call_abc123", "type": "function", "function": { "name": "get_gamma_spectrum", "arguments": "{\"energy_range\": [0, 3000], \"detector\": \"HPGe\"}" }

  // Pre-regex: try to extract tool calls from markdown JSON code blocks.
  //  Some providers (e.g., GPT-OSS) return tool calls inside ```json ... ``` blocks
  //  with a "tool" or "name" key.
  {
    const std::string stripped = extractFromMarkdown( content );
    if( stripped != content )
    {
      try
      {
        json obj = lenientlyParseJson( stripped );

        // Handle a single JSON object with "tool" or "name" key
        auto extractToolCallFromJson = [&]( const json &item ) -> bool
        {
          if( !item.is_object() )
            return false;

          std::string toolName;
          if( item.contains("name") && item["name"].is_string() )
            toolName = item["name"].get<std::string>();
          else if( item.contains("tool") && item["tool"].is_string() )
            toolName = item["tool"].get<std::string>();
          else
            return false;

          json arguments = json::object();
          if( item.contains("arguments") )
            arguments = item["arguments"];
          else if( item.contains("params") )
            arguments = item["params"];
          else if( item.contains("parameters") )
            arguments = item["parameters"];

          toolName = sanitizeToolName( toolName );
          if( m_tool_registry )
            toolName = resolveToolNameCaseInsensitive( toolName, *m_tool_registry );

          const auto epoch_ticks = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch() ).count();
          const std::string invocationId = "md_call_" + std::to_string( epoch_ticks );

          json function_def = json::object();
          function_def["name"] = toolName;
          function_def["arguments"] = arguments.dump();

          json call_def = json::object();
          call_def["id"] = invocationId;
          call_def["type"] = "function";
          call_def["function"] = std::move( function_def );

          tool_calls.push_back( std::move( call_def ) );

          if( m_debug_stream )
            (*m_debug_stream) << "Extracted tool call from markdown code block: " << toolName << endl;

          return true;
        };//extractToolCallFromJson lambda

        if( obj.is_object() )
        {
          extractToolCallFromJson( obj );
        }
        else if( obj.is_array() )
        {
          for( const auto &item : obj )
            extractToolCallFromJson( item );
        }
      }catch( const std::exception & )
      {
        // Not valid JSON in code block - fall through to regex patterns
      }
    }//if( stripped != content ) -- i.e., content had a markdown code fence
  }

  // If we found tool calls from markdown code blocks, skip regex patterns
  if( !tool_calls.empty() )
  {
    if( m_debug_stream )
      (*m_debug_stream) << "=== Found " << tool_calls.size() << " tool calls from markdown code blocks ===" << endl;

    return executeToolCallsAndSendResults( tool_calls, convo, requestId, rawResponseContent,
      thinkingContent, thinkingSignature, reasoningContent, reasoningDetails );
  }

  // Look for various tool call patterns
  std::vector<std::regex> toolCallPatterns = {
    // Pattern: [TOOL_REQUEST] {"name": "tool_name", "arguments": {...}} [END_TOOL_REQUEST]
    std::regex("\\[TOOL_REQUEST\\]\\s*([^\\[]+)\\s*\\[END_TOOL_REQUEST\\]"),
    // Pattern: [TOOL_REQUEST] tool_name {"param": "value"} [/TOOL_REQUEST]
    std::regex("\\[TOOL_REQUEST\\]\\s*(\\w+)\\s*(\\{[^}]*\\})\\s*\\[/TOOL_REQUEST\\]"),
    // Pattern: <tool_call name="tool_name">{"param": "value"}</tool_call>
    std::regex("<tool_call\\s+name=\"([^\"]+)\">([^<]*)</tool_call>"),
    // Pattern: <tool_call>{"name": "tool_name", "arguments": {...}}</tool_call>
    std::regex("<tool_call>\\s*([^<]+)\\s*</tool_call>"),
    // Pattern: Tool: tool_name Arguments: {"param": "value"}
    std::regex("Tool:\\s*(\\w+)\\s*Arguments:\\s*(\\{[^}]*\\})"),
    // Pattern: detected_peaks({"specType": "Foreground"})
    std::regex("(\\w+)\\s*\\(\\s*(\\{[^}]*\\})\\s*\\)")
  };
  
  std::smatch match;
  for( size_t patIdx = 0; patIdx < toolCallPatterns.size(); ++patIdx )
  {
    const std::regex &pattern = toolCallPatterns[patIdx];
    std::sregex_iterator iter(content.begin(), content.end(), pattern);
    std::sregex_iterator end;
    
    for( ; iter != end; ++iter )
    {
      const std::smatch& match = *iter;
      
      if( match.size() < 2 )
        continue;
      
      // Generate a unique invocation ID for this tool call (OpenAI-compatible format)
      const string invocationId = generateToolCallId();
      
      string toolName;
      json arguments = json::object();
      
      try
      {
        if( match.size() == 2 )
        {
          // Single capture group - JSON object with name and arguments fields
          string toolCallJson = match[1].str();
          if( m_debug_stream )
            (*m_debug_stream) << "Found JSON tool call: " << toolCallJson << " -- and about to parse" << endl;

          json toolCallObj = lenientlyParseJson( toolCallJson );
          if( m_debug_stream )
            (*m_debug_stream) << "Done parsing tool call JSON" << endl;

          if( toolCallObj.contains("name") )
            toolName = toolCallObj["name"];
          else if( toolCallObj.contains("tool") )
            toolName = toolCallObj["tool"];
          else
            throw std::runtime_error( "Invalid tool call JSON format - missing name" );

          if( toolCallObj.contains("arguments") )
            arguments = toolCallObj["arguments"];
          else if( toolCallObj.contains("params") )
            arguments = toolCallObj["params"];
          else if( toolCallObj.contains("parameters") )
            arguments = toolCallObj["parameters"];
        }else
        {
          assert( match.size() >= 3 );

          // Two capture groups - tool name and arguments separately
          toolName = match[1].str();
          string argumentsStr = match[2].str();

          // Parse arguments JSON
          if( !argumentsStr.empty() && (argumentsStr != "{}") )
          {
            if( m_debug_stream )
              (*m_debug_stream) << "About to parse alone tool call JSON" << endl;

            arguments = lenientlyParseJson( argumentsStr );

            if( m_debug_stream )
              (*m_debug_stream) << "Done parsing alone tool call JSON" << endl;
          }
        }

        if( m_debug_stream )
          (*m_debug_stream) << "Found text-based tool call: " << toolName << " with args: " << arguments.dump() << endl;
      }catch( const std::exception &e )
      {
        if( m_debug_stream )
          (*m_debug_stream) << "Error parsing text-based tool call: " << e.what() << endl;

        // Try to extract tool name from the matched content as a fallback
        if( toolName.empty() )
        {
          // Try to extract name from JSON-like content in the match
          // Pattern: look for "name": "tool_name" or \"name\": \"tool_name\"
          const string matchStr = match[0].str();
          std::regex namePattern( R"(["']?name["']?\s*:\s*["']([^"']+)["'])" );
          std::smatch nameMatch;

          if( std::regex_search( matchStr, nameMatch, namePattern ) && nameMatch.size() > 1 )
          {
            toolName = nameMatch[1].str();
            if( m_debug_stream )
              (*m_debug_stream) << "Extracted tool name from fallback pattern: " << toolName << endl;
          }else
          {
            // Last resort - use a descriptive error placeholder
            toolName = "unknown_tool_parse_error";
            if( m_debug_stream )
              (*m_debug_stream) << "Failed to extract tool name from match, using placeholder: " << toolName << endl;
          }
        }
      }//try / catch

      // Sanitize and resolve tool name case-insensitively
      toolName = sanitizeToolName( toolName );
      if( m_tool_registry )
        toolName = resolveToolNameCaseInsensitive( toolName, *m_tool_registry );

      // For aggressive patterns (function-call syntax), validate that the tool name
      //  is actually a known tool to avoid false positives
      if( (patIdx >= 5) && m_tool_registry && !m_tool_registry->getTool( toolName ) )
      {
        if( m_debug_stream )
          (*m_debug_stream) << "Skipping false-positive match for unknown tool: " << toolName << endl;
        continue;
      }

      // Add an entry to `tool_calls` that look like:
      // { "id": "call_abc123", "type": "function", "function": { "name": "get_gamma_spectrum", "arguments": "{...}" }
      // (note, adding function call, even if we got an exception extracting the tool call - we will let
      //  `executeToolCallsAndSendResults(...)` deal with errors).
      //
      // IMPORTANT: arguments must be stored as a STRING (JSON serialized), not as a JSON object,
      // because executeToolCallsAndSendResults expects to parse it from a string
      nlohmann::json function_def = nlohmann::json::object();
      function_def["name"] = toolName;
      function_def["arguments"] = arguments.dump();  // Serialize to string, not move as object

      nlohmann::json call_def = nlohmann::json::object();
      call_def["id"] = invocationId;
      call_def["type"] = "function";
      call_def["function"] = std::move(function_def);

      tool_calls.push_back( std::move(call_def) );
    }//for( ; iter != end; ++iter )
  }//for( size_t patIdx = 0; patIdx < toolCallPatterns.size(); ++patIdx )

  // Post-regex: if no tool calls found yet, check if the entire content is a bare
  //  tool call like "tool_name{json}" (e.g., Devstral returning
  //  `add_analysis_peaks_for_source{"source": "Eu152"}`).
  //  Only match when this IS the entire content to avoid false positives.
  if( tool_calls.empty() && m_tool_registry )
  {
    std::regex entireContentPattern( R"(^\s*(\w+)\s*(\{[\s\S]*\})\s*$)" );
    std::smatch entireMatch;
    if( std::regex_match( content, entireMatch, entireContentPattern ) && (entireMatch.size() >= 3) )
    {
      std::string toolName = sanitizeToolName( entireMatch[1].str() );
      toolName = resolveToolNameCaseInsensitive( toolName, *m_tool_registry );

      if( m_tool_registry->getTool( toolName ) )
      {
        try
        {
          const json arguments = lenientlyParseJson( entireMatch[2].str() );
          const std::string invocationId = generateToolCallId();

          json function_def = json::object();
          function_def["name"] = toolName;
          function_def["arguments"] = arguments.dump();

          json call_def = json::object();
          call_def["id"] = invocationId;
          call_def["type"] = "function";
          call_def["function"] = std::move( function_def );

          tool_calls.push_back( std::move( call_def ) );

          if( m_debug_stream )
            (*m_debug_stream) << "Extracted bare tool call from entire content: " << toolName << endl;
        }catch( const std::exception &e )
        {
          if( m_debug_stream )
            (*m_debug_stream) << "Failed to parse bare tool call JSON: " << e.what() << endl;
        }
      }//if( tool name is a known tool )
    }//if( entire content matches tool_name{json} pattern )
  }//if( no tool calls found yet and we have a tool registry )

  if( m_debug_stream )
    (*m_debug_stream) << "=== Complete extracting " << tool_calls.size() << " text-based tool calls ===" << endl;

  pair<shared_ptr<LlmToolRequest>,shared_ptr<LlmToolResults>> result
        = executeToolCallsAndSendResults( tool_calls, convo, requestId, rawResponseContent,
            thinkingContent, thinkingSignature, reasoningContent, reasoningDetails );

  return result;
}//parseContentForToolCallsAndSendResults



/** Sanitize a string to ensure it is valid UTF-8.

 Invalid UTF-8 sequences are replaced with the Unicode replacement character (U+FFFD).
 This is necessary because nlohmann::json requires valid UTF-8 strings.

 @param str The string to sanitize
 @return A valid UTF-8 string
 */
static std::string sanitizeUtf8( const std::string &str )
{
  std::string result;
  result.reserve( str.size() );

  size_t invalidCount = 0;

  for( size_t i = 0; i < str.size(); )
  {
    const unsigned char byte = static_cast<unsigned char>(str[i]);

    // Single-byte character (0xxxxxxx)
    if( byte < 0x80 )
    {
      result += str[i];
      ++i;
      continue;
    }

    // Determine the number of bytes in this UTF-8 sequence
    size_t seqLen = 0;
    unsigned char mask = 0;

    if( (byte & 0xE0) == 0xC0 )      // 2-byte sequence (110xxxxx)
    {
      seqLen = 2;
      mask = 0x1F;
    }
    else if( (byte & 0xF0) == 0xE0 ) // 3-byte sequence (1110xxxx)
    {
      seqLen = 3;
      mask = 0x0F;
    }
    else if( (byte & 0xF8) == 0xF0 ) // 4-byte sequence (11110xxx)
    {
      seqLen = 4;
      mask = 0x07;
    }
    else
    {
      // Invalid start byte - replace with replacement character
      result += "\xEF\xBF\xBD"; // U+FFFD in UTF-8
      ++invalidCount;
      ++i;
      continue;
    }

    // Check if we have enough bytes
    if( i + seqLen > str.size() )
    {
      // Truncated sequence - replace with replacement character
      result += "\xEF\xBF\xBD"; // U+FFFD in UTF-8
      ++invalidCount;
      break;
    }

    // Validate continuation bytes (10xxxxxx)
    bool valid = true;
    for( size_t j = 1; j < seqLen; ++j )
    {
      const unsigned char contByte = static_cast<unsigned char>(str[i + j]);
      if( (contByte & 0xC0) != 0x80 )
      {
        valid = false;
        break;
      }
    }

    if( valid )
    {
      // Copy the valid UTF-8 sequence
      for( size_t j = 0; j < seqLen; ++j )
        result += str[i + j];
      i += seqLen;
    }
    else
    {
      // Invalid sequence - replace with replacement character
      result += "\xEF\xBF\xBD"; // U+FFFD in UTF-8
      ++invalidCount;
      ++i;
    }
  }

#if( PERFORM_DEVELOPER_CHECKS )
  if( invalidCount > 0 )
  {
    std::cerr << "sanitizeUtf8: Replaced " << invalidCount << " invalid UTF-8 sequence(s)" << std::endl;
    std::cerr << "  Original length: " << str.size() << " bytes" << std::endl;
    std::cerr << "  Sanitized length: " << result.size() << " bytes" << std::endl;

    // Show a snippet of the original string with hex values for debugging
    std::cerr << "  First 200 bytes (hex): ";
    for( size_t i = 0; i < std::min(size_t(200), str.size()); ++i )
    {
      std::cerr << std::hex << std::setw(2) << std::setfill('0')
                << static_cast<int>(static_cast<unsigned char>(str[i])) << " ";
    }
    std::cerr << std::dec << std::endl;
  }
#endif

  return result;
}

std::string LlmInterface::stripThinkingContent(const std::string& content) {
  std::string result = content;
  
  // Remove <think>...</think> blocks (case insensitive, supports nested and multiline)
  // Use [\s\S] instead of . to match any character including newlines
  std::regex thinkRegex("<think[^>]*>[\\s\\S]*?</think>", 
    std::regex_constants::icase | std::regex_constants::ECMAScript);
  
  // Keep removing think blocks until no more are found (handles nested cases)
  std::string prevResult;
  do {
    prevResult = result;
    result = std::regex_replace(result, thinkRegex, "");
  } while (result != prevResult);
  
  // Clean up extra whitespace that may be left after removing think blocks
  // Replace multiple newlines with at most two newlines
  std::regex multiNewlineRegex("\n\n\n+");
  result = std::regex_replace(result, multiNewlineRegex, "\n\n");
  
  // Trim leading and trailing whitespace
  result = std::regex_replace(result, std::regex("^\\s+"), "");
  result = std::regex_replace(result, std::regex("\\s+$"), "");
  
  return result;
}

std::pair<std::string, std::string> LlmInterface::extractThinkingAndContent(const std::string& content) {
  std::string cleanContent = content;
  std::string thinkingContent;
  
  // Extract <think>...</think> blocks (case insensitive, supports nested and multiline)
  std::regex thinkRegex("<think[^>]*>([\\s\\S]*?)</think>", 
    std::regex_constants::icase | std::regex_constants::ECMAScript);
  
  std::smatch match;
  std::string::const_iterator searchStart(content.cbegin());
  
  // Collect all thinking content from properly formatted <think>...</think> blocks
  while (std::regex_search(searchStart, content.cend(), match, thinkRegex)) {
    if (!thinkingContent.empty()) {
      thinkingContent += "\n";
    }
    thinkingContent += match[1].str();
    searchStart = match.suffix().first;
  }
  
  // Check for orphaned </think> tags (closing tag without opening tag)
  std::regex orphanedCloseRegex("</think>", std::regex_constants::icase);
  std::smatch orphanedMatch;
  if (std::regex_search(content, orphanedMatch, orphanedCloseRegex)) {
    // Find the position of the first </think> tag
    size_t closePos = orphanedMatch.position();
    
    // Extract all text before the </think> tag as thinking content
    std::string orphanedThinking = content.substr(0, closePos);
    
    // Only add if we found some content and it's not just whitespace
    if (!orphanedThinking.empty() && orphanedThinking.find_first_not_of("\\s\\t\\n\\r") != std::string::npos) {
      if (!thinkingContent.empty()) {
        thinkingContent += "\n";
      }
      thinkingContent += orphanedThinking;
      
      // Remove the orphaned thinking content and </think> tag from clean content
      cleanContent = content.substr(closePos + 7); // 7 is length of "</think>"
    }
  } else {
    // No orphaned tags, use normal stripping
    cleanContent = stripThinkingContent(content);
  }
  
  return {cleanContent, thinkingContent};
}

std::string LlmInterface::getSystemPromptForAgent( const AgentType agentType ) const
{
  if( !m_config || !m_config->llmApi.enabled )
    return "";

  // Search for agent in config
  for( const LlmConfig::AgentConfig &agent : m_config->agents )
  {
    if( agent.type == agentType )
      return agent.systemPrompt;
  }

  assert( 0 );
  throw std::runtime_error( "Failed to find agent config for agent " + agentTypeToString(agentType) );
  
  return "";
}//getSystemPromptForAgent(...)


/*
void LlmInterface::initializeStateMachineForConversation( std::shared_ptr<LlmInteraction> convo ) const
{
  if( !convo || !m_config )
    return;

  // Find the agent config for this conversation
  for( const LlmConfig::AgentConfig &agent : m_config->agents )
  {
    if( agent.type == convo->agent_type )
    {
      // If this agent has a state machine defined, create a fresh copy for this conversation
      if( agent.state_machine )
      {
        // Create an independent copy of the state machine for this conversation
        convo->state_machine = agent.state_machine->copy();

        if( m_debug_stream )
          (*m_debug_stream) << "Initialized state machine for " << agentTypeToString(agent.type)
             << " conversation, starting in state: " << convo->state_machine->getCurrentState() << endl;
      }
      return;
    }
  }

  // If we get here, we didn't find the agent config (shouldn't happen)
  cerr << "Warning: Could not find agent config for " << agentTypeToString(convo->agent_type) << endl;
}//initializeStateMachineForConversation(...)
*/


std::shared_ptr<LlmInteraction> LlmInterface::getCurrentConversation() const
{
  return m_currentConversation.lock();
}//getCurrentConversation()


int LlmInterface::invokeSubAgent( std::shared_ptr<LlmInteraction> sub_agent_convo )
{
  assert( sub_agent_convo );
  if( !sub_agent_convo )
    throw std::logic_error( "LlmInterface::invokeSubAgent called with null conversation" );

  assert( sub_agent_convo->agent_type != AgentType::MainAgent );

  // Initialize state machine if the agent has one
  //initializeStateMachineForConversation( sub_agent_convo );

  // The initial request should be in the responses array as an InitialRequest turn
  assert( !sub_agent_convo->responses.empty() );
  assert( sub_agent_convo->responses.front()->type() == LlmInteractionTurn::Type::InitialRequest );

  // Get content from InitialRequest for debug output
  string contextTask;
  const LlmInteractionInitialRequest *initialReq = dynamic_cast<const LlmInteractionInitialRequest*>(sub_agent_convo->responses.front().get());
  if( initialReq )
    contextTask = initialReq->content();

  if( m_debug_stream )
  {
    (*m_debug_stream) << "=== Invoking sub-agent: " << agentTypeToString(sub_agent_convo->agent_type) << " (async with pause) ===" << endl;
    (*m_debug_stream) << "Context/Task: " << contextTask.substr(0, 100) << (contextTask.length() > 100 ? "..." : "") << endl;
    (*m_debug_stream) << "Invoke tool call ID: " << sub_agent_convo->conversationId << endl;
  }

  // Build API request for this specific agent - since the combined message is already in the
  const json requestJson = buildMessagesArray( sub_agent_convo );
  
  // Make tracked API call
  const int requestId = m_nextRequestId++;

  // Create pending request with sub-agent info
  PendingRequest pending;
  pending.requestId = requestId;
  pending.conversation = sub_agent_convo;
  pending.isSubAgentRequest = true;

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  pending.requestJson = requestJson;
#endif

  m_pendingRequests[requestId] = pending;

  if( m_debug_stream )
  {
    (*m_debug_stream) << "Main agent will pause until sub-agent completes" << endl;
  }

  // Make the actual API call
  std::string msg_content = makeApiCallWithId(requestJson, requestId);

  // Store raw JSON request in the InitialRequest turn's rawContent field
  if( !sub_agent_convo->responses.empty() )
  {
    std::shared_ptr<LlmInteractionTurn> firstTurn = sub_agent_convo->responses.front();
    if( firstTurn && firstTurn->type() == LlmInteractionTurn::Type::InitialRequest )
      firstTurn->setRawContent( std::move( msg_content ) );
  }

  return requestId;
}//invokeSubAgent(...)


nlohmann::json LlmInterface::buildMessagesArray( const std::shared_ptr<LlmInteraction> &convo )
{
  if( !m_config || !m_config->llmApi.enabled )
    throw std::logic_error( "LlmInterface: not configured" );
  
  assert( convo );
  if( !convo )
    throw std::logic_error( "LlmInterface::buildMessagesArray: null conversation history passed in." );
  
  json request;
  request["model"] = m_config->llmApi.model();
  
  // Use max_completion_tokens for newer OpenAI models, max_tokens for others
  if( m_config->llmApi.maxTokens() > 0 )
  {
    string modelName = m_config->llmApi.model();
    if( (modelName.find("gpt-") != string::npos)
       || (modelName.find("o1") != string::npos)
       || (modelName.find("o3") != string::npos)
       || (modelName.find("o4") != string::npos)
       || (modelName.find("gemini") != string::npos) )
    {
      request["max_completion_tokens"] = m_config->llmApi.maxTokens();
    } else
    {
      request["max_tokens"] = m_config->llmApi.maxTokens();
    }
  }//if( m_config->llmApi.maxTokens() > 0 )

  // Add reasoning configuration based on variant type
  if( std::holds_alternative<bool>(m_config->llmApi.reasoning()) )
  {
    // Boolean reasoning (OpenRouter style): reasoning: {enabled: true}
    if( std::get<bool>(m_config->llmApi.reasoning()) )
      request["reasoning"]["enabled"] = true;
  }else if( std::holds_alternative<LlmConfig::LlmApi::ReasoningEffort>(m_config->llmApi.reasoning()) )
  {
    // Effort-level reasoning (OpenAI o1/o3 style): reasoning_effort: "low"|"medium"|"high"
    const LlmConfig::LlmApi::ReasoningEffort effort = std::get<LlmConfig::LlmApi::ReasoningEffort>(m_config->llmApi.reasoning());
    switch( effort )
    {
      case LlmConfig::LlmApi::ReasoningEffort::low:
        request["reasoning_effort"] = "low";
        break;
      case LlmConfig::LlmApi::ReasoningEffort::medium:
        request["reasoning_effort"] = "medium";
        break;
      case LlmConfig::LlmApi::ReasoningEffort::high:
        request["reasoning_effort"] = "high";
        break;
    }
  }//if( std::holds_alternative<bool>(m_config->llmApi.reasoning()) ) / else

  // Add optional temperature parameter if configured
  if( m_config->llmApi.temperature().has_value() )
    request["temperature"] = m_config->llmApi.temperature().value();


  // Build system prompt with optional state machine overview
  string systemPrompt = getSystemPromptForAgent( convo->agent_type );

  if( convo->state_machine )
  {
    systemPrompt += "\n\n## Workflow State Machine\n\n";
    systemPrompt += "You will follow a structured workflow through these states:\n\n";
    systemPrompt += convo->state_machine->getFullStateMapSummary();
    systemPrompt += "\n\nState transitions are managed via the `set_workflow_state` tool. ";
    systemPrompt += "After each transition, you'll receive guidance for the new state.";
  }

  // Build messages array - delegate to LlmConversationHistory for MainAgent
  json messages;
  if( convo->agent_type == AgentType::MainAgent )
  {
    if( m_debug_stream )
    {
      const size_t nConvos = m_history->getConversations().size();
      if( nConvos > 0 )
        (*m_debug_stream) << "=== Including " << nConvos << " history messages in request ===" << endl;
      else
        (*m_debug_stream) << "=== No history to include ===" << endl;
    }

    messages = m_history->buildMessagesForMainAgent( convo, systemPrompt );
  }else
  {
    // For sub-agents, just include system prompt + current conversation
    messages = json::array();
    if( !systemPrompt.empty() )
    {
      json systemMsg;
      systemMsg["role"] = "system";
      systemMsg["content"] = systemPrompt;
      messages.push_back( systemMsg );
    }
    LlmConversationHistory::addConversationToLlmApiHistory( *convo, messages );
  }

  // Add ephemeral state machine reminder as the last message (not persisted in history).
  // Only inject after the first call (when there are already responses), to avoid consecutive
  // user messages on the initial request where the system prompt already has the state info.
  // Skip if the most recent interaction was a set_workflow_state call, since that tool's response
  // already contains all the state information.
  if( convo->state_machine && (convo->responses.size() > 1) )
  {
    bool skipEphemeral = false;
    if( !convo->responses.empty() )
    {
      const std::shared_ptr<LlmInteractionTurn> &lastResponse = convo->responses.back();
      if( lastResponse && lastResponse->type() == LlmInteractionTurn::Type::ToolResult )
      {
        const LlmToolResults *toolResults = dynamic_cast<const LlmToolResults *>( lastResponse.get() );
        if( toolResults )
        {
          for( const LlmToolCall &call : toolResults->toolCalls() )
          {
            if( call.toolName == "set_workflow_state" )
            {
              skipEphemeral = true;
              break;
            }
          }
        }
      }
    }

    if( !skipEphemeral )
    {
      const string &curState = convo->state_machine->getCurrentState();
      const vector<string> nextStates = convo->state_machine->getAllowedTransitions();

      if( !curState.empty() )
      {
        string ephemeralContent = "[System: Current workflow state is " + curState;

        if( !nextStates.empty() )
        {
          ephemeralContent += ". Allowed transitions: ";
          for( size_t i = 0; i < nextStates.size(); ++i )
          {
            if( i > 0 )
              ephemeralContent += ", ";
            ephemeralContent += nextStates[i];
          }
          ephemeralContent += ". Call set_workflow_state when ready to transition";
        }
        ephemeralContent += ".]";

        const string stateEphemeralTxt = convo->state_machine->getEphemeralMessageTxtForCurrentState();
        if( !stateEphemeralTxt.empty() )
          ephemeralContent += "\n\n" + stateEphemeralTxt;

        // Some providers (e.g., Mistral) require an assistant turn immediately after tool result
        // messages (role "tool"), so we can't insert a user message there.  Instead, append the
        // ephemeral state info to the last tool result message's content.
        if( !messages.empty() && (messages.back()["role"] == "tool") )
        {
          json &lastMsg = messages.back();
          if( lastMsg["content"].is_string() )
          {
            const string existing = lastMsg["content"].get<string>();
            lastMsg["content"] = existing + "\n\n" + ephemeralContent;
          }
          else
          {
            assert( lastMsg["content"].is_array() );
            if( lastMsg["content"].is_array() )
            {
              json textBlock;
              textBlock["type"] = "text";
              textBlock["text"] = ephemeralContent;
              lastMsg["content"].push_back( textBlock );
            }
          }
        }
        else
        {
          json stateMsg;
          stateMsg["role"] = "user";
          stateMsg["content"] = ephemeralContent;
          messages.push_back( stateMsg );
        }
      }
    }
  }//if( state machine ephemeral reminder )

  request["messages"] = messages;

  // Add tools
  json tools = json::array();
  
  const map<string, LlmTools::SharedTool> agentTools = m_tool_registry->getToolsForAgent(convo->agent_type);
  for( const auto &[name, tool] : agentTools )
  {
    json toolDef;
    toolDef["type"] = "function";
    toolDef["function"]["name"] = tool.name;
    toolDef["function"]["description"] = m_tool_registry->getDescriptionForAgent(tool.name, convo->agent_type);
    
    // Functions used to include "userSession" argument for the MCP server (now thats added in by the MCP server) - but we'll make sure of this here for the moment until we verify they have been totally removed.
    nlohmann::json par_schema = tool.parameters_schema;
    assert( !par_schema.contains("properties") || !par_schema["properties"].contains("userSession") );
    if( par_schema.contains("properties") && par_schema["properties"].contains("userSession") )
      par_schema["properties"].erase( "userSession" );
    
    toolDef["function"]["parameters"] = par_schema;
    tools.push_back(toolDef);
  }//for( const auto &[name, tool] : m_tool_registry->getTools() )
  
  assert( !tools.empty() );
  if( !tools.empty() )
  {
    request["tools"] = tools;

    // Use "required" when configured and in a non-final state machine state,
    // to force the model to make a tool call rather than just producing text.
    bool useRequired = false;
    if( m_config->llmApi.requireToolInStateMachine()
       && convo->state_machine
       && !convo->state_machine->isFinalState( convo->state_machine->getCurrentState() ) )
    {
      useRequired = true;
    }

    request["tool_choice"] = useRequired ? "required" : "auto";
  }

  return request;
}//nlohmann::json LlmInterface::buildMessagesArray( convo )


std::string LlmInterface::buildConfigFingerprint( const nlohmann::json &requestJson ) const
{
  string fp;
  if( requestJson.contains("model") )
    fp += requestJson["model"].dump();
  if( requestJson.contains("tools") )
    fp += "|tools=" + std::to_string( requestJson["tools"].size() );
  if( requestJson.contains("max_tokens") )
    fp += "|mt=" + requestJson["max_tokens"].dump();
  if( requestJson.contains("max_completion_tokens") )
    fp += "|mct=" + requestJson["max_completion_tokens"].dump();
  return fp;
}//string buildConfigFingerprint(...)


std::string LlmInterface::buildSerializedPrefix( const nlohmann::json &requestJson ) const
{
  string result = "{";
  bool first = true;

  auto addField = [&]( const string &key )
  {
    if( !requestJson.contains(key) )
      return;
    if( !first )
      result += ",";
    first = false;
    result += "\"" + key + "\":";
    const nlohmann::json sanitized = sanitizeJsonUtf8( requestJson[key] );
    result += sanitized.dump();
  };

  // Same priority order as serializeRequestForCaching
  addField( "model" );
  addField( "temperature" );
  addField( "max_completion_tokens" );
  addField( "max_tokens" );
  addField( "reasoning" );
  addField( "reasoning_effort" );
  addField( "tools" );
  addField( "tool_choice" );

  // Add any remaining non-messages fields
  for( auto it = requestJson.begin(); it != requestJson.end(); ++it )
  {
    const string &key = it.key();
    if( key == "model" || key == "temperature"
       || key == "max_completion_tokens" || key == "max_tokens"
       || key == "reasoning" || key == "reasoning_effort"
       || key == "tools" || key == "tool_choice"
       || key == "messages" )
      continue;

    if( !first )
      result += ",";
    first = false;
    result += "\"" + key + "\":";
    const nlohmann::json sanitized = sanitizeJsonUtf8( it.value() );
    result += sanitized.dump();
  }

  // Open the messages array
  if( !first )
    result += ",";
  result += "\"messages\":[";

  return result;
}//string buildSerializedPrefix(...)


std::string LlmInterface::serializeRequestIncremental( const nlohmann::json &requestJson )
{
  try
  {
    // Check if prefix cache is valid
    const string fingerprint = buildConfigFingerprint( requestJson );

    if( !m_serializedCache.isValid() || m_serializedCache.configFingerprint != fingerprint )
    {
      m_serializedCache.clear();
      m_serializedCache.configFingerprint = fingerprint;
      m_serializedCache.prefix = buildSerializedPrefix( requestJson );
    }

    // Get messages array
    assert( requestJson.contains("messages") && requestJson["messages"].is_array() );
    const nlohmann::json &messages = requestJson["messages"];
    const size_t totalMessages = messages.size();

    // If cached count > total (history was compacted or reset), rebuild
    if( m_serializedCache.serializedMessageCount > totalMessages )
    {
      m_serializedCache.serializedMessages.clear();
      m_serializedCache.serializedMessageCount = 0;
    }

    // Serialize only new messages
    for( size_t i = m_serializedCache.serializedMessageCount; i < totalMessages; ++i )
    {
      if( !m_serializedCache.serializedMessages.empty() )
        m_serializedCache.serializedMessages += ",";

      const nlohmann::json sanitized = sanitizeJsonUtf8( messages[i] );
      m_serializedCache.serializedMessages += sanitized.dump();
    }
    m_serializedCache.serializedMessageCount = totalMessages;

    // Assemble final string: prefix + messages + close
    return m_serializedCache.prefix + m_serializedCache.serializedMessages + "]}";
  }catch( const std::exception &e )
  {
    cerr << "ERROR in serializeRequestIncremental: " << e.what()
      << " - falling back to full serialization" << endl;
    m_serializedCache.clear();
    return serializeRequestForCaching( requestJson );
  }
}//string serializeRequestIncremental(...)


void LlmInterface::storeRawContentForToolResults( const std::shared_ptr<LlmInteraction> &convo,
                                                   const nlohmann::json &followupRequest )
{
  if( convo->responses.empty() )
    return;

  // Find the most recent ToolResult, then find the ToolRequest immediately preceding it
  for( auto it = convo->responses.rbegin(); it != convo->responses.rend(); ++it )
  {
    if( (*it)->type() == LlmInteractionTurn::Type::ToolResult )
    {
      std::shared_ptr<LlmToolResults> toolResults = std::dynamic_pointer_cast<LlmToolResults>( *it );

      // Find the ToolRequest that immediately precedes this ToolResult
      std::shared_ptr<LlmToolRequest> mostRecentToolRequest = nullptr;
      for( auto jt = std::next(it); jt != convo->responses.rend(); ++jt )
      {
        if( (*jt)->type() == LlmInteractionTurn::Type::ToolCall )
        {
          mostRecentToolRequest = std::dynamic_pointer_cast<LlmToolRequest>( *jt );
          break;
        }
      }

      // Verify invocationIds match between tool results and requests
      if( mostRecentToolRequest && toolResults )
      {
        const std::vector<LlmToolCall> &requestCalls = mostRecentToolRequest->toolCalls();
        const std::vector<LlmToolCall> &resultCalls = toolResults->toolCalls();

        if( requestCalls.size() == resultCalls.size() )
        {
          bool allMatch = true;
          for( const LlmToolCall &resultCall : resultCalls )
          {
            bool foundMatch = false;
            for( const LlmToolCall &requestCall : requestCalls )
            {
              if( resultCall.invocationId == requestCall.invocationId )
              {
                foundMatch = true;
                break;
              }
            }
            if( !foundMatch )
            {
              allMatch = false;
              cerr << "Warning: Tool result invocationId '" << resultCall.invocationId
                   << "' does not match any tool request invocationId" << endl;
              assert( foundMatch && "Tool result invocationId should match a tool request" );
              break;
            }
          }//for( result calls )

          if( allMatch || (requestCalls.size() == 1) )
            (*it)->setRawContent( followupRequest.dump() );
        }else
        {
          cerr << "Warning: Tool request count (" << requestCalls.size()
               << ") != tool result count (" << resultCalls.size() << ")" << endl;
          assert( requestCalls.size() == resultCalls.size() );
        }
      }else
      {
        // No tool request found or cast failed - store raw content anyway
        (*it)->setRawContent( followupRequest.dump() );
      }

      break; // Only update the most recent tool result
    }
  }//for( walk backwards through responses )
}//void storeRawContentForToolResults(...)


void LlmInterface::setupJavaScriptBridge() {
  auto app = Wt::WApplication::instance();
  if (!app) {
    cout << "Warning: No WApplication instance for JavaScript bridge setup" << endl;
    return;
  }
  
  //cout << "Setting up JavaScript bridge for HTTP requests..." << endl;
  
  // Set up shared JavaScript functions (only injected once per page)
  // The makeRequest/llmHttpRequest functions are shared across all LlmInterface instances;
  // each instance registers its own callback in window.llmResponseCallbacks[instanceId].
  string sharedJsCode = R"(
    if ( !window.llmHttpRequest ) {
      window.requestRetryState = {};
      window.llmResponseCallbacks = {};

      // Internal function to make a request (supports retry)
      function makeRequest(endpoint, requestJsonString, bearerToken, requestId, instanceId, isRetry, rateLimitRetryCount) {
        var MAX_RATE_LIMIT_RETRIES = 6;
        var DEFAULT_RATE_LIMIT_WAIT_MS = 5000;
        var MAX_RATE_LIMIT_WAIT_MS = 120000;
        var currentRateLimitRetry = rateLimitRetryCount || 0;

        var startTime = Date.now();
        var requestObj = null;

        try {
          requestObj = JSON.parse(requestJsonString);
        } catch(e) {
          console.error('Failed to parse request JSON:', e);
        }

        // Log request diagnostics
        if (!isRetry) {
          console.log('=== LLM Request ID', requestId, 'instance', instanceId, '===');
          console.log('Endpoint:', endpoint);
          console.log('Submit time:', new Date().toISOString());
          console.log('Request size:', requestJsonString.length, 'bytes');

          if (requestObj) {
            console.log('Model:', requestObj.model || 'not specified');
            console.log('Messages count:', requestObj.messages ? requestObj.messages.length : 0);
            console.log('Tools count:', requestObj.tools ? requestObj.tools.length : 0);

            // Calculate total conversation tokens (rough estimate)
            var totalChars = 0;
            if (requestObj.messages) {
              requestObj.messages.forEach(function(msg) {
                if (msg.content && typeof msg.content === 'string') {
                  totalChars += msg.content.length;
                }
              });
            }
            console.log('Estimated input chars:', totalChars);

            // Log truncated request JSON
            var jsonStr = JSON.stringify(requestObj);
            if (jsonStr.length > 80) {
              var truncated = jsonStr.substring(0, 70) + '...' + jsonStr.substring(jsonStr.length - 10);
              console.log('Request JSON (truncated):', truncated);
            } else {
              console.log('Request JSON:', jsonStr);
            }
          } else {
            var truncated = requestJsonString.length > 80
              ? requestJsonString.substring(0, 70) + '...' + requestJsonString.substring(requestJsonString.length - 10)
              : requestJsonString;
            console.log('Request JSON (raw, truncated):', truncated);
          }
        }

        var headers = {
          'Content-Type': 'application/json'
        };

        if (bearerToken && bearerToken.trim() !== '') {
          headers['Authorization'] = 'Bearer ' + bearerToken;

          // Some providers will reject the message if x-access-tokens is set, so only set it for the one I know needs it
          if (endpoint.indexOf('.ai/server/query') !== -1) {
            headers['x-access-tokens'] = bearerToken;
          }
        }

        // Create AbortController for timeout handling
        var controller = new AbortController();
        var timeoutMs = 300000; // 5 minutes
        var timeoutId = setTimeout(function() {
          var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
          console.error('=== LLM Request TIMEOUT (ID ' + requestId + ', instance ' + instanceId + (isRetry ? ', RETRY' : '') + ') ===');
          console.error('Elapsed time: ' + elapsed + 's');
          console.error('Timeout limit: ' + (timeoutMs / 1000) + 's');
          console.error('Endpoint: ' + endpoint);
          controller.abort();
        }, timeoutMs);

        if (!isRetry) {
          console.log('Timeout set for:', (timeoutMs / 1000) + 's');
          console.log('Initiating fetch...');
          //console.log('Sending request:', requestJsonString);
        }

        var httpStatus = 0;
        var retryAfterHeader = null;

        fetch(endpoint, {
          method: 'POST',
          headers: headers,
          body: requestJsonString,
          signal: controller.signal
        })
        .then(function(response) {
          var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
          console.log('=== LLM Response received (ID ' + requestId + ', instance ' + instanceId + (isRetry ? ', RETRY' : '') + ') ===');
          console.log('Status:', response.status, response.statusText);
          console.log('Elapsed time:', elapsed + 's');
          httpStatus = response.status;
          retryAfterHeader = response.headers.get('Retry-After');
          return response.text();
        })
        .then(function(responseText) {
          var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);
          console.log('Response parsed, size:', responseText.length, 'bytes');
          console.log('Total round-trip time:', elapsed + 's');
          //console.log('Response content', responseText );

          // Clear the timeout since we got a response
          clearTimeout(timeoutId);

          // Check for rate limit (HTTP 429 or error code) and auto-retry
          console.log('Rate limit check: httpStatus=' + httpStatus + ', rateLimitRetry=' + currentRateLimitRetry);
          var isRateLimited = (httpStatus === 429);
          var rateLimitBodyObj = null;
          if (!isRateLimited && currentRateLimitRetry < MAX_RATE_LIMIT_RETRIES) {
            try {
              rateLimitBodyObj = JSON.parse(responseText);
              if (rateLimitBodyObj && rateLimitBodyObj.error
                  && rateLimitBodyObj.error.code === 'rate_limit_exceeded') {
                isRateLimited = true;
              }
            } catch (e) { /* not JSON, not a rate limit error */ }
          }

          if (isRateLimited && currentRateLimitRetry < MAX_RATE_LIMIT_RETRIES) {
            var waitMs = DEFAULT_RATE_LIMIT_WAIT_MS;

            // Try Retry-After header
            if (retryAfterHeader) {
              var headerSeconds = parseFloat(retryAfterHeader);
              if (!isNaN(headerSeconds) && headerSeconds > 0) {
                waitMs = Math.ceil(headerSeconds * 1000);
              }
            }

            // Try extracting wait time from response body message (overrides header if found)
            // Handles both "in 3.363s" and "in 242ms" formats
            try {
              if (!rateLimitBodyObj) {
                rateLimitBodyObj = JSON.parse(responseText);
              }
              var errMsg = '';
              if (rateLimitBodyObj && rateLimitBodyObj.error && rateLimitBodyObj.error.message) {
                errMsg = rateLimitBodyObj.error.message;
              }
              var msMatch = errMsg.match(/try again in (\d+\.?\d*)\s*ms/i);
              var secMatch = errMsg.match(/try again in (\d+\.?\d*)\s*s(?!.*ms)/i);
              if (msMatch) {
                var parsedMs = Math.ceil(parseFloat(msMatch[1]));
                if (parsedMs > 0) {
                  waitMs = parsedMs;
                }
              } else if (secMatch) {
                var parsedMs = Math.ceil(parseFloat(secMatch[1]) * 1000);
                if (parsedMs > 0) {
                  waitMs = parsedMs;
                }
              }
            } catch (e) {
              // JSON parse failed, use header value or default
            }

            // Clamp to max wait time and add buffer (API-reported times are minimums)
            if (waitMs > MAX_RATE_LIMIT_WAIT_MS) {
              waitMs = MAX_RATE_LIMIT_WAIT_MS;
            }
            waitMs += 1500;

            var nextRetry = currentRateLimitRetry + 1;
            console.log('=== Rate limit hit (ID ' + requestId + ', instance ' + instanceId
              + ') - auto-retry ' + nextRetry + '/' + MAX_RATE_LIMIT_RETRIES
              + ' after ' + (waitMs / 1000).toFixed(1) + 's ===');

            setTimeout(function() {
              console.log('=== Launching rate-limit retry ' + nextRetry + ' for request ID ' + requestId + ' ===');
              makeRequest(endpoint, requestJsonString, bearerToken, requestId, instanceId, isRetry, nextRetry);
            }, waitMs);

            return;
          }

          // Check if we're a retry and original already completed
          if (isRetry && window.requestRetryState[requestId] && window.requestRetryState[requestId].originalCompleted) {
            console.log('Retry completed but original already succeeded - ignoring retry response');
            return;
          }

          // Mark as completed
          if (window.requestRetryState[requestId]) {
            window.requestRetryState[requestId].originalCompleted = true;
          }

          // Route response to the correct LlmInterface instance's callback
          var cb = window.llmResponseCallbacks[instanceId];
          if (cb) {
  console.log( 'Sent back for requestId:', requestId);
            cb(responseText, requestId);
          } else {
            console.error('No llmResponseCallback registered for instance', instanceId);
          }
        })
        .catch(function(error) {
          var elapsed = ((Date.now() - startTime) / 1000).toFixed(1);

          // Clear the timeout since we got an error
          clearTimeout(timeoutId);

          var errorType = error.name === 'AbortError' ? 'timeout_error' : 'network_error';

          // If this is the original request (not a retry) and it's a retryable error, attempt retry
          if (!isRetry && (errorType === 'timeout_error' || errorType === 'network_error')) {
            console.error('=== LLM Request ERROR (ID ' + requestId + ', instance ' + instanceId + ') - attempting retry ===');
            console.error('Error type:', errorType);

            // Initialize retry state
            if (!window.requestRetryState[requestId]) {
              window.requestRetryState[requestId] = { isRetrying: true, originalCompleted: false };
            }

            // Launch retry after small delay (100ms)
            setTimeout(function() {
              console.log('=== Launching retry for request ID ' + requestId + ' ===');
              makeRequest(endpoint, requestJsonString, bearerToken, requestId, instanceId, true);
            }, 100);

            // Don't send error to C++ yet - wait for retry
            return;
          }

          // This is a retry that also failed, OR original completed during retry
          if (window.requestRetryState[requestId] && window.requestRetryState[requestId].originalCompleted) {
            console.log('Retry failed but original already succeeded - ignoring');
            delete window.requestRetryState[requestId];
            return;
          }

          console.error('=== LLM Request ERROR (ID ' + requestId + ', instance ' + instanceId + (isRetry ? ', after retry' : '') + ') ===');
          console.error('Error type:', errorType);
          console.error('Error message:', error.message);
          console.error('Elapsed time:', elapsed + 's');

          // Clean up retry state
          delete window.requestRetryState[requestId];

          // Send error response to C++
          var errorResponse = JSON.stringify({
            error: {
              message: 'HTTP request failed: ' + error.message,
              type: errorType,
              elapsed_seconds: parseFloat(elapsed),
              request_id: requestId,
              retry_attempted: isRetry
            }
          });

          var cb = window.llmResponseCallbacks[instanceId];
          if (cb) {
            cb(errorResponse, requestId);
          } else {
            console.error('No llmResponseCallback registered for instance', instanceId);
          }
        });
      }

      // Public entry point - always starts as non-retry
      window.llmHttpRequest = function(endpoint, requestJsonString, bearerToken, requestId, instanceId) {
        makeRequest(endpoint, requestJsonString, bearerToken, requestId, instanceId, false);
      };

      console.log('LLM JavaScript bridge initialized with automatic retry support');
    }
  )";

  app->doJavaScript( sharedJsCode );

  // Register this instance's callback in the JavaScript registry
  const string callbackJs =
    "if (!window.llmResponseCallbacks) { window.llmResponseCallbacks = {}; }\n"
    "window.llmResponseCallbacks['" + m_instanceId + "'] = function(response, requestId) { "
    "" + m_responseSignal->createCall("response", "requestId") + ";"
    "};";

  app->doJavaScript( callbackJs );
  
  //cout << "JavaScript bridge setup complete" << endl;
}

void LlmInterface::handleJavaScriptResponse(std::string response, int requestId)
{
  if( m_debug_stream )
  {
    (*m_debug_stream) << "\n\n=== Recieved requestId=" << requestId << " ===" << endl;
    SpecUtils::trim( response );
    
    std::string responsePreview = response;
    SpecUtils::ireplace_all( responsePreview, "\n", " ");
    SpecUtils::ireplace_all( responsePreview, "\r", "");
    if( responsePreview.length() > 300 )
      responsePreview = responsePreview.substr(0, 300) + "...";
    (*m_debug_stream) << "Response: " << response << endl;
  }

  std::shared_ptr<LlmInteraction> convo;
  
  try
  {
    // Find and remove the pending request
    PendingRequest pendingRequest;
    if( m_pendingRequests.find(requestId) != m_pendingRequests.end() )
    {
      pendingRequest = std::move(m_pendingRequests[requestId]);
      m_pendingRequests.erase(requestId);
    }else
    {
      cerr << "Got response that didnt have pending request: requestId=" << requestId << ", response='" << response << "'" << endl;
      assert( 0 );
      return;
    }
    
    convo = pendingRequest.conversation.lock();
    if( !convo )
    {
      cerr << "For JavaScript response, found original request, but conversation is nullptr, so ending this conversation." << endl;
      assert( convo );
      
      return;
    }//if( !convo )
    
    // Check for errors first
    if( m_debug_stream )
      (*m_debug_stream) << "--- about to parse response to json --- " << endl;
    
    json responseJson = parseLenientJson( response );
    
    if( m_debug_stream )
      (*m_debug_stream) << "--- done parsing response to json --- " << endl;


    // In addition to checking for error, check status == 400...
    // We get here
    // "[{\n  \"error\": {\n    \"code\": 401,\n    \"message\": \"Request had invalid authentication credentials. Expected OAuth 2 access token, login cookie or other valid authentication credential. See https://developers.google.com/identity/sign-in/web/devconsole-project.\",\n    \"status\": \"UNAUTHENTICATED\",\n    \"details\": [\n      {\n        \"@type\": \"type.googleapis.com/google.rpc.ErrorInfo\",\n        \"reason\": \"ACCESS_TOKEN_TYPE_UNSUPPORTED\",\n        \"metadata\": {\n          \"method\": \"google.cloud.aiplatform.v1beta1.PredictionService.ChatCompletions\",\n          \"service\": \"aiplatform.googleapis.com\"\n        }\n      }\n    ]\n  }\n}\n]
    if( responseJson.is_array() && (responseJson.size() == 1) && responseJson[0].is_object() )
    {
      responseJson = responseJson[0];
    }//if( responseJson.is_array() && responseJson.size() > 0 )

    // TODO: we need to make sure `handleApiResponse` gives an error if it doesnt recognize the format of the result - it currently doesnt


    if( (responseJson.contains("error") && !responseJson["error"].is_null())
       || (responseJson.contains("choices") && responseJson["choices"].contains("error")) )
    {
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
      string now_str = SpecUtils::to_iso_string( now );
      const string::size_type period_pos = now_str.find('.');
      if( period_pos != string::npos )
        now_str = now_str.substr(0,period_pos);
        
      string debug_name = "llm_request_with_error_id" + std::to_string(requestId) + "_" + now_str + ".json";
      string debug_result = "llm_result_with_error_id" + std::to_string(requestId) + "_" + now_str + ".json";
#ifdef _WIN32
      const std::wstring wdebug_name = SpecUtils::convert_from_utf8_to_utf16(debug_name);
      std::ofstream output_request_json( wdebug_name.c_str(), ios::binary | ios::out );
      const std::wstring wdebug_result = SpecUtils::convert_from_utf8_to_utf16(debug_result);
      std::ofstream output_result_json( wdebug_result.c_str(), ios::binary | ios::out );
#else
      std::ofstream output_request_json( debug_name.c_str(), ios::binary | ios::out);
      std::ofstream output_result_json( debug_result.c_str(), ios::binary | ios::out );
#endif
      if( output_request_json )
        output_request_json << pendingRequest.requestJson.dump(2);
      if( output_result_json )
        output_result_json << response;
      cout << "\nLLM request error: wrote request input and output to '"
      << debug_name << "', and '" << debug_result << "', respectively."
      << endl << endl;
      if( m_debug_stream )
        (*m_debug_stream) << "\nLLM request error: wrote request input and output to '"
        << debug_name << "', and '" << debug_result << "', respectively."
        << endl << endl;
#endif //#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
      
      // Determine error type from JavaScript response
      LlmInteractionError::ErrorType errorType = LlmInteractionError::ErrorType::Unknown;
      bool retryAttempted = false;

      if( responseJson.contains("error") && responseJson["error"].is_object() )
      {
        const json &errorObj = responseJson["error"];
        /* Example `errorObj`:
         {
           "code": 401,
           "message": "Request had invalid authentication ...",
           "status": "UNAUTHENTICATED",
           "details": [{"@type": "type.googleapis.com/google.rpc.ErrorInfo", "reason": "ACCESS_TOKEN_TYPE_UNSUPPORTED", "metadata": {...}
         }
        */

        if( errorObj.contains("code")
           && ((errorObj["code"].is_number() && errorObj["code"] == 401)
               || (errorObj["code"].is_string() && errorObj["code"] == "401")) )
        {

          string msg = "Response code 401 (Unauthorized).";
          if( errorObj["message"] && errorObj["message"].is_string() )
            msg += " Message: " + errorObj["message"].get<string>();
          throw runtime_error( msg );
        }

        // Check for error type field from JavaScript
        if( errorObj.contains("type") && errorObj["type"].is_string() )
        {
          const std::string typeStr = errorObj["type"].get<std::string>();
          if( typeStr == "timeout_error" )
            errorType = LlmInteractionError::ErrorType::Timeout;
          else if( typeStr == "network_error" )
            errorType = LlmInteractionError::ErrorType::Network;
          else
            errorType = LlmInteractionError::ErrorType::LlmApi;
        }

        // Check if automatic retry was attempted
        if( errorObj.contains("retry_attempted") && errorObj["retry_attempted"].is_boolean() )
          retryAttempted = errorObj["retry_attempted"].get<bool>();
      }

      string errorMsg = "LLM API Error: " + responseJson["error"].dump(2);
      cerr << errorMsg << endl;
      if( m_debug_stream )
        (*m_debug_stream) << errorMsg << endl;

      // Add error to conversation history with type
      std::shared_ptr<LlmInteractionError> errorResponse = nullptr;
      if( m_history )
      {
        errorResponse = m_history->addErrorMessage( errorMsg, response, convo, errorType );
        if( errorResponse )
          errorResponse->setRetryAttempted( retryAttempted );
      }

      // IMPORTANT: Do NOT set finishTime for any errors - all errors pause the conversation
      // The user can either retry or continue anyway via UI buttons
      cerr << "Conversation PAUSED due to error (type: "
           << static_cast<int>(errorType) << ")" << endl;
      if( m_debug_stream )
        (*m_debug_stream) << "Conversation PAUSED due to error (type: " << static_cast<int>(errorType) << ")" << endl;

      // On error, clean up any deferred tool results that will never complete
      //  (since the response that would have triggered them never arrived).
      if( !m_deferredToolResults.empty() )
      {
        cerr << "Error response - clearing " << m_deferredToolResults.size()
             << " orphaned deferred tool results" << endl;
        if( m_debug_stream )
          (*m_debug_stream) << "Error response - clearing " << m_deferredToolResults.size()
            << " orphaned deferred tool results" << endl;
        m_deferredToolResults.clear();
      }

      // Signal that an error response was received
      if( m_pendingRequests.empty() )
      {
        cerr << "Error response - no pending requests, emitting error signal" << endl;
        if( m_debug_stream )
          (*m_debug_stream) << "Error response - no pending requests, emitting error signal" << endl;
        m_responseError.emit();
      }else
      {
        cerr << "Error response - still have " << m_pendingRequests.size()
             << " pending requests, not emitting signal yet" << endl;
        if( m_debug_stream )
          (*m_debug_stream) << "Error response - still have " << m_pendingRequests.size()
            << " pending requests, not emitting signal yet" << endl;
      }
      return;
    }
    
    // Handle summarization responses specially - apply summary, then send queued user message
    if( pendingRequest.isSummarizationRequest )
    {
      if( m_debug_stream )
        (*m_debug_stream) << "=== Processing summarization response ===" << endl;

      // Extract the summary text from the LLM response
      string summaryText;
      try
      {
        if( responseJson.contains("choices") && !responseJson["choices"].empty() )
        {
          const json &choice = responseJson["choices"][0];
          if( choice.contains("message") && choice["message"].contains("content") )
          {
            if( choice["message"]["content"].is_string() )
              summaryText = choice["message"]["content"].get<string>();
            else if( choice["message"]["content"].is_array() )
            {
              // Handle Claude content block array
              for( const json &block : choice["message"]["content"] )
              {
                if( block.value("type", "") == "text" )
                  summaryText += block.value("text", "");
              }
            }
          }
        }
      }catch( const std::exception &e )
      {
        cerr << "Error extracting summary from LLM response: " << e.what() << endl;
      }

      if( summaryText.empty() )
      {
        cerr << "Summarization response was empty, skipping summary application" << endl;
        if( m_debug_stream )
          (*m_debug_stream) << "Warning: Summarization response was empty" << endl;
      }else
      {
        if( m_debug_stream )
          (*m_debug_stream) << "Applying summary (" << summaryText.size() << " chars)" << endl;

        // Add the summary response to the summarization conversation (visible in the widget)
        m_history->addAssistantMessageWithThinking( summaryText, "", response, convo );

        // Apply the summary to replace older conversations
        m_history->applySummary( summaryText );

        // Invalidate serialization cache since history structure changed
        m_serializedCache.clear();
      }

      // Now send the queued user message
      shared_ptr<LlmInteraction> pendingUserConvo = m_summarizationPendingConvo.lock();
      m_summarizationPendingConvo.reset();

      if( pendingUserConvo )
      {
        if( m_debug_stream )
          (*m_debug_stream) << "=== Sending queued user message after summarization ===" << endl;

        json requestJson = buildMessagesArray( pendingUserConvo );
        std::pair<int,string> request_id_content = makeTrackedApiCall( requestJson, pendingUserConvo );

        if( m_debug_stream )
          (*m_debug_stream) << "Sent queued user message with request ID: " << request_id_content.first << endl;

        // Store raw JSON in the InitialRequest turn
        if( !pendingUserConvo->responses.empty() )
        {
          std::shared_ptr<LlmInteractionTurn> firstTurn = pendingUserConvo->responses.front();
          if( firstTurn && firstTurn->type() == LlmInteractionTurn::Type::InitialRequest )
            firstTurn->setRawContent( std::move(request_id_content.second) );
        }
      }else
      {
        // No pending user conversation - this is expected for manual compaction.
        // Emit conversationFinished so the UI re-enables input.
        if( m_debug_stream )
          (*m_debug_stream) << "Summarization completed (no queued user message - manual compaction)" << endl;

        if( m_pendingRequests.empty() && m_deferredToolResults.empty() )
          m_conversationFinished.emit();
      }

      return;
    }//if( summarization request )

    if( pendingRequest.isSubAgentRequest && m_debug_stream )
      (*m_debug_stream) << "=== Processing sub-agent response for: " << agentTypeToString(convo->agent_type) << " ===" << endl;

    handleApiResponse( response, convo, requestId );
  }catch( const json::parse_error &e )
  {
    string errorMsg = "Failed to parse LLM response as JSON: " + string(e.what());

    cerr << errorMsg << endl << "Raw response: " << response << endl;
    if( m_debug_stream )
      (*m_debug_stream) << errorMsg << endl << "Raw response: " << response << endl;

    if( m_history && convo )
      m_history->addErrorMessage( errorMsg, response, convo, LlmInteractionError::ErrorType::JsonParse );

    if( m_pendingRequests.empty() && m_deferredToolResults.empty() )
      m_responseError.emit();
  }catch( const std::exception &e )
  {
    string errorMsg = "Error processing LLM response: " + string(e.what());

    cout << errorMsg << endl;
    if( m_debug_stream )
      (*m_debug_stream) << errorMsg << endl;

    if( m_history && convo )
      m_history->addErrorMessage( errorMsg, response, convo, LlmInteractionError::ErrorType::Unknown );

    if( m_pendingRequests.empty() && m_deferredToolResults.empty() )
      m_responseError.emit();
  }
}

std::pair<int,std::string> LlmInterface::makeTrackedApiCall( const nlohmann::json& requestJson,
                                      std::shared_ptr<LlmInteraction> convo )
{
  assert( convo );
  
  const int requestId = m_nextRequestId++;  

  // Store the pending request
  PendingRequest pending;
  pending.requestId = requestId;
  pending.conversation = convo;
  pending.isSubAgentRequest = (convo->agent_type != AgentType::MainAgent);
#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  pending.requestJson = requestJson;
#endif

  m_pendingRequests[requestId] = pending;
  
  // Make the call with request ID tracking
  std::string request_content = makeApiCallWithId(requestJson, requestId);
  
  return make_pair(requestId,std::move(request_content));
}

int LlmInterface::sendToolResultsToLLM( std::shared_ptr<LlmInteraction> convo )
{
  // Build request with current conversation state (tool calls + results already in history)
  json followupRequest = buildMessagesArray( convo );

  // Make tracked API call to get LLM's response to the tool results
  std::pair<int,string> request_id_content = makeTrackedApiCall( followupRequest, convo );

  // Store the raw JSON being sent on the most recent ToolResult turn
  storeRawContentForToolResults( convo, followupRequest );

  return request_id_content.first;
}//int sendToolResultsToLLM( std::shared_ptr<LlmInteraction> convo )


#if( PERFORM_DEVELOPER_CHECKS )
// Test namespace to expose static functions for unit testing
namespace LlmInterfaceTests
{
  nlohmann::json lenientlyParseJson( const std::string &jsonStr )
  {
    return ::lenientlyParseJson( jsonStr );
  }

  std::string sanitizeJsonString( const std::string &jsonStr )
  {
    return ::sanitizeJsonString( jsonStr );
  }

  std::string repairIncompleteJson( const std::string &jsonStr,
                                    std::string *repairLog )
  {
    return ::repairIncompleteJson( jsonStr, repairLog );
  }
}//namespace LlmInterfaceTests
#endif


