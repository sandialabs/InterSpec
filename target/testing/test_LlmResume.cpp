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
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <Wt/WApplication.h>
#include <Wt/Test/WTestEnvironment.h>

#define BOOST_TEST_MODULE LlmResume_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/LlmConversationHistory.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

/* Tests that a conversation cut off mid-work - the LLM Assistant closed, or the session ended, while
 the LLM or its tools were still working - can be picked back up.  No LLM is called: requests go to
 the browser transport (a no-op here), and responses are fed in through
 LlmInterface::handleJavaScriptResponse(), as the browser bridge would.
 */

namespace
{
  // An enabled config with an OpenAI-compatible endpoint, so responses are in the chat-completions format.
  const char * const sm_llm_config =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<LlmConfig version=\"0\">\n"
    "  <LlmApi version=\"1\">\n"
    "    <Enabled>true</Enabled>\n"
    "    <ApiProvider active=\"true\">\n"
    "      <ApiEndpoint>http://localhost:1/v1/chat/completions</ApiEndpoint>\n"
    "      <BearerToken>x</BearerToken>\n"
    "      <Model active=\"true\"><Name>test-model</Name></Model>\n"
    "    </ApiProvider>\n"
    "  </LlmApi>\n"
    "  <McpServer version=\"1\"><Enabled>false</Enabled><BearerToken>x</BearerToken></McpServer>\n"
    "</LlmConfig>\n";


  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char ** const argv = boost::unit_test::framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
    }

    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
        {
          datadir = d;
          break;
        }
      }
    }//if( datadir.empty() )

    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path( datadir, "sandia.decay.xml" ) ),
                           "Could not find data directory; pass --datadir=" );
    InterSpec::setStaticDataDirectory( datadir );
    InterSpec::setWritableDataDirectory( datadir );

    // The agent and tool definitions come from the data directory; the provider settings from here.
    const string config_dir = SpecUtils::temp_file_name( "llm_resume_cfg", SpecUtils::temp_dir() );
    BOOST_REQUIRE( SpecUtils::create_directory( config_dir ) == 1 );
    {
      ofstream out( SpecUtils::append_path( config_dir, "llm_config.xml" ).c_str(), ios::binary | ios::out );
      out << sm_llm_config;
    }
    LlmConfig::setConfigDirectory( config_dir );
  }//void set_data_dir()


  class ResumeFixture
  {
  public:
    unique_ptr<Wt::Test::WTestEnvironment> m_env;
    InterSpecApp *m_app;
    unique_ptr<Wt::WApplication::UpdateLock> m_update_lock;
    InterSpec *m_interspec;

    /** @param create_tool Open the assistant; without it, as the stdio MCP server before its first
               prompt, or the desktop app with the assistant closed. */
    ResumeFixture( const bool create_tool = true )
      : m_env( nullptr ), m_app( nullptr ), m_interspec( nullptr )
    {
      set_data_dir();

      const string app_root = SpecUtils::lexically_normalize_path(
                                    SpecUtils::append_path( InterSpec::staticDataDirectory(), ".." ) );
      m_env.reset( new Wt::Test::WTestEnvironment( "", "", Wt::EntryPointType::Application ) );
      m_env->setAppRoot( app_root );

      m_app = new InterSpecApp( *m_env );
      m_update_lock.reset( new Wt::WApplication::UpdateLock( m_app ) );
      m_interspec = m_app->viewer();
      BOOST_REQUIRE( m_interspec );

      BOOST_REQUIRE( LlmToolGui::llmToolIsConfigured() );
      if( !create_tool )
        return;

      m_interspec->createLlmTool();
      BOOST_REQUIRE( tool() );
      BOOST_REQUIRE( tool()->llmInterface() );
    }

    ~ResumeFixture()
    {
      m_update_lock.reset();
      m_env.reset();   // deletes m_app
    }

    LlmToolGui *tool() const { return m_interspec->currentLlmTool(); }

    LlmInterface &llm() const
    {
      BOOST_REQUIRE( tool() && tool()->llmInterface() );
      return *tool()->llmInterface();
    }

    /** Shows `convos` in the tool, as when they are loaded from a file. */
    void setHistory( const vector<shared_ptr<LlmInteraction>> &convos ) const
    {
      tool()->setConversationHistory( make_shared<vector<shared_ptr<LlmInteraction>>>( convos ) );
    }
  };//class ResumeFixture


  /** The ids of the requests in flight; ids are handed out sequentially from 1. */
  vector<int> pending_request_ids( const LlmInterface &llm )
  {
    vector<int> ids;
    for( int id = 1; id < 10000; ++id )
    {
      if( llm.isRequestPending( id ) )
        ids.push_back( id );
    }
    return ids;
  }


  /** The single request in flight. */
  int pending_request_id( const LlmInterface &llm )
  {
    const vector<int> ids = pending_request_ids( llm );
    BOOST_REQUIRE_EQUAL( ids.size(), 1 );
    return ids.front();
  }


  /** A chat-completions response with the final answer `text`. */
  string final_response( const string &text )
  {
    json message;
    message["role"] = "assistant";
    message["content"] = text;

    json choice;
    choice["index"] = 0;
    choice["message"] = message;
    choice["finish_reason"] = "stop";

    json response;
    response["id"] = "test-response";
    response["choices"] = json::array( { choice } );
    response["usage"] = { {"prompt_tokens", 10}, {"completion_tokens", 5}, {"total_tokens", 15} };
    return response.dump();
  }


  LlmToolCall tool_call( const string &name, const string &id, const LlmToolCall::CallStatus status,
                         const string &content, const shared_ptr<LlmInteraction> &subAgent = nullptr )
  {
    LlmToolCall call( name, id, json::object() );
    call.status = status;
    call.content = content;
    call.sub_agent_conversation = subAgent;
    return call;
  }


  /** Adds a tool round of `call` to `convo`: its ToolCall turn, and `call` as its result. */
  void add_tool_round( LlmConversationHistory &history, const shared_ptr<LlmInteraction> &convo,
                       const LlmToolCall &call )
  {
    vector<LlmToolCall> requests{ LlmToolCall( call.toolName, call.invocationId, json::object() ) };
    history.addToolCalls( std::move( requests ), "", convo );

    vector<LlmToolCall> results{ call };
    history.addToolResults( std::move( results ), "", convo );
  }


  const LlmToolCall &last_result( const shared_ptr<LlmInteraction> &convo, const size_t index = 0 )
  {
    BOOST_REQUIRE( convo && !convo->responses.empty() );
    const LlmToolResults * const results = dynamic_cast<const LlmToolResults *>( convo->responses.back().get() );
    BOOST_REQUIRE( results && (results->toolCalls().size() > index) );
    return results->toolCalls()[index];
  }


  LlmInteraction::Interruption interruption( const chrono::system_clock::duration ago, const bool byUser )
  {
    LlmInteraction::Interruption answer;
    answer.time = chrono::system_clock::now() - ago;
    answer.cause = byUser ? LlmInteraction::Interruption::Cause::AssistantClosed
                          : LlmInteraction::Interruption::Cause::SessionEnded;
    return answer;
  }


  shared_ptr<LlmInteraction> sub_agent( const LlmInterface &llm, const AgentType type, const string &id )
  {
    shared_ptr<LlmInteraction> agent = LlmInteraction::createAgent( LlmInteraction::Type::User,
                                         "Context:\nA test.\n\nTask:\nDo the thing.", type, llm.config().get() );
    agent->conversationId = id;
    return agent;
  }

  const string sm_sub_agent_placeholder
    = "{\"status\":\"pending\",\"message\":\"Sub-agent is processing (will update when complete)\"}";


  /** A chat-completions response asking for one call of tool `name`. */
  string tool_call_response( const string &name, const string &id, const json &arguments )
  {
    json call;
    call["id"] = id;
    call["type"] = "function";
    call["function"] = { {"name", name}, {"arguments", arguments.dump()} };

    json message;
    message["role"] = "assistant";
    message["content"] = "";
    message["tool_calls"] = json::array( { call } );

    json choice;
    choice["index"] = 0;
    choice["message"] = message;
    choice["finish_reason"] = "tool_calls";

    json response;
    response["id"] = "test-tool-response";
    response["choices"] = json::array( { choice } );
    response["usage"] = { {"prompt_tokens", 10}, {"completion_tokens", 5}, {"total_tokens", 15} };
    return response.dump();
  }


  /** A one-spectrum file: a flat continuum with a peak at `peak_energy`. */
  shared_ptr<SpecMeas> make_spectrum_file( const string &filename, const double peak_energy )
  {
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );

    auto counts = make_shared<vector<float>>( 1024, 10.0f );
    const double peak_chan = cal->channel_for_energy( peak_energy );
    for( int ch = 0; ch < 1024; ++ch )
    {
      const double x = (ch + 0.5) - peak_chan;
      if( fabs(x) < 20 )
        (*counts)[ch] += static_cast<float>( 1000.0*exp( -0.5*x*x/4.0 ) );
    }

    auto m = make_shared<SpecUtils::Measurement>();
    m->set_gamma_counts( counts, 300.0f, 315.0f );
    m->set_energy_calibration( cal );
    m->set_sample_number( 1 );
    m->set_detector_name( "Ge" );
    m->set_source_type( SpecUtils::SourceType::Foreground );

    auto meas = make_shared<SpecMeas>();
    meas->add_measurement( m, false );
    meas->cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );
    meas->set_filename( filename );
    return meas;
  }//make_spectrum_file(...)


  bool history_shown( const LlmToolGui &tool, const shared_ptr<LlmInteraction> &convo )
  {
    const shared_ptr<vector<shared_ptr<LlmInteraction>>> shown = tool.getConversationHistory();
    return shown && (std::find( begin(*shown), end(*shown), convo ) != end(*shown));
  }
}//namespace


BOOST_AUTO_TEST_CASE( InterruptThenResumeRequest )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  const shared_ptr<LlmInteraction> convo = llm.sendUserMessage( "What is in my spectrum?" );
  BOOST_REQUIRE( convo );
  BOOST_CHECK( llm.isConversationActive( convo ) );

  // As when the user closes the assistant: marked to be resumed, not finished.
  llm.interruptInFlight( LlmInteraction::Interruption::Cause::AssistantClosed, false );
  BOOST_CHECK( !llm.hasActiveRequests() );
  BOOST_REQUIRE( convo->interruption.has_value() );
  BOOST_CHECK( convo->interruption->cause == LlmInteraction::Interruption::Cause::AssistantClosed );
  BOOST_CHECK( !convo->finishTime.has_value() );
  BOOST_CHECK( llm.canResumeInterrupted( convo ) );

  BOOST_CHECK( llm.resumeInterrupted( convo ) );
  BOOST_CHECK( !convo->interruption.has_value() );
  BOOST_CHECK( llm.isConversationActive( convo ) );

  llm.handleJavaScriptResponse( final_response( "Cs-137." ), pending_request_id( llm ) );
  BOOST_CHECK( !llm.hasActiveRequests() );
  BOOST_CHECK( convo->finishTime.has_value() );
  BOOST_REQUIRE( !convo->responses.empty() );
  BOOST_CHECK( convo->responses.back()->type() == LlmInteractionTurn::Type::FinalLlmResponse );
}//InterruptThenResumeRequest


BOOST_AUTO_TEST_CASE( AutoResumeOnlyWhenRecentAndNotClosedByUser )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  const auto make_interrupted = []( const LlmInteraction::Interruption &when ) -> shared_ptr<LlmInteraction> {
    LlmConversationHistory scratch;
    const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "Identify it" );
    convo->interruption = when;
    return convo;
  };

  // The session ended a moment ago (e.g., the phone app was briefly in the background): resumed.
  const shared_ptr<LlmInteraction> recent = make_interrupted( interruption( chrono::seconds( 10 ), false ) );
  fixture.setHistory( { recent } );
  fixture.tool()->resumeInterruptedConversation();
  BOOST_CHECK( !recent->interruption.has_value() );
  BOOST_CHECK( llm.isConversationActive( recent ) );
  llm.handleJavaScriptResponse( final_response( "Done." ), pending_request_id( llm ) );
  BOOST_CHECK( recent->finishTime.has_value() );

  // The user closed the assistant, just now: left for them to Resume or Dismiss.
  const shared_ptr<LlmInteraction> closed = make_interrupted( interruption( chrono::seconds( 10 ), true ) );
  fixture.setHistory( { closed } );
  fixture.tool()->resumeInterruptedConversation();
  BOOST_CHECK( closed->interruption.has_value() );
  BOOST_CHECK( !llm.hasActiveRequests() );

  // The session ended a while ago: also left for the user.
  const shared_ptr<LlmInteraction> old = make_interrupted( interruption( chrono::minutes( 10 ), false ) );
  fixture.setHistory( { old } );
  fixture.tool()->resumeInterruptedConversation();
  BOOST_CHECK( old->interruption.has_value() );
  BOOST_CHECK( !llm.hasActiveRequests() );

  // Only the newest conversation can be resumed; an older interrupted one is ended when shown.
  const shared_ptr<LlmInteraction> older = make_interrupted( interruption( chrono::seconds( 20 ), false ) );
  const shared_ptr<LlmInteraction> newer = make_interrupted( interruption( chrono::seconds( 10 ), false ) );
  newer->conversationId = older->conversationId + "_newer";
  fixture.setHistory( { older, newer } );
  BOOST_CHECK( !older->interruption.has_value() );
  BOOST_CHECK( older->finishTime.has_value() );
  BOOST_CHECK( newer->interruption.has_value() );
}//AutoResumeOnlyWhenRecentAndNotClosedByUser


BOOST_AUTO_TEST_CASE( ResumeToolRoundWithUnfinishedAsyncTool )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  LlmConversationHistory scratch;
  const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "Fit the peaks" );
  add_tool_round( scratch, convo, tool_call( "add_analysis_peak", "call_1", LlmToolCall::CallStatus::Pending,
                                             "{\"status\": \"pending\"}" ) );
  convo->interruption = interruption( chrono::seconds( 5 ), false );
  fixture.setHistory( { convo } );

  BOOST_CHECK( llm.resumeInterrupted( convo ) );

  // The tool's run was lost, so the LLM is told it was interrupted, and can call it again.
  const LlmToolCall &call = last_result( convo );
  BOOST_CHECK( call.status == LlmToolCall::CallStatus::Error );
  BOOST_CHECK( call.content.find( "interrupted" ) != string::npos );

  llm.handleJavaScriptResponse( final_response( "The fit was interrupted." ), pending_request_id( llm ) );
  BOOST_CHECK( convo->finishTime.has_value() );
  BOOST_CHECK( !llm.hasActiveRequests() );
}//ResumeToolRoundWithUnfinishedAsyncTool


BOOST_AUTO_TEST_CASE( ResumeWaitsOnAsyncToolStillRunning )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  // Interrupted within this session (another spectrum was shown), while the tool was still running.
  LlmConversationHistory scratch;
  const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "Fit the peaks" );
  LlmToolCall running = tool_call( "add_analysis_peak", "call_1", LlmToolCall::CallStatus::Pending,
                                   "{\"status\": \"pending\"}" );
  running.asyncRunning = true;
  add_tool_round( scratch, convo, running );
  convo->interruption = interruption( chrono::seconds( 5 ), true );
  convo->interruption->cause = LlmInteraction::Interruption::Cause::SpectrumChanged;
  fixture.setHistory( { convo } );

  // Resuming waits for the tool, rather than telling the LLM it was lost - so it would run it again.
  BOOST_CHECK( llm.resumeInterrupted( convo ) );
  BOOST_CHECK( last_result( convo ).status == LlmToolCall::CallStatus::Pending );
  BOOST_CHECK( llm.isConversationActive( convo ) );
  BOOST_CHECK( pending_request_ids( llm ).empty() );

  // As the tool's callback does: fill in the result, then hand it to the resumed round, which sends it.
  LlmToolResults * const results = dynamic_cast<LlmToolResults *>( convo->responses.back().get() );
  BOOST_REQUIRE( results );
  results->toolCalls()[0].status = LlmToolCall::CallStatus::Success;
  results->toolCalls()[0].content = "{\"success\": true}";
  results->toolCalls()[0].asyncRunning = false;
  llm.resolveResumedAsyncCall( convo, "call_1" );
  llm.resolveResumedAsyncCall( convo, "call_1" );  // A repeat does nothing

  llm.handleJavaScriptResponse( final_response( "The peak is fit." ), pending_request_id( llm ) );
  BOOST_CHECK( convo->finishTime.has_value() );
  BOOST_CHECK( !llm.hasActiveRequests() );
}//ResumeWaitsOnAsyncToolStillRunning


BOOST_AUTO_TEST_CASE( ResumeSubAgent )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  // The main agent's tool round waits on a sub-agent, which was waiting on the LLM.
  LlmConversationHistory scratch;
  const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "What activity?" );
  const shared_ptr<LlmInteraction> agent = sub_agent( llm, AgentType::ActivityFit, "subagent_ActivityFit_1" );
  add_tool_round( scratch, convo, tool_call( "invoke_ActivityFit", "call_1", LlmToolCall::CallStatus::Pending,
                                             sm_sub_agent_placeholder, agent ) );
  convo->interruption = interruption( chrono::seconds( 5 ), false );
  agent->interruption = convo->interruption;
  fixture.setHistory( { convo } );

  BOOST_CHECK( llm.resumeInterrupted( convo ) );
  BOOST_CHECK( !agent->interruption.has_value() );
  BOOST_CHECK( llm.isConversationActive( agent ) );
  BOOST_CHECK( llm.isConversationActive( convo ) );  // its round waits on the sub-agent

  // The sub-agent's answer is handed back to the main agent, whose results then go to the LLM.
  llm.handleJavaScriptResponse( final_response( "About 5 uCi." ), pending_request_id( llm ) );
  BOOST_CHECK( agent->finishTime.has_value() );
  const LlmToolCall &call = last_result( convo );
  BOOST_CHECK( call.status == LlmToolCall::CallStatus::Success );
  BOOST_CHECK( call.content.find( "About 5 uCi." ) != string::npos );
  BOOST_CHECK( !convo->finishTime.has_value() );

  llm.handleJavaScriptResponse( final_response( "It is about 5 uCi." ), pending_request_id( llm ) );
  BOOST_CHECK( convo->finishTime.has_value() );
  BOOST_CHECK( !llm.hasActiveRequests() );
}//ResumeSubAgent


BOOST_AUTO_TEST_CASE( ResumeWithPausedSubAgentStopsRound )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  // A sub-agent paused on an error, which the user had not yet chosen to Retry or Cancel.
  LlmConversationHistory scratch;
  const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "What activity?" );
  const shared_ptr<LlmInteraction> agent = sub_agent( llm, AgentType::ActivityFit, "subagent_ActivityFit_2" );
  scratch.addErrorMessage( "Rate limited", "", agent );
  add_tool_round( scratch, convo, tool_call( "invoke_ActivityFit", "call_1", LlmToolCall::CallStatus::Pending,
                                             sm_sub_agent_placeholder, agent ) );
  convo->interruption = interruption( chrono::seconds( 5 ), false );
  fixture.setHistory( { convo } );

  BOOST_CHECK( llm.resumeInterrupted( convo ) );

  // Back to waiting on the user, as before it was cut off.
  BOOST_CHECK( llm.isConversationStopped( convo ) );
  BOOST_CHECK( !llm.hasActiveRequests() );
  BOOST_CHECK( !convo->finishTime.has_value() );
  BOOST_CHECK( !agent->finishTime.has_value() );
}//ResumeWithPausedSubAgentStopsRound


BOOST_AUTO_TEST_CASE( ContinuingPausedSubAgentResumesItsConversation )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  // As above, but the user goes to the paused agent's Retry (or Cancel Agent), not the conversation's
  //  Resume - which, done after, would send the agent's request a second time.
  LlmConversationHistory scratch;
  const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "What activity?" );
  const shared_ptr<LlmInteraction> agent = sub_agent( llm, AgentType::ActivityFit, "subagent_ActivityFit_5" );
  scratch.addErrorMessage( "Rate limited", "", agent );
  add_tool_round( scratch, convo, tool_call( "invoke_ActivityFit", "call_1", LlmToolCall::CallStatus::Pending,
                                             sm_sub_agent_placeholder, agent ) );
  convo->interruption = interruption( chrono::seconds( 5 ), true );
  fixture.setHistory( { convo } );

  BOOST_CHECK( llm.resumeInterruptedFor( agent ) );
  BOOST_CHECK( !convo->interruption.has_value() );
  BOOST_CHECK( !llm.canResumeInterrupted( convo ) );
  BOOST_CHECK( llm.isConversationStopped( convo ) );
  BOOST_CHECK( !llm.hasActiveRequests() );
  BOOST_CHECK( llm.resumeInterruptedFor( agent ) );  // Already waited on

  // One no conversation is waiting on can not be continued.
  const shared_ptr<LlmInteraction> stray = sub_agent( llm, AgentType::ActivityFit, "subagent_ActivityFit_6" );
  BOOST_CHECK( !llm.resumeInterruptedFor( stray ) );
}//ContinuingPausedSubAgentResumesItsConversation


BOOST_AUTO_TEST_CASE( ResumeRestoresSubAgentStateMachine )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  const shared_ptr<LlmInteraction> agent = sub_agent( llm, AgentType::NuclideId, "subagent_NuclideId_1" );
  BOOST_REQUIRE_MESSAGE( agent->state_machine, "The NuclideId agent is expected to have a state machine" );

  // A state other than the one it starts in.
  string newState;
  for( const auto &state : agent->state_machine->allStates() )
  {
    if( state.first != agent->state_machine->getInitialState() )
    {
      newState = state.first;
      break;
    }
  }
  BOOST_REQUIRE( !newState.empty() );

  // Read back from a file, it has no state machine; it had moved to `newState`.
  agent->state_machine.reset();

  LlmConversationHistory scratch;
  json transition;
  transition["success"] = true;
  transition["new_state"] = newState;
  add_tool_round( scratch, agent, tool_call( "set_workflow_state", "call_s", LlmToolCall::CallStatus::Success,
                                             transition.dump() ) );

  const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "What is it?" );
  add_tool_round( scratch, convo, tool_call( "invoke_NuclideId", "call_1", LlmToolCall::CallStatus::Pending,
                                             sm_sub_agent_placeholder, agent ) );
  convo->interruption = interruption( chrono::seconds( 5 ), false );
  agent->interruption = convo->interruption;
  fixture.setHistory( { convo } );

  BOOST_CHECK( llm.resumeInterrupted( convo ) );
  BOOST_REQUIRE( agent->state_machine );
  BOOST_CHECK_EQUAL( agent->state_machine->getCurrentState(), newState );
  BOOST_CHECK( llm.isConversationActive( agent ) );
}//ResumeRestoresSubAgentStateMachine


BOOST_AUTO_TEST_CASE( AbandonInterrupted )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  LlmConversationHistory scratch;
  const shared_ptr<LlmInteraction> convo = scratch.addUserMessageToMainConversation( "What activity?" );
  const shared_ptr<LlmInteraction> agent = sub_agent( llm, AgentType::ActivityFit, "subagent_ActivityFit_3" );
  {
    vector<LlmToolCall> requests{ LlmToolCall( "invoke_ActivityFit", "call_1", json::object() ),
                                  LlmToolCall( "add_analysis_peak", "call_2", json::object() ) };
    scratch.addToolCalls( std::move( requests ), "", convo );
    vector<LlmToolCall> results{
      tool_call( "invoke_ActivityFit", "call_1", LlmToolCall::CallStatus::Pending, sm_sub_agent_placeholder, agent ),
      tool_call( "add_analysis_peak", "call_2", LlmToolCall::CallStatus::Pending, "{\"status\": \"pending\"}" )
    };
    scratch.addToolResults( std::move( results ), "", convo );
  }
  convo->interruption = interruption( chrono::seconds( 5 ), true );
  agent->interruption = convo->interruption;
  fixture.setHistory( { convo } );

  // Dismissed: no placeholder is left to be sent to the LLM with a later question.
  llm.abandonInterrupted( convo );
  BOOST_CHECK( !convo->interruption.has_value() );
  BOOST_CHECK( !agent->interruption.has_value() );
  BOOST_CHECK( convo->finishTime.has_value() );
  BOOST_CHECK( agent->finishTime.has_value() );
  BOOST_CHECK( last_result( convo, 0 ).status == LlmToolCall::CallStatus::Error );
  BOOST_CHECK( last_result( convo, 1 ).status == LlmToolCall::CallStatus::Error );
  BOOST_CHECK( !llm.canResumeInterrupted( convo ) );

  // Asking a new question also ends an interrupted conversation.
  LlmConversationHistory scratch2;
  const shared_ptr<LlmInteraction> left = scratch2.addUserMessageToMainConversation( "Left hanging" );
  left->interruption = interruption( chrono::seconds( 5 ), true );
  fixture.setHistory( { left } );

  const shared_ptr<LlmInteraction> next = llm.sendUserMessage( "Something else" );
  BOOST_CHECK( !left->interruption.has_value() );
  BOOST_CHECK( left->finishTime.has_value() );
  BOOST_CHECK( llm.isConversationActive( next ) );
}//AbandonInterrupted


BOOST_AUTO_TEST_CASE( CompactionIsNotResumed )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();

  // Its summarization prompt is not kept, so the request can not be sent again.
  LlmConversationHistory scratch;
  const shared_ptr<LlmInteraction> compaction
    = scratch.addSystemMessageToMainConversation( "[Context Compaction] Summarizing conversation history..." );
  compaction->interruption = interruption( chrono::seconds( 5 ), false );

  BOOST_CHECK( !llm.canResumeInterrupted( compaction ) );
  BOOST_CHECK( !llm.resumeInterrupted( compaction ) );
  BOOST_CHECK( !compaction->interruption.has_value() );
  BOOST_CHECK( compaction->finishTime.has_value() );
  BOOST_CHECK( !llm.hasActiveRequests() );
}//CompactionIsNotResumed


BOOST_AUTO_TEST_CASE( ClosingToolOrSessionInterrupts )
{
  shared_ptr<LlmInteraction> closedOn, sessionEnded;

  {// Begin the user closing the tool
    ResumeFixture fixture;
    closedOn = fixture.llm().sendUserMessage( "Closed on" );
    fixture.m_interspec->handleLlmToolClose();
    BOOST_CHECK( !fixture.tool() );
  }// End the user closing the tool

  BOOST_REQUIRE( closedOn->interruption.has_value() );
  BOOST_CHECK( closedOn->interruption->cause == LlmInteraction::Interruption::Cause::AssistantClosed );
  BOOST_CHECK( !closedOn->finishTime.has_value() );

  {// Begin the session ending with the tool open
    ResumeFixture fixture;
    sessionEnded = fixture.llm().sendUserMessage( "Session ended" );
  }// End the session ending with the tool open

  BOOST_REQUIRE( sessionEnded->interruption.has_value() );
  BOOST_CHECK( sessionEnded->interruption->cause == LlmInteraction::Interruption::Cause::SessionEnded );
  BOOST_CHECK( !sessionEnded->finishTime.has_value() );
}//ClosingToolOrSessionInterrupts


BOOST_AUTO_TEST_CASE( SpectrumSwitchInterruptsWork )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();
  LlmToolGui &tool = *fixture.tool();

  const shared_ptr<SpecMeas> first = make_spectrum_file( "first.n42", 661.7 );
  const shared_ptr<SpecMeas> second = make_spectrum_file( "second.n42", 1460.8 );
  fixture.m_interspec->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );

  // The main agent hands the work to a sub-agent, which is waiting on the LLM.
  const shared_ptr<LlmInteraction> convo = llm.sendUserMessage( "What is in my spectrum?" );
  BOOST_REQUIRE( convo );
  const json args = { {"context", "A test spectrum."}, {"task", "Fit the activity."} };
  llm.handleJavaScriptResponse( tool_call_response( "invoke_ActivityFit", "call_1", args ), pending_request_id( llm ) );
  BOOST_REQUIRE( llm.hasActiveRequests() );
  const shared_ptr<LlmInteraction> agent = last_result( convo ).sub_agent_conversation;
  BOOST_REQUIRE( agent );
  BOOST_CHECK( llm.isConversationActive( agent ) );

  // Showing the same conversations again (e.g., a time-range change) leaves the work alone.
  tool.setConversationHistory( tool.getConversationHistory() );
  BOOST_CHECK( llm.isConversationActive( agent ) );
  BOOST_CHECK( !convo->interruption.has_value() );

  // Showing another file interrupts it, rather than ending it.
  fixture.m_interspec->setSpectrum( second, {}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_CHECK( !llm.hasActiveRequests() );
  BOOST_CHECK( !history_shown( tool, convo ) );
  BOOST_REQUIRE( convo->interruption.has_value() );
  BOOST_CHECK( convo->interruption->cause == LlmInteraction::Interruption::Cause::SpectrumChanged );
  BOOST_CHECK( agent->interruption.has_value() );
  BOOST_CHECK( !convo->finishTime.has_value() );
  BOOST_CHECK( !agent->finishTime.has_value() );
  BOOST_CHECK( last_result( convo ).status == LlmToolCall::CallStatus::Pending );

  // Back on the first file, it is offered for resuming - not resumed on its own, as the user moved away.
  fixture.m_interspec->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_REQUIRE( history_shown( tool, convo ) );
  BOOST_CHECK( convo->interruption.has_value() );
  BOOST_CHECK( llm.canResumeInterrupted( convo ) );
  tool.resumeInterruptedConversation();
  BOOST_CHECK( !llm.hasActiveRequests() );

  // Resumed, the sub-agent picks back up, and hands its answer to the main agent.
  BOOST_CHECK( llm.resumeInterrupted( convo ) );
  BOOST_CHECK( llm.isConversationActive( agent ) );
  llm.handleJavaScriptResponse( final_response( "Cs-137." ), pending_request_id( llm ) );
  BOOST_CHECK( agent->finishTime.has_value() );
  BOOST_CHECK( last_result( convo ).status == LlmToolCall::CallStatus::Success );

  llm.handleJavaScriptResponse( final_response( "It is Cs-137." ), pending_request_id( llm ) );
  BOOST_CHECK( convo->finishTime.has_value() );
  BOOST_CHECK( !llm.hasActiveRequests() );
}//SpectrumSwitchInterruptsWork


BOOST_AUTO_TEST_CASE( SpectrumSwitchInterruptsStoppedRound )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();
  LlmToolGui &tool = *fixture.tool();

  const shared_ptr<SpecMeas> first = make_spectrum_file( "first.n42", 661.7 );
  const shared_ptr<SpecMeas> second = make_spectrum_file( "second.n42", 1460.8 );
  fixture.m_interspec->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );

  const shared_ptr<LlmInteraction> convo = llm.sendUserMessage( "What is in my spectrum?" );
  BOOST_REQUIRE( convo );
  const json args = { {"context", "A test spectrum."}, {"task", "Fit the activity."} };
  llm.handleJavaScriptResponse( tool_call_response( "invoke_ActivityFit", "call_1", args ), pending_request_id( llm ) );
  const shared_ptr<LlmInteraction> agent = last_result( convo ).sub_agent_conversation;
  BOOST_REQUIRE( agent );

  // The user pressed Stop: the sub-agent is paused, waiting on them to Retry or Cancel Agent.
  llm.cancelAll();
  BOOST_REQUIRE( llm.isConversationStopped( convo ) );

  // Showing another file, then this one again, keeps that choice open (it used to end the agent).
  fixture.m_interspec->setSpectrum( second, {}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_CHECK( !agent->finishTime.has_value() );
  BOOST_CHECK( !convo->finishTime.has_value() );
  fixture.m_interspec->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_REQUIRE( history_shown( tool, convo ) );
  BOOST_CHECK( llm.canResumeInterrupted( convo ) );

  BOOST_CHECK( llm.resumeInterrupted( convo ) );
  BOOST_CHECK( llm.isConversationStopped( convo ) );
}//SpectrumSwitchInterruptsStoppedRound


// An MCP client starts a new conversation with load_spectrum_file's conversation "new" and no path,
//  whether or not a foreground is loaded; without "new", unloading leaves the conversation alone.
BOOST_AUTO_TEST_CASE( McpUnloadWithNewConversationClearsIt )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();
  LlmToolGui &tool = *fixture.tool();
  InterSpec * const viewer = fixture.m_interspec;
  const shared_ptr<const LlmTools::ToolRegistry> registry = llm.toolRegistry();
  BOOST_REQUIRE( registry );

  const auto ask = [&llm]( const string &prompt ) -> shared_ptr<LlmInteraction> {
    const shared_ptr<LlmInteraction> convo = llm.sendUserMessage( prompt );
    BOOST_REQUIRE( convo );
    llm.handleJavaScriptResponse( final_response( "Cs-137." ), pending_request_id( llm ) );
    BOOST_REQUIRE( convo->finishTime.has_value() );
    return convo;
  };

  const json unload_new = { {"specType", "Foreground"}, {"conversation", "new"} };

  // No foreground: just the conversation is cleared.
  BOOST_REQUIRE( !viewer->measurment( SpecUtils::SpectrumType::Foreground ) );
  const shared_ptr<LlmInteraction> unattached = ask( "What is in my spectrum?" );
  BOOST_REQUIRE( history_shown( tool, unattached ) );
  json result;
  BOOST_REQUIRE_NO_THROW( result = registry->executeTool( "load_spectrum_file", unload_new, viewer ) );
  BOOST_CHECK_EQUAL( result.value( "conversation", string() ), "new" );
  BOOST_CHECK( !history_shown( tool, unattached ) );
  BOOST_CHECK( !tool.getConversationHistory() );

  // With a foreground: unloaded, the conversation is cleared, and the file keeps its conversation.
  const shared_ptr<SpecMeas> first = make_spectrum_file( "first.n42", 661.7 );
  viewer->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );
  const shared_ptr<LlmInteraction> attached = ask( "And now?" );
  BOOST_REQUIRE_NO_THROW( result = registry->executeTool( "load_spectrum_file", unload_new, viewer ) );
  BOOST_CHECK( !viewer->measurment( SpecUtils::SpectrumType::Foreground ) );
  BOOST_CHECK( !tool.getConversationHistory() );
  const shared_ptr<vector<shared_ptr<LlmInteraction>>> saved = first->llmConversationHistory( {1} );
  BOOST_CHECK( saved && (std::find( begin(*saved), end(*saved), attached ) != end(*saved)) );

  // Without "new", unloading the foreground keeps the conversation.
  viewer->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_REQUIRE( history_shown( tool, attached ) );
  BOOST_REQUIRE_NO_THROW( registry->executeTool( "load_spectrum_file", json{ {"specType", "Foreground"} }, viewer ) );
  BOOST_CHECK( !viewer->measurment( SpecUtils::SpectrumType::Foreground ) );
  BOOST_CHECK( history_shown( tool, attached ) );

  // "resume" needs a file, and is refused before anything is unloaded.
  viewer->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_CHECK_THROW( registry->executeTool( "load_spectrum_file",
                       json{ {"specType", "Foreground"}, {"conversation", "resume"} }, viewer ), std::exception );
  BOOST_CHECK( viewer->measurment( SpecUtils::SpectrumType::Foreground ) == first );

  // Loading the displayed foreground's path again with "new" keeps that copy of the file (and its
  //  work); as for any foreground load, its unfinished work is set aside with the file (interrupted).
  const string n42_path = SpecUtils::temp_file_name( "mcp_new_convo", SpecUtils::temp_dir() ) + ".n42";
  {
    ofstream out( n42_path.c_str(), ios::binary | ios::out );
    BOOST_REQUIRE( first->write_2012_N42( out ) );
  }
  const json load_new = { {"path", n42_path}, {"specType", "Foreground"}, {"conversation", "new"} };
  BOOST_REQUIRE_NO_THROW( registry->executeTool( "load_spectrum_file", load_new, viewer ) );
  const shared_ptr<SpecMeas> loaded = viewer->measurment( SpecUtils::SpectrumType::Foreground );
  BOOST_REQUIRE( loaded && (loaded != first) );

  const shared_ptr<LlmInteraction> again = ask( "What about now?" );
  BOOST_REQUIRE_NO_THROW( registry->executeTool( "load_spectrum_file", load_new, viewer ) );
  BOOST_CHECK( viewer->measurment( SpecUtils::SpectrumType::Foreground ) == loaded );
  BOOST_CHECK( !history_shown( tool, again ) );
  BOOST_CHECK( !tool.getConversationHistory() );
  SpecUtils::remove_file( n42_path );
}//McpUnloadWithNewConversationClearsIt


// Clearing the conversation (as "new" does) ends a turn stopped on a paused sub-agent, rather than
//  leaving it pointing at a removed conversation.
BOOST_AUTO_TEST_CASE( ClearConversationEndsStoppedRound )
{
  ResumeFixture fixture;
  LlmInterface &llm = fixture.llm();
  LlmToolGui &tool = *fixture.tool();

  const shared_ptr<SpecMeas> first = make_spectrum_file( "first.n42", 661.7 );
  fixture.m_interspec->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );

  const shared_ptr<LlmInteraction> convo = llm.sendUserMessage( "What is in my spectrum?" );
  BOOST_REQUIRE( convo );
  const json args = { {"context", "A test spectrum."}, {"task", "Fit the activity."} };
  llm.handleJavaScriptResponse( tool_call_response( "invoke_ActivityFit", "call_1", args ), pending_request_id( llm ) );
  llm.cancelAll();
  BOOST_REQUIRE( llm.isConversationStopped( convo ) );

  tool.clearConversationHistory();
  BOOST_CHECK( !llm.isConversationStopped( convo ) );
  BOOST_CHECK( convo->finishTime.has_value() );
  BOOST_CHECK( !tool.getConversationHistory() );
}//ClearConversationEndsStoppedRound


namespace
{
  /** Writes a file carrying a conversation, as an N42 saved by InterSpec (or made up by someone
   else), to a new temporary file; returns its path. */
  string write_n42_with_conversation( const string &question )
  {
    const shared_ptr<LlmInteraction> saved = LlmInteraction::createUserInteractionForMainConvo( question );
    const shared_ptr<SpecMeas> source = make_spectrum_file( "saved.n42", 661.7 );
    source->setLlmConversationHistory( {1}, make_shared<vector<shared_ptr<LlmInteraction>>>( 1, saved ) );
    const string n42_path = SpecUtils::temp_file_name( "mcp_saved_convo", SpecUtils::temp_dir() ) + ".n42";
    ofstream out( n42_path.c_str(), ios::binary | ios::out );
    BOOST_REQUIRE( source->write_2012_N42( out ) );
    return n42_path;
  }//write_n42_with_conversation(...)


  /** In a session with no assistant open, loads an N42 that carries a conversation with
   load_spectrum_file's conversation `mode`, then opens the assistant; returns how many conversations
   it shows (and checks the load's result reported the same).
   */
  size_t turns_opened_after_loading( const string &mode )
  {
    ResumeFixture fixture( false );
    InterSpec * const viewer = fixture.m_interspec;
    BOOST_REQUIRE( !fixture.tool() );

    const shared_ptr<const LlmConfig> config = LlmConfig::load();
    BOOST_REQUIRE( config );
    const LlmTools::ToolRegistry registry( *config );

    const string n42_path = write_n42_with_conversation( "Saved question" );
    const json load = { {"path", n42_path}, {"specType", "Foreground"}, {"conversation", mode} };
    json result;
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "load_spectrum_file", load, viewer ) );
    SpecUtils::remove_file( n42_path );
    BOOST_REQUIRE( result.contains( "restored" ) );
    const size_t reported = result["restored"].value( "conversation_turns", size_t(999) );

    viewer->createLlmTool();
    BOOST_REQUIRE( fixture.tool() );
    const shared_ptr<vector<shared_ptr<LlmInteraction>>> shown = fixture.tool()->getConversationHistory();
    const size_t turns = shown ? shown->size() : size_t(0);
    BOOST_CHECK_MESSAGE( reported == turns, mode << ": the load reported " << reported
                         << " conversations, the assistant shows " << turns );
    return turns;
  }//turns_opened_after_loading(...)
}//namespace


// Without an assistant open yet (the stdio MCP server before its first prompt, or the desktop app
//  with the assistant closed), "keep" and "new" decide the conversation it opens with, as they do
//  for an open one: a loaded file's saved conversation is only taken up with "resume".
BOOST_AUTO_TEST_CASE( ConversationModeWithoutAssistantOpen )
{
  BOOST_CHECK_EQUAL( turns_opened_after_loading( "resume" ), 1 );
  BOOST_CHECK_EQUAL( turns_opened_after_loading( "new" ), 0 );
  BOOST_CHECK_EQUAL( turns_opened_after_loading( "keep" ), 0 );
}//ConversationModeWithoutAssistantOpen


// With the assistant closed (as in the desktop app), "keep" carries the shown foreground's
//  conversation over to the new one, as it does with the assistant open.
BOOST_AUTO_TEST_CASE( KeepWithoutAssistantOpenCarriesConversation )
{
  ResumeFixture fixture( false );
  InterSpec * const viewer = fixture.m_interspec;
  const shared_ptr<const LlmConfig> config = LlmConfig::load();
  BOOST_REQUIRE( config );
  const LlmTools::ToolRegistry registry( *config );

  const shared_ptr<LlmInteraction> current = LlmInteraction::createUserInteractionForMainConvo( "Current question" );
  const shared_ptr<SpecMeas> first = make_spectrum_file( "first.n42", 1460.8 );
  first->setLlmConversationHistory( {1}, make_shared<vector<shared_ptr<LlmInteraction>>>( 1, current ) );
  viewer->setSpectrum( first, {}, SpecUtils::SpectrumType::Foreground, {} );

  const string n42_path = write_n42_with_conversation( "Saved question" );
  const json load = { {"path", n42_path}, {"specType", "Foreground"}, {"conversation", "keep"} };
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_spectrum_file", load, viewer ) );
  SpecUtils::remove_file( n42_path );
  BOOST_REQUIRE( viewer->measurment( SpecUtils::SpectrumType::Foreground ) != first );

  viewer->createLlmTool();
  BOOST_REQUIRE( fixture.tool() );
  const shared_ptr<vector<shared_ptr<LlmInteraction>>> shown = fixture.tool()->getConversationHistory();
  BOOST_REQUIRE( shown && (shown->size() == 1) );
  BOOST_CHECK( shown->front() == current );
}//KeepWithoutAssistantOpenCarriesConversation
