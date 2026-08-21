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

#include <string>
#include <vector>
#include <memory>
#include <iostream>

#define BOOST_TEST_MODULE LlmApiProtocol_suite
#include <boost/test/included/unit_test.hpp>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmApiProtocol.h"

using namespace std;
using json = nlohmann::json;

using ApiFormat = LlmConfig::LlmApi::ApiFormat;

/* These tests pin the Anthropic prompt-cache breakpoint placement.

 Anthropic caching is an exact prefix match terminating at an explicitly marked content block, so a
 marker placed on a block that the *next* request will not reproduce can never be read back - it
 only ever pays the cache-write premium.  The state-machine reminder appended by
 appendEphemeralUserNote() is exactly such a block: it is deliberately not persisted into history.
 */
namespace
{
  /** A minimal Anthropic-format LlmApi with one provider and one active model. */
  LlmConfig::LlmApi make_anthropic_api()
  {
    LlmConfig::LlmApi api;
    api.enabled = true;

    LlmConfig::LlmApi::ModelInfo m;
    m.name = "claude-opus-4-8";
    m.maxTokens = 32000;
    m.reasoning = false;  // keep thinking off so tool_choice/temperature stay deterministic

    LlmConfig::LlmApi::ApiProvider p;
    p.apiEndpoint = "https://api.anthropic.com/v1/messages";
    p.apiFormat = ApiFormat::Anthropic;
    p.models.push_back( m );

    api.providers.push_back( p );
    return api;
  }//make_anthropic_api()


  std::vector<NormalizedTool> make_tools()
  {
    NormalizedTool t;
    t.name = "get_peaks";
    t.description = "Return the currently fit peaks.";
    t.parameters = json::parse( R"({"type":"object","properties":{}})" );
    return { t };
  }//make_tools()


  /** The text of the not-persisted state-machine reminder, as LlmInterface builds it. */
  string state_note( const string &state )
  {
    return "[System: Current workflow state is " + state
           + ". Allowed transitions: done. Call set_workflow_state when ready to transition.]";
  }


  /** Flatten `messages` into one entry per content block, as "role|<block json>", with any
   `cache_control` marker stripped.  Markers legitimately move between requests, so they must not
   participate in the prefix comparison - only the block *content* is part of the cache key.
   */
  vector<string> flatten_blocks( const json &messages )
  {
    vector<string> out;
    if( !messages.is_array() )
      return out;

    for( const json &msg : messages )
    {
      const string role = msg.value( "role", "" );
      if( !msg.contains("content") )
        continue;

      const json &content = msg["content"];
      if( content.is_string() )
      {
        json b;
        b["type"] = "text";
        b["text"] = content.get<string>();
        out.push_back( role + "|" + b.dump() );
      }else if( content.is_array() )
      {
        for( json block : content )   // by value: we strip the marker below
        {
          if( block.is_object() )
            block.erase( "cache_control" );
          out.push_back( role + "|" + block.dump() );
        }
      }
    }//for( msg : messages )

    return out;
  }//flatten_blocks(...)


  /** Index (into flatten_blocks()) one past the last `cache_control` marker, i.e. the number of
   leading blocks that the provider will store as the cached prefix.  Zero if unmarked.
   */
  size_t cached_prefix_len( const json &messages )
  {
    size_t idx = 0, lastMarked = 0;
    if( !messages.is_array() )
      return 0;

    for( const json &msg : messages )
    {
      if( !msg.contains("content") )
        continue;

      const json &content = msg["content"];
      if( content.is_string() )
      {
        idx += 1;
      }else if( content.is_array() )
      {
        for( const json &block : content )
        {
          idx += 1;
          if( block.is_object() && block.contains("cache_control") )
            lastMarked = idx;
        }
      }
    }//for( msg : messages )

    return lastMarked;
  }//cached_prefix_len(...)


  size_t count_markers( const json &messages )
  {
    size_t n = 0;
    if( !messages.is_array() )
      return n;

    for( const json &msg : messages )
    {
      if( !msg.contains("content") || !msg["content"].is_array() )
        continue;
      for( const json &block : msg["content"] )
        n += (block.is_object() && block.contains("cache_control")) ? 1 : 0;
    }

    return n;
  }//count_markers(...)


  /** Whether any marked block's text is the state-machine reminder. */
  bool marker_on_state_note( const json &messages )
  {
    if( !messages.is_array() )
      return false;

    for( const json &msg : messages )
    {
      if( !msg.contains("content") || !msg["content"].is_array() )
        continue;

      for( const json &block : msg["content"] )
      {
        if( !block.is_object() || !block.contains("cache_control") )
          continue;

        const string txt = block.value( "text", "" );
        if( txt.rfind( "[System: Current workflow state is", 0 ) == 0 )
          return true;
      }
    }//for( msg : messages )

    return false;
  }//marker_on_state_note(...)
}//namespace


BOOST_AUTO_TEST_CASE( AnthropicCacheBreakpointNotOnEphemeralNote )
{
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  json messages = json::array();
  messages.push_back( json{ {"role","user"}, {"content","Context:\nTh232 spectrum\n\nTask:\nCalibrate energy."} } );

  proto->placeConversationCacheBreakpoints( messages );
  proto->appendEphemeralUserNote( messages, state_note("gather_peaks") );

  const json req = proto->buildRequestBody( make_anthropic_api(), "You are the energy calibration agent.",
                                            messages, make_tools(), ToolChoice::Auto );

  BOOST_REQUIRE( req.contains("messages") );

  // The reminder is rebuilt (or gone) next round, so a marker on it can never be read back.
  BOOST_CHECK_MESSAGE( !marker_on_state_note( req["messages"] ),
                       "cache_control must not land on the ephemeral state-machine reminder" );

  // Something must still be marked, or the conversation is not cached at all.
  BOOST_CHECK_GT( count_markers( req["messages"] ), 0u );
}


BOOST_AUTO_TEST_CASE( AnthropicCachedPrefixSurvivesToNextRequest )
{
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  const LlmConfig::LlmApi api = make_anthropic_api();
  const vector<NormalizedTool> tools = make_tools();
  const string sysPrompt = "You are the energy calibration agent.";
  const string initialTurn = "Context:\nTh232 spectrum\n\nTask:\nCalibrate energy.";

  // --- Round 1: the sub-agent's opening turn. -------------------------------------------------
  json msgs1 = json::array();
  msgs1.push_back( json{ {"role","user"}, {"content",initialTurn} } );
  proto->placeConversationCacheBreakpoints( msgs1 );
  proto->appendEphemeralUserNote( msgs1, state_note("gather_peaks") );
  const json req1 = proto->buildRequestBody( api, sysPrompt, msgs1, tools, ToolChoice::Auto );

  // --- Round 2: the model called a tool, we fed the result back, state advanced. ---------------
  json msgs2 = json::array();
  msgs2.push_back( json{ {"role","user"}, {"content",initialTurn} } );
  msgs2.push_back( json{ {"role","assistant"}, {"content", json::array({
      json{ {"type","tool_use"}, {"id","toolu_1"}, {"name","get_peaks"}, {"input",json::object()} } }) } } );
  msgs2.push_back( json{ {"role","user"}, {"content", json::array({
      json{ {"type","tool_result"}, {"tool_use_id","toolu_1"}, {"content","1460.8 keV, 2614.5 keV"} } }) } } );
  proto->placeConversationCacheBreakpoints( msgs2 );
  proto->appendEphemeralUserNote( msgs2, state_note("fit_calibration") );  // different state text
  const json req2 = proto->buildRequestBody( api, sysPrompt, msgs2, tools, ToolChoice::Auto );

  const size_t prefixLen = cached_prefix_len( req1["messages"] );
  BOOST_REQUIRE_MESSAGE( prefixLen > 0, "round 1 wrote no cache breakpoint" );

  const vector<string> blocks1 = flatten_blocks( req1["messages"] );
  const vector<string> blocks2 = flatten_blocks( req2["messages"] );

  BOOST_REQUIRE_GE( blocks1.size(), prefixLen );
  BOOST_REQUIRE_MESSAGE( blocks2.size() >= prefixLen,
                         "round 2 is shorter than round 1's cached prefix" );

  // THE invariant: everything round 1 asked the provider to cache must reappear, byte-identical
  // and in the same position, in round 2 - otherwise the entry is unreadable and we paid the
  // cache-write premium for nothing.
  for( size_t i = 0; i < prefixLen; ++i )
  {
    BOOST_CHECK_MESSAGE( blocks1[i] == blocks2[i],
        "cached prefix diverges at block " + std::to_string(i)
        + "\n  round 1: " + blocks1[i] + "\n  round 2: " + blocks2[i] );
  }

  BOOST_CHECK( !marker_on_state_note( req2["messages"] ) );
}


BOOST_AUTO_TEST_CASE( AnthropicCacheBreakpointCountWithinApiLimit )
{
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  json messages = json::array();
  for( int i = 0; i < 6; ++i )
  {
    messages.push_back( json{ {"role","user"}, {"content", json::array({
        json{ {"type","tool_result"}, {"tool_use_id","toolu_" + std::to_string(i)}, {"content","ok"} } }) } } );
    messages.push_back( json{ {"role","assistant"}, {"content", json::array({
        json{ {"type","tool_use"}, {"id","toolu_" + std::to_string(i+1)}, {"name","get_peaks"}, {"input",json::object()} } }) } } );
  }
  messages.push_back( json{ {"role","user"}, {"content","and now summarize"} } );

  proto->placeConversationCacheBreakpoints( messages );
  proto->appendEphemeralUserNote( messages, state_note("summarize") );

  const json req = proto->buildRequestBody( make_anthropic_api(), "sys", messages, make_tools(), ToolChoice::Auto );

  // Anthropic allows at most 4 breakpoints per request; the `system` block already uses one.
  BOOST_CHECK_LE( count_markers( req["messages"] ), 3u );
  BOOST_CHECK_GT( count_markers( req["messages"] ), 0u );
  BOOST_CHECK( !marker_on_state_note( req["messages"] ) );
}


BOOST_AUTO_TEST_CASE( AnthropicCacheBreakpointWithoutEphemeralNote )
{
  // The MainAgent has no state machine, so appendEphemeralUserNote() is never called for it.
  // buildRequestBody() must still place the conversation breakpoints itself.
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  json messages = json::array();
  messages.push_back( json{ {"role","user"}, {"content","Fit the 1460 keV peak"} } );
  messages.push_back( json{ {"role","assistant"}, {"content", json::array({
      json{ {"type","tool_use"}, {"id","toolu_1"}, {"name","get_peaks"}, {"input",json::object()} } }) } } );
  messages.push_back( json{ {"role","user"}, {"content", json::array({
      json{ {"type","tool_result"}, {"tool_use_id","toolu_1"}, {"content","fit ok"} } }) } } );

  proto->placeConversationCacheBreakpoints( messages );
  const json req = proto->buildRequestBody( make_anthropic_api(), "sys", messages, make_tools(), ToolChoice::Auto );

  BOOST_CHECK_GT( count_markers( req["messages"] ), 0u );
  BOOST_CHECK_LE( count_markers( req["messages"] ), 3u );

  // The newest user turn must carry a marker - that is the entry the next round reads back.
  const json &lastMsg = req["messages"].back();
  BOOST_REQUIRE( lastMsg["content"].is_array() );
  BOOST_CHECK( lastMsg["content"].back().contains("cache_control") );
}


BOOST_AUTO_TEST_CASE( AnthropicSystemBlockCarriesBreakpoint )
{
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  json messages = json::array();
  messages.push_back( json{ {"role","user"}, {"content","hello"} } );

  const json req = proto->buildRequestBody( make_anthropic_api(), "a stable system prompt",
                                            messages, make_tools(), ToolChoice::Auto );

  // Tools render before system, so this one marker caches the whole tools+system prefix.
  BOOST_REQUIRE( req.contains("system") );
  BOOST_REQUIRE( req["system"].is_array() );
  BOOST_REQUIRE( !req["system"].empty() );
  BOOST_CHECK( req["system"].back().contains("cache_control") );
  BOOST_CHECK_EQUAL( req["system"].back()["cache_control"].value("type",""), string("ephemeral") );
}


BOOST_AUTO_TEST_CASE( AnthropicUserTurnWireShapeStableAcrossRounds )
{
  // Marking is what promotes a user turn's string content to a block array.  If only the marked
  // turns were promoted, a message's wire shape would oscillate between rounds (array while it is
  // the head/anchor, bare string once it is neither), and whether the provider treats the two
  // shapes as the same cached content would become a silent, cost-only assumption.
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  const string opening = "Context:\nTh232 spectrum\n\nTask:\nCalibrate energy.";

  json round1 = json::array();
  round1.push_back( json{ {"role","user"}, {"content",opening} } );
  proto->placeConversationCacheBreakpoints( round1 );

  // Three user turns later the opening turn is neither head nor anchor.
  json round3 = json::array();
  round3.push_back( json{ {"role","user"}, {"content",opening} } );
  for( int i = 1; i <= 2; ++i )
  {
    round3.push_back( json{ {"role","assistant"}, {"content", json::array({
        json{ {"type","tool_use"}, {"id","toolu_" + std::to_string(i)}, {"name","get_peaks"}, {"input",json::object()} } }) } } );
    round3.push_back( json{ {"role","user"}, {"content", json::array({
        json{ {"type","tool_result"}, {"tool_use_id","toolu_" + std::to_string(i)}, {"content","ok"} } }) } } );
  }
  proto->placeConversationCacheBreakpoints( round3 );

  json shape1 = round1[0]["content"];
  json shape3 = round3[0]["content"];
  for( json &b : shape1 ) if( b.is_object() ) b.erase( "cache_control" );
  for( json &b : shape3 ) if( b.is_object() ) b.erase( "cache_control" );

  BOOST_CHECK_MESSAGE( shape1 == shape3,
      "opening user turn changed wire shape between rounds\n  round 1: " + shape1.dump()
      + "\n  round 3: " + shape3.dump() );
}


BOOST_AUTO_TEST_CASE( AnthropicNoBreakpointOnNoteWhenNothingElseIsMarkable )
{
  // A ToolResult turn carrying zero tool calls serializes to an EMPTY content array, which is not
  // markable.  Placement must then mark nothing at all rather than falling through and marking the
  // ephemeral note - the note is the one block guaranteed not to reappear next round.
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  json messages = json::array();
  messages.push_back( json{ {"role","user"}, {"content", json::array()} } );  // empty tool-result turn

  proto->placeConversationCacheBreakpoints( messages );
  proto->appendEphemeralUserNote( messages, state_note("gather_peaks") );

  const json req = proto->buildRequestBody( make_anthropic_api(), "sys", messages, make_tools(), ToolChoice::Auto );

  BOOST_CHECK_MESSAGE( !marker_on_state_note( req["messages"] ),
                       "cache_control fell through onto the ephemeral reminder" );
  BOOST_CHECK_EQUAL( count_markers( req["messages"] ), 0u );
}


BOOST_AUTO_TEST_CASE( AnthropicOneShotRequestGetsNoMessageBreakpoint )
{
  // The compaction request is one-shot: its prompt differs every time, so a message-level marker
  // could never be read back and would be a pure cache-write surcharge.  buildRequestBody() must
  // not place one on its own - only the caller knows a conversation will be continued.
  const unique_ptr<const LlmApiProtocol> proto = LlmApiProtocol::create( ApiFormat::Anthropic );
  BOOST_REQUIRE( !!proto );

  json messages = json::array();
  messages.push_back( json{ {"role","user"}, {"content","Summarize the following conversation: ..."} } );

  const std::vector<NormalizedTool> noTools;
  const json req = proto->buildRequestBody( make_anthropic_api(), "compaction system prompt",
                                            messages, noTools, ToolChoice::Auto );

  BOOST_CHECK_EQUAL( count_markers( req["messages"] ), 0u );

  // ...but the system prompt is reused across compactions, so it keeps its breakpoint.
  BOOST_REQUIRE( req["system"].is_array() && !req["system"].empty() );
  BOOST_CHECK( req["system"].back().contains("cache_control") );
}
