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
#include <iostream>

#define BOOST_TEST_MODULE LlmPromptTemplate_suite
#include <boost/test/included/unit_test.hpp>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmPromptTemplate.h"

using namespace std;
using json = nlohmann::json;

using Verbosity = LlmConfig::LlmApi::InstructionVerbosity;
using Surface = LlmPromptTemplate::Surface;

namespace
{
  // Build a minimal, valid in-memory LlmConfig with a single active provider + model.
  LlmConfig make_config( const bool supportsImages,
                         const Verbosity verbosity,
                         const string &endpoint = "https://api.anthropic.com/v1/messages" )
  {
    LlmConfig cfg;
    cfg.llmApi.enabled = true;

    LlmConfig::LlmApi::ModelInfo m;
    m.name = "test-model";
    m.supportsImages = supportsImages;
    m.reasoning = LlmConfig::LlmApi::ReasoningEffort::high;
    m.maxTokens = 1234;
    m.contextLengthLimit = 5678;
    m.instructionVerbosity = verbosity;

    LlmConfig::LlmApi::ApiProvider prov;
    prov.apiEndpoint = endpoint;
    prov.models.push_back( m );
    prov.activeModelIndex = 0;

    cfg.llmApi.providers.push_back( prov );
    cfg.llmApi.activeProviderIndex = 0;
    cfg.mcpServer.enabled = false;
    return cfg;
  }//make_config(...)


  /** Locate the shipped `data` directory, from `--datadir=` or by searching upward. */
  string data_dir()
  {
    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char ** const argv = boost::unit_test::framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
    }
    SpecUtils::ireplace_all( datadir, "%20", " " );

    if( datadir.empty() )
    {
      for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, "llm_agent_NuclideId.xml" ) ) )
        {
          datadir = d;
          break;
        }
      }
    }//if( datadir.empty() )

    return datadir;
  }//data_dir()
}//namespace


// Marker-free text (prose, markdown, single-brace JSON) must pass through byte-for-byte, so that
// existing (non-templated) instruction text is unchanged.
BOOST_AUTO_TEST_CASE( FastPathLeavesMarkerFreeTextUnchanged )
{
  const json ctx = LlmPromptTemplate::buildContext(
      make_config( true, Verbosity::Normal ), Surface::Agent, AgentType::MainAgent );

  const vector<string> unchanged = {
    "",
    "Identify peaks, then assign nuclides.",
    "## Heading\n\n1. Do X\n2. Do Y\n- bullet with { brace } and $var",
    "Example schema: { \"type\": \"object\", \"properties\": { \"AD\": 20.25 } }",
    "Use `set_workflow_state` to transition."
  };

  for( const string &s : unchanged )
    BOOST_CHECK_EQUAL( LlmPromptTemplate::render( s, ctx ), s );
}


BOOST_AUTO_TEST_CASE( RendersModelVariableAndConditionals )
{
  const json imgCtx = LlmPromptTemplate::buildContext(
      make_config( true, Verbosity::Normal ), Surface::Agent, AgentType::MainAgent );
  const json noImgCtx = LlmPromptTemplate::buildContext(
      make_config( false, Verbosity::Normal ), Surface::Agent, AgentType::MainAgent );

  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( "{{ model.supports_images }}", imgCtx ), "true" );
  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( "{{ model.supports_images }}", noImgCtx ), "false" );

  // Conditional inclusion (no whitespace adjacent to the tags, so it is unambiguous).
  const string tmplt = "X{% if model.supports_images %}IMG{% endif %}Y";
  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( tmplt, imgCtx ), "XIMGY" );
  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( tmplt, noImgCtx ), "XY" );

  // A realistic instruction appears only when the model supports images.
  const string real = "{% if model.supports_images %}You can view the spectrum image.{% endif %}";
  BOOST_CHECK( LlmPromptTemplate::render( real, imgCtx ).find( "spectrum image" ) != string::npos );
  BOOST_CHECK( LlmPromptTemplate::render( real, noImgCtx ).find( "spectrum image" ) == string::npos );
}


// The shipped NuclideId VALIDATE_CANDIDATE guidance carries a vision-only instruction to visually
// verify sub-50 keV peaks (they are easily faked by the detector turn-on edge, or by Ge/NaI x-ray
// escape peaks).  Render the real file so an unbalanced {% if %} or a stray raw '<' is caught here.
BOOST_AUTO_TEST_CASE( NuclideIdLowEnergyGuidanceIsVisionGated )
{
  const string datadir = data_dir();
  const string agent_file = SpecUtils::append_path( datadir, "llm_agent_NuclideId.xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( agent_file ),
                         "llm_agent_NuclideId.xml not at '" << agent_file << "'" );

  vector<LlmConfig::AgentConfig> agents;
  BOOST_REQUIRE_NO_THROW( agents = LlmConfig::loadAgentsFromFile( agent_file ) );
  BOOST_REQUIRE( !agents.empty() );
  BOOST_REQUIRE( agents[0].state_machine );
  BOOST_REQUIRE( agents[0].state_machine->hasState( "VALIDATE_CANDIDATE" ) );

  const string guidance
      = agents[0].state_machine->getStateDefinition( "VALIDATE_CANDIDATE" ).prompt_guidance;
  BOOST_REQUIRE( guidance.find( "{% if model.supports_images %}" ) != string::npos );

  const json imgCtx = LlmPromptTemplate::buildContext(
      make_config( true, Verbosity::Normal ), Surface::Agent, AgentType::NuclideId );
  const json noImgCtx = LlmPromptTemplate::buildContext(
      make_config( false, Verbosity::Normal ), Surface::Agent, AgentType::NuclideId );

  const string withImages = LlmPromptTemplate::render( guidance, imgCtx );
  const string withoutImages = LlmPromptTemplate::render( guidance, noImgCtx );

  BOOST_CHECK( withImages.find( "Peaks below ~50 keV" ) != string::npos );
  BOOST_CHECK( withImages.find( "get_spectrum_image" ) != string::npos );
  BOOST_CHECK( withoutImages.find( "Peaks below ~50 keV" ) == string::npos );
  BOOST_CHECK( withoutImages.find( "get_spectrum_image" ) == string::npos );

  // Neither rendering may leave template markers behind, and the guidance that is not gated
  // (e.g. the Ultimate Parent rule) must survive in both.
  BOOST_CHECK( withImages.find( "{%" ) == string::npos );
  BOOST_CHECK( withoutImages.find( "{%" ) == string::npos );
  BOOST_CHECK( withImages.find( "Ultimate Parent Rule" ) != string::npos );
  BOOST_CHECK( withoutImages.find( "Ultimate Parent Rule" ) != string::npos );
}


BOOST_AUTO_TEST_CASE( VerbosityGating )
{
  const json verbose = LlmPromptTemplate::buildContext(
      make_config( true, Verbosity::Verbose ), Surface::Agent, AgentType::NuclideId );
  const json terse = LlmPromptTemplate::buildContext(
      make_config( true, Verbosity::Terse ), Surface::Agent, AgentType::NuclideId );

  const string tmplt = "{% if verbosity == \"verbose\" %}LONG{% else %}SHORT{% endif %}";
  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( tmplt, verbose ), "LONG" );
  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( tmplt, terse ), "SHORT" );
}


// A malformed template must not throw or abort; it degrades to the raw (unrendered) text.
BOOST_AUTO_TEST_CASE( MalformedTemplateReturnsRawAndDoesNotThrow )
{
  const json ctx = LlmPromptTemplate::buildContext(
      make_config( true, Verbosity::Normal ), Surface::Agent, AgentType::MainAgent );

  const string bad = "{{ unclosed";
  string out;
  BOOST_CHECK_NO_THROW( out = LlmPromptTemplate::render( bad, ctx ) );
  BOOST_CHECK_EQUAL( out, bad );
}


BOOST_AUTO_TEST_CASE( BuildContextFieldValues )
{
  const LlmConfig cfg = make_config( true, Verbosity::Verbose );
  const json ctx = LlmPromptTemplate::buildContext( cfg, Surface::Agent, AgentType::NuclideId );

  BOOST_CHECK_EQUAL( ctx["model"]["name"].get<string>(), "test-model" );
  BOOST_CHECK_EQUAL( ctx["model"]["supports_images"].get<bool>(), true );
  BOOST_CHECK_EQUAL( ctx["model"]["reasoning_enabled"].get<bool>(), true );
  BOOST_CHECK_EQUAL( ctx["model"]["reasoning_effort"].get<string>(), "high" );
  BOOST_CHECK_EQUAL( ctx["model"]["api_format"].get<string>(), "anthropic" );
  BOOST_CHECK_EQUAL( ctx["model"]["max_tokens"].get<int>(), 1234 );
  BOOST_CHECK_EQUAL( ctx["model"]["context_length_limit"].get<int>(), 5678 );
  BOOST_CHECK_EQUAL( ctx["verbosity"].get<string>(), "verbose" );
  BOOST_CHECK_EQUAL( ctx["surface"].get<string>(), "agent" );
  BOOST_CHECK_EQUAL( ctx["agent"]["name"].get<string>(), "NuclideId" );
}


// The MCP surface is model-agnostic: even with a text-only, terse local model it exposes neutral
// capabilities and surface=="mcp", so authors branch on `surface` rather than local model caps.
BOOST_AUTO_TEST_CASE( McpSurfaceIsModelAgnostic )
{
  const LlmConfig cfg = make_config( false, Verbosity::Terse );
  const json ctx = LlmPromptTemplate::buildContext( cfg, Surface::Mcp, AgentType::MainAgent );

  BOOST_CHECK_EQUAL( ctx["surface"].get<string>(), "mcp" );
  BOOST_CHECK_EQUAL( ctx["model"]["supports_images"].get<bool>(), true );   // neutral default
  BOOST_CHECK_EQUAL( ctx["verbosity"].get<string>(), "normal" );            // neutral default
  BOOST_CHECK( !ctx.contains( "agent" ) );                                  // agent omitted for MCP

  const string tmplt = "{% if surface == \"mcp\" %}MCP{% else %}AGENT{% endif %}";
  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( tmplt, ctx ), "MCP" );
}


// Regression: Inja's default line-statement token is "##", which collides with the Markdown
// "##"/"###" headers used throughout the agent prompts.  We disable line statements, so headers must
// render as literal text even when the prompt also contains a real {% %} tag.
BOOST_AUTO_TEST_CASE( MarkdownHeadersAreNotLineStatements )
{
  const json imgCtx = LlmPromptTemplate::buildContext(
      make_config( true, Verbosity::Normal ), Surface::Agent, AgentType::MainAgent );
  const json noImgCtx = LlmPromptTemplate::buildContext(
      make_config( false, Verbosity::Normal ), Surface::Agent, AgentType::MainAgent );

  const string tmplt =
      "## Role and Goal\n\nYou are an expert.\n### Steps\n1. Do X{% if model.supports_images %} (and view the image){% endif %}.\n";

  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( tmplt, imgCtx ),
      "## Role and Goal\n\nYou are an expert.\n### Steps\n1. Do X (and view the image).\n" );
  BOOST_CHECK_EQUAL( LlmPromptTemplate::render( tmplt, noImgCtx ),
      "## Role and Goal\n\nYou are an expert.\n### Steps\n1. Do X.\n" );
}
