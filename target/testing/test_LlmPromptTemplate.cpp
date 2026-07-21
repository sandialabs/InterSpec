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
