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
#include <variant>
#include <iostream>
#include <exception>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmPromptTemplate.h"

using namespace std;

namespace
{
  const char *verbosityToString( const LlmConfig::LlmApi::InstructionVerbosity v )
  {
    switch( v )
    {
      case LlmConfig::LlmApi::InstructionVerbosity::Terse:   return "terse";
      case LlmConfig::LlmApi::InstructionVerbosity::Normal:  return "normal";
      case LlmConfig::LlmApi::InstructionVerbosity::Verbose: return "verbose";
    }
    return "normal";
  }//verbosityToString(...)


  const char *apiFormatToJsonString( const LlmConfig::LlmApi::ApiFormat fmt )
  {
    switch( fmt )
    {
      case LlmConfig::LlmApi::ApiFormat::OpenAiChat:      return "openai_chat";
      case LlmConfig::LlmApi::ApiFormat::OpenAiResponses: return "openai_responses";
      case LlmConfig::LlmApi::ApiFormat::Anthropic:       return "anthropic";
    }
    return "openai_chat";
  }//apiFormatToJsonString(...)
}//namespace


namespace LlmPromptTemplate
{

nlohmann::json buildContext( const LlmConfig &config, const Surface surface, const AgentType agent )
{
  const LlmConfig::LlmApi &api = config.llmApi;
  const bool isMcp = (surface == Surface::Mcp);

  // Flatten reasoning (variant<bool, ReasoningEffort>) into (enabled, effort-string) so templates can
  // branch on a single string.
  bool reasoning_enabled = false;
  string reasoning_effort = "off";
  const std::variant<bool, LlmConfig::LlmApi::ReasoningEffort> &reasoning = api.reasoning();
  if( std::holds_alternative<bool>( reasoning ) )
  {
    reasoning_enabled = std::get<bool>( reasoning );
    reasoning_effort = reasoning_enabled ? "on" : "off";
  }else
  {
    reasoning_enabled = true;
    switch( std::get<LlmConfig::LlmApi::ReasoningEffort>( reasoning ) )
    {
      case LlmConfig::LlmApi::ReasoningEffort::low:    reasoning_effort = "low";    break;
      case LlmConfig::LlmApi::ReasoningEffort::medium: reasoning_effort = "medium"; break;
      case LlmConfig::LlmApi::ReasoningEffort::high:   reasoning_effort = "high";   break;
    }
  }

  nlohmann::json model;
  model["name"] = api.model();
  // The MCP consuming model is remote/unknown: expose neutral defaults and always allow image tools
  // (the MCP server exposes them unconditionally); authors branch on `surface`, not these.
  model["supports_images"] = isMcp ? true : api.supportsImages();
  model["reasoning_enabled"] = isMcp ? false : reasoning_enabled;
  model["reasoning_effort"] = isMcp ? string("off") : reasoning_effort;
  model["api_format"] = apiFormatToJsonString( api.apiFormat() );
  model["max_tokens"] = api.maxTokens();
  model["context_length_limit"] = api.contextLengthLimit();
  if( api.temperature().has_value() )
    model["temperature"] = api.temperature().value();
  else
    model["temperature"] = nullptr;

  nlohmann::json ctx;
  ctx["model"] = std::move( model );
  ctx["verbosity"] = isMcp ? string("normal") : string( verbosityToString( api.instructionVerbosity() ) );
  ctx["surface"] = isMcp ? "mcp" : "agent";
  if( !isMcp )
    ctx["agent"] = { { "name", agentTypeToString( agent ) } };
  ctx["mcp_enabled"] = config.mcpServer.enabled;
  ctx["deep_research_enabled"] = !api.deep_research_url.empty();

  return ctx;
}//buildContext(...)


std::string render( const std::string &tmplt, const nlohmann::json &ctx )
{
  // Fast path: no Inja markers -> return unchanged.  This makes templating opt-in and guarantees all
  // existing marker-free instruction text is byte-for-byte identical and never touches Inja.
  if( (tmplt.find("{{") == string::npos)
      && (tmplt.find("{%") == string::npos)
      && (tmplt.find("{#") == string::npos) )
    return tmplt;

  try
  {
    inja::Environment env;
    // Instruction templates never pull other templates off disk; keep rendering hermetic.
    env.set_search_included_templates_in_files( false );
    // Disable Inja "line statements": its default line-statement token is "##", which collides with
    // the Markdown "##"/"###" headers used throughout the agent prompts.  Point it at a control-char
    // sentinel that never appears in prose, so only the explicit {{ }} / {% %} / {# #} forms are tags.
    env.set_line_statement( "\x01" );
    // Standard Jinja whitespace (no block trimming): wrapping a clause in {% if %}...{% endif %}
    // reproduces the original text byte-for-byte for the enabled branch.  Use {%- / -%} to trim.
    return env.render( tmplt, ctx );
  }catch( const std::exception &e )
  {
    // A malformed template must never break the assistant (these files are user-overridable), so log
    // and fall through to returning the raw text.  Broken *shipped* templates are caught by the unit
    // tests and the DebugFile request dump rather than by aborting here.
    cerr << "LlmPromptTemplate::render: failed to render instruction template; using raw text. Error: "
         << e.what() << endl;
  }//try / catch

  // Degrade gracefully to the raw (unrendered) template so a bad template never breaks the assistant.
  return tmplt;
}//render(...)

}//namespace LlmPromptTemplate
