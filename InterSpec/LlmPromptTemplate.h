#ifndef LLM_PROMPT_TEMPLATE_H
#define LLM_PROMPT_TEMPLATE_H
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

#include <nlohmann/json.hpp>

#include "InterSpec/LlmConfig.h"

/** Renders the LLM assistant's instruction text - agent system prompts, tool descriptions, per-state
 state-machine guidance, and the compaction prompt - as Inja templates, so the text can adapt to the
 currently-active model's capabilities and an author-selected verbosity level.

 Rendering is done once at model-setup time (the tool-registry build and
 `LlmInterface::prepareModelInstructions()`), not per request, and the finished strings are stored and
 reused by every turn.

 Templating is opt-in: any string containing none of the Inja markers "{{", "{%", "{#" is returned
 unchanged (see `render()`), so existing marker-free instruction text is byte-for-byte identical and
 never invokes Inja.
 */
namespace LlmPromptTemplate
{
  /** Which consumer the rendered text is destined for.

   - Agent: InterSpec's own LLM assistant, where the active model's real capabilities are known.
   - Mcp:   external MCP clients, whose model is remote/unknown - so `buildContext()` uses neutral
            capability defaults and sets `surface` = "mcp" (author templates should branch on that,
            not on `model.*`).
   */
  enum class Surface
  {
    Agent,
    Mcp
  };//enum class Surface

  /** Build the JSON "capability context" that instruction templates render against.

   Describes the active model (name, image support, reasoning, api-format, token limits), the
   author-selected `verbosity`, the `surface`, and a couple of feature flags (`mcp_enabled`,
   `deep_research_enabled`).  Pure; performs no I/O.  For `Surface::Mcp` the `agent` argument is
   ignored and `model.*` are neutral defaults.
   */
  nlohmann::json buildContext( const LlmConfig &config, Surface surface, AgentType agent );

  /** Render an Inja template string against `ctx`.

   Fast path: if `tmplt` contains none of "{{", "{%", "{#", it is returned unchanged.  On any Inja
   parse/render error it logs and returns the raw `tmplt` (these files are user-overridable, so a bad
   template must never break the assistant).  Never throws.
   */
  std::string render( const std::string &tmplt, const nlohmann::json &ctx );
}//namespace LlmPromptTemplate

#endif //LLM_PROMPT_TEMPLATE_H
