#ifndef LLM_KNOWLEDGE_BASE_H
#define LLM_KNOWLEDGE_BASE_H
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

#include <functional>
#include <string>
#include <vector>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

namespace LlmTools
{
struct SharedTool;
}

namespace LlmKnowledgeBase
{

/** A supplemental file in a skill's `references/` sub-directory; loaded only on request, so a
 skill's SKILL.md can stay a short index while detailed material is retrieved piecemeal.
 */
struct SkillReference
{
  /** File stem (e.g., "identifinder-family" for "references/identifinder-family.md"). */
  std::string name;
  /** From the reference's frontmatter `description`, or else its first Markdown heading. */
  std::string summary;
  std::string file_path;
};//struct SkillReference


struct SkillInfo
{
  std::string canonical_name;
  std::string skill_name;
  std::string description;
  std::string body;
  std::string file_path;
  std::vector<SkillReference> references;
};//struct SkillInfo


std::string skillToolNameFromCanonical( std::string canonical_name );

/** Loads all `SKILL.md` files (and their `references/*.md`) found recursively under the given
 `llm_knowledge` directory.  Skills are sorted by path; duplicate canonical names keep the first.
 */
std::vector<SkillInfo> loadSkillsFromKnowledgeDir( const std::string &knowledge_dir );

/** The existing directories skills are loaded from, in increasing order of precedence:
 `<static data>/llm_knowledge`, `<static data>/gamma_spec_skills/llm_knowledge`,
 `<writable data>/llm_knowledge`, `<writable data>/gamma_spec_skills/llm_knowledge`.
 The `gamma_spec_skills/llm_knowledge` locations let a checked-out knowledge repository be used in place.
 */
std::vector<std::string> knowledgeDirectories();

/** Loads skills from all of #knowledgeDirectories; a skill in a later directory overrides a same-named
 skill (by canonical name) from an earlier one - e.g., a writable-directory skill overrides a shipped one.
 */
std::vector<SkillInfo> loadSkills();

/** Text appended to the KnowledgeBase agent's `invoke_KnowledgeBase` tool description, so the calling agent
 knows what topics it can answer: a comma-separated list of the available skill names (kept brief, since
 every calling agent sees this on every turn), a note that it
 can run JavaScript computations, and - if `has_rag` - that it can also search reference documents.
 */
std::string invokeDescriptionSuffix( const std::vector<SkillInfo> &skills, const bool has_rag );

/** The JSON parameter schema for a skill's tool; includes an optional `reference` enum if the skill
 has any references.
 */
nlohmann::json skillToolParameterSchema( const SkillInfo &skill );

/** Executes a skill tool call.  Without a `reference` parameter, returns the SKILL.md body, plus
 the list of available references; with one, returns that reference's content.
 Throws if an unknown reference is requested, or the reference file can not be read.
 */
nlohmann::json executeSkillTool( const SkillInfo &skill, const nlohmann::json &params );

/** Registers the KnowledgeBase skill tools, and - if `rag_url` is non-empty - the tool that
 queries that endpoint.  An empty `rag_model` or `rag_corpora` means the built-in
 default is sent for that field of the request.
 */
void registerKnowledgeBaseTools( const std::vector<SkillInfo> &skills,
                                 const std::string &rag_url,
                                 const std::string &rag_model,
                                 const std::vector<std::string> &rag_corpora,
                                 const std::function<void(const LlmTools::SharedTool&)> &register_tool );

nlohmann::json executeQueryRagEndpoint( const nlohmann::json& params,
                                                 const std::string &rag_url,
                                                 const std::string &rag_model,
                                                 const std::vector<std::string> &rag_corpora );

}//namespace LlmKnowledgeBase

#endif // LLM_KNOWLEDGE_BASE_H
