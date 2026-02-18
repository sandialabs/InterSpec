#ifndef LLM_DEEP_RESEARCH_AGENT_H
#define LLM_DEEP_RESEARCH_AGENT_H
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

namespace LlmDeepResearch
{

struct SkillInfo
{
  std::string canonical_name;
  std::string skill_name;
  std::string description;
  std::string body;
  std::string file_path;
};//struct SkillInfo


std::string skillToolNameFromCanonical( std::string canonical_name );

std::vector<SkillInfo> loadSkills();

void registerDeepResearchTools( const std::string &deep_research_url,
                                const std::function<void(const LlmTools::SharedTool&)> &register_tool,
                                bool &loaded_any_skill );

nlohmann::json executeQueryDeepResearchEndpoint( const nlohmann::json& params, const std::string &deep_research_url );

}//namespace LlmDeepResearch

#endif // LLM_DEEP_RESEARCH_AGENT_H
