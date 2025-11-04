#ifndef LLM_TOOL_REGISTRY_H
#define LLM_TOOL_REGISTRY_H
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

#include <map>
#include <string>
#include <memory>

// Forward declarations
class InterSpec;
class DetectorPeakResponse;
class LlmConfig;

#include "InterSpec/LlmConfig.h"  // For AgentType enum
#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

namespace LlmTools {

/** Shared tool definition for the LLM tool registry.

 This is separate from the MCP Tool struct to allow for InterSpec* parameter.
 */
struct SharedTool {
    std::string name;
    std::string description;  // Default description
    nlohmann::json parameters_schema;
    // The executor takes parameters and InterSpec instance, returns result as JSON
    std::function<nlohmann::json(const nlohmann::json&, InterSpec*)> executor;

    // Agent-specific configurations
    std::vector<AgentType> availableForAgents;  // List of agent types that can use this tool (empty = all agents)
    std::map<AgentType, std::string> roleDescriptions;  // Role-specific descriptions (AgentType -> description)
};

/** Central registry for LLM tools that can be shared between LlmInterface and LlmMcpResource.
 
 We could implement a singleton pattern to manage tool registration and lookup, but currently arent to
 avoid any potential threading issues, or constraints as we implement further.
 */
class ToolRegistry {
public:
  ToolRegistry( const LlmConfig &config );
  ~ToolRegistry() = default;
  
  /** Register a new tool with the registry.
   @param tool The tool to register. If a tool with the same name exists, it will be replaced.
   */
  void registerTool(const SharedTool& tool);
  
  /** Register all default tools provided by InterSpec. */
  void registerDefaultTools( const LlmConfig &config );
  
  /** Get all registered tools.
   @return A map of tool name to SharedTool struct.
   */
  const std::map<std::string, SharedTool>& getTools() const;

  /** Get tools available for a specific agent.
   @param agentType The type of the agent
   @return Map of tool name to SharedTool struct, filtered for the agent
   */
  std::map<std::string, SharedTool> getToolsForAgent( AgentType agentType ) const;

  /** Get the description for a tool for a specific agent role.
   @param toolName The name of the tool
   @param agentType The type of the agent/role
   @return The role-specific description if available, otherwise the default description
   */
  std::string getDescriptionForAgent( const std::string &toolName, AgentType agentType ) const;
  
  /** Get a specific tool by name.
   @param name The name of the tool to look up.
   @return Pointer to the tool if found, nullptr otherwise.
   */
  const SharedTool* getTool(const std::string& name) const;
  
  /** Execute a tool by name with the given parameters and InterSpec context.
   @param toolName The name of the tool to execute.
   @param parameters JSON parameters for the tool.
   @param interspec The InterSpec instance to execute the tool against.
   @return JSON result from the tool execution.
   @throws std::runtime_error if tool not found or execution fails.
   */
  nlohmann::json executeTool(const std::string& toolName, 
                           const nlohmann::json& parameters, 
                           InterSpec* interspec) const;
  
  /** Clear all registered tools (mainly for testing). */
  void clearTools();
  
private:
  ToolRegistry() = delete;
  
  // Prevent copying
  ToolRegistry(const ToolRegistry&) = delete;
  ToolRegistry& operator=(const ToolRegistry&) = delete;
  
  std::map<std::string, SharedTool> m_tools;

  /** Factory function to create a SharedTool with executor based on tool name.
   *  Returns a tool with only the name and executor set - description and schema
   *  will be populated from the ToolConfig.
   */
  static SharedTool createToolWithExecutor( const std::string &toolName );

  // Helper functions for default tools
  //  The conventions that seem useful are:
  //  - If there is an error, throw an exception and the LLM executor will take care of catching and sending a response
  //  - Use "source", and not "nuclide" everywhere - to be consistent and inclusive of flourescent x-rays and reactions,
  //    unless it really is a nuclide-only thing, like decay chain.
  static nlohmann::json executePeakDetection(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetSpectrumInfo(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executePeakFit(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetUserPeaks(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetCharacteristicGammasForSource( const nlohmann::json& params );
  static nlohmann::json executeGetNuclidesWithCharacteristicsInEnergyRange( const nlohmann::json& params, InterSpec* interspec );
  static nlohmann::json executeGetLoadedSpectra(const nlohmann::json& params, InterSpec* interspec);
// I read that you should minimize the number of tool calls - so this next define makes it so the
// "analyst_notes_for_source" and "sources_associated_with_source" tool calls will not be defined, but rather this info
// is included with the "source_info" tool call
#define INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO 1
  static nlohmann::json executeGetAssociatedSources(const nlohmann::json& params );
  static nlohmann::json executeGetSourceAnalystNotes(const nlohmann::json& params );
  static nlohmann::json executeGetSourceInfo(const nlohmann::json& params, InterSpec* interspec );
  static nlohmann::json executeGetNuclideDecayChain(const nlohmann::json& params );
  static nlohmann::json executeGetAutomatedRiidId(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeFitPeaksForNuclide(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetCountsInEnergyRange(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetExpectedFwhm(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeCurrieMdaCalc(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetMaterials(InterSpec* interspec);
  static nlohmann::json executeGetMaterialInfo(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetAttenuationOfShielding(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetSourcePhotons(const nlohmann::json& params);
  
  static nlohmann::json executeAvailableDetectors(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeLoadDetectorEfficiency(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeGetDetectorInfo(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executePhotopeakDetectionCalc(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeSearchSourcesByEnergy(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeEditAnalysisPeak(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeEscapePeakCheck(const nlohmann::json& params, InterSpec* interspec);
  static nlohmann::json executeSumPeakCheck(const nlohmann::json& params, InterSpec* interspec);

  /** Helper function to find and load a detector by identifier.
   @param identifier The detector identifier (name, path, or URI)
   @param detectorName Optional specific detector name for multi-detector files
   @param sourceHint Source hint for where to search
   @param interspec InterSpec instance (may be null for some sources)
   @param loadedFrom Output parameter indicating where detector was loaded from
   @returns Shared pointer to loaded detector, or nullptr if not found
   */
  static std::shared_ptr<DetectorPeakResponse> findDetectorByIdentifier(
    const std::string& identifier,
    const std::string& detectorName,
    const std::string& sourceHint,
    InterSpec* interspec,
    std::string& loadedFrom
  );

};

} // namespace LlmTools

#endif // LLM_TOOL_REGISTRY_H
