#ifndef LlmActivityFitTool_h
#define LlmActivityFitTool_h
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

// Forward declarations
class InterSpec;

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"


namespace LlmTools
{

/** Namespace for LLM activity/shielding fit tool executors.

 This namespace contains all the tool executor functions for LLM-based activity and shielding
 source fitting, including peak marking, configuration management, and fit execution.
 */
namespace ActivityFitTool
{
  /** Close the activity/shielding fit display.

   Closes the activity/shielding source fit GUI window if it's currently displayed.

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" field
   */
  nlohmann::json executeCloseActivityShieldingDisplay(const nlohmann::json& params, InterSpec* interspec);

  /** Ask the user a question and wait for response.

   Displays a dialog to the user with a question and waits for their response.
   This is useful when the LLM needs clarification or additional information.

   @param params JSON object with "question" string field
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and "answer" fields
   */
  nlohmann::json executeAskUserQuestion(const nlohmann::json& params, InterSpec* interspec);

  /** Mark peaks for activity/shielding fitting.

   Marks or unmarks peaks as available for activity/shielding source fitting.
   Peaks can be specified by energy or index.

   @param params JSON object with "action" ("mark" or "unmark") and "peaks" array
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and number of peaks marked/unmarked
   */
  nlohmann::json executeMarkPeaksForActivityFit(const nlohmann::json& params, InterSpec* interspec);

  /** Get the current shielding/source configuration.

   Returns the current activity/shielding fit configuration including:
   - Source definitions (nuclides, activities, ages, fit settings)
   - Shielding definitions (materials, dimensions, fit settings)
   - Geometry and distance settings
   - Fit options

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with complete configuration state
   */
  nlohmann::json executeGetShieldingSourceConfig(const nlohmann::json& params, InterSpec* interspec);

  /** Modify the shielding/source configuration.

   Modifies the activity/shielding fit configuration. Supports operations:
   - add_source: Add a new source
   - remove_source: Remove a source by nuclide name
   - modify_source: Update source parameters (activity, age, fit settings)
   - add_shielding: Add shielding material
   - remove_shielding: Remove shielding by index
   - modify_shielding: Update shielding parameters
   - set_distance: Change source-detector distance
   - set_geometry: Change source geometry
   - set_options: Modify global fit options

   @param params JSON object with "operation" and operation-specific fields
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and updated configuration summary
   */
  nlohmann::json executeModifyShieldingSourceConfig(const nlohmann::json& params, InterSpec* interspec);

  /** Perform activity/shielding fit calculation using app state.

   Executes an activity and shielding fit using the current GUI configuration.
   Requires the Activity/Shielding fit GUI to be open.

   Returns fit results including activities, ages, shielding, and fit quality metrics.

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with fit results
   */
  nlohmann::json executeActivityFit(const nlohmann::json& params, InterSpec* interspec);

  /** Perform one-off activity/shielding fit calculation with full specification.

   Executes an activity and shielding fit with a fully-specified configuration.
   Does not use or modify the GUI state. This is ideal for one-off fits where you
   want complete control over all parameters.

   Required parameters:
   - distance: Source-to-detector distance with units (e.g., "1 m", "100 cm")
   - geometry: Geometry type ("Spherical", "CylinderSideOn", "CylinderEndOn")
   - peak_energies: Array of peak energies to use for fitting (in keV)

   Optional parameters:
   - sources: Array of source definitions (if omitted, auto-generated from peaks)
     Each source can have: nuclide, activity, fit_activity, age, fit_age
   - shielding: Array of shielding layer definitions
     Each layer can have: material, thickness/radius/length, fit flags, trace_sources
   - detector: Detector efficiency function name
   - options: Fit options (background_peak_subtract, same_age_isotopes, etc.)

   @param params JSON object with full fit specification
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with fit results
   */
  nlohmann::json executeActivityFitOneOff(const nlohmann::json& params, InterSpec* interspec);

}//namespace ActivityFitTool

}//namespace LlmTools

#endif //LlmActivityFitTool_h
