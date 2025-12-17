#ifndef LlmIsotopicsTool_h
#define LlmIsotopicsTool_h
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

/** Namespace for LLM isotopics analysis tool executors.

 This namespace contains all the tool executor functions for LLM-based isotopics analysis,
 including configuration management, nuclide manipulation, ROI management, and calculation execution.
 */
namespace IsotopicsTool
{
  /** List available isotopics preset configurations.

   Returns a list of available preset configurations that can be loaded via load_isotopics_preset.
   Presets are found in the user and static data directories under "auto_rel_eff_eqn_form/".

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and "presets" array
   */
  nlohmann::json executeListIsotopicsPresets(const nlohmann::json& params, InterSpec* interspec);

  /** Get the current isotopics configuration.

   Returns the complete current isotopics configuration including nuclides, ROIs,
   relative efficiency curve settings, and global options.

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with complete configuration state
   */
  nlohmann::json executeGetIsotopicsConfig(const nlohmann::json& params, InterSpec* interspec);

  /** Reset the isotopics configuration to defaults.

   Clears all nuclides, ROIs, and resets settings to default values.

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" field
   */
  nlohmann::json executeResetIsotopicsConfig(const nlohmann::json& params, InterSpec* interspec);

  /** Load an isotopics preset configuration.

   Loads a named preset configuration from the available presets.
   Use executeListIsotopicsPresets to see available preset names.

   @param params JSON object with "preset" string field
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and loaded configuration
   */
  nlohmann::json executeLoadIsotopicsPreset(const nlohmann::json& params, InterSpec* interspec);

  /** Perform the isotopics calculation.

   Executes the isotopics analysis using the current configuration and spectrum.
   Returns results including nuclide activities, chi-square, and fit quality metrics.

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with calculation results
   */
  nlohmann::json executePerformIsotopics(const nlohmann::json& params, InterSpec* interspec);

  /** Modify nuclides in the isotopics configuration.

   Add, remove, or update nuclides in the current configuration.
   Supports actions: "add", "remove", "update"

   @param params JSON object with "action" and "nuclides" array
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and updated configuration summary
   */
  nlohmann::json executeModifyIsotopicsNuclides(const nlohmann::json& params, InterSpec* interspec);

  /** Modify ROIs (regions of interest) in the isotopics configuration.

   Add, remove, update, or replace all ROIs in the current configuration.
   Supports actions: "add", "remove", "update", "replace_all"

   @param params JSON object with "action" and "rois" array
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and updated configuration summary
   */
  nlohmann::json executeModifyIsotopicsRois(const nlohmann::json& params, InterSpec* interspec);

  /** Modify relative efficiency curve settings.

   Change the equation type, order, shielding, and other curve-specific settings.

   @param params JSON object with curve settings fields to modify
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and updated configuration summary
   */
  nlohmann::json executeModifyIsotopicsCurveSettings(const nlohmann::json& params, InterSpec* interspec);

  /** Modify global isotopics analysis options.

   Change FWHM form, skew type, energy calibration fitting, background subtraction,
   and other global analysis settings.

   @param params JSON object with option fields to modify
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and updated configuration summary
   */
  nlohmann::json executeModifyIsotopicsOptions(const nlohmann::json& params, InterSpec* interspec);

  /** Modify activity ratio and mass fraction constraints.

   Add, remove, update, or clear all constraints for nuclides.
   Supports constraint types: "activity_ratio", "mass_fraction"
   Supports actions: "add", "remove", "update", "clear_all"

   @param params JSON object with "action", "constraint_type", and "constraints" array
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with "success" and updated configuration summary
   */
  nlohmann::json executeModifyIsotopicsConstraints(const nlohmann::json& params, InterSpec* interspec);

  /** Get schema and valid values for isotopics configuration fields.

   Returns detailed information about all configurable fields including valid enum values,
   field descriptions, constraints, and relationships between fields.

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with schema information
   */
  nlohmann::json executeGetIsotopicsConfigSchema(const nlohmann::json& params, InterSpec* interspec);

}//namespace IsotopicsTool

}//namespace LlmTools

#endif //LlmIsotopicsTool_h
