#ifndef LlmRelActManualTool_h
#define LlmRelActManualTool_h
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

/** Namespace for LLM peak-based relative efficiency tool executors.

 This namespace contains tool executor functions for LLM-based peak relative efficiency analysis,
 including performing relative efficiency calculations using existing fit peaks and retrieving
 the current RelActManual state.
 */
namespace RelActManualTool
{
  /** Perform peak-based relative efficiency calculation.

   Performs a relative efficiency analysis using existing fit peaks. Peaks can be specified
   either by energy list or by source nuclide list. This is useful for:
   - Validating nuclide ID (peaks way off the fit curve are likely not the assigned nuclide)
   - Performing enrichment checks for Uranium or Plutonium
   - Calculating relative activities of fission products when geometry is unknown

   @param params JSON object with the following fields:
     - peak_energies: Array of peak energies in keV (mutually exclusive with sources)
     - sources: Array of source names like "U235", "U238" (mutually exclusive with peak_energies)
     - eqn_form: Optional rel eff equation form ("LnX", "LnY", "LnXLnY", "FramEmpirical", 
                 "FramPhysicalModel"), default "LnY"
     - eqn_order: Optional polynomial order (default 3, ignored for FramPhysicalModel)
     - match_tolerance: Optional peak matching tolerance in sigma (default 1.5)
     - background_subtract: Optional boolean for background subtraction (default false)
     - additional_uncertainty: Optional enum string ("Unweighted", "StatOnly", "OnePercent", 
                              "FivePercent", "TenPercent", "TwentyFivePercent", "FiftyPercent",
                              "SeventyFivePercent", "OneHundredPercent"), default "StatOnly"
     - nuclide_ages: Optional object mapping nuclide names to age strings, e.g., 
                    {"U235": "20 years", "U238": "4.5e9 years"}
     - save_to_state: Optional boolean to save results to app state (default false)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with calculation results including:
     - sources: Array of sources with rel_activity, rel_activity_uncert, mass_fraction
     - activity_ratios: Array of activity ratios between source pairs
     - mass_ratios: Array of mass ratios between source pairs
     - peaks: Array of peak info with energy, counts, observed_efficiency, fit_efficiency,
              residual_sigma (0-5 good, 5-10 acceptable, >10 poor)
     - fit_quality: Object with chi2, dof, chi2_per_dof
     - rel_eff_equation: String representation of the fit equation
     - warnings: Array of warning messages
     - saved_to_state: Boolean indicating if state was saved
   */
  nlohmann::json executePeakBasedRelativeEfficiency(const nlohmann::json& params, InterSpec* interspec);

  /** Get the current RelActManual state.

   Retrieves the current RelActManual configuration state from either the GUI (if open)
   or from the SpecMeas. This allows inspection of the current settings without
   performing a calculation.

   @param params JSON object (currently no parameters required)
   @param interspec Pointer to the InterSpec instance
   @returns JSON object with:
     - success: Boolean indicating success
     - has_state: Boolean indicating if a state exists
     - state: Object with current configuration (eqn_form, eqn_order, match_tolerance,
              additional_uncertainty, background_subtract, nuclide_ages, nuclide_decay_correct)
     - source_from_gui: Boolean indicating if state came from open GUI
   */
  nlohmann::json executeGetRelActManualState(const nlohmann::json& params, InterSpec* interspec);

}//namespace RelActManualTool

}//namespace LlmTools

#endif //LlmRelActManualTool_h

