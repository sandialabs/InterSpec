#ifndef DecayBatchCalc_h
#define DecayBatchCalc_h
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

namespace SandiaDecay
{
  struct Nuclide;
}

/** Non-GUI core for the "Batch Decay" tool.

 Decays potentially many initial nuclides for one, or a number of, time points and produces a
 table of activities (and, optionally, particle-line rates).  This replaces the legacy command-line
 utility at external_libs/SandiaDecay/examples/batch_decay.cpp.

 The output table is oriented rows = items (nuclides or particle-lines), columns = time steps -
 matching the legacy tool.  The activity/particle units are carried in the row label so different
 row types (Bq/Ci vs counts/second) can coexist in one table.
 */
namespace DecayBatchCalc
{
  /** A single initial nuclide the user wants to decay. */
  struct BatchNuclide
  {
    /** Resolved nuclide; may be null if the input string could not be resolved (an error). */
    const SandiaDecay::Nuclide *nuclide = nullptr;

    /** The nuclide string as entered/resolved (e.g. "U238"). */
    std::string nuclide_str;

    /** Initial age of the nuclide, in SandiaDecay time units (seconds). */
    double age = 0.0;

    /** Initial activity, in SandiaDecay activity units (becquerel). */
    double activity = 0.0;

    /** Optional areal/volumetric suffix to carry through to output labels (e.g. "/m2").
     Empty for a normal activity input.  Purely a display label; the numeric decay is unaffected.
     */
    std::string unit_label;
  };//struct BatchNuclide


  /** Decay options; mirrors the legacy CLI arguments and the CSV-export dialog options. */
  struct BatchDecayOptions
  {
    /** Time to decay to, in SandiaDecay time units (seconds).  Must be > 0. */
    double time_span = 0.0;

    /** Number of time points to evaluate, from 0 to `time_span` (inclusive when > 1).
     A value of 1 evaluates only at `time_span`.  Must be >= 1.
     */
    std::size_t num_steps = 1;

    /** Sum all inputs into a single mixture and co-decay them (implies show_progeny). */
    bool mix_input = false;

    /** For each (un-mixed) input, also give the activity of every progeny nuclide. */
    bool show_progeny = false;

    /** Display activities in curie (true) or becquerel (false). */
    bool use_curie = true;

    bool include_activity = true;
    bool include_xrays = false;
    bool include_gammas = false;
    bool include_alphas = false;
    bool include_betas = false;
  };//struct BatchDecayOptions


  /** The computed result table. */
  struct BatchDecayResult
  {
    /** Header row; size == (1 + options.num_steps).  Element 0 is the label-column header. */
    std::vector<std::string> column_headers;

    /** Data rows; each row size == (1 + options.num_steps).  Element 0 is the row label. */
    std::vector<std::vector<std::string>> rows;

    /** Number of numeric data cells (rows * num_steps); used for the GUI preview size cap. */
    std::size_t num_data_cells = 0;

    /** Non-fatal notes accumulated while computing (e.g. skipped stable nuclides). */
    std::string warnings;
  };//struct BatchDecayResult


  /** Decays `inputs` per `opts` and returns the result table.

   Throws std::runtime_error on invalid options (e.g. no valid inputs, non-positive time span).
   Individual invalid/stable inputs are reported in `BatchDecayResult::warnings` rather than throwing.
   */
  BatchDecayResult decay( const std::vector<BatchNuclide> &inputs,
                          const BatchDecayOptions &opts );


  /** Serializes a result table to CSV text (fields comma separated, rows separated by "\r\n"). */
  std::string result_to_csv( const BatchDecayResult &result );


  /** Parses CSV/TSV text into a list of initial nuclides.

   Two formats are auto-detected:
     - Header-keyed: the first non-empty row contains columns named (case-insensitive) "Product",
       "Value", and "Unit" (any order; extra columns ignored).  Nuclide comes from the leading token
       of "Product", magnitude from "Value", and units from "Unit" (activity units are parsed; any
       areal/volumetric suffix such as "/m2" is carried through as `BatchNuclide::unit_label`).
     - Simple: each line is "nuclide, activity[units]" (comma or tab delimited).  Lines beginning
       with '#' and blank lines are ignored.  Units on the activity are optional (default becquerel).

   Throws std::runtime_error with a human-readable message on malformed input.
   */
  std::vector<BatchNuclide> parse_csv( const std::string &file_contents );

}//namespace DecayBatchCalc

#endif //DecayBatchCalc_h
