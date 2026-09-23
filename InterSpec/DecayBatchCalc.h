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
#include <utility>

namespace SandiaDecay
{
  struct Nuclide;
}

/** Non-GUI core for the "Batch Decay" tool.

 Decays potentially many initial nuclides for one, or a number of, time points and produces a
 table of activities (and, optionally, particle-line rates).  This replaces the legacy command-line
 utility at external_libs/SandiaDecay/examples/batch_decay.cpp.

 There are two output orientations:
   - "Wide" (the default): rows = items (nuclides or particle-lines), columns = time steps - matching
     the legacy tool.  The activity/particle units are carried in the row label so different row
     types (Bq/Ci vs counts/second) can coexist in one table.
   - "Long"/grouped: used when the inputs carry locations (see `BatchNuclide::location`), where each
     output row is one (location, nuclide) pair and mirrors the input file's own columns.  See
     #decay for why, and the "Multiple physical locations" section below.

 ## Multiple physical locations

 Some input files hold measurements for many physical locations at once (first column = location).
 Each location is an independent sample set: its nuclides are co-decayed together as one mixture, and
 values from different locations are never summed or mixed.  Setting `BatchNuclide::location`
 switches #decay into this grouped mode.
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

    /** Age of the nuclide at the measurement - time zero of the output - in SandiaDecay time units
     (seconds); i.e. how long its progeny have been growing in.  See #decay for how it is used when
     decaying backwards.
     */
    double age = 0.0;

    /** Initial activity, in SandiaDecay activity units (becquerel). */
    double activity = 0.0;

    /** Optional areal/volumetric suffix to carry through to output labels (e.g. "/m2").
     Empty for a normal activity input.  Purely a display label; the numeric decay is unaffected.
     */
    std::string unit_label;

    /** The physical location this measurement belongs to, or empty when the input is not grouped.

     All inputs sharing a location are co-decayed as one mixture, and never mixed with another
     location's.  A non-empty value on *any* input switches #decay to grouped/long output.
     */
    std::string location;

    /** For grouped input, the text following the nuclide token in the source "Product" cell (e.g.
     " Deposition at 28 hrs"), reproduced in the output label.  Empty when not applicable.
     */
    std::string product_suffix;

    /** For grouped input, extra source columns (name -> value, in source order; e.g. Latitude,
     Longitude) echoed verbatim into the output row.  This core does not interpret them.
     */
    std::vector<std::pair<std::string,std::string>> extra_columns;

    /** For grouped input, the activity unit token the value was given in (e.g. "uCi"), so grouped
     output can be written back in the input's own units rather than the user's Ci/Bq preference.
     Empty when not applicable.
     */
    std::string activity_unit;

    /** For grouped input, the source file's own names for the location, product, value and unit
     columns (in that order), so the output header echoes the input rather than hard-coded English.
     Empty when not applicable, in which case #decay falls back to the legacy English names.
     */
    std::vector<std::string> fixed_column_names;
  };//struct BatchNuclide


  /** Decay options; mirrors the legacy CLI arguments and the CSV-export dialog options. */
  struct BatchDecayOptions
  {
    /** Time to decay to, in SandiaDecay time units (seconds).  Must be non-zero.

     A *negative* value decays backwards in time: the inputs are taken as present-day measurements
     and the activities they must have had `|time_span|` ago are solved for.  See #decay for the
     caveats this brings (inputs decayed together that share an ancestor are coupled, and some past
     activities are not recoverable at all).
     */
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

    /** Drop results whose activity is below this value, as printed in its own displayed unit; zero
     (the default) applies no cut.  Unit-less on purpose: it is compared against the number the user
     sees, so a cut of 1e-20 against values printed in "uCi/m2" means 1e-20 uCi/m2.

     In grouped/long output an individual (location, nuclide) row is dropped, so a location whose
     every row is below the cut disappears entirely.  In wide output a row spans several time steps,
     so a row is dropped only when *every* step is below the cut (keeping the table rectangular).

     Applies to activity rows only: the particle-line rows are counts/second, which an activity
     threshold cannot meaningfully be compared against, so they are always kept.
     */
    double min_activity = 0.0;

    /** The decay time as the user typed it (e.g. "48h"), used verbatim to build grouped output
     labels like "Cs137 ...+48h".  Optional; when empty a compact form of `time_span` is used.
     */
    std::string time_span_str;
  };//struct BatchDecayOptions


  /** The computed result table.

   Two shapes, depending on whether the inputs carried locations (see `BatchNuclide::location`):
     - Wide (ungrouped): `column_headers.size() == (1 + options.num_steps)`, element 0 being the
       label column; likewise each row, whose element 0 is the row label.
     - Long (grouped): the columns mirror the input file's own (location, product, ...extras,
       value, unit), so there is a single value column and `num_steps > 1` emits one block of rows
       per time step.
   Either way element 0 of a row is textual and the trailing cells are numeric.
   */
  struct BatchDecayResult
  {
    /** Header row. */
    std::vector<std::string> column_headers;

    /** Data rows; every row has `column_headers.size()` entries. */
    std::vector<std::vector<std::string>> rows;

    /** Number of cells the GUI would render; used for its preview size cap. */
    std::size_t num_data_cells = 0;

    /** Non-fatal notes accumulated while computing (e.g. skipped stable nuclides). */
    std::string warnings;
  };//struct BatchDecayResult


  /** Decays `inputs` per `opts` and returns the result table.

   When any input carries a `location`, the inputs are grouped by it: each location is co-decayed as
   its own mixture and reported as its own block of rows, in the long/grouped output shape.  Values
   from different locations are never combined.

   A negative `opts.time_span` decays backwards in time, solving for the activities the inputs must
   have had `|time_span|` ago.  Because inputs decayed together (mixed, or of one location) that share
   an ancestor are coupled (e.g. Cs137 and its Ba137m progeny), this is a coupled inverse problem
   rather than a per-nuclide division, and two things follow: a short-lived nuclide's own past
   activity may be *unrecoverable* (nothing observable today depends on it), in which case it is
   reported as zero if an ancestor decayed with it accounts for the measurement and throws if none
   does; and if the inputs are not mutually consistent with having decayed from a common past state,
   the recovered state cannot reproduce them exactly - which is noted in `warnings`.

   An input's `age` is its age at the measurement, so looking back `|time_span|` its mixture is seeded
   `|time_span|` younger.  A (non-zero) age shorter than that is taken as freshly made at the past
   time, with a note in `warnings`; an age of zero (none given) means the same, without the note.

   Throws std::runtime_error on invalid options (e.g. no valid inputs, zero time span), and when a
   measured nuclide is so many half-lives old that its past activity cannot be recovered at all (the
   message names the nuclide and asks for a shorter time).
   Individual invalid/stable inputs are reported in `BatchDecayResult::warnings` rather than throwing.
   */
  BatchDecayResult decay( const std::vector<BatchNuclide> &inputs,
                          const BatchDecayOptions &opts );


  /** The location key for a source row's first-column value; see `BatchNuclide::location`.

   Rows of one location differ only by a trailing "_<something>" that indexes the row within the
   file, so that suffix is dropped when it is all digits, or when it contains a digit *and* the text
   before it looks like a structured ID (contains a <alnum>-<alnum>, e.g. "TeamA-01_run3" ->
   "TeamA-01").  Any other name is returned unchanged, so genuinely distinct names like "Site_North" -
   or "Site-A_North", whose suffix names a place rather than indexing one - are never merged.

   Exposed for testing; #parse_csv applies it, and falls back to the raw name if the result would
   merge rows that cannot belong to one location.
   */
  std::string location_key( const std::string &probe_name );


  /** Serializes a result table to CSV text (fields comma separated, rows separated by "\r\n"). */
  std::string result_to_csv( const BatchDecayResult &result );


  /** Parses CSV/TSV text into a list of initial nuclides.

   Three formats are auto-detected:
     - Multi-location: a header-keyed file that also has a leading name/location column (e.g.
       "Probe name,Product,Latitude,Longitude,Value,Unit").  Rows are grouped into locations by
       #location_key, with the remaining columns carried through for the output (see
       `BatchNuclide::location`, `::extra_columns`).  If grouping would merge rows that cannot belong
       to one location (a repeated nuclide, or disagreeing carried-through cells) the raw names are
       used instead; and if no location ends up with more than one nuclide the file is treated as
       ungrouped.
     - Header-keyed: the first non-empty row contains columns named (case-insensitive) "Product",
       "Value", and "Unit" (any order; extra columns ignored).  Nuclide comes from the leading token
       of "Product", magnitude from "Value", and units from "Unit" (activity units are parsed; any
       areal/volumetric suffix such as "/m2" is carried through as `BatchNuclide::unit_label`).
     - Simple: each line is "nuclide, activity[units]" (comma or tab delimited).  Lines beginning
       with '#' and blank lines are ignored.  Units on the activity are optional (default becquerel).

   A leading UTF-8 byte-order mark (as spreadsheet "CSV UTF-8" exports write) is ignored.
   Note `BatchNuclide::age` is not set by any of these formats.

   Throws std::runtime_error with a human-readable message on malformed input.
   */
  std::vector<BatchNuclide> parse_csv( const std::string &file_contents );


  /** Whether some text looks like a batch-decay input file, i.e. #parse_csv accepts it.

   Used to recognize these files when dropped on the main window, where only the start of the file may
   be at hand: with `is_whole_file` false, a trailing partial line is ignored (and text with no line
   end at all is rejected).  Text containing a NUL byte is rejected.
   */
  bool is_candidate_file( const std::string &start_of_file, const bool is_whole_file );

}//namespace DecayBatchCalc

#endif //DecayBatchCalc_h
