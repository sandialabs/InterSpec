#ifndef BatchSampleSelect_h
#define BatchSampleSelect_h
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

#include <set>
#include <memory>
#include <string>
#include <vector>

// Forward declarations
class SpecMeas;

namespace SpecUtils
{
  class SpecFile;
}//namespace SpecUtils


/** Deciding which sample number(s) of a spectrum file to analyze.

 This logic used to be duplicated, with subtly different tie-breaks, in `BatchActivity`,
 `BatchPeak`, and `BatchRelActAuto`; this namespace is the single implementation those now share.
 It also provides the expansion of input files into batch "work items", which is how a file
 holding several foreground records gets analyzed once per record.
 */
namespace BatchSampleSelect
{
  /** Classification of a single sample number of a file.

   Follows the convention of `SpecUtils::SpecFile::sum_measurements()`: the type of a sample is the
   common `SourceType` of the `SpecUtils::Measurement`s making it up, and is `Unknown` when those
   Measurements disagree with each other.  This keeps the answer independent of detector ordering.
   */
  enum class SampleClass : int
  {
    /** Every Measurement of this sample is marked Foreground. */
    Foreground,

    /** Marked Unknown, or the Measurements of this sample disagree with each other. */
    Unknown,

    /** Every Measurement of this sample is marked Background. */
    Background,

    /** Every Measurement of this sample is marked IntrinsicActivity and/or Calibration. */
    Other,

    /** Derived data, in a file that also holds non-derived data. */
    Derived,

    /** No Measurements for this sample number. */
    Empty
  };//enum class SampleClass


  /** All the sample numbers of a file, bucketed by `SampleClass`.

   Callers apply their own ladder to these; e.g. `BatchPeak::fit_peaks_in_file()` is deliberately
   more permissive than `BatchActivity`, since fitting peaks in a lone calibration spectrum is a
   sensible thing to do, while fitting an activity to one is not.

   `SampleClass::Derived` and `SampleClass::Empty` samples appear in none of these sets.
   */
  struct SampleBuckets
  {
    std::set<int> foreground;
    std::set<int> unknown;
    std::set<int> background;
    std::set<int> other;
  };//struct SampleBuckets


  /** Classifies a single sample number; see `SampleClass`. */
  InterSpec_API SampleClass classify_sample( const SpecUtils::SpecFile &file, const int sample_number );

  /** Classifies every sample number of `file`. */
  InterSpec_API SampleBuckets classify_samples( const SpecUtils::SpecFile &file );


  /** Options controlling `candidate_foreground_samples()`. */
  struct Options
  {
    /** If no sample classifies as Foreground or Unknown, but exactly one classifies as Background,
     return that background sample.  This preserves the long-standing behaviour of being able to
     analyze a file whose sole record happens to be marked background.
     */
    bool allow_lone_background_fallback = true;

    /** If true, Unknown samples are returned together with Foreground ones.  If false, Unknown
     samples are used only when there are no Foreground samples at all.
     */
    bool mix_unknown_with_foreground = false;
  };//struct Options


  /** The sample numbers that plausibly make up "the foreground" of a file.

   Returns an empty set when the answer is not determinable.

   IntrinsicActivity, Calibration, and derived-data samples are never returned; Background samples
   are returned only via `Options::allow_lone_background_fallback`.
   */
  InterSpec_API std::set<int> candidate_foreground_samples( const SpecUtils::SpecFile &file,
                                                            const Options &options = Options{} );

  /** Returns the single sample number to use as the foreground.

   Foreground and Unknown samples are considered together, and the lone-background fallback applies.

   Throws `std::runtime_error` if the answer is not unique.
   */
  InterSpec_API int single_foreground_sample( const SpecUtils::SpecFile &file );

  /** Returns the single sample number to use as the background.

   Throws `std::runtime_error` if the answer is not unique.
   */
  InterSpec_API int single_background_sample( const SpecUtils::SpecFile &file );


  /** How to treat an input file that holds more than one candidate foreground spectrum.

   Files that are passthrough/search-mode, or that hold a single candidate foreground sample, are
   never affected by this option.
   */
  enum class MultiSampleHandling : int
  {
    /** Historical behaviour: no sample numbers are passed down, and the fit function auto-detects.
     A file with more than one candidate foreground sample yields a single failed result, carrying
     the "ambiguous which spectrum to use" warning.  Costs no extra file parsing.
     */
    Auto,

    /** One independent analysis per candidate foreground sample; each gets its own entry in the
     summary JSON and its own complete set of output files, whose names gain a "_sampleN" infix.
     */
    EachSampleSeparately,

    /** Sum all candidate foreground samples into one spectrum and perform a single analysis.
     Output file names gain a "_summed" infix.
     */
    SumAllSamples
  };//enum class MultiSampleHandling


  /** Converts to/from the strings used for the command line and the report JSON. */
  InterSpec_API const char *to_str( const MultiSampleHandling handling );

  /** Parses "auto", "each"/"each-sample"/"separate", or "sum"/"sum-all" (case-insensitive).
   Throws `std::runtime_error` on anything else.
   */
  InterSpec_API MultiSampleHandling multi_sample_handling_from_str( std::string val );


  /** One unit of batch work: an input file, optionally narrowed to specific sample numbers. */
  struct InputWorkItem
  {
    /** Index into the caller's `files` vector this item came from. */
    size_t input_index = 0;

    /** Verbatim `files[input_index]` - never decorated.

     Used for the reports `Filepath`/`ParentDir`, as the load-on-demand path, and - critically -
     for the `options.background_subtract_file == filename` test that decides whether the
     background lives in the same file as the foreground.  Decorating this would silently break
     that test for every split work item.
     */
    std::string filename;

    /** Leaf name (with extension) that every output file name is derived from.
     Equals `SpecUtils::filename(filename)` unless this file was split or summed.
     */
    std::string output_base_name;

    /** Short human label used in warning messages, e.g. "foo.n42 (sample 3)". */
    std::string label;

    /** The parsed file to analyze; may be null, in which case the fit function loads `filename`.
     Shared between all work items that came from the same input file - see `needs_private_copy`.
     */
    std::shared_ptr<SpecMeas> source;

    /** True when more than one work item shares `source`.

     The fit functions destructively modify the file handed to them (energy calibration
     propagation, and hard background subtraction replaces all the Measurements), so a driver must
     hand each such work item its own `SpecMeas::uniqueCopyContents()` copy.  Copying is left to
     the driver, rather than done here, so only one copy is alive at a time.
     */
    bool needs_private_copy = false;

    /** The sample numbers to use as the foreground; empty means "auto-detect downstream". */
    std::set<int> foreground_sample_numbers;
  };//struct InputWorkItem


  /** Expands a single input file into the batch work items to analyze for it.

   With `MultiSampleHandling::Auto` this is a straight one-to-one pass-through, and the file is not
   parsed here.  Otherwise the file is parsed (if not already provided in `cached_file`) so its
   candidate foreground samples can be enumerated.

   Callers should prefer this over `expand_input_files()` when processing many files, so that only
   one input file is held parsed in memory at a time.

   @param filename The input file path.
   @param input_index Index of this file in the caller's list; copied into `InputWorkItem`.
   @param cached_file Already-parsed input file, or null.
   @param handling How to treat a file with more than one candidate foreground spectrum.
   */
  InterSpec_API std::vector<InputWorkItem>
    expand_input_file( const std::string &filename,
                       const size_t input_index,
                       const std::shared_ptr<SpecMeas> &cached_file,
                       const MultiSampleHandling handling );

  /** Expands a whole list of input files; see `expand_input_file()`.

   Note that for anything but `MultiSampleHandling::Auto` this parses, and keeps parsed, every
   input file - so it is only appropriate for small lists.

   @param files The input file paths.
   @param cached_files Already-parsed input files; if non-empty must be the same length as `files`.
   @param handling How to treat files with more than one candidate foreground spectrum.
   */
  InterSpec_API std::vector<InputWorkItem>
    expand_input_files( const std::vector<std::string> &files,
                        const std::vector<std::shared_ptr<SpecMeas>> &cached_files,
                        const MultiSampleHandling handling );

}//namespace BatchSampleSelect

#endif //BatchSampleSelect_h
