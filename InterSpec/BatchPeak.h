#ifndef BatchPeak_h
#define BatchPeak_h
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
#include <deque>
#include <memory>
#include <string>
#include <vector>

#include "InterSpec/BatchSampleSelect.h"
#include "InterSpec/DetectionLimitCalc.h"

// Forward declarations
class PeakDef;
class SpecMeas;
struct DetectorPeakResponse;

namespace SpecUtils
{
  class Measurement;
  class EnergyCalibration;
}//namespace SpecUtils

namespace SandiaDecay
{
  struct Nuclide;
}//namespace SandiaDecay


/** Whether the width of the peak region used for detection limits can be set by the user.

 Off by default, so batch analyses match the fixed 2.5 FWHM (i.e., the peak mean +-1.25 FWHM) that
 `DetectionLimitSimple` and `DetectionLimitTool` use, and that ISO 11929:2010 recommends - a limit
 from a batch report is then directly comparable to one from the GUI tools.

 Setting this to 1 exposes the width as a command-line option and a control in the batch GUI.  Note
 that the calculation assumes the entire peak is inside the region, so a narrower region gives an
 optimistically low limit unless that is corrected for; see
 `NotFitPeakMda::signal_fraction_in_roi`.

 This lives here, rather than in `BatchGuiAnaWidget.h`, because the command line
 (`BatchCommandLine.cpp`) and the report JSON (`BatchInfoLog.cpp`) implement the option too.
 */
#define ALLOW_SPECIFY_MDA_ROI_WIDTH 0


namespace BatchPeak
{
  /** The width of the peak region used for detection limits, in multiples of the peak FWHM.

   The value ISO 11929:2010 recommends, and what the GUI detection-limit tools use.
   */
  static const double sm_default_mda_roi_num_fwhm = 2.5;

  /** Fits the peaks in an 'exemplar' file, in a number of other similar file.
   A work in progress.
   
   Currently the results aren't too great, could use:
   - Limits, that are maybe adjustable, for how significant a peak has to be, before it is kept
   - The FWHM needs to be enforced to be reasonable - like after fitting the peaks, fit the FWHM, and then constrain the FWHM to be
      something reasonable around that (hopefully the higher-statistics peaks will dominate the FWHM - or instead could do auto
      peak search, and use those to get the FWHM - this is probably the better idea - or maybe some combination)
   - The FWHM within a single ROI should be constrained
   
   TODO:
   - specify a struct that contains all the options for the fit - we will likely be adding more options, like statistical significance
   - have a .ini file able to back the command line options, so users can specify their own default options
   - When peaks are not fit for, print out their Currie detection limit
   */
  
  /** How to compute the detection limit ("MDA") for exemplar peaks that could not be fit in a
   spectrum.

   The Currie-style (ISO 11929 gross-counts) limit is cheap, and only needs the spectrum, so it is
   the default.  The deconvolution-style limit additionally fits the peak shape at a range of
   source strengths, so it is considerably slower, and is opt-in; when selected, both limits are
   computed.
   */
  enum class NotFitPeakMdaMethod : int
  {
    /** Dont compute detection limits for peaks that were not fit. */
    None,

    /** Currie-style (ISO 11929) gross-counts limit only. */
    Currie,

    /** Currie-style limit, plus the (slower) deconvolution-style limit. */
    CurrieAndDecon
  };//enum class NotFitPeakMdaMethod

  /** Returns "none", "currie", or "currie-and-decon". */
  InterSpec_API const char *to_str( const NotFitPeakMdaMethod method );

  /** Converts "none", "currie", or "currie-and-decon" to a `NotFitPeakMdaMethod`.
   Throws exception if `str` isnt a valid value.
   */
  InterSpec_API NotFitPeakMdaMethod not_fit_peak_mda_method_from_str( const std::string &str );


  struct InterSpec_API BatchPeakFitOptions
  {
    /** If specified to be true, then instead of looking for just the peaks in the exemplar file, a search for all peaks in the specrtum will be performed. */
    bool fit_all_peaks;
    bool to_stdout;
    bool refit_energy_cal;
    bool use_exemplar_energy_cal;
    bool write_n42_with_results;
    bool show_nonfit_peaks;
    bool overwrite_output_files;
    bool create_csv_output;
    bool create_json_output;
    bool concatenate_to_n42;

    /** How to treat an input file that holds more than one candidate foreground spectrum.
     Defaults to the historical behaviour, so a caller that doesnt set it is unaffected.
     */
    BatchSampleSelect::MultiSampleHandling multi_sample_handling
                                    = BatchSampleSelect::MultiSampleHandling::Auto;

    std::string output_dir;
    std::string background_subtract_file;
    std::set<int> background_subtract_samples;
    std::shared_ptr<SpecMeas> cached_background_subtract_spec;
    bool use_existing_background_peaks;
    bool use_exemplar_energy_cal_for_background;
    // TODO: right now there is no option to refit energy calibration of background
    
    /** The improvement to the Chi2 of a peak fit required, over just fitting the continuum, to the ROI.
     
     A negative or zero value indicates no requirement (and default, since we are asserting peak
     is likely in the spectrum for batch analysis), and for general peak searching, reasonable
     values are between ~1 (a weak peak) and ~5 (a significant peak).
     */
    double peak_stat_threshold;
    
    /** Specifies how well the peak must match in shape to a gaussian in order to keep the peak.
     
     The higher this number, the more like a gaussian the fit peak is.
     It is the ratio of the null hypothesis chi2 (continuum only, no Gaussian),
     to the test hypothesis (continuum + Gaussian) chi2.
     A reasonable value for this seems to be ~4.
     A zero or negative value will mean no requirement, and also no
     `peak_stat_threshold` requirement.
     */
    double peak_hypothesis_threshold;

    /** How to compute the detection limit for exemplar peaks that could not be fit.

     Not applicable when `fit_all_peaks` is true (there are no exemplar peaks to miss).
     \sa BatchPeakFitResult::not_fit_peak_mdas
     */
    NotFitPeakMdaMethod not_fit_peak_mda = NotFitPeakMdaMethod::Currie;

    /** The confidence level, as a fraction in (0.5, 1.0), used for the not-fit peak detection
     limits; e.g., 0.95 for 95% CL.
     */
    double mda_confidence_level = 0.95;

    /** Number of channels on each side of the peak region used to estimate the continuum for the
     Currie-style detection limit.  Must be at least one; four is typical.
     */
    size_t mda_num_side_channels = 4;

    /** The width of the peak region used for the detection limit, in multiples of the peaks FWHM.

     Only user-settable when `ALLOW_SPECIFY_MDA_ROI_WIDTH` is enabled; otherwise it is left at
     `sm_default_mda_roi_num_fwhm`, so batch limits match the GUI detection-limit tools.

     A narrower region admits less continuum, and so can give a smaller limit, but note that the
     calculation assumes the whole peak is inside the region - see
     `NotFitPeakMda::signal_fraction_in_roi`.
     */
    double mda_roi_num_fwhm = sm_default_mda_roi_num_fwhm;

    /** The directory to allow report template to look in to include other templates.
     If specified, then the standard report directory cant be used.
     */
    std::string template_include_dir;
    
    
    /** File paths to report templates, that will be saved for each input files. 
     
      If the string contains the string ':--DisplayName--:', then everything before this string will
      be the path to the template, and everything after will be the display name of the template.
      \sa sm_report_display_name_marker
    */
    std::vector<std::string> report_templates;
    
    /** File path to report templates that summarizes all input files. 
     
     Similar to `report_templates`, the string may be delimited by ':--DisplayName--:', to specify 
     the display name of the template.
     \sa sm_report_display_name_marker
    */
    std::vector<std::string> summary_report_templates;

    /** Delimeter used within report template filesystem paths to seperate the filesystem path, and the display name of the template.
     If this delimeter is not present, then the full string is used for both these properties.
     
     This is used for report templates uploaded via HTML, and spooled to disk, so we can still generate reasonably named reports.
     */
    static const char * const sm_report_display_name_marker; // = ":--DisplayName--:"
  };//struct BatchPeakFitOptions

  /** Detection limit ("MDA") information for a single exemplar peak that could not be fit in the
   spectrum being analyzed.

   The Currie-style quantities are always in counts.  The activity quantities are only filled out
   for activity/shielding fits, where the peaks nuclide was one of the fitted sources - see
   #has_activity.
   */
  struct InterSpec_API NotFitPeakMda
  {
    /** The exemplar peak this limit is for; an entry of `BatchPeakFitResult::unfit_exemplar_peaks`. */
    std::shared_ptr<const PeakDef> exemplar_peak;

    /** Classification of the Currie-style result; the same classification used by the GUI
     detection-limit tools.
     */
    enum class MdaResultType : int
    {
      /** Signal is below the decision threshold (L_c) - only an upper limit can be given. */
      NotDetected,

      /** Signal is above the decision threshold, even though a peak could not be fit. */
      Detected,

      /** Significantly fewer counts observed than the neighboring continuum predicts. */
      Deficit,

      /** The limit could not be computed; see #currie_error. */
      Error
    };//enum class MdaResultType

    bool currie_computed = false;
    std::string currie_error;
    DetectionLimitCalc::CurrieMdaResult currie_result;
    MdaResultType result_type = MdaResultType::Error;

    /** Whether a peak that _was_ fit, or another exemplar peak that was not fit, falls within the
     energy range used to estimate the continuum - in either case the limit may be biased.
     */
    bool overlaps_fit_peak = false;
    bool overlaps_other_unfit_peak = false;

    /** The fraction of a Gaussian peak that lies inside the peak region the limit was computed
     over; 0.997 at the default region width of 2.5 FWHM.

     The Currie-style calculation assumes the entire peak is inside the region, so the limits are
     optimistic by roughly this factor.  At the default width the effect is negligible, but it
     grows quickly for narrower regions - e.g., a 1.0 FWHM wide region holds only 76% of the peak.
     \sa BatchPeakFitOptions::mda_roi_num_fwhm
     */
    double signal_fraction_in_roi = 1.0;

    /** A brief English phrase for the result, short enough for a table cell; e.g. "Less than Lc".

     Use `result_type` to branch on in templates, and this to display.
     */
    std::string short_description;

    /** A ready-to-include-in-a-report English description of the result; e.g., "Not detected...".

     Assembled from the three parts below by `update_description(...)`; set that way rather than
     appended to directly, so the caveats stay at the end of the paragraph even though the
     activity information is filled in later, by the activity/shielding fit.
     */
    std::string description;

    /** What the limit is, in counts. */
    std::string result_summary;

    /** What the limit is as an activity, or why there isnt one; empty for plain peak fits. */
    std::string activity_summary;

    /** Notes about anything that makes the limit less reliable; may be empty. */
    std::string caveats;

    /** The deconvolution-style limit; only computed when
     `BatchPeakFitOptions::not_fit_peak_mda == NotFitPeakMdaMethod::CurrieAndDecon`.

     For plain peak fits the quantity limited is peak counts; for activity/shielding fits it is
     the source activity - see #decon_quantity_is_counts.
     */
    bool decon_attempted = false;
    bool decon_computed = false;
    bool decon_quantity_is_counts = true;
    std::string decon_error;
    std::shared_ptr<const DetectionLimitCalc::DeconActivityOrDistanceLimitResult> decon_result;

    /** The deconvolution limit divided by the Currie-style limit; zero if not both available.

     The two methods look at the same data, so they should broadly agree.  When they do not, it is
     usually because the continuum under the peak is not straight, or another peak falls in the
     region - i.e. exactly where the gross-counts assumption breaks down.  A large disagreement is
     therefore worth surfacing rather than silently reporting two different numbers.
     \sa methods_disagree
     */
    double decon_over_currie_ratio = 0.0;

    /** Whether the two limits differ by more than a factor of two. */
    bool methods_disagree = false;

    /** The counts-to-activity conversion; only filled out by `BatchActivity`, and only when the
     peaks nuclide was one of the fitted sources.
     */
    bool has_activity = false;
    std::string no_activity_reason;
    const SandiaDecay::Nuclide *nuclide = nullptr;
    double branching_ratio = 0.0;
    double shield_transmission = 1.0;
    double air_transmission = 1.0;
    double det_efficiency = 0.0;
    double live_time = 0.0;

    /** Counts observed per becquerel of source activity; `activity = counts / gammas_per_bq`. */
    double gammas_per_bq = 0.0;

    /** Whether activities should be displayed in becquerel, rather than curie. */
    bool use_bq = false;

    /** Postfix to add to displayed activities, for fixed-geometry detector responses;
     e.g., "/cm2".  Empty for the usual case.
     */
    std::string activity_postfix;
  };//struct NotFitPeakMda


  struct InterSpec_API BatchPeakFitResult
  {
    std::string file_path;
    BatchPeakFitOptions options;

    std::shared_ptr<const SpecMeas> exemplar;
    std::set<int> exemplar_sample_nums;
    std::deque<std::shared_ptr<const PeakDef>> exemplar_peaks;
    std::shared_ptr<const SpecUtils::Measurement> exemplar_spectrum;
    std::vector<std::shared_ptr<const PeakDef>> unfit_exemplar_peaks;  //Exemplar peaks not found in the spectrum

    /** Detection limits for each entry of `unfit_exemplar_peaks`, in the same order.

     Empty if detection limits were not computed - e.g., `BatchPeakFitOptions::not_fit_peak_mda`
     is `None`, all exemplar peaks were fit, or the fit failed.
     */
    std::vector<NotFitPeakMda> not_fit_peak_mdas;

    std::shared_ptr<SpecMeas> measurement;
    std::shared_ptr<SpecUtils::Measurement> spectrum;
    std::set<int> sample_numbers;
    std::deque<std::shared_ptr<const PeakDef>> fit_peaks;

    /** Background spectrum that was subtracted from the foreground, to make `spectrum`, if any.

     The background subtraction can either be on a peak-by-peak basis, or a hard
     background subtraction, see `BatchPeakFitOptions::use_exemplar_energy_cal_for_background`.
     */
    std::shared_ptr<SpecUtils::Measurement> background;

    bool success;
    std::vector<std::string> warnings;

    /** The original energy calibration of the spectrum, before re-fitting it (if done). */
    std::shared_ptr<const SpecUtils::EnergyCalibration> original_energy_cal;

    /** The energy calibration after fitting for it - will only be non-null if energy calibration was performed */
    std::shared_ptr<const SpecUtils::EnergyCalibration> refit_energy_cal;
  };//struct BatchPeakFitResult


  struct BatchPeakFitSummary
  {
    BatchPeakFitOptions options;
    std::string exemplar_filename;
    std::shared_ptr<const SpecMeas> exemplar;
    std::set<int> exemplar_sample_nums;

    /** Each of these next four `file_*` variables will have the same number of entries as the number
     of analyses performed; this equals the number of input files, unless
     `BatchPeakFitOptions::multi_sample_handling` split an input file into per-sample analyses. */
    std::vector<BatchPeakFitResult> file_results;
    std::vector<std::string> file_json;
    std::vector<std::string> file_peak_csvs;
    std::vector<std::vector<std::string>> file_reports;

    std::string summary_json;
    std::vector<std::string> summary_reports;
    std::vector<std::string> warnings;
    std::shared_ptr<SpecMeas> concatenated_results;
  };//struct BatchPeakFitSummary


  /** Fits the peaks in a number of spectrum files, prodicing individual reports (if wanted) as well as a summary report.
   
   @param exemplar_filename The filesystem path to the spectrum file to use as the exemplar file.  This may be an N42 
          or a peaks CSV file.
   @param optional_parsed_exemplar_n42 If non-null, this object will be used as the exemplar file, with 
          `exemplar_filename` then serving as just the display name of the exemplar file to use in the reports.
   @param exemplar_sample_nums If a N42-2012 file is used for exemplar, and which peaks to use is ambiguous, these
          sample numbers specify which peaks to use.  Must be blank if exemplar is CSV file, or if N42-2012 file, this
          combination of sample numbers must specify peaks to use.  Can be empty in non-ambigous.
   @param files The list of spectrum files to fit peaks to.
   @param optional_cached_files If non-empty, then these objects will be used as the input files for analysis, and 
          `files` will only be used as the display names of the input files in the reports.  If non-empty, then
          must be exactly the same length as `files`.
   @param options The options to use for fitting peaks.
   @param results If non-null, then the results will be stored in this object.
   */
  InterSpec_API void fit_peaks_in_files( const std::string &exemplar_filename,
                          std::shared_ptr<const SpecMeas> optional_parsed_exemplar_n42,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          std::vector<std::shared_ptr<SpecMeas>> optional_cached_files,
                          const BatchPeakFitOptions &options,
                          BatchPeakFitSummary * const results = nullptr );


  
  
  /** Fits the exemplar peaks for a given file.
   
   @param exemplar_filename The file-path of the N42-2012 file with the example peaks, or the file-path of the CSV with peak info
   @param exemplar_sample_nums If a N42-2012 file is used for exemplar, and which peaks to use is ambiguous, these
          sample numbers specify which peaks to use.  Must be blank if exemplar is CSV file, or if N42-2012 file, this combination
          of sample numbers must specify peaks to use.
   @param cached_exemplar_n42 If non-null, then `exemplar_filename` will be ignored, and this file will be used; to avoid re-parsing
          of the exemplar file over-and-over again.
   @param filename The name of the spectrum file to fit peaks to.
   @param cached_spectrum If you have already parsed/opened the `filename` spectrum file, you can provide it here to
          avoid overhead of re-parsing it.
   @param foreground_sample_numbers The sample numbers to fit the peaks to.  If left empty, will try to automatically determine.
   @param options The options to use for fitting peaks; note, not all options are used, as some of them are only applicable to
          #fit_peaks_in_files
   */
  InterSpec_API BatchPeakFitResult fit_peaks_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          std::shared_ptr<SpecMeas> cached_spectrum,
                          std::set<int> foreground_sample_numbers,
                          const BatchPeakFitOptions &options );
  
  /** Computes the detection limit ("MDA"), in counts, for each exemplar peak that could not be fit.

   The Currie-style (ISO 11929 gross-counts) limit uses a peak region of the exemplar peaks mean
   +-1.25 FWHM, with `BatchPeakFitOptions::mda_num_side_channels` channels on either side of it to
   estimate the continuum.  If `BatchPeakFitOptions::not_fit_peak_mda` is `CurrieAndDecon`, then
   the deconvolution-style limit, also in counts, is computed as well.

   Never throws; a peak whose limit could not be computed (e.g., its region falls off the end of
   the spectrum) is returned with `NotFitPeakMda::result_type` of `Error`, and the reason in
   `NotFitPeakMda::currie_error`.

   @param unfit_exemplar_peaks The exemplar peaks that could not be fit; the returned vector has
          one entry per peak, in the same order.
   @param fit_peaks The peaks that _were_ fit; used only to flag limits whose continuum estimate
          may be contaminated by a nearby fit peak.
   @param spectrum The spectrum that was analyzed.
   @param exemplar_drf The detector response of the exemplar file; may be null.  Only used to
          provide the peak widths for the deconvolution-style limit, and even then only if it has
          resolution information (otherwise the exemplar peak widths are used).
   @param options The options for the analysis.
   */
  InterSpec_API std::vector<NotFitPeakMda> compute_not_fit_peak_mdas(
                          const std::vector<std::shared_ptr<const PeakDef>> &unfit_exemplar_peaks,
                          const std::deque<std::shared_ptr<const PeakDef>> &fit_peaks,
                          const std::shared_ptr<const SpecUtils::Measurement> &spectrum,
                          const std::shared_ptr<const DetectorPeakResponse> &exemplar_drf,
                          const BatchPeakFitOptions &options );

  /** Runs the deconvolution-style detection limit calculation, and stores the outcome into `mda`.

   Never throws; a failure is recorded in `NotFitPeakMda::decon_error`.

   @param mda The limit information to fill out.  `NotFitPeakMda::currie_result` must already have
          been computed, as it is used to seed the range of quantities searched over.
   @param input The calculation input; the ROI information, spectrum, and detector response must
          already be filled out.
   @param gammas_per_bq The number of counts observed per unit of the quantity being limited; use
          1.0 when the quantity being limited is counts.
   @param quantity_is_counts Whether the quantity being limited is peak counts (true), or source
          activity (false).
   @param confidence_level The confidence level, as a fraction; e.g., 0.95.
   @param use_curie Whether the text strings the calculation produces should use curie, rather
          than becquerel.  Only matters when the quantity being limited is an activity.
   */
  InterSpec_API void compute_decon_limit( NotFitPeakMda &mda,
                          const DetectionLimitCalc::DeconComputeInput &input,
                          const double gammas_per_bq,
                          const bool quantity_is_counts,
                          const double confidence_level,
                          const bool use_curie );

  /** Formats a confidence level as a string for display; e.g., 0.95 -> "95%", and
   0.999999426696856 -> "1-5.7E-07".
   */
  InterSpec_API std::string confidence_level_str( const double confidence_level );

  /** Returns the fraction of a Gaussian peak that lies within a region centered on the peak, and
   `num_fwhm` FWHMs wide; e.g., 2.5 -> 0.99676.
   */
  InterSpec_API double gaussian_fraction_in_roi( const double num_fwhm );

  /** Computes the deconvolution-style detection limit, in counts, for each entry that does not
   already have one.

   Entries whose `decon_attempted` is already true are left alone, so an activity/shielding fit can
   compute the (better) activity-space limit first, and then use this to fill in the peaks it could
   not convert to an activity - rather than leaving them with no limit at all.

   Never throws; per-entry failures are recorded in `NotFitPeakMda::decon_error`.
   */
  InterSpec_API void add_counts_decon_limits( std::vector<NotFitPeakMda> &mdas,
                          const std::shared_ptr<const SpecUtils::Measurement> &spectrum,
                          const std::shared_ptr<const DetectorPeakResponse> &drf,
                          const BatchPeakFitOptions &options );

  /** Rebuilds `NotFitPeakMda::description` from its parts; call after changing any of them. */
  InterSpec_API void update_description( NotFitPeakMda &mda );


  /** One analysis' contribution to the concatenated N42 file.
   \sa write_concatenated_n42
   */
  struct InterSpec_API ConcatRecord
  {
    /** The input file this record came from; used for the record title, and to detect input files
     that produced more than one analysis. */
    std::string source_file_path;

    /** Sample numbers analyzed; only spelled out in the title when the input file was split. */
    std::set<int> sample_numbers;

    std::shared_ptr<const SpecUtils::Measurement> spectrum;
    std::deque<std::shared_ptr<const PeakDef>> peaks;
  };//struct ConcatRecord


  /** Writes `<output_dir>/concatenated.n42` holding the spectrum and fit peaks of every analysis.

   Records are accumulated by the caller as it goes, rather than taken from a
   `BatchPeakFitSummary`, so that this works whether or not the caller asked for the full results
   to be retained.

   Does nothing if `records` is empty.  Any problem is reported by appending to `warnings` rather
   than throwing.

   @returns The concatenated file, or nullptr if nothing was written.
   */
  InterSpec_API std::shared_ptr<SpecMeas>
    write_concatenated_n42( const std::vector<ConcatRecord> &records,
                            const BatchPeakFitOptions &options,
                            std::vector<std::string> &warnings );


  /** Function that applies the energy calibration from the exemplar spectrum, to a spectrum from a different file.
   
   @param energy_cal The energy calibration to apply to `to_spectrum`, and optionally `to_specfile`
   @param to_spectrum The spectrum, which may or may not be in `to_specfile`, to apply the energy calibration from `from_spectrum`
   @param to_specfile The (optional) spectrum file to apply the energy calibration to; this will also take care of shifting peak energies
   @param used_sample_nums The sample numbers used to create the `to_spectrum` from the `to_specfile` - used to keep peaks from
          being moved twice.
   */
  void propagate_energy_cal( const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                                          std::shared_ptr<SpecUtils::Measurement> &to_spectrum,
                                          std::shared_ptr<SpecMeas> &to_specfile,
                                          const std::set<int> &used_sample_nums );
  
  /** Finds the spectrum, peaks, and sample numbers to use from the exemplar file.
   
   @param [out] exemplar_spectrum Will be set to the spectrum to use as the exemplar spectrum (may be a sum of
          multiple Measurements)
   @param [out] exemplar_peaks Will be set to the peaks to use from the exemplar file.
   @param [in|out] exemplar_sample_nums Sample numbers to use in the exemplar file.  If non-empty, these sample
          number will be used to retrieve the spectrum to use.  Contents will be set to the used sample numbers.
   @param [in] exemplar_n42 The spectrum file to retrieve the spectrum/peaks for
   @param [in] require_peaks If true, and there are no peaks for the specified measurement, an exception will be thrown.

   Throws exception on error (not a valid measurmeent specification, possibly no peaks, etc).
   */
  void get_exemplar_spectrum_and_peaks(
                                       std::shared_ptr<const SpecUtils::Measurement> &exemplar_spectrum,
                                       std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &exemplar_peaks,
                                       std::set<int> &exemplar_sample_nums,
                                       const std::shared_ptr<const SpecMeas> &exemplar_n42,
                                       const bool require_peaks );
}//namespace BatchPeak

#endif //BatchPeak_h
