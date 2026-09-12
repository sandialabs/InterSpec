#ifndef BatchInfoLog_h
#define BatchInfoLog_h
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
#include <set>
#include <deque>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"
#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

// Forward declarations
class PeakDef;
class SpecMeas;
class DetectorPeakResponse;


namespace ShieldingSourceFitCalc
{
  struct ModelFitResults;
  struct SupplementalPeakInfo;
  struct ShieldingSourceFitOptions;
}

namespace DetectionLimitCalc
{
  struct PeakCurrieCheck;
}

namespace GammaInteractionCalc
{
  struct PeakDetail;
  struct PeakDetailSrc;
  struct SourceDetails;
  struct ShieldingDetails;
  struct ShieldingSourceFitOptions;
  enum class GeometryType : int;
}

namespace BatchPeak
{
  struct NotFitPeakMda;
  struct BatchPeakFitResult;
  struct BatchPeakFitOptions;
}

namespace SpecUtils
{
  class Measurement;
  struct EnergyCalibration;
}

namespace BatchActivity
{
  struct BatchActivityFitOptions;
}


namespace BatchInfoLog
{
  /** Returns the equivalent of `InterSpec_resources/static_text/ShieldSourceFitLog/`, with trailing
   slash that inja requires. */
  std::string default_template_dir();
  
  /** Returns the user specified template include directory - or if none/default is specified, then `default_template_dir()`.
   
   Will include trailing path separator, as required by inja.
   */
  std::string template_include_dir( const BatchPeak::BatchPeakFitOptions &options );
  
  /** Returns the default inja environment, with include directory, options, and callbacks set, as well as default templates loaded.

   @param tmplt_dir The template include directory to use.  Should include trailing path separator (as required by inja).
                    If empty string, then no custom template directory will be used.
                    To get the directory from options, use `template_include_dir(options)`.
   */
  inja::Environment get_default_inja_env( const std::string &tmplt_dir );
  
  /** Returns key-value pairs of of the file contents of the JS and CSS files needed for SpectrumChar.  Specifically returns:

   "D3_JS":                        contents of `InterSpec_resources/d3.v3.min.js`
   "SpectrumChart_JS":    contents of `InterSpec_resources/SpectrumChartD3.js`
   "SpectrumChart_CSS": contents of `InterSpec_resources/SpectrumChartD3.css`

   */
  std::vector<std::pair<std::string,std::string>> load_spectrum_chart_js_and_css();

  /** Returns key-value pairs of the file contents of the JS and CSS files needed for ShieldingSourceFitPlot. Specifically returns:

   "D3_JS":                              contents of `InterSpec_resources/d3.v3.min.js`
   "ShieldingSourceFitPlot_JS":  contents of `InterSpec_resources/ShieldingSourceFitPlot.js`
   "ShieldingSourceFitPlot_CSS": contents of `InterSpec_resources/ShieldingSourceFitPlot.css`

   */
  std::vector<std::pair<std::string,std::string>> load_shielding_fit_plot_js_and_css();
  
  /** An enum to provide context of what default templates names "csv", "txt", and "html" refer to for `render_template(...)` */
  enum class TemplateRenderType : int
  {
    ActShieldIndividual,
    ActShieldSummary,
    PeakFitIndividual,
    PeakFitSummary
  };//enum class TemplateRenderType

  /** Report template strings may bundle a filesystem path with a display name, delimited by
   `BatchPeak::BatchPeakFitOptions::sm_report_display_name_marker` (":--DisplayName--:") - this is
   how templates uploaded via the GUI (and spooled to a temporary path) carry a human-readable name.

   `report_template_disk_path(...)` returns the part before the marker (the path to open), and
   `report_template_display_name(...)` returns the part after it (used for naming output files).
   If the marker is not present, the full string is returned unchanged by both.
   */
  std::string report_template_disk_path( const std::string &tmplt );
  std::string report_template_display_name( const std::string &tmplt );
  
  /** Renders data to a template, given the inputs.
   
   @param tmplt The name of the template file to use.  May be a default template given by values "txt" and "html"
          for ActShieldIndividual, or "csv" and "html" for ActShieldSummary.  If a file path, must either be a full
          path to file, or just the filename if in the specified or default template include directories.
   @param env The inja environment used to render the template.  It is expected to have been setup by `get_default_inja_env(...)`
   @param type The context for how to interpret the default template names
   @param options The user options, used to get template include directory.
   @param data The data to be rendered to the template.
   @returns The rendered template contents.
   
   Throws exception on error; `std::exception` on not finding template, `inja::InjaError` on templating error.
   */
  std::string render_template( std::string tmplt,
                              inja::Environment &env,
                              const TemplateRenderType type,
                              const BatchPeak::BatchPeakFitOptions &options,
                              const nlohmann::json &data );
  
  
  /** Combines the data filename, with report template name, and appends to report output directory, to return
   the suggested name and path for the output report file.
   */
  std::string suggested_output_report_filename( const std::string &filename,
                                               const std::string tmplt,
                                               const TemplateRenderType type,
                                               const BatchPeak::BatchPeakFitOptions &options );
  
  /** Callback from inja templating to print a floating point number to a specified number of decimals.
   Takes two arguments, the first is a `double` for the value to print, the second is an int, for the number of decimals to print.
   */
  std::string printFixed( std::vector<const nlohmann::json *> &args );
  
  /** Callback from inja templating to print a floating point number using `SpecUtils::printCompact(...)`. */
  std::string printCompact( std::vector<const nlohmann::json *> &args );

  /** Callback from inja templating to print a floating point number in scientific notation (e.g. "1.891E+04").
   Takes two arguments: the first is a `double` value, the second is an int giving the number of decimals.
   Matches the FRMAC/Genie report format `d.dddE+NN`.
   */
  std::string printExp( std::vector<const nlohmann::json *> &args );
  
  void add_basic_src_details( const GammaInteractionCalc::SourceDetails &src,
                            const std::shared_ptr<const DetectorPeakResponse> &drf,
                            const bool useBq,
                            const std::vector<GammaInteractionCalc::ShieldingDetails> *shield_details,
                             nlohmann::basic_json<> &src_json );
  
  void add_act_shield_fit_options_to_json( const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options,
                               const double distance,
                               const GammaInteractionCalc::GeometryType geometry,
                               const std::shared_ptr<const DetectorPeakResponse> &drf,
                                          nlohmann::basic_json<> &data );
  
  /** Adds basic information about a peak (energy, fwhm, counts, etc), but not any information
     about gammas that contribute to it, etc
   */
  void add_basic_peak_info( const GammaInteractionCalc::PeakDetail &peak, nlohmann::basic_json<> &peak_json );
  
  void shield_src_fit_results_to_json( const ShieldingSourceFitCalc::ModelFitResults &results,
                                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                                      const bool useBq,
                                      nlohmann::basic_json<> &data );
  
  void add_gamma_info_for_peak( const GammaInteractionCalc::PeakDetailSrc &ps,
                    const GammaInteractionCalc::SourceDetails * const src,
                    const std::shared_ptr<const DetectorPeakResponse> &drf,
                    const bool useBq,
                    const std::vector<GammaInteractionCalc::ShieldingDetails> * const shield_details,
                               nlohmann::basic_json<> &gamma_json );
  
  
  void add_hist_to_json( nlohmann::basic_json<> &data,
                       const bool is_background,
                       const double display_scale_factor,
                       const std::shared_ptr<const SpecUtils::Measurement> &spec_ptr,
                       const std::shared_ptr<const SpecMeas> &spec_file,
                       const std::set<int> &sample_numbers,
                       const std::string &filename,
                       const std::deque<std::shared_ptr<const PeakDef>> * const peak_fit );
  
  void add_energy_cal_json( nlohmann::basic_json<> &data,
                           const std::shared_ptr<const SpecUtils::EnergyCalibration> &cal );
  
  void add_activity_fit_options_to_json( nlohmann::basic_json<> &data,
                                        const BatchActivity::BatchActivityFitOptions &options );
  
  void add_exe_info_to_json( nlohmann::basic_json<> &data );
  
  void add_peak_fit_options_to_json( nlohmann::basic_json<> &data, const BatchPeak::BatchPeakFitOptions &options );
  
  void add_peak_fit_results_to_json( nlohmann::basic_json<> &data, const BatchPeak::BatchPeakFitResult &fit_results );

  /** Fills out `json` with the "Continua", "Peaks", and "PeakSortIndex_*" entries for `peaks`.

   @param mdas If non-null, then each peak gets a "HasMda" entry, and if it has a matching entry
          (matched by pointer), an "Mda" object.  Pass nullptr for collections of peaks that were
          fit, so that reports written before detection limits were added see no change.
   */
  void add_peaks_to_json( nlohmann::json &json,
                         std::deque<std::shared_ptr<const PeakDef>> peaks,
                         const std::shared_ptr<const SpecUtils::Measurement> &spectrum,
                         const std::vector<BatchPeak::NotFitPeakMda> * const mdas );

  /** Fills out `json` with the exemplar peaks that could not be fit, and their detection limits.

   Always sets `json["HasMdas"]`, so templates can rely on it being there.
   */
  void add_not_fit_peaks_to_json( nlohmann::json &json,
                                 const BatchPeak::BatchPeakFitResult &fit_results );

  /** Adds the detection limit information for a single peak that could not be fit. */
  void add_mda_to_json( nlohmann::basic_json<> &mda_json, const BatchPeak::NotFitPeakMda &mda );

  /** Adds the energy, width, area, and assigned source of a peak. */
  void add_peak_identity_to_json( const PeakDef &peak, nlohmann::basic_json<> &peak_json );

  /** Adds the assigned source ("SourceType", "SourceName", "SourceEnergy", ...) of a peak. */
  void add_peak_source_info_to_json( const PeakDef &peak, nlohmann::basic_json<> &peak_json );

  /** Adds the counts-space quantities of a Currie-style detection limit check. */
  void add_currie_check_to_json( nlohmann::basic_json<> &json,
                                const DetectionLimitCalc::PeakCurrieCheck &check );

  /** Adds the "SupplementalPeakInfo" object to an activity/shielding fit result: the per-peak
   detection limit checks, and the activities implied by peaks that were fit but not used in the
   model.

   Called by `shield_src_fit_results_to_json(...)`, so the GUIs fit log and the batch reports both
   get it.  Emits nothing at all when `supp_info` is empty, so reports written before this existed
   render unchanged.
   */
  void add_supplemental_peak_info_to_json( nlohmann::basic_json<> &data,
                  const std::vector<ShieldingSourceFitCalc::SupplementalPeakInfo> &supp_info,
                  const std::shared_ptr<const DetectorPeakResponse> &drf,
                  const bool useBq );

  /** Adds the "NotFitPeaks" and "AnyNotFitPeakMda" entries to an activity/shielding fit result.

   The peak-fit results are not otherwise included in activity/shielding fit reports, so this is
   what makes the detection limits available to those templates.
   */
  void add_not_fit_peaks_to_act_shield_json( nlohmann::basic_json<> &data,
                                            const BatchPeak::BatchPeakFitResult &peak_fit_results );

  /** Adds a per-source detection-limit rollup to the "Sources" of an activity/shielding fit result.
   A batch-only (e.g. FRMAC) convenience, so it is kept out of the GUI-shared
   `shield_src_fit_results_to_json(...)`.  Fields set on each source:
     - "HasDetectionLimit": bool; true when an Lc/MDA could be computed (guards the numeric fields).
     - "Lc_uCi"/"Lc_bq", "Mda_uCi"/"Mda_bq": the Currie decision threshold and detection limit,
       converted from counts to activity via the reporting line's `counts_per_bq`.
     - "RepresentativePeakEnergy": energy (keV) of the line the limit was taken from.
     - "UsedSubstitutePeak": bool; true when the exemplars designated line was not fit in this
       spectrum and a fallback line was used instead (see below).
     - "IsDetected": bool (source counts >= Currie decision threshold); report wording for detected
       vs. not (e.g. FRMAC "Approved" / "Less than Lc") is left to the template.
     - "DetectionLimitStatus": always-present short human-readable string explaining which line the
       limit came from, or why none could be computed - so a template can surface the reason rather
       than leaving a silent blank cell.

   `rep_energy_by_nuclide` maps a nuclide symbol to the energy (keV) of the representative gamma
   line chosen from the EXEMPLAR (its largest-amplitude peak used for the activity fit).  For each
   source this prefers that same line in `supp_info` (this spectrums peaks), so the reported Lc/MDA
   come from a consistent line across files.  When that designated line was not fit in a given
   spectrum, it falls back to the nuclides most prominent fitted line (flagged via
   "UsedSubstitutePeak") rather than reporting nothing; only when the nuclide has no usable fitted
   line at all is the limit omitted (with an explanatory "DetectionLimitStatus").
   */
  void add_exemplar_detection_limit_rollup_to_sources( nlohmann::json &data,
                  const std::vector<ShieldingSourceFitCalc::SupplementalPeakInfo> &supp_info,
                  const std::map<std::string,double> &rep_energy_by_nuclide );



  void write_json( const BatchPeak::BatchPeakFitOptions &options,
                  std::vector<std::string> &warnings,
                  const std::string &filename, nlohmann::json json_to_write );
}//namespace BatchInfoLog

#endif //BatchInfoLog_h
