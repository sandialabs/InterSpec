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

#include <set>
#include <memory>
#include <vector>
#include <utility>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"
#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

// Forward declarations
class SpecMeas;
class DetectorPeakResponse;


namespace ShieldingSourceFitCalc
{
  struct ModelFitResults;
  struct ShieldingSourceFitOptions;
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
  struct BatchPeakFitResult;
  struct BatchPeakFitOptions;
}

namespace SpecUtils
{
  class Measurement;
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
  
  /** Returns the default inja environment, with include directory, options, and callbacks set, as well as default templates loaded. */
  inja::Environment get_default_inja_env( const BatchPeak::BatchPeakFitOptions &options );
  
  /** Returns key-value pairs of of the file contents of the JS and CSS files needed for SpectrumChar.  Specifically returns:

   "D3_JS":                        contents of `InterSpec_resources/d3.v3.min.js`
   "SpectrumChart_JS":    contents of `InterSpec_resources/SpectrumChartD3.js`
   "SpectrumChart_CSS": contents of `InterSpec_resources/SpectrumChartD3.css`
   
   */
  std::vector<std::pair<std::string,std::string>> load_spectrum_chart_js_and_css();
  
  /** An enum to provide context of what default templates names "csv", "txt", and "html" refer to for `render_template(...)` */
  enum class TemplateRenderType : int
  {
    ActShieldIndividual,
    ActShieldSummary,
    PeakFitIndividual,
    PeakFitSummary
  };//enum class TemplateRenderType
  
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
  
  
  void write_json( const BatchPeak::BatchPeakFitOptions &options,
                  std::vector<std::string> &warnings,
                  const std::string &filename, nlohmann::json json_to_write );
}//namespace BatchInfoLog

#endif //BatchInfoLog_h
