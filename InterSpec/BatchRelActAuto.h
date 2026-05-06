#ifndef BatchRelActAuto_h
#define BatchRelActAuto_h
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
#include <optional>

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalcAuto.h"

class SpecMeas;
class DetectorPeakResponse;
namespace SpecUtils { class Measurement; }

/** Batch-mode driver for `RelActCalcAuto::solve`, parallel to `BatchActivity` and
 `BatchPeak`.  Reads an exemplar (an N42 with a stored `<RelActCalcAuto>` GUI
 state, or an external stand-alone XML config), runs `RelActCalcAuto::solve` on
 one or more input files, and renders per-file + summary reports via
 `RelActAutoReport`.
 */
namespace BatchRelActAuto
{
  /** Per-file options.  Does NOT inherit from `BatchPeak::BatchPeakFitOptions` —
   the peak-fit options (stat threshold, fit-all-peaks, refit-energy-cal, ...)
   don't apply to this analysis path; only the file/output and report bits are
   shared in spirit, declared here directly. */
  struct InterSpec_API Options
  {
    // -- Output bookkeeping (parallel to BatchPeakFitOptions output knobs) --
    std::string output_dir;
    bool overwrite_output_files = false;
    bool to_stdout = false;
    bool write_n42_with_results = false;
    bool create_json_output = false;

    // -- Reporting --
    /** One template per file.  Each entry is either a shorthand
     ("html", "txt", "json"), a filename to be resolved against
     `template_include_dir` (or the writable / bundled defaults), or an
     absolute path.  See `RelActAutoReport::render_template`. */
    std::vector<std::string> report_templates;
    std::vector<std::string> summary_report_templates;
    std::string template_include_dir;

    // -- Background plumbing --
    /** If non-empty, this background spectrum file overrides whatever the
     exemplar's stored `RelActAutoGuiState` says.  Forces
     `state.background_subtract = true`.  May be the same path as the input
     foreground file or the exemplar. */
    std::string background_subtract_file;
    std::set<int> background_subtract_samples;
    std::shared_ptr<SpecMeas> cached_background_subtract_spec;

    // -- Override slots (each empty / nullopt = "use exemplar value") --
    /** Detector response function override.  Required when the exemplar's
     state has any `RelEffEqnForm::FramPhysicalModel` rel-eff curve, OR when
     the FWHM-estimation method needs a DRF to seed from. */
    std::shared_ptr<const DetectorPeakResponse> drf_override;

    /** If set, fully replaces the `RelActAutoGuiState` that would otherwise be
     loaded from the exemplar's N42 — every field (rel-eff curves, ROIs, FWHM
     and skew options, energy-cal type, `background_subtract`, display ranges)
     comes from this state. */
    std::shared_ptr<const RelActCalcAuto::RelActAutoGuiState> state_override;

    std::optional<RelActCalcAuto::EnergyCalFitType> energy_cal_type;
    std::optional<RelActCalcAuto::FwhmForm> fwhm_form;
    std::optional<PeakDef::SkewType> skew_type;

    /** When set, forces `RelActAutoGuiState::background_subtract` regardless of
     what the exemplar / state-override said.  `false` here means
     "no background subtract"; `true` requires a usable background to be
     resolvable (either via `background_subtract_file` / cached spec, or via
     the exemplar's stored background filename). */
    std::optional<bool> background_subtract;
  };//struct Options


  enum class InterSpec_API ResultCode
  {
    Success,
    NoExemplar,
    CouldntOpenExemplar,
    ExemplarMissingRelActState,
    CouldntOpenStateOverride,
    CouldntOpenInputFile,
    CouldntOpenBackgroundFile,
    ForegroundSampleNumberUnderSpecified,
    BackgroundSampleNumberUnderSpecified,
    FwhmMethodNeedsDrfButNoneAvailable,
    SolveFailedToSetup,
    SolveFailedToSolve,
    SolveUserCanceled,
    SolveThrewException,
    UnknownStatus
  };//enum class ResultCode

  InterSpec_API const char *to_str( const ResultCode code );


  struct InterSpec_API Result
  {
    ResultCode m_result_code = ResultCode::UnknownStatus;
    std::string m_error_msg;
    std::vector<std::string> m_warnings;

    std::string m_filename;
    std::shared_ptr<const SpecMeas> m_foreground_file;
    std::set<int> m_foreground_sample_numbers;

    std::shared_ptr<const SpecMeas> m_exemplar_file;
    std::set<int> m_exemplar_sample_numbers;

    std::shared_ptr<const SpecMeas> m_background_file;
    std::set<int> m_background_sample_numbers;

    Options m_options;

    /** Will be default-constructed (Status::NotInitiated) on early-exit
     failures; on success, holds the full `RelActCalcAuto::solve` output. */
    RelActCalcAuto::RelActAutoSolution m_solution;
  };//struct Result


  struct InterSpec_API Summary
  {
    Options options;
    std::string exemplar_filename;
    std::shared_ptr<const SpecMeas> exemplar;
    std::set<int> exemplar_sample_nums;

    /** All four `file_*` vectors are kept the same length as the input file
     list. */
    std::vector<Result> file_results;
    std::vector<std::string> file_json;
    std::vector<std::vector<std::string>> file_reports;

    std::string summary_json;
    std::vector<std::string> summary_reports;
    std::vector<std::string> warnings;
  };//struct Summary


  /** Load a stand-alone `<RelActCalcAuto>` XML file (the exact format produced
   by `RelActAutoGui::guiStateToXml()` and read by
   `SpecMeas::getRelActAutoGuiState`) into an in-memory state struct.  Throws
   `std::runtime_error` on I/O / parse / schema failure. */
  InterSpec_API std::shared_ptr<RelActCalcAuto::RelActAutoGuiState>
    load_state_from_xml_file( const std::string &filename );


  /** Run a single-file analysis.  Mirrors
   `BatchActivity::fit_activities_in_file`.  Never throws on per-file errors —
   the failure is reported via `Result::m_result_code` / `m_error_msg`. */
  InterSpec_API Result run_on_file( const std::string &exemplar_filename,
                       std::set<int> exemplar_sample_nums,
                       std::shared_ptr<const SpecMeas> cached_exemplar,
                       const std::string &filename,
                       std::shared_ptr<SpecMeas> cached_file,
                       const Options &options );


  /** Run a multi-file batch.  Mirrors `BatchActivity::fit_activities_in_files`.
   Renders per-file and summary reports, writes them to disk if
   `options.output_dir` is non-empty, and (optionally) populates `summary` for
   programmatic consumers (the GUI builds its result dialog from this). */
  InterSpec_API void run_in_files( const std::string &exemplar_filename,
                      std::shared_ptr<const SpecMeas> cached_exemplar,
                      const std::set<int> &exemplar_sample_nums,
                      const std::vector<std::string> &files,
                      std::vector<std::shared_ptr<SpecMeas>> cached_files,
                      const Options &options,
                      Summary *summary = nullptr );

}//namespace BatchRelActAuto

#endif //BatchRelActAuto_h
