#ifndef PeakSearchEval_h
#define PeakSearchEval_h
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
#include <functional>

#include "FitPeaksCorpusScore.h"
#include "FitPeaksCorpusReport.h"

namespace PeakFitUtils
{
  enum class CoarseResolutionType : int;
}

/** The `--search-only` mode of fit_peaks_corpus_eval: scores the automated peak search
 (`ExperimentalAutomatedPeakSearch::search_for_peaks`) and background-peak recovery against the
 GADRAS-inject truth, with sweeps of the search's significance cuts.
 */
namespace FitPeaksCorpus
{

struct SearchEvalSettings
{
  std::string out;
  PeakFitUtils::CoarseResolutionType det_type;
  /** Also search each problem's background spectrum and run background-peak recovery. */
  bool recovery = true;
  bool plot_data = true;
  /** Score the search's own (chi2) decision fits, not its peaks as reported (sparse ROIs refit by
   `ExperimentalAutomatedPeakSearch::refit_sparse_rois`). */
  bool decision_peaks = false;
  /** `ExperimentalAutomatedPeakSearch::SearchCuts` fields, as name=value. */
  std::vector<std::pair<std::string,std::string>> sets;
  std::string sweep_name;
  std::vector<std::string> sweep_values;

  double min_energy = 20.0;       // truth and search peaks below this are ignored
  double truth_min_z = 1.5;       // truth photopeaks below this are not scored as found/missed
  double moderate_z = 3.0;
  double strong_z = 8.0;
  double match_num_fwhm = 1.0;    // a search peak within this many truth FWHM ...
  double match_min_kev = 0.5;     // ... (at least this many keV) of a truth photopeak finds it
};//struct SearchEvalSettings


/** Applies one `ExperimentalAutomatedPeakSearch::SearchCuts` setting; returns false for an unknown name or
 bad value.  Not synchronized: call before any search or fit. */
bool apply_search_setting( const std::string &name, const std::string &value );

/** The current search cuts as name=value lines. */
std::string search_settings_text();

/** Runs the evaluation (or the sweep) and writes its outputs to `settings.out`.  `parallel(n, fcn)`
 must call `fcn(i)` for every i < n.  Returns the process exit code. */
int run_search_eval( const std::vector<CorpusProblem> &problems,
                     const SearchEvalSettings &settings,
                     RunMeta meta,
                     const std::function<void(size_t, const std::function<void(size_t)> &)> &parallel );

}//namespace FitPeaksCorpus

#endif //PeakSearchEval_h
