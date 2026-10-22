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

#ifndef FitPeaksCorpusReport_h
#define FitPeaksCorpusReport_h

/** Output writers for fit_peaks_corpus_eval: TSV/JSON summaries, the HTML review gallery, and
 result N42 files.  See FitPeaksCorpusScore.h for the data they consume. */

#include <string>
#include <vector>

#include "FitPeaksCorpusScore.h"

namespace FitPeaksCorpus
{

/** Provenance written to run_meta.txt and the gallery header. */
struct RunMeta
{
  std::string command_line;
  std::string started;
  std::string git_head;
  std::string git_dirty;
  std::string binary_sha256;
  std::string corpus_dir;
  std::string mode;
  std::string config_text;    // PeakFitForNuclideConfig::to_string()
  std::string weights_text;   // ScoreWeights::to_string()
};

/** Corpus-wide totals of the per-problem counts and costs. */
struct AggregateCounts
{
  size_t problems = 0;
  size_t mechanical_failures = 0;
  size_t nondeterministic = 0;
  size_t dev_check_failures = 0;
  size_t truth_scored = 0, matched = 0;
  size_t missed_strong = 0, missed_moderate = 0, missed_weak = 0;
  size_t extra_significant = 0, extra_weak = 0, ghosts = 0;
  size_t pairs = 0, share_disagree = 0;
  size_t rois_compared = 0, family_disagree = 0;
  size_t extent_sides = 0, extent_gt1fwhm_sides = 0;
  size_t extra_legit = 0, extra_xray = 0, extra_bkg_line = 0, truth_on_bkg_line = 0, truth_merged_away = 0;
  size_t truth_clusters = 0, truth_xray_clusters = 0, truth_fit_graded = 0, truth_ref_graded = 0, truth_fit_bad = 0, truth_ref_bad = 0;
  double truth_fit_abs_pull = 0.0, truth_ref_abs_pull = 0.0;
  double sum_raw_cost = 0.0;
  double mean_norm_cost = 0.0;
  double cost_missed = 0.0, cost_extra = 0.0, cost_share = 0.0, cost_family = 0.0;
  double cost_extent = 0.0, cost_area = 0.0, cost_mean = 0.0, cost_failure = 0.0;
  double cost_truth_area = 0.0;
  double wall_seconds = 0.0;
};

AggregateCounts aggregate( const std::vector<ProblemResult> &results );

/** One-line, human-readable aggregate (also used for sweep_summary.tsv). */
std::string summary_line( const AggregateCounts &counts );
std::string summary_tsv_header();
std::string summary_tsv_row( const std::string &label, const AggregateCounts &counts );

void write_run_meta( const std::string &path, const RunMeta &meta );
void write_per_problem_tsv( const std::string &path, const std::vector<ProblemResult> &results );
void write_per_peak_tsv( const std::string &path, const std::vector<ProblemResult> &results );

/** ROIs the planner produced, one row per ROI: id, lower, upper.  Lets a planning change be scored
 without running the solve (see PeakFitForNuclideConfig::stop_after_plan), which is milliseconds
 against ~10 s per problem for a full fit. */
void write_planned_rois_tsv( const std::string &path, const std::vector<ProblemResult> &results );
void write_per_roi_tsv( const std::string &path, const std::vector<ProblemResult> &results );
void write_per_pair_tsv( const std::string &path, const std::vector<ProblemResult> &results );
/** The truth-area grading: one row per merged inject-truth photopeak (fitted and reference pooled
 areas and pulls); empty unless problems carried an inject truth. */
void write_per_truth_area_tsv( const std::string &path, const std::vector<ProblemResult> &results );

/** Per-problem plot data for plot_fit_peaks_rois.py: the spectrum, the (live-time scaled)
 background, and for the fitted and reference sets every ROI with its per-channel continuum and
 per-peak Gaussian counts, plus the inject truth photopeaks and the automated-search peaks. */
void write_plot_data_json( const std::string &path, const ProblemResult &result );
void write_summary_json( const std::string &path, const std::vector<ProblemResult> &results,
                         const RunMeta &meta );

/** Self-contained HTML gallery: sortable summary table plus one lazily initialised D3 spectrum
 chart per problem (fitted peaks colored by verdict, a toggle to the reference peaks, missed
 reference peaks as reference lines) and the per-peak / per-ROI tables.  `resources_dir` must
 hold d3.v3.min.js, SpectrumChartD3.js and SpectrumChartD3.css. */
void write_gallery_html( const std::string &path, const std::string &resources_dir,
                         const std::string &title, const RunMeta &meta,
                         const std::vector<ProblemResult> &results, const size_t max_charts,
                         const bool include_background );

/** Writes the problem's spectra with the fitted peaks as the user peak set and the reference peaks
 as the automated-search peak set, so both can be inspected in InterSpec. */
void write_result_n42( const CorpusProblem &problem, const ProblemResult &result,
                       const std::string &path );

/** The "re-mined corpus rules" report: family counts, share/separate vs separation, extent
 quantiles, step-vs-linear detection-z quantiles, and quadratic ROI widths, computed with the
 canonical statistic over the reference sets. */
void write_truth_stats( const std::string &path, const std::vector<ProblemResult> &results );

}//namespace FitPeaksCorpus

#endif //FitPeaksCorpusReport_h
