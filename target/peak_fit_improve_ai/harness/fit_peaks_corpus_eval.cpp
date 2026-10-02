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

/** fit_peaks_corpus_eval: scores FitPeaksForNuclides::fit_peaks_for_nuclides against manually fit
 reference N42 files (peak presence, ROI grouping, continuum family, ROI extent, parameters), with
 parameter sweeps, determinism repeats, an HTML review gallery and result N42s.

 Typical use (Release build; see the plan in ~/.claude/plans for the surrounding workflow):

   fit_peaks_corpus_eval --datadir=/path/to/InterSpec/data \
      --corpus=/path/to/manual_fits/Detective-X_300_seconds --out=/tmp/run1
   fit_peaks_corpus_eval ... --sweep=auto_keep_significance_z=2,2.5,3,4 --out=/tmp/sweep_keepz
   fit_peaks_corpus_eval ... --mode=norm-all --out=/tmp/norm
   fit_peaks_corpus_eval ... --dump-truth-stats --out=/tmp/truth

 Run with --help for the full option list.
 */
#include "InterSpec_config.h"

#include <map>
#include <set>
#include <mutex>
#include <cmath>
#include <ctime>
#include <atomic>
#include <chrono>
#include <thread>
#include <string>
#include <vector>
#include <cstdio>
#include <memory>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include <Wt/WFlags.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/FitPeaksForNuclides.h"

#include "FitPeaksCorpusScore.h"
#include "FitPeaksCorpusReport.h"

using namespace std;
using namespace FitPeaksCorpus;

namespace
{

struct Options
{
  string datadir;
  string corpus;
  string out;
  string problems;             // comma list of ids or id prefixes
  string problem_glob;         // simple wildcard on the id
  string mode = "source";      // source | norm-all | norm-each | all | sequential
  string engine = "legacy";    // legacy | policy
  string fit_options;          // comma list of FitSrcPeaksOptions names
  string det_type = "High";
  string background = "none";  // none | file
  string config_file;
  string resources_dir;
  string order;                // sequential mode source order
  string sources_override;     // force the requested sources for every selected problem
  string score_set = "observable";  // observable | fit | uncombined
  string corpus_format = "n42";     // n42 (InterSpec files with user peaks) | inject (GADRAS-inject truth csv + pcf)
  string truth_inject;              // inject directory whose <Source>_truth.csv files describe the n42 corpus spectra
  string debug_id;
  string sweep_name;
  vector<string> sweep_values;
  vector<pair<string,string>> sets;
  vector<pair<string,string>> weight_sets;
  int threads = 0;
  int solve_threads = 1;
  int repeat = 1;
  int gallery_max = 400;
  double fit_timeout = 900.0;   // seconds per fit before the cancel flag is raised
  bool html = true;
  bool plot_data = true;
  bool gallery_background = false;
  bool n42 = false;
  bool dump_truth_stats = false;
  bool list_only = false;
  bool force = false;
  bool help = false;
};//struct Options


void print_usage()
{
  cout <<
  "fit_peaks_corpus_eval - score fit_peaks_for_nuclides against manually fit reference N42 files\n\n"
  "Options (all --name=value or --name value):\n"
  "  --datadir DIR           InterSpec data directory (sandia.decay.xml); guessed if omitted\n"
  "  --corpus DIR            directory of reference N42 files (default: the Detective-X 300 s manual fits)\n"
  "  --out DIR               output directory (must not exist unless --force)\n"
  "  --problems a,b,...      only problems whose id equals or starts with one of these\n"
  "  --problem-glob PAT      only problems whose id matches the wildcard pattern (e.g. 'q1*')\n"
  "  --mode M                source | norm-all | norm-each | all | sequential (default source)\n"
  "  --engine E              legacy (default) | policy (adds DisableAutoInterfererFit)\n"
  "  --options a,b           FitSrcPeaksOptions by name, e.g. DoNotVaryEnergyCal,FitNormBkgrndPeaks\n"
  "  --sources a,b           force the requested sources for every selected problem\n"
  "  --order a,b             sequential mode: source order (default: the problem's sources)\n"
  "  --det-type T            PeakFitDetPrefs coarse type (default High)\n"
  "  --background none|file  supply the file's background spectrum to the fit (default none)\n"
  "  --set name=value        override a PeakFitForNuclideConfig field (repeatable)\n"
  "  --config FILE           name=value lines of config overrides\n"
  "  --sweep name=v1,v2,...  run once per value of one config field\n"
  "  --weight name=value     override a ScoreWeights field (repeatable)\n"
  "  --truth-min-z Z         don't-care threshold for reference peaks (default 1.5)\n"
  "  --score-set S           observable (default) | fit | uncombined\n"
  "  --truth-inject DIR      GADRAS-inject directory (<Det>/<Loc>/<time>_seconds) holding the truth of the\n"
  "                          n42 corpus spectra; unmatched fitted peaks on its photopeaks are 'legit', and\n"
  "                          fitted/reference areas are graded against it (default: the Detective-X\n"
  "                          Livermore 300 s directory when the default corpus is used; 'none' to skip)\n"
  "  --no-plot-data          skip the per-problem plot_data/<id>.json files (for plot_fit_peaks_rois.py)\n"
  "  --threads N             problems in parallel (default: hardware threads - 2)\n"
  "  --solve-threads N       threads inside each solve (default 1)\n"
  "  --repeat K              run K times and flag nondeterministic problems\n"
  "  --fit-timeout S         cancel a fit after S seconds and score it as a failure (default 900)\n"
  "  --n42                   write result N42 files (fitted peaks + reference peaks)\n"
  "  --no-html               skip the HTML gallery\n"
  "  --gallery-max N         charts in the gallery (default 400)\n"
  "  --gallery-background    also draw the background spectrum in the gallery charts (bigger file)\n"
  "  --resources-dir DIR     directory with d3.v3.min.js, SpectrumChartD3.js/.css\n"
  "  --dump-truth-stats      only write the reference-set statistics (no fitting)\n"
  "  --list                  list the selected problems and their sources, then exit\n"
  "  --debug ID              run one problem serially with the fitter's debug trace\n"
  "  --force                 allow writing into an existing output directory\n";
}


bool wildcard_match( const string &pattern, const string &text )
{
  size_t p = 0, t = 0, star = string::npos, match = 0;
  while( t < text.size() )
  {
    if( (p < pattern.size()) && ((pattern[p] == '?') || (pattern[p] == text[t])) )
    {
      ++p; ++t;
    }else if( (p < pattern.size()) && (pattern[p] == '*') )
    {
      star = p++;
      match = t;
    }else if( star != string::npos )
    {
      p = star + 1;
      t = ++match;
    }else
    {
      return false;
    }
  }
  while( (p < pattern.size()) && (pattern[p] == '*') )
    ++p;
  return p == pattern.size();
}


Options parse_options( const int argc, char **argv )
{
  Options opt;
  for( int i = 1; i < argc; ++i )
  {
    string arg = argv[i];
    if( !SpecUtils::istarts_with( arg, "--" ) )
      throw runtime_error( "Unexpected argument '" + arg + "'" );
    arg = arg.substr( 2 );
    string name = arg, value;
    bool have_value = false;
    const size_t eq = arg.find( '=' );
    if( eq != string::npos )
    {
      name = arg.substr( 0, eq );
      value = arg.substr( eq + 1 );
      have_value = true;
    }

    const auto need_value = [&]() -> string {
      if( have_value )
        return value;
      if( (i + 1) >= argc )
        throw runtime_error( "Option --" + name + " needs a value" );
      return argv[++i];
    };

    if( name == "help" )                  opt.help = true;
    else if( name == "datadir" )          opt.datadir = need_value();
    else if( name == "corpus" )           opt.corpus = need_value();
    else if( name == "out" )              opt.out = need_value();
    else if( name == "problems" )         opt.problems = need_value();
    else if( name == "problem-glob" )     opt.problem_glob = need_value();
    else if( name == "mode" )             opt.mode = need_value();
    else if( name == "engine" )           opt.engine = need_value();
    else if( name == "options" )          opt.fit_options = need_value();
    else if( name == "sources" )          opt.sources_override = need_value();
    else if( name == "order" )            opt.order = need_value();
    else if( name == "det-type" )         opt.det_type = need_value();
    else if( name == "background" )       opt.background = need_value();
    else if( name == "config" )           opt.config_file = need_value();
    else if( name == "resources-dir" )    opt.resources_dir = need_value();
    else if( name == "score-set" )        opt.score_set = need_value();
    else if( name == "corpus-format" )    opt.corpus_format = need_value();
    else if( name == "truth-inject" )     opt.truth_inject = need_value();
    else if( name == "no-plot-data" )     opt.plot_data = false;
    else if( name == "no-structure" )
    {
      for( const char *w : { "share_disagree", "family_step_disagree", "family_poly_disagree", "extent_per_fwhm" } )
        opt.weight_sets.emplace_back( w, "0" );
    }
    else if( name == "debug" )            opt.debug_id = need_value();
    else if( name == "threads" )          opt.threads = std::stoi( need_value() );
    else if( name == "solve-threads" )    opt.solve_threads = std::stoi( need_value() );
    else if( name == "repeat" )           opt.repeat = std::stoi( need_value() );
    else if( name == "gallery-max" )      opt.gallery_max = std::stoi( need_value() );
    else if( name == "fit-timeout" )      opt.fit_timeout = std::stod( need_value() );
    else if( name == "truth-min-z" )      opt.weight_sets.emplace_back( "truth_min_z", need_value() );
    else if( name == "n42" )              opt.n42 = true;
    else if( name == "gallery-background" ) opt.gallery_background = true;
    else if( name == "no-html" )          opt.html = false;
    else if( name == "html" )             opt.html = true;
    else if( name == "dump-truth-stats" ) opt.dump_truth_stats = true;
    else if( name == "list" )             opt.list_only = true;
    else if( name == "force" )            opt.force = true;
    else if( (name == "set") || (name == "weight") || (name == "sweep") )
    {
      const string v = need_value();
      const size_t veq = v.find( '=' );
      if( veq == string::npos )
        throw runtime_error( "--" + name + " expects name=value, got '" + v + "'" );
      const string field = v.substr( 0, veq );
      const string rhs = v.substr( veq + 1 );
      if( name == "set" )
        opt.sets.emplace_back( field, rhs );
      else if( name == "weight" )
        opt.weight_sets.emplace_back( field, rhs );
      else
      {
        if( !opt.sweep_name.empty() )
          throw runtime_error( "Only one --sweep is supported" );
        opt.sweep_name = field;
        SpecUtils::split( opt.sweep_values, rhs, "," );
      }
    }else
    {
      throw runtime_error( "Unknown option --" + name );
    }
  }//for( int i = 1; i < argc; ++i )

  return opt;
}//parse_options


bool set_weight( ScoreWeights &w, const string &name, const string &value )
{
  const map<string,double *> fields = {
    { "truth_min_z", &w.truth_min_z }, { "min_scored_energy", &w.min_scored_energy }, { "match_num_fwhm", &w.match_num_fwhm },
    { "match_min_kev", &w.match_min_kev }, { "weak_z", &w.weak_z }, { "strong_z", &w.strong_z },
    { "missed_weak", &w.missed_weak }, { "missed_moderate", &w.missed_moderate },
    { "missed_strong", &w.missed_strong }, { "ghost_z", &w.ghost_z }, { "extra_ghost", &w.extra_ghost },
    { "extra_weak", &w.extra_weak }, { "extra_significant", &w.extra_significant },
    { "share_max_sep_fwhm", &w.share_max_sep_fwhm }, { "share_disagree", &w.share_disagree },
    { "family_step_disagree", &w.family_step_disagree }, { "family_poly_disagree", &w.family_poly_disagree },
    { "extent_per_fwhm", &w.extent_per_fwhm }, { "extent_deadband_fwhm", &w.extent_deadband_fwhm },
    { "extent_cap_fwhm", &w.extent_cap_fwhm }, { "area_pull_weight", &w.area_pull_weight },
    { "area_rel_floor", &w.area_rel_floor }, { "area_pull_cap", &w.area_pull_cap },
    { "mean_offset_weight", &w.mean_offset_weight }, { "mean_offset_deadband_fwhm", &w.mean_offset_deadband_fwhm },
    { "mean_offset_cap_fwhm", &w.mean_offset_cap_fwhm }, { "failure", &w.failure },
    { "nondeterminism", &w.nondeterminism }, { "legit_truth_min_z", &w.legit_truth_min_z },
    { "truth_area_weight", &w.truth_area_weight }, { "truth_area_rel_floor", &w.truth_area_rel_floor },
    { "truth_area_cap", &w.truth_area_cap }, { "truth_area_min_z", &w.truth_area_min_z },
    { "legit_max_pull", &w.legit_max_pull }, { "merged_max_sep_fwhm", &w.merged_max_sep_fwhm }
  };
  const auto pos = fields.find( name );
  if( pos == end(fields) )
    return false;
  *pos->second = std::stod( value );
  return true;
}


string run_command( const string &command )
{
  string output;
  FILE *pipe = popen( command.c_str(), "r" );
  if( !pipe )
    return output;
  char buffer[512];
  while( fgets( buffer, sizeof(buffer), pipe ) )
    output += buffer;
  pclose( pipe );
  SpecUtils::trim( output );
  return output;
}


string now_string()
{
  const time_t now = time( nullptr );
  char buffer[64];
  strftime( buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", localtime( &now ) );
  return buffer;
}


Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> parse_fit_options( const string &csv )
{
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> flags;
  vector<string> names;
  SpecUtils::split( names, csv, "," );
  for( string &n : names )
  {
    SpecUtils::trim( n );
    if( n.empty() )
      continue;
    if( SpecUtils::iequals_ascii( n, "DoNotUseExistingRois" ) )
      flags |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois;
    else if( SpecUtils::iequals_ascii( n, "ExistingPeaksAsFreePeak" ) )
      flags |= FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak;
    else if( SpecUtils::iequals_ascii( n, "DoNotVaryEnergyCal" ) )
      flags |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotVaryEnergyCal;
    else if( SpecUtils::iequals_ascii( n, "DoNotRefineEnergyCal" ) )
      flags |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotRefineEnergyCal;
    else if( SpecUtils::iequals_ascii( n, "FitNormBkgrndPeaks" ) )
      flags |= FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks;
    else if( SpecUtils::iequals_ascii( n, "FitNormBkgrndPeaksDontUse" ) )
      flags |= FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse;
    else if( SpecUtils::iequals_ascii( n, "DisableAutoInterfererFit" ) )
      flags |= FitPeaksForNuclides::FitSrcPeaksOptions::DisableAutoInterfererFit;
    else
      throw runtime_error( "Unknown FitSrcPeaksOptions name '" + n + "'" );
  }
  return flags;
}


const char *status_name( const RelActCalcAuto::RelActAutoSolution::Status status )
{
  using Status = RelActCalcAuto::RelActAutoSolution::Status;
  switch( status )
  {
    case Status::Success:              return "Success";
    case Status::NotInitiated:         return "NotInitiated";
    case Status::FailedToSetupProblem: return "FailedToSetupProblem";
    case Status::FailToSolveProblem:   return "FailToSolveProblem";
    case Status::UserCanceled:         return "UserCanceled";
    case Status::UsableWithWarnings:   return "UsableWithWarnings";
  }
  return "Unknown";
}


/** One fit to perform: a problem, the sources, the option flags and the reference peaks to score against. */
struct Variant
{
  size_t problem = 0;
  string id;
  string mode;
  string sources_str;
  vector<RelActCalcAuto::SrcVariant> sources;
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options;
  vector<shared_ptr<const PeakDef>> truth;
  vector<shared_ptr<const PeakDef>> user_peaks;
};


void parallel_for( const size_t count, const int threads, const std::function<void(size_t)> &fcn )
{
  const int nthreads = std::max( 1, std::min( threads, static_cast<int>(count) ) );
  if( nthreads == 1 )
  {
    for( size_t i = 0; i < count; ++i )
      fcn( i );
    return;
  }

  std::atomic<size_t> next{ 0 };
  vector<std::thread> workers;
  for( int t = 0; t < nthreads; ++t )
  {
    workers.emplace_back( [&](){
      while( true )
      {
        const size_t i = next.fetch_add( 1 );
        if( i >= count )
          break;
        fcn( i );
      }
    } );
  }
  for( std::thread &w : workers )
    w.join();
}


vector<shared_ptr<const PeakDef>> to_shared( const vector<PeakDef> &peaks )
{
  vector<shared_ptr<const PeakDef>> answer;
  answer.reserve( peaks.size() );
  for( const PeakDef &p : peaks )
    answer.push_back( make_shared<PeakDef>( p ) );
  return answer;
}


struct FitContext
{
  const vector<CorpusProblem> *problems = nullptr;
  const vector<vector<shared_ptr<const PeakDef>>> *auto_peaks = nullptr;
  PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::High;
  bool supply_background = false;
  string score_set = "observable";
  double fit_timeout = 900.0;
  bool truth_is_inject = false;   // truth came from a GADRAS inject file, not from hand fits
  ScoreWeights weights;
};


/** Runs one fit and scores it. */
ProblemResult run_variant( const Variant &v, const FitPeaksForNuclides::PeakFitForNuclideConfig &config,
                           const FitContext &ctx )
{
  const CorpusProblem &problem = (*ctx.problems)[v.problem];
  ProblemResult r;
  r.id = v.id;
  r.mode = v.mode;
  r.sources = v.sources_str;
  r.spectrum_path = problem.path;
  r.foreground = problem.foreground;
  r.background = problem.background;
  r.truth_peaks = v.truth;

  auto prefs = make_shared<PeakFitDetPrefs>();
  prefs->m_det_type = ctx.det_type;

  const shared_ptr<const SpecUtils::Measurement> background = ctx.supply_background ? problem.background : nullptr;

  FitPeaksForNuclides::PeakFitResult fit;
  const auto start = chrono::steady_clock::now();

  // Watchdog: raise the cooperative cancel flag when the fit exceeds its time budget, so one
  // pathological solve cannot hold the whole corpus run.
  auto cancel_flag = make_shared<std::atomic_bool>( false );
  std::atomic<bool> fit_done{ false };
  std::atomic<bool> timed_out{ false };
  std::thread watchdog( [&](){
    const auto deadline = start + chrono::duration_cast<chrono::steady_clock::duration>( chrono::duration<double>( ctx.fit_timeout ) );
    while( !fit_done.load() )
    {
      if( chrono::steady_clock::now() > deadline )
      {
        timed_out.store( true );
        cancel_flag->store( true );
        break;
      }
      std::this_thread::sleep_for( chrono::milliseconds( 200 ) );
    }
  } );

  try
  {
    fit = FitPeaksForNuclides::fit_peaks_for_nuclides( (*ctx.auto_peaks)[v.problem], problem.foreground,
                                                       v.sources, v.user_peaks, background, nullptr,
                                                       v.options, config, prefs, cancel_flag );
    r.status = status_name( fit.status );
    for( const RelActCalcAuto::RoiRange &roi : fit.planned_rois )
      r.planned_rois.emplace_back( roi.lower_energy, roi.upper_energy );
    r.error = fit.error_message;
    r.warnings = fit.warnings;
    // The solver's own warnings do not reach PeakFitResult::warnings, but they are what turns a
    // fit into UsableWithWarnings (rank deficiency, an exhausted iteration budget); surface them
    // so a run can be asked which fits the data did not actually constrain.
    for( const std::string &w : fit.solution.m_warnings )
      r.warnings.push_back( "Solver: " + w );
    r.roi_plan_trace = fit.roi_plan_trace;
    r.mechanical_failure = !RelActCalcAuto::RelActAutoSolution::is_usable_status( fit.status );
  }catch( std::exception &e )
  {
    r.status = "Exception";
    r.error = e.what();
    r.mechanical_failure = true;
  }
  fit_done.store( true );
  watchdog.join();
  r.wall_seconds = chrono::duration<double>( chrono::steady_clock::now() - start ).count();
  if( timed_out.load() )
  {
    r.status = "Timeout";
    r.error = "Fit exceeded the " + std::to_string( ctx.fit_timeout ) + " s budget and was cancelled"
              + (r.error.empty() ? "" : ("; " + r.error));
    r.mechanical_failure = true;
  }

  for( const string &w : r.warnings )
  {
    if( SpecUtils::starts_with( w, "DevCheck:" ) )
      r.dev_check_failures += 1;
  }

  if( ctx.score_set == "fit" )
    r.fitted_peaks = to_shared( fit.fit_peaks );
  else if( ctx.score_set == "uncombined" )
    r.fitted_peaks = to_shared( fit.uncombined_fit_peaks );
  else
    r.fitted_peaks = to_shared( fit.observable_peaks );

  r.truth = make_peak_set( v.truth, problem.foreground, ctx.weights.truth_min_z );
  r.fitted = make_peak_set( r.fitted_peaks, problem.foreground, -1.0 );
  r.autosearch = make_peak_set( (*ctx.auto_peaks)[v.problem], problem.foreground, -1.0 );
  // The inject truth describes the source's own photopeaks; it applies to the source fits only
  // (not to NORM fits of the background files, whose truth would be the background photopeaks).
  r.inject_truth = (v.mode == "source") ? problem.inject_truth : nullptr;
  r.truth_is_inject = ctx.truth_is_inject;
  const std::vector<double> xrays = source_xray_energies( v.sources );
  const std::vector<double> src_lines = source_line_energies( v.sources );
  r.score = score_problem( r.truth, r.fitted, ctx.weights, r.inject_truth.get(),
                           xrays.empty() ? nullptr : &xrays, problem.foreground,
                           src_lines.empty() ? nullptr : &src_lines );
  if( r.mechanical_failure )
    r.score.cost_failure += ctx.weights.failure;

  return r;
}//run_variant


/** Applies a fit result to an accumulated user-peak list the way the GUI does. */
vector<shared_ptr<const PeakDef>> apply_fit_result( const vector<shared_ptr<const PeakDef>> &current,
                                                    const FitPeaksForNuclides::PeakFitResult &result )
{
  if( !RelActCalcAuto::RelActAutoSolution::is_usable_status( result.status ) || result.observable_peaks.empty() )
    return current;

  set<const PeakDef *> remove;
  for( const shared_ptr<const PeakDef> &p : result.original_peaks_to_remove )
    remove.insert( p.get() );

  vector<shared_ptr<const PeakDef>> updated;
  for( const shared_ptr<const PeakDef> &p : current )
  {
    if( p && !remove.count( p.get() ) )
      updated.push_back( p );
  }
  for( const PeakDef &p : result.observable_peaks )
    updated.push_back( make_shared<const PeakDef>( p ) );
  return updated;
}


/** Sequential-mode consistency: all-at-once vs one-at-a-time (given and reversed order). */
vector<ProblemResult> run_sequential( const size_t problem_index, const vector<RelActCalcAuto::SrcVariant> &order,
                                      const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options,
                                      const FitPeaksForNuclides::PeakFitForNuclideConfig &config,
                                      const FitContext &ctx )
{
  const CorpusProblem &problem = (*ctx.problems)[problem_index];
  vector<ProblemResult> results;

  Variant all;
  all.problem = problem_index;
  all.id = problem.id + "|all";
  all.mode = "sequential";
  all.sources = order;
  for( const RelActCalcAuto::SrcVariant &s : order )
    all.sources_str += (all.sources_str.empty() ? "" : ",") + RelActCalcAuto::to_name( s );
  all.options = options;
  all.truth = problem.truth_peaks;
  ProblemResult all_result = run_variant( all, config, ctx );
  const vector<shared_ptr<const PeakDef>> all_peaks = all_result.fitted_peaks;
  results.push_back( all_result );

  auto prefs = make_shared<PeakFitDetPrefs>();
  prefs->m_det_type = ctx.det_type;
  const shared_ptr<const SpecUtils::Measurement> background = ctx.supply_background ? problem.background : nullptr;

  for( const bool reversed : { false, true } )
  {
    vector<RelActCalcAuto::SrcVariant> seq = order;
    if( reversed )
      std::reverse( begin(seq), end(seq) );

    ProblemResult r;
    r.id = problem.id + (reversed ? "|seq_rev" : "|seq_fwd");
    r.mode = "sequential";
    r.foreground = problem.foreground;
    r.background = problem.background;
    r.truth_peaks = all_peaks;
    r.status = "Success";
    vector<shared_ptr<const PeakDef>> accumulated;
    const auto start = chrono::steady_clock::now();
    for( const RelActCalcAuto::SrcVariant &src : seq )
    {
      r.sources += (r.sources.empty() ? "" : ",") + RelActCalcAuto::to_name( src );
      try
      {
        const FitPeaksForNuclides::PeakFitResult fit = FitPeaksForNuclides::fit_peaks_for_nuclides(
          (*ctx.auto_peaks)[problem_index], problem.foreground, vector<RelActCalcAuto::SrcVariant>{ src },
          accumulated, background, nullptr, options, config, prefs );
        for( const string &w : fit.warnings )
        {
          r.warnings.push_back( RelActCalcAuto::to_name( src ) + ": " + w );
          if( SpecUtils::starts_with( w, "DevCheck:" ) )
            r.dev_check_failures += 1;
        }
        for( const string &line : fit.roi_plan_trace )
          r.roi_plan_trace.push_back( "[" + RelActCalcAuto::to_name( src ) + "] " + line );
        if( !RelActCalcAuto::RelActAutoSolution::is_usable_status( fit.status ) )
        {
          r.warnings.push_back( RelActCalcAuto::to_name( src ) + ": " + status_name( fit.status ) + " " + fit.error_message );
          r.mechanical_failure = true;
          r.status = status_name( fit.status );
        }
        accumulated = apply_fit_result( accumulated, fit );
      }catch( std::exception &e )
      {
        r.warnings.push_back( RelActCalcAuto::to_name( src ) + ": exception " + e.what() );
        r.mechanical_failure = true;
        r.status = "Exception";
      }
    }//for( const RelActCalcAuto::SrcVariant &src : seq )
    r.wall_seconds = chrono::duration<double>( chrono::steady_clock::now() - start ).count();

    r.fitted_peaks = accumulated;
    r.truth = make_peak_set( all_peaks, problem.foreground, ctx.weights.truth_min_z );
    r.fitted = make_peak_set( accumulated, problem.foreground, -1.0 );
    r.score = score_problem( r.truth, r.fitted, ctx.weights );
    if( r.mechanical_failure )
      r.score.cost_failure += ctx.weights.failure;
    results.push_back( r );
  }//for( const bool reversed : { false, true } )

  return results;
}//run_sequential


void ensure_directory( const string &dir, const bool force )
{
  if( SpecUtils::is_directory( dir ) )
  {
    if( !force && !SpecUtils::ls_files_in_directory( dir, "" ).empty() )
      throw runtime_error( "Output directory '" + dir + "' already has files (use --force to overwrite)" );
    return;
  }
  SpecUtils::create_directory( dir );
  if( !SpecUtils::is_directory( dir ) )
    throw runtime_error( "Could not create output directory '" + dir + "'" );
}


map<string,string> read_manifest( const string &corpus )
{
  map<string,string> answer;
  const string path = SpecUtils::append_path( corpus, "manifest.tsv" );
  ifstream in( path );
  if( !in.good() )
    return answer;
  string line;
  vector<string> header;
  int id_col = -1, src_col = -1;
  while( std::getline( in, line ) )
  {
    vector<string> cells;
    SpecUtils::split_no_delim_compress( cells, line, "\t" );
    if( header.empty() )
    {
      header = cells;
      for( size_t i = 0; i < header.size(); ++i )
      {
        if( (header[i] == "problem_id") || (header[i] == "id") )
          id_col = static_cast<int>( i );
        if( (header[i] == "requested_sources") || (header[i] == "sources") )
          src_col = static_cast<int>( i );
      }
      if( (id_col < 0) || (src_col < 0) )
      {
        cerr << "manifest.tsv has no problem_id/requested_sources columns; ignoring it" << endl;
        return answer;
      }
      continue;
    }
    if( (static_cast<int>(cells.size()) > id_col) && (static_cast<int>(cells.size()) > src_col) )
      answer[cells[id_col]] = cells[src_col];
  }
  return answer;
}


void write_outputs( const string &dir, const Options &opt, const vector<ProblemResult> &results,
                    const RunMeta &meta, const vector<CorpusProblem> &problems,
                    const map<string,size_t> &problem_index_by_result )
{
  ensure_directory( dir, true );
  write_run_meta( SpecUtils::append_path( dir, "run_meta.txt" ), meta );
  write_per_problem_tsv( SpecUtils::append_path( dir, "per_problem.tsv" ), results );
  write_per_peak_tsv( SpecUtils::append_path( dir, "per_peak.tsv" ), results );
  write_planned_rois_tsv( SpecUtils::append_path( dir, "planned_rois.tsv" ), results );
  write_per_roi_tsv( SpecUtils::append_path( dir, "per_roi.tsv" ), results );
  write_per_pair_tsv( SpecUtils::append_path( dir, "per_pair.tsv" ), results );
  write_per_truth_area_tsv( SpecUtils::append_path( dir, "per_truth_area.tsv" ), results );
  write_summary_json( SpecUtils::append_path( dir, "summary.json" ), results, meta );

  if( opt.plot_data )
  {
    const string plot_dir = SpecUtils::append_path( dir, "plot_data" );
    ensure_directory( plot_dir, true );
    for( const ProblemResult &r : results )
    {
      string name = r.id;
      for( char &c : name )
      {
        if( (c == '|') || (c == '/') || (c == ' ') )
          c = '_';
      }
      try
      {
        write_plot_data_json( SpecUtils::append_path( plot_dir, name + ".json" ), r );
      }catch( std::exception &e )
      {
        cerr << "Failed to write plot data for " << r.id << ": " << e.what() << endl;
      }
    }
  }//if( opt.plot_data )

  {
    ofstream timing( SpecUtils::append_path( dir, "timing.tsv" ) );
    ofstream exceptions( SpecUtils::append_path( dir, "exceptions.log" ) );
    ofstream devchecks( SpecUtils::append_path( dir, "dev_checks.log" ) );
    ofstream plan_trace( SpecUtils::append_path( dir, "roi_plan_trace.txt" ) );
    timing << "id\twall_seconds\tstatus\n";
    for( const ProblemResult &r : results )
    {
      timing << r.id << '\t' << r.wall_seconds << '\t' << r.status << '\n';
      if( r.mechanical_failure )
        exceptions << r.id << ": " << r.status << ": " << r.error << '\n';
      if( !r.roi_plan_trace.empty() )
      {
        plan_trace << "==== " << r.id << " (" << r.sources << ")\n";
        for( const string &line : r.roi_plan_trace )
          plan_trace << "  " << line << '\n';
      }
      for( const string &w : r.warnings )
      {
        if( SpecUtils::starts_with( w, "DevCheck:" ) )
          devchecks << r.id << ": " << w << '\n';
      }
    }
  }

  if( opt.html )
  {
    write_gallery_html( SpecUtils::append_path( dir, "gallery.html" ), opt.resources_dir,
                        "fit_peaks_corpus_eval " + meta.mode, meta, results,
                        static_cast<size_t>( std::max( 0, opt.gallery_max ) ), opt.gallery_background );
  }

  if( opt.n42 )
  {
    const string n42_dir = SpecUtils::append_path( dir, "n42" );
    ensure_directory( n42_dir, true );
    for( const ProblemResult &r : results )
    {
      const auto pos = problem_index_by_result.find( r.id );
      if( pos == end(problem_index_by_result) )
        continue;
      string name = r.id;
      for( char &c : name )
      {
        if( (c == '|') || (c == '/') || (c == ' ') )
          c = '_';
      }
      try
      {
        write_result_n42( problems[pos->second], r, SpecUtils::append_path( n42_dir, name + ".n42" ) );
      }catch( std::exception &e )
      {
        cerr << "Failed to write N42 for " << r.id << ": " << e.what() << endl;
      }
    }
  }
}//write_outputs

}//namespace


int main( int argc, char **argv )
{
  Options opt;
  try
  {
    opt = parse_options( argc, argv );
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
    print_usage();
    return 1;
  }

  if( opt.help )
  {
    print_usage();
    return 0;
  }

  try
  {
    if( opt.datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data", "../../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
        {
          opt.datadir = d;
          break;
        }
      }
    }
    if( opt.datadir.empty() || !SpecUtils::is_file( SpecUtils::append_path( opt.datadir, "sandia.decay.xml" ) ) )
      throw runtime_error( "--datadir must point at the InterSpec data directory" );

    if( opt.corpus.empty() )
      opt.corpus = "/Users/wcjohns/coding/InterSpec_peaks_for_source_opt/scratch/20260724_peaks_for_sources_manual_fits/Detective-X_300_seconds";
    if( !SpecUtils::is_directory( opt.corpus ) )
      throw runtime_error( "--corpus '" + opt.corpus + "' is not a directory" );

    if( opt.resources_dir.empty() )
      opt.resources_dir = SpecUtils::append_path( SpecUtils::append_path( SpecUtils::append_path( opt.datadir, ".." ), "external_libs" ), "SpecUtils/d3_resources" );

    if( opt.out.empty() && !opt.list_only )
      throw runtime_error( "--out is required" );

    if( opt.threads <= 0 )
      opt.threads = std::max( 1, static_cast<int>( std::thread::hardware_concurrency() ) - 2 );
    if( !opt.debug_id.empty() )
      opt.threads = 1;

    InterSpec::setStaticDataDirectory( opt.datadir );
    if( !DecayDataBaseServer::database() )
      throw runtime_error( "Could not load the decay database from " + opt.datadir );

    RelActCalc::set_max_solve_threads( static_cast<unsigned>( std::max( 1, opt.solve_threads ) ) );
    FitPeaksForNuclides::set_dev_checks_throw( true );
    FitPeaksForNuclides::set_debug_printout( !opt.debug_id.empty() );

    // ---- discover and load problems ----
    const bool inject_format = (opt.corpus_format == "inject");
    if( !inject_format && (opt.corpus_format != "n42") )
      throw runtime_error( "Unknown --corpus-format '" + opt.corpus_format + "' (n42 or inject)" );
    const map<string,string> manifest = inject_format ? map<string,string>{} : read_manifest( opt.corpus );
    vector<string> files = SpecUtils::ls_files_in_directory( opt.corpus, inject_format ? "_truth.csv" : ".n42" );
    std::sort( begin(files), end(files) );
    const auto id_of_file = [inject_format]( const string &f ) -> string {
      string id = problem_id_from_path( f );
      if( inject_format && SpecUtils::iends_with( id, "_truth" ) )
        id = id.substr( 0, id.size() - 6 );
      return id;
    };

    vector<string> id_filters;
    SpecUtils::split( id_filters, opt.problems, "," );
    const auto selected = [&]( const string &id ) -> bool {
      if( !opt.debug_id.empty() )
        return id == opt.debug_id;
      bool ok = id_filters.empty();
      for( const string &f : id_filters )
        ok = ok || (id == f) || SpecUtils::starts_with( id, f.c_str() );
      if( !ok )
        return false;
      if( !opt.problem_glob.empty() && !wildcard_match( opt.problem_glob, id ) )
        return false;
      return true;
    };

    vector<string> source_files, background_files;
    for( const string &f : files )
    {
      const string id = id_of_file( f );
      if( !inject_format && SpecUtils::istarts_with( id, "background" ) )
        background_files.push_back( f );
      else if( selected( id ) )
        source_files.push_back( f );
    }

    const bool want_source = (opt.mode == "source") || (opt.mode == "all") || (opt.mode == "sequential");
    const bool want_norm_all = (opt.mode == "norm-all") || (opt.mode == "all");
    const bool want_norm_each = (opt.mode == "norm-each") || (opt.mode == "all");
    if( !want_source && !want_norm_all && !want_norm_each )
      throw runtime_error( "Unknown --mode '" + opt.mode + "'" );

    vector<string> load_files;
    if( want_source )
      load_files = source_files;
    if( want_norm_all || want_norm_each )
      load_files.insert( end(load_files), begin(background_files), end(background_files) );
    if( load_files.empty() )
      throw runtime_error( "No problems selected" );

    vector<CorpusProblem> problems( load_files.size() );
    vector<string> load_errors( load_files.size() );
    parallel_for( load_files.size(), opt.threads, [&]( const size_t i ){
      try
      {
        problems[i] = inject_format ? load_inject_problem( load_files[i] )
                                    : load_problem( load_files[i], manifest );
      }catch( std::exception &e )
      {
        load_errors[i] = e.what();
      }
    } );

    {
      // Inject problems with an unknown source (pure beta emitters, unparsable names) are skipped
      // with a note; the hand-fit corpus must load completely.
      vector<CorpusProblem> loaded;
      for( size_t i = 0; i < problems.size(); ++i )
      {
        if( load_errors[i].empty() )
          loaded.push_back( std::move( problems[i] ) );
        else if( inject_format )
          cerr << "Skipping " << load_files[i] << ": " << load_errors[i] << endl;
        else
          throw runtime_error( "Problem load failed: " + load_errors[i] );
      }
      problems = std::move( loaded );
      if( problems.empty() )
        throw runtime_error( "No problems could be loaded" );
    }

    // Inject truth for the hand-fit corpus (the 300 s Detective-X hand fits are the Livermore inject
    // spectra; the truth files carry the expected photopeak areas that settle what the hand fits left
    // ambiguous).
    if( !inject_format )
    {
      const string default_corpus = "/Users/wcjohns/coding/InterSpec_peaks_for_source_opt/scratch/20260724_peaks_for_sources_manual_fits/Detective-X_300_seconds";
      const string default_truth = "/Users/wcjohns/coding/InterSpec_peak_fit_improve/peak_fit_accuracy_inject_compact/Detective-X/Livermore/300_seconds";
      if( opt.truth_inject.empty() && (opt.corpus == default_corpus) && SpecUtils::is_directory( default_truth ) )
        opt.truth_inject = default_truth;
      if( opt.truth_inject == "none" )
        opt.truth_inject.clear();
      if( !opt.truth_inject.empty() )
      {
        if( !SpecUtils::is_directory( opt.truth_inject ) )
          throw runtime_error( "--truth-inject '" + opt.truth_inject + "' is not a directory" );
        std::atomic<size_t> attached{ 0 }, unverified{ 0 };
        vector<string> attach_errors( problems.size() );
        parallel_for( problems.size(), opt.threads, [&]( const size_t i ){
          try
          {
            if( attach_inject_truth( problems[i], opt.truth_inject ) )
            {
              attached += 1;
              if( !problems[i].notes.empty() && SpecUtils::icontains( problems[i].notes.back(), "NOT VERIFIED" ) )
                unverified += 1;
            }
          }catch( std::exception &e )
          {
            attach_errors[i] = e.what();
          }
        } );
        for( size_t i = 0; i < problems.size(); ++i )
        {
          if( !attach_errors[i].empty() )
            cerr << "Inject truth for " << problems[i].id << ": " << attach_errors[i] << endl;
        }
        cout << "Inject truth attached to " << attached.load() << " of " << problems.size() << " problems ("
             << unverified.load() << " with a spectrum mismatch) from " << opt.truth_inject << endl;
      }
    }//if( !inject_format )

    const vector<RelActCalcAuto::SrcVariant> forced_sources = [&](){
      vector<RelActCalcAuto::SrcVariant> srcs;
      for( const string &n : split_source_list( opt.sources_override ) )
      {
        const RelActCalcAuto::SrcVariant s = RelActCalcAuto::source_from_string( n );
        if( RelActCalcAuto::is_null( s ) )
          throw runtime_error( "Unknown source '" + n + "' in --sources" );
        srcs.push_back( s );
      }
      return srcs;
    }();

    if( !forced_sources.empty() )
    {
      for( CorpusProblem &p : problems )
      {
        p.sources = forced_sources;
        p.requested_source_names = split_source_list( opt.sources_override );
        p.notes.push_back( "sources forced from --sources" );
      }
    }

    if( opt.list_only )
    {
      for( const CorpusProblem &p : problems )
      {
        cout << p.id << "\t" << SpecUtils::filename( p.path ) << "\tsources=";
        for( size_t i = 0; i < p.sources.size(); ++i )
          cout << (i ? "," : "") << RelActCalcAuto::to_name( p.sources[i] );
        cout << "\ttruth_peaks=" << p.truth_peaks.size() << "\tbackground=" << (p.background ? "yes" : "no")
             << "\tinject_truth=" << (p.inject_truth ? std::to_string( p.inject_truth->merged_signal.size() ) + " photopeaks" : string("none"));
        for( const string &n : p.notes )
          cout << "\t[" << n << "]";
        cout << "\n";
      }
      return 0;
    }

    // ---- weights, config ----
    ScoreWeights weights;
    for( const pair<string,string> &kv : opt.weight_sets )
    {
      if( !set_weight( weights, kv.first, kv.second ) )
        throw runtime_error( "Unknown --weight field '" + kv.first + "'" );
    }

    const PeakFitUtils::CoarseResolutionType det_type = PeakFitDetPrefs::coarse_res_from_str( opt.det_type );
    FitPeaksForNuclides::PeakFitForNuclideConfig config
      = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( det_type );

    if( !opt.config_file.empty() )
    {
      ifstream in( opt.config_file );
      if( !in.good() )
        throw runtime_error( "Could not open --config '" + opt.config_file + "'" );
      string line;
      while( std::getline( in, line ) )
      {
        const size_t hash = line.find( '#' );
        if( hash != string::npos )
          line = line.substr( 0, hash );
        SpecUtils::trim( line );
        if( line.empty() )
          continue;
        const size_t eq = line.find( '=' );
        if( eq == string::npos )
          throw runtime_error( "Bad config line '" + line + "'" );
        if( !config.set_field( line.substr( 0, eq ), line.substr( eq + 1 ) ) )
          throw runtime_error( "Bad config line '" + line + "' (unknown field or value)" );
      }
    }
    for( const pair<string,string> &kv : opt.sets )
    {
      if( !config.set_field( kv.first, kv.second ) )
        throw runtime_error( "--set " + kv.first + "=" + kv.second + ": unknown field or bad value" );
    }

    Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> base_options = parse_fit_options( opt.fit_options );
    if( opt.engine == "policy" )
      base_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DisableAutoInterfererFit;
    else if( opt.engine != "legacy" )
      throw runtime_error( "Unknown --engine '" + opt.engine + "'" );

    // ---- variants ----
    vector<Variant> variants;
    vector<size_t> sequential_problems;
    const vector<string> norm_nuclides = { "Th232", "U238", "Ra226", "K40", "U235" };
    for( size_t i = 0; i < problems.size(); ++i )
    {
      const CorpusProblem &p = problems[i];
      const bool is_background = SpecUtils::istarts_with( p.id, "background" );

      if( !is_background && want_source )
      {
        if( opt.mode == "sequential" )
        {
          sequential_problems.push_back( i );
          continue;
        }
        Variant v;
        v.problem = i;
        v.id = p.id;
        v.mode = "source";
        v.sources = p.sources;
        for( size_t k = 0; k < p.sources.size(); ++k )
          v.sources_str += (k ? "," : "") + RelActCalcAuto::to_name( p.sources[k] );
        v.options = base_options;
        v.truth = p.truth_peaks;
        variants.push_back( v );
      }

      if( is_background && want_norm_all )
      {
        Variant v;
        v.problem = i;
        v.id = p.id + "|norm-all";
        v.mode = "norm-all";
        v.sources_str = "NORM";
        v.options = base_options | FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks;
        v.truth = p.truth_peaks;
        variants.push_back( v );
      }

      if( is_background && want_norm_each )
      {
        for( const string &nuc : norm_nuclides )
        {
          Variant v;
          v.truth = filter_peaks_by_source( p.truth_peaks, { nuc } );
          if( v.truth.empty() )
            continue;
          v.problem = i;
          v.id = p.id + "|" + nuc;
          v.mode = "norm-each";
          v.sources_str = nuc;
          v.sources.push_back( RelActCalcAuto::source_from_string( nuc ) );
          v.options = base_options;
          variants.push_back( v );
        }
      }
    }//for( size_t i = 0; i < problems.size(); ++i )

    // ---- auto-search peaks, once per problem ----
    vector<vector<shared_ptr<const PeakDef>>> auto_peaks( problems.size() );
    {
      const auto start = chrono::steady_clock::now();
      parallel_for( problems.size(), opt.threads, [&]( const size_t i ){
        auto prefs = make_shared<PeakFitDetPrefs>();
        prefs->m_det_type = det_type;
        try
        {
          auto_peaks[i] = ExperimentalAutomatedPeakSearch::search_for_peaks( problems[i].foreground, nullptr, nullptr, true, prefs );
        }catch( std::exception &e )
        {
          cerr << "Auto-search failed for " << problems[i].id << ": " << e.what() << endl;
        }
      } );
      cout << "Auto-search of " << problems.size() << " spectra took "
           << chrono::duration<double>( chrono::steady_clock::now() - start ).count() << " s" << endl;
    }

    FitContext ctx;
    ctx.problems = &problems;
    ctx.auto_peaks = &auto_peaks;
    ctx.det_type = det_type;
    ctx.supply_background = (opt.background == "file");
    ctx.score_set = opt.score_set;
    ctx.fit_timeout = opt.fit_timeout;
    ctx.weights = weights;
    ctx.truth_is_inject = inject_format;

    map<string,size_t> problem_index_by_result;
    for( const Variant &v : variants )
      problem_index_by_result[v.id] = v.problem;
    for( const size_t pi : sequential_problems )
    {
      problem_index_by_result[problems[pi].id + "|all"] = pi;
      problem_index_by_result[problems[pi].id + "|seq_fwd"] = pi;
      problem_index_by_result[problems[pi].id + "|seq_rev"] = pi;
    }

    RunMeta meta;
    for( int i = 0; i < argc; ++i )
      meta.command_line += string( i ? " " : "" ) + argv[i];
    meta.started = now_string();
    const string repo_dir = SpecUtils::append_path( opt.datadir, ".." );
    meta.git_head = run_command( "git -C '" + repo_dir + "' rev-parse --short HEAD 2>/dev/null" );
    meta.git_dirty = run_command( "git -C '" + repo_dir + "' status --porcelain 2>/dev/null | head -c 1" ).empty() ? "clean" : "dirty";
    meta.binary_sha256 = run_command( string("shasum -a 256 '") + argv[0] + "' 2>/dev/null | cut -c1-16" );
    meta.corpus_dir = opt.corpus;
    meta.mode = opt.mode + " engine=" + opt.engine + " background=" + opt.background + " score_set=" + opt.score_set
                + (opt.truth_inject.empty() ? string() : (" truth_inject=" + opt.truth_inject));
    meta.weights_text = weights.to_string();

    ensure_directory( opt.out, opt.force );

    // ---- truth statistics only ----
    if( opt.dump_truth_stats )
    {
      vector<ProblemResult> truth_only;
      for( const Variant &v : variants )
      {
        ProblemResult r;
        r.id = v.id;
        r.mode = v.mode;
        r.sources = v.sources_str;
        r.status = "truth-only";
        r.foreground = problems[v.problem].foreground;
        r.truth_peaks = v.truth;
        r.truth = make_peak_set( v.truth, problems[v.problem].foreground, weights.truth_min_z );
        truth_only.push_back( r );
      }
      meta.config_text = config.to_string();
      write_truth_stats( SpecUtils::append_path( opt.out, "truth_stats.txt" ), truth_only );
      write_per_peak_tsv( SpecUtils::append_path( opt.out, "per_truth_peak.tsv" ), truth_only );
      write_per_roi_tsv( SpecUtils::append_path( opt.out, "per_truth_roi.tsv" ), truth_only );
      write_run_meta( SpecUtils::append_path( opt.out, "run_meta.txt" ), meta );
      cout << "Wrote reference statistics to " << opt.out << endl;
      return 0;
    }

    // ---- the run (with optional sweep and repeats) ----
    const auto run_all = [&]( const FitPeaksForNuclides::PeakFitForNuclideConfig &cfg, const string &progress_dir ) -> vector<ProblemResult> {
      vector<ProblemResult> results( variants.size() );
      std::mutex print_mutex;
      std::atomic<size_t> done{ 0 };
      ensure_directory( progress_dir, true );
      ofstream progress( SpecUtils::append_path( progress_dir, "progress.tsv" ), std::ios::app );
      parallel_for( variants.size(), opt.threads, [&]( const size_t i ){
        results[i] = run_variant( variants[i], cfg, ctx );
        const size_t n = ++done;
        std::lock_guard<std::mutex> lock( print_mutex );
        cout << "[" << n << "/" << variants.size() << "] " << results[i].id << " " << results[i].status
             << " raw=" << results[i].score.raw_cost() << " (" << results[i].wall_seconds << " s)" << endl;
        // Incremental record, so a killed or straggling run still leaves per-problem evidence.
        progress << results[i].id << '\t' << results[i].status << '\t' << results[i].score.raw_cost()
                 << '\t' << results[i].wall_seconds << '\t' << results[i].error << endl;
      } );

      vector<vector<ProblemResult>> seq_results( sequential_problems.size() );
      parallel_for( sequential_problems.size(), opt.threads, [&]( const size_t i ){
        seq_results[i] = run_sequential( sequential_problems[i], problems[sequential_problems[i]].sources,
                                         base_options, cfg, ctx );
        std::lock_guard<std::mutex> lock( print_mutex );
        for( const ProblemResult &r : seq_results[i] )
          cout << "[seq] " << r.id << " " << r.status << " raw=" << r.score.raw_cost() << endl;
      } );
      for( const vector<ProblemResult> &sr : seq_results )
        results.insert( end(results), begin(sr), end(sr) );

      if( opt.repeat > 1 )
      {
        for( int rep = 1; rep < opt.repeat; ++rep )
        {
          vector<ProblemResult> again( variants.size() );
          parallel_for( variants.size(), opt.threads, [&]( const size_t i ){
            again[i] = run_variant( variants[i], cfg, ctx );
          } );
          for( size_t i = 0; i < variants.size(); ++i )
          {
            if( fitted_fingerprint( again[i].fitted ) != fitted_fingerprint( results[i].fitted ) )
            {
              if( !results[i].nondeterministic )
              {
                results[i].nondeterministic = true;
                results[i].score.cost_failure += ctx.weights.nondeterminism;
                results[i].warnings.push_back( "Nondeterministic: repeat " + std::to_string( rep ) + " differs" );
              }
            }
          }
        }
      }//if( opt.repeat > 1 )

      return results;
    };//run_all

    if( !opt.sweep_name.empty() )
    {
      ofstream sweep_summary( SpecUtils::append_path( opt.out, "sweep_summary.tsv" ) );
      sweep_summary << summary_tsv_header() << "\n";
      for( const string &value : opt.sweep_values )
      {
        FitPeaksForNuclides::PeakFitForNuclideConfig cfg = config;
        if( !cfg.set_field( opt.sweep_name, value ) )
          throw runtime_error( "--sweep " + opt.sweep_name + "=" + value + ": unknown field or bad value" );
        cout << "==== sweep " << opt.sweep_name << "=" << value << " ====" << endl;
        const string sub_dir = SpecUtils::append_path( opt.out, "sweep_" + opt.sweep_name + "_" + value );
        const vector<ProblemResult> results = run_all( cfg, sub_dir );
        const AggregateCounts agg = aggregate( results );
        cout << summary_line( agg ) << endl;
        RunMeta sweep_meta = meta;
        sweep_meta.mode += " sweep " + opt.sweep_name + "=" + value;
        sweep_meta.config_text = cfg.to_string();
        write_outputs( sub_dir, opt, results, sweep_meta, problems, problem_index_by_result );
        sweep_summary << summary_tsv_row( opt.sweep_name + "=" + value, agg ) << "\n";
        sweep_summary.flush();
      }
      cout << "Sweep summary written to " << SpecUtils::append_path( opt.out, "sweep_summary.tsv" ) << endl;
      return 0;
    }

    const vector<ProblemResult> results = run_all( config, opt.out );
    const AggregateCounts agg = aggregate( results );
    meta.config_text = config.to_string();
    write_outputs( opt.out, opt, results, meta, problems, problem_index_by_result );
    cout << summary_line( agg ) << endl;
    cout << "Outputs written to " << opt.out << endl;
  }catch( std::exception &e )
  {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }

  return 0;
}//main
