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
#include <cmath>
#include <deque>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"

#include "PeakSearchEval.h"

using namespace std;

namespace FitPeaksCorpus
{
namespace
{

/** A truth photopeak of the whole spectrum: the source's and the background's rows, merged where the
 detector cannot resolve them. */
struct SpectrumTruth
{
  TruthPhotopeak peak;
  bool has_signal = false;   // at least one of the merged rows is a source line
  double z = 0.0;            // TruthPhotopeak::z_det()
};


vector<SpectrumTruth> spectrum_truth( const vector<TruthPhotopeak> &signal, const vector<TruthPhotopeak> &background )
{
  vector<TruthPhotopeak> rows = signal;
  rows.insert( end(rows), begin(background), end(background) );
  std::stable_sort( begin(rows), end(rows), []( const TruthPhotopeak &a, const TruthPhotopeak &b ){
    return a.energy < b.energy;
  } );

  vector<SpectrumTruth> answer;
  for( const TruthPhotopeak &m : merge_unresolved_photopeaks( rows ) )
  {
    SpectrumTruth t;
    t.peak = m;
    t.z = m.z_det();
    for( const TruthPhotopeak &s : signal )
      t.has_signal |= ((s.energy >= (m.energy_lo - 1.0E-6)) && (s.energy <= (m.energy_hi + 1.0E-6)));
    answer.push_back( t );
  }
  return answer;
}//spectrum_truth(...)


/** One search peak, as written to search_peaks.tsv. */
struct PeakRow
{
  string set;               // fg | bg | recovered
  double energy = 0.0, fwhm = 0.0, amplitude = 0.0, amplitude_uncert = 0.0;
  double marg_z = 0.0;      // amplitude / amplitude uncertainty, as the fit reports them
  double det_z = 0.0;       // PeakFitLM::peak_detection_significance(...)
  double chi2dof = 0.0;
  string verdict;           // signal | background | annihilation | escape | real_untabulated | unexplained | below_min_energy
  double truth_energy = 0.0, truth_z = 0.0;
  double reference_z = 0.0; // for a peak on no truth line: see reference_excess_z(...)
  // For a peak on no truth line, +-8 FWHM around it: channel lower energies, measured counts, the
  // reference spectrum, and the search fit's model of the peak's ROI (written to review_peaks.jsonl).
  vector<float> win_x, win_y, win_ref, win_model;
};


/** One truth photopeak, as written to search_truth.tsv. */
struct TruthRow
{
  string set;               // fg | bg (before recovery) | bg_recovered
  double energy = 0.0, fwhm = 0.0, area = 0.0, z = 0.0;
  bool has_signal = false;
  string truth_class;       // strong | moderate | weak
  bool found = false;
  double found_marg_z = 0.0, found_det_z = 0.0;
  double found_area = 0.0, found_area_uncert = 0.0;   // the matching search peak's
};


struct ClassCounts
{
  size_t strong = 0, strong_found = 0, moderate = 0, moderate_found = 0, weak = 0, weak_found = 0;
  void add( const TruthRow &t )
  {
    size_t &n = (t.truth_class == "strong") ? strong : ((t.truth_class == "moderate") ? moderate : weak);
    size_t &f = (t.truth_class == "strong") ? strong_found : ((t.truth_class == "moderate") ? moderate_found : weak_found);
    n += 1;
    f += t.found ? 1 : 0;
  }
  string str() const
  {
    ostringstream s;
    s << "strong=" << strong_found << "/" << strong << " moderate=" << moderate_found << "/" << moderate
      << " weak=" << weak_found << "/" << weak;
    return s.str();
  }
};//struct ClassCounts


struct VerdictCounts
{
  size_t total = 0, signal = 0, background = 0, annihilation = 0, escape = 0, real_untabulated = 0, unexplained = 0;
  void add( const PeakRow &p )
  {
    if( p.verdict == "below_min_energy" )
      return;
    total += 1;
    if( p.verdict == "signal" ) signal += 1;
    else if( p.verdict == "background" ) background += 1;
    else if( p.verdict == "annihilation" ) annihilation += 1;
    else if( p.verdict == "escape" ) escape += 1;
    else if( p.verdict == "real_untabulated" ) real_untabulated += 1;
    else unexplained += 1;
  }
};//struct VerdictCounts


struct ProblemOutcome
{
  string id;
  string error;
  string det_class;          // PeakFitUtils::coarse_resolution_from_peaks of the foreground search peaks
  double wall_seconds = 0.0;
  vector<PeakRow> peaks;
  vector<TruthRow> truth;
  vector<shared_ptr<const PeakDef>> fg_peaks;
  size_t num_recovered = 0;
};


/** A spectrum that shows the real features of a measurement better than the measurement itself: the
 inject PCF's noise-free record (3) for the foreground, its long background (record 2) for the background;
 `scale` takes it to the measurement's live time. */
struct ReferenceSpectrum
{
  std::shared_ptr<const SpecUtils::Measurement> spectrum;
  double scale = 1.0;
};


/** The net excess of the reference spectrum at `energy` - counts within +-1 FWHM over a linear continuum
 from sidebands 1.5-3 FWHM out on either side - in units of the measured spectrum's Poisson noise there.  A
 peak on no truth photopeak where this is at least 2 is a real feature the truth list does not hold (sum
 peaks, shield fluorescence, ...), not a false peak.  0 without a reference. */
double reference_excess_z( const ReferenceSpectrum &ref, const shared_ptr<const SpecUtils::Measurement> &data,
                           const double energy, const double fwhm )
{
  if( !ref.spectrum || !data || !(fwhm > 0.0) )
    return 0.0;
  const auto integral = [&ref]( const double lo, const double hi ) -> double {
    return ref.scale * ref.spectrum->gamma_integral( static_cast<float>(lo), static_cast<float>(hi) );
  };
  const double core = integral( energy - fwhm, energy + fwhm );
  const double left = integral( energy - 3.0*fwhm, energy - 1.5*fwhm );
  const double right = integral( energy + 1.5*fwhm, energy + 3.0*fwhm );
  const double continuum = (left + right) * (2.0*fwhm) / (3.0*fwhm);
  const double measured = data->gamma_integral( static_cast<float>(energy - fwhm), static_cast<float>(energy + fwhm) );
  return (core - continuum) / std::sqrt( std::max( measured, 1.0 ) );
}//reference_excess_z(...)


/** Fills `row.win_*` (see PeakRow) for `peak`, whose ROI is the peaks of `peaks` sharing its continuum. */
void fill_review_window( PeakRow &row, const PeakDef &peak, const vector<shared_ptr<const PeakDef>> &peaks,
                         const shared_ptr<const SpecUtils::Measurement> &data, const ReferenceSpectrum &reference )
{
  if( !data || !data->channel_energies() || !(peak.fwhm() > 0.0) )
    return;
  const size_t nchan = data->num_gamma_channels();
  const size_t lo = std::min( data->find_gamma_channel( static_cast<float>(peak.mean() - 8.0*peak.fwhm()) ), nchan - 1 );
  const size_t hi = std::min( data->find_gamma_channel( static_cast<float>(peak.mean() + 8.0*peak.fwhm()) ), nchan - 1 );
  if( hi <= lo )
    return;
  const size_t n = hi - lo + 1;
  const float * const energies = data->channel_energies()->data() + lo;

  vector<const PeakDef *> roi;
  for( const shared_ptr<const PeakDef> &p : peaks )
  {
    if( p && p->gausPeak() && (p->continuum() == peak.continuum()) )
      roi.push_back( p.get() );
  }
  vector<double> model( n, 0.0 );
  try
  {
    // Only within the ROI's range does its model mean anything.
    const shared_ptr<const PeakContinuum> cont = peak.continuum();
    const size_t roi_lo = data->find_gamma_channel( static_cast<float>( cont->lowerEnergy() ) );
    const size_t roi_hi = data->find_gamma_channel( static_cast<float>( cont->upperEnergy() ) );
    const size_t first = std::max( lo, roi_lo ), last = std::min( hi, roi_hi );
    if( last >= first )
    {
      cont->offset_integral( energies + (first - lo), model.data() + (first - lo), last - first + 1, data,
                             roi.data(), roi.size() );
      for( const PeakDef *p : roi )
        p->gauss_integral( energies + (first - lo), model.data() + (first - lo), last - first + 1 );
    }
  }catch( std::exception & )
  {
    std::fill( begin(model), end(model), 0.0 );
  }

  for( size_t i = 0; i < n; ++i )
  {
    row.win_x.push_back( energies[i] );
    row.win_y.push_back( data->gamma_channel_content( lo + i ) );
    row.win_ref.push_back( reference.spectrum
                             ? static_cast<float>( reference.scale * reference.spectrum->gamma_channel_content( lo + i ) )
                             : 0.0f );
    row.win_model.push_back( static_cast<float>( model[i] ) );
  }
}//fill_review_window(...)


/** Labels search peaks against `truth` and marks the truth photopeaks they find. */
void score_peaks( const string &set_name, const string &truth_set_name,
                  const vector<shared_ptr<const PeakDef>> &peaks,
                  const shared_ptr<const SpecUtils::Measurement> &data,
                  const vector<SpectrumTruth> &truth,
                  const ReferenceSpectrum &reference,
                  const SearchEvalSettings &settings, const bool chi2_fits,
                  vector<PeakRow> &peak_rows, vector<TruthRow> &truth_rows )
{
  const auto window = [&settings]( const double fwhm ) -> double {
    return std::max( settings.match_num_fwhm*fwhm, settings.match_min_kev );
  };
  const auto distance = []( const double energy, const TruthPhotopeak &t ) -> double {
    return std::max( 0.0, std::max( t.energy_lo - energy, energy - t.energy_hi ) );
  };

  vector<TruthRow> rows;
  for( const SpectrumTruth &t : truth )
  {
    TruthRow row;
    row.set = truth_set_name;
    row.energy = t.peak.energy;
    row.fwhm = t.peak.fwhm;
    row.area = t.peak.area;
    row.z = t.z;
    row.has_signal = t.has_signal;
    row.truth_class = (t.z >= settings.strong_z) ? "strong" : ((t.z >= settings.moderate_z) ? "moderate" : "weak");
    rows.push_back( row );
  }

  for( const shared_ptr<const PeakDef> &p : peaks )
  {
    if( !p || !p->gausPeak() )
      continue;

    PeakRow row;
    row.set = set_name;
    row.energy = p->mean();
    row.fwhm = p->fwhm();
    row.amplitude = p->amplitude();
    row.amplitude_uncert = p->amplitudeUncert();
    row.marg_z = (row.amplitude_uncert > 0.0) ? (row.amplitude / row.amplitude_uncert) : 0.0;
    row.det_z = PeakFitLM::peak_detection_significance( *p, peaks, data, chi2_fits );
    row.chi2dof = p->chi2dof();

    // The nearest truth photopeak in units of its match window.
    size_t best = truth.size();
    double best_ratio = 1.0;
    for( size_t i = 0; i < truth.size(); ++i )
    {
      const double ratio = distance( row.energy, truth[i].peak ) / window( truth[i].peak.fwhm );
      if( ratio <= best_ratio )
      {
        best = i;
        best_ratio = ratio;
      }
    }

    if( row.energy < settings.min_energy )
    {
      row.verdict = "below_min_energy";
    }else if( best < truth.size() )
    {
      row.verdict = truth[best].has_signal ? "signal" : "background";
      row.truth_energy = truth[best].peak.energy;
      row.truth_z = truth[best].z;
      TruthRow &t = rows[best];
      if( !t.found || (row.marg_z > t.found_marg_z) )
      {
        t.found_marg_z = row.marg_z;
        t.found_det_z = row.det_z;
        t.found_area = row.amplitude;
        t.found_area_uncert = row.amplitude_uncert;
      }
      t.found = true;
    }else if( std::fabs( row.energy - 510.999 ) <= window( row.fwhm ) )
    {
      row.verdict = "annihilation";
    }else
    {
      // Single and double escape peaks of a significant line.
      bool escape = false;
      for( const SpectrumTruth &t : truth )
      {
        if( (t.peak.energy < 1100.0) || (t.z < settings.moderate_z) )
          continue;
        escape |= (std::fabs( row.energy - (t.peak.energy - 510.999) ) <= window( t.peak.fwhm ))
                  || (std::fabs( row.energy - (t.peak.energy - 1021.998) ) <= window( t.peak.fwhm ));
      }
      if( escape )
      {
        row.verdict = "escape";
      }else
      {
        row.reference_z = reference_excess_z( reference, data, row.energy, row.fwhm );
        row.verdict = (row.reference_z >= 2.0) ? "real_untabulated" : "unexplained";
        fill_review_window( row, *p, peaks, data, reference );
      }
    }

    peak_rows.push_back( row );
  }//for( const shared_ptr<const PeakDef> &p : peaks )

  for( const TruthRow &t : rows )
  {
    if( (t.energy >= settings.min_energy) && (t.z >= settings.truth_min_z) )
      truth_rows.push_back( t );
  }
}//score_peaks(...)


bool apply_setting( const string &name, const string &value )
{
  double v = 0.0;
  try
  {
    size_t pos = 0;
    v = std::stod( value, &pos );
    if( pos != value.size() )
      return false;
  }catch( std::exception & )
  {
    return false;
  }

  ExperimentalAutomatedPeakSearch::SearchCuts &cuts = ExperimentalAutomatedPeakSearch::search_cuts();
  const map<string,double *> doubles = {
    { "highres_min_nsigma", &cuts.highres_min_nsigma },
    { "highres_low_stat_min_nsigma", &cuts.highres_low_stat_min_nsigma },
    { "highres_med_stat_min_nsigma", &cuts.highres_med_stat_min_nsigma },
    { "highres_chi2dof_cut_max_nsigma", &cuts.highres_chi2dof_cut_max_nsigma },
    { "lowres_single_min_nsigma", &cuts.lowres_single_min_nsigma },
    { "lowres_multi_min_nsigma", &cuts.lowres_multi_min_nsigma },
    { "highres_multi_min_nsigma", &cuts.highres_multi_min_nsigma },
    { "background_recovery_min_nsigma", &cuts.background_recovery_min_nsigma }
  };
  const auto d = doubles.find( name );
  if( d != end(doubles) )
  {
    *d->second = v;
    return true;
  }
  if( name == "detection_z_chi2_weights" )
  {
    cuts.detection_z_chi2_weights = (v != 0.0);
    return true;
  }
  return false;
}//apply_setting(...)


string settings_text()
{
  const ExperimentalAutomatedPeakSearch::SearchCuts &c = ExperimentalAutomatedPeakSearch::search_cuts();
  ostringstream out;
  out << "detection_z_chi2_weights=" << c.detection_z_chi2_weights
      << "\nhighres_min_nsigma=" << c.highres_min_nsigma
      << "\nhighres_low_stat_min_nsigma=" << c.highres_low_stat_min_nsigma
      << "\nhighres_med_stat_min_nsigma=" << c.highres_med_stat_min_nsigma
      << "\nhighres_chi2dof_cut_max_nsigma=" << c.highres_chi2dof_cut_max_nsigma
      << "\nlowres_single_min_nsigma=" << c.lowres_single_min_nsigma
      << "\nlowres_multi_min_nsigma=" << c.lowres_multi_min_nsigma
      << "\nhighres_multi_min_nsigma=" << c.highres_multi_min_nsigma
      << "\nbackground_recovery_min_nsigma=" << c.background_recovery_min_nsigma << "\n";
  return out.str();
}//settings_text()


ProblemOutcome run_problem( const CorpusProblem &problem, const SearchEvalSettings &settings )
{
  ProblemOutcome outcome;
  outcome.id = problem.id;
  const auto start = chrono::steady_clock::now();

  auto prefs = make_shared<PeakFitDetPrefs>();
  prefs->m_det_type = settings.det_type;

  try
  {
    if( !problem.inject_truth )
      throw runtime_error( "no inject truth" );

    const vector<shared_ptr<const PeakDef>> fg_decision_peaks
              = ExperimentalAutomatedPeakSearch::search_for_peaks( problem.foreground, nullptr, nullptr, true, prefs );
    outcome.fg_peaks = fg_decision_peaks;
    // Scored as reported to a user (the GUI's search and BatchPeak refit the search's sparse ROIs).
    if( !settings.decision_peaks )
    {
      outcome.fg_peaks = ExperimentalAutomatedPeakSearch::refit_sparse_rois( outcome.fg_peaks, problem.foreground,
                                                                              settings.det_type, {} );
    }
    outcome.det_class = PeakFitDetPrefs::to_str( PeakFitUtils::coarse_resolution_from_peaks( outcome.fg_peaks, problem.foreground ) );

    // The PCF's noise-free foreground (record 3) and long background (record 2).
    ReferenceSpectrum fg_reference, bg_reference;
    {
      SpecUtils::SpecFile pcf;
      if( pcf.load_file( problem.path, SpecUtils::ParserType::Auto ) && (pcf.num_measurements() >= 4) )
      {
        const shared_ptr<const SpecUtils::Measurement> noise_free = pcf.measurement_at_index( 3 );
        const shared_ptr<const SpecUtils::Measurement> long_bg = pcf.measurement_at_index( 2 );
        if( noise_free && (noise_free->num_gamma_channels() == problem.foreground->num_gamma_channels()) )
          fg_reference.spectrum = noise_free;
        if( long_bg && problem.background && (long_bg->live_time() > 0.0f)
            && (long_bg->num_gamma_channels() == problem.background->num_gamma_channels()) )
        {
          bg_reference.spectrum = long_bg;
          bg_reference.scale = problem.background->live_time() / long_bg->live_time();
        }
      }
    }

    const vector<SpectrumTruth> fg_truth = spectrum_truth( problem.inject_truth->signal, problem.inject_truth->background );
    score_peaks( "fg", "fg", outcome.fg_peaks, problem.foreground, fg_truth, fg_reference, settings,
                 settings.decision_peaks, outcome.peaks, outcome.truth );

    if( settings.recovery && problem.background )
    {
      // The background record has the foreground's duration, so the truth's background photopeaks
      // describe it.
      const vector<SpectrumTruth> bg_truth = spectrum_truth( {}, problem.inject_truth->background );
      const vector<shared_ptr<const PeakDef>> bg_peaks
                = ExperimentalAutomatedPeakSearch::search_for_peaks( problem.background, nullptr, nullptr, true, prefs );
      score_peaks( "bg", "bg", bg_peaks, problem.background, bg_truth, bg_reference, settings, true,
                   outcome.peaks, outcome.truth );

      // Run, as the GUI runs it, on the searches' decision fits (its hint peaks).
      const auto bg_deque = make_shared<const deque<shared_ptr<const PeakDef>>>( begin(bg_peaks), end(bg_peaks) );
      const shared_ptr<const deque<shared_ptr<const PeakDef>>> recovered
        = ExperimentalAutomatedPeakSearch::recover_background_peaks_under_foreground( fg_decision_peaks,
                                                    problem.background, bg_deque, nullptr, prefs );
      const vector<shared_ptr<const PeakDef>> recovered_peaks = recovered
          ? vector<shared_ptr<const PeakDef>>( begin(*recovered), end(*recovered) ) : bg_peaks;

      // Only the peaks recovery added (or replaced) get the "recovered" label.
      vector<PeakRow> recovered_rows;
      vector<TruthRow> recovered_truth;
      score_peaks( "recovered", "bg_recovered", recovered_peaks, problem.background, bg_truth, bg_reference, settings,
                   true, recovered_rows, recovered_truth );
      for( size_t i = 0; i < recovered_peaks.size(); ++i )
      {
        const bool is_new = std::none_of( begin(bg_peaks), end(bg_peaks), [&]( const shared_ptr<const PeakDef> &p ){
          return p.get() == recovered_peaks[i].get();
        } );
        if( !is_new )
          continue;
        outcome.num_recovered += 1;
        const auto row = std::find_if( begin(recovered_rows), end(recovered_rows), [&]( const PeakRow &r ){
          return r.energy == recovered_peaks[i]->mean();
        } );
        if( row != end(recovered_rows) )
          outcome.peaks.push_back( *row );
      }
      outcome.truth.insert( end(outcome.truth), begin(recovered_truth), end(recovered_truth) );
    }//if( settings.recovery && problem.background )
  }catch( std::exception &e )
  {
    outcome.error = e.what();
  }

  outcome.wall_seconds = chrono::duration<double>( chrono::steady_clock::now() - start ).count();
  return outcome;
}//run_problem(...)


struct SearchAggregate
{
  ClassCounts fg, bg, bg_recovered;
  VerdictCounts fg_peaks, bg_peaks, recovered_peaks;
  size_t problems = 0, failures = 0;
  double wall_seconds = 0.0;
  map<string,size_t> det_classes;

  void add( const ProblemOutcome &o )
  {
    problems += 1;
    failures += o.error.empty() ? 0 : 1;
    wall_seconds += o.wall_seconds;
    det_classes[o.det_class] += 1;
    for( const TruthRow &t : o.truth )
    {
      if( t.set == "fg" ) fg.add( t );
      else if( t.set == "bg" ) bg.add( t );
      else bg_recovered.add( t );
    }
    for( const PeakRow &p : o.peaks )
    {
      if( p.set == "fg" ) fg_peaks.add( p );
      else if( p.set == "bg" ) bg_peaks.add( p );
      else recovered_peaks.add( p );
    }
  }

  string summary_line() const
  {
    ostringstream s;
    s << "problems=" << problems << " fail=" << failures << " fg: " << fg.str()
      << " peaks=" << fg_peaks.total << " unexplained=" << fg_peaks.unexplained
      << " real_untabulated=" << fg_peaks.real_untabulated << " annih=" << fg_peaks.annihilation
      << " escape=" << fg_peaks.escape << " | bg: " << bg.str() << " unexplained=" << bg_peaks.unexplained
      << " | recovered=" << recovered_peaks.total << " (unexplained " << recovered_peaks.unexplained << "): "
      << bg_recovered.str() << " | det_class:";
    for( const auto &dc : det_classes )
      s << " " << (dc.first.empty() ? string("none") : dc.first) << "=" << dc.second;
    s << " | cpu=" << std::fixed << std::setprecision( 1 ) << wall_seconds << "s";
    return s.str();
  }

  static string tsv_header()
  {
    return "label\tproblems\tfail\tfg_strong\tfg_strong_found\tfg_moderate\tfg_moderate_found\tfg_weak\tfg_weak_found"
           "\tfg_peaks\tfg_unexplained\tfg_annihilation\tfg_escape\tbg_strong_found\tbg_moderate_found\tbg_weak_found"
           "\tbg_unexplained\trecovered\trecovered_unexplained\trec_strong_found\trec_moderate_found\trec_weak_found";
  }

  string tsv_row( const string &label ) const
  {
    ostringstream s;
    s << label << '\t' << problems << '\t' << failures << '\t' << fg.strong << '\t' << fg.strong_found
      << '\t' << fg.moderate << '\t' << fg.moderate_found << '\t' << fg.weak << '\t' << fg.weak_found
      << '\t' << fg_peaks.total << '\t' << fg_peaks.unexplained << '\t' << fg_peaks.annihilation
      << '\t' << fg_peaks.escape << '\t' << bg.strong_found << '\t' << bg.moderate_found << '\t' << bg.weak_found
      << '\t' << bg_peaks.unexplained << '\t' << recovered_peaks.total << '\t' << recovered_peaks.unexplained
      << '\t' << bg_recovered.strong_found << '\t' << bg_recovered.moderate_found << '\t' << bg_recovered.weak_found;
    return s.str();
  }
};//struct SearchAggregate


void write_outputs( const string &dir, const vector<CorpusProblem> &problems, const vector<ProblemOutcome> &outcomes,
                    const SearchEvalSettings &settings, RunMeta meta, const SearchAggregate &agg )
{
  if( !SpecUtils::is_directory( dir ) && SpecUtils::create_directory( dir ) != 1 )
    throw runtime_error( "Could not create directory '" + dir + "'" );

  const auto open = [&dir]( const string &name ) -> ofstream {
    ofstream out( SpecUtils::append_path( dir, name ) );
    if( !out.good() )
      throw runtime_error( "Could not write '" + name + "' in '" + dir + "'" );
    out << std::setprecision( 8 );
    return out;
  };

  {
    ofstream out = open( "search_peaks.tsv" );
    out << "problem\tset\tenergy\tfwhm\tamplitude\tamplitude_uncert\tmarg_z\tdet_z\tchi2dof\tverdict\ttruth_energy\ttruth_z"
           "\treference_z\n";
    for( const ProblemOutcome &o : outcomes )
    {
      for( const PeakRow &p : o.peaks )
        out << o.id << '\t' << p.set << '\t' << p.energy << '\t' << p.fwhm << '\t' << p.amplitude << '\t'
            << p.amplitude_uncert << '\t' << p.marg_z << '\t' << p.det_z << '\t' << p.chi2dof << '\t' << p.verdict
            << '\t' << p.truth_energy << '\t' << p.truth_z << '\t' << p.reference_z << '\n';
    }
  }

  {
    // One JSON object per peak on no truth line, for tools/search_review.py.
    ofstream out = open( "review_peaks.jsonl" );
    const auto array = []( const vector<float> &v ) -> string {
      string a = "[";
      for( size_t i = 0; i < v.size(); ++i )
        a += (i ? "," : "") + SpecUtils::printCompact( v[i], 6 );
      return a + "]";
    };
    for( const ProblemOutcome &o : outcomes )
    {
      for( const PeakRow &p : o.peaks )
      {
        if( p.win_x.empty() )
          continue;
        out << "{\"problem\":\"" << o.id << "\",\"set\":\"" << p.set << "\",\"energy\":" << p.energy
            << ",\"fwhm\":" << p.fwhm << ",\"amplitude\":" << p.amplitude << ",\"marg_z\":" << p.marg_z
            << ",\"det_z\":" << p.det_z << ",\"verdict\":\"" << p.verdict << "\",\"reference_z\":" << p.reference_z
            << ",\"x\":" << array( p.win_x ) << ",\"y\":" << array( p.win_y ) << ",\"ref\":" << array( p.win_ref )
            << ",\"model\":" << array( p.win_model ) << "}\n";
      }
    }
  }

  {
    ofstream out = open( "search_truth.tsv" );
    out << "problem\tset\tenergy\tfwhm\tarea\tz\tclass\thas_signal\tfound\tfound_marg_z\tfound_det_z"
           "\tfound_area\tfound_area_uncert\n";
    for( const ProblemOutcome &o : outcomes )
    {
      for( const TruthRow &t : o.truth )
        out << o.id << '\t' << t.set << '\t' << t.energy << '\t' << t.fwhm << '\t' << t.area << '\t' << t.z << '\t'
            << t.truth_class << '\t' << (t.has_signal ? 1 : 0) << '\t' << (t.found ? 1 : 0) << '\t' << t.found_marg_z
            << '\t' << t.found_det_z << '\t' << t.found_area << '\t' << t.found_area_uncert << '\n';
    }
  }

  {
    ofstream out = open( "per_problem.tsv" );
    out << "problem\tdet_class\tfg_peaks\tfg_unexplained\trecovered\tcpu_s\terror\n";
    for( const ProblemOutcome &o : outcomes )
    {
      VerdictCounts v;
      for( const PeakRow &p : o.peaks )
      {
        if( p.set == "fg" )
          v.add( p );
      }
      out << o.id << '\t' << o.det_class << '\t' << v.total << '\t' << v.unexplained << '\t' << o.num_recovered
          << '\t' << o.wall_seconds << '\t' << o.error << '\n';
    }
  }

  {
    ofstream out = open( "summary.tsv" );
    out << SearchAggregate::tsv_header() << "\n" << agg.tsv_row( "all" ) << "\n";
  }

  meta.config_text = settings_text();
  write_run_meta( SpecUtils::append_path( dir, "run_meta.txt" ), meta );

  if( settings.plot_data )
  {
    const string plot_dir = SpecUtils::append_path( dir, "plot_data" );
    if( !SpecUtils::is_directory( plot_dir ) && SpecUtils::create_directory( plot_dir ) != 1 )
      throw runtime_error( "Could not create '" + plot_dir + "'" );
    for( size_t i = 0; i < outcomes.size(); ++i )
    {
      ProblemResult r;
      r.id = outcomes[i].id;
      r.status = outcomes[i].error.empty() ? "search" : "search failed";
      r.foreground = problems[i].foreground;
      r.background = problems[i].background;
      r.truth = make_peak_set( problems[i].truth_peaks, problems[i].foreground, settings.truth_min_z );
      r.fitted = make_peak_set( outcomes[i].fg_peaks, problems[i].foreground, -1.0 );
      r.truth_is_inject = true;
      write_plot_data_json( SpecUtils::append_path( plot_dir, r.id + ".json" ), r );
    }
  }//if( settings.plot_data )
}//write_outputs(...)

}//namespace


bool apply_search_setting( const string &name, const string &value )
{
  return apply_setting( name, value );
}


string search_settings_text()
{
  return settings_text();
}


int run_search_eval( const vector<CorpusProblem> &problems,
                     const SearchEvalSettings &settings,
                     RunMeta meta,
                     const std::function<void(size_t, const std::function<void(size_t)> &)> &parallel )
{
  for( const pair<string,string> &kv : settings.sets )
  {
    if( !apply_setting( kv.first, kv.second ) )
      throw runtime_error( "--search-set " + kv.first + "=" + kv.second + ": unknown name or bad value" );
  }

  const auto run_all = [&]() -> pair<vector<ProblemOutcome>,SearchAggregate> {
    vector<ProblemOutcome> outcomes( problems.size() );
    parallel( problems.size(), [&]( const size_t i ){
      outcomes[i] = run_problem( problems[i], settings );
    } );
    SearchAggregate agg;
    for( const ProblemOutcome &o : outcomes )
    {
      agg.add( o );
      if( !o.error.empty() )
        cerr << "Search of " << o.id << " failed: " << o.error << endl;
    }
    return { std::move( outcomes ), agg };
  };

  meta.mode += " search-only";

  if( settings.sweep_name.empty() )
  {
    const pair<vector<ProblemOutcome>,SearchAggregate> result = run_all();
    write_outputs( settings.out, problems, result.first, settings, meta, result.second );
    cout << result.second.summary_line() << endl;
    cout << "Outputs written to " << settings.out << endl;
    return 0;
  }

  if( !SpecUtils::is_directory( settings.out ) && SpecUtils::create_directory( settings.out ) != 1 )
    throw runtime_error( "Could not create '" + settings.out + "'" );
  ofstream summary( SpecUtils::append_path( settings.out, "sweep_summary.tsv" ) );
  summary << SearchAggregate::tsv_header() << "\n";
  for( const string &value : settings.sweep_values )
  {
    if( !apply_setting( settings.sweep_name, value ) )
      throw runtime_error( "--search-sweep " + settings.sweep_name + "=" + value + ": unknown name or bad value" );
    const pair<vector<ProblemOutcome>,SearchAggregate> result = run_all();
    const string label = settings.sweep_name + "=" + value;
    SearchEvalSettings sub_settings = settings;
    sub_settings.plot_data = false;
    write_outputs( SpecUtils::append_path( settings.out, "sweep_" + settings.sweep_name + "_" + value ),
                   problems, result.first, sub_settings, meta, result.second );
    summary << result.second.tsv_row( label ) << "\n";
    summary.flush();
    cout << label << ": " << result.second.summary_line() << endl;
  }
  cout << "Sweep summary written to " << SpecUtils::append_path( settings.out, "sweep_summary.tsv" ) << endl;
  return 0;
}//run_search_eval(...)

}//namespace FitPeaksCorpus
