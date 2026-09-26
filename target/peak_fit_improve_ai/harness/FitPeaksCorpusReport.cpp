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
#include <string>
#include <vector>
#include <cstdio>
#include <memory>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>

#include <Wt/WColor.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"

#include "FitPeaksCorpusReport.h"

using namespace std;

namespace FitPeaksCorpus
{

namespace
{
  string fmt( const double v, const int precision = 4 )
  {
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%.*g", precision, v );
    return buffer;
  }

  string tsv_safe( string s )
  {
    for( char &c : s )
    {
      if( (c == '\t') || (c == '\n') || (c == '\r') )
        c = ' ';
    }
    return s;
  }

  string html_escape( const string &s )
  {
    string out;
    for( const char c : s )
    {
      switch( c )
      {
        case '&': out += "&amp;"; break;
        case '<': out += "&lt;"; break;
        case '>': out += "&gt;"; break;
        case '"': out += "&quot;"; break;
        default: out += c;
      }
    }
    return out;
  }

  string json_escape( const string &s )
  {
    string out;
    for( const char c : s )
    {
      switch( c )
      {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default: out += c;
      }
    }
    return out;
  }

  const char *gamma_type_name( const PeakDef::SourceGammaType type )
  {
    switch( type )
    {
      case PeakDef::NormalGamma:       return "gamma";
      case PeakDef::AnnihilationGamma: return "annihilation";
      case PeakDef::SingleEscapeGamma: return "single_escape";
      case PeakDef::DoubleEscapeGamma: return "double_escape";
      case PeakDef::XrayGamma:         return "xray";
    }
    return "unknown";
  }

  string join_warnings( const vector<string> &warnings )
  {
    string out;
    for( const string &w : warnings )
      out += (out.empty() ? "" : " | ") + tsv_safe( w );
    return out;
  }

  ofstream open_or_throw( const string &path )
  {
    ofstream out( path );
    if( !out.good() )
      throw runtime_error( "Could not open '" + path + "' for writing" );
    return out;
  }

  vector<double> quantiles( vector<double> values, const vector<double> &probs )
  {
    vector<double> answer;
    if( values.empty() )
      return vector<double>( probs.size(), 0.0 );
    std::sort( begin(values), end(values) );
    for( const double p : probs )
    {
      const double pos = p * static_cast<double>( values.size() - 1 );
      const size_t lo = static_cast<size_t>( std::floor( pos ) );
      const size_t hi = std::min( values.size() - 1, lo + 1 );
      const double frac = pos - static_cast<double>( lo );
      answer.push_back( values[lo] + frac*(values[hi] - values[lo]) );
    }
    return answer;
  }

  string quantile_str( const vector<double> &values )
  {
    const vector<double> q = quantiles( values, { 0.05, 0.25, 0.5, 0.75, 0.95 } );
    ostringstream out;
    out << "n=" << values.size() << " q5/25/50/75/95=";
    for( size_t i = 0; i < q.size(); ++i )
      out << (i ? "/" : "") << fmt( q[i], 3 );
    return out.str();
  }
}//namespace


AggregateCounts aggregate( const vector<ProblemResult> &results )
{
  AggregateCounts a;
  double norm_sum = 0.0;
  for( const ProblemResult &r : results )
  {
    const ProblemScore &s = r.score;
    a.problems += 1;
    a.mechanical_failures += r.mechanical_failure ? 1 : 0;
    a.nondeterministic += r.nondeterministic ? 1 : 0;
    a.dev_check_failures += r.dev_check_failures;
    a.truth_scored += s.n_truth_scored;
    a.matched += s.n_matched;
    a.missed_strong += s.missed_strong;
    a.missed_moderate += s.missed_moderate;
    a.missed_weak += s.missed_weak;
    a.extra_significant += s.extra_significant;
    a.extra_weak += s.extra_weak;
    a.ghosts += s.extra_ghost;
    a.pairs += s.pairs;
    a.share_disagree += s.share_disagree;
    a.rois_compared += s.rois_compared;
    a.family_disagree += s.family_disagree;
    a.extent_sides += s.extent_sides;
    a.extent_gt1fwhm_sides += s.extent_gt1fwhm_sides;
    a.extra_legit += s.extra_legit;
    a.extra_xray += s.extra_xray;
    a.extra_bkg_line += s.extra_bkg_line;
    a.truth_on_bkg_line += s.truth_on_bkg_line;
    a.truth_merged_away += s.truth_merged_away;
    a.truth_clusters += s.truth_clusters;
    a.truth_xray_clusters += s.truth_xray_clusters;
    a.truth_fit_graded += s.truth_fit_graded;
    a.truth_ref_graded += s.truth_ref_graded;
    a.truth_fit_bad += s.truth_fit_bad;
    a.truth_ref_bad += s.truth_ref_bad;
    a.truth_fit_abs_pull += s.truth_fit_abs_pull;
    a.truth_ref_abs_pull += s.truth_ref_abs_pull;
    a.cost_truth_area += s.cost_truth_area;
    a.sum_raw_cost += s.raw_cost();
    norm_sum += s.norm_cost();
    a.cost_missed += s.cost_missed;
    a.cost_extra += s.cost_extra;
    a.cost_share += s.cost_share;
    a.cost_family += s.cost_family;
    a.cost_extent += s.cost_extent;
    a.cost_area += s.cost_area;
    a.cost_mean += s.cost_mean;
    a.cost_failure += s.cost_failure;
    a.wall_seconds += r.wall_seconds;
  }
  a.mean_norm_cost = a.problems ? (norm_sum / static_cast<double>(a.problems)) : 0.0;
  return a;
}//aggregate


string summary_line( const AggregateCounts &a )
{
  ostringstream out;
  out << "problems=" << a.problems
      << " fail=" << a.mechanical_failures
      << " nondet=" << a.nondeterministic
      << " devchk=" << a.dev_check_failures
      << " raw=" << fmt( a.sum_raw_cost, 6 )
      << " norm=" << fmt( a.mean_norm_cost, 4 )
      << " matched=" << a.matched << "/" << a.truth_scored
      << " missed(s/m/w)=" << a.missed_strong << "/" << a.missed_moderate << "/" << a.missed_weak
      << " extra(sig/weak/ghost)=" << a.extra_significant << "/" << a.extra_weak << "/" << a.ghosts
      << " share=" << a.share_disagree << "/" << a.pairs
      << " family=" << a.family_disagree << "/" << a.rois_compared
      << " extent>1=" << a.extent_gt1fwhm_sides << "/" << a.extent_sides
      << " cost(miss/extra/share/fam/ext/area/mean/fail/truthA)="
      << fmt(a.cost_missed,4) << "/" << fmt(a.cost_extra,4) << "/" << fmt(a.cost_share,4) << "/"
      << fmt(a.cost_family,4) << "/" << fmt(a.cost_extent,4) << "/" << fmt(a.cost_area,4) << "/"
      << fmt(a.cost_mean,4) << "/" << fmt(a.cost_failure,4) << "/" << fmt(a.cost_truth_area,4);
  if( a.truth_clusters )
  {
    out << " legit=" << a.extra_legit << " xray=" << a.extra_xray << " bkgline=" << a.extra_bkg_line << " refbkg=" << a.truth_on_bkg_line << " refmerged=" << a.truth_merged_away
        << " truthA|pull|(fit/ref)=" << fmt( a.truth_fit_graded ? a.truth_fit_abs_pull/a.truth_fit_graded : 0.0, 3 )
        << "/" << fmt( a.truth_ref_graded ? a.truth_ref_abs_pull/a.truth_ref_graded : 0.0, 3 )
        << " xrayclusters=" << a.truth_xray_clusters << " bad(fit/ref)=" << a.truth_fit_bad << "/" << a.truth_ref_bad << " of " << a.truth_fit_graded << "/" << a.truth_ref_graded;
  }
  out << " wall=" << fmt( a.wall_seconds, 4 ) << "s";
  return out.str();
}//summary_line


string summary_tsv_header()
{
  return "label\tproblems\tmechanical_failures\tnondeterministic\tdev_check_failures\tsum_raw_cost"
         "\tmean_norm_cost\ttruth_scored\tmatched\tmissed_strong\tmissed_moderate\tmissed_weak"
         "\textra_significant\textra_weak\tghosts\tpairs\tshare_disagree\trois_compared"
         "\tfamily_disagree\textent_sides\textent_gt1fwhm_sides\tcost_missed\tcost_extra"
         "\tcost_share\tcost_family\tcost_extent\tcost_area\tcost_mean\tcost_failure\twall_seconds"
         "\textra_legit\textra_bkg_line\ttruth_clusters\ttruth_fit_graded\ttruth_ref_graded"
         "\ttruth_fit_abs_pull\ttruth_ref_abs_pull\ttruth_fit_bad\ttruth_ref_bad\tcost_truth_area";
}


string summary_tsv_row( const string &label, const AggregateCounts &a )
{
  ostringstream out;
  out << tsv_safe(label) << '\t' << a.problems << '\t' << a.mechanical_failures << '\t'
      << a.nondeterministic << '\t' << a.dev_check_failures << '\t' << fmt(a.sum_raw_cost,8) << '\t'
      << fmt(a.mean_norm_cost,6) << '\t' << a.truth_scored << '\t' << a.matched << '\t'
      << a.missed_strong << '\t' << a.missed_moderate << '\t' << a.missed_weak << '\t'
      << a.extra_significant << '\t' << a.extra_weak << '\t' << a.ghosts << '\t' << a.pairs << '\t'
      << a.share_disagree << '\t' << a.rois_compared << '\t' << a.family_disagree << '\t'
      << a.extent_sides << '\t' << a.extent_gt1fwhm_sides << '\t' << fmt(a.cost_missed,6) << '\t'
      << fmt(a.cost_extra,6) << '\t' << fmt(a.cost_share,6) << '\t' << fmt(a.cost_family,6) << '\t'
      << fmt(a.cost_extent,6) << '\t' << fmt(a.cost_area,6) << '\t' << fmt(a.cost_mean,6) << '\t'
      << fmt(a.cost_failure,6) << '\t' << fmt(a.wall_seconds,6) << '\t'
      << a.extra_legit << '\t' << a.extra_bkg_line << '\t' << a.truth_clusters << '\t'
      << a.truth_fit_graded << '\t' << a.truth_ref_graded << '\t' << fmt(a.truth_fit_abs_pull,6) << '\t'
      << fmt(a.truth_ref_abs_pull,6) << '\t' << a.truth_fit_bad << '\t' << a.truth_ref_bad << '\t'
      << fmt(a.cost_truth_area,6);
  return out.str();
}


void write_run_meta( const string &path, const RunMeta &meta )
{
  ofstream out = open_or_throw( path );
  out << "command_line: " << meta.command_line << "\n"
      << "started: " << meta.started << "\n"
      << "git_head: " << meta.git_head << "\n"
      << "git_dirty: " << meta.git_dirty << "\n"
      << "binary_sha256: " << meta.binary_sha256 << "\n"
      << "corpus_dir: " << meta.corpus_dir << "\n"
      << "mode: " << meta.mode << "\n"
      << "weights: " << meta.weights_text << "\n"
      << "config:\n" << meta.config_text << "\n";
}


void write_per_problem_tsv( const string &path, const vector<ProblemResult> &results )
{
  ofstream out = open_or_throw( path );
  out << "id\tmode\tsources\tstatus\tmechanical_failure\tnondeterministic\tdev_check_failures"
         "\traw_cost\tnorm_cost\tn_truth_scored\tn_truth_dontcare\tn_fitted\tn_matched"
         "\tmissed_strong\tmissed_moderate\tmissed_weak\textra_significant\textra_weak\tghost\tneutral"
         "\tpairs\tshare_disagree\trois_compared\tfamily_disagree\textent_sides\textent_gt1fwhm_sides"
         "\tcost_missed\tcost_extra\tcost_share\tcost_family\tcost_extent\tcost_area\tcost_mean"
         "\tcost_failure\twall_seconds\textra_legit\textra_bkg_line\ttruth_clusters\ttruth_fit_graded"
         "\ttruth_ref_graded\ttruth_fit_abs_pull\ttruth_ref_abs_pull\ttruth_fit_bad\ttruth_ref_bad"
         "\tcost_truth_area\terror\twarnings\n";
  for( const ProblemResult &r : results )
  {
    const ProblemScore &s = r.score;
    out << tsv_safe(r.id) << '\t' << r.mode << '\t' << tsv_safe(r.sources) << '\t' << r.status << '\t'
        << (r.mechanical_failure ? 1 : 0) << '\t' << (r.nondeterministic ? 1 : 0) << '\t'
        << r.dev_check_failures << '\t' << fmt(s.raw_cost(),8) << '\t' << fmt(s.norm_cost(),6) << '\t'
        << s.n_truth_scored << '\t' << s.n_truth_dontcare << '\t' << s.n_fitted << '\t' << s.n_matched << '\t'
        << s.missed_strong << '\t' << s.missed_moderate << '\t' << s.missed_weak << '\t'
        << s.extra_significant << '\t' << s.extra_weak << '\t' << s.extra_ghost << '\t' << s.neutral << '\t'
        << s.pairs << '\t' << s.share_disagree << '\t' << s.rois_compared << '\t' << s.family_disagree << '\t'
        << s.extent_sides << '\t' << s.extent_gt1fwhm_sides << '\t'
        << fmt(s.cost_missed,6) << '\t' << fmt(s.cost_extra,6) << '\t' << fmt(s.cost_share,6) << '\t'
        << fmt(s.cost_family,6) << '\t' << fmt(s.cost_extent,6) << '\t' << fmt(s.cost_area,6) << '\t'
        << fmt(s.cost_mean,6) << '\t' << fmt(s.cost_failure,6) << '\t' << fmt(r.wall_seconds,5) << '\t'
        << s.extra_legit << '\t' << s.extra_bkg_line << '\t' << s.truth_clusters << '\t'
        << s.truth_fit_graded << '\t' << s.truth_ref_graded << '\t' << fmt(s.truth_fit_abs_pull,6) << '\t'
        << fmt(s.truth_ref_abs_pull,6) << '\t' << s.truth_fit_bad << '\t' << s.truth_ref_bad << '\t'
        << fmt(s.cost_truth_area,6) << '\t'
        << tsv_safe(r.error) << '\t' << join_warnings(r.warnings) << '\n';
  }
}//write_per_problem_tsv


void write_planned_rois_tsv( const string &path, const vector<ProblemResult> &results )
{
  ofstream out( path.c_str() );
  if( !out )
    throw runtime_error( "Could not open '" + path + "'" );
  out << "id\tlower\tupper\n";
  for( const ProblemResult &r : results )
    for( const std::pair<double,double> &roi : r.planned_rois )
      out << tsv_safe(r.id) << '\t' << fmt(roi.first,6) << '\t' << fmt(roi.second,6) << '\n';
}//write_planned_rois_tsv


void write_per_peak_tsv( const string &path, const vector<ProblemResult> &results )
{
  ofstream out = open_or_throw( path );
  out << "id\tset\tindex\tenergy\tfwhm\tamplitude\tamplitude_uncert\tcontinuum_counts\tz_det"
         "\tcontinuum_type\tfamily\troi_index\troi_lower\troi_upper\troi_width_fwhm\troi_num_peaks"
         "\tsource\tgamma_type\tdont_care\tverdict\tmatch_energy\tmatch_amplitude\tauto_dist_fwhm\tauto_z"
         "\ttruth_energy\ttruth_area\ttruth_z\n";
  const auto write_set = [&out]( const string &id, const char *label, const PeakSet &set, const PeakSet &other,
                                 const PeakSet *autosearch, const InjectTruth *inject ){
    for( size_t i = 0; i < set.peaks.size(); ++i )
    {
      const ScoredPeak &p = set.peaks[i];
      const ScoredRoi &roi = set.rois[p.roi_index];
      // nearest automated-search peak, in FWHM of this peak (diagnostic for the admission gate)
      double auto_dist = -1.0, auto_z = 0.0;
      if( autosearch && (p.fwhm > 0.0) )
      {
        for( const ScoredPeak &a : autosearch->peaks )
        {
          const double d = std::fabs( a.energy - p.energy ) / p.fwhm;
          if( (auto_dist < 0.0) || (d < auto_dist) )
          {
            auto_dist = d;
            auto_z = a.z_det;
          }
        }
      }
      out << tsv_safe(id) << '\t' << label << '\t' << i << '\t' << fmt(p.energy,7) << '\t' << fmt(p.fwhm,5) << '\t'
          << fmt(p.amplitude,7) << '\t' << fmt(p.amplitude_uncert,5) << '\t' << fmt(p.continuum_counts,6) << '\t'
          << fmt(p.z_det,5) << '\t' << PeakContinuum::offset_type_str(p.continuum_type) << '\t'
          << family_name(p.family) << '\t' << p.roi_index << '\t' << fmt(roi.lower,7) << '\t'
          << fmt(roi.upper,7) << '\t' << fmt(roi.width_fwhm,4) << '\t' << roi.peaks.size() << '\t'
          << p.source_name << '\t' << gamma_type_name(p.gamma_type) << '\t' << (p.dont_care ? 1 : 0) << '\t'
          << p.verdict << '\t';
      if( (p.match >= 0) && (static_cast<size_t>(p.match) < other.peaks.size()) )
        out << fmt( other.peaks[p.match].energy, 7 ) << '\t' << fmt( other.peaks[p.match].amplitude, 7 );
      else
        out << '\t';
      out << '\t' << ((auto_dist >= 0.0) ? fmt( auto_dist, 4 ) : string()) << '\t' << ((auto_dist >= 0.0) ? fmt( auto_z, 4 ) : string());
      if( inject && (p.truth_cluster >= 0) && (static_cast<size_t>(p.truth_cluster) < inject->merged_signal.size()) )
      {
        const TruthPhotopeak &t = inject->merged_signal[p.truth_cluster];
        out << '\t' << fmt(t.energy,7) << '\t' << fmt(t.area,6) << '\t' << fmt(t.z_det(),4);
      }else
      {
        out << "\t\t\t";
      }
      out << '\n';
    }
  };
  for( const ProblemResult &r : results )
  {
    write_set( r.id, "truth", r.truth, r.fitted, &r.autosearch, r.inject_truth.get() );
    write_set( r.id, "fitted", r.fitted, r.truth, &r.autosearch, r.inject_truth.get() );
    write_set( r.id, "auto", r.autosearch, r.truth, nullptr, nullptr );
  }
}//write_per_peak_tsv


void write_per_truth_area_tsv( const string &path, const vector<ProblemResult> &results )
{
  ofstream out = open_or_throw( path );
  out << "id\tcluster\tenergy\ttruth_area\ttruth_z\tsigma\tfit_n\tfit_area\tfit_pull\tref_n\tref_area\tref_pull\n";
  for( const ProblemResult &r : results )
  {
    for( const TruthAreaRecord &rec : r.score.truth_area_records )
    {
      out << tsv_safe(r.id) << '\t' << rec.cluster << '\t' << fmt(rec.energy,7) << '\t' << fmt(rec.area,6) << '\t'
          << fmt(rec.z,4) << '\t' << fmt(rec.sigma,5) << '\t' << rec.fit_n << '\t' << fmt(rec.fit_area,6) << '\t'
          << (rec.fit_graded ? fmt(rec.fit_pull,4) : string()) << '\t' << rec.ref_n << '\t' << fmt(rec.ref_area,6) << '\t'
          << (rec.ref_graded ? fmt(rec.ref_pull,4) : string()) << '\n';
    }
  }
}//write_per_truth_area_tsv


void write_per_roi_tsv( const string &path, const vector<ProblemResult> &results )
{
  ofstream out = open_or_throw( path );
  out << "id\tset\troi_index\tlower\tupper\tfirst_channel\tlast_channel\tcontinuum_type\tfamily"
         "\tn_peaks\tdominant_energy\twidth_fwhm\tfit_roi\tfit_family\td_lower_fwhm\td_upper_fwhm"
         "\tcost_family\tcost_extent\n";
  for( const ProblemResult &r : results )
  {
    map<size_t,const RoiRecord *> records;
    for( const RoiRecord &rec : r.score.roi_records )
      records[rec.truth_roi] = &rec;

    const auto write_set = [&]( const char *label, const PeakSet &set, const bool is_truth ){
      for( size_t i = 0; i < set.rois.size(); ++i )
      {
        const ScoredRoi &roi = set.rois[i];
        out << tsv_safe(r.id) << '\t' << label << '\t' << i << '\t' << fmt(roi.lower,7) << '\t' << fmt(roi.upper,7)
            << '\t' << roi.first_channel << '\t' << roi.last_channel << '\t'
            << PeakContinuum::offset_type_str(roi.continuum_type) << '\t' << family_name(roi.family) << '\t'
            << roi.peaks.size() << '\t' << fmt(set.peaks[roi.dominant].energy,7) << '\t' << fmt(roi.width_fwhm,4);
        const auto pos = is_truth ? records.find( i ) : records.end();
        if( pos != records.end() )
        {
          const RoiRecord &rec = *pos->second;
          out << '\t' << rec.fit_roi << '\t' << family_name(rec.fit_family) << '\t' << fmt(rec.d_lower_fwhm,4)
              << '\t' << fmt(rec.d_upper_fwhm,4) << '\t' << fmt(rec.cost_family,4) << '\t' << fmt(rec.cost_extent,4);
        }else
        {
          out << "\t\t\t\t\t\t";
        }
        out << '\n';
      }
    };
    write_set( "truth", r.truth, true );
    write_set( "fitted", r.fitted, false );
  }
}//write_per_roi_tsv


void write_per_pair_tsv( const string &path, const vector<ProblemResult> &results )
{
  ofstream out = open_or_throw( path );
  out << "id\ttruth_energy_a\ttruth_energy_b\tseparation_fwhm\ttruth_share\tfit_share\tagree\n";
  for( const ProblemResult &r : results )
  {
    for( const PairRecord &rec : r.score.pair_records )
    {
      out << tsv_safe(r.id) << '\t' << fmt(r.truth.peaks[rec.truth_a].energy,7) << '\t'
          << fmt(r.truth.peaks[rec.truth_b].energy,7) << '\t' << fmt(rec.separation_fwhm,4) << '\t'
          << (rec.truth_share ? 1 : 0) << '\t' << (rec.fit_share ? 1 : 0) << '\t'
          << ((rec.truth_share == rec.fit_share) ? 1 : 0) << '\n';
    }
  }
}//write_per_pair_tsv


void write_summary_json( const string &path, const vector<ProblemResult> &results, const RunMeta &meta )
{
  const AggregateCounts a = aggregate( results );
  ofstream out = open_or_throw( path );
  out << "{\n  \"meta\": {\n"
      << "    \"command_line\": \"" << json_escape(meta.command_line) << "\",\n"
      << "    \"started\": \"" << json_escape(meta.started) << "\",\n"
      << "    \"git_head\": \"" << json_escape(meta.git_head) << "\",\n"
      << "    \"git_dirty\": \"" << json_escape(meta.git_dirty) << "\",\n"
      << "    \"binary_sha256\": \"" << json_escape(meta.binary_sha256) << "\",\n"
      << "    \"corpus_dir\": \"" << json_escape(meta.corpus_dir) << "\",\n"
      << "    \"mode\": \"" << json_escape(meta.mode) << "\",\n"
      << "    \"weights\": \"" << json_escape(meta.weights_text) << "\"\n  },\n"
      << "  \"aggregate\": {\n"
      << "    \"problems\": " << a.problems << ",\n"
      << "    \"mechanical_failures\": " << a.mechanical_failures << ",\n"
      << "    \"nondeterministic\": " << a.nondeterministic << ",\n"
      << "    \"dev_check_failures\": " << a.dev_check_failures << ",\n"
      << "    \"sum_raw_cost\": " << fmt(a.sum_raw_cost,8) << ",\n"
      << "    \"mean_norm_cost\": " << fmt(a.mean_norm_cost,6) << ",\n"
      << "    \"truth_scored\": " << a.truth_scored << ",\n"
      << "    \"matched\": " << a.matched << ",\n"
      << "    \"missed_strong\": " << a.missed_strong << ",\n"
      << "    \"missed_moderate\": " << a.missed_moderate << ",\n"
      << "    \"missed_weak\": " << a.missed_weak << ",\n"
      << "    \"extra_significant\": " << a.extra_significant << ",\n"
      << "    \"extra_weak\": " << a.extra_weak << ",\n"
      << "    \"ghosts\": " << a.ghosts << ",\n"
      << "    \"pairs\": " << a.pairs << ",\n"
      << "    \"share_disagree\": " << a.share_disagree << ",\n"
      << "    \"rois_compared\": " << a.rois_compared << ",\n"
      << "    \"family_disagree\": " << a.family_disagree << ",\n"
      << "    \"extent_sides\": " << a.extent_sides << ",\n"
      << "    \"extent_gt1fwhm_sides\": " << a.extent_gt1fwhm_sides << ",\n"
      << "    \"cost_missed\": " << fmt(a.cost_missed,6) << ",\n"
      << "    \"cost_extra\": " << fmt(a.cost_extra,6) << ",\n"
      << "    \"cost_share\": " << fmt(a.cost_share,6) << ",\n"
      << "    \"cost_family\": " << fmt(a.cost_family,6) << ",\n"
      << "    \"cost_extent\": " << fmt(a.cost_extent,6) << ",\n"
      << "    \"cost_area\": " << fmt(a.cost_area,6) << ",\n"
      << "    \"cost_mean\": " << fmt(a.cost_mean,6) << ",\n"
      << "    \"cost_failure\": " << fmt(a.cost_failure,6) << ",\n"
      << "    \"wall_seconds\": " << fmt(a.wall_seconds,6) << "\n  },\n"
      << "  \"problems\": [\n";
  for( size_t i = 0; i < results.size(); ++i )
  {
    const ProblemResult &r = results[i];
    out << "    {\"id\": \"" << json_escape(r.id) << "\", \"status\": \"" << json_escape(r.status)
        << "\", \"raw_cost\": " << fmt(r.score.raw_cost(),8) << ", \"norm_cost\": " << fmt(r.score.norm_cost(),6)
        << ", \"mechanical_failure\": " << (r.mechanical_failure ? "true" : "false") << "}"
        << ((i + 1 < results.size()) ? ",\n" : "\n");
  }
  out << "  ]\n}\n";
}//write_summary_json


namespace
{
  string peaks_json_colored( const vector<ScoredPeak> &peaks, const map<string,string> &verdict_colors,
                             const shared_ptr<const SpecUtils::Measurement> &spectrum )
  {
    vector<shared_ptr<const PeakDef>> colored;
    for( const ScoredPeak &p : peaks )
    {
      if( !p.peak )
        continue;
      shared_ptr<PeakDef> copy = make_shared<PeakDef>( *p.peak );
      const auto pos = verdict_colors.find( p.verdict );
      copy->setLineColor( Wt::WColor( (pos != end(verdict_colors)) ? pos->second : "#1f77b4" ) );
      colored.push_back( copy );
    }
    return PeakDef::peak_json( colored, spectrum, Wt::WColor(), -1 );
  }

  string verdict_css_class( const string &verdict )
  {
    if( verdict == "matched" ) return "matched";
    if( verdict == "missed" ) return "missed";
    if( (verdict == "extra") || (verdict == "extra_weak") || (verdict == "extra_bkg") ) return "extra";
    if( verdict == "ghost" ) return "ghost";
    if( verdict == "neutral" ) return "neutral";
    if( (verdict == "legit") || (verdict == "xray") ) return "legit";
    if( (verdict == "dontcare_merged") || (verdict == "dontcare_bkg") ) return "excused";
    return "dontcare";
  }

  string pull_css_class( const bool graded, const double pull )
  {
    if( !graded )
      return "dontcare";
    const double a = std::fabs( pull );
    return (a <= 2.0) ? "pullgood" : ((a <= 4.0) ? "pullwarn" : "pullbad");
  }
}//namespace


void write_gallery_html( const string &path, const string &resources_dir, const string &title,
                         const RunMeta &meta, const vector<ProblemResult> &results, const size_t max_charts,
                         const bool include_background )
{
  ofstream out = open_or_throw( path );
  D3SpectrumExport::write_html_page_header( out, title, resources_dir );

  out << "<body>\n<style>\n"
         "body{font-family:-apple-system,Helvetica,Arial,sans-serif;margin:12px;font-size:13px}\n"
         "table.summary{border-collapse:collapse;font-size:12px} table.summary th{cursor:pointer;background:#eee;position:sticky;top:0}\n"
         "table.summary td,table.summary th{border:1px solid #ccc;padding:2px 5px;text-align:right}\n"
         "table.summary td:first-child,table.summary td:nth-child(2){text-align:left}\n"
         "tr.fail{background:#fdd} tr.clickable{cursor:pointer}\n"
         "fieldset{margin:14px 0;border:1px solid #bbb} legend{font-weight:bold}\n"
         ".chart{width:96%;height:380px;margin:0 auto}\n"
         "table.peaks{font-size:11px;border-collapse:collapse;margin:6px 0} table.peaks td,table.peaks th{border:1px solid #ddd;padding:1px 5px;text-align:right}\n"
         "table.peaks td:first-child,table.peaks th{text-align:left}\n"
         ".matched{color:#1f77b4} .missed{color:#e6550d;font-weight:bold} .extra{color:#d62728;font-weight:bold} .ghost{color:#999} .neutral{color:#9467bd} .dontcare{color:#bbb} .legit{color:#17becf;font-weight:bold} .excused{color:#8c564b}\n"
         ".pullgood{color:#2ca02c} .pullwarn{color:#e6a000;font-weight:bold} .pullbad{color:#d62728;font-weight:bold}\n"
         ".meta{font-size:11px;color:#555;white-space:pre-wrap}\n"
         "</style>\n";

  out << "<h1>" << html_escape( title ) << "</h1>\n";
  const AggregateCounts agg = aggregate( results );
  out << "<p><b>" << html_escape( summary_line( agg ) ) << "</b></p>\n";
  out << "<div class=\"meta\">" << html_escape( "started " + meta.started + "  git " + meta.git_head + " " + meta.git_dirty
         + "\nbinary " + meta.binary_sha256 + "\n" + meta.command_line + "\nweights " + meta.weights_text ) << "</div>\n";

  out << "<p>Colors: <span class=\"matched\">matched</span>, <span class=\"missed\">missed reference</span>, "
         "<span class=\"extra\">extra fitted</span>, <span class=\"ghost\">ghost (z&lt;1)</span>, "
         "<span class=\"neutral\">neutral (matches a don't-care reference peak)</span>, "
         "<span class=\"legit\">legit (not in the reference, but on a truth photopeak)</span>, "
         "<span class=\"dontcare\">don't-care reference</span>.  Click a column header to sort; click a row to jump to its chart.  "
         "Truth-area pulls: <span class=\"pullgood\">|pull| &le; 2</span>, <span class=\"pullwarn\">2-4</span>, <span class=\"pullbad\">&gt; 4</span>.</p>\n";

  // Summary table
  out << "<button onclick=\"loadAll()\">Load all charts</button>\n";
  out << "<table class=\"summary\" id=\"summary\"><thead><tr>";
  const vector<string> columns = { "id", "sources", "status", "raw", "norm", "miss s/m/w", "extra sig/weak/ghost",
                                   "legit/bkg", "share dis/pairs", "family dis/rois", "extent>1/sides",
                                   "truth |pull| fit", "truth |pull| ref", "truth bad fit/ref", "devchk", "wall s" };
  for( size_t c = 0; c < columns.size(); ++c )
    out << "<th onclick=\"sortTable(" << c << ")\">" << html_escape( columns[c] ) << "</th>";
  out << "</tr></thead><tbody>\n";
  for( size_t i = 0; i < results.size(); ++i )
  {
    const ProblemResult &r = results[i];
    const ProblemScore &s = r.score;
    out << "<tr class=\"clickable" << (r.mechanical_failure ? " fail" : "") << "\" onclick=\"jumpTo('prob_" << i << "')\">"
        << "<td>" << html_escape(r.id) << "</td><td>" << html_escape(r.sources) << "</td><td>" << html_escape(r.status) << "</td>"
        << "<td data-v=\"" << fmt(s.raw_cost(),6) << "\">" << fmt(s.raw_cost(),4) << "</td>"
        << "<td data-v=\"" << fmt(s.norm_cost(),6) << "\">" << fmt(s.norm_cost(),3) << "</td>"
        << "<td data-v=\"" << (4*s.missed_strong + 2*s.missed_moderate + s.missed_weak) << "\">" << s.missed_strong << "/" << s.missed_moderate << "/" << s.missed_weak << "</td>"
        << "<td data-v=\"" << (4*s.extra_significant + 2*s.extra_weak + s.extra_ghost) << "\">" << s.extra_significant << "/" << s.extra_weak << "/" << s.extra_ghost << "</td>"
        << "<td data-v=\"" << (s.extra_legit + s.extra_bkg_line) << "\">" << s.extra_legit << "/" << s.extra_bkg_line << "</td>"
        << "<td data-v=\"" << s.share_disagree << "\">" << s.share_disagree << "/" << s.pairs << "</td>"
        << "<td data-v=\"" << s.family_disagree << "\">" << s.family_disagree << "/" << s.rois_compared << "</td>"
        << "<td data-v=\"" << s.extent_gt1fwhm_sides << "\">" << s.extent_gt1fwhm_sides << "/" << s.extent_sides << "</td>"
        << "<td data-v=\"" << fmt( s.truth_fit_graded ? s.truth_fit_abs_pull/s.truth_fit_graded : 0.0, 5 ) << "\">"
        << (s.truth_fit_graded ? fmt( s.truth_fit_abs_pull/s.truth_fit_graded, 3 ) : string("-")) << "</td>"
        << "<td data-v=\"" << fmt( s.truth_ref_graded ? s.truth_ref_abs_pull/s.truth_ref_graded : 0.0, 5 ) << "\">"
        << (s.truth_ref_graded ? fmt( s.truth_ref_abs_pull/s.truth_ref_graded, 3 ) : string("-")) << "</td>"
        << "<td data-v=\"" << s.truth_fit_bad << "\">" << s.truth_fit_bad << "/" << s.truth_ref_bad << "</td>"
        << "<td>" << r.dev_check_failures << "</td>"
        << "<td data-v=\"" << fmt(r.wall_seconds,5) << "\">" << fmt(r.wall_seconds,3) << "</td></tr>\n";
  }
  out << "</tbody></table>\n";

  out << "<script>\n"
         "function sortTable(col){const t=document.getElementById('summary');const rows=Array.from(t.tBodies[0].rows);"
         "const asc=(t.dataset.sortcol==col)&&(t.dataset.asc!='1');"
         "rows.sort((a,b)=>{const x=a.cells[col].dataset.v!==undefined?a.cells[col].dataset.v:a.cells[col].textContent;"
         "const y=b.cells[col].dataset.v!==undefined?b.cells[col].dataset.v:b.cells[col].textContent;"
         "const nx=parseFloat(x),ny=parseFloat(y);const c=(isNaN(nx)||isNaN(ny))?x.localeCompare(y):(nx-ny);return asc?c:-c;});"
         "rows.forEach(r=>t.tBodies[0].appendChild(r));t.dataset.sortcol=col;t.dataset.asc=asc?'1':'0';}\n"
         "const chartInits={};\n"
         "function ensureInit(id){const f=chartInits[id];if(f){delete chartInits[id];try{f();}catch(e){console.error(id,e);}}}\n"
         "function jumpTo(id){ensureInit(id);const el=document.getElementById(id);if(el)el.scrollIntoView();}\n"
         "function loadAll(){Object.keys(chartInits).forEach(ensureInit);}\n"
         "function showSet(div,which){const c=window['spec_chart_'+div];const d=window['data_'+div];if(!c||!d)return;"
         "d.spectra[0].peaks=(which==='truth')?window['truth_peaks_'+div]:window['fit_peaks_'+div];c.setData(d);}\n"
         "function toggleTruth(div,on){const c=window['spec_chart_'+div];if(!c)return;"
         "c.setReferenceLines(on?window['truth_reflines_'+div]:window['base_reflines_'+div]);}\n"
         "</script>\n";

  // Per-problem sections, worst raw cost first
  vector<size_t> order( results.size() );
  for( size_t i = 0; i < order.size(); ++i )
    order[i] = i;
  std::stable_sort( begin(order), end(order), [&results]( const size_t a, const size_t b ){
    return results[a].score.raw_cost() > results[b].score.raw_cost();
  } );

  const map<string,string> fitted_colors = {
    { "matched", "#1f77b4" }, { "extra", "#d62728" }, { "extra_weak", "#d62728" }, { "extra_bkg", "#d62728" },
    { "ghost", "#999999" }, { "neutral", "#9467bd" }, { "legit", "#17becf" }, { "xray", "#17becf" },
    { "extra_bkg", "#d62728" }
  };
  const map<string,string> truth_colors = {
    { "matched", "#2ca02c" }, { "missed", "#e6550d" }, { "dontcare", "#bbbbbb" }, { "matched_dontcare", "#bbbbbb" },
    { "dontcare_merged", "#8c564b" }, { "dontcare_bkg", "#8c564b" }
  };

  size_t charts_written = 0;
  for( const size_t i : order )
  {
    const ProblemResult &r = results[i];
    const ProblemScore &s = r.score;
    const string div = "chart_" + std::to_string( i );
    const bool write_chart = r.foreground && (charts_written < max_charts);

    out << "<fieldset id=\"prob_" << i << "\"" << (write_chart ? " data-chart=\"1\"" : "") << ">\n<legend>"
        << html_escape( r.id ) << " &mdash; " << html_escape( r.sources ) << " &mdash; " << html_escape( r.status )
        << " &mdash; raw " << fmt( s.raw_cost(), 4 ) << " (norm " << fmt( s.norm_cost(), 3 ) << ")</legend>\n";

    if( !r.spectrum_path.empty() )
      out << "<div class=\"meta\"><code>" << html_escape( r.spectrum_path ) << "</code></div>\n";

    if( !r.error.empty() )
      out << "<p class=\"extra\">" << html_escape( r.error ) << "</p>\n";
    for( const string &w : r.warnings )
      out << "<div class=\"meta\">" << html_escape( w ) << "</div>\n";

    if( write_chart )
    {
      charts_written += 1;

      double x_min = r.foreground->gamma_energy_min(), x_max = r.foreground->gamma_energy_max();
      double lo = 1.0e30, hi = -1.0e30;
      for( const PeakSet *set : { &r.truth, &r.fitted } )
      {
        for( const ScoredRoi &roi : set->rois )
        {
          lo = std::min( lo, roi.lower );
          hi = std::max( hi, roi.upper );
        }
      }
      if( lo < hi )
      {
        x_min = std::max( x_min, lo - 15.0 );
        x_max = std::min( x_max, hi + 15.0 );
      }

      // Missed reference peaks as reference lines (visible in both views); the inject truth
      // photopeaks (signal and background, heights by sqrt(area)) as a second, toggleable set.
      map<string,string> reference_lines;
      string truth_reflines_json = "[]";
      {
        string lines;
        for( const ScoredPeak &p : r.truth.peaks )
        {
          if( p.verdict != "missed" )
            continue;
          lines += (lines.empty() ? "" : ",") + string("{\"e\":") + fmt( p.energy, 7 ) + ",\"h\":1,\"particle\":\"gamma\"}";
        }
        if( !lines.empty() )
          reference_lines["Missed reference"] = "{\"color\":\"#e6550d\",\"parent\":\"Missed reference\",\"lines\":[" + lines + "]}";

        if( r.inject_truth )
        {
          // Heights are the expected areas scaled linearly onto 0-1, the strongest line being 1.
          const auto set_json = []( const vector<TruthPhotopeak> &peaks, const string &parent, const string &color ) -> string {
            double max_area = 0.0;
            for( const TruthPhotopeak &t : peaks )
              max_area = std::max( max_area, t.area );
            if( !(max_area > 0.0) )
              return string();
            string js;
            for( const TruthPhotopeak &t : peaks )
            {
              js += (js.empty() ? "" : ",") + string("{\"e\":") + fmt( t.energy, 7 ) + ",\"h\":" + fmt( std::max( 0.0, t.area ) / max_area, 5 )
                    + ",\"particle\":\"gamma\",\"decay\":\"A=" + fmt( t.area, 4 ) + " z=" + fmt( t.z_det(), 3 ) + "\"}";
            }
            return js.empty() ? string() : ("{\"color\":\"" + color + "\",\"parent\":\"" + parent + "\",\"lines\":[" + js + "]}");
          };
          // Only the source's own photopeaks: the background lines are not what the fit is asked to
          // find, and drawing them makes the source's lines hard to pick out.
          const string sig = set_json( r.inject_truth->merged_signal, "Truth " + r.inject_truth->source_name, "#2ca02c" );
          string all;
          for( const string &part : { reference_lines.count("Missed reference") ? reference_lines["Missed reference"] : string(), sig } )
          {
            if( !part.empty() )
              all += (all.empty() ? "" : ",") + part;
          }
          truth_reflines_json = "[" + all + "]";
        }
      }

      const D3SpectrumExport::D3SpectrumChartOptions chart_options(
        r.id, "Energy (keV)", "Counts/Channel", "", true, false, false, true, true,
        false, false, false, false, false, false, false, false, false,
        static_cast<float>(x_min), static_cast<float>(x_max), reference_lines );

      D3SpectrumExport::D3SpectrumOptions fg_opts;
      fg_opts.line_color = "black";
      fg_opts.title = r.id;
      fg_opts.display_scale_factor = 1.0;
      fg_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;
      fg_opts.peaks_json = peaks_json_colored( r.fitted.peaks, fitted_colors, r.foreground );

      const string truth_peaks_json = peaks_json_colored( r.truth.peaks, truth_colors, r.foreground );

      vector<pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions>> measurements;
      measurements.emplace_back( r.foreground.get(), fg_opts );
      if( include_background && r.background && (r.background->live_time() > 0.0f) )
      {
        D3SpectrumExport::D3SpectrumOptions bg_opts;
        bg_opts.line_color = "steelblue";
        bg_opts.title = "Background";
        bg_opts.display_scale_factor = r.foreground->live_time() / r.background->live_time();
        bg_opts.spectrum_type = SpecUtils::SpectrumType::Background;
        measurements.emplace_back( r.background.get(), bg_opts );
      }

      out << "<div id=\"" << div << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>\n"
          << "<div><button onclick=\"showSet('" << div << "','fit')\">Fitted peaks</button> "
          << "<button onclick=\"showSet('" << div << "','truth')\">"
          << (r.truth_is_inject ? "Truth photopeaks" : "Reference peaks") << "</button>";
      if( r.inject_truth && !r.inject_truth->merged_signal.empty() )
        out << " <label><input type=\"checkbox\" onchange=\"toggleTruth('" << div << "',this.checked)\"> Truth photopeak lines ("
            << html_escape( r.inject_truth->source_name ) << " only; height = area / strongest area)</label>";
      else if( r.inject_truth )
        out << " <span class=\"meta\">(the truth file for this problem lists no photopeaks)</span>";
      out << "</div>\n";

      out << "<script>\nchartInits['prob_" << i << "'] = function(){\n";
      D3SpectrumExport::write_js_for_chart( out, div, chart_options.m_dataTitle, chart_options.m_xAxisTitle, chart_options.m_yAxisTitle );
      D3SpectrumExport::write_and_set_data_for_chart( out, div, measurements );

      // The reference view reuses the same spectrum data and only swaps the peak list.
      out << "var truth_peaks_" << div << " = " << (truth_peaks_json.empty() ? "[]" : truth_peaks_json) << ";\n";

      D3SpectrumExport::write_set_options_for_chart( out, div, chart_options );
      out << "spec_chart_" << div << ".setReferenceLines( reference_lines_" << div << " );\n"
          << "window['base_reflines_" << div << "'] = reference_lines_" << div << ";\n"
          << "window['truth_reflines_" << div << "'] = " << truth_reflines_json << ";\n"
          << "window['spec_chart_" << div << "'] = spec_chart_" << div << ";\n"
          << "window['data_" << div << "'] = data_" << div << ";\n"
          << "window['fit_peaks_" << div << "'] = data_" << div << ".spectra[0].peaks;\n"
          << "window['truth_peaks_" << div << "'] = truth_peaks_" << div << ";\n"
          << "spec_chart_" << div << ".handleResize();\n"
          << "};\n</script>\n";
    }//if( write_chart )

    // Reference peaks table
    out << "<details><summary>" << (r.truth_is_inject ? "Truth photopeaks (" : "Reference peaks (")
        << r.truth.peaks.size() << ") and fitted extras</summary>\n"
        << "<table class=\"peaks\"><tr><th>" << (r.truth_is_inject ? "truth E" : "ref E") << "</th><th>source</th><th>type</th><th>A</th><th>z</th><th>family</th>"
           "<th>ROI</th><th>verdict</th><th>fit E</th><th>fit A</th><th>fit z</th><th>fit family</th><th>fit ROI</th></tr>\n";
    for( const ScoredPeak &p : r.truth.peaks )
    {
      const ScoredRoi &roi = r.truth.rois[p.roi_index];
      out << "<tr class=\"" << verdict_css_class( p.verdict ) << "\"><td>" << fmt(p.energy,6) << "</td><td>" << html_escape(p.source_name)
          << "</td><td>" << gamma_type_name(p.gamma_type) << "</td><td>" << fmt(p.amplitude,5) << "</td><td>" << fmt(p.z_det,3)
          << "</td><td>" << family_name(p.family) << "</td><td>[" << fmt(roi.lower,5) << ", " << fmt(roi.upper,5) << "]</td><td>"
          << p.verdict << "</td>";
      if( p.match >= 0 )
      {
        const ScoredPeak &f = r.fitted.peaks[p.match];
        const ScoredRoi &froi = r.fitted.rois[f.roi_index];
        out << "<td>" << fmt(f.energy,6) << "</td><td>" << fmt(f.amplitude,5) << "</td><td>" << fmt(f.z_det,3) << "</td><td>"
            << family_name(f.family) << "</td><td>[" << fmt(froi.lower,5) << ", " << fmt(froi.upper,5) << "]</td>";
      }else
      {
        out << "<td></td><td></td><td></td><td></td><td></td>";
      }
      out << "</tr>\n";
    }
    for( const ScoredPeak &f : r.fitted.peaks )
    {
      if( f.match >= 0 )
        continue;
      const ScoredRoi &froi = r.fitted.rois[f.roi_index];
      out << "<tr class=\"" << verdict_css_class( f.verdict ) << "\"><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td>"
          << f.verdict << "</td><td>" << fmt(f.energy,6) << "</td><td>" << fmt(f.amplitude,5) << "</td><td>" << fmt(f.z_det,3)
          << "</td><td>" << family_name(f.family) << "</td><td>[" << fmt(froi.lower,5) << ", " << fmt(froi.upper,5) << "] "
          << html_escape(f.source_name) << "</td></tr>\n";
    }
    out << "</table>\n";

    if( !s.truth_area_records.empty() )
    {
      out << "<details><summary>Truth areas (" << s.truth_area_records.size() << " photopeaks; fit |pull| "
          << (s.truth_fit_graded ? fmt( s.truth_fit_abs_pull/s.truth_fit_graded, 3 ) : string("-")) << ", reference |pull| "
          << (s.truth_ref_graded ? fmt( s.truth_ref_abs_pull/s.truth_ref_graded, 3 ) : string("-")) << ")</summary>\n"
          << "<table class=\"peaks\"><tr><th>truth E</th><th>truth A</th><th>truth z</th><th>sigma</th>"
             "<th>fit n</th><th>fit A</th><th>fit pull</th><th>ref n</th><th>ref A</th><th>ref pull</th></tr>\n";
      for( const TruthAreaRecord &rec : s.truth_area_records )
      {
        out << "<tr><td>" << fmt(rec.energy,6) << "</td><td>" << fmt(rec.area,5) << "</td><td>" << fmt(rec.z,3) << "</td><td>"
            << fmt(rec.sigma,4) << "</td><td>" << rec.fit_n << "</td><td>" << fmt(rec.fit_area,5) << "</td><td class=\""
            << pull_css_class( rec.fit_graded, rec.fit_pull ) << "\">" << (rec.fit_graded ? fmt(rec.fit_pull,3) : string("-"))
            << "</td><td>" << rec.ref_n << "</td><td>" << fmt(rec.ref_area,5) << "</td><td class=\""
            << pull_css_class( rec.ref_graded, rec.ref_pull ) << "\">" << (rec.ref_graded ? fmt(rec.ref_pull,3) : string("-"))
            << "</td></tr>\n";
      }
      out << "</table></details>\n";
    }

    if( !s.pair_records.empty() || !s.roi_records.empty() )
    {
      out << "<table class=\"peaks\"><tr><th>ROI / pair</th><th>reference</th><th>fitted</th><th>note</th></tr>\n";
      for( const PairRecord &rec : s.pair_records )
      {
        if( rec.truth_share == rec.fit_share )
          continue;
        out << "<tr class=\"extra\"><td>pair " << fmt(r.truth.peaks[rec.truth_a].energy,6) << " / "
            << fmt(r.truth.peaks[rec.truth_b].energy,6) << " (" << fmt(rec.separation_fwhm,3) << " FWHM)</td><td>"
            << (rec.truth_share ? "share" : "separate") << "</td><td>" << (rec.fit_share ? "share" : "separate")
            << "</td><td>share/separate disagreement</td></tr>\n";
      }
      for( const RoiRecord &rec : s.roi_records )
      {
        if( (rec.cost_family <= 0.0) && (rec.cost_extent <= 0.0) )
          continue;
        const ScoredRoi &troi = r.truth.rois[rec.truth_roi];
        const ScoredRoi &froi = r.fitted.rois[rec.fit_roi];
        out << "<tr><td>ROI [" << fmt(troi.lower,5) << ", " << fmt(troi.upper,5) << "]</td><td>" << family_name(rec.truth_family)
            << "</td><td>" << family_name(rec.fit_family) << " [" << fmt(froi.lower,5) << ", " << fmt(froi.upper,5) << "]</td><td>"
            << "dLower=" << fmt(rec.d_lower_fwhm,3) << " dUpper=" << fmt(rec.d_upper_fwhm,3) << " FWHM; cost family "
            << fmt(rec.cost_family,3) << ", extent " << fmt(rec.cost_extent,3) << "</td></tr>\n";
      }
      out << "</table>\n";
    }
    out << "</details>\n</fieldset>\n";
  }//for( const size_t i : order )

  out << "<script>\n"
         "const chartObserver=new IntersectionObserver((entries)=>{entries.forEach(en=>{if(en.isIntersecting)ensureInit(en.target.id);});},{rootMargin:'500px'});\n"
         "document.querySelectorAll('fieldset[data-chart]').forEach(el=>chartObserver.observe(el));\n"
         "</script>\n</body></html>\n";
}//write_gallery_html


void write_plot_data_json( const string &path, const ProblemResult &result )
{
  const shared_ptr<const SpecUtils::Measurement> &fg = result.foreground;
  if( !fg || !fg->num_gamma_channels() || !fg->channel_energies() )
    return;
  const size_t nchannel = fg->num_gamma_channels();
  const vector<float> &energies = *fg->channel_energies();
  const shared_ptr<const vector<float>> counts = fg->gamma_counts();
  if( !counts || (counts->size() < nchannel) || (energies.size() <= nchannel) )
    return;

  ofstream out = open_or_throw( path );
  out << std::setprecision( 8 );
  out << "{\"id\":\"" << json_escape( result.id ) << "\",\"sources\":\"" << json_escape( result.sources )
      << "\",\"status\":\"" << json_escape( result.status ) << "\",\"raw_cost\":" << fmt( result.score.raw_cost(), 6 )
      << ",\"live_time\":" << fg->live_time() << ",\n\"x\":[";
  for( size_t i = 0; i <= nchannel; ++i )
    out << (i ? "," : "") << energies[i];
  out << "],\n\"y\":[";
  for( size_t i = 0; i < nchannel; ++i )
    out << (i ? "," : "") << (*counts)[i];
  out << "]";

  const shared_ptr<const SpecUtils::Measurement> &bg = result.background;
  if( bg && bg->gamma_counts() && (bg->gamma_counts()->size() == nchannel) && (bg->live_time() > 0.0f) )
  {
    const double scale = fg->live_time() / bg->live_time();
    out << ",\n\"bg_scale\":" << scale << ",\"bg_y\":[";
    for( size_t i = 0; i < nchannel; ++i )
      out << (i ? "," : "") << (*bg->gamma_counts())[i];
    out << "]";
  }

  // ROIs of a set: per-channel continuum and per-peak Gaussian counts over the ROI's channels
  const auto write_set = [&]( const char *label, const PeakSet &set ){
    out << ",\n\"" << label << "\":[";
    for( size_t ri = 0; ri < set.rois.size(); ++ri )
    {
      const ScoredRoi &roi = set.rois[ri];
      size_t first = fg->find_gamma_channel( static_cast<float>(roi.lower) );
      size_t last = fg->find_gamma_channel( static_cast<float>(roi.upper) );
      first = std::min( first, nchannel - 1 );
      last = std::min( std::max( last, first ), nchannel - 1 );
      const size_t n = last - first + 1;
      vector<shared_ptr<const PeakDef>> roi_peaks;
      for( const size_t pi : roi.peaks )
        roi_peaks.push_back( set.peaks[pi].peak );
      vector<double> cont( n, 0.0 );
      const shared_ptr<const PeakContinuum> c = roi_peaks.empty() ? nullptr : roi_peaks.front()->continuum();
      if( c && c->parametersProbablySet() )
      {
        try
        {
          c->offset_integral( &energies[first], &cont[0], n, fg, roi_peaks );
        }catch( std::exception & )
        {
          std::fill( begin(cont), end(cont), 0.0 );
        }
      }
      out << (ri ? ",\n" : "\n") << "{\"lower\":" << roi.lower << ",\"upper\":" << roi.upper << ",\"first\":" << first
          << ",\"type\":\"" << PeakContinuum::offset_type_str( roi.continuum_type ) << "\",\"family\":\"" << family_name( roi.family )
          << "\",\"cont\":[";
      for( size_t i = 0; i < n; ++i )
        out << (i ? "," : "") << fmt( cont[i], 6 );
      out << "],\"peaks\":[";
      for( size_t k = 0; k < roi.peaks.size(); ++k )
      {
        const ScoredPeak &p = set.peaks[roi.peaks[k]];
        vector<double> gauss( n, 0.0 );
        if( p.peak )
          p.peak->gauss_integral( &energies[first], &gauss[0], n );
        out << (k ? "," : "") << "{\"mean\":" << p.energy << ",\"fwhm\":" << p.fwhm << ",\"amp\":" << p.amplitude
            << ",\"z\":" << fmt( p.z_det, 4 ) << ",\"source\":\"" << json_escape( p.source_name ) << "\",\"verdict\":\""
            << json_escape( p.verdict ) << "\",\"y\":[";
        for( size_t i = 0; i < n; ++i )
          out << (i ? "," : "") << fmt( gauss[i], 6 );
        out << "]}";
      }
      out << "]}";
    }
    out << "]";
  };
  write_set( "fit", result.fitted );
  write_set( "ref", result.truth );

  out << ",\n\"auto\":[";
  for( size_t i = 0; i < result.autosearch.peaks.size(); ++i )
  {
    const ScoredPeak &p = result.autosearch.peaks[i];
    out << (i ? "," : "") << "{\"mean\":" << p.energy << ",\"fwhm\":" << p.fwhm << ",\"amp\":" << p.amplitude << ",\"z\":" << fmt( p.z_det, 4 ) << "}";
  }
  out << "]";

  if( result.inject_truth )
  {
    out << ",\n\"truth_source\":\"" << json_escape( result.inject_truth->source_name ) << "\",\"truth\":[";
    bool first = true;
    for( const vector<TruthPhotopeak> *set : { &result.inject_truth->merged_signal, &result.inject_truth->merged_background } )
    {
      for( const TruthPhotopeak &t : *set )
      {
        out << (first ? "" : ",") << "{\"e\":" << t.energy << ",\"area\":" << fmt( t.area, 6 ) << ",\"z\":" << fmt( t.z_det(), 4 )
            << ",\"fwhm\":" << t.fwhm << ",\"signal\":" << (t.signal ? "true" : "false") << "}";
        first = false;
      }
    }
    out << "]";
    out << ",\n\"truth_area\":[";
    for( size_t i = 0; i < result.score.truth_area_records.size(); ++i )
    {
      const TruthAreaRecord &rec = result.score.truth_area_records[i];
      out << (i ? "," : "") << "{\"e\":" << rec.energy << ",\"area\":" << fmt( rec.area, 6 ) << ",\"z\":" << fmt( rec.z, 4 )
          << ",\"fit_area\":" << fmt( rec.fit_area, 6 ) << ",\"fit_pull\":" << (rec.fit_graded ? fmt( rec.fit_pull, 4 ) : string("null"))
          << ",\"ref_area\":" << fmt( rec.ref_area, 6 ) << ",\"ref_pull\":" << (rec.ref_graded ? fmt( rec.ref_pull, 4 ) : string("null")) << "}";
    }
    out << "]";
  }
  out << "}\n";
}//write_plot_data_json


void write_result_n42( const CorpusProblem &problem, const ProblemResult &result, const string &path )
{
  SpecMeas meas;
  if( !meas.load_N42_file( problem.path ) && !meas.load_file( problem.path, SpecUtils::ParserType::Auto ) )
    throw runtime_error( "Could not re-load '" + problem.path + "'" );

  const set<int> samples{ problem.foreground_sample };
  std::deque<shared_ptr<const PeakDef>> fitted( begin(result.fitted_peaks), end(result.fitted_peaks) );
  meas.setPeaks( fitted, samples );

  auto truth = make_shared<std::deque<shared_ptr<const PeakDef>>>( begin(result.truth_peaks), end(result.truth_peaks) );
  meas.setAutomatedSearchPeaks( samples, truth );

  if( !meas.save2012N42File( path ) )
    throw runtime_error( "Could not write '" + path + "'" );
}//write_result_n42


void write_truth_stats( const string &path, const vector<ProblemResult> &results )
{
  ofstream out = open_or_throw( path );

  map<string,size_t> family_rois, family_peaks, type_rois;
  vector<double> kept_z, ext_lo_all, ext_hi_all, ext_lo_single, ext_hi_single, width_all, width_single;
  vector<double> step_z, linear_z, poly_z, poly_width, poly_energy;
  map<int,size_t> peaks_per_roi;
  size_t total_peaks = 0, total_rois = 0, dontcare = 0;

  const vector<double> sep_bins = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 1.0e9 };
  vector<size_t> share_counts( sep_bins.size() - 1, 0 ), separate_counts( sep_bins.size() - 1, 0 );

  for( const ProblemResult &r : results )
  {
    const PeakSet &t = r.truth;
    total_peaks += t.peaks.size();
    total_rois += t.rois.size();
    for( const ScoredPeak &p : t.peaks )
    {
      family_peaks[family_name(p.family)] += 1;
      kept_z.push_back( p.z_det );
      dontcare += p.dont_care ? 1 : 0;
    }

    // adjacent pairs (all reference peaks)
    for( size_t k = 1; k < t.peaks.size(); ++k )
    {
      const ScoredPeak &a = t.peaks[k-1];
      const ScoredPeak &b = t.peaks[k];
      const double mid = 0.5*(a.fwhm + b.fwhm);
      if( mid <= 0.0 )
        continue;
      const double sep = (b.energy - a.energy) / mid;
      for( size_t bin = 0; bin + 1 < sep_bins.size(); ++bin )
      {
        if( (sep >= sep_bins[bin]) && (sep < sep_bins[bin+1]) )
        {
          if( a.roi_index == b.roi_index )
            share_counts[bin] += 1;
          else
            separate_counts[bin] += 1;
          break;
        }
      }
    }

    for( const ScoredRoi &roi : t.rois )
    {
      family_rois[family_name(roi.family)] += 1;
      type_rois[PeakContinuum::offset_type_str(roi.continuum_type)] += 1;
      peaks_per_roi[static_cast<int>(roi.peaks.size())] += 1;
      const ScoredPeak &first = t.peaks[roi.peaks.front()];
      const ScoredPeak &last = t.peaks[roi.peaks.back()];
      const ScoredPeak &dom = t.peaks[roi.dominant];
      if( (first.fwhm > 0.0) && (last.fwhm > 0.0) )
      {
        const double lo = (first.energy - roi.lower) / first.fwhm;
        const double hi = (roi.upper - last.energy) / last.fwhm;
        ext_lo_all.push_back( lo );
        ext_hi_all.push_back( hi );
        width_all.push_back( roi.width_fwhm );
        if( roi.peaks.size() == 1 )
        {
          ext_lo_single.push_back( lo );
          ext_hi_single.push_back( hi );
          width_single.push_back( roi.width_fwhm );
        }
      }
      if( is_step_family( roi.family ) )
        step_z.push_back( dom.z_det );
      else if( roi.family == ContinuumFamily::Linear )
        linear_z.push_back( dom.z_det );
      else if( roi.family == ContinuumFamily::Poly2Plus )
      {
        poly_z.push_back( dom.z_det );
        poly_width.push_back( roi.width_fwhm );
        poly_energy.push_back( dom.energy );
      }
    }
  }//for( const ProblemResult &r : results )

  out << "Reference-set statistics with the canonical detection statistic z = S/sqrt(S+B),\n"
         "S = 0.9815*A (Gaussian counts within +/-1 FWHM), B = the reference continuum over the same window.\n\n";
  out << "problems=" << results.size() << " peaks=" << total_peaks << " rois=" << total_rois
      << " dont_care_peaks=" << dontcare << "\n\n";

  out << "Continuum types by ROI:\n";
  for( const auto &kv : type_rois )
    out << "  " << kv.first << ": " << kv.second << "\n";
  out << "Continuum families by ROI:\n";
  for( const auto &kv : family_rois )
    out << "  " << kv.first << ": " << kv.second << "\n";
  out << "Continuum families by peak:\n";
  for( const auto &kv : family_peaks )
    out << "  " << kv.first << ": " << kv.second << "\n";
  out << "Peaks per ROI:\n";
  for( const auto &kv : peaks_per_roi )
    out << "  " << kv.first << ": " << kv.second << "\n";

  out << "\nAdjacent reference peak pairs: separation (FWHM) -> share ROI / separate ROIs\n";
  for( size_t bin = 0; bin + 1 < sep_bins.size(); ++bin )
    out << "  [" << fmt(sep_bins[bin],3) << ", " << (sep_bins[bin+1] > 1.0e8 ? string("inf") : fmt(sep_bins[bin+1],3))
        << "): share " << share_counts[bin] << "  separate " << separate_counts[bin] << "\n";

  out << "\nROI extent beyond the outermost peak mean (FWHM of that peak):\n"
      << "  all ROIs low side:    " << quantile_str( ext_lo_all ) << "\n"
      << "  all ROIs high side:   " << quantile_str( ext_hi_all ) << "\n"
      << "  single-peak low side: " << quantile_str( ext_lo_single ) << "\n"
      << "  single-peak high side:" << quantile_str( ext_hi_single ) << "\n"
      << "  ROI width (FWHM), all:    " << quantile_str( width_all ) << "\n"
      << "  ROI width (FWHM), single: " << quantile_str( width_single ) << "\n";

  out << "\nDominant-peak z by continuum family:\n"
      << "  step:    " << quantile_str( step_z ) << "\n"
      << "  linear:  " << quantile_str( linear_z ) << "\n"
      << "  poly2+:  " << quantile_str( poly_z ) << "\n";
  out << "  step-vs-linear separability by z threshold (fraction of ROIs at or above):\n";
  for( const double thr : { 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0 } )
  {
    const size_t ns = std::count_if( begin(step_z), end(step_z), [thr]( double z ){ return z >= thr; } );
    const size_t nl = std::count_if( begin(linear_z), end(linear_z), [thr]( double z ){ return z >= thr; } );
    out << "    z>=" << fmt(thr,3) << ": step " << ns << "/" << step_z.size() << "  linear " << nl << "/" << linear_z.size() << "\n";
  }
  out << "  poly2+ ROI widths (FWHM): " << quantile_str( poly_width ) << "\n"
      << "  poly2+ dominant energies (keV): " << quantile_str( poly_energy ) << "\n";

  out << "\nAll reference peaks z: " << quantile_str( kept_z ) << "\n";
  for( const double thr : { 1.0, 1.5, 2.0, 2.5, 3.0 } )
    out << "  peaks with z<" << fmt(thr,2) << ": " << std::count_if( begin(kept_z), end(kept_z), [thr]( double z ){ return z < thr; } ) << "\n";
}//write_truth_stats

}//namespace FitPeaksCorpus
