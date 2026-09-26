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
#include <fstream>
#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <cstdio>
#include <memory>
#include <set>
#include <cctype>
#include <limits>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/RelActCalcAuto.h"

#include "FitPeaksCorpusScore.h"

using namespace std;

namespace FitPeaksCorpus
{

ContinuumFamily continuum_family( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case PeakContinuum::OffsetType::NoOffset:
    case PeakContinuum::OffsetType::Constant:
    case PeakContinuum::OffsetType::Linear:
      return ContinuumFamily::Linear;

    case PeakContinuum::OffsetType::Quadratic:
    case PeakContinuum::OffsetType::Cubic:
      return ContinuumFamily::Poly2Plus;

    case PeakContinuum::OffsetType::FlatStep:
    case PeakContinuum::OffsetType::FlatStepCDF:
      return ContinuumFamily::FlatStep;

    case PeakContinuum::OffsetType::LinearStep:
    case PeakContinuum::OffsetType::LinearStepCDF:
      return ContinuumFamily::LinearStep;

    case PeakContinuum::OffsetType::BiLinearStep:
    case PeakContinuum::OffsetType::BiLinearStepCDF:
      return ContinuumFamily::BiLinearStep;

    case PeakContinuum::OffsetType::External:
      return ContinuumFamily::Other;
  }//switch( type )

  return ContinuumFamily::Other;
}//continuum_family


const char *family_name( const ContinuumFamily family )
{
  switch( family )
  {
    case ContinuumFamily::Linear:       return "Linear";
    case ContinuumFamily::Poly2Plus:    return "Poly2Plus";
    case ContinuumFamily::FlatStep:     return "FlatStep";
    case ContinuumFamily::LinearStep:   return "LinearStep";
    case ContinuumFamily::BiLinearStep: return "BiLinearStep";
    case ContinuumFamily::Other:        return "Other";
  }
  return "Other";
}//family_name


bool is_step_family( const ContinuumFamily family )
{
  return (family == ContinuumFamily::FlatStep)
         || (family == ContinuumFamily::LinearStep)
         || (family == ContinuumFamily::BiLinearStep);
}


double gaussian_fraction_within_num_fwhm( const double num_fwhm )
{
  // +/- n FWHM = +/- n*2.35482 sigma; the enclosed Gaussian fraction is erf( x / sqrt(2) )
  return std::erf( num_fwhm * PhysicalUnits::fwhm_nsigma / std::sqrt( 2.0 ) );
}


double detection_z( const PeakDef &peak,
                    const shared_ptr<const SpecUtils::Measurement> &data,
                    const vector<shared_ptr<const PeakDef>> &roi_peaks,
                    double *continuum_counts )
{
  if( continuum_counts )
    *continuum_counts = 0.0;

  if( !peak.gausPeak() )
    return 0.0;

  const double mean = peak.mean();
  const double fwhm = peak.fwhm();
  const double signal = gaussian_fraction_within_num_fwhm( 1.0 ) * peak.amplitude();
  const double x0 = mean - fwhm;
  const double x1 = mean + fwhm;

  double continuum = 0.0;
  const shared_ptr<const PeakContinuum> cont = peak.continuum();
  if( cont && (cont->type() != PeakContinuum::OffsetType::NoOffset) && cont->parametersProbablySet() )
  {
    try
    {
      continuum = std::max( 0.0, cont->offset_integral( x0, x1, data, roi_peaks ) );
    }catch( std::exception & )
    {
      continuum = 0.0;
    }
  }else if( data )
  {
    continuum = std::max( 0.0, data->gamma_integral( static_cast<float>(x0), static_cast<float>(x1) ) - signal );
  }

  if( continuum_counts )
    *continuum_counts = continuum;

  return signal / std::sqrt( std::max( 1.0, signal + continuum ) );
}//detection_z


namespace
{
  string source_name_of( const PeakDef &peak )
  {
    if( peak.parentNuclide() )
      return peak.parentNuclide()->symbol;
    if( peak.xrayElement() )
      return peak.xrayElement()->symbol;
    if( peak.reaction() )
      return peak.reaction()->name();
    return string();
  }
}//namespace


PeakSet make_peak_set( const vector<shared_ptr<const PeakDef>> &input_peaks,
                       const shared_ptr<const SpecUtils::Measurement> &data,
                       const double truth_min_z )
{
  PeakSet answer;

  vector<shared_ptr<const PeakDef>> peaks;
  for( const shared_ptr<const PeakDef> &p : input_peaks )
  {
    if( p && p->gausPeak() )
      peaks.push_back( p );
  }

  std::sort( begin(peaks), end(peaks),
    []( const shared_ptr<const PeakDef> &a, const shared_ptr<const PeakDef> &b ){
      return a->mean() < b->mean();
  } );

  // ROI membership: peaks sharing a PeakContinuum object share a ROI (as PeakModel defines it).
  map<const void *, size_t> roi_of_continuum;
  vector<vector<shared_ptr<const PeakDef>>> roi_peak_ptrs;

  for( size_t index = 0; index < peaks.size(); ++index )
  {
    const shared_ptr<const PeakDef> &p = peaks[index];
    const shared_ptr<const PeakContinuum> cont = p->continuum();
    const void * const key = cont ? static_cast<const void *>(cont.get()) : static_cast<const void *>(p.get());

    auto pos = roi_of_continuum.find( key );
    if( pos == end(roi_of_continuum) )
    {
      ScoredRoi roi;
      roi.continuum_type = cont ? cont->type() : PeakContinuum::OffsetType::NoOffset;
      roi.family = continuum_family( roi.continuum_type );
      if( cont && cont->energyRangeDefined() )
      {
        roi.lower = cont->lowerEnergy();
        roi.upper = cont->upperEnergy();
      }else
      {
        roi.lower = p->mean() - 3.0*p->fwhm();
        roi.upper = p->mean() + 3.0*p->fwhm();
      }
      pos = roi_of_continuum.emplace( key, answer.rois.size() ).first;
      answer.rois.push_back( roi );
      roi_peak_ptrs.emplace_back();
    }

    answer.rois[pos->second].peaks.push_back( index );
    roi_peak_ptrs[pos->second].push_back( p );
  }//for( size_t index = 0; index < peaks.size(); ++index )

  for( size_t index = 0; index < peaks.size(); ++index )
  {
    const shared_ptr<const PeakDef> &p = peaks[index];
    const shared_ptr<const PeakContinuum> cont = p->continuum();
    const void * const key = cont ? static_cast<const void *>(cont.get()) : static_cast<const void *>(p.get());
    const size_t roi_index = roi_of_continuum[key];

    ScoredPeak sp;
    sp.peak = p;
    sp.energy = p->mean();
    sp.sigma = p->sigma();
    sp.fwhm = p->fwhm();
    sp.amplitude = p->amplitude();
    sp.amplitude_uncert = p->amplitudeUncert();
    sp.z_det = detection_z( *p, data, roi_peak_ptrs[roi_index], &sp.continuum_counts );
    sp.continuum_type = cont ? cont->type() : PeakContinuum::OffsetType::NoOffset;
    sp.family = continuum_family( sp.continuum_type );
    sp.roi_index = roi_index;
    sp.source_name = source_name_of( *p );
    sp.gamma_type = p->sourceGammaType();
    sp.dont_care = (truth_min_z >= 0.0) && (sp.z_det < truth_min_z);
    answer.peaks.push_back( sp );
  }//for( size_t index = 0; index < peaks.size(); ++index )

  for( ScoredRoi &roi : answer.rois )
  {
    roi.dominant = roi.peaks.front();
    for( const size_t pi : roi.peaks )
    {
      if( answer.peaks[pi].amplitude > answer.peaks[roi.dominant].amplitude )
        roi.dominant = pi;
    }
    const double dom_fwhm = answer.peaks[roi.dominant].fwhm;
    roi.width_fwhm = (dom_fwhm > 0.0) ? ((roi.upper - roi.lower) / dom_fwhm) : 0.0;

    if( data && data->num_gamma_channels() )
    {
      roi.first_channel = data->find_gamma_channel( static_cast<float>(roi.lower) );
      roi.last_channel = data->find_gamma_channel( static_cast<float>(roi.upper) );
    }
  }//for( ScoredRoi &roi : answer.rois )

  return answer;
}//make_peak_set


vector<shared_ptr<const PeakDef>> filter_peaks_by_source(
  const vector<shared_ptr<const PeakDef>> &peaks,
  const vector<string> &sources )
{
  vector<shared_ptr<const PeakDef>> answer;
  for( const shared_ptr<const PeakDef> &p : peaks )
  {
    if( !p )
      continue;
    const string name = source_name_of( *p );
    if( std::find( begin(sources), end(sources), name ) != end(sources) )
      answer.push_back( p );
  }
  return answer;
}//filter_peaks_by_source


string ScoreWeights::to_string( const string &separator ) const
{
  ostringstream out;
  out << "min_scored_energy=" << min_scored_energy << separator
      << "truth_min_z=" << truth_min_z << separator
      << "match_num_fwhm=" << match_num_fwhm << separator
      << "match_min_kev=" << match_min_kev << separator
      << "weak_z=" << weak_z << separator
      << "strong_z=" << strong_z << separator
      << "missed_weak=" << missed_weak << separator
      << "missed_moderate=" << missed_moderate << separator
      << "missed_strong=" << missed_strong << separator
      << "ghost_z=" << ghost_z << separator
      << "extra_ghost=" << extra_ghost << separator
      << "extra_weak=" << extra_weak << separator
      << "extra_significant=" << extra_significant << separator
      << "share_max_sep_fwhm=" << share_max_sep_fwhm << separator
      << "share_disagree=" << share_disagree << separator
      << "family_step_disagree=" << family_step_disagree << separator
      << "family_poly_disagree=" << family_poly_disagree << separator
      << "extent_per_fwhm=" << extent_per_fwhm << separator
      << "extent_deadband_fwhm=" << extent_deadband_fwhm << separator
      << "extent_cap_fwhm=" << extent_cap_fwhm << separator
      << "area_pull_weight=" << area_pull_weight << separator
      << "area_rel_floor=" << area_rel_floor << separator
      << "area_pull_cap=" << area_pull_cap << separator
      << "mean_offset_weight=" << mean_offset_weight << separator
      << "mean_offset_deadband_fwhm=" << mean_offset_deadband_fwhm << separator
      << "mean_offset_cap_fwhm=" << mean_offset_cap_fwhm << separator
      << "failure=" << failure << separator
      << "nondeterminism=" << nondeterminism << separator
      << "legit_truth_min_z=" << legit_truth_min_z << separator
      << "truth_area_weight=" << truth_area_weight << separator
      << "truth_area_rel_floor=" << truth_area_rel_floor << separator
      << "truth_area_cap=" << truth_area_cap << separator
      << "truth_area_min_z=" << truth_area_min_z << separator
      << "legit_max_pull=" << legit_max_pull << separator
      << "merged_max_sep_fwhm=" << merged_max_sep_fwhm;
  return out.str();
}//ScoreWeights::to_string


double ProblemScore::raw_cost() const
{
  return cost_missed + cost_extra + cost_share + cost_family + cost_extent
         + cost_area + cost_mean + cost_failure + cost_truth_area;
}


double ProblemScore::norm_cost() const
{
  return raw_cost() / static_cast<double>( std::max( size_t(1), n_truth_scored ) );
}


namespace
{
  /** Index of the merged truth photopeak whose energy range (widened by the match window) holds
   `energy`, choosing the nearest when several do; -1 when none. */
  int nearest_truth_cluster( const vector<TruthPhotopeak> &clusters, const double energy,
                             const double fwhm, const ScoreWeights &w )
  {
    int best = -1;
    double best_dist = 0.0;
    for( size_t i = 0; i < clusters.size(); ++i )
    {
      const TruthPhotopeak &c = clusters[i];
      const double window = std::max( w.match_num_fwhm * std::max( c.fwhm, fwhm ), w.match_min_kev );
      const double dist = (energy < c.energy_lo) ? (c.energy_lo - energy)
                        : ((energy > c.energy_hi) ? (energy - c.energy_hi) : 0.0);
      if( (dist <= window) && ((best < 0) || (dist < best_dist)) )
      {
        best = static_cast<int>( i );
        best_dist = dist;
      }
    }
    return best;
  }
}//namespace


ProblemScore score_problem( PeakSet &truth, PeakSet &fitted, const ScoreWeights &w,
                            const InjectTruth *inject, const vector<double> *source_xrays,
                            const shared_ptr<const SpecUtils::Measurement> &foreground,
                            const vector<double> *source_lines )
{
  ProblemScore score;

  // Net counts above a sideband-median continuum over the peak core, as a Poisson significance.
  // Deliberately independent of both peak sets: it asks only whether the spectrum has an excess.
  const auto data_excess_z = [&foreground]( const double energy, const double fwhm ) -> double {
    if( !foreground || !foreground->num_gamma_channels() || !foreground->channel_energies()
        || !(fwhm > 0.0) )
      return std::numeric_limits<double>::infinity();   // no spectrum: never excuse a miss
    const vector<float> &counts = *foreground->gamma_counts();
    vector<float> lo_side, hi_side;
    double gross = 0.0;
    size_t nchan_core = 0;
    for( size_t ch = 0; ch < counts.size(); ++ch )
    {
      const double e = 0.5*( foreground->gamma_channel_lower(ch) + foreground->gamma_channel_upper(ch) );
      const double d = e - energy;
      if( std::fabs(d) <= fwhm ){ gross += counts[ch]; ++nchan_core; }
      else if( (d < -1.5*fwhm) && (d >= -3.0*fwhm) ) lo_side.push_back( counts[ch] );
      else if( (d > 1.5*fwhm) && (d <= 3.0*fwhm) ) hi_side.push_back( counts[ch] );
    }
    if( !nchan_core || (lo_side.size() < 2) || (hi_side.size() < 2) )
      return std::numeric_limits<double>::infinity();
    // Take the LOWER of the two side medians.  A single median over both sides sits on whatever
    // neighbours crowd the window, which in an x-ray multiplet cancels the very excess being tested
    // (it excused a 2678-count Np237 94.7 keV reference peak).  In a multiplet at least one side is
    // usually clean, and erring toward a larger net keeps a real peak scored as a miss.
    std::sort( begin(lo_side), end(lo_side) );
    std::sort( begin(hi_side), end(hi_side) );
    const double level = std::min( static_cast<double>( lo_side[lo_side.size()/2] ),
                                   static_cast<double>( hi_side[hi_side.size()/2] ) );
    const double net = gross - level*static_cast<double>( nchan_core );
    return net / std::sqrt( std::max( 1.0, gross ) );
  };

  // Below min_scored_energy nothing is judged: the planner has a physical low-energy floor, and a
  // detector that sees below it (a planar HPGe) would otherwise be charged for every line there.
  for( ScoredPeak &p : truth.peaks )
  {
    if( p.energy < w.min_scored_energy )
      p.dont_care = true;
  }
  for( ScoredPeak &p : fitted.peaks )
  {
    if( p.energy < w.min_scored_energy )
      p.verdict = "below_range";
  }

  for( ScoredPeak &p : truth.peaks )
  {
    p.match = -1;
    p.truth_cluster = -1;
    p.verdict = p.dont_care ? "dontcare" : "missed";
    // See ScoreWeights::miss_requires_data_z - a truth peak absent from the spectrum is not a miss.
    if( !p.dont_care && (w.miss_requires_data_z > 0.0)
        && (data_excess_z( p.energy, p.fwhm ) < w.miss_requires_data_z) )
    {
      p.dont_care = true;
      p.verdict = "dontcare_nodata";
      score.truth_absent_from_data += 1;
    }

    // See ScoreWeights::shield_fluorescence_num_fwhm - a shield x-ray no requested source can emit.
    if( !p.dont_care && (w.shield_fluorescence_num_fwhm > 0.0) && (p.fwhm > 0.0) )
    {
      // Pb K series; the shielding in these inject geometries is lead.
      static const double sm_shield_xrays[] = { 72.805, 74.969, 84.450, 84.936, 87.360 };
      const double window = w.shield_fluorescence_num_fwhm * p.fwhm;
      bool on_shield_xray = false;
      for( const double e : sm_shield_xrays )
        on_shield_xray = on_shield_xray || (std::fabs( e - p.energy ) <= window);

      bool source_explains = false;
      if( source_lines )
      {
        for( const double e : *source_lines )
          source_explains = source_explains || (std::fabs( e - p.energy ) <= window);
      }

      if( on_shield_xray && !source_explains )
      {
        p.dont_care = true;
        p.verdict = "dontcare_shield_xray";
        score.truth_shield_xray += 1;
      }
    }
    // A reference peak on a background photopeak where the source has none is a mislabelled
    // background line (a "Pu239 338.1 keV" peak on Ac228 338.3 keV): neither missed nor matched.
    // Requires a real signal list: ~12 % of the corpus attaches a truth whose signal section was
    // never computed, and without this guard the "source has none" test is vacuously true there.
    if( inject && !p.dont_care && !inject->merged_signal.empty()
        && (nearest_truth_cluster( inject->merged_signal, p.energy, p.fwhm, w ) < 0)
        && (nearest_truth_cluster( inject->merged_background, p.energy, p.fwhm, w ) >= 0) )
    {
      p.dont_care = true;
      p.verdict = "dontcare_bkg";
      score.truth_on_bkg_line += 1;
    }
  }
  for( ScoredPeak &p : fitted.peaks )
  {
    p.match = -1;
    p.truth_cluster = -1;
    p.verdict.clear();
  }

  // Greedy matching: scored truth peaks by descending amplitude first, then don't-care peaks.
  vector<size_t> truth_order;
  for( size_t i = 0; i < truth.peaks.size(); ++i )
    truth_order.push_back( i );
  std::stable_sort( begin(truth_order), end(truth_order), [&truth]( const size_t a, const size_t b ){
    const ScoredPeak &pa = truth.peaks[a];
    const ScoredPeak &pb = truth.peaks[b];
    if( pa.dont_care != pb.dont_care )
      return !pa.dont_care;
    return pa.amplitude > pb.amplitude;
  } );

  for( const size_t ti : truth_order )
  {
    ScoredPeak &t = truth.peaks[ti];
    const double window = std::max( w.match_num_fwhm * t.fwhm, w.match_min_kev );

    int best = -1;
    bool best_same_source = false;
    double best_offset = 0.0;
    for( size_t fi = 0; fi < fitted.peaks.size(); ++fi )
    {
      const ScoredPeak &f = fitted.peaks[fi];
      if( f.match >= 0 )
        continue;
      const double offset = std::fabs( f.energy - t.energy );
      if( offset > window )
        continue;
      const bool same_source = !t.source_name.empty() && (t.source_name == f.source_name);
      if( (best < 0)
          || (same_source && !best_same_source)
          || ((same_source == best_same_source) && (offset < best_offset)) )
      {
        best = static_cast<int>( fi );
        best_same_source = same_source;
        best_offset = offset;
      }
    }//for( size_t fi = 0; fi < fitted.peaks.size(); ++fi )

    if( best < 0 )
      continue;

    t.match = best;
    fitted.peaks[best].match = static_cast<int>( ti );
    if( t.dont_care )
    {
      t.verdict = "matched_dontcare";
      fitted.peaks[best].verdict = "neutral";
    }else
    {
      t.verdict = "matched";
      fitted.peaks[best].verdict = "matched";
    }
  }//for( const size_t ti : truth_order )

  // Truth-photopeak membership, needed before the presence terms: a reference peak the truth merges
  // into a photopeak the fit did report is not a miss (the truth says there is one detectable peak
  // there), it is the reference having split one peak in two.
  if( inject && !inject->merged_signal.empty() )
  {
    for( ScoredPeak &f : fitted.peaks )
      f.truth_cluster = nearest_truth_cluster( inject->merged_signal, f.energy, f.fwhm, w );
    for( ScoredPeak &t : truth.peaks )
      t.truth_cluster = nearest_truth_cluster( inject->merged_signal, t.energy, t.fwhm, w );

    std::set<int> clusters_with_a_fitted_peak;
    for( const ScoredPeak &f : fitted.peaks )
    {
      // A ghost must not buy an excuse: it costs 0.1 and could forgive a 4.0 miss.
      if( (f.truth_cluster >= 0) && (f.z_det >= w.ghost_z) )
        clusters_with_a_fitted_peak.insert( f.truth_cluster );
    }
    for( ScoredPeak &t : truth.peaks )
    {
      if( t.dont_care || (t.match >= 0) || (t.truth_cluster < 0)
          || !clusters_with_a_fitted_peak.count( t.truth_cluster ) )
        continue;
      // ... and the peak must be genuinely unresolvable from one the fit DID find: a matched
      // reference peak within merged_max_sep_fwhm.  The truth's merge is single-linkage and chains
      // (3.1 FWHM on this HPGe, 5.5 on NaI), so "same cluster" alone would excuse peaks a good fit
      // ought to separate.
      bool unresolvable = false;
      for( const ScoredPeak &o : truth.peaks )
      {
        if( (o.match < 0) || (o.truth_cluster != t.truth_cluster) || (t.fwhm <= 0.0) )
          continue;
        if( std::fabs( o.energy - t.energy ) <= (w.merged_max_sep_fwhm * t.fwhm) )
          unresolvable = true;
      }
      if( !unresolvable )
        continue;
      t.dont_care = true;
      t.verdict = "dontcare_merged";
      score.truth_merged_away += 1;
    }
  }//if( inject truth )

  // Presence terms
  for( const ScoredPeak &t : truth.peaks )
  {
    if( t.dont_care )
    {
      score.n_truth_dontcare += 1;
      continue;
    }
    score.n_truth_scored += 1;
    if( t.match >= 0 )
    {
      score.n_matched += 1;
      continue;
    }
    // Class a miss by the smallest significance any source claims for it, INCLUDING what the
    // spectrum itself shows.  A truth file computed from the source term rather than the recording
    // can be badly optimistic - Pd103's 39.76 keV line is claimed at 3869 counts (z=62) where the
    // Fulcrum data holds 36 counts (z~6) - and grading a miss by the claim rather than by the
    // evidence makes a marginal peak look like a headline failure.
    double class_z = t.z_det;
    if( inject )
    {
      const int c = nearest_truth_cluster( inject->merged_signal, t.energy, t.fwhm, w );
      if( c >= 0 )
        class_z = std::min( class_z, inject->merged_signal[c].z_det() );
    }
    if( w.class_miss_by_data_z )
    {
      const double observed_z = data_excess_z( t.energy, t.fwhm );
      if( std::isfinite(observed_z) )
        class_z = std::min( class_z, observed_z );
    }
    if( class_z >= w.strong_z )
    {
      score.missed_strong += 1;
      score.cost_missed += w.missed_strong;
    }else if( class_z >= w.weak_z )
    {
      score.missed_moderate += 1;
      score.cost_missed += w.missed_moderate;
    }else
    {
      score.missed_weak += 1;
      score.cost_missed += w.missed_weak;
    }
  }//for( const ScoredPeak &t : truth.peaks )

  score.n_fitted = fitted.peaks.size();
  for( ScoredPeak &f : fitted.peaks )
  {
    if( f.energy < w.min_scored_energy )
    {
      f.verdict = "below_range";
      continue;
    }
    if( f.match >= 0 )
    {
      if( f.verdict == "neutral" )
        score.neutral += 1;
      continue;
    }

    // Inject truth: a peak the reference did not fit is still right when the source really emits
    // a visible photopeak there (the hand fits skipped some weak lines and some 511 keV peaks).
    // A peak that also sits on a background photopeak is never "legit", whatever the source line
    // beside it: absorbing a background line into a source-labelled peak is the very thing
    // extra_bkg exists to catch, and the legit area guard scales with the background it swallowed.
    const bool on_background = inject && !inject->merged_background.empty()
                               && (nearest_truth_cluster( inject->merged_background, f.energy, f.fwhm, w ) >= 0);
    if( inject && !on_background && (f.z_det >= w.ghost_z) )
    {
      const int c = nearest_truth_cluster( inject->merged_signal, f.energy, f.fwhm, w );
      if( (c >= 0) && (inject->merged_signal[c].z_det() >= w.legit_truth_min_z) )
      {
        // ... provided the fitted area is what the line can supply; a 160-count peak on a
        // 20-count line is the background peak underneath it, wearing the source's label.
        const TruthPhotopeak &t = inject->merged_signal[c];
        const double sigma = std::sqrt( std::max( 1.0, t.area + t.continuum_within_fwhm() )
                                        + std::pow( w.truth_area_rel_floor * t.area, 2.0 ) );
        if( ((f.amplitude - t.area) / sigma) <= w.legit_max_pull )
        {
          f.verdict = "legit";
          score.extra_legit += 1;
          continue;
        }
      }
    }

    // Annihilation: real, present in essentially every positron-emitter and pair-production
    // spectrum, and absent from every GADRAS photopeak list, so it would otherwise be charged as an
    // extra on every such source (Br76's 511 keV peak is 11700 counts at z=107).  Free, like "xray".
    if( w.free_annihilation_peak && !on_background && (f.fwhm > 0.0)
        && (std::fabs( f.energy - 510.9989 ) <= std::max( w.match_num_fwhm * f.fwhm, w.match_min_kev )) )
    {
      f.verdict = "annih";
      score.fitted_annihilation += 1;
      continue;
    }

    // A characteristic x-ray of a requested source is a real peak, whatever the reference and the
    // GADRAS photopeak list say (neither enumerates x-rays).  Free, like "legit".
    if( source_xrays && !on_background && (f.z_det >= w.ghost_z) && (f.fwhm > 0.0) )
    {
      bool on_xray = false;
      for( const double e : *source_xrays )
        on_xray = on_xray || (std::fabs( e - f.energy ) <= std::max( w.match_num_fwhm * f.fwhm, w.match_min_kev ));
      if( on_xray )
      {
        f.verdict = "xray";
        score.extra_xray += 1;
        continue;
      }
    }

    if( f.z_det < w.ghost_z )
    {
      f.verdict = "ghost";
      score.extra_ghost += 1;
      score.cost_extra += w.extra_ghost;
    }else if( f.z_det < w.weak_z )
    {
      f.verdict = "extra_weak";
      score.extra_weak += 1;
      score.cost_extra += w.extra_weak;
    }else
    {
      f.verdict = "extra";
      score.extra_significant += 1;
      score.cost_extra += w.extra_significant;
    }

    // A real background photopeak carrying the requested source's label: costed as the extra it is,
    // but counted separately so mislabels can be told from invented peaks.
    if( on_background && (f.verdict != "ghost") )
    {
      f.verdict = "extra_bkg";
      score.extra_bkg_line += 1;
    }
  }//for( ScoredPeak &f : fitted.peaks )

  // Truth-area grading: pool the fitted and the reference peaks on each merged signal photopeak
  // and compare the pooled areas with the expected counts.
  if( inject && !inject->merged_signal.empty() )
  {
    const vector<TruthPhotopeak> &clusters = inject->merged_signal;
    score.truth_clusters = clusters.size();
    vector<TruthAreaRecord> records( clusters.size() );
    for( size_t c = 0; c < clusters.size(); ++c )
    {
      const TruthPhotopeak &t = clusters[c];
      records[c].cluster = c;
      records[c].energy = t.energy;
      records[c].area = t.area;
      records[c].z = t.z_det();
      records[c].sigma = std::sqrt( std::max( 1.0, t.area + t.continuum_within_fwhm() )
                                    + std::pow( w.truth_area_rel_floor * t.area, 2.0 ) );
      records[c].fit_pull = records[c].ref_pull = std::numeric_limits<double>::quiet_NaN();
    }

    for( const ScoredPeak &f : fitted.peaks )
    {
      if( f.truth_cluster >= 0 )
      {
        records[f.truth_cluster].fit_n += 1;
        records[f.truth_cluster].fit_area += f.amplitude;
      }
    }
    for( const ScoredPeak &t : truth.peaks )
    {
      if( t.truth_cluster >= 0 )
      {
        records[t.truth_cluster].ref_n += 1;
        records[t.truth_cluster].ref_area += t.amplitude;
      }
    }

    // Peaks that x-rays contribute to are not graded on area: the decay data's x-ray yields are
    // unreliable (for I123 the Te K-alpha lines come back 22,000x below their published per-decay
    // intensity, while the same mixture's gammas are exact), and a peak may in any case take a
    // fluorescence contribution the source model does not carry.
    for( TruthAreaRecord &rec : records )
    {
      if( source_xrays )
      {
        const double win = std::max( w.match_num_fwhm * clusters[rec.cluster].fwhm, w.match_min_kev );
        for( const double e : *source_xrays )
          rec.xray_affected = rec.xray_affected || (std::fabs( e - rec.energy ) <= win);
      }
      if( rec.xray_affected )
      {
        score.truth_xray_clusters += 1;
        continue;
      }
      const bool strong_enough = (rec.z >= w.truth_area_min_z);
      rec.fit_graded = strong_enough || (rec.fit_n > 0);
      rec.ref_graded = strong_enough || (rec.ref_n > 0);
      if( rec.fit_graded )
      {
        rec.fit_pull = (rec.fit_area - rec.area) / rec.sigma;
        const double capped = std::min( std::fabs( rec.fit_pull ), w.truth_area_cap );
        score.truth_fit_graded += 1;
        score.truth_fit_abs_pull += capped;
        score.truth_fit_bad += (std::fabs( rec.fit_pull ) > 3.0) ? 1 : 0;
        score.cost_truth_area += w.truth_area_weight * capped;
      }
      if( rec.ref_graded )
      {
        rec.ref_pull = (rec.ref_area - rec.area) / rec.sigma;
        score.truth_ref_graded += 1;
        score.truth_ref_abs_pull += std::min( std::fabs( rec.ref_pull ), w.truth_area_cap );
        score.truth_ref_bad += (std::fabs( rec.ref_pull ) > 3.0) ? 1 : 0;
      }
    }
    score.truth_area_records = std::move( records );
  }//if( inject truth )

  // Share/separate agreement on adjacent matched scored truth peaks (truth peaks are energy sorted)
  vector<size_t> matched_truth;
  for( size_t ti = 0; ti < truth.peaks.size(); ++ti )
  {
    const ScoredPeak &t = truth.peaks[ti];
    if( !t.dont_care && (t.match >= 0) )
      matched_truth.push_back( ti );
  }

  for( size_t k = 1; k < matched_truth.size(); ++k )
  {
    const ScoredPeak &a = truth.peaks[matched_truth[k-1]];
    const ScoredPeak &b = truth.peaks[matched_truth[k]];
    const double mid_fwhm = 0.5*(a.fwhm + b.fwhm);
    if( mid_fwhm <= 0.0 )
      continue;
    const double sep = (b.energy - a.energy) / mid_fwhm;
    if( sep > w.share_max_sep_fwhm )
      continue;

    PairRecord rec;
    rec.truth_a = matched_truth[k-1];
    rec.truth_b = matched_truth[k];
    rec.separation_fwhm = sep;
    rec.truth_share = (a.roi_index == b.roi_index);
    rec.fit_share = (fitted.peaks[a.match].roi_index == fitted.peaks[b.match].roi_index);
    score.pairs += 1;
    if( rec.truth_share != rec.fit_share )
    {
      score.share_disagree += 1;
      score.cost_share += w.share_disagree;
    }
    score.pair_records.push_back( rec );
  }//for( size_t k = 1; k < matched_truth.size(); ++k )

  // Continuum family and extent, per truth ROI whose dominant scored peak is matched
  for( size_t ri = 0; ri < truth.rois.size(); ++ri )
  {
    const ScoredRoi &troi = truth.rois[ri];
    int dominant = -1;
    for( const size_t pi : troi.peaks )
    {
      const ScoredPeak &p = truth.peaks[pi];
      if( p.dont_care )
        continue;
      if( (dominant < 0) || (p.amplitude > truth.peaks[dominant].amplitude) )
        dominant = static_cast<int>( pi );
    }
    if( (dominant < 0) || (truth.peaks[dominant].match < 0) )
      continue;

    const ScoredPeak &dom = truth.peaks[dominant];
    const ScoredPeak &fdom = fitted.peaks[dom.match];
    const ScoredRoi &froi = fitted.rois[fdom.roi_index];

    RoiRecord rec;
    rec.truth_roi = ri;
    rec.fit_roi = static_cast<int>( fdom.roi_index );
    rec.truth_family = troi.family;
    rec.fit_family = froi.family;
    score.rois_compared += 1;

    if( troi.family != froi.family )
    {
      score.family_disagree += 1;
      const bool step_involved = is_step_family( troi.family ) || is_step_family( froi.family );
      rec.cost_family = step_involved ? w.family_step_disagree : w.family_poly_disagree;
      score.cost_family += rec.cost_family;
    }

    const double fwhm = dom.fwhm;
    if( fwhm > 0.0 )
    {
      rec.d_lower_fwhm = (froi.lower - troi.lower) / fwhm;
      rec.d_upper_fwhm = (froi.upper - troi.upper) / fwhm;
      for( const double d : { rec.d_lower_fwhm, rec.d_upper_fwhm } )
      {
        score.extent_sides += 1;
        if( std::fabs(d) > 1.0 )
          score.extent_gt1fwhm_sides += 1;
        const double excess = std::min( w.extent_cap_fwhm, std::max( 0.0, std::fabs(d) - w.extent_deadband_fwhm ) );
        rec.cost_extent += w.extent_per_fwhm * excess;
      }
      score.cost_extent += rec.cost_extent;
    }

    score.roi_records.push_back( rec );
  }//for( size_t ri = 0; ri < truth.rois.size(); ++ri )

  // Parameter agreement on matched scored peaks
  for( const ScoredPeak &t : truth.peaks )
  {
    if( t.dont_care || (t.match < 0) )
      continue;
    const ScoredPeak &f = fitted.peaks[t.match];

    const double stat_uncert = (t.amplitude_uncert > 0.0) ? t.amplitude_uncert : std::sqrt( std::max( 1.0, t.amplitude ) );
    const double sigma_area = std::sqrt( stat_uncert*stat_uncert + std::pow( w.area_rel_floor * t.amplitude, 2.0 ) );
    const double pull = std::fabs( f.amplitude - t.amplitude ) / std::max( 1.0e-9, sigma_area );
    score.cost_area += w.area_pull_weight * std::min( pull, w.area_pull_cap );

    if( t.fwhm > 0.0 )
    {
      const double offset = std::fabs( f.energy - t.energy ) / t.fwhm;
      const double excess = std::min( w.mean_offset_cap_fwhm, std::max( 0.0, offset - w.mean_offset_deadband_fwhm ) );
      score.cost_mean += w.mean_offset_weight * excess;
    }
  }//for( const ScoredPeak &t : truth.peaks )

  return score;
}//score_problem


string problem_id_from_path( const string &path )
{
  string name = SpecUtils::filename( path );
  const size_t dot = name.rfind( '.' );
  if( dot != string::npos )
    name = name.substr( 0, dot );
  return name;
}


vector<string> split_source_list( const string &text )
{
  vector<string> parts;
  SpecUtils::split( parts, text, ",;" );
  vector<string> answer;
  for( string &p : parts )
  {
    SpecUtils::trim( p );
    if( !p.empty() )
      answer.push_back( p );
  }
  return answer;
}


double TruthPhotopeak::continuum_within_fwhm() const
{
  const double width = roi_upper - roi_lower;
  const double fraction = (width > 0.0) ? std::min( 1.0, (2.0 * fwhm) / width ) : 1.0;
  return std::max( 0.0, continuum_area ) * fraction;
}


double TruthPhotopeak::z_det() const
{
  const double signal = gaussian_fraction_within_num_fwhm( 1.0 ) * area;
  return signal / std::sqrt( std::max( 1.0, signal + continuum_within_fwhm() ) );
}


InjectTruth parse_inject_truth_csv( const string &path )
{
  InjectTruth truth;
  truth.path = path;

  ifstream input( path.c_str() );
  if( !input )
    throw runtime_error( "Could not open '" + path + "'" );

  string line;
  int section = 0;   // 1 = signal, 2 = background, 0/other = not a section we keep
  bool have_signal_section = false, past_header = false;
  while( SpecUtils::safe_get_line( input, line, 4096 ) )
  {
    SpecUtils::trim( line );
    if( line.empty() )
      continue;
    if( SpecUtils::istarts_with( line, "# Source:" ) )
    {
      truth.source_name = SpecUtils::trim_copy( line.substr( strlen("# Source:") ) );
      continue;
    }
    if( line.back() == ':' )
    {
      if( SpecUtils::istarts_with( line, "Expected Signal Photopeaks:" ) )
      {
        section = 1;
        have_signal_section = true;
      }else if( SpecUtils::istarts_with( line, "Expected Background Photopeaks:" ) )
      {
        section = 2;
      }else
      {
        section = 0;
      }
      past_header = false;
      continue;
    }
    if( section == 0 )
      continue;
    if( !past_header )
    {
      past_header = true;   // the column-name row
      continue;
    }
    vector<string> fields;
    SpecUtils::split( fields, line, "," );
    if( fields.size() != 8 )
      continue;   // per-gamma-line rows carry six fields
    TruthPhotopeak p;
    double num_lines = 0.0;
    if( !SpecUtils::parse_double( fields[0].c_str(), fields[0].size(), p.nsigma )
       || !SpecUtils::parse_double( fields[1].c_str(), fields[1].size(), p.roi_lower )
       || !SpecUtils::parse_double( fields[2].c_str(), fields[2].size(), p.roi_upper )
       || !SpecUtils::parse_double( fields[3].c_str(), fields[3].size(), p.energy )
       || !SpecUtils::parse_double( fields[4].c_str(), fields[4].size(), p.area )
       || !SpecUtils::parse_double( fields[5].c_str(), fields[5].size(), p.continuum_area )
       || !SpecUtils::parse_double( fields[6].c_str(), fields[6].size(), p.fwhm )
       || !SpecUtils::parse_double( fields[7].c_str(), fields[7].size(), num_lines ) )
      continue;
    if( !(p.fwhm > 0.0) || !(p.area > 0.0) || !(p.roi_upper > p.roi_lower) )
      continue;
    p.energy_lo = p.energy_hi = p.energy;
    p.num_lines = static_cast<size_t>( std::max( 0.0, num_lines ) );
    p.signal = (section == 1);
    (p.signal ? truth.signal : truth.background).push_back( p );
  }//while( lines )

  if( !have_signal_section )
    throw runtime_error( "no 'Expected Signal Photopeaks' section in '" + path + "'" );

  const auto by_energy = []( const TruthPhotopeak &a, const TruthPhotopeak &b ){ return a.energy < b.energy; };
  std::stable_sort( begin(truth.signal), end(truth.signal), by_energy );
  std::stable_sort( begin(truth.background), end(truth.background), by_energy );
  truth.merged_signal = merge_unresolved_photopeaks( truth.signal );
  truth.merged_background = merge_unresolved_photopeaks( truth.background );
  return truth;
}//parse_inject_truth_csv


vector<TruthPhotopeak> merge_unresolved_photopeaks( const vector<TruthPhotopeak> &rows )
{
  // The truth lists every emitted photopeak, but two lines closer than about a FWHM are one peak in
  // the data and one peak in any fit - scoring them separately counts a perfect fit as missing half
  // of them (a quarter of the Eu152 rows on a NaI detector).  Merge into the previous entry when
  // they overlap: areas add, the energy is area-weighted, the ROI is the union.
  vector<TruthPhotopeak> merged;
  for( const TruthPhotopeak &row : rows )
  {
    if( !merged.empty() )
    {
      TruthPhotopeak &prev = merged.back();
      const double sep = row.energy - prev.energy_hi;
      const double mean_fwhm = 0.5*(row.fwhm + prev.fwhm);
      if( (sep >= 0.0) && (sep < mean_fwhm) )
      {
        const double total = prev.area + row.area;
        prev.energy = (total > 0.0) ? ((prev.area*prev.energy + row.area*row.energy) / total) : row.energy;
        prev.energy_hi = row.energy;
        prev.area = total;
        prev.fwhm = std::max( prev.fwhm, row.fwhm );
        prev.roi_lower = std::min( prev.roi_lower, row.roi_lower );
        prev.roi_upper = std::max( prev.roi_upper, row.roi_upper );
        prev.continuum_area = std::max( prev.continuum_area, row.continuum_area );
        prev.nsigma = std::max( prev.nsigma, row.nsigma );
        prev.num_lines += row.num_lines;
        prev.num_rows += 1;
        continue;
      }
    }
    merged.push_back( row );
  }
  return merged;
}//merge_unresolved_photopeaks


string inject_name_for_problem_id( const string &problem_id )
{
  // "q115-Lu177m-Unsh" -> "Lu177m-Unsh" -> "Lu177m_Unsh"; "q098-K40-Sh-Point" -> "K40_Sh-Point"
  string name = problem_id;
  if( (name.size() > 1) && (name[0] == 'q') && isdigit( static_cast<unsigned char>(name[1]) ) )
  {
    const size_t dash = name.find( '-' );
    if( dash != string::npos )
      name = name.substr( dash + 1 );
  }
  const size_t first = name.find( '-' );
  if( first != string::npos )
    name[first] = '_';
  return name;
}//inject_name_for_problem_id


bool attach_inject_truth( CorpusProblem &problem, const string &inject_dir )
{
  const string name = inject_name_for_problem_id( problem.id );
  const string truth_path = SpecUtils::append_path( inject_dir, name + "_truth.csv" );
  if( !SpecUtils::is_file( truth_path ) )
    return false;

  auto truth = make_shared<InjectTruth>( parse_inject_truth_csv( truth_path ) );

  // The truth only applies if this really is the same spectrum: compare against the PCF's record 0.
  const string pcf_path = SpecUtils::append_path( inject_dir, name + ".pcf" );
  string verification;
  try
  {
    SpecMeas meas;
    if( !SpecUtils::is_file( pcf_path ) || !meas.load_file( pcf_path, SpecUtils::ParserType::Pcf )
        || (meas.num_measurements() < 1) )
      throw runtime_error( "could not load " + pcf_path );
    const shared_ptr<const SpecUtils::Measurement> fg = meas.measurement_at_index( 0 );
    const shared_ptr<const vector<float>> a = fg ? fg->gamma_counts() : nullptr;
    const shared_ptr<const vector<float>> b = problem.foreground ? problem.foreground->gamma_counts() : nullptr;
    if( !a || !b )
      throw runtime_error( "no channel data to compare" );
    if( a->size() != b->size() )
      throw runtime_error( "channel count differs (" + std::to_string(a->size()) + " vs " + std::to_string(b->size()) + ")" );
    size_t differing = 0;
    for( size_t i = 0; i < a->size(); ++i )
      differing += (std::fabs( (*a)[i] - (*b)[i] ) > 1.0e-3f * std::max( 1.0f, std::fabs( (*b)[i] ) )) ? 1 : 0;
    if( differing )
      throw runtime_error( std::to_string(differing) + " of " + std::to_string(a->size()) + " channels differ" );
  }catch( std::exception &e )
  {
    verification = e.what();
  }

  problem.inject_truth = truth;
  if( verification.empty() )
    problem.notes.push_back( "inject truth " + name + " (spectrum verified)" );
  else
    problem.notes.push_back( "inject truth " + name + " NOT VERIFIED: " + verification );
  return true;
}//attach_inject_truth


vector<double> source_xray_energies( const vector<RelActCalcAuto::SrcVariant> &sources,
                                     const double min_intensity )
{
  vector<double> answer;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    return answer;

  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    vector<SandiaDecay::EnergyRatePair> xrays;
    const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( src );
    if( nuc )
    {
      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( nuc, 1.0, 0.0 );
      xrays = mix.xrays( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy );
    }else
    {
      const SandiaDecay::Element * const el = RelActCalcAuto::element( src );
      if( el )
      {
        for( const SandiaDecay::EnergyIntensityPair &p : el->xrays )
          xrays.push_back( SandiaDecay::EnergyRatePair( p.intensity, p.energy ) );
      }
    }

    double strongest = 0.0;
    for( const SandiaDecay::EnergyRatePair &p : xrays )
      strongest = std::max( strongest, p.numPerSecond );
    for( const SandiaDecay::EnergyRatePair &p : xrays )
    {
      if( (strongest > 0.0) && (p.numPerSecond >= min_intensity*strongest) )
        answer.push_back( p.energy );
    }
  }//for( const RelActCalcAuto::SrcVariant &src : sources )

  std::sort( begin(answer), end(answer) );
  answer.erase( std::unique( begin(answer), end(answer) ), end(answer) );
  return answer;
}//source_xray_energies


vector<double> source_line_energies( const vector<RelActCalcAuto::SrcVariant> &sources,
                                     const double min_intensity )
{
  vector<double> answer = source_xray_energies( sources, min_intensity );
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    return answer;

  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( src );
    if( !nuc )
      continue;
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuc, 1.0 );
    const double age = PeakDef::defaultDecayTime( nuc, nullptr );
    const vector<SandiaDecay::EnergyRatePair> gammas
        = mix.gammas( age, SandiaDecay::NuclideMixture::OrderByEnergy, true );
    double strongest = 0.0;
    for( const SandiaDecay::EnergyRatePair &p : gammas )
      strongest = std::max( strongest, p.numPerSecond );
    for( const SandiaDecay::EnergyRatePair &p : gammas )
    {
      if( (strongest > 0.0) && (p.numPerSecond >= min_intensity*strongest) )
        answer.push_back( p.energy );
    }
  }//for( const RelActCalcAuto::SrcVariant &src : sources )

  std::sort( begin(answer), end(answer) );
  answer.erase( std::unique( begin(answer), end(answer) ), end(answer) );
  return answer;
}//source_line_energies


CorpusProblem load_inject_problem( const string &truth_csv_path )
{
  CorpusProblem problem;
  const string base = SpecUtils::filename( truth_csv_path );
  const string suffix = "_truth.csv";
  if( (base.size() <= suffix.size()) || !SpecUtils::iends_with( base, suffix ) )
    throw runtime_error( "Not an inject truth file: '" + truth_csv_path + "'" );
  problem.id = base.substr( 0, base.size() - suffix.size() );
  problem.path = SpecUtils::append_path( SpecUtils::parent_path( truth_csv_path ), problem.id + ".pcf" );

  // Requested source(s).  The file name carries one nuclide and a shielding suffix, but these
  // sources are mixtures; the mapping below is the user's, determined by manual inspection of the
  // set (target/peak_fit_improve/FitPeaksForNuclideDev.cpp).  Asking only for the named nuclide
  // makes several problems unanswerable - U233 without its U232 daughters cannot explain the
  // Pb212/Tl208 series that dominates its spectrum - which scores as a fitter failure when it is
  // really a badly posed question.
  const string token = problem.id.substr( 0, problem.id.find( '_' ) );
  string sources_text;
  if( (token == "Tl201woTl202") || (token == "Tl201wTl202") || (token == "Tl201") )
    sources_text = "Tl201,Tl202";
  else if( token == "I125" )
    sources_text = "I125,I126";
  // Uranium and plutonium items are always isotope mixtures, and the metal fluoresces its own K
  // x-rays, so the element goes in beside the isotopes (user, 2026-09-06).  Ages need no special
  // handling: PeakDef::defaultDecayTime already gives every U/Pu isotope with a half-life over two
  // years the 20 years these sources are assumed to have.
  else if( token == "U233" )
    sources_text = "U,U232,U233,U234,U235,U238";   // U233 fuel, plus the standard uranium suite
  else if( (token == "Pu238") || (token == "Pu239") )
    sources_text = "Pu,Pu238,Pu239,Pu240,Pu241";
  else if( token == "Uore" )
    sources_text = "U,U232,U234,U235,U238,Ra226";
  else if( (token == "U235") || (token == "U238") )
    sources_text = "U,U232,U234,U235,U238";
  else if( token == "Np237" )
    sources_text = "Np,Np237";
  else if( token == "Xe133" )
    sources_text = "Xe133,Xe133m";
  else if( (token == "Cf252") || (token == "Am241Li") )
    throw runtime_error( "source '" + token + "' is not a nuclide mixture this evaluation covers" );
  else
    sources_text = token;

  // The shielded configurations are lead, and it fluoresces: the Pb K series (72.8/75.0/84.9/87.4
  // keV) appears in these spectra at up to z=39 and no requested nuclide can produce it.  Requested
  // as an ELEMENT, so only its characteristic x-rays are modelled.
  if( !sources_text.empty() && SpecUtils::icontains( problem.id, "_Sh" ) )
    sources_text += ",Pb";

  for( const string &name : split_source_list( sources_text ) )
  {
    const RelActCalcAuto::SrcVariant src = RelActCalcAuto::source_from_string( name );
    if( RelActCalcAuto::is_null( src ) )
      throw runtime_error( "Unknown source '" + name + "' for inject problem " + problem.id );
    problem.requested_source_names.push_back( name );
    problem.sources.push_back( src );
  }

  // Spectra: record 0 = Poisson source + background, record 1 = background of the same duration.
  SpecMeas meas;
  if( !meas.load_file( problem.path, SpecUtils::ParserType::Auto ) || (meas.num_measurements() < 1) )
    throw runtime_error( "Failed to load '" + problem.path + "'" );
  const shared_ptr<const SpecUtils::Measurement> fg = meas.measurement_at_index( 0 );
  if( !fg || !fg->num_gamma_channels() )
    throw runtime_error( "No foreground gamma spectrum in '" + problem.path + "'" );
  problem.foreground = fg;
  problem.foreground_sample = fg->sample_number();
  if( meas.num_measurements() > 1 )
  {
    const shared_ptr<const SpecUtils::Measurement> bg = meas.measurement_at_index( 1 );
    if( bg && bg->num_gamma_channels() && SpecUtils::icontains( bg->title(), "Background" ) )
      problem.background = bg;
    else
      problem.notes.push_back( "record 1 of the PCF is not titled Background; no background used" );
  }

  // Truth: the resolution-merged "Expected Signal Photopeaks" rows as Gaussians on flat continua.
  auto truth = make_shared<InjectTruth>( parse_inject_truth_csv( truth_csv_path ) );
  problem.inject_truth = truth;
  const double fwhm_to_sigma = 1.0 / PhysicalUnits::fwhm_nsigma;
  for( const TruthPhotopeak &t : truth->merged_signal )
  {
    shared_ptr<PeakDef> peak = make_shared<PeakDef>( t.energy, t.fwhm * fwhm_to_sigma, t.area );
    shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
    cont->setType( PeakContinuum::OffsetType::Constant );
    cont->setRange( t.roi_lower, t.roi_upper );
    cont->setParameters( t.energy, vector<double>{ std::max( 0.0, t.continuum_area ) / (t.roi_upper - t.roi_lower) }, vector<double>{} );
    peak->setContinuum( cont );
    // Poisson-aware area uncertainty: the peak plus the continuum under +-1 FWHM.
    peak->setAmplitudeUncert( std::sqrt( std::max( 1.0, t.area + t.continuum_within_fwhm() ) ) );
    problem.truth_peaks.push_back( peak );
  }

  // Some inject truth files carry only the section headers (the expected-photopeak computation was
  // never run for that source); scoring against an empty truth would count every fitted peak as an
  // extra, so such problems are skipped rather than loaded.
  // An empty SIGNAL section is what makes a problem unscoreable (every fitted peak becomes an
  // extra with no truth to excuse it); a populated background section does not redeem it.
  if( truth->merged_signal.empty() )
    throw runtime_error( "no expected signal photopeaks in '" + truth_csv_path + "'" );
  return problem;
}//load_inject_problem(...)


CorpusProblem load_problem( const string &path, const map<string,string> &manifest_sources )
{
  CorpusProblem problem;
  problem.path = path;
  problem.id = problem_id_from_path( path );

  SpecMeas meas;
  if( !meas.load_N42_file( path ) )
    throw runtime_error( "Failed to load '" + path + "'" );

  for( const string &warning : meas.parse_warnings() )
    problem.notes.push_back( "parse warning: " + warning );

  int fg_sample = -1, bg_sample = -1;
  for( const int sample : meas.sample_numbers() )
  {
    const vector<shared_ptr<const SpecUtils::Measurement>> ms = meas.sample_measurements( sample );
    if( ms.empty() )
      continue;
    const SpecUtils::SourceType type = ms.front()->source_type();
    if( (type == SpecUtils::SourceType::Foreground) && (fg_sample < 0) )
      fg_sample = sample;
    else if( (type == SpecUtils::SourceType::Background) && (bg_sample < 0) )
      bg_sample = sample;
  }

  if( fg_sample < 0 )
  {
    // No Foreground-typed sample: a background-only reference file, or an unspecified type.
    for( const int sample : meas.sample_numbers() )
    {
      if( sample != bg_sample )
      {
        fg_sample = sample;
        break;
      }
    }
    if( fg_sample < 0 )
    {
      fg_sample = bg_sample;
      bg_sample = -1;
    }
  }

  if( fg_sample < 0 )
    throw runtime_error( "No spectra in '" + path + "'" );

  problem.foreground_sample = fg_sample;
  problem.foreground = meas.sum_measurements( set<int>{ fg_sample }, meas.detector_names(), nullptr );
  if( bg_sample >= 0 )
    problem.background = meas.sum_measurements( set<int>{ bg_sample }, meas.detector_names(), nullptr );

  if( !problem.foreground || !problem.foreground->num_gamma_channels() )
    throw runtime_error( "No foreground gamma spectrum in '" + path + "'" );

  // Reference peaks: the peak set covering the foreground sample, deep-copied with private continua
  // (fitting code mutates continua in place; see BatchPeak::get_exemplar_spectrum_and_peaks).
  shared_ptr<const std::deque<shared_ptr<const PeakDef>>> ref_peaks;
  for( const set<int> &samples : meas.sampleNumsWithPeaks() )
  {
    if( samples.count( fg_sample ) )
    {
      ref_peaks = meas.peaks( samples );
      break;
    }
  }

  if( ref_peaks )
  {
    map<shared_ptr<const PeakContinuum>, shared_ptr<PeakContinuum>> continua;
    for( const shared_ptr<const PeakDef> &p : *ref_peaks )
    {
      if( !p )
        continue;
      shared_ptr<PeakDef> copy = make_shared<PeakDef>( *p );
      const shared_ptr<const PeakContinuum> cont = p->continuum();
      if( cont )
      {
        auto pos = continua.find( cont );
        if( pos == end(continua) )
          pos = continua.emplace( cont, make_shared<PeakContinuum>( *cont ) ).first;
        copy->setContinuum( pos->second );
      }
      problem.truth_peaks.push_back( copy );
    }
  }//if( ref_peaks )

  // Requested sources: remark, manifest, then title
  string requested;
  const vector<shared_ptr<const SpecUtils::Measurement>> fg_meas = meas.sample_measurements( fg_sample );
  for( const shared_ptr<const SpecUtils::Measurement> &m : fg_meas )
  {
    for( const string &remark : m->remarks() )
    {
      if( SpecUtils::istarts_with( remark, "Requested sources:" ) )
      {
        requested = remark.substr( strlen("Requested sources:") );
        problem.notes.push_back( "sources from remark" );
        break;
      }
    }
    if( !requested.empty() )
      break;
  }

  if( requested.empty() )
  {
    const auto pos = manifest_sources.find( problem.id );
    if( pos != end(manifest_sources) )
    {
      requested = pos->second;
      problem.notes.push_back( "sources from manifest" );
    }
  }

  if( requested.empty() && !fg_meas.empty() && !fg_meas.front()->title().empty()
      && (fg_meas.front()->title().find( ',' ) != string::npos) )
  {
    requested = fg_meas.front()->title();
    problem.notes.push_back( "sources from title" );
  }

  problem.requested_source_names = split_source_list( requested );

  if( problem.requested_source_names.empty() && (problem.id.find( "Am241Li" ) != string::npos) )
  {
    problem.requested_source_names.push_back( "Am241" );
    problem.notes.push_back( "Am241Li problem: fitting Am241" );
  }

  for( const string &name : problem.requested_source_names )
  {
    const RelActCalcAuto::SrcVariant src = RelActCalcAuto::source_from_string( name );
    if( RelActCalcAuto::is_null( src ) )
    {
      problem.notes.push_back( "unrecognized source '" + name + "'" );
      continue;
    }
    if( std::find( begin(problem.sources), end(problem.sources), src ) == end(problem.sources) )
      problem.sources.push_back( src );
  }

  return problem;
}//load_problem


string fitted_fingerprint( const PeakSet &fitted )
{
  string answer;
  char buffer[128];
  for( const ScoredPeak &p : fitted.peaks )
  {
    const ScoredRoi &roi = fitted.rois[p.roi_index];
    snprintf( buffer, sizeof(buffer), "%.4f|%.6g|%.4f|%.4f|%d;",
              p.energy, p.amplitude, roi.lower, roi.upper, static_cast<int>(p.continuum_type) );
    answer += buffer;
  }
  return answer;
}

}//namespace FitPeaksCorpus
