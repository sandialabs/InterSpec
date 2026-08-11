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

#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <algorithm>
#include <optional>
#include <functional>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WString>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>

#include "Eigen/Dense"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/ExportSpecFileCAM.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;

namespace
{
  /** Merges lines the detector could not resolve from each other into a single yield-weighted
   line: the merged energy is the yield-weighted mean, the merged yield is the sum, and the
   merged energy uncertainty is a yield-weighted sample standard deviation of energy.

   Two lines are treated as unresolvable when they are within
   `ExportSpecFileCAM::sm_cluster_num_sigma` peak sigma of each other - the same criterion (and
   the same 1.25 sigma) the Activity/Shielding fit and manual rel-eff calculations use to decide
   which gamma lines contribute to a single peak (see
   `ShieldingSourceFitCalc::ShieldingSourceFitOptions::photopeak_cluster_sigma`).  Sigma comes
   from the `FWHM = fwhm_coeffs.first + fwhm_coeffs.second*sqrt(energy)` shape calibration being
   written into the CNF file, so the library ends up matching what Genie itself can resolve.

   Clusters grow greedily from low energy up, and each candidate is tested against the cluster's
   running yield-weighted mean, so a merged line is re-checked against the line that follows it.
   */
  void merge_unresolvable_lines( vector<ExportSpecFileCAM::GenieLibraryLine> &lines,
                                 const pair<float,float> &fwhm_coeffs )
  {
    if( lines.size() < 2 )
      return;

    std::sort( begin(lines), end(lines), []( const auto &lhs, const auto &rhs ){
      return lhs.energy < rhs.energy;
    } );

    const vector<float> coefs{ fwhm_coeffs.first, fwhm_coeffs.second };

    // Genie's FWHM equation can go negative at low energy (its own NaI defaults do); clamp to a
    //  small positive width so we never cluster on a nonsensical (or zero) resolution.
    const auto sigma_at = [&coefs]( const double energy ) -> double {
      const double sigma = DetectorPeakResponse::peakResolutionSigma( static_cast<float>(energy),
                                DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy,
                                coefs );
      return (sigma > 0.001) ? sigma : 0.001;
    };//sigma_at lambda

    vector<ExportSpecFileCAM::GenieLibraryLine> answer;

    size_t index = 0;
    while( index < lines.size() )
    {
      // Grow a cluster from lines[index], tracking the yield-weighted mean as we go, so each
      //  candidate is compared against where the merged line would actually sit.
      double yield_sum = lines[index].yield;
      double energy_num = lines[index].energy * static_cast<double>(lines[index].yield);
      double mean_energy = lines[index].energy;
      bool any_xray = lines[index].is_xray;

      size_t group_end = index + 1;
      while( group_end < lines.size() )
      {
        const double candidate_energy = lines[group_end].energy;
        if( (candidate_energy - mean_energy) > (ExportSpecFileCAM::sm_cluster_num_sigma * sigma_at(mean_energy)) )
          break;

        yield_sum += lines[group_end].yield;
        energy_num += candidate_energy * static_cast<double>(lines[group_end].yield);
        mean_energy = (yield_sum > 0.0) ? (energy_num / yield_sum) : candidate_energy;
        any_xray = any_xray || lines[group_end].is_xray;

        ++group_end;
      }//while( grow the cluster )

      if( (group_end - index) < 2 )
      {
        answer.push_back( lines[index] );
        index = group_end;
        continue;
      }

      // Weighted sample standard deviation of the merged energies.  For reliability weights the
      //  unbiased denominator is `V1 - V2/V1` (V1 = sum of weights, V2 = sum of squared weights);
      //  the `(n-1)/n * V1` form used before over-states the variance by 2x for a two-line merge.
      double energy_var_num = 0.0, weight_sq_sum = 0.0;
      for( size_t i = index; i < group_end; ++i )
      {
        energy_var_num += lines[i].yield * std::pow( lines[i].energy - mean_energy, 2.0 );
        weight_sq_sum += std::pow( static_cast<double>(lines[i].yield), 2.0 );
      }

      const double denom = yield_sum - (weight_sq_sum / yield_sum);
      const double energy_unc = ((yield_sum > 0.0) && (denom > 0.0))
                                 ? std::sqrt( energy_var_num / denom ) : 0.0;

      ExportSpecFileCAM::GenieLibraryLine merged;
      merged.energy = static_cast<float>( mean_energy );
      merged.energy_uncert = (energy_unc > 0.0) ? static_cast<float>(energy_unc) : -1.0f; //let Genie estimate it, if we couldn't compute one
      merged.yield = static_cast<float>( yield_sum );
      merged.yield_uncert = -1.0f; //let Genie estimate it
      merged.is_xray = any_xray;
      merged.included = true;

      answer.push_back( merged );
      index = group_end;
    }//while( index < lines.size() )

    lines.swap( answer );
  }//void merge_unresolvable_lines(...)


  /** Marks each source's highest-scoring line as the key line, approximating (for display
   purposes) what `CAMInputOutput::CAMIO::AssignKeyLines()` will independently compute when
   the file is written: score = energy_keV/1000 + yield_percent/10, and a candidate line is
   rejected - in favor of the source's next-best line - if any other line (of any source)
   falls within `key_line_intf_limit_kev` of it.  (`AssignKeyLines()` only actually checks the
   immediately-adjacent lines in a single, globally energy-sorted list, so this is not a
   byte-for-byte match, but should agree in the common case.)
   */
  void assign_key_lines( vector<ExportSpecFileCAM::GenieLibrarySource> &sources,
                         const float key_line_intf_limit_kev = 2.0f )
  {
    struct LineRef
    {
      size_t source_index;
      size_t line_index;
      float energy;
    };

    vector<LineRef> all_lines;
    for( size_t s = 0; s < sources.size(); ++s )
    {
      for( size_t l = 0; l < sources[s].lines.size(); ++l )
      {
        sources[s].lines[l].is_key_line = false;
        all_lines.push_back( LineRef{ s, l, sources[s].lines[l].energy } );
      }
    }

    const auto score = []( const ExportSpecFileCAM::GenieLibraryLine &ln ) -> float {
      return ln.energy/1000.0f + (100.0f*ln.yield)/10.0f;
    };

    for( size_t s = 0; s < sources.size(); ++s )
    {
      ExportSpecFileCAM::GenieLibrarySource &source = sources[s];
      if( source.lines.empty() )
        continue;

      if( source.lines.size() == 1 )
      {
        source.lines[0].is_key_line = true;
        continue;
      }

      vector<size_t> order( source.lines.size() );
      std::iota( begin(order), end(order), size_t(0) );
      std::sort( begin(order), end(order), [&source,&score]( const size_t lhs, const size_t rhs ){
        return score(source.lines[lhs]) > score(source.lines[rhs]);
      } );

      bool assigned = false;
      for( const size_t candidate_index : order )
      {
        const float energy = source.lines[candidate_index].energy;

        bool interferes = false;
        for( const LineRef &ref : all_lines )
        {
          if( (ref.source_index == s) && (ref.line_index == candidate_index) )
            continue;
          if( std::fabs(ref.energy - energy) <= key_line_intf_limit_kev )
          {
            interferes = true;
            break;
          }
        }//for( check every other line for interference )

        if( !interferes )
        {
          source.lines[candidate_index].is_key_line = true;
          assigned = true;
          break;
        }
      }//for( candidate_index : order )

      if( !assigned ) //every candidate had an interference; fall back to the highest-score line
        source.lines[order.front()].is_key_line = true;
    }//for( each source )
  }//void assign_key_lines(...)


  /** Returns the yield (photons, or x-rays, per decay) of the line in `photons` closest to
   `energy_kev`, if within `tolerance_kev`; otherwise returns a negative value.
   */
  double find_yield_near( const vector<SandiaDecay::EnergyRatePair> &photons,
                          const float energy_kev, const float tolerance_kev = 1.0f )
  {
    double best_yield = -1.0;
    double best_de = tolerance_kev;
    for( const SandiaDecay::EnergyRatePair &p : photons )
    {
      const double de = std::fabs( p.energy - energy_kev );
      if( de <= best_de )
      {
        best_de = de;
        best_yield = p.numPerSecond; // numPerSecond, for an activity of 1 becquerel, is photons/decay
      }
    }
    return best_yield;
  }//find_yield_near(...)
}//namespace


namespace ExportSpecFileCAM
{

vector<GenieLibrarySource> build_genie_library(
                                const deque<shared_ptr<const PeakDef>> &peaks,
                                const GenieLibraryLineMode mode,
                                const double yield_threshold_percent,
                                const bool combine_unresolvable_lines,
                                const pair<float,float> &fwhm_coeffs,
                                const pair<float,float> &energy_range,
                                const map<const SandiaDecay::Nuclide *, double> &nuclide_ages,
                                vector<string> *warnings )
{
  // Group the peaks that have a usable source (a nuclide) by that nuclide; peaks with a
  // reaction, or an x-ray with no parent nuclide, cant be represented in a Genie library.
  map<const SandiaDecay::Nuclide *, vector<shared_ptr<const PeakDef>>> peaks_by_nuclide;

  for( const shared_ptr<const PeakDef> &peak : peaks )
  {
    if( !peak || !peak->hasSourceGammaAssigned() )
      continue;

    // Escape peaks sit at E-511/E-1022, which is not where the nuclide emits anything - writing
    //  that energy into a library would attribute a real yield to an energy the source has no
    //  line at, so skip them rather than let the energy match find some unrelated neighbour.
    const PeakDef::SourceGammaType src_type = peak->sourceGammaType();
    if( (src_type == PeakDef::SourceGammaType::SingleEscapeGamma)
       || (src_type == PeakDef::SourceGammaType::DoubleEscapeGamma) )
    {
      if( warnings )
        warnings->push_back( WString::tr("esfcam-warn-escape-peak")
                             .arg( SpecUtils::printCompact(peak->mean(), 5) ).toUTF8() );
      continue;
    }

    const SandiaDecay::Nuclide *nuc = peak->parentNuclide();
    if( nuc )
    {
      peaks_by_nuclide[nuc].push_back( peak );
      continue;
    }

    if( peak->reaction() )
    {
      if( warnings )
        warnings->push_back( WString::tr("esfcam-warn-reaction-peak")
                             .arg( SpecUtils::printCompact(peak->mean(), 5) )
                             .arg( peak->reaction()->name() ).toUTF8() );
      continue;
    }

    // A peak assigned to a bare element (fluorescence x-ray) has no nuclide and hence no
    //  half-life, so it cannot be a Genie library entry.  Decay x-rays - which DO have a parent
    //  nuclide - are handled above, along with that nuclide's gammas.
    if( peak->xrayElement() )
    {
      if( warnings )
        warnings->push_back( WString::tr("esfcam-warn-fluorescence-xray-peak")
                             .arg( SpecUtils::printCompact(peak->mean(), 5) )
                             .arg( peak->xrayElement()->symbol ).toUTF8() );
      continue;
    }
  }//for( loop over peaks )

  vector<GenieLibrarySource> answer;

  for( const auto &nuc_peaks : peaks_by_nuclide )
  {
    const SandiaDecay::Nuclide * const nuc = nuc_peaks.first;

    GenieLibrarySource source;
    source.nuclide = nuc;
    source.name = nuc->symbol;
    source.half_life_seconds = static_cast<float>( nuc->halfLife / PhysicalUnits::second );
    source.half_life_uncert_seconds = -1.0f; //not tracked by SandiaDecay; let Genie estimate it
    source.included = true;
    source.is_ageable = !PeakDef::ageFitNotAllowed( nuc );

    // Even when the user cant usefully vary this nuclide's age (e.g. Cs137, whose spectrum
    // doesnt meaningfully change once its short-lived Ba137m daughter reaches equilibrium),
    // we still need a physically-reasonable default age - not 0 - so that daughter-product
    // line yields (like Cs137's 661.7 keV line) are computed at equilibrium, not at t=0 when
    // the daughter population (and its gamma yield) is still zero.
    double age = PeakDef::defaultDecayTime( nuc );
    if( source.is_ageable )
    {
      const auto age_iter = nuclide_ages.find( nuc );
      if( age_iter != end(nuclide_ages) )
        age = age_iter->second;
    }
    source.age_seconds = age / PhysicalUnits::second;

    SandiaDecay::NuclideMixture mixture;
    mixture.addNuclideByActivity( nuc, 1.0*PhysicalUnits::becquerel );
    const vector<SandiaDecay::EnergyRatePair> photons
                    = mixture.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );

    // `addNuclideByActivity(...)` sets the parent to 1 Bq at t=0, so `numPerSecond` at `age` is
    //  photons per second per *initial* becquerel - it still carries the parent's decay factor.
    //  A Genie library wants photons per decay of the parent, so divide the parent's activity at
    //  `age` back out.  (This is the same `age_sf` correction
    //  `GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities(...)` makes; without
    //  it, e.g. Cs137 at its default age comes out 5.66x too low, and any nuclide whose default
    //  age is 7 half-lives comes out 128x too low.)
    const double parent_activity = mixture.activity( age, nuc );
    if( parent_activity <= 0.0 )
    {
      if( warnings )
        warnings->push_back( WString::tr("esfcam-warn-decayed-away").arg( source.name ).toUTF8() );
      continue;
    }

    const double yield_sf = 1.0 / parent_activity;

    double max_yield = 0.0;
    for( const SandiaDecay::EnergyRatePair &p : photons )
      max_yield = std::max( max_yield, p.numPerSecond );

    if( mode == GenieLibraryLineMode::PeakLinesOnly )
    {
      for( const shared_ptr<const PeakDef> &peak : nuc_peaks.second )
      {
        // Note: use `gammaParticleEnergy()` for x-rays too.  `PeakDef::setNuclearTransition(...)`
        //  zeroes `m_xrayEnergy` whenever a transition is set, so `xrayEnergy()` is always 0 for a
        //  decay x-ray that has a parent nuclide - which is every peak that reaches here.
        const bool is_xray = (peak->sourceGammaType() == PeakDef::SourceGammaType::XrayGamma);
        const float energy = peak->gammaParticleEnergy();

        const double yield = find_yield_near( photons, energy );
        if( yield <= 0.0 )
        {
          if( warnings )
            warnings->push_back( WString::tr("esfcam-warn-no-yield")
                                 .arg( source.name )
                                 .arg( SpecUtils::printCompact(energy, 5) ).toUTF8() );
          continue;
        }

        // Two peaks can be assigned to the same gamma (e.g. in different sample sets, or a
        //  mis-assignment); writing the line twice would double its yield once the duplicates get
        //  merged as "unresolvable", making Genie report half the true activity.
        const bool already_have = std::any_of( begin(source.lines), end(source.lines),
                                    [energy]( const GenieLibraryLine &l ){
          return (std::fabs(l.energy - energy) < 0.001f);
        } );

        if( already_have )
          continue;

        GenieLibraryLine line;
        line.energy = energy;
        line.yield = static_cast<float>( yield_sf * yield );
        line.is_xray = is_xray;
        line.included = true;
        source.lines.push_back( line );
      }//for( loop over this nuclide's peaks )
    }else
    {
      assert( mode == GenieLibraryLineMode::AllLinesAboveThreshold );

      const double threshold = (yield_threshold_percent/100.0) * max_yield;

      // Figure out if a line is an x-ray by checking if it's also in xrays(...); gammas()
      // and xrays() are disjoint, and photons() is their union (plus annihilation gammas).
      const vector<SandiaDecay::EnergyRatePair> xrays = mixture.xrays( age );

      for( const SandiaDecay::EnergyRatePair &p : photons )
      {
        if( p.numPerSecond < threshold )
          continue;

        // Lines outside the measurement's energy range cannot correspond to a peak Genie could
        //  ever find, and just add noise to the library.
        if( (energy_range.second > energy_range.first)
           && ((p.energy < energy_range.first) || (p.energy > energy_range.second)) )
          continue;

        const bool is_xray = std::any_of( begin(xrays), end(xrays),
                                          [&p]( const SandiaDecay::EnergyRatePair &x ){
          return std::fabs(x.energy - p.energy) < 0.001;
        } );

        GenieLibraryLine line;
        line.energy = static_cast<float>( p.energy );
        line.yield = static_cast<float>( yield_sf * p.numPerSecond );
        line.is_xray = is_xray;
        line.included = true;
        source.lines.push_back( line );
      }//for( loop over all photons of this nuclide )
    }//if( PeakLinesOnly ) / else

    if( source.lines.empty() )
      continue;

    if( combine_unresolvable_lines )
      merge_unresolvable_lines( source.lines, fwhm_coeffs );

    answer.push_back( source );
  }//for( loop over peaks_by_nuclide )

  std::sort( begin(answer), end(answer), []( const GenieLibrarySource &lhs, const GenieLibrarySource &rhs ){
    return lhs.name < rhs.name;
  } );

  assign_key_lines( answer );

  return answer;
}//build_genie_library(...)


vector<CAMInputOutput::CnfGenieExtras::LibraryLine>
    to_library_lines( const vector<GenieLibrarySource> &sources )
{
  vector<CAMInputOutput::CnfGenieExtras::LibraryLine> answer;

  for( const GenieLibrarySource &source : sources )
  {
    if( !source.included )
      continue;

    for( const GenieLibraryLine &line : source.lines )
    {
      if( !line.included )
        continue;

      CAMInputOutput::CnfGenieExtras::LibraryLine ll;
      ll.nuclide_name = source.name;
      ll.half_life_seconds = source.half_life_seconds;
      ll.half_life_uncert_seconds = source.half_life_uncert_seconds;
      ll.energy = line.energy;
      ll.energy_uncert = line.energy_uncert;
      // GENIE library abundances are PERCENT, not a fraction - real Genie-produced libraries hold
      //  e.g. 85.6 for Ni57's 122 keV line.  `GenieLibraryLine::yield` is a fraction, so convert
      //  here, at the boundary.
      ll.yield = 100.0f * line.yield;
      ll.yield_uncert = (line.yield_uncert < 0.0f) ? -1.0f : (100.0f * line.yield_uncert);
      ll.no_weight_mean = line.is_xray;
      // Pass our key-line choice explicitly, so what the user previewed is what gets written
      //  (`CAMIO::AssignKeyLines()` would otherwise pick its own, and could disagree).
      ll.is_key_line = line.is_key_line;
      answer.push_back( ll );
    }//for( loop over source's lines )
  }//for( loop over sources )

  return answer;
}//to_library_lines(...)


namespace
{
  /** Least-squares fit of Genie's `FWHM = a + b*sqrt(energy)` to `fwhm_fcn` over
   [lower_energy, upper_energy].

   Samples are log-spaced and weighted by `1/FWHM`, so the fit minimizes *relative* error evenly
   across the range.  Sampling linearly and fitting absolute residuals instead puts almost every
   sample above a few hundred keV and lets the (much wider) high-energy end dominate, which
   noticeably under-predicts the width at low energy - where peaks are narrowest and a bad width
   costs the most.
   */
  pair<float,float> fit_genie_fwhm( const std::function<double(double)> &fwhm_fcn,
                                    const double lower_energy, const double upper_energy )
  {
    const int num_points = 16;
    Eigen::MatrixXd A( num_points, 2 );
    Eigen::VectorXd b( num_points );

    for( int i = 0; i < num_points; ++i )
    {
      const double frac = static_cast<double>(i) / (num_points - 1);
      const double energy = lower_energy * std::pow( upper_energy/lower_energy, frac );
      const double fwhm = fwhm_fcn( energy );
      const double weight = (fwhm > 0.0) ? (1.0 / fwhm) : 1.0;

      A(i,0) = weight;
      A(i,1) = weight * std::sqrt( energy );
      b(i) = weight * fwhm;
    }

    const Eigen::VectorXd sol = A.colPivHouseholderQr().solve( b );

    return make_pair( static_cast<float>(sol(0)), static_cast<float>(sol(1)) );
  }//fit_genie_fwhm(...)
}//namespace


double genie_low_tail( const PeakDef &peak )
{
  // Genie's low tail is a Gaussian core stitched to an exponential below `centroid - T`:
  //    P(x) = H*exp(-0.5*u^2)                       for u >= -T/sigma
  //    P(x) = H*exp(0.5*(T/sigma)^2 + (T/sigma)*u)  for u <  -T/sigma,   u = (x-mean)/sigma
  // which is exactly `PeakDists::gauss_exp_pdf(...)` with `skew == T/sigma`.  So a peak fit with
  // InterSpec's GaussExp skew converts exactly, and ExpGaussExp's lower-side parameter uses the
  // same convention.  Other skew forms (Bortel, Crystal Ball, ...) are different functions, and a
  // fitted-but-wrong tail would be worse for Genie than none, so those report no tail.
  //
  // T shares units with sigma (the tail exponent has to be dimensionless), and the CAM peak
  // record's FWHM is in keV - so T is in keV too.
  const double sigma = peak.sigma();
  if( sigma <= 0.0 )
    return -1.0;

  double skew = 0.0;   //T/sigma
  switch( peak.skewType() )
  {
    case PeakDef::SkewType::GaussExp:
    case PeakDef::SkewType::ExpGaussExp:
      skew = peak.coefficient( PeakDef::CoefficientType::SkewPar0 );
      break;

    default:
      break;
  }//switch( peak.skewType() )

  return (skew > 0.0) ? (skew * sigma) : -1.0;
}//genie_low_tail(...)


std::optional<pair<float,float>> fit_genie_low_tail_cal(
                                const deque<shared_ptr<const PeakDef>> &peaks )
{
  // Genie carries a low-tail curve `T(E) = B2 + B3*E` alongside the FWHM shape calibration, and
  //  uses it when fitting peaks itself.  Derive it from whichever peaks actually have a tail, so
  //  the curve agrees with the per-peak T values written into the PEAK block.
  vector<pair<double,double>> energy_tail;  //{energy, T}
  for( const shared_ptr<const PeakDef> &peak : peaks )
  {
    if( !peak || !peak->gausPeak() )
      continue;

    const double tail = genie_low_tail( *peak );
    if( tail > 0.0 )
      energy_tail.emplace_back( peak->mean(), tail );
  }

  if( energy_tail.empty() )
    return {};   //no tails - leave the calibration at zero, which is what "no tail" looks like

  if( energy_tail.size() == 1 )
    return make_pair( static_cast<float>(energy_tail[0].second), 0.0f );

  // Unweighted linear least squares of T against energy.
  double sum_e = 0.0, sum_t = 0.0, sum_ee = 0.0, sum_et = 0.0;
  for( const pair<double,double> &et : energy_tail )
  {
    sum_e += et.first;
    sum_t += et.second;
    sum_ee += et.first * et.first;
    sum_et += et.first * et.second;
  }

  const double n = static_cast<double>( energy_tail.size() );
  const double denom = (n*sum_ee) - (sum_e*sum_e);
  if( !(std::fabs(denom) > 0.0) )
    return make_pair( static_cast<float>(sum_t/n), 0.0f );

  const double slope = ((n*sum_et) - (sum_e*sum_t)) / denom;
  const double offset = (sum_t - (slope*sum_e)) / n;

  if( !std::isfinite(offset) || !std::isfinite(slope) )
    return {};

  return make_pair( static_cast<float>(offset), static_cast<float>(slope) );
}//fit_genie_low_tail_cal(...)


vector<CAMInputOutput::Peak> to_cam_peaks(
                        const deque<shared_ptr<const PeakDef>> &peaks,
                        const shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                        const float live_time_s,
                        const GenieEfficiencyResult * const efficiency,
                        const shared_ptr<const SpecUtils::Measurement> &data )
{
  vector<CAMInputOutput::Peak> answer;

  const bool have_cal = !!energy_cal && energy_cal->valid();
  const bool have_lt = (live_time_s > 0.0f);

  for( const shared_ptr<const PeakDef> &peak : peaks )
  {
    // A "data defined" peak has no fitted centroid or width to report.
    if( !peak || !peak->gausPeak() )
      continue;

    CAMInputOutput::Peak p{};

    p.Energy = static_cast<float>( peak->mean() );
    p.FullWidthAtHalfMaximum = static_cast<float>( peak->fwhm() );
    p.Area = static_cast<float>( peak->peakArea() );
    p.AreaUncertainty = static_cast<float>( peak->peakAreaUncert() );
    // Genie's low tail is a Gaussian core stitched to an exponential below `centroid - T`:
    //    P(x) = H*exp(-0.5*u^2)              for u >= -T/sigma
    //    P(x) = H*exp(0.5*(T/sigma)^2 + (T/sigma)*u)  for u <  -T/sigma,   u = (x-mean)/sigma
    // which is exactly `PeakDists::gauss_exp_pdf(...)` with `skew == T/sigma`.  So a peak fit with
    // InterSpec's GaussExp skew converts exactly; ExpGaussExp's lower-side parameter uses the same
    // convention.  Other skew forms (Bortel, Crystal Ball, ...) are different functions, and a
    // fitted-but-wrong tail would be worse for Genie than none, so those get the "no tail" value.
    p.LowTail = CAMInputOutput::Peak::sm_no_low_tail;

    const double low_tail = genie_low_tail( *peak );
    if( low_tail > 0.0 )
      p.LowTail = static_cast<float>( low_tail );  //T, in keV (see CAMInputOutput::Peak::LowTail)

    if( have_cal )
    {
      const double centroid_channel = energy_cal->channel_for_energy( peak->mean() );
      p.Centroid = static_cast<float>( centroid_channel );

      // Convert the mean's uncertainty to channels using the local gain.
      const double sigma_energy = peak->meanUncert();
      if( sigma_energy > 0.0 )
      {
        const double lower = energy_cal->channel_for_energy( peak->mean() - sigma_energy );
        const double upper = energy_cal->channel_for_energy( peak->mean() + sigma_energy );
        p.CentroidUncertainty = static_cast<float>( 0.5*std::fabs(upper - lower) );
      }

      const double left_channel = energy_cal->channel_for_energy( peak->lowerX() );
      const double right_channel = energy_cal->channel_for_energy( peak->upperX() );
      p.LeftChannel = static_cast<int>( std::floor(std::max(0.0, left_channel)) );
      p.RightChannel = static_cast<int>( std::ceil(std::max(0.0, right_channel)) );
      if( p.RightChannel < p.LeftChannel )
        p.RightChannel = p.LeftChannel;
    }//if( have_cal )

    // Counts in the continuum under the peak's ROI.
    //
    // A stepped continuum's shape comes from the data, and a CDF-step's also from the other peaks
    //  sharing the ROI, so both have to be supplied: `offset_integral` throws without the spectrum,
    //  which is why every stepped ROI used to be written with a continuum area of zero.
    const shared_ptr<const PeakContinuum> continuum = peak->continuum();
    if( continuum )
    {
      vector<shared_ptr<const PeakDef>> roi_peaks;
      for( const shared_ptr<const PeakDef> &other : peaks )
      {
        if( other && (other->continuum() == continuum) )
          roi_peaks.push_back( other );
      }

      try
      {
        p.Continuum = static_cast<float>( continuum->offset_integral( peak->lowerX(), peak->upperX(),
                                                                      data, roi_peaks ) );
      }catch( std::exception &e )
      {
        // Only reachable for a stepped continuum with no usable spectrum; a zero continuum area is
        //  at least honest, and better than refusing the export.
        cerr << "to_cam_peaks: continuum area for the peak at " << peak->mean()
             << " keV could not be computed: " << e.what() << endl;
        p.Continuum = 0.0f;
      }
    }//if( continuum )

    if( have_lt )
      p.CountRate = static_cast<float>( peak->peakArea() / live_time_s );

    // Note: this field is a RELATIVE uncertainty, as a percent - not counts per second (verified
    //  against a real Genie file; see `CAMInputOutput::Peak`).
    if( peak->peakArea() > 0.0 )
      p.CountRateUncertainty = static_cast<float>( 100.0 * peak->peakAreaUncert() / peak->peakArea() );

    // What Genie's report prints as "Peak significance"; see `CAMInputOutput::Peak` for why this
    //  is Area/AreaUncertainty rather than the peak-search statistic Genie itself puts here.
    if( peak->peakAreaUncert() > 0.0 )
      p.Significance = static_cast<float>( peak->peakArea() / peak->peakAreaUncert() );

    // Poisson sigma of the continuum counts under the peak, and the Currie decision limit built
    //  on it.  Genie's own files have `CriticalLevel/BackgroundSigma` around 2.13; 2.33 is the
    //  textbook coefficient, and near enough that the pair stays self-consistent.
    if( p.Continuum > 0.0f )
    {
      p.BackgroundSigma = std::sqrt( p.Continuum );
      p.CriticalLevel = 2.33f * p.BackgroundSigma;
    }

    // The efficiency curve's value at the peak, so the PEAK and GEOM blocks agree - this is what
    //  Genie fills in, and it is left zero for a file with no efficiency calibration.
    if( efficiency )
    {
      const double eff = genie_efficiency_at( *efficiency, peak->mean() );
      if( eff > 0.0 )
      {
        p.Efficiency = static_cast<float>( eff );
        p.EfficiencyUncertainty
              = static_cast<float>( eff * genie_default_eff_uncertainty( static_cast<float>(peak->mean()) ) );
      }
    }//if( efficiency )

    answer.push_back( p );
  }//for( loop over peaks )

  // Genie peak lists are energy ordered.
  std::sort( begin(answer), end(answer),
            []( const CAMInputOutput::Peak &lhs, const CAMInputOutput::Peak &rhs ){
    return lhs.Energy < rhs.Energy;
  } );

  return answer;
}//to_cam_peaks(...)


pair<float,float> default_genie_fwhm( const bool is_hpge )
{
  // For HPGe, use the same Ge default `CAMIO::AddDetectorType(...)` writes, so an export with
  //  FWHM writing turned off and one with it left at the default agree.
  if( is_hpge )
    return make_pair( 1.0f, 0.035f );

  // The Genie manual's own NaI defaults ({-7.0, 2.0}) go negative below ~12 keV and noticeably
  //  under-predict the width at low energy, so instead fit the Genie equation form to
  //  InterSpec's nominal "NaI 3x3" resolution over the range NaI spectra are actually used over.
  //  Fit from 20 keV (not 60) so the low-energy end is interpolated rather than extrapolated -
  //  fitting only 60-2614 keV drove the offset far enough negative that the equation went
  //  negative below ~18 keV and was ~50% low at 30 keV, which is worse than Genie's own default.
  static const pair<float,float> s_nai_coeffs
        = fit_genie_fwhm( []( const double energy ) -> double {
            return PeakFitUtils::nai_fwhm_fcn( static_cast<float>(energy) );
          }, 20.0, 2614.0 );

  return s_nai_coeffs;
}//default_genie_fwhm(...)


pair<float,float> fit_genie_fwhm_from_drf( const DetectorPeakResponse &drf )
{
  double lower_energy = drf.lowerEnergy();
  double upper_energy = drf.upperEnergy();
  if( !(upper_energy > lower_energy) )
  {
    lower_energy = 59.0;
    upper_energy = 2614.0;
  }

  // A DRF can legitimately report a lower energy of 0 (it is the default, and can also come from
  //  a tabulated efficiency whose first point is at 0 keV).  `fit_genie_fwhm` samples
  //  logarithmically, so a zero lower bound would make every sample energy 0*inf = NaN, and the
  //  NaN coefficients would go straight into the file.
  lower_energy = std::max( lower_energy, 10.0 );
  upper_energy = std::max( upper_energy, lower_energy + 1.0 );

  const pair<float,float> answer = fit_genie_fwhm( [&drf]( const double energy ) -> double {
    return drf.peakResolutionFWHM( static_cast<float>(energy) );
  }, lower_energy, upper_energy );

  if( !std::isfinite(answer.first) || !std::isfinite(answer.second) )
    throw runtime_error( "fit_genie_fwhm_from_drf: the detector response function's FWHM could not"
                         " be fit to Genie's FWHM = A0 + A1*sqrt(energy) form." );

  return answer;
}//fit_genie_fwhm_from_drf(...)


namespace
{
  /** Linear least-squares fit of `ln(y) = sum_i{ coeffs[i]*ln(x)^i }`.

   `rel_uncerts` gives each point's relative uncertainty on `y`, which is also the absolute
   uncertainty of `ln(y)`, so the fit is weighted by `1/rel_uncerts[i]^2` - the same weighting
   Genie uses (see `CAMIO::AddEfficiencyFit(...)`).  Pass an empty vector for an unweighted fit.
   */
  vector<double> fit_log_log_polynomial( const vector<double> &x, const vector<double> &y,
                                         const vector<double> &rel_uncerts, const size_t order )
  {
    assert( rel_uncerts.empty() || (rel_uncerts.size() == x.size()) );

    const size_t n = x.size();
    const size_t num_coeffs = order + 1;
    Eigen::MatrixXd A( n, num_coeffs );
    Eigen::VectorXd b( n );

    for( size_t row = 0; row < n; ++row )
    {
      const double sigma = rel_uncerts.empty() ? 1.0 : rel_uncerts[row];
      const double weight = (sigma > 0.0) ? (1.0 / sigma) : 1.0;

      const double lnx = std::log( x[row] );
      double pow_term = 1.0;
      for( size_t col = 0; col < num_coeffs; ++col )
      {
        A(row,col) = weight * pow_term;
        pow_term *= lnx;
      }
      b(row) = weight * std::log( y[row] );
    }

    const Eigen::VectorXd sol = A.colPivHouseholderQr().solve( b );

    vector<double> coeffs( num_coeffs );
    for( size_t i = 0; i < num_coeffs; ++i )
      coeffs[i] = sol(static_cast<int>(i));
    return coeffs;
  }//fit_log_log_polynomial(...)


  /** Evaluates `exp( sum_i{ coeffs[i]*ln(energy)^i } )`. */
  double eval_log_log_polynomial( const vector<double> &coeffs, const double energy )
  {
    double ln_eff = 0.0, pow_term = 1.0;
    const double lnE = std::log( energy );
    for( const double c : coeffs )
    {
      ln_eff += c * pow_term;
      pow_term *= lnE;
    }

    return std::exp( ln_eff );
  }//eval_log_log_polynomial(...)


  /** Genie's own default polynomial order, based on the number of calibration points, per the
   Genie 2000 Customization Tools Manual's "Efficiency Calibration" chapter.
   */
  size_t genie_default_poly_order( const size_t num_points )
  {
    if( num_points >= 10 )
      return 5;
    if( num_points >= 8 )
      return 4;
    if( num_points >= 6 )
      return 3;
    return 2;
  }
}//namespace


vector<double> fit_genie_efficiency_curve( const vector<double> &energies,
                                           const vector<double> &efficiencies,
                                           const vector<double> &rel_uncerts,
                                           const size_t order )
{
  if( energies.size() != efficiencies.size() )
    throw runtime_error( "fit_genie_efficiency_curve: energy and efficiency counts differ." );

  if( energies.size() < (order + 1) )
    throw runtime_error( "fit_genie_efficiency_curve: fewer points than coefficients." );

  return fit_log_log_polynomial( energies, efficiencies, rel_uncerts, order );
}//fit_genie_efficiency_curve(...)


double genie_efficiency_at( const GenieEfficiencyResult &efficiency, const double energy )
{
  if( !(energy > 0.0) )
    return 0.0;

  if( !efficiency.fit_coeffs.empty() )
  {
    const vector<double> coeffs( begin(efficiency.fit_coeffs), end(efficiency.fit_coeffs) );
    const double eff = eval_log_log_polynomial( coeffs, energy );
    return std::isfinite(eff) ? std::max( 0.0, eff ) : 0.0;
  }//if( we wrote a fitted curve )

  // No fit, so do what Genie's "Interpolated" model does with the points we did write.  They are
  //  energy ordered (see `convert_efficiency_to_genie(...)`), and outside their range there is
  //  nothing to interpolate - Genie leaves the peak's efficiency zero there rather than
  //  extrapolating, and so do we.
  const vector<CAMInputOutput::EfficiencyPoint> &pts = efficiency.points;
  if( (pts.size() < 2) || (energy < pts.front().Energy) || (energy > pts.back().Energy) )
    return 0.0;

  for( size_t i = 1; i < pts.size(); ++i )
  {
    if( energy > pts[i].Energy )
      continue;

    const double e0 = pts[i-1].Energy, e1 = pts[i].Energy;
    const double y0 = pts[i-1].Efficiency, y1 = pts[i].Efficiency;
    if( !(e0 > 0.0) || !(y0 > 0.0) || !(y1 > 0.0) || (e1 <= e0) )
      return y0;

    const double frac = (std::log(energy) - std::log(e0)) / (std::log(e1) - std::log(e0));
    return std::exp( std::log(y0) + frac*(std::log(y1) - std::log(y0)) );
  }//for( find the bracketing points )

  return pts.back().Efficiency;
}//genie_efficiency_at(...)


float genie_default_eff_uncertainty( const float energy )
{
  if( energy <= 55.0f )
    return 0.15f;
  if( energy <= 175.0f )
    return 0.10f;
  if( energy <= 400.0f )
    return 0.08f;
  if( energy <= 900.0f )
    return 0.06f;
  return 0.04f;
}//genie_default_eff_uncertainty(...)


GenieEfficiencyResult convert_efficiency_to_genie( const DetectorPeakResponse &drf,
                                                   const double distance,
                                                   vector<string> *warnings )
{
  GenieEfficiencyResult answer;

  // Points below this are dropped rather than clamped: a clamped point is a fabricated efficiency,
  //  and Genie would happily interpolate it and report activities orders of magnitude too high.
  const double min_useful_eff = 1.0E-9;
  // Two points is the minimum an interpolated curve needs; a DRF's own tabulated efficiency can
  //  legitimately have only a handful.
  const size_t min_useful_points = 2;

  // Genie wants absolute (per gamma emitted at the source) efficiency; InterSpec's
  //  `intrinsicEfficiency(...)` is per gamma striking the detector face, so for a far-field DRF
  //  the solid angle has to be folded back in - this is the same
  //  `fixed_geom ? intrinsicEfficiency(E) : efficiency(E,dist)` split the detection-limit and
  //  activity-fit code uses.
  const bool fixed_geom = drf.isFixedGeometry();

  if( !fixed_geom && !(distance > 0.0) )
    throw runtime_error( "convert_efficiency_to_genie: a positive source distance is required for"
                         " a non-fixed-geometry detector response function." );

  const auto absolute_eff = [&drf,fixed_geom,distance]( const double energy ) -> double {
    return fixed_geom ? drf.intrinsicEfficiency( static_cast<float>(energy) )
                      : drf.efficiency( static_cast<float>(energy), distance );
  };//absolute_eff lambda

  const DetectorPeakResponse::EfficiencyFnctForm form = drf.efficiencyFcnType();
  const bool from_pairs = (form == DetectorPeakResponse::EfficiencyFnctForm::kEnergyEfficiencyPairs);

  // A tabulated DRF contributes its own point energies; the analytic forms
  //  (kExpOfLogPowerSeries, kFunctialEfficienyForm) are sampled log-spaced across their valid
  //  range, since efficiency curves are fit in log-log space.
  vector<double> energies, effs;
  size_t num_tried = 0, num_dropped = 0;

  if( from_pairs )
  {
    const vector<DetectorPeakResponse::EnergyEfficiencyPair> &pairs = drf.getEnergyEfficiencyPair();
    num_tried = pairs.size();

    for( const DetectorPeakResponse::EnergyEfficiencyPair &pt : pairs )
    {
      // Note: use the DRF's own tabulated energies, but go back through the DRF for the value, so
      //  the solid-angle (and any absolute-to-intrinsic) correction is applied consistently.
      const double eff = absolute_eff( pt.energy );
      if( !(eff > min_useful_eff) )
      {
        num_dropped += 1;
        continue;
      }

      energies.push_back( pt.energy );
      effs.push_back( eff );
    }//for( loop over the DRF's tabulated points )
  }else
  {
    double lower_energy = drf.lowerEnergy();
    double upper_energy = drf.upperEnergy();
    if( !(upper_energy > lower_energy) )
    {
      lower_energy = 59.0;
      upper_energy = 2614.0;
    }
    lower_energy = std::max( lower_energy, 10.0 ); //avoid ln(0)/negative energies

    // Enough points that Genie's own point-count rule (`genie_default_poly_order`) reaches its
    //  highest order, with margin for any that get dropped.
    const size_t num_samples = 20;
    num_tried = num_samples;

    for( size_t i = 0; i < num_samples; ++i )
    {
      const double frac = static_cast<double>(i) / (num_samples - 1);
      const double energy = lower_energy * std::pow( upper_energy/lower_energy, frac );
      const double eff = absolute_eff( energy );

      // Drop, dont clamp - see `min_useful_eff`.
      if( !(eff > min_useful_eff) )
      {
        num_dropped += 1;
        continue;
      }

      energies.push_back( energy );
      effs.push_back( eff );
    }//for( loop over sample energies )
  }//if( from_pairs ) / else

  if( num_dropped && warnings )
    warnings->push_back( WString::tr("esfcam-warn-eff-points-dropped")
                         .arg( static_cast<int>(num_dropped) ).toUTF8() );

  if( energies.size() < min_useful_points )
    throw runtime_error( "convert_efficiency_to_genie: too few usable efficiency points" );

  // For the sampled branch we also know how many points we *tried* - if most of the range is
  //  negligible the curve is not worth writing, whatever the absolute count.  (A tabulated DRF's
  //  points are wherever its author put them, so the same test would be meaningless there.)
  if( !from_pairs && (energies.size() < (num_tried / 2)) )
    throw runtime_error( "convert_efficiency_to_genie: the detector response function's absolute"
                         " efficiency is negligible over much of its energy range at this"
                         " distance." );

  // Genie weights its efficiency fit by 1/sigma_rel^2, so these uncertainties shape the curve, and
  //  a zero would leave it undefined - which is why Genie could previously only offer its
  //  "Interpolated" model for our files.  No DRF form carries per-point uncertainties today
  //  (`EnergyEfficiencyPair` has no uncertainty member at all), so the defaults are used
  //  throughout; see `genie_default_eff_uncertainty(...)`.
  vector<double> rel_uncerts( energies.size() );
  for( size_t i = 0; i < energies.size(); ++i )
    rel_uncerts[i] = genie_default_eff_uncertainty( static_cast<float>(energies[i]) );

  // Fit the same curve Genie would, so the file carries it whether or not Genie recomputes it.
  //  The order must leave at least as many points as coefficients, and CAM has room for
  //  `sm_geom_max_fit_coeffs` of them.
  size_t order = genie_default_poly_order( energies.size() );
  order = std::min( order, energies.size() - 1 );
  order = std::min( order, CAMInputOutput::CAMIO::sm_geom_max_fit_coeffs - 1 );

  const vector<double> coeffs = fit_log_log_polynomial( energies, effs, rel_uncerts, order );

  // A fit is "good" if no point is off by more than ~5% in efficiency; a fitted model that
  //  misrepresents the DRF is worse for Genie than plainly interpolating between the points.
  double max_abs_ln_resid = 0.0;
  bool coeffs_finite = true;
  for( const double c : coeffs )
    coeffs_finite = (coeffs_finite && std::isfinite(c));

  for( size_t i = 0; coeffs_finite && (i < energies.size()); ++i )
  {
    const double model_eff = eval_log_log_polynomial( coeffs, energies[i] );
    max_abs_ln_resid = std::max( max_abs_ln_resid,
                                 std::fabs( std::log(model_eff) - std::log(effs[i]) ) );
  }

  const double max_acceptable_ln_resid = 0.05; //~5% relative efficiency error
  const bool fit_is_good = coeffs_finite && (max_abs_ln_resid <= max_acceptable_ln_resid);

  if( fit_is_good )
  {
    answer.model = CAMInputOutput::CAMIO::EfficiencyModel::EMPIRICAL;
    answer.fit_coeffs.reserve( coeffs.size() );
    for( const double c : coeffs )
      answer.fit_coeffs.push_back( static_cast<float>(c) );
    answer.fit_reference_energy = CAMInputOutput::CAMIO::sm_geom_default_fit_ref_energy;
    answer.detector_name = drf.name();

    // Reduced chi-square of the fit, in ln(efficiency) space where the uncertainties live.
    //
    // Floored at 1: the uncertainties this is divided by are defaults rather than measurements
    //  (see `genie_default_eff_uncertainty(...)`), so a vanishing chi-square says only that the
    //  DRF happened to be an exact member of the fit's own functional family - not that the
    //  calibration is that good.  Genie never leaves this field zero, and by now we know better
    //  than to be the first to.
    const size_t dof = (energies.size() > coeffs.size()) ? (energies.size() - coeffs.size()) : 1;
    double chi2 = 0.0;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double resid = std::log( eval_log_log_polynomial(coeffs, energies[i]) )
                            - std::log( effs[i] );
      chi2 += (resid*resid) / (rel_uncerts[i]*rel_uncerts[i]);
    }
    answer.fit_chi_square = std::max( 1.0f, static_cast<float>(chi2 / dof) );
  }else
  {
    answer.model = CAMInputOutput::CAMIO::EfficiencyModel::INTERPOL;
    if( warnings )
      warnings->push_back( WString::tr("esfcam-warn-eff-not-fit").toUTF8() );
  }//if( fit_is_good ) / else

  for( size_t i = 0; i < energies.size(); ++i )
  {
    CAMInputOutput::EfficiencyPoint eff_pt;
    eff_pt.Index = static_cast<int>( i );
    eff_pt.Energy = static_cast<float>( energies[i] );
    eff_pt.Efficiency = static_cast<float>( effs[i] );
    eff_pt.EfficiencyUncertainty = static_cast<float>( effs[i] * rel_uncerts[i] );
    answer.points.push_back( eff_pt );
  }

  return answer;
}//convert_efficiency_to_genie(...)


GenieLibraryModel::GenieLibraryModel( Wt::WObject *parent )
  : Wt::WAbstractItemModel( parent )
{
}


void GenieLibraryModel::setSources( vector<GenieLibrarySource> sources )
{
  layoutAboutToBeChanged().emit();
  m_sources = std::move( sources );
  layoutChanged().emit();
}//void setSources(...)


bool GenieLibraryModel::anySourceIsAgeable() const
{
  return std::any_of( begin(m_sources), end(m_sources),
                      []( const GenieLibrarySource &s ){ return s.is_ageable; } );
}//anySourceIsAgeable()


WModelIndex GenieLibraryModel::index( int row, int column, const WModelIndex &parent ) const
{
  if( (column < 0) || (column >= static_cast<int>(Column::NumColumns)) || (row < 0) )
    return WModelIndex();

  if( !parent.isValid() )
  {
    if( row >= static_cast<int>(m_sources.size()) )
      return WModelIndex();
    return createIndex( row, column, uint64_t(0) );
  }

  // `parent` must itself be a source row (id==0); line rows have no children.
  if( parent.internalId() != 0 )
    return WModelIndex();

  const int source_index = parent.row();
  if( (source_index < 0) || (source_index >= static_cast<int>(m_sources.size())) )
    return WModelIndex();

  if( !m_sources[source_index].included ) //see rowCount(...)
    return WModelIndex();

  if( row >= static_cast<int>(m_sources[source_index].lines.size()) )
    return WModelIndex();

  return createIndex( row, column, uint64_t(source_index + 1) );
}//index(...)


WModelIndex GenieLibraryModel::parent( const WModelIndex &index ) const
{
  if( !index.isValid() || (index.internalId() == 0) )
    return WModelIndex();

  const int source_index = static_cast<int>( index.internalId() - 1 );
  return createIndex( source_index, 0, uint64_t(0) );
}//parent(...)


int GenieLibraryModel::rowCount( const WModelIndex &parent ) const
{
  if( !parent.isValid() )
    return static_cast<int>( m_sources.size() );

  if( parent.internalId() != 0 ) //parent is a line row - no children
    return 0;

  const int source_index = parent.row();
  if( (source_index < 0) || (source_index >= static_cast<int>(m_sources.size())) )
    return 0;

  // A source that isnt being written has no lines to show - leaving them visible (and checked)
  //  would tell the user lines are going into the file when they are not.
  if( !m_sources[source_index].included )
    return 0;

  return static_cast<int>( m_sources[source_index].lines.size() );
}//rowCount(...)


int GenieLibraryModel::columnCount( const WModelIndex & ) const
{
  return static_cast<int>( Column::NumColumns );
}//columnCount(...)


boost::any GenieLibraryModel::data( const WModelIndex &index, int role ) const
{
  if( !index.isValid() )
    return boost::any();

  const bool is_source_row = (index.internalId() == 0);
  const int source_index = is_source_row ? index.row() : static_cast<int>( index.internalId() - 1 );
  if( (source_index < 0) || (source_index >= static_cast<int>(m_sources.size())) )
    return boost::any();

  const GenieLibrarySource &source = m_sources[source_index];

  if( is_source_row )
  {
    switch( Column(index.column()) )
    {
      case Column::Name:
        if( role != Wt::DisplayRole )
          return boost::any();
        return boost::any( WString::fromUTF8(source.name) );

      case Column::Info:
        //  Just the time - the column header says what it is.
        if( role != Wt::DisplayRole )
          return boost::any();
        return boost::any( WString::fromUTF8( PhysicalUnits::printToBestTimeUnits(source.half_life_seconds) ) );

      case Column::Age:
        if( !source.is_ageable || ((role != Wt::DisplayRole) && (role != Wt::EditRole)) )
          return boost::any();
        return boost::any( WString::fromUTF8( PhysicalUnits::printToBestTimeUnits(source.age_seconds) ) );

      case Column::Key: //the key line belongs to a line row, not the source row
        return boost::any();

      case Column::Include:
        if( role != Wt::CheckStateRole )
          return boost::any();
        return boost::any( source.included );

      case Column::NumColumns:
        break;
    }//switch( Column(index.column()) )

    return boost::any();
  }//if( is_source_row )

  const int line_index = index.row();
  if( (line_index < 0) || (line_index >= static_cast<int>(source.lines.size())) )
    return boost::any();

  const GenieLibraryLine &line = source.lines[line_index];

  switch( Column(index.column()) )
  {
    case Column::Name:
      if( role != Wt::DisplayRole )
        return boost::any();
      {
        //  Just the number - the column header carries the "keV".
        char buffer[32];
        snprintf( buffer, sizeof(buffer), "%.2f", line.energy );
        return boost::any( WString::fromUTF8(buffer) );
      }

    case Column::Info:
    {
      if( role != Wt::DisplayRole )
        return boost::any();
      char buffer[32];
      snprintf( buffer, sizeof(buffer), "%.4g", 100.0*line.yield );
      WString answer = WString::tr("esfcam-line-yield-value").arg(buffer);
      if( line.is_xray )
        answer += WString::tr("esfcam-line-is-xray");
      return boost::any( answer );
    }

    case Column::Age:
      return boost::any();

    case Column::Key:
      //  A line that isnt being written cannot be the key line, so it gets no checkbox at all
      //  rather than an unchecked one the user could click to no effect.
      if( (role != Wt::CheckStateRole) || !line.included )
        return boost::any();
      return boost::any( line.is_key_line );

    case Column::Include:
      if( role != Wt::CheckStateRole )
        return boost::any();
      return boost::any( line.included );

    case Column::NumColumns:
      break;
  }//switch( Column(index.column()) )

  return boost::any();
}//data(...)


bool GenieLibraryModel::setData( const WModelIndex &index, const boost::any &value, int role )
{
  if( !index.isValid() )
    return false;

  const bool is_source_row = (index.internalId() == 0);
  const int source_index = is_source_row ? index.row() : static_cast<int>( index.internalId() - 1 );
  if( (source_index < 0) || (source_index >= static_cast<int>(m_sources.size())) )
    return false;

  GenieLibrarySource &source = m_sources[source_index];

  if( is_source_row )
  {
    switch( Column(index.column()) )
    {
      case Column::Include:
      {
        if( role != Wt::CheckStateRole )
          return false;

        const bool checked = boost::any_cast<bool>( value );
        if( checked == source.included )
          return true;

        const int num_lines = static_cast<int>( source.lines.size() );

        // Unchecking a source hides its lines entirely (see `rowCount(...)`), so the rows have to
        //  be removed/inserted properly rather than just repainted.
        //  Note: the parent index handed to begin{Remove,Insert}Rows has to be the source row's
        //  COLUMN 0 index - that is what `parent()` returns and what the view has registered as
        //  the parent node.  Passing this column's index (Include) silently does nothing.
        const WModelIndex source_row = this->index( index.row(), 0 );

        if( !checked && (num_lines > 0) )
          beginRemoveRows( source_row, 0, num_lines - 1 );
        else if( checked && (num_lines > 0) )
          beginInsertRows( source_row, 0, num_lines - 1 );

        source.included = checked;
        if( checked )
        {
          // Re-checking a source brings all of its lines back, checked.
          for( GenieLibraryLine &line : source.lines )
            line.included = true;
        }

        if( num_lines > 0 )
        {
          if( checked )
            endInsertRows();
          else
            endRemoveRows();
        }

        dataChanged().emit( index, index );
        return true;
      }

      case Column::Age:
      {
        if( !source.is_ageable || (role != Wt::EditRole) )
          return false;
        const string str_value = boost::any_cast<WString>(value).toUTF8();
        try
        {
          const double age = PhysicalUnits::stringToTimeDurationPossibleHalfLife(
                                                str_value, source.half_life_seconds );
          source.age_seconds = age / PhysicalUnits::second;
        }catch( std::exception & )
        {
          return false;
        }
        dataChanged().emit( index, index );
        m_ageEdited.emit( source.nuclide, source.age_seconds );
        return true;
      }

      default:
        return false;
    }//switch( Column(index.column()) )
  }//if( is_source_row )

  const int line_index = index.row();
  if( (line_index < 0) || (line_index >= static_cast<int>(source.lines.size())) )
    return false;

  const int num_lines = static_cast<int>( source.lines.size() );

  if( (Column(index.column()) == Column::Include) && (role == Wt::CheckStateRole) )
  {
    source.lines[line_index].included = boost::any_cast<bool>( value );

    // Un-writing the key line (or bringing a line back when the source has none) has to move the
    //  key line, so the whole source is repainted rather than just this row.
    if( ensureSingleKeyLine( source_index ) )
      dataChanged().emit( this->index( 0, 0, index.parent() ),
                          this->index( num_lines - 1, static_cast<int>(Column::NumColumns) - 1, index.parent() ) );
    else
      dataChanged().emit( index, index );

    return true;
  }

  if( (Column(index.column()) == Column::Key) && (role == Wt::CheckStateRole) )
  {
    GenieLibraryLine &line = source.lines[line_index];
    if( !line.included )
      return false;

    // Every nuclide needs exactly one key line, so this column is a radio group: unchecking the
    //  current key line is refused (repainted back to checked), and checking a line moves the key
    //  line off of whichever line held it.
    if( !boost::any_cast<bool>( value ) )
    {
      dataChanged().emit( index, index );
      return true;
    }

    for( GenieLibraryLine &l : source.lines )
      l.is_key_line = false;
    line.is_key_line = true;

    dataChanged().emit( this->index( 0, static_cast<int>(Column::Key), index.parent() ),
                        this->index( num_lines - 1, static_cast<int>(Column::Key), index.parent() ) );
    return true;
  }//if( key-line checkbox )

  return false;
}//setData(...)


bool GenieLibraryModel::ensureSingleKeyLine( const int source_index )
{
  if( (source_index < 0) || (source_index >= static_cast<int>(m_sources.size())) )
    return false;

  vector<GenieLibraryLine> &lines = m_sources[source_index].lines;

  int current_key = -1, best_included = -1;
  bool changed = false;
  for( size_t i = 0; i < lines.size(); ++i )
  {
    if( lines[i].is_key_line )
    {
      if( lines[i].included && (current_key < 0) )
      {
        current_key = static_cast<int>( i );
      }else
      {
        lines[i].is_key_line = false; //drop duplicates, and key lines that arent being written
        changed = true;
      }
    }

    if( lines[i].included
       && ((best_included < 0) || (lines[i].yield > lines[best_included].yield)) )
    {
      best_included = static_cast<int>( i );
    }
  }//for( size_t i = 0; i < lines.size(); ++i )

  if( (current_key < 0) && (best_included >= 0) )
  {
    lines[best_included].is_key_line = true;
    changed = true;
  }

  return changed;
}//ensureSingleKeyLine(...)


WFlags<ItemFlag> GenieLibraryModel::flags( const WModelIndex &index ) const
{
  if( !index.isValid() )
    return WFlags<ItemFlag>();

  const bool is_source_row = (index.internalId() == 0);

  if( Column(index.column()) == Column::Include )
    return WFlags<ItemFlag>( ItemFlag::ItemIsUserCheckable );

  if( !is_source_row && (Column(index.column()) == Column::Key) )
  {
    const int source_index = static_cast<int>( index.internalId() - 1 );
    const int line_index = index.row();
    if( (source_index >= 0) && (source_index < static_cast<int>(m_sources.size()))
       && (line_index >= 0) && (line_index < static_cast<int>(m_sources[source_index].lines.size()))
       && m_sources[source_index].lines[line_index].included )
    {
      return WFlags<ItemFlag>( ItemFlag::ItemIsUserCheckable );
    }
  }//if( a line rows key-line checkbox )

  if( is_source_row && (Column(index.column()) == Column::Age) )
  {
    const int source_index = index.row();
    if( (source_index >= 0) && (source_index < static_cast<int>(m_sources.size()))
       && m_sources[source_index].is_ageable )
    {
      return WFlags<ItemFlag>( ItemFlag::ItemIsEditable );
    }
  }

  return WFlags<ItemFlag>();
}//flags(...)


boost::any GenieLibraryModel::headerData( int section, Orientation orientation, int role ) const
{
  if( orientation != Horizontal )
    return WAbstractItemModel::headerData( section, orientation, role );

  if( role == Wt::ToolTipRole )
  {
    if( Column(section) == Column::Key )
      return boost::any( WString::tr("esfcam-col-key-tt") );
    return boost::any();
  }

  if( role != Wt::DisplayRole )
    return WAbstractItemModel::headerData( section, orientation, role );

  switch( Column(section) )
  {
    case Column::Name:        return boost::any( WString::tr("esfcam-col-source-energy") );
    case Column::Info:        return boost::any( WString::tr("esfcam-col-halflife") );
    case Column::Age:         return boost::any( WString::tr("esfcam-col-age") );
    case Column::Key:         return boost::any( WString::tr("esfcam-col-key") );
    case Column::Include:     return boost::any( WString::tr("esfcam-col-write") );
    case Column::NumColumns:  break;
  }

  return boost::any();
}//headerData(...)


Wt::Signal<const SandiaDecay::Nuclide *, double> &GenieLibraryModel::ageEdited()
{
  return m_ageEdited;
}


GenieCnfOptionsWidget::GenieCnfOptionsWidget( InterSpec *viewer, WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_interspec( viewer ),
    m_writeSpectrumCb( nullptr ),
    m_writePeaksCb( nullptr ),
    m_writeLibraryCb( nullptr ),
    m_libraryModeCb( nullptr ),
    m_thresholdLabel( nullptr ),
    m_thresholdEdit( nullptr ),
    m_combineLinesCb( nullptr ),
    m_libraryTable( nullptr ),
    m_libraryModel( nullptr ),
    m_writeFwhmCb( nullptr ),
    m_fwhmSourceCb( nullptr ),
    m_fitFwhmFromPeaksBtn( nullptr ),
    m_fwhmFitTxt( nullptr ),
    m_have_manual_fwhm( false ),
    m_manual_fwhm_coeffs( 0.0f, 0.0f ),
    m_writeEfficiencyCb( nullptr ),
    m_effDistanceLabel( nullptr ),
    m_effDistanceEdit( nullptr ),
    m_writeEnergyCalCb( nullptr ),
    m_warningsTxt( nullptr )
{
  assert( m_interspec );
  if( m_interspec )
  {
    m_interspec->useMessageResourceBundle( "ExportSpecFile" );

    // The DRF can change while this dialog is open (the "Fit from peaks..." tool can even create
    //  one), which decides whether an efficiency curve can be written at all.
    m_interspec->detectorChanged().connect( this, &GenieCnfOptionsWidget::handleDetectorChanged );
    m_interspec->detectorModified().connect( this, &GenieCnfOptionsWidget::handleDetectorChanged );
  }//if( m_interspec )

  addStyleClass( "GenieCnfOptions" );

  WText *title = new WText( WString::tr("esfcam-title"), this );
  title->addStyleClass( "ExportColTitle" );

  // --- What to write ---
  m_writeSpectrumCb = new WCheckBox( WString::tr("esfcam-write-spectrum"), this );
  m_writeSpectrumCb->addStyleClass( "CbNoLineBreak" );
  m_writeSpectrumCb->setChecked( true );
  HelpSystem::attachToolTipOn( m_writeSpectrumCb, WString::tr("esfcam-write-spectrum-tt"),
                               true, HelpSystem::ToolTipPosition::Right,
                               HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_writeSpectrumCb->checked().connect( this, &GenieCnfOptionsWidget::handleWriteSpectrumChanged );
  m_writeSpectrumCb->unChecked().connect( this, &GenieCnfOptionsWidget::handleWriteSpectrumChanged );

  m_writePeaksCb = new WCheckBox( WString::tr("esfcam-write-peaks"), this );
  m_writePeaksCb->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_writePeaksCb, WString::tr("esfcam-write-peaks-tt"),
                               true, HelpSystem::ToolTipPosition::Right,
                               HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_writePeaksCb->changed().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );

  // --- Nuclide library ---
  m_writeLibraryCb = new WCheckBox( WString::tr("esfcam-write-library"), this );
  HelpSystem::attachToolTipOn( m_writeLibraryCb, WString::tr("esfcam-write-library-tt"),
                               true, HelpSystem::ToolTipPosition::Right,
                               HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_writeLibraryCb->addStyleClass( "CbNoLineBreak" );
  m_writeLibraryCb->checked().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );
  m_writeLibraryCb->unChecked().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );

  m_libraryModeCb = new WComboBox( this );
  m_libraryModeCb->addItem( WString::tr("esfcam-lines-from-peaks") );
  m_libraryModeCb->addItem( WString::tr("esfcam-lines-above-thresh") );
  m_libraryModeCb->setCurrentIndex( 0 );
  m_libraryModeCb->activated().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );

  // Keep the label and its input on one line, so hiding them together leaves no stray gap.
  WContainerWidget *thresholdRow = new WContainerWidget( this );
  thresholdRow->addStyleClass( "GenieCnfRow" );
  m_thresholdLabel = new WLabel( WString::tr("esfcam-threshold-label"), thresholdRow );
  m_thresholdEdit = new NativeFloatSpinBox( thresholdRow );
  m_thresholdEdit->setRange( 0.0f, 100.0f );
  m_thresholdEdit->setValue( 1.0f );
  m_thresholdLabel->setBuddy( m_thresholdEdit );
  m_thresholdEdit->valueChanged().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );

  m_combineLinesCb = new WCheckBox( WString::tr("esfcam-combine-lines"), this );
  m_combineLinesCb->addStyleClass( "CbNoLineBreak" );
  m_combineLinesCb->setChecked( true );
  HelpSystem::attachToolTipOn( m_combineLinesCb, WString::tr("esfcam-combine-lines-tt"),
                               true, HelpSystem::ToolTipPosition::Right,
                               HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_combineLinesCb->checked().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );
  m_combineLinesCb->unChecked().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );

  m_libraryModel = new GenieLibraryModel( this );
  m_libraryModel->ageEdited().connect( this, &GenieCnfOptionsWidget::handleSourceAgeEdited );

  m_libraryTable = new RowStretchTreeView( this );
  m_libraryTable->setRootIsDecorated( true );
  //  Sorting is not useful here (the rows are already energy ordered within each source), and the
  //  sort indicators eat header width we need for the column labels.
  m_libraryTable->setSortingEnabled( false );
  m_libraryTable->addStyleClass( "GenieLibraryTable" );
  m_libraryTable->setModel( m_libraryModel );
  // These have to add up to less than the dialog column the panel lives in (widened to ~340px by
  //  the `ExportSpecFileToolCnf` style class), or the "Write" checkbox column - the whole point of
  //  the table - ends up off-screen behind a scrollbar the user is unlikely to find.
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Name), 95 );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Info), 75 );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Age), 65 );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Key), 35 );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Include), 45 );
  m_libraryTable->expanded().connect( this, &GenieCnfOptionsWidget::handleSourceExpanded );
  m_libraryTable->collapsed().connect( this, &GenieCnfOptionsWidget::handleSourceCollapsed );

  //  Directly under the table, so warnings and "why you cannot export" messages are next to what
  //  they are about, rather than at the bottom of a panel the user has to scroll to reach.
  m_warningsTxt = new WText( this );
  m_warningsTxt->addStyleClass( "GenieCnfWarnings" );
  m_warningsTxt->setTextAlignment( Wt::AlignmentFlag::AlignLeft );
  m_warningsTxt->hide();
  // Note: dont call `setHeight(...)`/`resize(...)` on a RowStretchTreeView (see its header) - the
  //  table's height comes from the `.GenieLibraryTable` rule in ExportSpecFile.css.

  // --- FWHM ---
  m_writeFwhmCb = new WCheckBox( WString::tr("esfcam-write-fwhm"), this );
  m_writeFwhmCb->addStyleClass( "CbNoLineBreak" );
  m_writeFwhmCb->checked().connect( this, &GenieCnfOptionsWidget::handleFwhmSourceChanged );
  m_writeFwhmCb->unChecked().connect( this, &GenieCnfOptionsWidget::handleFwhmSourceChanged );

  //  Contents are filled in by `updateAvailableFwhmSources()`, so only sources that can actually
  //  produce a shape calibration for the current spectrum are ever offered.
  m_fwhmSourceCb = new WComboBox( this );
  m_fwhmSourceCb->activated().connect( this, &GenieCnfOptionsWidget::handleFwhmSourceChanged );

  //  Above the button, so the equation that will actually be written is visible right where the
  //  user picks/fits it.
  m_fwhmFitTxt = new WText( this );
  m_fwhmFitTxt->addStyleClass( "GenieCnfFitResult" );
  m_fwhmFitTxt->hide();

  m_fitFwhmFromPeaksBtn = new WPushButton( WString::tr("esfcam-fit-from-peaks-btn"), this );
  m_fitFwhmFromPeaksBtn->addStyleClass( "LightButton" );
  m_fitFwhmFromPeaksBtn->clicked().connect( this, &GenieCnfOptionsWidget::handleFitFwhmFromPeaksClicked );
  m_fitFwhmFromPeaksBtn->hide();

  // --- Efficiency / Energy cal ---
  //  Note: the "Export" button is a plain link to a `Wt::WResource` (see `ExportSpecFileTool`),
  //  so following it does not round-trip pending form state to the server.  Every input in this
  //  widget must therefore have a change handler connected, otherwise the user's selection can
  //  still be un-sent when the file gets written.
  m_writeEfficiencyCb = new WCheckBox( WString::tr("esfcam-write-eff"), this );
  m_writeEfficiencyCb->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_writeEfficiencyCb, WString::tr("esfcam-write-eff-tt"),
                               true, HelpSystem::ToolTipPosition::Right,
                               HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_writeEfficiencyCb->changed().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );

  // Genie efficiency curves are absolute, so a far-field DRF needs the source distance the
  //  exported curve is for; hidden for fixed-geometry DRFs, where there is no distance to give.
  WContainerWidget *effDistRow = new WContainerWidget( this );
  effDistRow->addStyleClass( "GenieCnfRow" );
  m_effDistanceLabel = new WLabel( WString::tr("esfcam-eff-distance-label"), effDistRow );
  m_effDistanceEdit = new WLineEdit( "1 m", effDistRow );
  m_effDistanceEdit->addStyleClass( "GenieEffDistance" );
  m_effDistanceLabel->setBuddy( m_effDistanceEdit );
  HelpSystem::attachToolTipOn( m_effDistanceEdit, WString::tr("esfcam-eff-distance-tt"),
                               true, HelpSystem::ToolTipPosition::Right,
                               HelpSystem::ToolTipPrefOverride::AlwaysShow );

  WRegExpValidator *distValidator
              = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, m_effDistanceEdit );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  distValidator->setMandatory( true );
  m_effDistanceEdit->setValidator( distValidator );

  m_effDistanceEdit->changed().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );
  m_effDistanceEdit->enterPressed().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );
  m_effDistanceEdit->blurred().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );

  m_writeEnergyCalCb = new WCheckBox( WString::tr("esfcam-write-energy-cal"), this );
  m_writeEnergyCalCb->addStyleClass( "CbNoLineBreak" );
  m_writeEnergyCalCb->setChecked( true );
  HelpSystem::attachToolTipOn( m_writeEnergyCalCb, WString::tr("esfcam-write-energy-cal-tt"),
                               true, HelpSystem::ToolTipPosition::Right,
                               HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_writeEnergyCalCb->changed().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );

  updateAvailableFwhmSources();
  handleWriteSpectrumChanged();
  handleLibraryOptionsChanged();
  handleFwhmSourceChanged();
  handleEfficiencyOrEnergyCalChanged();
}//GenieCnfOptionsWidget constructor


void GenieCnfOptionsWidget::updateForFile( const shared_ptr<const SpecMeas> &spec, const set<int> &samples )
{
  // `ExportSpecFileTool::refreshSampleAndDetectorOptions()` calls this on nearly every change in
  //  the export dialog, so bail early when nothing we care about actually changed - otherwise we
  //  would throw away the user's FWHM/age selections, and their per-source/line library choices,
  //  every time they touched an unrelated option.
  //
  //  Note: compare *contents*, not the `shared_ptr` - `currentlySelectedFile()` builds a brand new
  //  `SpecMeas` on every call when the "+ Back." combination is selected, so pointer identity
  //  never holds there and everything would reset constantly.
  const bool same_spectrum = (spec == m_spec)
        || (spec && m_spec
             && (spec->filename() == m_spec->filename())
             && (spec->uuid() == m_spec->uuid())
             && (spec->detector() == m_spec->detector())
             && (spec->num_measurements() == m_spec->num_measurements()));

  if( same_spectrum && (samples == m_samples) )
  {
    // Peaks or the DRF can still have changed under us; keep the options consistent with them.
    m_spec = spec;
    updateAvailableFwhmSources();
    return;
  }

  m_spec = spec;
  m_samples = samples;
  m_expanded_sources.clear();
  m_nuclide_ages.clear();
  m_have_manual_fwhm = false;

  const shared_ptr<const DetectorPeakResponse> drf = spec ? spec->detector() : nullptr;
  const bool have_drf = !!drf && drf->isValid();

  m_writeEfficiencyCb->setEnabled( have_drf );
  if( !have_drf )
    m_writeEfficiencyCb->setChecked( false );

  bool energy_cal_ok = false;
  bool is_hpge = true;
  if( spec )
  {
    for( const shared_ptr<const SpecUtils::Measurement> &m : spec->measurements() )
    {
      if( !m || !m->energy_calibration() )
        continue;

      switch( m->energy_calibration()->type() )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::FullRangeFraction:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          energy_cal_ok = true;
          break;
        case SpecUtils::EnergyCalType::LowerChannelEdge:
        case SpecUtils::EnergyCalType::InvalidEquationType:
          break;
      }//switch( energy cal type )

      if( energy_cal_ok )
      {
        is_hpge = PeakFitUtils::is_high_res( m );
        break;
      }
    }//for( loop over measurements, looking for one with a usable energy cal )
  }//if( spec )

  m_energy_cal_ok = energy_cal_ok;
  m_writeEnergyCalCb->setChecked( energy_cal_ok );

  // Only offer the FWHM sources that can actually produce a shape calibration for this spectrum.
  updateAvailableFwhmSources();

  // Pick a reasonable default: from the DRF if it can be fit, else HPGe/NaI by resolution.
  const GenieFwhmSource preferred = (have_drf && drf->hasResolutionInfo())
                                      ? GenieFwhmSource::FromDrf
                                      : (is_hpge ? GenieFwhmSource::DefaultHPGe
                                                 : GenieFwhmSource::DefaultNaI);
  for( size_t i = 0; i < m_available_fwhm_sources.size(); ++i )
  {
    if( m_available_fwhm_sources[i] == preferred )
      m_fwhmSourceCb->setCurrentIndex( static_cast<int>(i) );
  }

  handleWriteSpectrumChanged();
  handleFwhmSourceChanged();
  handleEfficiencyOrEnergyCalChanged();

  rebuildLibraryTable();
}//void updateForFile(...)


void GenieCnfOptionsWidget::rebuildLibraryTable()
{
  const shared_ptr<const SpecMeas> spec = m_spec;
  m_libraryModel->setSources( {} );
  m_library_warnings.clear();

  if( !spec )
  {
    updateWarningsAndExportable();
    return;
  }

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = spec->peaks( m_samples );
  if( !peaks || peaks->empty() )
  {
    m_library_warnings.push_back( WString::tr("esfcam-warn-no-peaks").toUTF8() );
    updateWarningsAndExportable();
    return;
  }

  const GenieLibraryLineMode mode = (m_libraryModeCb->currentIndex() == 0)
                                    ? GenieLibraryLineMode::PeakLinesOnly
                                    : GenieLibraryLineMode::AllLinesAboveThreshold;

  // Remember what the user had un-checked, so changing an unrelated option (the mode, the
  //  threshold, an age, the FWHM source, ...) doesn't silently re-check everything.
  map<string,bool> prev_source_included;
  map<pair<string,int>,bool> prev_line_included;
  for( const GenieLibrarySource &src : m_libraryModel->sources() )
  {
    prev_source_included[src.name] = src.included;
    for( const GenieLibraryLine &line : src.lines )
      prev_line_included[ {src.name, static_cast<int>(std::round(100.0*line.energy))} ] = line.included;
  }

  pair<float,float> energy_range{ 0.0f, 0.0f };
  for( const shared_ptr<const SpecUtils::Measurement> &m : spec->measurements() )
  {
    if( m && m->energy_calibration() && m->energy_calibration()->valid() )
    {
      energy_range.first = m->gamma_energy_min();
      energy_range.second = m->gamma_energy_max();
      break;
    }
  }

  vector<string> warnings;
  vector<GenieLibrarySource> sources = build_genie_library( *peaks, mode, m_thresholdEdit->value(),
                                                            m_combineLinesCb->isChecked(),
                                                            currentFwhmCoefficients(), energy_range,
                                                            m_nuclide_ages, &warnings );

  for( GenieLibrarySource &src : sources )
  {
    const auto src_pos = prev_source_included.find( src.name );
    if( src_pos != end(prev_source_included) )
      src.included = src_pos->second;

    for( GenieLibraryLine &line : src.lines )
    {
      const auto line_pos = prev_line_included.find( {src.name, static_cast<int>(std::round(100.0*line.energy))} );
      if( line_pos != end(prev_line_included) )
        line.included = line_pos->second;
    }
  }//for( restore the user's selections )

  m_libraryModel->setSources( sources );

  // Nothing in this library can be aged (e.g. every source is a Cs137-like "age doesn't matter"
  //  nuclide), so the Age column is dead space.
  m_libraryTable->setColumnHidden( static_cast<int>(GenieLibraryModel::Column::Age),
                                   !m_libraryModel->anySourceIsAgeable() );

  // `setSources()` resets the model, which drops WTreeView's expanded set - re-expand what the
  //  user had open.
  for( int row = 0; row < m_libraryModel->rowCount(); ++row )
  {
    const WModelIndex index = m_libraryModel->index( row, 0 );
    const GenieLibrarySource &src = m_libraryModel->sources()[row];
    if( m_expanded_sources.count(src.name) )
      m_libraryTable->expand( index );
  }

  m_library_warnings = warnings;
  updateWarningsAndExportable();
}//void rebuildLibraryTable()


void GenieCnfOptionsWidget::handleLibraryOptionsChanged()
{
  const bool write_library = m_writeLibraryCb->isChecked();
  m_libraryModeCb->setHidden( !write_library );
  m_combineLinesCb->setHidden( !write_library );
  m_libraryTable->setHidden( !write_library );

  const bool show_threshold = write_library && (m_libraryModeCb->currentIndex() != 0);
  m_thresholdLabel->setHidden( !show_threshold );
  m_thresholdEdit->setHidden( !show_threshold );

  if( write_library )
    rebuildLibraryTable();
  else
    updateWarningsAndExportable();
}//void handleLibraryOptionsChanged()


void GenieCnfOptionsWidget::handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> /*drf*/ )
{
  const shared_ptr<const DetectorPeakResponse> drf = m_spec ? m_spec->detector() : nullptr;
  const bool have_drf = !!drf && drf->isValid();

  m_writeEfficiencyCb->setEnabled( have_drf );
  if( !have_drf )
    m_writeEfficiencyCb->setChecked( false );

  updateAvailableFwhmSources();
  handleFwhmSourceChanged();
  handleEfficiencyOrEnergyCalChanged();
}//handleDetectorChanged(...)


void GenieCnfOptionsWidget::handleSourceExpanded( const Wt::WModelIndex &index )
{
  if( index.isValid() && (index.internalId() == 0)
     && (index.row() < static_cast<int>(m_libraryModel->sources().size())) )
    m_expanded_sources.insert( m_libraryModel->sources()[index.row()].name );
}//handleSourceExpanded(...)


void GenieCnfOptionsWidget::handleSourceCollapsed( const Wt::WModelIndex &index )
{
  if( index.isValid() && (index.internalId() == 0)
     && (index.row() < static_cast<int>(m_libraryModel->sources().size())) )
    m_expanded_sources.erase( m_libraryModel->sources()[index.row()].name );
}//handleSourceCollapsed(...)


void GenieCnfOptionsWidget::handleWriteSpectrumChanged()
{
  // Writing the spectrum always writes its energy calibration - channel counts are meaningless
  //  without it - so the option is only offered for a library/calibration-only file.
  const bool writing_spectrum = writeSpectrum();
  m_writeEnergyCalCb->setHidden( writing_spectrum || !m_energy_cal_ok );

  updateWarningsAndExportable();
}//handleWriteSpectrumChanged()


void GenieCnfOptionsWidget::updateAvailableFwhmSources()
{
  const shared_ptr<const DetectorPeakResponse> drf = m_spec ? m_spec->detector() : nullptr;

  // Only offer sources that can actually produce a `FWHM = A0 + A1*sqrt(E)` right now.
  vector<GenieFwhmSource> available{ GenieFwhmSource::DefaultHPGe, GenieFwhmSource::DefaultNaI };

  if( drf && drf->isValid() && drf->hasResolutionInfo() )
  {
    // ...and only if it can actually be fit to Genie's equation form.
    try
    {
      fit_genie_fwhm_from_drf( *drf );
      available.push_back( GenieFwhmSource::FromDrf );
    }catch( std::exception & )
    {
    }
  }//if( the DRF has usable resolution info )

  // "From peaks..." needs peaks to fit to.
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks
                            = m_spec ? m_spec->peaks( m_samples ) : nullptr;
  if( peaks && !peaks->empty() )
    available.push_back( GenieFwhmSource::FromPeaks );

  if( available == m_available_fwhm_sources )
    return;

  // Keep the user's choice if it is still available, else fall back to a sensible default.
  GenieFwhmSource wanted = GenieFwhmSource::None;
  const int current = m_fwhmSourceCb->currentIndex();
  if( (current >= 0) && (current < static_cast<int>(m_available_fwhm_sources.size())) )
    wanted = m_available_fwhm_sources[current];

  m_available_fwhm_sources = available;
  m_fwhmSourceCb->clear();

  int wanted_index = -1;
  for( size_t i = 0; i < m_available_fwhm_sources.size(); ++i )
  {
    switch( m_available_fwhm_sources[i] )
    {
      case GenieFwhmSource::DefaultHPGe:
        m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-default-hpge") );  break;
      case GenieFwhmSource::DefaultNaI:
        m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-default-nai") );   break;
      case GenieFwhmSource::FromDrf:
        m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-from-drf") );      break;
      case GenieFwhmSource::FromPeaks:
        m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-from-peaks") );    break;
      case GenieFwhmSource::None:
        continue;
    }//switch( source )

    if( m_available_fwhm_sources[i] == wanted )
      wanted_index = static_cast<int>( i );
  }//for( fill the combo )

  m_fwhmSourceCb->setCurrentIndex( (wanted_index >= 0) ? wanted_index : 0 );
}//updateAvailableFwhmSources()


void GenieCnfOptionsWidget::updateWarningsAndExportable()
{
  vector<string> messages = m_library_warnings;

  WString reason;
  if( !canExport( &reason ) && !reason.empty() )
    messages.push_back( reason.toUTF8() );

  if( messages.empty() )
  {
    m_warningsTxt->hide();
  }else
  {
    string txt;
    for( const string &m : messages )
      txt += (txt.empty() ? "" : "<br />") + m;
    m_warningsTxt->setText( WString::fromUTF8(txt) );
    m_warningsTxt->show();
  }

  m_exportableChanged.emit();
}//updateWarningsAndExportable()


bool GenieCnfOptionsWidget::writeSpectrum() const
{
  return !m_writeSpectrumCb || m_writeSpectrumCb->isChecked();
}//writeSpectrum()


bool GenieCnfOptionsWidget::canExport( WString *reason ) const
{
  // A file has to contain *something*.
  if( !writeSpectrum() && !writeLibrary() && !writePeaks() && !writeFwhm() && !writeEfficiency() )
  {
    if( reason )
      *reason = WString::tr("esfcam-err-nothing-to-write");
    return false;
  }

  // A checked option whose input is missing/invalid must not be silently dropped.
  const shared_ptr<const DetectorPeakResponse> drf = m_spec ? m_spec->detector() : nullptr;
  if( writeEfficiency() && drf && drf->isValid() && !drf->isFixedGeometry()
     && (currentEfficiencyDistance() <= 0.0) )
  {
    if( reason )
      *reason = WString::tr("esfcam-err-bad-distance");
    return false;
  }

  if( writeLibrary() && m_libraryModel->sources().empty() )
  {
    if( reason )
      *reason = WString::tr("esfcam-err-no-library-lines");
    return false;
  }

  // "From peaks..." with no fit yet would silently fall back to the detector-type default rather
  //  than write what the user asked for.
  if( m_writeFwhmCb->isChecked() && (currentFwhmSource() == GenieFwhmSource::FromPeaks)
     && !m_have_manual_fwhm )
  {
    if( reason )
      *reason = WString::tr("esfcam-err-fwhm-not-fit");
    return false;
  }

  return true;
}//canExport(...)


Wt::Signal<> &GenieCnfOptionsWidget::exportableChanged()
{
  return m_exportableChanged;
}


void GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged()
{
  // (this handler also exists to force a round-trip; see its declaration)
  const shared_ptr<const DetectorPeakResponse> drf = m_spec ? m_spec->detector() : nullptr;
  const bool have_drf = !!drf && drf->isValid();

  // A distance is only needed to turn the DRF's intrinsic efficiency into the absolute
  //  efficiency Genie wants; a fixed-geometry DRF is already absolute.
  // No peaks fitted means nothing to write; dont offer the option.
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> fit_peaks
                              = m_spec ? m_spec->peaks( m_samples ) : nullptr;
  const bool have_peaks = !!fit_peaks && !fit_peaks->empty();
  m_writePeaksCb->setEnabled( have_peaks );
  if( !have_peaks )
    m_writePeaksCb->setChecked( false );

  const bool need_distance = writeEfficiency() && have_drf && !drf->isFixedGeometry();
  m_effDistanceLabel->setHidden( !need_distance );
  m_effDistanceEdit->setHidden( !need_distance );

  if( need_distance && (currentEfficiencyDistance() <= 0.0) )
    m_effDistanceEdit->addStyleClass( "InvalidInput" );
  else
    m_effDistanceEdit->removeStyleClass( "InvalidInput" );

  updateWarningsAndExportable();
}//void handleEfficiencyOrEnergyCalChanged()


void GenieCnfOptionsWidget::handleSourceAgeEdited( const SandiaDecay::Nuclide *nuc, double age_seconds )
{
  if( !nuc )
    return;
  m_nuclide_ages[nuc] = age_seconds * PhysicalUnits::second;
  rebuildLibraryTable();
}//void handleSourceAgeEdited(...)


void GenieCnfOptionsWidget::handleFwhmSourceChanged()
{
  const bool write_fwhm = m_writeFwhmCb->isChecked();
  m_fwhmSourceCb->setHidden( !write_fwhm );

  const GenieFwhmSource source = currentFwhmSource();
  m_fitFwhmFromPeaksBtn->setHidden( !write_fwhm || (source != GenieFwhmSource::FromPeaks) );

  // Show the equation that will actually be written, for whichever source is selected - the user
  //  should not have to guess what "Default NaI" or a peak fit amounts to.
  m_fwhmFitTxt->setHidden( !write_fwhm );
  if( write_fwhm )
  {
    if( (source == GenieFwhmSource::FromPeaks) && !m_have_manual_fwhm )
    {
      m_fwhmFitTxt->setText( WString::tr("esfcam-fwhm-not-fit-yet") );
    }else
    {
      const pair<float,float> coeffs = currentFwhmCoefficients();

      // Build the equation itself in code so the sign reads correctly; the surrounding sentence
      //  stays localizable.
      const string a0 = SpecUtils::printCompact( coeffs.first, 4 );
      const string a1 = SpecUtils::printCompact( std::fabs(coeffs.second), 4 );
      const string sign = (coeffs.second < 0.0f) ? " - " : " + ";
      const string eqn = a0 + sign + a1 + "\xE2\x88\x9A" "E";  //U+221A SQUARE ROOT

      m_fwhmFitTxt->setText( WString::tr("esfcam-fwhm-equation").arg( WString::fromUTF8(eqn) ) );
    }
  }//if( write_fwhm )

  // The library's line clustering uses the FWHM being written, so it has to be re-done.
  if( m_writeLibraryCb->isChecked() && m_combineLinesCb->isChecked() )
    rebuildLibraryTable();

  updateWarningsAndExportable();
}//void handleFwhmSourceChanged()


void GenieCnfOptionsWidget::handleFitFwhmFromPeaksClicked()
{
  // Go through InterSpec rather than `new MakeFwhmForDrfWindow(...)` directly: it guards against
  //  opening a second window, hides the window once a fit is accepted, and - crucially - connects
  //  `finished()` so the window is actually deleted.  `AuxWindow::setHidden()` does not call
  //  `WDialog::setHidden()`, so without that connection Wt never pops the modal dialog cover and
  //  the whole application becomes unclickable once this window is closed.
  if( !m_interspec )
    return;

  MakeFwhmForDrfWindow *window = m_interspec->fwhmFromForegroundWindow( true );
  MakeFwhmForDrf *tool = window ? window->tool() : nullptr;
  if( !tool )
    return;

  tool->restrictToConstantPlusSqrtEnergy();
  tool->updatedDrf().connect( this, &GenieCnfOptionsWidget::handleFwhmFitFromToolUpdated );
}//void handleFitFwhmFromPeaksClicked()


void GenieCnfOptionsWidget::handleFwhmFitFromToolUpdated( shared_ptr<DetectorPeakResponse> new_drf )
{
  if( !new_drf || (new_drf->resolutionFcnType() != DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy) )
    return;

  const vector<float> &coefs = new_drf->resolutionFcnCoefficients();
  if( coefs.size() != 2 )
    return;

  m_manual_fwhm_coeffs = make_pair( coefs[0], coefs[1] );
  m_have_manual_fwhm = true;

  // Show the newly fit equation, and let the export become possible now that we have one.
  handleFwhmSourceChanged();

  if( m_writeLibraryCb->isChecked() && m_combineLinesCb->isChecked() )
    rebuildLibraryTable();
}//void handleFwhmFitFromToolUpdated(...)


/** The `GenieFwhmSource` currently selected in the (dynamically populated) combo. */
GenieFwhmSource GenieCnfOptionsWidget::currentFwhmSource() const
{
  const int index = m_fwhmSourceCb->currentIndex();
  if( (index < 0) || (index >= static_cast<int>(m_available_fwhm_sources.size())) )
    return GenieFwhmSource::None;
  return m_available_fwhm_sources[index];
}//currentFwhmSource()


bool GenieCnfOptionsWidget::writeLibrary() const
{
  return m_writeLibraryCb->isChecked();
}


bool GenieCnfOptionsWidget::writePeaks() const
{
  return m_writePeaksCb->isEnabled() && m_writePeaksCb->isChecked();
}


bool GenieCnfOptionsWidget::writeFwhm() const
{
  return m_writeFwhmCb->isChecked() && (currentFwhmSource() != GenieFwhmSource::None);
}


bool GenieCnfOptionsWidget::writeEfficiency() const
{
  return m_writeEfficiencyCb->isEnabled() && m_writeEfficiencyCb->isChecked();
}


bool GenieCnfOptionsWidget::writeEnergyCal() const
{
  // Writing the spectrum always writes its calibration; only a library/calibration-only file
  //  gives the user a choice (and the checkbox is only shown in that case).
  if( writeSpectrum() )
    return true;

  return m_writeEnergyCalCb->isChecked();
}


pair<float,float> GenieCnfOptionsWidget::currentFwhmCoefficients() const
{
  const shared_ptr<const DetectorPeakResponse> drf = m_spec ? m_spec->detector() : nullptr;

  const GenieFwhmSource source = m_writeFwhmCb->isChecked() ? currentFwhmSource()
                                                            : GenieFwhmSource::None;

  switch( source )
  {
    case GenieFwhmSource::DefaultHPGe:
      return default_genie_fwhm( true );

    case GenieFwhmSource::DefaultNaI:
      return default_genie_fwhm( false );

    case GenieFwhmSource::FromDrf:
      if( drf && drf->isValid() && drf->hasResolutionInfo() )
        return fit_genie_fwhm_from_drf( *drf );
      break;

    case GenieFwhmSource::FromPeaks:
      if( m_have_manual_fwhm )
        return m_manual_fwhm_coeffs;
      break;

    case GenieFwhmSource::None:
      break;
  }//switch( source )

  // Nothing chosen (or the choice isnt usable yet) - fall back to the detector-type default,
  //  which is what `CAMIO::AddDetectorType(...)` will end up writing anyway.
  const bool is_hpge = !drf || !drf->isValid() || !drf->hasResolutionInfo()
                        || (drf->peakResolutionFWHM(661.657f) < 10.0f);
  return default_genie_fwhm( is_hpge );
}//currentFwhmCoefficients()


double GenieCnfOptionsWidget::currentEfficiencyDistance() const
{
  if( !m_effDistanceEdit || m_effDistanceEdit->isHidden() )
    return -1.0;

  try
  {
    const double distance = PhysicalUnits::stringToDistance( m_effDistanceEdit->text().toUTF8() );
    return (distance > 0.0) ? distance : -1.0;
  }catch( std::exception & )
  {
    return -1.0;
  }
}//currentEfficiencyDistance()


CAMInputOutput::CnfGenieExtras GenieCnfOptionsWidget::currentExtras() const
{
  CAMInputOutput::CnfGenieExtras extras;

  const shared_ptr<const SpecMeas> spec = m_spec;

  extras.omit_spectrum = !writeSpectrum();
  extras.omit_energy_calibration = !writeEnergyCal();

  // If the peaks carry a low-energy tail, tell GENIE about it - both as the per-peak `T` values
  //  (in `to_cam_peaks(...)`) and as the shape calibration's low-tail curve.
  if( spec )
  {
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = spec->peaks( m_samples );
    if( peaks && !peaks->empty() )
    {
      const std::optional<pair<float,float>> tail_cal = fit_genie_low_tail_cal( *peaks );
      if( tail_cal )
        extras.low_tail_cal.set( tail_cal->first, tail_cal->second );
    }
  }

  if( writeFwhm() )
  {
    // The combo only offers sources that can produce coefficients (see
    //  `updateAvailableFwhmSources()`), and "From peaks..." is only exportable once a fit has been
    //  made (see `canExport()`), so this always has something real to write.
    const GenieFwhmSource source = currentFwhmSource();
    if( (source != GenieFwhmSource::FromPeaks) || m_have_manual_fwhm )
    {
      const pair<float,float> coefs = currentFwhmCoefficients();
      extras.shape_cal.set( coefs.first, coefs.second );
    }
  }//if( writeFwhm() )

  // Kept beyond the block below so the peaks can record the same curve's value at their energies.
  std::unique_ptr<GenieEfficiencyResult> efficiency;

  if( writeEfficiency() && spec )
  {
    const shared_ptr<const DetectorPeakResponse> drf = spec->detector();
    if( drf && drf->isValid() )
    {
      // Genie efficiency curves are absolute, so a far-field DRF needs the user's distance;
      //  silently skipping is better than writing an efficiency that is wrong by the solid angle.
      //  `canExport()` refuses the export when a far-field DRF has no usable distance, so by the
      //  time we get here the distance is either irrelevant (fixed geometry) or valid.
      const double distance = drf->isFixedGeometry() ? 1.0 : currentEfficiencyDistance();
      efficiency.reset( new GenieEfficiencyResult( convert_efficiency_to_genie(*drf, distance) ) );
      extras.eff_model = efficiency->model;
      extras.eff_points = efficiency->points;
      extras.eff_fit_coeffs = efficiency->fit_coeffs;
      extras.eff_fit_reference_energy = efficiency->fit_reference_energy;
      extras.eff_fit_chi_square = efficiency->fit_chi_square;
      extras.eff_detector_name = efficiency->detector_name;
    }//if( drf && drf->isValid() )
  }//if( writeEfficiency() )

  if( writeLibrary() )
    extras.library_lines = to_library_lines( m_libraryModel->sources() );

  if( writePeaks() && spec )
  {
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = spec->peaks( m_samples );
    if( peaks && !peaks->empty() )
    {
      // Use the calibration and live time of the record actually being written.
      shared_ptr<const SpecUtils::EnergyCalibration> cal;
      float live_time = 0.0f;
      for( const int sample : m_samples )
      {
        for( const string &det : spec->detector_names() )
        {
          const shared_ptr<const SpecUtils::Measurement> m = spec->measurement( sample, det );
          if( !m )
            continue;
          if( !cal && m->energy_calibration() && m->energy_calibration()->valid() )
            cal = m->energy_calibration();
          live_time += m->live_time();
        }
      }//for( find the calibration and total live time )

      // The summed spectrum the peaks were fit to; a stepped continuum's area cannot be computed
      //  without it.  Summing here matches what `write_cnf(...)` writes as the file's spectrum.
      shared_ptr<const SpecUtils::Measurement> data;
      try
      {
        data = spec->sum_measurements( m_samples, spec->detector_names(), nullptr );
      }catch( std::exception &e )
      {
        cerr << "currentExtras: could not sum the spectrum for the peak continuum areas: "
             << e.what() << endl;
      }

      extras.peaks = to_cam_peaks( *peaks, cal, live_time, efficiency.get(), data );
    }//if( there are peaks )
  }//if( writePeaks() )

  return extras;
}//currentExtras(...)


}//namespace ExportSpecFileCAM
