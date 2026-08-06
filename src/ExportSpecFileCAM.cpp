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

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
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

      double energy_var_num = 0.0;
      for( size_t i = index; i < group_end; ++i )
        energy_var_num += lines[i].yield * std::pow( lines[i].energy - mean_energy, 2.0 );

      const size_t n = group_end - index;
      const double energy_unc = (yield_sum > 0.0)
                          ? std::sqrt( energy_var_num / (((static_cast<double>(n) - 1.0) / static_cast<double>(n)) * yield_sum) )
                          : 0.0;

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

    if( peak->xrayElement() )
    {
      if( warnings )
        warnings->push_back( WString::tr("esfcam-warn-orphan-xray-peak")
                             .arg( SpecUtils::printCompact(peak->mean(), 5) ).toUTF8() );
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

    double max_yield = 0.0;
    for( const SandiaDecay::EnergyRatePair &p : photons )
      max_yield = std::max( max_yield, p.numPerSecond );

    if( mode == GenieLibraryLineMode::PeakLinesOnly )
    {
      for( const shared_ptr<const PeakDef> &peak : nuc_peaks.second )
      {
        const bool is_xray = (peak->sourceGammaType() == PeakDef::SourceGammaType::XrayGamma);
        const float energy = is_xray ? peak->xrayEnergy() : peak->gammaParticleEnergy();

        const double yield = find_yield_near( photons, energy );
        if( yield <= 0.0 )
        {
          if( warnings )
            warnings->push_back( WString::tr("esfcam-warn-no-yield")
                                 .arg( source.name )
                                 .arg( SpecUtils::printCompact(energy, 5) ).toUTF8() );
          continue;
        }

        GenieLibraryLine line;
        line.energy = energy;
        line.yield = static_cast<float>( yield );
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

        const bool is_xray = std::any_of( begin(xrays), end(xrays),
                                          [&p]( const SandiaDecay::EnergyRatePair &x ){
          return std::fabs(x.energy - p.energy) < 0.001;
        } );

        GenieLibraryLine line;
        line.energy = static_cast<float>( p.energy );
        line.yield = static_cast<float>( p.numPerSecond );
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
      ll.yield = line.yield;
      ll.yield_uncert = line.yield_uncert;
      ll.no_weight_mean = line.is_xray;
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


pair<float,float> default_genie_fwhm( const bool is_hpge )
{
  // For HPGe, use the same Ge default `CAMIO::AddDetectorType(...)` writes, so an export with
  //  FWHM writing turned off and one with it left at the default agree.
  if( is_hpge )
    return make_pair( 1.0f, 0.035f );

  // The Genie manual's own NaI defaults ({-7.0, 2.0}) go negative below ~12 keV and noticeably
  //  under-predict the width at low energy, so instead fit the Genie equation form to
  //  InterSpec's nominal "NaI 3x3" resolution over the range NaI spectra are actually used over.
  static const pair<float,float> s_nai_coeffs
        = fit_genie_fwhm( []( const double energy ) -> double {
            return PeakFitUtils::nai_fwhm_fcn( static_cast<float>(energy) );
          }, 60.0, 2614.0 );

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

  return fit_genie_fwhm( [&drf]( const double energy ) -> double {
    return drf.peakResolutionFWHM( static_cast<float>(energy) );
  }, lower_energy, upper_energy );
}//fit_genie_fwhm_from_drf(...)


namespace
{
  /** Unweighted linear least-squares fit of `ln(y) = sum_i{ coeffs[i]*ln(x)^i }`. */
  vector<double> fit_log_log_polynomial( const vector<double> &x, const vector<double> &y,
                                         const size_t order )
  {
    const size_t n = x.size();
    const size_t num_coeffs = order + 1;
    Eigen::MatrixXd A( n, num_coeffs );
    Eigen::VectorXd b( n );

    for( size_t row = 0; row < n; ++row )
    {
      const double lnx = std::log( x[row] );
      double pow_term = 1.0;
      for( size_t col = 0; col < num_coeffs; ++col )
      {
        A(row,col) = pow_term;
        pow_term *= lnx;
      }
      b(row) = std::log( y[row] );
    }

    const Eigen::VectorXd sol = A.colPivHouseholderQr().solve( b );

    vector<double> coeffs( num_coeffs );
    for( size_t i = 0; i < num_coeffs; ++i )
      coeffs[i] = sol(static_cast<int>(i));
    return coeffs;
  }//fit_log_log_polynomial(...)


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


GenieEfficiencyResult convert_efficiency_to_genie( const DetectorPeakResponse &drf,
                                                   const double distance,
                                                   vector<string> *warnings )
{
  GenieEfficiencyResult answer;

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

  if( form == DetectorPeakResponse::EfficiencyFnctForm::kEnergyEfficiencyPairs )
  {
    answer.model = CAMInputOutput::CAMIO::EfficiencyModel::SPLINE;
    for( const DetectorPeakResponse::EnergyEfficiencyPair &pt : drf.getEnergyEfficiencyPair() )
    {
      CAMInputOutput::EfficiencyPoint eff_pt;
      eff_pt.Index = static_cast<int>( answer.points.size() );
      eff_pt.Energy = pt.energy;
      // Note: use the DRF's own tabulated energies, but go back through the DRF for the value, so
      //  the solid-angle (and any absolute-to-intrinsic) correction is applied consistently.
      eff_pt.Efficiency = static_cast<float>( absolute_eff( pt.energy ) );
      eff_pt.EfficiencyUncertainty = 0.0f;
      answer.points.push_back( eff_pt );
    }
    return answer;
  }//if( kEnergyEfficiencyPairs )

  // kExpOfLogPowerSeries (an exact analytic formula) and kFunctialEfficienyForm (an arbitrary
  // formula) are both handled by sampling points across the DRF's valid range: for
  // kExpOfLogPowerSeries the sampled points are exact (no fitting error), and Genie's DUAL model
  // - its own default for Ge/NaI detectors, and the same ln(eff)-vs-ln(energy) polynomial family
  // - is always tagged; for kFunctialEfficienyForm, DUAL is only tagged if a polynomial fit of
  // the sampled points is actually good, otherwise SPLINE (plain point interpolation) is used.
  const bool is_exact_log_power_series =
                  (form == DetectorPeakResponse::EfficiencyFnctForm::kExpOfLogPowerSeries);

  double lower_energy = drf.lowerEnergy();
  double upper_energy = drf.upperEnergy();
  if( !(upper_energy > lower_energy) )
  {
    lower_energy = 59.0;
    upper_energy = 2614.0;
  }
  lower_energy = std::max( lower_energy, 10.0 ); //avoid ln(0)/negative energies

  const int num_points = 15;
  vector<double> energies( num_points ), effs( num_points );
  for( int i = 0; i < num_points; ++i )
  {
    const double frac = static_cast<double>(i) / (num_points - 1);
    // log-spaced points, since efficiency curves are typically fit in log-log space
    const double energy = lower_energy * std::pow( upper_energy/lower_energy, frac );
    energies[i] = energy;
    effs[i] = std::max<double>( absolute_eff( energy ), 1.0E-9 );
  }

  bool use_dual = is_exact_log_power_series;
  if( !is_exact_log_power_series )
  {
    const size_t order = genie_default_poly_order( energies.size() );
    const vector<double> coeffs = fit_log_log_polynomial( energies, effs, order );

    double max_abs_ln_resid = 0.0;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      double model_ln_eff = 0.0, pow_term = 1.0;
      const double lnE = std::log( energies[i] );
      for( const double c : coeffs )
      {
        model_ln_eff += c * pow_term;
        pow_term *= lnE;
      }
      max_abs_ln_resid = std::max( max_abs_ln_resid, std::fabs(model_ln_eff - std::log(effs[i])) );
    }

    const double max_acceptable_ln_resid = 0.05; //~5% relative efficiency error
    use_dual = (max_abs_ln_resid <= max_acceptable_ln_resid);
    if( !use_dual && warnings )
    {
      warnings->push_back( WString::tr("esfcam-warn-eff-not-dual").toUTF8() );
    }
  }//if( !is_exact_log_power_series )

  answer.model = use_dual ? CAMInputOutput::CAMIO::EfficiencyModel::DUAL
                           : CAMInputOutput::CAMIO::EfficiencyModel::SPLINE;
  for( size_t i = 0; i < energies.size(); ++i )
  {
    CAMInputOutput::EfficiencyPoint eff_pt;
    eff_pt.Index = static_cast<int>( i );
    eff_pt.Energy = static_cast<float>( energies[i] );
    eff_pt.Efficiency = static_cast<float>( effs[i] );
    eff_pt.EfficiencyUncertainty = 0.0f;
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
        if( role != Wt::DisplayRole )
          return boost::any();
        return boost::any( WString::tr("esfcam-half-life-value")
                             .arg( PhysicalUnits::printToBestTimeUnits(source.half_life_seconds) ) );

      case Column::Age:
        if( !source.is_ageable || ((role != Wt::DisplayRole) && (role != Wt::EditRole)) )
          return boost::any();
        return boost::any( WString::fromUTF8( PhysicalUnits::printToBestTimeUnits(source.age_seconds) ) );

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
        char buffer[32];
        snprintf( buffer, sizeof(buffer), "%.2f", line.energy );
        return boost::any( WString::tr("esfcam-line-energy-value").arg(buffer) );
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
      if( line.is_key_line )
        answer += WString::tr("esfcam-line-is-key");
      return boost::any( answer );
    }

    case Column::Age:
      return boost::any();

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
        source.included = checked;
        for( GenieLibraryLine &line : source.lines )
          line.included = checked;
        const WModelIndex last_line = this->index( static_cast<int>(source.lines.size()) - 1,
                                                    static_cast<int>(Column::Include), index );
        dataChanged().emit( index, source.lines.empty() ? index : last_line );
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

  if( (Column(index.column()) == Column::Include) && (role == Wt::CheckStateRole) )
  {
    source.lines[line_index].included = boost::any_cast<bool>( value );
    dataChanged().emit( index, index );
    return true;
  }

  return false;
}//setData(...)


WFlags<ItemFlag> GenieLibraryModel::flags( const WModelIndex &index ) const
{
  if( !index.isValid() )
    return WFlags<ItemFlag>();

  const bool is_source_row = (index.internalId() == 0);

  if( Column(index.column()) == Column::Include )
    return WFlags<ItemFlag>( ItemFlag::ItemIsUserCheckable );

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
  if( (orientation != Horizontal) || (role != Wt::DisplayRole) )
    return WAbstractItemModel::headerData( section, orientation, role );

  switch( Column(section) )
  {
    case Column::Name:        return boost::any( WString::tr("esfcam-col-source-energy") );
    case Column::Info:        return boost::any( WString::tr("esfcam-col-hl-yield") );
    case Column::Age:         return boost::any( WString::tr("esfcam-col-age") );
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
    m_interspec->useMessageResourceBundle( "ExportSpecFile" );

  addStyleClass( "GenieCnfOptions" );

  WText *title = new WText( WString::tr("esfcam-title"), this );
  title->addStyleClass( "ExportColTitle" );

  // --- Nuclide library ---
  m_writeLibraryCb = new WCheckBox( WString::tr("esfcam-write-library"), this );
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
  m_combineLinesCb->setToolTip( WString::tr("esfcam-combine-lines-tt") );
  m_combineLinesCb->checked().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );
  m_combineLinesCb->unChecked().connect( this, &GenieCnfOptionsWidget::handleLibraryOptionsChanged );

  m_libraryModel = new GenieLibraryModel( this );
  m_libraryModel->ageEdited().connect( this, &GenieCnfOptionsWidget::handleSourceAgeEdited );

  m_libraryTable = new RowStretchTreeView( this );
  m_libraryTable->setRootIsDecorated( true );
  m_libraryTable->addStyleClass( "GenieLibraryTable" );
  m_libraryTable->setModel( m_libraryModel );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Name), 140 );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Info), 160 );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Age), 90 );
  m_libraryTable->setColumnWidth( static_cast<int>(GenieLibraryModel::Column::Include), 50 );
  // Note: dont call `setHeight(...)`/`resize(...)` on a RowStretchTreeView (see its header) - the
  //  table's height comes from the `.GenieLibraryTable` rule in ExportSpecFile.css.

  // --- FWHM ---
  m_writeFwhmCb = new WCheckBox( WString::tr("esfcam-write-fwhm"), this );
  m_writeFwhmCb->addStyleClass( "CbNoLineBreak" );
  m_writeFwhmCb->checked().connect( this, &GenieCnfOptionsWidget::handleFwhmSourceChanged );
  m_writeFwhmCb->unChecked().connect( this, &GenieCnfOptionsWidget::handleFwhmSourceChanged );

  m_fwhmSourceCb = new WComboBox( this );
  m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-none") );
  m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-default-hpge") );
  m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-default-nai") );
  m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-from-drf") );
  m_fwhmSourceCb->addItem( WString::tr("esfcam-fwhm-from-peaks") );
  m_fwhmSourceCb->activated().connect( this, &GenieCnfOptionsWidget::handleFwhmSourceChanged );

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
  m_writeEfficiencyCb->changed().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );

  // Genie efficiency curves are absolute, so a far-field DRF needs the source distance the
  //  exported curve is for; hidden for fixed-geometry DRFs, where there is no distance to give.
  WContainerWidget *effDistRow = new WContainerWidget( this );
  effDistRow->addStyleClass( "GenieCnfRow" );
  m_effDistanceLabel = new WLabel( WString::tr("esfcam-eff-distance-label"), effDistRow );
  m_effDistanceEdit = new WLineEdit( "1 m", effDistRow );
  m_effDistanceEdit->addStyleClass( "GenieEffDistance" );
  m_effDistanceLabel->setBuddy( m_effDistanceEdit );
  m_effDistanceEdit->setToolTip( WString::tr("esfcam-eff-distance-tt") );

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
  m_writeEnergyCalCb->changed().connect( this, &GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged );

  m_warningsTxt = new WText( this );
  m_warningsTxt->addStyleClass( "GenieCnfWarnings" );
  m_warningsTxt->setTextAlignment( Wt::AlignmentFlag::AlignLeft );
  m_warningsTxt->hide();

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
  if( (spec == m_spec) && (samples == m_samples) )
    return;

  m_spec = spec;
  m_samples = samples;
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

  // `SpecUtils::SpecFile::write_cnf(...)` always writes the summed spectrum, and channel counts
  //  are not interpretable without their energy calibration - so whenever a spectrum is being
  //  written the calibration goes with it unconditionally, and there is no choice to offer.
  //  The option is only meaningful for a spectrum-less CNF (e.g. a library-only export), which
  //  this dialog cannot currently produce; see `writeEnergyCal()`.
  const bool writing_spectrum = true;
  m_writeEnergyCalCb->setHidden( writing_spectrum || !energy_cal_ok );
  m_writeEnergyCalCb->setChecked( energy_cal_ok );

  // Pick a reasonable default FWHM source: from the DRF if it has resolution info, else
  // Default HPGe/NaI based on the spectrum's resolution.
  const int fwhm_index = have_drf && drf->hasResolutionInfo()
                          ? static_cast<int>(GenieFwhmSource::FromDrf)
                          : (is_hpge ? static_cast<int>(GenieFwhmSource::DefaultHPGe)
                                     : static_cast<int>(GenieFwhmSource::DefaultNaI));
  m_fwhmSourceCb->setCurrentIndex( fwhm_index );
  handleFwhmSourceChanged();
  handleEfficiencyOrEnergyCalChanged();

  rebuildLibraryTable();
}//void updateForFile(...)


void GenieCnfOptionsWidget::rebuildLibraryTable()
{
  const shared_ptr<const SpecMeas> spec = m_spec;
  m_libraryModel->setSources( {} );
  m_warningsTxt->hide();

  if( !spec )
    return;

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = spec->peaks( m_samples );
  if( !peaks || peaks->empty() )
    return;

  const GenieLibraryLineMode mode = (m_libraryModeCb->currentIndex() == 0)
                                    ? GenieLibraryLineMode::PeakLinesOnly
                                    : GenieLibraryLineMode::AllLinesAboveThreshold;

  vector<string> warnings;
  vector<GenieLibrarySource> sources = build_genie_library( *peaks, mode, m_thresholdEdit->value(),
                                                            m_combineLinesCb->isChecked(),
                                                            currentFwhmCoefficients(),
                                                            m_nuclide_ages, &warnings );

  m_libraryModel->setSources( sources );

  if( !warnings.empty() )
  {
    string txt;
    for( const string &w : warnings )
      txt += (txt.empty() ? "" : "<br />") + w;
    m_warningsTxt->setText( WString::fromUTF8(txt) );
    m_warningsTxt->show();
  }
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
}//void handleLibraryOptionsChanged()


void GenieCnfOptionsWidget::handleEfficiencyOrEnergyCalChanged()
{
  // (this handler also exists to force a round-trip; see its declaration)
  const shared_ptr<const DetectorPeakResponse> drf = m_spec ? m_spec->detector() : nullptr;
  const bool have_drf = !!drf && drf->isValid();

  // A distance is only needed to turn the DRF's intrinsic efficiency into the absolute
  //  efficiency Genie wants; a fixed-geometry DRF is already absolute.
  const bool need_distance = writeEfficiency() && have_drf && !drf->isFixedGeometry();
  m_effDistanceLabel->setHidden( !need_distance );
  m_effDistanceEdit->setHidden( !need_distance );

  if( need_distance && (currentEfficiencyDistance() <= 0.0) )
    m_effDistanceEdit->addStyleClass( "InvalidInput" );
  else
    m_effDistanceEdit->removeStyleClass( "InvalidInput" );
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

  const GenieFwhmSource source = static_cast<GenieFwhmSource>( m_fwhmSourceCb->currentIndex() );
  m_fitFwhmFromPeaksBtn->setHidden( !write_fwhm || (source != GenieFwhmSource::FromPeaks) );

  // The library's line clustering uses the FWHM being written, so it has to be re-done.
  if( m_writeLibraryCb->isChecked() && m_combineLinesCb->isChecked() )
    rebuildLibraryTable();
}//void handleFwhmSourceChanged()


void GenieCnfOptionsWidget::handleFitFwhmFromPeaksClicked()
{
  MakeFwhmForDrfWindow *window = new MakeFwhmForDrfWindow( true );
  MakeFwhmForDrf *tool = window->tool();
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

  if( m_writeLibraryCb->isChecked() && m_combineLinesCb->isChecked() )
    rebuildLibraryTable();
}//void handleFwhmFitFromToolUpdated(...)


bool GenieCnfOptionsWidget::writeLibrary() const
{
  return m_writeLibraryCb->isChecked();
}


bool GenieCnfOptionsWidget::writeFwhm() const
{
  return m_writeFwhmCb->isChecked()
         && (static_cast<GenieFwhmSource>(m_fwhmSourceCb->currentIndex()) != GenieFwhmSource::None);
}


bool GenieCnfOptionsWidget::writeEfficiency() const
{
  return m_writeEfficiencyCb->isEnabled() && m_writeEfficiencyCb->isChecked();
}


bool GenieCnfOptionsWidget::writeEnergyCal() const
{
  // Hidden means a spectrum is being written, and then the calibration always goes with it;
  //  see `updateForFile(...)`.
  return m_writeEnergyCalCb->isHidden() || m_writeEnergyCalCb->isChecked();
}


pair<float,float> GenieCnfOptionsWidget::currentFwhmCoefficients() const
{
  const shared_ptr<const DetectorPeakResponse> drf = m_spec ? m_spec->detector() : nullptr;

  const GenieFwhmSource source = m_writeFwhmCb->isChecked()
                    ? static_cast<GenieFwhmSource>( m_fwhmSourceCb->currentIndex() )
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

  extras.omit_energy_calibration = !writeEnergyCal();

  if( writeFwhm() )
  {
    const pair<float,float> coeffs = currentFwhmCoefficients();
    const GenieFwhmSource source = static_cast<GenieFwhmSource>( m_fwhmSourceCb->currentIndex() );

    // `currentFwhmCoefficients()` falls back to a default when the chosen source cant supply
    //  coefficients (e.g. "From peaks..." before a fit has been made); only write a shape
    //  calibration when the user's actual choice produced one.
    const bool have_coeffs = (source == GenieFwhmSource::DefaultHPGe)
                             || (source == GenieFwhmSource::DefaultNaI)
                             || ((source == GenieFwhmSource::FromDrf)
                                  && spec && spec->detector() && spec->detector()->isValid()
                                  && spec->detector()->hasResolutionInfo())
                             || ((source == GenieFwhmSource::FromPeaks) && m_have_manual_fwhm);
    if( have_coeffs )
      extras.shape_cal = coeffs;
  }//if( writeFwhm() )

  if( writeEfficiency() && spec )
  {
    const shared_ptr<const DetectorPeakResponse> drf = spec->detector();
    if( drf && drf->isValid() )
    {
      // Genie efficiency curves are absolute, so a far-field DRF needs the user's distance;
      //  silently skipping is better than writing an efficiency that is wrong by the solid angle.
      const double distance = drf->isFixedGeometry() ? 1.0 : currentEfficiencyDistance();
      if( drf->isFixedGeometry() || (distance > 0.0) )
      {
        try
        {
          const GenieEfficiencyResult eff = convert_efficiency_to_genie( *drf, distance );
          extras.eff_model = eff.model;
          extras.eff_points = eff.points;
        }catch( std::exception & )
        {
          //Leave the efficiency out rather than write something meaningless.
        }
      }//if( we have a usable geometry )
    }//if( drf && drf->isValid() )
  }//if( writeEfficiency() )

  if( writeLibrary() )
    extras.library_lines = to_library_lines( m_libraryModel->sources() );

  return extras;
}//currentExtras(...)


}//namespace ExportSpecFileCAM
