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
#include <mutex>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include "SandiaDecay.h"

#include "cascade/AnalyticCascade.h"
#include "cascade/SandiaDecayCascade.h"
#include "cascade/CascadeSummingProgram.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/CascadeSummingCalc.h"
#include "InterSpec/RelActCalcCascade.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/RelActCalcCascade_imp.hpp"

using namespace std;

namespace
{
  /** A fit line and a cascade member are the same line when their energies agree to this (keV); both come
   from the same decay data, so they agree exactly - this only allows for round-off.
   */
  const double ns_line_match_kev = 1.0E-3;

  /** Sums of photon pairs closer than this (keV) are one sum peak. */
  const double ns_sum_merge_kev = 0.2;

  /** A coincident photon whose detection could not change any line's summing-out factor by this much, even
   with the source touching the detector and nothing attenuated, is left out of the evaluation.
   */
  const double ns_min_partner_effect = 1.0E-7;

  /** A sum peak is made of its strongest photon pairs, down to this fraction of its (upper bound) rate left
   out: a peak can have thousands of pairs (every x-ray line of every converting transition, with every
   gamma), nearly all of them negligible.
   */
  const double ns_sum_peak_dropped_fraction = 1.0E-2;

  /** Below this energy (keV), a detector response that does not state its own energy range is not trusted. */
  const double ns_default_lower_energy_kev = 40.0;

  /** A response's total efficiency is taken as unphysical when it is above this, or below this fraction of
   the full-energy efficiency.
   */
  const double ns_max_physical_total_eff = 1.02;
  const double ns_min_total_to_fep = 0.9;

  /** Total-to-peak ratio of a small coaxial HPGe crystal (50 mm diameter x 30 mm, an ORTEC Detective), from a
   CeeLo Monte-Carlo response of it: the small-crystal end of the fall-back used for a response without a total
   efficiency (the large-crystal end is InterSpec's example response).  Energies in keV.
   */
  const double ns_small_crystal_energy[] = { 60.0, 93.35, 143.76, 185.71, 258.26, 400.0, 583.19, 766.37, 1000.99,
                                             1460.0, 2614.5 };
  const double ns_small_crystal_total_to_peak[] = { 1.177, 1.193, 1.400, 1.620, 2.131, 3.214, 4.489, 5.597, 6.791,
                                                    8.760, 13.98 };
  static_assert( (sizeof(ns_small_crystal_energy)/sizeof(ns_small_crystal_energy[0]))
                 == (sizeof(ns_small_crystal_total_to_peak)/sizeof(ns_small_crystal_total_to_peak[0])),
                 "one total-to-peak ratio per energy" );

  /** The energies (keV) a response's total efficiency is checked at, to decide whether it is physical; those
   outside the response's stated range are not used.
   */
  const double ns_total_check_energies[] = { 45.0, 60.0, 80.0, 100.0, 130.0, 170.0, 220.0, 300.0, 400.0, 550.0,
                                             750.0, 1000.0, 1400.0, 2000.0, 2700.0 };

  /** The two energies (keV) whose full-energy efficiency ratio stands for the crystal's size, and that ratio's
   log for the small crystal above.
   */
  const double ns_size_low_energy = 185.71, ns_size_high_energy = 1000.99;
  const double ns_small_crystal_log_fep_ratio = -1.999;

  /** Reference energy (keV) a geometry-aware response's near-field total efficiency is read at. */
  const double ns_solid_angle_ref_energy = 600.0;


  /** Log-log interpolation of the small crystal's total-to-peak ratio, held constant past its ends. */
  double small_crystal_total_to_peak( const double energy )
  {
    const size_t n = sizeof(ns_small_crystal_energy)/sizeof(ns_small_crystal_energy[0]);
    if( energy <= ns_small_crystal_energy[0] )
      return ns_small_crystal_total_to_peak[0];
    if( energy >= ns_small_crystal_energy[n-1] )
      return ns_small_crystal_total_to_peak[n-1];
    size_t i = 1;
    while( ns_small_crystal_energy[i] < energy )
      ++i;
    const double frac = std::log( energy/ns_small_crystal_energy[i-1] )
                        / std::log( ns_small_crystal_energy[i]/ns_small_crystal_energy[i-1] );
    return std::exp( (1.0 - frac)*std::log( ns_small_crystal_total_to_peak[i-1] )
                     + frac*std::log( ns_small_crystal_total_to_peak[i] ) );
  }//small_crystal_total_to_peak(...)


  /** Fractional solid angle (of 4 pi) of a disk of radius `disk_radius`, from a point a height `height` above
   its plane and `offset` from its axis.

   Integrated over the direction, in the disk's plane, from the foot of the perpendicular: the disk spans
   distances rho_1 to rho_2 in a direction, which subtend `height/sqrt(height^2 + rho_1^2)` less the same for
   rho_2.  The integrand is smooth however near the plane the point is (a quadrature over the disk's area is
   not: its integrand peaks ever more sharply under the point), and the limit at the plane is right: a half
   over the disk, nothing beside it.
   */
  double off_axis_solid_angle( const double offset, const double height, const double disk_radius )
  {
    if( !(disk_radius > 0.0) || !(height >= 0.0) )
      return 0.0;

    const auto cosine = [height]( const double rho ) -> double {
      const double dist = std::sqrt( height*height + rho*rho );
      return (dist > 0.0) ? (height/dist) : 1.0;
    };

    const int num_steps = 256;
    double sum = 0.0;
    if( offset < disk_radius )
    {
      // Under the disk: every direction reaches its edge, from rho = 0.
      for( int i = 0; i < num_steps; ++i )
      {
        const double phi = (i + 0.5)*2.0*PhysicalUnits::pi/num_steps;
        const double edge = -offset*std::cos( phi )
                            + std::sqrt( std::max( 0.0, disk_radius*disk_radius - std::pow( offset*std::sin( phi ), 2 ) ) );
        sum += 1.0 - cosine( edge );
      }
      return sum*(2.0*PhysicalUnits::pi/num_steps) / (4.0*PhysicalUnits::pi);
    }//if( under the disk )

    // Beside the disk: only directions within `phi_max` of its centre cross it.  phi = phi_max*sin(t) takes
    //  out the square-root behaviour of the chord at the tangent directions.
    const double phi_max = std::asin( std::min( 1.0, disk_radius/offset ) );
    for( int i = 0; i < num_steps; ++i )
    {
      const double t = -0.5*PhysicalUnits::pi + (i + 0.5)*PhysicalUnits::pi/num_steps;
      const double phi = phi_max*std::sin( t );
      const double half_chord = std::sqrt( std::max( 0.0, disk_radius*disk_radius - std::pow( offset*std::sin( phi ), 2 ) ) );
      const double mid = offset*std::cos( phi );
      sum += (cosine( std::max( 0.0, mid - half_chord ) ) - cosine( mid + half_chord )) * phi_max*std::cos( t );
    }
    return sum*(PhysicalUnits::pi/num_steps) / (4.0*PhysicalUnits::pi);
  }//off_axis_solid_angle(...)


  /** The prepared cascade model of `nuc` at `age`, aged from a pure parent.  Building one walks the whole
   decay chain, so they are cached - by the exact age: rates per current decay are formed with the caller's
   age, so a model of a nearby age would give answers that depend on which age was asked for first.
   */
  shared_ptr<const ceelo::AnalyticCascadeModel> cascade_model( const SandiaDecay::Nuclide * const nuc,
                                                              const double age )
  {
    typedef pair<const SandiaDecay::Nuclide *,double> Key;
    static std::mutex s_mutex;
    static map<Key,shared_ptr<const ceelo::AnalyticCascadeModel>> s_models;
    const size_t max_cached = 64;

    const Key key( nuc, age );
    {
      std::lock_guard<std::mutex> lock( s_mutex );
      const map<Key,shared_ptr<const ceelo::AnalyticCascadeModel>>::const_iterator pos = s_models.find( key );
      if( pos != end(s_models) )
        return pos->second;
    }

    ceelo::cascade_adapter::CascadeOptions opts;
    opts.age_seconds = age / PhysicalUnits::second;
    opts.prompt_equilibrium = false;

    const shared_ptr<const ceelo::AnalyticCascadeModel> model
          = make_shared<const ceelo::AnalyticCascadeModel>( ceelo::cascade_adapter::build_cascades( nuc, opts ) );

    std::lock_guard<std::mutex> lock( s_mutex );
    if( s_models.size() >= max_cached )
      s_models.clear();
    return s_models.emplace( key, model ).first->second;
  }//cascade_model(...)


  /** The energies (keV) the per-response tables are made at, for interpolating to a model's own energies:
   20 keV to 3 MeV, 5% apart.
   */
  vector<double> response_grid_energies()
  {
    const size_t num_energies = 96;
    const double lowest = 20.0, highest = 3000.0;
    vector<double> energies( num_energies );
    for( size_t i = 0; i < num_energies; ++i )
      energies[i] = lowest * std::pow( highest/lowest, static_cast<double>(i)/(num_energies - 1) );
    return energies;
  }//response_grid_energies()


  /** The table of the photons that scatter in an attenuator and still reach the detector (see
   `ShieldScatterAugment`), for a detector response, at #response_grid_energies; null if GADRAS's scatter
   data can not be loaded.  Built once per response (a few tenths of a second).
   */
  shared_ptr<const GammaInteractionCalc::ShieldScatterAugment> scatter_grid( const shared_ptr<const DetectorPeakResponse> &drf )
  {
    typedef shared_ptr<const GammaInteractionCalc::ShieldScatterAugment> TablePtr;
    static std::mutex s_mutex;
    static map<uint64_t,TablePtr> s_tables;
    const size_t max_cached = 8;

    const uint64_t key = drf->hashValue();
    {
      std::lock_guard<std::mutex> lock( s_mutex );
      const map<uint64_t,TablePtr>::const_iterator pos = s_tables.find( key );
      if( pos != end(s_tables) )
        return pos->second;
    }

    const vector<double> energies = response_grid_energies();

    // The detector's intrinsic total efficiency at any energy: log-log between the grid's, and zero below
    //  the response's range, as everywhere here.
    const RelActCalcCascade::DetectorEfficiencies effs = RelActCalcCascade::detector_efficiencies( drf, energies );
    const std::function<double(double)> total_intrinsic = [&energies, &effs]( const double energy ) -> double {
      if( energy < effs.lower_energy )
        return 0.0;
      const vector<double>::const_iterator pos = std::lower_bound( begin(energies), end(energies), energy );
      if( pos == begin(energies) )
        return effs.total.front();
      if( pos == end(energies) )
        return effs.total.back();
      const size_t upper = static_cast<size_t>( pos - begin(energies) ), lower = upper - 1;
      if( !(effs.total[lower] > 0.0) || !(effs.total[upper] > 0.0) )
        return std::max( effs.total[lower], effs.total[upper] );
      const double frac = std::log( energy/energies[lower] ) / std::log( energies[upper]/energies[lower] );
      return std::exp( (1.0 - frac)*std::log( effs.total[lower] ) + frac*std::log( effs.total[upper] ) );
    };//total_intrinsic

    TablePtr table;
    try
    {
      const shared_ptr<GammaInteractionCalc::ShieldScatterAugment> built
                                            = make_shared<GammaInteractionCalc::ShieldScatterAugment>();
      // No hydrogen: the attenuators of a Physical model are metals and oxides, and it is second order.
      built->build( energies, total_intrinsic, 0.0, InterSpec::staticDataDirectory() );
      if( built->valid() )
        table = built;
    }catch( std::exception & )
    {
      // No scatter data.  Not remembered: the data directory may be set later.
      return nullptr;
    }

    std::lock_guard<std::mutex> lock( s_mutex );
    if( s_tables.size() >= max_cached )
      s_tables.clear();
    s_tables[key] = table;
    return table;
  }//scatter_grid(...)


  /** Index of `energy` in the sorted `energies`; throws if it is not there. */
  size_t index_of_energy( const vector<double> &energies, const double energy )
  {
    const vector<double>::const_iterator pos = std::lower_bound( begin(energies), end(energies), energy );
    if( (pos == end(energies)) || (*pos != energy) )
      throw std::logic_error( "RelActCalcCascade: energy " + std::to_string(energy) + " keV not tabulated." );
    return static_cast<size_t>( pos - begin(energies) );
  }//index_of_energy(...)


  void sort_unique( vector<double> &values )
  {
    std::sort( begin(values), end(values) );
    values.erase( std::unique( begin(values), end(values) ), end(values) );
  }
}//namespace


namespace RelActCalcCascade
{

const DepthNodes &depth_nodes()
{
  // Panels shrinking by a factor of four toward the front face, 4-point Gauss-Legendre on each: integrates
  //  exp(-mu*AD*x), and the ratios of such integrals the kernel forms, to 3 parts in 10^4 for any mu*AD
  //  (checked in test_RelActCalcCascade).  Three points per panel are six times worse; panels shrinking by
  //  three, with four points, six times better for four more nodes.
  static const DepthNodes s_nodes = []() -> DepthNodes {
    const int npoints = 4;
    const double gl_x[npoints] = { -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 };
    const double gl_w[npoints] = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };

    vector<double> edges( 1, 0.0 );
    for( int power = 5; power >= 0; --power )
      edges.push_back( std::pow( 0.25, power ) );

    DepthNodes nodes;
    for( size_t panel = 0; (panel + 1) < edges.size(); ++panel )
    {
      const double lower = edges[panel], upper = edges[panel + 1];
      for( int i = 0; i < npoints; ++i )
      {
        nodes.fraction.push_back( 0.5*(upper + lower) + 0.5*(upper - lower)*gl_x[i] );
        nodes.weight.push_back( 0.5*(upper - lower)*gl_w[i] );
      }
    }//for( loop over panels )

    return nodes;
  }();

  return s_nodes;
}//depth_nodes()


DetectorEfficiencies detector_efficiencies( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                            const std::vector<double> &energies_kev )
{
  if( !drf || !drf->isValid() )
    throw runtime_error( "Cascade summing needs a valid detector response." );

  if( drf->isFixedGeometry() )
    throw runtime_error( "Cascade summing is not supported for fixed-geometry detector responses." );

  DetectorEfficiencies answer;
  answer.detector_diameter = drf->detectorDiameter();
  answer.detector_setback = drf->detectorSetback();
  if( !(answer.detector_diameter > 0.0) )
    throw runtime_error( "Cascade summing needs a detector response with a detector diameter." );

  answer.crystal_diameter = drf->modeledCrystalDiameter();
  answer.crystal_setback = drf->modeledCrystalSetback();
  if( !(answer.crystal_diameter > 0.0) )
  {
    answer.crystal_diameter = answer.detector_diameter;
    answer.crystal_setback = answer.detector_setback;
  }

  const size_t num_energies = energies_kev.size();
  answer.fep.assign( num_energies, 0.0 );
  answer.total.assign( num_energies, 0.0 );

  // Below its characterized range a response is an extrapolation - commonly of an exponential of a polynomial
  //  in log-energy, which can be anything; photons that soft are taken as not detected.
  const double stated_lower = drf->lowerEnergy(), stated_upper = drf->upperEnergy();
  const bool states_range = (stated_upper > stated_lower);
  const double lower_energy = states_range ? stated_lower : ns_default_lower_energy_kev;
  answer.lower_energy = lower_energy;

  for( size_t i = 0; i < num_energies; ++i )
  {
    if( energies_kev[i] < lower_energy )
      continue;
    const double eff = drf->farFieldIntrinsicEfficiency( static_cast<float>(energies_kev[i]) );
    answer.fep[i] = (std::isfinite(eff) ? std::max( 0.0, std::min( 1.0, eff ) ) : 0.0);
  }

  using GammaInteractionCalc::CascadeSummingCalc;

  bool use_response_total = drf->hasAnyTotalEfficiencyInfo();
  if( !use_response_total )
  {
    answer.total_source_reason = "the detector response has no total efficiency";
  }else
  {
    // A Monte-Carlo response evaluates its total efficiency by quadrature, and there can be a thousand
    //  partner energies: go through the shared interpolation table (within 0.1% above 60 keV).
    const std::function<double(double)> total_fcn = CascadeSummingCalc::interpolatedEff(
      "RelActCalcCascade-total-" + std::to_string( drf->hashValue() ),
      [drf]( double energy ) -> double {
        try
        {
          return drf->totalIntrinsicEfficiencyAny( static_cast<float>(energy) );
        }catch( std::exception & )
        {
          return 0.0;
        }
      } );

    // Whether the response's total is physical is decided on a fixed set of energies, not the ones asked
    //  for, so that every use of a response decides the same way.
    for( const double energy : ns_total_check_energies )
    {
      if( (energy < lower_energy) || (states_range && (energy > stated_upper)) )
        continue;

      const double fep = drf->farFieldIntrinsicEfficiency( static_cast<float>(energy) );
      const double total = total_fcn( energy );
      if( !std::isfinite(fep) || (fep <= 0.0) )
        continue;

      if( !std::isfinite(total) || (total > ns_max_physical_total_eff) || (total < ns_min_total_to_fep*std::min( 1.0, fep )) )
      {
        use_response_total = false;
        answer.total_source_reason = "the detector response's total efficiency is not physical at "
                                     + std::to_string( static_cast<int>(energy) ) + " keV";
        break;
      }
    }//for( loop over the energies checked )

    for( size_t i = 0; use_response_total && (i < num_energies); ++i )
    {
      if( answer.fep[i] <= 0.0 )
        continue;

      const double total = total_fcn( energies_kev[i] );
      answer.total[i] = std::isfinite( total ) ? std::max( answer.fep[i], std::min( 1.0, total ) ) : answer.fep[i];
    }//for( loop over energies )
  }//if( !use_response_total ) / else

  if( use_response_total )
  {
    answer.total_source = TotalEffSource::DetectorResponse;
    return answer;
  }

  // Fall back to a reference total-to-peak ratio, on this response's own full-energy efficiency: between that
  //  of a small crystal and of the (large) example detector, by where this response's full-energy efficiency
  //  ratio between a low and a high energy - which is set by the crystal's size - falls between theirs.
  const shared_ptr<const DetectorPeakResponse> example = CascadeSummingCalc::exampleDetector();
  if( !example )
    throw runtime_error( "Cascade summing: " + answer.total_source_reason
                         + ", and the example detector response could not be loaded." );

  answer.total_source = TotalEffSource::ExampleDetectorRatio;
  const double far_distance = 100.0*PhysicalUnits::cm;
  const double example_lower_kev = 35.0;  //the example response is not characterized below this
  const std::function<double(double)> ratio_fcn = CascadeSummingCalc::interpolatedEff(
    "RelActCalcCascade-example-total-to-peak",
    [example, far_distance]( double energy ) -> double {
      const double ex_fep = CascadeSummingCalc::detectorFepEffAbs( example, energy, 0.0, 0.0, far_distance );
      const double ex_tot = CascadeSummingCalc::detectorTotEffAbs( example, energy, 0.0, 0.0, far_distance );
      return ((ex_fep > 0.0) && (ex_tot > ex_fep)) ? (ex_tot / ex_fep) : 1.0;
    } );

  // How large a crystal: 0 for the small reference, 1 for the example detector (or larger).
  double large_weight = 1.0;
  {
    const double example_low = CascadeSummingCalc::detectorFepEffAbs( example, ns_size_low_energy, 0.0, 0.0, far_distance );
    const double example_high = CascadeSummingCalc::detectorFepEffAbs( example, ns_size_high_energy, 0.0, 0.0, far_distance );
    const double drf_low = drf->farFieldIntrinsicEfficiency( static_cast<float>(ns_size_low_energy) );
    const double drf_high = drf->farFieldIntrinsicEfficiency( static_cast<float>(ns_size_high_energy) );
    if( (example_low > 0.0) && (example_high > 0.0) && (drf_low > 0.0) && (drf_high > 0.0)
       && (stated_lower <= ns_size_low_energy)
       && ((stated_upper <= stated_lower) || (stated_upper >= ns_size_high_energy)) )
    {
      const double example_log = std::log( example_high/example_low );
      const double drf_log = std::log( drf_high/drf_low );
      if( example_log > ns_small_crystal_log_fep_ratio )
        large_weight = (drf_log - ns_small_crystal_log_fep_ratio) / (example_log - ns_small_crystal_log_fep_ratio);
      large_weight = std::max( 0.0, std::min( 1.0, large_weight ) );
    }
  }

  for( size_t i = 0; i < num_energies; ++i )
  {
    if( answer.fep[i] <= 0.0 )
    {
      answer.total[i] = 0.0;
      continue;
    }

    const double large = std::max( 1.0, ratio_fcn( std::max( energies_kev[i], example_lower_kev ) ) );
    const double small_ratio = small_crystal_total_to_peak( energies_kev[i] );
    const double ratio = std::exp( (1.0 - large_weight)*std::log( small_ratio ) + large_weight*std::log( large ) );
    answer.total[i] = std::min( 1.0, answer.fep[i]*ratio );
  }//for( loop over energies )

  return answer;
}//detector_efficiencies(...)


struct SummingModel::MaterialMuCache
{
  std::mutex mutex;
  // Keyed on the material object: a shared_ptr is kept so the address can not be reused by another material.
  map<shared_ptr<const Material>,shared_ptr<const vector<double>>> values;
};//struct MaterialMuCache


SummingModel::SummingModel( const std::vector<SourceInput> &sources,
                            const std::shared_ptr<const DetectorPeakResponse> &drf,
                            const bool shield_scatter )
  : m_drf( drf ),
    m_material_mu( make_shared<MaterialMuCache>() )
{
  // The detector's face geometry and validity (throws if unusable), before any real work.
  m_detector = detector_efficiencies( drf, {} );

  // A geometry-aware response knows its own near-field efficiency; take distances from it.
  if( drf->ceeloResponse() )
  {
    const double ref_energy = std::max( ns_solid_angle_ref_energy, drf->lowerEnergy() );
    const DetectorEfficiencies ref = detector_efficiencies( drf, { ref_energy } );
    if( (ref.total_source == TotalEffSource::DetectorResponse) && (ref.total[0] > 0.0) )
    {
      m_ref_energy = ref_energy;
      m_ref_total_intrinsic = ref.total[0];
    }
  }//if( geometry-aware response )

  m_max_solid_angle = axis_solid_angle( 0.0 );

  struct PendingPair
  {
    ceelo::CoincidenceSumPair pair;
    double sum = 0.0, max_rate = 0.0;
  };

  struct PendingSource
  {
    shared_ptr<const ceelo::AnalyticCascadeModel> model;
    vector<vector<PendingPair>> peaks;   //the pairs of each kept sum peak
    vector<double> peak_energy, peak_max_rate;
    vector<int> peak_line;
    vector<ceelo::SummingOutProgram::Pair> program_pairs;  //the distinct transition pairs of the kept pairs
  };

  vector<PendingSource> pending( sources.size() );
  m_sources.resize( sources.size() );

  for( size_t src = 0; src < sources.size(); ++src )
  {
    const SourceInput &input = sources[src];
    if( !input.nuclide )
      throw runtime_error( "RelActCalcCascade::SummingModel: null nuclide." );

    Source &source = m_sources[src];
    source.input = input;

    PendingSource &pend = pending[src];
    pend.model = cascade_model( input.nuclide, input.age );

    // The cascade model's rates are per decay of the parent at age zero; callers (like RelActCalcAuto's
    //  yields) work per decay of the parent now.
    const double per_current_decay = std::exp( std::min( 50.0, input.nuclide->decayConstant() * input.age ) );

    // Merge overlapping ranges, so a pair is not found twice.
    vector<SumRange> ranges = input.sum_ranges;
    std::sort( begin(ranges), end(ranges), []( const SumRange &lhs, const SumRange &rhs ){
      return lhs.lower_energy < rhs.lower_energy;
    } );
    vector<SumRange> merged;
    for( const SumRange &range : ranges )
    {
      if( !(range.upper_energy > range.lower_energy) )
        continue;
      if( !merged.empty() && (range.lower_energy <= merged.back().upper_energy) )
      {
        merged.back().upper_energy = std::max( merged.back().upper_energy, range.upper_energy );
        merged.back().min_rate_per_decay = std::min( merged.back().min_rate_per_decay, range.min_rate_per_decay );
      }else
      {
        merged.push_back( range );
      }
    }//for( loop over ranges )

    for( const SumRange &range : merged )
    {
      const vector<ceelo::CoincidenceSumPair> found
                    = ceelo::coincidence_sum_pairs( *pend.model, range.lower_energy, range.upper_energy );

      // Each pair's rate with the source touching the detector, nothing attenuated; a pair with a member the
      //  detector can not see contributes nothing.
      vector<double> pair_energies;
      for( const ceelo::CoincidenceSumPair &pair : found )
      {
        pair_energies.push_back( pair.energy_a );
        pair_energies.push_back( pair.energy_b );
      }
      sort_unique( pair_energies );
      const DetectorEfficiencies pair_effs = detector_efficiencies( drf, pair_energies );

      vector<PendingPair> candidates;
      for( const ceelo::CoincidenceSumPair &pair : found )
      {
        const double fep_a = pair_effs.fep[index_of_energy( pair_energies, pair.energy_a )];
        const double fep_b = pair_effs.fep[index_of_energy( pair_energies, pair.energy_b )];
        PendingPair cand;
        cand.pair = pair;
        cand.pair.weight *= per_current_decay;
        cand.pair.joint *= per_current_decay;
        cand.pair.emission_a *= per_current_decay;
        cand.pair.emission_b *= per_current_decay;
        cand.sum = pair.energy_a + pair.energy_b;
        cand.max_rate = cand.pair.weight * fep_a * fep_b * m_max_solid_angle * m_max_solid_angle;
        if( cand.max_rate > 0.0 )
          candidates.push_back( cand );
      }//for( loop over found pairs )

      std::sort( begin(candidates), end(candidates), []( const PendingPair &lhs, const PendingPair &rhs ){
        return lhs.sum < rhs.sum;
      } );

      for( size_t first = 0; first < candidates.size(); )
      {
        size_t last = first + 1;
        while( (last < candidates.size()) && ((candidates[last].sum - candidates[first].sum) <= ns_sum_merge_kev) )
          ++last;

        vector<PendingPair> group( begin(candidates) + first, begin(candidates) + last );
        first = last;

        auto rate_and_energy = []( const vector<PendingPair> &pairs, double &rate, double &energy ){
          rate = energy = 0.0;
          for( const PendingPair &p : pairs )
          {
            rate += p.max_rate;
            energy += p.max_rate * p.sum;
          }
          energy = (rate > 0.0) ? (energy / rate) : 0.0;
        };

        double rate = 0.0, energy = 0.0;
        rate_and_energy( group, rate, energy );

        // Summing-in to one of the source's lines: a pair with a member that is that line is not a sum into
        //  it (the engine's rule; such a pair is the line plus a soft photon).
        int line_index = -1;
        double nearest = SummingModel::sm_sum_on_line_kev;
        for( size_t line = 0; line < input.line_energies.size(); ++line )
        {
          const double diff = fabs( input.line_energies[line] - energy );
          if( diff <= nearest )
          {
            nearest = diff;
            line_index = static_cast<int>( line );
          }
        }

        if( line_index >= 0 )
        {
          const double line_energy = input.line_energies[line_index];
          group.erase( std::remove_if( begin(group), end(group), [line_energy]( const PendingPair &p ){
            return (fabs( p.pair.energy_a - line_energy ) <= ns_line_match_kev)
                   || (fabs( p.pair.energy_b - line_energy ) <= ns_line_match_kev);
          } ), end(group) );
          rate_and_energy( group, rate, energy );
        }//if( line_index >= 0 )

        if( group.empty() || !(rate > 0.0) || (rate < range.min_rate_per_decay) )
          continue;

        // Keep the strongest pairs; the peak's energy and rate bound stay those of all its pairs.
        std::sort( begin(group), end(group), []( const PendingPair &lhs, const PendingPair &rhs ){
          return lhs.max_rate > rhs.max_rate;
        } );
        double kept_rate = 0.0;
        size_t num_keep = 0;
        while( (num_keep < group.size()) && (kept_rate < (1.0 - ns_sum_peak_dropped_fraction)*rate) )
          kept_rate += group[num_keep++].max_rate;
        group.resize( num_keep );

        // The kept pairs stand in for the dropped ones, so the peak is not biased low by what was dropped.
        if( kept_rate > 0.0 )
        {
          for( PendingPair &p : group )
            p.pair.weight *= (rate / kept_rate);
        }

        pend.peaks.push_back( std::move(group) );
        pend.peak_energy.push_back( energy );
        pend.peak_max_rate.push_back( rate );
        pend.peak_line.push_back( line_index );
      }//for( loop over groups of pairs with (nearly) the same sum )
    }//for( loop over sum ranges )

    // The summing-out program: the modeled lines, then the gamma members of the sum pairs whose coincident
    //  survival is estimated from their members' own summing-out (those not of two level-scheme transitions).
    vector<ceelo::SummingOutProgram::Line> program_lines;
    for( const double energy : input.line_energies )
      program_lines.push_back( { energy, ns_line_match_kev } );

    for( const vector<PendingPair> &peak : pend.peaks )
    {
      for( const PendingPair &p : peak )
      {
        if( (p.pair.transition_a >= 0) && (p.pair.transition_b >= 0) )
          continue;  //gets the engine's exact survival

        const pair<double,bool> members[2] = { { p.pair.energy_a, p.pair.a_is_gamma },
                                               { p.pair.energy_b, p.pair.b_is_gamma } };
        for( const pair<double,bool> &member : members )
        {
          if( !member.second )
            continue;
          const bool have = std::any_of( begin(program_lines), end(program_lines),
            [&member]( const ceelo::SummingOutProgram::Line &line ){
              return fabs( line.energy_keV - member.first ) <= ns_line_match_kev;
            } );
          if( !have )
            program_lines.push_back( { member.first, ns_line_match_kev } );
        }
      }//for( loop over pairs )
    }//for( loop over sum peaks )

    // The pairs of two level-scheme transitions get the engine's exact coincident survival, which depends only
    //  on the two transitions: a gamma with each x-ray line of a converting transition shares one.
    vector<ceelo::SummingOutProgram::Pair> &program_pairs = pend.program_pairs;
    for( const vector<PendingPair> &peak : pend.peaks )
    {
      for( const PendingPair &p : peak )
      {
        if( (p.pair.transition_a < 0) || (p.pair.transition_b < 0) )
          continue;
        const bool have = std::any_of( begin(program_pairs), end(program_pairs),
          [&p]( const ceelo::SummingOutProgram::Pair &pp ){
            return (pp.branch == p.pair.branch) && (pp.transition_a == p.pair.transition_a)
                   && (pp.transition_b == p.pair.transition_b);
          } );
        if( !have )
          program_pairs.push_back( { p.pair.branch, p.pair.transition_a, p.pair.transition_b } );
      }
    }

    source.program = make_shared<const ceelo::SummingOutProgram>( pend.model, program_lines, program_pairs );

    const size_t num_lines = input.line_energies.size();
    source.line_found.assign( num_lines, 0 );
    source.line_complete.assign( num_lines, 1 );
    source.line_emission.assign( num_lines, 0.0 );
    for( size_t line = 0; line < num_lines; ++line )
    {
      source.line_emission[line] = per_current_decay * source.program->line_emission( line );
      source.line_found[line] = (source.line_emission[line] > 0.0);
      source.line_complete[line] = source.program->line_complete( line );
    }
  }//for( loop over sources )

  // Which of each program's energies can matter: the first-order loss a photon of that energy causes any of
  //  the program's lines, at the largest total efficiency it could have.
  vector<vector<char>> energy_matters( m_sources.size() );
  vector<double> all_energies;
  for( size_t src = 0; src < m_sources.size(); ++src )
  {
    const ceelo::SummingOutProgram &program = *m_sources[src].program;
    const vector<double> &prog_energies = program.energies();
    const DetectorEfficiencies effs = detector_efficiencies( drf, prog_energies );
    const vector<vector<double>> coefs = program.first_order_coefficients();

    energy_matters[src].assign( prog_energies.size(), 0 );
    for( size_t i = 0; i < prog_energies.size(); ++i )
    {
      double max_coef = 0.0;
      for( const vector<double> &line_coefs : coefs )
        max_coef = std::max( max_coef, fabs( line_coefs[i] ) );
      if( (max_coef * effs.total[i] * m_max_solid_angle) >= ns_min_partner_effect )
      {
        energy_matters[src][i] = 1;
        all_energies.push_back( prog_energies[i] );
      }
    }//for( loop over the program's energies )

    for( size_t line = 0; line < program.num_lines(); ++line )
      all_energies.push_back( program.line( line ).energy_keV );

    for( const vector<PendingPair> &peak : pending[src].peaks )
    {
      for( const PendingPair &p : peak )
      {
        all_energies.push_back( p.pair.energy_a );
        all_energies.push_back( p.pair.energy_b );
      }
    }
  }//for( loop over sources )

  sort_unique( all_energies );
  m_energies = all_energies;
  m_detector = detector_efficiencies( drf, m_energies );
  build_near_field_tables();
  if( shield_scatter )
    build_scatter_table();

  for( size_t src = 0; src < m_sources.size(); ++src )
  {
    Source &source = m_sources[src];
    const ceelo::SummingOutProgram &program = *source.program;
    const vector<double> &prog_energies = program.energies();

    source.program_energy.assign( prog_energies.size(), -1 );
    for( size_t i = 0; i < prog_energies.size(); ++i )
    {
      if( energy_matters[src][i] )
        source.program_energy[i] = static_cast<int>( index_of_energy( m_energies, prog_energies[i] ) );
    }

    source.program_line_energy.resize( program.num_lines() );
    for( size_t line = 0; line < program.num_lines(); ++line )
      source.program_line_energy[line] = index_of_energy( m_energies, program.line( line ).energy_keV );

    auto program_line_of = [&program]( const double energy, const bool is_gamma ) -> int {
      if( !is_gamma )
        return -1;
      for( size_t line = 0; line < program.num_lines(); ++line )
      {
        if( fabs( program.line( line ).energy_keV - energy ) <= ns_line_match_kev )
          return (program.line_emission( line ) > 0.0) ? static_cast<int>( line ) : -1;
      }
      return -1;
    };//program_line_of

    const PendingSource &pend = pending[src];
    for( size_t p = 0; p < pend.peaks.size(); ++p )
    {
      SumPeak peak;
      peak.energy = pend.peak_energy[p];
      peak.line_index = pend.peak_line[p];
      peak.max_rate_per_decay = pend.peak_max_rate[p];
      for( const PendingPair &cand : pend.peaks[p] )
      {
        const ceelo::CoincidenceSumPair &cp = cand.pair;
        SumPair pair;
        pair.energy_a = index_of_energy( m_energies, cp.energy_a );
        pair.energy_b = index_of_energy( m_energies, cp.energy_b );
        pair.weight = cp.weight;
        // `weight` carries the angular factor of a both-full-energy pair; a partner's any-deposit does not.
        //  These are expected numbers of the partner per photon, so can pass one: an annihilation gives two.
        pair.prob_b_given_a = (cp.emission_a > 0.0) ? std::min( 2.0, cp.joint / cp.emission_a ) : 0.0;
        pair.prob_a_given_b = (cp.emission_b > 0.0) ? std::min( 2.0, cp.joint / cp.emission_b ) : 0.0;
        for( size_t i = 0; (pair.program_pair < 0) && (i < pend.program_pairs.size()); ++i )
        {
          const ceelo::SummingOutProgram::Pair &pp = pend.program_pairs[i];
          if( (pp.branch == cp.branch) && (pp.transition_a == cp.transition_a) && (pp.transition_b == cp.transition_b) )
            pair.program_pair = static_cast<int>( i );
        }
        assert( (pair.program_pair >= 0) == ((cp.transition_a >= 0) && (cp.transition_b >= 0)) );
        pair.program_line_a = program_line_of( cp.energy_a, cp.a_is_gamma );
        pair.program_line_b = program_line_of( cp.energy_b, cp.b_is_gamma );
        peak.pairs.push_back( pair );
      }//for( loop over the peak's pairs )

      source.sum_peaks.push_back( std::move(peak) );
    }//for( loop over sum peaks )

    std::sort( begin(source.sum_peaks), end(source.sum_peaks), []( const SumPeak &lhs, const SumPeak &rhs ){
      return lhs.energy < rhs.energy;
    } );
  }//for( loop over sources )

  // Whether there is anything for the geometry to change: a sum peak, or a line that loses counts with the
  //  source touching the detector and nothing attenuated.
  EvalInput<double> closest;
  closest.front_solid_angle = m_max_solid_angle;
  closest.external_transmission.assign( m_energies.size(), 1.0 );
  EvalResult<double> at_closest;
  evaluate( closest, at_closest );
  for( size_t src = 0; src < m_sources.size(); ++src )
  {
    m_has_effect = (m_has_effect || !m_sources[src].sum_peaks.empty());
    for( const double factor : at_closest.summing_out[src] )
      m_has_effect = (m_has_effect || (fabs( factor - 1.0 ) > 1.0E-9));
  }
}//SummingModel constructor


/** A geometry-aware response's on-axis near-field efficiencies, as ratios to "front solid angle times
 far-field intrinsic efficiency": [energy][solid angle], not-a-number where the response could not say.
 */
struct SummingModel::NearFieldGrid
{
  vector<double> energies;      //keV, ascending
  vector<double> solid_angle;   //ascending from zero (far away, where the ratios are one)
  vector<vector<double>> fep, total;
  /** Whether the response has a full-energy efficiency at each of `energies` (not below its range). */
  vector<char> detected;
};//struct NearFieldGrid


std::shared_ptr<const SummingModel::NearFieldGrid> SummingModel::near_field_grid() const
{
  typedef shared_ptr<const NearFieldGrid> GridPtr;
  static std::mutex s_mutex;
  static map<uint64_t,GridPtr> s_grids;
  const size_t max_cached = 8;

  const uint64_t key = m_drf->hashValue();
  {
    std::lock_guard<std::mutex> lock( s_mutex );
    const map<uint64_t,GridPtr>::const_iterator pos = s_grids.find( key );
    if( pos != end(s_grids) )
      return pos->second;
  }

  const shared_ptr<NearFieldGrid> answer = make_shared<NearFieldGrid>();
  answer->energies = response_grid_energies();
  const size_t num_energies = answer->energies.size();
  const DetectorEfficiencies detector = detector_efficiencies( m_drf, answer->energies );

  answer->solid_angle.assign( 1, 0.0 );
  answer->fep.assign( num_energies, vector<double>{ 1.0 } );
  answer->total.assign( num_energies, vector<double>{ 1.0 } );
  answer->detected.resize( num_energies );
  for( size_t i = 0; i < num_energies; ++i )
    answer->detected[i] = (detector.fep[i] > 0.0);

  // From far to touching; the front solid angle rises all the way.
  const double distances_cm[] = { 200.0, 80.0, 40.0, 25.0, 16.0, 11.0, 8.0, 6.0, 4.5, 3.25, 2.25, 1.5, 1.0,
                                  0.5, 0.2, 0.0 };
  for( const double distance_cm : distances_cm )
  {
    const double distance = distance_cm * PhysicalUnits::cm;
    const double solid_angle = axis_solid_angle( distance );
    if( !(solid_angle > (1.0001*answer->solid_angle.back() + 1.0E-9)) )
      continue;   //the cap at touching, or a response that is flat here

    const std::shared_ptr<const DetectorPeakResponse::PositionedQuadrature> quadrature
                                                        = m_drf->apertureQuadrature( 0.0, 0.0, distance );
    answer->solid_angle.push_back( solid_angle );
    for( size_t i = 0; i < num_energies; ++i )
    {
      double fep_ratio = 1.0, tot_ratio = 1.0;
      try
      {
        const float energy = static_cast<float>( answer->energies[i] );
        if( detector.fep[i] > 0.0 )
        {
          const DetectorPeakResponse::EffEval eval = m_drf->fepEfficiencyEval( energy, 0.0, 0.0, distance, quadrature );
          fep_ratio = eval.value / (solid_angle * detector.fep[i]);
        }
        // A total that is not the response's own has no near-field behaviour of the response's to follow.
        if( (detector.total[i] > 0.0) && (detector.total_source == TotalEffSource::DetectorResponse) )
        {
          const DetectorPeakResponse::EffEval eval = m_drf->totalEfficiencyEval( energy, 0.0, 0.0, distance );
          tot_ratio = eval.value / (solid_angle * detector.total[i]);
        }
      }catch( std::exception & )
      {
        fep_ratio = tot_ratio = std::numeric_limits<double>::quiet_NaN();
      }

      answer->fep[i].push_back( fep_ratio );
      answer->total[i].push_back( tot_ratio );
    }//for( loop over energies )
  }//for( loop over distances )

  std::lock_guard<std::mutex> lock( s_mutex );
  if( s_grids.size() >= max_cached )
    s_grids.clear();
  s_grids[key] = answer;
  return answer;
}//near_field_grid()


void SummingModel::build_near_field_tables()
{
  m_near_solid_angle.clear();
  m_near_fep_ratio.clear();
  m_near_fep_slope.clear();
  m_near_tot_ratio.clear();
  m_near_tot_slope.clear();

  if( !(m_ref_total_intrinsic > 0.0) || m_energies.empty() )
    return;   //not a geometry-aware response, or nothing to evaluate

  const shared_ptr<const NearFieldGrid> response = near_field_grid();
  const vector<double> &grid = response->solid_angle;
  if( grid.size() < 3 )
    return;

  // Each energy's ratios, between those of the response's grid energies on either side (in the log of the
  //  energy); from one side only where the response detects nothing at the other.
  const size_t num_energies = m_energies.size();
  const size_t num_grid_energies = response->energies.size();
  vector<vector<double>> fep( num_energies, vector<double>( grid.size(), 1.0 ) ), tot = fep;
  for( size_t i = 0; i < num_energies; ++i )
  {
    if( !(m_detector.fep[i] > 0.0) )
      continue;   //not detected: nothing to scale

    const double energy = m_energies[i];
    size_t upper = static_cast<size_t>( std::lower_bound( begin(response->energies), end(response->energies), energy )
                                        - begin(response->energies) );
    upper = std::min( upper, num_grid_energies - 1 );
    size_t lower = (upper > 0) ? (upper - 1) : upper;
    if( energy >= response->energies[upper] )
      lower = upper;
    if( !response->detected[lower] )
      lower = upper;
    if( !response->detected[upper] )
      upper = lower;

    const double frac = (upper == lower) ? 0.0
                        : (std::log( energy/response->energies[lower] )
                           / std::log( response->energies[upper]/response->energies[lower] ));
    for( size_t k = 0; k < grid.size(); ++k )
    {
      fep[i][k] = (1.0 - frac)*response->fep[lower][k] + frac*response->fep[upper][k];
      tot[i][k] = (1.0 - frac)*response->total[lower][k] + frac*response->total[upper][k];
    }
  }//for( loop over energies )

  // An energy the response can not say for (a ratio that is not a number, or nonsense) keeps the far-field
  //  scaling at every distance - resetting single grid points would put steps into the table.
  const auto sane_or_flat = []( vector<double> &row ){
    const bool sane = std::all_of( begin(row), end(row), []( const double ratio ) -> bool {
      return std::isfinite( ratio ) && (ratio >= 0.2) && (ratio <= 5.0);
    } );
    if( !sane )
      std::fill( begin(row), end(row), 1.0 );
  };
  for( size_t i = 0; i < num_energies; ++i )
  {
    sane_or_flat( fep[i] );
    sane_or_flat( tot[i] );
  }

  // Slopes of a monotone cubic interpolant (Fritsch-Carlson): no overshoot between grid points.
  auto slopes_of = [&grid]( const vector<double> &y ) -> vector<double> {
    const size_t n = grid.size();
    vector<double> secant( n - 1 ), slope( n, 0.0 );
    for( size_t k = 0; (k + 1) < n; ++k )
      secant[k] = (y[k + 1] - y[k]) / (grid[k + 1] - grid[k]);
    slope[0] = secant[0];
    slope[n - 1] = secant[n - 2];
    for( size_t k = 1; (k + 1) < n; ++k )
    {
      if( (secant[k - 1]*secant[k]) <= 0.0 )
      {
        slope[k] = 0.0;
      }else
      {
        const double w1 = 2.0*(grid[k + 1] - grid[k]) + (grid[k] - grid[k - 1]);
        const double w2 = (grid[k + 1] - grid[k]) + 2.0*(grid[k] - grid[k - 1]);
        slope[k] = (w1 + w2) / (w1/secant[k - 1] + w2/secant[k]);
      }
    }
    return slope;
  };//slopes_of

  m_near_solid_angle = grid;
  m_near_fep_ratio = fep;
  m_near_tot_ratio = tot;
  m_near_fep_slope.resize( num_energies );
  m_near_tot_slope.resize( num_energies );
  for( size_t i = 0; i < num_energies; ++i )
  {
    m_near_fep_slope[i] = slopes_of( fep[i] );
    m_near_tot_slope[i] = slopes_of( tot[i] );
  }
}//build_near_field_tables()


void SummingModel::build_scatter_table()
{
  m_scatter.reset();
  if( m_energies.empty() )
    return;

  // Interpolated from the response's table on a fixed grid of energies, which is built once: building one
  //  for a model's own energies takes about a second, and a fit makes several models.
  const shared_ptr<const GammaInteractionCalc::ShieldScatterAugment> grid = scatter_grid( m_drf );
  if( !grid )
    return;   //no scatter data: `models_shield_scatter()` says so

  const shared_ptr<GammaInteractionCalc::ShieldScatterAugment> scatter
                                          = make_shared<GammaInteractionCalc::ShieldScatterAugment>();
  scatter->buildFrom( *grid, m_energies );
  assert( scatter->valid() && (scatter->energies() == m_energies) );
  if( scatter->valid() && (scatter->energies() == m_energies) )
    m_scatter = scatter;
}//build_scatter_table()


bool SummingModel::models_shield_scatter() const
{
  return (m_scatter && m_scatter->valid());
}


double SummingModel::axis_solid_angle( const double distance ) const
{
  const double dist = std::max( 0.0, distance );
  if( m_ref_total_intrinsic > 0.0 )
  {
    try
    {
      const DetectorPeakResponse::EffEval eval
                    = m_drf->totalEfficiencyEval( static_cast<float>(m_ref_energy), 0.0, 0.0, dist );
      if( std::isfinite( eval.value ) && (eval.value > 0.0) )
        return std::min( 1.0, eval.value / m_ref_total_intrinsic );
    }catch( std::exception & )
    {
      // Fall through to the flat disk.
    }
  }//if( distances are taken from the response )

  return DetectorPeakResponse::fractionalSolidAngle( m_detector.detector_diameter, dist + m_detector.detector_setback );
}//axis_solid_angle(...)


double SummingModel::front_solid_angle( const double distance, const double source_radius ) const
{
  const double axis = axis_solid_angle( distance );
  if( !(source_radius > 0.0) )
    return axis;

  // How much less of the crystal's disk the rest of the source's face sees than its centre does: <g^2>/<g>
  //  over equal-area rings of the face, relative to the on-axis value.
  const double crystal_radius = 0.5*m_detector.crystal_diameter;
  const double height = std::max( 0.0, distance ) + m_detector.crystal_setback;

  const int num_rings = 16;
  double sum = 0.0, sum_sq = 0.0;
  for( int ring = 0; ring < num_rings; ++ring )
  {
    const double offset = source_radius*std::sqrt( (ring + 0.5)/num_rings );
    const double g = off_axis_solid_angle( offset, height, crystal_radius );
    sum += g;
    sum_sq += g*g;
  }
  const double on_axis = off_axis_solid_angle( 0.0, height, crystal_radius );
  if( !(sum > 0.0) || !(on_axis > 0.0) )
    return axis;

  return axis * (sum_sq/sum) / on_axis;
}//front_solid_angle(...)


double SummingModel::front_distance( const double solid_angle, const double source_radius ) const
{
  // The solid angle falls monotonically with distance: bisect.
  double near_distance = 0.0, far_distance = 10.0*PhysicalUnits::meter;
  if( solid_angle >= front_solid_angle( near_distance, source_radius ) )
    return near_distance;
  if( solid_angle <= front_solid_angle( far_distance, source_radius ) )
    return far_distance;

  for( int iter = 0; iter < 60; ++iter )
  {
    const double mid = 0.5*(near_distance + far_distance);
    if( front_solid_angle( mid, source_radius ) > solid_angle )
      near_distance = mid;
    else
      far_distance = mid;
    if( (far_distance - near_distance) < 1.0E-4*PhysicalUnits::cm )
      break;
  }
  return 0.5*(near_distance + far_distance);
}//front_distance(...)


std::shared_ptr<const std::vector<double>> SummingModel::material_mu( const std::shared_ptr<const Material> &material ) const
{
  if( !material )
    throw runtime_error( "SummingModel::material_mu: null material." );

  std::lock_guard<std::mutex> lock( m_material_mu->mutex );
  const map<shared_ptr<const Material>,shared_ptr<const vector<double>>>::const_iterator pos
                                                                    = m_material_mu->values.find( material );
  if( pos != end(m_material_mu->values) )
    return pos->second;

  const shared_ptr<vector<double>> mu = make_shared<vector<double>>( m_energies.size(), 0.0 );
  for( size_t i = 0; i < m_energies.size(); ++i )
    (*mu)[i] = GammaInteractionCalc::transmition_length_coefficient( material.get(), m_energies[i] ) / material->density;

  // A model can outlive many sets of options, each with its own material objects: do not keep them all.
  const size_t max_materials = 32;
  if( m_material_mu->values.size() >= max_materials )
    m_material_mu->values.clear();

  m_material_mu->values[material] = mu;
  return mu;
}//material_mu(...)

}//namespace RelActCalcCascade
