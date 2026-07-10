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

/* GADRAS-scatter effect on cascade corrections (quantification, per plan).
 Net summing factor C_net with vs without the ShieldScatterAugment term, for
 point sources at 10 cm from the Detective-X model grounded to the LANL curve
 (test_CascadeSummingFit's `CascadeScatterQuantification` case, analytic + DRF
 curves, no MC):

   Co-60  1173 keV behind Fe 10 mm:  no-scatter 0.9888  with 0.9831  (-0.57%)
   Co-60  1173 keV behind Pb  6 mm:  no-scatter 0.9884  with 0.9858  (-0.26%)
   Ba-133  356 keV behind Fe 10 mm:  ~1.0000          ~1.0000       (<0.01%)
   Ba-133  356 keV behind Pb  6 mm:  ~1.0000          ~1.0000       (<0.01%)

 The scatter continuum keeps the coincident partner's TOTAL efficiency higher
 than pure transmission predicts, so it slightly deepens summing-OUT (C_net a
 few tenths of a percent lower for Co-60).  For Ba-133 behind shielding the
 low-energy coincidence partners are removed by differential transmission, so
 there is little summing left and the scatter term is negligible.  In all cases
 the effect is <~0.6% on C_net: the shielding effect on summing is
 transmission-dominated and the scatter augmentation is a second-order restore
 (see mc_det_eff_plan/studies/cascade/cascade_shield_gadras_benchmark.md,
 Result 2).  The full shielded correction (transmission + this augmentation) is
 validated end-to-end by the shielded-point truth gates in
 test_CascadeSummingFit (Co-60/Ba-133 behind Al/Fe/Pb pass).
 */

#include "InterSpec_config.h"

#include <set>
#include <cmath>
#include <mutex>
#include <string>
#include <memory>
#include <vector>
#include <stdexcept>

#include "SandiaDecay.h"

#include "SpecUtils/Filesystem.h"

#include "cascade/CascadeTypes.h"
#include "cascade/AnalyticCascade.h"
#include "cascade/SandiaDecayCascade.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/CascadeSummingCalc.h"
#include "InterSpec/GadrasShieldScatter.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace GammaInteractionCalc
{

// ============================ ShieldScatterAugment ==========================

ShieldScatterAugment::ShieldScatterAugment()
{
}


void ShieldScatterAugment::build( const std::vector<double> &energies,
                                  const std::function<double(double)> &eps_totint,
                                  const double fracH,
                                  const std::string &data_directory )
{
  m_energies.clear();
  m_table.clear();

  if( energies.empty() )
    return;

  const string db_path = SpecUtils::append_path( data_directory,
                                                 "sandia.shieldscatter.db" );
  const GadrasShieldScatter scatter( db_path );  //throws on failure

  // Grid axes.  AD includes a zero anchor (no shield => no scatter); AN spans
  //  the tabulated range and clamps outside (see evaluate()).
  m_an_grid = { 4.0, 13.0, 26.0, 42.0, 60.0, 82.0,
                scatter.maxAtomicNumber() };
  m_ad_grid = { 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
                scatter.maxArealDensity() };

  m_energies = energies;
  std::sort( begin(m_energies), end(m_energies) );

  const double frac_h = std::max( 0.0, std::min( 1.0, fracH ) );

  m_table.assign( m_energies.size() * m_an_grid.size() * m_ad_grid.size(), 0.0 );

  vector<double> bounds, scatter_per_leak;
  for( size_t ei = 0; ei < m_energies.size(); ++ei )
  {
    const double energy = m_energies[ei];
    if( energy < 20.0 )
      continue;  //nothing meaningfully scatters through a shield at these energies

    const double eps_primary = eps_totint( energy );
    if( eps_primary <= 1.0E-8 )
      continue;

    scatter.groupBounds( energy, bounds );

    // Intrinsic total efficiency at the scatter-group centers; contributions
    //  below 15 keV are dropped (they do not reach the crystal, and the DRF
    //  curve is typically not defined there).
    vector<double> eps_at_center( scatter.groupCount(), 0.0 );
    for( int k = 0; k < scatter.groupCount(); ++k )
    {
      const double center = 0.5*(bounds[k] + bounds[k+1]);
      if( center >= 15.0 )
        eps_at_center[k] = eps_totint( center );
    }

    for( size_t ai = 0; ai < m_an_grid.size(); ++ai )
    {
      for( size_t di = 0; di < m_ad_grid.size(); ++di )
      {
        const double ad = m_ad_grid[di];
        if( ad <= 0.0 )
          continue;  //A(AD = 0) = 0 anchor

        // Per-uncollided-LEAKAGE-photon scatter groups; the augment term in
        //  eps_tot_shielded = eps_tot_bare*( T + A ) needs per-SOURCE-photon
        //  scatter, i.e. T * S_leak - so scale by GADRAS's own transmission.
        scatter.computeShieldScatter( energy, m_an_grid[ai], ad, frac_h,
                                      0.0 /*point-in-sphere*/, scatter_per_leak );

        vector<float> answer;
        vector<float> binning( bounds.size() - 1 );
        for( size_t k = 0; (k + 1) < bounds.size(); ++k )
          binning[k] = static_cast<float>( bounds[k] );
        const float transmitted = scatter.getContinuum( answer, static_cast<float>(energy),
                                        1.0f, static_cast<float>(m_an_grid[ai]),
                                        static_cast<float>(ad),
                                        static_cast<float>(frac_h), binning, 0.0f );
        (void)transmitted;  //getContinuum already scales `answer` per source photon

        double aug = 0.0;
        for( size_t k = 0; (k < answer.size()) && (k < eps_at_center.size()); ++k )
          aug += static_cast<double>(answer[k]) * eps_at_center[k] / eps_primary;

        m_table[ (ei*m_an_grid.size() + ai)*m_ad_grid.size() + di ] = aug;
      }//for( AD grid )
    }//for( AN grid )
  }//for( energies )
}//ShieldScatterAugment::build(...)


// ============================ CascadeSummingCalc ============================

namespace
{
  /** Age bucket: ~2% steps in log-age, so fitting the age does not force a
   cascade re-enumeration every iteration.
   */
  int age_bucket( const double age )
  {
    const double age_s = age / PhysicalUnits::second;
    if( age_s <= 1.0 )
      return -1;
    return static_cast<int>( std::floor( 115.0 * std::log10(age_s) ) );
  }
}//namespace


CascadeSummingCalc::CascadeSummingCalc(
              const std::vector<std::pair<const SandiaDecay::Nuclide *,double>> &nuclide_initial_ages,
              const std::vector<std::pair<double,double>> &peak_energy_widths,
              const double photopeak_cluster_sigma,
              const std::shared_ptr<const DetectorPeakResponse> &drf,
              const double shield_fracH,
              const bool shields_present,
              const std::string &data_directory )
{
  if( !drfHasNeededInfo(drf) )
    throw runtime_error( "Cascade summing correction requires the detector"
        " response to have total-efficiency information; add it via the"
        " detector editor (Detector Response Select -> Modify), e.g. by"
        " attaching a Monte-Carlo characterization." );

  // The fit's peak windows: clustering matches lines within
  //  photopeak_cluster_sigma * sigma of the peak's gamma energy.
  for( const pair<double,double> &ew : peak_energy_widths )
  {
    ceelo::PeakWindow win;
    win.energy_keV = ew.first;
    win.tolerance_keV = std::max( 0.25, photopeak_cluster_sigma * ew.second );
    m_windows.push_back( win );
  }

  // Pre-enumerate cascades at the starting ages, and collect every emission
  //  energy an evaluation can query (gammas, vacancy K/L x-ray lines, 511).
  set<double> energies;
  set<int> daughter_zs;
  for( const pair<const SandiaDecay::Nuclide *,double> &na : nuclide_initial_ages )
  {
    const shared_ptr<const vector<ceelo::DecayCascade>> casc
                                            = cascades( na.first, na.second );
    for( const ceelo::DecayCascade &dc : *casc )
    {
      for( const ceelo::CascadeMember &m : dc.members )
      {
        if( m.energy_keV >= 5.0 )
          energies.insert( m.energy_keV );
      }
      const int z = dc.daughter_Z ? dc.daughter_Z : dc.level_scheme.daughter_Z;
      if( z > 0 )
        daughter_zs.insert( z );
    }//for( branches )
  }//for( nuclides )

  for( const int z : daughter_zs )
  {
    for( const ceelo::detail::XrayLine &ln : ceelo::detail::k_lines(z) )
      if( ln.energy >= 5.0 )
        energies.insert( ln.energy );
    for( int s = 0; s < 3; ++s )
      for( const ceelo::detail::XrayLine &ln : ceelo::detail::l_lines(z, s) )
        if( ln.energy >= 5.0 )
          energies.insert( ln.energy );
  }//for( daughter Zs )

  energies.insert( 510.998950 );  //annihilation

  m_partner_energies.assign( begin(energies), end(energies) );

  if( shields_present )
  {
    const std::function<double(double)> eps_totint = [drf]( double energy ) -> double {
      try
      {
        return drf->totalEfficiencyEval( static_cast<float>(energy), 0.0, 0.0,
                                         1.0E6*PhysicalUnits::cm ).value
               / DetectorPeakResponse::fractionalSolidAngle(
                          drf->detectorDiameter(),
                          1.0E6*PhysicalUnits::cm + drf->detectorSetback() );
      }catch( std::exception & )
      {
        return 0.0;
      }
    };

    m_scatter.build( m_partner_energies, eps_totint, shield_fracH, data_directory );

#if( PERFORM_DEVELOPER_CHECKS )
    // The scatter term must demonstrably engage for a real shield.
    if( m_scatter.valid() && !m_partner_energies.empty() )
    {
      double test_e = 661.0;
      for( const double e : m_partner_energies )
        if( (e > 300.0) && (e < 1500.0) ){ test_e = e; break; }
      const double a = m_scatter.evaluate( test_e, 26.0, 8.0 );
      assert( (a > 0.0) || (eps_totint(test_e) <= 1.0E-8) );
    }
#endif
  }//if( shields_present )
}//CascadeSummingCalc constructor


bool CascadeSummingCalc::drfHasNeededInfo( const std::shared_ptr<const DetectorPeakResponse> &drf )
{
  return drf && drf->isValid() && drf->hasAnyTotalEfficiencyInfo();
}


double CascadeSummingCalc::estimateMaxSummingMagnitude(
                              const std::shared_ptr<const DetectorPeakResponse> &drf,
                              const double distance )
{
  if( !drfHasNeededInfo(drf) )
    return 0.0;

  double max_tot = 0.0;
  for( const double energy : { 150.0, 300.0, 600.0, 1000.0, 1500.0 } )
  {
    try
    {
      if( drf->isFixedGeometry() )
      {
        max_tot = std::max( max_tot,
            static_cast<double>( drf->totalIntrinsicEfficiency( static_cast<float>(energy) ) ) );
      }else
      {
        const DetectorPeakResponse::EffEval eval
              = drf->totalEfficiencyEval( static_cast<float>(energy), 0.0, 0.0, distance );
        max_tot = std::max( max_tot, eval.value );
      }
    }catch( std::exception & )
    {
    }
  }//for( sample energies )

  return max_tot;
}//estimateMaxSummingMagnitude(...)


double CascadeSummingCalc::detectorFepEffAbs( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                              const double energy, const double theta,
                                              const double phi, const double distance )
{
  if( !drf || !drf->isValid() )
    return 0.0;

  try
  {
    if( drf->isFixedGeometry() )
      return drf->intrinsicEfficiency( static_cast<float>(energy) );

    const DetectorPeakResponse::EffEval eval
          = drf->fepEfficiencyEval( static_cast<float>(energy), theta, phi, distance );
    return std::max( 0.0, eval.value );
  }catch( std::exception & )
  {
    return 0.0;
  }
}//detectorFepEffAbs(...)


double CascadeSummingCalc::detectorTotEffAbs( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                              const double energy, const double theta,
                                              const double phi, const double distance )
{
  if( !drf || !drf->isValid() )
    return 0.0;

  try
  {
    if( drf->isFixedGeometry() )
      return drf->totalIntrinsicEfficiency( static_cast<float>(energy) );

    const DetectorPeakResponse::EffEval eval
          = drf->totalEfficiencyEval( static_cast<float>(energy), theta, phi, distance );
    return std::max( 0.0, eval.value );
  }catch( std::exception & )
  {
    return 0.0;
  }
}//detectorTotEffAbs(...)


std::shared_ptr<const std::vector<ceelo::DecayCascade>> CascadeSummingCalc::cascades(
                            const SandiaDecay::Nuclide *nuc, const double age ) const
{
  if( !nuc )
    throw runtime_error( "CascadeSummingCalc::cascades: null nuclide" );

  const pair<const SandiaDecay::Nuclide *,int> key( nuc, age_bucket(age) );

  {
    std::lock_guard<std::mutex> lock( m_cascade_mutex );
    const auto pos = m_cascades.find( key );
    if( pos != end(m_cascades) )
      return pos->second;
  }

  ceelo::cascade_adapter::CascadeOptions opts;
  opts.age_seconds = age / PhysicalUnits::second;
  //defaults: prompt_equilibrium, include_xrays, vacancy_xray_model,
  // include_annihilation, angular_correlations all on.

  auto casc = make_shared<vector<ceelo::DecayCascade>>(
                          ceelo::cascade_adapter::build_cascades( nuc, opts ) );

  std::lock_guard<std::mutex> lock( m_cascade_mutex );
  m_cascades[key] = casc;
  return casc;
}//cascades(...)

}//namespace GammaInteractionCalc
