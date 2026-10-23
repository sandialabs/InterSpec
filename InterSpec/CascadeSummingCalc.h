#ifndef CascadeSummingCalc_h
#define CascadeSummingCalc_h
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
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include "cascade/AnalyticCascade.h"

class DetectorPeakResponse;
class GadrasShieldScatter;
namespace SandiaDecay{ struct Nuclide; }

namespace GammaInteractionCalc
{

/** GADRAS shield-scatter augmentation A(E; AN, AD) of the total efficiency.

 For a shielded source, the photons that scatter in the shield but still reach
 the detector keep the total (any-deposit) efficiency - and hence cascade
 summing-OUT - higher than pure transmission predicts.  Per the validated
 recipe (mc_det_eff_plan/studies/cascade/cascade_shield_gadras_benchmark.md):

     eps_fep_shielded(E) = eps_fep_bare(E) * T(E)                 (no scatter term)
     eps_tot_shielded(E) = eps_tot_bare(E) * ( T(E) + A(E) )
     A(E)  =  Sum_k S(E'_k) * eps_totint(E'_k) / eps_totint(E)

 where S(E'_k) is the GADRAS scattered-photon spectrum at the detector per
 source photon (GadrasShieldScatter::getContinuum with unit intensity - its
 output already includes the transmission scaling) and eps_totint the
 detector's intrinsic total efficiency.  The scatter part of the recipe is
 ~20%-accurate; the summing change is transmission-dominated, so A is a
 second-order restore (see the benchmark's Result 2).

 The table is prebuilt (session thread, at chi2-function creation) on an
 (atomic-number, areal-density) grid for every emission energy the cascade
 evaluation will query; evaluation is bilinear in (AN, AD), so an automatic-
 differentiation dual number (ceres::Jet areal density along the fit ray)
 differentiates through it.  Hydrogen mass-fraction is fixed at the shield
 stack's construction-time value (second-order).
 */
class ShieldScatterAugment
{
public:
  /** Empty/no-op table: evaluate() returns 0. */
  ShieldScatterAugment();

  /** Builds the grid.

   @param energies  Every energy (keV) evaluate() may be asked for.
   @param eps_totint  The detector's intrinsic total efficiency vs energy (keV).
   @param fracH  Hydrogen mass fraction of the shield stack (0..1).
   @param data_directory  Directory holding `sandia.shieldscatter.db`.

   Throws std::runtime_error if the scatter database cannot be loaded.
   */
  void build( const std::vector<double> &energies,
              const std::function<double(double)> &eps_totint,
              const double fracH,
              const std::string &data_directory );

  bool valid() const { return !m_energies.empty(); }

  /** A at (energy, AN, AD).  T is double or an AD dual number (e.g.
   ceres::Jet); AN/AD clamp to the grid (d/dparam = 0 outside, the correct
   derivative of a hard clamp).  `areal_density` in g/cm2.  Returns T(0) for
   an un-built table or an energy that was not in `energies`.
   */
  template <class T>
  T evaluate( const double energy, const T &atomic_number,
              const T &areal_density_gcm2 ) const
  {
    if( m_energies.empty() )
      return T(0.0);

    // Exact-set energy lookup (the partner energies are known up front): the first sorted entry
    //  within 0.25 keV.  Bisect to a start safely below any match, then scan forward with that test.
    size_t ei = m_energies.size();
    const size_t start = static_cast<size_t>( std::lower_bound( m_energies.begin(), m_energies.end(),
                                                                energy - 0.5 ) - m_energies.begin() );
    for( size_t i = start; (i < m_energies.size()) && (m_energies[i] <= (energy + 0.5)); ++i )
    {
      if( std::fabs(m_energies[i] - energy) < 0.25 )
      {
        ei = i;
        break;
      }
    }
    if( ei >= m_energies.size() )
      return T(0.0);

    // Clamp into the grid; constants zero the derivative lanes (hard clamp).
    T an = atomic_number;
    if( an < m_an_grid.front() )
      an = T( m_an_grid.front() );
    if( an > m_an_grid.back() )
      an = T( m_an_grid.back() );

    T ad = areal_density_gcm2;
    if( ad < m_ad_grid.front() )
      ad = T( m_ad_grid.front() );
    if( ad > m_ad_grid.back() )
      ad = T( m_ad_grid.back() );

    size_t ai = 0;
    while( (ai + 2) < m_an_grid.size() && !(an < m_an_grid[ai+1]) )
      ++ai;
    size_t di = 0;
    while( (di + 2) < m_ad_grid.size() && !(ad < m_ad_grid[di+1]) )
      ++di;

    const T u = (an - m_an_grid[ai]) / (m_an_grid[ai+1] - m_an_grid[ai]);
    const T v = (ad - m_ad_grid[di]) / (m_ad_grid[di+1] - m_ad_grid[di]);

    const double a00 = valueAt( ei, ai, di );
    const double a10 = valueAt( ei, ai+1, di );
    const double a01 = valueAt( ei, ai, di+1 );
    const double a11 = valueAt( ei, ai+1, di+1 );

    return (1.0 - u)*(1.0 - v)*a00 + u*(1.0 - v)*a10
           + (1.0 - u)*v*a01 + u*v*a11;
  }//evaluate(...)

private:
  double valueAt( const size_t ei, const size_t ai, const size_t di ) const
  {
    return m_table[ (ei*m_an_grid.size() + ai)*m_ad_grid.size() + di ];
  }

  std::vector<double> m_energies;   //every energy evaluate() may be asked for
  std::vector<double> m_an_grid;
  std::vector<double> m_ad_grid;
  std::vector<double> m_table;      //[energy][an][ad], flattened
};//class ShieldScatterAugment


/** True-coincidence (cascade) summing corrections for the Activity/Shielding
 fit.  A thin adapter around CeeLo's analytic summing engine
 (ceelo::compute_cascade_analytic - the exact level-scheme joint-probability
 treatment): this class owns the per-nuclide cascade enumeration (SandiaDecay
 adapter, cached per age bucket), the fit's peak windows, and the GADRAS
 shield-scatter augmentation table; the caller supplies the absolute
 full-energy-peak and total efficiencies at its geometry (which is where all
 fit-parameter dependence enters - as `ceres::Jet` on the Ceres path, so
 d(correction)/d(shield thickness) flows into the Jacobian).

 Built once per fit setup (on the session thread, in
 ShieldingSourceChi2Fcn::create); afterwards only const methods are called
 (from fit worker threads) - the cascade cache is mutex-guarded.
 */
class CascadeSummingCalc
{
public:
  /** @param nuclide_initial_ages  The fits point-source nuclides and their
          starting ages (PhysicalUnits) - used to pre-enumerate cascades and
          collect the partner-energy set.
   @param peak_energy_widths  ShieldingSourceChi2Fcn::observedPeakEnergyWidths
          output: (gamma energy, peak sigma) per fit peak.
   @param photopeak_cluster_sigma  The fits clustering half-width in sigma.
   @param drf  The detector response; must satisfy #drfHasNeededInfo.
   @param shield_fracH  AD-weighted hydrogen mass fraction of the shield stack.
   @param shields_present  Whether to build the scatter-augmentation table.
   @param data_directory  InterSpec static data dir (sandia.shieldscatter.db).
   @param age_fit_nuclides  Those of the nuclides whose age may change during the fit; their
          partner energies cover every member of their decay chain (see #partnerEnergies), the
          others' only what their cascades at the given age emit.

   Throws std::runtime_error on missing DRF info or scatter-db load failure.
   */
  CascadeSummingCalc( const std::vector<std::pair<const SandiaDecay::Nuclide *,double>> &nuclide_initial_ages,
                      const std::vector<std::pair<double,double>> &peak_energy_widths,
                      const double photopeak_cluster_sigma,
                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                      const double shield_fracH,
                      const bool shields_present,
                      const std::string &data_directory,
                      const std::vector<const SandiaDecay::Nuclide *> &age_fit_nuclides = {} );

  /** Re-windowing copy: same physics as `other` (cascade enumeration, partner energies and the
   shield-scatter table are all peak-independent), but for a different set of fit peaks.

   Every expensive part of the primary constructor - `cascades(...)` per nuclide, the daughter
   x-ray line enumeration, the partner-energy detector constants, and
   `ShieldScatterAugment::build(...)` - is reused, including `other`s memoized cascade cache.
   Only the peak windows are rebuilt.

   Use this rather than the shared_ptr when the peak list differs: #evaluate returns one result
   per *stored* window, so a calculator built for a smaller peak set silently applies no
   correction to the peaks it has no window for.
   */
  CascadeSummingCalc( const CascadeSummingCalc &other,
                      const std::vector<std::pair<double,double>> &peak_energy_widths,
                      const double photopeak_cluster_sigma );

  /** Every emission energy (keV, sorted) a correction for these nuclides can query efficiencies
   at, at ANY age - cascade gammas, daughter K/L x-rays, and 511 of every member of each nuclide's
   decay chain; a superset of #allPartnerEnergies for a calculator built for them.  The ages are
   ignored (kept for the call sites).  Enumerates the cascades (not cached), so call it once, off
   the per-evaluation path.
   */
  static std::vector<double> partnerEnergies(
              const std::vector<std::pair<const SandiaDecay::Nuclide *,double>> &nuclide_ages );

  /** Whether the DRF carries the total-efficiency info corrections need. */
  static bool drfHasNeededInfo( const std::shared_ptr<const DetectorPeakResponse> &drf );

  /** The age bucket #cascades() keys its enumeration on - ~2% steps in log-age, so fitting an age
   does not force a re-enumeration every iteration.  Exposed so callers that memoise anything
   derived from a cascade set key on the same bucket: two ages in one bucket produce an identical
   set and therefore an identical correction, while keying on the raw age both never hits and grows
   without bound during an age fit.
   */
  static int ageCacheBucket( const double age );

  /** Cheap pre-fit estimate of the maximum possible summing magnitude
   |1 - C|, for GUI gating: fractionalSolidAngle x max total intrinsic
   efficiency (sampled 150-1500 keV).  For fixed-geometry DRFs, just the max
   total efficiency.  Returns 0 if the DRF lacks total-efficiency info.
   */
  static double estimateMaxSummingMagnitude( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                             const double distance );

  /** The fits peak windows (energy + clustering tolerance), in the order
   #evaluate returns results.
   */
  const std::vector<ceelo::PeakWindow> &peakWindows() const { return m_windows; }

  /** Every emission energy (gammas, K/L x-ray lines, 511) an evaluation may
   query the efficiency functors / scatter table at.
   */
  const std::vector<double> &allPartnerEnergies() const { return m_partner_energies; }

  /** The detector's far-field intrinsic FEP and total efficiencies at each #allPartnerEnergies
   energy (same order), tabulated once at construction from the DRF given to it: each total-
   efficiency query of a CeeLo response ray-traces a fresh aperture quadrature, and these are
   needed for every partner on every volumetric evaluation.
   */
  const std::vector<double> &partnerFepIntrinsic() const { return m_partner_fep_int; }
  const std::vector<double> &partnerTotIntrinsic() const { return m_partner_tot_int; }

  /** Index of the entry of the sorted `energies` (e.g. #allPartnerEnergies) nearest `energy`, or -1 if none is
   within 0.25 keV.

   It must be the nearest entry: another source's line can sit just 0.25 keV away (U-238's 221.15 keV below
   U-235's 221.40 keV), so taking the first entry at or above `energy - 0.25` can land on that line, fail the
   tolerance, and report no efficiency at all for the exact match just after it.
   */
  static int nearestEnergyIndex( const std::vector<double> &energies, const double energy )
  {
    const std::vector<double>::const_iterator pos = std::lower_bound( energies.begin(), energies.end(), energy );
    int index = -1;
    double diff = 0.25;
    if( (pos != energies.end()) && (std::fabs(*pos - energy) < diff) )
    {
      index = static_cast<int>( pos - energies.begin() );
      diff = std::fabs(*pos - energy);
    }
    if( (pos != energies.begin()) && (std::fabs(*(pos - 1) - energy) < diff) )
      index = static_cast<int>( (pos - 1) - energies.begin() );
    return index;
  }

  const ShieldScatterAugment &scatterAugment() const { return m_scatter; }

  /** Detector absolute FEP / total efficiency at (energy, theta, phi,
   distance), CeeLo-response-aware (near-field / off-axis) with legacy-curve
   fallback - the "bare" part of the efficiency functors (no shields/air).
   */
  static double detectorFepEffAbs( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                   const double energy, const double theta,
                                   const double phi, const double distance );
  static double detectorTotEffAbs( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                   const double energy, const double theta,
                                   const double phi, const double distance );

  /** Wraps an absolute-efficiency function in a cached, log-spaced interpolation
   table, for callers that need an efficiency at many energies.

   A detector-response evaluation costs on the order of half a millisecond, and
   cascade summing asks for an efficiency at every emitted line and every cascade
   partner energy - about a thousand of each for a decay chain, i.e. roughly a
   second of pure response evaluation per update.  The table uses 256 log-spaced
   points between 5 keV and 3.5 MeV, interpolated in log-log; measured against the
   example HPGe response that is within 0.31% everywhere and 0.06% above 60 keV,
   which is far inside what a qualitative display or a summing correction needs.
   Energies outside the tabulated range fall through to `fcn`.

   Points are evaluated on first use, and tables are cached (the last few) under
   `key`, so repeated calls with the same detector and geometry are nearly free.
   `key` must therefore capture everything `fcn` depends on - detector identity,
   distance, and any fallback tier.  Thread-safe.
   */
  static std::function<double(double)> interpolatedEff( const std::string &key,
                                                        std::function<double(double)> fcn );

  /** The cascades for `nuc` at `age` (PhysicalUnits); enumerated on first use
   per (nuclide, age bucket) and cached (thread-safe).
   */
  std::shared_ptr<const std::vector<ceelo::DecayCascade>> cascades(
                    const SandiaDecay::Nuclide *nuc, const double age ) const;

  /** The prepared engine model of #cascades (level-scheme DAGs, occurrence lists, ...), cached
   alongside them, so evaluations do not re-prepare it.
   */
  std::shared_ptr<const ceelo::AnalyticCascadeModel> cascadeModel(
                    const SandiaDecay::Nuclide *nuc, const double age ) const;

  /** Per-window summing results for one nuclide, given the absolute
   efficiencies at the fit geometry.  T = double, or ceres::Jet on the Ceres
   path (the functors then carry the parameter derivatives).

   Windows with no emitted line of this nuclide come back `found == false`
   (c_net == 1) - skip those.  Throws if the engine reports unmatched
   energies.
   */
  template <class T>
  std::vector<ceelo::AnalyticPeakResultT<T>> evaluate(
                      const SandiaDecay::Nuclide *nuc,
                      const double age,
                      const std::function<T(double)> &eps_fep_abs,
                      const std::function<T(double)> &eps_tot_abs ) const
  {
    return evaluateWindows( nuc, age, m_windows, eps_fep_abs, eps_tot_abs );
  }

  /** Same as #evaluate, but for an explicit window subset - e.g. the volume
   integrand evaluates a single calculator's window per element.

   With `report_pure_sums`, a window holding no line of the nuclide also gets its pure coincidence-sum rate
   (`ceelo::AnalyticCascadeOptions::report_pure_sums`).
   */
  template <class T>
  std::vector<ceelo::AnalyticPeakResultT<T>> evaluateWindows(
                      const SandiaDecay::Nuclide *nuc,
                      const double age,
                      const std::vector<ceelo::PeakWindow> &windows,
                      const std::function<T(double)> &eps_fep_abs,
                      const std::function<T(double)> &eps_tot_abs,
                      const bool report_pure_sums = false ) const
  {
    struct FunctorProv final : public ceelo::EfficiencyProviderT<T>
    {
      const std::function<T(double)> *fep_f = nullptr;
      const std::function<T(double)> *tot_f = nullptr;
      T fep( double e ) const override { return (*fep_f)( e ); }
      T total( double e ) const override { return (*tot_f)( e ); }
      bool has( double ) const override { return true; }
    };

    const std::shared_ptr<const ceelo::AnalyticCascadeModel> model = cascadeModel( nuc, age );

    FunctorProv prov;
    prov.fep_f = &eps_fep_abs;
    prov.tot_f = &eps_tot_abs;

    ceelo::AnalyticCascadeOptions opts;  //defaults: triples + W(0) angular on
    opts.report_pure_sums = report_pure_sums;

    std::vector<ceelo::AnalyticPeakResultT<T>> results
              = ceelo::compute_cascade_analytic( *model, windows, prov, opts );

    for( const ceelo::AnalyticPeakResultT<T> &r : results )
    {
      if( !r.unmatched_energies.empty() )
        throw std::runtime_error( "CascadeSummingCalc: efficiency undefined at "
                    + std::to_string(r.unmatched_energies[0]) + " keV" );
    }

    return results;
  }//evaluate(...)

private:
  std::vector<ceelo::PeakWindow> m_windows;
  std::vector<double> m_partner_energies;
  std::vector<double> m_partner_fep_int, m_partner_tot_int;
  ShieldScatterAugment m_scatter;

  //Cascade cache: keyed by (nuclide, age bucket).  Age enters the branch
  // weights only through progeny in-growth, so ages are bucketed (~2% steps in
  // log age) to avoid re-enumerating per fit iteration when the age is being
  // fit.  The first age to reach a bucket defines its entry (short-lived
  // progeny grow in proportion to age, so weights can vary by ~2% within a
  // bucket); every evaluation in the bucket then uses that one entry.
  struct CachedCascades
  {
    std::shared_ptr<const std::vector<ceelo::DecayCascade>> cascades;
    std::shared_ptr<const ceelo::AnalyticCascadeModel> model;  //owns a reference to `cascades`
  };

  /** The cache entry for (nuc, bucket of age), enumerating and preparing it on first use. */
  CachedCascades cachedCascades( const SandiaDecay::Nuclide *nuc, const double age ) const;

  mutable std::mutex m_cascade_mutex;
  mutable std::map<std::pair<const SandiaDecay::Nuclide *,int>, CachedCascades> m_cascades;
};//class CascadeSummingCalc

}//namespace GammaInteractionCalc

#endif //CascadeSummingCalc_h
