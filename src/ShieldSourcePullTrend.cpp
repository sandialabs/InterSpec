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

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include <Eigen/Dense>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/CascadeSummingCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldSourcePullTrend.h"
#include "InterSpec/ShieldingSourceFitCalc.h"

using namespace std;

namespace
{
  // ---- Decision (profile likelihood-ratio) constants; calibrated by scratch/20261111_pull_trend_study
  //  (BASELINE.md): on real FitPeaksForNuclides fits these give ~2% false alarms on correct models,
  //  and catch ~75% (HPGe) / ~55% (NaI) of visible shielding-amount errors.

  /** Fractional floor added in quadrature to each peaks uncertainty (times its model counts), so a
   few very high-statistics peaks - whose areas are rarely good to better than a percent or two -
   cannot decide the diagnosis on their own.
   */
  const double sm_floor_frac = 0.02;

  /** Prior minimum residual scale: real peak fits scatter 2-5x their statistical uncertainty (the
   study measured ~1.3x for automated HPGe fits and ~4-5x for NaI), so the statistics are never
   divided by less than this, even when the fit leaves too few degrees of freedom to estimate it.
   */
  const double sm_min_residual_scale = 2.0;

  /** |t| (normal-equivalent sigma) needed for an amount or firm atomic-number call. */
  const double sm_call_t = 3.5;

  /** Significance needed for the non-specific "other" call.  Higher than the physical calls: with
   no physical direction to check it against, it otherwise fires on routine peak-fit systematics
   (CZT automated fits, whose areas run low by an energy-dependent amount, gave 11% false "other"
   at 3.5 sigma vs 4% at 5).
   */
  const double sm_other_t = 5.0;

  /** Soft ("possibly") atomic-number calls down to this |t|, only when the data have the power. */
  const double sm_soft_z_t = 2.5;

  /** Power gate for soft atomic-number calls; see TrendResult::anDiscriminationT. */
  const double sm_an_power_gate = 2.5;

  /** The reference effective-Z error the power is computed for (a factor, either direction). */
  const double sm_an_reference_z_factor = 1.5;

  /** Confidence tiers (|t|): weak below the first, moderate, then strong. */
  const double sm_moderate_t = 4.0;
  const double sm_strong_t = 6.0;

  /** Chi2 margin (in scaled units) by which the best reference-Z layer must beat every other one
   before a "model lacks shielding" message names a Z class.
   */
  const double sm_zclass_margin = 4.0;

  /** Reference layers tried when the model has no shielding: water/plastic, Al, Fe, Sn, Pb like. */
  const std::array<double,5> sm_reference_z = { 7.0, 13.0, 26.0, 50.0, 82.0 };

  /** A "too much shielding" call whose best-fit total areal density is below this says "little or
   no shielding" instead of quoting a near-zero areal density.
   */
  const double sm_little_or_no_shielding = 1.0 * PhysicalUnits::g / PhysicalUnits::cm2;

  /** Cascade-summing caveat: only for sources closer than this to the detector face. */
  const double sm_summing_max_distance = 15.0 * PhysicalUnits::cm;

  /** ...and only when the detector's maximum summing magnitude is at least this (when estimable). */
  const double sm_summing_min_magnitude = 0.01;


  // ---- Descriptive trend-curve (polynomial) constants ----

  /** |t| of the polynomial slope/curvature for the descriptive curve to be drawn when no
   conclusion was reached.
   */
  const double sm_curve_t_threshold = 3.5;

  /** Pulls beyond this many sigma are winsorized for the descriptive polynomial. */
  const double sm_max_abs_pull = 10.0;

  /** Minimum fraction of a peaks model counts from its assigned nuclide, for the polynomial. */
  const double sm_min_nuclide_fraction = 0.7;

  /** Number of samples of the fitted curve returned per-nuclide, for plotting. */
  const size_t sm_num_curve_samples = 40;


  /** A pull point accepted into the descriptive polynomial. */
  struct PullPoint
  {
    double energy;      //keV
    double pull;        //winsorized (observed-model)/uncert
    size_t nucIndex;    //index into the nuclide grouping
  };//struct PullPoint


  /** Mass attenuation coefficient (cm2/g) of a generic element at an energy (keV). */
  double mu_generic( const double z, const double energy )
  {
    const double one_gcm2 = PhysicalUnits::g / PhysicalUnits::cm2;
    return GammaInteractionCalc::transmition_coefficient_generic( std::min( 98.0, std::max( 1.0, z ) ),
                                                                  one_gcm2, energy );
  }


  /** Normal-equivalent sigma of a chi2 value with 2 degrees of freedom (survival exp(-x/2)). */
  double z_from_chi2_2dof( const double x )
  {
    if( !(x > 0.0) )
      return 0.0;
    const double p = std::exp( -0.5*x );
    if( !(p > 0.0) )
      return 38.0;
    double lo = 0.0, hi = 38.0;  //bisection on the normal upper tail
    for( int i = 0; i < 60; ++i )
    {
      const double mid = 0.5*(lo + hi);
      if( 0.5*std::erfc( mid/std::sqrt(2.0) ) > p )
        lo = mid;
      else
        hi = mid;
    }
    return 0.5*(lo + hi);
  }//z_from_chi2_2dof(...)


  /** The per-peak data of the profile likelihood, and the activities-profiled chi2 of any change
   of the per-peak shield optical depth.  The model of each peak is sum_k a_k * m_ik * exp(-dtau_i),
   with the activity multipliers a_k re-solved (linear least squares) for every dtau.
   */
  struct ProfileData
  {
    size_t n = 0, K = 0;
    std::vector<double> energy, w, y, tau;   //w = 1/sigma', y = w*observed
    std::vector<double> logm;                //n x K, log of the model contribution (-inf when 0)

    // Work buffers (chi2 is called a few thousand times).
    mutable Eigen::MatrixXd G;
    mutable Eigen::VectorXd yv;

    double chi2( const std::vector<double> &dtau, std::vector<double> *model = nullptr ) const
    {
      G.resize( n, K );
      yv.resize( n );
      for( size_t i = 0; i < n; ++i )
      {
        yv(i) = y[i];
        for( size_t j = 0; j < K; ++j )
        {
          const double lm = logm[i*K + j];
          G(i,j) = std::isfinite(lm) ? (w[i] * std::exp( std::min( 700.0, lm - dtau[i] ) )) : 0.0;
        }
      }

      // Pivoted QR, not the normal equations: when a nuclide's lines are attenuated to almost
      //  nothing its column is tiny, and squaring the condition number then corrupts the solution.
      const Eigen::VectorXd a = G.colPivHouseholderQr().solve( yv );
      double c2 = 0.0;
      if( model )
        model->resize( n );
      for( size_t i = 0; i < n; ++i )
      {
        double fit = 0.0;
        for( size_t j = 0; j < K; ++j )
          fit += G(i,j) * (std::isfinite(a(j)) ? a(j) : 0.0);
        c2 += (y[i] - fit)*(y[i] - fit);
        if( model )
          (*model)[i] = (w[i] > 0.0) ? (fit / w[i]) : 0.0;
      }
      return c2;
    }//chi2(...)

    std::vector<double> mu( const double z ) const
    {
      std::vector<double> v( n );
      for( size_t i = 0; i < n; ++i )
        v[i] = mu_generic( z, energy[i] );
      return v;
    }
  };//struct ProfileData


  /** Minimizes f over a grid, then refines between the best grid points neighbours. */
  template<class F>
  double minimize_1d( const F &f, const std::vector<double> &grid, double &xbest )
  {
    size_t ib = 0;
    double best = std::numeric_limits<double>::infinity();
    for( size_t i = 0; i < grid.size(); ++i )
    {
      const double v = f( grid[i] );
      if( v < best )
      {
        best = v;
        ib = i;
      }
    }
    xbest = grid[ib];
    double lo = grid[(ib > 0) ? (ib - 1) : 0];
    double hi = grid[std::min( grid.size() - 1, ib + 1 )];
    const double gr = 0.5*(std::sqrt(5.0) - 1.0);
    double c = hi - gr*(hi - lo), d = lo + gr*(hi - lo);
    double fc = f( c ), fd = f( d );
    for( int it = 0; it < 30; ++it )
    {
      if( fc < fd )
      {
        hi = d;
        d = c;
        fd = fc;
        c = hi - gr*(hi - lo);
        fc = f( c );
      }else
      {
        lo = c;
        c = d;
        fc = fd;
        d = lo + gr*(hi - lo);
        fd = f( d );
      }
    }//for( golden-section iterations )
    const double xm = 0.5*(lo + hi), fm = f( xm );
    if( fm < best )
    {
      best = fm;
      xbest = xm;
    }
    return best;
  }//minimize_1d(...)


  /** Outcome of the profile likelihood alternatives for one data set. */
  struct ProfileResult
  {
    bool ok = false;
    double chi2_0 = 0, chi2_a = 0, chi2_b = 0, chi2_best = 0, chi2_other = 0;
    bool added_layer = false;     //amount alternative was an added reference layer (no model shielding)
    double u_best = 1.0;          //amount alternative scale of the model optical depth
    double add_z = 0, add_ad = 0; //added layer Z and AD (g/cm2)
    double add_fe_ad = 0;         //added layer AD (g/cm2) if it is Fe-like (Z 26)
    bool add_z_resolved = false;
    double z_best = 0, ad_best = 0;   //free-Z alternative (AD in g/cm2)
    bool z_beat_amount = false;   //the free-Z scan improved on the amount alternative
    int p_best = 0;               //parameters of the best alternative
  };//struct ProfileResult


  /** Runs the amount alternative (always), the free-Z alternative (if `do_z`), and the "other"
   polynomial beyond the best of them.  `z0`/`ad0_gcm2` describe the model shielding (total areal
   density along the source-detector line); `model_has_shield` false means the model has no shielding,
   so the amount alternative adds a layer (any optical depth the peaks report, e.g. air folded into a
   beside-the-detector efficiency, is not shielding to scale).
   */
  ProfileResult run_profile( const ProfileData &P, const bool model_has_shield, const bool do_z,
                             const double z0, const double ad0_gcm2 )
  {
    ProfileResult res;
    if( (P.K < 1) || (P.n < P.K + 2) )
      return res;
    res.ok = true;

    const std::vector<double> zero( P.n, 0.0 );
    res.chi2_0 = P.chi2( zero );
    res.chi2_a = res.chi2_b = res.chi2_best = res.chi2_0;

    const double max_tau = *std::max_element( begin(P.tau), end(P.tau) );
    const bool have_tau = model_has_shield && (max_tau > 1.0E-6);

    static const std::vector<double> s_ad_grid = { 0.0, 0.1, 0.3, 0.6, 1, 1.5, 2, 3, 4.5, 6.5, 9, 12, 16, 21,
                                                   28, 37, 50, 65, 85, 110, 150, 200, 300, 400, 500 };
    std::vector<double> dtau( P.n, 0.0 );

    if( have_tau )
    {
      // Scale the model shielding's optical depth by u; for a thin model shield, far enough to reach
      //  ~500 g/cm2 (the largest areal density a fit allows).
      static const std::vector<double> s_ugrid = { 0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 1.0,
                                                   1.1, 1.25, 1.5, 1.8, 2.2, 2.7, 3.5, 5.0, 7.0, 10.0, 15.0, 25.0 };
      std::vector<double> ugrid = s_ugrid;
      const double u_max = (ad0_gcm2 > 0.0) ? (500.0 / ad0_gcm2) : 25.0;
      while( (ugrid.back() < u_max) && (ugrid.size() < 80) )
        ugrid.push_back( 1.6*ugrid.back() );

      res.chi2_a = minimize_1d( [&]( const double u ){
        for( size_t i = 0; i < P.n; ++i )
          dtau[i] = (u - 1.0)*P.tau[i];
        return P.chi2( dtau );
      }, ugrid, res.u_best );
    }else
    {
      // No model shielding: add a reference layer, AD >= 0.
      std::array<double,sm_reference_z.size()> c2s;
      double best = std::numeric_limits<double>::infinity();
      for( size_t iz = 0; iz < sm_reference_z.size(); ++iz )
      {
        const std::vector<double> mu = P.mu( sm_reference_z[iz] );
        double adb = 0.0;
        c2s[iz] = minimize_1d( [&]( const double ad ){
          for( size_t i = 0; i < P.n; ++i )
            dtau[i] = std::max( 0.0, ad )*mu[i];
          return P.chi2( dtau );
        }, s_ad_grid, adb );
        adb = std::max( 0.0, adb );
        if( sm_reference_z[iz] == 26.0 )
          res.add_fe_ad = adb;
        if( c2s[iz] < best )
        {
          best = c2s[iz];
          res.add_z = sm_reference_z[iz];
          res.add_ad = adb;
        }
      }
      res.chi2_a = std::min( best, res.chi2_0 );
      res.added_layer = true;

      // The Z class is only meaningful if the best layer beats every other one clearly.
      double runner_up = std::numeric_limits<double>::infinity();
      for( size_t iz = 0; iz < sm_reference_z.size(); ++iz )
        if( sm_reference_z[iz] != res.add_z )
          runner_up = std::min( runner_up, c2s[iz] );
      res.add_z_resolved = ((runner_up - best) >= sm_zclass_margin*sm_min_residual_scale*sm_min_residual_scale);
    }
    res.chi2_best = std::min( res.chi2_best, res.chi2_a );
    res.p_best = 1;

    // Per-peak ratio of the model optical depth to generic mu(z0)*ad0, so the free-Z alternative
    //  anchors exactly on the current model (FEP-removal coefficient and chord vs plain mu).
    std::vector<double> corr( P.n, 1.0 );
    if( do_z && have_tau && (z0 > 0.0) && (ad0_gcm2 > 0.0) )
    {
      const std::vector<double> mu0 = P.mu( z0 );
      for( size_t i = 0; i < P.n; ++i )
        corr[i] = (mu0[i]*ad0_gcm2 > 0.0) ? (P.tau[i] / (mu0[i]*ad0_gcm2)) : 1.0;

      static const std::vector<double> s_zgrid = { 4, 6, 7, 8, 10, 12, 13, 15, 18, 22, 26, 30, 35, 42, 50, 58,
                                                   66, 74, 82, 90, 94 };
      std::vector<double> zs = s_zgrid;
      zs.push_back( std::min( 94.0, std::max( 4.0, z0 ) ) );
      std::sort( begin(zs), end(zs) );

      double best = std::numeric_limits<double>::infinity();
      for( const double z : zs )
      {
        std::vector<double> mu = P.mu( z );
        for( size_t i = 0; i < P.n; ++i )
          mu[i] *= corr[i];
        double adb = 0.0;
        const double c2 = minimize_1d( [&]( const double ad ){
          for( size_t i = 0; i < P.n; ++i )
            dtau[i] = std::max( 0.0, ad )*mu[i] - P.tau[i];
          return P.chi2( dtau );
        }, s_ad_grid, adb );
        if( c2 < best )
        {
          best = c2;
          res.z_best = z;
          res.ad_best = std::max( 0.0, adb );
        }
      }
      res.z_beat_amount = (best < res.chi2_a);
      res.chi2_b = std::min( best, res.chi2_a );
      res.chi2_best = std::min( res.chi2_best, res.chi2_b );
      res.p_best = 2;
    }//if( free-Z alternative )

    // "Other": a shared quadratic-in-log(E) efficiency multiplier beyond the best physical alternative.
    std::fill( begin(dtau), end(dtau), 0.0 );
    if( res.z_beat_amount )
    {
      const std::vector<double> mu = P.mu( res.z_best );
      for( size_t i = 0; i < P.n; ++i )
        dtau[i] = res.ad_best*mu[i]*corr[i] - P.tau[i];
    }else if( res.added_layer )
    {
      const std::vector<double> mu = P.mu( (res.add_z > 0.0) ? res.add_z : 26.0 );
      for( size_t i = 0; i < P.n; ++i )
        dtau[i] = res.add_ad*mu[i];
    }else
    {
      for( size_t i = 0; i < P.n; ++i )
        dtau[i] = (res.u_best - 1.0)*P.tau[i];
    }

    res.chi2_other = res.chi2_best;
    const double emin = *std::min_element( begin(P.energy), end(P.energy) );
    const double emax = *std::max_element( begin(P.energy), end(P.energy) );
    const double xmid = 0.5*(std::log(emin) + std::log(emax));
    const double xh = 0.5*(std::log(emax) - std::log(emin));
    if( (xh > 0.0) && (P.n >= P.K + 3 + static_cast<size_t>(res.p_best)) )
    {
      std::vector<double> best_model;
      P.chi2( dtau, &best_model );
      const size_t ncol = P.K + 2;
      Eigen::MatrixXd G( P.n, ncol );
      Eigen::VectorXd yv( P.n );
      for( size_t i = 0; i < P.n; ++i )
      {
        const double x = (std::log(P.energy[i]) - xmid)/xh;
        for( size_t j = 0; j < P.K; ++j )
        {
          const double lm = P.logm[i*P.K + j];
          G(i,j) = std::isfinite(lm) ? (P.w[i] * std::exp( std::min( 700.0, lm - dtau[i] ) )) : 0.0;
        }
        G(i,P.K) = P.w[i] * best_model[i] * x;
        G(i,P.K+1) = P.w[i] * best_model[i] * x*x;
        yv(i) = P.y[i];
      }
      const Eigen::VectorXd a = G.colPivHouseholderQr().solve( yv );
      const double c2 = (yv - G*a).squaredNorm();
      if( std::isfinite(c2) )
        res.chi2_other = std::min( res.chi2_best, c2 );
    }

    return res;
  }//run_profile(...)


  ShieldSourcePullTrend::TrendConfidence confidence_from_t( const double abs_t )
  {
    using ShieldSourcePullTrend::TrendConfidence;
    if( abs_t >= sm_strong_t )
      return TrendConfidence::Strong;
    if( abs_t >= sm_moderate_t )
      return TrendConfidence::Moderate;
    if( abs_t >= sm_call_t )
      return TrendConfidence::Weak;
    return TrendConfidence::Tentative;
  }//confidence_from_t(...)


  /** Whether a source of nuclide `nuc` at age `age` emits coincident photons that can sum: some
   member of its decay chain, with at least 1% of the parent's activity at that age, has a decay mode
   (branch >= 1%) with two gammas above 30 keV (each >= 0.5% per decay) linked through a level of the
   daughter's level scheme (one feeds the level the other de-excites - the same topology the cascade
   engine uses), or a gamma together with a positron (gamma + annihilation-photon summing, e.g. Na22).
   Without level data the decay mode instead needs more than one gamma above 30 keV per decay on
   average.  Progeny behind a long-lived intermediate (e.g. Am241 -> Np237 -> ...) or fed only by a
   minor branch are not counted.  The level link, not the multiplicity, catches low-yield cascades
   such as Lu177's 208 -> 113 keV.
   */
  bool has_coincident_gammas( const SandiaDecay::Nuclide *nuc, const double age )
  {
    if( !nuc )
      return false;

    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuc, 1.0 );
    const double parent_act = mix.activity( std::max( 0.0, age ), nuc );

    std::vector<const SandiaDecay::Nuclide *> chain = nuc->descendants();
    if( std::find( begin(chain), end(chain), nuc ) == end(chain) )
      chain.push_back( nuc );

    for( const SandiaDecay::Nuclide *n : chain )
    {
      if( !n )
        continue;
      const double act = mix.activity( std::max( 0.0, age ), n );
      if( !(parent_act > 0.0) || !(act >= 0.01*parent_act) )
        continue;

      for( const SandiaDecay::Transition *trans : n->decaysToChildren )
      {
        if( !trans || (trans->branchRatio < 0.01) )
          continue;
        const std::vector<SandiaDecay::RadParticle> &products = trans->products;
        const auto is_summable_gamma = []( const SandiaDecay::RadParticle &p ) -> bool {
          return (p.type == SandiaDecay::GammaParticle) && (p.energy > 30.0) && (p.intensity >= 0.005);
        };

        double gammas = 0.0, positrons = 0.0;
        bool have_levels = false;
        for( const SandiaDecay::RadParticle &p : products )
        {
          if( (p.type == SandiaDecay::GammaParticle) && (p.energy > 30.0) )
          {
            gammas += p.intensity;
            have_levels = (have_levels || (p.from_level >= 0));
          }else if( p.type == SandiaDecay::PositronParticle )
          {
            positrons += p.intensity;
          }
        }//for( products )

        if( (gammas > 0.05) && (positrons > 0.05) )
          return true;

        if( !have_levels )
        {
          if( gammas > 1.05 )
            return true;
          continue;
        }

        for( const SandiaDecay::RadParticle &feeder : products )
        {
          if( !is_summable_gamma( feeder ) || (feeder.to_level < 0) )
            continue;
          for( const SandiaDecay::RadParticle &follower : products )
          {
            if( (&follower != &feeder) && is_summable_gamma( follower )
                && (follower.from_level == feeder.to_level) )
              return true;
          }
        }//for( feeder gammas )
      }//for( transitions )
    }//for( chain )
    return false;
  }//has_coincident_gammas(...)


  /** Areal density (PhysicalUnits) a shielding layer presents along the source-detector line, for a
   point source at the center: the thickness along the ray (spherical: [0]; cylinder end-on: the
   (half-)length [1]; cylinder side-on: the radius [0]; rectangular: the depth [2]) times the density,
   or the areal density of a generic shielding.
   */
  double along_ray_areal_density( const ShieldingSourceFitCalc::ShieldingInfo &s )
  {
    if( s.m_isGenericMaterial )
      return std::fabs( s.m_dimensions[1] );
    if( !s.m_material )
      return 0.0;

    double thickness = s.m_dimensions[0];
    switch( s.m_geometry )
    {
      case GammaInteractionCalc::GeometryType::CylinderEndOn: thickness = s.m_dimensions[1]; break;
      case GammaInteractionCalc::GeometryType::Rectangular:   thickness = s.m_dimensions[2]; break;
      case GammaInteractionCalc::GeometryType::Spherical:
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        break;
    }//switch( s.m_geometry )

    return std::fabs( thickness * s.m_material->density );
  }//along_ray_areal_density(...)


  /** Effective (areal-density-weighted) atomic number of the model shieldings, for the
   discrimination-power check.  Atomic number is unit-free, so for a single shielding this is
   robust to the AD unit conventions (the scale cancels in sum_z_ad/sum_ad).  Both entry points
   (final fit and live preview) store generic areal densities in PhysicalUnits, the same units as
   a material's thickness*density, so mixed generic + material models weight consistently too.
   Returns false if there is no usable shielding.
   */
  bool effective_atomic_number( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                double &eff_z )
  {
    double sum_ad = 0.0, sum_z_ad = 0.0;
    for( const ShieldingSourceFitCalc::ShieldingInfo &s : shieldings )
    {
      double z = 0.0;
      if( s.m_isGenericMaterial )
        z = s.m_dimensions[0];
      else if( s.m_material )
        z = s.m_material->massWeightedAtomicNumber();
      else
        continue;
      const double ad = along_ray_areal_density( s );
      if( (ad > 0.0) && (z > 0.0) )
      {
        sum_ad += ad;
        sum_z_ad += z*ad;
      }
    }
    if( sum_ad <= 0.0 )
      return false;
    eff_z = sum_z_ad / sum_ad;
    return true;
  }//effective_atomic_number(...)


  /** Evaluates the regressor transform g(E) for the requested basis.
   Returns NaN if the transform is not evaluable (caller falls back to LogEnergy).
   */
  double basis_transform( const ShieldSourcePullTrend::TrendBasis basis,
                          const double energy,
                          const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings )
  {
    using ShieldSourcePullTrend::TrendBasis;

    switch( basis )
    {
      case TrendBasis::Energy:
        return energy;

      case TrendBasis::LogEnergy:
        return (energy > 0.0) ? std::log(energy) : std::numeric_limits<double>::quiet_NaN();

      case TrendBasis::InvEnergy:
        return (energy > 0.0) ? (1.0/energy) : std::numeric_limits<double>::quiet_NaN();

      case TrendBasis::MuModel:
      {
        // Total optical depth (sum of mu*thickness) of the model shieldings at this energy.
        //  For non-spherical geometries the first dimension is used as the thickness - this is
        //  only a regressor, so the approximation just mildly distorts the basis.
        if( shieldings.empty() )
          return std::numeric_limits<double>::quiet_NaN();

        double optical_depth = 0.0;
        for( const ShieldingSourceFitCalc::ShieldingInfo &shield : shieldings )
        {
          try
          {
            if( shield.m_isGenericMaterial )
            {
              optical_depth += GammaInteractionCalc::transmition_coefficient_generic(
                                  shield.m_dimensions[0], shield.m_dimensions[1], energy );
            }else if( shield.m_material )
            {
              optical_depth += GammaInteractionCalc::transmition_coefficient_material(
                                  shield.m_material.get(), energy, shield.m_dimensions[0] );
            }
          }catch( std::exception & )
          {
            return std::numeric_limits<double>::quiet_NaN();
          }
        }//for( loop over shieldings )

        return optical_depth;
      }//case TrendBasis::MuModel
    }//switch( basis )

    assert( 0 );
    return std::numeric_limits<double>::quiet_NaN();
  }//basis_transform(...)


}//namespace


namespace ShieldSourcePullTrend
{

FitConfigSummary summarize_fit_config(
                    const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                    const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources )
{
  FitConfigSummary config;

  for( const ShieldingSourceFitCalc::ShieldingInfo &shield : shieldings )
  {
    if( shield.m_isGenericMaterial )
    {
      config.genericAnFit = (config.genericAnFit || shield.m_fitDimensions[0]);
      config.genericAdFit = (config.genericAdFit || shield.m_fitDimensions[1]);
    }else
    {
      for( size_t i = 0; i < 3; ++i )
        config.physicalDimensionFit = (config.physicalDimensionFit || shield.m_fitDimensions[i]);
    }

    if( !shield.m_nuclideFractions_.empty() || !shield.m_traceSources.empty() )
      config.hasSelfAttenOrTraceSrc = true;
  }//for( loop over shieldings )

  for( const ShieldingSourceFitCalc::SourceFitDef &src : sources )
  {
    if( src.sourceType != ShieldingSourceFitCalc::ModelSourceType::Point )
      config.hasSelfAttenOrTraceSrc = true;
    if( src.fitAge )
      config.anyAgeFit = true;
  }

  return config;
}//summarize_fit_config(...)


void set_summing_info( FitConfigSummary &config, const double distance, const bool summing_corrected,
                       const std::shared_ptr<const DetectorPeakResponse> &drf,
                       const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources )
{
  config.possibleCascadeSumming = false;
  try
  {
    if( summing_corrected )
      return;

    // A fixed-geometry DRF has no meaningful distance (the sample sits at its fixed, usually close,
    //  position); otherwise only close sources.
    const bool fixed_geom = (drf && drf->isFixedGeometry());
    if( !fixed_geom && (!(distance >= 0.0) || (distance >= sm_summing_max_distance)) )
      return;

    bool cascades = false;
    for( const ShieldingSourceFitCalc::SourceFitDef &src : sources )
      cascades = (cascades || has_coincident_gammas( src.nuclide, src.age ));
    if( !cascades )
      return;

    // When the DRF lets us estimate it, require a non-negligible summing magnitude; otherwise
    //  (no total-efficiency information) summing cannot be ruled out.
    if( drf && GammaInteractionCalc::CascadeSummingCalc::drfHasNeededInfo( drf ) )
    {
      const double magnitude = GammaInteractionCalc::CascadeSummingCalc::estimateMaxSummingMagnitude( drf, distance );
      if( magnitude < sm_summing_min_magnitude )
        return;
    }

    config.possibleCascadeSumming = true;
  }catch( std::exception & )
  {
    config.possibleCascadeSumming = false;
  }
}//set_summing_info(...)


TrendResult fit_pull_trend(
              const std::vector<GammaInteractionCalc::PeakResultPlotInfo> &peak_comparisons,
              const std::vector<GammaInteractionCalc::PeakDetail> &peak_details,
              const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
              const FitConfigSummary &config,
              const TrendBasis basis )
{
  TrendResult result;
  result.basis = basis;
  result.model = TrendModel::None;
  result.possibleCascadeSumming = config.possibleCascadeSumming;

  const size_t npoints = std::min( peak_comparisons.size(), peak_details.size() );

  // ---- Peaks for the profile likelihood (all valid peaks; per-nuclide model contributions) ----
  //  Two fit peaks assigned the same gamma share one model entry (the second gets no source or
  //  attenuation detail), so only the first is used.
  ProfileData prof;
  {
    vector<string> prof_nucs;
    vector<vector<double>> m_rows;
    vector<double> used_energies;
    for( size_t i = 0; i < npoints; ++i )
    {
      const GammaInteractionCalc::PeakResultPlotInfo &c = peak_comparisons[i];
      const GammaInteractionCalc::PeakDetail &d = peak_details[i];
      if( !std::isfinite(c.energy) || (c.energy <= 0.0) || !std::isfinite(c.observedCounts)
          || !std::isfinite(c.expectedCounts) || (c.expectedCounts < 0.0)
          || !std::isfinite(c.observedUncert) || (c.observedUncert <= 0.0) )
        continue;
      if( std::find( begin(used_energies), end(used_energies), c.energy ) != end(used_energies) )
        continue;
      used_energies.push_back( c.energy );

      const auto index_of = [&prof_nucs]( const string &name ) -> size_t {
        const vector<string>::const_iterator pos = std::find( prof_nucs.cbegin(), prof_nucs.cend(), name );
        if( pos != end(prof_nucs) )
          return static_cast<size_t>( pos - prof_nucs.cbegin() );
        prof_nucs.push_back( name );
        return prof_nucs.size() - 1;
      };

      vector<pair<size_t,double>> contribs;
      double sum = 0.0;
      for( const GammaInteractionCalc::PeakDetailSrc &src : d.m_sources )
      {
        if( !src.nuclide || !(src.modelContribToPeak > 0.0) )
          continue;
        contribs.emplace_back( index_of( src.nuclide->symbol ), src.modelContribToPeak );
        sum += src.modelContribToPeak;
      }
      if( contribs.empty() && (c.expectedCounts > 0.0) && !d.assignedNuclide.empty() )
      {
        contribs.emplace_back( index_of( d.assignedNuclide ), c.expectedCounts );
        sum = c.expectedCounts;
      }

      // Normalize to the reported model counts (they differ when e.g. summing corrections apply).
      //  The ratio first: a line attenuated to ~1e-270 counts would underflow as a product.
      vector<double> row( prof_nucs.size(), 0.0 );
      for( const pair<size_t,double> &cc : contribs )
        row[cc.first] += (sum > 0.0) ? (cc.second * (c.expectedCounts / sum)) : 0.0;
      m_rows.push_back( row );

      const double sigma2 = c.observedUncert*c.observedUncert
                            + std::pow( sm_floor_frac*c.expectedCounts, 2 );
      prof.energy.push_back( c.energy );
      prof.w.push_back( 1.0 / std::sqrt( sigma2 ) );
      prof.y.push_back( c.observedCounts / std::sqrt( sigma2 ) );
      const double atten = d.m_totalShieldAttenFactor;
      prof.tau.push_back( (atten > 0.0) ? std::max( 0.0, -std::log( std::min( 1.0, atten ) ) ) : 745.0 );
    }//for( peaks )

    // Activity columns only for nuclides that contribute somewhere.
    vector<size_t> cols;
    for( size_t k = 0; k < prof_nucs.size(); ++k )
    {
      double tot = 0.0;
      for( const vector<double> &row : m_rows )
        tot += (k < row.size()) ? row[k] : 0.0;
      if( tot > 0.0 )
        cols.push_back( k );
    }
    prof.n = m_rows.size();
    prof.K = cols.size();
    prof.logm.assign( prof.n*prof.K, -std::numeric_limits<double>::infinity() );
    for( size_t i = 0; i < prof.n; ++i )
      for( size_t j = 0; j < prof.K; ++j )
      {
        const double m = (cols[j] < m_rows[i].size()) ? m_rows[i][cols[j]] : 0.0;
        if( m > 0.0 )
          prof.logm[i*prof.K + j] = std::log( m );
      }
  }

  // ---- Model shielding summary ----
  double z0 = 0.0, ad0_gcm2 = 0.0;
  size_t nshields = 0;
  bool amount_at_bound = false;
  for( const ShieldingSourceFitCalc::ShieldingInfo &s : shieldings )
  {
    const double ad_gcm2 = along_ray_areal_density( s ) / (PhysicalUnits::g / PhysicalUnits::cm2);
    const bool amount_fit = s.m_isGenericMaterial
                              ? s.m_fitDimensions[1]
                              : (s.m_fitDimensions[0] || s.m_fitDimensions[1] || s.m_fitDimensions[2]);
    if( ad_gcm2 > 1.0E-6 )
      nshields += 1;
    if( amount_fit && (ad_gcm2 < 0.01) )
      amount_at_bound = true;   //e.g. a bare source with a thickness fit that went to zero
    ad0_gcm2 += ad_gcm2;
  }
  {
    double eff_z = 0.0;
    if( effective_atomic_number( shieldings, eff_z ) )
      z0 = eff_z;
  }
  result.modelHasShielding = (nshields > 0);

  const bool amount_was_fit = (config.physicalDimensionFit || config.genericAdFit) && !amount_at_bound;
  const bool an_was_fit = config.genericAnFit;
  const bool do_amount = !amount_was_fit;
  // Physical (amount / atomic-number) calls are not made when an age was fit: age changes line
  //  ratios in ways the shielding alternatives do not model.
  const bool physical_ok = !config.anyAgeFit;
  const bool do_z = !an_was_fit && (nshields == 1) && (z0 > 0.0);

  // ---- Profile likelihood-ratio tests ----
  result.numPeaksTested = prof.n;
  result.residualScale = 0.0;
  ProfileResult pr;
  if( !config.hasSelfAttenOrTraceSrc && (prof.n >= prof.K + 2) && (prof.K >= 1) )
  {
    try
    {
      pr = run_profile( prof, result.modelHasShielding, do_z, z0, ad0_gcm2 );
      if( pr.ok )
      {
        result.testsRun = true;
        const int nalt = std::max( 1, pr.p_best ) + 2;  //+2: the "other" polynomial
        const int dof = static_cast<int>( prof.n ) - static_cast<int>( prof.K ) - nalt;
        const double s2 = std::max( sm_min_residual_scale*sm_min_residual_scale,
                                    (dof > 0) ? (pr.chi2_other / dof) : 1.0 );
        result.residualScale = std::sqrt( s2 );

        const double d_amount = std::max( 0.0, pr.chi2_0 - pr.chi2_a );
        const double d_z = std::max( 0.0, pr.chi2_a - pr.chi2_b );
        const double d_other = std::max( 0.0, pr.chi2_best - pr.chi2_other );
        const double amount_sign = (pr.added_layer || (pr.u_best >= 1.0)) ? 1.0 : -1.0;
        result.amountT = (d_amount > 0.0) ? (amount_sign * std::sqrt( d_amount / s2 )) : 0.0;
        result.atomicNumberT = (d_z > 0.0) ? (((pr.z_best >= z0) ? 1.0 : -1.0) * std::sqrt( d_z / s2 )) : 0.0;
        result.otherT = z_from_chi2_2dof( d_other / s2 );

        // Atomic-number power: the atomicNumberT a 1.5x Z error (either way) would give, from
        //  synthetic data with the reference Z and AD matched at the most significant peak.  Only
        //  where an atomic-number call is possible, so "could resolve Z" is never reported for a fit
        //  whose Z cannot be judged.
        if( do_z && physical_ok && (ad0_gcm2 > 0.0) )
        {
          size_t jbest = 0;
          double best_signif = -1.0;
          vector<double> model_counts( prof.n, 0.0 );
          for( size_t i = 0; i < prof.n; ++i )
          {
            for( size_t j = 0; j < prof.K; ++j )
              if( std::isfinite( prof.logm[i*prof.K + j] ) )
                model_counts[i] += std::exp( prof.logm[i*prof.K + j] );
            if( model_counts[i]*prof.w[i] > best_signif )
            {
              best_signif = model_counts[i]*prof.w[i];
              jbest = i;
            }
          }
          const vector<double> mu0 = prof.mu( z0 );
          double power = 0.0;
          for( const double zr : { z0*sm_an_reference_z_factor, z0/sm_an_reference_z_factor } )
          {
            const double zref = std::min( 94.0, std::max( 3.0, zr ) );
            const vector<double> mur = prof.mu( zref );
            const double corr_j = (mu0[jbest]*ad0_gcm2 > 0.0) ? (prof.tau[jbest]/(mu0[jbest]*ad0_gcm2)) : 1.0;
            if( !(mur[jbest]*corr_j > 0.0) )
              continue;
            const double ad_ref = prof.tau[jbest] / (mur[jbest]*corr_j);
            ProfileData alt = prof;
            for( size_t i = 0; i < prof.n; ++i )
            {
              const double corr_i = (mu0[i]*ad0_gcm2 > 0.0) ? (prof.tau[i]/(mu0[i]*ad0_gcm2)) : 1.0;
              const double lc = (model_counts[i] > 0.0) ? std::log( model_counts[i] ) : -1000.0;
              alt.y[i] = prof.w[i] * std::exp( std::min( 700.0, lc - (ad_ref*mur[i]*corr_i - prof.tau[i]) ) );
            }
            const ProfileResult pa = run_profile( alt, true, true, z0, ad0_gcm2 );
            if( pa.ok )
              power = std::max( power, std::sqrt( std::max( 0.0, pa.chi2_a - pa.chi2_b ) / s2 ) );
          }//for( both directions )
          result.anDiscriminationT = power;
          result.anDiscriminable = (power >= sm_an_power_gate);
        }//if( compute power )
      }//if( pr.ok )
    }catch( std::exception & )
    {
      result.testsRun = false;
    }
  }//if( enough peaks to test )

  // ---- Decision ----
  //  AN and AD both fit: only "other" is possible.  Otherwise an amount problem (when the amount
  //  was not fit) is reported first - it is the most actionable ("fit the thickness") - then an
  //  atomic-number problem (soft calls only when the data have the power), then "other".
  if( result.testsRun )
  {
    if( an_was_fit && config.genericAdFit )
    {
      // Neither the amount nor the composition can be blamed.
    }else if( physical_ok && do_amount && (std::fabs(result.amountT) >= sm_call_t) )
    {
      if( pr.added_layer )
        result.conclusion = TrendConclusion::ModelLacksShielding;
      else
        result.conclusion = (result.amountT > 0.0) ? TrendConclusion::ShieldingAmountTooLow
                                                   : TrendConclusion::ShieldingAmountTooHigh;
      result.confidence = confidence_from_t( std::fabs(result.amountT) );
    }else if( physical_ok && do_z
              && (std::fabs(result.atomicNumberT) >= (result.anDiscriminable ? sm_soft_z_t : sm_call_t)) )
    {
      result.conclusion = (result.atomicNumberT > 0.0) ? TrendConclusion::AtomicNumberTooLow
                                                       : TrendConclusion::AtomicNumberTooHigh;
      result.confidence = confidence_from_t( std::fabs(result.atomicNumberT) );
    }

    if( (result.conclusion == TrendConclusion::NoConclusion) && (result.otherT >= sm_other_t) )
    {
      result.conclusion = TrendConclusion::OtherModelInconsistency;
      result.confidence = confidence_from_t( result.otherT );
    }
  }//if( tests ran )

  // ---- Hints: from the alternative that matches the conclusion ----
  const double g_cm2 = PhysicalUnits::g / PhysicalUnits::cm2;
  switch( result.conclusion )
  {
    case TrendConclusion::ModelLacksShielding:
      result.suggestedAtomicNumberResolved = pr.add_z_resolved;
      result.suggestedAtomicNumber = pr.add_z;
      result.suggestedArealDensity = (pr.add_z_resolved ? pr.add_ad : pr.add_fe_ad) * g_cm2;
      break;

    case TrendConclusion::ShieldingAmountTooLow:
    case TrendConclusion::ShieldingAmountTooHigh:
      result.suggestedArealDensity = pr.u_best * ad0_gcm2 * g_cm2;
      break;

    case TrendConclusion::AtomicNumberTooLow:
    case TrendConclusion::AtomicNumberTooHigh:
      result.suggestedAtomicNumber = pr.z_best;
      result.suggestedArealDensity = pr.ad_best * g_cm2;
      result.suggestedAtomicNumberResolved = true;
      break;

    case TrendConclusion::NoConclusion:
    case TrendConclusion::OtherModelInconsistency:
      break;
  }//switch( result.conclusion )


  // ---- Descriptive polynomial curve (pulls vs g(E)), for the chi chart only ----
  vector<string> nuc_names;   //order of first appearance
  vector<string> nuc_colors;
  vector<PullPoint> points;

  for( size_t i = 0; i < npoints; ++i )
  {
    const GammaInteractionCalc::PeakResultPlotInfo &comparison = peak_comparisons[i];
    const GammaInteractionCalc::PeakDetail &detail = peak_details[i];

    if( detail.assignedNuclide.empty() )
      continue;

    if( !std::isfinite(comparison.energy) || !std::isfinite(comparison.numSigmaOff)
        || (comparison.energy <= 0.0) )
      continue;

    // Ambiguity guard: require the assigned nuclide to dominate the peaks model counts.
    double total_at_src = 0.0, assigned_at_src = 0.0;
    for( const GammaInteractionCalc::PeakDetailSrc &src : detail.m_sources )
    {
      total_at_src += src.countsAtSource;
      if( src.nuclide && (src.nuclide->symbol == detail.assignedNuclide) )
        assigned_at_src += src.countsAtSource;
    }

    if( (total_at_src > 0.0) && (assigned_at_src < sm_min_nuclide_fraction*total_at_src) )
      continue;

    size_t nuc_index = nuc_names.size();
    for( size_t j = 0; j < nuc_names.size(); ++j )
    {
      if( nuc_names[j] == detail.assignedNuclide )
      {
        nuc_index = j;
        break;
      }
    }

    if( nuc_index == nuc_names.size() )
    {
      nuc_names.push_back( detail.assignedNuclide );
      nuc_colors.push_back( "" );
    }

    if( nuc_colors[nuc_index].empty() && !comparison.peakColor.isDefault() )
      nuc_colors[nuc_index] = comparison.peakColor.cssText( false );

    PullPoint point;
    point.energy = comparison.energy;
    point.pull = std::max( -sm_max_abs_pull, std::min( sm_max_abs_pull, comparison.numSigmaOff ) );
    point.nucIndex = nuc_index;
    points.push_back( point );
  }//for( loop over peaks )

  // Nuclides with a single point cant constrain the shared shape (their intercept just
  //  zeroes their residual), so drop them from the fit.
  {
    vector<size_t> nuc_counts( nuc_names.size(), 0 );
    for( const PullPoint &point : points )
      nuc_counts[point.nucIndex] += 1;

    vector<size_t> new_index( nuc_names.size(), 0 );
    vector<string> kept_names, kept_colors;
    for( size_t j = 0; j < nuc_names.size(); ++j )
    {
      new_index[j] = kept_names.size();
      if( nuc_counts[j] >= 2 )
      {
        kept_names.push_back( nuc_names[j] );
        kept_colors.push_back( nuc_colors[j] );
      }
    }

    vector<PullPoint> kept_points;
    for( const PullPoint &point : points )
    {
      if( nuc_counts[point.nucIndex] >= 2 )
      {
        PullPoint kept = point;
        kept.nucIndex = new_index[point.nucIndex];
        kept_points.push_back( kept );
      }
    }

    nuc_names.swap( kept_names );
    nuc_colors.swap( kept_colors );
    points.swap( kept_points );
  }

  const size_t num_nucs = nuc_names.size();
  result.numPointsUsed = points.size();

  if( num_nucs < 1 )
    return result;

  // ---- Regressor: g(E), centered and scaled to [-1,1] over the used points ----
  TrendBasis used_basis = basis;
  vector<double> gvals( points.size(), 0.0 );

  for( int attempt = 0; attempt < 2; ++attempt )
  {
    bool all_finite = true;
    for( size_t i = 0; i < points.size(); ++i )
    {
      gvals[i] = basis_transform( used_basis, points[i].energy, shieldings );
      all_finite = (all_finite && std::isfinite(gvals[i]));
    }

    if( all_finite )
      break;

    if( used_basis == TrendBasis::LogEnergy )
      return result; //energies already validated >0, shouldnt happen
    used_basis = TrendBasis::LogEnergy;
  }//for( basis fallback attempt )

  result.basis = used_basis;

  const auto minmax_g = std::minmax_element( begin(gvals), end(gvals) );
  const double g_mid = 0.5*((*minmax_g.first) + (*minmax_g.second));
  const double g_half_range = 0.5*((*minmax_g.second) - (*minmax_g.first));

  if( (g_half_range <= 0.0) || !std::isfinite(g_half_range) )
    return result;

  const auto scaled_x = [g_mid,g_half_range]( const double g ) -> double {
    return (g - g_mid) / g_half_range;
  };

  // Joint LLS: per-nuclide intercept columns + shared x (+ shared x^2); coefficient uncertainties
  //  scaled by the residual RMS sqrt(max(1, chi2/dof)).  Descriptive only.
  struct Fit
  {
    bool ok = false;
    std::vector<double> intercept, interceptUncert;
    double slope = 0.0, slopeUncert = 0.0, slopeT = 0.0;
    double curv = 0.0, curvUncert = 0.0, curvT = 0.0;
    double chi2 = 0.0, redChi2 = 0.0;
    unsigned int dof = 0;
  };//struct Fit

  const auto solve_trend = [&]( const bool quadratic ) -> Fit {
    Fit f;
    const size_t nshape = quadratic ? 2 : 1;
    const size_t ncols = num_nucs + nshape;
    const size_t nrows = points.size();
    if( nrows < (ncols + 1) )  //require at least 1 DOF
      return f;

    Eigen::MatrixXd design = Eigen::MatrixXd::Zero( nrows, ncols );
    Eigen::VectorXd rhs( nrows );
    for( size_t i = 0; i < nrows; ++i )
    {
      const double x = scaled_x( gvals[i] );
      design(i, points[i].nucIndex) = 1.0;
      design(i, num_nucs) = x;
      if( quadratic )
        design(i, num_nucs + 1) = x*x;
      rhs(i) = points[i].pull;
    }

    const Eigen::LDLT<Eigen::MatrixXd> ldlt( design.transpose() * design );
    if( ldlt.info() != Eigen::Success )
      return f;
    const Eigen::VectorXd solution = ldlt.solve( design.transpose() * rhs );
    const Eigen::MatrixXd covariance = ldlt.solve( Eigen::MatrixXd::Identity(ncols,ncols) );
    for( size_t i = 0; i < ncols; ++i )
    {
      if( !std::isfinite(solution(i)) || !std::isfinite(covariance(i,i)) || (covariance(i,i) <= 0.0) )
        return f;
    }

    f.chi2 = (rhs - design*solution).squaredNorm();
    f.dof = static_cast<unsigned int>( nrows - ncols );
    f.redChi2 = (f.dof > 0) ? (f.chi2 / f.dof) : 0.0;
    const double err_scale = std::sqrt( std::max( 1.0, f.redChi2 ) );

    f.intercept.resize( num_nucs );
    f.interceptUncert.resize( num_nucs );
    for( size_t j = 0; j < num_nucs; ++j )
    {
      f.intercept[j] = solution(j);
      f.interceptUncert[j] = err_scale * std::sqrt( covariance(j,j) );
    }
    f.slope = solution(num_nucs);
    f.slopeUncert = err_scale * std::sqrt( covariance(num_nucs,num_nucs) );
    f.slopeT = (f.slopeUncert > 0.0) ? (f.slope/f.slopeUncert) : 0.0;
    if( quadratic )
    {
      f.curv = solution(num_nucs+1);
      f.curvUncert = err_scale * std::sqrt( covariance(num_nucs+1,num_nucs+1) );
      f.curvT = (f.curvUncert > 0.0) ? (f.curv/f.curvUncert) : 0.0;
    }
    f.ok = true;
    return f;
  };//solve_trend lambda

  const Fit quad = solve_trend( true );
  const Fit lin  = solve_trend( false );
  if( !quad.ok && !lin.ok )
    return result;

  const auto use_quadratic = [&](){
    result.model = TrendModel::Quadratic;
    result.chi2 = quad.chi2; result.dof = quad.dof; result.redChi2 = quad.redChi2;
    result.slope = quad.slope; result.slopeUncert = quad.slopeUncert; result.slopeT = quad.slopeT;
    result.curvature = quad.curv; result.curvatureUncert = quad.curvUncert; result.curvatureT = quad.curvT;
  };
  const auto use_linear = [&](){
    result.model = TrendModel::Linear;
    result.chi2 = lin.chi2; result.dof = lin.dof; result.redChi2 = lin.redChi2;
    result.slope = lin.slope; result.slopeUncert = lin.slopeUncert; result.slopeT = lin.slopeT;
    result.curvature = 0.0; result.curvatureUncert = 0.0; result.curvatureT = 0.0;
  };

  // Draw a curve whenever a conclusion was reached (quadratic, which follows attenuation shapes
  //  best), else only when the polynomial itself is firm.
  const bool firm_curv  = quad.ok && (std::fabs(quad.curvT) >= sm_curve_t_threshold);
  const bool firm_slope = lin.ok  && (std::fabs(lin.slopeT) >= sm_curve_t_threshold);
  if( result.conclusion != TrendConclusion::NoConclusion )
  {
    if( quad.ok ) use_quadratic(); else use_linear();
  }else if( firm_curv )
  {
    use_quadratic();
  }else if( firm_slope )
  {
    use_linear();
  }else
  {
    // No curve: keep the quadratic's scatter for transparency, but no trend coefficients.
    const Fit &f = quad.ok ? quad : lin;
    result.chi2 = f.chi2; result.dof = f.dof; result.redChi2 = f.redChi2;
  }

  // ---- Per-nuclide curves (from the chosen model) ----
  const std::vector<double> &draw_intercept = (result.model == TrendModel::Quadratic)
                                                ? quad.intercept : lin.intercept;
  const std::vector<double> &draw_intercept_unc = (result.model == TrendModel::Quadratic)
                                                ? quad.interceptUncert : lin.interceptUncert;
  for( size_t j = 0; j < num_nucs; ++j )
  {
    NuclideTrend trend;
    trend.nuclide = nuc_names[j];
    trend.cssColor = nuc_colors[j];
    trend.intercept = (j < draw_intercept.size()) ? draw_intercept[j] : 0.0;
    trend.interceptUncert = (j < draw_intercept_unc.size()) ? draw_intercept_unc[j] : 0.0;
    trend.minEnergy = std::numeric_limits<double>::max();
    trend.maxEnergy = 0.0;

    for( const PullPoint &point : points )
    {
      if( point.nucIndex != j )
        continue;
      trend.numPeaksUsed += 1;
      trend.minEnergy = std::min( trend.minEnergy, point.energy );
      trend.maxEnergy = std::max( trend.maxEnergy, point.energy );
    }

    if( (result.model != TrendModel::None) && (trend.numPeaksUsed >= 2) )
    {
      for( size_t k = 0; k < sm_num_curve_samples; ++k )
      {
        const double frac = static_cast<double>(k) / (sm_num_curve_samples - 1);
        const double energy = trend.minEnergy + frac*(trend.maxEnergy - trend.minEnergy);
        const double g = basis_transform( result.basis, energy, shieldings );
        const double x = scaled_x( g );
        const double chi = trend.intercept + result.slope*x
                           + ((result.model == TrendModel::Quadratic) ? result.curvature*x*x : 0.0);
        if( std::isfinite(chi) )
          trend.curve.emplace_back( energy, chi );
      }
    }//if( draw a curve for this nuclide )

    result.nuclideTrends.push_back( trend );
  }//for( loop over nuclides )

  return result;
}//fit_pull_trend(...)


const char *conclusion_txt_key( const TrendConclusion conclusion )
{
  switch( conclusion )
  {
    case TrendConclusion::NoConclusion:            return nullptr;
    case TrendConclusion::ShieldingAmountTooLow:   return "x2g-trend-shield-too-low";
    case TrendConclusion::ShieldingAmountTooHigh:  return "x2g-trend-shield-too-high";
    case TrendConclusion::AtomicNumberTooLow:      return "x2g-trend-an-too-low";
    case TrendConclusion::AtomicNumberTooHigh:     return "x2g-trend-an-too-high";
    case TrendConclusion::OtherModelInconsistency: return "x2g-trend-other-inconsistency";
    case TrendConclusion::ModelLacksShielding:     return "x2g-trend-lacks-shielding";
  }//switch( conclusion )

  assert( 0 );
  return nullptr;
}//conclusion_txt_key(...)


const char *confidence_txt_key( const TrendConfidence confidence )
{
  switch( confidence )
  {
    case TrendConfidence::Tentative: return "x2g-trend-conf-tentative";
    case TrendConfidence::Weak:      return "x2g-trend-conf-weak";
    case TrendConfidence::Moderate:  return "x2g-trend-conf-moderate";
    case TrendConfidence::Strong:    return "x2g-trend-conf-strong";
  }//switch( confidence )

  assert( 0 );
  return "x2g-trend-conf-weak";
}//confidence_txt_key(...)


std::string areal_density_txt( const double areal_density )
{
  const double gcm2 = areal_density / (PhysicalUnits::g / PhysicalUnits::cm2);
  std::ostringstream strm;
  if( gcm2 >= 10.0 )
    strm << std::fixed << std::setprecision(0) << gcm2;
  else if( gcm2 >= 1.0 )
    strm << std::fixed << std::setprecision(1) << gcm2;
  else
    strm << std::setprecision(1) << gcm2;
  strm << " g/cm\u00B2";
  return strm.str();
}//areal_density_txt(...)


namespace
{
  /** Z-class key for the reference layers of the "lacks shielding" hint. */
  const char *zclass_txt_key( const double z )
  {
    if( z <= 8.0 )  return "x2g-trend-zclass-water";
    if( z <= 16.0 ) return "x2g-trend-zclass-al";
    if( z <= 35.0 ) return "x2g-trend-zclass-fe";
    if( z <= 65.0 ) return "x2g-trend-zclass-sn";
    return "x2g-trend-zclass-pb";
  }

  const char *zclass_english( const double z )
  {
    if( z <= 8.0 )  return "a very low-Z (water- or plastic-like)";
    if( z <= 16.0 ) return "a low-Z (aluminum-like)";
    if( z <= 35.0 ) return "a medium-Z (iron-like)";
    if( z <= 65.0 ) return "a high-Z (tin-like)";
    return "a very high-Z (lead-like)";
  }
}//namespace


MessageHint message_hint( const TrendResult &trend )
{
  MessageHint hint;
  switch( trend.conclusion )
  {
    case TrendConclusion::ShieldingAmountTooLow:
    case TrendConclusion::ShieldingAmountTooHigh:
      if( (trend.conclusion == TrendConclusion::ShieldingAmountTooHigh)
          && (trend.suggestedArealDensity < sm_little_or_no_shielding) )
      {
        hint.key = "x2g-trend-hint-no-shield";   //e.g., a bare source modeled behind shielding
      }else if( trend.suggestedArealDensity > 0.0 )
      {
        hint.key = "x2g-trend-hint-amount";
        hint.value_txt = areal_density_txt( trend.suggestedArealDensity );
      }
      break;

    case TrendConclusion::ModelLacksShielding:
      if( trend.suggestedArealDensity > 0.0 )
      {
        hint.value_txt = areal_density_txt( trend.suggestedArealDensity );
        if( trend.suggestedAtomicNumberResolved )
        {
          hint.key = "x2g-trend-hint-lacks-z";
          hint.zclass_key = zclass_txt_key( trend.suggestedAtomicNumber );
        }else
        {
          hint.key = "x2g-trend-hint-lacks-feeq";
        }
      }
      break;

    case TrendConclusion::AtomicNumberTooLow:
    case TrendConclusion::AtomicNumberTooHigh:
      if( trend.suggestedAtomicNumber > 0.0 )
      {
        hint.key = "x2g-trend-hint-z";
        hint.value_txt = std::to_string( static_cast<int>( std::lround( trend.suggestedAtomicNumber ) ) );
      }
      break;

    case TrendConclusion::NoConclusion:
    case TrendConclusion::OtherModelInconsistency:
      break;
  }//switch( trend.conclusion )

  return hint;
}//message_hint(...)


const char *summing_caveat_txt_key()
{
  return "x2g-trend-summing-caveat";
}


std::string conclusion_report_text( const TrendResult &trend )
{
  if( trend.conclusion == TrendConclusion::NoConclusion )
    return std::string();

  const char * const caveat = trend.possibleCascadeSumming
      ? " The source is close to the detector and cascade-summing correction is off, so summing may"
        " also contribute."
      : "";

  if( trend.conclusion == TrendConclusion::OtherModelInconsistency )
    return std::string( "A residual energy trend remains that changing the shielding does not explain;"
                        " the distance, geometry, source age, or detector response may be off." ) + caveat;

  const char *qualifier = "";
  switch( trend.confidence )
  {
    case TrendConfidence::Tentative: qualifier = "possibly indicate";  break;
    case TrendConfidence::Weak:      qualifier = "may indicate";       break;
    case TrendConfidence::Moderate:  qualifier = "probably indicate";  break;
    case TrendConfidence::Strong:    qualifier = "strongly indicate";  break;
  }

  const char *body = nullptr;
  switch( trend.conclusion )
  {
    case TrendConclusion::NoConclusion:            return std::string();
    case TrendConclusion::ShieldingAmountTooLow:   body = "too little modeled shielding"; break;
    case TrendConclusion::ShieldingAmountTooHigh:  body = "too much modeled shielding";   break;
    case TrendConclusion::AtomicNumberTooLow:      body = "the shielding's effective atomic number is too low";  break;
    case TrendConclusion::AtomicNumberTooHigh:     body = "the shielding's effective atomic number is too high"; break;
    case TrendConclusion::ModelLacksShielding:     body = "the model lacks shielding"; break;
    case TrendConclusion::OtherModelInconsistency: return std::string(); //handled above
  }

  std::string hint;
  const MessageHint h = message_hint( trend );
  if( h.key )
  {
    string ad = h.value_txt;
    const string superscript_two = "\u00B2";
    const size_t pos = ad.find( superscript_two );
    if( pos != string::npos )
      ad.replace( pos, superscript_two.size(), "2" );  //plain "g/cm2" in non-HTML reports
    if( string(h.key) == "x2g-trend-hint-no-shield" )
      hint = "; the peaks fit best with little or no shielding";
    else if( string(h.key) == "x2g-trend-hint-amount" )
      hint = "; about " + ad + " in total would fit best";
    else if( string(h.key) == "x2g-trend-hint-lacks-z" )
      hint = "; the peaks are consistent with roughly " + ad + " of " + zclass_english( trend.suggestedAtomicNumber ) + " material";
    else if( string(h.key) == "x2g-trend-hint-lacks-feeq" )
      hint = "; the peaks are consistent with roughly " + ad + " of iron-equivalent shielding";
    else if( string(h.key) == "x2g-trend-hint-z" )
      hint = "; the peaks favor an effective atomic number near " + h.value_txt;
  }

  return "Peak residuals " + std::string(qualifier) + " " + std::string(body) + hint + "." + caveat;
}//conclusion_report_text(...)

}//namespace ShieldSourcePullTrend
