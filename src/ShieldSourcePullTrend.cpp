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
#include <string>
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>

#include <Eigen/Dense>

#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldSourcePullTrend.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"

using namespace std;

namespace
{
  /** Thresholds and constants of the trend analysis.

   Calibrated by scratch/chi_trend_study (see RESULTS.md, 2026-07).  The reported t-statistics
   already divide the coefficient uncertainties by the residual RMS sqrt(max(1,chi2/dof)) (see
   solve_trend), which restores an approximately-N(0,1) null distribution even with realistic
   ~2x-Poisson peak-fit scatter: the calibrated null percentiles are p95~2.1, p99~3.8.  A
   threshold of 3.5 yields ~1.5% false-positive rate on correct models - a good operating point
   for a passive advisory that should err toward silence.  (Before the residual-variance scaling
   the null |t| tail ran to ~11 and forced a much higher, less powerful threshold.)
   */

  /** |t| of the quadratic coefficient required to keep the quadratic term. */
  const double sm_curvature_t_threshold = 3.5;

  /** |t| of the slope required to report any linear trend. */
  const double sm_slope_t_threshold = 3.5;

  /** Soft-advisory floor for atomic-number calls: when the data demonstrably has the power to
   resolve a Z error (see the discrimination-power check), even a weak observed curvature down to
   this |t| is worth a hedged "possibly" advisory.  Below the firm 3.5 threshold, so it only ever
   fires for atomic-number calls that pass the power gate - never for the amount diagnosis, which
   has no such physical detectability handle.  The null false-positive rate at |t|=1.5 is ~13%,
   which is why the wording is explicitly tentative.
   */
  const double sm_an_advisory_t = 1.5;

  /** |t| tier boundaries mapping to TrendConfidence::{Tentative,Weak,Moderate,Strong}.  Strong is
   near/above the null 99th percentile (~3.8) times a margin so "strongly" is rarely spurious.
   */
  const double sm_weak_t_threshold = 2.5;
  const double sm_moderate_t_threshold = 4.0;
  const double sm_strong_t_threshold = 6.0;

  /** Pulls beyond this many sigma are winsorized (clamped) so a single catastrophically
   mis-fit peak cannot dominate the unweighted least-squares trend.
   */
  const double sm_max_abs_pull = 10.0;

  /** Minimum fraction of a peaks model counts that must come from its assigned nuclide for
   the point to be used in the trend fit (guards against ambiguous multi-nuclide peaks).
   */
  const double sm_min_nuclide_fraction = 0.7;

  /** The reference effective-atomic-number error used for the discrimination-power check: the
   power is the curvature t-statistic that a Z error of this factor would imprint on the current
   peak set.  1.5 is "an error worth flagging" - e.g. mistaking steel (~26) for tin (~40).
   */
  const double sm_an_reference_z_factor = 1.5;

  /** Discrimination-power GATE: an atomic-number call is only made when anDiscriminationT reaches
   this value.  Calibrated by scratch/chi_trend_study --powercal (RESULTS.md): across realistic
   scenarios carrying a known 1.5x Z error, detection of that error becomes reliable once the
   predicted power reaches ~2.5 (caught >60% per power-bin, ~80% among gate-passers), and falls
   off below it.  Deliberately below the 3.5 curve-drawing threshold so the soft "possibly"
   advisory can fire on discriminable data.
   */
  const double sm_an_power_gate = 2.5;

  /** Number of samples of the fitted curve returned per-nuclide, for plotting. */
  const size_t sm_num_curve_samples = 40;


  /** A pull point accepted into the trend fit. */
  struct PullPoint
  {
    double energy;      //keV
    double pull;        //winsorized (observed-model)/uncert
    double significance;//observed counts / uncertainty (for the discrimination-power check)
    double tauModel;    //model optical depth -ln(total shield attenuation fraction) (dimensionless)
    size_t nucIndex;    //index into the nuclide grouping
  };//struct PullPoint


  /** Expected curvature t-statistic that a reference atomic-number error (Z_model -> Z_ref) would
   imprint on the given peaks, with the shielding amount re-fit.  Inputs per peak: energy,
   significance (area/uncert), and the model's actual optical depth tau_model = -ln(attenuation
   fraction) - a DIMENSIONLESS quantity taken from the fit, so this is robust to areal-density
   unit conventions and needs no separate areal density.

   Physics: the Z error scales each peak's optical depth by the ratio mu(Z_ref,E)/mu(Z_model,E)
   (the areal density cancels), so tau_ref(E) = tau_model(E) * mu(Z_ref,E)/mu(Z_model,E).  With the
   amount re-fit, the induced pull is significance(E) * [residual of tau_ref projected out of
   {1, tau_model}].  A pure design property (no observed pulls).  Returns 0 if not evaluable.
   */
  double expected_an_curvature_t( const std::vector<double> &energies,
                                  const std::vector<double> &signif,
                                  const std::vector<double> &tau_model,
                                  const double z_model_in, const double z_ref_in )
  {
    const size_t n = energies.size();
    if( (n < 4) || (z_model_in <= 0.0) || (z_ref_in <= 0.0) )
      return 0.0;

    // Need at least 3 distinct energies or the quadratic design {1,x,x^2} is rank-deficient.
    {
      std::vector<double> uniq( energies );
      std::sort( begin(uniq), end(uniq) );
      uniq.erase( std::unique( begin(uniq), end(uniq) ), end(uniq) );
      if( uniq.size() < 3 )
        return 0.0;
    }

    // The mass-attenuation table is only valid for Z in [1,98]; clamp so a heavy shield (Pb, W, U)
    //  whose reference Z would exceed 98 still gets a computed power instead of throwing (which
    //  would silently disable the atomic-number diagnosis for the most common heavy shields).
    const double z_model = std::min( 98.0, std::max( 1.0, z_model_in ) );
    const double z_ref = std::min( 98.0, std::max( 1.0, z_ref_in ) );

    Eigen::MatrixXd A( n, 2 );  //columns {1, tau_model}
    Eigen::VectorXd tau_ref( n ), w( n );
    double sum_tau = 0.0;
    for( size_t i = 0; i < n; ++i )
    {
      const double mu0 = GammaInteractionCalc::mass_attenuation_coef( static_cast<float>(z_model), static_cast<float>(energies[i]) );
      const double mu1 = GammaInteractionCalc::mass_attenuation_coef( static_cast<float>(z_ref),   static_cast<float>(energies[i]) );
      if( !std::isfinite(mu0) || !std::isfinite(mu1) || (mu0 <= 0.0) )
        return 0.0;
      A(i,0) = 1.0; A(i,1) = tau_model[i];
      tau_ref(i) = tau_model[i] * (mu1/mu0);
      w(i) = signif[i];
      sum_tau += tau_model[i];
    }
    if( sum_tau <= 0.0 )  //no attenuation in the model -> no Z signal possible
      return 0.0;

    const Eigen::MatrixXd Aw = w.asDiagonal() * A;
    const Eigen::VectorXd tw = w.asDiagonal() * tau_ref;
    const Eigen::LDLT<Eigen::MatrixXd> ldlt( Aw.transpose()*Aw );
    if( ldlt.info() != Eigen::Success )
      return 0.0;
    const Eigen::VectorXd ab = ldlt.solve( Aw.transpose()*tw );

    // Unweighted quadratic (in log-energy) of the induced pulls -> curvature t (design/unit variance).
    const double lmin = std::log(*std::min_element(begin(energies),end(energies)));
    const double lmax = std::log(*std::max_element(begin(energies),end(energies)));
    const double lmid = 0.5*(lmin+lmax), lhalf = 0.5*(lmax-lmin);
    if( lhalf <= 0.0 )
      return 0.0;

    Eigen::MatrixXd X( n, 3 );
    Eigen::VectorXd y( n );
    for( size_t i = 0; i < n; ++i )
    {
      const double xc = (std::log(energies[i]) - lmid)/lhalf;
      X(i,0) = 1.0; X(i,1) = xc; X(i,2) = xc*xc;
      y(i) = signif[i] * (tau_ref(i) - ab(0) - ab(1)*A(i,1));
    }
    const Eigen::MatrixXd cov = (X.transpose()*X).inverse();
    const Eigen::VectorXd b = cov * (X.transpose()*y);
    if( !std::isfinite(cov(2,2)) || (cov(2,2) <= 0.0) )
      return 0.0;
    return b(2) / std::sqrt( cov(2,2) );
  }//expected_an_curvature_t(...)


  /** Effective (areal-density-weighted) atomic number of the model shieldings, for the
   discrimination-power check.  Atomic number is unit-free, so for a single shielding this is
   robust to the AD unit conventions (the scale cancels in sum_z_ad/sum_ad).
   Caveat: when MIXING a generic shield (weight = areal density) with a material shield (weight =
   thickness*density) the two weights are in different units, and the fit_model vs interactive
   entry points store the generic AD in different units - so eff_z (which only selects a reference
   Z for the power gate) can differ slightly between the live preview and the final result for a
   multi-shield model.  The effect on the gate is minor; single-shield models are exact.
   Returns false if there is no usable shielding.
   */
  bool effective_atomic_number( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                double &eff_z )
  {
    double sum_ad = 0.0, sum_z_ad = 0.0;
    for( const ShieldingSourceFitCalc::ShieldingInfo &s : shieldings )
    {
      double ad = 0.0, z = 0.0;
      if( s.m_isGenericMaterial )
      {
        z = s.m_dimensions[0];
        ad = std::fabs( s.m_dimensions[1] );  //weight only; units irrelevant to the ratio
      }else if( s.m_material )
      {
        z = s.m_material->massWeightedAtomicNumber();
        ad = std::fabs( s.m_dimensions[0] * s.m_material->density );
      }else
      {
        continue;
      }
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
                                  static_cast<float>(shield.m_dimensions[0]),
                                  static_cast<float>(shield.m_dimensions[1]),
                                  static_cast<float>(energy) );
            }else if( shield.m_material )
            {
              optical_depth += GammaInteractionCalc::transmition_coefficient_material(
                                  shield.m_material.get(),
                                  static_cast<float>(energy),
                                  static_cast<float>(shield.m_dimensions[0]) );
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


  /** Maps |t| of the driving statistic to a confidence tier. */
  ShieldSourcePullTrend::TrendConfidence confidence_from_t( const double abs_t )
  {
    using ShieldSourcePullTrend::TrendConfidence;

    if( abs_t >= sm_strong_t_threshold )
      return TrendConfidence::Strong;
    if( abs_t >= sm_moderate_t_threshold )
      return TrendConfidence::Moderate;
    if( abs_t >= sm_weak_t_threshold )
      return TrendConfidence::Weak;
    return TrendConfidence::Tentative;
  }//confidence_from_t(...)
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
  }

  return config;
}//summarize_fit_config(...)


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

  // ---- Collect usable pull points, grouped by assigned nuclide ----
  const size_t npoints = std::min( peak_comparisons.size(), peak_details.size() );

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
    point.significance = (comparison.observedUncert > 0.0)
                            ? (comparison.observedCounts / comparison.observedUncert) : 0.0;
    // Model optical depth from the dimensionless total-shield attenuation fraction (1.0 = none).
    const double atten = detail.m_totalShieldAttenFactor;
    point.tauModel = ((atten > 0.0) && (atten < 1.0)) ? -std::log(atten) : 0.0;
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

  // ---- Joint LLS: per-nuclide intercept columns + shared x (+ shared x^2) ----
  //  Pulls are ~N(0,1) under a correct model, so the base covariance is (X^T X)^-1.  But real
  //  pulls are over-dispersed (2x-Poisson peak-fit systematics, and occasionally a poorly-
  //  converged shielding fit), so we scale the coefficient uncertainties by the residual RMS
  //  sqrt(max(1, chi2/dof)) - the standard unknown-noise-level correction.  This keeps the
  //  t-statistics calibrated (controlling false positives) and, crucially, suppresses garbage
  //  from poorly-converged fits: the calibration study (scratch/chi_trend_study, RESULTS.md)
  //  showed such fits otherwise flip the curvature sign and would break the atomic-number call.
  struct Fit
  {
    bool ok = false;
    std::vector<double> intercept, interceptUncert;
    double slope = 0.0, slopeUncert = 0.0, slopeT = 0.0;
    double curv = 0.0, curvUncert = 0.0, curvT = 0.0;
    double chi2 = 0.0, redChi2 = 0.0;
    unsigned int dof = 0;
  };//struct Fit

  // Joint LLS with residual-variance-scaled uncertainties (see comment above).  Returns a
  //  self-contained result so the quadratic and linear fits don't clobber each other.
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

  // ---- Atomic-number discrimination power (uses the quadratic residual scale) ----
  //  Whether these peaks could resolve a reference-sized effective-Z error at all - a design
  //  property of the data (energies + significances + model shielding), independent of the
  //  observed pulls.  Scaled by the same residual RMS so it is comparable to the observed t.
  try
  {
    double eff_z = 0.0;
    if( quad.ok && effective_atomic_number( shieldings, eff_z ) )
    {
      vector<double> pe, ps, ptau;
      for( const PullPoint &point : points )
      {
        pe.push_back( point.energy );
        ps.push_back( point.significance );
        ptau.push_back( point.tauModel );
      }
      // Probe a reference Z error in BOTH directions and take the larger response: for light/mid Z
      //  the up-shift (eff_z*factor) dominates (as in the calibration), while for heavy Z it
      //  saturates at the Z=98 table edge and the down-shift (eff_z/factor) carries the contrast.
      const double up_t = expected_an_curvature_t( pe, ps, ptau, eff_z, eff_z*sm_an_reference_z_factor );
      const double dn_t = expected_an_curvature_t( pe, ps, ptau, eff_z, eff_z/sm_an_reference_z_factor );
      const double design_t = (std::fabs(up_t) >= std::fabs(dn_t)) ? up_t : dn_t;
      const double err_scale = std::sqrt( std::max( 1.0, quad.redChi2 ) );
      result.anDiscriminationT = std::fabs( design_t ) / err_scale;
      result.anDiscriminable = (result.anDiscriminationT >= sm_an_power_gate);
    }
  }catch( std::exception & )
  {
    result.anDiscriminationT = 0.0;
    result.anDiscriminable = false;
  }

  // Default the reported statistics to the quadratic fit (a firm-model choice below overrides).
  //  This keeps redChi2/curvatureT populated for transparency even when no conclusion is drawn.
  if( quad.ok )
  {
    result.chi2 = quad.chi2; result.dof = quad.dof; result.redChi2 = quad.redChi2;
    result.slope = quad.slope; result.slopeUncert = quad.slopeUncert; result.slopeT = quad.slopeT;
    result.curvature = quad.curv; result.curvatureUncert = quad.curvUncert; result.curvatureT = quad.curvT;
  }else
  {
    result.chi2 = lin.chi2; result.dof = lin.dof; result.redChi2 = lin.redChi2;
    result.slope = lin.slope; result.slopeUncert = lin.slopeUncert; result.slopeT = lin.slopeT;
  }

  // ---- Conclusion + which model to draw ----
  //
  // Sign conventions, from the calibration study (scratch/chi_trend_study/RESULTS.md, 2026-07;
  //  the --zsweep controlled experiment made these unambiguous):
  //  - pull = (observed - model)/uncert.  A model with too LITTLE shielding (amount NOT fit)
  //    over-predicts the low-energy peaks -> POSITIVE slope vs log(energy) (mean slope t = +6.9
  //    for a 0.5x-shielding model; -6.4 for 2x).
  //  - A wrong effective atomic number (amount re-fit) shows up as CURVATURE whose SIGN gives the
  //    direction: model-Z-too-low => pulls concave-DOWN (negative curvature in log-energy),
  //    too-high => concave-UP (positive) - monotonic across nuclides/detectors, once the
  //    residual-variance scaling down-weights poorly-converged fits.
  //
  // Atomic-number calls are gated on the discrimination-power check (above): above ~250 keV all
  //  elements have the same attenuation SHAPE, so without low-energy leverage a Z error produces
  //  no curvature - a curvature that appears anyway cannot be a Z signature.  When the data DO
  //  have the power, we trust the observed curvature down to a soft "possibly" advisory (|t|>=1.5),
  //  since the physics guarantees a real effect would be visible; the confidence tier scales with
  //  the observed |t|.
  //
  //  amount was fit / generic-AD fit: amount is optimal, residual curvature -> effective Z.
  //  generic AN+AD both fit: attenuation freely shaped -> a residual points elsewhere (Other).
  //  activity only: slope -> amount direction; else curvature -> Z direction.

  const bool amount_was_fit = (config.physicalDimensionFit || config.genericAdFit);

  const bool firm_curv  = quad.ok && (std::fabs(quad.curvT) >= sm_curvature_t_threshold);
  const bool firm_slope = lin.ok  && (std::fabs(lin.slopeT) >= sm_slope_t_threshold);
  // Atomic-number advisory reachable (soft) only when the data can resolve Z.
  const bool an_advisory = quad.ok && result.anDiscriminable
                           && (std::fabs(quad.curvT) >= sm_an_advisory_t);

  // Choose the model to display/draw and populate result coefficients from that fit.
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

  const auto an_direction = [&]() -> TrendConclusion {
    return (quad.curv < 0.0) ? TrendConclusion::AtomicNumberTooLow
                             : TrendConclusion::AtomicNumberTooHigh;
  };

  // Decide conclusion and drawn model together.
  if( config.hasSelfAttenOrTraceSrc )
  {
    // No conclusion (physics differs), but still draw the strongest firm trend if present.
    if( firm_curv ) use_quadratic();
    else if( firm_slope ) use_linear();
  }else if( config.genericAnFit && config.genericAdFit )
  {
    if( firm_curv ) use_quadratic(); else if( firm_slope ) use_linear();
    if( firm_curv || firm_slope )
    {
      result.conclusion = TrendConclusion::OtherModelInconsistency;
      result.confidence = confidence_from_t( std::fabs( firm_curv ? quad.curvT : lin.slopeT ) );
    }
  }else if( amount_was_fit )
  {
    if( an_advisory )
    {
      use_quadratic();  //draw the curve even for a soft advisory
      result.conclusion = an_direction();
      result.confidence = confidence_from_t( std::fabs(quad.curvT) );
    }else if( firm_curv || firm_slope )
    {
      // A real trend, but not attributable to atomic number (no power) -> something else.
      if( firm_curv ) use_quadratic(); else use_linear();
      result.conclusion = TrendConclusion::OtherModelInconsistency;
      result.confidence = confidence_from_t( std::fabs( firm_curv ? quad.curvT : lin.slopeT ) );
    }
  }else
  {
    // Activity-only fit: a firm linear slope is the shielding amount; otherwise a Z curvature.
    if( firm_slope && !firm_curv )
    {
      use_linear();
      result.conclusion = (lin.slopeT > 0.0) ? TrendConclusion::ShieldingAmountTooLow
                                             : TrendConclusion::ShieldingAmountTooHigh;
      result.confidence = confidence_from_t( std::fabs(lin.slopeT) );
    }else if( an_advisory )
    {
      use_quadratic();
      result.conclusion = an_direction();
      result.confidence = confidence_from_t( std::fabs(quad.curvT) );
    }else if( firm_curv || firm_slope )
    {
      if( firm_curv ) use_quadratic(); else use_linear();
      result.conclusion = TrendConclusion::OtherModelInconsistency;
      result.confidence = confidence_from_t( std::fabs( firm_curv ? quad.curvT : lin.slopeT ) );
    }
  }

  // If no model was selected, the trend coefficients (defaulted to the quadratic fit above) are
  //  not meaningful - zero them so a reader never mistakes stale quad values for a real trend.
  if( result.model == TrendModel::None )
  {
    result.slope = result.slopeUncert = result.slopeT = 0.0;
    result.curvature = result.curvatureUncert = result.curvatureT = 0.0;
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


std::string conclusion_report_text( const TrendResult &trend )
{
  // The "other inconsistency" case does not use a confidence qualifier (matches the GUI message).
  if( trend.conclusion == TrendConclusion::OtherModelInconsistency )
    return "A residual energy trend remains despite fitting shielding; the distance, geometry,"
           " source age, or detector response may be slightly off.";

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
    case TrendConclusion::OtherModelInconsistency: return std::string(); //handled above
  }

  return "Peak residuals " + std::string(qualifier) + " " + std::string(body) + ".";
}//conclusion_report_text(...)

}//namespace ShieldSourcePullTrend
