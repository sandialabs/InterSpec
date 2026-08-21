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

// When enabled, observable peaks are computed on the original (non-energy-cal-adjusted)
// foreground with background subtraction, then the continuum is refit to the raw foreground.
// When disabled (0), observable peaks are computed on the fitted-cal foreground without
// background subtraction, then translated back to original energy cal.
#define OBSERVABLE_PEAKS_USING_ORIGINAL_CAL_WITH_BACK_SUB 0

#include <atomic>
#include <cstdint>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <array>
#include <tuple>
#include <vector>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"

#include "InterSpec/FitPeaksForNuclides.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"
#include "SpecUtils/SpecUtilsAsync.h"

// Note: PeakFitImprove::debug_printout is used in the original code but is in target/peak_fit_improve/
// For now, we'll conditionally compile debug output based on PERFORM_DEVELOPER_CHECKS
// This can be refined later if needed
#ifndef PERFORM_DEVELOPER_CHECKS
#define PERFORM_DEVELOPER_CHECKS 0
#endif

namespace
{
  // Local debug flag - can be set externally if needed
  // Note: In the original code, this was PeakFitImprove::debug_printout from target/peak_fit_improve/
  // For the moved code, we use a local flag that can be controlled
  bool local_debug_printout = false;
  
  // Helper to check if debug output should be printed
  // Replaces PeakFitImprove::debug_printout from the original code
  inline bool should_debug_print()
  {
#if PERFORM_DEVELOPER_CHECKS
    return local_debug_printout;
#else
    return false;
#endif
  }
  
  // Macro to replace PeakFitImprove::debug_printout throughout the code
  // This allows easy replacement of all instances
#define PEAK_FIT_DEBUG_PRINTOUT should_debug_print()
  
  // Helper functions that need access to should_use_desperation_shielding and create_desperation_shielding
  // These are used by fit_peaks_for_nuclide_relactauto but were in the anonymous namespace in the original
  bool should_use_desperation_shielding( const double atomic_number,
                                         const std::vector<RelActCalcAuto::RoiRange> &rois )
  {
    // Check 1: Valid atomic number
    if( atomic_number < 1.0 || atomic_number > 98.0 )
      return false;

    // Check 2: At least 2 ROIs required
    if( rois.size() < 2 )
      return false;

    // Check 3: Sufficient energy range
    // Find min lower_energy and max upper_energy across all ROIs
    double min_lower = std::numeric_limits<double>::max();
    double max_upper = std::numeric_limits<double>::lowest();
    bool has_low_energy_roi = false;

    for( const RelActCalcAuto::RoiRange &roi : rois )
    {
      min_lower = std::min( min_lower, roi.lower_energy );
      max_upper = std::max( max_upper, roi.upper_energy );

      if( roi.lower_energy < 120.0 )
        has_low_energy_roi = true;
    }

    const double energy_range = max_upper - min_lower;
    const double required_range = has_low_energy_roi ? 60.0 : 100.0;

    if( energy_range < required_range )
      return false;

    return true;
  }//should_use_desperation_shielding(...)


  std::shared_ptr<RelActCalc::PhysicalModelShieldInput>
  create_desperation_shielding( const double atomic_number, const double starting_areal_density )
  {
    if( atomic_number < 1.0 || atomic_number > 98.0 )
      throw std::invalid_argument( "create_desperation_shielding: atomic_number must be in range [1, 98]" );

    std::shared_ptr<RelActCalc::PhysicalModelShieldInput> shield
      = std::make_shared<RelActCalc::PhysicalModelShieldInput>();

    // Configure as generic element (not material)
    shield->atomic_number = atomic_number;
    shield->material = nullptr;

    // Set starting areal density (convert from g/cm2 to PhysicalUnits)
    shield->areal_density = starting_areal_density * PhysicalUnits::g_per_cm2;

    // Do NOT fit atomic number - only fit areal density
    shield->fit_atomic_number = false;

    // Configure areal density fitting with reasonable bounds
    shield->fit_areal_density = true;
    shield->lower_fit_areal_density = 0.1 * PhysicalUnits::g_per_cm2;
    shield->upper_fit_areal_density = 10.0 * PhysicalUnits::g_per_cm2;

    return shield;
  }//create_desperation_shielding(...)
}

using namespace std;

namespace FitPeaksForNuclides
{

// Development-harness hook to enable the verbose `should_debug_print()` tracing used throughout this
// file.  Not thread-safe; never enable during parallel GA optimization.
void set_debug_printout( bool enable )
{
  local_debug_printout = enable;
}

const char *automatic_roi_decision_name( const AutomaticRoiDecision decision )
{
  switch( decision )
  {
    case AutomaticRoiDecision::KeepSeparate: return "KeepSeparate";
    case AutomaticRoiDecision::MergeInseparable: return "MergeInseparable";
    case AutomaticRoiDecision::MergeInseparableWide: return "MergeInseparableWide";
    case AutomaticRoiDecision::UnmodeledFeatureBlocked: return "UnmodeledFeatureBlocked";
    case AutomaticRoiDecision::ProtectedGeometry: return "ProtectedGeometry";
    case AutomaticRoiDecision::R6LegacyBypass: return "R6LegacyBypass";
    case AutomaticRoiDecision::SourceBridgeRetained: return "SourceBridgeRetained";
    case AutomaticRoiDecision::SourceBridgeRejectedContinuum:
      return "SourceBridgeRejectedContinuum";
    case AutomaticRoiDecision::SourceBridgeRejectedFreeFeature:
      return "SourceBridgeRejectedFreeFeature";
    case AutomaticRoiDecision::InfeasiblePartition: return "InfeasiblePartition";
  }
  return "Unknown";
}

namespace detail
{
bool should_try_source_clean_recovery( const size_t num_source_anchors,
                                       const size_t num_preserved_anchors )
{
  return (num_source_anchors >= 2)
         && ((num_source_anchors - std::min(num_source_anchors, num_preserved_anchors)) >= 2);
}

bool should_accept_source_clean_challenger( const bool solve_succeeded,
                                            const size_t incumbent_preserved_anchors,
                                            const size_t candidate_preserved_anchors,
                                            const size_t incumbent_fitted_anchors,
                                            const size_t candidate_fitted_anchors,
                                            const double incumbent_score,
                                            const double candidate_score )
{
  const double unavailable_score = std::numeric_limits<double>::max();
  const bool score_improves = std::isfinite( candidate_score )
      && (candidate_score < unavailable_score)
      && ((incumbent_score == unavailable_score)
          || (std::isfinite(incumbent_score) && (candidate_score < incumbent_score)));
  return solve_succeeded
      && (candidate_preserved_anchors > incumbent_preserved_anchors)
      && (candidate_fitted_anchors >= incumbent_fitted_anchors)
      && score_improves;
}

double data_only_aicc( const double data_chi2,
                       const size_t num_data_rows,
                       const size_t num_parameters,
                       const double parameter_penalty )
{
  if( !std::isfinite(data_chi2) || (data_chi2 < 0.0) || !(parameter_penalty > 0.0)
      || (num_data_rows <= (num_parameters + 1)) )
    return std::numeric_limits<double>::max();
  return data_chi2 + parameter_penalty*num_parameters
      + parameter_penalty*num_parameters*(num_parameters + 1.0)
          / (num_data_rows - num_parameters - 1.0);
}
}//namespace detail

// Anonymous namespace for helper functions
namespace
{
  // Reduced chi-square (chi2 per degree of freedom) for a fit solution, treating a non-positive
  // number of degrees of freedom as an infinitely-bad fit.  Makes the auto-stage escalation and
  // acceptance gates deterministic instead of dividing by zero (inf/nan) when a solve is exactly-
  // or over-determined.  NOTE: only for the auto (RelActAuto) solution gates - the manual stage
  // legitimately produces dof==0 for single-peak sources and must not treat that as a failure.
  template<class SolutionType>
  double reduced_chi2( const SolutionType &sol )
  {
    return (sol.m_dof > 0) ? (sol.m_chi2 / sol.m_dof) : std::numeric_limits<double>::max();
  }

  // Forward declarations for structs used in helper functions
  struct RoiSignificanceResult
  {
    double chi2_with_peaks = 0.0;
    double chi2_continuum_only = 0.0;    // quadratic continuum-only null fit
    double chi2_reduction = 0.0;         // chi2_continuum_only - chi2_with_peaks
    double max_peak_significance = 0.0;  // Max per-peak S/sqrt(S+B); diagnostic output only
    size_t num_channels = 0;

    // Likelihood-ratio (Wilks) equivalent-z: chi2_reduction referred to a chi2 distribution with
    // dof = number of peaks in the ROI, converted to a one-sided normal quantile.  The quadratic
    // null absorbs smooth continuum curvature, so curvature cannot masquerade as peak significance.
    double equivalent_z = 0.0;

    bool has_significant_peaks = false;  // equivalent_z >= threshold
  };

  struct FixedRoiModelScore
  {
    double poisson_deviance = 0.0;
    size_t num_channels = 0;
    bool valid = false;
  };

  struct LocalMinimum
  {
    size_t channel;                  // Absolute channel number of minimum
    double synthetic_value;
    double depth_score;              // For tiebreaking (higher = better)
    double statistical_significance; // Primary criterion (lower = better breakpoint)
  };

  struct PredictedGamma
  {
    double energy;
    double expected_counts;
    RelActCalcAuto::SrcVariant source;
    size_t rel_eff_curve_index;
  };

  struct ClusteredGammaInfo
  {
    double lower;
    double upper;
    std::vector<double> gamma_energies;    // energies of gamma lines in this cluster
    std::vector<double> gamma_amplitudes;  // expected peak areas/amplitudes
    size_t joined_groups = 1;
    // Provisional source groups retain their source/curve provenance until the over-wide local
    // evidence gate completes.  Atoms are minted only afterward, so a rejected H0/Hf bridge was
    // never admitted and the exact-once atom transaction remains literal.
    std::vector<PredictedGamma> predicted_gammas;
    std::vector<detail::RoiAtom> atoms;
    uint64_t selected_partition_id = 0;
    double selected_parent_lower = std::numeric_limits<double>::quiet_NaN();
    double selected_parent_upper = std::numeric_limits<double>::quiet_NaN();
  };

  struct SelectedRoiComponentPartition
  {
    uint64_t id = 0;
    double parent_lower = 0.0;
    double parent_upper = 0.0;
    std::vector<RelActCalcAuto::RoiRange> children;
    std::shared_ptr<const SpecUtils::EnergyCalibration> calibration;
  };

  struct MarginalRejectedCluster
  {
    ClusteredGammaInfo cluster;
    std::vector<PredictedGamma> predicted_gammas;
    double expected_counts;
    double background_counts;
    double keep_significance;
  };
  /** Get photon energies and intensities for a source at a given age.

   For nuclides: Uses NuclideMixture with the specified age
   For elements: Returns xrays (age parameter ignored)
   For reactions: Returns gammas (age parameter ignored)

   \param src The source variant
   \param activity Activity in becquerels (for nuclides) or scaling factor (for elements/reactions)
   \param age Age in PhysicalUnits seconds (for nuclides only, ignored for others)
   \return Vector of EnergyRatePair with photon energies and rates
   \throws runtime_error if source is null or invalid
   */
  std::vector<SandiaDecay::EnergyRatePair> get_source_photons(
    const RelActCalcAuto::SrcVariant &src,
    const double activity,
    const double age )
  {
    std::vector<SandiaDecay::EnergyRatePair> result;

    const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( src );
    if( nuc )
    {
      // Use NuclideMixture for nuclides
      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( nuc, activity, age );
      result = mix.photons( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy );
      return result;
    }

    const SandiaDecay::Element *el = RelActCalcAuto::element( src );
    if( el )
    {
      // Convert xrays to EnergyRatePair format
      // Xrays have relative intensities, scale by activity for consistency
      for( const SandiaDecay::EnergyIntensityPair &xray : el->xrays )
      {
        // intensity is relative to strongest line (1.0 = strongest)
        // Scale by activity to get rate
        result.emplace_back( activity * xray.intensity, xray.energy );
      }

      return result;
    }

    const ReactionGamma::Reaction *rxn = RelActCalcAuto::reaction( src );
    if( rxn )
    {
      // Convert reaction gammas to EnergyRatePair format
      for( const ReactionGamma::Reaction::EnergyYield &gamma : rxn->gammas )
      {
        // abundance is yield per reaction
        // Scale by activity to get rate
        result.emplace_back( activity * gamma.abundance, gamma.energy );
      }

      return result;
    }

    throw std::runtime_error( "get_source_photons: invalid or null source" );
  }//get_source_photons(...)


  /** Get appropriate age for a source.

   For nuclides: Returns PeakDef::defaultDecayTime() if age < 0, else returns age
   For elements/reactions: Returns 0.0 (age not applicable)

   \param src The source variant
   \param age_input User-specified age (negative means use default for nuclides)
   \return Age in PhysicalUnits seconds
   */
  double get_source_age( const RelActCalcAuto::SrcVariant &src, const double age_input )
  {
    const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( src );
    if( nuc )
    {
      if( age_input >= 0.0 )
        return age_input;

      return PeakDef::defaultDecayTime( nuc );
    }

    // Elements and reactions don't have age
    return 0.0;
  }//get_source_age(...)


  /** Sets the default `useForShieldingSourceFit` and `useForManualRelEff` flags on
   each peak in `peaks`, matching the per-source-type behavior of `PeakModel::setNuclide`
   (src/PeakModel.cpp).  Used by the "Fit Source" workflow so peaks produced through
   FitPeaksForNuclides land in the user's model with the same default-on/off use flags
   as peaks created via the manual double-click + nuclide-assignment path.
   */
  void apply_default_use_flags( std::vector<PeakDef> &peaks )
  {
    for( PeakDef &peak : peaks )
    {
      bool useShield = false, useRelEff = false;
      switch( peak.sourceGammaType() )
      {
        case PeakDef::NormalGamma:
        {
          const SandiaDecay::Nuclide * const nuc = peak.parentNuclide();
          const SandiaDecay::Transition * const tr = peak.nuclearTransition();
          const int idx = peak.decayParticleIndex();
          if( nuc && tr && (idx >= 0) && (idx < static_cast<int>(tr->products.size())) )
          {
            const float energy = tr->products[idx].energy;
            useShield = PeakModel::recommendUseForFit( nuc, energy );
            useRelEff = PeakModel::recommendUseForManualRelEff( nuc, energy );
          }
          break;
        }

        case PeakDef::AnnihilationGamma:
          if( peak.parentNuclide() )
            useShield = PeakModel::recommendUseForFit( peak.parentNuclide(), 510.99891f );
          break;

        case PeakDef::XrayGamma:
        case PeakDef::SingleEscapeGamma:
        case PeakDef::DoubleEscapeGamma:
          break;
      }//switch( peak.sourceGammaType() )

      peak.useForShieldingSourceFit( useShield );
      peak.useForManualRelEff( useRelEff );
    }//for( PeakDef &peak : peaks )
  }//apply_default_use_flags(...)
}//namespace


// Fixed (non-GA-configurable) statistical safety rails.  These protect regimes where the
// significance-based tests below are not valid; they are deliberately NOT exposed to the
// genetic-algorithm config, since the GA would otherwise optimize away its own guard rails.
//
// Minimum expected peak counts for a gamma cluster to be kept: below roughly this many counts the
// Gaussian approximation to Poisson statistics (and hence the z-based keep gate) breaks down, and
// a fit is unlikely to converge meaningfully regardless of significance.
static constexpr double sm_keep_gate_min_est_counts = 15.0;

namespace
{
std::vector<RelActCalcManual::GenericPeakInfo> distinct_significant_source_anchors(
  const std::vector<RelActCalcManual::GenericPeakInfo> &anchors,
  const double minimum_z )
{
  std::vector<RelActCalcManual::GenericPeakInfo> candidates;
  for( const RelActCalcManual::GenericPeakInfo &anchor : anchors )
  {
    const double z = (anchor.m_counts_uncert > 0.0)
      ? (anchor.m_counts / anchor.m_counts_uncert) : 0.0;
    if( (z >= minimum_z) && (anchor.m_counts > 0.0)
        && std::isfinite(anchor.m_energy) && !anchor.m_source_gammas.empty() )
      candidates.push_back( anchor );
  }

  std::sort( std::begin(candidates), std::end(candidates),
    []( const RelActCalcManual::GenericPeakInfo &lhs,
        const RelActCalcManual::GenericPeakInfo &rhs ) {
      const double lhs_z = lhs.m_counts / lhs.m_counts_uncert;
      const double rhs_z = rhs.m_counts / rhs.m_counts_uncert;
      return lhs_z > rhs_z;
    } );

  std::vector<RelActCalcManual::GenericPeakInfo> result;
  for( const RelActCalcManual::GenericPeakInfo &candidate : candidates )
  {
    const bool overlaps_existing = std::any_of( std::begin(result), std::end(result),
      [&candidate]( const RelActCalcManual::GenericPeakInfo &existing ) {
        const double candidate_fwhm = std::max( 0.0, candidate.m_fwhm );
        const double existing_fwhm = std::max( 0.0, existing.m_fwhm );
        const double required_separation = 0.5*(candidate_fwhm + existing_fwhm);
        return std::fabs(candidate.m_energy - existing.m_energy) < required_separation;
      } );
    if( !overlaps_existing )
      result.push_back( candidate );
  }

  std::sort( std::begin(result), std::end(result),
    []( const RelActCalcManual::GenericPeakInfo &lhs,
        const RelActCalcManual::GenericPeakInfo &rhs ) {
      return lhs.m_energy < rhs.m_energy;
    } );
  return result;
}

double predicted_source_anchor_counts(
  const RelActCalcAuto::RelActAutoSolution &solution,
  const size_t rel_eff_curve_index,
  const RelActCalcManual::GenericPeakInfo &anchor,
  const double live_time )
{
  if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status)
      || (rel_eff_curve_index >= solution.m_rel_activities.size())
      || !(live_time > 0.0) )
    return 0.0;

  const double rel_eff = solution.relative_efficiency( anchor.m_energy, rel_eff_curve_index );
  if( !std::isfinite(rel_eff) || !(rel_eff > 0.0) )
    return 0.0;

  double predicted = 0.0;
  for( const RelActCalcManual::GenericLineInfo &line : anchor.m_source_gammas )
  {
    const std::vector<RelActCalcAuto::NuclideRelAct> &activities
      = solution.m_rel_activities[rel_eff_curve_index];
    const std::vector<RelActCalcAuto::NuclideRelAct>::const_iterator found = std::find_if(
      std::begin(activities), std::end(activities),
      [&line]( const RelActCalcAuto::NuclideRelAct &activity ) {
        return activity.name() == line.m_isotope;
      } );
    if( found != std::end(activities) )
      predicted += found->rel_activity * live_time * line.m_yield * rel_eff;
  }
  return std::isfinite(predicted) ? std::max( 0.0, predicted ) : 0.0;
}

size_t preserved_source_anchor_count(
  const RelActCalcAuto::RelActAutoSolution &solution,
  const size_t rel_eff_curve_index,
  const std::vector<RelActCalcManual::GenericPeakInfo> &anchors,
  const double live_time )
{
  return static_cast<size_t>( std::count_if( std::begin(anchors), std::end(anchors),
    [&solution, rel_eff_curve_index, live_time](
        const RelActCalcManual::GenericPeakInfo &anchor ) {
      return predicted_source_anchor_counts( solution, rel_eff_curve_index, anchor, live_time )
             > sm_keep_gate_min_est_counts;
    } ) );
}

std::vector<RelActCalcManual::GenericPeakInfo> preserved_source_anchors(
  const RelActCalcAuto::RelActAutoSolution &solution,
  const size_t rel_eff_curve_index,
  const std::vector<RelActCalcManual::GenericPeakInfo> &anchors,
  const double live_time )
{
  std::vector<RelActCalcManual::GenericPeakInfo> result;
  std::copy_if( std::begin(anchors), std::end(anchors), std::back_inserter(result),
    [&solution, rel_eff_curve_index, live_time](
        const RelActCalcManual::GenericPeakInfo &anchor ) {
      return predicted_source_anchor_counts( solution, rel_eff_curve_index, anchor, live_time )
             > sm_keep_gate_min_est_counts;
    } );
  return result;
}

}//namespace

// Minimum detection significance of an auto-search peak for a dropped energy-extreme ROI to be
// restored on the strength of the data (see the "data should override the model" block).  Replaces
// a former absolute 100-count floor.
static constexpr double sm_edge_roi_restore_min_z = 4.0;

// Adaptive ROI sideband extension (detail::extend_roi_by_sidebands) internals.
// Extension proceeds in blocks of this many FWHM (widened to at least 2 channels).
static constexpr double sm_extend_block_fwhm = 0.375;
// Cumulative slow-drift guard: extension stops once (mean block z)^2 x n_blocks exceeds this,
// catching gentle continuum curvature that individual blocks would pass.
static constexpr double sm_extend_drift_chi2 = 4.0;
// Extra low-side core allowance (in FWHM) when the fit will use a skewed peak shape, since skew
// puts appreciable peak area below the Gaussian core.
static constexpr double sm_skew_low_side_extra_fwhm = 0.75;
// Sideband retained beyond the core (in FWHM) when a ROI is shrunk after dropping edge peaks
// during the observable-peaks filter.
static constexpr double sm_post_drop_sideband_fwhm = 1.0;

// Width prior for the per-ROI continuum-order AICc selection: a quadratic candidate is only
// offered for ROIs at least this many FWHM wide (narrower windows cannot support curvature).
static constexpr double sm_quad_cont_min_roi_fwhm = 4.0;

// Loose sideband-asymmetry pre-filter for the step-continuum trial fit: the low-side continuum
// must exceed the high side by at least this many sigma before the (cheap, but not free) trial
// runs.  Deliberately permissive - it exists to skip clearly-symmetric strong peaks, not to make
// the step decision (the trial fit does that); see trial_step_continuum.
static constexpr double sm_step_trial_min_asym_z = 1.0;

// Half-width (in FWHM, each side of the peak) of a TIGHT ROI seeded directly for a found+matched
// auto-search peak - a data-confirmed source line that bypasses the predicted-signal keep-gate (see
// seed_tight_rois_for_found_peaks).  Deliberately tight (not the adaptive ~4-FWHM extent) so the
// output significance test is not diluted over a wide window.
static constexpr double sm_found_peak_roi_half_num_fwhm = 1.5;


// R6 auto co-fit of strong unmodeled interfering lines (detail::find_strong_unmodeled_interferers).
// A requested-source gamma this many FWHM (evaluated at the interfering line) from a strong
// unmodeled line triggers a co-fit check.
static constexpr double sm_interferer_near_num_fwhm = 2.0;
// A confirming foreground auto-search peak must lie within this many FWHM of the interfering line
// to count as data-confirmation.  Dedicated background discovery remains deferred.
static constexpr double sm_interferer_confirm_num_fwhm = 0.5;
// Minimum detection significance (area/uncert) of the confirming auto-search peak.
static constexpr double sm_interferer_min_detect_z = 5.0;
// If any requested-source gamma is within this many FWHM of the candidate line, the source's own
// chain is assumed to explain it and no interferer is added.
static constexpr double sm_interferer_source_owns_num_fwhm = 1.0;
// A single-line source vs a single-line interferer closer than this many FWHM is an unresolvable
// doublet: skip and warn rather than co-fit (co-fitting would just split one peak).
static constexpr double sm_interferer_doublet_min_fwhm = 1.0;
// Bound the number of automatically-added nuisance activities.  Allowing every detected NORM
// parent onto the extra curve can make a multi-source solve poorly conditioned.
static constexpr size_t sm_max_auto_interferer_nuclides = 2;

// R2 bounded fit-then-prune rescue.  These are statistical safety rails, not tuning genes.
static constexpr double sm_rescue_z_fraction = 0.7;
static constexpr double sm_rescue_guard_num_fwhm = 1.0;
static constexpr size_t sm_max_rescued_rois = 4;

#if( PERFORM_DEVELOPER_CHECKS )
bool sm_bounded_rescue_enabled_for_test = true;
bool sm_force_next_rescue_admission_failure_for_test = false;
bool sm_force_next_rescue_evaluation_failure_for_test = false;
#endif

bool bounded_rescue_enabled()
{
#if( PERFORM_DEVELOPER_CHECKS )
  return sm_bounded_rescue_enabled_for_test;
#else
  return true;
#endif
}


namespace detail
{

#if( PERFORM_DEVELOPER_CHECKS )
void set_bounded_rescue_enabled_for_test( const bool enabled )
{
  sm_bounded_rescue_enabled_for_test = enabled;
}


void force_next_bounded_rescue_admission_failure_for_test()
{
  sm_force_next_rescue_admission_failure_for_test = true;
}


void force_next_bounded_rescue_evaluation_failure_for_test()
{
  sm_force_next_rescue_evaluation_failure_for_test = true;
}
#endif

double LocalContinuumEstimate::integral( const double x0, const double x1 ) const
{
  if( !valid || (x1 <= x0) )
    return 0.0;

  const double area = PeakContinuum::offset_eqn_integral(
      coeffs, PeakContinuum::OffsetType::Linear, x0, x1, reference_energy );

  return std::max( 0.0, area );
}//LocalContinuumEstimate::integral


LocalContinuumEstimate estimate_local_continuum(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const double region_lower,
  const double region_upper,
  const double fwhm,
  const double sideband_num_fwhm,
  const std::function<double(double,double)> &predicted_signal,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks )
{
  LocalContinuumEstimate result;

  if( !foreground || !foreground->channel_energies() || (foreground->num_gamma_channels() < 4)
      || !std::isfinite(region_lower) || !std::isfinite(region_upper)
      || (region_upper <= region_lower) || !std::isfinite(fwhm) || (fwhm <= 0.0) )
  {
    return result;
  }

  const size_t nchannel = foreground->num_gamma_channels();
  const size_t lowchannel = foreground->find_gamma_channel( static_cast<float>(region_lower) );
  const size_t highchannel = foreground->find_gamma_channel( static_cast<float>(region_upper) );

  if( (lowchannel >= highchannel) || (highchannel >= nchannel) )
    return result;

  // Average at least 2 channels per side; more when the resolution spans many channels.
  const double channel_width = std::max( 1.0e-6,
      static_cast<double>( foreground->gamma_channel_width( lowchannel ) ) );
  const size_t sideband_channels = std::max( size_t(2),
      static_cast<size_t>( std::llround( (sideband_num_fwhm * fwhm) / channel_width ) ) );

  const std::shared_ptr<const SpecUtils::EnergyCalibration> cal = foreground->energy_calibration();
  if( !cal || !cal->valid() )
    return result;

  // Returns true when an unfit auto-search peak (within half its own FWHM of its mean) overlaps
  // the window - such structure would bias the sideband average away from the true continuum.
  const auto sideband_contaminated = [&unfit_auto_peaks]( const double w_lo, const double w_hi ) -> bool {
    for( const std::shared_ptr<const PeakDef> &p : unfit_auto_peaks )
    {
      if( !p )
        continue;
      const double half_w = 0.5 * p->fwhm();
      if( ((p->mean() + half_w) > w_lo) && ((p->mean() - half_w) < w_hi) )
        return true;
    }
    return false;
  };

  // Locate a usable sideband window on one side: start adjacent to the region edge and, when an
  // unfit auto-search peak sits in the window, slide one window-width further from the region
  // (up to a few tries) so the sample measures continuum rather than an unrelated peak.
  // Returns false when no clean in-spectrum window exists.
  const auto find_sideband = [&]( const bool low_side, size_t &first_ch, size_t &last_ch ) -> bool {
    const size_t max_shifts = 3;
    for( size_t shift = 0; shift <= max_shifts; ++shift )
    {
      if( low_side )
      {
        const size_t offset = (shift + 1) * sideband_channels;
        if( lowchannel < offset )
          return false;  // would extend past the first channel
        first_ch = lowchannel - offset;
        last_ch = first_ch + sideband_channels - 1;
      }else
      {
        first_ch = highchannel + 1 + shift*sideband_channels;
        last_ch = first_ch + sideband_channels - 1;
        if( last_ch >= nchannel )
          return false;  // would extend past the last channel
      }

      const double w_lo = cal->energy_for_channel( static_cast<double>(first_ch) );
      const double w_hi = cal->energy_for_channel( static_cast<double>(last_ch + 1) );
      if( !sideband_contaminated( w_lo, w_hi ) )
        return true;
    }
    return false;
  };//find_sideband lambda

  size_t low_first = 0, low_last = 0, up_first = 0, up_last = 0;
  if( !find_sideband( true, low_first, low_last ) || !find_sideband( false, up_first, up_last ) )
    return result;

  const double lower_low_energy = cal->energy_for_channel( static_cast<double>(low_first) );
  const double lower_up_energy = cal->energy_for_channel( static_cast<double>(low_last + 1) );
  const double upper_low_energy = cal->energy_for_channel( static_cast<double>(up_first) );
  const double upper_up_energy = cal->energy_for_channel( static_cast<double>(up_last + 1) );

  const double lower_dx = lower_up_energy - lower_low_energy;
  const double upper_dx = upper_up_energy - upper_low_energy;
  if( (lower_dx <= 0.0) || (upper_dx <= 0.0) )
    return result;

  const double lower_raw = foreground->gamma_channels_sum( low_first, low_last );
  const double upper_raw = foreground->gamma_channels_sum( up_first, up_last );

  // Subtract the caller-predicted signal (e.g., tails of the cluster's own gammas, or of known
  // neighboring peaks) so the sideband measures continuum rather than peak leakage; clamp at 0.
  double lower_net = lower_raw, upper_net = upper_raw;
  if( predicted_signal )
  {
    lower_net = std::max( 0.0, lower_raw - predicted_signal( lower_low_energy, lower_up_energy ) );
    upper_net = std::max( 0.0, upper_raw - predicted_signal( upper_low_energy, upper_up_energy ) );
  }

  result.reference_energy = region_lower;
  result.lower_sideband_counts = lower_net;
  result.upper_sideband_counts = upper_net;
  result.lower_sideband_raw_counts = lower_raw;
  result.upper_sideband_raw_counts = upper_raw;
  result.lower_sideband_lo = lower_low_energy;
  result.lower_sideband_hi = lower_up_energy;
  result.upper_sideband_lo = upper_low_energy;
  result.upper_sideband_hi = upper_up_energy;

  // Linear density y = m*x + b (x relative to reference_energy) whose integral over each sideband
  // window equals that window's (signal-subtracted) counts - the same algebra as
  // PeakContinuum::eqn_from_offsets, generalized to relocated windows and net counts:
  //   lower_net = b*lower_dx + 0.5*m*lower_sqr_diff
  //   upper_net = b*upper_dx + 0.5*m*upper_sqr_diff
  const double lower_x1_rel = lower_low_energy - result.reference_energy;
  const double lower_x2_rel = lower_up_energy - result.reference_energy;
  const double upper_x1_rel = upper_low_energy - result.reference_energy;
  const double upper_x2_rel = upper_up_energy - result.reference_energy;
  const double lower_sqr_diff = lower_x2_rel*lower_x2_rel - lower_x1_rel*lower_x1_rel;
  const double upper_sqr_diff = upper_x2_rel*upper_x2_rel - upper_x1_rel*upper_x1_rel;

  double m = 0.0, b = 0.0;
  const double denom = upper_sqr_diff - (upper_dx * lower_sqr_diff / lower_dx);
  if( std::fabs(denom) < FLT_EPSILON )
  {
    m = 0.0;
    b = (0.5 * lower_net / lower_dx) + (0.5 * upper_net / upper_dx);
  }else
  {
    m = 2.0 * (upper_net - (upper_dx * lower_net / lower_dx)) / denom;
    b = (lower_net - 0.5*m*lower_sqr_diff) / lower_dx;
  }

  result.coeffs[0] = b;
  result.coeffs[1] = m;
  result.valid = std::isfinite(b) && std::isfinite(m);

  return result;
}//estimate_local_continuum


double GlobalContinuumEstimate::integral( double x0, double x1 ) const
{
  if( !valid() || (x1 <= x0) )
    return 0.0;
  const double v = snip->gamma_integral( static_cast<float>(x0), static_cast<float>(x1) );
  return (v > 0.0) ? v : 0.0;
}//GlobalContinuumEstimate::integral


double GlobalContinuumEstimate::density_at( double E ) const
{
  // Counts/keV over a narrow window (gamma_integral clamps to the spectrum, so no edge throw).
  if( !valid() )
    return 0.0;
  return integral( E - 0.5, E + 0.5 );  // 1 keV window => integral == density
}//GlobalContinuumEstimate::density_at


double GlobalContinuumEstimate::integral_variance( double x0, double x1 ) const
{
  // Conservative Poisson variance: the DATA counts over [x0,x1] - an upper bound, largest exactly
  // where peaks sit (which is where the SNIP continuum is least trustworthy).
  if( !valid() || (x1 <= x0) )
    return 0.0;
  const double v = foreground->gamma_integral( static_cast<float>(x0), static_cast<float>(x1) );
  return (v > 0.0) ? v : 0.0;
}//GlobalContinuumEstimate::integral_variance


GlobalContinuumEstimate make_global_continuum(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::function<double(double)> &fwhm_at_energy,
  PeakFitUtils::CoarseResolutionType det_type,
  double restrict_lower_energy,
  double restrict_upper_energy )
{
  GlobalContinuumEstimate est;

  if( !foreground || !fwhm_at_energy )
    return est;

  try
  {
    // Per-class SNIP parameters settled during the SNIP work: HPGe = 2.0xFWHM / order 2 / 3-ch
    // presmooth / LLS on; NaI/LaBr/CZT = 1.5xFWHM / order 2 / 7-ch presmooth / LLS off.  Both
    // restricted to the valid extent so the low-energy detector turn-on cannot pull it up.
    const bool is_hpge = (det_type == PeakFitUtils::CoarseResolutionType::High);
    const double num_fwhm_window       = is_hpge ? 2.0 : 1.5;
    const int    filter_order          = 2;
    const int    presmooth_halfwidth   = is_hpge ? 1 : 3;   // 1 = 3-ch boxcar, 3 = 7-ch boxcar
    const bool   lls                   = is_hpge;

    est.snip = estimateContinuum( foreground, fwhm_at_energy, num_fwhm_window, filter_order,
                                  presmooth_halfwidth, lls, restrict_lower_energy, restrict_upper_energy );
    est.foreground = foreground;
    est.built = static_cast<bool>( est.snip );
  }catch( const std::exception & )
  {
    est = GlobalContinuumEstimate();  // invalid => callers fall back to local estimation
  }

  return est;
}//make_global_continuum


namespace
{
thread_local std::vector<RoiBoundaryShadowResult> s_roi_boundary_shadow_diagnostics;
thread_local std::vector<AutomaticRoiDecisionDiagnostic> s_automatic_roi_diagnostics;
}

void record_automatic_roi_diagnostic( const AutomaticRoiDecisionDiagnostic &diagnostic )
{
  s_automatic_roi_diagnostics.push_back( diagnostic );
}

std::vector<AutomaticRoiDecisionDiagnostic> take_automatic_roi_diagnostics()
{
  std::vector<AutomaticRoiDecisionDiagnostic> answer;
  answer.swap( s_automatic_roi_diagnostics );
  return answer;
}


RoiBoundaryShadowResult optimize_roi_boundaries_shadow(
  const std::vector<RoiBoundaryShadowGroup> &input_groups,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const GlobalContinuumEstimate &global_continuum,
  const std::function<double(double)> &fwhm_at_energy,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
  const double catastrophic_max_fwhm_width,
  const double roi_core_num_fwhm )
{
  RoiBoundaryShadowResult result;
  if( input_groups.empty() )
  {
    result.fallback_reason = "no accepted source-gamma groups";
    return result;
  }
  if( !foreground || !foreground->channel_energies() || !global_continuum.valid()
      || !fwhm_at_energy )
  {
    result.fallback_reason = "invalid foreground, FWHM, or shared SNIP continuum";
    return result;
  }

  std::vector<RoiBoundaryShadowGroup> groups = input_groups;
  groups.erase( std::remove_if( std::begin(groups), std::end(groups),
    []( const RoiBoundaryShadowGroup &group ) {
      return group.gamma_energies.empty();
    } ), std::end(groups) );
  std::sort( std::begin(groups), std::end(groups),
    []( const RoiBoundaryShadowGroup &lhs, const RoiBoundaryShadowGroup &rhs ) {
      const std::pair<std::vector<double>::const_iterator,
                      std::vector<double>::const_iterator> lhs_minmax
        = std::minmax_element( std::begin(lhs.gamma_energies), std::end(lhs.gamma_energies) );
      const std::pair<std::vector<double>::const_iterator,
                      std::vector<double>::const_iterator> rhs_minmax
        = std::minmax_element( std::begin(rhs.gamma_energies), std::end(rhs.gamma_energies) );
      if( *lhs_minmax.first != *rhs_minmax.first )
        return *lhs_minmax.first < *rhs_minmax.first;
      if( *lhs_minmax.second != *rhs_minmax.second )
        return *lhs_minmax.second < *rhs_minmax.second;
      if( lhs.legacy_lower != rhs.legacy_lower )
        return lhs.legacy_lower < rhs.legacy_lower;
      return lhs.legacy_upper < rhs.legacy_upper;
    } );
  if( groups.empty() )
  {
    result.fallback_reason = "accepted groups had no gamma energies";
    return result;
  }

  const double spectrum_lower = foreground->gamma_channel_lower( 0 );
  const double spectrum_upper = foreground->gamma_channel_upper(
      foreground->num_gamma_channels() - 1 );
  std::vector<double> core_lower, core_upper;
  std::set<double> candidate_set;
  for( RoiBoundaryShadowGroup &group : groups )
  {
    std::sort( std::begin(group.gamma_energies), std::end(group.gamma_energies) );
    const std::pair<std::vector<double>::const_iterator,
                    std::vector<double>::const_iterator> minmax
      = std::minmax_element( std::begin(group.gamma_energies), std::end(group.gamma_energies) );
    const double lo_fwhm = fwhm_at_energy( *minmax.first );
    const double hi_fwhm = fwhm_at_energy( *minmax.second );
    if( !std::isfinite(lo_fwhm) || !std::isfinite(hi_fwhm)
        || !(lo_fwhm > 0.0) || !(hi_fwhm > 0.0) )
    {
      result.fallback_reason = "invalid FWHM over a source-group core";
      return result;
    }
    core_lower.push_back( std::max(spectrum_lower,
        *minmax.first - roi_core_num_fwhm*lo_fwhm) );
    core_upper.push_back( std::min(spectrum_upper,
        *minmax.second + roi_core_num_fwhm*hi_fwhm) );
    candidate_set.insert( group.legacy_lower );
    candidate_set.insert( group.legacy_upper );
    candidate_set.insert( core_lower.back() );
    candidate_set.insert( core_upper.back() );
  }

  const std::vector<float> &channel_energies = *foreground->channel_energies();
  const std::shared_ptr<const SpecUtils::Measurement> snip = global_continuum.snip;
  struct UnmodeledPeakExclusion
  {
    double lower;
    double upper;
  };
  std::vector<UnmodeledPeakExclusion> unmodeled_exclusions;
  for( const std::shared_ptr<const PeakDef> &peak : unfit_auto_peaks )
  {
    if( !peak || (peak->mean() <= spectrum_lower) || (peak->mean() >= spectrum_upper) )
      continue;
    const double fwhm = fwhm_at_energy( peak->mean() );
    if( !std::isfinite(fwhm) || !(fwhm > 0.0) )
      continue;
    UnmodeledPeakExclusion exclusion;
    exclusion.lower = std::max( spectrum_lower, peak->mean() - fwhm );
    exclusion.upper = std::min( spectrum_upper, peak->mean() + fwhm );
    if( exclusion.upper > exclusion.lower )
    {
      unmodeled_exclusions.push_back( exclusion );
      candidate_set.insert( exclusion.lower );
      candidate_set.insert( exclusion.upper );
    }
  }

  for( size_t i = 0; i + 1 < groups.size(); ++i )
  {
    const double gap_lo = core_upper[i];
    const double gap_hi = core_lower[i + 1];
    if( !(gap_hi > gap_lo) )
      continue;
    candidate_set.insert( 0.5*(gap_lo + gap_hi) );

    const size_t first_channel = foreground->find_gamma_channel(
        static_cast<float>(gap_lo) );
    const size_t last_channel = foreground->find_gamma_channel(
        static_cast<float>(gap_hi) );
    if( last_channel > (first_channel + 2) )
    {
      size_t valley_channel = first_channel + 1;
      size_t curvature_channel = valley_channel;
      double valley_value = std::numeric_limits<double>::max();
      double max_curvature = -1.0;
      for( size_t channel = first_channel + 1; channel < last_channel; ++channel )
      {
        const double value = snip->gamma_channel_content( channel );
        if( value < valley_value )
        {
          valley_value = value;
          valley_channel = channel;
        }
        const double curvature = std::fabs(
            snip->gamma_channel_content(channel + 1) - 2.0*value
            + snip->gamma_channel_content(channel - 1) );
        if( curvature > max_curvature )
        {
          max_curvature = curvature;
          curvature_channel = channel;
        }
      }
      candidate_set.insert( channel_energies[valley_channel] );
      candidate_set.insert( channel_energies[curvature_channel] );
    }
  }

  std::vector<double> candidates;
  for( const double energy : candidate_set )
  {
    if( std::isfinite(energy) && (energy >= spectrum_lower) && (energy <= spectrum_upper) )
      candidates.push_back( energy );
  }
  std::sort( std::begin(candidates), std::end(candidates) );
  candidates.erase( std::unique(std::begin(candidates), std::end(candidates),
    []( const double lhs, const double rhs ) {
      return std::fabs(lhs - rhs) < 0.05;
    } ), std::end(candidates) );
  if( candidates.size() < 2 )
  {
    result.fallback_reason = "too few finite boundary candidates";
    return result;
  }

  struct IntervalScore
  {
    bool valid = false;
    double score = std::numeric_limits<double>::max();
    double mismatch = 0.0;
    PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;
    size_t num_channels = 0;
    size_t start_channel = 0;
    double width_fwhm = 0.0;
    std::vector<double> predictions;
  };
  std::map<std::pair<size_t,size_t>,IntervalScore> score_cache;
  const auto score_interval = [&]( const size_t lower_index, const size_t upper_index )
      -> IntervalScore {
    const std::pair<size_t,size_t> key( lower_index, upper_index );
    const auto cached = score_cache.find( key );
    if( cached != std::end(score_cache) )
      return cached->second;

    IntervalScore best;
    const double lower = candidates[lower_index];
    const double upper = candidates[upper_index];
    const double midpoint = 0.5*(lower + upper);
    const double midpoint_fwhm = fwhm_at_energy( midpoint );
    if( !(upper > lower) || !std::isfinite(midpoint_fwhm) || !(midpoint_fwhm > 0.0) )
      return score_cache.emplace( key, best ).first->second;
    best.width_fwhm = (upper - lower) / midpoint_fwhm;
    if( (catastrophic_max_fwhm_width > 0.0)
        && (best.width_fwhm > catastrophic_max_fwhm_width) )
      return score_cache.emplace( key, best ).first->second;

    const size_t start_channel = foreground->find_gamma_channel(
        static_cast<float>(lower) );
    const size_t end_channel = std::min( foreground->find_gamma_channel(
        static_cast<float>(upper) ), foreground->num_gamma_channels() - 1 );
    if( end_channel <= (start_channel + 4) )
      return score_cache.emplace( key, best ).first->second;
    const size_t nbin = end_channel - start_channel;
    best.num_channels = nbin;
    best.start_channel = start_channel;
    std::vector<float> snip_counts( nbin );
    std::vector<float> raw_variances( nbin );
    for( size_t i = 0; i < nbin; ++i )
    {
      snip_counts[i] = snip->gamma_channel_content( start_channel + i );
      raw_variances[i] = static_cast<float>( std::max( 1.0,
          static_cast<double>(foreground->gamma_channel_content(start_channel + i)) ) );
    }

    const PeakContinuum::OffsetType families[] = {
      PeakContinuum::OffsetType::Linear,
      PeakContinuum::OffsetType::Quadratic,
      PeakContinuum::OffsetType::FlatStep,
      PeakContinuum::OffsetType::LinearStep
    };
    const std::vector<double> no_means, no_sigmas;
    const std::vector<PeakDef> no_fixed_peaks;
    for( const PeakContinuum::OffsetType family : families )
    {
      const size_t num_parameters = PeakContinuum::num_parameters( family );
      if( nbin <= (num_parameters + 2) )
        continue;
      std::vector<double> amplitudes, coefficients, amplitude_uncerts, coefficient_uncerts;
      std::vector<double> predictions( nbin, 0.0 );
      try
      {
        static_cast<void>( PeakFit::fit_amp_and_offset_imp<PeakDef,double>(
            &channel_energies[start_channel], snip_counts.data(), raw_variances.data(), nbin,
            family, 0.0, midpoint, no_means, no_sigmas, no_fixed_peaks,
            PeakDef::SkewType::NoSkew, nullptr, amplitudes, coefficients,
            amplitude_uncerts, coefficient_uncerts, predictions.data() ) );
      }catch( const std::exception & )
      {
        continue;
      }
      if( coefficients.size() != num_parameters )
        continue;

      bool physically_valid = true;
      double mismatch = 0.0;
      for( size_t i = 0; i < nbin; ++i )
      {
        const double prediction = predictions[i];
        if( !std::isfinite(prediction) || (prediction < -1.0e-6) )
        {
          physically_valid = false;
          break;
        }
        const double residual = snip_counts[i] - prediction;
        const double variance = std::max( 1.0,
            static_cast<double>(foreground->gamma_channel_content(start_channel + i)) );
        mismatch += residual*residual / variance;
      }
      if( !physically_valid )
        continue;

      const double k = static_cast<double>(num_parameters);
      const double n = static_cast<double>(nbin);
      const double aicc = mismatch + 2.0*k
          + ((n > (k + 1.0)) ? (2.0*k*(k + 1.0)/(n - k - 1.0)) : 1.0e6);
      // Normalize the information score by the modeled channel count.  Boundary alternatives
      // legitimately cover different amounts of continuum; raw AICc would otherwise prefer a
      // short interval simply because it omits observations.  The DP sum still charges the
      // complexity term once per ROI, so gratuitous splitting is penalized.
      const double normalized_aicc = aicc / n;
      if( normalized_aicc < best.score )
      {
        best.valid = true;
        best.score = normalized_aicc;
        best.mismatch = mismatch / n;
        best.continuum_type = family;
        best.predictions = std::move( predictions );
      }
    }
    score_cache[key] = best;
    return best;
  };

  struct Path
  {
    bool valid = false;
    double cost = std::numeric_limits<double>::max();
    std::vector<RoiBoundaryShadowInterval> intervals;
  };
  const auto intersects_unmodeled_exclusion = [&unmodeled_exclusions](
      const double lower, const double upper ) -> bool {
    return std::any_of( std::begin(unmodeled_exclusions),
      std::end(unmodeled_exclusions), [lower, upper]( const UnmodeledPeakExclusion &exclusion ) {
        return (lower < exclusion.upper) && (upper > exclusion.lower);
      } );
  };
  const auto count_unmodeled_exclusions = [&unmodeled_exclusions](
      const double lower, const double upper ) -> size_t {
    return static_cast<size_t>( std::count_if( std::begin(unmodeled_exclusions),
      std::end(unmodeled_exclusions), [lower, upper]( const UnmodeledPeakExclusion &exclusion ) {
        return (lower < exclusion.upper) && (upper > exclusion.lower);
      } ) );
  };
  const auto permissible_next_starts = [&candidates, &unmodeled_exclusions](
      const size_t end_index ) {
    std::vector<size_t> indices( 1, end_index );
    for( const UnmodeledPeakExclusion &exclusion : unmodeled_exclusions )
    {
      if( std::fabs(candidates[end_index] - exclusion.lower) >= 0.05 )
        continue;
      const auto upper = std::lower_bound(
          std::begin(candidates), std::end(candidates), exclusion.upper - 0.05);
      if( (upper != std::end(candidates))
          && (std::fabs(*upper - exclusion.upper) < 0.05) )
        indices.push_back( static_cast<size_t>(std::distance(std::begin(candidates), upper)) );
    }
    std::sort( std::begin(indices), std::end(indices) );
    indices.erase( std::unique(std::begin(indices), std::end(indices)), std::end(indices) );
    return indices;
  };
  std::map<std::pair<size_t,size_t>,Path> memo;
  std::function<Path(size_t,size_t)> solve = [&]( const size_t group_index,
                                                   const size_t start_index ) -> Path {
    const std::pair<size_t,size_t> key( group_index, start_index );
    const auto found = memo.find( key );
    if( found != std::end(memo) )
      return found->second;

    Path best_path;
    if( (group_index >= groups.size()) || (candidates[start_index] > core_lower[group_index]) )
      return memo.emplace( key, best_path ).first->second;

    double required_core_upper = -std::numeric_limits<double>::max();
    for( size_t last_group = group_index; last_group < groups.size(); ++last_group )
    {
      required_core_upper = std::max( required_core_upper, core_upper[last_group] );
      for( size_t end_index = start_index + 1; end_index < candidates.size(); ++end_index )
      {
        const double end = candidates[end_index];
        if( end < required_core_upper )
          continue;
        if( (last_group + 1 < groups.size()) && (end > core_lower[last_group + 1]) )
          break;
        if( intersects_unmodeled_exclusion(candidates[start_index], end) )
          continue;
        const IntervalScore interval_score = score_interval( start_index, end_index );
        if( !interval_score.valid )
          continue;

        Path suffix;
        if( last_group + 1 < groups.size() )
        {
          for( const size_t next_start : permissible_next_starts(end_index) )
          {
            if( candidates[next_start] > core_lower[last_group + 1] )
              continue;
            const Path candidate_suffix = solve( last_group + 1, next_start );
            if( candidate_suffix.valid && (!suffix.valid || (candidate_suffix.cost < suffix.cost)) )
              suffix = candidate_suffix;
          }
          if( !suffix.valid )
            continue;
        }else
        {
          suffix.valid = true;
          suffix.cost = 0.0;
        }

        const double total_cost = interval_score.score + suffix.cost;
        if( total_cost >= best_path.cost )
          continue;

        RoiBoundaryShadowInterval interval;
        interval.lower = candidates[start_index];
        interval.upper = end;
        interval.legacy_lower = groups[group_index].legacy_lower;
        interval.legacy_upper = groups[last_group].legacy_upper;
        interval.width_fwhm = interval_score.width_fwhm;
        interval.num_channels = interval_score.num_channels;
        interval.continuum_type = interval_score.continuum_type;
        interval.normalized_continuum_mismatch = interval_score.mismatch;
        interval.interval_score = interval_score.score;
        interval.first_group = group_index;
        interval.last_group = last_group;
        for( size_t covered_group = group_index;
             covered_group <= last_group; ++covered_group )
        {
        interval.group_gamma_energies.insert( std::end(interval.group_gamma_energies),
              std::begin(groups[covered_group].gamma_energies),
              std::end(groups[covered_group].gamma_energies) );
        }
        const size_t sample_stride = std::max<size_t>(
            1, (interval_score.num_channels + 119) / 120);
        for( size_t sample = 0; sample < interval_score.num_channels;
             sample += sample_stride )
        {
          const size_t channel = interval_score.start_channel + sample;
          interval.profile_energies.push_back( 0.5 * (
              foreground->gamma_channel_lower(channel)
              + foreground->gamma_channel_upper(channel)) );
          interval.profile_foreground.push_back(
              foreground->gamma_channel_content(channel) );
          interval.profile_snip.push_back( snip->gamma_channel_content(channel) );
          interval.profile_continuum.push_back(
              sample < interval_score.predictions.size()
                ? interval_score.predictions[sample] : 0.0 );
        }
        for( const std::shared_ptr<const PeakDef> &peak : unfit_auto_peaks )
        {
          if( peak && (peak->mean() >= std::min(interval.lower, interval.legacy_lower))
              && (peak->mean() <= std::max(interval.upper, interval.legacy_upper)) )
            interval.unmodeled_peak_energies.push_back( peak->mean() );
        }
        interval.unmodeled_peak_conflicts = count_unmodeled_exclusions(
            interval.legacy_lower, interval.legacy_upper );
        if( last_group > group_index )
          interval.reason = "merge: joint continuum score favors adjacent source groups";
        else if( (std::fabs(interval.lower - interval.legacy_lower) > 0.05)
                 || (std::fabs(interval.upper - interval.legacy_upper) > 0.05) )
          interval.reason = "boundary adjustment: SNIP valley/curvature or group core";
        else
          interval.reason = "legacy bounds retained";

        best_path.valid = true;
        best_path.cost = total_cost;
        best_path.intervals.clear();
        best_path.intervals.push_back( interval );
        best_path.intervals.insert( std::end(best_path.intervals),
            std::begin(suffix.intervals), std::end(suffix.intervals) );
      }
    }
    memo[key] = best_path;
    return best_path;
  };

  Path best;
  for( size_t start_index = 0; start_index < candidates.size(); ++start_index )
  {
    if( candidates[start_index] > core_lower.front() )
      break;
    const Path candidate = solve( 0, start_index );
    if( candidate.valid && (candidate.cost < best.cost) )
      best = candidate;
  }
  if( !best.valid )
  {
    result.fallback_reason = "no feasible non-overlapping SNIP partition";
    return result;
  }

  const auto candidate_index_for = [&candidates]( const double energy ) -> size_t {
    const auto pos = std::lower_bound( std::begin(candidates), std::end(candidates), energy );
    if( (pos != std::end(candidates)) && (std::fabs(*pos - energy) < 0.05) )
      return static_cast<size_t>( std::distance(std::begin(candidates), pos) );
    if( (pos != std::begin(candidates)) && (std::fabs(*std::prev(pos) - energy) < 0.05) )
      return static_cast<size_t>( std::distance(std::begin(candidates), std::prev(pos)) );
    return candidates.size();
  };

  std::set<std::pair<double,double>> unique_legacy_intervals;
  for( const RoiBoundaryShadowGroup &group : groups )
    unique_legacy_intervals.insert( {group.legacy_lower, group.legacy_upper} );
  bool legacy_total_valid = true;
  double legacy_total_score = 0.0;
  for( const std::pair<double,double> &legacy : unique_legacy_intervals )
  {
    const size_t lower_index = candidate_index_for( legacy.first );
    const size_t upper_index = candidate_index_for( legacy.second );
    const IntervalScore score = ((lower_index < candidates.size())
        && (upper_index < candidates.size()) && (upper_index > lower_index))
      ? score_interval( lower_index, upper_index ) : IntervalScore();
    if( !score.valid )
    {
      legacy_total_valid = false;
      break;
    }
    legacy_total_score += score.score;
  }

  for( RoiBoundaryShadowInterval &interval : best.intervals )
  {
    std::set<std::pair<double,double>> covered_legacy_intervals;
    for( size_t group_index = interval.first_group;
         group_index <= interval.last_group; ++group_index )
      covered_legacy_intervals.insert( {groups[group_index].legacy_lower,
                                       groups[group_index].legacy_upper} );
    bool covered_legacy_valid = true;
    double covered_legacy_score = 0.0;
    double covered_legacy_mismatch = 0.0;
    PeakContinuum::OffsetType covered_legacy_type = PeakContinuum::OffsetType::Linear;
    for( const std::pair<double,double> &legacy : covered_legacy_intervals )
    {
      const size_t legacy_lower_index = candidate_index_for( legacy.first );
      const size_t legacy_upper_index = candidate_index_for( legacy.second );
      const IntervalScore legacy_score = ((legacy_lower_index < candidates.size())
          && (legacy_upper_index < candidates.size())
          && (legacy_upper_index > legacy_lower_index))
        ? score_interval( legacy_lower_index, legacy_upper_index ) : IntervalScore();
      if( !legacy_score.valid )
      {
        covered_legacy_valid = false;
        break;
      }
      covered_legacy_score += legacy_score.score;
      covered_legacy_mismatch += legacy_score.mismatch;
      covered_legacy_type = legacy_score.continuum_type;
    }
    if( covered_legacy_valid && !covered_legacy_intervals.empty() )
    {
      interval.legacy_score = covered_legacy_score;
      interval.legacy_normalized_continuum_mismatch
        = covered_legacy_mismatch / covered_legacy_intervals.size();
      interval.legacy_continuum_type = covered_legacy_type;
    }else
    {
      interval.legacy_score = std::numeric_limits<double>::quiet_NaN();
      interval.legacy_normalized_continuum_mismatch
        = std::numeric_limits<double>::quiet_NaN();
    }

    const bool split_legacy_roi = (interval.first_group == interval.last_group)
      && std::any_of( std::begin(groups), std::end(groups),
        [&interval, &groups]( const RoiBoundaryShadowGroup &group ) {
          const RoiBoundaryShadowGroup &covered = groups[interval.first_group];
          return (&group != &covered)
              && (std::fabs(group.legacy_lower - covered.legacy_lower) < 0.05)
              && (std::fabs(group.legacy_upper - covered.legacy_upper) < 0.05);
        } );
    if( split_legacy_roi )
      interval.reason = "split: joint SNIP continuum favors separate source groups";
  }

  result.valid = true;
  result.legacy_total_score = legacy_total_valid ? legacy_total_score
      : std::numeric_limits<double>::quiet_NaN();
  result.proposed_total_score = best.cost;
  result.intervals = std::move( best.intervals );
  return result;
}//optimize_roi_boundaries_shadow(...)


std::vector<RoiBoundaryShadowResult> take_roi_boundary_shadow_diagnostics()
{
  std::vector<RoiBoundaryShadowResult> answer;
  answer.swap( s_roi_boundary_shadow_diagnostics );
  return answer;
}//take_roi_boundary_shadow_diagnostics()


void record_roi_boundary_shadow_result( RoiBoundaryShadowResult result )
{
  s_roi_boundary_shadow_diagnostics.push_back( std::move(result) );
}//record_roi_boundary_shadow_result(...)


double LocalContinuumEstimate::sideband_asymmetry_z() const
{
  if( !valid )
    return 0.0;

  const double w_lo = lower_sideband_hi - lower_sideband_lo;
  const double w_up = upper_sideband_hi - upper_sideband_lo;
  if( (w_lo <= 0.0) || (w_up <= 0.0) )
    return 0.0;

  const double lo_dens = lower_sideband_counts / w_lo;
  const double up_dens = upper_sideband_counts / w_up;

  // Poisson variance of the density estimates, from the RAW (pre-subtraction) counts - the
  // subtracted prediction is treated as noise-free, which is slightly conservative.
  const double var = std::max( 1.0, lower_sideband_raw_counts ) / (w_lo * w_lo)
                     + std::max( 1.0, upper_sideband_raw_counts ) / (w_up * w_up);

  return (lo_dens - up_dens) / std::sqrt( var );
}//LocalContinuumEstimate::sideband_asymmetry_z


/** Expected counts over [x0,x1] from a set of Gaussian lines with the given total areas - the
 shared predicted-signal kernel for the sideband extension, clean-gap, keep-gate, and step-gate
 tests.  Lines with non-positive amplitude or invalid FWHM are skipped; result clamped >= 0. */
template <class FwhmFcn>
double predicted_gaussian_counts( const std::vector<double> &energies,
                                  const std::vector<double> &amplitudes,
                                  const FwhmFcn &fwhm_at_energy,
                                  const double x0, const double x1 )
{
  assert( energies.size() == amplitudes.size() );

  const double root_two = std::sqrt( 2.0 );
  double sum = 0.0;
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const double amp = amplitudes[i];
    if( amp <= 0.0 )
      continue;
    const double gamma_fwhm = fwhm_at_energy( energies[i] );
    if( !std::isfinite(gamma_fwhm) || (gamma_fwhm <= 0.0) )
      continue;
    const double gamma_sigma = gamma_fwhm / PhysicalUnits::fwhm_nsigma;
    const double t0 = (x0 - energies[i]) / (root_two * gamma_sigma);
    const double t1 = (x1 - energies[i]) / (root_two * gamma_sigma);
    sum += amp * 0.5 * (std::erf(t1) - std::erf(t0));
  }
  return std::max( 0.0, sum );
}//predicted_gaussian_counts


AdaptiveExtentResult extend_roi_by_sidebands(
  const std::vector<double> &gamma_energies,
  const std::vector<double> &gamma_amplitudes,
  const double effective_fwhm,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::function<double(double)> &fwhm_at_energy,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
  const double core_num_fwhm,
  const double extend_z,
  const double max_num_fwhm,
  const PeakDef::SkewType skew_type,
  const double lowest_energy,
  const double highest_energy )
{
  AdaptiveExtentResult result;

  assert( gamma_energies.size() == gamma_amplitudes.size() );

  if( gamma_energies.empty() || !std::isfinite(effective_fwhm) || (effective_fwhm <= 0.0) )
    return result;

  const auto minmax_gamma = std::minmax_element( std::begin(gamma_energies), std::end(gamma_energies) );
  const double e_lo = *minmax_gamma.first;
  const double e_hi = *minmax_gamma.second;
  const double fwhm = effective_fwhm;

  const double skew_extra = (skew_type == PeakDef::SkewType::NoSkew) ? 0.0 : sm_skew_low_side_extra_fwhm;

  const double limit_lo = std::max( lowest_energy, e_lo - (max_num_fwhm + skew_extra)*fwhm );
  const double limit_hi = std::min( highest_energy, e_hi + max_num_fwhm*fwhm );

  double cur_lo = std::clamp( e_lo - (core_num_fwhm + skew_extra)*fwhm, limit_lo, limit_hi );
  double cur_hi = std::clamp( e_hi + core_num_fwhm*fwhm, cur_lo, limit_hi );

  const double core_lo = cur_lo;
  const double core_hi = cur_hi;

  result.lower = cur_lo;
  result.upper = cur_hi;

  if( !foreground || !foreground->channel_energies() || (foreground->num_gamma_channels() < 8) )
    return result;  // no data to test against; keep the core extent

  // Expected counts from the cluster's own gammas over [x0,x1] (Gaussian model); this is the tail
  // leakage the extension test must not mistake for continuum inconsistency.
  const auto predicted_signal = [&]( const double x0, const double x1 ) -> double {
    return predicted_gaussian_counts( gamma_energies, gamma_amplitudes, fwhm_at_energy, x0, x1 );
  };//predicted_signal lambda

  // Interference veto: an unfit auto-search peak overlapping the candidate block (within half its
  // own FWHM) means structure that would contaminate the continuum - stop extending.
  const auto interfering_peak_near = [&unfit_auto_peaks]( const double b_lo, const double b_hi ) -> bool
  {
    for( const std::shared_ptr<const PeakDef> &p : unfit_auto_peaks )
    {
      if( !p )
        continue;
      const double half_w = 0.5 * p->fwhm();
      if( ((p->mean() + half_w) > b_lo) && ((p->mean() - half_w) < b_hi) )
        return true;
    }
    return false;
  };//interfering_peak_near lambda

  const double samp_num_fwhm = 0.5;  // sideband sample width for the continuum anchor

  // Per-block acceptance threshold, Bonferroni-calibrated so `extend_z` retains its meaning as
  // the FAMILY-wise z for a full side of extension: testing ~N independent ~0.375-FWHM blocks at
  // a fixed per-block z gives a family false-stop probability of ~N x alpha_block, so a genuinely
  // flat continuum would stop early ~28% of the time at extend_z=2 with N~7.  Splitting the
  // family alpha across the expected block count removes that data-length dependence.
  const double blocks_per_side = std::max( 1.0,
      std::ceil( (max_num_fwhm - core_num_fwhm) / sm_extend_block_fwhm ) );
  double block_z_thresh = extend_z;
  if( std::isfinite(extend_z) && (extend_z > 0.0) )
  {
    const boost::math::normal_distribution<double> gaus_dist;
    const double family_alpha = 2.0 * boost::math::cdf( gaus_dist, -extend_z );
    const double block_alpha = std::max( 1.0e-15, family_alpha / blocks_per_side );
    block_z_thresh = -boost::math::quantile( gaus_dist, 0.5 * block_alpha );
  }

  // Extend one side; dir = -1 for the low side, +1 for the high side.
  const auto extend_side = [&]( const int dir )
  {
    double z_sum = 0.0;
    size_t n_blocks = 0;

    while( true )
    {
      const double edge = (dir < 0) ? cur_lo : cur_hi;
      const double limit = (dir < 0) ? limit_lo : limit_hi;

      const double f_loc = fwhm_at_energy( edge );
      if( !std::isfinite(f_loc) || (f_loc <= 0.0) )
        break;

      // Block width: a fraction of the local FWHM, widened to at least 2 channels.
      const size_t edge_ch = foreground->find_gamma_channel( static_cast<float>(edge) );
      const double chan_w = std::max( 1.0e-6,
          static_cast<double>( foreground->gamma_channel_width( edge_ch ) ) );
      const double block_w = std::max( sm_extend_block_fwhm * f_loc, 2.0*chan_w );

      const double cand_lo = (dir < 0) ? (edge - block_w) : edge;
      const double cand_hi = (dir < 0) ? edge : (edge + block_w);

      if( (dir < 0) ? (cand_lo < limit) : (cand_hi > limit) )
        break;  // block would cross the extension cap - stop (no partial blocks)

      if( interfering_peak_near( cand_lo, cand_hi ) )
        break;

      // Continuum anchors: one sample-window at each edge, strictly INSIDE the accepted extent
      // (so the candidate block never helps anchor itself - that would make the test circular),
      // with the predicted signal of the cluster's own gammas subtracted (so strong-peak tails
      // at the core edges do not bias the line upward and spuriously stop extension).
      const double samp_w = samp_num_fwhm * f_loc;
      if( (cur_hi - cur_lo) <= 2.0*samp_w )
        break;

      const double lo_x0 = cur_lo, lo_x1 = cur_lo + samp_w;
      const double hi_x0 = cur_hi - samp_w, hi_x1 = cur_hi;

      const double lo_raw = foreground->gamma_integral( static_cast<float>(lo_x0), static_cast<float>(lo_x1) );
      const double hi_raw = foreground->gamma_integral( static_cast<float>(hi_x0), static_cast<float>(hi_x1) );
      const double lo_net = lo_raw - predicted_signal( lo_x0, lo_x1 );
      const double hi_net = hi_raw - predicted_signal( hi_x0, hi_x1 );

      const double lo_dens = std::max( 0.0, lo_net ) / samp_w;  // continuum counts per keV
      const double hi_dens = std::max( 0.0, hi_net ) / samp_w;
      const double lo_pos = 0.5*(lo_x0 + lo_x1);
      const double hi_pos = 0.5*(hi_x0 + hi_x1);

      // Linear density through the two anchors, integrated over the candidate block (the
      // integral of a linear density equals the midpoint density times the width).
      const double anchor_span = std::max( 1.0e-9, hi_pos - lo_pos );
      const double slope = (hi_dens - lo_dens) / anchor_span;
      const double cand_mid = 0.5*(cand_lo + cand_hi);
      const double cand_w = cand_hi - cand_lo;
      const double c_pred = std::max( 0.0, (lo_dens + slope*(cand_mid - lo_pos)) * cand_w );

      // Estimation variance of c_pred: it is a linear combination of the two noisy anchor
      // densities, c_pred = w*[(1-t)*lo_dens + t*hi_dens] with t the (possibly extrapolating,
      // |t|>1 leveraged) fractional position of the block midpoint between the anchors.  Using
      // the RAW anchor counts as the Poisson variance basis (the subtracted prediction is treated
      // as noise-free).  Omitting this term made the block z over-dispersed - extension stopped
      // earlier than the nominal threshold implied.
      const double t_lever = (cand_mid - lo_pos) / anchor_span;
      const double var_scale = (cand_w * cand_w) / (samp_w * samp_w);
      const double c_pred_var = var_scale * ( (1.0 - t_lever)*(1.0 - t_lever)*std::max( 1.0, lo_raw )
                                              + t_lever*t_lever*std::max( 1.0, hi_raw ) );

      const double s_pred = predicted_signal( cand_lo, cand_hi );
      const double d_obs = foreground->gamma_integral( static_cast<float>(cand_lo),
                                                       static_cast<float>(cand_hi) );

      const double z = (d_obs - c_pred - s_pred)
                       / std::sqrt( std::max( 1.0, c_pred + s_pred + c_pred_var ) );

      if( std::fabs(z) > block_z_thresh )
        break;

      // NOTE: successive block z's are positively correlated (they share the far anchor and the
      // fitted slope), so both this test and the drift guard below are approximate; the GA-tuned
      // extend_z absorbs the residual miscalibration.

      // Slow-drift guard: (mean z)^2 x n = z_sum^2 / n
      z_sum += z;
      const double n = static_cast<double>( n_blocks + 1 );
      if( ((z_sum * z_sum) / n) > sm_extend_drift_chi2 )
        break;  // candidate rejected: cumulative drift indicates curvature/structure

      n_blocks += 1;
      if( dir < 0 )
        cur_lo = cand_lo;
      else
        cur_hi = cand_hi;
    }//while( true )
  };//extend_side lambda

  extend_side( -1 );
  extend_side( +1 );

  result.lower = cur_lo;
  result.upper = cur_hi;
  result.sideband_lower_kev = std::max( 0.0, core_lo - cur_lo );
  result.sideband_upper_kev = std::max( 0.0, cur_hi - core_hi );

  return result;
}//extend_roi_by_sidebands


bool find_clean_gap_between(
  const std::vector<double> &left_energies,
  const std::vector<double> &left_amplitudes,
  const std::vector<double> &right_energies,
  const std::vector<double> &right_amplitudes,
  const double left_anchor,
  const double right_anchor,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::function<double(double)> &fwhm_at_energy,
  const double merge_tail_z,
  const double clean_gap_num_fwhm,
  double *clean_win_lo,
  double *clean_win_hi,
  const GlobalContinuumEstimate *global_continuum )
{
  if( clean_win_lo )
    *clean_win_lo = 0.0;
  if( clean_win_hi )
    *clean_win_hi = 0.0;

  assert( left_energies.size() == left_amplitudes.size() );
  assert( right_energies.size() == right_amplitudes.size() );

  if( !(right_anchor > left_anchor) )
    return false;

  const double mid_fwhm = fwhm_at_energy( 0.5*(left_anchor + right_anchor) );
  if( !std::isfinite(mid_fwhm) || (mid_fwhm <= 0.0) )
    return false;

  const double need = clean_gap_num_fwhm * mid_fwhm;
  if( (right_anchor - left_anchor) < need )
    return false;  // no room to anchor a continuum between the peaks - must merge

  // Predicted counts from BOTH groups' gammas over [x0,x1] (Gaussian model).
  const auto predicted_signal = [&]( const double x0, const double x1 ) -> double {
    return predicted_gaussian_counts( left_energies, left_amplitudes, fwhm_at_energy, x0, x1 )
           + predicted_gaussian_counts( right_energies, right_amplitudes, fwhm_at_energy, x0, x1 );
  };//predicted_signal lambda

  // Local continuum anchored outside the two anchor peaks, with the groups' own predicted tails
  // subtracted from the sideband samples: the sidebands sit only ~1-1.5 FWHM from the anchor
  // gammas (and other group members can sit right in them), so without subtraction c_est is
  // biased HIGH exactly when the peaks are strong and close - spuriously judging windows "clean"
  // and splitting ROIs that should merge.
  LocalContinuumEstimate cont;
  if( foreground )
    cont = estimate_local_continuum( foreground, left_anchor - mid_fwhm, right_anchor + mid_fwhm,
                                     mid_fwhm, 0.5, predicted_signal );

  const double step = 0.25 * mid_fwhm;
  bool found = false;
  double best_pred = std::numeric_limits<double>::max();

  for( double win_lo = left_anchor; (win_lo + need) <= (right_anchor + 1.0e-9); win_lo += step )
  {
    const double win_hi = win_lo + need;

    // Contamination-vs-noise at WINDOW level: could a continuum anchored on this whole window be
    // biased by the groups' tails beyond its own Poisson uncertainty?  (This was formerly tested
    // per ~0.25-FWHM block, but a per-block z understates the window-level contamination by
    // ~sqrt(block/window) - S scales with width while the noise scales with its square root -
    // making the old test anti-conservative, i.e. biased toward splitting.)
    const double s_pred = predicted_signal( win_lo, win_hi );

    // R1 step 2: anchor the clean-gap continuum on the single global SNIP continuum when available;
    // fall back to the tail-subtracted local estimate, then to gross counts.  The predicted-signal
    // tail term (s_pred) is unchanged.
    double c_est = 0.0;
    if( global_continuum && global_continuum->valid() )
      c_est = global_continuum->integral( win_lo, win_hi );
    else if( cont.valid )
      c_est = cont.integral( win_lo, win_hi );
    else if( foreground )
      c_est = foreground->gamma_integral( static_cast<float>(win_lo), static_cast<float>(win_hi) );

    const double z_window = s_pred / std::sqrt( std::max( 1.0, c_est ) );
    const bool clean = (z_window < merge_tail_z);

    if( clean && (s_pred < best_pred) )
    {
      found = true;
      best_pred = s_pred;
      if( clean_win_lo )
        *clean_win_lo = win_lo;
      if( clean_win_hi )
        *clean_win_hi = win_hi;
    }
  }//for( slide window across the gap )

  return found;
}//find_clean_gap_between


namespace
{
struct MeasuredRoiModelFit
{
  bool valid = false;
  double poisson_deviance = std::numeric_limits<double>::max();
  double aicc = std::numeric_limits<double>::max();
  size_t num_parameters = 0;
  PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;
  std::vector<double> amplitudes;
  std::vector<double> amplitude_uncertainties;
};


/** Fit one of the small structural hypotheses used below on an immutable foreground channel
 domain.  Peak means/widths are fixed; their areas are either linear fit parameters (`means`) or
 fixed by `fixed_peaks`.  The returned likelihood is always evaluated from the measured counts,
 even though the linear solve uses the usual Poisson-weighted least-squares initialization. */
MeasuredRoiModelFit fit_measured_roi_model(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const size_t first_channel,
  const size_t last_channel,
  const double reference_energy,
  const std::vector<double> &means,
  const std::vector<double> &sigmas,
  const std::vector<PeakDef> &fixed_peaks,
  const double aicc_penalty,
  std::vector<MeasuredRoiModelFit> *all_fits = nullptr )
{
  MeasuredRoiModelFit best;
  if( all_fits )
    all_fits->clear();
  if( !foreground || !foreground->channel_energies()
      || (first_channel >= foreground->num_gamma_channels())
      || (last_channel < first_channel) || !(aicc_penalty > 0.0)
      || (means.size() != sigmas.size()) )
    return best;

  const size_t nbin = 1 + last_channel - first_channel;
  const std::vector<float> &channel_energies = *foreground->channel_energies();
  if( channel_energies.size() <= (first_channel + nbin) )
    return best;

  std::vector<float> counts( nbin );
  std::vector<float> variances( nbin );
  for( size_t index = 0; index < nbin; ++index )
  {
    counts[index] = foreground->gamma_channel_content( first_channel + index );
    variances[index] = std::max( 1.0f, counts[index] );
  }

  const PeakContinuum::OffsetType families[] = {
    PeakContinuum::OffsetType::Linear,
    PeakContinuum::OffsetType::Quadratic,
    PeakContinuum::OffsetType::FlatStep,
    PeakContinuum::OffsetType::LinearStep
  };
  for( const PeakContinuum::OffsetType family : families )
  {
    const size_t num_parameters = PeakContinuum::num_parameters( family ) + means.size();
    if( nbin <= (num_parameters + 1) )
      continue;

    std::vector<double> amplitudes, coefficients, amplitude_uncertainties;
    std::vector<double> coefficient_uncertainties;
    std::vector<double> predictions( nbin, 0.0 );
    try
    {
      static_cast<void>( PeakFit::fit_amp_and_offset_imp<PeakDef,double>(
          &channel_energies[first_channel], counts.data(), variances.data(), nbin,
          family, 0.0, reference_energy, means, sigmas, fixed_peaks,
          PeakDef::SkewType::NoSkew, nullptr, amplitudes, coefficients,
          amplitude_uncertainties, coefficient_uncertainties, predictions.data() ) );
    }catch( const std::exception & )
    {
      continue;
    }
    if( (coefficients.size() != PeakContinuum::num_parameters(family))
        || (amplitudes.size() != means.size())
        || std::any_of( std::begin(amplitudes), std::end(amplitudes),
             []( const double amplitude ) {
               return !std::isfinite(amplitude) || (amplitude < 0.0);
             } ) )
      continue;

    bool physically_valid = true;
    double deviance = 0.0;
    for( size_t index = 0; index < nbin; ++index )
    {
      const double observed = std::max( 0.0, static_cast<double>(counts[index]) );
      double predicted = predictions[index];
      if( !std::isfinite(predicted) || (predicted < -1.0e-6) )
      {
        physically_valid = false;
        break;
      }
      predicted = std::max( 1.0e-9, predicted );
      deviance += (observed > 0.0)
        ? 2.0*(predicted - observed + observed*std::log(observed/predicted))
        : 2.0*predicted;
    }
    if( !physically_valid )
      continue;

    const double score = data_only_aicc( deviance, nbin, num_parameters, aicc_penalty );
    MeasuredRoiModelFit candidate;
    candidate.valid = true;
    candidate.poisson_deviance = deviance;
    candidate.aicc = score;
    candidate.num_parameters = num_parameters;
    candidate.continuum_type = family;
    candidate.amplitudes = amplitudes;
    candidate.amplitude_uncertainties = amplitude_uncertainties;
    if( all_fits )
      all_fits->push_back( candidate );
    if( score < best.aicc )
    {
      best = std::move( candidate );
    }
  }
  return best;
}//fit_measured_roi_model


/** Fit one nonnegative scale multiplying the exact fixed-ratio source-line template.  Each trial
 scale re-fits the production continuum families on the same channels; AICc counts the shared scale
 once, independent of the number of source lines. */
MeasuredRoiModelFit fit_measured_roi_scaled_source_model(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const size_t first_channel,
  const size_t last_channel,
  const double reference_energy,
  const std::vector<double> &means,
  const std::vector<double> &sigmas,
  const std::vector<double> &weights,
  const double aicc_penalty )
{
  MeasuredRoiModelFit best;
  if( !foreground || (means.size() != sigmas.size()) || (means.size() != weights.size())
      || means.empty() || (last_channel < first_channel) )
    return best;

  const double weight_sum = std::accumulate( std::begin(weights), std::end(weights), 0.0,
      []( const double sum, const double weight ) {
        return sum + ((std::isfinite(weight) && (weight > 0.0)) ? weight : 0.0);
      } );
  if( !(weight_sum > 0.0) )
    return best;
  const size_t num_channels = 1 + last_channel - first_channel;
  const double gross_counts = foreground->gamma_integral(
      foreground->gamma_channel_lower(first_channel),
      foreground->gamma_channel_upper(last_channel) );
  const double max_area = std::max( 1.0, 4.0*gross_counts );

  const auto evaluate = [&]( const double total_area ) {
    MeasuredRoiModelFit trial_best;
    std::vector<PeakDef> fixed_peaks;
    fixed_peaks.reserve( means.size() );
    for( size_t index = 0; index < means.size(); ++index )
    {
      if( !std::isfinite(means[index]) || !std::isfinite(sigmas[index])
          || !(sigmas[index] > 0.0) || !std::isfinite(weights[index])
          || !(weights[index] > 0.0) )
        continue;
      fixed_peaks.emplace_back( means[index], sigmas[index],
          total_area*weights[index]/weight_sum );
    }
    if( fixed_peaks.empty() )
      return trial_best;
    std::vector<MeasuredRoiModelFit> family_fits;
    const std::vector<double> no_means, no_sigmas;
    static_cast<void>( fit_measured_roi_model( foreground, first_channel, last_channel,
        reference_energy, no_means, no_sigmas, fixed_peaks, aicc_penalty, &family_fits ) );
    for( MeasuredRoiModelFit &fit : family_fits )
    {
      ++fit.num_parameters;  // the one common source-template scale
      fit.aicc = data_only_aicc( fit.poisson_deviance, num_channels,
          fit.num_parameters, aicc_penalty );
      fit.amplitudes = { total_area };
      if( fit.aicc < trial_best.aicc )
        trial_best = std::move( fit );
    }
    return trial_best;
  };

  // Coarse bracketing followed by a bounded golden-section refinement.  The best continuum family
  // may change with scale, so retain the best scored trial rather than assuming differentiability.
  constexpr size_t num_grid_steps = 16;
  size_t best_grid = 0;
  for( size_t index = 0; index <= num_grid_steps; ++index )
  {
    const double area = max_area*static_cast<double>(index)
        / static_cast<double>(num_grid_steps);
    MeasuredRoiModelFit candidate = evaluate( area );
    if( candidate.aicc < best.aicc )
    {
      best = std::move( candidate );
      best_grid = index;
    }
  }
  double lower = max_area*static_cast<double>(best_grid > 0 ? best_grid - 1 : 0)
      / static_cast<double>(num_grid_steps);
  double upper = max_area*static_cast<double>(std::min(num_grid_steps, best_grid + 1))
      / static_cast<double>(num_grid_steps);
  constexpr double inv_phi = 0.6180339887498948482;
  double left = upper - inv_phi*(upper - lower);
  double right = lower + inv_phi*(upper - lower);
  MeasuredRoiModelFit left_fit = evaluate( left );
  MeasuredRoiModelFit right_fit = evaluate( right );
  for( size_t iteration = 0; iteration < 14; ++iteration )
  {
    if( left_fit.aicc < right_fit.aicc )
    {
      upper = right;
      right = left;
      right_fit = left_fit;
      left = upper - inv_phi*(upper - lower);
      left_fit = evaluate( left );
    }else
    {
      lower = left;
      left = right;
      left_fit = right_fit;
      right = lower + inv_phi*(upper - lower);
      right_fit = evaluate( right );
    }
  }
  if( left_fit.aicc < best.aicc )
    best = std::move( left_fit );
  if( right_fit.aicc < best.aicc )
    best = std::move( right_fit );
  return best;
}//fit_measured_roi_scaled_source_model
}//namespace


SourceClusterEvidenceResult evaluate_source_cluster_evidence(
  const std::vector<double> &source_energies,
  const std::vector<double> &source_areas,
  const double lower_energy,
  const double upper_energy,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::function<double(double)> &fwhm_at_energy,
  const std::vector<std::shared_ptr<const PeakDef>> &found_peaks,
  const double significance_z,
  const double source_core_num_fwhm,
  const double aicc_penalty )
{
  SourceClusterEvidenceResult result;
  if( !foreground || !foreground->energy_calibration()
      || !foreground->energy_calibration()->valid() || !fwhm_at_energy
      || source_energies.empty() || (source_energies.size() != source_areas.size())
      || !(upper_energy > lower_energy) || !(significance_z > 0.0)
      || !(source_core_num_fwhm > 0.0)
      || !(aicc_penalty > 0.0) )
  {
    result.reason = "invalid local source-evidence input";
    return result;
  }

  const size_t first_channel = foreground->find_gamma_channel(
      static_cast<float>(lower_energy) );
  const size_t last_channel = std::min( foreground->find_gamma_channel(
      static_cast<float>(upper_energy) ), foreground->num_gamma_channels() - 1 );
  if( last_channel <= (first_channel + 4) )
  {
    result.reason = "too few channels for local source-evidence comparison";
    return result;
  }
  const double reference_energy = 0.5*(lower_energy + upper_energy);
  const std::vector<double> no_means, no_sigmas;
  const std::vector<PeakDef> no_fixed_peaks;
  const MeasuredRoiModelFit null_fit = fit_measured_roi_model(
      foreground, first_channel, last_channel, reference_energy,
      no_means, no_sigmas, no_fixed_peaks, aicc_penalty );

  // Hs has one local scale parameter multiplying the exact requested-source line mixture.  Its
  // line energies, widths, and relative areas are fixed, while the shared scale is allowed to
  // correct an imperfect provisional global activity.
  std::vector<double> source_means, source_sigmas, source_weights;
  for( size_t index = 0; index < source_energies.size(); ++index )
  {
    const double fwhm = fwhm_at_energy( source_energies[index] );
    if( !std::isfinite(fwhm) || !(fwhm > 0.0)
        || !std::isfinite(source_areas[index]) || !(source_areas[index] > 0.0) )
      continue;
    source_means.push_back( source_energies[index] );
    source_sigmas.push_back( fwhm / PhysicalUnits::fwhm_nsigma );
    source_weights.push_back( source_areas[index] );
  }
  if( source_means.empty() )
  {
    result.reason = "source hypothesis has no finite positive prediction";
    return result;
  }
  const MeasuredRoiModelFit source_fit = fit_measured_roi_scaled_source_model(
      foreground, first_channel, last_channel, reference_energy,
      source_means, source_sigmas, source_weights, aicc_penalty );

  result.null_aicc = null_fit.valid ? null_fit.aicc
                                    : std::numeric_limits<double>::quiet_NaN();
  result.source_aicc = source_fit.valid ? source_fit.aicc
                                        : std::numeric_limits<double>::quiet_NaN();
  if( !null_fit.valid )
  {
    result.reason = "local continuum/source comparison was ill-conditioned";
    return result;
  }
  if( !source_fit.valid )
  {
    result.decision = SourceClusterEvidenceDecision::RejectContinuumOnly;
    result.reason = "no nonnegative source-tied amplitude improves the continuum model";
    return result;
  }

  const double source_delta = null_fit.poisson_deviance - source_fit.poisson_deviance;
  result.source_likelihood_z = std::sqrt( std::max(0.0, source_delta) );

  std::vector<double> free_means, free_sigmas;
  for( const std::shared_ptr<const PeakDef> &peak : found_peaks )
  {
    if( !peak || !peak->gausPeak() || (peak->mean() < lower_energy)
        || (peak->mean() > upper_energy) )
      continue;
    const bool matches_source_core = std::any_of( std::begin(source_means),
        std::end(source_means), [&peak, &fwhm_at_energy, source_core_num_fwhm](
            const double source_energy ) {
          const double source_fwhm = fwhm_at_energy( source_energy );
          return std::isfinite(source_fwhm) && (source_fwhm > 0.0)
              && (std::fabs(peak->mean() - source_energy)
                  < (source_core_num_fwhm*source_fwhm));
        } );
    if( matches_source_core )
      continue;  // a requested line cannot compete against its own fixed-ratio Hs mixture as Hf
    const double amplitude = peak->amplitude();
    const double uncertainty = peak->amplitudeUncert();
    const double found_z = (uncertainty > 0.0) ? (amplitude / uncertainty)
        : ((amplitude > 0.0) ? std::sqrt(amplitude) : 0.0);
    if( found_z < significance_z )
      continue;
    const double sigma = (peak->sigma() > 0.0)
        ? peak->sigma() : (fwhm_at_energy(peak->mean()) / PhysicalUnits::fwhm_nsigma);
    if( !std::isfinite(sigma) || !(sigma > 0.0) )
      continue;
    const bool distinct = std::none_of( std::begin(free_means), std::end(free_means),
        [&peak, &fwhm_at_energy]( const double accepted_mean ) {
          const double fwhm = fwhm_at_energy( 0.5*(accepted_mean + peak->mean()) );
          return std::isfinite(fwhm) && (fwhm > 0.0)
              && (std::fabs(accepted_mean - peak->mean()) < fwhm);
        } );
    if( distinct )
    {
      free_means.push_back( peak->mean() );
      free_sigmas.push_back( sigma );
    }
  }
  const MeasuredRoiModelFit best_free_fit = free_means.empty()
      ? MeasuredRoiModelFit()
      : fit_measured_roi_model( foreground, first_channel, last_channel, reference_energy,
          free_means, free_sigmas, no_fixed_peaks, aicc_penalty );
  if( !free_means.empty() )
    result.free_feature_energy = free_means.front();
  result.free_feature_aicc = best_free_fit.valid ? best_free_fit.aicc
                                                  : std::numeric_limits<double>::quiet_NaN();

  // AICc is the actual transactional selector here: it already requires the one additional source
  // scale parameter to earn its place on the identical channels.  Requiring the peak-only z to
  // ALSO exceed the later observable threshold double-penalizes blended source groups.  In dense
  // forests a real shoulder can win the parameter-penalized comparison while its isolated z is
  // just below that downstream gate (the U235 90-keV complex is a representative case).  Such a
  // group is evidence, not a demonstrably empty bridge, so retain it and leave final significance
  // pruning to the existing observable stage.
  const bool source_beats_null = (source_fit.aicc < null_fit.aicc);
  if( !source_beats_null )
  {
    result.decision = SourceClusterEvidenceDecision::RejectContinuumOnly;
    result.reason = "continuum-only model explains the provisional source group";
    return result;
  }
  if( best_free_fit.valid && (best_free_fit.aicc < source_fit.aicc) )
  {
    result.decision = SourceClusterEvidenceDecision::RejectFreeFeature;
    result.reason = "free significant peak beats the source-tied prediction";
    return result;
  }

  result.decision = SourceClusterEvidenceDecision::RetainSource;
  result.reason = "source-tied prediction beats continuum and is not worse than a free feature";
  return result;
}//evaluate_source_cluster_evidence


AutomaticRoiPolicyResult evaluate_automatic_roi_boundary(
  const AutomaticRoiGroup &left,
  const AutomaticRoiGroup &right,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const GlobalContinuumEstimate *global_continuum,
  const std::function<double(double)> &fwhm_at_energy,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
  const AutomaticRoiPolicySettings &settings )
{
  AutomaticRoiPolicyResult result;
  AutomaticRoiDecisionDiagnostic &diag = result.diagnostic;
  diag.stage = settings.stage;
  diag.left_lower = left.lower;
  diag.left_upper = left.upper;
  diag.right_lower = right.lower;
  diag.right_upper = right.upper;
  diag.calibration_num_channels = foreground ? foreground->num_gamma_channels() : 0;
  diag.used_global_continuum = global_continuum && global_continuum->valid()
      && (global_continuum->foreground == foreground);

  if( left.protected_geometry || right.protected_geometry )
  {
    result.decision = diag.decision = AutomaticRoiDecision::ProtectedGeometry;
    diag.reason = "protected user/mixed ROI geometry";
    return result;
  }

  if( !foreground || !foreground->energy_calibration()
      || !foreground->energy_calibration()->valid() || !fwhm_at_energy
      || left.peak_energies.empty() || right.peak_energies.empty() )
  {
    result.decision = diag.decision = AutomaticRoiDecision::KeepSeparate;
    diag.reason = "insufficient calibrated evidence; conservative separation";
    return result;
  }

  const double left_anchor = *std::max_element(
      std::begin(left.peak_energies), std::end(left.peak_energies) );
  const double right_anchor = *std::min_element(
      std::begin(right.peak_energies), std::end(right.peak_energies) );
  const double midpoint = 0.5*(left_anchor + right_anchor);
  const double fwhm = fwhm_at_energy( midpoint );
  if( std::isfinite(fwhm) && (fwhm > 0.0) )
  {
    diag.combined_width_fwhm = (std::max(left.upper, right.upper)
        - std::min(left.lower, right.lower)) / fwhm;
    const double width_ratio = (settings.max_width_fwhm > 0.0)
        ? (diag.combined_width_fwhm / settings.max_width_fwhm) : 0.0;
    const double excess_width = std::max( 0.0, width_ratio - 1.0 );
    diag.width_pressure = settings.continuum_aicc_penalty * excess_width * excess_width
        * std::max<size_t>( 1, left.joined_groups + right.joined_groups - 1 );
  }
  if( !std::isfinite(fwhm) || !(fwhm > 0.0) )
  {
    result.decision = diag.decision = AutomaticRoiDecision::KeepSeparate;
    diag.reason = "invalid local resolution; conservative separation";
    return result;
  }
  if( !(right_anchor > left_anchor) )
  {
    result.decision = diag.decision = (diag.width_pressure > 0.0)
        ? AutomaticRoiDecision::MergeInseparableWide
        : AutomaticRoiDecision::MergeInseparable;
    diag.reason = "modeled peak cores overlap";
    return result;
  }
  const bool modeled_cores_overlap = (right_anchor - left_anchor)
      <= (2.0 * settings.peak_core_num_fwhm * fwhm);
  if( modeled_cores_overlap
      && !(settings.allow_overwide_overlap_partition && (diag.width_pressure > 0.0)) )
  {
    result.decision = diag.decision = (diag.width_pressure > 0.0)
        ? AutomaticRoiDecision::MergeInseparableWide
        : AutomaticRoiDecision::MergeInseparable;
    diag.reason = "modeled peak cores overlap";
    return result;
  }

  diag.separation_fwhm = (right_anchor - left_anchor) / fwhm;

  std::vector<double> left_areas = left.peak_areas;
  std::vector<double> right_areas = right.peak_areas;
  left_areas.resize( left.peak_energies.size(), 0.0 );
  right_areas.resize( right.peak_energies.size(), 0.0 );

  double clean_lo = 0.0;
  double clean_hi = 0.0;
  const bool have_clean_window = find_clean_gap_between(
      left.peak_energies, left_areas, right.peak_energies, right_areas,
      left_anchor, right_anchor, foreground, fwhm_at_energy, settings.merge_tail_z,
      settings.merge_clean_gap_fwhm, &clean_lo, &clean_hi,
      diag.used_global_continuum ? global_continuum : nullptr );
  result.boundary_energy = diag.boundary_energy = have_clean_window
      ? 0.5*(clean_lo + clean_hi) : midpoint;
  diag.boundary_channel = foreground->find_gamma_channel(
      static_cast<float>(result.boundary_energy) );

  // A boundary may not pass through a credible unmodeled +/-1-FWHM core.  Such a feature needs
  // its own excluded/transactional neighborhood; joining across it would hide that fact.
  for( const std::shared_ptr<const PeakDef> &peak : unfit_auto_peaks )
  {
    if( !peak || !peak->gausPeak() )
      continue;
    const double peak_fwhm = fwhm_at_energy( peak->mean() );
    if( !std::isfinite(peak_fwhm) || !(peak_fwhm > 0.0) )
      continue;
    if( (peak->mean() + peak_fwhm) > left_anchor
        && (peak->mean() - peak_fwhm) < right_anchor )
    {
      diag.unmodeled_core_blocked = true;
      result.exclusion_lower = peak->mean() - peak_fwhm;
      result.exclusion_upper = peak->mean() + peak_fwhm;
      result.decision = diag.decision = AutomaticRoiDecision::UnmodeledFeatureBlocked;
      diag.reason = "unmodeled peak core occupies the proposed gap";
      return result;
    }
  }

  const double valley_lo = have_clean_window ? clean_lo
      : std::max( left_anchor, midpoint - 0.5*settings.merge_clean_gap_fwhm*fwhm );
  const double valley_hi = have_clean_window ? clean_hi
      : std::min( right_anchor, midpoint + 0.5*settings.merge_clean_gap_fwhm*fwhm );
  if( valley_hi > valley_lo )
  {
    diag.observed_valley_counts = foreground->gamma_integral(
        static_cast<float>(valley_lo), static_cast<float>(valley_hi) );
    diag.snip_valley_counts = diag.used_global_continuum
        ? global_continuum->integral( valley_lo, valley_hi )
        : diag.observed_valley_counts;
    diag.modeled_tail_counts = predicted_gaussian_counts(
        left.peak_energies, left_areas, fwhm_at_energy, valley_lo, valley_hi )
      + predicted_gaussian_counts(
        right.peak_energies, right_areas, fwhm_at_energy, valley_lo, valley_hi );
    const double continuum_variance = diag.used_global_continuum
        ? global_continuum->integral_variance( valley_lo, valley_hi )
        : diag.observed_valley_counts;
    const double noise = std::sqrt( std::max( 1.0,
        continuum_variance + diag.modeled_tail_counts ) );
    diag.modeled_tail_significance = diag.modeled_tail_counts / noise;
    diag.unexplained_excess_significance = std::max( 0.0,
        diag.observed_valley_counts - diag.snip_valley_counts - diag.modeled_tail_counts ) / noise;
    diag.snip_mismatch_significance = std::fabs(
        diag.observed_valley_counts - diag.snip_valley_counts ) / noise;
  }

  const auto count_sideband_channels = [&]( const AutomaticRoiGroup &group,
                                             const double bound_lo,
                                             const double bound_hi ) -> size_t {
    size_t count = 0;
    const size_t first = foreground->find_gamma_channel( static_cast<float>(bound_lo) );
    const size_t last = std::min( foreground->find_gamma_channel(
        static_cast<float>(bound_hi) ), foreground->num_gamma_channels() - 1 );
    for( size_t channel = first; channel <= last; ++channel )
    {
      const double channel_lo = foreground->gamma_channel_lower( channel );
      const double channel_hi = foreground->gamma_channel_upper( channel );
      bool core = false;
      for( const double energy : group.peak_energies )
      {
        const double peak_fwhm = fwhm_at_energy( energy );
        core = core || ((channel_hi > (energy - settings.peak_core_num_fwhm*peak_fwhm))
            && (channel_lo < (energy + settings.peak_core_num_fwhm*peak_fwhm)));
      }
      if( !core )
        ++count;
      if( channel == last )
        break;
    }
    return count;
  };
  diag.left_sideband_channels = count_sideband_channels(
      left, left.lower, result.boundary_energy );
  diag.right_sideband_channels = count_sideband_channels(
      right, result.boundary_energy, right.upper );
  diag.sidebands_adequate = (diag.left_sideband_channels >= 8)
      && (diag.right_sideband_channels >= 8);

  struct ContinuumFit
  {
    bool valid = false;
    double mismatch = std::numeric_limits<double>::max();
    size_t num_parameters = 0;
  };
  const auto continuum_fits = [&]( const size_t first, const size_t last ) {
    std::vector<ContinuumFit> fits;
    if( !diag.used_global_continuum || (last < first) )
      return fits;
    const size_t nbin = last - first + 1;
    const std::vector<float> &channel_energies = *foreground->channel_energies();
    std::vector<float> snip_counts( nbin );
    std::vector<float> raw_variances( nbin );
    for( size_t index = 0; index < nbin; ++index )
    {
      snip_counts[index] = global_continuum->snip->gamma_channel_content( first + index );
      raw_variances[index] = static_cast<float>( std::max( 1.0,
          static_cast<double>(foreground->gamma_channel_content(first + index)) ) );
    }
    const PeakContinuum::OffsetType families[] = {
      PeakContinuum::OffsetType::Linear, PeakContinuum::OffsetType::Quadratic,
      PeakContinuum::OffsetType::FlatStep, PeakContinuum::OffsetType::LinearStep
    };
    const std::vector<double> no_means, no_sigmas;
    const std::vector<PeakDef> no_fixed_peaks;
    for( const PeakContinuum::OffsetType family : families )
    {
      const size_t k = PeakContinuum::num_parameters( family );
      if( nbin <= (k + 1) )
        continue;
      std::vector<double> amplitudes, coefficients, amplitude_uncerts, coefficient_uncerts;
      std::vector<double> predictions( nbin, 0.0 );
      try
      {
        static_cast<void>( PeakFit::fit_amp_and_offset_imp<PeakDef,double>(
            &channel_energies[first], snip_counts.data(), raw_variances.data(), nbin,
            family, 0.0, midpoint, no_means, no_sigmas, no_fixed_peaks,
            PeakDef::SkewType::NoSkew, nullptr, amplitudes, coefficients,
            amplitude_uncerts, coefficient_uncerts, predictions.data() ) );
      }catch( const std::exception & )
      {
        continue;
      }
      ContinuumFit fit;
      fit.valid = (coefficients.size() == k);
      fit.num_parameters = k;
      fit.mismatch = 0.0;
      for( size_t index = 0; fit.valid && (index < nbin); ++index )
      {
        const double prediction = predictions[index];
        if( !std::isfinite(prediction) || (prediction < -1.0e-6) )
        {
          fit.valid = false;
          break;
        }
        const double residual = snip_counts[index] - prediction;
        fit.mismatch += residual*residual / raw_variances[index];
      }
      if( fit.valid )
        fits.push_back( fit );
    }
    return fits;
  };
  const size_t union_first = foreground->find_gamma_channel(
      static_cast<float>(std::min(left.lower, right.lower)) );
  const size_t union_last = std::min( foreground->find_gamma_channel(
      static_cast<float>(std::max(left.upper, right.upper)) ),
      foreground->num_gamma_channels() - 1 );
  const size_t n = union_last - union_first + 1;
  // AICc is a sum over channel observations.  Treat the configured soft-width value as a
  // per-channel structural pressure too; adding the former dimensionless component-level value to
  // a channel-summed likelihood made width pressure vanish as ROIs gained channels.
  diag.width_pressure *= n;
  const auto aicc = [&]( const double mismatch, const size_t k ) -> double {
    if( n <= (k + 1) )
      return std::numeric_limits<double>::max();
    return mismatch + settings.continuum_aicc_penalty*k
        + settings.continuum_aicc_penalty*k*(k + 1.0)/(n - k - 1.0);
  };
  diag.one_roi_aicc = std::numeric_limits<double>::max();
  for( const ContinuumFit &fit : continuum_fits(union_first, union_last) )
    diag.one_roi_aicc = std::min( diag.one_roi_aicc,
        aicc(fit.mismatch, fit.num_parameters) + diag.width_pressure );
  diag.two_roi_aicc = std::numeric_limits<double>::max();
  const size_t first_candidate_channel = foreground->find_gamma_channel(
      static_cast<float>(have_clean_window ? clean_lo : result.boundary_energy) );
  const size_t last_candidate_channel = foreground->find_gamma_channel(
      static_cast<float>(have_clean_window ? clean_hi : result.boundary_energy) );
  for( size_t split_channel = first_candidate_channel;
       (split_channel <= last_candidate_channel) && (split_channel < union_last);
       ++split_channel )
  {
    if( split_channel < union_first )
      continue;
    const double candidate_energy = 0.5*(foreground->gamma_channel_lower(split_channel)
        + foreground->gamma_channel_upper(split_channel));
    if( (count_sideband_channels(left, left.lower,
             foreground->gamma_channel_lower(split_channel)) < 8)
        || (count_sideband_channels(right,
             foreground->gamma_channel_upper(split_channel), right.upper) < 8) )
      continue;
    const std::vector<ContinuumFit> left_fits = continuum_fits( union_first, split_channel );
    const std::vector<ContinuumFit> right_fits = continuum_fits( split_channel + 1, union_last );
    for( const ContinuumFit &left_fit : left_fits )
      for( const ContinuumFit &right_fit : right_fits )
      {
        const double candidate_aicc = aicc(
            left_fit.mismatch + right_fit.mismatch,
            left_fit.num_parameters + right_fit.num_parameters );
        if( candidate_aicc < diag.two_roi_aicc )
        {
          diag.two_roi_aicc = candidate_aicc;
          result.boundary_energy = diag.boundary_energy = candidate_energy;
          diag.boundary_channel = split_channel;
          diag.left_sideband_channels = count_sideband_channels(
              left, left.lower, foreground->gamma_channel_lower(split_channel) );
          diag.right_sideband_channels = count_sideband_channels(
              right, foreground->gamma_channel_upper(split_channel), right.upper );
          diag.sidebands_adequate = (diag.left_sideband_channels >= 8)
              && (diag.right_sideband_channels >= 8);
        }
      }
  }

  // This is evidence that no statistically material modeled tail or unexplained peak-like excess
  // connects the children; it does not require a morphological local minimum in the raw counts.
  const bool credible_unbridged_boundary = have_clean_window && diag.sidebands_adequate
      && (diag.modeled_tail_significance <= settings.merge_tail_z)
      && (diag.unexplained_excess_significance <= settings.merge_tail_z);
  const double unavailable_score = std::numeric_limits<double>::max();
  const bool one_score_valid = (diag.one_roi_aicc < unavailable_score);
  const bool two_score_valid = (diag.two_roi_aicc < unavailable_score);
  const bool two_no_worse = two_score_valid
      && (!one_score_valid || (diag.two_roi_aicc <= diag.one_roi_aicc));
  if( credible_unbridged_boundary && two_no_worse )
  {
    result.decision = diag.decision = AutomaticRoiDecision::KeepSeparate;
    diag.reason = "no significant peak bridge and two-continuum AICc is no worse";
  }else
  {
    const bool wide = (diag.width_pressure > 0.0);
    result.decision = diag.decision = wide
        ? AutomaticRoiDecision::MergeInseparableWide
        : AutomaticRoiDecision::MergeInseparable;
    if( !credible_unbridged_boundary )
      diag.reason = wide ? "no defensible peak-unbridged boundary; inseparable wide group"
                         : "no defensible peak-unbridged boundary; statistically inseparable";
    else
      diag.reason = wide ? "one-continuum AICc overcomes width pressure; inseparable wide group"
                         : "one-continuum AICc favors a shared continuum";
  }
  return result;
}//evaluate_automatic_roi_boundary


//=============================================================================================
// Atom-safe automatic ROI partition layer.
//
// evaluate_automatic_roi_boundary (above) remains the pure decision oracle.  The functions below
// own all geometry materialization for policy-mode ROI split/combine and carry stable atoms with
// the geometry so ownership is never reconstructed by energy containment.  Every operation is
// validated to preserve the admitted-atom multiset exactly-once before it can replace incumbent
// geometry; when no core-safe partition exists the incumbent is retained (merge or unchanged),
// never a dropped side.
//=============================================================================================

uint64_t next_roi_atom_id()
{
  static std::atomic<uint64_t> s_counter{ 1 };
  return s_counter.fetch_add( 1, std::memory_order_relaxed );
}


// Clamp an energy to a valid channel index for `fg`.
static size_t atomlayer_channel_for_energy(
    const std::shared_ptr<const SpecUtils::Measurement> &fg, const double energy )
{
  const size_t nchan = fg->num_gamma_channels();
  if( nchan == 0 )
    return 0;
  size_t ch = fg->find_gamma_channel( static_cast<float>( energy ) );
  if( ch >= nchan )
    ch = nchan - 1;
  return ch;
}


// Core half-width (keV) of an atom: peak_core_num_fwhm * FWHM(energy).  NaN if resolution invalid.
static double atomlayer_core_halfwidth( const double energy,
    const std::function<double(double)> &fwhm_at, const double core_num_fwhm )
{
  const double f = fwhm_at( energy );
  if( !std::isfinite(f) || (f <= 0.0) )
    return std::numeric_limits<double>::quiet_NaN();
  return core_num_fwhm * f;
}


// True iff the channel band [gap_first, gap_last] can be excluded from both children without
// cutting any atom core: every atom's +/- core lies wholly below the band's lower edge or wholly
// above its upper edge, and no atom energy lies within the band.  A core that cannot be evaluated
// (invalid resolution) makes the band unsafe.
static bool atomlayer_gap_core_safe( const size_t gap_first, const size_t gap_last,
    const std::vector<RoiAtom> &atoms,
    const std::shared_ptr<const SpecUtils::Measurement> &fg,
    const std::function<double(double)> &fwhm_at, const double core_num_fwhm )
{
  if( gap_last < gap_first )
    return false;
  const double glo = fg->gamma_channel_lower( gap_first );
  const double ghi = fg->gamma_channel_upper( gap_last );
  for( const RoiAtom &a : atoms )
  {
    const double hw = atomlayer_core_halfwidth( a.energy, fwhm_at, core_num_fwhm );
    if( !std::isfinite(hw) )
      return false;
    if( a.energy < glo )
    {
      if( (a.energy + hw) > glo )
        return false;
    }else if( a.energy > ghi )
    {
      if( (a.energy - hw) < ghi )
        return false;
    }else
    {
      return false;  // atom lies inside the excluded band
    }
  }
  return true;
}


// Build a component spanning channels [first_ch, last_ch] with the given atoms, copying
// pass-through metadata (continuum/range/protected/joined) from `meta`.
static AutomaticRoiComponent atomlayer_component_by_channels( size_t first_ch, size_t last_ch,
    std::vector<RoiAtom> atoms, const AutomaticRoiComponent &meta,
    const std::shared_ptr<const SpecUtils::Measurement> &fg )
{
  const size_t nchan = fg->num_gamma_channels();
  if( nchan )
  {
    if( first_ch >= nchan ) first_ch = nchan - 1;
    if( last_ch >= nchan ) last_ch = nchan - 1;
  }
  if( last_ch < first_ch )
    last_ch = first_ch;
  AutomaticRoiComponent c;
  c.first_channel = first_ch;
  c.last_channel = last_ch;
  c.lower = fg->gamma_channel_lower( first_ch );
  c.upper = fg->gamma_channel_upper( last_ch );
  std::sort( std::begin(atoms), std::end(atoms),
    []( const RoiAtom &a, const RoiAtom &b ){ return a.energy < b.energy; } );
  c.atoms = std::move( atoms );
  c.joined_groups = meta.joined_groups;
  c.protected_geometry = meta.protected_geometry;
  c.continuum_type = meta.continuum_type;
  c.range_limits_type = meta.range_limits_type;
  return c;
}


// Candidate single-channel boundaries in the anchor gap [search_lo, search_hi] that keep every
// `constraint_atoms` core intact, ordered nearest-first to `target_energy`.
static std::vector<size_t> atomlayer_core_safe_boundaries(
    const std::vector<RoiAtom> &constraint_atoms, const double target_energy,
    const double search_lo, const double search_hi,
    const std::shared_ptr<const SpecUtils::Measurement> &fg,
    const std::function<double(double)> &fwhm_at, const double core_num_fwhm )
{
  std::vector<size_t> out;
  if( !(search_hi > search_lo) )
    return out;
  const size_t lo_ch = atomlayer_channel_for_energy( fg, search_lo );
  const size_t hi_ch = atomlayer_channel_for_energy( fg, search_hi );
  if( hi_ch < lo_ch )
    return out;
  const size_t target_ch = atomlayer_channel_for_energy( fg, target_energy );
  for( size_t bc = lo_ch; bc <= hi_ch; ++bc )
  {
    if( atomlayer_gap_core_safe( bc, bc, constraint_atoms, fg, fwhm_at, core_num_fwhm ) )
      out.push_back( bc );
  }
  std::stable_sort( std::begin(out), std::end(out),
    [target_ch]( const size_t a, const size_t b ){
      const size_t da = (a > target_ch) ? (a - target_ch) : (target_ch - a);
      const size_t db = (b > target_ch) ? (b - target_ch) : (target_ch - b);
      return da < db;
    } );
  return out;
}


// Widen a child's OUTER edge to reach `min_width_fwhm` (never crossing the boundary, a barrier, or
// the spectrum edge).  Returns true if the child meets the width (or no minimum is imposed / the
// resolution is unavailable); false if it cannot reach the width outward.
static bool atomlayer_widen_child( AutomaticRoiComponent &child, const bool extend_lower_edge,
    const double min_width_fwhm, const double lowest, const double highest,
    const double left_barrier, const std::shared_ptr<const SpecUtils::Measurement> &fg,
    const std::function<double(double)> &fwhm_at )
{
  if( min_width_fwhm <= 0.0 )
    return true;
  double mid = 0.5 * (child.lower + child.upper);
  double f = fwhm_at( mid );
  if( !std::isfinite(f) || (f <= 0.0) )
    return true;
  const double need = min_width_fwhm * f;
  if( (child.upper - child.lower) >= need )
    return true;
  if( extend_lower_edge )
  {
    double target = child.upper - need;
    double floor_energy = lowest;
    if( std::isfinite(left_barrier) )
      floor_energy = std::max( floor_energy, left_barrier );
    target = std::max( target, floor_energy );
    size_t fc = atomlayer_channel_for_energy( fg, target );
    if( std::isfinite(left_barrier) )
    {
      const size_t bar_ch = atomlayer_channel_for_energy( fg, left_barrier );
      if( fc <= bar_ch )
        fc = bar_ch + 1;
    }
    if( fc > child.last_channel )
      return false;
    child.first_channel = fc;
    child.lower = fg->gamma_channel_lower( fc );
  }else
  {
    double target = std::min( child.lower + need, highest );
    size_t lc = atomlayer_channel_for_energy( fg, target );
    if( lc < child.first_channel )
      return false;
    child.last_channel = lc;
    child.upper = fg->gamma_channel_upper( lc );
  }
  mid = 0.5 * (child.lower + child.upper);
  f = fwhm_at( mid );
  if( !std::isfinite(f) || (f <= 0.0) )
    return true;
  return (child.upper - child.lower) >= (min_width_fwhm * f - 1.0e-9);
}


// Count atoms that ended up in a different original group than they started (spatial reassignment).
static size_t atomlayer_count_reassigned( const std::vector<RoiAtom> &orig_left,
    const std::vector<RoiAtom> &orig_right, const AutomaticRoiComponent &new_left,
    const AutomaticRoiComponent &new_right )
{
  std::set<uint64_t> left_ids, right_ids;
  for( const RoiAtom &a : orig_left ) left_ids.insert( a.id );
  for( const RoiAtom &a : orig_right ) right_ids.insert( a.id );
  size_t moved = 0;
  for( const RoiAtom &a : new_left.atoms )
    if( right_ids.count( a.id ) ) ++moved;
  for( const RoiAtom &a : new_right.atoms )
    if( left_ids.count( a.id ) ) ++moved;
  return moved;
}


// Handle a pair where at least one side is protected user/mixed geometry: pin the boundary to the
// protected edge, trim only the automatic side, book atoms whose cores fall inside the protected
// range to it (bounds untouched), and orphan any atom whose core straddles the pin.
static AutomaticRoiPartitionResult atomlayer_partition_protected(
    const AutomaticRoiComponent &left, const AutomaticRoiComponent &right,
    const std::shared_ptr<const SpecUtils::Measurement> &fg,
    const std::function<double(double)> &fwhm_at, const double core_num_fwhm,
    const std::string &stage )
{
  AutomaticRoiPartitionResult result;
  result.policy.decision = AutomaticRoiDecision::ProtectedGeometry;
  result.policy.diagnostic.decision = AutomaticRoiDecision::ProtectedGeometry;
  result.policy.diagnostic.stage = stage;
  result.policy.diagnostic.left_lower = left.lower;
  result.policy.diagnostic.left_upper = left.upper;
  result.policy.diagnostic.right_lower = right.lower;
  result.policy.diagnostic.right_upper = right.upper;

  // Both protected (should not overlap in practice): keep both unchanged.
  if( left.protected_geometry && right.protected_geometry )
  {
    result.outcome = AutomaticRoiPartitionOutcome::KeptSeparate;
    result.components = { left, right };
    result.policy.diagnostic.reason = "both sides protected; retained unchanged";
    return result;
  }

  const AutomaticRoiComponent &prot = left.protected_geometry ? left : right;
  const AutomaticRoiComponent &autom = left.protected_geometry ? right : left;
  const bool protected_is_left = left.protected_geometry;

  // Distribute the automatic side's atoms: inside the protected range -> booked to protected;
  // clear of the pin on the automatic side -> stays; straddling the pin -> orphaned.
  std::vector<RoiAtom> prot_atoms = prot.atoms;
  std::vector<RoiAtom> auto_atoms;
  const size_t nchan = fg->num_gamma_channels();

  // Automatic side's retained channel span (channel-disjoint from the protected range).
  size_t auto_first = autom.first_channel, auto_last = autom.last_channel;
  if( protected_is_left )
  {
    const size_t pin_ch = atomlayer_channel_for_energy( fg, prot.upper );
    auto_first = std::max( auto_first, pin_ch + 1 );
  }else
  {
    const size_t pin_ch = atomlayer_channel_for_energy( fg, prot.lower );
    if( pin_ch == 0 )
    {
      auto_first = 1;  // no room left of channel 0 -> force an empty automatic span (dissolve)
      auto_last = 0;
    }else
    {
      auto_last = std::min( auto_last, pin_ch - 1 );
    }
  }
  const bool auto_span_valid = (nchan > 0) && (auto_last >= auto_first) && (auto_first < nchan);
  const double auto_lower = auto_span_valid ? fg->gamma_channel_lower( auto_first ) : 0.0;
  const double auto_upper = auto_span_valid ? fg->gamma_channel_upper( auto_last ) : 0.0;

  for( const RoiAtom &a : autom.atoms )
  {
    const double hw = atomlayer_core_halfwidth( a.energy, fwhm_at, core_num_fwhm );
    const double clo = std::isfinite(hw) ? (a.energy - hw) : a.energy;
    const double chi = std::isfinite(hw) ? (a.energy + hw) : a.energy;
    if( (clo >= prot.lower) && (chi <= prot.upper) )
    {
      prot_atoms.push_back( a );  // core fully inside protected range
    }else if( auto_span_valid && (clo >= auto_lower) && (chi <= auto_upper) )
    {
      auto_atoms.push_back( a );  // core fully clear on the automatic side
    }else
    {
      result.orphaned_atoms.push_back( a );  // straddles the pin / spectrum edge
    }
  }

  AutomaticRoiComponent prot_out = atomlayer_component_by_channels(
      prot.first_channel, prot.last_channel, std::move(prot_atoms), prot, fg );
  // Protected bounds must remain bit-identical to the input.
  prot_out.lower = prot.lower;
  prot_out.upper = prot.upper;
  prot_out.first_channel = prot.first_channel;
  prot_out.last_channel = prot.last_channel;

  if( !result.orphaned_atoms.empty() )
    result.infeasible_reason = "automatic atom core straddles protected boundary";

  if( auto_span_valid && !auto_atoms.empty() )
  {
    AutomaticRoiComponent auto_out = atomlayer_component_by_channels(
        auto_first, auto_last, std::move(auto_atoms), autom, fg );
    result.outcome = AutomaticRoiPartitionOutcome::KeptSeparate;
    if( protected_is_left )
      result.components = { prot_out, auto_out };
    else
      result.components = { auto_out, prot_out };
    result.policy.diagnostic.reason = "protected geometry pinned; automatic side trimmed";
  }else
  {
    // The automatic side dissolved (its atoms went to the protected range or were orphaned).
    result.outcome = AutomaticRoiPartitionOutcome::Merged;
    result.components = { prot_out };
    result.policy.diagnostic.reason = "protected geometry absorbed automatic neighbor";
  }
  return result;
}


AutomaticRoiPartitionResult partition_automatic_roi_pair(
    const AutomaticRoiComponent &left,
    const AutomaticRoiComponent &right,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints )
{
  const bool have_cal = foreground && foreground->energy_calibration()
      && foreground->energy_calibration()->valid() && (foreground->num_gamma_channels() > 0);

  // Merge fallback: one component covering both, atoms concatenated - always atom-safe.
  const auto build_merged = [&]( const AutomaticRoiPolicyResult &policy, const std::string &reason,
      const bool infeasible ) -> AutomaticRoiPartitionResult
  {
    std::vector<RoiAtom> atoms = left.atoms;
    atoms.insert( std::end(atoms), std::begin(right.atoms), std::end(right.atoms) );
    AutomaticRoiPartitionResult r;
    r.outcome = AutomaticRoiPartitionOutcome::Merged;
    r.policy = policy;
    r.policy.diagnostic.stage = settings.stage;
    r.policy.diagnostic.reason = reason;
    r.policy.diagnostic.partition_infeasible = infeasible;
    if( have_cal )
    {
      const size_t uf = std::min( atomlayer_channel_for_energy( foreground, left.lower ),
                                  atomlayer_channel_for_energy( foreground, right.lower ) );
      const size_t ul = std::max( atomlayer_channel_for_energy( foreground, left.upper ),
                                  atomlayer_channel_for_energy( foreground, right.upper ) );
      AutomaticRoiComponent merged = atomlayer_component_by_channels( uf, ul, std::move(atoms),
          left, foreground );
      merged.joined_groups = left.joined_groups + right.joined_groups;
      merged.protected_geometry = left.protected_geometry || right.protected_geometry;
      r.components = { merged };
    }else
    {
      AutomaticRoiComponent merged = left;
      merged.lower = std::min( left.lower, right.lower );
      merged.upper = std::max( left.upper, right.upper );
      std::sort( std::begin(atoms), std::end(atoms),
        []( const RoiAtom &a, const RoiAtom &b ){ return a.energy < b.energy; } );
      merged.atoms = std::move( atoms );
      merged.joined_groups = left.joined_groups + right.joined_groups;
      merged.protected_geometry = left.protected_geometry || right.protected_geometry;
      r.components = { merged };
    }
    return r;
  };//build_merged

  // Self-validation wrapper: any result that fails the invariant collapses to the merge fallback
  // (dev builds assert loudly first).
  const auto finalize = [&]( AutomaticRoiPartitionResult r ) -> AutomaticRoiPartitionResult
  {
    const AutomaticRoiTransactionCheck chk = validate_automatic_roi_transaction(
        { left, right }, r.components, r.orphaned_atoms, foreground, fwhm_at_energy,
        constraints.peak_core_num_fwhm );
    if( chk.valid )
      return r;
#if( PERFORM_DEVELOPER_CHECKS )
    assert( 0 && "partition_automatic_roi_pair produced an invalid transaction" );
#endif
    if( left.protected_geometry || right.protected_geometry )
      return r;  // cannot merge across protected geometry; keep the (flagged) result
    return build_merged( r.policy, "partition validation failed; merged as fallback", true );
  };//finalize

  if( !have_cal || left.atoms.empty() || right.atoms.empty() )
  {
    AutomaticRoiPolicyResult policy;
    policy.decision = AutomaticRoiDecision::MergeInseparable;
    if( left.protected_geometry || right.protected_geometry )
      return finalize( atomlayer_partition_protected( left, right, foreground, fwhm_at_energy,
          constraints.peak_core_num_fwhm, settings.stage ) );
    return finalize( build_merged( policy,
        have_cal ? "a side carries no atoms; merged" : "no valid calibration; merged", !have_cal ) );
  }

  if( left.protected_geometry || right.protected_geometry )
    return finalize( atomlayer_partition_protected( left, right, foreground, fwhm_at_energy,
        constraints.peak_core_num_fwhm, settings.stage ) );

  // Ask the (unchanged) oracle whether these groups may be separated.
  AutomaticRoiGroup left_group;
  left_group.lower = left.lower;
  left_group.upper = left.upper;
  left_group.joined_groups = left.joined_groups;
  for( const RoiAtom &a : left.atoms ){ left_group.peak_energies.push_back( a.energy );
                                        left_group.peak_areas.push_back( a.area ); }
  AutomaticRoiGroup right_group;
  right_group.lower = right.lower;
  right_group.upper = right.upper;
  right_group.joined_groups = right.joined_groups;
  for( const RoiAtom &a : right.atoms ){ right_group.peak_energies.push_back( a.energy );
                                         right_group.peak_areas.push_back( a.area ); }

  const AutomaticRoiPolicyResult policy = evaluate_automatic_roi_boundary( left_group, right_group,
      foreground, global_continuum, fwhm_at_energy, unfit_auto_peaks, settings );

  if( (policy.decision == AutomaticRoiDecision::MergeInseparable)
      || (policy.decision == AutomaticRoiDecision::MergeInseparableWide) )
    return finalize( build_merged( policy, policy.diagnostic.reason, false ) );

  // KeepSeparate / UnmodeledFeatureBlocked -> materialize a core-safe partition.
  const size_t union_first = std::min( atomlayer_channel_for_energy( foreground, left.lower ),
                                       atomlayer_channel_for_energy( foreground, right.lower ) );
  const size_t union_last = std::max( atomlayer_channel_for_energy( foreground, left.upper ),
                                      atomlayer_channel_for_energy( foreground, right.upper ) );

  std::vector<RoiAtom> union_atoms = left.atoms;
  union_atoms.insert( std::end(union_atoms), std::begin(right.atoms), std::end(right.atoms) );

  double left_max = -std::numeric_limits<double>::infinity();
  for( const RoiAtom &a : left.atoms ) left_max = std::max( left_max, a.energy );
  double right_min = std::numeric_limits<double>::infinity();
  for( const RoiAtom &a : right.atoms ) right_min = std::min( right_min, a.energy );
  double target_boundary = policy.boundary_energy;
  if( !((target_boundary > left_max) && (target_boundary < right_min)) )
    target_boundary = 0.5 * (left_max + right_min);

  // Build the two children at a chosen gap band [gap_first, gap_last], assigning atoms spatially
  // and widening any under-width child outward.  Returns false if a child cannot be materialized.
  const auto try_build_children = [&]( const size_t gap_first, const size_t gap_last,
      AutomaticRoiPartitionResult &out ) -> bool
  {
    if( (gap_first <= union_first) || (gap_last >= union_last) || (gap_last < gap_first) )
      return false;
    const double gap_lo_edge = foreground->gamma_channel_lower( gap_first );
    const double gap_hi_edge = foreground->gamma_channel_upper( gap_last );
    std::vector<RoiAtom> left_atoms, right_atoms;
    for( const RoiAtom &a : union_atoms )
    {
      if( a.energy < gap_lo_edge )
        left_atoms.push_back( a );
      else if( a.energy > gap_hi_edge )
        right_atoms.push_back( a );
      else
        return false;  // atom inside the excluded band (should be pre-filtered by core-safety)
    }
    if( left_atoms.empty() || right_atoms.empty() )
      return false;
    AutomaticRoiComponent left_child = atomlayer_component_by_channels( union_first, gap_first - 1,
        std::move(left_atoms), left, foreground );
    AutomaticRoiComponent right_child = atomlayer_component_by_channels( gap_last + 1, union_last,
        std::move(right_atoms), right, foreground );
    if( !atomlayer_widen_child( left_child, /*extend_lower_edge=*/true, constraints.min_width_fwhm,
          constraints.lowest_energy, constraints.highest_energy, constraints.left_barrier,
          foreground, fwhm_at_energy ) )
      return false;
    if( !atomlayer_widen_child( right_child, /*extend_lower_edge=*/false, constraints.min_width_fwhm,
          constraints.lowest_energy, constraints.highest_energy,
          -std::numeric_limits<double>::infinity(), foreground, fwhm_at_energy ) )
      return false;
    out.outcome = AutomaticRoiPartitionOutcome::KeptSeparate;
    out.components = { left_child, right_child };
    out.policy = policy;
    out.policy.diagnostic.stage = settings.stage;
    out.policy.diagnostic.atoms_reassigned = atomlayer_count_reassigned( left.atoms, right.atoms,
        left_child, right_child );
    return true;
  };//try_build_children

  if( policy.decision == AutomaticRoiDecision::UnmodeledFeatureBlocked )
  {
    // Try the exclusion band first, but never carve through an admitted atom core.
    if( policy.exclusion_upper > policy.exclusion_lower )
    {
      const size_t band_first = atomlayer_channel_for_energy( foreground, policy.exclusion_lower );
      const size_t band_last = atomlayer_channel_for_energy( foreground, policy.exclusion_upper );
      AutomaticRoiPartitionResult carved;
      if( atomlayer_gap_core_safe( band_first, band_last, union_atoms, foreground, fwhm_at_energy,
            constraints.peak_core_num_fwhm )
          && try_build_children( band_first, band_last, carved ) )
      {
        carved.policy.diagnostic.reason = "unmodeled feature excluded between children";
        return finalize( carved );
      }
    }
    // Fallback: a plain single-channel boundary that clears both admitted cores AND the unmodeled
    // peaks (added as constraint pseudo-atoms so the boundary avoids their cores too).
    std::vector<RoiAtom> constraint_atoms = union_atoms;
    for( const std::shared_ptr<const PeakDef> &pk : unfit_auto_peaks )
    {
      if( !pk || !pk->gausPeak() )
        continue;
      const double m = pk->mean();
      if( (m > left_max) && (m < right_min) )
      {
        RoiAtom pseudo;
        pseudo.energy = m;
        pseudo.kind = RoiAtomKind::FloatingFeature;
        constraint_atoms.push_back( pseudo );
      }
    }
    const std::vector<size_t> candidates = atomlayer_core_safe_boundaries( constraint_atoms,
        target_boundary, left_max, right_min, foreground, fwhm_at_energy,
        constraints.peak_core_num_fwhm );
    for( const size_t bc : candidates )
    {
      AutomaticRoiPartitionResult split;
      if( try_build_children( bc, bc, split ) )
      {
        split.policy.diagnostic.reason = "core-safe boundary avoiding unmodeled feature";
        return finalize( split );
      }
    }
    return finalize( build_merged( policy,
        "unmodeled feature inseparable from admitted cores; merged", true ) );
  }

  // KeepSeparate: pick the core-safe boundary nearest the oracle's boundary energy.
  const std::vector<size_t> candidates = atomlayer_core_safe_boundaries( union_atoms,
      target_boundary, left_max, right_min, foreground, fwhm_at_energy,
      constraints.peak_core_num_fwhm );
  for( const size_t bc : candidates )
  {
    AutomaticRoiPartitionResult split;
    if( try_build_children( bc, bc, split ) )
      return finalize( split );
  }

  // No core-safe boundary exists (or none leaves both children wide enough) -> keep the incumbent
  // geometry by merging rather than dropping a requested line.
  return finalize( build_merged( policy, "no core-safe boundary; retained merged", true ) );
}//partition_automatic_roi_pair


AutomaticRoiComponentPartitionResult partition_overwide_automatic_component(
    const std::vector<AutomaticRoiComponent> &component,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints )
{
  AutomaticRoiComponentPartitionResult out;
  AutomaticRoiDecisionDiagnostic &diag = out.diagnostic;
  diag.stage = settings.stage;
  diag.calibration_num_channels = foreground ? foreground->num_gamma_channels() : 0;
  if( component.empty() || !foreground || !foreground->energy_calibration()
      || !foreground->energy_calibration()->valid() || !fwhm_at_energy )
  {
    out.failure_reason = diag.reason = "invalid whole-component partition input";
    return out;
  }
  if( std::any_of( std::begin(component), std::end(component),
        []( const AutomaticRoiComponent &c ) { return c.protected_geometry; } ) )
  {
    out.components = component;
    out.valid = true;
    diag.decision = AutomaticRoiDecision::ProtectedGeometry;
    diag.reason = "protected geometry bypasses whole-component partition";
    return out;
  }

  size_t union_first = component.front().first_channel;
  size_t union_last = component.front().last_channel;
  size_t joined_groups = 0;
  std::vector<RoiAtom> atoms;
  for( const AutomaticRoiComponent &c : component )
  {
    union_first = std::min( union_first, c.first_channel );
    union_last = std::max( union_last, c.last_channel );
    joined_groups += c.joined_groups;
    atoms.insert( std::end(atoms), std::begin(c.atoms), std::end(c.atoms) );
  }
  std::sort( std::begin(atoms), std::end(atoms),
    []( const RoiAtom &lhs, const RoiAtom &rhs ) { return lhs.energy < rhs.energy; } );
  AutomaticRoiComponent incumbent = atomlayer_component_by_channels(
      union_first, union_last, atoms, component.front(), foreground );
  incumbent.joined_groups = std::max<size_t>( 1, joined_groups );
  incumbent.protected_geometry = false;
  out.components = { incumbent };
  diag.left_lower = incumbent.lower;
  diag.left_upper = incumbent.upper;

  const double union_mid = 0.5*(incumbent.lower + incumbent.upper);
  const double union_fwhm = fwhm_at_energy( union_mid );
  if( !std::isfinite(union_fwhm) || !(union_fwhm > 0.0) || (atoms.size() < 2) )
  {
    out.valid = true;
    diag.decision = AutomaticRoiDecision::MergeInseparableWide;
    diag.reason = "too few finite anchors for whole-component partition";
    return out;
  }
  diag.combined_width_fwhm = (incumbent.upper - incumbent.lower) / union_fwhm;
  const size_t union_num_channels = 1 + union_last - union_first;
  const double union_width_ratio = (settings.max_width_fwhm > 0.0)
      ? (diag.combined_width_fwhm / settings.max_width_fwhm) : 0.0;
  const double union_excess = std::max( 0.0, union_width_ratio - 1.0 );
  diag.width_pressure = settings.continuum_aicc_penalty * union_excess * union_excess
      * std::max<size_t>( 1, incumbent.joined_groups - 1 ) * union_num_channels;
  if( !(diag.width_pressure > 0.0) )
  {
    out.valid = true;
    diag.decision = AutomaticRoiDecision::MergeInseparable;
    diag.reason = "component is below the soft-width onset";
    return out;
  }

  struct Anchor
  {
    std::vector<RoiAtom> atoms;
    double mean = 0.0;
    double sigma = 0.0;
    double core_lower = 0.0;
    double core_upper = 0.0;
  };
  std::vector<Anchor> anchors;
  for( const RoiAtom &atom : atoms )
  {
    const double half_width = atomlayer_core_halfwidth(
        atom.energy, fwhm_at_energy, constraints.peak_core_num_fwhm );
    if( !std::isfinite(half_width) )
      continue;
    const double core_lower = atom.energy - half_width;
    if( anchors.empty() || (core_lower > anchors.back().core_upper) )
    {
      Anchor anchor;
      anchor.atoms.push_back( atom );
      anchor.core_lower = core_lower;
      anchor.core_upper = atom.energy + half_width;
      anchors.push_back( std::move(anchor) );
    }else
    {
      anchors.back().atoms.push_back( atom );
      anchors.back().core_lower = std::min( anchors.back().core_lower, core_lower );
      anchors.back().core_upper = std::max(
          anchors.back().core_upper, atom.energy + half_width );
    }
  }
  for( Anchor &anchor : anchors )
  {
    double sum_weights = 0.0;
    double sum_energy = 0.0;
    for( const RoiAtom &atom : anchor.atoms )
    {
      const double weight = (std::isfinite(atom.area) && (atom.area > 0.0)) ? atom.area : 1.0;
      sum_weights += weight;
      sum_energy += weight*atom.energy;
    }
    anchor.mean = sum_energy / std::max( 1.0, sum_weights );
    double sum_variance = 0.0;
    for( const RoiAtom &atom : anchor.atoms )
    {
      const double weight = (std::isfinite(atom.area) && (atom.area > 0.0)) ? atom.area : 1.0;
      const double sigma = fwhm_at_energy(atom.energy) / PhysicalUnits::fwhm_nsigma;
      sum_variance += weight*(sigma*sigma + (atom.energy - anchor.mean)*(atom.energy - anchor.mean));
    }
    anchor.sigma = std::sqrt( sum_variance / std::max(1.0, sum_weights) );
  }
  if( anchors.size() < 2 )
  {
    out.valid = true;
    diag.decision = AutomaticRoiDecision::MergeInseparableWide;
    diag.reason = "all modeled atoms occupy one FWHM-connected core";
    return out;
  }
  struct SegmentCandidate
  {
    size_t first_anchor = 0;
    size_t last_anchor = 0;
    size_t first_channel = 0;
    size_t last_channel = 0;
    double objective = 0.0;
    size_t num_parameters = 0;
    PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;
  };
  const auto segment_candidates = [&]( const size_t first_anchor,
                                        const size_t last_anchor,
                                        const size_t first,
                                        const size_t last ) {
    std::vector<SegmentCandidate> result;
    if( last <= first )
      return result;
    const double lower = foreground->gamma_channel_lower( first );
    const double upper = foreground->gamma_channel_upper( last );
    const double midpoint = 0.5*(lower + upper);
    const double fwhm = fwhm_at_energy( midpoint );
    if( !std::isfinite(fwhm) || !(fwhm > 0.0) )
      return result;
    const double width_fwhm = (upper - lower) / fwhm;
    if( (constraints.min_width_fwhm > 0.0) && (width_fwhm < constraints.min_width_fwhm) )
      return result;

    std::vector<double> means, sigmas;
    for( size_t anchor_index = first_anchor; anchor_index <= last_anchor; ++anchor_index )
    {
      means.push_back( anchors[anchor_index].mean );
      sigmas.push_back( anchors[anchor_index].sigma );
    }
    const std::vector<PeakDef> no_fixed_peaks;
    std::vector<MeasuredRoiModelFit> fits;
    static_cast<void>( fit_measured_roi_model( foreground, first, last, midpoint,
        means, sigmas, no_fixed_peaks, settings.continuum_aicc_penalty, &fits ) );
    const double width_ratio = (settings.max_width_fwhm > 0.0)
        ? (width_fwhm / settings.max_width_fwhm) : 0.0;
    const double excess = std::max( 0.0, width_ratio - 1.0 );
    const size_t joined = 1 + last_anchor - first_anchor;
    const double width_pressure = settings.continuum_aicc_penalty * excess * excess
        * std::max<size_t>( 1, joined - 1 ) * (1 + last - first);
    for( const MeasuredRoiModelFit &fit : fits )
    {
      SegmentCandidate candidate;
      candidate.first_anchor = first_anchor;
      candidate.last_anchor = last_anchor;
      candidate.first_channel = first;
      candidate.last_channel = last;
      candidate.objective = fit.poisson_deviance + width_pressure;
      candidate.num_parameters = fit.num_parameters;
      candidate.continuum_type = fit.continuum_type;
      result.push_back( candidate );
    }
    return result;
  };
  double best_incumbent_score = std::numeric_limits<double>::max();
  double best_partition_score = std::numeric_limits<double>::max();
  double best_forced_partition_score = std::numeric_limits<double>::max();
  SegmentCandidate best_left, best_right;
  SegmentCandidate best_forced_left, best_forced_right;
  const std::vector<SegmentCandidate> incumbent_candidates
    = segment_candidates( 0, anchors.size() - 1, union_first, union_last );
  for( const SegmentCandidate &candidate : incumbent_candidates )
  {
    const double score = data_only_aicc( candidate.objective, union_num_channels,
        candidate.num_parameters, settings.continuum_aicc_penalty );
    best_incumbent_score = std::min( best_incumbent_score, score );
  }
  if( !(best_incumbent_score < std::numeric_limits<double>::max()) )
  {
    out.valid = true;
    diag.one_roi_aicc = best_incumbent_score;
    diag.two_roi_aicc = std::numeric_limits<double>::max();
    diag.decision = AutomaticRoiDecision::MergeInseparableWide;
    diag.partition_infeasible = true;
    diag.reason = "over-wide incumbent has no valid common-channel measured-data fit";
    return out;
  }

  const auto has_clean_gap = [&]( const size_t gap ) {
    std::vector<double> left_energies, left_areas, right_energies, right_areas;
    for( size_t index = 0; index <= gap; ++index )
      for( const RoiAtom &atom : anchors[index].atoms )
      {
        left_energies.push_back( atom.energy );
        left_areas.push_back( std::max( 0.0, atom.area ) );
      }
    for( size_t index = gap + 1; index < anchors.size(); ++index )
      for( const RoiAtom &atom : anchors[index].atoms )
      {
        right_energies.push_back( atom.energy );
        right_areas.push_back( std::max( 0.0, atom.area ) );
      }
    return find_clean_gap_between( left_energies, left_areas, right_energies, right_areas,
        anchors[gap].mean, anchors[gap + 1].mean, foreground, fwhm_at_energy,
        settings.merge_tail_z, settings.merge_clean_gap_fwhm, nullptr, nullptr,
        settings.global_continuum );
  };
  const auto has_residual_valley = [&]( const size_t gap ) {
    // This deliberately is not an ordinary merge decision.  It is only an opt-in late
    // sparse-component admission rail: a resolved modeled structure may have appreciable,
    // correctly modeled tails in its valley, while still having no evidence for an omitted peak
    // or continuum bridge.  Requiring the shared SNIP estimate prevents a local fitted
    // continuum from manufacturing that absence-of-excess evidence.
    if( !(settings.residual_valley_max_excess_z > 0.0)
        || !settings.global_continuum || !settings.global_continuum->valid()
        || (settings.global_continuum->foreground != foreground) )
      return false;
    std::vector<double> left_energies, left_areas, right_energies, right_areas;
    for( size_t index = 0; index <= gap; ++index )
      for( const RoiAtom &atom : anchors[index].atoms )
      {
        left_energies.push_back( atom.energy );
        left_areas.push_back( std::max( 0.0, atom.area ) );
      }
    for( size_t index = gap + 1; index < anchors.size(); ++index )
      for( const RoiAtom &atom : anchors[index].atoms )
      {
        right_energies.push_back( atom.energy );
        right_areas.push_back( std::max( 0.0, atom.area ) );
      }
    const double midpoint = 0.5*(anchors[gap].mean + anchors[gap + 1].mean);
    const double fwhm = fwhm_at_energy( midpoint );
    if( !std::isfinite(fwhm) || !(fwhm > 0.0) )
      return false;
    const double window_width = settings.merge_clean_gap_fwhm * fwhm;
    if( (anchors[gap + 1].mean - anchors[gap].mean) < window_width )
      return false;
    const auto tails = [&]( const double lo, const double hi ) {
      return predicted_gaussian_counts( left_energies, left_areas, fwhm_at_energy, lo, hi )
        + predicted_gaussian_counts( right_energies, right_areas, fwhm_at_energy, lo, hi );
    };
    const double step = 0.25 * fwhm;
    for( double lo = anchors[gap].mean;
         (lo + window_width) <= (anchors[gap + 1].mean + 1.0e-9); lo += step )
    {
      const double hi = lo + window_width;
      const double observed = foreground->gamma_integral(
          static_cast<float>(lo), static_cast<float>(hi) );
      const double continuum = settings.global_continuum->integral( lo, hi );
      const double modeled_tails = tails( lo, hi );
      const double variance = std::max( 1.0, observed + continuum + modeled_tails
          + settings.global_continuum->integral_variance( lo, hi ) );
      const double excess_z = std::max( 0.0, observed - continuum - modeled_tails )
          / std::sqrt( variance );
      if( excess_z <= settings.residual_valley_max_excess_z )
        return true;
    }
    return false;
  };
  std::vector<bool> clean_gaps( anchors.size() - 1, false );
  double strongest_clean_separation = 0.0;
  if( settings.allow_clean_gap_partition_override
      || (settings.residual_valley_max_excess_z > 0.0) )
  {
    for( size_t gap = 0; gap + 1 < anchors.size(); ++gap )
    {
      clean_gaps[gap] = has_clean_gap( gap ) || has_residual_valley( gap );
      if( clean_gaps[gap] )
        strongest_clean_separation = std::max( strongest_clean_separation,
            (anchors[gap + 1].mean - anchors[gap].mean) / union_fwhm );
    }
  }

  // A binary challenger cannot represent a broad sparse component with several independent
  // valleys: its global AICc winner can be an arbitrary central cut, while the physically clear
  // valleys remain untouched.  In explicitly enabled clean-gap mode, form a bounded partition
  // from the strongest non-overlapping clean valleys.  Each child is still fit on the original
  // channel bins and the caller later re-solves the complete transaction, so this only proposes
  // geometry; it never discards or silently reassigns modeled atoms.
  if( (settings.allow_clean_gap_partition_override
          || (settings.residual_valley_max_excess_z > 0.0))
      && (settings.max_partition_children > 2) )
  {
    struct CleanBoundary
    {
      size_t gap = 0;
      size_t channel = 0;
      double separation = 0.0;
      double core_gap_fwhm = 0.0;
    };
    std::vector<CleanBoundary> boundaries;
    for( size_t gap = 0; gap + 1 < anchors.size(); ++gap )
    {
      if( !clean_gaps[gap] )
        continue;
      const double core_gap_fwhm = (anchors[gap + 1].core_lower
          - anchors[gap].core_upper) / union_fwhm;
      if( (settings.minimum_partition_gap_fwhm > 0.0)
          && (core_gap_fwhm < settings.minimum_partition_gap_fwhm) )
        continue;
      const double target = 0.5*(anchors[gap].core_upper + anchors[gap + 1].core_lower);
      size_t best_channel = union_last;
      double best_distance = std::numeric_limits<double>::max();
      for( size_t channel = union_first; channel < union_last; ++channel )
      {
        const double edge = foreground->gamma_channel_upper( channel );
        if( (edge < anchors[gap].core_upper) || (edge > anchors[gap + 1].core_lower) )
          continue;
        const bool crosses_unmodeled_core = std::any_of( std::begin(unfit_auto_peaks),
            std::end(unfit_auto_peaks), [&fwhm_at_energy, edge](
                const std::shared_ptr<const PeakDef> &peak ) {
              if( !peak || !peak->gausPeak() )
                return false;
              const double peak_fwhm = fwhm_at_energy( peak->mean() );
              return std::isfinite(peak_fwhm) && (peak_fwhm > 0.0)
                  && (edge > (peak->mean() - peak_fwhm))
                  && (edge < (peak->mean() + peak_fwhm));
            } );
        if( crosses_unmodeled_core )
          continue;
        const double distance = std::fabs( edge - target );
        if( distance < best_distance )
        {
          best_distance = distance;
          best_channel = channel;
        }
      }
      if( best_channel < union_last )
      {
        CleanBoundary boundary;
        boundary.gap = gap;
        boundary.channel = best_channel;
        boundary.separation = (anchors[gap + 1].mean - anchors[gap].mean) / union_fwhm;
        boundary.core_gap_fwhm = core_gap_fwhm;
        boundaries.push_back( boundary );
      }
    }
    std::sort( std::begin(boundaries), std::end(boundaries),
      []( const CleanBoundary &lhs, const CleanBoundary &rhs ) {
        return lhs.separation > rhs.separation;
      } );
    const size_t max_boundaries = settings.max_partition_children - 1;
    if( boundaries.size() > max_boundaries )
      boundaries.resize( max_boundaries );
    std::sort( std::begin(boundaries), std::end(boundaries),
      []( const CleanBoundary &lhs, const CleanBoundary &rhs ) {
        return lhs.gap < rhs.gap;
      } );
    if( boundaries.size() >= 2 )
    {
      std::vector<SegmentCandidate> selected;
      size_t first_anchor = 0;
      size_t first_channel = union_first;
      bool feasible = true;
      for( size_t index = 0; index <= boundaries.size(); ++index )
      {
        const size_t last_anchor = (index < boundaries.size())
          ? boundaries[index].gap : (anchors.size() - 1);
        const size_t last_channel = (index < boundaries.size())
          ? boundaries[index].channel : union_last;
        std::vector<SegmentCandidate> candidates = segment_candidates(
            first_anchor, last_anchor, first_channel, last_channel );
        if( candidates.empty() )
        {
          feasible = false;
          break;
        }
        const auto best = std::min_element( std::begin(candidates), std::end(candidates),
          []( const SegmentCandidate &lhs, const SegmentCandidate &rhs ) {
            return lhs.objective < rhs.objective;
          } );
        selected.push_back( *best );
        first_anchor = last_anchor + 1;
        first_channel = last_channel + 1;
      }
      if( feasible )
      {
        double objective = 0.0;
        size_t num_parameters = 0;
        for( const SegmentCandidate &segment : selected )
        {
          objective += segment.objective;
          num_parameters += segment.num_parameters;
        }
        const double score = data_only_aicc( objective, union_num_channels, num_parameters,
            settings.continuum_aicc_penalty );
        const bool force = (settings.force_partition_gap_fwhm > 0.0)
          && std::all_of( std::begin(boundaries), std::end(boundaries),
            [&settings]( const CleanBoundary &boundary ) {
              return boundary.core_gap_fwhm >= settings.force_partition_gap_fwhm;
            } );
        if( (score < best_incumbent_score) || force )
        {
          std::vector<AutomaticRoiComponent> children;
          children.reserve( selected.size() );
          for( const SegmentCandidate &segment : selected )
          {
            std::vector<RoiAtom> child_atoms;
            for( size_t anchor_index = segment.first_anchor;
                 anchor_index <= segment.last_anchor; ++anchor_index )
              child_atoms.insert( std::end(child_atoms), std::begin(anchors[anchor_index].atoms),
                                  std::end(anchors[anchor_index].atoms) );
            AutomaticRoiComponent child = atomlayer_component_by_channels(
                segment.first_channel, segment.last_channel, std::move(child_atoms), incumbent,
                foreground );
            child.joined_groups = 1 + segment.last_anchor - segment.first_anchor;
            child.continuum_type = segment.continuum_type;
            children.push_back( std::move(child) );
          }
          const AutomaticRoiTransactionCheck check = validate_automatic_roi_transaction(
              component, children, {}, foreground, fwhm_at_energy, constraints.peak_core_num_fwhm );
          if( check.valid )
          {
            out.components = std::move(children);
            out.valid = true;
            out.changed = true;
            diag.decision = AutomaticRoiDecision::KeepSeparate;
            diag.two_roi_aicc = score;
            diag.reason = force
              ? "configured clean core gaps force a bounded multi-boundary partition"
              : "measured-data AICc favors a bounded clean-valley multi-boundary partition";
            diag.left_lower = out.components.front().lower;
            diag.left_upper = out.components.front().upper;
            diag.right_lower = out.components.back().lower;
            diag.right_upper = out.components.back().upper;
            diag.boundary_channel = boundaries.front().channel;
            diag.boundary_energy = foreground->gamma_channel_upper( diag.boundary_channel );
            return out;
          }
        }
      }
    }
  }

  // Search every channel edge between adjacent, FWHM-distinct source cores.  Candidate children
  // use exactly the same union channels as the incumbent and no edge may cross an unmodeled peak
  // core.  The selected structural alternative is exactly two continua, so diagnostics and the
  // later component transaction can carry it without lossy reconstruction.
  for( size_t gap = 0; gap + 1 < anchors.size(); ++gap )
  {
    const double anchor_separation = (anchors[gap + 1].mean - anchors[gap].mean) / union_fwhm;
    // In explicit clean-valley mode, split the strongest supported valley first.  AICc still
    // decides whether that physical proposal is worthwhile, but cannot substitute an arbitrary
    // central cut merely because it balances the two child continua.
    if( settings.allow_clean_gap_partition_override
        && (!clean_gaps[gap]
            || (anchor_separation + 1.0e-6 < strongest_clean_separation)) )
      continue;
    const double core_gap_fwhm = (anchors[gap + 1].core_lower
        - anchors[gap].core_upper) / union_fwhm;
    bool clean_gap_override = false;
    if( settings.allow_clean_gap_partition_override
        && (settings.minimum_partition_gap_fwhm > 0.0)
        && (core_gap_fwhm < settings.minimum_partition_gap_fwhm) )
    {
      clean_gap_override = clean_gaps[gap];
    }
    if( (settings.minimum_partition_gap_fwhm > 0.0)
        && (core_gap_fwhm < settings.minimum_partition_gap_fwhm)
        && !clean_gap_override )
      continue;
    for( size_t split_channel = union_first; split_channel < union_last; ++split_channel )
    {
      const double edge = foreground->gamma_channel_upper( split_channel );
      if( (edge < anchors[gap].core_upper) || (edge > anchors[gap + 1].core_lower) )
        continue;
      const bool crosses_unmodeled_core = std::any_of( std::begin(unfit_auto_peaks),
          std::end(unfit_auto_peaks), [&fwhm_at_energy, edge](
              const std::shared_ptr<const PeakDef> &peak ) {
            if( !peak || !peak->gausPeak() )
              return false;
            const double peak_fwhm = fwhm_at_energy( peak->mean() );
            return std::isfinite(peak_fwhm) && (peak_fwhm > 0.0)
                && (edge > (peak->mean() - peak_fwhm))
                && (edge < (peak->mean() + peak_fwhm));
          } );
      if( crosses_unmodeled_core )
        continue;
      const std::vector<SegmentCandidate> left_candidates
        = segment_candidates( 0, gap, union_first, split_channel );
      const std::vector<SegmentCandidate> right_candidates
        = segment_candidates( gap + 1, anchors.size() - 1,
            split_channel + 1, union_last );
      for( const SegmentCandidate &left : left_candidates )
      {
        for( const SegmentCandidate &right : right_candidates )
        {
          const double score = data_only_aicc( left.objective + right.objective,
              union_num_channels, left.num_parameters + right.num_parameters,
              settings.continuum_aicc_penalty );
          if( score < best_partition_score )
          {
            best_partition_score = score;
            best_left = left;
            best_right = right;
          }
          // The force-gap rail is intentionally distinct from AICc preference.  Retain the
          // best *eligible* clean-gap challenger, rather than testing the arbitrary global
          // AICc winner: a slightly lower score at a narrow safe edge must not hide a broader,
          // explicitly configured separation between the same two source cores.
          if( (settings.force_partition_gap_fwhm > 0.0)
              && (core_gap_fwhm >= settings.force_partition_gap_fwhm)
              && (score < best_forced_partition_score) )
          {
            best_forced_partition_score = score;
            best_forced_left = left;
            best_forced_right = right;
          }
        }
      }
    }
  }
  const bool aicc_prefers_partition = best_partition_score < best_incumbent_score;
  const bool forced_by_clean_gap = !aicc_prefers_partition
    && (best_forced_partition_score < std::numeric_limits<double>::max());
  if( forced_by_clean_gap )
  {
    best_partition_score = best_forced_partition_score;
    best_left = best_forced_left;
    best_right = best_forced_right;
  }
  diag.one_roi_aicc = best_incumbent_score;
  diag.two_roi_aicc = best_partition_score;
  const bool incumbent_valid = best_incumbent_score < std::numeric_limits<double>::max();
  const bool partition_valid = best_partition_score < std::numeric_limits<double>::max();
  if( !incumbent_valid || !partition_valid
      || (!aicc_prefers_partition && !forced_by_clean_gap) )
  {
    out.valid = true;
    diag.decision = AutomaticRoiDecision::MergeInseparableWide;
    diag.reason = (incumbent_valid && partition_valid)
        ? "measured-data whole-component AICc retains the over-wide union"
        : "no feasible measured-data whole-component partition";
    diag.partition_infeasible = !incumbent_valid || !partition_valid;
    return out;
  }

  std::vector<AutomaticRoiComponent> children;
  const SegmentCandidate selected_segments[] = { best_left, best_right };
  children.reserve( 2 );
  for( const SegmentCandidate &segment : selected_segments )
  {
    std::vector<RoiAtom> child_atoms;
    for( size_t anchor_index = segment.first_anchor;
         anchor_index <= segment.last_anchor; ++anchor_index )
      child_atoms.insert( std::end(child_atoms), std::begin(anchors[anchor_index].atoms),
                          std::end(anchors[anchor_index].atoms) );
    AutomaticRoiComponent child = atomlayer_component_by_channels(
        segment.first_channel, segment.last_channel, std::move(child_atoms), incumbent, foreground );
    child.joined_groups = 1 + segment.last_anchor - segment.first_anchor;
    child.continuum_type = segment.continuum_type;
    children.push_back( std::move(child) );
  }
  const AutomaticRoiTransactionCheck check = validate_automatic_roi_transaction(
      component, children, {}, foreground, fwhm_at_energy, constraints.peak_core_num_fwhm );
  if( !check.valid )
  {
    out.valid = true;
    out.changed = false;
    out.failure_reason = check.failure_reason;
    diag.decision = AutomaticRoiDecision::MergeInseparableWide;
    diag.partition_infeasible = true;
    diag.reason = "whole-component challenger failed atom transaction; retained union";
    return out;
  }

  out.components = std::move( children );
  out.valid = true;
  out.changed = true;
  diag.decision = AutomaticRoiDecision::KeepSeparate;
  diag.reason = forced_by_clean_gap
    ? "configured clean core gap overrides measured-data whole-component AICc"
    : "measured-data whole-component AICc and soft width favor a core-safe partition";
  if( out.components.size() > 1 )
  {
    diag.left_lower = out.components.front().lower;
    diag.left_upper = out.components.front().upper;
    diag.right_lower = out.components.back().lower;
    diag.right_upper = out.components.back().upper;
    diag.boundary_channel = out.components.front().last_channel;
    diag.boundary_energy = foreground->gamma_channel_upper( diag.boundary_channel );
  }
  return out;
}//partition_overwide_automatic_component


AutomaticRoiTransactionCheck validate_automatic_roi_transaction(
    const std::vector<AutomaticRoiComponent> &before,
    const std::vector<AutomaticRoiComponent> &after,
    const std::vector<RoiAtom> &reported_orphans,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    const double peak_core_num_fwhm )
{
  AutomaticRoiTransactionCheck chk;

  // 1. Atom-ID multiset preserved exactly, and exactly-once ownership after.
  std::multiset<uint64_t> before_ids, after_ids;
  for( const AutomaticRoiComponent &c : before )
    for( const RoiAtom &a : c.atoms ) before_ids.insert( a.id );
  for( const AutomaticRoiComponent &c : after )
    for( const RoiAtom &a : c.atoms ) after_ids.insert( a.id );
  for( const RoiAtom &a : reported_orphans ) after_ids.insert( a.id );
  if( before_ids != after_ids )
  {
    chk.failure_reason = "atom-ID multiset changed across the operation";
    return chk;
  }
  std::set<uint64_t> seen;
  for( const AutomaticRoiComponent &c : after )
    for( const RoiAtom &a : c.atoms )
      if( !seen.insert( a.id ).second )
      {
        chk.failure_reason = "an atom is owned by more than one component";
        return chk;
      }

  // 2. Components sorted and strictly channel-disjoint.
  for( size_t i = 1; i < after.size(); ++i )
  {
    if( after[i].first_channel <= after[i-1].last_channel )
    {
      std::ostringstream message;
      message << "components " << (i - 1) << " and " << i
              << " are not channel-disjoint: [" << after[i-1].lower << ", "
              << after[i-1].upper << "] channels " << after[i-1].first_channel << "-"
              << after[i-1].last_channel << " vs [" << after[i].lower << ", "
              << after[i].upper << "] channels " << after[i].first_channel << "-"
              << after[i].last_channel;
      chk.failure_reason = message.str();
      return chk;
    }
  }

  // 3. Atom energy and (clamped) core containment within the owning component.
  const bool have_cal = foreground && (foreground->num_gamma_channels() > 0);
  const double spec_lo = have_cal ? foreground->gamma_channel_lower( 0 ) : 0.0;
  const double spec_hi = have_cal
      ? foreground->gamma_channel_upper( foreground->num_gamma_channels() - 1 ) : 0.0;
  for( const AutomaticRoiComponent &c : after )
  {
    for( const RoiAtom &a : c.atoms )
    {
      if( (a.energy < c.lower - 1.0e-6) || (a.energy > c.upper + 1.0e-6) )
      {
        chk.failure_reason = "an atom energy lies outside its component";
        return chk;
      }
      const double hw = atomlayer_core_halfwidth( a.energy, fwhm_at_energy, peak_core_num_fwhm );
      if( std::isfinite(hw) && have_cal )
      {
        const double clo = std::max( a.energy - hw, spec_lo );
        const double chi = std::min( a.energy + hw, spec_hi );
        if( (clo < c.lower - 1.0e-6) || (chi > c.upper + 1.0e-6) )
        {
          std::ostringstream message;
          message << "atom " << a.id << " core [" << clo << ", " << chi
                  << "] lies outside component [" << c.lower << ", " << c.upper << "]";
          chk.failure_reason = message.str();
          return chk;
        }
      }
    }
  }

  // 4. Every protected component appears in `after` with bit-identical bounds/metadata.
  const auto matches = []( const AutomaticRoiComponent &a, const AutomaticRoiComponent &b ){
    return (std::fabs(a.lower - b.lower) < 1.0e-6) && (std::fabs(a.upper - b.upper) < 1.0e-6)
        && (a.continuum_type == b.continuum_type) && (a.range_limits_type == b.range_limits_type);
  };
  for( const AutomaticRoiComponent &b : before )
  {
    if( !b.protected_geometry )
      continue;
    const bool found = std::any_of( std::begin(after), std::end(after),
        [&]( const AutomaticRoiComponent &a ){ return a.protected_geometry && matches(a, b); } );
    if( !found )
    {
      chk.failure_reason = "a protected component's bounds/metadata changed";
      return chk;
    }
  }

  chk.valid = true;
  return chk;
}//validate_automatic_roi_transaction


static AutomaticRoiReconcileResult reconcile_automatic_components_impl(
    std::vector<AutomaticRoiComponent> components,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints,
    std::vector<AutomaticRoiDecisionDiagnostic> *diagnostics,
    const bool reconcile_transitive_collisions )
{
  AutomaticRoiReconcileResult out;
  const std::vector<AutomaticRoiComponent> before = components;

  const auto component_less = []( const AutomaticRoiComponent &a,
                                  const AutomaticRoiComponent &b ){
    if( a.lower != b.lower )
      return a.lower < b.lower;
    if( a.upper != b.upper )
      return a.upper < b.upper;
    return a.first_channel < b.first_channel;
  };

  std::sort( std::begin(components), std::end(components), component_less );

  // Preserve the established adjacent decisions first.  A globally sorted repair pass below is
  // entered only if pair materialization leaves a later channel collision.
  std::vector<AutomaticRoiComponent> resolved;
  std::vector<RoiAtom> orphans;
  for( const AutomaticRoiComponent &component : components )
  {
    if( resolved.empty()
        || (component.first_channel > resolved.back().last_channel) )
    {
      resolved.push_back( component );
      continue;
    }

    const bool back_empty = resolved.back().atoms.empty();
    const bool component_empty = component.atoms.empty();
    const bool drop_component = component_empty && !component.protected_geometry;
    const bool drop_back = !drop_component && back_empty && !resolved.back().protected_geometry;
    if( drop_component || drop_back )
    {
      if( diagnostics )
      {
        AutomaticRoiDecisionDiagnostic insufficient;
        insufficient.decision = AutomaticRoiDecision::KeepSeparate;
        insufficient.stage = settings.stage;
        insufficient.left_lower = resolved.back().lower;
        insufficient.left_upper = resolved.back().upper;
        insufficient.right_lower = component.lower;
        insufficient.right_upper = component.upper;
        insufficient.calibration_num_channels = foreground ? foreground->num_gamma_channels() : 0;
        insufficient.reason = "late ROI lacks modeled peak evidence; rejected automatic addition";
        diagnostics->push_back( insufficient );
      }
      if( drop_back )
        resolved.back() = component;
      continue;
    }

    AutomaticRoiPartitionConstraints pair_constraints = constraints;
    pair_constraints.left_barrier = (resolved.size() >= 2)
        ? resolved[resolved.size() - 2].upper : constraints.left_barrier;
    const AutomaticRoiPartitionResult partition = partition_automatic_roi_pair(
        resolved.back(), component, foreground, global_continuum, fwhm_at_energy,
        unfit_auto_peaks, settings, pair_constraints );
    if( diagnostics )
      diagnostics->push_back( partition.policy.diagnostic );
    orphans.insert( std::end(orphans), std::begin(partition.orphaned_atoms),
                    std::end(partition.orphaned_atoms) );
    if( (partition.outcome == AutomaticRoiPartitionOutcome::KeptSeparate)
        && (partition.components.size() == 2) )
    {
      resolved.back() = partition.components[0];
      resolved.push_back( partition.components[1] );
    }else if( !partition.components.empty() )
    {
      resolved.back() = partition.components[0];
    }
  }

  // A pair fold is useful for protected/late collision handling, but an over-wide merged result
  // gets one final whole-component transaction so the last pair visited cannot dictate the split.
  std::vector<AutomaticRoiComponent> jointly_partitioned;
  for( const AutomaticRoiComponent &component : resolved )
  {
    const double midpoint = 0.5*(component.lower + component.upper);
    const double fwhm = fwhm_at_energy ? fwhm_at_energy(midpoint)
                                       : std::numeric_limits<double>::quiet_NaN();
    const bool overwide = settings.allow_overwide_overlap_partition
        && !component.protected_geometry && std::isfinite(fwhm) && (fwhm > 0.0)
        && (settings.max_width_fwhm > 0.0)
        && (((component.upper - component.lower) / fwhm) > settings.max_width_fwhm);
    if( !overwide )
    {
      jointly_partitioned.push_back( component );
      continue;
    }
    AutomaticRoiPartitionConstraints local_constraints = constraints;
    local_constraints.left_barrier = jointly_partitioned.empty()
      ? constraints.left_barrier : jointly_partitioned.back().upper;
    const AutomaticRoiComponentPartitionResult partition
      = partition_overwide_automatic_component( { component }, foreground, fwhm_at_energy,
          unfit_auto_peaks, settings, local_constraints );
    if( diagnostics )
      diagnostics->push_back( partition.diagnostic );
    if( partition.valid && partition.changed )
      jointly_partitioned.insert( std::end(jointly_partitioned),
          std::begin(partition.components), std::end(partition.components) );
    else
      jointly_partitioned.push_back( component );
  }
  resolved = std::move( jointly_partitioned );

  // A split can expand its right child for the minimum-width requirement past the lower bound of
  // a component that the adjacent fold has not visited yet.  Re-sort and repeat globally until
  // every channel belongs to at most one component.
  bool converged = !reconcile_transitive_collisions;
  const size_t max_reconcile_steps = 16 * (resolved.size() + 1) * (resolved.size() + 1);
  for( size_t step = 0; reconcile_transitive_collisions && (step < max_reconcile_steps); ++step )
  {
    std::sort( std::begin(resolved), std::end(resolved), component_less );
    size_t right_index = resolved.size();
    for( size_t index = 1; index < resolved.size(); ++index )
    {
      if( resolved[index].first_channel <= resolved[index - 1].last_channel )
      {
        right_index = index;
        break;
      }
    }
    if( right_index >= resolved.size() )
    {
      converged = true;
      break;
    }

    const size_t left_index = right_index - 1;
    const AutomaticRoiComponent left = resolved[left_index];
    const AutomaticRoiComponent right = resolved[right_index];
    std::vector<AutomaticRoiComponent> replacement;
    const bool left_empty = left.atoms.empty();
    const bool right_empty = right.atoms.empty();
    const bool drop_right = right_empty && !right.protected_geometry;
    const bool drop_left = !drop_right && left_empty && !left.protected_geometry;
    if( drop_right || drop_left )
    {
      if( diagnostics )
      {
        AutomaticRoiDecisionDiagnostic insufficient;
        insufficient.decision = AutomaticRoiDecision::KeepSeparate;
        insufficient.stage = settings.stage;
        insufficient.left_lower = left.lower;
        insufficient.left_upper = left.upper;
        insufficient.right_lower = right.lower;
        insufficient.right_upper = right.upper;
        insufficient.calibration_num_channels = foreground ? foreground->num_gamma_channels() : 0;
        insufficient.reason = "late ROI lacks modeled peak evidence; rejected automatic addition";
        diagnostics->push_back( insufficient );
      }
      replacement.push_back( drop_left ? right : left );
    }else
    {
      AutomaticRoiPartitionConstraints pair_constraints = constraints;
      pair_constraints.left_barrier = (left_index >= 1)
          ? resolved[left_index - 1].upper : constraints.left_barrier;
      const AutomaticRoiPartitionResult partition = partition_automatic_roi_pair(
          left, right, foreground, global_continuum, fwhm_at_energy, unfit_auto_peaks,
          settings, pair_constraints );
      if( diagnostics )
        diagnostics->push_back( partition.policy.diagnostic );
      orphans.insert( std::end(orphans), std::begin(partition.orphaned_atoms),
                      std::end(partition.orphaned_atoms) );
      replacement = partition.components;

      const bool unchanged_overlap = (replacement.size() == 2)
          && (replacement[0].first_channel == left.first_channel)
          && (replacement[0].last_channel == left.last_channel)
          && (replacement[1].first_channel == right.first_channel)
          && (replacement[1].last_channel == right.last_channel);
      if( replacement.empty() || unchanged_overlap )
        break;
    }

    resolved.erase( std::begin(resolved) + static_cast<std::ptrdiff_t>(left_index),
                    std::begin(resolved) + static_cast<std::ptrdiff_t>(right_index + 1) );
    resolved.insert( std::begin(resolved) + static_cast<std::ptrdiff_t>(left_index),
                     std::begin(replacement), std::end(replacement) );
  }
  std::sort( std::begin(resolved), std::end(resolved), component_less );

  const AutomaticRoiTransactionCheck chk = validate_automatic_roi_transaction( before, resolved,
      orphans, foreground, fwhm_at_energy, constraints.peak_core_num_fwhm );
  out.components = std::move( resolved );
  out.orphaned_atoms = std::move( orphans );
  out.valid = converged && chk.valid;
  out.failure_reason = converged ? chk.failure_reason
                                  : "automatic ROI reconciliation did not converge";
#if( PERFORM_DEVELOPER_CHECKS )
  assert( out.valid );
#endif
  return out;
}//reconcile_automatic_components_impl


AutomaticRoiReconcileResult reconcile_automatic_components(
    std::vector<AutomaticRoiComponent> components,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints,
    std::vector<AutomaticRoiDecisionDiagnostic> *diagnostics )
{
  return reconcile_automatic_components_impl( std::move(components), foreground,
      global_continuum, fwhm_at_energy, unfit_auto_peaks, settings, constraints,
      diagnostics, true );
}


static AutomaticRoiReconcileResult reconcile_automatic_components_one_pass(
    std::vector<AutomaticRoiComponent> components,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints,
    std::vector<AutomaticRoiDecisionDiagnostic> *diagnostics )
{
  // The downstream resolver and private merge pass retain their established transaction; the
  // whole-list retry belongs only to initial automatic clustering.
  return reconcile_automatic_components_impl( std::move(components), foreground,
      global_continuum, fwhm_at_energy, unfit_auto_peaks, settings, constraints,
      diagnostics, false );
}


void assign_atoms_to_disjoint_rois(
    const std::vector<RoiAtom> &universe,
    const std::vector<RelActCalcAuto::RoiRange> &rois,
    std::vector<std::vector<RoiAtom>> &per_roi_atoms,
    std::vector<RoiAtom> &unowned_atoms )
{
  per_roi_atoms.assign( rois.size(), std::vector<RoiAtom>() );
  unowned_atoms.clear();

#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < rois.size(); ++i )
    assert( rois[i].lower_energy >= rois[i-1].upper_energy );
#endif

  for( const RoiAtom &a : universe )
  {
    long best = -1;
    double best_dist = std::numeric_limits<double>::infinity();
    for( size_t j = 0; j < rois.size(); ++j )
    {
      if( (a.energy >= rois[j].lower_energy) && (a.energy <= rois[j].upper_energy) )
      {
        const double mid = 0.5 * (rois[j].lower_energy + rois[j].upper_energy);
        const double dist = std::fabs( a.energy - mid );
        if( (best < 0) || (dist < best_dist) )
        {
          best = static_cast<long>( j );
          best_dist = dist;
        }
      }
    }
    if( best >= 0 )
      per_roi_atoms[static_cast<size_t>(best)].push_back( a );
    else
      unowned_atoms.push_back( a );
  }
}//assign_atoms_to_disjoint_rois


PeakContinuum::OffsetType select_continuum_order_by_sidebands(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const double roi_lower,
  const double roi_upper,
  const double core_lo,
  const double core_hi,
  const double aicc_penalty )
{
  if( !foreground || !foreground->channel_energies() || (foreground->num_gamma_channels() < 8)
     || !(roi_upper > roi_lower) )
    return PeakContinuum::OffsetType::Linear;

  // Gather sideband channels: inside the ROI but outside the peak core.
  std::vector<double> xs, ys;  // channel-center energy (relative to roi_lower), counts
  const size_t first_ch = foreground->find_gamma_channel( static_cast<float>(roi_lower) );
  const size_t last_ch = foreground->find_gamma_channel( static_cast<float>(roi_upper) );

  for( size_t ch = first_ch; (ch <= last_ch) && (ch < foreground->num_gamma_channels()); ++ch )
  {
    const double ch_lo = foreground->gamma_channel_lower( ch );
    const double ch_hi = foreground->gamma_channel_upper( ch );
    if( (ch_hi <= roi_lower) || (ch_lo >= roi_upper) )
      continue;
    if( (ch_hi > core_lo) && (ch_lo < core_hi) )
      continue;  // overlaps the peak core - not continuum

    xs.push_back( 0.5*(ch_lo + ch_hi) - roi_lower );
    ys.push_back( foreground->gamma_channel_content( ch ) );
  }

  const double num_data = static_cast<double>( xs.size() );
  if( xs.size() < 8 )
    return PeakContinuum::OffsetType::Linear;  // too few sideband channels to select on

  // Poisson-weighted least-squares chi2 of a polynomial (in x, counts-per-channel) of the given
  // parameter count, via normal equations solved by Gaussian elimination (max 3x3).
  const auto poly_fit_chi2 = [&xs, &ys]( const size_t num_par ) -> double
  {
    double ata[3][3] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
    double atb[3] = { 0.0, 0.0, 0.0 };

    for( size_t i = 0; i < xs.size(); ++i )
    {
      const double w = 1.0 / std::max( 1.0, ys[i] );  // Poisson variance, floored at 1 count
      double basis[3] = { 1.0, xs[i], xs[i]*xs[i] };
      for( size_t r = 0; r < num_par; ++r )
      {
        atb[r] += w * basis[r] * ys[i];
        for( size_t c = 0; c < num_par; ++c )
          ata[r][c] += w * basis[r] * basis[c];
      }
    }

    // Gaussian elimination with partial pivoting
    double coef[3] = { 0.0, 0.0, 0.0 };
    {
      double a[3][4];
      for( size_t r = 0; r < num_par; ++r )
      {
        for( size_t c = 0; c < num_par; ++c )
          a[r][c] = ata[r][c];
        a[r][num_par] = atb[r];
      }
      for( size_t col = 0; col < num_par; ++col )
      {
        size_t piv = col;
        for( size_t r = col + 1; r < num_par; ++r )
          if( std::fabs(a[r][col]) > std::fabs(a[piv][col]) )
            piv = r;
        if( std::fabs(a[piv][col]) < 1.0e-30 )
          return std::numeric_limits<double>::max();  // degenerate
        if( piv != col )
          for( size_t c = 0; c <= num_par; ++c )
            std::swap( a[piv][c], a[col][c] );
        for( size_t r = col + 1; r < num_par; ++r )
        {
          const double f = a[r][col] / a[col][col];
          for( size_t c = col; c <= num_par; ++c )
            a[r][c] -= f * a[col][c];
        }
      }
      for( size_t col = num_par; col-- > 0; )
      {
        double s = a[col][num_par];
        for( size_t c = col + 1; c < num_par; ++c )
          s -= a[col][c] * coef[c];
        coef[col] = s / a[col][col];
      }
    }

    double chi2 = 0.0;
    for( size_t i = 0; i < xs.size(); ++i )
    {
      const double pred = coef[0] + coef[1]*xs[i] + coef[2]*xs[i]*xs[i];
      const double resid = ys[i] - pred;
      chi2 += (resid * resid) / std::max( 1.0, ys[i] );
    }
    return chi2;
  };//poly_fit_chi2 lambda

  const auto aicc = [&num_data, aicc_penalty]( const double chi2, const double num_par ) -> double {
    if( num_data <= (num_par + 1.0) )
      return std::numeric_limits<double>::max();
    return chi2 + aicc_penalty*num_par
           + (aicc_penalty * num_par * (num_par + 1.0)) / (num_data - num_par - 1.0);
  };

  const double aicc_linear = aicc( poly_fit_chi2( 2 ), 2.0 );
  const double aicc_quad = aicc( poly_fit_chi2( 3 ), 3.0 );

  return (aicc_quad < aicc_linear) ? PeakContinuum::OffsetType::Quadratic
                                   : PeakContinuum::OffsetType::Linear;
}//select_continuum_order_by_sidebands

}//namespace detail


/** Returns auto-search peaks that do NOT correspond to any user peak.

 Matching is by energy proximity: an auto-search peak is considered to match a user peak if their
 means are within 0.15 * auto_peak->fwhm().
 */
std::vector<std::shared_ptr<const PeakDef>> compute_unfit_auto_peaks(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::vector<std::shared_ptr<const PeakDef>> &user_peaks )
{
  std::vector<std::shared_ptr<const PeakDef>> unfit_peaks;

  for( const std::shared_ptr<const PeakDef> &auto_peak : auto_search_peaks )
  {
    if( !auto_peak || !auto_peak->gausPeak() )
      continue;

    const double auto_mean = auto_peak->mean();
    const double tolerance = 0.15 * auto_peak->fwhm();

    bool matches_user_peak = false;
    for( const std::shared_ptr<const PeakDef> &user_peak : user_peaks )
    {
      if( !user_peak || !user_peak->gausPeak() )
        continue;

      if( std::fabs( user_peak->mean() - auto_mean ) < tolerance )
      {
        matches_user_peak = true;
        break;
      }
    }//for( const auto &user_peak : user_peaks )

    if( !matches_user_peak )
      unfit_peaks.push_back( auto_peak );
  }//for( const auto &auto_peak : auto_search_peaks )

  if( should_debug_print() )
  {
    std::cerr << "compute_unfit_auto_peaks: " << auto_search_peaks.size() << " auto peaks, "
         << user_peaks.size() << " user peaks -> " << unfit_peaks.size() << " unfit peaks" << std::endl;
  }

  return unfit_peaks;
}//compute_unfit_auto_peaks(...)


/** Returns the energy of the channel with fewest counts between two energies.

 Searches the foreground spectrum for the channel with the minimum count in the range
 [lower_energy, upper_energy].  Returns the center energy of that channel.
 If the range spans fewer than 2 channels, returns the midpoint energy.
 */
double find_min_counts_energy(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const double lower_energy,
  const double upper_energy )
{
  assert( foreground && foreground->energy_calibration() && foreground->energy_calibration()->valid() );
  assert( lower_energy < upper_energy );

  if( !foreground || !foreground->energy_calibration() || !foreground->energy_calibration()->valid() )
    return 0.5 * (lower_energy + upper_energy);

  const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal
    = foreground->energy_calibration();

  const size_t num_channels = foreground->num_gamma_channels();
  const size_t lower_ch = static_cast<size_t>( energy_cal->channel_for_energy( lower_energy ) );
  const size_t upper_ch = std::min( static_cast<size_t>( energy_cal->channel_for_energy( upper_energy ) ),
                                    num_channels - 1 );

  if( (upper_ch <= lower_ch) || ((upper_ch - lower_ch) < 2) )
    return 0.5 * (lower_energy + upper_energy);

  float min_counts = std::numeric_limits<float>::max();
  size_t min_channel = lower_ch;

  for( size_t ch = lower_ch; ch <= upper_ch; ++ch )
  {
    const float counts = foreground->gamma_channel_content( ch );
    if( counts < min_counts )
    {
      min_counts = counts;
      min_channel = ch;
    }
  }//for( size_t ch = lower_ch; ch <= upper_ch; ++ch )

  // Return center energy of the minimum-counts channel
  return energy_cal->energy_for_channel( static_cast<double>(min_channel) + 0.5 );
}//find_min_counts_energy(...)



/** Shrinks ROIs to avoid interference from auto-search peaks that do not correspond to source/NORM
 gammas.

 For each peak in unfit_auto_peaks:
 - If the peak matches a significant source/NORM gamma (within cluster_num_sigma sigma), skip it.
 - If the peak's ROI overlaps with one of our analysis ROIs, shrink the analysis ROI to the
   channel with fewest counts between the significant gamma energy and the interfering peak mean.
 - Minimum ROI extent from the nearest gamma is enforced: min_fwhm_roi_lower FWHM on the lower
   side, min_fwhm_roi_upper FWHM on the upper side.
 - For multi-gamma ROIs where the nearest gamma (edge gamma) is not the largest gamma in the ROI,
   more aggressive shrinking is allowed: down to edge_gamma ± 0.2*FWHM, while ensuring the
   largest gamma retains its full extent.  The ROI is shrunk until the interfering peak's
   Gaussian integral over the ROI is less than 20% of the edge gamma's expected area.

  TODO: This function is by no means optimized, or super well behaving - but its something that kinda covers some obvious cases
 
 \param rois_and_gammas  ROIs paired with their clustered gamma info (modified in-place)
 \param unfit_auto_peaks  Auto-search peaks not matching user peaks
 \param foreground  Foreground spectrum (may have adjusted energy calibration)
 \param orig_cal  Energy calibration of the original foreground (that auto_search_peaks are in)
 \param current_cal  Current energy calibration (that rois_and_gammas are in)
 */
void shrink_rois_for_interfering_peaks(
  std::vector<std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo>> &rois_and_gammas,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
  const std::vector<float> &fwhm_coefficients,
  const double fwhm_lower_energy,
  const double fwhm_upper_energy,
  const double cluster_num_sigma,
  const double min_fwhm_roi_lower,
  const double min_fwhm_roi_upper,
  const std::shared_ptr<const SpecUtils::EnergyCalibration> &orig_cal,
  const std::shared_ptr<const SpecUtils::EnergyCalibration> &current_cal )
{
  if( unfit_auto_peaks.empty() || rois_and_gammas.empty() )
    return;

  const bool have_fwhm_range = (fwhm_lower_energy > 0.0)
    && (fwhm_upper_energy > 0.0)
    && (fwhm_lower_energy < fwhm_upper_energy);

  // Determine if we need to translate peak energies from original to current calibration
  const bool need_cal_translation = (orig_cal && current_cal && (orig_cal != current_cal));

  // Collect all significant gamma energies across all ROIs, for checking if a peak matches a source gamma
  std::vector<double> all_gamma_energies;
  for( const std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &rg : rois_and_gammas )
  {
    for( const double gamma_e : rg.second.gamma_energies )
      all_gamma_energies.push_back( gamma_e );
  }

  for( const std::shared_ptr<const PeakDef> &peak : unfit_auto_peaks )
  {
    if( !peak || !peak->gausPeak() )
      continue;
    
    // Translate peak energies to current calibration if needed
    double peak_mean = peak->mean();
    double peak_lower = peak->lowerX();
    double peak_upper = peak->upperX();
    //const double peak_fwhm = peak->fwhm(); //TODO: we might be able to this instead of computing the FWHM for each gamma - but there would be edge-cases we should maybe consider.

    if( need_cal_translation )
    {
      const double ch_mean = orig_cal->channel_for_energy( peak_mean );
      peak_mean = current_cal->energy_for_channel( ch_mean );

      const double ch_lower = orig_cal->channel_for_energy( peak_lower );
      peak_lower = current_cal->energy_for_channel( ch_lower );

      const double ch_upper = orig_cal->channel_for_energy( peak_upper );
      peak_upper = current_cal->energy_for_channel( ch_upper );
    }//if( need_cal_translation )

    // Check if this peak matches any significant source/NORM gamma across all ROIs
    bool matches_source_gamma = false;
    for( const double gamma_energy : all_gamma_energies )
    {
      const double fwhm_eval_energy = have_fwhm_range
        ? std::clamp( gamma_energy, fwhm_lower_energy, fwhm_upper_energy )
        : gamma_energy;

      const float gamma_fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(fwhm_eval_energy), fwhm_form, fwhm_coefficients );

      if( !std::isfinite(gamma_fwhm) || (gamma_fwhm <= 0.0f) )
        continue;

      const double gamma_sigma = gamma_fwhm / PhysicalUnits::fwhm_nsigma;

      if( std::fabs( peak_mean - gamma_energy ) < (cluster_num_sigma * gamma_sigma) )
      {
        matches_source_gamma = true;
        break;
      }
    }//for( const double gamma_energy : all_gamma_energies )

    if( matches_source_gamma )
      continue;

    // This is an interfering peak.  Check each ROI for overlap.
    for( std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &rg : rois_and_gammas )
    {
      RelActCalcAuto::RoiRange &roi = rg.first;
      const ClusteredGammaInfo &gamma_info = rg.second;

      // Skip already-invalidated ROIs
      if( (roi.lower_energy < 0.0) || (roi.upper_energy < 0.0)
        || (roi.lower_energy >= roi.upper_energy) )
        continue;

      // Check for overlap between the interfering peak's ROI and our analysis ROI
      if( (peak_upper <= roi.lower_energy) || (peak_lower >= roi.upper_energy) )
        continue;

      // Find the nearest gamma to the interfering peak (edge gamma), and the largest gamma in the ROI
      double nearest_gamma_energy = 0.5 * (roi.lower_energy + roi.upper_energy); // fallback
      double nearest_gamma_amplitude = 0.0;
      double nearest_dist = std::numeric_limits<double>::max();
      size_t largest_gamma_index = 0;

      assert( gamma_info.gamma_energies.size() == gamma_info.gamma_amplitudes.size() );

      for( size_t gi = 0; gi < gamma_info.gamma_energies.size(); ++gi )
      {
        const double gamma_e = gamma_info.gamma_energies[gi];
        const double dist = std::fabs( gamma_e - peak_mean );
        if( dist < nearest_dist )
        {
          nearest_dist = dist;
          nearest_gamma_energy = gamma_e;
          nearest_gamma_amplitude = (gi < gamma_info.gamma_amplitudes.size())
            ? gamma_info.gamma_amplitudes[gi] : 0.0;
        }

        if( (gi < gamma_info.gamma_amplitudes.size())
          && (gamma_info.gamma_amplitudes[gi] > gamma_info.gamma_amplitudes[largest_gamma_index]) )
        {
          largest_gamma_index = gi;
        }
      }//for( size_t gi = 0; gi < gamma_info.gamma_energies.size(); ++gi )

      // Sometimes we'll run into the case where there are multiple gammas in the ROI, and the gamma near this "unfit"
      //  peak is actually a pretty small gamma, so in this case, we'll treat it a little different, and allow a more
      //  agressive shrinking of the ROI - this is so we dont totally mess up the primary line (because the continuum
      //  may fit really high in amplitude if the ROI extends well into a large "unfit" peak that we are otherwise not
      //  taking into account).
      const bool is_multi_gamma = (gamma_info.gamma_energies.size() > 1);
      const bool edge_gamma_is_largest = !is_multi_gamma
        || (std::fabs( nearest_gamma_energy - gamma_info.gamma_energies[largest_gamma_index] ) < 0.01);

      // Compute FWHM at the nearest gamma (edge gamma) for minimum ROI width enforcement
      const double fwhm_eval = have_fwhm_range
        ? std::clamp( nearest_gamma_energy, fwhm_lower_energy, fwhm_upper_energy )
        : nearest_gamma_energy;
      const float fwhm_at_gamma = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(fwhm_eval), fwhm_form, fwhm_coefficients );

      if( !std::isfinite( fwhm_at_gamma ) || (fwhm_at_gamma <= 0.0f) )
        continue;

      // For multi-gamma ROIs with non-dominant edge gamma, compute FWHM at the largest gamma
      double largest_gamma_energy = nearest_gamma_energy;
      float fwhm_at_largest = fwhm_at_gamma;
      if( is_multi_gamma && !edge_gamma_is_largest )
      {
        largest_gamma_energy = gamma_info.gamma_energies[largest_gamma_index];
        const double lg_fwhm_eval = have_fwhm_range
          ? std::clamp( largest_gamma_energy, fwhm_lower_energy, fwhm_upper_energy )
          : largest_gamma_energy;
        fwhm_at_largest = DetectorPeakResponse::peakResolutionFWHM(
          static_cast<float>(lg_fwhm_eval), fwhm_form, fwhm_coefficients );
        if( !std::isfinite( fwhm_at_largest ) || (fwhm_at_largest <= 0.0f) )
          fwhm_at_largest = fwhm_at_gamma;
      }//if( multi-gamma with non-dominant edge gamma )

      const double old_lower = roi.lower_energy;
      const double old_upper = roi.upper_energy;

      if( peak_mean > nearest_gamma_energy )
      {
        // Interfering peak is on the upper side of the nearest gamma
        const double min_count_energy = find_min_counts_energy( foreground, nearest_gamma_energy, peak_mean );

        // Enforce minimum ROI upper extent from the edge gamma
        double effective_min_upper = nearest_gamma_energy + min_fwhm_roi_upper * fwhm_at_gamma;

        // For multi-gamma ROIs where the edge gamma is not the largest, allow more aggressive
        // shrinking down to edge_gamma + 0.2*FWHM, but the largest gamma keeps its full extent.
        // Shrink until the interfering peak's integral in the ROI is < 20% of edge gamma's area.
        // (the 0.2*FWHM was fairly arbitrarily chosen, and not optimized)
        if( is_multi_gamma && !edge_gamma_is_largest
          && (nearest_gamma_amplitude > 0.0) && (peak->peakArea() > 0.0) )
        {
          const double hard_min = nearest_gamma_energy + 0.2 * fwhm_at_gamma;
          const double largest_min = largest_gamma_energy + min_fwhm_roi_upper * fwhm_at_largest;
          const double aggressive_min = std::max( hard_min, largest_min );

          // Only try aggressive shrinking if it allows a tighter boundary than the normal one
          if( aggressive_min < effective_min_upper )
          {
            // Check contribution at the normal minimum
            const double candidate = std::max( min_count_energy, effective_min_upper );
            const double contribution = peak->gauss_integral( roi.lower_energy, candidate );
            const double max_allowed = 0.20 * nearest_gamma_amplitude;

            if( contribution >= max_allowed )
            {
              // Use coverage_limits to find where the interfering peak's left-tail CDF equals
              // max_allowed/peakArea(), i.e., the lower quantile at fraction p/2 = max_allowed/peakArea()
              const double p = 2.0 * max_allowed / peak->peakArea();
              try
              {
                const double * const skew_pars = peak->coefficients() + static_cast<int>(PeakDef::SkewPar0);
                const double quantile_lower = PeakDists::coverage_limits( p, peak->skewType(),
                                                peak->mean(), peak->sigma(), skew_pars ).first;
                effective_min_upper = std::max( quantile_lower, aggressive_min );
              }catch( std::exception & )
              {
                // Fallback to candidate if coverage_limits fails
                effective_min_upper = std::max( candidate, aggressive_min );
              }

              if( should_debug_print() )
              {
                std::cerr << "shrink_rois_for_interfering_peaks: Aggressive shrink for multi-gamma ROI"
                     << " (edge gamma at " << nearest_gamma_energy << " keV is not the largest)"
                     << ", new effective_min_upper=" << effective_min_upper << " keV" << std::endl;
              }
            }//if( contribution >= max_allowed )
          }//if( aggressive_min < effective_min_upper )
        }//if( multi-gamma with non-dominant edge gamma )

        const double new_upper = std::max( min_count_energy, effective_min_upper );

        // Only shrink, never expand
        if( new_upper < roi.upper_energy )
          roi.upper_energy = new_upper;
      }
      else
      {
        // Interfering peak is on the lower side of the nearest gamma
        const double min_count_energy = find_min_counts_energy( foreground, peak_mean, nearest_gamma_energy );

        // Enforce minimum ROI lower extent from the edge gamma
        double effective_max_lower = nearest_gamma_energy - min_fwhm_roi_lower * fwhm_at_gamma;

        // For multi-gamma ROIs where the edge gamma is not the largest, allow more aggressive
        // shrinking down to edge_gamma - 0.2*FWHM, but the largest gamma keeps its full extent.
        if( is_multi_gamma && !edge_gamma_is_largest
          && (nearest_gamma_amplitude > 0.0) && (peak->peakArea() > 0.0) )
        {
          const double hard_max = nearest_gamma_energy - 0.2 * fwhm_at_gamma;
          const double largest_max = largest_gamma_energy - min_fwhm_roi_lower * fwhm_at_largest;
          const double aggressive_max = std::min( hard_max, largest_max );

          // Only try aggressive shrinking if it allows a tighter boundary than the normal one
          if( aggressive_max > effective_max_lower )
          {
            // Check contribution at the normal minimum
            const double candidate = std::min( min_count_energy, effective_max_lower );
            const double contribution = peak->gauss_integral( candidate, roi.upper_energy );
            const double max_allowed = 0.20 * nearest_gamma_amplitude;

            if( contribution >= max_allowed )
            {
              // Use coverage_limits to find where the interfering peak's right-tail CDF equals
              // max_allowed/peakArea(), i.e., the upper quantile at fraction p/2 = max_allowed/peakArea()
              const double p = 2.0 * max_allowed / peak->peakArea();
              try
              {
                const double * const skew_pars = peak->coefficients() + static_cast<int>(PeakDef::SkewPar0);
                const double quantile_upper = PeakDists::coverage_limits( p, peak->skewType(),
                                                peak->mean(), peak->sigma(), skew_pars ).second;
                effective_max_lower = std::min( quantile_upper, aggressive_max );
              }
              catch( std::exception & )
              {
                // Fallback to candidate if coverage_limits fails
                effective_max_lower = std::min( candidate, aggressive_max );
              }

              if( should_debug_print() )
              {
                std::cerr << "shrink_rois_for_interfering_peaks: Aggressive shrink for multi-gamma ROI"
                     << " (edge gamma at " << nearest_gamma_energy << " keV is not the largest)"
                     << ", new effective_max_lower=" << effective_max_lower << " keV" << std::endl;
              }
            }//if( contribution >= max_allowed )
          }//if( aggressive_max > effective_max_lower )
        }//if( multi-gamma with non-dominant edge gamma )

        const double new_lower = std::min( min_count_energy, effective_max_lower );

        // Only shrink, never expand
        if( new_lower > roi.lower_energy )
          roi.lower_energy = new_lower;
      }//if( peak on upper side ) / else

      if( should_debug_print()
        && ((roi.lower_energy != old_lower) || (roi.upper_energy != old_upper)) )
      {
        std::cerr << "shrink_rois_for_interfering_peaks: Shrunk ROI from ["
             << old_lower << ", " << old_upper << "] to ["
             << roi.lower_energy << ", " << roi.upper_energy << "] keV"
             << " due to interfering peak at " << peak_mean << " keV"
             << " (nearest gamma at " << nearest_gamma_energy << " keV)" << std::endl;
      }

      // Mark for removal if ROI collapsed
      if( roi.lower_energy >= roi.upper_energy )
      {
        if( should_debug_print() )
        {
          std::cerr << "shrink_rois_for_interfering_peaks: ROI collapsed and removed"
               << " (was [" << old_lower << ", " << old_upper << "] keV)" << std::endl;
        }
        roi.lower_energy = -1.0;
        roi.upper_energy = -1.0;
        break;
      }
    }//for( auto &rg : rois_and_gammas )
  }//for( const auto &peak : unfit_auto_peaks )

  // Remove any ROIs that were marked as invalid
  rois_and_gammas.erase(
    std::remove_if( rois_and_gammas.begin(), rois_and_gammas.end(),
      []( const std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &rg ) {
        return (rg.first.lower_energy < 0.0) || (rg.first.upper_energy < 0.0)
          || (rg.first.lower_energy >= rg.first.upper_energy);
      } ),
    rois_and_gammas.end()
  );

#if( PERFORM_DEVELOPER_CHECKS )
  for( const std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &rg : rois_and_gammas )
  {
    assert( rg.first.lower_energy < rg.first.upper_energy );
  }
#endif
}//shrink_rois_for_interfering_peaks(...)


// Returns the NORM background nuclides (U238, Ra226, U235, Th232, K40) as NucInputInfo entries,
// excluding any that already appear in `sources`.  Ages are set to prompt/secular equilibrium
// half-lives appropriate for each nuclide (see `getBackgroundRefLines()` in ReferenceLineInfo.cpp).
// If `color_css` is non-empty it is assigned to each entry's peak_color_css so the Rel. Eff.
// chart can render data points for NORM sources.
std::vector<RelActCalcAuto::NucInputInfo> get_norm_sources(
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::string &color_css = {} )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    return {};

  const SandiaDecay::Nuclide * const u238  = db->nuclide( "U238" );
  const SandiaDecay::Nuclide * const ra226 = db->nuclide( "Ra226" );
  const SandiaDecay::Nuclide * const u235  = db->nuclide( "U235" );
  const SandiaDecay::Nuclide * const th232 = db->nuclide( "Th232" );
  const SandiaDecay::Nuclide * const k40   = db->nuclide( "K40" );

  assert( u238 && ra226 && u235 && th232 && k40 );
  if( !u238 || !ra226 || !u235 || !th232 || !k40 )
    return {};

  // {nuclide, activity_scale, age}
  // Activity scales are chosen so representative peaks have specific amplitudes before normalisation;
  // denominators are integrals at the reference energy times the branching ratio.
  // Ages are prompt/secular equilibrium half-lives so all daughters are in equilibrium.
  const std::vector<std::tuple<const SandiaDecay::Nuclide *, double, double>> nuc_activity{
    { u238,   0.0004653/410.2892,  5.0*u238->promptEquilibriumHalfLife()  },
    { ra226,  0.02515/17990.5430,  5.0*ra226->promptEquilibriumHalfLife() },
    { u235,   0.001482/14603.0156, 5.0*u235->promptEquilibriumHalfLife()  },
    { th232,  0.02038/27897.2617,  5.0*th232->secularEquilibriumHalfLife() },
    { k40,    0.1066/6523.8994,    0.0 }
  };//nuc_activity

  std::vector<RelActCalcAuto::NucInputInfo> norm_sources;
  norm_sources.reserve( nuc_activity.size() );

  for( const std::tuple<const SandiaDecay::Nuclide *, double, double> &info : nuc_activity )
  {
    const SandiaDecay::Nuclide * const norm_nuc = std::get<0>( info );

    // Skip if this NORM nuclide is already one of the input sources
    bool already_in_sources = false;
    for( const RelActCalcAuto::NucInputInfo &src : sources )
    {
      if( RelActCalcAuto::nuclide( src.source ) == norm_nuc )
      {
        already_in_sources = true;
        break;
      }
    }
    if( already_in_sources )
      continue;

    RelActCalcAuto::NucInputInfo norm_src;
    norm_src.source        = norm_nuc;
    norm_src.age           = std::get<2>( info );
    norm_src.fit_age       = false;
    norm_src.peak_color_css = color_css;
    norm_sources.push_back( norm_src );
  }

  return norm_sources;
}//get_norm_sources(...)


// Nuclides that should be treated as NORM-like for the peak-fit-improve GA
// background-false-positive penalty, but are not the canonical five NORM
// nuclides handled by get_norm_sources().  These either share many of their
// peaks with NORM lines (e.g. U232/U233 overlap with U238/Th232 chain), or
// produce only the annihilation line (F18) which the broader 511 keV
// handling already controls.  This list is policy, not physics — extend it
// as the 1-generation diagnostic identifies more sources that systematically
// fit peaks on backgrounds.
static const std::array<const char *, 3> sk_norm_like_extras = {
  "U232", "U233", "F18"
};

bool is_norm_like_for_ga( const RelActCalcAuto::SrcVariant &src )
{
  // Elements (xrays) and reactions have no decay chain to test; treat as
  // not NORM-like.  If a future use case needs to exclude particular
  // elements/reactions, add them here.
  const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( src );
  if( !nuc )
    return false;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    return false;

  // Hand-curated extras (string compare against the nuclide symbol).
  for( const char * const sym : sk_norm_like_extras )
  {
    if( nuc->symbol == sym )
      return true;
  }

  // Decay-chain test: forebearers() is recursive and self-seeds with `this`,
  // so this also covers U238/U235/Th232 themselves.  K40 is included as a
  // separate ancestor check since it has no parents in SandiaDecay.
  const SandiaDecay::Nuclide * const u238  = db->nuclide( "U238" );
  const SandiaDecay::Nuclide * const u235  = db->nuclide( "U235" );
  const SandiaDecay::Nuclide * const th232 = db->nuclide( "Th232" );
  const SandiaDecay::Nuclide * const k40   = db->nuclide( "K40" );

  if( nuc == k40 )
    return true;

  const std::vector<const SandiaDecay::Nuclide *> ancestors = nuc->forebearers();
  for( const SandiaDecay::Nuclide * const ancestor : ancestors )
  {
    if( ancestor == u238 || ancestor == u235 || ancestor == th232 )
      return true;
  }

  return false;
}//is_norm_like_for_ga(...)


// Strong gamma + xray lines from NORM nuclides and their decay chains that commonly appear in
// ambient backgrounds, each attributed to the get_norm_sources() parent nuclide that emits it
// (that parent is added aged to equilibrium, so a Pb-214/Bi-214 daughter line is attributed to
// Ra-226, Th-234/Pa-234m to U-238, Ac-228/Tl-208/Bi-212 to Th-232 - matching the parentNuclide()
// of the peaks RelActCalcAuto produces).  parent_symbol is nullptr for the Pb/Bi/Tl K-xrays, which
// are element fluorescence with no single attributable NORM parent.
// Two consumers:
//   - is_near_strong_norm_gamma(): the GA's background-fit penalty, to suppress false-positive
//     penalties when a source's gamma coincides with a real NORM peak.  Uses ONLY the energies; its
//     behavior must stay identical (see NuclideConfig_GA.cpp).
//   - detail::find_strong_unmodeled_interferers(): R6 auto co-fit, which uses parent_symbol to add
//     the interfering NORM nuclide to the fit.
// Hand-curated from K-40 + U-238 chain (Pb-214/Bi-214) + Th-232 chain (Ac-228/Tl-208/Bi-212)
// + Pa-234m + Ra-226 + Th-234 + Pb/Bi/Tl K-xrays.
// Entries are GROUPED BY nuclide/decay-chain for readability, NOT sorted by energy.  Both consumers
// do a linear scan, so ordering is not relied upon; do not switch to binary_search without first
// sorting this array by energy.
struct StrongNormGammaLine
{
  double energy;              // keV
  const char *parent_symbol;  // get_norm_sources() parent symbol; nullptr => K-xray (no parent)
};

static const std::array<StrongNormGammaLine, 51> sk_strong_norm_gamma_lines = {{
  // Pb/Bi/Tl K x-rays (NORM-chain element fluorescence) - no single attributable parent
  {  70.83, nullptr }, {  72.81, nullptr }, {  72.87, nullptr }, {  74.81, nullptr },
  {  74.97, nullptr }, {  77.10, nullptr }, {  84.94, nullptr },
  // Th-234 (U-238 chain)
  {  63.29, "U238" }, {  92.38, "U238" }, {  92.80, "U238" },
  // Ac-228 (Th-232 chain), Pb-214/Bi-214 mix at low energy
  { 129.07, "Th232" },
  // U-235 (top lines - also in is_norm_like_for_ga decay test, listed for
  // mis-attribution overlap on non-U235 sources)
  { 143.76, "U235" }, { 163.36, "U235" }, { 185.72, "U235" }, { 205.31, "U235" },
  // Ra-226 itself
  { 186.21, "Ra226" },
  // Ac-228 (Th-232 chain)
  { 209.25, "Th232" }, { 270.24, "Th232" },
  // Tl-208 (Th-232 chain)
  { 277.36, "Th232" },
  // Pb-214 (Ra-226 chain)
  { 241.99, "Ra226" }, { 258.87, "Ra226" }, { 295.22, "Ra226" }, { 351.93, "Ra226" },
  // Ac-228
  { 328.00, "Th232" }, { 338.32, "Th232" }, { 463.00, "Th232" }, { 562.50, "Th232" },
  // Tl-208
  { 583.19, "Th232" },
  // Bi-214 (Ra-226 chain)
  { 609.31, "Ra226" }, { 768.36, "Ra226" }, { 806.17, "Ra226" }, { 934.06, "Ra226" },
  { 1120.29, "Ra226" }, { 1238.11, "Ra226" }, { 1377.67, "Ra226" }, { 1407.98, "Ra226" },
  { 1509.21, "Ra226" }, { 1661.27, "Ra226" }, { 1764.49, "Ra226" }, { 2204.21, "Ra226" },
  // Bi-212 (Th-232 chain)
  { 727.33, "Th232" }, { 785.37, "Th232" }, { 1620.50, "Th232" },
  // Ac-228
  { 755.32, "Th232" }, { 794.95, "Th232" }, { 911.20, "Th232" }, { 968.97, "Th232" },
  // Tl-208 high-energy
  { 860.56, "Th232" }, { 2614.51, "Th232" },
  // Pa-234m (U-238 chain)
  { 1001.03, "U238" },
  // K-40
  { 1460.82, "K40" }
}};

// Ambient (non-NORM) lines that commonly sit in a foreground and can steal counts from a requested
// source's nearby weak line, just like the NORM lines above.  Used ONLY by R6 interferer detection
// (NOT by is_near_strong_norm_gamma).  Cs-137 661 keV is single-line; Co-60 is the 1173/1332 pair.
static const std::array<StrongNormGammaLine, 3> sk_ambient_interferer_lines = {{
  {  661.657, "Cs137" }, { 1173.228, "Co60" }, { 1332.492, "Co60" }
}};

bool is_near_strong_norm_gamma( const double energy_kev, const double tolerance_kev )
{
  const double tol = std::max( 0.1, tolerance_kev );
  for( const StrongNormGammaLine &line : sk_strong_norm_gamma_lines )
  {
    if( std::fabs( line.energy - energy_kev ) < tol )
      return true;
  }
  return false;
}//is_near_strong_norm_gamma(...)


// R6 interferer detection.  Reopened here (after the strong-NORM table) so it can read
// sk_strong_norm_gamma_lines directly; declared in the header's detail namespace for unit testing.
namespace detail
{

bool is_marginal_keep_reject( const double expected_counts, const double significance,
                              const double keep_z )
{
  return (expected_counts > sm_keep_gate_min_est_counts)
      && !(significance > keep_z)
      && (significance >= (sm_rescue_z_fraction * keep_z));
}//is_marginal_keep_reject(...)

std::vector<InterfererCandidate> find_strong_unmodeled_interferers(
  const std::vector<RequestedSourceGammas> &source_gammas,
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::function<double(double)> &fwhm_at_energy,
  const bool fit_norm_peaks,
  const double min_valid_energy,
  const double max_valid_energy,
  const std::shared_ptr<const SpecUtils::Measurement> &background,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs,
  std::vector<std::string> *warnings,
  const GlobalContinuumEstimate *global_continuum,
  std::vector<double> *guard_energies )
{
  std::vector<InterfererCandidate> candidates;

  if( guard_energies )
    guard_energies->clear();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    return candidates;

  // Nuclides already accounted for by the model: the requested source nuclides, plus the five NORM
  // parents when NORM peaks are being fit (they are already on the NORM rel-eff curve).
  std::set<const SandiaDecay::Nuclide *> modeled_nucs;
  for( const RequestedSourceGammas &sg : source_gammas )
  {
    const SandiaDecay::Nuclide * const n = RelActCalcAuto::nuclide( sg.source );
    if( n )
      modeled_nucs.insert( n );
  }
  if( fit_norm_peaks )
  {
    for( const char * const sym : { "U238", "Ra226", "U235", "Th232", "K40" } )
    {
      const SandiaDecay::Nuclide * const n = db->nuclide( sym );
      if( n )
        modeled_nucs.insert( n );
    }
  }

  const auto fmt_kev = []( const double e ) -> std::string {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(1) << e;
    return ss.str();
  };

  // Max area/uncert of a foreground auto-search peak within the confirm-window of `energy`
  // (0 if none).  A missing amplitude uncert falls back to Poisson sqrt(area).
  const auto confirm_z = [&]( const double energy ) -> double {
    const double fwhm = std::max( 0.1, fwhm_at_energy( energy ) );
    const double tol = sm_interferer_confirm_num_fwhm * fwhm;
    double best_z = 0.0;
    for( const std::shared_ptr<const PeakDef> &p : auto_search_peaks )
    {
      if( !p || (std::fabs( p->mean() - energy ) > tol) )
        continue;
      const double area = p->amplitude();
      const double uncert = p->amplitudeUncert();
      const double z = (uncert > 0.0) ? (area / uncert) : ((area > 0.0) ? std::sqrt( area ) : 0.0);
      best_z = std::max( best_z, z );
    }
    return best_z;
  };

  // Number of significant (>= 10% of the strongest) in-range lines a requested source has; used to
  // decide whether the source is effectively single-line for the doublet guard.
  const auto num_significant_lines = [&]( const RequestedSourceGammas &s ) -> size_t {
    // Contract: energies and yields are parallel.  If a caller violates it, conservatively treat the
    // source as multi-line so the doublet guard never spuriously fires on a mis-sized input.
    if( s.yields.size() != s.energies.size() )
      return s.energies.size();
    double max_yield = 0.0;
    for( size_t i = 0; (i < s.energies.size()) && (i < s.yields.size()); ++i )
    {
      if( (s.energies[i] >= min_valid_energy) && (s.energies[i] <= max_valid_energy) )
        max_yield = std::max( max_yield, s.yields[i] );
    }
    if( max_yield <= 0.0 )
      return s.energies.size();  // no yields supplied: cannot judge, treat each as significant
    size_t n = 0;
    for( size_t i = 0; (i < s.energies.size()) && (i < s.yields.size()); ++i )
    {
      if( (s.energies[i] >= min_valid_energy) && (s.energies[i] <= max_valid_energy)
          && (s.yields[i] >= 0.1*max_yield) )
        ++n;
    }
    return n;
  };

  // Count of strong-NORM table lines in range attributed to `parent`; our proxy for whether the
  // interferer is effectively single-line (K-40) vs a multi-line chain (Ra-226, Th-232, ...).
  const auto interferer_is_single_line = [&]( const SandiaDecay::Nuclide * const parent ) -> bool {
    int n = 0;
    for( const StrongNormGammaLine &l : sk_strong_norm_gamma_lines )
    {
      if( l.parent_symbol && (l.energy >= min_valid_energy) && (l.energy <= max_valid_energy)
          && (db->nuclide( l.parent_symbol ) == parent) )
        ++n;
    }
    for( const StrongNormGammaLine &l : sk_ambient_interferer_lines )
    {
      if( (l.energy >= min_valid_energy) && (l.energy <= max_valid_energy)
          && (db->nuclide( l.parent_symbol ) == parent) )
        ++n;
    }
    return (n <= 1);
  };

  std::set<double> emitted_line_energies;  // de-dup: at most one candidate per interfering line

  // Per-(source gamma, strong line) interferer check, shared by the strong-NORM and ambient (Cs137/
  // Co60) sweeps below.  Emits a nuclide candidate when the line is near a requested-source gamma,
  // its parent is not already modeled, the source's own chain does not explain it, and a foreground
  // auto-search peak data-confirms it.  Skips K-xray entries (null parent).
  const auto check_interferer_line = [&]( const double es, const bool src_single_line,
                                          const double le, const char * const parent_symbol )
  {
    if( !parent_symbol )
      return;
    if( (le < min_valid_energy) || (le > max_valid_energy) || emitted_line_energies.count( le ) )
      return;

    const double fwhm = std::max( 0.1, fwhm_at_energy( le ) );

    // Trigger: a requested-source gamma within the near-window of the line.
    if( std::fabs( le - es ) >= (sm_interferer_near_num_fwhm * fwhm) )
      return;

    const SandiaDecay::Nuclide * const parent = db->nuclide( parent_symbol );
    if( !parent || modeled_nucs.count( parent ) )
      return;

    // Data confirmation: a foreground auto-search peak on the line at >= min detection z.
    const double z = confirm_z( le );
    if( z < sm_interferer_min_detect_z )
      return;

    // An unresolvable doublet is actionable only when the interfering line was actually observed.
    // Check this after foreground confirmation so an absent NORM line never produces a warning.
    if( src_single_line
        && (std::fabs( le - es ) < (sm_interferer_doublet_min_fwhm * fwhm))
        && interferer_is_single_line( parent ) )
    {
      if( guard_energies )
        guard_energies->push_back( le );
      if( warnings )
        warnings->push_back( "Requested source line at " + fmt_kev(es) + " keV overlaps a strong "
          + parent->symbol + " line at " + fmt_kev(le) + " keV within one FWHM; they are an"
          " unresolvable doublet and were not separately co-fit - the fitted area here may be"
          " contaminated." );
      emitted_line_energies.insert( le );
      return;
    }

    // Source-owns-it: for resolvable/multi-line cases, any requested-source gamma within one FWHM
    // means the source's own chain already explains the line.  The confirmed single-line-doublet
    // warning above is intentionally checked first so its ambiguity is not silently discarded.
    for( const RequestedSourceGammas &sg2 : source_gammas )
    {
      for( const double e2 : sg2.energies )
      {
        if( std::fabs( e2 - le ) < (sm_interferer_source_owns_num_fwhm * fwhm) )
          return;
      }
    }

    if( guard_energies )
      guard_energies->push_back( le );
    candidates.push_back( InterfererCandidate{ le, parent, z, /*from_background_search=*/false } );
    emitted_line_energies.insert( le );
  };//check_interferer_line

  // Sweep each requested-source gamma against the active strong-NORM table.  Candidates are
  // foreground-confirmed; the ambient Cs137/Co60 sweep remains intentionally disabled below.
  for( const RequestedSourceGammas &sg : source_gammas )
  {
    const bool src_single_line = (num_significant_lines( sg ) <= 1);
    for( const double es : sg.energies )
    {
      if( (es < min_valid_energy) || (es > max_valid_energy) )
        continue;
      for( const StrongNormGammaLine &line : sk_strong_norm_gamma_lines )
        check_interferer_line( es, src_single_line, line.energy, line.parent_symbol );
      // NOTE: ambient (Cs137/Co60) sweep temporarily disabled - it destabilized the {K40,Eu152}
      // joint fit (0 peaks); under investigation.  See below and scratch/peak_fit_improve/REVIEW_FINDINGS.md.
      // for( const StrongNormGammaLine &line : sk_ambient_interferer_lines )
      //   check_interferer_line( es, src_single_line, line.energy, line.parent_symbol );
    }
  }

  // ---- DEFERRED: a dedicated background auto-search (for disambiguating blended source+ambient
  // features), the ambient Cs137/Co60 sweep, and unattributable (source-less) floating interferer
  // peaks.  These remain future refinements until the multi-source conditioning is understood.
  (void)background;
  (void)drf;
  (void)peak_fit_prefs;
  (void)global_continuum;

  if( guard_energies )
  {
    std::sort( std::begin(*guard_energies), std::end(*guard_energies) );
    guard_energies->erase( std::unique(std::begin(*guard_energies), std::end(*guard_energies)),
                           std::end(*guard_energies) );
  }

  return candidates;
}//find_strong_unmodeled_interferers(...)

}//namespace detail


/** Add a floating peak at 511 keV if appropriate conditions are met.
 
 Physics reasoning for 511 keV floating peak:
 The 511 keV annihilation line can have enhanced intensity from cosmics, high-energy photons create e+e- pairs, etc
 
 For sources like F-18 where 511 keV IS a primary gamma (in top branching ratios),
 we should NOT add a floating peak since the source already accounts for it.
 But for sources with weak 511 lines (not in top ~5 BRs) that have other strong peaks,
 the measured 511 intensity may be dominated by background contributions that aren't
 well-modeled by the source's weak 511 BR. Adding a floating peak in this case allows
 the fit to accommodate excess 511 counts without distorting the relative efficiency curve.
 
 This feature is only enabled for high-resolution detectors (HPGe) since low-resolution
 detectors typically cannot resolve the 511 keV peak clearly enough to benefit from
 this treatment.
 
 \param options The RelActCalcAuto::Options to potentially add a floating peak to
 \param sources The input sources being fit
 \param fit_norm_peaks Whether NORM background peaks are being fit
 \param det_type The coarse resolution type of the detector
 \param min_valid_energy Minimum energy of the valid spectroscopic range (keV)
 \param max_valid_energy Maximum energy of the valid spectroscopic range (keV)
 */
void add_floating_511_peak_if_appropriate(
  RelActCalcAuto::Options &options,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const bool fit_norm_peaks,
  const PeakFitUtils::CoarseResolutionType det_type,
  const double min_valid_energy,
  const double max_valid_energy )
{
  // Only add 511 keV peaks for high-resolution detectors (HPGe)
  if( det_type != PeakFitUtils::CoarseResolutionType::High )
    return;
  
  // Check if there's a ROI covering 511 keV
  const double annihilation_energy = 510.9989;
  const double energy_tolerance = 2.0; // keV tolerance for ROI coverage
  
  bool have_511_roi = false;
  for( const RelActCalcAuto::RoiRange &roi : options.rois )
  {
    if( (annihilation_energy >= (roi.lower_energy - energy_tolerance))
        && (annihilation_energy <= (roi.upper_energy + energy_tolerance)) )
    {
      have_511_roi = true;
      break;
    }
  }
  
  if( !have_511_roi )
  {
    // Remove any existing floating 511 peak (e.g. copied from a previous iteration's options)
    // since there is no longer a ROI covering it - leaving it would cause an error in the solver.
    auto &fps = options.floating_peaks;
    fps.erase( std::remove_if( begin(fps), end(fps),
      [&]( const RelActCalcAuto::FloatingPeak &fp ){
        return std::fabs( fp.energy - annihilation_energy ) < 1.0;
      } ), end(fps) );
    return;
  }
  
  // Check if any source has 511 in top ~5 gamma BRs AND contributes to other ROIs
  bool should_add_floating_511 = false;
  
  // Combine all sources (input sources + NORM sources if fitting them)
  std::vector<RelActCalcAuto::NucInputInfo> all_sources = sources;
  if( fit_norm_peaks )
  {
    const std::vector<RelActCalcAuto::NucInputInfo> norm_sources = get_norm_sources( sources );
    all_sources.insert( all_sources.end(), norm_sources.begin(), norm_sources.end() );
  }
  
  for( const RelActCalcAuto::NucInputInfo &src : all_sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
      continue;
    
    // Get the source's age
    const double src_age = get_source_age( src.source, src.age );
    
    // Get photons for this source
    std::vector<SandiaDecay::EnergyRatePair> photons;
    try
    {
      photons = get_source_photons( src.source, 1.0, src_age );
    }catch( const std::exception & )
    {
      continue; // Skip invalid sources
    }
    
    if( photons.empty() )
      continue;
    
    // Filter photons to only include those within the valid energy range
    std::vector<SandiaDecay::EnergyRatePair> valid_photons;
    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      if( (photon.energy >= min_valid_energy) && (photon.energy <= max_valid_energy) )
        valid_photons.push_back( photon );
    }
    
    if( valid_photons.empty() )
      continue;
    
    // Sort valid photons by intensity (rate) to find top BRs
    std::vector<SandiaDecay::EnergyRatePair> sorted_photons = valid_photons;
    std::sort( sorted_photons.begin(), sorted_photons.end(),
               []( const SandiaDecay::EnergyRatePair &a, const SandiaDecay::EnergyRatePair &b ) {
                 return a.numPerSecond > b.numPerSecond;
               } );
    
    // Check if 511 keV is in the top 5 gamma lines by branching ratio
    // (considering only photons within the valid energy range)
    const size_t num_top_gammas = std::min( size_t(5), sorted_photons.size() );
    bool has_511_in_top_gammas = false;
    for( size_t i = 0; i < num_top_gammas; ++i )
    {
      if( std::fabs( sorted_photons[i].energy - annihilation_energy ) < 1.0 ) // 1 keV tolerance
      {
        has_511_in_top_gammas = true;
        break;
      }
    }
    
    // Check if this source contributes to other ROIs (has gammas in other energy ranges)
    bool contributes_to_other_rois = false;
    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      // Skip photons very close to 511 keV
      if( std::fabs( photon.energy - annihilation_energy ) < 2.0 )
        continue;
      
      // Check if this photon energy falls in any ROI
      for( const RelActCalcAuto::RoiRange &roi : options.rois )
      {
        if( (photon.energy >= roi.lower_energy) && (photon.energy <= roi.upper_energy) )
        {
          contributes_to_other_rois = true;
          break;
        }
      }
      
      if( contributes_to_other_rois )
        break;
    }
    
    // If this source does NOT have 511 in top gammas but DOES contribute to other ROIs,
    // we should add a floating peak at 511
    if( !has_511_in_top_gammas && contributes_to_other_rois )
    {
      should_add_floating_511 = true;
      if( PEAK_FIT_DEBUG_PRINTOUT )
      {
        std::cout << "Source '" << src.name() << "' contributes to other ROIs but does not have "
                  << "511 keV in top 5 gammas - will add floating 511 peak" << std::endl;
      }
      break; // Only need to find one qualifying source
    }
  }//for( loop over all sources )
  
  // Add the floating peak if conditions are met
  if( !should_add_floating_511 )
    return;
  
  // Check if we already have a floating peak at ~511 keV
  bool already_have_511_floating = false;
  for( const RelActCalcAuto::FloatingPeak &fp : options.floating_peaks )
  {
    if( std::fabs( fp.energy - annihilation_energy ) < 1.0 )
    {
      already_have_511_floating = true;
      break;
    }
  }
  
  if( already_have_511_floating )
    return;
  
  RelActCalcAuto::FloatingPeak floating_511;
  floating_511.energy = annihilation_energy; // 510.9989 keV
  
  // Allow FWHM to vary since annihilation peaks can have different widths
  // due to Doppler broadening from positron momentum distribution
  floating_511.release_fwhm = true;
  
  // 510.9989 keV is a true physical energy (electron rest mass), so always Known.
  floating_511.energy_origin = RelActCalcAuto::FloatingPeak::EnergyType::Known;

  options.floating_peaks.push_back( floating_511 );

  if( PEAK_FIT_DEBUG_PRINTOUT )
  {
    std::cout << "Added floating peak at " << annihilation_energy
              << " keV with release_fwhm=true and energy_origin=Known" << std::endl;
  }
}//add_floating_511_peak_if_appropriate(...)


/** Add floating peaks for single and double escape peaks of high-energy gammas if appropriate.
 
 This function checks if auto_search_peaks contains single escape peaks for high-energy gammas
 (like Th232's 2614 keV line). If found and fit_norm_peaks is enabled, it adds floating peaks
 for both single and double escape peaks to allow the fit to properly account for them.
 
 The implementation is kept general to support other high-energy lines in the future (e.g., Ra226).
 
 \param options The RelActCalcAuto::Options to potentially add floating peaks to
 \param auto_search_peaks Peaks found by automated peak search
 \param fit_norm_peaks Whether NORM background peaks are being fit
 \param min_valid_energy Minimum energy of the valid spectroscopic range (keV)
 \param max_valid_energy Maximum energy of the valid spectroscopic range (keV)
 */
void add_escape_peak_floating_peaks_if_appropriate(
  RelActCalcAuto::Options &options,
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const bool fit_norm_peaks,
  const PeakFitUtils::CoarseResolutionType det_type,
  const double min_valid_energy,
  const double max_valid_energy,
  const PeakFitForNuclideConfig &config )
{
  // Only add escape peaks for high-resolution detectors (HPGe)
  // Escape peaks are smeared out and not distinguishable in lower-resolution detectors
  if( !fit_norm_peaks || (det_type != PeakFitUtils::CoarseResolutionType::High) )
    return;
  
  const double electron_rest_mass = 510.9989; // keV
  const double single_escape_offset = electron_rest_mass;
  const double double_escape_offset = 2.0 * electron_rest_mass;
  
  // Define high-energy gamma lines that commonly have escape peaks
  // Escape peaks become significant above ~1.5 MeV
  struct EscapePeakCandidate
  {
    double parent_energy;
    std::string source_name;
  };
  
  const std::vector<EscapePeakCandidate> candidates = {
    { 2614.533, "Th232" },  // Th232 (Tl208) 2614 keV - most prominent
    // Future candidates (commented out for now, but structure supports them):
    // { 2204.21, "Ra226" },  // Ra226 (Bi214) 2204 keV
    // { 1764.49, "Ra226" },  // Ra226 (Bi214) 1764 keV
  };
  
  const double min_parent_energy_for_escape = 1600.0; // keV - escape peaks significant above ~1.5 MeV
  
  // For each high-energy candidate, first find the parent peak in auto_search_peaks
  for( const EscapePeakCandidate &candidate : candidates )
  {
    // Skip if parent energy is outside valid range or too low for escape peaks
    if( (candidate.parent_energy < min_parent_energy_for_escape)
       || (candidate.parent_energy > max_valid_energy) )
      continue;
    
    // Find the parent peak in auto_search_peaks using 0.75 FWHM threshold
    std::shared_ptr<const PeakDef> parent_peak;
    for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
    {
      if( !peak )
        continue;
      
      const double parent_tolerance = 0.75 * peak->fwhm();
      if( std::fabs( peak->mean() - candidate.parent_energy ) < parent_tolerance )
      {
        parent_peak = peak;
        break;
      }
    }
    
    // Calculate theoretical escape peak energies based on the candidate's nominal energy
    // Use theoretical energies for floating peaks, not fitted parent energy
    const double se_energy = candidate.parent_energy - single_escape_offset;
    const double de_energy = candidate.parent_energy - double_escape_offset;

    if( !parent_peak )
    {
      // No parent peak found in auto_search_peaks - remove any orphaned escape floating peaks
      // that may have been copied from a previous solution's options but now lack a ROI.
      // This mirrors how add_floating_511_peak_if_appropriate removes the 511 peak when
      // there is no ROI covering it.
      const double fp_tolerance = 1.0; // keV
      auto &fps = options.floating_peaks;
      fps.erase( std::remove_if( begin(fps), end(fps),
        [&]( const RelActCalcAuto::FloatingPeak &fp ) -> bool {
          if( std::fabs( fp.energy - se_energy ) > fp_tolerance
              && std::fabs( fp.energy - de_energy ) > fp_tolerance )
            return false;
          // Only remove if the floating peak has no covering ROI
          for( const RelActCalcAuto::RoiRange &roi : options.rois )
          {
            if( (fp.energy >= roi.lower_energy) && (fp.energy <= roi.upper_energy) )
              return false;
          }
          return true;
        } ), end(fps) );
      continue;
    }
    
    // Calculate expected positions based on fitted parent for checking auto_search_peaks
    const double se_expected_from_fit = parent_peak->mean() - single_escape_offset;
    const double de_expected_from_fit = parent_peak->mean() - double_escape_offset;
    
    // Skip if escape energies are outside valid range
    if( (se_energy < min_valid_energy) || (de_energy < min_valid_energy) )
      continue;
    
    // Check if single escape peak is in auto_search_peaks using 0.5 FWHM of parent
    // Use the expected position based on fitted parent peak
    const double escape_tolerance = 0.5 * parent_peak->fwhm();
    bool have_se_peak = false;
    for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
    {
      if( peak && (std::fabs( peak->mean() - se_expected_from_fit ) < escape_tolerance) )
      {
        have_se_peak = true;
        break;
      }
    }
    
    if( !have_se_peak )
      continue;  // No S.E. peak found
    
    // Check if there's a ROI covering the parent energy
    bool have_parent_roi = false;
    for( const RelActCalcAuto::RoiRange &roi : options.rois )
    {
      if( (parent_peak->mean() >= roi.lower_energy)
         && (parent_peak->mean() <= roi.upper_energy) )
      {
        have_parent_roi = true;
        break;
      }
    }
    
    if( !have_parent_roi )
      continue;
    
    // Check if double escape peak exists in auto_search_peaks (optional but nice to know)
    // Use the expected position based on fitted parent peak
    bool have_de_peak = false;
    for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
    {
      if( peak && (std::fabs( peak->mean() - de_expected_from_fit ) < escape_tolerance) )
      {
        have_de_peak = true;
        break;
      }
    }
    
    if( PEAK_FIT_DEBUG_PRINTOUT )
    {
      std::cout << "Found parent peak at " << parent_peak->mean() << " keV for " << candidate.source_name
                << " with S.E. at " << se_energy << " keV";
      if( have_de_peak )
        std::cout << " and D.E. at " << de_energy << " keV";
      std::cout << " - adding escape peak floating peaks and ROIs" << std::endl;
    }
    
    // Check if ROIs exist for S.E. and D.E. energies, add if missing.
    // Escape ROIs get the configured core width plus a fixed one-FWHM continuum sideband; they
    // are single isolated peaks so the full adaptive extension machinery is not warranted.
    const double escape_half_width_fwhm = config.auto_roi_core_num_fwhm + sm_post_drop_sideband_fwhm;
    const double roi_width_lower = escape_half_width_fwhm * parent_peak->fwhm();
    const double roi_width_upper = escape_half_width_fwhm * parent_peak->fwhm();
    
    // Check/add S.E. ROI
    bool have_se_roi = false;
    for( const RelActCalcAuto::RoiRange &roi : options.rois )
    {
      if( (se_energy >= roi.lower_energy) && (se_energy <= roi.upper_energy) )
      {
        have_se_roi = true;
        break;
      }
    }
    
    if( !have_se_roi )
    {
      RelActCalcAuto::RoiRange se_roi;
      // Set continuum/range-limits explicitly instead of leaving the RoiRange struct defaults.
      // The shared late-boundary policy decides any overlap; Fixed prevents a second generic
      // RelActAuto recombination after that decision.
      se_roi.continuum_type = PeakContinuum::OffsetType::Linear;
      se_roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
      se_roi.lower_energy = se_energy - roi_width_lower;
      se_roi.upper_energy = se_energy + roi_width_upper;
      
      // Make sure ROI is within valid energy range
      if( se_roi.lower_energy < min_valid_energy )
        se_roi.lower_energy = min_valid_energy;
      if( se_roi.upper_energy > max_valid_energy )
        se_roi.upper_energy = max_valid_energy;
      
      // Allow ROI to be as small as 0.5 FWHM on each side (1 FWHM total minimum)
      // If it overlaps with existing ROIs, we'll resolve that later before calling solve
      const double min_width_each_side = 0.5 * parent_peak->fwhm();
      
      // Shrink ROI if needed to maintain minimum width on each side
      if( (se_energy - se_roi.lower_energy) < min_width_each_side )
        se_roi.lower_energy = se_energy - min_width_each_side;
      if( (se_roi.upper_energy - se_energy) < min_width_each_side )
        se_roi.upper_energy = se_energy + min_width_each_side;
      
      // Constrain to valid energy range
      if( se_roi.lower_energy < min_valid_energy )
        se_roi.lower_energy = min_valid_energy;
      if( se_roi.upper_energy > max_valid_energy )
        se_roi.upper_energy = max_valid_energy;
      
      // Add the ROI - overlaps will be resolved later
      options.rois.push_back( se_roi );
      
      if( PEAK_FIT_DEBUG_PRINTOUT )
      {
        std::cout << "  Added ROI [" << se_roi.lower_energy << ", " << se_roi.upper_energy 
                  << "] keV for S.E. peak (may overlap, will be resolved)" << std::endl;
      }
    }
    
    // Check/add D.E. ROI
    bool have_de_roi = false;
    for( const RelActCalcAuto::RoiRange &roi : options.rois )
    {
      if( (de_energy >= roi.lower_energy) && (de_energy <= roi.upper_energy) )
      {
        have_de_roi = true;
        break;
      }
    }
    
    if( !have_de_roi )
    {
      RelActCalcAuto::RoiRange de_roi;
      // See the S.E. ROI above: set continuum/range-limits explicitly (not the struct defaults) so a
      // later merge with a source ROI inherits source-like settings.
      de_roi.continuum_type = PeakContinuum::OffsetType::Linear;
      de_roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
      de_roi.lower_energy = de_energy - roi_width_lower;
      de_roi.upper_energy = de_energy + roi_width_upper;
      
      // Make sure ROI is within valid energy range
      if( de_roi.lower_energy < min_valid_energy )
        de_roi.lower_energy = min_valid_energy;
      if( de_roi.upper_energy > max_valid_energy )
        de_roi.upper_energy = max_valid_energy;
      
      // Allow ROI to be as small as 0.5 FWHM on each side (1 FWHM total minimum)
      // If it overlaps with existing ROIs, we'll resolve that later before calling solve
      const double min_width_each_side = 0.5 * parent_peak->fwhm();
      
      // Shrink ROI if needed to maintain minimum width on each side
      if( (de_energy - de_roi.lower_energy) < min_width_each_side )
        de_roi.lower_energy = de_energy - min_width_each_side;
      if( (de_roi.upper_energy - de_energy) < min_width_each_side )
        de_roi.upper_energy = de_energy + min_width_each_side;
      
      // Constrain to valid energy range
      if( de_roi.lower_energy < min_valid_energy )
        de_roi.lower_energy = min_valid_energy;
      if( de_roi.upper_energy > max_valid_energy )
        de_roi.upper_energy = max_valid_energy;
      
      // Add the ROI - overlaps will be resolved later
      options.rois.push_back( de_roi );
      
      if( PEAK_FIT_DEBUG_PRINTOUT )
      {
        std::cout << "  Added ROI [" << de_roi.lower_energy << ", " << de_roi.upper_energy 
                  << "] keV for D.E. peak (may overlap, will be resolved)" << std::endl;
      }
    }
    
    // Add single escape floating peak if not already present
    // Use tighter tolerance based on parent FWHM
    const double fp_check_tolerance = 0.5 * parent_peak->fwhm();
    bool already_have_se = false;
    for( const RelActCalcAuto::FloatingPeak &fp : options.floating_peaks )
    {
      if( std::fabs( fp.energy - se_energy ) < fp_check_tolerance )
      {
        already_have_se = true;
        break;
      }
    }
    
    if( !already_have_se )
    {
      RelActCalcAuto::FloatingPeak se_fp;
      se_fp.energy = se_energy;
      se_fp.release_fwhm = false;
      se_fp.energy_origin = RelActCalcAuto::FloatingPeak::EnergyType::Known;
      options.floating_peaks.push_back( se_fp );
      
      if( PEAK_FIT_DEBUG_PRINTOUT )
      {
        std::cout << "  Added S.E. floating peak at " << se_energy << " keV" << std::endl;
      }
    }
    
    // Add double escape floating peak if not already present
    bool already_have_de = false;
    for( const RelActCalcAuto::FloatingPeak &fp : options.floating_peaks )
    {
      if( std::fabs( fp.energy - de_energy ) < fp_check_tolerance )
      {
        already_have_de = true;
        break;
      }
    }
    
    if( !already_have_de )
    {
      RelActCalcAuto::FloatingPeak de_fp;
      de_fp.energy = de_energy;
      de_fp.release_fwhm = false;
      de_fp.energy_origin = RelActCalcAuto::FloatingPeak::EnergyType::Known;
      options.floating_peaks.push_back( de_fp );
      
      if( PEAK_FIT_DEBUG_PRINTOUT )
      {
        std::cout << "  Added D.E. floating peak at " << de_energy << " keV" << std::endl;
      }
    }
  }//for( loop over escape peak candidates )
}//add_escape_peak_floating_peaks_if_appropriate(...)


/** Resolve any overlapping ROIs by merging them.

 This is called after escape peak ROIs are added, which may cause overlaps.  The input ROIs are
 otherwise already de-overlapped, so in practice every overlap involves an escape ROI.  Overlapping
 ROIs are merged (the lower-energy ROI is extended to cover both); there is intentionally no width
 cap and no split path (see the merge site for rationale).

 \param rois Vector of ROI ranges that may have overlaps - will be modified in place
 \param floating_peaks Vector of floating peaks; used by the developer-check to confirm overlaps are escape-related
 */
void resolve_overlapping_rois( std::vector<RelActCalcAuto::RoiRange> &rois,
                               const std::vector<RelActCalcAuto::FloatingPeak> &floating_peaks )
{
  if( rois.size() < 2 )
    return;
  
  // Sort by lower energy
  std::sort( rois.begin(), rois.end(), 
            []( const RelActCalcAuto::RoiRange &a, const RelActCalcAuto::RoiRange &b ) {
              return a.lower_energy < b.lower_energy;
            } );
  
  // Check for and resolve overlaps
  std::vector<RelActCalcAuto::RoiRange> resolved_rois;
  resolved_rois.reserve( rois.size() );
  
  for( size_t i = 0; i < rois.size(); ++i )
  {
    if( resolved_rois.empty() )
    {
      resolved_rois.push_back( rois[i] );
      continue;
    }
    
    RelActCalcAuto::RoiRange &last = resolved_rois.back();
    const RelActCalcAuto::RoiRange &current = rois[i];
    
    // Check for overlap
    if( current.lower_energy < last.upper_energy )
    {
      // We have an overlap
#if( PERFORM_DEVELOPER_CHECKS )
      // Sanity check: overlaps should only arise from escape-peak ROIs (the input ROIs are already
      // de-overlapped; add_escape_peak_floating_peaks_if_appropriate is the only thing that injects a
      // possibly-overlapping ROI before this call).  Escape ROIs are built around an escape floating
      // peak whose energy_origin is EnergyType::Known, so recognize the overlap as escape-related
      // when either ROI contains such a peak.  (This generalizes a previous check that hard-coded
      // only Th-232 2614's S.E./D.E. energies and would false-positive for any other escape parent.)
      const double last_width = last.upper_energy - last.lower_energy;
      const double current_width = current.upper_energy - current.lower_energy;
      const double overlap = last.upper_energy - current.lower_energy;

      auto roi_has_escape_floating_peak = [&]( const RelActCalcAuto::RoiRange &roi ) -> bool {
        for( const RelActCalcAuto::FloatingPeak &fp : floating_peaks )
        {
          if( (fp.energy >= roi.lower_energy) && (fp.energy <= roi.upper_energy)
              && (fp.energy_origin == RelActCalcAuto::FloatingPeak::EnergyType::Known) )
            return true;
        }
        return false;
      };

      const bool last_has_escape_peak = roi_has_escape_floating_peak( last );
      const bool current_has_escape_peak = roi_has_escape_floating_peak( current );

      if( !last_has_escape_peak && !current_has_escape_peak )
      {
        std::cerr << "WARNING: resolve_overlapping_rois found overlap NOT due to escape peaks:\n"
                  << "  ROI 1: [" << last.lower_energy << ", " << last.upper_energy
                  << "] keV (width=" << last_width << " keV)\n"
                  << "  ROI 2: [" << current.lower_energy << ", " << current.upper_energy
                  << "] keV (width=" << current_width << " keV)\n"
                  << "  Overlap: " << overlap << " keV\n"
                  << "  Neither ROI contains a Known-origin (escape) floating peak\n"
                  << "  This suggests ROI generation logic has a bug beyond escape peaks." << std::endl;
        assert( last_has_escape_peak || current_has_escape_peak );
      }
#endif

      // Merge the ROIs by extending the last one to include both.  We deliberately do NOT cap the
      // merged width or split the overlap: escape peaks are a rare exception (only added for
      // high-energy parents in the NORM-fit path) and occur at higher energies where ROIs tend to be
      // narrow, so an over-wide merge is not a practical concern here.
      last.upper_energy = std::max( last.upper_energy, current.upper_energy );
      
      if( PEAK_FIT_DEBUG_PRINTOUT )
      {
        std::cout << "Merged overlapping ROIs: [" << last.lower_energy << ", " 
                  << current.upper_energy << "] -> [" << last.lower_energy << ", " 
                  << last.upper_energy << "]" << std::endl;
      }
    }
    else
    {
      // No overlap - add as new ROI
      resolved_rois.push_back( current );
    }
  }
  
  rois = resolved_rois;
}//resolve_overlapping_rois(...)

void ensure_min_channel_gap(
    std::vector<RelActCalcAuto::RoiRange> &rois,
    const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal );


/** Central reconciliation for automatic ranges added after clustering.  The legacy resolver above
 remains available only to the transactional R6 path. */
void resolve_automatic_overlapping_rois(
    std::vector<RelActCalcAuto::RoiRange> &rois,
    const std::vector<RelActCalcAuto::FloatingPeak> &floating_peaks,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const detail::GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const PeakFitForNuclideConfig &config,
    const std::string &stage,
    std::vector<AutomaticRoiDecisionDiagnostic> *diagnostics,
    const std::vector<RelActCalcAuto::RoiRange> &protected_ranges = {},
    const std::vector<std::pair<double,double>> &modeled_peak_candidates = {},
    const bool use_automatic_roi_policy = true )
{
  if( rois.empty() )
    return;
  if( !use_automatic_roi_policy )
  {
    resolve_overlapping_rois( rois, floating_peaks );
    ensure_min_channel_gap( rois, foreground ? foreground->energy_calibration() : nullptr );
    return;
  }
  const auto is_protected = [&protected_ranges]( const RelActCalcAuto::RoiRange &roi ) {
    return std::any_of( std::begin(protected_ranges), std::end(protected_ranges),
      [&roi]( const RelActCalcAuto::RoiRange &protected_roi ) {
        return (std::fabs(roi.lower_energy - protected_roi.lower_energy) < 1.0e-6)
            && (std::fabs(roi.upper_energy - protected_roi.upper_energy) < 1.0e-6)
            && (roi.continuum_type == protected_roi.continuum_type);
      } );
  };
  std::sort( std::begin(rois), std::end(rois),
    []( const RelActCalcAuto::RoiRange &lhs, const RelActCalcAuto::RoiRange &rhs ) {
      return lhs.lower_energy < rhs.lower_energy;
    } );

  // POLICY MODE: reconcile the (possibly overlapping) ROIs through the atom-safe partition layer.
  // Each modeled candidate and floating feature is assigned EXACTLY ONCE to the input ROI that
  // contains it (nearest midpoint on ties) - replacing the former inclusive-containment scan that
  // could claim an atom in two overlapping ROIs - and the reconciler then folds adjacent ROIs
  // without any pop_back/skip drop path.
  const bool have_cal = foreground && foreground->energy_calibration()
      && foreground->energy_calibration()->valid() && (foreground->num_gamma_channels() > 0);
  if( !have_cal )
  {
    // No usable calibration for a channel-aligned partition; fall back to the legacy resolver.
    resolve_overlapping_rois( rois, floating_peaks );
    ensure_min_channel_gap( rois, foreground ? foreground->energy_calibration() : nullptr );
    return;
  }

  const size_t nroi = rois.size();
  const auto assign_to_roi = [&rois, nroi]( const double energy ) -> long {
    long best = -1;
    double best_dist = std::numeric_limits<double>::infinity();
    for( size_t j = 0; j < nroi; ++j )
    {
      if( (energy >= rois[j].lower_energy) && (energy <= rois[j].upper_energy) )
      {
        const double mid = 0.5 * (rois[j].lower_energy + rois[j].upper_energy);
        const double d = std::fabs( energy - mid );
        if( (best < 0) || (d < best_dist) ){ best = static_cast<long>(j); best_dist = d; }
      }
    }
    return best;
  };//assign_to_roi

  std::vector<std::vector<detail::RoiAtom>> roi_atoms( nroi );
  for( const std::pair<double,double> &peak : modeled_peak_candidates )
  {
    const long j = assign_to_roi( peak.first );
    if( j < 0 )
      continue;  // a candidate outside every ROI has no owner here (pre-existing floater)
    detail::RoiAtom atom;
    atom.id = detail::next_roi_atom_id();
    atom.energy = peak.first;
    atom.area = peak.second;
    atom.kind = detail::RoiAtomKind::ModeledGamma;
    roi_atoms[static_cast<size_t>(j)].push_back( atom );
  }
  for( const RelActCalcAuto::FloatingPeak &peak : floating_peaks )
  {
    const long j = assign_to_roi( peak.energy );
    if( j < 0 )
      continue;
    double area = 0.0;
    for( const std::shared_ptr<const PeakDef> &observed : unfit_auto_peaks )
    {
      if( !observed || !observed->gausPeak() )
        continue;
      const double local_fwhm = fwhm_at_energy( peak.energy );
      if( std::isfinite(local_fwhm) && (local_fwhm > 0.0)
          && (std::fabs(observed->mean() - peak.energy) <= 0.75*local_fwhm) )
        area = std::max( area, observed->amplitude() );
    }
    std::vector<detail::RoiAtom> &atoms = roi_atoms[static_cast<size_t>(j)];
    const std::vector<detail::RoiAtom>::iterator dup = std::find_if( std::begin(atoms),
        std::end(atoms), [&peak]( const detail::RoiAtom &a ){
          return std::fabs(a.energy - peak.energy) < 1.0e-6; } );
    if( dup == std::end(atoms) )
    {
      detail::RoiAtom atom;
      atom.id = detail::next_roi_atom_id();
      atom.energy = peak.energy;
      atom.area = area;
      atom.kind = detail::RoiAtomKind::FloatingFeature;
      atoms.push_back( atom );
    }else
    {
      dup->area = std::max( dup->area, area );  // keep the modeled atom identity, strongest area
    }
  }

  std::vector<detail::AutomaticRoiComponent> comps;
  comps.reserve( nroi );
  for( size_t j = 0; j < nroi; ++j )
  {
    detail::AutomaticRoiComponent c;
    c.lower = rois[j].lower_energy;
    c.upper = rois[j].upper_energy;
    c.first_channel = foreground->find_gamma_channel( static_cast<float>(rois[j].lower_energy) );
    c.last_channel = foreground->find_gamma_channel( static_cast<float>(rois[j].upper_energy) );
    c.continuum_type = rois[j].continuum_type;
    c.range_limits_type = rois[j].range_limits_type;
    c.protected_geometry = is_protected( rois[j] );
    c.joined_groups = 1;
    c.atoms = std::move( roi_atoms[j] );
    comps.push_back( std::move(c) );
  }

  // Unfit peaks coinciding with a modeled/floating atom are not interfering features.
  std::vector<std::shared_ptr<const PeakDef>> filtered_unfit;
  for( const std::shared_ptr<const PeakDef> &pk : unfit_auto_peaks )
  {
    if( !pk || !pk->gausPeak() )
      continue;
    bool matches = false;
    for( const detail::AutomaticRoiComponent &c : comps )
    {
      for( const detail::RoiAtom &a : c.atoms )
      {
        const double f = fwhm_at_energy( a.energy );
        if( std::isfinite(f) && (f > 0.0) && (std::fabs(pk->mean() - a.energy) <= 0.75*f) )
        {
          matches = true;
          break;
        }
      }
      if( matches )
        break;
    }
    if( !matches )
      filtered_unfit.push_back( pk );
  }

  detail::AutomaticRoiPolicySettings policy_settings;
  policy_settings.merge_tail_z = config.merge_tail_z;
  policy_settings.merge_clean_gap_fwhm = config.merge_clean_gap_fwhm;
  policy_settings.continuum_aicc_penalty = config.cont_order_aicc_penalty;
  policy_settings.peak_core_num_fwhm = config.auto_roi_core_num_fwhm;
  policy_settings.max_width_fwhm = config.auto_rel_eff_sol_max_fwhm;
  policy_settings.minimum_partition_gap_fwhm = config.auto_roi_partition_min_gap_fwhm;
  policy_settings.allow_clean_gap_partition_override
    = config.auto_roi_partition_allow_clean_gap_override;
  policy_settings.force_partition_gap_fwhm = config.auto_roi_partition_force_gap_fwhm;
  policy_settings.max_partition_children = config.auto_roi_partition_max_children;
  policy_settings.allow_overwide_overlap_partition = config.auto_roi_partition_overwide;
  policy_settings.stage = stage;

  detail::AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = foreground->gamma_channel_lower( 0 );
  cons.highest_energy = foreground->gamma_channel_upper( foreground->num_gamma_channels() - 1 );
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0;  // resolve imposes no minimum width (avoid geometry drift)
  cons.peak_core_num_fwhm = config.auto_roi_core_num_fwhm;

  const detail::AutomaticRoiReconcileResult rr = detail::reconcile_automatic_components_one_pass(
      std::move(comps), foreground, global_continuum, fwhm_at_energy, filtered_unfit,
      policy_settings, cons, diagnostics );

  if( !rr.valid )
  {
    // The transaction failed its own invariant check (should not happen; a dev build already
    // asserted).  Honor the all-or-nothing contract: retain incumbent geometry via the legacy
    // resolver, which still guarantees channel-disjoint output.
    resolve_overlapping_rois( rois, floating_peaks );
    ensure_min_channel_gap( rois, foreground ? foreground->energy_calibration() : nullptr );
    return;
  }

  std::vector<RelActCalcAuto::RoiRange> resolved;
  resolved.reserve( rr.components.size() );
  for( const detail::AutomaticRoiComponent &c : rr.components )
  {
    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = c.lower;
    roi.upper_energy = c.upper;
    roi.continuum_type = c.continuum_type;
    roi.range_limits_type = c.protected_geometry
        ? c.range_limits_type : RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    resolved.push_back( roi );
  }
  rois = std::move( resolved );

#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t index = 1; index < rois.size(); ++index )
  {
    assert( rois[index - 1].upper_energy <= rois[index].lower_energy );
    assert( is_protected(rois[index - 1])
        || (rois[index - 1].range_limits_type == RelActCalcAuto::RoiRange::RangeLimitsType::Fixed) );
  }
#endif
}//resolve_automatic_overlapping_rois(...)


/** Ensure at least one channel gap between consecutive ROIs.

 After resolving overlaps, adjacent ROIs may abut exactly (upper_energy of one equals
 lower_energy of the next), sharing the same channel boundary.  This function ensures
 at least one channel separates consecutive ROIs by adjusting their boundaries.

 ROIs must be sorted by lower_energy before calling this function.
 */
void ensure_min_channel_gap( std::vector<RelActCalcAuto::RoiRange> &rois,
                             const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal )
{
  if( !energy_cal || !energy_cal->valid() || (rois.size() < 2) )
    return;

  for( size_t i = 1; i < rois.size(); ++i )
  {
    const size_t prev_upper_ch = energy_cal->channel_for_energy( rois[i - 1].upper_energy );
    const size_t curr_lower_ch = energy_cal->channel_for_energy( rois[i].lower_energy );

    if( curr_lower_ch <= prev_upper_ch )
    {
      // Abutting or overlapping in channel space - create a 1-channel gap
      // Split at the midpoint energy, assigning mid_ch to prev and mid_ch+1 to current
      const size_t mid_ch = (prev_upper_ch + curr_lower_ch) / 2;
      const size_t num_channels = energy_cal->num_channels();

      if( (mid_ch + 1) < num_channels )
      {
        rois[i - 1].upper_energy = energy_cal->energy_for_channel( mid_ch );
        rois[i].lower_energy = energy_cal->energy_for_channel( mid_ch + 1 );
      }
    }
  }//for( size_t i = 1; i < rois.size(); ++i )
}//ensure_min_channel_gap(...)


/** Remove any floating peaks that are not covered by any ROI.

 This is called before each RelActCalcAuto::solve() invocation to guard against
 orphaned floating peaks inherited from a previous solution's options (e.g. when
 the iterative refinement loop replaces ROIs with re-clustered ones that no longer
 cover an earlier floating peak).  The solver throws if a floating peak has no
 covering ROI, so we clean up proactively here rather than duplicating the logic
 in every add_* helper.
 */
void remove_floating_peaks_without_roi( RelActCalcAuto::Options &options )
{
  auto &fps = options.floating_peaks;
  fps.erase( std::remove_if( begin(fps), end(fps),
    [&]( const RelActCalcAuto::FloatingPeak &fp ) -> bool {
      for( const RelActCalcAuto::RoiRange &roi : options.rois )
      {
        if( (fp.energy >= roi.lower_energy) && (fp.energy <= roi.upper_energy) )
          return false;
      }
      return true;  // no ROI covers this floating peak – remove it
    } ), end(fps) );
}//remove_floating_peaks_without_roi(...)


/** Assign escape peak relationships for high-energy gamma lines if appropriate.
 
 This function checks if fit peaks contain escape peaks that were added as floating peaks,
 and if they're significant, assigns them as single or double escape peaks of their parent.
 Currently handles Th232 2614 keV (and structured to support Ra226 lines in the future).
 
 \param fit_peaks Vector of fit peaks to potentially modify with escape peak assignments
 \param fit_norm_peaks Whether NORM background peaks were fit
 */
void assign_escape_peak_relationships(
  std::vector<PeakDef> &fit_peaks,
  const bool fit_norm_peaks,
  const PeakFitUtils::CoarseResolutionType det_type )
{
  // Only assign escape peaks for high-resolution detectors (HPGe)
  if( !fit_norm_peaks || (det_type != PeakFitUtils::CoarseResolutionType::High) )
    return;
  
  const double electron_rest_mass = 510.9989; // keV
  const double single_escape_offset = electron_rest_mass;
  const double double_escape_offset = 2.0 * electron_rest_mass;
  
  // Define high-energy gamma lines that commonly have escape peaks
  struct EscapePeakCandidate
  {
    double parent_energy;
    std::string parent_symbol;  // e.g., "Th232" for Tl208 in Th232 decay chain
  };
  
  const std::vector<EscapePeakCandidate> candidates = {
    { 2614.533, "Th232" },  // Th232 (Tl208) 2614 keV
    // Future candidates:
    // { 2204.21, "Ra226" },  // Ra226 (Bi214) 2204 keV
    // { 1764.49, "Ra226" },  // Ra226 (Bi214) 1764 keV
  };
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    return;
  
  // For each high-energy candidate, look for its parent peak and potential escape peaks
  for( const EscapePeakCandidate &candidate : candidates )
  {
    const double se_energy = candidate.parent_energy - single_escape_offset;
    const double de_energy = candidate.parent_energy - double_escape_offset;
    
    // Find the parent peak (2614 keV for Th232)
    PeakDef *parent_peak = nullptr;
    for( PeakDef &peak : fit_peaks )
    {
      if( std::fabs( peak.mean() - candidate.parent_energy ) < 2.0 ) // 2 keV tolerance
      {
        parent_peak = &peak;
        break;
      }
    }
    
    if( !parent_peak )
      continue;  // No parent peak found
    
    // Get the nuclide for this parent
    const SandiaDecay::Nuclide * const parent_nuclide = db->nuclide( candidate.parent_symbol );
    if( !parent_nuclide )
      continue;
    
    // Check if parent peak already has transition assigned - if so, use it
    // Otherwise, search for the transition
    const SandiaDecay::Transition *parent_transition = parent_peak->nuclearTransition();
    int parent_particle_index = parent_peak->decayParticleIndex();
    
    if( !parent_transition || (parent_particle_index < 0) )
    {
      // Need to find the transition manually by searching descendants
      const std::vector<const SandiaDecay::Nuclide *> descendants = parent_nuclide->descendants();
      
      for( const SandiaDecay::Nuclide * const nuc : descendants )
      {
        if( !nuc )
          continue;
        
        for( const SandiaDecay::Transition * const transition : nuc->decaysToChildren )
        {
          if( !transition )
            continue;
          
          for( size_t i = 0; i < transition->products.size(); ++i )
          {
            const SandiaDecay::RadParticle &particle = transition->products[i];
            if( (particle.type == SandiaDecay::GammaParticle)
               && (std::fabs( particle.energy - candidate.parent_energy ) < 0.5) ) // 0.5 keV tolerance
            {
              parent_transition = transition;
              parent_particle_index = static_cast<int>( i );
              break;
            }
          }
          if( parent_transition )
            break;
        }
        if( parent_transition )
          break;
      }
      
      if( !parent_transition || (parent_particle_index < 0) )
        continue;  // Couldn't find the transition
    }
    
    // Look for single escape peak
    for( PeakDef &peak : fit_peaks )
    {
      if( std::fabs( peak.mean() - se_energy ) < 2.0 ) // 2 keV tolerance
      {
        // Check if peak is significant (area > some threshold)
        const double peak_area = peak.peakArea();
        const double peak_area_uncert = peak.peakAreaUncert();
        // peakAreaUncert() returns the -1 sentinel when uncertainty was never computed; fall back to
        // a Poisson sqrt(area) so a real escape peak is not silently dropped.
        const double significance = (peak_area_uncert > 0.0)
          ? (peak_area / peak_area_uncert)
          : ((peak_area > 0.0) ? std::sqrt(peak_area) : 0.0);
        
        if( significance > 3.0 ) // At least 3-sigma significance
        {
          // Assign as single escape peak
          peak.setNuclearTransition( parent_nuclide, parent_transition, parent_particle_index,
                                     PeakDef::SourceGammaType::SingleEscapeGamma );
          
          if( PEAK_FIT_DEBUG_PRINTOUT )
          {
            std::cout << "Assigned peak at " << peak.mean() << " keV as S.E. of "
                      << candidate.parent_symbol << " " << candidate.parent_energy << " keV"
                      << std::endl;
          }
        }
        break;
      }
    }
    
    // Look for double escape peak
    for( PeakDef &peak : fit_peaks )
    {
      if( std::fabs( peak.mean() - de_energy ) < 2.0 ) // 2 keV tolerance
      {
        // Check if peak is significant (area > some threshold)
        const double peak_area = peak.peakArea();
        const double peak_area_uncert = peak.peakAreaUncert();
        // peakAreaUncert() returns the -1 sentinel when uncertainty was never computed; fall back to
        // a Poisson sqrt(area) so a real escape peak is not silently dropped.
        const double significance = (peak_area_uncert > 0.0)
          ? (peak_area / peak_area_uncert)
          : ((peak_area > 0.0) ? std::sqrt(peak_area) : 0.0);
        
        if( significance > 3.0 ) // At least 3-sigma significance
        {
          // Assign as double escape peak
          peak.setNuclearTransition( parent_nuclide, parent_transition, parent_particle_index,
                                     PeakDef::SourceGammaType::DoubleEscapeGamma );
          
          if( PEAK_FIT_DEBUG_PRINTOUT )
          {
            std::cout << "Assigned peak at " << peak.mean() << " keV as D.E. of "
                      << candidate.parent_symbol << " " << candidate.parent_energy << " keV"
                      << std::endl;
          }
        }
        break;
      }
    }
  }//for( loop over escape peak candidates )
}//assign_escape_peak_relationships(...)


std::pair<double,double> find_valid_energy_range( const std::shared_ptr<const SpecUtils::Measurement> &meas )
{
  // This implementation is an updated implementation of `find_spectroscopic_extent(...)`,
  //  and that function should probably be updated
  if( !meas 
    || !meas->energy_calibration() 
    || !meas->energy_calibration()->valid()
    || (meas->num_gamma_channels() < 7) )
  {
    return {0.0, 0.0};
  }
  
  size_t lower_channel = 0, upper_channel = 0;

  const std::vector<float> &channel_counts = *meas->gamma_counts();
  const size_t nbin = channel_counts.size();
  
  //First detect where spectrum begins
  const int side_bins = 3;
  const int order = 2;
  const int derivative = 2;
  std::vector<float> smoothed_2nd, smoothed_2nd_variance;
  SavitzyGolayCoeffs sgcoeffs( side_bins, side_bins, order, derivative );
  sgcoeffs.smooth_with_variance( channel_counts, smoothed_2nd, smoothed_2nd_variance );

  // Lower extent: delegate to the unified resolution-aware estimator ("C1") now living in
  //  ExperimentalPeakSearch::find_spectroscopic_extent(); see its definition for the method.
  //  This function keeps its own upper-extent logic (below), which callers rely on.
  {
    size_t lo_ch = 0, hi_ch_ignored = 0;
    if( ExperimentalPeakSearch::find_spectroscopic_extent( meas, lo_ch, hi_ch_ignored ) )
      lower_channel = lo_ch;
    else
      lower_channel = nbin;  // force the simple fallback below (lower_channel > nbin/3)
  }


  // Find the upper energy limit by looking for the last channel with meaningful signal
  // Start from the end and work backwards
  size_t upperlastchannel = nbin - 1;

  // First, skip trailing channels with very low counts (< 3 counts)
  while( (upperlastchannel > 0) && (channel_counts[upperlastchannel] < 3.0f) )
    --upperlastchannel;

  // Now increase upperlastchannel until we find at least two consecutive zero-count channels
  bool found_two_zeros = false;
  while( (upperlastchannel < nbin - 2) && !found_two_zeros )
  {
    ++upperlastchannel;

    // Check if we have two consecutive zeros starting at upperlastchannel
    if( (channel_counts[upperlastchannel] == 0.0f)
        && (channel_counts[upperlastchannel + 1] == 0.0f) )
    {
      found_two_zeros = true;
      ++upperlastchannel; // Set to one more than the second zero
    }
  }

  // Ensure upperlastchannel is valid
  upperlastchannel = std::min( upperlastchannel, nbin - 1 );


  //  upperlastchannel = (std::min)( upperlastchannel + 5, nbin - 1 );
  
  // Search backwards from the end to find the last peak-like structure
  // Look for significant negative curvature, then find where it becomes positive on the right side
  size_t lastchannel = upperlastchannel;

  bool found_last_signal = false;
  for( size_t i = upperlastchannel; i > lower_channel; --i )
  {
    // Look for statistically significant negative curvature (peak shape)
    if( (i >= smoothed_2nd.size()) || (smoothed_2nd_variance[i] <= 0.0f) )
      continue;

    const float sigma_2nd = std::sqrt( smoothed_2nd_variance[i] );
    if( smoothed_2nd[i] < (-2.0f * sigma_2nd) )
    {
      // Found significant negative curvature - this is a potential peak
      // Now find where the second derivative becomes positive on the right side
      size_t right_extent = i;
      for( size_t j = i + 1; j < smoothed_2nd.size(); ++j )
      {
        if( smoothed_2nd[j] > 0.0f )
        {
          right_extent = j;
          break;
        }
      }

      // Calculate peak width and add ~25% buffer
      const size_t peak_width = (right_extent > i) ? (right_extent - i) : 1;
      const size_t buffer = std::max( size_t(4), peak_width / 4 );
      lastchannel = std::min( right_extent + buffer, nbin - 1 );
      found_last_signal = true;
      break;
    }
  }

  // If we didn't find peak structure in the backwards search, use the simple upperlastchannel
  if( !found_last_signal )
    lastchannel = upperlastchannel;

  // If the sophisticated detection failed, fall back to a simple algorithm
  if( (lower_channel > (nbin/3)) || (lastchannel <= lower_channel) )
  {
    // Find first channel where this channel and the next three all have non-zero counts
    size_t first_nonzero = 0;
    while( first_nonzero + 3 < nbin )
    {
      if( (channel_counts[first_nonzero] > 0.0f)
         && (channel_counts[first_nonzero + 1] > 0.0f)
         && (channel_counts[first_nonzero + 2] > 0.0f)
         && (channel_counts[first_nonzero + 3] > 0.0f) )
      {
        break;
      }
      ++first_nonzero;
    }

    // Find last channel where this channel and the previous channel both have non-zero counts
    size_t last_nonzero = nbin - 1;
    while( last_nonzero > 0 )
    {
      if( (channel_counts[last_nonzero] > 0.0f) && (channel_counts[last_nonzero - 1] > 0.0f) )
        break;
      --last_nonzero;
    }

    // Make sure we found valid channels
    if( first_nonzero < nbin && last_nonzero > first_nonzero )
    {
      lower_channel = first_nonzero;
      upper_channel = last_nonzero;
    }
    else
    {
      // Ultimate fallback - use entire spectrum
      lower_channel = 0;
      upper_channel = nbin - 1;
    }

    return { meas->gamma_channel_lower(lower_channel), meas->gamma_channel_upper(upper_channel) };
  }//

  // Add some buffer beyond the last detected signal to ensure we don't cut off peak tails
  const size_t buffer_channels = std::max( size_t(10), size_t(0.005 * nbin) );
  upper_channel = std::min( lastchannel + buffer_channels, nbin - 1 );

  // Make sure we don't cut off too early - at minimum use upperlastchannel
  upper_channel = std::max( upper_channel, upperlastchannel );
  upper_channel = std::min( upper_channel, nbin - 1 );

  return { meas->gamma_channel_lower(lower_channel), meas->gamma_channel_upper(upper_channel) };
}//find_valid_energy_range


// Physically-meaningful low-energy floor for the analysis range.  find_valid_energy_range() reports
// where the spectrum has DATA, which for simulated spectra can extend down to a few keV - well below
// where the detector response / rel-eff are valid (e.g. the generic NaI DRF only spans 45-3000 keV).
// Analyzing there makes RelActCalcAuto evaluate an out-of-range rel-eff (inf/NaN) and, for the wide
// low-energy ROIs, dilutes peak significance.  Clamp the low bound to a per-detector floor (20 keV
// for HPGe, 25 keV otherwise), unless the supplied DRF is explicitly valid below that.
//
// Note this floor is the same for both rel-eff forms: for the PHYSICAL model, gammas below the DRF's
// valid range additionally yield a non-finite rel-eff and are skipped in the cost function (so the
// effective physical floor is the DRF's lower energy); for EMPIRICAL forms the polynomial rel-eff
// stays finite, so source lines down to this floor are kept.
double low_energy_analysis_floor( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                  const PeakFitUtils::CoarseResolutionType det_type )
{
  const double abs_floor = (det_type == PeakFitUtils::CoarseResolutionType::High) ? 20.0 : 25.0;
  if( drf && drf->isValid() && (drf->lowerEnergy() > 0.0) && (drf->lowerEnergy() < abs_floor) )
    return drf->lowerEnergy();
  return abs_floor;
}//low_energy_analysis_floor


bool rois_are_similar( const std::vector<RelActCalcAuto::RoiRange> &a,
                       const std::vector<RelActCalcAuto::RoiRange> &b,
                       const double energy_tolerance = 1.0 )
{
  if( a.size() != b.size() )
    return false;

  for( size_t i = 0; i < a.size(); ++i )
  {
    if( std::fabs(a[i].lower_energy - b[i].lower_energy) > energy_tolerance )
      return false;
    if( std::fabs(a[i].upper_energy - b[i].upper_energy) > energy_tolerance )
      return false;
  }
  return true;
}//rois_are_similar


RoiSignificanceResult compute_roi_chi2_significance(
  const RelActCalcAuto::RoiRange &roi,
  const std::vector<PeakDef> &all_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const double min_roi_significance_z,
  const bool include_peak_count_significance = true,
  const bool same_continuum_family_for_null = false )
{
  RoiSignificanceResult result;

  // Find peaks in this ROI
  std::vector<std::shared_ptr<const PeakDef>> peaks_in_roi;
  for( const PeakDef &peak : all_peaks )
  {
    if( peak.mean() >= roi.lower_energy && peak.mean() <= roi.upper_energy )
      peaks_in_roi.push_back( std::make_shared<PeakDef>( peak ) );
  }

  if( peaks_in_roi.empty() )
    return result;

  // Get channel range
  const size_t start_channel = data->find_gamma_channel( static_cast<float>( roi.lower_energy ) );
  const size_t end_channel = data->find_gamma_channel( static_cast<float>( roi.upper_energy ) );
  result.num_channels = (end_channel > start_channel) ? (end_channel - start_channel) : 0;

  if( result.num_channels < 3 )
    return result;


  // Get continuum parameters from first peak
  const std::shared_ptr<const PeakContinuum> continuum = peaks_in_roi[0]->continuum();
  if( !continuum )
    return result;

  const PeakContinuum::OffsetType cont_type = continuum->type();

  // Prepare data arrays for fit_amp_and_offset
  const std::vector<float> &channel_energies = *data->channel_energies();
  std::vector<float> channel_counts( result.num_channels );
  for( size_t i = 0; i < result.num_channels; ++i )
    channel_counts[i] = data->gamma_channel_content( start_channel + i );

  std::vector<double> empty_means, empty_sigmas;
  std::vector<double> continuum_coeffs, continuum_uncerts;
  std::vector<double> dummy_amps, dummy_amp_uncerts;

  // To compare apples to apples, call `fit_amp_and_offset(...)` for `peaks_in_roi`,
  // but with the amplitudes of the peaks fixed, so the continuum is re-fit, but not the peak amplitudes.
  // This way we can use the chi2 computed by fit_amp_and_offset, to consistently compare against
  // the chi2 computed for the continuum fit for no peaks.
  std::vector<PeakDef> fixed_peaks;
  for( const std::shared_ptr<const PeakDef> &peak : peaks_in_roi )
    fixed_peaks.push_back( *peak );

  result.chi2_with_peaks = fit_amp_and_offset(
    &channel_energies[start_channel],
    channel_counts.data(),
    result.num_channels,
    cont_type,
    continuum->referenceEnergy(),
    empty_means,
    empty_sigmas,
    fixed_peaks,  // Pass peaks as fixed - their amplitudes won't be fit
    PeakDef::SkewType::NoSkew,
    nullptr,
    dummy_amps,
    continuum_coeffs,
    dummy_amp_uncerts,
    continuum_uncerts
  );

  // Fit continuum only (no peaks) - use quadratic so the null hypothesis
  // has enough flexibility to model smooth curvature; otherwise a linear
  // continuum may make a continuum feature look like a significant peak.
  std::vector<PeakDef> empty_fixed_peaks;

  const PeakContinuum::OffsetType no_peak_cont_type = same_continuum_family_for_null
    ? cont_type : PeakContinuum::OffsetType::Quadratic;

  result.chi2_continuum_only = fit_amp_and_offset(
    &channel_energies[start_channel],
    channel_counts.data(),
    result.num_channels,
    no_peak_cont_type,
    continuum->referenceEnergy(),
    empty_means,
    empty_sigmas,
    empty_fixed_peaks,
    PeakDef::SkewType::NoSkew,
    nullptr,
    dummy_amps,
    continuum_coeffs,
    dummy_amp_uncerts,
    continuum_uncerts
  );

  // Compute chi2 reduction
  result.chi2_reduction = result.chi2_continuum_only - result.chi2_with_peaks;

  // Likelihood-ratio test (Wilks): refer the chi2 improvement from adding the ROI's peaks to a
  // chi2 distribution with dof = number of peaks (one amplitude each), and convert the survival
  // probability to a one-sided normal quantile so a single threshold works for any peak count.
  // The amplitudes were fit by the enclosing RelActAuto solve rather than per-ROI maximum
  // likelihood, and low-count channels are Poisson-not-Gaussian, so this z is a calibrated
  // ranking statistic rather than an exact p-value - the GA-tuned threshold absorbs the
  // miscalibration.
  {
    const size_t num_peak_dof = peaks_in_roi.size();
    if( (result.chi2_reduction > 0.0) && (num_peak_dof > 0) )
    {
      const boost::math::chi_squared_distribution<double> chi2_dist( static_cast<double>(num_peak_dof) );
      const double p_value = boost::math::cdf( boost::math::complement( chi2_dist, result.chi2_reduction ) );

      if( p_value < 1.0e-300 )
      {
        result.equivalent_z = 40.0;  // p underflows; overwhelmingly significant
      }else if( p_value > (1.0 - 1.0e-12) )
      {
        // chi2 improvement indistinguishable from zero: cdf rounds to exactly 1.0, where the
        // normal quantile throws an erfc_inv overflow (which would fail the whole fit).
        result.equivalent_z = -40.0;
      }else
      {
        const boost::math::normal_distribution<double> gaus_dist;
        result.equivalent_z = -boost::math::quantile( gaus_dist, p_value );
      }
    }
  }

  // Compute peak significance for each peak: peak_area / sqrt(peak_area + continuum)
  // This is a Poisson detection significance over ±1 FWHM that includes both the
  // continuum's noise and the peak's own Poisson noise.  A bare peak/sqrt(continuum)
  // form is overoptimistic when the peak is comparable to or larger than the
  // continuum (e.g. high-energy ROIs with very low background), where the peak's
  // own Poisson noise is the dominant uncertainty contribution.
  const double fwhm_coverage_fraction = 0.9793;

  for( const std::shared_ptr<const PeakDef> &peak : peaks_in_roi )
  {
    const double peak_mean = peak->mean();
    const double peak_fwhm = peak->fwhm();
    const double peak_lower = peak_mean - peak_fwhm;
    const double peak_upper = peak_mean + peak_fwhm;

    // Integrate continuum between peak_lower and peak_upper
    const std::shared_ptr<const PeakContinuum> peak_cont = peak->continuum();
    if( !peak_cont )
      continue;

    const double cont_integral_raw
      = peak_cont->offset_integral( peak_lower, peak_upper, data, peaks_in_roi );
    const double continuum_integral = std::max( 0.0, cont_integral_raw );

    // Peak amplitude in the ±1 FWHM range
    const double peak_amp_in_range = std::max( 0.0, peak->amplitude() * fwhm_coverage_fraction );

    // Total counts under +/- 1 FWHM (peak + continuum), floor 1.0 to avoid div-by-zero.
    const double total_in_range = std::max( 1.0, peak_amp_in_range + continuum_integral );
    const double peak_sig = peak_amp_in_range / std::sqrt( total_in_range );

    if( peak_sig > result.max_peak_significance )
      result.max_peak_significance = peak_sig;
  }//for( loop over peaks in ROI )

  // Keep an ROI if EITHER the whole-ROI Wilks significance (equivalent_z) OR its single strongest
  // peak's per-(+/-1 FWHM) data significance (max_peak_significance, computed just above but
  // historically used only for debug) clears the bar.  equivalent_z alone DILUTES a genuinely strong
  // peak when the ROI carries many rel-eff-tied peaks - its dof is the peak count, so a real 6-7 sigma
  // line inside a ~95-peak NORM-chain ROI (e.g. Th232) is referred to a chi^2_95 null and scores
  // NEGATIVE, discarding it.  max_peak_significance is the honest +/-1 FWHM Poisson significance (the
  // same statistic the injected truth reports as NSigmaOverBkg) and cannot pass junk: it keeps the ROI
  // only if its strongest peak clears the threshold against the fitted continuum, and the per-peak
  // observable refit still drops the weak siblings inside the ROI.  [architecture review, 2026-07-18]
  result.has_significant_peaks = (result.equivalent_z >= min_roi_significance_z)
      || (include_peak_count_significance
          && (result.max_peak_significance >= min_roi_significance_z));

#if( PERFORM_DEVELOPER_CHECKS )
  if( should_debug_print() )
  {
    std::cout << "compute_roi_chi2_significance: ROI [" << roi.lower_energy << ", " << roi.upper_energy << "] keV"
         << ", nch=" << result.num_channels
         << ", chi2_reduction=" << result.chi2_reduction
         << " (" << peaks_in_roi.size() << " peak dof)"
         << ", equivalent_z=" << result.equivalent_z
         << " (need " << min_roi_significance_z << ")"
         << ", max_peak_sig=" << result.max_peak_significance
         << ", significant=" << result.has_significant_peaks
         << std::endl;
  }
#endif

  return result;
}//compute_roi_chi2_significance


/** Evaluate one fitted ROI model over an explicit, immutable channel domain.

 `PeakFit::chi2_for_region` intentionally follows a peak continuum's fitted energy range and can
 therefore ignore caller-supplied channels.  That is useful for ordinary reporting, but invalid for
 the R6 nested comparison because the source-only and nuisance fits may have independently adjusted
 energy calibrations.  This helper selects the fitted continuum that best corresponds to the target
 ROI, then evaluates both its continuum and all peaks sharing it over exactly [lower_channel,
 upper_channel].
 */
FixedRoiModelScore fixed_roi_model_score(
  const std::vector<std::shared_ptr<const PeakDef>> &model_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const size_t lower_channel,
  const size_t upper_channel,
  const double target_lower_energy,
  const double target_upper_energy )
{
  FixedRoiModelScore result;
  if( !data || !data->channel_energies() || model_peaks.empty()
      || !(target_upper_energy > target_lower_energy)
      || (lower_channel >= data->num_gamma_channels())
      || (upper_channel < lower_channel) )
    return result;

  struct ContinuumGroup
  {
    std::shared_ptr<const PeakContinuum> continuum;
    std::vector<std::shared_ptr<const PeakDef>> peaks;
    double overlap = 0.0;
    double center_distance = std::numeric_limits<double>::max();
  };

  std::vector<ContinuumGroup> groups;
  for( const std::shared_ptr<const PeakDef> &peak : model_peaks )
  {
    if( !peak || !peak->continuum() )
      continue;
    const std::shared_ptr<const PeakContinuum> continuum = peak->continuum();
    auto pos = std::find_if( std::begin(groups), std::end(groups),
      [&continuum]( const ContinuumGroup &group ) {
        return group.continuum == continuum;
      } );
    if( pos == std::end(groups) )
    {
      ContinuumGroup group;
      group.continuum = continuum;
      groups.push_back( std::move(group) );
      pos = std::prev( std::end(groups) );
    }
    pos->peaks.push_back( peak );
  }

  const double target_center = 0.5 * (target_lower_energy + target_upper_energy);
  for( ContinuumGroup &group : groups )
  {
    double lower = target_lower_energy;
    double upper = target_upper_energy;
    if( group.continuum->energyRangeDefined() )
    {
      lower = group.continuum->lowerEnergy();
      upper = group.continuum->upperEnergy();
    }else
    {
      lower = std::numeric_limits<double>::max();
      upper = -std::numeric_limits<double>::max();
      for( const std::shared_ptr<const PeakDef> &peak : group.peaks )
      {
        lower = std::min( lower, peak->lowerX() );
        upper = std::max( upper, peak->upperX() );
      }
    }
    group.overlap = std::max( 0.0, std::min(upper, target_upper_energy)
                                   - std::max(lower, target_lower_energy) );
    group.center_distance = std::fabs( 0.5 * (lower + upper) - target_center );
  }

  const auto best = std::max_element( std::begin(groups), std::end(groups),
    []( const ContinuumGroup &lhs, const ContinuumGroup &rhs ) {
      if( lhs.overlap != rhs.overlap )
        return lhs.overlap < rhs.overlap;
      return lhs.center_distance > rhs.center_distance;
    } );
  if( (best == std::end(groups)) || !(best->overlap > 0.0) )
    return result;

  const size_t last_channel = std::min( upper_channel, data->num_gamma_channels() - 1 );
  const size_t num_channels = 1 + last_channel - lower_channel;
  const std::vector<float> &energies = *data->channel_energies();
  if( energies.size() <= last_channel + 1 )
    return result;

  std::vector<double> expected( num_channels, 0.0 );
  best->continuum->offset_integral( &energies[lower_channel], expected.data(), num_channels,
                                    data, best->peaks );
  for( const std::shared_ptr<const PeakDef> &peak : best->peaks )
    peak->gauss_integral( &energies[lower_channel], expected.data(), num_channels );

  double deviance = 0.0;
  for( size_t i = 0; i < num_channels; ++i )
  {
    const double observed = std::max( 0.0,
        static_cast<double>(data->gamma_channel_content(lower_channel + i)) );
    double predicted = expected[i];
    if( !std::isfinite(predicted) || (predicted < -1.0e-6) )
      return result;
    predicted = std::max( predicted, 1.0e-9 );
    deviance += (observed > 0.0)
      ? 2.0 * (predicted - observed + observed * std::log(observed / predicted))
      : 2.0 * predicted;
  }

  result.poisson_deviance = deviance;
  result.num_channels = num_channels;
  result.valid = std::isfinite(deviance);
  return result;
}//fixed_roi_model_score(...)


size_t solution_data_row_count( const RelActCalcAuto::RelActAutoSolution &solution )
{
  const std::shared_ptr<const SpecUtils::Measurement> &data = solution.m_foreground;
  if( !data || (data->num_gamma_channels() == 0) )
    return 0;
  const bool have_spec_cal_rois
    = (solution.m_final_roi_ranges_in_spectrum_cal.size() == solution.m_final_roi_ranges.size());
  const std::vector<RelActCalcAuto::RoiRange> &rois = have_spec_cal_rois
    ? solution.m_final_roi_ranges_in_spectrum_cal : solution.m_final_roi_ranges;
  size_t rows = 0;
  for( const RelActCalcAuto::RoiRange &roi : rois )
  {
    const size_t first = data->find_gamma_channel( static_cast<float>(roi.lower_energy) );
    const size_t last = std::min( data->find_gamma_channel(
        static_cast<float>(roi.upper_energy) ), data->num_gamma_channels() - 1 );
    if( last >= first )
      rows += last - first + 1;
  }
  return rows;
}


/** Data-only AICc on immutable, FWHM-local channel windows around the same source anchors.
 The effective parameter count is derived from the solution's data rows and m_dof_data. */
double source_anchor_data_aicc(
  const RelActCalcAuto::RelActAutoSolution &solution,
  const std::vector<RelActCalcManual::GenericPeakInfo> &anchors,
  const double aicc_penalty )
{
  const std::shared_ptr<const SpecUtils::Measurement> &data = solution.m_foreground;
  if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status)
      || !data || anchors.empty() || !(aicc_penalty > 0.0) )
    return std::numeric_limits<double>::max();

  std::vector<std::shared_ptr<const PeakDef>> model_peaks;
  model_peaks.reserve( solution.m_peaks_without_back_sub.size() );
  for( const PeakDef &peak : solution.m_peaks_without_back_sub )
    model_peaks.push_back( std::make_shared<const PeakDef>(peak) );

  double deviance = 0.0;
  size_t common_rows = 0;
  for( const RelActCalcManual::GenericPeakInfo &anchor : anchors )
  {
    const double half_width = 0.5 * anchor.m_fwhm;
    if( !(half_width > 0.0) )
      return std::numeric_limits<double>::max();
    const double lower = anchor.m_energy - half_width;
    const double upper = anchor.m_energy + half_width;
    const size_t first = data->find_gamma_channel( static_cast<float>(lower) );
    const size_t last = std::min( data->find_gamma_channel(
        static_cast<float>(upper) ), data->num_gamma_channels() - 1 );
    const FixedRoiModelScore score = fixed_roi_model_score(
        model_peaks, data, first, last, lower, upper );
    if( !score.valid )
      return std::numeric_limits<double>::max();
    deviance += score.poisson_deviance;
    common_rows += score.num_channels;
  }

  const size_t all_data_rows = solution_data_row_count( solution );
  if( (all_data_rows < solution.m_dof_data) || (common_rows == 0) )
    return std::numeric_limits<double>::max();
  const size_t num_parameters = all_data_rows - solution.m_dof_data;
  return detail::data_only_aicc(
      deviance, common_rows, num_parameters, aicc_penalty );
}


std::vector<PeakDef> distinct_significant_requested_source_peaks(
  const RelActCalcAuto::RelActAutoSolution &fit,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const double minimum_z,
  const std::function<double(double)> &fwhm_at_energy )
{
  std::vector<const PeakDef *> candidates;
  for( const PeakDef &peak : fit.m_peaks_without_back_sub )
  {
    const bool requested = std::any_of( std::begin(sources), std::end(sources),
      [&peak]( const RelActCalcAuto::NucInputInfo &source ) {
        return (peak.parentNuclide()
            && (peak.parentNuclide() == RelActCalcAuto::nuclide(source.source)))
          || (peak.xrayElement()
            && (peak.xrayElement() == RelActCalcAuto::element(source.source)))
          || (peak.reaction()
            && (peak.reaction() == RelActCalcAuto::reaction(source.source)));
      } );
    const double uncertainty = peak.amplitudeUncert();
    const double significance = (uncertainty > 0.0)
      ? (peak.amplitude() / uncertainty)
      : ((peak.amplitude() > 0.0) ? std::sqrt(peak.amplitude()) : 0.0);
    if( requested && (significance >= minimum_z) )
      candidates.push_back( &peak );
  }
  std::sort( std::begin(candidates), std::end(candidates),
    []( const PeakDef *lhs, const PeakDef *rhs ) {
      const double lhs_uncertainty = lhs->amplitudeUncert();
      const double rhs_uncertainty = rhs->amplitudeUncert();
      const double lhs_significance = (lhs_uncertainty > 0.0)
        ? (lhs->amplitude() / lhs_uncertainty) : std::sqrt(lhs->amplitude());
      const double rhs_significance = (rhs_uncertainty > 0.0)
        ? (rhs->amplitude() / rhs_uncertainty) : std::sqrt(rhs->amplitude());
      return lhs_significance > rhs_significance;
    } );
  std::vector<const PeakDef *> distinct;
  for( const PeakDef *candidate : candidates )
  {
    const bool overlaps = std::any_of( std::begin(distinct), std::end(distinct),
      [candidate, &fwhm_at_energy]( const PeakDef *existing ) {
        return std::fabs(candidate->mean() - existing->mean())
            < 0.5*(fwhm_at_energy(candidate->mean()) + fwhm_at_energy(existing->mean()));
      } );
    if( !overlaps )
      distinct.push_back( candidate );
  }
  std::vector<PeakDef> result;
  result.reserve( distinct.size() );
  for( const PeakDef *peak : distinct )
    result.push_back( *peak );
  return result;
}


size_t significant_requested_source_anchor_count(
  const RelActCalcAuto::RelActAutoSolution &candidate,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::vector<RelActCalcManual::GenericPeakInfo> &anchors,
  const double minimum_z,
  const std::function<double(double)> &fwhm_at_energy )
{
  const std::vector<PeakDef> candidate_lines
    = distinct_significant_requested_source_peaks(
        candidate, sources, minimum_z, fwhm_at_energy );
  return static_cast<size_t>( std::count_if( std::begin(anchors), std::end(anchors),
    [&candidate_lines, &fwhm_at_energy](
        const RelActCalcManual::GenericPeakInfo &anchor ) {
      const double anchor_center = std::isfinite(anchor.m_mean)
        ? anchor.m_mean : anchor.m_energy;
      return std::any_of( std::begin(candidate_lines), std::end(candidate_lines),
        [anchor_center, &anchor, &fwhm_at_energy]( const PeakDef &candidate_line ) {
          const double candidate_fwhm = fwhm_at_energy( candidate_line.mean() );
          if( !std::isfinite(candidate_fwhm) || !(candidate_fwhm > 0.0) )
            return false;
          const double anchor_fwhm = (std::isfinite(anchor.m_fwhm) && (anchor.m_fwhm > 0.0))
            ? anchor.m_fwhm : candidate_fwhm;
          return std::fabs(candidate_line.mean() - anchor_center)
            < (0.5 * (anchor_fwhm + candidate_fwhm));
        } );
    } ) );
}


bool significant_requested_source_anchors_preserved(
  const RelActCalcAuto::RelActAutoSolution &candidate,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::vector<RelActCalcManual::GenericPeakInfo> &anchors,
  const double minimum_z,
  const std::function<double(double)> &fwhm_at_energy )
{
  return significant_requested_source_anchor_count(
      candidate, sources, anchors, minimum_z, fwhm_at_energy ) == anchors.size();
}


/** Average chi2 per CHANNEL (not per fitted dof) over the solution's ROIs that contain
 significant peaks, filling `insignificant_roi_indices` with the ROIs that do not.

 Everything is evaluated in one calibration frame - the SPECTRUM cal of the foreground the
 solution was fit on: `m_final_roi_ranges_in_spectrum_cal` ROI bounds, `m_peaks_without_back_sub`
 peaks, and `solution.m_foreground` data.  (Formerly this paired TRUE-energy `m_final_roi_ranges`
 with spectrum-cal peaks against a caller-supplied - possibly cal-advanced - spectrum, so when the
 solve fit a non-trivial energy-cal adjustment, either the channel windows or the peaks were
 displaced by the shift, mis-scoring ROI significance on NaI/CZT.)
 */
double compute_filtered_chi2_per_channel(
  const RelActCalcAuto::RelActAutoSolution &solution,
  const double min_roi_significance_z,
  std::vector<size_t> &insignificant_roi_indices )
{
  insignificant_roi_indices.clear();

  double total_chi2 = 0.0;
  size_t total_channels = 0;

  const std::shared_ptr<const SpecUtils::Measurement> &data = solution.m_foreground;
  if( !data )
    return std::numeric_limits<double>::max();

  // Spectrum-cal ROI bounds match the spectrum-cal peaks and data; fall back to the true-energy
  // ranges only if the spectrum-cal vector was not populated (e.g., failed solve).
  const bool have_spec_cal_rois
    = (solution.m_final_roi_ranges_in_spectrum_cal.size() == solution.m_final_roi_ranges.size());
  assert( have_spec_cal_rois || solution.m_final_roi_ranges_in_spectrum_cal.empty() );

  const std::vector<RelActCalcAuto::RoiRange> &roi_ranges = have_spec_cal_rois
    ? solution.m_final_roi_ranges_in_spectrum_cal
    : solution.m_final_roi_ranges;

  for( size_t roi_idx = 0; roi_idx < roi_ranges.size(); ++roi_idx )
  {
    const RelActCalcAuto::RoiRange &roi = roi_ranges[roi_idx];

    const RoiSignificanceResult sig_result = compute_roi_chi2_significance(
      roi, solution.m_peaks_without_back_sub, data, min_roi_significance_z );

    if( sig_result.has_significant_peaks )
    {
      total_chi2 += sig_result.chi2_with_peaks;
      total_channels += sig_result.num_channels;
    }
    else
    {
      insignificant_roi_indices.push_back( roi_idx );
    }
  }//for( loop over ROIs )

  if( total_channels == 0 )
    return std::numeric_limits<double>::max();

  return total_chi2 / static_cast<double>( total_channels );
}//compute_filtered_chi2_per_channel


bool should_combine_peaks( const PeakDef &larger_peak,
                           const PeakDef &smaller_peak,
                           const double always_combine_nsigma )
{
  assert( larger_peak.peakArea() >= smaller_peak.peakArea() );

  const double dist = std::fabs( larger_peak.mean() - smaller_peak.mean() );
  const double sigma_large = larger_peak.sigma();

  // Criterion 1: Always combine if within threshold sigma
  if( dist < always_combine_nsigma * sigma_large )
    return true;

  // Criterion 2: Check if smaller peak is overwhelmed by larger peak's contribution
  // Calculate larger peak's contribution at smaller peak's location
  const double sigma_small = smaller_peak.sigma();
  const double mean_small = smaller_peak.mean();

  // Integrate larger peak over +/- 0.5 sigma of smaller peak's mean
  const double x0 = mean_small - 0.5 * sigma_small;
  const double x1 = mean_small + 0.5 * sigma_small;
  const double contribution = larger_peak.gauss_integral( x0, x1 );

  // If smaller peak area < 4 * sqrt(contribution), combine
  // This means the smaller peak is not statistically distinguishable from the tail of the larger peak
  const double smaller_area = smaller_peak.peakArea();
  if( contribution > 0.0 && smaller_area < 4.0 * std::sqrt( contribution ) )
    return true;

  return false;
}//should_combine_peaks


PeakDef combine_peaks( const std::vector<const PeakDef *> &peaks_to_combine )
{
  assert( !peaks_to_combine.empty() );
  if( peaks_to_combine.empty() )
    throw std::invalid_argument( "combine_peaks: empty peaks vector" );

  // Find dominant peak (largest area) and verify all share same continuum
  const PeakDef *dominant = peaks_to_combine[0];
  const std::shared_ptr<const PeakContinuum> cont = peaks_to_combine[0]->continuum();

  double total_area = 0.0;
  double sum_area_mean = 0.0;
  double sum_area_sigma = 0.0;
  double sum_uncert_sq = 0.0;
  double sum_means = 0.0, sum_sigmas = 0.0;

  for( const PeakDef *peak : peaks_to_combine )
  {
    assert( peak->continuum() == cont );
    if( peak->continuum() != cont )
      throw std::invalid_argument( "combine_peaks: peaks must share same continuum" );

    const double area = peak->peakArea();
    const double uncert = peak->peakAreaUncert();
    // Guard the -1 uncertainty sentinel before squaring (else (-1)*(-1)=+1 fabricates a bogus
    // uncertainty on the combined peak); fall back to a Poisson sqrt(area).
    const double eff_uncert = (uncert > 0.0) ? uncert : ((area > 0.0) ? std::sqrt(area) : 0.0);

    total_area += area;
    sum_area_mean += area * peak->mean();
    sum_area_sigma += area * peak->sigma();
    sum_uncert_sq += eff_uncert * eff_uncert;
    sum_means += peak->mean();
    sum_sigmas += peak->sigma();

    if( area > dominant->peakArea() )
      dominant = peak;
  }//for( const PeakDef *peak : peaks_to_combine )

  // Start with a copy of the dominant peak - this copies the continuum, source assignment,
  // skew type, line color, and other settings
  PeakDef combined = *dominant;

  // Update the gaussian parameters to the combined values
  combined.setPeakArea( total_area );
  combined.setPeakAreaUncert( std::sqrt( sum_uncert_sq ) );

  if( total_area > 0.0 )
  {
    combined.setMean( sum_area_mean / total_area );
    combined.setSigma( sum_area_sigma / total_area );
  }else
  {
    combined.setMean( sum_means / peaks_to_combine.size() );
    combined.setSigma( sum_sigmas / peaks_to_combine.size() );
  }


  return combined;
}//combine_peaks


std::vector<PeakDef> combine_overlapping_peaks_in_rois(
    const std::vector<PeakDef> &uncombined_peaks,
    const std::function<bool(const PeakDef &,const PeakDef &)> &may_combine = {} )
{
  if( uncombined_peaks.empty() )
    return {};

  // Phase 1: Group peaks by continuum (ROI), sorted by energy for deterministic order.
  std::vector<std::pair<const PeakContinuum *, std::vector<PeakDef>>> sorted_groups
    = group_peaks_by_roi( uncombined_peaks );

  std::vector<PeakDef> combined_peaks;
  combined_peaks.reserve( uncombined_peaks.size() );

  // Phase 2 & 3: Process each ROI independently
  for( const std::pair<const PeakContinuum *, std::vector<PeakDef>> &roi_group : sorted_groups )
  {
    const std::vector<PeakDef> &roi_peaks = roi_group.second;

    if( roi_peaks.size() == 1 )
    {
      // Single peak in ROI - no combination needed
      combined_peaks.push_back( roi_peaks[0] );
      continue;
    }

    // Sort by peak area (descending) - we'll process largest peaks first
    std::vector<size_t> sorted_indices( roi_peaks.size() );
    std::iota( sorted_indices.begin(), sorted_indices.end(), 0 );
    std::sort( sorted_indices.begin(), sorted_indices.end(),
      [&roi_peaks]( size_t a, size_t b ) {
        return roi_peaks[a].peakArea() > roi_peaks[b].peakArea();
      });

    // Greedy clustering: start with largest peak, cluster nearby peaks into it,
    // then move to next largest unclustered peak
    std::set<size_t> clustered;  // Tracks which indices have been assigned to a cluster
    std::vector<std::vector<const PeakDef *>> clusters;

    for( size_t i = 0; i < sorted_indices.size(); ++i )
    {
      const size_t idx_i = sorted_indices[i];

      // Skip if this peak has already been clustered with a larger peak
      if( clustered.count( idx_i ) )
        continue;

      // Start a new cluster with this peak (the largest remaining unclustered peak)
      const PeakDef &anchor_peak = roi_peaks[idx_i];
      std::vector<const PeakDef *> cluster;
      cluster.push_back( &anchor_peak );
      clustered.insert( idx_i );

      // Find all smaller peaks that should be combined with this anchor peak
      for( size_t j = i + 1; j < sorted_indices.size(); ++j )
      {
        const size_t idx_j = sorted_indices[j];

        // Skip if already in another cluster
        if( clustered.count( idx_j ) )
          continue;

        const PeakDef &candidate_peak = roi_peaks[idx_j];

        // Check if this smaller peak should be combined with the anchor peak
        if( (!may_combine || may_combine(anchor_peak, candidate_peak))
            && should_combine_peaks( anchor_peak, candidate_peak, 1.5 ) )
        {
          cluster.push_back( &candidate_peak );
          clustered.insert( idx_j );
        }
      }//for( size_t j = i + 1; j < sorted_indices.size(); ++j )

      clusters.push_back( std::move( cluster ) );
    }//for( size_t i = 0; i < sorted_indices.size(); ++i )

    // Create combined peaks for each cluster
    std::vector<PeakDef> combined_roi_peaks;
    for( const std::vector<const PeakDef *> &cluster : clusters )
    {
      if( cluster.size() == 1 )
      {
        combined_roi_peaks.push_back( *cluster[0] );
      }
      else
      {
        combined_roi_peaks.push_back( combine_peaks( cluster ) );
      }
    }//for( const std::vector<const PeakDef *> &cluster : clusters )

    // Make all peaks in this ROI share a new continuum
    if( !combined_roi_peaks.empty() )
    {
      std::shared_ptr<PeakContinuum> roi_continuum = std::make_shared<PeakContinuum>( *combined_roi_peaks.front().continuum() );
      for( PeakDef &peak : combined_roi_peaks )
        peak.setContinuum( roi_continuum );
    }

    // Add this ROI's peaks to the combined output
    combined_peaks.insert( combined_peaks.end(), combined_roi_peaks.begin(), combined_roi_peaks.end() );
  }//for( roi_group : sorted_groups )

  // Sort combined peaks by mean energy
  std::sort( combined_peaks.begin(), combined_peaks.end(),
    []( const PeakDef &lhs, const PeakDef &rhs ) {
      return lhs.mean() < rhs.mean();
    });

  return combined_peaks;
}//combine_overlapping_peaks_in_rois


std::vector<PeakDef> compute_observable_peaks(
  const std::vector<PeakDef> &fit_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const PeakFitUtils::CoarseResolutionType det_type,
  const PeakFitForNuclideConfig &config
#if( OBSERVABLE_PEAKS_USING_ORIGINAL_CAL_WITH_BACK_SUB )
  , const std::shared_ptr<const SpecUtils::Measurement> &background
#endif
  , const std::function<bool(const PeakDef &,const PeakDef &)> &may_combine_peaks = {}
  , const std::function<bool(const PeakDef &)> &must_refit_peak = {}
  , const std::function<bool(const PeakDef &)> &must_keep_peak = {}
  )
{
  
#if( PERFORM_DEVELOPER_CHECKS )
  // Debug: Check if input fit_peaks have overlapping ROIs
  if( should_debug_print() )
  {
    std::map<std::shared_ptr<const PeakContinuum>, std::vector<double>> continuum_to_peak_means_input;
    for( const PeakDef &peak : fit_peaks )
      continuum_to_peak_means_input[peak.continuum()].push_back( peak.mean() );

    struct RoiDebugInfoInput
    {
      double lower_energy;
      double upper_energy;
      std::vector<double> peak_means;
    };

    std::vector<RoiDebugInfoInput> rois_debug_input;
    for( const auto &entry : continuum_to_peak_means_input )
    {
      RoiDebugInfoInput info;
      info.lower_energy = entry.first->lowerEnergy();
      info.upper_energy = entry.first->upperEnergy();
      info.peak_means = entry.second;
      rois_debug_input.push_back( info );
    }

    std::sort( std::begin(rois_debug_input), std::end(rois_debug_input),
      []( const RoiDebugInfoInput &a, const RoiDebugInfoInput &b ) { return a.lower_energy < b.lower_energy; } );

    bool found_overlap_input = false;
    for( size_t i = 1; i < rois_debug_input.size(); ++i )
    {
      const RoiDebugInfoInput &prev_roi = rois_debug_input[i - 1];
      const RoiDebugInfoInput &curr_roi = rois_debug_input[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap_input )
          std::cerr << "compute_observable_peaks: INPUT fit_peaks ALREADY HAVE OVERLAPPING ROIs:" << std::endl;
        found_overlap_input = true;

        std::cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : prev_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;

        std::cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : curr_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;
        std::cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << std::endl;
      }
    }

    if( !found_overlap_input && should_debug_print() )
      std::cout << "compute_observable_peaks: Input fit_peaks have no overlapping ROIs" << std::endl;
  }//if( should_debug_print() )
#endif

  std::shared_ptr<const SpecUtils::Measurement> refit_data = foreground;
#if( OBSERVABLE_PEAKS_USING_ORIGINAL_CAL_WITH_BACK_SUB )
  // If background is provided, create a background-subtracted spectrum for refitting.
  // The bg-subtracted data gives Gaussian params that match the peak signal,
  // while the raw foreground is used afterwards to refit just the continuum for display.
  if( background && background->energy_calibration() && background->energy_calibration()->valid()
     && (background->live_time() > 0.0f) && (foreground->live_time() > 0.0f) )
  {
    const double lt_sf = foreground->live_time() / background->live_time();
    std::vector<float> channel_counts = *foreground->gamma_counts();
    const std::vector<float> &fore_energies = *foreground->energy_calibration()->channel_energies();
    const std::vector<float> &back_energies = *background->energy_calibration()->channel_energies();

    std::vector<float> bg_rebinned;
    SpecUtils::rebin_by_lower_edge( back_energies, *background->gamma_counts(),
                                     fore_energies, bg_rebinned );

    for( size_t i = 0; i < channel_counts.size(); ++i )
    {
      const double fg = std::max( static_cast<double>( channel_counts[i] ), 0.0 );
      const double bg = (i < bg_rebinned.size()) ? std::max( static_cast<double>( bg_rebinned[i] ), 0.0 ) : 0.0;
      channel_counts[i] = static_cast<float>( std::max( fg - lt_sf * bg, 0.0 ) );
    }

    std::shared_ptr<SpecUtils::Measurement> bg_sub = std::make_shared<SpecUtils::Measurement>( *foreground );
    bg_sub->set_gamma_counts( std::make_shared<std::vector<float>>( std::move(channel_counts) ),
                               foreground->live_time(), foreground->real_time() );
    refit_data = bg_sub;
  }//if( background provided )
#endif

  // Fraction of Gaussian area within +/-1 FWHM
  // erf(sqrt(ln(2))) = 0.7607
  const double fwhm_fraction = 0.7607;
  const double initial_significance_threshold = config.observable_peak_initial_significance_threshold;
  const double final_significance_threshold = config.observable_peak_final_significance_threshold;

  // ROI half-widths used when shrinking a ROI after edge peaks are dropped: the configured core
  // plus a fixed one-FWHM sideband (the dropped-peak region was continuum-consistent, but there
  // is no per-ROI adaptive information at this point).
  const double skew_low_extra = (config.skew_type == PeakDef::SkewType::NoSkew)
      ? 0.0 : sm_skew_low_side_extra_fwhm;
  const double roi_width_num_fwhm_lower = config.auto_roi_core_num_fwhm + skew_low_extra + sm_post_drop_sideband_fwhm;
  const double roi_width_num_fwhm_upper = config.auto_roi_core_num_fwhm + sm_post_drop_sideband_fwhm;

  // Lambda to adjust ROI bounds when edge peaks are removed.
  // Returns true if bounds were adjusted, false otherwise.
  // Takes peaks sorted by mean energy.
  const auto reduce_roi_bounds_if_needed = [roi_width_num_fwhm_lower, roi_width_num_fwhm_upper](
    std::vector<PeakDef> &peaks,
    const double orig_left_mean,
    const double orig_right_mean ) -> bool
  {
    if( peaks.empty() )
      return false;

    // Sort peaks by mean for edge detection
    std::sort( std::begin(peaks), std::end(peaks), &PeakDef::lessThanByMean );

    const double new_left_mean = peaks.front().mean();
    const double new_right_mean = peaks.back().mean();
    const bool left_edge_changed = (std::fabs(new_left_mean - orig_left_mean) > 0.1);
    const bool right_edge_changed = (std::fabs(new_right_mean - orig_right_mean) > 0.1);

    if( !left_edge_changed && !right_edge_changed )
      return false;

    // Get current continuum from first peak
    std::shared_ptr<const PeakContinuum> old_continuum = peaks.front().continuum();

    // Calculate new ROI bounds based on new edge peaks
    const double new_left_fwhm = peaks.front().fwhm();
    const double new_right_fwhm = peaks.back().fwhm();
    double new_lower_energy = new_left_mean - roi_width_num_fwhm_lower * new_left_fwhm;
    double new_upper_energy = new_right_mean + roi_width_num_fwhm_upper * new_right_fwhm;

    // Constrain new bounds to not expand beyond original ROI bounds
    // This function is only called when edge peaks are removed, so ROI should only shrink, never expand
    const double orig_lower = old_continuum->lowerEnergy();
    const double orig_upper = old_continuum->upperEnergy();
    new_lower_energy = std::max( new_lower_energy, orig_lower );
    new_upper_energy = std::min( new_upper_energy, orig_upper );

    // Create new continuum as copy with updated energy range
    std::shared_ptr<PeakContinuum> new_continuum = std::make_shared<PeakContinuum>( *old_continuum );
    new_continuum->setRange( new_lower_energy, new_upper_energy );

    // Set new continuum to all peaks
    for( PeakDef &p : peaks )
      p.setContinuum( new_continuum );

    if( should_debug_print() )
    {
      std::cout << "  Observable filter: adjusted ROI bounds from ["
           << old_continuum->lowerEnergy() << ", " << old_continuum->upperEnergy()
           << "] to [" << new_lower_energy << ", " << new_upper_energy << "] keV" << std::endl;
    }

    return true;
  };//reduce_roi_bounds_if_needed lambda

  // Step 1: Group input peaks by ROI (shared continuum), skipping peaks whose
  //  mean is outside their own continuum (tail contributions to adjacent ROIs).
  //  NOTE/TODO: if any peaks migrate out the ROI while fitting the observable peaks,
  //             we currently expand their ROI so this next filter wont delete them.
  //             This isnt the best solution, and should be improved.      
  std::vector<PeakDef> peaks_in_bounds;
  for( const PeakDef &peak : fit_peaks )
  {
    const double mean = peak.mean();
    if( peak.continuum()
      && (mean >= peak.continuum()->lowerEnergy())
      && (mean <= peak.continuum()->upperEnergy()) )
    {
      peaks_in_bounds.push_back( peak );
    }
  }

  std::vector<std::pair<const PeakContinuum *, std::vector<PeakDef>>> input_rois
    = group_peaks_by_roi( peaks_in_bounds );

  // Step 2: Initial significance filter and ROI adjustment per ROI
  std::vector<PeakDef> filtered_peaks;
  for( auto &roi_entry : input_rois )
  {
    std::vector<PeakDef> &roi_peaks = roi_entry.second;

    // Sort by mean for edge tracking
    std::sort( std::begin(roi_peaks), std::end(roi_peaks), &PeakDef::lessThanByMean );

    const double orig_left_mean = roi_peaks.front().mean();
    const double orig_right_mean = roi_peaks.back().mean();

    // Step-type continuum integrals need the ROI's peaks
    std::vector<std::shared_ptr<const PeakDef>> roi_peak_ptrs;
    roi_peak_ptrs.reserve( roi_peaks.size() );
    for( const PeakDef &p : roi_peaks )
      roi_peak_ptrs.push_back( std::make_shared<PeakDef>( p ) );

    // Filter peaks by initial significance
    std::vector<PeakDef> kept_peaks;
    for( const PeakDef &peak : roi_peaks )
    {
      const double mean = peak.mean();
      const double fwhm = peak.fwhm();
      const double lower_energy = mean - fwhm;
      const double upper_energy = mean + fwhm;

      // Initial significance z = S/sqrt(S+B) over +/-1 FWHM, with B from the FITTED continuum
      // rather than the gross data, so a neighboring peak's counts in the window no longer
      // dilute this peak's significance (and the test is invariant to live-time).
      const double cont_b = peak.continuum()
          ? std::max( 0.0, peak.continuum()->offset_integral( lower_energy, upper_energy,
                                                              foreground, roi_peak_ptrs ) )
          : foreground->gamma_integral( lower_energy, upper_energy );

      const double peak_contrib = peak.amplitude() * fwhm_fraction;
      const double significance = peak_contrib / std::sqrt( std::max( 1.0, peak_contrib + cont_b ) );

      if( (must_refit_peak && must_refit_peak(peak))
          || (significance >= initial_significance_threshold) )
      {
        kept_peaks.push_back( peak );
      }
      else if( should_debug_print() )
      {
        std::cout << "  Observable filter (initial sig=" << significance << " < " << initial_significance_threshold
             << "): peak at " << mean << " keV" << std::endl;
      }
    }//for( const PeakDef &peak : roi_peaks )

    // Adjust ROI bounds if edge peaks were removed
    if( !kept_peaks.empty() )
    {
      reduce_roi_bounds_if_needed( kept_peaks, orig_left_mean, orig_right_mean );
      filtered_peaks.insert( std::end(filtered_peaks), std::begin(kept_peaks), std::end(kept_peaks) );
    }
  }//for( auto &roi_entry : input_rois )

  if( filtered_peaks.empty() )
    return filtered_peaks;

  // Step 3: Group peaks by shared continuum (ROI) - may have new continuums from adjustment.
  std::vector<std::pair<const PeakContinuum *, std::vector<PeakDef>>> rois
    = group_peaks_by_roi( filtered_peaks );

  // Debug: Check if rois already have overlaps before parallel processing
#if( PERFORM_DEVELOPER_CHECKS )
  if( should_debug_print() )
  {
    struct RoiDebugInfo2
    {
      double lower_energy;
      double upper_energy;
      size_t num_peaks;
    };

    std::vector<RoiDebugInfo2> rois_debug2;
    for( const auto &entry : rois )
    {
      RoiDebugInfo2 info;
      info.lower_energy = entry.first->lowerEnergy();
      info.upper_energy = entry.first->upperEnergy();
      info.num_peaks = entry.second.size();
      rois_debug2.push_back( info );
    }

    std::sort( std::begin(rois_debug2), std::end(rois_debug2),
      []( const RoiDebugInfo2 &a, const RoiDebugInfo2 &b ) { return a.lower_energy < b.lower_energy; } );

    bool found_overlap2 = false;
    for( size_t i = 1; i < rois_debug2.size(); ++i )
    {
      const RoiDebugInfo2 &prev_roi = rois_debug2[i - 1];
      const RoiDebugInfo2 &curr_roi = rois_debug2[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap2 )
          std::cerr << "compute_observable_peaks: OVERLAPS BEFORE PARALLEL PROCESSING:" << std::endl;
        found_overlap2 = true;

        std::cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy
             << "] keV, " << prev_roi.num_peaks << " peaks" << std::endl;
        std::cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy
             << "] keV, " << curr_roi.num_peaks << " peaks" << std::endl;
        std::cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << std::endl;
      }
    }

    if( !found_overlap2 && should_debug_print() )
      std::cout << "compute_observable_peaks: No overlaps before parallel processing" << std::endl;
  }
#endif

  // Step 4: Iteratively refit each ROI and remove insignificant peaks - in parallel
  const size_t num_rois = rois.size();
  std::vector<std::vector<PeakDef>> roi_results( num_rois );

  // Copy ROI data to vector for parallel processing
  std::vector<std::pair<const PeakContinuum *, std::vector<PeakDef>>> roi_vec( rois.begin(), rois.end() );

  SpecUtilsAsync::ThreadPool pool;

  for( size_t roi_index = 0; roi_index < num_rois; ++roi_index )
  {
    pool.post( [roi_index, &roi_vec, &roi_results, &refit_data,
                final_significance_threshold, det_type, &reduce_roi_bounds_if_needed,
                &may_combine_peaks, &must_refit_peak, &must_keep_peak]()
    {
      std::vector<PeakDef> roi_peaks = roi_vec[roi_index].second;

      // Sort peaks by mean energy for edge detection
      std::sort( std::begin(roi_peaks), std::end(roi_peaks), &PeakDef::lessThanByMean );

      const size_t max_iterations = 3;
      size_t iteration = 0;
      bool changed = true;
      while( changed && !roi_peaks.empty() && (iteration < max_iterations) )
      {
        changed = false;

        // Track original edge peaks before filtering
        const double orig_left_mean = roi_peaks.front().mean();
        const double orig_right_mean = roi_peaks.back().mean();

        // Convert to shared_ptr for refitPeaksThatShareROI_LM
        std::vector<std::shared_ptr<const PeakDef>> input_peaks;
        for( const PeakDef &p : roi_peaks )
          input_peaks.push_back( std::make_shared<PeakDef>(p) );

#if( PERFORM_DEVELOPER_CHECKS )
        if( should_debug_print() )
        {
          std::cout << "compute_observable_peaks: refit input ROI ["
               << roi_peaks.front().continuum()->lowerEnergy() << ", "
               << roi_peaks.front().continuum()->upperEnergy() << "] keV, "
               << input_peaks.size() << " peaks:" << std::endl;
          for( const auto &p : input_peaks )
          {
            std::cout << "  mean=" << p->mean() << " keV, sigma=" << p->sigma()
                 << ", fwhm=" << p->fwhm() << ", amp=" << p->amplitude()
                 << ", area=" << p->peakArea() << ", areaUncert=" << p->peakAreaUncert()
                 << ", chi2dof=" << p->chi2dof()
                 << (p->parentNuclide() ? (std::string(", nuc=") + p->parentNuclide()->symbol) : "")
                 << std::endl;
          }
        }
#endif

        Wt::WFlags<PeakFitLM::PeakFitLMOptions> refine_amount;
        if( det_type == PeakFitUtils::CoarseResolutionType::High )
          refine_amount |= PeakFitLM::PeakFitLMOptions::SmallRefinementOnly;
        else
          refine_amount |= PeakFitLM::PeakFitLMOptions::MediumRefinementOnly;

        // Use fit_peaks_in_spectrum_LM rather than refitPeaksThatShareROI_LM: when a peak becomes
        // INSIGNIFICANT in the (honest) refit it is reported in `lost_peaks` and DROPPED here.
        // refitPeaksThatShareROI_LM instead returns an empty vector in that case, which is
        // indistinguishable from a numerical failure, so the caller used to KEEP the original peak
        // (often a continuum-curvature / degenerate "peak" with inflated significance).  Now only a
        // genuine numerical failure (status != Success) falls back to keeping the current peaks.
        const bool has_protected_peak = must_refit_peak && std::any_of(
          std::begin(roi_peaks), std::end(roi_peaks), must_refit_peak );
        const double refit_significance_threshold
          = has_protected_peak ? 0.0 : final_significance_threshold;
        const PeakFitLM::FitPeaksResults refit_res = PeakFitLM::fit_peaks_in_spectrum_LM(
            input_peaks, refit_data, /*stat_threshold=*/ refit_significance_threshold,
            /*hypothesis_threshold=*/ 0.0, det_type, std::nullopt /*keep peaks' own skew*/,
            refine_amount, may_combine_peaks );

        if( refit_res.status != PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Success )
          break;  // genuine refit failure - keep current peaks (defensive)

        // Survivors (refit_res.fit_peaks); peaks in refit_res.lost_peaks became insignificant and are
        // dropped by simply not carrying them forward.
        const std::vector<std::shared_ptr<const PeakDef>> &refit_result = refit_res.fit_peaks;

        // Dropping insignificant peaks changes the ROI, so iterate to refit the survivors.
        if( !refit_res.lost_peaks.empty() )
          changed = true;

#if( PERFORM_DEVELOPER_CHECKS )
        if( should_debug_print() )
        {
          std::cout << "  Refit produced " << refit_result.size() << " peaks ("
               << refit_res.lost_peaks.size() << " dropped as insignificant):" << std::endl;
          for( const auto &p : refit_result )
            std::cout << "    mean=" << p->mean() << " keV, area=" << p->peakArea()
                 << ", areaUncert=" << p->peakAreaUncert() << std::endl;
          for( const auto &p : refit_res.lost_peaks )
            std::cout << "    DROPPED(insignificant): mean=" << p->mean() << " keV, area=" << p->peakArea() << std::endl;
        }
#endif

        // Check significance and remove insignificant peaks
        vector<PeakDef> kept_peaks;
        for( const shared_ptr<const PeakDef> &peak : refit_result )
        {
          const double mean = peak->mean();
          const double area = peak->peakArea();
          const double area_uncert = peak->peakAreaUncert();
          // Guard the -1 uncertainty sentinel: fall back to Poisson sqrt(area).
          const double final_sig = (area_uncert > 0.0)
            ? (area / area_uncert)
            : ((area > 0.0) ? std::sqrt(area) : 0.0);

          if( (must_keep_peak && must_keep_peak(*peak))
              || (final_sig >= final_significance_threshold) )
          {
            kept_peaks.push_back( *peak );

#if( PERFORM_DEVELOPER_CHECKS )
            if( should_debug_print() )
            {
              // Compute chi2 with peak vs continuum-only
              const double roi_lower = peak->continuum()->lowerEnergy();
              const double roi_upper = peak->continuum()->upperEnergy();
              const int lch = static_cast<int>( refit_data->find_gamma_channel( roi_lower ) );
              const int uch = static_cast<int>( refit_data->find_gamma_channel( roi_upper ) );
              const int nch = std::max( uch - lch, 1 );

              std::vector<std::shared_ptr<const PeakDef>> peak_vec = { peak };
              const double chi2_with_peak = chi2_for_region( peak_vec, refit_data, lch, uch );

              // Continuum-only chi2: make a copy with zero amplitude
              PeakDef no_peak_copy = *peak;
              no_peak_copy.setAmplitude( 0.0 );
              std::vector<std::shared_ptr<const PeakDef>> no_peak_vec
                = { std::make_shared<PeakDef>( no_peak_copy ) };
              const double chi2_cont_only = chi2_for_region( no_peak_vec, refit_data, lch, uch );

              std::cout << "  Kept peak at " << mean << " keV: sig=" << final_sig
                   << ", chi2/dof=" << (chi2_with_peak / nch)
                   << ", chi2_cont_only/dof=" << (chi2_cont_only / nch)
                   << ", chi2_improvement=" << (chi2_cont_only - chi2_with_peak)
                   << ", nch=" << nch << std::endl;
            }
#endif
          }
          else
          {
            changed = true;
            if( should_debug_print() )
            {
              std::cout << "  Observable filter post-refit (final sig=" << final_sig << " < " << final_significance_threshold
                   << "): peak at " << mean << " keV" << std::endl;
            }
          }
        }//for( const std::shared_ptr<const PeakDef> &peak : refit_result )

        // Adjust ROI bounds if edge peaks were removed
        if( !kept_peaks.empty() )
        {
          const bool bounds_adjusted = reduce_roi_bounds_if_needed( kept_peaks, orig_left_mean, orig_right_mean );
          if( bounds_adjusted )
            changed = true;  // Need to refit with new continuum bounds
        }

        roi_peaks = std::move( kept_peaks );
        iteration += 1;
      }//while( changed && !roi_peaks.empty() && (iteration < max_iterations) )

      // The refit (using MediumRefinementOnly) allows peak means to shift up to ±0.5*sigma,
      // which can push a peak mean slightly outside the continuum bounds when the peak
      // starts near the ROI edge.  Rather than constraining the fitter (which could hurt
      // fit quality for edge peaks), we expand the ROI bounds to encompass all peak means.
      // TODO: revisit whether the refit should more tightly constrain means to stay within
      //       the ROI, or whether the ROI bounds from the RelActAuto solver should provide
      //       more margin around edge peaks to begin with.
      // NOTE: the `reduce_roi_bounds_if_needed` lambda will delete any pak whose mean is outside the ROI
      if( !roi_peaks.empty() )
      {
        std::shared_ptr<const PeakContinuum> continuum = roi_peaks.front().continuum();
        assert( continuum );

        double lower = continuum->lowerEnergy();
        double upper = continuum->upperEnergy();
        bool needs_expansion = false;

        for( const PeakDef &p : roi_peaks )
        {
          if( p.mean() < lower )
          {
            lower = p.mean() - 0.1;
            needs_expansion = true;
          }
          if( p.mean() > upper )
          {
            upper = p.mean() + 0.1;
            needs_expansion = true;
          }
        }//for( const PeakDef &p : roi_peaks )

        if( needs_expansion )
        {
          std::shared_ptr<PeakContinuum> new_continuum = std::make_shared<PeakContinuum>( *continuum );
          new_continuum->setRange( lower, upper );
          for( PeakDef &p : roi_peaks )
            p.setContinuum( new_continuum );
        }
      }//if( !roi_peaks.empty() )

      roi_results[roi_index] = std::move( roi_peaks );
    });//pool.post lambda
  }//for( size_t roi_index = 0; roi_index < num_rois; ++roi_index )

  pool.join();

  // Collect all results
  std::vector<PeakDef> observable_peaks;
  for( const std::vector<PeakDef> &roi_result : roi_results )
    observable_peaks.insert( std::end(observable_peaks), std::begin(roi_result), std::end(roi_result) );

  // Sort by energy
  std::sort( std::begin(observable_peaks), std::end(observable_peaks), &PeakDef::lessThanByMean );

#if( OBSERVABLE_PEAKS_USING_ORIGINAL_CAL_WITH_BACK_SUB )
  // If we refit on bg-subtracted data, the peak Gaussians are correct but the
  // continuum won't match the raw foreground for display.  Refit just the continuum
  // on the raw foreground while keeping peak Gaussians fixed.
  if( background && !observable_peaks.empty() )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    if( should_debug_print() )
    {
      std::cout << "compute_observable_peaks: before refit_roi_continuums (" << observable_peaks.size() << " peaks):" << std::endl;
      for( const PeakDef &p : observable_peaks )
      {
        std::cout << "  mean=" << p.mean() << " keV, sigma=" << p.sigma()
             << ", area=" << p.peakArea() << ", chi2dof=" << p.chi2dof()
             << ", ROI=[" << p.continuum()->lowerEnergy() << ", " << p.continuum()->upperEnergy() << "]"
             << std::endl;
      }
    }
#endif

    observable_peaks = RelActCalc::refit_roi_continuums( observable_peaks, foreground );

#if( PERFORM_DEVELOPER_CHECKS )
    if( should_debug_print() )
    {
      std::cout << "compute_observable_peaks: after refit_roi_continuums (" << observable_peaks.size() << " peaks):" << std::endl;
      for( const PeakDef &p : observable_peaks )
      {
        std::cout << "  mean=" << p.mean() << " keV, sigma=" << p.sigma()
             << ", area=" << p.peakArea() << ", chi2dof=" << p.chi2dof()
             << ", ROI=[" << p.continuum()->lowerEnergy() << ", " << p.continuum()->upperEnergy() << "]"
             << std::endl;
      }
    }
#endif
  }
#endif
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  // Debug: Check for overlapping ROI bounds
  if( should_debug_print() )
  {
    // Group peaks by continuum to identify unique ROIs
    std::map<std::shared_ptr<const PeakContinuum>, std::vector<double>> continuum_to_peak_means;
    for( const PeakDef &peak : observable_peaks )
      continuum_to_peak_means[peak.continuum()].push_back( peak.mean() );

    // Extract ROIs sorted by lower energy
    struct RoiDebugInfo
    {
      double lower_energy;
      double upper_energy;
      std::vector<double> peak_means;
    };

    std::vector<RoiDebugInfo> rois_debug;
    for( const auto &entry : continuum_to_peak_means )
    {
      RoiDebugInfo info;
      info.lower_energy = entry.first->lowerEnergy();
      info.upper_energy = entry.first->upperEnergy();
      info.peak_means = entry.second;
      rois_debug.push_back( info );
    }

    std::sort( std::begin(rois_debug), std::end(rois_debug),
      []( const RoiDebugInfo &a, const RoiDebugInfo &b ) { return a.lower_energy < b.lower_energy; } );

    // Check for overlaps
    bool found_overlap = false;
    for( size_t i = 1; i < rois_debug.size(); ++i )
    {
      const RoiDebugInfo &prev_roi = rois_debug[i - 1];
      const RoiDebugInfo &curr_roi = rois_debug[i];

      if( curr_roi.lower_energy < prev_roi.upper_energy )
      {
        if( !found_overlap )
          std::cerr << "compute_observable_peaks: FOUND OVERLAPPING ROIs:" << std::endl;
        found_overlap = true;

        std::cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : prev_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;

        std::cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
        for( double mean : curr_roi.peak_means )
          std::cerr << mean << " ";
        std::cerr << "keV" << std::endl;
        std::cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << std::endl;
      }
    }

    if( !found_overlap && should_debug_print() )
      std::cout << "compute_observable_peaks: No overlapping ROIs detected" << std::endl;
  }//if( should_debug_print() )
#endif

  // Drop orphan 511 keV peaks: an unattributed peak near the annihilation
  // energy whose ROI contains no source-attributed peak is almost certainly
  // ambient annihilation in the background and shouldn't be claimed.
  // Sources whose 511 peak IS attributed (e.g. F-18) and ROIs that also
  // contain another source-attributed peak (e.g. NORM Bi-214 at 511) are
  // unaffected.
  {
    const double annihilation_energy = 510.9989;
    const double energy_tolerance = 2.0;

    auto has_source = []( const PeakDef &p ) -> bool {
      return p.hasSourceGammaAssigned() || p.parentNuclide()
             || p.xrayElement() || p.reaction();
    };

    // Identify each ROI's peaks via shared continuum pointer identity.
    std::map<std::shared_ptr<const PeakContinuum>, std::vector<size_t>> roi_indices;
    for( size_t i = 0; i < observable_peaks.size(); ++i )
      roi_indices[ observable_peaks[i].continuum() ].push_back( i );

    std::vector<bool> to_remove( observable_peaks.size(), false );
    for( const auto &entry : roi_indices )
    {
      const std::vector<size_t> &idxs = entry.second;

      bool any_source_in_roi = false;
      for( size_t i : idxs )
        any_source_in_roi = any_source_in_roi || has_source( observable_peaks[i] );

      if( any_source_in_roi )
        continue;  // ROI has at least one attributed peak; keep everything

      for( size_t i : idxs )
      {
        const PeakDef &p = observable_peaks[i];
        if( std::fabs( p.mean() - annihilation_energy ) < energy_tolerance
            && !has_source(p) )
        {
          to_remove[i] = true;
          if( should_debug_print() )
          {
            std::cout << "compute_observable_peaks: dropping orphan 511 peak at "
                 << p.mean() << " keV (no source in ROI ["
                 << (p.continuum() ? p.continuum()->lowerEnergy() : 0.0) << ", "
                 << (p.continuum() ? p.continuum()->upperEnergy() : 0.0) << "])"
                 << std::endl;
          }
        }
      }
    }

    if( std::any_of( begin(to_remove), end(to_remove), []( bool b ){ return b; } ) )
    {
      std::vector<PeakDef> kept;
      kept.reserve( observable_peaks.size() );
      for( size_t i = 0; i < observable_peaks.size(); ++i )
      {
        if( !to_remove[i] )
          kept.push_back( std::move(observable_peaks[i]) );
      }
      observable_peaks = std::move(kept);
    }
  }

  return observable_peaks;
}//compute_observable_peaks


std::vector<double> build_synthetic_spectrum(
  const std::vector<double> &gamma_energies,
  const std::vector<double> &gamma_amplitudes,
  const std::function<double(double)> &fwhm_at_energy,
  const float *channel_energies,
  size_t start_channel,
  size_t num_channels )
{
  assert( gamma_energies.size() == gamma_amplitudes.size() );

  // Zero-initialize the synthetic spectrum (same binning as data)
  std::vector<double> synthetic( num_channels, 0.0 );

  // Add each gamma line's Gaussian contribution using photopeak_function_integral
  for( size_t g = 0; g < gamma_energies.size(); ++g )
  {
    const double mean = gamma_energies[g];
    const double amplitude = gamma_amplitudes[g];
    const double sigma = fwhm_at_energy( mean ) / 2.35482;  // FWHM to sigma

    // photopeak_function_integral adds the peak contribution to synthetic[]
    // Pass pointer to channel energies starting at start_channel
    PeakDists::photopeak_function_integral<double>(
      mean,
      sigma,
      amplitude,
      PeakDef::SkewType::NoSkew,  // No skew for synthetic spectrum
      nullptr,                     // No skew parameters
      num_channels,
      &channel_energies[start_channel],
      synthetic.data()
    );
  }

  return synthetic;
}//build_synthetic_spectrum


double compute_significance_in_region(
  const std::vector<double> &synthetic,
  size_t start_channel,
  size_t check_start_ch,
  size_t check_end_ch,
  const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  if( check_end_ch <= check_start_ch )
    return 0.0;  // Too few channels, return neutral

  // Ensure check range is within synthetic range
  if( check_start_ch < start_channel )
    check_start_ch = start_channel;
  if( check_end_ch > start_channel + synthetic.size() )
    check_end_ch = start_channel + synthetic.size();

  if( check_end_ch <= check_start_ch )
    return 0.0;

  // Sum the expected peak counts from synthetic and the data counts
  double sum_synthetic = 0.0;
  double sum_data = 0.0;

  for( size_t ch = check_start_ch; ch < check_end_ch; ++ch )
  {
    const size_t syn_idx = ch - start_channel;
    sum_synthetic += synthetic[syn_idx];
    sum_data += data->gamma_channel_content( ch );
  }

  // Significance = expected_signal / noise
  // where noise = sqrt(data_counts) (Poisson statistics)
  if( sum_data <= 0.0 )
    return 0.0;

  return sum_synthetic / std::sqrt( sum_data );
}//compute_significance_in_region


bool has_significant_peak_between(
  size_t lower_channel,
  size_t upper_channel,
  const std::vector<double> &synthetic,
  size_t start_channel,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const float *channel_energies,
  const std::function<double(double)> &fwhm_at_energy,
  double check_fwhm_fraction,
  double significance_threshold )
{
  const size_t num_channels = synthetic.size();

  // Convert to indices within synthetic
  if( lower_channel < start_channel )
    lower_channel = start_channel;
  if( upper_channel > start_channel + num_channels )
    upper_channel = start_channel + num_channels;

  // Find local maxima in synthetic between lower_channel and upper_channel
  for( size_t ch = lower_channel; ch < upper_channel; ++ch )
  {
    const size_t i = ch - start_channel;
    if( i == 0 || i >= num_channels - 1 )
      continue;

    // Check if this is a local maximum in the synthetic spectrum
    if( synthetic[i] > synthetic[i-1] && synthetic[i] > synthetic[i+1] )
    {
      // Get energy at this channel for FWHM calculation
      const double ch_energy = 0.5 * (channel_energies[ch] + channel_energies[ch + 1]);
      const double fwhm = fwhm_at_energy( ch_energy );

      // Convert FWHM fraction to channel range for significance check
      const double half_width = check_fwhm_fraction * fwhm;
      const size_t check_start = data->find_gamma_channel( static_cast<float>(ch_energy - half_width) );
      const size_t check_end = data->find_gamma_channel( static_cast<float>(ch_energy + half_width) );

      const double significance = compute_significance_in_region(
        synthetic, start_channel, check_start, check_end, data );

      if( significance >= significance_threshold )
        return true;  // Found a significant peak
    }
  }

  return false;  // No significant peak found
}//has_significant_peak_between


std::vector<LocalMinimum> find_synthetic_minima(
  const std::vector<double> &synthetic,
  size_t start_channel,
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  const float *channel_energies,
  const std::function<double(double)> &fwhm_at_energy,
  double check_fwhm_fraction )
{
  std::vector<LocalMinimum> minima;
  const size_t num_channels = synthetic.size();

  for( size_t i = 1; i < num_channels - 1; ++i )
  {
    // Check if this is a local minimum
    if( synthetic[i] < synthetic[i-1] && synthetic[i] < synthetic[i+1] )
    {
      const double min_value = synthetic[i];
      const size_t abs_channel = start_channel + i;

      // Find left maximum (scan left until we find a local max or boundary)
      double left_max = min_value;
      for( size_t j = i; j > 0; --j )
      {
        if( synthetic[j] > left_max )
          left_max = synthetic[j];
        if( j > 0 && synthetic[j] > synthetic[j-1] && synthetic[j] > synthetic[j+1] )
          break;  // Found left local maximum
      }

      // Find right maximum
      double right_max = min_value;
      for( size_t j = i; j < num_channels - 1; ++j )
      {
        if( synthetic[j] > right_max )
          right_max = synthetic[j];
        if( synthetic[j] > synthetic[j-1] && synthetic[j] > synthetic[j+1] )
          break;  // Found right local maximum
      }

      // Compute relative depth score (for tiebreaking)
      const double smaller_max = std::min( left_max, right_max );
      const double larger_max = std::max( left_max, right_max );
      const double depth = smaller_max - min_value;
      const double depth_score = (larger_max > 0) ? depth / larger_max : 0.0;

      // Get energy at this channel for FWHM-based check region
      const double ch_energy = 0.5 * (channel_energies[abs_channel] + channel_energies[abs_channel + 1]);
      const double fwhm = fwhm_at_energy( ch_energy );
      const double half_width = check_fwhm_fraction * fwhm;
      const size_t check_start = data->find_gamma_channel( static_cast<float>(ch_energy - half_width) );
      const size_t check_end = data->find_gamma_channel( static_cast<float>(ch_energy + half_width) );

      const double significance = compute_significance_in_region(
        synthetic, start_channel, check_start, check_end, data );

      LocalMinimum lm;
      lm.channel = abs_channel;
      lm.synthetic_value = min_value;
      lm.depth_score = depth_score;
      lm.statistical_significance = significance;
      minima.push_back( lm );
    }
  }

  return minima;
}//find_synthetic_minima


/** Decide between a polynomial continuum and its step-continuum variant for a ROI by directly
 trial-fitting BOTH against the data and comparing chi2 - "let the model choose."

 The step candidate is chosen with the SAME parameter count as the polynomial (Linear vs FlatStep,
 Quadratic vs LinearStep), so the AICc complexity penalties cancel and the comparison reduces to
 chi2, biased against the step by `chi2_margin` (GA-tunable per detector class).

 Both trials are linear least-squares (`fit_amp_and_offset` - amplitudes + continuum coefficients
 for fixed means/sigmas), NOT Ceres solves: the cost is a couple of small matrix solves, which is
 noise next to a RelActAuto iteration even inside the GA.  Callers gate this behind (a) the
 peak-dominance significance test and (b) a loose sideband-asymmetry pre-filter, so only a handful
 of ROIs per spectrum ever reach it.

 The trial peak list is the cluster's predicted gammas, merged within 1 sigma (a weighted mean)
 and capped at the several largest, so line-dense sources (Eu152/Pu/U swarms of insignificant
 lines) do not produce an ill-conditioned solve - per-line amplitudes are free parameters here,
 not predictions, so tiny lines add columns without adding information.

 Fits use NoSkew: the fitted skew parameters are not known at clustering time, and using the same
 shape for both candidates keeps any skew-tail mismatch from favoring either one.

 Returns `poly_type` unchanged when the trial cannot be run (too few channels, degenerate fit).

 This replaces `should_use_step_continuum`'s fixed +/-1.5-FWHM probe windows, which structurally
 self-vetoed (the probes needed 1.625 FWHM inside a 1.5-FWHM core, so a ROI with no sideband
 extension - e.g. a peak against a Compton edge, the canonical step case - could never get a
 step) and read neighboring cluster gammas as continuum.
 */
PeakContinuum::OffsetType trial_step_continuum(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const double roi_lower,
  const double roi_upper,
  const std::vector<double> &gamma_energies,
  const std::vector<double> &gamma_amplitudes,
  const std::function<double(double)> &fwhm_at,
  const PeakContinuum::OffsetType poly_type,
  const double chi2_margin )
{
  assert( gamma_energies.size() == gamma_amplitudes.size() );

  if( !foreground || !foreground->channel_energies() || gamma_energies.empty()
      || !(roi_upper > roi_lower) )
    return poly_type;

  const PeakContinuum::OffsetType step_type = (poly_type == PeakContinuum::OffsetType::Quadratic)
      ? PeakContinuum::OffsetType::LinearStep
      : PeakContinuum::OffsetType::FlatStep;

  // Merge the cluster's predicted lines into effective trial peaks: lines within 1 sigma of a
  // (stronger) anchor combine into its amplitude-weighted mean; keep at most the 8 largest.
  std::vector<std::pair<double,double>> lines;  // (amplitude, energy), for sorting
  lines.reserve( gamma_energies.size() );
  for( size_t i = 0; i < gamma_energies.size(); ++i )
  {
    if( gamma_amplitudes[i] > 0.0 )
      lines.emplace_back( gamma_amplitudes[i], gamma_energies[i] );
  }
  std::sort( std::begin(lines), std::end(lines), std::greater<std::pair<double,double>>() );

  std::vector<double> means, sigmas, weights;
  for( const std::pair<double,double> &line : lines )
  {
    const double gamma_fwhm = fwhm_at( line.second );
    if( !std::isfinite(gamma_fwhm) || (gamma_fwhm <= 0.0) )
      continue;
    const double gamma_sigma = gamma_fwhm / PhysicalUnits::fwhm_nsigma;

    bool absorbed = false;
    for( size_t j = 0; j < means.size(); ++j )
    {
      if( std::fabs( line.second - means[j] ) < sigmas[j] )
      {
        means[j] = (means[j]*weights[j] + line.second*line.first) / (weights[j] + line.first);
        weights[j] += line.first;
        absorbed = true;
        break;
      }
    }

    if( !absorbed && (means.size() < 8) )
    {
      means.push_back( line.second );
      sigmas.push_back( gamma_sigma );
      weights.push_back( line.first );
    }
  }//for( const auto &line : lines )

  if( means.empty() )
    return poly_type;

  // Gather the ROI's channel data.
  const size_t nchannel = foreground->num_gamma_channels();
  const size_t start_ch = foreground->find_gamma_channel( static_cast<float>(roi_lower) );
  const size_t end_ch = std::min( foreground->find_gamma_channel( static_cast<float>(roi_upper) ),
                                  nchannel - 1 );
  if( (end_ch <= start_ch) || ((end_ch - start_ch) < (means.size() + 4)) )
    return poly_type;  // too few channels to distinguish the candidates

  const size_t nbin = end_ch - start_ch;
  const std::vector<float> &channel_energies = *foreground->channel_energies();
  std::vector<float> channel_counts( nbin );
  for( size_t i = 0; i < nbin; ++i )
    channel_counts[i] = foreground->gamma_channel_content( start_ch + i );

  const auto trial_chi2 = [&]( const PeakContinuum::OffsetType cont_type ) -> double
  {
    const std::vector<PeakDef> no_fixed_peaks;
    std::vector<double> amps, cont_coeffs, amp_uncerts, cont_uncerts;
    try
    {
      return fit_amp_and_offset( &channel_energies[start_ch], channel_counts.data(), nbin,
                                 cont_type, roi_lower, means, sigmas, no_fixed_peaks,
                                 PeakDef::SkewType::NoSkew, nullptr,
                                 amps, cont_coeffs, amp_uncerts, cont_uncerts );
    }catch( const std::exception & )
    {
      return std::numeric_limits<double>::max();
    }
  };//trial_chi2 lambda

  const double chi2_poly = trial_chi2( poly_type );
  const double chi2_step = trial_chi2( step_type );

  if( (chi2_poly == std::numeric_limits<double>::max())
      || (chi2_step == std::numeric_limits<double>::max()) )
    return poly_type;

  return ((chi2_step + chi2_margin) < chi2_poly) ? step_type : poly_type;
}//trial_step_continuum


// Due to the large size of the remaining helper functions, they will be added
// by copying from FitPeaksForNuclideDev.cpp. For now, adding stubs that will
// be replaced with full implementations.


/** Find the energy of the minimum in a smoothed spectrum between two energies.
 Useful for determining the natural "valley" between two peaks for splitting ROIs.

 @param foreground The spectrum to search (if null, returns the midpoint of constraints)
 @param search_lower Lower bound of the search region (keV)
 @param search_upper Upper bound of the search region (keV)
 @param fwhm Approximate FWHM at the search region (used for smoothing window)
 @param constraint_lower The split point must not go below this energy (e.g., to keep
        the left ROI valid). Pass search_lower if no constraint.
 @param constraint_upper The split point must not go above this energy. Pass search_upper
        if no constraint.
 @returns Energy of the minimum in the smoothed spectrum within constraints.
*/
double find_spectrum_valley(
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const double search_lower,
  const double search_upper,
  const double fwhm,
  const double constraint_lower,
  const double constraint_upper )
{
  const double fallback = 0.5 * (std::max( search_lower, constraint_lower )
                                + std::min( search_upper, constraint_upper ));

  if( !foreground || (search_lower >= search_upper) || (fwhm <= 0.0) )
    return fallback;

  const size_t num_channels = foreground->num_gamma_channels();
  if( num_channels < 16 )
    return fallback;

  const size_t lower_ch = foreground->find_gamma_channel( static_cast<float>( search_lower ) );
  const size_t upper_ch = foreground->find_gamma_channel( static_cast<float>( search_upper ) );

  if( (lower_ch >= upper_ch) || (upper_ch >= num_channels) )
    return fallback;

  // Smoothing half-window: use roughly half-FWHM worth of channels
  const double channel_width = (foreground->gamma_channel_upper( upper_ch )
                              - foreground->gamma_channel_lower( lower_ch ))
                              / static_cast<double>( upper_ch - lower_ch );

  size_t smooth_half = (channel_width > 0.0)
    ? std::max( size_t(1), static_cast<size_t>( 0.5 * fwhm / channel_width ) )
    : size_t(1);
  smooth_half = std::min( smooth_half, (upper_ch - lower_ch) / 4 );

  if( lower_ch + smooth_half > upper_ch - smooth_half )
    return fallback;

  double min_smoothed = std::numeric_limits<double>::max();
  size_t min_channel = (lower_ch + upper_ch) / 2;

  for( size_t ch = lower_ch + smooth_half; ch <= upper_ch - smooth_half; ++ch )
  {
    // Check that the center of this channel is within constraints
    const double ch_center = 0.5 * (foreground->gamma_channel_lower( ch )
                                   + foreground->gamma_channel_upper( ch ));
    if( (ch_center < constraint_lower) || (ch_center > constraint_upper) )
      continue;

    double sum = 0.0;
    for( size_t k = ch - smooth_half; k <= ch + smooth_half; ++k )
      sum += foreground->gamma_channel_content( k );

    if( sum < min_smoothed )
    {
      min_smoothed = sum;
      min_channel = ch;
    }
  }//for( size_t ch ... )

  const double result = 0.5 * (foreground->gamma_channel_lower( min_channel )
                              + foreground->gamma_channel_upper( min_channel ));

  // Final safety: clamp to constraints
  return std::clamp( result, constraint_lower, constraint_upper );
}//find_spectrum_valley(...)


std::vector<std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo>> cluster_gammas_to_rois(
    const std::vector<std::function<double(double)>> &rel_eff_fcns,
    const std::vector<std::vector<std::tuple<RelActCalcAuto::SrcVariant, double /*age*/, double/*act*/>>> &sources_age_activity_sets,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
    const std::vector<float> &fwhm_coefficients,
    const double fwhm_lower_energy,
    const double fwhm_upper_energy,
    const double lowest_energy,
    const double highest_energy,
    const GammaClusteringSettings &settings,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks = {},
    std::vector<MarginalRejectedCluster> *marginal_rejects = nullptr,
    const std::vector<PredictedGamma> *supplied_predicted_gammas = nullptr,
    const std::function<double(double)> *fwhm_override = nullptr,
    std::vector<PredictedGamma> *all_predicted_gammas = nullptr,
    const std::string &shadow_stage = std::string(),
    const detail::GlobalContinuumEstimate *shadow_global_override = nullptr,
    std::vector<AutomaticRoiDecisionDiagnostic> *roi_policy_diagnostics = nullptr,
    const bool use_source_evidence_pruning = false )
{
  assert( rel_eff_fcns.size() == sources_age_activity_sets.size() );
  if( rel_eff_fcns.size() != sources_age_activity_sets.size() )
    throw runtime_error( "cluster_gammas_to_rois: there is a different number of relative efficiency functions and sets of sources" );

  vector<std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo>> result_rois;
  if( marginal_rejects )
    marginal_rejects->clear();

  // Keep source and curve provenance with every prediction.  R2 uses it to apply per-source
  // extrapolation guards; the accepted ROI construction below still consumes the same energies
  // and counts as before.
  vector<PredictedGamma> gammas_by_counts;

  if( supplied_predicted_gammas )
    gammas_by_counts = *supplied_predicted_gammas;

  for( size_t rel_eff_index = 0;
       !supplied_predicted_gammas && (rel_eff_index < rel_eff_fcns.size()); ++rel_eff_index )
  {
    const function<double(double)> &rel_eff_fcn = rel_eff_fcns[rel_eff_index];
    const vector<tuple<RelActCalcAuto::SrcVariant,double,double>> &src_age_and_activities = sources_age_activity_sets[rel_eff_index];

    for( const tuple<RelActCalcAuto::SrcVariant,double,double> &src_age_act : src_age_and_activities )
    {
      const RelActCalcAuto::SrcVariant &src = get<0>(src_age_act);
      const double age = get<1>(src_age_act);
      const double activity = get<2>(src_age_act);

      if( RelActCalcAuto::is_null(src) || (activity <= 0.0) )
        continue;

      if( should_debug_print() )
      {
        std::cerr << "cluster_gammas_to_rois: Source " << RelActCalcAuto::to_name( src ) << ", activity=" << activity
        << ", age=" << (age / PhysicalUnits::second) << " seconds ("
        << (age / PhysicalUnits::year) << " years)" << std::endl;
      }

      const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src, activity, age );

      if( should_debug_print() )
      {
        std::cerr << "  " << photons.size() << " photons from " << RelActCalcAuto::to_name( src ) << ", energy range ["
        << lowest_energy << ", " << highest_energy << "] keV" << std::endl;
      }

      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        if( (photon.energy < lowest_energy)
           || (photon.energy > highest_energy)
           || (photon.numPerSecond <= std::numeric_limits<double>::epsilon()) )
        {
          continue;
        }

        const double rel_eff = rel_eff_fcn( photon.energy );
        if( rel_eff <= 0.0 )
          continue;

        gammas_by_counts.push_back( PredictedGamma{ photon.energy,
            photon.numPerSecond * rel_eff, src, rel_eff_index } );
      }//for( const SandiaDecay::EnergyRatePair &photon : photons )
    }//for( const auto &src_act : sources_and_activities )
  }//for( size_t rel_eff_index = 0; rel_eff_index < rel_eff_fcns.size(); ++rel_eff_index )

  if( all_predicted_gammas )
    *all_predicted_gammas = gammas_by_counts;

  if( gammas_by_counts.empty() )
    return result_rois;

  // Create a copy sorted by energy for efficient lookup
  std::vector<PredictedGamma> gammas_by_energy = gammas_by_counts;
  const auto lessThanByEnergy = []( const PredictedGamma &lhs, const PredictedGamma &rhs ) {
    if( lhs.energy != rhs.energy )
      return lhs.energy < rhs.energy;
    if( lhs.rel_eff_curve_index != rhs.rel_eff_curve_index )
      return lhs.rel_eff_curve_index < rhs.rel_eff_curve_index;
    return RelActCalcAuto::to_name(lhs.source) < RelActCalcAuto::to_name(rhs.source);
  };
  std::sort( std::begin(gammas_by_energy), std::end(gammas_by_energy), lessThanByEnergy );

  // Sort gammas by expected counts (highest first), with energy as tiebreaker
  // for determinism when two gammas have equal expected counts.
  std::sort( std::begin(gammas_by_counts), std::end(gammas_by_counts),
    []( const PredictedGamma &lhs, const PredictedGamma &rhs ) {
      if( lhs.expected_counts != rhs.expected_counts )
        return lhs.expected_counts > rhs.expected_counts;
      if( lhs.energy != rhs.energy )
        return lhs.energy < rhs.energy;
      if( lhs.rel_eff_curve_index != rhs.rel_eff_curve_index )
        return lhs.rel_eff_curve_index < rhs.rel_eff_curve_index;
      return RelActCalcAuto::to_name(lhs.source) < RelActCalcAuto::to_name(rhs.source);
  } );

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: Input gammas (" << gammas_by_counts.size() << " total):" << std::endl;
    std::cerr << "  Top 20 by expected counts:" << std::endl;
    for( size_t i = 0; i < std::min( gammas_by_counts.size(), size_t(20) ); ++i )
      std::cerr << "    " << gammas_by_counts[i].energy << " keV, est_counts="
                << gammas_by_counts[i].expected_counts << std::endl;
  }

  // Cluster gamma lines
  std::vector<ClusteredGammaInfo> clustered_gammas;

  // Check if we have a valid FWHM energy range for clamping
  const bool have_fwhm_range = ((fwhm_lower_energy > 0.0)
                                && (fwhm_upper_energy > 0.0)
                                && (fwhm_lower_energy < fwhm_upper_energy));

  // FWHM at an energy, clamped to the valid FWHM-fit range - shared by the adaptive-extent and
  // clean-gap-merge code below.
  const auto fwhm_at = [have_fwhm_range, fwhm_lower_energy, fwhm_upper_energy, fwhm_form,
                        &fwhm_coefficients, fwhm_override]( const double energy ) -> double
  {
    if( fwhm_override )
      return (*fwhm_override)( energy );
    const double e = have_fwhm_range
        ? std::clamp( energy, fwhm_lower_energy, fwhm_upper_energy ) : energy;
    return DetectorPeakResponse::peakResolutionFWHM( static_cast<float>(e), fwhm_form, fwhm_coefficients );
  };//fwhm_at lambda

  // NOTE: the keep-gate below no longer has a "rescue" fallback for the all-rejected case.  A source
  // whose every predicted cluster is sub-threshold simply yields no clustered ROIs here; the
  // data-confirmed found+matched auto-search peaks are given tight ROIs by the caller
  // (seed_tight_rois_for_found_peaks), and a genuinely-empty ROI set is a valid empty result rather
  // than a setup failure (see fit_peaks_for_nuclide_relactauto).  [architecture review 2026-07-18]
  for( const PredictedGamma &energy_counts : gammas_by_counts )
  {
    auto ene_pos = std::lower_bound( std::begin(gammas_by_energy), std::end(gammas_by_energy),
                                    energy_counts, lessThanByEnergy );
    if( ene_pos == std::end(gammas_by_energy) )
      continue;

    if( ene_pos->energy != energy_counts.energy )
    {
      // Already removed from gammas_by_energy (absorbed into another cluster)
      if( should_debug_print() && (energy_counts.energy >= 800.0) && (energy_counts.energy <= 820.0) )
      {
        std::cerr << "cluster_gammas_to_rois: Gamma at " << energy_counts.energy
             << " keV already absorbed into another cluster" << std::endl;
      }
      continue;
    }

    const double energy = energy_counts.energy;
    const double counts = energy_counts.expected_counts;

    const double fwhm = fwhm_at( energy );

    // Check for invalid FWHM (NaN, zero, or negative)
    if( !std::isfinite(fwhm) || (fwhm <= 0.0f) )
    {
      if( should_debug_print() )
        std::cerr << "Warning: Invalid FWHM=" << fwhm << " at energy=" << energy << " keV, skipping gamma line" << std::endl;
      continue;
    }

    const double sigma = fwhm / PhysicalUnits::fwhm_nsigma;

    const double lower = std::max( lowest_energy, energy - settings.cluster_num_sigma * sigma );
    const double upper = std::min( highest_energy, energy + settings.cluster_num_sigma * sigma );

    // Find all gammas in this range
    const auto start_remove = std::lower_bound( std::begin(gammas_by_energy),
      std::end(gammas_by_energy), lower,
      []( const PredictedGamma &gamma, const double value ) {
        return gamma.energy < value;
      } );
    const auto end_remove = std::upper_bound( std::begin(gammas_by_energy),
      std::end(gammas_by_energy), upper,
      []( const double value, const PredictedGamma &gamma ) {
        return value < gamma.energy;
      } );

    const double counts_in_region = std::accumulate( start_remove, end_remove, 0.0,
        []( const double sum, const PredictedGamma &gamma ) {
          return sum + gamma.expected_counts;
    } );

    // Capture the gamma lines before erasing them - as separate energy and amplitude arrays
    std::vector<double> gamma_energies_in_cluster;
    std::vector<double> gamma_amplitudes_in_cluster;
    std::vector<PredictedGamma> predicted_gammas_in_cluster;
    for( auto it = start_remove; it != end_remove; ++it )
    {
      gamma_energies_in_cluster.push_back( it->energy );
      gamma_amplitudes_in_cluster.push_back( it->expected_counts );
      predicted_gammas_in_cluster.push_back( *it );
    }

    // Keep-gate significance: z = S / sqrt(S + B).  S is the cluster's total expected counts; B
    // is the sideband-estimated continuum over the cluster's CORE extent - the outermost gammas
    // +/- roi_core_num_fwhm x FWHM, i.e. the always-included part of the ROI the fit will actually
    // see (formerly the +/- cluster_num_sigma seed window, which is neither the fit window nor
    // guaranteed to contain all the signal).  The sideband samples have the cluster's own
    // predicted tails subtracted and slide away from interfering unfit auto-search peaks, so a
    // busy continuum no longer biases B high and suppresses genuine weak lines.  A fixed minimum
    // expected-count floor protects the Gaussian-statistics regime.
    const double core_lo = std::max( lowest_energy,
        gamma_energies_in_cluster.front() - settings.roi_core_num_fwhm * fwhm );
    const double core_hi = std::min( highest_energy,
        gamma_energies_in_cluster.back() + settings.roi_core_num_fwhm * fwhm );

    const auto cluster_predicted_signal = [&]( const double x0, const double x1 ) -> double {
      return detail::predicted_gaussian_counts( gamma_energies_in_cluster,
                                                gamma_amplitudes_in_cluster, fwhm_at, x0, x1 );
    };//cluster_predicted_signal lambda

    const double data_area = foreground->gamma_integral( static_cast<float>(core_lo),
                                                         static_cast<float>(core_hi) );

    const detail::LocalContinuumEstimate local_cont = detail::estimate_local_continuum(
        foreground, core_lo, core_hi, fwhm, 0.5, cluster_predicted_signal, unfit_auto_peaks );
    // R1 step 2: prefer the shared SNIP global continuum for the keep-gate B; fall back to the local
    // two-sideband estimate (then gross data area) when the global provider is absent/invalid.
    const double b_est = (settings.global_continuum && settings.global_continuum->valid())
                         ? settings.global_continuum->integral( core_lo, core_hi )
                         : (local_cont.valid ? local_cont.integral( core_lo, core_hi ) : data_area);

    gammas_by_energy.erase( start_remove, end_remove );

    const double signif = counts_in_region / std::sqrt( std::max( 1.0, counts_in_region + b_est ) );

    // Additional safety check - ensure lower and upper are finite and valid
    if( !std::isfinite(lower) || !std::isfinite(upper) || (lower >= upper) )
    {
      if( should_debug_print() )
        std::cerr << "Warning: Invalid cluster bounds [" << lower << ", " << upper << "], skipping" << std::endl;
      continue;
    }

    const bool passes_counts = (counts_in_region > sm_keep_gate_min_est_counts);
    const bool passes_signif = (signif > settings.keep_significance_z);
    // Marginal R2 candidates are already past their provisional gate and may mint immediately.
    // Normally accepted clusters defer minting until the over-wide H0/Hs/Hf bridge check below;
    // this makes "rejected before admission" an exact atom-ledger statement rather than a waiver.
    const auto make_cluster_atoms = [&settings]( const std::vector<PredictedGamma> &predictions )
        -> std::vector<detail::RoiAtom> {
      std::vector<detail::RoiAtom> atoms;
      atoms.reserve( predictions.size() );
      for( const PredictedGamma &pg : predictions )
      {
        detail::RoiAtom atom;
        atom.id = detail::next_roi_atom_id();
        atom.energy = pg.energy;
        atom.area = pg.expected_counts;
        atom.kind = detail::RoiAtomKind::ModeledGamma;
        atom.source = pg.source;
        atom.rel_eff_curve_index = pg.rel_eff_curve_index;
        atom.admission = settings.cluster_admission_stage;
        atoms.push_back( atom );
      }
      return atoms;
    };//make_cluster_atoms

    if( passes_counts && passes_signif )
    {
      ClusteredGammaInfo cluster_info;
      cluster_info.lower = lower;
      cluster_info.upper = upper;
      cluster_info.predicted_gammas = std::move( predicted_gammas_in_cluster );
      cluster_info.gamma_energies = std::move( gamma_energies_in_cluster );
      cluster_info.gamma_amplitudes = std::move( gamma_amplitudes_in_cluster );
      clustered_gammas.push_back( std::move( cluster_info ) );
    }else if( marginal_rejects && detail::is_marginal_keep_reject(
              counts_in_region, signif, settings.keep_significance_z ) )
    {
      MarginalRejectedCluster marginal;
      marginal.cluster.lower = lower;
      marginal.cluster.upper = upper;
      marginal.cluster.atoms = make_cluster_atoms( predicted_gammas_in_cluster );
      marginal.cluster.gamma_energies = gamma_energies_in_cluster;
      marginal.cluster.gamma_amplitudes = gamma_amplitudes_in_cluster;
      marginal.predicted_gammas = std::move( predicted_gammas_in_cluster );
      marginal.expected_counts = counts_in_region;
      marginal.background_counts = b_est;
      marginal.keep_significance = signif;
      marginal_rejects->push_back( std::move(marginal) );
    }

    if( should_debug_print() )
    {
      const std::string status_str = (passes_counts && passes_signif) ? "Accepted" : "Rejected";
      std::cerr << "cluster_gammas_to_rois: " << status_str << " [" << std::fixed << std::setprecision(1) << lower << ", " << upper << "] keV (e="
           << energy << " keV): "
           << "est_counts=" << counts_in_region << (passes_counts ? " > " : " < ") << sm_keep_gate_min_est_counts << "; "
         << "z=" << signif << (passes_signif ? " > " : " < ") << settings.keep_significance_z
         << " (b_est=" << b_est << (local_cont.valid ? "" : " gross-fallback") << ", gross=" << data_area << ")"
         << std::endl;
    }
  }//for( const PredictedGamma &energy_counts : gammas_by_counts )

  // Sort by lower energy
  std::sort( std::begin(clustered_gammas), std::end(clustered_gammas),
    []( const ClusteredGammaInfo &lhs, const ClusteredGammaInfo &rhs ) {
      return lhs.lower < rhs.lower;
  } );

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: Initial " << clustered_gammas.size() << " clustered ROIs:" << std::endl;
    for( size_t i = 0; i < clustered_gammas.size(); ++i )
    {
      const ClusteredGammaInfo &c = clustered_gammas[i];
      std::cerr << "  [" << i << "] range=[" << c.lower << ", " << c.upper << "] keV, "
           << c.gamma_energies.size() << " gammas: ";
      for( size_t j = 0; j < std::min( c.gamma_energies.size(), size_t(5) ); ++j )
        std::cerr << c.gamma_energies[j] << (j + 1 < c.gamma_energies.size() ? ", " : "");
      if( c.gamma_energies.size() > 5 )
        std::cerr << "... (" << c.gamma_energies.size() - 5 << " more)";
      std::cerr << std::endl;
    }
  }

  // Update cluster bounds based on weighted mean and effective FWHM
  // The effective FWHM accounts for both individual peak widths and the spread of gamma line energies
  // Uses the law of total variance for a mixture of Gaussians:
  //   σ²_total = (weighted avg of individual variances) + (weighted variance of the means)
  for( ClusteredGammaInfo &cluster : clustered_gammas )
  {
    if( cluster.gamma_energies.empty() )
      continue;
    
    // Calculate weighted mean energy (weighted by expected amplitude)
    double sum_weighted_energy = 0.0;
    double sum_weights = 0.0;
    for( size_t i = 0; i < cluster.gamma_energies.size(); ++i )
    {
      sum_weighted_energy += cluster.gamma_energies[i] * cluster.gamma_amplitudes[i];
      sum_weights += cluster.gamma_amplitudes[i];
    }
    
    if( sum_weights <= 0.0 )
      continue;
    
    const double weighted_mean = sum_weighted_energy / sum_weights;
    
    // Calculate effective variance using law of total variance for mixture of Gaussians
    // σ²_total = E[Var(X|I)] + Var(E[X|I])
    //          = (weighted avg of individual σ²) + (weighted variance of means)
    double sum_weighted_var = 0.0;        // For E[Var(X|I)] - weighted average of individual variances
    double sum_weighted_sq_dev = 0.0;     // For Var(E[X|I]) - weighted variance of the means
    
    for( size_t i = 0; i < cluster.gamma_energies.size(); ++i )
    {
      const double energy = cluster.gamma_energies[i];
      const double amplitude = cluster.gamma_amplitudes[i];
      
      const double fwhm_i = fwhm_at( energy );
      
      if( !std::isfinite(fwhm_i) || (fwhm_i <= 0.0) )
        continue;
      
      const double sigma_i = fwhm_i / PhysicalUnits::fwhm_nsigma;
      const double var_i = sigma_i * sigma_i;
      
      // Weighted average of individual variances
      sum_weighted_var += amplitude * var_i;
      
      // Weighted squared deviation from weighted mean
      const double dev = energy - weighted_mean;
      sum_weighted_sq_dev += amplitude * dev * dev;
    }
    
    // Total variance = weighted avg variance + weighted variance of means
    const double total_var = (sum_weighted_var + sum_weighted_sq_dev) / sum_weights;
    const double effective_sigma = std::sqrt( total_var );
    const double effective_fwhm = effective_sigma * PhysicalUnits::fwhm_nsigma;
    
    if( !std::isfinite(effective_fwhm) || (effective_fwhm <= 0.0) )
      continue;
    
    // Data-driven extent: an always-included core around the outermost gammas, then sideband
    // extension while the data stays statistically consistent with a local linear continuum
    // (see detail::extend_roi_by_sidebands).  Replaces the former fixed weighted_mean +/- k*FWHM.
    const double old_lower = cluster.lower;
    const double old_upper = cluster.upper;

    const detail::AdaptiveExtentResult extent = detail::extend_roi_by_sidebands(
        cluster.gamma_energies, cluster.gamma_amplitudes, effective_fwhm, foreground, fwhm_at,
        unfit_auto_peaks, settings.roi_core_num_fwhm, settings.roi_extend_z,
        settings.roi_max_num_fwhm, settings.skew_type, lowest_energy, highest_energy );

    cluster.lower = extent.lower;
    cluster.upper = extent.upper;

    if( should_debug_print() )
    {
      std::cerr << "cluster_gammas_to_rois: Adaptive extent for cluster with "
      << cluster.gamma_energies.size() << " gammas:" << std::endl
      << "  weighted_mean=" << weighted_mean << " keV, effective_fwhm=" << effective_fwhm << " keV" << std::endl
      << "  old range=[" << old_lower << ", " << old_upper << "] keV"
      << " -> new range=[" << cluster.lower << ", " << cluster.upper << "] keV"
      << " (sidebands " << extent.sideband_lower_kev << " / " << extent.sideband_upper_kev << " keV)" << std::endl;
    }
  }//for( ClusteredGammaInfo &cluster : clustered_gammas )

  // Adaptive sideband extension changes both bounds independently, so the ordering established
  // before extension is no longer an invariant.  The evidence-component scan and every subsequent
  // adjacent reconciliation require increasing lower bounds.
  if( settings.use_automatic_roi_policy )
  {
    std::sort( std::begin(clustered_gammas), std::end(clustered_gammas),
      []( const ClusteredGammaInfo &lhs, const ClusteredGammaInfo &rhs ) {
        return lhs.lower < rhs.lower;
    } );
  }

  // Predicted significance admits provisional source groups, but it cannot establish that the
  // measured spectrum contains them.  After the source-collapse gate has fired, and in an
  // over-wide transitive overlap component only, compare each provisional group transactionally on
  // its own fixed channels: continuum-only (H0), one local scale multiplying the exact source-line
  // mixture (Hs), and all FWHM-distinct significant found features together (Hf).  Rejected groups
  // have not minted atoms yet, so they cannot become artificial bridges and exact-once ownership
  // begins only after this evidence gate.  At least two supported groups must survive or the whole
  // local rejection set is rolled back.
  if( settings.use_automatic_roi_policy && use_source_evidence_pruning
      && (clustered_gammas.size() >= 2) )
  {
    std::vector<bool> retain( clustered_gammas.size(), true );
    size_t component_first = 0;
    while( component_first < clustered_gammas.size() )
    {
      size_t component_last = component_first;
      double component_upper = clustered_gammas[component_first].upper;
      while( (component_last + 1 < clustered_gammas.size())
             && (clustered_gammas[component_last + 1].lower <= component_upper) )
      {
        ++component_last;
        component_upper = std::max( component_upper, clustered_gammas[component_last].upper );
      }
      const double component_lower = clustered_gammas[component_first].lower;
      const double midpoint = 0.5*(component_lower + component_upper);
      const double fwhm = fwhm_at( midpoint );
      const bool overwide = (component_last > component_first) && std::isfinite(fwhm)
          && (fwhm > 0.0) && (settings.max_fwhm_width > 0.0)
          && (((component_upper - component_lower) / fwhm) > settings.max_fwhm_width);
      if( overwide )
      {
        struct PendingEvidence
        {
          size_t cluster_index = 0;
          detail::SourceClusterEvidenceResult evidence;
        };
        std::vector<PendingEvidence> pending;
        size_t supported = 0;
        for( size_t index = component_first; index <= component_last; ++index )
        {
          PendingEvidence item;
          item.cluster_index = index;
          const ClusteredGammaInfo &cluster = clustered_gammas[index];
          item.evidence = detail::evaluate_source_cluster_evidence(
              cluster.gamma_energies, cluster.gamma_amplitudes, cluster.lower, cluster.upper,
              foreground, fwhm_at, unfit_auto_peaks, settings.keep_significance_z,
              settings.roi_core_num_fwhm,
              settings.cont_order_aicc_penalty );
          const bool rejected = (item.evidence.decision
                  == detail::SourceClusterEvidenceDecision::RejectContinuumOnly)
              || (item.evidence.decision
                  == detail::SourceClusterEvidenceDecision::RejectFreeFeature);
          if( !rejected )
            ++supported;
          pending.push_back( std::move(item) );
        }
        const bool apply_rejections = (supported >= 2);
        for( const PendingEvidence &item : pending )
        {
          const detail::SourceClusterEvidenceResult &evidence = item.evidence;
          const bool rejected_continuum = (evidence.decision
              == detail::SourceClusterEvidenceDecision::RejectContinuumOnly);
          const bool rejected_free = (evidence.decision
              == detail::SourceClusterEvidenceDecision::RejectFreeFeature);
          const bool rejected = apply_rejections && (rejected_continuum || rejected_free);
          retain[item.cluster_index] = !rejected;

          AutomaticRoiDecisionDiagnostic diagnostic;
          diagnostic.stage = shadow_stage.empty() ? "automatic clustering source evidence"
                                                   : shadow_stage + " source evidence";
          diagnostic.left_lower = clustered_gammas[item.cluster_index].lower;
          diagnostic.left_upper = clustered_gammas[item.cluster_index].upper;
          diagnostic.calibration_num_channels = foreground->num_gamma_channels();
          diagnostic.combined_width_fwhm = (component_upper - component_lower) / fwhm;
          diagnostic.source_evidence_tested = true;
          diagnostic.source_null_aicc = evidence.null_aicc;
          diagnostic.source_tied_aicc = evidence.source_aicc;
          diagnostic.free_feature_aicc = evidence.free_feature_aicc;
          diagnostic.source_likelihood_z = evidence.source_likelihood_z;
          diagnostic.free_feature_energy = evidence.free_feature_energy;
          if( rejected_continuum && apply_rejections )
            diagnostic.decision = AutomaticRoiDecision::SourceBridgeRejectedContinuum;
          else if( rejected_free && apply_rejections )
            diagnostic.decision = AutomaticRoiDecision::SourceBridgeRejectedFreeFeature;
          else
            diagnostic.decision = AutomaticRoiDecision::SourceBridgeRetained;
          diagnostic.reason = (!apply_rejections && (rejected_continuum || rejected_free))
              ? "local source-evidence rejection rolled back: fewer than two supported anchors"
              : evidence.reason;
          if( roi_policy_diagnostics )
            roi_policy_diagnostics->push_back( diagnostic );
          else
            detail::record_automatic_roi_diagnostic( diagnostic );

          if( should_debug_print() )
          {
            std::cerr << "cluster_gammas_to_rois: "
              << automatic_roi_decision_name(diagnostic.decision)
              << " provisional group [" << diagnostic.left_lower << ", "
              << diagnostic.left_upper << "] keV: H0=" << diagnostic.source_null_aicc
              << ", Hs=" << diagnostic.source_tied_aicc
              << ", Hf=" << diagnostic.free_feature_aicc
              << ", source_z=" << diagnostic.source_likelihood_z
              << "; " << diagnostic.reason << std::endl;
          }
        }
      }
      component_first = component_last + 1;
    }
    std::vector<ClusteredGammaInfo> evidence_filtered;
    evidence_filtered.reserve( clustered_gammas.size() );
    for( size_t index = 0; index < clustered_gammas.size(); ++index )
      if( retain[index] )
        evidence_filtered.push_back( std::move(clustered_gammas[index]) );
    clustered_gammas = std::move( evidence_filtered );
  }

  // Atom admission starts here, after every provisional bridge decision.  Preserve source and
  // rel-eff-curve provenance exactly for the retained groups.
  if( settings.use_automatic_roi_policy )
  {
    for( ClusteredGammaInfo &cluster : clustered_gammas )
    {
      cluster.atoms.clear();
      cluster.atoms.reserve( cluster.predicted_gammas.size() );
      for( const PredictedGamma &prediction : cluster.predicted_gammas )
      {
        detail::RoiAtom atom;
        atom.id = detail::next_roi_atom_id();
        atom.energy = prediction.energy;
        atom.area = prediction.expected_counts;
        atom.kind = detail::RoiAtomKind::ModeledGamma;
        atom.source = prediction.source;
        atom.rel_eff_curve_index = prediction.rel_eff_curve_index;
        atom.admission = settings.cluster_admission_stage;
        cluster.atoms.push_back( std::move(atom) );
      }
    }
  }
  
  // Collect all source gamma energies for checking if an unfit peak matches a source gamma;
  // used during merge prevention below.
  std::vector<double> all_source_gamma_energies;
  if( !unfit_auto_peaks.empty() )
  {
    for( const ClusteredGammaInfo &c : clustered_gammas )
      all_source_gamma_energies.insert( all_source_gamma_energies.end(),
        c.gamma_energies.begin(), c.gamma_energies.end() );
  }

  // Atom-safe conversions between a cluster and a partition component (policy mode only).  The
  // cluster's atom ledger is authoritative; gamma_energies/gamma_amplitudes are regenerated from it
  // after every partition so the downstream emit sees a consistent view.
  const auto cluster_to_component = [&]( const ClusteredGammaInfo &c ) -> detail::AutomaticRoiComponent {
    detail::AutomaticRoiComponent comp;
    comp.lower = c.lower;
    comp.upper = c.upper;
    comp.first_channel = foreground->find_gamma_channel( static_cast<float>(c.lower) );
    const float upper_inside = std::nextafter(
        static_cast<float>(c.upper), static_cast<float>(c.lower) );
    comp.last_channel = foreground->find_gamma_channel( upper_inside );
    comp.joined_groups = c.joined_groups;
    comp.protected_geometry = false;
    comp.continuum_type = PeakContinuum::OffsetType::Linear;
    comp.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    comp.atoms = c.atoms;
    if( comp.atoms.empty() )  // defensive: synthesize provenance-less atoms from the flat view
    {
      for( size_t i = 0; i < c.gamma_energies.size(); ++i )
      {
        detail::RoiAtom atom;
        atom.id = detail::next_roi_atom_id();
        atom.energy = c.gamma_energies[i];
        atom.area = (i < c.gamma_amplitudes.size()) ? c.gamma_amplitudes[i] : 0.0;
        atom.admission = settings.cluster_admission_stage;
        comp.atoms.push_back( atom );
      }
    }
    // Materialize the transaction on closed channel ownership and include every admitted core.
    // This also removes sub-channel roundoff between adaptive floating bounds and core validation.
    for( const detail::RoiAtom &atom : comp.atoms )
    {
      const double core_half_width = detail::atomlayer_core_halfwidth(
          atom.energy, fwhm_at, settings.roi_core_num_fwhm );
      if( !std::isfinite(core_half_width) )
        continue;
      comp.first_channel = std::min( comp.first_channel,
          detail::atomlayer_channel_for_energy(foreground, atom.energy - core_half_width) );
      comp.last_channel = std::max( comp.last_channel,
          detail::atomlayer_channel_for_energy(foreground, atom.energy + core_half_width) );
    }
    comp.lower = foreground->gamma_channel_lower( comp.first_channel );
    comp.upper = foreground->gamma_channel_upper( comp.last_channel );
    return comp;
  };//cluster_to_component
  const auto component_to_cluster = []( const detail::AutomaticRoiComponent &comp ) -> ClusteredGammaInfo {
    ClusteredGammaInfo c;
    c.lower = comp.lower;
    c.upper = comp.upper;
    c.joined_groups = comp.joined_groups;
    c.atoms = comp.atoms;
    for( const detail::RoiAtom &a : comp.atoms )
    {
      c.gamma_energies.push_back( a.energy );
      c.gamma_amplitudes.push_back( a.area );
    }
    return c;
  };//component_to_cluster

  // Merge overlapping clusters, but prevent merging when an unfit auto-search peak lies between
  // the largest gammas of each cluster (the interfering peak would contaminate the combined ROI).
  std::vector<ClusteredGammaInfo> merged_clusters;
  if( settings.use_automatic_roi_policy )
  {
    std::vector<detail::AutomaticRoiComponent> components;
    components.reserve( clustered_gammas.size() );
    for( const ClusteredGammaInfo &cluster : clustered_gammas )
      components.push_back( cluster_to_component(cluster) );

    // Auto-search peaks coincident with any requested-source gamma are modeled evidence, not
    // unmodeled barriers.  Preserve the existing exemption while reconciling the whole list.
    std::vector<std::shared_ptr<const PeakDef>> blocking_unfit_peaks;
    for( const std::shared_ptr<const PeakDef> &peak : unfit_auto_peaks )
    {
      bool matches_source = false;
      for( const double gamma_energy : all_source_gamma_energies )
      {
        const double gamma_fwhm = fwhm_at( gamma_energy );
        if( !std::isfinite(gamma_fwhm) || (gamma_fwhm <= 0.0) )
          continue;
        const double gamma_sigma = gamma_fwhm / PhysicalUnits::fwhm_nsigma;
        if( std::fabs(peak->mean() - gamma_energy) < (settings.cluster_num_sigma * gamma_sigma) )
        {
          matches_source = true;
          break;
        }
      }
      if( !matches_source )
        blocking_unfit_peaks.push_back( peak );
    }

    detail::AutomaticRoiPolicySettings policy_settings;
    policy_settings.merge_tail_z = settings.merge_tail_z;
    policy_settings.merge_clean_gap_fwhm = settings.merge_clean_gap_fwhm;
    policy_settings.continuum_aicc_penalty = settings.cont_order_aicc_penalty;
    policy_settings.peak_core_num_fwhm = settings.roi_core_num_fwhm;
    policy_settings.max_width_fwhm = settings.max_fwhm_width;
    policy_settings.allow_overwide_overlap_partition = false;
    policy_settings.stage = shadow_stage.empty() ? "automatic clustering" : shadow_stage;

    detail::AutomaticRoiPartitionConstraints constraints;
    constraints.lowest_energy = lowest_energy;
    constraints.highest_energy = highest_energy;
    constraints.left_barrier = -std::numeric_limits<double>::infinity();
    constraints.min_width_fwhm = settings.min_fwhm_roi;
    constraints.peak_core_num_fwhm = settings.roi_core_num_fwhm;

    std::vector<AutomaticRoiDecisionDiagnostic> local_diagnostics;
    std::vector<AutomaticRoiDecisionDiagnostic> *diagnostic_sink
      = roi_policy_diagnostics ? roi_policy_diagnostics : &local_diagnostics;
    const detail::AutomaticRoiReconcileResult reconciliation
      = detail::reconcile_automatic_components( std::move(components), foreground,
          settings.global_continuum, fwhm_at, blocking_unfit_peaks, policy_settings,
          constraints, diagnostic_sink );
    if( !roi_policy_diagnostics )
      for( const AutomaticRoiDecisionDiagnostic &diagnostic : local_diagnostics )
        detail::record_automatic_roi_diagnostic( diagnostic );

    if( !reconciliation.valid || !reconciliation.orphaned_atoms.empty() )
    {
      std::ostringstream message;
      message << "cluster_gammas_to_rois: invalid automatic ROI reconciliation";
      if( !reconciliation.failure_reason.empty() )
        message << ": " << reconciliation.failure_reason;
      if( !reconciliation.orphaned_atoms.empty() )
        message << " (" << reconciliation.orphaned_atoms.size() << " orphaned atoms)";
      std::cerr << "ERROR: " << message.str() << std::endl;
      throw std::runtime_error( message.str() );
    }

    merged_clusters.reserve( reconciliation.components.size() );
    for( const detail::AutomaticRoiComponent &component : reconciliation.components )
      merged_clusters.push_back( component_to_cluster(component) );
  }else
  {
  for( const ClusteredGammaInfo &cluster : clustered_gammas )
  {
    if( merged_clusters.empty() || (cluster.lower > merged_clusters.back().upper) )
    {
      merged_clusters.push_back( cluster );
    }
    else
    {
      // Clusters overlap - find the dominant gamma in each cluster for comparison
      const ClusteredGammaInfo &prev = merged_clusters.back();
      assert( !prev.gamma_amplitudes.empty() && !cluster.gamma_amplitudes.empty() );

      const size_t prev_max_idx = static_cast<size_t>(
        std::max_element( prev.gamma_amplitudes.begin(), prev.gamma_amplitudes.end() )
        - prev.gamma_amplitudes.begin() );
      const double prev_largest_energy = prev.gamma_energies[prev_max_idx];

      const size_t curr_max_idx = static_cast<size_t>(
        std::max_element( cluster.gamma_amplitudes.begin(), cluster.gamma_amplitudes.end() )
        - cluster.gamma_amplitudes.begin() );
      const double curr_largest_energy = cluster.gamma_energies[curr_max_idx];

      // Check if an unfit peak lies between the largest gammas of each cluster
      bool unfit_peak_between = false;
      std::vector<std::shared_ptr<const PeakDef>> blocking_unfit_peaks;
      if( !unfit_auto_peaks.empty() )
      {
        const double between_lo = std::min( prev_largest_energy, curr_largest_energy );
        const double between_hi = std::max( prev_largest_energy, curr_largest_energy );

        for( const std::shared_ptr<const PeakDef> &peak : unfit_auto_peaks )
        {
          const double peak_mean = peak->mean();

          // Check if this unfit peak matches any source gamma - if so, it's not interfering
          bool matches_source = false;
          for( const double gamma_energy : all_source_gamma_energies )
          {
            const double gamma_fwhm = fwhm_at( gamma_energy );
            if( !std::isfinite( gamma_fwhm ) || (gamma_fwhm <= 0.0f) )
              continue;
            const double gamma_sigma = gamma_fwhm / PhysicalUnits::fwhm_nsigma;
            if( std::fabs( peak_mean - gamma_energy ) < (settings.cluster_num_sigma * gamma_sigma) )
            {
              matches_source = true;
              break;
            }
          }//for( const double gamma_energy : all_source_gamma_energies )

          if( matches_source )
            continue;

          blocking_unfit_peaks.push_back( peak );
          if( (peak_mean <= between_lo) || (peak_mean >= between_hi) )
            continue;
          unfit_peak_between = true;

          if( should_debug_print() )
          {
            std::cerr << "cluster_gammas_to_rois: NOT merging clusters ["
                 << prev.lower << ", " << prev.upper << "] and ["
                 << cluster.lower << ", " << cluster.upper
                 << "] due to unfit peak at " << peak_mean
                 << " keV between largest gammas at " << prev_largest_energy
                 << " and " << curr_largest_energy << " keV" << std::endl;
          }
        }//for( const auto &peak : unfit_auto_peaks )
      }//if( !unfit_auto_peaks.empty() )

      // Clean-gap test: keep the clusters separate only when a continuum-anchoring window exists
      // between the dominant gammas where the predicted tail contamination from BOTH clusters is
      // statistically negligible vs the local continuum noise.  Replaces the amplitude-relative
      // tail-fraction test, which ignored counting statistics.
      bool have_clean_gap = false;
      double clean_win_lo = 0.0, clean_win_hi = 0.0;
      if( !unfit_peak_between )
      {
        have_clean_gap = detail::find_clean_gap_between(
            prev.gamma_energies, prev.gamma_amplitudes,
            cluster.gamma_energies, cluster.gamma_amplitudes,
            std::min( prev_largest_energy, curr_largest_energy ),
            std::max( prev_largest_energy, curr_largest_energy ),
            foreground, fwhm_at, settings.merge_tail_z, settings.merge_clean_gap_fwhm,
            &clean_win_lo, &clean_win_hi, settings.global_continuum );

        if( have_clean_gap && should_debug_print() )
        {
          std::cerr << "cluster_gammas_to_rois: NOT merging clusters ["
               << prev.lower << ", " << prev.upper << "] and ["
               << cluster.lower << ", " << cluster.upper
               << "] - clean gap [" << clean_win_lo << ", " << clean_win_hi << "] keV"
               << " (merge_tail_z=" << settings.merge_tail_z
               << ", clean_gap_fwhm=" << settings.merge_clean_gap_fwhm << ")" << std::endl;
        }
      }//if( !unfit_peak_between )

      if( settings.use_automatic_roi_policy )
      {
        // POLICY MODE: replace the greedy split/merge below with ONE atom-safe partition of this
        // adjacent pair.  No admitted gamma can be clipped-and-dropped; when no core-safe boundary
        // exists the pair merges (never a dropped side).  The legacy split/merge code below runs
        // only for R6-enabled fits (use_automatic_roi_policy == false).
        detail::AutomaticRoiPolicySettings policy_settings;
        policy_settings.merge_tail_z = settings.merge_tail_z;
        policy_settings.merge_clean_gap_fwhm = settings.merge_clean_gap_fwhm;
        policy_settings.continuum_aicc_penalty = settings.cont_order_aicc_penalty;
        policy_settings.peak_core_num_fwhm = settings.roi_core_num_fwhm;
        policy_settings.max_width_fwhm = settings.max_fwhm_width;
        policy_settings.stage = shadow_stage.empty() ? "automatic clustering" : shadow_stage;

        detail::AutomaticRoiPartitionConstraints cons;
        cons.lowest_energy = lowest_energy;
        cons.highest_energy = highest_energy;
        cons.left_barrier = (merged_clusters.size() >= 2)
            ? merged_clusters[merged_clusters.size() - 2].upper
            : -std::numeric_limits<double>::infinity();
        cons.min_width_fwhm = settings.min_fwhm_roi;
        cons.peak_core_num_fwhm = settings.roi_core_num_fwhm;

        const detail::AutomaticRoiComponent left_comp = cluster_to_component( merged_clusters.back() );
        const detail::AutomaticRoiComponent right_comp = cluster_to_component( cluster );
        const detail::AutomaticRoiPartitionResult pr = detail::partition_automatic_roi_pair(
            left_comp, right_comp, foreground, settings.global_continuum, fwhm_at,
            blocking_unfit_peaks, policy_settings, cons );
        if( roi_policy_diagnostics )
          roi_policy_diagnostics->push_back( pr.policy.diagnostic );
        else
          detail::record_automatic_roi_diagnostic( pr.policy.diagnostic );
        assert( pr.orphaned_atoms.empty() );  // clustering has no protected geometry -> no orphans

        if( (pr.outcome == detail::AutomaticRoiPartitionOutcome::KeptSeparate)
            && (pr.components.size() == 2) )
        {
          merged_clusters.back() = component_to_cluster( pr.components[0] );
          merged_clusters.push_back( component_to_cluster( pr.components[1] ) );
        }else if( !pr.components.empty() )
        {
          merged_clusters.back() = component_to_cluster( pr.components[0] );  // Merged
        }
        continue;
      }//if( settings.use_automatic_roi_policy )

      detail::AutomaticRoiGroup left_group;
      left_group.lower = prev.lower;
      left_group.upper = prev.upper;
      left_group.peak_energies = prev.gamma_energies;
      left_group.peak_areas = prev.gamma_amplitudes;
      left_group.joined_groups = prev.joined_groups;
      detail::AutomaticRoiGroup right_group;
      right_group.lower = cluster.lower;
      right_group.upper = cluster.upper;
      right_group.peak_energies = cluster.gamma_energies;
      right_group.peak_areas = cluster.gamma_amplitudes;
      right_group.joined_groups = cluster.joined_groups;
      detail::AutomaticRoiPolicySettings policy_settings;
      policy_settings.merge_tail_z = settings.merge_tail_z;
      policy_settings.merge_clean_gap_fwhm = settings.merge_clean_gap_fwhm;
      policy_settings.continuum_aicc_penalty = settings.cont_order_aicc_penalty;
      policy_settings.peak_core_num_fwhm = settings.roi_core_num_fwhm;
      policy_settings.max_width_fwhm = settings.max_fwhm_width;
      policy_settings.allow_overwide_overlap_partition = use_source_evidence_pruning;
      policy_settings.stage = shadow_stage.empty() ? "automatic clustering" : shadow_stage;
      detail::AutomaticRoiPolicyResult policy;
      if( settings.use_automatic_roi_policy )
      {
        policy = detail::evaluate_automatic_roi_boundary(
            left_group, right_group, foreground, settings.global_continuum, fwhm_at,
            blocking_unfit_peaks, policy_settings );
        if( roi_policy_diagnostics )
          roi_policy_diagnostics->push_back( policy.diagnostic );
        else
          detail::record_automatic_roi_diagnostic( policy.diagnostic );
      }
      const bool keep_separate = settings.use_automatic_roi_policy
        ? ((policy.decision == AutomaticRoiDecision::KeepSeparate)
            || (policy.decision == AutomaticRoiDecision::UnmodeledFeatureBlocked)
            || (policy.decision == AutomaticRoiDecision::ProtectedGeometry))
        : (unfit_peak_between || have_clean_gap);

      if( keep_separate )
      {
        // Don't merge - find the natural valley in the overlap for the split point.
        // Constrain so both clusters still contain their dominant gamma, and to the clean
        // window when one was found (that is where a continuum can actually be anchored).
        const double overlap_lo = cluster.lower;
        const double overlap_hi = merged_clusters.back().upper;
        double split_lo = std::max( overlap_lo, prev_largest_energy );
        double split_hi = std::min( overlap_hi, curr_largest_energy );
        if( settings.use_automatic_roi_policy && (policy.boundary_energy > 0.0) )
        {
          split_lo = std::max( split_lo, policy.boundary_energy );
          split_hi = std::min( split_hi, policy.boundary_energy );
        }

        double split_point;
        if( settings.use_automatic_roi_policy
            && (policy.boundary_energy >= split_lo) && (policy.boundary_energy <= split_hi) )
        {
          split_point = policy.boundary_energy;
        }
        else if( split_lo >= split_hi )
        {
          split_point = 0.5 * (overlap_lo + overlap_hi);
        }
        else
        {
          const double gap_mid = 0.5 * (prev_largest_energy + curr_largest_energy);
          const double gap_fwhm = fwhm_at( gap_mid );
          const double search_fwhm = (std::isfinite( gap_fwhm ) && (gap_fwhm > 0.0))
              ? gap_fwhm : (overlap_hi - overlap_lo);

          split_point = find_spectrum_valley( foreground,
              overlap_lo, overlap_hi, search_fwhm,
              split_lo, split_hi );
        }

        ClusteredGammaInfo adjusted_cluster = cluster;
        if( settings.use_automatic_roi_policy
            && (policy.decision == AutomaticRoiDecision::UnmodeledFeatureBlocked)
            && (policy.exclusion_upper > policy.exclusion_lower) )
        {
          const size_t lower_channel = foreground->find_gamma_channel(
              static_cast<float>(policy.exclusion_lower) );
          const size_t upper_channel = foreground->find_gamma_channel(
              static_cast<float>(policy.exclusion_upper) );
          merged_clusters.back().upper = foreground->gamma_channel_lower( lower_channel );
          adjusted_cluster.lower = foreground->gamma_channel_upper( upper_channel );
        }else if( settings.use_automatic_roi_policy )
        {
          const size_t split_channel = foreground->find_gamma_channel(
              static_cast<float>(split_point) );
          merged_clusters.back().upper = foreground->gamma_channel_lower( split_channel );
          adjusted_cluster.lower = foreground->gamma_channel_upper( split_channel );
        }else
        {
          merged_clusters.back().upper = split_point;
          adjusted_cluster.lower = split_point;
        }
        const bool left_valid = (merged_clusters.back().upper > merged_clusters.back().lower)
            && (prev_largest_energy <= merged_clusters.back().upper);
        const bool right_valid = (adjusted_cluster.upper > adjusted_cluster.lower)
            && (curr_largest_energy >= adjusted_cluster.lower);
        if( !left_valid )
        {
          merged_clusters.pop_back();
          if( right_valid )
            merged_clusters.push_back( adjusted_cluster );
          continue;
        }
        if( !right_valid )
          continue;  // retain the anchored left child and reject the new automatic group
        merged_clusters.push_back( adjusted_cluster );
      }
      else
      {
        if( should_debug_print() )
        {
          std::cerr << "cluster_gammas_to_rois: Merging cluster [" << cluster.lower << ", " << cluster.upper
               << "] into [" << merged_clusters.back().lower << ", " << merged_clusters.back().upper << "]"
               << " -> new upper=" << std::max( merged_clusters.back().upper, cluster.upper ) << std::endl;
        }

        merged_clusters.back().upper = std::min( highest_energy, std::max( merged_clusters.back().upper, cluster.upper ) );
        merged_clusters.back().gamma_energies.insert(
          std::end(merged_clusters.back().gamma_energies),
          std::begin(cluster.gamma_energies),
          std::end(cluster.gamma_energies)
        );
        merged_clusters.back().gamma_amplitudes.insert(
          std::end(merged_clusters.back().gamma_amplitudes),
          std::begin(cluster.gamma_amplitudes),
          std::end(cluster.gamma_amplitudes)
        );
        merged_clusters.back().joined_groups += cluster.joined_groups;
      }//if( keep_separate ) / else
    }
  }//for( const ClusteredGammaInfo &cluster : clustered_gammas )
  }//if( settings.use_automatic_roi_policy ) / else

  // Validate merged clusters - ensure merge didn't introduce NaN values
  merged_clusters.erase(
    std::remove_if( std::begin(merged_clusters), std::end(merged_clusters),
      [&]( const ClusteredGammaInfo &cluster ) {
        const bool invalid = !std::isfinite(cluster.lower) || !std::isfinite(cluster.upper) || (cluster.lower >= cluster.upper);
        if( invalid && should_debug_print() )
          std::cerr << "Warning: Removing invalid merged cluster with bounds [" << cluster.lower << ", " << cluster.upper << "]" << std::endl;
        // A validated atom-safe partition never yields NaN/inverted bounds, so an erased cluster in
        // policy mode would be a silent atom loss - assert it carries no atoms.
        assert( !invalid || !settings.use_automatic_roi_policy || cluster.atoms.empty() );
        return invalid;
      } ),
    std::end(merged_clusters)
  );

  // Greedy overlap collection establishes the transitive component, but after a measured
  // source-anchor collapse it is not allowed to choose that component's final boundaries.  Keep
  // this geometry fallback inside the same evidence-gated recovery trial; applying it during the
  // initial ordinary solve changed unrelated dense/unshielded spectra before any collapse existed.
  // For every over-wide union in that trial, score the full component and transactionally replace
  // it only when the best core-safe two-ROI partition beats the union plus soft-width pressure.
  if( settings.use_automatic_roi_policy && use_source_evidence_pruning )
  {
    std::vector<ClusteredGammaInfo> jointly_partitioned;
    for( const ClusteredGammaInfo &cluster : merged_clusters )
    {
      const double cluster_midpoint = 0.5*(cluster.lower + cluster.upper);
      const double cluster_fwhm = fwhm_at( cluster_midpoint );
      const bool overwide = std::isfinite(cluster_fwhm) && (cluster_fwhm > 0.0)
          && (settings.max_fwhm_width > 0.0)
          && (((cluster.upper - cluster.lower) / cluster_fwhm) > settings.max_fwhm_width);
      if( !overwide )
      {
        jointly_partitioned.push_back( cluster );
        continue;
      }
      const detail::AutomaticRoiComponent component = cluster_to_component( cluster );
      detail::AutomaticRoiPolicySettings policy_settings;
      policy_settings.merge_tail_z = settings.merge_tail_z;
      policy_settings.merge_clean_gap_fwhm = settings.merge_clean_gap_fwhm;
      policy_settings.continuum_aicc_penalty = settings.cont_order_aicc_penalty;
      policy_settings.peak_core_num_fwhm = settings.roi_core_num_fwhm;
      policy_settings.max_width_fwhm = settings.max_fwhm_width;
      policy_settings.stage = shadow_stage.empty() ? "automatic clustering whole component"
                                                   : shadow_stage + " whole component";
      detail::AutomaticRoiPartitionConstraints constraints;
      constraints.lowest_energy = lowest_energy;
      constraints.highest_energy = highest_energy;
      constraints.left_barrier = jointly_partitioned.empty()
        ? -std::numeric_limits<double>::infinity() : jointly_partitioned.back().upper;
      constraints.min_width_fwhm = settings.min_fwhm_roi;
      constraints.peak_core_num_fwhm = settings.roi_core_num_fwhm;
      const detail::AutomaticRoiComponentPartitionResult partition
        = detail::partition_overwide_automatic_component( { component }, foreground, fwhm_at,
            unfit_auto_peaks, policy_settings, constraints );
      if( roi_policy_diagnostics )
        roi_policy_diagnostics->push_back( partition.diagnostic );
      else
        detail::record_automatic_roi_diagnostic( partition.diagnostic );
      if( partition.valid && partition.changed )
      {
        const uint64_t partition_id = detail::next_roi_atom_id();
        for( const detail::AutomaticRoiComponent &child : partition.components )
        {
          ClusteredGammaInfo materialized = component_to_cluster( child );
          materialized.selected_partition_id = partition_id;
          materialized.selected_parent_lower = cluster.lower;
          materialized.selected_parent_upper = cluster.upper;
          jointly_partitioned.push_back( std::move(materialized) );
        }
      }else
      {
        jointly_partitioned.push_back( cluster );
      }
      if( should_debug_print() && (partition.diagnostic.width_pressure > 0.0) )
      {
        std::cerr << "cluster_gammas_to_rois: whole-component "
          << automatic_roi_decision_name(partition.diagnostic.decision)
          << " [" << cluster.lower << ", " << cluster.upper << "] keV: union_AICc="
          << partition.diagnostic.one_roi_aicc << ", partition_AICc="
          << partition.diagnostic.two_roi_aicc << ", width_pressure="
          << partition.diagnostic.width_pressure << "; "
          << partition.diagnostic.reason << std::endl;
      }
    }
    merged_clusters = std::move( jointly_partitioned );
  }

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: After merging, " << merged_clusters.size() << " clusters:" << std::endl;
    for( size_t i = 0; i < merged_clusters.size(); ++i )
    {
      const ClusteredGammaInfo &c = merged_clusters[i];
      std::cerr << "  [" << i << "] range=[" << c.lower << ", " << c.upper << "] keV ("
           << (c.upper - c.lower) << " keV wide), " << c.gamma_energies.size() << " gammas" << std::endl;
    }
  }

  std::vector<ClusteredGammaInfo> final_clusters;
  if( settings.use_automatic_roi_policy )
  {
    // Width is a soft pressure in the shared policy; an explicitly inseparable wide group is not
    // overridden by the former synthetic-minimum post-hoc splitter.
    final_clusters = merged_clusters;
  }else
  {
  for( ClusteredGammaInfo &cluster : merged_clusters )
  {
    const double mid_energy = 0.5 * (cluster.lower + cluster.upper);

    const double mid_fwhm = fwhm_at( mid_energy );
    const double max_width = settings.max_fwhm_width * mid_fwhm;
    const double current_width = cluster.upper - cluster.lower;

    if( should_debug_print() )
    {
      std::cerr << "cluster_gammas_to_rois: Cluster [" << cluster.lower << ", " << cluster.upper
           << "] keV: width=" << current_width << " keV, mid_fwhm=" << mid_fwhm
           << " keV, max_fwhm_width=" << settings.max_fwhm_width
           << ", max_width=" << max_width << " keV, "
           << (current_width <= max_width ? "NOT breaking" : "NEEDS breaking") << std::endl;
    }

    if( current_width <= max_width )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    // Need to break up this cluster using synthetic spectrum-based breakpoint selection
    // Sort gamma energies and amplitudes together by energy
    {
      // Create index array for sorting
      std::vector<size_t> indices( cluster.gamma_energies.size() );
      std::iota( std::begin(indices), std::end(indices), 0 );
      std::sort( std::begin(indices), std::end(indices),
        [&cluster]( size_t a, size_t b ) {
          return cluster.gamma_energies[a] < cluster.gamma_energies[b];
        } );

      // Reorder both arrays based on sorted indices
      std::vector<double> sorted_energies( cluster.gamma_energies.size() );
      std::vector<double> sorted_amplitudes( cluster.gamma_amplitudes.size() );
      for( size_t i = 0; i < indices.size(); ++i )
      {
        sorted_energies[i] = cluster.gamma_energies[indices[i]];
        sorted_amplitudes[i] = cluster.gamma_amplitudes[indices[i]];
      }
      cluster.gamma_energies = std::move( sorted_energies );
      cluster.gamma_amplitudes = std::move( sorted_amplitudes );
    }

    if( cluster.gamma_energies.size() <= 1 )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    // Get channel range for this cluster
    const size_t start_channel = foreground->find_gamma_channel( static_cast<float>(cluster.lower) );
    const size_t end_channel = foreground->find_gamma_channel( static_cast<float>(cluster.upper) );
    const size_t num_channels = end_channel - start_channel;

    if( num_channels < 3 )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }

    const float *channel_energies = foreground->channel_energies()->data();

    // Build synthetic spectrum from expected Gaussians (same binning as data)
    std::vector<double> synthetic = build_synthetic_spectrum(
      cluster.gamma_energies,
      cluster.gamma_amplitudes,
      fwhm_at,
      channel_energies,
      start_channel,
      num_channels );

    // Find local minima in synthetic spectrum (with significance computed)
    const double start_energy = channel_energies[start_channel];
    const double end_energy = channel_energies[end_channel];
    const bool debug_this_cluster = should_debug_print() &&
      ((start_energy >= 60.0 && end_energy <= 140.0) ||
       (start_energy >= 140.0 && start_energy <= 160.0));

    if( debug_this_cluster )
    {
      std::cerr << "DEBUG: Breaking cluster [" << cluster.lower << ", " << cluster.upper
           << "] keV: start_channel=" << start_channel
           << ", end_channel=" << end_channel
           << ", num_channels=" << num_channels << std::endl;
    }

    std::vector<LocalMinimum> minima = find_synthetic_minima(
      synthetic,
      start_channel,
      foreground,
      channel_energies,
      fwhm_at,
      settings.break_check_fwhm_fraction );

    if( debug_this_cluster )
    {
      std::cerr << "DEBUG: Found " << minima.size() << " minima in cluster ["
           << cluster.lower << ", " << cluster.upper << "] keV:" << std::endl;
      for( size_t i = 0; i < minima.size(); ++i )
      {
        const LocalMinimum &m = minima[i];
        std::cerr << "  [" << i << "] channel=" << m.channel
             << ", energy=" << channel_energies[m.channel] << " keV"
             << ", synthetic_value=" << m.synthetic_value
             << ", depth_score=" << m.depth_score
             << ", stat_sig=" << m.statistical_significance << std::endl;
      }
    }

    // If no minima found, don't break the ROI
    if( minima.empty() )
    {
      final_clusters.push_back( std::move( cluster ) );
      continue;
    }
    
    // Sort by statistical_significance (ascending - least significant first)
    // Use depth_score as tiebreaker (descending - deepest first)
    std::sort( std::begin(minima), std::end(minima),
      []( const LocalMinimum &a, const LocalMinimum &b ) {
        if( a.statistical_significance != b.statistical_significance )
          return a.statistical_significance < b.statistical_significance;
        return a.depth_score > b.depth_score;
      });

    // Select breakpoints, validating significant peaks exist between them
    std::vector<size_t> selected_breakpoint_channels;
    size_t current_lower_ch = start_channel;

    for( const LocalMinimum &candidate : minima )
    {
      // Check if significant peak exists between current_lower_ch and this candidate
      const bool has_left_peak = has_significant_peak_between(
        current_lower_ch,
        candidate.channel,
        synthetic,
        start_channel,
        foreground,
        channel_energies,
        fwhm_at,
        settings.break_check_fwhm_fraction,
        settings.break_peak_significance_threshold );

      // Check if significant peak exists between this candidate and ROI upper
      const bool has_right_peak = has_significant_peak_between(
        candidate.channel,
        end_channel,
        synthetic,
        start_channel,
        foreground,
        channel_energies,
        fwhm_at,
        settings.break_check_fwhm_fraction,
        settings.break_peak_significance_threshold );

      if( debug_this_cluster )
      {
        std::cerr << "DEBUG: Evaluating breakpoint candidate at channel " << candidate.channel
             << " (" << channel_energies[candidate.channel] << " keV): "
             << "has_left_peak=" << (has_left_peak ? "YES" : "NO")
             << ", has_right_peak=" << (has_right_peak ? "YES" : "NO")
             << ", accepted=" << ((has_left_peak && has_right_peak) ? "YES" : "NO") << std::endl;
      }

      if( has_left_peak && has_right_peak )
      {
        selected_breakpoint_channels.push_back( candidate.channel );
        current_lower_ch = candidate.channel;

        // Check if we've satisfied the max_fwhm_width constraint for all sub-ROIs
        bool all_rois_ok = true;
        size_t prev_ch = start_channel;
        for( size_t bp_ch : selected_breakpoint_channels )
        {
          const double sub_lower = channel_energies[prev_ch];
          const double sub_upper = channel_energies[bp_ch];
          const double sub_width = sub_upper - sub_lower;
          const double sub_center = 0.5 * (sub_lower + sub_upper);
          const double sub_fwhm = fwhm_at( sub_center );
          if( sub_width > settings.max_fwhm_width * sub_fwhm )
            all_rois_ok = false;
          prev_ch = bp_ch;
        }
        // Check last sub-ROI
        {
          const double sub_lower = channel_energies[prev_ch];
          const double sub_upper = channel_energies[end_channel];
          const double sub_width = sub_upper - sub_lower;
          const double sub_center = 0.5 * (sub_lower + sub_upper);
          const double sub_fwhm = fwhm_at( sub_center );
          if( sub_width > settings.max_fwhm_width * sub_fwhm )
            all_rois_ok = false;
        }

        if( all_rois_ok )
          break;
      }
    }//for( const LocalMinimum &candidate : minima )

    // Create sub-clusters based on selected breakpoint channels
    if( selected_breakpoint_channels.empty() )
    {
      final_clusters.push_back( std::move( cluster ) );
    }
    else
    {
      std::sort( std::begin(selected_breakpoint_channels), std::end(selected_breakpoint_channels) );

      // Convert breakpoint channels to energies
      std::vector<double> breakpoint_energies;
      for( size_t bp_ch : selected_breakpoint_channels )
        breakpoint_energies.push_back( channel_energies[bp_ch] );

      if( should_debug_print() )
      {
        std::cerr << "cluster_gammas_to_rois: Breaking up cluster [" << cluster.lower << ", " << cluster.upper
             << "] at " << breakpoint_energies.size() << " breakpoints: ";
        for( size_t i = 0; i < breakpoint_energies.size(); ++i )
          std::cerr << breakpoint_energies[i] << (i + 1 < breakpoint_energies.size() ? ", " : "");
        std::cerr << " keV" << std::endl;
      }

      double seg_start = cluster.lower;
      size_t gamma_start_idx = 0;

      for( const double bp_energy : breakpoint_energies )
      {
        ClusteredGammaInfo sub_cluster;
        sub_cluster.lower = std::max( lowest_energy, seg_start );
        sub_cluster.upper = std::min( highest_energy, bp_energy );

        // Add gamma lines that fall within this segment
        for( size_t j = gamma_start_idx; j < cluster.gamma_energies.size(); ++j )
        {
          if( cluster.gamma_energies[j] >= bp_energy )
          {
            gamma_start_idx = j;
            break;
          }
          sub_cluster.gamma_energies.push_back( cluster.gamma_energies[j] );
          sub_cluster.gamma_amplitudes.push_back( cluster.gamma_amplitudes[j] );
          if( j == cluster.gamma_energies.size() - 1 )
            gamma_start_idx = j + 1;
        }

        if( !sub_cluster.gamma_energies.empty() )
          final_clusters.push_back( std::move( sub_cluster ) );

        seg_start = bp_energy;
      }

      // Add the last segment
      ClusteredGammaInfo sub_cluster;
      sub_cluster.lower = std::max( lowest_energy, seg_start );
      sub_cluster.upper = std::min( highest_energy, cluster.upper );
      for( size_t j = gamma_start_idx; j < cluster.gamma_energies.size(); ++j )
      {
        sub_cluster.gamma_energies.push_back( cluster.gamma_energies[j] );
        sub_cluster.gamma_amplitudes.push_back( cluster.gamma_amplitudes[j] );
      }

      if( !sub_cluster.gamma_energies.empty() )
        final_clusters.push_back( std::move( sub_cluster ) );
    }
  }//legacy for( ClusteredGammaInfo &cluster : merged_clusters )
  }

  // Policy-mode finalization is an exact channel transaction, not a best-effort energy clamp.
  // Keep the diagnostic assertion in developer builds and fail closed in Release if a later stage
  // ever violates the reconciler's postcondition.
  for( size_t i = 1; i < final_clusters.size(); ++i )
  {
    const ClusteredGammaInfo &prev = final_clusters[i - 1];
    const ClusteredGammaInfo &curr = final_clusters[i];
    const float prev_upper_inside = std::nextafter(
        static_cast<float>(prev.upper), static_cast<float>(prev.lower) );
    const size_t prev_last_channel = foreground->find_gamma_channel(
        prev_upper_inside );
    const size_t curr_first_channel = foreground->find_gamma_channel(
        static_cast<float>(curr.lower) );
    const bool energy_overlap = (curr.lower < prev.upper);
    const bool channel_overlap = (curr_first_channel <= prev_last_channel);
    if( energy_overlap || (settings.use_automatic_roi_policy && channel_overlap) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      std::cerr << "ERROR: final_clusters[" << (i-1) << "] and [" << i << "] overlap: "
           << "[" << prev.lower << ", " << prev.upper << "] vs "
           << "[" << curr.lower << ", " << curr.upper << "], channels "
           << prev_last_channel << " and " << curr_first_channel << std::endl;
      assert( !energy_overlap );
      assert( !settings.use_automatic_roi_policy || !channel_overlap );
#endif
      if( settings.use_automatic_roi_policy )
        throw std::runtime_error( "cluster_gammas_to_rois produced overlapping final clusters" );
    }
  }

  // Create ROIs from final clusters
  double previous_roi_upper = 0.0;

  for( const ClusteredGammaInfo &cluster : final_clusters )
  {
    RelActCalcAuto::RoiRange roi;

    // Use the pre-calculated bounds from cluster (based on weighted mean and effective FWHM)
    // Constrain lower bound to not overlap with previous ROI, and clamp to valid energy range
    const double unclamped_lower = std::max( lowest_energy, cluster.lower );
#if( PERFORM_DEVELOPER_CHECKS )
    assert( !settings.use_automatic_roi_policy || (unclamped_lower >= previous_roi_upper) );
#endif
    roi.lower_energy = std::max( unclamped_lower, previous_roi_upper );
    roi.upper_energy = std::min( highest_energy, cluster.upper );

    const double mid_energy = 0.5 * (roi.upper_energy + roi.lower_energy);
    const double mid_fwhm = fwhm_at( mid_energy );
    const double num_fwhm_wide = (roi.upper_energy - roi.lower_energy) / mid_fwhm;

    if( num_fwhm_wide < settings.min_fwhm_roi )
    {
      // The cluster passed the keep-gate, so its gammas should not just vanish.  Under-width ROIs
      // arise from clamps, not from the cluster itself (the adaptive-extent core alone is
      // ~2*roi_core_num_fwhm wide): either the previous ROI's upper bound shaved this one's lower
      // bound, or a clean-gap/valley split left a sliver.  When the ROI abuts the previously
      // emitted ROI, fold the cluster back into it (mirroring merge_rois' under-width fold-back
      // guards) instead of silently dropping gammas the keep decision said were significant.
      const double fold_tol = 0.1 * mid_fwhm;
      if( !result_rois.empty() && ((roi.lower_energy - previous_roi_upper) <= fold_tol) )
      {
        std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &previous = result_rois.back();
        previous.first.upper_energy = std::max( previous.first.upper_energy, roi.upper_energy );
        previous.second.upper = std::max( previous.second.upper, cluster.upper );
        previous.second.gamma_energies.insert( std::end(previous.second.gamma_energies),
            std::begin(cluster.gamma_energies), std::end(cluster.gamma_energies) );
        previous.second.gamma_amplitudes.insert( std::end(previous.second.gamma_amplitudes),
            std::begin(cluster.gamma_amplitudes), std::end(cluster.gamma_amplitudes) );
        // Carry the atom ledger too, so a policy-mode fold never loses admitted provenance
        // (a no-op in legacy mode, where clusters carry no atoms).
        previous.second.atoms.insert( std::end(previous.second.atoms),
            std::begin(cluster.atoms), std::end(cluster.atoms) );
        previous_roi_upper = previous.first.upper_energy;
        continue;
      }
      if( settings.use_automatic_roi_policy )
      {
        // The atom-safe partition already enforced the minimum width and channel-disjointness, so a
        // residual sliver here is a rare clamp/rounding artifact.  Emit it (fall through) rather
        // than drop the admitted gammas - dropping a requested line is exactly the regression this
        // work removes.
        if( should_debug_print() )
          std::cerr << "cluster_gammas_to_rois: emitting narrow policy ROI [" << roi.lower_energy
               << ", " << roi.upper_energy << "] keV rather than dropping its atoms" << std::endl;
      }else
      {
        if( should_debug_print() )
        {
          std::cerr << "cluster_gammas_to_rois: Rejected ROI [" << roi.lower_energy << ", " << roi.upper_energy
               << "] keV: too narrow (" << num_fwhm_wide << " FWHM < " << settings.min_fwhm_roi << " min)" << std::endl;
        }
        // Don't update previous_roi_upper since we're not adding this ROI
        continue;
      }
    }

    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;

    // Continuum polynomial order by AICc over the ROI's continuum sidebands (channels outside
    // the peak core) - whether the continuum actually curves is a property of the data, not of
    // the window width.  A fixed width prior still gates the quadratic candidate: very narrow
    // ROIs cannot support curvature.
    if( (num_fwhm_wide >= sm_quad_cont_min_roi_fwhm) && !cluster.gamma_energies.empty() )
    {
      const auto minmax_g = std::minmax_element( std::begin(cluster.gamma_energies),
                                                 std::end(cluster.gamma_energies) );
      const double skew_extra = (settings.skew_type == PeakDef::SkewType::NoSkew)
          ? 0.0 : sm_skew_low_side_extra_fwhm;
      const double core_lo = *minmax_g.first
          - (settings.roi_core_num_fwhm + skew_extra) * fwhm_at( *minmax_g.first );
      const double core_hi = *minmax_g.second
          + settings.roi_core_num_fwhm * fwhm_at( *minmax_g.second );

      roi.continuum_type = detail::select_continuum_order_by_sidebands(
          foreground, roi.lower_energy, roi.upper_energy, core_lo, core_hi,
          settings.cont_order_aicc_penalty );
    }

    // Step-continuum decision: a gated "let the model choose" ladder.
    //   Gate 1 - peak dominance: a step riser is generated by the peak's own scattered photons,
    //     so it only matters for peaks that tower over the continuum (z = S/sqrt(S+B), GA gene).
    //   Gate 2 - loose sideband asymmetry: the low-side continuum must sit above the high side at
    //     >= sm_step_trial_min_asym_z (deliberately permissive - a "worth investigating" hint that
    //     cannot self-veto like the former fixed probe geometry).
    //   Stage 3 - trial fit: fit the ROI with the selected polynomial continuum and with its
    //     equal-parameter-count step variant (linear least squares, not Ceres) and keep the step
    //     only if it beats the polynomial chi2 by settings.step_trial_chi2_margin.
    if( !cluster.gamma_amplitudes.empty() )
    {
      const size_t max_amp_idx = static_cast<size_t>(
          std::max_element( std::begin(cluster.gamma_amplitudes), std::end(cluster.gamma_amplitudes) )
          - std::begin(cluster.gamma_amplitudes) );
      const double max_amplitude = cluster.gamma_amplitudes[max_amp_idx];
      const double dominant_energy = cluster.gamma_energies[max_amp_idx];

      // Significance of the dominant peak over the local continuum: z = S/sqrt(S+B), with B the
      // sideband-estimated continuum over a ~95% Gaussian window (1.665 FWHM) at the dominant
      // gamma.  The sideband estimate has the cluster's own predicted tails subtracted so tall
      // peaks do not inflate their own B.
      const double win_lo = dominant_energy - 0.5*1.665*mid_fwhm;
      const double win_hi = dominant_energy + 0.5*1.665*mid_fwhm;

      const auto step_predicted_signal = [&]( const double x0, const double x1 ) -> double {
        return detail::predicted_gaussian_counts( cluster.gamma_energies, cluster.gamma_amplitudes,
                                                  fwhm_at, x0, x1 );
      };//step_predicted_signal lambda

      const detail::LocalContinuumEstimate step_cont = detail::estimate_local_continuum(
          foreground, roi.lower_energy, roi.upper_energy, mid_fwhm, 0.5,
          step_predicted_signal, unfit_auto_peaks );

      // Gate 2 needs a VALID sideband estimate, so there is no gross-count fallback for Gate 1's
      // B: a ROI at a spectrum edge, or with hopelessly contaminated sidebands, never gets a step.
      if( step_cont.valid )
      {
        // R1 step 2: take Gate-1's dominance background from the single global SNIP continuum when
        // available (falling back to the local step_cont, still gated on step_cont.valid so a step
        // still needs a valid sideband estimate - no gross fallback here).  Gate-2's asymmetry test
        // below stays local: it needs the raw sideband counts SNIP does not provide.
        const double b_est = (settings.global_continuum && settings.global_continuum->valid())
            ? settings.global_continuum->integral( win_lo, win_hi )
            : step_cont.integral( win_lo, win_hi );
        const double est_significance = max_amplitude
            / std::sqrt( std::max( 1.0, max_amplitude + b_est ) );

        if( (est_significance >= settings.step_cont_min_peak_significance)
            && (step_cont.sideband_asymmetry_z() >= sm_step_trial_min_asym_z) )
        {
          roi.continuum_type = trial_step_continuum( foreground, roi.lower_energy, roi.upper_energy,
              cluster.gamma_energies, cluster.gamma_amplitudes, fwhm_at,
              roi.continuum_type, settings.step_trial_chi2_margin );
        }
      }//if( step_cont.valid )
    }//if( !cluster.gamma_amplitudes.empty() )

    result_rois.push_back( std::make_pair( roi, cluster ) );
    previous_roi_upper = roi.upper_energy;
  }//for( const ClusteredGammaInfo &cluster : final_clusters )

  if( should_debug_print() )
  {
    std::cerr << "cluster_gammas_to_rois: Final " << result_rois.size() << " ROIs:" << std::endl;
    for( size_t i = 0; i < result_rois.size(); ++i )
    {
      const RelActCalcAuto::RoiRange &roi = result_rois[i].first;
      const double width = roi.upper_energy - roi.lower_energy;
      const double mid = 0.5 * (roi.lower_energy + roi.upper_energy);
      const double fwhm = fwhm_at( mid );
      const char *cont_str = "Unknown";
      switch( roi.continuum_type )
      {
        case PeakContinuum::OffsetType::Linear:     cont_str = "Linear";     break;
        case PeakContinuum::OffsetType::Quadratic:  cont_str = "Quadratic";  break;
        case PeakContinuum::OffsetType::FlatStep:   cont_str = "FlatStep";   break;
        case PeakContinuum::OffsetType::LinearStep: cont_str = "LinearStep"; break;
        default: break;
      }
      std::cerr << "  [" << i << "] range=[" << roi.lower_energy << ", " << roi.upper_energy << "] keV ("
           << width << " keV, " << (width / fwhm) << " FWHM), cont=" << cont_str
           << ", " << result_rois[i].second.gamma_energies.size() << " gammas" << std::endl;
    }
  }

  // R4 shadow mode: jointly propose boundaries from the pre-merge source groups and shared SNIP,
  // but deliberately leave `result_rois` untouched.  The evaluator drains the diagnostics from
  // this thread after each fit for paired old/proposed review.
  if( !supplied_predicted_gammas
      && (std::getenv("PEAKFIT_ROI_SHADOW_TSV")
          || std::getenv("INTERSPEC_ROI_BOUNDARY_SHADOW")) )
  {
    std::vector<detail::RoiBoundaryShadowGroup> shadow_groups;
    for( const ClusteredGammaInfo &cluster : clustered_gammas )
    {
      for( const double gamma_energy : cluster.gamma_energies )
      {
        const auto legacy = std::find_if( std::begin(result_rois), std::end(result_rois),
          [gamma_energy]( const std::pair<RelActCalcAuto::RoiRange,ClusteredGammaInfo> &entry ) {
            return (gamma_energy >= entry.first.lower_energy)
                && (gamma_energy <= entry.first.upper_energy);
          } );
        if( legacy == std::end(result_rois) )
          continue;
        detail::RoiBoundaryShadowGroup group;
        group.legacy_lower = legacy->first.lower_energy;
        group.legacy_upper = legacy->first.upper_energy;
        group.gamma_energies.push_back( gamma_energy );
        shadow_groups.push_back( std::move(group) );
      }
    }
    std::sort( std::begin(shadow_groups), std::end(shadow_groups),
      []( const detail::RoiBoundaryShadowGroup &lhs,
          const detail::RoiBoundaryShadowGroup &rhs ) {
        if( lhs.gamma_energies.front() != rhs.gamma_energies.front() )
          return lhs.gamma_energies.front() < rhs.gamma_energies.front();
        if( lhs.legacy_lower != rhs.legacy_lower )
          return lhs.legacy_lower < rhs.legacy_lower;
        return lhs.legacy_upper < rhs.legacy_upper;
      } );
    shadow_groups.erase( std::unique( std::begin(shadow_groups), std::end(shadow_groups),
      []( const detail::RoiBoundaryShadowGroup &lhs,
          const detail::RoiBoundaryShadowGroup &rhs ) {
        return (std::fabs(lhs.gamma_energies.front() - rhs.gamma_energies.front()) < 0.01)
            && (std::fabs(lhs.legacy_lower - rhs.legacy_lower) < 0.05)
            && (std::fabs(lhs.legacy_upper - rhs.legacy_upper) < 0.05);
      } ), std::end(shadow_groups) );

    const detail::GlobalContinuumEstimate invalid_global;
    const detail::GlobalContinuumEstimate &shadow_global = shadow_global_override
      ? *shadow_global_override
      : (settings.global_continuum ? *settings.global_continuum : invalid_global);
    detail::RoiBoundaryShadowResult shadow_result = detail::optimize_roi_boundaries_shadow(
        shadow_groups, foreground, shadow_global, fwhm_at, unfit_auto_peaks,
        settings.max_fwhm_width, settings.roi_core_num_fwhm );
    shadow_result.stage = shadow_stage.empty() ? "unspecified clustering" : shadow_stage;
    detail::record_roi_boundary_shadow_result( std::move(shadow_result) );
  }

  // Developer check: Validate result ROIs don't overlap
#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < result_rois.size(); ++i )
  {
    const RelActCalcAuto::RoiRange &prev_roi = result_rois[i - 1].first;
    const RelActCalcAuto::RoiRange &curr_roi = result_rois[i].first;
    if( curr_roi.lower_energy < prev_roi.upper_energy )
    {
      std::cerr << "ERROR: result_rois[" << (i-1) << "] and [" << i << "] overlap: "
           << "[" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] vs "
           << "[" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "]" << std::endl;
      assert( curr_roi.lower_energy >= prev_roi.upper_energy );
    }
  }
#endif

  return result_rois;
}//cluster_gammas_to_rois


struct InitialRoi
{
  RelActCalcAuto::RoiRange roi;
  double center_energy;
  double fwhm;

  /** Estimated peak area, in COUNTS, for the clean-gap tail check in merge_rois (which compares
   predicted tail counts against the Poisson noise of the local continuum - see
   detail::find_clean_gap_between).  Pass 0 when no counts-scale estimate exists (e.g. a bare
   br*efficiency yield): the tail check then degenerates to a pure gap-width test, which is the
   documented "unknown amplitude" behavior rather than a silently-always-passing one.
  */
  double estimated_amplitude = 0.0;  // expected peak area in counts (0 = unknown)
  std::vector<double> modeled_energies;
  std::vector<double> modeled_areas;
  size_t joined_groups = 1;
};

std::vector<RelActCalcAuto::RoiRange> merge_rois(
    std::vector<InitialRoi> initial_rois,
    const PeakFitForNuclideConfig &config,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks = {},
    const std::shared_ptr<const SpecUtils::Measurement> &foreground = {},
    const detail::GlobalContinuumEstimate *global_continuum = nullptr,
    std::vector<AutomaticRoiDecisionDiagnostic> *roi_policy_diagnostics = nullptr,
    const std::string &policy_stage = "initial ROI merge",
    const bool use_automatic_roi_policy = true,
    std::vector<std::pair<double,double>> *merged_modeled_peaks = nullptr )
{
  // Sort by lower_energy for merging
  std::sort( initial_rois.begin(), initial_rois.end(), [](const InitialRoi &a, const InitialRoi &b){
      return a.roi.lower_energy < b.roi.lower_energy;
  } );

  // Step 4: Merge/split overlapping ROIs
  // Critical: ensure NO overlapping ROIs in final result
  std::vector<RelActCalcAuto::RoiRange> merged_rois;
  std::vector<std::vector<double>> merged_centers;  // Track center energies for each merged ROI
  std::vector<std::vector<double>> merged_amplitudes;  // Track amplitude per center for tail check
  std::vector<double> merged_fwhms;  // Track FWHM for validation
  std::vector<size_t> merged_joined_groups;

  const auto current_centers = []( const InitialRoi &roi ) {
    return roi.modeled_energies.empty() ? std::vector<double>{roi.center_energy}
                                        : roi.modeled_energies;
  };
  const auto current_areas = []( const InitialRoi &roi ) {
    std::vector<double> areas = roi.modeled_areas.empty()
      ? std::vector<double>{roi.estimated_amplitude} : roi.modeled_areas;
    areas.resize( roi.modeled_energies.empty() ? size_t(1) : roi.modeled_energies.size(), 0.0 );
    return areas;
  };

  // POLICY MODE: reconcile ALL ROIs at once through the atom-safe partition layer, which cannot
  // drop or duplicate a modeled line (the legacy per-pair split/merge below has pop_back/skip drop
  // paths).  Atoms are minted from each ROI's modeled lines for this call; the reconciler carries
  // them exactly-once.  The legacy loop still runs for R6-enabled fits, or when there is no usable
  // calibration for the channel-aligned partition.
  const bool policy_ok = use_automatic_roi_policy && foreground
      && foreground->energy_calibration() && foreground->energy_calibration()->valid()
      && (foreground->num_gamma_channels() > 0);
  if( policy_ok )
  {
    // FWHM(E) from the nearest ROI center (FWHM varies slowly; the legacy path used a single
    // per-pair mid-FWHM, which this strictly improves on).
    const auto fwhm_at = [&initial_rois]( const double energy ) -> double {
      double best_fwhm = std::numeric_limits<double>::quiet_NaN();
      double best_dist = std::numeric_limits<double>::infinity();
      for( const InitialRoi &ir : initial_rois )
      {
        if( !std::isfinite(ir.fwhm) || (ir.fwhm <= 0.0) )
          continue;
        const double d = std::fabs( energy - ir.center_energy );
        if( d < best_dist ){ best_dist = d; best_fwhm = ir.fwhm; }
      }
      return best_fwhm;
    };//fwhm_at

    std::vector<detail::AutomaticRoiComponent> comps;
    comps.reserve( initial_rois.size() );
    for( const InitialRoi &ir : initial_rois )
    {
      detail::AutomaticRoiComponent c;
      c.lower = ir.roi.lower_energy;
      c.upper = ir.roi.upper_energy;
      c.first_channel = foreground->find_gamma_channel( static_cast<float>(ir.roi.lower_energy) );
      c.last_channel = foreground->find_gamma_channel( static_cast<float>(ir.roi.upper_energy) );
      c.continuum_type = ir.roi.continuum_type;
      c.range_limits_type = ir.roi.range_limits_type;
      c.joined_groups = ir.joined_groups;
      c.protected_geometry = false;
      const std::vector<double> centers = current_centers( ir );
      const std::vector<double> areas = current_areas( ir );
      const bool is_modeled = !ir.modeled_energies.empty();
      for( size_t i = 0; i < centers.size(); ++i )
      {
        detail::RoiAtom atom;
        atom.id = detail::next_roi_atom_id();
        atom.energy = centers[i];
        atom.area = (i < areas.size()) ? areas[i] : 0.0;
        atom.kind = is_modeled ? detail::RoiAtomKind::ModeledGamma
                               : detail::RoiAtomKind::FoundPeakEvidence;
        c.atoms.push_back( atom );
      }
      comps.push_back( std::move(c) );
    }

    // Unfit auto-search peaks that coincide with a modeled line are not interfering features; drop
    // them so the policy's unmodeled-feature test only sees genuine unmodeled peaks (matches the
    // resolve-stage filter).
    std::vector<std::shared_ptr<const PeakDef>> filtered_unfit;
    for( const std::shared_ptr<const PeakDef> &pk : unfit_auto_peaks )
    {
      if( !pk || !pk->gausPeak() )
        continue;
      bool matches_modeled = false;
      for( const detail::AutomaticRoiComponent &c : comps )
      {
        for( const detail::RoiAtom &a : c.atoms )
        {
          const double f = fwhm_at( a.energy );
          if( std::isfinite(f) && (f > 0.0) && (std::fabs(pk->mean() - a.energy) <= 0.75*f) )
          {
            matches_modeled = true;
            break;
          }
        }
        if( matches_modeled )
          break;
      }
      if( !matches_modeled )
        filtered_unfit.push_back( pk );
    }

    detail::AutomaticRoiPolicySettings policy_settings;
    policy_settings.merge_tail_z = config.merge_tail_z;
    policy_settings.merge_clean_gap_fwhm = config.merge_clean_gap_fwhm;
    policy_settings.continuum_aicc_penalty = config.cont_order_aicc_penalty;
    policy_settings.peak_core_num_fwhm = config.auto_roi_core_num_fwhm;
    policy_settings.max_width_fwhm = config.auto_rel_eff_sol_max_fwhm;
    policy_settings.minimum_partition_gap_fwhm = config.auto_roi_partition_min_gap_fwhm;
    policy_settings.allow_clean_gap_partition_override
      = config.auto_roi_partition_allow_clean_gap_override;
    policy_settings.force_partition_gap_fwhm = config.auto_roi_partition_force_gap_fwhm;
    policy_settings.max_partition_children = config.auto_roi_partition_max_children;
    policy_settings.allow_overwide_overlap_partition = config.auto_roi_partition_overwide;
    // This is the final automatic merge of independently constructed source/NORM ranges.  The
    // measured-data partition is opt-in through the configuration gene, so its interaction with
    // the maximum-width and continuum penalties remains experimentally controllable.
    policy_settings.stage = policy_stage;

    detail::AutomaticRoiPartitionConstraints cons;
    cons.lowest_energy = foreground->gamma_channel_lower( 0 );
    cons.highest_energy = foreground->gamma_channel_upper( foreground->num_gamma_channels() - 1 );
    cons.left_barrier = -std::numeric_limits<double>::infinity();
    cons.min_width_fwhm = 1.0;  // matches the legacy "child width >= FWHM" validity check
    cons.peak_core_num_fwhm = config.auto_roi_core_num_fwhm;

    // Route diagnostics to the caller's vector, or to the thread-local sink when it passed none
    // (matching the legacy path, which the whole-fit result later drains).
    std::vector<AutomaticRoiDecisionDiagnostic> local_diags;
    std::vector<AutomaticRoiDecisionDiagnostic> *diag_sink
        = roi_policy_diagnostics ? roi_policy_diagnostics : &local_diags;
    const detail::AutomaticRoiReconcileResult rr = detail::reconcile_automatic_components_one_pass(
        std::move(comps), foreground, global_continuum, fwhm_at, filtered_unfit,
        policy_settings, cons, diag_sink );
    if( !roi_policy_diagnostics )
      for( const AutomaticRoiDecisionDiagnostic &d : local_diags )
        detail::record_automatic_roi_diagnostic( d );

    for( const detail::AutomaticRoiComponent &c : rr.components )
    {
      RelActCalcAuto::RoiRange roi;
      roi.lower_energy = c.lower;
      roi.upper_energy = c.upper;
      roi.continuum_type = c.continuum_type;
      roi.range_limits_type = c.range_limits_type;
      merged_rois.push_back( roi );
      if( merged_modeled_peaks )
        for( const detail::RoiAtom &a : c.atoms )
          merged_modeled_peaks->emplace_back( a.energy, a.area );
    }
    return merged_rois;
  }//if( policy_ok )

  for( size_t roi_idx = 0; roi_idx < initial_rois.size(); ++roi_idx )
  {
    const InitialRoi &current = initial_rois[roi_idx];
    if( merged_rois.empty() )
    {
      merged_rois.push_back( current.roi );
      merged_centers.push_back( current_centers(current) );
      merged_amplitudes.push_back( current_areas(current) );
      merged_fwhms.push_back( current.fwhm );
      merged_joined_groups.push_back( current.joined_groups );
      continue;
    }

    RelActCalcAuto::RoiRange &last = merged_rois.back();
    std::vector<double> &last_centers = merged_centers.back();
    std::vector<double> &last_amps = merged_amplitudes.back();
    const double last_fwhm = merged_fwhms.back();

    // Check if ROIs overlap
    const bool overlaps = (current.roi.lower_energy < last.upper_energy);

    if( !overlaps )
    {
      // No overlap - add new ROI
      merged_rois.push_back( current.roi );
      merged_centers.push_back( current_centers(current) );
      merged_amplitudes.push_back( current_areas(current) );
      merged_fwhms.push_back( current.fwhm );
      merged_joined_groups.push_back( current.joined_groups );
      continue;
    }

    // ROIs overlap - check width constraint
    const double combined_upper = std::max( last.upper_energy, current.roi.upper_energy );
    const double combined_width = combined_upper - last.lower_energy;
    const double mid_energy = 0.5 * (last.lower_energy + combined_upper);
    const double last_mid = 0.5 * (last.lower_energy + last.upper_energy);
    const double current_mid = 0.5 * (current.roi.lower_energy + current.roi.upper_energy);
    const double mid_dist = (fabs(current_mid - last_mid) < 0.1) ? 0.5 : ((mid_energy - last_mid) / (current_mid - last_mid));
    const double mid_fwhm = last_fwhm + mid_dist*(current.fwhm - last_fwhm);
    assert( ((mid_fwhm >= current.fwhm) && (mid_fwhm <= last_fwhm)) || ((mid_fwhm <= current.fwhm) && (mid_fwhm >= last_fwhm)) );
    //const double mid_fwhm = DetectorPeakResponse::peakResolutionFWHM(
    //  static_cast<float>(mid_energy), fwhmFnctnlForm, fwhm_coefficients );

    const bool width_ok = (combined_width <= config.auto_rel_eff_sol_max_fwhm * mid_fwhm);

    // Clean-gap test: keep the ROIs separate only when a continuum-anchoring window exists
    // between their nearest peak centers (see detail::find_clean_gap_between).  Replaces the
    // amplitude-relative tail-fraction test, which ignored counting statistics.  When
    // amplitudes are unknown (0) the test degenerates to a pure gap-width check.
    bool have_clean_gap = false;
    double clean_win_lo = 0.0, clean_win_hi = 0.0;
    if( width_ok && !last_centers.empty() )
    {
      // FWHM varies little between adjacent ROIs; the interpolated mid-point value suffices.
      const auto fwhm_at = [mid_fwhm]( const double ) -> double { return mid_fwhm; };
      const std::vector<double> right_centers = current_centers( current );
      const std::vector<double> right_amps = current_areas( current );

      have_clean_gap = detail::find_clean_gap_between(
          last_centers, last_amps, right_centers, right_amps,
          std::min( last_centers.back(), right_centers.front() ),
          std::max( last_centers.back(), right_centers.front() ),
          foreground, fwhm_at, config.merge_tail_z, config.merge_clean_gap_fwhm,
          &clean_win_lo, &clean_win_hi, global_continuum );

      if( have_clean_gap && should_debug_print() )
      {
        std::cerr << "merge_rois: NOT merging ROIs ["
             << last.lower_energy << ", " << last.upper_energy << "] and ["
             << current.roi.lower_energy << ", " << current.roi.upper_energy
             << "] - clean gap [" << clean_win_lo << ", " << clean_win_hi << "] keV" << std::endl;
      }
    }//if( width_ok && !last_centers.empty() )

    // Check if an unfit peak lies between the center energies of the two ROIs.
    // Skip unfit peaks that match a source gamma (center energy) within clustering tolerance.
    bool unfit_peak_between = false;
    std::vector<std::shared_ptr<const PeakDef>> blocking_unfit_peaks;
    if( !unfit_auto_peaks.empty()
        && (use_automatic_roi_policy || (width_ok && !have_clean_gap)) )
    {
      // Use the last center energy of the merged ROI and the current center energy
      const double last_center = last_centers.back();
      const std::vector<double> right_centers = current_centers( current );
      const double curr_center = right_centers.front();
      const double between_lo = std::min( last_center, curr_center );
      const double between_hi = std::max( last_center, curr_center );

      // Collect all center energies for source-gamma matching
      std::vector<double> all_centers = last_centers;
      all_centers.insert( std::end(all_centers), std::begin(right_centers),
                          std::end(right_centers) );

      for( const std::shared_ptr<const PeakDef> &peak : unfit_auto_peaks )
      {
        const double peak_mean = peak->mean();

        // Check if this unfit peak matches any source gamma (center energy)
        const double peak_sigma = peak->sigma();
        const double tolerance = config.auto_rel_eff_cluster_num_sigma * peak_sigma;
        bool matches_source = false;
        for( const double center : all_centers )
        {
          if( std::fabs( peak_mean - center ) < tolerance )
          {
            matches_source = true;
            break;
          }
        }

        if( matches_source )
          continue;

        // The shared policy owns the +/-FWHM core-intersection test.  Passing only peaks whose
        // means are strictly between the anchors misses a core that crosses into the gap from
        // just outside it.  Keep the historical mean-only Boolean for the legacy path.
        if( use_automatic_roi_policy )
          blocking_unfit_peaks.push_back( peak );
        if( (peak_mean <= between_lo) || (peak_mean >= between_hi) )
          continue;
        if( !use_automatic_roi_policy )
          blocking_unfit_peaks.push_back( peak );
        unfit_peak_between = true;

        if( should_debug_print() )
        {
          std::cerr << "merge_rois: NOT merging ROIs ["
               << last.lower_energy << ", " << last.upper_energy << "] and ["
               << current.roi.lower_energy << ", " << current.roi.upper_energy
               << "] due to unfit peak at " << peak_mean
               << " keV between centers at " << last_center
               << " and " << curr_center << " keV" << std::endl;
        }
      }//for( const auto &peak : unfit_auto_peaks )
    }//if( width_ok && !unfit_auto_peaks.empty() )

    detail::AutomaticRoiPolicyResult policy;
    if( use_automatic_roi_policy )
    {
      const auto fwhm_at = [mid_fwhm]( const double ) -> double { return mid_fwhm; };
      detail::AutomaticRoiGroup left_group;
      left_group.lower = last.lower_energy;
      left_group.upper = last.upper_energy;
      left_group.peak_energies = last_centers;
      left_group.peak_areas = last_amps;
      left_group.joined_groups = merged_joined_groups.back();
      detail::AutomaticRoiGroup right_group;
      right_group.lower = current.roi.lower_energy;
      right_group.upper = current.roi.upper_energy;
      right_group.peak_energies = current_centers( current );
      right_group.peak_areas = current_areas( current );
      right_group.joined_groups = current.joined_groups;
      detail::AutomaticRoiPolicySettings policy_settings;
      policy_settings.merge_tail_z = config.merge_tail_z;
      policy_settings.merge_clean_gap_fwhm = config.merge_clean_gap_fwhm;
      policy_settings.continuum_aicc_penalty = config.cont_order_aicc_penalty;
      policy_settings.peak_core_num_fwhm = config.auto_roi_core_num_fwhm;
      policy_settings.max_width_fwhm = config.auto_rel_eff_sol_max_fwhm;
      policy_settings.allow_overwide_overlap_partition = config.auto_roi_partition_overwide;
      policy_settings.stage = policy_stage;
      policy = detail::evaluate_automatic_roi_boundary(
          left_group, right_group, foreground, global_continuum, fwhm_at,
          blocking_unfit_peaks, policy_settings );
      if( roi_policy_diagnostics )
        roi_policy_diagnostics->push_back( policy.diagnostic );
      else
        detail::record_automatic_roi_diagnostic( policy.diagnostic );
    }
    const bool policy_allows_merge = use_automatic_roi_policy
      ? ((policy.decision == AutomaticRoiDecision::MergeInseparable)
          || (policy.decision == AutomaticRoiDecision::MergeInseparableWide))
      : (width_ok && !have_clean_gap && !unfit_peak_between);

    if( policy_allows_merge )
    {
      // MERGE: Extend last ROI to encompass both
      last.upper_energy = combined_upper;
      const std::vector<double> centers = current_centers( current );
      const std::vector<double> areas = current_areas( current );
      last_centers.insert( std::end(last_centers), std::begin(centers), std::end(centers) );
      last_amps.insert( std::end(last_amps), std::begin(areas), std::end(areas) );
      merged_joined_groups.back() += current.joined_groups;
      // Update FWHM to use the average
      merged_fwhms.back() = 0.5 * (last_fwhm + current.fwhm);

      if( should_debug_print() )
      {
        std::cerr << "Merged overlapping ROI [" << current.roi.lower_energy << ", "
             << current.roi.upper_energy << "] into [" << last.lower_energy
             << ", " << last.upper_energy << "]" << std::endl;
      }
    }
    else
    {
      // SPLIT: find the natural valley between the peaks for the split point.
      // Constrain the split so both ROIs still contain their peak centers.
      const double overlap_lower = current.roi.lower_energy;
      const double overlap_upper = last.upper_energy;

      // The split point must be between the rightmost peak of the left ROI
      // and the center of the right ROI, so both ROIs still contain their peaks;
      // additionally constrained to the clean window when one was found.
      double split_constraint_lower = last_centers.empty()
        ? overlap_lower : std::max( overlap_lower, last_centers.back() );
      const std::vector<double> right_centers = current_centers( current );
      const double right_anchor = *std::min_element(
          std::begin(right_centers), std::end(right_centers) );
      double split_constraint_upper = std::min( overlap_upper, right_anchor );
      if( use_automatic_roi_policy && (policy.boundary_energy > 0.0) )
      {
        split_constraint_lower = std::max( split_constraint_lower, policy.boundary_energy );
        split_constraint_upper = std::min( split_constraint_upper, policy.boundary_energy );
      }
      else if( have_clean_gap )
      {
        split_constraint_lower = std::max( split_constraint_lower, clean_win_lo );
        split_constraint_upper = std::min( split_constraint_upper, clean_win_hi );
      }

      double split_point;
      if( use_automatic_roi_policy
          && (policy.boundary_energy >= split_constraint_lower)
          && (policy.boundary_energy <= split_constraint_upper) )
      {
        split_point = policy.boundary_energy;
      }
      else if( split_constraint_lower >= split_constraint_upper )
      {
        // Constraints don't allow a valid range - fall back to overlap midpoint
        split_point = 0.5 * (overlap_lower + overlap_upper);
      }
      else
      {
        // Search for the spectrum valley between the peak centers
        split_point = find_spectrum_valley( foreground,
            overlap_lower, overlap_upper, mid_fwhm,
            split_constraint_lower, split_constraint_upper );
      }

      // Adjust boundaries at the split point, leaving at least ~1 channel gap
      // so the two ROIs don't share a boundary and the continuum fits are
      // independent.
      double half_gap = 0.5;  // keV, default fallback
      if( foreground && foreground->energy_calibration() && foreground->energy_calibration()->valid() )
      {
        const double ch = foreground->find_gamma_channel( split_point );
        half_gap = 0.5 * foreground->gamma_channel_width( static_cast<size_t>( ch ) );
      }

      const double original_last_upper = last.upper_energy;
      RelActCalcAuto::RoiRange adjusted_current = current.roi;
      if( use_automatic_roi_policy
          && (policy.decision == AutomaticRoiDecision::UnmodeledFeatureBlocked)
          && (policy.exclusion_upper > policy.exclusion_lower) && foreground )
      {
        const size_t lower_channel = foreground->find_gamma_channel(
            static_cast<float>(policy.exclusion_lower) );
        const size_t upper_channel = foreground->find_gamma_channel(
            static_cast<float>(policy.exclusion_upper) );
        last.upper_energy = foreground->gamma_channel_lower( lower_channel );
        adjusted_current.lower_energy = foreground->gamma_channel_upper( upper_channel );
      }else
      {
        last.upper_energy = split_point - half_gap;
        adjusted_current.lower_energy = split_point + half_gap;
      }

      // Validate adjusted last ROI still contains all its center energies and is wide enough
      const double last_width = last.upper_energy - last.lower_energy;
      bool last_valid = (last_width >= last_fwhm);
      for( const double center : last_centers )
      {
        if( center < last.lower_energy || center > last.upper_energy )
        {
          last_valid = false;
          break;
        }
      }

      // Validate adjusted current ROI is still useful
      const double adjusted_width = adjusted_current.upper_energy - adjusted_current.lower_energy;
      const bool current_contains_center = std::all_of(
          std::begin(right_centers), std::end(right_centers),
          [&adjusted_current]( const double center ) {
            return (center >= adjusted_current.lower_energy)
                && (center <= adjusted_current.upper_energy);
          } );
      const bool current_wide_enough = (adjusted_width >= current.fwhm);
      const bool current_valid = current_contains_center && current_wide_enough;

      if( !last_valid )
      {
        // Last ROI is no longer valid after split - remove it and re-process current
        // against the new back of merged_rois (which may also overlap).
        // Using --roi_idx so the for-loop increment brings us back to the same index.
        merged_rois.pop_back();
        merged_centers.pop_back();
        merged_amplitudes.pop_back();
        merged_fwhms.pop_back();
        merged_joined_groups.pop_back();

        if( should_debug_print() )
        {
          std::cerr << "Removed last ROI (invalid after split), will re-process current ROI at "
               << current.center_energy << " keV against prior ROI" << std::endl;
        }

        --roi_idx;
        continue;
      }
      else if( current_valid )
      {
        // Both ROIs valid - add the split current ROI
        merged_rois.push_back( adjusted_current );
        merged_centers.push_back( current_centers(current) );
        merged_amplitudes.push_back( current_areas(current) );
        merged_fwhms.push_back( current.fwhm );
        merged_joined_groups.push_back( current.joined_groups );

        if( should_debug_print() )
        {
          std::cerr << "Split overlapping ROIs: split_point=" << split_point
               << " keV, last=[" << last.lower_energy << ", " << last.upper_energy
               << "], current=[" << adjusted_current.lower_energy << ", "
               << adjusted_current.upper_energy << "]" << std::endl;
        }
      }
      else if( !use_automatic_roi_policy && width_ok )
      {
        // Preserve the legacy R6-enabled behavior exactly: if a split would invalidate the new
        // child but the combined range is under the old width limit, fold it back into the left.
        last.upper_energy = combined_upper;
        const std::vector<double> centers = current_centers( current );
        const std::vector<double> areas = current_areas( current );
        last_centers.insert( std::end(last_centers), std::begin(centers), std::end(centers) );
        last_amps.insert( std::end(last_amps), std::begin(areas), std::end(areas) );
        merged_joined_groups.back() += current.joined_groups;
        merged_fwhms.back() = 0.5 * (last_fwhm + current.fwhm);
      }
      else
      {
        // The policy rejected a join, and the new child cannot retain its modeled core.  Reject
        // the new child instead of silently bypassing the boundary decision.
        last.upper_energy = original_last_upper;

        if( should_debug_print() )
        {
          std::cerr << "Skipping adjusted current ROI at " << current.center_energy << " keV: ";
          if( !current_contains_center )
            std::cerr << "doesn't contain source energy";
          else
            std::cerr << "too narrow (" << adjusted_width << " keV < " << current.fwhm << " keV FWHM)";
          std::cerr << " (automatic ROI policy forbids folding it through the boundary)" << std::endl;
        }
      }
    }
  }

  // Validate no overlaps (developer check)
#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < merged_rois.size(); ++i )
  {
    if( merged_rois[i].lower_energy < merged_rois[i-1].upper_energy )
    {
      std::cerr << "merge_rois OVERLAP at index " << i << " of " << merged_rois.size()
           << ": ROI[" << (i-1) << "]=[" << merged_rois[i-1].lower_energy
           << ", " << merged_rois[i-1].upper_energy
           << "], ROI[" << i << "]=[" << merged_rois[i].lower_energy
           << ", " << merged_rois[i].upper_energy << "]" << std::endl;
      std::cerr << "  Input had " << initial_rois.size() << " ROIs:" << std::endl;
      for( size_t j = 0; j < initial_rois.size(); ++j )
      {
        std::cerr << "    input[" << j << "]: roi=[" << initial_rois[j].roi.lower_energy
             << ", " << initial_rois[j].roi.upper_energy
             << "], center=" << initial_rois[j].center_energy
             << ", fwhm=" << initial_rois[j].fwhm << std::endl;
      }
      std::cerr << "  Output has " << merged_rois.size() << " ROIs:" << std::endl;
      for( size_t j = 0; j < merged_rois.size(); ++j )
      {
        std::cerr << "    output[" << j << "]: [" << merged_rois[j].lower_energy
             << ", " << merged_rois[j].upper_energy << "]" << std::endl;
      }
    }
    assert( merged_rois[i].lower_energy >= merged_rois[i-1].upper_energy );
  }
#endif

  if( should_debug_print() )
  {
    std::cerr << "estimate_initial_rois_without_peaks: Created " << merged_rois.size()
         << " final ROIs" << std::endl;
  }

  if( merged_modeled_peaks )
  {
    for( size_t roi_index = 0; roi_index < merged_centers.size(); ++roi_index )
    {
      for( size_t peak_index = 0; peak_index < merged_centers[roi_index].size(); ++peak_index )
      {
        const double area = (peak_index < merged_amplitudes[roi_index].size())
          ? merged_amplitudes[roi_index][peak_index] : 0.0;
        merged_modeled_peaks->emplace_back( merged_centers[roi_index][peak_index], area );
      }
    }
  }

  return merged_rois;
}//std::vector<InitialRoi> merge_rois(...)


std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_without_peaks(
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitUtils::CoarseResolutionType det_type,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const double min_valid_energy,
  const double max_valid_energy,
  const PeakFitForNuclideConfig &config,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground = {},
  const bool use_automatic_roi_policy = true,
  std::vector<std::pair<double,double>> *modeled_peak_candidates = nullptr )
{
  // Step 1: Get or create valid DRF (use generic if nullptr)
  std::shared_ptr<const DetectorPeakResponse> drf_to_use = drf;
  if( !drf_to_use || !drf_to_use->isValid() )
  {
    switch( det_type )
    {
      case PeakFitUtils::CoarseResolutionType::High:
        drf_to_use = DetectorPeakResponse::getGenericHPGeDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::LaBr:
      case PeakFitUtils::CoarseResolutionType::MedRes:
        drf_to_use = DetectorPeakResponse::getGenericLaBrDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::CZT:
        drf_to_use = DetectorPeakResponse::getGenericCZTGeneralDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::Low:
      case PeakFitUtils::CoarseResolutionType::LowOrMedRes:
      case PeakFitUtils::CoarseResolutionType::Unknown:
      default:
        drf_to_use = DetectorPeakResponse::getGenericNaIDetector();
        break;
    }//switch( det_type )
  }

  if( !drf_to_use || !drf_to_use->isValid() )
    return {};

  // Step 2: Collect top gammas per source
  // Data structure to hold gamma info
  struct GammaInfo
  {
    double energy;
    double br_times_eff;
    RelActCalcAuto::SrcVariant source;  // For debugging
  };

  std::vector<GammaInfo> selected_gammas;

#if( PERFORM_DEVELOPER_CHECKS )
  // Verify no duplicate sources
  for( size_t i = 0; i < sources.size(); ++i )
  {
    for( size_t j = i + 1; j < sources.size(); ++j )
    {
      if( sources[i].source == sources[j].source )
      {
        std::cerr << "estimate_initial_rois_without_peaks: Duplicate source at indices "
             << i << " and " << j << ": " << RelActCalcAuto::to_name( sources[i].source ) << std::endl;
      }
      assert( sources[i].source != sources[j].source );
    }
  }
#endif

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
      continue;

    // Get source age and photons
    const double age = src.age;
    const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src.source, 1.0, age );

    // Compute BR*eff scores for valid gammas
    std::vector<GammaInfo> candidates;
    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      if( photon.energy < min_valid_energy || photon.energy > max_valid_energy )
        continue;

      const double br = photon.numPerSecond;  // BR since we used unit activity
      const double eff = drf_to_use->intrinsicEfficiency( static_cast<float>(photon.energy) );
      const double score = br * eff;

      if( score > 0.0 )
        candidates.push_back( {photon.energy, score, src.source} );
    }

    // Sort by score (descending) and take top 4
    std::sort( candidates.begin(), candidates.end(),
      [](const GammaInfo &a, const GammaInfo &b) { return a.br_times_eff > b.br_times_eff; } );

    const size_t num_to_take = std::min( candidates.size(), size_t(4) );
    for( size_t i = 0; i < num_to_take; ++i )
      selected_gammas.push_back( candidates[i] );

    // Debug output
    if( should_debug_print() )
    {
      std::cerr << "estimate_initial_rois_without_peaks: Source "
           << RelActCalcAuto::to_name( src.source ) << " - selected " << num_to_take << " gammas" << std::endl;
      for( size_t i = 0; i < num_to_take; ++i )
      {
        std::cerr << "  " << candidates[i].energy << " keV, BR*eff=" << candidates[i].br_times_eff << std::endl;
      }
    }
  }

  if( selected_gammas.empty() )
  {
    if( should_debug_print() )
      std::cerr << "estimate_initial_rois_without_peaks: No valid gammas found" << std::endl;
    return {};
  }

  // Step 3: Create initial ROIs
  std::vector<InitialRoi> initial_rois;

  const bool have_fwhm_range = ((lower_fwhm_energy > 0.0)
                                && (upper_fwhm_energy > 0.0)
                                && (lower_fwhm_energy < upper_fwhm_energy));

  for( const GammaInfo &gamma : selected_gammas )
  {
    // Compute FWHM (clamp energy to valid range to avoid extrapolation)
    const double fwhm_eval_energy = have_fwhm_range
        ? std::clamp( gamma.energy, lower_fwhm_energy, upper_fwhm_energy )
        : gamma.energy;
    const double fwhm = DetectorPeakResponse::peakResolutionFWHM(
      static_cast<float>(fwhm_eval_energy), fwhmFnctnlForm, fwhm_coefficients );

    // Validate FWHM
    if( !std::isfinite( fwhm ) || fwhm <= 0.0 )
    {
      if( should_debug_print() )
        std::cerr << "Warning: Invalid FWHM at " << gamma.energy << " keV, skipping" << std::endl;
      continue;
    }

    // Create ROI: energy ± 2.5 FWHM, clamped to valid range
    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = std::max( min_valid_energy, gamma.energy - 2.5 * fwhm );
    roi.upper_energy = std::min( max_valid_energy, gamma.energy + 2.5 * fwhm );
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = use_automatic_roi_policy
      ? RelActCalcAuto::RoiRange::RangeLimitsType::Fixed
      : RelActCalcAuto::RoiRange::RangeLimitsType::CanBeBrokenUp;

    // No activity estimate exists in this no-matched-peaks path, so there is no counts-scale
    // amplitude; pass 0 ("unknown") rather than the dimensionless br*eff yield, which
    // find_clean_gap_between would misinterpret as a (near-zero) count and silently pass the
    // tail-contamination check for every window.
    InitialRoi initial{roi, gamma.energy, fwhm, 0.0};
    initial.modeled_energies.push_back( gamma.energy );
    initial.modeled_areas.push_back( 0.0 );
    initial_rois.push_back( std::move(initial) );
  }

  if( initial_rois.empty() )
  {
    if( should_debug_print() )
      std::cerr << "estimate_initial_rois_without_peaks: No valid ROIs created" << std::endl;
    return {};
  }

  // With no source match, every auto-search core is unmodeled evidence and must remain available
  // to the shared boundary policy; it may exclude a gap but cannot authorize a merge-through.
  const auto policy_fwhm = [=, &fwhm_coefficients]( const double energy ) {
    const double eval_energy = have_fwhm_range
      ? std::clamp(energy, lower_fwhm_energy, upper_fwhm_energy) : energy;
    return static_cast<double>( DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(eval_energy), fwhmFnctnlForm, fwhm_coefficients) );
  };
  detail::GlobalContinuumEstimate initial_continuum;
  if( use_automatic_roi_policy )
  {
    initial_continuum = detail::make_global_continuum(
        foreground, policy_fwhm, det_type, min_valid_energy, max_valid_energy );
  }
  const std::vector<std::shared_ptr<const PeakDef>> no_unfit_peaks;
  const std::vector<std::shared_ptr<const PeakDef>> &merge_unfit_peaks
    = use_automatic_roi_policy ? unfit_auto_peaks : no_unfit_peaks;
  return merge_rois( initial_rois, config, merge_unfit_peaks, foreground,
      initial_continuum.valid() ? &initial_continuum : nullptr, nullptr,
      "no matched peaks", use_automatic_roi_policy, modeled_peak_candidates );
}//estimate_initial_rois_without_peaks

// Forward declarations (defined below) - used by the shielding-robust fallback rel-eff curve.
std::shared_ptr<const DetectorPeakResponse> generic_drf_for_rel_eff_extrap(
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitUtils::CoarseResolutionType det_type );
double shape_rel_eff_above_boundary( const double re_hi, const double energy, const double e_hi,
                                     const std::shared_ptr<const DetectorPeakResponse> &drf );

std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_fallback(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::vector<RelActCalcManual::GenericPeakInfo> &matched_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitUtils::CoarseResolutionType det_type,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const GammaClusteringSettings &settings,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
  std::vector<std::pair<double,double>> *modeled_peak_candidates = nullptr )
{
  // Step 1: Get or create valid DRF
  std::shared_ptr<const DetectorPeakResponse> drf_to_use = drf;
  if( !drf_to_use || !drf_to_use->isValid() )
  {
    switch( det_type )
    {
      case PeakFitUtils::CoarseResolutionType::High:
        drf_to_use = DetectorPeakResponse::getGenericHPGeDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::LaBr:
      case PeakFitUtils::CoarseResolutionType::MedRes:
        drf_to_use = DetectorPeakResponse::getGenericLaBrDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::CZT:
        drf_to_use = DetectorPeakResponse::getGenericCZTGeneralDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::Low:
      case PeakFitUtils::CoarseResolutionType::LowOrMedRes:
      case PeakFitUtils::CoarseResolutionType::Unknown:
      default:
        drf_to_use = DetectorPeakResponse::getGenericNaIDetector();
        break;
    }//switch( det_type )
  }

  if( !drf_to_use || !drf_to_use->isValid() )
    return {};

  const double live_time = foreground->live_time();
  if( live_time <= 0.0 )
    return {};

  // Find valid energy range, clamped to a physically-valid low-energy floor (see low_energy_analysis_floor).
  const std::pair<double,double> raw_valid_range = find_valid_energy_range( foreground );
  const double low_e_floor = low_energy_analysis_floor( drf, det_type );
  const double min_valid_energy = (low_e_floor < raw_valid_range.second)
                                  ? std::max( raw_valid_range.first, low_e_floor ) : raw_valid_range.first;
  const double max_valid_energy = raw_valid_range.second;

  // Step 2: Build a shielding-robust relative-efficiency seed.
  //
  // The original fallback calibrated each source's activity from its single highest-branching-ratio
  // gamma and predicted every other line with the generic (un-shielded) DRF efficiency.  For a
  // shielded source that inverts reality: it predicts the strongly-attenuated low-energy lines as
  // dominant and drops the high-energy lines that actually carry the spectrum (e.g. shielded Ir-192,
  // where 588-885 keV dominate but 296-468 keV are attenuated).  Instead, derive the efficiency
  // *shape* from the OBSERVED matched-peak areas (area / yield) - that is measured, so it captures
  // whatever shielding is present.  Only sources with no matched peak fall back to the generic-DRF
  // brightest-gamma estimate.
  const auto eff_drf = [&drf_to_use]( double energy ) -> double {
    return std::max( 1.0e-12, static_cast<double>( drf_to_use->intrinsicEfficiency( static_cast<float>(energy) ) ) );
  };

  std::vector<tuple<RelActCalcAuto::SrcVariant, double, double>> source_age_and_acts;

  // Empirical (energy, rel_eff) points pooled across sources.  Each source's points are normalised by
  // that source's scalar activity, so they share a common scale and trace out the true (shielded)
  // efficiency curve.
  std::vector<std::pair<double,double>> empirical_pts;

  const bool have_fwhm_range = ((lower_fwhm_energy > 0.0)
                                && (upper_fwhm_energy > 0.0)
                                && (lower_fwhm_energy < upper_fwhm_energy));

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
      continue;

    const double age = src.age;
    const std::string src_name = RelActCalcAuto::to_name( src.source );

    // This source's matched peaks: (energy, observed area, dominant yield for this source).
    // matched_peaks is the contaminant-pruned set from the manual stage (peaks the source's OWN fitted
    // curve could account for), so no extra NORM/contaminant filtering is done here - filtering by
    // energy (e.g. is_near_strong_norm_gamma) would wrongly drop a source's own lines when the source
    // itself is a NORM nuclide, and would over-reach on low-resolution detectors.
    std::vector<std::tuple<double,double,double>> matched_for_src;
    for( const RelActCalcManual::GenericPeakInfo &peak : matched_peaks )
    {
      double yield = 0.0;
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
        if( (line.m_isotope == src_name) && (line.m_yield > yield) )
          yield = line.m_yield;
      if( (yield > 0.0) && (peak.m_counts > 0.0) )
        matched_for_src.emplace_back( peak.m_energy, peak.m_counts, yield );
    }

    if( !matched_for_src.empty() )
    {
      // Scalar activity from the generic DRF.  It cancels at matched energies (where the empirical
      // curve forces predicted == observed area) and only sets the absolute scale for this source's
      // UN-matched lines, so an un-shielded DRF here is acceptable.
      std::vector<double> acts;
      for( const std::tuple<double,double,double> &t : matched_for_src )
        acts.push_back( std::get<1>(t) / (std::get<2>(t) * eff_drf( std::get<0>(t) )) );
      std::sort( begin(acts), end(acts) );
      const double activity = acts[acts.size()/2];  // median

      if( activity > 0.0 )
      {
        for( const std::tuple<double,double,double> &t : matched_for_src )
          empirical_pts.emplace_back( std::get<0>(t), std::get<1>(t) / (activity * std::get<2>(t)) );
        source_age_and_acts.emplace_back( src.source, age, activity );

        if( should_debug_print() )
          std::cout << "Fallback(empirical): " << src_name << " from " << matched_for_src.size()
                    << " matched peaks, activity=" << activity << std::endl;
        continue;
      }
    }//if( this source has matched peaks )

    // No matched peak for this source: estimate from the single brightest gamma using the generic
    // (un-shielded) DRF, matching an auto-search peak if one is near, else integrating (original
    // behaviour, retained for the rare no-match source).
    const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src.source, 1.0, age );
    double best_yield = 0.0, best_energy = 0.0, best_br = 0.0, best_eff = 0.0;
    for( const SandiaDecay::EnergyRatePair &photon : photons )
    {
      if( photon.energy < min_valid_energy || photon.energy > max_valid_energy )
        continue;
      const double br = photon.numPerSecond;
      const double eff = eff_drf( photon.energy );
      if( (br*eff) > best_yield ){ best_yield = br*eff; best_energy = photon.energy; best_br = br; best_eff = eff; }
    }
    if( best_yield <= 0.0 || best_br <= 0.0 || best_eff <= 0.0 )
      continue;

    const double fwhm_eval_energy = have_fwhm_range
        ? std::clamp( best_energy, lower_fwhm_energy, upper_fwhm_energy ) : best_energy;
    const double fwhm_at_energy = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(fwhm_eval_energy), fwhmFnctnlForm, fwhm_coefficients );

    std::shared_ptr<const PeakDef> matched_peak = nullptr;
    double min_distance = 0.5 * fwhm_at_energy;
    for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
    {
      if( !peak || !peak->gausPeak() )
        continue;
      const double distance = std::abs( peak->mean() - best_energy );
      if( distance < min_distance ){ min_distance = distance; matched_peak = peak; }
    }

    // Do NOT divide by live_time (matches cluster_gammas_to_rois' activity = counts/(br*eff) convention).
    double estimated_activity = 0.0;
    if( matched_peak )
    {
      estimated_activity = matched_peak->peakArea() / (best_br * best_eff);
    }else
    {
      const double total_counts = foreground->gamma_integral( static_cast<float>(best_energy - 0.5*fwhm_at_energy),
                                                              static_cast<float>(best_energy + 0.5*fwhm_at_energy) );
      estimated_activity = (total_counts / 4.0) / (best_br * best_eff);
    }

    if( estimated_activity > 0.0 )
    {
      source_age_and_acts.emplace_back( src.source, age, estimated_activity );
      if( should_debug_print() )
        std::cout << "Fallback(generic-DRF): " << src_name << " no matched peak; brightest gamma "
                  << best_energy << " keV, est activity=" << estimated_activity << std::endl;
    }
  }//for( loop over sources )

  if( source_age_and_acts.empty() )
  {
    if( should_debug_print() )
      std::cerr << "Fallback: Could not estimate activity for any source" << std::endl;
    return {};
  }

  // Step 3: Relative-efficiency function.  Prefer the empirical (measured, shielding-aware) curve;
  // fall back to the generic DRF only if fewer than two empirical points are available.
  std::sort( begin(empirical_pts), end(empirical_pts),
             []( const std::pair<double,double> &a, const std::pair<double,double> &b ){ return a.first < b.first; } );
  const std::shared_ptr<const DetectorPeakResponse> rel_eff_extrap_drf
    = generic_drf_for_rel_eff_extrap( drf, det_type );

  std::function<double(double)> fallback_rel_eff;
  if( empirical_pts.size() >= 2 )
  {
    fallback_rel_eff = [empirical_pts, rel_eff_extrap_drf]( double energy ) -> double {
      // Below the lowest measured point: flat-clamp (downward extrapolation is unreliable; for a
      // shielded source the true efficiency keeps dropping, so flat is a conservative over-estimate
      // that errs toward keeping a low-energy ROI rather than spuriously dropping one).
      if( energy <= empirical_pts.front().first )
        return empirical_pts.front().second;
      // Above the highest measured point: DRF-shaped extrapolation anchored to the boundary value.
      if( energy >= empirical_pts.back().first )
        return shape_rel_eff_above_boundary( empirical_pts.back().second, energy,
                                             empirical_pts.back().first, rel_eff_extrap_drf );
      // Between points: log-linear interpolation.
      size_t i = 1;
      while( (i + 1 < empirical_pts.size()) && (empirical_pts[i].first < energy) )
        ++i;
      const double e0 = empirical_pts[i-1].first, e1 = empirical_pts[i].first;
      const double r0 = empirical_pts[i-1].second, r1 = empirical_pts[i].second;
      const double t = (e1 > e0) ? ((energy - e0) / (e1 - e0)) : 0.0;
      if( (r0 > 0.0) && (r1 > 0.0) )
        return r0 * std::pow( r1 / r0, t );
      return r0 + t * (r1 - r0);
    };
  }else
  {
    fallback_rel_eff = [drf_to_use]( double energy ) -> double {
      return drf_to_use->intrinsicEfficiency( static_cast<float>(energy) );
    };
  }

  // Step 4: Call cluster_gammas_to_rois with estimated activities
  const auto policy_fwhm = [=, &fwhm_coefficients]( const double energy ) {
    const double eval_energy = have_fwhm_range
      ? std::clamp(energy, lower_fwhm_energy, upper_fwhm_energy) : energy;
    return static_cast<double>( DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(eval_energy), fwhmFnctnlForm, fwhm_coefficients) );
  };
  detail::GlobalContinuumEstimate initial_continuum;
  if( settings.use_automatic_roi_policy )
  {
    initial_continuum = detail::make_global_continuum(
        foreground, policy_fwhm, det_type, min_valid_energy, max_valid_energy );
  }
  const detail::GlobalContinuumEstimate * const initial_continuum_ptr
    = initial_continuum.valid() ? &initial_continuum : nullptr;
  GammaClusteringSettings policy_settings = settings;
  if( settings.use_automatic_roi_policy )
    policy_settings.global_continuum = initial_continuum_ptr;
  const std::vector<std::shared_ptr<const PeakDef>> no_unfit_peaks;
  const std::vector<std::shared_ptr<const PeakDef>> &clustering_unfit_peaks
    = settings.use_automatic_roi_policy ? unfit_auto_peaks : no_unfit_peaks;
  const std::vector<std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo>> rois_and_gammas
    = cluster_gammas_to_rois(
      {fallback_rel_eff},
      {source_age_and_acts},
      foreground,
      fwhmFnctnlForm,
      fwhm_coefficients,
      lower_fwhm_energy,
      upper_fwhm_energy,
      min_valid_energy,
      max_valid_energy,
      policy_settings, clustering_unfit_peaks, nullptr, nullptr, nullptr, nullptr,
      "initial fallback", initial_continuum_ptr
    );

  std::vector<RelActCalcAuto::RoiRange> result;
  result.reserve( rois_and_gammas.size() );
  for( const std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &p : rois_and_gammas )
  {
    result.push_back( p.first );
    if( modeled_peak_candidates )
    {
      for( size_t peak_index = 0; peak_index < p.second.gamma_energies.size(); ++peak_index )
      {
        const double area = (peak_index < p.second.gamma_amplitudes.size())
          ? p.second.gamma_amplitudes[peak_index] : 0.0;
        modeled_peak_candidates->emplace_back( p.second.gamma_energies[peak_index], area );
      }
    }
  }

  return result;
}//estimate_initial_rois_fallback


// Resolve a non-null, valid detector response for SHAPING rel-eff extrapolation above the fitted
// span.  Falls back to a generic detector for the resolution class when `drf` is null/invalid;
// returns null only if even the generic detector is unavailable (caller then holds flat).
std::shared_ptr<const DetectorPeakResponse> generic_drf_for_rel_eff_extrap(
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitUtils::CoarseResolutionType det_type )
{
  if( drf && drf->isValid() )
    return drf;

  std::shared_ptr<const DetectorPeakResponse> generic;
  switch( det_type )
  {
    case PeakFitUtils::CoarseResolutionType::High:
      generic = DetectorPeakResponse::getGenericHPGeDetector();
      break;
    case PeakFitUtils::CoarseResolutionType::LaBr:
    case PeakFitUtils::CoarseResolutionType::MedRes:
      generic = DetectorPeakResponse::getGenericLaBrDetector();
      break;
    case PeakFitUtils::CoarseResolutionType::CZT:
      generic = DetectorPeakResponse::getGenericCZTGeneralDetector();
      break;
    case PeakFitUtils::CoarseResolutionType::Low:
    case PeakFitUtils::CoarseResolutionType::LowOrMedRes:
    case PeakFitUtils::CoarseResolutionType::Unknown:
    default:
      generic = DetectorPeakResponse::getGenericNaIDetector();
      break;
  }//switch( det_type )

  if( generic && !generic->isValid() )
    generic.reset();
  return generic;
}//generic_drf_for_rel_eff_extrap


// Shape a boundary rel-eff value to a higher energy using the DRF intrinsic-efficiency falloff:
//   rel_eff(E) = rel_eff(E_hi) * drf_eff(E)/drf_eff(E_hi),  for E > E_hi.
// DetectorPeakResponse::intrinsicEfficiency() self-clamps at the DRF's own energy-range edges, so
// the multiplier is bounded (no blow-up).  Degrades to a flat hold (returns re_hi) when no usable
// DRF is available or the boundary efficiency is non-positive.  Used for UPPER-side extrapolation
// only; the lower side keeps a flat clamp (extrapolating rel-eff downward is unreliable too, and
// the steep low-energy DRF rise would over-predict efficiency below the fitted span).
double shape_rel_eff_above_boundary( const double re_hi, const double energy, const double e_hi,
                                     const std::shared_ptr<const DetectorPeakResponse> &drf )
{
  if( !drf )
    return re_hi;
  const double eff_hi = drf->intrinsicEfficiency( static_cast<float>(e_hi) );
  if( eff_hi <= 0.0 )
    return re_hi;
  const double eff_e = drf->intrinsicEfficiency( static_cast<float>(energy) );
  return re_hi * (eff_e / eff_hi);
}//shape_rel_eff_above_boundary


/** Guarantee an ROI for every data-confirmed source line the auto-search already found.

 For each (energy, FWHM) in `found_energy_fwhm` (auto-search peaks matched to a requested source)
 whose energy is not already inside an ROI in `rois`, append a TIGHT ROI of +/- half_num_fwhm x FWHM.
 These peaks are directly confirmed in the data, so they intentionally bypass the predicted-signal
 keep-gate in cluster_gammas_to_rois - that gate uses z = S_pred/sqrt(S_pred + B) with B integrated
 over a wide core, which deflates on wide-FWHM low-count NaI and drops strong lines the search found
 (e.g. Co58 810, Zr89 909).  Seeding tight (not the adaptive ~4-FWHM extent) also keeps the downstream
 whole-ROI significance test from being diluted over a wide window.

 Self-limiting by construction: the auto-search surfaces only a handful of real peaks even on busy
 NORM chains (measured: 3 on a 300 s NaI Th232), so this cannot flood ROIs the way lowering the
 predicted keep-gate threshold does.  Returns the number of ROIs seeded.  [architecture review 2026-07-18]
 */
size_t seed_tight_rois_for_found_peaks(
    std::vector<RelActCalcAuto::RoiRange> &rois,
    const std::vector<std::pair<double,double>> &found_energy_fwhm,
    const double half_num_fwhm,
    const double lowest_energy,
    const double highest_energy,
    const bool use_automatic_roi_policy = true )
{
  size_t num_seeded = 0;
  for( const std::pair<double,double> &ef : found_energy_fwhm )
  {
    const double energy = ef.first;
    const double fwhm = ef.second;
    if( !std::isfinite(energy) || !std::isfinite(fwhm) || (fwhm <= 0.0)
        || (energy < lowest_energy) || (energy > highest_energy) )
      continue;

    bool covered = false;
    for( const RelActCalcAuto::RoiRange &r : rois )
      covered = covered || ((energy >= r.lower_energy) && (energy <= r.upper_energy));
    if( covered )
      continue;

    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = std::max( lowest_energy, energy - half_num_fwhm*fwhm );
    roi.upper_energy = std::min( highest_energy, energy + half_num_fwhm*fwhm );
    if( !(roi.lower_energy < roi.upper_energy) )
      continue;
    // A tight single-line window: a linear continuum is the honest null (a quadratic can start to
    // absorb the peak); the downstream continuum-order selection may still upgrade it after merging.
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    // FitPeaksForNuclides owns this automatic boundary; keeping it Fixed prevents RelActAuto's
    // generic significant-range recombination from silently bypassing the policy.
    roi.range_limits_type = use_automatic_roi_policy
      ? RelActCalcAuto::RoiRange::RangeLimitsType::Fixed
      : RelActCalcAuto::RoiRange::RangeLimitsType::CanBeBrokenUp;
    rois.push_back( roi );
    ++num_seeded;
  }//for( const std::pair<double,double> &ef : found_energy_fwhm )

  return num_seeded;
}//seed_tight_rois_for_found_peaks


std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_using_relactmanual(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::shared_ptr<const SpecUtils::Measurement> &background,
  const std::vector<std::shared_ptr<const PeakDef>> &background_auto_search_peaks,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitUtils::CoarseResolutionType det_type,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const double min_valid_energy,
  const double max_valid_energy,
  const GammaClusteringSettings &manual_settings,
  const PeakFitForNuclideConfig &config,
  const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
  std::string &fallback_warning,
  std::vector<std::pair<double,double>> *modeled_peak_candidates = nullptr,
  std::vector<RelActCalcManual::GenericPeakInfo> *source_anchor_candidates = nullptr,
  std::vector<RelActCalcAuto::RoiRange> *clean_source_rois = nullptr,
  bool *has_provisional_fallback_source_anchors = nullptr )
{
  std::vector<RelActCalcAuto::RoiRange> initial_rois;
  // A manual rel-eff fit can fail before it has established which matched peaks its fitted curve
  // accounts for.  Keep separately-qualified matches for the later, transactional source-clean
  // challenger; they must never change the generic fallback incumbent itself.
  std::vector<RelActCalcManual::GenericPeakInfo> provisional_fallback_source_anchors;
  if( has_provisional_fallback_source_anchors )
    *has_provisional_fallback_source_anchors = false;

  // Step 1: Convert auto_search_peaks to RelActManual format
  std::vector<RelActCalcManual::GenericPeakInfo> rel_act_manual_peaks;
  const double foreground_live_time = foreground ? foreground->live_time() : 0.0;
  const double background_live_time = background ? background->live_time() : 0.0;
  const double background_scale = (foreground_live_time > 0.0) && (background_live_time > 0.0)
    ? (foreground_live_time / background_live_time) : 0.0;

  for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
  {
    assert( peak && peak->gausPeak() );
    if( !peak || !peak->gausPeak() )
      continue;

    // RelActCalcManual hard-rejects the ENTIRE fit if any peak has counts <= 0 ("peak counts must be
    // >0.0").  The auto-search can occasionally return a non-positive (or non-finite) amplitude - e.g.
    // a peak fit on a steeply-falling continuum - so skip such peaks here rather than let a single one
    // abort the whole manual rel-eff seed and force the (degraded) fallback.  Observed for shielded
    // Eu-152 on NaI: one matched line had amplitude <= 0, so no peaks were fit at all.
    const double amp = peak->amplitude();
    if( !(amp > 0.0) || !std::isfinite(amp) )
    {
      if( should_debug_print() )
        std::cout << "Skipping auto-search peak at " << peak->mean()
                  << " keV with non-positive/non-finite amplitude (" << amp << ")" << std::endl;
      continue;
    }

    RelActCalcManual::GenericPeakInfo peak_info;
    peak_info.m_energy = peak_info.m_mean = peak->mean();
    peak_info.m_fwhm = peak->fwhm();
    peak_info.m_counts = amp;
    // RelActCalcManual rejects a zero counts-uncert unless m_base_rel_eff_uncert == -1, and the HPGe
    // default base uncert is 0, so a peak whose amplitudeUncert() was never computed (returns <= 0)
    // would hard-fail the entire manual rel-eff seed.  Supply a Poisson counting-stats estimate.
    const double amp_uncert = peak->amplitudeUncert();
    peak_info.m_counts_uncert = ( amp_uncert > 0.0 ) ? amp_uncert : std::sqrt( std::max( 1.0, amp ) );
    peak_info.m_base_rel_eff_uncert = config.rel_eff_manual_base_rel_eff_uncert;

    if( background_scale > 0.0 && !background_auto_search_peaks.empty() )
    {
      const PeakDef *matched_background_peak = nullptr;
      double smallest_delta = std::numeric_limits<double>::infinity();
      for( const std::shared_ptr<const PeakDef> &background_peak : background_auto_search_peaks )
      {
        if( !background_peak || !background_peak->gausPeak()
            || !(background_peak->amplitude() > 0.0) )
          continue;
        const double delta = std::fabs( peak->mean() - background_peak->mean() );
        const double match_width = 0.75 * 0.5 * (peak->fwhm() + background_peak->fwhm());
        if( (delta < match_width) && (delta < smallest_delta) )
        {
          matched_background_peak = background_peak.get();
          smallest_delta = delta;
        }
      }//for( background_auto_search_peaks )

      if( matched_background_peak )
      {
        const double background_amp = matched_background_peak->amplitude();
        const double background_uncert = matched_background_peak->amplitudeUncert();
        const double safe_background_uncert = (background_uncert > 0.0)
          ? background_uncert : std::sqrt( std::max( 1.0, background_amp ) );
        const double net_counts = amp - background_scale*background_amp;
        const double net_uncert = std::sqrt( peak_info.m_counts_uncert*peak_info.m_counts_uncert
                                              + background_scale*background_scale
                                                *safe_background_uncert*safe_background_uncert );
        const double net_significance = net_counts / net_uncert;
        if( !(net_counts > 0.0) || !std::isfinite(net_counts)
            || !(net_significance >= 2.25) || !std::isfinite(net_significance) )
        {
          if( should_debug_print() )
            std::cout << "Skipping background-net manual point at " << peak_info.m_energy
                      << " keV: net=" << net_counts << " +/- " << net_uncert
                      << " (z=" << net_significance << ")" << std::endl;
          continue;
        }
        peak_info.m_counts = net_counts;
        peak_info.m_counts_uncert = net_uncert;
        if( should_debug_print() )
          std::cout << "Using background-net manual point at " << peak_info.m_energy
                    << " keV: gross=" << amp << ", scaled background="
                    << background_scale*background_amp << ", net=" << net_counts
                    << " +/- " << net_uncert << std::endl;
      }
    }

    rel_act_manual_peaks.push_back( peak_info );
  }

  // Step 2: Match peaks to source nuclides
  std::vector<RelActCalcManual::PeakCsvInput::NucAndAge> rel_act_manual_srcs;
  for( const RelActCalcAuto::NucInputInfo &src : sources )
    rel_act_manual_srcs.emplace_back( RelActCalcAuto::to_name( src.source ), -1.0, false );

  const std::vector<std::pair<float,float>> energy_ranges{};
  const std::vector<float> excluded_peak_energies{};
  const float real_time = foreground->real_time();

  const RelActCalcManual::PeakCsvInput::NucMatchResults peak_match_results
   = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( rel_act_manual_peaks,
                                                          RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay,
                                                          energy_ranges,
                                                          rel_act_manual_srcs,
                                                          config.initial_nuc_match_cluster_num_sigma, excluded_peak_energies,
                                                          real_time );

  std::vector<RelActCalcManual::GenericPeakInfo> peaks_matched = peak_match_results.peaks_matched;
  std::sort( std::begin(peaks_matched), std::end(peaks_matched),
    []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ){
      return lhs.m_energy < rhs.m_energy;
  });

  if( should_debug_print() )
  {
    const std::vector<std::string> &used_isotopes = peak_match_results.used_isotopes;
    const std::vector<std::string> &unused_isotopes = peak_match_results.unused_isotopes;
    if( unused_isotopes.empty() )
      std::cout << "Matched up all source nuclides to initial peak fit" << std::endl;
    std::cout << "Failed to match up nuclides: {";
    for( const std::string &nuc : unused_isotopes )
      std::cout << nuc << ", ";
    std::cout << "} to initial auto-fit peaks" << std::endl;
    if( !used_isotopes.empty() )
    {
      std::cout << "Matched up nuclides: {";
      for( const std::string &nuc : used_isotopes )
        std::cout << nuc << ", ";
      std::cout << "} to initial auto-fit peaks to a total of " << peaks_matched.size() << " peaks." << std::endl;

      if( !peaks_matched.empty() )
      {
        std::cout << "Matched peak energies (keV): {";
        for( size_t i = 0; i < peaks_matched.size(); ++i )
        {
          std::cout << peaks_matched[i].m_energy;
          if( i < peaks_matched.size() - 1 )
            std::cout << ", ";
        }
        std::cout << "}" << std::endl;
      }
    }
  }

  // If no matched peaks, fall back to estimate_initial_rois_without_peaks
  if( peaks_matched.empty() )
  {
    return estimate_initial_rois_without_peaks(
      sources, drf, det_type,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      min_valid_energy, max_valid_energy, config, unfit_auto_peaks, foreground,
      manual_settings.use_automatic_roi_policy, modeled_peak_candidates );
  }

  // Step 3: The rel-eff equation form/order is chosen per spectrum by small-sample-corrected AIC
  // (AICc) over a bounded candidate ladder (see the solve block below), instead of a fixed GA-tuned
  // form/order per matched-peak count.  Here we just bound the ladder by what the data supports:
  // RelActCalcManual requires (num_fit_activities + eqn_order) <= num_peaks.
  const size_t num_distinct_nuclides = std::max<size_t>( 1, peak_match_results.used_isotopes.size() );
  const size_t max_eqn_order = std::min<size_t>( 4,
      (peaks_matched.size() > num_distinct_nuclides) ? (peaks_matched.size() - num_distinct_nuclides)
                                                     : size_t(0) );

  // Detector for the physical-model candidate: the supplied DRF when available, else the
  // per-class generic (same fallback the old physical retry used).
  std::shared_ptr<const DetectorPeakResponse> phys_model_det = drf;
  if( !phys_model_det )
  {
    switch( det_type )
    {
      case PeakFitUtils::CoarseResolutionType::High:
        phys_model_det = DetectorPeakResponse::getGenericHPGeDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::LaBr:
      case PeakFitUtils::CoarseResolutionType::MedRes:
        phys_model_det = DetectorPeakResponse::getGenericLaBrDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::CZT:
        phys_model_det = DetectorPeakResponse::getGenericCZTGeneralDetector();
        break;
      case PeakFitUtils::CoarseResolutionType::Low:
      case PeakFitUtils::CoarseResolutionType::LowOrMedRes:
      case PeakFitUtils::CoarseResolutionType::Unknown:
      default:
        phys_model_det = DetectorPeakResponse::getGenericNaIDetector();
        break;
    }//switch( det_type )
  }//if( !phys_model_det )

  // Judge a manual rel-eff solution using the systematic uncertainty the curve was actually fit with,
  // rather than RelActCalcManual's reported m_chi2 (which is statistical-only - it explicitly ignores
  // m_base_rel_eff_uncert, see RelActCalcManual.cpp).  For strong/many-count peaks the statistical bars
  // are tiny, so the stat-only chi2/dof is astronomically large and the accept gates below would ALWAYS
  // reject a perfectly good curve and force the (less-informed) estimate_initial_rois_fallback path.
  // Add the configured base rel-eff uncertainty (or a 10% floor when it is zero, as for the HPGe default)
  // in quadrature onto each peak.  This mirrors the internal chi2 loop in solve_relative_efficiency but
  // with the inflated denominator.
  const double judge_sys_frac = (config.rel_eff_manual_base_rel_eff_uncert > 1.0e-6)
                                  ? config.rel_eff_manual_base_rel_eff_uncert : 0.10;
  // A matched peak whose PREDICTED source counts are below this fraction of the OBSERVED counts is a
  // peak the source(s) demonstrably cannot account for - in a NORM background (e.g. trinitite's U/Th
  // chain) the matcher attaches strong background lines (Pb214 295/352, Pb212 238, 511 annih., Bi212
  // 727, Ac228 969, ...) to faint source gammas with ~0 branching ratio.  Such peaks are uninformative
  // about the rel-eff curve (curve value x ~0 yield) yet dominate a stat-only chi2/dof and force an
  // unnecessary fallback.  Exclude them from the acceptance judgment (this is the robust generalization
  // of "exclude peaks in the background / matching NORM").
  const double judge_min_source_frac = 0.25;
  const auto source_accounts_for_peak = [judge_min_source_frac](
      const double expected, const RelActCalcManual::GenericPeakInfo &peak ) -> bool {
    const double boundary = judge_min_source_frac * peak.m_counts;
    const double boundary_uncert = judge_min_source_frac * peak.m_counts_uncert;
    return expected >= (boundary - boundary_uncert);
  };

  // Raw judged chi2 plus the included/excluded peak counts - the AICc model selection below needs
  // the raw terms, while acceptance gates use the chi2/dof wrapper.
  struct JudgedChi2Terms
  {
    double chi2 = std::numeric_limits<double>::max();  // max() means "not judgeable" (failed fit)
    int num_included = 0;
    int num_excluded = 0;
  };//struct JudgedChi2Terms

  const auto judgment_chi2_terms = [judge_sys_frac, &source_accounts_for_peak](
      const RelActCalcManual::RelEffSolution &sol ) -> JudgedChi2Terms
  {
    JudgedChi2Terms terms;

    // Use the solution's own internally-computed predicted counts (m_predicted_peak_counts): the
    // activity/efficiency normalization is not consistent across equation forms, so re-deriving the
    // prediction externally from relative_activity()*relative_efficiency() is unreliable.
    if( (sol.m_status != RelActCalcManual::ManualSolutionStatus::Success)
       || (sol.m_predicted_peak_counts.size() != sol.m_input.peaks.size()) )
      return terms;

    terms.chi2 = 0.0;
    for( size_t i = 0; i < sol.m_input.peaks.size(); ++i )
    {
      const RelActCalcManual::GenericPeakInfo &peak = sol.m_input.peaks[i];
      const double expected = sol.m_predicted_peak_counts[i];

      if( !source_accounts_for_peak(expected, peak) )
      {
        terms.num_excluded += 1;  // background/mismatch peak the source cannot account for - uninformative
        continue;
      }

      terms.num_included += 1;
      const double sys = judge_sys_frac * peak.m_counts;
      const double unc2 = (peak.m_counts_uncert * peak.m_counts_uncert) + (sys * sys);
      if( unc2 > 0.0 )
        terms.chi2 += ((expected - peak.m_counts) * (expected - peak.m_counts)) / unc2;
    }

    return terms;
  };//judgment_chi2_terms

  const auto judgment_chi2_dof = [&judgment_chi2_terms]( const RelActCalcManual::RelEffSolution &sol ) -> double {
    const JudgedChi2Terms terms = judgment_chi2_terms( sol );
    if( terms.chi2 == std::numeric_limits<double>::max() )
      return std::numeric_limits<double>::max();
    // If peaks were excluded as contaminants AND fewer than two informative peaks remain, the curve
    // has nothing to be judged against - reject (fall back) rather than vacuously "accepting" with ~0
    // chi2.  But a legitimately exactly-determined fit with no exclusions (e.g. a source matched to a
    // single peak - very common) must NOT be rejected here; it should be accepted (chi2 ~ 0) so the
    // proper joint rel-act/rel-eff manual solution is used instead of diverting to the fallback.
    if( (terms.num_included < 2) && (terms.num_excluded > 0) )
      return std::numeric_limits<double>::max();
    // Each excluded peak removes a data point but not a fit parameter, so the effective dof drops by
    // the number excluded.  dof<=0 is a legitimate exactly-determined fit (e.g. single-peak source);
    // treat its (near-zero) chi2 as-is rather than rejecting, matching the prior std::max(m_dof,1).
    const int eff_dof = sol.m_dof - terms.num_excluded;
    return (eff_dof > 0) ? (terms.chi2 / eff_dof) : terms.chi2;
  };//judgment_chi2_dof

  // Contaminant-pruned matched peaks (those the source's own fitted curve can account for) - handed to
  // the fallback if the manual fit cannot be made acceptable, so it builds its empirical rel-eff from a
  // clean set without any energy-based (NORM-unsafe) filtering.  Stays empty if the manual solve never
  // converges, in which case the fallback uses its generic-DRF brightest-gamma estimate.
  std::vector<RelActCalcManual::GenericPeakInfo> clean_matched_peaks;

  // Step 2: Solve for relative efficiency, choosing the equation form/order per spectrum by
  // small-sample-corrected AIC (AICc) over a bounded candidate ladder.  The judged chi2
  // (systematic-inflated, contaminant-excluded) is the goodness term; the GA tunes only the
  // complexity-penalty scale config.manual_releff_aicc_penalty (kappa; 2.0 = textbook AIC), which
  // also absorbs the judged chi2's non-textbook scale.  This adapts the model to each spectrum -
  // a 3-peak Cs137+Ba133 spectrum and a 3-peak Eu152 spectrum get different answers - replacing
  // the former fixed form/order-per-peak-count genes, whose single global choice could not be
  // right for all spectra at once, plus the separate physical-model retry (now just a candidate).
  try
  {
    struct ManualRelEffCandidate
    {
      RelActCalc::RelEffEqnForm form;
      size_t order;
    };//struct ManualRelEffCandidate

    // Preference-ordered (strict-improvement selection means earlier candidates win exact ties):
    // LnXLnY first - the historical default, stable at low energy; the physical model last (the
    // only Ceres-fit candidate).
    const RelActCalc::RelEffEqnForm empirical_forms[] = {
      RelActCalc::RelEffEqnForm::LnXLnY, RelActCalc::RelEffEqnForm::LnX,
      RelActCalc::RelEffEqnForm::LnY, RelActCalc::RelEffEqnForm::FramEmpirical
    };

    std::vector<ManualRelEffCandidate> candidates;
    for( size_t order = 0; order <= max_eqn_order; ++order )
    {
      for( const RelActCalc::RelEffEqnForm form : empirical_forms )
        candidates.push_back( { form, order } );
    }
    candidates.push_back( { RelActCalc::RelEffEqnForm::FramPhysicalModel, 0 } );

    const double kappa = config.manual_releff_aicc_penalty;

    RelActCalcManual::RelEffSolution manual_solution;  // m_status defaults to NotInitialized
    double best_metric = std::numeric_limits<double>::max();
    bool best_is_aicc_valid = false;

    for( const ManualRelEffCandidate &cand : candidates )
    {
      RelActCalcManual::RelEffInput cand_input;
      cand_input.peaks = peaks_matched;
      cand_input.eqn_form = cand.form;
      cand_input.eqn_order = cand.order;
      cand_input.use_ceres_to_fit_eqn = false;
      cand_input.phys_model_use_hoerl = false;

      const bool is_physical = (cand.form == RelActCalc::RelEffEqnForm::FramPhysicalModel);
      if( is_physical )
      {
        cand_input.eqn_order = 0;
        cand_input.use_ceres_to_fit_eqn = true;
        cand_input.phys_model_use_hoerl = (peaks_matched.size() > 3);
        cand_input.phys_model_detector = phys_model_det;
        cand_input.phys_model_self_atten = std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>{};
        cand_input.phys_model_external_attens = std::vector<std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>>{};
      }//if( is_physical )

      RelActCalcManual::RelEffSolution sol;
      try
      {
        sol = RelActCalcManual::solve_relative_efficiency( cand_input );
      }catch( const std::exception & )
      {
        continue;
      }

      if( sol.m_status != RelActCalcManual::ManualSolutionStatus::Success )
        continue;

      const JudgedChi2Terms terms = judgment_chi2_terms( sol );
      if( (terms.chi2 == std::numeric_limits<double>::max())
         || ((terms.num_included < 2) && (terms.num_excluded > 0)) )
        continue;  // not judgeable - nothing informative to select on

      // Fitted-parameter count: curve coefficients plus one activity per matched nuclide (the
      // Hoerl modifier adds two).  kappa absorbs any small miscount.
      const double num_par = is_physical
          ? (static_cast<double>(num_distinct_nuclides)
             + (cand_input.phys_model_use_hoerl ? 2.0 : 0.0))
          : (static_cast<double>(cand.order) + 1.0
             + static_cast<double>(num_distinct_nuclides));
      const double num_data = terms.num_included;

      // AICc = chi2 + kappa*k + kappa*k*(k+1)/(n-k-1); the correction term requires n > k+1.
      // When no candidate satisfies that (very few matched peaks), order-0 candidates - which
      // have no flexibility to overfit - compete on chi2 + kappa*k alone.
      const bool aicc_valid = (num_data > (num_par + 1.0));
      double metric;
      if( aicc_valid )
        metric = terms.chi2 + kappa*num_par + (kappa * num_par * (num_par + 1.0)) / (num_data - num_par - 1.0);
      else if( cand.order == 0 )
        metric = terms.chi2 + kappa*num_par;
      else
        continue;

      if( should_debug_print() )
      {
        std::cout << "Manual rel-eff ladder: form=" << RelActCalc::to_str( cand.form )
                  << ", order=" << cand.order << ": judged_chi2=" << terms.chi2
                  << " (" << terms.num_included << " incl, " << terms.num_excluded << " excl)"
                  << ", k=" << num_par << ", metric=" << metric
                  << (aicc_valid ? "" : " [no AICc correction]") << std::endl;
      }

      // Any AICc-valid candidate beats every non-valid one; within a tier, lower metric wins.
      const bool better = best_is_aicc_valid
          ? (aicc_valid && (metric < best_metric))
          : (aicc_valid || (metric < best_metric));
      if( better )
      {
        best_metric = metric;
        best_is_aicc_valid = aicc_valid;
        manual_solution = std::move( sol );
      }
    }//for( const ManualRelEffCandidate &cand : candidates )

    double chi2_dof = judgment_chi2_dof( manual_solution );

    if( should_debug_print()
       && (manual_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success) )
    {
      std::cout << "Manual rel-eff ladder winner: form=" << RelActCalc::to_str( manual_solution.m_input.eqn_form )
                << ", order=" << manual_solution.m_input.eqn_order
                << (manual_solution.m_input.phys_model_external_attens.empty()
                    ? "" : ", external attenuation")
                << ", judged chi2/dof=" << chi2_dof << std::endl;
    }

    if( manual_solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
      throw std::runtime_error( "Failed to fit initial RelActCalcManual::RelEffSolution: " + manual_solution.m_error_message );

    // Outlier-removal-and-refit recovery.  A few matched peaks that don't lie on any smooth rel-eff
    // curve - a handful of lines a too-stiff empirical form can't follow through a steep shielding
    // rolloff (e.g. shielded Ir-192's 468-489 keV region) - can push the judged chi2/dof over the
    // acceptance gate even though the bulk of the source's lines are mutually consistent.  Rather than
    // discard the whole (proper, joint rel-act/rel-eff) solution for the cruder fallback, iteratively
    // drop the single worst-fitting INCLUDED peak and re-solve until the fit is acceptable or too few
    // peaks remain.  The outlier test is curve-relative (per-peak chi2 against the source's OWN fitted
    // curve) and never keys on absolute line energies, so it is safe when the source itself is NORM.
    if( chi2_dof > config.initial_manual_rel_eff_max_chi2_dof )
    {
      const size_t num_src = std::max<size_t>( 1, manual_solution.m_rel_activities.size() );
      const size_t min_peaks = num_src + manual_solution.m_input.eqn_order + 2;
      const int max_removals = 4;

      for( int iter = 0;
           (iter < max_removals)
             && (chi2_dof > config.initial_manual_rel_eff_max_chi2_dof)
             && (manual_solution.m_input.peaks.size() > min_peaks)
             && (manual_solution.m_predicted_peak_counts.size() == manual_solution.m_input.peaks.size());
           ++iter )
      {
        // Worst INCLUDED peak by curve-relative judged chi2.  Contaminants (predicted << observed)
        // carry little curve leverage so removing them rarely lowers chi2 - they are dropped from
        // clean_matched_peaks below instead.
        double worst_c2 = -1.0;
        size_t worst_idx = manual_solution.m_input.peaks.size();
        for( size_t i = 0; i < manual_solution.m_input.peaks.size(); ++i )
        {
          const RelActCalcManual::GenericPeakInfo &pk = manual_solution.m_input.peaks[i];
          const double expected = manual_solution.m_predicted_peak_counts[i];
          if( !source_accounts_for_peak(expected, pk) )
            continue;  // contaminant - not a curve-misfit candidate

          // Don't orphan a nuclide: keep at least one peak for this peak's dominant isotope.
          std::string dom; double dy = 0.0;
          for( const RelActCalcManual::GenericLineInfo &ln : pk.m_source_gammas )
            if( ln.m_yield > dy ){ dy = ln.m_yield; dom = ln.m_isotope; }
          int dom_count = 0;
          for( const RelActCalcManual::GenericPeakInfo &pk2 : manual_solution.m_input.peaks )
          {
            std::string d2; double y2 = 0.0;
            for( const RelActCalcManual::GenericLineInfo &ln2 : pk2.m_source_gammas )
              if( ln2.m_yield > y2 ){ y2 = ln2.m_yield; d2 = ln2.m_isotope; }
            if( d2 == dom ) ++dom_count;
          }
          if( dom.empty() || (dom_count <= 1) )
            continue;

          const double resid = pk.m_counts - expected;
          const double sys = judge_sys_frac * pk.m_counts;
          const double unc2 = (pk.m_counts_uncert*pk.m_counts_uncert) + (sys*sys);
          const double c2 = (unc2 > 0.0) ? (resid*resid/unc2) : 0.0;
          if( c2 > worst_c2 ){ worst_c2 = c2; worst_idx = i; }
        }

        if( worst_idx >= manual_solution.m_input.peaks.size() )
          break;  // nothing removable without orphaning a nuclide

        RelActCalcManual::RelEffInput trial_input = manual_solution.m_input;
        const double dropped_energy = trial_input.peaks[worst_idx].m_energy;
        trial_input.peaks.erase( trial_input.peaks.begin() + static_cast<long>(worst_idx) );

        RelActCalcManual::RelEffSolution trial_sol = RelActCalcManual::solve_relative_efficiency( trial_input );
        const double trial_chi2 = judgment_chi2_dof( trial_sol );
        if( (trial_sol.m_status != RelActCalcManual::ManualSolutionStatus::Success) || (trial_chi2 >= chi2_dof) )
        {
          if( should_debug_print() )
            std::cout << "Outlier removal: dropping " << dropped_energy << " keV did not help (chi2/dof "
                      << chi2_dof << " -> " << trial_chi2 << "); stopping." << std::endl;
          break;
        }

        if( should_debug_print() )
          std::cout << "Outlier removal: dropped " << dropped_energy << " keV (chi2_judged=" << worst_c2
                    << "); chi2/dof " << chi2_dof << " -> " << trial_chi2 << std::endl;
        manual_solution = std::move( trial_sol );
        chi2_dof = trial_chi2;
      }//for( outlier-removal iterations )

      if( should_debug_print() && (chi2_dof <= config.initial_manual_rel_eff_max_chi2_dof) )
        std::cout << "Outlier removal recovered an acceptable manual solution (chi2/dof=" << chi2_dof
                  << ", " << manual_solution.m_input.peaks.size() << " peaks)." << std::endl;
    }//if( chi2_dof > gate ) - outlier removal

    // Contaminant-pruned set for the fallback: keep peaks the (possibly pruned) fitted curve can
    // account for (predicted >= 25% of observed).  Curve-relative, so NORM-safe.
    if( manual_solution.m_predicted_peak_counts.size() == manual_solution.m_input.peaks.size() )
    {
      clean_matched_peaks.clear();
      for( size_t i = 0; i < manual_solution.m_input.peaks.size(); ++i )
      {
        if( source_accounts_for_peak( manual_solution.m_predicted_peak_counts[i],
                                      manual_solution.m_input.peaks[i] ) )
          clean_matched_peaks.push_back( manual_solution.m_input.peaks[i] );
      }
    }

    if( should_debug_print() )
    {
      // Per-peak rel-eff residuals, so we can see whether a bad judged chi2/dof is driven by a few
      // outliers (e.g. peaks contaminated by NORM/background lines) or a broad curve misfit.
      std::cout << "Per-peak rel-eff residuals (judged with " << (100.0*judge_sys_frac) << "% systematic floor), form="
                << RelActCalc::to_str( manual_solution.m_input.eqn_form ) << ":" << std::endl;
      double dbg_chi2_stat = 0.0, dbg_chi2_judged = 0.0;
      int dbg_num_excluded = 0;
      const bool dbg_have_pred = (manual_solution.m_predicted_peak_counts.size() == manual_solution.m_input.peaks.size());
      for( size_t pi = 0; pi < manual_solution.m_input.peaks.size(); ++pi )
      {
        const RelActCalcManual::GenericPeakInfo &peak = manual_solution.m_input.peaks[pi];
        std::string iso;
        for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
          if( iso.empty() ){ iso = line.m_isotope; break; }
        const double expected = dbg_have_pred ? manual_solution.m_predicted_peak_counts[pi] : 0.0;
        const double resid = peak.m_counts - expected;
        const double pct = (expected > 0.0) ? (100.0 * resid / expected) : 0.0;
        const double sys = judge_sys_frac * peak.m_counts;
        const double unc2_judged = (peak.m_counts_uncert*peak.m_counts_uncert) + (sys*sys);
        const double c2_judged = (unc2_judged > 0.0) ? (resid*resid/unc2_judged) : 0.0;
        const bool excluded = !source_accounts_for_peak(expected, peak);
        if( !excluded )
          dbg_chi2_judged += c2_judged;
        else
          ++dbg_num_excluded;
        dbg_chi2_stat += (peak.m_counts_uncert > 0.0) ? (resid*resid/(peak.m_counts_uncert*peak.m_counts_uncert)) : 0.0;
        const bool near_norm = is_near_strong_norm_gamma( peak.m_energy, std::max(1.0, peak.m_fwhm) );
        std::cout << "    " << peak.m_energy << " keV " << iso
                  << ": data=" << peak.m_counts << " pred=" << expected
                  << " resid=" << pct << "% chi2_judged=" << c2_judged
                  << (excluded ? "   [EXCLUDED: source<25% of peak]" : "")
                  << (near_norm ? " (near NORM)" : "") << std::endl;
      }
      const int dbg_eff_dof = manual_solution.m_dof - dbg_num_excluded;
      std::cout << "  Sum chi2: stat-only(all)=" << dbg_chi2_stat << ", judged(included)=" << dbg_chi2_judged
                << ", excluded " << dbg_num_excluded << " peaks, eff_dof=" << dbg_eff_dof
                << ", judged chi2/dof=" << (dbg_eff_dof > 0 ? dbg_chi2_judged/dbg_eff_dof : dbg_chi2_judged) << std::endl;
    }//if( should_debug_print() )

    // Reject extremely bad fits.  If matched peaks can't be reconciled with any
    // smooth relative-efficiency curve, the source(s) probably aren't really
    // present and we should treat this the same as a hard failure so the caller
    // falls back to `estimate_initial_rois_fallback`.
    const double final_chi2_dof = judgment_chi2_dof( manual_solution );
    if( final_chi2_dof > config.initial_manual_rel_eff_max_chi2_dof )
    {
      throw std::runtime_error(
        "Initial RelActCalcManual::RelEffSolution chi2/dof="
        + std::to_string(final_chi2_dof)
        + " exceeds threshold "
        + std::to_string(config.initial_manual_rel_eff_max_chi2_dof)
        + " - matched peaks are not consistent with the source(s) being present" );
    }

    if( should_debug_print() )
    {
      std::cout << "Successfully fitted initial RelActCalcManual::RelEffSolution: chi2/dof="
      << manual_solution.m_chi2 << "/" << manual_solution.m_dof << "="
      << (manual_solution.m_chi2 / (std::max)(manual_solution.m_dof, 1))
      << " using " << peaks_matched.size() << " peaks"
      << std::endl;
      std::cout << std::endl;
    }

    // Step 3: Create rel_eff lambda from manual solution - handles extrapolation clamping.
    // Lower side: flat-clamp to the lowest matched-peak energy (downward extrapolation unreliable).
    // Upper side: DRF-shaped extrapolation anchored to the boundary rel-eff.
    const std::shared_ptr<const DetectorPeakResponse> rel_eff_extrap_drf
      = generic_drf_for_rel_eff_extrap( drf, det_type );
    const auto manual_rel_eff = [&manual_solution, &peaks_matched, rel_eff_extrap_drf]( double energy ) -> double {
      if( energy < peaks_matched.front().m_energy )
        return manual_solution.relative_efficiency( peaks_matched.front().m_energy );
      const double e_hi = peaks_matched.back().m_energy;
      if( energy <= e_hi )
        return manual_solution.relative_efficiency( energy );
      const double re_hi = manual_solution.relative_efficiency( e_hi );
      return shape_rel_eff_above_boundary( re_hi, energy, e_hi, rel_eff_extrap_drf );
    };

    // Step 4: Collect sources and activities from the manual solution
    vector<tuple<RelActCalcAuto::SrcVariant, double, double>> source_age_and_acts;
    for( const RelActCalcManual::IsotopeRelativeActivity &rel_act : manual_solution.m_rel_activities )
    {
      const RelActCalcAuto::SrcVariant src = RelActCalcAuto::source_from_string( rel_act.m_isotope );
      if( RelActCalcAuto::is_null( src ) )
        throw std::logic_error( "Failed to get source from RelAct isotope '" + rel_act.m_isotope + "'" );

      double age = 0.0;
      if( RelActCalcAuto::nuclide(src) )
      {
        bool found_src = false;
        for( const RelActCalcAuto::NucInputInfo &input_src : sources )
        {
          if( RelActCalcAuto::to_name(input_src.source) == rel_act.m_isotope )
          {
            age = input_src.age;
            found_src = true;
            break;
          }
        }//for( const RelActCalcAuto::NucInputInfo &input_src : sources )

        assert( found_src );
      }//if( RelActCalcAuto::nuclide(src) )

      const double act = manual_solution.relative_activity( rel_act.m_isotope );
      source_age_and_acts.emplace_back( src, age, act );
    }//for( loop over manual_solution.m_rel_activities )


    // Step 5: Use the reusable clustering function to create ROIs
    {
      const bool have_policy_fwhm_range = (lower_fwhm_energy > 0.0)
        && (upper_fwhm_energy > lower_fwhm_energy);
      const auto policy_fwhm = [=, &fwhm_coefficients]( const double energy ) {
        const double eval_energy = have_policy_fwhm_range
          ? std::clamp(energy, lower_fwhm_energy, upper_fwhm_energy) : energy;
        return static_cast<double>( DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(eval_energy), fwhmFnctnlForm, fwhm_coefficients) );
      };
      detail::GlobalContinuumEstimate initial_continuum;
      if( manual_settings.use_automatic_roi_policy )
      {
        initial_continuum = detail::make_global_continuum(
            foreground, policy_fwhm, det_type, min_valid_energy, max_valid_energy );
      }
      const detail::GlobalContinuumEstimate * const initial_continuum_ptr
        = initial_continuum.valid() ? &initial_continuum : nullptr;
      GammaClusteringSettings policy_settings = manual_settings;
      if( manual_settings.use_automatic_roi_policy )
        policy_settings.global_continuum = initial_continuum_ptr;
      const std::vector<std::shared_ptr<const PeakDef>> no_unfit_peaks;
      const std::vector<std::shared_ptr<const PeakDef>> &clustering_unfit_peaks
        = manual_settings.use_automatic_roi_policy ? unfit_auto_peaks : no_unfit_peaks;
      const std::vector<std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo>> rois_and_gammas
        = cluster_gammas_to_rois( {manual_rel_eff}, {source_age_and_acts}, foreground,
                                  fwhmFnctnlForm, fwhm_coefficients,
                                  lower_fwhm_energy, upper_fwhm_energy,
                                  min_valid_energy, max_valid_energy,
                                  policy_settings, clustering_unfit_peaks,
                                  nullptr, nullptr, nullptr, nullptr,
                                  "initial manual", initial_continuum_ptr );

      initial_rois.clear();
      initial_rois.reserve( rois_and_gammas.size() );
      for( const std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &p : rois_and_gammas )
      {
        initial_rois.push_back( p.first );
        if( modeled_peak_candidates )
        {
          for( size_t peak_index = 0; peak_index < p.second.gamma_energies.size(); ++peak_index )
          {
            const double area = (peak_index < p.second.gamma_amplitudes.size())
              ? p.second.gamma_amplitudes[peak_index] : 0.0;
            modeled_peak_candidates->emplace_back( p.second.gamma_energies[peak_index], area );
          }
        }
      }
    }

    if( should_debug_print() )
    {
      std::cout << "Initial ROIs from RelActManual: ";
      for( const RelActCalcAuto::RoiRange &roi : initial_rois )
        std::cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
      std::cout << std::endl;
    }
  }catch( std::exception &e )
  {
    if( should_debug_print() )
    {
      std::cerr << "Error trying to fit initial manual rel-eff solution: " << e.what() << std::endl;
      std::cerr << "Using fallback activity estimation..." << std::endl;
    }

    initial_rois = estimate_initial_rois_fallback(
      auto_search_peaks, clean_matched_peaks, foreground, sources, drf, det_type,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      manual_settings, unfit_auto_peaks, modeled_peak_candidates );

    // Do not hand every matcher coincidence to the generic fallback: it remains deliberately
    // unchanged when the manual curve is unavailable.  A later source-clean challenger may use
    // only independently significant, distinct requested-source matches.  This preserves a real
    // measured anchor through a manual-solve failure without treating it as a trusted rel-eff
    // calibration point or protecting it from the ordinary final fit/refit significance checks.
    if( manual_settings.use_automatic_roi_policy )
    {
      provisional_fallback_source_anchors = distinct_significant_source_anchors(
          peaks_matched, config.roi_significance_z );
      if( provisional_fallback_source_anchors.size() < 2 )
        provisional_fallback_source_anchors.clear();
      else
      {
        if( has_provisional_fallback_source_anchors )
          *has_provisional_fallback_source_anchors = true;
        if( should_debug_print() )
        {
          std::cout << "Manual rel-eff fallback retained "
                    << provisional_fallback_source_anchors.size()
                    << " data-supported matched source anchors for the transactional challenger: ";
          for( const RelActCalcManual::GenericPeakInfo &anchor : provisional_fallback_source_anchors )
            std::cout << anchor.m_energy << " keV, ";
          std::cout << std::endl;
        }
      }
    }

    fallback_warning = "RelActManual fitting failed (" + std::string(e.what())
                     + "); used empirical peak-area activity estimation fallback";

    if( should_debug_print() )
    {
      std::cout << "Fallback ROIs: ";
      for( const RelActCalcAuto::RoiRange &roi : initial_rois )
        std::cout << "[" << roi.lower_energy << ", " << roi.upper_energy << "], ";
      std::cout << std::endl;
    }
  }

  // Change 3: guarantee a tight ROI for every found+matched auto-search peak the clustering did not
  // already cover.  These are data-confirmed source lines (matched to a requested source by
  // fill_in_nuclide_info above); the predicted-signal keep-gate can wrongly drop them on wide-FWHM
  // low-count NaI, and losing them was the root of the "No ROIs"/silent-empty failures.  Runs on both
  // the manual-success and fallback paths.  peaks_matched is non-empty here (the empty case returned
  // early via estimate_initial_rois_without_peaks).  [architecture review 2026-07-18]
  {
    const std::vector<RelActCalcAuto::RoiRange> clustered_rois = initial_rois;
    std::vector<std::pair<double,double>> found_energy_fwhm;
    found_energy_fwhm.reserve( peaks_matched.size() );
    for( const RelActCalcManual::GenericPeakInfo &pk : peaks_matched )
    {
      found_energy_fwhm.emplace_back( pk.m_energy, pk.m_fwhm );
      if( modeled_peak_candidates )
      {
        const std::vector<std::pair<double,double>>::iterator existing = std::find_if(
            std::begin(*modeled_peak_candidates), std::end(*modeled_peak_candidates),
            [&pk]( const std::pair<double,double> &candidate ) {
              return std::fabs(candidate.first - pk.m_energy) < 1.0e-6;
            } );
        if( existing == std::end(*modeled_peak_candidates) )
          modeled_peak_candidates->emplace_back( pk.m_energy, pk.m_counts );
        else
          existing->second = std::max( existing->second, pk.m_counts );
      }
    }

    const size_t num_seeded = seed_tight_rois_for_found_peaks(
        initial_rois, found_energy_fwhm, sm_found_peak_roi_half_num_fwhm,
        min_valid_energy, max_valid_energy, manual_settings.use_automatic_roi_policy );

    if( should_debug_print() && num_seeded )
      std::cout << "Seeded " << num_seeded << " tight ROI(s) for found+matched auto-search peak(s)"
                << " not covered by clustering" << std::endl;

    const std::vector<RelActCalcManual::GenericPeakInfo> &recovery_source_anchors
      = clean_matched_peaks.empty() ? provisional_fallback_source_anchors : clean_matched_peaks;
    if( clean_source_rois && manual_settings.use_automatic_roi_policy )
    {
      // Provisional anchors have not survived a manual curve fit.  Their challenger must start
      // from only their tight, measured ROIs rather than inherit a generic-fallback ROI that was
      // admitted by an unconstrained extrapolation.
      *clean_source_rois = provisional_fallback_source_anchors.empty()
        ? clustered_rois : std::vector<RelActCalcAuto::RoiRange>();
      std::vector<std::pair<double,double>> clean_energy_fwhm;
      clean_energy_fwhm.reserve( recovery_source_anchors.size() );
      for( const RelActCalcManual::GenericPeakInfo &pk : recovery_source_anchors )
        clean_energy_fwhm.emplace_back( pk.m_energy, pk.m_fwhm );
      seed_tight_rois_for_found_peaks( *clean_source_rois, clean_energy_fwhm,
          sm_found_peak_roi_half_num_fwhm, min_valid_energy, max_valid_energy, true );
    }

    if( source_anchor_candidates )
      *source_anchor_candidates = recovery_source_anchors;
  }

  if( source_anchor_candidates && source_anchor_candidates->empty() )
    *source_anchor_candidates = clean_matched_peaks;

  return initial_rois;
}//estimate_initial_rois_using_relactmanual

PeakFitResult fit_peaks_for_nuclide_relactauto(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &orig_foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::vector<RelActCalcAuto::RoiRange> &input_rois,
  const std::vector<RelActCalcAuto::RoiRange> &clean_source_rois,
  const std::vector<std::pair<double,double>> &initial_modeled_peak_candidates,
  const std::vector<RelActCalcManual::GenericPeakInfo> &source_anchor_candidates,
  const bool has_provisional_fallback_source_anchors,
  const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &orig_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const Wt::WFlags<FitSrcPeaksOptions> user_options,
  const PeakFitForNuclideConfig &config,
  const DetectorPeakResponse::ResolutionFnctForm fwhm_form,
  const std::vector<float> &fwhm_coefficients,
  const PeakFitUtils::CoarseResolutionType det_type,
  const double fwhm_lower_energy,
  const double fwhm_upper_energy,
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs )
{
  PeakFitResult result;
  result.automatic_roi_diagnostics = detail::take_automatic_roi_diagnostics();

  const bool fit_norm_peaks = user_options.test(FitSrcPeaksOptions::FitNormBkgrndPeaks)
                              || user_options.test(FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse);
  const bool norm_peaks_dont_use = user_options.test(FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse);
  const bool disable_auto_interferer_fit
    = user_options.test(FitSrcPeaksOptions::DisableAutoInterfererFit);
  const bool use_automatic_roi_policy = fit_norm_peaks || disable_auto_interferer_fit;
  const bool apply_energy_cal_between = config.fit_energy_cal
                                        && !user_options.test(FitSrcPeaksOptions::DoNotVaryEnergyCal)
                                        && !user_options.test(FitSrcPeaksOptions::DoNotRefineEnergyCal);
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );

  // Validate sources
  if( sources.empty() && !fit_norm_peaks )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "No sources provided";
    return result;
  }

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Null source in sources vector";
      return result;
    }
  }


  // Use skew type from config
  const PeakDef::SkewType skew_type = config.skew_type;

  // Create relative efficiency curve from config parameters
  RelActCalcAuto::RelEffCurveInput rel_eff_curve;
  rel_eff_curve.rel_eff_eqn_type = config.rel_eff_eqn_type;
  rel_eff_curve.rel_eff_eqn_order = config.rel_eff_eqn_order;

  // FramPhysicalModel requires rel_eff_eqn_order to be 0
  assert( (config.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel)
         || (rel_eff_curve.rel_eff_eqn_order == 0) );
  if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    rel_eff_curve.rel_eff_eqn_order = 0;

  // Determine nucs_of_el_same_age and potentially synchronize ages for nuclides of the same element.
  // If no ages are being fit and we have multiple nuclides of the same element, set to false.
  // If ages are being fit, synchronize starting ages: for U use U235's age, for Pu use Pu239's age,
  // otherwise use the median age of isotopes for that element.
  bool any_age_fit = false;
  std::map<int, std::vector<size_t>> element_to_source_indices; // atomic number -> indices in sources
  
  for( size_t i = 0; i < sources.size(); ++i )
  {
    const RelActCalcAuto::NucInputInfo &src = sources[i];
    if( src.fit_age )
      any_age_fit = true;
    
    const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( src.source );
    if( nuc )
      element_to_source_indices[nuc->atomicNumber].push_back( i );
  }//for( loop over sources )
  
  // Check if any element has multiple nuclides
  bool has_multi_nuc_element = false;
  for( const std::pair<const int, std::vector<size_t>> &el_indices : element_to_source_indices )
  {
    if( el_indices.second.size() > 1 )
    {
      has_multi_nuc_element = true;
      break;
    }
  }//for( loop over element_to_source_indices )
  
  // Make a mutable copy of sources for potential age synchronization
  std::vector<RelActCalcAuto::NucInputInfo> synced_sources = sources;
  
  if( has_multi_nuc_element )
  {
    if( !any_age_fit )
    {
      // No ages being fit and multiple nuclides of same element - set to false
      rel_eff_curve.nucs_of_el_same_age = false;
    }else
    {
      // Ages are being fit - synchronize starting ages for nuclides of the same element
      rel_eff_curve.nucs_of_el_same_age = config.nucs_of_el_same_age;
      
      for( const std::pair<const int, std::vector<size_t>> &el_indices : element_to_source_indices )
      {
        const int atomic_num = el_indices.first;
        const std::vector<size_t> &indices = el_indices.second;
        
        if( indices.size() <= 1 )
          continue;
        
        // Determine the age to use for this element
        double age_to_use = -1.0;
        
        // Special case for U (Z=92): use U235's age if present
        if( atomic_num == 92 )
        {
          for( const size_t idx : indices )
          {
            const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( synced_sources[idx].source );
            if( nuc && (nuc->massNumber == 235) && (nuc->isomerNumber == 0) )
            {
              age_to_use = synced_sources[idx].age;
              break;
            }
          }//for( loop over indices )
        }else if( atomic_num == 94 )
        {
          // Special case for Pu (Z=94): use Pu239's age if present
          for( const size_t idx : indices )
          {
            const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( synced_sources[idx].source );
            if( nuc && (nuc->massNumber == 239) && (nuc->isomerNumber == 0) )
            {
              age_to_use = synced_sources[idx].age;
              break;
            }
          }//for( loop over indices )
        }//if( U ) else if( Pu )
        
        // If no special case matched, use the median age
        if( age_to_use < 0.0 )
        {
          std::vector<double> ages;
          ages.reserve( indices.size() );
          for( const size_t idx : indices )
          {
            if( synced_sources[idx].age >= 0.0 )
              ages.push_back( synced_sources[idx].age );
          }
          
          if( !ages.empty() )
          {
            std::sort( ages.begin(), ages.end() );
            const size_t mid = ages.size() / 2;
            if( ages.size() % 2 == 0 )
              age_to_use = 0.5 * (ages[mid - 1] + ages[mid]);
            else
              age_to_use = ages[mid];
          }
        }//if( age_to_use < 0.0 )
        
        // Apply the synchronized age to all nuclides of this element
        if( age_to_use >= 0.0 )
        {
          for( const size_t idx : indices )
            synced_sources[idx].age = age_to_use;
        }
      }//for( loop over element_to_source_indices )
    }//if( !any_age_fit ) / else
  }else
  {
    // No element has multiple nuclides - use config value
    rel_eff_curve.nucs_of_el_same_age = config.nucs_of_el_same_age;
  }//if( has_multi_nuc_element ) / else
  
  
  rel_eff_curve.phys_model_corr.corr_fcn = config.phys_model_use_hoerl
                            ? RelActCalc::PhysModelCorrFcn::Hoerl : RelActCalc::PhysModelCorrFcn::None;
  rel_eff_curve.nuclides = synced_sources;

  // Copy shielding inputs from config (converting to const shared_ptr)
  if( !config.phys_model_self_atten.empty() )
    rel_eff_curve.phys_model_self_atten = config.phys_model_self_atten.front();

  rel_eff_curve.phys_model_external_atten.clear();
  for( const auto &shield : config.phys_model_external_atten )
    rel_eff_curve.phys_model_external_atten.push_back( shield );

  // Generate name based on equation type and shielding (just for informational purposes)
  if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    rel_eff_curve.name = "Physical Model";
    if( config.phys_model_external_atten.empty() )
    {
      rel_eff_curve.name += " (no shielding) Peak Fit";
    }else
    {
      const int z = static_cast<int>( config.phys_model_external_atten.front()->atomic_number + 0.5 );
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      const SandiaDecay::Element *el = db->element(z);
      assert( el );
      rel_eff_curve.name += " " + (el ? el->symbol : ("z=" + std::to_string(z)) ) + std::string(" Peak Fit");
    }
  }else
  {
    rel_eff_curve.name = RelActCalc::to_str(config.rel_eff_eqn_type) + std::string(" Peak Fit");
  }
  
  // Track which rel_eff_curves index corresponds to sources vs NORM;
  //  either or both may be valid (-1 means not present).
  int sources_rel_eff_index = -1, norm_rel_eff_index = -1;

  // R6 auto co-fit of strong unmodeled interferers: nuclides added on an extra rel-eff curve (their
  // peaks are dropped from the returned results), and the energies of any auto-added floating
  // interferer peaks.  Both are consulted in the result-filter block after the solve.
  std::set<const SandiaDecay::Nuclide *> auto_interferer_nucs;
  std::vector<double> auto_interferer_float_energies;
  std::vector<double> auto_interferer_lines;
  std::vector<double> interferer_guard_energies;
  std::vector<detail::InterfererCandidate> interferer_candidates;

  // Energies (original spectrum cal) of existing user peaks carried into the fit as bystander
  // FloatingPeaks.  Their identity is reported by the bystander itself after the solve (either the
  // FloatingPeakResult-updated copy we append to `observable_peaks`, or the users untouched
  // original when it could not be updated), so the corresponding unassigned model peak must not
  // ALSO be reported - see `is_carried_bystander_float` in the result-filter block.
  std::vector<double> carried_bystander_float_energies;

  // Create RelActAuto options from config
  RelActCalcAuto::Options options;
  if( !rel_eff_curve.nuclides.empty() )
  {
    sources_rel_eff_index = 0;
    options.rel_eff_curves.push_back( rel_eff_curve );
  }
  
  options.rois = input_rois;
  //options.energy_cal_type = config.fit_energy_cal
  //  ? RelActCalcAuto::EnergyCalFitType::NonLinearFit
  //  : RelActCalcAuto::EnergyCalFitType::NoFit;
  options.energy_cal_type = user_options.test(FitSrcPeaksOptions::DoNotVaryEnergyCal)
                              ? RelActCalcAuto::EnergyCalFitType::NoFit
                              : RelActCalcAuto::EnergyCalFitType::NonLinearFit;
  
  options.fwhm_form = config.fwhm_form;
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
  options.skew_type = skew_type;
  options.additional_br_uncert = config.rel_eff_auto_base_rel_eff_uncert;

  // Copy fixed skew parameter values from PeakFitDetPrefs, if available and not ROI-independent
  if( peak_fit_prefs && !peak_fit_prefs->m_roi_independent_skew )
  {
    for( size_t i = 0; i < 4; ++i )
    {
      options.fixed_lower_skew[i] = peak_fit_prefs->m_lower_energy_skew[i];
      options.fixed_upper_skew[i] = peak_fit_prefs->m_upper_energy_skew[i];
    }
  }//if( peak_fit_prefs && !roi_independent )

  // Find valid energy range, clamped to a physically-valid low-energy floor (see low_energy_analysis_floor).
  const std::pair<double,double> raw_valid_range = find_valid_energy_range( orig_foreground );
  const double low_e_floor = low_energy_analysis_floor( drf, det_type );
  const double min_valid_energy = (low_e_floor < raw_valid_range.second)
                                  ? std::max( raw_valid_range.first, low_e_floor ) : raw_valid_range.first;
  const double max_valid_energy = raw_valid_range.second;

  // R1 step 2: build the shared SNIP global continuum ONCE over the valid extent, so all the gating
  // B(E) estimates below (keep gate, step-gate dominance, edge-ROI restore) reason about the same
  // B(E) instead of each re-estimating an unreliable local two-sideband line.  Per-class SNIP params
  // live in make_global_continuum; if it fails to build, `valid()` is false and every consumer
  // transparently falls back to its prior local estimate.
  const auto global_fwhm_at = [&]( double e ) -> double {
    const bool have = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0)
                      && (fwhm_lower_energy < fwhm_upper_energy);
    const double ee = have ? std::clamp( e, fwhm_lower_energy, fwhm_upper_energy ) : e;
    return DetectorPeakResponse::peakResolutionFWHM( static_cast<float>(ee), fwhm_form, fwhm_coefficients );
  };
  // Not const: when the refinement loop advances the energy calibration, we re-stamp this estimate's
  // calibration to match the working foreground (see the cal-advance block below).
  detail::GlobalContinuumEstimate global_cont = detail::make_global_continuum(
      orig_foreground, global_fwhm_at, det_type, min_valid_energy, max_valid_energy );

  const bool do_not_use_existing_rois = user_options.test( FitSrcPeaksOptions::DoNotUseExistingRois );
  const bool existing_peaks_as_free   = user_options.test( FitSrcPeaksOptions::ExistingPeaksAsFreePeak );

  assert( !(do_not_use_existing_rois && existing_peaks_as_free) );

  // Helper: returns true if a peak's assigned source matches any of the input sources (or NORM nuclides
  // when fit_norm_peaks is requested).  Used by ExistingPeaksAsFreePeak to decide whether an existing
  // peak is "owned" by this fit or is an unrelated bystander.
  const auto peak_source_is_in_fit = [&]( const std::shared_ptr<const PeakDef> &peak ) -> bool {
    const SandiaDecay::Nuclide    *peak_nuc  = peak->parentNuclide();
    const SandiaDecay::Element    *peak_el   = peak->xrayElement();
    const ReactionGamma::Reaction *peak_rxn  = peak->reaction();

    if( !peak_nuc && !peak_el && !peak_rxn )
      return false;  // no source assigned

    // Check against input sources
    for( const RelActCalcAuto::NucInputInfo &src : sources )
    {
      if( peak_nuc && (RelActCalcAuto::nuclide(src.source) == peak_nuc) )
        return true;
      if( peak_el  && (RelActCalcAuto::element(src.source) == peak_el) )
        return true;
      if( peak_rxn && (RelActCalcAuto::reaction(src.source) == peak_rxn) )
        return true;
    }

    // Check against NORM nuclides if we are fitting NORM
    if( fit_norm_peaks && peak_nuc )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      assert( db );
      if( db )
      {
        const SandiaDecay::Nuclide *norm_nucs[] = {
          db->nuclide("U238"), db->nuclide("Ra226"), db->nuclide("U235"),
          db->nuclide("Th232"), db->nuclide("K40")
        };
        for( const SandiaDecay::Nuclide *norm_nuc : norm_nucs )
        {
          assert( norm_nuc );
          if( norm_nuc && (peak_nuc == norm_nuc) )
            return true;
        }//
      }//if( db )
    }//if( fit_norm_peaks && peak_nuc )

    return false;
  };// peak_source_is_in_fit


  // Build existing ROI ranges at function scope so the refinement loop can also filter against them.
  // Each entry stores the ROI bounds plus the peak means within each ROI, so that the
  // trimming logic can enforce FWHM-based margins from existing peak means.
  // Populated in both DoNotUseExistingRois mode (all user ROIs) and default mode (other-source ROIs only).
  struct ExistingRoiInfo
  {
    double lower_energy;
    double upper_energy;
    vector<double> peak_means;
  };

  vector<ExistingRoiInfo> existing_roi_ranges;
  std::vector<RelActCalcAuto::RoiRange> protected_mixed_rois;

  // In default mode, tracks bystander peaks from mixed ROIs (existing ROIs containing both
  // same-source and other-source peaks). After the fit, these bystanders are updated with
  // FloatingPeakResult parameters and added to observable_peaks.
  std::vector<std::pair<std::shared_ptr<const PeakDef>, double>> default_mode_bystander_peaks;

  if( !existing_peaks_as_free && !user_peaks.empty() )
  {
    if( do_not_use_existing_rois )
    {
      // DoNotUseExistingRois mode: all user peak ROIs are excluded from the fit
      map<const PeakContinuum *, size_t> continuum_to_index;
      for( const std::shared_ptr<const PeakDef> &peak : user_peaks )
      {
        const PeakContinuum * const cont = peak->continuum().get();
        const map<const PeakContinuum *, size_t>::const_iterator it = continuum_to_index.find( cont );
        if( it == continuum_to_index.end() )
        {
          continuum_to_index[cont] = existing_roi_ranges.size();
          ExistingRoiInfo info;
          info.lower_energy = peak->lowerX();
          info.upper_energy = peak->upperX();
          info.peak_means.push_back( peak->mean() );
          existing_roi_ranges.push_back( info );
        }
        else
        {
          existing_roi_ranges[it->second].peak_means.push_back( peak->mean() );
        }
      }
    }else
    {
      // Default mode: group user peaks by ROI (shared PeakContinuum), classify each ROI as
      // "all other-source" or "mixed" (contains peaks for both the source being fit and other sources).
      // All-other-source ROIs get added to existing_roi_ranges for trimming.
      // Mixed ROIs have their other-source peaks added as bystander FloatingPeaks, and the ROI
      // is re-created using the existing user peak continuum bounds (not auto-search estimates).
      // Note: In the future, we may also trim existing ROIs so they don't overlap the new peaks,
      // but for now we only trim the new ROIs.
      std::map<const PeakContinuum *, std::vector<std::shared_ptr<const PeakDef>>> roi_groups;
      for( const std::shared_ptr<const PeakDef> &peak : user_peaks )
      {
        if( peak )
          roi_groups[peak->continuum().get()].push_back( peak );
      }

      for( const std::pair<const PeakContinuum *const, std::vector<std::shared_ptr<const PeakDef>>> &roi_group : roi_groups )
      {
        const std::vector<std::shared_ptr<const PeakDef>> &peaks = roi_group.second;
        assert( !peaks.empty() );

        bool has_fit_source = false;
        bool has_other_source = false;
        for( const std::shared_ptr<const PeakDef> &p : peaks )
        {
          if( peak_source_is_in_fit( p ) )
            has_fit_source = true;
          else
            has_other_source = true;
        }

        if( !has_fit_source )
        {
          // All other-source: add to existing_roi_ranges for trimming
          ExistingRoiInfo info;
          info.lower_energy = peaks.front()->lowerX();
          info.upper_energy = peaks.front()->upperX();
          for( const std::shared_ptr<const PeakDef> &p : peaks )
            info.peak_means.push_back( p->mean() );
          existing_roi_ranges.push_back( info );
        }else if( has_other_source )
        {
          // Mixed ROI: contains both same-source and other-source peaks.
          // Add other-source peaks as bystander FloatingPeaks so the fit accounts for them.
          for( const std::shared_ptr<const PeakDef> &p : peaks )
          {
            if( !peak_source_is_in_fit( p ) )
            {
              RelActCalcAuto::FloatingPeak fp;
              fp.energy = p->mean();
              fp.release_fwhm = false;
              // Known (not ObservedInSpectrum) so the bystander tracks the fitted energy-cal
              //  adjustment like the source peaks - see the ExistingPeaksAsFreePeak path below.
              fp.energy_origin = RelActCalcAuto::FloatingPeak::EnergyType::Known;
              options.floating_peaks.push_back( fp );

              default_mode_bystander_peaks.emplace_back( p, p->mean() );
              carried_bystander_float_energies.push_back( p->mean() );
            }
          }

          // Create an explicit RoiRange using the existing user ROI bounds, preserving
          // the existing ROI extent rather than using auto-search estimates.
          const double existing_lower = peaks.front()->lowerX();
          const double existing_upper = peaks.front()->upperX();

          RelActCalcAuto::RoiRange mixed_roi;
          mixed_roi.lower_energy = existing_lower;
          mixed_roi.upper_energy = existing_upper;
          mixed_roi.continuum_type = peaks.front()->continuum()->type();
          mixed_roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;

          // Remove any auto-search-estimated ROI that overlaps this mixed ROI
          // (the existing ROI takes precedence)
          options.rois.erase(
            std::remove_if( options.rois.begin(), options.rois.end(),
              [existing_lower, existing_upper]( const RelActCalcAuto::RoiRange &r ) -> bool {
                return (r.lower_energy < existing_upper) && (r.upper_energy > existing_lower);
              } ),
            options.rois.end()
          );

          options.rois.push_back( mixed_roi );
          protected_mixed_rois.push_back( mixed_roi );
        }
        // else: all same-source, no trimming needed (fit will replace these peaks)
      }//for( roi_groups )
    }//if( do_not_use_existing_rois ) / else

    // Sort existing_roi_ranges by energy for deterministic trimming order
    std::sort( existing_roi_ranges.begin(), existing_roi_ranges.end(),
      []( const ExistingRoiInfo &a, const ExistingRoiInfo &b ) {
        return a.lower_energy < b.lower_energy;
      } );

    // Filter options.rois: reduce, or remove ROIs that overlap existing other-source ROIs
    if( !existing_roi_ranges.empty() )
    {
      std::vector<RelActCalcAuto::RoiRange> filtered_rois;
      filtered_rois.reserve( options.rois.size() );

      const double min_roi_channels = 5.0;    // Minimum ROI size in channels

      const std::shared_ptr<const SpecUtils::EnergyCalibration> energy_cal = orig_foreground->energy_calibration();
      assert( energy_cal && energy_cal->valid() );

      // Helper to find the nearest source gamma energy in [lo, hi] for the sources being fit
      //  (including NORM if applicable). Returns 0.0 if none found.
      const auto nearest_source_gamma_in_range = [&]( const double lo, const double hi,
                                                      const bool near_lower ) -> double
      {
        double best_energy = 0.0;
        double best_dist = std::numeric_limits<double>::max();
        const double ref = near_lower ? lo : hi;

        const auto check_photons = [&]( const vector<SandiaDecay::EnergyRatePair> &photons )
        {
          for( const SandiaDecay::EnergyRatePair &ph : photons )
          {
            if( (ph.energy >= lo) && (ph.energy <= hi)
               && (ph.numPerSecond > std::numeric_limits<float>::epsilon()) )
            {
              const double dist = std::fabs( ph.energy - ref );
              if( dist < best_dist )
              {
                best_dist = dist;
                best_energy = ph.energy;
              }
            }
          }
        };

        for( const RelActCalcAuto::NucInputInfo &src : sources )
        {
          if( RelActCalcAuto::is_null( src.source ) )
            continue;
          check_photons( get_source_photons( src.source, 1.0, src.age ) );
        }

        if( (best_energy <= 0.0) && fit_norm_peaks )
        {
          const SandiaDecay::Nuclide *norm_nucs[] = {
            db->nuclide("U238"), db->nuclide("Ra226"), db->nuclide("U235"),
            db->nuclide("Th232"), db->nuclide("K40")
          };
          for( const SandiaDecay::Nuclide *norm_nuc : norm_nucs )
          {
            if( norm_nuc )
              check_photons( get_source_photons( norm_nuc, 1.0, 0.0 ) );
          }
        }

        return best_energy;
      };

      for( RelActCalcAuto::RoiRange roi : options.rois )
      {
        // Check for overlaps with existing ROIs and reduce bounds if needed
        for( const ExistingRoiInfo &existing : existing_roi_ranges )
        {
          // Check if there's an overlap
          if( (roi.lower_energy < existing.upper_energy) && (roi.upper_energy > existing.lower_energy) )
          {
            // Four mutually exclusive overlap cases:
            const bool new_contains_existing = (roi.lower_energy < existing.lower_energy) && (roi.upper_energy > existing.upper_energy);
            const bool new_contained_by_existing = (roi.lower_energy >= existing.lower_energy) && (roi.upper_energy <= existing.upper_energy);
            const bool lower_edge_inside = (roi.lower_energy >= existing.lower_energy) && (roi.lower_energy < existing.upper_energy);
            const bool upper_edge_inside = (roi.upper_energy > existing.lower_energy) && (roi.upper_energy <= existing.upper_energy);

            if( new_contained_by_existing )
            {
              // New ROI is completely inside the existing ROI - discard it entirely
              roi.lower_energy = -1.0;
              roi.upper_energy = -1.0;
              break;
            }else if( new_contains_existing )
            {
              // New ROI completely contains the existing ROI - keep the larger segment
              const double left_width = existing.lower_energy - roi.lower_energy;
              const double right_width = roi.upper_energy - existing.upper_energy;
              if( left_width > right_width )
                roi.upper_energy = existing.lower_energy;
              else
                roi.lower_energy = existing.upper_energy;
            }else if( lower_edge_inside )
            {
              // New ROI's lower edge is inside the existing ROI - allow abutting, but ensure
              //  the nearest source gamma still has at least 1 FWHM of extent below it.
              const double nearest_gamma = nearest_source_gamma_in_range(
                existing.upper_energy, roi.upper_energy, true );
              if( nearest_gamma > 0.0 )
              {
                const bool have_fwhm_range = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0)
                                             && (fwhm_lower_energy < fwhm_upper_energy);
                const double fwhm_eval = have_fwhm_range
                  ? std::clamp( nearest_gamma, fwhm_lower_energy, fwhm_upper_energy ) : nearest_gamma;
                const double gamma_fwhm = DetectorPeakResponse::peakResolutionFWHM(
                  static_cast<float>(fwhm_eval), fwhm_form, fwhm_coefficients );
                const double min_lower = nearest_gamma - std::max( gamma_fwhm, 1.0 );
                roi.lower_energy = std::max( existing.upper_energy, min_lower );
              }else
              {
                roi.lower_energy = existing.upper_energy;
              }
            }else if( upper_edge_inside )
            {
              // New ROI's upper edge is inside the existing ROI - allow abutting, but ensure
              //  the nearest source gamma still has at least 1 FWHM of extent above it.
              const double nearest_gamma = nearest_source_gamma_in_range(
                roi.lower_energy, existing.lower_energy, false );
              if( nearest_gamma > 0.0 )
              {
                const bool have_fwhm_range = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0)
                                             && (fwhm_lower_energy < fwhm_upper_energy);
                const double fwhm_eval = have_fwhm_range
                  ? std::clamp( nearest_gamma, fwhm_lower_energy, fwhm_upper_energy ) : nearest_gamma;
                const double gamma_fwhm = DetectorPeakResponse::peakResolutionFWHM(
                  static_cast<float>(fwhm_eval), fwhm_form, fwhm_coefficients );
                const double max_upper = nearest_gamma + std::max( gamma_fwhm, 1.0 );
                roi.upper_energy = std::min( existing.lower_energy, max_upper );
              }else
              {
                roi.upper_energy = existing.lower_energy;
              }
            }
          }//if( overlap )

          // Enforce FWHM-based margin from existing peak means.
          // New ROI edges must not be closer than 0.5*auto_rel_eff_sol_min_fwhm_roi * FWHM
          // to any peak mean in the existing ROI.
          for( const double peak_mean : existing.peak_means )
          {
            if( (peak_mean <= roi.lower_energy) || (peak_mean >= roi.upper_energy) )
              continue;  // peak mean is outside this (possibly trimmed) ROI - no issue

            // Peak mean is inside the new ROI - this means the new ROI extends past
            // an existing peak mean. Pull back whichever edge is closer.
            const bool have_fwhm_range = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0)
                                         && (fwhm_lower_energy < fwhm_upper_energy);
            const double fwhm_eval = have_fwhm_range
              ? std::clamp( peak_mean, fwhm_lower_energy, fwhm_upper_energy ) : peak_mean;
            const double peak_fwhm = DetectorPeakResponse::peakResolutionFWHM(
              static_cast<float>(fwhm_eval), fwhm_form, fwhm_coefficients );
            const double min_margin = 0.5 * config.auto_rel_eff_sol_min_fwhm_roi * peak_fwhm;

            const double dist_from_lower = peak_mean - roi.lower_energy;
            const double dist_from_upper = roi.upper_energy - peak_mean;

            if( dist_from_lower < dist_from_upper )
            {
              // Peak mean is closer to lower edge - pull lower edge away
              roi.lower_energy = peak_mean + min_margin;
            }
            else
            {
              // Peak mean is closer to upper edge - pull upper edge away
              roi.upper_energy = peak_mean - min_margin;
            }
          }//for( peak means )
        }//for( const ExistingRoiInfo &existing : existing_roi_ranges )

        // Skip if ROI was marked for removal
        if( (roi.lower_energy < 0.0) || (roi.upper_energy < 0.0) )
          continue;

        // Skip if ROI bounds are invalid after trimming
        if( roi.lower_energy >= roi.upper_energy )
          continue;

        // Check if ROI is large enough (at least ~5 channels)
        const size_t lower_channel = energy_cal->channel_for_energy( roi.lower_energy );
        const size_t upper_channel = energy_cal->channel_for_energy( roi.upper_energy );
        if( (upper_channel - lower_channel) < min_roi_channels )
          continue;

        // Check if any source/NORM gamma is still present in the reduced ROI
        bool has_source_gamma = false;

        // Check source gammas
        for( size_t src_index = 0; !has_source_gamma && (src_index < sources.size()); ++src_index )
        {
          const RelActCalcAuto::NucInputInfo &src = sources[src_index];
          assert( !RelActCalcAuto::is_null(src.source) );

          if( RelActCalcAuto::is_null( src.source ) )
            continue;

          const vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( src.source, 1.0, src.age );
          for( size_t photon_index = 0; !has_source_gamma && (photon_index < photons.size()); ++photon_index )
          {
            const SandiaDecay::EnergyRatePair &photon = photons[photon_index];
            has_source_gamma = ( (photon.energy >= roi.lower_energy)
                                && (photon.energy <= roi.upper_energy)
                                && (photon.numPerSecond > std::numeric_limits<float>::epsilon()) );
          }
        }//for( const RelActCalcAuto::NucInputInfo &src : sources )

        // Check NORM gammas if fitting NORM peaks
        if( !has_source_gamma && fit_norm_peaks )
        {
          std::array<const SandiaDecay::Nuclide *,5> norm_nucs{
            db->nuclide("U238"), db->nuclide("Ra226"), db->nuclide("U235"),
            db->nuclide("Th232"), db->nuclide("K40")
          };

          for( size_t norm_index = 0; norm_index < norm_nucs.size(); ++norm_index )
          {
            const SandiaDecay::Nuclide *norm_nuc = norm_nucs[norm_index];
            assert( norm_nuc );
            if( !norm_nuc )
              continue;

            const std::vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( norm_nuc, 1.0, 0.0 );
            for( size_t photon_index = 0; !has_source_gamma && (photon_index < photons.size()); ++photon_index )
            {
              const SandiaDecay::EnergyRatePair &photon = photons[photon_index];
              has_source_gamma = ((photon.energy >= roi.lower_energy)
                                  && (photon.energy <= roi.upper_energy)
                                  && (photon.numPerSecond > std::numeric_limits<float>::epsilon()) );
            }
          }//for( size_t norm_index = 0; norm_index < norm_nucs.size(); ++norm_index )
        }//if( !has_source_gamma && fit_norm_peaks )

        // Only keep ROI if it still contains a source/NORM gamma
        if( has_source_gamma )
          filtered_rois.push_back( roi );
      }//for( RelActCalcAuto::RoiRange roi : options.rois )

      options.rois = filtered_rois;

#if( PERFORM_DEVELOPER_CHECKS )
      // Verify no filtered ROI overlaps with any existing other-source ROI
      for( const RelActCalcAuto::RoiRange &roi : filtered_rois )
      {
        for( const ExistingRoiInfo &existing : existing_roi_ranges )
        {
          const bool overlaps = (roi.lower_energy < existing.upper_energy) && (roi.upper_energy > existing.lower_energy);
          if( overlaps )
          {
            std::cerr << "Existing ROI trimming: filtered ROI [" << roi.lower_energy << ", " << roi.upper_energy
                      << "] still overlaps existing ROI [" << existing.lower_energy << ", " << existing.upper_energy << "]" << std::endl;
            assert( !overlaps );
          }
        }
      }
#endif
    }//if( !existing_roi_ranges.empty() )
  }// if( !existing_peaks_as_free && !user_peaks.empty() )

  // Lambda to filter a set of ROIs against existing_roi_ranges (populated in DoNotUseExistingRois
  // mode with all user ROIs, or in default mode with other-source-only ROIs).  Used for each
  // iteration of the refinement loop, so that refined ROIs also avoid existing user ROIs.
  const auto filter_rois_for_existing = [&](
      std::vector<RelActCalcAuto::RoiRange> rois,
      const std::shared_ptr<const SpecUtils::Measurement> &working_foreground )
    -> std::vector<RelActCalcAuto::RoiRange>
  {
    if( existing_roi_ranges.empty() )
      return rois;

    std::vector<ExistingRoiInfo> working_existing_roi_ranges = existing_roi_ranges;
    const std::shared_ptr<const SpecUtils::EnergyCalibration> original_cal
      = orig_foreground->energy_calibration();
    const std::shared_ptr<const SpecUtils::EnergyCalibration> working_cal
      = working_foreground ? working_foreground->energy_calibration() : original_cal;
    if( original_cal && working_cal && (original_cal != working_cal) )
    {
      for( ExistingRoiInfo &existing : working_existing_roi_ranges )
      {
        existing.lower_energy = working_cal->energy_for_channel(
            original_cal->channel_for_energy(existing.lower_energy) );
        existing.upper_energy = working_cal->energy_for_channel(
            original_cal->channel_for_energy(existing.upper_energy) );
        for( double &peak_mean : existing.peak_means )
          peak_mean = working_cal->energy_for_channel(
              original_cal->channel_for_energy(peak_mean) );
      }
    }

    const double overlap_buffer_kev = 1.0;
    const double min_roi_channels = 5.0;
    const std::shared_ptr<const SpecUtils::EnergyCalibration> energy_cal = working_cal;

    std::vector<RelActCalcAuto::RoiRange> result_rois;
    result_rois.reserve( rois.size() );

    for( RelActCalcAuto::RoiRange roi : rois )
    {
      for( const ExistingRoiInfo &existing : working_existing_roi_ranges )
      {
        if( !((roi.lower_energy < existing.upper_energy) && (roi.upper_energy > existing.lower_energy)) )
          continue;

        const bool new_contained_by_existing = (roi.lower_energy >= existing.lower_energy) && (roi.upper_energy <= existing.upper_energy);
        const bool new_contains_existing = (roi.lower_energy < existing.lower_energy) && (roi.upper_energy > existing.upper_energy);
        const bool lower_edge_inside = (roi.lower_energy >= existing.lower_energy) && (roi.lower_energy < existing.upper_energy);
        const bool upper_edge_inside = (roi.upper_energy > existing.lower_energy) && (roi.upper_energy <= existing.upper_energy);

        if( new_contained_by_existing )
        {
          roi.lower_energy = -1.0;
          roi.upper_energy = -1.0;
          break;
        }else if( new_contains_existing )
        {
          const double left_width = existing.lower_energy - roi.lower_energy;
          const double right_width = roi.upper_energy - existing.upper_energy;
          if( left_width > right_width )
            roi.upper_energy = existing.lower_energy - overlap_buffer_kev;
          else
            roi.lower_energy = existing.upper_energy + overlap_buffer_kev;
        }else if( lower_edge_inside )
        {
          roi.lower_energy = existing.upper_energy + overlap_buffer_kev;
        }else if( upper_edge_inside )
        {
          roi.upper_energy = existing.lower_energy - overlap_buffer_kev;
        }

        // Enforce FWHM-based margin from existing peak means
        for( const double peak_mean : existing.peak_means )
        {
          if( (peak_mean <= roi.lower_energy) || (peak_mean >= roi.upper_energy) )
            continue;

          const bool have_fwhm_range = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0)
                                       && (fwhm_lower_energy < fwhm_upper_energy);
          const double fwhm_eval = have_fwhm_range
            ? std::clamp( peak_mean, fwhm_lower_energy, fwhm_upper_energy ) : peak_mean;
          const double fwhm_at_peak = DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>( fwhm_eval ),
            fwhm_form, fwhm_coefficients );
          const double min_margin = 0.5 * config.auto_rel_eff_sol_min_fwhm_roi * fwhm_at_peak;

          if( (peak_mean - roi.lower_energy) < (roi.upper_energy - peak_mean) )
            roi.lower_energy = peak_mean + min_margin;
          else
            roi.upper_energy = peak_mean - min_margin;
        }
      }//for( existing_roi_ranges )

      if( (roi.lower_energy < 0.0) || (roi.upper_energy < 0.0) )
        continue;
      if( roi.lower_energy >= roi.upper_energy )
        continue;

      if( energy_cal && energy_cal->valid() )
      {
        const size_t lower_ch = energy_cal->channel_for_energy( roi.lower_energy );
        const size_t upper_ch = energy_cal->channel_for_energy( roi.upper_energy );
        if( (upper_ch - lower_ch) < static_cast<size_t>(min_roi_channels) )
          continue;
      }

      // Only keep if it still contains a source gamma
      bool has_source_gamma = false;
      for( size_t src_idx = 0; !has_source_gamma && (src_idx < sources.size()); ++src_idx )
      {
        const vector<SandiaDecay::EnergyRatePair> photons
          = get_source_photons( sources[src_idx].source, 1.0, sources[src_idx].age );
        for( const SandiaDecay::EnergyRatePair &ph : photons )
        {
          if( (ph.energy >= roi.lower_energy) && (ph.energy <= roi.upper_energy)
             && (ph.numPerSecond > std::numeric_limits<float>::epsilon()) )
          {
            has_source_gamma = true;
            break;
          }
        }
      }

      if( !has_source_gamma && fit_norm_peaks )
      {
        const std::array<const SandiaDecay::Nuclide *,5> norm_nucs{
          db->nuclide("U238"), db->nuclide("Ra226"), db->nuclide("U235"),
          db->nuclide("Th232"), db->nuclide("K40")
        };
        for( const SandiaDecay::Nuclide *norm_nuc : norm_nucs )
        {
          if( !norm_nuc )
            continue;
          const vector<SandiaDecay::EnergyRatePair> photons = get_source_photons( norm_nuc, 1.0, 0.0 );
          for( const SandiaDecay::EnergyRatePair &ph : photons )
          {
            if( (ph.energy >= roi.lower_energy) && (ph.energy <= roi.upper_energy)
               && (ph.numPerSecond > std::numeric_limits<float>::epsilon()) )
            {
              has_source_gamma = true;
              break;
            }
          }
          if( has_source_gamma )
            break;
        }
      }

      // A transactionally accepted R6 nuisance line is part of the fitted model even though it is
      // not a requested source and will be hidden from the public peak vectors.  Preserve its ROI
      // during later existing-ROI filtering; the initial R6 transaction has already rejected any
      // new coverage window that conflicts with a protected user ROI.
      if( !has_source_gamma )
      {
        for( const double energy : auto_interferer_lines )
        {
          if( (energy >= roi.lower_energy) && (energy <= roi.upper_energy) )
          {
            has_source_gamma = true;
            break;
          }
        }
      }

      if( has_source_gamma )
        result_rois.push_back( roi );
    }//for( roi : rois )

    return result_rois;
  };//filter_rois_for_existing lambda


  // Tracks which input user_peaks will be replaced when the result is accepted.
  // Populated only when ExistingPeaksAsFreePeak is set; each entry is {original peak ptr, energy
  // of the floating peak it maps to}.  After the fit, this is used to populate
  // result.original_peaks_to_remove based on whether the ROI ended up in the solution.
  std::vector<std::pair<std::shared_ptr<const PeakDef>, double>> existing_peaks_added_as_floating;

  if( existing_peaks_as_free && !user_peaks.empty() )
  {
    assert( !do_not_use_existing_rois );
    
    // Treat a shared user continuum atomically.  A policy boundary may touch only one peak in a
    // mixed ROI; enrolling that peak alone would resize the continuum out from under its sibling.
    std::map<const PeakContinuum *, std::vector<std::shared_ptr<const PeakDef>>> roi_groups;
    for( const std::shared_ptr<const PeakDef> &peak : user_peaks )
    {
      assert( peak );
      if( peak )
        roi_groups[peak->continuum().get()].push_back( peak );
    }

    for( const std::pair<const PeakContinuum *const,
             std::vector<std::shared_ptr<const PeakDef>>> &roi_group : roi_groups )
    {
      const std::vector<std::shared_ptr<const PeakDef>> &peaks = roi_group.second;
      assert( !peaks.empty() );
      const bool touched = std::any_of( std::begin(peaks), std::end(peaks),
        [&options]( const std::shared_ptr<const PeakDef> &peak ) {
          return std::any_of( std::begin(options.rois), std::end(options.rois),
            [&peak]( const RelActCalcAuto::RoiRange &roi ) {
              return (peak->mean() >= roi.lower_energy) && (peak->mean() <= roi.upper_energy);
            } );
        } );

      if( !touched )
      {
        ExistingRoiInfo info;
        info.lower_energy = peaks.front()->lowerX();
        info.upper_energy = peaks.front()->upperX();
        for( const std::shared_ptr<const PeakDef> &peak : peaks )
          info.peak_means.push_back( peak->mean() );
        existing_roi_ranges.push_back( info );
        continue;
      }

      bool has_bystander = false;
      for( const std::shared_ptr<const PeakDef> &peak : peaks )
      {
        const double peak_energy = peak->mean();
        if( !peak_source_is_in_fit(peak) )
        {
          has_bystander = true;
          RelActCalcAuto::FloatingPeak fp;
          fp.energy = peak_energy;
          fp.release_fwhm = false;
          // Existing bystanders use their observed mean as a Known reference energy so they track
          // the fitted energy-cal adjustment with source peaks instead of forming a shifted doublet.
          fp.energy_origin = RelActCalcAuto::FloatingPeak::EnergyType::Known;
          options.floating_peaks.push_back( fp );
          carried_bystander_float_energies.push_back( peak_energy );
        }
        existing_peaks_added_as_floating.emplace_back( peak, peak_energy );
      }

      if( has_bystander )
      {
        RelActCalcAuto::RoiRange protected_roi;
        protected_roi.lower_energy = peaks.front()->lowerX();
        protected_roi.upper_energy = peaks.front()->upperX();
        protected_roi.continuum_type = peaks.front()->continuum()->type();
        protected_roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
        options.rois.erase( std::remove_if( std::begin(options.rois), std::end(options.rois),
          [&protected_roi]( const RelActCalcAuto::RoiRange &roi ) {
            return (roi.lower_energy < protected_roi.upper_energy)
                && (roi.upper_energy > protected_roi.lower_energy);
          } ), std::end(options.rois) );
        options.rois.push_back( protected_roi );
        protected_mixed_rois.push_back( protected_roi );
      }
    }//for( roi_groups )

    std::sort( existing_roi_ranges.begin(), existing_roi_ranges.end(),
      []( const ExistingRoiInfo &a, const ExistingRoiInfo &b ) {
        return a.lower_energy < b.lower_energy;
      } );

    if( !existing_roi_ranges.empty() )
      options.rois = filter_rois_for_existing( options.rois, orig_foreground );
  }// if( existing_peaks_as_free && !user_peaks.empty() )


  if( fit_norm_peaks )
  {
    try
    {
      vector<RelActCalcAuto::NucInputInfo> norm_sources = get_norm_sources( sources, config.norm_css_color );

      // Self-attenuating "soil" shielding for the NORM Physical-Model rel-eff curve.
      auto self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>();
      //const char *soil_chem_formula = "H0.022019C0.009009O0.593577Al0.066067Si0.272289K0.01001Fe0.027029 d=1.6";
      //self_atten->material = make_share<Material>( MaterialDB::materialFromChemicalFormula( soil_chem_formula, db ) );
      self_atten->atomic_number = 10.4;
      self_atten->areal_density = (1.6 * PhysicalUnits::g / PhysicalUnits::cm3) * (100 * PhysicalUnits::cm);  //Transmission frac @2614 keV, though 100 cm soil is 0.2%
      self_atten->fit_atomic_number = false;
      self_atten->fit_areal_density = false;

      // NOTE: an external attenuator (~1 mm concrete, fit areal density) used to be configured
      // here, but only ever fed a since-removed (dead) manual NORM pre-solve - it was never
      // applied to norm_rel_eff_curve below.  Add it to phys_model_external_atten there if the
      // NORM curve should model external attenuation.

      // Note: input_rois already contains NORM ROIs (merged by the caller in fit_peaks_for_nuclides).
      // We rebuild the InitialRoi metadata from the CURRENT options.rois - NOT from the raw
      // input_rois - so the existing-ROI trimming and mixed-ROI fixed bounds applied above are
      // preserved (rebuilding from input_rois silently discarded them).
      vector<InitialRoi> initial_src_norm_info_rois;
      for( const RelActCalcAuto::RoiRange &roi : options.rois )
      {
        const bool protected_geometry = std::any_of(
            std::begin(protected_mixed_rois), std::end(protected_mixed_rois),
            [&roi]( const RelActCalcAuto::RoiRange &protected_roi ) {
              return (std::fabs(roi.lower_energy - protected_roi.lower_energy) < 1.0e-6)
                  && (std::fabs(roi.upper_energy - protected_roi.upper_energy) < 1.0e-6);
            } );
        if( protected_geometry )
          continue;  // reinsert unchanged after automatic source/NORM reconciliation
        InitialRoi roi_info;
        roi_info.roi = roi;
        roi_info.center_energy = 0.5*(roi.upper_energy + roi.lower_energy);
        const bool have_fwhm_range = ((fwhm_lower_energy > 0.0)
                                      && (fwhm_upper_energy > 0.0)
                                      && (fwhm_lower_energy < fwhm_upper_energy));
        const double fwhm_eval_energy = have_fwhm_range
            ? std::clamp( roi_info.center_energy, fwhm_lower_energy, fwhm_upper_energy )
            : roi_info.center_energy;
        roi_info.fwhm = DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(fwhm_eval_energy), fwhm_form, fwhm_coefficients );

        float min_sigma_width, max_sigma_width;
        expected_peak_width_limits( roi_info.fwhm,
          det_type,
          orig_foreground, min_sigma_width, max_sigma_width );

        if( roi_info.fwhm < (min_sigma_width*PhysicalUnits::fwhm_nsigma) )
          roi_info.fwhm = min_sigma_width*PhysicalUnits::fwhm_nsigma;
        if( roi_info.fwhm > (max_sigma_width*PhysicalUnits::fwhm_nsigma) )
          roi_info.fwhm = max_sigma_width*PhysicalUnits::fwhm_nsigma;

        // Find the largest auto_search_peak within the ROI for the amplitude estimate; when one
        // exists, also anchor center_energy on its mean rather than the ROI midpoint - the
        // clean-gap merge test and split-point constraints model the ROI's signal as a Gaussian
        // at center_energy, and a real peak position grounds that far better than the geometric
        // midpoint of a (possibly asymmetric) ROI.
        for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
        {
          if( (peak->mean() >= roi_info.roi.lower_energy)
            && (peak->mean() <= roi_info.roi.upper_energy)
            && (peak->amplitude() > roi_info.estimated_amplitude) )
          {
            roi_info.estimated_amplitude = peak->amplitude();
            roi_info.center_energy = peak->mean();
          }
        }

        for( const std::pair<double,double> &candidate : initial_modeled_peak_candidates )
        {
          if( (candidate.first >= roi.lower_energy) && (candidate.first <= roi.upper_energy) )
          {
            roi_info.modeled_energies.push_back( candidate.first );
            roi_info.modeled_areas.push_back( candidate.second );
          }
        }
        if( !roi_info.modeled_energies.empty() )
        {
          const std::vector<double>::const_iterator dominant = std::max_element(
              std::begin(roi_info.modeled_areas), std::end(roi_info.modeled_areas) );
          const size_t dominant_index = static_cast<size_t>(
              dominant - std::begin(roi_info.modeled_areas) );
          roi_info.center_energy = roi_info.modeled_energies[dominant_index];
          roi_info.estimated_amplitude = *dominant;
        }

        initial_src_norm_info_rois.push_back( roi_info );
      }//for( const RelActCalcAuto::RoiRange &roi : options.rois )

      options.rois = merge_rois( initial_src_norm_info_rois, config, {}, orig_foreground,
                                 global_cont.valid() ? &global_cont : nullptr,
                                 &result.automatic_roi_diagnostics, "source/NORM merge" );
      options.rois.insert( std::end(options.rois), std::begin(protected_mixed_rois),
                           std::end(protected_mixed_rois) );

      RelActCalcAuto::RelEffCurveInput norm_rel_eff_curve;
      norm_rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
      norm_rel_eff_curve.rel_eff_eqn_order = 0;
      norm_rel_eff_curve.nucs_of_el_same_age = false;
      norm_rel_eff_curve.phys_model_corr.corr_fcn = RelActCalc::PhysModelCorrFcn::None;
      norm_rel_eff_curve.phys_model_self_atten = self_atten;
      norm_rel_eff_curve.nuclides = norm_sources;
      norm_rel_eff_curve.name = "NORM curve";

      if( !norm_sources.empty() )
      {
        norm_rel_eff_index = static_cast<int>( options.rel_eff_curves.size() );
        options.rel_eff_curves.push_back( norm_rel_eff_curve );
      }
    }catch( const std::exception &e )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
      result.error_message = "Error performing initial estimation of background peaks: " + string(e.what());
    }
  }//if( fit_norm_peaks )

  // R6 detection is deliberately side-effect free here.  The requested-source model is solved
  // first; only a successful, non-empty incumbent may later be augmented transactionally with the
  // bounded nuisance curve.  This prevents a nuisance ROI from turning an honestly empty source
  // request into a non-empty fit, and gives every augmented failure a safe source-only fallback.
  if( !fit_norm_peaks && !disable_auto_interferer_fit )
  {
    try
    {
      const auto fwhm_at_energy = [&]( double e ) -> double {
        const bool have = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0)
                          && (fwhm_lower_energy < fwhm_upper_energy);
        const double ee = have ? std::clamp( e, fwhm_lower_energy, fwhm_upper_energy ) : e;
        return DetectorPeakResponse::peakResolutionFWHM( static_cast<float>(ee), fwhm_form, fwhm_coefficients );
      };

      // Expand requested sources to their in-range gammas (mirrors add_floating_511's expansion).
      std::vector<detail::RequestedSourceGammas> source_gammas;
      for( const RelActCalcAuto::NucInputInfo &src : sources )
      {
        if( RelActCalcAuto::is_null( src.source ) )
          continue;
        detail::RequestedSourceGammas sg;
        sg.source = src.source;
        std::vector<SandiaDecay::EnergyRatePair> photons;
        try{
          photons = get_source_photons( src.source, 1.0, get_source_age( src.source, src.age ) );
        }catch( const std::exception & ){
          continue;
        }
        for( const SandiaDecay::EnergyRatePair &p : photons )
        {
          if( (p.energy >= min_valid_energy) && (p.energy <= max_valid_energy) )
          {
            sg.energies.push_back( p.energy );
            sg.yields.push_back( p.numPerSecond );
          }
        }
        if( !sg.energies.empty() )
          source_gammas.push_back( std::move(sg) );
      }

      std::vector<std::string> interferer_warnings;
      interferer_candidates = detail::find_strong_unmodeled_interferers(
            source_gammas, auto_search_peaks, fwhm_at_energy, fit_norm_peaks,
            min_valid_energy, max_valid_energy, orig_background, drf, peak_fit_prefs,
            &interferer_warnings, global_cont.valid() ? &global_cont : nullptr,
            &interferer_guard_energies );

      for( const std::string &w : interferer_warnings )
        result.warnings.push_back( w );
    }catch( const std::exception &e )
    {
      if( should_debug_print() )
        std::cerr << "R6 interferer detection failed (continuing source-only): " << e.what() << std::endl;
    }
  }//if( !fit_norm_peaks && !disable_auto_interferer_fit ) -- R6 interferer detection

  // Add a floating peak at 511 keV if appropriate (see function documentation for physics reasoning)
  add_floating_511_peak_if_appropriate( options, sources, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy );

  // Add floating peaks for escape peaks of high-energy gammas if appropriate
  add_escape_peak_floating_peaks_if_appropriate( options, auto_search_peaks, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy, config );

  // Compute auto-search peaks that don't correspond to user peaks -- these are peaks in the
  // auto-search that may interfere with our source/NORM ROIs.  Used during the iterative
  // refinement loop to shrink ROIs away from interfering peaks.
  const std::vector<std::shared_ptr<const PeakDef>> unfit_auto_peaks
    = compute_unfit_auto_peaks( auto_search_peaks, user_peaks );

  try
  {
    // Resolve any overlapping ROIs (may occur from escape peak ROIs)
    resolve_automatic_overlapping_rois( options.rois, options.floating_peaks,
        orig_foreground, global_cont.valid() ? &global_cont : nullptr, global_fwhm_at,
        unfit_auto_peaks, config, "initial escape/overlap finalization",
        &result.automatic_roi_diagnostics, protected_mixed_rois,
        initial_modeled_peak_candidates,
        use_automatic_roi_policy );
    remove_floating_peaks_without_roi( options );

    // Change 2: a genuinely-empty ROI set is a valid empty result, NOT a setup error.  With the
    // keep-gate rescue removed and found+matched peaks seeded upstream, reaching this point with no
    // ROIs means every predicted source line was sub-threshold AND the search found no confirming
    // peak - i.e. there is honestly nothing to fit.  Return Success with zero peaks and a warning
    // rather than letting RelActCalcAuto::solve throw "No ROIs are defined" (which the caller would
    // report as FailedToSetupProblem, conflating "nothing to fit" with "we dropped everything").
    // The throw in solve() is left in place as a guard for other callers.  [architecture review 2026-07-18]
    if( options.rois.empty() )
    {
      if( should_debug_print() )
        std::cout << "fit_peaks_for_nuclide_relactauto: no ROIs to fit - returning valid empty"
                     " (Success, 0 peaks)" << std::endl;
      result.status = RelActCalcAuto::RelActAutoSolution::Status::Success;
      result.solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::Success;
      result.warnings.push_back( "No significant peaks were found for the requested source(s)." );
      return result;
    }//if( options.rois.empty() )

    // Call RelActAuto::solve with provided options
    RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
      options, orig_foreground, orig_background, drf, auto_search_peaks, det_type
    );

    const std::vector<RelActCalcManual::GenericPeakInfo> significant_source_anchors
      = use_automatic_roi_policy
        ? distinct_significant_source_anchors(
            source_anchor_candidates, config.roi_significance_z )
        : std::vector<RelActCalcManual::GenericPeakInfo>();
    const size_t initial_source_anchors_preserved
      = (use_automatic_roi_policy && (sources_rel_eff_index >= 0)
          && RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status))
        ? preserved_source_anchor_count( solution,
            static_cast<size_t>(sources_rel_eff_index), significant_source_anchors,
            orig_foreground->live_time() )
        : significant_source_anchors.size();
    const bool source_anchor_collapse_detected = use_automatic_roi_policy
        && detail::should_try_source_clean_recovery(
            significant_source_anchors.size(), initial_source_anchors_preserved );
    std::vector<RelActCalcManual::GenericPeakInfo> recovered_source_anchors;

    // The manual matcher can associate a real found peak with the requested source even when the
    // fitted source curve says that source explains only a small fraction of its area.  Those
    // contaminant-like matches used to receive unconditional tight ROIs, and a dense set of them
    // could pull the first empirical curve into the collapsed solution.  Retry transactionally with
    // only the manual stage's source-accounted-for seeds.  Ordinary spectra retain the incumbent;
    // the challenger is considered only after the source-evidence collapse gate fires.
    if( source_anchor_collapse_detected && !clean_source_rois.empty()
        && !fit_norm_peaks && user_peaks.empty()
        && (sources_rel_eff_index >= 0) )
    {
      const size_t diagnostic_checkpoint = result.automatic_roi_diagnostics.size();
      bool retained_challenger = false;
      try
      {
        RelActCalcAuto::Options clean_options = options;
        clean_options.rois = clean_source_rois;
        if( has_provisional_fallback_source_anchors
            && (sources_rel_eff_index >= 0) )
        {
          // With only pre-manual, data-supported anchors there is no defensible basis for the
          // configured high-order curve.  Use its order-zero member as a conservative challenger;
          // it can be retained only through the same source-anchor and common-channel score gates.
          clean_options.rel_eff_curves[sources_rel_eff_index].rel_eff_eqn_order = 0;
        }
        add_floating_511_peak_if_appropriate( clean_options, sources, fit_norm_peaks,
            det_type, min_valid_energy, max_valid_energy );
        add_escape_peak_floating_peaks_if_appropriate( clean_options, auto_search_peaks,
            fit_norm_peaks, det_type, min_valid_energy, max_valid_energy, config );
        resolve_automatic_overlapping_rois( clean_options.rois,
            clean_options.floating_peaks, orig_foreground,
            global_cont.valid() ? &global_cont : nullptr, global_fwhm_at,
            unfit_auto_peaks, config, "source-clean seed challenger finalization",
            &result.automatic_roi_diagnostics, protected_mixed_rois,
            initial_modeled_peak_candidates, use_automatic_roi_policy );
        remove_floating_peaks_without_roi( clean_options );

        RelActCalcAuto::RelActAutoSolution clean_solution = RelActCalcAuto::solve(
            clean_options, orig_foreground, orig_background, drf,
            auto_search_peaks, det_type );
        const size_t clean_preserved = preserved_source_anchor_count(
            clean_solution, static_cast<size_t>(sources_rel_eff_index),
            significant_source_anchors, orig_foreground->live_time() );
        const double incumbent_score = source_anchor_data_aicc(
            solution, significant_source_anchors, config.manual_releff_aicc_penalty );
        const double clean_score = RelActCalcAuto::RelActAutoSolution::is_usable_status(clean_solution.m_status)
            ? source_anchor_data_aicc( clean_solution, significant_source_anchors,
                config.manual_releff_aicc_penalty )
            : std::numeric_limits<double>::max();
        const size_t incumbent_fit_anchors
          = significant_requested_source_anchor_count( solution, sources,
              significant_source_anchors, config.roi_significance_z, global_fwhm_at );
        const size_t clean_fit_anchors
          = significant_requested_source_anchor_count( clean_solution, sources,
              significant_source_anchors, config.roi_significance_z, global_fwhm_at );
        const bool accept = detail::should_accept_source_clean_challenger(
            RelActCalcAuto::RelActAutoSolution::is_usable_status(clean_solution.m_status),
            initial_source_anchors_preserved, clean_preserved,
            incumbent_fit_anchors, clean_fit_anchors, incumbent_score, clean_score );
        if( should_debug_print() )
        {
          std::cerr << "Source-clean seed challenger: anchors="
                    << significant_source_anchors.size() << ", preserved="
                    << initial_source_anchors_preserved << "->" << clean_preserved
                    << ", fitted anchors=" << incumbent_fit_anchors
                    << "->" << clean_fit_anchors
                    << ", common-anchor data AICc=" << incumbent_score << "->" << clean_score
                    << ", accepted=" << accept << std::endl;
        }
        if( accept )
        {
          options = std::move( clean_options );
          solution = std::move( clean_solution );
          recovered_source_anchors = preserved_source_anchors( solution,
              static_cast<size_t>(sources_rel_eff_index), significant_source_anchors,
              orig_foreground->live_time() );
          retained_challenger = true;
          result.warnings.push_back( "Replaced contaminant-like found-peak seeds after they"
            " collapsed multiple significant requested-source anchors." );
        }
      }catch( const std::exception &error )
      {
        result.warnings.push_back( "The source-clean seed challenger failed; retained the"
          " successful incumbent fit: " + std::string(error.what()) );
      }
      if( !retained_challenger )
        result.automatic_roi_diagnostics.resize( diagnostic_checkpoint );
    }

    // As of 20260103, energy calibration adjustments may cause failure to fit the correct solution sometimes,
    //  so if our current solution failed, or is really bad, we'll try without fitting energy cal
    if( ( options.energy_cal_type != RelActCalcAuto::EnergyCalFitType::NoFit )
      && ((!RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status))
        || (reduced_chi2(solution) > 10.0)) ) //10.0 arbitrary - and un-explored
    {
      RelActCalcAuto::Options no_ecal_opts = options;
      no_ecal_opts.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;

      add_floating_511_peak_if_appropriate( no_ecal_opts, sources, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy );
      add_escape_peak_floating_peaks_if_appropriate( no_ecal_opts, auto_search_peaks, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy, config );
      resolve_automatic_overlapping_rois( no_ecal_opts.rois, no_ecal_opts.floating_peaks,
          orig_foreground, global_cont.valid() ? &global_cont : nullptr, global_fwhm_at,
          unfit_auto_peaks, config, "no-ecal escape/overlap finalization",
          &result.automatic_roi_diagnostics, protected_mixed_rois,
          initial_modeled_peak_candidates,
          use_automatic_roi_policy );
      remove_floating_peaks_without_roi( no_ecal_opts );

      RelActCalcAuto::RelActAutoSolution no_ecal_solution = RelActCalcAuto::solve(
        no_ecal_opts, orig_foreground, orig_background, drf, auto_search_peaks, det_type
      );

      // If the solution is still really bad - we'll try a Physical Model solution
      // Optionally with external shielding if configured.
      // NOTE: gate on the just-computed no-ecal solution (was erroneously testing the
      // superseded `solution`), so escalation reflects the fit we are about to keep.
      // The no-ecal solution is kept in its own variable so a successful-but-mediocre no-ecal
      // fit is not lost when the desperation solve fails or is worse; the final selection below
      // takes the best of {original, no-ecal, desperation}.
      RelActCalcAuto::RelActAutoSolution desperation_solution;
      desperation_solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::NotInitiated;

      if( ((!RelActCalcAuto::RelActAutoSolution::is_usable_status(no_ecal_solution.m_status))
          || (reduced_chi2(no_ecal_solution) > 10.0))
         && (sources_rel_eff_index >= 0) )
      {
        // Base the desperation solve on the no-ecal options: the retry exists because energy-cal
        // fitting may be destabilizing the solution, so the desperation attempt must not quietly
        // re-enable it (it previously inherited energy_cal_type = NonLinearFit from `options`).
        RelActCalcAuto::Options desperation_opts = options;
        desperation_opts.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
        RelActCalcAuto::RelEffCurveInput &curve = desperation_opts.rel_eff_curves[sources_rel_eff_index];
        curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
        curve.rel_eff_eqn_order = 0;
        curve.phys_model_self_atten = nullptr;
        curve.phys_model_external_atten.clear();
        
        // Apply external shielding if conditions are met
        if( should_use_desperation_shielding( config.desperation_phys_model_atomic_number, options.rois ) )
        {
          try
          {
            std::shared_ptr<RelActCalc::PhysicalModelShieldInput> shield
            = create_desperation_shielding( config.desperation_phys_model_atomic_number,
                                           config.desperation_phys_model_areal_density_g_per_cm2 );
            curve.phys_model_external_atten.push_back( shield );
            
            if( should_debug_print() )
            {
              std::cerr << "First desperation attempt: using external shielding with AN="
              << config.desperation_phys_model_atomic_number
              << ", starting AD=" << config.desperation_phys_model_areal_density_g_per_cm2
              << " g/cm2" << std::endl;
            }
          }
          catch( const std::exception &e )
          {
            if( should_debug_print() )
              std::cerr << "Failed to create desperation shielding: " << e.what() << std::endl;
            // Continue without shielding
          }
        }
        else
        {
          if( should_debug_print() )
            std::cerr << "First desperation attempt: not using external shielding" << std::endl;
        }
        
        curve.phys_model_corr.corr_fcn = (options.rois.size() > 2)
                            ? RelActCalc::PhysModelCorrFcn::Hoerl : RelActCalc::PhysModelCorrFcn::None;
        
        add_floating_511_peak_if_appropriate( desperation_opts, sources, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy );
        add_escape_peak_floating_peaks_if_appropriate( desperation_opts, auto_search_peaks, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy, config );

        resolve_automatic_overlapping_rois( desperation_opts.rois, desperation_opts.floating_peaks,
            orig_foreground, global_cont.valid() ? &global_cont : nullptr, global_fwhm_at,
            unfit_auto_peaks, config, "desperation escape/overlap finalization",
            &result.automatic_roi_diagnostics, protected_mixed_rois,
            initial_modeled_peak_candidates,
            use_automatic_roi_policy );
        remove_floating_peaks_without_roi( desperation_opts );

        desperation_solution = RelActCalcAuto::solve( desperation_opts, orig_foreground, orig_background, drf, auto_search_peaks, det_type );
      }//If( still a bad solution )

      // Keep the best of {original (ecal), no-ecal, desperation}: a successful candidate replaces
      // `solution` when `solution` failed or the candidate has a lower reduced chi2.
      const auto preserves_recovered_evidence
        = [&]( const RelActCalcAuto::RelActAutoSolution &candidate ) -> bool {
        if( recovered_source_anchors.empty() )
          return true;
        if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(candidate.m_status)
            || (sources_rel_eff_index < 0) )
          return false;
        const bool predicted_preserved = std::all_of(
            std::begin(recovered_source_anchors), std::end(recovered_source_anchors),
            [&candidate, sources_rel_eff_index, &orig_foreground](
                const RelActCalcManual::GenericPeakInfo &anchor ) {
              return predicted_source_anchor_counts( candidate,
                  static_cast<size_t>(sources_rel_eff_index), anchor,
                  orig_foreground->live_time() ) > sm_keep_gate_min_est_counts;
            } );
        return predicted_preserved && significant_requested_source_anchors_preserved(
            candidate, sources, recovered_source_anchors,
            config.roi_significance_z, global_fwhm_at );
      };
      const auto is_better_than_solution
        = [&solution, &preserves_recovered_evidence](
            const RelActCalcAuto::RelActAutoSolution &cand ) -> bool {
        return RelActCalcAuto::RelActAutoSolution::is_usable_status(cand.m_status)
               && preserves_recovered_evidence( cand )
               && ( (!RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status))
                   || (reduced_chi2(solution) > reduced_chi2(cand)) );
      };

      if( is_better_than_solution(no_ecal_solution) )
      {
        if( should_debug_print() )
          std::cerr << "Abandoning fitting e-cal for nuclide" << std::endl;
        solution = std::move( no_ecal_solution );
      }

      if( is_better_than_solution(desperation_solution) )
      {
        if( should_debug_print() )
          std::cerr << "Using Physical-Model desperation solution for nuclide" << std::endl;
        solution = std::move( desperation_solution );
      }
    }

    // Check if initial solve failed
    if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status) )
    {
      result.status = solution.m_status;
      result.error_message = solution.m_error_message;
      return result;
    }

    const auto initial_fwhm_at_energy = [&]( const double energy ) -> double {
      const bool have = (fwhm_lower_energy > 0.0) && (fwhm_upper_energy > 0.0)
                        && (fwhm_lower_energy < fwhm_upper_energy);
      const double eval_energy = have
          ? std::clamp( energy, fwhm_lower_energy, fwhm_upper_energy ) : energy;
      return DetectorPeakResponse::peakResolutionFWHM(
          static_cast<float>(eval_energy), fwhm_form, fwhm_coefficients );
    };
    // All peaks in a successful RelActAuto solution share its fitted resolution model.  Interpolate
    // the resulting fitted peak widths so R6 guards and the R2 rebuild follow the current solution,
    // rather than the pre-solve DRF/auto-search seed coefficients.
    const std::function<double(double)> solution_fwhm_at_energy
      = [&]( const double energy ) -> double {
      const PeakDef *below = nullptr;
      const PeakDef *above = nullptr;
      for( const PeakDef &peak : solution.m_peaks_without_back_sub )
      {
        const double fitted_fwhm = PhysicalUnits::fwhm_nsigma * peak.sigma();
        if( !std::isfinite(fitted_fwhm) || !(fitted_fwhm > 0.0) )
          continue;
        if( (peak.mean() <= energy) && (!below || (peak.mean() > below->mean())) )
          below = &peak;
        if( (peak.mean() >= energy) && (!above || (peak.mean() < above->mean())) )
          above = &peak;
      }
      if( below && above && (above != below)
          && ((above->mean() - below->mean()) > 1.0e-6) )
      {
        const double fraction = (energy - below->mean()) / (above->mean() - below->mean());
        const double below_fwhm = PhysicalUnits::fwhm_nsigma * below->sigma();
        const double above_fwhm = PhysicalUnits::fwhm_nsigma * above->sigma();
        return below_fwhm + fraction * (above_fwhm - below_fwhm);
      }
      if( below )
        return PhysicalUnits::fwhm_nsigma * below->sigma();
      if( above )
        return PhysicalUnits::fwhm_nsigma * above->sigma();
      return initial_fwhm_at_energy( energy );
    };
    const auto is_requested_peak = [&sources]( const PeakDef &peak ) -> bool {
      for( const RelActCalcAuto::NucInputInfo &source : sources )
      {
        if( peak.parentNuclide()
            && (peak.parentNuclide() == RelActCalcAuto::nuclide(source.source)) )
          return true;
        if( peak.xrayElement()
            && (peak.xrayElement() == RelActCalcAuto::element(source.source)) )
          return true;
        if( peak.reaction()
            && (peak.reaction() == RelActCalcAuto::reaction(source.source)) )
          return true;
      }
      return false;
    };
    const auto peak_fit_significance = []( const PeakDef &peak ) -> double {
      const double amplitude = peak.amplitude();
      const double uncertainty = peak.amplitudeUncert();
      return (uncertainty > 0.0) ? (amplitude / uncertainty)
          : ((amplitude > 0.0) ? std::sqrt(amplitude) : 0.0);
    };
    const auto same_gamma_identity = [&]( const PeakDef &lhs, const PeakDef &rhs ) -> bool {
      if( (lhs.parentNuclide() != rhs.parentNuclide())
          || (lhs.xrayElement() != rhs.xrayElement())
          || (lhs.reaction() != rhs.reaction()) )
        return false;
      if( lhs.hasSourceGammaAssigned() && rhs.hasSourceGammaAssigned() )
        return std::fabs(lhs.gammaParticleEnergy() - rhs.gammaParticleEnergy()) < 0.05;
      return std::fabs(lhs.mean() - rhs.mean())
          < (0.25 * solution_fwhm_at_energy(lhs.mean()));
    };
    const auto requested_anchors_preserved = [&]( const RelActCalcAuto::RelActAutoSolution &incumbent,
                                                   const RelActCalcAuto::RelActAutoSolution &challenger,
                                                   const std::vector<double> &affected_energies ) -> bool {
      for( const PeakDef &anchor : incumbent.m_peaks_without_back_sub )
      {
        if( !is_requested_peak(anchor)
            || (peak_fit_significance(anchor) < config.roi_significance_z) )
          continue;

        bool affected = false;
        for( const double energy : affected_energies )
        {
          if( std::fabs(anchor.mean() - energy)
              < (2.0 * solution_fwhm_at_energy(energy)) )
          {
            affected = true;
            break;
          }
        }
        if( affected )
          continue;

        bool found = false;
        for( const PeakDef &candidate : challenger.m_peaks_without_back_sub )
        {
          if( same_gamma_identity(anchor, candidate)
              && (peak_fit_significance(candidate) >= config.roi_significance_z) )
          {
            found = true;
            break;
          }
        }
        if( !found )
          return false;
      }
      return true;
    };
    const auto observable_requested_peaks
      = [&]( const RelActCalcAuto::RelActAutoSolution &candidate ) {
        std::vector<size_t> insignificant_indices;
        compute_filtered_chi2_per_channel(
            candidate, config.roi_significance_z, insignificant_indices );
        const bool have_spectrum_ranges
          = (candidate.m_final_roi_ranges_in_spectrum_cal.size()
              == candidate.m_final_roi_ranges.size());
        const std::vector<RelActCalcAuto::RoiRange> &ranges = have_spectrum_ranges
          ? candidate.m_final_roi_ranges_in_spectrum_cal : candidate.m_final_roi_ranges;

        std::vector<PeakDef> fitted;
        for( const PeakDef &peak : candidate.m_peaks_without_back_sub )
        {
          if( !is_requested_peak(peak) || !peak.continuum() )
            continue;
          const double lower = peak.continuum()->lowerEnergy();
          const double upper = peak.continuum()->upperEnergy();
          bool in_insignificant_roi = false;
          for( const size_t index : insignificant_indices )
          {
            if( index >= ranges.size() )
              continue;
            const RelActCalcAuto::RoiRange &range = ranges[index];
            if( (std::fabs(lower - range.lower_energy) < 1.0)
                && (std::fabs(upper - range.upper_energy) < 1.0) )
            {
              in_insignificant_roi = true;
              break;
            }
          }
          if( in_insignificant_roi )
            continue;
          bool mean_in_significant_roi = false;
          for( size_t index = 0; index < ranges.size(); ++index )
          {
            if( std::find( std::begin(insignificant_indices),
                    std::end(insignificant_indices), index )
                != std::end(insignificant_indices) )
              continue;
            const RelActCalcAuto::RoiRange &range = ranges[index];
            if( (peak.mean() >= range.lower_energy)
                && (peak.mean() <= range.upper_energy) )
            {
              mean_in_significant_roi = true;
              break;
            }
          }
          if( mean_in_significant_roi && (peak.mean() >= lower) && (peak.mean() <= upper) )
            fitted.push_back( peak );
        }
        std::vector<PeakDef> combined = combine_overlapping_peaks_in_rois( fitted );
        return compute_observable_peaks(
            combined, candidate.m_foreground, det_type, config );
      };
    const auto observable_requested_anchors_preserved
      = [&]( const RelActCalcAuto::RelActAutoSolution &incumbent,
             const RelActCalcAuto::RelActAutoSolution &challenger ) {
        const std::vector<PeakDef> incumbent_observable
          = observable_requested_peaks( incumbent );
        const std::vector<PeakDef> challenger_observable
          = observable_requested_peaks( challenger );
        for( const PeakDef &anchor : incumbent_observable )
        {
          const bool found = std::any_of( std::begin(challenger_observable),
              std::end(challenger_observable),
              [&anchor, &same_gamma_identity, &solution_fwhm_at_energy](
                  const PeakDef &candidate_peak ) {
                if( same_gamma_identity(anchor, candidate_peak) )
                  return true;
                if( anchor.hasSourceGammaAssigned()
                    || candidate_peak.hasSourceGammaAssigned() )
                  return false;
                const double anchor_fwhm = solution_fwhm_at_energy( anchor.mean() );
                const double candidate_fwhm
                  = solution_fwhm_at_energy( candidate_peak.mean() );
                return std::isfinite(anchor_fwhm) && (anchor_fwhm > 0.0)
                    && std::isfinite(candidate_fwhm) && (candidate_fwhm > 0.0)
                    && (std::fabs(anchor.mean() - candidate_peak.mean())
                      < 0.5*(anchor_fwhm + candidate_fwhm));
              } );
          if( !found )
          {
            if( should_debug_print() )
              std::cerr << "Observable-anchor guard: lost requested peak at "
                        << anchor.mean() << " keV" << std::endl;
            return false;
          }
        }
        return true;
      };
    const auto requested_anchors_catastrophically_removed
      = [&]( const RelActCalcAuto::RelActAutoSolution &incumbent,
             const RelActCalcAuto::RelActAutoSolution &challenger ) -> bool {
        size_t num_incumbent = 0;
        size_t num_matched = 0;
        for( const PeakDef &anchor : incumbent.m_peaks_without_back_sub )
        {
          const double anchor_z = peak_fit_significance( anchor );
          if( !is_requested_peak(anchor) || (anchor_z < config.roi_significance_z) )
            continue;
          ++num_incumbent;
          const bool matched = std::any_of(
            std::begin(challenger.m_peaks_without_back_sub),
            std::end(challenger.m_peaks_without_back_sub),
            [&]( const PeakDef &candidate ) {
              return same_gamma_identity(anchor, candidate)
                  && (peak_fit_significance(candidate) >= config.roi_significance_z);
            } );
          num_matched += matched ? 1u : 0u;
          if( !matched && (anchor_z >= std::max(8.0, 2.0*config.roi_significance_z)) )
          {
            if( should_debug_print() )
              std::cerr << "Rescue anchor guard: lost strong requested peak at " << anchor.mean()
                        << " keV (z=" << anchor_z << ")" << std::endl;
            return true;
          }
        }
        if( num_incumbent == 0 )
          return false;
        const bool catastrophic = (num_matched == 0) || (5*num_matched < 4*num_incumbent);
        if( catastrophic && should_debug_print() )
          std::cerr << "Rescue anchor guard: matched " << num_matched << " of "
                    << num_incumbent << " requested anchors" << std::endl;
        return catastrophic;
      };

    // R6 transaction: augment a successful source-only incumbent with at most two foreground-NORM
    // nuisance nuclides.  A supplied background already models/subtracts these lines, so raw-
    // foreground confirmation is not a residual test and must not add another curve in that mode.
    // Any failure, source-anchor loss, or insufficient nested-likelihood gain keeps the incumbent.
    if( !fit_norm_peaks && !orig_background && !interferer_candidates.empty()
        && !solution.m_peaks_without_back_sub.empty() )
    {
      try
      {
        struct RankedInterferer
        {
          const SandiaDecay::Nuclide *nuclide = nullptr;
          double max_detection_z = 0.0;
          std::vector<double> line_energies;
        };

        std::vector<RankedInterferer> ranked;
        for( const detail::InterfererCandidate &candidate : interferer_candidates )
        {
          if( !candidate.nuclide )
            continue;  // the unattributable floating path remains disabled

          auto pos = std::find_if( std::begin(ranked), std::end(ranked),
            [&candidate]( const RankedInterferer &entry ) {
              return entry.nuclide == candidate.nuclide;
            } );
          if( pos == std::end(ranked) )
          {
            RankedInterferer entry;
            entry.nuclide = candidate.nuclide;
            entry.max_detection_z = candidate.detection_z;
            entry.line_energies.push_back( candidate.energy );
            ranked.push_back( std::move(entry) );
          }else
          {
            pos->max_detection_z = std::max( pos->max_detection_z, candidate.detection_z );
            pos->line_energies.push_back( candidate.energy );
          }
        }

        std::sort( std::begin(ranked), std::end(ranked),
          []( const RankedInterferer &lhs, const RankedInterferer &rhs ) {
            if( lhs.max_detection_z != rhs.max_detection_z )
              return lhs.max_detection_z > rhs.max_detection_z;
            if( lhs.nuclide->atomicNumber != rhs.nuclide->atomicNumber )
              return lhs.nuclide->atomicNumber < rhs.nuclide->atomicNumber;
            if( lhs.nuclide->massNumber != rhs.nuclide->massNumber )
              return lhs.nuclide->massNumber < rhs.nuclide->massNumber;
            return lhs.nuclide->isomerNumber < rhs.nuclide->isomerNumber;
          } );
        RelActCalcAuto::Options augmented_options = solution.m_options;
        const std::vector<RelActCalcAuto::NucInputInfo> norm_all
          = get_norm_sources( sources, config.norm_css_color );
        std::vector<RankedInterferer> selected;
        for( RankedInterferer entry : ranked )
        {
          if( selected.size() >= sm_max_auto_interferer_nuclides )
            break;
          double nuisance_age = PeakDef::defaultDecayTime( entry.nuclide );
          for( const RelActCalcAuto::NucInputInfo &norm_source : norm_all )
          {
            if( RelActCalcAuto::nuclide(norm_source.source) == entry.nuclide )
            {
              nuisance_age = norm_source.age;
              break;
            }
          }

          // A nuisance curve emits every parent line in every fitted ROI, not only the line that
          // confirmed discovery.  If any such modeled line duplicates an explicit floating peak,
          // exclude the whole parent rather than allowing a second candidate line to reintroduce it.
          const std::vector<SandiaDecay::EnergyRatePair> parent_photons
            = get_source_photons( entry.nuclide, 1.0, nuisance_age );
          const bool parent_duplicates_floating_peak = std::any_of(
            std::begin(parent_photons), std::end(parent_photons),
            [&]( const SandiaDecay::EnergyRatePair &photon ) {
              const bool in_fitted_roi = std::any_of(
                std::begin(augmented_options.rois), std::end(augmented_options.rois),
                [&]( const RelActCalcAuto::RoiRange &roi ) {
                  return (photon.energy >= roi.lower_energy)
                      && (photon.energy <= roi.upper_energy);
                } );
              if( !in_fitted_roi )
                return false;
              const double line_fwhm = solution_fwhm_at_energy( photon.energy );
              return std::any_of(
                std::begin(augmented_options.floating_peaks),
                std::end(augmented_options.floating_peaks),
                [&]( const RelActCalcAuto::FloatingPeak &peak ) {
                  return std::fabs( peak.energy - photon.energy ) < line_fwhm;
                } );
            } );
          if( parent_duplicates_floating_peak )
            continue;

          std::vector<double> usable_lines;
          for( const double energy : entry.line_energies )
          {
            const double line_fwhm = solution_fwhm_at_energy( energy );
            const bool duplicates_floating_peak = std::any_of(
              std::begin(augmented_options.floating_peaks),
              std::end(augmented_options.floating_peaks),
              [&]( const RelActCalcAuto::FloatingPeak &peak ) {
                return std::fabs( peak.energy - energy ) < line_fwhm;
              } );
            if( duplicates_floating_peak )
              continue;

            bool covered = false;
            for( const RelActCalcAuto::RoiRange &roi : augmented_options.rois )
            {
              if( (energy >= roi.lower_energy) && (energy <= roi.upper_energy) )
              {
                covered = true;
                break;
              }
            }

            if( !covered )
            {
              const double half = std::max( 1.0,
                  config.auto_roi_core_num_fwhm * line_fwhm );
              const double lo = std::max( min_valid_energy, energy - half );
              const double hi = std::min( max_valid_energy, energy + half );
              bool conflicts_with_existing = false;
              for( const ExistingRoiInfo &existing : existing_roi_ranges )
              {
                if( (lo < existing.upper_energy) && (hi > existing.lower_energy) )
                {
                  conflicts_with_existing = true;
                  break;
                }
              }
              if( conflicts_with_existing )
                continue;
            }

            usable_lines.push_back( energy );
          }

          if( !usable_lines.empty() )
          {
            entry.line_energies = std::move( usable_lines );
            selected.push_back( std::move(entry) );
          }
        }

        if( !selected.empty() )
        {
          std::vector<RelActCalcAuto::NucInputInfo> interferer_nucs;
          std::vector<double> selected_lines;
          for( const RankedInterferer &entry : selected )
          {
            bool from_norm = false;
            for( const RelActCalcAuto::NucInputInfo &norm_source : norm_all )
            {
              if( RelActCalcAuto::nuclide( norm_source.source ) == entry.nuclide )
              {
                interferer_nucs.push_back( norm_source );
                from_norm = true;
                break;
              }
            }
            if( !from_norm )
            {
              RelActCalcAuto::NucInputInfo info;
              info.source = entry.nuclide;
              info.age = 0.0;
              info.fit_age = false;
              info.peak_color_css = config.norm_css_color;
              interferer_nucs.push_back( info );
            }
            selected_lines.insert( std::end(selected_lines),
                std::begin(entry.line_energies), std::end(entry.line_energies) );
          }

          std::shared_ptr<RelActCalc::PhysicalModelShieldInput> self_atten
            = std::make_shared<RelActCalc::PhysicalModelShieldInput>();
          self_atten->atomic_number = 10.4;
          self_atten->areal_density
            = (1.6 * PhysicalUnits::g / PhysicalUnits::cm3) * (100 * PhysicalUnits::cm);
          self_atten->fit_atomic_number = false;
          self_atten->fit_areal_density = false;

          RelActCalcAuto::RelEffCurveInput interferer_curve;
          interferer_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
          interferer_curve.rel_eff_eqn_order = 0;
          interferer_curve.nucs_of_el_same_age = false;
          interferer_curve.phys_model_corr.corr_fcn = RelActCalc::PhysModelCorrFcn::None;
          interferer_curve.phys_model_self_atten = self_atten;
          interferer_curve.nuclides = interferer_nucs;
          interferer_curve.name = "Interfering-line curve";
          augmented_options.rel_eff_curves.push_back( interferer_curve );

          for( const double energy : selected_lines )
          {
            bool covered = false;
            for( const RelActCalcAuto::RoiRange &roi : augmented_options.rois )
            {
              if( (energy >= roi.lower_energy) && (energy <= roi.upper_energy) )
              {
                covered = true;
                break;
              }
            }
            if( covered )
              continue;

            const double half = std::max( 1.0,
                config.auto_roi_core_num_fwhm * solution_fwhm_at_energy( energy ) );
            double lo = std::max( min_valid_energy, energy - half );
            double hi = std::min( max_valid_energy, energy + half );
            std::vector<RelActCalcAuto::RoiRange> kept;
            PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;
            bool merged_any = false;
            for( const RelActCalcAuto::RoiRange &roi : augmented_options.rois )
            {
              if( (roi.upper_energy >= lo) && (roi.lower_energy <= hi) )
              {
                if( !merged_any )
                  continuum_type = roi.continuum_type;
                lo = std::min( lo, roi.lower_energy );
                hi = std::max( hi, roi.upper_energy );
                merged_any = true;
              }else
              {
                kept.push_back( roi );
              }
            }

            RelActCalcAuto::RoiRange coverage;
            coverage.lower_energy = lo;
            coverage.upper_energy = hi;
            coverage.continuum_type = continuum_type;
            coverage.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
            kept.push_back( coverage );
            augmented_options.rois = std::move( kept );
          }

          std::sort( std::begin(augmented_options.rois), std::end(augmented_options.rois),
            []( const RelActCalcAuto::RoiRange &lhs, const RelActCalcAuto::RoiRange &rhs ) {
              return lhs.lower_energy < rhs.lower_energy;
            } );
          resolve_overlapping_rois( augmented_options.rois, augmented_options.floating_peaks );
          ensure_min_channel_gap( augmented_options.rois, orig_foreground->energy_calibration() );
          AutomaticRoiDecisionDiagnostic r6_bypass;
          r6_bypass.decision = AutomaticRoiDecision::R6LegacyBypass;
          r6_bypass.stage = "R6 transactional coverage union";
          r6_bypass.reason = "R6 legacy transactional geometry intentionally preserved";
          result.automatic_roi_diagnostics.push_back( r6_bypass );
          remove_floating_peaks_without_roi( augmented_options );

          RelActCalcAuto::Options common_domain_source_options = augmented_options;
          common_domain_source_options.rel_eff_curves.pop_back();
          const RelActCalcAuto::RelActAutoSolution common_domain_source_solution
            = RelActCalcAuto::solve( common_domain_source_options, orig_foreground,
                orig_background, drf, auto_search_peaks, det_type );
          RelActCalcAuto::RelActAutoSolution augmented_solution = RelActCalcAuto::solve(
              augmented_options, orig_foreground, orig_background, drf, auto_search_peaks, det_type );

          size_t incumbent_requested_count = 0;
          for( const PeakDef &peak : solution.m_peaks_without_back_sub )
            incumbent_requested_count += is_requested_peak( peak ) ? 1u : 0u;

          size_t augmented_requested_count = 0;
          for( const PeakDef &peak : augmented_solution.m_peaks_without_back_sub )
            augmented_requested_count += is_requested_peak( peak ) ? 1u : 0u;

          const auto is_selected_interferer = [&selected]( const PeakDef &peak ) -> bool {
            const SandiaDecay::Nuclide * const parent = peak.parentNuclide();
            if( !parent )
              return false;
            return std::any_of( std::begin(selected), std::end(selected),
              [parent]( const RankedInterferer &entry ) {
                return entry.nuclide == parent;
              } );
          };
          const auto is_confirming_interferer_peak = [&]( const PeakDef &peak ) -> bool {
            if( !is_selected_interferer(peak) || !peak.hasSourceGammaAssigned() )
              return false;
            return std::any_of( std::begin(selected_lines), std::end(selected_lines),
              [&peak]( const double energy ) {
                return std::fabs(peak.gammaParticleEnergy() - energy) < 0.1;
              } );
          };

          // Evaluate the nested source-only and source+nuisance solves on the identical affected
          // ROI channels.  The source-only control uses the augmented ROI geometry, so its tied
          // rel-eff/activity model cannot freely assign the interferer counts to one weak source
          // line, while the comparison remains independent of any added non-affected channels.
          size_t num_affected_rois = 0;
          bool affected_domains_valid = true;
          double affected_deviance_source = 0.0;
          double affected_deviance_augmented = 0.0;
          std::set<const SandiaDecay::Nuclide *> contributing_interferers;
          if( RelActCalcAuto::RelActAutoSolution::is_usable_status(
                common_domain_source_solution.m_status)
              && RelActCalcAuto::RelActAutoSolution::is_usable_status(
                augmented_solution.m_status) )
          {
            for( const RelActCalcAuto::RoiRange &roi : augmented_options.rois )
            {
              const bool affected = std::any_of(
                std::begin(selected_lines), std::end(selected_lines),
                [&]( const double energy ) {
                  if( (energy < roi.lower_energy) || (energy > roi.upper_energy) )
                    return false;
                  return std::any_of(
                    std::begin(augmented_solution.m_peaks_without_back_sub),
                    std::end(augmented_solution.m_peaks_without_back_sub),
                    [&]( const PeakDef &peak ) {
                      return is_confirming_interferer_peak(peak)
                          && (std::fabs(peak.gammaParticleEnergy() - energy) < 0.1);
                    } );
                } );
              if( !affected )
                continue;

              std::vector<std::shared_ptr<const PeakDef>> source_model_peaks;
              std::vector<std::shared_ptr<const PeakDef>> augmented_model_peaks;
              const auto append_roi_peaks = [&roi]( const std::vector<PeakDef> &peaks,
                  std::vector<std::shared_ptr<const PeakDef>> &output ) {
                for( const PeakDef &peak : peaks )
                {
                  const std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
                  if( !continuum || (continuum->upperEnergy() <= roi.lower_energy)
                      || (continuum->lowerEnergy() >= roi.upper_energy) )
                    continue;
                  output.push_back( std::make_shared<PeakDef>(peak) );
                }
              };
              append_roi_peaks( common_domain_source_solution.m_peaks_without_back_sub,
                                source_model_peaks );
              append_roi_peaks( augmented_solution.m_peaks_without_back_sub,
                                augmented_model_peaks );
              if( source_model_peaks.empty() || augmented_model_peaks.empty() )
              {
                affected_domains_valid = false;
                break;
              }

              const size_t lower_channel
                = orig_foreground->find_gamma_channel(roi.lower_energy);
              const size_t upper_channel = std::min(
                  orig_foreground->find_gamma_channel(roi.upper_energy),
                  orig_foreground->num_gamma_channels() - 1 );
              if( upper_channel <= lower_channel )
              {
                affected_domains_valid = false;
                break;
              }
              const std::shared_ptr<const SpecUtils::Measurement> source_model_data
                = common_domain_source_solution.m_foreground
                  ? common_domain_source_solution.m_foreground : orig_foreground;
              const std::shared_ptr<const SpecUtils::Measurement> augmented_model_data
                = augmented_solution.m_foreground
                  ? augmented_solution.m_foreground : orig_foreground;
              const FixedRoiModelScore source_score = fixed_roi_model_score(
                  source_model_peaks, source_model_data, lower_channel, upper_channel,
                  source_model_data->gamma_channel_lower(lower_channel),
                  source_model_data->gamma_channel_upper(upper_channel) );
              const FixedRoiModelScore augmented_score = fixed_roi_model_score(
                  augmented_model_peaks, augmented_model_data, lower_channel, upper_channel,
                  augmented_model_data->gamma_channel_lower(lower_channel),
                  augmented_model_data->gamma_channel_upper(upper_channel) );
              if( !source_score.valid || !augmented_score.valid
                  || (source_score.num_channels != augmented_score.num_channels) )
              {
                affected_domains_valid = false;
                break;
              }
              ++num_affected_rois;
              affected_deviance_source += source_score.poisson_deviance;
              affected_deviance_augmented += augmented_score.poisson_deviance;
              for( const std::shared_ptr<const PeakDef> &peak : augmented_model_peaks )
              {
                if( peak && is_confirming_interferer_peak(*peak)
                    && (peak->gammaParticleEnergy() >= roi.lower_energy)
                    && (peak->gammaParticleEnergy() <= roi.upper_energy)
                    && (peak_fit_significance(*peak) >= sm_interferer_min_detect_z) )
                  contributing_interferers.insert( peak->parentNuclide() );
              }
            }
          }

          const bool have_affected_roi = affected_domains_valid && (num_affected_rois > 0);
          const bool every_selected_parent_contributed
            = (contributing_interferers.size() == selected.size());
          const double affected_delta_chi2
            = affected_deviance_source - affected_deviance_augmented;
          double affected_nested_z = -40.0;
          if( have_affected_roi && (affected_delta_chi2 > 0.0)
              && !contributing_interferers.empty() )
          {
            const boost::math::chi_squared_distribution<double> chi2_dist(
                static_cast<double>(contributing_interferers.size()) );
            const double p_value = boost::math::cdf(
                boost::math::complement( chi2_dist, affected_delta_chi2 ) );
            if( p_value < 1.0e-300 )
              affected_nested_z = 40.0;
            else if( p_value < (1.0 - 1.0e-12) )
            {
              const boost::math::normal_distribution<double> normal_dist;
              affected_nested_z = -boost::math::quantile( normal_dist, p_value );
            }
          }
          const bool affected_rois_significant
            = have_affected_roi && (affected_nested_z >= config.roi_significance_z);

          const bool anchors_preserved = requested_anchors_preserved(
              solution, augmented_solution, selected_lines );
          const bool accept = RelActCalcAuto::RelActAutoSolution::is_usable_status(
                augmented_solution.m_status)
              && (incumbent_requested_count > 0)
              && (augmented_requested_count > 0)
              && anchors_preserved
              && have_affected_roi
              && every_selected_parent_contributed
              && affected_rois_significant;

          if( accept )
          {
            solution = std::move( augmented_solution );
            auto_interferer_lines = selected_lines;
            std::string names;
            for( const RankedInterferer &entry : selected )
            {
              auto_interferer_nucs.insert( entry.nuclide );
              names += (names.empty() ? "" : ", ") + entry.nuclide->symbol;
            }
            result.warnings.push_back( "Auto co-fit interfering line(s) from " + names
              + " after preserving the successful source-only fit; these nuisance peaks are hidden"
                " from the returned peak vectors." );
          }else
          {
            std::ostringstream warning;
            warning << "Rejected the automatic interfering-line augmentation because it failed,"
              " removed a requested-source anchor, or lacked significant affected-ROI likelihood"
              " gain; retained the successful source-only fit";
            warning << " (solve=" << RelActCalcAuto::RelActAutoSolution::is_usable_status(
                                          augmented_solution.m_status)
                    << ", requested=" << incumbent_requested_count << "->"
                    << augmented_requested_count << ", anchors=" << anchors_preserved
                    << ", affected_roi=" << have_affected_roi
                    << ", every_parent_contributed=" << every_selected_parent_contributed
                    << ", affected_significant=" << affected_rois_significant;
            if( have_affected_roi && std::isfinite(affected_nested_z) )
              warning << ", affected nested z=" << affected_nested_z << ", required "
                      << config.roi_significance_z;
            warning << ").";
            result.warnings.push_back( warning.str() );
          }
        }//if( !selected.empty() )
      }catch( const std::exception &e )
      {
        result.warnings.push_back( "Automatic interfering-line augmentation failed; retained the"
          " successful source-only fit: " + std::string(e.what()) );
      }
    }//transactional R6 augmentation

    const size_t print_curve_idx = (sources_rel_eff_index >= 0) ? static_cast<size_t>( sources_rel_eff_index ) : 0u;
    //std::cout << "Initial RelActAuto solution (" << options.rel_eff_curves[print_curve_idx].name << "):" << std::endl;
    //solution.print_summary( std::cout );
    //std::cout << "Chi2/DOF = " << solution.m_chi2 << "/" << solution.m_dof << " = " << (solution.m_chi2 / solution.m_dof) << std::endl;

    // We may adjust energy calibration - sow we'll use `background` and `foreground` for this.
    shared_ptr<const SpecUtils::Measurement> foreground = orig_foreground;
    shared_ptr<const SpecUtils::Measurement> background = orig_background;

    // Iteratively refine ROIs using RelActAuto solutions
    // The idea is that each iteration provides a better relative efficiency estimate,
    // which allows us to better identify significant gamma lines and create better ROIs.
    std::vector<RelActCalcAuto::RoiRange> rescued_roi_ranges;
    {
      const size_t max_iterations = 3;
      size_t num_extra_allowed = 0; //If we switch to our "desperation" model type retry - we will increment this to 1.
      bool rescue_attempted = false;
      std::vector<RelActCalcAuto::RoiRange> rescue_rejected_ranges;
      for( size_t iter = 0; iter < (max_iterations + num_extra_allowed); ++iter )
      {
        bool rescue_solve_this_iteration = false;
        std::vector<std::shared_ptr<const PeakDef>> current_unfit_auto_peaks;
        std::vector<RelActCalcAuto::RoiRange> current_protected_mixed_rois;
        bool calibration_guards_valid = true;
        const auto refresh_calibrated_guards = [&]() {
          current_unfit_auto_peaks = unfit_auto_peaks;
          current_protected_mixed_rois.clear();
          const std::shared_ptr<const SpecUtils::EnergyCalibration> original_cal
            = orig_foreground->energy_calibration();
          const std::shared_ptr<const SpecUtils::EnergyCalibration> working_cal
            = foreground->energy_calibration();
          if( original_cal && working_cal && (original_cal != working_cal)
              && !unfit_auto_peaks.empty() )
          {
            try
            {
              const std::deque<std::shared_ptr<const PeakDef>> original_peaks(
                  std::begin(unfit_auto_peaks), std::end(unfit_auto_peaks) );
              const std::deque<std::shared_ptr<const PeakDef>> translated
                = EnergyCal::translatePeaksForCalibrationChange(
                    original_peaks, original_cal, working_cal );
              current_unfit_auto_peaks.assign( std::begin(translated), std::end(translated) );
            }catch( const std::exception &error )
            {
              calibration_guards_valid = false;
              current_unfit_auto_peaks.clear();
              result.warnings.push_back( "Could not transform automatic ROI guards into the"
                " current calibration; stopped ROI refinement and retained the successful"
                " incumbent (" + std::string(error.what()) + ")." );
            }
          }
          for( RelActCalcAuto::RoiRange protected_roi : protected_mixed_rois )
          {
            if( original_cal && working_cal && (original_cal != working_cal) )
            {
              const double lower_channel = original_cal->channel_for_energy(
                  protected_roi.lower_energy );
              const double upper_channel = original_cal->channel_for_energy(
                  protected_roi.upper_energy );
              protected_roi.lower_energy = working_cal->energy_for_channel( lower_channel );
              protected_roi.upper_energy = working_cal->energy_for_channel( upper_channel );
            }
            current_protected_mixed_rois.push_back( protected_roi );
          }
        };
        refresh_calibrated_guards();
        bool rescue_transaction_failed = false;
        std::vector<RelActCalcAuto::RoiRange> proposed_rescued_rois;
        if( apply_energy_cal_between && config.fit_energy_cal )
        {
          const shared_ptr<SpecUtils::EnergyCalibration> fitted_cal = solution.get_adjusted_energy_cal();
          const shared_ptr<const SpecUtils::EnergyCalibration> orig_fg_cal = foreground->energy_calibration();

#if( PERFORM_DEVELOPER_CHECKS )
          if( should_debug_print() )
          {
            std::cout << "Energy cal iteration " << iter << ":" << std::endl;
            std::cout << "  Linear adjustments:";
            for( size_t ei = 0; ei < solution.m_energy_cal_adjustments.size(); ++ei )
            {
              const double adj_keV = (solution.m_energy_cal_adjustments[ei]
                / RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                * RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
              std::cout << " [" << ei << "]=" << adj_keV << " keV"
                   << (solution.m_fit_energy_cal[ei] ? " (fit)" : " (fixed)");
            }
            std::cout << std::endl;

            if( !solution.m_deviation_pair_offsets.empty() )
            {
              std::cout << "  Deviation pairs (" << solution.m_deviation_pair_offsets.size()
                   << ", " << solution.m_num_deviations_fit << " fit):" << std::endl;
              for( const auto &dp : solution.m_deviation_pair_offsets )
                std::cout << "    energy=" << dp.first << " keV, offset=" << dp.second << " keV" << std::endl;
            }

            // Show a few sample energy differences between original and fitted cal
            std::cout << "  Sample energy shifts (orig -> fitted):" << std::endl;
            for( const double e : { 50.0, 81.0, 160.0, 233.0, 303.0, 384.0, 500.0, 1000.0 } )
            {
              const double orig_ch = orig_fg_cal->channel_for_energy( e );
              const double fitted_e = fitted_cal->energy_for_channel( orig_ch );
              std::cout << "    " << e << " keV: ch=" << std::setprecision(2) << orig_ch
                   << ", fitted_e=" << std::setprecision(4) << fitted_e
                   << ", shift=" << std::setprecision(4) << (fitted_e - e) << " keV" << std::endl;
            }
          }//if( should_debug_print() )
#endif

          shared_ptr<SpecUtils::Measurement> new_foreground = make_shared<SpecUtils::Measurement>( *foreground );
          new_foreground->set_energy_calibration( fitted_cal );
          foreground = new_foreground;

          // R1 step 2: keep the shared global SNIP continuum in the SAME (fitted) energy frame as the
          // working foreground, so its energy->channel lookups match the data the gates now see.  The
          // fit only re-labels the calibration (never re-bins), so the SNIP's per-channel counts are
          // invariant - we just re-stamp the calibration instead of recomputing SNIP (cheap + exact).
          if( global_cont.valid() )
          {
            shared_ptr<SpecUtils::Measurement> recal_snip = make_shared<SpecUtils::Measurement>( *global_cont.snip );
            recal_snip->set_energy_calibration( fitted_cal );
            global_cont.snip = recal_snip;
            global_cont.foreground = foreground;  // same original counts, now on the fitted cal
          }

          if( background )
          {
            // Propagate the foreground cal change (orig_fg_cal -> fitted_cal) to
            // the background.  propogate_energy_cal_change requires orig and new
            // cals to be Polynomial or FRF; if LowerChannelEdge just use fitted_cal directly.
            shared_ptr<const SpecUtils::EnergyCalibration> new_bkg_cal;
            if( orig_fg_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge )
              new_bkg_cal = fitted_cal;
            else
              new_bkg_cal = EnergyCal::propogate_energy_cal_change( orig_fg_cal, fitted_cal, background->energy_calibration() );

            shared_ptr<SpecUtils::Measurement> new_background = make_shared<SpecUtils::Measurement>( *background );
            new_background->set_energy_calibration( new_bkg_cal );
            background = new_background;
          }
        }//if( apply_energy_cal_between && config.fit_energy_cal )

        // Energy-cal refinement above re-labels foreground/SNIP channels.  Re-materialize every
        // guard only after that advance so no policy input is left in the prior frame.
        refresh_calibrated_guards();
        if( !calibration_guards_valid )
          break;


        vector<function<double(double)>> auto_rel_effs;
        vector<vector<tuple<RelActCalcAuto::SrcVariant,double,double>>> source_age_and_acts;

        // Energy span of the ROIs actually used in the most recent fit.  When clustering candidate
        // gammas we CLAMP the rel-eff evaluation to this span: beyond the fitted ROIs the curve is
        // being extrapolated, and polynomial/physical extrapolation is unreliable (it can droop toward
        // zero or blow up), which would wrongly starve real lines or invent ROIs at energies the fit
        // had no data for.  The LOWER side holds the boundary efficiency flat; the UPPER side is
        // shaped by the DRF intrinsic-efficiency falloff so far high-energy gammas get
        // a realistic declining efficiency instead of an over-optimistic flat hold.
        double rel_eff_clamp_lo = std::numeric_limits<double>::max();
        double rel_eff_clamp_hi = std::numeric_limits<double>::lowest();
        for( const RelActCalcAuto::RoiRange &roi : solution.m_final_roi_ranges )
        {
          rel_eff_clamp_lo = std::min( rel_eff_clamp_lo, roi.lower_energy );
          rel_eff_clamp_hi = std::max( rel_eff_clamp_hi, roi.upper_energy );
        }
        const bool have_rel_eff_clamp = (rel_eff_clamp_lo < rel_eff_clamp_hi);

        // DRF used to shape rel-eff extrapolation ABOVE the fitted span.
        const std::shared_ptr<const DetectorPeakResponse> rel_eff_extrap_drf
          = generic_drf_for_rel_eff_extrap( drf, det_type );

        for( size_t rel_eff_index = 0; rel_eff_index < solution.m_rel_activities.size(); ++rel_eff_index )
        {
          source_age_and_acts.push_back( {} );

          // rel_eff lambda from the current RelActAuto solution.  Below the fitted-ROI span hold the
          // boundary efficiency flat; above it, shape the boundary value by the DRF intrinsic-
          // efficiency falloff so far high-energy gammas get a realistic declining efficiency.
          auto_rel_effs.push_back(
            [&solution, rel_eff_index, rel_eff_clamp_lo, rel_eff_clamp_hi, have_rel_eff_clamp, rel_eff_extrap_drf]( double energy ) -> double {
              if( !have_rel_eff_clamp )
                return solution.relative_efficiency( energy, rel_eff_index );
              if( energy <= rel_eff_clamp_hi )
                return solution.relative_efficiency( std::max( energy, rel_eff_clamp_lo ), rel_eff_index );
              const double re_hi = solution.relative_efficiency( rel_eff_clamp_hi, rel_eff_index );
              return shape_rel_eff_above_boundary( re_hi, energy, rel_eff_clamp_hi, rel_eff_extrap_drf );
          } );

          // Collect sources and activities from the current solution
          // m_rel_activities is a 2D vector: [rel_eff_curve_index][source_index]

          if( should_debug_print() )
            std::cout << "Collecting " << solution.m_rel_activities[rel_eff_index].size() << " sources from RelActAuto solution:" << std::endl;

          for( const RelActCalcAuto::NuclideRelAct &nuc_act : solution.m_rel_activities[rel_eff_index] )
          {
            if( RelActCalcAuto::is_null( nuc_act.source ) )
              continue;
            const SandiaDecay::Nuclide * const activity_parent
              = RelActCalcAuto::nuclide( nuc_act.source );
            if( activity_parent && auto_interferer_nucs.count(activity_parent) )
              continue;  // R6 nuisances may refine only inside their confirmed incumbent ROIs.

            const double live_time_seconds = foreground->live_time();
            // RelActAuto's rel_activity is per second, need to multiply by live time for clustering
            const double activity_for_clustering = nuc_act.rel_activity * live_time_seconds;

            if( should_debug_print() )
              std::cout << "  " << nuc_act.name() << ": rel_activity=" << nuc_act.rel_activity
                   << ", live_time=" << live_time_seconds
                   << "s, activity_for_clustering=" << activity_for_clustering << std::endl;

            source_age_and_acts.back().emplace_back( nuc_act.source, nuc_act.age, activity_for_clustering );
          }//for( loop over solution.m_rel_activities[rel_eff_index] )
        }//for( size_t rel_eff_index = 0; rel_eff_index < ; ++rel_eff_index )

        // Use the FWHM coefficients passed to this function (computed from auto-search peaks
        // or DRF in fit_peaks_for_nuclides), rather than relying on solution.m_drf which may
        // not have valid FWHM info or may have incorrect values

        // Get auto clustering settings from config
        GammaClusteringSettings auto_settings = config.get_auto_clustering_settings();
        auto_settings.use_automatic_roi_policy = use_automatic_roi_policy;
        auto_settings.global_continuum = global_cont.valid() ? &global_cont : nullptr;  // R1 step 2

        // Cluster gammas using current solution's relative efficiency
        std::vector<MarginalRejectedCluster> marginal_rejects;
        std::vector<PredictedGamma> fitted_cluster_predictions;
        const std::string refinement_stage = "refinement " + std::to_string(iter);
        std::vector<std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo>> refined_rois_and_gammas
          = cluster_gammas_to_rois(
              auto_rel_effs, source_age_and_acts, foreground,
              fwhm_form, fwhm_coefficients,
              fwhm_lower_energy, fwhm_upper_energy,
              min_valid_energy, max_valid_energy,
              auto_settings, current_unfit_auto_peaks,
              (!rescue_attempted && (iter == 0)) ? &marginal_rejects : nullptr,
              nullptr, nullptr,
              (!rescue_attempted && (iter == 0)) ? &fitted_cluster_predictions : nullptr,
              refinement_stage, nullptr,
              &result.automatic_roi_diagnostics,
              !recovered_source_anchors.empty() );
        std::vector<std::pair<double,double>> current_modeled_peak_candidates;
        for( const std::pair<RelActCalcAuto::RoiRange,ClusteredGammaInfo> &entry
             : refined_rois_and_gammas )
        {
          for( size_t peak_index = 0;
               peak_index < entry.second.gamma_energies.size(); ++peak_index )
          {
            const double area = (peak_index < entry.second.gamma_amplitudes.size())
                ? entry.second.gamma_amplitudes[peak_index] : 0.0;
            current_modeled_peak_candidates.emplace_back(
                entry.second.gamma_energies[peak_index], area );
          }
        }
        for( const std::pair<double,double> &candidate : initial_modeled_peak_candidates )
        {
          const bool duplicate = std::any_of(
              std::begin(current_modeled_peak_candidates),
              std::end(current_modeled_peak_candidates),
              [&candidate]( const std::pair<double,double> &current_candidate ) {
                return std::fabs(current_candidate.first - candidate.first) < 1.0e-6;
              } );
          if( !duplicate )
            current_modeled_peak_candidates.push_back( candidate );
        }

        // Carry the exact measured-data partition that was selected in this calibration frame.
        // Capture it before the ordinary interference-shrinking pass below: the local transaction
        // must solve the identical channel bins that won the measured-data AICc comparison, not a
        // later edge refinement that was never scored by that comparison.  This is deliberately
        // independent of diagnostics: rollback must not reconstruct geometry from a reason string
        // or one representative boundary.
        std::vector<SelectedRoiComponentPartition> selected_component_partitions;
        std::map<uint64_t,size_t> selected_partition_indices;
        for( const std::pair<RelActCalcAuto::RoiRange,ClusteredGammaInfo> &entry
             : refined_rois_and_gammas )
        {
          const ClusteredGammaInfo &cluster = entry.second;
          if( cluster.selected_partition_id == 0 )
            continue;
          const auto inserted = selected_partition_indices.emplace(
              cluster.selected_partition_id, selected_component_partitions.size() );
          if( inserted.second )
          {
            SelectedRoiComponentPartition record;
            record.id = cluster.selected_partition_id;
            record.parent_lower = cluster.selected_parent_lower;
            record.parent_upper = cluster.selected_parent_upper;
            record.calibration = foreground->energy_calibration();
            selected_component_partitions.push_back( std::move(record) );
          }
          selected_component_partitions[inserted.first->second].children.push_back( entry.first );
        }
        selected_component_partitions.erase( std::remove_if(
            std::begin(selected_component_partitions), std::end(selected_component_partitions),
            []( const SelectedRoiComponentPartition &record ) {
              return (record.children.size() < 2)
                  || !std::isfinite(record.parent_lower)
                  || !std::isfinite(record.parent_upper)
                  || !(record.parent_upper > record.parent_lower);
            } ), std::end(selected_component_partitions) );

        // Early predicted-gamma clustering cannot see every range that the solver will ultimately
        // construct.  In particular, a final fitted continuum can span widely separated, resolved
        // peaks even when no early transitive cluster exposed the split.  Propose a partition from
        // the accepted solver geometry itself, but do not edit that geometry in place: the existing
        // component-local transaction below re-solves the exact child ranges and retains them only
        // if its requested-source anchors survive.  This remains opt-in while its detector-specific
        // width and continuum controls are tuned.
        // This final-production challenger is intentionally dormant until its separate
        // force-gap experiment is enabled.  The normal over-wide switch controls the earlier
        // atom-safe policy only; coupling it to this unfinished post-solve path would confound
        // detector tuning and expose unreviewed transaction semantics.
        if( config.auto_roi_final_fitted_partition
            && auto_interferer_lines.empty()
            && !rescue_solve_this_iteration && solution.m_foreground
            && solution.m_foreground->energy_calibration()
            && solution.m_foreground->energy_calibration()->valid() )
        {
          struct RankedFinalPartition
          {
            double improvement = 0.0;
            SelectedRoiComponentPartition record;
          };
          std::vector<RankedFinalPartition> final_partition_candidates;
          const std::vector<RelActCalcAuto::RoiRange> &fitted_ranges
            = (solution.m_final_roi_ranges_in_spectrum_cal.size()
                == solution.m_final_roi_ranges.size())
              ? solution.m_final_roi_ranges_in_spectrum_cal : solution.m_final_roi_ranges;
          const std::shared_ptr<const SpecUtils::Measurement> fitted_foreground
            = solution.m_foreground;
          const std::shared_ptr<const SpecUtils::EnergyCalibration> fitted_cal
            = fitted_foreground->energy_calibration();
          detail::AutomaticRoiPolicySettings final_policy;
          final_policy.merge_tail_z = config.merge_tail_z;
          final_policy.merge_clean_gap_fwhm = config.merge_clean_gap_fwhm;
          final_policy.continuum_aicc_penalty = config.cont_order_aicc_penalty;
          final_policy.peak_core_num_fwhm = config.auto_roi_core_num_fwhm;
          final_policy.max_width_fwhm = config.auto_rel_eff_sol_max_fwhm;
          final_policy.minimum_partition_gap_fwhm = config.auto_roi_partition_min_gap_fwhm;
          final_policy.allow_clean_gap_partition_override
            = config.auto_roi_partition_allow_clean_gap_override;
          final_policy.residual_valley_max_excess_z
            = config.auto_roi_partition_residual_valley_max_excess_z;
          final_policy.global_continuum = global_cont.valid() ? &global_cont : nullptr;
          final_policy.force_partition_gap_fwhm = config.auto_roi_partition_force_gap_fwhm;
          final_policy.max_partition_children = config.auto_roi_partition_max_children;
          final_policy.stage = refinement_stage + " final fitted ROI partition";
          detail::AutomaticRoiPartitionConstraints final_constraints;
          final_constraints.lowest_energy = fitted_foreground->gamma_channel_lower( 0 );
          final_constraints.highest_energy = fitted_foreground->gamma_channel_upper(
              fitted_foreground->num_gamma_channels() - 1 );
          final_constraints.left_barrier = -std::numeric_limits<double>::infinity();
          final_constraints.min_width_fwhm = config.auto_rel_eff_sol_min_fwhm_roi;
          final_constraints.peak_core_num_fwhm = config.auto_roi_core_num_fwhm;

          for( const RelActCalcAuto::RoiRange &roi : fitted_ranges )
          {
            const double midpoint = 0.5*(roi.lower_energy + roi.upper_energy);
            const double fwhm = solution_fwhm_at_energy( midpoint );
            const double late_partition_min_width_fwhm = std::max(
                config.auto_rel_eff_sol_max_fwhm,
                config.auto_roi_final_partition_min_width_fwhm );
            if( !std::isfinite(fwhm) || !(fwhm > 0.0)
                || (((roi.upper_energy - roi.lower_energy) / fwhm)
                    <= late_partition_min_width_fwhm) )
              continue;

            detail::AutomaticRoiComponent component;
            component.lower = roi.lower_energy;
            component.upper = roi.upper_energy;
            component.first_channel = fitted_foreground->find_gamma_channel(
                static_cast<float>(component.lower) );
            component.last_channel = fitted_foreground->find_gamma_channel(
                static_cast<float>(component.upper) );
            component.continuum_type = roi.continuum_type;
            component.range_limits_type = roi.range_limits_type;
            component.joined_groups = 1;
            // `fitted_ranges` and `fitted_foreground` are in the spectrum calibration used by
            // the solve.  The uncombined no-background-subtraction peaks use that same frame;
            // `m_fit_peaks` are true-energy/public coordinates and must not be compared to these
            // bounds after an energy-calibration refinement.
            for( const PeakDef &peak : solution.m_peaks_without_back_sub )
            {
              const std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
              // The final spectrum-cal ROI list owns membership.  A fitted continuum can carry
              // slightly different fitted-object edges after the solver's calibration/refinement
              // bookkeeping, so equality of its displayed bounds is not a reliable ownership key.
              if( !continuum || (peak.mean() < roi.lower_energy)
                  || (peak.mean() > roi.upper_energy) )
                continue;
              detail::RoiAtom atom;
              atom.id = detail::next_roi_atom_id();
              atom.energy = peak.mean();
              atom.area = std::max( 0.0, peak.peakArea() );
              atom.kind = detail::RoiAtomKind::ModeledGamma;
              component.atoms.push_back( atom );
            }
            if( (config.auto_roi_final_partition_max_atoms > 0)
                && (component.atoms.size()
                    > config.auto_roi_final_partition_max_atoms) )
            {
              AutomaticRoiDecisionDiagnostic dense;
              dense.decision = AutomaticRoiDecision::MergeInseparableWide;
              dense.stage = final_policy.stage;
              dense.left_lower = roi.lower_energy;
              dense.left_upper = roi.upper_energy;
              dense.reason = "late fitted ROI partition skipped: modeled atom count exceeds configured sparse-ROI limit";
              result.automatic_roi_diagnostics.push_back( std::move(dense) );
              continue;
            }
            // Auto-search peaks that coincide with a fitted modeled atom are evidence for that
            // atom, not an intervening unmodeled feature.  Treating them as the latter vetoes
            // every candidate boundary in a final solver ROI (the same filtering is used by the
            // earlier automatic merge/reconciliation stages).
            std::vector<std::shared_ptr<const PeakDef>> final_unfit;
            for( const std::shared_ptr<const PeakDef> &unfit : current_unfit_auto_peaks )
            {
              if( !unfit || !unfit->gausPeak() )
                continue;
              const bool matches_modeled = std::any_of( std::begin(component.atoms),
                  std::end(component.atoms), [&solution_fwhm_at_energy, &unfit](
                      const detail::RoiAtom &atom ) {
                    const double atom_fwhm = solution_fwhm_at_energy( atom.energy );
                    return std::isfinite(atom_fwhm) && (atom_fwhm > 0.0)
                        && (std::fabs(unfit->mean() - atom.energy) <= 0.75*atom_fwhm);
                  } );
              if( !matches_modeled )
                final_unfit.push_back( unfit );
            }
            const detail::AutomaticRoiComponentPartitionResult partition
              = detail::partition_overwide_automatic_component( { component }, fitted_foreground,
                  solution_fwhm_at_energy, final_unfit, final_policy,
                  final_constraints );
            result.automatic_roi_diagnostics.push_back( partition.diagnostic );
            if( !partition.valid || !partition.changed || (partition.components.size() < 2) )
              continue;

            SelectedRoiComponentPartition record;
            record.id = detail::next_roi_atom_id();
            record.parent_lower = roi.lower_energy;
            record.parent_upper = roi.upper_energy;
            record.calibration = fitted_cal;
            for( const detail::AutomaticRoiComponent &child : partition.components )
            {
              RelActCalcAuto::RoiRange child_roi;
              child_roi.lower_energy = child.lower;
              child_roi.upper_energy = child.upper;
              child_roi.continuum_type = child.continuum_type;
              child_roi.range_limits_type = child.range_limits_type;
              record.children.push_back( child_roi );
            }
            RankedFinalPartition candidate;
            candidate.improvement = partition.diagnostic.one_roi_aicc
                - partition.diagnostic.two_roi_aicc;
            candidate.record = std::move(record);
            final_partition_candidates.push_back( std::move(candidate) );
          }
          std::sort( std::begin(final_partition_candidates),
                     std::end(final_partition_candidates),
            []( const RankedFinalPartition &lhs, const RankedFinalPartition &rhs ) {
              if( lhs.improvement != rhs.improvement )
                return lhs.improvement > rhs.improvement;
              return lhs.record.parent_lower < rhs.record.parent_lower;
            } );
          const size_t max_final_partitions = std::min(
              config.auto_roi_final_partition_max_proposals,
              final_partition_candidates.size() );
          for( size_t index = 0; index < max_final_partitions; ++index )
            selected_component_partitions.push_back(
                std::move(final_partition_candidates[index].record) );
        }

        // Shrink ordinary refinement ROIs to avoid interference from unfit auto-search peaks.  The
        // exact selected-component challenger above intentionally retains its pre-shrink bounds.
        const double min_fwhm_above = 0.5*config.auto_rel_eff_sol_min_fwhm_roi;
        const double min_fwhm_below = 0.5*config.auto_rel_eff_sol_min_fwhm_roi;
        shrink_rois_for_interfering_peaks( refined_rois_and_gammas, unfit_auto_peaks,
            foreground, fwhm_form, fwhm_coefficients,
            fwhm_lower_energy, fwhm_upper_energy,
            config.auto_rel_eff_cluster_num_sigma,
            min_fwhm_below, min_fwhm_above,
            orig_foreground->energy_calibration(),
            foreground->energy_calibration() );

        // Extract just the ROIs for downstream use
        std::vector<RelActCalcAuto::RoiRange> refined_rois;
        refined_rois.reserve( refined_rois_and_gammas.size() );
        for( const std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &p : refined_rois_and_gammas )
          refined_rois.push_back( p.first );

        // The nuisance curve is deliberately excluded from re-clustering above: only its
        // foreground-confirmed line neighborhoods passed the R6 transaction.  Keep those admitted
        // neighborhoods represented without allowing unconfirmed sibling lines to seed new ROIs.
        for( const double interferer_energy : auto_interferer_lines )
        {
          const bool already_covered = std::any_of( std::begin(refined_rois),
            std::end(refined_rois), [interferer_energy]( const RelActCalcAuto::RoiRange &roi ) {
              return (interferer_energy >= roi.lower_energy)
                  && (interferer_energy <= roi.upper_energy);
            } );
          if( already_covered )
            continue;

          const auto incumbent = std::find_if( std::begin(solution.m_options.rois),
            std::end(solution.m_options.rois),
            [interferer_energy]( const RelActCalcAuto::RoiRange &roi ) {
              return (interferer_energy >= roi.lower_energy)
                  && (interferer_energy <= roi.upper_energy);
            } );
          if( incumbent == std::end(solution.m_options.rois) )
            continue;

          RelActCalcAuto::RoiRange retained = *incumbent;
          std::vector<RelActCalcAuto::RoiRange> nonoverlapping;
          for( const RelActCalcAuto::RoiRange &roi : refined_rois )
          {
            if( (roi.lower_energy < retained.upper_energy)
                && (roi.upper_energy > retained.lower_energy) )
            {
              retained.lower_energy = std::min( retained.lower_energy, roi.lower_energy );
              retained.upper_energy = std::max( retained.upper_energy, roi.upper_energy );
            }else
            {
              nonoverlapping.push_back( roi );
            }
          }
          nonoverlapping.push_back( retained );
          refined_rois = std::move( nonoverlapping );
        }

        // Restore initial edge ROIs that re-clustering dropped due to rel. eff.
        // extrapolation giving near-zero expected counts.  The rel. eff. curve is
        // fitted from the interior peaks and extrapolates poorly at the energy
        // extremes.  When there is a significant auto-search peak AND the source
        // has a major gamma line there, the data should override the model.
        // To guard against NORM contamination, we require the source has a
        // substantial branching ratio there, and that the observed peak area
        // is physically plausible relative to the nearest kept ROI's peak.
        if( !refined_rois.empty() && !input_rois.empty() )
        {
          const double refined_min_center = 0.5 * (refined_rois.front().lower_energy + refined_rois.front().upper_energy);
          const double refined_max_center = 0.5 * (refined_rois.back().lower_energy + refined_rois.back().upper_energy);

          for( const RelActCalcAuto::RoiRange &init_roi : input_rois )
          {
            const double init_center = 0.5 * (init_roi.lower_energy + init_roi.upper_energy);

            // Only consider ROIs at the energy extremes
            if( (init_center >= refined_min_center) && (init_center <= refined_max_center) )
              continue;

            // Check not already covered by a refined ROI
            bool covered = false;
            for( const RelActCalcAuto::RoiRange &ref_roi : refined_rois )
            {
              if( (init_center >= ref_roi.lower_energy) && (init_center <= ref_roi.upper_energy) )
              {
                covered = true;
                break;
              }
            }
            if( covered )
              continue;

            // Find the best auto-search peak in this dropped ROI
            std::shared_ptr<const PeakDef> best_peak;
            for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
            {
              if( (peak->mean() >= init_roi.lower_energy)
                && (peak->mean() <= init_roi.upper_energy)
                && (!best_peak || (peak->peakArea() > best_peak->peakArea())) )
              {
                best_peak = peak;
              }
            }

            if( !best_peak )
              continue;

            const double best_peak_area = best_peak->peakArea();

            // Require a statistically real detection rather than the former absolute 100-count
            // floor (which could not transfer across live-times).  Prefer the auto-search fit's
            // own area uncertainty; fall back to a sideband-continuum significance if absent.
            double detection_z = 0.0;
            if( best_peak->peakAreaUncert() > 0.0 )
            {
              detection_z = best_peak_area / best_peak->peakAreaUncert();
            }else
            {
              const double pk_mean = best_peak->mean();
              const double pk_fwhm = best_peak->fwhm();
              const double s_est = 0.7607 * best_peak_area;  // Gaussian fraction within +/-1 FWHM

              // R1 step 2: prefer the single global SNIP continuum for B here (this is the weakest
              // local estimator - no signal subtraction/relocation); fall back to the local estimate,
              // then to gross counts, whenever the global provider is unavailable.
              double b_est;
              if( global_cont.valid() )
              {
                b_est = global_cont.integral( pk_mean - pk_fwhm, pk_mean + pk_fwhm );
              }else
              {
                const detail::LocalContinuumEstimate edge_cont = detail::estimate_local_continuum(
                    foreground, pk_mean - 2.0*pk_fwhm, pk_mean + 2.0*pk_fwhm, pk_fwhm, 0.5 );
                b_est = edge_cont.valid
                    ? edge_cont.integral( pk_mean - pk_fwhm, pk_mean + pk_fwhm )
                    : foreground->gamma_integral( pk_mean - pk_fwhm, pk_mean + pk_fwhm );
              }
              detection_z = s_est / std::sqrt( std::max( 1.0, s_est + b_est ) );
            }

            if( detection_z < sm_edge_roi_restore_min_z )
              continue;

            // Find the nearest kept ROI (by center energy) that has an auto-search peak
            double nearest_dist = std::numeric_limits<double>::max();
            double nearest_roi_peak_area = 0.0;
            double nearest_roi_center = 0.0;

            for( const RelActCalcAuto::RoiRange &ref_roi : refined_rois )
            {
              const double ref_center = 0.5 * (ref_roi.lower_energy + ref_roi.upper_energy);
              const double dist = std::abs( ref_center - init_center );

              // Find the best auto-search peak in this refined ROI
              double ref_peak_area = 0.0;
              for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
              {
                if( (peak->mean() >= ref_roi.lower_energy)
                  && (peak->mean() <= ref_roi.upper_energy)
                  && (peak->peakArea() > ref_peak_area) )
                {
                  ref_peak_area = peak->peakArea();
                }
              }

              if( (ref_peak_area > 0.0) && (dist < nearest_dist) )
              {
                nearest_dist = dist;
                nearest_roi_peak_area = ref_peak_area;
                nearest_roi_center = ref_center;
              }
            }//for( each refined ROI )

            if( nearest_roi_peak_area <= 0.0 )
              continue;

            // Find the max source BR in the dropped ROI and in the nearest kept ROI
            double dropped_max_br = 0.0;
            double nearest_max_br = 0.0;

            for( const RelActCalcAuto::NucInputInfo &src : sources )
            {
              const std::vector<SandiaDecay::EnergyRatePair> photons
                = get_source_photons( src.source, 1.0, src.age );

              for( const SandiaDecay::EnergyRatePair &p : photons )
              {
                if( (p.energy >= init_roi.lower_energy) && (p.energy <= init_roi.upper_energy) )
                  dropped_max_br = std::max( dropped_max_br, p.numPerSecond );
              }

              // Find the best BR near the nearest ROI center (within half the ROI width)
              const double half_width = 0.25 * (refined_rois.front().upper_energy - refined_rois.front().lower_energy + 1.0);
              for( const SandiaDecay::EnergyRatePair &p : photons )
              {
                if( std::abs( p.energy - nearest_roi_center ) < half_width )
                  nearest_max_br = std::max( nearest_max_br, p.numPerSecond );
              }
            }//for( each source )

            // Require the source has a non-trivial BR in the dropped ROI
            if( dropped_max_br < 0.001 )
              continue;

            // Sanity check: the observed area ratio should not wildly exceed the BR ratio.
            // Shielding reduces the low-energy-to-high-energy ratio, never increases it.
            // A factor of 10 allows for efficiency and calibration variations.
            if( (nearest_max_br > 0.0) && (dropped_max_br > 0.0) )
            {
              const double area_ratio = best_peak_area / nearest_roi_peak_area;
              const double br_ratio = dropped_max_br / nearest_max_br;

              if( area_ratio > br_ratio * 10.0 )
              {
                if( should_debug_print() )
                {
                  std::cout << "Edge ROI [" << init_roi.lower_energy << ", " << init_roi.upper_energy
                       << "] keV not restored: area_ratio=" << area_ratio
                       << " > br_ratio*10=" << (br_ratio * 10.0) << " - likely NORM" << std::endl;
                }
                continue;
              }
            }

            refined_rois.push_back( init_roi );
            if( should_debug_print() )
            {
              std::cout << "Restored edge ROI [" << init_roi.lower_energy << ", " << init_roi.upper_energy
                   << "] keV: auto peak area=" << best_peak_area
                   << ", dropped_br=" << dropped_max_br << std::endl;
            }
          }//for( each initial ROI )

          std::sort( std::begin(refined_rois), std::end(refined_rois),
            []( const RelActCalcAuto::RoiRange &a, const RelActCalcAuto::RoiRange &b ){
              return a.lower_energy < b.lower_energy;
          } );

          // Any restored-edge overlap is reconciled by the shared policy after all late additions.
        }//if( check for dropped edge ROIs )

        // Change 3 (auto stage): carry the found-peak guarantee into refinement.  Re-clustering uses
        // the same predicted-signal keep-gate as the manual stage (at the stricter auto threshold),
        // so it can drop an ROI whose data-confirmed source line the search already found (e.g. Co58
        // 810, whose predicted z falls just under the auto bar).  For every INPUT ROI that still holds
        // a significant auto-search peak, re-seed a tight ROI if the refined set no longer covers it.
        // Deliberately outside the "dropped edge ROIs" block above so it also fires when re-clustering
        // returned nothing at all.  [architecture review 2026-07-18]
        {
          std::vector<std::pair<double,double>> found_energy_fwhm;
          for( const RelActCalcAuto::RoiRange &in_roi : input_rois )
          {
            std::shared_ptr<const PeakDef> best;
            for( const std::shared_ptr<const PeakDef> &pk : auto_search_peaks )
            {
              if( !pk || !pk->gausPeak()
                  || (pk->mean() < in_roi.lower_energy) || (pk->mean() > in_roi.upper_energy) )
                continue;
              const double z = (pk->amplitudeUncert() > 0.0)
                              ? (pk->amplitude() / pk->amplitudeUncert()) : 0.0;
              if( z < config.roi_significance_z )
                continue;  // only data-significant found peaks pin an ROI (self-limiting)
              if( !best || (pk->amplitude() > best->amplitude()) )
                best = pk;
            }//for( loop over auto_search_peaks )

            if( best )
              found_energy_fwhm.emplace_back( best->mean(), best->fwhm() );
          }//for( const RelActCalcAuto::RoiRange &in_roi : input_rois )

          const size_t num_seeded = seed_tight_rois_for_found_peaks(
              refined_rois, found_energy_fwhm, sm_found_peak_roi_half_num_fwhm,
              min_valid_energy, max_valid_energy, use_automatic_roi_policy );

          if( num_seeded )
          {
            // rois_are_similar (and the downstream overlap resolver) expect energy-sorted ROIs.
            std::sort( std::begin(refined_rois), std::end(refined_rois),
              []( const RelActCalcAuto::RoiRange &a, const RelActCalcAuto::RoiRange &b ){
                return a.lower_energy < b.lower_energy;
            } );

            if( should_debug_print() )
              std::cout << "Auto re-cluster: re-seeded " << num_seeded
                        << " tight ROI(s) for significant found peak(s)" << std::endl;
          }
        }

        // Keep previously admitted rescue ROIs present while the ordinary re-clustering converges.
        // They still pass through the final insignificant-ROI filter below.
        for( const RelActCalcAuto::RoiRange &rescued : rescued_roi_ranges )
        {
          const double center = 0.5 * (rescued.lower_energy + rescued.upper_energy);
          const bool covered = std::any_of( std::begin(refined_rois), std::end(refined_rois),
            [center]( const RelActCalcAuto::RoiRange &roi ) {
              return (center >= roi.lower_energy) && (center <= roi.upper_energy);
            } );
          if( !covered )
            refined_rois.push_back( rescued );
        }

        // R2: one bounded fit-then-prune admission pass, on the first fitted re-clustering only.
        // The normal accepted set above is untouched; this considers only provenance-preserving
        // clusters that missed the keep threshold by less than the fixed rescue fraction.
        if( bounded_rescue_enabled() && !rescue_attempted && (iter == 0) )
        {
          rescue_attempted = true;
          const std::vector<RelActCalcAuto::RoiRange> pre_rescue_refined_rois = refined_rois;
          try
          {
#if( PERFORM_DEVELOPER_CHECKS )
          if( sm_force_next_rescue_admission_failure_for_test )
          {
            sm_force_next_rescue_admission_failure_for_test = false;
            throw std::runtime_error( "forced rescue-admission failure for developer test" );
          }
#endif

          const auto source_is_requested = [&sources]( const RelActCalcAuto::SrcVariant &source ) {
            return std::any_of( std::begin(sources), std::end(sources),
              [&source]( const RelActCalcAuto::NucInputInfo &input ) {
                return input.source == source;
              } );
          };

          // Re-cluster and reclassify only the requested-source portion of the rejected
          // predictions.  The fitted solution may now contain an R6 nuisance curve; retaining the
          // mixed cluster's old counts/z would let nuisance counts promote a sub-marginal source
          // line into rescue.  This second pass uses the identical production clustering and keep
          // gate with the fitted resolution, but cannot change the normal accepted ROI set.
          std::vector<PredictedGamma> requested_rejected_predictions;
          for( const PredictedGamma &gamma : fitted_cluster_predictions )
          {
            if( source_is_requested(gamma.source) )
              requested_rejected_predictions.push_back( gamma );
          }
          std::vector<MarginalRejectedCluster> requested_only_marginals;
          if( !requested_rejected_predictions.empty() )
          {
            GammaClusteringSettings r2_settings = auto_settings;
            r2_settings.use_automatic_roi_policy = false;
            static_cast<void>( cluster_gammas_to_rois( {}, {}, foreground,
                fwhm_form, fwhm_coefficients, fwhm_lower_energy, fwhm_upper_energy,
                min_valid_energy, max_valid_energy, r2_settings, current_unfit_auto_peaks,
                &requested_only_marginals, &requested_rejected_predictions,
                &solution_fwhm_at_energy, nullptr, "R2 requested-only reclassification",
                nullptr, nullptr ) );
          }
          marginal_rejects = std::move( requested_only_marginals );

          // Each curve needs at least two distinct fitted requested-source energies before its
          // shape has a defensible interpolation span for rescue.
          std::vector<std::vector<double>> fitted_curve_energies(
              solution.m_options.rel_eff_curves.size() );
          for( const PeakDef &peak : solution.m_peaks_without_back_sub )
          {
            if( !peak.continuum() || (peak.mean() < peak.continuum()->lowerEnergy())
                || (peak.mean() > peak.continuum()->upperEnergy()) )
              continue;
            const std::vector<RelActCalcAuto::RoiRange> &fitted_ranges
              = solution.m_final_roi_ranges_in_spectrum_cal.empty()
                ? solution.m_final_roi_ranges : solution.m_final_roi_ranges_in_spectrum_cal;
            const bool covered_by_fitted_data = std::any_of( std::begin(fitted_ranges),
              std::end(fitted_ranges), [&peak]( const RelActCalcAuto::RoiRange &roi ) {
                return (peak.mean() >= roi.lower_energy) && (peak.mean() <= roi.upper_energy);
              } );
            if( !covered_by_fitted_data )
              continue;

            RelActCalcAuto::SrcVariant peak_source;
            if( peak.parentNuclide() )
              peak_source = peak.parentNuclide();
            else if( peak.xrayElement() )
              peak_source = peak.xrayElement();
            else if( peak.reaction() )
              peak_source = peak.reaction();
            else
              continue;
            if( !source_is_requested(peak_source) )
              continue;

            const double energy = peak.hasSourceGammaAssigned()
                ? peak.gammaParticleEnergy() : peak.mean();
            for( size_t curve_index = 0;
                 curve_index < solution.m_options.rel_eff_curves.size(); ++curve_index )
            {
              const RelActCalcAuto::RelEffCurveInput &curve
                = solution.m_options.rel_eff_curves[curve_index];
              const bool on_curve = std::any_of( std::begin(curve.nuclides),
                std::end(curve.nuclides),
                [&peak_source]( const RelActCalcAuto::NucInputInfo &input ) {
                  return input.source == peak_source;
                } );
              if( on_curve )
                fitted_curve_energies[curve_index].push_back( energy );
            }
          }
          for( std::vector<double> &energies : fitted_curve_energies )
          {
            std::sort( std::begin(energies), std::end(energies) );
            energies.erase( std::unique(std::begin(energies), std::end(energies),
              []( const double lhs, const double rhs ) {
                return std::fabs(lhs - rhs) < 0.05;
              } ), std::end(energies) );
          }

          struct RankedMarginal
          {
            size_t index;
            double significance;
            double energy;
          };
          std::vector<RankedMarginal> ranked_marginals;
          const SandiaDecay::SandiaDecayDataBase * const decay_db
            = DecayDataBaseServer::database();
          const std::shared_ptr<const SpecUtils::EnergyCalibration> original_cal
            = orig_foreground->energy_calibration();
          const std::shared_ptr<const SpecUtils::EnergyCalibration> working_cal
            = foreground->energy_calibration();
          const auto unfit_peak_in_working_cal
            = [&original_cal, &working_cal]( const PeakDef &peak ) -> double {
              if( !original_cal || !working_cal || (original_cal == working_cal) )
                return peak.mean();
              const double channel = original_cal->channel_for_energy( peak.mean() );
              return working_cal->energy_for_channel( channel );
            };

          for( size_t marginal_index = 0;
               marginal_index < marginal_rejects.size(); ++marginal_index )
          {
            MarginalRejectedCluster &marginal = marginal_rejects[marginal_index];
            bool guarded = false;
            for( const PredictedGamma &gamma : marginal.predicted_gammas )
            {
              assert( source_is_requested(gamma.source) );

              const double fwhm = solution_fwhm_at_energy( gamma.energy );
              if( !std::isfinite(fwhm) || !(fwhm > 0.0)
                  || (gamma.rel_eff_curve_index >= fitted_curve_energies.size()) )
              {
                guarded = true;
                break;
              }

              const std::vector<double> &span
                = fitted_curve_energies[gamma.rel_eff_curve_index];
              if( (span.size() < 2) || (gamma.energy < span.front())
                  || (gamma.energy > span.back()) )
              {
                guarded = true;
                break;
              }

              if( std::any_of( std::begin(refined_rois), std::end(refined_rois),
                    [&gamma]( const RelActCalcAuto::RoiRange &roi ) {
                      return (gamma.energy >= roi.lower_energy) && (gamma.energy <= roi.upper_energy);
                    } ) )
              {
                guarded = true;
                break;
              }

              const auto near_energy = [&]( const double energy ) {
                return std::fabs(energy - gamma.energy)
                    < (sm_rescue_guard_num_fwhm * fwhm);
              };
              if( std::any_of( std::begin(unfit_auto_peaks), std::end(unfit_auto_peaks),
                    [&]( const std::shared_ptr<const PeakDef> &peak ) {
                      return peak && near_energy(unfit_peak_in_working_cal(*peak));
                    } )
                  || std::any_of( std::begin(interferer_guard_energies),
                    std::end(interferer_guard_energies), near_energy ) )
              {
                guarded = true;
                break;
              }

              const SandiaDecay::Nuclide * const source_nuclide
                = RelActCalcAuto::nuclide( gamma.source );
              for( const StrongNormGammaLine &line : sk_strong_norm_gamma_lines )
              {
                if( !near_energy(line.energy) )
                  continue;
                const SandiaDecay::Nuclide * const norm_parent
                  = (decay_db && line.parent_symbol) ? decay_db->nuclide(line.parent_symbol) : nullptr;
                if( !source_nuclide || (norm_parent != source_nuclide) )
                {
                  guarded = true;
                  break;
                }
              }
              if( guarded )
                break;
            }//for( marginal predicted gammas )

            if( guarded || marginal.predicted_gammas.empty() )
              continue;
            const double energy = marginal.predicted_gammas.front().energy;
            ranked_marginals.push_back( RankedMarginal{
                marginal_index, marginal.keep_significance, energy } );
          }//for( marginal rejects )

          std::sort( std::begin(ranked_marginals), std::end(ranked_marginals),
            []( const RankedMarginal &lhs, const RankedMarginal &rhs ) {
              if( lhs.significance != rhs.significance )
                return lhs.significance > rhs.significance;
              return lhs.energy < rhs.energy;
            } );

          GammaClusteringSettings rescue_settings = auto_settings;
          rescue_settings.keep_significance_z = 0.0;
          rescue_settings.use_automatic_roi_policy = false;
          size_t inspected = 0;
          for( const RankedMarginal &ranked : ranked_marginals )
          {
            if( proposed_rescued_rois.size() >= sm_max_rescued_rois )
              break;

            const MarginalRejectedCluster &marginal = marginal_rejects[ranked.index];
            const std::vector<std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo>> candidates
              = cluster_gammas_to_rois( {}, {}, foreground, fwhm_form, fwhm_coefficients,
                  fwhm_lower_energy, fwhm_upper_energy, min_valid_energy, max_valid_energy,
                  rescue_settings, current_unfit_auto_peaks, nullptr,
                  &marginal.predicted_gammas, &solution_fwhm_at_energy, nullptr,
                  "R2 bounded candidate construction", nullptr, nullptr );

            for( const std::pair<RelActCalcAuto::RoiRange, ClusteredGammaInfo> &candidate : candidates )
            {
              if( (inspected >= sm_max_rescued_rois)
                  || (proposed_rescued_rois.size() >= sm_max_rescued_rois) )
                break;
              const RelActCalcAuto::RoiRange &roi = candidate.first;
              const bool overlaps = std::any_of( std::begin(refined_rois),
                std::end(refined_rois), [&roi]( const RelActCalcAuto::RoiRange &accepted ) {
                  return (roi.lower_energy < accepted.upper_energy)
                      && (roi.upper_energy > accepted.lower_energy);
                } ) || std::any_of( std::begin(proposed_rescued_rois),
                std::end(proposed_rescued_rois), [&roi]( const RelActCalcAuto::RoiRange &accepted ) {
                  return (roi.lower_energy < accepted.upper_energy)
                      && (roi.upper_energy > accepted.lower_energy);
                } ) || std::any_of( std::begin(rescue_rejected_ranges),
                std::end(rescue_rejected_ranges), [&roi]( const RelActCalcAuto::RoiRange &rejected ) {
                  return (roi.lower_energy < rejected.upper_energy)
                      && (roi.upper_energy > rejected.lower_energy);
                } );
              if( overlaps )
                continue;

              const double roi_center = 0.5 * (roi.lower_energy + roi.upper_energy);
              const bool incumbent_contains_center = std::any_of(
                std::begin(solution.m_options.rois), std::end(solution.m_options.rois),
                [roi_center]( const RelActCalcAuto::RoiRange &incumbent ) {
                  return (roi_center >= incumbent.lower_energy)
                      && (roi_center <= incumbent.upper_energy);
                } );
              const bool partially_overlaps_incumbent = !incumbent_contains_center
                && std::any_of( std::begin(solution.m_options.rois),
                  std::end(solution.m_options.rois), [&roi]( const RelActCalcAuto::RoiRange &incumbent ) {
                    return (roi.lower_energy < incumbent.upper_energy)
                        && (roi.upper_energy > incumbent.lower_energy);
                  } );
              if( partially_overlaps_incumbent )
              {
                rescue_rejected_ranges.push_back( roi );
                continue;
              }

              ++inspected;

              std::shared_ptr<PeakContinuum> continuum = std::make_shared<PeakContinuum>();
              continuum->setType( roi.continuum_type );
              continuum->setRange( roi.lower_energy, roi.upper_energy );
              continuum->setParameters( roi_center,
                  std::vector<double>(PeakContinuum::num_parameters(roi.continuum_type), 0.0), {} );
              std::vector<PeakDef> provisional_peaks;
              for( size_t i = 0; i < candidate.second.gamma_energies.size(); ++i )
              {
                const double energy = candidate.second.gamma_energies[i];
                const double sigma = solution_fwhm_at_energy(energy)
                    / PhysicalUnits::fwhm_nsigma;
                PeakDef peak( energy, sigma, candidate.second.gamma_amplitudes[i] );
                peak.setContinuum( continuum );
                provisional_peaks.push_back( std::move(peak) );
              }

              const RoiSignificanceResult significance = compute_roi_chi2_significance(
                  roi, provisional_peaks, foreground, config.roi_significance_z,
                  /*include_peak_count_significance=*/false,
                  /*same_continuum_family_for_null=*/true );
              if( significance.has_significant_peaks )
                proposed_rescued_rois.push_back( roi );
              else
                rescue_rejected_ranges.push_back( roi );
            }//for( rebuilt candidate ROIs )
          }//for( ranked marginal clusters )

          if( !proposed_rescued_rois.empty() )
          {
            // Make rescue a narrow transaction on the successful incumbent geometry.  Using the
            // ordinary re-clustered set here can simultaneously discard unrelated, strongly fitted
            // source anchors; starting from the incumbent changes only the admitted ranges and
            // makes rollback/anchor comparison meaningful.
            refined_rois = solution.m_options.rois;
            for( const RelActCalcAuto::RoiRange &rescued : proposed_rescued_rois )
            {
              const double center = 0.5 * (rescued.lower_energy + rescued.upper_energy);
              const bool covered = std::any_of( std::begin(refined_rois),
                std::end(refined_rois), [center]( const RelActCalcAuto::RoiRange &roi ) {
                  return (center >= roi.lower_energy) && (center <= roi.upper_energy);
                } );
              if( !covered )
                refined_rois.push_back( rescued );
            }
          }
          }
          catch( const std::exception &error )
          {
            refined_rois = pre_rescue_refined_rois;
            proposed_rescued_rois.clear();
            rescue_transaction_failed = true;
            result.warnings.push_back( "The bounded marginal-line rescue admission failed; retained"
              " the successful incumbent source fit (" + std::string(error.what()) + ")." );
          }
          catch( ... )
          {
            refined_rois = pre_rescue_refined_rois;
            proposed_rescued_rois.clear();
            rescue_transaction_failed = true;
            result.warnings.push_back( "The bounded marginal-line rescue admission failed; retained"
              " the successful incumbent source fit." );
          }
        }//one R2 rescue admission pass

        if( rescue_transaction_failed )
          break;

        if( !refined_rois.empty() )
        {
          std::sort( std::begin(refined_rois), std::end(refined_rois),
            []( const RelActCalcAuto::RoiRange &lhs, const RelActCalcAuto::RoiRange &rhs ) {
              return lhs.lower_energy < rhs.lower_energy;
            } );
        }

        // Retained mixed user ROIs are materialized by channel in the current solve calibration
        // on every pass, even when automatic re-clustering found nothing.
        for( const RelActCalcAuto::RoiRange &protected_roi : current_protected_mixed_rois )
        {
          const bool present = std::any_of( std::begin(refined_rois), std::end(refined_rois),
            [&protected_roi]( const RelActCalcAuto::RoiRange &roi ) {
              return (std::fabs(roi.lower_energy - protected_roi.lower_energy) < 1.0e-6)
                  && (std::fabs(roi.upper_energy - protected_roi.upper_energy) < 1.0e-6);
            } );
          if( !present )
            refined_rois.push_back( protected_roi );
        }

        if( refined_rois.empty() )
        {
          // If we lost all ROIs and are not already using a PhysicalModel on the sources
          // curve, try re-fitting with one.  This is similar to the "desperation" approach above.
          const bool using_physical_model = (sources_rel_eff_index >= 0)
            && (solution.m_options.rel_eff_curves[sources_rel_eff_index].rel_eff_eqn_type
                == RelActCalc::RelEffEqnForm::FramPhysicalModel);

          if( !using_physical_model && (sources_rel_eff_index >= 0) )
          {
            if( should_debug_print() )
              std::cerr << "Lost all ROIs, trying PhysicalModel as desperation attempt..." << std::endl;

            RelActCalcAuto::Options desperation_opts = options;
            RelActCalcAuto::RelEffCurveInput &curve = desperation_opts.rel_eff_curves[sources_rel_eff_index];
            curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
            curve.rel_eff_eqn_order = 0;
            curve.phys_model_self_atten = nullptr;
            curve.phys_model_external_atten.clear();

            // Apply external shielding if conditions are met
            // Note: We use desperation_opts.rois (not options.rois) since they may differ at this point
            if( should_use_desperation_shielding( config.desperation_phys_model_atomic_number, desperation_opts.rois ) )
            {
              try
              {
                std::shared_ptr<RelActCalc::PhysicalModelShieldInput> shield
                  = create_desperation_shielding( config.desperation_phys_model_atomic_number,
                                                  config.desperation_phys_model_areal_density_g_per_cm2 );
                curve.phys_model_external_atten.push_back( shield );

                if( should_debug_print() )
                {
                  std::cerr << "Second desperation attempt: using external shielding with AN="
                            << config.desperation_phys_model_atomic_number
                            << ", starting AD=" << config.desperation_phys_model_areal_density_g_per_cm2
                            << " g/cm2" << std::endl;
                }
              }
              catch( const std::exception &e )
              {
                if( should_debug_print() )
                  std::cerr << "Failed to create desperation shielding: " << e.what() << std::endl;
                // Continue without shielding
              }
            }
            else
            {
              if( should_debug_print() )
                std::cerr << "Second desperation attempt: not using external shielding" << std::endl;
            }

            curve.phys_model_corr.corr_fcn = (desperation_opts.rois.size() > 2)
                            ? RelActCalc::PhysModelCorrFcn::Hoerl : RelActCalc::PhysModelCorrFcn::None;

            add_floating_511_peak_if_appropriate( desperation_opts, sources, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy );
            add_escape_peak_floating_peaks_if_appropriate( desperation_opts, auto_search_peaks, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy, config );

            if( auto_interferer_lines.empty() )
            {
              resolve_automatic_overlapping_rois( desperation_opts.rois,
                  desperation_opts.floating_peaks, foreground,
                  global_cont.valid() ? &global_cont : nullptr, solution_fwhm_at_energy,
                  current_unfit_auto_peaks, config,
                  "refinement desperation escape/overlap finalization",
                  &result.automatic_roi_diagnostics, current_protected_mixed_rois,
                  current_modeled_peak_candidates, use_automatic_roi_policy );
            }else
            {
              resolve_overlapping_rois( desperation_opts.rois, desperation_opts.floating_peaks );
              ensure_min_channel_gap( desperation_opts.rois, foreground->energy_calibration() );
            }
            remove_floating_peaks_without_roi( desperation_opts );

            RelActCalcAuto::RelActAutoSolution desperation_solution = RelActCalcAuto::solve(
              desperation_opts, foreground, background, drf, auto_search_peaks, det_type
            );

            if( RelActCalcAuto::RelActAutoSolution::is_usable_status(desperation_solution.m_status)
                && (reduced_chi2(desperation_solution) < reduced_chi2(solution)) )
            {
              if( should_debug_print() )
                std::cerr << "PhysicalModel desperation solution succeeded and improved chi2/dof" << std::endl;
              
              if( !num_extra_allowed ) //Allow an extra iteration since we changed
                num_extra_allowed += 1;
              
              solution = desperation_solution;
              // Continue with another iteration using the new solution
              continue;
            }
            else
            {
              if( should_debug_print() )
                std::cerr << "PhysicalModel desperation solution did not improve result" << std::endl;
            }
          }//if( !using_physical_model )

          std::cerr << "Have lost all ROIs!  Halting iterations to refine solution." << std::endl;
          break;
        }//if( refined_rois.empty() )

        if( auto_interferer_lines.empty() )
        {
          resolve_automatic_overlapping_rois( refined_rois, solution.m_options.floating_peaks,
              foreground, global_cont.valid() ? &global_cont : nullptr,
              solution_fwhm_at_energy, current_unfit_auto_peaks, config,
              "refinement edge/found/rescue reconciliation",
              &result.automatic_roi_diagnostics, current_protected_mixed_rois,
              current_modeled_peak_candidates, use_automatic_roi_policy );
        }else
        {
          resolve_overlapping_rois( refined_rois, solution.m_options.floating_peaks );
          ensure_min_channel_gap( refined_rois, foreground->energy_calibration() );
          AutomaticRoiDecisionDiagnostic r6_bypass;
          r6_bypass.decision = AutomaticRoiDecision::R6LegacyBypass;
          r6_bypass.stage = "R6 refinement neighborhood retention";
          r6_bypass.reason = "accepted R6 nuisance geometry remains legacy transactional behavior";
          result.automatic_roi_diagnostics.push_back( r6_bypass );
        }
        
        // Debug output: print refined ROIs with expected counts
        if( should_debug_print() )
        {
          std::cout << "Iteration " << iter << " refined ROIs:" << std::endl;
          for( size_t roi_idx = 0; roi_idx < refined_rois.size(); ++roi_idx )
          {
            const RelActCalcAuto::RoiRange &roi = refined_rois[roi_idx];
            const double roi_data_counts = foreground->gamma_integral(
                static_cast<float>(roi.lower_energy), static_cast<float>(roi.upper_energy) );
            std::cout << "  ROI " << roi_idx << ": [" << roi.lower_energy << " - " << roi.upper_energy
                 << "] keV, width=" << (roi.upper_energy - roi.lower_energy)
                 << " keV, data_counts=" << roi_data_counts << std::endl;
          }
        }

        // Check if ROIs changed significantly - if not, stop iterating
        // A final fitted-ROI partition is intentionally proposed from the *incumbent* solver
        // geometry.  It commonly leaves the ordinary predicted-gamma re-clustering unchanged,
        // which used to take this early convergence exit before the component-local transaction
        // below could test the proposal.  Keep the established early exit for ordinary
        // convergence, but let an explicit bounded partition challenger reach its isolated
        // accept-or-rollback solve.  That solve starts from the incumbent and changes only the
        // selected component, so it cannot make an otherwise-similar re-cluster mutate globally.
        if( rois_are_similar( refined_rois, solution.m_options.rois )
            && selected_component_partitions.empty() )
        {
          if( !proposed_rescued_rois.empty() )
          {
            rescued_roi_ranges.insert( std::end(rescued_roi_ranges),
                std::begin(proposed_rescued_rois), std::end(proposed_rescued_rois) );
            result.warnings.push_back( "Retained "
              + std::to_string(proposed_rescued_rois.size())
              + " marginal source ROI(s) through the bounded fit-then-prune rescue; the"
                " successful incumbent already used the same ROI geometry, and final ROI"
                " significance filtering remains authoritative." );
          }
          if( should_debug_print() )
            std::cout << "Iteration " << iter << ": ROIs are similar, stopping refinement" << std::endl;
          break;
        }

        if( should_debug_print() )
          std::cout << "Iteration " << iter << ": trying " << refined_rois.size() << " refined ROIs" << std::endl;

        // Re-run RelActAuto with refined ROIs.  When a rescue challenger is present, everything
        // from option preparation through post-solve evaluation is transactional: no exception in
        // challenger-only machinery may replace the already-successful incumbent.
        const bool rescue_transaction_requested = !proposed_rescued_rois.empty();
        try
        {
        RelActCalcAuto::Options refined_options = solution.m_options;
        // Apply DoNotUseExistingRois filtering to refined ROIs as well, so iterative
        // re-clustering can't produce ROIs that land on existing user peaks' locations.
        refined_options.rois = filter_rois_for_existing( refined_rois, foreground );

        add_floating_511_peak_if_appropriate( refined_options, sources, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy );
        add_escape_peak_floating_peaks_if_appropriate( refined_options, auto_search_peaks, fit_norm_peaks, det_type, min_valid_energy, max_valid_energy, config );

        if( auto_interferer_lines.empty() )
        {
          resolve_automatic_overlapping_rois( refined_options.rois,
              refined_options.floating_peaks, foreground,
              global_cont.valid() ? &global_cont : nullptr, solution_fwhm_at_energy,
              current_unfit_auto_peaks, config,
              "refined solve escape/overlap finalization",
              &result.automatic_roi_diagnostics, current_protected_mixed_rois,
              current_modeled_peak_candidates, use_automatic_roi_policy );
        }else
        {
          resolve_overlapping_rois( refined_options.rois, refined_options.floating_peaks );
          ensure_min_channel_gap( refined_options.rois, foreground->energy_calibration() );
        }
        remove_floating_peaks_without_roi( refined_options );

        if( !proposed_rescued_rois.empty() )
        {
          proposed_rescued_rois.erase( std::remove_if(
              std::begin(proposed_rescued_rois), std::end(proposed_rescued_rois),
              [&refined_options]( const RelActCalcAuto::RoiRange &rescued ) {
                const double center = 0.5 * (rescued.lower_energy + rescued.upper_energy);
                return !std::any_of( std::begin(refined_options.rois),
                  std::end(refined_options.rois), [center]( const RelActCalcAuto::RoiRange &roi ) {
                    return (center >= roi.lower_energy) && (center <= roi.upper_energy);
                  } );
              } ), std::end(proposed_rescued_rois) );
          rescue_solve_this_iteration = !proposed_rescued_rois.empty();
        }

        RelActCalcAuto::RelActAutoSolution refined_solution;
        try
        {
          refined_solution = RelActCalcAuto::solve(
              refined_options, foreground, background, drf, auto_search_peaks, det_type );
        }
        catch( const std::exception &error )
        {
          if( !rescue_solve_this_iteration )
            throw;
          result.warnings.push_back( "The bounded marginal-line rescue solve threw; retained"
            " the successful incumbent source fit (" + std::string(error.what()) + ")." );
          break;
        }
        catch( ... )
        {
          if( !rescue_solve_this_iteration )
            throw;
          result.warnings.push_back( "The bounded marginal-line rescue solve threw; retained"
            " the successful incumbent source fit." );
          break;
        }

        if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(refined_solution.m_status) )
        {
          if( rescue_solve_this_iteration )
          {
            result.warnings.push_back( "The bounded marginal-line rescue solve failed; retained"
              " the successful incumbent source fit." );
          }
          if( should_debug_print() )
            std::cout << "Iteration " << iter << " failed: " << refined_solution.m_error_message << std::endl;
          break;
        }

        const std::vector<double> anchor_exclusions = rescue_solve_this_iteration
            ? std::vector<double>() : auto_interferer_lines;
        bool recovered_anchor_failure = false;
        if( !recovered_source_anchors.empty() )
        {
          const bool predicted_preserved = (sources_rel_eff_index >= 0)
            && std::all_of( std::begin(recovered_source_anchors),
                std::end(recovered_source_anchors),
                [&refined_solution, sources_rel_eff_index, &orig_foreground](
                    const RelActCalcManual::GenericPeakInfo &anchor ) {
                  return predicted_source_anchor_counts( refined_solution,
                      static_cast<size_t>(sources_rel_eff_index), anchor,
                      orig_foreground->live_time() ) > sm_keep_gate_min_est_counts;
                } );
          const bool fitted_preserved = significant_requested_source_anchors_preserved(
              refined_solution, sources, recovered_source_anchors,
              config.roi_significance_z, solution_fwhm_at_energy );
          if( should_debug_print() && (!predicted_preserved || !fitted_preserved) )
          {
            const std::vector<PeakDef> fitted_lines
              = distinct_significant_requested_source_peaks( refined_solution, sources,
                  config.roi_significance_z, solution_fwhm_at_energy );
            std::cerr << "Recovered source-anchor check: predicted=" << predicted_preserved
                      << ", observed=" << fitted_preserved << std::endl;
            for( const RelActCalcManual::GenericPeakInfo &anchor : recovered_source_anchors )
            {
              const double center = std::isfinite(anchor.m_mean) ? anchor.m_mean : anchor.m_energy;
              const bool observed = std::any_of( std::begin(fitted_lines), std::end(fitted_lines),
                [center, &anchor, &solution_fwhm_at_energy]( const PeakDef &line ) {
                  const double line_fwhm = solution_fwhm_at_energy( line.mean() );
                  const double anchor_fwhm = (std::isfinite(anchor.m_fwhm)
                      && (anchor.m_fwhm > 0.0)) ? anchor.m_fwhm : line_fwhm;
                  return std::isfinite(line_fwhm) && (line_fwhm > 0.0)
                      && (std::fabs(line.mean() - center)
                          < 0.5*(anchor_fwhm + line_fwhm));
                } );
              const double predicted = (sources_rel_eff_index >= 0)
                ? predicted_source_anchor_counts( refined_solution,
                    static_cast<size_t>(sources_rel_eff_index), anchor,
                    orig_foreground->live_time() ) : 0.0;
              if( !observed || !(predicted > sm_keep_gate_min_est_counts) )
                std::cerr << "  missing clean anchor at " << anchor.m_energy
                          << " keV (observed mean " << center << ", predicted "
                          << predicted << " counts, observed=" << observed << ")" << std::endl;
            }
          }
          recovered_anchor_failure = !predicted_preserved || !fitted_preserved;
        }
        const bool refinement_observable_failure = !recovered_source_anchors.empty()
            && !observable_requested_anchors_preserved( solution, refined_solution );
        bool anchor_failure = recovered_anchor_failure || refinement_observable_failure
          || (rescue_solve_this_iteration
            ? requested_anchors_catastrophically_removed(solution, refined_solution)
            : (!auto_interferer_lines.empty()
                && !requested_anchors_preserved(solution, refined_solution, anchor_exclusions)));

        // Apply a locally selected structural change as its own transaction.  The exact parent and
        // two child ROIs are carried from the measured-data decision, in its calibration frame;
        // unrelated re-clustering geometry is never imported.  Translating by channel into the
        // original spectrum frame makes the local solve use the identical bins that won AICc.
        bool retained_component_transaction = false;
        bool attempted_component_transaction = false;
        const auto label_component_proposal
          = [&]( const double affected_lower, const double affected_upper,
                 const bool accepted, const std::string &reason ) {
            const std::string stage_prefix = refinement_stage;
            const double midpoint = 0.5*(affected_lower + affected_upper);
            const double local_fwhm = solution_fwhm_at_energy( midpoint );
            const double tolerance = std::max( 1.0,
                (std::isfinite(local_fwhm) && (local_fwhm > 0.0)) ? local_fwhm : 0.0 );
            for( AutomaticRoiDecisionDiagnostic &diagnostic
                 : result.automatic_roi_diagnostics )
            {
              if( diagnostic.stage.compare(0, stage_prefix.size(), stage_prefix) != 0 )
                continue;
              double diagnostic_lower = diagnostic.left_lower;
              double diagnostic_upper = diagnostic.left_upper;
              if( diagnostic.right_upper > diagnostic.right_lower )
              {
                diagnostic_lower = std::min( diagnostic_lower, diagnostic.right_lower );
                diagnostic_upper = std::max( diagnostic_upper, diagnostic.right_upper );
              }
              const bool overlaps = (diagnostic_lower < (affected_upper + tolerance))
                  && (diagnostic_upper > (affected_lower - tolerance));
              if( !overlaps )
                continue;
              const std::string status = accepted
                  ? "accepted component transaction"
                  : "rolled-back component proposal";
              if( diagnostic.stage.find(status) == std::string::npos )
                diagnostic.stage += " [" + status + "]";
              const std::string prior_decision
                = automatic_roi_decision_name( diagnostic.decision );
              const std::string prior_reason = diagnostic.reason;
              diagnostic.reason = reason + "; proposed " + prior_decision
                  + (prior_reason.empty() ? std::string() : ": " + prior_reason);
              if( !accepted )
                diagnostic.decision = AutomaticRoiDecision::MergeInseparableWide;
            }
            if( should_debug_print() )
            {
              std::cerr << "Component-local diagnostic status [" << affected_lower
                        << ", " << affected_upper << "] keV: "
                        << (accepted ? "accepted" : "rolled back")
                        << "; " << reason << std::endl;
            }
          };
        if( use_automatic_roi_policy && auto_interferer_lines.empty()
            && !rescue_solve_this_iteration && !selected_component_partitions.empty() )
        {
          attempted_component_transaction = true;
          try
          {
            const std::shared_ptr<const SpecUtils::EnergyCalibration> original_cal
              = orig_foreground->energy_calibration();
            if( !original_cal || !original_cal->valid() )
              throw std::runtime_error( "original calibration is unavailable" );
            const auto to_original_energy = [&original_cal](
                const std::shared_ptr<const SpecUtils::EnergyCalibration> &from,
                const double energy ) {
              if( !from || !from->valid() )
                throw std::runtime_error( "partition calibration is unavailable" );
              return original_cal->energy_for_channel( from->channel_for_energy( energy ) );
            };

            std::vector<SelectedRoiComponentPartition> translated_partitions;
            translated_partitions.reserve( selected_component_partitions.size() );
            for( const SelectedRoiComponentPartition &record : selected_component_partitions )
            {
              SelectedRoiComponentPartition translated = record;
              translated.parent_lower = to_original_energy(
                  record.calibration, record.parent_lower );
              translated.parent_upper = to_original_energy(
                  record.calibration, record.parent_upper );
              translated.calibration = original_cal;
              for( RelActCalcAuto::RoiRange &child : translated.children )
              {
                child.lower_energy = to_original_energy(
                    record.calibration, child.lower_energy );
                child.upper_energy = to_original_energy(
                    record.calibration, child.upper_energy );
              }
              if( (translated.children.size() < 2)
                  || !(translated.parent_upper > translated.parent_lower)
                  || (std::adjacent_find( std::begin(translated.children),
                      std::end(translated.children),
                      []( const RelActCalcAuto::RoiRange &lhs,
                          const RelActCalcAuto::RoiRange &rhs ) {
                        return rhs.lower_energy < lhs.upper_energy;
                      } ) != std::end(translated.children)) )
                throw std::runtime_error( "translated partition is not a disjoint child transaction" );
              translated_partitions.push_back( std::move(translated) );
            }

            RelActCalcAuto::RelActAutoSolution working_solution = solution;
            size_t retained_partitions = 0;
            const auto mark_partition_rolled_back
              = [&]( const double affected_lower, const double affected_upper,
                     const std::string &reason ) {
                label_component_proposal(
                    affected_lower, affected_upper, false, reason );
              };
            for( size_t partition_index = 0;
                 partition_index < translated_partitions.size(); ++partition_index )
            {
              const SelectedRoiComponentPartition &record
                = translated_partitions[partition_index];
              // The proposal's child ranges are the geometry being installed.  After an
              // energy-calibration refinement their translated outer edges can differ slightly
              // from the bookkeeping parent range captured before translation.  Matching only
              // the stale parent can leave an incumbent ROI overlapping a real child and falsely
              // reject an otherwise valid exact component transaction.
              double replacement_lower = std::numeric_limits<double>::max();
              double replacement_upper = -std::numeric_limits<double>::max();
              for( const RelActCalcAuto::RoiRange &child : record.children )
              {
                replacement_lower = std::min( replacement_lower, child.lower_energy );
                replacement_upper = std::max( replacement_upper, child.upper_energy );
              }
              if( !(replacement_upper > replacement_lower) )
                throw std::runtime_error( "partition children have invalid replacement extent" );
              RelActCalcAuto::Options local_options = working_solution.m_options;
              std::vector<RelActCalcAuto::RoiRange> local_rois;
              bool matched_partition = false;
              double affected_lower = std::numeric_limits<double>::max();
              double affected_upper = -std::numeric_limits<double>::max();
              for( const RelActCalcAuto::RoiRange &incumbent : local_options.rois )
              {
                const bool overlaps = (incumbent.lower_energy < replacement_upper)
                    && (incumbent.upper_energy > replacement_lower);
                if( overlaps )
                {
                  matched_partition = true;
                  affected_lower = std::min( affected_lower, incumbent.lower_energy );
                  affected_upper = std::max( affected_upper, incumbent.upper_energy );
                }else
                {
                  local_rois.push_back( incumbent );
                }
              }
              if( !matched_partition )
              {
                if( should_debug_print() )
                  std::cerr << "Component-local partition did not match an incumbent component"
                            << std::endl;
                mark_partition_rolled_back( record.parent_lower, record.parent_upper,
                    "whole-component partition rolled back: no matching incumbent component" );
                continue;
              }
              local_rois.insert( std::end(local_rois), std::begin(record.children),
                  std::end(record.children) );
              std::sort( std::begin(local_rois), std::end(local_rois),
                []( const RelActCalcAuto::RoiRange &lhs,
                    const RelActCalcAuto::RoiRange &rhs ) {
                  return lhs.lower_energy < rhs.lower_energy;
                } );
              const bool disjoint = std::adjacent_find( std::begin(local_rois),
                  std::end(local_rois),
                  []( const RelActCalcAuto::RoiRange &lhs,
                      const RelActCalcAuto::RoiRange &rhs ) {
                    return rhs.lower_energy < lhs.upper_energy;
                  } ) == std::end(local_rois);
              if( !disjoint )
              {
                if( should_debug_print() )
                  std::cerr << "Component-local partition overlapped an unchanged incumbent ROI"
                            << std::endl;
                mark_partition_rolled_back( affected_lower, affected_upper,
                    "whole-component partition rolled back: child overlapped incumbent geometry" );
                continue;
              }

              local_options.rois = std::move( local_rois );
              remove_floating_peaks_without_roi( local_options );
              RelActCalcAuto::RelActAutoSolution local_solution = RelActCalcAuto::solve(
                  local_options, orig_foreground, orig_background,
                  drf, auto_search_peaks, det_type );
              const bool solved = RelActCalcAuto::RelActAutoSolution::is_usable_status(
                                    local_solution.m_status);
              const bool local_predicted = solved && (recovered_source_anchors.empty()
                  || ((sources_rel_eff_index >= 0)
                    && std::all_of( std::begin(recovered_source_anchors),
                        std::end(recovered_source_anchors),
                        [&local_solution, sources_rel_eff_index, &orig_foreground](
                            const RelActCalcManual::GenericPeakInfo &anchor ) {
                          return predicted_source_anchor_counts( local_solution,
                              static_cast<size_t>(sources_rel_eff_index), anchor,
                              orig_foreground->live_time() ) > sm_keep_gate_min_est_counts;
                        } )));
              const bool local_observed = local_predicted
                  && (recovered_source_anchors.empty()
                    || significant_requested_source_anchors_preserved( local_solution, sources,
                        recovered_source_anchors, config.roi_significance_z,
                        solution_fwhm_at_energy ));
              const bool incumbent_observed = local_observed
                  && observable_requested_anchors_preserved(
                      working_solution, local_solution );
              // A split must never exchange a significant requested fitted line for a lower
              // local continuum score.  Observable-peak preservation is deliberately retained
              // as the public-output contract; this stricter fitted-anchor check additionally
              // protects dense, unresolved source multiplets before observability filtering.
              const bool incumbent_requested = incumbent_observed
                  && requested_anchors_preserved( working_solution, local_solution, {} );
              if( should_debug_print() )
              {
                std::cerr << "Component-local challenger [" << record.parent_lower
                          << ", " << record.parent_upper << "] keV: status="
                          << static_cast<int>(local_solution.m_status)
                          << ", predicted=" << local_predicted
                          << ", clean anchors observed=" << local_observed
                          << ", incumbent observable anchors=" << incumbent_observed
                          << ", incumbent fitted anchors=" << incumbent_requested << std::endl;
              }
              if( local_predicted && local_observed && incumbent_requested )
              {
                working_solution = std::move( local_solution );
                ++retained_partitions;
                label_component_proposal( affected_lower, affected_upper, true,
                    "whole-component partition retained by its component-local solve" );
              }else
              {
                mark_partition_rolled_back( affected_lower, affected_upper,
                    "whole-component partition rolled back by its component-local solve" );
              }
            }
            if( retained_partitions > 0 )
            {
              solution = std::move( working_solution );
              retained_component_transaction = true;
              result.warnings.push_back( "Retained "
                  + std::to_string(retained_partitions)
                  + " measured-data ROI component partition(s) while transactionally"
                    " retaining all unrelated incumbent geometry." );
            }
          }
          catch( const std::exception &error )
          {
            if( should_debug_print() )
              std::cerr << "Component-local partition challenger failed: "
                        << error.what() << std::endl;
          }
        }
        const auto mark_component_partitions_rolled_back
          = [&]( const std::string &reason ) {
            label_component_proposal(
                min_valid_energy, max_valid_energy, false, reason );
          };
        if( retained_component_transaction )
        {
          // A committed local transaction is the final result for this refinement pass.  A later
          // ordinary all-ROI re-cluster has no knowledge of the exact child transaction and can
          // silently re-merge it, leaving an "accepted" diagnostic while returning the old wide
          // public ROI.  All selected components for this pass have already been applied through
          // the same isolated transaction above; a subsequent geometry change must be proposed
          // by a fresh top-level fit, not by importing global re-clustering state here.
          break;
        }
        if( attempted_component_transaction )
        {
          // A selected component is evaluated only through the exact local transaction above.
          // Falling through would let the enclosing all-ROI refinement import that same rejected
          // geometry (and unrelated re-clustering changes), defeating transactional rollback.
          mark_component_partitions_rolled_back(
              "whole-component partition rolled back because its local transaction was rejected" );
          result.warnings.push_back( "Rejected a measured-data ROI component partition because"
            " its component-local solve failed or removed an observable requested-source anchor;"
            " retained the prior accepted solution." );
          break;
        }
        if( anchor_failure )
        {
          mark_component_partitions_rolled_back(
              "whole-component partition rolled back with rejected refinement" );
          result.warnings.push_back( refinement_observable_failure
            ? "Rejected ROI refinement because it removed an observable requested-source"
              " anchor; retained the prior accepted solution."
            : (recovered_anchor_failure
              ? "Rejected ROI refinement because it removed source evidence recovered by the"
                " source-clean challenger; retained the prior accepted solution."
            : (rescue_solve_this_iteration
              ? "Rejected the bounded marginal-line rescue because it removed a significant"
                " requested-source anchor; retained the incumbent source fit."
              : "Rejected a post-interferer ROI refinement because it removed a significant"
                " requested-source anchor; retained the prior accepted solution.") ) );
          break;
        }

        // Developer check: Validate refined solution's final ROIs don't overlap
#if( PERFORM_DEVELOPER_CHECKS )
        for( size_t i = 1; i < refined_solution.m_final_roi_ranges.size(); ++i )
        {
          const RelActCalcAuto::RoiRange &prev_roi = refined_solution.m_final_roi_ranges[i - 1];
          const RelActCalcAuto::RoiRange &curr_roi = refined_solution.m_final_roi_ranges[i];
          if( curr_roi.lower_energy < prev_roi.upper_energy )
          {
            std::cerr << "ERROR: RelActAuto returned overlapping ROIs[" << (i-1) << "] and [" << i << "]: "
                 << "[" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] vs "
                 << "[" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "]" << std::endl;
            assert( curr_roi.lower_energy >= prev_roi.upper_energy );
          }
        }
#endif

#if( PERFORM_DEVELOPER_CHECKS )
        if( rescue_solve_this_iteration
            && sm_force_next_rescue_evaluation_failure_for_test )
        {
          sm_force_next_rescue_evaluation_failure_for_test = false;
          throw std::runtime_error( "forced bounded-rescue post-solve evaluation failure" );
        }
#endif

        // Compute filtered chi2-per-channel that only includes ROIs with significant peaks.
        // This avoids the problem where adding a ROI in a flat region (with no real peaks)
        // would artificially reduce the average chi2.  Each solution is evaluated against its
        // own m_foreground (see compute_filtered_chi2_per_channel), so the incumbent and the
        // challenger are each scored in their own consistent calibration frame - the loop-local
        // `foreground` may be one cal-step ahead of the incumbent solution.
        std::vector<size_t> old_insignificant_rois, new_insignificant_rois;
        const double old_chi2_dof = compute_filtered_chi2_per_channel(
          solution, config.roi_significance_z, old_insignificant_rois );
        const double new_chi2_dof = compute_filtered_chi2_per_channel(
          refined_solution, config.roi_significance_z, new_insignificant_rois );

        // Check if chi2/channel improved
        if( !rescue_solve_this_iteration && (new_chi2_dof >= old_chi2_dof) )
        {
          mark_component_partitions_rolled_back(
              "whole-component partition rolled back because the enclosing solve did not improve" );
          if( should_debug_print() )
          {
            std::cout << "Iteration " << iter << " did not improve filtered chi2/channel ("
                 << old_chi2_dof << " -> " << new_chi2_dof << "), stopping" << std::endl;
            if( !new_insignificant_rois.empty() )
              std::cout << "  (" << new_insignificant_rois.size() << " ROIs had insignificant peaks)" << std::endl;
          }
          break;
        }

        solution = std::move( refined_solution );

        if( rescue_solve_this_iteration )
        {
          for( const RelActCalcAuto::RoiRange &rescued : proposed_rescued_rois )
          {
            const double center = 0.5 * (rescued.lower_energy + rescued.upper_energy);
            const bool retained = std::any_of( std::begin(solution.m_options.rois),
              std::end(solution.m_options.rois), [center]( const RelActCalcAuto::RoiRange &roi ) {
                return (center >= roi.lower_energy) && (center <= roi.upper_energy);
              } );
            if( retained )
              rescued_roi_ranges.push_back( rescued );
          }
          result.warnings.push_back( "Admitted " + std::to_string(rescued_roi_ranges.size())
            + " marginal source ROI(s) through the bounded fit-then-prune rescue; final ROI"
              " significance filtering remains authoritative." );
        }

        if( should_debug_print() )
          std::cout << "Iteration " << iter << " improved: chi2/dof=" << new_chi2_dof
               << " (was " << old_chi2_dof << ")" << std::endl;
        }
        catch( const std::exception &error )
        {
          if( !rescue_transaction_requested )
            throw;
          result.warnings.push_back( "The bounded marginal-line rescue challenger threw during"
            " preparation or evaluation; retained the successful incumbent source fit ("
            + std::string(error.what()) + ")." );
          break;
        }
        catch( ... )
        {
          if( !rescue_transaction_requested )
            throw;
          result.warnings.push_back( "The bounded marginal-line rescue challenger threw during"
            " preparation or evaluation; retained the successful incumbent source fit." );
          break;
        }
      }//for( size_t iter = 0; iter < max_iterations; ++iter )

      if( should_debug_print() )
      {
        std::cout << "Final solution after refinement:" << std::endl;
        solution.print_summary( std::cout );
        std::cout << std::endl;

        // Print ROIs and sum fit peak areas for each ROI
        std::cout << "Solution ROIs and fit peak areas:" << std::endl;
      }

      for( size_t roi_index = 0; roi_index < solution.m_final_roi_ranges.size(); ++roi_index )
      {
        const RelActCalcAuto::RoiRange &roi = solution.m_final_roi_ranges[roi_index];
        double sum_peak_area = 0.0;
        size_t num_peaks_in_roi = 0;

        for( const PeakDef &peak : solution.m_fit_peaks )
        {
          const double peak_roi_lower = peak.continuum()->lowerEnergy();
          const double peak_roi_upper = peak.continuum()->upperEnergy();
          // Match if peak's ROI bounds are within 1 keV of the solution ROI bounds
          if( (std::fabs( peak_roi_lower - roi.lower_energy ) < 1.0)
             && (std::fabs( peak_roi_upper - roi.upper_energy ) < 1.0) )
          {
            sum_peak_area += peak.peakArea();
            ++num_peaks_in_roi;
          }
        }//for( loop over fit peaks )

        if( should_debug_print() )
          std::cout << "  ROI " << roi_index << ": [" << roi.lower_energy << ", " << roi.upper_energy << "] keV"
               << ", " << num_peaks_in_roi << " peaks, sum area = " << sum_peak_area << std::endl;
      }//for( loop over ROIs )

      if( should_debug_print() )
        std::cout << std::endl;
    }//iterative refinement

    // Identify ROIs without significant peaks for filtering.  The significance is evaluated
    // against solution.m_foreground (see compute_filtered_chi2_per_channel), so the result is
    // calibration-consistent regardless of how the refinement loop exited (the loop-local
    // `foreground` can be one cal-step past the accepted solution).
    std::vector<size_t> final_insignificant_rois;
    compute_filtered_chi2_per_channel( solution,
      config.roi_significance_z, final_insignificant_rois );

    // Build set of insignificant ROI ranges for filtering
    std::vector<std::pair<double,double>> insignificant_roi_ranges;
    for( const size_t roi_idx : final_insignificant_rois )
    {
      // Use the spectrum-cal ROI ranges (not the true-energy m_final_roi_ranges): these bounds are
      // compared below against peak.continuum() bounds, which are in spectrum cal.  When energy cal is
      // fit (NaI/CZT) the two cals can differ by >1 keV, so true-energy ranges would mis-match.
      // (Fall back to the true-energy ranges if the spectrum-cal vector was not populated,
      // matching the index source used inside compute_filtered_chi2_per_channel.)
      const RelActCalcAuto::RoiRange &roi
        = (roi_idx < solution.m_final_roi_ranges_in_spectrum_cal.size())
          ? solution.m_final_roi_ranges_in_spectrum_cal[roi_idx]
          : solution.m_final_roi_ranges[roi_idx];
      insignificant_roi_ranges.emplace_back( roi.lower_energy, roi.upper_energy );
    }

    // Populate result, filtering out peaks from insignificant ROIs
    result.status = solution.m_status;
    result.error_message = solution.m_error_message;

    // Returns true if peak_src is one of the input sources being fit.
    auto is_input_source = [&sources]( const PeakDef &p ) -> bool {
      RelActCalcAuto::SrcVariant psrc;
      if( p.parentNuclide() )
        psrc = p.parentNuclide();
      if( p.xrayElement() )
        psrc = p.xrayElement();
      if( p.reaction() )
        psrc = p.reaction();
      
      if( RelActCalcAuto::is_null( psrc ) )
        return false;
      for( const RelActCalcAuto::NucInputInfo &src : sources )
      {
        if( src.source == psrc )
          return true;
      }
      return false;
    };//is_input_source lambda

    const auto is_auto_interferer = [&]( const PeakDef &peak ) -> bool {
      const SandiaDecay::Nuclide * const parent = peak.parentNuclide();
      if( parent && auto_interferer_nucs.count(parent) )
        return true;
      if( peak.hasSourceGammaAssigned() )
        return false;
      for( const double energy : auto_interferer_float_energies )
      {
        if( std::fabs(peak.mean() - energy) < 1.0 )
          return true;
      }
      return false;
    };//is_auto_interferer lambda

    /** True if `peak` is just the model's stand-in for an existing user peak we carried into the
     fit as a bystander FloatingPeak.  Such a peak has no source assigned (the fit does not know
     what the bystander was), and the bystander is reported separately - with its source, color and
     label intact - by the bystander blocks further below.  Reporting both would put two peaks a
     fraction of a FWHM apart into the user's peak list for one physical line.
     */
    const auto is_carried_bystander_float = [&]( const PeakDef &peak ) -> bool {
      if( peak.hasSourceGammaAssigned() )
        return false;
      // Wide enough to cover the fitted mean drifting off the enrolled energy, but still well
      //  inside one peak width, at any detector resolution.
      const double tol = std::max( 1.0, 0.5*peak.fwhm() );
      for( const double energy : carried_bystander_float_energies )
      {
        if( std::fabs(peak.mean() - energy) < tol )
          return true;
      }
      return false;
    };//is_carried_bystander_float lambda

    const auto is_rescued_source_peak = [&]( const PeakDef &peak ) -> bool {
      if( !is_input_source(peak) )
        return false;
      const double energy = peak.hasSourceGammaAssigned()
        ? peak.gammaParticleEnergy() : peak.mean();
      return std::any_of( std::begin(rescued_roi_ranges), std::end(rescued_roi_ranges),
        [energy]( const RelActCalcAuto::RoiRange &roi ) {
          return (energy >= roi.lower_energy) && (energy <= roi.upper_energy);
        } );
    };//is_rescued_source_peak lambda

    const auto hide_from_public_results = [&]( const PeakDef &peak ) -> bool {
      return is_auto_interferer( peak )
          || (norm_peaks_dont_use && !is_input_source(peak));
    };//hide_from_public_results lambda

    // Keep accepted automatic R6 nuisance peaks in this private model through combining and the
    // observable LM refit.  Other FitNormBkgrndPeaksDontUse peaks retain their legacy early-filter
    // behavior; broadening filter-late beyond the transactional R6 path would alter unrelated modes.
    std::vector<PeakDef> full_model_peaks;
    result.fit_peaks.clear();
    if( should_debug_print() && !insignificant_roi_ranges.empty() )
      std::cout << "Peak filtering by ROI significance:" << std::endl;

    for( const PeakDef &peak : solution.m_peaks_without_back_sub )
    {
      const double mean = peak.mean();

      const bool hidden = hide_from_public_results( peak );
      const bool keep_in_private_model = !hidden || is_auto_interferer(peak);
      if( hidden && norm_peaks_dont_use && !is_input_source(peak) )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        // Sanity: filtered peak should either have no source assigned (free peak)
        // or its parentNuclide should be one of the known NORM nuclides.
        if( peak.hasSourceGammaAssigned() )
        {
          const SandiaDecay::Nuclide * const pnuc = peak.parentNuclide();
          const bool is_known_norm = pnuc
            && ( (pnuc->symbol == "U238") || (pnuc->symbol == "Ra226")
              || (pnuc->symbol == "U235") || (pnuc->symbol == "Th232")
              || (pnuc->symbol == "K40") );
          if( !is_known_norm )
            std::cerr << "WARNING: FitNormBkgrndPeaksDontUse filtered peak at " << mean
                      << " keV with unexpected source: " << peak.sourceName() << std::endl;
        }
#endif
        if( should_debug_print() )
          std::cout << "  Filtered before observable refit (NORM/background): " << mean << " keV" << std::endl;
      }
      else if( hidden && should_debug_print() )
      {
        std::cout << "  Hidden until observable refit (auto interferer): " << mean << " keV" << std::endl;
      }

      const double peak_roi_lower = peak.continuum()->lowerEnergy();
      const double peak_roi_upper = peak.continuum()->upperEnergy();
      bool in_insignificant_roi = false;

      for( const std::pair<double,double> &roi_range : insignificant_roi_ranges )
      {
        // Match if peak's ROI bounds are within 1 keV of the insignificant ROI bounds
        if( (std::fabs( peak_roi_lower - roi_range.first ) < 1.0)
           && (std::fabs( peak_roi_upper - roi_range.second ) < 1.0) )
        {
          in_insignificant_roi = true;
          break;
        }
      }

      if( in_insignificant_roi )
      {
        if( should_debug_print() )
        {
          std::cout << "  Filtered (insignificant ROI [" << peak.continuum()->lowerEnergy() << ", " << peak.continuum()->upperEnergy()
               << "] keV): peak at " << mean << " keV, area = " << peak.peakArea() << std::endl;
        }
      }else
      {
        bool mean_in_roi = false;
        for( size_t roi_index = 0; !mean_in_roi && (roi_index < solution.m_final_roi_ranges_in_spectrum_cal.size()); ++roi_index )
        {
          const auto pos = std::find( std::begin(final_insignificant_rois), std::end(final_insignificant_rois), roi_index );
          if( pos != std::end(final_insignificant_rois) )
            continue;
          const RelActCalcAuto::RoiRange &roi = solution.m_final_roi_ranges_in_spectrum_cal[roi_index];
          mean_in_roi = ((mean >= roi.lower_energy) && (mean <= roi.upper_energy));
        }

        // Also require the peak mean is within its own continuum/ROI range;
        //  peaks whose mean is outside their continuum are just tail contributions
        //  to adjacent ROIs and should not be included as separate peaks.
        if( mean_in_roi )
          mean_in_roi = (mean >= peak_roi_lower) && (mean <= peak_roi_upper);

        if( should_debug_print() && !insignificant_roi_ranges.empty() )
        {
          std::cout << "  Kept (significant ROI [" << peak.continuum()->lowerEnergy() << ", " << peak.continuum()->upperEnergy()
               << "] keV): peak at " << mean << " keV, area = " << peak.peakArea() << (mean_in_roi ? " (was in ROI)" : " (skipping peak, not in a ROI)") << std::endl;
        }
        
        if( mean_in_roi )
        {
          if( keep_in_private_model )
            full_model_peaks.push_back( peak );
          if( !hidden )
            result.fit_peaks.push_back( peak );
        }
      }
    }//for( const PeakDef &peak : solution.m_peaks_without_back_sub )

    if( !insignificant_roi_ranges.empty() )
    {
      const size_t num_filtered = solution.m_peaks_without_back_sub.size() - result.fit_peaks.size();
      if( should_debug_print() )
        std::cout << "Filtered out " << num_filtered << " peaks from "
             << insignificant_roi_ranges.size() << " ROIs without significant chi2 improvement" << std::endl;
    }

    // With FitNormBkgrndPeaksDontUse and no requested sources, every fit peak is a NORM peak and
    // gets filtered, leaving a Success status with zero peaks - warn so callers can tell this from
    // "no peaks were observable".
    if( norm_peaks_dont_use && sources.empty() && result.fit_peaks.empty()
        && !solution.m_peaks_without_back_sub.empty() )
    {
      result.warnings.push_back( "All fit peaks were NORM background peaks, which"
        " FitNormBkgrndPeaksDontUse excludes from the returned results." );
    }

    

    result.solution = std::move( solution );

    // The peaks in result.fit_peaks come from the final accepted solution and are expressed in
    // `result.solution.m_foreground`'s energy calibration.  On the common loop-exit paths the
    // loop-local `foreground` has been advanced one cal-step beyond the accepted solution (by that
    // solution's own fitted energy-cal adjustment), so using it here would mis-translate the peaks and
    // refit observable peaks against mismatched data.  Use the solution's own foreground as the
    // canonical "fitted" spectrum/cal for both observable-peak refitting and the translate-back below.
    const shared_ptr<const SpecUtils::Measurement> solution_foreground
      = result.solution.m_foreground ? result.solution.m_foreground : foreground;

    std::vector<PeakDef> full_uncombined_peaks = std::move( full_model_peaks );
    const auto may_combine_model_peaks = [&]( const PeakDef &lhs, const PeakDef &rhs ) -> bool {
      // Never let a hidden nuisance/source peak donate its area and provenance to the other class.
      // Combining within the public class or within the nuisance class remains unchanged.
      return is_auto_interferer(lhs) == is_auto_interferer(rhs);
    };
    const auto must_refit_model_peak = [&]( const PeakDef &peak ) -> bool {
      // Nuisances must remain in the private model to prevent re-absorption.  A rescued source peak
      // already passed the ROI Wilks gate, so let it reach the honest LM refit instead of allowing
      // the coarse S/sqrt(S+B) prefilter to undo R2; unlike nuisances, it must still pass the final
      // post-refit source significance threshold.
      return is_auto_interferer(peak) || is_rescued_source_peak(peak);
    };
    std::vector<PeakDef> full_combined_peaks
      = combine_overlapping_peaks_in_rois( full_uncombined_peaks, may_combine_model_peaks );

    if( should_debug_print() && (full_combined_peaks.size() != full_uncombined_peaks.size()) )
    {
      std::cout << "Combined " << full_uncombined_peaks.size() << " full-model peaks into "
           << full_combined_peaks.size() << " peaks" << std::endl;
    }

    // Translate peaks from fitted energy cal back to original energy cal (if needed)
    auto translate_peaks_to_orig_cal = [&]( vector<PeakDef> &peaks )
    {
      if( !apply_energy_cal_between || !config.fit_energy_cal )
        return;

      const shared_ptr<const SpecUtils::EnergyCalibration> fitted_cal = solution_foreground->energy_calibration();
      const shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = orig_foreground->energy_calibration();

      if( !fitted_cal || !orig_cal || (*fitted_cal == *orig_cal) )
        return;

      deque<shared_ptr<const PeakDef>> tmp_peaks;
      for( const PeakDef &p : peaks )
        tmp_peaks.push_back( make_shared<const PeakDef>( p ) );

      const deque<shared_ptr<const PeakDef>> translated
        = EnergyCal::translatePeaksForCalibrationChange( tmp_peaks, fitted_cal, orig_cal );

      peaks.clear();
      for( const shared_ptr<const PeakDef> &p : translated )
        peaks.push_back( *p );
    };

#if( OBSERVABLE_PEAKS_USING_ORIGINAL_CAL_WITH_BACK_SUB )
    // Translate the complete model to original energy cal first, then compute observable peaks on
    // the original foreground with background
    // subtraction.  This avoids the poor continuum fits that result from refitting
    // on the energy-cal-adjusted spectrum and then translating peaks back.
    translate_peaks_to_orig_cal( full_combined_peaks );
    translate_peaks_to_orig_cal( full_uncombined_peaks );

    std::vector<PeakDef> full_observable_peaks = compute_observable_peaks(
      full_combined_peaks, orig_foreground, det_type, config, orig_background,
      may_combine_model_peaks, must_refit_model_peak, is_auto_interferer );
#else
    // Existing path: refit the complete model on fitted-cal foreground, then translate all peaks.
    std::vector<PeakDef> full_observable_peaks
      = compute_observable_peaks( full_combined_peaks, solution_foreground, det_type, config,
          may_combine_model_peaks, must_refit_model_peak, is_auto_interferer );

    translate_peaks_to_orig_cal( full_combined_peaks );
    translate_peaks_to_orig_cal( full_uncombined_peaks );
    translate_peaks_to_orig_cal( full_observable_peaks );
#endif

    const auto public_peaks_only = [&]( const std::vector<PeakDef> &peaks ) {
      std::vector<PeakDef> public_peaks;
      public_peaks.reserve( peaks.size() );
      for( const PeakDef &peak : peaks )
      {
        if( !hide_from_public_results(peak) )
          public_peaks.push_back( peak );
      }
      return public_peaks;
    };

    result.uncombined_fit_peaks = public_peaks_only( full_uncombined_peaks );
    result.fit_peaks = public_peaks_only( full_combined_peaks );
    result.observable_peaks = public_peaks_only( full_observable_peaks );

    // Drop the model stand-ins for carried bystanders; the bystander blocks below append the real
    //  thing (or deliberately leave the user's original peak in place).
    result.observable_peaks.erase(
      std::remove_if( std::begin(result.observable_peaks), std::end(result.observable_peaks),
                      is_carried_bystander_float ),
      std::end(result.observable_peaks) );

    // Sort observable_peaks by mean energy for deterministic ordering.
    // compute_observable_peaks processes ROIs via a map keyed by continuum
    // pointer, whose iteration order is non-deterministic across runs.
    std::sort( result.observable_peaks.begin(), result.observable_peaks.end(),
      &PeakDef::lessThanByMean );

    if( should_debug_print() && (result.observable_peaks.size() != result.fit_peaks.size()) )
    {
      std::cout << "Observable peaks: " << result.observable_peaks.size() << " of "
           << result.fit_peaks.size() << " fit_peaks" << std::endl;
    }

    // Developer check: Look for duplicate peaks (same mean, different ROI)
#if( PERFORM_DEVELOPER_CHECKS )
    for( size_t i = 0; i < result.observable_peaks.size(); ++i )
    {
      for( size_t j = i + 1; j < result.observable_peaks.size(); ++j )
      {
        const PeakDef &peak_i = result.observable_peaks[i];
        const PeakDef &peak_j = result.observable_peaks[j];
        const double mean_diff = std::fabs( peak_i.mean() - peak_j.mean() );
        if( mean_diff < 0.5 )  // Same energy within 0.5 keV
        {
          const bool same_continuum = (peak_i.continuum() == peak_j.continuum());
          if( !same_continuum )
          {
            std::cerr << "WARNING: Duplicate peaks at " << peak_i.mean() << " keV with different ROIs: "
                 << "[" << peak_i.continuum()->lowerEnergy() << ", " << peak_i.continuum()->upperEnergy() << "] vs "
                 << "[" << peak_j.continuum()->lowerEnergy() << ", " << peak_j.continuum()->upperEnergy() << "]" << std::endl;
          }
        }
      }
    }
#endif

    // Populate original_peaks_to_remove for ExistingPeaksAsFreePeak.
    //
    // Source-matched peaks (assigned to one of our fit sources) are unconditionally removed
    // since the fit produced a replacement.
    //
    // Bystander peaks (assigned to a different source) were added as FloatingPeaks so the
    // RelActAuto fit could account for their signal.  We reconstruct updated bystander PeakDefs
    // from the FloatingPeakResult (which has the fit-updated amplitude/fwhm), preserving the
    // original source attribution and color.  The updated bystander is added to observable_peaks
    // and the original is marked for removal.
    //
    // A bystander is only updated/removed if:
    //   1. Its energy falls within an observable peak's ROI (the ROI was used in the solution), AND
    //   2. The fit actually produced source-attributed peaks in that ROI (so the removal is justified), AND
    //   3. A matching FloatingPeakResult exists with positive amplitude.
    // If any condition fails, the bystander is left untouched.
    // Helper: returns true if any peak in result.fit_peaks that is attributed to one of
    // the requested fit sources (or NORM nuclides when fitting NORM) has its mean within
    // [lower_energy, upper_energy].  This tells us whether the fit produced useful results
    // for a given ROI, justifying the replacement of bystander peaks.
    const auto fit_has_source_peak_in_range = [&]( const double lower_energy,
                                                    const double upper_energy ) -> bool
    {
      for( const PeakDef &fp : result.fit_peaks )
      {
        if( (fp.mean() < lower_energy) || (fp.mean() > upper_energy) )
          continue;

        const SandiaDecay::Nuclide *fp_nuc = fp.parentNuclide();
        const SandiaDecay::Element *fp_el = fp.xrayElement();
        const ReactionGamma::Reaction *fp_rxn = fp.reaction();

        if( !fp_nuc && !fp_el && !fp_rxn )
          continue;

        for( const RelActCalcAuto::NucInputInfo &src : sources )
        {
          if( fp_nuc && (RelActCalcAuto::nuclide( src.source ) == fp_nuc) )
            return true;
          if( fp_el && (RelActCalcAuto::element( src.source ) == fp_el) )
            return true;
          if( fp_rxn && (RelActCalcAuto::reaction( src.source ) == fp_rxn) )
            return true;
        }

        if( fit_norm_peaks && fp_nuc && db )
        {
          const SandiaDecay::Nuclide *norm_nucs[] = {
            db->nuclide("U238"), db->nuclide("Ra226"), db->nuclide("U235"),
            db->nuclide("Th232"), db->nuclide("K40")
          };
          for( const SandiaDecay::Nuclide *norm_nuc : norm_nucs )
          {
            if( norm_nuc && (fp_nuc == norm_nuc) )
              return true;
          }
        }
      }//for( result.fit_peaks )

      return false;
    };//fit_has_source_peak_in_range

    // Helper: finds the FloatingPeakResult in result.solution.m_floating_peaks whose energy
    // matches the given bystander energy (in original calibration).  For ObservedInSpectrum
    // floating peaks, FloatingPeakResult::energy equals the input energy.
    // Each result is consumed at most once so two bystanders closer together than the match
    // tolerance cannot both bind to the same (nearest) result - each input floating peak
    // corresponds to exactly one bystander.
    std::set<const RelActCalcAuto::FloatingPeakResult *> consumed_floating_results;
    const auto find_floating_peak_result = [&]( const double bystander_energy )
      -> const RelActCalcAuto::FloatingPeakResult *
    {
      const RelActCalcAuto::FloatingPeakResult *best = nullptr;
      double best_diff = std::numeric_limits<double>::max();

      for( const RelActCalcAuto::FloatingPeakResult &fpr : result.solution.m_floating_peaks )
      {
        if( consumed_floating_results.count( &fpr ) )
          continue;

        const double diff = std::fabs( fpr.energy - bystander_energy );
        if( diff < best_diff )
        {
          best_diff = diff;
          best = &fpr;
        }
      }//for( result.solution.m_floating_peaks )

      // Use a tight tolerance since these should match very closely
      const double match_tol = 0.5; // keV
      if( best && (best_diff < match_tol) )
      {
        consumed_floating_results.insert( best );
        return best;
      }
      return nullptr;
    };//find_floating_peak_result

    if( existing_peaks_as_free
       && !existing_peaks_added_as_floating.empty()
       && RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
    {
      // Where the solver ended up putting each bystander that we manage to match to a
      //  FloatingPeakResult below, keyed by the original peak.  `find_floating_peak_result` consumes
      //  each result at most once, so the de-duplication pass afterwards must re-use what the match
      //  here found rather than ask for it a second time.
      std::map<const PeakDef *, double> bystander_fit_energies;

      // Everything the solver produced sits in [0, num_solver_observables); the loop below appends
      //  our own updated bystanders after it.  The de-duplication pass needs that split: it removes
      //  solver stand-ins by matching energy among *unattributed* peaks, and an updated bystander
      //  built from a user peak with no source assigned is unattributed and sits at exactly the
      //  stand-in's energy - so attribution cannot tell them apart, and matching by energy alone
      //  would delete the replacement we just made.
      const size_t num_solver_observables = result.observable_peaks.size();

      for( const std::pair<std::shared_ptr<const PeakDef>, double> &orig_and_energy : existing_peaks_added_as_floating )
      {
        const shared_ptr<const PeakDef> &orig_peak = orig_and_energy.first;
        const double orig_energy = orig_and_energy.second;

        if( peak_source_is_in_fit( orig_peak ) )
        {
          // Source-matched peak: always remove (fit produced replacement)
          result.original_peaks_to_remove.push_back( orig_peak );
        }else
        {
          // Bystander floating peak: only update/remove if justified.

          // Step 1: Find the observable ROI containing the bystander energy.
          std::shared_ptr<PeakContinuum> roi_continuum;
          double roi_lower = -1.0, roi_upper = -1.0;
          for( PeakDef &obs_peak : result.observable_peaks )
          {
            if( obs_peak.continuum()
               && (orig_energy >= obs_peak.continuum()->lowerEnergy())
               && (orig_energy <= obs_peak.continuum()->upperEnergy()) )
            {
              roi_continuum = obs_peak.continuum();
              roi_lower = obs_peak.continuum()->lowerEnergy();
              roi_upper = obs_peak.continuum()->upperEnergy();
              break;
            }
          }

          if( !roi_continuum )
          {
#if( PERFORM_DEVELOPER_CHECKS )
            std::cerr << "ExistingPeaksAsFreePeak: bystander at " << orig_energy
                 << " keV - NOT in any observable ROI. Observable ROIs:";
            std::set<const PeakContinuum *> printed;
            for( const PeakDef &obs : result.observable_peaks )
            {
              if( obs.continuum() && printed.insert( obs.continuum().get() ).second )
                std::cerr << " [" << obs.continuum()->lowerEnergy() << ", " << obs.continuum()->upperEnergy() << "]";
            }
            std::cerr << std::endl;
#endif
            continue;  // Bystander not in any observable ROI — leave untouched
          }

          // Step 2: Check if the fit produced source-attributed peaks in this ROI.
          //  If the fit added no peaks for the requested sources here, removing the
          //  bystander would be purely destructive — leave it alone.
          if( !fit_has_source_peak_in_range( roi_lower, roi_upper ) )
          {
#if( PERFORM_DEVELOPER_CHECKS )
            if( should_debug_print() )
            {
              std::cout << "ExistingPeaksAsFreePeak: bystander at " << orig_energy
                   << " keV retained — no source peaks in ROI ["
                   << roi_lower << ", " << roi_upper << "] keV" << std::endl;
            }
#endif
            continue;  // No source peaks in this ROI — leave bystander untouched
          }

          // Step 3: Reconstruct updated bystander from FloatingPeakResult.
          //  The fit updated the bystander's amplitude/fwhm to account for the influence
          //  of the new source's gammas.  We preserve the original source attribution.
          const RelActCalcAuto::FloatingPeakResult *fpr = find_floating_peak_result( orig_energy );

          if( !fpr )
          {
#if( PERFORM_DEVELOPER_CHECKS )
            std::cerr << "WARNING: ExistingPeaksAsFreePeak: bystander at " << orig_energy
                 << " keV has no matching FloatingPeakResult — keeping original." << std::endl;
#endif
            continue;  // No matching fit result — leave untouched
          }

          // Recorded before the retain-original bail-outs below: a retained bystander stays in the
          //  peak list, so the de-duplication pass needs its fitted energy just as much as an
          //  updated one does - more so, since that is the case that would otherwise show twice.
          bystander_fit_energies[orig_peak.get()] = fpr->original_spectrum_cal_energy;

          if( fpr->amplitude <= 0.0 )
          {
            // The solver gave the floating peak zero or negative amplitude.  We can't build a
            //  zero/negative-amplitude replacement, and silently dropping the user's existing peak
            //  would degrade their data when fitting a new source.  So retain the original bystander
            //  untouched (it keeps its own ROI) - the same conservative choice we make above when the
            //  fit produced no source peaks in the bystander's ROI.
#if( PERFORM_DEVELOPER_CHECKS )
            std::cerr << "ExistingPeaksAsFreePeak: bystander at " << orig_energy
                 << " keV has non-positive amplitude (" << fpr->amplitude
                 << ") — retaining original." << std::endl;
#endif
            continue;
          }

          // The solver can also come back with an amplitude it could not actually determine: the new
          //  source's gammas sitting on top of the bystander leave the pair degenerate, and the fit
          //  reports an uncertainty as large as, or larger than, the amplitude itself.  That is not a
          //  re-measurement of the user's peak, it is the fit saying it cannot tell, and swapping it
          //  in trades a measured peak for a meaningless one.  A missing or non-finite uncertainty is
          //  the same situation with even less to go on - and worse, the update below would then pair
          //  the fit's new amplitude with the ORIGINAL peak's (small) uncertainty, reporting a
          //  confident peak that nothing ever measured.  Both retain the original, for the same
          //  reason as the non-positive-amplitude case above.
          //
          //  The test is deliberately scale-free rather than a significance threshold: this
          //  uncertainty carries the solve's sqrt(chi2/dof) covariance inflation (RelActCalcAuto's
          //  `m_cov_scale`), while the thresholds used to filter observable peaks are calibrated on
          //  un-inflated ROI refits.  Comparing across those two scales would reject sound updates on
          //  any fit with a large chi2/dof; "the uncertainty exceeds the amplitude" needs no
          //  calibration to mean what it says.
          if( !(fpr->amplitude_uncert > 0.0)   //`!(x > 0)` so a NaN uncertainty lands here too
             || (fpr->amplitude_uncert >= fpr->amplitude) )
          {
#if( PERFORM_DEVELOPER_CHECKS )
            if( should_debug_print() )
            {
              std::cerr << "ExistingPeaksAsFreePeak: bystander at " << orig_energy
                   << " keV was not determined by the fit (amplitude " << fpr->amplitude
                   << " +- " << fpr->amplitude_uncert << ") - retaining original." << std::endl;
            }
#endif
            continue;
          }//if( the fit did not determine this bystander's amplitude )

          // Create updated bystander: copy original to preserve source info, color, labels, etc.
          PeakDef updated_bystander( *orig_peak );

          // Update Gaussian parameters from the fit result.
          // FloatingPeakResult::original_spectrum_cal_energy is in the original spectrum
          // calibration, same as the observable_peaks after translate_peaks_to_orig_cal.
          updated_bystander.setMean( fpr->original_spectrum_cal_energy );
          const double sigma = fpr->fwhm / (2.0 * std::sqrt( 2.0 * std::log( 2.0 ) ));
          updated_bystander.setSigma( sigma );
          updated_bystander.setAmplitude( fpr->amplitude );
          if( fpr->amplitude_uncert > 0.0 )
            updated_bystander.setAmplitudeUncert( fpr->amplitude_uncert );
          if( fpr->fwhm_uncert > 0.0 )
          {
            const double sigma_uncert = fpr->fwhm_uncert / (2.0 * std::sqrt( 2.0 * std::log( 2.0 ) ));
            updated_bystander.setSigmaUncert( sigma_uncert );
          }

          // Share the continuum with the observable peaks in the same ROI so the
          // bystander is part of the same ROI as the source peaks.
          updated_bystander.setContinuum( roi_continuum );

          if( should_debug_print() )
          {
            std::cout << "ExistingPeaksAsFreePeak: updating bystander at " << orig_energy
                 << " keV -> " << fpr->original_spectrum_cal_energy << " keV"
                 << " (amp: " << orig_peak->amplitude() << " -> " << fpr->amplitude
                 << ", fwhm: " << orig_peak->fwhm() << " -> " << fpr->fwhm << ")"
                 << ", source: " << orig_peak->sourceName() << std::endl;
          }

          result.observable_peaks.push_back( updated_bystander );
          result.original_peaks_to_remove.push_back( orig_peak );
        }
      }// for( existing_peaks_added_as_floating )

      // Remove unattributed observable peaks that originated from bystander FloatingPeaks.
      //
      // The RelActCalcAuto solver includes floating peaks in m_peaks_without_back_sub,
      // which flow through to observable_peaks without source attribution.  These must be
      // removed to prevent duplication with the original user peaks, which are either
      // retained as-is (when the fit produced no source peaks in the bystander's ROI) or
      // replaced by source-attributed updated bystanders (created above from FloatingPeakResult).
      {
        std::vector<double> bystander_fp_energies;
        for( const auto &orig_and_energy : existing_peaks_added_as_floating )
        {
          if( peak_source_is_in_fit( orig_and_energy.first ) )
            continue;

          bystander_fp_energies.push_back( orig_and_energy.second );

          // The observable being matched here is where the SOLVER put this floating peak, which the
          //  energy-cal fit can move well away from where it went in - so also match on the fitted
          //  energy.  Matching only the pre-fit energy leaves the tolerance below load-bearing for
          //  how far the calibration moved: a bystander retained above stays in the peak list, so a
          //  missed match here shows the same line twice, once from each.
          //
          //  Use the energy recorded when this bystander was matched above rather than looking it up
          //  again - `find_floating_peak_result` consumes each result at most once, so a repeat call
          //  for an already-matched bystander finds nothing, or worse binds to a neighbour's result.
          const std::map<const PeakDef *, double>::const_iterator fit_energy_pos
                                   = bystander_fit_energies.find( orig_and_energy.first.get() );
          if( fit_energy_pos != bystander_fit_energies.end() )
          {
            bystander_fp_energies.push_back( fit_energy_pos->second );
          }else
          {
            // Never matched above (e.g. the bystander sat in no observable ROI), so its result is
            //  still unconsumed, and here is the only chance to pick the fitted energy up.
            const RelActCalcAuto::FloatingPeakResult *fpr
                                       = find_floating_peak_result( orig_and_energy.second );
            if( fpr )
              bystander_fp_energies.push_back( fpr->original_spectrum_cal_energy );
          }
        }

        if( !bystander_fp_energies.empty() )
        {
          // Tolerance for matching floating peak observables to bystander energies.
          // Both are in the original spectrum calibration at this point (after
          // translate_peaks_to_orig_cal), so the difference should be small.
          const double match_tol = 1.0; // keV

          assert( num_solver_observables <= result.observable_peaks.size() );
          const std::vector<PeakDef>::iterator solver_end
                         = result.observable_peaks.begin() + static_cast<long>(num_solver_observables);

          // Only the solver's own peaks are candidates - see `num_solver_observables` above.  The
          //  updated bystanders appended after `solver_end` are what we are de-duplicating *for*,
          //  and are left alone by construction rather than by inspecting their source.
          const auto new_end = std::remove_if( result.observable_peaks.begin(),
            solver_end,
            [&bystander_fp_energies, match_tol]( const PeakDef &peak ) -> bool
            {
              // Only remove unattributed peaks; a stand-in the solver made for a floating peak
              //  carries no source, while a real fit result does.
              if( peak.parentNuclide() || peak.xrayElement() || peak.reaction() )
                return false;

              const double peak_mean = peak.mean();
              for( const double fp_energy : bystander_fp_energies )
              {
                if( std::fabs( peak_mean - fp_energy ) < match_tol )
                  return true;
              }
              return false;
            } );

#if( PERFORM_DEVELOPER_CHECKS )
          if( should_debug_print() )
          {
            const size_t num_removed = std::distance( new_end, solver_end );
            if( num_removed > 0 )
            {
              std::cout << "ExistingPeaksAsFreePeak: removed " << num_removed
                   << " unattributed floating peak observable(s) to prevent"
                   << " duplication with retained/updated bystander user peaks."
                   << std::endl;
            }
          }
#endif

          result.observable_peaks.erase( new_end, solver_end );
        }//if( !bystander_fp_energies.empty() )
      }

      // Developer check: warn if any observable ROI overlaps an active existing ROI.
      // The pre-fit filter_rois_for_existing should prevent this, but RelActAuto may
      // adjust ROI bounds during solving.  We do not modify ROI bounds post-fit since
      // that would make the displayed bounds inconsistent with the actual fit.
#if( PERFORM_DEVELOPER_CHECKS )
      {
        std::set<const PeakContinuum *> removed_continuums;
        for( const std::shared_ptr<const PeakDef> &peak : result.original_peaks_to_remove )
          removed_continuums.insert( peak->continuum().get() );

        std::set<const PeakContinuum *> seen_obs;
        for( const PeakDef &obs_peak : result.observable_peaks )
        {
          if( !obs_peak.continuum() || !seen_obs.insert( obs_peak.continuum().get() ).second )
            continue;

          const double obs_lower = obs_peak.continuum()->lowerEnergy();
          const double obs_upper = obs_peak.continuum()->upperEnergy();

          std::set<const PeakContinuum *> seen_user;
          for( const std::shared_ptr<const PeakDef> &up : user_peaks )
          {
            if( !up || removed_continuums.count( up->continuum().get() ) )
              continue;
            if( !seen_user.insert( up->continuum().get() ).second )
              continue;
            const double user_lower = up->lowerX();
            const double user_upper = up->upperX();
            if( (obs_lower < user_upper) && (obs_upper > user_lower) )
            {
              std::cerr << "WARNING: ExistingPeaksAsFreePeak: observable ROI ["
                   << obs_lower << ", " << obs_upper
                   << "] overlaps active existing ROI [" << user_lower << ", " << user_upper
                   << "] keV" << std::endl;
            }
          }
        }
      }
#endif
    }// if( existing_peaks_as_free ... )

    // Default mode (neither DoNotUseExistingRois nor ExistingPeaksAsFreePeak):
    // Populate original_peaks_to_remove with user peaks whose source matches the fit sources
    // and whose energy falls within a fitted observable ROI.  This ensures the caller removes
    // old peaks before adding new fit results, preventing duplicate peaks at the same location.
    // DoNotUseExistingRois is excluded: existing ROIs are intentionally left untouched.
    if( !do_not_use_existing_rois && !existing_peaks_as_free
       && !user_peaks.empty()
       && RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
    {
      for( const std::shared_ptr<const PeakDef> &peak : user_peaks )
      {
        if( !peak || !peak_source_is_in_fit(peak) )
          continue;

        // Only remove if the fit actually covered this peak's location.
        const double peak_energy = peak->mean();
        for( const PeakDef &obs_peak : result.observable_peaks )
        {
          if( obs_peak.continuum()
             && (peak_energy >= obs_peak.continuum()->lowerEnergy())
             && (peak_energy <= obs_peak.continuum()->upperEnergy()) )
          {
            result.original_peaks_to_remove.push_back( peak );
            break;
          }
        }
      }//for( user_peaks )

      // Process default mode bystander peaks from mixed ROIs.
      // These were added as FloatingPeaks during ROI setup; now update them from FloatingPeakResults
      // and add to observable_peaks (preserving original source attribution).
      for( const std::pair<std::shared_ptr<const PeakDef>, double> &orig_and_energy : default_mode_bystander_peaks )
      {
        const std::shared_ptr<const PeakDef> &orig_peak = orig_and_energy.first;
        const double orig_energy = orig_and_energy.second;

        // Find the observable ROI containing the bystander energy.
        std::shared_ptr<PeakContinuum> roi_continuum;
        double roi_lower = -1.0, roi_upper = -1.0;
        for( PeakDef &obs_peak : result.observable_peaks )
        {
          if( obs_peak.continuum()
             && (orig_energy >= obs_peak.continuum()->lowerEnergy())
             && (orig_energy <= obs_peak.continuum()->upperEnergy()) )
          {
            roi_continuum = obs_peak.continuum();
            roi_lower = obs_peak.continuum()->lowerEnergy();
            roi_upper = obs_peak.continuum()->upperEnergy();
            break;
          }
        }

        if( !roi_continuum )
          continue;  // Bystander not in any observable ROI — leave untouched

        // Only update if the fit produced source peaks in this ROI
        if( !fit_has_source_peak_in_range( roi_lower, roi_upper ) )
          continue;

        // Reconstruct updated bystander from FloatingPeakResult, preserving source info and color.
        const RelActCalcAuto::FloatingPeakResult *fpr = find_floating_peak_result( orig_energy );

        if( !fpr || (fpr->amplitude <= 0.0) )
        {
#if( PERFORM_DEVELOPER_CHECKS )
          if( !fpr )
          {
            std::cerr << "WARNING: Default mode bystander at " << orig_energy
                 << " keV has no matching FloatingPeakResult — keeping original." << std::endl;
          }else
          {
            std::cerr << "WARNING: Default mode bystander at " << orig_energy
                 << " keV has non-positive amplitude (" << fpr->amplitude
                 << ") in FloatingPeakResult — keeping original." << std::endl;
          }
#endif
          continue;
        }

        PeakDef updated_bystander( *orig_peak );

        updated_bystander.setMean( fpr->original_spectrum_cal_energy );
        const double sigma = fpr->fwhm / (2.0 * std::sqrt( 2.0 * std::log( 2.0 ) ));
        updated_bystander.setSigma( sigma );
        updated_bystander.setAmplitude( fpr->amplitude );
        if( fpr->amplitude_uncert > 0.0 )
          updated_bystander.setAmplitudeUncert( fpr->amplitude_uncert );
        if( fpr->fwhm_uncert > 0.0 )
        {
          const double sigma_uncert = fpr->fwhm_uncert / (2.0 * std::sqrt( 2.0 * std::log( 2.0 ) ));
          updated_bystander.setSigmaUncert( sigma_uncert );
        }

        updated_bystander.setContinuum( roi_continuum );

        if( should_debug_print() )
        {
          std::cout << "Default mode: updating bystander at " << orig_energy
               << " keV -> " << fpr->original_spectrum_cal_energy << " keV"
               << " (amp: " << orig_peak->amplitude() << " -> " << fpr->amplitude
               << ", fwhm: " << orig_peak->fwhm() << " -> " << fpr->fwhm << ")"
               << ", source: " << orig_peak->sourceName() << std::endl;
        }

        result.observable_peaks.push_back( updated_bystander );
        result.original_peaks_to_remove.push_back( orig_peak );
      }//for( default_mode_bystander_peaks )

      // Bug 3 fix: Expand original_peaks_to_remove to include all sibling peaks sharing a
      // PeakContinuum with any removed peak.  Without this, sibling peaks from other sources
      // would keep stale ROI bounds that may overlap with new fit results.
      // For each removed same-source peak, map its old continuum to the replacement observable
      // ROI's continuum, then use that mapping for sibling peaks.
      if( !result.original_peaks_to_remove.empty() )
      {
        // Map old user-peak continuum → replacement observable peak continuum
        std::map<const PeakContinuum *, std::shared_ptr<PeakContinuum>> continuum_replacement_map;
        for( const std::shared_ptr<const PeakDef> &removed : result.original_peaks_to_remove )
        {
          const PeakContinuum *old_cont = removed->continuum().get();
          if( continuum_replacement_map.count( old_cont ) )
            continue;

          // Find the observable peak whose ROI covers this removed peak's energy
          const double removed_energy = removed->mean();
          for( PeakDef &obs_peak : result.observable_peaks )
          {
            if( obs_peak.continuum()
               && (removed_energy >= obs_peak.continuum()->lowerEnergy())
               && (removed_energy <= obs_peak.continuum()->upperEnergy()) )
            {
              continuum_replacement_map[old_cont] = obs_peak.continuum();
              break;
            }
          }
        }//for( removed peaks )

        // Find sibling user peaks sharing the same continuum as removed peaks
        std::set<const PeakContinuum *> removed_continuums;
        for( const std::shared_ptr<const PeakDef> &rm : result.original_peaks_to_remove )
          removed_continuums.insert( rm->continuum().get() );

        for( const std::shared_ptr<const PeakDef> &peak : user_peaks )
        {
          if( !peak )
            continue;
          if( removed_continuums.count( peak->continuum().get() ) == 0 )
            continue;

          // Check not already in original_peaks_to_remove (pointer identity)
          const bool already_marked = std::any_of(
            result.original_peaks_to_remove.begin(), result.original_peaks_to_remove.end(),
            [&peak]( const std::shared_ptr<const PeakDef> &rm ) -> bool { return rm == peak; }
          );
          if( already_marked )
            continue;

          // Sibling shares a continuum with a removed peak, so it must be removed too (all peaks
          // of a ROI are removed/replaced together).
          result.original_peaks_to_remove.push_back( peak );

          // Only re-add a copy for OTHER-source siblings (bystanders whose peak we must preserve
          // across the ROI replacement).  A same-source sibling was already evaluated by the fit:
          // if it were observable it would be in observable_peaks by now, so re-adding a stale
          // copy would resurrect a peak the fit intentionally dropped as insignificant (with an
          // amplitude that is inconsistent with the replacement continuum to boot).
          if( peak_source_is_in_fit( peak ) )
          {
            if( should_debug_print() )
            {
              std::cout << "Default mode: removing same-source sibling at " << peak->mean()
                   << " keV without replacement (fit found it unobservable)." << std::endl;
            }
            continue;
          }

          const PeakContinuum *old_cont = peak->continuum().get();
          const std::map<const PeakContinuum *, std::shared_ptr<PeakContinuum>>::const_iterator replacement_it
            = continuum_replacement_map.find( old_cont );

          if( replacement_it != continuum_replacement_map.end() )
          {
            PeakDef sibling_copy( *peak );
            sibling_copy.setContinuum( replacement_it->second );
            result.observable_peaks.push_back( sibling_copy );
          }else
          {
            // No replacement ROI found; add the peak as-is (shouldn't normally happen)
            PeakDef sibling_copy( *peak );
            result.observable_peaks.push_back( sibling_copy );
#if( PERFORM_DEVELOPER_CHECKS )
            std::cerr << "WARNING: sibling peak at " << peak->mean()
                 << " keV has no replacement continuum in observable_peaks." << std::endl;
#endif
          }
        }//for( user_peaks )

#if( PERFORM_DEVELOPER_CHECKS )
        // Verify: no user peak shares a continuum with a removed peak unless it is also removed
        for( const std::shared_ptr<const PeakDef> &peak : user_peaks )
        {
          if( !peak )
            continue;

          const bool is_removed = std::any_of(
            result.original_peaks_to_remove.begin(), result.original_peaks_to_remove.end(),
            [&peak]( const std::shared_ptr<const PeakDef> &rm ) -> bool { return rm == peak; }
          );
          if( is_removed )
            continue;

          for( const std::shared_ptr<const PeakDef> &rm : result.original_peaks_to_remove )
          {
            if( peak->continuum().get() == rm->continuum().get() )
            {
              std::cerr << "BUG: user peak at " << peak->mean()
                   << " keV shares continuum with removed peak at " << rm->mean()
                   << " keV but was not also removed." << std::endl;
              assert( false );
            }
          }
        }
#endif
      }//if( !result.original_peaks_to_remove.empty() )
    }// default mode

    // Assign escape peak relationships for high-energy gammas if appropriate
    // This checks fit peaks (including observable_peaks) and assigns S.E. and D.E. relationships
    // to significant escape peaks of high-energy lines like Th232 2614 keV
    assign_escape_peak_relationships( result.fit_peaks, fit_norm_peaks, det_type );
    assign_escape_peak_relationships( result.observable_peaks, fit_norm_peaks, det_type );
    assign_escape_peak_relationships( result.uncombined_fit_peaks, fit_norm_peaks, det_type );

    // Mirror the use-flag defaults that PeakModel::setNuclide applies on the manual
    // double-click path, so "Fit Source" peaks are pre-selected for shielding/source
    // and manual rel-eff fits where appropriate.
    apply_default_use_flags( result.fit_peaks );
    apply_default_use_flags( result.observable_peaks );
    apply_default_use_flags( result.uncombined_fit_peaks );

  }catch( const std::exception &e )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    result.error_message = e.what();
  }

  return result;
}//fit_peaks_for_nuclide_relactauto


const PeakFitForNuclideConfig &PeakFitForNuclideConfig::default_config( const PeakFitUtils::CoarseResolutionType det_type )
{
  static const PeakFitForNuclideConfig s_default_hpge_config;

  if( det_type == PeakFitUtils::CoarseResolutionType::High )
    return s_default_hpge_config;

  // TODO: add a constructor for PeakFitForNuclideConfig, so we can avoid all this mutex stuff
  static std::mutex s_have_inited_non_hpge_config_mutex;
  static bool s_have_inited_non_hpge_config = false;
  static PeakFitForNuclideConfig s_default_non_hpge_config;

  std::unique_lock<std::mutex> lock( s_have_inited_non_hpge_config_mutex );
  if( s_have_inited_non_hpge_config )
    return s_default_non_hpge_config;

  // Settings from a single Genetic Algorithm optimization, using only 50 individuals, and 8 sources - so could be a lot better!
  s_default_non_hpge_config.fwhm_functional_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
  s_default_non_hpge_config.rel_eff_manual_base_rel_eff_uncert=0.308957;
  s_default_non_hpge_config.initial_nuc_match_cluster_num_sigma=1.41369;
  s_default_non_hpge_config.manual_eff_cluster_num_sigma=3.80554;
  // Manual rel-eff form/order is now selected per spectrum by AICc (kappa below); the former
  // per-peak-count form/order fields are gone.
  s_default_non_hpge_config.manual_releff_aicc_penalty=2.0;
  s_default_non_hpge_config.cont_order_aicc_penalty=2.0;
  s_default_non_hpge_config.manual_keep_significance_z=5.26621;
  s_default_non_hpge_config.manual_rel_eff_sol_min_fwhm_roi=2.0;
  s_default_non_hpge_config.manual_rel_eff_sol_max_fwhm=19.9579;
  // Adaptive-extent seeds: the old fixed half-widths (~2.1-2.75 FWHM) split into an always-kept
  // core plus data-driven extension up to the cap.  Re-tuned by the GA.
  s_default_non_hpge_config.manual_roi_core_num_fwhm=1.25;
  s_default_non_hpge_config.fwhm_form = RelActCalcAuto::FwhmForm::Berstein_3;
  s_default_non_hpge_config.rel_eff_auto_base_rel_eff_uncert=0.191199;
  s_default_non_hpge_config.auto_rel_eff_cluster_num_sigma=4.0;
  s_default_non_hpge_config.auto_keep_significance_z=6.38959;
  s_default_non_hpge_config.auto_roi_core_num_fwhm=1.25;
  s_default_non_hpge_config.roi_extend_z=2.0;
  s_default_non_hpge_config.roi_max_num_fwhm=4.0;
  s_default_non_hpge_config.auto_rel_eff_sol_max_fwhm=12.2638;
  s_default_non_hpge_config.merge_tail_z=2.0;
  s_default_non_hpge_config.merge_clean_gap_fwhm=1.0;
  s_default_non_hpge_config.auto_rel_eff_sol_min_fwhm_roi=0.691922;
  // LnXLnY (pure ln(x) polynomial, no 1/x term) is used instead of FramEmpirical order 1: the latter
  // is exp(a + b/x^2), whose lone 1/x^2 term runs away at low energy (e.g. on NaI it explodes below
  // ~100 keV and drives the source activity to 0, so no peaks are fit).  LnXLnY order 2 matches the
  // HPGe / struct default, is stable, and does not depend on the detector efficiency.
  s_default_non_hpge_config.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnXLnY;
  s_default_non_hpge_config.rel_eff_eqn_order=2;
  s_default_non_hpge_config.desperation_phys_model_atomic_number=41.2303;
  s_default_non_hpge_config.desperation_phys_model_areal_density_g_per_cm2=13.7506;
  s_default_non_hpge_config.nucs_of_el_same_age = true;
  s_default_non_hpge_config.phys_model_use_hoerl = false;
  s_default_non_hpge_config.fit_energy_cal = true;  //manually changed from `false`
  // Equivalent-z seed for the unified LR test: the old delta-chi2 gate of 24.235 with one peak
  // dof corresponds to z = sqrt(24.235) ~ 4.9.  Re-tuned by the GA.
  s_default_non_hpge_config.roi_significance_z=4.9;
  s_default_non_hpge_config.observable_peak_initial_significance_threshold=4.2546;
  s_default_non_hpge_config.observable_peak_final_significance_threshold=3.73082;
  s_default_non_hpge_config.step_cont_min_peak_significance=61.9867;
  // step_trial_chi2_margin left at the struct default: the GA has not yet tuned the new
  // step-trial parameterization (the former step_cont_left_right_nsigma gene it replaces was
  // tuned against the old probe-window test).

  s_have_inited_non_hpge_config = true;

  return s_default_non_hpge_config;
}

PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const Wt::WFlags<FitSrcPeaksOptions> options,
  const PeakFitForNuclideConfig &config,
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs )
{

  std::vector<RelActCalcAuto::NucInputInfo> base_nuclides;
  base_nuclides.reserve( sources.size() );

  for( const RelActCalcAuto::SrcVariant &src : sources )
  {
    RelActCalcAuto::NucInputInfo nuc_info;
    nuc_info.age = get_source_age( src, -1.0 );
    nuc_info.source = src;
    nuc_info.fit_age = false;  // not currently exposed in UI
    base_nuclides.push_back( nuc_info );
  }

  return fit_peaks_for_nuclides( auto_search_peaks, foreground, base_nuclides,
                                user_peaks, background, drf_input, options, config, peak_fit_prefs );
}
  
  
PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
  std::shared_ptr<const SpecUtils::Measurement> long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const Wt::WFlags<FitSrcPeaksOptions> options,
  const PeakFitForNuclideConfig &config,
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs )
{
  assert( peak_fit_prefs );
  // Diagnostics are per top-level fit; discard any abandoned thread-local construction from a
  // prior exception before beginning this request.
  static_cast<void>( detail::take_automatic_roi_diagnostics() );
  const PeakFitUtils::CoarseResolutionType det_type = peak_fit_prefs
    ? peak_fit_prefs->m_det_type
    : PeakFitUtils::coarse_det_type( foreground, nullptr );

  PeakFitResult result;

  // NOTE: `local_debug_printout` is left at its ambient value (false by default) - callers such as
  // the standalone debug harnesses set it true explicitly.  It must NOT be force-enabled here: this
  // function runs on many parallel worker threads during GA optimization, where an unsynchronized
  // write to the shared flag is a data race and floods stdout/stderr.

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "Failed to open SandiaDecayDataBase";
    return result;
  }

  // Validate sources
  if( sources.empty() && !options.test(FitSrcPeaksOptions::FitNormBkgrndPeaks) )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
    result.error_message = "No sources provided";
    return result;
  }

  for( const RelActCalcAuto::NucInputInfo &src : sources )
  {
    if( RelActCalcAuto::is_null( src.source ) )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Null source in sources vector";
      return result;
    }
  }

  // Use input DRF or create a copy we can modify
  std::shared_ptr<const DetectorPeakResponse> drf = drf_input;

  vector<string> local_warnings;  // Set if we use the fallback activity estimation

  try
  {
    // Step 1: Determine FWHM functional form from auto-search peaks or DRF
    DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm = config.fwhm_functional_form;
    double lower_fwhm_energy = -1.0, upper_fwhm_energy = -1.0;
    std::vector<float> fwhm_coefficients, fwhm_uncerts;

    if( !drf || !drf->isValid() || !drf->hasResolutionInfo() || (auto_search_peaks.size() > 6) )
    {
      bool got_fwhm_fcn = false;
      
      // If we have peaks, estimate FWHM from peak widths. Otherwise fall back to a generic detector.
      if( !auto_search_peaks.empty() )
      {
        try
        {
          const int num_auto_peaks = static_cast<int>(auto_search_peaks.size());
          int sqrtEqnOrder = (std::min)( 6, num_auto_peaks / (1 + (num_auto_peaks > 3)) );
          if( auto_search_peaks.size() < 3 )
            sqrtEqnOrder = static_cast<int>( auto_search_peaks.size() );
          
          std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> auto_search_peaks_dq
          = std::make_shared<const std::deque<std::shared_ptr<const PeakDef>>>( begin(auto_search_peaks), end(auto_search_peaks) );
          
          MakeDrfFit::performResolutionFit( auto_search_peaks_dq, fwhmFnctnlForm, sqrtEqnOrder, fwhm_coefficients, fwhm_uncerts );
          auto_search_peaks_dq = MakeDrfFit::removeOutlyingWidthPeaks( auto_search_peaks_dq, fwhmFnctnlForm, fwhm_coefficients );
          MakeDrfFit::performResolutionFit( auto_search_peaks_dq, fwhmFnctnlForm, sqrtEqnOrder, fwhm_coefficients, fwhm_uncerts );
          
          // Set energy range based on peaks used for FWHM fit
          if( !auto_search_peaks_dq->empty() )
          {
            lower_fwhm_energy = auto_search_peaks_dq->front()->mean();
            upper_fwhm_energy = auto_search_peaks_dq->back()->mean();
          }
          
          got_fwhm_fcn = true;
        }catch( std::exception &e )
        {
          cerr << "Failed to fit FWHM functional form in fit_peaks_for_nuclides: " << e.what() << ". Will use generic FWHM." << endl;
          got_fwhm_fcn = false;  // fall through to the generic-detector FWHM fallback below
        }
      }//if( !auto_search_peaks.empty() )
      
      if( !got_fwhm_fcn )
      {
        // Prefer the caller-supplied DRF's own resolution info when it has some (we may only be
        // here because the peak-based FWHM fit threw - the DRF was valid but >6 auto-search peaks
        // made us try the peak fit first); replacing a real DRF with a generic detector would also
        // corrupt the efficiency/rel-eff usage of `drf` downstream.
        if( drf && drf->isValid() && drf->hasResolutionInfo() )
        {
          fwhmFnctnlForm = drf->resolutionFcnType();
          fwhm_coefficients = drf->resolutionFcnCoefficients();
          lower_fwhm_energy = drf->lowerEnergy();
          upper_fwhm_energy = drf->upperEnergy();
          local_warnings.push_back( "Estimating resolution from peaks failed; using the detector response function's FWHM parameters." );
        }else
        {
          // With no peaks and no DRF resolution info, use generic detector resolution coefficients.
          switch( det_type )
          {
            case PeakFitUtils::CoarseResolutionType::High:
              drf = DetectorPeakResponse::getGenericHPGeDetector();
              break;
            case PeakFitUtils::CoarseResolutionType::LaBr:
            case PeakFitUtils::CoarseResolutionType::MedRes:
              drf = DetectorPeakResponse::getGenericLaBrDetector();
              break;
            case PeakFitUtils::CoarseResolutionType::CZT:
              drf = DetectorPeakResponse::getGenericCZTGeneralDetector();
              break;
            case PeakFitUtils::CoarseResolutionType::Low:
            case PeakFitUtils::CoarseResolutionType::LowOrMedRes:
            case PeakFitUtils::CoarseResolutionType::Unknown:
            default:
              drf = DetectorPeakResponse::getGenericNaIDetector();
              break;
          }//switch( det_type )

          if( drf && drf->isValid() && drf->hasResolutionInfo() )
          {
            fwhmFnctnlForm = drf->resolutionFcnType();
            fwhm_coefficients = drf->resolutionFcnCoefficients();
            lower_fwhm_energy = drf->lowerEnergy();
            upper_fwhm_energy = drf->upperEnergy();
            local_warnings.push_back( "No peaks were available to estimate resolution; using generic detector FWHM parameters." );
          }
        }//if( input DRF has resolution info ) / else
      }//if( !auto_search_peaks.empty() ) / else
    }else
    {
      fwhmFnctnlForm = drf->resolutionFcnType();
      fwhm_coefficients = drf->resolutionFcnCoefficients();

      // Get energy range from detector response function
      lower_fwhm_energy = drf->lowerEnergy();
      upper_fwhm_energy = drf->upperEnergy();
    }

    // Validate that we have valid FWHM coefficients
    if( fwhm_coefficients.empty() )
    {
      result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      result.error_message = "Failed to determine FWHM coefficients - unable to proceed with peak fitting";
      return result;
    }

    // Check that coefficients are finite
    for( size_t i = 0; i < fwhm_coefficients.size(); ++i )
    {
      if( !std::isfinite( fwhm_coefficients[i] ) )
      {
        result.status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
        result.error_message = "FWHM coefficient[" + std::to_string(i) + "] is not finite (value="
                               + std::to_string(fwhm_coefficients[i]) + ")";
        return result;
      }
    }

    // Find valid energy range, clamped to a physically-valid low-energy floor (see low_energy_analysis_floor).
    const std::pair<double,double> raw_valid_range = find_valid_energy_range( foreground );
    const double low_e_floor = low_energy_analysis_floor( drf_input, det_type );
    const double min_valid_energy = (low_e_floor < raw_valid_range.second)
                                    ? std::max( raw_valid_range.first, low_e_floor ) : raw_valid_range.first;
    const double max_valid_energy = raw_valid_range.second;

    if( should_debug_print() )
      std::cout << "fit_peaks_for_nuclides: valid energy range = ["
                << min_valid_energy << ", " << max_valid_energy << "] keV" << std::endl;

    /*
    // Determine energy range for gamma lines
    double highest_energy_gamma = 0.0, lowest_energy_gamma = std::numeric_limits<double>::max();

    
    for( const RelActCalcAuto::NucInputInfo &src : sources )
    {
      const std::vector<SandiaDecay::EnergyRatePair> photons
          = get_source_photons( src.source, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits, src.age );
      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        highest_energy_gamma = (std::max)( highest_energy_gamma, photon.energy );
        lowest_energy_gamma = (std::min)( lowest_energy_gamma, photon.energy );
      }
    }

    const double energy_pad = (det_type == PeakFitUtils::CoarseResolutionType::High) ? 5.0 : 25.0;
    lowest_energy_gamma = (std::max)( lowest_energy_gamma - energy_pad, (double)foreground->gamma_energy_min() );
    highest_energy_gamma = (std::min)( highest_energy_gamma + energy_pad, (double)foreground->gamma_energy_max() );
    */

    // Step 2, 3 & 4: Estimate initial ROIs using RelActManual with multiple fallbacks
    // This function internally:
    // - Converts auto_search_peaks to RelActManual format and matches to sources
    // - Falls back to estimate_initial_rois_without_peaks() if no peaks match
    // - Fits relative efficiency curve and clusters gammas into ROIs
    // - Falls back to estimate_initial_rois_fallback() if RelActManual fails
    GammaClusteringSettings manual_settings = config.get_manual_clustering_settings();
    const bool r6_enabled = !options.test(FitSrcPeaksOptions::FitNormBkgrndPeaks)
        && !options.test(FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse)
        && !options.test(FitSrcPeaksOptions::DisableAutoInterfererFit);
    manual_settings.use_automatic_roi_policy = !r6_enabled;
    const std::vector<std::shared_ptr<const PeakDef>> local_unfit_auto_peaks
      = compute_unfit_auto_peaks( auto_search_peaks, user_peaks );
    std::vector<std::shared_ptr<const PeakDef>> background_auto_search_peaks;
    if( long_background && !sources.empty() )
    {
      try
      {
        background_auto_search_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks(
          long_background, nullptr, nullptr, true, peak_fit_prefs );
      }catch( const std::exception &e )
      {
        local_warnings.push_back( "Unable to search supplied background for RelActManual seeding: "
                                  + std::string(e.what()) );
      }
    }
    const std::vector<std::shared_ptr<const PeakDef>> no_background_auto_search_peaks;

    string fallback_warning;
    std::vector<std::pair<double,double>> initial_modeled_peak_candidates;
    std::vector<RelActCalcManual::GenericPeakInfo> source_anchor_candidates;
    std::vector<RelActCalcAuto::RoiRange> clean_source_rois;
    bool has_provisional_fallback_source_anchors = false;
    const vector<RelActCalcAuto::RoiRange> source_rois = sources.empty()
      ? vector<RelActCalcAuto::RoiRange>{}
      : estimate_initial_rois_using_relactmanual(
          auto_search_peaks, foreground, long_background, background_auto_search_peaks,
          sources, drf, det_type,
          fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
          min_valid_energy, max_valid_energy, manual_settings,
          config, local_unfit_auto_peaks, fallback_warning,
          &initial_modeled_peak_candidates, &source_anchor_candidates, &clean_source_rois,
          &has_provisional_fallback_source_anchors
        );
    
    if( !fallback_warning.empty() )
      local_warnings.push_back( fallback_warning );

    // We need to fit NORM peaks on a different Rel Eff curve than source nuclides (since they will have two differnt
    //  efficiency curves), and then combine the ROIs together.
    vector<RelActCalcAuto::RoiRange> norm_rois;
    if( options.test(FitSrcPeaksOptions::FitNormBkgrndPeaks)
       || options.test(FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse) )
    {
      long_background = nullptr; //We dont want to do background subtraction if we are trying to fit NORM peaks.
      
      const vector<RelActCalcAuto::NucInputInfo> norm_sources = get_norm_sources( sources, config.norm_css_color );
      if( !norm_sources.empty() )
      {
        try
        {
          fallback_warning.clear();
          norm_rois = estimate_initial_rois_using_relactmanual(
            auto_search_peaks, foreground, nullptr, no_background_auto_search_peaks,
            norm_sources, drf, det_type,
            fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
            min_valid_energy, max_valid_energy, manual_settings,
            config, local_unfit_auto_peaks, fallback_warning,
            &initial_modeled_peak_candidates, nullptr
          );
          
          if( !fallback_warning.empty() )
            local_warnings.push_back( fallback_warning );
        }catch( std::exception &e )
        {
          local_warnings.push_back( "Unable to estimate initial ROIs for NORM peaks: " + std::string(e.what()) );
        }//try / catch to get norm ROIs
      }//if( !norm_sources.empty() )
    }//if( use NORM peaks )
    
    
    vector<RelActCalcAuto::RoiRange> initial_rois;
    {// Begin combine `source_rois` and `norm_rois`
      // Combine source and NORM ROIs, then merge any that overlap
      vector<RelActCalcAuto::RoiRange> all_rois = source_rois;
      all_rois.insert( end(all_rois), begin(norm_rois), end(norm_rois) );


      const bool have_fwhm_range = (lower_fwhm_energy > 0.0) && (upper_fwhm_energy > 0.0) && (lower_fwhm_energy < upper_fwhm_energy);

      vector<InitialRoi> all_roi_infos;
      for( const RelActCalcAuto::RoiRange &roi : all_rois )
      {
        InitialRoi roi_info;
        roi_info.roi = roi;
        roi_info.center_energy = 0.5 * (roi.upper_energy + roi.lower_energy);
        const double fwhm_eval_energy = have_fwhm_range
            ? std::clamp( roi_info.center_energy, lower_fwhm_energy, upper_fwhm_energy )
            : roi_info.center_energy;
        roi_info.fwhm = DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(fwhm_eval_energy), fwhmFnctnlForm, fwhm_coefficients );

        float min_sigma_width, max_sigma_width;
        expected_peak_width_limits( roi_info.fwhm,
          det_type,
          foreground, min_sigma_width, max_sigma_width );

        if( roi_info.fwhm < (min_sigma_width * PhysicalUnits::fwhm_nsigma) )
          roi_info.fwhm = min_sigma_width * PhysicalUnits::fwhm_nsigma;
        if( roi_info.fwhm > (max_sigma_width * PhysicalUnits::fwhm_nsigma) )
          roi_info.fwhm = max_sigma_width * PhysicalUnits::fwhm_nsigma;

        // Find the largest auto_search_peak within the ROI for the amplitude estimate; when one
        // exists, also anchor center_energy on its mean rather than the ROI midpoint - the
        // clean-gap merge test and split-point constraints model the ROI's signal as a Gaussian
        // at center_energy, and a real peak position grounds that far better than the geometric
        // midpoint of a (possibly asymmetric) ROI.
        for( const std::shared_ptr<const PeakDef> &peak : auto_search_peaks )
        {
          if( (peak->mean() >= roi_info.roi.lower_energy)
            && (peak->mean() <= roi_info.roi.upper_energy)
            && (peak->amplitude() > roi_info.estimated_amplitude) )
          {
            roi_info.estimated_amplitude = peak->amplitude();
            roi_info.center_energy = peak->mean();
          }
        }

        if( manual_settings.use_automatic_roi_policy )
        {
          for( const std::pair<double,double> &candidate : initial_modeled_peak_candidates )
          {
            if( (candidate.first >= roi.lower_energy) && (candidate.first <= roi.upper_energy) )
            {
              roi_info.modeled_energies.push_back( candidate.first );
              roi_info.modeled_areas.push_back( candidate.second );
            }
          }
          if( !roi_info.modeled_energies.empty() )
          {
            const std::vector<double>::const_iterator dominant = std::max_element(
                std::begin(roi_info.modeled_areas), std::end(roi_info.modeled_areas) );
            const size_t dominant_index = static_cast<size_t>(
                dominant - std::begin(roi_info.modeled_areas) );
            roi_info.center_energy = roi_info.modeled_energies[dominant_index];
            roi_info.estimated_amplitude = *dominant;
          }
        }

        all_roi_infos.push_back( roi_info );
      }//for( const RelActCalcAuto::RoiRange &roi : all_rois )

      const auto initial_policy_fwhm = [=, &fwhm_coefficients]( const double energy ) {
        const double eval_energy = have_fwhm_range
          ? std::clamp(energy, lower_fwhm_energy, upper_fwhm_energy) : energy;
        return static_cast<double>( DetectorPeakResponse::peakResolutionFWHM(
            static_cast<float>(eval_energy), fwhmFnctnlForm, fwhm_coefficients) );
      };
      detail::GlobalContinuumEstimate initial_merge_continuum;
      if( manual_settings.use_automatic_roi_policy )
      {
        initial_merge_continuum = detail::make_global_continuum(
            foreground, initial_policy_fwhm, det_type, min_valid_energy, max_valid_energy );
      }
      const std::vector<std::shared_ptr<const PeakDef>> no_unfit_peaks;
      const std::vector<std::shared_ptr<const PeakDef>> &initial_merge_unfit
        = manual_settings.use_automatic_roi_policy ? local_unfit_auto_peaks : no_unfit_peaks;
      initial_rois = merge_rois( all_roi_infos, config, initial_merge_unfit, foreground,
          initial_merge_continuum.valid() ? &initial_merge_continuum : nullptr, nullptr,
          "initial source/NORM merge", manual_settings.use_automatic_roi_policy );
    }// End combine `source_rois` and `norm_rois`
    
    
    // Call RelActAuto with initial_rois
    result = fit_peaks_for_nuclide_relactauto(
      auto_search_peaks, foreground, sources,
      initial_rois, norm_rois.empty() ? clean_source_rois : initial_rois,
      initial_modeled_peak_candidates, source_anchor_candidates,
      has_provisional_fallback_source_anchors,
      user_peaks, long_background,
      drf, options, config,
      fwhmFnctnlForm, fwhm_coefficients, det_type,
      lower_fwhm_energy, upper_fwhm_energy,
      peak_fit_prefs
    );

  }catch( const std::exception &e )
  {
    result.status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    result.error_message = e.what();
  }

  // Add any local warnings
  result.warnings.insert( end(result.warnings), begin(local_warnings), end(local_warnings) );
  std::vector<AutomaticRoiDecisionDiagnostic> pending_roi_diagnostics
    = detail::take_automatic_roi_diagnostics();
  result.automatic_roi_diagnostics.insert( std::end(result.automatic_roi_diagnostics),
      std::begin(pending_roi_diagnostics), std::end(pending_roi_diagnostics) );

  return result;
}//fit_peaks_for_nuclides


int debug_fit_peaks_for_nuclides()
{
  using namespace std;

  // Enable internal debug traces throughout the pipeline
  local_debug_printout = true;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
  {
    cerr << "debug_fit_peaks_for_nuclides: Failed to open SandiaDecayDataBase" << endl;
    local_debug_printout = false;
    return 1;
  }

  const string data_dir = InterSpec::staticDataDirectory();
  const string spec_dir = SpecUtils::append_path( data_dir,
    "reference_spectra/Common_Field_Nuclides/RadEagle NaI 3x1" );

  // Load foreground
  const string fg_path = SpecUtils::append_path( spec_dir, "Xe133_Shielded.txt" );
  SpecUtils::SpecFile fg_file;
  if( !fg_file.load_file( fg_path, SpecUtils::ParserType::Auto ) )
  {
    cerr << "debug_fit_peaks_for_nuclides: Failed to load foreground: " << fg_path << endl;
    local_debug_printout = false;
    return 1;
  }

  shared_ptr<const SpecUtils::Measurement> foreground = fg_file.measurement_at_index( size_t(0) );
  if( !foreground || !foreground->num_gamma_channels() )
  {
    cerr << "debug_fit_peaks_for_nuclides: Foreground has no gamma data" << endl;
    local_debug_printout = false;
    return 1;
  }

  cout << "\n=== Foreground Spectrum ===" << endl;
  cout << "  File: " << fg_path << endl;
  cout << "  Channels: " << foreground->num_gamma_channels() << endl;
  cout << "  Energy range: " << foreground->gamma_energy_min()
       << " - " << foreground->gamma_energy_max() << " keV" << endl;
  cout << "  Live time: " << foreground->live_time() << " s" << endl;
  cout << "  Real time: " << foreground->real_time() << " s" << endl;

  // Load background
  const string bg_path = SpecUtils::append_path( spec_dir, "background.txt" );
  SpecUtils::SpecFile bg_file;
  shared_ptr<const SpecUtils::Measurement> background;
  if( bg_file.load_file( bg_path, SpecUtils::ParserType::Auto ) )
  {
    background = bg_file.measurement_at_index( size_t(0) );
    if( background && background->num_gamma_channels() )
    {
      cout << "\n=== Background Spectrum ===" << endl;
      cout << "  Channels: " << background->num_gamma_channels() << endl;
      cout << "  Live time: " << background->live_time() << " s" << endl;
    }else
    {
      cerr << "  Warning: Background loaded but has no gamma data" << endl;
      background = nullptr;
    }
  }else
  {
    cerr << "  Warning: Failed to load background: " << bg_path << endl;
  }

  // Run auto peak search
  // TODO: derive det_type from the spectrum, instead of hardcoding
  const PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::Low;
  std::shared_ptr<PeakFitDetPrefs> peak_fit_prefs = std::make_shared<PeakFitDetPrefs>();
  peak_fit_prefs->m_det_type = det_type;

  cout << "\n=== Auto Peak Search ===" << endl;
  const vector<shared_ptr<const PeakDef>> auto_peaks
    = ExperimentalAutomatedPeakSearch::search_for_peaks( foreground, nullptr, nullptr, true, peak_fit_prefs );

  cout << "  Found " << auto_peaks.size() << " auto-search peaks:" << endl;
  for( size_t i = 0; i < auto_peaks.size(); ++i )
  {
    const PeakDef &p = *auto_peaks[i];
    cout << "    [" << i << "] mean=" << std::fixed << std::setprecision(1) << p.mean()
         << " keV, sigma=" << std::setprecision(2) << p.sigma()
         << " keV, amplitude=" << std::setprecision(1) << p.amplitude()
         << ", area=" << std::setprecision(1) << p.peakArea() << endl;
  }

  // Set up sources
  const SandiaDecay::Nuclide * const xe133 = db->nuclide( "Xe133" );
  const SandiaDecay::Nuclide * const xe133m = db->nuclide( "Xe133m" );
  if( !xe133 || !xe133m )
  {
    cerr << "debug_fit_peaks_for_nuclides: Failed to find Xe133 or Xe133m nuclide" << endl;
    local_debug_printout = false;
    return 1;
  }

  vector<RelActCalcAuto::NucInputInfo> sources( 2 );
  sources[0].source = xe133;
  sources[0].age = get_source_age( RelActCalcAuto::SrcVariant{xe133}, -1.0 );
  sources[0].fit_age = false;
  sources[1].source = xe133m;
  sources[1].age = get_source_age( RelActCalcAuto::SrcVariant{xe133m}, -1.0 );
  sources[1].fit_age = false;

  // Print gamma lines for each source
  const double activity = GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits;
  for( size_t s = 0; s < sources.size(); ++s )
  {
    const string src_name = RelActCalcAuto::to_name( sources[s].source );
    const double age = get_source_age( sources[s].source, sources[s].age );
    cout << "\n=== Source: " << src_name << " (age=" << std::setprecision(1)
         << age / PhysicalUnits::second << " s) ===" << endl;

    const vector<SandiaDecay::EnergyRatePair> photons
      = get_source_photons( sources[s].source, activity, age );

    // Sort by rate descending for display
    vector<SandiaDecay::EnergyRatePair> sorted_photons = photons;
    sort( begin(sorted_photons), end(sorted_photons),
      []( const SandiaDecay::EnergyRatePair &a, const SandiaDecay::EnergyRatePair &b ){
        return a.numPerSecond > b.numPerSecond;
    } );

    const size_t num_to_print = min( sorted_photons.size(), size_t(20) );
    cout << "  Top " << num_to_print << " gammas (of " << photons.size() << " total):" << endl;
    for( size_t i = 0; i < num_to_print; ++i )
    {
      const double br = sorted_photons[i].numPerSecond / activity;
      cout << "    " << std::setprecision(2) << sorted_photons[i].energy
           << " keV, BR=" << std::scientific << std::setprecision(4) << br
           << std::fixed << endl;
    }
  }//for( each source )

  // Recompute chi2 against the original foreground and print peaks
  const auto print_peaks = [&foreground]( const string &label, vector<PeakDef> peaks ){
    // Recompute chi2/dof against original foreground so values are trustworthy
    if( !peaks.empty() )
      set_chi2_dof( foreground, peaks, 0, peaks.size() );

    cout << "\n  " << label << ": " << peaks.size() << " peaks" << endl;
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      const PeakDef &p = peaks[i];
      cout << "    [" << i << "] mean=" << std::fixed << std::setprecision(2) << p.mean()
           << " keV, sigma=" << std::setprecision(2) << p.sigma()
           << ", chi2dof=" << std::setprecision(3) << p.chi2dof()
           << ", area=" << std::setprecision(1) << p.peakArea()
           << " +/- " << std::setprecision(1) << p.peakAreaUncert();

      if( p.parentNuclide() )
        cout << ", nuc=" << p.parentNuclide()->symbol;
      if( p.hasSourceGammaAssigned() )
        cout << ", gamma=" << std::setprecision(2) << p.gammaParticleEnergy() << " keV";

      if( p.continuum() )
      {
        cout << ", ROI=[" << std::setprecision(1) << p.continuum()->lowerEnergy()
             << ", " << p.continuum()->upperEnergy() << "]";
      }
      cout << endl;
    }
  };

  const vector<shared_ptr<const PeakDef>> user_peaks; // empty

  // Run with energy calibration fitting enabled
  {
    cout << "\n=== Fit WITH energy cal ===" << endl;
    PeakFitForNuclideConfig config = PeakFitForNuclideConfig::default_config( det_type );
    config.fit_energy_cal = true;
    const Wt::WFlags<FitSrcPeaksOptions> options; // no flags = allow energy cal

    const PeakFitResult result = fit_peaks_for_nuclides(
      auto_peaks, foreground, sources, user_peaks,
      background, nullptr, options, config, peak_fit_prefs
    );

    cout << "  Status: " << static_cast<int>( result.status ) << endl;
    if( !result.error_message.empty() )
      cout << "  Error: " << result.error_message << endl;
    for( const string &w : result.warnings )
      cout << "  Warning: " << w << endl;

    print_peaks( "fit_peaks (WITH energy cal)", result.fit_peaks );
    print_peaks( "observable_peaks (WITH energy cal)", result.observable_peaks );
  }

  // Run without energy calibration fitting
  {
    cout << "\n=== Fit WITHOUT energy cal ===" << endl;
    PeakFitForNuclideConfig config = PeakFitForNuclideConfig::default_config( det_type );
    config.fit_energy_cal = false;
    const Wt::WFlags<FitSrcPeaksOptions> options = FitSrcPeaksOptions::DoNotVaryEnergyCal;

    const PeakFitResult result = fit_peaks_for_nuclides(
      auto_peaks, foreground, sources, user_peaks,
      background, nullptr, options, config, peak_fit_prefs
    );

    cout << "  Status: " << static_cast<int>( result.status ) << endl;
    if( !result.error_message.empty() )
      cout << "  Error: " << result.error_message << endl;
    for( const string &w : result.warnings )
      cout << "  Warning: " << w << endl;

    print_peaks( "fit_peaks (WITHOUT energy cal)", result.fit_peaks );
    print_peaks( "observable_peaks (WITHOUT energy cal)", result.observable_peaks );
  }

  cout << "\n=== Done ===" << endl;

  local_debug_printout = false;
  return 0;
}//debug_fit_peaks_for_nuclides

}//namespace FitPeaksForNuclides
