/*
  Shared types and deterministic helpers for one-dimensional profile-likelihood scans.

  This header deliberately has no dependency on the RelActAuto solver.  It holds the evaluation
  and result types consumed and produced by the vertex-anchored engine in
  RelActCalcAuto_ProfileFit.h, the exact feasible-range projections used to place a scan, the
  automatic weak-quantity trigger, and the deferred better-baseline reselection machinery.  The
  production profile code supplies conditional-fit evaluations, while unit tests supply analytic
  objectives, under the same rules.
*/
#ifndef InterSpec_RelActCalcAuto_ProfileScan_h
#define InterSpec_RelActCalcAuto_ProfileScan_h

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace RelActCalcAutoImp
{
namespace ProfileLikelihood
{

enum class Direction : int { Lower, Upper };

enum class EvaluationStatus : int
{
  Success,
  Failed,
  Canceled,
  BetterBaseline,
  FitCapReached
};

struct Evaluation
{
  EvaluationStatus status = EvaluationStatus::Failed;
  double reported_fraction = 0.0;
  double delta_chi2 = std::numeric_limits<double>::infinity();
  /** The requested outer control point lay beyond a structural feasible limit; the returned
   reported fraction and objective are the independently evaluated limiting point. */
  bool reached_feasible_bound = false;
  std::string diagnostic;
};

enum class ScanStatus : int
{
  Complete,
  BoundaryLimited,
  NonIdentifiable,
  Failed,
  Canceled,
  BetterBaseline,
  FitCapReached
};

struct Endpoint
{
  double reported_fraction = 0.0;
  bool likelihood_crossing = false;
};

struct Interval
{
  double delta_chi2 = 0.0;
  Endpoint lower;
  Endpoint upper;
};

struct Sample
{
  double control_fraction = 0.0;
  double reported_fraction = 0.0;
  double delta_chi2 = 0.0;
  bool reached_feasible_bound = false;
};

struct ScanResult
{
  ScanStatus status = ScanStatus::Failed;
  std::array<Interval,2> intervals;
  std::vector<Sample> samples;
  std::size_t num_fits = 0;
  std::string diagnostic;
};

using Evaluator = std::function<Evaluation(double,Direction,std::size_t)>;

/** Project component and sibling box bounds through `sum(fractions) == 1`.

 Unconstrained siblings are represented by `[0,1]`.  The returned interval is the exact projection
 of the box/simplex intersection onto the selected component (assuming the input windows are
 jointly feasible).
 */
inline std::pair<double,double> project_simplex_component_bounds(
    const double component_lower, const double component_upper,
    const std::vector<std::pair<double,double>> &sibling_bounds )
{
  double sibling_lower_sum = 0.0;
  double sibling_upper_sum = 0.0;
  for( const std::pair<double,double> &bounds : sibling_bounds )
  {
    sibling_lower_sum += bounds.first;
    sibling_upper_sum += bounds.second;
  }
  return {
    (std::max)(component_lower,(std::max)(0.0,1.0-sibling_upper_sum)),
    (std::min)(component_upper,(std::min)(1.0,1.0-sibling_lower_sum))
  };
}

/** Whether a symmetric local-Gaussian band covers the quantity's actual feasible interval.

 This deliberately accepts the already-projected input/simplex/ratio bounds instead of assuming
 every mass fraction can range over [0,1].  It is the automatic-profile trigger, not an interval
 reported to users. */
inline bool local_gaussian_band_spans_feasible_range(
    const double fraction, const double one_sigma,
    const double feasible_lower, const double feasible_upper,
    const double gaussian_scale = 1.959963984540054 )
{
  if( !std::isfinite(fraction) || !std::isfinite(one_sigma)
      || !std::isfinite(feasible_lower) || !std::isfinite(feasible_upper)
      || !std::isfinite(gaussian_scale) || (one_sigma < 0.0)
      || (gaussian_scale < 0.0) || (feasible_lower > feasible_upper) )
    return false;
  const double half_width = gaussian_scale*one_sigma;
  return ((fraction-half_width) <= feasible_lower)
      && ((fraction+half_width) >= feasible_upper);
}


/** A fixed-composition activity-ratio group contributing to one element's mass denominator.

 `root_activity` is the sole free activity for the group.  Every member's fixed activity ratio and
 activity-per-gram factor has already been absorbed into the two mass coefficients. */
struct ActivityBoxMassGroup
{
  double total_mass_per_root_activity = 0.0;
  double target_mass_per_root_activity = 0.0;
  double lower_root_activity = 0.0;
  double upper_root_activity = std::numeric_limits<double>::infinity();
};


/** Exact target-fraction projection for independent activity boxes and fixed ratio groups.

 A nuclide belongs to exactly one ratio group, so precisely one group can contain the selected
 target.  Its fraction increases monotonically with that root activity and decreases monotonically
 with every competing root.  The extrema therefore use the corresponding activity-box corners;
 unbounded boxes are handled as limits.  `nullopt` denotes malformed or infeasible groups. */
inline std::optional<std::pair<double,double>> activity_box_fraction_bounds(
    const std::vector<ActivityBoxMassGroup> &groups )
{
  if( groups.empty() )
    return std::nullopt;

  size_t target_group = groups.size();
  for( size_t index = 0; index < groups.size(); ++index )
  {
    const ActivityBoxMassGroup &group = groups[index];
    if( !std::isfinite(group.total_mass_per_root_activity)
        || !(group.total_mass_per_root_activity > 0.0)
        || !std::isfinite(group.target_mass_per_root_activity)
        || (group.target_mass_per_root_activity < 0.0)
        || !std::isfinite(group.lower_root_activity)
        || (group.lower_root_activity < 0.0)
        || std::isnan(group.upper_root_activity)
        || (group.upper_root_activity < group.lower_root_activity) )
      return std::nullopt;
    if( group.target_mass_per_root_activity > 0.0 )
    {
      if( target_group != groups.size() )
        return std::nullopt;
      target_group = index;
    }
  }
  if( target_group == groups.size() )
    return std::nullopt;

  const ActivityBoxMassGroup &selected = groups[target_group];
  const double within_group_fraction = selected.target_mass_per_root_activity
                                     / selected.total_mass_per_root_activity;

  // A target root fixed at zero contributes no mass.  It has fraction zero whenever any competing
  // root can be positive; if every root is fixed at zero the normalized composition is infeasible.
  if( selected.upper_root_activity == 0.0 )
  {
    const bool competitor_can_be_positive = std::any_of(
        begin(groups),end(groups),[&selected]( const ActivityBoxMassGroup &group ){
          return (&group != &selected) && (group.upper_root_activity > 0.0);
        } );
    if( !competitor_can_be_positive )
      return std::nullopt;
    return std::pair<double,double>{0.0,0.0};
  }

  double competitor_max_mass = 0.0;
  double competitor_min_mass = 0.0;
  bool competitor_max_unbounded = false;
  for( size_t index = 0; index < groups.size(); ++index )
  {
    if( index == target_group )
      continue;
    const ActivityBoxMassGroup &group = groups[index];
    competitor_min_mass += group.total_mass_per_root_activity*group.lower_root_activity;
    if( std::isinf(group.upper_root_activity) )
      competitor_max_unbounded = true;
    else
      competitor_max_mass += group.total_mass_per_root_activity*group.upper_root_activity;
  }

  double lower = 0.0;
  if( !competitor_max_unbounded && (competitor_max_mass == 0.0) )
  {
    // Identically-zero competitors do not create a zero-fraction limit.  Every positive-total
    // feasible point is wholly in the selected fixed-ratio group.
    lower = within_group_fraction;
  }else if( !competitor_max_unbounded && (selected.lower_root_activity > 0.0) )
  {
    const double target_mass
        = selected.total_mass_per_root_activity*selected.lower_root_activity;
    lower = within_group_fraction*target_mass/(target_mass+competitor_max_mass);
  }

  double upper = 0.0;
  if( std::isinf(selected.upper_root_activity) )
  {
    upper = within_group_fraction;
  }else if( selected.upper_root_activity > 0.0 )
  {
    const double target_mass
        = selected.total_mass_per_root_activity*selected.upper_root_activity;
    upper = within_group_fraction*target_mass/(target_mass+competitor_min_mass);
  }

  if( !std::isfinite(lower) || !std::isfinite(upper) || (lower > upper) )
    return std::nullopt;
  return std::pair<double,double>{
      (std::max)(0.0,(std::min)(1.0,lower)),
      (std::max)(0.0,(std::min)(1.0,upper)) };
}

/** Minimal rank information for a profile point which improved the retained baseline.

 The physical objective is deliberately kept separate from the profile's reported coordinate and
 scan order.  `semantic_key` must identify the curve and target independently of caller order.
 */
struct PendingBaselineDiscovery
{
  double full_objective = std::numeric_limits<double>::infinity();
  std::string semantic_key;
};


/** Policy for a conditional point which improves the retained physical objective.

 The first profile pass must defer every such point until all eligible targets have been visited;
 only then is the deterministic best seed selected.  Profiles restarted from that selected
 baseline are the one permitted retry, so another discovery is an explicit baseline-not-optimal
 failure rather than a second restart.

 Note on the conditional-fit budget: the cap passed to `scan` is per scan on a given baseline, so a
 baseline reselection deliberately grants the same semantic target a fresh budget rather than
 carrying the first pass's spend forward.  Sharing one budget across both baselines would starve the
 restarted scan and turn a legitimately bounded profile into a `FitCapReached` failure, which is
 reported as `Failed`; the restarted scan is a scan of a *different* physical baseline and is
 budgeted as such.
 */
enum class BaselineDiscoveryDisposition
{
  DeferUntilPassComplete,
  RejectAfterRestart
};

inline BaselineDiscoveryDisposition baseline_discovery_disposition(
    const unsigned completed_baseline_restarts )
{
  return completed_baseline_restarts == 0
      ? BaselineDiscoveryDisposition::DeferUntilPassComplete
      : BaselineDiscoveryDisposition::RejectAfterRestart;
}


/** The smallest objective decrease that may be called a genuinely better baseline.

 A conditional fit which lands microscopically below the retained solution has not found another
 basin; it has resampled the same one.  Two independent scales bound how small "microscopic" is:

 - **What the optimizer can deliver.**  These fits are routinely ill-conditioned - the 610-775 keV
   free-age plutonium problem solves at a condition number near 1e7 with a genuinely rank-deficient
   direction - so the attainable objective accuracy is far coarser than any solver tolerance
   suggests.  `solve_ceres` documents the measured figure beside `parameter_tolerance`: about 0.1 in
   chi-square for an example physical model.  A gate below that fires on round-off, not on physics.
 - **What the result can express.**  Interval endpoints are crossings of `cov_scale*1` and
   `cov_scale*3.841458820694124`.  An improvement that is a thousandth of the 68.27% threshold
   cannot move a reported endpoint by anything a reader could detect.

 Chasing a decrease below both scales is worse than useless: along a numerically flat direction the
 lowest sampled point is chosen by round-off, so re-seeding to it makes the reported nominal
 arbitrary rather than more accurate.  The honest report of a flat direction is the bounded or
 non-identifiable interval that spans it, which is exactly what this scanner produces.

 The guard this feeds is unchanged - one deferred reselection, then a visible failure - and a real
 alternative basin clears these scales by orders of magnitude.
 */
inline double baseline_improvement_tolerance( const double objective, const double cov_scale )
{
  const double magnitude = std::isfinite(objective) ? std::fabs(objective) : 0.0;
  const double scale = std::isfinite(cov_scale) ? std::fabs(cov_scale) : 0.0;
  return (std::max)( {1.0e-6, 1.0e-5*(1.0 + magnitude), 1.0e-3*scale} );
}


/** Select the deterministic seed for a single deferred baseline reselection.

 Only finite full-objective values are eligible.  A strictly lower objective wins; an exactly
 equal objective is broken by the lexicographically smaller semantic key.  If both fields are
 equal, the earlier entry is retained because the entries are semantically indistinguishable.
 `std::nullopt` means that no finite discovery was supplied.
 */
inline std::optional<std::size_t> best_pending_baseline_discovery_index(
    const std::vector<PendingBaselineDiscovery> &discoveries )
{
  std::optional<std::size_t> best;
  for( std::size_t index = 0; index < discoveries.size(); ++index )
  {
    const PendingBaselineDiscovery &candidate = discoveries[index];
    if( !std::isfinite(candidate.full_objective) )
      continue;
    if( !best )
    {
      best = index;
      continue;
    }

    const PendingBaselineDiscovery &incumbent = discoveries[*best];
    if( (candidate.full_objective < incumbent.full_objective)
        || ((candidate.full_objective == incumbent.full_objective)
            && (candidate.semantic_key < incumbent.semantic_key)) )
      best = index;
  }
  return best;
}


/** Order every finite deferred discovery best-first using exactly the rule of
 `best_pending_baseline_discovery_index`.

 One warm seed cannot cover basins that different profile targets entered independently, so the
 single permitted baseline reselection is given several deterministic seeds.  Non-finite entries are
 dropped rather than sorted to the end: they can never seed anything.  Entries whose objective and
 semantic key are both equal keep their original relative order, so the result is a total order that
 never depends on caller/source ordering.

 The returned indices refer to `discoveries`.  `front()` is always the same entry
 `best_pending_baseline_discovery_index` would return, which keeps the one-seed and many-seed paths
 provably consistent.
 */
inline std::vector<std::size_t> ordered_pending_baseline_discoveries(
    const std::vector<PendingBaselineDiscovery> &discoveries )
{
  std::vector<std::size_t> order;
  order.reserve(discoveries.size());
  for( std::size_t index = 0; index < discoveries.size(); ++index )
    if( std::isfinite(discoveries[index].full_objective) )
      order.push_back(index);

  std::stable_sort( order.begin(), order.end(),
    [&discoveries]( const std::size_t lhs, const std::size_t rhs ) {
      const PendingBaselineDiscovery &a = discoveries[lhs], &b = discoveries[rhs];
      if( a.full_objective < b.full_objective )
        return true;
      if( b.full_objective < a.full_objective )
        return false;
      return a.semantic_key < b.semantic_key;
    } );
  return order;
}


/** De-duplicate an ordered discovery list by semantic key, keeping the best entry per key.

 Two profile targets in the same curve routinely fall into the same basin, and seeding the search
 twice from it wastes a candidate slot that a genuinely different basin could use.  Input must
 already be ordered by `ordered_pending_baseline_discoveries`, so the retained entry per key is that
 key's lowest objective.
 */
inline std::vector<std::size_t> unique_pending_baseline_discoveries(
    const std::vector<PendingBaselineDiscovery> &discoveries,
    const std::vector<std::size_t> &order )
{
  std::vector<std::size_t> unique;
  unique.reserve(order.size());
  for( const std::size_t index : order )
  {
    bool duplicate = false;
    for( const std::size_t retained : unique )
      duplicate = duplicate
                || (discoveries[index].semantic_key == discoveries[retained].semantic_key);
    if( !duplicate )
      unique.push_back(index);
  }
  return unique;
}

namespace detail
{
inline bool finite_sample( const Evaluation &evaluation )
{
  return std::isfinite(evaluation.reported_fraction)
         && std::isfinite(evaluation.delta_chi2);
}

inline bool at_control_bound( const Sample &sample, const double bound,
                              const double full_span )
{
  const double tolerance = (std::max)(1.0e-13, 1.0e-11*full_span);
  return std::fabs(sample.control_fraction - bound) <= tolerance;
}


}//namespace detail


}//namespace ProfileLikelihood
}//namespace RelActCalcAutoImp

#endif //InterSpec_RelActCalcAuto_ProfileScan_h
