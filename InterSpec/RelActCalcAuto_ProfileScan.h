/*
  Deterministic one-dimensional profile-likelihood scanning utilities.

  This header deliberately has no dependency on the RelActAuto solver.  The production profile
  code supplies conditional-fit evaluations, while unit tests can supply analytic objectives.  In
  both cases the same bracketing, fit-budget, endpoint, and non-identifiability rules are used.
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
  /** Number of actual conditional optimizations used to produce this evaluation. */
  std::size_t num_fits = 1;
  /** Absolute accuracy with which the conditional solve reached its requested control coordinate.
   Brent must not refine below this scale when the objective was evaluated at a nearby reported
   coordinate. */
  double control_tolerance = 0.0;
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
  double control_tolerance = 0.0;
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


/** Coordinate used for a conditional solve at a requested profile point.

 Every interior coordinate is preserved exactly, including sub-ppm input limits.  Only the open
 production parameterizations at the mathematical zero/one limits need an inward representative;
 its displacement scales with the actual feasible span instead of imposing an absolute 2e-6 floor. */
inline double conditional_solve_control( const double requested,
                                         const double baseline,
                                         const double feasible_lower,
                                         const double feasible_upper,
                                         const double inward_fraction = 2.0e-6,
                                         const double minimum_upper_inset = 0.0 )
{
  if( !std::isfinite(requested) || !std::isfinite(baseline)
      || !std::isfinite(feasible_lower)
      || !std::isfinite(feasible_upper) || !(feasible_lower <= feasible_upper)
      || !(inward_fraction > 0.0) || (minimum_upper_inset < 0.0) )
    return requested;

  const double span = feasible_upper-feasible_lower;
  if( !(span > 0.0) )
    return requested;
  if( (requested == 0.0) && (feasible_lower == 0.0) )
  {
    // Scale against the distance to the baseline, not the whole feasible range: for a trace
    // fraction near 1e-9 in [0,1], a fixed 2e-6 inset crosses to the wrong side of the optimum.
    const double directional_span = (std::max)(0.0,baseline-feasible_lower);
    if( !(directional_span > 0.0) )
      return requested;
    const double inward = feasible_lower + inward_fraction*directional_span;
    return (std::min)(feasible_upper,
        (std::max)(std::nextafter(feasible_lower,feasible_upper),inward));
  }
  if( (requested == 1.0) && (feasible_upper == 1.0) )
  {
    const double directional_span = (std::max)(0.0,feasible_upper-baseline);
    if( !(directional_span > 0.0) )
      return requested;
    const double inset = (std::max)(minimum_upper_inset,
                                    inward_fraction*directional_span);
    const double inward = feasible_upper-inset;
    return (std::max)(feasible_lower,
        (std::min)(std::nextafter(feasible_upper,feasible_lower),inward));
  }
  return requested;
}


/** Scale-aware reported-coordinate equality tolerance for an AL conditional solve. */
inline double conditional_equality_tolerance( const double requested,
                                              const double baseline,
                                              const double feasible_lower,
                                              const double feasible_upper )
{
  if( !std::isfinite(requested) || !std::isfinite(baseline)
      || !std::isfinite(feasible_lower)
      || !std::isfinite(feasible_upper) || !(feasible_lower <= feasible_upper) )
    return std::numeric_limits<double>::quiet_NaN();
  const bool lower_direction = (requested <= baseline);
  const double directional_span = lower_direction
      ? (std::max)(0.0,baseline-feasible_lower)
      : (std::max)(0.0,feasible_upper-baseline);
  const double scale = (std::max)(directional_span,
                                  std::numeric_limits<double>::min());
  const double roundoff_floor = 64.0*std::numeric_limits<double>::epsilon()
                              * (std::max)(1.0,std::fabs(requested));
  // Near a trace baseline, a conditional solve can be accurate to roughly a per-mille of the
  // requested coordinate while still missing an absolute tolerance derived from the entire
  // direction by orders of magnitude.  Since the scanner records the independently evaluated
  // coordinate (rather than pretending it hit the request exactly), accepting that precision is
  // both honest and ample for endpoint refinement.  Cap this allowance at 1e-8 so ordinary
  // fractions retain the tighter directional tolerance, and scale from the smaller of the
  // baseline/request so mathematical-zero surrogates are not mistaken for a nearby trace value.
  const double trace_coordinate_tolerance = (std::min)(
      1.0e-8,1.0e-3*(std::min)(std::fabs(requested),std::fabs(baseline)) );
  return (std::max)({roundoff_floor,2.0e-5*scale,trace_coordinate_tolerance});
}


/** Whether repeated reported-coordinate AL fits have reached an interior structural limit. */
inline bool stable_reported_structural_bound( const double requested,
                                              const double actual,
                                              const double baseline,
                                              const double last_reported_change,
                                              const double equality_tolerance,
                                              const double directional_span,
                                              const std::size_t num_fits,
                                              const Direction direction )
{
  if( !std::isfinite(requested) || !std::isfinite(actual) || !std::isfinite(baseline)
      || !std::isfinite(last_reported_change) || !std::isfinite(equality_tolerance)
      || !std::isfinite(directional_span) || (equality_tolerance < 0.0)
      || (directional_span < 0.0) || (num_fits < 2) )
    return false;
  const bool on_requested_side = (direction == Direction::Lower)
      ? ((actual >= requested) && (actual <= baseline))
      : ((actual <= requested) && (actual >= baseline));
  const double stability_tolerance = (std::max)(
      5.0*equality_tolerance,2.0e-6*directional_span );
  return on_requested_side
      && (std::fabs(actual-requested) > equality_tolerance)
      && (last_reported_change <= stability_tolerance);
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

/** Return the stable index of the cached control point nearest `requested_control`.

 Non-finite cache entries are ignored and ties retain the earlier entry.  `controls.size()` means
 no usable cache point.  Production uses this for semantic conditional-solve warm starts; exposing
 the tiny policy here keeps the non-chronological Brent/68-to-95 traversal independently testable.
 */
inline std::size_t nearest_control_index( const std::vector<double> &controls,
                                         const double requested_control )
{
  std::size_t nearest = controls.size();
  double nearest_distance = std::numeric_limits<double>::infinity();
  if( !std::isfinite(requested_control) )
    return nearest;
  for( std::size_t index = 0; index < controls.size(); ++index )
  {
    if( !std::isfinite(controls[index]) )
      continue;
    const double distance = std::fabs(controls[index]-requested_control);
    if( distance < nearest_distance )
    {
      nearest = index;
      nearest_distance = distance;
    }
  }
  return nearest;
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
         && std::isfinite(evaluation.delta_chi2)
         && std::isfinite(evaluation.control_tolerance)
         && (evaluation.control_tolerance >= 0.0);
}

inline bool at_control_bound( const Sample &sample, const double bound,
                              const double full_span )
{
  const double tolerance = (std::max)(1.0e-13, 1.0e-11*full_span);
  return std::fabs(sample.control_fraction - bound) <= tolerance;
}


/** Select an independently evaluated crossing endpoint at the true global fit cap.

 This is deliberately more restrictive than ordinary Brent convergence: it is available only for
 the final bracketed crossing after the complete optimization budget has been spent.  A bracket is
 terminally usable when its control width is at most ten ordinary control tolerances (capped at
 1e-4 of the directional span), or when one endpoint is already within 0.5% of the requested
 likelihood threshold.  The evaluated outside endpoint is returned so the reported confidence
 interval never truncates parameter values which are independently known to remain below the
 threshold.  A broad, poorly localized bracket remains a fit-cap failure.
 */
inline std::optional<Sample> terminal_cap_crossing_sample(
    const Sample &inside, const Sample &outside, const double threshold,
    const double normal_control_tolerance, const double directional_span,
    const bool global_fit_cap_exhausted, const bool final_bracketed_crossing )
{
  if( !global_fit_cap_exhausted || !final_bracketed_crossing
      || !std::isfinite(threshold) || !(threshold > 0.0)
      || !std::isfinite(normal_control_tolerance) || !(normal_control_tolerance > 0.0)
      || !std::isfinite(directional_span) || !(directional_span > 0.0)
      || !std::isfinite(inside.control_fraction)
      || !std::isfinite(outside.control_fraction)
      || !std::isfinite(inside.reported_fraction)
      || !std::isfinite(outside.reported_fraction)
      || !std::isfinite(inside.delta_chi2) || !std::isfinite(outside.delta_chi2)
      || !(inside.delta_chi2 < threshold) || !(outside.delta_chi2 >= threshold) )
    return std::nullopt;

  const double width = std::fabs(outside.control_fraction-inside.control_fraction);
  const double roundoff_floor = 64.0*std::numeric_limits<double>::epsilon()
      *(std::max)({1.0,std::fabs(inside.control_fraction),
                   std::fabs(outside.control_fraction)});
  const double tight_width_limit = (std::max)(roundoff_floor,
      (std::min)(10.0*normal_control_tolerance,1.0e-4*directional_span));
  const double outside_residual = std::fabs(outside.delta_chi2-threshold);
  const double terminal_objective_tolerance = (std::max)(1.0e-5,5.0e-3*threshold);
  if( (width > tight_width_limit)
      && (outside_residual > terminal_objective_tolerance) )
    return std::nullopt;

  return outside;
}

}//namespace detail


/** Scan a bounded scalar profile at two nested likelihood thresholds.

 The scan brackets geometrically at 1/16, 1/8, ..., 1 of each available side.  Crossings are refined
 with the bracket-preserving Brent-Dekker algorithm: inverse quadratic interpolation when possible,
 a secant step with two distinct values, and bisection whenever interpolation is unsafe.  Every
 crossing shares one global conditional-optimization budget.  Successful refinement samples
 are retained for the nested threshold, and no callback may take the total beyond `max_fits`.

 `evaluate` must return the objective difference from the same independently evaluated baseline.
 It may map the control coordinate to a different reported coordinate, as is needed for Pu-242
 correlation/renormalization.
 */
inline ScanResult scan( const double baseline_control,
                        const double baseline_reported,
                        const double lower_control,
                        const double upper_control,
                        const std::array<double,2> thresholds,
                        const std::size_t max_fits,
                        const Evaluator &evaluate )
{
  ScanResult answer;
  answer.intervals[0].delta_chi2 = thresholds[0];
  answer.intervals[1].delta_chi2 = thresholds[1];

  if( !evaluate || !std::isfinite(baseline_control) || !std::isfinite(baseline_reported)
      || !std::isfinite(lower_control) || !std::isfinite(upper_control)
      || (lower_control > baseline_control) || (baseline_control > upper_control)
      || !std::isfinite(thresholds[0]) || !std::isfinite(thresholds[1])
      || (max_fits == 0)
      || !(thresholds[0] > 0.0) || !(thresholds[1] > thresholds[0]) )
  {
    answer.status = ScanStatus::Failed;
    answer.diagnostic = "Invalid profile bounds, baseline, thresholds, or evaluator.";
    return answer;
  }

  const Sample baseline{baseline_control,baseline_reported,0.0,false};
  answer.samples.push_back(baseline);
  std::vector<Sample> lower_samples{baseline};
  std::vector<Sample> upper_samples{baseline};
  const double full_span = upper_control - lower_control;
  const bool have_lower_span = baseline_control > lower_control + 1.0e-14;
  const bool have_upper_span = baseline_control < upper_control - 1.0e-14;
  bool lower_done = !have_lower_span;
  bool upper_done = !have_upper_span;
  bool used_terminal_cap_crossing = false;

  auto request = [&]( const double control, const Direction direction,
                      Sample &sample, const std::size_t local_fit_cap ) -> bool {
    if( answer.num_fits >= max_fits )
    {
      answer.status = ScanStatus::FitCapReached;
      answer.diagnostic = "The conditional-fit cap was reached before the profile was classified.";
      return false;
    }
    const std::size_t remaining_fits
        = (std::min)(max_fits - answer.num_fits,local_fit_cap);
    if( remaining_fits == 0 )
    {
      answer.status = ScanStatus::FitCapReached;
      answer.diagnostic = "The conditional-fit allocation was exhausted before a profile point.";
      return false;
    }
    const Evaluation evaluation = evaluate(control,direction,remaining_fits);
    if( evaluation.num_fits > remaining_fits )
    {
      answer.status = ScanStatus::FitCapReached;
      answer.num_fits = max_fits;
      answer.diagnostic = "A conditional evaluation exceeded the remaining profile-fit budget.";
      return false;
    }
    answer.num_fits += evaluation.num_fits;
    switch( evaluation.status )
    {
      case EvaluationStatus::Canceled:
        answer.status = ScanStatus::Canceled;
        answer.diagnostic = evaluation.diagnostic.empty()
                          ? "Mass-fraction profiling was canceled." : evaluation.diagnostic;
        return false;
      case EvaluationStatus::BetterBaseline:
        answer.status = ScanStatus::BetterBaseline;
        answer.diagnostic = evaluation.diagnostic.empty()
                          ? "A conditional fit found a better baseline." : evaluation.diagnostic;
        return false;
      case EvaluationStatus::FitCapReached:
        if( answer.num_fits >= max_fits )
        {
          answer.status = ScanStatus::FitCapReached;
          answer.diagnostic = evaluation.diagnostic.empty()
                            ? "The conditional-fit cap was reached." : evaluation.diagnostic;
        }else
        {
          answer.status = ScanStatus::Failed;
          answer.diagnostic = evaluation.diagnostic.empty()
              ? "A conditional evaluator reported a fit cap before exhausting the global budget."
              : evaluation.diagnostic;
        }
        return false;
      case EvaluationStatus::Failed:
        answer.status = ScanStatus::Failed;
        answer.diagnostic = evaluation.diagnostic.empty()
                          ? "A conditional profile fit failed." : evaluation.diagnostic;
        return false;
      case EvaluationStatus::Success:
        break;
    }
    if( !detail::finite_sample(evaluation) )
    {
      answer.status = ScanStatus::Failed;
      answer.diagnostic = "A conditional profile fit returned a non-finite result.";
      return false;
    }
    sample.control_fraction = evaluation.reached_feasible_bound
                            ? evaluation.reported_fraction : control;
    sample.reported_fraction = evaluation.reported_fraction;
    sample.delta_chi2 = (std::max)(0.0,evaluation.delta_chi2);
    sample.reached_feasible_bound = evaluation.reached_feasible_bound;
    sample.control_tolerance = evaluation.control_tolerance;
    answer.samples.push_back(sample);
    return true;
  };

  // Geometric bracketing, interleaved so a small cap cannot silently classify just one side.
  for( int exponent = -4; exponent <= 0 && (!lower_done || !upper_done); ++exponent )
  {
    const double outward_fraction = std::ldexp(1.0,exponent);
    if( !lower_done )
    {
      const double control = (exponent == 0) ? lower_control
                         : baseline_control
                           - outward_fraction*(baseline_control-lower_control);
      Sample sample;
      if( !request(control,Direction::Lower,sample,max_fits) )
        return answer;
      if( sample.reached_feasible_bound
          && (sample.control_fraction > lower_samples.back().control_fraction + 1.0e-12) )
      {
        answer.status = ScanStatus::Failed;
        answer.diagnostic = "A reported lower structural bound was not outward from the baseline.";
        return answer;
      }
      lower_samples.push_back(sample);
      lower_done = (sample.delta_chi2 >= thresholds[1])
                   || sample.reached_feasible_bound
                   || detail::at_control_bound(sample,lower_control,full_span);
    }
    if( !upper_done )
    {
      const double control = (exponent == 0) ? upper_control
                         : baseline_control
                           + outward_fraction*(upper_control-baseline_control);
      Sample sample;
      if( !request(control,Direction::Upper,sample,max_fits) )
        return answer;
      if( sample.reached_feasible_bound
          && (sample.control_fraction < upper_samples.back().control_fraction - 1.0e-12) )
      {
        answer.status = ScanStatus::Failed;
        answer.diagnostic = "A reported upper structural bound was not outward from the baseline.";
        return answer;
      }
      upper_samples.push_back(sample);
      upper_done = (sample.delta_chi2 >= thresholds[1])
                   || sample.reached_feasible_bound
                   || detail::at_control_bound(sample,upper_control,full_span);
    }
  }

  if( !lower_done || !upper_done )
  {
    answer.status = ScanStatus::Failed;
    answer.diagnostic = "Geometric profile bracketing ended before reaching a crossing or bound.";
    return answer;
  }

  const auto crossing_index = []( const std::vector<Sample> &side, const double threshold ) {
    for( std::size_t i = 1; i < side.size(); ++i )
      if( side[i].delta_chi2 >= threshold )
        return i;
    return side.size();
  };
  std::size_t crossing_endpoints_remaining = 0;
  for( const double threshold : thresholds )
  {
    crossing_endpoints_remaining += have_lower_span
          && (crossing_index(lower_samples,threshold) != lower_samples.size());
    crossing_endpoints_remaining += have_upper_span
          && (crossing_index(upper_samples,threshold) != upper_samples.size());
  }

  auto refine_endpoint = [&]( std::vector<Sample> &side,
                              const bool have_span,
                              const double bound,
                              const Direction direction,
                              const double threshold,
                              Endpoint &endpoint ) -> bool {
    if( !have_span )
    {
      endpoint.reported_fraction = baseline_reported;
      endpoint.likelihood_crossing = false;
      return true;
    }

    const std::size_t crossing = crossing_index(side,threshold);

    if( crossing == side.size() )
    {
      if( !side.back().reached_feasible_bound
          && !detail::at_control_bound(side.back(),bound,full_span) )
      {
        answer.status = ScanStatus::Failed;
        answer.diagnostic = "A likelihood direction reached neither its threshold nor its bound.";
        return false;
      }
      endpoint.reported_fraction = side.back().reported_fraction;
      endpoint.likelihood_crossing = false;
      return true;
    }

    Sample inside = side[crossing-1];
    Sample outside = side[crossing];

    // All endpoints share one hard global budget.  A formerly equal local split could stop a
    // difficult early crossing while many fits remained reserved for endpoints which were never
    // attempted.  In particular, a callback may consume up to three augmented-Lagrangian solves,
    // so a seven-fit local share was only two useful profile points.  Let a crossing borrow unused
    // global fits; successful samples are inserted into `side` below and normally make the nested
    // threshold substantially cheaper.
    if( crossing_endpoints_remaining == 0 )
    {
      answer.status = ScanStatus::Failed;
      answer.diagnostic = "Internal profile crossing-count inconsistency.";
      return false;
    }
    --crossing_endpoints_remaining;

    const auto accept_terminal_cap_crossing = [&]() {
      const double directional_span = std::fabs(bound-baseline_control);
      const double normal_control_tolerance = (std::max)({
          64.0*std::numeric_limits<double>::epsilon()
            *(std::max)({1.0,std::fabs(inside.control_fraction),
                         std::fabs(outside.control_fraction)}),
          1.0e-5*directional_span,
          inside.control_tolerance,
          outside.control_tolerance});
      const std::optional<Sample> terminal_sample = detail::terminal_cap_crossing_sample(
          inside,outside,threshold,normal_control_tolerance,directional_span,
          (answer.num_fits >= max_fits),(crossing_endpoints_remaining == 0));
      if( !terminal_sample )
        return false;
      endpoint.reported_fraction = terminal_sample->reported_fraction;
      endpoint.likelihood_crossing = true;
      used_terminal_cap_crossing = true;
      return true;
    };

    if( answer.num_fits >= max_fits )
    {
      // The outer threshold may already have tightened a shared bracket before consuming the
      // final fit.  Let the nested final endpoint use that independently evaluated bracket under
      // the same conservative terminal-cap rule instead of failing without entering Brent.
      if( accept_terminal_cap_crossing() )
        return true;
      answer.status = ScanStatus::FitCapReached;
      answer.diagnostic = "The conditional-fit cap left no budget for a bracketed crossing.";
      return false;
    }

    // Classic van Wijngaarden-Dekker-Brent root iteration.  `a` is the previous best point, `b`
    // the current best, and `c` the opposite-sign bracket endpoint.  The interpolation path uses
    // inverse quadratic interpolation when three distinct function values are available, secant
    // otherwise; Brent's acceptance inequalities fall back to bisection when interpolation is not
    // safely inside the shrinking bracket.
    Sample sample_a = inside;
    Sample sample_b = outside;
    Sample sample_c = outside;
    double a = sample_a.control_fraction;
    double b = sample_b.control_fraction;
    double c = sample_c.control_fraction;
    double fa = sample_a.delta_chi2-threshold;
    double fb = sample_b.delta_chi2-threshold;
    double fc = fb;
    double d = b-a;
    double e = d;
    bool converged = false;

    for( std::size_t iteration = 0; iteration < 64; ++iteration )
    {
      const double objective_tolerance = (std::max)(1.0e-5,5.0e-4*threshold);
      // `inside` and `outside` are the tightest independently evaluated opposite-sign samples.
      // Brent's historical b/c bookkeeping can temporarily retain a wider pair after an
      // interpolation swap, so honor the tighter bracket directly before proposing another
      // expensive conditional solve.
      const double tight_bracket_width
          = std::fabs(outside.control_fraction-inside.control_fraction);
      const double tight_control_tolerance = (std::max)({
          64.0*std::numeric_limits<double>::epsilon()
            *(std::max)({1.0,std::fabs(inside.control_fraction),
                         std::fabs(outside.control_fraction)}),
          1.0e-5*std::fabs(bound-baseline_control),
          inside.control_tolerance,
          outside.control_tolerance});
      const double tight_inside_residual = std::fabs(inside.delta_chi2-threshold);
      const double tight_outside_residual = std::fabs(outside.delta_chi2-threshold);
      if( (tight_bracket_width <= tight_control_tolerance)
          || ((std::min)(tight_inside_residual,tight_outside_residual)
                <= objective_tolerance) )
      {
        sample_b = (tight_inside_residual <= tight_outside_residual) ? inside : outside;
        b = sample_b.control_fraction;
        fb = sample_b.delta_chi2-threshold;
        converged = true;
        break;
      }

      if( ((fb > 0.0) && (fc > 0.0)) || ((fb < 0.0) && (fc < 0.0)) )
      {
        c = a;
        sample_c = sample_a;
        fc = fa;
        d = e = b-a;
      }
      if( std::fabs(fc) < std::fabs(fb) )
      {
        const double old_b = b, old_fb = fb;
        const Sample old_sample_b = sample_b;
        a = b; sample_a = sample_b; fa = fb;
        b = c; sample_b = sample_c; fb = fc;
        c = old_b; sample_c = old_sample_b; fc = old_fb;
      }

      // Scale root accuracy to this direction.  A trace fraction near zero can have a complete
      // lower bracket narrower than 1e-5 of [0,1]; using the full feasible span would incorrectly
      // accept the baseline itself as both lower crossings.
      const double directional_span = std::fabs(bound-baseline_control);
      const double control_tolerance = (std::max)({
          64.0*std::numeric_limits<double>::epsilon()
            *(std::max)(1.0,std::fabs(b)),
          1.0e-5*directional_span,
          sample_b.control_tolerance,
          sample_c.control_tolerance});
      const double tolerance = 2.0*std::numeric_limits<double>::epsilon()*std::fabs(b)
                               + 0.5*control_tolerance;
      const double midpoint = 0.5*(c-b);
      if( (std::fabs(midpoint) <= tolerance) || (std::fabs(fb) <= objective_tolerance) )
      {
        converged = true;
        break;
      }

      if( answer.num_fits >= max_fits )
        break;

      if( (std::fabs(e) >= tolerance) && (std::fabs(fa) > std::fabs(fb)) )
      {
        const double s = fb/fa;
        double p = 0.0, q = 0.0;
        if( a == c )
        {
          p = 2.0*midpoint*s;
          q = 1.0-s;
        }else
        {
          q = fa/fc;
          const double r = fb/fc;
          p = s*(2.0*midpoint*q*(q-r) - (b-a)*(r-1.0));
          q = (q-1.0)*(r-1.0)*(s-1.0);
        }
        if( p > 0.0 )
          q = -q;
        p = std::fabs(p);
        const double interpolation_limit
            = (std::min)(3.0*midpoint*q-std::fabs(tolerance*q),std::fabs(e*q));
        if( (q != 0.0) && (2.0*p < interpolation_limit) )
        {
          e = d;
          d = p/q;
        }else
        {
          d = midpoint;
          e = d;
        }
      }else
      {
        d = midpoint;
        e = d;
      }

      a = b;
      sample_a = sample_b;
      fa = fb;
      const double step = (std::fabs(d) > tolerance)
                        ? d : std::copysign(tolerance,midpoint);
      double next_control = b+step;
      const double tight_low
          = (std::min)(inside.control_fraction,outside.control_fraction);
      const double tight_high
          = (std::max)(inside.control_fraction,outside.control_fraction);
      // On a flat profile followed by a very steep rise, a secant step can hug the flat endpoint;
      // classic Brent then spends the next call bisecting and repeats this two-call cycle.  Reject
      // interpolation which contracts less than 10% of the known sign bracket and bisect now.
      // This is the same progress safeguard Brent applies to its internal history, but expressed
      // against the tighter physical bracket retained above.
      constexpr double minimum_bracket_contraction = 0.10;
      const double minimum_progress = minimum_bracket_contraction*(tight_high-tight_low);
      const bool endpoint_hugging = (next_control <= tight_low+minimum_progress)
                                 || (next_control >= tight_high-minimum_progress);
      const double smaller_residual
          = (std::min)(tight_inside_residual,tight_outside_residual);
      const double larger_residual
          = (std::max)(tight_inside_residual,tight_outside_residual);
      const bool flat_to_steep_imbalance
          = (smaller_residual > 0.05*threshold)
            && (larger_residual > 100.0*smaller_residual);
      if( !std::isfinite(next_control)
          || (endpoint_hugging && flat_to_steep_imbalance) )
        next_control = 0.5*(tight_low+tight_high);
      Sample trial;
      // `request` enforces the sole hard global cap.  Passing the true remainder lets the callback
      // complete its bounded augmented-Lagrangian updates instead of reporting a local allocation
      // failure while global fits remain.
      if( !request(next_control,direction,trial,max_fits-answer.num_fits) )
        return false;
      b = trial.control_fraction;
      sample_b = trial;
      fb = trial.delta_chi2-threshold;
      if( fb < 0.0 )
        inside = trial;
      else
        outside = trial;

      // Retain refinement work for the nested threshold.  Directional sample vectors are ordered
      // from the baseline outwards, just like the original geometric grid.
      side.push_back(trial);
      std::stable_sort(begin(side),end(side),[direction]( const Sample &lhs, const Sample &rhs ){
        return (direction == Direction::Lower)
             ? (lhs.control_fraction > rhs.control_fraction)
             : (lhs.control_fraction < rhs.control_fraction);
      });
    }//Brent iterations

    if( !converged )
    {
      const bool exhausted = (answer.num_fits >= max_fits);
      if( accept_terminal_cap_crossing() )
        return true;
      answer.status = exhausted ? ScanStatus::FitCapReached : ScanStatus::Failed;
      answer.diagnostic = exhausted
          ? "The conditional-fit cap was reached before a Brent crossing converged."
          : "The Brent iteration limit was reached before a profile crossing converged.";
      return false;
    }
    endpoint.reported_fraction = sample_b.reported_fraction;
    endpoint.likelihood_crossing = true;
    return true;
  };

  // Refine the outer threshold first.  Its bracket encloses the inner crossing, so the retained
  // Brent samples tighten the later inner problem; the reverse order cannot similarly narrow the
  // still-outward outer bracket.
  for( const std::size_t level : std::array<std::size_t,2>{{1,0}} )
  {
    if( !refine_endpoint(lower_samples,have_lower_span,lower_control,Direction::Lower,
                         thresholds[level],answer.intervals[level].lower) )
      return answer;
    if( !refine_endpoint(upper_samples,have_upper_span,upper_control,Direction::Upper,
                         thresholds[level],answer.intervals[level].upper) )
      return answer;
    if( answer.intervals[level].lower.reported_fraction
          > answer.intervals[level].upper.reported_fraction )
      std::swap(answer.intervals[level].lower,answer.intervals[level].upper);
  }

  const auto whole_side_inside = [&]( const std::vector<Sample> &side,
                                      const bool have_span, const double bound ) {
    if( !have_span )
      return true;
    if( !side.back().reached_feasible_bound
        && !detail::at_control_bound(side.back(),bound,full_span) )
      return false;
    return std::all_of(side.begin(),side.end(),[&]( const Sample &sample ){
      return sample.delta_chi2 < thresholds[1];
    });
  };
  const bool entire_range_inside = whole_side_inside(lower_samples,have_lower_span,lower_control)
                                && whole_side_inside(upper_samples,have_upper_span,upper_control);
  if( entire_range_inside )
  {
    answer.status = ScanStatus::NonIdentifiable;
    answer.diagnostic = "The entire feasible mass-fraction range remains inside the 95% likelihood threshold.";
  }else
  {
    bool has_boundary = false;
    for( const Interval &interval : answer.intervals )
      has_boundary = has_boundary || !interval.lower.likelihood_crossing
                                  || !interval.upper.likelihood_crossing;
    answer.status = has_boundary ? ScanStatus::BoundaryLimited : ScanStatus::Complete;
    if( used_terminal_cap_crossing )
      answer.diagnostic = "The final likelihood crossing used the independently evaluated outside"
                          " endpoint of a terminal sign bracket at the conditional-fit cap.";
    else
      answer.diagnostic = has_boundary
        ? "At least one confidence endpoint is set by a physical or input bound."
        : "Both likelihood directions were bracketed.";
  }
  return answer;
}

}//namespace ProfileLikelihood
}//namespace RelActCalcAutoImp

#endif //InterSpec_RelActCalcAuto_ProfileScan_h
