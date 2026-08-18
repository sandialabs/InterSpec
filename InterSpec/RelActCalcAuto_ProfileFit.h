/*
  Vertex-anchored profile-likelihood fitting.

  Like `RelActCalcAuto_ProfileScan.h`, this header deliberately has no dependency on the RelActAuto
  solver: the production code supplies conditional-fit evaluations, while unit tests supply analytic
  objectives.  It replaces the Brent root-finder used previously, which needed a function it could
  bracket and could not cope with the several-chi-squared scatter a conditional re-solve produces.

  The model exploits the fact that the vertex is known a priori.  The unconstrained optimum IS the
  minimum of the profile in its own coordinate, so `dchi2(q0) == 0` with zero slope by construction.
  What is left is a one-, two-, or three-parameter model

      dchi2(d) = c2*d^2 + c3*d^3 + c4*d^4,     d = |q - q0|

  fitted independently on each side (the profile is generically asymmetric, and forcing symmetry is
  a real modelling error near a bound).  A single off-optimum point already determines `c2`;
  everything beyond the first point per side exists to TEST the model, not to bracket anything.

  Two rules govern which samples the model may use, and both are load-bearing:

   - Rejection is ONE-SIDED.  A conditional fit is a minimization, so a failed, under-converged, or
     wrong-basin solve can only push dchi2 too HIGH, never too low; the true profile is the lower
     envelope of the samples.  What separates noise from curvature is not where a point sits but
     whether it is ISOLATED - real curvature lifts a point and its neighbours, a failed solve lifts
     one point alone - so isolation is tested directly, by leave-one-out.

   - Reported crossings INTERPOLATE, never extrapolate.  A crossing must lie inside the span of
     evaluated samples.  The vertex is exempt and serves as the inboard anchor: `(0,0)` is not a
     sample that might be wrong, it is known analytically.  So the requirement reduces to needing an
     evaluated sample at or above the threshold outboard of the crossing.
*/
#ifndef InterSpec_RelActCalcAuto_ProfileFit_h
#define InterSpec_RelActCalcAuto_ProfileFit_h

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include <limits>
#include <optional>
#include <algorithm>
#include <functional>

#include "InterSpec/RelActCalcAuto_ProfileScan.h"

namespace RelActCalcAutoImp
{
namespace ProfileLikelihood
{

/** Point budget.  These are scheduling numbers to be measured and adjusted, not derived constants;
 they live together so a measurement can move them without touching any logic.

 Only `sm_profile_max_points_per_quantity` is a hard cap.  The others are per-side allowances, and
 the arithmetic below deliberately closes with room to spare: with both sides' ladders and re-probes
 reserved up front, and the worst-case extension draw fitting inside the surplus, neither side can
 starve the other IN ANY ORDER.  That is why no priority rule is stated - no contention can arise.
 A profile whose answer depended on which direction happened to be scanned first would be exactly
 the defect the sibling machinery is built to avoid.
 */
constexpr std::size_t sm_profile_points_per_side         = 5;  //!< placement ladder (probe + 4)
constexpr std::size_t sm_profile_reprobe_per_side        = 1;  //!< rescale-failure re-probe
constexpr std::size_t sm_profile_max_extensions_per_side = 3;  //!< outward extensions
constexpr std::size_t sm_profile_max_points_per_quantity = 20; //!< the ONLY hard cap

/** Model-order ceiling for the REJECTION diagnostic, deliberately below the ceiling used for the
 crossing model.

 These are two different questions.  The crossing model wants enough freedom to follow the genuine
 shape of a clean profile; the rejection diagnostic asks the strictly local question "is this point
 isolated from its neighbours", for which a low-order smoother is the right tool precisely because
 it cannot chase the point under test.

 Measured, not assumed.  At a ceiling of 3, leaving one point out of five leaves four points for
 three parameters, which very nearly interpolates and then swings wildly on extrapolation: on a
 vector carrying one real outlier the diagnostic flagged a perfectly good interior point at +4.5
 alongside the true outlier at +35.  At a ceiling of 2 the same vector flags only the outlier, the
 recorded Pu-240 vector flags only its 23.129, and a clean quadratic flags nothing.
 */
constexpr std::size_t sm_rejection_model_order_ceiling = 2;

static_assert( 2*(sm_profile_points_per_side + sm_profile_reprobe_per_side
                  + sm_profile_max_extensions_per_side) <= sm_profile_max_points_per_quantity,
               "The point pool must be able to honour both sides' reserved allocations plus the"
               " worst-case extension draw simultaneously, or the budget becomes order-dependent." );


/** The vertex-anchored model.  `order` selects how many terms are active: 1 -> c2 only,
 2 -> c2 and c4 (biquadratic, still closed form), 3 -> all three.
 */
struct VertexModel
{
  double c2 = 0.0;
  double c3 = 0.0;
  double c4 = 0.0;
  std::size_t order = 0;

  double value( const double d ) const
  {
    return ((c4*d + c3)*d + c2)*d*d;
  }

  double slope( const double d ) const
  {
    return ((4.0*c4*d + 3.0*c3)*d + 2.0*c2)*d;
  }

  double curvature( const double d ) const
  {
    return (12.0*c4*d + 6.0*c3)*d + 2.0*c2;
  }

  /** `c2` is the curvature at the optimum; a non-positive one means the "optimum" is not a minimum
   in this coordinate, which is a failure and a strong hint the baseline is wrong.
   */
  bool usable() const
  {
    return (order >= 1) && (order <= 3) && std::isfinite(c2) && (c2 > 0.0)
           && std::isfinite(c3) && std::isfinite(c4);
  }
};


/** One sample reduced to what the model sees: a distance from the vertex, and a delta chi2. */
struct FitPoint
{
  double d = 0.0;
  double delta_chi2 = 0.0;
  /** Index into the caller's sample vector, so a rejected point can be reported by name. */
  std::size_t sample_index = 0;
};


/** Rejection tolerance of section 2.3.  It must sit above the genuine point-to-point scatter and
 below the outlier gap; on the recorded Pu-240 vector those are 1.08 and 12.26 respectively, and
 this returns 3.39.
 */
inline double rejection_tolerance( const double threshold_95 )
{
  return (std::max)( 1.0, 0.2*threshold_95 );
}


namespace detail
{

/** Solve a small symmetric normal-equation system by Gaussian elimination with partial pivoting,
 rejecting an ill-conditioned system rather than returning noise.

 Returns false when the pivots span more than `1e8`, which is what happens when the sampled points
 are clustered too tightly in `d` to identify the requested number of parameters.  Stepping the
 model order down on that signal is what makes clustering self-correcting, and it is why no ad-hoc
 "spread" heuristic is needed: the recorded Pu-240 vector (four points spanning 0.4% in `d`) is
 rejected at orders 3 and 2 and accepted at order 1, which is exactly the right answer.
 */
inline bool solve_small_system( std::vector<std::vector<double>> matrix,
                                std::vector<double> rhs,
                                std::vector<double> &solution )
{
  const std::size_t n = rhs.size();
  if( (n == 0) || (matrix.size() != n) )
    return false;

  double largest_pivot = 0.0;
  double smallest_pivot = std::numeric_limits<double>::infinity();

  for( std::size_t col = 0; col < n; ++col )
  {
    std::size_t pivot_row = col;
    for( std::size_t row = col + 1; row < n; ++row )
    {
      if( std::fabs(matrix[row][col]) > std::fabs(matrix[pivot_row][col]) )
        pivot_row = row;
    }

    const double pivot = matrix[pivot_row][col];
    if( !std::isfinite(pivot) || (pivot == 0.0) )
      return false;

    std::swap( matrix[col], matrix[pivot_row] );
    std::swap( rhs[col], rhs[pivot_row] );

    largest_pivot = (std::max)( largest_pivot, std::fabs(pivot) );
    smallest_pivot = (std::min)( smallest_pivot, std::fabs(pivot) );

    for( std::size_t row = col + 1; row < n; ++row )
    {
      const double factor = matrix[row][col]/pivot;
      if( !std::isfinite(factor) )
        return false;
      for( std::size_t k = col; k < n; ++k )
        matrix[row][k] -= factor*matrix[col][k];
      rhs[row] -= factor*rhs[col];
    }
  }

  if( !(smallest_pivot > 0.0) || !(largest_pivot > 0.0)
      || (largest_pivot/smallest_pivot > 1.0e8) )
    return false;

  solution.assign( n, 0.0 );
  for( std::size_t back = n; back > 0; --back )
  {
    const std::size_t row = back - 1;
    double accumulated = rhs[row];
    for( std::size_t k = row + 1; k < n; ++k )
      accumulated -= matrix[row][k]*solution[k];
    solution[row] = accumulated/matrix[row][row];
    if( !std::isfinite(solution[row]) )
      return false;
  }

  return true;
}

}//namespace detail


/** Least-squares fit of the vertex-anchored model at a requested order.

 The fit is done in a scaled coordinate `u = d/d_max` so the monomial powers stay near unity; on raw
 `d` of order 1e-4 the order-3 normal matrix entries would span `d^4` to `d^8`, i.e. 1e-16 to 1e-32,
 and would be numerically meaningless long before they were statistically meaningless.

 Returns a model with `order == 0` when the requested order cannot be identified from these points.
 */
inline VertexModel fit_vertex_model( const std::vector<FitPoint> &points, const std::size_t order )
{
  VertexModel model;
  if( (order < 1) || (order > 3) || (points.size() < order) )
    return model;

  double d_max = 0.0;
  for( const FitPoint &point : points )
  {
    if( !std::isfinite(point.d) || !std::isfinite(point.delta_chi2) || (point.d <= 0.0) )
      return model;
    d_max = (std::max)( d_max, point.d );
  }

  if( !(d_max > 0.0) )
    return model;

  // Active monomial exponents in the scaled coordinate, by order.
  std::vector<int> exponents;
  switch( order )
  {
    case 1:  exponents = {2};       break;
    case 2:  exponents = {2,4};     break;
    default: exponents = {2,3,4};   break;
  }

  const std::size_t n = exponents.size();
  std::vector<std::vector<double>> normal( n, std::vector<double>(n,0.0) );
  std::vector<double> rhs( n, 0.0 );

  for( const FitPoint &point : points )
  {
    const double u = point.d/d_max;
    std::vector<double> basis( n, 0.0 );
    for( std::size_t j = 0; j < n; ++j )
      basis[j] = std::pow( u, exponents[j] );

    for( std::size_t j = 0; j < n; ++j )
    {
      for( std::size_t k = 0; k < n; ++k )
        normal[j][k] += basis[j]*basis[k];
      rhs[j] += basis[j]*point.delta_chi2;
    }
  }

  std::vector<double> scaled;
  if( !detail::solve_small_system(normal,rhs,scaled) )
    return model;

  // Undo the scaling: a coefficient multiplying u^p multiplies d^p once divided by d_max^p.
  for( std::size_t j = 0; j < n; ++j )
  {
    const double coefficient = scaled[j]/std::pow( d_max, exponents[j] );
    if( !std::isfinite(coefficient) )
      return VertexModel();
    switch( exponents[j] )
    {
      case 2: model.c2 = coefficient; break;
      case 3: model.c3 = coefficient; break;
      default: model.c4 = coefficient; break;
    }
  }

  model.order = order;
  if( !model.usable() )
    return VertexModel();

  return model;
}


/** Solve `model(d) == threshold` for the smallest positive `d`, in closed form where possible.

 Returns nothing when the model never reaches the threshold, which for a concave quartic is a real
 possibility and must not be answered with the square root of a negative number.
 */
inline std::optional<double> solve_crossing( const VertexModel &model, const double threshold )
{
  if( !model.usable() || !std::isfinite(threshold) || !(threshold > 0.0) )
    return std::nullopt;

  const double quadratic_root = std::sqrt( threshold/model.c2 );
  if( !std::isfinite(quadratic_root) || !(quadratic_root > 0.0) )
    return std::nullopt;

  if( model.order == 1 )
    return quadratic_root;

  // Biquadratic: c4*d^4 + c2*d^2 - T == 0, solved for d^2.  The `+` branch is the smaller positive
  // root for BOTH signs of c4.
  const auto biquadratic_root = [&]() -> std::optional<double> {
    if( model.c4 == 0.0 )
      return quadratic_root;
    const double discriminant = model.c2*model.c2 + 4.0*model.c4*threshold;
    if( !std::isfinite(discriminant) || (discriminant < 0.0) )
      return std::nullopt;   // The model turns over below T and never reaches it.
    const double d_squared = (-model.c2 + std::sqrt(discriminant))/(2.0*model.c4);
    if( !std::isfinite(d_squared) || !(d_squared > 0.0) )
      return std::nullopt;
    return std::sqrt( d_squared );
  };

  if( model.order == 2 )
    return biquadratic_root();

  // Cubic term present: Newton from the quadratic root, which is arithmetic on an analytic model
  // rather than extra conditional fits, so it costs nanoseconds.
  double d = quadratic_root;
  for( std::size_t iteration = 0; iteration < 4; ++iteration )
  {
    const double slope = model.slope(d);
    if( !std::isfinite(slope) || !(slope > 0.0) )
    {
      // Stepping through a non-increasing region would leave the branch we want; bisect back
      // toward the quadratic root instead.
      d = 0.5*(d + quadratic_root);
      continue;
    }

    const double step = (model.value(d) - threshold)/slope;
    if( !std::isfinite(step) )
      break;

    const double next = d - step;
    if( !std::isfinite(next) || !(next > 0.0) )
    {
      d = 0.5*(d + quadratic_root);
      continue;
    }
    d = next;
  }

  if( std::isfinite(d) && (d > 0.0)
      && (std::fabs(model.value(d) - threshold) <= 1.0e-6*threshold) )
    return d;

  // Newton did not land on the threshold; the closed-form biquadratic is the honest fallback.
  return biquadratic_root();
}


/** Acceptance guard 1: the model must be CONVEX and increasing on `(0,d_root]`.

 An earlier formulation required only that the model be increasing, which is vacuous: for `c4 < 0`
 the selected branch always lies before the turn-over, so an increasing model is guaranteed and the
 test passes for EVERY concave fit - and a concave fit is precisely the failure this guards (a
 rejection that ate real curvature flattens the model and widens the interval).  A guard that cannot
 fail on the failure it targets is not a guard.

 A genuine profile likelihood is convex near its minimum for any well-posed problem, so local
 concavity in the FITTED model is a statement about the fit, not about the physics.

 `curvature()` is a quadratic in `d`, so checking its endpoints plus its interior stationary point
 (present only when `c4 < 0`) is the complete set.
 */
inline bool model_is_convex_and_increasing( const VertexModel &model, const double d_root )
{
  if( !model.usable() || !std::isfinite(d_root) || !(d_root > 0.0) )
    return false;

  if( !(model.slope(d_root) > 0.0) )
    return false;

  if( !(model.curvature(0.0) >= 0.0) || !(model.curvature(d_root) >= 0.0) )
    return false;

  if( model.c4 < 0.0 )
  {
    const double stationary = -model.c3/(4.0*model.c4);
    if( std::isfinite(stationary) && (stationary > 0.0) && (stationary < d_root)
        && !(model.curvature(stationary) >= 0.0) )
      return false;
  }

  return true;
}


/** Result of the one-sided rejection pass of section 2.3. */
struct RejectionResult
{
  std::vector<FitPoint> kept;
  /** Set when a point was dropped as an isolated outlier. */
  std::optional<FitPoint> dropped;
  /** The dropped point was the only evaluated sample at or above the 95% threshold, so the drop
   must not become final until an extension has tested the region.  Section 2.3 step 6. */
  bool drop_removed_sole_anchor = false;
};


/** One-sided leave-one-out outlier rejection.

 At most ONE point is dropped per side, in total - not per pass.  With five points and three
 parameters a single drop already leaves one degree of freedom; a second would leave the model
 interpolating, after which every residual is identically zero and the rejection loop would be
 certifying its own output.

 A point is never dropped for lying BELOW the model.  A spuriously low delta chi2 is not something a
 minimization can produce; a genuinely low one means the baseline was not the optimum, which the
 caller routes to the deferred-rebaseline flow instead.
 */
inline RejectionResult reject_one_outlier( const std::vector<FitPoint> &points,
                                           const double tau,
                                           const double threshold_95,
                                           const std::size_t order
                                                   = sm_rejection_model_order_ceiling )
{
  RejectionResult result;
  result.kept = points;

  if( (points.size() < 2) || !std::isfinite(tau) || !(tau > 0.0) )
    return result;

  // Leave-one-out, not in-sample: an outlier drags the in-sample fit toward itself and so hides its
  // own residual.
  std::size_t worst_index = points.size();
  double worst_excess = 0.0;
  for( std::size_t candidate = 0; candidate < points.size(); ++candidate )
  {
    std::vector<FitPoint> without;
    without.reserve( points.size() - 1 );
    for( std::size_t i = 0; i < points.size(); ++i )
    {
      if( i != candidate )
        without.push_back( points[i] );
    }

    std::size_t trial_order = (std::min)( order, without.size() );
    VertexModel model;
    while( (trial_order >= 1) && !model.usable() )
    {
      model = fit_vertex_model( without, trial_order );
      if( !model.usable() )
        --trial_order;
    }

    if( !model.usable() )
      continue;

    const double excess = points[candidate].delta_chi2 - model.value( points[candidate].d );
    if( (excess > tau) && (excess > worst_excess) )
    {
      worst_excess = excess;
      worst_index = candidate;
    }
  }

  if( worst_index >= points.size() )
    return result;

  std::vector<FitPoint> kept;
  kept.reserve( points.size() - 1 );
  for( std::size_t i = 0; i < points.size(); ++i )
  {
    if( i != worst_index )
      kept.push_back( points[i] );
  }

  // Survivors-consistency test.  If two or more points are still high after the drop, that is
  // curvature rather than a failed solve - real curvature lifts a point AND its neighbours - so
  // restore the candidate and fit all of them.  This is the guard that protects against rejection
  // eating genuine curvature, and it does so without any reference to ordering.
  std::size_t consistency_order = (std::min)( order, kept.size() );
  VertexModel kept_model;
  while( (consistency_order >= 1) && !kept_model.usable() )
  {
    kept_model = fit_vertex_model( kept, consistency_order );
    if( !kept_model.usable() )
      --consistency_order;
  }

  if( !kept_model.usable() )
    return result;

  double largest_residual = 0.0;
  for( const FitPoint &point : kept )
    largest_residual = (std::max)( largest_residual,
                                   std::fabs(point.delta_chi2 - kept_model.value(point.d)) );

  if( largest_residual > 0.5*tau )
    return result;   // Two or more points still high: this is curvature, so keep everything.

  // The drop stands.  Flag it when it removes the only evaluated anchor for the 95% crossing: the
  // rejection rule drops a point precisely when it is isolatedly high, which in a one-point outboard
  // region is the same sample the bracketing rule depends on.
  if( std::isfinite(threshold_95) && (threshold_95 > 0.0)
      && (points[worst_index].delta_chi2 >= threshold_95) )
  {
    const bool another_anchor = std::any_of( begin(kept), end(kept),
        [threshold_95]( const FitPoint &point ){ return point.delta_chi2 >= threshold_95; } );
    result.drop_removed_sole_anchor = !another_anchor;
  }

  result.dropped = points[worst_index];
  result.kept = std::move( kept );
  return result;
}


/** Fit the best model these points support, stepping the order down until one is identifiable.

 The order ceiling comes from the number of DISTINCT points, but the binding constraint in practice
 is conditioning, which `fit_vertex_model` reports by returning an unusable model.
 */
inline VertexModel fit_best_model( const std::vector<FitPoint> &points )
{
  const std::size_t ceiling = (std::min)( std::size_t(3), points.size() );
  for( std::size_t order = ceiling; order >= 1; --order )
  {
    const VertexModel model = fit_vertex_model( points, order );
    if( model.usable() )
      return model;
  }
  return VertexModel();
}


/** Solve for a crossing that satisfies the convexity guard, stepping the model order down on
 violation - which is what an earlier draft's "report a bound" would have silently converted into a
 `BoundaryLimited` answer for a side that was simply under-sampled.
 */
inline std::optional<double> guarded_crossing( const VertexModel &model, const double threshold,
                                               VertexModel &accepted )
{
  VertexModel trial = model;
  while( trial.usable() )
  {
    const std::optional<double> root = solve_crossing( trial, threshold );
    if( root && model_is_convex_and_increasing(trial,*root) )
    {
      accepted = trial;
      return root;
    }

    if( trial.order <= 1 )
      break;

    if( trial.order == 3 )
    {
      trial.c3 = 0.0;
      trial.order = 2;
    }else
    {
      trial.c4 = 0.0;
      trial.order = 1;
    }
  }

  return std::nullopt;
}

/** Per-side hygiene of section 2.5: order the samples by achieved reported distance, collapse
 near-duplicates, and detect a folded mapping.

 Samples are requested outward in the CONTROL coordinate, but the model is fitted in the REPORTED
 one, and the map between them is where all the nonlinearity lives.  Two things can therefore go
 wrong and neither is a modelling problem - both are data hygiene:

  - Two requests can land on the same reported value (the pin was silently clamped by the bounds
    projection, or the map is locally flat).  Keep the LOWER delta chi2, by the one-sided rule.
  - The reported ordering can invert relative to the control ordering, i.e. the map has folded.  Do
    not fit through a fold, and do not extend past one either: a sample beyond a fold is not
    outboard in the reported coordinate and cannot anchor a bracket.
 */
struct SideHygiene
{
  std::vector<FitPoint> points;      //!< sorted by `d` ascending, near-duplicates collapsed
  bool folded = false;
  double fold_first = 0.0;           //!< the two colliding reported values, for the diagnostic
  double fold_second = 0.0;
};

inline SideHygiene side_hygiene( const std::vector<Sample> &samples,
                                 const double baseline_reported,
                                 const double baseline_control )
{
  SideHygiene answer;
  if( samples.empty() )
    return answer;

  // Two samples are "the same measurement" only if they are close RELATIVE TO THE RANGE BEING
  // EXPLORED.  A fixed `1e-4*|q0|` is right when the scan spans much more than that - it is what
  // collapses the recorded Pu-240 near-duplicate pair - but it is catastrophic when the scan spans
  // LESS: where the reported quantity responds only weakly to the pinned coordinate, the whole
  // profile can live inside 3e-5 of fraction while the fixed tolerance is 9e-5, and then every
  // sample collapses into one and a side with a perfectly clear crossing reports that it never
  // bracketed anything.  Scaling by the observed span keeps both cases right.
  double reported_span = 0.0;
  for( const Sample &sample : samples )
  {
    const double distance = std::fabs( sample.reported_fraction - baseline_reported );
    if( std::isfinite(distance) )
      reported_span = (std::max)( reported_span, distance );
  }

  double duplicate_tolerance = (std::max)( 1.0e-12, 1.0e-4*std::fabs(baseline_reported) );
  if( reported_span > 0.0 )
    duplicate_tolerance = (std::max)( 1.0e-12,
                                      (std::min)(duplicate_tolerance,0.05*reported_span) );

  // A fold is defined against the CONTROL ordering, not the order the points happened to be
  // requested in: a probe that overshoots is rescaled correctly and the ladder then lands INSIDE
  // it, so request order is not outward order, and inferring a fold from it would misfire on the
  // cheapest and most desirable outcome there is.
  std::vector<std::size_t> by_control( samples.size() );
  for( std::size_t i = 0; i < by_control.size(); ++i )
    by_control[i] = i;
  std::stable_sort( begin(by_control), end(by_control),
      [&]( const std::size_t lhs, const std::size_t rhs ){
        return std::fabs(samples[lhs].control_fraction - baseline_control)
             < std::fabs(samples[rhs].control_fraction - baseline_control);
      } );

  // Walking outward in the control coordinate, the reported distance must not decrease.  Where it
  // does, the map has folded; truncate the side there rather than fitting through it.
  std::vector<FitPoint> ordered;
  double previous_distance = 0.0;
  for( std::size_t position = 0; position < by_control.size(); ++position )
  {
    const std::size_t index = by_control[position];
    const double distance = std::fabs( samples[index].reported_fraction - baseline_reported );
    if( !std::isfinite(distance) )
      break;

    if( (position > 0) && (distance < previous_distance - duplicate_tolerance) )
    {
      answer.folded = true;
      answer.fold_first = samples[by_control[position-1]].reported_fraction;
      answer.fold_second = samples[index].reported_fraction;
      break;
    }

    FitPoint point;
    point.d = distance;
    point.delta_chi2 = samples[index].delta_chi2;
    point.sample_index = index;
    ordered.push_back( point );
    previous_distance = distance;
  }

  std::stable_sort( begin(ordered), end(ordered),
                    []( const FitPoint &lhs, const FitPoint &rhs ){ return lhs.d < rhs.d; } );

  for( const FitPoint &point : ordered )
  {
    if( !answer.points.empty()
        && ((point.d - answer.points.back().d) <= duplicate_tolerance) )
    {
      // Same reported coordinate reached twice: keep the lower objective.
      if( point.delta_chi2 < answer.points.back().delta_chi2 )
        answer.points.back() = point;
      continue;
    }
    answer.points.push_back( point );
  }

  // A point at the vertex itself carries no information and would divide by zero in the scaled fit.
  answer.points.erase( std::remove_if( begin(answer.points), end(answer.points),
                                       []( const FitPoint &point ){ return !(point.d > 0.0); } ),
                       end(answer.points) );

  return answer;
}


/** Fit a bounded scalar profile at two nested likelihood thresholds by pinning and modelling.

 Drop-in for `scan()`, with one added argument: the local one-sigma in the CONTROL coordinate, used
 only to place the first probe.  Zero or non-finite means unavailable - the common case for an
 automatically-profiled weak quantity, whose trigger is literally "the 95% band already covers the
 whole feasible interval" - and the probe then falls back to a fraction of the directional span.
 The placement is self-calibrating either way: it rescales off the probe's measured delta chi2 and
 so never assumes how sigma relates to the threshold.

 `evaluate` must return the objective difference from the same independently evaluated baseline, and
 may map the control coordinate to a different reported coordinate.
 */
inline ScanResult fit_profile( const double baseline_control,
                               const double baseline_reported,
                               const double lower_control,
                               const double upper_control,
                               const std::array<double,2> thresholds,
                               const std::size_t max_points,
                               const Evaluator &evaluate,
                               const double control_one_sigma = 0.0 )
{
  ScanResult answer;
  answer.intervals[0].delta_chi2 = thresholds[0];
  answer.intervals[1].delta_chi2 = thresholds[1];

  if( !evaluate || !std::isfinite(baseline_control) || !std::isfinite(baseline_reported)
      || !std::isfinite(lower_control) || !std::isfinite(upper_control)
      || (lower_control > baseline_control) || (baseline_control > upper_control)
      || !std::isfinite(thresholds[0]) || !std::isfinite(thresholds[1])
      || (max_points == 0)
      || !(thresholds[0] > 0.0) || !(thresholds[1] > thresholds[0]) )
  {
    answer.status = ScanStatus::Failed;
    answer.diagnostic = "Invalid profile bounds, baseline, thresholds, or evaluator.";
    return answer;
  }

  const double full_span = upper_control - lower_control;
  const double threshold_95 = thresholds[1];
  const double tau = rejection_tolerance( threshold_95 );
  const std::size_t point_pool = (std::min)( max_points, sm_profile_max_points_per_quantity );

  std::size_t points_used = 0;
  std::size_t discarded_points = 0;
  std::string last_discard_diagnostic;
  bool hard_failure = false;

  /** Evaluate one control value, or report why the whole scan must stop. */
  const auto request = [&]( const double control, const Direction direction,
                            Sample &out ) -> bool {
    if( points_used >= point_pool )
      return false;
    ++points_used;

    const Evaluation evaluation = evaluate( control, direction, points_used );
    switch( evaluation.status )
    {
      case EvaluationStatus::Canceled:
        answer.status = ScanStatus::Canceled;
        answer.diagnostic = evaluation.diagnostic.empty()
                          ? "Mass-fraction profiling was canceled." : evaluation.diagnostic;
        hard_failure = true;
        return false;

      case EvaluationStatus::BetterBaseline:
        answer.status = ScanStatus::BetterBaseline;
        answer.diagnostic = evaluation.diagnostic.empty()
            ? "A conditional profile fit found a better unconstrained point."
            : evaluation.diagnostic;
        hard_failure = true;
        return false;

      case EvaluationStatus::FitCapReached:
        // The pool is this function's own accounting; an evaluator reporting its own cap is a
        // budget the caller imposed, and is not recoverable here.
        answer.status = ScanStatus::FitCapReached;
        answer.diagnostic = evaluation.diagnostic.empty()
                          ? "The conditional-fit cap was reached." : evaluation.diagnostic;
        hard_failure = true;
        return false;

      case EvaluationStatus::Failed:
        // DISCARD the point, do not fail the scan.  This design assumes some conditionals will come
        // back unusable - that assumption is the entire basis of the one-sided rejection rule - so
        // letting one of them terminate the whole profile would throw away every good sample
        // alongside it.  The point is charged to the pool (it cost a solve) and the placement moves
        // on; a side that ends with too few usable samples fails on its own, with a diagnostic that
        // says so.
        ++discarded_points;
        if( last_discard_diagnostic.empty() )
          last_discard_diagnostic = evaluation.diagnostic;
        return false;

      case EvaluationStatus::Success:
        break;
    }//switch( evaluation.status )

    if( !detail::finite_sample(evaluation) )
    {
      answer.status = ScanStatus::Failed;
      answer.diagnostic = "A conditional profile fit returned a non-finite result.";
      hard_failure = true;
      return false;
    }

    // POINTS, not `ceres::Solve` calls.  The pool this function enforces is a point budget, and
    // every consumer downstream reports this field as "points"; the solve count is reported
    // separately by the host, which is what makes a slow scan diagnosable as "more points" versus
    // "more solves per point" rather than one conflated number.
    ++answer.num_fits;
    // Always the REQUESTED control coordinate.  The scanner this replaces overwrote it with the
    // reported value at a structural bound, which is harmless when both are mass fractions but
    // actively wrong here: the control coordinate is a fit parameter on a completely different
    // scale, and mixing the two would corrupt placement, fold detection and duplicate rejection.
    // `reached_feasible_bound` carries that information instead.
    out.control_fraction = control;
    out.reported_fraction = evaluation.reported_fraction;
    // A conditional cannot fall below its own unconstrained minimum; a negative value is solver
    // noise about zero, not information.
    out.delta_chi2 = (std::max)( 0.0, evaluation.delta_chi2 );
    out.reached_feasible_bound = evaluation.reached_feasible_bound;
    answer.samples.push_back( out );
    return true;
  };

  struct SideResult
  {
    std::array<Endpoint,2> endpoints;
    bool failed = false;
    bool entirely_inside = false;    //!< reached the bound with every sample below the 95% threshold
    std::string diagnostic;
  };

  const auto scan_side = [&]( const Direction direction ) -> SideResult {
    SideResult side;
    const double sign = (direction == Direction::Lower) ? -1.0 : 1.0;
    const double bound = (direction == Direction::Lower) ? lower_control : upper_control;
    const double span = (direction == Direction::Lower)
                      ? (baseline_control - lower_control)
                      : (upper_control - baseline_control);

    for( Endpoint &endpoint : side.endpoints )
    {
      endpoint.reported_fraction = baseline_reported;
      endpoint.likelihood_crossing = false;
    }

    if( !(span > 0.0) )
    {
      // No room on this side at all: the baseline IS the bound.
      side.entirely_inside = true;
      return side;
    }

    std::vector<Sample> samples;
    std::size_t ladder_used = 0;
    std::size_t reprobe_used = 0;
    std::size_t extensions_used = 0;

    const auto control_at = [&]( const double distance ){
      return baseline_control + sign*(std::min)( distance, span );
    };

    const auto already_evaluated = [&]( const double distance ){
      const double tolerance = (std::max)( 1.0e-12, 1.0e-9*span );
      for( const Sample &sample : samples )
      {
        if( std::fabs(std::fabs(sample.control_fraction - baseline_control) - distance)
              <= tolerance )
          return true;
      }
      return false;
    };

    const auto take = [&]( const double distance ) -> bool {
      if( already_evaluated(distance) )
        return false;
      Sample sample;
      if( !request(control_at(distance),direction,sample) )
        return false;
      samples.push_back( sample );
      return true;
    };

    // --- Probe -------------------------------------------------------------------------------
    //
    // The probe exists to measure delta chi2 at a KNOWN distance so the rescale can place the
    // ladder; any interior point does that job.  It is therefore held to a quarter of the span,
    // which matters more than it looks: where the reported quantity responds only weakly to the
    // pinned coordinate, the one-sigma distance carried into that coordinate can exceed the whole
    // span, and an unclamped probe then lands exactly ON the bound - an activity of zero, say -
    // where the model cannot be evaluated at all.  Every ladder point would clip to the same
    // unusable extreme and the side would produce no samples whatsoever.
    double probe_distance = (std::isfinite(control_one_sigma) && (control_one_sigma > 0.0))
                          ? (std::max)( 0.98*control_one_sigma, 0.02*span )
                          : 0.02*span;
    probe_distance = (std::min)( probe_distance, 0.25*span );
    const bool probe_usable = take( probe_distance );
    ++ladder_used;
    if( hard_failure )
    {
      side.failed = true;
      return side;
    }

    // A discarded probe costs the rescale, not the side: fall back to spreading over what is
    // reachable, which is the same placement the clipping branch uses.
    double best_probe_distance = probe_distance;
    double best_probe_chi2 = probe_usable ? samples.back().delta_chi2 : 0.0;

    // A probe that barely moved the objective gives a rescale factor that is all noise; spend the
    // reserved re-probe far enough out to measure something.  The opposite case needs nothing: a
    // probe already past the 95% threshold has bracketed the crossing by itself, which is the
    // cheapest possible outcome.
    if( probe_usable && (best_probe_chi2 < 0.05*thresholds[0])
        && (reprobe_used < sm_profile_reprobe_per_side) && (probe_distance < span) )
    {
      const double retry = (std::min)( 10.0*probe_distance, span );
      if( take(retry) )
      {
        ++reprobe_used;
        if( samples.back().delta_chi2 > best_probe_chi2 )
        {
          best_probe_chi2 = samples.back().delta_chi2;
          best_probe_distance = retry;
        }
      }
    }

    // --- Ladder ------------------------------------------------------------------------------
    // Under the quadratic model the 95% crossing sits at d1*sqrt(T95/dchi2_1).
    const double predicted_95 = (probe_usable && (best_probe_chi2 > 0.0))
                              ? best_probe_distance*std::sqrt(threshold_95/best_probe_chi2)
                              : std::numeric_limits<double>::infinity();

    // Two of the four ladder points are deliberately OUTBOARD of the predicted crossing.  That is
    // required, not incidental: the bracketing rule needs an evaluated sample at or above the
    // threshold, and the rejection rule drops a point precisely when it is isolatedly high - so with
    // only one point out there, the anchor would be the sample most likely to be discarded, and the
    // restore-on-consistency test could never fire in a one-point region.
    const std::array<double,4> crossing_ladder{{0.55, 0.85, 1.15, 1.45}};
    // When the predicted crossing is past the feasible bound - the COMMON case, since the automatic
    // trigger is literally "the 95% band already covers the whole feasible interval" - abandon
    // crossing placement and spread over what is actually reachable.
    const std::array<double,4> reachable_ladder{{0.25, 0.50, 0.75, 1.00}};
    const bool crossing_is_reachable = std::isfinite(predicted_95) && (predicted_95 <= span);
    const std::array<double,4> &ladder = crossing_is_reachable ? crossing_ladder : reachable_ladder;
    const double ladder_scale = crossing_is_reachable ? predicted_95 : span;

    for( const double multiplier : ladder )
    {
      if( ladder_used >= sm_profile_points_per_side )
        break;
      if( take((std::min)(multiplier*ladder_scale,span)) )
        ++ladder_used;
      if( hard_failure )
      {
        side.failed = true;
        return side;
      }
    }

    if( samples.empty() )
    {
      side.failed = true;
      side.diagnostic = "No usable profile sample could be evaluated on this side"
          + (last_discard_diagnostic.empty() ? std::string(".")
                                             : (": " + last_discard_diagnostic));
      return side;
    }

    // --- Fit, extend, classify ---------------------------------------------------------------
    bool rejection_done = false;
    std::optional<FitPoint> dropped_point;
    SideHygiene hygiene;
    VertexModel model;
    std::optional<double> root_95;
    VertexModel accepted;

    const auto reached_bound = [&](){
      for( const Sample &sample : samples )
      {
        if( sample.reached_feasible_bound
            || detail::at_control_bound(sample,bound,full_span) )
          return true;
      }
      return false;
    };

    // The innermost EVALUATED sample at or above a threshold - the outboard bracket anchor.  Stated
    // on evaluated values, not on the model's values at those points: a model sitting above the data
    // could otherwise satisfy the test while no measurement ever reached the threshold, which is
    // extrapolation wearing a bracket's clothing.
    const auto anchor_distance = [&]( const double threshold ) -> std::optional<double> {
      std::optional<double> innermost;
      for( const FitPoint &point : hygiene.points )
      {
        if( (point.delta_chi2 >= threshold) && (!innermost || (point.d < *innermost)) )
          innermost = point.d;
      }
      return innermost;
    };

    /** Is the crossing bracketed by an evaluated sample?

     The comparison carries a relative hair of 1e-9 because a sample can land exactly ON the
     threshold - a probe placed at one sigma does so almost by construction at the 68% level - and
     then the root and the anchor are the same number computed two different ways.  Without it a
     perfect bracket is converted into a failure by round-off.  1e-9 relative is many orders below
     any resolution these samples carry, so it permits no extrapolation in any meaningful sense.
     */
    const auto is_bracketed = []( const double root, const double anchor ){
      return root <= anchor*(1.0 + 1.0e-9) + 1.0e-15;
    };

    // Outermost EVALUATED control distance.  Deliberately NOT the outermost survivor: a drop must
    // not move the next measurement back inside a point already discarded.  Rejection governs what
    // the fit sees; it must not govern where the next measurement goes.
    const auto outermost_evaluated_distance = [&](){
      double furthest = 0.0;
      for( const Sample &sample : samples )
        furthest = (std::max)( furthest,
                               std::fabs(sample.control_fraction - baseline_control) );
      return furthest;
    };

    // Measured, not assumed: the reported-per-control ratio comes from the outermost sample, so the
    // extension step needs no analytic derivative of a map we only ever sample.
    const auto reported_per_control = [&]() -> double {
      const double control_distance = outermost_evaluated_distance();
      if( !(control_distance > 0.0) )
        return 0.0;
      double reported_distance = 0.0;
      for( const Sample &sample : samples )
        reported_distance = (std::max)( reported_distance,
                                        std::fabs(sample.reported_fraction - baseline_reported) );
      return reported_distance/control_distance;
    };

    const auto extend_outward = [&]() -> bool {
      if( (extensions_used >= sm_profile_max_extensions_per_side)
          || (points_used >= point_pool) )
        return false;

      const double out = outermost_evaluated_distance();
      if( !(out < span) )
        return false;   // Already at the bound; there is nowhere further to go.

      // 1.3*out is a floor guaranteeing geometric progress, so a badly-conditioned prediction
      // cannot produce an infinitesimal step and the bound stays reachable.  1.15*predicted aims
      // deliberately PAST the crossing: undershooting costs a whole conditional solve, while
      // overshooting costs only a slightly wider bracket that the interpolating fit absorbs.
      double next = 1.3*out;
      const double ratio = reported_per_control();
      if( root_95 && (ratio > 0.0) )
        next = (std::max)( next, 1.15*(*root_95)/ratio );
      next = (std::min)( next, span );

      if( !take(next) )
        return false;
      ++extensions_used;
      return true;
    };

    for( std::size_t pass = 0; pass < (2*sm_profile_max_extensions_per_side + 4); ++pass )
    {
      hygiene = side_hygiene( samples, baseline_reported, baseline_control );

      if( hygiene.points.empty() )
      {
        side.failed = true;
        side.diagnostic = hygiene.folded
            ? "The control-to-reported mapping folded before any usable profile sample."
            : "No usable profile samples remain on this side.";
        return side;
      }

      std::vector<FitPoint> fitted_points = hygiene.points;
      if( !rejection_done )
      {
        const RejectionResult rejection
            = reject_one_outlier( hygiene.points, tau, threshold_95 );

        // A drop that removes the only bracketing anchor is not final: spend one point outward and
        // re-run the whole rejection over the enlarged set.  If the new sample confirms the region
        // is high the two are mutually consistent, the candidate is restored, and the "outlier" is
        // revealed as real curvature.  If it agrees with the survivors' trend the drop stands and
        // the new sample is the anchor.  Either way the bracket survives and the decision rests on
        // a measurement rather than on a rule.
        if( rejection.drop_removed_sole_anchor && !hygiene.folded )
        {
          if( extend_outward() )
            continue;
          if( hard_failure )
          {
            side.failed = true;
            return side;
          }
        }

        fitted_points = rejection.kept;
        dropped_point = rejection.dropped;
        rejection_done = true;
      }else
      {
        if( dropped_point )
        {
          fitted_points.erase(
              std::remove_if( begin(fitted_points), end(fitted_points),
                  [&]( const FitPoint &point ){
                    return point.sample_index == dropped_point->sample_index;
                  } ),
              end(fitted_points) );
        }
      }

      if( fitted_points.empty() )
      {
        side.failed = true;
        side.diagnostic = "Every profile sample on this side was rejected.";
        return side;
      }

      model = fit_best_model( fitted_points );
      if( !model.usable() )
      {
        // No usable curvature.  WHY matters, and conflating the two cases misreports a legitimate
        // answer as a failure:
        //
        //  - Nothing on this side ever approached the threshold.  The profile is FLAT along this
        //    coordinate - the objective genuinely does not respond to it - which is what
        //    non-identifiability IS, not a failure to measure it.  Fall through and let the bound
        //    and interior classification below say so.  Observed on Pu-241 in 610-775 keV, whose
        //    upper side returns delta chi2 of exactly zero at every sampled point.
        //  - Some sample DID cross the threshold, so a crossing exists and the model cannot
        //    describe it.  That is a real failure, and a strong hint the baseline is wrong.
        const bool any_sample_reached_threshold = std::any_of( begin(hygiene.points),
            end(hygiene.points), [threshold_95]( const FitPoint &point ){
              return point.delta_chi2 >= threshold_95;
            } );

        if( any_sample_reached_threshold )
        {
          side.failed = true;
          side.diagnostic = "The profile has no positive curvature at the optimum in this"
                            " direction even though a sample crossed the threshold, so the reported"
                            " solution is not a minimum along it.";
          return side;
        }
      }

      root_95 = model.usable() ? guarded_crossing( model, threshold_95, accepted )
                               : std::nullopt;
      const std::optional<double> anchor_95 = anchor_distance( threshold_95 );
      const bool bracketed = root_95 && anchor_95 && is_bracketed(*root_95,*anchor_95);

      // An anchor is enough to stop travelling outward, even when the model's own root sits beyond
      // it.  That happens when the true profile is steeper than the fitted polynomial, and no
      // amount of further extension fixes it - the measurements ALREADY bracket the threshold, so
      // spending points chasing the model would burn the budget and then report failure for a
      // quantity whose data answer the question.
      if( bracketed || anchor_95 || hygiene.folded )
        break;

      if( reached_bound() )
        break;   // The likelihood genuinely never reaches the threshold inside the feasible set.

      if( !extend_outward() )
        break;
      if( hard_failure )
      {
        side.failed = true;
        return side;
      }
    }//for( extension passes )

    const bool at_bound = reached_bound();

    // Every sample below the 95% threshold with the bound reached is a legitimate, correct answer -
    // the likelihood really does stay inside the threshold across the whole feasible set - and must
    // never be conflated with "we did not push far enough".
    side.entirely_inside = at_bound
        && std::all_of( begin(hygiene.points), end(hygiene.points),
                        [threshold_95]( const FitPoint &point ){
                          return point.delta_chi2 < threshold_95;
                        } );

    // The bound's reported value comes from the sample furthest out in CONTROL, which after a
    // rescaled ladder is not necessarily the last one requested.
    double bound_reported = baseline_reported;
    double furthest_control = -1.0;
    for( const Sample &sample : samples )
    {
      const double distance = std::fabs( sample.control_fraction - baseline_control );
      if( distance > furthest_control )
      {
        furthest_control = distance;
        bound_reported = sample.reported_fraction;
      }
    }

    for( std::size_t level = 0; level < 2; ++level )
    {
      const double threshold = thresholds[level];
      const std::optional<double> anchor = anchor_distance( threshold );

      VertexModel level_model;
      const std::optional<double> root = model.usable()
                                       ? guarded_crossing( model, threshold, level_model )
                                       : std::nullopt;

      if( root && anchor && is_bracketed(*root,*anchor) && !hygiene.folded )
      {
        side.endpoints[level].reported_fraction = baseline_reported + sign*(*root);
        side.endpoints[level].likelihood_crossing = true;
        continue;
      }

      // The model underpredicts at the innermost evaluated sample that reached the threshold.  The
      // crossing is bracketed all the same - the vertex is below it by construction and this
      // MEASUREMENT is at or above it - so report the measurement rather than a model root that
      // would extrapolate past it.  Being an evaluated point at or beyond the true crossing, it
      // errs WIDE, which is the safe direction and the same reasoning the terminal-bracket rule
      // has always used.  Note this necessarily precedes the bound test: a sample that reached the
      // threshold proves a crossing exists inside the feasible set.
      if( anchor && !hygiene.folded )
      {
        side.endpoints[level].reported_fraction = baseline_reported + sign*(*anchor);
        side.endpoints[level].likelihood_crossing = true;
        continue;
      }

      if( at_bound )
      {
        // The bound itself, exactly as before: a physical or input limit, not a crossing.
        side.endpoints[level].reported_fraction = bound_reported;
        side.endpoints[level].likelihood_crossing = false;
        continue;
      }

      // Nothing above reported an endpoint, so no evaluated sample ever reached the threshold and
      // the feasible bound was not arrived at either.  That is an honest "we never got near it
      // within budget", and it must never be dressed up as a crossing or as a bound.
      //
      // Note there is deliberately no separate straddling-pair fallback here.  The anchor clamp
      // above already returns what the deleted Brent scanner's `terminal_cap_crossing_sample`
      // did - the innermost EVALUATED sample at or above the threshold - and it does so without
      // that function's acceptance rule, which was written for a Brent bracket between two samples
      // and rejected precisely the case this engine must serve: a vertex-to-anchor bracket is wide
      // by construction, so the rule would refuse an endpoint whose data plainly bracket the
      // threshold.  Erring wide is the safe direction and is the same reasoning that function was
      // written on.
      side.failed = true;
      if( hygiene.folded )
      {
        side.diagnostic = "The control-to-reported mapping folded between "
            + std::to_string(hygiene.fold_first) + " and " + std::to_string(hygiene.fold_second)
            + ", so no crossing outboard of the fold can be interpolated.";
      }else
      {
        double best = 0.0;
        for( const FitPoint &point : hygiene.points )
          best = (std::max)( best, point.delta_chi2 );
        side.diagnostic = "The profile point budget was exhausted without bracketing the threshold:"
            " " + std::to_string(samples.size()) + " points and "
            + std::to_string(extensions_used) + " extensions reached delta chi2 "
            + std::to_string(best) + " against a threshold of " + std::to_string(threshold)
            + " (model order " + std::to_string(model.order)
            + ", c2=" + std::to_string(model.c2)
            + ", root " + (root ? std::to_string(*root) : std::string("none"))
            + ", anchor " + (anchor ? std::to_string(*anchor) : std::string("none"))
            + ", distinct samples " + std::to_string(hygiene.points.size()) + ").";
      }
      return side;
    }//for( each threshold level )

    return side;
  };

  const SideResult lower_side = scan_side( Direction::Lower );
  if( hard_failure )
    return answer;
  const SideResult upper_side = scan_side( Direction::Upper );
  if( hard_failure )
    return answer;

  for( std::size_t level = 0; level < 2; ++level )
  {
    answer.intervals[level].lower = lower_side.endpoints[level];
    answer.intervals[level].upper = upper_side.endpoints[level];
    if( answer.intervals[level].lower.reported_fraction
          > answer.intervals[level].upper.reported_fraction )
      std::swap( answer.intervals[level].lower, answer.intervals[level].upper );
  }

  if( lower_side.failed || upper_side.failed )
  {
    answer.status = ScanStatus::Failed;
    answer.diagnostic = lower_side.failed ? lower_side.diagnostic : upper_side.diagnostic;
    return answer;
  }

  if( lower_side.entirely_inside && upper_side.entirely_inside )
  {
    answer.status = ScanStatus::NonIdentifiable;
    answer.diagnostic = "The entire feasible mass-fraction range remains inside the 95% likelihood"
                        " threshold.";
    return answer;
  }

  bool has_boundary = false;
  for( const Interval &interval : answer.intervals )
    has_boundary = has_boundary || !interval.lower.likelihood_crossing
                                || !interval.upper.likelihood_crossing;

  answer.status = has_boundary ? ScanStatus::BoundaryLimited : ScanStatus::Complete;
  answer.diagnostic = has_boundary
      ? "At least one confidence endpoint is set by a physical or input bound."
      : "Both likelihood directions were bracketed by evaluated samples.";
  if( discarded_points > 0 )
    answer.diagnostic += "  " + std::to_string(discarded_points) + " conditional point(s) were"
        " discarded as unconverged and excluded from the fit.";
  return answer;
}


}//namespace ProfileLikelihood
}//namespace RelActCalcAutoImp

#endif //InterSpec_RelActCalcAuto_ProfileFit_h
