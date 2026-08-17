/* InterSpec: quantity-specific diagnostics for a truncated fit covariance. */
#ifndef InterSpec_RelActCalcAuto_CovarianceQuality_h
#define InterSpec_RelActCalcAuto_CovarianceQuality_h

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace RelActCalcAutoImp
{

/** Result of projecting one derived quantity onto the natural, column-equilibrated
    tangent coordinates used by the RelActAuto covariance SVD.

    Squared sensitivities are expressed in those dimensionless coordinates.  This
    makes the decision invariant to changing a fit parameter from, for example,
    years to seconds or activity to kilo-activity. */
struct QuantityCovarianceDiagnostic
{
  bool dimensions_valid = false;
  bool finite = false;
  bool depends_on_dropped_direction = false;
  bool depends_on_active_bound = false;
  double natural_sensitivity_squared = 0.0;
  double dropped_sensitivity_squared = 0.0;
  double active_bound_sensitivity_squared = 0.0;
  double dropped_fraction = 0.0;
  double active_bound_fraction = 0.0;
};


/** Assess whether a derived quantity's local Gaussian covariance omits a material
    part of that quantity's sensitivity.

    `ambient_from_equilibrated_tangent` is B = PlusJacobian * D, where D contains
    the inverse objective-Jacobian column norms (and a finite characteristic scale
    for a truly zero column).  Each row is an ambient parameter and each column is
    one tangent coordinate.  `dropped_ambient_directions` contains B*v for every
    right-singular direction discarded by the covariance pseudo-inverse.

    Two different pathologies are measured, and they deserve different thresholds:

    - A *dropped* direction was removed from the covariance pseudo-inverse, so its
      variance contribution is exactly zero.  Even a small projection onto it can
      make a reported sigma meaningless, hence the strict 1e-6 default.
    - An *active bound* merely truncates: Ceres reports ~0 variance for a pinned
      parameter, which understates - but does not invalidate - a quantity that only
      touches it.  A pinned parameter is very often perfectly well determined (an
      isotope that is genuinely absent, or a constrained element total sitting at
      its cap by construction), so requiring a *material* share here matches
      `sensitivity_rides_on_bound`, the sibling test used by the uncertainty-widening
      path, rather than flagging anything that touches a bound at all.

    Neither flag means "this uncertainty is wrong" on its own - they say the local
    covariance omits something.  Callers must combine them with a plausibility test on
    the reported uncertainty (see `reliability_floored_uncert`): a credible uncertainty
    is not made untrustworthy by an unrelated parameter sitting on a bound. */
inline QuantityCovarianceDiagnostic assess_quantity_covariance(
    const std::vector<double> &ambient_gradient,
    const std::vector<std::vector<double>> &ambient_from_equilibrated_tangent,
    const std::vector<std::vector<double>> &dropped_ambient_directions,
    const std::vector<char> &ambient_parameter_at_bound,
    const double material_squared_fraction = 1.0e-6,
    const double bound_squared_fraction = 0.25 )
{
  QuantityCovarianceDiagnostic answer;
  const size_t num_ambient = ambient_gradient.size();
  if( (num_ambient == 0)
      || (ambient_from_equilibrated_tangent.size() != num_ambient)
      || (!ambient_parameter_at_bound.empty()
          && (ambient_parameter_at_bound.size() != num_ambient))
      || !(material_squared_fraction >= 0.0)
      || !std::isfinite(material_squared_fraction)
      || !(bound_squared_fraction >= 0.0)
      || !std::isfinite(bound_squared_fraction) )
    return answer;

  const size_t num_tangent = ambient_from_equilibrated_tangent.front().size();
  if( num_tangent == 0 )
    return answer;
  for( const std::vector<double> &row : ambient_from_equilibrated_tangent )
    if( row.size() != num_tangent )
      return answer;
  for( const std::vector<double> &direction : dropped_ambient_directions )
    if( direction.size() != num_ambient )
      return answer;

  answer.dimensions_valid = true;

  std::vector<double> tangent_gradient( num_tangent, 0.0 );
  for( size_t ambient = 0; ambient < num_ambient; ++ambient )
  {
    const double gradient = ambient_gradient[ambient];
    if( !std::isfinite(gradient) )
      return answer;
    for( size_t tangent = 0; tangent < num_tangent; ++tangent )
    {
      const double mapping = ambient_from_equilibrated_tangent[ambient][tangent];
      if( !std::isfinite(mapping) )
        return answer;
      tangent_gradient[tangent] += gradient*mapping;
    }
  }

  for( const double derivative : tangent_gradient )
    answer.natural_sensitivity_squared += derivative*derivative;

  for( const std::vector<double> &direction : dropped_ambient_directions )
  {
    double projection = 0.0;
    for( size_t ambient = 0; ambient < num_ambient; ++ambient )
    {
      if( !std::isfinite(direction[ambient]) )
        return answer;
      projection += ambient_gradient[ambient]*direction[ambient];
    }
    answer.dropped_sensitivity_squared += projection*projection;
  }

  // The live Ceres manifold is a SubsetManifold today, for which a row of B has
  // at most one nonzero entry.  Using its row norm nevertheless keeps this
  // diagnostic correct and scale-aware if another orthonormal manifold is used.
  if( !ambient_parameter_at_bound.empty() )
  {
    for( size_t ambient = 0; ambient < num_ambient; ++ambient )
    {
      if( !ambient_parameter_at_bound[ambient] )
        continue;
      double row_scale_squared = 0.0;
      for( const double mapping : ambient_from_equilibrated_tangent[ambient] )
        row_scale_squared += mapping*mapping;
      answer.active_bound_sensitivity_squared
          += ambient_gradient[ambient]*ambient_gradient[ambient]*row_scale_squared;
    }
  }

  if( !std::isfinite(answer.natural_sensitivity_squared)
      || !std::isfinite(answer.dropped_sensitivity_squared)
      || !std::isfinite(answer.active_bound_sensitivity_squared) )
    return answer;

  answer.finite = true;
  if( !(answer.natural_sensitivity_squared
        > 64.0*std::numeric_limits<double>::min()) )
    return answer;

  answer.dropped_fraction = std::clamp(
      answer.dropped_sensitivity_squared/answer.natural_sensitivity_squared,0.0,1.0);
  answer.active_bound_fraction = std::clamp(
      answer.active_bound_sensitivity_squared/answer.natural_sensitivity_squared,0.0,1.0);
  answer.depends_on_dropped_direction
      = answer.dropped_sensitivity_squared
          > material_squared_fraction*answer.natural_sensitivity_squared;
  answer.depends_on_active_bound
      = answer.active_bound_sensitivity_squared
          > bound_squared_fraction*answer.natural_sensitivity_squared;
  return answer;
}

}//namespace RelActCalcAutoImp

#endif //InterSpec_RelActCalcAuto_CovarianceQuality_h
