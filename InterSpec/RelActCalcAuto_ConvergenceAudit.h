/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
*/

#ifndef InterSpec_RelActCalcAuto_ConvergenceAudit_h
#define InterSpec_RelActCalcAuto_ConvergenceAudit_h

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace RelActCalcAutoImp
{

/** Thresholds used to decide whether a Ceres NO_CONVERGENCE result is nevertheless a usable,
 * stationary budget-limited point.  This type is public only to the focused numerical oracle; it
 * is not part of the user-facing RelActAuto API. */
struct ConvergenceAuditThresholds
{
  std::size_t envelope_tail_length = 10;
  double stagnant_envelope_change = 1.0e-7;
  double maximum_tail_rise = 1.0e-7;
  double stationary_projected_gradient = 1.0e-4;
  double small_relative_step = 1.0e-5;
  double collapsed_trust_region_radius = 1.0e-12;
};

/** Independently measured inputs to the budget-exhaustion decision.
 *
 * `successful_full_objectives` must be in chronological order and in the same sum-of-squares
 * convention as `fresh_full_objective` (Ceres' one-half cost is multiplied by two by the caller).
 * The two fresh objectives are deliberately explicit: a finite Solver::Summary cost is not a
 * substitute for evaluating the retained parameter vector. */
struct ConvergenceAuditInput
{
  std::vector<double> successful_full_objectives;
  double fresh_full_objective = std::numeric_limits<double>::infinity();
  double fresh_data_objective = std::numeric_limits<double>::infinity();
  double projected_gradient = std::numeric_limits<double>::infinity();
  double relative_step = std::numeric_limits<double>::infinity();
  double trust_region_radius = std::numeric_limits<double>::infinity();
  bool parameters_finite_and_feasible = false;
};

enum class BudgetExhaustionDisposition
{
  Reject,
  UsableWithWarnings
};

struct ConvergenceAuditResult
{
  BudgetExhaustionDisposition disposition = BudgetExhaustionDisposition::Reject;
  bool usable_with_warnings = false;
  bool objectives_finite = false;
  bool history_finite = false;
  bool envelope_stagnant = false;
  bool tail_not_rising = false;
  bool stationary = false;
  bool trust_region_valid = false;
  bool step_or_radius_small = false;

  double envelope_relative_improvement = std::numeric_limits<double>::infinity();
  double tail_relative_rise = std::numeric_limits<double>::infinity();
};

/** Pure decision oracle for classifying a budget-exhausted solve.
 *
 * A cumulative best-cost envelope alone is flat during a sequence of increasingly worse accepted
 * nonmonotonic steps.  Therefore the envelope test is paired with a tail-rise guard.  Small
 * nonmonotonic numerical jitter is allowed up to `maximum_tail_rise`; a material rising tail is not
 * promoted even when the separately measured gradient and step happen to be small. */
inline ConvergenceAuditResult audit_budget_exhausted_solution(
    const ConvergenceAuditInput &input,
    const ConvergenceAuditThresholds &thresholds = ConvergenceAuditThresholds{} )
{
  ConvergenceAuditResult result;

  result.objectives_finite = std::isfinite(input.fresh_full_objective)
                          && std::isfinite(input.fresh_data_objective);
  result.history_finite = input.successful_full_objectives.size() >= 2;
  for( const double cost : input.successful_full_objectives )
    result.history_finite = result.history_finite && std::isfinite(cost) && (cost >= 0.0);

  if( result.history_finite )
  {
    const std::size_t count = input.successful_full_objectives.size();
    const std::size_t wanted_tail = (std::max)(std::size_t(2),thresholds.envelope_tail_length);
    const std::size_t tail_begin = (count > wanted_tail) ? (count - wanted_tail) : 0;

    double best_at_tail_begin = input.successful_full_objectives.front();
    for( std::size_t i = 1; i <= tail_begin; ++i )
      best_at_tail_begin = (std::min)(best_at_tail_begin,input.successful_full_objectives[i]);

    double best_final = best_at_tail_begin;
    double best_in_tail = input.successful_full_objectives[tail_begin];
    for( std::size_t i = tail_begin + 1; i < count; ++i )
    {
      best_final = (std::min)(best_final,input.successful_full_objectives[i]);
      best_in_tail = (std::min)(best_in_tail,input.successful_full_objectives[i]);
    }

    result.envelope_relative_improvement = (best_at_tail_begin - best_final)
        / (std::max)(1.0,std::fabs(best_at_tail_begin));
    result.tail_relative_rise = (input.successful_full_objectives.back() - best_in_tail)
        / (std::max)(1.0,std::fabs(best_in_tail));

    result.envelope_stagnant = std::isfinite(result.envelope_relative_improvement)
        && (result.envelope_relative_improvement >= 0.0)
        && (result.envelope_relative_improvement <= thresholds.stagnant_envelope_change);
    result.tail_not_rising = std::isfinite(result.tail_relative_rise)
        && (result.tail_relative_rise <= thresholds.maximum_tail_rise);
  }

  result.stationary = input.parameters_finite_and_feasible
                   && result.objectives_finite
                   && std::isfinite(input.projected_gradient)
                   && (input.projected_gradient <= thresholds.stationary_projected_gradient);
  result.trust_region_valid = std::isfinite(input.trust_region_radius)
                           && (input.trust_region_radius > 0.0);
  const bool relative_step_small = std::isfinite(input.relative_step)
                                && (input.relative_step <= thresholds.small_relative_step);
  const bool trust_region_collapsed = result.trust_region_valid
                                   && (input.trust_region_radius
                                       <= thresholds.collapsed_trust_region_radius);
  result.step_or_radius_small = relative_step_small || trust_region_collapsed;

  result.usable_with_warnings = result.history_finite
                             && result.envelope_stagnant
                             && result.tail_not_rising
                             && result.stationary
                             && result.trust_region_valid
                             && result.step_or_radius_small;
  result.disposition = result.usable_with_warnings
      ? BudgetExhaustionDisposition::UsableWithWarnings
      : BudgetExhaustionDisposition::Reject;
  return result;
}

} // namespace RelActCalcAutoImp

#endif // InterSpec_RelActCalcAuto_ConvergenceAudit_h
