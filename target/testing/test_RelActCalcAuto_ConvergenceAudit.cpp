/* InterSpec: focused RelActAuto budget-exhaustion classification tests. */

#include "InterSpec_config.h"

#include <limits>
#include <utility>
#include <vector>

#define BOOST_TEST_MODULE RelActCalcAuto_ConvergenceAudit_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/RelActCalcAuto_ConvergenceAudit.h"

using namespace RelActCalcAutoImp;

namespace
{
ConvergenceAuditInput otherwise_stationary_input( std::vector<double> history )
{
  ConvergenceAuditInput input;
  input.successful_full_objectives = std::move(history);
  input.fresh_full_objective = 90.0;
  input.fresh_data_objective = 87.0;
  input.projected_gradient = 2.0e-6;
  input.relative_step = 2.0e-7;
  input.trust_region_radius = 25.0;
  input.parameters_finite_and_feasible = true;
  return input;
}
} // namespace

BOOST_AUTO_TEST_CASE( rising_nonmonotonic_tail_is_rejected )
{
  // The global best (90) predates a rising accepted-step tail.  Its chronological best envelope is
  // therefore perfectly flat, reproducing the case that a signed oldest/newest test promoted.
  ConvergenceAuditInput input = otherwise_stationary_input(
      {90.0, 91.0, 100.0, 99.0, 101.0, 102.0, 103.0,
       104.0, 105.0, 106.0, 107.0, 108.0, 109.0} );
  input.fresh_full_objective = 109.0;

  const ConvergenceAuditResult result = audit_budget_exhausted_solution(input);
  BOOST_CHECK( result.envelope_stagnant );
  BOOST_CHECK( !result.tail_not_rising );
  BOOST_CHECK( result.stationary );
  BOOST_CHECK( result.step_or_radius_small );
  BOOST_CHECK( !result.usable_with_warnings );
  BOOST_CHECK( result.disposition == BudgetExhaustionDisposition::Reject );
}

BOOST_AUTO_TEST_CASE( flat_tail_with_large_projected_gradient_is_rejected )
{
  ConvergenceAuditInput input = otherwise_stationary_input(
      {90.0000004, 90.0000003, 90.0000002, 90.0000001, 90.0} );
  input.projected_gradient = 3.0e-2;

  const ConvergenceAuditResult result = audit_budget_exhausted_solution(input);
  BOOST_CHECK( result.envelope_stagnant );
  BOOST_CHECK( result.tail_not_rising );
  BOOST_CHECK( !result.stationary );
  BOOST_CHECK( !result.usable_with_warnings );
  BOOST_CHECK( result.disposition == BudgetExhaustionDisposition::Reject );
}

BOOST_AUTO_TEST_CASE( stationary_budget_limited_point_is_usable_with_warning )
{
  const ConvergenceAuditInput input = otherwise_stationary_input(
      {90.0000004, 90.0000003, 90.00000035, 90.0000002, 90.0000001, 90.0} );

  const ConvergenceAuditResult result = audit_budget_exhausted_solution(input);
  BOOST_CHECK( result.objectives_finite );
  BOOST_CHECK( result.history_finite );
  BOOST_CHECK( result.envelope_stagnant );
  BOOST_CHECK( result.tail_not_rising );
  BOOST_CHECK( result.stationary );
  BOOST_CHECK( result.trust_region_valid );
  BOOST_CHECK( result.step_or_radius_small );
  BOOST_CHECK( result.usable_with_warnings );
  BOOST_CHECK( result.disposition == BudgetExhaustionDisposition::UsableWithWarnings );
}

BOOST_AUTO_TEST_CASE( collapsed_trust_region_can_document_a_stationary_small_move )
{
  ConvergenceAuditInput input = otherwise_stationary_input(
      {90.0000002, 90.0000001, 90.0} );
  input.relative_step = 1.0e-2; //not independently small
  input.trust_region_radius = 1.0e-14;

  const ConvergenceAuditResult result = audit_budget_exhausted_solution(input);
  BOOST_CHECK( result.stationary );
  BOOST_CHECK( result.step_or_radius_small );
  BOOST_CHECK( result.disposition == BudgetExhaustionDisposition::UsableWithWarnings );
}

BOOST_AUTO_TEST_CASE( fresh_objective_or_feasibility_failure_is_rejected )
{
  ConvergenceAuditInput input = otherwise_stationary_input({90.0000001,90.0});
  input.fresh_data_objective = std::numeric_limits<double>::infinity();
  BOOST_CHECK( audit_budget_exhausted_solution(input).disposition
               == BudgetExhaustionDisposition::Reject );

  input.fresh_data_objective = 87.0;
  input.parameters_finite_and_feasible = false;
  BOOST_CHECK( audit_budget_exhausted_solution(input).disposition
               == BudgetExhaustionDisposition::Reject );
}
