/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; either version 2.1 of the License, or (at your option)
 any later version.
 */

#include "InterSpec_config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE RelActCalcAuto_ProfileScan_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/RelActCalcAuto_ProfileScan.h"

namespace PL = RelActCalcAutoImp::ProfileLikelihood;

namespace
{
constexpr double sm_threshold95 = 3.841458820694124;

PL::Evaluation success( const double reported, const double delta )
{
  PL::Evaluation result;
  result.status = PL::EvaluationStatus::Success;
  result.reported_fraction = reported;
  result.delta_chi2 = delta;
  return result;
}
}//namespace


BOOST_AUTO_TEST_CASE( quadratic_profile_has_symmetric_68_and_95_crossings )
{
  const auto result = PL::scan(0.4,0.4,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double x, const PL::Direction, const std::size_t ) {
      return success(x,std::pow((x-0.4)/0.1,2.0));
    });

  BOOST_REQUIRE( result.status == PL::ScanStatus::Complete );
  BOOST_REQUIRE_LE( result.num_fits,32u );
  BOOST_CHECK( result.intervals[0].lower.likelihood_crossing );
  BOOST_CHECK( result.intervals[0].upper.likelihood_crossing );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[0].lower.reported_fraction,0.3,0.01 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[0].upper.reported_fraction,0.5,0.01 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[1].lower.reported_fraction,
                             0.4-0.1*std::sqrt(sm_threshold95),0.01 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[1].upper.reported_fraction,
                             0.4+0.1*std::sqrt(sm_threshold95),0.01 );
}


BOOST_AUTO_TEST_CASE( asymmetric_profile_and_reported_coordinate_mapping )
{
  // The objective is narrower below the optimum.  The affine output mapping stands in for a
  // post-correlation reported fraction whose scan control is a fitted pre-correction coordinate.
  const auto result = PL::scan(0.3,0.34,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double x, const PL::Direction, const std::size_t ) {
      const double sigma = (x < 0.3) ? 0.05 : 0.15;
      return success(0.1+0.8*x,std::pow((x-0.3)/sigma,2.0));
    });

  BOOST_REQUIRE( result.status == PL::ScanStatus::Complete );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[0].lower.reported_fraction,
                             0.1+0.8*0.25,0.015 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[0].upper.reported_fraction,
                             0.1+0.8*0.45,0.015 );
  BOOST_CHECK_LT( 0.34-result.intervals[0].lower.reported_fraction,
                  result.intervals[0].upper.reported_fraction-0.34 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[1].lower.reported_fraction,
                             0.1+0.8*(0.3-0.05*std::sqrt(sm_threshold95)),0.015 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[1].upper.reported_fraction,
                             0.1+0.8*(0.3+0.15*std::sqrt(sm_threshold95)),0.015 );
}


BOOST_AUTO_TEST_CASE( nonchronological_profile_requests_choose_nearest_warm_point )
{
  // Mimic geometric bracketing followed by Brent and then the nested 68% pass: the most recently
  // evaluated point (0.18) is not the best warm start for a request back near the baseline.
  const std::vector<double> cached_controls{0.4,0.3,0.1,0.18};
  BOOST_CHECK_EQUAL( PL::nearest_control_index(cached_controls,0.34),1u );
  BOOST_CHECK_EQUAL( PL::nearest_control_index(cached_controls,0.12),2u );
  BOOST_CHECK_EQUAL( PL::nearest_control_index(cached_controls,0.19),3u );

  const std::vector<double> with_invalid{
      std::numeric_limits<double>::quiet_NaN(),0.25,0.75};
  BOOST_CHECK_EQUAL( PL::nearest_control_index(with_invalid,0.5),1u ); // stable tie
  BOOST_CHECK_EQUAL( PL::nearest_control_index(with_invalid,
                     std::numeric_limits<double>::quiet_NaN()),with_invalid.size() );
}


BOOST_AUTO_TEST_CASE( pending_baseline_selection_is_independent_of_discovery_order )
{
  const std::array<PL::PendingBaselineDiscovery,4> candidates{{
      {3045.0,"0:Pu239"},
      {3044.75,"0:Pu240"},
      {3044.5,"0:Pu238"},
      {3046.0,"0:Pu241"}
  }};
  std::array<std::size_t,4> order{{0,1,2,3}};
  do
  {
    std::vector<PL::PendingBaselineDiscovery> permuted;
    for( const std::size_t index : order )
      permuted.push_back(candidates[index]);
    const auto best = PL::best_pending_baseline_discovery_index(permuted);
    BOOST_REQUIRE( best );
    BOOST_CHECK_EQUAL( permuted[*best].semantic_key,"0:Pu238" );
    BOOST_CHECK_EQUAL( permuted[*best].full_objective,3044.5 );
  }while( std::next_permutation(order.begin(),order.end()) );
}


BOOST_AUTO_TEST_CASE( baseline_discoveries_are_deferred_for_exactly_one_profile_pass )
{
  using Disposition = PL::BaselineDiscoveryDisposition;
  BOOST_CHECK( PL::baseline_discovery_disposition(0)
               == Disposition::DeferUntilPassComplete );
  BOOST_CHECK( PL::baseline_discovery_disposition(1)
               == Disposition::RejectAfterRestart );
  BOOST_CHECK( PL::baseline_discovery_disposition(2)
               == Disposition::RejectAfterRestart );

  // Mimic the production first pass: every eligible target can contribute before selection.  The
  // lowest later discovery must win even though an earlier one already improved the baseline.
  const std::vector<PL::PendingBaselineDiscovery> completed_first_pass{
      {3045.0,"curve-a|target=Pu239"},
      {3044.5,"curve-a|target=Pu238"},
      {3044.8,"curve-a|target=Pu242"}
  };
  const auto selected = PL::best_pending_baseline_discovery_index(completed_first_pass);
  BOOST_REQUIRE( selected );
  BOOST_CHECK_EQUAL( *selected,1u );
}


BOOST_AUTO_TEST_CASE( pending_baseline_selection_uses_exact_semantic_tie_break )
{
  const double objective = 3045.177;
  const double strictly_lower = std::nextafter(objective,
                                  -std::numeric_limits<double>::infinity());
  std::vector<PL::PendingBaselineDiscovery> candidates{
      {objective,"0:Pu241"},
      {objective,"0:Pu238"},
      {strictly_lower,"0:Pu242"},
      {std::numeric_limits<double>::quiet_NaN(),"0:Pu000"},
      {-std::numeric_limits<double>::infinity(),"0:Pu001"}
  };

  auto best = PL::best_pending_baseline_discovery_index(candidates);
  BOOST_REQUIRE( best );
  // No objective tolerance is applied: even one representable step lower wins first.
  BOOST_CHECK_EQUAL( candidates[*best].semantic_key,"0:Pu242" );

  candidates.erase(candidates.begin()+2);
  best = PL::best_pending_baseline_discovery_index(candidates);
  BOOST_REQUIRE( best );
  BOOST_CHECK_EQUAL( candidates[*best].semantic_key,"0:Pu238" );

  std::reverse(candidates.begin(),candidates.end());
  best = PL::best_pending_baseline_discovery_index(candidates);
  BOOST_REQUIRE( best );
  BOOST_CHECK_EQUAL( candidates[*best].semantic_key,"0:Pu238" );

  const std::vector<PL::PendingBaselineDiscovery> no_finite{
      {std::numeric_limits<double>::infinity(),"0:Pu238"},
      {std::numeric_limits<double>::quiet_NaN(),"0:Pu239"}
  };
  BOOST_CHECK( !PL::best_pending_baseline_discovery_index(no_finite) );
}


BOOST_AUTO_TEST_CASE( better_baseline_gate_clears_optimizer_noise_but_not_a_real_basin )
{
  // Values measured on the JRC CBNM Pu70 free-age 610-775 keV fit: an ill-conditioned (kappa~1.7e7)
  // rank-deficient solve whose conditional fits resample the same flat direction.
  const double objective = 3045.1595741668584;
  const double cov_scale = 4.4068879510374215;
  const double tolerance = PL::baseline_improvement_tolerance(objective,cov_scale);

  // Both improvements the profile pass actually reported are round-off on this problem, and both
  // must stay below the gate: 5.7e-6 and 1.0e-6 of the objective respectively.
  BOOST_CHECK_GT( tolerance,0.017433 );   //first pass, Pu240
  BOOST_CHECK_GT( tolerance,0.003164 );   //post-restart, Pu238

  // ...yet the gate stays far below anything that could move a reported endpoint.  Interval
  // endpoints are crossings of cov_scale*1 and cov_scale*3.841458820694124.
  BOOST_CHECK_LT( tolerance,0.01*cov_scale );
  // ...and below the cost change solve_ceres documents as achievable (about 0.1 in chi-square), so
  // the gate remains conservative rather than merely permissive.
  BOOST_CHECK_LT( tolerance,0.1 );

  // A genuinely different basin clears it by orders of magnitude.
  BOOST_CHECK_LT( tolerance,1.0 );

  // Small-objective problems keep a usable absolute floor, and the gate is never zero or negative.
  BOOST_CHECK_GE( PL::baseline_improvement_tolerance(0.0,0.0),1.0e-6 );
  BOOST_CHECK_GT( PL::baseline_improvement_tolerance(0.0,0.0),0.0 );
  // Sign of the objective is irrelevant; only its magnitude sets the scale.
  BOOST_CHECK_EQUAL( PL::baseline_improvement_tolerance(-objective,cov_scale),tolerance );
  // A non-finite input must not produce a non-finite gate, which would disable the guard entirely.
  const double nan = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK( std::isfinite(PL::baseline_improvement_tolerance(nan,cov_scale)) );
  BOOST_CHECK( std::isfinite(PL::baseline_improvement_tolerance(objective,nan)) );

  // The interval-threshold term governs when an overdispersed fit needs a coarser gate than its
  // objective magnitude alone implies.
  BOOST_CHECK_EQUAL( PL::baseline_improvement_tolerance(1.0,1.0e4),10.0 );
  // The objective term governs the ordinary case.
  BOOST_CHECK_CLOSE_FRACTION( PL::baseline_improvement_tolerance(objective,1.0),
                              1.0e-5*(1.0 + objective),1.0e-12 );
}


BOOST_AUTO_TEST_CASE( deferred_discoveries_order_best_first_independently_of_caller_order )
{
  // One warm seed cannot cover basins that different profile targets entered independently, so the
  // single permitted reselection is offered several deterministic seeds.  The ordering must agree
  // with the single-best selector on its first element, and must not depend on caller order.
  const std::array<PL::PendingBaselineDiscovery,5> candidates{{
      {3045.0,"curve-a|target=Pu239"},
      {3044.5,"curve-a|target=Pu238"},
      {std::numeric_limits<double>::quiet_NaN(),"curve-a|target=Pu241"},
      {3044.8,"curve-a|target=Pu242"},
      {-std::numeric_limits<double>::infinity(),"curve-a|target=Pu240"}
  }};

  std::array<std::size_t,5> permutation{{0,1,2,3,4}};
  do
  {
    std::vector<PL::PendingBaselineDiscovery> permuted;
    for( const std::size_t index : permutation )
      permuted.push_back(candidates[index]);

    const std::vector<std::size_t> order = PL::ordered_pending_baseline_discoveries(permuted);
    // Both non-finite entries are dropped: neither can seed anything.
    BOOST_REQUIRE_EQUAL( order.size(),3u );
    BOOST_CHECK_EQUAL( permuted[order[0]].semantic_key,"curve-a|target=Pu238" );
    BOOST_CHECK_EQUAL( permuted[order[1]].semantic_key,"curve-a|target=Pu242" );
    BOOST_CHECK_EQUAL( permuted[order[2]].semantic_key,"curve-a|target=Pu239" );

    // The many-seed and one-seed paths must never disagree about the winner.
    const auto best = PL::best_pending_baseline_discovery_index(permuted);
    BOOST_REQUIRE( best );
    BOOST_CHECK_EQUAL( order.front(),*best );

    for( std::size_t i = 1; i < order.size(); ++i )
      BOOST_CHECK_LE( permuted[order[i-1]].full_objective,permuted[order[i]].full_objective );
  }while( std::next_permutation(permutation.begin(),permutation.end()) );

  BOOST_CHECK( PL::ordered_pending_baseline_discoveries({}).empty() );
}


BOOST_AUTO_TEST_CASE( deferred_discovery_ordering_breaks_exact_ties_semantically )
{
  const double objective = 3045.177;
  const double strictly_lower = std::nextafter(objective,
                                  -std::numeric_limits<double>::infinity());
  const std::vector<PL::PendingBaselineDiscovery> candidates{
      {objective,"0:Pu241"},
      {objective,"0:Pu238"},
      {strictly_lower,"0:Pu242"}
  };
  const std::vector<std::size_t> order = PL::ordered_pending_baseline_discoveries(candidates);
  BOOST_REQUIRE_EQUAL( order.size(),3u );
  // One representable step lower still beats every semantic key.
  BOOST_CHECK_EQUAL( candidates[order[0]].semantic_key,"0:Pu242" );
  BOOST_CHECK_EQUAL( candidates[order[1]].semantic_key,"0:Pu238" );
  BOOST_CHECK_EQUAL( candidates[order[2]].semantic_key,"0:Pu241" );
}


BOOST_AUTO_TEST_CASE( deferred_discoveries_deduplicate_by_semantic_key_keeping_the_best )
{
  // Two profile targets routinely fall into the same basin; seeding the reselection twice from it
  // wastes a bounded candidate slot that a genuinely different basin could have used.
  const std::vector<PL::PendingBaselineDiscovery> candidates{
      {3045.0,"curve-a|target=Pu239"},
      {3044.5,"curve-a|target=Pu239"},
      {3044.8,"curve-b|target=Pu239"},
      {3044.9,"curve-a|target=Pu239"}
  };
  const std::vector<std::size_t> order = PL::ordered_pending_baseline_discoveries(candidates);
  const std::vector<std::size_t> unique
      = PL::unique_pending_baseline_discoveries(candidates,order);

  BOOST_REQUIRE_EQUAL( order.size(),4u );
  BOOST_REQUIRE_EQUAL( unique.size(),2u );
  // The retained entry for a repeated key is that key's lowest objective, and overall order holds.
  BOOST_CHECK_EQUAL( unique[0],1u );
  BOOST_CHECK_EQUAL( candidates[unique[0]].full_objective,3044.5 );
  BOOST_CHECK_EQUAL( candidates[unique[1]].semantic_key,"curve-b|target=Pu239" );

  BOOST_CHECK( PL::unique_pending_baseline_discoveries(candidates,{}).empty() );
}


BOOST_AUTO_TEST_CASE( brent_dekker_refines_nonlinear_asymmetric_crossings )
{
  std::size_t calls = 0;
  const auto result = PL::scan(0.45,0.45,0.0,1.0,{{1.0,sm_threshold95}},32,
    [&calls]( const double x, const PL::Direction, const std::size_t remaining ) {
      ++calls;
      BOOST_REQUIRE_GT( remaining,0u );
      BOOST_REQUIRE_GE( x,0.0 );
      BOOST_REQUIRE_LE( x,1.0 );
      const double scale = (x < 0.45) ? 0.04 : 0.08;
      return success(x,std::expm1(std::fabs(x-0.45)/scale));
    });

  BOOST_REQUIRE( result.status == PL::ScanStatus::Complete );
  BOOST_CHECK_EQUAL( calls,result.num_fits );
  BOOST_CHECK_LE( result.num_fits,32u );
  const std::array<double,2> expected_lower{{
      0.45-0.04*std::log1p(1.0),
      0.45-0.04*std::log1p(sm_threshold95)
  }};
  const std::array<double,2> expected_upper{{
      0.45+0.08*std::log1p(1.0),
      0.45+0.08*std::log1p(sm_threshold95)
  }};
  for( std::size_t level = 0; level < 2; ++level )
  {
    BOOST_CHECK_SMALL( result.intervals[level].lower.reported_fraction
                         - expected_lower[level],2.0e-4 );
    BOOST_CHECK_SMALL( result.intervals[level].upper.reported_fraction
                         - expected_upper[level],2.0e-4 );
  }
}


BOOST_AUTO_TEST_CASE( steep_crossing_rejects_endpoint_hugging_interpolation )
{
  // A flat likelihood followed by a sub-resolution rise is the limiting case seen in the Pu-241
  // upper profile.  Unsafeguarded secant steps alternate between a negligible move off the flat
  // endpoint and a bisection, exhausting 32 fits before all four nested endpoints are classified.
  constexpr double baseline = 0.4;
  constexpr double steep_crossing = 0.7;
  std::size_t calls = 0;
  const auto result = PL::scan(baseline,baseline,0.0,1.0,{{1.0,sm_threshold95}},32,
    [&calls]( const double x, const PL::Direction direction, const std::size_t remaining ) {
      ++calls;
      BOOST_REQUIRE_GT( remaining,0u );
      const double delta = (direction == PL::Direction::Lower)
          ? std::pow((x-baseline)/0.05,2.0)
          : ((x < steep_crossing) ? 0.4 : 10000.0);
      return success(x,delta);
    });

  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::Complete,
                         result.diagnostic << ", fits=" << result.num_fits );
  BOOST_CHECK_EQUAL( calls,result.num_fits );
  BOOST_CHECK_LE( result.num_fits,32u );
  BOOST_CHECK_SMALL( result.intervals[0].upper.reported_fraction-steep_crossing,2.0e-5 );
  BOOST_CHECK_SMALL( result.intervals[1].upper.reported_fraction-steep_crossing,2.0e-5 );
  BOOST_CHECK( result.intervals[0].lower.likelihood_crossing );
  BOOST_CHECK( result.intervals[1].lower.likelihood_crossing );
}


BOOST_AUTO_TEST_CASE( trace_profile_borrows_global_budget_and_refines_nested_crossings )
{
  // This reproduces the geometry and fit accounting of the Pu-241 failure which originally
  // stopped at 16/32 fits.  The lower direction is a trace-scale steep profile; the upper side has
  // a broad shelf followed by a very steep rise.  The first difficult solve on each side consumes
  // three augmented-Lagrangian fits, while subsequent warm solves consume one.
  constexpr double baseline = 9.2053390775993291e-6;
  constexpr double first_lower = 8.6383401282942523e-6;
  const double lower_sigma = (baseline-first_lower)/std::sqrt(87.5044656717522);
  constexpr double upper_plateau = 0.438;
  constexpr double upper_rise_scale = 0.0114;
  constexpr double upper_rise_rate = 82.4;
  constexpr double structural_upper = 0.43796841532903286;
  bool charged_lower_al = false;
  bool charged_upper_al = false;

  const auto upper_delta = [&]( const double x ) {
    const double shelf = upper_plateau
                       * (1.0-std::exp(-(x-baseline)/0.01));
    const double rise = (x > 0.25)
        ? upper_rise_scale*std::expm1(upper_rise_rate*(x-0.25)) : 0.0;
    return shelf+rise;
  };

  const auto result = PL::scan(baseline,baseline,0.0,1.0,
      {{1.2575068653963255,4.83066}},32,
    [&]( const double requested, const PL::Direction direction, const std::size_t remaining ) {
      BOOST_REQUIRE_GT( remaining,0u );
      PL::Evaluation evaluation;
      evaluation.status = PL::EvaluationStatus::Success;
      evaluation.reported_fraction = requested;
      if( direction == PL::Direction::Lower )
      {
        evaluation.delta_chi2 = std::pow((baseline-requested)/lower_sigma,2.0);
        // Match the finite accuracy of the reported-coordinate augmented-Lagrangian solve.  Root
        // refinement should use this scale, but must not use 1e-5 of the full [0,1] interval.
        evaluation.control_tolerance = 8.0e-9;
        if( !charged_lower_al )
        {
          evaluation.num_fits = 3;
          charged_lower_al = true;
        }
      }else
      {
        const double evaluated = (std::min)(requested,structural_upper);
        evaluation.reported_fraction = evaluated;
        evaluation.delta_chi2 = upper_delta(evaluated);
        evaluation.control_tolerance = 1.0e-8;
        if( requested > structural_upper )
        {
          evaluation.reached_feasible_bound = true;
          if( !charged_upper_al )
          {
            evaluation.num_fits = 3;
            charged_upper_al = true;
          }
        }
      }
      return evaluation;
    });

  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::Complete,
                         result.diagnostic << ", fits=" << result.num_fits );
  BOOST_CHECK_LE( result.num_fits,32u );
  const double expected_lower68
      = baseline-lower_sigma*std::sqrt(result.intervals[0].delta_chi2);
  const double expected_lower95
      = baseline-lower_sigma*std::sqrt(result.intervals[1].delta_chi2);
  BOOST_CHECK_SMALL( result.intervals[0].lower.reported_fraction-expected_lower68,1.7e-8 );
  BOOST_CHECK_SMALL( result.intervals[1].lower.reported_fraction-expected_lower95,1.7e-8 );
  BOOST_CHECK_LT( result.intervals[1].lower.reported_fraction,
                  result.intervals[0].lower.reported_fraction );
  BOOST_CHECK_GT( result.intervals[0].upper.reported_fraction,0.25 );
  BOOST_CHECK_GT( result.intervals[1].upper.reported_fraction,
                  result.intervals[0].upper.reported_fraction );
}


BOOST_AUTO_TEST_CASE( physical_boundary_is_an_endpoint_below_threshold )
{
  const auto result = PL::scan(0.05,0.05,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double x, const PL::Direction, const std::size_t ) {
      return success(x,std::pow((x-0.05)/0.1,2.0));
    });

  BOOST_REQUIRE( result.status == PL::ScanStatus::BoundaryLimited );
  BOOST_CHECK_SMALL( result.intervals[0].lower.reported_fraction,1.0e-12 );
  BOOST_CHECK( !result.intervals[0].lower.likelihood_crossing );
  BOOST_CHECK( result.intervals[0].upper.likelihood_crossing );
  BOOST_CHECK_SMALL( result.intervals[1].lower.reported_fraction,1.0e-12 );
  BOOST_CHECK( !result.intervals[1].lower.likelihood_crossing );
}


BOOST_AUTO_TEST_CASE( entire_feasible_range_below_95_is_non_identifiable )
{
  const auto result = PL::scan(0.0,0.0,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double x, const PL::Direction, const std::size_t ) {
      return success(x,0.1*x*x);
    });

  BOOST_REQUIRE( result.status == PL::ScanStatus::NonIdentifiable );
  BOOST_CHECK_SMALL( result.intervals[1].lower.reported_fraction,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[1].upper.reported_fraction,1.0,1.0e-12 );
  BOOST_CHECK( !result.intervals[1].lower.likelihood_crossing );
  BOOST_CHECK( !result.intervals[1].upper.likelihood_crossing );
}


BOOST_AUTO_TEST_CASE( transformed_structural_bounds_classify_full_range )
{
  const auto result = PL::scan(0.4,0.4,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double x, const PL::Direction, const std::size_t ) {
      PL::Evaluation evaluation = success((std::max)(0.2,(std::min)(0.8,x)),0.1);
      evaluation.reached_feasible_bound = (x < 0.2) || (x > 0.8);
      return evaluation;
    });

  BOOST_REQUIRE( result.status == PL::ScanStatus::NonIdentifiable );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[1].lower.reported_fraction,0.2,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( result.intervals[1].upper.reported_fraction,0.8,1.0e-12 );
  BOOST_CHECK( !result.intervals[1].lower.likelihood_crossing );
  BOOST_CHECK( !result.intervals[1].upper.likelihood_crossing );
}


BOOST_AUTO_TEST_CASE( stable_augmented_lagrangian_limit_can_be_interior )
{
  constexpr double requested = 0.50000460267;
  constexpr double actual = 0.437968415302;
  constexpr double baseline = 9.2053390776e-6;
  constexpr double tolerance = 1.99998158932e-5;
  BOOST_CHECK( PL::stable_reported_structural_bound(
      requested,actual,baseline,1.11e-16,tolerance,1.0-baseline,3,
      PL::Direction::Upper) );
  BOOST_CHECK( !PL::stable_reported_structural_bound(
      requested,actual,baseline,1.0e-2,tolerance,1.0-baseline,3,
      PL::Direction::Upper) );
  BOOST_CHECK( !PL::stable_reported_structural_bound(
      requested,0.6,baseline,1.11e-16,tolerance,1.0-baseline,3,
      PL::Direction::Upper) );
  BOOST_CHECK( !PL::stable_reported_structural_bound(
      requested,actual,baseline,1.11e-16,tolerance,1.0-baseline,1,
      PL::Direction::Upper) );
}


BOOST_AUTO_TEST_CASE( fit_cap_failure_is_explicit_and_never_exceeded )
{
  std::size_t callback_count = 0;
  const auto result = PL::scan(0.5,0.5,0.0,1.0,{{1.0,sm_threshold95}},1,
    [&callback_count]( const double x, const PL::Direction, const std::size_t ) {
      ++callback_count;
      return success(x,std::pow((x-0.5)/0.1,2.0));
    });

  BOOST_CHECK( result.status == PL::ScanStatus::FitCapReached );
  BOOST_CHECK_EQUAL( result.num_fits,1u );
  BOOST_CHECK_EQUAL( callback_count,1u );
}


BOOST_AUTO_TEST_CASE( fit_cap_during_brent_is_not_reported_as_a_complete_interval )
{
  const auto result = PL::scan(0.4,0.4,0.0,1.0,{{1.0,sm_threshold95}},9,
    []( const double x, const PL::Direction, const std::size_t ) {
      return success(x,std::pow((x-0.4)/0.1,2.0));
    });
  BOOST_CHECK( result.status == PL::ScanStatus::FitCapReached );
  BOOST_CHECK_EQUAL( result.num_fits,9u );
}


BOOST_AUTO_TEST_CASE( terminal_fit_cap_uses_only_a_well_localized_sign_bracket )
{
  // Exact final Pu-239/68% bracket from the Pu70 free-age profile.  All 32 conditional fits had
  // been spent, but the two independently evaluated endpoints already localized the crossing to
  // 1.96e-5 in control space and straddled the threshold by about 0.3%.
  constexpr double baseline = 0.64620916490723379;
  constexpr double threshold = 4.4133000106711426;
  const PL::Sample inside{
      0.79698045369919035,0.79698045369919035,4.3995358928086716,false,0.0};
  const PL::Sample outside{
      0.79700009881129485,0.79700009881129485,4.4265069835273607,false,0.0};
  const double directional_span = 1.0-baseline;
  const double normal_control_tolerance = 1.0e-5*directional_span;

  const auto accepted = PL::detail::terminal_cap_crossing_sample(
      inside,outside,threshold,normal_control_tolerance,directional_span,true,true);
  BOOST_REQUIRE( accepted );
  // The outside endpoint is marginally closer to the requested likelihood threshold.
  BOOST_CHECK_EQUAL( accepted->control_fraction,outside.control_fraction );
  BOOST_CHECK_EQUAL( accepted->delta_chi2,outside.delta_chi2 );

  // Even when the inside residual is closer, returning it would truncate values independently
  // known to remain in the confidence set.  A terminal endpoint is always conservative.
  PL::Sample closer_inside = inside;
  closer_inside.delta_chi2 = threshold*0.999;
  PL::Sample farther_outside = outside;
  farther_outside.delta_chi2 = threshold*1.004;
  const auto conservative = PL::detail::terminal_cap_crossing_sample(
      closer_inside,farther_outside,threshold,normal_control_tolerance,
      directional_span,true,true);
  BOOST_REQUIRE( conservative );
  BOOST_CHECK_EQUAL( conservative->control_fraction,farther_outside.control_fraction );
  BOOST_CHECK_GE( conservative->delta_chi2,threshold );

  // This relaxation is unavailable before the true global cap and while another bracketed
  // endpoint remains to be refined.
  BOOST_CHECK( !PL::detail::terminal_cap_crossing_sample(
      inside,outside,threshold,normal_control_tolerance,directional_span,false,true) );
  BOOST_CHECK( !PL::detail::terminal_cap_crossing_sample(
      inside,outside,threshold,normal_control_tolerance,directional_span,true,false) );

  // A broad bracket remains an explicit fit-cap failure unless one independently evaluated end
  // is already within the separate 0.5%-of-threshold terminal objective tolerance.
  const PL::Sample broad_inside{0.70,0.70,1.0,false,0.0};
  const PL::Sample broad_outside{0.85,0.85,20.0,false,0.0};
  BOOST_CHECK( !PL::detail::terminal_cap_crossing_sample(
      broad_inside,broad_outside,threshold,normal_control_tolerance,
      directional_span,true,true) );

  PL::Sample broad_inside_close = broad_inside;
  broad_inside_close.delta_chi2 = threshold*0.999;
  BOOST_CHECK( !PL::detail::terminal_cap_crossing_sample(
      broad_inside_close,broad_outside,threshold,normal_control_tolerance,
      directional_span,true,true) );

  PL::Sample close_outside = broad_outside;
  close_outside.delta_chi2 = threshold*1.004;
  const auto objective_close = PL::detail::terminal_cap_crossing_sample(
      broad_inside,close_outside,threshold,normal_control_tolerance,
      directional_span,true,true);
  BOOST_REQUIRE( objective_close );
  BOOST_CHECK_EQUAL( objective_close->control_fraction,close_outside.control_fraction );

  close_outside.delta_chi2 = threshold*1.006;
  BOOST_CHECK( !PL::detail::terminal_cap_crossing_sample(
      broad_inside,close_outside,threshold,normal_control_tolerance,
      directional_span,true,true) );
}


BOOST_AUTO_TEST_CASE( final_nested_crossing_reuses_a_tight_bracket_at_exact_fit_cap )
{
  // The outer threshold tightens a discontinuous upper bracket.  Its last inside evaluation uses
  // the remaining conditional-fit allocation exactly; the nested 68% endpoint must reuse that
  // inherited sign bracket without another callback.  This exercises the pre-Brent no-budget
  // path, not merely the terminal fallback after a Brent iteration.
  std::size_t callback_count = 0;
  const auto result = PL::scan(0.0,0.0,0.0,1.0,{{1.0,sm_threshold95}},32,
    [&callback_count]( const double x, const PL::Direction, const std::size_t remaining ) {
      ++callback_count;
      PL::Evaluation evaluation = success(x,(x < 0.5) ? 0.0 : 10000.0);
      if( (x > 0.49999) && (x < 0.5) )
        evaluation.num_fits = remaining;
      return evaluation;
    });

  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::BoundaryLimited,
                         result.diagnostic << ", fits=" << result.num_fits );
  BOOST_CHECK_EQUAL( result.num_fits,32u );
  BOOST_CHECK_LT( callback_count,result.num_fits );
  BOOST_CHECK_NE( result.diagnostic.find("terminal sign bracket"),std::string::npos );
  BOOST_CHECK( result.intervals[0].upper.likelihood_crossing );
  BOOST_CHECK_EQUAL( result.intervals[0].upper.reported_fraction,0.5 );
  BOOST_CHECK_GE( result.intervals[0].upper.reported_fraction,0.5 );
  BOOST_CHECK( !result.intervals[0].lower.likelihood_crossing );
  BOOST_CHECK_EQUAL( result.intervals[0].lower.reported_fraction,0.0 );
}


BOOST_AUTO_TEST_CASE( brent_callback_receives_only_the_remaining_global_fit_allocation )
{
  std::vector<std::size_t> remaining_seen;
  const auto result = PL::scan(0.4,0.4,0.0,1.0,{{1.0,sm_threshold95}},32,
    [&remaining_seen]( const double x, const PL::Direction, const std::size_t remaining ) {
      remaining_seen.push_back(remaining);
      PL::Evaluation evaluation = success(x,std::pow((x-0.4)/0.1,2.0));
      // Deliberately consume the complete global remainder on the first Brent call.  This is a
      // pathological evaluator, but it verifies that the scanner accounts the hard cap exactly.
      if( remaining_seen.size() == 9 )
        evaluation.num_fits = remaining;
      return evaluation;
    });

  BOOST_REQUIRE_GE( remaining_seen.size(),9u );
  // Eight bracket fits leave the true global remainder of 24.  No local endpoint allocation may
  // masquerade as the hard cap; if the evaluator spends all 24, FitCapReached is reported at 32.
  BOOST_CHECK_EQUAL( remaining_seen[8],24u );
  BOOST_CHECK( result.status == PL::ScanStatus::FitCapReached );
  BOOST_CHECK_EQUAL( result.num_fits,32u );
}


BOOST_AUTO_TEST_CASE( sibling_windows_project_to_exact_simplex_bounds )
{
  const auto bounded = PL::project_simplex_component_bounds(
      0.0,1.0,{{0.10,0.25},{0.20,0.35},{0.05,0.10}} );
  BOOST_CHECK_CLOSE_FRACTION( bounded.first,0.30,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( bounded.second,0.65,1.0e-12 );

  const auto own_window = PL::project_simplex_component_bounds(
      0.4,0.6,{{0.0,1.0}} );
  BOOST_CHECK_CLOSE_FRACTION( own_window.first,0.4,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( own_window.second,0.6,1.0e-12 );
}


BOOST_AUTO_TEST_CASE( automatic_weak_trigger_uses_projected_feasible_interval )
{
  // This local 95% band is far too narrow to span the generic physical [0,1] domain, but it
  // covers the complete input-constrained interval and must therefore trigger profiling.
  BOOST_CHECK( !PL::local_gaussian_band_spans_feasible_range(0.50,0.04,0.0,1.0) );
  BOOST_CHECK( PL::local_gaussian_band_spans_feasible_range(0.50,0.04,0.45,0.55) );

  const auto projected = PL::project_simplex_component_bounds(
      0.0,1.0,{{0.20,0.25},{0.25,0.30}} );
  BOOST_CHECK_CLOSE_FRACTION( projected.first,0.45,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( projected.second,0.55,1.0e-12 );
  BOOST_CHECK( PL::local_gaussian_band_spans_feasible_range(
      0.50,0.04,projected.first,projected.second) );

  BOOST_CHECK( !PL::local_gaussian_band_spans_feasible_range(
      std::numeric_limits<double>::quiet_NaN(),0.04,0.45,0.55) );
  BOOST_CHECK( !PL::local_gaussian_band_spans_feasible_range(0.50,-0.04,0.45,0.55) );
  BOOST_CHECK( !PL::local_gaussian_band_spans_feasible_range(0.50,0.04,0.55,0.45) );
}


BOOST_AUTO_TEST_CASE( activity_boxes_and_fixed_ratio_groups_project_exact_fraction_bounds )
{
  // Target activity in [2,4], competing activity in [6,8].  With equal mass/activity
  // coefficients, the exact corner projection is [2/(2+8),4/(4+6)].
  const std::vector<PL::ActivityBoxMassGroup> boxed{
      {1.0,1.0,2.0,4.0},
      {1.0,0.0,6.0,8.0}
  };
  const auto bounds = PL::activity_box_fraction_bounds(boxed);
  BOOST_REQUIRE( bounds );
  BOOST_CHECK_CLOSE_FRACTION( bounds->first,0.20,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( bounds->second,0.40,1.0e-12 );

  // A target tied by activity to another isotope in the same sole root group has a fixed 1/2
  // within-group composition, independent of the roots activity box.
  const auto ratio_fixed = PL::activity_box_fraction_bounds({{2.0,1.0,3.0,7.0}});
  BOOST_REQUIRE( ratio_fixed );
  BOOST_CHECK_CLOSE_FRACTION( ratio_fixed->first,0.5,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( ratio_fixed->second,0.5,1.0e-12 );

  const auto zero_competitor = PL::activity_box_fraction_bounds({
      {1.0,1.0,0.0,std::numeric_limits<double>::infinity()},
      {1.0,0.0,0.0,0.0}
  });
  BOOST_REQUIRE( zero_competitor );
  BOOST_CHECK_CLOSE_FRACTION( zero_competitor->first,1.0,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( zero_competitor->second,1.0,1.0e-12 );

  const auto zero_target = PL::activity_box_fraction_bounds({
      {1.0,1.0,0.0,0.0}, {1.0,0.0,0.0,4.0}
  });
  BOOST_REQUIRE( zero_target );
  BOOST_CHECK_SMALL( zero_target->first,1.0e-15 );
  BOOST_CHECK_SMALL( zero_target->second,1.0e-15 );

  const std::vector<PL::ActivityBoxMassGroup> malformed{
      {1.0,1.0,2.0,1.0}, {1.0,0.0,0.0,1.0}
  };
  BOOST_CHECK( !PL::activity_box_fraction_bounds(malformed) );
}


BOOST_AUTO_TEST_CASE( trace_profile_controls_preserve_interior_and_scale_endpoint_surrogates )
{
  constexpr double upper = 9.0e-10;
  constexpr double interior = 4.0e-10;
  BOOST_CHECK_EQUAL( PL::conditional_solve_control(interior,interior,0.0,upper),interior );
  BOOST_CHECK_EQUAL( PL::conditional_solve_control(upper,interior,0.0,upper),upper );

  const double zero_surrogate = PL::conditional_solve_control(0.0,interior,0.0,upper);
  BOOST_CHECK_GT( zero_surrogate,0.0 );
  BOOST_CHECK_LT( zero_surrogate,1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( zero_surrogate,2.0e-6*interior,1.0e-12 );

  const double tolerance = PL::conditional_equality_tolerance(interior,interior,0.0,upper);
  BOOST_CHECK_GT( tolerance,0.0 );
  BOOST_CHECK_LT( tolerance,1.0e-12 );
  BOOST_CHECK_LT( tolerance,0.001*upper );

  // The upper input boundary is an ordinary interior probability and must be evaluated exactly;
  // the old absolute clamp incorrectly replaced it with 2e-6.
  bool saw_exact_upper = false;
  const auto scan = PL::scan(interior,interior,0.0,upper,{{1.0,sm_threshold95}},32,
    [&saw_exact_upper]( const double requested, const PL::Direction, const std::size_t ) {
      constexpr double bound = 9.0e-10;
      constexpr double nominal = 4.0e-10;
      const double evaluated = PL::conditional_solve_control(requested,nominal,0.0,bound);
      saw_exact_upper = saw_exact_upper || (requested == bound && evaluated == bound);
      return success(evaluated,0.1*std::pow((evaluated-4.0e-10)/bound,2.0));
    });
  BOOST_CHECK( saw_exact_upper );
  BOOST_CHECK( scan.status == PL::ScanStatus::NonIdentifiable );
  BOOST_CHECK_CLOSE_FRACTION( scan.intervals[1].upper.reported_fraction,upper,1.0e-12 );

  // A trace baseline in the full physical domain must keep its lower endpoint surrogate on the
  // lower side of the optimum and use a trace-scale AL tolerance, not the full unit span.
  constexpr double trace_baseline = 1.0e-9;
  const double full_range_zero
      = PL::conditional_solve_control(0.0,trace_baseline,0.0,1.0);
  BOOST_CHECK_GT( full_range_zero,0.0 );
  BOOST_CHECK_LT( full_range_zero,trace_baseline );
  BOOST_CHECK_CLOSE_FRACTION( full_range_zero,2.0e-6*trace_baseline,1.0e-12 );
  const double full_range_lower_tolerance
      = PL::conditional_equality_tolerance(full_range_zero,trace_baseline,0.0,1.0);
  BOOST_CHECK_GT( full_range_lower_tolerance,0.0 );
  BOOST_CHECK_LT( full_range_lower_tolerance,1.0e-12 );

  // A lower-side trace request is accepted at honest per-mille coordinate precision.  Production
  // stores the actual independently fitted coordinate, so this avoids a false profile failure
  // without labeling the requested point as exact.
  constexpr double pu241_baseline = 9.2053390775993291e-6;
  constexpr double pu241_request = 8.63000538525e-6;
  const double pu241_tolerance
      = PL::conditional_equality_tolerance(pu241_request,pu241_baseline,0.0,1.0);
  BOOST_CHECK_GT( pu241_tolerance,8.33474304e-9 );
  BOOST_CHECK_LE( pu241_tolerance,1.0e-8 );

  // Exact sigma-block constraints need at least 1e-6 remainder for an unconstrained sibling even
  // when the baseline is already close to one.  Reported-coordinate AL callers pass no minimum.
  const double exact_upper
      = PL::conditional_solve_control(1.0,0.9,0.0,1.0,2.0e-6,2.0e-6);
  BOOST_CHECK_LE( exact_upper,1.0-2.0e-6 );
  BOOST_CHECK_GT( exact_upper,0.9 );

  const auto narrow_crossing = PL::scan(4.5e-10,4.5e-10,0.0,9.0e-10,
      {{1.0,sm_threshold95}},32,
      []( const double x, const PL::Direction, const std::size_t ) {
        return success(x,std::pow((x-4.5e-10)/1.0e-10,2.0));
      });
  BOOST_REQUIRE( narrow_crossing.status == PL::ScanStatus::Complete );
  BOOST_CHECK_CLOSE_FRACTION(
      narrow_crossing.intervals[0].lower.reported_fraction,3.5e-10,0.01 );
  BOOST_CHECK_CLOSE_FRACTION(
      narrow_crossing.intervals[0].upper.reported_fraction,5.5e-10,0.01 );
}


BOOST_AUTO_TEST_CASE( cancellation_failure_and_better_baseline_are_distinct )
{
  const auto canceled = PL::scan(0.5,0.5,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double, const PL::Direction, const std::size_t ) {
      PL::Evaluation result;
      result.status = PL::EvaluationStatus::Canceled;
      return result;
    });
  BOOST_CHECK( canceled.status == PL::ScanStatus::Canceled );

  const auto failed = PL::scan(0.5,0.5,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double, const PL::Direction, const std::size_t ) {
      PL::Evaluation result;
      result.status = PL::EvaluationStatus::Failed;
      result.diagnostic = "synthetic conditional failure";
      return result;
    });
  BOOST_CHECK( failed.status == PL::ScanStatus::Failed );
  BOOST_CHECK_EQUAL( failed.diagnostic,"synthetic conditional failure" );

  const auto better = PL::scan(0.5,0.5,0.0,1.0,{{1.0,sm_threshold95}},32,
    []( const double, const PL::Direction, const std::size_t ) {
      PL::Evaluation result;
      result.status = PL::EvaluationStatus::BetterBaseline;
      result.diagnostic = "synthetic better point";
      return result;
    });
  BOOST_CHECK( better.status == PL::ScanStatus::BetterBaseline );
  BOOST_CHECK_EQUAL( better.diagnostic,"synthetic better point" );
}
