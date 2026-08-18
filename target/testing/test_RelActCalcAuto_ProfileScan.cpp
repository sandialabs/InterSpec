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




