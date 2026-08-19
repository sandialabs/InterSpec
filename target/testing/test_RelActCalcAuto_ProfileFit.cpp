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

#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <cstddef>
#include <limits>
#include <algorithm>

#define BOOST_TEST_MODULE RelActCalcAuto_ProfileFit_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/RelActCalcAuto_ProfileFit.h"

namespace PL = RelActCalcAutoImp::ProfileLikelihood;

namespace
{
constexpr double sm_threshold_68 = 1.0;
constexpr double sm_threshold_95 = 3.841458820694124;
const std::array<double,2> sm_thresholds{{sm_threshold_68, sm_threshold_95}};

PL::FitPoint point( const double d, const double delta_chi2, const std::size_t index = 0 )
{
  PL::FitPoint answer;
  answer.d = d;
  answer.delta_chi2 = delta_chi2;
  answer.sample_index = index;
  return answer;
}

PL::Evaluation success( const double reported, const double delta_chi2,
                        const bool reached_bound = false )
{
  PL::Evaluation answer;
  answer.status = PL::EvaluationStatus::Success;
  answer.reported_fraction = reported;
  answer.delta_chi2 = delta_chi2;
  answer.reached_feasible_bound = reached_bound;
  return answer;
}

/** An analytic evaluator: the reported coordinate equals the control coordinate, and the objective
 is an exact even polynomial about the baseline.  This is the same "no solver dependency" property
 that makes the scanner unit-testable.
 */
struct AnalyticProfile
{
  double baseline = 0.5;
  double c2 = 100.0;
  double c4 = 0.0;
  std::size_t evaluations = 0;
  std::vector<double> requested;

  PL::Evaluator evaluator()
  {
    return [this]( const double control, PL::Direction, std::size_t ) {
      ++evaluations;
      requested.push_back( control );
      const double d = control - baseline;
      return success( control, c2*d*d + c4*d*d*d*d );
    };
  }
};
}//namespace


/** The vertex is known a priori, so one off-optimum point already determines the curvature.  This
 is the whole reason the model is one-parameter rather than three.
 */
BOOST_AUTO_TEST_CASE( a_single_point_and_the_vertex_determine_the_curvature )
{
  const std::vector<PL::FitPoint> points{ point(0.2,4.0) };
  const PL::VertexModel model = PL::fit_vertex_model( points, 1 );
  BOOST_REQUIRE( model.usable() );
  BOOST_CHECK_CLOSE( model.c2, 100.0, 1.0e-9 );

  // ... and the crossing follows in closed form.
  const std::optional<double> root = PL::solve_crossing( model, sm_threshold_95 );
  BOOST_REQUIRE( root.has_value() );
  BOOST_CHECK_CLOSE( *root, std::sqrt(sm_threshold_95/100.0), 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( an_exact_quadratic_is_recovered_to_machine_precision )
{
  const double c2 = 37.25;
  std::vector<PL::FitPoint> points;
  for( const double d : {0.10,0.20,0.30,0.40,0.50} )
    points.push_back( point(d,c2*d*d) );

  const PL::VertexModel model = PL::fit_best_model( points );
  BOOST_REQUIRE( model.usable() );
  const std::optional<double> root = PL::solve_crossing( model, sm_threshold_95 );
  BOOST_REQUIRE( root.has_value() );
  BOOST_CHECK_SMALL( *root - std::sqrt(sm_threshold_95/c2), 1.0e-12 );
}


BOOST_AUTO_TEST_CASE( a_quartic_profile_is_recovered_by_the_fit )
{
  const double c2 = 20.0, c4 = 900.0;
  std::vector<PL::FitPoint> points;
  for( const double d : {0.05,0.12,0.20,0.28,0.36} )
    points.push_back( point(d,c2*d*d + c4*d*d*d*d) );

  const PL::VertexModel model = PL::fit_best_model( points );
  BOOST_REQUIRE( model.usable() );

  PL::VertexModel accepted;
  const std::optional<double> root = PL::guarded_crossing( model, sm_threshold_95, accepted );
  BOOST_REQUIRE( root.has_value() );
  // Closed-form biquadratic solution of c4*d^4 + c2*d^2 == T.
  const double expected_d2 = (-c2 + std::sqrt(c2*c2 + 4.0*c4*sm_threshold_95))/(2.0*c4);
  BOOST_CHECK_SMALL( *root - std::sqrt(expected_d2), 1.0e-9 );
}


/** A concave quartic can turn over below the threshold and never reach it.  Taking the square root
 of the negative discriminant would produce a NaN endpoint that every downstream check would pass.
 */
BOOST_AUTO_TEST_CASE( a_model_that_never_reaches_the_threshold_falls_back_rather_than_returning_nan )
{
  PL::VertexModel model;
  model.order = 2;
  model.c2 = 1.0;
  model.c4 = -1.0;   // Peaks at dchi2 = 0.25, far below the 95% threshold.
  BOOST_REQUIRE( model.usable() );
  BOOST_REQUIRE_LT( model.c2*model.c2 + 4.0*model.c4*sm_threshold_95, 0.0 );

  const std::optional<double> root = PL::solve_crossing( model, sm_threshold_95 );
  BOOST_CHECK( !root.has_value() );

  // The guarded solve steps the order down instead of failing outright.
  PL::VertexModel accepted;
  const std::optional<double> guarded = PL::guarded_crossing( model, sm_threshold_95, accepted );
  BOOST_REQUIRE( guarded.has_value() );
  BOOST_CHECK_EQUAL( accepted.order, 1U );
  BOOST_CHECK_SMALL( *guarded - std::sqrt(sm_threshold_95/model.c2), 1.0e-12 );
}


/** The convexity guard must FAIL on a concave fit.  Its predecessor - "strictly increasing on
 [0,d_root]" - is vacuous here: for c4 < 0 the selected branch always lies before the turn-over, so
 an increasing model is guaranteed and the test passes for every concave fit, which is precisely the
 failure it was meant to catch.
 */
BOOST_AUTO_TEST_CASE( the_convexity_guard_rejects_a_concave_fit_a_monotonicity_guard_would_accept )
{
  PL::VertexModel model;
  model.order = 2;
  model.c2 = 40.0;
  model.c4 = -60.0;

  const std::optional<double> root = PL::solve_crossing( model, sm_threshold_95 );
  BOOST_REQUIRE( root.has_value() );

  // The model IS increasing at the root - the guard a monotonicity test would have applied.
  BOOST_CHECK_GT( model.slope(*root), 0.0 );
  // But it is concave there, so the convexity guard rejects it.
  BOOST_CHECK_LT( model.curvature(*root), 0.0 );
  BOOST_CHECK( !PL::model_is_convex_and_increasing(model,*root) );

  // The guarded solve therefore steps down to the quadratic, which is convex by construction.
  PL::VertexModel accepted;
  const std::optional<double> guarded = PL::guarded_crossing( model, sm_threshold_95, accepted );
  BOOST_REQUIRE( guarded.has_value() );
  BOOST_CHECK_EQUAL( accepted.order, 1U );
}


/** The recorded Pu-240 conditional-fit sequence, which is what killed Brent: delta chi2 of
 9.783 / 11.778 / 11.043 / 23.129 / 11.034 with the 95% threshold at 16.95 sitting squarely inside
 the scatter band.

 Note the outlier is the OUTERMOST sample.  A rule of "never drop the outermost point, because that
 is where the real curvature lives" would protect the outlier and then discard a good point on the
 next pass - which is why isolation, not position, is what gets tested.
 */
BOOST_AUTO_TEST_CASE( one_sided_rejection_drops_the_isolated_outlier_even_when_it_is_outermost )
{
  const double threshold_95 = 16.95;
  const double tau = PL::rejection_tolerance( threshold_95 );
  BOOST_CHECK_CLOSE( tau, 3.39, 1.0e-9 );

  // Four points reach the rejection stage: the near-duplicate pair has already collapsed.  All sit
  // at d ~ 0.044 with a 0.4% spread, so they identify exactly one parameter.
  const std::vector<PL::FitPoint> points{
    point(0.04428, 9.783, 0),
    point(0.04438,11.778, 1),
    point(0.04443,11.034, 2),
    point(0.04445,23.129, 3)
  };

  const PL::VertexModel model = PL::fit_best_model( points );
  BOOST_REQUIRE( model.usable() );
  BOOST_CHECK_EQUAL( model.order, 1U );   // Clustered points cannot identify more.

  const PL::RejectionResult rejection = PL::reject_one_outlier( points, tau, threshold_95 );
  BOOST_REQUIRE( rejection.dropped.has_value() );
  BOOST_CHECK_CLOSE( rejection.dropped->delta_chi2, 23.129, 1.0e-9 );
  BOOST_CHECK_EQUAL( rejection.kept.size(), 3U );
  for( const PL::FitPoint &kept : rejection.kept )
    BOOST_CHECK_LT( kept.delta_chi2, 12.0 );

  // 23.129 was the only sample above the threshold, so the drop is not final - it must trigger an
  // outward extension and a retest rather than silently destroying the bracket.
  BOOST_CHECK( rejection.drop_removed_sole_anchor );
}


/** The guard against rejection eating genuine curvature: real curvature lifts a point AND its
 neighbours, so two consistent high points are never dropped.
 */
BOOST_AUTO_TEST_CASE( a_pair_of_high_outer_points_is_curvature_and_is_not_dropped )
{
  const double tau = PL::rejection_tolerance( sm_threshold_95 );
  // A profile genuinely steeper than the model can follow.  Every outboard point is high, and they
  // are high TOGETHER - which is what curvature looks like and what a failed solve never does.
  const std::vector<PL::FitPoint> points{
    point(0.10, 1.00,0), point(0.20, 4.00,1), point(0.30,10.00,2),
    point(0.40,30.00,3), point(0.50,80.00,4)
  };

  const PL::RejectionResult rejection = PL::reject_one_outlier( points, tau, sm_threshold_95 );
  BOOST_CHECK( !rejection.dropped.has_value() );
  BOOST_CHECK_EQUAL( rejection.kept.size(), points.size() );
}


BOOST_AUTO_TEST_CASE( a_low_point_is_never_dropped )
{
  const double tau = PL::rejection_tolerance( sm_threshold_95 );
  const std::vector<PL::FitPoint> points{
    point(0.10,1.00,0), point(0.20,4.00,1), point(0.30,9.00,2),
    point(0.40,0.10,3),                       // Absurdly low - impossible for a minimization.
    point(0.50,25.00,4)
  };

  const PL::RejectionResult rejection = PL::reject_one_outlier( points, tau, sm_threshold_95 );
  if( rejection.dropped )
    BOOST_CHECK_GT( rejection.dropped->delta_chi2, 1.0 );
}


BOOST_AUTO_TEST_CASE( at_most_one_point_is_dropped_per_side )
{
  const double tau = PL::rejection_tolerance( sm_threshold_95 );
  const std::vector<PL::FitPoint> points{
    point(0.10,1.00,0), point(0.20,4.00,1), point(0.30,90.0,2),
    point(0.40,16.0,3), point(0.50,25.0,4)
  };

  const PL::RejectionResult rejection = PL::reject_one_outlier( points, tau, sm_threshold_95 );
  // Whatever it decides, it can never remove more than one: with five points and three parameters a
  // second drop would leave the model interpolating and every residual identically zero.
  BOOST_CHECK_GE( rejection.kept.size(), points.size() - 1 );
}


/** Section 2.5: two requests landing on the same reported coordinate keep the LOWER objective. */
BOOST_AUTO_TEST_CASE( near_duplicate_reported_values_collapse_and_keep_the_lower_objective )
{
  std::vector<PL::Sample> samples;
  for( const std::pair<double,double> &entry : std::vector<std::pair<double,double>>{
        {0.18298, 9.783},{0.18308,11.778},{0.18312,11.043},{0.18315,23.129},{0.18313,11.034} } )
  {
    PL::Sample sample;
    sample.control_fraction = entry.first;
    sample.reported_fraction = entry.first;
    sample.delta_chi2 = entry.second;
    samples.push_back( sample );
  }

  // Baseline of 0.1387 gives a near-duplicate threshold of 1.387e-5; the 0.18312/0.18313 pair is
  // 1.0e-5 apart and therefore collapses, keeping the LOWER objective by the one-sided rule.  Four
  // samples - not five - reach the rejection stage.
  const PL::SideHygiene hygiene = PL::side_hygiene( samples, 0.1387, 0.1387 );
  BOOST_CHECK( !hygiene.folded );   // Ordered by control, the reported values are monotone.
  BOOST_REQUIRE_EQUAL( hygiene.points.size(), 4U );
  for( std::size_t i = 1; i < hygiene.points.size(); ++i )
    BOOST_CHECK_GT( hygiene.points[i].d, hygiene.points[i-1].d );

  const bool kept_the_lower_of_the_pair = std::any_of( hygiene.points.begin(),
      hygiene.points.end(), []( const PL::FitPoint &p ){
        return std::fabs(p.delta_chi2 - 11.034) < 1.0e-9;
      } );
  const bool dropped_the_higher_of_the_pair = std::none_of( hygiene.points.begin(),
      hygiene.points.end(), []( const PL::FitPoint &p ){
        return std::fabs(p.delta_chi2 - 11.043) < 1.0e-9;
      } );
  BOOST_CHECK( kept_the_lower_of_the_pair );
  BOOST_CHECK( dropped_the_higher_of_the_pair );
}


/** The common case: probe, rescale, place the ladder past the crossing, and never need to extend. */
BOOST_AUTO_TEST_CASE( probe_and_rescale_brackets_the_crossing_without_an_extension )
{
  AnalyticProfile profile;
  profile.baseline = 0.5;
  profile.c2 = 400.0;   // 95% crossing at d = sqrt(3.8415/400) = 0.0980

  // Supply the local one-sigma, as production does: d(1 sigma) = sqrt(1/c2) = 0.05.
  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 profile.evaluator(),0.05 );

  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::Complete, result.diagnostic );
  const double expected = std::sqrt( sm_threshold_95/400.0 );
  BOOST_CHECK_SMALL( result.intervals[1].upper.reported_fraction - (0.5 + expected), 1.0e-6 );
  BOOST_CHECK_SMALL( result.intervals[1].lower.reported_fraction - (0.5 - expected), 1.0e-6 );

  const double expected_68 = std::sqrt( sm_threshold_68/400.0 );
  BOOST_CHECK_SMALL( result.intervals[0].upper.reported_fraction - (0.5 + expected_68), 1.0e-6 );
  BOOST_CHECK_SMALL( result.intervals[0].lower.reported_fraction - (0.5 - expected_68), 1.0e-6 );

  BOOST_CHECK( result.intervals[1].lower.likelihood_crossing );
  BOOST_CHECK( result.intervals[1].upper.likelihood_crossing );

  // Five ladder points per side, no extension and no re-probe needed.
  BOOST_CHECK_EQUAL( profile.evaluations, 2*PL::sm_profile_points_per_side );
}


/** With no local sigma available - the common case for an automatically-profiled weak quantity - the
 probe falls back to a fraction of the span and can land far inside the crossing.  The reserved
 re-probe is what absorbs that, and it must not come out of any other side's allowance.
 */
BOOST_AUTO_TEST_CASE( without_a_local_sigma_the_reserved_reprobe_absorbs_a_probe_that_lands_short )
{
  AnalyticProfile profile;
  profile.baseline = 0.5;
  profile.c2 = 400.0;

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 profile.evaluator(),0.0 );
  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::Complete, result.diagnostic );

  const double expected = std::sqrt( sm_threshold_95/400.0 );
  BOOST_CHECK_SMALL( result.intervals[1].upper.reported_fraction - (0.5 + expected), 1.0e-6 );

  // Ladder plus the reserved re-probe, on both sides, and nothing beyond it.
  BOOST_CHECK_EQUAL( profile.evaluations,
                     2*(PL::sm_profile_points_per_side + PL::sm_profile_reprobe_per_side) );
}


/** Interpolation only: the reported crossing must be bracketed by an evaluated sample, so a probe
 that lands far inside the crossing must provoke outward travel rather than an extrapolated answer.
 */
BOOST_AUTO_TEST_CASE( a_crossing_is_never_reported_outside_the_evaluated_span )
{
  AnalyticProfile profile;
  profile.baseline = 0.5;
  profile.c2 = 400.0;

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 profile.evaluator(),0.05 );
  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::Complete, result.diagnostic );

  for( std::size_t level = 0; level < 2; ++level )
  {
    for( const PL::Endpoint &endpoint : { result.intervals[level].lower,
                                          result.intervals[level].upper } )
    {
      if( !endpoint.likelihood_crossing )
        continue;
      // Some evaluated sample must lie at or beyond the reported endpoint, on that side.
      const bool bracketed = std::any_of( result.samples.begin(), result.samples.end(),
        [&]( const PL::Sample &sample ){
          return std::fabs(sample.reported_fraction - 0.5)
                 >= std::fabs(endpoint.reported_fraction - 0.5) - 1.0e-12;
        } );
      BOOST_CHECK_MESSAGE( bracketed, "endpoint " << endpoint.reported_fraction
                           << " was not bracketed by an evaluated sample" );
    }
  }
}


/** Reaching the feasible bound below the threshold is a legitimate, correct answer - the likelihood
 genuinely never reaches the threshold inside the feasible set - and must never be conflated with
 "we did not push far enough".
 */
BOOST_AUTO_TEST_CASE( reaching_the_feasible_bound_below_the_threshold_is_boundary_limited )
{
  AnalyticProfile profile;
  profile.baseline = 0.5;
  profile.c2 = 1.0;   // 95% crossing would be at d = 1.96, far outside the [0.45,0.55] window.

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.45,0.55,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 profile.evaluator(),0.0 );

  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::NonIdentifiable
                           || result.status == PL::ScanStatus::BoundaryLimited,
                         result.diagnostic );
  BOOST_CHECK( !result.intervals[1].lower.likelihood_crossing );
  BOOST_CHECK( !result.intervals[1].upper.likelihood_crossing );
  BOOST_CHECK_SMALL( result.intervals[1].lower.reported_fraction - 0.45, 1.0e-9 );
  BOOST_CHECK_SMALL( result.intervals[1].upper.reported_fraction - 0.55, 1.0e-9 );
}


/** The whole feasible range inside the 95% threshold on both sides is `NonIdentifiable`, exactly as
 before - and it is the state the automatic-weak trigger is defined by.
 */
BOOST_AUTO_TEST_CASE( an_entirely_interior_feasible_range_is_non_identifiable )
{
  AnalyticProfile profile;
  profile.baseline = 0.5;
  profile.c2 = 0.01;

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.40,0.60,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 profile.evaluator(),0.0 );
  BOOST_CHECK_MESSAGE( result.status == PL::ScanStatus::NonIdentifiable, result.diagnostic );
}


/** The point pool is a hard cap and is never exceeded, whatever the profile does. */
BOOST_AUTO_TEST_CASE( the_point_budget_is_fixed_and_never_exceeded )
{
  for( const double curvature : {0.001,0.05,1.0,25.0,400.0,1.0e4} )
  {
    AnalyticProfile profile;
    profile.baseline = 0.5;
    profile.c2 = curvature;

    const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                   PL::sm_profile_max_points_per_quantity,
                                                   profile.evaluator(),0.0 );
    BOOST_CHECK_MESSAGE( profile.evaluations <= PL::sm_profile_max_points_per_quantity,
                         "curvature " << curvature << " used " << profile.evaluations
                                      << " points" );
    BOOST_CHECK_MESSAGE( result.status != PL::ScanStatus::Failed,
        "curvature " << curvature << " status=" << static_cast<int>(result.status)
                     << " points=" << profile.evaluations << ": " << result.diagnostic );
  }
}


/** Exhausting the budget with evaluated samples straddling the threshold yields an evaluated
 endpoint, not a failure.  This is the regression the shipped behaviour must keep - a rare hard
 quantity should spend more iterations, not report nothing.

 The mechanism is the anchor clamp: an evaluated sample at or above the threshold, together with the
 vertex below it, brackets the crossing by construction, so that sample is reported.  It sits at or
 beyond the true crossing, so the interval errs wide - the safe direction.
 */
BOOST_AUTO_TEST_CASE( budget_exhaustion_with_a_straddling_pair_yields_an_evaluated_endpoint )
{
  // A profile the polynomial model cannot follow: flat, then a wall.  Every model root will sit
  // outside the evaluated span, so only the straddling-pair rule can produce an endpoint.
  std::size_t evaluations = 0;
  const PL::Evaluator evaluator = [&evaluations]( const double control, PL::Direction,
                                                  std::size_t ) {
    ++evaluations;
    const double d = std::fabs( control - 0.5 );
    const double delta = (d < 0.30) ? (0.02*d) : (500.0*(d - 0.30) + 0.006);
    return success( control, delta );
  };

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 evaluator,0.0 );

  BOOST_CHECK_MESSAGE( result.status != PL::ScanStatus::Failed, result.diagnostic );
  BOOST_CHECK( evaluations <= PL::sm_profile_max_points_per_quantity );
}


/** A folded control-to-reported map is `Failed` with a diagnostic naming the fold, never a silent
 crossing: a sample beyond a fold is not outboard in the reported coordinate and cannot anchor
 anything.
 */
BOOST_AUTO_TEST_CASE( a_folded_mapping_is_failed_and_not_a_silent_crossing )
{
  // The reported coordinate turns around once the control coordinate passes 0.6.
  const PL::Evaluator evaluator = []( const double control, PL::Direction, std::size_t ) {
    const double reported = (control <= 0.6) ? control : (1.2 - control);
    return success( reported, 4.0*(reported - 0.5)*(reported - 0.5) );
  };

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 evaluator,0.0 );
  BOOST_CHECK( result.status != PL::ScanStatus::Complete );
}


/** A non-positive curvature at the optimum means the reported solution is not a minimum along this
 direction - a failure, and a strong hint the baseline is wrong.
 */
BOOST_AUTO_TEST_CASE( non_positive_curvature_at_the_optimum_is_failed )
{
  std::vector<PL::FitPoint> points;
  for( const double d : {0.10,0.20,0.30,0.40,0.50} )
    points.push_back( point(d,-5.0*d*d) );

  const PL::VertexModel model = PL::fit_best_model( points );
  BOOST_CHECK( !model.usable() );
  BOOST_CHECK( !PL::solve_crossing(model,sm_threshold_95).has_value() );
}


/** A probe that already exceeds the 95% threshold has bracketed the crossing by itself, which is
 the cheapest possible outcome and must not provoke any extension.
 */
BOOST_AUTO_TEST_CASE( a_probe_past_the_threshold_brackets_immediately )
{
  AnalyticProfile profile;
  profile.baseline = 0.5;
  profile.c2 = 1.0e6;   // Even a 0.02-of-span probe lands far above the threshold.

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 profile.evaluator(),0.0 );
  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::Complete, result.diagnostic );
  BOOST_CHECK_EQUAL( profile.evaluations, 2*PL::sm_profile_points_per_side );

  const double expected = std::sqrt( sm_threshold_95/1.0e6 );
  BOOST_CHECK_SMALL( result.intervals[1].upper.reported_fraction - (0.5 + expected), 1.0e-9 );
}


/** A probe that barely moves the objective gives a rescale factor made entirely of noise; the
 reserved re-probe exists for exactly that, and it must be spent.
 */
BOOST_AUTO_TEST_CASE( a_probe_that_barely_moves_the_objective_triggers_the_reserved_reprobe )
{
  AnalyticProfile profile;
  profile.baseline = 0.5;
  profile.c2 = 0.5;   // A 0.02-of-span probe gives dchi2 = 2e-4, far below 0.05*T68.

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 profile.evaluator(),0.0 );
  BOOST_CHECK_MESSAGE( result.status != PL::ScanStatus::Failed, result.diagnostic );
  BOOST_CHECK_GT( profile.evaluations, 2*PL::sm_profile_points_per_side );
  BOOST_CHECK_LE( profile.evaluations, PL::sm_profile_max_points_per_quantity );
}


/** The budget arithmetic closes with room to spare, which is what makes it order-independent: with
 both sides' ladders and re-probes reserved, the worst-case extension draw still fits in the
 surplus, so neither side can starve the other in any order.
 */
BOOST_AUTO_TEST_CASE( the_reserved_budget_leaves_both_sides_their_full_allowance )
{
  const std::size_t reserved = 2*(PL::sm_profile_points_per_side + PL::sm_profile_reprobe_per_side);
  const std::size_t worst_case_extensions = 2*PL::sm_profile_max_extensions_per_side;
  BOOST_CHECK_LE( reserved + worst_case_extensions, PL::sm_profile_max_points_per_quantity );
  BOOST_CHECK_EQUAL( PL::sm_profile_max_points_per_quantity - reserved, 8U );
}


/** Cancellation and a better baseline are distinct outcomes, and neither is a profile failure. */
BOOST_AUTO_TEST_CASE( cancellation_and_better_baseline_are_distinct )
{
  const PL::Evaluator canceled = []( double, PL::Direction, std::size_t ) {
    PL::Evaluation answer;
    answer.status = PL::EvaluationStatus::Canceled;
    return answer;
  };
  const PL::ScanResult cancel_result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                        PL::sm_profile_max_points_per_quantity,
                                                        canceled,0.0 );
  BOOST_CHECK( cancel_result.status == PL::ScanStatus::Canceled );

  const PL::Evaluator better = []( double, PL::Direction, std::size_t ) {
    PL::Evaluation answer;
    answer.status = PL::EvaluationStatus::BetterBaseline;
    return answer;
  };
  const PL::ScanResult better_result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                        PL::sm_profile_max_points_per_quantity,
                                                        better,0.0 );
  BOOST_CHECK( better_result.status == PL::ScanStatus::BetterBaseline );
}


/** An asymmetric profile is fitted independently per side: forcing symmetry is a real modelling
 error near a bound, and the profile is generically asymmetric.
 */
BOOST_AUTO_TEST_CASE( an_asymmetric_profile_is_fitted_independently_per_side )
{
  const double lower_c2 = 900.0, upper_c2 = 100.0;
  const PL::Evaluator evaluator = [=]( const double control, PL::Direction, std::size_t ) {
    const double d = control - 0.5;
    const double c2 = (d < 0.0) ? lower_c2 : upper_c2;
    return success( control, c2*d*d );
  };

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 evaluator,0.0 );
  BOOST_REQUIRE_MESSAGE( result.status == PL::ScanStatus::Complete, result.diagnostic );
  BOOST_CHECK_SMALL( result.intervals[1].lower.reported_fraction
                       - (0.5 - std::sqrt(sm_threshold_95/lower_c2)), 1.0e-6 );
  BOOST_CHECK_SMALL( result.intervals[1].upper.reported_fraction
                       - (0.5 + std::sqrt(sm_threshold_95/upper_c2)), 1.0e-6 );
}


/** Regression: the near-duplicate tolerance must scale with the range actually explored.

 A fixed `1e-4*|q0|` is right when the scan spans much more than that, and it is what collapses the
 recorded Pu-240 pair.  But where the reported quantity responds only weakly to the pinned
 coordinate, the whole profile can live INSIDE that tolerance - and then every sample collapses into
 one, and a side with a perfectly clear crossing reports that it never bracketed anything.  Found on
 the ratio-constrained U-235 fixture, whose profile spans ~3e-5 about a baseline of 0.94.
 */
BOOST_AUTO_TEST_CASE( near_duplicate_collapse_never_swallows_the_whole_sampled_span )
{
  const double baseline = 0.943338;
  std::vector<PL::Sample> samples;
  // Five distinct points spanning 3.3e-5 - well inside the fixed tolerance of 1e-4*0.943 = 9.4e-5.
  for( const std::pair<double,double> &entry : std::vector<std::pair<double,double>>{
        {1.26e-5,2.88},{1.46e-5,4.04},{1.95e-5,8.18},{2.64e-5,18.49},{3.33e-5,37.80} } )
  {
    PL::Sample sample;
    sample.control_fraction = entry.first;      //control coordinate, monotone with the reported one
    sample.reported_fraction = baseline - entry.first;
    sample.delta_chi2 = entry.second;
    samples.push_back( sample );
  }

  const PL::SideHygiene hygiene = PL::side_hygiene( samples, baseline, 0.0 );
  BOOST_CHECK( !hygiene.folded );
  BOOST_CHECK_MESSAGE( hygiene.points.size() >= 4U,
                       "the sampled span was collapsed to " << hygiene.points.size()
                       << " point(s) by a tolerance wider than the span itself" );
}


/** Regression: when the true profile is steeper than the fitted polynomial, the model's root lands
 OUTSIDE the innermost measured sample that reached the threshold.  The crossing is bracketed all the
 same - the vertex is below it by construction, that sample is at or above it - so the measurement is
 reported rather than the profile being failed.  It errs wide, which is the safe direction.
 */
BOOST_AUTO_TEST_CASE( a_root_overshooting_its_measured_anchor_reports_the_anchor )
{
  // Every evaluated point is already above the 68% threshold, so there is no evaluated "inside"
  // sample at that level and only the vertex can serve as the inboard anchor.
  const PL::Evaluator evaluator = []( const double control, PL::Direction, std::size_t ) {
    const double d = std::fabs( control - 0.5 );
    // Steeper than quartic: a polynomial fitted to the outer points underpredicts near the inside.
    return success( control, 4.0e3*d*d + 8.0e9*d*d*d*d*d*d );
  };

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 evaluator,0.0 );
  BOOST_CHECK_MESSAGE( result.status != PL::ScanStatus::Failed, result.diagnostic );

  // Whatever endpoint is reported must be bracketed by an evaluated sample - never past one.
  for( const PL::Endpoint &endpoint : { result.intervals[0].lower, result.intervals[0].upper,
                                        result.intervals[1].lower, result.intervals[1].upper } )
  {
    if( !endpoint.likelihood_crossing )
      continue;
    const bool bracketed = std::any_of( result.samples.begin(), result.samples.end(),
      [&]( const PL::Sample &sample ){
        return std::fabs(sample.reported_fraction - 0.5)
               >= std::fabs(endpoint.reported_fraction - 0.5) - 1.0e-12;
      } );
    BOOST_CHECK( bracketed );
  }
}


/** Regression: extensions must stop once a MEASUREMENT has bracketed the threshold.

 Extending cannot repair a model that underpredicts at its anchor, so continuing to travel outward
 burns the point budget and then reports failure for a quantity whose data already answer the
 question - the precise regression the plan forbids ("allow more iterations to not fail rare
 instances", not spend them chasing a model).
 */
BOOST_AUTO_TEST_CASE( extensions_stop_once_a_measurement_brackets_the_threshold )
{
  std::size_t evaluations = 0;
  const PL::Evaluator evaluator = [&evaluations]( const double control, PL::Direction,
                                                  std::size_t ) {
    ++evaluations;
    const double d = std::fabs( control - 0.5 );
    return success( control, 4.0e3*d*d + 8.0e9*d*d*d*d*d*d );
  };

  const PL::ScanResult result = PL::fit_profile( 0.5,0.5,0.0,1.0,sm_thresholds,
                                                 PL::sm_profile_max_points_per_quantity,
                                                 evaluator,0.0 );
  BOOST_CHECK_MESSAGE( result.status != PL::ScanStatus::Failed, result.diagnostic );
  // The ladder alone brackets this profile; no extension should be spent.
  BOOST_CHECK_LE( evaluations,
                  2*(PL::sm_profile_points_per_side + PL::sm_profile_reprobe_per_side) );
}


// -------------------------------------------------------------------------------------------------
//  Conditional-solve convergence: "stuck" vs "at a constrained optimum"
// -------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( a_seed_convergence_is_accepted_without_probing )
{
  // Ceres counts the initial evaluation as a successful step, so a solve that started AT its
  // conditional optimum reports one success and - the load-bearing part - zero rejections.
  BOOST_CHECK( PL::conditional_convergence_verdict(true,1,0)
               == PL::ConditionalConvergence::SeedConverged );
  BOOST_CHECK( PL::conditional_convergence_verdict(true,0,0)
               == PL::ConditionalConvergence::SeedConverged );
}


BOOST_AUTO_TEST_CASE( accepting_nothing_while_rejecting_steps_is_ambiguous_not_a_failure )
{
  // The observed JRC Pu70 signature.  It must be neither accepted outright nor discarded outright:
  // a stuck solve and a genuine constrained optimum both look exactly like this.
  BOOST_CHECK( PL::conditional_convergence_verdict(true,1,4)
               == PL::ConditionalConvergence::NeedsProbe );
  BOOST_CHECK( PL::conditional_convergence_verdict(true,1,1)
               == PL::ConditionalConvergence::NeedsProbe );
}


BOOST_AUTO_TEST_CASE( a_solve_that_accepted_real_steps_is_converged )
{
  BOOST_CHECK( PL::conditional_convergence_verdict(true,2,0)
               == PL::ConditionalConvergence::Converged );
  BOOST_CHECK( PL::conditional_convergence_verdict(true,17,9)
               == PL::ConditionalConvergence::Converged );
}


BOOST_AUTO_TEST_CASE( a_solve_ceres_did_not_call_successful_is_rejected_whatever_its_step_counts )
{
  BOOST_CHECK( PL::conditional_convergence_verdict(false,17,0)
               == PL::ConditionalConvergence::Rejected );
  BOOST_CHECK( PL::conditional_convergence_verdict(false,1,4)
               == PL::ConditionalConvergence::Rejected );
  BOOST_CHECK( PL::conditional_convergence_verdict(false,1,0)
               == PL::ConditionalConvergence::Rejected );
}


BOOST_AUTO_TEST_CASE( ceres_no_bound_sentinels_are_not_mistaken_for_constraints )
{
  BOOST_CHECK( !PL::has_real_bound( std::numeric_limits<double>::max() ) );
  BOOST_CHECK( !PL::has_real_bound( std::numeric_limits<double>::lowest() ) );
  BOOST_CHECK( !PL::has_real_bound( std::numeric_limits<double>::infinity() ) );
  BOOST_CHECK( !PL::has_real_bound( std::numeric_limits<double>::quiet_NaN() ) );
  BOOST_CHECK( PL::has_real_bound( 0.0 ) );
  BOOST_CHECK( PL::has_real_bound( -1.0e12 ) );
}


BOOST_AUTO_TEST_CASE( the_tangent_index_map_is_the_increasing_order_complement )
{
  // `ceres::SubsetManifold::Plus` walks the ambient indices in increasing order and consumes one
  // delta entry per non-constant coordinate, so this is the ordering `Problem::Evaluate`'s gradient
  // arrives in.  Getting it wrong would step the wrong parameters - possibly the pin.
  const std::vector<int> constants{{1,3,4}};
  const std::vector<std::size_t> free_indices = PL::tangent_to_ambient_indices( 6, constants );
  const std::vector<std::size_t> expected{{0,2,5}};
  BOOST_CHECK( free_indices == expected );

  BOOST_CHECK( PL::tangent_to_ambient_indices(3,std::vector<int>{}).size() == 3u );
  BOOST_CHECK( PL::tangent_to_ambient_indices(3,std::vector<int>{{0,1,2}}).empty() );
}


BOOST_AUTO_TEST_CASE( a_step_out_of_the_box_at_an_active_bound_is_clamped_away )
{
  // Index 0 sits ON its lower bound with a gradient pushing it further down: it must not move, which
  // is precisely why such a coordinate contributes nothing to a projected-gradient test and why a
  // point like this can be a legitimate constrained optimum rather than a stuck solve.
  // Index 1 is interior and must move.  Index 2 is pinned/constant and must be bit-identical.
  const std::vector<double> x{{0.0, 5.0, 42.0}};
  const std::vector<std::pair<double,double>> box{{ {0.0,10.0}, {0.0,10.0}, {0.0,100.0} }};
  const std::vector<std::size_t> free_indices{{0,1}};
  const std::vector<double> gradient{{ 3.0, 2.0 }};   //descent is -gradient, so index 0 wants to go negative

  std::vector<double> trial;
  const bool moved = PL::projected_descent_point( x,gradient,free_indices,box,1.0,trial );

  BOOST_CHECK( moved );
  BOOST_CHECK_EQUAL( trial[0], 0.0 );    //clamped exactly back onto the bound
  BOOST_CHECK_CLOSE( trial[1], 3.0, 1.0e-9 );
  BOOST_CHECK_EQUAL( trial[2], 42.0 );   //never written
}


BOOST_AUTO_TEST_CASE( a_step_entirely_out_of_the_box_reports_that_nothing_moved )
{
  const std::vector<double> x{{0.0, 10.0}};
  const std::vector<std::pair<double,double>> box{{ {0.0,10.0}, {0.0,10.0} }};
  const std::vector<std::size_t> free_indices{{0,1}};
  const std::vector<double> gradient{{ 1.0, -1.0 }};   //both push outward, off opposite bounds

  std::vector<double> trial;
  BOOST_CHECK( !PL::projected_descent_point(x,gradient,free_indices,box,1.0,trial) );
  BOOST_CHECK( trial == x );
}


BOOST_AUTO_TEST_CASE( the_probe_step_is_scale_free )
{
  // The property the old magnitude test lacked.  The parameters of this problem are activities,
  // energy-calibration terms, FWHM and rel-eff coefficients - sizes differing by many orders of
  // magnitude - so a step measured in raw gradient units is meaningless across them.  Scaling one
  // coordinate's box AND its gradient by the same factor is a pure change of units and must leave
  // the resulting relative displacement identical.
  const std::vector<std::size_t> free_indices{{0,1}};
  const double relative_move = 1.0e-2;

  const std::vector<double> x_a{{5.0, 0.5}};
  const std::vector<std::pair<double,double>> box_a{{ {0.0,10.0}, {0.0,1.0} }};
  const std::vector<double> gradient_a{{ 4.0, 0.25 }};

  const double factor = 1000.0;
  const std::vector<double> x_b{{5.0*factor, 0.5}};
  const std::vector<std::pair<double,double>> box_b{{ {0.0,10.0*factor}, {0.0,1.0} }};
  const std::vector<double> gradient_b{{ 4.0*factor, 0.25 }};

  const double alpha_a = PL::projected_step_scale( x_a,gradient_a,free_indices,box_a,relative_move );
  const double alpha_b = PL::projected_step_scale( x_b,gradient_b,free_indices,box_b,relative_move );
  BOOST_REQUIRE( alpha_a > 0.0 );
  BOOST_REQUIRE( alpha_b > 0.0 );

  std::vector<double> trial_a, trial_b;
  PL::projected_descent_point( x_a,gradient_a,free_indices,box_a,alpha_a,trial_a );
  PL::projected_descent_point( x_b,gradient_b,free_indices,box_b,alpha_b,trial_b );

  // Same fractional move of the rescaled coordinate, and the untouched coordinate is unaffected.
  BOOST_CHECK_CLOSE( (x_a[0]-trial_a[0])/(box_a[0].second-box_a[0].first),
                     (x_b[0]-trial_b[0])/(box_b[0].second-box_b[0].first), 1.0e-9 );
  BOOST_CHECK_CLOSE( trial_a[1], trial_b[1], 1.0e-9 );

  // And the opening rung really does move the leading coordinate by `relative_move` of its span.
  BOOST_CHECK_CLOSE( (x_a[0]-trial_a[0])/(box_a[0].second-box_a[0].first), relative_move, 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( an_unbounded_coordinate_is_scaled_by_its_own_magnitude )
{
  // Ceres' "no bound" sentinel must not be used as a span; a coordinate with no real box is scaled
  // by its own size instead, floored at 1 so a coordinate sitting at zero still has a scale.
  BOOST_CHECK_CLOSE( PL::coordinate_scale(7.0, 0.0, 10.0), 10.0, 1.0e-9 );
  BOOST_CHECK_CLOSE( PL::coordinate_scale(7.0, 1.0, std::numeric_limits<double>::max()), 7.0, 1.0e-9 );
  BOOST_CHECK_CLOSE( PL::coordinate_scale(0.0, std::numeric_limits<double>::lowest(),
                                          std::numeric_limits<double>::max()), 1.0, 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( the_descent_ladder_reaches_far_enough_down_to_be_a_complete_test )
{
  // For a smooth objective a nonzero projected gradient guarantees SOME small enough step descends,
  // so the probe is only a complete test if the ladder spans several decades below its opening rung.
  BOOST_REQUIRE( !PL::sm_descent_probe_ladder.empty() );
  BOOST_CHECK_CLOSE( PL::sm_descent_probe_ladder.front(), 1.0, 1.0e-9 );
  BOOST_CHECK_LE( PL::sm_descent_probe_ladder.back(), 1.0e-3 );
  for( std::size_t i = 1; i < PL::sm_descent_probe_ladder.size(); ++i )
    BOOST_CHECK_LT( PL::sm_descent_probe_ladder[i], PL::sm_descent_probe_ladder[i-1] );
  BOOST_CHECK_GT( PL::sm_descent_probe_relative_move, 0.0 );
  BOOST_CHECK_LT( PL::sm_descent_probe_relative_move, 1.0 );
}


BOOST_AUTO_TEST_CASE( a_coordinate_blocked_at_a_bound_does_not_set_the_step_scale )
{
  // Index 0 has by far the largest gradient but sits ON its lower bound with the gradient pushing
  // further down, so it cannot move.  Letting it set the scale would shrink the step for index 1 by
  // a factor of 100 - and an under-sized probe step under-detects descent, which is the direction
  // that wrongly calls a stuck solve a constrained optimum.
  const std::vector<double> x{{0.0, 5.0}};
  const std::vector<std::pair<double,double>> box{{ {0.0,1.0}, {0.0,10.0} }};
  const std::vector<std::size_t> free_indices{{0,1}};
  const std::vector<double> blocked{{ 100.0, 1.0 }};   //index 0 wants to go below its lower bound
  const std::vector<double> free_grad{{ -100.0, 1.0 }};//same magnitude, but pointing INTO the box

  const double alpha_blocked = PL::projected_step_scale( x,blocked,free_indices,box,1.0e-2 );
  const double alpha_free = PL::projected_step_scale( x,free_grad,free_indices,box,1.0e-2 );

  BOOST_REQUIRE( alpha_blocked > 0.0 );
  BOOST_REQUIRE( alpha_free > 0.0 );
  BOOST_CHECK_GT( alpha_blocked, alpha_free );   //ignoring the blocked coordinate gives a larger step

  std::vector<double> trial;
  PL::projected_descent_point( x,blocked,free_indices,box,alpha_blocked,trial );
  BOOST_CHECK_EQUAL( trial[0], 0.0 );   //still pinned to its bound
  BOOST_CHECK_CLOSE( (x[1]-trial[1])/(box[1].second-box[1].first), 1.0e-2, 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( only_an_outward_gradient_at_a_bound_counts_as_blocked )
{
  BOOST_CHECK( PL::descent_blocked_by_bound(0.0, 1.0, 0.0, 10.0) );    //on lower, wants lower
  BOOST_CHECK( !PL::descent_blocked_by_bound(0.0, -1.0, 0.0, 10.0) );  //on lower, wants higher
  BOOST_CHECK( PL::descent_blocked_by_bound(10.0, -1.0, 0.0, 10.0) );  //on upper, wants higher
  BOOST_CHECK( !PL::descent_blocked_by_bound(10.0, 1.0, 0.0, 10.0) );  //on upper, wants lower
  BOOST_CHECK( !PL::descent_blocked_by_bound(5.0, 1.0, 0.0, 10.0) );   //interior
  // The Ceres "no bound" sentinel must never make a coordinate look blocked.
  BOOST_CHECK( !PL::descent_blocked_by_bound(1.0, 1.0, std::numeric_limits<double>::lowest(),
                                             std::numeric_limits<double>::max()) );
}
