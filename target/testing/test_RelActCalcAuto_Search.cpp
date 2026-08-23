/* InterSpec: focused deterministic RelActAuto search/model-selection tests. */

#include "InterSpec_config.h"

#include <algorithm>
#include <optional>
#include <set>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE RelActCalcAuto_Search_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/RelActCalcAuto_Search_imp.hpp"

using namespace RelActCalcAutoImp;

BOOST_AUTO_TEST_CASE( semantic_candidate_matrix_is_fixed_and_bounded )
{
  const auto variants = search_seed_variants();
  BOOST_REQUIRE_EQUAL( variants.size(), 16U );
  BOOST_CHECK( variants.front() == SearchSeedVariant::Default );

  std::set<std::string> names;
  for( const SearchSeedVariant variant : variants )
    names.insert( search_seed_variant_name(variant) );
  BOOST_CHECK_EQUAL( names.size(), variants.size() );
}


BOOST_AUTO_TEST_CASE( frozen_branching_ratio_hash_covers_exact_ranges_and_membership )
{
  const std::vector<std::pair<double,double>> ranges{{99.0,101.0},{199.0,201.0}};
  std::vector<FrozenPeakRangeGamma> gammas{
    { "2:Pu240", 200.0, 1 },
    { "1:Pu239", 150.0, 1 },
    { "1:Pu239", 100.0, 2 }
  };
  const std::vector<size_t> membership = frozen_peak_range_membership(ranges,gammas);
  BOOST_REQUIRE_EQUAL( membership.size(), 3U );
  BOOST_CHECK_EQUAL( membership[0], 2U );
  BOOST_CHECK_EQUAL( membership[1], 0U );
  BOOST_CHECK_EQUAL( membership[2], 1U );

  const std::uint64_t baseline_hash = frozen_peak_range_policy_hash(ranges,gammas);
  std::reverse(gammas.begin(),gammas.end());
  BOOST_CHECK_EQUAL( frozen_peak_range_policy_hash(ranges,gammas), baseline_hash );

  // The membership is unchanged, but a one-ULP endpoint shift is still a different frozen
  // objective and therefore must produce a different layout fingerprint.
  std::vector<std::pair<double,double>> endpoint_shift = ranges;
  endpoint_shift[0].first = std::nextafter(endpoint_shift[0].first,0.0);
  BOOST_CHECK( frozen_peak_range_membership(endpoint_shift,gammas)
               == frozen_peak_range_membership(ranges,gammas) );
  BOOST_CHECK_NE( frozen_peak_range_policy_hash(endpoint_shift,gammas), baseline_hash );

  // The range count stays fixed while line-to-nuisance assignment changes.
  const std::vector<std::pair<double,double>> membership_shift{{149.0,151.0},{199.0,201.0}};
  BOOST_REQUIRE_EQUAL( membership_shift.size(), ranges.size() );
  BOOST_CHECK_NE( frozen_peak_range_policy_hash(membership_shift,gammas), baseline_hash );
  BOOST_CHECK( frozen_peak_range_membership(membership_shift,gammas)
               != frozen_peak_range_membership(ranges,gammas) );
}


BOOST_AUTO_TEST_CASE( semantic_candidates_copy_the_outer_branching_ratio_policy )
{
  const std::vector<std::pair<double,double>> outer_ranges{{99.0,101.0},{199.0,201.0}};
  const std::vector<std::pair<double,double>> candidate_a_seed{{89.0,91.0},{209.0,211.0}};
  const std::vector<std::pair<double,double>> candidate_b_seed{{109.0,111.0}};
  const std::vector<FrozenPeakRangeGamma> gammas{
    { "1:Pu239", 100.0, 1 }, { "2:Pu240", 200.0, 1 }
  };

  const std::vector<std::pair<double,double>> candidate_a
      = peak_ranges_for_frozen_solve(candidate_a_seed,&outer_ranges);
  const std::vector<std::pair<double,double>> candidate_b
      = peak_ranges_for_frozen_solve(candidate_b_seed,&outer_ranges);
  BOOST_CHECK( candidate_a == outer_ranges );
  BOOST_CHECK( candidate_b == outer_ranges );
  BOOST_CHECK_EQUAL( frozen_peak_range_policy_hash(candidate_a,gammas),
                     frozen_peak_range_policy_hash(candidate_b,gammas) );
  BOOST_CHECK( frozen_peak_range_membership(candidate_a,gammas)
               == frozen_peak_range_membership(candidate_b,gammas) );

  // A non-null empty policy freezes zero BR rows; it is not equivalent to "classify this seed".
  const std::vector<std::pair<double,double>> no_ranges;
  BOOST_CHECK( peak_ranges_for_frozen_solve(candidate_a_seed,&no_ranges).empty() );
  BOOST_CHECK( peak_ranges_for_frozen_solve(candidate_a_seed,nullptr) == candidate_a_seed );

  BOOST_CHECK( valid_frozen_peak_ranges(outer_ranges) );
  BOOST_CHECK( !valid_frozen_peak_ranges({{99.0,102.0},{101.0,103.0}}) );
}


BOOST_AUTO_TEST_CASE( selected_model_policy_is_exact_canonical_and_replayable )
{
  FrozenModelPolicy policy{
    { "Skew_1", 0.0 },
    { "REQ0_2", -0.0 },
    { "EAtt0(AD)", 10.0 }
  };
  const FrozenModelPolicy canonical = canonical_frozen_model_policy(policy);
  BOOST_REQUIRE_EQUAL( canonical.size(),3U );
  BOOST_CHECK_EQUAL( canonical[0].semantic_name,"EAtt0(AD)" );
  BOOST_CHECK_EQUAL( canonical[1].semantic_name,"REQ0_2" );
  BOOST_CHECK_EQUAL( canonical[2].semantic_name,"Skew_1" );

  const std::uint64_t baseline_hash = frozen_model_policy_hash(policy);
  std::reverse(policy.begin(),policy.end());
  BOOST_CHECK_EQUAL( frozen_model_policy_hash(policy),baseline_hash );

  FrozenModelPolicy changed_value = policy;
  changed_value[1].fixed_value = std::nextafter(changed_value[1].fixed_value,1.0);
  BOOST_CHECK_NE( frozen_model_policy_hash(changed_value),baseline_hash );
  FrozenModelPolicy changed_identity = policy;
  changed_identity[1].semantic_name += "-different";
  BOOST_CHECK_NE( frozen_model_policy_hash(changed_identity),baseline_hash );

  const std::vector<std::string> rebuilt_names{
      "EneOffset","Skew_1","EAtt0(AD)","REQ0_2"};
  const std::vector<size_t> indices
      = frozen_model_policy_indices(policy,rebuilt_names);
  BOOST_REQUIRE_EQUAL( indices.size(),3U );
  // Resolution follows canonical semantic order, never caller policy order.
  BOOST_CHECK_EQUAL( indices[0],2U );
  BOOST_CHECK_EQUAL( indices[1],3U );
  BOOST_CHECK_EQUAL( indices[2],1U );

  BOOST_CHECK_THROW( canonical_frozen_model_policy({
      {"REQ0_2",0.0},{"REQ0_2",1.0}}),std::invalid_argument );
  BOOST_CHECK_THROW( frozen_model_policy_indices(policy,{"Skew_1","EAtt0(AD)"}),
                     std::invalid_argument );
  BOOST_CHECK_THROW( frozen_model_policy_indices({{"Skew_1",0.0}},
                     {"Skew_1","Skew_1"}),std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( full_objective_ranking_is_order_independent )
{
  std::vector<DeterministicSearchScore> scores{
    { "z-default", 100.0, {1.0,2.0}, true },
    { "b-basin",    90.0, {3.0,4.0}, true },
    { "a-basin",    90.0 + 1.0e-12, {5.0,6.0}, true },
    { "failed",      1.0, {7.0,8.0}, false }
  };

  auto winner_name = []( const std::vector<DeterministicSearchScore> &input ){
    const std::vector<size_t> order = deterministic_search_order(input);
    BOOST_REQUIRE( !order.empty() );
    return input[order.front()].semantic_name;
  };

  BOOST_CHECK_EQUAL( winner_name(scores), "b-basin" );
  std::reverse(scores.begin(),scores.end());
  BOOST_CHECK_EQUAL( winner_name(scores), "b-basin" );
  std::rotate(scores.begin(),scores.begin()+1,scores.end());
  BOOST_CHECK_EQUAL( winner_name(scores), "b-basin" );
}

BOOST_AUTO_TEST_CASE( objective_order_is_strict_at_large_scale )
{
  // The former relative "tie" tolerance was 1.0 at this scale, so the semantic name could defeat
  // a genuinely lower full objective and the comparator's equivalence relation was non-transitive.
  std::vector<DeterministicSearchScore> scores{
    { "a-higher", 1000000001.0, {1.0}, true },
    { "z-lower",  1000000000.0, {2.0}, true },
    { "m-middle", 1000000000.5, {3.0}, true }
  };
  const std::vector<size_t> order = deterministic_search_order(scores);
  BOOST_REQUIRE_EQUAL( order.size(), 3U );
  BOOST_CHECK_EQUAL( scores[order[0]].semantic_name, "z-lower" );
  BOOST_CHECK_EQUAL( scores[order[1]].semantic_name, "m-middle" );
  BOOST_CHECK_EQUAL( scores[order[2]].semantic_name, "a-higher" );

  const std::vector<BackwardEliminationScore> removals{
    { "a-higher", 1000000001.0, true },
    { "z-lower",  1000000000.0, true }
  };
  const std::optional<size_t> selected
      = select_backward_elimination_removal(removals,1000000000.0,1.0);
  BOOST_REQUIRE( selected );
  BOOST_CHECK_EQUAL( removals[*selected].semantic_key, "z-lower" );
}

BOOST_AUTO_TEST_CASE( duplicate_basins_are_collapsed_after_ranking )
{
  const std::vector<DeterministicSearchScore> scores{
    { "better-copy", 10.0, {1.0,2.0,3.0}, true },
    { "same-basin",  10.1, {1.0 + 1.0e-8,2.0,3.0}, true },
    { "other-basin", 11.0, {8.0,9.0,10.0}, true }
  };
  const std::vector<size_t> order = deterministic_search_order(scores);
  BOOST_REQUIRE_EQUAL( order.size(), 2U );
  BOOST_CHECK_EQUAL( scores[order[0]].semantic_name, "better-copy" );
  BOOST_CHECK_EQUAL( scores[order[1]].semantic_name, "other-basin" );
}

BOOST_AUTO_TEST_CASE( backward_elimination_evaluates_all_removals_from_parent )
{
  const double parent = 50.0;
  std::vector<BackwardEliminationScore> trials{
    { "z-first-in-generation-order", 50.8, true },
    { "a-best-later-candidate",      49.9, true },
    { "outside-delta",               51.1, true },
    { "failed",                      1.0, false }
  };

  std::optional<size_t> selected = select_backward_elimination_removal(trials,parent,1.0);
  BOOST_REQUIRE( selected );
  BOOST_CHECK_EQUAL( trials[*selected].semantic_key, "a-best-later-candidate" );

  std::reverse(trials.begin(),trials.end());
  selected = select_backward_elimination_removal(trials,parent,1.0);
  BOOST_REQUIRE( selected );
  BOOST_CHECK_EQUAL( trials[*selected].semantic_key, "a-best-later-candidate" );
}

BOOST_AUTO_TEST_CASE( backward_elimination_uses_semantics_only_for_exact_ties )
{
  const std::vector<BackwardEliminationScore> distinct_trials{
    { "curve:Pu:term:2", 20.0, true },
    { "curve:Am:term:2", 20.0 + 1.0e-12, true }
  };
  std::optional<size_t> selected
      = select_backward_elimination_removal(distinct_trials,20.0,1.0);
  BOOST_REQUIRE( selected );
  BOOST_CHECK_EQUAL( distinct_trials[*selected].semantic_key, "curve:Pu:term:2" );

  const std::vector<BackwardEliminationScore> tied_trials{
    { "curve:Pu:term:2", 20.0, true },
    { "curve:Am:term:2", 20.0, true }
  };
  selected = select_backward_elimination_removal(tied_trials,20.0,1.0);
  BOOST_REQUIRE( selected );
  BOOST_CHECK_EQUAL( tied_trials[*selected].semantic_key, "curve:Am:term:2" );
}

BOOST_AUTO_TEST_CASE( weakest_svd_direction_is_semantically_remapped )
{
  const std::vector<std::string> previous_names{"Act:Pu239","Age:Pu239","Act:Pu240"};
  const std::vector<double> previous_direction{0.25,-2.0,1.5};
  const std::vector<std::string> current_names{"Act:Pu240","Act:Pu239","Age:Pu239"};
  const std::vector<double> mapped = remap_semantic_direction(
      previous_names,previous_direction,current_names );
  BOOST_REQUIRE_EQUAL( mapped.size(), 3U );
  BOOST_CHECK_EQUAL( mapped[0], 1.5 );
  BOOST_CHECK_EQUAL( mapped[1], 0.25 );
  BOOST_CHECK_EQUAL( mapped[2], -2.0 );

  BOOST_CHECK( remap_semantic_direction(previous_names,{1.0},current_names).empty() );
}

BOOST_AUTO_TEST_CASE( weakest_svd_seed_moves_collinearly_and_inward )
{
  const std::vector<double> seed{0.95,0.80,7.0};
  const std::vector<double> direction{1.0,2.0,100.0};
  const std::vector<std::optional<double>> lower{0.0,0.0,std::nullopt};
  const std::vector<std::optional<double>> upper{1.0,1.0,std::nullopt};
  const std::vector<char> constant{0,0,1};

  const WeakDirectionSeed moved = perturb_seed_along_weak_direction(
      seed,direction,lower,upper,constant );
  BOOST_REQUIRE( moved.applied );
  BOOST_CHECK_LT( moved.signed_step, 0.0 ); //negative has more room than the upper-bound direction
  BOOST_CHECK_GE( moved.values[0], 0.0 );
  BOOST_CHECK_LE( moved.values[0], 1.0 );
  BOOST_CHECK_GE( moved.values[1], 0.0 );
  BOOST_CHECK_LE( moved.values[1], 1.0 );
  BOOST_CHECK_EQUAL( moved.values[2], seed[2] );
  BOOST_CHECK_CLOSE_FRACTION( (moved.values[1] - seed[1])/(moved.values[0] - seed[0]),
                             2.0, 1.0e-12 );
}


BOOST_AUTO_TEST_CASE( weak_direction_sign_can_be_forced_against_the_room_heuristic )
{
  // Deliberately asymmetric room: the seed sits near its upper bounds, so the room heuristic walks
  // downward.  A profile likelihood that found a better point has already shown the upward half of
  // this direction is worth spending a candidate slot on, so the sign must be forceable.
  const std::vector<double> seed{0.95,0.80};
  const std::vector<double> direction{1.0,2.0};
  const std::vector<std::optional<double>> lower{0.0,0.0};
  const std::vector<std::optional<double>> upper{1.0,1.0};
  const std::vector<char> constant{0,0};

  const WeakDirectionSeed automatic = perturb_seed_along_weak_direction(
      seed,direction,lower,upper,constant );
  const WeakDirectionSeed defaulted = perturb_seed_along_weak_direction(
      seed,direction,lower,upper,constant,WeakDirectionSign::MoreFeasibleRoom );
  const WeakDirectionSeed positive = perturb_seed_along_weak_direction(
      seed,direction,lower,upper,constant,WeakDirectionSign::Positive );
  const WeakDirectionSeed negative = perturb_seed_along_weak_direction(
      seed,direction,lower,upper,constant,WeakDirectionSign::Negative );

  // Omitting the policy must remain byte-identical to the historical behavior.
  BOOST_REQUIRE( automatic.applied && defaulted.applied );
  BOOST_CHECK_EQUAL( automatic.signed_step, defaulted.signed_step );
  BOOST_CHECK_EQUAL( automatic.values[0], defaulted.values[0] );
  BOOST_CHECK_EQUAL( automatic.values[1], defaulted.values[1] );
  BOOST_CHECK_LT( automatic.signed_step, 0.0 );
  BOOST_CHECK_EQUAL( automatic.signed_step, negative.signed_step );

  BOOST_REQUIRE( positive.applied && negative.applied );
  BOOST_CHECK_GT( positive.signed_step, 0.0 );
  BOOST_CHECK_LT( negative.signed_step, 0.0 );
  BOOST_CHECK_GT( positive.values[0], seed[0] );
  BOOST_CHECK_LT( negative.values[0], seed[0] );

  // Both signs stay strictly feasible and collinear with the requested direction.
  for( const WeakDirectionSeed &result : {positive,negative} )
  {
    for( size_t i = 0; i < result.values.size(); ++i )
    {
      BOOST_CHECK_GE( result.values[i], *lower[i] );
      BOOST_CHECK_LE( result.values[i], *upper[i] );
    }
    BOOST_CHECK_CLOSE_FRACTION( (result.values[1] - seed[1])/(result.values[0] - seed[0]),
                                2.0, 1.0e-12 );
  }
}


BOOST_AUTO_TEST_CASE( opposite_weak_direction_is_named_but_outside_the_native_matrix )
{
  // The native matrix is the hard cap of 16.  The Opposite variant exists only for the one permitted
  // profile-triggered baseline reselection, and may occupy only a slot the applicability screen left
  // unused; letting it into search_seed_variants() would silently raise the cap to 17.
  const auto variants = search_seed_variants();
  BOOST_CHECK( std::find(variants.begin(),variants.end(),
                         SearchSeedVariant::WeakDirectionCheckerboardOpposite)
               == variants.end() );

  const std::string opposite
      = search_seed_variant_name(SearchSeedVariant::WeakDirectionCheckerboardOpposite);
  BOOST_CHECK( !opposite.empty() );
  BOOST_CHECK( opposite != "unknown" );
  for( const SearchSeedVariant variant : variants )
    BOOST_CHECK( opposite != search_seed_variant_name(variant) );

  BOOST_CHECK( weak_direction_sign(SearchSeedVariant::WeakDirectionCheckerboard)
               == WeakDirectionSign::MoreFeasibleRoom );
  BOOST_CHECK( weak_direction_sign(SearchSeedVariant::WeakDirectionCheckerboardOpposite)
               == WeakDirectionSign::Negative );
  // Every other variant leaves the historical heuristic alone.
  for( const SearchSeedVariant variant : variants )
    BOOST_CHECK( weak_direction_sign(variant) == WeakDirectionSign::MoreFeasibleRoom );
}
