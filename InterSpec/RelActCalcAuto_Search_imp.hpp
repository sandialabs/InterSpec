#ifndef InterSpec_RelActCalcAuto_Search_imp_hpp
#define InterSpec_RelActCalcAuto_Search_imp_hpp

/* InterSpec: deterministic basin-search and backward-elimination helpers.

   This header intentionally contains only small, data-independent ordering utilities.  Keeping the
   comparison policy here makes the production solver and focused unit tests exercise exactly the
   same NaN handling, objective tie-breaking, and basin de-duplication rules.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace RelActCalcAutoImp
{

/** A nominal gamma identity used to fingerprint which lines receive each frozen branching-ratio
    nuisance parameter.  The source id must be semantic (never a pointer or caller-local index). */
struct FrozenPeakRangeGamma
{
  std::string semantic_source;
  double energy = 0.0;
  int gamma_type = 0;
};


/** One degree of freedom removed by selected-basin backward elimination.

 The name is the canonical optimizer-parameter identity (after source canonicalization), and the
 value is the exact identity value at which the selected model holds that parameter.  Keeping the
 value in the policy is important: a semantic warm start is only a seed and must never be relied on
 to reproduce a fixed model coefficient. */
struct FrozenModelParameter
{
  std::string semantic_name;
  double fixed_value = 0.0;
};

using FrozenModelPolicy = std::vector<FrozenModelParameter>;


/** Validate and canonicalize a selected-model policy.  Policy order is presentation state; exact
 names and value bits define the mathematical model. */
inline FrozenModelPolicy canonical_frozen_model_policy( FrozenModelPolicy policy )
{
  for( const FrozenModelParameter &parameter : policy )
    if( parameter.semantic_name.empty() || !std::isfinite(parameter.fixed_value) )
      throw std::invalid_argument( "Invalid frozen selected-model parameter" );

  std::sort( policy.begin(),policy.end(),
    []( const FrozenModelParameter &lhs, const FrozenModelParameter &rhs ) {
      return lhs.semantic_name < rhs.semantic_name;
    } );
  for( size_t index = 1; index < policy.size(); ++index )
    if( policy[index - 1].semantic_name == policy[index].semantic_name )
      throw std::invalid_argument( "Duplicate frozen selected-model parameter '"
                                   + policy[index].semantic_name + "'" );
  return policy;
}


/** Resolve the canonical policy into a rebuilt optimizer layout.  Every selected-model parameter
 must exist exactly once; silently dropping a removed degree of freedom would change the profile
 likelihood's nuisance model. */
inline std::vector<size_t> frozen_model_policy_indices(
    const FrozenModelPolicy &policy, const std::vector<std::string> &parameter_names )
{
  const FrozenModelPolicy canonical = canonical_frozen_model_policy(policy);
  std::vector<size_t> indices;
  indices.reserve(canonical.size());
  for( const FrozenModelParameter &parameter : canonical )
  {
    const auto first = std::find( parameter_names.begin(),parameter_names.end(),
                                  parameter.semantic_name );
    if( first == parameter_names.end() )
      throw std::invalid_argument( "Frozen selected-model parameter '" + parameter.semantic_name
                                   + "' is absent from the rebuilt layout" );
    const auto duplicate = std::find( std::next(first),parameter_names.end(),
                                      parameter.semantic_name );
    if( duplicate != parameter_names.end() )
      throw std::invalid_argument( "Frozen selected-model parameter '" + parameter.semantic_name
                                   + "' is ambiguous in the rebuilt layout" );
    indices.push_back( static_cast<size_t>(first - parameter_names.begin()) );
  }
  return indices;
}


/** Stable exact-bit fingerprint of the selected model. */
inline std::uint64_t frozen_model_policy_hash( const FrozenModelPolicy &input_policy )
{
  const FrozenModelPolicy policy = canonical_frozen_model_policy(input_policy);
  std::uint64_t hash = UINT64_C(14695981039346656037);
  const auto add_byte = [&hash]( const std::uint8_t byte ) {
    hash ^= byte;
    hash *= UINT64_C(1099511628211);
  };
  const auto add_u64 = [&add_byte]( const std::uint64_t value ) {
    for( unsigned int shift = 0; shift < 64; shift += 8 )
      add_byte( static_cast<std::uint8_t>((value >> shift) & UINT64_C(0xff)) );
  };
  const auto add_string = [&add_u64,&add_byte]( const std::string &value ) {
    add_u64(value.size());
    for( const unsigned char byte : value )
      add_byte(byte);
  };

  add_string( "RelActAuto/frozen-selected-model/v1" );
  add_u64(policy.size());
  for( const FrozenModelParameter &parameter : policy )
  {
    add_string(parameter.semantic_name);
    std::uint64_t bits = 0;
    static_assert( sizeof(bits) == sizeof(parameter.fixed_value),
                   "Unexpected double representation" );
    std::memcpy( &bits,&parameter.fixed_value,sizeof(bits) );
    add_u64(bits);
  }
  return hash;
}


/** Validate the exact interval policy consumed by the branching-ratio nuisance parameters.  The
    generated policy is sorted and disjoint; recursively rebuilt candidates must preserve those
    same properties as well as every endpoint bit. */
inline bool valid_frozen_peak_ranges(
    const std::vector<std::pair<double,double>> &ranges )
{
  for( size_t index = 0; index < ranges.size(); ++index )
  {
    const std::pair<double,double> &range = ranges[index];
    if( !std::isfinite(range.first) || !std::isfinite(range.second)
        || !(range.first < range.second) )
      return false;
    if( index && (ranges[index - 1].second > range.first) )
      return false;
  }
  return true;
}


/** Select the peak-range policy for one solve.  Once an outer problem supplies a frozen policy,
    candidate-specific seed classification is deliberately ignored and the exact endpoint values
    are copied.  A non-null empty policy is meaningful: it freezes the absence of nuisance rows. */
inline std::vector<std::pair<double,double>> peak_ranges_for_frozen_solve(
    const std::vector<std::pair<double,double>> &seed_classified_ranges,
    const std::vector<std::pair<double,double>> *frozen_ranges )
{
  return frozen_ranges ? *frozen_ranges : seed_classified_ranges;
}


/** Return the first matching range plus one for every gamma (zero means no range).  First-match
    behavior intentionally mirrors the production evaluator, including a line exactly on a shared
    endpoint of two touching ranges. */
inline std::vector<size_t> frozen_peak_range_membership(
    const std::vector<std::pair<double,double>> &ranges,
    const std::vector<FrozenPeakRangeGamma> &gammas )
{
  if( !valid_frozen_peak_ranges(ranges) )
    throw std::invalid_argument( "Invalid frozen branching-ratio peak ranges" );

  std::vector<size_t> membership( gammas.size(), 0 );
  for( size_t gamma = 0; gamma < gammas.size(); ++gamma )
  {
    if( !std::isfinite(gammas[gamma].energy) )
      throw std::invalid_argument( "Non-finite gamma energy in frozen peak-range policy" );
    for( size_t range = 0; range < ranges.size(); ++range )
    {
      if( (gammas[gamma].energy >= ranges[range].first)
          && (gammas[gamma].energy <= ranges[range].second) )
      {
        membership[gamma] = range + 1;
        break;
      }
    }
  }
  return membership;
}


/** Stable fingerprint of exact range endpoints and semantic gamma membership.  Gamma inputs are
    sorted internally, so caller/source order cannot affect the result. */
inline std::uint64_t frozen_peak_range_policy_hash(
    const std::vector<std::pair<double,double>> &ranges,
    std::vector<FrozenPeakRangeGamma> gammas )
{
  if( !valid_frozen_peak_ranges(ranges) )
    throw std::invalid_argument( "Invalid frozen branching-ratio peak ranges" );
  for( const FrozenPeakRangeGamma &gamma : gammas )
    if( !std::isfinite(gamma.energy) )
      throw std::invalid_argument( "Non-finite gamma energy in frozen peak-range policy" );

  std::sort( gammas.begin(), gammas.end(),
    []( const FrozenPeakRangeGamma &lhs, const FrozenPeakRangeGamma &rhs ) {
      return std::tie(lhs.semantic_source,lhs.energy,lhs.gamma_type)
           < std::tie(rhs.semantic_source,rhs.energy,rhs.gamma_type);
    } );
  const std::vector<size_t> membership = frozen_peak_range_membership(ranges,gammas);

  std::uint64_t hash = UINT64_C(14695981039346656037);
  const auto add_byte = [&hash]( const std::uint8_t byte ) {
    hash ^= byte;
    hash *= UINT64_C(1099511628211);
  };
  const auto add_u64 = [&add_byte]( const std::uint64_t value ) {
    for( unsigned int shift = 0; shift < 64; shift += 8 )
      add_byte( static_cast<std::uint8_t>((value >> shift) & UINT64_C(0xff)) );
  };
  const auto add_double = [&add_u64]( const double value ) {
    std::uint64_t bits = 0;
    static_assert( sizeof(bits) == sizeof(value), "Unexpected double representation" );
    std::memcpy( &bits, &value, sizeof(bits) );
    add_u64(bits);
  };
  const auto add_string = [&add_u64,&add_byte]( const std::string &value ) {
    add_u64(value.size());
    for( const unsigned char byte : value )
      add_byte(byte);
  };

  add_string( "RelActAuto/frozen-branching-ratio-ranges/v1" );
  add_u64(ranges.size());
  for( const std::pair<double,double> &range : ranges )
  {
    add_double(range.first);
    add_double(range.second);
  }
  add_u64(gammas.size());
  for( size_t gamma = 0; gamma < gammas.size(); ++gamma )
  {
    add_string(gammas[gamma].semantic_source);
    add_double(gammas[gamma].energy);
    add_u64(static_cast<std::uint64_t>(gammas[gamma].gamma_type));
    add_u64(membership[gamma]);
  }
  return hash;
}

/** Which way to travel along a weak singular direction.

 `MoreFeasibleRoom` is the historical behavior and stays the default: it is the right choice for an
 ordinary exploratory candidate, because the sign with more room usually yields the larger useful
 displacement.  Explicit signs exist because a profile likelihood which discovered a lower objective
 has already demonstrated that *both* signs of the weak direction are physically interesting, and the
 room heuristic then silently discards half of the search space. */
enum class WeakDirectionSign : int
{
  MoreFeasibleRoom = 0,
  Positive,
  Negative
};


/** Named, semantic starting points used by the deterministic basin search.  The numeric values are
    persisted nowhere; their order is nevertheless fixed so a search is reproducible. */
enum class SearchSeedVariant : int
{
  Default = 0,
  EmMedianAttribution,
  PreEmUnweighted,
  AgeLowerQuartile,
  AgeUpperQuartile,
  CalibrationLowerQuartile,
  CalibrationUpperQuartile,
  FwhmLowerQuartile,
  FwhmUpperQuartile,
  ShieldingLowerQuartile,
  ShieldingUpperQuartile,
  EfficiencyPositiveShape,
  EfficiencyNegativeShape,
  ActivitySemanticSplitA,
  ActivitySemanticSplitB,
  WeakDirectionCheckerboard,

  /** The mirror of `WeakDirectionCheckerboard`.

   Deliberately absent from `search_seed_variants()`: the native matrix is capped at 16 and this
   variant is only justified once a profile likelihood has demonstrated that the retained baseline
   is not optimal.  It is offered to the one permitted baseline reselection through the restart
   context, and only into a candidate slot the applicable native matrix left unused. */
  WeakDirectionCheckerboardOpposite
};

/** The native candidate matrix.  Exactly 16 entries, which is also the hard cap on the number of
    candidates any single search may polish; extra restart-context seeds may only occupy slots the
    applicability screen left unused. */
inline constexpr std::array<SearchSeedVariant,16> search_seed_variants()
{
  return {{
    SearchSeedVariant::Default,
    SearchSeedVariant::EmMedianAttribution,
    SearchSeedVariant::PreEmUnweighted,
    SearchSeedVariant::AgeLowerQuartile,
    SearchSeedVariant::AgeUpperQuartile,
    SearchSeedVariant::CalibrationLowerQuartile,
    SearchSeedVariant::CalibrationUpperQuartile,
    SearchSeedVariant::FwhmLowerQuartile,
    SearchSeedVariant::FwhmUpperQuartile,
    SearchSeedVariant::ShieldingLowerQuartile,
    SearchSeedVariant::ShieldingUpperQuartile,
    SearchSeedVariant::EfficiencyPositiveShape,
    SearchSeedVariant::EfficiencyNegativeShape,
    SearchSeedVariant::ActivitySemanticSplitA,
    SearchSeedVariant::ActivitySemanticSplitB,
    SearchSeedVariant::WeakDirectionCheckerboard
  }};
}

inline const char *search_seed_variant_name( const SearchSeedVariant variant )
{
  switch( variant )
  {
    case SearchSeedVariant::Default:                    return "default";
    case SearchSeedVariant::EmMedianAttribution:        return "em-median-attribution";
    case SearchSeedVariant::PreEmUnweighted:            return "pre-em-unweighted";
    case SearchSeedVariant::AgeLowerQuartile:           return "age-lower-quartile";
    case SearchSeedVariant::AgeUpperQuartile:           return "age-upper-quartile";
    case SearchSeedVariant::CalibrationLowerQuartile:   return "calibration-lower-quartile";
    case SearchSeedVariant::CalibrationUpperQuartile:   return "calibration-upper-quartile";
    case SearchSeedVariant::FwhmLowerQuartile:          return "fwhm-lower-quartile";
    case SearchSeedVariant::FwhmUpperQuartile:          return "fwhm-upper-quartile";
    case SearchSeedVariant::ShieldingLowerQuartile:     return "shielding-lower-quartile";
    case SearchSeedVariant::ShieldingUpperQuartile:     return "shielding-upper-quartile";
    case SearchSeedVariant::EfficiencyPositiveShape:   return "efficiency-positive-shape";
    case SearchSeedVariant::EfficiencyNegativeShape:   return "efficiency-negative-shape";
    case SearchSeedVariant::ActivitySemanticSplitA:    return "activity-semantic-split-a";
    case SearchSeedVariant::ActivitySemanticSplitB:    return "activity-semantic-split-b";
    case SearchSeedVariant::WeakDirectionCheckerboard: return "weak-direction-checkerboard";
    case SearchSeedVariant::WeakDirectionCheckerboardOpposite:
                                                       return "weak-direction-checkerboard-opposite";
  }
  return "unknown";
}


/** The travel direction each weak-direction variant asks for. */
inline WeakDirectionSign weak_direction_sign( const SearchSeedVariant variant )
{
  return (variant == SearchSeedVariant::WeakDirectionCheckerboardOpposite)
       ? WeakDirectionSign::Negative : WeakDirectionSign::MoreFeasibleRoom;
}

struct DeterministicSearchScore
{
  std::string semantic_name;
  double full_objective = std::numeric_limits<double>::infinity();
  std::vector<double> basin_signature;
  bool usable = false;
};


/** Map an ambient direction between identical semantic layouts without relying on
    caller/source order.  Repeated names are matched by occurrence, which mirrors
    the deterministic canonical parameter construction.  An empty result means
    the layouts are not semantically equivalent. */
inline std::vector<double> remap_semantic_direction(
    const std::vector<std::string> &previous_names,
    const std::vector<double> &previous_direction,
    const std::vector<std::string> &current_names )
{
  if( previous_names.size() != previous_direction.size() )
    return {};

  std::vector<double> answer( current_names.size(), 0.0 );
  std::vector<char> used( previous_names.size(), 0 );
  for( size_t current = 0; current < current_names.size(); ++current )
  {
    size_t found = previous_names.size();
    for( size_t previous = 0; previous < previous_names.size(); ++previous )
    {
      if( !used[previous] && (previous_names[previous] == current_names[current]) )
      {
        found = previous;
        break;
      }
    }
    if( (found == previous_names.size()) || !std::isfinite(previous_direction[found]) )
      return {};
    used[found] = 1;
    answer[current] = previous_direction[found];
  }
  return answer;
}


struct WeakDirectionSeed
{
  std::vector<double> values;
  bool applied = false;
  double signed_step = 0.0;
};


/** Move a seed along one actual SVD direction while respecting every parameter
    bound.  By default the sign with more feasible room is selected
    deterministically; `sign_policy` can instead force either sign.  The
    intended displacement is one quarter of the largest parameter's natural
    bound/value scale, and the result remains 20% of the available room inside a
    blocking bound so the optimizer sees an inward derivative. */
inline WeakDirectionSeed perturb_seed_along_weak_direction(
    const std::vector<double> &seed,
    const std::vector<double> &direction,
    const std::vector<std::optional<double>> &lower_bounds,
    const std::vector<std::optional<double>> &upper_bounds,
    const std::vector<char> &constant_parameters,
    const WeakDirectionSign sign_policy = WeakDirectionSign::MoreFeasibleRoom )
{
  WeakDirectionSeed answer;
  answer.values = seed;
  if( seed.empty() || (direction.size() != seed.size())
      || (lower_bounds.size() != seed.size()) || (upper_bounds.size() != seed.size())
      || (!constant_parameters.empty() && (constant_parameters.size() != seed.size())) )
    return answer;

  std::vector<double> scales( seed.size(), 1.0 );
  double max_scaled_direction = 0.0;
  for( size_t i = 0; i < seed.size(); ++i )
  {
    if( !std::isfinite(seed[i]) || !std::isfinite(direction[i]) )
      return answer;
    if( !constant_parameters.empty() && constant_parameters[i] )
      continue;

    if( lower_bounds[i] && upper_bounds[i] && (*upper_bounds[i] > *lower_bounds[i]) )
      scales[i] = *upper_bounds[i] - *lower_bounds[i];
    else
    {
      scales[i] = 1.0 + std::fabs(seed[i]);
      if( lower_bounds[i] ) scales[i] = (std::max)(scales[i],std::fabs(*lower_bounds[i]));
      if( upper_bounds[i] ) scales[i] = (std::max)(scales[i],std::fabs(*upper_bounds[i]));
    }
    max_scaled_direction = (std::max)(max_scaled_direction,
                                      std::fabs(direction[i])/scales[i]);
  }
  if( !(max_scaled_direction > 0.0) || !std::isfinite(max_scaled_direction) )
    return answer;

  const double target_step = 0.25/max_scaled_direction;
  const double component_floor = 1.0e-12*max_scaled_direction;
  const auto feasible_step = [&]( const double sign ){
    double maximum = std::numeric_limits<double>::infinity();
    for( size_t i = 0; i < seed.size(); ++i )
    {
      if( (!constant_parameters.empty() && constant_parameters[i])
          || (std::fabs(direction[i])/scales[i] <= component_floor) )
        continue;
      const double delta = sign*direction[i];
      if( (delta > 0.0) && upper_bounds[i] )
      {
        const double room = (std::max)(0.0,*upper_bounds[i] - seed[i]);
        maximum = (std::min)(maximum,0.8*room/delta);
      }else if( (delta < 0.0) && lower_bounds[i] )
      {
        const double room = (std::max)(0.0,seed[i] - *lower_bounds[i]);
        maximum = (std::min)(maximum,0.8*room/(-delta));
      }
    }
    return (std::min)(target_step,maximum);
  };

  const double positive_step = feasible_step(1.0);
  const double negative_step = feasible_step(-1.0);
  double sign = 1.0;
  switch( sign_policy )
  {
    case WeakDirectionSign::Positive: sign = 1.0; break;
    case WeakDirectionSign::Negative: sign = -1.0; break;
    case WeakDirectionSign::MoreFeasibleRoom:
      sign = (negative_step > positive_step*(1.0 + 1.0e-12)) ? -1.0 : 1.0;
      break;
  }
  const double step = (sign > 0.0) ? positive_step : negative_step;
  if( !(step > 1.0e-12*target_step) || !std::isfinite(step) )
    return answer;

  bool changed = false;
  for( size_t i = 0; i < seed.size(); ++i )
  {
    if( !constant_parameters.empty() && constant_parameters[i] )
      continue;
    double value = seed[i] + sign*step*direction[i];
    if( lower_bounds[i] ) value = (std::max)(value,*lower_bounds[i]);
    if( upper_bounds[i] ) value = (std::min)(value,*upper_bounds[i]);
    if( !std::isfinite(value) )
      return WeakDirectionSeed{seed,false,0.0};
    changed = changed || (value != seed[i]);
    answer.values[i] = value;
  }
  answer.applied = changed;
  answer.signed_step = changed ? sign*step : 0.0;
  return answer;
}

inline bool search_objectives_tied( const double lhs, const double rhs )
{
  // Approximate equality is not transitive and therefore cannot safely be used by std::sort.  More
  // importantly, a relative tolerance at a large objective (for example 1e9) can call a real
  // delta-chi2 of one a "tie" and let a semantic name defeat the lower physical objective.  Search
  // candidates are evaluated by the same frozen scalar functor, so only numerical equality is a
  // tie; every representable decrease wins before semantic ordering is considered.
  return lhs == rhs;
}

inline bool same_search_basin( const std::vector<double> &lhs, const std::vector<double> &rhs )
{
  if( lhs.size() != rhs.size() )
    return false;
  double max_scaled_delta = 0.0;
  for( size_t i = 0; i < lhs.size(); ++i )
  {
    if( !std::isfinite(lhs[i]) || !std::isfinite(rhs[i]) )
      return false;
    const double scale = 1.0 + (std::max)(std::fabs(lhs[i]),std::fabs(rhs[i]));
    max_scaled_delta = (std::max)(max_scaled_delta, std::fabs(lhs[i] - rhs[i])/scale);
  }
  return max_scaled_delta <= 1.0e-6;
}

/** Return candidate indices best-first, with unusable/non-finite candidates omitted and duplicate
    basins collapsed.  Meaningful objective differences rank first; numerical ties use the semantic
    name, never generation or caller order. */
inline std::vector<size_t> deterministic_search_order( const std::vector<DeterministicSearchScore> &scores )
{
  std::vector<size_t> order;
  for( size_t i = 0; i < scores.size(); ++i )
    if( scores[i].usable && std::isfinite(scores[i].full_objective) )
      order.push_back(i);

  std::sort( order.begin(), order.end(), [&]( const size_t lhs, const size_t rhs ){
    const DeterministicSearchScore &a = scores[lhs], &b = scores[rhs];
    if( a.full_objective < b.full_objective )
      return a.full_objective < b.full_objective;
    if( b.full_objective < a.full_objective )
      return false;
    if( a.semantic_name != b.semantic_name )
      return a.semantic_name < b.semantic_name;
    return lhs < rhs;
  } );

  std::vector<size_t> unique;
  for( const size_t candidate : order )
  {
    bool duplicate = false;
    for( const size_t retained : unique )
      duplicate = duplicate || same_search_basin(scores[candidate].basin_signature,
                                                  scores[retained].basin_signature);
    if( !duplicate )
      unique.push_back(candidate);
  }
  return unique;
}

struct BackwardEliminationScore
{
  std::string semantic_key;
  double full_objective = std::numeric_limits<double>::infinity();
  bool solved = false;
};

/** Select one accepted one-DOF removal after every removal from the same parent was evaluated.
    Lower objective wins; objective ties use the stable semantic key. */
inline std::optional<size_t> select_backward_elimination_removal(
    const std::vector<BackwardEliminationScore> &trials,
    const double parent_full_objective,
    const double max_delta )
{
  std::optional<size_t> best;
  for( size_t i = 0; i < trials.size(); ++i )
  {
    const BackwardEliminationScore &trial = trials[i];
    if( !trial.solved || !std::isfinite(trial.full_objective)
        || (trial.full_objective - parent_full_objective > max_delta) )
      continue;
    if( !best )
    {
      best = i;
      continue;
    }
    const BackwardEliminationScore &incumbent = trials[*best];
    if( (trial.full_objective < incumbent.full_objective)
        || (search_objectives_tied(trial.full_objective,incumbent.full_objective)
            && (trial.semantic_key < incumbent.semantic_key)) )
      best = i;
  }
  return best;
}

}//namespace RelActCalcAutoImp

#endif //InterSpec_RelActCalcAuto_Search_imp_hpp
