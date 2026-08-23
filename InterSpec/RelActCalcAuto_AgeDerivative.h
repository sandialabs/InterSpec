/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
*/

#ifndef InterSpec_RelActCalcAuto_AgeDerivative_h
#define InterSpec_RelActCalcAuto_AgeDerivative_h

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <deque>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace RelActCalcAutoImp
{

/** A pointer-free identity and yield for one line in an aged source.  Gamma energies come from
 * immutable decay-database records, so their exact IEEE bits remain stable when the transition
 * that supplies the largest share of a combined line changes with age. */
struct AgeGammaYield
{
  std::uint64_t identity = 0;
  double yield = 0.0;
};

inline std::uint64_t exact_double_identity( const double value )
{
  static_assert( sizeof(double) == sizeof(std::uint64_t), "Unexpected double representation." );
  std::uint64_t bits = 0;
  std::memcpy( &bits, &value, sizeof(bits) );
  return bits;
}

inline std::uint64_t age_gamma_energy_identity( const double energy )
{
  return exact_double_identity(energy);
}

/** Complete semantic key for the production fitted-age derivative cache.  Exact bits are required:
 * adjacent representable ages must not share a derivative, and neither may analyses with different
 * physical/input bounds or exclusion sets. */
struct AgeDerivativeCacheKey
{
  const void *source = nullptr;
  std::uint64_t age_bits = 0;
  std::vector<std::uint64_t> excluded_gamma_bits;
  std::uint64_t lower_bound_bits = 0;
  std::uint64_t upper_bound_bits = 0;
  std::uint64_t characteristic_time_bits = 0;

  AgeDerivativeCacheKey( const void * const source_in, const double age,
                         const std::vector<double> &excluded_gammas,
                         const double lower_bound, const double upper_bound,
                         const double characteristic_time )
    : source( source_in ),
      age_bits( exact_double_identity(age) ),
      excluded_gamma_bits{},
      lower_bound_bits( exact_double_identity(lower_bound) ),
      upper_bound_bits( exact_double_identity(upper_bound) ),
      characteristic_time_bits( exact_double_identity(characteristic_time) )
  {
    excluded_gamma_bits.reserve( excluded_gammas.size() );
    for( const double energy : excluded_gammas )
      excluded_gamma_bits.push_back( exact_double_identity(energy) );
    std::sort( excluded_gamma_bits.begin(), excluded_gamma_bits.end() );
    excluded_gamma_bits.erase(
        std::unique(excluded_gamma_bits.begin(), excluded_gamma_bits.end()),
        excluded_gamma_bits.end() );
  }
};

struct AgeDerivativeCacheKeyLess
{
  bool operator()( const AgeDerivativeCacheKey &lhs, const AgeDerivativeCacheKey &rhs ) const
  {
    const std::less<const void *> pointer_less;
    if( pointer_less(lhs.source, rhs.source) )
      return true;
    if( pointer_less(rhs.source, lhs.source) )
      return false;
    if( lhs.age_bits != rhs.age_bits )
      return lhs.age_bits < rhs.age_bits;
    if( lhs.excluded_gamma_bits != rhs.excluded_gamma_bits )
      return lhs.excluded_gamma_bits < rhs.excluded_gamma_bits;
    if( lhs.lower_bound_bits != rhs.lower_bound_bits )
      return lhs.lower_bound_bits < rhs.lower_bound_bits;
    if( lhs.upper_bound_bits != rhs.upper_bound_bits )
      return lhs.upper_bound_bits < rhs.upper_bound_bits;
    return lhs.characteristic_time_bits < rhs.characteristic_time_bits;
  }
};

/** Equality for the derivative cache key.  Every field is an exact bit identity, so this is exact
    equality by construction - see AgeDerivativeCacheKeyLess for why each field must participate. */
inline bool operator==( const AgeDerivativeCacheKey &lhs, const AgeDerivativeCacheKey &rhs )
{
  return (lhs.source == rhs.source)
      && (lhs.age_bits == rhs.age_bits)
      && (lhs.lower_bound_bits == rhs.lower_bound_bits)
      && (lhs.upper_bound_bits == rhs.upper_bound_bits)
      && (lhs.characteristic_time_bits == rhs.characteristic_time_bits)
      && (lhs.excluded_gamma_bits == rhs.excluded_gamma_bits);
}


/** A bounded most-recent-first memoization cache.

 Decay/aging results were previously memoized in an unbounded `std::map`.  Because the key holds the
 *exact* fitted age, a fit that varies an age mints a new entry on essentially every objective
 evaluation, and with a 50000-iteration budget the map grows without limit.  That is survivable for
 one solve, but the deterministic candidate search keeps up to sixteen candidate solutions alive at
 once and profiling retains dozens of conditional solutions, each holding its own cost functor and
 therefore its own cache - which is how a 30 MB cache became tens of gigabytes and took the machine
 down.

 The access pattern is strong temporal locality rather than long-term reuse: a DynamicAutoDiff pass
 asks for the same physical age across many parameter lanes and ROIs in immediate succession.  A
 handful of recent entries therefore captures nearly all of the hit rate that the old unbounded map
 achieved.  Note the capacity cannot be as small as the per-source locality suggests: the key
 includes the source, so one evaluation round queries one entry per nuclide before any repeats, and
 the capacity must cover a full round or the cache thrashes.

 Eviction can only ever cost a recomputation - values are immutable and keyed exactly - so results
 are unchanged.  Size-triggered eviction is also deterministic in a way a time- or memory-triggered
 policy would not be.
 */
template<class KeyT, class ValueT>
class RecentCache
{
public:
  explicit RecentCache( const std::size_t capacity ) : m_capacity( capacity ) {}

  /** Newest-first lookup; null when absent.  Linear, which for a small capacity beats a tree.
   A zero capacity disables the cache outright, so "is this cache worth its complexity" can be
   measured rather than assumed. */
  std::shared_ptr<const ValueT> find( const KeyT &key ) const
  {
    if( !m_capacity )
      return nullptr;
    for( auto iter = m_entries.rbegin(); iter != m_entries.rend(); ++iter )
      if( iter->first == key )
        return iter->second;
    return nullptr;
  }

  /** Store and return {retained value, true if this call was the one that stored it}.

   A second thread can race the computation of the same key.  Returning the already-stored object
   preserves the previous map-based contract that every racer observes one immutable value.

   Callers that count diagnostics off the bool should note that an evicted-then-recomputed key
   reports true again.  That is the honest reading - the work really was done a second time - but it
   is a change from the old unbounded map, which could only ever report it once. */
  std::pair<std::shared_ptr<const ValueT>,bool> insert( const KeyT &key,
                                                       std::shared_ptr<const ValueT> value )
  {
    if( const std::shared_ptr<const ValueT> existing = find(key) )
      return { existing, false };
    if( !m_capacity )
      return { value, true };   //disabled: hand back the caller's value, retain nothing
    m_entries.emplace_back( key, std::move(value) );
    while( m_entries.size() > m_capacity )
      m_entries.pop_front();
    return { m_entries.back().second, true };
  }

  /** Resize the retained window.  Sized once from the problem's source count before any solve, so
   the shrink path below is defensive rather than exercised in production. */
  void set_capacity( const std::size_t capacity )
  {
    m_capacity = capacity;
    while( m_entries.size() > m_capacity )
      m_entries.pop_front();
  }

  std::size_t size() const { return m_entries.size(); }
  std::size_t capacity() const { return m_capacity; }

  /** Lookup tallies, for measuring whether the retained window is actually big enough.  A window
   smaller than the live working set does not fail, it silently degrades to recomputing the decay
   every time, so hit rate is the only way to tell a working cache from an expensive no-op.  Not
   synchronized: callers hold the same mutex they use for the entries. */
  std::size_t hits() const { return m_hits; }
  std::size_t misses() const { return m_misses; }
  void note_hit() const { ++m_hits; }
  void note_miss() const { ++m_misses; }

private:
  std::size_t m_capacity;
  mutable std::size_t m_hits = 0;
  mutable std::size_t m_misses = 0;
  std::deque<std::pair<KeyT,std::shared_ptr<const ValueT>>> m_entries;
};


/** Retained window for the aging/derivative caches, in entries.

 The useful locality is per source: a DynamicAutoDiff pass asks for a handful of recent ages of a
 given nuclide across many parameter lanes and ROIs.  But the cache key includes the source, so one
 evaluation round queries an entry per nuclide before anything repeats - the window has to hold a
 few ages for *every* source at once, not just a few entries overall.

 At roughly 27 kB per entry, even a 40-source problem bounds each cache near 17 MB, so the sixteen
 cost functors the candidate search keeps alive stay in the tens of megabytes rather than the tens
 of gigabytes an unbounded map reached.

 MEASURED 2026-08-22 on JRC Detective_Pu CBNM_Pu/4h with the shipped Pu presets (all four ages
 fitted), chi2 bit-identical across every arm - this cache only ever costs or saves time:

   preset            cache off     5 ages/src      16 ages/src
   610-775 keV       244.6 s       172.3 s (0.70x) 157.4 s (0.64x)
   120-780 keV      1248.6 s       198.8 s (0.16x) 186.3 s (0.15x)

 So the cache earns its keep decisively - 6.3x on the wide preset - and the original 5 was simply
 undersized: the age-derivative sampler routes its probes through this same cache, so one source
 needs its reference age plus a couple per refinement level at once.  Hit rate saturates near 27%
 (610-775) / 61% (120-780) by 16-64 ages per source, so a larger window buys nothing - which also
 means the unbounded map this replaced was paying its memory for almost no extra hits.
 Sweep with INTERSPEC_REL_ACT_DECAY_CACHE_AGES (0 disables) before changing this. */
inline std::size_t recent_decay_cache_capacity( const std::size_t num_sources )
{
  // `INTERSPEC_REL_ACT_DECAY_CACHE_AGES` overrides the ages-per-source window, so the window can be
  // swept against CPU/hit-rate on a real corpus instead of being argued about.  The shipped value
  // of 5 has never been measured, and the live working set is plausibly larger: the age-derivative
  // sampler routes its probes through this same cache, so a source can need its reference age plus
  // two per refinement level at once.
  static const std::size_t ages_per_source = []() -> std::size_t {
    const char * const value = std::getenv( "INTERSPEC_REL_ACT_DECAY_CACHE_AGES" );
    if( !value )
      return 16;
    const int parsed = std::atoi( value );
    if( parsed == 0 )
      return 0;   //explicitly disable, for the with/without measurement
    return (parsed > 0) ? static_cast<std::size_t>(parsed) : 16;
  }();

  if( !ages_per_source )
    return 0;
  return ages_per_source * (num_sources ? num_sources : 1);
}


enum class AgeDerivativeStencil
{
  Centered,
  Forward,
  Backward
};

enum class AgeDerivativeStatus
{
  Converged,
  BestFiniteEstimate,
  Failed
};

enum class AgeDerivativeFailure
{
  None,
  InvalidInput,
  NoFeasibleStep,
  DuplicateIdentity,
  IdentityMismatch,
  NonFiniteYield,
  NonFiniteDerivative
};

struct AgeDerivativeOptions
{
  /** The first step is this fraction of max(characteristic time, |age|, 1 second). */
  double initial_step_fraction = 0.02;

  /** Richardson compares derivatives after multiplying the absolute tolerance by yield/time. */
  double relative_tolerance = 1.0e-6;
  double absolute_scaled_tolerance = 1.0e-9;

  std::size_t minimum_refinements = 3;
  std::size_t maximum_refinements = 7;
};

struct AgeDerivativeResult
{
  /** Derivatives in the same order as the reference vector passed to the calculator. */
  std::vector<double> derivatives;

  AgeDerivativeStatus status = AgeDerivativeStatus::Failed;
  AgeDerivativeFailure failure = AgeDerivativeFailure::None;
  AgeDerivativeStencil stencil = AgeDerivativeStencil::Centered;
  std::size_t refinements = 0;
  std::size_t unconverged_lines = 0;
  double initial_step = 0.0;
  double final_step = 0.0;
  std::string diagnostic;
};

/** Attach a scalar decay-engine derivative to every active automatic-differentiation lane. */
template<typename JetType>
void attach_age_yield_derivative( JetType &yield, const double derivative,
                                  const JetType &physical_age )
{
  yield.v = derivative * physical_age.v;
}

namespace AgeDerivativeDetail
{

struct IndexedYield
{
  std::uint64_t identity = 0;
  double yield = 0.0;
  std::size_t original_index = 0;
};

inline bool canonical_reference( const std::vector<AgeGammaYield> &input,
                                 std::vector<IndexedYield> &answer,
                                 AgeDerivativeResult &result )
{
  answer.clear();
  answer.reserve( input.size() );
  for( std::size_t i = 0; i < input.size(); ++i )
  {
    if( !std::isfinite(input[i].yield) )
    {
      result.failure = AgeDerivativeFailure::NonFiniteYield;
      result.diagnostic = "The reference gamma-yield vector contains a non-finite value.";
      return false;
    }
    answer.push_back( { input[i].identity, input[i].yield, i } );
  }

  std::sort( answer.begin(), answer.end(), []( const IndexedYield &lhs, const IndexedYield &rhs ) {
    return lhs.identity < rhs.identity;
  } );
  for( std::size_t i = 1; i < answer.size(); ++i )
  {
    if( answer[i - 1].identity == answer[i].identity )
    {
      result.failure = AgeDerivativeFailure::DuplicateIdentity;
      result.diagnostic = "The reference gamma-yield vector contains a duplicate stable identity.";
      return false;
    }
  }
  return true;
}

inline bool align_sample( const std::vector<AgeGammaYield> &sample,
                          const std::vector<IndexedYield> &reference,
                          std::vector<double> &answer,
                          AgeDerivativeResult &result )
{
  if( sample.size() != reference.size() )
  {
    result.failure = AgeDerivativeFailure::IdentityMismatch;
    result.diagnostic = "An aged gamma sample has a different number of stable line identities.";
    return false;
  }

  std::vector<AgeGammaYield> sorted = sample;
  std::sort( sorted.begin(), sorted.end(), []( const AgeGammaYield &lhs, const AgeGammaYield &rhs ) {
    return lhs.identity < rhs.identity;
  } );

  answer.assign( reference.size(), 0.0 );
  for( std::size_t i = 0; i < sorted.size(); ++i )
  {
    if( (i > 0) && (sorted[i - 1].identity == sorted[i].identity) )
    {
      result.failure = AgeDerivativeFailure::DuplicateIdentity;
      result.diagnostic = "An aged gamma sample contains a duplicate stable line identity.";
      return false;
    }
    if( sorted[i].identity != reference[i].identity )
    {
      result.failure = AgeDerivativeFailure::IdentityMismatch;
      result.diagnostic = "An aged gamma sample does not match the reference line identities.";
      return false;
    }
    if( !std::isfinite(sorted[i].yield) )
    {
      result.failure = AgeDerivativeFailure::NonFiniteYield;
      result.diagnostic = "An aged gamma sample contains a non-finite yield.";
      return false;
    }
    answer[reference[i].original_index] = sorted[i].yield;
  }
  return true;
}

} // namespace AgeDerivativeDetail

/** Differentiate every gamma yield with respect to physical source age.

 The same second-order stencil is retained throughout a halving ladder.  Centered and second-order
 one-sided formulas therefore all admit the same Richardson factor of 1/3.  The stencil with the
 largest useful in-bounds step is selected, avoiding tiny centered probes immediately beside a
 physical or user-specified bound.  Returned derivatives are matched by stable gamma identity, not
 by vector position or by the age-dependent dominant transition metadata.

 `sample_at_age` must return finite yields with the same set of identities as `reference`.
 */
template<typename SampleAtAge>
AgeDerivativeResult richardson_age_yield_derivative(
    const double age,
    const double lower_bound,
    const double upper_bound,
    const double characteristic_time,
    const std::vector<AgeGammaYield> &reference,
    SampleAtAge sample_at_age,
    const AgeDerivativeOptions &options = AgeDerivativeOptions{} )
{
  AgeDerivativeResult result;
  result.derivatives.assign( reference.size(), 0.0 );

  const bool valid_options = std::isfinite(options.initial_step_fraction)
                          && (options.initial_step_fraction > 0.0)
                          && std::isfinite(options.relative_tolerance)
                          && (options.relative_tolerance >= 0.0)
                          && std::isfinite(options.absolute_scaled_tolerance)
                          && (options.absolute_scaled_tolerance >= 0.0)
                          && (options.minimum_refinements >= 2)
                          && (options.maximum_refinements >= options.minimum_refinements);
  if( !std::isfinite(age) || !std::isfinite(lower_bound) || !std::isfinite(upper_bound)
      || !std::isfinite(characteristic_time) || (characteristic_time <= 0.0)
      || (lower_bound >= upper_bound) || !valid_options )
  {
    result.failure = AgeDerivativeFailure::InvalidInput;
    result.diagnostic = "Invalid age, bounds, characteristic time, or Richardson options.";
    return result;
  }

  std::vector<AgeDerivativeDetail::IndexedYield> sorted_reference;
  if( !AgeDerivativeDetail::canonical_reference(reference, sorted_reference, result) )
    return result;

  if( reference.empty() )
  {
    result.status = AgeDerivativeStatus::Converged;
    result.diagnostic = "The empty gamma-yield vector has an empty derivative.";
    return result;
  }

  // Scaling a Ceres age parameter to physical seconds can put an otherwise exact bound one ULP
  // outside its input value.  Admit only that roundoff-sized excursion and treat the evaluated age
  // itself as the effective edge; materially infeasible probes remain an explicit failure.
  const double bound_scale = (std::max)(
      1.0, (std::max)(std::fabs(age), (std::max)(std::fabs(lower_bound),std::fabs(upper_bound))) );
  const double bound_roundoff = 128.0*std::numeric_limits<double>::epsilon()*bound_scale;
  if( (age < lower_bound - bound_roundoff) || (age > upper_bound + bound_roundoff) )
  {
    result.failure = AgeDerivativeFailure::InvalidInput;
    result.diagnostic = "The evaluated age lies materially outside its derivative bounds.";
    return result;
  }
  const double effective_lower_bound = (std::min)(lower_bound,age);
  const double effective_upper_bound = (std::max)(upper_bound,age);

  const double time_scale = (std::max)( 1.0, (std::max)(std::fabs(age), characteristic_time) );
  const double nominal_step = options.initial_step_fraction * time_scale;
  const double lower_room = (std::max)(0.0, age - effective_lower_bound);
  const double upper_room = (std::max)(0.0, effective_upper_bound - age);

  // A centered stencil is more accurate for a given h, but directly beside a bound its available
  // h can be orders of magnitude smaller than an inward one-sided step.  Compare their capacities
  // and use centered unless it is less than one quarter of the best one-sided capacity.
  const double centered_capacity = (std::min)(nominal_step, (std::min)(lower_room, upper_room));
  const double forward_capacity = (std::min)(nominal_step, 0.5*upper_room);
  const double backward_capacity = (std::min)(nominal_step, 0.5*lower_room);
  const double best_one_sided = (std::max)(forward_capacity, backward_capacity);

  double step = 0.0;
  if( (centered_capacity > 0.0) && (centered_capacity >= 0.25*best_one_sided) )
  {
    result.stencil = AgeDerivativeStencil::Centered;
    step = centered_capacity;
  }else if( forward_capacity >= backward_capacity )
  {
    result.stencil = AgeDerivativeStencil::Forward;
    step = forward_capacity;
  }else
  {
    result.stencil = AgeDerivativeStencil::Backward;
    step = backward_capacity;
  }

  const double minimum_step = 128.0*std::numeric_limits<double>::epsilon()*time_scale;
  if( !std::isfinite(step) || (step <= minimum_step) )
  {
    result.failure = AgeDerivativeFailure::NoFeasibleStep;
    result.diagnostic = "No numerically resolvable age-difference step fits inside the age bounds.";
    return result;
  }

  result.initial_step = step;

  std::map<std::uint64_t,std::vector<double>> sample_cache;
  const auto sample = [&]( const double sample_age, std::vector<double> &aligned ) -> bool {
    const std::uint64_t key = exact_double_identity(sample_age);
    const auto cached = sample_cache.find(key);
    if( cached != sample_cache.end() )
    {
      aligned = cached->second;
      return true;
    }

    const std::vector<AgeGammaYield> raw = sample_at_age(sample_age);
    if( !AgeDerivativeDetail::align_sample(raw, sorted_reference, aligned, result) )
      return false;
    sample_cache.emplace( key, aligned );
    return true;
  };

  std::vector<double> previous_raw, previous_richardson;
  std::vector<double> best_derivative( reference.size(), 0.0 );
  std::vector<double> best_change( reference.size(), std::numeric_limits<double>::infinity() );
  std::vector<char> has_converged( reference.size(), 0 );
  bool have_finite_estimate = false;

  for( std::size_t level = 0; level < options.maximum_refinements; ++level )
  {
    if( step <= minimum_step )
      break;

    std::vector<double> first, second;
    bool sampled = false;
    switch( result.stencil )
    {
      case AgeDerivativeStencil::Centered:
        sampled = sample(age + step, first) && sample(age - step, second);
        break;
      case AgeDerivativeStencil::Forward:
        sampled = sample(age + step, first) && sample(age + 2.0*step, second);
        break;
      case AgeDerivativeStencil::Backward:
        sampled = sample(age - step, first) && sample(age - 2.0*step, second);
        break;
    }
    if( !sampled )
      return result;

    std::vector<double> raw( reference.size(), 0.0 );
    for( std::size_t i = 0; i < raw.size(); ++i )
    {
      switch( result.stencil )
      {
        case AgeDerivativeStencil::Centered:
          raw[i] = (first[i] - second[i])/(2.0*step);
          break;
        case AgeDerivativeStencil::Forward:
          raw[i] = (-3.0*reference[i].yield + 4.0*first[i] - second[i])/(2.0*step);
          break;
        case AgeDerivativeStencil::Backward:
          raw[i] = (3.0*reference[i].yield - 4.0*first[i] + second[i])/(2.0*step);
          break;
      }
      if( !std::isfinite(raw[i]) )
      {
        result.failure = AgeDerivativeFailure::NonFiniteDerivative;
        result.diagnostic = "An age-difference stencil produced a non-finite derivative.";
        return result;
      }
    }

    ++result.refinements;
    result.final_step = step;

    if( !previous_raw.empty() )
    {
      std::vector<double> richardson( raw.size(), 0.0 );
      for( std::size_t i = 0; i < raw.size(); ++i )
      {
        // Each raw stencil is second-order, so D(h/2) + (D(h/2)-D(h))/(2^2-1).
        richardson[i] = raw[i] + (raw[i] - previous_raw[i])/3.0;
        if( !std::isfinite(richardson[i]) )
        {
          result.failure = AgeDerivativeFailure::NonFiniteDerivative;
          result.diagnostic = "Richardson extrapolation produced a non-finite derivative.";
          return result;
        }
      }

      if( previous_richardson.empty() )
      {
        best_derivative = richardson;
        have_finite_estimate = true;
      }else
      {
        for( std::size_t i = 0; i < richardson.size(); ++i )
        {
          const double change = std::fabs(richardson[i] - previous_richardson[i]);
          const double derivative_scale = (std::max)(std::fabs(richardson[i]),
                                                      std::fabs(previous_richardson[i]));
          const double yield_scale = (std::max)(1.0, std::fabs(reference[i].yield));
          const double tolerance = options.relative_tolerance*derivative_scale
                                 + options.absolute_scaled_tolerance*yield_scale/time_scale;
          if( change < best_change[i] )
          {
            best_change[i] = change;
            best_derivative[i] = richardson[i];
          }
          if( change <= tolerance )
            has_converged[i] = 1;
        }
        have_finite_estimate = true;

        const bool all_converged = std::find(has_converged.begin(), has_converged.end(), char(0))
                                   == has_converged.end();
        if( all_converged && (result.refinements >= options.minimum_refinements) )
        {
          result.derivatives = best_derivative;
          result.status = AgeDerivativeStatus::Converged;
          result.diagnostic = "The bounds-aware Richardson ladder converged for every gamma line.";
          return result;
        }
      }
      previous_richardson.swap(richardson);
    }else
    {
      best_derivative = raw;
      have_finite_estimate = true;
    }

    previous_raw.swap(raw);
    step *= 0.5;
  }

  if( have_finite_estimate )
  {
    result.derivatives = best_derivative;
    result.status = AgeDerivativeStatus::BestFiniteEstimate;
    result.unconverged_lines = static_cast<std::size_t>(
        std::count(has_converged.begin(), has_converged.end(), char(0)) );
    result.diagnostic = "The Richardson ladder did not satisfy its convergence tolerance for "
                      + std::to_string(result.unconverged_lines) + " gamma line(s); the most stable"
                        " finite estimates were retained.";
    return result;
  }

  result.failure = AgeDerivativeFailure::NoFeasibleStep;
  result.diagnostic = "The Richardson ladder could not form a finite derivative estimate.";
  return result;
}

} // namespace RelActCalcAutoImp

#endif // InterSpec_RelActCalcAuto_AgeDerivative_h
