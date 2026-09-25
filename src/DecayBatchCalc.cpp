/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <map>
#include <set>
#include <tuple>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include <boost/tokenizer.hpp>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayBatchCalc.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;

namespace
{
  /** Smallest surviving fraction a nuclide may have and still have its past activity solved for; see
   the use in #back_decay_activities.  2^-47 matches the ~47 half-lives SandiaDecay allows in
   NuclideMixture::addAgedNuclideByNumAtoms, i.e. an amplification of at most ~1.4E+14.
   */
  const double sm_min_back_decay_diagonal = 7.1E-15;

  /** The longest initial age accepted, in half-lives of the nuclide.  An aged seed is scaled up by
   2^(age/half-life) (see NuclideMixture::addAgedNuclideByActivity), which underflows - silently
   giving zero - past ~1000 half-lives; nothing that old could be measured, so leave a wide margin.
   */
  const double sm_max_age_half_lives = 500.0;

  /** Ends the error for an input activity that could not be read. */
  const std::string sm_activity_requirement = "; activities must be non-negative numbers.";

  /** The time of a given step, in seconds.  With `num_steps == 1` returns `time_span`; otherwise
   returns `step*time_span/(num_steps-1)` so the points run 0..time_span inclusive.  Matches the
   legacy batch_decay.cpp behavior.  A negative `time_span` walks backwards from 0 to it.
   */
  double time_for_step( const size_t step, const size_t num_steps, const double time_span )
  {
    if( num_steps <= 1 )
      return time_span;
    return (step * time_span) / static_cast<double>(num_steps - 1);
  }


  /** Formats an activity in the user-selected units, appending the unit label (e.g. "uCi/m2"). */
  string activity_unit_suffix( const bool use_curie, const string &extra_label )
  {
    return string(use_curie ? "Ci" : "Bq") + extra_label;
  }


  /** Formats a value the way the legacy multi-location output did: "%.6E", so zero prints as
   "0.000000E+00".  (SpecUtils::printCompact would give "0", and drops trailing zeros.)
   */
  string print_scientific( const double value )
  {
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.6E", value );
    return string( buffer );
  }


  /** A compact rendering of a time span for the grouped-output "+<time>" label suffix, e.g. "48h".
   `printToBestTimeUnits` gives a fixed-precision value and a space ("48.000000 h"); the legacy
   output has neither, so drop redundant trailing zeros and the space.
   */
  string compact_time_str( const double time_span )
  {
    const string txt = PhysicalUnits::printToBestTimeUnits( time_span, 6 );

    const string::size_type sp = txt.find( ' ' );
    string value = (sp == string::npos) ? txt : txt.substr( 0, sp );
    const string unit = (sp == string::npos) ? string() : txt.substr( sp + 1 );

    if( value.find( '.' ) != string::npos )
    {
      value.erase( value.find_last_not_of( '0' ) + 1 );
      if( !value.empty() && (value.back() == '.') )
        value.pop_back();
    }

    return value + unit;
  }


  /** What #back_decay_activities could not do, accumulated across however many mixtures were solved
   (one per location, so potentially dozens), and formatted once by #back_decay_warnings.  Collecting
   rather than appending text keeps a file of 35 locations from repeating the same warning 35 times.
   */
  struct BackDecayNotes
  {
    /** Nuclides whose own past activity is unidentifiable, but which have an ancestor in the set whose
     in-growth accounts for the measurement; reported as zero.  (With no such ancestor there is nothing
     to report at all, and #back_decay_activities throws instead of noting it here.)
     */
    std::set<const SandiaDecay::Nuclide *> unrecoverable;

    /** Nuclides whose measurement the recovered past state cannot reproduce, with the worst-case
     measured and ancestor-implied present activities seen for each.
     */
    std::map<const SandiaDecay::Nuclide *,std::pair<double,double>> inconsistent;

    /** Inputs whose stated (present) age is shorter than how far back the solve looked, so they were
     taken as freshly made at the past time; see #decay.
     */
    std::set<const SandiaDecay::Nuclide *> age_too_short;

    /** How far back the solve looked, for the message text. */
    double age = 0.0;

    /** The units the messages quote activities in - the divisor from SandiaDecay units, and the text
     to print - so a message never says "Bq" about a table of "uCi/m2".  Set by #decay.
     */
    double act_unit = PhysicalUnits::becquerel;
    std::string unit_str = "Bq";
  };//struct BackDecayNotes


  /** Thrown by #back_decay_activities when a measured nuclide's past activity cannot be determined, and
   no parent solved together with it accounts for the measurement.
   */
  struct UnrecoverablePastError : public std::runtime_error
  {
    using std::runtime_error::runtime_error;
  };


  /** Solves for the activities a set of nuclides must have had `age` ago to give the measured
   activities now.  `age` is a positive magnitude (how far back to look).

   Why this is not just `A_now / exp(-lambda*age)` per nuclide: if any input is a descendant of
   another (Cs137 -> Ba137m, Zr95 -> Nb95, ...) then part of its present activity grew in from that
   ancestor, so the nuclides are coupled and must be solved together.

   The forward map from past amounts to present amounts is linear, and lower-triangular when the
   nuclides are ordered ancestor-before-descendant, so the past amounts follow by forward
   substitution.  Each column of that matrix is obtained from SandiaDecay itself - seed a mixture with
   one unit of nuclide j, evaluate at `age`, and read off what it produced in every i - which avoids
   relying on the "barely tested" branchRatioFromForebear()/branchRatioToDecendant().

   Everything is done in numbers of atoms, which keeps the matrix well conditioned (activity would
   scale rows by wildly different decay constants) and lets stable nuclides participate.

   `seed_ages` are the ages each nuclide will be seeded with at the past time (see #decay), so each
   column carries that nuclide's in-grown progeny exactly as the seeded mixture will; a fresh seed here
   would have the solve miss (and so double count) an aged parent's daughters.

   Rows of one nuclide share a single unknown: in a mixture they add, and nothing tells them apart.  So
   each distinct nuclide is solved once, and its past amount is shared among its rows in proportion to
   their measurements.  A nuclide's surviving fraction does not depend on a row's seed age, so this is
   exact; each row still carries its own column (and so its own in-grown progeny) to the descendants.

   Two things can go wrong, and both are physical rather than numerical:
     - Over enough half-lives a nuclide's surviving fraction `exp(-lambda*age)` becomes too small to
       divide by, or underflows to exactly zero (Ba137m is ~1129 half-lives in 48 h), so its past
       amount has no recoverable effect on anything measurable today.  If an ancestor is among the
       inputs then its in-growth accounts for the measurement, and the nuclide itself is reported as
       zero; if not, there is nothing to report and this throws so the user can look back less far.
     - The measurements may not be mutually consistent with any past state (e.g. a measured
       Ba137m/Cs137 ratio that secular equilibrium forbids).  The solve then cannot reproduce them;
       `warnings` says so rather than silently returning a set that does not decay back to the input.

   Returns past activities in the same (SandiaDecay) units as `activities`, indexed like `nuclides`.
   Throws UnrecoverablePastError when a measured nuclide's past activity is wholly unrecoverable.
   */
  vector<double> back_decay_activities( const vector<const SandiaDecay::Nuclide *> &nuclides,
                                        const vector<double> &activities,
                                        const vector<double> &seed_ages,
                                        const double age,
                                        BackDecayNotes &notes )
  {
    const size_t nnuc = nuclides.size();
    assert( activities.size() == nnuc );
    assert( seed_ages.size() == nnuc );
    assert( age > 0.0 );

    notes.age = age;

    // Order ancestor-before-descendant, so the forward map is lower-triangular.  Counting, for each
    //  nuclide, how many of the other input nuclides are its ancestors gives such an order directly:
    //  if j is an ancestor of i then every ancestor of j is also one of i, plus j itself, so j's count
    //  is strictly smaller than i's.  Sorting by that count ascending therefore always places an
    //  ancestor before its descendants.
    vector<size_t> num_ancestors( nnuc, 0 );
    for( size_t i = 0; i < nnuc; ++i )
    {
      if( !nuclides[i] )
        continue;

      // forebearers() includes the nuclide itself, so skip any input that is the same nuclide.
      const vector<const SandiaDecay::Nuclide *> forebearers = nuclides[i]->forebearers();
      for( size_t j = 0; j < nnuc; ++j )
      {
        if( (j != i) && nuclides[j] && (nuclides[j] != nuclides[i])
           && (std::find( begin(forebearers), end(forebearers), nuclides[j] ) != end(forebearers)) )
        {
          num_ancestors[i] += 1;
        }
      }
    }//for( each nuclide )

    vector<size_t> order( nnuc );
    for( size_t i = 0; i < nnuc; ++i )
      order[i] = i;

    std::stable_sort( begin(order), end(order), [&num_ancestors]( const size_t a, const size_t b ){
      return num_ancestors[a] < num_ancestors[b];
    } );

    // Column j of the forward map: atoms of each nuclide at `age` from one atom of j, seeded the same
    //  way #decay seeds its mixtures (addAgedNuclideByNumAtoms would also refuse ages past ~45
    //  half-lives, which the mixtures take fine).  Also gives the diagonal, j's own surviving
    //  fraction - which an aged seed does not change, since its atom count is the parent's.
    vector<vector<double>> forward( nnuc, vector<double>( nnuc, 0.0 ) );
    for( size_t j = 0; j < nnuc; ++j )
    {
      const SandiaDecay::Nuclide * const src = nuclides[j];
      if( !src )
        continue;

      SandiaDecay::NuclideMixture mix;
      if( (seed_ages[j] > 0.0) && !src->isStable() )
        mix.addAgedNuclideByActivity( src, src->decayConstant(), seed_ages[j] );  // one atom's activity
      else
        mix.addNuclideByAbundance( src, 1.0 );

      const vector<SandiaDecay::NuclideNumAtomsPair> atoms = mix.numAtoms( age );
      for( const SandiaDecay::NuclideNumAtomsPair &nap : atoms )
      {
        for( size_t i = 0; i < nnuc; ++i )
        {
          // Another row of the same nuclide is not something j produces (see the substitution below).
          if( (nuclides[i] == nap.nuclide) && ((i == j) || (nuclides[i] != src)) )
            forward[i][j] = nap.numAtoms;
        }
      }
    }//for( each column )

    // Present-day atoms from the measured activities (a stable nuclide has no decay constant, so it
    //  cannot be expressed as an activity at all - it contributes nothing to solve for).
    vector<double> now_atoms( nnuc, 0.0 );
    for( size_t i = 0; i < nnuc; ++i )
    {
      const SandiaDecay::Nuclide * const nuc = nuclides[i];
      if( nuc && !nuc->isStable() && (nuc->decayConstant() > 0.0) )
        now_atoms[i] = activities[i] / nuc->decayConstant();
    }

#if( PERFORM_DEVELOPER_CHECKS )
    // The whole method rests on `order` making the forward map lower-triangular: a nuclide must not
    //  receive atoms from one solved after it.  Cheap to check, and silently wrong if it ever fails.
    for( size_t a = 0; a < nnuc; ++a )
    {
      for( size_t b = a + 1; b < nnuc; ++b )
        assert( forward[order[a]][order[b]] <= 0.0 );
    }
#endif

    // Forward substitution in ancestor-first order, one distinct nuclide at a time (see above).
    vector<double> past_atoms( nnuc, 0.0 );
    vector<bool> solved( nnuc, false );

    for( size_t pos = 0; pos < nnuc; ++pos )
    {
      const size_t i = order[pos];
      if( solved[i] )
        continue;

      // All rows of this nuclide; they have the same row of the forward map, so `i` stands for them.
      vector<size_t> rows;
      double now_total = 0.0;
      for( size_t k = pos; k < nnuc; ++k )
      {
        if( nuclides[order[k]] == nuclides[i] )
        {
          rows.push_back( order[k] );
          now_total += now_atoms[order[k]];
        }
      }

      // Subtract what the already-solved rows - which by construction include every ancestor of i -
      //  grow into i by now.
      double from_ancestors = 0.0;
      for( size_t k = 0; k < nnuc; ++k )
      {
        if( solved[k] )
          from_ancestors += forward[i][k] * past_atoms[k];
      }

      const double residual = now_total - from_ancestors;
      const double diagonal = forward[i][i];
      double past_total = 0.0;

      // Past enough half-lives a nuclide's own past amount stops being recoverable from today's
      //  measurement: the diagonal underflows to exactly zero (Ba137m over 48 h), or gets so small
      //  that dividing by it amplifies the measurement - and its noise - without bound (at 48 h Pr144
      //  has a diagonal of 7E-51, which would turn a real 1.4E-03 uCi/m2 into 2E+47 uCi/m2).  Both
      //  are the same situation, so treat them the same; SandiaDecay draws this same line at ~47
      //  half-lives in NuclideMixture::addAgedNuclideByNumAtoms.
      if( diagonal < sm_min_back_decay_diagonal )
      {
        // With an ancestor among the inputs, its in-growth accounts for the measurement and zero is
        //  the right answer for i itself - only a note is owed.  With none there is nothing to fall
        //  back on and no meaningful number to report, so have the user look back a shorter way.
        if( (from_ancestors <= 0.0) && (now_total > 0.0) )
        {
          const double half_life = (nuclides[i] && (nuclides[i]->halfLife > 0.0))
                                     ? nuclides[i]->halfLife : 0.0;
          char buffer[512] = { '\0' };
          snprintf( buffer, sizeof(buffer), "Cannot look back %s: that is %.0f half-lives of %s, and"
                    " none of its parents are decayed together with it, so its past activity cannot"
                    " be determined - recovering it would mean scaling the measurement up by more"
                    " than %.0G.  Please use a shorter time.",
                    compact_time_str( age ).c_str(),
                    (half_life > 0.0) ? (age / half_life) : 0.0,
                    (nuclides[i] ? nuclides[i]->symbol.c_str() : "?"),
                    1.0 / sm_min_back_decay_diagonal );
          throw UnrecoverablePastError( buffer );
        }

        if( now_total > 0.0 )
          notes.unrecoverable.insert( nuclides[i] );
      }else if( residual > 0.0 )
      {
        // A negative residual means the ancestors alone already over-produce i; no (non-negative) past
        //  amount of i can fix that, so it stays at zero and the disagreement is reported below.
        past_total = residual / diagonal;
      }

      for( const size_t r : rows )
      {
        past_atoms[r] = (now_total > 0.0) ? (past_total * now_atoms[r] / now_total) : 0.0;
        solved[r] = true;
      }

      // Whenever i's own past amount could not absorb the difference - because it is unidentifiable, or
      //  because the ancestors already over-produce i - the past state cannot reproduce the input.
      const bool absorbed = (diagonal >= sm_min_back_decay_diagonal) && (residual > 0.0);
      if( !absorbed && (now_total > 0.0) && (from_ancestors > 0.0) )
      {
        const double diff = fabs( from_ancestors - now_total ) / now_total;
        if( diff > 0.01 )
        {
          const double lambda = nuclides[i] ? nuclides[i]->decayConstant() : 0.0;

          // Keep the worst disagreement seen for this nuclide across all the mixtures solved.
          std::pair<double,double> &worst = notes.inconsistent[nuclides[i]];
          const double prev = (worst.first > 0.0) ? fabs(worst.second - worst.first)/worst.first : -1.0;
          if( diff > prev )
            worst = std::make_pair( now_total * lambda, from_ancestors * lambda );
        }
      }//if( the measurement could not be matched )
    }//for( each nuclide, ancestors first )

    // Back to activities.
    vector<double> past_activities( nnuc, 0.0 );
    for( size_t i = 0; i < nnuc; ++i )
    {
      const SandiaDecay::Nuclide * const nuc = nuclides[i];
      if( nuc && !nuc->isStable() )
        past_activities[i] = past_atoms[i] * nuc->decayConstant();
    }

    return past_activities;
  }//back_decay_activities(...)


  /** How far apart two activities are, as text: a percentage while that stays readable, and a plain
   multiple once it does not (a nearly-stable nuclide's ancestors can imply many orders of magnitude
   more than was measured, and "+2.5E+50%" tells the reader nothing).
   */
  string difference_str( const double measured, const double implied )
  {
    char buffer[64] = { '\0' };

    if( (measured > 0.0) && (implied > (2.0 * measured)) )
      snprintf( buffer, sizeof(buffer), "%.3G times more", implied / measured );
    else if( measured > 0.0 )
      snprintf( buffer, sizeof(buffer), "%+.1f%%", 100.0 * (implied - measured) / measured );
    else
      snprintf( buffer, sizeof(buffer), "measured as zero" );

    return string( buffer );
  }//difference_str(...)


  /** Renders the notes accumulated over every back-decay solve into user-facing warning text; one
   message per affected nuclide no matter how many locations hit it, and only the worst few spelled
   out - a whole fission-product mix can be inconsistent in dozens of nuclides at once.
   */
  string back_decay_warnings( const BackDecayNotes &notes )
  {
    const size_t max_detailed = 5;
    string warnings;

    // Nothing to report unless a solve ran, and a solve always sets `age` - so a note without one
    //  would render as "over 0s it decays away entirely".
    assert( (notes.age > 0.0)
           || (notes.unrecoverable.empty() && notes.inconsistent.empty()
               && notes.age_too_short.empty()) );

    if( !notes.age_too_short.empty() )
    {
      string names;
      for( const SandiaDecay::Nuclide * const nuc : notes.age_too_short )
        names += (names.empty() ? "" : ", ") + (nuc ? nuc->symbol : string("?"));

      warnings += "The initial age of " + names + " is shorter than " + compact_time_str( notes.age )
                  + ", so it did not exist that long ago; it is taken as freshly made "
                  + compact_time_str( notes.age ) + " ago instead.\n";
    }

    if( !notes.unrecoverable.empty() )
    {
      string names;
      for( const SandiaDecay::Nuclide * const nuc : notes.unrecoverable )
        names += (names.empty() ? "" : ", ") + (nuc ? nuc->symbol : string("?"));

      warnings += "Past activity of " + names + " cannot be determined: over "
                  + compact_time_str( notes.age ) + " it decays away entirely, so today's measurement"
                    " reflects only in-growth from its parent.  Reported as zero.\n";
    }

    if( notes.inconsistent.empty() )
      return warnings;

    // Worst disagreement first, so the detailed lines are the ones worth reading.
    vector<const SandiaDecay::Nuclide *> worst;
    for( const std::pair<const SandiaDecay::Nuclide * const,std::pair<double,double>> &nv
        : notes.inconsistent )
    {
      worst.push_back( nv.first );
    }

    std::sort( begin(worst), end(worst), [&notes]( const SandiaDecay::Nuclide *l,
                                                  const SandiaDecay::Nuclide *r ){
      const std::pair<double,double> &a = notes.inconsistent.at( l );
      const std::pair<double,double> &b = notes.inconsistent.at( r );
      const double la = (a.first > 0.0) ? fabs(a.second - a.first)/a.first : 0.0;
      const double lb = (b.first > 0.0) ? fabs(b.second - b.first)/b.first : 0.0;
      return la > lb;
    } );

    warnings += "Input activities are not self-consistent with having decayed from a state "
                + compact_time_str( notes.age ) + " ago, for " + std::to_string( worst.size() )
                + (worst.size() == 1 ? string(" nuclide") : string(" nuclides"))
                + ".  The reported past activities will not decay forward to exactly the input"
                  " values.\n";

    for( size_t i = 0; (i < worst.size()) && (i < max_detailed); ++i )
    {
      const SandiaDecay::Nuclide * const nuc = worst[i];
      const std::pair<double,double> &mi = notes.inconsistent.at( nuc );

      char buffer[512] = { '\0' };
      snprintf( buffer, sizeof(buffer), "  %s is measured at %.6G %s, but its parents imply"
                " %.6G %s (%s).\n",
                (nuc ? nuc->symbol.c_str() : "?"),
                mi.first / notes.act_unit, notes.unit_str.c_str(),
                mi.second / notes.act_unit, notes.unit_str.c_str(),
                difference_str( mi.first, mi.second ).c_str() );
      warnings += buffer;
    }//for( the worst few )

    if( worst.size() > max_detailed )
    {
      warnings += "  ...and " + std::to_string( worst.size() - max_detailed )
                  + " others.\n";
    }

    return warnings;
  }//back_decay_warnings(...)
}//anonymous namespace


namespace DecayBatchCalc
{

/** The grouped/long-format half of #decay, for input that carries locations.

 Each location is co-decayed as its own mixture and gets its own block of output rows; nothing is ever
 summed across locations.  The row set is the union of every location's nuclides and progeny, so all
 locations list the same nuclides in the same order (even those that are zero there) - which matches
 the legacy output this format reproduces, and makes the blocks directly comparable.

 Unlike the wide format, values are written in the input's own units and in "%.6E", and stable
 nuclides are listed too (as an activity, hence always zero).

 `eval_time` maps a step index to the mixture time to evaluate, `seed_age` gives the age to seed an
 input with, and `past_activities` back-decays one location's inputs when looking backwards; all come
 from #decay.
 */
static void decay_grouped( const vector<BatchNuclide> &inputs,
                           const BatchDecayOptions &opts,
                           const std::function<double(size_t)> &eval_time,
                           const std::function<double(const BatchNuclide &)> &seed_age,
                           const std::function<vector<double>(const vector<BatchNuclide> &)> &past_activities,
                           BatchDecayResult &result )
{
  // Group inputs by location, keeping first-seen order (so output follows the input file's order).
  vector<string> locations;
  std::map<string,vector<BatchNuclide>> by_location;
  for( const BatchNuclide &in : inputs )
  {
    if( by_location.find(in.location) == end(by_location) )
      locations.push_back( in.location );
    by_location[in.location].push_back( in );
  }

  // Header: mirror the source columns.  Rows take their extra cells from their own location's first
  //  input, so the header must name the columns of whichever input has the most of them - a row that
  //  was short a trailing cell would otherwise shift every later column under the wrong header.
  const BatchNuclide *header_src = &inputs.front();
  for( const BatchNuclide &in : inputs )
  {
    if( in.extra_columns.size() > header_src->extra_columns.size() )
      header_src = &in;
  }

  // Prefer the source file's own names (see BatchNuclide::fixed_column_names); the literals are the
  //  legacy names, used only for input that did not carry a header.
  const vector<string> &src_names = header_src->fixed_column_names;
  auto fixed_name = [&src_names]( const size_t i, const char * const fallback ) -> string {
    return ((i < src_names.size()) && !src_names[i].empty()) ? src_names[i] : string( fallback );
  };

  result.column_headers.push_back( fixed_name( 0, "Probe name" ) );
  result.column_headers.push_back( fixed_name( 1, "Product" ) );
  for( const std::pair<string,string> &nv : header_src->extra_columns )
    result.column_headers.push_back( nv.first );
  result.column_headers.push_back( fixed_name( 2, "Value" ) );
  result.column_headers.push_back( fixed_name( 3, "Unit" ) );

  const size_t num_extra_cols = header_src->extra_columns.size();

  // The nuclide list, unioned over all locations so every block lists the same rows.  A mixture's
  //  decayedToNuclidesEvolutions() is already sorted by (mass number, atomic number, isomer), and a
  //  union of such lists needs that same ordering - which is what Nuclide::lessThanForOrdering gives.
  vector<const SandiaDecay::Nuclide *> all_nuclides;
  for( const BatchNuclide &in : inputs )
  {
    for( const SandiaDecay::Nuclide * const nuc : in.nuclide->descendants() )
    {
      if( std::find( begin(all_nuclides), end(all_nuclides), nuc ) == end(all_nuclides) )
        all_nuclides.push_back( nuc );
    }
  }

  std::sort( begin(all_nuclides), end(all_nuclides),
            []( const SandiaDecay::Nuclide *l, const SandiaDecay::Nuclide *r ){
    return SandiaDecay::Nuclide::lessThanForOrdering( l, r );
  } );

  // Each step's "+48h" label for the Product column; the same for every location, so built once.
  vector<string> time_suffixes( opts.num_steps );
  for( size_t step = 0; step < opts.num_steps; ++step )
  {
    // A single-step run echoes the user's own notation (giving the legacy "+48h" for a typed "48h");
    //  multiple steps each need their own time, so they are formatted compactly.
    string suffix = ((opts.num_steps == 1) && !opts.time_span_str.empty())
                      ? opts.time_span_str
                      : compact_time_str( time_for_step( step, opts.num_steps, opts.time_span ) );
    SpecUtils::trim( suffix );

    // An all-whitespace user string trims to nothing, which would leave a bare "+".
    if( suffix.empty() )
      suffix = compact_time_str( time_for_step( step, opts.num_steps, opts.time_span ) );

    // The suffix reads as an offset ("Cs137 ... +48h"), so supply the sign unless the text already
    //  carries one - a user may type either "48h" or "+48h", and a negative time signs itself.
    if( (suffix.front() != '-') && (suffix.front() != '+') )
      suffix = "+" + suffix;

    time_suffixes[step] = suffix;
  }//for( each step )

  vector<string> locations_mixed_units;

  for( const string &location : locations )
  {
    const vector<BatchNuclide> &group = by_location[location];

    // One mixture per location; when looking backwards the location's nuclides are solved together,
    //  since an ancestor among them feeds its descendants.
    const vector<double> acts = past_activities( group );

    SandiaDecay::NuclideMixture mix;
    for( size_t i = 0; i < group.size(); ++i )
    {
      // A zero (or negative) activity must not be seeded: addAgedNuclideByActivity() forms
      //  `activity/aged_activity`, which is 0/0 = NaN for a zero input, and NuclideMixture then
      //  *sums* the seeded atoms of a repeated nuclide - so a zero row that is an ancestor of a
      //  non-zero row would turn that sibling's real activity into NaN, which
      //  NuclideTimeEvolution::numAtoms() silently clamps to zero via max(0.0,NaN).  Leaving the
      //  nuclide out of the mixture instead reports it as zero (see the `in_mix` test below),
      //  which is the right answer and keeps every location listing the same rows.
      if( acts[i] > 0.0 )
        mix.addAgedNuclideByActivity( group[i].nuclide, acts[i], seed_age( group[i] ) );
    }

    // What this mixture can actually be asked about (see the `in_mix` check below).
    vector<const SandiaDecay::Nuclide *> mix_nuclides;
    for( int i = 0; i < mix.numSolutionNuclides(); ++i )
      mix_nuclides.push_back( mix.solutionNuclide(i) );

    // Progeny rows have no input row of their own, so they inherit the location's passthrough cells.
    const BatchNuclide &first = group.front();

    // A row is written in its *own* input row's units; a progeny row has no input of its own, so it
    //  takes the location's first.  Using `first`'s unit for everything would silently scale a row
    //  whose input used a different one (e.g. mCi/m2 among uCi/m2 rows) while labelling it correctly.
    std::map<const SandiaDecay::Nuclide *,const BatchNuclide *> input_of;
    bool mixed_units = false;
    for( const BatchNuclide &in : group )
    {
      input_of[in.nuclide] = &in;
      mixed_units |= ((in.activity_unit != first.activity_unit)
                      || (in.unit_label != first.unit_label));
    }

    if( mixed_units )
      locations_mixed_units.push_back( location );

    // Each step is its own block of rows, labeled with the time it is at.
    for( size_t step = 0; step < opts.num_steps; ++step )
    {
      const string &time_suffix = time_suffixes[step];
      const double t = eval_time( step );

      for( const SandiaDecay::Nuclide * const nuc : all_nuclides )
      {
        const std::map<const SandiaDecay::Nuclide *,const BatchNuclide *>::const_iterator pos
                                                                              = input_of.find( nuc );
        const BatchNuclide &src = (pos == end(input_of)) ? first : *(pos->second);
        const double act_unit = src.activity_unit.empty()
                                 ? PhysicalUnits::becquerel
                                 : PhysicalUnits::stringToActivity( "1" + src.activity_unit );

        // The row set is the union over *every* location, so a nuclide may be absent from this
        //  location's mixture (it has none of it, and no ancestor of it) - which NuclideMixture
        //  reports by throwing rather than returning zero.  Zero is exactly the right answer, and
        //  keeps every block listing the same rows.
        const bool in_mix = (std::find( begin(mix_nuclides), end(mix_nuclides), nuc ) != end(mix_nuclides));
        const double act = in_mix ? (mix.activity( t, nuc ) / act_unit) : 0.0;
        assert( !IsNan(act) );
        if( (opts.min_activity > 0.0) && (!(act >= opts.min_activity)) )
          continue;

        vector<string> row;
        row.push_back( location );
        row.push_back( nuc->symbol + first.product_suffix + time_suffix );

        // Pad (or truncate) to the header's width, so a source row that was missing a trailing cell
        //  cannot shift Value/Unit under the wrong heading.
        for( size_t c = 0; c < num_extra_cols; ++c )
          row.push_back( (c < first.extra_columns.size()) ? first.extra_columns[c].second : string() );

        row.push_back( print_scientific( act ) );
        row.push_back( src.activity_unit + src.unit_label );
        assert( row.size() == result.column_headers.size() );

        result.rows.push_back( std::move(row) );
      }//for( each nuclide )
    }//for( each step )
  }//for( each location )

  if( !locations_mixed_units.empty() )
  {
    // Each row is still correct in its own stated unit, but the progeny rows had to pick one, so say so.
    string names;
    for( size_t i = 0; (i < locations_mixed_units.size()) && (i < 5); ++i )
      names += (names.empty() ? "" : ", ") + locations_mixed_units[i];
    if( locations_mixed_units.size() > 5 )
      names += ", ...";

    result.warnings += "Activity units differ between rows of " + names + "; each row is written in"
                       " its own unit, and progeny rows use the first unit of their location.\n";
  }

  // Long format has one value per row, so the GUI preview renders rows*columns cells.
  result.num_data_cells = result.rows.size() * result.column_headers.size();
}//decay_grouped(...)


string location_key( const string &probe_name )
{
  const string::size_type us = probe_name.rfind( '_' );
  if( (us == string::npos) || (us == 0) || ((us + 1) == probe_name.size()) )
    return probe_name;

  const string head = probe_name.substr( 0, us );
  const string tail = probe_name.substr( us + 1 );

  // A purely numeric suffix is a row index (e.g. "TeamA-01_39", "Grid_7").
  if( tail.find_first_not_of( "0123456789" ) == string::npos )
    return head;

  // Otherwise only strip when the name looks like a structured ID *and* the suffix looks like an index
  //  rather than a place name: the ID test alone would merge "Site-A_North" with "Site-A_South", and a
  //  location that a name distinguishes must never be merged with another.  So require the suffix to
  //  contain a digit (e.g. "TeamA-01_run3", "TeamA-01_r3") - a purely alphabetic suffix like "North"
  //  is taken to name the location.
  if( tail.find_first_of( "0123456789" ) == string::npos )
    return probe_name;

  // The '-' must be in the last segment of `head`, so the ID is what the suffix indexes (this keeps
  //  "Site-A_North2"-style names from being stripped on the strength of an unrelated earlier hyphen).
  const string::size_type seg = head.find_last_of( '_' );
  const string::size_type from = (seg == string::npos) ? 0 : (seg + 1);
  for( string::size_type i = from + 1; (i + 1) < head.size(); ++i )
  {
    if( (head[i] == '-') && isalnum( static_cast<unsigned char>(head[i-1]) )
       && isalnum( static_cast<unsigned char>(head[i+1]) ) )
    {
      return head;
    }
  }

  return probe_name;
}//location_key(...)


BatchDecayResult decay( const vector<BatchNuclide> &inputs, const BatchDecayOptions &opts )
{
  BatchDecayResult result;

  if( opts.time_span == 0.0 )
    throw runtime_error( "A non-zero decay time must be given." );

  if( opts.num_steps < 1 )
    throw runtime_error( "At least one time step is required." );

  if( inputs.empty() )
    throw runtime_error( "No input nuclides were provided." );

  const bool show_progeny = (opts.show_progeny || opts.mix_input);
  const double act_unit = opts.use_curie ? PhysicalUnits::curie : PhysicalUnits::becquerel;

  // A negative time span means "what were the activities |time_span| ago".  We solve for that past
  //  state once (per mixture - see back_decay_activities), seed the mixtures with it, and then decay
  //  *forwards* from there, so every step is a real point on one trajectory: the last step is the past
  //  state itself, and step 0 lands back on (approximately) the measured input.
  const bool back_decay = (opts.time_span < 0.0);
  const double back_age = -opts.time_span;

  // The forward time at which to evaluate a mixture for a given step.
  auto eval_time = [&opts,back_decay,back_age]( const size_t step ) -> double {
    const double t = time_for_step( step, opts.num_steps, opts.time_span );
    return back_decay ? (back_age + t) : t;
  };//eval_time

  // The age to seed an input's mixture with.  An input's age is its age at the measurement, i.e. at
  //  time zero of the table; looking backwards the mixture starts `back_age` earlier, when the sample
  //  was that much younger.  An age shorter than that (including the zero of "no age given") means
  //  the sample did not exist yet, so it is taken as freshly made at the past time.
  auto seed_age = [back_decay,back_age]( const BatchNuclide &in ) -> double {
    return back_decay ? std::max( 0.0, in.age - back_age ) : in.age;
  };//seed_age

  // Replaces each input's activity with the activity it must have had `back_age` ago, for one set of
  //  nuclides that share a mixture (so ancestor/descendant coupling is accounted for).  A no-op when
  //  not back-decaying.  Whatever could not be recovered accumulates into `notes`, which is turned into
  //  warning text once every mixture has been solved (there is one per location).
  BackDecayNotes notes;
  notes.act_unit = act_unit;
  notes.unit_str = activity_unit_suffix( opts.use_curie, string() );
  if( back_decay )
    notes.age = back_age;   // also set by the solve, but the warning text must not depend on that
  auto past_activities = [back_decay,back_age,&notes,&seed_age]( const vector<BatchNuclide> &group ) -> vector<double> {
    vector<double> acts( group.size(), 0.0 );
    for( size_t i = 0; i < group.size(); ++i )
      acts[i] = group[i].activity;

    if( !back_decay )
      return acts;

    vector<const SandiaDecay::Nuclide *> nucs( group.size(), nullptr );
    vector<double> ages( group.size(), 0.0 );
    for( size_t i = 0; i < group.size(); ++i )
    {
      nucs[i] = group[i].nuclide;
      ages[i] = seed_age( group[i] );
      if( (group[i].age > 0.0) && (group[i].age < back_age) )
        notes.age_too_short.insert( group[i].nuclide );
    }

    return back_decay_activities( nucs, acts, ages, back_age, notes );
  };//past_activities

  // The inputs we can actually decay.
  vector<BatchNuclide> valid_inputs;
  for( const BatchNuclide &in : inputs )
  {
    if( !in.nuclide )
      result.warnings += "Skipped invalid nuclide '" + in.nuclide_str + "'.\n";
    else if( in.nuclide->isStable() )
      result.warnings += "Skipped stable nuclide '" + in.nuclide_str + "'.\n";
    else if( IsNan(in.activity) || IsInf(in.activity) || (in.activity < 0.0) )
      throw runtime_error( "The activity of " + in.nuclide->symbol + " is not valid" + sm_activity_requirement );
    else if( !(in.age >= 0.0) || ((in.age / in.nuclide->halfLife) > sm_max_age_half_lives) )
    {
      char buffer[256] = { '\0' };
      snprintf( buffer, sizeof(buffer), "The initial age of %s, %.4G half-lives, cannot be used: ages"
                " from zero up to %.0f half-lives are supported.", in.nuclide->symbol.c_str(),
                in.age / in.nuclide->halfLife, sm_max_age_half_lives );
      throw runtime_error( buffer );
    }else
      valid_inputs.push_back( in );
  }//for( each input )

  if( valid_inputs.empty() )
    throw runtime_error( "No valid, unstable nuclides to decay." );

  // Grouped (multiple physical locations) input gets its own, quite different, output shape.
  bool any_location = false;
  for( const BatchNuclide &in : valid_inputs )
    any_location |= !in.location.empty();

  if( any_location )
  {
    // Grouped output stays in the input's own unit (`use_curie` does not apply), so the warnings must
    //  quote that unit too.
    const BatchNuclide &first_in = valid_inputs.front();
    if( !first_in.activity_unit.empty() )
    {
      notes.act_unit = PhysicalUnits::stringToActivity( "1" + first_in.activity_unit );
      notes.unit_str = first_in.activity_unit + first_in.unit_label;
    }else
    {
      // Grouped values are written as-given (parse_csv leaves `activity_unit` empty for a bare number,
      //  or for a unit cell of just "/m2"), so quoting "Bq"/"Ci" from `use_curie` would be wrong.
      notes.act_unit = PhysicalUnits::becquerel;
      notes.unit_str = first_in.unit_label;
    }

    decay_grouped( valid_inputs, opts, eval_time, seed_age, past_activities, result );
    result.warnings += back_decay_warnings( notes );  // after every location has been solved
    return result;
  }

  // Build the header row.
  result.column_headers.push_back( "Nuclide" );
  for( size_t step = 0; step < opts.num_steps; ++step )
  {
    const double t = time_for_step( step, opts.num_steps, opts.time_span );
    result.column_headers.push_back( PhysicalUnits::printToBestTimeUnits( t, 3 ) );
  }

  // A "source" is one decaying mixture: the single combined mixture when mixing inputs, or one
  //  per input otherwise.  `parent` is the input nuclide for the un-mixed case (used to label and to
  //  identify the initial nuclide); it is null when mixing (the mixture has several initial nuclides).
  struct Source
  {
    std::unique_ptr<SandiaDecay::NuclideMixture> mix;
    const SandiaDecay::Nuclide *parent = nullptr;
    string unit_label;
  };//struct Source

  vector<Source> sources;

  if( opts.mix_input )
  {
    // Sum all valid inputs into a single mixture and co-decay.  When looking backwards the whole set
    //  is solved together, since they share one mixture and so may feed each other.
    const vector<double> acts = past_activities( valid_inputs );

    Source src;
    src.mix.reset( new SandiaDecay::NuclideMixture() );
    bool consistent_unit = true;
    for( size_t i = 0; i < valid_inputs.size(); ++i )
    {
      const BatchNuclide &in = valid_inputs[i];
      if( acts[i] > 0.0 )                    // see the note in decay_grouped(): 0 activity -> NaN
        src.mix->addAgedNuclideByActivity( in.nuclide, acts[i], seed_age( in ) );
      if( i == 0 )
        src.unit_label = in.unit_label;      // capture the first added input's label
      else if( in.unit_label != src.unit_label )
        consistent_unit = false;             // subsequent labels must match
    }//for( each input )

    if( !src.mix->numInitialNuclides() )
      throw runtime_error( "All input activities are zero; there is nothing to decay." );

    if( !consistent_unit )
    {
      result.warnings += "Mixed inputs had differing unit labels; dropped from output.\n";
      src.unit_label.clear();
    }

    sources.push_back( std::move(src) );
  }else
  {
    // Each input is decayed independently, so looking backwards each is solved on its own - which
    //  reduces to A_past = A_now/exp(-lambda*|t|), matching DecayActivityDiv.
    for( const BatchNuclide &in : valid_inputs )
    {
      vector<double> acts;
      try
      {
        acts = past_activities( vector<BatchNuclide>{ in } );
      }catch( UnrecoverablePastError &e )
      {
        // Solved on its own, a parent among the other inputs cannot account for it - but mixing would.
        const vector<const SandiaDecay::Nuclide *> forebearers = in.nuclide->forebearers();
        for( const BatchNuclide &other : valid_inputs )
        {
          if( (other.nuclide != in.nuclide)
             && (std::find( begin(forebearers), end(forebearers), other.nuclide ) != end(forebearers)) )
          {
            throw runtime_error( string( e.what() ) + "  Its parent " + other.nuclide->symbol
                                 + " is among the inputs; check \"Mix inputs\" to decay them together." );
          }
        }//for( each other input )

        throw;
      }//try / catch

      // A zero activity would seed NaN (see decay_grouped()), and on its own it has no progeny to
      //  report either, so the input simply contributes no rows.
      if( acts[0] <= 0.0 )
      {
        result.warnings += "Skipped '" + in.nuclide->symbol + "': activity is zero.\n";
        continue;
      }

      Source src;
      src.mix.reset( new SandiaDecay::NuclideMixture() );
      src.parent = in.nuclide;
      src.unit_label = in.unit_label;
      src.mix->addAgedNuclideByActivity( in.nuclide, acts[0], seed_age( in ) );
      sources.push_back( std::move(src) );
    }//for( each input )
  }//if( mix_input ) / else

  if( sources.empty() )
    throw runtime_error( "No valid, unstable nuclides to decay." );

  // Whether a nuclide is one of a source's initial (input) nuclides.
  auto is_initial_nuclide = [&]( const Source &src, const SandiaDecay::Nuclide *nuc ) -> bool {
    if( src.parent )
      return (nuc == src.parent);
    for( int i = 0; i < src.mix->numInitialNuclides(); ++i )
    {
      if( src.mix->initialNuclide(i) == nuc )
        return true;
    }
    return false;
  };//is_initial_nuclide

  // --- Section 1: activities (initial nuclides, plus progeny when requested), grouped together. ---
  if( opts.include_activity )
  {
    for( const Source &src : sources )
    {
      const vector<SandiaDecay::NuclideTimeEvolution> &evos = src.mix->decayedToNuclidesEvolutions();
      for( const SandiaDecay::NuclideTimeEvolution &evo : evos )
      {
        const SandiaDecay::Nuclide * const nuc = evo.nuclide;
        if( !nuc || IsInf(nuc->halfLife) )  // skip stable nuclides
          continue;

        const bool is_initial = is_initial_nuclide( src, nuc );
        if( !show_progeny && !is_initial )
          continue;

        // "I135 activity (Ci)"; for an un-mixed progeny, note its originating input to disambiguate
        //  shared progeny that are reported separately (e.g. "Xe135 activity (Ci) (from I135)").
        string label = nuc->symbol + " activity ("
                       + activity_unit_suffix(opts.use_curie, src.unit_label) + ")";
        if( !is_initial && src.parent )
          label += " (from " + src.parent->symbol + ")";

        vector<string> row;
        row.push_back( std::move(label) );
        bool any_above_cut = false;
        for( size_t step = 0; step < opts.num_steps; ++step )
        {
          const double act = src.mix->activity( eval_time(step), nuc ) / act_unit;
          assert( !IsNan(act) );
          any_above_cut |= (act >= opts.min_activity);
          row.push_back( SpecUtils::printCompact( act, 6 ) );
        }

        // A row spans several time steps, so only drop it when no step makes the cut (which keeps the
        //  table rectangular).
        if( (opts.min_activity > 0.0) && !any_above_cut )
          continue;

        result.rows.push_back( std::move(row) );
      }//for( each solution nuclide )
    }//for( each source )
  }//if( include_activity )

  // A single emission line, attributed to the (parent -> child) transition that produced it.
  struct AttribLine
  {
    double energy;
    const SandiaDecay::Nuclide *parent;
    const SandiaDecay::Nuclide *child;
    double rate;
  };//struct AttribLine

  // Enumerates the `want` emission lines of `mix` at time `t`, keeping the (parent -> child)
  //  transition context that NuclideMixture::decayParticle discards.  Reproduces that function's rate
  //  (activity * intensity * branchRatio) and its vetoes: pure-E0 gammas and per-nuclide duplicate
  //  x-rays (Br83-style GS/isomer repeats).
  auto enumerate_lines = [&]( const SandiaDecay::NuclideMixture &mix, const double t,
                              const SandiaDecay::ProductType want ) -> vector<AttribLine> {
    vector<AttribLine> out;
    const vector<SandiaDecay::NuclideActivityPair> acts = mix.activity( t );
    for( const SandiaDecay::NuclideActivityPair &ap : acts )
    {
      const SandiaDecay::Nuclide * const nuc = ap.nuclide;
      if( !nuc )
        continue;
      const double activity = ap.activity;
      const size_t ndecays = nuc->decaysToChildren.size();
      for( size_t di = 0; di < ndecays; ++di )
      {
        const SandiaDecay::Transition * const tr = nuc->decaysToChildren[di];
        if( !tr )
          continue;

        // Reproduce decayParticle's x-ray de-duplication across a nuclide's transitions.
        if( want == SandiaDecay::XrayParticle )
        {
          bool duplicate = false;
          for( size_t pj = 0; !duplicate && (pj < di); ++pj )
          {
            const SandiaDecay::Transition * const pt = nuc->decaysToChildren[pj];
            duplicate = ( tr->parent && pt->parent && tr->child && pt->child
                          && (tr->parent->atomicNumber == pt->parent->atomicNumber)
                          && (tr->parent->massNumber == pt->parent->massNumber)
                          && (tr->parent->isomerNumber == pt->parent->isomerNumber)
                          && (tr->child->atomicNumber == pt->child->atomicNumber) );
          }
          if( duplicate )
            continue;
        }//if( want == XrayParticle )

        for( const SandiaDecay::RadParticle &p : tr->products )
        {
          if( p.type != want )
            continue;
          if( (want == SandiaDecay::GammaParticle) && p.e0_verified )  // pure E0 emits no photon
            continue;
          out.push_back( AttribLine{ p.energy, tr->parent, tr->child,
                                     activity * p.intensity * tr->branchRatio } );
        }//for( each product )
      }//for( each transition )
    }//for( each nuclide )
    return out;
  };//enumerate_lines

  // Appends one section of particle-line rows (one line per distinct energy+transition), across all
  //  sources.  Energies are keyed to a millieV grid; the set of (energy,transition) keys is unioned
  //  across time steps so an in-growth line absent at one step is still reported.  For gammas,
  //  positron transitions additionally contribute the 511 keV annihilation line (2 photons each),
  //  matching NuclideMixture::gammas(..., includeAnnihilation=true).
  auto add_particle_section = [&]( const char *type_word, const SandiaDecay::ProductType want,
                                   const bool annihilation ){
    using Key = std::tuple<long long, const SandiaDecay::Nuclide *, const SandiaDecay::Nuclide *>;
    for( const Source &src : sources )
    {
      vector<std::map<Key,double>> step_maps( opts.num_steps );
      std::map<Key,double> key_energy;  // representative energy per key (ordered by energy, then ptr)

      auto accumulate = [&]( const size_t step, const double energy,
                             const SandiaDecay::Nuclide *par, const SandiaDecay::Nuclide *chld,
                             const double rate ){
        const Key key{ std::llround( energy * 1.0E6 ), par, chld };
        step_maps[step][key] += rate;
        key_energy.emplace( key, energy );
      };//accumulate

      for( size_t step = 0; step < opts.num_steps; ++step )
      {
        const double t = eval_time( step );
        for( const AttribLine &ln : enumerate_lines( *src.mix, t, want ) )
          accumulate( step, ln.energy, ln.parent, ln.child, ln.rate );

        if( annihilation )
        {
          for( const AttribLine &ln : enumerate_lines( *src.mix, t, SandiaDecay::PositronParticle ) )
            accumulate( step, 510.998910, ln.parent, ln.child, 2.0 * ln.rate );
        }
      }//for( each step )

      for( const std::pair<const Key,double> &ke : key_energy )
      {
        const SandiaDecay::Nuclide * const par = std::get<1>( ke.first );
        const SandiaDecay::Nuclide * const chld = std::get<2>( ke.first );

        stringstream lbl;
        lbl << SpecUtils::printCompact( ke.second, 6 ) << " keV " << type_word << "/s ("
            << (par ? par->symbol : string("?"));
        if( chld )
          lbl << " -> " << chld->symbol;
        lbl << ")";

        vector<string> row;
        row.push_back( lbl.str() );
        for( size_t step = 0; step < opts.num_steps; ++step )
        {
          const std::map<Key,double>::const_iterator it = step_maps[step].find( ke.first );
          row.push_back( SpecUtils::printCompact( (it != step_maps[step].end()) ? it->second : 0.0, 6 ) );
        }
        result.rows.push_back( std::move(row) );
      }//for( each distinct energy+transition )
    }//for( each source )
  };//add_particle_section

  // --- Section 2+: particle lines, each type grouped together (gammas, x-rays, alphas, betas). ---
  if( opts.include_gammas )
    add_particle_section( "gamma", SandiaDecay::GammaParticle, true );
  if( opts.include_xrays )
    add_particle_section( "xray",  SandiaDecay::XrayParticle,  false );
  if( opts.include_alphas )
    add_particle_section( "alpha", SandiaDecay::AlphaParticle, false );
  if( opts.include_betas )
    add_particle_section( "beta",  SandiaDecay::BetaParticle,  false );

  // Cells the GUI would render - the label column included, so this means the same thing here as in
  //  the grouped path (which has several label columns).
  result.num_data_cells = result.rows.size() * result.column_headers.size();
  result.warnings += back_decay_warnings( notes );  // after every source has been solved

  return result;
}//decay(...)


string result_to_csv( const BatchDecayResult &result )
{
  auto escape = []( const string &field ) -> string {
    if( field.find_first_of( ",\"\r\n" ) == string::npos )
      return field;
    string out = "\"";
    for( const char c : field )
    {
      if( c == '"' )
        out += "\"\"";
      else
        out += c;
    }
    out += "\"";
    return out;
  };//escape

  stringstream out;
  const string eol = "\r\n";

  for( size_t i = 0; i < result.column_headers.size(); ++i )
    out << (i ? "," : "") << escape(result.column_headers[i]);
  out << eol;

  for( const vector<string> &row : result.rows )
  {
    for( size_t i = 0; i < row.size(); ++i )
      out << (i ? "," : "") << escape(row[i]);
    out << eol;
  }

  return out.str();
}//result_to_csv(...)


/** Reads the activity of one input row: a number followed by an activity unit (e.g. "3.2 uCi"), or a
 bare number, which is becquerel.  Returns false for anything else, including trailing text, and for a
 value that is not a finite, non-negative activity - no measurement gives one, and the decay would
 otherwise silently treat it as zero.  `has_unit` says whether a unit was given.
 */
static bool read_activity( const string &txt, double &activity, bool &has_unit )
{
  has_unit = false;
  try
  {
    activity = PhysicalUnits::stringToActivity( txt );
    has_unit = true;
  }catch( std::exception & )
  {
    try
    {
      size_t end_pos = 0;
      activity = std::stod( txt, &end_pos ) * PhysicalUnits::becquerel;
      if( txt.find_first_not_of( " \t", end_pos ) != string::npos )
        return false;
    }catch( std::exception & )
    {
      return false;
    }
  }//try / catch

  return !IsNan(activity) && !IsInf(activity) && (activity >= 0.0);
}//read_activity(...)


/** #parse_csv, also saying whether the text identifies itself as batch-decay input (see
 #is_candidate_file): it is in one of the column-keyed formats, or gives at least one activity with a
 unit.  A bare "nuclide, number" list could just as well be nuclides and energies.
 */
static vector<BatchNuclide> parse_csv_imp( const string &file_contents, bool &self_identifying )
{
  self_identifying = false;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Nuclear decay database is not available." );

  // Split into lines, tolerating \r\n, \r, and \n.
  string contents = file_contents;

  // Spreadsheet "CSV UTF-8" exports start with a byte-order mark, which would otherwise be read as
  //  part of the first nuclide or column name (and is invisible in any error message).
  if( SpecUtils::starts_with( contents, "\xEF\xBB\xBF" ) )
    contents.erase( 0, 3 );

  SpecUtils::ireplace_all( contents, "\r\n", "\n" );
  SpecUtils::ireplace_all( contents, "\r", "\n" );

  vector<string> raw_lines;
  SpecUtils::split( raw_lines, contents, "\n" );

  // Collect non-blank, non-comment lines.
  vector<string> lines;
  for( string line : raw_lines )
  {
    SpecUtils::trim( line );
    if( line.empty() || line[0] == '#' )
      continue;
    lines.push_back( line );
  }

  if( lines.empty() )
    throw runtime_error( "No data rows found in the file." );

  // Splits a line into trimmed cells, keeping empty ones.  A cell may be double-quoted, as spreadsheets
  //  write one holding the delimiter; the quotes are dropped (a doubled "" too - a literal quote is
  //  rare here).  There are no escape characters, since notes or file paths may hold a backslash.
  auto split_line = []( const string &line, const char * const delims ) -> vector<string> {
    const string no_escapes, quote = "\"";
    const boost::escaped_list_separator<char> separator( no_escapes, delims, quote );
    const boost::tokenizer<boost::escaped_list_separator<char>> tokens( line, separator );

    vector<string> fields;
    for( string field : tokens )
    {
      SpecUtils::trim( field );
      fields.push_back( std::move(field) );
    }
    return fields;
  };

  // A "nuclide, activity" line splits on comma or tab, ignoring empty cells.
  auto split_fields = [&split_line]( const string &line ) -> vector<string> {
    vector<string> fields = split_line( line, ",\t" );
    fields.erase( std::remove( begin(fields), end(fields), string() ), end(fields) );
    return fields;
  };

  // A column-keyed file splits on tab if its header line has one, else comma, keeping empty cells so
  //  a blank one can't shift the later columns - nor a comma within a cell of a tab-separated file.
  const char * const keyed_delim = (lines[0].find( '\t' ) != string::npos) ? "\t" : ",";
  auto split_fields_keep_empty = [&split_line,keyed_delim]( const string &line ) -> vector<string> {
    return split_line( line, keyed_delim );
  };

  // Detect the header-keyed ("Product"/"Value"/"Unit") format from the first line.
  const vector<string> header = split_fields_keep_empty( lines[0] );
  int product_col = -1, value_col = -1, unit_col = -1;
  for( size_t i = 0; i < header.size(); ++i )
  {
    if( SpecUtils::iequals_ascii(header[i], "product") ) product_col = static_cast<int>(i);
    else if( SpecUtils::iequals_ascii(header[i], "value") ) value_col = static_cast<int>(i);
    else if( SpecUtils::iequals_ascii(header[i], "unit") ) unit_col = static_cast<int>(i);
  }

  vector<BatchNuclide> answer;

  // A leading column before "Product" names a physical location, so each location's rows are one
  //  sample set; see the multi-location notes in DecayBatchCalc.h.
  const bool multi_location = (product_col > 0) && (value_col >= 0);

  if( multi_location )
  {

    // Columns other than the location, product, value and unit ride along to the output verbatim.
    vector<int> extra_cols;
    for( size_t i = 0; i < header.size(); ++i )
    {
      const int col = static_cast<int>( i );
      if( (col != 0) && (col != product_col) && (col != value_col) && (col != unit_col) )
        extra_cols.push_back( col );
    }

    for( size_t li = 1; li < lines.size(); ++li )
    {
      const vector<string> fields = split_fields_keep_empty( lines[li] );
      if( static_cast<int>(fields.size()) <= std::max(product_col, value_col) )
        throw runtime_error( "Row '" + lines[li] + "' has too few columns." );

      // Nuclide is the leading token of "Product"; the rest is a description carried to the output.
      const string &product = fields[product_col];
      const string::size_type sp = product.find_first_of( " \t" );
      const string nuc_str = (sp == string::npos) ? product : product.substr( 0, sp );

      BatchNuclide bn;
      bn.nuclide_str = nuc_str;
      bn.nuclide = db->nuclide( nuc_str );
      if( !bn.nuclide )
        throw runtime_error( "'" + nuc_str + "' (from '" + product + "') is not a valid nuclide." );

      bn.location = fields[0];
      bn.product_suffix = (sp == string::npos) ? string() : product.substr( sp );

      // Echo the source file's own column names, so the output header matches the input (and is in
      //  whatever language the input used) rather than being hard-coded English.
      bn.fixed_column_names.push_back( header.empty() ? string() : header[0] );
      bn.fixed_column_names.push_back( header[product_col] );
      bn.fixed_column_names.push_back( header[value_col] );
      bn.fixed_column_names.push_back( (unit_col >= 0 && unit_col < static_cast<int>(header.size()))
                                        ? header[unit_col] : string() );

      for( const int col : extra_cols )
      {
        if( col < static_cast<int>(fields.size()) )
          bn.extra_columns.emplace_back( header[col], fields[col] );
      }

      const string &value_str = fields[value_col];
      const string unit_str = (unit_col >= 0 && unit_col < static_cast<int>(fields.size()))
                                ? fields[unit_col] : string();

      // Separate any areal/volumetric suffix (e.g. "uCi/m2") from the activity unit.
      string act_unit = unit_str;
      const string::size_type slash = unit_str.find( '/' );
      if( slash != string::npos )
      {
        act_unit = unit_str.substr( 0, slash );
        bn.unit_label = unit_str.substr( slash ); // includes the leading '/'
      }
      bn.activity_unit = act_unit;

      // The Value cell is just the number; its unit, if any, is the Unit column's.
      bool has_unit = false;
      if( !read_activity( act_unit.empty() ? value_str : (value_str + " " + act_unit), bn.activity, has_unit )
         || (has_unit == act_unit.empty()) )
      {
        throw runtime_error( "Could not interpret activity '" + value_str
                             + (unit_str.empty() ? string() : (" " + unit_str)) + "' for '" + nuc_str
                             + "'" + sm_activity_requirement );
      }

      answer.push_back( bn );
    }//for( each data line )

    if( answer.empty() )
      throw runtime_error( "No nuclides were parsed from the file." );

    self_identifying = true;

    // Rows of one location differ only by a row-index suffix on the name, so strip it - but only if
    //  doing so can't merge rows that cannot belong to one location (a repeated nuclide, or
    //  disagreeing carried-through cells).  Otherwise keep the names as given.
    std::map<string,vector<size_t>> grouped;
    for( size_t i = 0; i < answer.size(); ++i )
      grouped[ location_key( answer[i].location ) ].push_back( i );

    bool key_is_valid = true;
    for( std::map<string,vector<size_t>>::const_iterator it = begin(grouped);
        key_is_valid && (it != end(grouped)); ++it )
    {
      const vector<size_t> &indices = it->second;
      const BatchNuclide &first = answer[indices.front()];
      std::set<const SandiaDecay::Nuclide *> seen;

      for( size_t i = 0; key_is_valid && (i < indices.size()); ++i )
      {
        const BatchNuclide &cur = answer[indices[i]];
        key_is_valid = seen.insert( cur.nuclide ).second      // no nuclide twice in one location
                       && (cur.extra_columns == first.extra_columns); // and one lat/lon etc per location
      }
    }//for( each candidate group )

    if( key_is_valid )
    {
      for( BatchNuclide &bn : answer )
        bn.location = location_key( bn.location );
    }

    // If no location holds more than one nuclide there is no grouping information in the file, so
    //  treat it as ungrouped (rather than showing one group per row).
    std::map<string,size_t> counts;
    size_t biggest = 0;
    for( const BatchNuclide &bn : answer )
      biggest = std::max( biggest, ++counts[bn.location] );

    if( biggest < 2 )
    {
      for( BatchNuclide &bn : answer )
      {
        bn.location.clear();
        bn.product_suffix.clear();
        bn.extra_columns.clear();
        bn.activity_unit.clear();
      }
    }

    return answer;
  }//if( multi_location )

  if( (product_col >= 0) && (value_col >= 0) )
  {
    // Header-keyed format; data starts at the second line.
    for( size_t li = 1; li < lines.size(); ++li )
    {
      const vector<string> fields = split_fields_keep_empty( lines[li] );
      if( static_cast<int>(fields.size()) <= std::max(product_col, value_col) )
        throw runtime_error( "Row '" + lines[li] + "' has too few columns." );

      // Nuclide is the leading token of the "Product" cell (e.g. "AM-241 Deposition ..." -> "AM-241").
      string product = fields[product_col];
      const string::size_type sp = product.find_first_of( " \t" );
      const string nuc_str = (sp == string::npos) ? product : product.substr(0, sp);

      BatchNuclide bn;
      bn.nuclide_str = nuc_str;
      bn.nuclide = db->nuclide( nuc_str );
      if( !bn.nuclide )
        throw runtime_error( "'" + nuc_str + "' (from '" + product + "') is not a valid nuclide." );

      const string value_str = fields[value_col];
      string unit_str = (unit_col >= 0 && unit_col < static_cast<int>(fields.size()))
                          ? fields[unit_col] : string();

      // Separate any areal/volumetric suffix (e.g. "uCi/m2") from the activity unit.
      string act_unit = unit_str;
      const string::size_type slash = unit_str.find( '/' );
      if( slash != string::npos )
      {
        act_unit = unit_str.substr( 0, slash );
        bn.unit_label = unit_str.substr( slash ); // includes the leading '/'
      }

      // As above: the Value cell is just the number, and with no unit it is becquerel.
      bool has_unit = false;
      if( !read_activity( act_unit.empty() ? value_str : (value_str + " " + act_unit), bn.activity, has_unit )
         || (has_unit == act_unit.empty()) )
      {
        throw runtime_error( "Could not interpret activity '" + value_str
                             + (unit_str.empty() ? string() : (" " + unit_str)) + "' for '" + nuc_str
                             + "'" + sm_activity_requirement );
      }

      answer.push_back( bn );
    }//for( each data line )

    self_identifying = true;
  }else
  {
    // Simple "nuclide, activity[units]" format.
    for( const string &line : lines )
    {
      const vector<string> fields = split_fields( line );
      if( fields.size() < 2 )
        throw runtime_error( "Line '" + line + "' does not have a nuclide and an activity." );

      BatchNuclide bn;
      bn.nuclide_str = fields[0];
      bn.nuclide = db->nuclide( fields[0] );
      if( !bn.nuclide )
        throw runtime_error( "'" + fields[0] + "' is not a valid nuclide." );

      // Accept a value with explicit activity units, or (matching the legacy CLI) a bare number,
      //  which is interpreted as becquerel.
      bool has_unit = false;
      if( !read_activity( fields[1], bn.activity, has_unit ) )
        throw runtime_error( "Could not interpret activity '" + fields[1] + "' for '" + fields[0]
                             + "'" + sm_activity_requirement );
      self_identifying |= has_unit;

      answer.push_back( bn );
    }//for( each line )
  }//if( header-keyed ) / else

  if( answer.empty() )
    throw runtime_error( "No nuclides were parsed from the file." );

  return answer;
}//parse_csv_imp(...)


vector<BatchNuclide> parse_csv( const string &file_contents )
{
  bool self_identifying = false;
  return parse_csv_imp( file_contents, self_identifying );
}//parse_csv(...)


bool is_candidate_file( const string &start_of_file, const bool is_whole_file )
{
  // Binary content is never ours (and a NUL would also end the text early).
  if( start_of_file.find( '\0' ) != string::npos )
    return false;

  // Only complete lines are judged, so a line cut off at the end of the sample cannot fail the file.
  string text = start_of_file;
  if( !is_whole_file )
  {
    const string::size_type last_eol = text.find_last_of( "\r\n" );
    if( last_eol == string::npos )
      return false;
    text.erase( last_eol );
  }

  try
  {
    // It must also say what it is, and hold something to decay.
    bool self_identifying = false;
    const vector<BatchNuclide> inputs = parse_csv_imp( text, self_identifying );
    return self_identifying
           && std::any_of( begin(inputs), end(inputs), []( const BatchNuclide &in ){
             return in.nuclide && !in.nuclide->isStable();
           } );
  }catch( std::exception & )
  {
  }

  return false;
}//is_candidate_file(...)

}//namespace DecayBatchCalc
