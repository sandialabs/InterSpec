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
#include <tuple>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <functional>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayBatchCalc.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;

namespace
{
  /** The time of a given step, in seconds.  With `num_steps == 1` returns `time_span`; otherwise
   returns `step*time_span/(num_steps-1)` so the points run 0..time_span inclusive.  Matches the
   legacy batch_decay.cpp behavior.
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
}//anonymous namespace


namespace DecayBatchCalc
{

BatchDecayResult decay( const vector<BatchNuclide> &inputs, const BatchDecayOptions &opts )
{
  BatchDecayResult result;

  if( opts.time_span <= 0.0 )
    throw runtime_error( "A positive decay time must be given." );

  if( opts.num_steps < 1 )
    throw runtime_error( "At least one time step is required." );

  if( inputs.empty() )
    throw runtime_error( "No input nuclides were provided." );

  const bool show_progeny = (opts.show_progeny || opts.mix_input);
  const double act_unit = opts.use_curie ? PhysicalUnits::curie : PhysicalUnits::becquerel;

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
    // Sum all valid inputs into a single mixture and co-decay.
    Source src;
    src.mix.reset( new SandiaDecay::NuclideMixture() );
    bool any_added = false;
    bool consistent_unit = true;
    for( const BatchNuclide &in : inputs )
    {
      if( !in.nuclide )
      {
        result.warnings += "Skipped invalid nuclide '" + in.nuclide_str + "'.\n";
        continue;
      }
      if( in.nuclide->isStable() )
      {
        result.warnings += "Skipped stable nuclide '" + in.nuclide_str + "'.\n";
        continue;
      }
      src.mix->addAgedNuclideByActivity( in.nuclide, in.activity, in.age );
      if( !any_added )
        src.unit_label = in.unit_label;      // capture the first added input's label
      else if( in.unit_label != src.unit_label )
        consistent_unit = false;             // subsequent labels must match
      any_added = true;
    }//for( each input )

    if( !any_added )
      throw runtime_error( "No valid, unstable nuclides to decay." );

    if( !consistent_unit )
    {
      result.warnings += "Mixed inputs had differing unit labels; dropped from output.\n";
      src.unit_label.clear();
    }

    sources.push_back( std::move(src) );
  }else
  {
    // Each input is decayed independently.
    for( const BatchNuclide &in : inputs )
    {
      if( !in.nuclide )
      {
        result.warnings += "Skipped invalid nuclide '" + in.nuclide_str + "'.\n";
        continue;
      }
      if( in.nuclide->isStable() )
      {
        result.warnings += "Skipped stable nuclide '" + in.nuclide_str + "'.\n";
        continue;
      }

      Source src;
      src.mix.reset( new SandiaDecay::NuclideMixture() );
      src.parent = in.nuclide;
      src.unit_label = in.unit_label;
      src.mix->addAgedNuclideByActivity( in.nuclide, in.activity, in.age );
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
        for( size_t step = 0; step < opts.num_steps; ++step )
        {
          const double t = time_for_step( step, opts.num_steps, opts.time_span );
          const double act = src.mix->activity( t, nuc );
          row.push_back( SpecUtils::printCompact( act / act_unit, 6 ) );
        }
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
        const double t = time_for_step( step, opts.num_steps, opts.time_span );
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

  result.num_data_cells = result.rows.size() * opts.num_steps;

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


vector<BatchNuclide> parse_csv( const string &file_contents )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Nuclear decay database is not available." );

  // Split into lines, tolerating \r\n, \r, and \n.
  string contents = file_contents;
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

  // Split a line on comma or tab.
  auto split_fields = []( const string &line ) -> vector<string> {
    vector<string> fields;
    SpecUtils::split( fields, line, ",\t" );
    for( string &f : fields )
      SpecUtils::trim( f );
    return fields;
  };

  // Detect the header-keyed ("Product"/"Value"/"Unit") format from the first line.
  const vector<string> header = split_fields( lines[0] );
  int product_col = -1, value_col = -1, unit_col = -1;
  for( size_t i = 0; i < header.size(); ++i )
  {
    if( SpecUtils::iequals_ascii(header[i], "product") ) product_col = static_cast<int>(i);
    else if( SpecUtils::iequals_ascii(header[i], "value") ) value_col = static_cast<int>(i);
    else if( SpecUtils::iequals_ascii(header[i], "unit") ) unit_col = static_cast<int>(i);
  }

  vector<BatchNuclide> answer;

  if( (product_col >= 0) && (value_col >= 0) )
  {
    // Header-keyed format; data starts at the second line.
    for( size_t li = 1; li < lines.size(); ++li )
    {
      const vector<string> fields = split_fields( lines[li] );
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

      try
      {
        if( act_unit.empty() )
        {
          // No activity unit supplied; interpret the bare value as becquerel.
          bn.activity = std::stod( value_str ) * PhysicalUnits::becquerel;
        }else
        {
          bn.activity = PhysicalUnits::stringToActivity( value_str + " " + act_unit );
        }
      }catch( std::exception & )
      {
        throw runtime_error( "Could not interpret activity '" + value_str + " " + unit_str
                             + "' for '" + nuc_str + "'." );
      }

      answer.push_back( bn );
    }//for( each data line )
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
      const string act_str = fields[1];
      try
      {
        bn.activity = PhysicalUnits::stringToActivity( act_str );
      }catch( std::exception & )
      {
        try
        {
          size_t end_pos = 0;
          const double val = std::stod( act_str, &end_pos );
          if( act_str.find_first_not_of( " \t", end_pos ) != string::npos )
            throw runtime_error( "trailing characters" );
          bn.activity = val * PhysicalUnits::becquerel;
        }catch( std::exception & )
        {
          throw runtime_error( "Could not interpret activity '" + act_str
                               + "' for '" + fields[0] + "'." );
        }
      }

      answer.push_back( bn );
    }//for( each line )
  }//if( header-keyed ) / else

  if( answer.empty() )
    throw runtime_error( "No nuclides were parsed from the file." );

  return answer;
}//parse_csv(...)

}//namespace DecayBatchCalc
