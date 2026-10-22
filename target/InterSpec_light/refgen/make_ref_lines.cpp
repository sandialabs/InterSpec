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

// Generates the reference-line library (ref_lines.json) embedded in InterSpec-light.
//
// Usage: make_ref_lines [--transition-children] <sandia.decay.xml> <sandia.reactiongamma.xml> <sources.txt>
//                       <out.json> [min_rel_importance]
//
// Each output entry is in the SpectrumChartD3 reference-line format, plus a few fields the light
//  app uses to assign sources to peaks:
//  { "parent": "Cs137", "kind": "nuclide"|"reaction"|"xray"|"background", "age": "30.1 y",
//    "desc_strs": [...],
//    "desc_child": [...] (with --transition-children: for each desc_strs entry, the daughter of the
//                         decay it describes, or ""; InterSpec needs it to read peak sources from files),
//    "lines": [ { "e": 661.66, "h": 1, "particle": "gamma"|"xray", "desc_ind": 0, "major": 1,
//                 "t": "SE"|"DE"|"annih" (if not a normal photon), "pe": photon energy (if != e),
//                 "src": nuclide name (background set only), "dp": decay parent (if != source),
//                 "na": 1 (not for assigning to peaks) }, ... ] }
//
// Nuclides are decayed to InterSpec's default age, with all progeny; photons below 10 keV are
//  dropped, as are lines that are neither important nor prominent (see `make_entry`).

#include <map>
#include <set>
#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "nlohmann/json.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/ReactionGamma.h"

using namespace std;
using json = nlohmann::json;

namespace
{
#include "soil_transport.inc"

const double sm_min_photon_energy = 10.0;

// A line not meeting `min_rel_importance` is still kept if its importance is at least this fraction
//  of the largest importance within max(10 keV, 10% of its energy), and at least
//  `sm_min_prominent_importance` of the largest importance of the source.
const double sm_prominence_fraction = 0.1;
const double sm_min_prominent_importance = 3.0E-6;

struct Line
{
  double e = 0.0, h = 0.0, pe = 0.0;
  bool xray = false, no_assign = false;
  string type;  // "", "SE", "DE", or "annih"
  string src, dp, dc, desc;  //dp and dc: parent and daughter of the decay that made the photon
};

// Whether to write each description's decay daughter ("desc_child")
bool sm_write_children = false;


/** Port of InterSpec's PeakDef::defaultDecayTime. */
double default_decay_time( const SandiaDecay::Nuclide *nuclide, string &label )
{
  label.clear();
  double decaytime = nuclide->canObtainSecularEquilibrium() ? 10.0*nuclide->secularEquilibriumHalfLife()
                                                             : 7.0*nuclide->promptEquilibriumHalfLife();
  if( (decaytime <= 0.0) || (decaytime < 0.5*nuclide->halfLife) )
  {
    decaytime = 2.5*nuclide->halfLife;
    label = "2.5 HL";
  }

  if( nuclide->decaysToStableChildren() )
  {
    decaytime = 0.0;
    label = "0 s";
  }

  if( (nuclide->halfLife > 100.0*SandiaDecay::year)
     || ((nuclide->halfLife > 2.0*SandiaDecay::year)
         && ((nuclide->atomicNumber == 92) || (nuclide->atomicNumber == 94))) )
  {
    decaytime = 20.0*SandiaDecay::year;
    label = "20 y";
  }

  if( decaytime > 100.0*nuclide->halfLife )
  {
    decaytime = 7.0*nuclide->halfLife;
    label = "7 HL";
  }

  if( label.empty() )
  {
    const pair<double,const char *> units[] = { {SandiaDecay::year, "y"}, {SandiaDecay::day, "d"},
      {SandiaDecay::hour, "h"}, {60.0*SandiaDecay::second, "m"}, {SandiaDecay::second, "s"} };
    for( const auto &u : units )
    {
      if( (decaytime >= u.first) || (u.first == SandiaDecay::second) )
      {
        char buffer[64];
        snprintf( buffer, sizeof(buffer), "%.3g %s", decaytime/u.first, u.second );
        label = buffer;
        break;
      }
    }
  }//if( label.empty() )

  return decaytime;
}//default_decay_time(...)


string transition_desc( const SandiaDecay::Transition *trans )
{
  string desc = trans->parent ? trans->parent->symbol : string();
  if( trans->child )
    desc += " to " + trans->child->symbol;
  desc += string(" via ") + SandiaDecay::to_str( trans->mode );
  return desc;
}


/** Merges lines at the same (0.01 keV rounded) energy, type, and source; the description of the
 largest contributor is kept. */
vector<Line> merge_lines( const vector<Line> &lines )
{
  map<tuple<long long,bool,string,string>,Line> merged;
  map<tuple<long long,bool,string,string>,double> largest;
  for( const Line &l : lines )
  {
    const auto key = make_tuple( llround(100.0*l.e), l.xray, l.type, l.src );
    Line &m = merged[key];
    if( m.h == 0.0 )
    {
      m = l;
      largest[key] = l.h;
      continue;
    }
    m.h += l.h;
    if( l.h > largest[key] )
    {
      largest[key] = l.h;
      m.desc = l.desc;
      m.dp = l.dp;
      m.dc = l.dc;
    }
  }//for( loop over lines )

  vector<Line> answer;
  for( const auto &kv : merged )
    answer.push_back( kv.second );
  std::sort( begin(answer), end(answer), []( const Line &a, const Line &b ){ return a.e < b.e; } );
  return answer;
}//merge_lines(...)


/** Normalizes to the strongest line, drops weak and low-energy lines, and writes the entry.

 Lines are selected by importance, h*sqrt(E) (as InterSpec's RefLineDynamic uses), so strong
 low-energy lines, which shielding often removes, do not crowd out weaker high-energy lines.  A line
 is kept if its importance is at least `min_rel` of the largest, or if it is prominent: the largest,
 or near the largest, in its energy neighborhood (e.g., Pu239's 375 and 414 keV lines, which are
 ~1E-3 of the U L x-rays).
 */
json make_entry( const string &parent, const string &kind, const string &age, vector<Line> lines,
                 const double min_rel )
{
  lines.erase( std::remove_if( begin(lines), end(lines), []( const Line &l ){
    return (l.e < sm_min_photon_energy) || !(l.h > 0.0);
  } ), end(lines) );
  lines = merge_lines( lines );

  double max_h = 0.0, sum_h = 0.0;
  for( const Line &l : lines )
  {
    max_h = std::max( max_h, l.h );
    sum_h += l.h;
  }
  if( max_h <= 0.0 )
    return nullptr;

  vector<double> importance( lines.size() );
  double max_importance = 0.0;
  for( size_t i = 0; i < lines.size(); ++i )
  {
    importance[i] = lines[i].h * std::sqrt( lines[i].e );
    max_importance = std::max( max_importance, importance[i] );
  }

  vector<bool> keep( lines.size(), false );
  for( size_t i = 0; i < lines.size(); ++i )
  {
    if( importance[i] >= min_rel*max_importance )
    {
      keep[i] = true;
      continue;
    }

    if( importance[i] < sm_min_prominent_importance*max_importance )
      continue;

    const double window = std::max( 10.0, 0.1*lines[i].e );
    double nearby_max = 0.0;
    for( size_t j = 0; j < lines.size(); ++j )
    {
      if( fabs(lines[j].e - lines[i].e) <= window )
        nearby_max = std::max( nearby_max, importance[j] );
    }
    keep[i] = (importance[i] >= sm_prominence_fraction*nearby_max);
  }//for( loop over lines )

  // "Major" lines: at least 5% of the strongest, or among the lines making up 75% of the intensity
  vector<double> sorted_h;
  for( const Line &l : lines )
    sorted_h.push_back( l.h );
  std::sort( begin(sorted_h), end(sorted_h), std::greater<double>() );
  double cumulative = 0.0, major_threshold = sorted_h.front();
  for( const double h : sorted_h )
  {
    cumulative += h;
    major_threshold = h;
    if( cumulative >= 0.75*sum_h )
      break;
  }
  major_threshold = std::min( major_threshold, 0.05*max_h );

  json entry;
  entry["parent"] = parent;
  entry["kind"] = kind;
  if( !age.empty() )
    entry["age"] = age;

  vector<string> descs, children;
  json jlines = json::array();
  for( size_t i = 0; i < lines.size(); ++i )
  {
    if( !keep[i] )
      continue;

    const Line &l = lines[i];
    const double h = l.h / max_h;

    json jl;
    jl["e"] = std::round( 100.0*l.e ) / 100.0;
    char hbuf[32];
    snprintf( hbuf, sizeof(hbuf), "%.4g", h );
    jl["h"] = std::stod( hbuf );
    jl["particle"] = l.xray ? "xray" : "gamma";

    if( !l.desc.empty() )
    {
      auto pos = std::find( begin(descs), end(descs), l.desc );
      if( pos == end(descs) )
      {
        pos = descs.insert( end(descs), l.desc );
        children.push_back( l.dc );
      }
      const size_t index = static_cast<size_t>( pos - begin(descs) );
      if( children[index] != l.dc )
        throw logic_error( "Description '" + l.desc + "' used for decays to " + children[index] + " and " + l.dc );
      jl["desc_ind"] = static_cast<int>( index );
    }

    if( l.h >= major_threshold )
      jl["major"] = 1;
    if( !l.type.empty() )
      jl["t"] = l.type;
    if( (l.pe > 0.0) && (fabs(l.pe - l.e) > 0.001) )
      jl["pe"] = std::round( 1000.0*l.pe ) / 1000.0;
    if( !l.src.empty() )
      jl["src"] = l.src;
    if( !l.dp.empty() && (l.dp != parent) && (l.dp != l.src) )
      jl["dp"] = l.dp;
    if( l.no_assign )
      jl["na"] = 1;
    jlines.push_back( jl );
  }//for( loop over lines )

  if( jlines.empty() )
    return nullptr;

  entry["desc_strs"] = descs;
  if( sm_write_children && std::any_of( begin(children), end(children), []( const string &c ){ return !c.empty(); } ) )
    entry["desc_child"] = children;
  entry["lines"] = jlines;
  return entry;
}//make_entry(...)


/** All photon lines of the decay chain of `activities`; `src` is set on every line if non-empty. */
vector<Line> decay_lines( const vector<SandiaDecay::NuclideActivityPair> &activities, const string &src )
{
  vector<Line> lines;
  double annihilation = 0.0;
  for( const SandiaDecay::NuclideActivityPair &nap : activities )
  {
    for( const SandiaDecay::Transition *trans : nap.nuclide->decaysToChildren )
    {
      for( const SandiaDecay::RadParticle &particle : trans->products )
      {
        const double br = nap.activity * trans->branchRatio * particle.intensity;
        switch( particle.type )
        {
          case SandiaDecay::GammaParticle:
          case SandiaDecay::XrayParticle:
          {
            Line l;
            l.e = l.pe = particle.energy;
            l.h = br;
            l.xray = (particle.type == SandiaDecay::XrayParticle);
            l.src = src;
            l.dp = trans->parent ? trans->parent->symbol : string();
            l.dc = trans->child ? trans->child->symbol : string();
            l.desc = l.xray ? ((trans->child ? trans->child->symbol : l.dp) + " x-ray") : transition_desc( trans );
            lines.push_back( l );
            break;
          }

          case SandiaDecay::PositronParticle:
            annihilation += 2.0*br;
            break;

          default:
            break;
        }//switch( particle.type )
      }//for( loop over products )
    }//for( loop over transitions )
  }//for( loop over nuclides in chain )

  if( annihilation > 0.0 )
  {
    Line l;
    l.e = l.pe = 510.9989;
    l.h = annihilation;
    l.type = "annih";
    l.src = src;
    l.desc = "Annihilation";
    lines.push_back( l );
  }

  return lines;
}//decay_lines(...)


json nuclide_entry( const SandiaDecay::Nuclide *nuc, const double min_rel )
{
  string age_label;
  const double age = default_decay_time( nuc, age_label );

  SandiaDecay::NuclideMixture mix;
  mix.addNuclideByActivity( nuc, 1.0E-3*SandiaDecay::curie );
  return make_entry( nuc->symbol, "nuclide", age_label, decay_lines( mix.activity( age ), "" ), min_rel );
}//nuclide_entry(...)


json background_entry( const SandiaDecay::SandiaDecayDataBase &db, const double min_rel )
{
  // Port of InterSpec's getBackgroundRefLines(): relative activities tuned to a representative
  //  background, with photons transported out of a 1 m radius soil sphere (precomputed table).
  struct BackSrc { const char *name; double rel_act; bool secular; };
  const BackSrc srcs[] = {
    { "U238",  0.0004653/410.2892,  false },
    { "Ra226", 0.02515/17990.5430,  false },
    { "U235",  0.001482/14603.0156, false },
    { "Th232", 0.02038/27897.2617,  true  },
    { "K40",   0.1066/6523.8994,    false }
  };

  vector<Line> lines;
  for( const BackSrc &s : srcs )
  {
    const SandiaDecay::Nuclide *nuc = db.nuclide( s.name );
    if( !nuc )
      throw runtime_error( string("Missing background nuclide ") + s.name );
    const double age = (string(s.name) == "K40") ? 0.0
            : 5.0*(s.secular ? nuc->secularEquilibriumHalfLife() : nuc->promptEquilibriumHalfLife());

    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc, s.rel_act, age );
    for( Line l : decay_lines( mix.activity( 0.0 ), s.name ) )
    {
      if( l.type == "annih" ) //InterSpec doesnt include positrons for background
        continue;
      const auto pos = std::lower_bound( begin(integral_energies), end(integral_energies), l.e );
      l.h *= (pos == end(integral_energies)) ? integral_values.back()
                                             : integral_values[pos - begin(integral_energies)];
      lines.push_back( l );
    }
  }//for( loop over background sources )

  double th232_2614 = 0.0, max_h = 0.0;
  for( const Line &l : lines )
  {
    max_h = std::max( max_h, l.h );
    if( (l.src == "Th232") && (fabs(l.e - 2614.53) < 0.05) )
      th232_2614 += l.h;
  }

  Line se, de, annih;
  se.e = 2614.533 - 510.9989;
  se.pe = de.pe = 2614.533;
  se.h = 0.144*th232_2614;
  se.type = "SE";
  se.src = de.src = "Th232";
  se.dp = de.dp = "Tl208";
  se.dc = de.dc = "Pb208";
  se.desc = "Th232 S.E. 2614 keV";
  de.e = 2614.533 - 2.0*510.9989;
  de.h = 0.082*th232_2614;
  de.type = "DE";
  de.desc = "Th232 D.E. 2614 keV";
  annih.e = annih.pe = 510.9989;
  annih.h = 0.052*max_h;
  annih.type = "annih";
  annih.no_assign = true;
  annih.desc = "Annihilation radiation (beta+)";
  lines.push_back( se );
  lines.push_back( de );
  lines.push_back( annih );

  return make_entry( "Background", "background", "", lines, min_rel );
}//background_entry(...)


json xray_entry( const SandiaDecay::Element *el, const double min_rel )
{
  vector<Line> lines;
  for( const SandiaDecay::EnergyIntensityPair &x : el->xrays )
  {
    Line l;
    l.e = l.pe = x.energy;
    l.h = x.intensity;
    l.xray = true;
    l.desc = el->symbol + " x-ray";
    lines.push_back( l );
  }
  return make_entry( el->symbol, "xray", "", lines, min_rel );
}//xray_entry(...)


json reaction_entry( const ReactionGamma &rg, const string &name, const double min_rel )
{
  vector<ReactionGamma::ReactionPhotopeak> peaks;
  const string label = rg.gammas( name, peaks );
  vector<Line> lines;
  for( const ReactionGamma::ReactionPhotopeak &p : peaks )
  {
    Line l;
    l.e = l.pe = p.energy;
    l.h = p.abundance;
    l.desc = p.reaction ? p.reaction->name() : name;
    lines.push_back( l );
  }
  return make_entry( label.empty() ? name : label, "reaction", "", lines, min_rel );
}//reaction_entry(...)
}//namespace


int main( int argc, char **argv )
{
  if( (argc > 1) && (string(argv[1]) == "--transition-children") )
  {
    sm_write_children = true;
    ++argv;
    --argc;
  }

  if( (argc != 5) && (argc != 6) )
  {
    cerr << "Usage: " << argv[0] << " [--transition-children] <sandia.decay.xml> <sandia.reactiongamma.xml>"
         << " <sources.txt> <out.json> [min_rel_importance]" << endl;
    return 1;
  }

  const double min_rel = (argc == 6) ? std::stod( argv[5] ) : 1.0E-3;

  try
  {
    SandiaDecay::SandiaDecayDataBase db( argv[1] );
    const ReactionGamma rg( argv[2], &db );

    ifstream sources( argv[3] );
    if( !sources )
      throw runtime_error( string("Could not open ") + argv[3] );

    json lib = json::array();
    set<string> seen;
    size_t nskipped = 0;
    string line;
    while( std::getline( sources, line ) )
    {
      const size_t pos = line.find_first_not_of( " \t" );
      if( (pos == string::npos) || (line[pos] == '#') )
        continue;
      const size_t space = line.find_first_of( " \t", pos );
      if( space == string::npos )
        continue;
      string kind = line.substr( pos, space - pos );
      string name = line.substr( line.find_first_not_of( " \t", space ) );
      while( !name.empty() && isspace( static_cast<unsigned char>(name.back()) ) )
        name.pop_back();

      if( (kind == "nuclide") && (name.find( '(' ) != string::npos) )
        kind = "reaction";

      json entry;
      try
      {
        if( kind == "nuclide" )
        {
          const SandiaDecay::Nuclide *nuc = db.nuclide( name );
          if( !nuc )
            throw runtime_error( "unknown nuclide" );
          name = nuc->symbol;
          if( !seen.count( name ) )
            entry = nuclide_entry( nuc, min_rel );
        }else if( kind == "xray" )
        {
          const SandiaDecay::Element *el = db.element( name );
          if( !el )
            throw runtime_error( "unknown element" );
          name = el->symbol;
          if( !seen.count( name ) )
            entry = xray_entry( el, min_rel );
        }else if( kind == "reaction" )
        {
          if( !seen.count( name ) )
            entry = reaction_entry( rg, name, min_rel );
        }else if( kind == "background" )
        {
          if( !seen.count( name ) )
            entry = background_entry( db, min_rel );
        }else
        {
          throw runtime_error( "unknown kind '" + kind + "'" );
        }
      }catch( std::exception &e )
      {
        cerr << "Skipping '" << line << "': " << e.what() << endl;
        ++nskipped;
        continue;
      }

      if( !entry.is_null() && !seen.count( entry["parent"].get<string>() ) )
      {
        seen.insert( entry["parent"].get<string>() );
        lib.push_back( entry );
      }
      seen.insert( name );
    }//while( loop over lines in sources.txt )

    ofstream out( argv[4] );
    out << lib.dump();
    if( !out )
      throw runtime_error( string("Failed writing ") + argv[4] );

    size_t nlines = 0;
    for( const json &e : lib )
      nlines += e["lines"].size();
    cout << "Wrote " << lib.size() << " sources (" << nlines << " lines) to " << argv[4]
         << "; skipped " << nskipped << endl;
  }catch( std::exception &e )
  {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }

  return 0;
}//main(...)
