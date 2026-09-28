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

#include <cmath>
#include <cctype>
#include <cstdio>
#include <limits>
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakFitUtils.h"

#include "RefLib.h"

using namespace std;
using json = nlohmann::json;

namespace
{
  const double sm_511 = 510.998950;

  PeakDef::SourceGammaType gamma_type_from_str( const string &t )
  {
    if( t == "SE" )    return PeakDef::SingleEscapeGamma;
    if( t == "DE" )    return PeakDef::DoubleEscapeGamma;
    if( t == "annih" ) return PeakDef::AnnihilationGamma;
    return PeakDef::NormalGamma;
  }

  bool same_source( const PeakDef::Source &a, const PeakDef::Source &b )
  {
    return RefLib::same_nuclide( a, b ) && (fabs( a.particle_energy - b.particle_energy ) < 0.01)
           && (a.gamma_type == b.gamma_type);
  }

  /** Lower-cased, without spaces or dashes (as `App.ref.normalize` in web/reflines.js). */
  string normalize_name( const string &name )
  {
    string answer;
    for( const char c : name )
    {
      if( !std::isspace( static_cast<unsigned char>(c) ) && (c != '-') )
        answer += static_cast<char>( std::tolower( static_cast<unsigned char>(c) ) );
    }
    return answer;
  }

  /** How well `line` matches a peak at `mean` of width `sigma` (smaller is better), as InterSpec's
   `assign_nuc_from_ref_lines` scores it: (0.25*sigma + |mean - E|)/intensity, with x-rays
   de-weighted, and if `escapes`, for high-energy gammas the best of the full-energy, S.E., and D.E.
   peaks.  Sets the source (with the gamma type used), and the expected peak energy and intensity.
   Returns infinity for a line too weak to use.
   */
  double score_line( const RefLib::Line &line, const double mean, const double sigma, const bool isHPGe,
                     const bool escapes, PeakDef::Source &src, double &expected_energy, double &intensity )
  {
    // Efficiency of S.E. and D.E. peaks, relative to F.E., for a generic 20% HPGe (energy in keV)
    const auto single_escape_sf = []( const double x ) -> double {
      return std::max( 0.0, (1.8768E-11 *x*x*x) - (9.1467E-08 *x*x) + (2.1565E-04 *x) - 0.16367 );
    };
    const auto double_escape_sf = []( const double x ) -> double {
      return std::max( 0.0, (1.8575E-11 *x*x*x) - (9.0329E-08 *x*x) + (2.1302E-04 *x) - 0.16176 );
    };

    const double pair_prod_thresh = 1255.0; //escape scale factors are negative below this
    const double always_check_escape_thresh = 4000.0;
    const double escape_suppression_factor = 0.5;
    const double xray_suppression_factor = 0.2;

    src = line.source;
    expected_energy = line.energy;
    intensity = line.intensity * (line.is_xray ? xray_suppression_factor : 1.0);
    if( !std::isfinite(intensity) || (intensity < std::numeric_limits<float>::min()) )
      return std::numeric_limits<double>::infinity();

    double dist = (0.25*sigma + fabs(mean - line.energy)) / intensity;

    const double e = line.energy;
    if( escapes && (src.gamma_type == PeakDef::NormalGamma)
       && ((isHPGe && (e > pair_prod_thresh)) || (e > always_check_escape_thresh)) )
    {
      const double se_abundance = intensity * single_escape_sf(e) * escape_suppression_factor;
      const double se_dist = (0.25*sigma + fabs(mean + sm_511 - e)) / se_abundance;
      const double de_abundance = intensity * double_escape_sf(e) * escape_suppression_factor;
      const double de_dist = (0.25*sigma + fabs(mean + 2.0*sm_511 - e)) / de_abundance;

      if( se_dist < dist )
      {
        dist = se_dist;
        intensity = se_abundance;
        expected_energy = e - sm_511;
        src.gamma_type = PeakDef::SingleEscapeGamma;
      }

      if( de_dist < dist )
      {
        dist = de_dist;
        intensity = de_abundance;
        expected_energy = e - 2.0*sm_511;
        src.gamma_type = PeakDef::DoubleEscapeGamma;
      }
    }//if( check for S.E. or D.E. )

    return dist;
  }//score_line(...)
}//namespace


namespace RefLib
{

void Library::load( const json &lib )
{
  if( !lib.is_array() )
    throw runtime_error( "Reference line library must be an array" );

  map<string,Source> sources;
  for( const json &entry : lib )
  {
    Source src;
    src.parent = entry.at( "parent" ).get<string>();
    const string kind = entry.value( "kind", string("nuclide") );
    const vector<string> descs = entry.value( "desc_strs", vector<string>() );
    const vector<string> children = entry.value( "desc_child", vector<string>() );

    for( const json &l : entry.at( "lines" ) )
    {
      if( l.value( "na", 0 ) )  //e.g., annihilation line of the combined background set
        continue;

      Line line;
      line.energy = l.at( "e" ).get<double>();
      line.intensity = l.at( "h" ).get<double>();
      line.is_xray = (l.value( "particle", string() ) == "xray");

      PeakDef::Source &s = line.source;
      s.gamma_type = line.is_xray ? PeakDef::XrayGamma : gamma_type_from_str( l.value( "t", string() ) );
      s.particle_energy = static_cast<float>( l.value( "pe", line.energy ) );
      s.decay_parent = l.value( "dp", string() );

      const size_t desc_ind = l.value( "desc_ind", descs.size() );
      if( desc_ind < children.size() )
        s.decay_child = children[desc_ind];

      if( kind == "xray" )
      {
        s.kind = PeakDef::Source::Kind::Xray;
        s.name = src.parent;
      }else if( kind == "reaction" )
      {
        // The line's description is InterSpec's name of the reaction that made it, e.g. "Fe(n,n)"
        s.kind = PeakDef::Source::Kind::Reaction;
        s.name = (desc_ind < descs.size()) ? descs[desc_ind] : src.parent;
      }else
      {
        // Nuclide x-rays are attributed to the nuclide (as InterSpec does); lines in a combined
        //  set (e.g., "Background") name their actual nuclide.
        s.kind = PeakDef::Source::Kind::Nuclide;
        s.name = l.value( "src", src.parent );
      }

      src.lines.push_back( line );
    }//for( loop over lines )

    const string key = normalize_name( src.parent );
    sources[key] = std::move( src );
  }//for( loop over entries )

  m_sources = std::move( sources );
}//void Library::load( const json &lib )


const Source *Library::find( const string &parent ) const
{
  const auto pos = m_sources.find( normalize_name( parent ) );
  return (pos == end(m_sources)) ? nullptr : &pos->second;
}//const Source *Library::find(...)


unique_ptr<pair<shared_ptr<const PeakDef>,PeakDef::Source>>
  assign_source( PeakDef &peak,
                 const vector<shared_ptr<const PeakDef>> &previous_peaks,
                 const vector<Shown> &shown,
                 const bool isHPGe )
{
  unique_ptr<pair<shared_ptr<const PeakDef>,PeakDef::Source>> other_change;
  if( shown.empty() )
    return other_change;

  const double minx = peak.lowerX();
  const double maxx = peak.upperX();
  const double mean = peak.mean();
  const double sigma = peak.gausPeak() ? peak.sigma() : peak.roiWidth();

  double mindist = std::numeric_limits<double>::max();
  bool found = false;
  PeakDef::Source best_src;
  double best_intensity = 0.0;
  string best_color;

  shared_ptr<const PeakDef> prevpeak;
  double prev_peak_dist = std::numeric_limits<double>::max(), prev_intensity = 0.0;
  PeakDef::Source prev_line_src;

  for( const Shown &sh : shown )
  {
    if( !sh.source )
      continue;

    for( const Line &line : sh.source->lines )
    {
      PeakDef::Source cand;
      double expected_energy = 0.0, intensity = 0.0;
      const double dist = score_line( line, mean, sigma, isHPGe, true, cand, expected_energy, intensity );
      if( (dist >= mindist) || (expected_energy < minx) || (expected_energy > maxx) )
        continue;

      bool currently_used = false;
      for( const shared_ptr<const PeakDef> &pp : previous_peaks )
      {
        if( pp && same_source( pp->source(), cand ) )
        {
          currently_used = true;
          if( dist < prev_peak_dist )
          {
            prevpeak = pp;
            prev_peak_dist = dist;
            prev_intensity = intensity;
            prev_line_src = cand;
          }
          break;
        }
      }//for( loop over previous peaks )

      if( !currently_used )
      {
        found = true;
        prevpeak.reset();
        mindist = dist;
        best_src = cand;
        best_intensity = intensity;
        best_color = sh.color;
      }
    }//for( loop over lines )
  }//for( loop over shown sources )

  if( !found )
    return other_change;

  // An existing peak claimed a line that scored better for this peak; if the relative amplitudes
  //  and intensities disagree, swap the two assignments.
  if( prevpeak )
  {
    const bool prev_amp_smaller = (prevpeak->amplitude() < peak.amplitude());
    const bool prev_intensity_smaller = (prev_intensity < best_intensity);
    if( prev_amp_smaller != prev_intensity_smaller )
    {
      other_change = make_unique<pair<shared_ptr<const PeakDef>,PeakDef::Source>>( prevpeak, best_src );
      best_src = prevpeak->source();
    }
  }//if( prevpeak )

  peak.setSource( best_src );
  if( !best_color.empty() )
    peak.setLineColor( Wt::WColor(best_color) );

  return other_change;
}//assign_source(...)


vector<PeakDef::Source> suggest_sources( const PeakDef &peak, const vector<Shown> &shown )
{
  const double mean = peak.mean();
  const double sigma = peak.gausPeak() ? peak.sigma() : 0.125*peak.roiWidth();

  vector<pair<double,PeakDef::Source>> best; //the best line of each nuclide, element, or reaction
  for( const Shown &sh : shown )
  {
    if( !sh.source )
      continue;

    for( const Line &line : sh.source->lines )
    {
      // Only the lines as shown (as InterSpec's right-click menu): escape peaks of every
      //  high-energy line would add weak lines' escapes
      PeakDef::Source src;
      double expected_energy = 0.0, intensity = 0.0;
      const double dist = score_line( line, mean, sigma, false, false, src, expected_energy, intensity );
      if( !std::isfinite(dist) || (fabs(expected_energy - mean) > 4.0*sigma) || same_nuclide( src, peak.source() ) )
        continue;

      const auto pos = std::find_if( begin(best), end(best), [&src]( const pair<double,PeakDef::Source> &b ){
        return same_nuclide( b.second, src );
      } );
      if( pos == end(best) )
        best.emplace_back( dist, src );
      else if( dist < pos->first )
        *pos = { dist, src };
    }//for( loop over lines )
  }//for( loop over shown sources )

  // Best first (there are only a few, and this is much smaller code than a sort)
  vector<PeakDef::Source> answer;
  while( !best.empty() )
  {
    const auto pos = std::min_element( begin(best), end(best), []( const pair<double,PeakDef::Source> &a, const pair<double,PeakDef::Source> &b ){
      return a.first < b.first;
    } );
    answer.push_back( pos->second );
    best.erase( pos );
  }
  return answer;
}//suggest_sources(...)


optional<PeakDef::Source> source_from_text( const Library &lib, const PeakDef &peak, const string &text,
                                            const Source *&from, string &note )
{
  from = nullptr;
  note.clear();

  string label = SpecUtils::to_lower_ascii_copy( SpecUtils::trim_copy( text ) );
  if( label.empty() || SpecUtils::contains( label, "none" ) || SpecUtils::contains( label, "undef" )
     || (label == "na") )
    return nullopt;

  // Labels mark annihilation lines with "annih."
  bool annihilation = false;
  for( const char *marker : { "annihilation", "annih.", "annih" } )
  {
    if( SpecUtils::contains( label, marker ) )
    {
      SpecUtils::ireplace_all( label, marker, " " );
      annihilation = true;
      break;
    }
  }//for( loop over annihilation markers )

  PeakDef::SourceGammaType type = PeakDef::NormalGamma;
  PeakDef::gammaTypeFromUserInput( label, type );
  if( annihilation && (type == PeakDef::NormalGamma) )
    type = PeakDef::AnnihilationGamma;

  double given_energy = PeakDef::extract_energy_from_peak_source_string( label ); //removes it from label

  // The name: a reaction through its ")"; else, as InterSpec, through the mass number, plus a
  //  following "m", "m2", "meta 2", etc.
  string name = label;
  const size_t paren = label.find( ')' );
  if( paren != string::npos )
  {
    name = label.substr( 0, paren + 1 );
  }else
  {
    const size_t num_start = label.find_first_of( "0123456789" );
    const size_t num_end = (num_start == string::npos) ? string::npos : label.find_first_not_of( "0123456789", num_start );
    if( num_end != string::npos )
    {
      vector<string> after;
      SpecUtils::split( after, label.substr( num_end ), " \t-\n\r" );
      name = label.substr( 0, num_end );
      const string first = after.empty() ? string() : after[0];
      string level;
      if( (first == "m") || (first == "meta") )
        level = (after.size() > 1) ? after[1] : string("1");
      else if( ((first.size() == 2) && (first[0] == 'm')) || ((first.size() == 5) && SpecUtils::starts_with( first, "meta" )) )
        level = first.substr( first.size() - 1 );

      // The library names metastable states e.g. "Tc99m" and "Hf178m2"
      if( level == "1" )
        name += "m";
      else if( (level == "2") || (level == "3") )
        name += "m" + level;
    }//if( there is something after the mass number )
  }//if( reaction ) / else

  const string key = normalize_name( name );
  if( key.empty() )
    throw runtime_error( "No source name in \"" + SpecUtils::trim_copy( text ) + "\"." );

  // The lines of the library source with this name, else of this nuclide/element/reaction in any source
  const auto lines_named = [&lib]( const string &key ) -> vector<pair<const Line *,const Source *>> {
    vector<pair<const Line *,const Source *>> answer;
    if( const Source *src = lib.find( key ) )
    {
      for( const Line &line : src->lines )
        answer.emplace_back( &line, src );
      return answer;
    }

    for( const auto &key_src : lib.sources() )
    {
      for( const Line &line : key_src.second.lines )
      {
        if( normalize_name( line.source.name ) == key )
          answer.emplace_back( &line, &key_src.second );
      }
    }
    return answer;
  };//lines_named lambda

  vector<pair<const Line *,const Source *>> lines = lines_named( key );

  // A whole number under 295 after a name is not read as an energy, as it may be a mass number (e.g.,
  //  "U 235"), so if that is not a source, try it as an energy (e.g., "Pb 75")
  const size_t space = label.find_last_of( ' ' );
  if( lines.empty() && !(given_energy > 0.0) && (space != string::npos) && ((space + 1) < label.size())
     && (label.find_first_of( "0123456789" ) == (space + 1))
     && (label.find_first_not_of( "0123456789", space + 1 ) == string::npos) )
  {
    lines = lines_named( normalize_name( label.substr( 0, space ) ) );
    if( !lines.empty() )
      given_energy = std::stod( label.substr( space + 1 ) );
  }

  if( lines.empty() )
    throw runtime_error( "\"" + SpecUtils::trim_copy( text ) + "\" is not a source in the reference-line library." );

  lines.erase( std::remove_if( begin(lines), end(lines), [type]( const pair<const Line *,const Source *> &l ) -> bool {
    const PeakDef::SourceGammaType line_type = l.first->source.gamma_type;
    switch( type )
    {
      case PeakDef::NormalGamma:
        return (line_type == PeakDef::SingleEscapeGamma) || (line_type == PeakDef::DoubleEscapeGamma);
      case PeakDef::XrayGamma:
        return !l.first->is_xray;
      case PeakDef::AnnihilationGamma:
        return (line_type != PeakDef::AnnihilationGamma);
      case PeakDef::SingleEscapeGamma:
      case PeakDef::DoubleEscapeGamma:
        return l.first->is_xray || (l.first->source.particle_energy <= 2.0*sm_511)
               || ((line_type != PeakDef::NormalGamma) && (line_type != type));
    }//switch( type )
    return true;
  } ), end(lines) );

  if( lines.empty() )
  {
    const char *what = "lines";
    switch( type )
    {
      case PeakDef::NormalGamma:       break;
      case PeakDef::XrayGamma:         what = "x-rays"; break;
      case PeakDef::AnnihilationGamma: what = "annihilation line"; break;
      case PeakDef::SingleEscapeGamma:
      case PeakDef::DoubleEscapeGamma: what = "lines above 1022 keV (which escape peaks need)"; break;
    }
    throw runtime_error( string("The reference-line library has no ") + what + " for \""
                         + SpecUtils::trim_copy( text ) + "\"." );
  }//if( lines.empty() )

  // Energies are of the photon, so an escape peak's is 511 or 1022 keV above the peak
  double energy = given_energy;
  if( !(energy > 0.0) )
  {
    energy = peak.mean();
    if( type == PeakDef::SingleEscapeGamma )
      energy += sm_511;
    else if( type == PeakDef::DoubleEscapeGamma )
      energy += 2.0*sm_511;
  }//if( no energy given )

  const pair<const Line *,const Source *> *best = nullptr;

  // An energy given to the 0.01 keV labels show picks that line (without "x-ray", a gamma if there is one)
  if( given_energy > 0.0 )
  {
    const auto rank = []( const pair<const Line *,const Source *> &l ){
      return make_pair( !l.first->is_xray, l.first->intensity );
    };
    for( const pair<const Line *,const Source *> &l : lines )
    {
      if( (fabs( l.first->source.particle_energy - energy ) < 0.006) && (!best || (rank( l ) > rank( *best ))) )
        best = &l;
    }
  }//if( given_energy > 0.0 )

  if( !best )
  {
    // As PeakDef::findNearestPhotopeak: of the lines within 4 sigma, the one with the smallest
    //  (0.1*window + |dE|)/(relative intensity); if there are none, the nearest line.
    const double window = 4.0 * (peak.gausPeak() ? peak.sigma() : 0.125*peak.roiWidth());
    double max_intensity = 0.0;
    for( const pair<const Line *,const Source *> &l : lines )
    {
      if( fabs( l.first->source.particle_energy - energy ) <= window )
        max_intensity = std::max( max_intensity, l.first->intensity );
    }

    double best_score = std::numeric_limits<double>::infinity();
    for( const pair<const Line *,const Source *> &l : lines )
    {
      const double de = fabs( l.first->source.particle_energy - energy );
      double score = de;
      if( max_intensity > 0.0 )
        score = ((de <= window) && (l.first->intensity > 0.0)) ? ((0.1*window + de) * max_intensity / l.first->intensity)
                                                               : std::numeric_limits<double>::infinity();
      if( score < best_score )
      {
        best_score = score;
        best = &l;
      }
    }//for( loop over lines )

    if( best && !(max_intensity > 0.0) )
    {
      char buffer[256];
      snprintf( buffer, sizeof(buffer), "%s has no line within %.1f keV of %.1f keV in the reference-line"
                " library; assigned its nearest, %.2f keV.", best->first->source.name.c_str(), window,
                energy, best->first->source.particle_energy );
      note = buffer;
    }
  }//if( !best )

  if( !best )
    throw runtime_error( "Could not find a line for \"" + SpecUtils::trim_copy( text ) + "\"." );

  PeakDef::Source answer = best->first->source;
  if( (type == PeakDef::SingleEscapeGamma) || (type == PeakDef::DoubleEscapeGamma) )
    answer.gamma_type = type;
  from = best->second;
  return answer;
}//source_from_text(...)


string source_label( const PeakDef::Source &src )
{
  if( src.kind == PeakDef::Source::Kind::None )
    return "";

  const char *prefix = "", *postfix = "";
  switch( src.gamma_type )
  {
    case PeakDef::SingleEscapeGamma: prefix = "S.E. "; break;
    case PeakDef::DoubleEscapeGamma: prefix = "D.E. "; break;
    case PeakDef::XrayGamma:         postfix = " x-ray"; break;
    case PeakDef::AnnihilationGamma: postfix = " annih."; break;
    case PeakDef::NormalGamma:       break;
  }
  char buffer[128];
  snprintf( buffer, sizeof(buffer), "%s%s %.2f keV%s", prefix, src.name.c_str(), src.particle_energy, postfix );
  return buffer;
}//source_label(...)


bool same_nuclide( const PeakDef::Source &a, const PeakDef::Source &b )
{
  return (a.kind == b.kind) && (a.kind != PeakDef::Source::Kind::None) && SpecUtils::iequals_ascii( a.name, b.name );
}

}//namespace RefLib
