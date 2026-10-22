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
#include <limits>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/RapidXmlUtils.hpp"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitUtils.h"

#include "RefLib.h"
#include "PeakCsv.h"
#include "PeakFileIo.h"

using namespace std;

namespace
{
  const char * const sm_n42_marker = "<DHS:InterSpec";
  const char * const sm_n42_end_marker = "</DHS:InterSpec>";
  const char * const sm_spe_marker = "$PEAK_INFO_CSV:";
  const double sm_511 = 510.9989;

  /** The version of the `<DHS:InterSpec>` element (InterSpec's `sm_specMeasSerializationVersion`),
   and of its `<Peaks>` element (`sm_peakXmlSerializationVersion`). */
  const char * const sm_interspec_node_version = "2";
  const char * const sm_peaks_node_version = "2";


  /** Finds the first of `markers` in the file (in blocks, as N42 files can be large); returns false if
   none are present, otherwise the marker's index and file offset.
   */
  /** Splits a CSV line the way InterSpec's peak CSV reader does (boost::escaped_list_separator with
   escape "\\", separators ",\t", and quote "\""): quotes toggle quoting and are removed, and a
   backslash escapes "n", a quote, a separator, or a backslash (anything else throws).  An empty line
   has no fields; a trailing separator gives an empty last field.
   */
  vector<string> split_csv_line( const string &line )
  {
    vector<string> fields;
    if( line.empty() )
      return fields;

    string field;
    bool in_quote = false;
    for( size_t i = 0; i < line.size(); ++i )
    {
      const char c = line[i];
      if( c == '\\' )
      {
        if( ++i == line.size() )
          throw runtime_error( "cannot end with escape" );
        const char next = line[i];
        if( next == 'n' )
          field += '\n';
        else if( (next == '"') || (next == ',') || (next == '\t') || (next == '\\') )
          field += next;
        else
          throw runtime_error( "unknown escape sequence" );
      }else if( (c == ',') || (c == '\t') )
      {
        if( in_quote )
        {
          field += c;
        }else
        {
          fields.push_back( field );
          field.clear();
        }
      }else if( c == '"' )
      {
        in_quote = !in_quote;
      }else
      {
        field += c;
      }
    }//for( loop over characters )

    fields.push_back( field );
    return fields;
  }//split_csv_line(...)


  bool find_first_marker( const string &path, const vector<string> &markers, size_t &which, size_t &offset )
  {
    ifstream in( path.c_str(), ios::in | ios::binary );
    if( !in )
      return false;

    size_t max_len = 0;
    for( const string &m : markers )
      max_len = std::max( max_len, m.size() );

    string buffer;
    size_t buffer_start = 0;
    vector<char> block( 1024*1024 );
    while( in )
    {
      in.read( block.data(), static_cast<streamsize>(block.size()) );
      const streamsize nread = in.gcount();
      if( nread <= 0 )
        break;
      buffer.append( block.data(), static_cast<size_t>(nread) );

      size_t best = string::npos;
      for( size_t i = 0; i < markers.size(); ++i )
      {
        const size_t pos = buffer.find( markers[i] );
        if( pos < best )
        {
          best = pos;
          which = i;
        }
      }

      if( best != string::npos )
      {
        offset = buffer_start + best;
        return true;
      }

      // Keep the tail, in case a marker spans two blocks
      if( buffer.size() >= max_len )
      {
        const size_t ndrop = buffer.size() - max_len + 1;
        buffer.erase( 0, ndrop );
        buffer_start += ndrop;
      }
    }//while( in )

    return false;
  }//find_first_marker(...)


  string read_file_from( const string &path, const size_t offset )
  {
    ifstream in( path.c_str(), ios::in | ios::binary );
    in.seekg( static_cast<streamoff>(offset) );
    stringstream strm;
    strm << in.rdbuf();
    return strm.str();
  }


  PeakFileIo::SampleNumsToPeaks peaks_from_xml( const rapidxml::xml_node<char> *peaks_node,
                                                vector<string> &warnings )
  {
    using namespace rapidxml;

    int version = 1;
    const string version_str = SpecUtils::xml_value_str( peaks_node->first_attribute( "version", 7 ) );
    if( !version_str.empty() && !(stringstream(version_str) >> version) )
      throw runtime_error( "Peaks element has invalid version '" + version_str + "'" );

    map<int,shared_ptr<PeakContinuum>> continuums;
    for( const xml_node<char> *node = peaks_node->first_node( "PeakContinuum", 13 );
         node; node = node->next_sibling( "PeakContinuum", 13 ) )
    {
      auto continuum = make_shared<PeakContinuum>();
      int cont_id;
      continuum->fromXml( node, cont_id );
      continuums[cont_id] = continuum;
    }//for( loop over PeakContinuum nodes )

    const auto sample_nums = []( const xml_node<char> *set_node ) -> set<int> {
      const xml_node<char> *samples_node = set_node->first_node( "SampleNumbers", 13 );
      vector<int> nums;
      if( !samples_node || !SpecUtils::split_to_ints( samples_node->value(), samples_node->value_size(), nums ) )
        throw runtime_error( "Did not find SampleNumbers under PeakSet" );
      return set<int>( begin(nums), end(nums) );
    };

    PeakFileIo::SampleNumsToPeaks answer;
    if( version < 2 )
    {
      // Version 1: each <PeakSet> holds its <Peak> elements
      for( const xml_node<char> *set_node = peaks_node->first_node( "PeakSet", 7 );
           set_node; set_node = set_node->next_sibling( "PeakSet", 7 ) )
      {
        PeakFileIo::PeakDeque peaks;
        for( const xml_node<char> *peak_node = set_node->first_node( "Peak", 4 );
             peak_node; peak_node = peak_node->next_sibling( "Peak", 4 ) )
        {
          auto peak = make_shared<PeakDef>();
          peak->fromXml( peak_node, continuums, &warnings );
          peaks.push_back( peak );
        }
        answer[sample_nums( set_node )] = peaks;
      }//for( loop over PeakSet nodes )
    }else if( version == 2 )
    {
      // Version 2: <Peak> elements have IDs, which each <PeakSet> lists
      map<int,shared_ptr<const PeakDef>> id_to_peak;
      for( const xml_node<char> *peak_node = peaks_node->first_node( "Peak", 4 );
           peak_node; peak_node = peak_node->next_sibling( "Peak", 4 ) )
      {
        auto peak = make_shared<PeakDef>();
        peak->fromXml( peak_node, continuums, &warnings );

        int id = 0;
        const string id_str = SpecUtils::xml_value_str( peak_node->first_attribute( "id", 2 ) );
        if( !(stringstream(id_str) >> id) || id_to_peak.count(id) )
          throw runtime_error( "Peak elements must have a unique, numeric, \"id\" attribute." );
        id_to_peak[id] = peak;
      }//for( loop over Peak nodes )

      for( const xml_node<char> *set_node = peaks_node->first_node( "PeakSet", 7 );
           set_node; set_node = set_node->next_sibling( "PeakSet", 7 ) )
      {
        // InterSpec also saves the peaks its automated search found, as hints for later fits
        const string source = SpecUtils::xml_value_str( set_node->first_attribute( "source", 6 ) );
        if( source != "UserPeaks" )
          continue;

        const xml_node<char> *ids_node = set_node->first_node( "PeakIds", 7 );
        vector<int> ids;
        if( !ids_node || !SpecUtils::split_to_ints( ids_node->value(), ids_node->value_size(), ids ) )
          throw runtime_error( "Did not find PeakIds under PeakSet" );

        PeakFileIo::PeakDeque peaks;
        for( const int id : ids )
        {
          const auto pos = id_to_peak.find( id );
          if( pos == end(id_to_peak) )
            throw runtime_error( "Could not find peak with id " + std::to_string(id) );
          peaks.push_back( pos->second );
        }
        answer[sample_nums( set_node )] = peaks;
      }//for( loop over PeakSet nodes )
    }else
    {
      throw runtime_error( "Unsupported Peaks element version " + std::to_string(version) );
    }//if( version < 2 ) / else

    for( auto &sp : answer )
      std::sort( begin(sp.second), end(sp.second), &PeakDef::lessThanByMeanShrdPtr );

    return answer;
  }//peaks_from_xml(...)


  /** Reads the `<DHS:InterSpec>` element that starts at `offset` in the file. */
  void read_n42_peaks( const string &path, const size_t offset, const SpecUtils::SpecFile &spec,
                       PeakFileIo::FilePeaks &answer )
  {
    using namespace rapidxml;

    string xml = read_file_from( path, offset );
    const size_t end_pos = xml.find( sm_n42_end_marker );
    if( end_pos == string::npos )
      throw runtime_error( "InterSpec element in N42 file is not closed" );
    xml.resize( end_pos + strlen(sm_n42_end_marker) );

    xml_document<char> doc;
    doc.parse<parse_trim_whitespace | allow_sloppy_parse>( &xml[0] );
    const xml_node<char> *interspec_node = doc.first_node();
    if( !interspec_node )
      throw runtime_error( "Could not parse InterSpec element in N42 file" );

    const xml_node<char> *node = interspec_node->first_node( "DisplayedSampleNumbers", 22 );
    vector<int> samples;
    if( node && SpecUtils::split_to_ints( node->value(), node->value_size(), samples ) )
    {
      for( const int sample : samples )
      {
        if( spec.sample_numbers().count( sample ) )
          answer.displayed_samples.insert( sample );
      }
    }//if( DisplayedSampleNumbers )

    answer.display_type = SpecUtils::xml_value_str( interspec_node->first_node( "DisplayType", 11 ) );

    node = interspec_node->first_node( "Peaks", 5 );
    if( node )
      answer.peaks = peaks_from_xml( node, answer.warnings );
  }//read_n42_peaks(...)


  /** Reads the `$PEAK_INFO_CSV:` section of the SPE file that starts at `offset`. */
  void read_spe_peaks( const string &path, const size_t offset, const SpecUtils::SpecFile &spec,
                       const RefLib::Library &lib, PeakFileIo::FilePeaks &answer )
  {
    stringstream section( read_file_from( path, offset ) );
    string line, csv;
    SpecUtils::safe_get_line( section, line ); //the "$PEAK_INFO_CSV:" line
    while( SpecUtils::safe_get_line( section, line ) )
    {
      SpecUtils::trim( line );
      if( SpecUtils::istarts_with( line, "$" ) )
        break;
      csv += line + "\r\n";
    }

    // Like InterSpec, the peaks are for the sum of everything in the file (normally one spectrum)
    const vector<shared_ptr<const SpecUtils::Measurement>> meass = spec.measurements();
    const shared_ptr<const SpecUtils::Measurement> summed = (meass.size() == 1) ? meass[0]
                  : spec.sum_measurements( spec.sample_numbers(), spec.detector_names(), nullptr );
    if( !summed )
      throw runtime_error( "could not sum the spectra" );

    stringstream csvstrm( csv );
    PeakFileIo::PeakDeque peaks;
    for( const PeakDef &p : PeakFileIo::read_peak_csv( csvstrm, summed, lib ) )
      peaks.push_back( make_shared<PeakDef>( p ) );

    answer.peaks[spec.sample_numbers()] = peaks;
  }//read_spe_peaks(...)


  /** The library line best matching the source of a peak CSV row (InterSpec writes e.g. "Cs137",
   "Ba133 (x-ray)", "Th232 (S.E.)", "Pb-xray", or "H(n g)", and the photon energy - for escape peaks,
   the escape-peak energy).  Returns a source without the decay info if the library does not have it.
   Fields have been lower-cased.
   */
  bool source_from_csv( const RefLib::Library &lib, string name, const string &energy_str, PeakDef::Source &answer )
  {
    float energy = 0.0f;
    if( !SpecUtils::parse_float( energy_str.c_str(), energy_str.size(), energy ) || !(energy > 0.0f) )
      return false;

    PeakDef::Source src;
    src.kind = PeakDef::Source::Kind::Nuclide;
    src.gamma_type = PeakDef::NormalGamma;

    const auto remove_suffix = [&name]( const string &suffix ) -> bool {
      if( !SpecUtils::iends_with( name, suffix ) )
        return false;
      name.erase( name.size() - suffix.size() );
      SpecUtils::trim( name );
      return true;
    };

    if( remove_suffix( "(x-ray)" ) )
      src.gamma_type = PeakDef::XrayGamma;
    else if( remove_suffix( "(s.e.)" ) )
      src.gamma_type = PeakDef::SingleEscapeGamma;
    else if( remove_suffix( "(d.e.)" ) )
      src.gamma_type = PeakDef::DoubleEscapeGamma;
    else if( remove_suffix( "-xray" ) )
    {
      src.kind = PeakDef::Source::Kind::Xray;
      src.gamma_type = PeakDef::XrayGamma;
    }

    const size_t paren = name.find( '(' );
    if( (src.kind == PeakDef::Source::Kind::Nuclide) && (paren != string::npos) )
    {
      // The CSV writes a reaction's comma as a space, e.g. "Fe(n n)"
      src.kind = PeakDef::Source::Kind::Reaction;
      const size_t space = name.find( ' ', paren );
      if( space != string::npos )
        name[space] = ',';
    }//if( a reaction )

    if( name.empty() )
      return false;

    // The photon energy; the CSV gives an escape peak's own energy
    double particle_energy = energy;
    if( src.gamma_type == PeakDef::SingleEscapeGamma )
      particle_energy += sm_511;
    else if( src.gamma_type == PeakDef::DoubleEscapeGamma )
      particle_energy += 2.0*sm_511;

    // The CSV and library energies agree to ~0.01 keV; a larger window could pick a different line
    //  when the file's line is not in the (filtered) library.
    const RefLib::Line *best = nullptr;
    double best_de = 0.05;  //keV
    for( const auto &name_src : lib.sources() )
    {
      for( const RefLib::Line &line : name_src.second.lines )
      {
        const PeakDef::Source &ls = line.source;
        if( (ls.kind != src.kind) || !SpecUtils::iequals_ascii( ls.name, name ) )
          continue;

        bool type_ok = false;
        switch( src.gamma_type )
        {
          case PeakDef::XrayGamma:
            type_ok = line.is_xray;
            break;
          case PeakDef::NormalGamma:  //An annihilation peak is written without a marker
            type_ok = !line.is_xray && ((ls.gamma_type == PeakDef::NormalGamma)
                                        || (ls.gamma_type == PeakDef::AnnihilationGamma));
            break;
          case PeakDef::SingleEscapeGamma:
          case PeakDef::DoubleEscapeGamma:
            type_ok = !line.is_xray && (ls.gamma_type != PeakDef::AnnihilationGamma);
            break;
          case PeakDef::AnnihilationGamma:
            type_ok = (ls.gamma_type == PeakDef::AnnihilationGamma);
            break;
        }//switch( src.gamma_type )

        const double de = fabs( ls.particle_energy - particle_energy );
        if( type_ok && (de < best_de) )
        {
          best = &line;
          best_de = de;
        }
      }//for( loop over lines )
    }//for( loop over library sources )

    if( best )
    {
      answer = best->source;
      if( (src.gamma_type == PeakDef::SingleEscapeGamma) || (src.gamma_type == PeakDef::DoubleEscapeGamma) )
        answer.gamma_type = src.gamma_type;
      return true;
    }//if( best )

    // Not in the library: keep the name (capitalized as InterSpec writes it) and energy, without the
    //  decay information InterSpec would need to identify the nuclear transition.
    name[0] = static_cast<char>( std::toupper( static_cast<unsigned char>(name[0]) ) );
    src.name = name;
    src.particle_energy = static_cast<float>( particle_energy );
    answer = src;
    return true;
  }//source_from_csv(...)
}//namespace


namespace PeakFileIo
{

SavedPeaksLocation locate_saved_peaks( const std::string &path )
{
  SavedPeaksLocation loc;
  size_t which = 0;
  if( find_first_marker( path, { sm_n42_marker, sm_spe_marker }, which, loc.offset ) )
    loc.type = (which == 0) ? SavedPeaksLocation::Type::N42 : SavedPeaksLocation::Type::Spe;
  return loc;
}//SavedPeaksLocation locate_saved_peaks(...)


void SampleKeepingSpecFile::cleanup_after_load( const unsigned int flags )
{
  // As InterSpec's SpecMeas::cleanup_after_load does for files it wrote
  const bool keep = m_keep_sample_numbers && !(flags & ReorderSamplesByTime);
  SpecUtils::SpecFile::cleanup_after_load( keep ? (flags | DontChangeOrReorderSamples) : flags );
}


FilePeaks read_file_peaks( const std::string &path, const SavedPeaksLocation &loc,
                           const SpecUtils::SpecFile &spec, const RefLib::Library &lib )
{
  FilePeaks answer;
  switch( loc.type )
  {
    case SavedPeaksLocation::Type::None:
      break;
    case SavedPeaksLocation::Type::N42:
      read_n42_peaks( path, loc.offset, spec, answer );
      break;
    case SavedPeaksLocation::Type::Spe:
      read_spe_peaks( path, loc.offset, spec, lib, answer );
      break;
  }//switch( loc.type )

  return answer;
}//FilePeaks read_file_peaks(...)


void write_n42( std::ostream &out, const SpecUtils::SpecFile &spec, const SampleNumsToPeaks &peaks,
                const std::set<int> &displayed_samples, const std::vector<std::string> &displayed_detectors )
{
  using namespace rapidxml;

  const shared_ptr<xml_document<char>> doc = spec.create_2012_N42_xml();
  xml_node<char> * const RadInstrumentData = doc ? doc->first_node( "RadInstrumentData", 17 ) : nullptr;
  if( !RadInstrumentData )
    throw runtime_error( "Failed to create N42 file." );

  // As InterSpec's SpecMeas::appendSpecMeasStuffToXml; the namespace URI is only an identifier.
  RadInstrumentData->append_attribute( doc->allocate_attribute( "xmlns:DHS", "https://github.com/sandialabs/InterSpec" ) );

  xml_node<char> *interspec_node = doc->allocate_node( node_element, "DHS:InterSpec" );
  RadInstrumentData->append_node( interspec_node );
  interspec_node->append_attribute( doc->allocate_attribute( "version", sm_interspec_node_version ) );

  string samples_str;
  for( const int sample : displayed_samples )
    samples_str += (samples_str.empty() ? "" : " ") + std::to_string( sample );
  interspec_node->append_node( doc->allocate_node( node_element, "DisplayedSampleNumbers",
                                                   doc->allocate_string( samples_str.c_str() ) ) );

  // Names are quoted, so blank names, or leading/trailing spaces, survive the XML
  xml_node<char> *dets_node = doc->allocate_node( node_element, "DisplayedDetectors" );
  interspec_node->append_node( dets_node );
  for( const string &name : displayed_detectors )
  {
    const string quoted = "\"" + name + "\"";
    dets_node->append_node( doc->allocate_node( node_element, "DetectorName", doc->allocate_string( quoted.c_str() ) ) );
  }

  interspec_node->append_node( doc->allocate_node( node_element, "DisplayType",
                                   SpecUtils::descriptionText( SpecUtils::SpectrumType::Foreground ) ) );

  const bool have_peaks = std::any_of( begin(peaks), end(peaks), []( const auto &sp ){ return !sp.second.empty(); } );
  if( !have_peaks )
  {
    rapidxml::print( out, *doc );
    return;
  }

  xml_node<char> *peaks_node = doc->allocate_node( node_element, "Peaks" );
  interspec_node->append_node( peaks_node );
  peaks_node->append_attribute( doc->allocate_attribute( "version", sm_peaks_node_version ) );

  // Port of SpecMeas::addPeaksToXmlHelper: each peak and continuum is written once, and each set of
  //  samples lists the IDs of its peaks.
  map<shared_ptr<const PeakDef>,int> peak_ids;
  map<shared_ptr<PeakContinuum>,int> continuum_ids;
  for( const auto &sp : peaks )
  {
    if( sp.second.empty() )
      continue;

    string ids_str;
    for( const shared_ptr<const PeakDef> &peak : sp.second )
    {
      auto pos = peak_ids.find( peak );
      if( pos == end(peak_ids) )
      {
        const int id = static_cast<int>( peak_ids.size() + 1 );
        pos = peak_ids.insert( { peak, id } ).first;
        xml_node<char> *peak_node = peak->toXml( peaks_node, peaks_node, continuum_ids );
        peak_node->append_attribute( doc->allocate_attribute( "id", doc->allocate_string( std::to_string(id).c_str() ) ) );
      }
      ids_str += (ids_str.empty() ? "" : " ") + std::to_string( pos->second );
    }//for( loop over peaks )

    string nums_str;
    for( const int sample : sp.first )
      nums_str += (nums_str.empty() ? "" : " ") + std::to_string( sample );

    xml_node<char> *set_node = doc->allocate_node( node_element, "PeakSet" );
    peaks_node->append_node( set_node );
    set_node->append_attribute( doc->allocate_attribute( "source", "UserPeaks" ) );
    set_node->append_node( doc->allocate_node( node_element, "SampleNumbers", doc->allocate_string( nums_str.c_str() ) ) );
    set_node->append_node( doc->allocate_node( node_element, "PeakIds", doc->allocate_string( ids_str.c_str() ) ) );
  }//for( loop over sets of samples )

  rapidxml::print( out, *doc );
}//void write_n42(...)


std::string add_spe_peaks( const std::string &spe, const PeakDeque &peaks,
                           const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  const string end_record = "$ENDRECORD:\r\n";
  string answer = spe;
  if( SpecUtils::iends_with( answer, end_record ) )
    answer.erase( answer.size() - end_record.size() );

  if( peaks.empty() || !data )
    return answer + end_record;

  // "$PEAKLABELS:" is (probably) for PeakEasy; InterSpec doesnt read it back.
  stringstream out;
  out << "$PEAKLABELS:\r\n";
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
  for( const shared_ptr<const PeakDef> &peak : peaks )
  {
    const PeakDef::Source &src = peak->source();
    if( peak->userLabel().empty() && !peak->hasSourceGammaAssigned() )
      continue;
    if( !cal || !cal->valid() )
      continue;

    double channel = 0.0;
    try
    {
      channel = cal->channel_for_energy( peak->hasSourceGammaAssigned() ? peak->gammaParticleEnergy() : peak->mean() );
    }catch( std::exception & )
    {
      continue;
    }

    string label = peak->userLabel();
    switch( src.kind )
    {
      case PeakDef::Source::Kind::None:     break;
      case PeakDef::Source::Kind::Nuclide:  label += " " + src.name; break;
      case PeakDef::Source::Kind::Xray:     label += " " + src.name + " xray"; break;
      case PeakDef::Source::Kind::Reaction: label += " " + src.name; break;
    }

    if( src.kind != PeakDef::Source::Kind::Xray )
    {
      switch( src.gamma_type )
      {
        case PeakDef::NormalGamma:                            break;
        case PeakDef::AnnihilationGamma: label += " annih."; break;
        case PeakDef::SingleEscapeGamma: label += " s.e.";   break;
        case PeakDef::DoubleEscapeGamma: label += " d.e.";   break;
        case PeakDef::XrayGamma:         label += " xray";   break;
      }
    }//if( not an element x-ray )

    SpecUtils::ireplace_all( label, "\r\n", " " );
    SpecUtils::ireplace_all( label, "\r", " " );
    SpecUtils::ireplace_all( label, "\n", " " );
    SpecUtils::ireplace_all( label, "\"", "&quot;" );
    SpecUtils::ireplace_all( label, "  ", " " );
    SpecUtils::trim( label );

    if( !label.empty() )
      out << channel << ", \"" << label << "\"\r\n";
  }//for( loop over peaks )

  // "$PEAK_INFO_CSV:" is InterSpec's peak CSV; characters that could confuse SPE parsers are removed
  stringstream csv;
  PeakCsv::write( csv, peaks, data );
  string csv_str = csv.str();
  SpecUtils::ireplace_all( csv_str, "$", " " );
  SpecUtils::ireplace_all( csv_str, ":", " " );
  SpecUtils::ireplace_all( csv_str, "\t", " " );
  out << "$PEAK_INFO_CSV:\r\n" << csv_str << "\r\n";

  return answer + out.str() + end_record;
}//std::string add_spe_peaks(...)


std::vector<PeakDef> read_peak_csv( std::istream &csv,
                                    const std::shared_ptr<const SpecUtils::Measurement> &meas,
                                    const RefLib::Library &lib )
{
  using SpecUtils::trim_copy;
  using SpecUtils::to_lower_ascii_copy;

  if( !meas || !meas->gamma_counts() || (meas->gamma_counts()->size() < 7) )
    throw runtime_error( "input data invalid" );

  const float minenergy = meas->gamma_energy_min();
  const float maxenergy = meas->gamma_energy_max();

  // The header is the first non-empty, non-comment line
  string line;
  while( SpecUtils::safe_get_line( csv, line, 2048 ) )
  {
    SpecUtils::trim( line );
    if( !line.empty() && (line[0] != '#') )
      break;
  }

  if( line.empty() || !csv )
    throw runtime_error( "Failed to get first line" );

  // Columns that may or may not be present are -1 if not
  int mean_index = -1, area_index = -1, fwhm_index = -1;
  int roi_lower_index = -1, roi_upper_index = -1, nuc_index = -1, nuc_energy_index = -1;
  int color_index = -1, label_index = -1, cont_type_index = -1, skew_type_index = -1;
  int cont_coef_index = -1, skew_coef_index = -1, area_uncert_index = -1, peak_type_index = -1;

  {//begin get column indexes
    vector<string> headers;
    for( const string &header : split_csv_line( line ) )
      headers.push_back( to_lower_ascii_copy( trim_copy( header ) ) );

    const auto index_of = [&headers]( const char *name, const vector<string>::const_iterator start ) -> int {
      const auto pos = std::find( start, headers.cend(), name );
      return (pos == headers.cend()) ? -1 : static_cast<int>( pos - headers.cbegin() );
    };

    mean_index = index_of( "centroid", headers.cbegin() );
    area_index = index_of( "net_area", headers.cbegin() );
    fwhm_index = index_of( "fwhm", headers.cbegin() );   //the first FWHM is in keV; the second is %
    if( mean_index < 0 )
      throw runtime_error( "Header did not contain 'Centroid'" );
    if( area_index < 0 )
      throw runtime_error( "Header did not contain 'Net_Area'" );
    if( fwhm_index < 0 )
      throw runtime_error( "Header did not contain 'FWHM'" );

    // The second "Net_Area" column is the uncertainty
    area_uncert_index = index_of( "net_area", headers.cbegin() + area_index + 1 );
    nuc_index = index_of( "nuclide", headers.cbegin() );
    nuc_energy_index = index_of( "photopeak_energy", headers.cbegin() );
    roi_lower_index = index_of( "roi_lower_energy", headers.cbegin() );
    roi_upper_index = index_of( "roi_upper_energy", headers.cbegin() );
    color_index = index_of( "color", headers.cbegin() );
    label_index = index_of( "user_label", headers.cbegin() );
    cont_type_index = index_of( "continuum_type", headers.cbegin() );
    skew_type_index = index_of( "skew_type", headers.cbegin() );
    cont_coef_index = index_of( "continuum_coefficients", headers.cbegin() );
    skew_coef_index = index_of( "skew_coefficients", headers.cbegin() );
    peak_type_index = index_of( "peak_type", headers.cbegin() );
  }//end get column indexes

  vector<PeakDef> answer;

  // Continua read with the pre-version-3 BiLinearStepCDF convention; converted once all the peaks
  //  sharing them are read.
  set<shared_ptr<PeakContinuum>> legacy_bilinear_cdf_continua;

  while( SpecUtils::safe_get_line( csv, line, 2048 ) )
  {
    SpecUtils::trim( line );

    if( SpecUtils::istarts_with( line, "#END " ) || SpecUtils::istarts_with( line, "# END " ) )
      break;

    if( line.empty() || (line[0] == '#') || (!isdigit( static_cast<unsigned char>(line[0]) ) && (line[0] != '+') && (line[0] != '-')) )
      continue;

    vector<string> fields;
    int token_col = 0;
    for( const string &token : split_csv_line( line ) )
    {
      // Like InterSpec, fields are lower-cased, except the label (and here, the color)
      string val = trim_copy( token );
      if( (token_col != label_index) && (token_col != color_index) )
        val = to_lower_ascii_copy( val );
      fields.push_back( val );
      ++token_col;
    }

    const int nfields = static_cast<int>( fields.size() );
    if( (nfields <= mean_index) || (nfields <= area_index) || (nfields <= fwhm_index) )
      continue;

    const auto have = [nfields]( const int index ){ return (index >= 0) && (index < nfields); };

    try
    {
      const float centroid = std::stof( fields[mean_index] );
      const float fwhm = std::stof( fields[fwhm_index] );
      const float area = std::stof( fields[area_index] );

      if( (centroid <= minenergy) || (centroid >= maxenergy) || (fwhm <= 0.0) || (area <= 0.0) )
        continue;

      PeakDef peak( centroid, fwhm/2.35482, area );

      if( have(area_uncert_index) )
      {
        const float area_uncert = std::stof( fields[area_uncert_index] );
        if( area_uncert > 0.0 )
          peak.setPeakAreaUncert( area_uncert );
      }

      if( have(roi_lower_index) && have(roi_upper_index) )
      {
        const float roi_lower = std::max( minenergy, std::stof( fields[roi_lower_index] ) );
        const float roi_upper = std::min( maxenergy, std::stof( fields[roi_upper_index] ) );
        if( (roi_lower >= roi_upper) || (centroid < roi_lower) || (centroid > roi_upper) )
          throw runtime_error( "ROI range invalid." );

        peak.continuum()->setRange( roi_lower, roi_upper );
        peak.continuum()->calc_linear_continuum_eqn( meas, centroid, roi_lower, roi_upper, 3, 3 );
      }else
      {
        const vector<shared_ptr<const PeakDef>> peakv( 1, make_shared<const PeakDef>( peak ) );
        const bool isHPGe = (PeakFitUtils::coarse_resolution_from_peaks( peakv ) == PeakFitUtils::CoarseResolutionType::High);

        double lower_energy, upper_energy;
        findROIEnergyLimits( lower_energy, upper_energy, peak, meas, isHPGe );
        peak.continuum()->setRange( lower_energy, upper_energy );
        peak.continuum()->calc_linear_continuum_eqn( meas, centroid, lower_energy, upper_energy, 3, 3 );
      }//if( the CSV gives the ROI extent ) / else

      if( have(cont_type_index) )
      {
        const string &strval = fields[cont_type_index];

        try
        {
          // A trailing "(vN)" gives the serialization convention the coefficients follow; its
          //  absence means version 2.
          string type_str = strval;
          int cont_version = 2;
          const size_t paren_pos = type_str.rfind( "(v" );
          if( !type_str.empty() && (type_str.back() == ')') && (paren_pos != string::npos) )
          {
            const string ver_str = type_str.substr( paren_pos + 2, type_str.size() - paren_pos - 3 );
            if( !ver_str.empty() && (ver_str.find_first_not_of( "0123456789" ) == string::npos) )
            {
              cont_version = std::stoi( ver_str );
              type_str = type_str.substr( 0, paren_pos );
            }
          }//if( the type may carry a version tag )

          if( cont_version > PeakContinuum::sm_xmlSerializationVersion )
            throw runtime_error( "continuum convention version " + std::to_string(cont_version) + " is too new" );

          const PeakContinuum::OffsetType type = PeakContinuum::str_to_offset_type_str( type_str.c_str(), type_str.size() );
          peak.continuum()->setType( type );

          if( have(cont_coef_index) )
          {
            vector<float> values;
            SpecUtils::split_to_floats( fields[cont_coef_index].c_str(), values, " ,\r\n\t;", false );

            const size_t num_cont_par = PeakContinuum::num_parameters( type );
            if( values.size() == (1 + num_cont_par) )
            {
              const vector<double> dvalues( begin(values) + 1, end(values) );
              peak.continuum()->setParameters( values[0], dvalues, {} );

              static_assert( PeakContinuum::sm_xmlSerializationVersion == 3, "Check legacy conversion" );
              if( (type == PeakContinuum::BiLinearStepCDF) && (cont_version < 3) )
                legacy_bilinear_cdf_continua.insert( peak.continuum() );
            }//if( reference energy, plus the right number of coefficients )
          }//if( have(cont_coef_index) )
        }catch( std::exception & )
        {
          // Continuum left as fit from the data, like InterSpec
        }
      }//if( have(cont_type_index) )

      if( have(skew_type_index) )
      {
        try
        {
          const PeakDef::SkewType type = PeakDef::skew_from_string( fields[skew_type_index] );
          peak.setSkewType( type );

          if( have(skew_coef_index) )
          {
            vector<float> values;
            SpecUtils::split_to_floats( fields[skew_coef_index].c_str(), values, " ,\r\n\t;", false );
            if( values.size() == PeakDef::num_skew_parameters(type) )
            {
              for( size_t i = 0; i < values.size(); ++i )
                peak.set_coefficient( values[i], PeakDef::CoefficientType(PeakDef::SkewPar0 + i) );
            }
          }//if( have(skew_coef_index) )
        }catch( std::exception & )
        {
        }
      }//if( have(skew_type_index) )

      if( have(nuc_index) && have(nuc_energy_index) && !fields[nuc_index].empty() && !fields[nuc_energy_index].empty() )
      {
        PeakDef::Source src;
        if( source_from_csv( lib, fields[nuc_index], fields[nuc_energy_index], src ) )
          peak.setSource( src );
      }

      if( have(color_index) && !fields[color_index].empty() )
        peak.setLineColor( Wt::WColor( fields[color_index] ) );

      if( have(label_index) && !fields[label_index].empty() )
        peak.setUserLabel( fields[label_index] );

      if( have(peak_type_index)
         && (PeakDef::peak_type_from_str( fields[peak_type_index].c_str() ) == PeakDef::DataDefined) )
      {
        // Area is the data, minus the continuum, in the ROI
        peak.setPeakType( PeakDef::DataDefined );
        const shared_ptr<const SpecUtils::EnergyCalibration> cal = meas->energy_calibration();
        const double lx = peak.lowerX(), ux = peak.upperX();

        if( cal && cal->valid() && (peak.continuum()->type() == PeakContinuum::Linear) )
        {
          try
          {
            const double ref_energy = 0.5*(lx + ux);
            double coefficients[2] = { 0.0, 0.0 };
            PeakContinuum::eqn_from_offsets( meas->find_gamma_channel( lx ), meas->find_gamma_channel( ux ),
                                             ref_energy, meas, 3, 3, coefficients[1], coefficients[0] );
            peak.continuum()->setParameters( ref_energy, coefficients, nullptr );
          }catch( std::exception & )
          {
          }
        }//if( linear continuum )

        double datasum = 0.0, continuumsum = 0.0;
        if( cal && cal->valid() && !PeakContinuum::is_peak_cdf_step_continuum( peak.continuum()->type() ) )
        {
          datasum = meas->gamma_integral( lx, ux );
          continuumsum = peak.continuum()->offset_integral( lx, ux, meas, nullptr, 0 );
        }
        double peaksum = datasum - continuumsum;
        peaksum = ((peaksum >= 0.0) && std::isfinite(peaksum)) ? peaksum : 0.0;
        datasum = ((datasum >= 0.0) && std::isfinite(datasum)) ? datasum : 0.0;
        peak.set_coefficient( peaksum, PeakDef::GaussAmplitude );
        peak.set_uncertainty( sqrt(datasum), PeakDef::GaussAmplitude );
      }//if( a data-defined peak )

      // Peaks with the same ROI share a continuum
      if( peak.continuum()->energyRangeDefined() )
      {
        for( PeakDef &p : answer )
        {
          if( p.continuum()->energyRangeDefined()
             && (fabs(p.continuum()->lowerEnergy() - peak.continuum()->lowerEnergy()) < 0.0001)
             && (fabs(p.continuum()->upperEnergy() - peak.continuum()->upperEnergy()) < 0.0001) )
          {
            peak.setContinuum( p.continuum() );
          }
        }
      }//if( new peak has a defined energy range )

      answer.push_back( peak );

      if( answer.size() > 25000 )
        throw runtime_error( "Too many peaks in CSV file (max 25000)." );
    }catch( std::exception &e )
    {
      throw runtime_error( "Invalid value on line '" + line + "', " + string(e.what()) );
    }//try / catch to parse a line into a peak
  }//while( loop over lines )

  if( answer.empty() )
    throw runtime_error( "No peak rows found." );

  // Bring pre-version-3 BiLinearStepCDF continua to the current convention; this needs the total
  //  peak area of each ROI, so is done after all rows are read.
  map<shared_ptr<PeakContinuum>,pair<double,double>> legacy_roi_sums;  //{sum amp, sum amp*CDF(lower)}
  for( const PeakDef &p : answer )
  {
    const shared_ptr<PeakContinuum> cont = std::const_pointer_cast<PeakContinuum>( p.continuum() );
    if( !legacy_bilinear_cdf_continua.count(cont) )
      continue;

    const double amp = std::max( p.amplitude(), 0.0 );
    if( !std::isfinite(amp) )
      continue;

    pair<double,double> &sums = legacy_roi_sums[cont];
    sums.first += amp;
    if( !p.gausPeak() || (p.sigma() <= 0.0) )
      continue;

    const double cdf0 = PeakDists::peak_cdf( cont->lowerEnergy(), p.mean(), p.sigma(), p.skewType(),
                                            p.coefficients() + PeakDef::CoefficientType::SkewPar0 );
    if( std::isfinite(cdf0) )
      sums.second += amp * cdf0;
  }//for( loop over peaks )

  for( const auto &roi : legacy_roi_sums )
  {
    vector<double> values = roi.first->parameters();
    vector<double> uncerts = roi.first->uncertainties();
    PeakContinuum::convert_legacy_bilinear_step_cdf( values, uncerts, roi.second.first, roi.second.second );
    roi.first->setParameters( roi.first->referenceEnergy(), values, uncerts );
  }

  return answer;
}//read_peak_csv(...)

}//namespace PeakFileIo
