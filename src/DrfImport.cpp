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

#include <string>
#include <vector>
#include <memory>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/DrfImport.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/AngleOutxImport.h"
#include "InterSpec/DetectorEffG2kPar.h"
#include "InterSpec/GadrasDetectorDat.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;


namespace
{
  // Case-insensitive search of the header bytes.
  bool header_contains( const string &header, const string &term )
  {
    return SpecUtils::icontains( header, term );
  }


  // Whether `term1` and `term2` both start fields of one comma/tab-delimited line.
  bool header_line_has_both( const string &header, const string &term1, const string &term2 )
  {
    string text = header.substr( 0, strnlen( header.c_str(), header.size() ) );
    if( SpecUtils::istarts_with( text, "\xEF\xBB\xBF" ) )  //UTF-8 BOM
      text = text.substr( 3 );

    vector<string> lines;
    SpecUtils::split( lines, text, "\r\n" );

    for( const string &line : lines )
    {
      vector<string> fields;
      SpecUtils::split( fields, line, ",\t" );
      if( fields.size() < 2 )
        continue;

      bool has_term1 = false, has_term2 = false;
      for( string field : fields )
      {
        SpecUtils::trim( field );
        if( (field.size() > 1) && (field.front() == '"') )
          field = field.substr( 1 );
        has_term1 |= SpecUtils::istarts_with( field, term1 );
        has_term2 |= SpecUtils::istarts_with( field, term2 );
      }

      if( has_term1 && has_term2 )
        return true;
    }//for( const string &line : lines )

    return false;
  }//header_line_has_both(...)


  /** A legacy GADRAS Detector.dat has no signature - it is a bare table of numbered parameter
   lines, "<index> <value> <fit-flag>  <label>" - so it is recognized by that table's shape: after
   any leading '!' or '#' comments, three or more consecutively-numbered lines (from 0 or 1), each
   carrying at least two more numbers.  Requiring the run to be consecutive AND to start the file
   is what keeps a CSV with a stray numeric row from matching.  The XML variant announces itself.
   */
  bool looks_like_gadras_dat( const string &header, const size_t fileSize )
  {
    if( header_contains( header, "<gamma_detector" ) )
      return true;

    if( (fileSize <= 256) || (fileSize >= 256*1024) )
      return false;

    int run = 0, expect = -1;
    size_t pos = 0;
    bool past_comments = false;
    while( (pos < header.size()) && (run < 3) )
    {
      const size_t eol = std::min( header.find('\n', pos), header.size() );
      string line = header.substr( pos, eol - pos );
      pos = eol + 1;
      SpecUtils::trim( line );

      if( !past_comments )
      {
        if( line.empty() || (line[0] == '!') || (line[0] == '#') )
          continue;
        past_comments = true;
      }

      int idx = 0;
      double value = 0.0, flag = 0.0;
      const bool parsed = (sscanf( line.c_str(), "%d %lf %lf", &idx, &value, &flag ) == 3);

      if( !parsed || ((run == 0) ? ((idx != 0) && (idx != 1)) : (idx != expect)) )
        break;   //the table has to start the file, not appear somewhere in it

      run += 1;
      expect = idx + 1;
    }//while( looking for the parameter table )

    return (run >= 3);
  }//looks_like_gadras_dat(...)


  string file_stem( const string &displayName )
  {
    string stem = SpecUtils::filename( displayName );
    const string ext = SpecUtils::file_extension( stem );
    if( ext.size() < stem.size() )
      stem = stem.substr( 0, stem.size() - ext.size() );
    return stem;
  }


  // Parsers' placeholder names, which say less than the file name does.
  bool is_placeholder_name( const string &name )
  {
    return name.empty() || (name == "DetectorPeakResponse") || (name == "Efficiency CSV")
           || (name == "ANGLE detector");
  }


  void parse_as( DrfImport::ParsedFile &out, const DrfImport::FileKind kind )
  {
    using DrfImport::FileKind;

    const string &data = *out.data;
    istringstream strm( data );

    out.kind = kind;

    switch( kind )
    {
      case FileKind::DrfXml:
      {
        vector<char> buffer( begin(data), end(data) );
        buffer.push_back( '\0' );

        rapidxml::xml_document<char> doc;
        doc.parse<rapidxml::parse_default>( buffer.data() );
        const rapidxml::xml_node<char> * const node = doc.first_node( "DetectorPeakResponse" );
        if( !node )
          throw runtime_error( "No DetectorPeakResponse XML element." );

        auto drf = make_shared<DetectorPeakResponse>();
        drf->fromXml( node );
        if( !drf->isValid() )
          throw runtime_error( "DRF XML does not define a valid detector." );
        out.drfs.push_back( drf );
        break;
      }//case FileKind::DrfXml:

      case FileKind::MakeDrfCsv:
      {
        const shared_ptr<DetectorPeakResponse> drf = DetectorPeakResponse::parseInterSpecRelEffCsv( strm );
        if( !drf || !drf->isValid() )
          throw runtime_error( "Not a Make Detector Response CSV file." );
        out.drfs.push_back( drf );
        break;
      }//case FileKind::MakeDrfCsv:

      case FileKind::MultiDrfCsv:
      case FileKind::GammaQuantCsv:
      {
        vector<string> credits, warnings;
        vector<shared_ptr<DetectorPeakResponse>> drfs;
        if( kind == FileKind::MultiDrfCsv )
          DetectorPeakResponse::parseMultipleRelEffDrfCsv( strm, credits, drfs );
        else
          DetectorPeakResponse::parseGammaQuantRelEffDrfCsv( strm, drfs, credits, warnings );

        for( const shared_ptr<DetectorPeakResponse> &drf : drfs )
        {
          if( drf && drf->isValid() )
            out.drfs.push_back( drf );
        }

        if( out.drfs.empty() )
          throw runtime_error( "No detector efficiency functions found in file." );

        for( const string &s : credits )
        {
          if( !s.empty() )
            out.credits.push_back( s );
        }
        for( const string &s : warnings )
        {
          if( !s.empty() )
            out.notes.push_back( s );
        }
        break;
      }//case FileKind::MultiDrfCsv / GammaQuantCsv:

      case FileKind::IsocsEcc:
      {
        const DetectorPeakResponse::EccParseResult ecc = DetectorPeakResponse::parseEccFile( strm );
        if( !ecc.drf || !ecc.drf->isValid() )
          throw runtime_error( "Invalid ECC file." );
        out.drfs.push_back( ecc.drf );
        out.sourceArea = ecc.sourceArea;
        out.sourceMass = ecc.sourceMass;
        out.uncertEnergies = ecc.uncertEnergies;
        out.baselineFrac = ecc.baselineFrac;
        out.convergenceFrac = ecc.convergenceFrac;
        break;
      }//case FileKind::IsocsEcc:

      case FileKind::Angle:
      {
        // A file with results gives a fixed-geometry curve; one without (e.g., a .detx) may still
        //  describe the detector well enough to characterize.
        shared_ptr<AngleOutxContents> contents;
        try
        {
          contents = make_shared<AngleOutxContents>( DetectorPeakResponse::parseAngleOutxFileFull( strm ) );
          if( contents->fixedGeomDrf && contents->fixedGeomDrf->isValid() )
            out.drfs.push_back( contents->fixedGeomDrf );
        }catch( std::exception & )
        {
          istringstream geomstrm( data );
          contents = make_shared<AngleOutxContents>( AngleOutx::parse( geomstrm ) );
          if( !contents->hasGeometry && !contents->hasReference )
            throw;
        }//try / catch

        out.notes.insert( end(out.notes), begin(contents->parseNotes), end(contents->parseNotes) );
        out.angle = contents;
        break;
      }//case FileKind::Angle:

      case FileKind::EfficiencyCsv:
      case FileKind::GadrasEfficiencyCsv:
      {
        const DetectorPeakResponse::EffCsvParseResult result
                                            = DetectorPeakResponse::parseEfficiencyCsvFile( strm );
        if( !result.drf || !result.drf->isValid() )
          throw runtime_error( "Not an efficiency CSV file." );
        out.kind = result.is_gadras_format ? FileKind::GadrasEfficiencyCsv : FileKind::EfficiencyCsv;
        out.drfs.push_back( result.drf );
        break;
      }//case FileKind::EfficiencyCsv / GadrasEfficiencyCsv:

      case FileKind::GadrasDetectorDat:
      {
        if( !GadrasDetectorDat::isCandidateDetectorDat( strm ) )
          throw runtime_error( "Not a GADRAS Detector.dat file." );

        strm.clear();
        strm.seekg( 0, ios::beg );
        const GadrasDetectorDat dat = GadrasDetectorDat::fromStream( strm );

        // Everything the file defines except an efficiency; applyGadrasDat attaches the geometry,
        //  but swallows why one could not be built, so ask for that here.
        strm.clear();
        strm.seekg( 0, ios::beg );
        auto seed = make_shared<DetectorPeakResponse>();
        seed->fromGadrasDatOnly( strm );

        try
        {
          vector<string> warnings;
          CeeLoUtils::buildGadrasGeometry( dat, warnings );
          out.notes.insert( end(out.notes), begin(warnings), end(warnings) );
        }catch( std::exception &e )
        {
          out.notes.push_back( e.what() );
        }

        out.drfs.push_back( seed );
        break;
      }//case FileKind::GadrasDetectorDat:

      case FileKind::ParGrid:
      {
        const vector<uint8_t> bytes( begin(data), end(data) );
        out.par = make_shared<const DetEffG2kPar::ParFile>( DetEffG2kPar::parseParFile( bytes ) );
        break;
      }//case FileKind::ParGrid:

      case FileKind::ParDetectorTxt:
      {
        for( const DetEffG2kPar::DetectorDef &def : DetEffG2kPar::parseDetectorTxt( strm ) )
        {
          if( (def.d1_crystal_diam_mm > 0.0) && (def.d2_crystal_len_mm > 0.0) )
            out.detectorDefs.push_back( def );
        }

        if( out.detectorDefs.empty() )
          throw runtime_error( "No detector records found in file." );
        break;
      }//case FileKind::ParDetectorTxt:
    }//switch( kind )
  }//parse_as(...)
}//namespace


namespace DrfImport
{

vector<FileKind> candidateKinds( const uint8_t *header, const size_t headerLen, const size_t fileSize )
{
  vector<FileKind> answer;
  if( !header || !headerLen )
    return answer;

  // Binary first: a .par grid has no signature, only a plausible header.
  if( DetEffG2kPar::isCandidateParFile( header, headerLen, fileSize ) )
    answer.push_back( FileKind::ParGrid );

  string text( reinterpret_cast<const char *>(header), headerLen );
  if( SpecUtils::istarts_with( text, "\xEF\xBB\xBF" ) )  //UTF-8 BOM
    text = text.substr( 3 );

  if( header_contains( text, "<DetectorPeakResponse" ) )
    answer.push_back( FileKind::DrfXml );

  if( header_contains( text, "# Detector Response Function" ) )
    answer.push_back( FileKind::MakeDrfCsv );

  if( header_contains( text, "Relative Eff" ) )
    answer.push_back( FileKind::MultiDrfCsv );

  if( header_contains( text, "Detector ID" ) )
    answer.push_back( FileKind::GammaQuantCsv );

  if( header_contains( text, "SGI_template" ) || header_contains( text, "ISOCS_file_name" ) )
    answer.push_back( FileKind::IsocsEcc );

  if( header_contains( text, "<angle" ) )
    answer.push_back( FileKind::Angle );

  // parseFile tells a GADRAS Efficiency.csv from other energy/efficiency CSVs
  if( header_line_has_both( text, "en", "eff" ) || header_contains( text, "energy,peak,pcom" ) )
    answer.push_back( FileKind::EfficiencyCsv );

  if( looks_like_gadras_dat( text, fileSize ) )
    answer.push_back( FileKind::GadrasDetectorDat );

  if( DetEffG2kPar::isCandidateDetectorTxt( text ) )
    answer.push_back( FileKind::ParDetectorTxt );

  return answer;
}//candidateKinds(...)


bool isCompleteDrfKind( const FileKind kind )
{
  switch( kind )
  {
    case FileKind::DrfXml:
    case FileKind::MakeDrfCsv:
    case FileKind::MultiDrfCsv:
    case FileKind::GammaQuantCsv:
      return true;

    case FileKind::IsocsEcc:
    case FileKind::Angle:
    case FileKind::EfficiencyCsv:
    case FileKind::GadrasEfficiencyCsv:
    case FileKind::GadrasDetectorDat:
    case FileKind::ParGrid:
    case FileKind::ParDetectorTxt:
      break;
  }//switch( kind )

  return false;
}//isCompleteDrfKind(...)


bool areCompanions( const FileKind a, const FileKind b )
{
  const auto is_pair = [a,b]( const FileKind x, const FileKind y ){
    return ((a == x) && (b == y)) || ((a == y) && (b == x));
  };

  return is_pair( FileKind::GadrasEfficiencyCsv, FileKind::GadrasDetectorDat )
         || is_pair( FileKind::EfficiencyCsv, FileKind::GadrasDetectorDat )
         || is_pair( FileKind::ParGrid, FileKind::ParDetectorTxt );
}//areCompanions(...)


shared_ptr<const string> readFile( const string &path )
{
#ifdef _WIN32
  const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16( path );
  ifstream input( wpath.c_str(), ios::in | ios::binary );
#else
  ifstream input( path.c_str(), ios::in | ios::binary );
#endif

  if( !input )
    throw runtime_error( "Could not open file." );

  input.seekg( 0, ios::end );
  const streamoff size = input.tellg();
  input.seekg( 0, ios::beg );

  if( size < 0 )
    throw runtime_error( "Could not read file." );
  if( static_cast<size_t>(size) > sm_maxFileBytes )
    throw runtime_error( "File is too large to be a detector efficiency file." );

  auto data = make_shared<string>( static_cast<size_t>(size), '\0' );
  if( size && !input.read( &((*data)[0]), size ) )
    throw runtime_error( "Could not read file." );

  return data;
}//readFile(...)


shared_ptr<const ParsedFile> parseFile( const string &displayName,
                                        shared_ptr<const string> data,
                                        const bool anyTextAsCsv )
{
  if( !data || data->empty() )
    throw runtime_error( "Empty file." );

  if( data->size() > sm_maxFileBytes )
    throw runtime_error( "File is too large to be a detector efficiency file." );

  // The parsers don't expect a UTF-8 BOM (a binary .par could start with these bytes, but then
  //  its header would not be plausible anyway)
  if( SpecUtils::istarts_with( *data, "\xEF\xBB\xBF" ) )
    data = make_shared<const string>( data->substr( 3 ) );

  const size_t headerLen = std::min( data->size(), size_t(1024) );
  vector<FileKind> kinds = candidateKinds( reinterpret_cast<const uint8_t *>( data->data() ),
                                           headerLen, data->size() );

  // The header only nominates; a CSV's column headings can be further in than it looks.
  if( anyTextAsCsv && (std::find( begin(kinds), end(kinds), FileKind::ParGrid ) == end(kinds)) )
  {
    for( const FileKind kind : { FileKind::MakeDrfCsv, FileKind::EfficiencyCsv, FileKind::MultiDrfCsv } )
    {
      if( std::find( begin(kinds), end(kinds), kind ) == end(kinds) )
        kinds.push_back( kind );
    }
  }

  string first_error;
  for( const FileKind kind : kinds )
  {
    auto answer = make_shared<ParsedFile>();
    answer->displayName = displayName;
    answer->data = data;

    try
    {
      parse_as( *answer, kind );
      return answer;
    }catch( std::exception &e )
    {
      if( first_error.empty() )
        first_error = e.what();
    }
  }//for( const FileKind kind : kinds )

  if( kinds.empty() || first_error.empty() )
    throw runtime_error( "Not a recognized detector efficiency file." );

  throw runtime_error( first_error );
}//parseFile(...)


bool isSlow( const ParsedFile *a, const ParsedFile *b )
{
  return a && b && areCompanions( a->kind, b->kind )
         && ((a->kind == FileKind::ParGrid) || (b->kind == FileKind::ParGrid));
}//isSlow(...)


Source makeSource( shared_ptr<const ParsedFile> a, shared_ptr<const ParsedFile> b, const size_t record )
{
  Source src;

  if( !a )
  {
    src.error = "No file.";
    return src;
  }

  if( b && !areCompanions( a->kind, b->kind ) )
    b.reset();

  // `a` is the efficiency-bearing half of a pair.
  if( b && ((a->kind == FileKind::GadrasDetectorDat) || (a->kind == FileKind::ParDetectorTxt)) )
    std::swap( a, b );

  src.notes = a->notes;
  src.credits = a->credits;
  if( b )
    src.notes.insert( end(src.notes), begin(b->notes), end(b->notes) );

  const vector<Interpretation> curve_interps{ Interpretation::FarFieldIntrinsic,
                                              Interpretation::FarFieldAbsolute,
                                              Interpretation::FixedTotal };

  try
  {
    switch( a->kind )
    {
      case FileKind::DrfXml:
      case FileKind::MakeDrfCsv:
      case FileKind::MultiDrfCsv:
      case FileKind::GammaQuantCsv:
      {
        if( a->drfs.size() > 1 )
        {
          for( const shared_ptr<const DetectorPeakResponse> &drf : a->drfs )
            src.recordNames.push_back( drf->name() );
        }

        src.record = std::min( record, a->drfs.size() - 1 );
        src.base = a->drfs.at( src.record );
        src.status = Status::Ready;
        break;
      }//case complete DRF(s)

      case FileKind::IsocsEcc:
      {
        src.base = a->drfs.at( 0 );
        src.interpretations = curve_interps;
        if( a->sourceArea > 0.0 )
        {
          src.interpretations.push_back( Interpretation::FixedPerCm2 );
          src.interpretations.push_back( Interpretation::FixedPerM2 );
        }
        if( a->sourceMass > 0.0 )
          src.interpretations.push_back( Interpretation::FixedPerGram );
        src.defaultInterpretation = Interpretation::FixedTotal;
        src.sourceArea = a->sourceArea;
        src.sourceMass = a->sourceMass;
        src.uncertEnergies = a->uncertEnergies;
        src.baselineFrac = a->baselineFrac;
        src.convergenceFrac = a->convergenceFrac;
        src.status = Status::Ready;
        break;
      }//case FileKind::IsocsEcc:

      case FileKind::Angle:
      {
        if( !a->angle )
          throw runtime_error( "ANGLE file was not parsed." );
        const AngleOutxContents &contents = *a->angle;

        // A file describing the whole detector is best used as one - geometry-modeled, answering
        //  any source position - so that is the default when it can be built.
        shared_ptr<const ceelo::GeometryDescriptor> geometry;
        if( contents.hasGeometry && contents.modeASupported )
        {
          vector<string> warnings;
          try
          {
            geometry = make_shared<const ceelo::GeometryDescriptor>(
                                        CeeLoUtils::buildAngleGeometry( contents, warnings ) );
          }catch( std::exception &e )
          {
            src.notes.push_back( e.what() );
          }
          src.notes.insert( end(src.notes), begin(warnings), end(warnings) );
        }else if( contents.hasGeometry && !contents.modeAObstruction.empty() )
        {
          src.notes.push_back( contents.modeAObstruction );
        }

        if( geometry && !contents.hasReference )
          src.notes.push_back( "The file has no measured reference efficiency curve for a detector"
                               " model to be anchored on, so it cannot be used as a generic detector." );

        if( geometry && contents.hasReference )
        {
          try
          {
            const shared_ptr<DetectorPeakResponse> generic = CeeLoUtils::buildAngleSeedDrf( contents );
            generic->setGeometry( geometry );
            CeeLoUtils::attachCurveTransferResponse( *generic );
            src.generic = generic;
          }catch( std::exception &e )
          {
            src.notes.push_back( e.what() );
          }
        }//if( geometry && contents.hasReference )

        if( !a->drfs.empty() )
        {
          src.base = a->drfs[0];
          src.interpretations = curve_interps;
        }

        if( src.generic )
        {
          src.interpretations.push_back( Interpretation::GenericDetector );
          src.defaultInterpretation = Interpretation::GenericDetector;
          src.status = Status::Ready;
        }else if( src.base )
        {
          src.defaultInterpretation = Interpretation::FixedTotal;
          src.status = Status::Ready;
        }else if( geometry )
        {
          // Geometry only (e.g., a .detx): a seed for Monte Carlo characterization.
          const string name = contents.detName.empty() ? file_stem( a->displayName ) : contents.detName;
          auto seed = make_shared<DetectorPeakResponse>( name, contents.detDescription );
          seed->setDetectorDiameter( static_cast<float>( 2.0 * contents.crystalRadius ) );
          seed->setGeometry( geometry );
          src.base = seed;
          src.status = Status::NeedsCharacterization;
        }else
        {
          throw runtime_error( "The ANGLE file has neither efficiency results nor a usable detector"
                               " description." );
        }
        break;
      }//case FileKind::Angle:

      case FileKind::EfficiencyCsv:
      case FileKind::GadrasEfficiencyCsv:
      {
        if( b )
        {
          // With its Detector.dat: FWHM, peak shape, crystal geometry and a curve-transfer
          //  response all come along (see DetectorPeakResponse::applyGadrasDat), and - as GADRAS
          //  does - the efficiencies are taken as intrinsic.
          istringstream csv( *a->data ), dat( *b->data );
          auto drf = make_shared<DetectorPeakResponse>();
          drf->fromGadrasDefinition( csv, dat );
          drf->setDrfSource( DetectorPeakResponse::DrfSource::UserImportedGadrasDrf );
          src.base = drf;
          if( a->kind == FileKind::EfficiencyCsv )
            src.notes.push_back( "The efficiencies are taken to be intrinsic, as they are for a"
                                 " GADRAS Efficiency.csv." );
        }else
        {
          src.base = a->drfs.at( 0 );
          src.interpretations = curve_interps;
          if( a->kind == FileKind::GadrasEfficiencyCsv )
          {
            // GADRAS efficiencies are intrinsic, but the diameter is in the Detector.dat
            src.defaultInterpretation = Interpretation::FarFieldIntrinsic;
          }else
          {
            // The file's efficiencies may already be per unit area or mass of source
            src.interpretations.push_back( Interpretation::FixedPerCm2 );
            src.interpretations.push_back( Interpretation::FixedPerM2 );
            src.interpretations.push_back( Interpretation::FixedPerGram );
            src.defaultInterpretation = Interpretation::FixedTotal;
          }
        }//if( b ) / else
        src.status = Status::Ready;
        break;
      }//case FileKind::EfficiencyCsv / GadrasEfficiencyCsv:

      case FileKind::GadrasDetectorDat:
      {
        src.base = a->drfs.at( 0 );
        if( !src.base->geometry() )
          throw runtime_error( "The Detector.dat does not describe a detector geometry that can be"
                               " characterized; its Efficiency.csv is needed." );
        src.status = Status::NeedsCharacterization;
        break;
      }//case FileKind::GadrasDetectorDat:

      case FileKind::ParGrid:
      case FileKind::ParDetectorTxt:
      {
        if( !b )
        {
          src.status = Status::NeedsCompanion;
          break;
        }

        if( (a->kind != FileKind::ParGrid) || !a->par || (b->kind != FileKind::ParDetectorTxt) )
          throw runtime_error( "The .par grid or its DETECTOR.txt was not parsed." );
        const vector<DetEffG2kPar::DetectorDef> &defs = b->detectorDefs;

        DetEffG2kPar::DetectorDef def;
        try
        {
          def = DetEffG2kPar::selectDetectorDef( defs, a->displayName );
        }catch( std::exception & )
        {
          // No record names this .par; let the user pick.
          for( const DetEffG2kPar::DetectorDef &d : defs )
            src.recordNames.push_back( d.name + " (" + d.parFile + ")" );
          src.record = std::min( record, defs.size() - 1 );
          def = defs.at( src.record );
          src.notes.push_back( "None of the records in '" + b->displayName + "' are for '"
                               + SpecUtils::filename( a->displayName ) + "'; choose which to use." );
        }//try / catch

        src.base = DetEffG2kPar::makeDrf( *a->par, def );
        src.status = Status::Ready;
        break;
      }//case FileKind::ParGrid / ParDetectorTxt:
    }//switch( a->kind )
  }catch( std::exception &e )
  {
    src.status = Status::Error;
    src.error = e.what();
    src.base.reset();
    src.generic.reset();
  }//try / catch

  const shared_ptr<const DetectorPeakResponse> named = src.generic ? src.generic : src.base;
  src.name = named ? named->name() : string();
  if( is_placeholder_name( src.name ) )
  {
    src.name = file_stem( a->displayName );
    src.nameIsFileStem = true;
  }

  return src;
}//makeSource(...)


Result build( const Source &src, const Options &options )
{
  Result result;
  result.status = src.status;
  result.error = src.error;

  if( src.status != Status::Ready )
    return result;

  Interpretation interp = Interpretation::AsIs;
  if( !src.interpretations.empty() )
  {
    const bool allowed = (std::find( begin(src.interpretations), end(src.interpretations),
                                     options.interpretation ) != end(src.interpretations));
    interp = allowed ? options.interpretation : src.defaultInterpretation;
  }

  const bool far_field = ((interp == Interpretation::FarFieldIntrinsic)
                          || (interp == Interpretation::FarFieldAbsolute));
  if( far_field && !(options.diameter > 0.0) )
  {
    result.status = Status::NeedsDiameter;
    return result;
  }

  if( (interp == Interpretation::FarFieldAbsolute) && !(options.distance > 0.0) )
  {
    result.status = Status::NeedsDistance;
    return result;
  }

  try
  {
    const shared_ptr<const DetectorPeakResponse> base = (interp == Interpretation::GenericDetector)
                                                        ? src.generic : src.base;
    if( !base )
      throw runtime_error( "No detector to build from." );

    shared_ptr<DetectorPeakResponse> drf;
    switch( interp )
    {
      case Interpretation::AsIs:
      case Interpretation::GenericDetector:
        drf = make_shared<DetectorPeakResponse>( *base );
        break;

      case Interpretation::FarFieldIntrinsic:
        drf = base->reinterpretAsFarFieldIntrinsicEfficiency( options.diameter );
        break;

      case Interpretation::FarFieldAbsolute:
        drf = base->reinterpretAsFarFieldAbsEfficiency( options.diameter, options.distance, true );
        break;

      case Interpretation::FixedTotal:
        drf = base->reinterpretAsFixedGeom( DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
        break;

      case Interpretation::FixedPerCm2:
      case Interpretation::FixedPerM2:
      case Interpretation::FixedPerGram:
      {
        const DetectorPeakResponse::EffGeometryType type
          = (interp == Interpretation::FixedPerCm2) ? DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2
          : (interp == Interpretation::FixedPerM2)  ? DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2
                                                    : DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram;
        const double quantity = (interp == Interpretation::FixedPerGram) ? src.sourceMass : src.sourceArea;

        // A file stating its source's area/mass (an ISOCS .ecc) is for the total activity, so is
        //  converted; otherwise the efficiencies are taken to already be per unit area/mass.
        if( quantity > 0.0 )
          drf = base->reinterpretAsFixedGeom( DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct )
                    ->convertFixedGeometryType( quantity, type );
        else
          drf = base->reinterpretAsFixedGeom( type );
        break;
      }//case per unit area or mass
    }//switch( interp )

    if( !drf )
      throw runtime_error( "Failed to create detector." );

    // The reinterpret* functions may hand back `base` itself when there is nothing to change.
    if( drf.get() == base.get() )
      drf = make_shared<DetectorPeakResponse>( *base );

    if( far_field && (options.setback >= 0.0) )
      drf->setDetectorSetback( options.setback );

    if( options.setUncert )
      drf->setEfficiencyUncert( options.uncert );

    drf->setName( options.name.empty() ? src.name : options.name );

    result.drf = drf;
    result.status = Status::Ready;
  }catch( std::exception &e )
  {
    result.status = Status::Error;
    result.error = e.what();
  }//try / catch

  return result;
}//build(...)

}//namespace DrfImport
