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

/* Tests of `DrfImport` - identifying detector-efficiency files, pairing a file with its companion
 (Efficiency.csv + Detector.dat, .par + DETECTOR.txt), and building DRFs from the user's choices -
 which is everything the DRF import GUI does except the widgets.
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE TestDrfImport
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "io/DetectorResponse.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfImport.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorEffG2kPar.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using DrfImport::FileKind;
using DrfImport::Status;
using DrfImport::Interpretation;

namespace
{
string g_data_dir, g_test_file_dir;

void set_dirs()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;
  s_have_set = true;

  const int argc = boost::unit_test::framework::master_test_suite().argc;
  char ** const argv = boost::unit_test::framework::master_test_suite().argv;

  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      g_data_dir = arg.substr( 10 );
    else if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }
  SpecUtils::ireplace_all( g_data_dir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

  if( g_data_dir.empty() )
  {
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        g_data_dir = d;
        break;
      }
    }
  }//if( g_data_dir.empty() )

  if( g_test_file_dir.empty() )
  {
    for( const char * const d : { "test_data", "../test_data", "../../test_data" } )
    {
      if( SpecUtils::is_directory( SpecUtils::append_path(d, "det_eff") ) )
      {
        g_test_file_dir = d;
        break;
      }
    }
  }//if( g_test_file_dir.empty() )

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(g_data_dir, "sandia.decay.xml") ),
                         "sandia.decay.xml not in '" << g_data_dir << "'; pass --datadir=" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory( SpecUtils::append_path(g_test_file_dir, "det_eff") ),
                         "No det_eff in '" << g_test_file_dir << "'; pass --testfiledir=" );

  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path( g_data_dir, "sandia.decay.xml" ) );
}//set_dirs()


shared_ptr<const DrfImport::ParsedFile> parse_path( const string &path )
{
  return DrfImport::parseFile( SpecUtils::filename( path ), DrfImport::readFile( path ), false );
}


shared_ptr<const DrfImport::ParsedFile> parse_text( const string &name, const string &contents )
{
  return DrfImport::parseFile( name, make_shared<const string>( contents ), false );
}


bool offers( const DrfImport::Source &src, const Interpretation interp )
{
  return std::find( begin(src.interpretations), end(src.interpretations), interp )
         != end(src.interpretations);
}


void append_u16( string &out, const uint16_t v )
{
  out.push_back( static_cast<char>( v & 0xFF ) );
  out.push_back( static_cast<char>( (v >> 8) & 0xFF ) );
}

template<class T>
void append_le( string &out, const T value )
{
  static_assert( (sizeof(T) == 4) || (sizeof(T) == 8), "float or double only" );
  uint8_t bytes[sizeof(T)];
  memcpy( bytes, &value, sizeof(T) );  // all supported platforms are little-endian
  out.append( reinterpret_cast<const char *>( bytes ), sizeof(T) );
}


/** The bytes of a small .par grid, in the layout `DetEffG2kPar::parseParFile` decodes: an energy
 header, then one record per energy - a 28-byte record header (the marker double at offset 16),
 then the uint16 cells.
 */
string synthetic_par_bytes( DetEffG2kPar::ParFile &par )
{
  par.emin_keV = 60.0;
  par.emax_keV = 1332.0;
  par.energies_keV = { 60.0, 300.0, 1332.0 };
  par.grids.clear();

  string out;
  append_le<double>( out, par.emin_keV );
  append_le<double>( out, par.emax_keV );
  append_u16( out, static_cast<uint16_t>( par.energies_keV.size() ) );
  for( const double energy : par.energies_keV )
    append_le<double>( out, energy );

  for( size_t e = 0; e < par.energies_keV.size(); ++e )
  {
    DetEffG2kPar::ParGrid g;
    g.ncols = 19;                                    // 0..180 deg in 10 deg steps
    g.nrows = 40;
    g.theta_step_rad = static_cast<float>( 10.0 * 3.14159265358979323846 / 180.0 );
    g.r_step = 0.2f;                                 // ln(mm)
    for( int r = 0; r < g.nrows; ++r )
    {
      for( int c = 0; c < g.ncols; ++c )  // efficiency falls with distance, angle and energy
        g.V.push_back( static_cast<uint16_t>( 2000 + 400*e + 40*r + 15*c ) );
    }

    const size_t rec_start = out.size();
    append_u16( out, g.ncols );
    append_u16( out, g.nrows );
    out.append( 8, '\0' );
    append_le<float>( out, static_cast<float>( g.theta_step_rad ) );
    append_le<double>( out, 4707532.0 );
    append_le<float>( out, static_cast<float>( g.r_step ) );
    BOOST_REQUIRE_EQUAL( out.size() - rec_start, 28 );
    for( const uint16_t v : g.V )
      append_u16( out, v );

    par.grids.push_back( g );
  }//for( each energy )

  return out;
}//synthetic_par_bytes(...)


const char * const sm_detector_txt =
  "# 12345 - SYNTHETIC - S/N-TEST-1\r\n"
  "SynthDet,70.0,60.0,0,80.0,150.0,4.0,4.0,26,synth.par,4, #\r\n"
  "ge,0.7,5.35, #\r\n"
  "al,0.5,2.7, #\r\n"
  ",,, #\r\n"
  "ge,0.7,5.35, #\r\n"
  "al,1.0,2.7, #\r\n"
  ",,, #\r\n"
  "ge,0.5,5.35, #\r\n"
  "al,3.0,2.7, #\r\n"
  ",,\r\n";

const char * const sm_two_record_detector_txt =
  "FirstDet,70.0,60.0,0,80.0,150.0,4.0,4.0,26,first.par,4, #\r\n"
  "ge,0.7,5.35, #\r\n"
  "al,0.5,2.7, #\r\n"
  "\r\n"
  "SecondDet,50.0,40.0,0,60.0,120.0,3.0,3.0,26,second.par,4, #\r\n"
  "ge,0.7,5.35, #\r\n"
  "al,0.5,2.7, #\r\n";
}//namespace


BOOST_AUTO_TEST_CASE( IdentifiesTestFiles )
{
  set_dirs();

  const string det_eff = SpecUtils::append_path( g_test_file_dir, "det_eff" );
  const string gadras_dir = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/HPGe 40%" );

  const vector<pair<string,FileKind>> expected{
    { SpecUtils::append_path( det_eff, "Detective-X_in-situ.ecc" ), FileKind::IsocsEcc },
    { SpecUtils::append_path( det_eff, "Angle-example-efficiency.outx" ), FileKind::Angle },
    { SpecUtils::append_path( det_eff, "Angle-detector-only.detx" ), FileKind::Angle },
    { SpecUtils::append_path( det_eff, "example_extractedfromgamanal_gamEff.csv" ), FileKind::EfficiencyCsv },
    { SpecUtils::append_path( det_eff, "Run_effoutput_BIG8_Foil_86mm_110p_Abs-1.csv" ), FileKind::EfficiencyCsv },
    { SpecUtils::append_path( gadras_dir, "Efficiency.csv" ), FileKind::GadrasEfficiencyCsv },
    { SpecUtils::append_path( gadras_dir, "Detector.dat" ), FileKind::GadrasDetectorDat },
    { SpecUtils::append_path( g_test_file_dir, "gadras_detectors/NaI_3x3_text/Detector.dat" ), FileKind::GadrasDetectorDat },
    { SpecUtils::append_path( g_test_file_dir, "gadras_detectors/Detective_X_xml/Detector.dat" ), FileKind::GadrasDetectorDat },
    { SpecUtils::append_path( g_test_file_dir, "RelActAutoBatch/Eu152/other_det_eff.drf.xml" ), FileKind::DrfXml },
  };

  for( const pair<string,FileKind> &file : expected )
  {
    BOOST_TEST_CONTEXT( file.first )
    {
      shared_ptr<const DrfImport::ParsedFile> parsed;
      BOOST_REQUIRE_NO_THROW( parsed = parse_path( file.first ) );
      BOOST_CHECK( parsed->kind == file.second );
    }
  }//for( const pair<string,FileKind> &file : expected )

  // A multi-DRF list offers a choice of detector.  This one's column headings are past the part of
  //  the file an app-wide drop looks at, but the "Import" tab still accepts it.
  const string list_path = SpecUtils::append_path( g_data_dir, "common_drfs.tsv" );
  shared_ptr<const DrfImport::ParsedFile> list_file;
  BOOST_REQUIRE_NO_THROW( list_file = DrfImport::parseFile( "common_drfs.tsv",
                                                    DrfImport::readFile( list_path ), true ) );
  BOOST_CHECK( list_file->kind == FileKind::MultiDrfCsv );
  const DrfImport::Source list = DrfImport::makeSource( list_file, nullptr, 1 );
  BOOST_CHECK( list.status == Status::Ready );
  BOOST_CHECK_GT( list.recordNames.size(), 1 );
  BOOST_CHECK_EQUAL( list.record, 1 );

  // A spectrum file is not a DRF file, and would not open the import dialog.
  {
    const string n42 = SpecUtils::append_path( g_test_file_dir, "AnalystTests/Ba133_Cs137_RandomSummingOnly.n42" );
    const shared_ptr<const string> data = DrfImport::readFile( n42 );
    const size_t len = (std::min)( data->size(), size_t(1024) );
    for( const FileKind kind : DrfImport::candidateKinds( reinterpret_cast<const uint8_t *>( data->data() ),
                                                          len, data->size() ) )
      BOOST_CHECK( DrfImport::isCompleteDrfKind( kind ) );
    BOOST_CHECK_THROW( DrfImport::parseFile( "spec.n42", data, false ), std::exception );
  }

  // Oversized input is refused before any parsing.
  BOOST_CHECK_THROW( DrfImport::parseFile( "big.csv",
                        make_shared<const string>( DrfImport::sm_maxFileBytes + 1, ' ' ), true ),
                     std::exception );
}//BOOST_AUTO_TEST_CASE( IdentifiesTestFiles )


BOOST_AUTO_TEST_CASE( EccInterpretations )
{
  set_dirs();

  const string ecc_path = SpecUtils::append_path( g_test_file_dir, "det_eff/Detective-X_in-situ.ecc" );
  const DrfImport::Source src = DrfImport::makeSource( parse_path( ecc_path ), nullptr, 0 );
  BOOST_REQUIRE( src.status == Status::Ready );
  BOOST_REQUIRE( src.base );
  BOOST_CHECK( src.defaultInterpretation == Interpretation::FixedTotal );
  BOOST_CHECK( offers( src, Interpretation::FixedPerCm2 ) );   //the file states a source area
  BOOST_CHECK( offers( src, Interpretation::FixedPerGram ) );  //and a source mass
  BOOST_CHECK( !offers( src, Interpretation::GenericDetector ) );
  BOOST_CHECK_GE( src.uncertEnergies.size(), 2 );

  const uint64_t base_hash = src.base->hashValue();

  DrfImport::Options options;
  options.interpretation = Interpretation::FixedTotal;
  const DrfImport::Result total = DrfImport::build( src, options );
  BOOST_REQUIRE( total.status == Status::Ready );
  BOOST_CHECK( total.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );

  options.interpretation = Interpretation::FixedPerCm2;
  const DrfImport::Result per_cm2 = DrfImport::build( src, options );
  BOOST_REQUIRE( per_cm2.status == Status::Ready );
  BOOST_CHECK( per_cm2.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2 );

  options.interpretation = Interpretation::FarFieldIntrinsic;
  BOOST_CHECK( DrfImport::build( src, options ).status == Status::NeedsDiameter );
  options.diameter = 5.0 * PhysicalUnits::cm;
  const DrfImport::Result intrinsic = DrfImport::build( src, options );
  BOOST_REQUIRE( intrinsic.status == Status::Ready );
  BOOST_CHECK( intrinsic.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  options.interpretation = Interpretation::FarFieldAbsolute;
  BOOST_CHECK( DrfImport::build( src, options ).status == Status::NeedsDistance );
  options.distance = 25.0 * PhysicalUnits::cm;
  BOOST_CHECK( DrfImport::build( src, options ).status == Status::Ready );

  options.interpretation = Interpretation::FixedPerGram;
  const DrfImport::Result per_gram = DrfImport::build( src, options );
  BOOST_REQUIRE( per_gram.status == Status::Ready );
  BOOST_CHECK( per_gram.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram );

  // An interpretation the file does not offer falls back to its default.
  options.interpretation = Interpretation::GenericDetector;
  const DrfImport::Result fallback = DrfImport::build( src, options );
  BOOST_REQUIRE( fallback.status == Status::Ready );
  BOOST_CHECK( fallback.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );

  // Every build is a new object, and none touches the source - the GUI relies on this, since
  //  the DRFs it has already handed out are in use elsewhere.
  options.interpretation = Interpretation::FixedTotal;
  options.name = "Renamed";
  options.setUncert = true;   //with a null uncertainty: clears it
  const DrfImport::Result renamed = DrfImport::build( src, options );
  BOOST_REQUIRE( renamed.status == Status::Ready );
  BOOST_CHECK( renamed.drf != total.drf );
  BOOST_CHECK( renamed.drf.get() != src.base.get() );
  BOOST_CHECK_EQUAL( renamed.drf->name(), "Renamed" );
  BOOST_CHECK( !renamed.drf->efficiencyUncert() );
  BOOST_CHECK_EQUAL( src.base->hashValue(), base_hash );
  BOOST_CHECK( src.base->name() != "Renamed" );
}//BOOST_AUTO_TEST_CASE( EccInterpretations )


BOOST_AUTO_TEST_CASE( GadrasPair )
{
  set_dirs();

  const string dir = SpecUtils::append_path( g_test_file_dir, "gadras_detectors/Detective_X_xml" );
  const shared_ptr<const DrfImport::ParsedFile> csv = parse_path( SpecUtils::append_path( dir, "Efficiency.csv" ) );
  const shared_ptr<const DrfImport::ParsedFile> dat = parse_path( SpecUtils::append_path( dir, "Detector.dat" ) );
  BOOST_REQUIRE( DrfImport::areCompanions( csv->kind, dat->kind ) );

  // Efficiency.csv alone: intrinsic, but the diameter is in the Detector.dat
  const DrfImport::Source csv_only = DrfImport::makeSource( csv, nullptr, 0 );
  BOOST_REQUIRE( csv_only.status == Status::Ready );
  BOOST_CHECK( csv_only.defaultInterpretation == Interpretation::FarFieldIntrinsic );
  DrfImport::Options options;
  options.interpretation = Interpretation::FarFieldIntrinsic;
  BOOST_CHECK( DrfImport::build( csv_only, options ).status == Status::NeedsDiameter );

  // Together, in either order, they give the same complete detector - the way GADRAS itself reads
  //  a detector directory.
  const DrfImport::Source ab = DrfImport::makeSource( csv, dat, 0 );
  const DrfImport::Source ba = DrfImport::makeSource( dat, csv, 0 );
  BOOST_REQUIRE_MESSAGE( ab.status == Status::Ready, ab.error );
  BOOST_REQUIRE_MESSAGE( ba.status == Status::Ready, ba.error );
  BOOST_CHECK( ab.interpretations.empty() );
  BOOST_CHECK( ab.nameIsFileStem );

  const DrfImport::Result a = DrfImport::build( ab, DrfImport::Options() );
  const DrfImport::Result b = DrfImport::build( ba, DrfImport::Options() );
  BOOST_REQUIRE( a.drf && b.drf );
  BOOST_CHECK( a.drf->isValid() );
  BOOST_CHECK( a.drf->geometry() );
  BOOST_CHECK( a.drf->drfSource() == DetectorPeakResponse::DrfSource::UserImportedGadrasDrf );
  BOOST_CHECK_GT( a.drf->detectorDiameter(), 0.0f );
  BOOST_CHECK_EQUAL( a.drf->detectorDiameter(), b.drf->detectorDiameter() );
  for( const float energy : { 60.0f, 661.0f, 1332.0f } )
    BOOST_CHECK_CLOSE( a.drf->farFieldIntrinsicEfficiency( energy ), b.drf->farFieldIntrinsicEfficiency( energy ), 1.0E-4 );

  // The reference: the same pair read the way the GADRAS tab does.
  auto reference = make_shared<DetectorPeakResponse>();
  reference->fromGadrasDirectory( dir );
  for( const float energy : { 60.0f, 661.0f, 1332.0f } )
    BOOST_CHECK_CLOSE( a.drf->farFieldIntrinsicEfficiency( energy ), reference->farFieldIntrinsicEfficiency( energy ), 1.0E-4 );

  // A Detector.dat alone: a geometry to characterize, not a usable detector
  const string lone_dat = SpecUtils::append_path( g_test_file_dir, "gadras_detectors/NaI_3x3_text/Detector.dat" );
  const DrfImport::Source dat_only = DrfImport::makeSource( parse_path( lone_dat ), nullptr, 0 );
  BOOST_REQUIRE_MESSAGE( dat_only.status == Status::NeedsCharacterization, dat_only.error );
  BOOST_REQUIRE( dat_only.base );
  BOOST_CHECK( dat_only.base->geometry() );
  BOOST_CHECK( !dat_only.base->isValid() );
  BOOST_CHECK( DrfImport::build( dat_only, options ).status == Status::NeedsCharacterization );

  // ... unless it states no usable crystal (these generic ones give a zero length), which is said
  //  rather than offering a characterization that cannot run.
  const string no_geom_dat = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors/HPGe 40%/Detector.dat" );
  const DrfImport::Source no_geom = DrfImport::makeSource( parse_path( no_geom_dat ), nullptr, 0 );
  BOOST_CHECK( no_geom.status == Status::Error );
  BOOST_CHECK( !no_geom.error.empty() );
  BOOST_CHECK( !no_geom.notes.empty() );

  // Two halves of different pairs, or two of the same half, are not companions.
  BOOST_CHECK( !DrfImport::areCompanions( FileKind::GadrasEfficiencyCsv, FileKind::ParDetectorTxt ) );
  BOOST_CHECK( !DrfImport::areCompanions( FileKind::GadrasDetectorDat, FileKind::GadrasDetectorDat ) );
}//BOOST_AUTO_TEST_CASE( GadrasPair )


BOOST_AUTO_TEST_CASE( AngleFiles )
{
  set_dirs();

  const string det_eff = SpecUtils::append_path( g_test_file_dir, "det_eff" );

  // A file with the detector model and a measured reference curve defaults to a geometry-modeled
  //  detector, with its curve-transfer response already attached.
  size_t num_generic = 0;
  for( const string &path : SpecUtils::recursive_ls( det_eff, ".outx" ) )
  {
    shared_ptr<const DrfImport::ParsedFile> parsed;
    try
    {
      parsed = parse_path( path );
    }catch( std::exception & )
    {
      continue;  //some fixtures are deliberately malformed
    }

    const DrfImport::Source src = DrfImport::makeSource( parsed, nullptr, 0 );
    if( !src.generic )
      continue;

    ++num_generic;
    BOOST_TEST_CONTEXT( path )
    {
      BOOST_CHECK( src.defaultInterpretation == Interpretation::GenericDetector );

      DrfImport::Options options;
      options.interpretation = Interpretation::GenericDetector;
      const DrfImport::Result result = DrfImport::build( src, options );
      BOOST_REQUIRE( result.status == Status::Ready );
      BOOST_CHECK( result.drf->geometry() );
      BOOST_CHECK( result.drf->ceeloResponse() );
      BOOST_CHECK( !result.drf->isFixedGeometry() );
    }
  }//for( each .outx file )
  BOOST_CHECK_GE( num_generic, 1 );

  // A well detector has no geometry model, but its results are still usable, and it says why.
  {
    const DrfImport::Source well = DrfImport::makeSource(
                                parse_path( SpecUtils::append_path( det_eff, "Angle-well.outx" ) ), nullptr, 0 );
    BOOST_CHECK( well.status == Status::Ready );
    BOOST_CHECK( !offers( well, Interpretation::GenericDetector ) );
    BOOST_CHECK( !well.notes.empty() );
  }

  // Geometry only: something to characterize.
  {
    const DrfImport::Source detx = DrfImport::makeSource(
                        parse_path( SpecUtils::append_path( det_eff, "Angle-detector-only.detx" ) ), nullptr, 0 );
    BOOST_CHECK( detx.status == Status::NeedsCharacterization );
    BOOST_REQUIRE( detx.base );
    BOOST_CHECK( detx.base->geometry() );
    BOOST_CHECK_GT( detx.base->detectorDiameter(), 0.0f );
  }
}//BOOST_AUTO_TEST_CASE( AngleFiles )


BOOST_AUTO_TEST_CASE( ParAndDetectorTxt )
{
  set_dirs();

  DetEffG2kPar::ParFile expected;
  const string par_bytes = synthetic_par_bytes( expected );

  // The serializer writes what the reader reads.
  const DetEffG2kPar::ParFile decoded = DetEffG2kPar::parseParFile(
                                          vector<uint8_t>( begin(par_bytes), end(par_bytes) ) );
  BOOST_REQUIRE_EQUAL( decoded.energies_keV.size(), expected.energies_keV.size() );
  for( size_t e = 0; e < decoded.grids.size(); ++e )
    BOOST_CHECK( decoded.grids[e].V == expected.grids[e].V );

  BOOST_CHECK( DetEffG2kPar::isCandidateParFile( reinterpret_cast<const uint8_t *>( par_bytes.data() ),
                                                 (std::min)( par_bytes.size(), size_t(1024) ),
                                                 par_bytes.size() ) );
  BOOST_CHECK( DetEffG2kPar::isCandidateDetectorTxt( sm_detector_txt ) );

  // Neither looks like the other, nor like text that is not a DETECTOR.txt
  BOOST_CHECK( !DetEffG2kPar::isCandidateDetectorTxt( par_bytes.substr( 0, 1024 ) ) );
  BOOST_CHECK( !DetEffG2kPar::isCandidateParFile( reinterpret_cast<const uint8_t *>( sm_detector_txt ),
                                                  strlen( sm_detector_txt ), strlen( sm_detector_txt ) ) );
  const string csv_text = "Energy (keV), Efficiency (%)\n10, 0.0\n20,10.3\n25,28.3\n30,43.6\n";
  BOOST_CHECK( !DetEffG2kPar::isCandidateParFile( reinterpret_cast<const uint8_t *>( csv_text.data() ),
                                                  csv_text.size(), csv_text.size() ) );
  BOOST_CHECK( !DetEffG2kPar::isCandidateDetectorTxt( csv_text ) );

  const shared_ptr<const DrfImport::ParsedFile> par = parse_text( "synth.par", par_bytes );
  const shared_ptr<const DrfImport::ParsedFile> txt = parse_text( "DETECTOR.TXT", sm_detector_txt );
  BOOST_REQUIRE( par->kind == FileKind::ParGrid );
  BOOST_REQUIRE( txt->kind == FileKind::ParDetectorTxt );

  // Each needs the other
  BOOST_CHECK( DrfImport::makeSource( par, nullptr, 0 ).status == Status::NeedsCompanion );
  BOOST_CHECK( DrfImport::makeSource( txt, nullptr, 0 ).status == Status::NeedsCompanion );
  BOOST_CHECK( DrfImport::isSlow( par.get(), txt.get() ) );
  BOOST_CHECK( !DrfImport::isSlow( par.get(), nullptr ) );

  // Building from a grid is slow (tens of seconds in a Debug build), so only two builds: the pair
  //  in reverse order (the .par first is the usual one, and is covered by the record choice below),
  //  and a record chosen from a DETECTOR.txt with several.
  {
    const DrfImport::Source src = DrfImport::makeSource( txt, par, 0 );
    BOOST_REQUIRE_MESSAGE( src.status == Status::Ready, src.error );
    BOOST_CHECK( src.recordNames.empty() );
    BOOST_CHECK_EQUAL( src.name, "SynthDet" );

    const DrfImport::Result result = DrfImport::build( src, DrfImport::Options() );
    BOOST_REQUIRE( result.drf );
    BOOST_CHECK( result.drf->drfSource() == DetectorPeakResponse::DrfSource::CharacterizationParFile );
    BOOST_CHECK( result.drf->ceeloResponse() );
    BOOST_CHECK( result.drf->geometry() );

    // Within the grid's energies and at ordinary distances, nothing is out of range
    for( const float energy : { 60.0f, 300.0f, 661.0f, 1332.0f } )
    {
      for( const double dist_cm : { 5.0, 25.0, 100.0 } )
      {
        const DetectorPeakResponse::EffEval ev = result.drf->efficiencyEval( energy, dist_cm*PhysicalUnits::cm );
        BOOST_CHECK_MESSAGE( ev.flag == DetectorPeakResponse::EffFlag::Ok,
                             energy << " keV at " << dist_cm << " cm flagged "
                             << DetectorPeakResponse::effFlagName( ev.flag ) );
      }
    }

    // The float-stored 10 deg step still yields a 90 deg node at exactly cos = 0, so a source in
    //  the face plane is not "clamped".
    const shared_ptr<const ceelo::DetectorResponse> resp = result.drf->ceeloResponse();
    BOOST_REQUIRE( resp );
    BOOST_CHECK_EQUAL( resp->eta_fep.cos_thetas.size(), 10 );
    BOOST_CHECK_EQUAL( resp->eta_fep.cos_thetas.front(), 0.0 );
    const ceelo::EffResult side_on = resp->eps_fep_at( 661.0, ceelo::source_position( 30.0, 0.0, 0.0 ) );
    BOOST_CHECK( side_on.flag != ceelo::ResponseFlag::OutOfRangeClamped );
  }

  // A DETECTOR.txt with several records, none naming this .par: the user picks one.
  const shared_ptr<const DrfImport::ParsedFile> two = parse_text( "DETECTOR.TXT", sm_two_record_detector_txt );
  BOOST_REQUIRE_EQUAL( two->detectorDefs.size(), 2 );
  const DrfImport::Source pick = DrfImport::makeSource( par, two, 1 );
  BOOST_REQUIRE_MESSAGE( pick.status == Status::Ready, pick.error );
  BOOST_CHECK_EQUAL( pick.recordNames.size(), 2 );
  BOOST_CHECK_EQUAL( pick.record, 1 );
  BOOST_CHECK_EQUAL( pick.name, "SecondDet" );
  BOOST_CHECK( !pick.notes.empty() );

  // ... unless one of them names it (case-insensitively), when there is no choice to make
  BOOST_CHECK_EQUAL( DetEffG2kPar::selectDetectorDef( two->detectorDefs, "second.PAR" ).name, "SecondDet" );
  BOOST_CHECK_THROW( DetEffG2kPar::selectDetectorDef( two->detectorDefs, "synth.par" ), std::exception );
}//BOOST_AUTO_TEST_CASE( ParAndDetectorTxt )


BOOST_AUTO_TEST_CASE( ParFilesAreNotSpectra )
{
  set_dirs();

  // Dropped on the app, a file is first tried as a spectrum; neither of these may be taken for one.
  DetEffG2kPar::ParFile par;
  const vector<pair<string,string>> files{
    { "synth.par", synthetic_par_bytes( par ) },
    { "DETECTOR.txt", sm_detector_txt }
  };

  for( const pair<string,string> &file : files )
  {
    const string tmpname = SpecUtils::temp_file_name( "drf_import_" + file.first, SpecUtils::temp_dir() );
    {
      ofstream out( tmpname.c_str(), ios::out | ios::binary );
      out.write( file.second.data(), file.second.size() );
      BOOST_REQUIRE( out.good() );
    }
    BOOST_REQUIRE_EQUAL( SpecUtils::file_size( tmpname ), file.second.size() );

    SpecUtils::SpecFile spec;
    const bool loaded = spec.load_file( tmpname, SpecUtils::ParserType::Auto,
                                        SpecUtils::file_extension( file.first ) );
    SpecUtils::remove_file( tmpname );
    BOOST_CHECK_MESSAGE( !loaded, file.first << " was read as a spectrum file" );
  }//for( both files )
}//BOOST_AUTO_TEST_CASE( ParFilesAreNotSpectra )


BOOST_AUTO_TEST_CASE( MakeDrfCsv )
{
  set_dirs();

  // As written by MakeDrf::writeCsvSummary (trimmed): the coefficient line, its uncertainties, and
  //  the geometry type column.
  const string csv =
    "# Detector Response Function generated by InterSpec\n"
    "# Name,SynthMakeDrf\n"
    "# Description,A test detector\n"
    "# Intrinsic Efficiency Coefficients for equation of form Eff(x) = exp( C_0 + C_1*log(x) + C_2*log(x)^2 + ...) where x is energy in MeV\n"
    "#  i.e. equation for probability of gamma that hits the face of the detector being detected in the full energy photopeak.\n"
    "# Name,Relative Eff @ 1332 keV,eff.c name,c0,c1,c2,c3,c4,c5,c6,c7,p0,p1,p2,Calib Distance,Radius (cm),G factor,GeometryType\n"
    "SynthMakeDrf Intrinsic,50%,,-1.5,-0.6,-0.05,,,,,,,,,0.0,3.81,0.5,FixedPerGram\n"
    "# 1 sigma Uncertainties,1%,,0.01,0.02,0.003\n"
    "# Chi2 / DOF = 5 / 4 = 1.25\n"
    "Valid energy range: 50 keV to 3000 keV.\n";

  const shared_ptr<const DrfImport::ParsedFile> parsed = parse_text( "synth.csv", csv );
  BOOST_REQUIRE( parsed->kind == FileKind::MakeDrfCsv );
  BOOST_REQUIRE_EQUAL( parsed->drfs.size(), 1 );

  const shared_ptr<const DetectorPeakResponse> drf = parsed->drfs.front();
  BOOST_CHECK( drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram );
  BOOST_CHECK_CLOSE( drf->lowerEnergy(), 50.0, 1.0E-6 );
  BOOST_CHECK_CLOSE( drf->upperEnergy(), 3000.0, 1.0E-6 );

  const shared_ptr<const DetectorEfficiencyCurve> curve = drf->efficiencyCurve();
  BOOST_REQUIRE( curve );
  const vector<float> &uncerts = curve->expOfLogPowerSeriesUncerts();
  BOOST_REQUIRE_GE( uncerts.size(), 3 );
  BOOST_CHECK_CLOSE( uncerts[0], 0.01f, 1.0E-3 );
  BOOST_CHECK_CLOSE( uncerts[1], 0.02f, 1.0E-3 );
  BOOST_CHECK_CLOSE( uncerts[2], 0.003f, 1.0E-3 );

  // Older exports had a yes/no fixed-geometry column
  string old_csv = csv;
  SpecUtils::ireplace_all( old_csv, "GeometryType", "FixedGeometry" );
  SpecUtils::ireplace_all( old_csv, "FixedPerGram", "yes" );
  const shared_ptr<const DrfImport::ParsedFile> old_parsed = parse_text( "old.csv", old_csv );
  BOOST_REQUIRE( old_parsed->kind == FileKind::MakeDrfCsv );
  BOOST_CHECK( old_parsed->drfs.front()->geometryType()
               == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
}//BOOST_AUTO_TEST_CASE( MakeDrfCsv )


BOOST_AUTO_TEST_CASE( MakeDrfCsvFullExport )
{
  set_dirs();

  // A far-field detector with a geometry, to write the geometry line MakeDrf writes.
  const string gadras_dir = SpecUtils::append_path( g_test_file_dir, "gadras_detectors/Detective_X_xml" );
  auto gadras = make_shared<DetectorPeakResponse>();
  gadras->fromGadrasDirectory( gadras_dir );
  BOOST_REQUIRE( gadras->geometry() );
  string geometry_xml = gadras->geometry()->to_xml_string();
  SpecUtils::ireplace_all( geometry_xml, "\r", "" );
  SpecUtils::ireplace_all( geometry_xml, "\n", "" );

  // The layout MakeDrf::writeCsvSummary writes for a far-field DRF with an FWHM fit.  The
  //  description mentions "MeV", but the equation is in keV.
  const double c0 = -20.0, c1 = 5.0, c2 = -0.5;
  const string header =
    "# Detector Response Function generated by InterSpec 20260930T120000\r\n"
    "\r\n"
    "# Name,FullExport\r\n"
    "# Description,Good from 50 keV to 3 MeV\r\n"
    "\r\n"
    "# Intrinsic Efficiency Coefficients for equation of form Eff(x) = exp( C_0 + C_1*log(x) + C_2*log(x)^2 + ...) where x is energy in keV\r\n"
    "#  i.e. equation for probability of gamma that hits the face of the detector being detected in the full energy photopeak.\r\n"
    "# Name,Relative Eff @ 1332 keV,eff.c name,c0,c1,c2,c3,c4,c5,c6,c7,p0,p1,p2,Calib Distance,Radius (cm),G factor,GeometryType\r\n"
    "FullExport Intrinsic,50%,,-20,5,-0.5,,,,,,,,,0.0,3.81,0.5,FarField\r\n";
  const string uncerts = "# 1 sigma Uncertainties,1%,,0.1,0.02,0.003\r\n";
  const string rest =
    "# Chi2 / DOF = 5 / 4 = 1.25\r\n"
    "\r\n"
    "# Absolute Efficiency Coefficients (i.e. probability of gamma emitted from source at 25cm being detected in the full energy photopeak) for equation of form Eff(x) = exp( C_0 + C_1*log(x) + C_2*log(x)^2 + ...) where x is energy in keV and at distance of 25 cm\r\n"
    "#  i.e. equation for probability of gamma emitted from source at 25cm being detected in the full energy photopeak.\r\n"
    "# Name,Relative Eff @ 1332 keV,eff.c name,c0,c1,c2,c3,c4,c5,c6,c7,p0,p1,p2,Calib Distance,Radius (cm),G factor\r\n"
    "FullExport Absolute,50%,,-26,5,-0.5,,,,,,,,,25,3.81,0.0023\r\n"
    "# 1 sigma Uncertainties,1%,,0.1,0.02,0.003\r\n"
    "\r\n";
  const string fwhm =
    "# Full width half maximum (FWHM) follows equation: FWHM = sqrt( A0 + A1*energy + A2/energy )\r\n"
    "# Energy in keV\r\n"
    "# ,A0,A1,A2\r\n"
    "Values,1.5,0.002,3\r\n"
    "Uncertainties,0.1,0.0001,0.5\r\n"
    "# Chi2 / DOF = 3 / 2 = 1.5\r\n"
    "\r\n"
    "Detector diameter = 7.62 cm.\r\n"
    "Detector setback = 2.5 cm.\r\n"
    "# Detector geometry (CeeLo XML): " + geometry_xml + "\r\n"
    "Valid energy range: 50 keV to 3000 keV.\r\n"
    "\r\n"
    "# Peaks used to create DRF\r\n";

  const auto check_drf = [=]( const shared_ptr<const DetectorPeakResponse> &drf, const bool has_uncerts ){
    BOOST_REQUIRE( drf );
    BOOST_CHECK( drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    BOOST_CHECK_CLOSE( drf->detectorDiameter(), 7.62*PhysicalUnits::cm, 1.0E-3 );

    // The equation is in keV, whatever the description says
    const double lnE = std::log( 661.0 );
    const double expected = std::exp( c0 + c1*lnE + c2*lnE*lnE );
    BOOST_CHECK_CLOSE( drf->farFieldIntrinsicEfficiency( 661.0f ), expected, 1.0E-3 );

    // Everything after the FWHM block used to be skipped
    BOOST_CHECK_CLOSE( drf->lowerEnergy(), 50.0, 1.0E-6 );
    BOOST_CHECK_CLOSE( drf->upperEnergy(), 3000.0, 1.0E-6 );
    BOOST_CHECK_CLOSE( drf->detectorSetback(), 2.5*PhysicalUnits::cm, 1.0E-3 );
    BOOST_CHECK( drf->geometry() );

    BOOST_CHECK( drf->resolutionFcnType() == DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse );
    BOOST_REQUIRE_EQUAL( drf->resolutionFcnCoefficients().size(), 3 );
    BOOST_CHECK_CLOSE( drf->resolutionFcnCoefficients()[2], 3.0f, 1.0E-4 );
    BOOST_REQUIRE_EQUAL( drf->resolutionFcnUncertainties().size(), 3 );
    BOOST_CHECK_CLOSE( drf->resolutionFcnUncertainties()[2], 0.5f, 1.0E-4 );

    BOOST_CHECK_EQUAL( drf->efficiencyCurve()->expOfLogPowerSeriesUncerts().size(), (has_uncerts ? 3 : 0) );
  };//check_drf

  {
    std::istringstream strm( header + uncerts + rest + fwhm );
    check_drf( DetectorPeakResponse::parseInterSpecRelEffCsv( strm ), true );
  }

  // Without the uncertainties line, the line after the coefficients is still looked at - here,
  //  the start of the FWHM block.
  {
    std::istringstream strm( header + fwhm );
    check_drf( DetectorPeakResponse::parseInterSpecRelEffCsv( strm ), false );
  }

  // An unusable intrinsic table must not fall through to the absolute one after it
  {
    string bad = header + uncerts + rest + fwhm;
    SpecUtils::ireplace_all( bad, "FullExport Intrinsic,50%,,-20,", "FullExport Intrinsic,50%,,abc," );
    std::istringstream strm( bad );
    BOOST_CHECK( !DetectorPeakResponse::parseInterSpecRelEffCsv( strm ) );
  }
}//BOOST_AUTO_TEST_CASE( MakeDrfCsvFullExport )


BOOST_AUTO_TEST_CASE( ParFileHeaderGuards )
{
  set_dirs();

  // A record marker at the very start used to wrap the framing arithmetic around, and read before
  //  the buffer; every placement near the start must be rejected cleanly.
  const double marker = 4707532.0;
  for( size_t offset = 0; offset <= 40; ++offset )
  {
    vector<uint8_t> bytes( 64, 0 );
    memcpy( &bytes[offset], &marker, sizeof(marker) );
    BOOST_CHECK_THROW( DetEffG2kPar::parseParFile( bytes ), std::exception );
  }

  // ... and a header whose energy count does not fit before the first record
  DetEffG2kPar::ParFile par;
  string good = synthetic_par_bytes( par );
  good[0x10] = static_cast<char>( 200 );  //claims 200 energies
  BOOST_CHECK_THROW( DetEffG2kPar::parseParFile( vector<uint8_t>( begin(good), end(good) ) ), std::exception );
}//BOOST_AUTO_TEST_CASE( ParFileHeaderGuards )


BOOST_AUTO_TEST_CASE( MoreImportChoices )
{
  set_dirs();

  const string det_eff = SpecUtils::append_path( g_test_file_dir, "det_eff" );
  const shared_ptr<const DrfImport::ParsedFile> gameff
                  = parse_path( SpecUtils::append_path( det_eff, "example_extractedfromgamanal_gamEff.csv" ) );
  BOOST_REQUIRE( gameff->kind == FileKind::EfficiencyCsv );

  // A plain efficiency CSV may already be per unit area or mass: relabelled, values unchanged
  const DrfImport::Source csv_only = DrfImport::makeSource( gameff, nullptr, 0 );
  BOOST_REQUIRE( csv_only.status == Status::Ready );
  BOOST_CHECK( offers( csv_only, Interpretation::FixedPerCm2 ) );
  DrfImport::Options options;
  options.interpretation = Interpretation::FixedTotal;
  const DrfImport::Result total = DrfImport::build( csv_only, options );
  options.interpretation = Interpretation::FixedPerCm2;
  const DrfImport::Result per_cm2 = DrfImport::build( csv_only, options );
  BOOST_REQUIRE( total.drf && per_cm2.drf );
  BOOST_CHECK( per_cm2.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2 );
  BOOST_CHECK_CLOSE( per_cm2.drf->efficiencyCurve()->efficiency( 661.0f ),
                     total.drf->efficiencyCurve()->efficiency( 661.0f ), 1.0E-4 );

  // With no name given, the suggested one is used (a copy is re-hashed on the way)
  BOOST_CHECK_EQUAL( total.drf->name(), csv_only.name );

  // A GADRAS Detector.dat pairs with any efficiency CSV, as intrinsic efficiencies
  const shared_ptr<const DrfImport::ParsedFile> dat
      = parse_path( SpecUtils::append_path( g_test_file_dir, "gadras_detectors/NaI_3x3_text/Detector.dat" ) );
  BOOST_REQUIRE( DrfImport::areCompanions( gameff->kind, dat->kind ) );
  const DrfImport::Source paired = DrfImport::makeSource( dat, gameff, 0 );
  BOOST_REQUIRE_MESSAGE( paired.status == Status::Ready, paired.error );
  BOOST_CHECK( !paired.notes.empty() );
  const DrfImport::Result paired_drf = DrfImport::build( paired, DrfImport::Options() );
  BOOST_REQUIRE( paired_drf.drf );
  BOOST_CHECK( paired_drf.drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  BOOST_CHECK( paired_drf.drf->geometry() );
  BOOST_CHECK_GT( paired_drf.drf->detectorDiameter(), 0.0f );

  // An ANGLE detector model with no reference curve says why it isn't offered as a generic detector
  const DrfImport::Source no_ref = DrfImport::makeSource(
                        parse_path( SpecUtils::append_path( det_eff, "Angle-no-refcurve.outx" ) ), nullptr, 0 );
  BOOST_CHECK( !offers( no_ref, Interpretation::GenericDetector ) );
  bool explained = false;
  for( const string &note : no_ref.notes )
    explained |= SpecUtils::icontains( note, "reference efficiency curve" );
  BOOST_CHECK( explained );

  // A UTF-8 BOM, or spaces and quotes around CSV column names, don't hide a file's kind
  {
    const shared_ptr<const string> dat_bytes = DrfImport::readFile(
                SpecUtils::append_path( g_test_file_dir, "gadras_detectors/NaI_3x3_text/Detector.dat" ) );
    const shared_ptr<const DrfImport::ParsedFile> with_bom
               = DrfImport::parseFile( "Detector.dat", make_shared<const string>( "\xEF\xBB\xBF" + *dat_bytes ), false );
    BOOST_CHECK( with_bom->kind == FileKind::GadrasDetectorDat );

    const string spaced = "\"Energy (keV)\", \"Efficiency (%)\"\n60, 10.5\n122, 20.1\n662, 5.2\n1332, 2.1\n";
    const vector<FileKind> kinds = DrfImport::candidateKinds( reinterpret_cast<const uint8_t *>( spaced.data() ),
                                                              spaced.size(), spaced.size() );
    BOOST_CHECK( std::find( begin(kinds), end(kinds), FileKind::EfficiencyCsv ) != end(kinds) );
  }
}//BOOST_AUTO_TEST_CASE( MoreImportChoices )
