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
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE test_GadrasDetectorDat_suite
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GadrasDetectorDat.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

// Data directory globals (set from command line).
std::string g_data_dir = "";
std::string g_test_data_dir = "";

namespace
{
  bool close_enough( const double a, const double b, const double abs_tol = 1.0e-4 )
  {
    return std::fabs(a - b) <= abs_tol;
  }

  // Directory holding the copied GADRAS example detectors.
  string gadras_dir()
  {
    return SpecUtils::append_path( g_test_data_dir, "gadras_detectors" );
  }

  string dat_path( const string &det_name )
  {
    return SpecUtils::append_path( SpecUtils::append_path( gadras_dir(), det_name ), "Detector.dat" );
  }

  GadrasDetectorDat load( const string &det_name )
  {
    const string path = dat_path( det_name );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(path), "Missing fixture: " + path );
    return GadrasDetectorDat::fromFile( path );
  }
}//namespace


BOOST_AUTO_TEST_CASE( test_parse_text_nai )
{
  const GadrasDetectorDat dat = load( "NaI_3x3_text" );

  BOOST_CHECK( dat.m_format == GadrasDetectorDat::SourceFormat::LegacyText );
  BOOST_CHECK( close_enough( dat.length(), 7.62 ) );
  BOOST_CHECK( close_enough( dat.width(), 6.66 ) );
  BOOST_CHECK( close_enough( dat.heightToWidth(), 1.0 ) );
  BOOST_CHECK( close_enough( dat.shapeFactor(), 100.0 ) );
  BOOST_CHECK( close_enough( dat.resOffset(), -6.0 ) );
  BOOST_CHECK( close_enough( dat.resFWHM661(), 7.3 ) );
  BOOST_CHECK( close_enough( dat.resPower(), 0.55 ) );
  BOOST_CHECK( close_enough( dat.solidAnglePercent(), 0.03505 ) );
  BOOST_CHECK( close_enough( dat.distanceCm(), 100.0 ) );
  BOOST_CHECK( close_enough( dat.effScalar(), 1.0 ) );
  BOOST_CHECK( close_enough( dat.lldKeV(), 27.0 ) );
  BOOST_CHECK( close_enough( dat.setbackCm(), 0.30 ) );
  BOOST_CHECK_EQUAL( dat.materialIndex(), 5 );
  BOOST_CHECK_EQUAL( dat.materialName(), string("NaI") );
  BOOST_CHECK_EQUAL( dat.m_response_version, string("18.7.4") );
  BOOST_CHECK_EQUAL( dat.m_highest_param, 80 );
}


BOOST_AUTO_TEST_CASE( test_parse_text_czt )
{
  const GadrasDetectorDat dat = load( "CZT_1cm_text" );

  BOOST_CHECK( dat.m_format == GadrasDetectorDat::SourceFormat::LegacyText );
  BOOST_CHECK_EQUAL( dat.materialIndex(), 8 );
  BOOST_CHECK_EQUAL( dat.materialName(), string("CZT") );
  BOOST_CHECK( close_enough( dat.length(), 1.07268 ) );
  BOOST_CHECK( close_enough( dat.width(), 0.94284 ) );
  BOOST_CHECK( close_enough( dat.shapeFactor(), 10.0 ) );
  BOOST_CHECK( close_enough( dat.lldKeV(), 25.0 ) );
  BOOST_CHECK( close_enough( dat.deadLayerMm(), 0.05029 ) );
}


BOOST_AUTO_TEST_CASE( test_parse_text_hpge_shields )
{
  const GadrasDetectorDat dat = load( "HPGe_Planar50_text" );

  BOOST_CHECK_EQUAL( dat.materialIndex(), 3 );
  BOOST_CHECK_EQUAL( dat.materialName(), string("HPGe") );
  BOOST_CHECK( dat.width() > 0.0f );
  BOOST_CHECK( close_enough( dat.length(), 2.97 ) );
  BOOST_CHECK( close_enough( dat.setbackCm(), 2.0 ) );

  const GadrasDetectorDat::Attenuator inner = dat.innerAttenuator();
  BOOST_CHECK( close_enough( inner.atomicNumber, 5.0 ) );
  BOOST_CHECK( close_enough( inner.arealDensity, 0.5 ) );
  BOOST_CHECK( close_enough( inner.porosityPercent, 0.0 ) );

  const GadrasDetectorDat::Attenuator outer = dat.outerAttenuator();
  BOOST_CHECK( close_enough( outer.atomicNumber, 4.0 ) );
}


BOOST_AUTO_TEST_CASE( test_parse_xml_identifinder_ngh )
{
  // Real-XML fields the previous InterSpec reader missed (setback nested in <dimensions>, and
  //  the LLD/energy-range blocks) must now parse.
  const GadrasDetectorDat dat = load( "IdentiFINDER_NGH_xml" );

  BOOST_CHECK( dat.m_format == GadrasDetectorDat::SourceFormat::Xml );
  BOOST_CHECK_EQUAL( dat.materialName(), string("NaI") );
  BOOST_CHECK_EQUAL( dat.materialIndex(), 5 );
  BOOST_CHECK( close_enough( dat.setbackCm(), 2.0 ) );
  BOOST_CHECK( close_enough( dat.length(), 5.1 ) );
  BOOST_CHECK( close_enough( dat.width(), 3.05 ) );
  BOOST_CHECK( close_enough( dat.shapeFactor(), 100.0 ) );
  BOOST_CHECK( close_enough( dat.resOffset(), -6.5 ) );
  BOOST_CHECK( close_enough( dat.resFWHM661(), 7.5 ) );
  BOOST_CHECK( close_enough( dat.lldKeV(), 14.0 ) );
  BOOST_CHECK( close_enough( dat.lowerEnergyKeV(), 25.0 ) );
  BOOST_CHECK( close_enough( dat.upperEnergyKeV(), 3000.0 ) );

  const GadrasDetectorDat::Attenuator inner = dat.innerAttenuator();
  BOOST_CHECK( close_enough( inner.atomicNumber, 13.6 ) );
  BOOST_CHECK( close_enough( inner.arealDensity, 1.09 ) );

  BOOST_CHECK_EQUAL( dat.m_file_version, string("1.2.0") );
  BOOST_CHECK_EQUAL( dat.m_response_version, string("19.2.3") );
}


BOOST_AUTO_TEST_CASE( test_parse_xml_czt_scientific )
{
  // Values in this file are in scientific notation (e.g. 8.060000e-01).
  const GadrasDetectorDat dat = load( "MikesCzt_xml" );

  BOOST_CHECK( dat.m_format == GadrasDetectorDat::SourceFormat::Xml );
  BOOST_CHECK_EQUAL( dat.materialName(), string("CZT") );
  BOOST_CHECK( close_enough( dat.setbackCm(), 0.806 ) );
  BOOST_CHECK( close_enough( dat.length(), 0.8 ) );
  BOOST_CHECK( close_enough( dat.width(), 1.95 ) );
  BOOST_CHECK( close_enough( dat.shapeFactor(), 100.0 ) );
  BOOST_CHECK( close_enough( dat.effScalar(), 4.0 ) );
  BOOST_CHECK( close_enough( dat.deadLayerMm(), 0.05065652, 1.0e-6 ) );
}


BOOST_AUTO_TEST_CASE( test_equivalent_diameter )
{
  const GadrasDetectorDat dat = load( "NaI_3x3_text" );
  const float w = dat.width(), h2w = dat.heightToWidth();
  const float expected = 2.0f * std::sqrt( w*w*h2w / 3.14159265359f );
  BOOST_CHECK( close_enough( dat.equivalentCircularDiameterCm(), expected ) );
}


BOOST_AUTO_TEST_CASE( test_infer_shape )
{
  // Scintillator, square cross-section -> cylinder.
  const GadrasDetectorDat nai = load( "NaI_3x3_text" );
  GadrasDetectorDat::InferredShape nai_shape = nai.inferShape();
  BOOST_CHECK( nai_shape.shape == GadrasDetectorDat::Shape::Cylinder );
  BOOST_CHECK( close_enough( nai_shape.dimA, 0.5*nai.width() ) );
  BOOST_CHECK( close_enough( nai_shape.dimB, nai.length() ) );

  // HPGe -> coaxial cylinder.
  const GadrasDetectorDat hpge = load( "Detective_X_xml" );
  BOOST_CHECK( hpge.inferShape().shape == GadrasDetectorDat::Shape::CoaxialCylinder );

  // CZT -> box.
  const GadrasDetectorDat czt = load( "CZT_1cm_text" );
  BOOST_CHECK( czt.inferShape().shape == GadrasDetectorDat::Shape::Box );

  // Fabricated non-unity height/width -> rectangular.
  GadrasDetectorDat rect = load( "NaI_3x3_text" );
  rect.setValue( 12, 4.0f );
  GadrasDetectorDat::InferredShape rect_shape = rect.inferShape();
  BOOST_CHECK( rect_shape.shape == GadrasDetectorDat::Shape::Rectangular );
  BOOST_CHECK( close_enough( rect_shape.dimA, rect.length() ) );
  BOOST_CHECK( close_enough( rect_shape.dimB, rect.width() ) );
  BOOST_CHECK( close_enough( rect_shape.dimC, 4.0*rect.width() ) );
}


// Round-trip a parsed file through toXml and back, asserting the named fields are preserved.
void check_round_trip( const string &det_name )
{
  const GadrasDetectorDat a = load( det_name );

  stringstream xml;
  a.toXml( xml );

  const GadrasDetectorDat b = GadrasDetectorDat::fromStream( xml );

  BOOST_CHECK_MESSAGE( close_enough( a.length(), b.length() ), det_name + ": length" );
  BOOST_CHECK_MESSAGE( close_enough( a.width(), b.width() ), det_name + ": width" );
  BOOST_CHECK_MESSAGE( close_enough( a.heightToWidth(), b.heightToWidth() ), det_name + ": h2w" );
  BOOST_CHECK_MESSAGE( close_enough( a.shapeFactor(), b.shapeFactor() ), det_name + ": shapeFactor" );
  BOOST_CHECK_MESSAGE( close_enough( a.resOffset(), b.resOffset() ), det_name + ": resOffset" );
  BOOST_CHECK_MESSAGE( close_enough( a.resFWHM661(), b.resFWHM661() ), det_name + ": resFWHM661" );
  BOOST_CHECK_MESSAGE( close_enough( a.resPower(), b.resPower() ), det_name + ": resPower" );
  BOOST_CHECK_MESSAGE( close_enough( a.setbackCm(), b.setbackCm() ), det_name + ": setback" );
  BOOST_CHECK_MESSAGE( close_enough( a.lldKeV(), b.lldKeV() ), det_name + ": lld" );
  BOOST_CHECK_MESSAGE( close_enough( a.effScalar(), b.effScalar() ), det_name + ": scalar" );
  BOOST_CHECK_MESSAGE( close_enough( a.deadLayerMm(), b.deadLayerMm(), 1.0e-4 ), det_name + ": deadLayer" );
  BOOST_CHECK_MESSAGE( close_enough( a.lowerEnergyKeV(), b.lowerEnergyKeV() ), det_name + ": lowerE" );
  BOOST_CHECK_MESSAGE( close_enough( a.upperEnergyKeV(), b.upperEnergyKeV() ), det_name + ": upperE" );
  BOOST_CHECK_MESSAGE( a.materialName() == b.materialName(), det_name + ": material" );
}


BOOST_AUTO_TEST_CASE( test_round_trip_xml )
{
  check_round_trip( "NaI_3x3_text" );        // text -> xml -> struct
  check_round_trip( "CZT_1cm_text" );
  check_round_trip( "HPGe_Planar50_text" );
  check_round_trip( "IdentiFINDER_NGH_xml" );  // xml -> xml -> struct
  check_round_trip( "MikesCzt_xml" );
}


BOOST_AUTO_TEST_CASE( test_material_table )
{
  BOOST_CHECK_EQUAL( GadrasDetectorDat::materialTable().size(), 36u );

  const GadrasDetectorDat::MaterialInfo * const hpge = GadrasDetectorDat::materialByName( "HPGe" );
  BOOST_REQUIRE( hpge );
  BOOST_CHECK_EQUAL( hpge->index, 3 );
  BOOST_CHECK( close_enough( hpge->density, 5.33, 1.0e-3 ) );
  BOOST_CHECK( hpge->specialCased );

  // Case-insensitive lookup.
  const GadrasDetectorDat::MaterialInfo * const czt = GadrasDetectorDat::materialByName( "czt" );
  BOOST_REQUIRE( czt );
  BOOST_CHECK_EQUAL( czt->index, 8 );
  BOOST_CHECK( czt->specialCased );

  const GadrasDetectorDat::MaterialInfo * const nai = GadrasDetectorDat::materialByName( "NaI" );
  BOOST_REQUIRE( nai );
  BOOST_CHECK( nai->specialCased );

  const GadrasDetectorDat::MaterialInfo * const yap = GadrasDetectorDat::materialByName( "YAP" );
  BOOST_REQUIRE( yap );
  BOOST_CHECK( !yap->specialCased );

  BOOST_CHECK( GadrasDetectorDat::materialByIndex( 0 ) == nullptr );
  BOOST_CHECK( GadrasDetectorDat::materialByIndex( 5 ) != nullptr );
  BOOST_CHECK( GadrasDetectorDat::materialByName( "not-a-material" ) == nullptr );
}


BOOST_AUTO_TEST_CASE( test_integrated_from_directory )
{
  // Regression: the refactored DetectorPeakResponse must find the setback for a real XML file
  //  (nested in <dimensions>), which the previous reader missed.
  const string det_dir = SpecUtils::append_path( gadras_dir(), "IdentiFINDER_NGH_xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(det_dir, "Efficiency.csv") ),
                         "Missing Efficiency.csv for integrated test" );

  DetectorPeakResponse drf;
  BOOST_REQUIRE_NO_THROW( drf.fromGadrasDirectory( det_dir ) );

  BOOST_CHECK( close_enough( drf.detectorSetback() / PhysicalUnits::cm, 2.0, 0.01 ) );

  const float expected_diam = 2.0f * std::sqrt( 3.05f*3.05f*1.0f / 3.14159265359f );
  BOOST_CHECK( close_enough( drf.detectorDiameter() / PhysicalUnits::cm, expected_diam, 0.01 ) );
}


// Parse the command-line data directories (mirrors other InterSpec Boost tests).
struct GlobalFixture
{
  GlobalFixture()
  {
    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( arg.find("--datadir=") == 0 )
        g_data_dir = arg.substr( 10 );
      else if( arg.find("--testfiledir=") == 0 )
        g_test_data_dir = arg.substr( 14 );
    }

    cout << "  Data directory: " << (g_data_dir.empty() ? "(not set)" : g_data_dir) << endl;
    cout << "  Test data directory: " << (g_test_data_dir.empty() ? "(not set)" : g_test_data_dir) << endl;

    if( g_test_data_dir.empty() )
      cerr << "Warning: test data directory not set (use --testfiledir=...)" << endl;

    if( !g_data_dir.empty() )
      InterSpec::setStaticDataDirectory( g_data_dir );
  }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture );
