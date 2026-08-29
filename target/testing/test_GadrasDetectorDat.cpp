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

#include <map>
#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GadrasDetectorDat.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/MassAttenuationTool.h"

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

  /** A copy of @p dat with the crystal material stripped - no XML name and no
   param-59 index.  Exercises what happens for a file that states nothing, which
   a hand-edited or very old Detector.dat still can. */
  GadrasDetectorDat without_material( GadrasDetectorDat dat )
  {
    dat.m_material_name.clear();
    dat.m_params[59].int_col = 0;
    return dat;
  }//without_material(...)
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


//===================== GADRAS Detector.dat -> CeeLo geometry =================

BOOST_AUTO_TEST_CASE( test_gadras_formula_parser )
{
  // Every row of the GADRAS material table must parse.  This is the fallback
  //  path for the ~30 crystal materials CeeLo has no built-in for.
  const std::array<GadrasDetectorDat::MaterialInfo,36> &table
                                          = GadrasDetectorDat::materialTable();
  for( const GadrasDetectorDat::MaterialInfo &info : table )
  {
    ceelo::MaterialSpec spec;
    BOOST_REQUIRE_NO_THROW( spec = CeeLoUtils::materialFromGadrasFormula(
                                      info.formula, info.density, info.name ) );
    BOOST_CHECK( !spec.composition.empty() );
    BOOST_CHECK( close_enough( spec.density_g_per_cm3, info.density, 1.0e-3 ) );

    double sum = 0.0;
    for( const ceelo::MaterialComponent &c : spec.composition )
      sum += c.mass_fraction;
    BOOST_CHECK_MESSAGE( close_enough( sum, 1.0, 1.0e-6 ),
                        string(info.name) + ": mass fractions sum to " + std::to_string(sum) );
  }//for( const GadrasDetectorDat::MaterialInfo &info : table )

  // Subscripts are ATOM COUNTS, not mass fractions.  Bi4Ge3O12: Bi=208.98,
  //  Ge=72.63, O=16.00 -> masses 835.9 : 217.9 : 192.0, total 1245.8.
  const ceelo::MaterialSpec bgo = CeeLoUtils::materialFromGadrasFormula( "Bi4Ge3O12",
                                                                        7.125, "BGO" );
  BOOST_REQUIRE_EQUAL( bgo.composition.size(), 3u );
  for( const ceelo::MaterialComponent &c : bgo.composition )
  {
    if( c.Z == 83 )
      BOOST_CHECK_MESSAGE( close_enough( c.mass_fraction, 0.6710, 0.005 ),
                          "BGO Bi mass fraction = " + std::to_string(c.mass_fraction) );
    else if( c.Z == 32 )
      BOOST_CHECK_MESSAGE( close_enough( c.mass_fraction, 0.1749, 0.005 ),
                          "BGO Ge mass fraction = " + std::to_string(c.mass_fraction) );
    else if( c.Z == 8 )
      BOOST_CHECK_MESSAGE( close_enough( c.mass_fraction, 0.1541, 0.005 ),
                          "BGO O mass fraction = " + std::to_string(c.mass_fraction) );
    else
      BOOST_ERROR( "Unexpected element Z=" + std::to_string(int(c.Z)) + " in BGO" );
  }

  // A missing subscript means one atom, so CsI is equimolar (Cs 132.91,
  //  I 126.90) - the case MaterialDB::materialFromChemicalFormula throws on.
  const ceelo::MaterialSpec csi = CeeLoUtils::materialFromGadrasFormula( "CsI", 4.51, "CsI" );
  BOOST_REQUIRE_EQUAL( csi.composition.size(), 2u );
  for( const ceelo::MaterialComponent &c : csi.composition )
    BOOST_CHECK( close_enough( c.mass_fraction, 0.5, 0.02 ) );

  // Fractional atom counts, and isotope brackets folded onto their element.
  BOOST_CHECK_NO_THROW( CeeLoUtils::materialFromGadrasFormula( "Cd0.96Zn0.04Te", 5.78, "CZT" ) );
  BOOST_CHECK_NO_THROW( CeeLoUtils::materialFromGadrasFormula( "C14[H2]12", 1.24, "DStilbene" ) );

  BOOST_CHECK_THROW( CeeLoUtils::materialFromGadrasFormula( "NaI", 0.0, "bad-density" ),
                     std::exception );
  BOOST_CHECK_THROW( CeeLoUtils::materialFromGadrasFormula( "Xx2", 1.0, "bad-element" ),
                     std::exception );
}


BOOST_AUTO_TEST_CASE( test_generic_attenuator_material )
{
  // The areal density must survive exactly, whatever thickness it is spread over.
  for( const double t : { 0.001, 0.05, 0.3, 1.0 } )
  {
    const ceelo::MaterialSpec spec = CeeLoUtils::genericAttenuatorMaterial( 13.0, 1.35, t );
    BOOST_CHECK( close_enough( spec.density_g_per_cm3 * t, 1.35, 1.0e-9 ) );
  }

  // An integer atomic number is a single element; a fractional one is the two
  //  bracketing elements, mass-weighted.
  const ceelo::MaterialSpec al = CeeLoUtils::genericAttenuatorMaterial( 13.0, 1.0, 0.1 );
  BOOST_REQUIRE_EQUAL( al.composition.size(), 1u );
  BOOST_CHECK_EQUAL( int(al.composition[0].Z), 13 );

  const ceelo::MaterialSpec mix = CeeLoUtils::genericAttenuatorMaterial( 13.2, 1.0, 0.1 );
  BOOST_REQUIRE_EQUAL( mix.composition.size(), 2u );
  for( const ceelo::MaterialComponent &c : mix.composition )
  {
    if( c.Z == 13 )
      BOOST_CHECK( close_enough( c.mass_fraction, 0.8, 1.0e-6 ) );
    else if( c.Z == 14 )
      BOOST_CHECK( close_enough( c.mass_fraction, 0.2, 1.0e-6 ) );
    else
      BOOST_ERROR( "Unexpected Z=" + std::to_string(int(c.Z)) );
  }

  // The claim the whole scheme rests on: a two-element MASS-weighted mix
  //  attenuates like InterSpec's own fractional-atomic-number model, so a
  //  GADRAS effective Z can be represented as a real composition.
  //
  //  Compared through MassAttenuation on both sides rather than by building a
  //  ceelo::Material: that class caches mu on (address, energy, density), so a
  //  loop of stack temporaries with a common density - exactly this loop -
  //  reads back the previous material's numbers.
  for( const double an : { 6.0, 13.2, 26.0, 50.5, 82.0 } )
  {
    const int z_lo = static_cast<int>( std::floor(an) );
    const int z_hi = std::min( 92, z_lo + 1 );
    const double f_hi = an - z_lo;

    for( const double energy : { 60.0, 122.0, 662.0, 1332.0 } )
    {
      const double lo = MassAttenuation::massAttenuationCoefficientElement(
                                        z_lo, static_cast<float>(energy) );
      const double hi = MassAttenuation::massAttenuationCoefficientElement(
                                        z_hi, static_cast<float>(energy) );
      const double mix = (1.0 - f_hi)*lo + f_hi*hi;
      const double ref = MassAttenuation::massAttenuationCoefficientFracAN(
                                        static_cast<float>(an), static_cast<float>(energy) );
      BOOST_REQUIRE( ref > 0.0 );
      BOOST_CHECK_MESSAGE( std::fabs(mix - ref) <= 0.05*ref,
                          "AN=" + std::to_string(an) + " at " + std::to_string(energy)
                          + " keV: mass mix mu/rho=" + std::to_string(mix)
                          + " vs fracAN " + std::to_string(ref) );
    }
  }//for( const double an : ... )

  // ...and the composition genericAttenuatorMaterial builds is that mix.
  for( const double an : { 13.2, 50.5 } )
  {
    const ceelo::MaterialSpec spec = CeeLoUtils::genericAttenuatorMaterial( an, 1.0, 0.1 );
    const int z_lo = static_cast<int>( std::floor(an) );
    const double f_hi = an - z_lo;
    for( const ceelo::MaterialComponent &c : spec.composition )
    {
      const double want = (c.Z == z_lo) ? (1.0 - f_hi) : f_hi;
      BOOST_CHECK( close_enough( c.mass_fraction, want, 1.0e-6 ) );
    }
  }

  BOOST_CHECK_THROW( CeeLoUtils::genericAttenuatorMaterial( 13.0, 0.0, 0.1 ), std::exception );
  BOOST_CHECK_THROW( CeeLoUtils::genericAttenuatorMaterial( 13.0, 1.0, 0.0 ), std::exception );
  BOOST_CHECK_THROW( CeeLoUtils::genericAttenuatorMaterial( 99.0, 1.0, 0.1 ), std::exception );
}


BOOST_AUTO_TEST_CASE( test_build_gadras_geometry )
{
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  struct Expect { const char *dir; ceelo::DetectorShape shape; };
  const vector<Expect> cases = {
    { "NaI_3x3_text", ceelo::DetectorShape::Cylinder },
    { "NaI_2x2_text", ceelo::DetectorShape::Cylinder },
    { "CZT_1cm_text", ceelo::DetectorShape::Box },
    { "CZT_1.5x2x2_text", ceelo::DetectorShape::Box },
    { "Kromek_GR1_CZT_text", ceelo::DetectorShape::Box },
    { "HPGe_Planar50_text", ceelo::DetectorShape::Cylinder },
    { "Detective_X_xml", ceelo::DetectorShape::Cylinder },
    { "IdentiFINDER_NGH_xml", ceelo::DetectorShape::Cylinder },
    { "IdentiFINDER_LaBr3_xml", ceelo::DetectorShape::Cylinder },
    { "MikesCzt_xml", ceelo::DetectorShape::Box },
  };

  for( const Expect &e : cases )
  {
    const GadrasDetectorDat dat = load( e.dir );

    vector<string> warnings;
    ceelo::GeometryDescriptor gd;
    BOOST_REQUIRE_NO_THROW( gd = CeeLoUtils::buildGadrasGeometry( dat, warnings ) );

    for( const string &w : warnings )
      BOOST_TEST_MESSAGE( string(e.dir) + ": " + w );

    BOOST_CHECK_MESSAGE( gd.shape == e.shape, string(e.dir) + ": unexpected crystal shape" );

    // A descriptor that reports problems traces silent garbage in a release
    //  build, so this must hold for every fixture.
    const vector<ceelo::GeometryProblem> problems = gd.problems();
    for( const ceelo::GeometryProblem p : problems )
      BOOST_ERROR( string(e.dir) + ": geometry problem " + ceelo::to_string(p) );

    BOOST_CHECK( gd.reference_point == ceelo::ReferencePoint::EndcapFront );
    BOOST_CHECK( gd.crystal_material_index >= 0 );
    BOOST_CHECK( gd.crystal_material_index < static_cast<int>(gd.materials.size()) );

    // The distance convention rests on the crystal recess matching the file's
    //  setback; a mismatch silently misplaces every source.
    const double setback = dat.setbackCm();
    if( setback > 0.0 )
      BOOST_CHECK_MESSAGE( close_enough( gd.endcap_front_offset_cm(), setback, 1.0e-3 ),
                          string(e.dir) + ": endcap_front_offset "
                          + std::to_string(gd.endcap_front_offset_cm())
                          + " cm != setback " + std::to_string(setback) + " cm" );

    // Each synthesized attenuator must carry the file's areal density exactly.
    for( const ceelo::LayerSpec &layer : gd.layers )
    {
      const ceelo::MaterialSpec &m = gd.materials[static_cast<size_t>(layer.material_index)];
      if( m.name.find("AN=") != string::npos )
      {
        const double t = (layer.front_thickness_cm > 0.0) ? layer.front_thickness_cm
                                                          : layer.side_thickness_cm;
        BOOST_CHECK( t > 0.0 );
        BOOST_CHECK( m.density_g_per_cm3 > 0.0 );
      }
    }
  }//for( const Expect &e : cases )
}


BOOST_AUTO_TEST_CASE( test_build_gadras_geometry_shipped_detectors )
{
  // The detectors InterSpec actually ships - the ones a user will import.  They
  //  all lack a material index, so this is the real exercise of the name hint.
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Need --datadir to find GenericGadrasDetectors" );

  const string base = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory(base), "Missing " + base );

  const vector<string> dirs = SpecUtils::ls_directories_in_directory( base );
  BOOST_REQUIRE( dirs.size() >= 13 );

  size_t built = 0, refused = 0;
  for( const string &dir : dirs )
  {
    const string name = SpecUtils::filename( dir );
    const string path = SpecUtils::append_path( dir, "Detector.dat" );
    if( !SpecUtils::is_file(path) )
      continue;

    const GadrasDetectorDat dat = GadrasDetectorDat::fromFile( path );

    vector<string> warnings;
    try
    {
      const ceelo::GeometryDescriptor gd = CeeLoUtils::buildGadrasGeometry( dat, warnings );

      const vector<ceelo::GeometryProblem> problems = gd.problems();
      for( const ceelo::GeometryProblem p : problems )
        BOOST_ERROR( name + ": geometry problem " + ceelo::to_string(p) );

      // Every one of these is a right circular cylinder in reality; a Box here
      //  means the material hint failed and the shape-factor fallback ran.
      BOOST_CHECK_MESSAGE( gd.shape == ceelo::DetectorShape::Cylinder,
                          name + ": expected a cylinder, got a box - the crystal"
                          " material was probably not resolved" );

      if( dat.setbackCm() > 0.0 )
        BOOST_CHECK_MESSAGE( close_enough( gd.endcap_front_offset_cm(), dat.setbackCm(), 1.0e-3 ),
                            name + ": endcap_front_offset "
                            + std::to_string(gd.endcap_front_offset_cm())
                            + " != setback " + std::to_string(dat.setbackCm()) );
      built += 1;
    }catch( std::exception &e )
    {
      // "HPGe 40%" has a crystal length of zero, so it cannot be modeled; it
      //  must fail cleanly rather than build a zero-length crystal.
      BOOST_TEST_MESSAGE( name + ": Mode A unavailable - " + e.what() );
      BOOST_CHECK_MESSAGE( (dat.length() <= 0.0f) || (dat.width() <= 0.0f),
                          name + ": refused for a reason other than bad dimensions: "
                          + e.what() );
      refused += 1;
    }
  }//for( const string &dir : dirs )

  BOOST_TEST_MESSAGE( "Shipped GADRAS detectors: " + std::to_string(built)
                      + " modeled, " + std::to_string(refused) + " refused." );
  BOOST_CHECK( built >= 12 );
}


BOOST_AUTO_TEST_CASE( test_gadras_solid_angle_convention )
{
  // GADRAS's own fractional solid angle (param 9) is computed at the source
  //  distance PLUS the setback, over the round-face-equivalent diameter.  That
  //  is the convention the Efficiency.csv cross-validation compares through, so
  //  it is pinned here.
  //
  //  Measured agreement (|computed/stated - 1|), which is what the 5% gate is
  //  set from - and note how much the setback term carries:
  //      NaI_3x3_text          0.00%   (0.60% without the setback)
  //      Kromek_GR1_CZT_text   0.02%   (no setback in that file)
  //      CZT_1cm_text          2.15%   (no setback in that file)
  //      CZT_1.5x2x2_text      3.07%   (15.81% without the setback)
  //  Param 9 is a stored, fitted GADRAS parameter rather than a recomputed
  //  disk solid angle, so a couple of percent of file-specific slop is
  //  expected; dropping the setback is what produces the double-digit errors.
const vector<string> dets = { "NaI_3x3_text", "NaI_2x2_text", "HPGe_Planar50_text",
                                "CZT_1cm_text", "Kromek_GR1_CZT_text" };
  for( const string &det : dets )
  {
    const GadrasDetectorDat dat = load( det );

    const double stated = 0.01 * dat.solidAnglePercent();
    if( stated <= 0.0 )
      continue;

    const double diam = dat.equivalentCircularDiameterCm() * PhysicalUnits::cm;
    const double dist = (dat.distanceCm() + dat.setbackCm()) * PhysicalUnits::cm;
    BOOST_REQUIRE( (diam > 0.0) && (dist > 0.0) );

    const double computed = DetectorPeakResponse::fractionalSolidAngle( diam, dist );

    BOOST_CHECK_MESSAGE( std::fabs(computed - stated) <= 0.05*stated,
                        det + ": GADRAS solid angle " + std::to_string(stated)
                        + " vs fractionalSolidAngle(D_eq, dist+setback) "
                        + std::to_string(computed) );

    // ...and without the setback it should be measurably worse, which is what
    //  makes the setback term load-bearing rather than cosmetic.
    if( dat.setbackCm() > 0.05f )
    {
      const double no_setback = DetectorPeakResponse::fractionalSolidAngle(
                                    diam, dat.distanceCm()*PhysicalUnits::cm );
      BOOST_TEST_MESSAGE( det + ": solid angle w/ setback "
                          + std::to_string(computed) + ", w/o "
                          + std::to_string(no_setback) + ", GADRAS "
                          + std::to_string(stated) );
    }
  }//for( const string &det : dets )
}


/** Every DRF InterSpec ships must come out of a refactor byte-identical.

 `fromGadrasDefinition` is about to be split into a Detector.dat half and an
 Efficiency.csv half so a `.dat` can be imported on its own; the hash covers the
 efficiency curve, FWHM, diameter, setback, energy range and geometry type, so a
 mismatch here means the split changed a shipped detector.

 Recording rather than asserting fixed values: the point is that the number does
 not MOVE, and pinning literals would just have to be rewritten whenever the
 hash's own definition changes.  Run with --record-gadras-hashes to print the
 current set.
 */
BOOST_AUTO_TEST_CASE( test_shipped_gadras_drfs_unchanged )
{
  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Need --datadir" );
  const string base = SpecUtils::append_path( g_data_dir, "GenericGadrasDetectors" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory(base), "Missing " + base );

  // name -> hash, re-recorded 2026-08-29 (--record-gadras-hashes), when `applyGadrasDat` began
  //  recording the crystal geometry on the DRF and attaching a measured-curve transfer response to
  //  it.  Both feed computeHash(), so every shipped detector that HAS an Efficiency.csv and a
  //  usable geometry changed hash - deliberately: they answer off-axis and near-field now.
  //
  //  "HPGe 40%" is the one entry that did NOT change: its Detector.dat states a zero crystal length
  //  (parameter 10), so no geometry can be built for it and it keeps the flat-disk treatment.  If
  //  that data file is ever corrected, this hash moves too.
  static const std::map<string,uint64_t> sm_expected = {
    { "HPGe 10%", 12999987909897955040ull },
    { "HPGe 20%", 15942574301746852451ull },
    { "HPGe 40%", 11793998863736797277ull },
    { "LaBr 10%", 3405271502809736741ull },
    { "LaBr 5%", 9533508627378423206ull },
    { "NaI 10%", 12710450218362030861ull },
    { "NaI 12%", 13838579886097982910ull },
    { "NaI 1x1", 9444249551207657343ull },
    { "NaI 25%", 14359049774689294627ull },
    { "NaI 2x2", 17708187484627074873ull },
    { "NaI 30%", 3282691297936797345ull },
    { "NaI 3x3", 12212646696940897187ull },
    { "NaI 5%", 7940542305173954851ull }
  };

  bool record = false;
  {
    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;
    for( int i = 1; i < argc; ++i )
      if( string(argv[i]) == "--record-gadras-hashes" )
        record = true;
  }

  vector<string> dirs = SpecUtils::ls_directories_in_directory( base );
  std::sort( begin(dirs), end(dirs) );
  BOOST_REQUIRE( dirs.size() >= 13 );

  size_t checked = 0;
  for( const string &dir : dirs )
  {
    const string name = SpecUtils::filename( dir );
    if( !SpecUtils::is_file( SpecUtils::append_path(dir, "Efficiency.csv") ) )
      continue;

    DetectorPeakResponse drf;
    BOOST_REQUIRE_NO_THROW( drf.fromGadrasDirectory( dir ) );
    BOOST_REQUIRE_MESSAGE( drf.isValid(), name + ": DRF is not valid" );

    const uint64_t hash = drf.hashValue();
    if( record )
    {
      cout << "    { \"" << name << "\", " << hash << "ull }," << endl;
      continue;
    }

    const auto pos = sm_expected.find( name );
    BOOST_REQUIRE_MESSAGE( pos != end(sm_expected),
                           name + " is not in the recorded set - add it" );
    BOOST_CHECK_MESSAGE( hash == pos->second,
                        name + ": hash changed, " + std::to_string(pos->second)
                        + " -> " + std::to_string(hash)
                        + ".  A shipped detector response is no longer what it was." );
    checked += 1;
  }//for( const string &dir : dirs )

  if( !record )
    BOOST_CHECK_MESSAGE( checked >= 13, "Only checked " + std::to_string(checked)
                         + " shipped detectors" );
}//test_shipped_gadras_drfs_unchanged


/** A Detector.dat with no Efficiency.csv beside it.

 `IdentiFINDER_LaBr3_xml` and `NaI_3x3_text` are the two fixtures shipped without
 an Efficiency.csv, so they are what a geometry-only import actually looks like.
 */
BOOST_AUTO_TEST_CASE( test_gadras_dat_only_import )
{
  const vector<string> dat_only = { "IdentiFINDER_LaBr3_xml", "NaI_3x3_text" };

  for( const string &det : dat_only )
  {
    const string dir = SpecUtils::append_path( gadras_dir(), det );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory(dir), "Missing fixture " + dir );
    BOOST_REQUIRE_MESSAGE( !SpecUtils::is_file( SpecUtils::append_path(dir,"Efficiency.csv") ),
                           det + " unexpectedly has an Efficiency.csv" );

    // The default still refuses, so no existing caller can be handed a detector
    //  that has no efficiency.
    {
      DetectorPeakResponse strict;
      BOOST_CHECK_THROW( strict.fromGadrasDirectory( dir ), std::exception );
    }

    DetectorPeakResponse drf;
    BOOST_REQUIRE_NO_THROW( drf.fromGadrasDirectory( dir, true ) );

    // Deliberately not valid: it has a shape but no sensitivity yet.
    BOOST_CHECK_MESSAGE( !drf.isValid(),
                         det + ": a geometry-only DRF must not report itself valid" );
    BOOST_CHECK( drf.drfSource() == DetectorPeakResponse::DrfSource::GadrasDetectorDatOnly );

    // ...but everything the file DOES define must have come across.
    const GadrasDetectorDat dat = GadrasDetectorDat::fromFile(
                                    SpecUtils::append_path(dir, "Detector.dat") );
    BOOST_CHECK_MESSAGE( drf.detectorDiameter() > 0.0f, det + ": no crystal diameter" );
    BOOST_CHECK( close_enough( drf.detectorDiameter() / PhysicalUnits::cm,
                               dat.equivalentCircularDiameterCm(), 0.01 ) );
    BOOST_CHECK( close_enough( drf.detectorSetback() / PhysicalUnits::cm,
                               dat.setbackCm(), 0.01 ) );

    BOOST_CHECK( drf.resolutionFcnType() == DetectorPeakResponse::kGadrasResolutionFcn );
    BOOST_CHECK_MESSAGE( drf.peakResolutionFWHM( 661.0f ) > 0.0f,
                         det + ": no usable FWHM at 661 keV" );

    // A GADRAS characterization always states a peak shape - "no skew" is one of
    //  the shapes it can state - so the preferences come across either way.
    const shared_ptr<const PeakFitDetPrefs> prefs = drf.peakFitDetPrefs();
    BOOST_REQUIRE_MESSAGE( !!prefs, det + ": no GADRAS peak-shape preferences" );
    BOOST_CHECK( (prefs->m_peak_skew_type == PeakDef::SkewType::GadrasGeneric)
                 || (prefs->m_peak_skew_type == PeakDef::SkewType::GadrasCZT) );

    // GADRAS fit this shape against the real detector, so all six coefficients
    //  are FIXED - a nullopt would mean "fit it per-ROI" and would throw that
    //  measurement away.
    for( int i = 0; i < 6; ++i )
      BOOST_CHECK_MESSAGE( prefs->m_lower_energy_skew[i].has_value(),
                          det + ": skew parameter " + std::to_string(i)
                          + " was left to be fit rather than taken from the file" );

    // The description carries the provenance and the low-energy caveat, since
    //  there is no measured curve to speak for this detector.
    BOOST_CHECK( SpecUtils::icontains( drf.description(), "Monte-Carlo" ) );
    BOOST_CHECK( SpecUtils::icontains( drf.description(), "60 keV" ) );

    BOOST_TEST_MESSAGE( "  " << det << ": " << drf.description() );
  }//for( const string &det : dat_only )
}



/** The sniffer must not call a spectrum or an efficiency CSV a Detector.dat -
 the text parser itself is far too permissive to use for classification. */
BOOST_AUTO_TEST_CASE( test_is_candidate_detector_dat )
{
  for( const string &det : { "NaI_3x3_text", "Detective_X_xml", "CZT_1cm_text" } )
  {
    ifstream in( dat_path(det).c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( in.is_open() );
    BOOST_CHECK_MESSAGE( GadrasDetectorDat::isCandidateDetectorDat( in ),
                         det + ": a real Detector.dat was not recognized" );
    // The stream must be left usable.
    BOOST_CHECK( !!GadrasDetectorDat::fromStream( in ).width() );
  }

  const string effcsv = SpecUtils::append_path(
        SpecUtils::append_path( gadras_dir(), "Detective_X_xml" ), "Efficiency.csv" );
  if( SpecUtils::is_file(effcsv) )
  {
    ifstream in( effcsv.c_str(), ios::in | ios::binary );
    BOOST_CHECK_MESSAGE( !GadrasDetectorDat::isCandidateDetectorDat( in ),
                         "an Efficiency.csv was mistaken for a Detector.dat" );
  }

  {
    stringstream empty( "" );
    BOOST_CHECK( !GadrasDetectorDat::isCandidateDetectorDat( empty ) );

    stringstream prose( "This is just some text.\n1 2 3\n" );
    BOOST_CHECK( !GadrasDetectorDat::isCandidateDetectorDat( prose ) );
  }

  // Leading '!' / '#' comment lines must not hide the table from the sniffer -
  //  the drag-and-drop classifier applies the same rule on its header buffer.
  {
    string commented = "! a comment GADRAS might write\n# another\n";
    ifstream real( dat_path("NaI_3x3_text").c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( real.is_open() );
    stringstream body;
    body << real.rdbuf();
    commented += body.str();

    stringstream withcomments( commented );
    BOOST_CHECK_MESSAGE( GadrasDetectorDat::isCandidateDetectorDat( withcomments ),
                         "leading comments hid a real Detector.dat from the sniffer" );
  }
}


BOOST_AUTO_TEST_CASE( test_generic_attenuator_name_round_trip )
{
  // The geometry form shows a GADRAS attenuator as this text and rebuilds it
  //  from whatever the user leaves there, so the two must agree exactly.
  for( const double an : { 6.0, 13.0, 13.2, 26.0, 82.0 } )
  {
    for( const double ad : { 0.05, 0.5, 1.35, 12.0 } )
    {
      const string text = CeeLoUtils::genericAttenuatorName( an, ad );
      double got_an = 0.0, got_ad = 0.0;
      BOOST_REQUIRE_MESSAGE( CeeLoUtils::parseGenericAttenuatorName( text, got_an, got_ad ),
                             "could not read back '" + text + "'" );
      BOOST_CHECK_CLOSE( got_an, an, 0.1 );
      BOOST_CHECK_CLOSE( got_ad, ad, 0.1 );

      // Whatever thickness it is given, the areal density must survive.
      const ceelo::MaterialSpec spec
              = CeeLoUtils::genericAttenuatorMaterial( got_an, got_ad, 0.3 );
      BOOST_CHECK_CLOSE( spec.density_g_per_cm3 * 0.3, ad, 0.01 );
    }
  }

  // Tolerant of how a user might retype it.
  double an = 0.0, ad = 0.0;
  BOOST_CHECK( CeeLoUtils::parseGenericAttenuatorName( "an=13, ad=0.5", an, ad ) );
  BOOST_CHECK_CLOSE( an, 13.0, 0.1 );
  BOOST_CHECK( CeeLoUtils::parseGenericAttenuatorName( "AN = 26  AD = 2.5 g/cm2", an, ad ) );
  BOOST_CHECK_CLOSE( ad, 2.5, 0.1 );

  // ...and must NOT swallow ordinary material names, which fall through to the
  //  MaterialDB lookup.
  for( const char *notgeneric : { "Aluminum", "Fe", "C0.5H0.5 d=1.2", "", "AN=0, AD=1" } )
    BOOST_CHECK_MESSAGE( !CeeLoUtils::parseGenericAttenuatorName( notgeneric, an, ad ),
                         string("'") + notgeneric + "' was mistaken for a generic attenuator" );
}


/** Every GADRAS attenuator must survive the geometry form: the descriptor built
 from a Detector.dat, rendered into the form, and read back out again has to
 carry the same areal densities and the same crystal.  This is the round trip
 that used to turn an imported BGO detector into NaI. */
BOOST_AUTO_TEST_CASE( test_gadras_geometry_names_are_recoverable )
{
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  for( const string &det : { "NaI_3x3_text", "Detective_X_xml", "CZT_1cm_text" } )
  {
    const GadrasDetectorDat dat = load( det );
    vector<string> warnings;
    const ceelo::GeometryDescriptor gd = CeeLoUtils::buildGadrasGeometry( dat, warnings );

    // The crystal must be a material the geometry form actually offers, or it
    //  will be substituted when the user opens the form.
    BOOST_REQUIRE( gd.crystal_material_index >= 0 );
    const string crystal = gd.materials[static_cast<size_t>(gd.crystal_material_index)].name;
    BOOST_CHECK_MESSAGE( CeeLoUtils::gadrasCrystalMaterial( crystal ).name == crystal,
                        det + ": crystal '" + crystal + "' does not resolve back to itself" );

    // Every synthesized attenuator must read back as the AN/AD it was built
    //  from, at whatever thickness it was given.
    for( const ceelo::LayerSpec &layer : gd.layers )
    {
      const ceelo::MaterialSpec &m = gd.materials[static_cast<size_t>(layer.material_index)];
      double an = 0.0, ad = 0.0;
      if( !CeeLoUtils::parseGenericAttenuatorName( m.name, an, ad ) )
        continue;   //the transparent standoff spacer

      const double t = (layer.front_thickness_cm > 0.0) ? layer.front_thickness_cm
                                                        : layer.side_thickness_cm;
      BOOST_REQUIRE( t > 0.0 );

      // The text is what the user reads and edits, so it carries four
      //  significant figures rather than the full double.  That is only the
      //  precision of an EDITED row: a row the user leaves alone keeps the exact
      //  composition it was seeded with (DetectorGeometryInput's `use_seeded`
      //  path), so nothing is lost by opening the form and closing it.
      BOOST_CHECK_MESSAGE( std::fabs( m.density_g_per_cm3*t - ad ) <= 1.0e-3*ad,
                          det + ": layer '" + m.name + "' carries areal density "
                          + std::to_string(m.density_g_per_cm3*t) + ", not " + std::to_string(ad) );
    }
  }//for( const string &det : ... )
}


/** Every material GADRAS knows must classify, and the classes that are not a
 question of resolution must not be decided by one.

 The corpus only exercises four materials; this covers all 36 table entries, so
 a new one (or a re-ordered table) cannot quietly fall through to Unknown. */
BOOST_AUTO_TEST_CASE( test_coarse_type_for_every_gadras_material )
{
  const std::array<GadrasDetectorDat::MaterialInfo,36> &table
                                          = GadrasDetectorDat::materialTable();

  // Named expectations for the ones whose class is a physics statement rather
  //  than a resolution: semiconductors tail, germanium is HPGe, the lanthanum
  //  and cerium halides are their own class.
  const std::map<string,PeakFitUtils::CoarseResolutionType> expect = {
    { "CZT",   PeakFitUtils::CoarseResolutionType::CZT  },
    { "CdTe",  PeakFitUtils::CoarseResolutionType::CZT  },
    { "TlBr",  PeakFitUtils::CoarseResolutionType::CZT  },
    { "HgI",   PeakFitUtils::CoarseResolutionType::CZT  },
    { "HPGe",  PeakFitUtils::CoarseResolutionType::High },
    { "Si",    PeakFitUtils::CoarseResolutionType::High },
    { "LaBr3", PeakFitUtils::CoarseResolutionType::LaBr },
    { "LaCl3", PeakFitUtils::CoarseResolutionType::LaBr },
    { "CeBr3", PeakFitUtils::CoarseResolutionType::LaBr },
    { "NaI",   PeakFitUtils::CoarseResolutionType::Low  },
    { "CsI",   PeakFitUtils::CoarseResolutionType::Low  },
    { "BGO",   PeakFitUtils::CoarseResolutionType::Low  },
  };

  size_t unknown = 0;
  for( const GadrasDetectorDat::MaterialInfo &info : table )
  {
    DetectorPeakResponse drf;
    // Drive it through the real path: a minimal Detector.dat naming this
    //  material, with a resolution that would classify it differently, so the
    //  material must be what decides.
    GadrasDetectorDat dat;
    dat.setValue( 10, 5.0f );          //length
    dat.setValue( 11, 5.0f );          //width
    dat.setValue( 12, 1.0f );          //height/width
    dat.setValue( 7, info.resolution661 );
    dat.m_material_name = info.name;
    dat.m_format = GadrasDetectorDat::SourceFormat::Xml;

    stringstream xml;
    dat.toXml( xml );
    DetectorPeakResponse built;
    BOOST_REQUIRE_NO_THROW( built.fromGadrasDatOnly( xml ) );

    const shared_ptr<const PeakFitDetPrefs> prefs = built.peakFitDetPrefs();
    BOOST_REQUIRE_MESSAGE( !!prefs, string(info.name) + ": no preferences" );

    const auto pos = expect.find( info.name );
    if( pos != end(expect) )
      BOOST_CHECK_MESSAGE( prefs->m_det_type == pos->second,
                          string(info.name) + ": classified as "
                          + PeakFitDetPrefs::to_str(prefs->m_det_type) );

    if( prefs->m_det_type == PeakFitUtils::CoarseResolutionType::Unknown )
    {
      unknown += 1;
      BOOST_TEST_MESSAGE( "  unclassified: " << info.name
                          << " (nominal " << info.resolution661 << "% FWHM)" );
    }

    // A semiconductor gets the CZT tail construction; nothing else does.
    const bool want_czt = (prefs->m_det_type == PeakFitUtils::CoarseResolutionType::CZT);
    BOOST_CHECK_MESSAGE( (prefs->m_peak_skew_type == PeakDef::SkewType::GadrasCZT) == want_czt,
                        string(info.name) + ": peak-shape family does not follow its class" );
  }//for( every material )

  BOOST_CHECK_MESSAGE( unknown == 0, std::to_string(unknown)
                       + " GADRAS materials could not be classified" );
}


/** Every crystal a Detector.dat can name has to be one the geometry form offers,
 or opening the form would substitute a different crystal for it.

 The form's list is built from GadrasDetectorDat::materialTable(), so this holds
 by construction - but only as long as every table entry actually resolves.  A
 new table row whose formula does not parse would silently drop off the list and
 take its detectors' crystals with it. */
BOOST_AUTO_TEST_CASE( test_every_gadras_crystal_is_offered )
{
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  size_t offered = 0;
  for( const GadrasDetectorDat::MaterialInfo &info : GadrasDetectorDat::materialTable() )
  {
    ceelo::MaterialSpec spec;
    BOOST_REQUIRE_NO_THROW( spec = CeeLoUtils::gadrasCrystalMaterial( info.name ) );

    // The name the import stores must resolve back to the same material, which
    //  is what lets the form match it rather than substitute.
    ceelo::MaterialSpec again;
    BOOST_REQUIRE_NO_THROW( again = CeeLoUtils::gadrasCrystalMaterial( spec.name ) );
    BOOST_CHECK_MESSAGE( again.name == spec.name,
                        string(info.name) + " stores as '" + spec.name
                        + "', which resolves to '" + again.name + "'" );
    BOOST_CHECK( !spec.composition.empty() );
    offered += 1;
  }

  BOOST_CHECK_EQUAL( offered, GadrasDetectorDat::materialTable().size() );
  BOOST_TEST_MESSAGE( "  all " << offered << " GADRAS crystals resolve" );
}


/** A drag-and-dropped Detector.dat.

 The browser sends only the base name, which for this format is the useless
 "Detector.dat" - so there is no directory to recover the crystal material from.
 That must not be fatal: the geometry form has a crystal dropdown the user can
 set, so an assumed crystal plus a loud note beats refusing the import.
 (Reported from a real drag-and-drop of "LaBr 10%/Detector.dat", which failed
 with "the detector's geometry cannot be modeled".)
 */
BOOST_AUTO_TEST_CASE( test_dropped_detector_dat_has_no_directory )
{
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Need --datadir" );

  for( const char * const det : { "LaBr 10%", "NaI 3x3", "HPGe 20%" } )
  {
    const string path = SpecUtils::append_path(
        SpecUtils::append_path( SpecUtils::append_path(g_data_dir, "GenericGadrasDetectors"), det ),
        "Detector.dat" );
    if( !SpecUtils::is_file(path) )
      continue;

    ifstream in( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( in.is_open() );

    // Exactly what the drop path has to work with: the bytes, and the name
    //  "Detector.dat" - and, for this half of the test, a file that names no
    //  crystal, so there is genuinely nothing to go on.
    const GadrasDetectorDat dat = without_material(
                              GadrasDetectorDat::fromStream( in ) );
    BOOST_CHECK_MESSAGE( dat.materialName().empty(),
                        string(det) + ": nothing should identify this crystal" );

    vector<string> warnings;
    ceelo::GeometryDescriptor gd;
    BOOST_REQUIRE_NO_THROW( gd = CeeLoUtils::buildGadrasGeometry( dat, warnings ) );

    // A usable geometry, not an exception.
    BOOST_CHECK( gd.problems().empty() );
    BOOST_REQUIRE( gd.crystal_material_index >= 0 );
    BOOST_CHECK( !gd.dimensions_cm.empty() );

    // ...and the assumption has to be stated, not silent.
    bool said_so = false;
    for( const string &w : warnings )
      said_so |= SpecUtils::icontains( w, "was assumed" );
    BOOST_CHECK_MESSAGE( said_so, string(det) + ": assumed a crystal without saying so" );

    const string crystal = gd.materials[static_cast<size_t>(gd.crystal_material_index)].name;
    BOOST_TEST_MESSAGE( "  " << det << " dropped with no directory -> assumed "
                        << crystal << " (" << dat.resFWHM661() << "% FWHM@661)" );

    // The resolution at least has to separate germanium from the scintillators.
    if( dat.resFWHM661() < 1.5f )
      BOOST_CHECK_EQUAL( crystal, string("HPGe") );
    else
      BOOST_CHECK_EQUAL( crystal, string("NaI") );
  }//for( each detector )

  // ...but assuming is strictly the LAST resort.  A file that states its
  //  material - as a param-59 index or an XML <material> - must be believed even
  //  when dropped with no directory to go on.  (Only 2 of the 13 detectors
  //  InterSpec ships state one; the fixtures, which came from real GADRAS
  //  installs, mostly do.)
  struct Stated { const char *fixture; const char *material; };
  const vector<Stated> stated = {
    { "NaI_3x3_text",           "NaI"   },   //param-59 index 5
    { "CZT_1cm_text",           "CZT"   },   //param-59 index 8
    { "HPGe_Planar50_text",     "HPGe"  },   //param-59 index 3
    { "IdentiFINDER_LaBr3_xml", "LaBr3" },   //XML <material>
    { "MikesCzt_xml",           "CZT"   },   //XML <material>
  };

  for( const Stated &e : stated )
  {
    const string path = dat_path( e.fixture );
    (void)path;
    if( !SpecUtils::is_file(path) )
      continue;

    ifstream in( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( in.is_open() );
    const GadrasDetectorDat dat = GadrasDetectorDat::fromStream( in );

    vector<string> notes;
    BOOST_CHECK_MESSAGE( dat.materialName() == string(e.material),
                        string(e.fixture) + ": states '" + e.material
                        + "' but resolved to '" + dat.materialName() + "'" );
    BOOST_CHECK_MESSAGE( notes.empty(),
                        string(e.fixture) + ": a stated material must not be reported"
                        " as inferred" );

    vector<string> warnings;
    ceelo::GeometryDescriptor gd;
    BOOST_REQUIRE_NO_THROW( gd = CeeLoUtils::buildGadrasGeometry( dat, warnings ) );
    for( const string &w : warnings )
      BOOST_CHECK_MESSAGE( !SpecUtils::icontains( w, "was assumed" ),
                          string(e.fixture) + ": assumed a crystal despite the file"
                          " stating one" );

    BOOST_TEST_MESSAGE( "  " << e.fixture << " dropped with no directory -> "
                        << gd.materials[static_cast<size_t>(gd.crystal_material_index)].name
                        << " (from the file)" );
  }//for( const Stated &e : stated )

  // The detectors InterSpec ships now state their material too, so dropping one
  //  identifies its crystal from the file rather than assuming - which is what a
  //  real user's Detector.dat does, and what a dropped "LaBr 10%" failed to do
  //  before the material indices were filled in.
  for( const char * const det : { "LaBr 10%", "LaBr 5%", "NaI 3x3", "HPGe 20%" } )
  {
    const string path = SpecUtils::append_path(
        SpecUtils::append_path( SpecUtils::append_path(g_data_dir, "GenericGadrasDetectors"), det ),
        "Detector.dat" );
    if( !SpecUtils::is_file(path) )
      continue;

    ifstream in( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( in.is_open() );
    const GadrasDetectorDat dat = GadrasDetectorDat::fromStream( in );

    vector<string> notes;
    const string crystal = dat.materialName();
    BOOST_CHECK_MESSAGE( !crystal.empty(),
                        string(det) + ": the shipped file no longer identifies its crystal" );
    BOOST_CHECK_MESSAGE( notes.empty(),
                        string(det) + ": a stated material was reported as inferred" );
    BOOST_TEST_MESSAGE( "  " << det << " dropped with no directory -> " << crystal
                        << " (from the file)" );
  }
}


/** An Efficiency.csv imported with its Detector.dat beside it must end up with
 the .dat's FWHM, peak shape and crystal type on the finished DRF - the
 "Modify Detector Response" dialog reads the crystal back off those preferences,
 so losing them shows the detector as NaI whatever it really is.

 Reproduces what SpecMeasManager's efficiency-CSV dialog does on accept.
 */
BOOST_AUTO_TEST_CASE( test_eff_csv_plus_dat_keeps_detector_type )
{
  BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Need --datadir" );

  struct Expect { const char *dir; PeakFitUtils::CoarseResolutionType type; };
  const vector<Expect> cases = {
    { "LaBr 10%", PeakFitUtils::CoarseResolutionType::LaBr },
    { "HPGe 20%", PeakFitUtils::CoarseResolutionType::High },
    { "NaI 3x3",  PeakFitUtils::CoarseResolutionType::Low  },
  };

  for( const Expect &e : cases )
  {
    const string dir = SpecUtils::append_path(
                        SpecUtils::append_path(g_data_dir, "GenericGadrasDetectors"), e.dir );
    const string csvpath = SpecUtils::append_path( dir, "Efficiency.csv" );
    const string datpath = SpecUtils::append_path( dir, "Detector.dat" );
    if( !SpecUtils::is_file(csvpath) || !SpecUtils::is_file(datpath) )
      continue;

    // 1. The CSV alone, as the dialog first parses it.
    ifstream csv( csvpath.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( csv.is_open() );
    DetectorPeakResponse::EffCsvParseResult parsed;
    BOOST_REQUIRE_NO_THROW( parsed = DetectorPeakResponse::parseEfficiencyCsvFile( csv ) );
    BOOST_REQUIRE( parsed.drf );
    BOOST_CHECK_MESSAGE( !parsed.drf->peakFitDetPrefs(),
                        string(e.dir) + ": an efficiency CSV alone cannot know the crystal" );

    // 2. The Detector.dat uploaded alongside it.
    ifstream datf( datpath.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( datf.is_open() );
    auto enriched = make_shared<DetectorPeakResponse>();
    BOOST_REQUIRE_NO_THROW( enriched->fromGadrasDatOnly( datf ) );

    const shared_ptr<const PeakFitDetPrefs> from_dat = enriched->peakFitDetPrefs();
    BOOST_REQUIRE_MESSAGE( !!from_dat, string(e.dir) + ": the .dat gave no preferences" );
    BOOST_CHECK_MESSAGE( from_dat->m_det_type == e.type,
                        string(e.dir) + ": the .dat itself has the wrong detector type" );

    // 3. What the dialog builds on accept: reinterpret the CSV curve at the
    //    .dat's diameter, then carry the .dat's own values across.
    shared_ptr<DetectorPeakResponse> final_drf;
    BOOST_REQUIRE_NO_THROW( final_drf = parsed.drf->reinterpretAsFarFieldIntrinsicEfficiency(
                                                    enriched->detectorDiameter() ) );
    BOOST_REQUIRE( final_drf );
    final_drf->setPeakFitDetPrefs( from_dat );
    if( enriched->detectorSetback() > 0.0 )
      final_drf->setDetectorSetback( enriched->detectorSetback() );

    // 4. This is what the Modify dialog reads to choose the crystal.
    const shared_ptr<const PeakFitDetPrefs> prefs = final_drf->peakFitDetPrefs();
    BOOST_REQUIRE_MESSAGE( !!prefs, string(e.dir) + ": preferences lost building the final DRF" );
    BOOST_CHECK_MESSAGE( prefs->m_det_type == e.type,
                        string(e.dir) + ": final DRF has detector type "
                        + PeakFitDetPrefs::to_str(prefs->m_det_type)
                        + ", so the geometry form would show the wrong crystal" );

    BOOST_TEST_MESSAGE( "  " << e.dir << ": CSV+dat -> "
                        << PeakFitDetPrefs::to_str(prefs->m_det_type) );
  }//for( const Expect &e : cases )
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
