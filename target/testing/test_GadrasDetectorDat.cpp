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


BOOST_AUTO_TEST_CASE( test_resolve_crystal_from_name_hint )
{
  // Every detector shipped in data/GenericGadrasDetectors has a material index
  //  of 0, so without a name hint inferShape() falls through to its
  //  shape-factor heuristic and calls a 3x3 NaI a box.  This is the guard on
  //  that: the hint must recover the material, and the material must fix the
  //  shape.
  GadrasDetectorDat bare;
  bare.setValue( 10, 7.62f );    //length
  bare.setValue( 11, 6.35f );    //width
  bare.setValue( 12, 1.0f );     //height/width
  bare.setValue( 13, 100.0f );   //shape factor -> "Box" without a material
  BOOST_REQUIRE( bare.materialName().empty() );
  BOOST_CHECK( bare.inferShape().shape == GadrasDetectorDat::Shape::Box );

  vector<string> notes;
  const string nai = CeeLoUtils::resolveGadrasCrystalName( bare, "NaI 3x3", notes );
  BOOST_CHECK_EQUAL( nai, string("NaI") );
  BOOST_CHECK( !notes.empty() );   //a guessed material must always be reported
  BOOST_CHECK( bare.inferShape( nai ).shape == GadrasDetectorDat::Shape::Cylinder );

  notes.clear();
  BOOST_CHECK_EQUAL( CeeLoUtils::resolveGadrasCrystalName( bare, "HPGe 40%", notes ),
                     string("HPGe") );
  notes.clear();
  BOOST_CHECK_EQUAL( CeeLoUtils::resolveGadrasCrystalName( bare, "LaBr 10%", notes ),
                     string("LaBr3") );
  notes.clear();
  BOOST_CHECK_EQUAL( CeeLoUtils::resolveGadrasCrystalName( bare, "Kromek GR1 CZT", notes ),
                     string("CZT") );

  // Nothing to go on -> empty, so the caller can refuse the geometry import
  //  rather than model the wrong crystal.
  notes.clear();
  BOOST_CHECK( CeeLoUtils::resolveGadrasCrystalName( bare, "Detector 17", notes ).empty() );
  BOOST_CHECK( CeeLoUtils::resolveGadrasCrystalName( bare, "", notes ).empty() );

  // A material stated in the file always wins over the name.
  GadrasDetectorDat stated = bare;
  stated.m_material_name = "CsI";
  notes.clear();
  BOOST_CHECK_EQUAL( CeeLoUtils::resolveGadrasCrystalName( stated, "NaI 3x3", notes ),
                     string("CsI") );
  BOOST_CHECK( notes.empty() );
}


BOOST_AUTO_TEST_CASE( test_build_gadras_geometry )
{
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );

  struct Expect { const char *dir; const char *hint; ceelo::DetectorShape shape; };
  const vector<Expect> cases = {
    { "NaI_3x3_text",          "NaI 3x3",        ceelo::DetectorShape::Cylinder },
    { "NaI_2x2_text",          "NaI 2x2",        ceelo::DetectorShape::Cylinder },
    { "CZT_1cm_text",          "CZT 1cm",        ceelo::DetectorShape::Box },
    { "CZT_1.5x2x2_text",      "CZT 1.5x2x2",    ceelo::DetectorShape::Box },
    { "Kromek_GR1_CZT_text",   "Kromek GR1 CZT", ceelo::DetectorShape::Box },
    { "HPGe_Planar50_text",    "HPGe Planar50",  ceelo::DetectorShape::Cylinder },
    { "Detective_X_xml",       "Detective X",    ceelo::DetectorShape::Cylinder },
    { "IdentiFINDER_NGH_xml",  "IdentiFINDER NGH", ceelo::DetectorShape::Cylinder },
    { "IdentiFINDER_LaBr3_xml","IdentiFINDER LaBr3", ceelo::DetectorShape::Cylinder },
    { "MikesCzt_xml",          "MikesCzt",       ceelo::DetectorShape::Box },
  };

  for( const Expect &e : cases )
  {
    const GadrasDetectorDat dat = load( e.dir );

    vector<string> warnings;
    ceelo::GeometryDescriptor gd;
    BOOST_REQUIRE_NO_THROW( gd = CeeLoUtils::buildGadrasGeometry( dat, e.hint, warnings ) );

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
      const ceelo::GeometryDescriptor gd = CeeLoUtils::buildGadrasGeometry( dat, name, warnings );

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
